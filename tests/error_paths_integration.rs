//! Integration tests for inputs alevin-fry must reject with an error rather
//! than a panic.
//!
//! Release builds set `panic = "abort"`, so every panic on user input is a
//! SIGABRT with a core dump and a "panicked at" frame instead of a one-line
//! error and exit status 1. Each test below feeds one such input: a RAD file
//! that is too new, truncated or empty, a resolution or bootstrap request USA
//! mode cannot honour, a multi-sample knee run where a sample has only one or
//! two cell barcodes, and `view` / `collate` invocations that used to end in an
//! `unwrap` or `expect`. The library entry points are called directly where
//! they return the error; the binary is run where the error has to make it out
//! of `main`. Test builds unwind instead of aborting, so a regression shows up
//! as exit status 101 there.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use libradicl::chunk::Chunk;
use libradicl::codec::ChunkCodec;
use libradicl::header::{RadHeader, RadPrelude, SpecVersion};
use libradicl::rad_types::{
    RadIntId, RadType, TagDesc, TagMap, TagSection, TagSectionLabel, TagValue,
};
use libradicl::record::{
    AlevinFryReadRecord, AlevinFryRecordContext, MultiBarcodeReadRecord, MultiBarcodeRecordContext,
    RecordContext,
};
use libradicl::writers::RadFileWriter;
use smallvec::smallvec;

use alevin_fry::cellfilter::{CellFilterMethod, generate_permit_list};
use alevin_fry::collate::collate;
use alevin_fry::prog_opts::{GenPermitListOpts, QuantOpts, SampleCorrectionMode};
use alevin_fry::quant::{ResolutionStrategy, SplicedAmbiguityModel, quantify};
use bio_types::strand::Strand;

/// Spliced/unspliced reference pairs: `tx{2g}` is gene `g` spliced and
/// `tx{2g+1}` its unspliced counterpart.
const NUM_GENES: u64 = 3;
const NUM_REFS: u64 = 2 * NUM_GENES;
const CELL_BC_LEN: u16 = 16;
const UMI_LEN: u16 = 12;
const SAMPLE_BC_LEN: u16 = 8;
// The real version: quant rejects multi-sample collations written before 0.15.0.
const TEST_VERSION: &str = env!("CARGO_PKG_VERSION");

fn packed(idx: u64, len: u16) -> u64 {
    let mask = (1u64 << (2 * len as u64)) - 1;
    idx.wrapping_mul(2654435761) & mask
}

/// Render a packed barcode as nucleotides, first base in the most significant
/// bit pair (the `needletail::bitkmer` layout the barcode lists are read with).
fn packed_to_nuc(packed: u64, len: usize) -> String {
    const NUCS: [char; 4] = ['A', 'C', 'G', 'T'];
    (0..len)
        .map(|i| NUCS[((packed >> (2 * (len - 1 - i))) & 3) as usize])
        .collect()
}

fn make_test_logger() -> slog::Logger {
    if std::env::var("AF_TEST_VERBOSE").is_ok() {
        let decorator = slog_term::PlainDecorator::new(slog_term::TestStdoutWriter);
        let drain = slog_term::CompactFormat::new(decorator).build();
        let drain = std::sync::Mutex::new(drain);
        let drain = slog::Fuse::new(drain);
        return slog::Logger::root(drain, slog::o!());
    }
    slog::Logger::root(slog::Discard, slog::o!())
}

/// The full error chain, as `main` prints it.
fn chain(err: &anyhow::Error) -> String {
    format!("{err:#}")
}

/// Write a versioned (spec major 2) single-barcode RAD with `b`/`u` read tags.
/// `reads_per_cell[i]` reads carry cell barcode `i`, cycling over the references.
fn write_single_bc_rad(path: &Path, reads_per_cell: &[usize]) {
    let hdr = RadHeader {
        version: SpecVersion::current(),
        is_paired: 0,
        ref_count: NUM_REFS,
        ref_names: (0..NUM_REFS).map(|i| format!("tx{i}")).collect(),
        num_chunks: 0, // backpatched by RadFileWriter
    };
    let mut file_tags = TagSection::new_with_label(TagSectionLabel::FileTags);
    file_tags.add_tag_desc(TagDesc::new("cblen", RadType::Int(RadIntId::U16)));
    file_tags.add_tag_desc(TagDesc::new("ulen", RadType::Int(RadIntId::U16)));
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    read_tags.add_tag_desc(TagDesc::new("b", RadType::Int(RadIntId::U32)));
    read_tags.add_tag_desc(TagDesc::new("u", RadType::Int(RadIntId::U32)));
    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    aln_tags.add_tag_desc(TagDesc::new(
        "compressed_ori_refid",
        RadType::Int(RadIntId::U32),
    ));
    let prelude = RadPrelude::from_header_and_tag_sections(hdr, file_tags, read_tags, aln_tags);
    let mut file_tag_map = TagMap::with_keyset(&prelude.file_tags.tags);
    file_tag_map.add(TagValue::U16(CELL_BC_LEN));
    file_tag_map.add(TagValue::U16(UMI_LEN));

    let ctx = AlevinFryRecordContext::get_context_from_tag_section(
        &prelude.file_tags,
        &prelude.read_tags,
        &prelude.aln_tags,
    )
    .unwrap();
    let mut fw = RadFileWriter::new(File::create(path).unwrap(), &prelude, &file_tag_map).unwrap();
    for (cell, &nreads) in reads_per_cell.iter().enumerate() {
        let reads: Vec<_> = (0..nreads)
            .map(|r| AlevinFryReadRecord {
                bc: packed(cell as u64 + 1, CELL_BC_LEN),
                umi: packed((cell * 1000 + r) as u64, UMI_LEN),
                dirs: vec![true],
                refs: vec![(r as u64 % NUM_REFS) as u32],
            })
            .collect();
        let chunk = Chunk::<AlevinFryReadRecord> {
            nbytes: 0,
            nrec: reads.len() as u32,
            reads,
        };
        fw.write_chunk(&chunk, &ctx).unwrap();
    }
    fw.finalize().unwrap();
}

/// Write a two-level (sample, cell) multi-barcode RAD. `samples[s]` lists the
/// read count of each cell of sample `s`; sample `s` has sample barcode
/// `sample_bc(s)`.
fn write_multi_bc_rad(path: &Path, samples: &[Vec<usize>]) {
    let hdr = RadHeader {
        version: SpecVersion::Legacy,
        is_paired: 0,
        ref_count: NUM_REFS,
        ref_names: (0..NUM_REFS).map(|i| format!("tx{i}")).collect(),
        num_chunks: 0,
    };
    let mut file_tags = TagSection::new_with_label(TagSectionLabel::FileTags);
    for name in ["num_barcodes", "b0len", "b1len", "ulen"] {
        file_tags.add_tag_desc(TagDesc::new(name, RadType::Int(RadIntId::U16)));
    }
    file_tags.add_tag_desc(TagDesc::new("known_rad_type", RadType::String));
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    for name in ["b0", "b1", "u"] {
        read_tags.add_tag_desc(TagDesc::new(name, RadType::Int(RadIntId::U32)));
    }
    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    aln_tags.add_tag_desc(TagDesc::new(
        "compressed_ori_refid",
        RadType::Int(RadIntId::U32),
    ));
    let prelude = RadPrelude::from_header_and_tag_sections(hdr, file_tags, read_tags, aln_tags);
    let mut file_tag_map = TagMap::with_keyset(&prelude.file_tags.tags);
    file_tag_map.add(TagValue::U16(2));
    file_tag_map.add(TagValue::U16(SAMPLE_BC_LEN));
    file_tag_map.add(TagValue::U16(CELL_BC_LEN));
    file_tag_map.add(TagValue::U16(UMI_LEN));
    file_tag_map.add(TagValue::String("sc_rna_multi_bc".to_string()));

    let ctx = MultiBarcodeRecordContext::get_context_from_tag_section(
        &prelude.file_tags,
        &prelude.read_tags,
        &prelude.aln_tags,
    )
    .unwrap();
    let mut fw = RadFileWriter::new(File::create(path).unwrap(), &prelude, &file_tag_map).unwrap();
    for (s, cells) in samples.iter().enumerate() {
        for (cell, &nreads) in cells.iter().enumerate() {
            let cell_bc = packed((s * 1000 + cell) as u64 + 1, CELL_BC_LEN);
            let reads: Vec<_> = (0..nreads)
                .map(|r| MultiBarcodeReadRecord {
                    barcodes: smallvec![sample_bc(s), cell_bc],
                    umi: packed((s * 100_000 + cell * 100 + r) as u64, UMI_LEN),
                    dirs: vec![true],
                    refs: vec![(r as u64 % NUM_REFS) as u32],
                })
                .collect();
            let chunk = Chunk::<MultiBarcodeReadRecord> {
                nbytes: 0,
                nrec: reads.len() as u32,
                reads,
            };
            fw.write_chunk(&chunk, &ctx).unwrap();
        }
    }
    fw.finalize().unwrap();
}

fn sample_bc(s: usize) -> u64 {
    packed(100 * (s as u64 + 1), SAMPLE_BC_LEN)
}

/// A named sample list (`barcode<TAB>name`) for samples `0..n`.
fn write_sample_list(path: &Path, n: usize) {
    let mut f = BufWriter::new(File::create(path).unwrap());
    for s in 0..n {
        writeln!(
            f,
            "{}\tsample{s}",
            packed_to_nuc(sample_bc(s), SAMPLE_BC_LEN as usize)
        )
        .unwrap();
    }
}

fn single_gpl_opts<'a>(
    rad_dir: &'a PathBuf,
    out_dir: &'a PathBuf,
    fmeth: CellFilterMethod,
    log: &'a slog::Logger,
) -> GenPermitListOpts<'a, 'a, 'a, 'a, 'a> {
    GenPermitListOpts::builder()
        .input_dir(rad_dir)
        .output_dir(out_dir)
        .fmeth(fmeth)
        .expected_ori(Strand::Unknown)
        .version(TEST_VERSION)
        .threads(2)
        .velo_mode(false)
        .cmdline("test")
        .log(log)
        .build()
}

fn multi_gpl_opts<'a>(
    rad_dir: &'a PathBuf,
    out_dir: &'a PathBuf,
    fmeth: CellFilterMethod,
    sample_list: PathBuf,
    log: &'a slog::Logger,
) -> GenPermitListOpts<'a, 'a, 'a, 'a, 'a> {
    GenPermitListOpts::builder()
        .input_dir(rad_dir)
        .output_dir(out_dir)
        .fmeth(fmeth)
        .expected_ori(Strand::Unknown)
        .version(TEST_VERSION)
        .threads(2)
        .velo_mode(false)
        .cmdline("test")
        .log(log)
        .sample_bc_list(Some(sample_list))
        .sample_names(None)
        .sample_correction_mode(SampleCorrectionMode::Exact)
        .build()
}

/// A RAD written by a newer spec major, a truncated one and an empty one are
/// all errors from generate-permit-list that name the file, and the version
/// error keeps libradicl's "newer than supported" wording.
#[test]
fn gpl_rejects_unreadable_rad_header_with_an_error() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();
    let good = tmp.path().join("good.rad");
    write_single_bc_rad(&good, &[5, 5, 5]);
    let bytes = std::fs::read(&good).unwrap();
    assert_eq!(&bytes[..8], b"RAD_FILE");

    let mut too_new = bytes.clone();
    too_new[8] = 3; // spec major
    let cases: [(&str, Vec<u8>, &str); 3] = [
        (
            "too_new",
            too_new,
            "RAD spec major version 3 is newer than supported",
        ),
        (
            "truncated",
            bytes[..12].to_vec(),
            "failed to fill whole buffer",
        ),
        ("empty", Vec::new(), "failed to fill whole buffer"),
    ];
    for (name, data, cause) in cases {
        let rad_dir = tmp.path().join(name);
        std::fs::create_dir_all(&rad_dir).unwrap();
        std::fs::write(rad_dir.join("map.rad"), data).unwrap();
        let out_dir = tmp.path().join(format!("{name}_out"));
        let err = generate_permit_list(single_gpl_opts(
            &rad_dir,
            &out_dir,
            CellFilterMethod::ForceCells(3),
            &log,
        ))
        .expect_err(name);
        let msg = chain(&err);
        assert!(
            msg.contains(&format!("{}", rad_dir.join("map.rad").display())),
            "{name}: {msg}"
        );
        assert!(msg.contains(cause), "{name}: {msg}");
    }
}

/// USA mode (a 3-column tg-map) supports neither trivial resolution nor
/// bootstrapping; both are refused before any output is written, and an
/// allowed resolution still quantifies.
#[test]
fn usa_mode_rejects_trivial_resolution_and_bootstraps() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();
    let rad_dir = tmp.path().join("rad");
    let gpl_dir = tmp.path().join("gpl");
    std::fs::create_dir_all(&rad_dir).unwrap();
    write_single_bc_rad(&rad_dir.join("map.rad"), &[12, 9, 7]);
    generate_permit_list(single_gpl_opts(
        &rad_dir,
        &gpl_dir,
        CellFilterMethod::ForceCells(3),
        &log,
    ))
    .unwrap();
    collate(
        gpl_dir.clone(),
        &rad_dir,
        2,
        1_000,
        ChunkCodec::None,
        "test",
        TEST_VERSION,
        &log,
    )
    .unwrap();

    let usa_tg_map = tmp.path().join("t2g_3col.tsv");
    let mut f = BufWriter::new(File::create(&usa_tg_map).unwrap());
    for g in 0..NUM_GENES {
        writeln!(f, "tx{}\tgene{g}\tS", 2 * g).unwrap();
        writeln!(f, "tx{}\tgene{g}\tU", 2 * g + 1).unwrap();
    }
    drop(f);

    let quant = |resolution: ResolutionStrategy, num_bootstraps: u32, out: &PathBuf| {
        quantify(
            QuantOpts::builder()
                .input_dir(&gpl_dir)
                .tg_map(&usa_tg_map)
                .output_dir(out)
                .num_threads(2)
                .num_bootstraps(num_bootstraps)
                .init_uniform(false)
                .summary_stat(false)
                .dump_eq(false)
                .resolution(resolution)
                .pug_exact_umi(false)
                .sa_model(SplicedAmbiguityModel::WinnerTakeAll)
                .small_thresh(0)
                .large_graph_thresh(0)
                .filter_list(None)
                .cmdline("test")
                .version(TEST_VERSION)
                .log(&log)
                .build(),
        )
    };

    let out = tmp.path().join("q_trivial");
    let err = quant(ResolutionStrategy::Trivial, 0, &out).expect_err("trivial in USA mode");
    assert!(
        chain(&err).contains("--resolution trivial is not supported in USA mode"),
        "{err:#}"
    );
    assert!(!out.exists(), "a refused run must not leave partial output");

    let out = tmp.path().join("q_boot");
    let err =
        quant(ResolutionStrategy::CellRangerLikeEm, 4, &out).expect_err("bootstraps in USA mode");
    assert!(chain(&err).contains("bootstrapping"), "{err:#}");
    assert!(!out.exists(), "a refused run must not leave partial output");

    let out = tmp.path().join("q_crlike");
    quant(ResolutionStrategy::CellRangerLike, 0, &out).unwrap();
    assert!(out.join("alevin").join("quants_mat.mtx").exists());
}

/// In multi-sample knee mode a sample with only one or two distinct cell
/// barcodes used to trip knee finding's assertion and abort the whole run. It
/// is now skipped like an empty sample (no cells, no permit files), the other
/// samples are filtered as before, and other filter methods still keep those
/// barcodes.
#[test]
fn multi_sample_knee_skips_samples_with_too_few_barcodes() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();
    let rad_dir = tmp.path().join("rad");
    std::fs::create_dir_all(&rad_dir).unwrap();
    // sample0: a clear knee (8 large cells, then a tail of small ones);
    // sample1: a single cell barcode; sample2: two.
    let mut sample0 = vec![200; 8];
    sample0.extend([3, 2, 2, 1, 1, 1, 1, 1]);
    write_multi_bc_rad(&rad_dir.join("map.rad"), &[sample0, vec![40], vec![30, 20]]);
    let sample_list = tmp.path().join("samples.tsv");
    write_sample_list(&sample_list, 3);

    let read_info = |dir: &PathBuf| -> serde_json::Value {
        serde_json::from_reader(File::open(dir.join("sample_info.json")).unwrap()).unwrap()
    };
    let sample = |info: &serde_json::Value, name: &str| -> serde_json::Value {
        info["samples"]
            .as_array()
            .unwrap()
            .iter()
            .find(|s| s["name"] == name)
            .unwrap_or_else(|| panic!("no {name} in sample_info.json"))
            .clone()
    };

    let knee_dir = tmp.path().join("knee");
    generate_permit_list(multi_gpl_opts(
        &rad_dir,
        &knee_dir,
        CellFilterMethod::KneeFinding,
        sample_list.clone(),
        &log,
    ))
    .unwrap();
    let info = read_info(&knee_dir);
    let kept = sample(&info, "sample0")["num_cells"].as_u64().unwrap();
    assert!(
        (1..16).contains(&kept),
        "knee finding should filter sample0, kept {kept}"
    );
    assert!(
        knee_dir
            .join("sample_sample0")
            .join("permit_freq.bin")
            .exists()
    );
    for name in ["sample1", "sample2"] {
        let s = sample(&info, name);
        assert_eq!(s["num_cells"], 0, "{name}");
        assert_eq!(s["num_reads"], 0, "{name}");
        assert!(
            !knee_dir
                .join(format!("sample_{name}"))
                .join("permit_freq.bin")
                .exists(),
            "{name} should have no permit files"
        );
    }

    let force_dir = tmp.path().join("force");
    generate_permit_list(multi_gpl_opts(
        &rad_dir,
        &force_dir,
        CellFilterMethod::ForceCells(2),
        sample_list,
        &log,
    ))
    .unwrap();
    let info = read_info(&force_dir);
    assert_eq!(sample(&info, "sample1")["num_cells"], 1);
    assert_eq!(sample(&info, "sample2")["num_cells"], 2);
}

/// Single-sample knee finding on one or two distinct barcodes is an error
/// rather than an assertion failure; with none, it keeps none as before.
#[test]
fn single_sample_knee_with_too_few_barcodes_is_an_error() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();
    for cells in [vec![10], vec![10, 4]] {
        let name = format!("cells{}", cells.len());
        let rad_dir = tmp.path().join(&name);
        std::fs::create_dir_all(&rad_dir).unwrap();
        write_single_bc_rad(&rad_dir.join("map.rad"), &cells);
        let out_dir = tmp.path().join(format!("{name}_out"));
        let err = generate_permit_list(single_gpl_opts(
            &rad_dir,
            &out_dir,
            CellFilterMethod::KneeFinding,
            &log,
        ))
        .expect_err(&name);
        assert!(
            chain(&err).contains("knee finding needs at least 3 distinct cell barcodes"),
            "{name}: {err:#}"
        );
    }
}

/// `view` prints only the single-barcode layout; on a multi-barcode RAD it
/// returns an error that says so, and a too-new RAD is an error too.
#[test]
fn view_rejects_multi_barcode_and_unreadable_rad() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();
    let multi = tmp.path().join("multi.rad");
    write_multi_bc_rad(&multi, &[vec![3, 2]]);
    let err = alevin_fry::convert::view(&multi, false, &log).expect_err("multi-barcode view");
    assert!(chain(&err).contains("multi-barcode"), "{err:#}");

    let single = tmp.path().join("single.rad");
    write_single_bc_rad(&single, &[2]);
    let mut bytes = std::fs::read(&single).unwrap();
    bytes[8] = 3;
    std::fs::write(&single, bytes).unwrap();
    let err = alevin_fry::convert::view(&single, false, &log).expect_err("too-new view");
    assert!(
        chain(&err).contains("RAD spec major version 3 is newer than supported"),
        "{err:#}"
    );
}

/// Errors from `collate` and `view` leave `main` as an error (exit status 1)
/// instead of an `expect`/`unwrap` panic. A multi-sample permit list with no
/// cells is the real-world trigger for collate: a sample barcode list that
/// matches nothing.
#[test]
fn binary_exits_with_an_error_for_zero_cell_collate_and_multi_barcode_view() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();
    let rad_dir = tmp.path().join("rad");
    std::fs::create_dir_all(&rad_dir).unwrap();
    write_multi_bc_rad(&rad_dir.join("map.rad"), &[vec![5, 5, 5]]);
    // A sample list whose only entry is a barcode no read carries.
    let sample_list = tmp.path().join("samples.tsv");
    std::fs::write(&sample_list, "AAAAAAAA\tnone\n").unwrap();
    let gpl_dir = tmp.path().join("gpl");
    let cells = generate_permit_list(multi_gpl_opts(
        &rad_dir,
        &gpl_dir,
        CellFilterMethod::ForceCells(3),
        sample_list,
        &log,
    ))
    .unwrap();
    assert_eq!(cells, 0);

    let af = env!("CARGO_BIN_EXE_alevin-fry");
    let output = Command::new(af)
        .arg("collate")
        .arg("-i")
        .arg(&gpl_dir)
        .arg("-r")
        .arg(&rad_dir)
        .args(["-t", "2"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("could not collate"), "{stderr}");
    assert!(!stderr.contains("panicked"), "{stderr}");

    let output = Command::new(af)
        .arg("view")
        .arg("-r")
        .arg(rad_dir.join("map.rad"))
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("multi-barcode"), "{stderr}");
    assert!(!stderr.contains("panicked"), "{stderr}");
}
