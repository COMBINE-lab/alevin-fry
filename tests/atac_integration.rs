//! Integration tests for the scATAC `generate-permit-list -> collate -> deduplicate`
//! path.
//!
//! This path had no coverage at all, which is how two separate defects survived
//! from the commit that introduced the ATAC module until they were found by
//! hand:
//!
//!   * `collate` aborted on any input containing unmapped or multi-mapping
//!     reads, because each temp bucket's record budget came from the permit
//!     list (which counts every read) while the scatter phase skips those
//!     records.
//!   * `deduplicate` never started its producer thread, so its consumers span
//!     forever on a queue nothing would ever fill.
//!
//! Both are cheap to re-introduce and expensive to notice, because `simpleaf`
//! drives scATAC through `atac sort` and never calls either command. The
//! fixtures below therefore deliberately include unmapped and multi-mapping
//! records — a fixture of cleanly, uniquely mapped reads passes even against
//! the unfixed code.

use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::time::Duration;

use libradicl::chunk::Chunk;
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{
    RadIntId, RadType, TagDesc, TagMap, TagSection, TagSectionLabel, TagValue,
};
use libradicl::record::{AtacSeqReadRecord, AtacSeqRecordContext, RecordContext};
use libradicl::writers::RadFileWriter;

use alevin_fry::atac::cellfilter::{CellFilterMethod, generate_permit_list};
use alevin_fry::atac::collate::collate;
use alevin_fry::atac::deduplicate::deduplicate;
use alevin_fry::atac::prog_opts::{DeduplicateOpts, GenPermitListOpts};
use alevin_fry::atac::sort::sort;
use alevin_fry::barcode_correction::{BarcodeNeighborhood, BarcodeResolution};
use alevin_fry::correction_plan::CorrectionPlan;

const TEST_VERSION: &str = "0.17.0";
const CELL_BC_LEN: u16 = 16;
const REF_NAMES: [&str; 2] = ["chr1", "chr2"];
const REF_LENGTHS: [u32; 2] = [1_000_000, 800_000];

/// A hard ceiling on any single pipeline stage.
///
/// `deduplicate`'s failure mode was an infinite spin, not a crash, so a plain
/// call would hang the whole test binary rather than fail it. Every stage runs
/// under this deadline so a regression reports itself.
const STAGE_TIMEOUT: Duration = Duration::from_secs(120);

/// Pack `idx` into a 2-bit-per-base barcode of length `len`.
///
/// Mixed so consecutive cells do not differ by a single base — otherwise the
/// permit list's 1-mismatch recovery would treat them as neighbours.
fn make_packed_bc(idx: u64, len: u16) -> u64 {
    let mask = (1u64 << (2 * len as u64)) - 1;
    idx.wrapping_mul(2654435761) & mask
}

/// Render a packed barcode as nucleotides.
///
/// The first base occupies the *most* significant bit pair, matching
/// `needletail::bitkmer::BitNuclKmer`, which is what reads the whitelist back
/// in. Getting this backwards makes every barcode silently fail to match.
fn packed_to_nuc(packed: u64, len: usize) -> String {
    const NUCS: [char; 4] = ['A', 'C', 'G', 'T'];
    (0..len)
        .map(|i| NUCS[((packed >> (2 * (len - 1 - i))) & 3) as usize])
        .collect()
}

/// Build the prelude piscem-rs emits for scATAC: paired flag set, the `cblen`
/// / `known_rad_type` / `ref_lengths` file tags, a `b` read tag, and the four
/// alignment tags.
fn make_atac_prelude() -> (RadPrelude, TagMap) {
    let hdr = RadHeader {
        is_paired: 1,
        ref_count: REF_NAMES.len() as u64,
        ref_names: REF_NAMES.iter().map(|s| s.to_string()).collect(),
        num_chunks: 0, // backpatched by RadFileWriter
    };

    let mut file_tags = TagSection::new_with_label(TagSectionLabel::FileTags);
    file_tags.add_tag_desc(TagDesc {
        name: "cblen".to_string(),
        typeid: RadType::Int(RadIntId::U16),
    });
    file_tags.add_tag_desc(TagDesc {
        name: "known_rad_type".to_string(),
        typeid: RadType::String,
    });
    file_tags.add_tag_desc(TagDesc {
        name: "ref_lengths".to_string(),
        typeid: RadType::Array(
            RadIntId::U32,
            libradicl::rad_types::RadAtomicId::Int(RadIntId::U32),
        ),
    });

    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    read_tags.add_tag_desc(TagDesc {
        name: "b".to_string(),
        typeid: RadType::Int(RadIntId::U32),
    });

    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    for (name, typeid) in [
        ("ref", RadType::Int(RadIntId::U32)),
        ("type", RadType::Int(RadIntId::U8)),
        ("start_pos", RadType::Int(RadIntId::U32)),
        ("frag_len", RadType::Int(RadIntId::U16)),
    ] {
        aln_tags.add_tag_desc(TagDesc {
            name: name.to_string(),
            typeid,
        });
    }

    let prelude = RadPrelude {
        hdr,
        file_tags,
        read_tags,
        aln_tags,
    };

    let mut file_tag_map = TagMap::with_keyset(&prelude.file_tags.tags);
    file_tag_map.add(TagValue::U16(CELL_BC_LEN));
    file_tag_map.add(TagValue::String("sc_atac".to_string()));
    file_tag_map.add(TagValue::ArrayU32(REF_LENGTHS.to_vec()));

    (prelude, file_tag_map)
}

/// How many records of each kind each synthetic cell contributes.
struct CellMix {
    /// Uniquely mapped, properly paired — the only kind `collate` keeps.
    good: usize,
    /// No alignments at all. `AtacSeqReadRecord::is_empty` is true, so the
    /// scatter phase drops these while the permit list still counted them.
    unmapped: usize,
    /// Two alignments. Dropped by the scatter phase for `naln > 1`.
    multimapped: usize,
}

impl CellMix {
    fn total(&self) -> usize {
        self.good + self.unmapped + self.multimapped
    }
}

/// Write a synthetic scATAC RAD file, one chunk per cell.
///
/// Returns the packed cell barcodes in the order they were written.
fn create_synthetic_atac_rad(path: &Path, num_cells: usize, mix: &CellMix) -> Vec<u64> {
    let (prelude, file_tag_map) = make_atac_prelude();
    let ctx = AtacSeqRecordContext::get_context_from_tag_section(
        &prelude.file_tags,
        &prelude.read_tags,
        &prelude.aln_tags,
    )
    .expect("could not build ATAC record context");

    let file = File::create(path).expect("could not create synthetic RAD");
    let mut fw =
        RadFileWriter::new(file, &prelude, &file_tag_map).expect("could not write prelude");

    let mut barcodes = Vec::with_capacity(num_cells);

    for cell_idx in 0..num_cells {
        let bc = make_packed_bc(cell_idx as u64 + 1, CELL_BC_LEN);
        barcodes.push(bc);

        let mut reads = Vec::with_capacity(mix.total());

        for r in 0..mix.good {
            let ref_id = (r % REF_NAMES.len()) as u32;
            // Spread positions out so fragments do not all collapse into one
            // another during deduplication.
            let start = 1_000 + (cell_idx as u32 * 977) + (r as u32 * 131);
            reads.push(AtacSeqReadRecord {
                bc,
                start_pos: vec![start],
                refs: vec![ref_id],
                frag_lengths: vec![120],
                // 4 is the "properly mapped pair" type that `deduplicate` keeps.
                map_type: vec![4],
            });
        }

        for _ in 0..mix.unmapped {
            reads.push(AtacSeqReadRecord {
                bc,
                start_pos: vec![],
                refs: vec![],
                frag_lengths: vec![],
                map_type: vec![],
            });
        }

        for r in 0..mix.multimapped {
            let start = 500_000 + (cell_idx as u32 * 13) + (r as u32 * 7);
            reads.push(AtacSeqReadRecord {
                bc,
                start_pos: vec![start, start + 250],
                refs: vec![0, 1],
                frag_lengths: vec![110, 115],
                map_type: vec![4, 4],
            });
        }

        let chunk = Chunk::<AtacSeqReadRecord> {
            nbytes: 0, // computed by write_chunk
            nrec: reads.len() as u32,
            reads,
        };
        fw.write_chunk(&chunk, &ctx).expect("could not write chunk");
    }

    fw.finalize().expect("could not finalize synthetic RAD");
    barcodes
}

/// Write the barcode whitelist `generate-permit-list` consumes.
fn write_whitelist(path: &Path, barcodes: &[u64]) {
    let mut f = File::create(path).expect("could not create whitelist");
    for bc in barcodes {
        writeln!(f, "{}", packed_to_nuc(*bc, CELL_BC_LEN as usize))
            .expect("could not write whitelist entry");
    }
}

fn make_test_logger() -> slog::Logger {
    use slog::Drain;
    let decorator = slog_term::PlainSyncDecorator::new(slog_term::TestStdoutWriter);
    let drain = slog_term::FullFormat::new(decorator).build().fuse();
    slog::Logger::root(drain, slog::o!())
}

/// Run `f` on its own thread and fail if it outlives [`STAGE_TIMEOUT`].
///
/// A hung stage leaks its thread; that is deliberate. The alternative is
/// blocking the test binary forever, which is exactly how the `deduplicate`
/// spin went unnoticed.
fn run_with_timeout<F>(what: &str, f: F)
where
    F: FnOnce() -> anyhow::Result<()> + Send + 'static,
{
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(f());
    });
    match rx.recv_timeout(STAGE_TIMEOUT) {
        Ok(Ok(())) => {}
        Ok(Err(e)) => panic!("{what} failed: {e:?}"),
        Err(mpsc::RecvTimeoutError::Timeout) => {
            panic!("{what} did not finish within {STAGE_TIMEOUT:?}; it is most likely spinning")
        }
        Err(mpsc::RecvTimeoutError::Disconnected) => panic!("{what} panicked"),
    }
}

/// Build a fixture and run generate-permit-list, returning (rad_dir, gpl_dir).
fn stage_permit_list(
    tmp: &Path,
    num_cells: usize,
    mix: &CellMix,
    log: &slog::Logger,
) -> (PathBuf, PathBuf) {
    let rad_dir = tmp.join("rad");
    let gpl_dir = tmp.join("gpl");
    std::fs::create_dir_all(&rad_dir).unwrap();
    std::fs::create_dir_all(&gpl_dir).unwrap();

    let barcodes = create_synthetic_atac_rad(&rad_dir.join("map.rad"), num_cells, mix);
    let wl_path = tmp.join("whitelist.txt");
    write_whitelist(&wl_path, &barcodes);

    // piscem always emits this alongside map.rad, and `collate` opens it
    // unconditionally. A zero-length file is the legitimate "no unmapped
    // barcodes" case: it is read as a stream of (u64 barcode, u32 count) pairs.
    File::create(rad_dir.join("unmapped_bc_count.bin")).expect("could not create unmapped counts");

    let gpl_opts = GenPermitListOpts::builder()
        .input_dir(&rad_dir)
        .output_dir(&gpl_dir)
        .fmeth(CellFilterMethod::UnfilteredExternalList(wl_path.clone(), 1))
        .threads(2)
        .rc(false)
        .cmdline("atac_integration_test")
        .version(TEST_VERSION)
        .log(log)
        .build();

    // The return value is the number of barcodes *rescued* by 1-mismatch
    // correction, not the number retained. Every synthetic barcode is in the
    // whitelist verbatim, so nothing needs rescuing.
    let num_corrected = generate_permit_list(gpl_opts).expect("atac generate-permit-list failed");
    assert_eq!(
        num_corrected, 0,
        "synthetic barcodes match the whitelist exactly, so none should need correction"
    );

    let correction_plan = CorrectionPlan::read_from(&gpl_dir.join("correction_plan.bin"))
        .expect("ATAC GPL wrote no valid correction plan");
    assert_eq!(correction_plan.cell_barcode_len, CELL_BC_LEN as u8);
    assert_eq!(correction_plan.cell_scopes.len(), 1);
    assert_eq!(
        correction_plan.cell_scopes[0].spec.neighborhood,
        BarcodeNeighborhood::HammingOne
    );
    assert_eq!(
        correction_plan.cell_scopes[0].spec.resolution,
        BarcodeResolution::Unique
    );
    assert_eq!(
        correction_plan.cell_scopes[0].corrections.len(),
        num_cells,
        "the compact ATAC plan must include retained identities"
    );

    // What we actually care about is the retained set, which lands in
    // permit_freq.bin. That file opens with two u64s — the format version and
    // the barcode length — before the bincode-encoded barcode -> count map.
    let freq_file = File::open(gpl_dir.join("permit_freq.bin"))
        .expect("generate-permit-list wrote no permit_freq.bin");
    let mut freq_reader = std::io::BufReader::new(freq_file);
    let mut header = [0u8; 16];
    freq_reader
        .read_exact(&mut header)
        .expect("permit_freq.bin is truncated");
    assert_eq!(
        u64::from_le_bytes(header[8..16].try_into().unwrap()),
        CELL_BC_LEN as u64,
        "permit_freq.bin records an unexpected barcode length"
    );
    let freq: std::collections::HashMap<u64, u64> =
        bincode::deserialize_from(&mut freq_reader).expect("could not read permit_freq.bin");
    assert_eq!(
        freq.len(),
        num_cells,
        "every synthetic barcode is in the whitelist, so all should be retained"
    );
    for count in freq.values() {
        assert_eq!(
            *count,
            mix.total() as u64,
            "the permit list counts every read of a retained barcode, including the \
             unmapped and multi-mapping ones that collate will later drop"
        );
    }

    (rad_dir, gpl_dir)
}

/// `collate` used to abort whenever the input held records the scatter phase
/// drops, because the temp-bucket record budget came from the permit list.
///
/// The budget is the loop bound the gather phase reads each bucket with, so
/// getting it wrong is not merely a failed assertion.
#[test]
fn atac_collate_tolerates_unmapped_and_multimapping_records() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();

    // Every cell contributes records of all three kinds, so the permit-list
    // count and the number of records actually scattered cannot agree.
    let mix = CellMix {
        good: 6,
        unmapped: 3,
        multimapped: 2,
    };
    let num_cells = 24;

    let (rad_dir, gpl_dir) = stage_permit_list(tmp.path(), num_cells, &mix, &log);

    // Prove this run consumes the compiled artifact rather than silently
    // relying on the compatibility map.
    std::fs::remove_file(gpl_dir.join("permit_map.bin")).unwrap();

    let gpl_for_collate = gpl_dir.clone();
    let rad_for_collate = rad_dir.clone();
    let collate_log = log.clone();
    run_with_timeout("atac collate", move || {
        collate(
            gpl_for_collate,
            rad_for_collate,
            2,
            10_000,
            false,
            "atac_integration_test",
            TEST_VERSION,
            &collate_log,
        )
    });

    let collated = gpl_dir.join("map.collated.rad");
    let len = std::fs::metadata(&collated)
        .expect("collate produced no output file")
        .len();
    assert!(
        len > 0,
        "collated RAD is empty; collate reported success without writing records"
    );

    // One chunk per cell: no cell here loses *all* of its records.
    let f = File::open(&collated).unwrap();
    let mut br = std::io::BufReader::new(f);
    let prelude = RadPrelude::from_bytes(&mut br).expect("could not parse collated prelude");
    assert_eq!(
        prelude.hdr.num_chunks, num_cells as u64,
        "expected one collated chunk per cell"
    );
}

#[test]
fn atac_sort_uses_compiled_plan_without_legacy_map() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();
    let mix = CellMix {
        good: 3,
        unmapped: 0,
        multimapped: 0,
    };
    let (rad_dir, gpl_dir) = stage_permit_list(tmp.path(), 4, &mix, &log);
    std::fs::remove_file(gpl_dir.join("permit_map.bin")).unwrap();

    let sort_input = gpl_dir.clone();
    let sort_rad = rad_dir.clone();
    let sort_log = log.clone();
    run_with_timeout("atac sort with compiled correction plan", move || {
        sort(
            sort_input,
            sort_rad,
            2,
            10_000,
            false,
            "atac_integration_test",
            TEST_VERSION,
            &sort_log,
        )
    });
    assert!(
        std::fs::metadata(gpl_dir.join("map.bed")).is_ok_and(|metadata| metadata.len() > 0),
        "ATAC sort produced no BED output"
    );
}

#[test]
fn atac_sort_uses_legacy_map_without_compiled_plan() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();
    let mix = CellMix {
        good: 3,
        unmapped: 0,
        multimapped: 0,
    };
    let (rad_dir, gpl_dir) = stage_permit_list(tmp.path(), 4, &mix, &log);
    std::fs::remove_file(gpl_dir.join("correction_plan.bin")).unwrap();

    let sort_input = gpl_dir.clone();
    let sort_rad = rad_dir.clone();
    let sort_log = log.clone();
    run_with_timeout("atac sort with legacy permit map", move || {
        sort(
            sort_input,
            sort_rad,
            2,
            10_000,
            false,
            "atac_integration_test",
            TEST_VERSION,
            &sort_log,
        )
    });
    assert!(
        std::fs::metadata(gpl_dir.join("map.bed")).is_ok_and(|metadata| metadata.len() > 0),
        "ATAC legacy fallback produced no BED output"
    );
}

/// A cell whose records are *all* dropped by the scatter phase contributes no
/// chunk, so the collated file legitimately holds fewer chunks than there are
/// retained barcodes. That used to be a hard assertion failure.
#[test]
fn atac_collate_handles_cells_with_no_surviving_records() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();

    // Not one usable record anywhere: every read is unmapped or multi-mapping.
    let mix = CellMix {
        good: 0,
        unmapped: 4,
        multimapped: 2,
    };
    let num_cells = 8;

    let (rad_dir, gpl_dir) = stage_permit_list(tmp.path(), num_cells, &mix, &log);

    let gpl_for_collate = gpl_dir.clone();
    let collate_log = log.clone();
    run_with_timeout("atac collate (all records dropped)", move || {
        collate(
            gpl_for_collate,
            rad_dir,
            2,
            10_000,
            false,
            "atac_integration_test",
            TEST_VERSION,
            &collate_log,
        )
    });

    let f = File::open(gpl_dir.join("map.collated.rad")).unwrap();
    let mut br = std::io::BufReader::new(f);
    let prelude = RadPrelude::from_bytes(&mut br).expect("could not parse collated prelude");
    assert_eq!(
        prelude.hdr.num_chunks, 0,
        "no record survives the scatter phase, so no chunk should be written"
    );
}

/// `deduplicate` built a reader, spawned consumers, and never started the
/// producer, so it span forever. This drives the whole path and asserts the
/// BED actually describes the fragments we put in.
#[test]
fn atac_deduplicate_terminates_and_writes_fragments() {
    let tmp = tempfile::tempdir().unwrap();
    let log = make_test_logger();

    let mix = CellMix {
        good: 5,
        unmapped: 2,
        multimapped: 1,
    };
    let num_cells = 16;

    let (rad_dir, gpl_dir) = stage_permit_list(tmp.path(), num_cells, &mix, &log);

    let gpl_for_collate = gpl_dir.clone();
    let collate_log = log.clone();
    run_with_timeout("atac collate", move || {
        collate(
            gpl_for_collate,
            rad_dir,
            2,
            10_000,
            false,
            "atac_integration_test",
            TEST_VERSION,
            &collate_log,
        )
    });

    let gpl_for_dedup = gpl_dir.clone();
    let dedup_log = log.clone();
    run_with_timeout("atac deduplicate", move || {
        let opts = DeduplicateOpts::builder()
            .input_dir(&gpl_for_dedup)
            .num_threads(4)
            .rev(false)
            .cmdline("atac_integration_test")
            .version(TEST_VERSION)
            .log(&dedup_log)
            .build();
        deduplicate(opts)
    });

    let bed = std::fs::read_to_string(gpl_dir.join("map.bed")).expect("deduplicate wrote no BED");
    let lines: Vec<&str> = bed.lines().collect();
    assert_eq!(
        lines.len(),
        num_cells * mix.good,
        "every uniquely-mapped fragment should appear exactly once"
    );

    let valid_refs: Vec<String> = REF_NAMES.iter().map(|s| s.to_string()).collect();
    for line in &lines {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 5, "malformed BED line: {line:?}");
        assert!(valid_refs.contains(&f[0].to_string()), "bad ref: {line:?}");
        let start: u64 = f[1].parse().expect("non-numeric start");
        let end: u64 = f[2].parse().expect("non-numeric end");
        assert!(end > start, "empty or inverted interval: {line:?}");
        assert_eq!(
            f[3].len(),
            CELL_BC_LEN as usize,
            "barcode has the wrong length: {line:?}"
        );
        let count: u64 = f[4].parse().expect("non-numeric count");
        assert!(count >= 1, "non-positive count: {line:?}");
    }

    // The multi-mapping records must not reach the BED: the scatter phase
    // drops them, and they were placed far from the uniquely-mapped ones.
    assert!(
        lines.iter().all(|l| {
            let start: u64 = l.split('\t').nth(1).unwrap().parse().unwrap();
            start < 500_000
        }),
        "a multi-mapping record leaked into the deduplicated output"
    );
}
