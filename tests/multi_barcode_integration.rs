//! Integration tests for the multi-barcode (10x Flex) pipeline.
//!
//! These tests create synthetic multi-barcode RAD files programmatically,
//! then exercise the generate-permit-list → collate → quant pipeline.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

use libradicl::chunk::Chunk;
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{
    RadIntId, RadType, TagDesc, TagMap, TagSection, TagSectionLabel, TagValue,
};
use libradicl::record::{MultiBarcodeReadRecord, MultiBarcodeRecordContext, RecordContext};
use libradicl::writers::RadFileWriter;
use smallvec::smallvec;

// -----------------------------------------------------------------------
// Synthetic RAD file builder
// -----------------------------------------------------------------------

/// Number of reference targets in test data
const NUM_REFS: u64 = 10;
const TEST_VERSION: &str = "0.12.0";
/// Sample barcode length (bases)
const SAMPLE_BC_LEN: u16 = 8;
/// Cell barcode length (bases)
const CELL_BC_LEN: u16 = 16;
/// UMI length (bases)
const UMI_LEN: u16 = 12;

/// Simple deterministic barcode generator from an index.
fn make_packed_bc(idx: u64, len: u16) -> u64 {
    // Use low bits, masked to the barcode length
    let mask = (1u64 << (2 * len as u64)) - 1;
    // Mix the index a bit so different barcodes don't collide trivially
    idx.wrapping_mul(2654435761) & mask
}

/// Build a multi-barcode RadPrelude and TagMap for test RAD files.
fn make_multi_bc_prelude() -> (RadPrelude, TagMap) {
    let mut ref_names = Vec::new();
    for i in 0..NUM_REFS {
        ref_names.push(format!("gene_{}", i));
    }

    let hdr = RadHeader {
        is_paired: 0,
        ref_count: NUM_REFS,
        ref_names,
        num_chunks: 0, // backpatched by RadFileWriter
    };

    // File-level tags: num_barcodes, b0len, b1len, ulen, known_rad_type
    let mut file_tags = TagSection::new_with_label(TagSectionLabel::FileTags);
    file_tags.add_tag_desc(TagDesc {
        name: "num_barcodes".to_string(),
        typeid: RadType::Int(RadIntId::U16),
    });
    file_tags.add_tag_desc(TagDesc {
        name: "b0len".to_string(),
        typeid: RadType::Int(RadIntId::U16),
    });
    file_tags.add_tag_desc(TagDesc {
        name: "b1len".to_string(),
        typeid: RadType::Int(RadIntId::U16),
    });
    file_tags.add_tag_desc(TagDesc {
        name: "ulen".to_string(),
        typeid: RadType::Int(RadIntId::U16),
    });
    file_tags.add_tag_desc(TagDesc {
        name: "known_rad_type".to_string(),
        typeid: RadType::String,
    });

    // Read-level tags: b0 (sample), b1 (cell), u (UMI)
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    read_tags.add_tag_desc(TagDesc {
        name: "b0".to_string(),
        typeid: RadType::Int(RadIntId::U32),
    });
    read_tags.add_tag_desc(TagDesc {
        name: "b1".to_string(),
        typeid: RadType::Int(RadIntId::U32),
    });
    read_tags.add_tag_desc(TagDesc {
        name: "u".to_string(),
        typeid: RadType::Int(RadIntId::U32),
    });

    // Alignment-level tags
    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    aln_tags.add_tag_desc(TagDesc {
        name: "compressed_ori_refid".to_string(),
        typeid: RadType::Int(RadIntId::U32),
    });

    let prelude = RadPrelude {
        hdr,
        file_tags,
        read_tags,
        aln_tags,
    };

    // File-level tag values
    let mut file_tag_map = TagMap::with_keyset(&prelude.file_tags.tags);
    file_tag_map.add(TagValue::U16(2)); // num_barcodes
    file_tag_map.add(TagValue::U16(SAMPLE_BC_LEN)); // b0len
    file_tag_map.add(TagValue::U16(CELL_BC_LEN)); // b1len
    file_tag_map.add(TagValue::U16(UMI_LEN)); // ulen
    file_tag_map.add(TagValue::String("sc_rna_multi_bc".to_string())); // known_rad_type

    (prelude, file_tag_map)
}

/// Create a synthetic multi-barcode RAD file.
///
/// Creates `num_samples` samples, each with `cells_per_sample` cells,
/// each with `reads_per_cell` reads. All reads map to a single reference.
fn create_synthetic_multi_bc_rad(
    path: &Path,
    num_samples: usize,
    cells_per_sample: usize,
    reads_per_cell: usize,
    sample_barcodes: &[u64],
) -> anyhow::Result<MultiBarcodeRecordContext> {
    create_synthetic_multi_bc_rad_with_shared_cells(
        path,
        num_samples,
        cells_per_sample,
        reads_per_cell,
        sample_barcodes,
        false,
    )
}

/// Create a synthetic multi-barcode RAD file with optional shared cell barcodes across samples.
fn create_synthetic_multi_bc_rad_with_shared_cells(
    path: &Path,
    num_samples: usize,
    cells_per_sample: usize,
    reads_per_cell: usize,
    sample_barcodes: &[u64],
    share_cell_barcodes_across_samples: bool,
) -> anyhow::Result<MultiBarcodeRecordContext> {
    let (prelude, file_tag_map) = make_multi_bc_prelude();

    let ctx = MultiBarcodeRecordContext::get_context_from_tag_section(
        &prelude.file_tags,
        &prelude.read_tags,
        &prelude.aln_tags,
    )?;

    let file = File::create(path)?;
    let mut fw = RadFileWriter::new(file, &prelude, &file_tag_map)?;

    // Generate chunks. Each chunk contains all reads for one cell.
    // We interleave samples so records are NOT pre-sorted by sample.
    for cell_idx in 0..cells_per_sample {
        for (sample_idx, sample_bc) in sample_barcodes
            .iter()
            .copied()
            .enumerate()
            .take(num_samples)
        {
            let cell_seed = if share_cell_barcodes_across_samples {
                cell_idx as u64
            } else {
                (sample_idx * 1000 + cell_idx) as u64
            };
            let cell_bc = make_packed_bc(cell_seed, CELL_BC_LEN);

            let mut reads = Vec::with_capacity(reads_per_cell);
            for read_idx in 0..reads_per_cell {
                let umi = make_packed_bc(
                    (sample_idx * 100000 + cell_idx * 100 + read_idx) as u64,
                    UMI_LEN,
                );
                let ref_id = (read_idx % NUM_REFS as usize) as u32;
                reads.push(MultiBarcodeReadRecord {
                    barcodes: smallvec![sample_bc, cell_bc],
                    umi,
                    dirs: vec![true],
                    refs: vec![ref_id],
                });
            }

            let chunk = Chunk::<MultiBarcodeReadRecord> {
                nbytes: 0, // computed by write_chunk
                nrec: reads.len() as u32,
                reads,
            };
            fw.write_chunk(&chunk, &ctx)?;
        }
    }

    fw.finalize()?;
    Ok(ctx)
}

/// Write a sample barcode list file (one nucleotide sequence per line).
fn write_sample_bc_list(path: &Path, barcodes: &[u64], bc_len: u16) -> anyhow::Result<()> {
    let mut f = BufWriter::new(File::create(path)?);
    for &bc in barcodes {
        // Convert packed barcode back to nucleotide string
        let seq = packed_to_nuc(bc, bc_len as usize);
        writeln!(f, "{}", seq)?;
    }
    Ok(())
}

fn write_named_sample_bc_list(
    path: &Path,
    entries: &[(&str, u64)],
    bc_len: u16,
) -> anyhow::Result<()> {
    let mut f = BufWriter::new(File::create(path)?);
    for (name, bc) in entries {
        let seq = packed_to_nuc(*bc, bc_len as usize);
        writeln!(f, "{}\t{}", seq, name)?;
    }
    Ok(())
}

fn write_tg_map(path: &Path) -> anyhow::Result<()> {
    let mut f = BufWriter::new(File::create(path)?);
    for i in 0..NUM_REFS {
        writeln!(f, "gene_{}\tgene_{}", i, i)?;
    }
    Ok(())
}

fn make_test_logger() -> slog::Logger {
    let decorator = slog_term::PlainDecorator::new(slog_term::TestStdoutWriter);
    let drain = slog_term::CompactFormat::new(decorator).build();
    let drain = std::sync::Mutex::new(drain);
    let drain = slog::Fuse::new(drain);
    slog::Logger::root(drain, slog::o!())
}

/// Convert a 2-bit packed barcode to a nucleotide string.
fn packed_to_nuc(packed: u64, len: usize) -> String {
    let bases = ['A', 'C', 'G', 'T'];
    let mut s = String::with_capacity(len);
    for i in 0..len {
        let shift = 2 * (len - 1 - i);
        let code = ((packed >> shift) & 3) as usize;
        s.push(bases[code]);
    }
    s
}

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

/// Test that we can create and read back a synthetic multi-barcode RAD file.
#[test]
fn test_synthetic_rad_roundtrip() {
    let tmp = tempfile::tempdir().unwrap();
    let rad_path = tmp.path().join("map.rad");

    let sample_bcs = vec![
        make_packed_bc(100, SAMPLE_BC_LEN),
        make_packed_bc(200, SAMPLE_BC_LEN),
    ];
    let num_samples = 2;
    let cells_per_sample = 3;
    let reads_per_cell = 5;

    let ctx = create_synthetic_multi_bc_rad(
        &rad_path,
        num_samples,
        cells_per_sample,
        reads_per_cell,
        &sample_bcs,
    )
    .unwrap();

    // Read it back
    let file = File::open(&rad_path).unwrap();
    let mut reader = BufReader::new(file);
    let prelude = RadPrelude::from_bytes(&mut reader).unwrap();
    let file_tag_map = prelude
        .file_tags
        .parse_tags_from_bytes(&mut reader)
        .unwrap();

    // Verify header
    let expected_chunks = (num_samples * cells_per_sample) as u64;
    assert_eq!(prelude.hdr.num_chunks, expected_chunks);
    assert_eq!(prelude.hdr.ref_count, NUM_REFS);

    // Verify file tags
    let num_bc: u16 = file_tag_map
        .get("num_barcodes")
        .unwrap()
        .try_into()
        .unwrap();
    assert_eq!(num_bc, 2);
    let b0len: u16 = file_tag_map.get("b0len").unwrap().try_into().unwrap();
    assert_eq!(b0len, SAMPLE_BC_LEN);
    let b1len: u16 = file_tag_map.get("b1len").unwrap().try_into().unwrap();
    assert_eq!(b1len, CELL_BC_LEN);

    // Read all chunks and verify record counts
    let mut total_reads = 0u64;
    let mut sample_bc_counts: HashMap<u64, u64> = HashMap::new();
    for _ in 0..expected_chunks {
        let c = Chunk::<MultiBarcodeReadRecord>::from_bytes(&mut reader, &ctx);
        assert_eq!(c.nrec, reads_per_cell as u32);
        for read in &c.reads {
            total_reads += 1;
            let sample_bc = read.barcodes[0];
            *sample_bc_counts.entry(sample_bc).or_insert(0) += 1;
        }
    }

    assert_eq!(
        total_reads,
        (num_samples * cells_per_sample * reads_per_cell) as u64,
    );
    // Each sample should have cells_per_sample * reads_per_cell reads
    for &bc in &sample_bcs {
        assert_eq!(
            *sample_bc_counts.get(&bc).unwrap_or(&0),
            (cells_per_sample * reads_per_cell) as u64,
        );
    }
}

/// Test the generate-permit-list step for multi-barcode data.
#[test]
fn test_multi_bc_generate_permit_list() {
    use alevin_fry::cellfilter::{CellFilterMethod, generate_permit_list};
    use alevin_fry::prog_opts::{GenPermitListOpts, SampleCorrectionMode};
    use bio_types::strand::Strand;

    let tmp = tempfile::tempdir().unwrap();
    let rad_dir = tmp.path().join("rad");
    std::fs::create_dir_all(&rad_dir).unwrap();
    let output_dir = tmp.path().join("output");
    std::fs::create_dir_all(&output_dir).unwrap();

    let sample_bcs = vec![
        make_packed_bc(100, SAMPLE_BC_LEN),
        make_packed_bc(200, SAMPLE_BC_LEN),
    ];
    let num_samples = 2;
    let cells_per_sample = 5;
    let reads_per_cell = 10;

    // Create RAD file
    let rad_path = rad_dir.join("map.rad");
    create_synthetic_multi_bc_rad(
        &rad_path,
        num_samples,
        cells_per_sample,
        reads_per_cell,
        &sample_bcs,
    )
    .unwrap();

    // Write sample barcode list
    let sample_list_path = tmp.path().join("sample_barcodes.txt");
    write_sample_bc_list(&sample_list_path, &sample_bcs, SAMPLE_BC_LEN).unwrap();

    // Setup logger
    let log = make_test_logger();

    // Run generate-permit-list
    let gpl_opts = GenPermitListOpts::builder()
        .input_dir(&rad_dir)
        .output_dir(&output_dir)
        .fmeth(CellFilterMethod::ForceCells(cells_per_sample))
        .expected_ori(Strand::Unknown)
        .version("test")
        .threads(1)
        .velo_mode(false)
        .cmdline("test")
        .log(&log)
        .sample_bc_list(Some(sample_list_path))
        .sample_names(None)
        .sample_correction_mode(SampleCorrectionMode::Exact)
        .build();

    let total_cells = generate_permit_list(gpl_opts).unwrap();
    assert!(
        total_cells > 0,
        "expected at least some cells to be retained"
    );

    // Verify outputs
    assert!(output_dir.join("sample_permit_map.bin").exists());
    assert!(output_dir.join("sample_info.json").exists());
    assert!(output_dir.join("generate_permit_list.json").exists());

    // Check sample_info.json
    let info_file = File::open(output_dir.join("sample_info.json")).unwrap();
    let info: serde_json::Value = serde_json::from_reader(info_file).unwrap();
    assert_eq!(info["num_samples"].as_u64().unwrap(), num_samples as u64);
    assert_eq!(info["num_barcodes"].as_u64().unwrap(), 2);
    assert!(info["matched_reads"].as_u64().unwrap() > 0);

    // Check per-sample permit maps exist
    let samples = info["samples"].as_array().unwrap();
    for entry in samples {
        let name = entry["name"].as_str().unwrap();
        let sample_dir = output_dir.join(format!("sample_{}", name));
        // At least one sample should have permit maps
        if entry["num_cells"].as_u64().unwrap() > 0 {
            assert!(
                sample_dir.join("permit_map.bin").exists(),
                "permit_map.bin should exist for sample {} with cells",
                name
            );
            assert!(
                sample_dir.join("permit_freq.bin").exists(),
                "permit_freq.bin should exist for sample {} with cells",
                name
            );
        }
    }
}

#[test]
fn test_multi_bc_collate_and_quant_preserve_sample_cell_identity() {
    use alevin_fry::cellfilter::{CellFilterMethod, generate_permit_list};
    use alevin_fry::collate::collate;
    use alevin_fry::prog_opts::{GenPermitListOpts, QuantOpts, SampleCorrectionMode};
    use alevin_fry::quant::{ResolutionStrategy, SplicedAmbiguityModel, quantify};
    use bio_types::strand::Strand;
    use std::collections::HashSet;

    let tmp = tempfile::tempdir().unwrap();
    let rad_dir = tmp.path().join("rad");
    std::fs::create_dir_all(&rad_dir).unwrap();
    let output_dir = tmp.path().join("output");
    std::fs::create_dir_all(&output_dir).unwrap();
    let quant_dir = tmp.path().join("quant");

    let sample_bcs = vec![
        make_packed_bc(100, SAMPLE_BC_LEN),
        make_packed_bc(200, SAMPLE_BC_LEN),
    ];
    let sample_entries = [("sample_a", sample_bcs[0]), ("sample_b", sample_bcs[1])];
    let num_samples = sample_bcs.len();
    let cells_per_sample = 4;
    let reads_per_cell = 8;

    let rad_path = rad_dir.join("map.rad");
    create_synthetic_multi_bc_rad_with_shared_cells(
        &rad_path,
        num_samples,
        cells_per_sample,
        reads_per_cell,
        &sample_bcs,
        true,
    )
    .unwrap();

    let sample_list_path = tmp.path().join("sample_barcodes.tsv");
    write_named_sample_bc_list(&sample_list_path, &sample_entries, SAMPLE_BC_LEN).unwrap();

    let tg_map_path = tmp.path().join("tg_map.tsv");
    write_tg_map(&tg_map_path).unwrap();

    let log = make_test_logger();
    let gpl_opts = GenPermitListOpts::builder()
        .input_dir(&rad_dir)
        .output_dir(&output_dir)
        .fmeth(CellFilterMethod::ForceCells(cells_per_sample))
        .expected_ori(Strand::Unknown)
        .version(TEST_VERSION)
        .threads(2)
        .velo_mode(false)
        .cmdline("test")
        .log(&log)
        .sample_bc_list(Some(sample_list_path))
        .sample_names(None)
        .sample_correction_mode(SampleCorrectionMode::Exact)
        .build();
    let total_cells = generate_permit_list(gpl_opts).unwrap();
    assert_eq!(total_cells, (num_samples * cells_per_sample) as u64);

    collate(
        output_dir.clone(),
        &rad_dir,
        2,
        1_000,
        false,
        "test",
        TEST_VERSION,
        &log,
    )
    .unwrap();

    let collated_path = output_dir.join("map.collated.rad");
    let collated_file = File::open(&collated_path).unwrap();
    let mut collated_reader = BufReader::new(collated_file);
    let collated_prelude = RadPrelude::from_bytes(&mut collated_reader).unwrap();
    let _collated_ftm = collated_prelude
        .file_tags
        .parse_tags_from_bytes(&mut collated_reader)
        .unwrap();
    assert_eq!(
        collated_prelude.hdr.num_chunks,
        (num_samples * cells_per_sample) as u64
    );

    let quant_opts = QuantOpts::builder()
        .input_dir(&output_dir)
        .tg_map(&tg_map_path)
        .output_dir(&quant_dir)
        .num_threads(2)
        .num_bootstraps(0)
        .init_uniform(false)
        .summary_stat(false)
        .dump_eq(false)
        .resolution(ResolutionStrategy::Trivial)
        .pug_exact_umi(false)
        .sa_model(SplicedAmbiguityModel::WinnerTakeAll)
        .small_thresh(0)
        .large_graph_thresh(0)
        .filter_list(None)
        .cmdline("test")
        .version(TEST_VERSION)
        .log(&log)
        .build();
    quantify(quant_opts).unwrap();

    let rows =
        std::fs::read_to_string(quant_dir.join("alevin").join("quants_mat_rows.txt")).unwrap();
    let row_labels: Vec<&str> = rows.lines().collect();
    let unique_rows: HashSet<&str> = row_labels.iter().copied().collect();
    assert_eq!(row_labels.len(), num_samples * cells_per_sample);
    assert_eq!(unique_rows.len(), row_labels.len());
    assert!(row_labels.iter().any(|r| r.starts_with("sample_a_")));
    assert!(row_labels.iter().any(|r| r.starts_with("sample_b_")));
}

/// Regression test for COMBINE-lab/simpleaf#195.
///
/// The bug: collate.rs used to write the *sparse plate index* `sample_idx`
/// (a position in the chemistry's full sample-BC list, e.g. 0..384 for 10x
/// Flex v2) into each record's `barcodes[0]`.  But quant.rs builds a
/// *densely packed* `sample_names: Vec<String>` of length = number of
/// present samples and indexes that Vec directly with `barcodes[0]`.  When
/// the wells actually used in a run don't occupy contiguous positions
/// 0..N, the high sidxs miss the dense Vec, `sample_name` becomes None,
/// the row is written with 9 fields instead of 10, and downstream Polars
/// parsing in af-anndata infers `sample_name` as Int64 and crashes on
/// the first real string value.
///
/// This test wires up a sample BC list with **8 entries** but generates
/// reads for only positions **{0, 3, 7}** (non-contiguous), then runs the
/// full generate_permit_list → collate → quantify pipeline and asserts
/// three invariant families:
///
/// 1. Field-shape:   every row of featureDump.txt has exactly 10 fields.
/// 2. Identity:      all three expected sample names appear in
///                   quants_mat_rows.txt.
/// 3. Polars-shape:  the file parses with af-anndata's CSV schema
///                   (CB+sample_name as String) without error, and the
///                   sample_name column dtype is String with exactly the
///                   three expected values.
#[test]
fn test_multi_bc_quant_handles_sparse_sample_positions() {
    use alevin_fry::cellfilter::{CellFilterMethod, generate_permit_list};
    use alevin_fry::collate::collate;
    use alevin_fry::prog_opts::{GenPermitListOpts, QuantOpts, SampleCorrectionMode};
    use alevin_fry::quant::{ResolutionStrategy, SplicedAmbiguityModel, quantify};
    use bio_types::strand::Strand;
    use polars::prelude::{DataType, Field, Schema};
    use polars_io::csv::read::{CsvParseOptions, CsvReadOptions};
    use polars_io::SerReader;
    use std::collections::HashSet;
    use std::sync::Arc;

    let tmp = tempfile::tempdir().unwrap();
    let rad_dir = tmp.path().join("rad");
    std::fs::create_dir_all(&rad_dir).unwrap();
    let output_dir = tmp.path().join("output");
    std::fs::create_dir_all(&output_dir).unwrap();
    let quant_dir = tmp.path().join("quant");

    // 8-entry "plate" — only positions {0, 3, 7} will actually have reads,
    // mimicking a user-multiplexed run on a small subset of chemistry wells.
    let all_bcs: Vec<u64> = (0..8u64)
        .map(|i| make_packed_bc(100 + i * 37, SAMPLE_BC_LEN))
        .collect();
    let names: Vec<String> = (0..8).map(|i| format!("sample_{:02}", i)).collect();
    let sample_entries: Vec<(&str, u64)> = names
        .iter()
        .map(|s| s.as_str())
        .zip(all_bcs.iter().copied())
        .collect();

    // Positions used for the synthetic reads.  These MUST be non-contiguous
    // and include at least one value > num_present_samples - 1 (here, > 2)
    // to actually trigger the bug.
    let used_positions: [usize; 3] = [0, 3, 7];
    let used_bcs: Vec<u64> = used_positions.iter().map(|&p| all_bcs[p]).collect();
    let expected_sample_names: HashSet<String> = used_positions
        .iter()
        .map(|&p| names[p].clone())
        .collect();
    let num_used = used_bcs.len();
    let cells_per_sample = 4;
    let reads_per_cell = 8;

    // RAD records exist only for the used positions' BCs.
    let rad_path = rad_dir.join("map.rad");
    create_synthetic_multi_bc_rad_with_shared_cells(
        &rad_path,
        num_used,
        cells_per_sample,
        reads_per_cell,
        &used_bcs,
        true,
    )
    .unwrap();

    // Sample BC list includes ALL 8 entries — generate_permit_list will
    // iterate over the full plate, find reads only for positions 0/3/7,
    // and emit sample_info.json with sparse coverage.
    let sample_list_path = tmp.path().join("sample_barcodes.tsv");
    write_named_sample_bc_list(&sample_list_path, &sample_entries, SAMPLE_BC_LEN).unwrap();

    let tg_map_path = tmp.path().join("tg_map.tsv");
    write_tg_map(&tg_map_path).unwrap();

    let log = make_test_logger();
    let gpl_opts = GenPermitListOpts::builder()
        .input_dir(&rad_dir)
        .output_dir(&output_dir)
        .fmeth(CellFilterMethod::ForceCells(cells_per_sample))
        .expected_ori(Strand::Unknown)
        .version(TEST_VERSION)
        .threads(2)
        .velo_mode(false)
        .cmdline("test")
        .log(&log)
        .sample_bc_list(Some(sample_list_path))
        .sample_names(None)
        .sample_correction_mode(SampleCorrectionMode::Exact)
        .build();
    let total_cells = generate_permit_list(gpl_opts).unwrap();
    assert_eq!(total_cells, (num_used * cells_per_sample) as u64);

    collate(
        output_dir.clone(),
        &rad_dir,
        2,
        1_000,
        false,
        "test",
        TEST_VERSION,
        &log,
    )
    .unwrap();

    let quant_opts = QuantOpts::builder()
        .input_dir(&output_dir)
        .tg_map(&tg_map_path)
        .output_dir(&quant_dir)
        .num_threads(2)
        .num_bootstraps(0)
        .init_uniform(false)
        .summary_stat(false)
        .dump_eq(false)
        .resolution(ResolutionStrategy::Trivial)
        .pug_exact_umi(false)
        .sa_model(SplicedAmbiguityModel::WinnerTakeAll)
        .small_thresh(0)
        .large_graph_thresh(0)
        .filter_list(None)
        .cmdline("test")
        .version(TEST_VERSION)
        .log(&log)
        .build();
    quantify(quant_opts).unwrap();

    // ---- Assertion family 1: field-shape ----
    // Every featureDump.txt row (including header) must have 10 tab-separated
    // fields.  Pre-fix, ~2/3 of rows had 9 fields.
    let feat_dump_path = quant_dir.join("featureDump.txt");
    let feat_dump_text = std::fs::read_to_string(&feat_dump_path).unwrap();
    let mut field_counts: HashSet<usize> = HashSet::new();
    let mut total_data_rows = 0usize;
    for (i, line) in feat_dump_text.lines().enumerate() {
        let nf = line.split('\t').count();
        field_counts.insert(nf);
        if i > 0 {
            total_data_rows += 1;
        }
    }
    assert_eq!(
        field_counts,
        HashSet::from([10]),
        "featureDump.txt has rows with non-10 field counts (field counts seen: {:?})",
        field_counts
    );
    assert_eq!(
        total_data_rows,
        num_used * cells_per_sample,
        "expected {} featureDump rows, found {}",
        num_used * cells_per_sample,
        total_data_rows
    );

    // ---- Assertion family 2: identity (row labels match expected samples) ----
    let rows =
        std::fs::read_to_string(quant_dir.join("alevin").join("quants_mat_rows.txt")).unwrap();
    let row_labels: Vec<&str> = rows.lines().collect();
    assert_eq!(row_labels.len(), num_used * cells_per_sample);
    for n in &expected_sample_names {
        let prefix = format!("{}_", n);
        assert!(
            row_labels.iter().any(|r| r.starts_with(&prefix)),
            "no rows in quants_mat_rows.txt prefixed with {:?} (rows: {:?})",
            prefix,
            row_labels
        );
    }

    // ---- Assertion family 3: polars-shape ----
    // Parse featureDump.txt with af-anndata's schema (CB+sample_name as
    // String) — pre-fix this used to crash with a dtype-inference error.
    let feat_dump_schema = Arc::new(Schema::from_iter([
        Field::new("CB".into(), DataType::String),
        Field::new("sample_name".into(), DataType::String),
        Field::new("CorrectedReads".into(), DataType::Int64),
        Field::new("MappedReads".into(), DataType::Int64),
        Field::new("DeduplicatedReads".into(), DataType::Float64),
        Field::new("MappingRate".into(), DataType::Float64),
        Field::new("DedupRate".into(), DataType::Float64),
        Field::new("MeanByMax".into(), DataType::Float64),
        Field::new("NumGenesExpressed".into(), DataType::Int64),
        Field::new("NumGenesOverMean".into(), DataType::Int64),
    ]));
    let parse_options = CsvParseOptions::default().with_separator(b'\t');
    let df = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .with_schema_overwrite(Some(feat_dump_schema))
        .with_raise_if_empty(true)
        .try_into_reader_with_file_path(Some(feat_dump_path.clone()))
        .unwrap()
        .finish()
        .expect("polars CSV read of featureDump.txt must succeed");

    assert_eq!(df.width(), 10, "expected 10 columns, got {}", df.width());

    let sample_name_col = df.column("sample_name").expect("sample_name column");
    assert_eq!(
        sample_name_col.dtype(),
        &DataType::String,
        "sample_name dtype should be String, got {:?}",
        sample_name_col.dtype()
    );

    let observed_names: HashSet<String> = sample_name_col
        .str()
        .expect("sample_name should be a string series")
        .into_iter()
        .filter_map(|opt| opt.map(|s| s.to_string()))
        .collect();
    assert_eq!(
        observed_names, expected_sample_names,
        "sample_name values mismatch — expected {:?}, got {:?}",
        expected_sample_names, observed_names
    );
}

/// End-to-end Flex v2 smoke test for issue #195 against real subset data.
///
/// This test is `#[ignore]` by default — it requires pre-staged inputs
/// (a mapped multi-barcode RAD directory, the cached Flex v2 sample BC
/// list, and a t2g map).  Run it with:
///
/// ```text
/// AF_TEST_FLEXV2_RAD=/path/to/rad_dir \
/// AF_TEST_FLEXV2_SAMPLE_BC_LIST=/path/to/737K-flex-v2.txt \
/// AF_TEST_FLEXV2_TG_MAP=/path/to/t2g.tsv \
/// cargo test -p alevin-fry --test multi_barcode_integration -- \
///     --ignored flexv2_real_data --nocapture
/// ```
///
/// Staging the RAD directory (one-time, outside this test):
/// 1. Subset the Flex v2 FASTQs (e.g. with `zcat | head -n 4000000`) so
///    only reads from a small number of wells appear.  The fastest way is
///    to use a sample of `/scratch2/rob/read_data/10x_flexv2/`.
/// 2. Run `simpleaf workflow` (or piscem directly) to produce a
///    multi-barcode `map.rad` against a Flex v2 index.
/// 3. Point AF_TEST_FLEXV2_RAD at the directory containing that map.rad.
///
/// The test does NOT hard-code which wells are present — it reads
/// sample_info.json after generate_permit_list to learn which wells got
/// non-zero reads, then asserts those wells all appear as sample_name
/// values in featureDump.txt.  Bug pre-fix: only the *first* such well
/// appears; this test fails with a "missing sample names" assertion.
#[test]
#[ignore]
fn test_multi_bc_quant_flexv2_real_data() {
    use alevin_fry::cellfilter::{CellFilterMethod, generate_permit_list};
    use alevin_fry::collate::collate;
    use alevin_fry::prog_opts::{GenPermitListOpts, QuantOpts, SampleBarcodeOri, SampleCorrectionMode};
    use alevin_fry::quant::{ResolutionStrategy, SplicedAmbiguityModel, quantify};
    use bio_types::strand::Strand;
    use polars::prelude::{DataType, Field, Schema};
    use polars_io::csv::read::{CsvParseOptions, CsvReadOptions};
    use polars_io::SerReader;
    use std::collections::HashSet;
    use std::path::PathBuf;
    use std::sync::Arc;

    let rad_dir = match std::env::var("AF_TEST_FLEXV2_RAD") {
        Ok(v) => PathBuf::from(v),
        Err(_) => {
            eprintln!(
                "Skipping flexv2_real_data: AF_TEST_FLEXV2_RAD unset \
                 (point it at a pre-staged Flex v2 map.rad directory)"
            );
            return;
        }
    };
    let sample_bc_list = match std::env::var("AF_TEST_FLEXV2_SAMPLE_BC_LIST") {
        Ok(v) => PathBuf::from(v),
        Err(_) => {
            eprintln!("Skipping flexv2_real_data: AF_TEST_FLEXV2_SAMPLE_BC_LIST unset");
            return;
        }
    };
    let tg_map_path = match std::env::var("AF_TEST_FLEXV2_TG_MAP") {
        Ok(v) => PathBuf::from(v),
        Err(_) => {
            eprintln!("Skipping flexv2_real_data: AF_TEST_FLEXV2_TG_MAP unset");
            return;
        }
    };
    if !rad_dir.join("map.rad").exists() {
        eprintln!(
            "Skipping flexv2_real_data: {} does not contain map.rad",
            rad_dir.display()
        );
        return;
    }
    // 10x Flex v2 sample barcodes appear reverse-complemented relative to the
    // published whitelist, so the orientation must be "reverse" for that
    // chemistry.  Controlled by an env var (default "forward") so the test
    // also works for forward-orientation chemistries.
    let sample_bc_ori = match std::env::var("AF_TEST_FLEXV2_SAMPLE_BC_ORI")
        .unwrap_or_else(|_| "forward".to_string())
        .to_ascii_lowercase()
        .as_str()
    {
        "reverse" | "rc" => SampleBarcodeOri::Reverse,
        _ => SampleBarcodeOri::Forward,
    };

    let tmp = tempfile::tempdir().unwrap();
    let output_dir = tmp.path().join("output");
    std::fs::create_dir_all(&output_dir).unwrap();
    let quant_dir = tmp.path().join("quant");

    let log = make_test_logger();
    let gpl_opts = GenPermitListOpts::builder()
        .input_dir(&rad_dir)
        .output_dir(&output_dir)
        // ForceCells with a generous cap — the data drives the actual count.
        .fmeth(CellFilterMethod::ForceCells(10_000))
        .expected_ori(Strand::Unknown)
        .version(TEST_VERSION)
        .threads(4)
        .velo_mode(false)
        .cmdline("test_flexv2_real_data")
        .log(&log)
        .sample_bc_list(Some(sample_bc_list))
        .sample_names(None)
        .sample_correction_mode(SampleCorrectionMode::Exact)
        .sample_bc_ori(sample_bc_ori)
        .build();
    generate_permit_list(gpl_opts).expect("generate_permit_list failed");

    // Discover the wells that actually have reads from sample_info.json.
    // For #195 to manifest, the cached BC list must contain entries for
    // wells the data does NOT use — typical for the 384-entry Flex v2 file.
    let sample_info: serde_json::Value =
        serde_json::from_reader(std::fs::File::open(output_dir.join("sample_info.json")).unwrap())
            .unwrap();
    let expected_sample_names: HashSet<String> = sample_info["samples"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|e| e["num_reads"].as_u64().unwrap_or(0) > 0)
        .map(|e| e["name"].as_str().unwrap().to_string())
        .collect();
    eprintln!(
        "Flexv2 real-data test: {} wells have reads ({:?})",
        expected_sample_names.len(),
        expected_sample_names
    );
    assert!(
        !expected_sample_names.is_empty(),
        "no wells with reads — the input RAD likely has no sample matches",
    );

    collate(
        output_dir.clone(),
        &rad_dir,
        4,
        100_000,
        false,
        "test_flexv2_real_data",
        TEST_VERSION,
        &log,
    )
    .expect("collate failed");

    let quant_opts = QuantOpts::builder()
        .input_dir(&output_dir)
        .tg_map(&tg_map_path)
        .output_dir(&quant_dir)
        .num_threads(4)
        .num_bootstraps(0)
        .init_uniform(false)
        .summary_stat(false)
        .dump_eq(false)
        .resolution(ResolutionStrategy::CellRangerLike)
        .pug_exact_umi(false)
        .sa_model(SplicedAmbiguityModel::WinnerTakeAll)
        .small_thresh(0)
        .large_graph_thresh(0)
        .filter_list(None)
        .cmdline("test_flexv2_real_data")
        .version(TEST_VERSION)
        .log(&log)
        .build();
    quantify(quant_opts).expect("quantify failed");

    // ---- Assertion family 1: field-shape ----
    let feat_dump_path = quant_dir.join("featureDump.txt");
    let feat_dump_text = std::fs::read_to_string(&feat_dump_path).unwrap();
    let mut field_counts: HashSet<usize> = HashSet::new();
    for line in feat_dump_text.lines() {
        field_counts.insert(line.split('\t').count());
    }
    assert_eq!(
        field_counts,
        HashSet::from([10]),
        "featureDump.txt has mixed-width rows (field counts: {:?})",
        field_counts
    );

    // ---- Assertion family 2: identity ----
    // Every well that had reads in sample_info.json should appear as a
    // row prefix in quants_mat_rows.txt (multi-sample naming is
    // "{sample_name}_{cell_barcode}").
    let rows =
        std::fs::read_to_string(quant_dir.join("alevin").join("quants_mat_rows.txt")).unwrap();
    let row_labels: Vec<&str> = rows.lines().collect();
    for n in &expected_sample_names {
        let prefix = format!("{}_", n);
        assert!(
            row_labels.iter().any(|r| r.starts_with(&prefix)),
            "no rows in quants_mat_rows.txt prefixed with {:?} \
             (sample present in sample_info.json with non-zero reads)",
            prefix
        );
    }

    // ---- Assertion family 3: polars-shape ----
    let feat_dump_schema = Arc::new(Schema::from_iter([
        Field::new("CB".into(), DataType::String),
        Field::new("sample_name".into(), DataType::String),
        Field::new("CorrectedReads".into(), DataType::Int64),
        Field::new("MappedReads".into(), DataType::Int64),
        Field::new("DeduplicatedReads".into(), DataType::Float64),
        Field::new("MappingRate".into(), DataType::Float64),
        Field::new("DedupRate".into(), DataType::Float64),
        Field::new("MeanByMax".into(), DataType::Float64),
        Field::new("NumGenesExpressed".into(), DataType::Int64),
        Field::new("NumGenesOverMean".into(), DataType::Int64),
    ]));
    let parse_options = CsvParseOptions::default().with_separator(b'\t');
    let df = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .with_schema_overwrite(Some(feat_dump_schema))
        .with_raise_if_empty(true)
        .try_into_reader_with_file_path(Some(feat_dump_path))
        .unwrap()
        .finish()
        .expect("polars CSV read of featureDump.txt must succeed");

    assert_eq!(df.width(), 10);
    let sample_name_col = df.column("sample_name").unwrap();
    assert_eq!(sample_name_col.dtype(), &DataType::String);

    let observed: HashSet<String> = sample_name_col
        .str()
        .unwrap()
        .into_iter()
        .filter_map(|opt| opt.map(|s| s.to_string()))
        .collect();
    assert_eq!(
        observed, expected_sample_names,
        "featureDump.txt sample_name values do not match sample_info.json non-zero wells"
    );
}

/// Test the collation manifest roundtrip.
#[test]
fn test_collation_manifest_roundtrip() {
    use libradicl::collation::{CollationManifest, SampleGroup};

    let tmp = tempfile::tempdir().unwrap();
    let manifest_path = tmp.path().join("collation_manifest.bin");

    let mut manifest = CollationManifest::new(vec!["sample".to_string(), "cell".to_string()]);
    manifest.add_sample_group(SampleGroup {
        key: 0x1234,
        name: Some("sample_A".to_string()),
        chunk_start: 0,
        num_chunks: 100,
        num_records: 50000,
    });
    manifest.add_sample_group(SampleGroup {
        key: 0x5678,
        name: Some("sample_B".to_string()),
        chunk_start: 100,
        num_chunks: 80,
        num_records: 40000,
    });

    manifest.write_to_file(&manifest_path).unwrap();

    let restored = CollationManifest::read_from_file(&manifest_path).unwrap();
    assert_eq!(restored.levels.len(), 2);
    assert_eq!(restored.sample_groups.len(), 2);
    assert_eq!(restored.sample_groups[0].name, Some("sample_A".to_string()));
    assert_eq!(restored.sample_groups[1].num_chunks, 80);
    assert_eq!(restored.total_chunks(), 180);
    assert_eq!(restored.total_records(), 90000);
}

/// Test reading from the real Flex RAD file (if available)
#[test]
fn test_read_real_flex_rad() {
    let rad_path = std::path::Path::new("flex_data/map_output3/map.rad");
    if !rad_path.exists() {
        eprintln!("Skipping test_read_real_flex_rad: RAD file not found");
        return;
    }

    let f = std::fs::File::open(rad_path).unwrap();
    let mut reader = std::io::BufReader::new(f);

    let prelude = libradicl::header::RadPrelude::from_bytes(&mut reader).unwrap();
    let _ftm = prelude
        .file_tags
        .parse_tags_from_bytes(&mut reader)
        .unwrap();

    let ctx = prelude
        .get_record_context::<libradicl::record::MultiBarcodeRecordContext>()
        .unwrap();
    eprintln!("Context: {:?}", ctx);

    // Read first chunk
    let chunk = libradicl::chunk::Chunk::<MultiBarcodeReadRecord>::from_bytes(&mut reader, &ctx);
    eprintln!(
        "Chunk: nbytes={}, nrec={}, reads.len()={}",
        chunk.nbytes,
        chunk.nrec,
        chunk.reads.len()
    );

    assert!(chunk.nrec > 0, "First chunk should have records");
    assert_eq!(
        chunk.reads.len(),
        chunk.nrec as usize,
        "reads.len() should match nrec"
    );

    let r = &chunk.reads[0];
    eprintln!(
        "First read: barcodes={:?}, umi={}, refs={:?}",
        r.barcodes.as_slice(),
        r.umi,
        &r.refs
    );
    assert_eq!(
        r.barcodes.len(),
        2,
        "Should have 2 barcodes (sample + cell)"
    );
}
