//! Integration tests for the multi-barcode (10x Flex) pipeline.
//!
//! These tests create synthetic multi-barcode RAD files programmatically,
//! then exercise the generate-permit-list → collate → quant pipeline.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Cursor, Read, Write};
use std::path::{Path, PathBuf};

use libradicl::chunk::Chunk;
use libradicl::collation::CollationManifest;
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{RadIntId, RadType, TagDesc, TagMap, TagSection, TagSectionLabel, TagValue};
use libradicl::record::{
    MultiBarcodeReadRecord, MultiBarcodeReadRecordT, MultiBarcodeRecordContext, RecordContext,
};
use libradicl::writers::RadFileWriter;
use smallvec::smallvec;

// -----------------------------------------------------------------------
// Synthetic RAD file builder
// -----------------------------------------------------------------------

/// Number of reference targets in test data
const NUM_REFS: u64 = 10;
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
    let mixed = idx.wrapping_mul(2654435761) & mask;
    mixed
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
    file_tag_map.add(TagValue::U16(2));                    // num_barcodes
    file_tag_map.add(TagValue::U16(SAMPLE_BC_LEN));        // b0len
    file_tag_map.add(TagValue::U16(CELL_BC_LEN));          // b1len
    file_tag_map.add(TagValue::U16(UMI_LEN));              // ulen
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
        for sample_idx in 0..num_samples {
            let sample_bc = sample_barcodes[sample_idx];
            let cell_bc = make_packed_bc((sample_idx * 1000 + cell_idx) as u64, CELL_BC_LEN);

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

/// Write a transcript-to-gene map file.
fn write_tg_map(path: &Path) -> anyhow::Result<()> {
    let mut f = BufWriter::new(File::create(path)?);
    for i in 0..NUM_REFS {
        // Each transcript maps to a gene with the same name
        writeln!(f, "gene_{}\tgene_{}", i, i)?;
    }
    Ok(())
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
        &rad_path, num_samples, cells_per_sample, reads_per_cell, &sample_bcs,
    )
    .unwrap();

    // Read it back
    let file = File::open(&rad_path).unwrap();
    let mut reader = BufReader::new(file);
    let prelude = RadPrelude::from_bytes(&mut reader).unwrap();
    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut reader).unwrap();

    // Verify header
    let expected_chunks = (num_samples * cells_per_sample) as u64;
    assert_eq!(prelude.hdr.num_chunks, expected_chunks);
    assert_eq!(prelude.hdr.ref_count, NUM_REFS);

    // Verify file tags
    let num_bc: u16 = file_tag_map.get("num_barcodes").unwrap().try_into().unwrap();
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
    use alevin_fry::cellfilter::{generate_permit_list, CellFilterMethod};
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
        &rad_path, num_samples, cells_per_sample, reads_per_cell, &sample_bcs,
    )
    .unwrap();

    // Write sample barcode list
    let sample_list_path = tmp.path().join("sample_barcodes.txt");
    write_sample_bc_list(&sample_list_path, &sample_bcs, SAMPLE_BC_LEN).unwrap();

    // Setup logger
    let decorator = slog_term::PlainDecorator::new(slog_term::TestStdoutWriter);
    let drain = slog_term::CompactFormat::new(decorator).build();
    let drain = std::sync::Mutex::new(drain);
    let drain = slog::Fuse::new(drain);
    let log = slog::Logger::root(drain, slog::o!());

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
    assert!(total_cells > 0, "expected at least some cells to be retained");

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
            assert!(sample_dir.join("permit_map.bin").exists(),
                "permit_map.bin should exist for sample {} with cells", name);
            assert!(sample_dir.join("permit_freq.bin").exists(),
                "permit_freq.bin should exist for sample {} with cells", name);
        }
    }
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
