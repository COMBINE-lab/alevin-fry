//! Integration tests for the self-describing "tag roles" (#64/#66) read paths.
//!
//! These build a **role-only** single-barcode RAD — the barcode and UMI read
//! tags carry `TagRole::Barcode`/`TagRole::Umi` but use non-conventional names
//! (not `b`/`u`), and there is no `cblen`/`num_barcodes` file tag — and drive it
//! through generate-permit-list and collate. Before PR #195 such a file was
//! rejected by generate-permit-list (which keyed on the `b`/`u` tag NAMES), which
//! also made collate's role auto-route (`!has_bridge && has_barcode_role`)
//! unreachable via the real pipeline. This asserts both now work end to end.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use libradicl::chunk::Chunk;
use libradicl::codec::ChunkCodec;
use libradicl::header::{RadHeader, RadPrelude, SpecVersion};
use libradicl::rad_types::{
    RadIntId, RadType, TagDesc, TagMap, TagRole, TagSection, TagSectionLabel,
};
use libradicl::record::{AlevinFryReadRecord, AlevinFryRecordContext, RecordContext};
use libradicl::writers::RadFileWriter;

const NUM_REFS: u64 = 4;
const CELL_BC_LEN: u16 = 16;
const UMI_LEN: u16 = 12;
const TEST_VERSION: &str = "0.12.0";

fn packed(idx: u64, len: u16) -> u64 {
    let mask = (1u64 << (2 * len as u64)) - 1;
    idx.wrapping_mul(2654435761) & mask
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

/// Build a role-only single-barcode prelude: read tags `cell_bc` (Barcode role,
/// len 16) and `umi_tag` (Umi role, len 12) — deliberately NOT named `b`/`u` — and
/// an alignment tag with no declared role. Spec major 2 so the roles round-trip.
/// No `cblen`/`num_barcodes` file tags: the lengths live only on the roles.
fn make_role_only_prelude() -> (RadPrelude, TagMap) {
    let ref_names: Vec<String> = (0..NUM_REFS).map(|i| format!("gene_{i}")).collect();
    let hdr = RadHeader {
        version: SpecVersion::current(),
        is_paired: 0,
        ref_count: NUM_REFS,
        ref_names,
        num_chunks: 0, // backpatched by RadFileWriter
    };

    let file_tags = TagSection::new_with_label(TagSectionLabel::FileTags);

    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    read_tags.add_tag_desc(
        TagDesc::new("cell_bc", RadType::Int(RadIntId::U32)).with_role(TagRole::Barcode {
            level: 0,
            len: CELL_BC_LEN as u8,
        }),
    );
    read_tags.add_tag_desc(
        TagDesc::new("umi_tag", RadType::Int(RadIntId::U32))
            .with_role(TagRole::Umi { len: UMI_LEN as u8 }),
    );

    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    aln_tags.add_tag_desc(TagDesc::new(
        "compressed_ori_refid",
        RadType::Int(RadIntId::U32),
    ));

    let prelude = RadPrelude::from_header_and_tag_sections(hdr, file_tags, read_tags, aln_tags);
    let file_tag_map = TagMap::with_keyset(&prelude.file_tags.tags); // empty (no file tags)
    (prelude, file_tag_map)
}

/// Write a role-only single-barcode RAD: `num_cells` cells, each one chunk of
/// `reads_per_cell` reads mapping to a single reference.
fn create_role_only_rad(
    path: &Path,
    num_cells: usize,
    reads_per_cell: usize,
) -> anyhow::Result<()> {
    let (prelude, file_tag_map) = make_role_only_prelude();
    // The parsing context is built from the declared roles (non-`b`/`u` names).
    let ctx = AlevinFryRecordContext::get_context_prefer_roles(
        &prelude.file_tags,
        &prelude.read_tags,
        &prelude.aln_tags,
    )?;

    let file = File::create(path)?;
    let mut fw = RadFileWriter::new(file, &prelude, &file_tag_map)?;

    for cell_idx in 0..num_cells {
        let cell_bc = packed(cell_idx as u64, CELL_BC_LEN);
        let mut reads = Vec::with_capacity(reads_per_cell);
        for read_idx in 0..reads_per_cell {
            let umi = packed((cell_idx * 1000 + read_idx) as u64, UMI_LEN);
            let ref_id = (read_idx as u64 % NUM_REFS) as u32;
            reads.push(AlevinFryReadRecord {
                bc: cell_bc,
                umi,
                dirs: vec![true],
                refs: vec![ref_id],
            });
        }
        let chunk = Chunk::<AlevinFryReadRecord> {
            nbytes: 0,
            nrec: reads.len() as u32,
            reads,
        };
        fw.write_chunk(&chunk, &ctx)?;
    }
    fw.finalize()?;
    Ok(())
}

/// A role-only single-barcode RAD flows through generate-permit-list (which now
/// accepts the Barcode/Umi ROLES rather than only `b`/`u` names) and then through
/// collate, which auto-routes to the tag-driven generic gather because the file
/// declares a barcode role but lacks the `b`/`u` name bridge.
#[test]
fn role_only_single_barcode_flows_through_gpl_and_collate() {
    use alevin_fry::cellfilter::{CellFilterMethod, generate_permit_list};
    use alevin_fry::collate::collate;
    use alevin_fry::prog_opts::GenPermitListOpts;
    use bio_types::strand::Strand;

    let tmp = tempfile::tempdir().unwrap();
    let rad_dir = tmp.path().join("rad");
    std::fs::create_dir_all(&rad_dir).unwrap();
    let output_dir = tmp.path().join("output");
    std::fs::create_dir_all(&output_dir).unwrap();

    let num_cells = 6usize;
    let reads_per_cell = 8usize;
    create_role_only_rad(&rad_dir.join("map.rad"), num_cells, reads_per_cell).unwrap();

    let log = make_test_logger();

    // generate-permit-list: expected_ori Unknown so no Orientation role is needed.
    let gpl_opts = GenPermitListOpts::builder()
        .input_dir(&rad_dir)
        .output_dir(&output_dir)
        .fmeth(CellFilterMethod::ForceCells(num_cells))
        .expected_ori(Strand::Unknown)
        .version(TEST_VERSION)
        .threads(2)
        .velo_mode(false)
        .cmdline("test")
        .log(&log)
        .build();
    // The single-barcode filtered path returns the count of *corrected* (error)
    // barcodes, which is legitimately 0 for these exact barcodes; the retained
    // cells are proved below by the collated chunk count. The key assertion for
    // fix 2 is that this call SUCCEEDS at all — before it, a role-only RAD was
    // rejected because `validate_tag_types` keyed on the `b`/`u` names.
    generate_permit_list(gpl_opts).unwrap();
    assert!(output_dir.join("generate_permit_list.json").exists());
    assert!(output_dir.join("permit_freq.bin").exists());

    // collate: auto-routes to the tag-driven generic gather (role key), since the
    // read tags declare a Barcode role but are not named `b`/`u`.
    collate(
        output_dir.clone(),
        &rad_dir,
        2,
        1_000,
        ChunkCodec::None,
        "test",
        TEST_VERSION,
        &log,
    )
    .unwrap();

    // The collated output exists and holds exactly one chunk per retained cell.
    let collated_path = output_dir.join("map.collated.rad");
    assert!(
        collated_path.exists(),
        "collate must produce map.collated.rad for the role-only RAD"
    );
    let f = File::open(&collated_path).unwrap();
    let mut reader = BufReader::new(f);
    let collated_prelude = RadPrelude::from_bytes(&mut reader).unwrap();
    assert_eq!(
        collated_prelude.hdr.num_chunks, num_cells as u64,
        "collated RAD should have one chunk per retained cell"
    );
    // The collated prelude carries the same declared Barcode role.
    assert!(
        collated_prelude
            .read_tags
            .tags
            .iter()
            .any(|t| matches!(t.role, TagRole::Barcode { .. })),
        "collated prelude should preserve the declared Barcode role"
    );
}
