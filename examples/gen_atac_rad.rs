// Generate a scaled synthetic scATAC RAD (+ whitelist.txt + empty
// unmapped_bc_count.bin) for benchmarking collate/deduplicate. Mirrors the
// fixture in tests/atac_integration.rs but parameterized and larger.
//
// Usage: gen_atac_rad <rad_dir> <num_cells> <good> <unmapped> <multi>
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{
    RadAtomicId, RadIntId, RadType, TagDesc, TagMap, TagSection, TagSectionLabel, TagValue,
};
use libradicl::record::{AtacSeqReadRecord, AtacSeqRecordContext, RecordContext};
use libradicl::{RadFileWriter, chunk::Chunk};
use std::fs::File;
use std::io::Write;

const CELL_BC_LEN: u16 = 16;
const REF_NAMES: [&str; 4] = ["chr1", "chr2", "chr3", "chr4"];
const REF_LENGTHS: [u32; 4] = [250_000_000, 240_000_000, 200_000_000, 190_000_000];

fn make_packed_bc(idx: u64, len: u16) -> u64 {
    let mask = (1u64 << (2 * len as u64)) - 1;
    idx.wrapping_mul(2654435761) & mask
}
fn packed_to_nuc(packed: u64, len: usize) -> String {
    const NUCS: [char; 4] = ['A', 'C', 'G', 'T'];
    (0..len)
        .map(|i| NUCS[((packed >> (2 * (len - 1 - i))) & 3) as usize])
        .collect()
}

fn make_prelude() -> (RadPrelude, TagMap) {
    let hdr = RadHeader {
        version: libradicl::header::SpecVersion::Legacy,
        is_paired: 1,
        ref_count: REF_NAMES.len() as u64,
        ref_names: REF_NAMES.iter().map(|s| s.to_string()).collect(),
        num_chunks: 0,
    };
    let mut file_tags = TagSection::new_with_label(TagSectionLabel::FileTags);
    file_tags.add_tag_desc(TagDesc::new("cblen", RadType::Int(RadIntId::U16)));
    file_tags.add_tag_desc(TagDesc::new("known_rad_type", RadType::String));
    file_tags.add_tag_desc(TagDesc::new(
        "ref_lengths",
        RadType::Array(RadIntId::U32, RadAtomicId::Int(RadIntId::U32)),
    ));
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    read_tags.add_tag_desc(TagDesc::new("b", RadType::Int(RadIntId::U32)));
    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    for (name, typeid) in [
        ("ref", RadType::Int(RadIntId::U32)),
        ("type", RadType::Int(RadIntId::U8)),
        ("start_pos", RadType::Int(RadIntId::U32)),
        ("frag_len", RadType::Int(RadIntId::U16)),
    ] {
        aln_tags.add_tag_desc(TagDesc::new(name, typeid));
    }
    let prelude = RadPrelude {
        hdr,
        file_tags,
        read_tags,
        aln_tags,
    };
    let mut m = TagMap::with_keyset(&prelude.file_tags.tags);
    m.add(TagValue::U16(CELL_BC_LEN));
    m.add(TagValue::String("sc_atac".to_string()));
    m.add(TagValue::ArrayU32(REF_LENGTHS.to_vec()));
    (prelude, m)
}

fn main() -> anyhow::Result<()> {
    let mut a = std::env::args().skip(1);
    let dir = a.next().expect("rad_dir");
    let num_cells: usize = a.next().expect("num_cells").parse()?;
    let good: usize = a.next().expect("good").parse()?;
    let unmapped: usize = a.next().expect("unmapped").parse()?;
    let multi: usize = a.next().expect("multi").parse()?;
    std::fs::create_dir_all(&dir)?;

    let (prelude, ftm) = make_prelude();
    let ctx = AtacSeqRecordContext::get_context_from_tag_section(
        &prelude.file_tags,
        &prelude.read_tags,
        &prelude.aln_tags,
    )?;
    let f = File::create(format!("{dir}/map.rad"))?;
    let mut fw = RadFileWriter::new(f, &prelude, &ftm)?;
    let mut wl = std::io::BufWriter::new(File::create(format!("{dir}/whitelist.txt"))?);

    for cell in 0..num_cells {
        let bc = make_packed_bc(cell as u64 + 1, CELL_BC_LEN);
        writeln!(wl, "{}", packed_to_nuc(bc, CELL_BC_LEN as usize))?;
        let mut reads = Vec::with_capacity(good + unmapped + multi);
        for r in 0..good {
            let ref_id = (r % REF_NAMES.len()) as u32;
            // Keep positions well within the reference (gpl bins by position).
            let span = REF_LENGTHS[ref_id as usize] as u64 - 200_000;
            let start = 1_000
                + (((cell as u64).wrapping_mul(977) + (r as u64).wrapping_mul(131)) % span) as u32;
            reads.push(AtacSeqReadRecord {
                bc,
                start_pos: vec![start],
                refs: vec![ref_id],
                frag_lengths: vec![120],
                map_type: vec![4],
            });
        }
        for _ in 0..unmapped {
            reads.push(AtacSeqReadRecord {
                bc,
                start_pos: vec![],
                refs: vec![],
                frag_lengths: vec![],
                map_type: vec![],
            });
        }
        for r in 0..multi {
            let span = REF_LENGTHS[0] as u64 - 200_000;
            let start = 50_000 + (((cell as u64).wrapping_mul(13) + (r as u64) * 7) % span) as u32;
            reads.push(AtacSeqReadRecord {
                bc,
                start_pos: vec![start, start + 250],
                refs: vec![0, 1],
                frag_lengths: vec![110, 115],
                map_type: vec![4, 4],
            });
        }
        let chunk = Chunk::<AtacSeqReadRecord> {
            nbytes: 0,
            nrec: reads.len() as u32,
            reads,
        };
        fw.write_chunk(&chunk, &ctx)?;
    }
    fw.finalize()?;
    wl.flush()?;
    // piscem always emits this; collate opens it unconditionally.
    File::create(format!("{dir}/unmapped_bc_count.bin"))?;
    eprintln!(
        "wrote {num_cells} cells (good={good} unmapped={unmapped} multi={multi}) -> {dir}/map.rad"
    );
    Ok(())
}
