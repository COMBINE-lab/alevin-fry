// Generate a synthetic multi-barcode (Flex-like) scRNA RAD, either in the legacy
// name-convention layout (b0/b1/u + num_barcodes/b0len/b1len/ulen file tags) or a
// role-only layout (sample_bc/cell_bc/umi carrying Barcode{level,len}/Umi roles,
// spec major 2, and NO num_barcodes/bNlen tags). Same records either way, so the
// two can be run through gpl->collate->quant and compared. Also writes the
// matching sample barcode list. For validating role-driven multi gpl (#66/B6).
//
// Usage: gen_multi_rad <out_dir> <names|roles> <num_samples> <cells_per_sample> <reads_per_cell>
//   also writes <out_dir>/../sample_bc.txt (the sample barcodes, one per line)
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{
    RadIntId, RadType, TagDesc, TagMap, TagRole, TagSection, TagSectionLabel, TagValue,
};
use libradicl::record::{MultiBarcodeReadRecord, MultiBarcodeRecordContext, RecordContext};
use libradicl::{chunk::Chunk, RadFileWriter};
use smallvec::smallvec;
use std::fs::File;
use std::io::Write;

const SAMPLE_LEN: u16 = 8;
const CELL_LEN: u16 = 16;
const UMI_LEN: u16 = 12;
const REF_NAMES: [&str; 3] = ["gene0", "gene1", "gene2"];

fn packed(idx: u64, len: u16) -> u64 {
    let mask = (1u64 << (2 * len as u64)) - 1;
    idx.wrapping_mul(2654435761) & mask
}
fn to_nuc(p: u64, len: usize) -> String {
    const N: [char; 4] = ['A', 'C', 'G', 'T'];
    (0..len)
        .map(|i| N[((p >> (2 * (len - 1 - i))) & 3) as usize])
        .collect()
}

fn desc(name: &str, typeid: RadType, role: TagRole) -> TagDesc {
    TagDesc::new(name, typeid).with_role(role)
}

fn main() -> anyhow::Result<()> {
    let mut a = std::env::args().skip(1);
    let dir = a.next().expect("out_dir");
    let mode = a.next().expect("names|roles");
    let ns: usize = a.next().expect("num_samples").parse()?;
    let cps: usize = a.next().expect("cells_per_sample").parse()?;
    let rpc: usize = a.next().expect("reads_per_cell").parse()?;
    let roles = mode == "roles";
    std::fs::create_dir_all(&dir)?;

    let u32i = || RadType::Int(RadIntId::U32);
    let u16i = || RadType::Int(RadIntId::U16);

    let mut file_tags = TagSection::new_with_label(TagSectionLabel::FileTags);
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    if roles {
        // Role-only: no num_barcodes/bNlen conventions; roles carry everything.
        file_tags.add_tag_desc(desc("known_rad_type", RadType::String, TagRole::None));
        read_tags.add_tag_desc(desc(
            "sample_bc",
            u32i(),
            TagRole::Barcode {
                level: 0,
                len: SAMPLE_LEN as u8,
            },
        ));
        read_tags.add_tag_desc(desc(
            "cell_bc",
            u32i(),
            TagRole::Barcode {
                level: 1,
                len: CELL_LEN as u8,
            },
        ));
        read_tags.add_tag_desc(desc("umi", u32i(), TagRole::Umi { len: 12 }));
        aln_tags.add_tag_desc(desc("compressed_ori_refid", u32i(), TagRole::Orientation));
    } else {
        for (n, t) in [
            ("num_barcodes", u16i()),
            ("b0len", u16i()),
            ("b1len", u16i()),
            ("ulen", u16i()),
            ("known_rad_type", RadType::String),
        ] {
            file_tags.add_tag_desc(desc(n, t, TagRole::None));
        }
        read_tags.add_tag_desc(desc("b0", u32i(), TagRole::None));
        read_tags.add_tag_desc(desc("b1", u32i(), TagRole::None));
        read_tags.add_tag_desc(desc("u", u32i(), TagRole::None));
        aln_tags.add_tag_desc(desc("compressed_ori_refid", u32i(), TagRole::None));
    }

    let hdr = RadHeader {
        version: if roles {
            libradicl::header::SpecVersion::current()
        } else {
            libradicl::header::SpecVersion::Legacy
        },
        is_paired: 0,
        ref_count: REF_NAMES.len() as u64,
        ref_names: REF_NAMES.iter().map(|s| s.to_string()).collect(),
        num_chunks: 0,
    };
    let prelude = RadPrelude {
        hdr,
        file_tags,
        read_tags,
        aln_tags,
    };
    let mut ftm = TagMap::with_keyset(&prelude.file_tags.tags);
    if roles {
        ftm.add(TagValue::String("sc_rna_multi_bc".to_string()));
    } else {
        ftm.add(TagValue::U16(2));
        ftm.add(TagValue::U16(SAMPLE_LEN));
        ftm.add(TagValue::U16(CELL_LEN));
        ftm.add(TagValue::U16(UMI_LEN));
        ftm.add(TagValue::String("sc_rna_multi_bc".to_string()));
    }

    let ctx = MultiBarcodeRecordContext::get_context_from_tag_section(
        &prelude.file_tags,
        &prelude.read_tags,
        &prelude.aln_tags,
    )
    .unwrap_or_else(|_| {
        // role-only prelude: build the context from roles instead of names
        MultiBarcodeRecordContext::from_roles(&prelude.read_tags)
            .unwrap()
            .expect("role context")
    });

    let f = File::create(format!("{dir}/map.rad"))?;
    let mut fw = RadFileWriter::new(f, &prelude, &ftm)?;

    // sample barcodes: a fixed known set
    let sample_bcs: Vec<u64> = (0..ns).map(|s| packed(s as u64 + 1, SAMPLE_LEN)).collect();
    for (s, &sbc) in sample_bcs.iter().enumerate() {
        for c in 0..cps {
            let cbc = packed((s * cps + c) as u64 + 1, CELL_LEN);
            let mut reads = Vec::with_capacity(rpc);
            for r in 0..rpc {
                let refid = (r % REF_NAMES.len()) as u32;
                reads.push(MultiBarcodeReadRecord {
                    barcodes: smallvec![sbc, cbc],
                    umi: packed((c * rpc + r) as u64 + 1, UMI_LEN),
                    dirs: vec![true],
                    refs: vec![refid],
                });
            }
            let chunk = Chunk::<MultiBarcodeReadRecord> {
                nbytes: 0,
                nrec: reads.len() as u32,
                reads,
            };
            fw.write_chunk(&chunk, &ctx)?;
        }
    }
    fw.finalize()?;

    // sample barcode list (shared by both variants; write once)
    let sbl = format!("{dir}/../sample_bc.txt");
    let mut w = std::io::BufWriter::new(File::create(&sbl)?);
    for s in &sample_bcs {
        writeln!(w, "{}", to_nuc(*s, SAMPLE_LEN as usize))?;
    }
    w.flush()?;
    eprintln!(
        "wrote {mode} multi RAD: {ns} samples x {cps} cells x {rpc} reads -> {dir}/map.rad (+ sample_bc.txt)"
    );
    Ok(())
}
