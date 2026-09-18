// Rewrite a legacy multi-barcode (Flex) RAD's prelude into a versioned
// (spec-major-2) one, stamping the two barcode read tags with Barcode roles
// (levels 0/1), the UMI tag with the Umi role, and the orientation alignment
// field with the Orientation role -- and renaming them away from the b0/b1/u
// bridge so the collate dispatch must use the declared roles. Chunk data is
// copied verbatim. For validating the role-driven composite collation path.
//
// Usage: stamp_multi_roles <in.rad> <out.rad>
//   (assumes read tags b0, b1, u and alignment tag compressed_ori_refid)
use libradicl::header::RadPrelude;
use libradicl::rad_types::TagRole;
use std::io::{BufReader, Seek, SeekFrom, Write};

fn main() -> anyhow::Result<()> {
    let mut a = std::env::args().skip(1);
    let inp = a.next().expect("in.rad");
    let outp = a.next().expect("out.rad");

    let mut br = BufReader::new(std::fs::File::open(&inp)?);
    let mut prelude = RadPrelude::from_bytes(&mut br)?;
    let desc_end = br.stream_position()?;

    prelude.hdr.major_version = libradicl::constants::RAD_SPEC_MAJOR;
    prelude.hdr.minor_version = libradicl::constants::RAD_SPEC_MINOR;

    for t in &mut prelude.read_tags.tags {
        // Carry the nucleotide length in the role (10x Flex: 8bp sample, 16bp cell),
        // so a role-only RAD needs no b0len/b1len file tags.
        match t.name.as_str() {
            "b0" => {
                t.role = TagRole::Barcode { level: 0, len: 8 };
                t.name = "sample_bc".to_string();
            }
            "b1" => {
                t.role = TagRole::Barcode { level: 1, len: 16 };
                t.name = "cell_bc".to_string();
            }
            "u" => {
                t.role = TagRole::Umi;
                t.name = "umi".to_string();
            }
            _ => {}
        }
    }
    for t in &mut prelude.aln_tags.tags {
        if t.name == "compressed_ori_refid" {
            t.role = TagRole::Orientation;
        }
    }

    let mut out = std::io::BufWriter::new(std::fs::File::create(&outp)?);
    prelude.write(&mut out)?;
    let mut rest = std::fs::File::open(&inp)?;
    rest.seek(SeekFrom::Start(desc_end))?;
    std::io::copy(&mut rest, &mut out)?;
    out.flush()?;
    eprintln!("stamped multi-barcode roles (sample_bc/cell_bc/umi + Orientation) -> {outp}");
    Ok(())
}
