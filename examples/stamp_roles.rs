// Rewrite a legacy RAD's prelude into a versioned (spec-major-2) one, stamping a
// read-level tag with a Barcode collation role — for validating the role-driven
// generic collation path end-to-end. Chunk data is copied verbatim.
//
// Usage: stamp_roles <in.rad> <out.rad> <read_tag_name> <barcode_level> [ori_aln_tag]
use libradicl::header::RadPrelude;
use libradicl::rad_types::TagRole;
use std::io::{BufReader, Seek, SeekFrom, Write};

fn main() -> anyhow::Result<()> {
    let mut a = std::env::args().skip(1);
    let inp = a.next().expect("in.rad");
    let outp = a.next().expect("out.rad");
    let tag = a.next().expect("read_tag_name");
    let level: u8 = a.next().expect("barcode_level").parse()?;
    let ori_tag = a.next(); // optional alignment tag to stamp as Orientation

    let mut br = BufReader::new(std::fs::File::open(&inp)?);
    let mut prelude = RadPrelude::from_bytes(&mut br)?;
    // byte offset in the original file just past the tag *descriptor* sections
    // (i.e. where the file-tag values begin) — everything from here is copied.
    // BufReader::stream_position() already reports the logical (consumed) position.
    let desc_end = br.stream_position()?;

    // stamp: versioned header + the named read tag gets a Barcode role
    prelude.hdr.major_version = libradicl::constants::RAD_SPEC_MAJOR;
    prelude.hdr.minor_version = libradicl::constants::RAD_SPEC_MINOR;
    // Optionally rename the barcode tag (STAMP_RENAME_BC) so the file no longer
    // carries the `b` bridge name and is thus unknown to the fast engine —
    // exercising the role-driven auto-routing of unknown record types.
    let rename_bc = std::env::var("STAMP_RENAME_BC").ok();
    // Barcode nucleotide length carried by the role (STAMP_BC_LEN, default 16).
    let bc_len: u8 = std::env::var("STAMP_BC_LEN")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(16);
    let mut stamped = false;
    for t in &mut prelude.read_tags.tags {
        if t.name == tag {
            t.role = TagRole::Barcode {
                level,
                len: bc_len,
            };
            if let Some(new_name) = &rename_bc {
                t.name = new_name.clone();
            }
            stamped = true;
        }
    }
    anyhow::ensure!(stamped, "read tag `{tag}` not found");
    // Optionally stamp a Umi role on a read tag (STAMP_UMI=<name>) and rename it
    // (STAMP_RENAME_UMI=<new>), so the single-barcode context can be built from
    // roles alone (a Barcode role requires a matching Umi role).
    if let Ok(umi_name) = std::env::var("STAMP_UMI") {
        let rename_umi = std::env::var("STAMP_RENAME_UMI").ok();
        let mut ok = false;
        for t in &mut prelude.read_tags.tags {
            if t.name == umi_name {
                t.role = TagRole::Umi { len: 12 };
                if let Some(n) = &rename_umi {
                    t.name = n.clone();
                }
                ok = true;
            }
        }
        anyhow::ensure!(ok, "umi read tag `{umi_name}` not found");
    }
    if let Some(ot) = ori_tag {
        let mut ok = false;
        for t in &mut prelude.aln_tags.tags {
            if t.name == ot {
                t.role = TagRole::Orientation;
                ok = true;
            }
        }
        anyhow::ensure!(ok, "alignment tag `{ot}` not found");
    }

    let mut out = std::io::BufWriter::new(std::fs::File::create(&outp)?);
    prelude.write(&mut out)?; // v2 header (magic+version) + descriptors (with roles)

    // copy the remainder of the original file (file-tag values + all chunks)
    let mut rest = std::fs::File::open(&inp)?;
    rest.seek(SeekFrom::Start(desc_end))?;
    std::io::copy(&mut rest, &mut out)?;
    out.flush()?;
    eprintln!("stamped `{tag}` as Barcode{{level:{level}}} → {outp}");
    Ok(())
}
