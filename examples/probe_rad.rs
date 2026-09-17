// Probe a RAD's tag sections and reverse-engineer the per-alignment stride.
use libradicl::header::RadPrelude;
use libradicl::rad_types::RadIntId;
use std::io::{BufReader, Read, Seek, SeekFrom};

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: probe_rad <rad>");
    let mut br = BufReader::new(std::fs::File::open(&path)?);
    let prelude = RadPrelude::from_bytes(&mut br)?;
    let ftm = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    println!("num_chunks = {}", prelude.hdr.num_chunks);
    for (label, sec) in [
        ("FILE", &prelude.file_tags),
        ("READ", &prelude.read_tags),
        ("ALN", &prelude.aln_tags),
    ] {
        println!("--- {label} tags ---");
        for t in &sec.tags {
            println!("   {} : {:?}  role={:?}", t.name, t.typeid, t.role);
        }
    }
    println!("file tag values: {:?}", ftm);

    // position is now at first chunk start
    let mut u32b = [0u8; 4];
    br.read_exact(&mut u32b)?;
    let nbytes = u32::from_le_bytes(u32b);
    br.read_exact(&mut u32b)?;
    let nrec = u32::from_le_bytes(u32b);
    println!("chunk0: nbytes={nbytes} nrec={nrec}");

    // header widths: na(u32) + bc + umi. Guess bc/umi from read tag int types.
    let int_bytes = |t: &RadIntId| t.bytes_for_type();
    // read-level tags are typically [b (bc), u (umi)]
    let read_ints: Vec<usize> = prelude
        .read_tags
        .tags
        .iter()
        .filter_map(|t| match &t.typeid {
            libradicl::rad_types::RadType::Int(i) => Some(int_bytes(i)),
            _ => None,
        })
        .collect();
    println!("read-level int widths (bytes): {read_ints:?}");
    let hdr_after_na: usize = read_ints.iter().sum();
    println!(
        "record header = na(4) + {hdr_after_na} (bc+umi) = {}",
        4 + hdr_after_na
    );

    // Walk the chunk body brute-forcing aln stride. Body = nbytes - 8.
    let body = (nbytes as usize) - 8;
    let chunk_start = br.stream_position()?;

    // Byte-level dump of the first 220 bytes of the body: offset, hex, ascii.
    {
        let mut buf = [0u8; 220];
        br.read_exact(&mut buf)?;
        for (off, row) in buf.chunks(16).enumerate() {
            let hex: String = row.iter().map(|b| format!("{b:02x} ")).collect();
            let asc: String = row
                .iter()
                .map(|&b| {
                    if (0x20..0x7f).contains(&b) {
                        b as char
                    } else {
                        '.'
                    }
                })
                .collect();
            println!("  {:>4}: {:<48} {}", off * 16, hex, asc);
        }
        br.seek(SeekFrom::Start(chunk_start))?;
    }
    // Test: record = na(4) + bc+umi(hdr_after_na) + NAME(fixed) + na*aln_stride.
    let aln_stride = 20usize;
    for name in [0usize, 64, 68, 69, 72] {
        br.seek(SeekFrom::Start(chunk_start))?;
        let mut ok = true;
        let mut consumed = 0usize;
        let mut recs = 0usize;
        for _ in 0..nrec {
            let mut nb = [0u8; 4];
            if br.read_exact(&mut nb).is_err() {
                ok = false;
                break;
            }
            let na = u32::from_le_bytes(nb) as usize;
            let rec_bytes = 4 + hdr_after_na + name + na * aln_stride;
            let skip = rec_bytes - 4;
            if consumed + rec_bytes > body || br.seek(SeekFrom::Current(skip as i64)).is_err() {
                ok = false;
                break;
            }
            consumed += rec_bytes;
            recs += 1;
        }
        println!(
            "  name_field={name}: recs={recs}/{nrec} consumed={consumed} body={body} exact={}",
            ok && consumed == body && recs == nrec as usize
        );
    }
    Ok(())
}
