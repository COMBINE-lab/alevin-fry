// Convert a divergent ScLong RAD (zzare's ScLongRead/develop-refactor format,
// which stores a fixed-width read-name block per record) into the current
// libradicl ScLong format by stripping that name block. The prelude/tag sections
// and every numeric field (na, bc, umi, and each [ori+ref_id, AS, start, end,
// tlen] alignment) are preserved verbatim; only the name block is removed and each
// chunk's nbytes recomputed.
//
// Usage: convert_sclong <in.rad> <out.rad> [name_width=69]
use libradicl::header::RadPrelude;
use libradicl::rad_types::RadType;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let in_path = args
        .next()
        .expect("usage: convert_sclong <in> <out> [name_width]");
    let out_path = args
        .next()
        .expect("usage: convert_sclong <in> <out> [name_width]");
    let name_width: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(69);
    let aln_stride = 20usize; // [ori+ref_id, AS, start, end, tlen] = 5 * u32

    let mut br = BufReader::with_capacity(1 << 20, std::fs::File::open(&in_path)?);
    let prelude = RadPrelude::from_bytes(&mut br)?;
    let _ftm = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    let num_chunks = prelude.hdr.num_chunks as usize;
    // header (na) + bc + umi bytes, from the read-level int tags
    let hdr_after_na: usize = prelude
        .read_tags
        .tags
        .iter()
        .filter_map(|t| match &t.typeid {
            RadType::Int(i) => Some(i.bytes_for_type()),
            _ => None,
        })
        .sum();
    let rec_hdr = 4 + hdr_after_na; // na + bc + umi

    // prelude+file-tags occupy [0, prelude_end); copy them verbatim.
    let prelude_end = br.stream_position()?;
    let mut header = vec![0u8; prelude_end as usize];
    {
        let mut f = std::fs::File::open(&in_path)?;
        f.read_exact(&mut header)?;
    }
    let mut bw = BufWriter::with_capacity(1 << 20, std::fs::File::create(&out_path)?);
    bw.write_all(&header)?;

    br.seek(SeekFrom::Start(prelude_end))?;
    let mut body = Vec::new();
    let mut out = Vec::new();
    let mut u32b = [0u8; 4];
    let mut total_recs = 0u64;
    for c in 0..num_chunks {
        br.read_exact(&mut u32b)?;
        let nbytes = u32::from_le_bytes(u32b) as usize;
        br.read_exact(&mut u32b)?;
        let nrec = u32::from_le_bytes(u32b);
        body.resize(nbytes - 8, 0);
        br.read_exact(&mut body)?;

        out.clear();
        let mut p = 0usize;
        for _ in 0..nrec {
            let na = u32::from_le_bytes(body[p..p + 4].try_into().unwrap()) as usize;
            let aln = na * aln_stride;
            // copy na + bc + umi
            out.extend_from_slice(&body[p..p + rec_hdr]);
            // skip the name block, copy the alignment bytes
            let aln_start = p + rec_hdr + name_width;
            out.extend_from_slice(&body[aln_start..aln_start + aln]);
            p = aln_start + aln;
            total_recs += 1;
        }
        assert_eq!(
            p,
            body.len(),
            "chunk {c}: record walk did not consume body exactly"
        );

        let new_nbytes = (out.len() + 8) as u32;
        bw.write_all(&new_nbytes.to_le_bytes())?;
        bw.write_all(&nrec.to_le_bytes())?;
        bw.write_all(&out)?;
    }
    bw.flush()?;
    eprintln!("converted {num_chunks} chunks, {total_recs} records → {out_path}");
    Ok(())
}
