// Canonical dumper for a collated RAD (standard short-read sc-RNA, u64 barcode).
// For each chunk (one cell) prints: <barcode> <nrec> <hash-of-sorted-records>.
// Lines are sorted, so two collations that regroup the SAME records into the
// SAME cells produce identical output regardless of chunk / within-cell order.
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::io::{BufReader, Cursor, Read};

use libradicl::codec::{chunk_codec_from_tag_map, decompress_payload};
use libradicl::header::RadPrelude;
use libradicl::record::{AlevinFryReadRecordT, AlevinFryRecordContext, MappedRecord};

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: dump_collated <rad>");
    let f = std::fs::File::open(&path)?;
    let mut br = BufReader::new(f);

    let prelude = RadPrelude::from_bytes(&mut br)?;
    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    let codec = chunk_codec_from_tag_map(&file_tag_map)?;
    let ctx = prelude.get_record_context::<AlevinFryRecordContext>()?;
    let num_chunks = prelude.hdr.num_chunks as usize;

    let mut lines: Vec<(u64, u32, u64)> = Vec::with_capacity(num_chunks);
    let mut hdr = [0u8; 8];
    for _ in 0..num_chunks {
        br.read_exact(&mut hdr)?;
        let disk_nbytes = u32::from_le_bytes([hdr[0], hdr[1], hdr[2], hdr[3]]) as usize;
        let nrec = u32::from_le_bytes([hdr[4], hdr[5], hdr[6], hdr[7]]);
        let mut payload = vec![0u8; disk_nbytes - 8];
        br.read_exact(&mut payload)?;
        let records = decompress_payload(codec, &payload)?;

        let mut cur = Cursor::new(records);
        let mut per_rec: Vec<u64> = Vec::with_capacity(nrec as usize);
        let mut bc0: u64 = 0;
        for i in 0..nrec {
            let r = AlevinFryReadRecordT::<u64>::from_bytes_with_context(&mut cur, &ctx);
            if i == 0 {
                bc0 = r.bc;
            }
            // canonicalize within-record: sort (ref,dir) pairs
            let mut pairs: Vec<(u32, bool)> = r
                .refs
                .iter()
                .cloned()
                .zip(r.dirs.iter().cloned().chain(std::iter::repeat(false)))
                .collect();
            pairs.sort();
            let mut h = DefaultHasher::new();
            r.bc.hash(&mut h);
            r.umi.hash(&mut h);
            pairs.hash(&mut h);
            per_rec.push(h.finish());
        }
        per_rec.sort_unstable();
        let mut ch = DefaultHasher::new();
        per_rec.hash(&mut ch);
        lines.push((bc0, nrec, ch.finish()));
    }
    lines.sort_unstable();
    for (bc, n, h) in &lines {
        println!("{bc}\t{n}\t{h:016x}");
    }
    eprintln!("chunks={} codec={:?}", lines.len(), codec);
    Ok(())
}
