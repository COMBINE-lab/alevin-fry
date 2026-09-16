// Canonical dumper for a collated multi-barcode (Flex) RAD.
// Per chunk (one sample+cell group) prints: <sample> <cell> <nrec> <hash>.
// Sorted output ⇒ two collations that regroup the same records into the same
// (sample,cell) groups match regardless of chunk / within-cell order.
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::io::{BufReader, Cursor, Read};

use libradicl::codec::{chunk_codec_from_tag_map, decompress_payload};
use libradicl::header::RadPrelude;
use libradicl::record::{MappedRecord, MultiBarcodeReadRecordT, MultiBarcodeRecordContext};

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: dump_collated_multi <rad>");
    let f = std::fs::File::open(&path)?;
    let mut br = BufReader::new(f);

    let prelude = RadPrelude::from_bytes(&mut br)?;
    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    let codec = chunk_codec_from_tag_map(&file_tag_map)?;
    let ctx = prelude.get_record_context::<MultiBarcodeRecordContext>()?;
    let num_chunks = prelude.hdr.num_chunks as usize;

    let mut lines: Vec<(u64, u64, u32, u64)> = Vec::with_capacity(num_chunks);
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
        let (mut s0, mut c0): (u64, u64) = (0, 0);
        for i in 0..nrec {
            let r = MultiBarcodeReadRecordT::<u64>::from_bytes_with_context(&mut cur, &ctx);
            let sample = r.barcodes[0];
            let cell = *r.barcodes.last().unwrap();
            if i == 0 {
                s0 = sample;
                c0 = cell;
            }
            let mut pairs: Vec<(u32, bool)> = r
                .refs
                .iter()
                .cloned()
                .zip(r.dirs.iter().cloned().chain(std::iter::repeat(false)))
                .collect();
            pairs.sort();
            let mut h = DefaultHasher::new();
            sample.hash(&mut h);
            cell.hash(&mut h);
            r.umi.hash(&mut h);
            pairs.hash(&mut h);
            per_rec.push(h.finish());
        }
        per_rec.sort_unstable();
        let mut ch = DefaultHasher::new();
        per_rec.hash(&mut ch);
        lines.push((s0, c0, nrec, ch.finish()));
    }
    lines.sort_unstable();
    for (s, c, n, h) in &lines {
        println!("{s}\t{c}\t{n}\t{h:016x}");
    }
    eprintln!("chunks={} codec={:?}", lines.len(), codec);
    Ok(())
}
