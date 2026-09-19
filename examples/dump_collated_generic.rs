// Record-type-agnostic canonical dumper for a collated RAD with a fixed-width
// tag layout (any of the fixed record types: basic, position, multi-barcode…).
// Groups records by the collation key (the Barcode-role read tag, or the first
// read tag if none is annotated) and, per group, hashes the sorted set of raw
// record byte-strings. Sorted output ⇒ two collations that regroup the same
// records into the same cells match regardless of chunk / within-cell order.
//
// Usage: dump_collated_generic <rad>
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::io::{BufReader, Read};

use libradicl::codec::{chunk_codec_from_tag_map, decompress_payload};
use libradicl::header::RadPrelude;
use libradicl::rad_types::{RadType, TagRole};

fn int_bytes(t: &RadType) -> Option<usize> {
    match t {
        RadType::Int(i) => Some(i.bytes_for_type()),
        _ => None,
    }
}

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: dump_collated_generic <rad>");
    let mut br = BufReader::new(std::fs::File::open(&path)?);
    let prelude = RadPrelude::from_bytes(&mut br)?;
    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    let codec = chunk_codec_from_tag_map(&file_tag_map)?;

    // Fixed read header (na + read tags) + key field offset/width.
    let mut read_off = vec![]; // (offset, width) per read tag
    let mut off = 4usize; // past na
    for td in &prelude.read_tags.tags {
        let w = int_bytes(&td.typeid).expect("read tag must be fixed-width int");
        read_off.push((off, w));
        off += w;
    }
    let read_header = off;
    // key = first Barcode-role read tag, else first read tag.
    let key_idx = prelude
        .read_tags
        .tags
        .iter()
        .position(|t| matches!(t.role, TagRole::Barcode { .. }))
        .unwrap_or(0);
    let (key_off, key_w) = read_off[key_idx];
    let aln_stride: usize = prelude
        .aln_tags
        .tags
        .iter()
        .map(|t| int_bytes(&t.typeid).expect("aln tag must be fixed-width int"))
        .sum();

    let read_le = |b: &[u8], o: usize, w: usize| -> u128 {
        let mut v = 0u128;
        for i in 0..w {
            v |= (b[o + i] as u128) << (8 * i);
        }
        v
    };

    let num_chunks = prelude.hdr.num_chunks as usize;
    let mut lines: Vec<(u128, u32, u64)> = Vec::with_capacity(num_chunks);
    let mut hdr = [0u8; 8];
    for _ in 0..num_chunks {
        br.read_exact(&mut hdr)?;
        let disk_nbytes = u32::from_le_bytes(hdr[0..4].try_into().unwrap()) as usize;
        let nrec = u32::from_le_bytes(hdr[4..8].try_into().unwrap());
        let mut payload = vec![0u8; disk_nbytes - 8];
        br.read_exact(&mut payload)?;
        let recs = decompress_payload(codec, &payload)?;

        let mut pos = 0usize;
        let mut key0 = 0u128;
        let mut per_rec: Vec<u64> = Vec::with_capacity(nrec as usize);
        for i in 0..nrec {
            let na = u32::from_le_bytes(recs[pos..pos + 4].try_into().unwrap()) as usize;
            let rec_len = read_header + na * aln_stride;
            let rec = &recs[pos..pos + rec_len];
            if i == 0 {
                key0 = read_le(rec, key_off, key_w);
            }
            let mut h = DefaultHasher::new();
            rec.hash(&mut h);
            per_rec.push(h.finish());
            pos += rec_len;
        }
        per_rec.sort_unstable();
        let mut ch = DefaultHasher::new();
        per_rec.hash(&mut ch);
        lines.push((key0, nrec, ch.finish()));
    }
    lines.sort_unstable();
    for (k, n, h) in &lines {
        println!("{k}\t{n}\t{h:016x}");
    }
    eprintln!("chunks={} codec={:?}", lines.len(), codec);
    Ok(())
}
