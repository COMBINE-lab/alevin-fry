/*
 * Copyright (c) 2020-2026 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Benchmarks for the two CPU costs `collate` pays on top of its I/O.
//!
//! `collate` writes each temporary bucket through a Snappy frame encoder and
//! reads it back the same way, and it deserializes the barcode-correction map
//! produced by `generate-permit-list` before it can start. Neither is the
//! dominant cost of the stage -- collate is I/O bound -- but both are the parts
//! a dependency swap would move, so they are what a dependency benchmark must
//! measure. The wall-clock effect on the stage as a whole has to be read off a
//! full pipeline run, not off these numbers.
//!
//! Run with real data:
//! ```text
//! AF_BENCH_DATA=/path/to/fixtures cargo bench --bench collate_io
//! ```

mod common;

use std::collections::HashMap;
use std::io::{Read, Write};

/// Size of one temporary-bucket write in `collate`. The buckets are sized from
/// the input chunk offsets; a few MiB is the common case on 10x data.
const CHUNK_SIZE: usize = 4 << 20;

fn main() {
    divan::main();
}

/// A buffer of collated-RAD bytes: real record bytes when a collated RAD file
/// is available, otherwise a synthetic stand-in with a similar entropy profile
/// (dense u32 reference ids, small counts, packed barcodes).
fn chunk_bytes() -> Vec<u8> {
    if let Some(p) = common::fixture("map.collated.rad") {
        let mut f = std::fs::File::open(p).expect("could not open collated rad");
        // skip the header; the tail of the file is pure record payload, which
        // is what the bucket writer actually compresses.
        let len = f.metadata().unwrap().len();
        let skip = (len / 2).min(len.saturating_sub(CHUNK_SIZE as u64));
        use std::io::Seek;
        f.seek(std::io::SeekFrom::Start(skip)).unwrap();
        let mut buf = vec![0u8; CHUNK_SIZE];
        let mut read = 0usize;
        while read < CHUNK_SIZE {
            match f.read(&mut buf[read..]) {
                Ok(0) => break,
                Ok(n) => read += n,
                Err(e) => panic!("could not read collated rad: {e}"),
            }
        }
        buf.truncate(read);
        return buf;
    }

    common::note_synthetic("collate_io");
    let mut rng = common::Rng::new(0xC0_11A7E);
    let mut buf = Vec::with_capacity(CHUNK_SIZE);
    while buf.len() < CHUNK_SIZE {
        buf.extend_from_slice(&(rng.below(375_446) as u32).to_le_bytes());
        buf.extend_from_slice(&(1u32 + rng.below(4) as u32).to_le_bytes());
        buf.extend_from_slice(&rng.next_u64().to_le_bytes());
    }
    buf.truncate(CHUNK_SIZE);
    buf
}

fn snap_compress(src: &[u8]) -> Vec<u8> {
    let mut enc = snap::write::FrameEncoder::new(Vec::with_capacity(src.len()));
    enc.write_all(src).unwrap();
    enc.into_inner().unwrap()
}

#[divan::bench]
fn snap_encode(bencher: divan::Bencher) {
    let src = chunk_bytes();
    bencher
        .counter(divan::counter::BytesCount::of_slice(&src))
        .bench_local(|| snap_compress(&src));
}

#[divan::bench]
fn snap_decode(bencher: divan::Bencher) {
    let src = chunk_bytes();
    let comp = snap_compress(&src);
    bencher
        .counter(divan::counter::BytesCount::of_slice(&src))
        .bench_local(|| {
            let mut out = Vec::with_capacity(src.len());
            snap::read::FrameDecoder::new(&comp[..])
                .read_to_end(&mut out)
                .unwrap();
            out
        });
}

/// Compression ratio, reported once so that a throughput win can be weighed
/// against the extra bytes it writes.
#[divan::bench]
fn snap_ratio(bencher: divan::Bencher) {
    let src = chunk_bytes();
    let comp = snap_compress(&src);
    eprintln!(
        "snappy: {} -> {} bytes ({:.3}x)",
        src.len(),
        comp.len(),
        src.len() as f64 / comp.len() as f64
    );
    bencher.bench_local(|| comp.len());
}

/// The barcode-correction map, at the size the toy dataset produces
/// (~80k entries). `collate` deserializes one of these per invocation, and
/// every worker thread clones the `Arc` to it.
fn correction_map(n: usize) -> HashMap<u64, u64> {
    let bcs = common::synthetic_barcodes(n, 16, 0x5EED_5EED);
    bcs.iter().map(|b| (*b, b ^ 0x3)).collect()
}

#[divan::bench(args = [80_000])]
fn bincode_serialize_correction_map(bencher: divan::Bencher, n: usize) {
    let m = correction_map(n);
    bencher.bench_local(|| bincode::serialize(&m).unwrap());
}

#[divan::bench(args = [80_000])]
fn bincode_deserialize_correction_map(bencher: divan::Bencher, n: usize) {
    let m = correction_map(n);
    let bytes = bincode::serialize(&m).unwrap();
    bencher
        .counter(divan::counter::BytesCount::of_slice(&bytes))
        .bench_local(|| bincode::deserialize::<HashMap<u64, u64>>(&bytes).unwrap());
}
