/*
 * Copyright (c) 2020-2026 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Barcode-correction benchmarks.
//!
//! `generate_permitlist_map` builds the map from every one-edit neighbour of a
//! retained barcode back to that barcode; it is the correction step of the
//! filtered `generate-permit-list` path and is almost entirely hash-table work
//! (roughly 10 insertions per retained barcode). `get_all_one_edit_neighbors`
//! is the pure bit-twiddling underneath it, with no hashing beyond the small
//! scratch set, which separates the two costs.
//!
//! Run with real data:
//! ```text
//! AF_BENCH_DATA=/path/to/fixtures cargo bench --bench barcode
//! ```

mod common;

use alevin_fry::utils::{
    generate_permitlist_map, generate_whitelist_set, get_all_one_edit_neighbors,
};

const BC_LEN: usize = 16; // 10x v3 cell barcode

fn main() {
    divan::main();
}

fn barcodes(n: usize) -> Vec<u64> {
    match common::fixture("permit_list.txt") {
        Some(p) => common::load_barcodes(&p, n),
        None => {
            common::note_synthetic("barcode");
            common::synthetic_barcodes(n, BC_LEN, 0xB0_1CE5)
        }
    }
}

/// Build the full one-edit correction map. This is the routine
/// `process_filtered` calls once per run over the retained barcodes.
#[divan::bench(args = [1_000, 10_000])]
fn permitlist_map(bencher: divan::Bencher, n: usize) {
    let bcs = barcodes(n);
    bencher.bench_local(|| generate_permitlist_map(&bcs, BC_LEN).unwrap());
}

/// The set-valued variant, which drops the correction target and keeps only
/// membership.
#[divan::bench(args = [10_000])]
fn whitelist_set(bencher: divan::Bencher, n: usize) {
    let bcs = barcodes(n);
    bencher.bench_local(|| generate_whitelist_set(&bcs, BC_LEN).unwrap());
}

/// Neighbour enumeration for a single barcode, reusing the scratch set exactly
/// as the callers do.
#[divan::bench]
fn one_edit_neighbors(bencher: divan::Bencher) {
    let bcs = barcodes(1_000);
    // left to inference so the bench keeps compiling whichever `BuildHasher`
    // `get_all_one_edit_neighbors` is declared with.
    let mut neighbors = Default::default();
    get_all_one_edit_neighbors(bcs[0], BC_LEN, &mut neighbors).unwrap();
    neighbors.reserve(3 * BC_LEN + 8 * (BC_LEN - 1));
    let mut i = 0usize;
    bencher.bench_local(|| {
        let bc = bcs[i % bcs.len()];
        i += 1;
        get_all_one_edit_neighbors(bc, BC_LEN, &mut neighbors).unwrap();
        neighbors.len()
    });
}
