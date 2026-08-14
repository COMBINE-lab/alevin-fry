/*
 * Copyright (c) 2020-2026 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Fixture loading shared by the benchmarks.
//!
//! The benchmarks prefer *real* data. Point `AF_BENCH_DATA` at a directory
//! holding the artifacts of an actual run and every bench below switches from
//! its synthetic fallback to the real thing:
//!
//! ```text
//! $AF_BENCH_DATA/
//!   gene_eqclass.txt.gz   # from `quant -r cr-like-em --dump-eqclasses`
//!   geqc_counts.mtx       # ditto
//!   permit_list.txt       # a barcode-per-line whitelist (e.g. 10x_v3_permit.txt)
//!   permit_map.bin        # from `generate-permit-list`
//!   map.collated.rad      # from `collate`
//! ```
//!
//! No fixture is committed to the repository: the smallest useful ones are
//! hundreds of megabytes. When a fixture is absent the bench falls back to a
//! deterministic synthetic generator and says so on stderr, so that a number
//! produced without real data is never mistaken for one produced with it.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// Root of the real-data fixtures, if the user provided one.
pub fn data_dir() -> Option<PathBuf> {
    let d = PathBuf::from(std::env::var_os("AF_BENCH_DATA")?);
    if d.is_dir() { Some(d) } else { None }
}

/// Path of `name` under the fixture root, if both exist.
pub fn fixture(name: &str) -> Option<PathBuf> {
    let p = data_dir()?.join(name);
    if p.is_file() { Some(p) } else { None }
}

/// Announce, once per bench binary, that we are running on synthetic input.
pub fn note_synthetic(what: &str) {
    use std::sync::Once;
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        eprintln!(
            "note: AF_BENCH_DATA is unset or incomplete; `{what}` runs on synthetic input. \
             Numbers are comparable between branches but not representative of real data."
        );
    });
}

/// A logger that drops everything; the EM routines want one but the benchmarks
/// do not want the I/O.
pub fn null_logger() -> slog::Logger {
    slog::Logger::root(slog::Discard, slog::o!())
}

/// xorshift64*: a tiny deterministic PRNG so the synthetic fallbacks are
/// reproducible without pulling the `rand` seeding machinery into the benches.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed | 1)
    }
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    pub fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

/// Per-cell equivalence-class workload: `(eq_class_id, umi_count)` pairs, the
/// exact shape `em_optimize_subset` consumes.
pub type CellData = Vec<(u32, u32)>;

/// Read `geqc_counts.mtx` (cells x eq-classes) into one `CellData` per cell,
/// keeping the `max_cells` cells with the most equivalence classes. EM cost is
/// dominated by the large cells, and taking the head of the file instead would
/// make the bench a measure of whichever cells happened to be written first.
pub fn load_cell_data(path: &Path, max_cells: usize) -> Vec<CellData> {
    let tri: sprs::TriMatI<f32, u32> =
        sprs::io::read_matrix_market(path).expect("could not read geqc_counts.mtx");
    let csr = tri.to_csr::<usize>();

    let mut cells: Vec<CellData> = csr
        .outer_iterator()
        .map(|row| {
            row.iter()
                .map(|(eqid, count)| (eqid as u32, *count as u32))
                .collect()
        })
        .filter(|c: &CellData| !c.is_empty())
        .collect();

    cells.sort_unstable_by_key(|c| std::cmp::Reverse(c.len()));
    cells.truncate(max_cells);
    cells
}

/// Build a synthetic `IndexedEqList` plus matching per-cell workloads.
///
/// This is a deterministic stress workload, not a model of any specific real
/// dataset: most equivalence classes carry a single gene label, a long tail
/// carries 2-8, and a handful carry up to 64.
pub fn synthetic_workload(
    num_genes: usize,
    num_eqc: usize,
    num_cells: usize,
    eqc_per_cell: usize,
    seed: u64,
) -> (alevin_fry::eq_class::IndexedEqList, Vec<CellData>) {
    let mut rng = Rng::new(seed);
    let mut eql = alevin_fry::eq_class::IndexedEqList::new();
    eql.eq_label_starts.push(0);
    eql.num_genes = num_genes;

    let mut labels: Vec<u32> = Vec::with_capacity(64);
    for _ in 0..num_eqc {
        let r = rng.below(100);
        let len = match r {
            0..=64 => 1,
            65..=94 => 2 + rng.below(7) as usize,
            _ => 9 + rng.below(56) as usize,
        };
        labels.clear();
        for _ in 0..len {
            labels.push(rng.below(num_genes as u64) as u32);
        }
        labels.sort_unstable();
        labels.dedup();
        eql.add_label_vec(&labels);
    }

    let cells = (0..num_cells)
        .map(|_| {
            let mut c: CellData = (0..eqc_per_cell)
                .map(|_| (rng.below(num_eqc as u64) as u32, 1 + rng.below(8) as u32))
                .collect();
            c.sort_unstable();
            c.dedup_by_key(|(e, _)| *e);
            c
        })
        .collect();

    (eql, cells)
}

/// Synthetic 2-bit packed barcodes of length `bc_len`, deterministic in `seed`.
pub fn synthetic_barcodes(n: usize, bc_len: usize, seed: u64) -> Vec<u64> {
    let mut rng = Rng::new(seed);
    let mask = if bc_len >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * bc_len)) - 1
    };
    let mut v: Vec<u64> = (0..n).map(|_| rng.next_u64() & mask).collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// Read a barcode-per-line whitelist and 2-bit pack it, keeping the first `n`.
pub fn load_barcodes(path: &Path, n: usize, bc_len: usize) -> Vec<u64> {
    use std::io::BufRead;
    let f = std::fs::File::open(path).expect("could not open permit list");
    let mut v = Vec::with_capacity(n);
    for line in std::io::BufReader::new(f).lines() {
        let line = line.expect("could not read permit list");
        let l = line.trim();
        if l.is_empty() {
            continue;
        }
        if l.len() != bc_len {
            continue;
        }
        if let Some((_, km, _)) =
            needletail::bitkmer::BitNuclKmer::new(l.as_bytes(), bc_len as u8, false).next()
        {
            v.push(km.0);
        }
        if v.len() == n {
            break;
        }
    }
    assert_eq!(
        v.len(),
        n,
        "permit_list.txt contains fewer than {n} valid {bc_len}-base barcodes"
    );
    v
}
