/*
 * Copyright (c) 2020-2026 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! EM benchmarks.
//!
//! `em_optimize_subset` is the per-cell EM driver used by the `*-em` resolution
//! strategies; it runs the M-step (`em_update_subset`) to convergence. On the
//! project's toy dataset `quant -r cr-like-em --small-thresh 0` spends the bulk
//! of its ~50 CPU-seconds inside it, so it is the single hottest routine in the
//! quantification stage.
//!
//! Run with real data:
//! ```text
//! AF_BENCH_DATA=/path/to/fixtures cargo bench --bench em
//! ```

mod common;

use alevin_fry::em::{EmInitType, em_optimize_subset};
use alevin_fry::eq_class::IndexedEqList;
use common::CellData;

fn main() {
    divan::main();
}

/// The eq-class list and the per-cell workloads, loaded once for the whole run.
struct Workload {
    eql: IndexedEqList,
    cells: Vec<CellData>,
    num_alphas: usize,
}

fn workload(max_cells: usize) -> Workload {
    match (
        common::fixture("gene_eqclass.txt.gz"),
        common::fixture("geqc_counts.mtx"),
    ) {
        (Some(eqc), Some(counts)) => {
            let eql = IndexedEqList::init_from_eqc_file(&eqc);
            let cells = common::load_cell_data(&counts, max_cells);
            let num_alphas = eql.num_genes;
            Workload {
                eql,
                cells,
                num_alphas,
            }
        }
        _ => {
            common::note_synthetic("em");
            let (eql, cells) =
                common::synthetic_workload(20_000, 40_000, max_cells, 400, 0xA1E7_1F2E);
            let num_alphas = eql.num_genes;
            Workload {
                eql,
                cells,
                num_alphas,
            }
        }
    }
}

/// Full EM to convergence over a batch of cells, single-threaded.
///
/// `quant` runs one of these per cell across a thread pool; benchmarking the
/// batch serially removes the pool's scheduling noise while keeping the work
/// identical.
#[divan::bench(args = [16, 128])]
fn em_optimize_subset_batch(bencher: divan::Bencher, num_cells: usize) {
    let w = workload(num_cells);
    let log = common::null_logger();
    let cells = &w.cells[..w.cells.len().min(num_cells)];

    bencher
        .with_inputs(|| (vec![false; w.num_alphas], vec![true; w.num_alphas]))
        .bench_local_values(|(mut unique_evidence, mut no_ambiguity)| {
            let mut acc = 0.0f32;
            for cell in cells {
                let alphas = em_optimize_subset(
                    &w.eql,
                    cell,
                    &mut unique_evidence,
                    &mut no_ambiguity,
                    EmInitType::Informative,
                    w.num_alphas,
                    false,
                    None,
                    &log,
                );
                acc += alphas.iter().sum::<f32>();
            }
            acc
        });
}

/// The unique-counts-only path: no EM iteration, just the scatter that seeds
/// `alphas_in`. Isolates the cost of walking the CSR label list from the cost
/// of the iteration loop on top of it.
#[divan::bench(args = [128])]
fn em_unique_only(bencher: divan::Bencher, num_cells: usize) {
    let w = workload(num_cells);
    let log = common::null_logger();
    let cells = &w.cells[..w.cells.len().min(num_cells)];

    bencher
        .with_inputs(|| (vec![false; w.num_alphas], vec![true; w.num_alphas]))
        .bench_local_values(|(mut unique_evidence, mut no_ambiguity)| {
            let mut acc = 0.0f32;
            for cell in cells {
                let alphas = em_optimize_subset(
                    &w.eql,
                    cell,
                    &mut unique_evidence,
                    &mut no_ambiguity,
                    EmInitType::Informative,
                    w.num_alphas,
                    true,
                    None,
                    &log,
                );
                acc += alphas.iter().sum::<f32>();
            }
            acc
        });
}

/// Loading the dumped equivalence classes back into the flat CSR layout.
/// Only meaningful with real data; skipped otherwise.
#[divan::bench]
fn init_from_eqc_file(bencher: divan::Bencher) {
    match common::fixture("gene_eqclass.txt.gz") {
        Some(p) => {
            bencher.bench_local(|| IndexedEqList::init_from_eqc_file(&p));
        }
        None => {
            common::note_synthetic("em::init_from_eqc_file (skipped)");
        }
    }
}
