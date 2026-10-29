// coverage_diag.rs
//
// Dumps per-cell diagnostics for the coverage model to TSV files, so the
// effect of the coverage weights can be attributed in Python.
// Read-only: it looks at the resolved molecules and the per-transcript
// coverage state, and never changes what is quantified.
//
// Usage:
//   COVERAGE_DIAG_DIR=diag_out target/release/alevin-fry quant ... --coverage-model
//   COVERAGE_DIAG_STRIDE=10   optional; write the per-transcript rows only for
//                             cells with cell_num % 10 == 0 (default 1 = all)
//
// Output files (in COVERAGE_DIAG_DIR):
//   coverage_diag_cells.tsv         one row per cell (every cell)
//   coverage_diag_txps.tsv.gz       one row per (cell, candidate transcript)
//   coverage_diag_profiles.tsv.gz   one row per (cell, eligible transcript):
//                                   binned coverage and per-bin probability
//
// A transcript is "eligible" in a cell when the model actually fits a coverage
// profile for it, i.e. total_weight > MIN_ALIGNMENTS_FOR_COVERAGE; otherwise
// every bin has probability 1.0 and the transcript is coverage-neutral.
//
// Molecule categories (cells file):
//   single        one candidate transcript; its weight is 1.0 regardless
//   multi_neutral >1 candidates, none eligible; all coverage weights are 1.0
//   multi_same    >1 candidates, >=1 eligible, but no weight moved
//   changed       the coverage-adjusted weights differ from the base weights
//     changed_mixed     ... candidates are a mix of eligible and neutral
//     changed_all_elig  ... every candidate is eligible
//
// The .gz files are written as one gzip member per cell, so they stay valid
// if the run is interrupted and need no finalisation step.

use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::{Mutex, OnceLock};

use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::bitkmer::bitmer_to_bytes;

use crate::eq_class::{merge_intervals, CellTxpCoverage, MIN_ALIGNMENTS_FOR_COVERAGE};
use crate::pugutils_dp_new::ResolvedMolecule;

// a molecule weight that moves by less than this is treated as unchanged
const WEIGHT_EPS: f64 = 1e-6;

pub struct CoverageDiag {
    stride: usize,
    cells: Mutex<File>,
    txps: Mutex<File>,
    profiles: Mutex<File>,
}

#[derive(Default)]
struct TxpStat {
    n_mol: u32,
    n_mol_multi: u32,
    n_mol_multi_elig: u32,
    n_mol_changed: u32,
    w_base: f64,
    w_cov: f64,
    cov_w_sum: f64,
    max_abs_dw: f64,
}

#[derive(Default)]
struct CellStat {
    n_mol: u64,
    n_mol_single: u64,
    n_mol_multi_neutral: u64,
    n_mol_multi_same: u64,
    n_mol_changed: u64,
    n_mol_changed_mixed: u64,
    n_mol_changed_all_elig: u64,
    moved_weight: f64,
    moved_weight_mixed: f64,
    moved_weight_all_elig: f64,
    n_txp: u64,
    n_txp_eligible: u64,
    n_txp_only_single: u64,
    n_txp_neutral: u64,
    n_txp_same: u64,
    n_txp_changed: u64,
}

fn gz_member(text: &str) -> Vec<u8> {
    let mut enc = GzEncoder::new(Vec::new(), Compression::fast());
    enc.write_all(text.as_bytes()).expect("gzip write failed");
    enc.finish().expect("gzip finish failed")
}

fn is_eligible(txp: &CellTxpCoverage) -> bool {
    !txp.coverage_bins.is_empty() && txp.total_weight > MIN_ALIGNMENTS_FOR_COVERAGE
}

fn join_floats(vals: &[f64], precision: usize) -> String {
    vals.iter()
        .map(|v| format!("{:.*}", precision, v))
        .collect::<Vec<_>>()
        .join(",")
}

fn init() -> Option<CoverageDiag> {
    let dir = std::env::var("COVERAGE_DIAG_DIR").ok()?;
    let stride = std::env::var("COVERAGE_DIAG_STRIDE")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|v| *v > 0)
        .unwrap_or(1);

    let dir = Path::new(&dir);
    std::fs::create_dir_all(dir).expect("could not create COVERAGE_DIAG_DIR");
    let create = |name: &str| {
        File::create(dir.join(name)).unwrap_or_else(|e| panic!("could not create {name}: {e}"))
    };

    let mut cells = create("coverage_diag_cells.tsv");
    writeln!(
        cells,
        "cell_num\tbarcode\tn_mol\tn_mol_single\tn_mol_multi_neutral\tn_mol_multi_same\t\
         n_mol_changed\tn_mol_changed_mixed\tn_mol_changed_all_elig\tmoved_weight\t\
         moved_weight_mixed\tmoved_weight_all_elig\tn_txp\tn_txp_eligible\t\
         n_txp_only_single\tn_txp_neutral\tn_txp_same\tn_txp_changed"
    )
    .unwrap();

    let mut txps = create("coverage_diag_txps.tsv.gz");
    txps.write_all(&gz_member(
        "cell_num\tbarcode\ttxp_id\tcol_id\ttlen\teligible\tcov_n_mol\tn_mol\tn_mol_multi\t\
         n_mol_multi_elig\tn_mol_changed\tw_base\tw_cov\tmax_abs_dw\tcov_w_mean\t\
         cov_prob_min\tcov_prob_max\tcov_prob_mean\n",
    ))
    .unwrap();

    let mut profiles = create("coverage_diag_profiles.tsv.gz");
    profiles
        .write_all(&gz_member(
            "cell_num\tbarcode\ttxp_id\tcol_id\ttlen\tbin_width\tcov_n_mol\tcov_bins\tcov_probs\n",
        ))
        .unwrap();

    eprintln!(
        "coverage diagnostics -> {} (per-transcript rows every {} cell(s))",
        dir.display(),
        stride
    );
    Some(CoverageDiag {
        stride,
        cells: Mutex::new(cells),
        txps: Mutex::new(txps),
        profiles: Mutex::new(profiles),
    })
}

/// The diagnostics sink, or `None` when COVERAGE_DIAG_DIR is not set.
pub fn get() -> Option<&'static CoverageDiag> {
    static DIAG: OnceLock<Option<CoverageDiag>> = OnceLock::new();
    DIAG.get_or_init(init).as_ref()
}

impl CoverageDiag {
    /// Call once per cell, after `compute_cell_coverage_probs`.
    ///
    /// `coverage[t].coverage_bins` already includes the `total_weight / 100`
    /// floor added by `compute_cell_coverage_probs`.
    pub fn dump_cell(
        &self,
        cell_num: usize,
        bc: u64,
        barcode_len: u16,
        molecules: &[ResolvedMolecule],
        coverage: &[CellTxpCoverage],
        tid_to_gid: &[u32],
    ) {
        let bc_bytes = bitmer_to_bytes((bc, barcode_len as u8));
        let barcode = String::from_utf8_lossy(&bc_bytes);

        let mut cell = CellStat::default();
        let mut txp_stats: HashMap<u32, TxpStat> = HashMap::new();
        let mut cov_w: Vec<f64> = Vec::new();

        for molecule in molecules {
            let cands = &molecule.transcript_probs;
            if cands.is_empty() {
                continue;
            }

            // NOTE: mirrors the weight computation of the `--coverage-model`
            // branch in quant.rs; keep the two in sync.
            cov_w.clear();
            cov_w.extend(cands.iter().map(|&(transcript_id, _)| {
                molecule
                    .intervals
                    .get(&transcript_id)
                    .map(|intervals| {
                        let merged =
                            merge_intervals(intervals.iter().map(|i| (i.start, i.end)).collect());
                        coverage[transcript_id as usize].coverage_probability_for_intervals(&merged)
                    })
                    .unwrap_or(1.0)
            }));
            let total: f64 = cands.iter().zip(cov_w.iter()).map(|(&(_, w), c)| w * c).sum();
            let weights_ok = total.is_finite() && total > 0.0;

            let n_elig = cands
                .iter()
                .filter(|&&(t, _)| is_eligible(&coverage[t as usize]))
                .count();
            let multi = cands.len() > 1;

            let mut max_dw = 0.0_f64;
            let mut moved = 0.0_f64;
            for (i, &(transcript_id, base_w)) in cands.iter().enumerate() {
                let new_w = if weights_ok { base_w * cov_w[i] / total } else { base_w };
                let dw = (new_w - base_w).abs();
                max_dw = max_dw.max(dw);
                moved += dw;

                let s = txp_stats.entry(transcript_id).or_default();
                s.n_mol += 1;
                if multi {
                    s.n_mol_multi += 1;
                    if n_elig > 0 {
                        s.n_mol_multi_elig += 1;
                    }
                }
                if dw > WEIGHT_EPS {
                    s.n_mol_changed += 1;
                }
                s.w_base += base_w;
                s.w_cov += new_w;
                s.cov_w_sum += cov_w[i];
                s.max_abs_dw = s.max_abs_dw.max(dw);
            }
            // half the L1 distance = fraction of this molecule that is reassigned
            let moved = 0.5 * moved;

            cell.n_mol += 1;
            if !multi {
                cell.n_mol_single += 1;
            } else if n_elig == 0 {
                cell.n_mol_multi_neutral += 1;
            } else if max_dw <= WEIGHT_EPS {
                cell.n_mol_multi_same += 1;
            } else {
                cell.n_mol_changed += 1;
                cell.moved_weight += moved;
                if n_elig == cands.len() {
                    cell.n_mol_changed_all_elig += 1;
                    cell.moved_weight_all_elig += moved;
                } else {
                    cell.n_mol_changed_mixed += 1;
                    cell.moved_weight_mixed += moved;
                }
            }
        }

        let detail = cell_num % self.stride == 0;
        let mut txp_rows = String::new();
        let mut profile_rows = String::new();

        let mut txp_ids: Vec<u32> = txp_stats.keys().copied().collect();
        txp_ids.sort_unstable();
        for transcript_id in txp_ids {
            let s = &txp_stats[&transcript_id];
            let txp = &coverage[transcript_id as usize];
            let eligible = is_eligible(txp);

            cell.n_txp += 1;
            if eligible {
                cell.n_txp_eligible += 1;
            }
            if s.n_mol_multi == 0 {
                cell.n_txp_only_single += 1;
            } else if s.n_mol_multi_elig == 0 {
                cell.n_txp_neutral += 1;
            } else if s.n_mol_changed == 0 {
                cell.n_txp_same += 1;
            } else {
                cell.n_txp_changed += 1;
            }

            if !detail {
                continue;
            }

            let probs = &txp.coverage_prob;
            let (p_min, p_max, p_mean) = if probs.is_empty() {
                (1.0, 1.0, 1.0)
            } else {
                (
                    probs.iter().copied().fold(f64::INFINITY, f64::min),
                    probs.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                    probs.iter().sum::<f64>() / probs.len() as f64,
                )
            };
            let col_id = tid_to_gid[transcript_id as usize];
            txp_rows.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\n",
                cell_num,
                barcode,
                transcript_id,
                col_id,
                txp.len,
                eligible as u8,
                txp.total_weight,
                s.n_mol,
                s.n_mol_multi,
                s.n_mol_multi_elig,
                s.n_mol_changed,
                s.w_base,
                s.w_cov,
                s.max_abs_dw,
                s.cov_w_sum / s.n_mol as f64,
                p_min,
                p_max,
                p_mean,
            ));
            if eligible {
                profile_rows.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    cell_num,
                    barcode,
                    transcript_id,
                    col_id,
                    txp.len,
                    txp.bin_width,
                    txp.total_weight,
                    join_floats(&txp.coverage_bins, 3),
                    join_floats(probs, 4),
                ));
            }
        }

        let cell_row = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            cell_num,
            barcode,
            cell.n_mol,
            cell.n_mol_single,
            cell.n_mol_multi_neutral,
            cell.n_mol_multi_same,
            cell.n_mol_changed,
            cell.n_mol_changed_mixed,
            cell.n_mol_changed_all_elig,
            cell.moved_weight,
            cell.moved_weight_mixed,
            cell.moved_weight_all_elig,
            cell.n_txp,
            cell.n_txp_eligible,
            cell.n_txp_only_single,
            cell.n_txp_neutral,
            cell.n_txp_same,
            cell.n_txp_changed,
        );
        self.cells
            .lock()
            .unwrap()
            .write_all(cell_row.as_bytes())
            .expect("can't write coverage_diag_cells.tsv");

        if !txp_rows.is_empty() {
            let member = gz_member(&txp_rows);
            self.txps
                .lock()
                .unwrap()
                .write_all(&member)
                .expect("can't write coverage_diag_txps.tsv.gz");
        }
        if !profile_rows.is_empty() {
            let member = gz_member(&profile_rows);
            self.profiles
                .lock()
                .unwrap()
                .write_all(&member)
                .expect("can't write coverage_diag_profiles.tsv.gz");
        }
    }
}
