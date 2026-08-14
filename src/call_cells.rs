/*
 * Copyright (c) 2020-2026 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Cell calling by order of magnitude plus EmptyDrops, following Cell Ranger.
//!
//! This is a port of the two stages Cell Ranger runs after counting, from
//! `lib/python/cellranger/cell_calling_helpers.py` and
//! `lib/python/cellranger/cell_calling.py` at tag `cellranger-10.0.0`:
//!
//! 1. **ordmag**: keep every barcode within an order of magnitude of a
//!    high-count barcode. The threshold is bootstrapped, and the number of
//!    recovered cells is itself estimated by searching for a fixed point of
//!    that procedure.
//! 2. **EmptyDrops**: estimate an ambient expression profile from barcodes in
//!    a known-empty count range, smoothed by Simple Good-Turing, then test
//!    every remaining candidate barcode against a multinomial draw from that
//!    profile by Monte Carlo, and keep the ones that are significantly
//!    unlikely to be ambient after Benjamini-Hochberg correction.
//!
//! It lives in its own subcommand rather than in `generate-permit-list` for a
//! structural reason: EmptyDrops needs a barcode-by-feature count matrix, and
//! `generate-permit-list` runs before `collate` and `quant`, when only per
//! barcode *read* counts exist. Cell Ranger calls cells after counting, so this
//! does too.

use anyhow::{Context, bail};
use num_format::{Locale, ToFormattedString};
use serde_json::json;
use slog::info;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::prog_opts::CallCellsOpts;

/// Number of bootstrap samples in the ordmag stage
/// (`ORDMAG_NUM_BOOTSTRAP_SAMPLES`).
const ORDMAG_NUM_BOOTSTRAP_SAMPLES: usize = 100;
/// Quantile used to turn a recovered-cell count into a baseline index
/// (`ORDMAG_RECOVERED_CELLS_QUANTILE`).
const ORDMAG_RECOVERED_CELLS_QUANTILE: f64 = 0.99;
/// Floor on the estimated number of recovered cells
/// (`MIN_RECOVERED_CELLS_PER_GEM_GROUP`).
const MIN_RECOVERED_CELLS_PER_GEM_GROUP: usize = 50;
/// Ceiling on the estimated number of recovered cells
/// (`MAX_RECOVERED_CELLS_PER_GEM_GROUP`).
const MAX_RECOVERED_CELLS_PER_GEM_GROUP: usize = 1 << 18;
/// Threshold on the Turing/linear-Good-Turing switch, from the original
/// Simple Good-Turing paper and reproduced verbatim by Cell Ranger.
const SGT_SWITCH_THRESHOLD: f64 = 1.65;

/// A barcode-by-feature count matrix, held one barcode at a time.
struct BarcodeMatrix {
    /// Barcode sequences, in matrix row order.
    barcodes: Vec<String>,
    /// `(feature, count)` for each barcode.
    rows: Vec<Vec<(u32, f64)>>,
    /// Number of features (matrix columns).
    num_features: usize,
}

impl BarcodeMatrix {
    /// Total counts for each barcode.
    fn counts_per_barcode(&self) -> Vec<f64> {
        self.rows
            .iter()
            .map(|r| r.iter().map(|&(_, c)| c).sum())
            .collect()
    }
}

/// Read `quants_mat.mtx` and `quants_mat_rows.txt` from a `quant` output
/// directory.
fn read_matrix(quant_dir: &Path) -> anyhow::Result<BarcodeMatrix> {
    let alevin_dir = quant_dir.join("alevin");
    let mtx_path = alevin_dir.join("quants_mat.mtx");
    let rows_path = alevin_dir.join("quants_mat_rows.txt");

    let barcodes: Vec<String> = BufReader::new(
        std::fs::File::open(&rows_path)
            .with_context(|| format!("could not open {}", rows_path.display()))?,
    )
    .lines()
    .collect::<Result<Vec<_>, _>>()?
    .into_iter()
    .map(|l| l.trim().to_owned())
    .filter(|l| !l.is_empty())
    .collect();

    let f = std::fs::File::open(&mtx_path)
        .with_context(|| format!("could not open {}", mtx_path.display()))?;
    let mut lines = BufReader::new(f).lines();

    // banner, then any number of comment lines, then the dimension line
    let mut header = lines
        .next()
        .transpose()?
        .context("matrix market file was empty")?;
    while header.starts_with('%') {
        header = lines
            .next()
            .transpose()?
            .context("matrix market file had no dimension line")?;
    }
    let dims: Vec<usize> = header
        .split_whitespace()
        .map(|t| t.parse::<usize>())
        .collect::<Result<Vec<_>, _>>()
        .context("could not parse the matrix market dimension line")?;
    if dims.len() != 3 {
        bail!("expected 3 fields on the matrix market dimension line, found {dims:?}");
    }
    let (num_barcodes, num_features) = (dims[0], dims[1]);

    if num_barcodes != barcodes.len() {
        bail!(
            "the matrix has {} rows but quants_mat_rows.txt lists {} barcodes",
            num_barcodes,
            barcodes.len()
        );
    }

    let mut rows: Vec<Vec<(u32, f64)>> = vec![Vec::new(); num_barcodes];
    for line in lines {
        let line = line?;
        let mut it = line.split_whitespace();
        let (Some(i), Some(j), Some(v)) = (it.next(), it.next(), it.next()) else {
            continue;
        };
        let i: usize = i.parse()?;
        let j: u32 = j.parse()?;
        let v: f64 = v.parse()?;
        if v != 0.0 {
            rows[i - 1].push((j - 1, v));
        }
    }

    Ok(BarcodeMatrix {
        barcodes,
        rows,
        num_features,
    })
}

/// xorshift64*, so the bootstrap and the Monte Carlo are reproducible without
/// depending on how any particular `rand` version happens to seed itself.
/// Cell Ranger pins its own seed for the same reason (`np.random.RandomState(0)`).
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed | 1)
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    /// Uniform on [0, 1).
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

// ---------------------------------------------------------------------------
// Simple Good-Turing
// ---------------------------------------------------------------------------

/// `_averaging_transform`: divide each frequency-of-frequency by the width of
/// the gap around it.
fn averaging_transform(r: &[f64], nr: &[f64]) -> Vec<f64> {
    let n = r.len();
    let mut d = Vec::with_capacity(n);
    d.push(1.0);
    for i in 1..n {
        d.push(r[i] - r[i - 1]);
    }
    let mut dr = Vec::with_capacity(n);
    for i in 0..n - 1 {
        dr.push(0.5 * (d[i + 1] + d[i]));
    }
    dr.push(d[n - 1]);
    (0..n).map(|i| nr[i] / dr[i]).collect()
}

/// Least-squares slope of `y` on `x`, which is all Cell Ranger uses from
/// `scipy.stats.linregress`.
fn linregress_slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..x.len() {
        num += (x[i] - mx) * (y[i] - my);
        den += (x[i] - mx) * (x[i] - mx);
    }
    if den == 0.0 { 0.0 } else { num / den }
}

/// `simple_good_turing`: smoothed frequencies and the mass to assign to unseen
/// items. `xr` are the distinct non-zero frequencies, `xnr` their counts.
fn simple_good_turing(xr: &[f64], xnr: &[f64]) -> anyhow::Result<(Vec<f64>, f64)> {
    let n = xr.len();
    let x_total: f64 = (0..n).map(|i| xr[i] * xnr[i]).sum();

    let xnrz = averaging_transform(xr, xnr);
    let logx: Vec<f64> = xr.iter().map(|v| v.ln()).collect();
    let logz: Vec<f64> = xnrz.iter().map(|v| v.ln()).collect();
    let mut slope = linregress_slope(&logx, &logz);
    if slope >= -1.0 {
        // Cell Ranger clamps rather than failing, and says so on stdout.
        slope = -1.0;
    }

    // linear Good-Turing estimate
    let xrst: Vec<f64> = xr
        .iter()
        .map(|&r| r * (1.0 + 1.0 / r).powf(1.0 + slope))
        .collect();
    let xrstrel: Vec<f64> = (0..n).map(|i| xrst[i] / xr[i]).collect();

    // Turing estimate, defined only where r+1 is also observed
    let xrtry: Vec<bool> = (0..n)
        .map(|i| {
            let next = if i + 1 < n { xr[i + 1] - 1.0 } else { 0.0 };
            xr[i] == next
        })
        .collect();
    let mut xrstarel = vec![0.0f64; n];
    for i in 0..n {
        if xrtry[i] {
            xrstarel[i] = (xr[i] + 1.0) / xr[i] * xnr[i + 1] / xnr[i];
        }
    }

    let mut tursd = vec![1.0f64; n];
    for i in 0..n {
        if xrtry[i] {
            tursd[i] =
                (i as f64 + 2.0) / xnr[i] * (xnr[i + 1] * (1.0 + xnr[i + 1] / xnr[i])).sqrt();
        }
    }

    // Switch from the Turing estimate to the linear one at the first r where
    // they agree to within 1.65 standard deviations, and never switch back.
    let mut xrstcmbrel = vec![0.0f64; n];
    let mut use_turing = true;
    for r in 0..n {
        if !use_turing {
            xrstcmbrel[r] = xrstrel[r];
        } else if (xrstrel[r] - xrstarel[r]).abs() * (1.0 + r as f64) / tursd[r]
            > SGT_SWITCH_THRESHOLD
        {
            xrstcmbrel[r] = xrstarel[r];
        } else {
            use_turing = false;
            xrstcmbrel[r] = xrstrel[r];
        }
    }

    let sumpraw: f64 = (0..n)
        .map(|i| xrstcmbrel[i] * xr[i] * xnr[i] / x_total)
        .sum();
    if sumpraw == 0.0 {
        bail!("Simple Good-Turing produced a zero total probability");
    }
    let p0 = xnr[0] / x_total;
    for v in xrstcmbrel.iter_mut() {
        *v *= (1.0 - p0) / sumpraw;
    }

    Ok(((0..n).map(|i| xr[i] * xrstcmbrel[i]).collect(), p0))
}

/// `sgt_proportions`: adjusted proportions for the observed items, and the
/// total probability left for unobserved ones.
fn sgt_proportions(frequencies: &[f64]) -> anyhow::Result<(Vec<f64>, f64)> {
    if frequencies.is_empty() {
        bail!("Simple Good-Turing was given an empty frequency vector");
    }
    let max_f = frequencies.iter().cloned().fold(0.0f64, f64::max) as usize;
    let mut freqfreqs = vec![0.0f64; max_f + 1];
    for &f in frequencies {
        if f <= 0.0 {
            bail!("Simple Good-Turing requires strictly positive frequencies");
        }
        freqfreqs[f as usize] += 1.0;
    }
    let use_freqs: Vec<usize> = (1..=max_f).filter(|&i| freqfreqs[i] > 0.0).collect();
    if use_freqs.len() < 10 {
        bail!(
            "too few distinct non-zero frequencies ({}) to run Simple Good-Turing; \
             the ambient barcode range is probably empty or nearly so",
            use_freqs.len()
        );
    }

    let xr: Vec<f64> = use_freqs.iter().map(|&i| i as f64).collect();
    let xnr: Vec<f64> = use_freqs.iter().map(|&i| freqfreqs[i]).collect();
    let (rstar, p0) = simple_good_turing(&xr, &xnr)?;

    let rstar_sum: f64 = (0..xr.len()).map(|i| xnr[i] * rstar[i]).sum();
    let mut rstar_of = vec![0.0f64; max_f + 1];
    for (i, &f) in use_freqs.iter().enumerate() {
        rstar_of[f] = rstar[i];
    }
    let pstar: Vec<f64> = frequencies
        .iter()
        .map(|&f| (1.0 - p0) * (rstar_of[f as usize] / rstar_sum))
        .collect();
    Ok((pstar, p0))
}

// ---------------------------------------------------------------------------
// ordmag
// ---------------------------------------------------------------------------

/// `find_within_ordmag`: number of barcodes whose count is at least a tenth of
/// the count at `baseline_idx` places from the top.
fn find_within_ordmag(sorted_ascending: &[f64], baseline_idx: usize) -> usize {
    let n = sorted_ascending.len();
    if n == 0 {
        return 0;
    }
    let baseline = sorted_ascending[n - 1 - baseline_idx.min(n - 1)];
    let cutoff = (0.1 * baseline).round().max(1.0);
    // index of the first element >= cutoff, counted from the bottom
    let pos = sorted_ascending.partition_point(|&v| v < cutoff);
    n - pos
}

/// `estimate_recovered_cells_ordmag`: search a log2-spaced grid of candidate
/// recovered-cell counts for the one whose ordmag result reproduces it.
fn estimate_recovered_cells_ordmag(sample_sorted: &[f64], max_expected_cells: usize) -> (f64, f64) {
    let n = sample_sorted.len();
    let hi = (max_expected_cells as f64).log2();
    let mut candidates: Vec<usize> = (0..2000)
        .map(|i| {
            let t = 1.0 + (hi - 1.0) * (i as f64) / 1999.0;
            t.exp2().round() as usize
        })
        .collect();
    candidates.dedup();

    let mut best = (candidates[0] as f64, f64::INFINITY);
    for &rc in &candidates {
        let baseline_idx =
            ((rc as f64 * (1.0 - ORDMAG_RECOVERED_CELLS_QUANTILE)).round() as usize).min(n - 1);
        let filtered = find_within_ordmag(sample_sorted, baseline_idx) as f64;
        let loss = (filtered - rc as f64).powi(2) / rc as f64;
        if loss < best.1 {
            best = (rc as f64, loss);
        }
    }
    best
}

/// A bootstrap resample (with replacement) of `counts`, returned sorted
/// ascending, which is the only form the ordmag routines use.
fn bootstrap_sorted(counts: &[f64], rng: &mut Rng) -> Vec<f64> {
    let mut s: Vec<f64> = (0..counts.len())
        .map(|_| counts[rng.below(counts.len())])
        .collect();
    s.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    s
}

/// `filter_cellular_barcodes_ordmag`: the number of barcodes to call as cells.
fn ordmag_num_cells(
    counts: &[f64],
    recovered_cells: Option<usize>,
    max_expected_cells: usize,
    log: &slog::Logger,
) -> usize {
    let nonzero: Vec<f64> = counts.iter().cloned().filter(|&c| c > 0.0).collect();
    if nonzero.is_empty() {
        return 0;
    }
    let mut rng = Rng::new(0);

    let recovered = match recovered_cells {
        Some(rc) => rc.max(MIN_RECOVERED_CELLS_PER_GEM_GROUP),
        None => {
            let mut sum = 0.0;
            for _ in 0..ORDMAG_NUM_BOOTSTRAP_SAMPLES {
                let sample = bootstrap_sorted(&nonzero, &mut rng);
                sum += estimate_recovered_cells_ordmag(&sample, max_expected_cells).0;
            }
            let mean = sum / ORDMAG_NUM_BOOTSTRAP_SAMPLES as f64;
            (mean.round() as usize).max(MIN_RECOVERED_CELLS_PER_GEM_GROUP)
        }
    };
    info!(log, "ordmag: recovered_cells = {}", recovered);

    let baseline_idx = ((recovered as f64 * (1.0 - ORDMAG_RECOVERED_CELLS_QUANTILE)).round()
        as usize)
        .min(nonzero.len() - 1);

    let mut top_n_boot = Vec::with_capacity(ORDMAG_NUM_BOOTSTRAP_SAMPLES);
    for _ in 0..ORDMAG_NUM_BOOTSTRAP_SAMPLES {
        let sample = bootstrap_sorted(&nonzero, &mut rng);
        top_n_boot.push(find_within_ordmag(&sample, baseline_idx) as f64);
    }
    let mean = top_n_boot.iter().sum::<f64>() / top_n_boot.len() as f64;
    let mut nbcs = mean.round() as usize;
    if nbcs == 0 {
        return 0;
    }

    // `summarize_bootstrapped_top_n`: if the barcode at the cutoff has ties,
    // take them all, unless that would inflate the call by more than 20%.
    let mut sorted_desc = nonzero.clone();
    sorted_desc.sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());
    let cutoff = sorted_desc[nbcs - 1];
    let mut index = nbcs - 1;
    while index + 1 < sorted_desc.len() && sorted_desc[index] == cutoff {
        index += 1;
        if (index + 1 - nbcs) as f64 > 0.20 * nbcs as f64 {
            return nbcs;
        }
        nbcs = index + 1;
    }
    nbcs
}

// ---------------------------------------------------------------------------
// multinomial likelihoods
// ---------------------------------------------------------------------------

/// `ln(Gamma(x))` by the Lanczos approximation, used only for `ln(n!)`.
fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    let x = x - 1.0;
    let mut a = G[0];
    let t = x + 7.5;
    for (i, &g) in G.iter().enumerate().skip(1) {
        a += g / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// `ln(n!)`.
fn ln_factorial(n: f64) -> f64 {
    ln_gamma(n + 1.0)
}

/// Multinomial log-likelihood of one barcode's counts under `logp`.
///
/// `l = ln(N!) - sum_i ln(x_i!) + sum_i x_i * ln(p_i)`, with zero counts
/// dropped since they contribute nothing.
fn observed_loglik(row: &[(u32, f64)], logp: &[f64], total: f64) -> f64 {
    let mut l = ln_factorial(total);
    for &(feat, x) in row {
        l -= ln_factorial(x);
        l += x * logp[feat as usize];
    }
    l
}

/// One Monte Carlo replicate: draw `max_n` items from the ambient profile and
/// record the cumulative log-likelihood at each of the `distinct_n` sizes.
///
/// The cumulative form is Cell Ranger's: the increment contributed by the
/// `i`-th draw of feature `f`, when that draw is the `R`-th of `f`, is
/// `ln(i) - ln(R) + ln(p_f)`, which telescopes to the multinomial
/// log-likelihood above.
#[allow(clippy::too_many_arguments)]
fn simulate_one(
    p_cumulative: &[f64],
    logp: &[f64],
    max_n: usize,
    distinct_n: &[usize],
    seen: &mut [u32],
    touched: &mut Vec<u32>,
    rng: &mut Rng,
    out: &mut [f64],
) {
    touched.clear();
    let mut cum = 0.0f64;
    let mut next_idx = 0usize;
    for i in 1..=max_n {
        let u = rng.next_f64();
        let feat = p_cumulative.partition_point(|&c| c < u).min(logp.len() - 1);
        let seen_f = &mut seen[feat];
        if *seen_f == 0 {
            touched.push(feat as u32);
        }
        *seen_f += 1;
        cum += (i as f64).ln() - (*seen_f as f64).ln() + logp[feat];
        while next_idx < distinct_n.len() && distinct_n[next_idx] == i {
            out[next_idx] = cum;
            next_idx += 1;
        }
    }
    for &f in touched.iter() {
        seen[f as usize] = 0;
    }
}

/// Benjamini-Hochberg adjusted p-values, in the original order.
fn adjust_pvalues_bh(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by(|&a, &b| p[b].partial_cmp(&p[a]).unwrap());
    let mut adj = vec![0.0f64; n];
    let mut running = 1.0f64;
    for (k, &i) in order.iter().enumerate() {
        let rank = (n - k) as f64;
        running = running.min(p[i] * n as f64 / rank);
        adj[i] = running.min(1.0);
    }
    adj
}

// ---------------------------------------------------------------------------
// driver
// ---------------------------------------------------------------------------

/// Run ordmag then EmptyDrops over a `quant` output directory.
pub fn call_cells(opts: CallCellsOpts) -> anyhow::Result<usize> {
    let log = opts.log;
    let matrix = read_matrix(opts.input_dir)?;
    let umis_per_bc = matrix.counts_per_barcode();
    let num_bcs = matrix.barcodes.len();

    info!(
        log,
        "read a {} barcode x {} feature matrix",
        num_bcs.to_formatted_string(&Locale::en),
        matrix.num_features.to_formatted_string(&Locale::en)
    );

    // --- stage 1: ordmag -------------------------------------------------
    let max_expected_cells = (opts.n_partitions / 2).min(MAX_RECOVERED_CELLS_PER_GEM_GROUP);
    let n_ordmag = ordmag_num_cells(&umis_per_bc, opts.expect_cells, max_expected_cells, log);

    let mut by_count: Vec<usize> = (0..num_bcs).collect();
    by_count.sort_unstable_by(|&a, &b| {
        umis_per_bc[b]
            .partial_cmp(&umis_per_bc[a])
            .unwrap()
            .then(a.cmp(&b))
    });
    let mut is_cell = vec![false; num_bcs];
    for &i in by_count.iter().take(n_ordmag) {
        is_cell[i] = true;
    }
    info!(
        log,
        "ordmag called {} cells",
        n_ordmag.to_formatted_string(&Locale::en)
    );

    // --- stage 2: EmptyDrops ---------------------------------------------
    // Barcodes ranked [n_partitions/2, n_partitions) by descending count are
    // taken to be empty partitions; the non-zero ones among them form the
    // ambient pool.
    let (lo, hi) = (opts.n_partitions / 2, opts.n_partitions);
    let empty_bcs: Vec<usize> = by_count
        .iter()
        .skip(lo)
        .take(hi.saturating_sub(lo))
        .cloned()
        .collect();
    let ambient_bcs: Vec<usize> = empty_bcs
        .iter()
        .cloned()
        .filter(|&i| umis_per_bc[i] > 0.0)
        .collect();
    let max_background_umis = empty_bcs
        .iter()
        .map(|&i| umis_per_bc[i])
        .fold(0.0f64, f64::max);

    info!(
        log,
        "empty barcode range {}-{}: {} barcodes, {} with non-zero counts, max background UMIs {}",
        lo,
        hi,
        empty_bcs.len(),
        ambient_bcs.len(),
        max_background_umis
    );

    let mut num_nonambient = 0usize;
    let mut num_candidates = 0usize;
    let mut emptydrops_ran = false;

    if ambient_bcs.is_empty() {
        info!(
            log,
            "no ambient barcodes in the expected range, skipping EmptyDrops. If the input matrix \
             only holds barcodes that already passed a filter, re-run generate-permit-list with a \
             lower --min-reads so the empty partitions are present."
        );
    } else {
        // ambient profile over the features that are non-zero somewhere
        let mut feature_totals = vec![0.0f64; matrix.num_features];
        for &b in &ambient_bcs {
            for &(f, c) in &matrix.rows[b] {
                feature_totals[f as usize] += c;
            }
        }
        let mut any_nonzero = vec![false; matrix.num_features];
        for row in &matrix.rows {
            for &(f, _) in row {
                any_nonzero[f as usize] = true;
            }
        }
        let eval_features: Vec<usize> = (0..matrix.num_features)
            .filter(|&f| any_nonzero[f])
            .collect();

        let observed: Vec<f64> = eval_features
            .iter()
            .map(|&f| feature_totals[f])
            .filter(|&v| v > 0.0)
            .collect();

        match sgt_proportions(&observed) {
            Err(e) => {
                info!(log, "skipping EmptyDrops: {}", e);
            }
            Ok((pstar, p0)) => {
                emptydrops_ran = true;
                let n_zero = eval_features.len() - observed.len();
                let p0_each = if n_zero == 0 { 0.0 } else { p0 / n_zero as f64 };

                // profile over eval_features, in that order
                let mut profile = Vec::with_capacity(eval_features.len());
                let mut k = 0usize;
                for &f in &eval_features {
                    if feature_totals[f] > 0.0 {
                        profile.push(pstar[k]);
                        k += 1;
                    } else {
                        profile.push(p0_each);
                    }
                }
                // renormalize in the absence of a zero class, as Cell Ranger does
                if n_zero == 0 {
                    let s: f64 = profile.iter().sum();
                    for v in profile.iter_mut() {
                        *v /= s;
                    }
                }

                // dense feature -> position in `profile`
                let mut feat_pos = vec![usize::MAX; matrix.num_features];
                for (i, &f) in eval_features.iter().enumerate() {
                    feat_pos[f] = i;
                }

                let logp: Vec<f64> = profile
                    .iter()
                    .map(|&v| v.max(f64::MIN_POSITIVE).ln())
                    .collect();
                let mut p_cumulative = Vec::with_capacity(profile.len());
                let mut acc = 0.0;
                for &v in &profile {
                    acc += v;
                    p_cumulative.push(acc);
                }

                // candidates: not already cells, not in the empty range, and
                // above both the minimum and the largest background barcode
                let mut in_empty = vec![false; num_bcs];
                for &i in &empty_bcs {
                    in_empty[i] = true;
                }
                let min_umis = (opts.min_umis as f64).max(1.0 + max_background_umis);
                let mut eval_bcs: Vec<usize> = (0..num_bcs)
                    .filter(|&i| !is_cell[i] && !in_empty[i] && umis_per_bc[i] >= min_umis)
                    .collect();
                eval_bcs.sort_unstable();
                num_candidates = eval_bcs.len();

                info!(
                    log,
                    "EmptyDrops: {} candidate barcodes above {} UMIs",
                    num_candidates.to_formatted_string(&Locale::en),
                    min_umis
                );

                if !eval_bcs.is_empty() {
                    // observed log-likelihoods
                    let mut obs_loglk = Vec::with_capacity(eval_bcs.len());
                    for &b in &eval_bcs {
                        let row: Vec<(u32, f64)> = matrix.rows[b]
                            .iter()
                            .filter(|&&(f, _)| feat_pos[f as usize] != usize::MAX)
                            .map(|&(f, c)| (feat_pos[f as usize] as u32, c))
                            .collect();
                        let total: f64 = row.iter().map(|&(_, c)| c).sum();
                        obs_loglk.push(observed_loglik(&row, &logp, total));
                    }

                    // simulated log-likelihoods, one column per replicate
                    let mut distinct_n: Vec<usize> = eval_bcs
                        .iter()
                        .map(|&b| umis_per_bc[b].round() as usize)
                        .collect();
                    distinct_n.sort_unstable();
                    distinct_n.dedup();
                    let max_n = *distinct_n.last().unwrap();

                    info!(
                        log,
                        "EmptyDrops: {} simulations of up to {} draws over {} distinct sizes",
                        opts.num_sims.to_formatted_string(&Locale::en),
                        max_n.to_formatted_string(&Locale::en),
                        distinct_n.len().to_formatted_string(&Locale::en)
                    );

                    let num_threads = opts.threads.max(1);
                    let chunk = opts.num_sims.div_ceil(num_threads);
                    let sim_loglk: Vec<Vec<f64>> = std::thread::scope(|scope| {
                        let mut handles = Vec::new();
                        for t in 0..num_threads {
                            let start = t * chunk;
                            let end = ((t + 1) * chunk).min(opts.num_sims);
                            if start >= end {
                                continue;
                            }
                            let (p_cumulative, logp, distinct_n) =
                                (&p_cumulative, &logp, &distinct_n);
                            handles.push(scope.spawn(move || {
                                // seeded per thread so the run is reproducible
                                // regardless of how the work is split
                                let mut rng = Rng::new(0x5EED_0000 + t as u64);
                                let mut seen = vec![0u32; logp.len()];
                                let mut touched: Vec<u32> = Vec::new();
                                let mut out = vec![0.0f64; distinct_n.len()];
                                let mut cols = Vec::with_capacity(end - start);
                                for _ in start..end {
                                    simulate_one(
                                        p_cumulative,
                                        logp,
                                        max_n,
                                        distinct_n,
                                        &mut seen,
                                        &mut touched,
                                        &mut rng,
                                        &mut out,
                                    );
                                    cols.push(out.clone());
                                }
                                cols
                            }));
                        }
                        handles
                            .into_iter()
                            .flat_map(|h| h.join().expect("simulation thread panicked"))
                            .collect()
                    });

                    // p-values: fraction of simulated log-likelihoods, at the
                    // same N, that are below the observed one
                    let num_sims = sim_loglk.len();
                    let mut pvalues = Vec::with_capacity(eval_bcs.len());
                    for (i, &b) in eval_bcs.iter().enumerate() {
                        let n = umis_per_bc[b].round() as usize;
                        let idx = distinct_n.partition_point(|&d| d < n);
                        let mut lower = 0usize;
                        for col in &sim_loglk {
                            if col[idx] < obs_loglk[i] {
                                lower += 1;
                            }
                        }
                        pvalues.push((1 + lower) as f64 / (1 + num_sims) as f64);
                    }

                    let padj = adjust_pvalues_bh(&pvalues);
                    for (i, &b) in eval_bcs.iter().enumerate() {
                        if padj[i] <= opts.fdr {
                            is_cell[b] = true;
                            num_nonambient += 1;
                        }
                    }
                    info!(
                        log,
                        "EmptyDrops rescued {} barcodes at FDR {}",
                        num_nonambient.to_formatted_string(&Locale::en),
                        opts.fdr
                    );
                }
            }
        }
    }

    // --- output ----------------------------------------------------------
    let out_dir: &PathBuf = opts.output_dir;
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("could not create {}", out_dir.display()))?;

    let bc_path = out_dir.join("filtered_barcodes.txt");
    {
        let mut w = BufWriter::new(
            std::fs::File::create(&bc_path)
                .with_context(|| format!("could not create {}", bc_path.display()))?,
        );
        for &i in by_count.iter() {
            if is_cell[i] {
                writeln!(w, "{}", matrix.barcodes[i])?;
            }
        }
    }

    let num_cells = is_cell.iter().filter(|&&b| b).count();
    let meta = json!({
        "cmd": opts.cmdline,
        "version_str": opts.version,
        "num_barcodes": num_bcs,
        "num_features": matrix.num_features,
        "ordmag_cells": n_ordmag,
        "emptydrops_ran": emptydrops_ran,
        "emptydrops_candidates": num_candidates,
        "emptydrops_rescued": num_nonambient,
        "num_cells": num_cells,
        "fdr": opts.fdr,
        "num_sims": opts.num_sims,
        "n_partitions": opts.n_partitions,
        "min_umis": opts.min_umis,
    });
    let json_path = out_dir.join("call_cells.json");
    std::fs::write(&json_path, serde_json::to_string_pretty(&meta)?)
        .with_context(|| format!("could not write {}", json_path.display()))?;

    info!(
        log,
        "called {} cells in total, written to {}",
        num_cells.to_formatted_string(&Locale::en),
        bc_path.display()
    );
    Ok(num_cells)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ln_factorial_matches_known_values() {
        assert!((ln_factorial(0.0) - 0.0).abs() < 1e-9);
        assert!((ln_factorial(1.0) - 0.0).abs() < 1e-9);
        // ln(10!) = ln(3628800)
        assert!((ln_factorial(10.0) - 3_628_800f64.ln()).abs() < 1e-9);
        // ln(170!) stays finite where 170! itself would not overflow f64
        assert!(ln_factorial(170.0).is_finite());
    }

    #[test]
    fn bh_adjustment_is_monotone_and_bounded() {
        let p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205];
        let adj = adjust_pvalues_bh(&p);
        for (a, b) in adj.iter().zip(p.iter()) {
            assert!(*a >= *b - 1e-12, "adjusted p must not fall below raw p");
            assert!(*a <= 1.0);
        }
        // sorted order is preserved, which is what BH guarantees
        let mut sorted = adj.clone();
        sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(adj, sorted);
    }

    #[test]
    fn find_within_ordmag_takes_everything_above_a_tenth() {
        // ascending counts; baseline is the top element, cutoff = 10
        let counts = [1.0, 5.0, 9.0, 10.0, 50.0, 100.0];
        assert_eq!(find_within_ordmag(&counts, 0), 3);
    }

    #[test]
    fn sgt_rejects_too_few_distinct_frequencies() {
        let freqs = [1.0, 1.0, 2.0, 2.0, 3.0];
        assert!(sgt_proportions(&freqs).is_err());
    }

    #[test]
    fn sgt_proportions_sum_with_p0_to_one() {
        // 12 distinct frequencies, enough for the estimator to apply
        let mut freqs = Vec::new();
        for r in 1..=12u32 {
            for _ in 0..(20 - r) {
                freqs.push(r as f64);
            }
        }
        let (pstar, p0) = sgt_proportions(&freqs).expect("SGT should apply here");
        let total: f64 = pstar.iter().sum::<f64>() + p0;
        assert!((total - 1.0).abs() < 1e-6, "total was {total}");
    }
}
