/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Groups of cells that share isoform usage, found from the molecules with one candidate
//! transcript (`--truncation-discover-groups`), so that abundances can be pooled within groups
//! (`--truncation-groups`) without a file of labels.
//!
//! Per gene with several isoforms, a cell's one-candidate counts over the gene's transcripts are
//! multinomial around the shares of its group. The cells are split by bisecting. Each leaf's
//! molecules are thinned at random into four parts: the first two give the residual features
//! `(x - n phi) / sqrt(n phi (1 - phi))` on the genes whose usage varies between the leaf's cells
//! beyond sampling, and the candidate splits (k-means and Ward clusters of the leading singular
//! vectors of the residuals for several numbers of clusters, and the tails of each component; each
//! cluster against the rest); the third ranks the candidates by the held-out likelihood of its
//! molecules; the fourth confirms the best one against random sides of the same size. A side's
//! shares are its counts shrunk toward the leaf's with a strength per gene estimated from how much
//! the two sides differ beyond sampling (as `--truncation-group-shrink auto-gene` does), so that a
//! small side is penalized only on the genes where it does differ. A split is kept when its
//! confirmation gain beats the best random side by `min_gain` nats in total and
//! `min_gain_per_cell` per cell. Finally each cell moves to the group whose shares fit all its
//! molecules best, until none moves. Without any accepted split every cell is in one group: the
//! pooled model. `scripts/truncation_groups.py` runs the same procedure on a `--truncation-dump`.

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use nalgebra::{DMatrix, SymmetricEigen};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Binomial, Distribution, StandardNormal};

/// The discovery's settings; the defaults are those of `scripts/truncation_groups.py`.
#[derive(Clone, Debug)]
pub struct DiscoverySettings {
    /// A leaf below this many cells is not split, and no side may be smaller.
    pub min_cells: usize,
    /// Nats by which a split's confirmation gain must beat the best random side of the same size.
    pub min_gain: f64,
    /// Random sides of the same size the confirmation gain is compared with.
    pub null_sides: usize,
    /// The gain above the best random side must also exceed this per cell of the leaf.
    pub min_gain_per_cell: f64,
    /// Leading singular vectors of the residuals the candidates are found in.
    pub components: usize,
    /// A gene takes part when at least this many one-candidate molecules fall in it.
    pub min_gene_molecules: f64,
    /// Null standard deviations a gene's between-side deviation must exceed to get a strength of
    /// its own in the shrunk shares.
    pub gate_z: f64,
    /// Pairs of pseudo-molecules pulling a gene's strength toward the common one.
    pub gene_prior_weight: f64,
    /// Depth of the bisection.
    pub max_depth: usize,
    pub seed: u64,
    pub threads: usize,
}

impl Default for DiscoverySettings {
    fn default() -> Self {
        DiscoverySettings {
            min_cells: 8,
            min_gain: 20.0,
            null_sides: 5,
            min_gain_per_cell: 0.02,
            components: 20,
            min_gene_molecules: 50.0,
            gate_z: 4.0,
            gene_prior_weight: 250.0,
            max_depth: 6,
            seed: 1,
            threads: 1,
        }
    }
}

/// One attempted split of a leaf, as the discovery logs it.
#[derive(Clone, Debug)]
pub struct SplitRecord {
    pub depth: usize,
    pub cells: usize,
    pub variable_genes: usize,
    pub side: usize,
    pub rest: usize,
    pub how: String,
    pub selection_gain: f64,
    pub gain: f64,
    pub null_max: f64,
    pub per_cell: f64,
    pub genes_finite: usize,
    pub accepted: bool,
}

impl fmt::Display for SplitRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "depth {}: {} cells, {} variable genes -> {} + {} ({}; selection gain {:.1}), confirmation gain {:.1} nats against random sides {:.1} at most ({:.3}/cell above), {} genes with a finite strength: {}",
            self.depth,
            self.cells,
            self.variable_genes,
            self.side,
            self.rest,
            self.how,
            self.selection_gain,
            self.gain,
            self.null_max,
            self.per_cell,
            self.genes_finite,
            if self.accepted { "ACCEPT" } else { "reject" }
        )
    }
}

/// The groups found ([`discover_groups`]).
#[derive(Clone, Debug)]
pub struct Discovery {
    /// Each cell's group, `0..n_groups`, in the order of the input.
    pub labels: Vec<u32>,
    pub n_groups: usize,
    /// Cells per group.
    pub sizes: Vec<usize>,
    /// Genes with at least two isoforms seen and enough molecules, and their transcripts.
    pub genes: usize,
    pub transcripts: usize,
    /// One-candidate molecules in those genes.
    pub molecules: f64,
    pub splits: Vec<SplitRecord>,
    /// Cells moved in each refinement pass.
    pub refinement: Vec<usize>,
}

/// Runs `f(i)` for `i in 0..n` on `threads` threads; the results are in index order whatever the
/// scheduling, so the discovery does not depend on the number of threads.
fn par_map<T: Send, F: Fn(usize) -> T + Sync>(n: usize, threads: usize, f: F) -> Vec<T> {
    if threads <= 1 || n <= 1 {
        return (0..n).map(&f).collect();
    }
    let next = AtomicUsize::new(0);
    let slots: Vec<Mutex<Option<T>>> = (0..n).map(|_| Mutex::new(None)).collect();
    std::thread::scope(|s| {
        for _ in 0..threads.min(n) {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= n {
                    break;
                }
                let r = f(i);
                *slots[i].lock().unwrap() = Some(r);
            });
        }
    });
    slots
        .into_iter()
        .map(|m| m.into_inner().unwrap().expect("every task ran"))
        .collect()
}

/// Ranges that split `0..n` into about `parts` pieces.
fn chunks(n: usize, parts: usize) -> Vec<(usize, usize)> {
    let parts = parts.clamp(1, n.max(1));
    let size = n.div_ceil(parts);
    (0..n).step_by(size.max(1)).map(|a| (a, (a + size).min(n))).collect()
}

/// A sparse matrix in compressed-row form, with the columns of each row in increasing order.
#[derive(Clone)]
struct Csr {
    rows: usize,
    cols: usize,
    indptr: Vec<usize>,
    indices: Vec<u32>,
    data: Vec<f64>,
}

impl Csr {
    fn row(&self, r: usize) -> (&[u32], &[f64]) {
        let (a, b) = (self.indptr[r], self.indptr[r + 1]);
        (&self.indices[a..b], &self.data[a..b])
    }

    fn nnz(&self) -> usize {
        self.indices.len()
    }

    /// Column sums over the rows where `mask` is `value` (all rows without a mask).
    fn colsum(&self, mask: Option<(&[bool], bool)>) -> Vec<f64> {
        let mut out = vec![0.0; self.cols];
        for r in 0..self.rows {
            if let Some((m, v)) = mask {
                if m[r] != v {
                    continue;
                }
            }
            let (idx, dat) = self.row(r);
            for (&c, &x) in idx.iter().zip(dat) {
                out[c as usize] += x;
            }
        }
        out
    }

    fn select_rows(&self, rows: &[usize]) -> Csr {
        let mut indptr = Vec::with_capacity(rows.len() + 1);
        let mut indices = Vec::new();
        let mut data = Vec::new();
        indptr.push(0);
        for &r in rows {
            let (idx, dat) = self.row(r);
            indices.extend_from_slice(idx);
            data.extend_from_slice(dat);
            indptr.push(indices.len());
        }
        Csr {
            rows: rows.len(),
            cols: self.cols,
            indptr,
            indices,
            data,
        }
    }

    /// Elementwise sum of two matrices with the same pattern.
    fn add(&self, o: &Csr) -> Csr {
        debug_assert_eq!(self.indices, o.indices);
        Csr {
            rows: self.rows,
            cols: self.cols,
            indptr: self.indptr.clone(),
            indices: self.indices.clone(),
            data: self.data.iter().zip(&o.data).map(|(a, b)| a + b).collect(),
        }
    }

    fn transpose(&self) -> Csr {
        let mut count = vec![0_usize; self.cols + 1];
        for &c in &self.indices {
            count[c as usize + 1] += 1;
        }
        for c in 0..self.cols {
            count[c + 1] += count[c];
        }
        let indptr = count.clone();
        let mut pos = count;
        let mut indices = vec![0_u32; self.nnz()];
        let mut data = vec![0.0; self.nnz()];
        for r in 0..self.rows {
            let (idx, dat) = self.row(r);
            for (&c, &x) in idx.iter().zip(dat) {
                let p = pos[c as usize];
                indices[p] = r as u32;
                data[p] = x;
                pos[c as usize] += 1;
            }
        }
        Csr {
            rows: self.cols,
            cols: self.rows,
            indptr,
            indices,
            data,
        }
    }

    /// `self * dense`, where `dense` is `cols x l` in row-major order; the result is `rows x l`.
    fn mul_dense(&self, dense: &[f64], l: usize, threads: usize) -> Vec<f64> {
        let parts = chunks(self.rows, threads * 4);
        let out = par_map(parts.len(), threads, |p| {
            let (a, b) = parts[p];
            let mut out = vec![0.0; (b - a) * l];
            for i in a..b {
                let (idx, dat) = self.row(i);
                let o = &mut out[(i - a) * l..(i - a + 1) * l];
                for (&c, &v) in idx.iter().zip(dat) {
                    let w = &dense[c as usize * l..(c as usize + 1) * l];
                    for j in 0..l {
                        o[j] += v * w[j];
                    }
                }
            }
            out
        });
        out.concat()
    }
}

/// Multinomial thinning of integer counts into `parts` equal parts, keeping the pattern.
fn thin(x: &Csr, parts: usize, rng: &mut StdRng) -> Vec<Csr> {
    let mut data: Vec<Vec<f64>> = (0..parts).map(|_| vec![0.0; x.nnz()]).collect();
    for (i, &v) in x.data.iter().enumerate() {
        let mut remaining = v.round().max(0.0) as u64;
        for p in 0..parts {
            let draw = if p + 1 == parts || remaining == 0 {
                remaining
            } else {
                Binomial::new(remaining, 1.0 / (parts - p) as f64)
                    .expect("valid binomial")
                    .sample(rng)
            };
            data[p][i] = draw as f64;
            remaining -= draw;
        }
    }
    data.into_iter()
        .map(|d| Csr {
            rows: x.rows,
            cols: x.cols,
            indptr: x.indptr.clone(),
            indices: x.indices.clone(),
            data: d,
        })
        .collect()
}

/// The columns' genes: the columns of gene `g` are `start[g]..start[g + 1]`, in order.
struct Genes {
    gi: Vec<u32>,
    start: Vec<usize>,
}

impl Genes {
    fn n(&self) -> usize {
        self.start.len() - 1
    }

    fn cols(&self, g: usize) -> std::ops::Range<usize> {
        self.start[g]..self.start[g + 1]
    }

    /// Calls `f(gene, n, cols, values)` for each gene a row has entries in, with `n` the sum of
    /// the row's entries in the gene (the row's columns are sorted, and a gene's are contiguous).
    fn genes_of_row(&self, idx: &[u32], dat: &[f64], mut f: impl FnMut(usize, f64, &[u32], &[f64])) {
        let mut i = 0;
        while i < idx.len() {
            let g = self.gi[idx[i] as usize] as usize;
            let mut j = i;
            let mut n = 0.0;
            while j < idx.len() && self.gi[idx[j] as usize] as usize == g {
                n += dat[j];
                j += 1;
            }
            f(g, n, &idx[i..j], &dat[i..j]);
            i = j;
        }
    }

    /// Each transcript's share of its gene from summed counts `x` (uniform where the gene has none).
    fn plain_shares(&self, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; x.len()];
        for g in 0..self.n() {
            let r = self.cols(g);
            let n: f64 = x[r.clone()].iter().sum();
            if n > 0.0 {
                for c in r {
                    out[c] = x[c] / n;
                }
            } else {
                let u = 1.0 / (r.end - r.start) as f64;
                for c in r {
                    out[c] = u;
                }
            }
        }
        out
    }

    /// Shares of each side shrunk toward `phi` with a strength per gene from how much the sides
    /// differ beyond sampling: method of moments over the two sides with the pooled share plugged
    /// in (so that a unit's deviation counts with weight `1 - 2 n/N + S2/N^2`), gated at `gate_z`
    /// null standard deviations and shrunk toward the common strength with `w` pairs; returns the
    /// two sides' shares and the strength per gene (infinite where the gene is pooled).
    fn eb_shares(
        &self,
        x_side: &[f64],
        x_rest: &[f64],
        phi: &[f64],
        gate_z: f64,
        w: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let n_genes = self.n();
        let mut n_s = vec![0.0; n_genes];
        let mut n_r = vec![0.0; n_genes];
        for g in 0..n_genes {
            let r = self.cols(g);
            n_s[g] = x_side[r.clone()].iter().sum();
            n_r[g] = x_rest[r].iter().sum();
        }
        let mut num = vec![0.0; n_genes];
        let mut den = vec![0.0; n_genes];
        let mut null_var = vec![0.0; n_genes];
        for t in 0..phi.len() {
            let g = self.gi[t] as usize;
            let (ns, nr, p) = (n_s[g], n_r[g], phi[t]);
            let v = p * (1.0 - p);
            let total = ns + nr;
            let s2_all = ns * ns + nr * nr;
            let (mut dev, mut sn, mut sn2, mut sn3) = (0.0, 0.0, 0.0, 0.0);
            for (x, n) in [(x_side[t], ns), (x_rest[t], nr)] {
                if n >= 2.0 {
                    dev += (x - n * p) * (x - n * p);
                    sn += n;
                    sn2 += n * n;
                    sn3 += n * n * n;
                }
            }
            let npos = total.max(1e-300);
            let e = sn - sn2 / npos;
            let d = if total > 0.0 {
                (sn2 - sn) - 2.0 * (sn3 - sn2) / npos + sn2 * (s2_all - total) / (npos * npos)
            } else {
                0.0
            };
            let exc = dev - v * e;
            if d > 0.0 && total > 0.0 {
                num[g] += d * v;
                den[g] += exc;
                null_var[g] += 2.0 * d * v * v;
            }
        }
        let use_g: Vec<bool> = (0..n_genes)
            .map(|g| num[g] > 0.0 && den[g] > gate_z * null_var[g].sqrt())
            .collect();
        let (mut sd, mut sn) = (0.0, 0.0);
        for g in 0..n_genes {
            if use_g[g] {
                sd += den[g];
                sn += num[g];
            }
        }
        let common = if use_g.iter().any(|&u| u) {
            sd / sn.max(1e-300)
        } else {
            0.0
        };
        let beta: Vec<f64> = (0..n_genes)
            .map(|g| {
                if !use_g[g] {
                    return f64::INFINITY;
                }
                let r = (den[g] + w * common) / (num[g] + w);
                if r > 0.0 {
                    (1.0 / r.max(1e-300) - 1.0).max(1.0)
                } else {
                    f64::INFINITY
                }
            })
            .collect();
        let mut s_side = vec![0.0; phi.len()];
        let mut s_rest = vec![0.0; phi.len()];
        for t in 0..phi.len() {
            let g = self.gi[t] as usize;
            let b = beta[g].min(1e9);
            s_side[t] = (x_side[t] + b * phi[t]) / (n_s[g] + b).max(1e-300);
            s_rest[t] = (x_rest[t] + b * phi[t]) / (n_r[g] + b).max(1e-300);
        }
        (s_side, s_rest, beta)
    }

    /// Genes whose isoform shares vary across the rows of `x` beyond sampling: the excess of
    /// `sum_c (x - n phi)^2` over its multinomial expectation as a z-score against the spread of
    /// the per-cell terms; returns a mask over the columns (true: keep) and the genes kept.
    fn variable_genes(&self, x: &Csr, phi: &[f64]) -> (Vec<bool>, usize) {
        let mut exc = vec![0.0; x.cols];
        let mut sd2 = vec![0.0; x.cols];
        for r in 0..x.rows {
            let (idx, dat) = x.row(r);
            self.genes_of_row(idx, dat, |g, n, cols, vals| {
                if n <= 0.0 {
                    return;
                }
                let mut k = 0;
                for c in self.cols(g) {
                    let xv = if k < cols.len() && cols[k] as usize == c {
                        k += 1;
                        vals[k - 1]
                    } else {
                        0.0
                    };
                    let e = n * phi[c];
                    let v = e * (1.0 - phi[c]);
                    let term = (xv - e) * (xv - e) - v;
                    exc[c] += term;
                    sd2[c] += term * term;
                }
            });
        }
        let n_genes = self.n();
        let mut keep_g = vec![false; n_genes];
        let mut kept = 0;
        for g in 0..n_genes {
            let r = self.cols(g);
            let e: f64 = exc[r.clone()].iter().sum();
            let s: f64 = r.map(|c| sd2[c] + 1e-12).sum::<f64>().sqrt();
            keep_g[g] = e / s.max(1e-12) > 2.0;
            kept += keep_g[g] as usize;
        }
        if kept < 10 {
            keep_g.iter_mut().for_each(|k| *k = true);
            kept = n_genes;
        }
        ((0..x.cols).map(|c| keep_g[self.gi[c] as usize]).collect(), kept)
    }

    /// Pearson residuals of the counts around shares `phi`, over every transcript of each gene a
    /// cell has molecules in, clipped to ±10 and restricted to the columns in `keep`.
    fn residuals(&self, x: &Csr, phi: &[f64], keep: &[bool]) -> Csr {
        let mut new_col = vec![u32::MAX; x.cols];
        let mut n_keep = 0;
        for c in 0..x.cols {
            if keep[c] {
                new_col[c] = n_keep as u32;
                n_keep += 1;
            }
        }
        let mut indptr = Vec::with_capacity(x.rows + 1);
        let mut indices = Vec::new();
        let mut data = Vec::new();
        indptr.push(0);
        for r in 0..x.rows {
            let (idx, dat) = x.row(r);
            self.genes_of_row(idx, dat, |g, n, cols, vals| {
                if n <= 0.0 {
                    return;
                }
                let mut k = 0;
                for c in self.cols(g) {
                    let xv = if k < cols.len() && cols[k] as usize == c {
                        k += 1;
                        vals[k - 1]
                    } else {
                        0.0
                    };
                    if new_col[c] == u32::MAX {
                        continue;
                    }
                    let e = n * phi[c];
                    let v = e * (1.0 - phi[c]);
                    let res = ((xv - e) / (v + 0.25).sqrt()).clamp(-10.0, 10.0);
                    indices.push(new_col[c]);
                    data.push(res);
                }
            });
            indptr.push(indices.len());
        }
        Csr {
            rows: x.rows,
            cols: n_keep,
            indptr,
            indices,
            data,
        }
    }
}

/// Log-likelihood of each row's counts under shares `phi`.
fn loglik_rows(x: &Csr, phi: &[f64]) -> Vec<f64> {
    let lp: Vec<f64> = phi.iter().map(|p| p.max(1e-300).ln()).collect();
    (0..x.rows)
        .map(|r| {
            let (idx, dat) = x.row(r);
            idx.iter().zip(dat).map(|(&c, &v)| v * lp[c as usize]).sum()
        })
        .collect()
}

/// The thin Q of a QR decomposition of the row-major `n x l` matrix `y`, row-major.
fn qr_q(y: &[f64], n: usize, l: usize) -> Vec<f64> {
    let q = DMatrix::from_row_slice(n, l, y).qr().q();
    let mut out = vec![0.0; n * l];
    for i in 0..n {
        for j in 0..l.min(q.ncols()) {
            out[i * l + j] = q[(i, j)];
        }
    }
    out
}

/// The leading `k` left singular vectors of `r`, scaled by their singular values (randomized
/// subspace iteration), row-major `rows x k`.
fn embed(r: &Csr, k: usize, seed: u64, threads: usize) -> Vec<f64> {
    let (n, m) = (r.rows, r.cols);
    let l = (k + 10).min(n.min(m));
    let mut rng = StdRng::seed_from_u64(seed);
    let omega: Vec<f64> = (0..m * l).map(|_| rng.sample(StandardNormal)).collect();
    let rt = r.transpose();
    let mut q = qr_q(&r.mul_dense(&omega, l, threads), n, l);
    for _ in 0..7 {
        let q2 = qr_q(&rt.mul_dense(&q, l, threads), m, l);
        q = qr_q(&r.mul_dense(&q2, l, threads), n, l);
    }
    // B^T = R^T Q (m x l); the eigenvectors of B B^T rotate Q onto the singular vectors
    let bt = rt.mul_dense(&q, l, threads);
    let mut c: DMatrix<f64> = DMatrix::zeros(l, l);
    for row in bt.chunks(l) {
        for j in 0..l {
            for h in 0..l {
                c[(j, h)] += row[j] * row[h];
            }
        }
    }
    let eig = SymmetricEigen::new(c);
    let mut order: Vec<usize> = (0..l).collect();
    order.sort_by(|&a, &b| eig.eigenvalues[b].partial_cmp(&eig.eigenvalues[a]).unwrap());
    let mut z = vec![0.0; n * k];
    for (j, &e) in order.iter().take(k).enumerate() {
        let s = eig.eigenvalues[e].max(0.0).sqrt();
        for i in 0..n {
            let mut u = 0.0;
            for t in 0..l {
                u += q[i * l + t] * eig.eigenvectors[(t, e)];
            }
            z[i * k + j] = u * s;
        }
    }
    z
}

fn sqdist(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum()
}

/// One k-means run: greedy k-means++ seeding, then Lloyd's iterations to sklearn's tolerance.
fn kmeans_once(z: &[f64], n: usize, d: usize, k: usize, rng: &mut StdRng) -> (Vec<u32>, f64) {
    let pt = |i: usize| &z[i * d..(i + 1) * d];
    let mut centers = vec![0.0; k * d];
    let first = rng.random_range(0..n);
    centers[..d].copy_from_slice(pt(first));
    let mut closest: Vec<f64> = (0..n).map(|i| sqdist(pt(i), pt(first))).collect();
    let trials = 2 + (k as f64).ln() as usize;
    for c in 1..k {
        let pot: f64 = closest.iter().sum();
        let mut best: Option<(usize, f64, Vec<f64>)> = None;
        for _ in 0..trials {
            let cand = if pot > 0.0 {
                let mut u = rng.random::<f64>() * pot;
                let mut cand = n - 1;
                for (i, &w) in closest.iter().enumerate() {
                    if u < w {
                        cand = i;
                        break;
                    }
                    u -= w;
                }
                cand
            } else {
                rng.random_range(0..n)
            };
            let newd: Vec<f64> = (0..n)
                .map(|i| closest[i].min(sqdist(pt(i), pt(cand))))
                .collect();
            let p: f64 = newd.iter().sum();
            if best.as_ref().is_none_or(|b| p < b.1) {
                best = Some((cand, p, newd));
            }
        }
        let (cand, _, newd) = best.unwrap();
        centers[c * d..(c + 1) * d].copy_from_slice(pt(cand));
        closest = newd;
    }
    // sklearn's tolerance: 1e-4 of the mean per-dimension variance
    let mut tol = 0.0;
    for t in 0..d {
        let mean = (0..n).map(|i| z[i * d + t]).sum::<f64>() / n as f64;
        tol += (0..n).map(|i| (z[i * d + t] - mean).powi(2)).sum::<f64>() / n as f64;
    }
    let tol = 1e-4 * tol / d as f64;
    let mut labels = vec![u32::MAX; n];
    for _ in 0..300 {
        let mut changed = false;
        for i in 0..n {
            let (mut bj, mut bd) = (0, f64::INFINITY);
            for j in 0..k {
                let dd = sqdist(pt(i), &centers[j * d..(j + 1) * d]);
                if dd < bd {
                    bd = dd;
                    bj = j;
                }
            }
            if labels[i] != bj as u32 {
                labels[i] = bj as u32;
                changed = true;
            }
        }
        if !changed {
            break;
        }
        let mut sums = vec![0.0; k * d];
        let mut cnt = vec![0_usize; k];
        for i in 0..n {
            let j = labels[i] as usize;
            cnt[j] += 1;
            for t in 0..d {
                sums[j * d + t] += z[i * d + t];
            }
        }
        let mut shift = 0.0;
        for j in 0..k {
            if cnt[j] == 0 {
                // an empty cluster moves to the point farthest from its center
                let far = (0..n)
                    .max_by(|&a, &b| {
                        let da = sqdist(pt(a), &centers[labels[a] as usize * d..(labels[a] as usize + 1) * d]);
                        let db = sqdist(pt(b), &centers[labels[b] as usize * d..(labels[b] as usize + 1) * d]);
                        da.partial_cmp(&db).unwrap()
                    })
                    .unwrap();
                sums[j * d..(j + 1) * d].copy_from_slice(pt(far));
                cnt[j] = 1;
                labels[far] = j as u32;
            }
            for t in 0..d {
                let nc = sums[j * d + t] / cnt[j] as f64;
                shift += (nc - centers[j * d + t]).powi(2);
                centers[j * d + t] = nc;
            }
        }
        if shift <= tol {
            break;
        }
    }
    let inertia = (0..n)
        .map(|i| sqdist(pt(i), &centers[labels[i] as usize * d..(labels[i] as usize + 1) * d]))
        .sum();
    (labels, inertia)
}

/// k-means with `n_init` seedings in parallel; the labels of the run with the lowest inertia.
fn kmeans(z: &[f64], n: usize, d: usize, k: usize, n_init: usize, seed: u64, threads: usize) -> Vec<u32> {
    let runs = par_map(n_init, threads, |i| {
        let mut rng = StdRng::seed_from_u64(seed.wrapping_add(0x9e37_79b9_7f4a_7c15u64.wrapping_mul(i as u64 + 1)));
        kmeans_once(z, n, d, k, &mut rng)
    });
    runs.into_iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .map(|(labels, _)| labels)
        .unwrap()
}

/// Ward's agglomerative clustering by the nearest-neighbor chain: the `n - 1` merges, each as
/// the two cluster ids merged (points are `0..n`, the i-th merge creates id `n + i`) and the
/// Ward distance `|A||B| / (|A| + |B|) * |c_A - c_B|^2`.
fn ward_merges(z: &[f64], n: usize, d: usize) -> Vec<(usize, usize, f64)> {
    let total = 2 * n - 1;
    let mut cent = z.to_vec();
    cent.resize(total * d, 0.0);
    let mut size = vec![1.0_f64; total];
    let mut active: Vec<usize> = (0..n).collect();
    let mut pos_in_active = vec![usize::MAX; total];
    for (p, &a) in active.iter().enumerate() {
        pos_in_active[a] = p;
    }
    let mut merges = Vec::with_capacity(n - 1);
    let mut chain: Vec<usize> = Vec::new();
    let dist = |a: usize, b: usize, cent: &[f64], size: &[f64]| {
        size[a] * size[b] / (size[a] + size[b]) * sqdist(&cent[a * d..(a + 1) * d], &cent[b * d..(b + 1) * d])
    };
    while active.len() > 1 {
        if chain.is_empty() {
            chain.push(active[0]);
        }
        loop {
            let a = *chain.last().unwrap();
            let prev = if chain.len() >= 2 { Some(chain[chain.len() - 2]) } else { None };
            let (mut best, mut bd) = match prev {
                Some(p) => (p, dist(a, p, &cent, &size)),
                None => (usize::MAX, f64::INFINITY),
            };
            for &b in &active {
                if b != a {
                    let dd = dist(a, b, &cent, &size);
                    if dd < bd {
                        bd = dd;
                        best = b;
                    }
                }
            }
            if Some(best) == prev {
                chain.pop();
                chain.pop();
                let new = n + merges.len();
                let (sa, sb) = (size[a], size[best]);
                for t in 0..d {
                    cent[new * d + t] = (sa * cent[a * d + t] + sb * cent[best * d + t]) / (sa + sb);
                }
                size[new] = sa + sb;
                for old in [a, best] {
                    let p = pos_in_active[old];
                    let last = *active.last().unwrap();
                    active.swap_remove(p);
                    if p < active.len() {
                        pos_in_active[last] = p;
                    }
                    pos_in_active[old] = usize::MAX;
                }
                pos_in_active[new] = active.len();
                active.push(new);
                merges.push((a, best, bd));
                break;
            }
            chain.push(best);
        }
    }
    merges
}

/// Cuts the dendrogram of `merges` into `k` clusters for each `k` in `ks`: labels per point.
fn ward_cuts(merges: &[(usize, usize, f64)], n: usize, ks: &[usize]) -> Vec<(usize, Vec<u32>)> {
    let mut order: Vec<usize> = (0..merges.len()).collect();
    order.sort_by(|&a, &b| merges[a].2.partial_cmp(&merges[b].2).unwrap().then(a.cmp(&b)));
    let mut parent: Vec<usize> = (0..2 * n - 1).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let mut ks_sorted = ks.to_vec();
    ks_sorted.sort_unstable_by(|a, b| b.cmp(a));
    let mut applied = 0;
    let mut out = Vec::new();
    for k in ks_sorted {
        let k = k.clamp(1, n);
        while applied < n - k {
            let (a, b, _) = merges[order[applied]];
            let new = n + order[applied];
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            parent[rb] = ra;
            let rn = find(&mut parent, new);
            parent[rn] = ra;
            applied += 1;
        }
        let mut compact: HashMap<usize, u32> = HashMap::new();
        let labels = (0..n)
            .map(|i| {
                let r = find(&mut parent, i);
                let next = compact.len() as u32;
                *compact.entry(r).or_insert(next)
            })
            .collect();
        out.push((k, labels));
    }
    out
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n == 0 {
        0.0
    } else if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// Candidate sides from an embedding: k-means clusters for several k, Ward clusters for the same
/// k (up to 6,000 cells), and the tails of each component, each as a mask (cluster vs rest).
fn candidates(z: &[f64], n: usize, d: usize, s: &DiscoverySettings, seed: u64) -> Vec<(Vec<bool>, String)> {
    let mut seen: HashSet<Vec<u64>> = HashSet::new();
    let mut out: Vec<(Vec<bool>, String)> = Vec::new();
    let mut add = |mask: Vec<bool>, how: String| {
        let inside = mask.iter().filter(|&&m| m).count();
        if inside.min(n - inside) < s.min_cells {
            return;
        }
        let mut key = vec![0_u64; n.div_ceil(64)];
        for (i, &m) in mask.iter().enumerate() {
            if m {
                key[i / 64] |= 1 << (i % 64);
            }
        }
        if seen.insert(key) {
            out.push((mask, how));
        }
    };
    let ks: Vec<usize> = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]
        .into_iter()
        .take_while(|&kk| kk <= n / s.min_cells)
        .collect();
    let ward = if n <= 6000 && !ks.is_empty() {
        ward_cuts(&ward_merges(z, n, d), n, &ks)
    } else {
        Vec::new()
    };
    for (i, &kk) in ks.iter().enumerate() {
        let lab = kmeans(z, n, d, kk, 10, seed.wrapping_add(i as u64), s.threads);
        for c in 0..kk {
            add(lab.iter().map(|&l| l as usize == c).collect(), format!("k-means k={kk}"));
        }
        if let Some((_, lab)) = ward.iter().find(|(k, _)| *k == kk) {
            for c in 0..kk {
                add(lab.iter().map(|&l| l as usize == c).collect(), format!("ward k={kk}"));
            }
        }
    }
    for j in 0..d {
        let mut col: Vec<f64> = (0..n).map(|i| z[i * d + j]).collect();
        let med = median(&mut col.clone());
        let mut dev: Vec<f64> = col.iter().map(|v| (v - med).abs()).collect();
        let mad = median(&mut dev);
        for v in col.iter_mut() {
            *v = (*v - med) / (1.4826 * mad + 1e-12);
        }
        for sign in [1.0, -1.0] {
            for cut in [3.0, 5.0] {
                add(
                    col.iter().map(|&v| sign * v > cut).collect(),
                    format!("component {j} tail {:+.0}", sign * cut),
                );
            }
        }
    }
    out
}

struct Discoverer<'a> {
    x: Csr,
    genes: Genes,
    s: &'a DiscoverySettings,
    rng: StdRng,
    splits: Vec<SplitRecord>,
}

impl Discoverer<'_> {
    /// Held-out gain of splitting a leaf into side / rest: shares from `xa` (shrunk toward the
    /// leaf's), likelihood on `xb`, against the leaf's shares from `xa`; and the genes with a
    /// finite strength.
    fn gain_of(&self, side: &[bool], xa: &Csr, xb: &Csr) -> (f64, usize) {
        let phi = self.genes.plain_shares(&xa.colsum(None));
        let (s_in, s_out, beta) = self.genes.eb_shares(
            &xa.colsum(Some((side, true))),
            &xa.colsum(Some((side, false))),
            &phi,
            self.s.gate_z,
            self.s.gene_prior_weight,
        );
        let ll_leaf = loglik_rows(xb, &phi);
        let ll_in = loglik_rows(xb, &s_in);
        let ll_out = loglik_rows(xb, &s_out);
        let gain = (0..xb.rows)
            .map(|r| if side[r] { ll_in[r] } else { ll_out[r] } - ll_leaf[r])
            .sum();
        (gain, beta.iter().filter(|b| b.is_finite()).count())
    }

    /// Tries to split the cells `rows`; the two sides when a split is accepted.
    fn split(&mut self, rows: &[usize], depth: usize) -> Option<(Vec<usize>, Vec<usize>)> {
        let n = rows.len();
        if n < 2 * self.s.min_cells || depth >= self.s.max_depth {
            return None;
        }
        let parts = thin(&self.x.select_rows(rows), 4, &mut self.rng);
        let (x1, x2, x3, x4) = (&parts[0], &parts[1], &parts[2], &parts[3]);
        let xf = x1.add(x2);
        let phi_f = self.genes.plain_shares(&xf.colsum(None));
        let (keep, n_var) = self.genes.variable_genes(&xf, &phi_f);
        let r = self.genes.residuals(&xf, &phi_f, &keep);
        let k = self.s.components.min(r.rows.min(r.cols).saturating_sub(1));
        if k < 2 || r.nnz() == 0 {
            return None;
        }
        let seed_svd: u64 = self.rng.random();
        let seed_km: u64 = self.rng.random();
        let z = embed(&r, k, seed_svd, self.s.threads);
        let cands = candidates(&z, n, k, self.s, seed_km);
        if cands.is_empty() {
            return None;
        }
        let gains = par_map(cands.len(), self.s.threads, |i| self.gain_of(&cands[i].0, &xf, x3).0);
        let best = (0..cands.len())
            .max_by(|&a, &b| gains[a].partial_cmp(&gains[b]).unwrap().then(b.cmp(&a)))
            .unwrap();
        let (side, how) = &cands[best];
        // confirmation on molecules that neither made nor chose the split: shares from parts 1-3,
        // likelihood on part 4, against random sides of the same size scored the same way
        let xa = xf.add(x3);
        let (gain, genes_finite) = self.gain_of(side, &xa, x4);
        let inside = side.iter().filter(|&&m| m).count();
        let fakes: Vec<Vec<bool>> = (0..self.s.null_sides)
            .map(|_| {
                let mut fake = vec![false; n];
                for i in rand::seq::index::sample(&mut self.rng, n, inside) {
                    fake[i] = true;
                }
                fake
            })
            .collect();
        let null = par_map(fakes.len(), self.s.threads, |i| self.gain_of(&fakes[i], &xa, x4).0);
        let null_max = null.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let per_cell = (gain - null_max) / n as f64;
        let accepted = gain > null_max + self.s.min_gain && per_cell > self.s.min_gain_per_cell;
        self.splits.push(SplitRecord {
            depth,
            cells: n,
            variable_genes: n_var,
            side: inside,
            rest: n - inside,
            how: how.clone(),
            selection_gain: gains[best],
            gain,
            null_max,
            per_cell,
            genes_finite,
            accepted,
        });
        if !accepted {
            return None;
        }
        let (mut left, mut right) = (Vec::new(), Vec::new());
        for (i, &r) in rows.iter().enumerate() {
            if side[i] {
                left.push(r);
            } else {
                right.push(r);
            }
        }
        Some((left, right))
    }
}

/// Finds groups of cells that share isoform usage from their one-candidate molecules.
/// `counts[c]` lists `(transcript id, molecules whose only candidate it is)` for cell `c`;
/// `gene_of[t]` is transcript `t`'s gene.
pub fn discover_groups(counts: &[Vec<(u32, u32)>], gene_of: &[u32], s: &DiscoverySettings) -> Discovery {
    let n_cells = counts.len();
    let one_group = |genes, transcripts, molecules| Discovery {
        labels: vec![0; n_cells],
        n_groups: 1,
        sizes: vec![n_cells],
        genes,
        transcripts,
        molecules,
        splits: Vec::new(),
        refinement: Vec::new(),
    };
    // genes with at least two isoforms seen and enough molecules; columns sorted by gene
    let mut total: HashMap<u32, f64> = HashMap::new();
    for cell in counts {
        for &(t, c) in cell {
            if c > 0 {
                *total.entry(t).or_default() += c as f64;
            }
        }
    }
    let mut per_gene: HashMap<u32, (usize, f64)> = HashMap::new();
    for (&t, &c) in &total {
        let e = per_gene.entry(gene_of[t as usize]).or_default();
        e.0 += 1;
        e.1 += c;
    }
    let mut cols: Vec<(u32, u32)> = total
        .keys()
        .filter(|&&t| {
            let (iso, mols) = per_gene[&gene_of[t as usize]];
            iso >= 2 && mols >= s.min_gene_molecules
        })
        .map(|&t| (gene_of[t as usize], t))
        .collect();
    cols.sort_unstable();
    if cols.is_empty() || n_cells == 0 {
        return one_group(0, 0, 0.0);
    }
    let col_of: HashMap<u32, u32> = cols.iter().enumerate().map(|(i, &(_, t))| (t, i as u32)).collect();
    let mut gi = Vec::with_capacity(cols.len());
    let mut start = vec![0];
    for (i, &(g, _)) in cols.iter().enumerate() {
        if i > 0 && cols[i - 1].0 != g {
            start.push(i);
        }
        gi.push((start.len() - 1) as u32);
    }
    start.push(cols.len());
    let genes = Genes { gi, start };
    let mut indptr = Vec::with_capacity(n_cells + 1);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    indptr.push(0);
    let mut molecules = 0.0;
    for cell in counts {
        let mut row: Vec<(u32, f64)> = cell
            .iter()
            .filter(|&&(_, c)| c > 0)
            .filter_map(|&(t, c)| col_of.get(&t).map(|&col| (col, c as f64)))
            .collect();
        row.sort_unstable_by_key(|&(c, _)| c);
        for (c, v) in row {
            indices.push(c);
            data.push(v);
            molecules += v;
        }
        indptr.push(indices.len());
    }
    let x = Csr {
        rows: n_cells,
        cols: cols.len(),
        indptr,
        indices,
        data,
    };
    let (n_genes, n_txps) = (genes.n(), cols.len());
    let mut d = Discoverer {
        x,
        genes,
        s,
        rng: StdRng::seed_from_u64(s.seed),
        splits: Vec::new(),
    };
    let mut leaves: Vec<(Vec<usize>, usize)> = vec![((0..n_cells).collect(), 0)];
    let mut finals: Vec<Vec<usize>> = Vec::new();
    while let Some((rows, depth)) = leaves.pop() {
        match d.split(&rows, depth) {
            Some((left, right)) => {
                leaves.push((left, depth + 1));
                leaves.push((right, depth + 1));
            }
            None => finals.push(rows),
        }
    }
    let mut labels = vec![0_u32; n_cells];
    for (k, rows) in finals.iter().enumerate() {
        for &r in rows {
            labels[r] = k as u32;
        }
    }
    let n_groups = finals.len();
    let mut refinement = Vec::new();
    if n_groups > 1 {
        // each cell to the group whose shrunk shares fit all its molecules best
        let phi_all = d.genes.plain_shares(&d.x.colsum(None));
        for _ in 0..15 {
            let ll: Vec<Vec<f64>> = par_map(n_groups, s.threads, |k| {
                let mask: Vec<bool> = labels.iter().map(|&l| l as usize == k).collect();
                let (s_in, _, _) = d.genes.eb_shares(
                    &d.x.colsum(Some((&mask, true))),
                    &d.x.colsum(Some((&mask, false))),
                    &phi_all,
                    s.gate_z,
                    s.gene_prior_weight,
                );
                loglik_rows(&d.x, &s_in)
            });
            let mut new: Vec<u32> = (0..n_cells)
                .map(|c| {
                    (0..n_groups)
                        .max_by(|&a, &b| ll[a][c].partial_cmp(&ll[b][c]).unwrap().then(b.cmp(&a)))
                        .unwrap() as u32
                })
                .collect();
            for k in 0..n_groups as u32 {
                if new.iter().filter(|&&l| l == k).count() < s.min_cells {
                    for (c, &old) in labels.iter().enumerate() {
                        if old == k {
                            new[c] = k;
                        }
                    }
                }
            }
            let moved = new.iter().zip(&labels).filter(|(a, b)| a != b).count();
            labels = new;
            refinement.push(moved);
            if moved == 0 {
                break;
            }
        }
    }
    let mut sizes = vec![0_usize; n_groups];
    for &l in &labels {
        sizes[l as usize] += 1;
    }
    Discovery {
        labels,
        n_groups,
        sizes,
        genes: n_genes,
        transcripts: n_txps,
        molecules,
        splits: d.splits,
        refinement,
    }
}
