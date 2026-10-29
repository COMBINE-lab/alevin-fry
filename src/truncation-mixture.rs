//! Per-cell positional (truncation) model for long-read molecules.
//!
//! The model is the one of oarfish's `truncation.rs` (`truncation-mixture` branch), whose
//! likelihood, laws and M-step are reproduced here unchanged: a mixture of the ways a long read
//! can sit on a transcript.
//!
//! Long-read protocols sequence one read per captured molecule, so read counts need no
//! fragmentation-style length correction. What does depend on the transcript is where a read
//! can land on it, and that differs by protocol: direct RNA and oligo(dT)-primed cDNA give reads
//! anchored at the 3' end and truncated toward the 5' end; some library preparations keep the
//! 5' end and lose the 3' end; full-length selection gives reads spanning the whole molecule;
//! and fragmentation (or a simulator that places reads uniformly) gives reads anywhere on it.
//! Which of these a sample contains is learned from the sample, as a mixture.
//!
//! For an alignment covering `[s, e)` of a transcript of length `L`, with read length
//! `l = e - s` and gaps `d5 = s`, `d3 = L - e`: each end has a kernel `K3(d3)`, `K5(d5)`, the
//! distribution of a read's gap to that end when the read sits at it. Every component draws a
//! read length from its own law, censored at the length available to it (a read that would be
//! longer runs to the far end), and places the read by its own rule. With `f_c`, `S_c` the
//! density and survival of component `c`'s read-length law (a histogram, see [`LenLaw`]):
//!
//! ```text
//! 3'-anchored  K3(d3) ( f_3(l)  +  S_3(L - d3) K5(d5) )
//! 5'-anchored  K5(d5) ( f_5(l)  +  S_5(L - d5) K3(d3) )
//! full length  K3(d3) K5(d5)
//! uniform      f_u(l) / (L - l + 1)  +  S_u(L) K3(d3) K5(d5)
//! P(s, e | t) = p3 · 3'-anchored + p5 · 5'-anchored + pf · full length + pu · uniform
//! ```
//!
//! Here the model is fitted independently to each cell: [`fit_molecules`] runs an EM over one
//! cell's deduplicated molecules, estimating the transcript abundances together with the
//! mixture weights, read-length laws and end kernels (the latter by the truncated- and
//! censored-data M-steps), and then replaces each molecule's transcript probabilities by its
//! posterior. A molecule's likelihood under a transcript is its prior probability for that
//! transcript times `P(s, e | t)`. Cells are processed one at a time, so the EM here is
//! sequential (the caller parallelizes across cells).
//!
//! Alternatively ([`TruncationParamScope::Global`]) the parameters are estimated once, by
//! [`fit_global`], from the molecules of all cells (each cell with its own abundances), and
//! every cell's EM ([`fit_molecules_with`]) then estimates only its abundances with them fixed,
//! as oarfish does after estimating the parameters on a sample of the reads.

// The statistics merge helpers are kept as in oarfish, where they combine per-thread results.
#![allow(dead_code)]

use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

use serde::Serialize;

use crate::eq_class::merge_intervals;
use crate::pugutils_dp_new::ResolvedMolecule;

/// What the positional model needs of one alignment: its target and score probability, the
/// target's length, the aligned length, and the gaps to the transcript's 5' and 3' ends (the
/// span clamped to the transcript).
#[derive(Clone, Copy)]
struct Geo {
    t: u32,
    p: f64,
    len: u32,
    l: u32,
    d5: u32,
    d3: u32,
}

/// An end kernel: the distribution of the gap `d` (nt) between a read's end and the transcript's
/// end, for reads placed at that end. Nonparametric, on log-spaced bins `{0}`, `{1}`, `[2, 3]`,
/// `[4, 7]`, ... (uniform within a bin), so it can represent both reads flush with the end and
/// ends spread hundreds of bases from it (alternative polyadenylation, annotation mismatch).
/// For a transcript of length `L` it is normalized over the gaps the transcript allows,
/// `0..L`, and fitted by the truncated-data EM.
#[derive(Clone)]
struct GapLaw {
    h: Vec<f64>,
}

impl std::fmt::Debug for GapLaw {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.summary();
        write!(
            f,
            "GapLaw(0: {:.2}, 1-15: {:.2}, 16-255: {:.2}, 256+: {:.2})",
            s[0], s[1], s[2], s[3]
        )
    }
}

impl GapLaw {
    #[inline]
    fn bin(d: u32) -> usize {
        (32 - d.leading_zeros()) as usize
    }

    #[inline]
    fn width(b: usize) -> f64 {
        if b == 0 {
            1.0
        } else {
            (1u64 << (b - 1)) as f64
        }
    }

    /// Discretized exponential law of mean 100 nt over gaps up to `max_len`.
    fn initial(max_len: u32) -> Self {
        let nb = Self::bin(max_len) + 1;
        let cdf = |x: f64| 1.0 - (-x / 100.0).exp();
        let h: Vec<f64> = (0..nb)
            .map(|b| {
                let lo = if b == 0 { 0.0 } else { Self::width(b) };
                cdf(lo + Self::width(b)) - cdf(lo)
            })
            .collect();
        let s: f64 = h.iter().sum();
        GapLaw {
            h: h.into_iter().map(|x| x / s).collect(),
        }
    }

    /// Per-nt density `dens[d]` and distribution function `cdf[x] = P(gap <= x)` for gaps up
    /// to `max_len`.
    fn tables(&self, max_len: u32) -> (Vec<f64>, Vec<f64>) {
        let dens: Vec<f64> = (0..=max_len)
            .map(|d| {
                let b = Self::bin(d);
                self.h.get(b).map_or(0.0, |p| p / Self::width(b))
            })
            .collect();
        let mut acc = 0.0;
        let cdf = dens
            .iter()
            .map(|x| {
                acc += x;
                acc
            })
            .collect();
        (dens, cdf)
    }

    /// Mass at gap 0, in 1-15, 16-255 and 256 or more.
    fn summary(&self) -> [f64; 4] {
        let mut s = [0.0; 4];
        for (b, p) in self.h.iter().enumerate() {
            s[match b {
                0 => 0,
                1..=4 => 1,
                5..=8 => 2,
                _ => 3,
            }] += p;
        }
        s
    }

    fn moved(&self, o: &GapLaw) -> f64 {
        0.5 * self
            .h
            .iter()
            .zip(&o.h)
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
    }

    /// Truncated-data EM step from `old`: `s.stop[b]` is the observed weight in bin `b`,
    /// `s.cens[c]` the weight of observations truncated at a gap in bin `c` (the transcript
    /// allowed no larger gap); each implies unobserved draws beyond its truncation point,
    /// spread over the bins above `c` in proportion to `old.h`.
    fn law(s: &LStats, old: &GapLaw) -> GapLaw {
        let nb = old.h.len();
        let mut n = vec![0.0; nb];
        let (mut f, mut acc) = (0.0, 0.0);
        for (b, (nb_, &hb)) in n.iter_mut().zip(&old.h).enumerate() {
            *nb_ = s.stop.get(b).copied().unwrap_or(0.0) + acc * hb;
            f += hb;
            let c = s.cens.get(b).copied().unwrap_or(0.0);
            if c > 0.0 && f > 0.0 {
                // (1 - F) / F unobserved draws per truncated one, spread as h_b / (1 - F)
                acc += c / f;
            }
        }
        let tot: f64 = n.iter().sum();
        if tot <= 0.0 {
            return old.clone();
        }
        let eps = 1e-9 * tot;
        let z = tot + eps * nb as f64;
        GapLaw {
            h: n.into_iter().map(|x| (x + eps) / z).collect(),
        }
    }
}

/// Bins of the relative part of an end kernel (gap as a fraction of transcript length).
const RB: usize = 20;

/// An end kernel: how far a read placed at an end sits from it. A mixture of an absolute part
/// (a [`GapLaw`] in nucleotides, for reads flush with the end or a fixed distance from it) and a
/// relative part (a histogram over the gap as a fraction of the transcript's length, for ends
/// that fall at a relative position, as poly(A) sites in long 3' UTRs do), each a distribution
/// over the gaps `0..L` a transcript of length `L` allows: `K(d | L) = (1 - rho) K_abs(d) /
/// C_abs(L - 1) + rho g(d / L) / L`.
#[derive(Clone)]
struct EndLaw {
    abs: GapLaw,
    rel: Vec<f64>,
    rho: f64,
}

impl std::fmt::Debug for EndLaw {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let near: f64 = self.rel[..RB / 10].iter().sum();
        let far: f64 = self.rel[RB / 2..].iter().sum();
        write!(
            f,
            "{:?} + relative {:.2} (gap/L < 0.1: {:.2}, >= 0.5: {:.2})",
            self.abs, self.rho, near, far
        )
    }
}

impl EndLaw {
    fn initial(max_len: u32) -> Self {
        EndLaw {
            abs: GapLaw::initial(max_len),
            rel: vec![1.0 / RB as f64; RB],
            rho: 0.5,
        }
    }

    #[inline]
    fn rbin(d: u32, len: u32) -> usize {
        ((d as u64 * RB as u64) / len.max(1) as u64).min(RB as u64 - 1) as usize
    }

    /// The absolute and relative parts of `K(d | len)`, given the absolute part's tables.
    #[inline]
    fn parts(&self, d: u32, len: u32, dens: &[f64], cdf: &[f64]) -> (f64, f64) {
        let last = len.saturating_sub(1) as usize;
        let a = (1.0 - self.rho) * dens[d as usize] / cdf[last].max(1e-300);
        let r = self.rho * self.rel[Self::rbin(d, len)] * RB as f64 / len.max(1) as f64;
        (a, r)
    }

    fn moved(&self, o: &EndLaw) -> f64 {
        let tv: f64 = 0.5
            * self
                .rel
                .iter()
                .zip(&o.rel)
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>();
        self.abs.moved(&o.abs).max(tv).max((self.rho - o.rho).abs())
    }

    fn law(s: &EStats, old: &EndLaw) -> EndLaw {
        let n = s.n_abs + s.n_rel;
        if n <= 0.0 {
            return old.clone();
        }
        let tot: f64 = s.rel.iter().sum();
        let eps = 1e-9 * tot.max(1e-300);
        let z = tot + eps * RB as f64;
        EndLaw {
            abs: GapLaw::law(&s.abs, &old.abs),
            rel: if tot > 0.0 {
                (0..RB)
                    .map(|b| (s.rel.get(b).copied().unwrap_or(0.0) + eps) / z)
                    .collect()
            } else {
                old.rel.clone()
            },
            rho: (s.n_rel / n).clamp(1e-6, 1.0 - 1e-6),
        }
    }
}

/// Expected sufficient statistics of an end kernel.
#[derive(Default, Clone)]
struct EStats {
    abs: LStats,
    rel: Vec<f64>,
    n_abs: f64,
    n_rel: f64,
}

impl EStats {
    fn merge(mut self, o: EStats) -> EStats {
        self.abs = self.abs.merge(o.abs);
        if self.rel.len() < o.rel.len() {
            self.rel.resize(o.rel.len(), 0.0);
        }
        self.rel.iter_mut().zip(o.rel).for_each(|(a, b)| *a += b);
        self.n_abs += o.n_abs;
        self.n_rel += o.n_rel;
        self
    }

    /// Weight `w` of a gap `d` on a transcript of length `len`, split between the kernel's
    /// parts in proportion to their values `a`, `r` there.
    fn add(&mut self, w: f64, d: u32, len: u32, a: f64, r: f64) {
        let z = a + r;
        if w <= 0.0 || z <= 0.0 {
            return;
        }
        let (wa, wr) = (w * a / z, w * r / z);
        LStats::add(&mut self.abs.stop, GapLaw::bin(d), wa);
        LStats::add(&mut self.abs.cens, GapLaw::bin(len.saturating_sub(1)), wa);
        LStats::add(&mut self.rel, EndLaw::rbin(d, len), wr);
        self.n_abs += wa;
        self.n_rel += wr;
    }
}

/// A read-length law: a histogram over lengths in bins of [`LenLaw::W`] nt, uniform within a
/// bin. Lengths are right-censored at the length a transcript makes available (a read that
/// would be longer runs to the transcript's end), so the law is fitted by the censored-data EM
/// (each censored read is spread over the bins at or above its censoring length).
#[derive(Clone)]
struct LenLaw {
    h: Vec<f64>,
}

impl std::fmt::Debug for LenLaw {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "LenLaw(mean {:.0} nt)", self.mean())
    }
}

impl LenLaw {
    const W: f64 = 50.0;

    /// Discretized exponential law with the given mean, over lengths up to `max_len`.
    fn exponential(mean: f64, max_len: u32) -> Self {
        let nb = (max_len as f64 / Self::W) as usize + 1;
        let h: Vec<f64> = (0..nb)
            .map(|b| {
                let (a, z) = (b as f64 * Self::W, (b + 1) as f64 * Self::W);
                (-a / mean).exp() - (-z / mean).exp()
            })
            .collect();
        let s: f64 = h.iter().sum();
        LenLaw {
            h: h.into_iter().map(|x| x / s).collect(),
        }
    }

    #[inline]
    fn bin(x: f64) -> usize {
        (x / Self::W) as usize
    }

    /// [`LenLaw::bin`] of a whole number of bases, by integer division (the same bin: `x / W`
    /// is at least `1 / W` from the next whole number, far beyond rounding).
    #[inline]
    fn bin_nt(x: u32) -> usize {
        const W_NT: u32 = LenLaw::W as u32;
        (x / W_NT) as usize
    }

    fn mean(&self) -> f64 {
        self.h
            .iter()
            .enumerate()
            .map(|(b, p)| p * (b as f64 + 0.5) * Self::W)
            .sum()
    }

    /// Per-nt density `dens[x]` and survival `surv[x] = P(length >= x)` for `x` up to `max_len`.
    fn tables(&self, max_len: u32) -> (Vec<f64>, Vec<f64>) {
        let nb = self.h.len();
        let mut tail = vec![0.0; nb + 1];
        for b in (0..nb).rev() {
            tail[b] = tail[b + 1] + self.h[b];
        }
        let dens = (0..=max_len)
            .map(|x| self.h.get(Self::bin(x as f64)).map_or(0.0, |p| p / Self::W))
            .collect();
        let surv = (0..=max_len)
            .map(|x| {
                let b = Self::bin(x as f64);
                if b >= nb {
                    return 0.0;
                }
                let within = (x as f64 - b as f64 * Self::W) / Self::W;
                (tail[b + 1] + self.h[b] * (1.0 - within)).max(0.0)
            })
            .collect();
        (dens, surv)
    }

    /// Total-variation distance between two laws.
    fn moved(&self, o: &LenLaw) -> f64 {
        0.5 * self
            .h
            .iter()
            .zip(&o.h)
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
    }
}

/// Expected sufficient statistics of a length law: observed lengths and censoring lengths,
/// binned.
#[derive(Default, Clone)]
struct LStats {
    stop: Vec<f64>,
    cens: Vec<f64>,
}

impl LStats {
    #[inline]
    fn add(v: &mut Vec<f64>, b: usize, w: f64) {
        if w <= 0.0 {
            return;
        }
        if v.len() <= b {
            v.resize(b + 1, 0.0);
        }
        v[b] += w;
    }

    fn merge(mut self, o: LStats) -> LStats {
        for (dst, src) in [(&mut self.stop, o.stop), (&mut self.cens, o.cens)] {
            if dst.len() < src.len() {
                dst.resize(src.len(), 0.0);
            }
            dst.iter_mut().zip(src).for_each(|(a, b)| *a += b);
        }
        self
    }

    /// Censored-histogram EM step from `old`: an observed length counts in its bin; a read
    /// censored in bin `c` is spread over the bins `b >= c` in proportion to `old.h[b]`.
    fn law(&self, old: &LenLaw) -> LenLaw {
        let nb = old.h.len();
        let mut tail = vec![0.0; nb + 1];
        for b in (0..nb).rev() {
            tail[b] = tail[b + 1] + old.h[b];
        }
        let mut n = vec![0.0; nb];
        // running sum over c <= b of cens[c] / P(length in bins >= c)
        let mut acc = 0.0;
        for b in 0..nb {
            let c = self.cens.get(b).copied().unwrap_or(0.0);
            if c > 0.0 && tail[b] > 0.0 {
                acc += c / tail[b];
            }
            n[b] = self.stop.get(b).copied().unwrap_or(0.0) + acc * old.h[b];
        }
        let tot: f64 = n.iter().sum();
        if tot <= 0.0 {
            return old.clone();
        }
        // A small floor keeps every length possible.
        let eps = 1e-9 * tot;
        let z = tot + eps * nb as f64;
        LenLaw {
            h: n.into_iter().map(|x| (x + eps) / z).collect(),
        }
    }
}

/// The positional model's parameters.
#[derive(Clone, Debug)]
struct Params {
    /// Mixture weights: 3'-anchored, 5'-anchored, full length, uniformly placed.
    p3: f64,
    p5: f64,
    pf: f64,
    pu: f64,
    /// Read-length laws of the 3'-anchored, 5'-anchored and uniform components.
    len3: LenLaw,
    len5: LenLaw,
    lenu: LenLaw,
    /// End kernels at the 3' and 5' ends.
    k3: EndLaw,
    k5: EndLaw,
}

/// Per-alignment terms of the positional likelihood, each already weighted by its mixture
/// weight. `*_stop` is the part where the read's length was observed, `*_cens` the part where
/// it spans the transcript (length censored at `L`).
struct Terms {
    a3_stop: f64,
    a3_cens: f64,
    a5_stop: f64,
    a5_cens: f64,
    full: f64,
    u_stop: f64,
    u_cens: f64,
    /// Absolute and relative parts of the end kernels at this alignment's gaps.
    k3a: f64,
    k3r: f64,
    k5a: f64,
    k5r: f64,
}

impl Terms {
    #[inline]
    fn total(&self) -> f64 {
        self.a3_stop
            + self.a3_cens
            + self.a5_stop
            + self.a5_cens
            + self.full
            + self.u_stop
            + self.u_cens
    }
}

impl Params {
    fn initial(max_len: u32) -> Self {
        Params {
            p3: 0.25,
            p5: 0.25,
            pf: 0.25,
            pu: 0.25,
            len3: LenLaw::exponential(1000.0, max_len),
            len5: LenLaw::exponential(1000.0, max_len),
            lenu: LenLaw::exponential(1000.0, max_len),
            k3: EndLaw::initial(max_len),
            k5: EndLaw::initial(max_len),
        }
    }

    #[inline]
    fn terms(&self, g: &Geo, e: &Exps) -> Terms {
        let (len, l, d5, d3) = (g.len, g.l, g.d5, g.d3);
        let (l_us, len_us) = (l as usize, len as usize);
        // end kernels, each a distribution over the gaps this transcript allows (0..len)
        let (k3a, k3r) = self.k3.parts(d3, len, &e.g3d, &e.g3c);
        let (k5a, k5r) = self.k5.parts(d5, len, &e.g5d, &e.g5c);
        let (kk3, kk5) = (k3a + k3r, k5a + k5r);
        Terms {
            a3_stop: self.p3 * kk3 * e.d3[l_us],
            a3_cens: self.p3 * kk3 * e.s3[(len - d3) as usize] * kk5,
            a5_stop: self.p5 * kk5 * e.d5[l_us],
            a5_cens: self.p5 * kk5 * e.s5[(len - d5) as usize] * kk3,
            full: self.pf * kk3 * kk5,
            u_stop: self.pu * e.du[l_us] / ((len - l) as f64 + 1.0),
            u_cens: self.pu * e.su[len_us] * kk3 * kk5,
            k3a,
            k3r,
            k5a,
            k5r,
        }
    }
}

/// `exp(-rate x)` for each rate, and `(1 - rho)^x` for each kernel scale, for every whole number
/// of bases `x` up to the longest transcript, built once per iteration so that no alignment
/// computes an exponential.
struct Exps {
    /// Length-law densities and survivals of the three components.
    d3: Vec<f64>,
    s3: Vec<f64>,
    d5: Vec<f64>,
    s5: Vec<f64>,
    du: Vec<f64>,
    su: Vec<f64>,
    /// End-kernel densities and distribution functions.
    g3d: Vec<f64>,
    g3c: Vec<f64>,
    g5d: Vec<f64>,
    g5c: Vec<f64>,
}

impl Exps {
    fn new(par: &Params, max_len: u32) -> Self {
        let (d3, s3) = par.len3.tables(max_len);
        let (d5, s5) = par.len5.tables(max_len);
        let (du, su) = par.lenu.tables(max_len);
        let (g3d, g3c) = par.k3.abs.tables(max_len);
        let (g5d, g5c) = par.k5.abs.tables(max_len);
        Exps {
            d3,
            s3,
            d5,
            s5,
            du,
            su,
            g3d,
            g3c,
            g5d,
            g5c,
        }
    }
}

/// Expected sufficient statistics for the M-step, accumulated per thread.
#[derive(Default)]
struct Stats {
    /// Posterior mass, and its share in each component.
    n: f64,
    n3: f64,
    n5: f64,
    nf: f64,
    nu: f64,
    /// Length laws.
    l3: LStats,
    l5: LStats,
    lu: LStats,
    /// End kernels: observed gaps (`stop`) and truncation points (`cens`), binned.
    k3: EStats,
    k5: EStats,
    /// Log-likelihood of the reads at the step's input (abundances unnormalized).
    ll: f64,
}

impl Stats {
    fn merge(mut self, o: Stats) -> Stats {
        self.ll += o.ll;
        self.n += o.n;
        self.n3 += o.n3;
        self.n5 += o.n5;
        self.nf += o.nf;
        self.nu += o.nu;
        self.l3 = self.l3.merge(o.l3);
        self.l5 = self.l5.merge(o.l5);
        self.lu = self.lu.merge(o.lu);
        self.k3 = self.k3.merge(o.k3);
        self.k5 = self.k5.merge(o.k5);
        self
    }

    fn add(&mut self, g: &Geo, t: &Terms, post: f64) {
        let tot = t.total();
        if tot <= 0.0 {
            return;
        }
        let q = post / tot;
        let (a3s, a3c) = (q * t.a3_stop, q * t.a3_cens);
        let (a5s, a5c) = (q * t.a5_stop, q * t.a5_cens);
        let (us, uc) = (q * t.u_stop, q * t.u_cens);
        let fl = q * t.full;
        self.n += post;
        self.n3 += a3s + a3c;
        self.n5 += a5s + a5c;
        self.nf += fl;
        self.nu += us + uc;
        // A read reaching the far end is right-censored at the length available to it; any
        // other read is an observed stop after `l` bases.
        LStats::add(&mut self.l3.stop, LenLaw::bin_nt(g.l), a3s);
        LStats::add(&mut self.l3.cens, LenLaw::bin_nt(g.len - g.d3), a3c);
        LStats::add(&mut self.l5.stop, LenLaw::bin_nt(g.l), a5s);
        LStats::add(&mut self.l5.cens, LenLaw::bin_nt(g.len - g.d5), a5c);
        LStats::add(&mut self.lu.stop, LenLaw::bin_nt(g.l), us);
        LStats::add(&mut self.lu.cens, LenLaw::bin_nt(g.len), uc);
        // Every part that places the read at an end through that end's kernel observes its gap.
        let (w3, w5) = (a3s + a3c + a5c + fl + uc, a5s + a5c + a3c + fl + uc);
        let lu = g.len;
        self.k3.add(w3, g.d3, lu, t.k3a, t.k3r);
        self.k5.add(w5, g.d5, lu, t.k5a, t.k5r);
    }

    fn params(&self, old: &Params) -> Params {
        if self.n <= 0.0 {
            return old.clone();
        }
        let frac = |x: f64| (x / self.n).clamp(1e-6, 1.0);
        let (p3, p5, pf, pu) = (frac(self.n3), frac(self.n5), frac(self.nf), frac(self.nu));
        let z = p3 + p5 + pf + pu;
        Params {
            p3: p3 / z,
            p5: p5 / z,
            pf: pf / z,
            pu: pu / z,
            len3: self.l3.law(&old.len3),
            len5: self.l5.law(&old.len5),
            lenu: self.lu.law(&old.lenu),
            k3: EndLaw::law(&self.k3, &old.k3),
            k5: EndLaw::law(&self.k5, &old.k5),
        }
    }
}

/// How far the parameters moved: absolute change for the mixture weights; for each length law,
/// its change scaled by the weight of the one component that uses it, as that bounds its effect
/// on the likelihood (the law of a component holding a negligible share of the reads drifts for
/// as long as the EM runs, and must not hold the parameters open). The end kernels enter
/// several components and count in full.
fn moved(a: &Params, b: &Params) -> f64 {
    (a.p3 - b.p3)
        .abs()
        .max((a.p5 - b.p5).abs())
        .max((a.pf - b.pf).abs())
        .max((a.pu - b.pu).abs())
        .max(a.p3.max(b.p3) * a.len3.moved(&b.len3))
        .max(a.p5.max(b.p5) * a.len5.moved(&b.len5))
        .max(a.pu.max(b.pu) * a.lenu.moved(&b.lenu))
        .max(a.k3.moved(&b.k3))
        .max(a.k5.moved(&b.k5))
}

/// The EM's limits and thresholds, at oarfish's defaults (`--max-em-iter`,
/// `--convergence-thresh`, and its `util::constants`).
const MAX_EM_ITER: u32 = 1000;
const CONVERGENCE_THRESH: f64 = 1e-3;
const MIN_READ_THRESH: f64 = 1e-5;
const EM_DENOM_THRESH: f64 = 1e-30;

/// Parameters stop being re-estimated once none moves by more than `FREEZE_TOL` (see [`moved`])
/// for `FREEZE_AFTER` consecutive iterations, or once an iteration raises the log-likelihood by
/// less than `FREEZE_LL_GAIN` nats that many times in a row; the remaining iterations use fixed
/// weights and SQUAREM acceleration.
const FREEZE_TOL: f64 = 1e-6;
const FREEZE_AFTER: u32 = 3;
const FREEZE_LL_GAIN: f64 = 1.0;

/// Plain EM steps before SQUAREM when a cell's EM starts with the parameters fixed (estimated
/// from all cells), as oarfish does after estimating them on a read sample: without them the
/// estimates of transcripts with few molecules have not settled within the iteration limit.
const PLAIN_FIRST: u32 = 100;

/// How the several intervals a molecule has on one transcript enter the model.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub enum TruncationIntervalMode {
    /// Every interval is an alignment of its own, sharing the transcript's prior equally.
    #[default]
    All,
    /// The intervals are merged into one, from the first start to the last end.
    Union,
}

impl fmt::Display for TruncationIntervalMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::All => "all",
                Self::Union => "union",
            }
        )
    }
}

impl FromStr for TruncationIntervalMode {
    type Err = &'static str;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "all" => Ok(Self::All),
            "union" => Ok(Self::Union),
            _ => Err("expected all or union"),
        }
    }
}

/// Where the positional parameters are estimated.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub enum TruncationParamScope {
    /// Separately for each cell, together with its abundances.
    #[default]
    Cell,
    /// Once from the molecules of all cells ([`fit_global`]), then fixed in every cell's EM.
    Global,
}

impl fmt::Display for TruncationParamScope {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Cell => "cell",
                Self::Global => "global",
            }
        )
    }
}

impl FromStr for TruncationParamScope {
    type Err = &'static str;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "cell" => Ok(Self::Cell),
            "global" => Ok(Self::Global),
            _ => Err("expected cell or global"),
        }
    }
}

/// Which test froze the parameters, over the [`FREEZE_AFTER`] consecutive iterations that did.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FreezeReason {
    /// The parameters moved by less than [`FREEZE_TOL`] in every one of them.
    Params,
    /// The log-likelihood rose by less than [`FREEZE_LL_GAIN`] in every one of them.
    LogLik,
    /// Both held in every one of them.
    Both,
    /// Each passed one test or the other, but neither test held throughout.
    Mixed,
    /// Estimated from all cells and fixed from the start (`--truncation-params global`).
    Global,
}

impl fmt::Display for FreezeReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Params => "params",
                Self::LogLik => "loglik",
                Self::Both => "both",
                Self::Mixed => "mixed",
                Self::Global => "global",
            }
        )
    }
}

/// The freeze rule: counts the consecutive parameter updates that pass either test, and which
/// test has held throughout the current run of them. The log-likelihood test's threshold is
/// [`FREEZE_LL_GAIN`] times the share of the data the fit sees (1 unless it is a sample).
struct Freeze {
    ll_gain_tol: f64,
    still: u32,
    by_params: bool,
    by_loglik: bool,
}

impl Freeze {
    fn new(share: f64) -> Self {
        Freeze {
            ll_gain_tol: FREEZE_LL_GAIN * share,
            still: 0,
            by_params: true,
            by_loglik: true,
        }
    }

    /// Records one parameter update that moved the parameters by `moved` (see [`moved`]) and
    /// raised the log-likelihood by `ll_gain`; returns the reason once they should freeze.
    fn update(&mut self, moved: f64, ll_gain: f64) -> Option<FreezeReason> {
        let small_move = moved < FREEZE_TOL;
        let small_gain = (0.0..self.ll_gain_tol).contains(&ll_gain);
        if !(small_move || small_gain) {
            self.still = 0;
            return None;
        }
        if self.still == 0 {
            (self.by_params, self.by_loglik) = (true, true);
        }
        self.by_params &= small_move;
        self.by_loglik &= small_gain;
        self.still += 1;
        (self.still >= FREEZE_AFTER).then_some(match (self.by_params, self.by_loglik) {
            (true, true) => FreezeReason::Both,
            (true, false) => FreezeReason::Params,
            (false, true) => FreezeReason::LogLik,
            (false, false) => FreezeReason::Mixed,
        })
    }
}

/// A summary of one cell's fitted model, for reporting.
#[derive(Clone, Copy, Debug)]
pub struct TruncationParams {
    /// Mixture weights: 3'-anchored, 5'-anchored, full length, uniformly placed.
    pub p3: f64,
    pub p5: f64,
    pub pf: f64,
    pub pu: f64,
    /// Mean read length (nt) of the 3'-anchored, 5'-anchored and uniform components.
    pub len3_mean: f64,
    pub len5_mean: f64,
    pub lenu_mean: f64,
    /// Weight of the relative (fraction of transcript length) part of each end kernel.
    pub k3_rho: f64,
    pub k5_rho: f64,
    /// Absolute part of each end kernel: mass at gap 0, in 1-15, 16-255 and 256 or more.
    pub k3_gap: [f64; 4],
    pub k5_gap: [f64; 4],
}

impl TruncationParams {
    fn new(p: &Params) -> Self {
        Self {
            p3: p.p3,
            p5: p.p5,
            pf: p.pf,
            pu: p.pu,
            len3_mean: p.len3.mean(),
            len5_mean: p.len5.mean(),
            lenu_mean: p.lenu.mean(),
            k3_rho: p.k3.rho,
            k5_rho: p.k5.rho,
            k3_gap: p.k3.abs.summary(),
            k5_gap: p.k5.abs.summary(),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct TruncationFit {
    pub params: TruncationParams,
    pub iterations: usize,
    /// The iteration after which the parameters stopped being updated (0 if they were fixed
    /// from the start); `None` if they never froze.
    pub frozen_at: Option<usize>,
    /// Which test froze the parameters; `None` if they never froze.
    pub freeze_reason: Option<FreezeReason>,
    pub converged: bool,
    pub molecule_count: usize,
    pub has_positional_evidence: bool,
    /// Molecules with candidate transcripts both with and without a usable interval. The latter
    /// carry no positional factor, so they are not on the scale of the former and take nearly
    /// all of the molecule's mass; nonzero means the input lacks intervals it should have.
    pub mixed_interval_molecules: usize,
}

/// Positional parameters estimated once from the molecules of many cells ([`fit_global`]), for
/// every cell's EM to use as they are ([`fit_molecules_with`]).
#[derive(Clone, Debug)]
pub struct GlobalFit {
    par: Params,
    pub params: TruncationParams,
    pub cells: usize,
    pub molecules: usize,
    pub iterations: usize,
    /// The iteration after which the parameters stopped being updated; `None` if they were still
    /// moving at the iteration limit (they are then those of the last iteration).
    pub frozen_at: Option<usize>,
    pub freeze_reason: Option<FreezeReason>,
}

/// One molecule's candidate placement on one transcript. `span` is `(start, end, length)`;
/// without it (no interval is known) the placement carries no positional information.
#[derive(Clone, Copy)]
struct Alignment {
    target: usize,
    prior: f64,
    span: Option<(u32, u32, u32)>,
}

impl Alignment {
    fn geo(&self) -> Option<Geo> {
        self.span.map(|(start, end, len)| Geo {
            t: self.target as u32,
            p: self.prior,
            len,
            l: end - start,
            d5: start,
            d3: len - end,
        })
    }

    /// The abundance-free part of the weight: prior times positional likelihood.
    fn fixed_weight(&self, par: &Params, e: &Exps) -> f64 {
        match self.geo() {
            Some(g) => self.prior * par.terms(&g, e).total(),
            None => self.prior,
        }
    }

    /// The weight of this alignment under abundances `ab` and the current parameters, with the
    /// terms that produced it.
    fn weight(&self, ab: &[f64], par: &Params, e: &Exps) -> (f64, Option<(Geo, Terms)>) {
        match self.geo() {
            Some(g) => {
                let t = par.terms(&g, e);
                (ab[self.target] * self.prior * t.total(), Some((g, t)))
            }
            None => (ab[self.target] * self.prior, None),
        }
    }
}

/// One cell's molecules as the EM reads them: molecule `m` owns alignments `off[m]..off[m + 1]`,
/// whose targets index the cell's transcripts `0..n_targets`.
pub struct CellMolecules {
    off: Vec<u32>,
    alns: Vec<Alignment>,
    n_targets: usize,
}

impl CellMolecules {
    pub fn molecule_count(&self) -> usize {
        self.off.len() - 1
    }

    /// The longest transcript any of the cell's intervals lies on (1 if there are none).
    pub fn max_len(&self) -> u32 {
        self.alns
            .iter()
            .filter_map(|a| a.span.map(|(_, _, len)| len))
            .max()
            .unwrap_or(1)
    }

    fn molecule(&self, m: usize) -> &[Alignment] {
        &self.alns[self.off[m] as usize..self.off[m + 1] as usize]
    }
}

/// One cell's molecules as the EM reads them, its transcript ids (target `i` is `ids[i]`), and
/// the number of molecules mixing candidates with and without an interval (see
/// [`TruncationFit::mixed_interval_molecules`]).
fn build(
    molecules: &[ResolvedMolecule],
    interval_mode: TruncationIntervalMode,
) -> (CellMolecules, Vec<u32>, usize) {
    let mut target_ids = Vec::new();
    let mut target_indices = HashMap::new();
    let mut off = Vec::with_capacity(molecules.len() + 1);
    let mut alns = Vec::new();
    let mut mixed_interval_molecules = 0;
    off.push(0);

    for molecule in molecules {
        let (mut with_span, mut without_span) = (false, false);
        for &(target_id, prior) in &molecule.transcript_probs {
            assert!(
                prior.is_finite() && prior >= 0.0,
                "truncation model received an invalid transcript probability {prior} for transcript {target_id}"
            );
            // The softmax over log-likelihoods underflows to exactly zero for candidates far
            // below the best one; such a candidate can never receive mass in the EM.
            if prior == 0.0 {
                continue;
            }
            let target = *target_indices.entry(target_id).or_insert_with(|| {
                let index = target_ids.len();
                target_ids.push(target_id);
                index
            });
            if let Some(intervals) = molecule.intervals.get(&target_id) {
                let valid_intervals: Vec<_> = intervals
                    .iter()
                    .filter(|interval| interval.transcript_len > 0)
                    .collect();
                if !valid_intervals.is_empty() {
                    let transcript_len = valid_intervals[0].transcript_len;
                    let positions = match interval_mode {
                        TruncationIntervalMode::All => valid_intervals
                            .iter()
                            .map(|interval| (interval.start, interval.end))
                            .collect(),
                        TruncationIntervalMode::Union => {
                            let union = merge_intervals(
                                valid_intervals
                                    .iter()
                                    .map(|interval| (interval.start, interval.end))
                                    .collect(),
                            );
                            let start = union.first().map_or(0, |interval| interval.0);
                            let end = union.last().map_or(start, |interval| interval.1);
                            vec![(start, end)]
                        }
                    };
                    let interval_prior = prior / positions.len() as f64;
                    for (start, end) in positions {
                        let len = transcript_len;
                        let start = start.min(len.saturating_sub(1));
                        let end = end.clamp(start, len);
                        alns.push(Alignment {
                            target,
                            prior: interval_prior,
                            span: Some((start, end, len)),
                        });
                    }
                    with_span = true;
                    continue;
                }
            }
            without_span = true;
            alns.push(Alignment {
                target,
                prior,
                span: None,
            });
        }
        if with_span && without_span {
            mixed_interval_molecules += 1;
        }
        off.push(u32::try_from(alns.len()).expect("too many alignments in one cell"));
    }

    // Cells are held by the thousands for a global fit, so no spare capacity is kept.
    off.shrink_to_fit();
    alns.shrink_to_fit();
    let n_targets = target_ids.len();
    (
        CellMolecules {
            off,
            alns,
            n_targets,
        },
        target_ids,
        mixed_interval_molecules,
    )
}

/// One cell's molecules in the form [`fit_global`] takes.
pub fn cell_molecules(
    molecules: &[ResolvedMolecule],
    interval_mode: TruncationIntervalMode,
) -> CellMolecules {
    build(molecules, interval_mode).0
}

/// The frozen-parameter EM's inputs: each alignment's target and fixed weight (prior times
/// positional likelihood). Molecule `m` owns entries `off[m]..off[m + 1]`.
struct Fixed {
    off: Vec<usize>,
    t: Vec<usize>,
    w: Vec<f64>,
}

impl Fixed {
    fn new(cell: &CellMolecules, par: &Params, e: &Exps) -> Self {
        Fixed {
            off: cell.off.iter().map(|&o| o as usize).collect(),
            t: cell.alns.iter().map(|a| a.target).collect(),
            w: cell.alns.iter().map(|a| a.fixed_weight(par, e)).collect(),
        }
    }

    fn molecules(&self) -> usize {
        self.off.len() - 1
    }
}

/// One EM round with the parameters re-estimated: distributes each molecule over its
/// alignments, writing the abundances into `curr` and returning the positional statistics.
fn step_with(
    cell: &CellMolecules,
    par: &Params,
    e: &Exps,
    prev: &[f64],
    curr: &mut [f64],
) -> Stats {
    curr.fill(0.0);
    let mut st = Stats::default();
    let mut scored: Vec<(f64, Option<(Geo, Terms)>)> = Vec::new();
    for m in 0..cell.molecule_count() {
        let molecule = cell.molecule(m);
        scored.clear();
        let mut denom = 0.0;
        for a in molecule {
            let entry = a.weight(prev, par, e);
            denom += entry.0;
            scored.push(entry);
        }
        if denom > EM_DENOM_THRESH {
            st.ll += denom.ln();
            for (a, (w, positional)) in molecule.iter().zip(&scored) {
                let post = w / denom;
                curr[a.target] += post;
                if let Some((g, t)) = positional {
                    st.add(g, t, post);
                }
            }
        }
    }
    st
}

/// One EM round with frozen parameters: no terms or statistics are computed.
fn step_static(fx: &Fixed, prev: &[f64], curr: &mut [f64]) {
    curr.fill(0.0);
    for m in 0..fx.molecules() {
        let (lo, hi) = (fx.off[m], fx.off[m + 1]);
        let denom: f64 = (lo..hi).map(|k| prev[fx.t[k]] * fx.w[k]).sum();
        if denom > EM_DENOM_THRESH {
            for k in lo..hi {
                curr[fx.t[k]] += prev[fx.t[k]] * fx.w[k] / denom;
            }
        }
    }
}

/// Log-likelihood of abundances `theta` (normalized to proportions) under the fixed weights,
/// up to a constant: the objective the frozen-weight EM ascends.
fn loglik_static(fx: &Fixed, theta: &[f64]) -> f64 {
    let total: f64 = theta.iter().sum();
    let (mut ll, mut n) = (0.0, 0.0);
    for m in 0..fx.molecules() {
        let (lo, hi) = (fx.off[m], fx.off[m + 1]);
        let d: f64 = (lo..hi).map(|k| theta[fx.t[k]] * fx.w[k]).sum();
        if d > EM_DENOM_THRESH {
            ll += d.ln();
            n += 1.0;
        }
    }
    ll - n * total.ln()
}

/// State of the SQUAREM acceleration (Varadhan & Roland 2008, scheme S3) of the frozen-weight
/// EM: buffers for the iterates of a cycle, and the bound on the extrapolation step length,
/// which grows while long steps pay off and shrinks when one is rejected.
struct Squarem {
    t0: Vec<f64>,
    t1: Vec<f64>,
    t2: Vec<f64>,
    step_max: f64,
}

impl Squarem {
    const STEP_GROW: f64 = 4.0;

    fn new(n: usize) -> Self {
        Squarem {
            t0: vec![0.0; n],
            t1: vec![0.0; n],
            t2: vec![0.0; n],
            step_max: 1.0,
        }
    }
}

/// One SQUAREM cycle from the abundances in `prev`: two EM steps (θ1, θ2), a step along the
/// extrapolated path from θ0, and one EM step from there to stabilize it. The result is kept
/// only if its likelihood is at least that of θ2, the result of plain EM; otherwise θ2 is kept
/// and the step bound shrinks. On return `curr` holds the new abundances and `prev` the input
/// of the cycle's last EM step, so their difference is a one-step change as in plain EM.
fn squarem_cycle(fx: &Fixed, prev: &mut [f64], curr: &mut [f64], sq: &mut Squarem) {
    sq.t0.copy_from_slice(prev);
    step_static(fx, prev, curr);
    sq.t1.copy_from_slice(curr);
    prev.copy_from_slice(&sq.t1);
    step_static(fx, prev, curr);
    sq.t2.copy_from_slice(curr);
    let (mut rr, mut vv) = (0.0, 0.0);
    for i in 0..sq.t0.len() {
        let (a, b, c) = (sq.t0[i], sq.t1[i], sq.t2[i]);
        let (r, v) = (b - a, c - 2.0 * b + a);
        rr += r * r;
        vv += v * v;
    }
    // Step length -alpha in [1, step_max]; 1 reproduces θ2 exactly.
    let mut s = if vv > 0.0 { (rr / vv).sqrt() } else { 1.0 };
    if !s.is_finite() || s < 1.0 {
        s = 1.0;
    }
    let hit = s > 1.0 && s >= sq.step_max;
    s = s.min(sq.step_max);
    if s == 1.0 {
        // Plain EM: θ2 stands; report its last step. A step held back by the bound lets the
        // bound grow, as after an accepted long step.
        if hit {
            sq.step_max *= Squarem::STEP_GROW;
        }
        prev.copy_from_slice(&sq.t1);
        return;
    }
    for i in 0..prev.len() {
        let (a, b, c) = (sq.t0[i], sq.t1[i], sq.t2[i]);
        let (r, v) = (b - a, c - 2.0 * b + a);
        prev[i] = (a + 2.0 * s * r + s * s * v).max(0.0);
    }
    step_static(fx, prev, curr);
    let (l3, l2) = (loglik_static(fx, curr), loglik_static(fx, &sq.t2));
    if l3.is_finite() && l3 >= l2 {
        if hit {
            sq.step_max *= Squarem::STEP_GROW;
        }
    } else {
        sq.step_max = (sq.step_max / Squarem::STEP_GROW).max(1.0);
        prev.copy_from_slice(&sq.t1);
        curr.copy_from_slice(&sq.t2);
    }
}

/// Estimate the positional parameters once from the molecules of many cells: one set of
/// parameters shared by all cells, each cell with its own transcript abundances. The EM is the
/// parameter-estimating EM of [`fit_molecules`], each round's statistics summed over the cells,
/// and it stops when the parameters freeze (as oarfish's EM on a read sample does) or at
/// [`MAX_EM_ITER`]. The read-length laws and end kernels are tabulated up to `max_len` (or the
/// longest transcript in `cells`, if longer), which must cover every transcript of every cell
/// the parameters are then used in. When `cells` are a sample holding `share` of the data
/// (by reads), the log-likelihood freeze threshold shrinks by that share, as oarfish's does on
/// a read sample, so the sample freezes where a fit to all the cells would. The cells are split
/// over `n_threads` threads; the result depends on that number only through the order of
/// floating-point sums.
pub fn fit_global(
    cells: &[CellMolecules],
    max_len: u32,
    share: f64,
    n_threads: usize,
) -> GlobalFit {
    let max_len = cells
        .iter()
        .map(CellMolecules::max_len)
        .fold(max_len.max(1), u32::max);
    let mut par = Params::initial(max_len);
    let mut prev: Vec<Vec<f64>> = cells
        .iter()
        .map(|c| vec![c.molecule_count() as f64 / c.n_targets.max(1) as f64; c.n_targets])
        .collect();
    let mut curr: Vec<Vec<f64>> = cells.iter().map(|c| vec![0.0; c.n_targets]).collect();
    let per_thread = cells.len().div_ceil(n_threads.max(1)).max(1);
    let mut freeze = Freeze::new(share);
    let (mut frozen_at, mut freeze_reason) = (None, None);
    let mut last_ll = f64::NEG_INFINITY;
    let mut iterations = 0_usize;

    while iterations < MAX_EM_ITER as usize {
        let e = Exps::new(&par, max_len);
        let (par_ref, e_ref) = (&par, &e);
        let st = std::thread::scope(|scope| {
            let handles: Vec<_> = cells
                .chunks(per_thread)
                .zip(prev.chunks(per_thread))
                .zip(curr.chunks_mut(per_thread))
                .map(|((cs, ps), out)| {
                    scope.spawn(move || {
                        cs.iter().zip(ps).zip(out.iter_mut()).fold(
                            Stats::default(),
                            |st, ((cell, p), c)| st.merge(step_with(cell, par_ref, e_ref, p, c)),
                        )
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("truncation EM thread panicked"))
                .fold(Stats::default(), Stats::merge)
        });
        std::mem::swap(&mut prev, &mut curr);
        iterations += 1;
        let ll_gain = st.ll - last_ll;
        last_ll = st.ll;
        let next = st.params(&par);
        let reason = freeze.update(moved(&next, &par), ll_gain);
        par = next;
        if reason.is_some() {
            (frozen_at, freeze_reason) = (Some(iterations), reason);
            break;
        }
    }

    GlobalFit {
        params: TruncationParams::new(&par),
        par,
        cells: cells.len(),
        molecules: cells.iter().map(CellMolecules::molecule_count).sum(),
        iterations,
        frozen_at,
        freeze_reason,
    }
}

/// Re-estimate transcript probabilities for one cell's deduplicated molecules using the
/// truncation mixture model. The returned parameters are specific to this cell.
///
/// The fit follows oarfish's EM driver: the abundances and the positional parameters are
/// estimated together until the parameters freeze (see [`FREEZE_TOL`], [`FREEZE_LL_GAIN`]);
/// the abundances then continue under the fixed weights with SQUAREM acceleration until their
/// largest relative change is below [`CONVERGENCE_THRESH`] or [`MAX_EM_ITER`] is reached.
/// Abundances below [`MIN_READ_THRESH`] are then zeroed and each molecule is assigned its
/// posterior over its transcripts.
pub fn fit_molecules(
    molecules: &mut [ResolvedMolecule],
    interval_mode: TruncationIntervalMode,
) -> TruncationFit {
    fit_molecules_with(molecules, interval_mode, None)
}

/// As [`fit_molecules`], or, with `global`, with the positional parameters fixed to those
/// estimated from all cells ([`fit_global`]): the EM then estimates only the abundances, as
/// oarfish's does after estimating the parameters on a read sample ([`PLAIN_FIRST`] plain
/// steps, then SQUAREM), and the reported parameters are the global ones.
pub fn fit_molecules_with(
    molecules: &mut [ResolvedMolecule],
    interval_mode: TruncationIntervalMode,
    global: Option<&GlobalFit>,
) -> TruncationFit {
    let (cell, target_ids, mixed_interval_molecules) = build(molecules, interval_mode);
    let max_len = cell.max_len();
    let mut par = match global {
        Some(g) => g.par.clone(),
        None => Params::initial(max_len),
    };

    if molecules.is_empty() || target_ids.is_empty() {
        return TruncationFit {
            params: TruncationParams::new(&par),
            iterations: 0,
            frozen_at: None,
            freeze_reason: None,
            converged: false,
            molecule_count: molecules.len(),
            has_positional_evidence: false,
            mixed_interval_molecules,
        };
    }

    let n_targets = target_ids.len();
    let mut prev = vec![molecules.len() as f64 / n_targets as f64; n_targets];
    let mut curr = vec![0.0; n_targets];
    let mut sq = Squarem::new(n_targets);
    let mut frozen: Option<Fixed> =
        global.map(|_| Fixed::new(&cell, &par, &Exps::new(&par, max_len)));
    let (mut frozen_at, mut freeze_reason) = match global {
        Some(_) => (Some(0), Some(FreezeReason::Global)),
        None => (None, None),
    };
    let plain_first = if global.is_some() { PLAIN_FIRST } else { 0 };
    let mut freeze = Freeze::new(1.0);
    let mut niter = 0_u32;
    let mut steps = 0_usize;
    let mut last_ll = f64::NEG_INFINITY;
    let mut converged = false;
    let mut has_positional_evidence =
        global.is_some() && cell.alns.iter().any(|a| a.span.is_some());

    while niter < MAX_EM_ITER {
        match &frozen {
            Some(fx) if niter < plain_first => {
                step_static(fx, &prev, &mut curr);
                steps += 1;
            }
            Some(fx) => {
                squarem_cycle(fx, &mut prev, &mut curr, &mut sq);
                niter += 2;
                steps += 2;
            }
            None => {
                let e = Exps::new(&par, max_len);
                let st = step_with(&cell, &par, &e, &prev, &mut curr);
                steps += 1;
                has_positional_evidence |= st.n > 0.0;
                let ll_gain = st.ll - last_ll;
                last_ll = st.ll;
                let next = st.params(&par);
                let reason = freeze.update(moved(&next, &par), ll_gain);
                par = next;
                if reason.is_some() {
                    let e = Exps::new(&par, max_len);
                    frozen = Some(Fixed::new(&cell, &par, &e));
                    frozen_at = Some(steps);
                    freeze_reason = reason;
                }
            }
        }
        // Largest relative change in either direction (after a SQUAREM extrapolation most of
        // the motion is in decreases).
        let mut rel_diff = 0.0_f64;
        for (p, c) in prev.iter().zip(&curr) {
            if *p > MIN_READ_THRESH {
                rel_diff = rel_diff.max((c - p).abs() / p);
            }
        }
        std::mem::swap(&mut prev, &mut curr);
        if rel_diff < CONVERGENCE_THRESH && niter > 1 {
            converged = true;
            break;
        }
        niter += 1;
    }

    // Zero very small abundances, then assign each molecule its posterior with the
    // parameters fixed.
    for x in prev.iter_mut() {
        if *x < MIN_READ_THRESH {
            *x = 0.0;
        }
    }
    let exps = Exps::new(&par, max_len);
    for (m, molecule) in molecules.iter_mut().enumerate() {
        let alignments = cell.molecule(m);
        let mut weights: Vec<f64> = match &frozen {
            Some(fx) => (fx.off[m]..fx.off[m + 1])
                .map(|k| prev[fx.t[k]] * fx.w[k])
                .collect(),
            None => alignments
                .iter()
                .map(|a| a.weight(&prev, &par, &exps).0)
                .collect(),
        };
        let mut denominator: f64 = weights.iter().sum();
        if !(denominator.is_finite() && denominator > EM_DENOM_THRESH) {
            // The molecule is already deduplicated and must keep its counts, so it falls
            // back to its prior probabilities.
            weights = alignments.iter().map(|a| a.prior).collect();
            denominator = weights.iter().sum();
        }
        assert!(
            denominator.is_finite() && denominator > 0.0,
            "truncation model produced invalid final molecule probabilities"
        );
        let mut probabilities = vec![0.0; n_targets];
        for (alignment, weight) in alignments.iter().zip(weights) {
            probabilities[alignment.target] += weight / denominator;
        }
        molecule.transcript_probs = probabilities
            .into_iter()
            .enumerate()
            .filter_map(|(target, probability)| {
                (probability > 0.0).then_some((target_ids[target], probability))
            })
            .collect();
    }
    TruncationFit {
        params: TruncationParams::new(&par),
        iterations: steps,
        frozen_at,
        freeze_reason,
        converged,
        molecule_count: molecules.len(),
        has_positional_evidence,
        mixed_interval_molecules,
    }
}
