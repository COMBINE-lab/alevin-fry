//! Per-cell positional (truncation) model for long-read molecules.
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

use crate::eq_class::{merge_intervals, CellTxpCoverage};
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
    /// `relative`: whether the kernel has a relative part at all (without one, `rho` stays 0).
    fn initial(max_len: u32, relative: bool) -> Self {
        EndLaw {
            abs: GapLaw::initial(max_len),
            rel: vec![1.0 / RB as f64; RB],
            rho: if relative { 0.5 } else { 0.0 },
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
            rho: if old.rho == 0.0 {
                0.0
            } else {
                (s.n_rel / n).clamp(1e-6, 1.0 - 1e-6)
            },
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

/// Upper edges (nt) of the transcript-length classes whose read-length laws are estimated
/// separately under [`TruncationLengthLaw::Conditional`]. How long a read is depends on the
/// transcript it comes from: on short transcripts reads cover a fraction of the transcript,
/// on long ones they have a length of their own. With a single law, reads that stop short of
/// a short transcript's 5' end can only be explained as censored full-length reads with a lost
/// 5' end, which widens the 5' kernel and favors shorter isoforms for every read.
const LEN_CLASS_EDGES: [u32; 13] = [
    300, 500, 700, 900, 1100, 1400, 1700, 2000, 2500, 3000, 4000, 5500, 8000,
];
const N_LEN_CLASSES: usize = LEN_CLASS_EDGES.len() + 1;

/// How a read's length law depends on the length of its transcript.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub enum TruncationLengthLaw {
    /// One law per component for every transcript (the oarfish model).
    #[default]
    Marginal,
    /// One law per component and transcript-length class ([`LEN_CLASS_EDGES`]).
    Conditional,
}

impl fmt::Display for TruncationLengthLaw {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Marginal => "marginal",
                Self::Conditional => "conditional",
            }
        )
    }
}

impl FromStr for TruncationLengthLaw {
    type Err = &'static str;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "marginal" => Ok(Self::Marginal),
            "conditional" => Ok(Self::Conditional),
            _ => Err("expected marginal or conditional"),
        }
    }
}

/// Whether the end kernels have a relative part (a gap proportional to the transcript's length).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub enum TruncationEndKernels {
    /// Absolute and relative parts (the oarfish model).
    #[default]
    Mixed,
    /// Absolute part only. The relative part gives a short transcript a density of order
    /// `rho / len` at every gap; at the 5' end of 3'-anchored reads it turns truncated reads
    /// into censored full-length ones and favors shorter isoforms.
    Absolute,
}

impl fmt::Display for TruncationEndKernels {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Mixed => "mixed",
                Self::Absolute => "absolute",
            }
        )
    }
}

impl FromStr for TruncationEndKernels {
    type Err = &'static str;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "mixed" => Ok(Self::Mixed),
            "absolute" => Ok(Self::Absolute),
            _ => Err("expected mixed or absolute"),
        }
    }
}

/// Which positional law the model uses.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub enum TruncationPositionalLaw {
    /// oarfish's four-component mixture (3'-anchored, 5'-anchored, full length, uniformly placed).
    #[default]
    Mixture,
    /// The earlier four-parameter model ([`FourParams`]), run with this module's EM, freezing,
    /// pooling and shrinkage instead of its own per-cell EM.
    FourParam,
}

impl fmt::Display for TruncationPositionalLaw {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Mixture => "mixture",
                Self::FourParam => "four-param",
            }
        )
    }
}

impl FromStr for TruncationPositionalLaw {
    type Err = &'static str;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "mixture" => Ok(Self::Mixture),
            "four-param" | "four_param" | "fourparam" => Ok(Self::FourParam),
            _ => Err("expected mixture or four-param"),
        }
    }
}

/// Where the coverage models (`--coverage-model`, `--coverage-dedup`) take their coverage from.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub enum CoverageSource {
    /// Each cell's own molecules, as before.
    #[default]
    Cell,
    /// The molecules of all cells, accumulated in a pass before quantifying: transcript-specific
    /// coverage that a handful of molecules per cell cannot show.
    Global,
}

impl fmt::Display for CoverageSource {
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

impl FromStr for CoverageSource {
    type Err = &'static str;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "cell" => Ok(Self::Cell),
            "global" => Ok(Self::Global),
            _ => Err("expected cell or global"),
        }
    }
}

/// The shape of the positional model.
#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub struct TruncationModelSpec {
    pub length_law: TruncationLengthLaw,
    pub end_kernels: TruncationEndKernels,
    /// Power the positional likelihood is raised to in each alignment's weight: 1 is the
    /// model as it stands, 0 ignores positions, values between temper it.
    pub positional_power: f64,
    /// With pooled abundances: choose the power from the data ([`calibrate_power`]) instead of
    /// using `positional_power`.
    pub calibrate_power: bool,
    /// A molecule's final shares below this are dropped and the rest renormalized
    /// ([`prune_shares`]); 0 keeps every share.
    pub min_share: f64,
    /// A molecule's final shares below this give their combined mass to one of them, drawn by
    /// share ([`sample_small_shares`]); 1 draws one transcript per molecule, 0 keeps every share.
    pub sample_below: f64,
    /// With [`TruncationAbundance::Shrink`]: estimate the strength of the pooled prior from the
    /// data ([`prior_strength`]) instead of taking it from `--truncation-shrink K`.
    pub shrink_auto: bool,
    /// With `shrink_auto`: one strength per locus (gene, with [`GeneMap`]), shrunk toward the
    /// common one ([`GENE_PRIOR_WEIGHT`]).
    pub shrink_per_gene: bool,
    /// Which positional law the model uses.
    pub positional_law: TruncationPositionalLaw,
    /// With `--truncation-groups`: pseudo-molecules per locus pulling each group's pooled
    /// isoform shares toward the shares of all cells ([`fit_pooled`]); 0 pools the groups
    /// independently.
    pub group_shrink: f64,
    /// Estimate `group_shrink` from how much isoform usage varies between the groups
    /// ([`prior_strength_counts`] over groups) instead of taking it as given.
    pub group_shrink_auto: bool,
    /// With `group_shrink_auto`: the loci are genes ([`GeneMap`]); without it they are the
    /// connected sets of candidate transcripts.
    pub group_shrink_per_gene: bool,
}

impl Default for TruncationModelSpec {
    fn default() -> Self {
        Self {
            length_law: TruncationLengthLaw::default(),
            end_kernels: TruncationEndKernels::default(),
            positional_power: 1.0,
            calibrate_power: false,
            min_share: 0.0,
            sample_below: 0.0,
            shrink_auto: false,
            shrink_per_gene: false,
            positional_law: TruncationPositionalLaw::default(),
            group_shrink: 0.0,
            group_shrink_auto: false,
            group_shrink_per_gene: false,
        }
    }
}

/// A read-length law's density and tail sums by bin, for constant-time density and survival
/// lookups (the same values as [`LenLaw::tables`], without a table per base).
struct LenTab {
    h: Vec<f64>,
    tail: Vec<f64>,
}

impl LenTab {
    fn new(law: &LenLaw) -> Self {
        let nb = law.h.len();
        let mut tail = vec![0.0; nb + 1];
        for b in (0..nb).rev() {
            tail[b] = tail[b + 1] + law.h[b];
        }
        LenTab {
            h: law.h.clone(),
            tail,
        }
    }

    #[inline]
    fn density(&self, x: u32) -> f64 {
        self.h.get(LenLaw::bin_nt(x)).map_or(0.0, |p| p / LenLaw::W)
    }

    /// `P(length >= x)`.
    #[inline]
    fn survival(&self, x: u32) -> f64 {
        let b = LenLaw::bin_nt(x);
        if b >= self.h.len() {
            return 0.0;
        }
        let within = (x as f64 - b as f64 * LenLaw::W) / LenLaw::W;
        (self.tail[b + 1] + self.h[b] * (1.0 - within)).max(0.0)
    }
}

/// The window (nt) at a transcript's end within which the four-parameter law counts a read end
/// as at the end, as in the earlier model.
const FOUR_END_WINDOW: f64 = 50.0;

/// The earlier four-parameter law ([`TruncationPositionalLaw::FourParam`]): a read reaches the
/// transcript's 5' end (within [`FOUR_END_WINDOW`]) with probability `full`, otherwise it stopped
/// at a geometric distance from the 3' end with rate `stop` per base; its 3' end lies within the
/// window of the transcript's 3' end with probability `three_prime`, otherwise at an exponential
/// distance with rate `decay`. The two ends are independent.
#[derive(Clone, Copy, Debug)]
struct FourParams {
    full: f64,
    stop: f64,
    three_prime: f64,
    decay: f64,
}

/// The four-parameter law's terms for one alignment: the ways its 5' end (reached the end and
/// kept it, reached it as a survivor of the stopping process, or stopped) and its 3' end (at the
/// end, or off it) can be explained, each already weighted by its probability.
#[derive(Clone, Copy)]
struct FourTerms {
    full_w: f64,
    full_surv: f64,
    stopped: f64,
    at3: f64,
    off3: f64,
    ext: f64,
    d3: f64,
    len: f64,
}

impl FourTerms {
    fn start(&self) -> f64 {
        self.full_w + self.full_surv + self.stopped
    }

    fn end(&self) -> f64 {
        self.at3 + self.off3
    }

    fn total(&self) -> f64 {
        self.start() * self.end()
    }
}

impl FourParams {
    fn initial() -> Self {
        FourParams {
            full: 0.5,
            stop: 1e-3,
            three_prime: 0.5,
            decay: 2e-3,
        }
    }

    fn terms(&self, g: &Geo) -> FourTerms {
        let len = g.len as f64;
        let ext = len - g.d5 as f64;
        let d3 = g.d3 as f64;
        let full = ext > len - FOUR_END_WINDOW;
        FourTerms {
            full_w: if full { self.full / FOUR_END_WINDOW } else { 0.0 },
            full_surv: if full {
                (1.0 - self.full) * (-self.stop * len).exp() / FOUR_END_WINDOW
            } else {
                0.0
            },
            stopped: (1.0 - self.full) * self.stop * (-self.stop * ext).exp(),
            at3: if d3 < FOUR_END_WINDOW { self.three_prime / FOUR_END_WINDOW } else { 0.0 },
            off3: (1.0 - self.three_prime) * self.decay * (-self.decay * d3).exp(),
            ext,
            d3,
            len,
        }
    }

    /// The M-step of the earlier model.
    fn updated(&self, s: &FourStats) -> Self {
        if s.n <= 0.0 {
            return *self;
        }
        let fraction = |v: f64| v.clamp(1e-6, 1.0 - 1e-6);
        FourParams {
            full: fraction(s.n_full_w / s.n),
            stop: if s.exposure_stop > 0.0 {
                (s.n_stopped / s.exposure_stop).clamp(1e-7, 1.0)
            } else {
                self.stop
            },
            three_prime: fraction(s.n_at3 / s.n),
            decay: if s.exposure_three_prime > 0.0 {
                (s.n_off3 / s.exposure_three_prime).clamp(1e-7, 1.0)
            } else {
                self.decay
            },
        }
    }

    fn moved(&self, o: &Self) -> f64 {
        let relative = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-12);
        (self.full - o.full)
            .abs()
            .max((self.three_prime - o.three_prime).abs())
            .max(relative(self.stop, o.stop))
            .max(relative(self.decay, o.decay))
    }
}

/// Expected sufficient statistics of the four-parameter law.
#[derive(Clone, Copy, Default)]
struct FourStats {
    n: f64,
    n_full_w: f64,
    n_stopped: f64,
    n_at3: f64,
    n_off3: f64,
    exposure_stop: f64,
    exposure_three_prime: f64,
}

impl FourStats {
    fn add(&mut self, t: &FourTerms, post: f64) {
        let (start, end) = (t.start(), t.end());
        // an alignment whose terms underflowed carries no information about the parameters
        if !(start > 0.0 && end > 0.0) {
            return;
        }
        let p_full_w = post * t.full_w / start;
        let p_survived = post * t.full_surv / start;
        let p_stopped = post * t.stopped / start;
        let p_at3 = post * t.at3 / end;
        let p_off3 = post * t.off3 / end;
        self.n += post;
        self.n_full_w += p_full_w;
        self.n_stopped += p_stopped;
        self.n_at3 += p_at3;
        self.n_off3 += p_off3;
        self.exposure_stop += p_survived * t.len + p_stopped * t.ext;
        self.exposure_three_prime += p_off3 * t.d3;
    }

    fn merge(mut self, o: FourStats) -> FourStats {
        self.n += o.n;
        self.n_full_w += o.n_full_w;
        self.n_stopped += o.n_stopped;
        self.n_at3 += o.n_at3;
        self.n_off3 += o.n_off3;
        self.exposure_stop += o.exposure_stop;
        self.exposure_three_prime += o.exposure_three_prime;
        self
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
    /// Read-length laws of the 3'-anchored, 5'-anchored and uniform components, one per
    /// transcript-length class (a single one under [`TruncationLengthLaw::Marginal`]).
    len3: Vec<LenLaw>,
    len5: Vec<LenLaw>,
    lenu: Vec<LenLaw>,
    /// End kernels at the 3' and 5' ends.
    k3: EndLaw,
    k5: EndLaw,
    /// See [`TruncationModelSpec::positional_power`].
    pos_power: f64,
    /// The four-parameter law's parameters, under [`TruncationPositionalLaw::FourParam`].
    four: Option<FourParams>,
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
    /// The four-parameter law's terms instead, under that law.
    four: Option<FourTerms>,
}

impl Terms {
    fn four(f: FourTerms) -> Terms {
        Terms {
            a3_stop: 0.0,
            a3_cens: 0.0,
            a5_stop: 0.0,
            a5_cens: 0.0,
            full: 0.0,
            u_stop: 0.0,
            u_cens: 0.0,
            k3a: 0.0,
            k3r: 0.0,
            k5a: 0.0,
            k5r: 0.0,
            four: Some(f),
        }
    }

    #[inline]
    fn total(&self) -> f64 {
        if let Some(f) = &self.four {
            return f.total();
        }
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
    fn initial(max_len: u32, spec: TruncationModelSpec) -> Self {
        let relative = spec.end_kernels == TruncationEndKernels::Mixed;
        let classes = match spec.length_law {
            TruncationLengthLaw::Marginal => 1,
            TruncationLengthLaw::Conditional => N_LEN_CLASSES,
        };
        let laws = || vec![LenLaw::exponential(1000.0, max_len); classes];
        Params {
            p3: 0.25,
            p5: 0.25,
            pf: 0.25,
            pu: 0.25,
            len3: laws(),
            len5: laws(),
            lenu: laws(),
            k3: EndLaw::initial(max_len, relative),
            k5: EndLaw::initial(max_len, relative),
            pos_power: spec.positional_power,
            four: (spec.positional_law == TruncationPositionalLaw::FourParam)
                .then(FourParams::initial),
        }
    }

    /// The positional likelihood `p` as it enters an alignment's weight.
    #[inline]
    fn positional(&self, p: f64) -> f64 {
        if self.pos_power == 1.0 {
            p
        } else {
            p.powf(self.pos_power)
        }
    }

    /// The read-length class of a transcript of length `len`.
    #[inline]
    fn class_of(&self, len: u32) -> usize {
        if self.len3.len() == 1 {
            0
        } else {
            LEN_CLASS_EDGES.partition_point(|&e| e <= len)
        }
    }

    #[inline]
    fn terms(&self, g: &Geo, e: &Exps) -> Terms {
        if let Some(fp) = self.four {
            return Terms::four(fp.terms(g));
        }
        let (len, l, d5, d3) = (g.len, g.l, g.d5, g.d3);
        let c = self.class_of(len);
        let (l3, l5, lu) = (&e.l3[c], &e.l5[c], &e.lu[c]);
        // end kernels, each a distribution over the gaps this transcript allows (0..len)
        let (k3a, k3r) = self.k3.parts(d3, len, &e.g3d, &e.g3c);
        let (k5a, k5r) = self.k5.parts(d5, len, &e.g5d, &e.g5c);
        let (kk3, kk5) = (k3a + k3r, k5a + k5r);
        Terms {
            a3_stop: self.p3 * kk3 * l3.density(l),
            a3_cens: self.p3 * kk3 * l3.survival(len - d3) * kk5,
            a5_stop: self.p5 * kk5 * l5.density(l),
            a5_cens: self.p5 * kk5 * l5.survival(len - d5) * kk3,
            full: self.pf * kk3 * kk5,
            u_stop: self.pu * lu.density(l) / ((len - l) as f64 + 1.0),
            u_cens: self.pu * lu.survival(len) * kk3 * kk5,
            k3a,
            k3r,
            k5a,
            k5r,
            four: None,
        }
    }
}

/// The laws in the form the likelihood reads them, built once per iteration: the length laws
/// by bin (one per length class), the end kernels' absolute parts for every whole number of
/// bases up to the longest transcript.
struct Exps {
    /// Length laws of the three components, per length class.
    l3: Vec<LenTab>,
    l5: Vec<LenTab>,
    lu: Vec<LenTab>,
    /// End-kernel densities and distribution functions.
    g3d: Vec<f64>,
    g3c: Vec<f64>,
    g5d: Vec<f64>,
    g5c: Vec<f64>,
}

impl Exps {
    fn new(par: &Params, max_len: u32) -> Self {
        let tabs = |laws: &[LenLaw]| laws.iter().map(LenTab::new).collect();
        let (g3d, g3c) = par.k3.abs.tables(max_len);
        let (g5d, g5c) = par.k5.abs.tables(max_len);
        Exps {
            l3: tabs(&par.len3),
            l5: tabs(&par.len5),
            lu: tabs(&par.lenu),
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
    /// Length laws, per length class.
    l3: Vec<LStats>,
    l5: Vec<LStats>,
    lu: Vec<LStats>,
    /// End kernels: observed gaps (`stop`) and truncation points (`cens`), binned.
    k3: EStats,
    k5: EStats,
    /// Log-likelihood of the reads at the step's input (abundances unnormalized).
    ll: f64,
    /// The four-parameter law's statistics, under that law.
    four: FourStats,
}

impl Stats {
    fn merge(mut self, o: Stats) -> Stats {
        self.ll += o.ll;
        self.n += o.n;
        self.n3 += o.n3;
        self.n5 += o.n5;
        self.nf += o.nf;
        self.nu += o.nu;
        let merge = |a: Vec<LStats>, b: Vec<LStats>| -> Vec<LStats> {
            if a.is_empty() {
                return b;
            }
            if b.is_empty() {
                return a;
            }
            a.into_iter().zip(b).map(|(x, y)| x.merge(y)).collect()
        };
        self.l3 = merge(std::mem::take(&mut self.l3), o.l3);
        self.l5 = merge(std::mem::take(&mut self.l5), o.l5);
        self.lu = merge(std::mem::take(&mut self.lu), o.lu);
        self.k3 = self.k3.merge(o.k3);
        self.k5 = self.k5.merge(o.k5);
        self.four = self.four.merge(o.four);
        self
    }

    fn add(&mut self, g: &Geo, t: &Terms, post: f64, class: usize, n_classes: usize) {
        let tot = t.total();
        if tot <= 0.0 {
            return;
        }
        if let Some(f) = &t.four {
            self.n += post;
            self.four.add(f, post);
            return;
        }
        if self.l3.is_empty() {
            self.l3 = vec![LStats::default(); n_classes];
            self.l5 = vec![LStats::default(); n_classes];
            self.lu = vec![LStats::default(); n_classes];
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
        let (l3, l5, lu) = (&mut self.l3[class], &mut self.l5[class], &mut self.lu[class]);
        LStats::add(&mut l3.stop, LenLaw::bin_nt(g.l), a3s);
        LStats::add(&mut l3.cens, LenLaw::bin_nt(g.len - g.d3), a3c);
        LStats::add(&mut l5.stop, LenLaw::bin_nt(g.l), a5s);
        LStats::add(&mut l5.cens, LenLaw::bin_nt(g.len - g.d5), a5c);
        LStats::add(&mut lu.stop, LenLaw::bin_nt(g.l), us);
        LStats::add(&mut lu.cens, LenLaw::bin_nt(g.len), uc);
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
        if let Some(fp) = old.four {
            return Params {
                four: Some(fp.updated(&self.four)),
                ..old.clone()
            };
        }
        // a class without statistics keeps its law
        let laws = |st: &[LStats], old: &[LenLaw]| -> Vec<LenLaw> {
            old.iter()
                .enumerate()
                .map(|(c, law)| st.get(c).map_or_else(|| law.clone(), |s| s.law(law)))
                .collect()
        };
        let frac = |x: f64| (x / self.n).clamp(1e-6, 1.0);
        let (p3, p5, pf, pu) = (frac(self.n3), frac(self.n5), frac(self.nf), frac(self.nu));
        let z = p3 + p5 + pf + pu;
        Params {
            p3: p3 / z,
            p5: p5 / z,
            pf: pf / z,
            pu: pu / z,
            len3: laws(&self.l3, &old.len3),
            len5: laws(&self.l5, &old.len5),
            lenu: laws(&self.lu, &old.lenu),
            k3: EndLaw::law(&self.k3, &old.k3),
            k5: EndLaw::law(&self.k5, &old.k5),
            pos_power: old.pos_power,
            four: old.four,
        }
    }
}

/// How far the parameters moved: absolute change for the mixture weights; for each length law,
/// its change scaled by the weight of the one component that uses it, as that bounds its effect
/// on the likelihood (the law of a component holding a negligible share of the reads drifts for
/// as long as the EM runs, and must not hold the parameters open). The end kernels enter
/// several components and count in full.
fn moved(a: &Params, b: &Params) -> f64 {
    if let (Some(x), Some(y)) = (a.four, b.four) {
        return x.moved(&y);
    }
    (a.p3 - b.p3)
        .abs()
        .max((a.p5 - b.p5).abs())
        .max((a.pf - b.pf).abs())
        .max((a.pu - b.pu).abs())
        .max(a.p3.max(b.p3) * laws_moved(&a.len3, &b.len3))
        .max(a.p5.max(b.p5) * laws_moved(&a.len5, &b.len5))
        .max(a.pu.max(b.pu) * laws_moved(&a.lenu, &b.lenu))
        .max(a.k3.moved(&b.k3))
        .max(a.k5.moved(&b.k5))
}

fn laws_moved(a: &[LenLaw], b: &[LenLaw]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x.moved(y)).fold(0.0, f64::max)
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
    /// Intervals are ignored: no positional likelihood, only the alignment-score prior and the
    /// abundances (a control for what the positional model adds).
    Ignore,
}

impl fmt::Display for TruncationIntervalMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::All => "all",
                Self::Union => "union",
                Self::Ignore => "none",
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
            "none" => Ok(Self::Ignore),
            _ => Err("expected all, union or none"),
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
    /// Under the four-parameter law: full-length fraction, 5' stop rate, 3'-anchored fraction, 3' decay.
    pub four: Option<[f64; 4]>,
}

/// The mean read length over the length classes' laws (their plain average: a summary only).
fn mean_of(laws: &[LenLaw]) -> f64 {
    laws.iter().map(LenLaw::mean).sum::<f64>() / laws.len().max(1) as f64
}

impl TruncationParams {
    fn new(p: &Params) -> Self {
        Self {
            p3: p.p3,
            p5: p.p5,
            pf: p.pf,
            pu: p.pu,
            len3_mean: mean_of(&p.len3),
            len5_mean: mean_of(&p.len5),
            lenu_mean: mean_of(&p.lenu),
            k3_rho: p.k3.rho,
            k5_rho: p.k5.rho,
            k3_gap: p.k3.abs.summary(),
            k5_gap: p.k5.abs.summary(),
            four: p.four.map(|f| [f.full, f.stop, f.three_prime, f.decay]),
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

/// Where each cell's transcript abundances come from.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub enum TruncationAbundance {
    /// Estimated by each cell's own EM (maximum likelihood from that cell's molecules alone).
    #[default]
    Cell,
    /// One abundance per transcript estimated from the molecules of all cells; each cell's
    /// molecules are assigned by their posteriors under it. Isoform usage is then shared by all
    /// cells, which per-cell estimates cannot resolve from a handful of molecules per gene.
    Global,
    /// Each cell's own EM with a Dirichlet prior centered on the pooled isoform fractions of
    /// each locus ([`fit_molecules_with`]), so cells keep their own isoform usage where their
    /// molecules show it.
    Shrink,
}

impl fmt::Display for TruncationAbundance {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Cell => "cell",
                Self::Global => "global",
                Self::Shrink => "shrink",
            }
        )
    }
}

impl FromStr for TruncationAbundance {
    type Err = &'static str;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "cell" => Ok(Self::Cell),
            "global" => Ok(Self::Global),
            "shrink" => Ok(Self::Shrink),
            _ => Err("expected cell, global or shrink"),
        }
    }
}

/// The strength of the pooled prior under [`TruncationAbundance::Shrink`] (`--truncation-shrink`).
#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub enum TruncationShrink {
    /// `k` pseudo-molecules per molecule the cell has in the locus ([`assign_pooled`]).
    Fixed(f64),
    /// A number of pseudo-molecules per locus estimated from the data ([`prior_strength`]).
    Auto,
    /// As `Auto`, but one number per gene, shrunk toward the common one: a gene whose isoform
    /// shares vary between cells more than sampling allows gets a weak prior, the rest a strong one.
    AutoGene,
}

impl fmt::Display for TruncationShrink {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Fixed(k) => write!(f, "{k}"),
            Self::Auto => write!(f, "auto"),
            Self::AutoGene => write!(f, "auto-gene"),
        }
    }
}

impl FromStr for TruncationShrink {
    type Err = &'static str;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value.eq_ignore_ascii_case("auto") {
            return Ok(Self::Auto);
        }
        if ["auto-gene", "autogene", "auto_gene"].iter().any(|v| value.eq_ignore_ascii_case(v)) {
            return Ok(Self::AutoGene);
        }
        match value.parse::<f64>() {
            Ok(k) if k.is_finite() && k >= 0.0 => Ok(Self::Fixed(k)),
            _ => Err("expected a finite, non-negative number, auto or auto-gene"),
        }
    }
}

/// A named configuration of the truncation model's options (`--truncation-preset`); an option
/// given explicitly keeps its value.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum TruncationPreset {
    /// The configuration that did best across simulated and experimental long-read single-cell
    /// data: positional parameters pooled over all cells, absolute end kernels, the positional
    /// likelihood at power 0.5, each cell's abundances shrunk toward its group's pooled isoform
    /// shares with one pseudo-molecule per molecule (`--truncation-shrink 1`), the groups found
    /// from the data (`--truncation-discover-groups`) unless given, and the groups' own shares
    /// kept only in the loci where they differ beyond sampling (`--truncation-group-shrink auto`).
    Recommended,
}

impl fmt::Display for TruncationPreset {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Recommended => write!(f, "recommended"),
        }
    }
}

impl FromStr for TruncationPreset {
    type Err = &'static str;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "recommended" => Ok(Self::Recommended),
            _ => Err("expected recommended"),
        }
    }
}

/// Groups of cells that share pooled abundances (`--truncation-groups`): a barcode -> group map;
/// barcodes it does not list form one more group of their own.
#[derive(Clone, Debug, Default)]
pub struct TruncationGroups {
    map: HashMap<String, u32>,
    names: Vec<String>,
}

impl TruncationGroups {
    /// Reads `barcode<TAB>group` lines (no header).
    pub fn load(path: &std::path::Path) -> std::io::Result<Self> {
        let mut g = TruncationGroups::default();
        let mut index = HashMap::new();
        for line in std::fs::read_to_string(path)?.lines() {
            let mut f = line.split('\t');
            let (Some(bc), Some(name)) = (f.next(), f.next()) else {
                continue;
            };
            let next = g.names.len() as u32;
            let gi = *index.entry(name.trim().to_string()).or_insert_with(|| {
                g.names.push(name.trim().to_string());
                next
            });
            g.map.insert(bc.trim().to_string(), gi);
        }
        Ok(g)
    }

    /// The group of a barcode: its listed group, else the group of the unlisted barcodes.
    pub fn group_of(&self, barcode: &str) -> u32 {
        self.map.get(barcode).copied().unwrap_or(self.names.len() as u32)
    }

    /// Number of groups, the unlisted barcodes' group included.
    pub fn n_groups(&self) -> usize {
        self.names.len() + 1
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn n_barcodes(&self) -> usize {
        self.map.len()
    }

    /// The group of `barcode`, adding a group named after it when the barcode is not listed yet.
    pub fn own_group(&mut self, barcode: &str) -> u32 {
        if let Some(&g) = self.map.get(barcode) {
            return g;
        }
        let g = self.names.len() as u32;
        self.names.push(barcode.to_string());
        self.map.insert(barcode.to_string(), g);
        g
    }

    /// Groups from discovered labels (`--truncation-discover-groups`): group `k` is named `L{k}`.
    pub fn from_labels(barcodes: &[String], labels: &[u32]) -> Self {
        let n = labels.iter().map(|&l| l as usize + 1).max().unwrap_or(0);
        TruncationGroups {
            map: barcodes.iter().cloned().zip(labels.iter().copied()).collect(),
            names: (0..n).map(|k| format!("L{k}")).collect(),
        }
    }
}

/// Transcript-to-gene assignment for the pooled prior (`--truncation-genes`): with it the prior
/// acts within genes instead of within loci, which merge genes through molecules whose candidates
/// span several (a locus can hold thousands of transcripts), so that differences in gene
/// expression between cells are not mistaken for differences in isoform usage.
#[derive(Clone, Debug)]
pub struct GeneMap {
    /// Each transcript's gene as an index; a transcript not in the file gets an index of its own.
    pub gene_of: Vec<u32>,
    /// Genes named in the file.
    pub genes: usize,
    /// Transcripts of the reference found in the file.
    pub matched: usize,
}

impl GeneMap {
    /// Reads `transcript<TAB>gene` lines (no header; further columns are ignored). A transcript is
    /// matched by its name as written, else by its name without a version suffix (`.N`).
    pub fn load(path: &std::path::Path, ref_names: &[String]) -> std::io::Result<Self> {
        fn strip(s: &str) -> &str {
            match s.rsplit_once('.') {
                Some((stem, v)) if !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()) => stem,
                _ => s,
            }
        }
        let mut gene_index: HashMap<String, u32> = HashMap::new();
        let mut exact: HashMap<String, u32> = HashMap::new();
        let mut stripped: HashMap<String, u32> = HashMap::new();
        for line in std::fs::read_to_string(path)?.lines() {
            let mut f = line.split('\t');
            let (Some(t), Some(g)) = (f.next(), f.next()) else {
                continue;
            };
            let next = gene_index.len() as u32;
            let gi = *gene_index.entry(g.trim().to_string()).or_insert(next);
            exact.entry(t.trim().to_string()).or_insert(gi);
            stripped.entry(strip(t.trim()).to_string()).or_insert(gi);
        }
        let mut next_own = gene_index.len() as u32;
        let mut matched = 0;
        let gene_of = ref_names
            .iter()
            .map(|n| match exact.get(n).or_else(|| stripped.get(strip(n))) {
                Some(&g) => {
                    matched += 1;
                    g
                }
                None => {
                    next_own += 1;
                    next_own - 1
                }
            })
            .collect();
        Ok(GeneMap {
            gene_of,
            genes: gene_index.len(),
            matched,
        })
    }
}

/// The strength of the pooled prior estimated from the data ([`prior_strength`]).
#[derive(Clone, Debug)]
pub struct PriorStrength {
    /// Pseudo-molecules per locus; infinite when the cells show no variation in isoform usage
    /// beyond sampling, so that every cell uses the pooled shares as they are.
    pub alpha: f64,
    /// Loci with at least two transcripts that have one-candidate molecules.
    pub loci: usize,
    /// Cells with at least two one-candidate molecules in some such locus.
    pub cells: usize,
    /// One-candidate molecules in those loci.
    pub molecules: usize,
    /// Whether the loci were split by gene ([`GeneMap`]).
    pub within_genes: bool,
    /// Under `auto-gene`: the strength of each locus, indexed by its label (infinite where the
    /// locus shows no variation beyond sampling or has no data).
    pub per_locus: Option<Vec<f64>>,
    /// Under `auto-gene`: loci with a finite strength, and the median of those strengths.
    pub loci_finite: usize,
    pub alpha_median: f64,
}

impl PriorStrength {
    /// The weight a cell's own `n` molecules of a locus get against the pooled shares.
    pub fn own_weight(&self, n: f64) -> f64 {
        if self.alpha.is_infinite() {
            0.0
        } else {
            n / (n + self.alpha)
        }
    }
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
    /// With pooled abundances: each transcript's abundance (molecules) estimated from the
    /// molecules of all the cells of a group, per group, indexed by transcript id.
    alpha: Option<Vec<Vec<f64>>>,
    /// With pooled abundances: each transcript's locus (transcripts are in one locus when a
    /// molecule has both as candidates), indexed by transcript id.
    locus: Option<Vec<u32>>,
    /// EM steps of the pooled-abundance fit after the parameters froze.
    pub abundance_iterations: usize,
    /// The calibrated positional power ([`calibrate_power`]), when it was estimated.
    pub calibration: Option<Calibration>,
    /// The strength of the pooled prior ([`prior_strength`]), when it was estimated
    /// (`--truncation-shrink auto`).
    pub prior: Option<PriorStrength>,
    /// The strength of the group-level prior (`--truncation-group-shrink auto`), when it was
    /// estimated: how much isoform usage varies between the groups.
    pub group_prior: Option<PriorStrength>,
    /// The junction law ([`JunctionLaw`]) learned from the cells, when asked for.
    pub junction: Option<std::sync::Arc<JunctionLaw>>,
    locus_mass_cache: std::sync::OnceLock<Vec<Vec<f64>>>,
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
    fn fixed_weight(&self, par: &Params, e: &Exps, positional: bool) -> f64 {
        match self.geo().filter(|_| positional) {
            Some(g) => self.prior * par.positional(par.terms(&g, e).total()),
            None => self.prior,
        }
    }

    /// The weight of this alignment under abundances `ab` and the current parameters, with the
    /// terms that produced it.
    fn weight(
        &self,
        ab: &[f64],
        par: &Params,
        e: &Exps,
        positional: bool,
    ) -> (f64, Option<(Geo, Terms)>) {
        self.weight_at(ab[self.target], par, e, positional)
    }

    /// As [`Alignment::weight`], for its target's abundance `ab`.
    fn weight_at(
        &self,
        ab: f64,
        par: &Params,
        e: &Exps,
        positional: bool,
    ) -> (f64, Option<(Geo, Terms)>) {
        match self.geo().filter(|_| positional) {
            Some(g) => {
                let t = par.terms(&g, e);
                (ab * self.prior * par.positional(t.total()), Some((g, t)))
            }
            None => (ab * self.prior, None),
        }
    }
}

/// Per-alignment factors multiplied into the alignment-score prior before the positional model.
#[derive(Clone, Copy, Default)]
pub struct PriorFactors<'a> {
    /// A coverage model, this cell's own or one learned from all cells (`--coverage-source`):
    /// each candidate's prior is multiplied by the coverage probability of the molecule's
    /// intervals on it, as `--coverage-model` does before its EM.
    pub coverage: Option<&'a [CellTxpCoverage]>,
    /// A junction law (`--junction-likelihood`): each interval's prior is multiplied by the law's
    /// probability of its start.
    pub junction: Option<&'a JunctionLaw>,
}

/// Inputs of [`fit_global`] beyond the cells.
#[derive(Clone, Copy, Default)]
pub struct FitExtras<'a> {
    /// Each transcript's gene ([`GeneMap`]), to split the loci the pooled prior acts within.
    pub genes: Option<&'a [u32]>,
    /// Each transcript's splice-junction positions (ascending), to learn a [`JunctionLaw`].
    pub junctions: Option<&'a [Vec<u32>]>,
}

/// The coverage weight of a molecule's intervals on `target`, as `--coverage-model` computes it.
fn coverage_weight(molecule: &ResolvedMolecule, target: u32, coverage: &[CellTxpCoverage]) -> f64 {
    molecule.intervals.get(&target).map_or(1.0, |intervals| {
        let merged = merge_intervals(intervals.iter().map(|i| (i.start, i.end)).collect());
        coverage
            .get(target as usize)
            .map_or(1.0, |c| c.coverage_probability_for_intervals(&merged))
    })
}

/// A law for where a read starts on a transcript, relative to the splice junction before it
/// (`--junction-likelihood`): the probability that a start lies `d` bases after the previous
/// junction (the transcript's start counting as one), learned from the one-candidate molecules
/// of all cells on log-spaced bins of `d` ([`GapLaw::bin`]) and normalized over every start
/// position of each transcript, so that it compares like a likelihood across isoforms of
/// different structure and length. The earlier `--junction-model` weight is not normalized and
/// favours the isoform with the shortest extension before the read.
#[derive(Clone, Debug)]
pub struct JunctionLaw {
    junctions: Vec<Vec<u32>>,
    /// Density per base by bin of the distance.
    dens: Vec<f64>,
    /// `cum[d]`: the density summed over the distances below `d`.
    cum: Vec<f64>,
    /// Per transcript, the density summed over all its start positions (0: length unknown).
    z: Vec<f64>,
    /// One-candidate molecules the law was learned from.
    pub molecules: usize,
}

impl JunctionLaw {
    /// The largest junction at or before `start`, or 0 (the transcript's start).
    fn previous_junction(junctions: &[u32], start: u32) -> u32 {
        match junctions.partition_point(|&j| j <= start) {
            0 => 0,
            k => junctions[k - 1],
        }
    }

    /// Learns the law from the one-candidate molecules of `cells` (each counting once, over its
    /// intervals).
    pub fn learn(cells: &[CellMolecules], junctions: &[Vec<u32>]) -> Self {
        let n_refs = junctions.len();
        let mut len_of = vec![0u32; n_refs];
        let mut max_len = 1u32;
        for cell in cells {
            for a in &cell.alns {
                if let Some((_, _, len)) = a.span {
                    if let Some(l) = len_of.get_mut(cell.target_ids[a.target] as usize) {
                        *l = (*l).max(len);
                    }
                    max_len = max_len.max(len);
                }
            }
        }
        let nb = GapLaw::bin(max_len) + 1;
        let mut counts = vec![0.0; nb];
        let mut molecules = 0usize;
        for cell in cells {
            for m in 0..cell.molecule_count() {
                let molecule = cell.molecule(m);
                let target = molecule[0].target;
                if !molecule.iter().all(|a| a.target == target) {
                    continue;
                }
                let t = cell.target_ids[target] as usize;
                if t >= n_refs {
                    continue;
                }
                let k = molecule.iter().filter(|a| a.span.is_some()).count();
                for a in molecule {
                    if let Some((start, _, _)) = a.span {
                        let d = start - Self::previous_junction(&junctions[t], start);
                        counts[GapLaw::bin(d).min(nb - 1)] += 1.0 / k as f64;
                    }
                }
                molecules += usize::from(k > 0);
            }
        }
        let total: f64 = counts.iter().sum();
        let dens: Vec<f64> = (0..nb)
            .map(|b| (counts[b] + 0.5) / ((total + 0.5 * nb as f64) * GapLaw::width(b)))
            .collect();
        let mut cum = vec![0.0; max_len as usize + 2];
        for d in 0..=max_len as usize {
            cum[d + 1] = cum[d] + dens[GapLaw::bin(d as u32).min(nb - 1)];
        }
        let z = (0..n_refs)
            .map(|t| {
                let len = len_of[t];
                if len == 0 {
                    return 0.0;
                }
                let (mut z, mut prev) = (0.0, 0u32);
                for &j in &junctions[t] {
                    if j >= len {
                        break;
                    }
                    if j > prev {
                        z += cum[(j - prev) as usize];
                        prev = j;
                    }
                }
                z + cum[(len - prev) as usize]
            })
            .collect();
        JunctionLaw {
            junctions: junctions.to_vec(),
            dens,
            cum,
            z,
            molecules,
        }
    }

    /// The law's probability per base that a read on transcript `t` starts at `start`; 1 when
    /// the transcript's length is unknown to the law.
    pub fn factor(&self, t: u32, start: u32) -> f64 {
        let Some(&z) = self.z.get(t as usize) else {
            return 1.0;
        };
        if z <= 0.0 {
            return 1.0;
        }
        let d = start - Self::previous_junction(&self.junctions[t as usize], start);
        self.dens[GapLaw::bin(d).min(self.dens.len() - 1)] / z
    }

    /// Multiplies every alignment's prior in `cells` by its factor.
    fn apply(&self, cells: &mut [CellMolecules]) {
        for cell in cells.iter_mut() {
            for a in cell.alns.iter_mut() {
                if let Some((start, _, _)) = a.span {
                    a.prior *= self.factor(cell.target_ids[a.target], start);
                }
            }
        }
    }
}

/// Learns a [`JunctionLaw`] from `cells` and applies it to them.
fn learn_junction(
    cells: &mut [CellMolecules],
    junctions: Option<&[Vec<u32>]>,
) -> Option<std::sync::Arc<JunctionLaw>> {
    let law = JunctionLaw::learn(cells, junctions?);
    law.apply(cells);
    Some(std::sync::Arc::new(law))
}

/// One cell's molecules as the EM reads them: molecule `m` owns alignments `off[m]..off[m + 1]`,
/// whose targets index the cell's transcripts `0..n_targets` (transcript ids `target_ids`).
pub struct CellMolecules {
    off: Vec<u32>,
    alns: Vec<Alignment>,
    n_targets: usize,
    target_ids: Vec<u32>,
    /// The cell's group for pooled abundances (0 without `--truncation-groups`).
    group: u32,
    /// Whether the positional likelihood uses the alignments' spans; false under interval mode
    /// `none`, where the spans are kept for the coverage and junction factors only.
    positional: bool,
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

    /// The cell's group for pooled abundances.
    pub fn set_group(&mut self, group: u32) {
        self.group = group;
    }

    /// Molecules with one candidate transcript, per transcript: `(transcript id, molecules)`,
    /// sorted by transcript id, for transcripts with at least one such molecule.
    pub fn one_candidate_counts(&self) -> Vec<(u32, u32)> {
        let mut unique = vec![0_u32; self.n_targets];
        for m in 0..self.molecule_count() {
            let mol = self.molecule(m);
            let first = mol[0].target;
            if mol.iter().all(|a| a.target == first) {
                unique[first] += 1;
            }
        }
        let mut out: Vec<(u32, u32)> = unique
            .iter()
            .enumerate()
            .filter(|&(_, &u)| u > 0)
            .map(|(i, &u)| (self.target_ids[i], u))
            .collect();
        out.sort_unstable();
        out
    }
}

/// One cell's molecules as the EM reads them, its transcript ids (target `i` is `ids[i]`), and
/// the number of molecules mixing candidates with and without an interval (see
/// [`TruncationFit::mixed_interval_molecules`]).
fn build(
    molecules: &[ResolvedMolecule],
    interval_mode: TruncationIntervalMode,
    factors: Option<&PriorFactors>,
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
            // a coverage model multiplies the candidate's prior by how well its coverage fits
            let prior = match factors.and_then(|f| f.coverage) {
                Some(cov) => prior * coverage_weight(molecule, target_id, cov),
                None => prior,
            };
            // spans are kept under every mode: the coverage and junction factors use them even where
            // the positional likelihood ignores them
            let intervals = molecule.intervals.get(&target_id);
            if let Some(intervals) = intervals {
                let valid_intervals: Vec<_> = intervals
                    .iter()
                    .filter(|interval| interval.transcript_len > 0)
                    .collect();
                if !valid_intervals.is_empty() {
                    let transcript_len = valid_intervals[0].transcript_len;
                    let positions = match interval_mode {
                        TruncationIntervalMode::All | TruncationIntervalMode::Ignore => valid_intervals
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
                        // The record's start is 1-based (its end, 1-based inclusive, is already
                        // the 0-based exclusive end); the model works on 0-based gaps.
                        let start = start.saturating_sub(1).min(len.saturating_sub(1));
                        let end = end.clamp(start, len);
                        alns.push(Alignment {
                            target,
                            prior: match factors.and_then(|f| f.junction) {
                                Some(law) => interval_prior * law.factor(target_id, start),
                                None => interval_prior,
                            },
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
    target_ids.shrink_to_fit();
    let n_targets = target_ids.len();
    (
        CellMolecules {
            off,
            alns,
            n_targets,
            target_ids: target_ids.clone(),
            group: 0,
            positional: interval_mode != TruncationIntervalMode::Ignore,
        },
        target_ids,
        mixed_interval_molecules,
    )
}

/// One cell's molecules in the form [`fit_global`] takes.
pub fn cell_molecules(
    molecules: &[ResolvedMolecule],
    interval_mode: TruncationIntervalMode,
    group: u32,
    factors: Option<&PriorFactors>,
) -> CellMolecules {
    let mut cell = build(molecules, interval_mode, factors).0;
    cell.group = group;
    cell
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
            w: cell
                .alns
                .iter()
                .map(|a| a.fixed_weight(par, e, cell.positional))
                .collect(),
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
            let entry = a.weight(prev, par, e, cell.positional);
            denom += entry.0;
            scored.push(entry);
        }
        if denom > EM_DENOM_THRESH {
            st.ll += denom.ln();
            for (a, (w, positional)) in molecule.iter().zip(&scored) {
                let post = w / denom;
                curr[a.target] += post;
                if let Some((g, t)) = positional {
                    st.add(g, t, post, par.class_of(g.len), par.len3.len());
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
#[allow(clippy::too_many_arguments)]
pub fn fit_global(
    cells: &mut [CellMolecules],
    max_len: u32,
    share: f64,
    n_threads: usize,
    spec: TruncationModelSpec,
    abundance: TruncationAbundance,
    n_refs: usize,
    n_groups: usize,
    extras: FitExtras,
) -> GlobalFit {
    let junction = learn_junction(cells, extras.junctions);
    if abundance != TruncationAbundance::Cell {
        return fit_pooled(
            cells,
            max_len,
            share,
            n_threads,
            spec,
            n_refs,
            n_groups.max(1),
            extras.genes,
            junction,
        );
    }
    let max_len = cells
        .iter()
        .map(CellMolecules::max_len)
        .fold(max_len.max(1), u32::max);
    let mut par = Params::initial(max_len, spec);
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

    if std::env::var("TRUNC_PRINT_LAWS").is_ok() {
        eprintln!("LAWS K5 rho {:.4}", par.k5.rho);
        eprintln!("LAWS K5 abs {:?}", par.k5.abs.h);
        eprintln!("LAWS K5 rel {:?}", par.k5.rel);
        eprintln!("LAWS K3 rho {:.4}", par.k3.rho);
        eprintln!("LAWS K3 abs {:?}", par.k3.abs.h);
        eprintln!("LAWS K3 rel {:?}", par.k3.rel);
        eprintln!("LAWS p3 {:.5} p5 {:.5} pf {:.5} pu {:.5}", par.p3, par.p5, par.pf, par.pu);
    }
    GlobalFit {
        params: TruncationParams::new(&par),
        par,
        cells: cells.len(),
        molecules: cells.iter().map(CellMolecules::molecule_count).sum(),
        iterations,
        frozen_at,
        freeze_reason,
        alpha: None,
        locus: None,
        abundance_iterations: 0,
        calibration: None,
        prior: None,
        junction,
        group_prior: None,
        locus_mass_cache: std::sync::OnceLock::new(),
    }
}

/// Upper bound on the EM steps of the pooled-abundance fit once the parameters are fixed, and
/// its convergence threshold (largest relative change of an abundance above
/// [`POOLED_MIN_ABUNDANCE`] molecules).
const MAX_POOLED_ITER: usize = 2000;
const POOLED_CONVERGENCE: f64 = 1e-4;
const POOLED_MIN_ABUNDANCE: f64 = 1e-2;

/// Splits `cells` into chunks for `n_threads` threads.
fn chunk_len(n: usize, n_threads: usize) -> usize {
    n.div_ceil(n_threads.max(1)).max(1)
}

/// One step of the parameter-estimating EM with pooled abundances: every cell's molecules are
/// distributed over their alignments under the one abundance vector `alpha` (indexed by
/// transcript id). Returns the positional statistics and the next abundances.
fn pooled_step(
    cells: &[CellMolecules],
    par: &Params,
    e: &Exps,
    alpha: &[Vec<f64>],
    n_threads: usize,
) -> (Stats, Vec<Vec<f64>>) {
    let (n_groups, n_refs) = (alpha.len(), alpha[0].len());
    let n_classes = par.len3.len();
    std::thread::scope(|scope| {
        let handles: Vec<_> = cells
            .chunks(chunk_len(cells.len(), n_threads))
            .map(|cs| {
                scope.spawn(move || {
                    let mut st = Stats::default();
                    let mut next_all = vec![vec![0.0; n_refs]; n_groups];
                    let mut scored: Vec<(f64, Option<(Geo, Terms)>)> = Vec::new();
                    for cell in cs {
                        let gi = (cell.group as usize).min(n_groups - 1);
                        let alpha = &alpha[gi];
                        let next = &mut next_all[gi];
                        for m in 0..cell.molecule_count() {
                            let molecule = cell.molecule(m);
                            scored.clear();
                            let mut denom = 0.0;
                            for a in molecule {
                                let ab = alpha[cell.target_ids[a.target] as usize];
                                let entry = a.weight_at(ab, par, e, cell.positional);
                                denom += entry.0;
                                scored.push(entry);
                            }
                            if denom > EM_DENOM_THRESH {
                                st.ll += denom.ln();
                                for (a, (w, positional)) in molecule.iter().zip(&scored) {
                                    let post = w / denom;
                                    next[cell.target_ids[a.target] as usize] += post;
                                    if let Some((g, t)) = positional {
                                        st.add(g, t, post, par.class_of(g.len), n_classes);
                                    }
                                }
                            }
                        }
                    }
                    (st, next_all)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("truncation EM thread panicked"))
            .fold(
                (Stats::default(), vec![vec![0.0; n_refs]; n_groups]),
                |(st, mut acc), (s, v)| {
                    for (a, b) in acc.iter_mut().zip(v) {
                        a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
                    }
                    (st.merge(s), acc)
                },
            )
    })
}

/// The fixed (abundance-free) weights of every cell's alignments, with transcript ids for
/// targets, for the pooled-abundance EM once the parameters are fixed.
fn pooled_fixed(cells: &[CellMolecules], par: &Params, e: &Exps, n_threads: usize) -> Vec<Fixed> {
    std::thread::scope(|scope| {
        let handles: Vec<_> = cells
            .chunks(chunk_len(cells.len(), n_threads))
            .map(|cs| {
                scope.spawn(move || {
                    cs.iter()
                        .map(|cell| {
                            let mut fx = Fixed::new(cell, par, e);
                            fx.t.iter_mut().for_each(|t| *t = cell.target_ids[*t] as usize);
                            fx
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("truncation EM thread panicked"))
            .collect()
    })
}

/// One pooled EM step with fixed weights: the next abundances, and (with `ll`) the
/// log-likelihood of `alpha` normalized to proportions, up to a constant.
fn pooled_static(fx: &[&Fixed], alpha: &[f64], n_threads: usize, ll: bool) -> (Vec<f64>, f64) {
    let n_refs = alpha.len();
    let total: f64 = alpha.iter().sum();
    std::thread::scope(|scope| {
        let handles: Vec<_> = fx
            .chunks(chunk_len(fx.len(), n_threads))
            .map(|fs| {
                scope.spawn(move || {
                    let mut next = vec![0.0; n_refs];
                    let (mut l, mut n) = (0.0, 0.0);
                    for f in fs {
                        for m in 0..f.molecules() {
                            let (lo, hi) = (f.off[m], f.off[m + 1]);
                            let denom: f64 = (lo..hi).map(|k| alpha[f.t[k]] * f.w[k]).sum();
                            if denom > EM_DENOM_THRESH {
                                if ll {
                                    l += denom.ln();
                                    n += 1.0;
                                }
                                for k in lo..hi {
                                    next[f.t[k]] += alpha[f.t[k]] * f.w[k] / denom;
                                }
                            }
                        }
                    }
                    (next, l - n * total.ln())
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("truncation EM thread panicked"))
            .fold((vec![0.0; n_refs], 0.0), |(mut acc, l), (v, lv)| {
                acc.iter_mut().zip(v).for_each(|(a, b)| *a += b);
                (acc, l + lv)
            })
    })
}

/// Pooled abundances with the parameters fixed: EM over the molecules of all cells, one
/// abundance per transcript, accelerated with SQUAREM (as [`squarem_cycle`]), until the largest
/// relative change is below [`POOLED_CONVERGENCE`] or [`MAX_POOLED_ITER`] steps.
fn pooled_abundance(
    fx: &[&Fixed],
    mut alpha: Vec<f64>,
    n_threads: usize,
    pseudo: Option<&[f64]>,
    locus: &[u32],
) -> (Vec<f64>, usize) {
    // One EM step, or with `pseudo` (pseudo-molecules per transcript, a Dirichlet prior on each
    // locus's shares) one MAP step: the expected molecules plus the pseudo-molecules, scaled back
    // to the locus's own molecules so that every locus keeps the molecule scale and a molecule
    // whose candidates span two loci is not drawn to the one with the stronger prior; with `ll`
    // the log-likelihood, plus the log-prior, up to a constant.
    let n_refs = alpha.len();
    let step = |a: &[f64], ll: bool| -> (Vec<f64>, f64) {
        let (mut next, mut l) = pooled_static(fx, a, n_threads, ll);
        if let Some(p) = pseudo {
            if ll {
                let total: f64 = a.iter().sum::<f64>().max(1e-300);
                l += p
                    .iter()
                    .zip(a)
                    .filter(|(p, _)| **p > 0.0)
                    .map(|(p, x)| p * ((x + POOLED_FLOOR) / total).ln())
                    .sum::<f64>();
            }
            let (mut n_l, mut b_l) = (vec![0.0; n_refs], vec![0.0; n_refs]);
            for t in 0..n_refs {
                n_l[locus[t] as usize] += next[t];
                b_l[locus[t] as usize] += p[t];
            }
            for t in 0..n_refs {
                let l = locus[t] as usize;
                let s = n_l[l] + b_l[l];
                if s > 0.0 {
                    next[t] = (next[t] + p[t]) * n_l[l] / s;
                }
            }
        }
        (next, l)
    };
    let mut steps = 0;
    let mut step_max = 1.0_f64;
    while steps < MAX_POOLED_ITER {
        let t0 = alpha.clone();
        let (t1, _) = step(&t0, false);
        let (t2, l2) = step(&t1, true);
        steps += 2;
        let (mut rr, mut vv) = (0.0, 0.0);
        for i in 0..t0.len() {
            let (r, v) = (t1[i] - t0[i], t2[i] - 2.0 * t1[i] + t0[i]);
            rr += r * r;
            vv += v * v;
        }
        let mut s = if vv > 0.0 { (rr / vv).sqrt() } else { 1.0 };
        if !s.is_finite() || s < 1.0 {
            s = 1.0;
        }
        let hit = s > 1.0 && s >= step_max;
        s = s.min(step_max);
        let (prev, next) = if s > 1.0 {
            let x: Vec<f64> = (0..t0.len())
                .map(|i| {
                    let (r, v) = (t1[i] - t0[i], t2[i] - 2.0 * t1[i] + t0[i]);
                    (t0[i] + 2.0 * s * r + s * s * v).max(0.0)
                })
                .collect();
            let (t3, _) = step(&x, false);
            let (_, l3) = step(&t3, true);
            steps += 2;
            if l3.is_finite() && l3 >= l2 {
                if hit {
                    step_max *= Squarem::STEP_GROW;
                }
                (x, t3)
            } else {
                step_max = (step_max / Squarem::STEP_GROW).max(1.0);
                (t1, t2)
            }
        } else {
            if hit {
                step_max *= Squarem::STEP_GROW;
            }
            (t1, t2)
        };
        let rel = prev
            .iter()
            .zip(&next)
            .filter(|(p, _)| **p > POOLED_MIN_ABUNDANCE)
            .map(|(p, c)| (c - p).abs() / p)
            .fold(0.0, f64::max);
        alpha = next;
        if rel < POOLED_CONVERGENCE {
            break;
        }
    }
    (alpha, steps)
}

/// Each transcript's locus: the connected components of the graph joining the candidates of
/// every molecule of every cell (a transcript no molecule has is a locus of its own).
fn loci(cells: &[CellMolecules], n_refs: usize) -> Vec<u32> {
    let mut parent: Vec<u32> = (0..n_refs as u32).collect();
    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            parent[x as usize] = parent[parent[x as usize] as usize];
            x = parent[x as usize];
        }
        x
    }
    for cell in cells {
        for m in 0..cell.molecule_count() {
            let molecule = cell.molecule(m);
            if let Some(first) = molecule.first() {
                let a = find(&mut parent, cell.target_ids[first.target]);
                for al in &molecule[1..] {
                    let b = find(&mut parent, cell.target_ids[al.target]);
                    if a != b {
                        parent[b as usize] = a;
                    }
                }
            }
        }
    }
    (0..n_refs as u32).map(|x| find(&mut parent, x)).collect()
}

/// Splits each locus by gene: the transcripts of one locus and one gene form a group, labelled by
/// the smallest transcript id in it (a valid locus label).
fn split_by_genes(locus: &[u32], genes: &[u32]) -> Vec<u32> {
    let mut label: HashMap<(u32, u32), u32> = HashMap::new();
    (0..locus.len())
        .map(|t| *label.entry((locus[t], genes[t])).or_insert(t as u32))
        .collect()
}

/// The strength of the pooled prior that the data support, from the molecules whose transcript
/// is certain: those with one candidate transcript. (Molecules with several candidates settled
/// by their alignment scores are not used: their labels are wrong often enough, and unevenly
/// enough across cells, to look like cell-specific isoform usage.)
///
/// Each cell's shares of a locus are taken as drawn from a Dirichlet around the pooled shares
/// `phi` with concentration `alpha`, so that `Var(share_t) = phi_t (1 - phi_t) / (alpha + 1)`.
/// For a cell with `n` one-candidate molecules in the locus, `x` of them on transcript `t`,
///
/// ```text
/// E[(x - n phi_t)^2] = n phi_t (1 - phi_t) + n (n - 1) Var(share_t)
/// ```
///
/// which gives one `alpha` for all loci by the method of moments, each cell weighted by its
/// number of molecule pairs `n (n - 1)`:
///
/// ```text
/// alpha + 1 = sum_locus sum_t D phi_t (1 - phi_t) / sum_locus sum_t N_t,
/// N_t = sum_cells [(x - n phi_t)^2 - n phi_t (1 - phi_t)],   D = sum_cells n (n - 1),
/// ```
///
/// over the cells with `n >= 2`. A denominator of zero or less means no variation between
/// cells beyond sampling: `alpha` is infinite and every cell uses the pooled shares. Per-locus
/// estimates are too noisy at a few molecules per cell, which is why one value is estimated.
/// `locus` should be split by gene ([`split_by_genes`]) when genes are known: within a locus
/// that merges genes, the shares vary between cells with gene expression as well.
/// Weight of the common relative variance in a locus's own estimate under `auto-gene`, in the
/// units of the estimator's numerator (molecule pairs times `phi (1 - phi)`): about 1,000 pairs at
/// `phi = 1/2`. A locus with far more pairs follows its own data, one with fewer the common value.
const GENE_PRIOR_WEIGHT: f64 = 250.0;
/// Pseudo-molecules standing in for an infinite group-level strength: the group's shares are
/// then the shares of all cells to within one part in a billion.
const GROUP_SHRINK_MAX: f64 = 1e9;
/// Under `auto-gene`, a locus gets a strength of its own only when its deviation between the units
/// exceeds this many null standard deviations; the rest take the pooled shares.
const LOCUS_Z_MIN: f64 = 4.0;

fn prior_strength(
    cells: &[CellMolecules],
    locus: &[u32],
    n_refs: usize,
    within_genes: bool,
    per_gene: bool,
) -> PriorStrength {
    let units: Vec<Vec<(u32, f64)>> = cells.iter().map(unique_counts).collect();
    prior_strength_counts(&units, locus, n_refs, within_genes, per_gene)
}

/// A cell's one-candidate molecules per transcript, as `(transcript id, count)` pairs.
fn unique_counts(cell: &CellMolecules) -> Vec<(u32, f64)> {
    let mut count = vec![0.0; cell.n_targets];
    for m in 0..cell.molecule_count() {
        let molecule = cell.molecule(m);
        let target = molecule[0].target;
        if molecule.iter().all(|a| a.target == target) {
            count[target] += 1.0;
        }
    }
    count
        .iter()
        .enumerate()
        .filter(|&(_, &c)| c > 0.0)
        .map(|(i, &c)| (cell.target_ids[i], c))
        .collect()
}

/// [`prior_strength`] over any units of pooling: each unit's one-candidate molecules per
/// transcript (cells for the cell-level prior, groups of cells for the group-level one).
fn prior_strength_counts(
    units: &[Vec<(u32, f64)>],
    locus: &[u32],
    n_refs: usize,
    within_genes: bool,
    per_gene: bool,
) -> PriorStrength {
    // per transcript, over the units with n >= 2 in its locus: sum x^2, sum x n; over all units: sum x
    let mut sx = vec![0.0; n_refs];
    let mut sx2 = vec![0.0; n_refs];
    let mut sxn = vec![0.0; n_refs];
    // per locus: sum n and sum n^2 over all units; sum n, n^2 and n^3 over the units with n >= 2
    let mut sn_all = vec![0.0; n_refs];
    let mut sn2_all = vec![0.0; n_refs];
    let mut sn = vec![0.0; n_refs];
    let mut sn2 = vec![0.0; n_refs];
    let mut sn3 = vec![0.0; n_refs];
    let (mut used_units, mut used_molecules) = (0_usize, 0_usize);
    let mut per_locus: HashMap<u32, f64> = HashMap::new();
    for unit in units {
        per_locus.clear();
        for &(t, c) in unit {
            *per_locus.entry(locus[t as usize]).or_insert(0.0) += c;
        }
        let mut counted = false;
        for &(t, c) in unit {
            let t = t as usize;
            let n = per_locus[&locus[t]];
            sx[t] += c;
            if n >= 2.0 {
                sx2[t] += c * c;
                sxn[t] += c * n;
                used_molecules += c as usize;
                counted = true;
            }
        }
        for (&l, &n) in &per_locus {
            sn_all[l as usize] += n;
            sn2_all[l as usize] += n * n;
            if n >= 2.0 {
                sn[l as usize] += n;
                sn2[l as usize] += n * n;
                sn3[l as usize] += n * n * n;
            }
        }
        used_units += usize::from(counted);
    }
    // loci with at least two transcripts that have one-candidate molecules
    let mut with_counts = vec![0_u32; n_refs];
    for t in 0..n_refs {
        if sx[t] > 0.0 {
            with_counts[locus[t] as usize] += 1;
        }
    }
    // The pooled share is estimated from the same units, which pulls every unit's deviation
    // toward zero, the more the larger the unit's part of the locus. With N and S2 the locus's
    // molecules and squared molecules over all units, the deviations of the units with n >= 2
    // have expectation v (sum n - sum n^2 / N) + Var(share) D, with
    // D = sum n (n - 1) - 2 sum n^2 (n - 1) / N + sum n^2 (S2 - N) / N^2.
    // For many small units this is the usual v sum n + Var(share) sum n (n - 1); for one small
    // unit against one large one, D is about four times the small unit's squared size instead
    // of the large one's, which is what makes the group-level estimate usable.
    let (mut num, mut den) = (0.0, 0.0);
    let (mut lnum, mut lden) = (vec![0.0; n_refs], vec![0.0; n_refs]);
    // per locus: sum over its transcripts of v^2, for the null variance of its deviation
    let mut lv2 = vec![0.0; n_refs];
    let pairs = |l: usize| -> f64 {
        let n_tot = sn_all[l];
        (sn2[l] - sn[l]) - 2.0 * (sn3[l] - sn2[l]) / n_tot
            + sn2[l] * (sn2_all[l] - n_tot) / (n_tot * n_tot)
    };
    for t in 0..n_refs {
        let l = locus[t] as usize;
        if sx[t] <= 0.0 || with_counts[l] < 2 || sn_all[l] <= 0.0 {
            continue;
        }
        let d = pairs(l);
        if d <= 0.0 {
            continue;
        }
        let phi = sx[t] / sn_all[l];
        let v = phi * (1.0 - phi);
        let e = sn[l] - sn2[l] / sn_all[l];
        let excess = sx2[t] - 2.0 * phi * sxn[t] + phi * phi * sn2[l] - v * e;
        num += d * v;
        den += excess;
        lnum[l] += d * v;
        lden[l] += excess;
        lv2[l] += v * v;
    }
    // never below one pseudo-molecule: a strength near zero would make a unit's shares its raw
    // counts, with nothing for the isoforms it happens to have no molecule of
    let alpha = if den > 0.0 && num > 0.0 {
        (num / den - 1.0).max(1.0)
    } else {
        f64::INFINITY
    };
    // per locus: its relative between-unit variance, for the loci whose deviation is beyond what
    // sampling explains (under no variation a locus's deviation has mean 0 and variance about
    // 2 D sum_t v_t^2), shrunk toward the common one of those loci; the rest pool fully. Without
    // the gate, a strength estimated from two or three units is noise on every locus that does
    // not vary, and the units overfit their own sampling there.
    let per_locus: Option<Vec<f64>> = per_gene.then(|| {
        let gated: Vec<bool> = (0..n_refs)
            .map(|l| lnum[l] > 0.0 && lden[l] > LOCUS_Z_MIN * (2.0 * pairs(l) * lv2[l]).sqrt())
            .collect();
        let (gn, gd) = (0..n_refs)
            .filter(|&l| gated[l])
            .fold((0.0, 0.0), |(a, b), l| (a + lnum[l], b + lden[l]));
        let common = if gn > 0.0 { (gd / gn).max(0.0) } else { 0.0 };
        (0..n_refs)
            .map(|l| {
                if !gated[l] {
                    return f64::INFINITY;
                }
                let r = (lden[l] + GENE_PRIOR_WEIGHT * common) / (lnum[l] + GENE_PRIOR_WEIGHT);
                if r <= 0.0 {
                    f64::INFINITY
                } else {
                    (1.0 / r - 1.0).max(1.0)
                }
            })
            .collect()
    });
    let mut finite: Vec<f64> = per_locus
        .as_ref()
        .map(|a| a.iter().copied().filter(|x| x.is_finite()).collect())
        .unwrap_or_default();
    finite.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let alpha_median = finite.get(finite.len() / 2).copied().unwrap_or(f64::INFINITY);
    PriorStrength {
        alpha,
        loci: (0..n_refs)
            .filter(|&l| locus[l] as usize == l && with_counts[l] >= 2 && sn_all[l] > 0.0 && pairs(l) > 0.0)
            .count(),
        cells: used_units,
        molecules: used_molecules,
        within_genes,
        loci_finite: finite.len(),
        alpha_median,
        per_locus,
    }
}

/// [`fit_global`] with pooled abundances: the parameter-estimating EM distributes every cell's
/// molecules under one abundance vector shared by all cells (instead of one per cell), and once
/// the parameters freeze the abundances are re-estimated with them fixed, to convergence. With
/// `genes` (each transcript's gene, [`GeneMap`]) the loci the prior acts within are split by gene.
#[allow(clippy::too_many_arguments)]
fn fit_pooled(
    cells: &[CellMolecules],
    max_len: u32,
    share: f64,
    n_threads: usize,
    spec: TruncationModelSpec,
    n_refs: usize,
    n_groups: usize,
    genes: Option<&[u32]>,
    junction: Option<std::sync::Arc<JunctionLaw>>,
) -> GlobalFit {
    let max_len = cells
        .iter()
        .map(CellMolecules::max_len)
        .fold(max_len.max(1), u32::max);
    let mut par = Params::initial(max_len, spec);
    // With a group prior the parameters and the all-cell abundances are fitted as one pool, exactly
    // as a run without groups fits them; the groups come after.
    let group_prior_on = n_groups > 1 && (spec.group_shrink > 0.0 || spec.group_shrink_auto);
    let n_param_groups = if group_prior_on { 1 } else { n_groups };
    let mut alpha = vec![vec![0.0; n_refs]; n_param_groups];
    for cell in cells {
        let each = cell.molecule_count() as f64 / cell.n_targets.max(1) as f64;
        let g = if group_prior_on { 0 } else { cell.group as usize };
        for &t in &cell.target_ids {
            alpha[g][t as usize] += each;
        }
    }
    let mut freeze = Freeze::new(share);
    let (mut frozen_at, mut freeze_reason) = (None, None);
    let mut last_ll = f64::NEG_INFINITY;
    let mut iterations = 0_usize;
    while iterations < MAX_EM_ITER as usize {
        let e = Exps::new(&par, max_len);
        let (st, next_alpha) = pooled_step(cells, &par, &e, &alpha, n_threads);
        alpha = next_alpha;
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
    let mut locus = loci(cells, n_refs);
    if let Some(g) = genes {
        locus = split_by_genes(&locus, g);
    }
    // Without a group prior: each group's pooled EM on its own. With one: the all-cell fit as a
    // run without groups makes it, each group's own fit started from the all-cell shares (where
    // the likelihood is flat between near-identical isoforms the group then stays where all
    // cells are instead of wandering along the ridge), the deviance of those fits for the
    // strength ([`group_deviance`], [`group_strength`]), and each group's pooled EM once more as
    // a MAP with `beta` pseudo-molecules per locus drawn from the all-cell shares (given, or
    // estimated: infinite, and the group's shares those of all cells, where the groups do not
    // differ beyond sampling).
    let mut group_pseudo: Option<Vec<f64>> = None;
    let (mut alpha, mut abundance_iterations, group_prior) = if !group_prior_on {
        let (a, it) = group_abundances(cells, &par, max_len, alpha, n_threads, None, &locus);
        (a, it, None)
    } else {
        let (all, it0) = group_abundances(cells, &par, max_len, alpha, n_threads, None, &locus);
        let global = all.into_iter().next().expect("the all-cell abundances");
        let mut mass = vec![0.0; n_refs];
        for (t, a) in global.iter().enumerate() {
            mass[locus[t] as usize] += a + POOLED_FLOOR;
        }
        // each group's molecules per locus: its starting point is those molecules at the
        // all-cell shares
        let mut n_gl = vec![vec![0.0; n_refs]; n_groups];
        for cell in cells {
            let g = (cell.group as usize).min(n_groups - 1);
            for m in 0..cell.molecule_count() {
                if let Some(first) = cell.molecule(m).first() {
                    let l = locus[cell.target_ids[first.target] as usize] as usize;
                    n_gl[g][l] += 1.0;
                }
            }
        }
        let start: Vec<Vec<f64>> = (0..n_groups)
            .map(|g| {
                (0..n_refs)
                    .map(|t| {
                        let l = locus[t] as usize;
                        n_gl[g][l] * (global[t] + POOLED_FLOOR) / mass[l]
                    })
                    .collect()
            })
            .collect();
        let (own, it1) = group_abundances(cells, &par, max_len, start.clone(), n_threads, None, &locus);
        let e = Exps::new(&par, max_len);
        let fx = pooled_fixed(cells, &par, &e, n_threads);
        let prior = spec.group_shrink_auto.then(|| {
            let (dev, n, touched) = group_deviance(cells, &fx, &own, &global, &locus, n_threads);
            group_strength(&dev, &n, &touched, &global, &locus, genes.is_some(), spec.group_shrink_per_gene)
        });
        let pseudo: Vec<f64> = (0..n_refs)
            .map(|t| {
                let beta = match prior.as_ref() {
                    Some(p) => p.per_locus.as_ref().map_or(p.alpha, |pl| pl[locus[t] as usize]),
                    None => spec.group_shrink,
                };
                beta.min(GROUP_SHRINK_MAX) * (global[t] + POOLED_FLOOR) / mass[locus[t] as usize]
            })
            .collect();
        let (map, it2) = group_abundances(cells, &par, max_len, start, n_threads, Some(&pseudo), &locus);
        group_pseudo = Some(pseudo);
        (map, it0 + it1 + it2, prior)
    };
    let prior = spec
        .shrink_auto
        .then(|| prior_strength(cells, &locus, n_refs, genes.is_some(), spec.shrink_per_gene));
    let mut calibration = None;
    if spec.calibrate_power {
        // The power depends on the abundances and the abundances on the power: two rounds
        // settle both (the second changes the power by far less than the first).
        for round in 1..=CALIBRATION_ROUNDS {
            let mut cal = calibrate_power(cells, &par, max_len, &alpha, n_threads);
            cal.rounds = round;
            par.pos_power = cal.p;
            (alpha, abundance_iterations) =
                group_abundances(cells, &par, max_len, alpha, n_threads, group_pseudo.as_deref(), &locus);
            calibration = Some(cal);
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
        locus: Some(locus),
        alpha: Some(alpha),
        abundance_iterations,
        calibration,
        prior,
        junction,
        group_prior,
        locus_mass_cache: std::sync::OnceLock::new(),
    }
}

/// With the parameters fixed, each group's pooled abundances (independent across groups).
fn group_abundances(
    cells: &[CellMolecules],
    par: &Params,
    max_len: u32,
    alpha: Vec<Vec<f64>>,
    n_threads: usize,
    pseudo: Option<&[f64]>,
    locus: &[u32],
) -> (Vec<Vec<f64>>, usize) {
    let e = Exps::new(par, max_len);
    let fx = pooled_fixed(cells, par, &e, n_threads);
    let mut iterations = 0;
    let n_alpha = alpha.len();
    let alpha = alpha
        .into_iter()
        .enumerate()
        .map(|(gi, a)| {
            let sub: Vec<&Fixed> = fx
                .iter()
                .zip(cells)
                .filter(|(_, c)| (c.group as usize).min(n_alpha - 1) == gi)
                .map(|(f, _)| f)
                .collect();
            if sub.is_empty() {
                return a;
            }
            let (a, it) = pooled_abundance(&sub, a, n_threads, pseudo, locus);
            iterations = iterations.max(it);
            a
        })
        .collect();
    (alpha, iterations)
}

/// Per group and locus: the deviance of the group's own pooled shares against the shares of all
/// cells on the group's molecules, `2 sum_m [ln P(m | group) - ln P(m | all cells)]` with
/// `P(m | a) = sum_t a_t w_m,t / sum_(t in locus) a_t`, and the group's molecules in the locus.
/// Every molecule counts, ambiguous ones included, so a group that lacks an isoform whose reads
/// are never unique still shows it. Also returned: which transcripts each group's molecules have
/// among their candidates, since a group's fit can only move on those. `fx` are the cells'
/// fixed weights ([`pooled_fixed`]).
fn group_deviance(
    cells: &[CellMolecules],
    fx: &[Fixed],
    alpha: &[Vec<f64>],
    global: &[f64],
    locus: &[u32],
    n_threads: usize,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<bool>>) {
    let n_refs = global.len();
    let n_groups = alpha.len();
    let mass_of = |a: &[f64]| {
        let mut m = vec![0.0; n_refs];
        for (t, x) in a.iter().enumerate() {
            m[locus[t] as usize] += x + POOLED_FLOOR;
        }
        m
    };
    let mass_g: Vec<Vec<f64>> = alpha.iter().map(|a| mass_of(a)).collect();
    let mass_0 = mass_of(global);
    let chunk = chunk_len(cells.len(), n_threads);
    std::thread::scope(|scope| {
        let handles: Vec<_> = fx
            .chunks(chunk)
            .zip(cells.chunks(chunk))
            .map(|(fs, cs)| {
                let (mass_g, mass_0) = (&mass_g, &mass_0);
                scope.spawn(move || {
                    let mut dev = vec![vec![0.0; n_refs]; n_groups];
                    let mut n = vec![vec![0.0; n_refs]; n_groups];
                    let mut touched = vec![vec![false; n_refs]; n_groups];
                    for (f, cell) in fs.iter().zip(cs) {
                        let g = (cell.group as usize).min(n_groups - 1);
                        let a = &alpha[g];
                        for m in 0..f.molecules() {
                            let (lo, hi) = (f.off[m], f.off[m + 1]);
                            if lo == hi {
                                continue;
                            }
                            let l = locus[f.t[lo]] as usize;
                            let (mut sg, mut s0) = (0.0, 0.0);
                            for k in lo..hi {
                                let t = f.t[k];
                                touched[g][t] = true;
                                sg += (a[t] + POOLED_FLOOR) / mass_g[g][l] * f.w[k];
                                s0 += (global[t] + POOLED_FLOOR) / mass_0[l] * f.w[k];
                            }
                            if sg > EM_DENOM_THRESH && s0 > EM_DENOM_THRESH {
                                dev[g][l] += 2.0 * (sg.ln() - s0.ln());
                                n[g][l] += 1.0;
                            }
                        }
                    }
                    (dev, n, touched)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().expect("truncation EM thread panicked")).fold(
            (
                vec![vec![0.0; n_refs]; n_groups],
                vec![vec![0.0; n_refs]; n_groups],
                vec![vec![false; n_refs]; n_groups],
            ),
            |(mut d, mut n, mut tch), (dv, nv, tv)| {
                for g in 0..n_groups {
                    d[g].iter_mut().zip(&dv[g]).for_each(|(a, b)| *a += b);
                    n[g].iter_mut().zip(&nv[g]).for_each(|(a, b)| *a += b);
                    tch[g].iter_mut().zip(&tv[g]).for_each(|(a, b)| *a |= b);
                }
                (d, n, tch)
            },
        )
    })
}

/// The group-level prior strength from the groups' deviances ([`group_deviance`]). A group's
/// fit can move only on the `k_g` transcripts of a locus its molecules have as candidates, so
/// with `k` the locus's transcripts of nonzero abundance the deviances of `G` groups sum to about
/// `df = sum_g (k_g - 1) - (k - 1)` when the groups share the locus's isoform usage (the
/// likelihood-ratio test of homogeneity), and to that plus `sum_g n_g (k_g - 1) / (beta + 1)`
/// when their shares spread around the all-cell shares as a Dirichlet of concentration `beta`,
/// so that
///
/// ```text
/// beta + 1 = sum_g n_g (k_g - 1) / (sum_g 2 dev_g - df),
/// ```
///
/// per locus. A locus gets a strength of its own only when its deviance exceeds the chi-square
/// null by [`LOCUS_Z_MIN`] standard deviations (Wilson-Hilferty); the rest are infinite, and
/// their groups take the shares of all cells. (One strength for all loci as a ratio of sums is
/// not usable: the thousands of loci with few molecules per group sit below their nominal
/// degrees of freedom and swamp the few hundred that differ.) Strengths are floored at one
/// pseudo-molecule; `alpha` reports the median of the finite ones, infinite when there are none.
fn group_strength(
    dev: &[Vec<f64>],
    n: &[Vec<f64>],
    touched: &[Vec<bool>],
    global: &[f64],
    locus: &[u32],
    within_genes: bool,
    per_gene: bool,
) -> PriorStrength {
    let n_refs = global.len();
    let n_groups = dev.len();
    let mut k = vec![0_u32; n_refs];
    let mut kg = vec![vec![0_u32; n_refs]; n_groups];
    for t in 0..n_refs {
        if global[t] > 1e-12 {
            let l = locus[t] as usize;
            k[l] += 1;
            for g in 0..n_groups {
                if touched[g][t] {
                    kg[g][l] += 1;
                }
            }
        }
    }
    let (mut used_loci, mut used_molecules) = (0_usize, 0.0);
    let mut per_locus_beta = vec![f64::INFINITY; n_refs];
    let mut finite: Vec<f64> = Vec::new();
    for l in 0..n_refs {
        if locus[l] as usize != l || k[l] < 2 {
            continue;
        }
        let (mut sd, mut sn, mut scale, mut free, mut groups) = (0.0, 0.0, 0.0, 0.0, 0_usize);
        for g in 0..n_groups {
            if n[g][l] > 0.0 && kg[g][l] >= 1 {
                sd += dev[g][l];
                sn += n[g][l];
                scale += n[g][l] * (kg[g][l] as f64 - 1.0);
                free += kg[g][l] as f64 - 1.0;
                groups += 1;
            }
        }
        let df = free - (k[l] as f64 - 1.0);
        if groups < 2 || df < 1.0 || scale <= 0.0 {
            continue;
        }
        let excess = sd - df;
        used_loci += 1;
        used_molecules += sn;
        // Wilson-Hilferty: (sd / df)^(1/3) is normal with mean 1 - 2 / (9 df), var 2 / (9 df)
        let z = ((sd / df).max(0.0).cbrt() - (1.0 - 2.0 / (9.0 * df))) / (2.0 / (9.0 * df)).sqrt();
        if z > LOCUS_Z_MIN && excess > 0.0 {
            let beta = (scale / excess - 1.0).max(1.0);
            per_locus_beta[l] = beta;
            finite.push(beta);
        }
    }
    finite.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let alpha_median = finite.get(finite.len() / 2).copied().unwrap_or(f64::INFINITY);
    PriorStrength {
        alpha: alpha_median,
        loci: used_loci,
        cells: n_groups,
        molecules: used_molecules as usize,
        within_genes: within_genes || per_gene,
        loci_finite: finite.len(),
        alpha_median,
        per_locus: Some(per_locus_beta),
    }
}

/// Molecules whose alignment-score prior puts at least this share on one transcript are taken
/// as labelled by their splice structure in [`calibrate_power`].
const CALIBRATION_MIN_PRIOR: f64 = 0.99;
const CALIBRATION_ROUNDS: usize = 2;
const CALIBRATION_MAX_POWER: f64 = 2.0;

/// The calibrated positional power and how it was found.
#[derive(Clone, Debug)]
pub struct Calibration {
    /// The power maximizing the calibration log-likelihood.
    pub p: f64,
    /// Molecules used: at least two candidate transcripts, all with intervals, and a prior of at
    /// least [`CALIBRATION_MIN_PRIOR`] on one of them.
    pub molecules: usize,
    /// Share of them whose positional likelihood alone ranks the prior's transcript first.
    pub positional_agreement: f64,
    /// The calibration log-likelihood per molecule at a few powers, for reporting.
    pub curve: Vec<(f64, f64)>,
    pub rounds: usize,
}

/// The positional power `p` that best predicts, from abundance and position alone, the
/// transcript of molecules whose transcript is settled by their alignments (a prior of at least
/// [`CALIBRATION_MIN_PRIOR`] on one candidate):
///
/// ```text
/// p = argmax_p  sum_m [ log a_w(m) + p log L_m,w(m) - log sum_t a_t exp(p log L_m,t) ]
/// ```
///
/// with `a_t` the pooled abundance of candidate `t` (in the molecule's group), `L_m,t` the
/// positional likelihood of the molecule on `t` (the mean over its intervals) and `w(m)` the
/// candidate its prior settles on. This is temperature scaling of the positional evidence with
/// the splice structure as the label; the objective is concave in `p`, so the maximum is unique.
/// A correct positional model gives p near 1; a model whose likelihood ratios overstate the
/// evidence (real reads that the model's shapes describe only roughly) gives p below 1.
fn calibrate_power(
    cells: &[CellMolecules],
    par: &Params,
    max_len: u32,
    alpha: &[Vec<f64>],
    n_threads: usize,
) -> Calibration {
    let e = Exps::new(par, max_len);
    // per thread: molecule offsets, winner index, and per candidate (log a, log L)
    let parts: Vec<(Vec<u32>, Vec<u32>, Vec<f32>, Vec<f32>, usize)> = std::thread::scope(|scope| {
        let handles: Vec<_> = cells
            .chunks(chunk_len(cells.len(), n_threads))
            .map(|cs| {
                let e = &e;
                scope.spawn(move || {
                    let (mut off, mut win, mut la, mut ll) = (vec![0u32], Vec::new(), Vec::new(), Vec::new());
                    let mut agree = 0usize;
                    let mut per: Vec<(usize, f64, f64)> = Vec::new();
                    for cell in cs {
                        let a = &alpha[cell.group as usize];
                        for m in 0..cell.molecule_count() {
                            let mol = cell.molecule(m);
                            if mol.iter().any(|x| x.span.is_none()) {
                                continue;
                            }
                            // aggregate a molecule's intervals per target: prior and mean likelihood
                            per.clear();
                            for x in mol {
                                let g = x.geo().expect("span checked");
                                let lik = par.terms(&g, e).total();
                                match per.iter_mut().find(|(t, _, _)| *t == x.target) {
                                    Some(p) => {
                                        p.1 += x.prior;
                                        p.2 += x.prior * lik;
                                    }
                                    None => per.push((x.target, x.prior, x.prior * lik)),
                                }
                            }
                            if per.len() < 2 {
                                continue;
                            }
                            let tot: f64 = per.iter().map(|p| p.1).sum();
                            let (wi, wp) = per
                                .iter()
                                .enumerate()
                                .fold((0, f64::MIN), |b, (i, p)| if p.1 > b.1 { (i, p.1) } else { b });
                            if tot <= 0.0 || wp / tot < CALIBRATION_MIN_PRIOR {
                                continue;
                            }
                            let lw = per[wi].2 / per[wi].1;
                            if per.iter().all(|p| p.2 / p.1.max(1e-300) <= lw) {
                                agree += 1;
                            }
                            for p in &per {
                                let id = cell.target_ids[p.0] as usize;
                                la.push((a[id] + POOLED_FLOOR).ln() as f32);
                                ll.push((p.2 / p.1.max(1e-300)).max(1e-300).ln() as f32);
                            }
                            win.push(wi as u32);
                            off.push(la.len() as u32);
                        }
                    }
                    (off, win, la, ll, agree)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("calibration thread panicked"))
            .collect()
    });
    let objective = |p: f64| -> f64 {
        let mut f = 0.0;
        for (off, win, la, ll, _) in &parts {
            for m in 0..win.len() {
                let (lo, hi) = (off[m] as usize, off[m + 1] as usize);
                let mut mx = f64::NEG_INFINITY;
                for k in lo..hi {
                    mx = mx.max(la[k] as f64 + p * ll[k] as f64);
                }
                let lse = mx
                    + (lo..hi)
                        .map(|k| (la[k] as f64 + p * ll[k] as f64 - mx).exp())
                        .sum::<f64>()
                        .ln();
                let w = lo + win[m] as usize;
                f += la[w] as f64 + p * ll[w] as f64 - lse;
            }
        }
        f
    };
    let molecules: usize = parts.iter().map(|x| x.1.len()).sum();
    let agree: usize = parts.iter().map(|x| x.4).sum();
    // golden-section search on [0, CALIBRATION_MAX_POWER] (the objective is concave)
    let phi = 0.5 * (5f64.sqrt() - 1.0);
    let (mut a, mut b) = (0.0, CALIBRATION_MAX_POWER);
    let (mut c, mut d) = (b - phi * (b - a), a + phi * (b - a));
    let (mut fc, mut fd) = (objective(c), objective(d));
    while b - a > 1e-3 {
        if fc > fd {
            (b, d, fd) = (d, c, fc);
            c = b - phi * (b - a);
            fc = objective(c);
        } else {
            (a, c, fc) = (c, d, fd);
            d = a + phi * (b - a);
            fd = objective(d);
        }
    }
    let p = 0.5 * (a + b);
    let n = molecules.max(1) as f64;
    let curve = [0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, p]
        .iter()
        .map(|&q| (q, objective(q) / n))
        .collect();
    Calibration {
        p,
        molecules,
        positional_agreement: agree as f64 / n,
        curve,
        rounds: 0,
    }
}

/// Each molecule's positional likelihood on each of its candidate transcripts (the mean over
/// its intervals, untempered) under the parameters of `g`, keyed by transcript id; a candidate
/// without an interval is absent. For diagnostics.
pub fn molecule_positional_likelihoods(
    molecules: &[ResolvedMolecule],
    interval_mode: TruncationIntervalMode,
    g: &GlobalFit,
) -> Vec<Vec<(u32, f64)>> {
    let (cell, target_ids, _) = build(molecules, interval_mode, None);
    let e = Exps::new(&g.par, cell.max_len().max(1));
    (0..cell.molecule_count())
        .map(|m| {
            let mut per: Vec<(u32, f64, f64)> = Vec::new();
            for x in cell.molecule(m) {
                if let Some(geo) = x.geo() {
                    let lik = g.par.terms(&geo, &e).total();
                    let id = target_ids[x.target];
                    match per.iter_mut().find(|(t, _, _)| *t == id) {
                        Some(p) => {
                            p.1 += x.prior;
                            p.2 += x.prior * lik;
                        }
                        None => per.push((id, x.prior, x.prior * lik)),
                    }
                }
            }
            per.into_iter().map(|(t, w, s)| (t, s / w.max(1e-300))).collect()
        })
        .collect()
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
    spec: TruncationModelSpec,
) -> TruncationFit {
    fit_molecules_with(molecules, interval_mode, None, spec, None, 0, None)
}

/// As [`fit_molecules`], or, with `global`, with the positional parameters fixed to those
/// estimated from all cells ([`fit_global`]): the EM then estimates only the abundances, as
/// oarfish's does after estimating the parameters on a read sample ([`PLAIN_FIRST`] plain
/// steps, then SQUAREM), and the reported parameters are the global ones.
#[allow(clippy::too_many_arguments)]
pub fn fit_molecules_with(
    molecules: &mut [ResolvedMolecule],
    interval_mode: TruncationIntervalMode,
    global: Option<&GlobalFit>,
    spec: TruncationModelSpec,
    shrink: Option<f64>,
    group: u32,
    coverage: Option<&[CellTxpCoverage]>,
) -> TruncationFit {
    let factors = PriorFactors {
        coverage,
        junction: global.and_then(|g| g.junction.as_deref()),
    };
    let (cell, target_ids, mixed_interval_molecules) =
        build(molecules, interval_mode, Some(&factors));
    let max_len = cell.max_len();
    let mut par = match global {
        Some(g) => g.par.clone(),
        None => Params::initial(max_len, spec),
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
    if let Some(g) = global.filter(|g| g.alpha.is_some()) {
        return assign_pooled(
            molecules,
            &cell,
            g,
            shrink,
            max_len,
            mixed_interval_molecules,
            group,
            spec.min_share,
            spec.sample_below,
        );
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
                .map(|a| a.weight(&prev, &par, &exps, cell.positional).0)
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
        prune_shares(&mut probabilities, spec.min_share);
        sample_small_shares(&mut probabilities, spec.sample_below, molecule_uniform(&molecule.vertices));
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

/// Drops the shares of one molecule's transcript probabilities below `min_share` and renormalizes
/// the rest; the largest share is always kept. A molecule is one count, so splitting it into
/// many small shares spreads counts over transcripts the cell most likely does not express.
fn prune_shares(probabilities: &mut [f64], min_share: f64) {
    if min_share <= 0.0 {
        return;
    }
    let max = probabilities.iter().copied().fold(0.0, f64::max);
    let mut kept = 0.0;
    for p in probabilities.iter_mut() {
        if *p < min_share && *p < max {
            *p = 0.0;
        }
        kept += *p;
    }
    if kept > 0.0 {
        for p in probabilities.iter_mut() {
            *p /= kept;
        }
    }
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// A number in [0, 1) fixed by the molecule's (equivalence class, UMI) vertices, so that a draw
/// does not depend on thread scheduling and a rerun gives the same counts.
fn molecule_uniform(vertices: &[(u32, u64)]) -> f64 {
    let mut h: u64 = 0x2545_F491_4F6C_DD1D;
    for &(e, u) in vertices {
        h = splitmix64(h ^ u64::from(e) ^ u.rotate_left(17));
    }
    (splitmix64(h) >> 11) as f64 / (1u64 << 53) as f64
}

/// Gives the combined mass of a molecule's shares below `below` to one of them, drawn with
/// probability proportional to its share (`u` uniform in [0, 1)); `below >= 1` draws one
/// transcript for the whole molecule. Each transcript's expected share is unchanged, so sums
/// over many cells keep the soft counts in expectation, while each cell's counts stop spreading
/// small fractions of a molecule over many transcripts.
fn sample_small_shares(probabilities: &mut [f64], below: f64, u: f64) {
    if below <= 0.0 {
        return;
    }
    let small = |p: f64| p > 0.0 && (below >= 1.0 || p < below);
    let mass: f64 = probabilities.iter().copied().filter(|&p| small(p)).sum();
    if mass <= 0.0 {
        return;
    }
    let target = u * mass;
    let (mut acc, mut picked, mut last, mut whole) = (0.0, None, 0, true);
    for (i, p) in probabilities.iter_mut().enumerate() {
        if small(*p) {
            acc += *p;
            if picked.is_none() && target < acc {
                picked = Some(i);
            }
            last = i;
            *p = 0.0;
        } else if *p > 0.0 {
            whole = false;
        }
    }
    // a molecule drawn whole is exactly one count, not the rounded sum of its shares
    probabilities[picked.unwrap_or(last)] = if whole { 1.0 } else { mass };
}

/// Floor (molecules) added to every pooled abundance a cell uses, so that a transcript absent
/// from the cells the pooled fit saw (with `--truncation-global-cells`) keeps some weight.
const POOLED_FLOOR: f64 = 1e-3;

/// A cell's molecules assigned under the pooled fit of [`fit_pooled`]: with `shrink` `None`
/// ([`TruncationAbundance::Global`]) each molecule takes its posterior under the pooled
/// abundances; with `shrink` `Some(k)` ([`TruncationAbundance::Shrink`]) under the cell's own
/// MAP-EM abundances, with Dirichlet pseudo-counts `k * n(locus) * theta(t)` on each candidate
/// `t`, where `theta(t)` is `t`'s pooled share of its locus and `n(locus)` the cell's molecules
/// in that locus: `k = 0` is the cell's own maximum-likelihood estimate, a large `k` the pooled
/// isoform fractions. When the fit carries an estimated prior strength (`--truncation-shrink
/// auto`, [`GlobalFit::prior`]) the pseudo-counts are `alpha * theta(t)` instead, the same
/// number of pseudo-molecules in every cell, so that a cell with many molecules of a locus
/// keeps its own isoform usage and one with few takes the pooled shares; an infinite `alpha` is
/// the pooled assignment itself.
fn assign_pooled(
    molecules: &mut [ResolvedMolecule],
    cell: &CellMolecules,
    g: &GlobalFit,
    shrink: Option<f64>,
    max_len: u32,
    mixed_interval_molecules: usize,
    group: u32,
    min_share: f64,
    sample_below: f64,
) -> TruncationFit {
    let alphas = g.alpha.as_ref().expect("pooled fit has abundances");
    let group = (group as usize).min(alphas.len() - 1);
    let alpha = &alphas[group];
    let fx = Fixed::new(cell, &g.par, &Exps::new(&g.par, max_len));
    let targets = &cell.target_ids;
    let pooled: Vec<f64> = targets
        .iter()
        .map(|&t| alpha[t as usize] + POOLED_FLOOR)
        .collect();
    let mut ab = pooled.clone();
    let (mut steps, mut converged) = (0_usize, true);
    let auto = g.prior.as_ref().map(|p| p.alpha);
    let locus_alpha = g.prior.as_ref().and_then(|p| p.per_locus.as_deref());
    // an infinite common strength with no per-locus strengths is the pooled assignment itself
    let pooled_only = auto.is_some_and(f64::is_infinite) && locus_alpha.is_none();
    if let Some(k) = shrink.filter(|_| !pooled_only) {
        let locus = g.locus.as_ref().expect("pooled fit has loci");
        let mass = &g.locus_mass()[group];
        // each target's pooled share of its locus, and its locus as an index over the cell's loci
        let theta: Vec<f64> = targets
            .iter()
            .zip(&pooled)
            .map(|(&t, a)| a / mass[locus[t as usize] as usize].max(1e-300))
            .collect();
        let mut local = HashMap::new();
        let lidx: Vec<usize> = targets
            .iter()
            .map(|&t| {
                let n = local.len();
                *local.entry(locus[t as usize]).or_insert(n)
            })
            .collect();
        let mut cnt = vec![0.0; targets.len()];
        let mut per_locus = vec![0.0; local.len()];
        converged = false;
        while steps < MAX_EM_ITER as usize {
            step_static(&fx, &ab, &mut cnt);
            steps += 1;
            per_locus.fill(0.0);
            for (c, &l) in cnt.iter().zip(&lidx) {
                per_locus[l] += c;
            }
            let mut rel = 0.0_f64;
            for i in 0..ab.len() {
                let pseudo = match (locus_alpha, auto) {
                    (Some(la), _) => la[locus[targets[i] as usize] as usize],
                    (None, Some(alpha)) => alpha,
                    (None, None) => k * per_locus[lidx[i]],
                };
                // a locus with no variation between cells keeps the pooled shares
                let next = if pseudo.is_infinite() {
                    pooled[i]
                } else {
                    cnt[i] + pseudo * theta[i]
                };
                if ab[i] > MIN_READ_THRESH {
                    rel = rel.max((next - ab[i]).abs() / ab[i]);
                }
                ab[i] = next;
            }
            if rel < CONVERGENCE_THRESH {
                converged = true;
                break;
            }
        }
    }
    for (m, molecule) in molecules.iter_mut().enumerate() {
        let (lo, hi) = (fx.off[m], fx.off[m + 1]);
        let mut weights: Vec<f64> = (lo..hi).map(|i| ab[fx.t[i]] * fx.w[i]).collect();
        let mut denominator: f64 = weights.iter().sum();
        if !(denominator.is_finite() && denominator > EM_DENOM_THRESH) {
            weights = cell.molecule(m).iter().map(|a| a.prior).collect();
            denominator = weights.iter().sum();
        }
        let mut probabilities = vec![0.0; targets.len()];
        for (i, w) in (lo..hi).zip(weights) {
            probabilities[fx.t[i]] += w / denominator;
        }
        prune_shares(&mut probabilities, min_share);
        sample_small_shares(&mut probabilities, sample_below, molecule_uniform(&molecule.vertices));
        molecule.transcript_probs = probabilities
            .into_iter()
            .enumerate()
            .filter_map(|(t, p)| (p > 0.0).then_some((targets[t], p)))
            .collect();
    }
    TruncationFit {
        params: g.params,
        iterations: steps,
        frozen_at: Some(0),
        freeze_reason: Some(FreezeReason::Global),
        converged,
        molecule_count: molecules.len(),
        has_positional_evidence: cell.positional && cell.alns.iter().any(|a| a.span.is_some()),
        mixed_interval_molecules,
    }
}

impl GlobalFit {
    /// Total pooled abundance of each locus (indexed by the locus id, a transcript id), per group.
    fn locus_mass(&self) -> &[Vec<f64>] {
        self.locus_mass_cache.get_or_init(|| {
            let locus = self.locus.as_ref().expect("pooled fit has loci");
            self.alpha
                .as_ref()
                .expect("pooled fit has abundances")
                .iter()
                .map(|alpha| {
                    let mut mass = vec![0.0; alpha.len()];
                    for (a, &l) in alpha.iter().zip(locus) {
                        mass[l as usize] += a + POOLED_FLOOR;
                    }
                    mass
                })
                .collect()
        })
    }

    /// Number of groups with their own pooled abundances (1 without `--truncation-groups`).
    pub fn n_groups(&self) -> usize {
        self.alpha.as_ref().map_or(1, Vec::len)
    }
}

impl GlobalFit {
    /// The pooled abundances of the first group (all cells without `--truncation-groups`),
    /// indexed by transcript id; empty without pooled abundances.
    pub fn pooled_abundances(&self) -> &[f64] {
        self.alpha
            .as_ref()
            .and_then(|a| a.first())
            .map_or(&[], Vec::as_slice)
    }

    /// This fit with `alpha` as its per-group abundances in place of the pooled ones
    /// (`--truncation-oracle`): every cell's molecules are then assigned under its group's
    /// vector exactly as under pooled abundances, with no prior strength.
    pub fn with_abundances(mut self, alpha: Vec<Vec<f64>>) -> Self {
        self.alpha = Some(alpha);
        self.prior = None;
        self.locus_mass_cache = std::sync::OnceLock::new();
        self
    }
}

/// Fixed per-cell abundances read for `--truncation-oracle` ([`load_oracle`]).
pub struct OracleAbundances {
    /// One group per barcode the file lists, in order of first appearance, named after it;
    /// unlisted barcodes fall in the extra group as usual.
    pub groups: TruncationGroups,
    /// The abundance vector of each group, indexed by transcript id, and last the pooled
    /// abundances for the unlisted barcodes.
    pub alpha: Vec<Vec<f64>>,
    /// Distinct transcripts of the reference the file named.
    pub transcripts: usize,
    /// Lines whose transcript is not in the reference.
    pub unmatched: usize,
}

/// A transcript's name without a version suffix (`.N`).
fn strip_version(s: &str) -> &str {
    match s.rsplit_once('.') {
        Some((stem, v)) if !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()) => stem,
        _ => s,
    }
}

/// Reads `transcript<TAB>barcode<TAB>count` lines (a header line, whose count does not parse,
/// is skipped) into one abundance vector per barcode, transcripts matched by name as written or
/// without a version suffix; the barcodes are written as alevin-fry writes them. The vectors are
/// dense: the file's cells times the reference's transcripts. For evaluation on simulated data,
/// where the counts are the truth: the assignment under them bounds what any choice of
/// abundances can achieve with the same positional model.
pub fn load_oracle(
    path: &std::path::Path,
    ref_names: &[String],
    pooled: &[f64],
) -> std::io::Result<OracleAbundances> {
    use std::io::BufRead;
    let mut exact: HashMap<&str, u32> = HashMap::new();
    let mut stripped: HashMap<&str, u32> = HashMap::new();
    for (i, n) in ref_names.iter().enumerate() {
        exact.entry(n.as_str()).or_insert(i as u32);
        stripped.entry(strip_version(n)).or_insert(i as u32);
    }
    let n_refs = ref_names.len();
    let mut groups = TruncationGroups::default();
    let mut alpha: Vec<Vec<f64>> = Vec::new();
    let mut seen = vec![false; n_refs];
    let mut unmatched = 0_usize;
    for line in std::io::BufReader::new(std::fs::File::open(path)?).lines() {
        let line = line?;
        let mut f = line.split('\t');
        let (Some(t), Some(bc), Some(c)) = (f.next(), f.next(), f.next()) else {
            continue;
        };
        let Ok(count) = c.trim().parse::<f64>() else {
            continue;
        };
        let t = t.trim();
        let Some(&tid) = exact.get(t).or_else(|| stripped.get(strip_version(t))) else {
            unmatched += 1;
            continue;
        };
        let gi = groups.own_group(bc.trim()) as usize;
        if gi == alpha.len() {
            alpha.push(vec![0.0; n_refs]);
        }
        alpha[gi][tid as usize] += count;
        seen[tid as usize] = true;
    }
    let mut last = pooled.to_vec();
    last.resize(n_refs, 0.0);
    alpha.push(last);
    Ok(OracleAbundances {
        groups,
        alpha,
        transcripts: seen.iter().filter(|&&s| s).count(),
        unmatched,
    })
}

/// Writes, for `--truncation-dump`, every cell's molecule counts per candidate transcript:
/// gzip-compressed `barcode<TAB>transcript<TAB>unique<TAB>candidates` lines, `unique` the
/// molecules whose only candidate is the transcript and `candidates` the molecules that have it
/// among their candidates. Returns the lines written and the one-candidate molecules counted.
pub fn write_unique_counts(
    cells: &[CellMolecules],
    barcodes: &[String],
    ref_names: &[String],
    path: &std::path::Path,
) -> std::io::Result<(usize, usize)> {
    use std::io::Write;
    let mut out = std::io::BufWriter::new(flate2::write::GzEncoder::new(
        std::fs::File::create(path)?,
        flate2::Compression::fast(),
    ));
    writeln!(out, "barcode\ttranscript\tunique\tcandidates")?;
    let (mut lines, mut molecules) = (0_usize, 0_usize);
    let (mut unique, mut any) = (Vec::new(), Vec::new());
    let mut targets = Vec::new();
    for (cell, bc) in cells.iter().zip(barcodes) {
        unique.clear();
        unique.resize(cell.n_targets, 0_u32);
        any.clear();
        any.resize(cell.n_targets, 0_u32);
        for m in 0..cell.molecule_count() {
            let mol = cell.molecule(m);
            targets.clear();
            targets.extend(mol.iter().map(|a| a.target));
            targets.sort_unstable();
            targets.dedup();
            if targets.len() == 1 {
                unique[targets[0]] += 1;
                molecules += 1;
            }
            for &t in &targets {
                any[t] += 1;
            }
        }
        for i in 0..cell.n_targets {
            if any[i] > 0 {
                writeln!(
                    out,
                    "{bc}\t{}\t{}\t{}",
                    ref_names[cell.target_ids[i] as usize],
                    unique[i],
                    any[i]
                )?;
                lines += 1;
            }
        }
    }
    out.flush()?;
    out.into_inner().map_err(|e| e.into_error())?.finish()?;
    Ok((lines, molecules))
}
