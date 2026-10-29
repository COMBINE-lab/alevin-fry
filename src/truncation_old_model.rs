use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

use serde::Serialize;

use crate::eq_class::merge_intervals;
use crate::pugutils_dp_new::ResolvedMolecule;

const END_WINDOW: f64 = 50.0;
const MAX_ITERATIONS: usize = 100;
const CONVERGENCE_TOLERANCE: f64 = 1e-6;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub enum TruncationIntervalMode {
    #[default]
    All,
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

#[derive(Clone, Copy, Debug)]
pub struct TruncationParams {
    pub full_length_fraction: f64,
    pub stop_rate: f64,
    pub three_prime_fraction: f64,
    pub three_prime_decay: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct TruncationFit {
    pub params: TruncationParams,
    pub iterations: usize,
    pub converged: bool,
    pub molecule_count: usize,
    pub has_positional_evidence: bool,
}

impl TruncationParams {
    fn initial() -> Self {
        Self {
            full_length_fraction: 0.5,
            stop_rate: 1e-3,
            three_prime_fraction: 0.5,
            three_prime_decay: 2e-3,
        }
    }

    fn terms(self, alignment: &Alignment) -> Terms {
        if !alignment.has_position {
            return Terms::neutral();
        }

        let len = alignment.transcript_len as f64;
        let start = alignment.start as f64;
        let ext = len - start;
        let d3 = (len - alignment.end as f64).max(0.0);
        let full = ext > len - END_WINDOW;
        let full_w = if full {
            self.full_length_fraction / END_WINDOW
        } else {
            0.0
        };
        let full_surv = if full {
            (1.0 - self.full_length_fraction) * (-self.stop_rate * len).exp() / END_WINDOW
        } else {
            0.0
        };
        let stopped =
            (1.0 - self.full_length_fraction) * self.stop_rate * (-self.stop_rate * ext).exp();
        let at3 = if d3 < END_WINDOW {
            self.three_prime_fraction / END_WINDOW
        } else {
            0.0
        };
        let off3 = (1.0 - self.three_prime_fraction)
            * self.three_prime_decay
            * (-self.three_prime_decay * d3).exp();

        Terms {
            full_w,
            full_surv,
            stopped,
            at3,
            off3,
        }
    }

    fn updated(self, stats: &SufficientStats) -> Self {
        if stats.n == 0.0 {
            return self;
        }
        let fraction = |value: f64| value.clamp(1e-6, 1.0 - 1e-6);
        Self {
            full_length_fraction: fraction(stats.n_full_w / stats.n),
            stop_rate: if stats.exposure_stop > 0.0 {
                (stats.n_stopped / stats.exposure_stop).clamp(1e-7, 1.0)
            } else {
                self.stop_rate
            },
            three_prime_fraction: fraction(stats.n_at3 / stats.n),
            three_prime_decay: if stats.exposure_three_prime > 0.0 {
                (stats.n_off3 / stats.exposure_three_prime).clamp(1e-7, 1.0)
            } else {
                self.three_prime_decay
            },
        }
    }

    fn movement(self, next: Self) -> f64 {
        let relative = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-12);
        (self.full_length_fraction - next.full_length_fraction)
            .abs()
            .max((self.three_prime_fraction - next.three_prime_fraction).abs())
            .max(relative(self.stop_rate, next.stop_rate))
            .max(relative(self.three_prime_decay, next.three_prime_decay))
    }
}

#[derive(Clone, Copy)]
struct Alignment {
    target: usize,
    prior: f64,
    start: u32,
    end: u32,
    transcript_len: u32,
    has_position: bool,
}

struct Terms {
    full_w: f64,
    full_surv: f64,
    stopped: f64,
    at3: f64,
    off3: f64,
}

impl Terms {
    fn neutral() -> Self {
        Self {
            full_w: 0.0,
            full_surv: 0.0,
            stopped: 1.0,
            at3: 1.0,
            off3: 0.0,
        }
    }

    fn start(&self) -> f64 {
        self.full_w + self.full_surv + self.stopped
    }

    fn end(&self) -> f64 {
        self.at3 + self.off3
    }
}

#[derive(Default)]
struct SufficientStats {
    n: f64,
    n_full_w: f64,
    n_stopped: f64,
    n_at3: f64,
    n_off3: f64,
    exposure_stop: f64,
    exposure_three_prime: f64,
}

impl SufficientStats {
    fn add(&mut self, alignment: &Alignment, terms: &Terms, posterior: f64) {
        if !alignment.has_position {
            return;
        }
        let start = terms.start();
        let end = terms.end();
        // An alignment whose start or end likelihood underflowed to zero carries no information
        // about the parameters; skipping it matches oarfish.
        if !(start > 0.0 && end > 0.0) {
            return;
        }
        let p_full_w = posterior * terms.full_w / start;
        let p_survived = posterior * terms.full_surv / start;
        let p_stopped = posterior * terms.stopped / start;
        let p_at3 = posterior * terms.at3 / end;
        let p_off3 = posterior * terms.off3 / end;

        self.n += posterior;
        self.n_full_w += p_full_w;
        self.n_stopped += p_stopped;
        self.n_at3 += p_at3;
        self.n_off3 += p_off3;
        self.exposure_stop += p_survived * alignment.transcript_len as f64
            + p_stopped * (alignment.transcript_len - alignment.start) as f64;
        self.exposure_three_prime += p_off3 * (alignment.transcript_len - alignment.end) as f64;
    }
}

/// Weight without positional terms: abundance times prior, or the prior alone if the abundance
/// is also unusable.
fn fallback_weight(alignment: &Alignment, abundances: &[f64]) -> f64 {
    let weight = abundances[alignment.target] * alignment.prior;
    if weight.is_finite() && weight > 0.0 {
        weight
    } else {
        alignment.prior
    }
}

/// Re-estimate transcript probabilities for one cell's deduplicated molecules using the
/// 3'-anchored, 5'-truncation model. The returned parameters are specific to this cell.
pub fn fit_molecules(
    molecules: &mut [ResolvedMolecule],
    interval_mode: TruncationIntervalMode,
) -> TruncationFit {
    let mut target_ids = Vec::new();
    let mut target_indices = HashMap::new();
    let mut alignments_by_molecule = Vec::<Vec<Alignment>>::with_capacity(molecules.len());

    for molecule in molecules.iter() {
        let mut alignments = Vec::new();
        for &(target_id, prior) in &molecule.transcript_probs {
            assert!(
                prior.is_finite() && prior > 0.0,
                "truncation model received a non-positive transcript probability"
            );
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
                        alignments.push(Alignment {
                            target,
                            prior: interval_prior,
                            start,
                            end,
                            transcript_len: len,
                            has_position: true,
                        });
                    }
                    continue;
                }
            }
            alignments.push(Alignment {
                target,
                prior,
                start: 0,
                end: 0,
                transcript_len: 0,
                has_position: false,
            });
        }
        alignments_by_molecule.push(alignments);
    }

    if molecules.is_empty() || target_ids.is_empty() {
        return TruncationFit {
            params: TruncationParams::initial(),
            iterations: 0,
            converged: false,
            molecule_count: molecules.len(),
            has_positional_evidence: false,
        };
    }

    let initial_abundance = molecules.len() as f64 / target_ids.len() as f64;
    let mut abundances = vec![initial_abundance; target_ids.len()];
    let mut next_abundances = vec![0.0; target_ids.len()];
    let mut params = TruncationParams::initial();
    let mut weighted_terms = Vec::new();
    let mut iterations = 0;
    let mut converged = false;
    let mut has_positional_evidence = false;

    for iteration in 0..MAX_ITERATIONS {
        next_abundances.fill(0.0);
        let mut stats = SufficientStats::default();
        for alignments in &alignments_by_molecule {
            weighted_terms.clear();
            let mut denominator = 0.0;
            for alignment in alignments {
                let terms = params.terms(alignment);
                let weight =
                    abundances[alignment.target] * alignment.prior * terms.start() * terms.end();
                denominator += weight;
                weighted_terms.push((weight, terms));
            }
            let positional = denominator.is_finite() && denominator > 0.0;
            if !positional {
                // The molecule is already deduplicated and must keep its counts, so it falls
                // back to abundance and prior without positional terms or parameter statistics.
                denominator = 0.0;
                for (alignment, entry) in alignments.iter().zip(weighted_terms.iter_mut()) {
                    entry.0 = fallback_weight(alignment, &abundances);
                    denominator += entry.0;
                }
            }
            for (alignment, (weight, terms)) in alignments.iter().zip(&weighted_terms) {
                let posterior = *weight / denominator;
                next_abundances[alignment.target] += posterior;
                if positional {
                    stats.add(alignment, terms, posterior);
                }
            }
        }

        has_positional_evidence |= stats.n > 0.0;
        let next_params = params.updated(&stats);
        let mut abundance_change = 0.0_f64;
        for (previous, next) in abundances.iter().zip(&next_abundances) {
            abundance_change =
                abundance_change.max((next - previous).abs() / previous.abs().max(1e-12));
        }
        let params_change = params.movement(next_params);
        std::mem::swap(&mut abundances, &mut next_abundances);
        params = next_params;
        iterations = iteration + 1;

        if iteration >= 9
            && abundance_change < CONVERGENCE_TOLERANCE
            && params_change < CONVERGENCE_TOLERANCE
        {
            converged = true;
            break;
        }
    }

    for (molecule, alignments) in molecules.iter_mut().zip(&alignments_by_molecule) {
        let mut weights: Vec<_> = alignments
            .iter()
            .map(|alignment| {
                let terms = params.terms(alignment);
                abundances[alignment.target] * alignment.prior * terms.start() * terms.end()
            })
            .collect();
        let mut denominator: f64 = weights.iter().sum();
        if !(denominator.is_finite() && denominator > 0.0) {
            weights = alignments
                .iter()
                .map(|alignment| fallback_weight(alignment, &abundances))
                .collect();
            denominator = weights.iter().sum();
        }
        assert!(
            denominator.is_finite() && denominator > 0.0,
            "truncation model produced invalid final molecule probabilities"
        );
        let mut probabilities = vec![0.0; target_ids.len()];
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
        params,
        iterations,
        converged,
        molecule_count: molecules.len(),
        has_positional_evidence,
    }
}
