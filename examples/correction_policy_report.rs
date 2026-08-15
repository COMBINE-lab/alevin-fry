/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Compare Hamming-one and historical shift-aware barcode correction.
//!
//! This is a developer measurement utility, not a supported user command:
//!
//! ```text
//! cargo run --release --example correction_policy_report -- \
//!   <all_freq.bin> <permit_freq.bin>
//! ```

use alevin_fry::barcode_correction::{
    BarcodeNeighborhood, BarcodeResolution, Confidence, CorrectionDecision, CorrectionIndex,
    CorrectionSpec, RetainedSource,
};
use anyhow::{Context, bail};
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Debug, Default, Serialize)]
struct WeightedCount {
    distinct: u64,
    reads: u64,
}

impl WeightedCount {
    fn add(&mut self, reads: u64) {
        self.distinct += 1;
        self.reads += reads;
    }
}

#[derive(Debug, Default, Serialize)]
struct Comparison {
    identical_target: WeightedCount,
    identical_rejection: WeightedCount,
    shift_only_rescue: WeightedCount,
    extra_shift_induced_rejection: WeightedCount,
    changed_target: WeightedCount,
    hamming_ambiguous_shift_not_found: WeightedCount,
    hamming_not_found_shift_ambiguous: WeightedCount,
}

#[derive(Debug, Serialize)]
struct TargetDelta {
    targets_with_changed_frequency: u64,
    targets_gaining_reads: u64,
    targets_losing_reads: u64,
    total_absolute_read_delta: u128,
    maximum_absolute_read_delta: u64,
    top_absolute_deltas: Vec<TargetDeltaEntry>,
}

#[derive(Debug, Serialize)]
struct TargetDeltaEntry {
    target: u64,
    hamming_reads: u64,
    shift_reads: u64,
    signed_delta: i128,
}

#[derive(Debug, Serialize)]
struct PolicyReport {
    comparison: Comparison,
    hamming_accepted: WeightedCount,
    shift_accepted: WeightedCount,
    target_frequency_delta: TargetDelta,
    elapsed_seconds: f64,
}

#[derive(Debug, Serialize)]
struct Report {
    all_freq: PathBuf,
    permit_freq: PathBuf,
    barcode_length: u8,
    observed_distinct: u64,
    observed_reads: u64,
    retained_targets: u64,
    frequency_confidence: String,
    unique: PolicyReport,
    frequency: PolicyReport,
}

fn read_frequency_map(path: &Path) -> anyhow::Result<(u64, u8, HashMap<u64, u64>)> {
    let file = File::open(path)
        .with_context(|| format!("could not open frequency file {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut header = [0_u8; 16];
    reader
        .read_exact(&mut header)
        .with_context(|| format!("could not read frequency header from {}", path.display()))?;
    let version = u64::from_le_bytes(header[..8].try_into().unwrap());
    let barcode_len_u64 = u64::from_le_bytes(header[8..].try_into().unwrap());
    let barcode_len = u8::try_from(barcode_len_u64)
        .with_context(|| format!("barcode length {barcode_len_u64} does not fit in u8"))?;
    if !(1..=32).contains(&barcode_len) {
        bail!(
            "frequency file {} declares invalid barcode length {barcode_len}",
            path.display()
        );
    }
    let counts = bincode::deserialize_from(reader)
        .with_context(|| format!("could not decode frequency map from {}", path.display()))?;
    Ok((version, barcode_len, counts))
}

fn compare_policies(
    observed: &HashMap<u64, u64>,
    retained: &HashMap<u64, u64>,
    barcode_len: u8,
    resolution: BarcodeResolution,
) -> anyhow::Result<PolicyReport> {
    let started = Instant::now();
    let retained_sources = || {
        retained.keys().map(|&barcode| RetainedSource {
            source: barcode,
            target: barcode,
            exact_count: observed.get(&barcode).copied().unwrap_or(0),
        })
    };
    let hamming = CorrectionIndex::new(
        CorrectionSpec {
            barcode_len,
            neighborhood: BarcodeNeighborhood::HammingOne,
            resolution,
        },
        retained_sources(),
    )?;
    let shift = CorrectionIndex::new(
        CorrectionSpec {
            barcode_len,
            neighborhood: BarcodeNeighborhood::SubstitutionOrShiftOne,
            resolution,
        },
        retained_sources(),
    )?;

    let mut comparison = Comparison::default();
    let mut hamming_accepted = WeightedCount::default();
    let mut shift_accepted = WeightedCount::default();
    let mut hamming_targets = BTreeMap::<u64, u64>::new();
    let mut shift_targets = BTreeMap::<u64, u64>::new();

    for (&observed_barcode, &reads) in observed {
        let hamming_decision = hamming.resolve(observed_barcode);
        let shift_decision = shift.resolve(observed_barcode);
        let hamming_target = hamming_decision.target();
        let shift_target = shift_decision.target();

        if let Some(target) = hamming_target {
            hamming_accepted.add(reads);
            *hamming_targets.entry(target).or_default() += reads;
        }
        if let Some(target) = shift_target {
            shift_accepted.add(reads);
            *shift_targets.entry(target).or_default() += reads;
        }

        match (hamming_target, shift_target) {
            (Some(lhs), Some(rhs)) if lhs == rhs => comparison.identical_target.add(reads),
            (Some(_), Some(_)) => comparison.changed_target.add(reads),
            (None, Some(_)) => comparison.shift_only_rescue.add(reads),
            (Some(_), None) => comparison.extra_shift_induced_rejection.add(reads),
            (None, None) => comparison.identical_rejection.add(reads),
        }
        match (hamming_decision, shift_decision) {
            (CorrectionDecision::Ambiguous, CorrectionDecision::NotFound) => {
                comparison.hamming_ambiguous_shift_not_found.add(reads);
            }
            (CorrectionDecision::NotFound, CorrectionDecision::Ambiguous) => {
                comparison.hamming_not_found_shift_ambiguous.add(reads);
            }
            _ => {}
        }
    }

    let mut all_targets = hamming_targets.keys().copied().collect::<Vec<_>>();
    all_targets.extend(shift_targets.keys().copied());
    all_targets.sort_unstable();
    all_targets.dedup();
    let mut deltas = Vec::new();
    for target in all_targets {
        let hamming_reads = hamming_targets.get(&target).copied().unwrap_or(0);
        let shift_reads = shift_targets.get(&target).copied().unwrap_or(0);
        if hamming_reads != shift_reads {
            deltas.push(TargetDeltaEntry {
                target,
                hamming_reads,
                shift_reads,
                signed_delta: i128::from(shift_reads) - i128::from(hamming_reads),
            });
        }
    }
    deltas.sort_unstable_by(|lhs, rhs| {
        rhs.signed_delta
            .unsigned_abs()
            .cmp(&lhs.signed_delta.unsigned_abs())
            .then_with(|| lhs.target.cmp(&rhs.target))
    });
    let target_frequency_delta = TargetDelta {
        targets_with_changed_frequency: deltas.len() as u64,
        targets_gaining_reads: deltas.iter().filter(|delta| delta.signed_delta > 0).count() as u64,
        targets_losing_reads: deltas.iter().filter(|delta| delta.signed_delta < 0).count() as u64,
        total_absolute_read_delta: deltas
            .iter()
            .map(|delta| delta.signed_delta.unsigned_abs())
            .sum(),
        maximum_absolute_read_delta: deltas
            .first()
            .map_or(0, |delta| delta.signed_delta.unsigned_abs() as u64),
        top_absolute_deltas: deltas.into_iter().take(10).collect(),
    };

    Ok(PolicyReport {
        comparison,
        hamming_accepted,
        shift_accepted,
        target_frequency_delta,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    })
}

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args_os().skip(1);
    let all_freq = PathBuf::from(args.next().context("missing all_freq.bin argument")?);
    let permit_freq = PathBuf::from(args.next().context("missing permit_freq.bin argument")?);
    if args.next().is_some() {
        bail!("usage: correction_policy_report <all_freq.bin> <permit_freq.bin>");
    }

    let (all_version, all_len, observed) = read_frequency_map(&all_freq)?;
    let (permit_version, permit_len, retained) = read_frequency_map(&permit_freq)?;
    if all_version != permit_version {
        bail!(
            "frequency-file versions differ: all_freq={all_version}, permit_freq={permit_version}"
        );
    }
    if all_len != permit_len {
        bail!("barcode lengths differ: all_freq={all_len}, permit_freq={permit_len}");
    }
    let confidence = Confidence::RNA;
    let unique = compare_policies(&observed, &retained, all_len, BarcodeResolution::Unique)?;
    let frequency = compare_policies(
        &observed,
        &retained,
        all_len,
        BarcodeResolution::Frequency {
            confidence,
            pseudocount: 1,
        },
    )?;
    let report = Report {
        all_freq,
        permit_freq,
        barcode_length: all_len,
        observed_distinct: observed.len() as u64,
        observed_reads: observed
            .values()
            .map(|&count| u128::from(count))
            .sum::<u128>() as u64,
        retained_targets: retained.len() as u64,
        frequency_confidence: confidence.to_string(),
        unique,
        frequency,
    };
    serde_json::to_writer_pretty(std::io::stdout().lock(), &report)?;
    println!();
    Ok(())
}
