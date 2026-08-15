/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Compare two compiled multi-barcode correction plans against a RAD file.
//!
//! This is a developer measurement utility, not a supported user command:
//!
//! ```text
//! cargo run --release --example multi_correction_plan_report -- \
//!   <map.rad> <hamming-output-dir> <shift-output-dir> [threads]
//! ```

use ahash::AHashMap;
use alevin_fry::correction_plan::{BarcodeCorrection, CORRECTION_PLAN_FILENAME, CorrectionPlan};
use anyhow::{Context, bail};
use libradicl::header::RadPrelude;
use libradicl::record::{
    CollatableMappedRecord, HierarchicallyCollatable, MultiBarcodeReadRecord,
    MultiBarcodeRecordContext,
};
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufReader, Read};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

const IDENTICAL_TARGET: u8 = 0;
const SHIFT_ONLY_RESCUE: u8 = 1;
const EXTRA_SHIFT_REJECTION: u8 = 2;
const CHANGED_TARGET: u8 = 3;

#[derive(Debug, Default, Clone, Copy, Serialize)]
struct OutcomeCounts {
    identical_target: u64,
    shift_only_rescue: u64,
    extra_shift_induced_rejection: u64,
    changed_target: u64,
    identical_rejection: u64,
}

#[derive(Debug, Serialize)]
struct TargetFrequencyDelta {
    targets_with_changed_frequency: u64,
    targets_gaining_reads: u64,
    targets_losing_reads: u64,
    total_absolute_read_delta: u128,
    maximum_absolute_read_delta: u64,
    top_absolute_deltas: Vec<TargetFrequencyDeltaEntry>,
}

#[derive(Debug, Serialize)]
struct TargetFrequencyDeltaEntry {
    sample: u64,
    cell: u64,
    hamming_reads: u64,
    shift_reads: u64,
    signed_delta: i128,
}

impl OutcomeCounts {
    fn record(&mut self, outcome: Option<u8>, count: u64) {
        match outcome {
            Some(IDENTICAL_TARGET) => self.identical_target += count,
            Some(SHIFT_ONLY_RESCUE) => self.shift_only_rescue += count,
            Some(EXTRA_SHIFT_REJECTION) => self.extra_shift_induced_rejection += count,
            Some(CHANGED_TARGET) => self.changed_target += count,
            None => self.identical_rejection += count,
            Some(_) => unreachable!("invalid correction outcome"),
        }
    }
}

#[derive(Debug, Serialize)]
struct Report {
    rad: PathBuf,
    hamming_plan: PathBuf,
    shift_plan: PathBuf,
    total_records: u64,
    identically_sample_routed_records: u64,
    sample_routing_difference_records: u64,
    distinct_accepted_outcomes: OutcomeCounts,
    read_weighted_outcomes: OutcomeCounts,
    target_frequency_delta: TargetFrequencyDelta,
    classified_observed_pairs: u64,
    elapsed_seconds: f64,
}

fn plan_path(directory: &Path) -> PathBuf {
    directory.join(CORRECTION_PLAN_FILENAME)
}

fn sample_map(plan: &CorrectionPlan) -> AHashMap<u64, u64> {
    plan.sample_corrections
        .iter()
        .map(|entry| (entry.observed, entry.corrected))
        .collect()
}

fn read_frequency_map(path: &Path) -> anyhow::Result<HashMap<u64, u64>> {
    let file = File::open(path)
        .with_context(|| format!("could not open frequency file {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut header = [0_u8; 16];
    reader
        .read_exact(&mut header)
        .with_context(|| format!("could not read frequency header from {}", path.display()))?;
    bincode::deserialize_from(reader)
        .with_context(|| format!("could not decode frequency map from {}", path.display()))
}

fn sample_names(directory: &Path) -> anyhow::Result<BTreeMap<u64, String>> {
    let file = File::open(directory.join("sample_info.json"))?;
    let document: serde_json::Value = serde_json::from_reader(BufReader::new(file))?;
    document["samples"]
        .as_array()
        .context("sample_info.json has no samples array")?
        .iter()
        .map(|sample| {
            let barcode = sample["barcode"]
                .as_str()
                .context("sample has no barcode")?
                .strip_prefix("0x")
                .context("sample barcode does not have a 0x prefix")?;
            Ok((
                u64::from_str_radix(barcode, 16)?,
                sample["name"]
                    .as_str()
                    .context("sample has no name")?
                    .to_owned(),
            ))
        })
        .collect()
}

fn target_frequency_delta(
    hamming_dir: &Path,
    shift_dir: &Path,
) -> anyhow::Result<TargetFrequencyDelta> {
    let hamming_names = sample_names(hamming_dir)?;
    let shift_names = sample_names(shift_dir)?;
    if hamming_names != shift_names {
        bail!("the two sample_info.json files declare different samples");
    }

    let mut deltas = Vec::new();
    for (sample, name) in hamming_names {
        let hamming_path = hamming_dir
            .join(format!("sample_{name}"))
            .join("permit_freq.bin");
        let shift_path = shift_dir
            .join(format!("sample_{name}"))
            .join("permit_freq.bin");
        let hamming = if hamming_path.exists() {
            read_frequency_map(&hamming_path)?
        } else {
            HashMap::new()
        };
        let shift = if shift_path.exists() {
            read_frequency_map(&shift_path)?
        } else {
            HashMap::new()
        };
        let mut cells = hamming.keys().copied().collect::<Vec<_>>();
        cells.extend(shift.keys().copied());
        cells.sort_unstable();
        cells.dedup();
        for cell in cells {
            let hamming_reads = hamming.get(&cell).copied().unwrap_or(0);
            let shift_reads = shift.get(&cell).copied().unwrap_or(0);
            if hamming_reads != shift_reads {
                deltas.push(TargetFrequencyDeltaEntry {
                    sample,
                    cell,
                    hamming_reads,
                    shift_reads,
                    signed_delta: i128::from(shift_reads) - i128::from(hamming_reads),
                });
            }
        }
    }
    deltas.sort_unstable_by(|lhs, rhs| {
        rhs.signed_delta
            .unsigned_abs()
            .cmp(&lhs.signed_delta.unsigned_abs())
            .then_with(|| lhs.sample.cmp(&rhs.sample))
            .then_with(|| lhs.cell.cmp(&rhs.cell))
    });
    Ok(TargetFrequencyDelta {
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
    })
}

fn merge_scope(
    sample: u64,
    cell_bits: u32,
    hamming: &[BarcodeCorrection],
    shift: &[BarcodeCorrection],
    outcomes: &mut AHashMap<u64, u8>,
    distinct: &mut OutcomeCounts,
) {
    let mut h = 0;
    let mut s = 0;
    while h < hamming.len() || s < shift.len() {
        let (observed, outcome) = match (hamming.get(h), shift.get(s)) {
            (Some(lhs), Some(rhs)) if lhs.observed == rhs.observed => {
                let outcome = if lhs.corrected == rhs.corrected {
                    IDENTICAL_TARGET
                } else {
                    CHANGED_TARGET
                };
                h += 1;
                s += 1;
                (lhs.observed, outcome)
            }
            (Some(lhs), Some(rhs)) if lhs.observed < rhs.observed => {
                h += 1;
                (lhs.observed, EXTRA_SHIFT_REJECTION)
            }
            (Some(_), Some(rhs)) => {
                s += 1;
                (rhs.observed, SHIFT_ONLY_RESCUE)
            }
            (Some(lhs), None) => {
                h += 1;
                (lhs.observed, EXTRA_SHIFT_REJECTION)
            }
            (None, Some(rhs)) => {
                s += 1;
                (rhs.observed, SHIFT_ONLY_RESCUE)
            }
            (None, None) => unreachable!(),
        };
        let key = (sample << cell_bits) | observed;
        let previous = outcomes.insert(key, outcome);
        assert!(previous.is_none(), "duplicate sample/cell correction key");
        distinct.record(Some(outcome), 1);
    }
}

fn build_outcome_index(
    hamming: &CorrectionPlan,
    shift: &CorrectionPlan,
) -> anyhow::Result<(AHashMap<u64, u8>, OutcomeCounts)> {
    if hamming.cell_barcode_len != shift.cell_barcode_len
        || hamming.sample_barcode_len != shift.sample_barcode_len
    {
        bail!("correction plans declare different barcode lengths");
    }
    let sample_len = hamming
        .sample_barcode_len
        .context("the first correction plan is not multi-barcode")?;
    let cell_bits = 2 * u32::from(hamming.cell_barcode_len);
    if 2 * u32::from(sample_len) + cell_bits > 64 {
        bail!("combined sample/cell barcode does not fit in this report's packed u64 key");
    }

    let hamming_scopes: BTreeMap<_, _> = hamming
        .cell_scopes
        .iter()
        .map(|scope| {
            Ok((
                scope
                    .sample_barcode
                    .context("ordinary scope in multi-barcode correction plan")?,
                scope.corrections.as_slice(),
            ))
        })
        .collect::<anyhow::Result<_>>()?;
    let shift_scopes: BTreeMap<_, _> = shift
        .cell_scopes
        .iter()
        .map(|scope| {
            Ok((
                scope
                    .sample_barcode
                    .context("ordinary scope in multi-barcode correction plan")?,
                scope.corrections.as_slice(),
            ))
        })
        .collect::<anyhow::Result<_>>()?;
    if hamming_scopes.keys().ne(shift_scopes.keys()) {
        bail!("correction plans have different canonical sample scopes");
    }

    let capacity = hamming_scopes
        .values()
        .map(|scope| scope.len())
        .sum::<usize>()
        .saturating_add(shift_scopes.values().map(|scope| scope.len()).sum());
    let mut outcomes = AHashMap::with_capacity(capacity);
    let mut distinct = OutcomeCounts::default();
    for (&sample, hamming_scope) in &hamming_scopes {
        merge_scope(
            sample,
            cell_bits,
            hamming_scope,
            shift_scopes[&sample],
            &mut outcomes,
            &mut distinct,
        );
    }
    Ok((outcomes, distinct))
}

fn main() -> anyhow::Result<()> {
    let started = Instant::now();
    let mut args = std::env::args_os().skip(1);
    let rad = PathBuf::from(args.next().context("missing map.rad argument")?);
    let hamming_dir = PathBuf::from(args.next().context("missing Hamming output directory")?);
    let shift_dir = PathBuf::from(args.next().context("missing shift output directory")?);
    let threads = args
        .next()
        .map(|value| value.to_string_lossy().parse::<usize>())
        .transpose()?
        .unwrap_or(16)
        .max(1);
    if args.next().is_some() {
        bail!(
            "usage: multi_correction_plan_report <map.rad> <hamming-output-dir> <shift-output-dir> [threads]"
        );
    }

    let hamming_path = plan_path(&hamming_dir);
    let shift_path = plan_path(&shift_dir);
    let hamming = CorrectionPlan::read_from(&hamming_path)?;
    let shift = CorrectionPlan::read_from(&shift_path)?;
    let hamming_samples = Arc::new(sample_map(&hamming));
    let shift_samples = Arc::new(sample_map(&shift));
    let cell_bits = 2 * u32::from(hamming.cell_barcode_len);
    let (outcomes, distinct_outcomes) = build_outcome_index(&hamming, &shift)?;
    let outcomes = Arc::new(outcomes);
    drop(hamming);
    drop(shift);

    let file = File::open(&rad).with_context(|| format!("could not open {}", rad.display()))?;
    let mut reader = BufReader::new(file);
    let prelude = RadPrelude::from_bytes(&mut reader)?;
    prelude.file_tags.parse_tags_from_bytes(&mut reader)?;
    prelude.get_record_context::<MultiBarcodeRecordContext>()?;
    let mut chunk_reader = libradicl::readers::ParallelChunkReader::<MultiBarcodeReadRecord>::new(
        &prelude,
        NonZeroUsize::new(threads).unwrap(),
    );

    let worker_results = std::thread::scope(|scope| -> anyhow::Result<Vec<_>> {
        let mut handles = Vec::with_capacity(threads);
        for _ in 0..threads {
            let chunks = chunk_reader.chunk_iter();
            let outcomes = outcomes.clone();
            let hamming_samples = hamming_samples.clone();
            let shift_samples = shift_samples.clone();
            handles.push(scope.spawn(move || {
                let mut total = 0_u64;
                let mut routed = 0_u64;
                let mut sample_difference = 0_u64;
                let mut counts = OutcomeCounts::default();
                for meta_chunk in chunks {
                    for chunk in meta_chunk.iter() {
                        for read in &chunk.reads {
                            total += 1;
                            let raw_sample = read.collation_key_at_level(0);
                            match (
                                hamming_samples.get(&raw_sample),
                                shift_samples.get(&raw_sample),
                            ) {
                                (Some(&h_sample), Some(&s_sample)) if h_sample == s_sample => {
                                    routed += 1;
                                    let observed_cell = read.collate_key();
                                    let key = (h_sample << cell_bits) | observed_cell;
                                    counts.record(outcomes.get(&key).copied(), 1);
                                }
                                _ => sample_difference += 1,
                            }
                        }
                    }
                }
                (total, routed, sample_difference, counts)
            }));
        }
        chunk_reader.start(&mut reader, None::<fn(u64, u64)>)?;
        handles
            .into_iter()
            .map(|handle| Ok(handle.join().expect("report worker panicked")))
            .collect()
    })?;

    let mut total_records = 0;
    let mut identically_sample_routed_records = 0;
    let mut sample_routing_difference_records = 0;
    let mut read_weighted_outcomes = OutcomeCounts::default();
    for (total, routed, sample_difference, counts) in worker_results {
        total_records += total;
        identically_sample_routed_records += routed;
        sample_routing_difference_records += sample_difference;
        read_weighted_outcomes.identical_target += counts.identical_target;
        read_weighted_outcomes.shift_only_rescue += counts.shift_only_rescue;
        read_weighted_outcomes.extra_shift_induced_rejection +=
            counts.extra_shift_induced_rejection;
        read_weighted_outcomes.changed_target += counts.changed_target;
        read_weighted_outcomes.identical_rejection += counts.identical_rejection;
    }
    let report = Report {
        rad,
        hamming_plan: hamming_path,
        shift_plan: shift_path,
        total_records,
        identically_sample_routed_records,
        sample_routing_difference_records,
        distinct_accepted_outcomes: distinct_outcomes,
        read_weighted_outcomes,
        target_frequency_delta: target_frequency_delta(&hamming_dir, &shift_dir)?,
        classified_observed_pairs: outcomes.len() as u64,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    };
    serde_json::to_writer_pretty(std::io::stdout().lock(), &report)?;
    println!();
    Ok(())
}
