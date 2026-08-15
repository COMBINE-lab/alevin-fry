/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use crate::barcode_correction::{
    BarcodeResolution, CorrectionDecision, CorrectionIndex, CorrectionSpec, CorrectionStats,
    RetainedSource,
};
use crate::correction_plan::{BarcodeCorrection, CellCorrectionScope, CorrectionPlan};
use crate::correction_spool::{DeferredPairSpool, SpoolStats};
use crate::diagnostics;
use crate::knee_finding;
use crate::prog_opts::GenPermitListOpts;
use crate::prog_opts::SampleBarcodeOri;
use crate::utils as afutils;
use crate::utils::KnownRecordType;
#[allow(unused_imports)]
use ahash::{AHashMap, AHashSet, AHasher, RandomState};
use anyhow::{Context, anyhow, bail};
use bio_types::strand::Strand;
use bstr::io::BufReadExt;
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use libradicl::exit_codes;
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{self, RadType, TagMap, TagSection};
use libradicl::record::{
    AlevinFryReadRecordWithPosition, CollatableMappedRecord, ConvertiblePrimitiveInteger,
    HierarchicallyCollatable, KnownSize, MappedRecord, MultiBarcodeReadRecord,
    MultiBarcodeRecordContext, RecordContext, ScLongReadRecord,
};
use libradicl::{chunk, record::AlevinFryReadRecord};
use needletail::bitkmer::*;
use num_format::{Locale, ToFormattedString};
use serde::Serialize;
use serde_json::json;
use slog::crit;
use slog::{info, warn};
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::hash_map::Entry;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::{Arc, atomic::Ordering};
use std::time::Instant;

#[derive(Clone, Debug, Serialize)]
pub enum CellFilterMethod {
    // cut off at this cell in
    // the frequency sorted list
    ForceCells(usize),
    // use this cell as a hint in
    // the frequency sorted list
    ExpectCells(usize),
    // correct all cells in an
    // edit distance of 1 of these
    // barcodes
    ExplicitList(PathBuf),
    // barcodes will be provided in the
    // form of an *unfiltered* external
    // permit list
    UnfilteredExternalList(PathBuf, usize),
    // use the distance method to
    // automatically find the knee
    // in the curve
    KneeFinding,
}

fn populate_unfiltered_barcode_map<T: Read>(
    br: BufReader<T>,
    first_bclen: &mut usize,
) -> DashMap<u64, u64, ahash::RandomState> {
    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let hm = DashMap::with_hasher(s);

    // read through the external unfiltered barcode list
    // and generate a vector of encoded barcodes
    // let mut kv = Vec::<u64>::new();
    for l in br.byte_lines().flatten() {
        if *first_bclen == 0 {
            *first_bclen = l.len();
        } else {
            assert_eq!(
                *first_bclen,
                l.len(),
                "found barcodes of different lengths {} and {}",
                *first_bclen,
                l.len()
            );
        }
        if let Some((_, km, _)) =
            needletail::bitkmer::BitNuclKmer::new(&l[..], l.len() as u8, false).next()
        {
            hm.insert(km.0, 0);
        }
    }
    hm
}

type BarcodeCountMap = HashMap<u64, u64, ahash::RandomState>;

/// The policy-free sample routing information needed by the RAD scan.
///
/// GPL has already resolved `observed -> canonical sample` before this is
/// built.  Compiling the canonical target to its dense sample index here
/// lets each aligned record use one fast lookup instead of separately
/// consulting the permit map, exact-source map, and canonical-index map.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ImmediateSampleRoute {
    sample_index: usize,
    exact_source: bool,
}

fn compile_immediate_sample_routes(
    sample_permit_map: &HashMap<u64, u64>,
    canonical_to_index: &HashMap<u64, usize>,
    exact_sources: &HashMap<u64, u64>,
) -> anyhow::Result<AHashMap<u64, ImmediateSampleRoute>> {
    let mut routes = AHashMap::with_capacity(sample_permit_map.len());
    for (&observed, corrected) in sample_permit_map {
        let &sample_index = canonical_to_index.get(corrected).ok_or_else(|| {
            anyhow!("sample correction target 0x{corrected:x} has no canonical sample index")
        })?;
        routes.insert(
            observed,
            ImmediateSampleRoute {
                sample_index,
                exact_source: exact_sources.contains_key(&observed),
            },
        );
    }
    Ok(routes)
}

fn flush_unmatched_counts(
    local: &mut [BarcodeCountMap],
    shared: &[std::sync::Mutex<BarcodeCountMap>],
) {
    for (sample_idx, buffer) in local.iter_mut().enumerate() {
        if buffer.is_empty() {
            continue;
        }
        let mut target = shared[sample_idx].lock().unwrap();
        for (cell, count) in buffer.drain() {
            *target.entry(cell).or_insert(0) += count;
        }
    }
}

fn flush_cell_counts(
    local: &mut [BarcodeCountMap],
    shared: &[DashMap<u64, u64, ahash::RandomState>],
) {
    for (sample_idx, buffer) in local.iter_mut().enumerate() {
        for (cell, count) in buffer.drain() {
            *shared[sample_idx].entry(cell).or_insert(0) += count;
        }
    }
}

#[inline]
fn buffer_count(
    buffer: &mut BarcodeCountMap,
    barcode: u64,
    count: u64,
    buffered_distinct: &mut usize,
) {
    match buffer.entry(barcode) {
        Entry::Occupied(mut entry) => *entry.get_mut() += count,
        Entry::Vacant(entry) => {
            entry.insert(count);
            *buffered_distinct += 1;
        }
    }
}

struct MultiSampleCellOutput {
    sample_idx: usize,
    num_cells: u64,
    info: serde_json::Value,
    scope: CellCorrectionScope,
}

#[allow(clippy::too_many_arguments)]
fn compile_multi_sample_cell_output(
    sample_idx: usize,
    sample_bc: u64,
    sample_name: &str,
    cell_hist: DashMap<u64, u64, ahash::RandomState>,
    unmatched_cells: std::sync::Mutex<BarcodeCountMap>,
    parent: &Path,
    cell_spec: CorrectionSpec,
    filter_method: &CellFilterMethod,
    min_reads: u64,
    cell_bc_len: u16,
    log: &slog::Logger,
) -> anyhow::Result<MultiSampleCellOutput> {
    let sample_dir = parent.join(format!("sample_{sample_name}"));
    std::fs::create_dir_all(&sample_dir)?;

    // Consume the sample histogram rather than retaining the DashMap and a
    // cloned HashMap simultaneously across the whole experiment.
    let cell_hist: BarcodeCountMap = cell_hist.into_iter().collect();
    info!(
        log,
        "Sample '{}' (bc=0x{:x}): {} distinct cell barcodes observed",
        sample_name,
        sample_bc,
        cell_hist.len().to_formatted_string(&Locale::en),
    );

    // Preserve the historical multi-sample output contract: unused samples
    // have metadata/directories but no empty permit files. The compiled plan
    // remains total so its empty lookup rejects any such records directly.
    if cell_hist.is_empty() {
        warn!(
            log,
            "Sample '{}' has no reads — skipping permit list generation", sample_name
        );
        return Ok(MultiSampleCellOutput {
            sample_idx,
            num_cells: 0,
            info: serde_json::json!({
                "name": sample_name,
                "barcode": format!("0x{:x}", sample_bc),
                "num_reads": 0_u64,
                "num_cells": 0_u64,
                "collatable_reads": 0_u64,
                "cell_correction": CorrectionStats::default(),
            }),
            scope: CellCorrectionScope {
                sample_barcode: Some(sample_bc),
                spec: cell_spec,
                corrections: Vec::new(),
            },
        });
    }

    let mut kept_bcs =
        select_retained_barcodes(&cell_hist, filter_method, min_reads, cell_bc_len, log);
    let mut observed_counts = cell_hist;
    for (barcode, count) in unmatched_cells.into_inner().unwrap() {
        *observed_counts.entry(barcode).or_insert(0) += count;
    }

    kept_bcs.sort_unstable();
    kept_bcs.dedup();
    // Historically, filtered multi-barcode modes reported only exact reads
    // assigned to the selected cells in sample_info.json.  Unfiltered mode
    // already included rescued reads. Keep that field stable while exposing
    // the complete compiled-correction total separately below.
    let exact_retained_reads = kept_bcs
        .iter()
        .map(|barcode| observed_counts.get(barcode).copied().unwrap_or(0))
        .sum::<u64>();
    info!(
        log,
        "  {} retained cell-barcode targets for sample '{}'",
        kept_bcs.len(),
        sample_name,
    );

    let correction_index = CorrectionIndex::new(
        cell_spec,
        kept_bcs.iter().map(|&barcode| RetainedSource {
            source: barcode,
            target: barcode,
            exact_count: observed_counts.get(&barcode).copied().unwrap_or(0),
        }),
    )?;
    let (observed_table, target_counts) = correction_index
        .compile_distinct_observed_with_target_counts(
            observed_counts
                .iter()
                .map(|(&barcode, &count)| (barcode, count)),
        );
    let sample_freq_map: BarcodeCountMap = target_counts.into_iter().collect();

    // Keep legacy permit_map.bin with its theoretical-neighbour semantics,
    // while compiling it under the same policy as the compact artifact.
    let full_table = (!matches!(
        filter_method,
        CellFilterMethod::UnfilteredExternalList(_, _)
    ))
    .then(|| correction_index.compile_full_neighborhood());
    let legacy_entries = full_table
        .as_ref()
        .map_or(observed_table.entries(), |table| table.entries());
    let full_permit_list: HashMap<u64, u64> = legacy_entries
        .iter()
        .copied()
        .filter(|(_, target)| sample_freq_map.contains_key(target))
        .collect();
    let map_file = File::create(sample_dir.join("permit_map.bin"))?;
    bincode::serialize_into(BufWriter::new(map_file), &full_permit_list)
        .map_err(|e| anyhow!("failed to serialize permit map: {}", e))?;

    let correction_stats = observed_table.stats();
    info!(
        log,
        "  Cell correction for sample '{}': {} exact reads, {} corrected, {} ambiguous, {} without a candidate",
        sample_name,
        correction_stats.exact_reads,
        correction_stats.corrected_reads,
        correction_stats.ambiguous_reads,
        correction_stats.not_found_reads,
    );

    let mut accepted_entries = observed_table.into_entries();
    accepted_entries.retain(|(_, target)| sample_freq_map.contains_key(target));
    let scope = CellCorrectionScope {
        sample_barcode: Some(sample_bc),
        spec: cell_spec,
        corrections: accepted_entries
            .into_iter()
            .map(|(observed, corrected)| BarcodeCorrection {
                observed,
                corrected,
            })
            .collect(),
    };

    afutils::write_permit_list_freq(
        &sample_dir.join("permit_freq.bin"),
        cell_bc_len,
        &sample_freq_map,
    )
    .map_err(|e| {
        anyhow!(
            "failed to write permit freq for sample '{}': {}",
            sample_name,
            e
        )
    })?;

    let collatable_reads: u64 = sample_freq_map.values().sum();
    let historical_num_reads = if matches!(
        filter_method,
        CellFilterMethod::UnfilteredExternalList(_, _)
    ) {
        collatable_reads
    } else {
        exact_retained_reads
    };
    Ok(MultiSampleCellOutput {
        sample_idx,
        num_cells: sample_freq_map.len() as u64,
        info: serde_json::json!({
            "name": sample_name,
            "barcode": format!("0x{:x}", sample_bc),
            "num_reads": historical_num_reads,
            "num_cells": sample_freq_map.len(),
            "collatable_reads": collatable_reads,
            "cell_correction": correction_stats,
        }),
        scope,
    })
}

#[allow(clippy::unnecessary_unwrap, clippy::too_many_arguments)]
fn process_unfiltered(
    hm: DashMap<u64, u64, ahash::RandomState>,
    unmatched_bc: HashMap<u64, u64, ahash::RandomState>,
    file_tag_map: &rad_types::TagMap,
    filter_meth: &CellFilterMethod,
    expected_ori: Strand,
    output_dir: &PathBuf,
    version: &str,
    max_ambiguity_read: usize,
    velo_mode: bool,
    cmdline: &str,
    log: &slog::Logger,
    gpl_opts: &GenPermitListOpts,
) -> anyhow::Result<u64> {
    let parent = std::path::Path::new(output_dir);
    std::fs::create_dir_all(parent)
        .with_context(|| format!("couldn't create directory path {}", parent.display()))?;
    let min_freq = match filter_meth {
        CellFilterMethod::UnfilteredExternalList(_, min_reads) => {
            info!(log, "minimum num reads for barcode pass = {}", *min_reads);
            *min_reads as u64
        }
        _ => {
            unimplemented!();
        }
    };
    let barcode_tag = file_tag_map
        .get("cblen")
        .expect("tag map must contain cblen");
    let barcode_len: u16 = barcode_tag.try_into()?;

    let mut observed_counts: HashMap<u64, u64, ahash::RandomState> = HashMap::default();
    let mut kept_barcodes = Vec::new();
    for entry in &hm {
        if *entry.value() > 0 {
            observed_counts.insert(*entry.key(), *entry.value());
        }
        if *entry.value() >= min_freq {
            kept_barcodes.push(*entry.key());
        }
    }
    for (barcode, count) in unmatched_bc {
        *observed_counts.entry(barcode).or_insert(0) += count;
    }
    kept_barcodes.sort_unstable();

    let spec = gpl_opts.cell_correction_spec(barcode_len as u8, false);
    warn_for_shift_frequency(spec, "cell", log);
    let correction_index = CorrectionIndex::new(
        spec,
        kept_barcodes.iter().map(|&barcode| RetainedSource {
            source: barcode,
            target: barcode,
            exact_count: observed_counts.get(&barcode).copied().unwrap_or(0),
        }),
    )?;
    info!(
        log,
        "found {} cells with non-trivial number of reads by exact barcode match",
        kept_barcodes.len().to_formatted_string(&Locale::en)
    );

    let started = Instant::now();
    let (correction_table, target_counts) = correction_index
        .compile_distinct_observed_with_target_counts(
            observed_counts
                .iter()
                .map(|(&barcode, &count)| (barcode, count)),
        );
    let permitted_map: HashMap<u64, u64, ahash::RandomState> = target_counts.into_iter().collect();
    let stats = correction_table.stats();
    info!(
        log,
        "Cell correction took {:?}: {} exact reads, {} corrected, {} ambiguous, {} without a candidate",
        started.elapsed(),
        stats.exact_reads,
        stats.corrected_reads,
        stats.ambiguous_reads,
        stats.not_found_reads,
    );

    afutils::write_permit_list_freq(&parent.join("permit_freq.bin"), barcode_len, &permitted_map)
        .map_err(|error| anyhow!("could not write permit frequencies: {error}"))?;
    let permit_map: HashMap<u64, u64, ahash::RandomState> = correction_table
        .entries()
        .iter()
        .copied()
        .filter(|(_, target)| permitted_map.contains_key(target))
        .collect();
    let pm_file = File::create(parent.join("permit_map.bin"))?;
    bincode::serialize_into(BufWriter::new(pm_file), &permit_map)
        .context("couldn't serialize permit list mapping")?;

    let compact_entries = correction_table
        .entries()
        .iter()
        .copied()
        .filter(|(_, target)| permitted_map.contains_key(target))
        .map(|(observed, corrected)| BarcodeCorrection {
            observed,
            corrected,
        })
        .collect();
    let mut correction_plan = CorrectionPlan {
        sample_barcode_len: None,
        cell_barcode_len: barcode_len as u8,
        sample_spec: None,
        sample_corrections: Vec::new(),
        cell_scopes: vec![CellCorrectionScope {
            sample_barcode: None,
            spec,
            corrections: compact_entries,
        }],
    };
    correction_plan.write_to(parent)?;

    let meta_info = json!({
    "velo_mode" : velo_mode,
    "expected_ori" : *expected_ori.strand_symbol(),
    "version_str" : version,
    "max-ambig-record" : max_ambiguity_read,
    "cmd" : cmdline,
    "permit-list-type" : "unfiltered",
    "gpl_options" : &gpl_opts,
    "resolved_cell_bc_neighborhood": spec.neighborhood.to_string(),
    "resolved_cell_bc_confidence": gpl_opts.cell_bc_confidence.to_string(),
    "correction_stats": stats,
    });

    let m_path = parent.join("generate_permit_list.json");
    let mut m_file = std::fs::File::create(m_path).context("could not create metadata file.")?;

    let meta_info_string =
        serde_json::to_string_pretty(&meta_info).context("could not format json.")?;
    m_file
        .write_all(meta_info_string.as_bytes())
        .context("cannot write to generate_permit_list.json file")?;

    info!(
        log,
        "total number of distinct corrected barcodes : {}",
        stats.corrected_distinct.to_formatted_string(&Locale::en)
    );

    Ok(stats.corrected_distinct)
}

#[allow(clippy::unnecessary_unwrap, clippy::too_many_arguments)]
fn process_filtered(
    hm: DashMap<u64, u64, ahash::RandomState>,
    file_tag_map: &rad_types::TagMap,
    filter_meth: &CellFilterMethod,
    expected_ori: Strand,
    output_dir: &PathBuf,
    version: &str,
    max_ambiguity_read: usize,
    velo_mode: bool,
    cmdline: &str,
    log: &slog::Logger,
    gpl_opts: &GenPermitListOpts,
) -> anyhow::Result<u64> {
    let hm: HashMap<u64, u64, ahash::RandomState> = hm.into_iter().collect();
    let barcode_tag = file_tag_map
        .get("cblen")
        .expect("tag map must contain cblen");
    let barcode_len: u16 = barcode_tag.try_into()?;
    let mut valid_barcodes = select_retained_barcodes(&hm, filter_meth, 0, barcode_len, log);
    valid_barcodes.sort_unstable();
    valid_barcodes.dedup();
    info!(
        log,
        "filtering selected {} retained barcode targets",
        valid_barcodes.len()
    );

    let spec = gpl_opts.cell_correction_spec(barcode_len as u8, false);
    warn_for_shift_frequency(spec, "cell", log);
    let correction_index = CorrectionIndex::new(
        spec,
        valid_barcodes.iter().map(|&barcode| RetainedSource {
            source: barcode,
            target: barcode,
            exact_count: hm.get(&barcode).copied().unwrap_or(0),
        }),
    )?;
    let (correction_table, target_counts) = correction_index
        .compile_distinct_observed_with_target_counts(
            hm.iter().map(|(&barcode, &count)| (barcode, count)),
        );
    let permitted_map: HashMap<u64, u64, ahash::RandomState> = target_counts.into_iter().collect();
    let full_permit_list: HashMap<u64, u64> = correction_index
        .compile_full_neighborhood()
        .entries()
        .iter()
        .copied()
        .filter(|(_, target)| permitted_map.contains_key(target))
        .collect();
    let stats = correction_table.stats();

    let parent = std::path::Path::new(output_dir);
    std::fs::create_dir_all(parent).with_context(|| {
        format!(
            "failed to create path to output location {}",
            parent.display()
        )
    })?;
    let o_path = parent.join("permit_freq.bin");

    match afutils::write_permit_list_freq(&o_path, barcode_len, &permitted_map) {
        Ok(_) => {}
        Err(error) => {
            panic!("Error: {}", error);
        }
    };

    let o_path = parent.join("all_freq.bin");

    match afutils::write_permit_list_freq(&o_path, barcode_len, &hm) {
        Ok(_) => {}
        Err(error) => {
            panic!("Error: {}", error);
        }
    };

    let s_path = parent.join("permit_map.bin");
    let s_file = std::fs::File::create(s_path).context("could not create serialization file.")?;
    let mut s_writer = BufWriter::new(&s_file);
    bincode::serialize_into(&mut s_writer, &full_permit_list)
        .context("couldn't serialize permit list.")?;

    let compact_entries = correction_table
        .entries()
        .iter()
        .copied()
        .filter(|(_, target)| permitted_map.contains_key(target))
        .map(|(observed, corrected)| BarcodeCorrection {
            observed,
            corrected,
        })
        .collect();
    let mut correction_plan = CorrectionPlan {
        sample_barcode_len: None,
        cell_barcode_len: barcode_len as u8,
        sample_spec: None,
        sample_corrections: Vec::new(),
        cell_scopes: vec![CellCorrectionScope {
            sample_barcode: None,
            spec,
            corrections: compact_entries,
        }],
    };
    correction_plan.write_to(parent)?;

    let meta_info = json!({
    "velo_mode" : velo_mode,
    "expected_ori" : *expected_ori.strand_symbol(),
    "version_str" : version,
    "max-ambig-record" : max_ambiguity_read,
    "cmd" : cmdline,
    "permit-list-type" : "filtered",
    "gpl_options" : &gpl_opts,
    "resolved_cell_bc_neighborhood": spec.neighborhood.to_string(),
    "resolved_cell_bc_confidence": gpl_opts.cell_bc_confidence.to_string(),
    "correction_stats": stats,
    });

    let m_path = parent.join("generate_permit_list.json");
    let mut m_file = std::fs::File::create(m_path).context("could not create metadata file.")?;

    let meta_info_string =
        serde_json::to_string_pretty(&meta_info).context("could not format json.")?;
    m_file
        .write_all(meta_info_string.as_bytes())
        .context("cannot write to generate_permit_list.json file")?;

    info!(
        log,
        "total number of distinct corrected barcodes : {}",
        stats.corrected_distinct.to_formatted_string(&Locale::en)
    );

    Ok(stats.corrected_distinct)
}

/// Given the input RAD file `input_file`, compute
/// and output (in `output_dir`) the list of valid
/// (i.e. "permitted") barcode values, as well as
/// a map from each correctable barcode to the
/// permitted barcode to which it maps.
pub fn generate_permit_list(gpl_opts: GenPermitListOpts) -> anyhow::Result<u64> {
    let rad_dir = gpl_opts.input_dir.clone();
    let log = gpl_opts.log;

    let i_dir = std::path::Path::new(&rad_dir);
    let i_file = File::open(i_dir.join("map.rad")).context("could not open input rad file")?;
    let mut ifile = BufReader::new(i_file);

    // should we assume this condition was already checked
    // during parsing?
    if !i_dir.exists() {
        crit!(
            log,
            "the input RAD path {} does not exist",
            rad_dir.display()
        );
        // std::process::exit(1);
        anyhow::bail!("execution terminated because input RAD path does not exist.");
    }

    let prelude = RadPrelude::from_bytes(&mut ifile).unwrap();
    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut ifile).unwrap();
    let rec_type = afutils::get_record_type_from_prelude(&prelude, &file_tag_map);

    match rec_type {
        KnownRecordType::RnaLong(_bc_len) => {
            info!(log, "record type is long read single-cell RNA-seq");
            do_generate_permit_list::<u64, ScLongReadRecord>(gpl_opts, ifile, prelude, file_tag_map)
        }
        KnownRecordType::AtacSeq(_bc_len) => {
            info!(log, "record type is short read single-cell ATAC-seq");
            anyhow::bail!("To process atac-seq data, you should use the \"atac\" sub-command");
        }
        KnownRecordType::RnaShortPos(_bc_len) => {
            info!(
                log,
                "record type is short read single-cell RNA-seq with positions"
            );
            do_generate_permit_list::<u64, AlevinFryReadRecordWithPosition>(
                gpl_opts,
                ifile,
                prelude,
                file_tag_map,
            )
        }
        KnownRecordType::RnaShort(_bc_len) => {
            info!(
                log,
                "record type is standard short read single-cell RNA-seq"
            );
            do_generate_permit_list::<u64, AlevinFryReadRecord>(
                gpl_opts,
                ifile,
                prelude,
                file_tag_map,
            )
        }
        KnownRecordType::RnaShortMultiBC(cell_bc_len, num_bc) => {
            info!(
                log,
                "record type is multi-barcode single-cell RNA-seq ({} barcode levels, cell BC len = {})",
                num_bc,
                cell_bc_len,
            );
            do_generate_permit_list_multi_bc(gpl_opts, ifile, prelude, file_tag_map, num_bc)
        }
    }
}

fn warn_for_shift_frequency(spec: CorrectionSpec, kind: &str, log: &slog::Logger) {
    if matches!(spec.resolution, BarcodeResolution::Frequency { .. })
        && spec.neighborhood
            == crate::barcode_correction::BarcodeNeighborhood::SubstitutionOrShiftOne
    {
        warn!(
            log,
            "{}-barcode Frequency correction is using substitution-or-shift-1. RAD stores no barcode base qualities or indel-error model, so this is an abundance heuristic rather than a calibrated indel posterior.",
            kind
        );
    }
}

fn select_retained_barcodes(
    histogram: &HashMap<u64, u64, ahash::RandomState>,
    method: &CellFilterMethod,
    unfiltered_min_reads: u64,
    barcode_len: u16,
    log: &slog::Logger,
) -> Vec<u64> {
    if let CellFilterMethod::ExplicitList(path) = method {
        return permit_list_from_file(path, barcode_len);
    }
    if histogram.is_empty() {
        return Vec::new();
    }

    let mut frequencies: Vec<u64> = histogram.values().copied().collect();
    frequencies.sort_unstable_by(|lhs, rhs| rhs.cmp(lhs));
    let threshold = match method {
        CellFilterMethod::UnfilteredExternalList(_, _) => unfiltered_min_reads,
        CellFilterMethod::KneeFinding => {
            let knee = knee_finding::get_knee(&frequencies, 100, log);
            frequencies[knee.min(frequencies.len() - 1)]
        }
        CellFilterMethod::ForceCells(number) => {
            if *number == 0 {
                return Vec::new();
            }
            frequencies[number.saturating_sub(1).min(frequencies.len() - 1)]
        }
        CellFilterMethod::ExpectCells(expected) => {
            let robust_index =
                ((*expected as f64 * 0.99).round() as usize).min(frequencies.len() - 1);
            std::cmp::max(1, (frequencies[robust_index] as f64 / 10.0).round() as u64)
        }
        CellFilterMethod::ExplicitList(_) => unreachable!(),
    };

    histogram
        .iter()
        .filter_map(|(&barcode, &count)| (count >= threshold).then_some(barcode))
        .collect()
}

/// Multi-barcode generate-permit-list implementation.
///
/// For protocols like 10x Flex with multiple barcodes per read (e.g., sample + cell):
/// 1. Reads all records from the multi-barcode RAD file
/// 2. Corrects sample barcodes against the provided known list
/// 3. Per-sample: counts cell barcode frequencies and generates permit lists
/// 4. Outputs: sample_permit_map.bin, per-sample permit_map.bin/permit_freq.bin, sample_info.json
fn do_generate_permit_list_multi_bc(
    gpl_opts: GenPermitListOpts,
    mut ifile: BufReader<File>,
    prelude: RadPrelude,
    file_tag_map: TagMap,
    num_barcodes: u16,
) -> anyhow::Result<u64> {
    let log = gpl_opts.log;
    let output_dir = gpl_opts.output_dir;
    let expected_ori = gpl_opts.expected_ori;

    // Require sample barcode list for multi-barcode mode
    let sample_bc_list_path = gpl_opts.sample_bc_list.as_ref().ok_or_else(|| {
        anyhow!(
            "Multi-barcode RAD file detected ({} barcode levels), but --sample-bc-list was not provided. \
             A known sample barcode list is required for multi-barcode processing.",
            num_barcodes,
        )
    })?;

    // Parse the record context for multi-barcode records
    let rec_ctx = prelude.get_record_context::<MultiBarcodeRecordContext>()?;
    info!(
        log,
        "Multi-barcode record context: {} barcode levels",
        rec_ctx.num_barcodes()
    );

    // Load known sample barcodes (with rotation → canonical mapping)
    let sample_info = load_sample_barcode_list(sample_bc_list_path, gpl_opts.sample_bc_ori, log)?;

    // Build sample barcode correction map (rotation → canonical)
    let sample_correction_spec = gpl_opts.sample_correction_spec(sample_info.barcode_len as u8);
    let sample_routing = build_sample_permit_map(&sample_info, sample_correction_spec, log)?;
    let sample_permit_map = sample_routing.permit_map;
    let sample_bc_to_idx = sample_routing.canonical_to_index;
    let ambiguous_sample_barcodes = sample_routing.ambiguous_sample_barcodes;
    let sample_frequency = sample_correction_spec
        .is_some_and(|spec| matches!(spec.resolution, BarcodeResolution::Frequency { .. }));
    let immediate_sample_routes = compile_immediate_sample_routes(
        &sample_permit_map,
        &sample_bc_to_idx,
        &sample_info.rotation_to_canonical,
    )?;

    // Get sample names from the barcode file (uses canonical → name mapping)
    let sample_names = get_sample_names(&sample_info);

    // Load external cell barcode whitelist if provided (e.g. 737K 10x list).
    // This is used to filter cell barcodes — only whitelist BCs are counted directly;
    // non-matching BCs are collected for 1-edit correction later.
    let (cell_bc_whitelist, min_reads) = match &gpl_opts.fmeth {
        CellFilterMethod::UnfilteredExternalList(wl_path, mr) => {
            info!(
                log,
                "Loading external cell barcode whitelist from {}",
                wl_path.display()
            );
            let wl_file = File::open(wl_path)
                .with_context(|| format!("couldn't open whitelist {}", wl_path.display()))?;
            let mut first_bclen = 0usize;
            let wl = populate_unfiltered_barcode_map(BufReader::new(wl_file), &mut first_bclen);
            let wl: AHashSet<u64> = wl.iter().map(|entry| *entry.key()).collect();
            info!(
                log,
                "Loaded {} cell barcodes from whitelist (bclen={})",
                wl.len(),
                first_bclen
            );
            (Some(wl), *mr as u64)
        }
        _ => (None, 0u64),
    };

    // First pass: read all records in parallel, correct sample BCs, count cell BCs per sample.
    // When a whitelist is present, only count cell BCs that are in the whitelist;
    // collect non-matching BCs per-sample for later 1-edit correction.
    info!(
        log,
        "First pass: counting cell barcodes per sample ({} threads)...", gpl_opts.threads
    );
    let num_samples = sample_info.canonical_barcodes.len();
    // Per-sample cell barcode frequency (only whitelist-matching BCs) — thread-safe
    let per_sample_cell_hist: Vec<DashMap<u64, u64, ahash::RandomState>> =
        (0..num_samples).map(|_| DashMap::default()).collect();
    // Per-sample unmatched cell barcode counts. Counts avoid the historical
    // occurrence-sized vectors while retaining exact read totals.
    let per_sample_unmatched: Vec<std::sync::Mutex<HashMap<u64, u64, ahash::RandomState>>> = (0
        ..num_samples)
        .map(|_| std::sync::Mutex::new(HashMap::default()))
        .collect();

    let total_reads = std::sync::atomic::AtomicU64::new(0);
    let matched_reads = std::sync::atomic::AtomicU64::new(0);
    let unmatched_reads = std::sync::atomic::AtomicU64::new(0);
    let exact_sample_reads = std::sync::atomic::AtomicU64::new(0);
    let structurally_unique_sample_reads = std::sync::atomic::AtomicU64::new(0);
    let deferred_sample_reads = std::sync::atomic::AtomicU64::new(0);

    let nworkers = gpl_opts.threads.max(1);
    let mut chunk_reader = libradicl::readers::ParallelChunkReader::<MultiBarcodeReadRecord>::new(
        &prelude,
        std::num::NonZeroUsize::new(nworkers).unwrap(),
    );

    let num_chunks = prelude.hdr.num_chunks;
    let pbar =
        ProgressBar::with_draw_target(Some(num_chunks), ProgressDrawTarget::stderr_with_hz(5));
    pbar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} chunks (ETA: {eta})",
        )
        .unwrap()
        .progress_chars("##-"),
    );

    // Wrap shared state in Arcs for the worker threads
    let per_sample_cell_hist = std::sync::Arc::new(per_sample_cell_hist);
    let per_sample_unmatched = std::sync::Arc::new(per_sample_unmatched);
    let immediate_sample_routes = std::sync::Arc::new(immediate_sample_routes);
    let ambiguous_sample_barcodes = std::sync::Arc::new(ambiguous_sample_barcodes);
    let cell_bc_whitelist = std::sync::Arc::new(cell_bc_whitelist);
    let total_reads_arc = &total_reads;
    let matched_reads_arc = &matched_reads;
    let unmatched_reads_arc = &unmatched_reads;
    let tmp_dir = gpl_opts
        .tmp_dir
        .clone()
        .unwrap_or_else(|| output_dir.clone());
    let spool_pairs_per_worker = if sample_frequency {
        let bytes = gpl_opts.effective_memory_limit() / nworkers as u64;
        usize::try_from(bytes / std::mem::size_of::<(u64, u64)>() as u64)
            .unwrap_or(usize::MAX)
            .max(1)
    } else {
        1
    };

    let (mut correction_spools, sample_exact_counts) = std::thread::scope(
        |s| -> anyhow::Result<(Vec<DeferredPairSpool>, BarcodeCountMap)> {
            // Spawn worker threads
            let mut handles = Vec::new();
            for worker in 0..nworkers {
                let chunks = chunk_reader.chunk_iter();
                let hist = per_sample_cell_hist.clone();
                let unmatched = per_sample_unmatched.clone();
                let sample_routes = immediate_sample_routes.clone();
                let ambiguous_samples = ambiguous_sample_barcodes.clone();
                let wl = cell_bc_whitelist.clone();
                let spool_dir = tmp_dir.clone();

                let handle = s.spawn(move || -> anyhow::Result<_> {
                    let mut local_total = 0u64;
                    let mut local_matched = 0u64;
                    let mut local_unmatched_sample = 0u64;
                    let mut local_exact_sample = 0u64;
                    let mut local_structurally_unique_sample = 0u64;
                    let mut local_deferred_sample = 0u64;
                    let mut local_sample_exact_counts = BarcodeCountMap::default();
                    let num_s = hist.len();
                    let mut local_cell_bufs: Vec<BarcodeCountMap> =
                        (0..num_s).map(|_| HashMap::default()).collect();
                    let mut local_unmatched_bufs: Vec<BarcodeCountMap> =
                        (0..num_s).map(|_| HashMap::default()).collect();
                    let mut local_buffered_distinct = 0usize;
                    // Bound the combined exact and unmatched per-worker
                    // histograms. Local aggregation avoids a concurrent map
                    // operation for every read without allowing every worker
                    // to retain a full copy of the barcode universe.
                    const CELL_COUNT_FLUSH_ENTRIES: usize = 1 << 16;
                    let mut spool = if sample_frequency {
                        Some(DeferredPairSpool::new(
                            &spool_dir,
                            worker,
                            spool_pairs_per_worker,
                        )?)
                    } else {
                        None
                    };

                    for meta_chunk in chunks {
                        for c in meta_chunk.iter() {
                            for read in &c.reads {
                                local_total += 1;

                                if !read.has_alignment_on_strand(expected_ori) {
                                    continue;
                                }

                                let sample_bc: u64 = read.collation_key_at_level(0);
                                let cell_bc: u64 = read.collate_key();

                                if let Some(route) = sample_routes.get(&sample_bc) {
                                    if route.exact_source {
                                        if sample_frequency {
                                            *local_sample_exact_counts
                                                .entry(sample_bc)
                                                .or_insert(0) += 1;
                                        }
                                        local_exact_sample += 1;
                                    } else {
                                        local_structurally_unique_sample += 1;
                                    }
                                    let sample_idx = route.sample_index;
                                    local_matched += 1;
                                    if let Some(ref wl_map) = *wl {
                                        if wl_map.contains(&cell_bc) {
                                            buffer_count(
                                                &mut local_cell_bufs[sample_idx],
                                                cell_bc,
                                                1,
                                                &mut local_buffered_distinct,
                                            );
                                        } else {
                                            buffer_count(
                                                &mut local_unmatched_bufs[sample_idx],
                                                cell_bc,
                                                1,
                                                &mut local_buffered_distinct,
                                            );
                                        }
                                    } else {
                                        buffer_count(
                                            &mut local_cell_bufs[sample_idx],
                                            cell_bc,
                                            1,
                                            &mut local_buffered_distinct,
                                        );
                                    }
                                    if local_buffered_distinct >= CELL_COUNT_FLUSH_ENTRIES {
                                        flush_cell_counts(&mut local_cell_bufs, &hist);
                                        flush_unmatched_counts(
                                            &mut local_unmatched_bufs,
                                            &unmatched,
                                        );
                                        local_buffered_distinct = 0;
                                    }
                                } else if sample_frequency
                                    && ambiguous_samples
                                        .as_ref()
                                        .as_ref()
                                        .expect("Frequency routing has an ambiguity index")
                                        .contains(&sample_bc)
                                {
                                    spool
                                        .as_mut()
                                        .expect("Frequency routing has a spool")
                                        .push(sample_bc, cell_bc)?;
                                    local_deferred_sample += 1;
                                } else {
                                    local_unmatched_sample += 1;
                                }
                            }
                        }
                    }

                    flush_cell_counts(&mut local_cell_bufs, &hist);
                    flush_unmatched_counts(&mut local_unmatched_bufs, &unmatched);

                    Ok((
                        local_total,
                        local_matched,
                        local_unmatched_sample,
                        local_exact_sample,
                        local_structurally_unique_sample,
                        local_deferred_sample,
                        local_sample_exact_counts,
                        spool,
                    ))
                });
                handles.push(handle);
            }

            // Start the chunk reader (this feeds the queue from the BufReader)
            let cb = |_new_bytes: u64, new_rec: u64| {
                pbar.inc(new_rec);
            };
            chunk_reader.start(&mut ifile, Some(cb))?;

            // Join workers and aggregate counts
            let mut spools = Vec::new();
            let mut exact_counts = BarcodeCountMap::default();
            for handle in handles {
                let (lt, lm, lu, le, ls, ld, local_exact_counts, spool) =
                    handle.join().expect("worker thread panicked")?;
                total_reads_arc.fetch_add(lt, Ordering::SeqCst);
                matched_reads_arc.fetch_add(lm, Ordering::SeqCst);
                unmatched_reads_arc.fetch_add(lu, Ordering::SeqCst);
                exact_sample_reads.fetch_add(le, Ordering::SeqCst);
                structurally_unique_sample_reads.fetch_add(ls, Ordering::SeqCst);
                deferred_sample_reads.fetch_add(ld, Ordering::SeqCst);
                for (barcode, count) in local_exact_counts {
                    *exact_counts.entry(barcode).or_insert(0) += count;
                }
                if let Some(spool) = spool {
                    spools.push(spool);
                }
            }
            Ok((spools, exact_counts))
        },
    )?;
    pbar.finish_and_clear();

    let total_reads = total_reads.load(Ordering::SeqCst);

    let mut sample_permit_map = sample_permit_map;

    let mut replay_time = std::time::Duration::ZERO;
    let mut spool_stats = SpoolStats::default();
    let immediate_rejected_sample_reads = unmatched_reads.load(Ordering::SeqCst);
    let mut deferred_accepted_sample_reads = 0_u64;
    let mut deferred_rejected_sample_reads = 0_u64;
    if let Some(spec) = sample_correction_spec
        .filter(|spec| matches!(spec.resolution, BarcodeResolution::Frequency { .. }))
    {
        let frequency_index = CorrectionIndex::new(
            spec,
            sample_info
                .rotation_to_canonical
                .iter()
                .map(|(&source, &target)| RetainedSource {
                    source,
                    target,
                    exact_count: sample_exact_counts.get(&source).copied().unwrap_or(0),
                }),
        )?;
        let table = frequency_index.compile_full_neighborhood();
        sample_permit_map = table.entries().iter().copied().collect();

        let replay_start = Instant::now();
        for spool in &mut correction_spools {
            spool.replay(|raw_sample, cell, count| {
                match frequency_index.resolve(raw_sample) {
                    CorrectionDecision::Exact(corrected)
                    | CorrectionDecision::Corrected(corrected) => {
                        deferred_accepted_sample_reads += count;
                        let sample_idx = sample_bc_to_idx[&corrected];
                        matched_reads.fetch_add(count, Ordering::SeqCst);
                        if let Some(ref whitelist) = *cell_bc_whitelist {
                            if whitelist.contains(&cell) {
                                *per_sample_cell_hist[sample_idx].entry(cell).or_insert(0) += count;
                            } else {
                                *per_sample_unmatched[sample_idx]
                                    .lock()
                                    .unwrap()
                                    .entry(cell)
                                    .or_insert(0) += count;
                            }
                        } else {
                            *per_sample_cell_hist[sample_idx].entry(cell).or_insert(0) += count;
                        }
                    }
                    CorrectionDecision::Ambiguous | CorrectionDecision::NotFound => {
                        deferred_rejected_sample_reads += count;
                        unmatched_reads.fetch_add(count, Ordering::SeqCst);
                    }
                }
                Ok(())
            })?;
            let stats = spool.stats();
            spool_stats.runs += stats.runs;
            spool_stats.pairs += stats.pairs;
            spool_stats.encoded_pairs += stats.encoded_pairs;
            spool_stats.uncompressed_bytes += stats.uncompressed_bytes;
            spool_stats.compressed_bytes += stats.compressed_bytes;
        }
        replay_time = replay_start.elapsed();
        info!(
            log,
            "Sample-Frequency replay: {} deferred reads in {} runs, {} compressed bytes from {} encoded bytes; replay took {:?}",
            spool_stats.pairs,
            spool_stats.runs,
            spool_stats.compressed_bytes,
            spool_stats.uncompressed_bytes,
            replay_time,
        );
    }

    let matched_reads = matched_reads.load(Ordering::SeqCst);
    let unmatched_reads = unmatched_reads.load(Ordering::SeqCst);

    // Unwrap the Arcs (single owner again after threads join)
    let per_sample_cell_hist = std::sync::Arc::into_inner(per_sample_cell_hist).unwrap();
    let per_sample_unmatched = std::sync::Arc::into_inner(per_sample_unmatched).unwrap();

    info!(
        log,
        "First pass complete: {} total reads, {} matched to samples, {} unmatched",
        total_reads.to_formatted_string(&Locale::en),
        matched_reads.to_formatted_string(&Locale::en),
        unmatched_reads.to_formatted_string(&Locale::en),
    );
    info!(
        log,
        "Sample routing: {} exact, {} structurally unique, {} immediately rejected, {} deferred ({} accepted and {} rejected after replay)",
        exact_sample_reads.load(Ordering::SeqCst),
        structurally_unique_sample_reads.load(Ordering::SeqCst),
        immediate_rejected_sample_reads,
        deferred_sample_reads.load(Ordering::SeqCst),
        deferred_accepted_sample_reads,
        deferred_rejected_sample_reads,
    );

    // Create output directory
    let parent = std::path::Path::new(output_dir);
    std::fs::create_dir_all(parent)
        .with_context(|| format!("couldn't create output directory {}", parent.display()))?;

    // Write sample_permit_map.bin
    let sample_map_path = parent.join("sample_permit_map.bin");
    let sample_map_file = File::create(&sample_map_path)?;
    bincode::serialize_into(BufWriter::new(sample_map_file), &sample_permit_map)
        .map_err(|e| anyhow!("failed to serialize sample permit map: {}", e))?;
    info!(
        log,
        "Wrote sample permit map to {}",
        sample_map_path.display()
    );

    // Per-sample: generate cell barcode permit lists
    let mut total_cells: u64 = 0;
    let mut sample_info_entries = Vec::new();
    let mut cell_correction_scopes = Vec::with_capacity(num_samples);

    // Get the cell barcode length from file tags
    let cell_bc_tag = format!("b{}len", num_barcodes - 1);
    let cell_bc_len: u16 = file_tag_map
        .get(&cell_bc_tag)
        .unwrap_or_else(|| panic!("expected '{}' file-level tag", cell_bc_tag))
        .try_into()
        .unwrap_or_else(|_| panic!("couldn't parse '{}' as u16", cell_bc_tag));
    let cell_spec = gpl_opts.cell_correction_spec(cell_bc_len as u8, true);
    warn_for_shift_frequency(cell_spec, "cell", log);
    let cell_correction_start = Instant::now();

    // Sample scopes are independent once exact priors are frozen. Use a small
    // worker cap: this removes the serial tail without multiplying the large
    // per-sample histograms by the user's full RAD-reader thread count.
    let correction_workers = gpl_opts.threads.min(4).min(num_samples).max(1);
    let mut partitions: Vec<Vec<_>> = (0..correction_workers).map(|_| Vec::new()).collect();
    for (sample_idx, input) in per_sample_cell_hist
        .into_iter()
        .zip(per_sample_unmatched)
        .enumerate()
    {
        partitions[sample_idx % correction_workers].push((sample_idx, input));
    }

    let canonical_barcodes = &sample_info.canonical_barcodes;
    let sample_names_ref = &sample_names;
    let filter_method = &gpl_opts.fmeth;
    let mut sample_outputs = std::thread::scope(|scope| -> anyhow::Result<Vec<_>> {
        let mut handles = Vec::with_capacity(correction_workers);
        for partition in partitions {
            let worker_log = log.clone();
            let handle = scope.spawn(move || -> anyhow::Result<Vec<_>> {
                partition
                    .into_iter()
                    .map(|(sample_idx, (cell_hist, unmatched_cells))| {
                        compile_multi_sample_cell_output(
                            sample_idx,
                            canonical_barcodes[sample_idx],
                            &sample_names_ref[sample_idx],
                            cell_hist,
                            unmatched_cells,
                            parent,
                            cell_spec,
                            filter_method,
                            min_reads,
                            cell_bc_len,
                            &worker_log,
                        )
                    })
                    .collect()
            });
            handles.push(handle);
        }

        let mut outputs = Vec::with_capacity(num_samples);
        for handle in handles {
            outputs.extend(handle.join().expect("cell-correction worker panicked")?);
        }
        Ok(outputs)
    })?;
    sample_outputs.sort_unstable_by_key(|output| output.sample_idx);
    for output in sample_outputs {
        total_cells += output.num_cells;
        sample_info_entries.push(output.info);
        cell_correction_scopes.push(output.scope);
    }

    let cell_correction_time = cell_correction_start.elapsed();
    let sample_corrections = sample_permit_map
        .iter()
        .map(|(&observed, &corrected)| BarcodeCorrection {
            observed,
            corrected,
        })
        .collect();
    let mut correction_plan = CorrectionPlan {
        sample_barcode_len: Some(sample_info.barcode_len as u8),
        cell_barcode_len: cell_bc_len as u8,
        sample_spec: sample_correction_spec,
        sample_corrections,
        cell_scopes: cell_correction_scopes,
    };
    let artifact_write_start = Instant::now();
    let correction_plan_path = correction_plan.write_to(parent)?;
    let artifact_write_time = artifact_write_start.elapsed();
    let correction_plan_bytes = std::fs::metadata(&correction_plan_path)?.len();
    info!(
        log,
        "Compiled cell corrections in {:?}; wrote {}-byte correction plan in {:?}",
        cell_correction_time,
        correction_plan_bytes,
        artifact_write_time,
    );

    // Write sample_info.json
    let resolved_sample_bc_correction = match gpl_opts.sample_correction_mode {
        crate::prog_opts::SampleCorrectionMode::Exact => "exact",
        crate::prog_opts::SampleCorrectionMode::Unique
        | crate::prog_opts::SampleCorrectionMode::OneEdit => "unique",
        crate::prog_opts::SampleCorrectionMode::Frequency => "frequency",
    };
    let resolved_sample_bc_neighborhood =
        sample_correction_spec.map(|spec| spec.neighborhood.to_string());
    let sample_info = serde_json::json!({
        "num_samples": num_samples,
        "num_barcodes": num_barcodes,
        "total_cells": total_cells,
        "total_reads": total_reads,
        "matched_reads": matched_reads,
        "unmatched_reads": unmatched_reads,
        "sample_correction_mode": format!("{:?}", gpl_opts.sample_correction_mode),
        "sample_bc_ori": format!("{:?}", gpl_opts.sample_bc_ori),
        "cell_bc_correction": gpl_opts.cell_bc_correction.to_string(),
        "cell_bc_neighborhood": cell_spec.neighborhood.to_string(),
        "cell_bc_confidence": gpl_opts.cell_bc_confidence.to_string(),
        "sample_bc_correction": resolved_sample_bc_correction,
        "sample_bc_neighborhood": resolved_sample_bc_neighborhood,
        "sample_bc_confidence": gpl_opts.sample_bc_confidence.to_string(),
        "cell_correction_seconds": cell_correction_time.as_secs_f64(),
        "correction_plan_bytes": correction_plan_bytes,
        "correction_plan_write_seconds": artifact_write_time.as_secs_f64(),
        "sample_routing": {
            "exact": exact_sample_reads.load(Ordering::SeqCst),
            "structurally_unique": structurally_unique_sample_reads.load(Ordering::SeqCst),
            "immediate_rejected": immediate_rejected_sample_reads,
            "deferred": deferred_sample_reads.load(Ordering::SeqCst),
            "deferred_accepted": deferred_accepted_sample_reads,
            "deferred_rejected": deferred_rejected_sample_reads,
            "spill_runs": spool_stats.runs,
            "spill_pairs_after_run_length_encoding": spool_stats.encoded_pairs,
            "spill_compressed_bytes": spool_stats.compressed_bytes,
            "spill_uncompressed_bytes": spool_stats.uncompressed_bytes,
            "replay_seconds": replay_time.as_secs_f64(),
        },
        "samples": sample_info_entries,
    });
    let info_path = parent.join("sample_info.json");
    let info_file = File::create(&info_path)?;
    serde_json::to_writer_pretty(BufWriter::new(info_file), &sample_info)?;

    // Write generate_permit_list.json (standard metadata)
    let gpl_meta = json!({
        "velo_mode": gpl_opts.velo_mode,
        "expected_ori": expected_ori.strand_symbol(),
        "version_str": gpl_opts.version,
        "cmd": gpl_opts.cmdline,
        "permit-list-type": format!("{:?}", gpl_opts.fmeth),
        "multi_barcode": true,
        "num_barcodes": num_barcodes,
        "cell_bc_correction": gpl_opts.cell_bc_correction.to_string(),
        "cell_bc_neighborhood": cell_spec.neighborhood.to_string(),
        "cell_bc_confidence": gpl_opts.cell_bc_confidence.to_string(),
        "sample_bc_correction": resolved_sample_bc_correction,
        "sample_bc_neighborhood": resolved_sample_bc_neighborhood,
        "sample_bc_confidence": gpl_opts.sample_bc_confidence.to_string(),
    });
    let gpl_meta_path = parent.join("generate_permit_list.json");
    let gpl_meta_file = File::create(&gpl_meta_path)?;
    serde_json::to_writer_pretty(BufWriter::new(gpl_meta_file), &gpl_meta)?;

    info!(
        log,
        "Multi-barcode permit list generation complete: {} samples, {} total cells",
        num_samples,
        total_cells,
    );

    Ok(total_cells)
}

/// Result of loading a sample barcode file with rotation/canonical structure.
struct SampleBarcodeInfo {
    /// The canonical (deduplicated) sample barcodes — one per true sample.
    canonical_barcodes: Vec<u64>,
    /// Maps every observed rotation barcode to its canonical barcode.
    rotation_to_canonical: HashMap<u64, u64>,
    /// Maps canonical barcode to sample name.
    canonical_to_name: HashMap<u64, String>,
    /// Number of nucleotides in each packed sample barcode.
    barcode_len: usize,
}

/// Load sample barcodes from a file.
///
/// Supports two formats:
/// 1. Simple: one barcode per line (each line is a separate sample)
/// 2. TSV with rotation mapping: `observed_bc  canonical_bc  sample_name`
///    Multiple observed barcodes can map to the same canonical/sample.
///    This is the standard 10x Flex probe barcode format where each sample
///    has 8 rotation variants.
///
/// Returns `SampleBarcodeInfo` with canonical barcodes and rotation mapping.
fn load_sample_barcode_list(
    path: &PathBuf,
    ori: SampleBarcodeOri,
    log: &slog::Logger,
) -> anyhow::Result<SampleBarcodeInfo> {
    let file = File::open(path)
        .with_context(|| format!("couldn't open sample barcode list: {}", path.display()))?;
    let reader = BufReader::new(file);

    if ori == SampleBarcodeOri::Reverse {
        info!(
            log,
            "Sample barcode whitelist orientation: reverse — reverse-complementing entries before lookup"
        );
    }

    let mut rotation_to_canonical: HashMap<u64, u64> = HashMap::new();
    let mut canonical_to_name: HashMap<u64, String> = HashMap::new();
    let mut canonical_order: Vec<u64> = Vec::new(); // preserves order, deduplicated
    let mut has_tsv_columns = false;
    let mut barcode_len = None;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = trimmed.split('\t').collect();

        let (observed_seq, canonical_seq, name) = if parts.len() >= 3 {
            // TSV format: observed_bc  canonical_bc  sample_name
            has_tsv_columns = true;
            (parts[0], parts[1], parts[2].to_string())
        } else if parts.len() == 2 {
            // Two-column: barcode  name (no rotation mapping)
            (parts[0], parts[0], parts[1].to_string())
        } else {
            // Single column: barcode only (each is its own sample)
            (parts[0], parts[0], parts[0].to_string())
        };

        if observed_seq.len() != canonical_seq.len() {
            bail!(
                "sample barcode '{}' and its canonical barcode '{}' have different lengths",
                observed_seq,
                canonical_seq
            );
        }
        if !(1..=32).contains(&observed_seq.len()) {
            bail!(
                "sample barcode length {} is unsupported; expected 1 through 32 bases",
                observed_seq.len()
            );
        }
        match barcode_len {
            Some(expected) if expected != observed_seq.len() => bail!(
                "sample barcode file mixes lengths {} and {}",
                expected,
                observed_seq.len()
            ),
            None => barcode_len = Some(observed_seq.len()),
            _ => {}
        }

        let rc = |seq: &str| -> String {
            seq.bytes()
                .rev()
                .map(|b| match b {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'a' => b't',
                    b't' => b'a',
                    b'c' => b'g',
                    b'g' => b'c',
                    other => other,
                })
                .map(|b| b as char)
                .collect()
        };

        let (observed_owned, canonical_owned);
        let (observed_seq, canonical_seq) = if ori == SampleBarcodeOri::Reverse {
            observed_owned = rc(observed_seq);
            canonical_owned = rc(canonical_seq);
            (observed_owned.as_str(), canonical_owned.as_str())
        } else {
            (observed_seq, canonical_seq)
        };

        let pack = |seq: &str| -> anyhow::Result<u64> {
            needletail::bitkmer::BitNuclKmer::new(seq.as_bytes(), seq.len() as u8, false)
                .next()
                .map(|(_, kmer, _)| kmer.0)
                .ok_or_else(|| anyhow!("couldn't pack sample barcode: {}", seq))
        };

        let obs_packed = pack(observed_seq)?;
        let canon_packed = pack(canonical_seq)?;

        match rotation_to_canonical.entry(obs_packed) {
            Entry::Occupied(entry) if *entry.get() != canon_packed => {
                bail!(
                    "sample construction barcode '{}' is assigned to conflicting canonical barcodes",
                    observed_seq
                );
            }
            Entry::Occupied(_) => {}
            Entry::Vacant(entry) => {
                entry.insert(canon_packed);
            }
        }
        match canonical_to_name.entry(canon_packed) {
            Entry::Occupied(entry) if entry.get() != &name => {
                bail!(
                    "canonical sample barcode '{}' is assigned conflicting names '{}' and '{}'",
                    canonical_seq,
                    entry.get(),
                    name
                );
            }
            Entry::Occupied(_) => {}
            Entry::Vacant(entry) => {
                canonical_order.push(canon_packed);
                entry.insert(name);
            }
        }
    }

    let num_rotations = rotation_to_canonical.len();
    let num_canonical = canonical_order.len();
    let barcode_len = barcode_len.ok_or_else(|| {
        anyhow!(
            "sample barcode list {} contains no barcodes",
            path.display()
        )
    })?;

    if has_tsv_columns && num_rotations > num_canonical {
        info!(
            log,
            "Loaded {} rotation barcodes mapping to {} canonical samples from {}",
            num_rotations,
            num_canonical,
            path.display(),
        );
    } else {
        info!(
            log,
            "Loaded {} sample barcodes from {}",
            num_canonical,
            path.display(),
        );
    }

    Ok(SampleBarcodeInfo {
        canonical_barcodes: canonical_order,
        rotation_to_canonical,
        canonical_to_name,
        barcode_len,
    })
}

/// Build the sample barcode permit map (observed -> canonical).
/// Returns (permit_map, bc_to_idx) where:
///   - permit_map maps every observed rotation barcode to its canonical barcode
///   - bc_to_idx maps canonical barcode to sample index
struct SampleCorrectionRouting {
    permit_map: HashMap<u64, u64>,
    canonical_to_index: HashMap<u64, usize>,
    ambiguous_sample_barcodes: Option<HashSet<u64, ahash::RandomState>>,
}

fn build_sample_permit_map(
    sample_info: &SampleBarcodeInfo,
    correction_spec: Option<CorrectionSpec>,
    log: &slog::Logger,
) -> anyhow::Result<SampleCorrectionRouting> {
    let mut permit_map: HashMap<u64, u64> = HashMap::new();
    let mut bc_to_idx: HashMap<u64, usize> = HashMap::new();

    // Map every rotation barcode to its canonical form
    for (&obs, &canon) in &sample_info.rotation_to_canonical {
        permit_map.insert(obs, canon);
    }

    // Index canonical barcodes
    for (idx, &canon) in sample_info.canonical_barcodes.iter().enumerate() {
        bc_to_idx.insert(canon, idx);
    }

    match correction_spec {
        None => {
            info!(
                log,
                "Sample barcode correction: exact match ({} rotation entries -> {} canonical samples)",
                permit_map.len(),
                bc_to_idx.len(),
            );
            Ok(SampleCorrectionRouting {
                permit_map,
                canonical_to_index: bc_to_idx,
                ambiguous_sample_barcodes: None,
            })
        }
        Some(spec) => {
            let topology_spec = CorrectionSpec {
                resolution: BarcodeResolution::Unique,
                ..spec
            };
            let topology = CorrectionIndex::new(
                topology_spec,
                sample_info
                    .rotation_to_canonical
                    .iter()
                    .map(|(&source, &target)| RetainedSource {
                        source,
                        target,
                        exact_count: 0,
                    }),
            )?;
            let structural = topology.compile_structural_neighborhood();
            permit_map.extend(structural.entries().iter().copied());
            info!(
                log,
                "Built deterministic sample correction topology with {} entries using {}",
                permit_map.len(),
                spec.neighborhood,
            );
            let ambiguous_sample_barcodes = if matches!(
                spec.resolution,
                BarcodeResolution::Frequency { .. }
            ) {
                let ambiguous = structural.ambiguous().iter().copied().collect();
                info!(
                    log,
                    "{} structurally ambiguous sample barcodes will be deferred until exact sample priors are frozen",
                    structural.ambiguous().len(),
                );
                Some(ambiguous)
            } else {
                info!(
                    log,
                    "Built sample permit map with {} entries ({} exact + {} corrected) -> {} canonical samples",
                    permit_map.len(),
                    sample_info.rotation_to_canonical.len(),
                    permit_map.len() - sample_info.rotation_to_canonical.len(),
                    bc_to_idx.len(),
                );
                None
            };
            Ok(SampleCorrectionRouting {
                permit_map,
                canonical_to_index: bc_to_idx,
                ambiguous_sample_barcodes,
            })
        }
    }
}

/// Get sample names for canonical barcodes.
/// Uses names from the barcode file, falling back to hex-encoded canonical barcode.
fn get_sample_names(sample_info: &SampleBarcodeInfo) -> Vec<String> {
    sample_info
        .canonical_barcodes
        .iter()
        .map(|&canon| {
            sample_info
                .canonical_to_name
                .get(&canon)
                .cloned()
                .unwrap_or_else(|| format!("{:x}", canon))
        })
        .collect()
}

/// update teh counts in the barcode histogram for those reads matching
/// the prescribed orientation.
pub fn update_barcode_hist_unfiltered<B, R>(
    hist: &DashMap<u64, u64, ahash::RandomState>,
    unmatched_bc: &mut HashMap<u64, u64, ahash::RandomState>,
    max_ambiguity_read: &mut usize,
    chunk: &chunk::Chunk<R>,
    expected_ori: &Strand,
) -> usize
where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B>,
{
    let mut num_strand_compat_reads = 0usize;
    for r in &chunk.reads {
        if r.has_alignment_on_strand(*expected_ori) {
            num_strand_compat_reads += 1;
            *max_ambiguity_read = r.num_aln().max(*max_ambiguity_read);
            // lookup the barcode in the map of unfiltered known
            // barcodes
            match hist.get_mut(&(r.collate_key().into())) {
                // if we find a match, increment the count
                Some(mut c) => *c += 1,
                // otherwise, push this into the unmatched list
                None => {
                    *unmatched_bc.entry(r.collate_key().into()).or_insert(0) += 1;
                }
            }
        }
    }

    /*
    match expected_ori {
        Strand::Unknown => {
            for r in &chunk.reads {
                num_strand_compat_reads += 1;
                *max_ambiguity_read = r.num_aln().max(*max_ambiguity_read);
                // lookup the barcode in the map of unfiltered known
                // barcodes
                match hist.get_mut(&(r.collate_key().into())) {
                    // if we find a match, increment the count
                    Some(mut c) => *c += 1,
                    // otherwise, push this into the unmatched list
                    None => {
                        *unmatched_bc.entry(r.collate_key().into()).or_insert(0) += 1;
                    }
                }
            }
        }
        Strand::Forward => {
            for r in &chunk.reads {
                if r.dirs.iter().any(|&x| x) {
                    num_strand_compat_reads += 1;
                    *max_ambiguity_read = r.num_aln().max(*max_ambiguity_read);
                    // lookup the barcode in the map of unfiltered known
                    // barcodes
                    match hist.get_mut(&(r.collate_key().into())) {
                        // if we find a match, increment the count
                        Some(mut c) => *c += 1,
                        // otherwise, push this into the unmatched list
                        None => {
                            *unmatched_bc.entry(r.collate_key().into()).or_insert(0) += 1;
                        }
                    }
                }
            }
        }
        Strand::Reverse => {
            for r in &chunk.reads {
                if r.dirs.iter().any(|&x| !x) {
                    num_strand_compat_reads += 1;
                    *max_ambiguity_read = r.num_aln().max(*max_ambiguity_read);
                    // lookup the barcode in the map of unfiltered known
                    // barcodes
                    match hist.get_mut(&(r.collate_key().into())) {
                        // if we find a match, increment the count
                        Some(mut c) => *c += 1,
                        // otherwise, push this into the unmatched list
                        None => {
                            *unmatched_bc.entry(r.collate_key().into()).or_insert(0) += 1;
                        }
                    }
                }
            }
        }
    }
    */
    num_strand_compat_reads
}

pub fn update_barcode_hist<B, R>(
    hist: &DashMap<u64, u64, ahash::RandomState>,
    max_ambiguity_read: &mut usize,
    chunk: &chunk::Chunk<R>,
    expected_ori: &Strand,
) -> usize
where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B>,
{
    let mut num_orientation_compat_reads = 0usize;
    for r in &chunk.reads {
        if r.has_alignment_on_strand(*expected_ori) {
            num_orientation_compat_reads += 1;
            *max_ambiguity_read = r.num_aln().max(*max_ambiguity_read);
            *hist.entry(r.collate_key().into()).or_insert(0) += 1;
        }
    }
    num_orientation_compat_reads
    /*
    match expected_ori {
        Strand::Unknown => {
            for r in &chunk.reads {
                *max_ambiguity_read = r.num_aln().max(*max_ambiguity_read);
                *hist.entry(r.collate_key().into()).or_insert(0) += 1;
            }
        }
        Strand::Forward => {
            for r in &chunk.reads {
                if r.dirs.iter().any(|&x| x) {
                    *max_ambiguity_read = r.num_aln().max(*max_ambiguity_read);
                    *hist.entry(r.collate_key().into()).or_insert(0) += 1;
                }
            }
        }
        Strand::Reverse => {
            for r in &chunk.reads {
                if r.dirs.iter().any(|&x| !x) {
                    *max_ambiguity_read = r.num_aln().max(*max_ambiguity_read);
                    *hist.entry(r.collate_key().into()).or_insert(0) += 1;
                }
            }
        }
    }
    */
}

pub fn permit_list_from_threshold(
    hist: &HashMap<u64, u64, ahash::RandomState>,
    min_freq: u64,
) -> Vec<u64> {
    hist.iter()
        .filter_map(|(k, v)| if v >= &min_freq { Some(*k) } else { None })
        .collect()
}

pub fn permit_list_from_file<P>(ifile: P, bclen: u16) -> Vec<u64>
where
    P: AsRef<Path>,
{
    let f = File::open(ifile).expect("couldn't open input barcode file.");
    let br = BufReader::new(f);
    let mut bc = Vec::<u64>::with_capacity(10_000);

    for l in br.lines() {
        let line = l.expect("couldn't read line from barcode file.");
        let mut bnk = BitNuclKmer::new(line.as_bytes(), bclen as u8, false);
        let (_, k, _) = bnk.next().expect("can't extract kmer");
        bc.push(k.0);
    }
    bc
}

struct UnfilteredBarcodeHist {
    pub unfiltered_bc_counts: Option<DashMap<u64, u64, ahash::RandomState>>,
    pub first_bclen: usize,
}

// Main entry point - now much cleaner
pub fn do_generate_permit_list<B, R>(
    gpl_opts: GenPermitListOpts,
    ifile: BufReader<File>,
    prelude: RadPrelude,
    file_tag_map: TagMap,
) -> anyhow::Result<u64>
where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize,
    <R as MappedRecord>::ParsingContext: RecordContext + Clone + Send,
{
    let UnfilteredBarcodeHist {
        unfiltered_bc_counts,
        first_bclen,
    } = load_unfiltered_barcodes(&gpl_opts, gpl_opts.log)?;

    if let Some(tag_val) = file_tag_map.get("cblen") {
        let cb_len: u64 = tag_val.try_into().expect("should be able to fit");
        if cb_len as usize != first_bclen {
            warn!(
                &gpl_opts.log,
                "The provided permit list had barcodes of length {}, but the mapped reads have barcodes of length {}",
                first_bclen,
                cb_len
            );
        }
    } else {
        warn!(
            &gpl_opts.log,
            "Expected \"cblen\" file-level tag, but it was not found."
        );
    }

    let rad_reader = setup_rad_reader::<R>(ifile, prelude, file_tag_map, gpl_opts.threads);

    log_rad_header_info(&rad_reader, gpl_opts.log);
    validate_tag_types(&rad_reader.prelude.read_tags, gpl_opts.log)?;

    let num_chunks = validate_chunks(&rad_reader.prelude.hdr, gpl_opts.log)?;
    let pbar = create_progress_bar(num_chunks);

    match &gpl_opts.fmeth {
        CellFilterMethod::UnfilteredExternalList(_, _) => {
            process_unfiltered_workflow(rad_reader, unfiltered_bc_counts, gpl_opts, pbar)
        }
        _ => process_filtered_workflow(rad_reader, gpl_opts, pbar),
    }
}

// Load and validate unfiltered barcode list if provided
fn load_unfiltered_barcodes(
    gpl_opts: &GenPermitListOpts,
    log: &slog::Logger,
) -> anyhow::Result<UnfilteredBarcodeHist> {
    let mut first_bclen = 0usize;
    let mut unfiltered_bc_counts = None;

    if let CellFilterMethod::UnfilteredExternalList(fname, _) = &gpl_opts.fmeth {
        let (reader, compression) = niffler::from_path(fname)
            .with_context(|| format!("could not open input file {}", fname.display()))?;
        let br = BufReader::new(reader);

        info!(
            log,
            "reading permit list from {}; inferred format {:#?}",
            fname.display(),
            compression
        );

        unfiltered_bc_counts = Some(populate_unfiltered_barcode_map(br, &mut first_bclen));
        info!(
            log,
            "number of unfiltered bcs read = {}",
            unfiltered_bc_counts
                .as_ref()
                .unwrap()
                .len()
                .to_formatted_string(&Locale::en)
        );
    }

    Ok(UnfilteredBarcodeHist {
        unfiltered_bc_counts,
        first_bclen,
    })
}

// Initialize the RAD reader with parallel processing
fn setup_rad_reader<R>(
    ifile: BufReader<File>,
    prelude: RadPrelude,
    file_tag_map: TagMap,
    nworkers: usize,
) -> libradicl::readers::ParallelRadReader<R, BufReader<File>>
where
    R: MappedRecord,
{
    libradicl::readers::ParallelRadReader::<R, _>::from_prelude_and_file_tag_map(
        ifile,
        prelude,
        file_tag_map,
        NonZeroUsize::new(nworkers).unwrap(),
    )
}

// Log information about the RAD file header
fn log_rad_header_info<R: MappedRecord, F: std::io::BufRead + std::io::Seek>(
    rad_reader: &libradicl::readers::ParallelRadReader<R, F>,
    log: &slog::Logger,
) {
    let hdr = &rad_reader.prelude.hdr;
    info!(
        log,
        "paired: {:?}, ref_count: {}, num_chunks: {}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en)
    );

    let fl_tags = &rad_reader.prelude.file_tags;
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());

    let rl_tags = &rad_reader.prelude.read_tags;
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());

    let al_tags = &rad_reader.prelude.aln_tags;
    info!(log, "read {:?} alignment-level tags", al_tags.tags.len());

    info!(log, "File-level tag values {:?}", &rad_reader.file_tag_map);
}

// Validate that barcode and UMI tags are present and of correct type
fn validate_tag_types(rl_tags: &TagSection, log: &slog::Logger) -> anyhow::Result<()> {
    const BNAME: &str = "b";
    const UNAME: &str = "u";

    let mut bct: Option<RadType> = None;
    let mut umit: Option<RadType> = None;

    for rt in &rl_tags.tags {
        if rt.name == BNAME || rt.name == UNAME {
            if !rt.typeid.is_int_type() {
                crit!(
                    log,
                    "currently only RAD types 1--4 are supported for 'b' and 'u' tags."
                );
                std::process::exit(exit_codes::EXIT_UNSUPPORTED_TAG_TYPE);
            }

            if rt.name == BNAME {
                bct = Some(rt.typeid);
            }
            if rt.name == UNAME {
                umit = Some(rt.typeid);
            }
        }
    }

    anyhow::ensure!(bct.is_some(), "barcode type tag must be present");
    anyhow::ensure!(umit.is_some(), "umi type tag must be present");

    Ok(())
}

// Validate that chunks are present in the RAD file
fn validate_chunks(hdr: &RadHeader, _log: &slog::Logger) -> anyhow::Result<u64> {
    hdr.num_chunks()
        .ok_or_else(|| anyhow!(
            "The RAD file appears to have no chunks; this most commonly occurs when no reads are mapped due to an incorrect chemistry being set. Please ensure that you have set the correct chemistry"
        ))
        .map(|nc| nc.get() as u64)
}

// Create progress bar for chunk processing
fn create_progress_bar(num_chunks: u64) -> ProgressBar {
    let pbar = ProgressBar::new(num_chunks);
    pbar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-"),
    );
    pbar.set_draw_target(ProgressDrawTarget::stderr_with_hz(5));
    pbar
}

// Process unfiltered barcode workflow
fn process_unfiltered_workflow<R, B>(
    mut rad_reader: libradicl::readers::ParallelRadReader<R, BufReader<File>>,
    unfiltered_bc_counts: Option<DashMap<u64, u64, ahash::RandomState>>,
    gpl_opts: GenPermitListOpts,
    pbar: ProgressBar,
) -> anyhow::Result<u64>
where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize,
    <R as MappedRecord>::ParsingContext: RecordContext + Clone + Send,
{
    let nworkers = gpl_opts.threads;
    let expected_ori = gpl_opts.expected_ori;

    let cb = create_progress_callback(pbar.clone());

    let Some(hmu) = unfiltered_bc_counts else {
        return Ok(0);
    };

    let hmu_arc = Arc::new(hmu);

    let (num_reads, num_orientation_compat_reads, unmatched_bc, max_ambiguity_read) =
        parse_chunks_unfiltered(
            &mut rad_reader,
            hmu_arc.clone(),
            nworkers,
            &expected_ori,
            cb,
        );

    pbar.finish_with_message("finished parsing RAD file\n");

    log_parsing_stats(
        gpl_opts.log,
        num_reads,
        num_orientation_compat_reads,
        rad_reader.prelude.hdr.num_chunks().unwrap(),
        max_ambiguity_read,
    );

    let unmatched_reads = unmatched_bc.values().sum::<u64>() as usize;
    check_permit_list_validity(gpl_opts.log, unmatched_reads, num_reads)?;

    let hmu = Arc::try_unwrap(hmu_arc).map_err(|_| anyhow!("Failed to unwrap Arc"))?;

    process_unfiltered(
        hmu,
        unmatched_bc,
        &rad_reader.file_tag_map,
        &gpl_opts.fmeth,
        expected_ori,
        gpl_opts.output_dir,
        gpl_opts.version,
        max_ambiguity_read,
        gpl_opts.velo_mode,
        gpl_opts.cmdline,
        gpl_opts.log,
        &gpl_opts,
    )
}

// Process filtered barcode workflow
fn process_filtered_workflow<R, B>(
    mut rad_reader: libradicl::readers::ParallelRadReader<R, BufReader<File>>,
    gpl_opts: GenPermitListOpts,
    pbar: ProgressBar,
) -> anyhow::Result<u64>
where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize,
    <R as MappedRecord>::ParsingContext: RecordContext + Clone + Send,
{
    let nworkers = gpl_opts.threads;
    let expected_ori = gpl_opts.expected_ori;

    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let hm = Arc::new(DashMap::with_hasher(s));

    let cb = create_progress_callback(pbar.clone());

    let (num_reads, num_orientation_compat_reads, max_ambiguity_read) =
        parse_chunks_filtered(&mut rad_reader, hm.clone(), nworkers, &expected_ori, cb);

    pbar.finish_with_message("finished parsing RAD file\n");

    info!(
        gpl_opts.log,
        "observed {} reads ({} orientation consistent) in {} chunks --- max ambiguity read occurs in {} refs",
        num_reads.to_formatted_string(&Locale::en),
        num_orientation_compat_reads.to_formatted_string(&Locale::en),
        rad_reader
            .prelude
            .hdr
            .num_chunks()
            .unwrap()
            .to_formatted_string(&Locale::en),
        max_ambiguity_read.to_formatted_string(&Locale::en)
    );

    let hm = Arc::into_inner(hm).ok_or_else(|| anyhow!("Failed to extract hash map"))?;

    process_filtered(
        hm,
        &rad_reader.file_tag_map,
        &gpl_opts.fmeth,
        expected_ori,
        gpl_opts.output_dir,
        gpl_opts.version,
        max_ambiguity_read,
        gpl_opts.velo_mode,
        gpl_opts.cmdline,
        gpl_opts.log,
        &gpl_opts,
    )
}

// Parse chunks with unfiltered barcode tracking
fn parse_chunks_unfiltered<R, B>(
    rad_reader: &mut libradicl::readers::ParallelRadReader<R, BufReader<File>>,
    hmu: Arc<DashMap<u64, u64, ahash::RandomState>>,
    nworkers: usize,
    expected_ori: &Strand,
    cb: impl Fn(u64, u64) + Send + Sync,
) -> (usize, usize, HashMap<u64, u64, ahash::RandomState>, usize)
where
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize,
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    <R as MappedRecord>::ParsingContext: RecordContext + Clone + Send,
{
    std::thread::scope(|s| {
        let mut handles = Vec::new();

        for _ in 0..nworkers {
            let chunks = rad_reader.chunk_iter();
            let hmu = hmu.clone();

            let handle = s.spawn(move || {
                let mut unmatched_bc: HashMap<u64, u64, ahash::RandomState> = HashMap::default();
                let mut max_ambiguity_read = 0usize;
                let mut num_reads = 0;
                let mut num_orientation_compat_reads = 0;

                for meta_chunk in chunks {
                    for c in meta_chunk.iter() {
                        num_orientation_compat_reads += update_barcode_hist_unfiltered(
                            &hmu,
                            &mut unmatched_bc,
                            &mut max_ambiguity_read,
                            &c,
                            expected_ori,
                        );
                        num_reads += c.reads.len();
                    }
                }
                (
                    num_reads,
                    num_orientation_compat_reads,
                    unmatched_bc,
                    max_ambiguity_read,
                )
            });
            handles.push(handle);
        }

        let _ = rad_reader.start_chunk_parsing(Some(cb));

        let mut total_reads = 0;
        let mut total_compat = 0;
        let mut all_unmatched: HashMap<u64, u64, ahash::RandomState> = HashMap::default();
        let mut max_ambig = 0;

        for handle in handles {
            let (nr, nocr, ubc, mar) = handle.join().expect("The parsing thread panicked");
            total_reads += nr;
            total_compat += nocr;
            for (barcode, count) in ubc {
                *all_unmatched.entry(barcode).or_insert(0) += count;
            }
            max_ambig = max_ambig.max(mar);
        }

        (total_reads, total_compat, all_unmatched, max_ambig)
    })
}

// Parse chunks with filtered barcode tracking
fn parse_chunks_filtered<R, B>(
    rad_reader: &mut libradicl::readers::ParallelRadReader<R, BufReader<File>>,
    hm: Arc<DashMap<u64, u64, ahash::RandomState>>,
    nworkers: usize,
    expected_ori: &Strand,
    cb: impl Fn(u64, u64) + Send + Sync,
) -> (usize, usize, usize)
where
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize,
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    <R as MappedRecord>::ParsingContext: RecordContext + Clone + Send,
{
    std::thread::scope(|s| {
        let mut handles = Vec::new();

        for _ in 0..nworkers {
            let chunks = rad_reader.chunk_iter();
            let hm = hm.clone();

            let handle = s.spawn(move || {
                let mut max_ambiguity_read = 0usize;
                let mut num_reads = 0;
                let mut num_orientation_compat_reads = 0;

                for meta_chunk in chunks {
                    for c in meta_chunk.iter() {
                        num_orientation_compat_reads +=
                            update_barcode_hist(&hm, &mut max_ambiguity_read, &c, expected_ori);
                        num_reads += c.reads.len();
                    }
                }
                (num_reads, num_orientation_compat_reads, max_ambiguity_read)
            });
            handles.push(handle);
        }

        let _ = rad_reader.start_chunk_parsing(Some(cb));

        let mut total_reads = 0;
        let mut total_compat = 0;
        let mut max_ambig = 0;

        for handle in handles {
            let (nr, nocr, mar) = handle.join().expect("The parsing thread panicked");
            total_reads += nr;
            total_compat += nocr;
            max_ambig = max_ambig.max(mar);
        }

        (total_reads, total_compat, max_ambig)
    })
}

// Create callback for progress bar updates
fn create_progress_callback(pbar: ProgressBar) -> impl Fn(u64, u64) + Send + Sync {
    move |_new_bytes: u64, new_chunks: u64| {
        pbar.inc(new_chunks);
    }
}

// Log parsing statistics
fn log_parsing_stats(
    log: &slog::Logger,
    num_reads: usize,
    num_orientation_compat_reads: usize,
    num_chunks: NonZeroUsize,
    max_ambiguity_read: usize,
) {
    info!(
        log,
        "observed {} reads ({} orientation consistent) in {} chunks --- max ambiguity read occurs in {} refs",
        num_reads.to_formatted_string(&Locale::en),
        num_orientation_compat_reads.to_formatted_string(&Locale::en),
        num_chunks.to_formatted_string(&Locale::en),
        max_ambiguity_read.to_formatted_string(&Locale::en)
    );
}

// Check if permit list appears valid
fn check_permit_list_validity(
    log: &slog::Logger,
    num_unmatched: usize,
    num_reads: usize,
) -> anyhow::Result<()> {
    let valid_thresh = 0.3f64;
    match diagnostics::likely_valid_permit_list(num_unmatched, num_reads, valid_thresh) {
        Ok(f) => {
            info!(
                log,
                "The percentage of mapped reads not matching a known barcode exactly is {:.3}%, which is < the warning threshold {:.3}%",
                f * 100f64,
                valid_thresh * 100f64
            );
            Ok(())
        }
        Err(e) => {
            warn!(log, "{:?}", e);
            Ok(())
        }
    }
}

/*
/// Dispatched by `generate_permit_list` with the appropriate generic type for the
/// read record.
pub fn do_generate_permit_list<B, R>(
    gpl_opts: GenPermitListOpts,
    ifile: BufReader<File>,
    prelude: RadPrelude,
    file_tag_map: TagMap) -> anyhow::Result<u64>
where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize,
       <R as MappedRecord>::ParsingContext: RecordContext,
       <R as MappedRecord>::ParsingContext: Clone,
       <R as MappedRecord>::ParsingContext: Send
{
    let output_dir = gpl_opts.output_dir;
    let filter_meth = gpl_opts.fmeth.clone();
    let expected_ori = gpl_opts.expected_ori;
    let version = gpl_opts.version;
    let velo_mode = gpl_opts.velo_mode;
    let cmdline = gpl_opts.cmdline;
    let log = gpl_opts.log;

    let mut first_bclen = 0usize;
    let mut unfiltered_bc_counts = None;
    if let CellFilterMethod::UnfilteredExternalList(fname, _) = &filter_meth {
        let (reader, compression) = niffler::from_path(fname)
            .with_context(|| format!("coult not open input file {}", fname.display()))?;
        let br = BufReader::new(reader);

        info!(
            log,
            "reading permit list from {}; inferred format {:#?}",
            fname.display(),
            compression
        );

        unfiltered_bc_counts = Some(populate_unfiltered_barcode_map(br, &mut first_bclen));
        info!(
            log,
            "number of unfiltered bcs read = {}",
            unfiltered_bc_counts
                .as_ref()
                .unwrap()
                .len()
                .to_formatted_string(&Locale::en)
        );
    }

    let nworkers: usize = gpl_opts.threads;

    let mut rad_reader = libradicl::readers::ParallelRadReader::<R, _,>::from_prelude_and_file_tag_map(ifile, prelude, file_tag_map, NonZeroUsize::new(nworkers).unwrap());

    let hdr = &rad_reader.prelude.hdr;
    info!(
        log,
        "paired : {:?}, ref_count : {}, num_chunks : {}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en)
    );
    let num_chunks = hdr.num_chunks();

    // file-level
    let fl_tags = &rad_reader.prelude.file_tags;
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = &rad_reader.prelude.read_tags;
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());

    // right now, we only handle BC and UMI types of U8—U64, so validate that
    const BNAME: &str = "b";
    const UNAME: &str = "u";

    let mut bct: Option<RadType> = None;
    let mut umit: Option<RadType> = None;

    for rt in &rl_tags.tags {
        // if this is one of our tags
        if rt.name == BNAME || rt.name == UNAME {
            if !rt.typeid.is_int_type() {
                crit!(
                    log,
                    "currently only RAD types 1--4 are supported for 'b' and 'u' tags."
                );
                std::process::exit(exit_codes::EXIT_UNSUPPORTED_TAG_TYPE);
            }

            if rt.name == BNAME {
                bct = Some(rt.typeid);
            }
            if rt.name == UNAME {
                umit = Some(rt.typeid);
            }
        }
    }
    assert!(bct.is_some(), "barcode type tag must be present.");
    assert!(umit.is_some(), "umi type tag must be present.");

    // alignment-level
    let al_tags = &rad_reader.prelude.aln_tags;
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    {
        let file_tag_map = &rad_reader.file_tag_map;
        info!(log, "File-level tag values {:?}", file_tag_map);
    }

    let mut num_reads: usize = 0;

    // if dealing with filtered type
    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let hm = std::sync::Arc::new(DashMap::with_hasher(s));

    // if dealing with the unfiltered type
    // the set of barcodes that are not an exact match for any known barcodes
    let mut unmatched_bc: Vec<u64>;
    let mut num_orientation_compat_reads = 0usize;
    let mut max_ambiguity_read = 0usize;

    let nc = num_chunks.ok_or(
        anyhow!("The RAD file appears to have no chunks; this most commonly occurs when no reads are mapped due to an incorrect chemistry being set. Please ensure that you have set the correct chemistry"))?
        .get() as u64;
    let pbar = ProgressBar::new(nc);
    pbar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-"),
    );
    pbar.set_draw_target(ProgressDrawTarget::stderr_with_hz(5));
    let cb = |_new_bytes: u64, new_chunks: u64| {
        pbar.inc(new_chunks);
    };

    match filter_meth {
        CellFilterMethod::UnfilteredExternalList(_, _min_reads) => {
            unmatched_bc = Vec::with_capacity(10_000_000);
            // the unfiltered_bc_count map must be valid in this branch
            if unfiltered_bc_counts.is_some() {
                let hmu = std::thread::scope(|s| {
                    let hmu = std::sync::Arc::new(unfiltered_bc_counts.unwrap());
                    let mut handles = Vec::<
                        std::thread::ScopedJoinHandle<(usize, usize, Vec<u64>, usize)>,
                    >::new();
                    for _ in 0..nworkers {
                        let rd = rad_reader.is_done();
                        let q = rad_reader.get_queue();
                        let hmu = hmu.clone();
                        let handle = s.spawn(move || {
                            let mut unmatched_bc = Vec::<u64>::new();
                            let mut max_ambiguity_read = 0usize;
                            let mut num_reads = 0;
                            let mut num_orientation_compat_reads = 0;
                            while !rd.load(Ordering::SeqCst) || !q.is_empty() {
                                while let Some(meta_chunk) = q.pop() {
                                    for c in meta_chunk.iter() {
                                        num_orientation_compat_reads +=
                                            update_barcode_hist_unfiltered(
                                                &hmu,
                                                &mut unmatched_bc,
                                                &mut max_ambiguity_read,
                                                &c,
                                                &expected_ori,
                                            );
                                        num_reads += c.reads.len();
                                    }
                                }
                            }
                            (
                                num_reads,
                                num_orientation_compat_reads,
                                unmatched_bc,
                                max_ambiguity_read,
                            )
                        });
                        handles.push(handle);
                    }
                    let _ = rad_reader.start_chunk_parsing(Some(cb)); //libradicl::readers::EMPTY_METACHUNK_CALLBACK);
                    for handle in handles {
                        let (nr, nocr, ubc, mar) =
                            handle.join().expect("The parsing thread panicked");
                        num_reads += nr;
                        num_orientation_compat_reads += nocr;
                        unmatched_bc.extend_from_slice(&ubc);
                        max_ambiguity_read = max_ambiguity_read.max(mar);
                    }
                    pbar.finish_with_message("finished parsing RAD file\n");
                    // return the hash map we no longer need
                    std::sync::Arc::<DashMap<u64, u64, ahash::RandomState>>::into_inner(hmu)
                });
                /*
                for _ in 0..(hdr.num_chunks as usize) {
                let c =
                chunk::Chunk::<AlevinFryReadRecord>::from_bytes(&mut br, &record_context);
                num_orientation_compat_reads += update_barcode_hist_unfiltered(
                &mut hmu,
                &mut unmatched_bc,
                &mut max_ambiguity_read,
                &c,
                &expected_ori,
                );
                num_reads += c.reads.len();
                }
                */
                info!(
                        log,
                        "observed {} reads ({} orientation consistent) in {} chunks --- max ambiguity read occurs in {} refs",
                        num_reads.to_formatted_string(&Locale::en),
                        num_orientation_compat_reads.to_formatted_string(&Locale::en),
                        num_chunks.expect("nonzero").to_formatted_string(&Locale::en),
                        max_ambiguity_read.to_formatted_string(&Locale::en)
                    );
                let valid_thresh = 0.3f64;
                match diagnostics::likely_valid_permit_list(
                    unmatched_bc.len(),
                    num_reads,
                    valid_thresh,
                ) {
                    Ok(f) => {
                        info!(log,
                        "The percentage of mapped reads not matching a known barcode exactly is {:.3}%, which is < the warning threshold {:.3}%",
                        f * 100f64, valid_thresh * 100f64);
                    }
                    Err(e) => {
                        warn!(log, "{:?}", e);
                    }
                }

                process_unfiltered(
                    hmu.unwrap(),
                    unmatched_bc,
                    &rad_reader.file_tag_map,
                    &filter_meth,
                    expected_ori,
                    output_dir,
                    version,
                    max_ambiguity_read,
                    velo_mode,
                    cmdline,
                    log,
                    &gpl_opts,
                )
            } else {
                Ok(0)
            }
        }
        _ => {
            let hm = std::thread::scope(|s| {
                let mut handles =
                    Vec::<std::thread::ScopedJoinHandle<(usize, usize, usize)>>::new();
                for _ in 0..nworkers {
                    let rd = rad_reader.is_done();
                    let q = rad_reader.get_queue();
                    let hm = hm.clone();
                    let handle = s.spawn(move || {
                        let mut max_ambiguity_read = 0usize;
                        let mut num_reads = 0;
                        while !rd.load(Ordering::SeqCst) || !q.is_empty() {
                            while let Some(meta_chunk) = q.pop() {
                                for c in meta_chunk.iter() {
                                    update_barcode_hist(
                                        &hm,
                                        &mut max_ambiguity_read,
                                        &c,
                                        &expected_ori,
                                    );
                                    num_reads += c.reads.len();
                                }
                            }
                        }
                        (num_reads, num_orientation_compat_reads, max_ambiguity_read)
                    });
                    handles.push(handle);
                }
                let _ = rad_reader.start_chunk_parsing(Some(cb));
                for handle in handles {
                    let (nr, nocr, mar) = handle.join().expect("The parsing thread panicked");
                    num_reads += nr;
                    num_orientation_compat_reads += nocr;
                    max_ambiguity_read = max_ambiguity_read.max(mar);
                }
                pbar.finish_with_message("finished parsing RAD file\n");
                // return the hash map we no longer need
                Arc::<DashMap<u64, u64, ahash::RandomState>>::into_inner(hm)
                    .expect("unique reference to DashMap")
            });
            info!(
                log,
                "observed {} reads in {} chunks --- max ambiguity read occurs in {} refs",
                num_reads.to_formatted_string(&Locale::en),
                num_chunks.unwrap().to_formatted_string(&Locale::en),
                max_ambiguity_read.to_formatted_string(&Locale::en)
            );
            process_filtered(
                hm,
                &rad_reader.file_tag_map,
                &filter_meth,
                expected_ori,
                output_dir,
                version,
                max_ambiguity_read,
                velo_mode,
                cmdline,
                log,
                &gpl_opts,
            )
        }
    }
}
*/

#[cfg(test)]
mod cell_barcode_correction_tests {
    use super::*;
    use crate::prog_opts::CellBarcodeCorrectionStrategy;

    #[test]
    fn correction_strategy_has_stable_metadata_names_and_default() {
        assert_eq!(
            CellBarcodeCorrectionStrategy::default(),
            CellBarcodeCorrectionStrategy::Unique
        );
        assert_eq!(
            serde_json::to_string(&CellBarcodeCorrectionStrategy::Frequency).unwrap(),
            "\"frequency\""
        );
    }

    #[test]
    fn immediate_sample_routes_fuse_target_index_and_exactness() {
        let permit_map = HashMap::from([(10, 100), (11, 100), (20, 200)]);
        let canonical_to_index = HashMap::from([(100, 0), (200, 1)]);
        let exact_sources = HashMap::from([(10, 100), (20, 200)]);

        let routes =
            compile_immediate_sample_routes(&permit_map, &canonical_to_index, &exact_sources)
                .unwrap();
        assert_eq!(
            routes.get(&10),
            Some(&ImmediateSampleRoute {
                sample_index: 0,
                exact_source: true,
            })
        );
        assert_eq!(
            routes.get(&11),
            Some(&ImmediateSampleRoute {
                sample_index: 0,
                exact_source: false,
            })
        );
        assert!(
            compile_immediate_sample_routes(
                &HashMap::from([(30, 300)]),
                &canonical_to_index,
                &exact_sources,
            )
            .is_err()
        );
    }

    #[test]
    fn buffered_cell_counts_reduce_exactly_and_clear_local_state() {
        let mut local = vec![BarcodeCountMap::default(), BarcodeCountMap::default()];
        let mut buffered_distinct = 0;
        buffer_count(&mut local[0], 7, 2, &mut buffered_distinct);
        buffer_count(&mut local[0], 7, 3, &mut buffered_distinct);
        buffer_count(&mut local[1], 9, 4, &mut buffered_distinct);
        assert_eq!(buffered_distinct, 2);

        let shared = vec![DashMap::default(), DashMap::default()];
        flush_cell_counts(&mut local, &shared);
        assert!(local.iter().all(HashMap::is_empty));
        assert_eq!(shared[0].get(&7).map(|count| *count), Some(5));
        assert_eq!(shared[1].get(&9).map(|count| *count), Some(4));

        buffer_count(&mut local[0], 7, 6, &mut buffered_distinct);
        flush_cell_counts(&mut local, &shared);
        assert_eq!(shared[0].get(&7).map(|count| *count), Some(11));
    }

    #[test]
    fn sample_one_edit_correction_drops_cross_sample_collisions() {
        let sample_info = SampleBarcodeInfo {
            canonical_barcodes: vec![0, 5],
            rotation_to_canonical: HashMap::from([(0, 0), (5, 5)]),
            canonical_to_name: HashMap::new(),
            barcode_len: 2,
        };
        let log = slog::Logger::root(slog::Discard, slog::o!());
        let spec = CorrectionSpec {
            barcode_len: 2,
            neighborhood: crate::barcode_correction::BarcodeNeighborhood::SubstitutionOrShiftOne,
            resolution: BarcodeResolution::Unique,
        };
        let routing = build_sample_permit_map(&sample_info, Some(spec), &log).unwrap();
        let permit_map = routing.permit_map;
        let sample_indices = routing.canonical_to_index;

        assert_eq!(permit_map.get(&0), Some(&0));
        assert_eq!(permit_map.get(&5), Some(&5));
        assert!(!permit_map.contains_key(&1));
        assert_eq!(sample_indices, HashMap::from([(0, 0), (5, 1)]));

        let rotations_of_one_sample = SampleBarcodeInfo {
            canonical_barcodes: vec![9],
            rotation_to_canonical: HashMap::from([(0, 9), (5, 9)]),
            canonical_to_name: HashMap::new(),
            barcode_len: 2,
        };
        let permit_map = build_sample_permit_map(&rotations_of_one_sample, Some(spec), &log)
            .unwrap()
            .permit_map;
        assert_eq!(permit_map.get(&1), Some(&9));
    }

    #[test]
    fn sample_corrections_are_whitelist_order_invariant_and_conflicts_fail() {
        let dir = tempfile::tempdir().unwrap();
        let forward_path = dir.path().join("forward.tsv");
        let reverse_path = dir.path().join("reverse.tsv");
        std::fs::write(&forward_path, "AA\tAA\ta\nAC\tAA\ta\nCC\tCC\tb\n").unwrap();
        std::fs::write(&reverse_path, "CC\tCC\tb\nAC\tAA\ta\nAA\tAA\ta\n").unwrap();
        let log = slog::Logger::root(slog::Discard, slog::o!());
        let spec = CorrectionSpec {
            barcode_len: 2,
            neighborhood: crate::barcode_correction::BarcodeNeighborhood::HammingOne,
            resolution: BarcodeResolution::Unique,
        };
        let forward =
            load_sample_barcode_list(&forward_path, SampleBarcodeOri::Forward, &log).unwrap();
        let reverse =
            load_sample_barcode_list(&reverse_path, SampleBarcodeOri::Forward, &log).unwrap();
        assert_eq!(
            build_sample_permit_map(&forward, Some(spec), &log)
                .unwrap()
                .permit_map,
            build_sample_permit_map(&reverse, Some(spec), &log)
                .unwrap()
                .permit_map
        );

        let conflict_path = dir.path().join("conflict.tsv");
        std::fs::write(&conflict_path, "AA\tAA\ta\nAA\tCC\tb\n").unwrap();
        assert!(load_sample_barcode_list(&conflict_path, SampleBarcodeOri::Forward, &log).is_err());
    }
}
