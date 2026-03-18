/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use crate::diagnostics;
use crate::knee_finding;
use crate::prog_opts::GenPermitListOpts;
use crate::utils as afutils;
use crate::utils::KnownRecordType;
#[allow(unused_imports)]
use ahash::{AHasher, RandomState};
use anyhow::{Context, anyhow};
use bio_types::strand::Strand;
use bstr::io::BufReadExt;
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::Itertools;
use libradicl::BarcodeLookupMap;
use libradicl::exit_codes;
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{self, RadType, TagMap, TagSection};
use libradicl::record::{
    AlevinFryReadRecordWithPosition, CollatableMappedRecord, ConvertiblePrimitiveInteger,
    HierarchicallyCollatable, KnownSize, MappedRecord, MultiBarcodeReadRecord,
    MultiBarcodeRecordContext, RecordContext, ScLongReadRecord,
};
use libradicl::{chunk, record::AlevinFryReadRecord};
use crate::prog_opts::SampleCorrectionMode;
use needletail::bitkmer::*;
use num_format::{Locale, ToFormattedString};
use serde::Serialize;
use serde_json::json;
use slog::crit;
use slog::{info, warn};
use std::cmp;
use std::collections::HashMap;
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

#[allow(clippy::unnecessary_unwrap, clippy::too_many_arguments)]
fn process_unfiltered(
    hm: DashMap<u64, u64, ahash::RandomState>,
    mut unmatched_bc: Vec<u64>,
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

    // the smallest number of reads we'll allow per barcode
    let min_freq = match filter_meth {
        CellFilterMethod::UnfilteredExternalList(_, min_reads) => {
            info!(log, "minimum num reads for barcode pass = {}", *min_reads);
            *min_reads as u64
        }
        _ => {
            unimplemented!();
        }
    };

    // the set of barcodes we'll keep
    let mut kept_bc = Vec::<u64>::new();

    // iterate over the count map
    for mut kvp in hm.iter_mut() {
        // if this satisfies our requirement for the minimum count
        // then keep this barcode
        if *kvp.value() >= min_freq {
            kept_bc.push(*kvp.key());
        } else {
            // otherwise, we have to add this barcode's
            // counts to our unmatched list
            for _ in 0..*kvp.value() {
                unmatched_bc.push(*kvp.key());
            }
            // and then reset the counter for this barcode to 0
            *kvp.value_mut() = 0u64;
        }
    }

    // drop the absent barcodes from hm
    hm.retain(|_, &mut v| v > 0);

    // how many we will keep
    let num_passing = kept_bc.len();
    info!(
        log,
        "num_passing = {}",
        num_passing.to_formatted_string(&Locale::en)
    );

    let barcode_tag = file_tag_map
        .get("cblen")
        .expect("tag map must contain cblen");
    let barcode_len: u16 = barcode_tag.try_into()?;

    // now, we create a second barcode map with just the barcodes
    // for cells we will keep / rescue.
    let bcmap2 = BarcodeLookupMap::new(kept_bc, barcode_len as u32);
    info!(
        log,
        "found {} cells with non-trivial number of reads by exact barcode match",
        bcmap2.barcodes.len().to_formatted_string(&Locale::en)
    );

    // finally, we'll go through the set of unmatched barcodes
    // and try to rescue those that have a *unique* neighbor in the
    // list of retained barcodes.

    //let mut found_exact = 0usize;
    let mut found_approx = 0usize;
    let mut ambig_approx = 0usize;
    let mut not_found = 0usize;

    let start_unmatched_time = Instant::now();

    unmatched_bc.sort_unstable();

    let mut distinct_unmatched_bc = 0usize;
    let mut distinct_recoverable_bc = 0usize;

    // mapping the uncorrected barcode to what it corrects to
    let mut corrected_list = Vec::<(u64, u64)>::with_capacity(1_000_000);

    for (count, ubc) in unmatched_bc.iter().dedup_with_count() {
        // try to find the unmatched barcode, but
        // look up to 1 edit away
        match bcmap2.find_neighbors(*ubc, false) {
            // if we have a match
            (Some(x), n) => {
                let cbc = bcmap2.barcodes[x];
                // if the uncorrected barcode had a
                // single, unique retained neighbor
                if cbc != *ubc && n == 1 {
                    // then increment the count of this
                    // barcode by 1 (because we'll correct to it)
                    if let Some(mut c) = hm.get_mut(&cbc) {
                        *c += count as u64;
                        corrected_list.push((*ubc, cbc));
                    }
                    // this counts as an approximate find
                    found_approx += count;
                    distinct_recoverable_bc += 1;
                }
                // if we had > 1 single-mismatch neighbor
                // then don't keep the barcode, but remember
                // the count of such events
                if n > 1 {
                    ambig_approx += count;
                }
            }
            // if we had no single-mismatch neighbor
            // then this barcode is not_found and gets
            // dropped.
            (None, _) => {
                not_found += count;
            }
        }
        distinct_unmatched_bc += 1;
    }
    let unmatched_duration = start_unmatched_time.elapsed();
    let num_corrected = distinct_recoverable_bc as u64;

    info!(
        log,
        "There were {} distinct unmatched barcodes, and {} that can be recovered",
        distinct_unmatched_bc.to_formatted_string(&Locale::en),
        distinct_recoverable_bc.to_formatted_string(&Locale::en)
    );
    info!(
        log,
        "Matching unmatched barcodes to retained barcodes took {:?}", unmatched_duration
    );
    info!(log, "Of the unmatched barcodes\n============");
    info!(
        log,
        "\t{} had exactly 1 single-edit neighbor in the retained list",
        found_approx.to_formatted_string(&Locale::en)
    );
    info!(
        log,
        "\t{} had >1 single-edit neighbor in the retained list",
        ambig_approx.to_formatted_string(&Locale::en)
    );
    info!(
        log,
        "\t{} had no neighbor in the retained list",
        not_found.to_formatted_string(&Locale::en)
    );

    let parent = std::path::Path::new(output_dir);
    std::fs::create_dir_all(parent).with_context(|| {
        format!(
            "couldn't create path to output directory {}",
            parent.display()
        )
    })?;
    let o_path = parent.join("permit_freq.bin");

    // convert the DashMap to a HashMap
    let mut hm: HashMap<u64, u64, ahash::RandomState> = hm.into_iter().collect();
    match afutils::write_permit_list_freq(&o_path, barcode_len, &hm) {
        Ok(_) => {}
        Err(error) => {
            panic!("Error: {}", error);
        }
    };

    /*
    // don't need this right now
    let s_path = parent.join("bcmap.bin");
    let s_file = std::fs::File::create(&s_path).expect("could not create serialization file.");
    let mut s_writer = BufWriter::new(&s_file);
    bincode::serialize_into(&mut s_writer, &bcmap2).expect("couldn't serialize barcode list.");
    */

    // now that we are done with using hm to count, we can repurpose it as
    // the correction map.
    for (k, v) in hm.iter_mut() {
        // each present barcode corrects to itself
        *v = *k;
        //*kvp.value_mut() = *kvp.key();
    }
    for (uncorrected, corrected) in corrected_list.iter() {
        hm.insert(*uncorrected, *corrected);
    }

    let pm_path = parent.join("permit_map.bin");
    let pm_file = std::fs::File::create(pm_path).context("could not create serialization file.")?;
    let mut pm_writer = BufWriter::new(&pm_file);
    bincode::serialize_into(&mut pm_writer, &hm)
        .context("couldn't serialize permit list mapping.")?;

    let meta_info = json!({
    "velo_mode" : velo_mode,
    "expected_ori" : *expected_ori.strand_symbol(),
    "version_str" : version,
    "max-ambig-record" : max_ambiguity_read,
    "cmd" : cmdline,
    "permit-list-type" : "unfiltered",
    "gpl_options" : &gpl_opts
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
        num_corrected.to_formatted_string(&Locale::en)
    );

    Ok(num_corrected)
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
    let valid_bc: Vec<u64>;
    let hm: HashMap<u64, u64, ahash::RandomState> = hm.into_iter().collect();
    let mut freq: Vec<u64> = hm.values().cloned().collect();
    freq.sort_unstable();
    freq.reverse();

    let barcode_tag = file_tag_map
        .get("cblen")
        .expect("tag map must contain cblen");
    let barcode_len: u16 = barcode_tag.try_into()?;

    // select from among supported filter methods
    match filter_meth {
        CellFilterMethod::KneeFinding => {
            let num_bc = knee_finding::get_knee(&freq[..], 100, log);
            let min_freq = freq[num_bc];

            // collect all of the barcodes that have a frequency
            // >= to min_thresh.
            valid_bc = permit_list_from_threshold(&hm, min_freq);
            info!(
                log,
                "knee distance method resulted in the selection of {} permitted barcodes.",
                valid_bc.len()
            );
        }
        CellFilterMethod::ForceCells(top_k) => {
            let num_bc = if freq.len() < *top_k {
                freq.len() - 1
            } else {
                top_k - 1
            };

            let min_freq = freq[num_bc];

            // collect all of the barcodes that have a frequency
            // >= to min_thresh.
            valid_bc = permit_list_from_threshold(&hm, min_freq);
        }
        CellFilterMethod::ExplicitList(valid_bc_file) => {
            valid_bc = permit_list_from_file(valid_bc_file, barcode_len);
        }
        CellFilterMethod::ExpectCells(expected_num_cells) => {
            let robust_quantile = 0.99f64;
            let robust_div = 10.0f64;
            let robust_ind = (*expected_num_cells as f64 * robust_quantile).round() as u64;
            // the robust ind must be valid
            let ind = cmp::min(freq.len() - 1, robust_ind as usize);
            let robust_freq = freq[ind];
            let min_freq = std::cmp::max(1u64, (robust_freq as f64 / robust_div).round() as u64);
            valid_bc = permit_list_from_threshold(&hm, min_freq);
        }
        CellFilterMethod::UnfilteredExternalList(_, _min_reads) => {
            unimplemented!();
        }
    }

    // generate the map from each permitted barcode to all barcodes within
    // edit distance 1 of it.
    let full_permit_list =
        afutils::generate_permitlist_map(&valid_bc, barcode_len as usize).unwrap();

    let s2 = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut permitted_map = HashMap::with_capacity_and_hasher(valid_bc.len(), s2);

    let mut num_corrected = 0;
    for (k, v) in hm.iter() {
        if let Some(&valid_key) = full_permit_list.get(k) {
            *permitted_map.entry(valid_key).or_insert(0u64) += *v;
            num_corrected += 1;
            //println!("{} was a neighbor of {}, with count {}", k, valid_key, v);
        }
    }

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

    let meta_info = json!({
    "velo_mode" : velo_mode,
    "expected_ori" : *expected_ori.strand_symbol(),
    "version_str" : version,
    "max-ambig-record" : max_ambiguity_read,
    "cmd" : cmdline,
    "permit-list-type" : "filtered",
    "gpl_options" : &gpl_opts
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
        num_corrected.to_formatted_string(&Locale::en)
    );

    Ok(num_corrected)
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
            do_generate_permit_list_multi_bc(
                gpl_opts,
                ifile,
                prelude,
                file_tag_map,
                num_bc,
            )
        }
    }
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
    info!(log, "Multi-barcode record context: {} barcode levels", rec_ctx.num_barcodes());

    // Load known sample barcodes (with rotation → canonical mapping)
    let sample_info = load_sample_barcode_list(sample_bc_list_path, log)?;

    // Build sample barcode correction map (rotation → canonical)
    let (sample_permit_map, sample_bc_to_idx) = build_sample_permit_map(
        &sample_info,
        &gpl_opts.sample_correction_mode,
        log,
    )?;

    // Get sample names from the barcode file (uses canonical → name mapping)
    let sample_names = get_sample_names(&sample_info);

    // Load external cell barcode whitelist if provided (e.g. 737K 10x list).
    // This is used to filter cell barcodes — only whitelist BCs are counted directly;
    // non-matching BCs are collected for 1-edit correction later.
    let (cell_bc_whitelist, min_reads) = match &gpl_opts.fmeth {
        CellFilterMethod::UnfilteredExternalList(wl_path, mr) => {
            info!(log, "Loading external cell barcode whitelist from {}", wl_path.display());
            let wl_file = File::open(wl_path)
                .with_context(|| format!("couldn't open whitelist {}", wl_path.display()))?;
            let mut first_bclen = 0usize;
            let wl = populate_unfiltered_barcode_map(BufReader::new(wl_file), &mut first_bclen);
            info!(log, "Loaded {} cell barcodes from whitelist (bclen={})", wl.len(), first_bclen);
            (Some(wl), *mr as u64)
        }
        _ => (None, 0u64),
    };

    // First pass: read all records in parallel, correct sample BCs, count cell BCs per sample.
    // When a whitelist is present, only count cell BCs that are in the whitelist;
    // collect non-matching BCs per-sample for later 1-edit correction.
    info!(log, "First pass: counting cell barcodes per sample ({} threads)...", gpl_opts.threads);
    let num_samples = sample_info.canonical_barcodes.len();
    // Per-sample cell barcode frequency (only whitelist-matching BCs) — thread-safe
    let per_sample_cell_hist: Vec<DashMap<u64, u64, ahash::RandomState>> =
        (0..num_samples).map(|_| DashMap::default()).collect();
    // Per-sample unmatched cell barcodes (for 1-edit rescue) — per-worker, merged after
    let per_sample_unmatched: Vec<std::sync::Mutex<Vec<u64>>> =
        (0..num_samples).map(|_| std::sync::Mutex::new(Vec::new())).collect();

    let total_reads = std::sync::atomic::AtomicU64::new(0);
    let matched_reads = std::sync::atomic::AtomicU64::new(0);
    let unmatched_reads = std::sync::atomic::AtomicU64::new(0);

    let nworkers = gpl_opts.threads.max(1);
    let mut chunk_reader = libradicl::readers::ParallelChunkReader::<MultiBarcodeReadRecord>::new(
        &prelude,
        std::num::NonZeroUsize::new(nworkers).unwrap(),
    );

    let num_chunks = prelude.hdr.num_chunks;
    let pbar = ProgressBar::with_draw_target(
        Some(num_chunks),
        ProgressDrawTarget::stderr_with_hz(5),
    );
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
    let sample_permit_map = std::sync::Arc::new(sample_permit_map);
    let sample_bc_to_idx = std::sync::Arc::new(sample_bc_to_idx);
    let cell_bc_whitelist = std::sync::Arc::new(cell_bc_whitelist);
    let total_reads_arc = &total_reads;
    let matched_reads_arc = &matched_reads;
    let unmatched_reads_arc = &unmatched_reads;

    std::thread::scope(|s| {
        // Spawn worker threads
        let mut handles = Vec::new();
        for _ in 0..nworkers {
            let q = chunk_reader.get_queue();
            let done = chunk_reader.is_done();
            let hist = per_sample_cell_hist.clone();
            let unmatched = per_sample_unmatched.clone();
            let spm = sample_permit_map.clone();
            let s2i = sample_bc_to_idx.clone();
            let wl = cell_bc_whitelist.clone();

            let handle = s.spawn(move || {
                let mut local_total = 0u64;
                let mut local_matched = 0u64;
                let mut local_unmatched_sample = 0u64;
                // Per-sample local unmatched buffers to reduce lock contention
                let num_s = hist.len();
                let mut local_unmatched_bufs: Vec<Vec<u64>> =
                    (0..num_s).map(|_| Vec::new()).collect();

                while !done.load(Ordering::SeqCst) || !q.is_empty() {
                    while let Some(meta_chunk) = q.pop() {
                        for c in meta_chunk.iter() {
                            for read in &c.reads {
                                local_total += 1;

                                if !read.has_alignment_on_strand(expected_ori) {
                                    continue;
                                }

                                let sample_bc: u64 = read.collation_key_at_level(0).into();
                                let cell_bc: u64 = read.collate_key().into();

                                if let Some(&corrected_sample) = spm.get(&sample_bc) {
                                    if let Some(&sample_idx) = s2i.get(&corrected_sample) {
                                        local_matched += 1;
                                        if let Some(ref wl_map) = *wl {
                                            if wl_map.contains_key(&cell_bc) {
                                                *hist[sample_idx].entry(cell_bc).or_insert(0) += 1;
                                            } else {
                                                local_unmatched_bufs[sample_idx].push(cell_bc);
                                            }
                                        } else {
                                            *hist[sample_idx].entry(cell_bc).or_insert(0) += 1;
                                        }
                                    }
                                } else {
                                    local_unmatched_sample += 1;
                                }
                            }
                        }
                    }
                }

                // Flush local unmatched buffers
                for (idx, buf) in local_unmatched_bufs.into_iter().enumerate() {
                    if !buf.is_empty() {
                        unmatched[idx].lock().unwrap().extend(buf);
                    }
                }

                (local_total, local_matched, local_unmatched_sample)
            });
            handles.push(handle);
        }

        // Start the chunk reader (this feeds the queue from the BufReader)
        let cb = |_new_bytes: u64, new_rec: u64| {
            pbar.inc(new_rec);
        };
        let _ = chunk_reader.start(&mut ifile, Some(cb));

        // Join workers and aggregate counts
        for handle in handles {
            let (lt, lm, lu) = handle.join().expect("worker thread panicked");
            total_reads_arc.fetch_add(lt, Ordering::SeqCst);
            matched_reads_arc.fetch_add(lm, Ordering::SeqCst);
            unmatched_reads_arc.fetch_add(lu, Ordering::SeqCst);
        }
    });
    pbar.finish_and_clear();

    let total_reads = total_reads.load(Ordering::SeqCst);
    let matched_reads = matched_reads.load(Ordering::SeqCst);
    let unmatched_reads = unmatched_reads.load(Ordering::SeqCst);

    // Unwrap the Arcs (single owner again after threads join)
    let sample_permit_map = std::sync::Arc::into_inner(sample_permit_map).unwrap();
    let sample_bc_to_idx = std::sync::Arc::into_inner(sample_bc_to_idx).unwrap();
    let per_sample_cell_hist = std::sync::Arc::into_inner(per_sample_cell_hist).unwrap();
    let per_sample_unmatched = std::sync::Arc::into_inner(per_sample_unmatched).unwrap();

    info!(
        log,
        "First pass complete: {} total reads, {} matched to samples, {} unmatched",
        total_reads.to_formatted_string(&Locale::en),
        matched_reads.to_formatted_string(&Locale::en),
        unmatched_reads.to_formatted_string(&Locale::en),
    );

    // Create output directory
    let parent = std::path::Path::new(output_dir);
    std::fs::create_dir_all(parent).with_context(|| {
        format!("couldn't create output directory {}", parent.display())
    })?;

    // Write sample_permit_map.bin
    let sample_map_path = parent.join("sample_permit_map.bin");
    let sample_map_file = File::create(&sample_map_path)?;
    bincode::serialize_into(BufWriter::new(sample_map_file), &sample_permit_map)
        .map_err(|e| anyhow!("failed to serialize sample permit map: {}", e))?;
    info!(log, "Wrote sample permit map to {}", sample_map_path.display());

    // Per-sample: generate cell barcode permit lists
    let mut total_cells: u64 = 0;
    let mut sample_info_entries = Vec::new();

    // Get the cell barcode length from file tags
    let cell_bc_tag = format!("b{}len", num_barcodes - 1);
    let cell_bc_len: u16 = file_tag_map
        .get(&cell_bc_tag)
        .unwrap_or_else(|| panic!("expected '{}' file-level tag", cell_bc_tag))
        .try_into()
        .unwrap_or_else(|_| panic!("couldn't parse '{}' as u16", cell_bc_tag));

    for (sample_idx, sample_bc) in sample_info.canonical_barcodes.iter().enumerate() {
        let sample_name = &sample_names[sample_idx];
        let sample_dir = parent.join(format!("sample_{}", sample_name));
        std::fs::create_dir_all(&sample_dir)?;

        // Convert DashMap to HashMap for this sample
        let cell_hist: HashMap<u64, u64, ahash::RandomState> =
            per_sample_cell_hist[sample_idx].clone().into_iter().collect();
        let num_cells_observed = cell_hist.len();

        info!(
            log,
            "Sample '{}' (bc=0x{:x}): {} distinct cell barcodes observed",
            sample_name,
            sample_bc,
            num_cells_observed.to_formatted_string(&Locale::en),
        );

        if cell_hist.is_empty() {
            warn!(log, "Sample '{}' has no reads — skipping permit list generation", sample_name);
            sample_info_entries.push(serde_json::json!({
                "name": sample_name,
                "barcode": format!("0x{:x}", sample_bc),
                "num_reads": 0_u64,
                "num_cells": 0_u64,
            }));
            continue;
        }

        // Determine which cell BCs to keep, applying the same logic as single-barcode.
        // For UnfilteredExternalList: keep all whitelist BCs with >= min_reads,
        // then rescue unmatched BCs via 1-edit correction.
        let mut kept_bcs: Vec<u64> = Vec::new();
        let mut sample_freq_map: HashMap<u64, u64, ahash::RandomState> = HashMap::default();

        match &gpl_opts.fmeth {
            CellFilterMethod::UnfilteredExternalList(_, _) => {
                // Filter by min_reads: keep whitelist BCs that pass threshold
                let mut below_threshold_bcs = Vec::new();
                for (bc, count) in cell_hist.iter() {
                    if *count >= min_reads {
                        kept_bcs.push(*bc);
                        sample_freq_map.insert(*bc, *count);
                    } else {
                        // BCs below threshold: add their reads to unmatched for rescue
                        for _ in 0..*count {
                            below_threshold_bcs.push(*bc);
                        }
                    }
                }

                info!(
                    log,
                    "  {} whitelist BCs pass min_reads={} for sample '{}'",
                    kept_bcs.len(), min_reads, sample_name,
                );

                // Build BarcodeLookupMap from kept BCs for 1-edit correction
                let bcmap = BarcodeLookupMap::new(kept_bcs.clone(), cell_bc_len as u32);

                // Rescue unmatched cell BCs (not in whitelist) via 1-edit correction
                let mut unmatched = per_sample_unmatched[sample_idx].lock().unwrap();
                unmatched.append(&mut below_threshold_bcs);
                unmatched.sort_unstable();

                let mut found_approx = 0usize;
                let mut ambig_approx = 0usize;
                let mut not_found = 0usize;
                let mut corrected_list = Vec::<(u64, u64)>::new();

                for (count, ubc) in unmatched.iter().dedup_with_count() {
                    match bcmap.find_neighbors(*ubc, false) {
                        (Some(x), n) => {
                            let cbc = bcmap.barcodes[x];
                            if cbc != *ubc && n == 1 {
                                *sample_freq_map.entry(cbc).or_insert(0) += count as u64;
                                corrected_list.push((*ubc, cbc));
                                found_approx += count;
                            }
                            if n > 1 {
                                ambig_approx += count;
                            }
                        }
                        (None, _) => {
                            not_found += count;
                        }
                    }
                }

                info!(
                    log,
                    "  {} unmatched BCs for sample '{}': {} rescued (1-edit), {} ambiguous, {} not found",
                    unmatched.len(), sample_name,
                    found_approx.to_formatted_string(&Locale::en),
                    ambig_approx.to_formatted_string(&Locale::en),
                    not_found.to_formatted_string(&Locale::en),
                );

                // Write the corrected_list for this sample
                let corr_path = sample_dir.join("permit_map.bin");
                // Build the full permit map: identity for kept + correction for rescued
                let mut full_permit_list: HashMap<u64, u64> = HashMap::with_capacity(
                    kept_bcs.len() + corrected_list.len(),
                );
                for &bc in &kept_bcs {
                    full_permit_list.insert(bc, bc);
                }
                for (uncorrected, corrected) in &corrected_list {
                    full_permit_list.entry(*uncorrected).or_insert(*corrected);
                }
                let map_file = File::create(&corr_path)?;
                bincode::serialize_into(BufWriter::new(map_file), &full_permit_list)
                    .map_err(|e| anyhow!("failed to serialize permit map: {}", e))?;
            }
            CellFilterMethod::KneeFinding => {
                let mut freq: Vec<u64> = cell_hist.values().copied().collect();
                freq.sort_unstable();
                freq.reverse();
                let knee = knee_finding::get_knee(&freq, 100, log);
                let threshold = if knee > 0 { freq[knee.saturating_sub(1)] } else { 0 };
                for (bc, count) in cell_hist.iter() {
                    if *count >= threshold {
                        kept_bcs.push(*bc);
                        sample_freq_map.insert(*bc, *count);
                    }
                }
                info!(log, "  Knee at {}, {} cells retained for sample '{}'", knee, kept_bcs.len(), sample_name);

                let full_permit_list = afutils::generate_permitlist_map(&kept_bcs, cell_bc_len as usize)
                    .map_err(|e| anyhow!("failed to generate permit list map: {}", e))?;
                let map_path = sample_dir.join("permit_map.bin");
                let map_file = File::create(&map_path)?;
                bincode::serialize_into(BufWriter::new(map_file), &full_permit_list)
                    .map_err(|e| anyhow!("failed to serialize permit map: {}", e))?;
            }
            _ => {
                // ForceCells, ExpectCells, ExplicitList: use frequency-based selection
                let mut freq: Vec<u64> = cell_hist.values().copied().collect();
                freq.sort_unstable();
                freq.reverse();
                let num_cells = match &gpl_opts.fmeth {
                    CellFilterMethod::ForceCells(n) => (*n).min(freq.len()),
                    CellFilterMethod::ExpectCells(n) => {
                        let threshold = freq[0] / *n as u64;
                        let idx = freq.iter().position(|&f| f < threshold).unwrap_or(freq.len());
                        (idx as f64 * 10.0) as usize
                    }
                    _ => freq.len(),
                }.min(freq.len());
                let threshold = if num_cells > 0 { freq[num_cells.saturating_sub(1)] } else { 0 };
                for (bc, count) in cell_hist.iter() {
                    if *count >= threshold {
                        kept_bcs.push(*bc);
                        sample_freq_map.insert(*bc, *count);
                    }
                }
                info!(log, "  {} cells retained for sample '{}'", kept_bcs.len(), sample_name);

                let full_permit_list = afutils::generate_permitlist_map(&kept_bcs, cell_bc_len as usize)
                    .map_err(|e| anyhow!("failed to generate permit list map: {}", e))?;
                let map_path = sample_dir.join("permit_map.bin");
                let map_file = File::create(&map_path)?;
                bincode::serialize_into(BufWriter::new(map_file), &full_permit_list)
                    .map_err(|e| anyhow!("failed to serialize permit map: {}", e))?;
            }
        }

        // Write per-sample permit_freq.bin
        let freq_path = sample_dir.join("permit_freq.bin");
        afutils::write_permit_list_freq(&freq_path, cell_bc_len, &sample_freq_map)
            .map_err(|e| anyhow!("failed to write permit freq for sample '{}': {}", sample_name, e))?;

        total_cells += kept_bcs.len() as u64;

        let sample_total_reads: u64 = sample_freq_map.values().sum();
        sample_info_entries.push(serde_json::json!({
            "name": sample_name,
            "barcode": format!("0x{:x}", sample_bc),
            "num_reads": sample_total_reads,
            "num_cells": kept_bcs.len(),
        }));
    }

    // Write sample_info.json
    let sample_info = serde_json::json!({
        "num_samples": num_samples,
        "num_barcodes": num_barcodes,
        "total_cells": total_cells,
        "total_reads": total_reads,
        "matched_reads": matched_reads,
        "unmatched_reads": unmatched_reads,
        "sample_correction_mode": format!("{:?}", gpl_opts.sample_correction_mode),
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
    });
    let gpl_meta_path = parent.join("generate_permit_list.json");
    let gpl_meta_file = File::create(&gpl_meta_path)?;
    serde_json::to_writer_pretty(BufWriter::new(gpl_meta_file), &gpl_meta)?;

    info!(
        log,
        "Multi-barcode permit list generation complete: {} samples, {} total cells",
        num_samples, total_cells,
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
fn load_sample_barcode_list(path: &PathBuf, log: &slog::Logger) -> anyhow::Result<SampleBarcodeInfo> {
    let file = File::open(path)
        .with_context(|| format!("couldn't open sample barcode list: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut rotation_to_canonical: HashMap<u64, u64> = HashMap::new();
    let mut canonical_to_name: HashMap<u64, String> = HashMap::new();
    let mut canonical_order: Vec<u64> = Vec::new(); // preserves order, deduplicated
    let mut has_tsv_columns = false;

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

        let pack = |seq: &str| -> anyhow::Result<u64> {
            needletail::bitkmer::BitNuclKmer::new(seq.as_bytes(), seq.len() as u8, false)
                .next()
                .map(|(_, kmer, _)| kmer.0)
                .ok_or_else(|| anyhow!("couldn't pack sample barcode: {}", seq))
        };

        let obs_packed = pack(observed_seq)?;
        let canon_packed = pack(canonical_seq)?;

        rotation_to_canonical.insert(obs_packed, canon_packed);
        if !canonical_to_name.contains_key(&canon_packed) {
            canonical_order.push(canon_packed);
        }
        canonical_to_name.insert(canon_packed, name);
    }

    let num_rotations = rotation_to_canonical.len();
    let num_canonical = canonical_order.len();

    if has_tsv_columns && num_rotations > num_canonical {
        info!(
            log,
            "Loaded {} rotation barcodes mapping to {} canonical samples from {}",
            num_rotations, num_canonical, path.display(),
        );
    } else {
        info!(
            log,
            "Loaded {} sample barcodes from {}",
            num_canonical, path.display(),
        );
    }

    Ok(SampleBarcodeInfo {
        canonical_barcodes: canonical_order,
        rotation_to_canonical,
        canonical_to_name,
    })
}

/// Build the sample barcode permit map (observed -> canonical).
/// Returns (permit_map, bc_to_idx) where:
///   - permit_map maps every observed rotation barcode to its canonical barcode
///   - bc_to_idx maps canonical barcode to sample index
fn build_sample_permit_map(
    sample_info: &SampleBarcodeInfo,
    correction_mode: &SampleCorrectionMode,
    log: &slog::Logger,
) -> anyhow::Result<(HashMap<u64, u64>, HashMap<u64, usize>)> {
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

    match correction_mode {
        SampleCorrectionMode::Exact => {
            info!(
                log,
                "Sample barcode correction: exact match ({} rotation entries -> {} canonical samples)",
                permit_map.len(),
                bc_to_idx.len(),
            );
        }
        SampleCorrectionMode::OneEdit => {
            info!(log, "Sample barcode correction mode: 1-edit distance");
            // Generate 1-edit neighbors for ALL observed rotation barcodes
            let all_observed: Vec<u64> = sample_info.rotation_to_canonical.keys().copied().collect();
            let bc_len = if let Some(&first) = all_observed.first() {
                // Estimate barcode length from packed value
                // For 8bp barcodes: 2*8 = 16 bits used
                8usize // TODO: derive from actual barcode length
            } else {
                8usize
            };
            let full_map = afutils::generate_permitlist_map(&all_observed, bc_len)
                .map_err(|e| anyhow!("failed to generate sample permit map: {}", e))?;
            for (raw, corrected) in &full_map {
                if !permit_map.contains_key(raw) {
                    // Map the 1-edit neighbor to the same canonical as its corrected barcode
                    if let Some(&canon) = sample_info.rotation_to_canonical.get(corrected) {
                        permit_map.insert(*raw, canon);
                    }
                }
            }
            info!(
                log,
                "Built sample permit map with {} entries ({} exact + {} corrected) -> {} canonical samples",
                permit_map.len(),
                sample_info.rotation_to_canonical.len(),
                permit_map.len() - sample_info.rotation_to_canonical.len(),
                bc_to_idx.len(),
            );
        }
    }

    Ok((permit_map, bc_to_idx))
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
    unmatched_bc: &mut Vec<u64>,
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
                    unmatched_bc.push(r.collate_key().into());
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
                        unmatched_bc.push(r.collate_key().into());
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
                            unmatched_bc.push(r.collate_key().into());
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
                            unmatched_bc.push(r.collate_key().into());
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
    let valid_bc: Vec<u64> = hist
        .iter()
        .filter_map(|(k, v)| if v >= &min_freq { Some(*k) } else { None })
        .collect();
    valid_bc
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

    check_permit_list_validity(gpl_opts.log, unmatched_bc.len(), num_reads)?;

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
) -> (usize, usize, Vec<u64>, usize)
where
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize,
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    <R as MappedRecord>::ParsingContext: RecordContext + Clone + Send,
{
    std::thread::scope(|s| {
        let mut handles = Vec::new();

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
        let mut all_unmatched = Vec::new();
        let mut max_ambig = 0;

        for handle in handles {
            let (nr, nocr, ubc, mar) = handle.join().expect("The parsing thread panicked");
            total_reads += nr;
            total_compat += nocr;
            all_unmatched.extend_from_slice(&ubc);
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
            let rd = rad_reader.is_done();
            let q = rad_reader.get_queue();
            let hm = hm.clone();

            let handle = s.spawn(move || {
                let mut max_ambiguity_read = 0usize;
                let mut num_reads = 0;
                let mut num_orientation_compat_reads = 0;

                while !rd.load(Ordering::SeqCst) || !q.is_empty() {
                    while let Some(meta_chunk) = q.pop() {
                        for c in meta_chunk.iter() {
                            num_orientation_compat_reads +=
                                update_barcode_hist(&hm, &mut max_ambiguity_read, &c, expected_ori);
                            num_reads += c.reads.len();
                        }
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
