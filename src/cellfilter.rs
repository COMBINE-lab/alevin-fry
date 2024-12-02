/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use anyhow::Context;
use dashmap::DashMap;
use slog::crit;
use slog::{info, warn};

use crate::diagnostics;
use crate::prog_opts::GenPermitListOpts;
use crate::utils as afutils;
#[allow(unused_imports)]
use ahash::{AHasher, RandomState};
use bio_types::strand::Strand;
use bstr::io::BufReadExt;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::Itertools;
use libradicl::exit_codes;
use libradicl::rad_types::{self, RadType};
use libradicl::BarcodeLookupMap;
use libradicl::{chunk, record::AlevinFryReadRecord};
use needletail::bitkmer::*;
use num_format::{Locale, ToFormattedString};
use serde::Serialize;
use serde_json::json;
use std::cmp;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::{atomic::Ordering, Arc};
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

struct Point {
    x: f64,
    y: f64,
}

/// compute the distance between the query point `Q`
/// and the line defined by points `P1` and `P2`.  The
/// formula used here is taken from :
/// https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
fn distance_to_line(p1: &Point, p2: &Point, q: &Point) -> f64 {
    let x_0 = q.x;
    let y_0 = q.y;

    let x_1 = p1.x;
    let y_1 = p1.y;

    let x_2 = p2.x;
    let y_2 = p2.y;

    let numer = ((y_2 - y_1) * x_0 - (x_2 - x_1) * y_0 + x_2 * y_1 - y_2 * x_1).abs();
    let denom = ((y_2 - y_1).powi(2) + (x_2 - x_1).powi(2)).sqrt();
    assert!(denom > 0.0f64);
    numer / denom
}

/// This method is a implementation of the distance method
/// used in umi_tools :
///     Smith, Tom, Andreas Heger, and Ian Sudbery.
///     "UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy."
///     Genome research 27.3 (2017): 491-499.
///
/// though this is a re-implementation and uses the same basic algorithm
/// the result may not be identical.
///
/// Given a list of cumulative frequencies, where the index of a
/// point is interpreted as its x-coordinate and its frequency
/// is interpreted as its y-coordinate, define the line L based on
/// `sorted_frequences.first()` and `sorted_frequencies.last()`.  Compute
/// the distance of each point from L, and return the index of the point
/// having the maximum distance.
fn get_max_distance_index(sorted_frequencies: &[u64], is_cumulative: bool) -> usize {
    assert!(sorted_frequencies.len() >= 2,
        "ERROR: when attempting to find a knee-distance threshold, the list of putative cells is only of length {}. Cannot proceed. Please check the mapping rate.",
        sorted_frequencies.len());
    let first = sorted_frequencies
        .first()
        .expect("cannot process empty frequency list.");
    let last = sorted_frequencies
        .last()
        .expect("cannot process empty frequency list.");

    // length as a float
    let max_x = sorted_frequencies.len() as f64;

    // if the distribution is cumulative, then the smallest y coordinate is
    // f, otherewise it is l
    let max_y = if is_cumulative {
        *last as f64
    } else {
        *first as f64
    };

    let p1 = Point {
        x: 0.0f64,
        y: (*first as f64) / max_y,
    };
    let p2 = Point {
        x: 1.0f64,
        y: (*last as f64) / max_y,
    };

    let mut max_d: f64 = -1.0;
    let mut max_ind: usize = 0;

    for (ind, freq) in sorted_frequencies.iter().enumerate() {
        let x = ind as f64 / max_x;
        let y = *freq as f64 / max_y;
        let q = Point { x, y };
        let d = distance_to_line(&p1, &p2, &q);
        if d >= max_d {
            max_d = d;
            max_ind = ind;
        }
    }
    max_ind
}

/// Get the knee of the cure using the `distance` method as described
/// in the [UMI-tools documentation](https://github.com/CGATOxford/UMI-tools).
/// This method takes a reverse-sorted (sorted in descending order) llist of
/// frequencies, and a maximum number of iterations to run the algorithm.  It
/// returns the point on the CDF of the reverse-sorted frequency vector that is
/// farthest from the line defined by the end-points.  The algorithm is taken from
/// [here](https://github.com/CGATOxford/UMI-tools/blob/master/umi_tools/whitelist_methods.py#L248).
fn get_knee(freq: &[u64], max_iterations: usize, log: &slog::Logger) -> usize {
    // get the cumulative frequency from the frequency
    let cfreq: Vec<u64> = freq
        .iter()
        .scan(0u64, |acc, &num| {
            *acc += num;
            Some(*acc)
        })
        .collect();
    // get the guess about the max distance point
    let mut prev_max = 0;
    let mut max_idx = get_max_distance_index(&cfreq[..], true);

    // if we think we should include no cells, something is probably wrong.
    assert_ne!(max_idx, 0,
	       "get_knee determined a knee index of 0. This probably should not happen with valid input data.");

    let mut iterations = 0;
    let iter_slack = 5;
    // while our algorithm hasn't converged
    while max_idx - prev_max != 0 {
        info!(log, "max_idx = {}", max_idx);
        prev_max = max_idx;
        iterations += 1;
        if iterations % 10 == 0 {
            info!(log, "knee-finding iter = {}", iterations);
        }
        if iterations > max_iterations {
            break;
        }
        let last_idx = std::cmp::min(cfreq.len() - 1, max_idx * iter_slack);
        max_idx = get_max_distance_index(&cfreq[0..last_idx], true);
        assert_ne!(max_idx, 0,
	       "get_knee determined a knee index of 0. This probably should not happen with valid input data.");
    }
    max_idx
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
            let num_bc = get_knee(&freq[..], 100, log);
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
    let rad_dir = gpl_opts.input_dir;
    let output_dir = gpl_opts.output_dir;
    let filter_meth = gpl_opts.fmeth.clone();
    let expected_ori = gpl_opts.expected_ori;
    let version = gpl_opts.version;
    let velo_mode = gpl_opts.velo_mode;
    let cmdline = gpl_opts.cmdline;
    let log = gpl_opts.log;

    let i_dir = std::path::Path::new(&rad_dir);

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

    let mut first_bclen = 0usize;
    let mut unfiltered_bc_counts = None;
    if let CellFilterMethod::UnfilteredExternalList(fname, _) = &filter_meth {
        let i_file = File::open(fname).context("could not open input file")?;
        let br = BufReader::new(i_file);
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
    let i_file = File::open(i_dir.join("map.rad")).context("could not open input rad file")?;
    let ifile = BufReader::new(i_file);
    let mut rad_reader = libradicl::readers::ParallelRadReader::<
        AlevinFryReadRecord,
        BufReader<File>,
    >::new(ifile, NonZeroUsize::new(nworkers).unwrap());

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

    let nc = num_chunks.expect("unknwon number of chunks").get() as u64;
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
            unmatched_bc = Vec::with_capacity(10000000);
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
                            while !rd.load(Ordering::SeqCst) {
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
                    pbar.finish_with_message(format!("finished parsing RAD file\n",));
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
                        while !rd.load(Ordering::SeqCst) {
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
                // return the hash map we no longer need
                Arc::<DashMap<u64, u64, ahash::RandomState>>::into_inner(hm)
                    .expect("unique reference to DashMap")
            });
            pbar.finish_with_message(format!("finished parsing RAD file\n",));
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

    /*
    let valid_bc: Vec<u64>;
    let mut freq: Vec<u64> = hm.values().cloned().collect();
    freq.sort_unstable();
    freq.reverse();

    // select from among supported filter methods
    match filter_meth {
    CellFilterMethod::KneeFinding => {
        let num_bc = get_knee(&freq[..], 100, &log);
        let min_freq = freq[num_bc];

        // collect all of the barcodes that have a frequency
        // >= to min_thresh.
        valid_bc = libradicl::permit_list_from_threshold(&hm, min_freq);
        info!(
        log,
        "knee distance method resulted in the selection of {} permitted barcodes.",
        valid_bc.len()
        );
    }
    CellFilterMethod::ForceCells(top_k) => {
        let num_bc = if freq.len() < top_k {
        freq.len() - 1
        } else {
        top_k - 1
        };

        let min_freq = freq[num_bc];

        // collect all of the barcodes that have a frequency
        // >= to min_thresh.
        valid_bc = libradicl::permit_list_from_threshold(&hm, min_freq);
    }
    CellFilterMethod::ExplicitList(valid_bc_file) => {
        valid_bc = libradicl::permit_list_from_file(valid_bc_file, ft_vals.bclen);
    }
    CellFilterMethod::ExpectCells(expected_num_cells) => {
        let robust_quantile = 0.99f64;
        let robust_div = 10.0f64;
        let robust_ind = (expected_num_cells as f64 * robust_quantile).round() as u64;
        // the robust ind must be valid
        let ind = cmp::min(freq.len() - 1, robust_ind as usize);
        let robust_freq = freq[ind];
        let min_freq = std::cmp::max(1u64, (robust_freq as f64 / robust_div).round() as u64);
        valid_bc = libradicl::permit_list_from_threshold(&hm, min_freq);
    }
    CellFilterMethod::UnfilteredExternalList(_, min_reads) => {
        unimplemented!();
    }
    }

    // generate the map from each permitted barcode to all barcodes within
    // edit distance 1 of it.
    let full_permit_list =
    libradicl::utils::generate_permitlist_map(&valid_bc, ft_vals.bclen as usize).unwrap();

    let s2 = RandomState::<Hash64>::new();
    let mut permitted_map = HashMap::with_capacity_and_hasher(valid_bc.len(), s2);

    let mut num_corrected = 0;
    for (k, v) in hm.iter() {
    if let Some(&valid_key) = full_permit_list.get(k) {
        *permitted_map.entry(valid_key).or_insert(0u64) += *v;
        num_corrected += 1;
        //println!("{} was a neighbor of {}, with count {}", k, valid_key, v);
    }
    }

    let parent = std::path::Path::new(&output_dir);
    std::fs::create_dir_all(&parent).unwrap();
    let o_path = parent.join("permit_freq.tsv");
    let output = std::fs::File::create(&o_path).expect("could not create output.");
    let mut writer = BufWriter::new(&output);

    for (k, v) in permitted_map {
    //let bc_mer: BitKmer = (k, ft_vals.bclen as u8);
    writeln!(
        &mut writer,
        "{}\t{}",
        k,
        //from_utf8(&bitmer_to_bytes(bc_mer)[..]).unwrap(),
        v
    )
    .expect("couldn't write to output file.");
    }

    let o_path = parent.join("all_freq.tsv");
    let output = std::fs::File::create(&o_path).expect("could not create output.");
    let mut writer = BufWriter::new(&output);
    for (k, v) in hm {
    let bc_mer: BitKmer = (k, ft_vals.bclen as u8);
    writeln!(
        &mut writer,
        "{}\t{}",
        from_utf8(&bitmer_to_bytes(bc_mer)[..]).unwrap(),
        v
    )
    .expect("couldn't write to output file.");
    }

    let s_path = parent.join("permit_map.bin");
    let s_file = std::fs::File::create(&s_path).expect("could not create serialization file.");
    let mut s_writer = BufWriter::new(&s_file);
    bincode::serialize_into(&mut s_writer, &full_permit_list)
    .expect("couldn't serialize permit list.");

    let meta_info = json!({
    "expected_ori" : expected_ori.strand_symbol()
    });

    let m_path = parent.join("generate_permit_list.json");
    let mut m_file = std::fs::File::create(&m_path).expect("could not create metadata file.");

    let meta_info_string =
    serde_json::to_string_pretty(&meta_info).expect("could not format json.");
    m_file
    .write_all(meta_info_string.as_bytes())
    .expect("cannot write to generate_permit_list.json file");

    info!(
    log,
    "total number of corrected barcodes : {}",
    num_corrected.to_formatted_string(&Locale::en)
    );
    Ok(num_corrected)
    */
}

pub fn update_barcode_hist_unfiltered(
    hist: &DashMap<u64, u64, ahash::RandomState>,
    unmatched_bc: &mut Vec<u64>,
    max_ambiguity_read: &mut usize,
    chunk: &chunk::Chunk<AlevinFryReadRecord>,
    expected_ori: &Strand,
) -> usize {
    let mut num_strand_compat_reads = 0usize;
    match expected_ori {
        Strand::Unknown => {
            for r in &chunk.reads {
                num_strand_compat_reads += 1;
                *max_ambiguity_read = r.refs.len().max(*max_ambiguity_read);
                // lookup the barcode in the map of unfiltered known
                // barcodes
                match hist.get_mut(&r.bc) {
                    // if we find a match, increment the count
                    Some(mut c) => *c += 1,
                    // otherwise, push this into the unmatched list
                    None => {
                        unmatched_bc.push(r.bc);
                    }
                }
            }
        }
        Strand::Forward => {
            for r in &chunk.reads {
                if r.dirs.iter().any(|&x| x) {
                    num_strand_compat_reads += 1;
                    *max_ambiguity_read = r.refs.len().max(*max_ambiguity_read);
                    // lookup the barcode in the map of unfiltered known
                    // barcodes
                    match hist.get_mut(&r.bc) {
                        // if we find a match, increment the count
                        Some(mut c) => *c += 1,
                        // otherwise, push this into the unmatched list
                        None => {
                            unmatched_bc.push(r.bc);
                        }
                    }
                }
            }
        }
        Strand::Reverse => {
            for r in &chunk.reads {
                if r.dirs.iter().any(|&x| !x) {
                    num_strand_compat_reads += 1;
                    *max_ambiguity_read = r.refs.len().max(*max_ambiguity_read);
                    // lookup the barcode in the map of unfiltered known
                    // barcodes
                    match hist.get_mut(&r.bc) {
                        // if we find a match, increment the count
                        Some(mut c) => *c += 1,
                        // otherwise, push this into the unmatched list
                        None => {
                            unmatched_bc.push(r.bc);
                        }
                    }
                }
            }
        }
    }
    num_strand_compat_reads
}

pub fn update_barcode_hist(
    hist: &DashMap<u64, u64, ahash::RandomState>,
    max_ambiguity_read: &mut usize,
    chunk: &chunk::Chunk<AlevinFryReadRecord>,
    expected_ori: &Strand,
) {
    match expected_ori {
        Strand::Unknown => {
            for r in &chunk.reads {
                *max_ambiguity_read = r.refs.len().max(*max_ambiguity_read);
                *hist.entry(r.bc).or_insert(0) += 1;
            }
        }
        Strand::Forward => {
            for r in &chunk.reads {
                if r.dirs.iter().any(|&x| x) {
                    *max_ambiguity_read = r.refs.len().max(*max_ambiguity_read);
                    *hist.entry(r.bc).or_insert(0) += 1;
                }
            }
        }
        Strand::Reverse => {
            for r in &chunk.reads {
                if r.dirs.iter().any(|&x| !x) {
                    *max_ambiguity_read = r.refs.len().max(*max_ambiguity_read);
                    *hist.entry(r.bc).or_insert(0) += 1;
                }
            }
        }
    }
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
