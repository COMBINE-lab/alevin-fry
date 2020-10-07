// Copyright 2020 Rob Patro, Avi Srivastava. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

extern crate fasthash;
extern crate slog;
use self::slog::crit;
use self::slog::info;

use crate as libradicl;
use bio_types::strand::Strand;
use fasthash::sea::Hash64;
use fasthash::RandomState;
use libradicl::exit_codes;
use needletail::bitkmer::*;
use num_format::{Locale, ToFormattedString};
use serde_json::json;
use std::cmp;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::{BufWriter, Write};
use std::str::from_utf8;

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
    ExplicitList(String),
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
    assert!(sorted_frequencies.len() >= 2);
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

/// Given the input RAD file `input_file`, compute
/// and output (in `output_dir`) the list of valid
/// (i.e. "permitted") barcode values, as well as
/// a map from each correctable barcode to the
/// permitted barcode to which it maps.
pub fn generate_permit_list(
    rad_dir: String,
    output_dir: String,
    filter_meth: CellFilterMethod,
    expected_ori: Strand,
    //top_k: Option<usize>,
    //valid_bc_file: Option<String>,
    //use_knee_distance: bool,
    log: &slog::Logger,
) -> Result<u64, Box<dyn std::error::Error>> {
    let i_dir = std::path::Path::new(&rad_dir);

    if !i_dir.exists() {
        crit!(log, "the input RAD path {} does not exist", rad_dir);
        std::process::exit(1);
    }

    let i_file = File::open(i_dir.join("map.rad")).expect("could not open input rad file");
    let mut br = BufReader::new(i_file);
    let hdr = libradicl::RADHeader::from_bytes(&mut br);
    info!(
        log,
        "paired : {:?}, ref_count : {}, num_chunks : {}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en)
    );
    // file-level
    let fl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());

    // right now, we only handle BC and UMI types of U8—U64, so validate that
    const BNAME: &str = "b";
    const UNAME: &str = "u";

    let mut bct: Option<u8> = None;
    let mut umit: Option<u8> = None;

    for rt in &rl_tags.tags {
        // if this is one of our tags
        if rt.name == BNAME || rt.name == UNAME {
            if libradicl::decode_int_type_tag(rt.typeid).is_none() {
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

    // alignment-level
    let al_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    let ft_vals = libradicl::FileTags::from_bytes(&mut br);
    info!(log, "File-level tag values {:?}", ft_vals);

    let mut num_reads: usize = 0;

    let s = RandomState::<Hash64>::new();
    let mut hm = HashMap::with_hasher(s);
    let bc_type = libradicl::decode_int_type_tag(bct.expect("no barcode tag description present."))
        .expect("unknown barcode type id.");
    let umi_type = libradicl::decode_int_type_tag(umit.expect("no umi tag description present"))
        .expect("unknown barcode type id.");

    for _ in 0..(hdr.num_chunks as usize) {
        let c = libradicl::Chunk::from_bytes(&mut br, &bc_type, &umi_type);
        libradicl::update_barcode_hist(&mut hm, &c, &expected_ori);
        num_reads += c.reads.len();
    }

    info!(
        log,
        "observed {} reads in {} chunks",
        num_reads.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en)
    );

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
}
