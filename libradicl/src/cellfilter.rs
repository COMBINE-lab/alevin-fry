// Copyright 2020 Rob Patro, Avi Srivastava. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

extern crate fasthash;
extern crate slog;
use self::slog::info;

use crate as libradicl;
use fasthash::sea::Hash64;
use fasthash::RandomState;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::{BufWriter, Write};

pub fn generate_permit_list(
    input_file: String,
    output_dir: String,
    top_k: Option<usize>,
    valid_bc_file: Option<String>,
    log: &slog::Logger,
) -> Result<u64, Box<dyn std::error::Error>> {
    let i_file = File::open(input_file).unwrap();
    let mut br = BufReader::new(i_file);
    let hdr = libradicl::RADHeader::from_bytes(&mut br);
    info!(
        log,
        "paired : {:?}, ref_count : {:?}, num_chunks : {:?}",
        hdr.is_paired,
        hdr.ref_count,
        hdr.num_chunks
    );
    // file-level
    let fl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    let ft_vals = libradicl::FileTags::from_bytes(&mut br);
    info!(log, "File-level tag values {:?}", ft_vals);

    let bct = rl_tags.tags[0].typeid;
    let umit = rl_tags.tags[1].typeid;

    let mut num_reads: usize = 0;

    let s = RandomState::<Hash64>::new();
    let mut hm = HashMap::with_hasher(s);
    let bc_type = libradicl::decode_int_type_tag(bct).expect("unknown barcode type id.");
    let umi_type = libradicl::decode_int_type_tag(umit).expect("unknown barcode type id.");

    for _ in 0..(hdr.num_chunks as usize) {
        let c = libradicl::Chunk::from_bytes(&mut br, &bc_type, &umi_type);
        libradicl::update_barcode_hist(&mut hm, &c);
        num_reads += c.reads.len();
    }

    info!(
        log,
        "observed {:?} reads in {:?} chunks", num_reads, hdr.num_chunks
    );

    let valid_bc: Vec<u64>;

    if top_k.is_some() {
        let top_k = top_k.unwrap() as usize;
        let mut freq: Vec<u64> = hm.values().cloned().collect();
        freq.sort_unstable();
        freq.reverse();

        let num_bc = if freq.len() < top_k {
            freq.len() - 1
        } else {
            top_k - 1
        };

        let min_freq = freq[num_bc];

        // collect all of the barcodes that have a frequency
        // >= to min_thresh.
        valid_bc = libradicl::permit_list_from_threshold(&hm, min_freq);
    } else {
        let valid_bc_file = valid_bc_file.expect("couldn't extract --valid-bc option.");
        valid_bc = libradicl::permit_list_from_file(valid_bc_file, ft_vals.bclen);
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
        writeln!(&mut writer, "{:?}\t{:?}", k, v).expect("couldn't write to output file.");
    }

    let s_path = parent.join("permit_map.bin");
    let s_file = std::fs::File::create(&s_path).expect("could not create serialization file.");
    let mut s_writer = BufWriter::new(&s_file);
    bincode::serialize_into(&mut s_writer, &full_permit_list)
        .expect("couldn't serialize permit list.");

    Ok(num_corrected)
}
