/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */
use crate::constants as afconst;
use crate::eq_class::IndexedEqList;
use anyhow::{anyhow, Context};
use bstr::io::BufReadExt;
use core::fmt;
use libradicl::utils::SPLICE_MASK_U32;
use needletail::bitkmer::*;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::str::FromStr;
use thiserror::Error;

/*
struct QuantArguments {
    num_threads: u64,
    num_bootstraps: u64,
    init_uniform: bool,
    summary_stat: bool,
    dump_eq: bool,
    use_mtx: bool,
    input_dir: String,
    output_dir: String,
    tg_map: String,
    resolution: ResolutionStrategy,
    sa_model: SplicedAmbiguityModel,
    small_thresh: u64,
    filter_list: String
}
*/

/// FROM https://github.com/10XGenomics/rust-debruijn/blob/master/src/dna_string.rs
/// count Hamming distance between 2 2-bit DNA packed u64s
pub(super) fn count_diff_2_bit_packed(a: u64, b: u64) -> usize {
    let bit_diffs = a ^ b;
    let two_bit_diffs = (bit_diffs | bit_diffs >> 1) & 0x5555555555555555;
    two_bit_diffs.count_ones() as usize
}

#[inline(always)]
fn unspliced_of(gid: u32) -> u32 {
    gid + 1
}

/// should always compile to no-op
#[inline(always)]
fn spliced_of(gid: u32) -> u32 {
    gid
}

// given a spliced or unspliced gene id, return
// the spliced (canonical) id for this gene.
#[inline(always)]
fn spliced_id(gid: u32) -> u32 {
    gid & SPLICE_MASK_U32
}

#[inline(always)]
pub fn same_gene(g1: u32, g2: u32, with_unspliced: bool) -> bool {
    (g1 == g2) || (with_unspliced && (spliced_id(g1) == spliced_id(g2)))
}

#[inline(always)]
pub fn is_spliced(gid: u32) -> bool {
    // if the id is even, then it's spliced
    (0x1 & gid) == 0
}

#[inline(always)]
pub fn is_unspliced(gid: u32) -> bool {
    // if it's not spliced, then it is unspliced
    !is_spliced(gid)
}

/// Write the permit_freq.bin and all_freq.bin files
pub fn write_permit_list_freq(
    o_path: &std::path::Path,
    bclen: u16,
    permit_freq_map: &HashMap<u64, u64, ahash::RandomState>,
) -> Result<(), Box<dyn std::error::Error>> {
    let output = std::fs::File::create(o_path)?;
    let mut writer = BufWriter::new(&output);

    {
        // the first u64 represents file format version.
        writer
            .write_all(&afconst::PERMIT_FILE_VER.to_le_bytes())
            .unwrap();

        // the second u64 represents barcode length
        writer.write_all(&(u64::from(bclen)).to_le_bytes()).unwrap();

        // the rest records the permitted barcode:freq hashmap
        bincode::serialize_into(&mut writer, &permit_freq_map)?;
    }
    Ok(())
}

/// Parse a 3 column tsv of the format
/// transcript_name gene_name   status
/// where status is one of S or U each gene will be allocated both a spliced and
/// unspliced variant, the spliced index will always be even and the unspliced odd,
/// and they will always be adjacent ids.  For example, if gene A is present in
/// the sample and it's spliced variant is assigned id i,  then it will always be true that
/// i % 2 == 0
/// and
/// (i+1) will be the id for the unspliced version of gene A
fn parse_tg_spliced_unspliced(
    rdr: &mut csv::Reader<File>,
    ref_count: usize,
    rname_to_id: &HashMap<String, u32, ahash::RandomState>,
    gene_names: &mut Vec<String>,
    gene_name_to_id: &mut HashMap<String, u32, ahash::RandomState>,
) -> anyhow::Result<(Vec<u32>, bool)> {
    // map each transcript id to the corresponding gene id
    // the transcript name can be looked up from the id in the RAD header,
    // and the gene name can be looked up from the id in the gene_names
    // vector.
    let mut tid_to_gid = vec![u32::MAX; ref_count];

    // Record will be transcript, gene, splicing status
    type TsvRec = (String, String, String);

    // the transcripts for which we've found a gene mapping
    let mut found = 0usize;

    // starting from 0, we assign each gene 2 ids (2 consecutive integers),
    // the even ids are for spliced txps, the odd ids are for unspliced txps
    // for convenience, we define a gid helper, next_gid
    let mut next_gid = 0u32;
    // apparently the "header" (first row) will be included
    // in the iterator returned by `deserialize` anyway
    /*let hdr = rdr.headers()?;
    let hdr_vec : Vec<Result<TsvRec,csv::Error>> = vec![hdr.deserialize(None)];
    */
    for result in rdr.deserialize() {
        let record: TsvRec = result?;
        // first, get the first id for this gene
        let gene_id = *gene_name_to_id.entry(record.1.clone()).or_insert_with(|| {
            // as we need to return the current next_gid if we run this code
            // we add by two and then return current gene id.
            let cur_gid = next_gid;
            next_gid += 2;
            // we haven't added this gene name already,
            // we append it now to the list of gene names.
            gene_names.push(record.1.clone());
            cur_gid
        });

        // get the transcript id
        if let Some(transcript_id) = rname_to_id.get(&record.0) {
            found += 1;
            if record.2.eq_ignore_ascii_case("U") {
                // This is an unspliced txp
                // we link it to the second gid of this gene
                tid_to_gid[*transcript_id as usize] = unspliced_of(gene_id);
            } else if record.2.eq_ignore_ascii_case("S") {
                // This is a spliced txp, we link it to the
                // first gid of this gene
                tid_to_gid[*transcript_id as usize] = spliced_of(gene_id);
            } else {
                return Err(anyhow!(
                    "Third column in 3 column txp-to-gene file must be S or U"
                ));
            }
        }
    }

    assert_eq!(
        found, ref_count,
        "The tg-map must contain a gene mapping for all transcripts in the header"
    );

    Ok((tid_to_gid, true))
}

fn parse_tg_spliced(
    rdr: &mut csv::Reader<File>,
    ref_count: usize,
    rname_to_id: &HashMap<String, u32, ahash::RandomState>,
    gene_names: &mut Vec<String>,
    gene_name_to_id: &mut HashMap<String, u32, ahash::RandomState>,
) -> anyhow::Result<(Vec<u32>, bool)> {
    // map each transcript id to the corresponding gene id
    // the transcript name can be looked up from the id in the RAD header,
    // and the gene name can be looked up from the id in the gene_names
    // vector.
    let mut tid_to_gid = vec![u32::MAX; ref_count];
    // now read in the transcript to gene map
    type TsvRec = (String, String);
    // now, map each transcript index to it's corresponding gene index
    let mut found = 0usize;
    // apparently the "header" (first row) will be included
    // in the iterator returned by `deserialize` anyway
    /*let hdr = rdr.headers()?;
    let hdr_vec : Vec<Result<TsvRec,csv::Error>> = vec![hdr.deserialize(None)];
    */
    for result in rdr.deserialize() {
        match result {
            Ok(record_in) => {
                let record: TsvRec = record_in;
                //let record: TSVRec = result?;
                // first, get the id for this gene
                let next_id = gene_name_to_id.len() as u32;
                let gene_id = *gene_name_to_id.entry(record.1.clone()).or_insert(next_id);
                // if we haven't added this gene name already, then
                // append it now to the list of gene names.
                if gene_id == next_id {
                    gene_names.push(record.1.clone());
                }
                // get the transcript id
                if let Some(transcript_id) = rname_to_id.get(&record.0) {
                    found += 1;
                    tid_to_gid[*transcript_id as usize] = gene_id;
                }
            }
            Err(e) => {
                /*
                crit!(
                log,
                "Encountered error [{}] when reading the transcript-to-gene map. Please make sure the transcript-to-gene mapping is a 2 column, tab separated file.",
                e
                );
                */
                return Err(anyhow!(
                    "failed to parse the transcript-to-gene map : {}.",
                    e
                ));
            }
        }
    }

    assert_eq!(
        found, ref_count,
        "The tg-map must contain a gene mapping for all transcripts in the header"
    );

    Ok((tid_to_gid, false))
}

pub fn parse_tg_map(
    tg_map: &PathBuf,
    ref_count: usize,
    rname_to_id: &HashMap<String, u32, ahash::RandomState>,
    gene_names: &mut Vec<String>,
    gene_name_to_id: &mut HashMap<String, u32, ahash::RandomState>,
) -> anyhow::Result<(Vec<u32>, bool)> {
    let t2g_file = std::fs::File::open(tg_map).context("couldn't open file")?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(t2g_file);

    let headers = rdr.headers()?;
    match headers.len() {
        2 => {
            // parse the 2 column format
            parse_tg_spliced(
                &mut rdr,
                ref_count,
                rname_to_id,
                gene_names,
                gene_name_to_id,
            )
        }
        3 => {
            // parse the 3 column format
            parse_tg_spliced_unspliced(
                &mut rdr,
                ref_count,
                rname_to_id,
                gene_names,
                gene_name_to_id,
            )
        }
        _ => {
            // not supported
            Err(anyhow!(
                "Transcript-gene mapping must have either 2 or 3 columns."
            ))
        }
    }
}

/// Extracts UMI counts from the `gene_eqc` HashMap.
/// This function is to be used when we are counting UMIs in
/// USA mode, and when we do not wish to consider gene-ambiguous
/// reads.
/// UMIs will be assigned to the spliced, unspliced, or ambiguous
/// version of their gene.  If a UMI is compatible with more than
/// one gene, but only one *spliced* gene, then it is assigned to
/// the spliced gene, unless there is too much multimapping
/// (i.e. it is compatible with > 10 different loci).
pub fn extract_counts(
    gene_eqc: &HashMap<Vec<u32>, u32, ahash::RandomState>,
    num_counts: usize,
) -> Vec<f32> {
    // the number of genes not considering status
    // i.e. spliced, unspliced, ambiguous
    let unspliced_offset = num_counts / 3;
    let ambig_offset = 2 * unspliced_offset;
    let mut counts = vec![0_f32; num_counts];

    for (labels, count) in gene_eqc {
        // the length of the label will tell us if this is a
        // splicing-unique, gene-unique (but splicing ambiguous).
        // or gene-ambiguous equivalence class label.
        match labels.len() {
            1 => {
                // determine if spliced or unspliced
                if let Some(gid) = labels.first() {
                    let idx = if is_spliced(*gid) {
                        (*gid >> 1) as usize
                    } else {
                        unspliced_offset + (*gid >> 1) as usize
                    };
                    counts[idx] += *count as f32;
                }
            }
            2 => {
                // spliced & unspliced of the same gene, or something differnet?
                if let (Some(g1), Some(g2)) = (labels.first(), labels.last()) {
                    if same_gene(*g1, *g2, true) {
                        let idx = ambig_offset + (*g1 >> 1) as usize;
                        //eprintln!("ambig count {} at {}!", *count, idx);
                        counts[idx] += *count as f32;
                    } else {
                        // report spliced if we can
                        match (is_spliced(*g1), is_spliced(*g2)) {
                            (true, false) => {
                                counts[(*g1 >> 1) as usize] += *count as f32;
                            }
                            (false, true) => {
                                counts[(*g2 >> 1) as usize] += *count as f32;
                            }
                            _ => { /* do nothing */ }
                        }
                    }
                }
            }
            3..=10 => {
                // if we don't have *too* many distinct genes matching this UMI
                // then apply the prefer-spliced rule.

                // See if there is precisely 1 spliced gene, and if so take it
                // but assign the read as ambiguous if it is for this gene
                let mut iter = labels.iter();
                // search for the first spliced index
                if let Some(sidx) = iter.position(|&x| is_spliced(x)) {
                    // if we found a spliced gene, check if there are any more
                    if let Some(_sidx2) = iter.position(|&x| is_spliced(x)) {
                        // in this case we had 2 spliced genes, so this is
                        // gene ambiguous and we just drop it.
                    } else {
                        // we only had one spliced gene.  Check to see if the
                        // index following the spliced gene we found is its
                        // unspliced variant or not.  If so, add it as ambiguous
                        // otherwise, add it as spliced
                        if let Some(sg) = labels.get(sidx) {
                            if let Some(ng) = labels.get(sidx + 1) {
                                if same_gene(*sg, *ng, true) {
                                    let idx = ambig_offset + (*sg >> 1) as usize;
                                    counts[idx] += *count as f32;
                                    continue;
                                }
                            }
                            counts[(*sg >> 1) as usize] += *count as f32;
                        }
                    }
                }
            }
            _ => {}
        }
    }
    counts
}

/// Extracts UMI counts from the `gene_eqc` HashMap.
/// This function is to be used when we are counting UMIs in
/// USA mode.  Multimappers will be uniformly allocated to the
/// genes to which they map.
pub fn extract_counts_mm_uniform(
    gene_eqc: &HashMap<Vec<u32>, u32, ahash::RandomState>,
    num_counts: usize,
) -> Vec<f32> {
    // the number of genes not considering status
    // i.e. spliced, unspliced, ambiguous
    let unspliced_offset = num_counts / 3;
    let ambig_offset = 2 * unspliced_offset;
    let mut counts = vec![0_f32; num_counts];
    let mut tvec = Vec::<usize>::with_capacity(16);

    for (labels, count) in gene_eqc {
        // the length of the label will tell us if this is a
        // splicing-unique, gene-unique (but splicing ambiguous).
        // or gene-ambiguous equivalence class label.
        match labels.len() {
            1 => {
                // determine if spliced or unspliced
                if let Some(gid) = labels.first() {
                    let idx = if is_spliced(*gid) {
                        (*gid >> 1) as usize
                    } else {
                        unspliced_offset + (*gid >> 1) as usize
                    };
                    counts[idx] += *count as f32;
                }
            }
            _ => {
                // iterate over all of the genes
                let mut iter = labels.iter().peekable();
                tvec.clear();
                while let Some(gn) = iter.next() {
                    // the base index of this gene
                    let mut idx = (gn >> 1) as usize;
                    // if the current gene is spliced
                    // check if the next item is the unspliced version
                    // of this gene.
                    if is_spliced(*gn) {
                        if let Some(ng) = iter.peek() {
                            // if this is the unspliced version
                            // of the same gene, then the count allocation
                            // goes to the ambiguous label
                            if same_gene(*gn, **ng, true) {
                                idx += ambig_offset;
                                // advance the iterator so we don't see
                                // this again.
                                iter.next();
                            }
                            // if it's not the same gene then add the
                            // contribution to the spliced molecule
                            // so do nothing here
                        }
                    } else {
                        // this is unspliced, so even if there is a next element
                        // it cannot belong to the same gene.
                        // modify the index so the contribution is
                        // to the unspliced gene index.
                        idx += unspliced_offset;
                    }
                    tvec.push(idx)
                }
                let fcount = (*count as f32) / (tvec.len() as f32);
                for g in &tvec {
                    counts[*g] += fcount;
                }
            }
        }
    }
    counts
}

/// Extracts an `IndexedEqList` and equivalence class ID / count
/// vector from the `gene_eqc` HashMap.
/// This function is used in USA-mode when we wish to resolve
/// multi-mapping UMIs via an EM algorithm. Equivalence class
/// labels (stored in `idx_eq_list`) will contain
/// spliced, unspliced and ambiguous gene IDs based on UMI mappings,
/// and `eq_id_count` will enumerate the count of UMIs for each
/// observed equivalence class.
pub fn extract_usa_eqmap(
    gene_eqc: &HashMap<Vec<u32>, u32, ahash::RandomState>,
    num_counts: usize,
    idx_eq_list: &mut IndexedEqList,
    eq_id_count: &mut Vec<(u32, u32)>,
) {
    // We use a little trick here.  Even though the resulting
    // USA-mode equivalence classes will not be over the same set
    // of gene IDs as the input list, we *do* know there will be
    // a 1-1 correspondence, such that each equivalence class label
    // in `gene_eqc` will produce exactly one USA-mode equivalence
    // class label, and that each USA-mode equivalence class label
    // will be unique.  This means we can just clear out our
    // `idx_eq_list` and add the new class labels and counts as we
    // encounter them.
    idx_eq_list.clear();
    eq_id_count.clear();

    // i.e. spliced, unspliced, ambiguous
    let unspliced_offset = num_counts / 3;
    let ambig_offset = 2 * unspliced_offset;
    let mut tvec = Vec::<u32>::with_capacity(16);

    for (ctr, (labels, count)) in gene_eqc.iter().enumerate() {
        // the length of the label will tell us if this is a
        // splicing-unique, gene-unique (but splicing ambiguous).
        // or gene-ambiguous equivalence class label.
        match labels.len() {
            1 => {
                // determine if spliced or unspliced
                if let Some(gid) = labels.first() {
                    let idx = if is_spliced(*gid) {
                        (*gid >> 1) as usize
                    } else {
                        unspliced_offset + (*gid >> 1) as usize
                    };
                    idx_eq_list.add_single_label(idx as u32);
                    eq_id_count.push((ctr as u32, *count));
                }
            }
            _ => {
                // iterate over all of the genes
                let mut iter = labels.iter().peekable();
                tvec.clear();
                while let Some(gn) = iter.next() {
                    // the base index of this gene
                    let mut idx = (gn >> 1) as usize;
                    // if the current gene is spliced
                    // check if the next item is the unspliced version
                    // of this gene.
                    if is_spliced(*gn) {
                        if let Some(ng) = iter.peek() {
                            // if this is the unspliced version
                            // of the same gene, then the count allocation
                            // goes to the ambiguous label
                            if same_gene(*gn, **ng, true) {
                                idx += ambig_offset;
                                // advance the iterator so we don't see
                                // this again.
                                iter.next();
                            }
                            // if it's not the same gene then add the
                            // contribution to the spliced molecule
                            // so do nothing here
                        }
                    } else {
                        // this is unspliced, so even if there is a next element
                        // it cannot belong to the same gene.
                        // modify the index so the contribution is
                        // to the unspliced gene index.
                        idx += unspliced_offset;
                    }
                    tvec.push(idx as u32);
                }
                // NOTE: the tvec won't necessarily be in sorted order
                // however, because we know the original eqc labels
                // and the USA mode labels are 1-1, we don't need this
                // so avoid the sort here.
                idx_eq_list.add_label_vec(tvec.as_slice());
                eq_id_count.push((ctr as u32, *count));
            }
        }
    }
}

pub fn get_bit_mask(nt_index: usize, fill_with: u64) -> u64 {
    let mut mask: u64 = fill_with;
    mask <<= 2 * (nt_index - 1);
    mask
}

pub fn get_all_snps(bc: u64, bc_length: usize) -> Vec<u64> {
    assert!(
        bc <= 2u64.pow(2 * bc_length as u32),
        "the barcode id is larger than possible (based on barcode length)"
    );
    assert!(
        bc_length <= 32,
        "barcode length greater than 32 not supported"
    );

    let mut snps: Vec<u64> = Vec::with_capacity(3 * bc_length);

    for nt_index in 1..=bc_length {
        // clearing the two relevant bits based on nucleotide position
        let bit_mask = bc & !get_bit_mask(nt_index, 3);

        // iterating over all 4 choices of the nucleotide
        for i in 0..=3 {
            let new_bc = bit_mask | get_bit_mask(nt_index, i);
            if new_bc != bc {
                snps.push(new_bc);
            }
        }
    }

    snps
}

pub fn get_all_indels(bc: u64, bc_length: usize) -> Vec<u64> {
    assert!(
        bc <= 2u64.pow(2 * bc_length as u32),
        "the barcode id is larger than possible (based on barcode length)"
    );
    assert!(
        bc_length <= 32,
        "barcode length greater than 32 not supported"
    );

    let mut indels: Vec<u64> = Vec::with_capacity(8 * (bc_length - 1));

    for nt_index in 1..bc_length {
        let mut bit_mask = 1 << (2 * nt_index);
        bit_mask -= 1;

        let upper_half = bc & !bit_mask;
        let lower_half = bc & bit_mask;

        // iterating over all 4 choices of the nucleotide
        for i in 0..=3 {
            let new_insertion_bc = upper_half | get_bit_mask(nt_index, i) | (lower_half >> 2);
            let new_deletion_bc = upper_half
                | get_bit_mask(1, i)
                | ((lower_half & !get_bit_mask(nt_index + 1, 3)) << 2);

            if new_insertion_bc != bc {
                indels.push(new_insertion_bc);
            }
            if new_deletion_bc != bc {
                indels.push(new_deletion_bc);
            }
        }
    }

    indels
}

pub fn get_all_one_edit_neighbors(
    bc: u64,
    bc_length: usize,
    neighbors: &mut HashSet<u64>,
) -> Result<(), Box<dyn Error>> {
    neighbors.clear();

    let snps: Vec<u64> = get_all_snps(bc, bc_length);
    let indels: Vec<u64> = get_all_indels(bc, bc_length);

    neighbors.extend(&snps);
    neighbors.extend(&indels);

    Ok(())
}

pub fn generate_whitelist_set(
    whitelist_bcs: &[u64],
    bc_length: usize,
) -> Result<HashSet<u64>, Box<dyn Error>> {
    let num_bcs = whitelist_bcs.len();

    let mut one_edit_barcode_hash: HashSet<u64> = HashSet::new();
    let mut neighbors: HashSet<u64> = HashSet::new();
    one_edit_barcode_hash.reserve(10 * num_bcs);
    // reserved space for 3*length SNP
    // + 4 * (length -1) insertion
    // + 4 * (length -1) deletion
    neighbors.reserve(3 * bc_length + 8 * (bc_length - 1));

    for bc in whitelist_bcs {
        get_all_one_edit_neighbors(*bc, bc_length, &mut neighbors)?;
        one_edit_barcode_hash.extend(&neighbors);
    }

    Ok(one_edit_barcode_hash)
}

/**
 * generates a map that contains all one edit distance neighbors
 * of the permitted barcodes.  The key is the neighbor and the value
 * is the original permitted barcode to which it maps.
 **/
pub fn generate_permitlist_map(
    permit_bcs: &[u64],
    bc_length: usize,
) -> Result<HashMap<u64, u64>, Box<dyn Error>> {
    let num_bcs = permit_bcs.len();

    let mut one_edit_barcode_map: HashMap<u64, u64> = HashMap::with_capacity(10 * num_bcs);
    // first insert everything already in the explicit permitlist
    for bc in permit_bcs {
        one_edit_barcode_map.insert(*bc, *bc);
    }

    // reserved space for 3*length SNP
    // + 4 * (length -1) insertion
    // + 4 * (length -1) deletion
    let mut neighbors: HashSet<u64> = HashSet::with_capacity(3 * bc_length + 8 * (bc_length - 1));

    for bc in permit_bcs {
        get_all_one_edit_neighbors(*bc, bc_length, &mut neighbors)?;
        for n in &neighbors {
            one_edit_barcode_map.entry(*n).or_insert(*bc);
        }
    }

    Ok(one_edit_barcode_map)
}

/// Reads the contents of the file `flist`, which should contain
/// a single barcode per-line, and returns a Result that is either
/// a HashSet containing the k-mer encoding of all barcodes or
/// the Error that was encountered parsing the file.
pub fn read_filter_list(
    flist: &PathBuf,
    bclen: u16,
) -> anyhow::Result<HashSet<u64, ahash::RandomState>> {
    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut fset = HashSet::<u64, ahash::RandomState>::with_hasher(s);

    let filt_file = std::fs::File::open(flist).context("couldn't open file")?;
    let mut reader = BufReader::new(filt_file);

    // Read the file line by line using the lines() iterator from std::io::BufRead.
    reader
        .for_byte_line(|line| {
            let mut bnk = BitNuclKmer::new(line, bclen as u8, false);
            let (_, k, _) = bnk.next().unwrap();
            fset.insert(k.0);
            Ok(true)
        })
        .unwrap();

    Ok(fset)
}

pub fn is_velo_mode(input_dir: &PathBuf) -> bool {
    let parent = std::path::Path::new(input_dir);
    // open the metadata file and read the json
    let meta_data_file = File::open(parent.join("generate_permit_list.json"))
        .expect("could not open the generate_permit_list.json file.");
    let mdata: serde_json::Value = serde_json::from_reader(meta_data_file)
        .expect("could not deseralize generate_permit_list.json");
    let vm = mdata.get("velo_mode");
    match vm {
        Some(v) => v.as_bool().unwrap_or(false),
        None => false,
    }
}

#[allow(dead_code)]
#[derive(Debug, PartialEq, Eq)]
pub struct InternalVersionInfo {
    pub major: u32,
    pub minor: u32,
    pub patch: u32,
}

impl InternalVersionInfo {
    pub fn is_compatible_with(&self, other: &InternalVersionInfo) -> Result<(), String> {
        if self.major == other.major && self.minor == other.minor {
            Ok(())
        } else {
            let s = format!(
                "running alevin-fry {} on {} results, please regenerate the results using alevin-fry {} or greater",
                self, other, self
            );
            Err(s)
        }
    }
}

impl fmt::Display for InternalVersionInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "v{}.{}.{}", self.major, self.minor, self.patch)
    }
}

#[derive(Error, Debug)]
pub enum VersionParseError {
    #[error("The version string should be of the format x.y.z; it was `{0}`")]
    IncorrectFormat(String),
}

impl FromStr for InternalVersionInfo {
    type Err = VersionParseError;

    fn from_str(vs: &str) -> Result<Self, Self::Err> {
        let versions: Vec<u32> = vs.split('.').map(|s| s.parse::<u32>().unwrap()).collect();
        if versions.len() != 3 {
            return Err(VersionParseError::IncorrectFormat(vs.to_string()));
        }
        Ok(Self {
            major: versions[0],
            minor: versions[1],
            patch: versions[2],
        })
    }
}

pub fn get_bc_string(
    kmerseq: &needletail::bitkmer::BitKmerSeq,
    reverse_barcode: bool,
    bc_len: u8,
) -> String {
    let kmseq = *kmerseq;
    let mut km: needletail::bitkmer::BitKmer = (kmseq, bc_len);
    if reverse_barcode {
        km = needletail::bitkmer::reverse_complement(km);
    }
    let bytes = needletail::bitkmer::bitmer_to_bytes(km);
    let seq: String = String::from_utf8(bytes).expect("Invalid barcode");
    seq
}

#[cfg(test)]
mod tests {
    use crate::utils::generate_whitelist_set;
    use crate::utils::get_all_indels;
    use crate::utils::get_all_one_edit_neighbors;
    use crate::utils::get_all_snps;
    use crate::utils::get_bit_mask;
    use crate::utils::InternalVersionInfo;
    use std::collections::HashSet;
    use std::str::FromStr;

    #[test]
    fn test_version_info() {
        let vi = InternalVersionInfo::from_str("1.2.3").unwrap();
        assert_eq!(
            vi,
            InternalVersionInfo {
                major: 1,
                minor: 2,
                patch: 3
            }
        );
    }

    #[test]
    fn test_get_bit_mask() {
        let mut output = Vec::new();
        for i in 0..=3 {
            let mask = get_bit_mask(2, i);
            output.push(mask);
        }
        assert_eq!(output, vec![0, 4, 8, 12]);
    }

    #[test]
    fn test_get_all_snps() {
        let mut output: Vec<u64> = get_all_snps(7, 3).into_iter().collect();
        output.sort_unstable();

        assert_eq!(output, vec![3, 4, 5, 6, 11, 15, 23, 39, 55]);
    }

    #[test]
    fn test_get_all_indels() {
        let mut output: Vec<u64> = get_all_indels(7, 3).into_iter().collect();
        output.sort_unstable();
        output.dedup();

        assert_eq!(output, vec![1, 4, 5, 6, 9, 12, 13, 14, 15, 28, 29, 30, 31]);
    }

    #[test]
    fn test_get_all_one_edit_neighbors() {
        let mut neighbors: HashSet<u64> = HashSet::new();
        get_all_one_edit_neighbors(7, 3, &mut neighbors).unwrap();

        let mut output: Vec<u64> = neighbors.into_iter().collect();

        output.sort_unstable();
        output.dedup();

        assert_eq!(
            output,
            vec![1, 3, 4, 5, 6, 9, 11, 12, 13, 14, 15, 23, 28, 29, 30, 31, 39, 55]
        );
    }

    #[test]
    fn test_generate_whitelist_hash() {
        let neighbors: HashSet<u64> = generate_whitelist_set(&[7], 3).unwrap();
        let mut output: Vec<u64> = neighbors.into_iter().collect();

        output.sort_unstable();
        output.dedup();

        assert_eq!(
            output,
            vec![1, 3, 4, 5, 6, 9, 11, 12, 13, 14, 15, 23, 28, 29, 30, 31, 39, 55]
        );
    }
}
