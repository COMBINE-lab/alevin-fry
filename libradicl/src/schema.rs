// Copyright 2020 Rob Patro, Avi Srivastava. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

extern crate bio_types;
extern crate fasthash;
extern crate quickersort;

use crate as libradicl;
use bio_types::strand::Strand;
use fasthash::{sea, sea::Hash64, RandomState};
use std::collections::HashMap;
use std::fmt;
use std::io::BufRead;
use std::str::FromStr;

/**
* Single-cell equivalence class
**/
#[derive(Debug)]
pub struct CellEQClass<'a> {
    // transcripts defining this eq. class
    pub transcripts: &'a Vec<u32>,
    // umis with multiplicities
    // the k-mer should be a k-mer class eventually
    pub umis: Vec<(u64, u32)>,
}

#[derive(Debug)]
pub(super) struct EqMapEntry {
    pub umis: Vec<(u64, u32)>,
    pub eq_num: u32,
}

#[derive(Debug)]
pub struct ProtocolInfo {
    // TODO: only makes sense
    // for single-strand protocols
    // right now.  Expand to be generic.
    pub expected_ori: Strand,
}

#[derive(PartialEq, Debug, Clone, Copy)]
pub enum ResolutionStrategy {
    Trivial,
    CellRangerLike,
    CellRangerLikeEM,
    Full,
    Parsimony,
}

impl fmt::Display for ResolutionStrategy {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

// Implement the trait
impl FromStr for ResolutionStrategy {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "trivial" => Ok(ResolutionStrategy::Trivial),
            "cr-like" => Ok(ResolutionStrategy::CellRangerLike),
            "cr-like-em" => Ok(ResolutionStrategy::CellRangerLikeEM),
            "full" => Ok(ResolutionStrategy::Full),
            "parsimony" => Ok(ResolutionStrategy::Parsimony),
            _ => Err("no match"),
        }
    }
}

// NOTE: this is _clearly_ redundant with the EqMap below.
// we should re-factor so that EqMap makes use of this class
// rather than replicates its members
pub(super) struct IndexedEqList {
    pub num_genes: usize,
    // the total size of the list of all reference
    // ids over all equivalence classes that occur in this
    // cell
    pub label_list_size: usize,
    // concatenated lists of the labels of all equivalence classes
    pub eq_labels: Vec<u32>,
    // vector that deliniates where each equivalence class label
    // begins and ends.  The label for equivalence class i begins
    // at offset eq_label_starts[i], and it ends at
    // eq_label_starts[i+1].  The length of this vector is 1 greater
    // than the number of equivalence classes.
    pub eq_label_starts: Vec<u32>,
}

impl IndexedEqList {
    /// returns the number of equivalence classes represented in this `IndexedEqList`
    pub(super) fn num_eq_classes(&self) -> usize {
        if !self.eq_label_starts.is_empty() {
            self.eq_label_starts.len() - 1
        } else {
            0usize
        }
    }

    //pub(super) fn num_genes(&self) -> usize {
    //    self.num_genes;
    //}

    /// clears out the contents of this `IndexedEqList`
    #[allow(dead_code)]
    pub(super) fn clear(&mut self) {
        self.label_list_size = 0usize;
        self.eq_labels.clear();
        self.eq_label_starts.clear();
    }

    /// Creates an `IndexedEqList` from a HashMap of eq labels to counts
    pub(super) fn init_from_hash(
        eqclasses: &HashMap<Vec<u32>, u32, fasthash::RandomState<Hash64>>,
        num_genes: usize,
    ) -> IndexedEqList {
        let num_eqc = eqclasses.len();

        let mut label_list_size = 0usize;

        // concatenated lists of the labels of all equivalence classes
        let mut eq_labels = Vec::<u32>::with_capacity(2 * num_eqc);

        // vector that deliniates where each equivalence class label
        // begins and ends.  The label for equivalence class i begins
        // at offset eq_label_starts[i], and it ends at
        // eq_label_starts[i+1].  The length of this vector is 1 greater
        // than the number of equivalence classes.
        let mut eq_label_starts = Vec::<u32>::with_capacity(num_eqc + 1);

        eq_label_starts.push(0);
        for (labels, _count) in eqclasses.iter() {
            label_list_size += labels.len();
            eq_labels.extend(labels.clone());
            eq_label_starts.push(eq_labels.len() as u32);
        }

        IndexedEqList {
            num_genes,
            label_list_size,
            eq_labels,
            eq_label_starts,
        }
    }

    /// Loads the `IndexedEqList` from a gzip compressed file
    pub(super) fn init_from_eqc_file<P: AsRef<std::path::Path>>(eqc_path: P) -> IndexedEqList {
        let file = flate2::read::GzDecoder::new(std::fs::File::open(eqc_path).unwrap());
        let reader = std::io::BufReader::new(file);

        let mut lit = reader.lines();

        let num_genes = lit.next().unwrap().unwrap().parse::<u32>().unwrap() as usize;
        let num_eqc = lit.next().unwrap().unwrap().parse::<u32>().unwrap();

        let mut eqid_map: Vec<Vec<u32>> = vec![vec![]; num_eqc as usize];

        let mut label_list_size = 0usize;

        // go over all other lines and fill in our equivalence class map
        for l in lit {
            let v = l
                .unwrap()
                .split_whitespace()
                .map(|y| y.parse::<u32>().unwrap())
                .collect::<Vec<u32>>();
            if let Some((id, ev)) = v.split_last() {
                label_list_size += ev.len();
                eqid_map[*id as usize] = ev.to_vec();
            }
        }

        // concatenated lists of the labels of all equivalence classes
        let mut eq_labels = Vec::<u32>::with_capacity(label_list_size);

        // vector that deliniates where each equivalence class label
        // begins and ends.  The label for equivalence class i begins
        // at offset eq_label_starts[i], and it ends at
        // eq_label_starts[i+1].  The length of this vector is 1 greater
        // than the number of equivalence classes.
        let mut eq_label_starts = Vec::<u32>::with_capacity(eqid_map.len() + 1);

        eq_label_starts.push(0);
        // now we have everything in order, we need to flatten it
        for v in eqid_map.iter() {
            eq_labels.extend(v);
            eq_label_starts.push(eq_labels.len() as u32);
        }

        IndexedEqList {
            num_genes,
            label_list_size,
            eq_labels,
            eq_label_starts,
        }
    }

    /// Returns a slice of gene IDs corresponding to the
    /// label for equivalence class `idx`.
    pub(super) fn refs_for_eqc(&self, idx: u32) -> &[u32] {
        &self.eq_labels[(self.eq_label_starts[idx as usize] as usize)
            ..(self.eq_label_starts[(idx + 1) as usize] as usize)]
    }
}

pub(super) struct EqMap {
    // for each equivalence class, holds the (umi, freq) pairs
    // and the id of that class
    pub eqc_info: Vec<EqMapEntry>,
    // the total number of refrence targets
    pub nref: u32,
    // the total size of the list of all reference
    // ids over all equivalence classes that occur in this
    // cell
    pub label_list_size: usize,
    // concatenated lists of the labels of all equivalence classes
    pub eq_labels: Vec<u32>,
    // vector that deliniates where each equivalence class label
    // begins and ends.  The label for equivalence class i begins
    // at offset eq_label_starts[i], and it ends at
    // eq_label_starts[i+1].  The length of this vector is 1 greater
    // than the number of equivalence classes.
    pub eq_label_starts: Vec<u32>,
    // the number of equivalence class labels in which each reference
    // appears
    pub label_counts: Vec<u32>,
    // the offset into the global list of equivalence class ids
    // where the equivalence class list for each reference starts
    pub ref_offsets: Vec<u32>,
    // the concatenated list of all equivalence class labels for
    // all transcripts
    pub ref_labels: Vec<u32>,
}

impl EqMap {
    pub(super) fn num_eq_classes(&self) -> usize {
        self.eqc_info.len()
    }

    #[allow(dead_code)]
    pub(super) fn clear(&mut self) {
        self.eqc_info.clear();
        // keep nref
        self.label_list_size = 0usize;
        self.eq_labels.clear();
        self.eq_label_starts.clear();
        // clear the label_counts, but resize
        // and fill with 0
        self.label_counts.clear();
        self.label_counts.resize(self.nref as usize, 0u32);

        self.ref_offsets.clear();
        self.ref_labels.clear();
    }

    pub(super) fn new(nref_in: u32) -> EqMap {
        EqMap {
            eqc_info: vec![], //HashMap::with_hasher(rs),
            nref: nref_in,
            label_list_size: 0usize,
            eq_labels: vec![],
            eq_label_starts: vec![],
            label_counts: vec![0; nref_in as usize],
            ref_offsets: vec![],
            ref_labels: vec![],
        }
    }

    pub(super) fn fill_ref_offsets(&mut self) {
        self.ref_offsets = self
            .label_counts
            .iter()
            .scan(0, |sum, i| {
                *sum += i;
                Some(*sum)
            })
            .collect::<Vec<_>>();
        self.ref_offsets.push(*self.ref_offsets.last().unwrap());
    }

    pub(super) fn fill_label_sizes(&mut self) {
        self.ref_labels = vec![u32::MAX; self.label_list_size + 1];
    }

    pub(super) fn init_from_chunk(&mut self, cell_chunk: &mut libradicl::Chunk) {
        // temporary map of equivalence class label to assigned
        // index.
        let s = RandomState::<sea::Hash64>::new();
        let mut eqid_map: HashMap<Vec<u32>, u32, fasthash::RandomState<fasthash::sea::Hash64>> =
            HashMap::with_hasher(s);

        // gather the equivalence class info
        for r in &mut cell_chunk.reads {
            // TODO: ensure this is done upstream so we
            // don't have to do it here.
            // NOTE: should be done if collate was run.
            // r.refs.sort();

            match eqid_map.get_mut(&r.refs) {
                // if we've seen this equivalence class before, just add the new
                // umi.
                Some(v) => {
                    self.eqc_info[*v as usize].umis.push((r.umi, 1));
                }
                // otherwise, add the new umi, but we also have some extra bookkeeping
                None => {
                    // each reference in this equivalence class label
                    // will have to point to this equivalence class id
                    let eq_num = self.eqc_info.len() as u32;
                    self.label_list_size += r.refs.len();
                    for r in r.refs.iter() {
                        let ridx = *r as usize;
                        self.label_counts[ridx] += 1;
                        //ref_to_eqid[*r as usize].push(eq_num);
                    }
                    self.eq_label_starts.push(self.eq_labels.len() as u32);
                    self.eq_labels.extend(&r.refs);
                    self.eqc_info.push(EqMapEntry {
                        umis: vec![(r.umi, 1)],
                        eq_num,
                    });
                    eqid_map.insert(r.refs.clone(), eq_num);
                    //self.eqc_map.insert(r.refs.clone(), EqMapEntry { umis : vec![(r.umi,1)], eq_num });
                }
            }
        }
        // final value to avoid special cases
        self.eq_label_starts.push(self.eq_labels.len() as u32);

        self.fill_ref_offsets();
        self.fill_label_sizes();

        // initially we inserted duplicate UMIs
        // here, collapse them and keep track of their count
        for idx in 0..self.num_eq_classes() {
            //} self.eqc_info.iter_mut().enumerate() {

            // for each reference in this
            // label, put it in the next free spot
            // TODO: @k3yavi, can we avoid this copy?
            let label = self.refs_for_eqc(idx as u32).to_vec();
            //println!("{:?}", label);
            for r in label {
                // k.iter() {
                self.ref_offsets[r as usize] -= 1;
                self.ref_labels[self.ref_offsets[r as usize] as usize] = self.eqc_info[idx].eq_num;
            }

            let v = &mut self.eqc_info[idx];
            // sort so dups are adjacent
            quickersort::sort(&mut v.umis[..]);
            // we need a copy of the vector b/c we
            // can't easily modify it in place
            // at least I haven't seen how (@k3yavi, help here if you can).
            let cv = v.umis.clone();
            // since we have a copy, clear the original to fill it
            // with the new contents.
            v.umis.clear();

            let mut count = 1;
            let mut cur_elem = cv.first().unwrap().0;
            for e in cv.iter().skip(1) {
                if e.0 == cur_elem {
                    count += 1;
                } else {
                    v.umis.push((cur_elem, count));
                    cur_elem = e.0;
                    count = 1;
                }
            }

            // remember to push the last element, since we
            // won't see a subsequent "different" element.
            v.umis.push((cur_elem, count));
        }
    }

    pub(super) fn eq_classes_containing(&self, r: u32) -> &[u32] {
        &self.ref_labels
            [(self.ref_offsets[r as usize] as usize)..(self.ref_offsets[(r + 1) as usize] as usize)]
    }

    pub(super) fn refs_for_eqc(&self, idx: u32) -> &[u32] {
        &self.eq_labels[(self.eq_label_starts[idx as usize] as usize)
            ..(self.eq_label_starts[(idx + 1) as usize] as usize)]
    }
}

#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub(super) enum SplicedStatus {
    Unspliced,
    Spliced,
    Ambiguous,
}

impl fmt::Display for SplicedStatus {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            SplicedStatus::Unspliced => write!(f, "U"),
            SplicedStatus::Spliced => write!(f, "S"),
            SplicedStatus::Ambiguous => write!(f, "A"),
        }
    }
}

#[derive(Debug)]
pub(super) enum PUGEdgeType {
    NoEdge,
    BiDirected,
    XToY,
    YToX,
}

#[derive(Debug)]
pub(super) struct PUGResolutionStatistics {
    pub used_alternative_strategy: bool,
    pub total_mccs: u64,
    pub ambiguous_mccs: u64,
    pub trivial_mccs: u64,
}

pub struct VeloCounts {
    pub spliced: Vec<f32>,
    pub unspliced: Vec<f32>,
}

impl VeloCounts {
    pub(super) fn new(num_genes: usize, init_value: f32) -> VeloCounts {
        VeloCounts {
            spliced: vec![init_value; num_genes],
            unspliced: vec![init_value; num_genes],
        }
    }

    pub(super) fn scale(&mut self) {
        self.spliced.iter_mut().for_each(|alpha| *alpha *= 1e-3);
        self.unspliced.iter_mut().for_each(|alpha| *alpha *= 1e-3);
    }

    pub(super) fn only_unique(&mut self) {
        self.spliced.iter_mut().for_each(|alpha| {
            *alpha -= 0.5;
        });
        self.unspliced.iter_mut().for_each(|alpha| {
            *alpha -= 0.5;
        });
    }

    #[allow(dead_code)]
    pub(super) fn unique_initialize(&mut self, gid: usize, count: f32, cp: f32) {
        self.spliced[gid] += count * cp;
        self.unspliced[gid] += count * (1.0 - cp);
    }

    pub(super) fn update_small_alpha(&mut self, min_alpha: f32) {
        self.spliced.iter_mut().for_each(|alpha| {
            if *alpha < min_alpha {
                *alpha = 0.0_f32;
            }
        });
        self.unspliced.iter_mut().for_each(|alpha| {
            if *alpha < min_alpha {
                *alpha = 0.0_f32;
            }
        });
    }

    pub(super) fn sum(&self) -> f32 {
        let spliced_sum: f32 = self.spliced.iter().sum();
        let unspliced_sum: f32 = self.unspliced.iter().sum();
        spliced_sum + unspliced_sum
    }

    pub(super) fn multigene_add(
        &mut self,
        alphas_in: &VeloCounts,
        gid: usize,
        inv_denominator: f32,
        cp: f32,
    ) {
        const ERROR_MARGIN: f32 = f32::EPSILON; // Use an epsilon for comparison
        let quota = alphas_in.gene_alpha(gid) * inv_denominator;
        if (cp - 1.0).abs() < ERROR_MARGIN {
            self.spliced[gid] += quota;
        } else if cp == 0.0 {
            self.unspliced[gid] += quota;
        } else {
            let spliced_proportion = (alphas_in.spliced[gid] * cp)
                / (alphas_in.spliced[gid] * cp + alphas_in.unspliced[gid] * (1.0 - cp));
            if !spliced_proportion.is_nan() {
                self.spliced[gid] += quota * spliced_proportion;
                self.unspliced[gid] += quota * (1.0 - spliced_proportion);
            }
        }
        // let spliced_proportion = (alphas_in.spliced[gid] * cp) / (alphas_in.spliced[gid] * cp +  alphas_in.unspliced[gid] * (1.0 - cp));
        // let quota = alphas_in.gene_alpha(gid) * inv_denominator;
        // if (quota * spliced_proportion).is_nan() {
        //     println!{"\n\nspliced: {}, unspliced: {}, gene_alpha:{} cp: {}, inv_den: {}, quota: {}", alphas_in.spliced[gid], alphas_in.unspliced[gid], alphas_in.gene_alpha(gid), cp, inv_denominator, quota};
        // }
        // self.spliced[gid] += quota * spliced_proportion;
        // self.unspliced[gid] += quota * (1.0 - spliced_proportion);
    }

    pub(super) fn unigene_add(&mut self, alphas_in: &VeloCounts, gid: usize, count: f32, cp: f32) {
        const ERROR_MARGIN: f32 = f32::EPSILON; // Use an epsilon for comparison
        if (cp - 1.0).abs() < ERROR_MARGIN {
            self.spliced[gid] += count;
        } else if cp == 0.0 {
            self.unspliced[gid] += count;
        } else {
            let spliced_proportion = alphas_in.spliced[gid] * cp
                / (alphas_in.spliced[gid] * cp + alphas_in.unspliced[gid] * (1.0 - cp));
            if !spliced_proportion.is_nan() {
                self.spliced[gid] += count * spliced_proportion;
                self.unspliced[gid] += count * (1.0 - spliced_proportion);
            }
        }

        // let spliced_proportion = alphas_in.spliced[gid] * cp / (alphas_in.spliced[gid] * cp +  alphas_in.unspliced[gid] * (1.0 - cp));
        // if (count * spliced_proportion).is_nan() {
        //     println!{"\n\nspliced: {}, unspliced: {}, gene_alpha:{} cp: {}, count: {}", alphas_in.spliced[gid], alphas_in.unspliced[gid], alphas_in.gene_alpha(gid), cp, count};
        // }
        // self.spliced[gid] += count * spliced_proportion;
        // self.unspliced[gid] += count * (1.0 - spliced_proportion);
    }

    pub(super) fn update(&mut self, alphas_out: &mut VeloCounts, gid: usize) {
        self.spliced[gid] = alphas_out.spliced[gid];
        self.unspliced[gid] = alphas_out.unspliced[gid];
        alphas_out.spliced[gid] = 0.0;
        alphas_out.unspliced[gid] = 0.0;
    }

    pub(super) fn diff(&self, alphas_out: &VeloCounts, gid: usize) -> f32 {
        self.gene_alpha(gid) - alphas_out.gene_alpha(gid)
    }

    #[allow(dead_code)]
    pub(super) fn len(&self) -> usize {
        self.spliced.len()
    }
    pub(super) fn gene_alpha(&self, gid: usize) -> f32 {
        self.spliced[gid] + self.unspliced[gid]
    }
    #[allow(dead_code)]
    pub(super) fn clone(&self) -> VeloCounts {
        VeloCounts {
            spliced: self.spliced.clone(),
            unspliced: self.unspliced.clone(),
        }
    }
    // pub(super) fn iter(&self) -> Counts {
    //     Counts {
    //         spliced: self.spliced.clone(),
    //         unspliced: self.unspliced.clone(),
    //     }
    // }
}
