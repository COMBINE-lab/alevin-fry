/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */
use std::collections::HashMap;
use std::hash::{BuildHasher, Hasher};
use std::io::BufRead;

use libradicl::chunk;
use libradicl::record::{MappedRecord, UmiTaggedRecord};

use crate::utils::OptionalAlignmentScores;

fn score_probabilities(scores: &[i32]) -> Vec<f64> {
    const DENOM: f64 = 10.0;
    let max_score = *scores.iter().max().unwrap();

    let mut out: Vec<f64> = scores
        .iter()
        .map(|&s| (((s - max_score) as f64) / DENOM).exp())
        .collect();
    let sum: f64 = out.iter().sum();
    if sum > 0.0 {
        for p in out.iter_mut() {
            *p /= sum;
        }
    } else {
        // raise some error!
    }
    out
}

// Modified from https://stackoverflow.com/questions/69764050/how-to-get-the-indices-that-would-sort-a-vec
// kmdreko
fn argsort<T: Ord>(data: &[T]) -> Vec<usize> {
    let mut indices = (0..data.len()).collect::<Vec<_>>();
    indices.sort_unstable_by_key(|&i| &data[i]);
    indices
}
/// initially suggested by Claude
fn argsort_by<T, F>(data: &[T], mut compare: F) -> Vec<usize>
where
    F: FnMut(&T, &T) -> std::cmp::Ordering,
{
    let mut indices: Vec<usize> = (0..data.len()).collect();
    indices.sort_unstable_by(|&i, &j| compare(&data[i], &data[j]));
    indices
}

/// Reorder a vector in-place using the given permutation indices.
/// This is more memory-efficient but modifies the original vector.
/// Time: O(n), Space: O(n) for tracking visited indices.
fn reorder_in_place<T>(data: &mut [T], indices: &[usize]) {
    let mut visited = vec![false; data.len()];

    for start in 0..data.len() {
        if visited[start] {
            continue;
        }

        let mut current = start;
        let mut next = indices[current];

        while next != start {
            visited[current] = true;
            data.swap(current, next);
            current = next;
            next = indices[next];
        }
        visited[current] = true;
    }
}

/**
* Single-cell equivalence class
**/
#[derive(Debug)]
pub struct CellEqClass<'a> {
    // transcripts defining this eq. class
    pub transcripts: &'a Vec<u32>,
    // umis with multiplicities
    // the k-mer should be a k-mer class eventually
    pub umis: Vec<(u64, u32)>,
}

#[derive(Debug, Clone, Copy)]
pub enum EqMapType {
    TranscriptLevel,
    GeneLevel,
}

#[derive(Debug)]
pub struct EqMapEntry {
    pub umis: Vec<(u64, u32)>,
    pub eq_num: u32,
}
// NOTE: this is _clearly_ redundant with the EqMap below.
// we should re-factor so that EqMap makes use of this class
// rather than replicates its members
pub struct IndexedEqList {
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
    pub fn new() -> Self {
        IndexedEqList {
            num_genes: 0usize,
            label_list_size: 0usize,
            eq_labels: Vec::<u32>::new(),
            eq_label_starts: Vec::<u32>::new(),
        }
    }

    /// returns the number of equivalence classes represented in this `IndexedEqList`
    pub fn num_eq_classes(&self) -> usize {
        if !self.eq_label_starts.is_empty() {
            self.eq_label_starts.len() - 1
        } else {
            0usize
        }
    }

    pub fn add_single_label(&mut self, lab: u32) {
        self.label_list_size += 1;
        self.eq_labels.push(lab);
        self.eq_label_starts.push(self.eq_labels.len() as u32);
    }

    pub fn add_label_vec(&mut self, lab: &[u32]) {
        self.label_list_size += lab.len();
        self.eq_labels.extend_from_slice(lab);
        self.eq_label_starts.push(self.eq_labels.len() as u32);
    }

    //pub(super) fn num_genes(&self) -> usize {
    //    self.num_genes;
    //}

    /// clears out the contents of this `IndexedEqList`
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.label_list_size = 0usize;
        self.eq_labels.clear();
        self.eq_label_starts.clear();
        self.eq_label_starts.push(0);
    }

    /// Creates an `IndexedEqList` from a HashMap of eq labels to counts
    pub fn init_from_hash(
        eqclasses: &HashMap<Vec<u32>, u32, ahash::RandomState>,
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
    pub fn init_from_eqc_file<P: AsRef<std::path::Path>>(eqc_path: P) -> IndexedEqList {
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
    pub fn refs_for_eqc(&self, idx: u32) -> &[u32] {
        &self.eq_labels[(self.eq_label_starts[idx as usize] as usize)
            ..(self.eq_label_starts[(idx + 1) as usize] as usize)]
    }
}

impl Default for IndexedEqList {
    fn default() -> Self {
        Self::new()
    }
}

pub struct EqMap {
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
    eqid_map: HashMap<Vec<u32>, u32, ahash::RandomState>,
    pub map_type: EqMapType,
    pub prob_map: Option<ProbMap>,
}

pub(crate) struct TempProbMap {
    pub eq_ids: Vec<u32>, // for each score chunk the equivalence class ID to which it belongs
    pub probs: Vec<f64>,  // the concatenated list of probabilities
    pub lengths: Vec<usize>, // the lengths of each prob slice
}

pub(crate) struct ProbMap {
    pub probs: Vec<f64>, // the concatenated list of probabilities
    /// access the probabilities for all reads (can be more than 1) corresponding
    /// to the i'th global UMI (which belongs to some equivalence class determined
    /// by the rank in eq_indices).
    pub umi_offsets: Vec<usize>,
    /// the cumulative indices in the umi_offsets vector that start each new
    /// equivalence classes' set of umis
    pub eq_indices: Vec<usize>,
}

impl ProbMap {
    pub fn new() -> Self {
        Self {
            probs: Vec::new(),
            umi_offsets: vec![0],
            eq_indices: vec![0],
        }
    }

    pub fn mark_umi_end(&mut self) {
        self.umi_offsets.push(self.probs.len());
    }

    pub fn mark_eq_class_end(&mut self) {
        self.eq_indices.push(self.umi_offsets.len());
    }

    pub fn num_umis_for_eq(&self, eq_id: usize) -> Option<usize> {
        if eq_id + 1 < self.eq_indices.len() {
            let end_offset = self.eq_indices[eq_id + 1];
            let start_offset = self.eq_indices[eq_id];
            Some(end_offset - start_offset)
        } else {
            None
        }
    }

    pub fn probs_for_eq_id_umi_rank(&self, eq_id: usize, umi_rank: usize) -> &[f64] {
        let eq_offset = self.eq_indices[eq_id];
        let start_prob_offsets = self.umi_offsets[eq_offset + umi_rank];
        let end_prob_offsets = self.umi_offsets[eq_offset + umi_rank + 1];
        &self.probs[start_prob_offsets..end_prob_offsets]
    }
}

pub(crate) struct EqClassAlnProbView<'a> {
    pub eq_id: u32,               // the id of this equivalence class
    pub probs: &'a [f64],         // the probabilities
    cumulative_offsets: Vec<u32>, // where the per-read offsets start
}

impl<'a> EqClassAlnProbView<'a> {
    // the number of reads in this eq class for which we have
    // information
    pub fn num_reads(&self) -> usize {
        self.cumulative_offsets.len() - 1
    }

    // get the probabilities for the read of the associated rank within
    // this equivalence class
    pub fn get_probs_for_read_rank(&self, r: usize) -> Option<&[f64]> {
        // if the rank
        if r + 1 > self.num_reads() {
            None
        } else {
            let start = self.cumulative_offsets[r] as usize;
            let end = self.cumulative_offsets[r + 1] as usize;
            Some(&self.probs[start..end])
        }
    }

    pub fn get_eq_id(&self) -> u32 {
        self.eq_id
    }
}

pub(crate) struct ProbMapEqClassIter<'a> {
    pub current_eq_id: u32,
    pub current_idx: usize,
    pub current_probs_offset: usize,
    prob_map: &'a TempProbMap,
}

impl<'a> Iterator for ProbMapEqClassIter<'a> {
    type Item = EqClassAlnProbView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        // if we haven't yet exhausted the probability map
        if self.current_probs_offset < self.prob_map.probs.len() {
            // the current equivalence class id
            let current_eq_id = self.prob_map.eq_ids[self.current_idx];
            // where the set of probabilities for alignments for all reads in
            // this equivalence class starts
            let start_prob_slice = self.current_probs_offset;
            // the cimulative length array to keep track of the probabilities
            // for the alignments for each individual read of this equivalence
            // class
            let mut cumulative_length = 0_u32;
            let mut cumulative_offsets = vec![cumulative_length];

            // the first iterator that will yield alignments for the first
            // read of this eq class
            let mut it = self
                .prob_map
                .eq_ids
                .iter()
                .skip(self.current_idx)
                .peekable();

            // we should have at least one entry
            let next_eq_id = it
                .next()
                .expect("each eq class should have at least one read");
            // the number of alignments for this read (should be the cardinality
            // of the label of this equivalence class).
            let eq_class_len = self.prob_map.lengths[self.current_idx] as u32;

            // update the cumulative length
            cumulative_length += eq_class_len;
            // and push it on the vector
            cumulative_offsets.push(eq_class_len);

            // update where we point in the global probability array
            self.current_probs_offset += eq_class_len as usize;
            // update our index in the eq_id array
            self.current_idx += 1;

            while let Some(next_eq_id) = it.peek() {
                // if the next read belongs to the same eq class
                if **next_eq_id == current_eq_id {
                    let eq_class_len = self.prob_map.lengths[self.current_idx] as u32;
                    cumulative_length += eq_class_len;
                    cumulative_offsets.push(cumulative_length);
                    self.current_probs_offset += eq_class_len as usize;
                    self.current_idx += 1;
                    it.next();
                } else {
                    break;
                }
            }

            let end_prob_slice = self.current_probs_offset;
            Some(EqClassAlnProbView {
                eq_id: current_eq_id,
                probs: &self.prob_map.probs[start_prob_slice..end_prob_slice],
                cumulative_offsets,
            })
        } else {
            None
        }
    }
}

impl TempProbMap {
    pub fn new() -> Self {
        Self {
            eq_ids: Vec::new(),
            probs: Vec::new(),
            lengths: Vec::new(),
        }
    }

    pub fn push(&mut self, eq_id: u32, probs: &[f64]) {
        self.probs.extend_from_slice(probs);
        self.lengths.push(probs.len());
        self.eq_ids.push(eq_id);
    }

    pub fn len(&self) -> usize {
        self.eq_ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Reorder this probability map's entries by the equivalence
    /// class IDs of the underlying map
    pub fn order_by_eq_id(&mut self) {
        let perm = argsort(&self.eq_ids);

        // get cumulative sums to re-order
        let mut cumsum = Vec::with_capacity(self.lengths.len() + 1);
        let mut tot = 0;
        cumsum.push(0);
        for l in self.lengths.iter() {
            tot += l;
            cumsum.push(tot);
        }

        // reorder the probability slices
        let mut new_probs = Vec::with_capacity(self.probs.len());
        for ind in perm.iter() {
            let start = cumsum[*ind];
            let end = cumsum[*ind + 1];
            new_probs.extend_from_slice(&self.probs[start..end]);
        }
        // reorder the equivalence class IDs
        reorder_in_place(&mut self.eq_ids, &perm);
        // reorder the lengths
        reorder_in_place(&mut self.lengths, &perm);
        std::mem::swap(&mut new_probs, &mut self.probs);
    }

    /// Get an iterator over the alignments for each equivalence class
    pub(crate) fn eq_class_aln_view_iter(&self) -> ProbMapEqClassIter<'_> {
        ProbMapEqClassIter {
            current_eq_id: 0,
            current_idx: 0,
            current_probs_offset: 0,
            prob_map: self,
        }
    }
}

impl EqMap {
    pub fn num_eq_classes(&self) -> usize {
        self.eqc_info.len()
    }

    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.eqc_info.clear();
        // keep nref
        self.label_list_size = 0usize;
        self.eq_labels.clear();
        self.eq_label_starts.clear();
        // clear the label_counts, but resize
        // and fill with 0
        self.label_counts.fill(0u32);

        self.ref_offsets.clear();
        self.ref_labels.clear();
        // keep map_type
        //self.eqid_map.clear();
    }

    pub fn new(nref_in: u32, map_type: EqMapType) -> EqMap {
        let rand_state = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
        EqMap {
            eqc_info: vec![], //HashMap::with_hasher(rs),
            nref: nref_in,
            label_list_size: 0usize,
            eq_labels: vec![],
            eq_label_starts: vec![],
            label_counts: vec![0; nref_in as usize],
            ref_offsets: vec![],
            ref_labels: vec![],
            eqid_map: HashMap::with_hasher(rand_state),
            map_type,
            prob_map: None,
        }
    }

    pub fn fill_ref_offsets(&mut self) {
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

    pub fn fill_label_sizes(&mut self) {
        self.ref_labels = vec![u32::MAX; self.label_list_size + 1];
    }

    #[allow(dead_code)]
    fn init_from_small_chunk<R>(&mut self, cell_chunk: &mut chunk::Chunk<R>)
    where
        R: MappedRecord + UmiTaggedRecord,
    {
        //let rand_state = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
        let mut hasher = self.eqid_map.hasher().build_hasher();
        //let mut hasher = rand_state.build_hasher();

        let mut hash_vec: Vec<(u64, u64, usize)> = cell_chunk
            .reads
            .iter()
            .enumerate()
            .map(|(idx, r)| -> (u64, u64, usize) {
                hasher.write(libradicl::as_u8_slice(r.refs()));
                (hasher.finish(), r.umi(), idx)
            })
            .collect();

        hash_vec.sort_unstable();

        let mut prev_hash = 0u64;
        let mut prev_umi = u64::MAX;
        let mut eq_num = 0usize;
        for (chash, cumi, idx) in hash_vec {
            if let Some(r) = cell_chunk.reads.get(idx) {
                // if the label is the same
                if chash == prev_hash {
                    // and the umi is the same
                    if cumi == prev_umi {
                        // increment the prev umi count
                        self.eqc_info[eq_num].umis.last_mut().unwrap().1 += 1;
                    } else {
                        // if the umi is different
                        self.eqc_info[eq_num].umis.push((cumi, 1));
                    }
                } else {
                    // new class
                    self.label_list_size += r.refs().len();
                    for r in r.refs().iter() {
                        let ridx = *r as usize;
                        self.label_counts[ridx] += 1;
                    }
                    eq_num = self.eq_label_starts.len();
                    self.eq_label_starts.push(self.eq_labels.len() as u32);
                    self.eq_labels.extend(r.refs());
                    self.eqc_info.push(EqMapEntry {
                        umis: vec![(cumi, 1)],
                        eq_num: eq_num as u32,
                    });
                    prev_hash = chash;
                }
                prev_umi = cumi;
            }
        }
        // final value to avoid special cases
        self.eq_label_starts.push(self.eq_labels.len() as u32);
        self.fill_ref_offsets();
        self.fill_label_sizes();

        for idx in 0..self.num_eq_classes() {
            // for each reference in this
            // label, put it in the next free spot
            let label = self.refs_for_eqc(idx as u32).to_vec();
            for r in label {
                self.ref_offsets[r as usize] -= 1;
                self.ref_labels[self.ref_offsets[r as usize] as usize] = self.eqc_info[idx].eq_num;
            }
        }
    }

    pub fn init_from_chunk_gene_level<R>(
        &mut self,
        cell_chunk: &mut chunk::Chunk<R>,
        tid_to_gid: &[u32],
    ) where
        R: MappedRecord + UmiTaggedRecord,
    {
        self.eqid_map.clear();

        let mut gvec: Vec<u32> = vec![];
        // gather the equivalence class info
        for r in &mut cell_chunk.reads {
            // project from txp-level to gene-level
            gvec.extend(r.refs().iter().map(|x| tid_to_gid[*x as usize]));
            gvec.sort_unstable();
            gvec.dedup();

            match self.eqid_map.get_mut(&gvec) {
                // if we've seen this equivalence class before, just add the new
                // umi.
                Some(v) => {
                    self.eqc_info[*v as usize].umis.push((r.umi(), 1));
                }
                // otherwise, add the new umi, but we also have some extra bookkeeping
                None => {
                    // each reference in this equivalence class label
                    // will have to point to this equivalence class id
                    let eq_num = self.eqc_info.len() as u32;
                    self.label_list_size += gvec.len();
                    for r in gvec.iter() {
                        let ridx = *r as usize;
                        self.label_counts[ridx] += 1;
                        //ref_to_eqid[*r as usize].push(eq_num);
                    }
                    self.eq_label_starts.push(self.eq_labels.len() as u32);
                    self.eq_labels.extend(&gvec);
                    self.eqc_info.push(EqMapEntry {
                        umis: vec![(r.umi(), 1)],
                        eq_num,
                    });
                    self.eqid_map.insert(gvec.clone(), eq_num);
                }
            }

            gvec.clear();
        }
        // final value to avoid special cases
        self.eq_label_starts.push(self.eq_labels.len() as u32);

        self.fill_ref_offsets();
        self.fill_label_sizes();

        // initially we inserted duplicate UMIs
        // here, collapse them and keep track of their count
        for idx in 0..self.num_eq_classes() {
            // for each reference in this
            // label, put it in the next free spot
            // TODO: @k3yavi, can we avoid this copy?
            let label = self.refs_for_eqc(idx as u32).to_vec();
            //println!("{:?}", label);
            for r in label {
                self.ref_offsets[r as usize] -= 1;
                self.ref_labels[self.ref_offsets[r as usize] as usize] = self.eqc_info[idx].eq_num;
            }

            let v = &mut self.eqc_info[idx];
            // sort so dups are adjacent
            v.umis.sort_unstable();
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

    pub(crate) fn init_from_chunk<R>(&mut self, cell_chunk: &mut chunk::Chunk<R>)
    where
        R: MappedRecord + UmiTaggedRecord + OptionalAlignmentScores,
    {
        /*
        if cell_chunk.reads.len() < 10 {
        self.init_from_small_chunk(cell_chunk);
        return;
        }
        */
        // temporary map of equivalence class label to assigned
        // index.
        // let mut eqid_map:
        // let rand_state = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
        // let mut eqid_map  = HashMap::with_hasher(rand_state);

        // a temporary probability map as we will need to rearrange things
        let mut temp_prob_map = TempProbMap::new();
        // the final map that will be added to the EQMap
        let mut prob_map = ProbMap::new();

        self.eqid_map.clear();
        // gather the equivalence class info
        for r in &mut cell_chunk.reads {
            // Take what we need from r up-front
            let refs = r.refs();
            let _umi = r.umi();

            // if the underlying record type provides scores, then we
            // get some.
            let maybe_scores = r.maybe_scores();

            match self.eqid_map.get_mut(refs) {
                // if we've seen this equivalence class before, just add the new
                // umi.
                Some(v) => {
                    self.eqc_info[*v as usize].umis.push((r.umi(), 1));
                    // if we have scores, add them labeled with this equivalence class
                    if let Some(scores) = maybe_scores {
                        let score_probs = score_probabilities(scores);
                        // push score probs with associated eq_id of *v
                        temp_prob_map.push(*v, &score_probs);
                    };
                }
                // otherwise, add the new umi, but we also have some extra bookkeeping
                None => {
                    // each reference in this equivalence class label
                    // will have to point to this equivalence class id
                    let eq_num = self.eqc_info.len() as u32;
                    self.label_list_size += r.refs().len();
                    for r in r.refs().iter() {
                        let ridx = *r as usize;
                        self.label_counts[ridx] += 1;
                        //ref_to_eqid[*r as usize].push(eq_num);
                    }
                    self.eq_label_starts.push(self.eq_labels.len() as u32);
                    self.eq_labels.extend(r.refs());
                    self.eqc_info.push(EqMapEntry {
                        umis: vec![(r.umi(), 1)],
                        eq_num,
                    });
                    self.eqid_map.insert(r.refs().to_vec(), eq_num);
                    // if we have scores, add them labeled with this equivalence class
                    if let Some(scores) = maybe_scores {
                        let score_probs = score_probabilities(scores);
                        // push score probs with associated eq_id of *v
                        temp_prob_map.push(eq_num, &score_probs);
                    }
                    //self.eqc_map.insert(r.refs.clone(), EqMapEntry { umis : vec![(r.umi,1)], eq_num });
                }
            }
        }
        // final value to avoid special cases
        self.eq_label_starts.push(self.eq_labels.len() as u32);

        self.fill_ref_offsets();
        self.fill_label_sizes();

        // sort the probability map by equivalence class ids
        // if we have them
        let have_probs = !temp_prob_map.is_empty();
        if have_probs {
            temp_prob_map.order_by_eq_id();
            //prob_map.umi_offsets.push(0);
            //prob_map.eq_indices.push(0);
        }
        // get an iterartor to iterate over alignment probabilities to the
        // reads belonging to each equivalence class (in order)
        let mut aln_view_iter = temp_prob_map.eq_class_aln_view_iter();

        // initially we inserted duplicate UMIs
        // here, collapse them and keep track of their count
        // iterate over the equivalence classes in order.
        for idx in 0..self.num_eq_classes() {
            // for each reference in this
            // label, put it in the next free spot
            // TODO: @k3yavi, can we avoid this copy?
            let label = self.refs_for_eqc(idx as u32).to_vec();
            //println!("{:?}", label);
            for r in label {
                self.ref_offsets[r as usize] -= 1;
                self.ref_labels[self.ref_offsets[r as usize] as usize] = self.eqc_info[idx].eq_num;
            }

            let v = &mut self.eqc_info[idx];

            // get the permutation that puts the UMIs in order by their value
            let perm = argsort_by(&v.umis, |a, b| a.0.cmp(&b.0));

            // if we have probabilites, reorder them as we have the UMIs
            if have_probs {
                // get the alignment view for the next equivalence class
                // (i.e. alignments for all reads with this mapping pattern)
                if let Some(aln_view) = aln_view_iter.next() {
                    // the first UMI in this equivalence class
                    let mut curr_umi = v.umis[perm[0]].0;

                    // for each read, in the appropriate order
                    for p in perm.iter() {
                        // if the current UMI changed, then update the current
                        // umi and mark then boundary offset for the last probability
                        // belonging to the previous umi.
                        if v.umis[*p].0 != curr_umi {
                            curr_umi = v.umis[*p].0;
                            prob_map.mark_umi_end();
                        }
                        // get the probability vector for this read
                        let probs = aln_view
                            .get_probs_for_read_rank(*p)
                            .expect("missing read rank");
                        // add them to the probability vector we are building
                        // in the final probability map
                        prob_map.probs.extend_from_slice(probs);
                    }
                    // make sure we mark the end offset for the last UMI
                    prob_map.mark_umi_end();

                    // once we've visited all reads for this equivalence class
                    // mark the end position in the eq_indices vector
                    prob_map.mark_eq_class_end();
                } else {
                    eprintln!(
                        "Number of probability vectors should equal total number of equivalence classes"
                    );
                }
            }

            // sort so dups are adjacent
            reorder_in_place(&mut v.umis, &perm);
            // v.umis.sort_unstable();

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
        // if we had read probabilities set them
        if !temp_prob_map.is_empty() {
            self.prob_map = Some(prob_map);
        } else {
            self.prob_map = None;
        }
    }

    /// Returns the probability row (length = label_len) for a specific graph node (eqid, umi_idx).
    pub fn probs_for_eq_umi_tx(
        &self,
        eqid: u32,
        umi_idx: u32,
        tx_index: usize,
    ) -> Option<Vec<f64>> {
        if let Some(ref pmap) = self.prob_map {
            let label_len = self.refs_for_eqc(eqid).len();
            // the index in umi_offsets where eqid starts
            let pmap_eq_idx = pmap.eq_indices[eqid as usize];
            // the umi_idx umi for this eq class
            let poffset = pmap.umi_offsets[pmap_eq_idx + umi_idx as usize];
            // number of occurrences of this UMI
            let freq = self.eqc_info[eqid as usize].umis[umi_idx as usize].1 as usize;
            let probs: Vec<f64> = pmap
                .probs
                .iter()
                .skip(poffset + tx_index)
                .step_by(label_len)
                .take(freq)
                .copied()
                .collect();
            Some(probs)
        } else {
            None
        }
    }

    pub fn eq_classes_containing(&self, r: u32) -> &[u32] {
        &self.ref_labels
            [(self.ref_offsets[r as usize] as usize)..(self.ref_offsets[(r + 1) as usize] as usize)]
    }

    pub fn refs_for_eqc(&self, idx: u32) -> &[u32] {
        &self.eq_labels[(self.eq_label_starts[idx as usize] as usize)
            ..(self.eq_label_starts[(idx + 1) as usize] as usize)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn temp_prob_map_add() {
        let mut tpm = TempProbMap::new();
        tpm.push(0, &[0.5, 0.5]);
        tpm.push(1, &[0.25, 0.25, 0.5]);
        tpm.push(0, &[0.75, 0.2]);
        assert_eq!(tpm.probs.len(), 7);
        assert_eq!(tpm.lengths, vec![2, 3, 2]);
        assert_eq!(tpm.eq_ids, vec![0_u32, 1, 0]);

        tpm.order_by_eq_id();
        assert_eq!(tpm.probs.len(), 7);
        assert_eq!(tpm.lengths, vec![2, 2, 3]);
        assert_eq!(tpm.eq_ids, vec![0_u32, 0, 1]);
    }

    #[test]
    fn temp_prob_map_iter() {
        let mut tpm = TempProbMap::new();
        tpm.push(0, &[0.5, 0.5]);
        tpm.push(1, &[0.25, 0.25, 0.5]);
        tpm.push(0, &[0.75, 0.25]);
        assert_eq!(tpm.probs.len(), 7);
        assert_eq!(tpm.lengths, vec![2, 3, 2]);
        assert_eq!(tpm.eq_ids, vec![0_u32, 1, 0]);

        tpm.order_by_eq_id();
        assert_eq!(tpm.probs.len(), 7);
        assert_eq!(tpm.lengths, vec![2, 2, 3]);
        assert_eq!(tpm.eq_ids, vec![0_u32, 0, 1]);

        let mut tpm_iter = tpm.eq_class_aln_view_iter();
        let eq_cards = vec![2, 3];
        let mut ctr = 0;
        while let Some(ti) = tpm_iter.next() {
            for rank in 0..ti.num_reads() {
                let p = ti.get_probs_for_read_rank(rank).expect("have read of rank");
                if ctr == 0 {
                    assert_eq!(ti.num_reads(), 2);
                    assert_eq!(&[0.5, 0.5], p);
                } else if ctr == 1 {
                    assert_eq!(&[0.75, 0.25], p);
                } else if ctr == 2 {
                    assert_eq!(ti.num_reads(), 1);
                    assert_eq!(&[0.25, 0.25, 0.5], p);
                }
                ctr += 1;
            }
        }
    }

    #[test]
    fn prob_map_access() {
        let mut pm = ProbMap::new();

        // eq 0 has 2 UMIs, one of frequency 3 the other of frequency 1
        pm.probs.extend_from_slice(&[0.1, 0.7, 0.1, 0.1]);
        pm.probs.extend_from_slice(&[0.05, 0.7, 0.15, 0.1]);
        pm.probs.extend_from_slice(&[0.15, 0.7, 0.05, 0.1]);
        pm.mark_umi_end();
        pm.probs.extend_from_slice(&[0.2, 0.2, 0.2, 0.4]);
        pm.mark_umi_end();
        pm.mark_eq_class_end();

        // eq 1 has 1 UMI of frequency 2
        pm.probs.extend_from_slice(&[0.25, 0.75]);
        pm.probs.extend_from_slice(&[0.2, 0.8]);
        pm.mark_umi_end();
        pm.mark_eq_class_end();

        // eq 2 has 3 UMIs each of frequency 1
        pm.probs.extend_from_slice(&[0.2, 0.6, 0.2]);
        pm.mark_umi_end();
        pm.probs.extend_from_slice(&[0.15, 0.35, 0.5]);
        pm.mark_umi_end();
        pm.probs.extend_from_slice(&[0.5, 0.2, 0.3]);
        pm.mark_umi_end();
        pm.mark_eq_class_end();

        assert_eq!(pm.num_umis_for_eq(0), Some(2));
        assert_eq!(pm.num_umis_for_eq(1), Some(1));
        assert_eq!(pm.num_umis_for_eq(2), Some(3));
        assert_eq!(pm.num_umis_for_eq(3), None);
    }
}
