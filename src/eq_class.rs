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
use libradicl::record::AlevinFryReadRecord;

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
    fn init_from_small_chunk(&mut self, cell_chunk: &mut chunk::Chunk<AlevinFryReadRecord>) {
        //let rand_state = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
        let mut hasher = self.eqid_map.hasher().build_hasher();
        //let mut hasher = rand_state.build_hasher();

        let mut hash_vec: Vec<(u64, u64, usize)> = cell_chunk
            .reads
            .iter()
            .enumerate()
            .map(|(idx, r)| -> (u64, u64, usize) {
                hasher.write(libradicl::as_u8_slice(&r.refs[..]));
                (hasher.finish(), r.umi, idx)
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
                    self.label_list_size += r.refs.len();
                    for r in r.refs.iter() {
                        let ridx = *r as usize;
                        self.label_counts[ridx] += 1;
                    }
                    eq_num = self.eq_label_starts.len();
                    self.eq_label_starts.push(self.eq_labels.len() as u32);
                    self.eq_labels.extend(&r.refs);
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

    pub fn init_from_chunk_gene_level(
        &mut self,
        cell_chunk: &mut chunk::Chunk<AlevinFryReadRecord>,
        tid_to_gid: &[u32],
    ) {
        self.eqid_map.clear();

        let mut gvec: Vec<u32> = vec![];
        // gather the equivalence class info
        for r in &mut cell_chunk.reads {
            // project from txp-level to gene-level
            gvec.extend(r.refs.iter().map(|x| tid_to_gid[*x as usize]));
            gvec.sort_unstable();
            gvec.dedup();

            match self.eqid_map.get_mut(&gvec) {
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
                    self.label_list_size += gvec.len();
                    for r in gvec.iter() {
                        let ridx = *r as usize;
                        self.label_counts[ridx] += 1;
                        //ref_to_eqid[*r as usize].push(eq_num);
                    }
                    self.eq_label_starts.push(self.eq_labels.len() as u32);
                    self.eq_labels.extend(&gvec);
                    self.eqc_info.push(EqMapEntry {
                        umis: vec![(r.umi, 1)],
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

    pub fn init_from_chunk(&mut self, cell_chunk: &mut chunk::Chunk<AlevinFryReadRecord>) {
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

        self.eqid_map.clear();
        // gather the equivalence class info
        for r in &mut cell_chunk.reads {
            // TODO: ensure this is done upstream so we
            // don't have to do it here.
            // NOTE: should be done if collate was run.
            // r.refs.sort();

            match self.eqid_map.get_mut(&r.refs) {
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
                    self.eqid_map.insert(r.refs.clone(), eq_num);
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

    pub fn eq_classes_containing(&self, r: u32) -> &[u32] {
        &self.ref_labels
            [(self.ref_offsets[r as usize] as usize)..(self.ref_offsets[(r + 1) as usize] as usize)]
    }

    pub fn refs_for_eqc(&self, idx: u32) -> &[u32] {
        &self.eq_labels[(self.eq_label_starts[idx as usize] as usize)
            ..(self.eq_label_starts[(idx + 1) as usize] as usize)]
    }
}
