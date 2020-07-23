extern crate fasthash;

use fasthash::{sea, RandomState};
use std::collections::HashMap;

use crate as libradicl;

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

    pub(super) fn new(rs: RandomState<sea::Hash64>, nref_in: u32) -> EqMap {
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

    pub(super) fn init_from_chunk(&mut self, cell_chunk: &libradicl::Chunk) {

        // temporary map of equivalence class label to assigned
        // index.
        let s = RandomState::<sea::Hash64>::new();
        let mut eqid_map: HashMap<Vec<u32>, u32, fasthash::RandomState<fasthash::sea::Hash64>> =
            HashMap::with_hasher(s);

        // gather the equivalence class info
        for r in &cell_chunk.reads {
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
            v.umis.sort();
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

#[derive(Debug)]
pub(super) enum PUGEdgeType {
    NoEdge,
    BiDirected,
    XToY,
    YToX,
}
