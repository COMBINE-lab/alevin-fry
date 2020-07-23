use std::collections::HashMap;
use fasthash::{sea, RandomState};

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
    pub umis : Vec<(u64, u32)>,
    pub eq_num : u32
}

pub(super) struct EqMap {
    pub eqc_map : HashMap<Vec<u32>, EqMapEntry, fasthash::RandomState<fasthash::sea::Hash64>>,
    pub nref : u32,
    pub label_list_size : usize,
    // concatenated lists of the labels of all equivalence classes
    pub eq_labels : Vec<u32>, //= Vec::new();
    // vector that deliniates where each equivalence class label 
    // begins and ends.  The label for equivalence class i begins 
    // at offset eq_label_starts[i], and it ends at 
    // eq_label_starts[i+1].  The length of this vector is 1 greater
    // than the number of equivalence classes.
    pub eq_label_starts : Vec<u32>, //= Vec::new();
    pub label_counts : Vec<u32>,// = vec![0; nref as usize];
    pub ref_offsets : Vec<u32>,
    pub ref_labels : Vec<u32> //= vec![ u32::MAX; label_list_size];
}

impl EqMap {

    pub(super) fn clear(&mut self) {
        self.eqc_map.clear();
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

    pub(super) fn new(rs : RandomState::<sea::Hash64>, nref_in : u32) -> EqMap {
        EqMap {
            eqc_map : HashMap::with_hasher(rs),
            nref : nref_in,
            label_list_size : 0usize,
            eq_labels : vec![],
            eq_label_starts : vec![],
            label_counts : vec![0; nref_in as usize],
            ref_offsets : vec![],
            ref_labels : vec![]
        }
    }

    pub(super) fn fill_ref_offsets(&mut self) {
        self.ref_offsets = self.label_counts.iter().scan(0, |sum, i| {*sum += i; Some(*sum)}).collect::<Vec<_>>();
        self.ref_offsets.push(*self.ref_offsets.last().unwrap());
    }

    pub(super) fn fill_label_sizes(&mut self)  {
        self.ref_labels = vec![ u32::MAX; self.label_list_size];
    }

    //pub(super) 
}


#[derive(Debug)]
pub(super) enum PUGEdgeType {
   NoEdge,
   BiDirected,
   XToY,
   YToX,
}