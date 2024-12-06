use needletail::bitkmer::*;

pub fn get_bin_id(pos: u32, ref_id: usize, size_range: u32, blens: &[u64]) -> usize {
    let bid = pos / size_range;
    let ind: usize = (blens[ref_id] + bid as u64) as usize;
    ind
}

pub fn get_bc_string(kmerseq: &BitKmerSeq, reverse_barcode: bool, bc_len: u8) -> String {
    let kmseq = *kmerseq;
    let mut km: BitKmer = (kmseq, bc_len);
    if reverse_barcode {
        km = needletail::bitkmer::reverse_complement(km);
    }
    let bytes = needletail::bitkmer::bitmer_to_bytes(km);
    let seq: String = String::from_utf8(bytes).expect("Invalid barcode");
    seq
}
