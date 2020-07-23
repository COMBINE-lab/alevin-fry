use std::collections::{HashMap, HashSet};
use std::error::Error;

/// FROM https://github.com/10XGenomics/rust-debruijn/blob/master/src/dna_string.rs
/// count Hamming distance between 2 2-bit DNA packed u64s
pub(super) fn count_diff_2_bit_packed(a: u64, b: u64) -> usize {
    let bit_diffs = a ^ b;
    let two_bit_diffs = (bit_diffs | bit_diffs >> 1) & 0x5555555555555555;
    let total_diffs = two_bit_diffs.count_ones() as usize;
    return total_diffs;
}

fn get_bit_mask(nt_index: usize, fill_with: u64) -> Result<u64, Box<dyn Error>> {
    let mut mask: u64 = fill_with;
    mask <<= 2 * (nt_index - 1);
    Ok(mask)
}

fn get_all_snps(bc: u64, bc_length: usize) -> Result<Vec<u64>, Box<dyn Error>> {
    assert!(
        bc <= 2u64.pow(2 * bc_length as u32),
        "the barcode id is larger than possible (based on barcode length)"
    );
    assert!(
        bc_length <= 32,
        "barcode length greater than 32 not supported"
    );

    let mut snps: Vec<u64> = Vec::new();
    snps.reserve(3 * bc_length);

    for nt_index in 1..=bc_length {
        // clearing the two relevant bits based on nucleotide position
        let bit_mask = bc & !get_bit_mask(nt_index, 3)?;

        // iterating over all 4 choices of the nucleotide
        for i in 0..=3 {
            let new_bc = bit_mask | get_bit_mask(nt_index, i)?;
            if new_bc != bc {
                snps.push(new_bc);
            }
        }
    }

    Ok(snps)
}

fn get_all_indels(bc: u64, bc_length: usize) -> Result<Vec<u64>, Box<dyn Error>> {
    assert!(
        bc <= 2u64.pow(2 * bc_length as u32),
        "the barcode id is larger than possible (based on barcode length)"
    );
    assert!(
        bc_length <= 32,
        "barcode length greater than 32 not supported"
    );

    let mut indels: Vec<u64> = Vec::new();
    indels.reserve(8 * (bc_length - 1));

    for nt_index in 1..bc_length {
        let mut bit_mask = 1 << (2 * nt_index);
        bit_mask -= 1;

        let upper_half = bc & !bit_mask;
        let lower_half = bc & bit_mask;

        // iterating over all 4 choices of the nucleotide
        for i in 0..=3 {
            let new_insertion_bc = upper_half | get_bit_mask(nt_index, i)? | (lower_half >> 2);
            let new_deletion_bc = upper_half
                | get_bit_mask(1, i)?
                | ((lower_half & !get_bit_mask(nt_index + 1, 3)?) << 2);

            if new_insertion_bc != bc {
                indels.push(new_insertion_bc);
            }
            if new_deletion_bc != bc {
                indels.push(new_deletion_bc);
            }
        }
    }

    Ok(indels)
}

fn get_all_one_edit_neighbors(
    bc: u64,
    bc_length: usize,
    neighbors: &mut HashSet<u64>,
) -> Result<(), Box<dyn Error>> {
    neighbors.clear();

    let snps: Vec<u64> = get_all_snps(bc, bc_length)?;
    let indels: Vec<u64> = get_all_indels(bc, bc_length)?;

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

#[cfg(test)]
mod tests {
    use utils::*;

    #[test]
    fn test_get_bit_mask() {
        let mut output = Vec::new();
        for i in 0..=3 {
            let mask = get_bit_mask(2, i).unwrap();
            output.push(mask);
        }
        assert_eq!(output, vec![0, 4, 8, 12]);
    }

    #[test]
    fn test_get_all_snps() {
        let mut output: Vec<u64> = get_all_snps(7, 3).unwrap().into_iter().collect();
        output.sort();

        assert_eq!(output, vec![3, 4, 5, 6, 11, 15, 23, 39, 55]);
    }

    #[test]
    fn test_get_all_indels() {
        let mut output: Vec<u64> = get_all_indels(7, 3).unwrap().into_iter().collect();
        output.sort();
        output.dedup();

        assert_eq!(output, vec![1, 4, 5, 6, 9, 12, 13, 14, 15, 28, 29, 30, 31]);
    }

    #[test]
    fn test_get_all_one_edit_neighbors() {
        let mut neighbors: HashSet<u64> = HashSet::new();
        get_all_one_edit_neighbors(7, 3, &mut neighbors).unwrap();

        let mut output: Vec<u64> = neighbors.into_iter().collect();

        output.sort();
        output.dedup();

        assert_eq!(
            output,
            vec![1, 3, 4, 5, 6, 9, 11, 12, 13, 14, 15, 23, 28, 29, 30, 31, 39, 55]
        );
    }

    #[test]
    fn test_generate_whitelist_hash() {
        let neighbors: HashSet<u64> = generate_whitelist_hash(&vec![7], 3).unwrap();
        let mut output: Vec<u64> = neighbors.into_iter().collect();

        output.sort();
        output.dedup();

        assert_eq!(
            output,
            vec![1, 3, 4, 5, 6, 9, 11, 12, 13, 14, 15, 23, 28, 29, 30, 31, 39, 55]
        );
    }
}
