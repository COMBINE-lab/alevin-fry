use std::collections::HashSet;
use std::error::Error;

fn get_bit_mask(nt_index: usize, fill_with: u64) -> Result<u64, Box<dyn Error>> {
    let mut mask: u64 = fill_with;
    mask = mask << (2 * (nt_index - 1));

    Ok(mask)
}

fn get_all_snps(bc: u64, bc_length: usize) -> Result<Vec<u64>, Box<dyn Error>> {
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
            let new_deletion_bc =
                upper_half | get_bit_mask(1, i)? | (lower_half & !get_bit_mask(nt_index, 3)?);

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

pub fn generate_whitelist_hash(
    whitelist_bcs: &Vec<u64>,
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
        get_all_one_edit_neighbors(bc.clone(), bc_length, &mut neighbors)?;
        one_edit_barcode_hash.extend(&neighbors);
    }

    Ok(one_edit_barcode_hash)
}
