use crate::mlp_spline::predict_with_tch;
use anyhow::Result;
use ndarray::prelude::*;
use ndarray::{Array1, Array2};
use ndarray::s;
use serde_json::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use tch::nn;
type MCCIndex = usize;
type CoveringTxpId = u32;
/// RAD alignment tuple: (is_reverse, ref_start)
type AlgnTuple = (bool, u32);
type ForsetiCheckingList = HashMap<(MCCIndex, CoveringTxpId), Vec<AlgnTuple>>;
use memchr::memmem;

fn reverse_complement(seq: &str) -> Result<String, String> {
    let complement = |base: char| match base {
        'A' | 'a' => Ok('T'),
        'T' | 't' => Ok('A'),
        'C' | 'c' => Ok('G'),
        'G' | 'g' => Ok('C'),
        'N' | 'n' => Ok('N'),
        _ => Err(format!("Invalid nucleotide found: {}", base)),
    };

    seq.chars().rev().map(complement).collect() // Collect into a Result<String, _>
}

fn build_kmers(sequence: &str, ksize: usize) -> Vec<String> {
    // If the sequence is shorter than k, there are no valid k-mers.
    // IMPORTANT: using saturating_sub here is *not* sufficient because it would yield 1
    // and then attempt to slice `sequence[0..ksize]`, which panics.
    if sequence.len() < ksize {
        return Vec::new();
    }
    let n_kmers = sequence.len() - ksize + 1;
    (0..n_kmers)
        .map(|i| sequence[i..i + ksize].to_string())
        .collect()
}

// fn __work_one_hot_encoder(kmer_list: &[String], has_enough_a: &[bool]) -> Array2<f32> {
//     let nucleotides = ['A', 'C', 'G', 'T', 'N'];
//     let num_classes = nucleotides.len();

//     // Create a mapping from ASCII codes to indices
//     let mut code_to_idx = [-1i32; 256];
//     for (i, &nuc) in nucleotides.iter().enumerate() {
//         code_to_idx[nuc as usize] = i as i32;
//     }

//     // Filter valid kmers
//     let valid_kmers: Vec<&String> = kmer_list
//         .iter()
//         .zip(has_enough_a)
//         .filter_map(|(kmer, &has_a)| if has_a { Some(kmer) } else { None })
//         .collect();

//     if valid_kmers.is_empty() {
//         return Array2::<f32>::zeros((0, 0));
//     }

//     let num_sequences = valid_kmers.len();
//     let sequence_length = valid_kmers[0].len();

//     // Prepare the output array
//     let mut one_hot = Array2::<f32>::zeros((num_sequences, sequence_length * num_classes));

//     for (i, seq) in valid_kmers.iter().enumerate() {
//         let seq_bytes = seq.as_bytes();
//         for (j, &byte) in seq_bytes.iter().enumerate() {
//             let idx = code_to_idx[byte as usize];
//             let idx = idx as usize;
//             one_hot[(i, j * num_classes + idx)] = 1.0;
//         }
//     }

//     one_hot
// }
fn one_hot_encoder(kmer_list: &[String], has_enough_a: &Array1<bool>) -> Array2<f32> {
    let nucleotides = ['A', 'C', 'G', 'T', 'N'];
    let num_classes = nucleotides.len();

    // Create a mapping from ASCII codes to indices
    let mut code_to_idx = [-1i32; 256];
    for (i, &nuc) in nucleotides.iter().enumerate() {
        code_to_idx[nuc as usize] = i as i32;
    }

    // Filter valid kmers using ndarray boolean masking
    let valid_kmers: Vec<&String> = kmer_list
        .iter()
        .zip(has_enough_a.iter())
        .filter_map(|(kmer, &has_a)| if has_a { Some(kmer) } else { None })
        .collect();

    if valid_kmers.is_empty() {
        return Array2::<f32>::zeros((0, 0));
    }

    let num_sequences = valid_kmers.len();
    let sequence_length = valid_kmers[0].len();

    // Prepare the output array
    let mut one_hot = Array2::<f32>::zeros((num_sequences, sequence_length * num_classes));

    // Process each kmer
    for (i, seq) in valid_kmers.iter().enumerate() {
        let seq_bytes = seq.as_bytes();
        let indices: Vec<usize> = seq_bytes
            .iter()
            .map(|&byte| code_to_idx[byte as usize] as usize)
            .collect();

        // Vectorized one-hot encoding for the current sequence
        for (j, &idx) in indices.iter().enumerate() {
            one_hot[(i, j * num_classes + idx)] = 1.0;
        }
    }

    one_hot
}

fn load_spline_lookup_table(file_path: &str) -> Result<Array1<f64>> {
    // Open the file
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    // Parse the JSON
    let json_data: Value = serde_json::from_reader(reader)?;

    // Extract the "y" field as an array
    if let Some(y_values) = json_data["y"].as_array() {
        // Convert JSON array to Vec<f64>
        let y_vec: Vec<f64> = y_values
            .iter()
            .map(|v| v.as_f64().unwrap_or(0.0)) // Ensure conversion
            .collect();

        // Convert Vec<f64> to ndarray::Array1
        Ok(Array1::from(y_vec))
    } else {
        Err(anyhow::anyhow!("Missing 'y' field in JSON"))
    }
}
fn process_binding_affinity(
    all_30mers: &[String],
    has_enough_a: Array1<bool>,
    mlp: &nn::Sequential,
    discount_perc: f64,
    binding_affinity_threshold: f64,
) -> Result<Array1<f64>> {
    // Collect indices where `has_enough_a` is true
    let indices: Vec<usize> = has_enough_a
        .indexed_iter()
        .filter_map(|(i, &val)| if val { Some(i) } else { None })
        .collect();
    // Filter `all_30mers` using the collected indices
    let filtered_kmers: Vec<&String> = indices.iter().map(|&i| &all_30mers[i]).collect();

    let encoded = one_hot_encoder(all_30mers, &has_enough_a);

    let mut nonzero_binding_affinity = predict_with_tch(&mlp, encoded)?;
    // Ensure dimensions match between `nonzero_binding_affinity` and filtered `has_enough_a`
    if nonzero_binding_affinity.len() != has_enough_a.iter().filter(|&&val| val).count() {
        return Err(anyhow::anyhow!(
            "Mismatched dimensions between binding affinity and has_enough_a"
        ));
    }

    // Apply discount percentage
    nonzero_binding_affinity *= discount_perc;

    // Handle "all A" patterns directly using indices
    let pattern = "A".repeat(30);
    for (i, kmer) in filtered_kmers.iter().enumerate() {
        if kmer == &&pattern {
            nonzero_binding_affinity[i] = 1.0;
        } else if nonzero_binding_affinity[i] <= binding_affinity_threshold {
            nonzero_binding_affinity[i] = 0.0;
        }
    }

    // Get the length of downstream 30-mers
    let mut downstream_binding_affinity = Array1::<f64>::zeros(has_enough_a.len());

    // Assign the affinity scores to the corresponding indices
    for (idx, &affinity) in indices.iter().zip(nonzero_binding_affinity.iter()) {
        downstream_binding_affinity[*idx] = affinity;
    }

    Ok(downstream_binding_affinity)
}
pub fn forseti_for_multi_best(
    forseti_checking_list: &ForsetiCheckingList,
    ref_names: &[String],
    spliceu_txome: &HashMap<u32, Vec<u8>>,
    spline_lookup: &Array1<f64>,
    mlp: &nn::Sequential,
    read_length: u16,
) -> Result<Vec<u16>> {
    // Set up parameters
    let snr_min_size = 6;
    let discount_perc = 1.0_f64;
    let polya_tail_len = 200;
    let max_frag_len = 1000;
    let binding_affinity_threshold = 0.0;
    let pattern = "A".repeat(snr_min_size);
    let pattern_bytes = pattern.as_bytes();

    // Avoid turning stderr into the bottleneck when many MCCs are skipped.
    const INVALID_POS_PRINT_LIMIT: u64 = 2;
    let mut invalid_pos_skipped: u64 = 0;

    // Reusable buffers to reduce per-MCC allocations.
    let mut ref_start_list: Vec<usize> = Vec::new();
    let mut ref_end_list: Vec<usize> = Vec::new();
    // reusable sum buffers for accumulating joint probabilities across alignments
    let mut sum_log_probs: Vec<f64> = Vec::new();
    let mut tail_sum_log_probs: Vec<f64> = Vec::new();
    const EPS: f64 = 1e-12;

    let mut mcc_to_tx_prob_dict: HashMap<MCCIndex, f64> = HashMap::new();
    // algn_tuple_list is direction(fw, reverse) and ref_start
    for ((mcc_idx, covering_txp_id), algn_tuple_list) in forseti_checking_list {
        let mut norm_sum_joint_prob = f64::NEG_INFINITY;

        let ref_seq_bytes = match spliceu_txome.get(covering_txp_id) {
            Some(seq) => seq.as_slice(),
            None => {
                let nm = ref_names
                    .get(*covering_txp_id as usize)
                    .map(|s| s.as_str())
                    .unwrap_or("<unknown>");
                eprintln!("Error: ref_id {} (name {}) not found in spliceu_txome.", covering_txp_id, nm);
                continue;
            }
        };
        let ref_seq = match std::str::from_utf8(ref_seq_bytes) {
            Ok(s) => s,
            Err(_) => {
                eprintln!("Error: ref_id {} sequence bytes were not valid utf-8.", covering_txp_id);
                continue;
            }
        };
        let tx_ref_end = ref_seq.len();
        let all_forward = algn_tuple_list.iter().all(|algn_tuple| algn_tuple.0);
        let all_reverse = algn_tuple_list.iter().all(|algn_tuple| !algn_tuple.0);

        let invalid_pos_limit: u32 = u32::MAX - read_length as u32 - 50;
        let mut mcc_invalid = false;

        if all_forward {
            ref_start_list.clear();
            ref_start_list.reserve(algn_tuple_list.len());
            for &(_is_forward, ref_start_u32) in algn_tuple_list.iter() {
                // NOTE: Mimic softclip
                // TODO:we could have a better way to handle this. if we can use signed int for the ref_start, we could know the length of overhang/clipping.
                if ref_start_u32 > invalid_pos_limit {
                    ref_start_list.push(0 as usize);
                }else{
                    let rs = ref_start_u32 as usize;
                    if rs > tx_ref_end {
                        println!("ref_start: {}", rs);
                        eprintln!("Error: ref_start is not close to 2^32, but exceed ref length");
                        mcc_invalid = true;
                        break;
                    }
                    ref_start_list.push(rs);
                }
            }
            if mcc_invalid{
                continue;
            }
            if ref_start_list.is_empty() {
                eprintln!("ref_start_list is empty! Should not happen.");
                continue;
            }

            let overlap_wdow_start = *ref_start_list.iter().max().unwrap();
            let overlap_wdow_end = *ref_start_list.iter().min().unwrap() + max_frag_len;
            if overlap_wdow_start >= overlap_wdow_end {
                continue;
            }

            // Downstream processing
            // if we have >30 bases downstream, we want to consider internal polyA sites
            // we use > 30 because the polyA should start one base after the window range start
            if tx_ref_end - overlap_wdow_start > 30 {
                let downstream_seq = &ref_seq
                    [(overlap_wdow_start + 1)..usize::min(overlap_wdow_end + 30, tx_ref_end)];
                let downstream_30mers = build_kmers(downstream_seq, 30);
                if downstream_30mers.is_empty() {
                    continue;
                }
                let has_enough_a: Array1<bool> = Array1::from(
                    downstream_30mers
                        .iter()
                        .map(|kmer| memmem::find(kmer.as_bytes(), pattern_bytes).is_some())
                        .collect::<Vec<bool>>(),
                );
                // if any of the 30mers has enough A, we can process the binding affinity for this mcc
                if has_enough_a.iter().any(|&x| x) {
                    let downstream_binding_affinity = process_binding_affinity(
                        &downstream_30mers,
                        has_enough_a,
                        mlp,
                        discount_perc,
                        binding_affinity_threshold,
                    )?;

                    sum_log_probs.clear();
                    sum_log_probs.resize(downstream_30mers.len(), 0.0);

                    // for each alignment, we compute the prob. distanse = (each algn's start to the overlapped window end), while poly A range is overlap wdow start to end.
                    for &ref_start in &ref_start_list {
                        let prefix_dis = overlap_wdow_start - ref_start;
                        let start_idx = prefix_dis + 1;
                        let end_idx = prefix_dis + downstream_30mers.len() + 1;
                        let downstream_frag_len_prob = spline_lookup.slice(s![start_idx..end_idx]);
                        for (i, (&affinity, &frag_prob)) in downstream_binding_affinity
                            .iter()
                            .zip(downstream_frag_len_prob.iter())
                            .enumerate()
                        {
                            // avoid log(0), + EPS
                            sum_log_probs[i] += (affinity * frag_prob + EPS).ln();
                        }
                    }

                    let max_sum_log = sum_log_probs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    // make it comparable across different n_alns
                    norm_sum_joint_prob = max_sum_log / ref_start_list.len() as f64;
                }
            }
            // the downstream window end is beyond the last 30 mer of the ref, we could borrow A from the poly A tail.
            // And, we only compute the case that we consider borrowed A from tail(avoid duplicate computation for cases above, pure internal polyA)
            if overlap_wdow_end > tx_ref_end - 30 {
                let needed_extra_a_len = 30 + overlap_wdow_end - tx_ref_end;
                // Build tail 30-mers with extra added "A"s from polyA tail
                // at most we add 15 extra A, since when 30mer has >15A, we will consider this as polyA tail. no need to add more& we can be efficient.
                let tail_seq = format!(
                    "{}{}",
                    &ref_seq[tx_ref_end.saturating_sub(30 - 1)..],
                    "A".repeat(needed_extra_a_len.min(15))
                );
                let tail_30mers = build_kmers(&tail_seq, 30);

                let has_enough_a: Array1<bool> = Array1::from(
                    tail_30mers
                        .iter()
                        .map(|kmer| memmem::find(kmer.as_bytes(), pattern_bytes).is_some())
                        .collect::<Vec<bool>>(),
                );

                let tail_binding_affinity = if has_enough_a.iter().any(|&x| x) {
                    process_binding_affinity(
                        &tail_30mers,
                        has_enough_a,
                        mlp,
                        discount_perc,
                        binding_affinity_threshold,
                    )?
                } else {
                    Array1::zeros(tail_30mers.len())
                };

                // Compute joint probabilities for the tail
                let all_a_prob = 1.0;

                // ovlp_wdow_dis_to_tx_end_30mer is the region we already computedin above downstream branch; (overlap_wdow_end  - overlap_wdow_start + 1)is the length of the shared sliding window.
                // let ovlp_wdow_dis_to_tx_end_30mer = tx_ref_end - 30 + 1 - overlap_wdow_start;
                let tx_end_30mer_start = tx_ref_end.saturating_sub(29); // = tx_ref_end - 30 + 1
                if overlap_wdow_start >= tx_end_30mer_start {
                    continue; // skip the overhanging cases
                }
                let ovlp_wdow_dis_to_tx_end_30mer = tx_end_30mer_start - overlap_wdow_start;

                let ovlp_wdow_length = (overlap_wdow_end  - overlap_wdow_start + 1)
                .min(ovlp_wdow_dis_to_tx_end_30mer + 1 + polya_tail_len);

                tail_sum_log_probs.clear();
                tail_sum_log_probs.resize(ovlp_wdow_length, 0.0);
                // tail_joint_prob length is determined per-alignment
                for &ref_start in &ref_start_list {
                    // here we compute the ovlp_start to the last 30 mer of the ref
                    // because this is the cases we did not covered by the above arm(downstream window end is within the ref)

                    let prefix_dis = overlap_wdow_start - ref_start;
                    let start_idx = prefix_dis + ovlp_wdow_dis_to_tx_end_30mer;
                    // NOTE: predfix_dis is the per algn dis;
                    let end_idx = prefix_dis+ ovlp_wdow_length;

                    let tail_frag_len_prob = spline_lookup.slice(s![start_idx..end_idx]);

                    // Update joint probabilities with tail_binding_affinity
                    for (i, frag_prob) in tail_frag_len_prob.iter().enumerate() {
                        let affinity = if i < tail_30mers.len() {
                             tail_binding_affinity[i]
                        } else {
                             all_a_prob // All A probability; if we got >15 A, also apply all A probability, as this is more likely to be polyA tail mode, not internal polyA mode.
                        };
                        let prob = affinity * frag_prob + EPS;
                        tail_sum_log_probs[i] += prob.ln();
                    }
                }

                // Sum and normalize joint probabilities for the tail
                let max_log = tail_sum_log_probs
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max);

                let norm_max_tail_log = max_log / ref_start_list.len() as f64;

                if norm_max_tail_log > norm_sum_joint_prob {
                    norm_sum_joint_prob = norm_max_tail_log;
                }
            }

            // Update mcc_to_tx_prob_dict
            mcc_to_tx_prob_dict.insert(*mcc_idx, norm_sum_joint_prob);
        } else if all_reverse {
            ref_end_list.clear();
            ref_end_list.reserve(algn_tuple_list.len());

            for &(_is_reverse, ref_start_u32) in algn_tuple_list.iter() {
                let mut ref_start = ref_start_u32 as usize;
                if ref_start_u32 > invalid_pos_limit {
                // mimic softclip
                    ref_start = 0 as usize;
                }else if ref_start > tx_ref_end{
                    eprintln!("Error: ref_start is not close to 2^32, but exceed ref length. This should not happen.");
                    mcc_invalid = true;
                    break;
                }
                // clamp to transcript end to avoid panics near the end (or for clipped alignments)
                let ref_end = ref_start
                    .saturating_add(read_length as usize)
                    .min(tx_ref_end);
                ref_end_list.push(ref_end);
            }

            if mcc_invalid{
                continue;
            }

            let overlap_wdow_end = *ref_end_list.iter().min().unwrap();
            let overlap_wdow_start = ref_end_list
                .iter()
                .cloned()
                .max()
                .unwrap()
                .saturating_sub(max_frag_len);
            if overlap_wdow_start >= overlap_wdow_end {
                continue;
            }

            if overlap_wdow_end > 30 {
                let start_pos = overlap_wdow_start.max(30);
                let end_pos = overlap_wdow_end.min(tx_ref_end);
                if end_pos <= start_pos {
                    continue;
                }
                let seq_slice = &ref_seq[start_pos..end_pos];
                let rev_comp_seq = reverse_complement(seq_slice).unwrap();

                // Build kmers
                let upstream_30mers = build_kmers(&rev_comp_seq, 30);

                // Check for has_enough_a
                let has_enough_a: Array1<bool> = Array1::from(
                    upstream_30mers
                        .iter()
                        .map(|kmer| memmem::find(kmer.as_bytes(), pattern_bytes).is_some())
                        .collect::<Vec<bool>>(),
                );

                if has_enough_a.iter().any(|&x| x) {
                    let upstream_binding_affinity = process_binding_affinity(
                        &upstream_30mers,
                        has_enough_a,
                        mlp,
                        discount_perc,
                        binding_affinity_threshold,
                    )?;
                    // New: add penalty for antisense reads
                    let anti_sense_penalty = 0.8;
                    // For each alignment, compute fragment length probabilities
                    sum_log_probs.clear();
                    sum_log_probs.resize(upstream_30mers.len(), 0.0);
                    for &ref_end in &ref_end_list {
                        let suffix_dis = ref_end - overlap_wdow_end;
                        let start_idx = suffix_dis + 1;
                        let end_idx = upstream_30mers.len() + suffix_dis + 1;
                        let upstream_frag_len_prob = spline_lookup.slice(s![start_idx..end_idx]);
                        for (i, (&affinity, &frag_prob)) in upstream_binding_affinity
                            .iter()
                            .zip(upstream_frag_len_prob.iter())
                            .enumerate()
                        {
                            sum_log_probs[i] += (affinity * frag_prob * anti_sense_penalty + EPS).ln();
                        }
                    }
                    let max_sum = sum_log_probs.iter().cloned().fold(0.0_f64, f64::max);
                    let norm_sum_joint_prob = max_sum / ref_end_list.len() as f64;
                    // Update mcc_to_tx_prob_dict
                    mcc_to_tx_prob_dict.insert(*mcc_idx, norm_sum_joint_prob);
                }
            }else{
                continue;
            }
        }else{
            eprintln!("Error: algn_tuple_list is not all forward or all reverse. This should not happen.");
            continue;
        }
    }

    // if invalid_pos_skipped > 0 {
    //     eprintln!(
    //         "Forseti summary: skipped {} MCC(s) due to invalid alignment start pos (printed first {}).",
    //         invalid_pos_skipped, INVALID_POS_PRINT_LIMIT,
    //     );
    // }

    // Process mcc_to_tx_prob_dict to find best indices
    let mcc_pred_score_list: Vec<(MCCIndex, f64)> = mcc_to_tx_prob_dict
        .iter()
        .map(|(&mcc_idx, &score)| (mcc_idx, score))
        .collect();

    if mcc_pred_score_list.is_empty() {
        return Ok(Vec::new());
    }

    let (_best_score, best_mcc_indices): (f64, Vec<u16>) = mcc_pred_score_list.iter().fold(
        (f64::MIN, Vec::new()),
        |(max_score, mut indices), &(mcc_idx, score)| {
            if score > max_score {
                (score, vec![mcc_idx.try_into().unwrap()]) // New max score found, reset indices
            } else if (score - max_score).abs() < f64::EPSILON {
                indices.push(mcc_idx.try_into().unwrap()); // Tie with current max, add index
                (max_score, indices)
            } else {
                (max_score, indices) // No change
            }
        },
    );
    // println!("Best score: {}", best_score);

    Ok(best_mcc_indices)
}
