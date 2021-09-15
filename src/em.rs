/*
 * Copyright (c) 2020-2021 Rob Patro, Avi Srivastava, Hirak Sarkar, Dongze He, Mohsen Zakeri.
 *
 * This file is part of alevin-fry
 * (see https://github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

#[allow(unused_imports)]
use crate::eq_class::IndexedEqList;
#[allow(unused_imports)]
use ahash::{AHasher, RandomState};
use rand::{thread_rng, Rng};
#[allow(unused_imports)]
use slog::info;
use statrs::distribution::Multinomial;
use std::collections::HashMap;
use std::f32;

//#[derive(Clone, Debug)]
//pub struct SalmonEQClass {
//    pub labels: Vec<u32>,
//    pub counts: u32,
//}

const MIN_ALPHA: f32 = 1e-8;
const MIN_OUTPUT_ALPHA: f32 = 0.01;
const ALPHA_CHECK_CUTOFF: f32 = 1e-2;

const MIN_ITER: u32 = 2;
const MAX_ITER: u32 = 100;
const REL_DIFF_TOLERANCE: f32 = 1e-2;

#[derive(Copy, Clone)]
pub enum EmInitType {
    Informative,
    Uniform,
    Random,
}

#[allow(dead_code)]
fn mean(data: &[f64]) -> Option<f64> {
    let sum = data.iter().sum::<f64>() as f64;
    let count = data.len();

    match count {
        positive if positive > 0 => Some(sum / count as f64),
        _ => None,
    }
}

#[allow(dead_code)]
fn std_deviation(data: &[f64]) -> Option<f64> {
    match (mean(data), data.len()) {
        (Some(data_mean), count) if count > 0 => {
            let variance = data
                .iter()
                .map(|value| {
                    let diff = data_mean - (*value as f64);

                    diff * diff
                })
                .sum::<f64>()
                / count as f64;

            Some(variance.sqrt())
        }
        _ => None,
    }
}

pub(crate) fn em_update_subset(
    alphas_in: &[f32],
    alphas_out: &mut Vec<f32>,
    eqclasses: &IndexedEqList,
    cell_data: &[(u32, u32)], // indices into eqclasses relevant for this cell
) {
    for (i, count) in cell_data {
        let labels = eqclasses.refs_for_eqc(*i);

        if labels.len() > 1 {
            let mut denominator: f32 = 0.0;
            for label in labels {
                denominator += alphas_in[*label as usize];
            }

            if denominator > 0.0 {
                let inv_denominator = *count as f32 / denominator;
                for label in labels {
                    let index = *label as usize;
                    let count = alphas_in[index] * inv_denominator;
                    alphas_out[index] += count;
                }
            }
        } else {
            let tidx = labels.get(0).expect("can't extract labels");
            alphas_out[*tidx as usize] += *count as f32;
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn em_optimize_subset(
    eqclasses: &IndexedEqList,
    cell_data: &[(u32, u32)], // indices into eqclasses relevant for this cell
    unique_evidence: &mut Vec<bool>,
    no_ambiguity: &mut Vec<bool>,
    init_type: EmInitType,
    num_alphas: usize,
    only_unique: bool,
    _log: &slog::Logger,
) -> Vec<f32> {
    let mut alphas_in: Vec<f32> = vec![0.0; num_alphas];
    let mut alphas_out: Vec<f32> = vec![0.0; num_alphas];

    let mut needs_em = false;
    for (i, count) in cell_data {
        let labels = eqclasses.refs_for_eqc(*i);
        if labels.len() == 1 {
            let idx = labels.get(0).expect("can't extract labels");
            alphas_in[*idx as usize] += *count as f32;
            unique_evidence[*idx as usize] = true;
        } else {
            for idx in labels {
                no_ambiguity[*idx as usize] = false;
            }
            needs_em = true;
        }
    }

    // if we are just pulling out unique counts
    // or there were no multi-mapping reads, then
    // we're done
    if only_unique || !needs_em {
        return alphas_in;
    }

    // fill in the alphas based on the initialization strategy
    let mut rng = rand::thread_rng();
    let uni_prior = 1.0 / (num_alphas as f32);
    for item in alphas_in.iter_mut().take(num_alphas) {
        match init_type {
            EmInitType::Uniform => {
                *item = uni_prior;
            }
            EmInitType::Informative => {
                *item = (*item + 0.5) * 1e-3;
            }
            EmInitType::Random => {
                *item = rng.gen::<f32>() + 1e-5;
            }
        }
    }

    let mut it_num: u32 = 0;
    let mut converged: bool = true;
    // allow one last round of the EM after thresholding
    // very small counts to 0.
    let mut last_round: bool = false;

    while it_num < MIN_ITER || (it_num < MAX_ITER && !converged) || last_round {
        // perform one round of em update
        em_update_subset(&alphas_in, &mut alphas_out, eqclasses, cell_data);

        converged = true;
        let mut max_rel_diff = -f32::INFINITY;

        for index in 0..num_alphas {
            if alphas_out[index] > ALPHA_CHECK_CUTOFF {
                let diff = alphas_in[index] - alphas_out[index];
                let rel_diff = diff.abs();

                max_rel_diff = max_rel_diff.max(rel_diff);

                if rel_diff > REL_DIFF_TOLERANCE {
                    converged = false;
                }
            } // end- in>out if

            alphas_in[index] = alphas_out[index];
            alphas_out[index] = 0.0_f32;
        } //end-for

        it_num += 1;

        // if this was the last round
        // then break the loop.
        if last_round {
            break;
        }

        // if we've run for at least the required number
        // of iterations, and if we are converged
        // then do one last round after filtering
        // very small values.
        if it_num >= MIN_ITER && converged {
            alphas_in.iter_mut().for_each(|alpha| {
                if *alpha < MIN_OUTPUT_ALPHA {
                    *alpha = 0.0_f32;
                }
            });
            last_round = true;
        }
    }

    // update too small alphas
    alphas_in.iter_mut().for_each(|alpha| {
        if *alpha < MIN_OUTPUT_ALPHA {
            *alpha = 0.0_f32;
        }
    });

    //let alphas_sum: f32 = alphas_in.iter().sum();
    //assert!(alphas_sum > 0.0, "Alpha Sum too small");
    alphas_in
}

pub fn em_update(
    alphas_in: &[f32],
    alphas_out: &mut Vec<f32>,
    eqclasses: &HashMap<Vec<u32>, u32, ahash::RandomState>,
) {
    // loop over all the eqclasses
    for (labels, count) in eqclasses {
        if labels.len() > 1 {
            let mut denominator: f32 = 0.0;
            for label in labels {
                denominator += alphas_in[*label as usize];
            }

            if denominator > 0.0 {
                let inv_denominator = *count as f32 / denominator;
                for label in labels {
                    let index = *label as usize;
                    let count = alphas_in[index] * inv_denominator;
                    alphas_out[index] += count;
                }
            }
        } else {
            let tidx = labels.get(0).expect("can't extract labels");
            alphas_out[*tidx as usize] += *count as f32;
        }
    }
}

pub fn em_optimize(
    eqclasses: &HashMap<Vec<u32>, u32, ahash::RandomState>,
    unique_evidence: &mut Vec<bool>,
    no_ambiguity: &mut Vec<bool>,
    init_type: EmInitType,
    num_alphas: usize,
    only_unique: bool,
    _log: &slog::Logger,
) -> Vec<f32> {
    let mut alphas_in: Vec<f32> = vec![0.0; num_alphas];
    let mut alphas_out: Vec<f32> = vec![0.0; num_alphas];

    for (labels, count) in eqclasses {
        if labels.len() == 1 {
            let idx = labels.get(0).expect("can't extract labels");
            alphas_in[*idx as usize] += *count as f32;
            unique_evidence[*idx as usize] = true;
        } else {
            for idx in labels {
                no_ambiguity[*idx as usize] = false;
            }
        }
    }

    if only_unique {
        return alphas_in;
    }

    // fill in the alphas based on the initialization strategy
    let mut rng = rand::thread_rng();
    let uni_prior = 1.0 / (num_alphas as f32);
    for item in alphas_in.iter_mut().take(num_alphas) {
        match init_type {
            EmInitType::Uniform => {
                *item = uni_prior;
            }
            EmInitType::Informative => {
                *item = (*item + 0.5) * 1e-3;
            }
            EmInitType::Random => {
                *item = rng.gen::<f32>() + 1e-5;
            }
        }
    }

    // TODO: is it even necessary?
    //alphas_in.iter_mut().for_each(|alpha| *alpha *= 1e-3);

    let mut it_num: u32 = 0;
    let mut converged: bool = true;
    while it_num < MIN_ITER || (it_num < MAX_ITER && !converged) {
        // perform one round of em update
        em_update(&alphas_in, &mut alphas_out, eqclasses);

        converged = true;
        let mut max_rel_diff = -f32::INFINITY;

        for index in 0..num_alphas {
            if alphas_out[index] > ALPHA_CHECK_CUTOFF {
                let diff = alphas_in[index] - alphas_out[index];
                let rel_diff = diff.abs();

                max_rel_diff = match rel_diff > max_rel_diff {
                    true => rel_diff,
                    false => max_rel_diff,
                };

                if rel_diff > REL_DIFF_TOLERANCE {
                    converged = false;
                }
            } // end- in>out if

            alphas_in[index] = alphas_out[index];
            alphas_out[index] = 0.0_f32;
        } //end-for

        it_num += 1;
    }

    // update too small alphas
    alphas_in.iter_mut().for_each(|alpha| {
        if *alpha < MIN_OUTPUT_ALPHA {
            *alpha = 0.0_f32;
        }
    });
    //let alphas_sum: f32 = alphas_in.iter().sum();
    //assert!(alphas_sum > 0.0, "Alpha Sum too small");
    /*
    info!(log,
    "Total Molecules after EM {}",
    alphas_sum
    );
    */
    alphas_in
}

pub(crate) fn run_bootstrap_subset(
    eqclasses: &IndexedEqList,
    cell_data: &[(u32, u32)], // (eq_id, count) vec for classes relevant for this cell
    num_alphas: u32,          // number of genes
    num_bootstraps: u32,      // number of bootstraps to draw
    _init_uniform: bool,
    summary_stat: bool, // if true, the output will simply be a vector of means and variances
    _log: &slog::Logger,
) -> Vec<Vec<f32>> {
    // the population sample size
    let total_fragments: u32 = cell_data.iter().map(|x| x.1).sum();
    assert!(
        total_fragments > 0,
        "Cannot bootstrap from a sample with 0 counts."
    );

    let num_alphas_us = num_alphas as usize;
    let mut unique_evidence = vec![false; num_alphas_us];
    let mut no_ambiguity = vec![false; num_alphas_us];

    let mut alphas_sum: Vec<f32> = vec![0.0; num_alphas_us];
    let mut alphas_square_sum: Vec<f32> = vec![0.0; num_alphas_us];
    let mut sample_mean: Vec<f32> = vec![0.0; num_alphas_us];
    let mut sample_var: Vec<f32> = vec![0.0; num_alphas_us];

    // define a multinomial with the probabilities given by the
    // original equivalence class counts
    let eq_counts: Vec<f64> = cell_data
        .iter()
        .map(|x| (x.1 as f64) / (total_fragments as f64))
        .collect();
    let dist = Multinomial::new(&eq_counts[..], total_fragments as u64).unwrap();

    // store bootstraps
    let mut bootstrap_counts = Vec::with_capacity(cell_data.len());

    // store bootstraps
    let num_output_bs = if summary_stat {
        2usize
    } else {
        num_bootstraps as usize
    };
    let mut bootstraps = Vec::with_capacity(num_output_bs);

    // bootstrap loop starts
    // let mut old_resampled_counts = Vec::new();
    for _bs_num in 0..num_bootstraps {
        // resample from multinomial
        let resampled_counts = thread_rng().sample(dist.clone());
        for (idx, (eq_id, _orig_count)) in cell_data.iter().enumerate() {
            bootstrap_counts.push((*eq_id, resampled_counts[idx].round() as u32));
        }

        let alphas = em_optimize_subset(
            &eqclasses,
            &bootstrap_counts[..], // indices into eqclasses relevant for this cell
            &mut unique_evidence,
            &mut no_ambiguity,
            EmInitType::Random,
            num_alphas_us,
            false, // only unique
            &_log,
        );

        // clear out for the next iteration.
        bootstrap_counts.clear();

        let est_frags: f32 = alphas.iter().sum();
        assert!(est_frags > 0.0, "Alpha sum is too small");
        // if we are collecting summary stats, we'll need these
        if summary_stat {
            for i in 0..num_alphas_us {
                alphas_sum[i] += alphas[i];
                alphas_square_sum[i] += alphas[i] * alphas[i];
            }
        } else {
            // otherwise, just push the bootstrap
            bootstraps.push(alphas.clone());
        }
    }

    // if we are only providing summary stats, then
    // do that computation here.
    if summary_stat {
        for i in 0..num_alphas_us {
            let mean_alpha = alphas_sum[i] / num_bootstraps as f32;
            sample_mean[i] = mean_alpha;
            sample_var[i] =
                (alphas_square_sum[i] / num_bootstraps as f32) - (mean_alpha * mean_alpha);
        }

        bootstraps.push(sample_mean);
        bootstraps.push(sample_var);
    }

    bootstraps
}

pub fn run_bootstrap(
    eqclasses: &HashMap<Vec<u32>, u32, ahash::RandomState>,
    num_bootstraps: u32,
    gene_alpha: &[f32],
    // unique_evidence: &mut Vec<bool>,
    // no_ambiguity: &mut Vec<bool>,
    // num_alphas: usize,
    // only_unique: bool,
    _init_uniform: bool,
    summary_stat: bool,
    _log: &slog::Logger,
) -> Vec<Vec<f32>> {
    // This function is just a thin wrapper around run_bootstrap_subset.
    // Here, we convert the hashmap to an `IndexedEqList`, we map the
    // counts to a `Vec<(u32, u32)>` representing the equivalence class
    // indices and counts for this cell's data.

    // NOTE: This working properly relies on iteration over a static hash
    // yielding the key value pairs in the same order.  This isn't generally
    // true between runs (b/c of how rust initalizes hashes), but within the
    // same run and over the same hash table (as here), it should always hold.
    let eql = IndexedEqList::init_from_hash(eqclasses, gene_alpha.len());

    // since this is a local (not global) eq-list, the cell data is just
    // the indices of all equivalence classes and their corresponding counts.
    let cell_data: Vec<(u32, u32)> = eqclasses
        .iter()
        .enumerate()
        .map(|(idx, (_labels, count))| (idx as u32, *count))
        .collect();

    // now that we have the `IndexedEqList` representation of this data, just
    // run that version of the bootstrap function and return the result.
    run_bootstrap_subset(
        &eql,
        &cell_data[..],
        gene_alpha.len() as u32,
        num_bootstraps,
        _init_uniform,
        summary_stat,
        _log,
    )
}

#[allow(dead_code)]
pub fn run_bootstrap_old(
    eqclasses: &HashMap<Vec<u32>, u32, ahash::RandomState>,
    num_bootstraps: u32,
    gene_alpha: &[f32],
    // unique_evidence: &mut Vec<bool>,
    // no_ambiguity: &mut Vec<bool>,
    // num_alphas: usize,
    // only_unique: bool,
    init_uniform: bool,
    summary_stat: bool,
    _log: &slog::Logger,
) -> Vec<Vec<f32>> {
    let mut total_fragments = 0u64;

    // println!("In bootstrapping");

    let mut alphas: Vec<f32> = vec![0.0; gene_alpha.len()];
    let mut alphas_mean: Vec<f32> = vec![0.0; gene_alpha.len()];
    let mut alphas_square: Vec<f32> = vec![0.0; gene_alpha.len()];
    let mut sample_mean: Vec<f32> = vec![0.0; gene_alpha.len()];
    let mut sample_var: Vec<f32> = vec![0.0; gene_alpha.len()];
    let mut alphas_prime: Vec<f32> = vec![0.0; gene_alpha.len()];
    // let mut means: Vec<f32> = vec![0.0; gene_alpha.len()];
    // let mut square_means: Vec<f32> = vec![0.0; gene_alpha.len()];

    // make discrete distribution of the eqclass counts
    // hash map to serialize the eqclasses
    // eqclasses_serialize : id -> label
    // eqclasses : label -> count
    let mut eq_counts = Vec::new();
    let mut eqclasses_serialize: HashMap<usize, Vec<u32>> = HashMap::new();

    for (idx, (labels, count)) in eqclasses.iter().enumerate() {
        eq_counts.push(*count as f64);
        total_fragments += *count as u64;
        eqclasses_serialize
            .entry(idx)
            .or_insert_with(|| labels.to_vec());
    }

    // println!("total fragments {:?}", total_fragments);

    // a new hashmap to be updated in each bootstrap
    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut eqclass_bootstrap: HashMap<Vec<u32>, u32, ahash::RandomState> = HashMap::with_hasher(s);
    // define a multinomial
    let dist = Multinomial::new(&eq_counts, total_fragments).unwrap();

    // store bootstraps
    let mut bootstraps = Vec::new();

    // bootstrap loop starts
    // let mut old_resampled_counts = Vec::new();
    for _bs_num in 0..num_bootstraps {
        // resample from multinomial
        let resampled_counts = thread_rng().sample(dist.clone());

        for (eq_id, labels) in &eqclasses_serialize {
            eqclass_bootstrap.insert(labels.to_vec(), resampled_counts[*eq_id].round() as u32);
            // eqclass_bootstrap
            //     .entry(labels.to_vec())
            //     .or_insert(resampled_counts[*eq_id].round() as u32);
        }

        // fill new alpha
        for i in 0..gene_alpha.len() {
            if init_uniform {
                alphas[i] = 1.0 / gene_alpha.len() as f32;
            } else {
                alphas[i] = (gene_alpha[i] + 0.5) * 1e-3;
            }
        }

        // let alpha_sum : f32 = alphas.iter().sum();
        // println!("Bootstrap num ... {:?}, alpha_sum ... {:?}", _bs_num, alpha_sum);

        let mut it_num: u32 = 0;
        let mut converged: bool = false;
        while it_num < MIN_ITER || (it_num < MAX_ITER && !converged) {
            // perform one round of em update
            em_update(&alphas, &mut alphas_prime, &eqclass_bootstrap);

            converged = true;
            let mut max_rel_diff = -f32::INFINITY;

            for index in 0..gene_alpha.len() {
                if alphas_prime[index] > ALPHA_CHECK_CUTOFF {
                    let diff = alphas[index] - alphas_prime[index];
                    let rel_diff = diff.abs();

                    max_rel_diff = match rel_diff > max_rel_diff {
                        true => rel_diff,
                        false => max_rel_diff,
                    };

                    if rel_diff > REL_DIFF_TOLERANCE {
                        converged = false;
                    }
                } // end- in>out if

                alphas[index] = alphas_prime[index];
                alphas_prime[index] = 0.0_f32;
            } //end-for

            it_num += 1;
        }

        // update too small alphas
        alphas.iter_mut().for_each(|alpha| {
            if *alpha < MIN_ALPHA {
                *alpha = 0.0_f32;
            }
        });

        let alphas_sum: f32 = alphas.iter().sum();
        assert!(alphas_sum > 0.0, "Alpha Sum too small");
        if summary_stat {
            for i in 0..gene_alpha.len() {
                alphas_mean[i] += alphas[i];
                alphas_square[i] += alphas[i] * alphas[i];
            }
        } else {
            bootstraps.push(alphas.clone());
        }
        // println!("After alpha sum: {:?}, it_num: {:?}", alphas_sum, it_num);
        // old_resampled_counts = resampled_counts.clone();
    }
    if summary_stat {
        for i in 0..gene_alpha.len() {
            let mean_alpha = alphas_mean[i] / num_bootstraps as f32;
            sample_mean[i] = mean_alpha;
            sample_var[i] = (alphas_square[i] / num_bootstraps as f32) - (mean_alpha * mean_alpha);
        }

        bootstraps.push(sample_mean);
        bootstraps.push(sample_var);
    }

    bootstraps
}
