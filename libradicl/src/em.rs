// Copyright 2020 Rob Patro, Avi Srivastava. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

extern crate fasthash;
extern crate slog;

#[allow(unused_imports)]
use self::slog::info;
use fasthash::sea::Hash64;
#[allow(unused_imports)]
use fasthash::RandomState;
use std::collections::HashMap;
use std::f32;

//#[derive(Clone, Debug)]
//pub struct SalmonEQClass {
//    pub labels: Vec<u32>,
//    pub counts: u32,
//}

const MIN_ALPHA: f32 = 1e-8;
const ALPHA_CHECK_CUTOFF: f32 = 1e-2;

const MIN_ITER: u32 = 50;
const MAX_ITER: u32 = 10_000;
const REL_DIFF_TOLERANCE: f32 = 1e-2;

pub fn em_update(
    alphas_in: &[f32],
    alphas_out: &mut Vec<f32>,
    eqclasses: &HashMap<Vec<u32>, u32, fasthash::RandomState<Hash64>>,
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
    eqclasses: &HashMap<Vec<u32>, u32, fasthash::RandomState<Hash64>>,
    unique_evidence: &mut Vec<bool>,
    no_ambiguity: &mut Vec<bool>,
    num_alphas: usize,
    only_unique: bool,
    _log: &slog::Logger,
) -> Vec<f32> {
    // set up starting alphas as 0.5
    let mut alphas_in: Vec<f32> = vec![0.5; num_alphas];
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
        alphas_in.iter_mut().for_each(|alpha| {
            *alpha -= 0.5;
        });

        return alphas_in;
    }

    // TODO: is it even necessary?
    alphas_in.iter_mut().for_each(|alpha| *alpha *= 1e-3);

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
            alphas_out[index] = 0.0 as f32;
        } //end-for

        it_num += 1;
    }

    // update too small alphas
    alphas_in.iter_mut().for_each(|alpha| {
        if *alpha < MIN_ALPHA {
            *alpha = 0.0 as f32;
        }
    });

    let alphas_sum: f32 = alphas_in.iter().sum();
    assert!(alphas_sum > 0.0, "Alpha Sum too small");
    /*
    info!(log,
        "Total Molecules after EM {}",
        alphas_sum
    );
    */
    alphas_in
}
