/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

#[allow(unused_imports)]
use crate::eq_class::IndexedEqList;
use crate::multinomial::Multinomial;
use crate::utils::EqClassPayload;
#[allow(unused_imports)]
use ahash::{AHasher, RandomState};
use rand::RngExt;
#[allow(unused_imports)]
use slog::info;
use std::collections::HashMap;
use std::f32;

//#[derive(Clone, Debug)]
//pub struct SalmonEQClass {
//    pub labels: Vec<u32>,
//    pub counts: u32,
//}

const _MIN_ALPHA: f32 = 1e-8;
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

/// Reusable storage for [`em_optimize_subset_with_scratch`].
///
/// Create one value per worker thread and reuse it across cells. The
/// bit-packed membership vector lets support construction append each alpha
/// index at most once without first collecting and sorting duplicate labels.
/// Dense alpha workspaces deliberately remain invocation-local so idle workers
/// do not retain a full gene-axis allocation.
#[derive(Debug, Default)]
pub(crate) struct EmSubsetScratch {
    support: Vec<u32>,
    support_membership: Vec<bool>,
}

impl EmSubsetScratch {
    fn new(num_alphas: usize) -> Self {
        Self {
            support: Vec::new(),
            support_membership: vec![false; num_alphas],
        }
    }

    fn clear_support(&mut self) {
        for &index in &self.support {
            self.support_membership[index as usize] = false;
        }
        self.support.clear();
    }

    #[inline]
    fn mark_supported(&mut self, index: usize, num_alphas: usize) {
        assert!(
            index < num_alphas,
            "equivalence-class label {index} is outside {num_alphas} alphas"
        );
        if !self.support_membership[index] {
            self.support_membership[index] = true;
            self.support
                .push(u32::try_from(index).expect("alpha index does not fit in u32"));
        }
    }

    /// Populate the possible support for a cell. Counts are intentionally not
    /// consulted: the resulting superset can be reused when bootstrap
    /// replicates assign zero observations to some equivalence classes.
    fn prepare_support(
        &mut self,
        eqclasses: &IndexedEqList,
        cell_data: &[(u32, u32)],
        num_alphas: usize,
        usa_offsets: Option<(usize, usize)>,
    ) {
        self.clear_support();
        self.support_membership.resize(num_alphas, false);

        for (eq_id, _) in cell_data {
            for &label in eqclasses.refs_for_eqc(*eq_id) {
                let index = label as usize;
                self.mark_supported(index, num_alphas);
                if let Some((unspliced_offset, ambig_offset)) = usa_offsets {
                    if index >= ambig_offset {
                        self.mark_supported(index - unspliced_offset, num_alphas);
                        self.mark_supported(index - ambig_offset, num_alphas);
                    } else if index >= unspliced_offset {
                        self.mark_supported(index + unspliced_offset, num_alphas);
                    } else {
                        self.mark_supported(index + ambig_offset, num_alphas);
                    }
                }
            }
        }
    }
}

/// The subset EM's dense storage, retained by one inference worker. Only
/// indices in the previous cell's possible support need clearing. This keeps
/// `infer` from allocating, zeroing and scanning the entire feature axis for
/// every cell; the EM arithmetic and stopping criteria remain unchanged.
#[derive(Default)]
pub(crate) struct SubsetInferenceCounts {
    support: EmSubsetScratch,
    values: Vec<f32>,
    next: Vec<f32>,
    unique_evidence: Vec<bool>,
    no_ambiguity: Vec<bool>,
}

impl SubsetInferenceCounts {
    pub(crate) fn optimize(
        &mut self,
        eqclasses: &IndexedEqList,
        cell_data: &[(u32, u32)],
        num_alphas: usize,
        usa_offsets: Option<(usize, usize)>,
    ) {
        // Clear before resizing: the axis can shrink between invocations.
        for &index in &self.support.support {
            let i = index as usize;
            self.values[i] = 0.0;
            self.next[i] = 0.0;
            self.unique_evidence[i] = false;
            self.no_ambiguity[i] = false;
        }
        self.values.resize(num_alphas, 0.0);
        self.next.resize(num_alphas, 0.0);
        self.unique_evidence.resize(num_alphas, false);
        self.no_ambiguity.resize(num_alphas, false);
        // Include the sibling statuses read by USA denominators, even if a
        // sibling is absent from every label or has zero count in this cell.
        self.support
            .prepare_support(eqclasses, cell_data, num_alphas, usa_offsets);
        em_optimize_subset_into(
            eqclasses,
            cell_data,
            &mut self.unique_evidence,
            &mut self.no_ambiguity,
            EmInitType::Informative,
            num_alphas,
            false,
            usa_offsets,
            &mut self.support,
            true,
            &mut self.values,
            &mut self.next,
        );
    }

    pub(crate) fn write_nonzero(&mut self, indices: &mut Vec<usize>, values: &mut Vec<f32>) {
        indices.clear();
        values.clear();
        // Preserve the original dense scan's increasing feature order.
        self.support.support.sort_unstable();
        for &index in &self.support.support {
            let value = self.values[index as usize];
            if value > 0.0 {
                indices.push(index as usize);
                values.push(value);
            }
        }
    }
}

#[allow(dead_code)]
fn mean(data: &[f64]) -> Option<f64> {
    let sum = data.iter().sum::<f64>();
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
                    let diff = data_mean - *value;

                    diff * diff
                })
                .sum::<f64>()
                / count as f64;

            Some(variance.sqrt())
        }
        _ => None,
    }
}

#[allow(dead_code)]
#[inline(always)]
fn spliced_idx_of(idx: u32, usa_offset: (usize, usize)) -> u32 {
    // if this index is ambiguous, subtract off
    // the ambig offset
    if idx as usize >= usa_offset.1 {
        idx - (usa_offset.1 as u32)
    } else if idx as usize >= usa_offset.0 {
        // otherwise, if it is unspliced, subtract
        // off the unspliced offset
        idx - (usa_offset.0 as u32)
    } else {
        // otherwise it is spliced, return it as is
        idx
    }
}

/// Based on the type of the index `idx` (spliced, unspliced or ambiguous)
/// get the appropriate count for this gene, status pair to use in the EM.
#[inline(always)]
fn get_abundance_for(idx: u32, alphas_in: &[f32], usa_offset: (usize, usize)) -> f32 {
    let index = idx as usize;
    // if this index is ambiguous, subtract off
    // the ambig offset
    if index >= usa_offset.1 {
        // in the ambiguous case, get the abundance
        // of spliced / unspliced / and ambiguous
        alphas_in[index - usa_offset.0] + // unspliced
        alphas_in[index - usa_offset.1] + // spliced
        alphas_in[index] // ambig
    } else if index >= usa_offset.0 {
        // in the unspliced case, get the abundance of
        // unspliced and ambiguous
        alphas_in[index + usa_offset.0] + // ambiguous
        alphas_in[index] // unspliced
    } else {
        // otherwise it is spliced, return it as is
        alphas_in[index + usa_offset.1] + // ambiguous
        alphas_in[index] // spliced
    }
}

pub(crate) fn em_update_subset(
    alphas_in: &[f32],
    alphas_out: &mut [f32],
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
            let tidx = labels.first().expect("can't extract labels");
            alphas_out[*tidx as usize] += *count as f32;
        }
    }
}

pub(crate) fn em_update_subset_usa(
    alphas_in: &[f32],
    alphas_out: &mut [f32],
    eqclasses: &IndexedEqList,
    cell_data: &[(u32, u32)], // indices into eqclasses relevant for this cell
    usa_offsets: (usize, usize),
) {
    for (i, count) in cell_data {
        let labels = eqclasses.refs_for_eqc(*i);

        if labels.len() > 1 {
            let mut denominator: f32 = 0.0;
            for label in labels {
                denominator += get_abundance_for(*label, alphas_in, usa_offsets);
            }

            if denominator > 0.0 {
                let inv_denominator = *count as f32 / denominator;
                for label in labels {
                    let count = get_abundance_for(*label, alphas_in, usa_offsets) * inv_denominator;
                    let index = *label as usize;
                    alphas_out[index] += count;
                }
            }
        } else {
            let tidx = labels.first().expect("can't extract labels");
            alphas_out[*tidx as usize] += *count as f32;
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn em_optimize_subset(
    eqclasses: &IndexedEqList,
    cell_data: &[(u32, u32)], // indices into eqclasses relevant for this cell
    unique_evidence: &mut [bool],
    no_ambiguity: &mut [bool],
    init_type: EmInitType,
    num_alphas: usize,
    only_unique: bool,
    usa_offsets: Option<(usize, usize)>,
    _log: &slog::Logger,
) -> Vec<f32> {
    let mut scratch = EmSubsetScratch::new(num_alphas);
    em_optimize_subset_with_scratch(
        eqclasses,
        cell_data,
        unique_evidence,
        no_ambiguity,
        init_type,
        num_alphas,
        only_unique,
        usa_offsets,
        _log,
        &mut scratch,
    )
}

/// Optimize one cell while reusing worker-local allocation buffers.
#[allow(clippy::too_many_arguments)]
pub(crate) fn em_optimize_subset_with_scratch(
    eqclasses: &IndexedEqList,
    cell_data: &[(u32, u32)],
    unique_evidence: &mut [bool],
    no_ambiguity: &mut [bool],
    init_type: EmInitType,
    num_alphas: usize,
    only_unique: bool,
    usa_offsets: Option<(usize, usize)>,
    _log: &slog::Logger,
    scratch: &mut EmSubsetScratch,
) -> Vec<f32> {
    em_optimize_subset_impl(
        eqclasses,
        cell_data,
        unique_evidence,
        no_ambiguity,
        init_type,
        num_alphas,
        only_unique,
        usa_offsets,
        scratch,
        false,
    )
}

#[allow(clippy::too_many_arguments)]
fn em_optimize_subset_impl(
    eqclasses: &IndexedEqList,
    cell_data: &[(u32, u32)],
    unique_evidence: &mut [bool],
    no_ambiguity: &mut [bool],
    init_type: EmInitType,
    num_alphas: usize,
    only_unique: bool,
    usa_offsets: Option<(usize, usize)>,
    scratch: &mut EmSubsetScratch,
    support_is_prepared: bool,
) -> Vec<f32> {
    let mut alphas_in = vec![0.0; num_alphas];
    let mut alphas_out = vec![0.0; num_alphas];
    em_optimize_subset_into(
        eqclasses,
        cell_data,
        unique_evidence,
        no_ambiguity,
        init_type,
        num_alphas,
        only_unique,
        usa_offsets,
        scratch,
        support_is_prepared,
        &mut alphas_in,
        &mut alphas_out,
    );
    alphas_in
}

#[allow(clippy::too_many_arguments)]
fn em_optimize_subset_into(
    eqclasses: &IndexedEqList,
    cell_data: &[(u32, u32)],
    unique_evidence: &mut [bool],
    no_ambiguity: &mut [bool],
    init_type: EmInitType,
    num_alphas: usize,
    only_unique: bool,
    usa_offsets: Option<(usize, usize)>,
    scratch: &mut EmSubsetScratch,
    support_is_prepared: bool,
    alphas_in: &mut [f32],
    alphas_out: &mut [f32],
) {
    let mut needs_em = false;
    for (i, count) in cell_data {
        let labels = eqclasses.refs_for_eqc(*i);
        if labels.len() == 1 {
            let idx = labels.first().expect("can't extract labels");
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
        return;
    }

    // The alpha vectors are dense over every gene (or every gene x status in
    // USA mode), but a single cell only ever touches the handful of them that
    // appear in its own equivalence classes. Walking all `num_alphas` entries
    // once per EM iteration made the per-cell cost independent of the size of
    // the cell: in USA mode, a run with 78,899 genes has 236,697 alpha values
    // to compare per iteration even when a cell has only five equivalence
    // classes.
    //
    // `support` is the set of indices this cell can touch. `em_update_subset`
    // only writes to the labels of the cell's equivalence classes, and in USA
    // mode `get_abundance_for` additionally *reads* the two sibling statuses of
    // each label, so those are collected too: an index that is read but never
    // written still has to be reset each iteration for the result to be
    // identical to the dense version.
    if !support_is_prepared {
        scratch.prepare_support(eqclasses, cell_data, num_alphas, usa_offsets);
    }
    let support = &scratch.support;

    // fill in the alphas based on the initialization strategy
    //
    // All initialization modes need only the possible support. In particular,
    // random values outside this set cannot occur in an update denominator or
    // receive output mass, so omitting those independent draws does not change
    // bootstrap semantics.
    let mut rng = rand::rng();
    let uni_prior = 1.0 / (num_alphas as f32);
    for &index in support {
        let item = &mut alphas_in[index as usize];
        match init_type {
            EmInitType::Uniform => {
                *item = uni_prior;
            }
            EmInitType::Informative => {
                *item = (*item + 0.5) * 1e-3;
            }
            EmInitType::Random => {
                *item = rng.random::<f32>() + 1e-5;
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
        match usa_offsets {
            Some(otup) => {
                em_update_subset_usa(alphas_in, alphas_out, eqclasses, cell_data, otup);
            }
            None => {
                em_update_subset(alphas_in, alphas_out, eqclasses, cell_data);
            }
        }

        converged = true;
        let mut max_rel_diff = -f64::INFINITY;

        for &index in support {
            let index = index as usize;
            if alphas_out[index] > ALPHA_CHECK_CUTOFF {
                let diff = alphas_in[index] - alphas_out[index];
                let rel_diff = diff.abs();

                max_rel_diff = max_rel_diff.max(rel_diff as f64);

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
            for &index in support {
                let alpha = &mut alphas_in[index as usize];
                if *alpha < MIN_OUTPUT_ALPHA {
                    *alpha = 0.0_f32;
                }
            }
            last_round = true;
        }
    }

    // update too small alphas; entries outside the support are still zero
    for &index in support {
        let alpha = &mut alphas_in[index as usize];
        if *alpha < MIN_OUTPUT_ALPHA {
            *alpha = 0.0_f32;
        }
    }

    //let alphas_sum: f32 = alphas_in.iter().sum();
    //assert!(alphas_sum > 0.0, "Alpha Sum too small");
}

pub fn em_update<P: EqClassPayload>(
    alphas_in: &[f32],
    alphas_out: &mut [f32],
    eqclasses: &HashMap<Vec<u32>, P, ahash::RandomState>,
) {
    // loop over all the eqclasses
    for (labels, payload) in eqclasses {
        let count = payload.count();
        if labels.len() > 1 {
            let mut denominator: f32 = 0.0;
            for label in labels {
                denominator += alphas_in[*label as usize];
            }

            if denominator > 0.0 {
                let inv_denominator = count as f32 / denominator;
                for label in labels {
                    let index = *label as usize;
                    let count = alphas_in[index] * inv_denominator;
                    alphas_out[index] += count;
                }
            }
        } else {
            let tidx = labels.first().expect("can't extract labels");
            alphas_out[*tidx as usize] += count as f32;
        }
    }
}

pub fn em_optimize<P: EqClassPayload>(
    eqclasses: &HashMap<Vec<u32>, P, ahash::RandomState>,
    unique_evidence: &mut [bool],
    no_ambiguity: &mut [bool],
    init_type: EmInitType,
    num_alphas: usize,
    only_unique: bool,
    _log: &slog::Logger,
) -> Vec<f32> {
    let mut alphas_in: Vec<f32> = vec![0.0; num_alphas];
    let mut alphas_out: Vec<f32> = vec![0.0; num_alphas];

    for (labels, payload) in eqclasses {
        let count = payload.count();
        if labels.len() == 1 {
            let idx = labels.first().expect("can't extract labels");
            alphas_in[*idx as usize] += count as f32;
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
    let mut rng = rand::rng();
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
                *item = rng.random::<f32>() + 1e-5;
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
        let mut max_rel_diff = -f64::INFINITY;

        for index in 0..num_alphas {
            if alphas_out[index] > ALPHA_CHECK_CUTOFF {
                let diff = alphas_in[index] - alphas_out[index];
                let rel_diff = diff.abs();

                max_rel_diff = match rel_diff > max_rel_diff as f32 {
                    true => rel_diff as f64,
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

/// One cell's counts on a feature axis (genes, or probes at probe level),
/// held by a worker and reused for every cell it resolves.
///
/// `em_optimize` allocates, initializes and scans vectors over the whole axis
/// for every cell, and the caller then scans the result again to find the
/// nonzero entries, although a cell only touches the features its own
/// equivalence classes name. Here `values` (and, once the EM has run, `scratch`)
/// span the axis once per worker, and every per-cell pass walks only `support`:
/// the features this cell's classes name. Between cells, both vectors are zero
/// outside `support`. Every method that starts a cell (`reset`, `unique_counts`,
/// `em`, and the resolvers that fill a `SparseCellCounts`) forgets the previous
/// one first.
#[derive(Default)]
pub struct SparseCellCounts {
    values: Vec<f32>,
    scratch: Vec<f32>,
    in_support: Vec<bool>,
    support: Vec<u32>,
}

impl SparseCellCounts {
    /// Forgets the previous cell and sizes the axis to `num_features`.
    pub fn reset(&mut self, num_features: usize) {
        for &index in &self.support {
            let index = index as usize;
            self.values[index] = 0.0;
            self.in_support[index] = false;
        }
        self.support.clear();
        self.values.resize(num_features, 0.0);
        self.in_support.resize(num_features, false);
    }

    #[inline]
    fn mark(&mut self, index: usize) {
        if !self.in_support[index] {
            self.in_support[index] = true;
            self.support
                .push(u32::try_from(index).expect("feature index does not fit in u32"));
        }
    }

    /// Adds `count` to feature `index`.
    #[inline]
    pub fn add(&mut self, index: usize, count: f32) {
        self.mark(index);
        self.values[index] += count;
    }

    /// Forgets the previous cell and counts the classes that name a single
    /// feature, over an axis of `num_features`, as
    /// `em_optimize(.., only_unique = true)` does.
    pub fn unique_counts<P: EqClassPayload>(
        &mut self,
        eqclasses: &HashMap<Vec<u32>, P, ahash::RandomState>,
        num_features: usize,
    ) {
        self.reset(num_features);
        for (labels, payload) in eqclasses {
            if let [index] = labels[..] {
                self.add(index as usize, payload.count() as f32);
            }
        }
    }

    /// Forgets the previous cell and runs `em_optimize(.., only_unique = false)`
    /// over this cell's support, on an axis of `num_alphas` features.
    ///
    /// A feature outside the support appears in no class: `em_update` never
    /// reads or writes it, so the dense version zeroes it in the first
    /// iteration (there are always at least `MIN_ITER`) whatever its initial
    /// value. Within the support, every step is the dense version's, in the
    /// same order, so for the uniform and informative initializations the
    /// result is bit-identical. A random initialization draws one value per
    /// supported feature instead of one per feature.
    pub fn em<P: EqClassPayload>(
        &mut self,
        eqclasses: &HashMap<Vec<u32>, P, ahash::RandomState>,
        init_type: EmInitType,
        num_alphas: usize,
    ) {
        self.reset(num_alphas);
        // Only the EM needs the second vector; it stays zero between cells.
        self.scratch.resize(num_alphas, 0.0);
        for (labels, payload) in eqclasses {
            if let [index] = labels[..] {
                self.add(index as usize, payload.count() as f32);
            } else {
                for &index in labels {
                    self.mark(index as usize);
                }
            }
        }

        let Self {
            values: alphas_in,
            scratch: alphas_out,
            support,
            ..
        } = self;

        let mut rng = rand::rng();
        let uni_prior = 1.0 / (num_alphas as f32);
        for &index in support.iter() {
            let item = &mut alphas_in[index as usize];
            match init_type {
                EmInitType::Uniform => {
                    *item = uni_prior;
                }
                EmInitType::Informative => {
                    *item = (*item + 0.5) * 1e-3;
                }
                EmInitType::Random => {
                    *item = rng.random::<f32>() + 1e-5;
                }
            }
        }

        let mut it_num: u32 = 0;
        let mut converged: bool = true;
        while it_num < MIN_ITER || (it_num < MAX_ITER && !converged) {
            em_update(alphas_in, alphas_out, eqclasses);

            converged = true;
            for &index in support.iter() {
                let index = index as usize;
                if alphas_out[index] > ALPHA_CHECK_CUTOFF {
                    let rel_diff = (alphas_in[index] - alphas_out[index]).abs();
                    if rel_diff > REL_DIFF_TOLERANCE {
                        converged = false;
                    }
                }
                alphas_in[index] = alphas_out[index];
                alphas_out[index] = 0.0_f32;
            }

            it_num += 1;
        }

        for &index in support.iter() {
            let alpha = &mut alphas_in[index as usize];
            if *alpha < MIN_OUTPUT_ALPHA {
                *alpha = 0.0_f32;
            }
        }
    }

    /// Writes this cell's nonzero counts to `ind`/`vals` in increasing
    /// feature order, and returns their maximum and sum. These are exactly
    /// the entries, order and floating-point sum of a scan over the dense
    /// vector: the entries skipped here are all zero.
    pub fn write_nonzero(&mut self, ind: &mut Vec<usize>, vals: &mut Vec<f32>) -> (f32, f32) {
        self.support.sort_unstable();
        ind.clear();
        vals.clear();
        let mut max = 0.0f32;
        let mut sum = 0.0f32;
        for &index in &self.support {
            let c = self.values[index as usize];
            if c > max {
                max = c;
            }
            sum += c;
            if c > 0.0 {
                ind.push(index as usize);
                vals.push(c);
            }
        }
        (max, sum)
    }

    /// The counts as a dense vector over the axis `reset` last sized.
    pub fn to_dense(&self) -> Vec<f32> {
        let mut dense = vec![0.0f32; self.values.len()];
        for &index in &self.support {
            dense[index as usize] = self.values[index as usize];
        }
        dense
    }
}

#[allow(clippy::too_many_arguments)]
fn run_bootstrap_subset_with_scratch(
    eqclasses: &IndexedEqList,
    cell_data: &[(u32, u32)],
    num_alphas: u32,
    num_bootstraps: u32,
    _init_uniform: bool,
    summary_stat: bool,
    _log: &slog::Logger,
    scratch: &mut EmSubsetScratch,
) -> Vec<Vec<f32>> {
    // the population sample size
    let total_fragments: usize = cell_data.iter().map(|x| x.1 as usize).sum();
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

    // define the sampling weights for our multinomial
    let eq_counts: Vec<u32> = cell_data.iter().map(|x| x.1).collect();
    let mut dist = Multinomial::new(eq_counts, total_fragments).unwrap();

    // store bootstraps
    let mut bootstrap_counts = Vec::with_capacity(cell_data.len());

    // store bootstraps
    let num_output_bs = if summary_stat {
        2usize
    } else {
        num_bootstraps as usize
    };
    let mut bootstraps = Vec::with_capacity(num_output_bs);

    // A bootstrap changes equivalence-class counts, not the classes or their
    // labels. Compute the cell's possible support once and reuse that safe
    // superset even when a replicate assigns zero observations to every class
    // containing a particular gene.
    scratch.prepare_support(eqclasses, cell_data, num_alphas_us, None);

    // bootstrap loop starts
    // let mut old_resampled_counts = Vec::new();
    let mut rnd = rand::rng();
    for _bs_num in 0..num_bootstraps {
        // resample from multinomial
        let resampled_counts = dist.sample_u32(&mut rnd);
        for (idx, (eq_id, _orig_count)) in cell_data.iter().enumerate() {
            bootstrap_counts.push((*eq_id, resampled_counts[idx]));
        }

        let alphas = em_optimize_subset_impl(
            eqclasses,
            &bootstrap_counts[..], // indices into eqclasses relevant for this cell
            &mut unique_evidence,
            &mut no_ambiguity,
            EmInitType::Random,
            num_alphas_us,
            false, // only unique
            None,
            scratch,
            true,
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

pub fn run_bootstrap<P: EqClassPayload>(
    eqclasses: &HashMap<Vec<u32>, P, ahash::RandomState>,
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
    let mut scratch = EmSubsetScratch::new(gene_alpha.len());
    run_bootstrap_with_scratch(
        eqclasses,
        num_bootstraps,
        gene_alpha,
        _init_uniform,
        summary_stat,
        _log,
        &mut scratch,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_bootstrap_with_scratch<P: EqClassPayload>(
    eqclasses: &HashMap<Vec<u32>, P, ahash::RandomState>,
    num_bootstraps: u32,
    gene_alpha: &[f32],
    _init_uniform: bool,
    summary_stat: bool,
    _log: &slog::Logger,
    scratch: &mut EmSubsetScratch,
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
        .map(|(idx, (_labels, payload))| {
            let count = payload.count();
            (idx as u32, count)
        })
        .collect();

    // now that we have the `IndexedEqList` representation of this data, just
    // run that version of the bootstrap function and return the result.
    run_bootstrap_subset_with_scratch(
        &eql,
        &cell_data[..],
        gene_alpha.len() as u32,
        num_bootstraps,
        _init_uniform,
        summary_stat,
        _log,
        scratch,
    )
}

/*
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
        eq_counts.push(*count);
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
    let mut dist = Multinomial::new(eq_counts, total_fragments as usize).unwrap();

    // store bootstraps
    let mut bootstraps = Vec::new();

    // bootstrap loop starts
    // let mut old_resampled_counts = Vec::new();
    let mut rnd = rand::rng();
    for _bs_num in 0..num_bootstraps {
        // resample from multinomial
        let resampled_counts = dist.sample_u32(&mut rnd);

        for (eq_id, labels) in &eqclasses_serialize {
            eqclass_bootstrap.insert(labels.to_vec(), resampled_counts[*eq_id]);
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
            let mut max_rel_diff = -f64::INFINITY;

            for index in 0..gene_alpha.len() {
                if alphas_prime[index] > ALPHA_CHECK_CUTOFF {
                    let diff = alphas[index] - alphas_prime[index];
                    let rel_diff = diff.abs();

                    max_rel_diff = match rel_diff > max_rel_diff as f32 {
                        true => rel_diff,
                        false => max_rel_diff as f32,
                    } as f64;

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
*/

//em tailored for long read dataset
pub fn em_update_long_read<P: EqClassPayload>(
    alphas_in: &[f32],
    alphas_out: &mut [f32],
    eqclasses_prob: &HashMap<Vec<u32>, P, ahash::RandomState>,
) {
    // loop over all the eqclasses
    for (labels, payload) in eqclasses_prob {
        let count = payload.count();
        if labels.len() > 1 {
            let probs = payload.probs();
            let mut denominator: f32 = 0.0;
            let mut avg_probs: Vec<f32> = Vec::with_capacity(labels.len());
            for (tx_idx, label) in labels.iter().enumerate() {
                let sum_prob: f64 = (0..count as usize).map(|i| probs[i][tx_idx]).sum();
                let avg_prob: f32 = sum_prob as f32 / count as f32;
                avg_probs.push(avg_prob);
                denominator += alphas_in[*label as usize] * avg_prob;
            }

            if denominator > 0.0 {
                let inv_denominator = (count as f32) / denominator;
                for (tx_idx, label) in labels.iter().enumerate() {
                    let index = *label as usize;
                    let count = alphas_in[index] * avg_probs[tx_idx] * inv_denominator;
                    alphas_out[index] += count;
                }
            }
        } else {
            let tidx = labels.first().expect("can't extract labels");
            alphas_out[*tidx as usize] += count as f32;
        }
    }
}

pub fn em_optimize_long_read<P: EqClassPayload>(
    eqclasses_prob: &HashMap<Vec<u32>, P, ahash::RandomState>,
    unique_evidence: &mut [bool],
    no_ambiguity: &mut [bool],
    init_type: EmInitType,
    num_alphas: usize,
    only_unique: bool,
    _log: &slog::Logger,
) -> Vec<f32> {
    let mut alphas_in: Vec<f32> = vec![0.0; num_alphas];
    let mut alphas_out: Vec<f32> = vec![0.0; num_alphas];

    for (labels, payload) in eqclasses_prob {
        if labels.len() == 1 {
            let idx = labels.first().expect("can't extract labels");
            alphas_in[*idx as usize] += payload.count() as f32;
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
    let mut rng = rand::rng();
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
                *item = rng.random::<f32>() + 1e-5;
            }
        }
    }

    // TODO: is it even necessary?
    //alphas_in.iter_mut().for_each(|alpha| *alpha *= 1e-3);

    let mut it_num: u32 = 0;
    let mut converged: bool = true;
    while it_num < MIN_ITER || (it_num < MAX_ITER && !converged) {
        // perform one round of em update
        em_update_long_read(&alphas_in, &mut alphas_out, eqclasses_prob);

        converged = true;
        let mut max_rel_diff = -f64::INFINITY;

        for index in 0..num_alphas {
            if alphas_out[index] > ALPHA_CHECK_CUTOFF {
                let diff = alphas_in[index] - alphas_out[index];
                let rel_diff = diff.abs();

                max_rel_diff = match rel_diff > max_rel_diff as f32 {
                    true => rel_diff as f64,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inference_workspace_clears_counts_between_usa_cells() {
        let eql = eq_list(2, &[&[0], &[2], &[4], &[1], &[3], &[5], &[0, 4], &[1, 5]]);
        // Reuse after an EM cell, a unique-only early return, another EM,
        // and an empty cell. Old counts and USA sibling states must not leak.
        let cells: &[&[(u32, u32)]] = &[
            &[(0, 20), (1, 4), (2, 3), (6, 8)],
            &[(3, 7)],
            &[(7, 5)],
            &[],
        ];
        let mut workspace = SubsetInferenceCounts::default();
        for cell in cells {
            let expected =
                assert_sparse_matches_dense(&eql, cell, EmInitType::Informative, 6, Some((2, 4)));
            workspace.optimize(&eql, cell, 6, Some((2, 4)));
            assert_eq!(bits(&workspace.values), bits(&expected));
        }
    }

    fn eq_list(num_genes: usize, labels: &[&[u32]]) -> IndexedEqList {
        let mut eqclasses = IndexedEqList::new();
        eqclasses.num_genes = num_genes;
        eqclasses.eq_label_starts.push(0);
        for label in labels {
            eqclasses.add_label_vec(label);
        }
        eqclasses
    }

    #[allow(clippy::too_many_arguments)]
    fn dense_reference(
        eqclasses: &IndexedEqList,
        cell_data: &[(u32, u32)],
        unique_evidence: &mut [bool],
        no_ambiguity: &mut [bool],
        init_type: EmInitType,
        num_alphas: usize,
        only_unique: bool,
        usa_offsets: Option<(usize, usize)>,
    ) -> Vec<f32> {
        let mut alphas_in = vec![0.0; num_alphas];
        let mut alphas_out = vec![0.0; num_alphas];

        let mut needs_em = false;
        for (eq_id, count) in cell_data {
            let labels = eqclasses.refs_for_eqc(*eq_id);
            if labels.len() == 1 {
                let index = labels[0] as usize;
                alphas_in[index] += *count as f32;
                unique_evidence[index] = true;
            } else {
                for &index in labels {
                    no_ambiguity[index as usize] = false;
                }
                needs_em = true;
            }
        }

        if only_unique || !needs_em {
            return alphas_in;
        }

        let uniform_prior = 1.0 / num_alphas as f32;
        for alpha in &mut alphas_in {
            match init_type {
                EmInitType::Uniform => *alpha = uniform_prior,
                EmInitType::Informative => *alpha = (*alpha + 0.5) * 1e-3,
                EmInitType::Random => unreachable!("random output is not reproducible"),
            }
        }

        let mut iteration = 0;
        let mut converged = true;
        let mut last_round = false;
        while iteration < MIN_ITER || (iteration < MAX_ITER && !converged) || last_round {
            match usa_offsets {
                Some(offsets) => {
                    em_update_subset_usa(&alphas_in, &mut alphas_out, eqclasses, cell_data, offsets)
                }
                None => em_update_subset(&alphas_in, &mut alphas_out, eqclasses, cell_data),
            }

            converged = true;
            for index in 0..num_alphas {
                if alphas_out[index] > ALPHA_CHECK_CUTOFF
                    && (alphas_in[index] - alphas_out[index]).abs() > REL_DIFF_TOLERANCE
                {
                    converged = false;
                }
                alphas_in[index] = alphas_out[index];
                alphas_out[index] = 0.0;
            }
            iteration += 1;

            if last_round {
                break;
            }
            if iteration >= MIN_ITER && converged {
                for alpha in &mut alphas_in {
                    if *alpha < MIN_OUTPUT_ALPHA {
                        *alpha = 0.0;
                    }
                }
                last_round = true;
            }
        }

        for alpha in &mut alphas_in {
            if *alpha < MIN_OUTPUT_ALPHA {
                *alpha = 0.0;
            }
        }
        alphas_in
    }

    fn assert_sparse_matches_dense(
        eqclasses: &IndexedEqList,
        cell_data: &[(u32, u32)],
        init_type: EmInitType,
        num_alphas: usize,
        usa_offsets: Option<(usize, usize)>,
    ) -> Vec<f32> {
        let mut sparse_unique = vec![false; num_alphas];
        let mut sparse_no_ambiguity = vec![true; num_alphas];
        let sparse = em_optimize_subset(
            eqclasses,
            cell_data,
            &mut sparse_unique,
            &mut sparse_no_ambiguity,
            init_type,
            num_alphas,
            false,
            usa_offsets,
            &slog::Logger::root(slog::Discard, slog::o!()),
        );

        let mut dense_unique = vec![false; num_alphas];
        let mut dense_no_ambiguity = vec![true; num_alphas];
        let dense = dense_reference(
            eqclasses,
            cell_data,
            &mut dense_unique,
            &mut dense_no_ambiguity,
            init_type,
            num_alphas,
            false,
            usa_offsets,
        );

        assert_eq!(sparse, dense);
        assert_eq!(sparse_unique, dense_unique);
        assert_eq!(sparse_no_ambiguity, dense_no_ambiguity);
        sparse
    }

    #[test]
    fn sparse_matches_dense_for_empty_singleton_mixed_and_converged_cells() {
        let eqclasses = eq_list(8, &[&[0], &[1], &[0, 1], &[1, 2], &[2, 3, 4]]);
        let cases: &[&[(u32, u32)]] = &[&[], &[(0, 7)], &[(0, 20), (1, 4), (2, 8), (3, 1), (4, 2)]];
        for &init_type in &[EmInitType::Uniform, EmInitType::Informative] {
            for cell_data in cases {
                assert_sparse_matches_dense(&eqclasses, cell_data, init_type, 8, None);
            }
        }
    }

    #[test]
    fn sparse_matches_dense_for_all_usa_states_and_mixed_labels() {
        // Three genes occupy spliced [0,3), unspliced [3,6), and ambiguous
        // [6,9). Each case has a multi-label class so it exercises the EM.
        let eqclasses = eq_list(3, &[&[0, 1], &[3, 4], &[6, 7], &[0, 4, 8], &[2]]);
        let cases: &[&[(u32, u32)]] = &[
            &[(0, 5)],
            &[(1, 5)],
            &[(2, 5)],
            &[(0, 3), (1, 4), (2, 5), (3, 7), (4, 2)],
        ];
        for &init_type in &[EmInitType::Uniform, EmInitType::Informative] {
            for cell_data in cases {
                assert_sparse_matches_dense(&eqclasses, cell_data, init_type, 9, Some((3, 6)));
            }
        }
    }

    #[test]
    fn sparse_matches_dense_at_output_threshold() {
        let eqclasses = eq_list(6, &[&[0], &[0, 1], &[1, 2], &[2, 3, 4]]);
        let output = assert_sparse_matches_dense(
            &eqclasses,
            &[(0, 10_000), (1, 1), (2, 1), (3, 1)],
            EmInitType::Informative,
            6,
            None,
        );
        assert!(output.contains(&0.0));
    }

    #[test]
    fn random_initialization_uses_only_possible_cell_support() {
        let eqclasses = eq_list(8, &[&[1, 2], &[2, 3]]);
        let mut scratch = EmSubsetScratch::new(8);
        scratch.prepare_support(&eqclasses, &[(0, 4), (1, 0)], 8, None);
        assert_eq!(scratch.support, vec![1, 2, 3]);

        let mut unique = vec![false; 8];
        let mut no_ambiguity = vec![true; 8];
        let output = em_optimize_subset_impl(
            &eqclasses,
            &[(0, 4), (1, 0)],
            &mut unique,
            &mut no_ambiguity,
            EmInitType::Random,
            8,
            false,
            None,
            &mut scratch,
            true,
        );
        assert!(output[1] + output[2] > 0.0);
        assert_eq!(output[0], 0.0);
        assert_eq!(output[3..], [0.0; 5]);
    }

    #[test]
    fn bootstrap_reuses_possible_support_and_never_emits_outside_it() {
        let eqclasses = eq_list(8, &[&[1, 2], &[2, 3], &[3]]);
        let cell_data = [(0, 5), (1, 3), (2, 2)];
        let mut scratch = EmSubsetScratch::new(8);
        let bootstraps = run_bootstrap_subset_with_scratch(
            &eqclasses,
            &cell_data,
            8,
            8,
            false,
            false,
            &slog::Logger::root(slog::Discard, slog::o!()),
            &mut scratch,
        );

        assert_eq!(scratch.support, vec![1, 2, 3]);
        assert_eq!(bootstraps.len(), 8);
        for output in bootstraps {
            assert!(output[1] + output[2] + output[3] > 0.0);
            assert_eq!(output[0], 0.0);
            assert_eq!(output[4..], [0.0; 4]);
        }
    }

    #[test]
    fn scratch_support_deduplicates_clears_and_reuses_membership() {
        let eqclasses = eq_list(5, &[&[0, 1, 1], &[1, 2], &[3, 4]]);
        let mut scratch = EmSubsetScratch::new(5);
        let membership_capacity = scratch.support_membership.capacity();

        scratch.prepare_support(&eqclasses, &[(0, 2), (1, 3)], 5, None);
        assert_eq!(scratch.support, vec![0, 1, 2]);
        scratch.prepare_support(&eqclasses, &[(2, 1)], 5, None);
        assert_eq!(scratch.support, vec![3, 4]);
        assert!(!scratch.support_membership[0]);
        assert!(!scratch.support_membership[1]);
        assert!(!scratch.support_membership[2]);
        assert_eq!(scratch.support_membership.capacity(), membership_capacity);
    }

    fn cell_eqc(
        classes: &[(&[u32], u32)],
    ) -> HashMap<Vec<u32>, crate::utils::BasicEqClassPayload, ahash::RandomState> {
        let mut eqc = HashMap::with_hasher(ahash::RandomState::with_seeds(2, 7, 1, 8));
        for &(labels, ct) in classes {
            eqc.insert(labels.to_vec(), crate::utils::BasicEqClassPayload { ct });
        }
        eqc
    }

    fn bits(v: &[f32]) -> Vec<u32> {
        v.iter().map(|x| x.to_bits()).collect()
    }

    /// The dense path this replaces: `em_optimize`, then a scan of its output.
    fn dense_cell(
        eqc: &HashMap<Vec<u32>, crate::utils::BasicEqClassPayload, ahash::RandomState>,
        init_type: EmInitType,
        num_alphas: usize,
        only_unique: bool,
    ) -> (Vec<f32>, Vec<usize>, Vec<f32>, f32, f32) {
        let mut unique_evidence = vec![false; num_alphas];
        let mut no_ambiguity = vec![false; num_alphas];
        let counts = em_optimize(
            eqc,
            &mut unique_evidence,
            &mut no_ambiguity,
            init_type,
            num_alphas,
            only_unique,
            &slog::Logger::root(slog::Discard, slog::o!()),
        );
        let (mut ind, mut vals, mut max, mut sum) = (Vec::new(), Vec::new(), 0.0f32, 0.0f32);
        for (i, &c) in counts.iter().enumerate() {
            max = if c > max { c } else { max };
            sum += c;
            if c > 0.0 {
                ind.push(i);
                vals.push(c);
            }
        }
        (counts, ind, vals, max, sum)
    }

    fn assert_cell_counts_match_dense(
        sparse: &mut SparseCellCounts,
        classes: &[(&[u32], u32)],
        init_type: EmInitType,
        num_alphas: usize,
        only_unique: bool,
    ) {
        let eqc = cell_eqc(classes);
        let (dense, d_ind, d_vals, d_max, d_sum) =
            dense_cell(&eqc, init_type, num_alphas, only_unique);
        if only_unique {
            sparse.unique_counts(&eqc, num_alphas);
        } else {
            sparse.em(&eqc, init_type, num_alphas);
        }
        assert_eq!(bits(&sparse.to_dense()), bits(&dense), "{classes:?}");
        let (mut ind, mut vals) = (vec![usize::MAX], vec![-1.0]);
        let (max, sum) = sparse.write_nonzero(&mut ind, &mut vals);
        assert_eq!(ind, d_ind, "{classes:?}");
        assert_eq!(bits(&vals), bits(&d_vals), "{classes:?}");
        assert_eq!(max.to_bits(), d_max.to_bits(), "{classes:?}");
        assert_eq!(sum.to_bits(), d_sum.to_bits(), "{classes:?}");
    }

    #[test]
    fn cell_counts_match_em_optimize_and_the_dense_scan() {
        let cells: &[&[(&[u32], u32)]] = &[
            &[],
            &[(&[3], 5)],
            &[
                (&[0], 20),
                (&[1], 4),
                (&[0, 1], 8),
                (&[1, 2], 1),
                (&[2, 3, 4], 2),
            ],
            // only ambiguous classes: every feature starts from its prior
            &[(&[5, 9], 3), (&[9, 11], 7)],
            // a feature whose EM share falls under MIN_OUTPUT_ALPHA
            &[(&[2], 1000), (&[2, 7], 1)],
            // slow convergence: a long chain of overlapping classes
            &[
                (&[0, 1], 9),
                (&[1, 2], 9),
                (&[2, 3], 9),
                (&[3, 4], 9),
                (&[4], 1),
                (&[0], 30),
            ],
        ];
        let mut sparse = SparseCellCounts::default();
        for &init_type in &[EmInitType::Uniform, EmInitType::Informative] {
            for only_unique in [true, false] {
                for cell in cells {
                    assert_cell_counts_match_dense(&mut sparse, cell, init_type, 12, only_unique);
                }
            }
        }
    }

    #[test]
    fn cell_counts_match_dense_on_generated_cells_with_one_reused_scratch() {
        // A worker reuses one scratch for every cell: each result must not
        // depend on what the previous cells touched.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = move |n: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % n
        };
        let num_alphas = 200;
        let mut sparse = SparseCellCounts::default();
        for cell_num in 0..500 {
            let num_classes = next(40) as usize;
            let mut owned: Vec<(Vec<u32>, u32)> = Vec::new();
            for _ in 0..num_classes {
                let len = 1 + next(4) as usize;
                let mut labels: Vec<u32> =
                    (0..len).map(|_| next(num_alphas as u64) as u32).collect();
                labels.sort_unstable();
                labels.dedup();
                owned.push((labels, 1 + next(50) as u32));
            }
            let classes: Vec<(&[u32], u32)> = owned.iter().map(|(l, c)| (&l[..], *c)).collect();
            let init_type = if cell_num % 2 == 0 {
                EmInitType::Uniform
            } else {
                EmInitType::Informative
            };
            assert_cell_counts_match_dense(
                &mut sparse,
                &classes,
                init_type,
                num_alphas,
                cell_num % 3 == 0,
            );
        }
    }

    #[test]
    fn usa_unique_counts_match_extract_counts_with_one_reused_scratch() {
        // USA labels are `gene << 1 | unspliced`; the axis is genes x {S, U, A}.
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move |n: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % n
        };
        let num_genes = 60u64;
        let num_counts = 3 * num_genes as usize;
        let mut sparse = SparseCellCounts::default();
        for _ in 0..500 {
            let mut owned: Vec<(Vec<u32>, u32)> = Vec::new();
            for _ in 0..next(30) {
                let len = 1 + next(4) as usize;
                let mut labels: Vec<u32> = (0..len)
                    .map(|_| {
                        let g = next(num_genes) as u32;
                        // often both statuses of one gene, the A column's case
                        (g << 1) | next(2) as u32
                    })
                    .collect();
                labels.sort_unstable();
                labels.dedup();
                owned.push((labels, 1 + next(50) as u32));
            }
            let classes: Vec<(&[u32], u32)> = owned.iter().map(|(l, c)| (&l[..], *c)).collect();
            let eqc = cell_eqc(&classes);
            let dense = crate::utils::extract_counts(&eqc, num_counts);
            crate::utils::extract_counts_sparse(&eqc, num_counts, &mut sparse);
            assert_eq!(bits(&sparse.to_dense()), bits(&dense), "{classes:?}");
        }
    }
}
