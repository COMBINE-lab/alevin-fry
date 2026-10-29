/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//use derive_builder::Builder;
use bio_types::strand::Strand;
use serde::Serialize;
use slog;
use typed_builder::TypedBuilder;

use crate::cellfilter::CellFilterMethod;
use crate::quant::{
    EMMode, EditDistanceMode, EndDistribution, EndModel, ResolutionStrategy, SplicedAmbiguityModel,
    TrueUmiOracleMode,
};
use crate::truncation::{
    CoverageSource, TruncationAbundance, TruncationEndKernels, TruncationIntervalMode,
    TruncationLengthLaw, TruncationParamScope, TruncationPositionalLaw, TruncationPreset,
};

use std::path::PathBuf;

#[derive(TypedBuilder, Debug, Serialize)]
//#[builder(name = "QuantOptsBuilder")]
pub struct QuantOpts<'a, 'b, 'c, 'd, 'e, 'f, 'g> {
    pub input_dir: &'a PathBuf,
    pub tg_map: &'b PathBuf,
    pub output_dir: &'c PathBuf,
    pub junction_model: bool,
    pub junction_gtf: Option<&'a PathBuf>,
    pub num_threads: u32,
    pub num_bootstraps: u32,
    pub init_uniform: bool,
    pub summary_stat: bool,
    pub dump_eq: bool,
    pub use_mtx: bool,
    pub resolution: ResolutionStrategy,
    pub pug_exact_umi: bool,
    pub sa_model: SplicedAmbiguityModel,
    pub small_thresh: usize,
    pub lambda_size: f64,
    pub tau_delta: f64,
    pub large_graph_thresh: usize,
    pub filter_list: Option<&'d PathBuf>,
    pub cmdline: &'e str,
    pub version: &'f str,
    pub true_umi_mode: TrueUmiOracleMode,
    pub true_txp_mode: Option<bool>,
    pub true_umi_txp_file: Option<PathBuf>,
    pub ed_distance_model: EditDistanceMode,
    pub ed_upper_bound: usize,
    pub em_prob_model: EMMode,
    pub end_model: EndModel,
    pub end_distribution: EndDistribution,
    pub end_threshold: f64,
    pub end_scale: f64,
    pub coverage_model: bool,
    pub coverage_dedup_model: bool,
    pub truncation_model: bool,
    pub truncation_preset: Option<TruncationPreset>,
    pub truncation_discover_groups: bool,
    pub truncation_discover_genes: Option<PathBuf>,
    pub truncation_interval_mode: TruncationIntervalMode,
    pub truncation_params: TruncationParamScope,
    pub truncation_length_law: TruncationLengthLaw,
    pub truncation_end_kernels: TruncationEndKernels,
    pub truncation_abundance: TruncationAbundance,
    pub truncation_shrink: f64,
    pub truncation_shrink_auto: bool,
    pub truncation_shrink_per_gene: bool,
    pub truncation_positional_power: f64,
    pub truncation_calibrate_power: bool,
    pub truncation_min_share: f64,
    pub truncation_sample_below: f64,
    pub truncation_global_cells: usize,
    pub truncation_groups: Option<PathBuf>,
    pub truncation_genes: Option<PathBuf>,
    pub truncation_dump: Option<PathBuf>,
    pub truncation_oracle: Option<PathBuf>,
    pub truncation_group_shrink: f64,
    pub truncation_group_shrink_auto: bool,
    pub truncation_group_shrink_per_gene: bool,
    pub truncation_positional_law: TruncationPositionalLaw,
    pub coverage_source: CoverageSource,
    pub junction_likelihood: bool,
    pub coverage_bin_width: u32,
    pub coverage_growth_rate: f64,
    #[serde(skip_serializing)]
    pub log: &'g slog::Logger,
}

#[derive(TypedBuilder, Debug, Serialize)]
pub struct GenPermitListOpts<'a, 'b, 'c, 'd, 'e> {
    pub input_dir: &'a PathBuf,
    pub output_dir: &'b PathBuf,
    pub fmeth: CellFilterMethod,
    pub expected_ori: Strand,
    pub velo_mode: bool,
    pub threads: usize,
    pub cmdline: &'c str,
    pub version: &'d str,
    #[serde(skip_serializing)]
    pub log: &'e slog::Logger,
}
