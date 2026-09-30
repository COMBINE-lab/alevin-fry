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
use crate::quant::{ResolutionStrategy, SplicedAmbiguityModel, TrueUmiOracleMode, EditDistanceMode, EMMode, EndModel, EndDistribution};

use std::path::PathBuf;

#[derive(TypedBuilder, Debug, Serialize)]
//#[builder(name = "QuantOptsBuilder")]
pub struct QuantOpts<'a, 'b, 'c, 'd, 'e, 'f, 'g> {
    pub input_dir: &'a PathBuf,
    pub tg_map: &'b PathBuf,
    pub output_dir: &'c PathBuf,
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
