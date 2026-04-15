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
use crate::quant::{ResolutionStrategy, SplicedAmbiguityModel};

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
    pub resolution: ResolutionStrategy,
    pub pug_exact_umi: bool,
    pub sa_model: SplicedAmbiguityModel,
    pub small_thresh: usize,
    pub large_graph_thresh: usize,
    pub filter_list: Option<&'d PathBuf>,
    pub cmdline: &'e str,
    pub version: &'f str,
    #[serde(skip_serializing)]
    pub log: &'g slog::Logger,
}

/// Correction mode for sample barcodes in multi-barcode protocols.
#[derive(Debug, Clone, Serialize)]
pub enum SampleCorrectionMode {
    /// Exact match only — no error correction
    Exact,
    /// Allow single-edit correction using BarcodeLookupMap
    OneEdit,
}

/// Orientation of the sample/probe barcodes in the whitelist relative to
/// how they appear in the read. `Forward` means the whitelist is already in
/// read-orientation; `Reverse` means each whitelist entry must be
/// reverse-complemented before lookup (e.g., 10x Flex v2, where the sample
/// BC on R1 downstream of the TTGCTAGGACCG anchor is the RC of the
/// vendor-published list).
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq, Default)]
pub enum SampleBarcodeOri {
    #[default]
    Forward,
    Reverse,
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
    /// Path to known sample barcode list (one per line). When present,
    /// triggers multi-barcode mode (e.g., 10x Flex).
    #[builder(default)]
    pub sample_bc_list: Option<PathBuf>,
    /// Path to sample name mapping file (TSV: barcode\tname).
    #[builder(default)]
    pub sample_names: Option<PathBuf>,
    /// Correction mode for sample barcodes.
    #[builder(default = SampleCorrectionMode::Exact)]
    pub sample_correction_mode: SampleCorrectionMode,
    /// Orientation of sample barcodes in the whitelist relative to the read.
    #[builder(default = SampleBarcodeOri::Forward)]
    pub sample_bc_ori: SampleBarcodeOri,
}
