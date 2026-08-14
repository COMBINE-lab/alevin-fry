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

/// Resolves a raw cell barcode that is one edit away from multiple retained
/// cell barcodes in a filtered RNA permit list.
///
/// This policy does not affect exact matches, unfiltered RNA workflows, sample
/// barcodes, or ATAC barcode correction.
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq, Default)]
pub enum BarcodeCollisionPolicy {
    /// Do not correct a barcode when more than one retained barcode is a
    /// one-edit neighbor. This is the conservative default.
    #[default]
    DropAmbiguous,
    /// Correct to the retained neighbor with the largest observed count. Ties
    /// are resolved by choosing the lowest packed-barcode value.
    PreferMostFrequent,
    /// Correct to the lowest packed-barcode neighbor.
    AcceptFirstNeighbor,
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
    /// Policy for ambiguous one-edit cell-barcode corrections in filtered RNA
    /// permit-list workflows. This is serialized into the run metadata.
    #[builder(default)]
    pub cell_barcode_collision_policy: BarcodeCollisionPolicy,
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
    /// How unmatched cell barcodes are corrected to retained barcodes when an
    /// unfiltered external permit list is used.
    #[builder(default)]
    pub cell_bc_correction: CellBarcodeCorrectionStrategy,
}

/// How to resolve an unmatched barcode that is one substitution away from more
/// than one retained barcode.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum CellBarcodeCorrectionStrategy {
    /// Correct only when there is exactly one such neighbour, and drop the
    /// barcode otherwise. This is what alevin-fry has always done.
    #[default]
    Unique,
    /// Rank retained neighbours by their Laplace-smoothed pre-correction
    /// frequencies and accept the best one when it carries at least 97.5% of
    /// the total weight. This is Cell Ranger's posterior decision rule under
    /// uniform base qualities, applied to alevin-fry's retained candidates.
    Frequency,
}

impl std::fmt::Display for CellBarcodeCorrectionStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Unique => write!(f, "unique"),
            Self::Frequency => write!(f, "frequency"),
        }
    }
}
