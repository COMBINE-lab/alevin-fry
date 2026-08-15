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

use crate::barcode_correction::{
    BarcodeNeighborhood, BarcodeResolution, Confidence, CorrectionSpec,
};
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum SampleCorrectionMode {
    /// Exact match only — no error correction
    Exact,
    /// Accept a non-exact barcode only when it has one canonical target.
    Unique,
    /// Resolve competing canonical samples using frozen exact counts.
    Frequency,
    /// Deprecated compatibility spelling for Unique plus the historical
    /// substitution-or-shift neighbourhood.
    OneEdit,
}

impl std::fmt::Display for SampleCorrectionMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Exact => f.write_str("exact"),
            Self::Unique => f.write_str("unique"),
            Self::Frequency => f.write_str("frequency"),
            Self::OneEdit => f.write_str("1-edit"),
        }
    }
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
    /// How unmatched cell barcodes are corrected to retained barcodes when an
    /// unfiltered external permit list is used.
    #[builder(default)]
    pub cell_bc_correction: CellBarcodeCorrectionStrategy,
    /// Explicit cell neighbourhood. `None` retains protocol/filter-specific
    /// historical defaults.
    #[builder(default)]
    pub cell_bc_neighborhood: Option<BarcodeNeighborhood>,
    /// Sample neighbourhood for the new Unique/Frequency modes.
    #[builder(default = BarcodeNeighborhood::HammingOne)]
    pub sample_bc_neighborhood: BarcodeNeighborhood,
    #[builder(default = Confidence::RNA)]
    pub cell_bc_confidence: Confidence,
    #[builder(default = Confidence::RNA)]
    pub sample_bc_confidence: Confidence,
    /// Total memory available to deferred sample-Frequency buffers.
    #[builder(default = 512_u64 << 20)]
    pub memory_limit: u64,
    /// Directory for correction spools; defaults to `output_dir`.
    #[builder(default)]
    pub tmp_dir: Option<PathBuf>,
}

impl GenPermitListOpts<'_, '_, '_, '_, '_> {
    pub fn resolved_cell_neighborhood(&self, multi_barcode: bool) -> BarcodeNeighborhood {
        self.cell_bc_neighborhood.unwrap_or({
            if multi_barcode || matches!(self.fmeth, CellFilterMethod::UnfilteredExternalList(_, _))
            {
                BarcodeNeighborhood::HammingOne
            } else {
                BarcodeNeighborhood::SubstitutionOrShiftOne
            }
        })
    }

    pub fn cell_correction_spec(&self, barcode_len: u8, multi_barcode: bool) -> CorrectionSpec {
        let resolution = match self.cell_bc_correction {
            CellBarcodeCorrectionStrategy::Unique => BarcodeResolution::Unique,
            CellBarcodeCorrectionStrategy::Frequency => BarcodeResolution::Frequency {
                confidence: self.cell_bc_confidence,
                pseudocount: 1,
            },
        };
        CorrectionSpec {
            barcode_len,
            neighborhood: self.resolved_cell_neighborhood(multi_barcode),
            resolution,
        }
    }

    pub fn sample_correction_spec(&self, barcode_len: u8) -> Option<CorrectionSpec> {
        let (neighborhood, resolution) = match self.sample_correction_mode {
            SampleCorrectionMode::Exact => return None,
            SampleCorrectionMode::Unique => {
                (self.sample_bc_neighborhood, BarcodeResolution::Unique)
            }
            SampleCorrectionMode::Frequency => (
                self.sample_bc_neighborhood,
                BarcodeResolution::Frequency {
                    confidence: self.sample_bc_confidence,
                    pseudocount: 1,
                },
            ),
            SampleCorrectionMode::OneEdit => (
                BarcodeNeighborhood::SubstitutionOrShiftOne,
                BarcodeResolution::Unique,
            ),
        };
        Some(CorrectionSpec {
            barcode_len,
            neighborhood,
            resolution,
        })
    }

    /// Enforce the documented practical lower bound inside the library API,
    /// not just in the CLI parser.
    pub fn effective_memory_limit(&self) -> u64 {
        const MINIMUM: u64 = 256_u64 << 20;
        if self.memory_limit < MINIMUM {
            slog::warn!(
                self.log,
                "Requested a {} byte GPL memory limit, but correction is not designed for less than 256 MiB; proceeding with 256 MiB.",
                self.memory_limit
            );
            MINIMUM
        } else {
            self.memory_limit
        }
    }
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
