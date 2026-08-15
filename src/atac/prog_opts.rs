use crate::atac::cellfilter::CellFilterMethod;
use crate::barcode_correction::{
    BarcodeNeighborhood, BarcodeResolution, Confidence, CorrectionSpec,
};
use crate::prog_opts::CellBarcodeCorrectionStrategy;
use serde::Serialize;
use slog;
use std::path::PathBuf;
use typed_builder::TypedBuilder;

#[derive(TypedBuilder, Debug, Serialize)]
pub struct GenPermitListOpts<'a, 'b, 'c, 'd, 'e> {
    pub input_dir: &'a PathBuf,
    pub output_dir: &'b PathBuf,
    pub fmeth: CellFilterMethod,
    pub threads: usize,
    pub rc: bool,
    pub cmdline: &'c str,
    pub version: &'d str,
    #[serde(skip_serializing)]
    pub log: &'e slog::Logger,
    #[builder(default)]
    pub cell_bc_correction: CellBarcodeCorrectionStrategy,
    #[builder(default = BarcodeNeighborhood::HammingOne)]
    pub cell_bc_neighborhood: BarcodeNeighborhood,
    #[builder(default = Confidence::ATAC)]
    pub cell_bc_confidence: Confidence,
}

impl GenPermitListOpts<'_, '_, '_, '_, '_> {
    pub fn correction_spec(&self, barcode_len: u8) -> CorrectionSpec {
        CorrectionSpec {
            barcode_len,
            neighborhood: self.cell_bc_neighborhood,
            resolution: match self.cell_bc_correction {
                CellBarcodeCorrectionStrategy::Unique => BarcodeResolution::Unique,
                CellBarcodeCorrectionStrategy::Frequency => BarcodeResolution::Frequency {
                    confidence: self.cell_bc_confidence,
                    pseudocount: 1,
                },
            },
        }
    }
}

#[derive(TypedBuilder, Debug, Serialize)]
pub struct DeduplicateOpts<'a, 'b, 'c, 'd> {
    pub input_dir: &'a PathBuf,
    pub num_threads: u32,
    pub rev: bool,
    pub cmdline: &'b str,
    pub version: &'c str,
    #[serde(skip_serializing)]
    pub log: &'d slog::Logger,
}
