use crate::atac::cellfilter::CellFilterMethod;
use serde::Serialize;
use slog;
use std::path::PathBuf;
use typed_builder::TypedBuilder;

#[derive(TypedBuilder, Debug, Serialize)]
pub struct GenPermitListOpts<'a, 'b, 'c, 'd, 'e> {
    pub input_dir: &'a PathBuf,
    pub output_dir: &'b PathBuf,
    pub fmeth: CellFilterMethod,
    pub rc: bool,
    pub cmdline: &'c str,
    pub version: &'d str,
    #[serde(skip_serializing)]
    pub log: &'e slog::Logger,
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
