//use derive_builder::Builder;
use crate::cellfilter::CellFilterMethod;
use serde::Serialize;
use slog;
use std::path::PathBuf;
use typed_builder::TypedBuilder;
// use crate::cellfilter::CellFilterMethod;
// use crate::quant::{ResolutionStrategy, SplicedAmbiguityModel};

// #[derive(TypedBuilder, Debug, Serialize)]
//#[builder(name = "QuantOptsBuilder")]
// pub struct QuantOpts<'a, 'b, 'c, 'd, 'e, 'f, 'g> {
//     pub input_dir: &'a PathBuf,
//     pub tg_map: &'b PathBuf,
//     pub output_dir: &'c PathBuf,
//     pub num_threads: u32,
//     pub num_bootstraps: u32,
//     pub init_uniform: bool,
//     pub summary_stat: bool,
//     pub dump_eq: bool,
//     pub use_mtx: bool,
//     pub resolution: ResolutionStrategy,
//     pub pug_exact_umi: bool,
//     pub sa_model: SplicedAmbiguityModel,
//     pub small_thresh: usize,
//     pub large_graph_thresh: usize,
//     pub filter_list: Option<&'d PathBuf>,
//     pub cmdline: &'e str,
//     pub version: &'f str,
//     #[serde(skip_serializing)]
//     pub log: &'g slog::Logger,
// }

// #[derive(TypedBuilder, Debug, Serialize)]
// pub struct GenPermitListOpts<'a, 'b, 'c, 'd, 'e> {
//     pub input_dir: &'a PathBuf,
//     pub output_dir: &'b PathBuf,
//     pub fmeth: CellFilterMethod,
//     pub expected_ori: Strand,
//     pub velo_mode: bool,
//     pub cmdline: &'c str,
//     pub version: &'d str,
//     #[serde(skip_serializing)]
//     pub log: &'e slog::Logger,
// }
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
