//use derive_builder::Builder;
use bio_types::strand::Strand;
use slog;
use typed_builder::TypedBuilder;

use crate::cellfilter::CellFilterMethod;
use crate::quant::{ResolutionStrategy, SplicedAmbiguityModel};

#[derive(TypedBuilder, Debug)]
//#[builder(name = "QuantOptsBuilder")]
pub struct QuantOpts<'b, 'c, 'd> {
    pub input_dir: String,
    pub tg_map: String,
    pub output_dir: String,
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
    pub large_graph_thresh: usize,
    pub filter_list: Option<String>,
    pub cmdline: &'b str,
    pub version: &'c str,
    pub log: &'d slog::Logger,
}

#[derive(TypedBuilder, Debug)]
pub struct GenPermitListOpts<'a, 'b, 'c> {
    pub input_dir: String,
    pub output_dir: String,
    pub fmeth: CellFilterMethod,
    pub expected_ori: Strand,
    pub velo_mode: bool,
    pub cmdline: &'a str,
    pub version: &'b str,
    pub log: &'c slog::Logger,
}
