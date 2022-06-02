use crate::quant::{ResolutionStrategy, SplicedAmbiguityModel};
//use derive_builder::Builder;
use typed_builder::TypedBuilder;

#[derive(TypedBuilder, Debug)]
//#[builder(name = "QuantOptsBuilder")]
pub struct QuantOpts<'a, 'b, 'c, 'd> {
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
    pub filter_list: Option<&'a str>,
    pub cmdline: &'b str,
    pub version: &'c str,
    pub log: &'d slog::Logger,
}
