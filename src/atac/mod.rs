// These mirror the top-level scRNA modules (`cellfilter`, `collate`, `quant`,
// ...), which lib.rs exports publicly. Keeping the scATAC equivalents private
// meant the whole path could only be reached through `run`, i.e. by driving the
// CLI — which is part of why it went untested for so long.
pub mod cellfilter;
pub mod collate;
pub mod deduplicate;
pub mod prog_opts;
pub mod run;
pub mod sort;
mod utils;
