//! scATAC-seq processing.
//!
//! The supported pipeline is [`cellfilter`] (`generate-permit-list`) followed
//! by [`sort`], which bins records by genomic position, deduplicates them, and
//! writes the BED. That is what `simpleaf atac process` drives.
//!
//! [`collate`] and [`deduplicate`] are a second, barcode-oriented route to a
//! similar end result. They are **provisional**: nothing in the supported
//! pipeline calls them, and their subcommands are hidden from `--help`. See
//! each module's own comment before building on them.

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
