/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use anyhow::Context;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};

#[allow(unused_imports)]
use slog::{crit, info, warn};

use needletail::bitkmer::*;
use num_format::{Locale, ToFormattedString};
use serde::Serialize;
use serde_json::json;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::str::FromStr;
use std::string::ToString;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use libradicl::collation::CollationManifest;
use libradicl::header::RadPrelude;
use libradicl::rad_types::TagMap;
use libradicl::record::{
    AlevinFryReadRecord, AlevinFryReadRecordWithPosition, CollatableMappedRecord,
    CollatableRecordHeader, ConvertiblePrimitiveInteger, KnownSize, MappedRecord,
    MultiBarcodeReadRecord, RecordContext, ScLongReadRecord, UmiTaggedRecord,
};

use std::fmt;
//use std::ptr;

use flate2::Compression;
use flate2::write::GzEncoder;

use crate::em::{
    EmInitType, em_optimize, em_optimize_long_read, em_optimize_subset, run_bootstrap,
};
use crate::eq_class::{EqMap, EqMapType, IndexedEqList};
use crate::prog_opts::QuantOpts;
use crate::pugutils;
use crate::utils as afutils;
use crate::utils::{
    BasicEqClassPayload, EqClassPayload, KnownRecordType, LongReadEqClassPayload,
    OptionalAlignmentExtras,
};

#[derive(PartialEq, Eq, Debug, Default, Clone, Copy, Serialize)]
pub enum SplicedAmbiguityModel {
    PreferAmbiguity,
    #[default]
    WinnerTakeAll,
}

// Implement the trait
impl FromStr for SplicedAmbiguityModel {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "prefer-ambig" => Ok(SplicedAmbiguityModel::PreferAmbiguity),
            "winner-take-all" => Ok(SplicedAmbiguityModel::WinnerTakeAll),
            _ => Err("no match"),
        }
    }
}

#[derive(PartialEq, Eq, Debug, Default, Clone, Copy, Serialize)]
pub enum ResolutionStrategy {
    Trivial,
    #[default]
    CellRangerLike,
    CellRangerLikeEm,
    ParsimonyEm,
    Parsimony,
    ParsimonyGeneEm,
    ParsimonyGene,
}

impl fmt::Display for ResolutionStrategy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

// Implement the trait
impl FromStr for ResolutionStrategy {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "trivial" => Ok(ResolutionStrategy::Trivial),
            "cr-like" => Ok(ResolutionStrategy::CellRangerLike),
            "cr-like-em" => Ok(ResolutionStrategy::CellRangerLikeEm),
            "parsimony-em" => Ok(ResolutionStrategy::ParsimonyEm),
            "parsimony" => Ok(ResolutionStrategy::Parsimony),
            "parsimony-gene" => Ok(ResolutionStrategy::ParsimonyGene),
            "parsimony-gene-em" => Ok(ResolutionStrategy::ParsimonyGeneEm),
            _ => Err("no match"),
        }
    }
}

/// Bootstrap output configuration.
///
/// Bootstrap summary statistics (mean and variance across replicates) are
/// written as sparse MTX files: `bootstraps_mean.mtx` and `bootstraps_var.mtx`,
/// with the same dimensions as the main count matrix (cells × genes).
///
/// When `--summary-stat` is NOT set, the full per-replicate counts are computed
/// but currently only the summary stats are written. Full per-replicate output
/// will be supported via HDF5/AnnData layers in a future release, which
/// provides efficient storage for the 3D (cells × genes × replicates) array
/// via chunked, compressed datasets in `.h5ad` format.
struct BootstrapHelper {
    num_bootstraps: u32,
    summary_stat: bool,
    /// Thread-local bootstrap mean triplets: (row, col, val)
    mean_triplets: Vec<(usize, usize, f32)>,
    /// Thread-local bootstrap variance triplets: (row, col, val)
    var_triplets: Vec<(usize, usize, f32)>,
}

impl BootstrapHelper {
    fn new(
        _output_path: &std::path::Path,
        num_bootstraps: u32,
        summary_stat: bool,
    ) -> BootstrapHelper {
        if num_bootstraps > 0 && !summary_stat {
            eprintln!(
                "NOTE: Full per-replicate bootstrap output is not yet supported in MTX format. \
                 Summary statistics (mean, variance) will be written instead. \
                 Full replicate output will be available in a future release via AnnData/h5ad."
            );
        }
        BootstrapHelper {
            num_bootstraps,
            summary_stat,
            mean_triplets: Vec::new(),
            var_triplets: Vec::new(),
        }
    }

    /// Record bootstrap results for one cell. Computes mean and variance
    /// across replicates and stores nonzero entries as sparse triplets.
    fn record_cell(&mut self, row_index: usize, bootstraps: &[Vec<f32>]) {
        if bootstraps.is_empty() {
            return;
        }

        let (mean_vec, var_vec) = if self.summary_stat && bootstraps.len() == 2 {
            // run_bootstrap already returned [mean, var]
            (&bootstraps[0], &bootstraps[1])
        } else {
            // Shouldn't happen in current flow, but handle gracefully
            return;
        };

        for (col, &val) in mean_vec.iter().enumerate() {
            if val != 0.0 {
                self.mean_triplets.push((row_index, col, val));
            }
        }
        for (col, &val) in var_vec.iter().enumerate() {
            if val != 0.0 {
                self.var_triplets.push((row_index, col, val));
            }
        }
    }

    /// Compute mean and variance from full bootstrap replicates for one cell.
    fn record_cell_from_replicates(&mut self, row_index: usize, bootstraps: &[Vec<f32>]) {
        if bootstraps.is_empty() {
            return;
        }
        let n = bootstraps.len() as f32;
        let num_genes = bootstraps[0].len();

        for col in 0..num_genes {
            let mean: f32 = bootstraps.iter().map(|b| b[col]).sum::<f32>() / n;
            if mean != 0.0 {
                self.mean_triplets.push((row_index, col, mean));
                let var: f32 = bootstraps
                    .iter()
                    .map(|b| {
                        let d = b[col] - mean;
                        d * d
                    })
                    .sum::<f32>()
                    / (n - 1.0).max(1.0);
                if var != 0.0 {
                    self.var_triplets.push((row_index, col, var));
                }
            }
        }
    }
}

struct QuantOutputInfo {
    barcode_file: BufWriter<fs::File>,
    feature_file: BufWriter<fs::File>,
    row_index: usize,
}

struct EqcMap {
    // the *global* gene-level equivalence class map
    global_eqc: HashMap<Vec<u32>, u64, ahash::RandomState>,
    // the list of equivalence classes (and corresponding umi count)
    // that occurs in each cell.
    cell_level_count: Vec<(u64, u32)>,
    // a vector of tuples that contains pairs of the form
    // (row offset, number of equivalence classes in this cell)
    cell_offset: Vec<(usize, usize)>,
}

fn write_eqc_counts(
    eqid_map_lock: &Arc<Mutex<EqcMap>>,
    num_genes: usize,
    usa_mode: bool,
    output_path: &std::path::Path,
    log: &slog::Logger,
) -> anyhow::Result<bool> {
    let eqmap_deref = eqid_map_lock.lock();
    let geqmap = eqmap_deref.unwrap();
    let num_eqclasses = geqmap.global_eqc.len();

    info!(
        log,
        "Writing gene-level equivalence class with {:?} classes",
        geqmap.global_eqc.len()
    );

    // the sparse matrix that will hold the equivalence class counts
    let mut eqmat = sprs::TriMatI::<f32, u32>::with_capacity(
        (geqmap.cell_offset.len(), num_eqclasses), // cells x eq-classes
        geqmap.cell_level_count.len(),             // num non-zero entries
    );

    // fill in the matrix
    let mut global_offset = 0usize;
    for (row_index, num_cell_eqs) in geqmap.cell_offset.iter() {
        let slice = global_offset..(global_offset + num_cell_eqs);
        for (eqid, umi_count) in geqmap.cell_level_count[slice].iter() {
            eqmat.add_triplet(*row_index, *eqid as usize, *umi_count as f32);
        }
        global_offset += num_cell_eqs;
    }

    // and write it to file.
    let mtx_path = output_path.join("geqc_counts.mtx");
    sprs::io::write_matrix_market(mtx_path, &eqmat).context("could not write geqc_counts.mtx")?;

    // write the sets of genes that define each eqc
    let gn_eq_path = output_path.join("gene_eqclass.txt.gz");
    let mut gn_eq_writer = BufWriter::new(GzEncoder::new(
        fs::File::create(gn_eq_path).unwrap(),
        Compression::default(),
    ));

    // number of genes
    gn_eq_writer
        .write_all(format!("{}\n", num_genes).as_bytes())
        .context("could not write to gene_eqclass.txt.gz")?;

    // number of classes
    gn_eq_writer
        .write_all(format!("{}\n", num_eqclasses).as_bytes())
        .context("could not write to gene_eqclass.txt.gz")?;

    // each line describes a class in terms of
    // the tab-separated tokens
    // g_1 g_2 ... g_k eqid

    if usa_mode {
        // if we are running in USA mode, then:
        // spliced (even) IDs get divided by 2
        // and odd IDs get divided by 2 and added to the
        // unspliced offset.

        // offset for unspliced gene ids
        let unspliced_offset = (num_genes / 3) as u32;
        // offset for ambiguous gene ids
        let ambig_offset = 2 * unspliced_offset;
        // to hold the gene labels as we write them.
        let mut gl;

        // if we are running the *standard* mode, then the gene_id
        // mapping is unaltered
        for (gene_list, eqid) in geqmap.global_eqc.iter() {
            // strategy for peeking ahead as needed derived from
            // https://sts10.github.io/2020/10/06/peeking-the-pivot.html
            let mut peekable_arr = gene_list.iter().peekable();
            // get the current gene label
            while let Some(cg) = peekable_arr.next() {
                // get the next gene label in the eq class
                if let Some(ng) = peekable_arr.peek() {
                    // if the gene label belongs to the same gene
                    // then it must be splicing ambiguous (because exact
                    // duplicate IDs can't occur in eq class labels).
                    if afutils::same_gene(*cg, **ng, true) {
                        gl = (cg >> 1) + ambig_offset;
                        gn_eq_writer
                            .write_all(format!("{}\t", gl).as_bytes())
                            .expect("could not write to gene_eqclass.txt.gz");
                        // we covered the next element here, so skip it in the
                        // next iteration.
                        peekable_arr.next();
                        continue;
                    }
                }
                // either the next element does *not* belong to the same
                // gene, or there is no next element.  In either case, deal
                // with this gene label individually.
                if afutils::is_spliced(*cg) {
                    gl = cg >> 1;
                } else {
                    gl = (cg >> 1) + unspliced_offset;
                }
                gn_eq_writer
                    .write_all(format!("{}\t", gl).as_bytes())
                    .context("could not write to gene_eqclass.txt.gz")?;
            }
            gn_eq_writer
                .write_all(format!("{}\n", eqid).as_bytes())
                .context("could not write to gene_eqclass.txt.gz")?;
        }
    } else {
        // if we are running the *standard* mode, then the gene_id
        // mapping is unaltered
        for (gene_list, eqid) in geqmap.global_eqc.iter() {
            for g in gene_list.iter() {
                gn_eq_writer
                    .write_all(format!("{}\t", g).as_bytes())
                    .context("could not write to gene_eqclass.txt.gz")?;
            }
            gn_eq_writer
                .write_all(format!("{}\n", eqid).as_bytes())
                .context("could not write to gene_eqclass.txt.gz")?;
        }
    }
    Ok(true)
}

// TODO: see if we'd rather pass an structure
// with these options
pub fn quantify(quant_opts: QuantOpts) -> anyhow::Result<()> {
    let parent = std::path::Path::new(quant_opts.input_dir);
    let log = quant_opts.log;

    // read the collate metadata
    let collate_md_file =
        File::open(parent.join("collate.json")).context("could not open the collate.json file.")?;
    let collate_md: serde_json::Value = serde_json::from_reader(&collate_md_file)?;

    // is the collated RAD file compressed?
    let compressed_input = collate_md["compressed_output"]
        .as_bool()
        .context("could not read compressed_output field from collate metadata.")?;

    if compressed_input {
        let i_file =
            File::open(parent.join("map.collated.rad.sz")).context("run collate before quant")?;
        let br = BufReader::new(snap::read::FrameDecoder::new(&i_file));

        info!(
            log,
            "quantifying from compressed, collated RAD file {:?}", i_file
        );

        do_quantify_dispatch(br, quant_opts)
    } else {
        let i_file =
            File::open(parent.join("map.collated.rad")).context("run collate before quant")?;
        let br = BufReader::new(&i_file);

        info!(
            log,
            "quantifying from uncompressed, collated RAD file {:?}", i_file
        );

        do_quantify_dispatch(br, quant_opts)
    }
}

struct WorkerConfig {
    resolution: ResolutionStrategy,
    usa_mode: bool,
    usa_offsets: Option<(usize, usize)>,
    em_init_type: EmInitType,
    large_graph_thresh: usize,
    pug_exact_umi: bool,
    sa_model: SplicedAmbiguityModel,
    num_bootstraps: u32,
    init_uniform: bool,
    summary_stat: bool,
    dump_eq: bool,
    num_genes: usize,
    num_rows: usize,
    barcode_len: u16,
}

struct WorkerSharedState<R: MappedRecord> {
    in_q: Arc<crossbeam_queue::ArrayQueue<libradicl::readers::MetaChunk<R>>>,
    is_done: Arc<std::sync::atomic::AtomicBool>,
    tid_to_gid: Arc<Vec<u32>>,
    cells_remaining: Arc<AtomicUsize>,
    bcout: Arc<Mutex<QuantOutputInfo>>,
    eqid_map_lock: Arc<Mutex<EqcMap>>,
    alt_res_cells: Arc<Mutex<Vec<u64>>>,
    empty_resolved_cells: Arc<Mutex<Vec<u64>>>,
    unmapped_count: Arc<libradicl::unmapped::CollatedUnmappedCounts>,
    mmrate: Arc<Mutex<Vec<f64>>>,
    /// Sample names indexed by sample index. None for single-barcode.
    sample_names: Option<Arc<Vec<String>>>,
    /// Extracts the sample index from a record. None for single-barcode types.
    sample_idx_extractor: Option<Arc<dyn Fn(&R) -> usize + Send + Sync>>,
}

/// Threshold (in number of records) below which a cell uses the fast path
/// that avoids HashMap-based equivalence class construction entirely.
const SMALL_CELL_FAST_THRESHOLD: usize = 100;

/// Sparse fast path for quantifying small cells without any HashMap or dense
/// vector overhead.
///
/// Implements the same cr-like (winner-take-all) UMI resolution as
/// `resolve_num_molecules_crlike_from_vec` in pugutils, but outputs sparse
/// `(slot_index, umi)` pairs in sorted order so the caller can extract counts
/// via run-length counting — avoiding both O(num_genes) dense vector zeroing
/// and the O(num_genes) scan to find non-zero entries.
///
/// For each read, all transcripts are projected to gene IDs (with dedup).
/// For each UMI, gene votes are tallied across all reads. If there is a
/// unique gene with the maximum vote count, that UMI is attributed to that
/// gene. Tied UMIs (multi-gene) are discarded in non-USA mode.
///
/// In USA mode, tied UMIs undergo splicing-aware resolution matching the
/// logic in `extract_counts`: same-gene S+U → ambiguous slot, cross-gene
/// prefer-spliced, etc. The output indices are USA slot offsets (spliced,
/// unspliced, or ambiguous) rather than raw gene IDs.
fn quantify_small_cell_sparse<B, R>(
    chunk: &mut libradicl::chunk::Chunk<R>,
    tid_to_gid: &[u32],
    usa_mode: bool,
    num_rows: usize, // num_genes in non-USA, 3*num_genes in USA
    // Output: sorted (slot_index, umi) pairs for resolved UMIs
    gene_umi_buf: &mut Vec<(u32, u64)>,
    // Reusable scratch buffer for (umi, gene, count) triplets
    umi_gene_triplets: &mut Vec<(u64, u32, u32)>,
) where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize + UmiTaggedRecord + 'static,
{
    gene_umi_buf.clear();
    umi_gene_triplets.clear();

    // Phase 1: For each read, project transcript refs to gene IDs (deduped)
    // and emit (umi, gene, 1) triplets — same as get_num_molecules_cell_ranger_like_small.
    for read in &chunk.reads {
        let umi = read.umi();
        let refs = read.refs();

        // Collect unique gene IDs for this read's reference set.
        // Use a small inline approach: for most reads the gene set is 1-3 genes.
        let first_ref = refs.first().copied().unwrap_or(0) as usize;
        if refs.is_empty() {
            continue;
        }
        let first_gid = tid_to_gid[first_ref];
        let mut is_single_gene = true;
        // Check if all refs map to the same gene (common case)
        for &ref_id in &refs[1..] {
            if tid_to_gid[ref_id as usize] != first_gid {
                is_single_gene = false;
                break;
            }
        }
        if is_single_gene {
            umi_gene_triplets.push((umi, first_gid, 1));
        } else {
            // Multi-gene read: collect unique gene IDs
            // For < 100 reads, this branch is rare and the allocation is tiny
            let mut gset: smallvec::SmallVec<[u32; 4]> =
                refs.iter().map(|&tid| tid_to_gid[tid as usize]).collect();
            gset.sort_unstable();
            gset.dedup();
            for &g in &gset {
                umi_gene_triplets.push((umi, g, 1));
            }
        }
    }

    if umi_gene_triplets.is_empty() {
        return;
    }

    // USA mode slot offsets
    let unspliced_offset = if usa_mode { num_rows / 3 } else { 0 };
    let ambig_offset = if usa_mode { 2 * unspliced_offset } else { 0 };

    // Map a raw USA gene ID to its output slot index.
    // In non-USA mode, the gene ID is used directly.
    let to_slot = |gid: u32| -> u32 {
        if !usa_mode {
            gid
        } else if afutils::is_spliced(gid) {
            gid >> 1
        } else {
            (unspliced_offset as u32) + (gid >> 1)
        }
    };

    // Phase 2: cr-like resolution — same algorithm as
    // resolve_num_molecules_crlike_from_vec in pugutils.
    // Sort by (umi, gene, count).
    umi_gene_triplets.sort_unstable();

    let mut curr_umi = umi_gene_triplets[0].0;
    let mut curr_gn = umi_gene_triplets[0].1;
    let mut max_count = 0u32;
    let mut count_aggr = 0u32;
    // Track the full set of best genes (needed for USA splicing resolution)
    let mut best_genes: smallvec::SmallVec<[u32; 4]> = smallvec::smallvec![curr_gn];

    // Commit a resolved UMI to gene_umi_buf. For single-winner UMIs,
    // maps to slot and pushes. For multi-gene ties in USA mode, applies
    // splicing-aware resolution matching extract_counts in utils.rs.
    let commit_umi = |best: &smallvec::SmallVec<[u32; 4]>, umi: u64, buf: &mut Vec<(u32, u64)>| {
        match best.len() {
            0 => {}
            1 => {
                buf.push((to_slot(best[0]), umi));
            }
            _ if !usa_mode => {
                // Non-USA: multi-gene tie → discard
            }
            2 => {
                // USA: same logic as extract_counts len==2
                let (g1, g2) = (best[0], best[1]);
                if afutils::same_gene(g1, g2, true) {
                    // S+U of same gene → ambiguous slot
                    buf.push(((ambig_offset as u32) + (g1 >> 1), umi));
                } else {
                    // Different genes: prefer spliced
                    match (afutils::is_spliced(g1), afutils::is_spliced(g2)) {
                        (true, false) => buf.push(((g1 >> 1), umi)),
                        (false, true) => buf.push(((g2 >> 1), umi)),
                        _ => {} // both spliced or both unspliced → discard
                    }
                }
            }
            n if n <= 10 => {
                // USA: same logic as extract_counts len==3..10
                // Find spliced genes
                let mut spliced_iter = best.iter().filter(|&&x| afutils::is_spliced(x));
                if let Some(&first_spliced) = spliced_iter.next() {
                    if spliced_iter.next().is_some() {
                        // 2+ spliced genes → gene-ambiguous, discard
                    } else {
                        // Exactly 1 spliced gene. Check if its unspliced
                        // counterpart is also in the set.
                        let has_unspliced_partner = best.iter().any(|&g| {
                            g != first_spliced && afutils::same_gene(first_spliced, g, true)
                        });
                        if has_unspliced_partner {
                            buf.push(((ambig_offset as u32) + (first_spliced >> 1), umi));
                        } else {
                            buf.push(((first_spliced >> 1), umi));
                        }
                    }
                }
                // No spliced genes at all → discard
            }
            _ => {} // >10 labels → discard
        }
    };

    for idx in 0..umi_gene_triplets.len() {
        let (umi, gn, ct) = umi_gene_triplets[idx];

        if umi != curr_umi {
            // Commit previous UMI
            commit_umi(&best_genes, curr_umi, gene_umi_buf);

            // Reset for new UMI
            curr_umi = umi;
            curr_gn = gn;
            count_aggr = ct;
            max_count = ct;
            best_genes.clear();
            best_genes.push(gn);
        } else {
            // Same UMI
            if gn == curr_gn {
                count_aggr += ct;
            } else {
                count_aggr = ct;
                curr_gn = gn;
            }

            use std::cmp::Ordering;
            match count_aggr.cmp(&max_count) {
                Ordering::Greater => {
                    max_count = count_aggr;
                    match &best_genes[..] {
                        [x] if *x == gn => {}
                        _ => {
                            best_genes.clear();
                            best_genes.push(gn);
                        }
                    }
                }
                Ordering::Equal => {
                    best_genes.push(gn);
                }
                Ordering::Less => {}
            }
        }

        // Commit last UMI
        if idx == umi_gene_triplets.len() - 1 {
            commit_umi(&best_genes, curr_umi, gene_umi_buf);
        }
    }

    // Sort output by (slot, umi) for the caller's run-length counting
    gene_umi_buf.sort_unstable();
}

fn run_worker_thread<B, R, P>(
    _worker_num: usize,
    config: WorkerConfig,
    shared: WorkerSharedState<R>,
    log: slog::Logger,
    num_eq_targets: u32,
    eq_map_type: EqMapType,
) -> (usize, Vec<(usize, usize, f32)>, BootstrapHelper)
where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord
        + CollatableMappedRecord<B>
        + KnownSize
        + UmiTaggedRecord
        + OptionalAlignmentExtras
        + 'static,
    <R as MappedRecord>::ParsingContext: RecordContext,
    <R as MappedRecord>::ParsingContext: Clone,
    <R as MappedRecord>::ParsingContext: Send,
    P: EqClassPayload,
{
    // these can be created once and cleared after processing
    // each cell.
    let mut unique_evidence = vec![false; config.num_rows];
    let mut no_ambiguity = vec![false; config.num_rows];
    let mut eq_map = EqMap::new(num_eq_targets, eq_map_type, P::HAS_PROBS);
    let mut expressed_vec = Vec::<f32>::with_capacity(config.num_genes);
    let mut expressed_ind = Vec::<usize>::with_capacity(config.num_genes);
    // Thread-local triplet buffer for MTX output. Accumulating triplets locally
    // avoids holding the global mutex during add_triplet, which was the primary
    // bottleneck serializing all workers.
    let mut local_triplets: Vec<(usize, usize, f32)> = Vec::new();
    // Thread-local bootstrap helper for accumulating summary stat triplets
    let mut boot_helper = BootstrapHelper::new(
        std::path::Path::new(""),
        config.num_bootstraps,
        config.summary_stat,
    );
    // Reusable buffers for the small-cell sparse fast path
    let mut gene_umi_buf: Vec<(u32, u64)> = Vec::new();
    let mut umi_gene_triplets: Vec<(u64, u32, u32)> = Vec::new();

    // the variable we will use to bind the *cell-specific* gene-level
    // equivalence class table.
    // Make gene-level eqclasses.
    // This is a map of gene ids to the count of
    // _de-duplicated_ reads observed for that set of genes.
    // For every gene set (label) of length 1, these are gene
    // unique reads.  Standard scRNA-seq counting results
    // can be obtained by simply discarding all equivalence
    // classes of size greater than 1, and probabilistic results
    // will attempt to resolve gene multi-mapping reads by
    // running an EM algorithm.
    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut gene_eqc: HashMap<Vec<u32>, P, ahash::RandomState> = HashMap::with_hasher(s);

    // If we are operating in USA-mode with an EM capable resolution
    // method, we'll use (re-use) these variables to hold the USA-mode
    // equivalence class information.
    let mut idx_eq_list = IndexedEqList::new();
    let mut eq_id_count = Vec::<(u32, u32)>::new();

    let mut local_nrec = 0usize;
    // pop MetaChunks from the work queue until everything is
    // processed
    while !shared.is_done.load(Ordering::SeqCst) || !shared.in_q.is_empty() {
        while let Some(meta_chunk) = shared.in_q.pop() {
            let first_cell_in_chunk = meta_chunk.first_chunk_index;
            for (cn, mut c) in meta_chunk.iter().enumerate() {
                shared.cells_remaining.fetch_sub(1, Ordering::SeqCst);
                let cell_num = first_cell_in_chunk + cn;

                let nbytes = c.nbytes;
                let nrec = c.nrec;
                local_nrec += nrec as usize;

                if c.reads.is_empty() {
                    warn!(
                        log,
                        "Discovered empty chunk; should not happen! cell_num = {}, nbytes = {}, nrec = {}",
                        cell_num,
                        nbytes,
                        nrec
                    );
                }

                // TODO: Clean up the expect() and merge with the check above
                // the expect shouldn't happen, but the message is redundant with
                // the above.  Plus, this would panic if it actually occurred.
                let first_rec = c.reads.first().expect("chunk with no reads");
                let bc = first_rec.collate_key();
                // For multi-barcode records, extract the sample index
                // from the first record (stored in barcodes[0] by the
                // scatter phase) for direct Vec lookup of the sample name.
                let sample_idx_from_rec: Option<usize> = shared
                    .sample_idx_extractor
                    .as_ref()
                    .map(|ext| ext(first_rec));

                // The structures we'll need to hold our output for this
                // cell.
                let mut counts: Vec<f32>;
                let mut alt_resolution = false;
                let mut max_umi = 0.0f32;
                let mut sum_umi = 0.0f32;
                let mut num_expr: u32 = 0;

                let mut bootstraps: Vec<Vec<f32>> = Vec::new();

                // Fast path for very small cells: bypass all HashMap-based
                // equivalence class machinery. Directly compute gene counts
                // via sorted Vec dedup, and build the sparse output
                // (expressed_vec/expressed_ind) directly — avoiding both the
                // O(num_genes) dense counts vector allocation AND the
                // O(num_genes) scan to extract non-zero entries.
                let used_fast_path;
                if c.reads.len() < SMALL_CELL_FAST_THRESHOLD {
                    used_fast_path = true;
                    quantify_small_cell_sparse::<B, R>(
                        &mut c,
                        &shared.tid_to_gid,
                        config.usa_mode,
                        config.num_rows,
                        &mut gene_umi_buf,
                        &mut umi_gene_triplets,
                    );
                    // gene_umi_buf is now sorted and deduped (gene_id, umi) pairs.
                    // Build expressed_vec/expressed_ind by run-length counting genes.
                    expressed_vec.clear();
                    expressed_ind.clear();
                    let mut sum_umi_local = 0.0f32;
                    let mut max_umi_local = 0.0f32;
                    if !gene_umi_buf.is_empty() {
                        let mut cur_gene = gene_umi_buf[0].0;
                        let mut cur_count = 0u32;
                        for &(gid, _) in gene_umi_buf.iter() {
                            if gid == cur_gene {
                                cur_count += 1;
                            } else {
                                // emit previous gene
                                let c = cur_count as f32;
                                expressed_ind.push(cur_gene as usize);
                                expressed_vec.push(c);
                                sum_umi_local += c;
                                if c > max_umi_local {
                                    max_umi_local = c;
                                }
                                cur_gene = gid;
                                cur_count = 1;
                            }
                        }
                        // emit last gene
                        let c = cur_count as f32;
                        expressed_ind.push(cur_gene as usize);
                        expressed_vec.push(c);
                        sum_umi_local += c;
                        if c > max_umi_local {
                            max_umi_local = c;
                        }
                    }
                    // Set the output variables that the rest of the function expects
                    max_umi = max_umi_local;
                    sum_umi = sum_umi_local;
                    num_expr = expressed_vec.len() as u32;
                    // counts is only needed for EDS output; build lazily below
                    counts = Vec::new();
                } else
                // Original path for larger cells
                {
                    used_fast_path = false;
                    let non_trivial = true; // all cells here have >= SMALL_CELL_FAST_THRESHOLD reads
                    if non_trivial {
                        // TODO: some testing was done, but see if there is a better way to set this value.
                        let small_cell = c.reads.len() <= 250;

                        // TODO: Is there an easy / clean way to have similar
                        // optimized code paths for other resolution methods?

                        match config.resolution {
                            ResolutionStrategy::CellRangerLike
                            | ResolutionStrategy::CellRangerLikeEm => {
                                if small_cell {
                                    pugutils::get_num_molecules_cell_ranger_like_small::<B, R, P>(
                                        &mut c,
                                        &shared.tid_to_gid,
                                        config.num_genes,
                                        &mut gene_eqc,
                                        config.sa_model,
                                        &log,
                                    );
                                } else {
                                    eq_map.init_from_chunk::<R>(&mut c);
                                    pugutils::get_num_molecules_cell_ranger_like(
                                        &eq_map,
                                        &shared.tid_to_gid,
                                        config.num_genes,
                                        &mut gene_eqc,
                                        config.sa_model,
                                        &log,
                                    );
                                    eq_map.clear();
                                }
                                let only_unique =
                                    config.resolution == ResolutionStrategy::CellRangerLike;

                                // NOTE: This configuration seems overly complicated
                                // see if we can simplify it.
                                match (config.usa_mode, only_unique) {
                                    (true, true) => {
                                        // USA mode, only gene-unqique reads
                                        counts =
                                            afutils::extract_counts(&gene_eqc, config.num_rows);
                                    }
                                    (true, false) => {
                                        // USA mode, use EM
                                        afutils::extract_usa_eqmap(
                                            &gene_eqc,
                                            config.num_rows,
                                            &mut idx_eq_list,
                                            &mut eq_id_count,
                                        );
                                        counts = em_optimize_subset(
                                            &idx_eq_list,
                                            &eq_id_count,
                                            &mut unique_evidence,
                                            &mut no_ambiguity,
                                            config.em_init_type,
                                            config.num_rows,
                                            only_unique,
                                            config.usa_offsets,
                                            &log,
                                        );
                                    }
                                    (false, _) => {
                                        // not USA-mode
                                        counts = em_optimize(
                                            &gene_eqc,
                                            &mut unique_evidence,
                                            &mut no_ambiguity,
                                            config.em_init_type,
                                            config.num_genes,
                                            only_unique,
                                            &log,
                                        );
                                    }
                                }
                            }
                            ResolutionStrategy::Trivial => {
                                eq_map.init_from_chunk(&mut c);
                                let ct = pugutils::get_num_molecules_trivial_discard_all_ambig(
                                    &eq_map,
                                    &shared.tid_to_gid,
                                    config.num_genes,
                                    &log,
                                );
                                counts = ct.0;
                                shared.mmrate.lock().unwrap()[cell_num] = ct.1;
                                eq_map.clear();
                            }
                            ResolutionStrategy::Parsimony
                            | ResolutionStrategy::ParsimonyEm
                            | ResolutionStrategy::ParsimonyGene
                            | ResolutionStrategy::ParsimonyGeneEm => {
                                if (config.resolution == ResolutionStrategy::ParsimonyGene)
                                    || (config.resolution == ResolutionStrategy::ParsimonyGeneEm)
                                {
                                    eq_map.init_from_chunk_gene_level(&mut c, &shared.tid_to_gid);
                                } else {
                                    //eprintln!("before the init from chunk");
                                    eq_map.init_from_chunk(&mut c);
                                    //eprintln!("after the init from chunk");
                                }

                                let g =
                                    pugutils::extract_graph(&eq_map, config.pug_exact_umi, &log);
                                // for the PUG resolution algorithm, set the hasher
                                // that will be used based on the cell barcode.
                                let s = ahash::RandomState::with_seeds(bc.into(), 7u64, 1u64, 8u64);
                                let pug_stats = pugutils::get_num_molecules::<P>(
                                    &g,
                                    &eq_map,
                                    &shared.tid_to_gid,
                                    &mut gene_eqc,
                                    &s,
                                    config.large_graph_thresh,
                                    &log,
                                );
                                alt_resolution = pug_stats.used_alternative_strategy; // alt_res;
                                eq_map.clear();

                                let only_unique = (config.resolution
                                    == ResolutionStrategy::Parsimony)
                                    || (config.resolution == ResolutionStrategy::ParsimonyGene);

                                // NOTE: This configuration seems overly complicated
                                // see if we can simplify it.
                                match (config.usa_mode, only_unique, P::HAS_PROBS) {
                                    (true, true, _) => {
                                        // USA mode, only gene-unqique reads
                                        counts =
                                            afutils::extract_counts(&gene_eqc, config.num_rows);
                                    }
                                    (true, false, _) => {
                                        // USA mode, use EM
                                        afutils::extract_usa_eqmap(
                                            &gene_eqc,
                                            config.num_rows,
                                            &mut idx_eq_list,
                                            &mut eq_id_count,
                                        );
                                        counts = em_optimize_subset(
                                            &idx_eq_list,
                                            &eq_id_count,
                                            &mut unique_evidence,
                                            &mut no_ambiguity,
                                            config.em_init_type,
                                            config.num_rows,
                                            only_unique,
                                            config.usa_offsets,
                                            &log,
                                        );
                                    }
                                    (false, _, false) => {
                                        // not USA-mode
                                        counts = em_optimize(
                                            &gene_eqc,
                                            &mut unique_evidence,
                                            &mut no_ambiguity,
                                            config.em_init_type,
                                            config.num_genes,
                                            only_unique,
                                            &log,
                                        );
                                    }
                                    (false, _, true) => {
                                        // not USA-mode
                                        counts = em_optimize_long_read(
                                            &gene_eqc,
                                            &mut unique_evidence,
                                            &mut no_ambiguity,
                                            config.em_init_type,
                                            config.num_genes,
                                            only_unique,
                                            &log,
                                        );
                                    }
                                }
                            }
                        }

                        if config.num_bootstraps > 0 {
                            bootstraps = run_bootstrap(
                                &gene_eqc,
                                config.num_bootstraps,
                                &counts,
                                config.init_uniform,
                                config.summary_stat,
                                &log,
                            );
                        }

                        // clear our local variables
                        // eq_map.clear();

                        // fill requires >= 1.50.0
                        unique_evidence.fill(false);
                        no_ambiguity.fill(false);

                        // done clearing
                    } else {
                        // very small number of reads, avoid data structure
                        // overhead and resolve looking at the actual records
                        pugutils::get_num_molecules_cell_ranger_like_small(
                            &mut c,
                            &shared.tid_to_gid,
                            config.num_genes,
                            &mut gene_eqc,
                            config.sa_model,
                            &log,
                        );
                        // USA-mode
                        if config.usa_mode {
                            // here, just like for non-USA mode,
                            // we substitute EM with uniform allocation in
                            // this special case
                            match config.resolution {
                                ResolutionStrategy::CellRangerLike
                                | ResolutionStrategy::Parsimony
                                | ResolutionStrategy::ParsimonyGene => {
                                    counts = afutils::extract_counts(&gene_eqc, config.num_rows);
                                }
                                ResolutionStrategy::CellRangerLikeEm
                                | ResolutionStrategy::ParsimonyEm
                                | ResolutionStrategy::ParsimonyGeneEm => {
                                    counts = afutils::extract_counts_mm_uniform(
                                        &gene_eqc,
                                        config.num_rows,
                                    );
                                }
                                _ => {
                                    counts = vec![0f32; config.num_genes];
                                    warn!(
                                        log,
                                        "Should not reach here, only cr-like, cr-like-em, parsimony(-gene) and parsimony(-gene)-em are supported in USA-mode."
                                    );
                                }
                            }
                        } else {
                            // non USA-mode
                            counts = vec![0f32; config.num_genes];
                            for (k, payload) in gene_eqc.iter() {
                                let v = payload.count();
                                if k.len() == 1 {
                                    counts[*k.first().unwrap() as usize] += v as f32;
                                } else {
                                    match config.resolution {
                                        ResolutionStrategy::CellRangerLikeEm
                                        | ResolutionStrategy::ParsimonyEm => {
                                            let contrib = 1.0 / (k.len() as f32);
                                            for g in k.iter() {
                                                counts[*g as usize] += contrib;
                                            }
                                        }
                                        _ => {
                                            // otherwise discard gene multimappers
                                        }
                                    }
                                }
                            }
                        }
                        // if the user requested bootstraps
                        // NOTE: we check that the specified resolution method
                        // is conceptually compatible with bootstrapping before
                        // invoking `quant`, so we don't bother checking that
                        // here.
                        if config.num_bootstraps > 0 {
                            // TODO: should issue a warning here,
                            // bootstrapping doesn't make sense for
                            // unfiltered data.
                            if config.summary_stat {
                                // sample mean = quant
                                bootstraps.push(counts.clone());
                                // sample var = 0
                                bootstraps.push(vec![0f32; config.num_genes]);
                            } else {
                                // no variation
                                for _ in 0..config.num_bootstraps {
                                    bootstraps.push(counts.clone());
                                }
                            }
                        } // if the user requested bootstraps
                    } // end of else branch for trivial size cells
                } // end of else block for non-fast-path cells

                if alt_resolution {
                    shared.alt_res_cells.lock().unwrap().push(cell_num as u64);
                }

                //
                // featuresStream << "\t" << numRawReads
                //   << "\t" << numMappedReads
                //
                // For fast-path cells, max_umi/sum_umi/num_expr/expressed_*
                // are already computed sparsely above. For non-fast-path
                // cells, extract them from the dense counts vector.
                if !used_fast_path {
                    let mut max_umi_local = 0.0f32;
                    let mut sum_umi_local = 0.0f32;
                    expressed_vec.clear();
                    expressed_ind.clear();

                    for (gn, c) in counts.iter().enumerate() {
                        max_umi_local = if *c > max_umi_local {
                            *c
                        } else {
                            max_umi_local
                        };
                        sum_umi_local += *c;
                        if *c > 0.0 {
                            num_expr += 1;
                            expressed_vec.push(*c);
                            expressed_ind.push(gn);
                        }
                    }
                    max_umi = max_umi_local;
                    sum_umi = sum_umi_local;
                }

                if num_expr == 0 {
                    shared
                        .empty_resolved_cells
                        .lock()
                        .unwrap()
                        .push(cell_num as u64);
                }

                let num_mapped = nrec;
                let dedup_rate = sum_umi / num_mapped as f32;

                let bcint = bc.into();
                let num_unmapped = shared.unmapped_count.get_single(bcint);

                let mapping_rate = num_mapped as f32 / (num_mapped + num_unmapped) as f32;

                // mean of the "expressed" genes
                let mean_expr = sum_umi / num_expr as f32;
                // number of genes with expression > expressed mean
                let num_genes_over_mean = expressed_vec
                    .iter()
                    .fold(0u32, |acc, x| if x > &mean_expr { acc + 1u32 } else { acc });
                // expressed mean / max expression
                let mean_by_max = mean_expr / max_umi;

                let row_index: usize; // the index for this row (cell)
                {
                    // writing the files
                    let bc_mer: BitKmer = (bc.into(), config.barcode_len as u8);

                    // Scope the lock to minimize hold time — triplet accumulation
                    // happens after the lock is released.
                    {
                        let writer_deref = shared.bcout.lock();
                        let writer = &mut *writer_deref.unwrap();

                        // get the row index and then increment it
                        row_index = writer.row_index;
                        writer.row_index += 1;

                        // write to barcode file
                        let bc_bytes = &bitmer_to_bytes(bc_mer)[..];
                        let bc_str = unsafe { std::str::from_utf8_unchecked(bc_bytes) };

                        // For multi-barcode data, prefix with sample name.
                        // The sample index was read from the first record's
                        // barcodes[0] (written by the scatter phase), so
                        // assignment is correct regardless of chunk order.
                        let sample_name = sample_idx_from_rec.and_then(|si| {
                            shared.sample_names.as_ref()?.get(si).map(|s| s.as_str())
                        });

                        if let Some(sn) = sample_name {
                            writeln!(&mut writer.barcode_file, "{}_{}", sn, bc_str)
                        } else {
                            writeln!(&mut writer.barcode_file, "{}", bc_str)
                        }
                        .expect("can't write to barcode file.");

                        // write to feature dump file
                        if let Some(sn) = sample_name {
                            writeln!(
                                &mut writer.feature_file,
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                bc_str,
                                sn,
                                (num_mapped + num_unmapped),
                                num_mapped,
                                sum_umi,
                                mapping_rate,
                                dedup_rate,
                                mean_by_max,
                                num_expr,
                                num_genes_over_mean
                            )
                        } else {
                            writeln!(
                                &mut writer.feature_file,
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                bc_str,
                                (num_mapped + num_unmapped),
                                num_mapped,
                                sum_umi,
                                mapping_rate,
                                dedup_rate,
                                mean_by_max,
                                num_expr,
                                num_genes_over_mean
                            )
                        }
                        .expect("can't write to feature file");
                    } // lock on bc_writer released here (end of scope)

                    // Accumulate MTX triplets in thread-local buffer (no lock held)
                    for (&ind, &val) in expressed_ind.iter().zip(expressed_vec.iter()) {
                        local_triplets.push((row_index, ind, val));
                    }

                    // Record bootstrap summary stats (mean/var) as sparse triplets
                    if config.num_bootstraps > 0 && !bootstraps.is_empty() {
                        if config.summary_stat {
                            boot_helper.record_cell(row_index, &bootstraps);
                        } else {
                            boot_helper.record_cell_from_replicates(row_index, &bootstraps);
                        }
                    }
                } // end of cell processing

                // if we are dumping the equivalence class output, fill in
                // the in-memory representation here.
                if config.dump_eq {
                    let eqmap_deref = shared.eqid_map_lock.lock();
                    let geqmap = &mut *eqmap_deref.unwrap();
                    // the next available global id for a gene-level
                    // equivalence class
                    let mut next_id = geqmap.global_eqc.len() as u64;
                    for (labels, payload) in gene_eqc.iter() {
                        let count = payload.count();
                        let mut found = true;
                        match geqmap.global_eqc.get(labels) {
                            Some(eqid) => {
                                geqmap.cell_level_count.push((*eqid, count));
                            }
                            None => {
                                found = false;
                                geqmap.cell_level_count.push((next_id, count));
                            }
                        }
                        if !found {
                            geqmap.global_eqc.insert(labels.to_vec().clone(), next_id);
                            next_id += 1;
                        }
                    }
                    //let bc_mer: BitKmer = (bc, bclen as u8);
                    geqmap.cell_offset.push((row_index, gene_eqc.len()));
                }
                // Reset the gene eqc map only if it was used (i.e., NOT the fast path).
                // For fast-path cells (nrec < SMALL_CELL_FAST_THRESHOLD), gene_eqc was
                // never touched, so clearing it is pure waste — and the O(capacity) clear
                // was the 95% bottleneck when processing millions of tiny cells.
                if nrec >= SMALL_CELL_FAST_THRESHOLD as u32 {
                    if gene_eqc.capacity() > 256 {
                        gene_eqc = HashMap::with_hasher(ahash::RandomState::with_seeds(
                            2u64, 7u64, 1u64, 8u64,
                        ));
                    } else {
                        gene_eqc.clear();
                    }
                }
            } // for all cells in this meta chunk
        } // while we can get work
    } // while cells remain
    (local_nrec, local_triplets, boot_helper)
}

pub(crate) fn do_quantify<T: BufRead, B, R, P>(
    mut br: T,
    quant_opts: QuantOpts,
    prelude: RadPrelude,
    file_tag_map: TagMap,
    sample_bc_extractor: Option<Arc<dyn Fn(&R) -> usize + Send + Sync>>,
) -> anyhow::Result<()>
where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord
        + CollatableMappedRecord<B>
        + KnownSize
        + UmiTaggedRecord
        + OptionalAlignmentExtras
        + 'static,
    <R as MappedRecord>::ParsingContext: RecordContext,
    <R as MappedRecord>::ParsingContext: Clone,
    <R as MappedRecord>::ParsingContext: Send,
    P: EqClassPayload,
{
    let parent = std::path::Path::new(quant_opts.input_dir);
    // Load collation manifest if present (multi-barcode data).
    // Build a lookup from cell_num (chunk index) → sample index,
    // plus a list of sample names. This is used to write sample-prefixed
    // barcodes and the sample_name column in featureDump.txt directly
    // during quantification, avoiding a non-deterministic post-processing step.
    let manifest_path = parent.join("collation_manifest.bin");
    // For multi-barcode data, build a Vec of sample names indexed by sample index.
    // The collation stores the integer sample index in barcodes[0] of each record,
    // so quant reads it directly and indexes into this Vec — no HashMap needed.
    let sample_names = if manifest_path.exists() {
        let manifest = CollationManifest::read_from_file(&manifest_path)?;
        let names: Vec<String> = manifest
            .sample_groups
            .iter()
            .map(|g| {
                g.name
                    .as_deref()
                    .unwrap_or(&format!("{:x}", g.key))
                    .to_string()
            })
            .collect();
        Some(Arc::new(names))
    } else {
        None
    };

    let init_uniform = quant_opts.init_uniform;
    let summary_stat = quant_opts.summary_stat;
    let dump_eq = quant_opts.dump_eq;
    let resolution = quant_opts.resolution;
    let pug_exact_umi = quant_opts.pug_exact_umi;
    let mut sa_model = quant_opts.sa_model;
    let _small_thresh = quant_opts.small_thresh;
    let large_graph_thresh = quant_opts.large_graph_thresh;
    let filter_list = quant_opts.filter_list;
    let log = quant_opts.log;
    let num_threads = quant_opts.num_threads;
    let num_bootstraps = quant_opts.num_bootstraps;

    let hdr = &prelude.hdr;
    // in the collated rad file, we have 1 cell per chunk.
    // we make this value `mut` since, if we have a non-empty
    // filter list, the number of cells will be dictated by
    // it's length.
    let mut num_cells = hdr.num_chunks;

    info!(
        log,
        "paired : {:?}, ref_count : {}, num_chunks : {}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en)
    );

    // now that we have the header, parse and convert the
    // tgmap.

    // first, build a hash of each transcript to it's index
    let rnhasher = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut rname_to_id: HashMap<String, u32, ahash::RandomState> =
        HashMap::with_capacity_and_hasher(hdr.ref_count as usize, rnhasher);
    for (i, n) in hdr.ref_names.iter().enumerate() {
        rname_to_id.insert(n.clone(), i as u32);
    }
    //println!("{:?}", hdr.ref_names);

    // will hold the unique gene names in the order they are encountered
    let mut gene_names: Vec<String> = Vec::with_capacity((hdr.ref_count / 2) as usize);
    let gnhasher = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut gene_name_to_id: HashMap<String, u32, ahash::RandomState> =
        HashMap::with_hasher(gnhasher);

    let usa_mode;
    let tid_to_gid;
    // parse the tg-map; this is expected to be a 2-column
    // tsv file if we are dealing with one status of transcript
    // e.g. just spliced, or 3-column tsv if we are dealing with
    // both spliced and unspliced.  The type will be automatically
    // determined.
    match afutils::parse_tg_map(
        quant_opts.tg_map,
        hdr.ref_count as usize,
        &rname_to_id,
        &mut gene_names,
        &mut gene_name_to_id,
    ) {
        Ok((v, us)) => {
            tid_to_gid = v;
            usa_mode = us;
            if usa_mode {
                assert_eq!(
                    num_bootstraps, 0,
                    "currently USA-mode (all-in-one unspliced/spliced/ambiguous) analysis cannot be used with bootstrapping."
                );

                match resolution {
                    ResolutionStrategy::Parsimony
                    | ResolutionStrategy::ParsimonyEm
                    | ResolutionStrategy::ParsimonyGene
                    | ResolutionStrategy::ParsimonyGeneEm => {
                        info!(
                            log,
                            "currently USA-mode (all-in-one unspliced/spliced/ambiguous) analysis using parsimony(-gene) or parsimony(-gene)-em resolution is EXPERIMENTAL."
                        );
                    }
                    _ => {}
                }
            } else {
                // the SplicedAmbiguityModel of PreferAmbiguity only makes sense when we are
                // operating `usa_mode`, so if the user has set that here, inform them
                // it will be changed back to winner-take-all
                match sa_model {
                    SplicedAmbiguityModel::WinnerTakeAll => {}
                    _ => {
                        info!(
                            log,
                            "When not operating in USA-mode (all-in-one unspliced/spliced/ambiguous), the SplicedAmbiguityModel will be ignored."
                        );
                        sa_model = SplicedAmbiguityModel::WinnerTakeAll;
                    }
                }
            }
        }
        Err(e) => {
            return Err(e);
        }
    }

    info!(
        log,
        "tg-map contained {} genes mapping to {} transcripts.",
        gene_names.len().to_formatted_string(&Locale::en),
        tid_to_gid.len().to_formatted_string(&Locale::en)
    );

    // read the map for the number of unmapped reads per corrected barcode
    let bc_unmapped_map: Arc<libradicl::unmapped::CollatedUnmappedCounts> = Arc::new(
        libradicl::unmapped::CollatedUnmappedCounts::read_from_file(
            &parent.join("unmapped_bc_count_collated.bin"),
        )
        .unwrap_or_else(|_| {
            libradicl::unmapped::CollatedUnmappedCounts::new_single(
                libradicl::rad_types::RadIntId::U32,
            )
        }),
    );

    // file-level
    let fl_tags = &prelude.file_tags;
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = &prelude.read_tags;
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = &prelude.aln_tags;
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    //let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    info!(log, "File-level tag values {:?}", file_tag_map);

    // Get barcode length: standard files use "cblen", multi-barcode files
    // use "b{N-1}len" where N is the number of barcodes. Try both.
    let barcode_tag = file_tag_map
        .get("cblen")
        .or_else(|| {
            // Multi-barcode: try b1len, b0len, etc.
            file_tag_map
                .get("b1len")
                .or_else(|| file_tag_map.get("b0len"))
        })
        .expect("tag map must contain cblen or bNlen for barcode length");
    let barcode_len: u16 = barcode_tag.try_into()?;

    // if we have a filter list, extract it here
    let mut retained_bc: Option<HashSet<u64, ahash::RandomState>> = None;
    if let Some(fname) = filter_list {
        match afutils::read_filter_list(fname, barcode_len) {
            Ok(fset) => {
                // the number of cells we expect to
                // actually process
                num_cells = fset.len() as u64;
                retained_bc = Some(fset);
            }
            Err(e) => {
                return Err(e);
            }
        }
    }

    let mut _num_reads: usize = 0;

    let pbar = ProgressBar::with_draw_target(
        Some(num_cells),
        ProgressDrawTarget::stderr_with_hz(5u8), // update at most 5 times/sec.
    );
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .expect("ProgressStyle template was invalid.")
            .progress_chars("╢▌▌░╟"),
    );

    // Trying this parallelization strategy to avoid
    // many temporary data structures.

    // We have a work queue that contains collated chunks
    // (i.e. data at the level of each cell).  One thread
    // populates the queue, and the remaining worker threads
    // pop a chunk, perform the quantification, and update the
    // output.  The updating of the output requires acquiring
    // two locks (1) to update the data in the matrix and
    // (2) to write to the barcode file.  We also have to
    // decrement an atomic coutner for the numebr of cells that
    // remain to be processed.

    // create a thread-safe queue based on the number of worker threads
    let n_workers = if num_threads > 1 {
        (num_threads - 1) as usize
    } else {
        1
    };
    let mut chunk_reader = libradicl::readers::ParallelChunkReader::<R>::new(
        &prelude,
        std::num::NonZeroUsize::new(n_workers).unwrap(),
    );

    // the number of cells left to process
    let cells_to_process = Arc::new(AtomicUsize::new(num_cells as usize));
    // each thread needs a *read-only* copy of this transcript <-> gene map
    let tid_to_gid_shared = std::sync::Arc::new(tid_to_gid);
    // the number of reference sequences
    let ref_count = hdr.ref_count as u32;

    // the number of genes (different than the number of reference sequences, which are transcripts)
    let num_genes = gene_name_to_id.len();

    // create our output directory
    let output_path = std::path::Path::new(quant_opts.output_dir);
    fs::create_dir_all(output_path)?;

    // create sub-directory for matrix
    let output_matrix_path = output_path.join("alevin");
    fs::create_dir_all(&output_matrix_path)?;

    // well need a protected handle to write out the barcode
    let bc_path = output_matrix_path.join("quants_mat_rows.txt");
    let bc_file = fs::File::create(bc_path)?;

    let _boot_helper = BootstrapHelper::new(output_path, num_bootstraps, summary_stat);

    let ff_path = output_path.join("featureDump.txt");
    let mut ff_file = fs::File::create(ff_path)?;
    if sample_names.is_some() {
        writeln!(
            ff_file,
            "CB\tsample_name\tCorrectedReads\tMappedReads\tDeduplicatedReads\tMappingRate\tDedupRate\tMeanByMax\tNumGenesExpressed\tNumGenesOverMean"
        )?;
    } else {
        writeln!(
            ff_file,
            "CB\tCorrectedReads\tMappedReads\tDeduplicatedReads\tMappingRate\tDedupRate\tMeanByMax\tNumGenesExpressed\tNumGenesOverMean"
        )?;
    }
    let alt_res_cells = Arc::new(Mutex::new(Vec::<u64>::new()));
    let empty_resolved_cells = Arc::new(Mutex::new(Vec::<u64>::new()));

    // Estimate initial triplet capacity for MTX output. The 10% density assumption
    // is far too high for large multiplexed datasets (actual density is often <0.5%).
    // Cap at 256M entries (~3GB) to avoid overallocation; the vec will grow as needed.
    let _tmcap = {
        let estimate = (0.1f64 * num_genes as f64 * num_cells as f64).round() as usize;
        estimate.min(256_000_000)
    };

    // the length of the vector of gene counts we'll use
    let num_rows = if usa_mode {
        // the number of genes should be the max gene id + 1
        // over the gene ids in gene_name_to_id.  The +2 is
        // because the ids in gene_name_to_id are only for
        // the spliced genes, leaving a space in between for each
        // unspliced variant; so +1 to get the largest valid index and
        // another +1 to get the size (ids are 0 based).
        let mid = (gene_name_to_id
            .values()
            .max()
            .expect("gene name to id map should not be empty.")
            + 2) as usize;

        // spliced, unspliced, ambiguous for each gene
        // but num genes already accounts for spliced & unspliced
        mid + (mid / 2)
    } else {
        num_genes
    };

    let usa_offsets = if usa_mode {
        Some((num_rows / 3, (2 * num_rows / 3)))
    } else {
        None
    };

    let bc_writer = Arc::new(Mutex::new(QuantOutputInfo {
        barcode_file: BufWriter::new(bc_file),
        feature_file: BufWriter::new(ff_file),
        row_index: 0usize,
    }));

    let mmrate = Arc::new(Mutex::new(vec![0f64; num_cells as usize]));

    type WorkerResult = (usize, Vec<(usize, usize, f32)>, BootstrapHelper);
    let mut thread_handles: Vec<thread::JoinHandle<WorkerResult>> = Vec::with_capacity(n_workers);

    // This is the hash table that will hold the global
    // (i.e. across all cells) gene-level equivalence
    // classes.  We have _not_ yet figured out an efficient
    // way to do this in a lock-free manner, so this
    // structure is protected by a lock for now.
    // This will only be used if the `dump_eq` paramater is true.
    let so = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let eqid_map_lock = Arc::new(Mutex::new(EqcMap {
        global_eqc: HashMap::with_hasher(so),
        cell_level_count: Vec::new(),
        cell_offset: Vec::new(),
    }));

    // for each worker, spawn off a thread
    for worker in 0..n_workers {
        // each thread will need to access the work queue
        //let in_q = q.clone();
        let in_q = chunk_reader.get_queue();
        let is_done = chunk_reader.is_done();
        // and the logger
        let log = log.clone();
        // the shared tid_to_gid map
        let tid_to_gid = tid_to_gid_shared.clone();
        // and the atomic counter of remaining work
        let cells_remaining = cells_to_process.clone();

        // and the file writer
        let bcout = bc_writer.clone();
        // global gene-level eqc map
        let eqid_map_lockc = eqid_map_lock.clone();
        // and will need to know the barcode length
        let alt_res_cells = alt_res_cells.clone();
        let empty_resolved_cells = empty_resolved_cells.clone();
        let unmapped_count = bc_unmapped_map.clone();
        let mmrate = mmrate.clone();

        // if we are performing parsimony-gene or parsimony-gene-em
        // resolution, then the equivalence classes will be immediately
        // projected to the gene level.
        let eq_map_type = match resolution {
            ResolutionStrategy::ParsimonyGene | ResolutionStrategy::ParsimonyGeneEm => {
                EqMapType::GeneLevel
            }
            _ => EqMapType::TranscriptLevel,
        };

        let num_eq_targets = match eq_map_type {
            EqMapType::TranscriptLevel => ref_count,
            EqMapType::GeneLevel => {
                // get the max spliced gene ID and add 1 to get the unspliced ID
                // and another 1 to get the size.
                gene_name_to_id
                    .values()
                    .max()
                    .expect("gene name to id map should not be empty.")
                    + 2
            }
        };

        let config = WorkerConfig {
            resolution,
            usa_mode,
            usa_offsets,
            em_init_type: if init_uniform {
                EmInitType::Uniform
            } else {
                EmInitType::Informative
            },
            large_graph_thresh,
            pug_exact_umi,
            sa_model,
            num_bootstraps,
            init_uniform,
            summary_stat,
            dump_eq,
            num_genes,
            num_rows,
            barcode_len,
        };

        let shared = WorkerSharedState {
            in_q,
            is_done,
            tid_to_gid,
            cells_remaining,
            bcout,
            eqid_map_lock: eqid_map_lockc,
            alt_res_cells,
            empty_resolved_cells,
            unmapped_count,
            mmrate,
            sample_names: sample_names.clone(),
            sample_idx_extractor: sample_bc_extractor.clone(),
        };

        // now, make the worker thread
        let handle = std::thread::spawn(move || {
            run_worker_thread::<_, _, P>(worker, config, shared, log, num_eq_targets, eq_map_type)
        });

        thread_handles.push(handle);
    }

    let cb = |_new_bytes: u64, new_rec: u64| {
        pbar.inc(new_rec);
    };

    // push the work onto the queue for the worker threads
    // we spawned above.
    let _ = if let Some(ret_bc) = retained_bc {
        let filter_fn =
            |buf: &[u8], record_context: &<R as MappedRecord>::ParsingContext| -> bool {
                let ch = R::peek_collatable_header(&buf[8..], record_context)
                    .expect("at least one record");
                let ck: u64 = ch.collate_key().into();
                ret_bc.contains(&ck)
            };
        chunk_reader.start_filtered(&mut br, filter_fn, Some(cb))
    } else {
        chunk_reader.start(&mut br, Some(cb))
    };

    let gn_path = output_matrix_path.join("quants_mat_cols.txt");
    let gn_file = File::create(gn_path).expect("couldn't create gene name file.");
    let mut gn_writer = BufWriter::new(gn_file);

    // if we are not using unspliced then just write the gene names
    if !usa_mode {
        for g in gene_names {
            gn_writer.write_all(format!("{}\n", g).as_bytes())?;
        }
    } else {
        // otherwise, we write the spliced names, the unspliced names, and then
        // the ambiguous names
        for g in gene_names.iter() {
            gn_writer.write_all(format!("{}\n", *g).as_bytes())?;
        }
        // unspliced
        for g in gene_names.iter() {
            gn_writer.write_all(format!("{}-U\n", *g).as_bytes())?;
        }
        // ambiguous
        for g in gene_names.iter() {
            gn_writer.write_all(format!("{}-A\n", *g).as_bytes())?;
        }
    }

    let mut total_records = 0usize;
    let mut all_triplets: Vec<(usize, usize, f32)> = Vec::new();
    let mut all_boot_mean_triplets: Vec<(usize, usize, f32)> = Vec::new();
    let mut all_boot_var_triplets: Vec<(usize, usize, f32)> = Vec::new();
    for h in thread_handles {
        match h.join() {
            Ok((rc, triplets, boot)) => {
                total_records += rc;
                all_triplets.extend(triplets);
                if boot.num_bootstraps > 0 {
                    all_boot_mean_triplets.extend(boot.mean_triplets);
                    all_boot_var_triplets.extend(boot.var_triplets);
                }
            }
            Err(_e) => {
                info!(log, "thread panicked");
            }
        }
    }

    // Write the MTX output
    info!(
        log,
        "building triplet matrix from {} entries",
        all_triplets.len().to_formatted_string(&Locale::en),
    );

    let mut trimat = sprs::TriMatI::<f32, u32>::with_capacity(
        (num_cells as usize, num_rows),
        all_triplets.len(),
    );
    for (row, col, val) in all_triplets {
        trimat.add_triplet(row, col, val);
    }

    let mtx_path = output_matrix_path.join("quants_mat.mtx");
    sprs::io::write_matrix_market(mtx_path, &trimat)?;

    // Write bootstrap summary stat matrices if bootstraps were computed
    if num_bootstraps > 0 && !all_boot_mean_triplets.is_empty() {
        let mut mean_trimat = sprs::TriMatI::<f32, u32>::with_capacity(
            (num_cells as usize, num_rows),
            all_boot_mean_triplets.len(),
        );
        for (row, col, val) in all_boot_mean_triplets {
            mean_trimat.add_triplet(row, col, val);
        }
        let mean_path = output_matrix_path.join("bootstraps_mean.mtx");
        sprs::io::write_matrix_market(mean_path, &mean_trimat)?;

        let mut var_trimat = sprs::TriMatI::<f32, u32>::with_capacity(
            (num_cells as usize, num_rows),
            all_boot_var_triplets.len(),
        );
        for (row, col, val) in all_boot_var_triplets {
            var_trimat.add_triplet(row, col, val);
        }
        let var_path = output_matrix_path.join("bootstraps_var.mtx");
        sprs::io::write_matrix_market(var_path, &var_trimat)?;

        info!(
            log,
            "wrote bootstrap summary statistics: mean ({} entries), variance ({} entries)",
            mean_trimat.nnz(),
            var_trimat.nnz(),
        );
    }

    let pb_msg = format!(
        "finished quantifying {} cells.",
        num_cells.to_formatted_string(&Locale::en)
    );
    pbar.finish_with_message(pb_msg);

    info!(
        log,
        "processed {} total read records",
        total_records.to_formatted_string(&Locale::en)
    );

    if dump_eq {
        write_eqc_counts(&eqid_map_lock, num_rows, usa_mode, &output_matrix_path, log)?;
    }

    let meta_info = json!({
    "cmd" : quant_opts.cmdline,
    "version_str": quant_opts.version,
    "resolution_strategy" : resolution.to_string(),
    "num_quantified_cells" : num_cells,
    "num_genes" : num_rows,
    "dump_eq" : dump_eq,
    "usa_mode" : usa_mode,
    "alt_resolved_cell_numbers" : *alt_res_cells.lock().unwrap(),
    "empty_resolved_cell_numbers" : *empty_resolved_cells.lock().unwrap(),
    "quant_options" : quant_opts
    });

    let mut meta_info_file =
        File::create(output_path.join("quant.json")).expect("couldn't create quant.json file.");
    let aux_info_str = serde_json::to_string_pretty(&meta_info).expect("could not format json.");
    meta_info_file
        .write_all(aux_info_str.as_bytes())
        .expect("cannot write to quant.json file");

    // k3yavi: Todo delete after api stability
    // creating a dummy cmd_info.json for R compatibility
    /*
    let cmd_info = json!({
     "salmon_version": "1.4.0",
     "auxDir": "aux_info"
    });
    let mut cmd_info_file = File::create(output_path.join("quant_cmd_info.json"))
    .expect("couldn't create quant_cmd_info.json file.");
    let cmd_info_str =
    serde_json::to_string_pretty(&cmd_info).expect("could not format quant_cmd_info json.");
    cmd_info_file
    .write_all(cmd_info_str.as_bytes())
    .expect("cannot write to quant_cmd_info.json file");
    */
    Ok(())
}

// TODO: see if we'd rather pass an structure
// with these options
pub fn do_quantify_dispatch<T: BufRead>(mut br: T, quant_opts: QuantOpts) -> anyhow::Result<()> {
    let log = quant_opts.log;
    let prelude = RadPrelude::from_bytes(&mut br)?;
    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br).unwrap();

    let rec_type = afutils::get_record_type_from_prelude(&prelude, &file_tag_map);

    match rec_type {
        KnownRecordType::RnaLong(_bc_len) => {
            info!(log, "record type is long read single-cell RNA-seq");
            do_quantify::<_, u64, ScLongReadRecord, LongReadEqClassPayload>(
                br,
                quant_opts,
                prelude,
                file_tag_map,
                None,
            )
        }
        KnownRecordType::AtacSeq(_bc_len) => {
            info!(log, "record type is short read single-cell ATAC-seq");
            anyhow::bail!("To process atac-seq data, you should use the \"atac\" sub-command");
        }
        KnownRecordType::RnaShortPos(_bc_len) => {
            info!(
                log,
                "record type is short read single-cell RNA-seq with positions"
            );
            do_quantify::<_, u64, AlevinFryReadRecordWithPosition, BasicEqClassPayload>(
                br,
                quant_opts,
                prelude,
                file_tag_map,
                None,
            )
        }
        KnownRecordType::RnaShort(_bc_len) => {
            info!(
                log,
                "record type is standard short read single-cell RNA-seq"
            );
            do_quantify::<_, u64, AlevinFryReadRecord, BasicEqClassPayload>(
                br,
                quant_opts,
                prelude,
                file_tag_map,
                None,
            )
        }
        KnownRecordType::RnaShortMultiBC(cell_bc_len, num_bc) => {
            info!(
                log,
                "record type is multi-barcode single-cell RNA-seq ({} barcode levels, cell BC len = {})",
                num_bc,
                cell_bc_len,
            );

            // Run quantification using MultiBarcodeReadRecord.
            // The collation has grouped records by (sample, cell) and each chunk
            // is one cell's data. collate_key() returns the cell barcode.
            // The sample index extractor reads barcodes[0] (the integer sample
            // index written by the scatter phase) from each record.
            let extractor: Arc<dyn Fn(&MultiBarcodeReadRecord) -> usize + Send + Sync> =
                Arc::new(|rec: &MultiBarcodeReadRecord| -> usize { rec.barcodes[0] as usize });
            do_quantify::<_, u64, MultiBarcodeReadRecord, BasicEqClassPayload>(
                br,
                quant_opts,
                prelude,
                file_tag_map,
                Some(extractor),
            )?;

            Ok(())
        }
    }
}

// TODO: see if we'd rather pass an structure
// with these options
pub fn velo_quantify(_quant_opts: QuantOpts) -> anyhow::Result<()> {
    unimplemented!("not implemented on this branch yet");
    //Ok(())
}
