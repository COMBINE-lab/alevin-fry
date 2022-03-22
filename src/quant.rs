/*
 * Copyright (c) 2020-2022 Rob Patro, Avi Srivastava, Hirak Sarkar, Dongze He, Mohsen Zakeri.
 *
 * This file is part of alevin-fry
 * (see https://github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use crossbeam_queue::ArrayQueue;
use indicatif::{ProgressBar, ProgressStyle};

#[allow(unused_imports)]
use slog::{crit, info, warn};

use needletail::bitkmer::*;
use num_format::{Locale, ToFormattedString};
use scroll::Pread;
use serde_json::json;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::str::FromStr;
use std::string::ToString;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use std::fmt;
//use std::ptr;

use flate2::write::GzEncoder;
use flate2::Compression;

use crate::em::{em_optimize, em_optimize_subset, run_bootstrap, EmInitType};
use crate::eq_class::{EqMap, IndexedEqList};
use crate::io_utils;
use crate::pugutils;
use crate::utils as afutils;
use libradicl::rad_types;

type BufferedGzFile = BufWriter<GzEncoder<fs::File>>;

#[derive(PartialEq, Debug, Clone, Copy)]
pub enum SplicedAmbiguityModel {
    PreferAmbiguity,
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

#[derive(PartialEq, Debug, Clone, Copy)]
pub enum ResolutionStrategy {
    Trivial,
    CellRangerLike,
    CellRangerLikeEm,
    Full,
    Parsimony,
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
            "parsimony-em" | "full" => Ok(ResolutionStrategy::Full),
            "parsimony" => Ok(ResolutionStrategy::Parsimony),
            _ => Err("no match"),
        }
    }
}

struct BootstrapHelper {
    bsfile: Option<BufferedGzFile>,
    mean_var_files: Option<(BufferedGzFile, BufferedGzFile)>,
}

impl BootstrapHelper {
    fn new(
        output_path: &std::path::Path,
        num_bootstraps: u32,
        summary_stat: bool,
    ) -> BootstrapHelper {
        if num_bootstraps > 0 {
            if summary_stat {
                let bootstrap_mean_path = output_path.join("bootstraps_mean.eds.gz");
                let bootstrap_var_path = output_path.join("bootstraps_var.eds.gz");
                let bt_mean_buffered = GzEncoder::new(
                    fs::File::create(bootstrap_mean_path).unwrap(),
                    Compression::default(),
                );
                let bt_var_buffered = GzEncoder::new(
                    fs::File::create(bootstrap_var_path).unwrap(),
                    Compression::default(),
                );
                BootstrapHelper {
                    bsfile: None,
                    mean_var_files: Some((
                        BufWriter::new(bt_mean_buffered),
                        BufWriter::new(bt_var_buffered),
                    )),
                }
            } else {
                let bootstrap_path = output_path.join("bootstraps.eds.gz");
                let bt_buffered = GzEncoder::new(
                    fs::File::create(bootstrap_path).unwrap(),
                    Compression::default(),
                );
                BootstrapHelper {
                    bsfile: Some(BufWriter::new(bt_buffered)),
                    mean_var_files: None,
                }
            }
        } else {
            BootstrapHelper {
                bsfile: None,
                mean_var_files: None,
            }
        }
    }
}

struct QuantOutputInfo {
    barcode_file: BufWriter<fs::File>,
    eds_file: BufWriter<GzEncoder<fs::File>>,
    feature_file: BufWriter<fs::File>,
    trimat: sprs::TriMatI<f32, u32>,
    row_index: usize,
    bootstrap_helper: BootstrapHelper, //sample_or_mean_and_var: (BufWriter<GzEncoder<fs::File>>)
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
) -> bool {
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
    sprs::io::write_matrix_market(&mtx_path, &eqmat).expect("could not write geqc_counts.mtx");

    // write the sets of genes that define each eqc
    let gn_eq_path = output_path.join("gene_eqclass.txt.gz");
    let mut gn_eq_writer = BufWriter::new(GzEncoder::new(
        fs::File::create(gn_eq_path).unwrap(),
        Compression::default(),
    ));

    // number of genes
    gn_eq_writer
        .write_all(format!("{}\n", num_genes).as_bytes())
        .expect("could not write to gene_eqclass.txt.gz");

    // number of classes
    gn_eq_writer
        .write_all(format!("{}\n", num_eqclasses).as_bytes())
        .expect("could not write to gene_eqclass.txt.gz");

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
        let ambig_offset = (2 * unspliced_offset) as u32;
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
                    .expect("could not write to gene_eqclass.txt.gz")
            }
            gn_eq_writer
                .write_all(format!("{}\n", eqid).as_bytes())
                .expect("could not write to gene_eqclass.txt.gz");
        }
    } else {
        // if we are running the *standard* mode, then the gene_id
        // mapping is unaltered
        for (gene_list, eqid) in geqmap.global_eqc.iter() {
            for g in gene_list.iter() {
                gn_eq_writer
                    .write_all(format!("{}\t", g).as_bytes())
                    .expect("could not write to gene_eqclass.txt.gz");
            }
            gn_eq_writer
                .write_all(format!("{}\n", eqid).as_bytes())
                .expect("could not write to gene_eqclass.txt.gz");
        }
    }
    true
}

// TODO: see if we'd rather pass an structure
// with these options
#[allow(clippy::too_many_arguments)]
pub fn quantify(
    input_dir: String,
    tg_map: String,
    output_dir: String,
    num_threads: u32,
    num_bootstraps: u32,
    init_uniform: bool,
    summary_stat: bool,
    dump_eq: bool,
    use_mtx: bool,
    resolution: ResolutionStrategy,
    sa_model: SplicedAmbiguityModel,
    small_thresh: usize,
    filter_list: Option<&str>,
    cmdline: &str,
    version: &str,
    log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {
    let parent = std::path::Path::new(&input_dir);

    // read the collate metadata
    let collate_md_file =
        File::open(parent.join("collate.json")).expect("could not open the collate.json file.");
    let collate_md: serde_json::Value = serde_json::from_reader(&collate_md_file)?;

    // is the collated RAD file compressed?
    let compressed_input = collate_md["compressed_output"].as_bool().unwrap();

    if compressed_input {
        let i_file =
            File::open(parent.join("map.collated.rad.sz")).expect("run collate before quant");
        let br = snap::read::FrameDecoder::new(BufReader::new(&i_file));

        info!(
            log,
            "quantifying from compressed, collated RAD file {:?}", i_file
        );

        do_quantify(
            input_dir,
            br,
            tg_map,
            output_dir,
            num_threads,
            num_bootstraps,
            init_uniform,
            summary_stat,
            dump_eq,
            use_mtx,
            resolution,
            sa_model,
            small_thresh,
            filter_list,
            cmdline,
            version,
            log,
        )
    } else {
        let i_file = File::open(parent.join("map.collated.rad")).expect("run collate before quant");
        let br = BufReader::new(&i_file);

        info!(
            log,
            "quantifying from uncompressed, collated RAD file {:?}", i_file
        );

        do_quantify(
            input_dir,
            br,
            tg_map,
            output_dir,
            num_threads,
            num_bootstraps,
            init_uniform,
            summary_stat,
            dump_eq,
            use_mtx,
            resolution,
            sa_model,
            small_thresh,
            filter_list,
            cmdline,
            version,
            log,
        )
    }
}

// TODO: see if we'd rather pass an structure
// with these options
#[allow(clippy::too_many_arguments)]
pub fn do_quantify<T: Read>(
    input_dir: String,
    mut br: T,
    tg_map: String,
    output_dir: String,
    num_threads: u32,
    num_bootstraps: u32,
    init_uniform: bool,
    summary_stat: bool,
    dump_eq: bool,
    use_mtx: bool,
    resolution: ResolutionStrategy,
    mut sa_model: SplicedAmbiguityModel,
    small_thresh: usize,
    filter_list: Option<&str>,
    cmdline: &str,
    version: &str,
    log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {
    let parent = std::path::Path::new(&input_dir);
    let hdr = rad_types::RadHeader::from_bytes(&mut br);

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

    let with_unspliced;
    let tid_to_gid;
    // parse the tg-map; this is expected to be a 2-column
    // tsv file if we are dealing with one status of transcript
    // e.g. just spliced, or 3-column tsv if we are dealing with
    // both spliced and unspliced.  The type will be automatically
    // determined.
    match afutils::parse_tg_map(
        &tg_map,
        hdr.ref_count as usize,
        &rname_to_id,
        &mut gene_names,
        &mut gene_name_to_id,
    ) {
        Ok((v, us)) => {
            tid_to_gid = v;
            with_unspliced = us;
            if with_unspliced {
                assert_eq!(
		     num_bootstraps, 0,
		     "currently USA-mode (all-in-one unspliced/spliced/ambiguous) analysis cannot be used with bootstrapping."
		 );
                assert!(
		     matches!(resolution,
			      ResolutionStrategy::CellRangerLike | ResolutionStrategy::CellRangerLikeEm),
		     "currently USA-mode (all-in-one unspliced/spliced/ambiguous) analysis can only be used with cr-like or cr-like-em resolution."
		 );
            } else {
                // the SplicedAmbiguityModel of PreferAmbiguity only makes sense when we are
                // operating `with_unspliced`, so if the user has set that here, inform them
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
    let bc_unmapped_file =
        std::fs::File::open(parent.join("unmapped_bc_count_collated.bin")).unwrap();
    let bc_unmapped_map: Arc<HashMap<u64, u32>> =
        Arc::new(bincode::deserialize_from(&bc_unmapped_file).unwrap());

    // file-level
    let fl_tags = rad_types::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = rad_types::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = rad_types::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    let ft_vals = rad_types::FileTags::from_bytes(&mut br);
    info!(log, "File-level tag values {:?}", ft_vals);

    let bct = rl_tags.tags[0].typeid;
    let umit = rl_tags.tags[1].typeid;

    // if we have a filter list, extract it here
    let mut retained_bc: Option<HashSet<u64, ahash::RandomState>> = None;
    if let Some(fname) = filter_list {
        match afutils::read_filter_list(fname, ft_vals.bclen) {
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

    let pbar = ProgressBar::new(num_cells);
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );
    let ddelta = 500_u64.min(num_cells / 10);
    pbar.set_draw_delta(ddelta);

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
    let q = Arc::new(ArrayQueue::<io_utils::MetaChunk>::new(4 * n_workers));

    // the number of cells left to process
    let cells_to_process = Arc::new(AtomicUsize::new(num_cells as usize));
    // each thread needs a *read-only* copy of this transcript <-> gene map
    let tid_to_gid_shared = std::sync::Arc::new(tid_to_gid);
    // the number of reference sequences
    let ref_count = hdr.ref_count as u32;
    // the types for the barcodes and umis
    let bc_type = rad_types::decode_int_type_tag(bct).expect("unsupported barcode type id.");
    let umi_type = rad_types::decode_int_type_tag(umit).expect("unsupported umi type id.");
    // the number of genes (different than the number of reference sequences, which are transcripts)
    let num_genes = gene_name_to_id.len();

    // create our output directory
    let output_path = std::path::Path::new(&output_dir);
    fs::create_dir_all(output_path)?;

    // create sub-directory for matrix
    let output_matrix_path = output_path.join("alevin");
    fs::create_dir_all(&output_matrix_path)?;

    // well need a protected handle to write out the barcode
    let bc_path = output_matrix_path.join("quants_mat_rows.txt");
    let bc_file = fs::File::create(bc_path)?;

    let mat_path = output_matrix_path.join("quants_mat.gz");
    let boot_helper = BootstrapHelper::new(output_path, num_bootstraps, summary_stat);
    let buffered = GzEncoder::new(fs::File::create(&mat_path)?, Compression::default());

    let ff_path = output_path.join("featureDump.txt");
    let mut ff_file = fs::File::create(ff_path)?;
    writeln!(
	 ff_file,
	 "CB\tCorrectedReads\tMappedReads\tDeduplicatedReads\tMappingRate\tDedupRate\tMeanByMax\tNumGenesExpressed\tNumGenesOverMean"
     )?;
    let alt_res_cells = Arc::new(Mutex::new(Vec::<u64>::new()));
    let empty_resolved_cells = Arc::new(Mutex::new(Vec::<u64>::new()));

    let tmcap = if use_mtx {
        (0.1f64 * num_genes as f64 * num_cells as f64).round() as usize
    } else {
        0usize
    };

    // the length of the vector of gene counts we'll use
    let num_rows = if with_unspliced {
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

    let usa_offsets = if with_unspliced {
        Some(((num_rows / 3) as usize, (2 * num_rows / 3) as usize))
    } else {
        None
    };

    let trimat =
        sprs::TriMatI::<f32, u32>::with_capacity((num_cells as usize, num_rows as usize), tmcap);

    let bc_writer = Arc::new(Mutex::new(QuantOutputInfo {
        barcode_file: BufWriter::new(bc_file),
        eds_file: BufWriter::new(buffered),
        feature_file: BufWriter::new(ff_file),
        trimat,
        row_index: 0usize,
        bootstrap_helper: boot_helper,
    }));

    let mmrate = Arc::new(Mutex::new(vec![0f64; num_cells as usize]));

    let mut thread_handles: Vec<thread::JoinHandle<usize>> = Vec::with_capacity(n_workers);

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
    for _worker in 0..n_workers {
        // each thread will need to access the work queue
        let in_q = q.clone();
        // and the logger
        let log = log.clone();
        // the shared tid_to_gid map
        let tid_to_gid = tid_to_gid_shared.clone();
        // and the atomic counter of remaining work
        let cells_remaining = cells_to_process.clone();
        // they will need to know the bc and umi type
        let bc_type = bc_type;
        let umi_type = umi_type;
        // and the file writer
        let bcout = bc_writer.clone();
        // global gene-level eqc map
        let eqid_map_lockc = eqid_map_lock.clone();
        // and will need to know the barcode length
        let bclen = ft_vals.bclen;
        let alt_res_cells = alt_res_cells.clone();
        let empty_resolved_cells = empty_resolved_cells.clone();
        let unmapped_count = bc_unmapped_map.clone();
        let mmrate = mmrate.clone();

        // now, make the worker thread
        let handle = std::thread::spawn(move || {
            // these can be created once and cleared after processing
            // each cell.
            let mut unique_evidence = vec![false; num_rows];
            let mut no_ambiguity = vec![false; num_rows];
            let mut eq_map = EqMap::new(ref_count);
            let mut expressed_vec = Vec::<f32>::with_capacity(num_genes);
            let mut expressed_ind = Vec::<usize>::with_capacity(num_genes);
            let mut eds_bytes = Vec::<u8>::new();
            let mut bt_eds_bytes: Vec<u8> = Vec::new();
            let mut eds_mean_bytes: Vec<u8> = Vec::new();
            let mut eds_var_bytes: Vec<u8> = Vec::new();

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
            let mut gene_eqc: HashMap<Vec<u32>, u32, ahash::RandomState> = HashMap::with_hasher(s);

            let em_init_type = if init_uniform {
                EmInitType::Uniform
            } else {
                EmInitType::Informative
            };

            // If we are operating in USA-mode with an EM capable resolution
            // method, we'll use (re-use) these variables to hold the USA-mode
            // equivalence class information.
            let mut idx_eq_list = IndexedEqList::new();
            let mut eq_id_count = Vec::<(u32, u32)>::new();

            let mut local_nrec = 0usize;
            // pop MetaChunks from the work queue until everything is
            // processed
            while cells_remaining.load(Ordering::SeqCst) > 0 {
                if let Some((
                    first_cell_in_chunk,
                    cells_in_chunk,
                    _nbytes_total,
                    _nrec_total,
                    buf,
                )) = in_q.pop()
                {
                    // for every cell (chunk) within this meta-chunk
                    let mut byte_offset = 0usize;
                    for cn in 0..cells_in_chunk {
                        cells_remaining.fetch_sub(1, Ordering::SeqCst);
                        let cell_num = first_cell_in_chunk + cn;
                        // nbytes for the current cell
                        let nbytes = buf[byte_offset..].pread::<u32>(0).unwrap();
                        let nrec = buf[byte_offset..].pread::<u32>(4).unwrap();
                        local_nrec += nrec as usize;
                        let mut nbr =
                            BufReader::new(&buf[byte_offset..(byte_offset + nbytes as usize)]);
                        byte_offset += nbytes as usize;

                        let mut c = rad_types::Chunk::from_bytes(&mut nbr, &bc_type, &umi_type);
                        if c.reads.is_empty() {
                            warn!(log, "Discovered empty chunk; should not happen! cell_num = {}, nbytes = {}, nrec = {}", cell_num, nbytes, nrec);
                        }

                        // TODO: Clean up the expect() and merge with the check above
                        // the expect shouldn't happen, but the message is redundant with
                        // the above.  Plus, this would panic if it actually occurred.
                        let bc = c.reads.first().expect("chunk with no reads").bc;

                        // The structures we'll need to hold our output for this
                        // cell.
                        let mut counts: Vec<f32>;
                        let mut alt_resolution = false;

                        let mut bootstraps: Vec<Vec<f32>> = Vec::new();

                        let non_trivial = c.reads.len() >= small_thresh;
                        if non_trivial {
                            // TODO: some testing was done, but see if there
                            // is a better way to set this value.
                            let small_cell = c.reads.len() <= 250;

                            // TODO: Is there an easy / clean way to have similar
                            // optimized code paths for other resolution methods?

                            match resolution {
                                ResolutionStrategy::CellRangerLike
                                | ResolutionStrategy::CellRangerLikeEm => {
                                    if small_cell {
                                        pugutils::get_num_molecules_cell_ranger_like_small(
                                            &mut c,
                                            &tid_to_gid,
                                            num_genes,
                                            &mut gene_eqc,
                                            with_unspliced,
                                            sa_model,
                                            &log,
                                        );
                                    } else {
                                        eq_map.init_from_chunk(&mut c);
                                        pugutils::get_num_molecules_cell_ranger_like(
                                            &eq_map,
                                            &tid_to_gid,
                                            num_genes,
                                            &mut gene_eqc,
                                            with_unspliced,
                                            sa_model,
                                            &log,
                                        );
                                        eq_map.clear();
                                    }
                                    let only_unique =
                                        resolution == ResolutionStrategy::CellRangerLike;

                                    // NOTE: This configuration seems overly complicated
                                    // see if we can simplify it.
                                    match (with_unspliced, only_unique) {
                                        (true, true) => {
                                            // USA mode, only gene-unqique reads
                                            counts = afutils::extract_counts(&gene_eqc, num_rows);
                                        }
                                        (true, false) => {
                                            // USA mode, use EM
                                            afutils::extract_usa_eqmap(
                                                &gene_eqc,
                                                num_rows,
                                                &mut idx_eq_list,
                                                &mut eq_id_count,
                                            );
                                            counts = em_optimize_subset(
                                                &idx_eq_list,
                                                &eq_id_count,
                                                &mut unique_evidence,
                                                &mut no_ambiguity,
                                                em_init_type,
                                                num_rows,
                                                only_unique,
                                                usa_offsets,
                                                &log,
                                            );
                                        }
                                        (false, _) => {
                                            // not USA-mode
                                            counts = em_optimize(
                                                &gene_eqc,
                                                &mut unique_evidence,
                                                &mut no_ambiguity,
                                                em_init_type,
                                                num_genes,
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
                                        &tid_to_gid,
                                        num_genes,
                                        &log,
                                    );
                                    counts = ct.0;
                                    mmrate.lock().unwrap()[cell_num] = ct.1;
                                    eq_map.clear();
                                }
                                ResolutionStrategy::Parsimony => {
                                    eq_map.init_from_chunk(&mut c);
                                    let g = pugutils::extract_graph(&eq_map, &log);
                                    let pug_stats = pugutils::get_num_molecules(
                                        &g,
                                        &eq_map,
                                        &tid_to_gid,
                                        num_genes,
                                        &mut gene_eqc,
                                        &log,
                                    );
                                    alt_resolution = pug_stats.used_alternative_strategy; // alt_res;
                                    counts = em_optimize(
                                        &gene_eqc,
                                        &mut unique_evidence,
                                        &mut no_ambiguity,
                                        em_init_type,
                                        num_genes,
                                        true, // only unqique evidence
                                        &log,
                                    );
                                    eq_map.clear();
                                }
                                ResolutionStrategy::Full => {
                                    eq_map.init_from_chunk(&mut c);
                                    let g = pugutils::extract_graph(&eq_map, &log);
                                    let pug_stats = pugutils::get_num_molecules(
                                        &g,
                                        &eq_map,
                                        &tid_to_gid,
                                        num_genes,
                                        &mut gene_eqc,
                                        &log,
                                    );
                                    alt_resolution = pug_stats.used_alternative_strategy; // alt_res;
                                    counts = em_optimize(
                                        &gene_eqc,
                                        &mut unique_evidence,
                                        &mut no_ambiguity,
                                        em_init_type,
                                        num_genes,
                                        false, // only unqique evidence
                                        &log,
                                    );
                                    eq_map.clear();
                                }
                            }

                            if num_bootstraps > 0 {
                                bootstraps = run_bootstrap(
                                    &gene_eqc,
                                    num_bootstraps,
                                    &counts,
                                    init_uniform,
                                    summary_stat,
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
                                &tid_to_gid,
                                num_genes,
                                &mut gene_eqc,
                                with_unspliced,
                                sa_model,
                                &log,
                            );
                            // USA-mode
                            if with_unspliced {
                                // here, just like for non-USA mode,
                                // we substitute EM with uniform allocation in
                                // this special case
                                match resolution {
                                    ResolutionStrategy::CellRangerLike => {
                                        counts = afutils::extract_counts(&gene_eqc, num_rows);
                                    }
                                    ResolutionStrategy::CellRangerLikeEm => {
                                        counts =
                                            afutils::extract_counts_mm_uniform(&gene_eqc, num_rows);
                                    }
                                    _ => {
                                        counts = vec![0f32; num_genes];
                                        warn!(log, "Should not reach here, only cr-like and cr-like-em are supported in USA-mode.");
                                    }
                                }
                            } else {
                                // non USA-mode
                                counts = vec![0f32; num_genes];
                                for (k, v) in gene_eqc.iter() {
                                    if k.len() == 1 {
                                        counts[*k.first().unwrap() as usize] += *v as f32;
                                    } else {
                                        match resolution {
                                            ResolutionStrategy::CellRangerLikeEm
                                            | ResolutionStrategy::Full => {
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
                            if num_bootstraps > 0 {
                                // TODO: should issue a warning here,
                                // bootstrapping doesn't make sense for
                                // unfiltered data.
                                if summary_stat {
                                    // sample mean = quant
                                    bootstraps.push(counts.clone());
                                    // sample var = 0
                                    bootstraps.push(vec![0f32; num_genes]);
                                } else {
                                    // no variation
                                    for _ in 0..num_bootstraps {
                                        bootstraps.push(counts.clone());
                                    }
                                }
                            } // if the user requested bootstraps
                        } // end of else branch for trivial size cells

                        if alt_resolution {
                            alt_res_cells.lock().unwrap().push(cell_num as u64);
                        }

                        //
                        // featuresStream << "\t" << numRawReads
                        //   << "\t" << numMappedReads
                        let mut max_umi = 0.0f32;
                        let mut sum_umi = 0.0f32;
                        let mut num_expr: u32 = 0;
                        expressed_vec.clear();
                        expressed_ind.clear();

                        for (gn, c) in counts.iter().enumerate() {
                            max_umi = if *c > max_umi { *c } else { max_umi };
                            sum_umi += *c;
                            if *c > 0.0 {
                                num_expr += 1;
                                expressed_vec.push(*c);
                                expressed_ind.push(gn);
                            }
                        }

                        if num_expr == 0 {
                            empty_resolved_cells.lock().unwrap().push(cell_num as u64);
                        }

                        let num_mapped = nrec;
                        let dedup_rate = sum_umi / num_mapped as f32;

                        let num_unmapped = match unmapped_count.get(&bc) {
                            Some(nu) => *nu,
                            None => 0u32,
                        };

                        let mapping_rate = num_mapped as f32 / (num_mapped + num_unmapped) as f32;

                        // mean of the "expressed" genes
                        let mean_expr = sum_umi / num_expr as f32;
                        // number of genes with expression > expressed mean
                        let num_genes_over_mean = expressed_vec.iter().fold(0u32, |acc, x| {
                            if x > &mean_expr {
                                acc + 1u32
                            } else {
                                acc
                            }
                        });
                        // expressed mean / max expression
                        let mean_by_max = mean_expr / max_umi;

                        let row_index: usize; // the index for this row (cell)
                        {
                            // writing the files
                            let bc_mer: BitKmer = (bc, bclen as u8);

                            if !use_mtx {
                                eds_bytes = sce::eds::as_bytes(&counts, num_rows)
                                    .expect("can't convert vector to eds");
                            }

                            // write bootstraps
                            if num_bootstraps > 0 {
                                // flatten the bootstraps
                                if summary_stat {
                                    eds_mean_bytes = sce::eds::as_bytes(&bootstraps[0], num_rows)
                                        .expect("can't convert vector to eds");
                                    eds_var_bytes = sce::eds::as_bytes(&bootstraps[1], num_rows)
                                        .expect("can't convert vector to eds");
                                } else {
                                    for i in 0..num_bootstraps {
                                        let bt_eds_bytes_slice =
                                            sce::eds::as_bytes(&bootstraps[i as usize], num_rows)
                                                .expect("can't convert vector to eds");
                                        bt_eds_bytes.append(&mut bt_eds_bytes_slice.clone());
                                    }
                                }
                            }

                            let writer_deref = bcout.lock();
                            let writer = &mut *writer_deref.unwrap();

                            // get the row index and then increment it
                            row_index = writer.row_index;
                            writer.row_index += 1;

                            // write to barcode file
                            let bc_bytes = &bitmer_to_bytes(bc_mer)[..];
                            writeln!(&mut writer.barcode_file, "{}", unsafe {
                                std::str::from_utf8_unchecked(bc_bytes)
                            })
                            .expect("can't write to barcode file.");

                            // write to matrix file
                            if !use_mtx {
                                // write in eds format
                                writer
                                    .eds_file
                                    .write_all(&eds_bytes)
                                    .expect("can't write to matrix file.");
                            } else {
                                // fill out the triplet matrix in memory
                                for (ind, val) in expressed_ind.iter().zip(expressed_vec.iter()) {
                                    writer.trimat.add_triplet(row_index as usize, *ind, *val);
                                }
                            }
                            writeln!(
                                &mut writer.feature_file,
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                unsafe { std::str::from_utf8_unchecked(bc_bytes) },
                                (num_mapped + num_unmapped),
                                num_mapped,
                                sum_umi,
                                mapping_rate,
                                dedup_rate,
                                mean_by_max,
                                num_expr,
                                num_genes_over_mean
                            )
                            .expect("can't write to feature file");

                            if num_bootstraps > 0 {
                                if summary_stat {
                                    if let Some((meanf, varf)) =
                                        &mut writer.bootstrap_helper.mean_var_files
                                    {
                                        meanf
                                            .write_all(&eds_mean_bytes)
                                            .expect("can't write to bootstrap mean file.");
                                        varf.write_all(&eds_var_bytes)
                                            .expect("can't write to bootstrap var file.");
                                    }
                                } else if let Some(bsfile) = &mut writer.bootstrap_helper.bsfile {
                                    bsfile
                                        .write_all(&bt_eds_bytes)
                                        .expect("can't write to bootstrap file");
                                }
                            } // done bootstrap writing
                        }

                        // if we are dumping the equivalence class output, fill in
                        // the in-memory representation here.
                        if dump_eq {
                            let eqmap_deref = eqid_map_lockc.lock();
                            let geqmap = &mut *eqmap_deref.unwrap();
                            // the next available global id for a gene-level
                            // equivalence class
                            let mut next_id = geqmap.global_eqc.len() as u64;
                            for (labels, count) in gene_eqc.iter() {
                                let mut found = true;
                                match geqmap.global_eqc.get(labels) {
                                    Some(eqid) => {
                                        geqmap.cell_level_count.push((*eqid, *count));
                                    }
                                    None => {
                                        found = false;
                                        geqmap.cell_level_count.push((next_id, *count));
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
                        // clear the gene eqc map
                        gene_eqc.clear();
                    } // for all cells in this meta chunk
                } // while we can get work
            } // while cells remain
            local_nrec
        });

        thread_handles.push(handle);
    }

    // push the work onto the queue for the worker threads
    // we spawned above.
    if let Some(ret_bc) = retained_bc {
        // we have a retained set
        io_utils::fill_work_queue_filtered(
            ret_bc,
            &rl_tags,
            q,
            br,
            hdr.num_chunks as usize,
            &pbar,
        )?;
    } else {
        // we're quantifying everything
        io_utils::fill_work_queue(q, br, hdr.num_chunks as usize, &pbar)?;
    }

    let gn_path = output_matrix_path.join("quants_mat_cols.txt");
    let gn_file = File::create(gn_path).expect("couldn't create gene name file.");
    let mut gn_writer = BufWriter::new(gn_file);

    // if we are not using unspliced then just write the gene names
    if !with_unspliced {
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
    for h in thread_handles {
        match h.join() {
            Ok(rc) => {
                total_records += rc;
            }
            Err(_e) => {
                info!(log, "thread panicked");
            }
        }
    }

    // write to matrix market if we are using it
    if use_mtx {
        let writer_deref = bc_writer.lock();
        let writer = &mut *writer_deref.unwrap();
        writer.eds_file.flush().unwrap();
        // now remove it
        fs::remove_file(&mat_path)?;
        let mtx_path = output_matrix_path.join("quants_mat.mtx");
        sprs::io::write_matrix_market(&mtx_path, &writer.trimat)?;
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
        write_eqc_counts(
            &eqid_map_lock,
            num_rows,
            with_unspliced,
            &output_matrix_path,
            log,
        );
    }

    let meta_info = json!({
    "cmd" : cmdline,
    "version_str": version,
    "resolution_strategy" : resolution.to_string(),
    "num_quantified_cells" : num_cells,
    "num_genes" : num_rows,
    "dump_eq" : dump_eq,
    "usa_mode" : with_unspliced,
    "alt_resolved_cell_numbers" : *alt_res_cells.lock().unwrap(),
    "empty_resolved_cell_numbers" : *empty_resolved_cells.lock().unwrap()
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
#[allow(clippy::too_many_arguments)]
pub fn velo_quantify(
    _input_dir: String,
    _tg_map: String,
    _output_dir: String,
    _num_threads: u32,
    _num_bootstraps: u32,
    _init_uniform: bool,
    _summary_stat: bool,
    _dump_eq: bool,
    _use_mtx: bool,
    _resolution: ResolutionStrategy,
    mut _sa_model: SplicedAmbiguityModel,
    _small_thresh: usize,
    _filter_list: Option<&str>,
    _cmdline: &str,
    _version: &str,
    _log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {
    unimplemented!("not implemented on this branch yet");
    //Ok(())
}
