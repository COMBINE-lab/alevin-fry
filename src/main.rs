/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use anyhow::{anyhow, bail};
use bio_types::strand::Strand;
use clap::{Command, arg, builder::ArgGroup, crate_authors, crate_version, value_parser};
use csv::Error as CSVError;
use csv::ErrorKind;
use itertools::Itertools;
use mimalloc::MiMalloc;
use rand::Rng;
use slog::{Drain, crit, info, o, warn};
use std::path::PathBuf;

use alevin_fry::cellfilter::{CellFilterMethod, generate_permit_list};
use alevin_fry::cmd_parse_utils::{
    pathbuf_directory_exists_validator, pathbuf_file_exists_validator,
};
use alevin_fry::prog_opts::{GenPermitListOpts, QuantOpts};
use alevin_fry::quant::{
    EMMode, EditDistanceMode, EndDistribution, EndModel, ResolutionStrategy, SplicedAmbiguityModel,
    TrueUmiOracleMode,
};
use alevin_fry::truncation::{
    CoverageSource, TruncationAbundance, TruncationEndKernels, TruncationIntervalMode,
    TruncationLengthLaw, TruncationParamScope, TruncationPositionalLaw, TruncationPreset,
    TruncationShrink,
};

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

// grab the version from the Cargo file.
const VERSION: &str = env!("CARGO_PKG_VERSION");

#[allow(dead_code)]
fn gen_random_kmer(k: usize) -> String {
    const CHARSET: &[u8] = b"ACGT";
    let mut rng = rand::rng();
    let s: String = (0..k)
        .map(|_| {
            let idx = rng.random_range(0..CHARSET.len());
            CHARSET[idx] as char
        })
        .collect();
    s
}

#[allow(clippy::manual_clamp)]
fn main() -> anyhow::Result<()> {
    let num_hardware_threads = num_cpus::get() as u32;
    let max_num_threads: String = (num_cpus::get() as u32).to_string();
    let max_num_collate_threads: String = (16_u32.min(num_hardware_threads).max(2_u32)).to_string();
    let max_num_gpl_threads: String = (8_u32.min(num_hardware_threads).max(2_u32)).to_string();

    let crate_authors = crate_authors!("\n");
    let version = crate_version!();

    // capture the entire command line as a string
    let cmdline = std::env::args().join(" ");

    let convert_app = Command::new("convert")
        .about("Convert a BAM file to a RAD file")
        .version(version)
        .author(crate_authors)
        .arg(
            arg!(-b --bam <BAMFILE> "input SAM/BAM file")
                .required(true)
                .value_parser(pathbuf_file_exists_validator),
        )
        .arg(
            arg!(-t --threads <THREADS> "number of threads to use for processing")
                .required(false)
                .value_parser(value_parser!(u32))
                .default_value(max_num_threads.clone()),
        )
        .arg(arg!(-f --filter_best "keep only the highest scoring alignments for each read"))
        .arg(
            arg!(-o --output <RADFILE> "output RAD file")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            arg!(--"score-threshold" <SCORETHRESH> "fraction of the best alignment score to retain for each read")
                .value_parser(value_parser!(f64))
                .default_value("0.95"),
        );

    let view_app = Command::new("view")
        .about("View a RAD file")
        .version(version)
        .author(crate_authors)
        .arg(
            arg!(-r --rad <RADFILE> "input RAD file")
                .required(true)
                .value_parser(pathbuf_file_exists_validator),
        )
        .arg(arg!(-H --header "flag for printing header"));

    let gen_app = Command::new("generate-permit-list")
        .about("Generate a permit list of barcodes from a RAD file")
        .version(version)
        .author(crate_authors)
        .arg(arg!(-i --input <INPUT>  "input directory containing the map.rad RAD file")
            .required(true)
            .value_parser(pathbuf_directory_exists_validator))
        .arg(arg!(-d --"expected-ori" <EXPECTEDORI> "the expected orientation of alignments")
             .required(true)
             .ignore_case(true)
             .value_parser(["fw", "rc", "both", "either"]))
        .arg(arg!(-o --"output-dir" <OUTPUTDIR>  "output directory").required(true).value_parser(value_parser!(PathBuf)))
        .arg(arg!(
            -k --"knee-distance"  "attempt to determine the number of barcodes to keep using the knee distance method."
            )
        )
        .arg(arg!(-t --threads <THREADS> "number of threads to use for the first phase of permit-list generation").value_parser(value_parser!(u32)).default_value(max_num_gpl_threads))
        .arg(arg!(-e --"expect-cells" <EXPECTCELLS> "defines the expected number of cells to use in determining the (read, not UMI) based cutoff")
             .value_parser(value_parser!(usize)))
        .arg(arg!(-f --"force-cells" <FORCECELLS>  "select the top-k most-frequent barcodes, based on read count, as valid (true)")
             .value_parser(value_parser!(usize)))
        .arg(
            arg!(-b --"valid-bc" <VALIDBC> "uses true barcode collected from a provided file")
            .value_parser(pathbuf_file_exists_validator)
        )
        .arg(
            arg!(-u --"unfiltered-pl" <UNFILTEREDPL> "uses an unfiltered external permit list")
            .value_parser(pathbuf_file_exists_validator)
        )
        .group(ArgGroup::new("filter-method")
               .args(["knee-distance", "expect-cells", "force-cells", "valid-bc", "unfiltered-pl"])
               .required(true)
               )
        .arg(
            arg!(-m --"min-reads" <MINREADS> "minimum read count threshold; only used with --unfiltered-pl")
                .value_parser(value_parser!(usize))
                .default_value("10"));

    let collate_app = Command::new("collate")
    .about("Collate a RAD file by corrected cell barcode")
    .version(version)
    .author(crate_authors)
    .arg(arg!(-i --"input-dir" <INPUTDIR> "input directory made by generate-permit-list")
        .required(true)
        .value_parser(pathbuf_directory_exists_validator))
    .arg(arg!(-r --"rad-dir" <RADFILE> "the directory containing the RAD file to be collated")
        .required(true)
        .value_parser(pathbuf_directory_exists_validator))
    .arg(arg!(-t --threads <THREADS> "number of threads to use for processing").value_parser(value_parser!(u32)).default_value(max_num_collate_threads))
    .arg(arg!(-c --compress "compress the output collated RAD file"))
    .arg(arg!(-m --"max-records" <MAXRECORDS> "the maximum number of read records to keep in memory at once")
         .value_parser(value_parser!(u32))
         .default_value("30000000"));

    let quant_app = Command::new("quant")
    .about("Quantify expression from a collated RAD file")
    .version(version)
    .author(crate_authors)
    .arg(arg!(-i --"input-dir" <INPUTDIR>  "input directory containing collated RAD file")
        .required(true)
        .value_parser(pathbuf_directory_exists_validator))
    .arg(arg!(-m --"tg-map" <TGMAP>  "transcript to gene map").required(true).value_parser(pathbuf_file_exists_validator))
    .arg(arg!(--"junction-model" "apply the preceding splice-junction distance model to alignment probabilities")
        .requires("junction-gtf"))
    .arg(arg!(--"junction-gtf" <GTF> "GTF annotation matching the transcript references; required with --junction-model or --junction-likelihood")
        .value_parser(pathbuf_file_exists_validator))
    .arg(arg!(--"junction-likelihood" "with --model-truncation --truncation-params global: learn from the one-candidate molecules of all cells where reads start relative to the previous splice junction, and multiply each interval's prior by that law, normalized over the transcript's start positions (a likelihood, unlike --junction-model's weight)")
        .requires("junction-gtf"))
    .arg(arg!(-o --"output-dir" <OUTPUTDIR> "output directory where quantification results will be written").required(true).value_parser(value_parser!(PathBuf)))
    .arg(arg!(-t --threads <THREADS> "number of threads to use for processing").value_parser(value_parser!(u32)).default_value(max_num_threads.clone()))
    .arg(arg!(-d --"dump-eqclasses" "flag for dumping equivalence classes"))
    .arg(arg!(-b --"num-bootstraps" <NUMBOOTSTRAPS> "number of bootstraps to use").value_parser(value_parser!(u32)).default_value("0"))
    .arg(arg!(--"init-uniform" "flag for uniform sampling").requires("num-bootstraps"))
    .arg(arg!(--"summary-stat" "flag for storing only summary statistics").requires("num-bootstraps"))
    .arg(arg!(--"use-mtx" "flag for writing output matrix in matrix market format (default)"))
    .arg(arg!(--"use-eds" "flag for writing output matrix in EDS format").conflicts_with("use-mtx"))
    .arg(arg!(--"quant-subset" <SFILE> "file containing list of barcodes to quantify, those not in this list will be ignored").value_parser(pathbuf_file_exists_validator))
    .arg(arg!(-r --resolution <RESOLUTION> "the resolution strategy by which molecules will be counted")
        .required(true)
        .ignore_case(true)
        .value_parser(value_parser!(ResolutionStrategy)))
    .arg(arg!(--"sa-model" <SAMODEL> "preferred model of splicing ambiguity")
        .ignore_case(true)
        .value_parser(value_parser!(SplicedAmbiguityModel))
        .default_value("winner-take-all")
        .hide(true))
    .arg(arg!(--"umi-edit-dist" <EDIST> "the Hamming distance within which potentially colliding UMIs will be considered for correction")
        .value_parser(value_parser!(u32))
        .default_value_ifs([
            ("resolution", "cr-like", Some("0")),
            ("resolution", "cr-like-em", Some("0")),
            ("resolution", "trivial", Some("0")),
            ("resolution", "parsimony", Some("1")),
            ("resolution", "parsimony-em", Some("1")),
            ("resolution", "parsimony-gene", Some("1")),
            ("resolution", "parsimony-gene-em", Some("1")),
        ])
        .hide(true))
    .arg(arg!(--"large-graph-thresh" <NVERT> "the order (number of nodes) of a PUG above which the alternative resolution strategy will be applied")
        .value_parser(value_parser!(usize))
        .default_value_ifs([
            ("resolution", "parsimony-gene-em", Some("1000")),
            ("resolution", "parsimony", Some("1000")),
            ("resolution", "parsimony-em", Some("1000")),
            ("resolution", "parsimony-gene", Some("1000")),
        ])
        .default_value("0") // for any other mode
        .hide(true))
    .arg(arg!(--"small-thresh" <SMALLTHRESH> "cells with fewer than these many reads will be resolved using a custom approach")
        .value_parser(value_parser!(usize))
        .default_value("10")
        .hide(true))
    .arg(arg!(--"lambda-size" <LAMBDASIZE> "reward for including a vertex in DP subtree selection")
        .value_parser(value_parser!(f64))
        .default_value("0.0"))
    .arg(arg!(--"tau-delta" <TAUDELTA> "minimum delta(v,t) = score(v,t) - best_alt(v) for transcript admissibility")
        .value_parser(value_parser!(f64))
        .default_value("0.0"))
    .arg(arg!(--"true-umi-mode" <MODE> "how (if at all) to use ground-truth UMIs to compute an oracle upper bound")
        .ignore_case(true)
        .value_parser(value_parser!(TrueUmiOracleMode))
        .default_value("off"))
    .arg(arg!(--"true-txp-mode" <TxpMODE> "if using the true transcript in the oracle mode")
        .ignore_case(true)
        .value_parser(value_parser!(bool)))
    .arg(arg!(--"true-umi-txp-file" <TRUEUMIFILE> "TSV file mapping read names to true UMI/gene ground truth; required unless --true-umi-mode=off")
        .value_parser(pathbuf_file_exists_validator))
    .arg(arg!(--"ed-distance-model" <EDMODE> "edit distance model to use for UMI correction")
        .ignore_case(true)
        .value_parser(value_parser!(EditDistanceMode))
        .default_value("general"))
    .arg(arg!(--"ed-upper-bound" <EDUPPERBOUND> "upper bound for edit distance")
        .value_parser(value_parser!(usize))
        .default_value("2"))
    .arg(arg!(--"em-prob-model" <EMMODE> "probability model to use for EM-based UMI correction")
        .ignore_case(true)
        .value_parser(value_parser!(EMMode))
        .default_value("uniform"))
    .arg(arg!(--"end-model" <ENDMODEL> "positional model for distance from the transcript 3' end")
        .ignore_case(true)
        .value_parser(value_parser!(EndModel))
        .default_value("off"))
    .arg(arg!(--"end-distribution" <ENDDISTRIBUTION> "distribution for distance from the transcript 3' end")
        .ignore_case(true)
        .value_parser(value_parser!(EndDistribution))
        .default_value("normal"))
    .arg(arg!(--"end-threshold" <ENDTHRESHOLD> "distance from the transcript 3' end before the penalty starts")
        .value_parser(value_parser!(f64))
        .default_value("10.0"))
    .arg(arg!(--"end-scale" <ENDSCALE> "scale of the transcript 3' end penalty")
        .value_parser(value_parser!(f64))
        .default_value("10.0"))
    .arg(arg!(--"coverage-model" "use per-cell coverage probabilities during long-read EM"))
    .arg(arg!(--"coverage-dedup" "estimate coverage after an initial deduplication, then repeat deduplication with coverage-adjusted alignment probabilities"))
    .arg(arg!(--"coverage-source" <SOURCE> "where the coverage models take their coverage from: cell (each cell's own molecules) or global (the molecules of all cells, in a pass before quantifying)")
        .ignore_case(true)
        .value_parser(value_parser!(CoverageSource))
        .default_value("cell"))
    .arg(arg!(--"model-truncation" "fit the positional (truncation) model to each cell's deduplicated molecules; see --truncation-params and --truncation-abundance for what is shared across cells"))
    .arg(arg!(--"truncation-preset" <NAME> "a named configuration of the truncation model (implies --model-truncation; a truncation option given explicitly keeps its value): recommended = the configuration that did best on simulated and experimental long-read single-cell data: --truncation-params global --truncation-abundance shrink --truncation-shrink 1 --truncation-end-kernels absolute --truncation-positional-power 0.5 --truncation-group-shrink auto, with --truncation-discover-groups unless --truncation-groups is given")
        .ignore_case(true)
        .value_parser(value_parser!(TruncationPreset)))
    .arg(arg!(--"truncation-discover-groups" [GENES] "with --truncation-abundance global or shrink, find the groups of cells that share isoform usage from the molecules with one candidate transcript (bisecting the cells with held-out likelihood tests, in the pass that fits the global parameters) in place of a --truncation-groups file; the groups are written to truncation_groups.tsv. GENES is a transcript<TAB>gene file the discovery takes its genes from, needed when the -m map does not group transcripts into genes (it does not change the loci the pooled prior acts within, unlike --truncation-genes); requires all cells in that pass (--truncation-global-cells 0)")
        .value_parser(value_parser!(String))
        .num_args(0..=1)
        .default_missing_value("")
        .conflicts_with("truncation-groups"))
    .arg(arg!(--"truncation-interval-mode" <MODE> "how multiple molecule intervals per transcript contribute to the truncation model: all, union or none (positions ignored)")
        .ignore_case(true)
        .value_parser(value_parser!(TruncationIntervalMode))
        .default_value("all"))
    .arg(arg!(--"truncation-params" <SCOPE> "where the truncation model's parameters are estimated: cell (separately in each cell) or global (once from all cells, then fixed in every cell's EM)")
        .ignore_case(true)
        .value_parser(value_parser!(TruncationParamScope))
        .default_value("cell"))
    .arg(arg!(--"truncation-length-law" <LAW> "how the truncation model's read-length laws depend on transcript length: marginal (one law) or conditional (one per transcript-length class)")
        .ignore_case(true)
        .value_parser(value_parser!(TruncationLengthLaw))
        .default_value("marginal"))
    .arg(arg!(--"truncation-positional-law" <LAW> "the truncation model's positional law: mixture (oarfish's four-component mixture) or four-param (the earlier four-parameter model, with the same EM, pooling and shrinkage options)")
        .ignore_case(true)
        .value_parser(value_parser!(TruncationPositionalLaw))
        .default_value("mixture"))
    .arg(arg!(--"truncation-end-kernels" <KIND> "the truncation model's end kernels: mixed (absolute and relative-to-length parts) or absolute (absolute part only)")
        .ignore_case(true)
        .value_parser(value_parser!(TruncationEndKernels))
        .default_value("mixed"))
    .arg(arg!(--"truncation-abundance" <MODE> "where each cell's transcript abundances come from in the truncation model: cell (the cell's own EM), global (pooled over all cells) or shrink (the cell's EM shrunk toward the pooled isoform fractions); global and shrink require --truncation-params global")
        .ignore_case(true)
        .value_parser(value_parser!(TruncationAbundance))
        .default_value("cell"))
    .arg(arg!(--"truncation-positional-power" <P> "power the truncation model's positional likelihood is raised to (1: as modelled, 0: positions ignored)")
        .value_parser(value_parser!(f64))
        .default_value("1.0"))
    .arg(arg!(--"truncation-calibrate-power" "with pooled abundances, choose the positional power from the data: the power that best predicts, from position and abundance, the transcript of molecules whose alignments settle it")
        .requires("truncation-abundance"))
    .arg(arg!(--"truncation-min-share" <EPS> "drop each molecule's final transcript shares below EPS and renormalize the rest (0: keep all)")
        .value_parser(value_parser!(f64))
        .default_value("0.0"))
    .arg(arg!(--"truncation-sample-below" <EPS> "give the combined mass of each molecule's final transcript shares below EPS to one of them, drawn by share (1: one transcript per molecule; 0: keep all)")
        .value_parser(value_parser!(f64))
        .default_value("0.0"))
    .arg(arg!(--"truncation-shrink" <K> "with --truncation-abundance shrink, the strength of the pooled prior: K pseudo-molecules per molecule the cell has in a locus, auto (a number of pseudo-molecules per locus estimated from the molecules with one candidate transcript) or auto-gene (the same per gene, shrunk toward the common value; needs --truncation-genes)")
        .ignore_case(true)
        .value_parser(value_parser!(TruncationShrink))
        .default_value("1.0"))
    .arg(arg!(--"truncation-groups" <FILE> "with --truncation-abundance global or shrink, pool abundances within groups of cells: a barcode<TAB>group file; unlisted barcodes form one more group")
        .value_parser(pathbuf_file_exists_validator))
    .arg(arg!(--"truncation-genes" <FILE> "with --truncation-abundance shrink, a transcript<TAB>gene file: the pooled prior then acts within genes instead of within loci (connected sets of candidate transcripts, which can merge genes), and --truncation-shrink auto estimates its strength from isoform shares within genes")
        .value_parser(pathbuf_file_exists_validator))
    .arg(arg!(--"truncation-group-shrink" <K> "with --truncation-groups, pull each group's pooled isoform shares toward those of all cells: K pseudo-molecules per locus; auto: the loci whose usage differs between the groups beyond sampling keep each group's own shares with a strength estimated from the groups' fits, the rest take the shares of all cells; auto-gene: the same with genes as the loci (needs --truncation-genes); 0: the groups are pooled independently")
        .ignore_case(true)
        .value_parser(value_parser!(TruncationShrink))
        .default_value("0"))
    .arg(arg!(--"truncation-dump" <DIR> "with --truncation-params global, write every cell's molecule counts per candidate transcript (barcode, transcript, molecules with that transcript as their only candidate, molecules with it among their candidates) to DIR/unique_counts.tsv.gz, for analyses of isoform usage across cells")
        .value_parser(value_parser!(PathBuf)))
    .arg(arg!(--"truncation-oracle" <FILE> "for evaluation on simulated data: a transcript<TAB>barcode<TAB>count file of true per-cell abundances, used as each listed cell's fixed abundances in place of the pooled ones (requires --truncation-abundance global; unlisted cells keep the pooled abundances)")
        .value_parser(pathbuf_file_exists_validator))
    .arg(arg!(--"truncation-global-cells" <N> "with --truncation-params global, estimate the parameters from about this many cells spread evenly through the input (0: all cells)")
        .value_parser(value_parser!(usize))
        .default_value("0"))
    .arg(arg!(--"coverage-bin-width" <BP> "transcript coverage bin width in bases")
        .value_parser(value_parser!(u32))
        .default_value("100"))
    .arg(arg!(--"coverage-growth-rate" <RATE> "logistic coverage-model growth rate")
        .value_parser(value_parser!(f64))
        .default_value("2.0"));

    let infer_app = Command::new("infer")
    .about("Perform inference on equivalence class count data")
    .version(version)
    .author(crate_authors)
    .arg(arg!(-c --"count-mat" <EQCMAT> "matrix of cells by equivalence class counts")
        .required(true)
        .value_parser(pathbuf_file_exists_validator))
    .arg(arg!(-e --"eq-labels" <EQLABELS> "file containing the gene labels of the equivalence classes")
        .required(true)
        .value_parser(pathbuf_file_exists_validator))
    .arg(arg!(-o --"output-dir" <OUTPUTDIR> "output directory where quantification results will be written").required(true).value_parser(value_parser!(PathBuf)))
    .arg(arg!(-t --threads <THREADS> "number of threads to use for processing").value_parser(value_parser!(u32)).default_value(max_num_threads))
    .arg(arg!(--usa "flag specifying that input equivalence classes were computed in USA mode"))
    .arg(arg!(--"quant-subset" <SFILE> "file containing list of barcodes to quantify, those not in this list will be ignored").value_parser(pathbuf_file_exists_validator))
    .arg(arg!(--"use-mtx" "flag for writing output matrix in matrix market format (default)"))
    .arg(arg!(--"use-eds" "flag for writing output matrix in EDS format").conflicts_with("use-mtx"));

    let atac_app = atac_sub_commands();

    let opts = Command::new("alevin-fry")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .version(version)
        .author(crate_authors)
        .about("Process RAD files from the command line")
        .subcommand(gen_app)
        .subcommand(collate_app)
        .subcommand(quant_app)
        .subcommand(infer_app)
        .subcommand(convert_app)
        .subcommand(view_app)
        .subcommand(atac_app)
        .get_matches();

    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator)
        .use_custom_timestamp(|out: &mut dyn std::io::Write| {
            write!(out, "{}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S")).unwrap();
            Ok(())
        })
        .build()
        .fuse();
    let drain = slog_async::Async::new(drain).build().fuse();
    let log = slog::Logger::root(drain, o!());

    use alevin_fry::atac;
    if let Some(t) = opts.subcommand_matches("atac") {
        info!(&log, "scATAC-seq mode");
        atac::run::run(t, VERSION, &cmdline, &log)?;
    }

    // You can handle information about subcommands by requesting their matches by name
    // (as below), requesting just the name used, or both at the same time
    if let Some(t) = opts.subcommand_matches("generate-permit-list") {
        let input_dir: &PathBuf = t.get_one("input").expect("no input directory specified");
        let output_dir: &PathBuf = t
            .get_one("output-dir")
            .expect("no output directory specified");

        let valid_ori: bool;
        let expected_ori = match t
            .get_one::<String>("expected-ori")
            .unwrap()
            .to_uppercase()
            .as_str()
        {
            "RC" => {
                valid_ori = true;
                Strand::Reverse
            }
            "FW" => {
                valid_ori = true;
                Strand::Forward
            }
            "BOTH" => {
                valid_ori = true;
                Strand::Unknown
            }
            "EITHER" => {
                valid_ori = true;
                Strand::Unknown
            }
            _ => {
                valid_ori = false;
                Strand::Unknown
            }
        };

        if !valid_ori {
            crit!(
                log,
                "{} is not a valid option for --expected-ori",
                expected_ori
            );
            std::process::exit(1);
        }

        let mut fmeth = CellFilterMethod::KneeFinding;

        let _expect_cells: Option<usize> = match t.get_one::<usize>("expect-cells") {
            Some(v) => {
                fmeth = CellFilterMethod::ExpectCells(*v);
                Some(*v)
            }
            None => None,
        };

        if t.get_flag("knee-distance") {
            fmeth = CellFilterMethod::KneeFinding;
        }

        let _force_cells = match t.get_one::<usize>("force-cells") {
            Some(v) => {
                fmeth = CellFilterMethod::ForceCells(*v);
                Some(*v)
            }
            None => None,
        };

        let _valid_bc = match t.get_one::<PathBuf>("valid-bc") {
            Some(v) => {
                fmeth = CellFilterMethod::ExplicitList(v.clone());
                Some(v)
            }
            None => None,
        };

        //let _unfiltered_pl = match t.get_one::<String>("unfiltered-pl") {
        if let Some(v) = t.get_one::<PathBuf>("unfiltered-pl") {
            let min_reads: usize = *t
                .get_one("min-reads")
                .expect("min-reads must be a valid integer");
            if min_reads < 1 {
                crit!(
                    log,
                    "min-reads < 1 is not supported, the value {} was provided",
                    min_reads
                );
                std::process::exit(1);
            }
            fmeth = CellFilterMethod::UnfilteredExternalList(v.clone(), min_reads);
        };

        // velo_mode --- currently, on this branch, it is always false
        let velo_mode = false; //t.get_flag("velocity-mode");
        let gpl_threads: usize = *t
            .get_one::<u32>("threads")
            .expect("valid integer number of threads") as usize;

        let gpl_opts = GenPermitListOpts::builder()
            .input_dir(input_dir)
            .output_dir(output_dir)
            .fmeth(fmeth)
            .expected_ori(expected_ori)
            .version(VERSION)
            .threads(gpl_threads)
            .velo_mode(velo_mode)
            .cmdline(&cmdline)
            .log(&log)
            .build();

        match generate_permit_list(gpl_opts) {
            Ok(0) => {
                warn!(log, "found 0 corrected barcodes; please check the input.");
            }
            Err(e) => return Err(e),
            _ => (),
        };
    }

    // convert a BAM file, in *transcriptomic coordinates*, with
    // the appropriate barcode and umi tags, into a RAD file
    if let Some(t) = opts.subcommand_matches("convert") {
        let input_file: &PathBuf = t.get_one("bam").unwrap();
        let rad_file: &PathBuf = t.get_one("output").unwrap();
        let num_threads: u32 = *t.get_one("threads").unwrap();
        let filter: bool = t.get_flag("filter_best");
        let score_threshold: f64 = *t.get_one("score-threshold").unwrap();

        alevin_fry::convert::bam2rad(
            input_file,
            rad_file,
            num_threads,
            filter,
            score_threshold,
            &log,
        )?
    }

    // convert a rad file to a textual representation and write to stdout
    if let Some(t) = opts.subcommand_matches("view") {
        let rad_file: &PathBuf = t.get_one("rad").unwrap();
        let print_header = t.get_flag("header");
        alevin_fry::convert::view(rad_file, print_header, &log)
    }

    // collate a rad file to group together all records corresponding
    // to the same corrected barcode.
    if let Some(t) = opts.subcommand_matches("collate") {
        let input_dir: &PathBuf = t.get_one("input-dir").unwrap();
        let rad_dir: &PathBuf = t.get_one("rad-dir").unwrap();
        let num_threads = *t.get_one("threads").unwrap();
        let compress_out = t.get_flag("compress");
        let max_records: u32 = *t.get_one("max-records").unwrap();
        alevin_fry::collate::collate(
            input_dir,
            rad_dir,
            num_threads,
            max_records,
            compress_out,
            &cmdline,
            VERSION,
            &log,
        )
        .expect("could not collate.");
    }

    // perform quantification of a collated rad file.
    if let Some(t) = opts.subcommand_matches("quant") {
        let num_threads = *t.get_one("threads").unwrap();
        let num_bootstraps = *t.get_one("num-bootstraps").unwrap();
        let init_uniform = t.get_flag("init-uniform");
        let summary_stat = t.get_flag("summary-stat");
        let dump_eq = t.get_flag("dump-eqclasses");
        let use_mtx = !t.get_flag("use-eds");
        let input_dir: &PathBuf = t.get_one("input-dir").unwrap();
        let output_dir: &PathBuf = t.get_one("output-dir").unwrap();
        let tg_map: &PathBuf = t.get_one("tg-map").unwrap();
        let junction_model = t.get_flag("junction-model");
        let junction_gtf: Option<&PathBuf> = t.get_one("junction-gtf");
        let junction_likelihood = t.get_flag("junction-likelihood");
        let resolution = *t.get_one::<ResolutionStrategy>("resolution").unwrap();
        let sa_model = *t.get_one::<SplicedAmbiguityModel>("sa-model").unwrap();
        let small_thresh = *t.get_one("small-thresh").unwrap();
        let lambda_size: f64 = *t.get_one("lambda-size").unwrap();
        let tau_delta: f64 = *t.get_one("tau-delta").unwrap();
        let filter_list: Option<&PathBuf> = t.get_one("quant-subset");
        let large_graph_thresh: usize = *t.get_one("large-graph-thresh").unwrap();
        let umi_edit_dist: u32 = *t.get_one("umi-edit-dist").unwrap();
        let true_umi_mode = *t.get_one::<TrueUmiOracleMode>("true-umi-mode").unwrap();
        let true_txp_mode: Option<bool> = t.get_one("true-txp-mode").cloned();
        let true_umi_txp_file: Option<&PathBuf> = t.get_one("true-umi-txp-file");
        let ed_distance_model = *t.get_one::<EditDistanceMode>("ed-distance-model").unwrap();
        let ed_upper_bound: usize = *t.get_one("ed-upper-bound").unwrap();
        let em_prob_model = *t.get_one::<EMMode>("em-prob-model").unwrap();
        let end_model = *t.get_one::<EndModel>("end-model").unwrap();
        let end_distribution = *t.get_one::<EndDistribution>("end-distribution").unwrap();
        let end_threshold: f64 = *t.get_one("end-threshold").unwrap();
        let end_scale: f64 = *t.get_one("end-scale").unwrap();
        let coverage_model = t.get_flag("coverage-model");
        let coverage_dedup_model = t.get_flag("coverage-dedup");
        let truncation_preset = t.get_one::<TruncationPreset>("truncation-preset").copied();
        let truncation_model = t.get_flag("model-truncation") || truncation_preset.is_some();
        // an optional value: the flag alone takes the genes from the -m map
        let truncation_discover_arg: Option<&String> = t.get_one("truncation-discover-groups");
        let mut truncation_discover_groups = truncation_discover_arg.is_some();
        let truncation_discover_genes: Option<PathBuf> =
            truncation_discover_arg.filter(|p| !p.is_empty()).map(PathBuf::from);
        if let Some(p) = truncation_discover_genes.as_ref() {
            if !p.is_file() {
                bail!("--truncation-discover-groups: {} is not a file", p.display());
            }
        }
        let truncation_interval_mode = *t
            .get_one::<TruncationIntervalMode>("truncation-interval-mode")
            .unwrap();
        let mut truncation_params = *t
            .get_one::<TruncationParamScope>("truncation-params")
            .unwrap();
        let truncation_global_cells: usize = *t.get_one("truncation-global-cells").unwrap();
        let truncation_groups: Option<&PathBuf> = t.get_one("truncation-groups");
        let truncation_genes: Option<&PathBuf> = t.get_one("truncation-genes");
        let truncation_dump: Option<&PathBuf> = t.get_one("truncation-dump");
        let truncation_oracle: Option<&PathBuf> = t.get_one("truncation-oracle");
        let (mut truncation_group_shrink, mut truncation_group_shrink_auto, truncation_group_shrink_per_gene) =
            match *t.get_one::<TruncationShrink>("truncation-group-shrink").unwrap() {
                TruncationShrink::Fixed(k) => (k, false, false),
                TruncationShrink::Auto => (0.0, true, false),
                TruncationShrink::AutoGene => (0.0, true, true),
            };
        let truncation_length_law = *t
            .get_one::<TruncationLengthLaw>("truncation-length-law")
            .unwrap();
        let truncation_positional_law = *t
            .get_one::<TruncationPositionalLaw>("truncation-positional-law")
            .unwrap();
        let mut truncation_end_kernels = *t
            .get_one::<TruncationEndKernels>("truncation-end-kernels")
            .unwrap();
        let mut truncation_abundance = *t
            .get_one::<TruncationAbundance>("truncation-abundance")
            .unwrap();
        let (truncation_shrink, truncation_shrink_auto, truncation_shrink_per_gene) =
            match *t.get_one::<TruncationShrink>("truncation-shrink").unwrap() {
                TruncationShrink::Fixed(k) => (k, false, false),
                TruncationShrink::Auto => (1.0, true, false),
                TruncationShrink::AutoGene => (1.0, true, true),
            };
        let mut truncation_positional_power: f64 = *t.get_one("truncation-positional-power").unwrap();
        let truncation_calibrate_power = t.get_flag("truncation-calibrate-power");
        // a preset sets every truncation option the command line left at its default
        if let Some(preset) = truncation_preset {
            let given = |id: &str| t.value_source(id) == Some(clap::parser::ValueSource::CommandLine);
            match preset {
                TruncationPreset::Recommended => {
                    if !given("truncation-params") {
                        truncation_params = TruncationParamScope::Global;
                    }
                    if !given("truncation-abundance") {
                        truncation_abundance = TruncationAbundance::Shrink;
                    }
                    if !given("truncation-end-kernels") {
                        truncation_end_kernels = TruncationEndKernels::Absolute;
                    }
                    if !given("truncation-positional-power") {
                        truncation_positional_power = 0.5;
                    }
                    if !given("truncation-group-shrink") {
                        truncation_group_shrink = 0.0;
                        truncation_group_shrink_auto = true;
                    }
                    if truncation_groups.is_none() && !given("truncation-discover-groups") {
                        truncation_discover_groups = true;
                    }
                }
            }
        }
        let truncation_min_share: f64 = *t.get_one("truncation-min-share").unwrap();
        let truncation_sample_below: f64 = *t.get_one("truncation-sample-below").unwrap();
        let coverage_bin_width: u32 = *t.get_one("coverage-bin-width").unwrap();
        let coverage_growth_rate: f64 = *t.get_one("coverage-growth-rate").unwrap();
        let coverage_source = *t.get_one::<CoverageSource>("coverage-source").unwrap();
        let coverage_enabled = coverage_model || coverage_dedup_model;

        if junction_model
            && matches!(
                resolution,
                ResolutionStrategy::ParsimonyGene | ResolutionStrategy::ParsimonyGeneEm
            )
        {
            bail!("--junction-model requires transcript-level resolution");
        }
        if junction_model && true_umi_mode == TrueUmiOracleMode::DirectCollapse {
            bail!("--junction-model is not supported with --true-umi-mode=direct");
        }

        if coverage_model && coverage_dedup_model {
            bail!("--coverage-model and --coverage-dedup are separate modes; select one");
        }
        if truncation_model && (end_model != EndModel::Off || num_bootstraps > 0) {
            bail!("--model-truncation cannot be combined with --end-model or bootstrapping");
        }
        if truncation_model && coverage_dedup_model && coverage_source == CoverageSource::Cell {
            bail!("--model-truncation with --coverage-dedup requires --coverage-source global: a cell's own coverage is only known after its deduplication");
        }
        if coverage_source == CoverageSource::Global && !(coverage_model || coverage_dedup_model) {
            bail!("--coverage-source global requires --coverage-model or --coverage-dedup");
        }
        if junction_likelihood && !(truncation_model && truncation_params == TruncationParamScope::Global) {
            bail!("--junction-likelihood requires --model-truncation --truncation-params global");
        }
        if junction_gtf.is_some() && !(junction_model || junction_likelihood) {
            bail!("--junction-gtf requires --junction-model or --junction-likelihood");
        }
        if truncation_model && resolution != ResolutionStrategy::ParsimonyEm {
            bail!("--model-truncation requires transcript-level --resolution parsimony-em");
        }
        if truncation_model && true_umi_mode == TrueUmiOracleMode::DirectCollapse {
            bail!("--model-truncation is not supported with --true-umi-mode=direct");
        }
        if !truncation_model && truncation_params == TruncationParamScope::Global {
            bail!("--truncation-params global requires --model-truncation");
        }
        if truncation_abundance != TruncationAbundance::Cell
            && truncation_params != TruncationParamScope::Global
        {
            bail!("--truncation-abundance global or shrink requires --truncation-params global");
        }
        if truncation_shrink_auto && truncation_abundance != TruncationAbundance::Shrink {
            bail!("--truncation-shrink auto requires --truncation-abundance shrink");
        }
        if truncation_shrink_per_gene && truncation_genes.is_none() {
            bail!("--truncation-shrink auto-gene requires --truncation-genes (without it the loci merge genes)");
        }
        if truncation_groups.is_some() && truncation_abundance == TruncationAbundance::Cell {
            bail!("--truncation-groups requires --truncation-abundance global or shrink");
        }
        if truncation_genes.is_some()
            && truncation_abundance != TruncationAbundance::Shrink
            && !truncation_group_shrink_auto
        {
            bail!("--truncation-genes requires --truncation-abundance shrink or --truncation-group-shrink auto");
        }
        if (truncation_group_shrink > 0.0 || truncation_group_shrink_auto)
            && truncation_groups.is_none()
            && !truncation_discover_groups
        {
            bail!("--truncation-group-shrink requires --truncation-groups or --truncation-discover-groups");
        }
        if truncation_discover_groups
            && !(truncation_model
                && truncation_params == TruncationParamScope::Global
                && truncation_abundance != TruncationAbundance::Cell)
        {
            bail!("--truncation-discover-groups requires --model-truncation --truncation-params global and --truncation-abundance global or shrink");
        }
        if truncation_discover_groups && truncation_global_cells > 0 {
            bail!("--truncation-discover-groups needs every cell in the global pass: it cannot be combined with --truncation-global-cells");
        }
        if truncation_discover_groups && truncation_oracle.is_some() {
            bail!("--truncation-oracle cannot be combined with --truncation-discover-groups");
        }
        if truncation_group_shrink_per_gene && truncation_genes.is_none() {
            bail!("--truncation-group-shrink auto-gene requires --truncation-genes");
        }
        if truncation_dump.is_some() && !(truncation_model && truncation_params == TruncationParamScope::Global) {
            bail!("--truncation-dump requires --model-truncation --truncation-params global");
        }
        if truncation_oracle.is_some()
            && !(truncation_model
                && truncation_params == TruncationParamScope::Global
                && truncation_abundance == TruncationAbundance::Global)
        {
            bail!("--truncation-oracle requires --model-truncation --truncation-params global --truncation-abundance global");
        }
        if truncation_oracle.is_some() && truncation_groups.is_some() {
            bail!("--truncation-oracle cannot be combined with --truncation-groups");
        }
        if truncation_global_cells > 0 && truncation_params != TruncationParamScope::Global {
            bail!("--truncation-global-cells requires --truncation-params global");
        }
        if coverage_enabled && (coverage_bin_width == 0 || coverage_growth_rate <= 0.0) {
            bail!("coverage bin width and growth rate must be positive");
        }
        if coverage_enabled && resolution != ResolutionStrategy::ParsimonyEm {
            bail!("coverage modes currently require transcript-level parsimony-em resolution");
        }
        if coverage_enabled && true_umi_mode == TrueUmiOracleMode::DirectCollapse {
            bail!("coverage modes are not supported with --true-umi-mode=direct");
        }
        if coverage_enabled && num_bootstraps > 0 {
            bail!("coverage modes cannot currently be combined with bootstrapping");
        }

        if end_threshold < 0.0 || end_scale <= 0.0 {
            bail!("--end-threshold must be non-negative and --end-scale must be positive");
        }

        if (true_umi_mode != TrueUmiOracleMode::Off || true_txp_mode == Some(true))
            && true_umi_txp_file.is_none()
        {
            crit!(
                log,
                "--true-umi-mode was set to {:?} but no --true-umi-txp-file was provided",
                true_umi_mode
            );
            bail!("Invalid command line option");
        }
        let mut pug_exact_umi = false;

        match umi_edit_dist {
            0 => {
                match resolution {
                    ResolutionStrategy::Trivial
                    | ResolutionStrategy::CellRangerLike
                    | ResolutionStrategy::CellRangerLikeEm => {
                        // already false, not pug_exact_umi because
                        // these methods don't use PUG
                    }
                    ResolutionStrategy::Parsimony
                    | ResolutionStrategy::ParsimonyEm
                    | ResolutionStrategy::ParsimonyGene
                    | ResolutionStrategy::ParsimonyGeneEm => {
                        pug_exact_umi = true;
                    }
                }
            }
            1 => {
                match resolution {
                    ResolutionStrategy::Trivial
                    | ResolutionStrategy::CellRangerLike
                    | ResolutionStrategy::CellRangerLikeEm => {
                        // these methods don't currently support 1 edit UMIs
                        crit!(
                            log,
                            "\n\nResolution strategy {:?} doesn't currently support 1-edit UMI resolution",
                            resolution
                        );
                        bail!("Invalid command line option");
                    }
                    ResolutionStrategy::Parsimony
                    | ResolutionStrategy::ParsimonyEm
                    | ResolutionStrategy::ParsimonyGene
                    | ResolutionStrategy::ParsimonyGeneEm => {
                        pug_exact_umi = false;
                    }
                }
            }
            j => {
                // no method currently supported edit distance 2 or greater correction
                crit!(
                    log,
                    "\n\nResolution strategy {:?} doesn't currently support {}-edit UMI resolution",
                    resolution,
                    j
                );
                bail!("Invalid command line option");
            }
        }

        if dump_eq && (resolution == ResolutionStrategy::Trivial) {
            crit!(
                log,
                "\n\nGene equivalence classes are not meaningful in case of Trivial resolution."
            );
            std::process::exit(1);
        }

        if num_bootstraps > 0 {
            match resolution {
                ResolutionStrategy::CellRangerLikeEm
                | ResolutionStrategy::ParsimonyEm
                | ResolutionStrategy::ParsimonyGeneEm => {
                    // sounds good
                }
                _ => {
                    eprintln!(
                        "\n\nThe num_bootstraps argument was set to {}, but bootstrapping can only be used with the cr-like-em, parsimony-em, or parsimony-gene-em resolution strategies",
                        num_bootstraps
                    );
                    std::process::exit(1);
                }
            }
        }

        // first make sure that the input direcory passed in has the
        // appropriate json file in it.
        // we should take care to document this workflow explicitly.
        let parent = std::path::Path::new(&input_dir);
        let json_path = parent.join("generate_permit_list.json");

        // build the QuantOpts structure
        let quant_opts = QuantOpts::builder()
            .input_dir(input_dir)
            .tg_map(tg_map)
            .output_dir(output_dir)
            .junction_model(junction_model)
            .junction_gtf(junction_gtf)
            .junction_likelihood(junction_likelihood)
            .num_threads(num_threads)
            .num_bootstraps(num_bootstraps)
            .init_uniform(init_uniform)
            .summary_stat(summary_stat)
            .dump_eq(dump_eq)
            .use_mtx(use_mtx)
            .resolution(resolution)
            .sa_model(sa_model)
            .small_thresh(small_thresh)
            .lambda_size(lambda_size)
            .tau_delta(tau_delta)
            .large_graph_thresh(large_graph_thresh)
            .filter_list(filter_list)
            .pug_exact_umi(pug_exact_umi)
            .true_umi_mode(true_umi_mode)
            .true_txp_mode(true_txp_mode)
            .true_umi_txp_file(true_umi_txp_file.cloned())
            .ed_distance_model(ed_distance_model)
            .ed_upper_bound(ed_upper_bound)
            .em_prob_model(em_prob_model)
            .end_model(end_model)
            .end_distribution(end_distribution)
            .end_threshold(end_threshold)
            .end_scale(end_scale)
            .coverage_model(coverage_model)
            .coverage_dedup_model(coverage_dedup_model)
            .truncation_model(truncation_model)
            .truncation_preset(truncation_preset)
            .truncation_discover_groups(truncation_discover_groups)
            .truncation_discover_genes(truncation_discover_genes)
            .truncation_interval_mode(truncation_interval_mode)
            .truncation_params(truncation_params)
            .truncation_length_law(truncation_length_law)
            .truncation_positional_law(truncation_positional_law)
            .truncation_end_kernels(truncation_end_kernels)
            .truncation_abundance(truncation_abundance)
            .truncation_shrink(truncation_shrink)
            .truncation_shrink_auto(truncation_shrink_auto)
            .truncation_shrink_per_gene(truncation_shrink_per_gene)
            .truncation_positional_power(truncation_positional_power)
            .truncation_calibrate_power(truncation_calibrate_power)
            .truncation_min_share(truncation_min_share)
            .truncation_sample_below(truncation_sample_below)
            .truncation_global_cells(truncation_global_cells)
            .truncation_groups(truncation_groups.cloned())
            .truncation_genes(truncation_genes.cloned())
            .truncation_dump(truncation_dump.cloned())
            .truncation_oracle(truncation_oracle.cloned())
            .truncation_group_shrink(truncation_group_shrink)
            .truncation_group_shrink_auto(truncation_group_shrink_auto)
            .truncation_group_shrink_per_gene(truncation_group_shrink_per_gene)
            .coverage_bin_width(coverage_bin_width)
            .coverage_growth_rate(coverage_growth_rate)
            .coverage_source(coverage_source)
            .cmdline(&cmdline)
            .version(VERSION)
            .log(&log)
            .build();

        // if the input directory contains the valid json file we want
        // then proceed.  otherwise print a critical error.
        if json_path.exists() {
            let velo_mode = alevin_fry::utils::is_velo_mode(quant_opts.input_dir);
            if velo_mode {
                if quant_opts.junction_model {
                    bail!("--junction-model is not supported for velocity quantification");
                }
                match alevin_fry::quant::velo_quantify(quant_opts) {
                    // if we're all good; then great!
                    Ok(_) => {}
                    // if we have an error, see if it's an error parsing
                    // the CSV or something else.
                    Err(e) => match e.downcast_ref::<CSVError>() {
                        Some(error) => {
                            match *error.kind() {
                                // if a deserialize error, we already complained about it
                                ErrorKind::Deserialize { .. } => {
                                    return Err(anyhow!("execution terminated unexpectedly"));
                                }
                                // if another type of error, just panic for now
                                _ => {
                                    panic!("could not quantify rad file: {e:#}");
                                }
                            }
                        }
                        // if something else, just panic
                        None => {
                            panic!("could not quantify rad file: {e:#}");
                        }
                    },
                }; // end match if
            } else {
                match alevin_fry::quant::quantify(quant_opts) {
                    // if we're all good; then great!
                    Ok(_) => {}
                    // if we have an error, see if it's an error parsing
                    // the CSV or something else.
                    Err(e) => match e.downcast_ref::<CSVError>() {
                        Some(error) => {
                            match *error.kind() {
                                // if a deserialize error, we already complained about it
                                ErrorKind::Deserialize { .. } => {
                                    return Err(anyhow!("execution terminated unexpectedly"));
                                }
                                // if another type of error, just panic for now
                                _ => {
                                    panic!("could not quantify rad file: {e:#}");
                                }
                            }
                        }
                        // if something else, just panic
                        None => {
                            panic!("could not quantify rad file: {e:#}");
                        }
                    },
                }; //end quant if
            }; // end velo_mode if
        } else {
            crit!(
                log,
                "The provided input directory lacks a generate_permit_list.json file; this should not happen."
            );
            bail!(
                "The provided input directory lacks a generate_permit_list.json file; this should not happen."
            );
        }
    } // end quant if

    // Given an input of equivalence class counts, perform inference
    // and output a target-by-cell count matrix.
    if let Some(t) = opts.subcommand_matches("infer") {
        let num_threads = *t.get_one("threads").unwrap();
        let use_mtx = !t.get_flag("use-eds");
        let output_dir = t.get_one("output-dir").unwrap();
        let count_mat = t.get_one("count-mat").unwrap();
        let eq_label_file = t.get_one("eq-labels").unwrap();
        let filter_list: Option<&PathBuf> = t.get_one("quant-subset");
        let usa_mode = t.get_flag("usa");

        alevin_fry::infer::infer(
            count_mat,
            eq_label_file,
            usa_mode,
            use_mtx,
            num_threads,
            filter_list,
            output_dir,
            &log,
        )
        .expect("could not perform inference from equivalence class counts.");
    }
    Ok(())
}

fn atac_sub_commands() -> Command {
    let num_hardware_threads = num_cpus::get() as u32;
    let max_num_threads: String = (num_cpus::get() as u32).to_string();
    let max_num_collate_threads: String = (16_u32.min(num_hardware_threads).max(2_u32)).to_string();
    let max_num_gpl_threads: String = (8_u32.min(num_hardware_threads).max(2_u32)).to_string();
    let max_num_sort_threads: String = (16_u32.min(num_hardware_threads).max(2_u32)).to_string();

    let crate_authors = crate_authors!("\n");
    let version = crate_version!();

    let gen_app = Command::new("generate-permit-list")
        .about("Generate a permit list of barcodes from a whitelist file")
        .version(version)
        .author(crate_authors)
        .arg(arg!(-i --input <INPUT>  "input directory containing the map.rad file")
            .required(true)
            .value_parser(pathbuf_directory_exists_validator))
        .arg(arg!(-o --"output-dir" <OUTPUTDIR>  "output directory")
            .required(true)
            .value_parser(value_parser!(PathBuf))
        )
        .arg(arg!(-t --threads <THREADS> "number of threads to use for the first phase of permit-list generation").value_parser(value_parser!(u32)).default_value(max_num_gpl_threads))
        .arg(
            arg!(-u --"unfiltered-pl" <UNFILTEREDPL> "uses an unfiltered external permit list")
                .value_parser(pathbuf_file_exists_validator)
        )
        .group(ArgGroup::new("filter-method")
            .args(["unfiltered-pl"])
            .required(true)
        )
        .arg(
            arg!(-m --"min-reads" <MINREADS> "minimum read count threshold; only used with --unfiltered-pl")
                .value_parser(value_parser!(usize))
                .default_value("10"))
        .arg(arg!(-d --"permit-bc-ori" <EXPECTEDORI> "the expected orientation of barcodes in the permit list")
             .ignore_case(true)
             .default_value("rc")
             .value_parser(["fw", "rc"]));

    let collate_app = Command::new("collate")
        .about("Collate a RAD file with corrected cell barcode")
        .version(version)
        .author(crate_authors)
        .arg(arg!(-i --"input-dir" <INPUTDIR> "output directory made by generate-permit-list")
            .required(true)
            .value_parser(pathbuf_directory_exists_validator))
        .arg(arg!(-r --"rad-dir" <RADDIR> "the directory containing the map.rad file which will be collated (typically produced as an output of the mapping)")
            .required(true)
            .value_parser(pathbuf_directory_exists_validator))
        .arg(arg!(-t --threads <THREADS> "number of threads to use for processing").value_parser(value_parser!(u32)).default_value(max_num_collate_threads.clone()))
        .arg(arg!(-c --compress "compress the output collated RAD file"))
        .arg(arg!(-m --"max-records" <MAXRECORDS> "the maximum number of read records to keep in memory at once")
            .value_parser(value_parser!(u32))
            .default_value("30000000"));

    let sort_app = Command::new("sort")
        .about("Produce coordinate sorted bed file")
        .version(version)
        .author(crate_authors)
        .arg(arg!(-i --"input-dir" <INPUTDIR> "output directory made by generate-permit-list")
            .required(true)
            .value_parser(pathbuf_directory_exists_validator))
        .arg(arg!(-r --"rad-dir" <RADDIR> "the directory containing the map.rad file which will be sorted (typically produced as an output of the mapping)")
            .required(true)
            .value_parser(pathbuf_directory_exists_validator))
        .arg(arg!(-t --threads <THREADS> "number of threads to use for processing").value_parser(value_parser!(u32)).default_value(max_num_sort_threads))
        .arg(arg!(-c --compress "compress the output of the sorted RAD file"))
        .arg(arg!(-m --"max-records" <MAXRECORDS> "the maximum number of read records to keep in memory at once")
            .value_parser(value_parser!(u32))
            .default_value("30000000"));

    let deduplicate_app = Command::new("deduplicate")
        .about("Deduplicate the RAD file and output a BED file")
        .version(version)
        .author(crate_authors)
        .arg(arg!(-i --"input-dir" <INPUTDIR> "input directory made by generate-permit-list that also contains the output of collate")
            .required(true)
            .value_parser(pathbuf_directory_exists_validator))
        .arg(arg!(-t --threads <THREADS> "number of threads to use for processing").value_parser(value_parser!(u32)).default_value(max_num_threads))
        .arg(arg!(-d --"permit-bc-ori" <EXPECTEDORI> "the expected orientation of barcodes in the permit list")
             .ignore_case(true)
             .default_value("rc")
             .value_parser(["fw", "rc"]));

    Command::new("atac")
        .about("subcommand for processing scATAC-seq RAD files")
        .version(version)
        .author(crate_authors)
        .subcommand(gen_app)
        .subcommand(sort_app)
        .subcommand(collate_app)
        .subcommand(deduplicate_app)
}
