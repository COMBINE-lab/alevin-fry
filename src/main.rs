/*
 * Copyright (c) 2020-2022 Rob Patro, Avi Srivastava, Hirak Sarkar, Dongze He, Mohsen Zakeri.
 *
 * This file is part of alevin-fry
 * (see https://github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use anyhow::{anyhow, bail};
use bio_types::strand::Strand;
use clap::{arg, builder::ArgGroup, crate_authors, crate_version, value_parser, Command};
use csv::Error as CSVError;
use csv::ErrorKind;
use itertools::Itertools;
use mimalloc::MiMalloc;
use rand::Rng;
use slog::{crit, o, warn, Drain};
use std::path::PathBuf;

use alevin_fry::cellfilter::{generate_permit_list, CellFilterMethod};
use alevin_fry::cmd_parse_utils::{
    pathbuf_directory_exists_validator, pathbuf_file_exists_validator,
};
use alevin_fry::prog_opts::{GenPermitListOpts, QuantOpts};
use alevin_fry::quant::{ResolutionStrategy, SplicedAmbiguityModel};

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

// grab the version from the Cargo file.
const VERSION: &str = env!("CARGO_PKG_VERSION");

#[allow(dead_code)]
fn gen_random_kmer(k: usize) -> String {
    const CHARSET: &[u8] = b"ACGT";
    let mut rng = rand::thread_rng();
    let s: String = (0..k)
        .map(|_| {
            let idx = rng.gen_range(0..CHARSET.len());
            CHARSET[idx] as char
        })
        .collect();
    s
}

fn main() -> anyhow::Result<()> {
    let num_hardware_threads = num_cpus::get() as u32;
    let max_num_threads: String = (num_cpus::get() as u32).to_string();
    let max_num_collate_threads: String = (16_u32.min(num_hardware_threads).max(2_u32)).to_string();

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
                .value_parser(pathbuf_file_exists_validator),
        )
        .arg(
            arg!(-t --threads [THREADS] "number of threads to use for processing")
                .value_parser(value_parser!(u32))
                .default_value(&max_num_threads),
        )
        .arg(arg!(-o --output <RADFILE> "output RAD file").value_parser(value_parser!(PathBuf)));

    let view_app = Command::new("view")
        .about("View a RAD file")
        .version(version)
        .author(crate_authors)
        .arg(arg!(-r --rad <RADFILE> "input RAD file").value_parser(pathbuf_file_exists_validator))
        .arg(
            arg!(-H --header "flag for printing header")
                .takes_value(false)
                .required(false),
        )
        .arg(arg!(-o --output [RADFILE] "output plain-text-file file"));

    let gen_app = Command::new("generate-permit-list")
        .about("Generate a permit list of barcodes from a RAD file")
        .version(version)
        .author(crate_authors)
        .arg(arg!(-i --input <INPUT>  "input directory containing the map.rad RAD file")
            .value_parser(pathbuf_directory_exists_validator))
        .arg(arg!(-d --"expected-ori" <EXPECTEDORI> "the expected orientation of alignments")
             .ignore_case(true)
             .value_parser(["fw", "rc", "both", "either"]))
        .arg(arg!(-o --"output-dir" <OUTPUTDIR>  "output directory").value_parser(value_parser!(PathBuf)))
        .arg(arg!(
            -k --"knee-distance"  "attempt to determine the number of barcodes to keep using the knee distance method."
            ).conflicts_with_all(&["force-cell", "valid-bc", "expect-cells", "unfiltered-pl"])
        )
        .arg(arg!(-e --"expect-cells" <EXPECTCELLS> "defines the expected number of cells to use in determining the (read, not UMI) based cutoff")
             .conflicts_with_all(&["force-cells", "valid-bc", "knee-distance", "unfiltered-pl"])
        )
        .arg(arg!(-f --"force-cells" <FORCECELLS>  "select the top-k most-frequent barcodes, based on read count, as valid (true)")
             .conflicts_with_all(&["expect-cells", "valid-bc", "knee-distance", "unfiltered-pl"])
        )
        .arg(
            arg!(-b --"valid-bc" <VALIDBC> "uses true barcode collected from a provided file")
            .conflicts_with_all(&["force-cells", "expect-cells", "knee-distance", "unfiltered-pl"])
            .value_parser(pathbuf_file_exists_validator)
        )
        .arg(
            arg!(-u --"unfiltered-pl" <UNFILTEREDPL> "uses an unfiltered external permit list")
            .conflicts_with_all(&["force-cells", "expect-cells", "knee-distance", "valid-bc"])
            .value_parser(pathbuf_file_exists_validator)
        )
        .group(ArgGroup::new("filter-method")
               .args(&["knee-distance", "expect-cells", "force-cells", "valid-bc", "unfiltered-pl"])
               .required(true)
               )
        .arg(
            arg!(-m --"min-reads" [MINREADS] "minimum read count threshold; only used with --unfiltered-pl")
                .value_parser(value_parser!(usize))
                .default_value("10"));

    let collate_app = Command::new("collate")
    .about("Collate a RAD file by corrected cell barcode")
    .version(version)
    .author(crate_authors)
    .arg(arg!(-i --"input-dir" <INPUTDIR> "input directory made by generate-permit-list")
        .value_parser(pathbuf_directory_exists_validator))
    .arg(arg!(-r --"rad-dir" <RADFILE> "the directory containing the RAD file to be collated")
        .value_parser(pathbuf_directory_exists_validator))
    .arg(arg!(-t --threads [THREADS] "number of threads to use for processing").value_parser(value_parser!(u32)).default_value(&max_num_collate_threads))
    .arg(arg!(-c --compress "compress the output collated RAD file").takes_value(false).required(false))
    .arg(arg!(-m --"max-records" [MAXRECORDS] "the maximum number of read records to keep in memory at once")
         .value_parser(value_parser!(u32))
         .default_value("30000000"));

    let quant_app = Command::new("quant")
    .about("Quantify expression from a collated RAD file")
    .version(version)
    .author(crate_authors)
    .arg(arg!(-i --"input-dir" <INPUTDIR>  "input directory containing collated RAD file")
        .value_parser(pathbuf_directory_exists_validator))
    .arg(arg!(-m --"tg-map" <TGMAP>  "transcript to gene map").value_parser(pathbuf_file_exists_validator))
    .arg(arg!(-o --"output-dir" <OUTPUTDIR> "output directory where quantification results will be written").value_parser(value_parser!(PathBuf)))
    .arg(arg!(-t --threads [THREADS] "number of threads to use for processing").value_parser(value_parser!(u32)).default_value(&max_num_threads))
    .arg(arg!(-d --"dump-eqclasses" "flag for dumping equivalence classes").takes_value(false).required(false))
    .arg(arg!(-b --"num-bootstraps" [NUMBOOTSTRAPS] "number of bootstraps to use").value_parser(value_parser!(u32)).default_value("0"))
    .arg(arg!(--"init-uniform" "flag for uniform sampling").requires("num-bootstraps").takes_value(false).required(false))
    .arg(arg!(--"summary-stat" "flag for storing only summary statistics").requires("num-bootstraps").takes_value(false).required(false))
    .arg(arg!(--"use-mtx" "flag for writing output matrix in matrix market format (default)").takes_value(false).required(false))
    .arg(arg!(--"use-eds" "flag for writing output matrix in EDS format").takes_value(false).required(false).conflicts_with("use-mtx"))
    .arg(arg!(--"quant-subset" [SFILE] "file containing list of barcodes to quantify, those not in this list will be ignored").value_parser(pathbuf_file_exists_validator))
    .arg(arg!(-r --resolution <RESOLUTION> "the resolution strategy by which molecules will be counted")
        .ignore_case(true)
        .value_parser(value_parser!(ResolutionStrategy)))
    .arg(arg!(--"sa-model" [SAMODEL] "preferred model of splicing ambiguity")
        .ignore_case(true)
        .value_parser(value_parser!(SplicedAmbiguityModel))
        .default_value("winner-take-all")
        .hide(true))
    .arg(arg!(--"umi-edit-dist" [EDIST] "the Hamming distance within which potentially colliding UMIs will be considered for correction")
        .value_parser(value_parser!(u32))
        .default_value_ifs(&[
            ("resolution", Some("cr-like"), Some("0")),
            ("resolution", Some("cr-like-em"), Some("0")),
            ("resolution", Some("trivial"), Some("0")),
            ("resolution", Some("parsimony"), Some("1")),
            ("resolution", Some("parsimony-em"), Some("1")),
            ("resolution", Some("parsimony-gene"), Some("1")),
            ("resolution", Some("parsimony-gene-em"), Some("1")),
        ])
        .hide(true))
    .arg(arg!(--"large-graph-thresh" [NVERT] "the order (number of nodes) of a PUG above which the alternative resolution strategy will be applied")
        .value_parser(value_parser!(usize))
        .default_value_ifs(&[
            ("resolution", Some("parsimony-gene-em"), Some("1000")),
            ("resolution", Some("parsimony"), Some("1000")),
            ("resolution", Some("parsimony-em"), Some("1000")),
            ("resolution", Some("parsimony-gene"), Some("1000")),
            ("resolution", Some("parsimony-gene"), Some("1000")),
        ])
        .default_value("0") // for any other mode
        .hide(true))
    .arg(arg!(--"small-thresh" [SMALLTHRESH] "cells with fewer than these many reads will be resolved using a custom approach")
        .value_parser(value_parser!(usize))
        .default_value("10")
        .hide(true));

    let infer_app = Command::new("infer")
    .about("Perform inference on equivalence class count data")
    .version(version)
    .author(crate_authors)
    .arg(arg!(-c --"count-mat" <EQCMAT> "matrix of cells by equivalence class counts")
        .value_parser(pathbuf_file_exists_validator).takes_value(true))
    //.arg(arg!(-b --barcodes=<barcodes> "file containing the barcodes labeling the matrix rows").takes_value(true).required(true))
    .arg(arg!(-e --"eq-labels" <EQLABELS> "file containing the gene labels of the equivalence classes")
        .value_parser(pathbuf_file_exists_validator).takes_value(true))
    .arg(arg!(-o --"output-dir" <OUTPUTDIR> "output directory where quantification results will be written").value_parser(value_parser!(PathBuf)).takes_value(true))
    .arg(arg!(-t --threads [THREADS] "number of threads to use for processing").value_parser(value_parser!(u32)).default_value(&max_num_threads))
    .arg(arg!(--usa "flag specifying that input equivalence classes were computed in USA mode").takes_value(false).required(false))
    .arg(arg!(--"quant-subset" [SFILE] "file containing list of barcodes to quantify, those not in this list will be ignored").value_parser(pathbuf_file_exists_validator))
    .arg(arg!(--"use-mtx" "flag for writing output matrix in matrix market format (default)").takes_value(false).required(false))
    .arg(arg!(--"use-eds" "flag for writing output matrix in EDS format").takes_value(false).required(false).conflicts_with("use-mtx"));

    let opts = Command::new("alevin-fry")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .version(version)
        .author(crate_authors)
        .about("Process RAD files from the command line")
        .subcommand(gen_app)
        //.subcommand(test_app)
        .subcommand(collate_app)
        .subcommand(quant_app)
        .subcommand(infer_app)
        .subcommand(convert_app)
        .subcommand(view_app)
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

        if t.is_present("knee-distance") {
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
        let velo_mode = false; //t.is_present("velocity-mode");

        let gpl_opts = GenPermitListOpts::builder()
            .input_dir(input_dir)
            .output_dir(output_dir)
            .fmeth(fmeth)
            .expected_ori(expected_ori)
            .version(VERSION)
            .velo_mode(velo_mode)
            .cmdline(&cmdline)
            .log(&log)
            .build();

        match generate_permit_list(gpl_opts) {
            Ok(nc) if nc == 0 => {
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
        alevin_fry::convert::bam2rad(input_file, rad_file, num_threads, &log)
    }

    // convert a rad file to a textual representation and write to stdout
    if let Some(t) = opts.subcommand_matches("view") {
        let rad_file: &PathBuf = t.get_one("rad").unwrap();
        let print_header = t.is_present("header");
        let mut out_file: String = String::from("");
        if t.is_present("output") {
            out_file = t.get_one::<String>("output").unwrap().clone();
        }
        alevin_fry::convert::view(rad_file, print_header, out_file, &log)
    }

    // collate a rad file to group together all records corresponding
    // to the same corrected barcode.
    if let Some(t) = opts.subcommand_matches("collate") {
        let input_dir: &PathBuf = t.get_one("input-dir").unwrap();
        let rad_dir: &PathBuf = t.get_one("rad-dir").unwrap();
        let num_threads = *t.get_one("threads").unwrap();
        let compress_out = t.is_present("compress");
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
        let init_uniform = t.is_present("init-uniform");
        let summary_stat = t.is_present("summary-stat");
        let dump_eq = t.is_present("dump-eqclasses");
        let use_mtx = !t.is_present("use-eds");
        let input_dir: &PathBuf = t.get_one("input-dir").unwrap();
        let output_dir: &PathBuf = t.get_one("output-dir").unwrap();
        let tg_map: &PathBuf = t.get_one("tg-map").unwrap();
        let resolution = t
            .get_one::<ResolutionStrategy>("resolution")
            .unwrap()
            .clone();
        let sa_model = t
            .get_one::<SplicedAmbiguityModel>("sa-model")
            .unwrap()
            .clone();
        let small_thresh = *t.get_one("small-thresh").unwrap();
        let filter_list: Option<&PathBuf> = t.get_one("quant-subset");
        let large_graph_thresh: usize = *t.get_one("large-graph-thresh").unwrap();
        let umi_edit_dist: u32 = *t.get_one("umi-edit-dist").unwrap();
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
            .num_threads(num_threads)
            .num_bootstraps(num_bootstraps)
            .init_uniform(init_uniform)
            .summary_stat(summary_stat)
            .dump_eq(dump_eq)
            .use_mtx(use_mtx)
            .resolution(resolution)
            .sa_model(sa_model)
            .small_thresh(small_thresh)
            .large_graph_thresh(large_graph_thresh)
            .filter_list(filter_list)
            .pug_exact_umi(pug_exact_umi)
            .cmdline(&cmdline)
            .version(VERSION)
            .log(&log)
            .build();

        // if the input directory contains the valid json file we want
        // then proceed.  otherwise print a critical error.
        if json_path.exists() {
            let velo_mode = alevin_fry::utils::is_velo_mode(quant_opts.input_dir);
            if velo_mode {
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
                                    panic!("could not quantify rad file.");
                                }
                            }
                        }
                        // if something else, just panic
                        None => {
                            panic!("could not quantify rad file.");
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
                                    panic!("could not quantify rad file.");
                                }
                            }
                        }
                        // if something else, just panic
                        None => {
                            panic!("could not quantify rad file.");
                        }
                    },
                }; //end quant if
            }; // end velo_mode if
        } else {
            crit!(log,
            "The provided input directory lacks a generate_permit_list.json file; this should not happen."
           );
        }
    } // end quant if

    // Given an input of equivalence class counts, perform inference
    // and output a target-by-cell count matrix.
    if let Some(t) = opts.subcommand_matches("infer") {
        let num_threads = *t.get_one("threads").unwrap();
        let use_mtx = !t.is_present("use-eds");
        let output_dir = t.get_one("output-dir").unwrap();
        let count_mat = t.get_one("count-mat").unwrap();
        let eq_label_file = t.get_one("eq-labels").unwrap();
        let filter_list: Option<&PathBuf> = t.get_one("quant-subset");
        let usa_mode = t.is_present("usa");

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
