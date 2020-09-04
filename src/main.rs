// Copyright 2020 Rob Patro, Avi Srivastava. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

extern crate bio_types;
extern crate chrono;
extern crate clap;
extern crate num_cpus;
extern crate rand;
extern crate slog;
extern crate slog_term;

use bio_types::strand::Strand;
use clap::{crate_authors, crate_version, App, Arg};
use libradicl::cellfilter::{generate_permit_list, CellFilterMethod};
use libradicl::schema::ResolutionStrategy;
use mimalloc::MiMalloc;
use rand::Rng;
use slog::{crit, o, warn, Drain};
use std::unimplemented;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

// static VERSION: &str = "0.0.1";
// static AUTHORS: &str = "Avi Srivastava, Rob Patro";

#[allow(dead_code)]
fn gen_random_kmer(k: usize) -> String {
    const CHARSET: &[u8] = b"ACGT";
    let mut rng = rand::thread_rng();
    let s: String = (0..k)
        .map(|_| {
            let idx = rng.gen_range(0, CHARSET.len());
            CHARSET[idx] as char
        })
        .collect();
    s
}

fn main() {
    let max_num_threads: String = (num_cpus::get() as u32).to_string();
    let crate_authors = crate_authors!("\n");
    let version = crate_version!();
    // [] add command for just counting barcode frequency
    // [] add other algorithms for determining barcode cutoff

    let convert_app = App::new("convert")
        .about("Convert a BAM file to a RAD file")
        .version(version)
        .author(crate_authors)
        .arg(Arg::from("-b, --bam=<bam-file> 'input SAM/BAM file'"))
        .arg(Arg::from("-t, --threads 'number of threads to use for processing'").default_value(&max_num_threads))
        .arg(Arg::from("-o, --output=<rad-file> 'output RAD file'"));

    let gen_app = App::new("generate-permit-list")
        .about("Generate a permit list of barcodes from a RAD file")
        .version(version)
        .author(crate_authors)
        .arg(Arg::from("-i, --input=<input>  'input RAD file'"))
        .arg(Arg::from("-d, --expected-ori=<expected-ori> 'the expected orientation of alignments'"))
        .arg(Arg::from(
            "-o, --output-dir=<output-dir>  'output directory'",
        ))
        .arg(Arg::from(
            "-k, --knee-distance  'attempt to determine the number of barcodes to keep using the knee distance method."
            ).conflicts_with_all(&["force-cell", "valid-bc", "expect-cells"])
        )
        .arg(Arg::from(
            "-e, --expect-cells=<expect-cells> 'defines the expected number of cells to use in determining the (read, not UMI) based cutoff'",
        ).conflicts_with_all(&["force-cells", "valid-bc", "knee-distance"])
        )
        .arg(Arg::from(
            "-f, --force-cells=<force-cells>  'select the top-k most-frequent barcodes, based on read count, as valid (true)'"
        ).conflicts_with_all(&["expect-cells", "valid-bc", "knee-distance"])
        )
        .arg(
            Arg::from(
                "-b, --valid-bc=<valid-bc> 'uses true barcode collected from a provided file'",
            )
            .conflicts_with_all(&["force-cells", "expect-cells", "knee-distance"]),
        );

    let collate_app = App::new("collate")
    .about("Collate a RAD file by corrected cell barcode")
    .version(version)
    .author(crate_authors)
    .arg(Arg::from("-i, --input-dir=<input-dir> 'input directory made by generate-permit-list'"))
    .arg(Arg::from("-r, --rad-file=<rad-file> 'the RAD file to be collated'"))
    .arg(Arg::from("-t, --threads 'number of threads to use for processing'").default_value(&max_num_threads))
    .arg(Arg::from("-m, --max-records=[max-records] 'the maximum number of read records to keep in memory at once'")
         .default_value("50000000"));
    //.arg(Arg::from("-e, --expected-ori=[expected-ori] 'the expected orientation of alignments'")
    //     .default_value("fw"));

    let quant_app = App::new("quant")
    .about("Quantify expression from a collated RAD file")
    .version(version)
    .author(crate_authors)
    .arg(Arg::from("-i, --input-dir=<input-dir>  'input directory containing collated RAD file'"))
    .arg(Arg::from("-m, --tg-map=<tg-map>  'transcript to gene map'"))
    .arg(Arg::from("-o, --output-dir=<output-dir> 'output directory where quantification results will be written'"))
    .arg(Arg::from("-t, --threads 'number of threads to use for processing'").default_value(&max_num_threads))
    .arg(Arg::from("-b, --num-bootstraps 'number of bootstraps to use'").default_value("0"))
    .arg(Arg::from("--init-uniform 'flag for uniform sampling'").requires("num-bootstraps").takes_value(false).required(false))
    .arg(Arg::from("--summary-stat 'flag for storing only summary statistics'").requires("num-bootstraps").takes_value(false).required(false))
    .arg(Arg::from("--use-mtx 'flag for writing output matrix in matrix market instead of EDS'").takes_value(false).required(false))
    .arg(Arg::from("-r, --resolution 'the resolution strategy by which molecules will be counted'")
        .possible_values(&["full", "trivial", "cr-like", "cr-like-em", "parsimony"])
        .default_value("full")
        .case_insensitive(true)
        .about("the resolution strategy by which molecules will be counted"));

    let opts = App::new("alevin-fry")
        .version(version)
        .author(crate_authors)
        .about("Process RAD files from the command line")
        .subcommand(gen_app)
        .subcommand(collate_app)
        .subcommand(quant_app)
        .subcommand(convert_app)
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
    if let Some(ref t) = opts.subcommand_matches("generate-permit-list") {
        let input_file: String = t.value_of_t("input").expect("no input directory specified");
        let output_dir: String = t
            .value_of_t("output-dir")
            .expect("no input directory specified");

        let valid_ori: bool;
        let expected_ori = match t.value_of("expected-ori").unwrap().to_uppercase().as_str() {
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

        let _expect_cells: Option<usize> = match t.value_of_t("expect-cells") {
            Ok(v) => {
                fmeth = CellFilterMethod::ExpectCells(v);
                Some(v)
            }
            Err(_) => None,
        };

        if t.is_present("knee-distance") {
            fmeth = CellFilterMethod::KneeFinding;
        }

        let _force_cells = match t.value_of_t("force-cells") {
            Ok(v) => {
                fmeth = CellFilterMethod::ForceCells(v);
                Some(v)
            }
            Err(_) => None,
        };

        let _valid_bc = match t.value_of_t::<String>("valid-bc") {
            Ok(v) => {
                fmeth = CellFilterMethod::ExplicitList(v.clone());
                Some(v)
            }
            Err(_) => None,
        };
        let nc = generate_permit_list(input_file, output_dir, fmeth, expected_ori, &log).unwrap();
        if nc == 0 {
            warn!(log, "found 0 corrected barcodes; please check the input.");
        }
    }

    if let Some(ref t) = opts.subcommand_matches("convert") {
        let input_file: String = t.value_of_t("bam").unwrap();
        let rad_file: String = t.value_of_t("output").unwrap();
        let num_threads: u32 = t.value_of_t("threads").unwrap();
        libradicl::convert::bam2rad(input_file, rad_file, num_threads, &log)
    }

    if let Some(ref t) = opts.subcommand_matches("collate") {
        let input_dir: String = t.value_of_t("input-dir").unwrap();
        let rad_file: String = t.value_of_t("rad-file").unwrap();
        let num_threads = t.value_of_t("threads").unwrap();
        let max_records: u32 = t.value_of_t("max-records").unwrap();
        libradicl::collate::collate(input_dir, rad_file, num_threads, max_records, &log)
            .expect("could not collate.");
    }

    if let Some(ref t) = opts.subcommand_matches("quant") {
        let num_threads = t.value_of_t("threads").unwrap();
        let num_bootstraps = t.value_of_t("num-bootstraps").unwrap();
        let init_uniform = t.is_present("init-uniform");
        let summary_stat = t.is_present("summary-stat");
        let use_mtx = t.is_present("use-mtx");
        let input_dir = t.value_of_t("input-dir").unwrap();
        let output_dir = t.value_of_t("output-dir").unwrap();
        let tg_map = t.value_of_t("tg-map").unwrap();
        let resolution: ResolutionStrategy = t.value_of_t("resolution").unwrap();
        libradicl::quant::quantify(
            input_dir,
            tg_map,
            output_dir,
            num_threads,
            num_bootstraps,
            init_uniform,
            summary_stat,
            use_mtx,
            resolution,
            &log,
        )
        .expect("could not quantify rad file.");
    }
}
