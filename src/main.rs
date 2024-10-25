use anyhow::bail;
use clap::{arg, builder::ArgGroup, crate_authors, crate_version, value_parser, Command};
use itertools::Itertools;
use slog::{crit, o, warn, Drain};
use std::path::PathBuf;

use piscem_atac::cellfilter::{generate_permit_list, CellFilterMethod};
use piscem_atac::cmd_parse_utils::{
    pathbuf_directory_exists_validator, pathbuf_file_exists_validator,
};
use piscem_atac::collate::collate;
use piscem_atac::deduplicate::deduplicate;
use piscem_atac::prog_opts::{DeduplicateOpts, GenPermitListOpts};

const VERSION: &str = env!("CARGO_PKG_VERSION");

fn main() -> anyhow::Result<()> {
    let num_hardware_threads = num_cpus::get() as u32;
    let max_num_threads: String = (num_cpus::get() as u32).to_string();
    let max_num_collate_threads: String = (16_u32.min(num_hardware_threads).max(2_u32)).to_string();

    let crate_authors = crate_authors!("\n");
    let version = crate_version!();
    let cmdline = std::env::args().join(" ");

    let gen_app = Command::new("generate-permit-list")
        .about("Generate a permit list of barcodes from a whitelist file")
        .version(version)
        .author(crate_authors)
        .arg(arg!(-i --input <INPUT>  "input directory containing the map.bed BED file")
            .required(true)
            .value_parser(pathbuf_directory_exists_validator))
        .arg(arg!(-o --"output-dir" <OUTPUTDIR>  "output directory")
            .required(true)
            .value_parser(value_parser!(PathBuf))
        )
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
        .arg(
            arg!(-r --"rev-comp" <REVERSECOMPLEMENT> "reverse complement")
            .value_parser(clap::builder::BoolishValueParser::new())
            .default_value("true")
        );

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

    let deduplicate_app = Command::new("deduplicate")
             .about("Collate a RAD file by corrected cell barcode")
             .version(version)
             .author(crate_authors)
             .arg(arg!(-i --"input-dir" <INPUTDIR> "input directory made by generate-permit-list that also contains the output of collate")
                 .required(true)
                 .value_parser(pathbuf_directory_exists_validator))
             .arg(arg!(-t --threads <THREADS> "number of threads to use for processing").value_parser(value_parser!(u32)).default_value(max_num_threads))
             .arg(
                arg!(-r --"rev-comp" <REVERSECOMPLEMENT> "reverse complement")
                .value_parser(clap::builder::BoolishValueParser::new())
                .default_value("true")
            );


    let opts = Command::new("piscem-atac")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .version(version)
        .author(crate_authors)
        .about("Process RAD files from the command line")
        .subcommand(gen_app)
        .subcommand(collate_app)
        .subcommand(deduplicate_app)
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

    if let Some(t) = opts.subcommand_matches("generate-permit-list") {
        let input_dir: &PathBuf = t.get_one("input").expect("no input directory specified");
        let output_dir: &PathBuf = t
            .get_one("output-dir")
            .expect("no output directory specified");
        let mut fmeth = CellFilterMethod::KneeFinding;

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
        let rc: bool = *t.get_one("rev-comp").expect("reverse comp must be boolean");

        let gpl_opts = GenPermitListOpts::builder()
            .input_dir(input_dir)
            .output_dir(output_dir)
            .fmeth(fmeth)
            .rc(rc)
            .version(VERSION)
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

    if let Some(t) = opts.subcommand_matches("collate") {
        let input_dir: &PathBuf = t.get_one("input-dir").unwrap();
        let rad_dir: &PathBuf = t.get_one("rad-dir").unwrap();
        let num_threads = *t.get_one("threads").unwrap();
        let compress_out = t.get_flag("compress");
        let max_records: u32 = *t.get_one("max-records").unwrap();
        collate(
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

    if let Some(t) = opts.subcommand_matches("deduplicate") {
        let input_dir: &PathBuf = t.get_one("input-dir").unwrap();
        let num_threads = *t.get_one("threads").unwrap();
        let rc: bool = *t.get_one("rev-comp").expect("reverse comp must be boolean");

        let dedup_opts = DeduplicateOpts::builder()
            .input_dir(input_dir)
            .num_threads(num_threads)
            .rev(rc)
            .cmdline(&cmdline)
            .version(version)
            .log(&log)
            .build();

        let parent = std::path::Path::new(&input_dir);
        let json_path = parent.join("generate_permit_list.json");
        let col_json_path = parent.join("collate.json");

        if json_path.exists() && col_json_path.exists() {
            match deduplicate(dedup_opts) {
                Ok(_) => {}
                Err(_e) => {
                    panic!("Could not dedupicate rad file");
                }
            };
        } else {
            crit!(log,
                "The provided input directory lacks a generate_permit_list.json or collate.json file; this should not happen."
            );
            bail!("The provided input directory lacks a generate_permit_list.json or collate.json file; this should not happen.");
        }
    }
    Ok(())
}
