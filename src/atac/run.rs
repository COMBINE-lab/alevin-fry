use crate::atac::cellfilter::generate_permit_list;
use crate::atac::cellfilter::CellFilterMethod;
use crate::atac::collate::collate;
use crate::atac::deduplicate::deduplicate;
use crate::atac::prog_opts::{DeduplicateOpts, GenPermitListOpts};
use crate::atac::sort::sort;
use anyhow::bail;
use clap::ArgMatches;
use slog::{crit, info, o, warn, Logger};
use std::path::PathBuf;

pub fn run(opts: &ArgMatches, version: &str, cmdline: &str, log: &Logger) -> anyhow::Result<()> {
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
            .version(version)
            .cmdline(cmdline)
            .log(log)
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
            cmdline,
            version,
            log,
        )
        .expect("could not collate.");
    }

    if let Some(t) = opts.subcommand_matches("sort") {
        let input_dir: &PathBuf = t.get_one("input-dir").unwrap();
        let rad_dir: &PathBuf = t.get_one("rad-dir").unwrap();
        let num_threads: u32 = *t.get_one("threads").unwrap();
        let compress_out = t.get_flag("compress");
        let max_records: u32 = *t.get_one("max-records").unwrap();

        sort(
            input_dir,
            rad_dir,
            num_threads,
            max_records,
            compress_out,
            cmdline,
            version,
            log,
        )
        .expect("could not sort.");
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
