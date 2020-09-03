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

use clap::{crate_authors, crate_version, App, Arg};
use mimalloc::MiMalloc;
use rand::Rng;
use slog::{o, Drain};

use libradicl::bulkrad::read_bulkRAD;

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

    let view_app = App::new("view")
        .about("View the RAD file")
        .version(version)
        .author(crate_authors)
        .arg(Arg::from("-i, --input=<input>  'input RAD file'"));

    let opts = App::new("rad")
        .version(version)
        .author(crate_authors)
        .about("Process RAD files from the command line")
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
    if let Some(ref t) = opts.subcommand_matches("view") {
        let input_file: String = t.value_of_t("input").expect("no input directory specified");
       
        read_bulkRAD(input_file, 16, &log).unwrap();
    }
}
