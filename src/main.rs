extern crate bincode;
extern crate bio_types;
extern crate chrono;
extern crate clap;
extern crate fasthash;
extern crate needletail;
extern crate num_cpus;
extern crate rand;
extern crate serde;
extern crate slog;
extern crate slog_term;

use bio_types::strand::Strand;
use fasthash::{sea, RandomState};
use mimalloc::MiMalloc;
use rand::Rng;
use slog::crit;
use slog::{info, o, Drain};
//use needletail::bitkmer::BitNuclKmer;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::{BufWriter, Write};

use clap::Clap;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

/// This doc string acts as a help message when the user runs '--help'
/// as do all doc strings on fields
#[derive(Clap)]
#[clap(version = "0.0.1", author = "Avi Srivastava, Rob Patro")]
struct Opts {
    #[clap(subcommand)]
    subcmd: SubCommand,
}

#[derive(Clap)]
enum SubCommand {
    #[clap(version = "0.0.1", author = "Avi Srivastava, Rob Patro")]
    GeneratePermitList(GeneratePermitList),
    #[clap(version = "0.0.1", author = "Avi Srivastava, Rob Patro")]
    Collate(Collate),
    #[clap(version = "0.0.1", author = "Avi Srivastava, Rob Patro")]
    Quant(Quant),
}

/// A subcommand for controlling testing
#[derive(Clap)]
struct GeneratePermitList {
    /// Print debug info
    #[clap(short, long)]
    input: String,
    #[clap(short, long)]
    output_dir: String,
    #[clap(short, long)]
    top_k: Option<u32>,
}

#[derive(Clap)]
struct Collate {
    ///
    #[clap(short, long)]
    input_dir: String,
    #[clap(short, long)]
    rad_file: String,
    #[clap(short, long, default_value = "10000000")]
    max_records: u32,
    #[clap(short, long, default_value = "fw")]
    expected_ori: String,
}

#[derive(Clap)]
struct Quant {
    ///
    #[clap(short, long)]
    input_dir: String,
    #[clap(short, long)]
    tg_map: String,
    #[clap(short, long)]
    output_dir: String,
    #[clap(short, long)]
    num_threads: Option<u32>,
    #[clap(short, long)]
    no_em: bool,
    #[clap(short, long)]
    naive: bool
}

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

fn generate_permit_list(
    t: GeneratePermitList,
    log: &slog::Logger,
) -> Result<u64, Box<dyn std::error::Error>> {
    let i_file = File::open(t.input).unwrap();
    let mut br = BufReader::new(i_file);
    let hdr = libradicl::RADHeader::from_bytes(&mut br);
    info!(
        log,
        "paired : {:?}, ref_count : {:?}, num_chunks : {:?}",
        hdr.is_paired,
        hdr.ref_count,
        hdr.num_chunks
    );
    // file-level
    let fl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    let ft_vals = libradicl::FileTags::from_bytes(&mut br);
    info!(log, "File-level tag values {:?}", ft_vals);

    let bct = rl_tags.tags[0].typeid;
    let umit = rl_tags.tags[1].typeid;

    let mut num_reads: usize = 0;

    let s = RandomState::<sea::Hash64>::new();
    let mut hm = HashMap::with_hasher(s);
    let bc_type = libradicl::decode_int_type_tag(bct).expect("unknown barcode type id.");
    let umi_type = libradicl::decode_int_type_tag(umit).expect("unknown barcode type id.");

    for _ in 0..(hdr.num_chunks as usize) {
        let c = libradicl::Chunk::from_bytes(&mut br, &bc_type, &umi_type);
        libradicl::update_barcode_hist(&mut hm, &c);
        num_reads += c.reads.len();
    }

    info!(
        log,
        "observed {:?} reads in {:?} chunks", num_reads, hdr.num_chunks
    );

    let mut freq: Vec<u64> = hm.values().cloned().collect();
    freq.sort_unstable();
    freq.reverse();

    let num_bc = t.top_k.unwrap_or((freq.len() - 1) as u32) as usize;
    let min_freq = freq[num_bc];

    // collect all of the barcodes that have a frequency
    // >= to min_thresh.
    let valid_bc = libradicl::permit_list_from_threshold(&hm, min_freq);

    // generate the map from each permitted barcode to all barcodes within
    // edit distance 1 of it.
    let full_permit_list =
        libradicl::utils::generate_permitlist_map(&valid_bc, ft_vals.bclen as usize).unwrap();

    let s2 = RandomState::<sea::Hash64>::new();
    let mut permitted_map = HashMap::with_capacity_and_hasher(valid_bc.len(), s2);

    let mut num_corrected = 0;
    for (k, v) in hm.iter() {
        if let Some(&valid_key) = full_permit_list.get(k) {
            *permitted_map.entry(valid_key).or_insert(0u64) += *v;
            num_corrected += 1;
            //println!("{} was a neighbor of {}, with count {}", k, valid_key, v);
        }
    }

    let parent = std::path::Path::new(&t.output_dir);
    std::fs::create_dir_all(&parent).unwrap();
    let o_path = parent.join("permit_freq.tsv");
    let output = std::fs::File::create(&o_path).expect("could not create output.");
    let mut writer = BufWriter::new(&output);

    for (k, v) in permitted_map {
        writeln!(&mut writer, "{:?}\t{:?}", k, v).expect("couldn't write to output file.");
    }

    let s_path = parent.join("permit_map.bin");
    let s_file = std::fs::File::create(&s_path).expect("could not create serialization file.");
    let mut s_writer = BufWriter::new(&s_file);
    bincode::serialize_into(&mut s_writer, &full_permit_list)
        .expect("couldn't serialize permit list.");

    Ok(num_corrected)
}

fn main() {
    let opts: Opts = Opts::parse();

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

    info!(log, "I'm using the library: {:?}", libradicl::lib_name());
    // You can handle information about subcommands by requesting their matches by name
    // (as below), requesting just the name used, or both at the same time
    match opts.subcmd {
        SubCommand::GeneratePermitList(t) => {
            let nc = generate_permit_list(t, &log).unwrap();
            info!(log, "total number of corrected barcodes : {}", nc);
        }
        SubCommand::Collate(t) => {
            let valid_ori: bool;
            let expected_ori = match t.expected_ori.to_uppercase().as_str() {
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
                    t.expected_ori
                );
                std::process::exit(1);
            }

            libradicl::collate::collate(t.input_dir, t.rad_file, t.max_records, expected_ori, &log)
                .expect("could not collate.");
        }
        SubCommand::Quant(t) => {
            let num_threads = match t.num_threads {
                Some(nt) => nt,
                None => num_cpus::get() as u32,
            };
            libradicl::quant::quantify(
                t.input_dir,
                t.tg_map,
                t.output_dir,
                num_threads,
                t.no_em,
                t.naive,
                &log,
            )
            .expect("could not quantify rad file.");
        }
    }
}
