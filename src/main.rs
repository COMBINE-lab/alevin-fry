extern crate bincode;
extern crate bio_types;
extern crate chrono;
extern crate clap;
extern crate fasthash;
extern crate num_cpus;
extern crate rand;
extern crate serde;
extern crate slog;
extern crate slog_term;

use clap::{Arg, App};
use bio_types::strand::Strand;
use fasthash::{sea, RandomState};
use mimalloc::MiMalloc;
use rand::Rng;
use slog::crit;
use slog::{info, o, Drain};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::{BufWriter, Write};
use libradicl::schema::ResolutionStrategy;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

/// This doc string acts as a help message when the user runs '--help'
/// as do all doc strings on fields
/*
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
    #[clap(short, long)]
    valid_bc: Option<String>,
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
    resolution: Option<u32>
}
*/

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
    input_file: String,
    output_dir: String,
    top_k: Option<usize>,
    valid_bc_file: Option<String>,
    log: &slog::Logger,
) -> Result<u64, Box<dyn std::error::Error>> {
    let i_file = File::open(input_file).unwrap();
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

    let valid_bc: Vec<u64>;

    if top_k.is_some() {
        let top_k = top_k.unwrap() as usize;
        let mut freq: Vec<u64> = hm.values().cloned().collect();
        freq.sort_unstable();
        freq.reverse();

        let num_bc = if freq.len() < top_k {
            freq.len() - 1
        } else {
            top_k - 1
        };

        let min_freq = freq[num_bc];

        // collect all of the barcodes that have a frequency
        // >= to min_thresh.
        valid_bc = libradicl::permit_list_from_threshold(&hm, min_freq);
    } else {
        let valid_bc_file = valid_bc_file.expect("couldn't extract --valid-bc option.");
        valid_bc = libradicl::permit_list_from_file(valid_bc_file, ft_vals.bclen);
    }

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

    let parent = std::path::Path::new(&output_dir);
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
    let max_num_threads : String = (num_cpus::get() as u32).to_string();

    let gen_app = App::new("generate-permit-list")
    .about("Generate a permit list of barcodes from a RAD file")
    .version("0.0.1")
    .author("Avi Srivastava, Rob Patro")
    .arg(Arg::from("-i, --input=<input>  'input RAD file'"))
    .arg(Arg::from("-o, --output-dir=<output-dir>  'output directory'"))
    .arg(Arg::from("-k, --top-k=<top-k>  'select the top-k most-frequent barcodes as valid (true)'"))
    .arg(Arg::from("-b, --valid-bc=<valid-bc> 'uses true barcode collected from a provided file'")
         .conflicts_with("top-k"));

    let collate_app = App::new("collate")
    .about("Collate a RAD file by corrected cell barcode")
    .version("0.0.1")
    .author("Avi Srivastava, Rob Patro")
    .arg(Arg::from("-i, --input-dir=<input-dir> 'input directory made by generate-permit-list'"))
    .arg(Arg::from("-r, --rad-file=<rad-file> 'the RAD file to be collated'"))
    .arg(Arg::from("-m, --max-records=[max-records] 'the maximum number of read records to keep in memory at once'")
         .default_value("10000000"))
    .arg(Arg::from("-e, --expected-ori=[expected-ori] 'the expected orientation of alignments'")
         .default_value("fw"));

    let quant_app = App::new("quant")
    .about("Quantify expression from a collated RAD file")
    .version("0.0.1")
    .author("Avi Srivastava, Rob Patro")
    .arg(Arg::from("-i, --input-dir=<input-dir>  'input directory containing collated RAD file'"))
    .arg(Arg::from("-m, --tg-map=<tg-map>  'transcript to gene map'"))
    .arg(Arg::from("-o, --output-dir=<output-dir> 'output directory where quantification results will be written'"))
    .arg(Arg::from("-t, --threads 'number of threads to use for processing'").default_value(&max_num_threads))
    .arg(Arg::from("-r, --resolution 'the resolution strategy by which molecules will be counted'")
        .possible_values(&["full", "trivial", "parsimony"])
        .default_value("full")
        .case_insensitive(true)
        .about("the resolution strategy by which molecules will be counted"));

    let opts = App::new("Radicl.")
                    .version("0.0.1")
                    .author("Avi Srivastava, Rob Patro")
                    .about("Process RAD files from the command line")
                    .subcommand(gen_app)
                    .subcommand(collate_app)
                    .subcommand(quant_app).get_matches();


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
            let input_file : String = t.value_of_t("input").expect("no input directory specified");
            let output_dir : String = t.value_of_t("output-dir").expect("no input directory specified");
            let top_k = match t.value_of_t("top-k") {
                Ok(v) => Some(v),
                Err(_) => None
            };
            let valid_bc = match t.value_of_t("valid-bc") {
                Ok(v) => Some(v),
                Err(_) => None
            };
            let nc = generate_permit_list(input_file, output_dir, top_k, valid_bc, &log).unwrap();
            info!(log, "total number of corrected barcodes : {}", nc);
    }

    if let Some(ref t) = opts.subcommand_matches("collate") {
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

        let input_dir : String = t.value_of_t("input-dir").unwrap();
        let rad_file : String = t.value_of_t("rad-file").unwrap();
        let max_records : u32 = t.value_of_t("max-records").unwrap();
        libradicl::collate::collate(input_dir, rad_file, max_records, expected_ori, &log)
            .expect("could not collate.");  
    }

    if let Some(ref t) = opts.subcommand_matches("quant") {
        let num_threads = t.value_of_t("threads").unwrap();
        let input_dir = t.value_of_t("input-dir").unwrap();
        let output_dir = t.value_of_t("output-dir").unwrap();
        let tg_map = t.value_of_t("tg-map").unwrap();
        let resolution : ResolutionStrategy = t.value_of_t("resolution").unwrap();
    libradicl::quant::quantify(
        input_dir,
        tg_map,
        output_dir,
        num_threads,
        resolution,
        &log).expect("could not quantify rad file.");
    }
}
