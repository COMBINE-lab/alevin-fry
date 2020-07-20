extern crate slog;
extern crate fasthash;
extern crate slog_term;
extern crate needletail;
extern crate rand;
extern crate clap;

use rand::Rng;
use fasthash::{sea, RandomState};
use slog::{Drain, o, info};
use needletail::bitkmer::BitNuclKmer;
use std::collections::HashMap;
use std::io::{BufReader};
use std::fs::File;


use clap::Clap;

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
    Read(Read),
}

/// A subcommand for controlling testing
#[derive(Clap)]
struct Read {
    /// Print debug info
    #[clap(short, long)]
    input: String,
    #[clap(short, long)]
    top_k: Option<u32>
}

fn gen_random_kmer(k : usize) -> String {
    use rand::Rng;
    const CHARSET: &[u8] = b"ACGT";
    let mut rng = rand::thread_rng();
    let s : String = (0..k).map(|_| {
            let idx = rng.gen_range(0, CHARSET.len());
            CHARSET[idx] as char
        })
        .collect();
    s
}

fn main() {
    let opts: Opts = Opts::parse();

    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator).build().fuse();
    let drain = slog_async::Async::new(drain).build().fuse();
    
    let log = slog::Logger::root(drain, o!());

    /*
    for _ in (0..1000) {
        let sk = gen_random_kmer(16);
        let (_, k, _) = BitNuclKmer::new(sk.as_bytes(), 16, false).next().unwrap();
        println!("{}\t{:?}", sk, k.0);
    }
    */


    info!(log, "I'm using the library: {:?}", libradicl::lib_name());
    // You can handle information about subcommands by requesting their matches by name
    // (as below), requesting just the name used, or both at the same time
    match opts.subcmd {
        SubCommand::Read(t) => {
            let f = File::open(t.input).unwrap();
            let mut br = BufReader::new(f);
            let h = libradicl::RADHeader::from_bytes(&mut br);
            info!(log, "paired : {:?}, ref_count : {:?}, num_chunks : {:?}", 
                      h.is_paired, h.ref_count, h.num_chunks);
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

            for _ in 0..(h.num_chunks as usize) {
                match (bct, umit) {
                    (3, 3) => {
                        let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U32, libradicl::RADIntID::U32);
                        libradicl::update_barcode_hist(&mut hm, &c);
                        num_reads += c.reads.len();
                        //info!(log, "{:?}", c)
                    },
                    (3, 4) => {
                        let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U32, libradicl::RADIntID::U64);
                        libradicl::update_barcode_hist(&mut hm, &c);
                        num_reads += c.reads.len();
                        //info!(log, "{:?}", c)
                    },
                    (4, 3) => {
                        let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U64, libradicl::RADIntID::U32);
                        libradicl::update_barcode_hist(&mut hm, &c);
                        num_reads += c.reads.len();
                        //info!(log, "{:?}", c)
                    },
                    (4, 4) => {
                        let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U64, libradicl::RADIntID::U64);
                        libradicl::update_barcode_hist(&mut hm, &c);
                        num_reads += c.reads.len();
                        //info!(log, "{:?}", c)
                    },
                    (_, _) => info!(log, "types not supported")
                }
            }

            info!(log, "observed {:?} reads in {:?} chunks", num_reads, h.num_chunks);

            let mut freq : Vec<u64> = hm.values().cloned().collect();
            freq.sort_unstable();
            freq.reverse();

            let num_bc = t.top_k.unwrap_or((freq.len() - 1)as u32) as usize;
            let min_freq = freq[num_bc];

            // collect all of the barcodes that have a frequency 
            // >= to min_thresh.
            let valid_bc = libradicl::permit_list_from_threshold(&hm, min_freq);
            // println!("{:?}", valid_bc);

            // generate the map from each permitted barcode to some barcode within
            // edit distance 1 of it.
            let full_permit_list = libradicl::utils::generate_permitlist_map(
                &valid_bc, 
                ft_vals.bclen as usize).unwrap();
            //println!("{:?}", full_permit_list);

            for (k,v) in hm.iter_mut() {
                //full_permit_list.get()
                match full_permit_list.get(k) {
                    Some(&valid_key) => { 
                        println!("{} was a neighbor of {}, with count {}", k, valid_key, v);
                    },
                    None => {}
                }
            }

        }
    }
}
