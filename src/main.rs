extern crate slog;
extern crate fasthash;
extern crate slog_term;

extern crate clap;
use fasthash::{sea, RandomState};
use slog::{Drain, o, info};
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
    #[clap(short)]
    input: String 
}


fn main() {
    let opts: Opts = Opts::parse();

    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator).build().fuse();
    let drain = slog_async::Async::new(drain).build().fuse();
    
    let log = slog::Logger::root(drain, o!());

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

            println!("{:?}", &freq[0..100] );
        }
    }
}
