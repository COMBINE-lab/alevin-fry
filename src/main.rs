extern crate clap;

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
    println!("I'm using the library: {:?}", libradicl::lib_name());


    // You can handle information about subcommands by requesting their matches by name
    // (as below), requesting just the name used, or both at the same time
    match opts.subcmd {
        SubCommand::Read(t) => {
            let f = File::open(t.input).unwrap();
            let mut br = BufReader::new(f);
            let h = libradicl::RADHeader::from_bytes(&mut br);
            println!("paired : {:?}, ref_count : {:?}, num_chunks : {:?}", 
                      h.is_paired, h.ref_count, h.num_chunks);
            // file-level
            let fl_tags = libradicl::TagSection::from_bytes(&mut br);
            println!("read {:?} file-level tags", fl_tags.tags.len());
            // read-level
            let rl_tags = libradicl::TagSection::from_bytes(&mut br);
            println!("read {:?} read-level tags", rl_tags.tags.len());
            // alignment-level
            let al_tags = libradicl::TagSection::from_bytes(&mut br);
            println!("read {:?} alignemnt-level tags", al_tags.tags.len());
        }
    }
}
