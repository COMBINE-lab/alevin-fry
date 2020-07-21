extern crate slog;
extern crate fasthash;
extern crate slog_term;
extern crate needletail;
extern crate rand;
extern crate clap;
extern crate serde;
extern crate bincode;

use bincode::{deserialize, serialize};
use serde::{Serialize, Deserialize};
use rand::Rng;
use fasthash::{sea, RandomState};
use slog::{Drain, o, info};
use needletail::bitkmer::BitNuclKmer;
use std::collections::HashMap;
use std::io::{BufReader};
use std::fs::File;
use std::io::{BufWriter, Write, Read, Seek, SeekFrom};

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
    GeneratePermitList(GeneratePermitList),
    #[clap(version = "0.0.1", author = "Avi Srivastava, Rob Patro")]
    Collate(Collate)
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
    top_k: Option<u32>
}

#[derive(Clap)]
struct Collate {
   /// 
   #[clap(short, long)]
   input_dir : String,
   #[clap(short, long)]
   rad_file : String,
   #[clap(short, long, default_value = "10000000")]
   max_records : u32
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


fn collate(t : Collate, log : &slog::Logger) -> Result<(), Box<dyn std::error::Error>> {

    let parent = std::path::Path::new(&t.input_dir);
    let mut ofile = File::create(parent.join("map.collated.rad")).unwrap();
    let mut f = File::open(&t.rad_file).unwrap();
    let f2 = f.try_clone().unwrap();

    let mut br = BufReader::new(f);
        
    let h = libradicl::RADHeader::from_bytes(&mut br);

    let end_header_pos = br.get_ref().seek(SeekFrom::Current(0)).unwrap() - (br.buffer().len() as u64);

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
    
    let pos = br.get_ref().seek(SeekFrom::Current(0)).unwrap() - (br.buffer().len() as u64);

    // copy the header 
    {
        br.get_mut().seek(SeekFrom::Start(0));
        let mut br2 = BufReader::new(br.get_ref());
        std::io::copy(&mut br2.by_ref().take(pos), &mut ofile);
    }
    
    let mut owriter = BufWriter::new(ofile);

    type TSVRec = (u64, u64);
    let mut tsv_map = Vec::<(u64,u64)>::new();//HashMap::<u64, u64>::new();

    let freq_file = std::fs::File::open(parent.join("permit_freq.tsv")).expect("couldn't open file");
    let mut rdr = csv::ReaderBuilder::new()
                  .has_headers(false)
                  .delimiter(b'\t').from_reader(freq_file);
    
    let mut total_to_collate = 0;
    for result in rdr.deserialize() {
        let record : TSVRec = result?;
        tsv_map.push(record);
        total_to_collate += record.1;
    }

    info!(log, "size of TSVMap = {:?}", tsv_map.len());
    info!(log, "tsv_map = {:?}", tsv_map);

    // get the correction map
    let cmfile = std::fs::File::open(parent.join("permit_map.bin")).unwrap();
    let correct_map : HashMap<u64,u64> = bincode::deserialize_from(&cmfile).unwrap();

    info!(log, "deserialized correction map of length : {:?}", correct_map.len());

    let mut pass_num = 1;
    let cc = libradicl::ChunkConfig{
                num_chunks : h.num_chunks, 
                bc_type : bct,
                umi_type : umit
            };

    let mut output_cache = HashMap::<u64, libradicl::CorrectedCBChunk>::new();
    let mut allocated_records = 0;
    let mut total_allocated_records = 0;
    let mut last_idx = 0;
    let mut num_output_chunks = 0;

    while last_idx < tsv_map.len() {
        allocated_records = 0;
        output_cache.clear();

        // The tsv_map tells us, for each "true" barcode
        // how many records belong to it.  We can scan this information
        // to determine what true barcodes we will keep in memory. 
        let init_offset = last_idx;
        for (i, rec) in tsv_map[init_offset..].iter().enumerate() {
            output_cache.insert(rec.0, libradicl::CorrectedCBChunk::from_counter(rec.1));
            allocated_records += rec.1;
            last_idx = i + 1;
            if allocated_records >= (t.max_records as u64){
                info!(log, "pass {:?} will collate {:?} records", pass_num, allocated_records);
                break;
            }
        }

        num_output_chunks += output_cache.len();
        total_allocated_records += allocated_records;
        last_idx += init_offset;

        // collect the output for the current barcode set
        libradicl::collect_records(&mut br, &cc, &correct_map, &mut output_cache);

        // dump the output we have
        libradicl::dump_output_cache(&mut owriter, &output_cache);

        // reset the reader to start of the chunks
        br.get_ref().seek(SeekFrom::Start(pos)).unwrap();
        pass_num += 1;
        info!(log, "total collated {:?} / {:?}", total_allocated_records, total_to_collate);
        info!(log, "last index processed {:?} / {:?}", last_idx, tsv_map.len());
    }

    info!(log, "writing num output chunks ({:?}) to header", num_output_chunks);
    owriter.get_ref().seek(SeekFrom::Start( end_header_pos - (std::mem::size_of::<u64>() as u64) ));
    owriter.write(&num_output_chunks.to_le_bytes());

    Ok(())
}

fn generate_permit_list(t : GeneratePermitList, log : &slog::Logger) -> Result<u64, Box<dyn std::error::Error>> {
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

    // generate the map from each permitted barcode to all barcodes within
    // edit distance 1 of it.
    let full_permit_list = libradicl::utils::generate_permitlist_map(
        &valid_bc, 
        ft_vals.bclen as usize).unwrap();

    let s2 = RandomState::<sea::Hash64>::new();
    let mut permitted_map = HashMap::with_capacity_and_hasher(valid_bc.len(), s2); 

    let mut num_corrected = 0;
    for (k,v) in hm.iter() {
        match full_permit_list.get(k) {
            Some(&valid_key) => { 
                *permitted_map.entry(valid_key).or_insert(0u64) += *v;
                num_corrected += 1;
                //println!("{} was a neighbor of {}, with count {}", k, valid_key, v);
            },
            None => {}
        }
    }

    let parent = std::path::Path::new(&t.output_dir);
    std::fs::create_dir_all(&parent).unwrap();
    let o_path = parent.join("permit_freq.tsv"); 
    let output = std::fs::File::create(&o_path).expect("could not create output.");
    let mut writer = BufWriter::new(&output);

    for (k,v) in permitted_map {
        writeln!(&mut writer, "{:?}\t{:?}", k, v);
    }

    let s_path = parent.join("permit_map.bin"); 
    let s_file = std::fs::File::create(&s_path).expect("could not create serialization file.");
    let mut s_writer = BufWriter::new(&s_file);
    bincode::serialize_into(&mut s_writer, &full_permit_list);

    Ok(num_corrected)
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
        SubCommand::GeneratePermitList(t) => {
            let nc = generate_permit_list(t, &log).unwrap();
            info!(log, "total number of corrected barcodes : {}", nc);
        },
        SubCommand::Collate(t) => {
            collate(t, &log).unwrap();
        }
    }
}
