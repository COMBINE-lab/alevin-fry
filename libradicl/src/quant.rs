extern crate slog;
extern crate serde;
extern crate bincode;
extern crate indicatif;

use self::indicatif::{ProgressBar, ProgressStyle};
use self::slog::{info};
use std::collections::HashMap;
use std::io::{BufReader};
use std::fs::File;
use std::io::{BufWriter, Write, Read, Seek, SeekFrom};

use crate as libradicl;


/*
fn eqc_from_chunk(cell_chunk : &libradicl::Chunk, 
                  eqc_map : &mut HashMap<u64, u64, fasthash::RandomState<fasthash::sea::Hash64>>) {
    
    // equivalence classes 
    for &r in cell_chunk.reads.iter() {
        eqc_map.entry(r.refs.as_slice()).or_insert(vec![r.umi]);
        
    }

    // transcript -> eq class map


}
*/

pub fn quantify(input_dir : String, log : &slog::Logger ) -> Result<(), Box<dyn std::error::Error>> {

    let i_file = File::open(input_dir).unwrap();
    let mut br = BufReader::new(i_file);
    let hdr = libradicl::RADHeader::from_bytes(&mut br);
    info!(log, "paired : {:?}, ref_count : {:?}, num_chunks : {:?}", 
              hdr.is_paired, hdr.ref_count, hdr.num_chunks);
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
    
    Ok(())
    /*
    let s = RandomState::<sea::Hash64>::new();
    let mut eq_map = HashMap::with_hasher(s);
    for _ in 0..(hdr.num_chunks as usize) {
        match (bct, umit) {
            (3, 3) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U32, libradicl::RADIntID::U32);
                let eq = eqc_from_chunk(&mut eq_map, &c);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (3, 4) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U32, libradicl::RADIntID::U64);
                let eq = eqc_from_chunk(&mut eq_map, &c);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (4, 3) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U64, libradicl::RADIntID::U32);
                let eq = eqc_from_chunk(&mut eq_map, &c);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (4, 4) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U64, libradicl::RADIntID::U64);
                let eq = eqc_from_chunk(&mut eq_map, &c);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (_, _) => info!(log, "types not supported")
        }
    }
    */


}