extern crate slog;
extern crate serde;
extern crate bincode;
extern crate indicatif;
extern crate fasthash;

use self::indicatif::{ProgressBar, ProgressStyle};
use self::slog::{info};
use std::collections::HashMap;
use fasthash::{sea, RandomState};
use std::io::{BufReader};
use std::fs::File;
use std::path;
use std::io::{BufWriter, Write, Read, Seek, SeekFrom};

use crate as libradicl;



fn eqc_from_chunk(cell_chunk : & libradicl::Chunk, 
                  eqc_map : &mut HashMap<Vec<u32>, Vec<(u64, u32)>, fasthash::RandomState<fasthash::sea::Hash64>>) {
    
    // equivalence classes 
    for r in &cell_chunk.reads {
        match eqc_map.get_mut(&r.refs) {
            Some(v) => { v.push((r.umi, 1)); },
            None => { eqc_map.insert(r.refs.clone(), vec![(r.umi,1)]); }
        }
        //eqc_map.entry(&r.refs).or_insert(vec![]).push((r.umi,1));
    }
    println!("for cell 1 eq_map has size {:?}", eqc_map.len());

}


pub fn quantify(input_dir : String, log : &slog::Logger ) -> Result<(), Box<dyn std::error::Error>> {
    let parent = std::path::Path::new(&input_dir);
    let i_file = File::open(parent.join("map.collated.rad")).unwrap();
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
    
    //Ok(())
    
    let s = RandomState::<sea::Hash64>::new();
    let mut eq_map = HashMap::with_hasher(s);
    for _ in 0..(hdr.num_chunks as usize) {
        match (bct, umit) {
            (3, 3) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U32, libradicl::RADIntID::U32);
                let eq = eqc_from_chunk(&c, &mut eq_map);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (3, 4) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U32, libradicl::RADIntID::U64);
                let eq = eqc_from_chunk(&c, &mut eq_map);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (4, 3) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U64, libradicl::RADIntID::U32);
                let eq = eqc_from_chunk(&c, &mut eq_map);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (4, 4) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U64, libradicl::RADIntID::U64);
                let eq = eqc_from_chunk(&c, &mut eq_map);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (_, _) => info!(log, "types not supported")
        }
    }
    Ok(())
}