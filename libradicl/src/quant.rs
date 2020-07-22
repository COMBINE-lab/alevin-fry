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

struct EqMapEntry {
    umis : Vec<(u64, u32)>,
    eq_num : u32
}

fn eqc_from_chunk(cell_chunk : & libradicl::Chunk, 
                  nref : u32,
                  eqc_map : &mut HashMap<Vec<u32>, EqMapEntry, fasthash::RandomState<fasthash::sea::Hash64>>) {
    
    let mut label_list_size = 0usize;
    let mut ref_to_eqid : Vec<Vec<u32>> = vec![ vec![]; nref as usize];
    let mut eq_label_starts : Vec<u32> = Vec::new();
    let mut eq_labels : Vec<u32> = Vec::new();

    // gather the equivalence class map 
    // which is a map from 
    // target set -> Vec< (UMI, count) >
    for r in &cell_chunk.reads {
        match eqc_map.get_mut(&r.refs) {
            Some(v) => { v.umis.push((r.umi, 1)); },
            None => { 
                // each reference in this equivalence class label 
                // will have to point to this equivalence class id
                let eq_num = eqc_map.len() as u32;
                label_list_size += eq_num as usize;
                for r in r.refs.iter() {
                    ref_to_eqid[*r as usize].push(eq_num);
                }
                eq_label_starts.push(eq_label_starts.len() as u32);
                eq_labels.extend(r.refs.iter());

                eqc_map.insert(r.refs.clone(), EqMapEntry { umis : vec![(r.umi,1)], eq_num }); 
            }
        }
    }
    eq_label_starts.push(eq_label_starts.len() as u32);

    // initially we inserted duplicate UMIs
    // here, collapse them and keep track of their count
    for (k, mut v) in eqc_map.iter_mut() {
        // sort so dups are adjacent
        v.umis.sort();
        // we need a copy of the vector b/c we 
        // can't easily modify it in place
        // at least I haven't seen how (@k3yavi, help here if you can).
        let cv = v.umis.clone();
        // since we have a copy, clear the original to fill it
        // with the new contents.
        v.umis.clear();
        
        let mut count = 1;
        let mut cur_elem = cv.first().unwrap().0;
        for e in cv.iter().skip(1) {
            if e.0 == cur_elem {
                count += 1;
            } else {
                v.umis.push((cur_elem, count));
                cur_elem = e.0;
                count = 1;
            }
        }
        // remember to push the last element, since we 
        // won't see a subsequent "different" element.
        v.umis.push((cur_elem, count));
    }

    println!("{:#?}", ref_to_eqid);
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
    let pbar = ProgressBar::new(hdr.num_chunks);
    pbar.set_style(ProgressStyle::default_bar()
    .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}")
    .progress_chars("╢▌▌░╟")); 

    let s = RandomState::<sea::Hash64>::new();
    let mut eq_map = HashMap::with_hasher(s);
    for _ in 0..(hdr.num_chunks as usize) {
        eq_map.clear();
        match (bct, umit) {
            (3, 3) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U32, libradicl::RADIntID::U32);
                let eq = eqc_from_chunk(&c, h.ref_count, &mut eq_map);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (3, 4) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U32, libradicl::RADIntID::U64);
                let eq = eqc_from_chunk(&c, h.ref_count, &mut eq_map);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (4, 3) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U64, libradicl::RADIntID::U32);
                let eq = eqc_from_chunk(&c, h.ref_count, &mut eq_map);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (4, 4) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U64, libradicl::RADIntID::U64);
                let eq = eqc_from_chunk(&c, h.ref_count, &mut eq_map);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (_, _) => info!(log, "types not supported")
        }

        pbar.inc(1);
    }

    pbar.finish_with_message("processed all cells.");
    Ok(())
}