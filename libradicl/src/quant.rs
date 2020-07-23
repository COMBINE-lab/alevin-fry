extern crate slog;
extern crate serde;
extern crate bincode;
extern crate indicatif;
extern crate fasthash;
extern crate petgraph;

use self::indicatif::{ProgressBar, ProgressStyle};
use self::slog::{info};
use std::collections::HashMap;
use fasthash::{sea, RandomState};
use std::io::{BufReader};
use std::fs::File;
use std::path;
use std::io::{BufWriter, Write, Read, Seek, SeekFrom};
use self::petgraph::prelude::*;

use crate as libradicl;

use self::libradicl::schema::{EqMap, EqMapEntry, PUGEdgeType};
use self::libradicl::utils::*;

fn eqc_from_chunk(cell_chunk : & libradicl::Chunk, 
                  eqmap : &mut EqMap)
                  /*
                  nref : u32,
                  EQMap
                  eqc_map : &mut HashMap<Vec<u32>, EqMapEntry, fasthash::RandomState<fasthash::sea::Hash64>>)*/ {
    
    //let mut label_list_size = 0usize;
    
    // concatenated lists of the labels of all equivalence classes
    //let mut eq_labels : Vec<u32> = Vec::new();
    // vector that deliniates where each equivalence class label 
    // begins and ends.  The label for equivalence class i begins 
    // at offset eq_label_starts[i], and it ends at 
    // eq_label_starts[i+1].  The length of this vector is 1 greater
    // than the number of equivalence classes.
    //let mut eq_label_starts : Vec<u32> = Vec::new();

    //let mut label_counts : Vec<u32> = vec![0; nref as usize];

    // gather the equivalence class map 
    // which is a map from 
    // target set -> Vec< (UMI, count) >
    for r in &cell_chunk.reads {
        match eqmap.eqc_map.get_mut(&r.refs) {
            Some(v) => { v.umis.push((r.umi, 1)); },
            None => { 
                // each reference in this equivalence class label 
                // will have to point to this equivalence class id
                let eq_num = eqmap.eqc_map.len() as u32;
                eqmap.label_list_size += eq_num as usize;
                for r in r.refs.iter() {
                    eqmap.label_counts[*r as usize] += 1;
                    //ref_to_eqid[*r as usize].push(eq_num);
                }
                eqmap.eq_label_starts.push(eqmap.eq_label_starts.len() as u32);
                eqmap.eq_labels.extend(r.refs.iter());

                eqmap.eqc_map.insert(r.refs.clone(), EqMapEntry { umis : vec![(r.umi,1)], eq_num }); 
            }
        }
    }
    // final value to avoid special cases
    eqmap.eq_label_starts.push(eqmap.eq_label_starts.len() as u32);

    eqmap.fill_ref_offsets();
    eqmap.fill_label_sizes();

    //let mut ref_offsets = label_counts.iter().scan(0, |sum, i| {*sum += i; Some(*sum)}).collect::<Vec<_>>();
    // final value to avoid special cases
    //ref_offsets.push(*ref_offsets.last().unwrap());
    //let mut ref_labels = vec![ u32::MAX; label_list_size];

    // initially we inserted duplicate UMIs
    // here, collapse them and keep track of their count
    for (k, v) in eqmap.eqc_map.iter_mut() {

        // for each reference in this 
        // label, put it in the next free spot
        for r in k.iter() {
            eqmap.ref_offsets[*r as usize] -= 1;
            eqmap.ref_labels[eqmap.ref_offsets[*r as usize] as usize] = v.eq_num;
        }

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

    println!("{:#?}", &eqmap.eq_labels[1..100]);
    //println!("{:#?}", ref_offsets);

}

/*
fn extract_graph() {

    // given 2 pairs (UMI, count), determine if an edge exists 
    // between them, and if so, what type.
    let has_edge = |x: &(u64, u32), y: &(u64, u32)| -> PUGEdgeType {
        if &x.0 == &y.0 {
            return PUGEdgeType::BiDirected;
        }
        if x.1 > (2 * y.1 - 1) {
            if count_diff_2_bit_packed(x.0, y.0) < 2 {
                return PUGEdgeType::XToY;
            } else {
                return PUGEdgeType::NoEdge;
            }
        } else if y.1 > (2 * x.1 - 1) {
            if count_diff_2_bit_packed(x.0, y.0) < 2 {
                return PUGEdgeType::YToX;
            } else {
                return PUGEdgeType::NoEdge;
            }
        }
        PUGEdgeType::NoEdge
    };

    let mut bidirected = 0u64;
    let mut unidirected = 0u64;

    let mut graph = DiGraphMap::<(u32, u32), ()>::new();
    for eqid in 0..ce.eq_classes.len() {
        if verbose && eqid % 1000 == 0 {
            print!("\rprocessed {:?} eq classes", eqid);
            io::stdout().flush().ok().expect("Could not flush stdout");
        }
        // for each equivalence class
        let eq = &ce.eq_classes[eqid as usize];
        let u1 = &eq.umis;
        for (xi, x) in u1.iter().enumerate() {
            graph.add_node((eqid as u32, xi as u32));
            for xi2 in (xi + 1)..u1.len() {
                let x2 = &u1[xi2];
                graph.add_node((eqid as u32, xi2 as u32));
                let et = has_edge(&x, &x2);
                match et {
                    EdgeType::BiDirected => {
                        graph.add_edge((eqid as u32, xi as u32), (eqid as u32, xi2 as u32), ());
                        graph.add_edge((eqid as u32, xi2 as u32), (eqid as u32, xi as u32), ());
                        bidirected += 1;
                        if multi_gene_vec[eqid] == true {
                            bidirected_in_multigene += 1;
                        }
                    }
                    EdgeType::XToY => {
                        graph.add_edge((eqid as u32, xi as u32), (eqid as u32, xi2 as u32), ());
                        unidirected += 1;
                        if multi_gene_vec[eqid] == true {
                            unidirected_in_multigene += 1;
                        }
                    }
                    EdgeType::YToX => {
                        graph.add_edge((eqid as u32, xi2 as u32), (eqid as u32, xi as u32), ());
                        unidirected += 1;
                        if multi_gene_vec[eqid] == true {
                            unidirected_in_multigene += 1;
                        }
                    }
                    EdgeType::NoEdge => {}
                }
            }
        }

        let mut hset = get_set();
        // for every transcript
        for t in eq.transcripts.iter() {
            // find the equivalence classes sharing this transcript
            for eq2id in tid_map[t].iter() {
                if (*eq2id as usize) <= eqid {
                    continue;
                }
                if hset.contains(eq2id) {
                    continue;
                }
                hset.insert(*eq2id);
                let eq2 = &ce.eq_classes[*eq2id as usize];
                // compare all the umis
                let u2 = &eq2.umis;
                for (xi, x) in u1.iter().enumerate() {
                    // Node for equiv : eqid and umi : xi
                    graph.add_node((eqid as u32, xi as u32));
                    for (yi, y) in u2.iter().enumerate() {
                        // Node for equiv : eq2id and umi : yi
                        graph.add_node((*eq2id as u32, yi as u32));
                        let et = has_edge(&x, &y);
                        match et {
                            EdgeType::BiDirected => {
                                graph.add_edge((eqid as u32, xi as u32), (*eq2id, yi as u32), ());
                                graph.add_edge((*eq2id, yi as u32), (eqid as u32, xi as u32), ());
                                bidirected += 1;
                                if multi_gene_vec[eqid] == true
                                    || multi_gene_vec[*eq2id as usize] == true
                                {
                                    bidirected_in_multigene += 1;
                                }
                            }
                            EdgeType::XToY => {
                                graph.add_edge((eqid as u32, xi as u32), (*eq2id, yi as u32), ());
                                unidirected += 1;
                                if multi_gene_vec[eqid] == true
                                    || multi_gene_vec[*eq2id as usize] == true
                                {
                                    unidirected_in_multigene += 1;
                                }
                            }
                            EdgeType::YToX => {
                                graph.add_edge((*eq2id, yi as u32), (eqid as u32, xi as u32), ());
                                unidirected += 1;
                                if multi_gene_vec[eqid] == true
                                    || multi_gene_vec[*eq2id as usize] == true
                                {
                                    unidirected_in_multigene += 1;
                                }
                            }
                            EdgeType::NoEdge => {}
                        }
                    }
                }
            }
        }
    }
    if verbose {
        info!("tid_map of size {:?}", tid_map.len());
        info!(
            "size of graph ({:?}, {:?})",
            graph.node_count(),
            graph.edge_count()
        );
    }

}
*/

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
    //let mut eq_map = HashMap::with_hasher(s);

    let mut eq_map = EqMap::new(s, hdr.ref_count as u32);

    for _ in 0..(hdr.num_chunks as usize) {
        eq_map.clear();
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

        pbar.inc(1);
    }

    pbar.finish_with_message("processed all cells.");
    Ok(())
}