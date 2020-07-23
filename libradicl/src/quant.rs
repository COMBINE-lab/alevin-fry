extern crate slog;
extern crate serde;
extern crate bincode;
extern crate indicatif;
extern crate fasthash;
extern crate petgraph;

use self::indicatif::{ProgressBar, ProgressStyle};
use self::slog::{info};
use std::io;
use std::collections::{HashMap, HashSet};
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

    // temporary map of equivalence class label to assigned 
    // index.
    let s = RandomState::<sea::Hash64>::new(); 
    let mut eqid_map : HashMap::<Vec<u32>, u32, fasthash::RandomState<fasthash::sea::Hash64>> = HashMap::with_hasher(s);

    // gather the equivalence class map 
    // which is a map from 
    // target set -> Vec< (UMI, count) >
    for r in &cell_chunk.reads {
        match eqid_map.get_mut(&r.refs) {
            Some(v) => { eqmap.eqc_info[*v as usize].umis.push((r.umi, 1)); },
            None => { 
                // each reference in this equivalence class label 
                // will have to point to this equivalence class id
                let eq_num = eqmap.eqc_info.len() as u32;
                eqmap.label_list_size += r.refs.len();
                for r in r.refs.iter() {
                    let ridx = *r as usize;
                    // if this is the first time we've seen
                    // this transcript, add it to the list of 
                    // active transcripts.
                    if eqmap.label_counts[ridx] == 0 {
                        eqmap.active_refs.push(*r);
                    }
                    eqmap.label_counts[ridx] += 1;
                    //ref_to_eqid[*r as usize].push(eq_num);
                }
                eqmap.eq_label_starts.push(eqmap.eq_labels.len() as u32);
                eqmap.eq_labels.extend(&r.refs);
                eqmap.eqc_info.push(EqMapEntry { umis : vec![(r.umi,1)], eq_num });
                eqid_map.insert(r.refs.clone(), eq_num);
                //eqmap.eqc_map.insert(r.refs.clone(), EqMapEntry { umis : vec![(r.umi,1)], eq_num }); 
            }
        }
    }
    // final value to avoid special cases
    eqmap.eq_label_starts.push(eqmap.eq_labels.len() as u32);

    eqmap.fill_ref_offsets();
    eqmap.fill_label_sizes();

    //println!("last offset = {}, num_eqc = {}\n", eqmap.eq_label_starts.last().unwrap(), eqid_map.len() );

    //let mut ref_offsets = label_counts.iter().scan(0, |sum, i| {*sum += i; Some(*sum)}).collect::<Vec<_>>();
    // final value to avoid special cases
    //ref_offsets.push(*ref_offsets.last().unwrap());
    //let mut ref_labels = vec![ u32::MAX; label_list_size];

    let mut written_labels = 0usize;

    // initially we inserted duplicate UMIs
    // here, collapse them and keep track of their count
    for idx in 0..eqmap.num_eq_classes() { //} eqmap.eqc_info.iter_mut().enumerate() {

        // for each reference in this 
        // label, put it in the next free spot
        // TODO: @k3yavi, can we avoid this copy?
        let label = eqmap.refs_for_eqc(idx as u32).to_vec();
        //println!("{:?}", label);
        for r in label {// k.iter() {
            eqmap.ref_offsets[r as usize] -= 1;
            eqmap.ref_labels[eqmap.ref_offsets[r as usize] as usize] = eqmap.eqc_info[idx].eq_num;
            written_labels += 1;
        }

        let v = &mut eqmap.eqc_info[idx];
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
    /*
    let mut cnt = 0usize;
    for e in eqmap.ref_labels.iter() {
        if *e != u32::MAX { cnt += 1;  }
    }
    
    println!("total number of labels = {}, expected = {}, written = {}", cnt, eqmap.label_list_size, written_labels);
    */
}


fn extract_graph(eqmap : &EqMap, 
                 log : &slog::Logger ) -> petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed> {
    let verbose = false;
    
    fn get_set() -> HashSet<u32, fasthash::sea::Hash64> {
        HashSet::with_hasher(sea::Hash64)
    }

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
    let mut ctr = 0;
    for eqid in 0..eqmap.num_eq_classes(){
        //let eqid = eq_val.eq_num;
        
        if verbose && ctr % 1000 == 0 {
            print!("\rprocessed {:?} eq classes", eqid);
            io::stdout().flush().ok().expect("Could not flush stdout");
        }
        ctr += 1;

        let eq = &eqmap.eqc_info[eqid];
        // for each equivalence class
        // let eq = &ce.eq_classes[eqid as usize];
        let u1 = &eq.umis;

        // for each (umi, count) pair and its index
        for (xi, x) in u1.iter().enumerate() {
            // add a node
            graph.add_node((eqid as u32, xi as u32));
            
            // for each (umi, count) pair and node after this one
            for xi2 in (xi + 1)..u1.len() {
                // x2 is the other (umi, count) pair
                let x2 = &u1[xi2];
                // add a node for it
                graph.add_node((eqid as u32, xi2 as u32));
                // determine if an edge exists between x and x2, and if so, what kind
                let et = has_edge(&x, &x2);
                // for each type of edge, add the appropriate edge in the graph
                match et {
                    PUGEdgeType::BiDirected => {
                        graph.add_edge((eqid as u32, xi as u32), (eqid as u32, xi2 as u32), ());
                        graph.add_edge((eqid as u32, xi2 as u32), (eqid as u32, xi as u32), ());
                        bidirected += 1;
                        //if multi_gene_vec[eqid] == true {
                        //    bidirected_in_multigene += 1;
                        //}
                    }
                    PUGEdgeType::XToY => {
                        graph.add_edge((eqid as u32, xi as u32), (eqid as u32, xi2 as u32), ());
                        unidirected += 1;
                        //if multi_gene_vec[eqid] == true {
                        //    unidirected_in_multigene += 1;
                        //}
                    }
                    PUGEdgeType::YToX => {
                        graph.add_edge((eqid as u32, xi2 as u32), (eqid as u32, xi as u32), ());
                        unidirected += 1;
                        //if multi_gene_vec[eqid] == true {
                        //    unidirected_in_multigene += 1;
                        //}
                    }
                    PUGEdgeType::NoEdge => {}
                }
            }
        }

        let mut hset = get_set();
        // for every reference id in this eq class 
        for r in eqmap.refs_for_eqc(eqid as u32) {
            // find the equivalence classes sharing this reference
            for eq2id in eqmap.eq_classes_containing(*r).iter() { 
                if (*eq2id as usize) <= eqid {
                    continue;
                }
                if hset.contains(eq2id) {
                    continue;
                }
                hset.insert(*eq2id);
                let eq2 = &eqmap.eqc_info[*eq2id as usize];
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
                            PUGEdgeType::BiDirected => {
                                graph.add_edge((eqid as u32, xi as u32), (*eq2id, yi as u32), ());
                                graph.add_edge((*eq2id, yi as u32), (eqid as u32, xi as u32), ());
                                bidirected += 1;
                                //if multi_gene_vec[eqid] == true
                                //    || multi_gene_vec[*eq2id as usize] == true
                                //{
                                //    bidirected_in_multigene += 1;
                                //}
                            }
                            PUGEdgeType::XToY => {
                                graph.add_edge((eqid as u32, xi as u32), (*eq2id, yi as u32), ());
                                unidirected += 1;
                                //if multi_gene_vec[eqid] == true
                                //    || multi_gene_vec[*eq2id as usize] == true
                                //{
                                //    unidirected_in_multigene += 1;
                                //}
                            }
                            PUGEdgeType::YToX => {
                                graph.add_edge((*eq2id, yi as u32), (eqid as u32, xi as u32), ());
                                unidirected += 1;
                                //if multi_gene_vec[eqid] == true
                                //    || multi_gene_vec[*eq2id as usize] == true
                                //{
                                //    unidirected_in_multigene += 1;
                                //}
                            }
                            PUGEdgeType::NoEdge => {}
                        }
                    }
                }
            }
        }
    }
    if verbose {
        //info!("tid_map of size {:?}", tid_map.len());
        info!(log, 
            "size of graph ({:?}, {:?})",
            graph.node_count(),
            graph.edge_count()
        );
    }
    graph
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
    
    let pbar = ProgressBar::new(hdr.num_chunks);
    pbar.set_style(ProgressStyle::default_bar()
    .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}")
    .progress_chars("╢▌▌░╟")); 

    let s = RandomState::<sea::Hash64>::new();
    let mut eq_map = EqMap::new(s, hdr.ref_count as u32);

    for _ in 0..(hdr.num_chunks as usize) {
        eq_map.clear();
        match (bct, umit) {
            (3, 3) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U32, libradicl::RADIntID::U32);
                eqc_from_chunk(&c, &mut eq_map);
                let g = extract_graph(&eq_map, log);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (3, 4) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U32, libradicl::RADIntID::U64);
                eqc_from_chunk(&c, &mut eq_map);
                let g = extract_graph(&eq_map, log);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (4, 3) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U64, libradicl::RADIntID::U32);
                eqc_from_chunk(&c, &mut eq_map);
                let g = extract_graph(&eq_map, log);
                num_reads += c.reads.len();
                //info!(log, "{:?}", c)
            },
            (4, 4) => {
                let c = libradicl::Chunk::from_bytes(&mut br, libradicl::RADIntID::U64, libradicl::RADIntID::U64);
                eqc_from_chunk(&c, &mut eq_map);
                let g = extract_graph(&eq_map, log);
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