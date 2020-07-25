
extern crate slog;
extern crate fasthash;
extern crate petgraph;

use std::collections::{HashMap, HashSet, VecDeque};
use std::iter::FromIterator;
use self::slog::{crit};
//use hashers::fnv::FNV1aHasher64;quant::
use fasthash::RandomState;
use fasthash::sea::Hash64;

use pugutils::petgraph::prelude::*;
use pugutils::petgraph::unionfind::*;
use pugutils::petgraph::visit::NodeIndexable;

use schema::{EqMap};

type GeneEqMapT = HashMap<u32, Vec<u32>, fasthash::RandomState<Hash64>>;

pub fn weakly_connected_components<G>(
    g: G,
) -> GeneEqMapT
where
    G: petgraph::visit::NodeCompactIndexable + petgraph::visit::IntoEdgeReferences,
{
    let mut vertex_sets = UnionFind::new(g.node_bound());
    for edge in g.edge_references() {
        let (a, b) = (edge.source(), edge.target());

        // union the two vertices of the edge
        vertex_sets.union(g.to_index(a), g.to_index(b));
    }
    let labels = vertex_sets.into_labeling();
    fn get_map() -> GeneEqMapT {
        let s = RandomState::<Hash64>::new();
        HashMap::with_hasher(s)
    }

    let mut components = get_map();
    for (i, v) in labels.iter().enumerate() {
        let ve = components.entry(*v as u32).or_insert_with(Vec::new);
        ve.push(i as u32);
    }
    components
}

fn collapse_vertices(
    v: u32,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap 
) -> (Vec<u32>, u32) {
    type VertexSet = HashSet<u32, fasthash::RandomState<Hash64>>;
    fn get_set(cap: u32) -> VertexSet {
        let s = RandomState::<Hash64>::new();
        VertexSet::with_capacity_and_hasher(cap as usize, s)
    }

    let mut largest_mcc: Vec<u32> = Vec::new();
    let mut chosen_txp = 0u32;
    let vert = g.from_index(v as usize);

    //unsafe {
        for txp in eqmap.refs_for_eqc(vert.0).iter() {
            let mut bfs_list = VecDeque::new();
            bfs_list.push_back(v);

            let mut visited_set = get_set(16);
            visited_set.insert(v);

            let mut current_mcc = Vec::new();

            while let Some(cv) = bfs_list.pop_front() {
                current_mcc.push(cv);

                for nv in g.neighbors_directed(g.from_index(cv as usize), Outgoing) {
                    let n = g.to_index(nv) as u32;
                    if visited_set.contains(&n) {
                        continue;
                    } else {
                        visited_set.insert(n);
                    }
                    let n_labels = eqmap.refs_for_eqc(nv.0);
                    if n_labels.contains(&txp) {
                        bfs_list.push_back(n);
                    }
                }
            }

            if largest_mcc.len() < current_mcc.len() {
                largest_mcc = current_mcc;
                chosen_txp = *txp;
            }
        }
    //}// unsafe

    (largest_mcc, chosen_txp)
}

pub(super) fn get_num_molecules(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    //ce: &CellExp,
    eqmap: &EqMap,
    tid_to_gid: &[u32],
    log: &slog::Logger
    //exp: &ProcessedExp,
    //gid_map: &HashMap<String, u32, BuildHasherDefault<FNV1aHasher64>>,
    //tgmap: &HashMap<String, String, BuildHasherDefault<FNV1aHasher64>>,
    //num_bootstraps: u32,
    //with_only_gene: Option<String>,
    //is_only_cell: bool,
    //only_unique: bool,
) {
    type U32Set = HashSet<u32, fasthash::RandomState<Hash64>>;
    fn get_set(cap: u32) -> U32Set {
        let s = RandomState::<Hash64>::new();
        U32Set::with_capacity_and_hasher(cap as usize, s)
    }


    // to process only specific genes
    /*
    let mut enable_gene_sampling = false;
    let mut subsample_gene_idx: u32 = <u32>::max_value();
    if !with_only_gene.clone().is_none() {
        enable_gene_sampling = true;
        subsample_gene_idx = *gid_map.get(&with_only_gene.expect("can't unwrap gene name"))
            .expect("can't found index of the given gene");

    }
    */

    let comps = weakly_connected_components(g);
    //println!("{:?}", comps.len());
    let mut one_vertex_components: Vec<usize> = vec![0, 0];
    //let mut identified_txps: Vec<u32> = Vec::with_capacity(comps.len());

    // code for dumping pre collapse graph for sanity check.
    /*
    let dump_no_collapse = false;
    if dump_no_collapse {
        unsafe {
            for (_comp_label, comp_verts) in comps.iter() {
                for vin in comp_verts.iter() {
                    let vert = g.from_index(*vin as usize);
                    for txp in ce
                        .eq_classes
                        .get_unchecked(vert.0 as usize)
                        .transcripts
                        .iter()
                    {
                        eprint!("{},", txp);
                    }
                }
                eprintln!();
            }
        }
    }
    */

    // Make salmon eqclasses for EM based counting
    let s = fasthash::RandomState::<Hash64>::new();
    let mut gene_eqclass_hash: HashMap<Vec<u32>, u32, fasthash::RandomState<Hash64>> =
        HashMap::with_hasher(s);


    // Get the genes that could potentially explain all
    // of the vertices in this mcc.
    // To do this, we first extract the set of _transcripts_
    // that label all vertices of the mcc, and then we project
    // the transcripts to their corresponding gene ids.
    //let mut global_txps : Vec<u32>;
    let mut global_txps  = get_set(16); 

    for (_comp_label, comp_verts) in comps.iter() {
        if comp_verts.len() > 1 {
            // vset will hold the set of vertices that are 
            // covered.
            let mut vset = HashSet::<u32>::from_iter(comp_verts.iter().cloned());

            // we will remove covered vertices from vset until they are 
            // all gone (until all vertices have been covered)
            while !vset.is_empty() {
                // will hold vertices in the best mcc
                let mut best_mcc: Vec<u32> = Vec::new();
                // the transcript that is responsible for the 
                // best mcc covering
                let mut best_covering_txp = std::u32::MAX;
                // for each vertex in the vertex set
                for v in vset.iter() {
                    // find the largest mcc starting from this vertex
                    // and the transcript that covers it
                    let (new_mcc, covering_txp) = collapse_vertices(*v, g, eqmap);

                    // if the new mcc is better than the current best, then 
                    // it becomes the new best
                    if best_mcc.len() < new_mcc.len() {
                        best_mcc = new_mcc;
                        best_covering_txp = covering_txp;
                    }
                }

                // get gene_id of best covering transcript
                let best_covering_gene = tid_to_gid[best_covering_txp as usize];

                if best_covering_txp == std::u32::MAX {
                    crit!(log, "Could not find a covering transcript");
                    std::process::exit(1);
                }

                //unsafe {
                    global_txps.clear();
                    // We iterate over all vertices in the mcc, and for 
                    // each one, we keep track of the (monotonically 
                    // non-increasing) set of transcripts that have appeared
                    // in all vertices.
                    for (index, vertex) in best_mcc.iter().enumerate() {
                        // get the underlying graph vertex for this 
                        // vertex in the mcc
                        let vert = g.from_index((*vertex) as usize);
                        
                        // the first element of the vertex tuple is the 
                        // equivalence class id
                        let eqid = vert.0 as usize;

                        // if this is the first vertex
                        if index == 0 {
                            for lt in eqmap.refs_for_eqc(eqid as u32) {
                                global_txps.insert(*lt);
                            }
                        } else {
                            for lt in eqmap.refs_for_eqc(eqid as u32) {
                                global_txps.remove(lt);
                            }
                        }
                    }
                //} // unsafe

                // project each coverting transcript to it's 
                // corresponding gene, and dedup the list
                let mut global_genes : Vec<u32> = global_txps.iter().cloned().map(|i| tid_to_gid[i as usize]).collect();
                global_genes.sort();
                global_genes.dedup();

                // extract gene level salmon eqclass and increment count by 1
                assert!(
                    !global_genes.is_empty(),
                    "can't find representative gene(s) for a molecule"
                );

                // assert the best covering gene in the global gene set
                //assert!(global_genes.contains(best_covering_gene), 
                //    "best gene not in covering set, shouldn't be possible");
                
                let counter = gene_eqclass_hash.entry(global_genes).or_insert(0);
                *counter += 1;

                // for every vertext that has been covered
                // remove it from vset
                for rv in best_mcc.iter() {
                    vset.remove(rv);
                }
            } //end-while
        } else { // this was a single-vertex subgraph
            let tv = comp_verts.first().expect("can't extract first vertex");
            let tl = eqmap.refs_for_eqc(g.from_index(*tv as usize).0);

            if tl.len() == 1 {
                one_vertex_components[0] += 1;
            } else {
                one_vertex_components[1] += 1;
            }

            let mut global_genes : Vec<u32> = tl.iter().map(|i| tid_to_gid[*i as usize]).collect();
            global_genes.sort();

            // extract gene level salmon eqclass and increment count by 1
            assert!(!global_genes.is_empty(),
                "can't find representative gene(s) for a molecule"
            );

            // incrementing the count of thie eqclass label by 1
            let counter = gene_eqclass_hash.entry(global_genes).or_insert(0);
            *counter += 1;
        }

        //let rand_cover = rand::thread_rng().choose(&tl)
        //    .expect("can;t get random cover");
        //identified_txps.push(*rand_cover as u32);
    }

    /*
    let mut salmon_eqclasses = Vec::<SalmonEQClass>::new();
    for (key, val) in salmon_eqclass_hash {
        salmon_eqclasses.push(SalmonEQClass {
            labels: key,
            counts: val,
        });
    }

    let mut unique_evidence: Vec<bool> = vec![false; gid_map.len()];
    let mut no_ambiguity: Vec<bool> = vec![true; gid_map.len()];
    if is_only_cell {
        info!("Total Networks: {}", comps.len());
        let num_txp_unique_networks: usize = one_vertex_components.iter().sum();
        let num_txp_ambiguous_networks: usize = comps.len() - num_txp_unique_networks;
        info!(
            ">1 vertices Network: {}, {}%",
            num_txp_ambiguous_networks,
            num_txp_ambiguous_networks as f32 * 100.0 / comps.len() as f32
        );
        info!(
            "1 vertex Networks w/ 1 txp: {}, {}%",
            one_vertex_components[0],
            one_vertex_components[0] as f32 * 100.0 / comps.len() as f32
        );
        info!(
            "1 vertex Networks w/ >1 txp: {}, {}%",
            one_vertex_components[1],
            one_vertex_components[1] as f32 * 100.0 / comps.len() as f32
        );

        //info!("Total Predicted Molecules {}", identified_txps.len());

        // iterate and extract gene names
        let mut gene_names: Vec<String> = vec!["".to_string(); gid_map.len()];
        for (gene_name, gene_idx) in gid_map {
            gene_names[*gene_idx as usize] = gene_name.clone();
        }

        if num_bootstraps > 0 {
            //entry point for bootstrapping
            let gene_counts: Vec<Vec<f32>> = do_bootstrapping(salmon_eqclasses, 
                &mut unique_evidence,
                &mut no_ambiguity,
                &num_bootstraps, 
                gid_map.len(),
                only_unique);
        
            write_bootstraps(gene_names, gene_counts, unique_evidence, 
                no_ambiguity, num_bootstraps);
            return None;
        }
        else{
            //entry point for EM
            //println!("{:?}", subsample_gene_idx);
            //println!("{:?}", &salmon_eqclasses);
            let gene_counts: Vec<f32> = optimize(salmon_eqclasses, &mut unique_evidence, 
                &mut no_ambiguity, gid_map.len(), only_unique);

            write_quants(gene_names, gene_counts, unique_evidence, no_ambiguity);
            return None;
        } // end-else
    }
    else {
        Some(optimize(salmon_eqclasses, 
            &mut unique_evidence,
            &mut no_ambiguity,
            gid_map.len(),
            only_unique
        ))
    }
    */
    //identified_txps
}