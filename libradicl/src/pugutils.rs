extern crate fasthash;
extern crate petgraph;
extern crate slog;
extern crate quickersort;

use self::slog::crit;
use fasthash::sea::Hash64;
use fasthash::RandomState;
use std::collections::{HashMap, HashSet, VecDeque};
use std::iter::FromIterator;

use petgraph::prelude::*;
use petgraph::unionfind::*;
use petgraph::visit::NodeIndexable;

use crate::schema::EqMap;

type CCMap = HashMap<u32, Vec<u32>, fasthash::RandomState<Hash64>>;

/// Extract the weakly connected components from the directed graph
/// G.  Interestingly, `petgraph` has a builtin algorithm for returning
/// the strongly-connected components of a digraph, and they have an
/// algorithm for returning the _number_ of connected components of an
/// undirected graph, but no algorithm for returning the actual
/// connected components.  So, we build our own using their union
/// find data structure.  This returns a HashMap, mapping each
/// connected component id (a u32) to the corresponding list of vertex
/// ids (also u32s) contained in the connected component.
pub fn weakly_connected_components<G>(g: G) -> CCMap
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
    fn get_map() -> CCMap {
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

/// Find the largest monochromatic spanning arboresence
/// in the graph `g` starting at vertex `v`.  The arboresence
/// is monochromatic if every vertex can be "covered" by a single
/// transcript (i.e. there exists a transcript that appears in the
/// equivalence class labels of all vertices in the arboresence).
fn collapse_vertices(
    v: u32,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
) -> (Vec<u32>, u32) {
    // get a new set to hold vertices
    type VertexSet = HashSet<u32, fasthash::RandomState<Hash64>>;
    fn get_set(cap: u32) -> VertexSet {
        let s = RandomState::<Hash64>::new();
        VertexSet::with_capacity_and_hasher(cap as usize, s)
    }

    // will hold the nodes in the largest arboresence found
    let mut largest_mcc: Vec<u32> = Vec::new();
    let mut chosen_txp = 0u32;
    let vert = g.from_index(v as usize);

    //unsafe {

    let nvert = g.node_count();

    // for every transcript in the equivalence class
    // label of the vertex
    for txp in eqmap.refs_for_eqc(vert.0).iter() {
        // start a bfs from this vertex
        let mut bfs_list = VecDeque::new();
        bfs_list.push_back(v);

        // the set to remember the nodes we've already
        // visited
        let mut visited_set = get_set(nvert as u32);
        visited_set.insert(v);

        // will hold the current arboresence we
        // are constructing
        let mut current_mcc = Vec::new();

        // get the next vertex in the BFS
        while let Some(cv) = bfs_list.pop_front() {
            // add it to the arboresence
            current_mcc.push(cv);

            // for all of the neighboring vertices that we can
            // reach (those with outgoing, or bidirected edges)
            for nv in g.neighbors_directed(g.from_index(cv as usize), Outgoing) {
                let n = g.to_index(nv) as u32;
                // check if we've seen n in this traversal
                // yet. The `insert()` method returns true
                // if the set didn't have the element, false
                // otherwise.
                if !visited_set.insert(n) {
                    continue;
                }
                // get the set of transcripts present in the
                // label of the current node.
                let n_labels = eqmap.refs_for_eqc(nv.0);
                if let Ok(_n) = n_labels.binary_search(&txp) {
                    bfs_list.push_back(n);
                }
            }
        }

        // if this arboresence is the largest we've yet
        // seen, then record it
        if largest_mcc.len() < current_mcc.len() {
            largest_mcc = current_mcc;
            chosen_txp = *txp;
        }
    }
    //}// unsafe

    (largest_mcc, chosen_txp)
}


pub(super) fn get_num_molecules_trivial(
    eq_map: &EqMap, 
    tid_to_gid: &[u32], 
    num_genes: usize,
    log: &slog::Logger) -> Vec<f32> {

    let mut counts = vec![0.0f32; num_genes];
    let s = RandomState::<Hash64>::new();
    let mut gene_map = HashMap::with_hasher(s);

    for eqinfo in &eq_map.eqc_info {
        let umis = &eqinfo.umis;
        let eqid = &eqinfo.eq_num;
        let tset = eq_map.refs_for_eqc(*eqid);
        let mut prev_gene_id = u32::MAX;
        let mut multi_gene = false;
        // if this is a single-gene equivalence class
        // then go ahead and assign the read
        for t in tset {
            let gid = tid_to_gid[*t as usize];
            if gid != prev_gene_id && prev_gene_id < u32::MAX {
                multi_gene = true;
                break;
            }
            prev_gene_id = gid;
        }
        
        // if the read is single-gene
        // then add this equivalence class' list 
        // of UMIs in the gene map
        if !multi_gene {
            gene_map.entry(prev_gene_id).or_insert(vec![]).extend(umis.iter().map(|x| x.0));
        }
    }

    // go over the map and merge umis from different 
    // equivalence classes that still map to the same
    // gene.
    for (k,v) in gene_map.iter_mut() {
        quickersort::sort(&mut v[..]); v.dedup();
        // the count is the number of distinct UMIs.
        counts[*k as usize] += v.len() as f32;
    }

    // return the counts
    counts
} 

/// Given the digraph `g` representing the PUGs within the current
/// cell, the EqMap `eqmap` to decode all equivalence classes
/// and the transcript-to-gene map `tid_to_gid`, apply the parsimonious
/// umi resolution algorithm.  Pass any relevant logging messages along to
/// `log`.
pub(super) fn get_num_molecules(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    tid_to_gid: &[u32],
    log: &slog::Logger,
) -> HashMap<Vec<u32>, u32, fasthash::RandomState<Hash64>> {
    type U32Set = HashSet<u32, fasthash::RandomState<Hash64>>;
    fn get_set(cap: u32) -> U32Set {
        let s = RandomState::<Hash64>::new();
        U32Set::with_capacity_and_hasher(cap as usize, s)
    }

    let comps = weakly_connected_components(g);
    // a vector of length 2 that records at index 0
    // the number of single-node subgraphs that are
    // transcript-unique and at index 1 the number of
    // single-node subgraphs that have more than one
    // associated transcript.
    let mut one_vertex_components: Vec<usize> = vec![0, 0];

    // Make gene-level eqclasses.
    // This is a map of gene ids to the count of
    // _de-duplicated_ reads observed for that set of genes.
    // For every gene set (label) of length 1, these are gene
    // unique reads.  Standard scRNA-seq counting results
    // can be obtained by simply discarding all equivalence
    // classes of size greater than 1, and probabilistic results
    // will attempt to resolve gene multi-mapping reads by
    // running and EM algorithm.
    let s = fasthash::RandomState::<Hash64>::new();
    let mut gene_eqclass_hash: HashMap<Vec<u32>, u32, fasthash::RandomState<Hash64>> =
        HashMap::with_hasher(s);

    // Get the genes that could potentially explain all
    // of the vertices in this mcc.
    // To do this, we first extract the set of _transcripts_
    // that label all vertices of the mcc, and then we project
    // the transcripts to their corresponding gene ids.
    //let mut global_txps : Vec<u32>;
    let mut global_txps = get_set(16);

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
                    // NOTE: what if there are multiple different mccs that
                    // are equally good? (@k3yavi — I don't think this case
                    // is even handled in the C++ code either).
                    let (new_mcc, covering_txp) = collapse_vertices(*v, g, eqmap);

                    // if the new mcc is better than the current best, then
                    // it becomes the new best
                    if best_mcc.len() < new_mcc.len() {
                        best_mcc = new_mcc;
                        best_covering_txp = covering_txp;
                    }
                }

                if best_covering_txp == std::u32::MAX {
                    crit!(log, "Could not find a covering transcript");
                    std::process::exit(1);
                }

                // get gene_id of best covering transcript
                let best_covering_gene = tid_to_gid[best_covering_txp as usize];

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
                        //crit!(log, "global txps = {:#?}\ncurr refs = {:#?}", global_txps, eqmap.refs_for_eqc(eqid as u32));
                        let txps_for_vert = eqmap.refs_for_eqc(eqid as u32);
                        global_txps.retain(|t| txps_for_vert.binary_search(t).is_ok());
                        //for lt in eqmap.refs_for_eqc(eqid as u32) {
                        //    global_txps.remove(lt);
                        //}
                    }
                }
                // at this point, whatever transcript ids remain in
                // global_txps appear in all vertices of the mcc

                //} // unsafe

                // project each coverting transcript to it's
                // corresponding gene, and dedup the list
                let mut global_genes: Vec<u32> = global_txps
                    .iter()
                    .cloned()
                    .map(|i| tid_to_gid[i as usize])
                    .collect();
                // sort since we will be hashing the ordered vector
                quickersort::sort(&mut global_genes[..]);
                // dedup as well since we don't care about duplicates
                global_genes.dedup();

                // assert the best covering gene in the global gene set
                assert!(
                    global_genes.contains(&best_covering_gene),
                    "best gene {} not in covering set, shouldn't be possible",
                    best_covering_gene
                );

                assert!(
                    !global_genes.is_empty(),
                    "can't find representative gene(s) for a molecule"
                );

                // in our hash, increment the count of this equivalence class
                // by 1 (and insert it if we've not seen it yet).
                let counter = gene_eqclass_hash.entry(global_genes).or_insert(0);
                *counter += 1;

                // for every vertext that has been covered
                // remove it from vset
                for rv in best_mcc.iter() {
                    vset.remove(rv);
                }
            } //end-while
        } else {
            // this was a single-vertex subgraph
            let tv = comp_verts.first().expect("can't extract first vertex");
            let tl = eqmap.refs_for_eqc(g.from_index(*tv as usize).0);

            if tl.len() == 1 {
                one_vertex_components[0] += 1;
            } else {
                one_vertex_components[1] += 1;
            }

            let mut global_genes: Vec<u32> = tl.iter().map(|i| tid_to_gid[*i as usize]).collect();
            quickersort::sort(&mut global_genes[..]);
            global_genes.dedup();

            // extract gene-level eqclass and increment count by 1
            assert!(
                !global_genes.is_empty(),
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

    gene_eqclass_hash
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
