/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

#[allow(unused_imports)]
use ahash::{AHasher, RandomState};
use arrayvec::ArrayVec;
use smallvec::SmallVec;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};
use std::io;
use std::io::Write;

use petgraph::prelude::*;
use petgraph::unionfind::*;
use petgraph::visit::NodeIndexable;

use libradicl::chunk;
use libradicl::record::{ConvertiblePrimitiveInteger, 
    MappedRecord, CollatableMappedRecord, KnownSize, 
    UmiTaggedRecord, RecordContext,AlevinFryReadRecordWithPosition,
};

use slog::{crit, info, warn};

use crate::eq_class::{EqMap, EqMapType};
use crate::quant::SplicedAmbiguityModel;
use crate::utils as afutils;

type CcMap = HashMap<u32, Vec<u32>, ahash::RandomState>;

// Forseti/RAD alignment tuple: (dir, ref_start)
type AlgnTuple = (bool, u32);

use crate::forseti::forseti_for_multi_best;
use ndarray::Array1;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use tch::nn;

#[derive(Debug)]
pub enum PugEdgeType {
    NoEdge,
    BiDirected,
    XToY,
    YToX,
}

#[derive(Debug)]
pub struct PugResolutionStatistics {
    pub used_alternative_strategy: bool,
    pub total_mccs: u64,
    pub ambiguous_mccs: u64,
    pub trivial_mccs: u64,
    // Forseti-specific / debug counters (kept as no-op for non-Forseti paths)
    pub non_triv_mccs_determined_ss: (u64, u64),
    pub non_triv_single_gene_ambiguous_ss: u64,
    pub large_mcc_forseti_needed: u64,
    pub triv_mcc_check: u64,
    pub triv_mccs_determined_ss: (u64, u64),
    pub triv_single_gene_ambiguous_ss: u64,
    pub triv_mcc_forseti_needed: u64,
}

/// Allowed suffix categories for transcript names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SuffixType {
    NoneDash, // No suffix tag (no "-X" at the end)
    T,        // "-T" (anothe type for unspliced)
    S,        // "-S" (spliced)
    U,        // "-U" (unspliced)
}

/// Extracts the parsimonious UMI graphs (PUGs) from the
/// equivalence class map for a given cell.
/// The returned graph is a directed graph (potentially with
/// bidirected edges) where each node consists of an (equivalence
/// class, UMI ID) pair.  Note, crucially, that the UMI ID is simply
/// the rank of the UMI in the list of all distinct UMIs for this
/// equivalence class.  There is a directed edge between any pair of
/// vertices whose set of transcripts overlap and whose UMIs are within
/// a Hamming distance of 1 of each other.  If one node has more than
/// twice the frequency of the other, the edge is directed from the
/// more frequent to the less freuqent node.  Otherwise, edges are
/// added in both directions.
pub fn extract_graph(
    eqmap: &EqMap,
    pug_exact_umi: bool, // true if only identical UMIs induce an edge
    log: &slog::Logger,
) -> petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed> {
    let verbose = false;
    let mut one_edit = 0u64;
    let mut zero_edit = 0u64;

    // given 2 pairs (UMI, count), determine if an edge exists
    // between them, and if so, what type.
    let mut has_edge = |x: &(u64, u32), y: &(u64, u32)| -> PugEdgeType {
        let hdist = if pug_exact_umi {
            if x.0 == y.0 {
                0
            } else {
                usize::MAX
            }
        } else {
            afutils::count_diff_2_bit_packed(x.0, y.0)
        };

        if hdist == 0 {
            zero_edit += 1;
            return PugEdgeType::BiDirected;
        }

        if hdist < 2 {
            one_edit += 1;
            return if x.1 > (2 * y.1 - 1) {
                PugEdgeType::XToY
            } else if y.1 > (2 * x.1 - 1) {
                PugEdgeType::YToX
            } else {
                PugEdgeType::BiDirected
            };
        }
        PugEdgeType::NoEdge
    };

    let mut _bidirected = 0u64;
    let mut _unidirected = 0u64;

    let mut graph = DiGraphMap::<(u32, u32), ()>::new();
    let mut hset = vec![0u8; eqmap.num_eq_classes()];
    let mut idxvec: SmallVec<[u32; 128]> = SmallVec::new();

    // insert all of the nodes up front to avoid redundant
    // checks later.
    for eqid in 0..eqmap.num_eq_classes() {
        // get the info Vec<(UMI, frequency)>
        let eq = &eqmap.eqc_info[eqid];
        let u1 = &eq.umis;
        for (xi, _x) in u1.iter().enumerate() {
            graph.add_node((eqid as u32, xi as u32));
        }
    }

    // for every equivalence class in this cell
    for eqid in 0..eqmap.num_eq_classes() {
        if verbose && eqid % 1000 == 0 {
            print!("\rprocessed {:?} eq classes", eqid);
            io::stdout().flush().expect("Could not flush stdout");
        }

        // get the info Vec<(UMI, frequency)>
        let eq = &eqmap.eqc_info[eqid];

        // for each (umi, count) pair and its index
        let u1 = &eq.umis;
        for (xi, x) in u1.iter().enumerate() {
            // add a node
            // graph.add_node((eqid as u32, xi as u32));

            // for each (umi, freq) pair and node after this one
            for (xi2, x2) in u1.iter().enumerate().skip(xi + 1) {
                //for xi2 in (xi + 1)..u1.len() {
                // x2 is the other (umi, freq) pair
                //let x2 = &u1[xi2];

                // add a node for it
                // graph.add_node((eqid as u32, xi2 as u32));

                // determine if an edge exists between x and x2, and if so, what kind
                let et = has_edge(x, x2);
                // for each type of edge, add the appropriate edge in the graph
                match et {
                    PugEdgeType::BiDirected => {
                        graph.add_edge((eqid as u32, xi as u32), (eqid as u32, xi2 as u32), ());
                        graph.add_edge((eqid as u32, xi2 as u32), (eqid as u32, xi as u32), ());
                        _bidirected += 1;
                        //if multi_gene_vec[eqid] == true {
                        //    bidirected_in_multigene += 1;
                        //}
                    }
                    PugEdgeType::XToY => {
                        graph.add_edge((eqid as u32, xi as u32), (eqid as u32, xi2 as u32), ());
                        _unidirected += 1;
                        //if multi_gene_vec[eqid] == true {
                        //    unidirected_in_multigene += 1;
                        //}
                    }
                    PugEdgeType::YToX => {
                        graph.add_edge((eqid as u32, xi2 as u32), (eqid as u32, xi as u32), ());
                        _unidirected += 1;
                        //if multi_gene_vec[eqid] == true {
                        //    unidirected_in_multigene += 1;
                        //}
                    }
                    PugEdgeType::NoEdge => {}
                }
            }
        }

        //hset.clear();
        //hset.resize(eqmap.num_eq_classes(), 0u8);
        for i in &idxvec {
            hset[*i as usize] = 0u8;
        }
        let stf = idxvec.len() > 128;
        idxvec.clear();
        if stf {
            idxvec.shrink_to_fit();
        }

        // for every reference id in this eq class
        for r in eqmap.refs_for_eqc(eqid as u32) {
            // find the equivalence classes sharing this reference
            for eq2id in eqmap.eq_classes_containing(*r).iter() {
                // if eq2id <= eqid, then we already observed the relevant edges
                // when we process eq2id
                if (*eq2id as usize) <= eqid {
                    continue;
                }
                // otherwise, if we have already processed this other equivalence
                // class because it shares _another_ reference (apart from r) with
                // the current equivalence class, then skip it.
                if hset[*eq2id as usize] > 0 {
                    continue;
                }

                // recall that we processed this eq class as a neighbor of eqid
                hset[*eq2id as usize] = 1;
                idxvec.push(*eq2id);
                let eq2 = &eqmap.eqc_info[*eq2id as usize];

                // compare all the umis between eqid and eq2id
                let u2 = &eq2.umis;
                for (xi, x) in u1.iter().enumerate() {
                    // Node for equiv : eqid and umi : xi
                    // graph.add_node((eqid as u32, xi as u32));

                    for (yi, y) in u2.iter().enumerate() {
                        // Node for equiv : eq2id and umi : yi
                        // graph.add_node((*eq2id as u32, yi as u32));

                        let et = has_edge(x, y);
                        match et {
                            PugEdgeType::BiDirected => {
                                graph.add_edge((eqid as u32, xi as u32), (*eq2id, yi as u32), ());
                                graph.add_edge((*eq2id, yi as u32), (eqid as u32, xi as u32), ());
                                _bidirected += 1;
                                //if multi_gene_vec[eqid] == true
                                //    || multi_gene_vec[*eq2id as usize] == true
                                //{
                                //    bidirected_in_multigene += 1;
                                //}
                            }
                            PugEdgeType::XToY => {
                                graph.add_edge((eqid as u32, xi as u32), (*eq2id, yi as u32), ());
                                _unidirected += 1;
                                //if multi_gene_vec[eqid] == true
                                //    || multi_gene_vec[*eq2id as usize] == true
                                //{
                                //    unidirected_in_multigene += 1;
                                //}
                            }
                            PugEdgeType::YToX => {
                                graph.add_edge((*eq2id, yi as u32), (eqid as u32, xi as u32), ());
                                _unidirected += 1;
                                //if multi_gene_vec[eqid] == true
                                //    || multi_gene_vec[*eq2id as usize] == true
                                //{
                                //    unidirected_in_multigene += 1;
                                //}
                            }
                            PugEdgeType::NoEdge => {}
                        }
                    }
                }
            }
        }
    }

    if verbose {
        info!(
            log,
            "\n\nsize of graph ({:?}, {:?})\n\n",
            graph.node_count(),
            graph.edge_count()
        );
        let total_edits = (one_edit + zero_edit) as f64;
        info!(log, "\n\n\n{}\n\n\n", one_edit as f64 / total_edits);
    }

    graph
}

/// Extract the weakly connected components from the directed graph
/// G.  Interestingly, `petgraph` has a builtin algorithm for returning
/// the strongly-connected components of a digraph, and they have an
/// algorithm for returning the _number_ of connected components of an
/// undirected graph, but no algorithm for returning the actual
/// connected components.  So, we build our own using their union
/// find data structure.  This returns a HashMap, mapping each
/// connected component id (a u32) to the corresponding list of vertex
/// ids (also u32s) contained in the connected component.
pub fn weakly_connected_components<G>(g: G) -> CcMap
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
    fn get_map() -> CcMap {
        let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
        HashMap::with_hasher(s)
    }

    let mut components = get_map();
    for (i, v) in labels.iter().enumerate() {
        let ve = components.entry(*v as u32).or_default();
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
    uncovered_vertices: &HashSet<u32, ahash::RandomState>, // the set of vertices already covered
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
) -> (Vec<u32>, u32) {
    // get a new set to hold vertices
    type VertexSet = HashSet<u32, ahash::RandomState>;
    let get_set =
        |cap: u32| VertexSet::with_capacity_and_hasher(cap as usize, hasher_state.clone());

    // will hold the nodes in the largest arboresence found
    let mut largest_mcc: Vec<u32> = Vec::new();
    let mut chosen_txp = 0u32;
    let vert = g.from_index(v as usize);
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

                // check if we should add this vertex or not:
                // uncovered_vertices contains the the current set of
                // *uncovered* vertices in this component (i.e. those
                // that we still need to explain by some molecule).
                // so, if n is *not* in uncovered_vertices, then it is not in the
                // uncovered set, and so it has already been
                // explained / covered.
                //
                // if n hasn't been covered yet, then
                // check if we've seen n in this traversal
                // yet. The `insert()` method returns true
                // if the set didn't have the element, false
                // otherwise.
                if !uncovered_vertices.contains(&n) || !visited_set.insert(n) {
                    continue;
                }

                // get the set of transcripts present in the
                // label of the current node.
                if let Ok(_n) = eqmap.refs_for_eqc(nv.0).binary_search(txp) {
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

fn collapse_vertices_keep_ties(
    v: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
) -> HashMap<Vec<u32>, Vec<u32>, ahash::RandomState> {
    type VertexSet = HashSet<u32, ahash::RandomState>;
    let get_set = |cap: u32| VertexSet::with_capacity_and_hasher(cap as usize, hasher_state.clone());

    let mut all_mcc_txp_map: HashMap<Vec<u32>, Vec<u32>, ahash::RandomState> = HashMap::with_capacity_and_hasher(g.node_count() as usize, hasher_state.clone());
    
    let mut max_len = 0usize;
    let vert = g.from_index(v as usize);

    for txp in eqmap.refs_for_eqc(vert.0).iter() {
        let mut bfs_list = VecDeque::new();
        bfs_list.push_back(v);
        let mut visited_set = get_set(g.node_count() as u32);
        visited_set.insert(v);

        let mut current_mcc = Vec::new();
        while let Some(cv) = bfs_list.pop_front() {
            current_mcc.push(cv);
            for nv in g.neighbors_directed(g.from_index(cv as usize), Outgoing) {
                let n = g.to_index(nv) as u32;
                if !uncovered_vertices.contains(&n) || !visited_set.insert(n) {
                    continue;
                }
                if let Ok(_) = eqmap.refs_for_eqc(nv.0).binary_search(txp) {
                    bfs_list.push_back(n);
                }
            }
        }

        let current_len = current_mcc.len();
        if current_len >= max_len && current_len > 0 {
            if current_len > max_len {
                max_len = current_len;
                all_mcc_txp_map.clear();
            }
            // sort the mcc to ensure the key is unique
            current_mcc.sort_unstable(); 
            all_mcc_txp_map.entry(current_mcc).or_default().push(*txp);
        }
    }

    // dedup the txp list in each mcc
    for txps in all_mcc_txp_map.values_mut() {
        txps.sort_unstable();
        txps.dedup();
    }

    all_mcc_txp_map
}

/// From a start vertex `v`, enumerate ALL maximum-size MCCs (ties kept).
/// Each returned MCC is sorted, so it can be used as a dedup key.
/// The txp_list is the set of transcripts that generate exactly that MCC
/// from this start vertex (i.e., compatible with all vertices in that MCC
/// under this BFS rule).
pub fn collapse_vertices_multi(
    v: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
) -> Vec<(Vec<u32>, Vec<u32>)> {
    type VertexSet = HashSet<u32, ahash::RandomState>;
    let get_set =
        |cap: usize| VertexSet::with_capacity_and_hasher(cap, hasher_state.clone());

    // (eqid, umi_idx) for this vertex
    let vert = g.from_index(v as usize);
    let nvert = g.node_count();
    let mut best_len: usize = 0;

    // Map sorted_mcc -> txp_list
    let mut best_map: HashMap<Vec<u32>, Vec<u32>, ahash::RandomState> =
        HashMap::with_capacity_and_hasher(8, hasher_state.clone());

    // Try each transcript in the starting vertex label
    for &txp in eqmap.refs_for_eqc(vert.0).iter() {
        let mut bfs = VecDeque::new();
        bfs.push_back(v);

        let mut visited = get_set(nvert);
        visited.insert(v);

        let mut current_mcc: Vec<u32> = Vec::new();

        while let Some(cv) = bfs.pop_front() {
            current_mcc.push(cv);

            // Traverse outgoing neighbors (same as original collapse_vertices)
            for nv in g.neighbors_directed(g.from_index(cv as usize), Outgoing) {
                let n = g.to_index(nv) as u32;

                // only consider uncovered vertices, and don't revisit
                if !uncovered_vertices.contains(&n) || !visited.insert(n) {
                    continue;
                }

                // Check if txp is present in this neighbor's eqc labels
                // NOTE: requires refs_for_eqc(...) sorted
                let n_labels = eqmap.refs_for_eqc(nv.0);
                if n_labels.binary_search(&txp).is_ok() {
                    bfs.push_back(n);
                }
            }
        }

        let cur_len = current_mcc.len();
        if cur_len == 0 {
            continue;
        }

        // Canonicalize MCC as sorted key
        current_mcc.sort_unstable();

        if cur_len > best_len {
            best_len = cur_len;
            best_map.clear();
            best_map.insert(current_mcc, vec![txp]);
        } else if cur_len == best_len {
            // Merge txp into the txp_list for this MCC key
            best_map
                .entry(current_mcc)
                .or_insert_with(|| Vec::new())
                .push(txp);
        }
        // if cur_len < best_len: ignore
    }

    // Convert map -> vec, and sort/dedup txp_list in each entry
    let mut out: Vec<(Vec<u32>, Vec<u32>)> = Vec::with_capacity(best_map.len());
    for (mcc, mut txps) in best_map {
        txps.sort_unstable();
        txps.dedup();
        out.push((mcc, txps));
    }

    // Optional: stable-ish ordering for deterministic behavior
    // (sort by mcc length desc (all same), then lexicographic)
    out.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    out
}

#[inline]
fn resolve_num_molecules_crlike_from_vec_prefer_ambig(
    umi_gene_count_vec: &mut [(u64, u32, u32)],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
) {
    // sort the triplets
    // first on umi
    // then on gene_id
    // then on count
    umi_gene_count_vec.sort_unstable();

    // hold the current umi and gene we are examining
    let mut curr_umi = umi_gene_count_vec.first().expect("cell with no UMIs").0;
    let first_gn = umi_gene_count_vec.first().expect("cell with no UMIs").1;
    // The capacity of curr_gn is 2 as it will be used to hold the
    // the spliced id of a gene, the unspliced id of a gene, or both
    let mut curr_gn = ArrayVec::<u32, 2>::new();
    curr_gn.push(first_gn);

    // hold the gene id having the max count for this umi
    // and the maximum count value itself
    let mut max_count = 0u32;
    // to aggregate the count should a (umi, gene) pair appear
    // more than once
    let mut count_aggr = 0u32;
    // the vector will hold the equivalent set of best genes
    let mut best_genes = Vec::<u32>::with_capacity(16);

    // look over all sorted triplets
    for (cidx, &(umi, gn, ct)) in umi_gene_count_vec.iter().enumerate() {
        // if this umi is different than
        // the one we are processing
        // then decide what action to take
        // on the previous umi
        if umi != curr_umi {
            // update the count of the equivalence class of genes
            // that gets this UMI

            *gene_eqclass_hash.entry(best_genes.clone()).or_insert(0) += 1;

            // the next umi and gene
            curr_umi = umi;
            curr_gn.clear();
            curr_gn.push(gn);

            // current gene is current best
            best_genes.clear();
            best_genes.push(gn);

            // count aggr = max count = ct
            count_aggr = ct;
            max_count = ct;
        } else {
            // the umi was the same

            let prev_gid = *curr_gn.last().expect("not empty");

            // if the gene is the same (modulo splicing), add the counts
            if afutils::same_gene(gn, prev_gid, true) {
                // if the gene is the same modulo splicing, but
                // curr_gn != gn, then this is gene will be set as ambiguous,
                // this is the transition from counting occurrences
                // of this UMI from the spliced to unspliced version
                // of the gene.
                if prev_gid != gn {
                    // mark the current gene as
                    // splicing ambiguous by pushing
                    // the unspliced id onto curr_gn
                    curr_gn.push(gn);
                }
                count_aggr += ct;
            } else {
                // if the gene is different, then restart the count_aggr
                // and set the current gene id
                count_aggr = ct;
                curr_gn.clear();
                curr_gn.push(gn);
            }

            // we have the following cases, consider we are
            // processing the records for gene g_i.  If
            // * curr_gn = [g_i^s], then we have so far only observed spliced reads for g_i
            // * curr_gn = [g_i^u], then we have only observed unspliced reads for g_i
            //   (and there are no spliced reads)
            // * curr_gn = [g_i^s, g_i^u], then we have observed both spliced and unspliced
            //   reads for g_i, and this gene will be considered splicing ambiguous.

            // the current best_genes is either empty or contains some
            // set of genes g_i-k, ..., g_i-1, and now,
            // g_i matches their count.  If g_i has only observed spliced reads,
            // then we will add g_i^s to the list.  If g_i has observed
            // only unspliced reads then we will add g_i^u, otherwise
            // we will add both g_i^s and g_i^u.

            // the current best_genes is either emtpy or contains some
            // set of genes g_i-k, ..., g_i-1, and now,
            // g_i *exceeds* their count.  If g_i has only observed spliced reads,
            // then we will *clear* best_genes, and replace it with g_i^s.
            // If g_i has only observed unspliced reads, then we will *clear*
            // best_genes and replace it with g_i^u.  Othewise we will *clear*
            // best_genes and replace it with [g_i^s, g_i^u];

            // if the count aggregator exceeded the max
            // then it is the new max, and this gene is
            // the new max gene.  Having a distinct max
            // also makes this UMI uniquely resolvable
            match count_aggr.cmp(&max_count) {
                Ordering::Greater => {
                    max_count = count_aggr;
                    best_genes.clear();
                    best_genes.extend(curr_gn.iter());
                }
                Ordering::Equal => {
                    // if we have a tie for the max count
                    // then the current UMI isn't uniquely-unresolvable
                    // it will stay this way unless we see a bigger
                    // count for this UMI.  We add the current
                    // "tied" gene to the equivalence class.
                    best_genes.extend(curr_gn.iter());
                }
                Ordering::Less => {
                    // we do nothing
                }
            }
        }

        // if this was the last UMI in the list
        if cidx == umi_gene_count_vec.len() - 1 {
            *gene_eqclass_hash.entry(best_genes.clone()).or_insert(0) += 1;
        }
    }
}

#[inline]
fn resolve_num_molecules_crlike_from_vec(
    umi_gene_count_vec: &mut [(u64, u32, u32)],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
) {
    // sort the triplets
    // first on umi
    // then on gene_id
    // then on count
    umi_gene_count_vec.sort_unstable();

    // hold the current umi and gene we are examining
    let mut curr_umi = umi_gene_count_vec.first().expect("cell with no UMIs").0;
    let mut curr_gn = umi_gene_count_vec.first().expect("cell with no UMIs").1;
    // hold the gene id having the max count for this umi
    // and the maximum count value itself
    // let mut max_count_gene = 0u32;
    let mut max_count = 0u32;
    // to aggregate the count should a (umi, gene) pair appear
    // more than once
    let mut count_aggr = 0u32;
    // the vector will hold the equivalent set of best genes
    let mut best_genes = Vec::<u32>::with_capacity(16);

    // look over all sorted triplets
    for (cidx, &(umi, gn, ct)) in umi_gene_count_vec.iter().enumerate() {
        // if this umi is different than
        // the one we are processing
        // then decide what action to take
        // on the previous umi
        if umi != curr_umi {
            // update the count of the equivalence class of genes
            // that gets this UMI
            *gene_eqclass_hash.entry(best_genes.clone()).or_insert(0) += 1;

            // the next umi and gene
            curr_umi = umi;
            curr_gn = gn;

            // current gene is current best
            // max_count_gene = gn;
            best_genes.clear();
            best_genes.push(gn);

            // count aggr = max count = ct
            count_aggr = ct;
            max_count = ct;
        } else {
            // the umi was the same

            // if the gene is the same, add the counts
            if gn == curr_gn {
                count_aggr += ct;
            } else {
                // if the gene is different, then restart the count_aggr
                // and set the current gene id
                count_aggr = ct;
                curr_gn = gn;
            }
            // if the count aggregator exceeded the max
            // then it is the new max, and this gene is
            // the new max gene.  Having a distinct max
            // also makes this UMI uniquely resolvable
            match count_aggr.cmp(&max_count) {
                Ordering::Greater => {
                    max_count = count_aggr;
                    // we want to avoid the case that we are just
                    // updating the count of the best gene above and
                    // here we clear out the vector and populate it
                    // with the same element again and again.  So
                    // if the current best_genes vector holds just
                    // gn, we do nothing.  Otherwise we clear it and
                    // add gn.
                    match &best_genes[..] {
                        [x] if *x == gn => { /* do nothing here */ }
                        _ => {
                            best_genes.clear();
                            best_genes.push(gn);
                        }
                    }
                }
                Ordering::Equal => {
                    // if we have a tie for the max count
                    // then the current UMI isn't uniquely-unresolvable
                    // it will stay this way unless we see a bigger
                    // count for this UMI.  We add the current
                    // "tied" gene to the equivalence class.
                    best_genes.push(gn);
                }
                Ordering::Less => {
                    // we do nothing
                }
            }
        }

        // if this was the last UMI in the list
        if cidx == umi_gene_count_vec.len() - 1 {
            *gene_eqclass_hash.entry(best_genes.clone()).or_insert(0) += 1;
        }
    }
}

pub fn get_num_molecules_cell_ranger_like_small<B, R>(
    cell_chunk: &mut chunk::Chunk<R>,
    tid_to_gid: &[u32],
    _num_genes: usize,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
    sa_model: SplicedAmbiguityModel,
    _log: &slog::Logger,
) where
B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize + UmiTaggedRecord, 
    <R as MappedRecord>::ParsingContext: RecordContext, 
    <R as MappedRecord>::ParsingContext: Clone,
    <R as MappedRecord>::ParsingContext: Send
{
    let mut umi_gene_count_vec: Vec<(u64, u32, u32)> = Vec::with_capacity(cell_chunk.nrec as usize);

    // for each record
    for rec in &cell_chunk.reads {
        // get the umi
        let umi = rec.umi();

        // project the transcript ids to gene ids
        let mut gset: Vec<u32> = rec
            .refs()
            .iter()
            .map(|tid| tid_to_gid[*tid as usize])
            .collect();
        // and make the gene ids unique
        gset.sort_unstable();
        gset.dedup();
        for g in &gset {
            umi_gene_count_vec.push((umi, *g, 1));
        }
    }
    match sa_model {
        SplicedAmbiguityModel::WinnerTakeAll => {
            resolve_num_molecules_crlike_from_vec(&mut umi_gene_count_vec, gene_eqclass_hash);
        }
        SplicedAmbiguityModel::PreferAmbiguity => {
            resolve_num_molecules_crlike_from_vec_prefer_ambig(
                &mut umi_gene_count_vec,
                gene_eqclass_hash,
            );
        }
    }
}

pub fn get_num_molecules_cell_ranger_like(
    eq_map: &EqMap,
    tid_to_gid: &[u32],
    _num_genes: usize,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
    sa_model: SplicedAmbiguityModel,
    _log: &slog::Logger,
) {
    // TODO: better capacity
    let mut umi_gene_count_vec: Vec<(u64, u32, u32)> = vec![];

    // for each equivalence class
    for eqinfo in &eq_map.eqc_info {
        // get the (umi, count) pairs
        let umis = &eqinfo.umis;
        let eqid = &eqinfo.eq_num;

        // project the transcript ids to gene ids
        let mut gset: Vec<u32> = eq_map
            .refs_for_eqc(*eqid)
            .iter()
            .map(|tid| tid_to_gid[*tid as usize])
            .collect();
        // and make the gene ids unique,
        // note, if we have both spliced and
        // unspliced gene ids, then they will
        // necessarily be adjacent here, since
        // they are always asigned adjacent ids
        // with spliced being even and unspliced odd.
        gset.sort_unstable();
        gset.dedup();

        // add every (umi, count), gene pair as a triplet
        // of (umi, gene_id, count) to the output vector
        for umi_ct in umis {
            for g in &gset {
                umi_gene_count_vec.push((umi_ct.0, *g, umi_ct.1));
            }
        }
    }
    match sa_model {
        SplicedAmbiguityModel::WinnerTakeAll => {
            resolve_num_molecules_crlike_from_vec(&mut umi_gene_count_vec, gene_eqclass_hash);
        }
        SplicedAmbiguityModel::PreferAmbiguity => {
            resolve_num_molecules_crlike_from_vec_prefer_ambig(
                &mut umi_gene_count_vec,
                gene_eqclass_hash,
            );
        }
    }
}

pub fn get_num_molecules_trivial_discard_all_ambig(
    eq_map: &EqMap,
    tid_to_gid: &[u32],
    num_genes: usize,
    _log: &slog::Logger,
) -> (Vec<f32>, f64) {
    let mut counts = vec![0.0f32; num_genes];
    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut gene_map: std::collections::HashMap<u32, Vec<u64>, ahash::RandomState> =
        HashMap::with_hasher(s);

    let mut total_umis = 0u64;
    let mut multi_gene_umis = 0u64;

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

        total_umis += umis.len() as u64;
        if multi_gene {
            multi_gene_umis += umis.len() as u64;
        }

        // if the read is single-gene
        // then add this equivalence class' list
        // of UMIs in the gene map
        if !multi_gene {
            gene_map
                .entry(prev_gene_id)
                .or_default()
                .extend(umis.iter().map(|x| x.0));
        }
    }

    // go over the map and merge umis from different
    // equivalence classes that still map to the same
    // gene.
    for (k, v) in gene_map.iter_mut() {
        v.sort_unstable();
        v.dedup();
        // the count is the number of distinct UMIs.
        counts[*k as usize] += v.len() as f32;
    }

    // return the counts
    (counts, multi_gene_umis as f64 / total_umis as f64)
}

/// given the connected component (subgraph) of `g` defined by the
/// vertices in `vertex_ids`, apply the cell-ranger-like algorithm
/// within this subgraph.
fn get_num_molecules_large_component(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eq_map: &EqMap,
    vertex_ids: &[u32],
    tid_to_gid: &[u32],
    hasher_state: &ahash::RandomState,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
    _log: &slog::Logger,
) {
    let gene_level_eq_map = match eq_map.map_type {
        EqMapType::GeneLevel => true,
        EqMapType::TranscriptLevel => false,
    };

    // TODO: better capacity
    let mut umi_gene_count_vec: Vec<(u64, u32, u32)> = vec![];

    // build a temporary hashmap from each
    // equivalence class id in the current subgraph
    // to the set of (UMI, frequency) pairs contained
    // in the subgraph
    //let ts = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut tmp_map =
        HashMap::<u32, Vec<(u64, u32)>, ahash::RandomState>::with_hasher(hasher_state.clone());

    // for each vertex id in the subgraph
    for vertex_id in vertex_ids {
        // get the corresponding vertex which is
        // an (eq_id, UMI index) pair
        let vert = g.from_index(*vertex_id as usize);
        // add the corresponding (UMI, frequency) pair to the map
        // for this eq_id
        let umis = tmp_map.entry(vert.0).or_default();
        umis.push(eq_map.eqc_info[vert.0 as usize].umis[vert.1 as usize]);
    }

    for (k, v) in tmp_map.iter() {
        // get the (umi, count) pairs
        let umis = v; //&eqinfo.umis;
        let eqid = k; //&eqinfo.eq_num;
                    // project the transcript ids to gene ids
        let mut gset: Vec<u32>;

        if gene_level_eq_map {
            gset = eq_map.refs_for_eqc(*eqid).to_vec();
        } else {
            gset = eq_map
                .refs_for_eqc(*eqid)
                .iter()
                .map(|tid| tid_to_gid[*tid as usize])
                .collect();
            // and make the gene ids unique
            gset.sort_unstable();
            gset.dedup();
        }

        // add every (umi, count), gene pair as a triplet
        // of (umi, gene_id, count) to the output vector
        for umi_ct in umis {
            for g in &gset {
                umi_gene_count_vec.push((umi_ct.0, *g, umi_ct.1));
            }
        }
    }

    resolve_num_molecules_crlike_from_vec(&mut umi_gene_count_vec, gene_eqclass_hash);
}

/// Given the digraph `g` representing the PUGs within the current
/// cell, the EqMap `eqmap` to decode all equivalence classes
/// and the transcript-to-gene map `tid_to_gid`, apply the parsimonious
/// umi resolution algorithm.  Pass any relevant logging messages along to
/// `log`.
pub fn get_num_molecules(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    tid_to_gid: &[u32],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
    hasher_state: &ahash::RandomState,
    large_graph_thresh: usize,
    log: &slog::Logger,
) -> PugResolutionStatistics
{
    type U32Set = HashSet<u32, ahash::RandomState>;
    let get_set = |cap: u32| {
        //let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
        U32Set::with_capacity_and_hasher(cap as usize, hasher_state.clone())
    };

    let gene_level_eq_map = match eqmap.map_type {
        EqMapType::GeneLevel => true,
        EqMapType::TranscriptLevel => false,
    };

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
    //let s = fasthash::RandomState::<Hash64>::new();
    //let mut gene_eqclass_hash: HashMap<Vec<u32>, u32, fasthash::RandomState<Hash64>> =
    //    HashMap::with_hasher(s);

    // Get the genes that could potentially explain all
    // of the vertices in this mcc.
    // To do this, we first extract the set of _transcripts_
    // that label all vertices of the mcc, and then we project
    // the transcripts to their corresponding gene ids.
    //let mut global_txps : Vec<u32>;
    let mut global_txps = get_set(16);
    let mut pug_stats = PugResolutionStatistics {
        used_alternative_strategy: false,
        total_mccs: 0u64,
        ambiguous_mccs: 0u64,
        trivial_mccs: 0u64,
        non_triv_mccs_determined_ss: (0u64, 0u64),
        non_triv_single_gene_ambiguous_ss: 0u64,
        large_mcc_forseti_needed: 0u64,
        triv_mcc_check: 0u64,
        triv_mccs_determined_ss: (0u64, 0u64),
        triv_single_gene_ambiguous_ss: 0u64,
        triv_mcc_forseti_needed: 0u64,
    };

    for (_comp_label, comp_verts) in comps.iter() {
        if comp_verts.len() > 1 {
            // the current parsimony resolution algorithm
            // can become slow for connected components that
            // are very large.  For components with > large_graph_thresh
            // vertices (this should be _very_ rare) we will instead
            // resolve the UMIs in the component using a simpler algorithm.
            if comp_verts.len() > large_graph_thresh {
                get_num_molecules_large_component(
                    g,
                    eqmap,
                    comp_verts,
                    tid_to_gid,
                    hasher_state,
                    gene_eqclass_hash,
                    log,
                );
                warn!(
                    log,
                    "found connected component with {} vertices; resolved with cr-like resolution.",
                    comp_verts.len(),
                );
                pug_stats.used_alternative_strategy = true;
                continue;
            }

            // uncovered_vertices will hold the set of vertices that are
            // *not yet* covered.
            //
            // Non-deterministic variant :
            // NOTE: The line below places the vertices into a HashSet that uses
            // a hasher with a RandomState which, by default, Rust will randomize between
            // runs.  That means that the output of the entire algorithm will, in general,
            // not be deterministic.  By using a RandomState with fixed seeds, this can
            // be made deterministic (see below), but it is unclear if this will increase
            // bias of resolving components in favor of certain transcripts (and therefore genes)
            // that tend to appear first in hash iteration order.
            // let mut uncovered_vertices = comp_verts.iter().cloned().collect::<HashSet<u32, ahash::RandomState>>();

            // Deterministic variant : replacing the above line with these two lines will
            // cause the parsimony resolution to be deterministic, but potentially at the
            // cost of increasing bias.
            let mut uncovered_vertices = get_set(comp_verts.len() as u32);
            for v in comp_verts.iter().cloned() {
                uncovered_vertices.insert(v);
            }

            // we will remove covered vertices from uncovered_vertices until they are
            // all gone (until all vertices have been covered)
            while !uncovered_vertices.is_empty() {
                let num_remaining = uncovered_vertices.len();
                // will hold vertices in the best mcc
                let mut best_mcc: Vec<u32> = Vec::new();
                // the transcript that is responsible for the
                // best mcc covering
                let mut best_covering_txp = u32::MAX;
                // for each vertex in the vertex set
                for v in uncovered_vertices.iter() {
                    // find the largest mcc starting from this vertex
                    // and the transcript that covers it
                    // NOTE: what if there are multiple different mccs that
                    // are equally good? (@k3yavi — I don't think this case
                    // is even handled in the C++ code either).
                    let (new_mcc, covering_txp) =
                        collapse_vertices(*v, &uncovered_vertices, g, eqmap, hasher_state);

                    let mcc_len = new_mcc.len();
                    // if the new mcc is better than the current best, then
                    // it becomes the new best
                    if best_mcc.len() < mcc_len {
                        best_mcc = new_mcc;
                        best_covering_txp = covering_txp;
                    }
                    // we can't do better than covering all
                    // remaining uncovered vertices.  So, if we
                    // accomplish that, then quit here.
                    if mcc_len == num_remaining {
                        break;
                    }
                }

                if best_covering_txp == u32::MAX {
                    crit!(log, "Could not find a covering transcript");
                    std::process::exit(1);
                }

                // get gene_id of best covering transcript
                let best_covering_gene = if gene_level_eq_map {
                    best_covering_txp
                } else {
                    tid_to_gid[best_covering_txp as usize]
                };

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

                // project each covering transcript to its
                // corresponding gene, and dedup the list
                let mut global_genes: Vec<u32> = if gene_level_eq_map {
                    global_txps.iter().cloned().collect()
                } else {
                    global_txps
                        .iter()
                        .cloned()
                        .map(|i| tid_to_gid[i as usize])
                        .collect()
                };
                // sort since we will be hashing the ordered vector
                global_genes.sort_unstable();
                // dedup as well since we don't care about duplicates
                global_genes.dedup();

                pug_stats.total_mccs += 1;
                if global_genes.len() > 1 {
                    pug_stats.ambiguous_mccs += 1;
                }

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
                // remove it from uncovered_vertices
                for rv in best_mcc.iter() {
                    uncovered_vertices.remove(rv);
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

            let mut global_genes: Vec<u32>;

            if gene_level_eq_map {
                global_genes = tl.to_vec();
            } else {
                global_genes = tl.iter().map(|i| tid_to_gid[*i as usize]).collect();
                global_genes.sort_unstable();
                global_genes.dedup();
            }

            // extract gene-level eqclass and increment count by 1
            assert!(
                !global_genes.is_empty(),
                "can't find representative gene(s) for a molecule"
            );

            pug_stats.total_mccs += 1;
            pug_stats.trivial_mccs += 1;
            if global_genes.len() > 1 {
                pug_stats.ambiguous_mccs += 1;
            }
            // incrementing the count of the eqclass label by 1
            let counter = gene_eqclass_hash.entry(global_genes).or_insert(0);
            *counter += 1;
        }

        //let rand_cover = rand::thread_rng().choose(&tl)
        //    .expect("can;t get random cover");
        //identified_txps.push(*rand_cover as u32);
    }

    /*(gene_eqclass_hash,*/
    pug_stats
    //)
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

pub fn get_num_molecules_forseti(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    tid_to_gid: &[u32],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
    hasher_state: &ahash::RandomState,
    large_graph_thresh: usize,
    cell_reads: &[AlevinFryReadRecordWithPosition],
    read_length: u16,
    max_frag_len: u16,
    ref_names: &[String],
    spliceu_txome: &HashMap<u32, Vec<u8>>,
    spline_lookup: &Array1<f64>,
    mlp: &nn::Sequential,
    tx_status_lookup: &Vec<u8>,
    log: &slog::Logger,
) -> PugResolutionStatistics{
    type U32Set = HashSet<u32, ahash::RandomState>;
    let get_set = |cap: u32| {
        //let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
        U32Set::with_capacity_and_hasher(cap as usize, hasher_state.clone())
    };
    let seed: u64 = 42;
    let gene_level_eq_map = match eqmap.map_type {
        EqMapType::GeneLevel => true,
        EqMapType::TranscriptLevel => false,
    };

    let comps = weakly_connected_components(g);
    // a vector of length 2 that records at index 0
    // the number of single-node subgraphs that are
    // transcript-unique and at index 1 the number of
    // single-node subgraphs that have more than one
    // associated transcript.
    let mut one_vertex_components: Vec<usize> = vec![0, 0];
    
    // Skip some previous commented code in this function.
    
    // Get the genes that could potentially explain all
    // of the vertices in this mcc.
    // To do this, we first extract the set of _transcripts_
    // that label all vertices of the mcc, and then we project
    // the transcripts to their corresponding gene ids.

    //let mut global_txps : Vec<u32>;
    // let mut global_txps = get_set(16);

    let mut pug_stats = PugResolutionStatistics {
        used_alternative_strategy: false,
        total_mccs: 0,
        ambiguous_mccs: 0,
        trivial_mccs: 0,
        non_triv_mccs_determined_ss: (0, 0),
        non_triv_single_gene_ambiguous_ss: 0,
        large_mcc_forseti_needed: 0,
        triv_mcc_check: 0,
        triv_mccs_determined_ss: (0, 0),
        triv_single_gene_ambiguous_ss: 0,
        triv_mcc_forseti_needed: 0,
    };

    for (_comp_label, comp_verts) in comps.iter() {
        let mut best_mcc: Vec<u32> = Vec::with_capacity(comp_verts.len());
        let mut predict_best_cvr_txp_list: Vec<u32> = Vec::with_capacity(16);
        let mut global_txps = Vec::with_capacity(16);
        let mut global_genes = Vec::with_capacity(16);
        let mut all_global_genes = Vec::with_capacity(16);
        let mut full_mcc_txp_pairs: Vec<(Vec<u32>, Vec<u32>)> = Vec::with_capacity(32);

        if comp_verts.len() > 1 {
            // if connected component is large, use existing function. the eq_map will safely use eqc_info rather than eqc_info_forseti.
            if comp_verts.len() > large_graph_thresh {
                get_num_molecules_large_component(
                    g,
                    eqmap,
                    comp_verts,
                    tid_to_gid,
                    hasher_state,
                    gene_eqclass_hash,
                    log,
                );
                warn!(
                    log,
                    "found connected component with {} vertices; resolved with cr-like resolution.",
                    comp_verts.len(),
                );
                pug_stats.used_alternative_strategy = true;
                continue;
            }

            let mut uncovered_vertices = get_set(comp_verts.len() as u32);
            for v in comp_verts.iter().cloned() {
                uncovered_vertices.insert(v);
            }

            // initial Vec buffers
            let mut best_mcc_list: Vec<Vec<u32>> = Vec::with_capacity(comp_verts.len());
            let mut predicted_best_mcc_list: Vec<Vec<u32>> = Vec::with_capacity(comp_verts.len());
            let mut current_global_txps = get_set(16);
            let mut all_global_txps = Vec::with_capacity(48);
            let mut full_mcc_txp_map: HashMap<Vec<u32>, Vec<u32>, ahash::RandomState> = HashMap::with_capacity_and_hasher(comp_verts.len() as usize, hasher_state.clone());
            let mut seen_mccs = HashSet::with_capacity_and_hasher(uncovered_vertices.len(), hasher_state.clone());
            
            // we will remove covered vertices from uncovered_vertices until they are
            // all gone (until all vertices have been covered)
            while !uncovered_vertices.is_empty() {
                let num_remaining = uncovered_vertices.len();
                let mut best_len: usize = 0;
                let mut current_covering_txp = u32::MAX;
                best_mcc.clear();
                global_txps.clear();
                best_mcc_list.clear();
                predicted_best_mcc_list.clear();
                seen_mccs.clear();
                full_mcc_txp_pairs.clear();
                full_mcc_txp_map.clear();
                // the transcript that is responsible for the
                // best mcc covering
                let mut find_single_best = false;

                // for each vertex in the vertex set
                for v in uncovered_vertices.iter() {
                    // find the largest mcc starting from this vertex
                    // and the transcript that covers it
                    // NOTE: mcc is the mcc consists of a set of vertices.(new_mcc.len() = size of the tree )
                    let keep_tie_mcc_txp_map = collapse_vertices_keep_ties(*v, &uncovered_vertices, g, eqmap, hasher_state);

                    let current_max_len = keep_tie_mcc_txp_map.keys().next().unwrap().len();
                    // if the new mcc is larger than the current best, then
                    // it becomes the new best; if the same length, append new ones
                    if current_max_len >= best_len && current_max_len > 0 {
                        if current_max_len > best_len {
                            // Found a new global maximum: reset everything
                            best_len = current_max_len;
                            best_mcc_list.clear();
                            full_mcc_txp_map.clear();
                            seen_mccs.clear();
                        }
                        // If current_max_len == best_len (or we just reset), merge these results
                        for (mcc, txps) in keep_tie_mcc_txp_map.into_iter() {
                            // Check if this specific set of vertices was already found via another starting vertex v
                            if seen_mccs.insert(mcc.clone()) {
                                full_mcc_txp_map.insert(mcc, txps);
                            }
                        }
                    }
                    if best_len == num_remaining { break; }
                }
                // Select initial representative for asserts/logic from the first winning MCC
                let (first_mcc, first_txp_list) = full_mcc_txp_map.iter().next().expect("No MCC found");
                let current_covering_txp = *first_txp_list.first().unwrap();
                // sanity check
                // only check the current_covering_txp, as a representative of the txps
                if current_covering_txp == u32::MAX {
                    crit!(log, "Could not find a covering transcript");
                    std::process::exit(1);
                }
                // get gene_id of best covering transcript
                let current_covering_gene = if gene_level_eq_map {
                    current_covering_txp
                } else {
                    tid_to_gid[current_covering_txp as usize]
                };

                // get all compatible txps for the best mcc list
                all_global_txps.clear();
                // We'll store the (mcc, txp) pairs here for Forseti scoring later
                // Reusing best_mcc_list to hold the final pairs

                for (mcc, txps) in full_mcc_txp_map.iter() {
                    all_global_txps.extend(txps.iter().copied());
                }
                all_global_txps.sort_unstable();
                all_global_txps.dedup();
                // project each covering transcript to its
                // corresponding gene, and dedup the list
                all_global_genes.clear();
                all_global_genes= if gene_level_eq_map {
                    all_global_txps.iter().cloned().collect()
                } else {
                    all_global_txps
                        .iter()
                        .cloned()
                        .map(|i| tid_to_gid[i as usize])
                        .collect()
                };
                // sort since we will be hashing the ordered vector
                all_global_genes.sort_unstable();
                // dedup as well since we don't care about duplicates
                all_global_genes.dedup();

                // assert the best covering gene in the global gene set
                assert!(
                    all_global_genes.contains(&current_covering_gene),
                    "best gene {} not in covering set, shouldn't be possible",
                    current_covering_gene
                );
                assert!(
                    !all_global_genes.is_empty(),
                    "can't find representative gene(s) for a molecule"
                );
                // check if the best_mcc_list (have 1 or more mcc) have multiple gene label/ amb splicing status(s.s.). If so, we will call forseti. Otherwise, we will randomly choose one mcc.

                // call forseti if we have multiple best mccs to determine the final best mcc and best covering txps
                let mut need_forseti = false;

                if all_global_txps.len() == 1 {
                    pug_stats.non_triv_mccs_determined_ss.0 += 1;
                } else if all_global_genes.len() == 1
                    && has_single_strict_suffix_type(&all_global_txps, tx_status_lookup)
                {
                    pug_stats.non_triv_mccs_determined_ss.1 += 1;
                } else {
                    need_forseti = true;
                }
                // as a default, use the first mcc and its global txps.
                best_mcc = full_mcc_txp_map.keys().next().unwrap().clone();
                global_txps = full_mcc_txp_map.values().next().unwrap().clone();
                global_genes = if gene_level_eq_map {
                    global_txps.iter().cloned().collect()
                } else {
                    global_txps
                        .iter()
                        .cloned()
                        .map(|i| tid_to_gid[i as usize])
                        .collect()
                };
                global_genes.sort_unstable();
                global_genes.dedup();

                let mut filtered_mcc_txp_pairs: HashMap<Vec<u32>, Vec<u32>, ahash::RandomState> = HashMap::with_capacity_and_hasher(full_mcc_txp_map.keys().len(), hasher_state.clone());

                // update best_mcc and global_txps if we have useful information from forseti.
                if need_forseti {
                    pug_stats.large_mcc_forseti_needed += 1;
                    // full_mcc_txp_map: mcc -->covering_txp_list
                    // check_mcc_list: [(mcc, covering_ref_id, rname_list1),(mcc, covering_ref_id2, rname_list2),...]
                    let check_mcc_list: Vec<(Vec<u32>, u32, Vec<u64>)> =
                    get_check_mcc_list_from_map(&full_mcc_txp_map, eqmap, g);
                    // get the algn tuple from align file;
                    // forseti_checking_list: {(iter_idx, covering_txp_id): [(dir, pos), (dir2, pos2), ...]}
                    let forseti_checking_list = get_forseti_check_list(cell_reads, &check_mcc_list);
                    // forseti_result returns the iter_idx and max score, allow multiple best score(to a list)

                    let forseti_result = forseti_for_multi_best(
                        &forseti_checking_list,
                        ref_names,
                        spliceu_txome,
                        spline_lookup,
                        mlp,
                        read_length as u16,
                        max_frag_len,
                    );

                    match forseti_result {
                        Ok(best_mcc_indices) => {
                            find_single_best = best_mcc_indices.len() == 1;

                            for mcc_idx in &best_mcc_indices {
                                let mcc_info_tuple = &check_mcc_list[*mcc_idx as usize];
                                let predict_best_cvr_txp = mcc_info_tuple.1;
                                let predict_best_mcc = mcc_info_tuple.0.clone();
                                filtered_mcc_txp_pairs.entry(predict_best_mcc)
                                    .or_insert(Vec::new())
                                    .push(predict_best_cvr_txp);
                            }
                        }
                        Err(e) => {
                            eprintln!(
                                "Error occurred when process the result from forseti: {:?}",
                                e
                            );
                            std::process::exit(1);
                        }
                    }
                    // after forseti
                    global_txps.clear(); //reset global_txps
                    if !filtered_mcc_txp_pairs.is_empty() {
                        // we have narrowed down mcc/txps from forseti. update best_mc and global_txps,  

                        //Even though forseti output all equally best mcc and txp. we will eventually choose just one "best_mcc" to proceed. here we use the first mcc.
                        let (bm, txps) = filtered_mcc_txp_pairs.iter().next().unwrap();
                        best_mcc = bm.clone();
                        global_txps.extend(txps.iter().copied());
                        //instead of uisng its original global_txp, we use its winner txps(narrowed down txps) from forseti.
                    }else{
                        best_mcc = full_mcc_txp_map.keys().next().unwrap().clone();
                        global_txps = full_mcc_txp_map.values().next().unwrap().clone();
                    }

                    // project each covering transcript to its
                    // corresponding gene, and dedup the list
                    global_genes.clear();
                    global_genes = if gene_level_eq_map {
                        global_txps.iter().cloned().collect()
                    } else {
                        global_txps
                            .iter()
                            .cloned()
                            .map(|i| tid_to_gid[i as usize])
                            .collect()
                    };
                    // sort since we will be hashing the ordered vector
                    global_genes.sort_unstable();
                    // dedup as well since we don't care about duplicates
                    global_genes.dedup();

                    if global_txps.len() == 1 {
                        // only 1 tx found, must have determined splicing status and unique gene
                        pug_stats.non_triv_mccs_determined_ss.0 += 1;
                    } else if global_genes.len() == 1 {
                        if has_single_strict_suffix_type(&global_txps, tx_status_lookup) {
                            pug_stats.non_triv_mccs_determined_ss.1 += 1;
                        } else {
                            pug_stats.non_triv_single_gene_ambiguous_ss += 1;
                        }
                    } else {
                        // multi global_genes {
                        pug_stats.ambiguous_mccs += 1;
                    }
                }
                pug_stats.total_mccs += 1;
                // in our hash, increment the count of this equivalence class
                // by 1 (and insert it if we've not seen it yet).
                // NOTE: we reuse the global_genes vector to same memo allocation but here we move it; so we need the std::mem::take to make the global_genes empty.
                assert!(
                    !global_genes.is_empty(),
                    " in the non-trivial mcc, can't find representative gene(s) for a molecule"
                );
                let counter = gene_eqclass_hash.entry(std::mem::take(&mut global_genes)).or_insert(0);
                *counter += 1;

                // for every vertext that has been covered
                // remove it from uncovered_vertices
                for rv in best_mcc.iter() {
                    uncovered_vertices.remove(rv);
                }
            } //end-while
        } else {
            pug_stats.triv_mcc_check += 1;
            pug_stats.total_mccs += 1;
            pug_stats.trivial_mccs += 1;

            // this was a single-vertex subgraph
            let tv = comp_verts.first().expect("can't extract first vertex");
            let tl = eqmap.refs_for_eqc(g.from_index(*tv as usize).0);
            let mut find_single_best: bool = false;
            global_genes.clear();
            // let best_mcc_winner: u32 = *tv;

            if gene_level_eq_map {
                global_genes = tl.to_vec();
            } else {
                global_genes = tl.iter().map(|i| tid_to_gid[*i as usize]).collect();
                global_genes.sort_unstable();
                global_genes.dedup();
            }

            // since we have only 1 mcc. only use forseti to narrow down the better coverting txp list.
            let mut need_forseti = false;
            // check if this mcc has multi-gene labels or has undetermined splicing status )
            if global_genes.len() == 1 {
                if tl.len() == 1 {
                    pug_stats.triv_mccs_determined_ss.0 += 1;
                    one_vertex_components[0] += 1;
                } else if has_single_strict_suffix_type(tl, tx_status_lookup) {
                    one_vertex_components[1] += 1;
                    pug_stats.triv_mccs_determined_ss.1 += 1;
                } else {
                    one_vertex_components[1] += 1;
                    need_forseti = true;
                }
            } else {
                one_vertex_components[1] += 1;
                // multi-genes
                need_forseti = true;
            }
            if need_forseti {
                pug_stats.triv_mcc_forseti_needed += 1;
                // need to call forseti to determine the best txps for splicing status or best gene, and we need to update the global_genes accordingly.

                // tv is the only vertex, tl is the txps of the eqclass
                full_mcc_txp_pairs.clear();
                full_mcc_txp_pairs.push((vec![*tv], tl.to_vec()));

                let check_mcc_list: Vec<(Vec<u32>, u32, Vec<u64>)> =
                    get_check_mcc_list(&full_mcc_txp_pairs, eqmap, g);

                let forseti_checking_list = get_forseti_check_list(cell_reads, &check_mcc_list);

                let forseti_result = forseti_for_multi_best(
                    &forseti_checking_list,
                    ref_names,
                    spliceu_txome,
                    spline_lookup,
                    mlp,
                    read_length as u16,
                    max_frag_len,
                );
                
                predict_best_cvr_txp_list.clear();
                match forseti_result {
                    Ok(best_mcc_indices) => {
                        // assign single best = 1(mcc,txp)
                        find_single_best = best_mcc_indices.len() == 1;

                        for mcc_idx in &best_mcc_indices {
                            let mcc_info_tuple = &check_mcc_list[*mcc_idx as usize];
                            predict_best_cvr_txp_list.push(mcc_info_tuple.1);
                        }
                    }
                    Err(e) => {
                        eprintln!(
                            "Error occurred when process the result from forseti: {:?}",
                            e
                        );
                        std::process::exit(1);
                    }
                }
                if predict_best_cvr_txp_list.len() > 0 {
                    // update the global genes only if the forseti has valid output.
                    global_genes.clear();
                    global_genes = predict_best_cvr_txp_list
                        .iter()
                        .map(|i| tid_to_gid[*i as usize])
                        .collect();
                    global_genes.sort_unstable();
                    global_genes.dedup();

                    //TODO: if single best txp is found, suggest that forseti solve the multi-gene trivial mcc
                    if find_single_best {
                        // single best txp after forseti. have determined splicing status and unique gene
                        pug_stats.triv_mccs_determined_ss.0 += 1;
                    } else if global_genes.len() == 1 {
                        if has_single_strict_suffix_type(&predict_best_cvr_txp_list, tx_status_lookup) {
                            pug_stats.triv_mccs_determined_ss.1 += 1;
                        } else {
                            pug_stats.triv_single_gene_ambiguous_ss += 1;
                        }
                    } else if global_genes.len() > 1 {
                        pug_stats.ambiguous_mccs += 1;
                    }
                }
            }
            assert!(
                !global_genes.is_empty(),
                " in the trivial mcc, can't find representative gene(s) for a molecule"
            );
            let counter = gene_eqclass_hash.entry(std::mem::take(&mut global_genes)).or_insert(0);
            *counter += 1;
        }
    }
    pug_stats
}

// check_mcc_list: [(mcc, covering_ref_id, read_idx_list1), (mcc, covering_ref_id2, read_idx_list2), ...]
fn get_check_mcc_list(
    full_mcc_txp_pairs: &[(Vec<u32>, Vec<u32>)],
    eqmap: &EqMap,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
) -> Vec<(Vec<u32>, u32, Vec<u64>)> {
    // Start from the full_mcc_txp_pairs, get the check_mcc_list(containing covering txp and read_idx_list) required by forseti
    // full_mcc_txp_pairs: [(mcc, covering_txp_list)], where mcc is a vector of vertex indices.
    let mut check_mcc_list: Vec<(Vec<u32>, u32, Vec<u64>)> = Vec::with_capacity(full_mcc_txp_pairs.len()*5);
    for (mcc, covering_txp_list) in full_mcc_txp_pairs.iter() {
        for covering_txp_id in covering_txp_list.iter() {
            let mut check_read_idx_list: Vec<u64> = Vec::new();
            for v in mcc.iter() {
                let (eq_id, umi_idx) = g.from_index(*v as usize);
                let eq = &eqmap.eqc_info_forseti[eq_id as usize];
                let (_umi_bit, read_idx_list) = &eq.umis[umi_idx as usize];
                check_read_idx_list.reserve(read_idx_list.len());
                check_read_idx_list.extend(read_idx_list.iter().copied());
            }
            check_mcc_list.push((mcc.clone(), *covering_txp_id, check_read_idx_list));
        }
    }
    return check_mcc_list;
}


fn get_check_mcc_list_from_map(
    full_mcc_txp_map: &HashMap<Vec<u32>, Vec<u32>, ahash::RandomState>,
    eqmap: &EqMap,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
) -> Vec<(Vec<u32>, u32, Vec<u64>)> {
    // Capacity: Total entries will be the sum of the lengths of all txp lists
    let total_entries: usize = full_mcc_txp_map.values().map(|v| v.len()).sum();
    let mut check_mcc_list: Vec<(Vec<u32>, u32, Vec<u64>)> = Vec::with_capacity(total_entries);
    let mut check_read_idx_list: Vec<u64> = Vec::new();
    for (mcc, covering_txp_list) in full_mcc_txp_map.iter() {
        check_read_idx_list.clear();
        // 1. Pre-collect all read indices for this MCC once
        // (Since all txps in the list cover the same MCC, the reads are the same)
        for v in mcc.iter() {
            let (eq_id, umi_idx) = g.from_index(*v as usize);
            let eq = &eqmap.eqc_info_forseti[eq_id as usize];
            let (_umi_bit, read_idx_list) = &eq.umis[umi_idx as usize];
            check_read_idx_list.extend(read_idx_list.iter().copied());
        }

        // 2. Sort and Dedup the read list to keep it clean and efficient for downstream
        check_read_idx_list.sort_unstable();
        check_read_idx_list.dedup();

        // 3. For every transcript that generated this MCC, push a record
        for &covering_txp_id in covering_txp_list.iter() {
            // We clone the MCC and the Read List for each transcript
            check_mcc_list.push((
                mcc.clone(), 
                covering_txp_id, 
                check_read_idx_list.clone()
            ));
        }
    }
    
    check_mcc_list
}

fn get_forseti_check_list(
    cell_reads: &[AlevinFryReadRecordWithPosition],
    check_mcc_list: &[(Vec<u32>, u32, Vec<u64>)],
) -> HashMap<(usize, u32), Vec<AlgnTuple>> {
    // Forseti only needs alignment tuples for a small subset of (read_idx, covering_txp_id)
    // pairs. Avoid building a large HashMap over all alignments for the cell; instead, look up
        // the tuple on-demand by scanning the corresponding read's alignments (na is typically small).
    #[inline]
    fn find_algn_tuple_for_read(
        cell_reads: &[AlevinFryReadRecordWithPosition],
        read_idx: u64,
        covering_txp_id: u32,
    ) -> Option<AlgnTuple> {
        let r = cell_reads
            .get(read_idx as usize)
            .expect("BUG: read_idx out of bounds when indexing cell_reads");
        for (aln_i, &tid) in r.refs().iter().enumerate() {
            if tid == covering_txp_id {
                let dir = *r.dirs.get(aln_i).expect("BUG: dirs index out of bounds");
                let pos = *r.pos.get(aln_i).expect("BUG: pos index out of bounds");
                return Some((dir, pos));
            }
        }
        None
    }

    // Pre-allocate roughly one entry per mcc
    let mut forseti_checking_list: HashMap<(usize, u32), Vec<AlgnTuple>> =
        HashMap::with_capacity(check_mcc_list.len());

    for (mcc_idx, (mcc, covering_txp_id, read_idx_list)) in check_mcc_list.iter().enumerate() {
        let key = (mcc_idx, *covering_txp_id);
        let vec_ref = forseti_checking_list
            .entry(key)
            .or_insert_with(|| Vec::with_capacity(read_idx_list.len()));
        for &read_idx in read_idx_list {
            if let Some(algn_tuple) = find_algn_tuple_for_read(cell_reads, read_idx, *covering_txp_id)
            {
                vec_ref.push(algn_tuple);
    }
        }
    }
    forseti_checking_list
}

/// Determine the suffix category of a transcript name.
/// Only "-T", "-S", "-U", or no suffix are allowed.
/// Any other trailing "-X" pattern will be treated as an error.
#[inline]
fn strict_suffix_type(name: &str) -> SuffixType {
    match name.rsplit_once('-') {
        None => SuffixType::NoneDash, // no dash at all
        Some((_, tag)) => match tag {
            "T" => SuffixType::T,
            "S" => SuffixType::S,
            "U" => SuffixType::U,
            _ => panic!("Unsupported transcript suffix: -{}", tag),
        },
    }
}
// if have single strict suffix type, return true, otherwise false.
#[inline]
fn has_single_strict_suffix_type(txp_ids: &[u32], status_lookup: &[u8]) -> bool {
    if txp_ids.is_empty() { return false; }

    let first_status = status_lookup[txp_ids[0] as usize];
    
    if first_status == 0 { return false; }

    // check if all the txps have the same status
    txp_ids[1..].iter().all(|&tid| status_lookup[tid as usize] == first_status)
}