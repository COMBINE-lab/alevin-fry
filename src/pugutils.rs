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
use libradicl::record::{
    CollatableMappedRecord, ConvertiblePrimitiveInteger, KnownSize, MappedRecord, RecordContext,
    UmiTaggedRecord,
};

use slog::{crit, info, warn};

use crate::eq_class::{EqMap, EqMapType};
use crate::quant::SplicedAmbiguityModel;
use crate::utils::{self as afutils, EqClassPayload};

type CcMap = HashMap<u32, Vec<u32>, ahash::RandomState>;

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
            if x.0 == y.0 { 0 } else { usize::MAX }
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
                let n_labels = eqmap.refs_for_eqc(nv.0);
                if let Ok(_n) = n_labels.binary_search(txp) {
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

/// Find the largest monochromatic spanning arboresence
/// in the graph `g` starting at vertex `v`.  The arboresence
/// is monochromatic if every vertex can be "covered" by a single
/// transcript (i.e. there exists a transcript that appears in the
/// equivalence class labels of all vertices in the arboresence).
fn collapse_vertices_weighted(
    v: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>, // the set of vertices already covered
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
) -> (Vec<u32>, u32, f64, Vec<(u32, f64)>) {
    // get a new set to hold vertices
    type VertexSet = HashSet<u32, ahash::RandomState>;
    let get_set =
        |cap: u32| VertexSet::with_capacity_and_hasher(cap as usize, hasher_state.clone());

    // will hold the nodes in the largest arboresence found
    let mut highest_prob_mcc: Vec<u32> = Vec::new();
    let mut highest_prob: f64 = 0.0;
    let mut eq_txps_prob: Vec<(u32, f64)> = Vec::new();
    let mut chosen_txp = 0u32;
    let vert = g.from_index(v as usize);

    //unsafe {

    let nvert = g.node_count();

    // for every transcript in the equivalence class
    for (tx_index, txp) in eqmap.refs_for_eqc(vert.0).iter().enumerate() {
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
        let mut current_prob = Vec::new();

        //obtain the average probabilities for this UMI
        let prob_vec = eqmap
            .probs_for_eq_umi_tx(vert.0, vert.1, tx_index)
            .expect("eq and umi should be valid");
        let avg_prob = prob_vec.iter().sum::<f64>() / prob_vec.len() as f64;
        current_prob.push(avg_prob);

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
                let n_labels = eqmap.refs_for_eqc(nv.0);
                if let Ok(_n) = n_labels.binary_search(txp) {
                    bfs_list.push_back(n);

                    //obtain the average probabilities for this UMI
                    let prob_vec = eqmap
                        .probs_for_eq_umi_tx(nv.0, nv.1, tx_index)
                        .expect("eq and umi should be valid");
                    let avg_prob = prob_vec.iter().sum::<f64>() / prob_vec.len() as f64;
                    current_prob.push(avg_prob);
                }
            }
        }

        //compute the average probabilities of the current prob
        let average_current_prob = current_prob.iter().sum::<f64>() / current_prob.len() as f64;
        // if this arboresence is the largest we've yet
        // seen, then record it
        if highest_prob < average_current_prob {
            highest_prob_mcc = current_mcc;
            chosen_txp = *txp;
            highest_prob = average_current_prob;
        }

        eq_txps_prob.push((*txp, average_current_prob));
    }
    //}// unsafe

    (highest_prob_mcc, chosen_txp, highest_prob, eq_txps_prob)
}

#[inline]
fn resolve_num_molecules_crlike_from_vec_prefer_ambig<P: EqClassPayload>(
    umi_gene_count_vec: &mut [(u64, u32, u32)],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
) {
    // A cell whose every (UMI, gene) key was dropped (all low-support) has nothing to resolve.
    if umi_gene_count_vec.is_empty() {
        return;
    }
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

            gene_eqclass_hash
                .entry(best_genes.clone())
                .or_insert(P::new(best_genes.len()))
                .inc();

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
            gene_eqclass_hash
                .entry(best_genes.clone())
                .or_insert(P::new(best_genes.len()))
                .inc();
        }
    }
}

#[inline]
fn resolve_num_molecules_crlike_from_vec<P: EqClassPayload>(
    umi_gene_count_vec: &mut [(u64, u32, u32)],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
) {
    // A cell whose every (UMI, gene) key was dropped (all low-support) has nothing to resolve.
    if umi_gene_count_vec.is_empty() {
        return;
    }
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
            gene_eqclass_hash
                .entry(best_genes.clone())
                .or_insert(P::new(best_genes.len()))
                .inc();

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
            gene_eqclass_hash
                .entry(best_genes.clone())
                .or_insert(P::new(best_genes.len()))
                .inc();
        }
    }
}

/// Reused working buffers for [`correct_umis_cellranger`], held per worker (on
/// [`CrLikeScratch`]) so the correction allocates nothing across cells. Every
/// buffer is `clear()`ed at the start of each call, keeping its capacity.
#[derive(Default)]
pub struct CorrScratch {
    /// `(umi, gene, count)` triplets after dedup-summing, sorted by `(umi, gene)`.
    raw: Vec<(u64, u32, u32)>,
    /// `(umi, gene) -> index into `raw``. Built lazily — only when some gene
    /// group is large enough to use the neighbour-enumeration probe.
    key_index: HashMap<(u64, u32), u32, RandomState>,
    /// gene id -> dense bucket id, for the O(n) counting sort that groups
    /// `idx_by_gene` by gene.
    gene_slot: HashMap<u32, u32, RandomState>,
    /// distinct gene ids in first-seen order (bucket id -> gene id).
    genes: Vec<u32>,
    /// per-bucket group size (number of distinct UMIs in that gene).
    counts: Vec<u32>,
    /// per-bucket running write offset during the counting-sort scatter.
    cursor: Vec<u32>,
    /// indices into `raw`, grouped by gene (ascending UMI within a gene, since
    /// `raw` is sorted by `(umi, gene)`).
    idx_by_gene: Vec<u32>,
    /// `dest[i]` = index into `raw` of `raw[i]`'s correction target (self if none).
    dest: Vec<u32>,
    /// intermediate one-read-moved counts (drop path only).
    inter: Vec<u32>,
    /// low-support flag per `raw` index (drop path only).
    low: Vec<bool>,
    /// merged read count per destination `raw` index.
    merged: Vec<u32>,
}

/// Cell Ranger-style Hamming-1 UMI correction within one cell, on `(umi, gene,
/// read count)` triplets (duplicate `(umi, gene)` pairs are summed first).
/// Mirrors cellranger's `mark_dups.rs` (`correct_umis` +
/// `determine_low_support_umigenes`).
///
/// UMIs are 2-bit packed MSB-first with A<C<G<T, so integer order equals
/// lexicographic order. `umi_len` MUST be the *packed* width (the RAD `ulen`
/// tag): the neighbour scan and the tie-break rely on it. Any producer padding of
/// short UMIs is a constant low-bit suffix at equal length, so it never creates a
/// spurious neighbour and never changes the relative order — no masking needed.
///
/// Steps:
/// 1. Correction: each `(umi, gene)` is relabelled to the best of {self} and its
///    Hamming-1 neighbours *of the same gene*, ordered by `(read count,
///    lexicographically larger UMI)`. Decided on raw counts (one step, no
///    chaining).
/// 2. Low-support ("chimeric") determination on an intermediate table where each
///    corrected UMI has moved exactly ONE read to its destination (Cell Ranger 3
///    quirk). Within a UMI, every gene below the per-UMI max count is
///    low-support; if the max is tied, every gene is.
/// 3. Relabel + merge into the destination bucket.
///
/// `drop_low_support` selects the downstream semantics:
/// - `true` (cr-like / winner-take-all): run all three steps and DROP low-support
///   keys, so at most one gene survives per UMI — output-identical to Cell Ranger.
/// - `false` (cr-like-em): run step 1 (correction) and step 3 (relabel + merge)
///   but SKIP the step-2 drop, so a multi-gene UMI survives as several
///   `(dest, gene)` keys and reaches the EM as an equivalence class. Only
///   Hamming-1 UMIs are merged; no molecule is discarded.
///
/// Not meaningful in USA mode (gene ids there are spliced/unspliced variants);
/// callers reject `--umi-edit-dist >= 1` there. Genes with a single distinct UMI
/// in the cell are skipped (no same-gene neighbour can exist) — the dominant
/// saving on real data.
// The all-pairs step-1 loops index `idx_by_gene` by position (needed for the
// self-skip and to map back into `raw`), so the range-loop form is intentional.
#[allow(clippy::needless_range_loop)]
pub fn correct_umis_cellranger(
    v: &mut Vec<(u64, u32, u32)>,
    umi_len: u32,
    drop_low_support: bool,
    cs: &mut CorrScratch,
) {
    if v.len() < 2 {
        return;
    }
    debug_assert!(
        umi_len <= 32,
        "umi_len ({umi_len}) exceeds the 32-base capacity of the u64 UMI packing; the neighbour-flip shift x << (2*pos) would overflow"
    );
    // Below this per-gene group size, all-pairs Hamming via 2-bit popcount beats
    // probing the 3*umi_len bit-flip neighbourhood against the key index.
    const SMALL_GENE_GROUP: usize = 64;

    let CorrScratch {
        raw,
        key_index,
        gene_slot,
        genes,
        counts,
        cursor,
        idx_by_gene,
        dest,
        inter,
        low,
        merged,
    } = cs;

    // Sort by (umi, gene). For the common case (umi_len <= 16, i.e. umi fits in
    // 32 bits) pack the key into a u64 and radix-sort — ~1.2-1.35x faster than
    // pdqsort on realistic per-cell sizes; the count field is payload and is
    // irrelevant post-dedup. Above 16 bases the packed key would overflow u64, so
    // fall back to the comparison sort (unreachable for real 10x/Flex UMIs).
    if umi_len <= 16 {
        radsort::sort_by_key(v, |&(umi, gene, _)| (umi << 32) | gene as u64);
    } else {
        v.sort_unstable();
    }
    // dedup-sum identical (umi, gene) into `raw` (stays sorted by (umi, gene)).
    raw.clear();
    for &t in v.iter() {
        if let Some(l) = raw.last_mut()
            && l.0 == t.0
            && l.1 == t.1
        {
            l.2 += t.2;
            continue;
        }
        raw.push(t);
    }
    let n = raw.len();

    dest.clear();
    dest.extend(0..n as u32);

    // O(n) counting sort that groups `idx_by_gene` by gene. Because `raw` is
    // sorted by (umi, gene), scattering in raw order leaves each gene's entries
    // in ascending-UMI order (identical to a (gene, umi) comparison sort) —
    // step 1 only needs the grouping, and steps 2/3 read `raw` order.
    gene_slot.clear();
    genes.clear();
    counts.clear();
    for &(_, g, _) in raw.iter() {
        let slot = *gene_slot.entry(g).or_insert_with(|| {
            genes.push(g);
            counts.push(0);
            (genes.len() - 1) as u32
        });
        counts[slot as usize] += 1;
    }
    // prefix-sum group sizes into write offsets, and note if any group is large.
    cursor.clear();
    cursor.resize(genes.len(), 0);
    let mut acc = 0u32;
    let mut any_large = false;
    for b in 0..genes.len() {
        cursor[b] = acc;
        acc += counts[b];
        if counts[b] as usize > SMALL_GENE_GROUP {
            any_large = true;
        }
    }
    idx_by_gene.clear();
    idx_by_gene.resize(n, 0);
    for i in 0..n {
        let b = gene_slot[&raw[i].1] as usize;
        idx_by_gene[cursor[b] as usize] = i as u32;
        cursor[b] += 1;
    }

    // `key_index` is only needed for the large-group neighbour-enumeration
    // probe; build it lazily. On Flex/scRNA every gene group is tiny, so this is
    // skipped entirely and the small-group path tracks the winner's index.
    if any_large {
        key_index.clear();
        for (i, &(u, g, _)) in raw.iter().enumerate() {
            key_index.insert((u, g), i as u32);
        }
    }

    // step 1: per-gene correction; groups are contiguous in `idx_by_gene`, with
    // boundaries straight from the prefix offsets (single-UMI genes skipped).
    let mut start = 0usize;
    for b in 0..genes.len() {
        let k = counts[b] as usize;
        let p = start;
        let q = start + k;
        start = q;
        if k < 2 {
            continue;
        }
        let g = genes[b];
        if k <= SMALL_GENE_GROUP {
            // all-pairs Hamming-1 via 2-bit popcount (no hashing); the winner is
            // an in-group entry whose index we already know.
            for a in p..q {
                let ia = idx_by_gene[a] as usize;
                let (ua, _, ca) = raw[ia];
                let mut best = (ca, ua);
                let mut best_idx = ia;
                for bb in p..q {
                    if bb == a {
                        continue;
                    }
                    let ib = idx_by_gene[bb] as usize;
                    let (ub, _, cb) = raw[ib];
                    if afutils::count_diff_2_bit_packed(ua, ub) == 1 && (cb, ub) > best {
                        best = (cb, ub);
                        best_idx = ib;
                    }
                }
                if best_idx != ia {
                    dest[ia] = best_idx as u32;
                }
            }
        } else {
            // probe the 3*umi_len bit-flip neighbourhood against the index
            for a in p..q {
                let ia = idx_by_gene[a] as usize;
                let (ua, _, ca) = raw[ia];
                let mut best = (ca, ua);
                let mut best_idx = ia;
                for pos in 0..umi_len {
                    for x in 1u64..=3 {
                        let neigh = ua ^ (x << (2 * pos));
                        if let Some(&j) = key_index.get(&(neigh, g)) {
                            let (uj, _, cj) = raw[j as usize];
                            if (cj, uj) > best {
                                best = (cj, uj);
                                best_idx = j as usize;
                            }
                        }
                    }
                }
                if best_idx != ia {
                    dest[ia] = best_idx as u32;
                }
            }
        }
    }

    // step 2: low-support determination (drop path only).
    low.clear();
    low.resize(n, false);
    if drop_low_support {
        inter.clear();
        inter.extend(raw.iter().map(|t| t.2));
        for i in 0..n {
            let d = dest[i] as usize;
            if d != i {
                inter[i] -= 1;
                inter[d] += 1;
            }
        }
        // UMI runs are contiguous in `raw` (sorted by (umi, gene)).
        let mut i = 0usize;
        while i < n {
            let u = raw[i].0;
            let mut j = i;
            while j < n && raw[j].0 == u {
                j += 1;
            }
            let maxc = (i..j).map(|k| inter[k]).max().unwrap();
            let tied = (i..j).filter(|&k| inter[k] == maxc).count() >= 2;
            for k in i..j {
                if tied || inter[k] < maxc {
                    low[k] = true;
                }
            }
            i = j;
        }
    }

    // step 3: relabel + merge; emit surviving buckets in (umi, gene) order.
    merged.clear();
    merged.resize(n, 0);
    for i in 0..n {
        merged[dest[i] as usize] += raw[i].2;
    }
    v.clear();
    for b in 0..n {
        if merged[b] > 0 && !low[b] {
            v.push((raw[b].0, raw[b].1, merged[b]));
        }
    }
}

/// Per-worker scratch buffers reused by the cr-like resolvers across cells.
///
/// The hot path is one call to a `get_num_molecules_cell_ranger_like*` per cell,
/// over millions of cells. Without reuse each call allocates a fresh working
/// vector (and each equivalence class a fresh gene-set vector), which the
/// profile shows dominating quant's CPU. Holding these buffers on the worker and
/// `clear()`ing them per use keeps the capacity and removes the churn. It is
/// purely an allocation optimisation: the contents built each call are identical
/// to the freshly-allocated version, so the output is unchanged.
#[derive(Default)]
pub struct CrLikeScratch {
    /// `(umi, gene_id, count)` triplets accumulated for one cell.
    umi_gene_count_vec: Vec<(u64, u32, u32)>,
    /// Deduplicated gene ids for one equivalence class / record.
    gset: Vec<u32>,
    /// Reused buffers for the cr-like Hamming-1 UMI correction.
    corr: CorrScratch,
}

#[allow(clippy::too_many_arguments)]
pub fn get_num_molecules_cell_ranger_like_small<B, R, P: EqClassPayload>(
    cell_chunk: &mut chunk::Chunk<R>,
    tid_to_gid: &[u32],
    _num_genes: usize,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
    sa_model: SplicedAmbiguityModel,
    scratch: &mut CrLikeScratch,
    umi_edit: u32,
    umi_len: u32,
    drop_low_support: bool,
    _log: &slog::Logger,
) where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize + UmiTaggedRecord,
    <R as MappedRecord>::ParsingContext: RecordContext,
    <R as MappedRecord>::ParsingContext: Clone,
    <R as MappedRecord>::ParsingContext: Send,
{
    // Disjoint borrows of the reusable buffers; cleared before use so prior
    // capacity is kept without carrying stale contents.
    let CrLikeScratch {
        umi_gene_count_vec,
        gset,
        corr,
    } = scratch;
    umi_gene_count_vec.clear();
    umi_gene_count_vec.reserve(cell_chunk.nrec as usize);

    // for each record
    for rec in &cell_chunk.reads {
        // get the umi
        let umi = rec.umi();

        // project the transcript ids to gene ids
        gset.clear();
        gset.extend(rec.refs().iter().map(|tid| tid_to_gid[*tid as usize]));
        // and make the gene ids unique
        gset.sort_unstable();
        gset.dedup();
        for g in gset.iter() {
            umi_gene_count_vec.push((umi, *g, 1));
        }
    }
    // Cell Ranger-style Hamming-1 UMI correction (per cell, per gene) before
    // resolution, when --umi-edit-dist >= 1. `drop_low_support` keeps the full
    // Cell Ranger chimera drop for cr-like (winner-take-all) and skips it for
    // cr-like-em so multi-gene UMIs survive as eqclasses for the EM.
    if umi_edit >= 1 {
        correct_umis_cellranger(umi_gene_count_vec, umi_len, drop_low_support, corr);
    }
    match sa_model {
        SplicedAmbiguityModel::WinnerTakeAll => {
            resolve_num_molecules_crlike_from_vec(umi_gene_count_vec, gene_eqclass_hash);
        }
        SplicedAmbiguityModel::PreferAmbiguity => {
            resolve_num_molecules_crlike_from_vec_prefer_ambig(
                umi_gene_count_vec,
                gene_eqclass_hash,
            );
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn get_num_molecules_cell_ranger_like<P: EqClassPayload>(
    eq_map: &EqMap,
    tid_to_gid: &[u32],
    _num_genes: usize,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
    sa_model: SplicedAmbiguityModel,
    scratch: &mut CrLikeScratch,
    umi_edit: u32,
    umi_len: u32,
    drop_low_support: bool,
    _log: &slog::Logger,
) {
    // Disjoint borrows of the reusable buffers; cleared before use so prior
    // capacity is kept without carrying stale contents.
    let CrLikeScratch {
        umi_gene_count_vec,
        gset,
        corr,
    } = scratch;
    umi_gene_count_vec.clear();

    // for each equivalence class
    for eqinfo in &eq_map.eqc_info {
        // get the (umi, count) pairs
        let umis = &eqinfo.umis;
        let eqid = &eqinfo.eq_num;

        // project the transcript ids to gene ids
        gset.clear();
        gset.extend(
            eq_map
                .refs_for_eqc(*eqid)
                .iter()
                .map(|tid| tid_to_gid[*tid as usize]),
        );
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
            for g in gset.iter() {
                umi_gene_count_vec.push((umi_ct.0, *g, umi_ct.1));
            }
        }
    }
    // Cell Ranger-style Hamming-1 UMI correction (per cell, per gene) before
    // resolution, when --umi-edit-dist >= 1. `drop_low_support` keeps the full
    // Cell Ranger chimera drop for cr-like (winner-take-all) and skips it for
    // cr-like-em so multi-gene UMIs survive as eqclasses for the EM.
    if umi_edit >= 1 {
        correct_umis_cellranger(umi_gene_count_vec, umi_len, drop_low_support, corr);
    }
    match sa_model {
        SplicedAmbiguityModel::WinnerTakeAll => {
            resolve_num_molecules_crlike_from_vec(umi_gene_count_vec, gene_eqclass_hash);
        }
        SplicedAmbiguityModel::PreferAmbiguity => {
            resolve_num_molecules_crlike_from_vec_prefer_ambig(
                umi_gene_count_vec,
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
fn get_num_molecules_large_component<P: EqClassPayload>(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eq_map: &EqMap,
    vertex_ids: &[u32],
    tid_to_gid: &[u32],
    hasher_state: &ahash::RandomState,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
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
pub fn get_num_molecules<P: EqClassPayload>(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    tid_to_gid: &[u32],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
    hasher_state: &ahash::RandomState,
    large_graph_thresh: usize,
    log: &slog::Logger,
) -> PugResolutionStatistics
//,)
{
    const EMPTY_VEC: Vec<(u32, f64)> = vec![];
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
    };

    for comp_verts in comps.values() {
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
            for v in comp_verts.iter() {
                uncovered_vertices.insert(*v);
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

                // in the long-read case
                let mut best_mcc_prob: f64 = 0.0;
                let mut best_mcc_txp_probs: Vec<(u32, f64)> = Vec::new();

                // for each vertex in the vertex set
                for v in uncovered_vertices.iter() {
                    // find the largest mcc starting from this vertex
                    // and the transcript that covers it
                    // NOTE: what if there are multiple different mccs that
                    // are equally good? (@k3yavi — I don't think this case
                    // is even handled in the C++ code either).
                    let (cand_mcc, cand_txp, cand_prob, eq_txs_prob) = if P::HAS_PROBS {
                        collapse_vertices_weighted(*v, &uncovered_vertices, g, eqmap, hasher_state)
                    } else {
                        let (new_mcc, covering_txp) =
                            collapse_vertices(*v, &uncovered_vertices, g, eqmap, hasher_state);
                        (new_mcc, covering_txp, 0_f64, EMPTY_VEC)
                    };

                    let mcc_len = cand_mcc.len();
                    if P::HAS_PROBS {
                        if best_mcc_prob < cand_prob {
                            best_mcc = cand_mcc;
                            best_mcc_prob = cand_prob;
                            best_covering_txp = cand_txp;
                            best_mcc_txp_probs = eq_txs_prob;
                        }
                    } else {
                        // if the new mcc is better than the current best, then
                        // it becomes the new best
                        if best_mcc.len() < mcc_len {
                            best_mcc = cand_mcc;
                            best_covering_txp = cand_txp;
                        }
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

                // for long reads obtain the probabiltiy for the global txps
                let mut global_txp_prob: Vec<f64> = vec![];
                if P::HAS_PROBS {
                    let mut txp_prob_temp: Vec<(u32, f64)> = best_mcc_txp_probs
                        .iter()
                        .filter(|(t, _)| global_txps.contains(t))
                        .cloned()
                        .collect();
                    txp_prob_temp.sort_unstable_by_key(|(t, _)| *t);
                    global_txp_prob = if txp_prob_temp.len() == 1 {
                        vec![1.0]
                    } else {
                        txp_prob_temp.iter().map(|(_, p)| *p).collect()
                    };
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
                let eq_label_len = global_genes.len();
                let payload = gene_eqclass_hash
                    .entry(global_genes)
                    .or_insert(P::new(eq_label_len));
                payload.inc();

                if P::HAS_PROBS {
                    payload.add_probs(&global_txp_prob);
                }

                // for every vertex that has been covered
                // remove it from uncovered_vertices
                for rv in best_mcc.iter() {
                    uncovered_vertices.remove(rv);
                }
            } //end-while
        } else {
            // this was a single-vertex subgraph
            let tv = comp_verts.first().expect("can't extract first vertex");
            let tl = eqmap.refs_for_eqc(g.from_index(*tv as usize).0);

            let mut global_txp_prob: Vec<f64> = vec![];

            let vcindex: usize = if tl.len() == 1 { 0 } else { 1 };
            one_vertex_components[vcindex] += 1;

            if P::HAS_PROBS {
                if tl.len() == 1 {
                    global_txp_prob = vec![1.0];
                } else {
                    let mut txp_prob_temp: Vec<(u32, f64)> = Vec::with_capacity(tl.len());
                    for (i, t) in tl.iter().enumerate() {
                        let (eq_id, umi_id) = g.from_index(*tv as usize);
                        let prob_vec = eqmap
                            .probs_for_eq_umi_tx(eq_id, umi_id, i)
                            .expect("should be a valid eq/umi pair");
                        let avg_prob = prob_vec.iter().sum::<f64>() / prob_vec.len() as f64;
                        txp_prob_temp.push((*t, avg_prob));
                    }

                    txp_prob_temp.sort_unstable_by_key(|(t, _)| *t);
                    global_txp_prob = txp_prob_temp.iter().map(|(_, p)| *p).collect();
                }
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
            let eq_label_len = global_genes.len();
            let payload = gene_eqclass_hash
                .entry(global_genes)
                .or_insert(P::new(eq_label_len));
            payload.inc();
            if P::HAS_PROBS {
                payload.add_probs(&global_txp_prob);
            }
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

#[cfg(test)]
mod cellranger_umi_tests {
    // Ported from ygao61's branch umi_ham_edit_1 (commit d92869f); these
    // include Cell Ranger's own mark_dups.rs test_correct_umis cases. They
    // exercise the drop path (`drop_low_support = true`, i.e. cr-like).
    use super::{CorrScratch, correct_umis_cellranger};
    // 4-nt UMIs, 2-bit MSB-first, A=0 C=1 G=2 T=3
    const AAAA: u64 = 0b0000_0000;
    const AAAT: u64 = 0b0000_0011;
    const AATT: u64 = 0b0000_1111;
    const CCCC: u64 = 0b0101_0101;
    const CGCC: u64 = 0b0110_0101;

    /// cellranger mark_dups.rs test_correct_umis case 1: AAAT(g0,2) -> AAAA(g0,3); then UMI AAAA
    /// is 4 reads in g0 vs 1 in g1 on the intermediate table, so (AAAA,g1) is low-support.
    #[test]
    fn greater_count_neighbour_and_low_support() {
        let mut v = vec![(AAAA, 0, 3), (AAAT, 0, 2), (AAAA, 1, 1), (AATT, 1, 1)];
        correct_umis_cellranger(&mut v, 32, true, &mut CorrScratch::default());
        assert_eq!(v, vec![(AAAA, 0, 5), (AATT, 1, 1)]);
    }

    /// case 2: equal counts -> merge into the lexicographically larger UMI (CCCC -> CGCC).
    #[test]
    fn equal_count_lexicographic_tiebreak() {
        let mut v = vec![(CCCC, 0, 1), (CGCC, 0, 1)];
        correct_umis_cellranger(&mut v, 32, true, &mut CorrScratch::default());
        assert_eq!(v, vec![(CGCC, 0, 2)]);
    }

    /// no chaining: A(1) -> B(2) and B(2) -> C(3) are decided on raw counts; A's reads land on B,
    /// B's reads on C, so both B and C survive as molecules (Cell Ranger one-step relabel).
    #[test]
    fn one_step_no_chaining() {
        let a = AAAA;
        let b = AAAT;
        let c = AATT; // a~b and b~c are Hamming-1, a~c is Hamming-2
        let mut v = vec![(a, 0, 1), (b, 0, 2), (c, 0, 3)];
        correct_umis_cellranger(&mut v, 32, true, &mut CorrScratch::default());
        assert_eq!(v, vec![(b, 0, 1), (c, 0, 5)]);
    }

    /// correction is per gene: a Hamming-1 neighbour in another gene does not attract reads.
    #[test]
    fn correction_is_within_gene() {
        let mut v = vec![(AAAA, 0, 1), (AAAT, 1, 5)];
        correct_umis_cellranger(&mut v, 32, true, &mut CorrScratch::default());
        assert_eq!(v, vec![(AAAA, 0, 1), (AAAT, 1, 5)]);
    }

    /// Packing sanity: the integer order of MSB-first 2-bit UMIs equals lexicographic order over
    /// ACGT. The tie-break in the correction ("merge into the lexicographically larger UMI")
    /// compares the packed integers directly and relies on this.
    #[test]
    fn packing_order_is_lexicographic() {
        fn pack(s: &str) -> u64 {
            s.bytes().fold(0u64, |k, b| {
                (k << 2)
                    | match b {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        b'T' => 3,
                        _ => panic!("bad base {b}"),
                    }
            })
        }
        let umis = [
            "AAAA", "AAAT", "ACGT", "CAAA", "GTTT", "TAAA", "TTTA", "TTTT", "CGCC", "CCCC",
        ];
        for a in umis {
            for b in umis {
                assert_eq!(pack(a).cmp(&pack(b)), a.cmp(b), "{a} vs {b}");
            }
        }
    }

    /// Every key low-support: the correction returns an empty vector, which the callers in
    /// `resolve_num_molecules_crlike_from_vec*` must tolerate. AAAA and AAAT are Hamming-1 with
    /// one read each in genes 0 and 1, so every UMI ends up with a tied per-gene maximum.
    #[test]
    fn all_low_support_gives_empty_vec() {
        let mut v = vec![(AAAA, 0, 1), (AAAT, 0, 1), (AAAA, 1, 1), (AAAT, 1, 1)];
        correct_umis_cellranger(&mut v, 32, true, &mut CorrScratch::default());
        assert!(v.is_empty());
    }
}

#[cfg(test)]
mod umi_ham_merge_tests {
    //! Merge-specific integration test (not in d92869f): the empty vector that
    //! `correct_umis_cellranger` can produce must flow through the resolver's
    //! empty guard without panicking on `.first().expect(...)`. This is the test
    //! that would catch a misplaced guard (advisor item 1).
    use super::*;
    use crate::utils::BasicEqClassPayload;

    #[test]
    fn all_low_support_result_resolves_without_panic() {
        // Same all-low-support input as ygao's `all_low_support_gives_empty_vec`.
        const AAAA: u64 = 0b0000_0000;
        const AAAT: u64 = 0b0000_0011;
        let mut v = vec![(AAAA, 0u32, 1u32), (AAAT, 0, 1), (AAAA, 1, 1), (AAAT, 1, 1)];
        correct_umis_cellranger(&mut v, 32, true, &mut CorrScratch::default());
        assert!(v.is_empty(), "all keys low-support -> empty");

        let mut h: HashMap<Vec<u32>, BasicEqClassPayload, ahash::RandomState> =
            HashMap::with_hasher(ahash::RandomState::with_seeds(2, 7, 1, 8));
        // Must return via the empty guard, leaving the eqclass hash empty.
        resolve_num_molecules_crlike_from_vec(&mut v, &mut h);
        assert!(h.is_empty(), "empty cell yields no molecules");
    }
}

#[cfg(test)]
mod umi_correction_optimization_tests {
    //! Guards for the gene-aware, allocation-free `correct_umis_cellranger`
    //! rewrite (advisor items C1 + performance): a naive reference of the exact
    //! three-step algorithm, a property test that the optimized version matches
    //! it in both modes, the cr-like-em multi-gene-survival semantics, and the
    //! needletail encoding-order invariant the tie-break relies on.
    use super::{CorrScratch, correct_umis_cellranger};
    use std::collections::{HashMap, HashSet};

    /// Straightforward, obviously-correct reference implementation of the three
    /// steps, parameterized by `drop_low_support`. No optimization; used only to
    /// validate the optimized `correct_umis_cellranger`.
    fn naive_correct(v: &mut Vec<(u64, u32, u32)>, umi_len: u32, drop_low_support: bool) {
        if v.len() < 2 {
            return;
        }
        v.sort_unstable();
        let mut raw: Vec<(u64, u32, u32)> = Vec::new();
        for &t in v.iter() {
            if let Some(l) = raw.last_mut()
                && l.0 == t.0
                && l.1 == t.1
            {
                l.2 += t.2;
                continue;
            }
            raw.push(t);
        }
        let mut counts: HashMap<(u64, u32), u32> = HashMap::new();
        for &(u, g, c) in &raw {
            counts.insert((u, g), c);
        }
        // step 1: corrections on raw counts
        let mut corr: HashMap<(u64, u32), u64> = HashMap::new();
        for &(u, g, c) in &raw {
            let mut best = (c, u);
            for pos in 0..umi_len {
                for x in 1u64..=3 {
                    let n = u ^ (x << (2 * pos));
                    if let Some(&nc) = counts.get(&(n, g))
                        && (nc, n) > best
                    {
                        best = (nc, n);
                    }
                }
            }
            if best.1 != u {
                corr.insert((u, g), best.1);
            }
        }
        // step 2: low support on the one-read-moved table (drop path only)
        let mut low: HashSet<(u64, u32)> = HashSet::new();
        if drop_low_support {
            let mut inter = counts.clone();
            for (&(u, g), &d) in &corr {
                *inter.get_mut(&(u, g)).unwrap() -= 1;
                *inter.get_mut(&(d, g)).unwrap() += 1;
            }
            let mut iv: Vec<(u64, u32, u32)> =
                inter.iter().map(|(&(u, g), &c)| (u, g, c)).collect();
            iv.sort_unstable();
            let mut i = 0;
            while i < iv.len() {
                let mut j = i;
                while j < iv.len() && iv[j].0 == iv[i].0 {
                    j += 1;
                }
                let maxc = iv[i..j].iter().map(|t| t.2).max().unwrap();
                let tied = iv[i..j].iter().filter(|t| t.2 == maxc).count() >= 2;
                for &(u, g, c) in &iv[i..j] {
                    if tied || c < maxc {
                        low.insert((u, g));
                    }
                }
                i = j;
            }
        }
        // step 3: relabel, merge, (drop)
        let mut fin: HashMap<(u64, u32), u32> = HashMap::new();
        for &(u, g, c) in &raw {
            let d = corr.get(&(u, g)).copied().unwrap_or(u);
            *fin.entry((d, g)).or_insert(0) += c;
        }
        let mut out: Vec<(u64, u32, u32)> = fin
            .into_iter()
            .filter(|(k, _)| !low.contains(k))
            .map(|((u, g), c)| (u, g, c))
            .collect();
        out.sort_unstable();
        v.clear();
        v.extend(out);
    }

    /// The optimized `correct_umis_cellranger` must match `naive_correct` for both
    /// `drop_low_support` values, over many random small cells (covers the
    /// all-pairs small-gene-group path) plus a large single-gene case (covers the
    /// neighbour-enumeration path).
    #[test]
    fn optimized_matches_naive_reference() {
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut rng = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut scratch = CorrScratch::default();
        // Two UMI-length regimes: umi_len=4 exercises the radix-sort path
        // (umi_len <= 16, packed u64 key); umi_len=20 exercises the
        // sort_unstable fallback (umi does not fit the packed u64 key).
        for &(umi_len, umi_space) in &[(4u32, 256u64), (20u32, 1u64 << 24)] {
            for _ in 0..3000 {
                let nkeys = (rng() % 38 + 2) as usize;
                let mut input: Vec<(u64, u32, u32)> = Vec::with_capacity(nkeys);
                for _ in 0..nkeys {
                    let umi = rng() % umi_space;
                    let gene = (rng() % 4) as u32; // few genes -> real per-gene groups
                    let cnt = (rng() % 4 + 1) as u32;
                    input.push((umi, gene, cnt));
                }
                for &drop in &[true, false] {
                    let mut opt = input.clone();
                    correct_umis_cellranger(&mut opt, umi_len, drop, &mut scratch);
                    let mut nai = input.clone();
                    naive_correct(&mut nai, umi_len, drop);
                    assert_eq!(
                        opt, nai,
                        "mismatch umi_len={umi_len} drop={drop} input={input:?}"
                    );
                }
            }
        }
        // 100 distinct UMIs in one gene forces the >SMALL_GENE_GROUP path.
        let big: Vec<(u64, u32, u32)> =
            (0u64..100).map(|u| (u, 0u32, (u as u32 % 5) + 1)).collect();
        for &drop in &[true, false] {
            let mut opt = big.clone();
            correct_umis_cellranger(&mut opt, 4, drop, &mut scratch);
            let mut nai = big.clone();
            naive_correct(&mut nai, 4, drop);
            assert_eq!(opt, nai, "large-group mismatch drop={drop}");
        }
    }

    /// cr-like drops a fully-chimeric UMI set to nothing (the winner-take-all
    /// behaviour), but cr-like-em must PRESERVE the corrected multi-gene UMI as
    /// keys in both genes so the EM can resolve it.
    #[test]
    fn crlike_em_preserves_multigene_umi_that_crlike_drops() {
        const AAAA: u64 = 0b0000_0000;
        const AAAT: u64 = 0b0000_0011; // Hamming-1 of AAAA
        let input = vec![(AAAA, 0u32, 1u32), (AAAT, 0, 1), (AAAA, 1, 1), (AAAT, 1, 1)];

        // cr-like (drop): every key is tied/chimeric on the intermediate table -> empty.
        let mut dropv = input.clone();
        correct_umis_cellranger(&mut dropv, 4, true, &mut CorrScratch::default());
        assert!(dropv.is_empty(), "cr-like drops the all-chimeric cell");

        // cr-like-em (no drop): AAAA corrects into AAAT in both genes, and the
        // resulting multi-gene UMI AAAT survives (genes 0 and 1) for the EM.
        let mut keepv = input;
        correct_umis_cellranger(&mut keepv, 4, false, &mut CorrScratch::default());
        assert_eq!(keepv, vec![(AAAT, 0, 2), (AAAT, 1, 2)]);
    }

    /// The tie-break compares packed UMIs as integers and relies on the *actual*
    /// needletail encoding (MSB-first, A<C<G<T) making integer order ==
    /// lexicographic order. This exercises `BitNuclKmer`, not an in-test packer.
    #[test]
    fn needletail_packing_order_is_lexicographic() {
        fn pack(s: &str) -> u64 {
            let (_, bk, _) =
                needletail::bitkmer::BitNuclKmer::new(s.as_bytes(), s.len() as u8, false)
                    .next()
                    .unwrap();
            bk.0
        }
        let umis = [
            "AAAA", "AAAT", "ACGT", "CAAA", "GTTT", "TAAA", "TTTA", "TTTT", "CGCC", "CCCC",
        ];
        for a in umis {
            for b in umis {
                assert_eq!(pack(a).cmp(&pack(b)), a.cmp(b), "{a} vs {b}");
            }
        }
    }
}
