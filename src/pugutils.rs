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
use crate::umi_general_distance::umi_edit_distance_from_packed_shifted;

use needletail::bitkmer::bitmer_to_bytes;
use triple_accel::levenshtein::levenshtein;

//#[inline]
//fn umi_edit_distance_from_packed(x: u64, y: u64, umi_len: u8) -> usize {
//    let xb = bitmer_to_bytes((x, umi_len));
//    let yb = bitmer_to_bytes((y, umi_len));
//    levenshtein(&xb, &yb) as usize
//}


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
    let umi_len: u8 = 12;
    let mut has_edge = |x: &(u64, u32), y: &(u64, u32)| -> PugEdgeType {
        let hdist = if pug_exact_umi {
            if x.0 == y.0 { 0 } else { usize::MAX }
        } else {
            umi_edit_distance_from_packed_shifted(x.0, y.0, umi_len)
            //umi_edit_distance_from_packed(x.0, y.0, umi_len)
            //afutils::count_diff_2_bit_packed(x.0, y.0)
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
                //PugEdgeType::NoEdge
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
pub fn collapse_vertices(
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

//new model?!
#[derive(Debug, Clone)]
struct CandidateComponent {
    vertices: Vec<u32>,
    score: f64,
    tx_weights: Vec<(u32, f64)>,
    tx_ids: Vec<u32>,
}

#[inline]
fn candidate_txps_from_root_with_gap(
    v: u32,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    root_tx_gap: f64,
) -> Vec<u32> {
    let vert = g.from_index(v as usize);
    let labels = eqmap.refs_for_eqc(vert.0);

    let mut tx_scores: Vec<(u32, f64)> = Vec::with_capacity(labels.len());
    for (tx_index, txp) in labels.iter().enumerate() {
        let ll = vertex_loglik_for_tx(vert.0, vert.1, tx_index, eqmap);
        tx_scores.push((*txp, ll));
    }

    if tx_scores.is_empty() {
        return Vec::new();
    }

    let best_score = tx_scores
        .iter()
        .map(|(_, s)| *s)
        .fold(f64::NEG_INFINITY, f64::max);

    let mut kept: Vec<u32> = tx_scores
        .into_iter()
        .filter(|(_, s)| *s >= best_score - root_tx_gap)
        .map(|(t, _)| t)
        .collect();

    kept.sort_unstable();
    kept.dedup();
    kept
}

#[inline]
fn component_soft_objective_with_filter(
    component: &[u32],
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    lambda: f64,
    tx_score_gap: f64,
) -> (f64, Vec<(u32, f64)>, Vec<u32>) {
    let tx_scores = component_tx_log_scores(component, g, eqmap, hasher_state);

    if tx_scores.is_empty() {
        return (f64::NEG_INFINITY, Vec::new(), Vec::new());
    }

    let best_score = tx_scores
        .iter()
        .map(|(_, s)| *s)
        .fold(f64::NEG_INFINITY, f64::max);

    let filtered: Vec<(u32, f64)> = tx_scores
        .into_iter()
        .filter(|(_, s)| *s >= best_score - tx_score_gap)
        .collect();

    if filtered.is_empty() {
        return (f64::NEG_INFINITY, Vec::new(), Vec::new());
    }

    let vals: Vec<f64> = filtered.iter().map(|(_, s)| *s).collect();
    let obj = logsumexp(&vals) - lambda;
    let tx_weights = softmax_log_scores(&filtered);
    let tx_ids = filtered.iter().map(|(t, _)| *t).collect();

    (obj, tx_weights, tx_ids)
}

#[inline]
fn objective_gain_of_adding(
    current_component: &[u32],
    candidate_vertex: u32,
    txp_seed: u32,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    lambda: f64,
    tx_score_gap: f64,
) -> Option<(f64, f64, Vec<(u32, f64)>, Vec<u32>)> {
    let vert = g.from_index(candidate_vertex as usize);
    let labels = eqmap.refs_for_eqc(vert.0);

    if labels.binary_search(&txp_seed).is_err() {
        return None;
    }

    let (curr_score, _, _) = component_soft_objective_with_filter(
        current_component,
        g,
        eqmap,
        hasher_state,
        lambda,
        tx_score_gap,
    );

    let mut trial = current_component.to_vec();
    trial.push(candidate_vertex);

    let (trial_score, trial_tx_weights, trial_tx_ids) = component_soft_objective_with_filter(
        &trial,
        g,
        eqmap,
        hasher_state,
        lambda,
        tx_score_gap,
    );

    Some((trial_score - curr_score, trial_score, trial_tx_weights, trial_tx_ids))
}

#[inline]
fn posthoc_prune_component(
    component: &[u32],
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    lambda: f64,
    tx_score_gap: f64,
) -> CandidateComponent {
    let mut best_vertices = component.to_vec();

    loop {
        let (curr_score, curr_tx_weights, curr_tx_ids) = component_soft_objective_with_filter(
            &best_vertices,
            g,
            eqmap,
            hasher_state,
            lambda,
            tx_score_gap,
        );

        if best_vertices.len() <= 1 {
            return CandidateComponent {
                vertices: best_vertices,
                score: curr_score,
                tx_weights: curr_tx_weights,
                tx_ids: curr_tx_ids,
            };
        }

        let mut improved = false;
        let mut best_trial_vertices = best_vertices.clone();
        let mut best_trial_score = curr_score;
        let mut best_trial_tx_weights = curr_tx_weights.clone();
        let mut best_trial_tx_ids = curr_tx_ids.clone();

        for idx in 0..best_vertices.len() {
            let mut trial = best_vertices.clone();
            trial.remove(idx);

            let (trial_score, trial_tx_weights, trial_tx_ids) = component_soft_objective_with_filter(
                &trial,
                g,
                eqmap,
                hasher_state,
                lambda,
                tx_score_gap,
            );

            if trial_score > best_trial_score {
                improved = true;
                best_trial_vertices = trial;
                best_trial_score = trial_score;
                best_trial_tx_weights = trial_tx_weights;
                best_trial_tx_ids = trial_tx_ids;
            }
        }

        if improved {
            best_vertices = best_trial_vertices;
        } else {
            return CandidateComponent {
                vertices: best_vertices,
                score: best_trial_score,
                tx_weights: best_trial_tx_weights,
                tx_ids: best_trial_tx_ids,
            };
        }
    }
}

fn collapse_vertices_log_weighted_greedy(
    v: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    lambda: f64,
) -> (Vec<u32>, f64, Vec<(u32, f64)>, Vec<u32>) {
    const ROOT_TX_GAP: f64 = 5.0;
    const TX_SCORE_GAP: f64 = 5.0;
    const MIN_GAIN_ADD: f64 = 1e-8;

    type VertexSet = HashSet<u32, ahash::RandomState>;
    let get_set =
        |cap: u32| VertexSet::with_capacity_and_hasher(cap as usize, hasher_state.clone());

    let mut best_component = vec![v];
    let mut best_score = f64::NEG_INFINITY;
    let mut best_tx_weights = Vec::new();
    let mut best_tx_ids = Vec::new();

    let root_txps = candidate_txps_from_root_with_gap(v, g, eqmap, ROOT_TX_GAP);

    if root_txps.is_empty() {
        let (score, tx_weights, tx_ids) = component_soft_objective_with_filter(
            &[v],
            g,
            eqmap,
            hasher_state,
            lambda,
            TX_SCORE_GAP,
        );
        return (vec![v], score, tx_weights, tx_ids);
    }

    for txp in root_txps.iter() {
        let mut current_component = vec![v];
        let (mut current_score, mut current_tx_weights, _) = component_soft_objective_with_filter(
            &current_component,
            g,
            eqmap,
            hasher_state,
            lambda,
            TX_SCORE_GAP,
        );

        let mut in_component = get_set(uncovered_vertices.len() as u32 + 1);
        in_component.insert(v);

        let mut frontier = get_set(uncovered_vertices.len() as u32 + 1);
        for nv in g.neighbors_directed(g.from_index(v as usize), Outgoing) {
            let n = g.to_index(nv) as u32;
            if uncovered_vertices.contains(&n) && !in_component.contains(&n) {
                frontier.insert(n);
            }
        }

        loop {
            let frontier_list: Vec<u32> = frontier.iter().cloned().collect();

            let mut best_neighbor = None;
            let mut best_gain = f64::NEG_INFINITY;
            let mut best_neighbor_score = f64::NEG_INFINITY;
            let mut best_neighbor_tx_weights = Vec::new();

            for n in frontier_list.iter() {
                if let Some((gain, trial_score, trial_tx_weights, _)) = objective_gain_of_adding(
                    &current_component,
                    *n,
                    *txp,
                    g,
                    eqmap,
                    hasher_state,
                    lambda,
                    TX_SCORE_GAP,
                ) {
                    if gain > best_gain {
                        best_gain = gain;
                        best_neighbor = Some(*n);
                        best_neighbor_score = trial_score;
                        best_neighbor_tx_weights = trial_tx_weights;
                    }
                }
            }

            if best_neighbor.is_none() || best_gain <= MIN_GAIN_ADD {
                break;
            }

            let n = best_neighbor.expect("checked above");
            frontier.remove(&n);
            in_component.insert(n);
            current_component.push(n);
            current_score = best_neighbor_score;
            current_tx_weights = best_neighbor_tx_weights;

            for nv in g.neighbors_directed(g.from_index(n as usize), Outgoing) {
                let m = g.to_index(nv) as u32;
                if uncovered_vertices.contains(&m) && !in_component.contains(&m) {
                    frontier.insert(m);
                }
            }
        }

        let pruned = posthoc_prune_component(
            &current_component,
            g,
            eqmap,
            hasher_state,
            lambda,
            TX_SCORE_GAP,
        );

        if pruned.score > best_score {
            best_component = pruned.vertices;
            best_score = pruned.score;
            best_tx_weights = pruned.tx_weights;
            best_tx_ids = pruned.tx_ids;
        } else if current_score > best_score {
            best_component = current_component;
            best_score = current_score;
            best_tx_weights = current_tx_weights;
            best_tx_ids = component_soft_objective_with_filter(
                &best_component,
                g,
                eqmap,
                hasher_state,
                lambda,
                TX_SCORE_GAP,
            ).2;
        }
    }

    (best_component, best_score, best_tx_weights, best_tx_ids)
}

pub fn get_num_molecules_log_greedy<P: EqClassPayload>(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    tid_to_gid: &[u32],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
    hasher_state: &ahash::RandomState,
    large_graph_thresh: usize,
    lambda: f64,
    log: &slog::Logger,
) -> PugResolutionStatistics {
    type U32Set = HashSet<u32, ahash::RandomState>;
    let get_set = |cap: u32| {
        U32Set::with_capacity_and_hasher(cap as usize, hasher_state.clone())
    };

    let gene_level_eq_map = match eqmap.map_type {
        EqMapType::GeneLevel => true,
        EqMapType::TranscriptLevel => false,
    };

    let comps = weakly_connected_components(g);

    let mut pug_stats = PugResolutionStatistics {
        used_alternative_strategy: false,
        total_mccs: 0u64,
        ambiguous_mccs: 0u64,
        trivial_mccs: 0u64,
    };

    for (_comp_label, comp_verts) in comps.iter() {
        if comp_verts.len() > 1 {
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
            for v in comp_verts.iter() {
                uncovered_vertices.insert(*v);
            }

            while !uncovered_vertices.is_empty() {
                let num_remaining = uncovered_vertices.len();

                let mut best_mcc: Vec<u32> = Vec::new();
                let mut best_mcc_score: f64 = f64::NEG_INFINITY;
                let mut best_mcc_txp_weights: Vec<(u32, f64)> = Vec::new();
                let mut best_mcc_tx_ids: Vec<u32> = Vec::new();
                let mut best_covering_txp = u32::MAX;

                for v in uncovered_vertices.iter() {
                    let (cand_mcc, cand_score, cand_tx_weights, cand_tx_ids) =
                        if P::HAS_PROBS {
                            collapse_vertices_log_weighted_greedy(
                                *v,
                                &uncovered_vertices,
                                g,
                                eqmap,
                                hasher_state,
                                lambda,
                            )
                        } else {
                            let (new_mcc, covering_txp) =
                                collapse_vertices(*v, &uncovered_vertices, g, eqmap, hasher_state);
                            (new_mcc, 0.0, vec![(covering_txp, 1.0)], vec![covering_txp])
                        };

                    let mcc_len = cand_mcc.len();

                    if P::HAS_PROBS {
                        if cand_score > best_mcc_score {
                            best_mcc = cand_mcc;
                            best_mcc_score = cand_score;
                            best_mcc_txp_weights = cand_tx_weights;
                            best_mcc_tx_ids = cand_tx_ids;

                            if let Some((top_txp, _)) = best_mcc_txp_weights
                                .iter()
                                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal))
                            {
                                best_covering_txp = *top_txp;
                            }
                        }
                    } else {
                        if best_mcc.len() < mcc_len {
                            best_mcc = cand_mcc;
                            best_covering_txp = cand_tx_ids[0];
                            best_mcc_tx_ids = cand_tx_ids;
                            best_mcc_txp_weights = vec![(best_covering_txp, 1.0)];
                        }
                    }

                    if mcc_len == num_remaining {
                        break;
                    }
                }

                if best_mcc.is_empty() {
                    crit!(log, "Could not find a valid component");
                    std::process::exit(1);
                }

                let mut global_txps = best_mcc_tx_ids.clone();
                global_txps.sort_unstable();
                global_txps.dedup();

                if global_txps.is_empty() {
                    crit!(log, "Selected component had no supported transcripts");
                    std::process::exit(1);
                }

                let mut global_txp_prob: Vec<f64> = Vec::new();
                if P::HAS_PROBS {
                    let mut txp_prob_temp: Vec<(u32, f64)> = best_mcc_txp_weights
                        .iter()
                        .filter(|(t, _)| global_txps.binary_search(t).is_ok())
                        .cloned()
                        .collect();

                    txp_prob_temp.sort_unstable_by_key(|(t, _)| *t);

                    if txp_prob_temp.is_empty() {
                        crit!(log, "Selected component had no transcript weights");
                        std::process::exit(1);
                    }

                    global_txp_prob = if txp_prob_temp.len() == 1 {
                        vec![1.0]
                    } else {
                        let sum_p: f64 = txp_prob_temp.iter().map(|(_, p)| *p).sum();
                        if !sum_p.is_finite() || sum_p <= 0.0 {
                            panic!("invalid transcript probability normalization");
                        }
                        txp_prob_temp.iter().map(|(_, p)| *p / sum_p).collect()
                    };
                }

                let mut global_genes: Vec<u32> = if gene_level_eq_map {
                    global_txps.clone()
                } else {
                    global_txps
                        .iter()
                        .cloned()
                        .map(|i| tid_to_gid[i as usize])
                        .collect()
                };

                global_genes.sort_unstable();
                global_genes.dedup();

                pug_stats.total_mccs += 1;
                if global_genes.len() > 1 {
                    pug_stats.ambiguous_mccs += 1;
                }

                if best_covering_txp != u32::MAX {
                    let best_covering_gene = if gene_level_eq_map {
                        best_covering_txp
                    } else {
                        tid_to_gid[best_covering_txp as usize]
                    };

                    assert!(
                        global_genes.contains(&best_covering_gene),
                        "best gene {} not in covering set",
                        best_covering_gene
                    );
                }

                assert!(!global_genes.is_empty(), "can't find representative gene(s) for a molecule");

                let eq_label_len = global_genes.len();
                let payload = gene_eqclass_hash
                    .entry(global_genes)
                    .or_insert(P::new(eq_label_len));
                payload.inc();

                if P::HAS_PROBS {
                    assert_eq!(global_txps.len(), global_txp_prob.len());
                    payload.add_probs(&global_txp_prob);
                }

                for rv in best_mcc.iter() {
                    uncovered_vertices.remove(rv);
                }
            }
        } else {
            let tv = comp_verts.first().expect("can't extract first vertex");
            let (eq_id, umi_id) = g.from_index(*tv as usize);
            let tl = eqmap.refs_for_eqc(eq_id);

            let mut tx_scores: Vec<(u32, f64)> = Vec::with_capacity(tl.len());
            for (i, t) in tl.iter().enumerate() {
                let ll = vertex_loglik_for_tx(eq_id, umi_id, i, eqmap);
                tx_scores.push((*t, ll));
            }

            tx_scores.sort_unstable_by_key(|(t, _)| *t);

            let global_txps: Vec<u32> = tx_scores.iter().map(|(t, _)| *t).collect();
            let global_txp_prob: Vec<f64> = if P::HAS_PROBS {
                if tx_scores.len() == 1 {
                    vec![1.0]
                } else {
                    let tx_weights = softmax_log_scores(&tx_scores);
                    tx_weights.iter().map(|(_, p)| *p).collect()
                }
            } else {
                Vec::new()
            };

            let mut global_genes: Vec<u32> = if gene_level_eq_map {
                global_txps.clone()
            } else {
                global_txps
                    .iter()
                    .map(|i| tid_to_gid[*i as usize])
                    .collect()
            };
            global_genes.sort_unstable();
            global_genes.dedup();

            assert!(!global_genes.is_empty(), "can't find representative gene(s) for a molecule");

            pug_stats.total_mccs += 1;
            pug_stats.trivial_mccs += 1;
            if global_genes.len() > 1 {
                pug_stats.ambiguous_mccs += 1;
            }

            let eq_label_len = global_genes.len();
            let payload = gene_eqclass_hash
                .entry(global_genes)
                .or_insert(P::new(eq_label_len));
            payload.inc();

            if P::HAS_PROBS {
                payload.add_probs(&global_txp_prob);
            }
        }
    }

    pug_stats
}

//new model based on the log-likelihood model
#[inline]
fn safe_ln_prob(p: f64) -> f64 {
    // avoid -inf from exact zero
    p.max(1e-12_f64).ln()
}

#[inline]
pub fn vertex_loglik_for_tx(
    eqid: u32,
    umi_idx: u32,
    tx_index: usize,
    eqmap: &EqMap,
) -> f64 {
    let prob_vec = eqmap
        .probs_for_eq_umi_tx(eqid, umi_idx, tx_index)
        .expect("eq and umi should be valid");

    if prob_vec.is_empty() {
        panic!(
            "empty prob_vec for eqid={}, umi_idx={}, tx_index={}",
            eqid, umi_idx, tx_index
        );
    }

    if let Some((_k, &_p)) = prob_vec.iter().enumerate().find(|(_, p)| !p.is_finite()) {
        panic!(
            "NON-FINITE prob_vec entry in vertex_loglik_for_tx for eqid={}, umi_idx={}, tx_index={}",
            eqid, umi_idx, tx_index
        );
    }

    prob_vec.iter().map(|p| safe_ln_prob(*p)).sum()
}

#[inline]
fn logsumexp(vals: &[f64]) -> f64 {
    if vals.is_empty() {
        return f64::NEG_INFINITY;
    }

    let max_val = vals
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    if !max_val.is_finite() {
        return max_val;
    }

    let sum_exp: f64 = vals.iter().map(|v| (*v - max_val).exp()).sum();
    max_val + sum_exp.ln()
}

#[inline]
pub fn softmax_log_scores(tx_scores: &[(u32, f64)]) -> Vec<(u32, f64)> {
    if tx_scores.is_empty() {
        return Vec::new();
    }

    let max_score = tx_scores
        .iter()
        .map(|(_, s)| *s)
        .fold(f64::NEG_INFINITY, f64::max);

    let denom: f64 = tx_scores.iter().map(|(_, s)| (*s - max_score).exp()).sum();

    tx_scores
        .iter()
        .map(|(t, s)| (*t, (*s - max_score).exp() / denom))
        .collect()
}

#[inline]
fn intersect_txps_for_component(
    component: &[u32],
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
) -> HashSet<u32, ahash::RandomState> {
    let mut global_txps =
        HashSet::<u32, ahash::RandomState>::with_hasher(hasher_state.clone());

    for (idx, vertex) in component.iter().enumerate() {
        let vert = g.from_index(*vertex as usize);
        let txps_for_vert = eqmap.refs_for_eqc(vert.0);

        if idx == 0 {
            for t in txps_for_vert {
                global_txps.insert(*t);
            }
        } else {
            global_txps.retain(|t| txps_for_vert.binary_search(t).is_ok());
        }
    }

    global_txps
}

#[inline]
pub fn component_tx_log_scores(
    component: &[u32],
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
) -> Vec<(u32, f64)> {
    let global_txps = intersect_txps_for_component(component, g, eqmap, hasher_state);

    let mut tx_scores: Vec<(u32, f64)> = Vec::with_capacity(global_txps.len());

    for txp in global_txps.iter() {
        let mut score_t = 0.0_f64;

        for vertex in component.iter() {
            let vert = g.from_index(*vertex as usize);
            let labels = eqmap.refs_for_eqc(vert.0);
            let tx_index = labels
                .binary_search(txp)
                .expect("transcript should exist in every vertex of the component");
            score_t += vertex_loglik_for_tx(vert.0, vert.1, tx_index, eqmap);
        }

        tx_scores.push((*txp, score_t));
    }

    tx_scores.sort_unstable_by_key(|(t, _)| *t);
    tx_scores
}

#[inline]
fn component_soft_objective(
    component: &[u32],
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    lambda: f64,
) -> (f64, Vec<(u32, f64)>) {
    let tx_scores = component_tx_log_scores(component, g, eqmap, hasher_state);

    if tx_scores.is_empty() {
        return (f64::NEG_INFINITY, Vec::new());
    }

    let score_vals: Vec<f64> = tx_scores.iter().map(|(_, s)| *s).collect();
    let obj = logsumexp(&score_vals) - lambda;
    let tx_weights = softmax_log_scores(&tx_scores);

    (obj, tx_weights)
}


/// Find a transcript-compatible candidate component starting from vertex `v`,
/// and score it using the soft probabilistic objective:
///
///   Phi(C) = log sum_t exp( Score(C,t) ) - lambda
///
/// where
///
///   Score(C,t) = sum_{vertex in C} sum_{reads in vertex} log p(read | t)
fn collapse_vertices_log_weighted(
    v: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    lambda: f64,
) -> (Vec<u32>, f64, Vec<(u32, f64)>) {
    type VertexSet = HashSet<u32, ahash::RandomState>;
    let get_set =
        |cap: u32| VertexSet::with_capacity_and_hasher(cap as usize, hasher_state.clone());

    let vert = g.from_index(v as usize);
    let nvert = g.node_count();

    let mut best_component: Vec<u32> = Vec::new();
    let mut best_component_score: f64 = f64::NEG_INFINITY;
    let mut best_tx_weights: Vec<(u32, f64)> = Vec::new();

    // For each transcript in the root node, grow a transcript-compatible BFS component.
    for txp in eqmap.refs_for_eqc(vert.0).iter() {
        let mut bfs_list = VecDeque::new();
        bfs_list.push_back(v);

        let mut visited_set = get_set(nvert as u32);
        visited_set.insert(v);

        let mut current_component: Vec<u32> = Vec::new();

        while let Some(cv) = bfs_list.pop_front() {
            current_component.push(cv);

            for nv in g.neighbors_directed(g.from_index(cv as usize), Outgoing) {
                let n = g.to_index(nv) as u32;

                if !uncovered_vertices.contains(&n) {
                    continue;
                }

                if !visited_set.insert(n) {
                    continue;
                }

                let n_labels = eqmap.refs_for_eqc(nv.0);
                if n_labels.binary_search(txp).is_ok() {
                    bfs_list.push_back(n);
                }
            }
        }

        if current_component.is_empty() {
            continue;
        }

        let (cand_score, cand_tx_weights) =
            component_soft_objective(&current_component, g, eqmap, hasher_state, lambda);

        if cand_score > best_component_score {
            best_component = current_component;
            best_component_score = cand_score;
            best_tx_weights = cand_tx_weights;
        }
    }

    (best_component, best_component_score, best_tx_weights)
}



pub fn get_num_molecules_log<P: EqClassPayload>(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    tid_to_gid: &[u32],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
    hasher_state: &ahash::RandomState,
    large_graph_thresh: usize,
    lambda: f64,
    log: &slog::Logger,
) -> PugResolutionStatistics {
    type U32Set = HashSet<u32, ahash::RandomState>;
    let get_set = |cap: u32| {
        U32Set::with_capacity_and_hasher(cap as usize, hasher_state.clone())
    };

    let gene_level_eq_map = match eqmap.map_type {
        EqMapType::GeneLevel => true,
        EqMapType::TranscriptLevel => false,
    };

    let comps = weakly_connected_components(g);

    let mut one_vertex_components: Vec<usize> = vec![0, 0];

    let mut pug_stats = PugResolutionStatistics {
        used_alternative_strategy: false,
        total_mccs: 0u64,
        ambiguous_mccs: 0u64,
        trivial_mccs: 0u64,
    };

    for (_comp_label, comp_verts) in comps.iter() {
        if comp_verts.len() > 1 {
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
            for v in comp_verts.iter() {
                uncovered_vertices.insert(*v);
            }

            while !uncovered_vertices.is_empty() {
                let num_remaining = uncovered_vertices.len();

                let mut best_mcc: Vec<u32> = Vec::new();
                let mut best_mcc_score: f64 = f64::NEG_INFINITY;
                let mut best_mcc_txp_weights: Vec<(u32, f64)> = Vec::new();

                // optional: for sanity/debug only
                let mut best_covering_txp = u32::MAX;

                for v in uncovered_vertices.iter() {
                    if P::HAS_PROBS {
                        let (cand_mcc, cand_score, cand_tx_weights) =
                            collapse_vertices_log_weighted(
                                *v,
                                &uncovered_vertices,
                                g,
                                eqmap,
                                hasher_state,
                                lambda,
                            );

                        let mcc_len = cand_mcc.len();

                        if cand_score > best_mcc_score {
                            best_mcc = cand_mcc;
                            best_mcc_score = cand_score;
                            best_mcc_txp_weights = cand_tx_weights;

                            if let Some((top_txp, _)) = best_mcc_txp_weights
                                .iter()
                                .max_by(|a, b| {
                                    a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal)
                                })
                            {
                                best_covering_txp = *top_txp;
                            }
                        }

                        if mcc_len == num_remaining {
                            break;
                        }
                    } else {
                        let (cand_mcc, cand_txp) =
                            collapse_vertices(*v, &uncovered_vertices, g, eqmap, hasher_state);

                        let mcc_len = cand_mcc.len();

                        if best_mcc.len() < mcc_len {
                            best_mcc = cand_mcc;
                            best_covering_txp = cand_txp;
                        }

                        if mcc_len == num_remaining {
                            break;
                        }
                    }
                }

                if best_mcc.is_empty() {
                    crit!(log, "Could not find a valid component / covering transcript");
                    std::process::exit(1);
                }

                // Intersect transcripts across the selected component
                let global_txps_set =
                    intersect_txps_for_component(&best_mcc, g, eqmap, hasher_state);

                let mut global_txps: Vec<u32> = global_txps_set.iter().cloned().collect();
                global_txps.sort_unstable();

                if global_txps.is_empty() {
                    crit!(log, "Selected component had no common transcripts");
                    std::process::exit(1);
                }

                // Transcript probabilities to pass to EM
                let mut global_txp_prob: Vec<f64> = Vec::new();

                if P::HAS_PROBS {
                    let mut txp_prob_temp: Vec<(u32, f64)> = best_mcc_txp_weights
                        .iter()
                        .filter(|(t, _)| global_txps.binary_search(t).is_ok())
                        .cloned()
                        .collect();

                    txp_prob_temp.sort_unstable_by_key(|(t, _)| *t);

                    if txp_prob_temp.is_empty() {
                        crit!(log, "Selected component had no transcript weights");
                        std::process::exit(1);
                    }

                    global_txp_prob = if txp_prob_temp.len() == 1 {
                        vec![1.0]
                    } else {
                        let sum_p: f64 = txp_prob_temp.iter().map(|(_, p)| *p).sum();
                        if !sum_p.is_finite() || sum_p <= 0.0 {
                            panic!("invalid transcript probability normalization");
                        }
                        txp_prob_temp.iter().map(|(_, p)| *p / sum_p).collect()
                    };
                }

                // Project transcripts to genes
                let mut global_genes: Vec<u32> = if gene_level_eq_map {
                    global_txps.clone()
                } else {
                    global_txps
                        .iter()
                        .cloned()
                        .map(|i| tid_to_gid[i as usize])
                        .collect()
                };

                global_genes.sort_unstable();
                global_genes.dedup();

                pug_stats.total_mccs += 1;
                if global_genes.len() > 1 {
                    pug_stats.ambiguous_mccs += 1;
                }

                if best_covering_txp != u32::MAX {
                    let best_covering_gene = if gene_level_eq_map {
                        best_covering_txp
                    } else {
                        tid_to_gid[best_covering_txp as usize]
                    };

                    assert!(
                        global_genes.contains(&best_covering_gene),
                        "best gene {} not in covering set, shouldn't be possible",
                        best_covering_gene
                    );
                }

                assert!(
                    !global_genes.is_empty(),
                    "can't find representative gene(s) for a molecule"
                );

                let eq_label_len = global_genes.len();
                let payload = gene_eqclass_hash
                    .entry(global_genes)
                    .or_insert(P::new(eq_label_len));
                payload.inc();

                if P::HAS_PROBS {
                    assert_eq!(
                        global_txps.len(),
                        global_txp_prob.len(),
                        "length of global_txps and global_txp_prob should match"
                    );
                    payload.add_probs(&global_txp_prob);
                }

                for rv in best_mcc.iter() {
                    uncovered_vertices.remove(rv);
                }
            }
        } else {
            // Single-vertex connected component
            let tv = comp_verts.first().expect("can't extract first vertex");
            let (eq_id, umi_id) = g.from_index(*tv as usize);
            let tl = eqmap.refs_for_eqc(eq_id);

            let vcindex: usize = if tl.len() == 1 { 0 } else { 1 };
            one_vertex_components[vcindex] += 1;

            let mut global_txps = tl.to_vec();
            global_txps.sort_unstable();

            let mut global_txp_prob: Vec<f64> = Vec::new();

            if P::HAS_PROBS {
                let mut tx_scores: Vec<(u32, f64)> = Vec::with_capacity(tl.len());

                for (i, t) in tl.iter().enumerate() {
                    let ll = vertex_loglik_for_tx(eq_id, umi_id, i, eqmap);
                    tx_scores.push((*t, ll));
                }

                tx_scores.sort_unstable_by_key(|(t, _)| *t);

                global_txp_prob = if tx_scores.len() == 1 {
                    vec![1.0]
                } else {
                    let tx_weights = softmax_log_scores(&tx_scores);
                    tx_weights.iter().map(|(_, p)| *p).collect()
                };
            }

            let mut global_genes: Vec<u32>;

            if gene_level_eq_map {
                global_genes = global_txps.clone();
            } else {
                global_genes = global_txps.iter().map(|i| tid_to_gid[*i as usize]).collect();
                global_genes.sort_unstable();
                global_genes.dedup();
            }

            assert!(
                !global_genes.is_empty(),
                "can't find representative gene(s) for a molecule"
            );

            pug_stats.total_mccs += 1;
            pug_stats.trivial_mccs += 1;
            if global_genes.len() > 1 {
                pug_stats.ambiguous_mccs += 1;
            }

            let eq_label_len = global_genes.len();
            let payload = gene_eqclass_hash
                .entry(global_genes)
                .or_insert(P::new(eq_label_len));
            payload.inc();

            if P::HAS_PROBS {
                payload.add_probs(&global_txp_prob);
            }
        }
    }

    pug_stats
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
        let mut avg_prob_by_vertex: HashMap<u32, f64> = HashMap::new();

        //obtain the average probabilities for this UMI
        let prob_vec = eqmap
            .probs_for_eq_umi_tx(vert.0, vert.1, tx_index)
            .expect("eq and umi should be valid");

        if let Some((k, &p)) = prob_vec.iter().enumerate().find(|(_, p)| !p.is_finite()) {
            panic!(
                "NON-FINITE prob_vec entry - 1"
            );
        }

        if prob_vec.is_empty() {
            panic!("empty prob_vec for eqid={}, umi_idx={}, tx_index={}", vert.0, vert.1, tx_index);
        }

        let avg_prob = prob_vec.iter().sum::<f64>() / prob_vec.len() as f64;
        current_prob.push(avg_prob);

        avg_prob_by_vertex.insert(v, avg_prob);

        // get the next vertex in the BFS
        while let Some(cv) = bfs_list.pop_front() {
            // add it to the arboresence
            current_mcc.push(cv);

            //fetch avg prob of current vertex
            let cv_avg = *avg_prob_by_vertex
                .get(&cv)
                .expect("missing avg prob for current vertex");

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
                if !uncovered_vertices.contains(&n) {
                    continue;
                }

                // get the set of transcripts present in the
                // label of the current node.
                let n_labels = eqmap.refs_for_eqc(nv.0);
                if let Ok(tx_index_nv) = n_labels.binary_search(txp) {
                    //bfs_list.push_back(n);

                    //obtain the average probabilities for this UMI
                    let prob_vec = eqmap
                        .probs_for_eq_umi_tx(nv.0, nv.1, tx_index_nv)
                        .expect("eq and umi should be valid");

                    if let Some((k, &p)) = prob_vec.iter().enumerate().find(|(_, p)| !p.is_finite()) {
                        panic!(
                            "NON-FINITE prob_vec entry - 2");
                    }

                    if prob_vec.is_empty() {
                        panic!("empty prob_vec for eqid={}, umi_idx={}, tx_index={}", nv.0, nv.1, tx_index_nv);
                    }

                    let avg_prob = prob_vec.iter().sum::<f64>() / prob_vec.len() as f64;

                    if !avg_prob.is_finite() {
                        panic!("the avg_prob is not finite in the collapse_vertices_weighted function: {:?}", prob_vec);
                    }

                    // only treat as connected if |cv_avg - avg_prob| < 0.5**
                    if (cv_avg - avg_prob).abs() >= 0.5 {
                        continue; 
                    }

                    // skip if already visited
                    if !visited_set.insert(n) {
                        continue;
                    }

                    // push only after passing threshold**
                    bfs_list.push_back(n);

                    current_prob.push(avg_prob);

                    // store neighbor avg prob**
                    avg_prob_by_vertex.insert(n, avg_prob);
                }
            }
        }

        if current_prob.len() == 0 {
            panic!("current_prob is zero in the collapse_vertices_weighted function");
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

    assert_eq!(eq_txps_prob.len(), eqmap.refs_for_eqc(vert.0).len(), "weigheted_collapsed_txp, length of eq_txps_prob should be the same as the number of transcripts in the equivalence class");

    (highest_prob_mcc, chosen_txp, highest_prob, eq_txps_prob)
}

#[inline]
fn resolve_num_molecules_crlike_from_vec_prefer_ambig<P: EqClassPayload>(
    umi_gene_count_vec: &mut [(u64, u32, u32)],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
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

pub fn get_num_molecules_cell_ranger_like_small<B, R, P: EqClassPayload>(
    cell_chunk: &mut chunk::Chunk<R>,
    tid_to_gid: &[u32],
    _num_genes: usize,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
    sa_model: SplicedAmbiguityModel,
    _log: &slog::Logger,
) where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize + UmiTaggedRecord,
    <R as MappedRecord>::ParsingContext: RecordContext,
    <R as MappedRecord>::ParsingContext: Clone,
    <R as MappedRecord>::ParsingContext: Send,
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

pub fn get_num_molecules_cell_ranger_like<P: EqClassPayload>(
    eq_map: &EqMap,
    tid_to_gid: &[u32],
    _num_genes: usize,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
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
pub fn get_num_molecules_large_component<P: EqClassPayload>(
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

                //eprintln!("@@@@@@@@@@@@@ global txp = {:?}", global_txps);
                //eprintln!("############# best_mcc_txp_probs = {:?}", best_mcc_txp_probs);

                // for long reads obtain the probabiltiy for the global txps
                let mut global_txp_prob: Vec<f64> = vec![];
                if P::HAS_PROBS {

                    if let Some((k, &p)) = best_mcc_txp_probs.iter().enumerate().find(|(_, p)| !p.1.is_finite()) {
                        panic!(
                            "NON-FINITE prob_vec entry -- 1 -- best_mcc_txp_prob"
                        );
                    }

                    let mut txp_prob_temp: Vec<(u32, f64)> = best_mcc_txp_probs
                        .iter()
                        .filter(|(t, _)| global_txps.contains(t))
                        .cloned()
                        .collect();

                    if global_txps.len() != txp_prob_temp.len() {
                        eprintln!("global_txp = {:?}, global_txp_len = {}, best_mcc = {:?}, best_mcc_len = {:?}, txp_prob_temp: {:?}, txp_prob_temp_len = {:?}, best_mcc_txp_probs = {:?}, best_mcc_txp_probs_len = {:?}", 
                        global_txps, global_txps.len(), best_mcc, best_mcc.len(), txp_prob_temp, txp_prob_temp.len(), best_mcc_txp_probs, best_mcc_txp_probs.len());
                        std::process::exit(1);
                    }
                    //assert_eq!(global_txps.len(), txp_prob_temp.len(), "after - length of global txps and their corresponding probs should be the same");

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
                    if let Some((k, &p)) = global_txp_prob.iter().enumerate().find(|(_, p)| !p.is_finite()) {
                        panic!(
                            "NON-FINITE prob_vec entry -- 1 -- at the end"
                        );
                    }
                    assert_eq!(global_txps.len(), global_txp_prob.len(), "length of global txps and their corresponding probs should be the same -- 2 -- at the end");
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
                if let Some((k, &p)) = global_txp_prob.iter().enumerate().find(|(_, p)| !p.is_finite()) {
                        panic!(
                            "NON-FINITE prob_vec entry -- 2 -- at the end"
                        );
                    }
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
