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
#[inline]
fn check_txps_suffix(txp_ids: &[u32], ref_names: &[String]) -> bool {
    txp_ids.iter().any(|&tid| {
        ref_names
            .get(tid as usize)
            .map(|nm| nm.ends_with("-U"))
            .unwrap_or(false)
    })
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

// Forseti path currently wants a list-of-candidates interface (even if length-1).
#[inline]
fn collapse_vertices_multi(
    v: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
) -> Vec<(Vec<u32>, u32)> {
    let (mcc, txp) = collapse_vertices(v, uncovered_vertices, g, eqmap, hasher_state);
    vec![(mcc, txp)]
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
    ref_names: &[String],
    spliceu_txome: &HashMap<u32, Vec<u8>>,
    spline_lookup: &Array1<f64>,
    mlp: &nn::Sequential,
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
    let mut global_txps = get_set(16);
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
            let mut best_mcc: Vec<u32> = Vec::with_capacity(comp_verts.len());
            let mut best_mcc_list: Vec<(Vec<u32>, u32)> = Vec::new();
            let mut dedup_best_mcc_list: Vec<(Vec<u32>, u32)> = Vec::new();
            let mut best_covering_txp_list: Vec<u32> = Vec::new();
            let mut predicted_best_mcc_list: Vec<Vec<u32>> = Vec::new();
            let mut all_global_txps = get_set(16);
            let mut global_txps_check: Vec<u32> = Vec::new();
            let mut sorted_verts: Vec<u32> = Vec::with_capacity(comp_verts.len());

            // we will remove covered vertices from uncovered_vertices until they are
            // all gone (until all vertices have been covered)
            let mut rng = StdRng::seed_from_u64(seed);
            while !uncovered_vertices.is_empty() {
                let num_remaining = uncovered_vertices.len();
                // will hold vertices in the best mcc
                best_mcc.clear();
                best_mcc_list.clear();
                dedup_best_mcc_list.clear();
                best_covering_txp_list.clear();
                predicted_best_mcc_list.clear();

                // the transcript that is responsible for the
                // best mcc covering
                let mut best_covering_txp = u32::MAX;
                let mut find_single_best = false;

                // for each vertex in the vertex set
                // TODO: should we change back later? YG. 
                // NOTE: to make the algorithm deterministic,
                sorted_verts.clear();
                sorted_verts.extend(uncovered_vertices.iter().copied());
                sorted_verts.sort_unstable(); // Sort the vertices to guarantee a fixed order
                for v in sorted_verts.iter() {
                    // find the largest mcc starting from this vertex
                    // and the transcript that covers it
                    // NOTE: mcc is the mcc consists of a set of vertices.(new_mcc.len() = size of the tree )
                    //largest_mcc_list = [(current_mcc, chosen_txp), ...]
                    let largest_mcc_list =
                        collapse_vertices_multi(*v, &uncovered_vertices, g, eqmap, hasher_state);

                    // though we have a list of largest mccs, they are all the same size, only compare the first one with the current best mcc
                    let first_mcc_tuple = &largest_mcc_list[0];
                    let mut new_mcc: Vec<u32> = first_mcc_tuple.0.clone();
                    let covering_txp: u32 = first_mcc_tuple.1;
                    let mcc_len = new_mcc.len();

                    match mcc_len.cmp(&best_mcc.len()) {
                        Ordering::Greater => {
                            // Move `new_mcc` into `best_mcc`
                            best_mcc.clear();
                            best_mcc.extend(new_mcc.iter().copied());
                            best_covering_txp = covering_txp;

                            best_mcc_list.clear();
                            // Use `extend` to add items from `largest_mcc_list`
                            best_mcc_list.extend(
                                largest_mcc_list
                                    .iter()
                                    .map(|(mcc, txp)| (mcc.clone(), *txp)),
                            );
                        }
                        Ordering::Equal => {
                            // Append items when lengths are equal
                            best_mcc_list.extend(
                                largest_mcc_list
                                    .iter()
                                    .map(|(mcc, txp)| (mcc.clone(), *txp)),
                            );
                        }
                        Ordering::Less => {
                            // Do nothing when `mcc_len` is less than `best_mcc.len()`
                        }
                    }
                    // YG. here the best mcc is the largest mcc after we loop all v as a start point.
                    // we can't do better than covering all
                    // remaining uncovered vertices.  So, if we
                    // accomplish that, then quit here.
                    // TODO: Now, if mcc_len == remaing & remaining == 1, we can call forseti to choose the better one. If we have time. we can do this.
                    if mcc_len == num_remaining {
                        break;
                    }
                }
                // YG-TODO: large_mcc is a vec of usize. we can sort it and recognize the same mcc that start from different v but have same combination of v. --> save time by repeat calculating the forseti for these case
                // best_mcc deduplication
                dedup_best_mcc_list.clear();
                for (mcc, txp) in best_mcc_list.iter() {
                    let mut sorted_mcc = mcc.clone();
                    sorted_mcc.sort_unstable();

                    // Only push if we don't have this exact (sorted_mcc, txp) pair yet.
                    if !dedup_best_mcc_list
                        .iter()
                        .any(|(existing_mcc, existing_txp)| {
                            existing_mcc == &sorted_mcc && existing_txp == txp
                        })
                    {
                        dedup_best_mcc_list.push((sorted_mcc, *txp));
                    }
                }

                all_global_txps.clear();
                for (_best_mcc, txp) in dedup_best_mcc_list.iter() {
                    //Hashset will ignore duplicates
                    all_global_txps.insert(*txp);
                }
                global_txps_check.clear();
                global_txps_check.extend(all_global_txps.iter().copied());

                // project each covering transcript to its
                // corresponding gene, and dedup the list
                let mut global_genes: Vec<u32> = if gene_level_eq_map {
                    global_txps_check.clone()
                } else {
                    all_global_txps
                        .iter()
                        .cloned()
                        .map(|i| tid_to_gid[i as usize])
                        .collect()
                };
                // sort since we will be hashing the ordered vector
                global_genes.sort_unstable();
                // dedup as well since we don't care about duplicates
                global_genes.dedup();
                assert!(
                    !global_genes.is_empty(),
                    "can't find representative gene(s) for a molecule"
                );

                // check if the best_mcc_list (have 1 or more mcc) have multiple gene label/ amb splicing status(s.s.). If so, we will call forseti. Otherwise, we will randomly choose one mcc.

                // call forseti if we have multiple best mccs to determine the final best mcc and best covering txps
                let mut need_forseti = false;

                if all_global_txps.len() == 1 {
                    pug_stats.non_triv_mccs_determined_ss.0 += 1;
                } else if global_genes.len() == 1
                    && check_txps_suffix(&global_txps_check, ref_names)
                {
                    pug_stats.non_triv_mccs_determined_ss.1 += 1;
                } else {
                    need_forseti = true;
                }

                if need_forseti {
                    pug_stats.large_mcc_forseti_needed += 1;
                    // best_mcc_list: [(mcc, covering_txp), (mcc2, covering_txp2), ...], where mcc is a vector of vertex indices.
                    // check_mcc_list: [(covering_ref_name, rname_list1),(covering_ref_name2, rname_list2),...]
                    let check_mcc_list: Vec<(u32, Vec<u64>)> =
                        get_check_mcc_list(&dedup_best_mcc_list, eqmap, g);
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
                    );

                    match forseti_result {
                        Ok(best_mcc_indices) => {
                            if best_mcc_indices.len() == 1 {
                                find_single_best = true;
                            }
                            for mcc_idx in &best_mcc_indices {
                                let mcc_info_tuple = &dedup_best_mcc_list[*mcc_idx as usize];
                                let best_covering_txp_element = mcc_info_tuple.1;
                                let predicted_best_mcc = mcc_info_tuple.0.clone();
                                predicted_best_mcc_list.push(predicted_best_mcc);
                                best_covering_txp_list.push(best_covering_txp_element);
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
                    if predicted_best_mcc_list.is_empty() || best_covering_txp_list.is_empty() {
                        // It is possible that the we cannot find a valid score and share txp for all reads in a mcc. In this case, we will have empty result from forseti.
                        //We will randomly choose the mcc from the best_mcc_list(before forseti) and the corresponding txp.)
                        // YG NOTE:
                        // we use the first mcc in the best_mcc_list, which is the same as the original af's implementation. Try to get a no-worse version, serve as a better baseline.
                        let chosen_mcc_info = &best_mcc_list[0];
                        let chosen_mcc = chosen_mcc_info.0.clone();
                        best_mcc = chosen_mcc.clone();
                        best_covering_txp = chosen_mcc_info.1;
                        // We iterate over all vertices in the mcc, and for each one, we keep track of the (monotonically non-increasing) set of transcripts that have appeared in all vertices.
                        for (index, vertex) in chosen_mcc.iter().enumerate() {
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
                    } else {
                        // We got valid better txp and mcc from forseti

                        //Even though forseti output all equally best mcc and txp. we will eventually choose just one "best_mcc" to proceed. here we randomly choose one best mcc.
                        let before_dedup_predicted_best_mcc_list = predicted_best_mcc_list.clone();
                        predicted_best_mcc_list.dedup();
                        if predicted_best_mcc_list.len() == 1 {
                            best_mcc = predicted_best_mcc_list[0].clone();
                            //NOTE: original af, will finally fix one mcc but the tx and corresponding gene could be multiple. then the so called 'global_genes' will be passed to the gene_eqclass_hash.
                            // NOTE: why change?[Jan6, 2025] original af will find and use all possible txp after choose a best mcc, but forseti might narrow down the txp set. we will use the narrowed txp set to update the gene_eqclass_hash.
                            for best_covering_txp in best_covering_txp_list.iter() {
                                global_txps.insert(*best_covering_txp);
                            }
                        } else {
                            // predicted multiple best mccs, randomly choose one
                            //generate random number as idx

                            // Logic: we random choose one mcc, then, instead of uisng its original global_txp, we'll use its winner txp from forseti.
                            let idx = rng.gen_range(0..before_dedup_predicted_best_mcc_list.len());
                            // Select elements at the same random index
                            best_mcc.clear();
                            best_mcc.extend(
                                before_dedup_predicted_best_mcc_list[idx]
                                    .iter()
                                    .copied(),
                            );

                            // Find all indices in predicted_best_mcc_list that match best_mcc
                            let matching_indices: Vec<usize> = before_dedup_predicted_best_mcc_list
                                .iter()
                                .enumerate()
                                .filter_map(
                                    |(i, mcc)| if mcc == &best_mcc { Some(i) } else { None },
                                )
                                .collect();

                            // Insert `best_covering_txp` values directly into `global_txps` using matching indices
                            for &i in &matching_indices {
                                global_txps.insert(best_covering_txp_list[i]);
                            }
                        }
                    }
                    // //check if this mcc have determined splicing status
                    let mut global_txps_check: Vec<u32> = global_txps.iter().cloned().collect();
                    global_txps_check.dedup();
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
                    if global_txps_check.len() == 1 {
                        // only 1 tx found, must have determined splicing status and unique gene
                        pug_stats.non_triv_mccs_determined_ss.0 += 1;
                    } else if global_genes.len() == 1 {
                        if check_txps_suffix(&global_txps_check, ref_names) {
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
                let counter = gene_eqclass_hash.entry(global_genes).or_insert(0);
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
            // let mut best_covering_txp = u32::MAX;
            // let mut best_mcc: Vec<u32> = Vec::new();
            let mut find_single_best: bool = false;
            let mut best_covering_txp_list: Vec<u32> = Vec::new();
            // let mut predicted_best_mcc_list: Vec<Vec<u32>> = Vec::new();
            let mut global_genes: Vec<u32>;

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
                    //YG_change. tx have determined splicing status and unique gene, as only 1 tx found
                    pug_stats.triv_mccs_determined_ss.0 += 1;
                    one_vertex_components[0] += 1;
                    // best_covering_txp_list.push(tl[0]); // no nned to update. we already have the global_genes above and only need the global_genes to be passed to gene_eqclass_hash.
                } else if check_txps_suffix(tl, ref_names) {
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
                let mut best_mcc_list: Vec<(Vec<u32>, u32)> = Vec::new();
                for txp in tl.iter() {
                    best_mcc_list.push((vec![*tv], *txp));
                }
                let check_mcc_list: Vec<(u32, Vec<u64>)> =
                    get_check_mcc_list(&best_mcc_list, eqmap, g);

                let forseti_checking_list = get_forseti_check_list(cell_reads, &check_mcc_list);

                let forseti_result = forseti_for_multi_best(
                    &forseti_checking_list,
                    ref_names,
                    spliceu_txome,
                    spline_lookup,
                    mlp,
                    read_length as u16,
                );
                match forseti_result {
                    Ok(best_mcc_indices) => {
                        //clear the best_covering_txp_list and insert based on the best_mcc_indices
                        best_covering_txp_list.clear();
                        // assign single best = 1(mcc,txp)
                        find_single_best = best_mcc_indices.len() == 1;

                        for mcc_idx in &best_mcc_indices {
                            let mcc_info_tuple = &best_mcc_list[*mcc_idx as usize];
                            let best_covering_txp_element = mcc_info_tuple.1;
                            // For single-vertex, we do not need to assign best_mcc. there is always only 1 mcc.
                            // let best_mcc = mcc_info_tuple.0.clone();
                            // predicted_best_mcc_list.push(best_mcc);
                            best_covering_txp_list.push(best_covering_txp_element);
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
                if !best_covering_txp_list.is_empty() {
                    // update the global genes only if the forseti has valid output.
                    // calc the global_genes
                    let mut global_txps_check: Vec<u32> = best_covering_txp_list.to_vec();
                    global_txps_check.dedup();
                    global_genes.clear();
                    global_genes = best_covering_txp_list
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
                        if check_txps_suffix(&global_txps_check, ref_names) {
                            pug_stats.triv_mccs_determined_ss.1 += 1;
                        } else {
                            pug_stats.triv_single_gene_ambiguous_ss += 1;
                        }
                    } else if global_genes.len() > 1 {
                        pug_stats.ambiguous_mccs += 1;
                    }
                }
            }
            let counter = gene_eqclass_hash.entry(global_genes).or_insert(0);
            *counter += 1;
        }
    }
    pug_stats
}

// check_mcc_list: [(covering_ref_id, read_idx_list1), (covering_ref_id, read_idx_list2), ...]
fn get_check_mcc_list(
    best_mcc_list: &[(Vec<u32>, u32)],
    eqmap: &EqMap,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
) -> Vec<(u32, Vec<u64>)> {
    // Start from the best_mcc_list, get the check_mcc_list(containing covering txp and read_idx_list) required by forseti
    // best_mcc_list: [(mcc, covering_txp)], where mcc is a vector of vertex indices.
    let mut check_mcc_list: Vec<(u32, Vec<u64>)> = Vec::with_capacity(best_mcc_list.len());
    for (mcc, covering_txp_id) in best_mcc_list {
        let mut check_read_idx_list: Vec<u64> = Vec::new();
        for v in mcc.iter() {
            let (eq_id, umi_idx) = g.from_index(*v as usize);
            let eq = &eqmap.eqc_info_forseti[eq_id as usize];
            let (_umi_bit, read_idx_list) = &eq.umis[umi_idx as usize];
            check_read_idx_list.reserve(read_idx_list.len());
            check_read_idx_list.extend(read_idx_list.iter().copied());
        }
        check_mcc_list.push((*covering_txp_id, check_read_idx_list));
    }
    return check_mcc_list;
}

fn get_forseti_check_list(
    cell_reads: &[AlevinFryReadRecordWithPosition],
    check_mcc_list: &[(u32, Vec<u64>)],
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

    for (mcc_idx, (covering_txp_id, read_idx_list)) in check_mcc_list.iter().enumerate() {
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