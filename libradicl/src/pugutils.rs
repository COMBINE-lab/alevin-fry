/*
 * Copyright (c) 2020-2021 Rob Patro, Avi Srivastava, Hirak Sarkar, Dongze He, Mohsen Zakeri.
 *
 * This file is part of alevin-fry
 * (see https://github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

extern crate petgraph;
extern crate quickersort;
extern crate slog;
use crate as libradicl;

use self::slog::{crit, warn};
#[allow(unused_imports)]
use ahash::{AHasher, RandomState};
use arrayvec::ArrayVec;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};

use petgraph::prelude::*;
use petgraph::unionfind::*;
use petgraph::visit::NodeIndexable;

use crate::schema::{EqMap, PugResolutionStatistics, SplicedAmbiguityModel};
use crate::utils;

type CcMap = HashMap<u32, Vec<u32>, ahash::RandomState>;

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
    uncovered_vertices: &HashSet<u32>, // the set of vertices already covered
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
) -> (Vec<u32>, u32) {
    // get a new set to hold vertices
    type VertexSet = HashSet<u32, ahash::RandomState>;
    fn get_set(cap: u32) -> VertexSet {
        let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
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

#[inline]
fn resolve_num_molecules_crlike_from_vec_prefer_ambig(
    umi_gene_count_vec: &mut [(u64, u32, u32)],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
) {
    // sort the triplets
    // first on umi
    // then on gene_id
    // then on count
    quickersort::sort(umi_gene_count_vec);

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
            if utils::same_gene(gn, prev_gid, true) {
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
    _with_unspliced: bool,
) {
    // sort the triplets
    // first on umi
    // then on gene_id
    // then on count
    quickersort::sort(umi_gene_count_vec);

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

pub(super) fn get_num_molecules_cell_ranger_like_small(
    cell_chunk: &mut libradicl::Chunk,
    tid_to_gid: &[u32],
    _num_genes: usize,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
    with_unspliced: bool,
    sa_model: SplicedAmbiguityModel,
    _log: &slog::Logger,
) {
    let mut umi_gene_count_vec: Vec<(u64, u32, u32)> = Vec::with_capacity(cell_chunk.nrec as usize);

    // for each record
    for rec in &cell_chunk.reads {
        // get the umi
        let umi = rec.umi;

        // project the transcript ids to gene ids
        let mut gset: Vec<u32> = rec
            .refs
            .iter()
            .map(|tid| tid_to_gid[*tid as usize])
            .collect();
        // and make the gene ids unique
        quickersort::sort(&mut gset[..]);
        gset.dedup();
        for g in &gset {
            umi_gene_count_vec.push((umi, *g, 1));
        }
    }
    match sa_model {
        SplicedAmbiguityModel::WinnerTakeAll => {
            resolve_num_molecules_crlike_from_vec(
                &mut umi_gene_count_vec,
                gene_eqclass_hash,
                with_unspliced,
            );
        }
        SplicedAmbiguityModel::PreferAmbiguity => {
            resolve_num_molecules_crlike_from_vec_prefer_ambig(
                &mut umi_gene_count_vec,
                gene_eqclass_hash,
            );
        }
    }
}

pub(super) fn get_num_molecules_cell_ranger_like(
    eq_map: &EqMap,
    tid_to_gid: &[u32],
    _num_genes: usize,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
    with_unspliced: bool,
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
        quickersort::sort(&mut gset[..]);
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
            resolve_num_molecules_crlike_from_vec(
                &mut umi_gene_count_vec,
                gene_eqclass_hash,
                with_unspliced,
            );
        }
        SplicedAmbiguityModel::PreferAmbiguity => {
            resolve_num_molecules_crlike_from_vec_prefer_ambig(
                &mut umi_gene_count_vec,
                gene_eqclass_hash,
            );
        }
    }
}

pub(super) fn get_num_molecules_trivial_discard_all_ambig(
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
        quickersort::sort(&mut v[..]);
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
    num_genes: usize,
    _log: &slog::Logger,
) -> Vec<u32> {
    let mut counts = vec![0u32; num_genes];

    // TODO: better capacity
    let mut umi_gene_count_vec: Vec<(u64, u32, u32)> = vec![];

    // build a temporary hashmap from each
    // equivalence class id in the current subgraph
    // to the set of (UMI, frequency) pairs contained
    // in the subgraph
    let ts = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut tmp_map = HashMap::<u32, Vec<(u64, u32)>, ahash::RandomState>::with_hasher(ts);

    // for each vertex id in the subgraph
    for vertex_id in vertex_ids {
        // get the corresponding vertex which is
        // an (eq_id, UMI index) pair
        let vert = g.from_index(*vertex_id as usize);
        // add the corresponding (UMI, frequency) pair to the map
        // for this eq_id
        let umis = tmp_map.entry(vert.0).or_insert_with(Vec::new);
        umis.push(eq_map.eqc_info[vert.0 as usize].umis[vert.1 as usize]);
    }

    for (k, v) in tmp_map.iter() {
        // get the (umi, count) pairs
        let umis = v; //&eqinfo.umis;
        let eqid = k; //&eqinfo.eq_num;
                      // project the transcript ids to gene ids
        let mut gset: Vec<u32> = eq_map
            .refs_for_eqc(*eqid)
            .iter()
            .map(|tid| tid_to_gid[*tid as usize])
            .collect();
        // and make the gene ids unique
        quickersort::sort(&mut gset[..]);
        gset.dedup();

        // add every (umi, count), gene pair as a triplet
        // of (umi, gene_id, count) to the output vector
        for umi_ct in umis {
            for g in &gset {
                umi_gene_count_vec.push((umi_ct.0, *g, umi_ct.1));
            }
        }
    }

    // sort the triplets
    // first on umi
    // then on gene_id
    // then on count
    quickersort::sort(&mut umi_gene_count_vec[..]);

    // hold the current umi and gene we are examining
    let mut curr_umi = umi_gene_count_vec.first().expect("cell with no UMIs").0;
    let mut curr_gn = umi_gene_count_vec.first().expect("cell with no UMIs").1;
    // hold the gene id having the max count for this umi
    // and the maximum count value itself
    let mut max_count_gene = 0u32;
    let mut max_count = 0u32;
    // to aggregate the count should a (umi, gene) pair appear
    // more than once
    let mut count_aggr = 0u32;
    // could this UMI be assigned toa best gene or not
    let mut unresolvable = false;
    // to keep track of the current index in the vector
    let mut cidx = 0usize;

    // look over all sorted triplets
    while cidx < umi_gene_count_vec.len() {
        let (umi, gn, ct) = umi_gene_count_vec[cidx];

        // if this umi is different than
        // the one we are processing
        // then decide what action to take
        // on the previous umi
        if umi != curr_umi {
            // if previous was resolvable, add it to the appropriate gene
            if !unresolvable {
                counts[max_count_gene as usize] += 1;
            }

            // the next umi and gene
            curr_umi = umi;
            curr_gn = gn;

            // the next umi will start as resolvable
            unresolvable = false;

            // current gene is current best
            max_count_gene = gn;

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
            // also makes this UMI resolvable
            match count_aggr.cmp(&max_count) {
                Ordering::Greater => {
                    max_count = count_aggr;
                    max_count_gene = gn;
                    unresolvable = false;
                }
                Ordering::Equal => {
                    // if we have a tie for the max count
                    // then the current UMI becomes unresolvable
                    // it will stay this way unless we see a bigger
                    // count for this UMI
                    unresolvable = true;
                }
                Ordering::Less => {
                    // we do nothing
                }
            }
        }

        // if this was the last UMI in the list
        if cidx == umi_gene_count_vec.len() - 1 && !unresolvable {
            counts[max_count_gene as usize] += 1;
        }
        cidx += 1;
    }

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
    num_genes: usize,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, u32, ahash::RandomState>,
    log: &slog::Logger,
) -> PugResolutionStatistics
//,)
{
    type U32Set = HashSet<u32, ahash::RandomState>;
    fn get_set(cap: u32) -> U32Set {
        let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
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
            // are very large.  For components with > 1000 vertices
            // (this should be _very_ rare) we will instead resolve
            // the UMIs in the component using a simpler algorithm.
            if comp_verts.len() > 1000 {
                let mut ng = 0u32;
                let mut numi = 0u32;
                let gene_increments = get_num_molecules_large_component(
                    g, eqmap, comp_verts, tid_to_gid, num_genes, log,
                );
                for (gn, val) in gene_increments.iter().enumerate() {
                    if *val > 0 {
                        let e = gene_eqclass_hash.entry(vec![gn as u32]).or_insert(0);
                        *e += *val;
                        ng += 1;
                        numi += *val;
                    }
                }

                warn!(
                    log,
                    "\n\nfound connected component with {} vertices, \
                    resolved into {} UMIs over {} genes with trivial resolution.\n\n",
                    comp_verts.len(),
                    numi,
                    ng
                );
                pug_stats.used_alternative_strategy = true;
                continue;
            }

            // uncovered_vertices will hold the set of vertices that are
            // *not yet* covered.
            let mut uncovered_vertices = comp_verts.iter().cloned().collect::<HashSet<u32>>();

            // we will remove covered vertices from uncovered_vertices until they are
            // all gone (until all vertices have been covered)
            while !uncovered_vertices.is_empty() {
                let num_remaining = uncovered_vertices.len();
                // will hold vertices in the best mcc
                let mut best_mcc: Vec<u32> = Vec::new();
                // the transcript that is responsible for the
                // best mcc covering
                let mut best_covering_txp = std::u32::MAX;
                // for each vertex in the vertex set
                for v in uncovered_vertices.iter() {
                    // find the largest mcc starting from this vertex
                    // and the transcript that covers it
                    // NOTE: what if there are multiple different mccs that
                    // are equally good? (@k3yavi — I don't think this case
                    // is even handled in the C++ code either).
                    let (new_mcc, covering_txp) =
                        collapse_vertices(*v, &uncovered_vertices, g, eqmap);

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

                // project each covering transcript to its
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

            let mut global_genes: Vec<u32> = tl.iter().map(|i| tid_to_gid[*i as usize]).collect();
            quickersort::sort(&mut global_genes[..]);
            global_genes.dedup();

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
