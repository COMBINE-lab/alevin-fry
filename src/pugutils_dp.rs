use ahash::RandomState;
use petgraph::prelude::*;
use petgraph::visit::NodeIndexable;
use slog::{crit, warn};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};

use crate::eq_class::{EqMap, EqMapType};
use crate::pugutils::{
    collapse_vertices, component_tx_log_scores, get_num_molecules_large_component,
    softmax_log_scores, vertex_loglik_for_tx, weakly_connected_components,
    PugResolutionStatistics,
};
use crate::utils::EqClassPayload;

#[derive(Debug, Clone)]
struct DpSubtreeResult {
    score: f64,
    vertices: Vec<u32>,
}

#[derive(Debug, Clone, Copy)]
pub struct DpConfig {
    /// Reward for including a vertex in the DP-selected component.
    /// Larger values encourage bigger components.
    pub lambda_size: f64,
    /// Require delta(v,t) = score(v,t) - best_alt(v) >= tau_delta
    /// for a vertex to be admissible for transcript t.
    pub tau_delta: f64,
}

#[inline]
fn vertex_loglik_for_named_txp(
    eqid: u32,
    umi_idx: u32,
    txp: u32,
    eqmap: &EqMap,
) -> Option<f64> {
    let labels = eqmap.refs_for_eqc(eqid);
    let tx_index = labels.binary_search(&txp).ok()?;
    Some(vertex_loglik_for_tx(eqid, umi_idx, tx_index, eqmap))
}

#[inline]
fn vertex_tx_delta(
    eqid: u32,
    umi_idx: u32,
    txp: u32,
    eqmap: &EqMap,
) -> Option<f64> {
    let labels = eqmap.refs_for_eqc(eqid);
    let tx_index = labels.binary_search(&txp).ok()?;

    let target_score = vertex_loglik_for_tx(eqid, umi_idx, tx_index, eqmap);

    let mut best_alt = f64::NEG_INFINITY;
    for (i, alt_txp) in labels.iter().enumerate() {
        if *alt_txp == txp {
            continue;
        }
        let s = vertex_loglik_for_tx(eqid, umi_idx, i, eqmap);
        if s > best_alt {
            best_alt = s;
        }
    }

    if best_alt == f64::NEG_INFINITY {
        Some(f64::INFINITY)
    } else {
        Some(target_score - best_alt)
    }
}

#[inline]
fn vertex_is_admissible_for_tx(
    eqid: u32,
    umi_idx: u32,
    txp: u32,
    eqmap: &EqMap,
    tau_delta: f64,
) -> bool {
    match vertex_tx_delta(eqid, umi_idx, txp, eqmap) {
        Some(delta) => delta >= tau_delta,
        None => false,
    }
}

#[inline]
fn dp_compatible_children(
    v: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    txp: u32,
    tau_delta: f64,
) -> Vec<u32> {
    let mut out = Vec::new();

    for nv in g.neighbors_directed(g.from_index(v as usize), Outgoing) {
        let n = g.to_index(nv) as u32;

        if !uncovered_vertices.contains(&n) {
            continue;
        }

        if vertex_is_admissible_for_tx(nv.0, nv.1, txp, eqmap, tau_delta) {
            out.push(n);
        }
    }

    out
}

/// Build a rooted traversal tree for transcript `txp`.
/// This converts the cyclic graph problem into a tree DP problem by
/// assigning each reachable compatible vertex a single parent.
fn build_dp_tree_for_tx(
    root: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    txp: u32,
    tau_delta: f64,
) -> HashMap<u32, Vec<u32>, ahash::RandomState> {
    type U32Set = HashSet<u32, ahash::RandomState>;
    let get_set = |cap: u32| {
        U32Set::with_capacity_and_hasher(cap as usize, hasher_state.clone())
    };

    let mut tree =
        HashMap::<u32, Vec<u32>, ahash::RandomState>::with_hasher(hasher_state.clone());

    let mut visited = get_set(uncovered_vertices.len() as u32 + 1);
    let mut bfs = VecDeque::new();

    visited.insert(root);
    bfs.push_back(root);
    tree.entry(root).or_default();

    while let Some(v) = bfs.pop_front() {
        let children = dp_compatible_children(v, uncovered_vertices, g, eqmap, txp, tau_delta);

        for c in children {
            if visited.insert(c) {
                tree.entry(v).or_default().push(c);
                tree.entry(c).or_default();
                bfs.push_back(c);
            }
        }
    }

    tree
}

/// Tree DP recurrence:
/// dp(v,t) = score_v(t) + lambda_size + sum_child max(0, dp(child,t))
fn dp_best_subtree_for_tx(
    v: u32,
    tree: &HashMap<u32, Vec<u32>, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    txp: u32,
    lambda_size: f64,
    memo: &mut HashMap<u32, DpSubtreeResult, ahash::RandomState>,
) -> DpSubtreeResult {
    if let Some(res) = memo.get(&v) {
        return res.clone();
    }

    let vert = g.from_index(v as usize);
    let score_v = vertex_loglik_for_named_txp(vert.0, vert.1, txp, eqmap)
        .expect("DP called with transcript absent from a tree vertex");

    let mut best_score = score_v + lambda_size;
    let mut best_vertices = vec![v];

    if let Some(children) = tree.get(&v) {
        for child in children.iter() {
            let child_res =
                dp_best_subtree_for_tx(*child, tree, g, eqmap, txp, lambda_size, memo);

            if child_res.score > 0.0 {
                best_score += child_res.score;
                best_vertices.extend(child_res.vertices);
            }
        }
    }

    let res = DpSubtreeResult {
        score: best_score,
        vertices: best_vertices,
    };
    memo.insert(v, res.clone());
    res
}

fn dp_component_from_root_tx(
    root: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    txp: u32,
    tau_delta: f64,
    lambda_size: f64,
) -> DpSubtreeResult {
    let tree = build_dp_tree_for_tx(
        root,
        uncovered_vertices,
        g,
        eqmap,
        hasher_state,
        txp,
        tau_delta,
    );

    let mut memo =
        HashMap::<u32, DpSubtreeResult, ahash::RandomState>::with_hasher(hasher_state.clone());

    dp_best_subtree_for_tx(root, &tree, g, eqmap, txp, lambda_size, &mut memo)
}

/// DP-based analogue of collapse_vertices_log_weighted / greedy.
/// For each transcript in the root, solve the best rooted subtree.
fn collapse_vertices_log_dp(
    v: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    cfg: DpConfig,
) -> (Vec<u32>, f64, Vec<(u32, f64)>, Vec<u32>) {
    let vert = g.from_index(v as usize);

    let mut best_component = vec![v];
    let mut best_score = f64::NEG_INFINITY;
    let mut best_tx_weights: Vec<(u32, f64)> = Vec::new();
    let mut best_tx_ids: Vec<u32> = Vec::new();

    for txp in eqmap.refs_for_eqc(vert.0).iter() {
        if !vertex_is_admissible_for_tx(vert.0, vert.1, *txp, eqmap, cfg.tau_delta) {
            continue;
        }

        let dp_res = dp_component_from_root_tx(
            v,
            uncovered_vertices,
            g,
            eqmap,
            hasher_state,
            *txp,
            cfg.tau_delta,
            cfg.lambda_size,
        );

        if dp_res.vertices.is_empty() {
            continue;
        }

        let tx_scores = component_tx_log_scores(&dp_res.vertices, g, eqmap, hasher_state);
        if tx_scores.is_empty() {
            continue;
        }

        let tx_weights = softmax_log_scores(&tx_scores);
        let tx_ids: Vec<u32> = tx_scores.iter().map(|(t, _)| *t).collect();

        if dp_res.score > best_score {
            best_component = dp_res.vertices;
            best_score = dp_res.score;
            best_tx_weights = tx_weights;
            best_tx_ids = tx_ids;
        }
    }

    // Fallback: keep root alone if every transcript got filtered out.
    if best_score == f64::NEG_INFINITY {
        let singleton = vec![v];
        let tx_scores = component_tx_log_scores(&singleton, g, eqmap, hasher_state);
        let tx_weights = softmax_log_scores(&tx_scores);
        let tx_ids: Vec<u32> = tx_scores.iter().map(|(t, _)| *t).collect();

        let mut singleton_best = f64::NEG_INFINITY;
        for (_, s) in tx_scores.iter() {
            if *s > singleton_best {
                singleton_best = *s;
            }
        }

        return (
            singleton,
            singleton_best + cfg.lambda_size,
            tx_weights,
            tx_ids,
        );
    }

    (best_component, best_score, best_tx_weights, best_tx_ids)
}

pub fn get_num_molecules_log_dp<P: EqClassPayload>(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    tid_to_gid: &[u32],
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
    hasher_state: &ahash::RandomState,
    large_graph_thresh: usize,
    cfg: DpConfig,
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
                let mut best_mcc_score = f64::NEG_INFINITY;
                let mut best_mcc_txp_weights: Vec<(u32, f64)> = Vec::new();
                let mut best_mcc_tx_ids: Vec<u32> = Vec::new();
                let mut best_covering_txp = u32::MAX;

                for v in uncovered_vertices.iter() {
                    let (cand_mcc, cand_score, cand_tx_weights, cand_tx_ids) = if P::HAS_PROBS {
                        collapse_vertices_log_dp(
                            *v,
                            &uncovered_vertices,
                            g,
                            eqmap,
                            hasher_state,
                            cfg,
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
                    crit!(log, "Could not find a valid DP component");
                    std::process::exit(1);
                }

                let mut global_txps = best_mcc_tx_ids.clone();
                global_txps.sort_unstable();
                global_txps.dedup();

                if global_txps.is_empty() {
                    crit!(log, "Selected DP component had no supported transcripts");
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
                        crit!(log, "Selected DP component had no transcript weights");
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