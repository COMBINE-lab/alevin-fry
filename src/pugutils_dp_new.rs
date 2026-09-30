// pugutils_dp_new.rs  —  BF (Bayes Factor) deduplication model
//
// Public API unchanged — quant.rs imports still work:
//   use crate::pugutils_dp_new::{get_num_molecules_log_dp_new, DpConfig};
//
// ── What changed from the previous mean_ll version and WHY ───────────────
//
//  CHANGE 1 — DpConfig: renamed lambda_size → log_prior_odds
//    WHY: lambda had no biological meaning and required manual tuning.
//         log_prior_odds has a clear meaning:
//           0.0  = no prior preference, let data decide (good default)
//           >0   = require stronger evidence to merge (conservative)
//           <0   = allow weaker evidence to merge (permissive)
//         NOTE: in quant.rs, change one line:
//           lambda_size: config.lambda_size  →  log_prior_odds: config.lambda_size
//
//  CHANGE 2 — DpDiagnostics: replaced score_v stats with BF stats
//    WHY: old diagnostics reported mean_ll scores which are no longer
//         used in the merge decision. New diagnostics report log_BF
//         values which are the actual merge criterion.
//
//  CHANGE 3 — REMOVED vertex_loglik_for_named_txp
//    WHY: only used by the old DP recurrence which is now gone.
//         (was already generating an "unused function" warning)
//
//  CHANGE 4 — REMOVED vertex_mean_loglik_for_named_txp
//    WHY: the mean_ll score was used in the old absorption condition
//         (mean_ll + lambda > 0). The BF model does not use this at all.
//         The BF is computed purely from sum_ll ratios, not absolute values.
//
//  CHANGE 5 — ADDED child_log_bf_for_txp
//    WHY: this is the core of the BF model.
//         log_BF = LL(child, t*) - LL(child, best_alt_transcript)
//         This asks: "do the child's reads prefer t* over all alternatives?"
//         Positive → child's reads are more consistent with t* → merge
//         Negative → child prefers a different transcript → reject
//         This replaces the lambda absorption condition entirely.
//
//  CHANGE 6 — CHANGED build_dp_tree_for_tx → build_bf_tree_for_tx
//    WHY: the old tree builder applied the tau gate to BOTH root and children.
//         In the BF model, children are filtered by their log_BF relative
//         to t*, NOT by their own tau gate. The tau gate now applies only
//         to the root (controlling which transcripts are tried as t*).
//         OLD: child enters tree if delta(child, t*) >= tau_delta
//         NEW: child enters tree if log_BF(child, t*) >= log_prior_odds
//
//  CHANGE 7 — REMOVED dp_best_subtree_for_tx (the DP recurrence)
//    WHY: the DP recurrence was needed because not every vertex in the
//         tree was guaranteed to be absorbed — it depended on the score.
//         In the BF model, every vertex that enters the tree has already
//         passed the BF test (log_BF >= log_prior_odds), so they are ALL
//         part of the component. No recursive scoring is needed.
//         This simplifies the algorithm from O(n²) to O(n).
//
//  CHANGE 8 — REMOVED dp_component_from_root_tx
//    WHY: this was just a wrapper calling build_dp_tree + dp_best_subtree.
//         Replaced by the simpler bf_component_from_root_tx.
//
//  CHANGE 9 — ADDED bf_component_from_root_tx
//    WHY: replaces dp_component_from_root_tx. Much simpler — just calls
//         build_bf_tree_for_tx and collects all vertices in the tree.
//         All vertices are already guaranteed to have passed the BF test.
//
//  CHANGE 10 — CHANGED collapse_vertices_log_dp → collapse_vertices_bf
//    WHY: updated to call bf_component_from_root_tx instead of
//         dp_component_from_root_tx. The score used to pick the best
//         component across transcripts is now the component SIZE (number
//         of vertices) rather than the DP score, since all included
//         vertices are equally valid (all passed BF test).
//
//  UNCHANGED: tau gate functions (vertex_tx_delta, vertex_is_admissible_for_tx)
//    These still apply to the ROOT vertex to control which transcripts
//    are tried as t*. This is desirable — we don't want to try a t* that
//    the root's own reads don't support.
//
//  UNCHANGED: entire outer loop, gene assignment, singleton handling,
//    component_tx_log_scores, softmax_log_scores, all public interface.
//
// ─────────────────────────────────────────────────────────────────────────

use petgraph::prelude::*;
use petgraph::visit::NodeIndexable;
use slog::{crit, warn};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};
use crate::umi_general_distance::umi_edit_distance_from_packed_shifted;

use crate::eq_class::{EqMap, EqMapType};
use crate::pugutils::{
    collapse_vertices, component_tx_log_scores, get_num_molecules_large_component,
    softmax_log_scores, vertex_loglik_for_tx, weakly_connected_components,
    PugResolutionStatistics,
};
use crate::utils::EqClassPayload;

use petgraph::unionfind::UnionFind;

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use crate::quant::EMMode;

static ORACLE_CELLS: AtomicU64 = AtomicU64::new(0);
static ORACLE_VERTICES: AtomicU64 = AtomicU64::new(0);
static ORACLE_VNODES: AtomicU64 = AtomicU64::new(0);
static ORACLE_MOLECULES: AtomicU64 = AtomicU64::new(0);
static ORACLE_READS: AtomicU64 = AtomicU64::new(0);
static ORACLE_CASE1_NODES: AtomicU64 = AtomicU64::new(0);
static ORACLE_CASE2_GROUPS: AtomicU64 = AtomicU64::new(0);
static ORACLE_MOL_CASE1: AtomicU64 = AtomicU64::new(0);
static ORACLE_MOL_CASE2: AtomicU64 = AtomicU64::new(0);
static ORACLE_MOL_BOTH: AtomicU64 = AtomicU64::new(0);
static ORACLE_MOL_EITHER: AtomicU64 = AtomicU64::new(0);
static ORACLE_MOL_NEITHER: AtomicU64 = AtomicU64::new(0);
static ORACLE_NODES_CASE1: AtomicU64 = AtomicU64::new(0);
static ORACLE_NODES_CASE2: AtomicU64 = AtomicU64::new(0);
static ORACLE_NODES_BOTH: AtomicU64 = AtomicU64::new(0);
static ORACLE_NODES_EITHER: AtomicU64 = AtomicU64::new(0);

// ─────────────────────────────────────────────────────────────────────────
// CHANGE 2 — Updated diagnostics
// Old version tracked: min_score_v, max_score_v, score_v_count, total_dp_calls
// New version tracks:  bf_min, bf_max, bf_sum, bf_count, total_children_merged
// WHY: the old stats were about mean_ll scores. The new stats are about
//      the actual BF values so you can see the evidence distribution.
// ─────────────────────────────────────────────────────────────────────────

fn diag_enabled() -> bool {
    std::env::var("DP_DIAG").map(|v| v == "1").unwrap_or(false)
}

#[derive(Default)]
struct DpDiagnostics {
    singleton_comps:           u64,
    multi_vertex_comps:        u64,
    max_comp_size:             usize,
    max_prob_vec_len:          usize,
    total_children_evaluated:  u64,
    total_children_merged:     u64,   // CHANGED: was total_children_absorbed
    total_children_rejected:   u64,   // NEW: explicit rejection count
    // CHANGED: was min_score_v / max_score_v / score_v_count
    // Now tracks log_BF distribution so you can see the evidence quality
    bf_min:   f64,
    bf_max:   f64,
    bf_sum:   f64,
    bf_count: u64,
}

impl DpDiagnostics {
    fn new() -> Self {
        Self {
            bf_min: f64::INFINITY,
            bf_max: f64::NEG_INFINITY,
            ..Default::default()
        }
    }

    fn record_bf(&mut self, bf: f64) {
        let capped = if bf.is_finite() { bf } else { 10.0 };
        self.bf_sum   += capped;
        self.bf_count += 1;
        if capped < self.bf_min { self.bf_min = capped; }
        if capped > self.bf_max { self.bf_max = capped; }
    }

    // CHANGED: report now shows BF stats instead of score_v stats
    fn report(&self, cfg: DpConfig) {
        let merge_rate = if self.total_children_evaluated > 0 {
            self.total_children_merged as f64 / self.total_children_evaluated as f64
        } else { 0.0 };
        let bf_mean = if self.bf_count > 0 {
            self.bf_sum / self.bf_count as f64
        } else { 0.0 };

        eprintln!("=== BF_DIAG report (log_prior_odds={:.3}, tau={:.3}) ===",
            cfg.log_prior_odds, cfg.tau_delta);
        eprintln!("  singleton_comps          = {}", self.singleton_comps);
        eprintln!("  multi_vertex_comps       = {}", self.multi_vertex_comps);
        eprintln!("  max_comp_size            = {}", self.max_comp_size);
        eprintln!("  max_prob_vec_len         = {}", self.max_prob_vec_len);
        eprintln!("  total_children_evaluated = {}", self.total_children_evaluated);
        eprintln!("  total_children_merged    = {}", self.total_children_merged);
        eprintln!("  total_children_rejected  = {}", self.total_children_rejected);
        eprintln!("  merge_rate               = {:.1}%", merge_rate * 100.0);
        // CHANGED: was "score_v range". Now shows BF range — more informative.
        // log_BF > 0 means child prefers t*. log_BF < 0 means child prefers alt.
        eprintln!("  log_BF range             = [{:.4}, {:.4}]  mean={:.4}  (n={})",
            self.bf_min, self.bf_max, bf_mean, self.bf_count);
        eprintln!("  merge condition          = log_BF >= {:.3}", cfg.log_prior_odds);
        if self.multi_vertex_comps == 0 {
            eprintln!("  >>> All components are singletons — model has no effect");
        } else if self.total_children_evaluated == 0 {
            eprintln!("  >>> Tau gate blocked all children — try tau_delta=f64::NEG_INFINITY");
        } else if self.total_children_merged == 0 {
            eprintln!("  >>> No merges: log_prior_odds ({:.3}) > max BF ({:.4})",
                cfg.log_prior_odds, self.bf_max);
            eprintln!("      Lower log_prior_odds to get merges");
        }
        eprintln!("======================================================");
    }
}

// ─────────────────────────────────────────────────────────────────────────
// CHANGE 1 — DpConfig: lambda_size renamed to log_prior_odds
//
// WHY: lambda_size was biologically meaningless. log_prior_odds is
//      interpretable: log(P(same molecule) / P(different molecule)).
//      Default 0.0 means equal prior — the data alone decides.
//
// quant.rs change needed (one line only):
//   OLD: let dp_cfg = DpConfig { lambda_size: config.lambda_size, tau_delta: ... };
//   NEW: let dp_cfg = DpConfig { log_prior_odds: config.lambda_size, tau_delta: ... };
// ─────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct DpSubtreeResult {
    score:    f64,
    vertices: Vec<u32>,
}

#[derive(Debug, Clone, Copy)]
pub struct DpConfig {
    /// Merge threshold.  Child merges if log_BF >= log_prior_odds.
    /// 0.0  = merge when child's reads prefer t* at least as much as any alt
    /// >0   = require positive BF evidence (conservative)
    /// <0   = allow child to slightly prefer an alt and still merge (permissive)
    pub log_prior_odds: f64,    // CHANGED: was lambda_size

    /// Admissibility gate for the ROOT vertex only.
    /// Controls which transcripts are tried as consensus t*.
    /// -1.0 = original default.  0.0 = strict.  NEG_INFINITY = try all.
    pub tau_delta: f64,         // UNCHANGED
}

// ─────────────────────────────────────────────────────────────────────────
// CHANGE 3 — REMOVED vertex_loglik_for_named_txp
// WHY: only used by dp_best_subtree_for_tx which is also removed.
//      Was already generating an "unused function" compiler warning.
//
// CHANGE 4 — REMOVED vertex_mean_loglik_for_named_txp
// WHY: the mean_ll score was used in the old absorption condition
//      (mean_ll + lambda > 0). The BF model does not use absolute
//      log-likelihood values — it only uses the RATIO between t* and
//      best_alt. So this function is simply not needed.
// ─────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────
// CHANGE 5 — ADDED child_log_bf_for_txp  (new function, replaces mean_ll)
//
// WHY: this is the BF model's core computation.
//
// The old absorption condition was:
//   mean_ll(child, t*) + lambda > 0
//   = "is the child's average alignment probability good enough in isolation?"
//   Problem: lambda is arbitrary, ignores whether other transcripts fit better.
//
// The new condition is:
//   log_BF = LL(child, t*) - LL(child, best_alt) >= log_prior_odds
//   = "do the child's reads prefer t* over all other transcripts?"
//   This is a RELATIVE test — it compares t* to the best alternative.
//   If t* is clearly the best explanation for the child's reads → merge.
//   If another transcript fits the child better → reject (keep separate).
//
// Special cases:
//   - child has only t* in its eq class → no alternative → BF = +inf → always merge
//   - child has t* but all alternatives have higher LL → BF < 0 → reject
//   - child has t* and all alternatives equally poor → BF > 0 → merge
// ─────────────────────────────────────────────────────────────────────────

#[inline]
fn child_log_bf_for_txp(
    eqid: u32,
    umi_idx: u32,
    txp: u32,           // t* — the parent component's consensus transcript
    eqmap: &EqMap,
    diag: &mut DpDiagnostics,
) -> Option<f64> {
    let labels = eqmap.refs_for_eqc(eqid);
    let tx_index = labels.binary_search(&txp).ok()?;  // None if child lacks t*

    // Record prob vec length for diagnostics
    if let Some(pv) = eqmap.probs_for_eq_umi_tx(eqid, umi_idx, tx_index) {
        if pv.len() > diag.max_prob_vec_len {
            diag.max_prob_vec_len = pv.len();
        }
    }

    // LL(child, t*) — sum of ln(p) across all reads for transcript t*
    let target_ll = vertex_loglik_for_tx(eqid, umi_idx, tx_index, eqmap);

    // LL(child, best_alt) — best LL among all other transcripts
    let mut best_alt_ll = f64::NEG_INFINITY;
    for (i, alt_txp) in labels.iter().enumerate() {
        if *alt_txp == txp { continue; }
        let s = vertex_loglik_for_tx(eqid, umi_idx, i, eqmap);
        if s > best_alt_ll { best_alt_ll = s; }
    }

    let log_bf = if best_alt_ll == f64::NEG_INFINITY {
        // t* is the only transcript for this child — merge is unambiguous
        f64::INFINITY
    } else {
        target_ll - best_alt_ll   // log_BF = LL(t*) - LL(best_alt)
    };

    diag.record_bf(log_bf);
    Some(log_bf)
}

// ─────────────────────────────────────────────────────────────────────────
// UNCHANGED — tau gate helpers (still used for ROOT vertex only)
//
// WHY kept: we still need to filter which transcripts are tried as t*
// for the root. Without this, the algorithm would try every transcript
// as a potential consensus, including ones the root's reads don't support.
// ─────────────────────────────────────────────────────────────────────────

#[inline]
fn vertex_tx_delta(eqid: u32, umi_idx: u32, txp: u32, eqmap: &EqMap) -> Option<f64> {
    let labels = eqmap.refs_for_eqc(eqid);
    let tx_index = labels.binary_search(&txp).ok()?;
    let target_score = vertex_loglik_for_tx(eqid, umi_idx, tx_index, eqmap);

    let mut best_alt = f64::NEG_INFINITY;
    for (i, alt_txp) in labels.iter().enumerate() {
        if *alt_txp == txp { continue; }
        let s = vertex_loglik_for_tx(eqid, umi_idx, i, eqmap);
        if s > best_alt { best_alt = s; }
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

// ─────────────────────────────────────────────────────────────────────────
// CHANGE 6 — REPLACED build_dp_tree_for_tx with build_bf_tree_for_tx
//
// WHY: the old tree builder applied the tau gate to children:
//   OLD: child enters tree if delta(child, t*) >= tau_delta
//        (tau gate = "does child's own best transcript agree with t*?")
//
//   NEW: child enters tree if log_BF(child, t*) >= log_prior_odds
//        (BF test = "do child's reads prefer t* over alternatives?")
//
// This is the RIGHT test for children because:
//   - The tau gate on children asked about the child's OWN preference.
//     A child could have delta >= tau but still prefer a different transcript.
//   - The BF test directly asks "is t* the best explanation for this child?"
//     which is exactly what we want to know for the merge decision.
//
// Also removed: counting total_children_tried (replaced by evaluated/merged/rejected)
// ─────────────────────────────────────────────────────────────────────────

fn build_bf_tree_for_tx(
    root: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    txp: u32,
    log_prior_odds: f64,    // CHANGED: was tau_delta applied to children
    diag: &mut DpDiagnostics,
) -> HashMap<u32, Vec<u32>, ahash::RandomState> {
    type U32Set = HashSet<u32, ahash::RandomState>;
    let get_set = |cap: u32| U32Set::with_capacity_and_hasher(cap as usize, hasher_state.clone());

    let mut tree = HashMap::<u32, Vec<u32>, ahash::RandomState>::with_hasher(hasher_state.clone());
    let mut visited = get_set(uncovered_vertices.len() as u32 + 1);
    let mut bfs = VecDeque::new();

    visited.insert(root);
    bfs.push_back(root);
    tree.entry(root).or_default();

    while let Some(v) = bfs.pop_front() {
        for nv in g.neighbors_directed(g.from_index(v as usize), Outgoing) {
            let n = g.to_index(nv) as u32;
            if !uncovered_vertices.contains(&n) { continue; }
            if visited.contains(&n)             { continue; }

            diag.total_children_evaluated += 1;

            // CHANGED: was vertex_is_admissible_for_tx (tau gate on child)
            // NOW: BF test — merge child iff its reads prefer t* over alternatives
            let should_merge = match child_log_bf_for_txp(nv.0, nv.1, txp, eqmap, diag) {
                Some(log_bf) => log_bf >= log_prior_odds,
                None         => false,  // child doesn't even have t* in its eq class
            };

            if should_merge {
                diag.total_children_merged += 1;
                visited.insert(n);
                tree.entry(v).or_default().push(n);
                tree.entry(n).or_default();
                bfs.push_back(n);
            } else {
                diag.total_children_rejected += 1;
            }
        }
    }

    tree
}

// ─────────────────────────────────────────────────────────────────────────
// CHANGE 7 — REMOVED dp_best_subtree_for_tx  (the recursive DP)
//
// WHY: the DP recurrence existed because not every vertex in the tree
//      was guaranteed to be part of the final component. A vertex was
//      only absorbed if its DP score > 0, i.e. mean_ll + lambda > 0.
//      So the tree could contain vertices that were explored but rejected.
//
//      In the BF model, vertices only enter the tree IF they pass the
//      BF test. So every vertex in the tree is already part of the
//      component. No further selection step is needed.
//      This simplifies the algorithm and removes the O(n) DP recurrence.
//
// CHANGE 8 — REMOVED dp_component_from_root_tx
// WHY: was just a thin wrapper around build_dp_tree + dp_best_subtree.
//      Replaced by the simpler bf_component_from_root_tx below.
// ─────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────
// CHANGE 9 — ADDED bf_component_from_root_tx
//
// WHY: replaces dp_component_from_root_tx. Much simpler because all
//      vertices in the BF tree already passed the merge test — we just
//      collect them all. No DP scoring needed.
// ─────────────────────────────────────────────────────────────────────────

fn bf_component_from_root_tx(
    root: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    txp: u32,
    log_prior_odds: f64,
    diag: &mut DpDiagnostics,
) -> DpSubtreeResult {
    let tree = build_bf_tree_for_tx(
        root, uncovered_vertices, g, eqmap, hasher_state,
        txp, log_prior_odds, diag,
    );

    // CHANGED: was dp_best_subtree_for_tx (recursive DP scoring).
    // NOW: simply collect all vertices — they all passed the BF test.
    let vertices: Vec<u32> = tree.keys().cloned().collect();
    let score    = vertices.len() as f64;  // larger groups preferred when BF ties

    DpSubtreeResult { score, vertices }
}

// ─────────────────────────────────────────────────────────────────────────
// CHANGE 10 — RENAMED collapse_vertices_log_dp → collapse_vertices_bf
//             and updated the internal call from dp_component → bf_component
//
// WHY: the function logic is structurally the same — try each admissible
//      transcript of the root as t*, pick the best component. But the
//      internal call now uses bf_component_from_root_tx instead of
//      dp_component_from_root_tx, which uses the BF test for children.
//
// The score used to compare across transcripts is now component SIZE
// rather than the DP score. This is valid because:
//   - All vertices in any BF component passed the same threshold
//   - A larger component means more evidence for that transcript
//   - Ties are broken by the first transcript tried (deterministic)
// ─────────────────────────────────────────────────────────────────────────

fn collapse_vertices_bf(
    v: u32,
    uncovered_vertices: &HashSet<u32, ahash::RandomState>,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    hasher_state: &ahash::RandomState,
    cfg: DpConfig,
    diag: &mut DpDiagnostics,
) -> (Vec<u32>, f64, Vec<(u32, f64)>, Vec<u32>) {
    let vert = g.from_index(v as usize);

    let mut best_component: Vec<u32> = vec![v];
    let mut best_score               = f64::NEG_INFINITY;
    let mut best_tx_weights          = Vec::new();
    let mut best_tx_ids              = Vec::new();

    for txp in eqmap.refs_for_eqc(vert.0).iter() {
        // Tau gate on ROOT — unchanged, still filters which t* are tried
        if !vertex_is_admissible_for_tx(vert.0, vert.1, *txp, eqmap, cfg.tau_delta) {
            continue;
        }

        // CHANGED: was dp_component_from_root_tx
        let bf_res = bf_component_from_root_tx(
            v, uncovered_vertices, g, eqmap, hasher_state,
            *txp, cfg.log_prior_odds, diag,  // CHANGED: log_prior_odds not lambda_size
        );

        if bf_res.vertices.is_empty() { continue; }

        let tx_scores = component_tx_log_scores(&bf_res.vertices, g, eqmap, hasher_state);
        if tx_scores.is_empty() { continue; }

        let tx_weights = softmax_log_scores(&tx_scores);
        let tx_ids: Vec<u32> = tx_scores.iter().map(|(t, _)| *t).collect();

        if bf_res.score > best_score {
            best_component = bf_res.vertices;
            best_score     = bf_res.score;
            best_tx_weights = tx_weights;
            best_tx_ids    = tx_ids;
        }
    }

    // Fallback — UNCHANGED: if all transcripts filtered by tau gate,
    // keep root as singleton (same behaviour as before)
    if best_score == f64::NEG_INFINITY {
        let singleton   = vec![v];
        let tx_scores   = component_tx_log_scores(&singleton, g, eqmap, hasher_state);
        let tx_weights  = softmax_log_scores(&tx_scores);
        let tx_ids: Vec<u32> = tx_scores.iter().map(|(t, _)| *t).collect();
        let singleton_best = tx_scores.iter()
            .map(|(_, s)| *s)
            .fold(f64::NEG_INFINITY, f64::max);
        return (singleton, singleton_best, tx_weights, tx_ids);
    }

    (best_component, best_score, best_tx_weights, best_tx_ids)
}

// ─────────────────────────────────────────────────────────────────────────
// UNCHANGED — public entry point
// Only changes: calls collapse_vertices_bf instead of collapse_vertices_log_dp
// ─────────────────────────────────────────────────────────────────────────

pub fn get_num_molecules_log_dp_new<P: EqClassPayload>(
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

    let do_diag = diag_enabled();
    let mut diag = DpDiagnostics::new();

    let gene_level_eq_map = match eqmap.map_type {
        EqMapType::GeneLevel  => true,
        EqMapType::TranscriptLevel => false,
    };

    let comps = weakly_connected_components(g);

    let mut pug_stats = PugResolutionStatistics {
        used_alternative_strategy: false,
        total_mccs:     0u64,
        ambiguous_mccs: 0u64,
        trivial_mccs:   0u64,
    };

    for (_comp_label, comp_verts) in comps.iter() {
        if comp_verts.len() > 1 {
            if do_diag {
                diag.multi_vertex_comps += 1;
                if comp_verts.len() > diag.max_comp_size {
                    diag.max_comp_size = comp_verts.len();
                }
            }

            if comp_verts.len() > large_graph_thresh {
                get_num_molecules_large_component(
                    g, eqmap, comp_verts, tid_to_gid, hasher_state,
                    gene_eqclass_hash, log,
                );
                warn!(log,
                    "found connected component with {} vertices; resolved with cr-like resolution.",
                    comp_verts.len()
                );
                pug_stats.used_alternative_strategy = true;
                continue;
            }

            let mut uncovered_vertices = get_set(comp_verts.len() as u32);
            for v in comp_verts.iter() { uncovered_vertices.insert(*v); }

            while !uncovered_vertices.is_empty() {
                let num_remaining = uncovered_vertices.len();

                let mut best_mcc:           Vec<u32>       = Vec::new();
                let mut best_mcc_score                     = f64::NEG_INFINITY;
                let mut best_mcc_txp_weights: Vec<(u32, f64)> = Vec::new();
                let mut best_mcc_tx_ids:    Vec<u32>       = Vec::new();
                let mut best_covering_txp                  = u32::MAX;

                for v in uncovered_vertices.iter() {
                    let (cand_mcc, cand_score, cand_tx_weights, cand_tx_ids) = if P::HAS_PROBS {
                        // CHANGED: was collapse_vertices_log_dp
                        collapse_vertices_bf(
                            *v, &uncovered_vertices, g, eqmap, hasher_state, cfg, &mut diag,
                        )
                    } else {
                        let (new_mcc, covering_txp) =
                            collapse_vertices(*v, &uncovered_vertices, g, eqmap, hasher_state);
                        (new_mcc, 0.0, vec![(covering_txp, 1.0)], vec![covering_txp])
                    };

                    let mcc_len = cand_mcc.len();

                    if P::HAS_PROBS {
                        if cand_score > best_mcc_score {
                            best_mcc             = cand_mcc;
                            best_mcc_score       = cand_score;
                            best_mcc_txp_weights = cand_tx_weights;
                            best_mcc_tx_ids      = cand_tx_ids;

                            if let Some((top_txp, _)) = best_mcc_txp_weights
                                .iter()
                                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal))
                            {
                                best_covering_txp = *top_txp;
                            }
                        }
                    } else {
                        if best_mcc.len() < mcc_len {
                            best_mcc          = cand_mcc;
                            best_covering_txp = cand_tx_ids[0];
                            best_mcc_tx_ids   = cand_tx_ids;
                            best_mcc_txp_weights = vec![(best_covering_txp, 1.0)];
                        }
                    }

                    if mcc_len == num_remaining { break; }
                }

                if best_mcc.is_empty() {
                    crit!(log, "Could not find a valid BF component");
                    std::process::exit(1);
                }

                let mut global_txps = best_mcc_tx_ids.clone();
                global_txps.sort_unstable();
                global_txps.dedup();

                if global_txps.is_empty() {
                    crit!(log, "Selected BF component had no supported transcripts");
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
                        crit!(log, "Selected BF component had no transcript weights");
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
                    global_txps.iter().cloned().map(|i| tid_to_gid[i as usize]).collect()
                };
                global_genes.sort_unstable();
                global_genes.dedup();

                pug_stats.total_mccs += 1;
                if global_genes.len() > 1 { pug_stats.ambiguous_mccs += 1; }

                if best_covering_txp != u32::MAX {
                    let best_covering_gene = if gene_level_eq_map {
                        best_covering_txp
                    } else {
                        tid_to_gid[best_covering_txp as usize]
                    };
                    assert!(
                        global_genes.contains(&best_covering_gene),
                        "best gene {} not in covering set", best_covering_gene
                    );
                }

                assert!(!global_genes.is_empty(),
                    "can't find representative gene(s) for a molecule");

                let eq_label_len = global_genes.len();
                let payload = gene_eqclass_hash
                    .entry(global_genes)
                    .or_insert(P::new(eq_label_len));
                payload.inc();

                if P::HAS_PROBS {
                    assert_eq!(global_txps.len(), global_txp_prob.len());
                    payload.add_probs(&global_txp_prob);
                }

                for rv in best_mcc.iter() { uncovered_vertices.remove(rv); }
            }

        } else {
            // ── UNCHANGED: singleton handling ─────────────────────────────
            if do_diag { diag.singleton_comps += 1; }

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
                global_txps.iter().map(|i| tid_to_gid[*i as usize]).collect()
            };
            global_genes.sort_unstable();
            global_genes.dedup();

            assert!(!global_genes.is_empty(),
                "can't find representative gene(s) for a molecule");

            pug_stats.total_mccs   += 1;
            pug_stats.trivial_mccs += 1;
            if global_genes.len() > 1 { pug_stats.ambiguous_mccs += 1; }

            let eq_label_len = global_genes.len();
            let payload = gene_eqclass_hash
                .entry(global_genes)
                .or_insert(P::new(eq_label_len));
            payload.inc();
            if P::HAS_PROBS { payload.add_probs(&global_txp_prob); }
        }
    }

    if do_diag {//&& (diag.singleton_comps + diag.multi_vertex_comps) < 5000 {
        diag.report(cfg);
    }

    pug_stats
}





pub fn get_num_molecules_true_umi_oracle<P: EqClassPayload>(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    tid_to_gid: &[u32],
    true_umi_table: &Vec<Vec<Vec<Vec<u8>>>>, // [eq_id][umi_idx] -> distinct true UMIs, sorted by frequency desc
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
    log: &slog::Logger,
) -> PugResolutionStatistics {
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

    // ── helper: transcript intersection across a set of vertices (global ids) ──
    let intersect_txps = |verts: &[u32]| -> Option<Vec<u32>> {
        let mut acc: Option<HashSet<u32>> = None;
        for v in verts.iter() {
            let (eq_id, _umi_idx) = g.from_index(*v as usize);
            let label: HashSet<u32> = eqmap.refs_for_eqc(eq_id).iter().cloned().collect();
            acc = Some(match acc {
                None => label,
                Some(mut a) => {
                    a.retain(|t| label.contains(t));
                    a
                }
            });
        }
        match acc {
            Some(s) if !s.is_empty() => {
                let mut v: Vec<u32> = s.into_iter().collect();
                v.sort_unstable();
                Some(v)
            }
            _ => None,
        }
    };

    // ── helper: best-covering-transcript fallback (never returns None) ──
    let best_covering_txps = |verts: &[u32]| -> Vec<u32> {
        let mut support: HashMap<u32, u32> = HashMap::new();
        for v in verts.iter() {
            let (eq_id, _umi_idx) = g.from_index(*v as usize);
            for t in eqmap.refs_for_eqc(eq_id).iter() {
                *support.entry(*t).or_insert(0) += 1;
            }
        }
        let max_support = support.values().max().cloned().unwrap_or(0);
        let mut txps: Vec<u32> = support
            .into_iter()
            .filter(|(_, c)| *c == max_support)
            .map(|(t, _)| t)
            .collect();
        txps.sort_unstable();
        txps
    };

    // ── helper: emit one molecule ──
    let mut emit = |verts_len: usize, global_txps: Vec<u32>, pug_stats: &mut PugResolutionStatistics| {
        let mut global_genes: Vec<u32> = if gene_level_eq_map {
            global_txps.clone()
        } else {
            let mut gs: Vec<u32> = global_txps.iter().map(|t| tid_to_gid[*t as usize]).collect();
            gs.sort_unstable();
            gs.dedup();
            gs
        };
        global_genes.sort_unstable();

        pug_stats.total_mccs += 1;
        if global_genes.len() > 1 {
            pug_stats.ambiguous_mccs += 1;
        }
        if verts_len == 1 {
            pug_stats.trivial_mccs += 1;
        }

        let eq_label_len = global_genes.len();
        let payload = gene_eqclass_hash
            .entry(global_genes)
            .or_insert(P::new(eq_label_len));
        payload.inc();

        if P::HAS_PROBS {
            let n_t = global_txps.len();
            let prob = if n_t == 1 { vec![1.0] } else { vec![1.0 / n_t as f64; n_t] };
            payload.add_probs(&prob);
        }
    };

    for (_comp_label, comp_verts) in comps.iter() {
        let n = comp_verts.len();
        let mut uf = UnionFind::new(n);

        let mut local_idx: HashMap<u32, usize> = HashMap::with_capacity(n);
        for (i, v) in comp_verts.iter().enumerate() {
            local_idx.insert(*v, i);
        }

        // per-node candidate UMI lists (already sorted by frequency desc)
        let node_candidates: Vec<&Vec<Vec<u8>>> = comp_verts
            .iter()
            .map(|v| {
                let (eq_id, umi_idx) = g.from_index(*v as usize);
                &true_umi_table[eq_id as usize][umi_idx as usize]
            })
            .collect();

        // ── pass 1: union using each node's TOP candidate only ──
        for (i, v) in comp_verts.iter().enumerate() {
            for nv in g.neighbors_directed(g.from_index(*v as usize), Outgoing) {
                let n_idx = g.to_index(nv) as u32;
                if let Some(&j) = local_idx.get(&n_idx) {
                    match (node_candidates[i].first(), node_candidates[j].first()) {
                        (Some(a), Some(b)) if a == b => {
                            uf.union(i, j);
                        }
                        _ => {}
                    }
                }
            }
        }

        let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..n {
            let root = uf.find(i);
            groups.entry(root).or_default().push(i);
        }

        for (_root, group_local) in groups.iter() {
            let group_verts: Vec<u32> = group_local.iter().map(|&i| comp_verts[i]).collect();

            if let Some(txps) = intersect_txps(&group_verts) {
                emit(group_verts.len(), txps, &mut pug_stats);
                continue;
            }

            // ── repair attempt: try swapping one ambiguous node at a time
            // to its second-choice UMI, re-cluster, and check if ALL
            // resulting sub-groups get a valid intersection ──
            let ambiguous: Vec<usize> = group_local
                .iter()
                .cloned()
                .filter(|&i| node_candidates[i].len() > 1)
                .collect();

            let mut repaired = false;

            for &amb in &ambiguous {
                // assign top candidate to everyone except `amb`, which gets its 2nd choice
                let assigned: HashMap<usize, &Vec<u8>> = group_local
                    .iter()
                    .filter_map(|&i| {
                        let cand = if i == amb {
                            node_candidates[i].get(1)
                        } else {
                            node_candidates[i].first()
                        };
                        cand.map(|c| (i, c))
                    })
                    .collect();

                // re-cluster this group only, using the new assignment
                let mut sub_uf = UnionFind::new(n);
                let group_set: HashSet<usize> = group_local.iter().cloned().collect();
                for &i in group_local.iter() {
                    let v = comp_verts[i];
                    for nv in g.neighbors_directed(g.from_index(v as usize), Outgoing) {
                        let n_idx = g.to_index(nv) as u32;
                        if let Some(&j) = local_idx.get(&n_idx) {
                            if !group_set.contains(&j) {
                                continue;
                            }
                            if let (Some(a), Some(b)) = (assigned.get(&i), assigned.get(&j)) {
                                if a == b {
                                    sub_uf.union(i, j);
                                }
                            }
                        }
                    }
                }

                let mut subgroups: HashMap<usize, Vec<u32>> = HashMap::new();
                for &i in group_local.iter() {
                    let root = sub_uf.find(i);
                    subgroups.entry(root).or_default().push(comp_verts[i]);
                }

                // check: do ALL resulting sub-groups have a valid intersection?
                let mut all_txps: Vec<(usize, Vec<u32>)> = Vec::new();
                let mut all_valid = true;
                for (root, sverts) in subgroups.iter() {
                    match intersect_txps(sverts) {
                        Some(t) => all_txps.push((*root, t)),
                        None => {
                            all_valid = false;
                            break;
                        }
                    }
                }

                if all_valid {
                    for (root, txps) in all_txps {
                        let sverts_len = subgroups[&root].len();
                        emit(sverts_len, txps, &mut pug_stats);
                    }
                    repaired = true;
                    break; // greedy: accept first swap that fully resolves the group
                }
            }

            if !repaired {
                // safety net: never silently drop a molecule. Fall back to
                // best-covering-transcript for the ORIGINAL (unsplit) group.
                warn!(
                    log,
                    "true-UMI group ({} vertices) could not be resolved by intersection or single-node repair; using best-covering-transcript fallback",
                    group_verts.len()
                );
                let txps = best_covering_txps(&group_verts);
                emit(group_verts.len(), txps, &mut pug_stats);
            }
        }
    }

    pug_stats
}






//#==============================
//pub fn get_num_molecules_true_umi_oracle_virtual_node<P: EqClassPayload>(
//    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
//    eqmap: &EqMap,
//    tid_to_gid: &[u32],
//    read_to_true_umi: &HashMap<String, Vec<u8>>, // read_name -> true UMI bytes
//    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
//    _log: &slog::Logger,
//) -> PugResolutionStatistics {
//    let gene_level_eq_map = match eqmap.map_type {
//        EqMapType::GeneLevel => true,
//        EqMapType::TranscriptLevel => false,
//    };
//
//    let comps = weakly_connected_components(g);
//
//    let mut pug_stats = PugResolutionStatistics {
//        used_alternative_strategy: false,
//        total_mccs: 0u64,
//        ambiguous_mccs: 0u64,
//        trivial_mccs: 0u64,
//    };
//
//    let mut emit = |num_reads: usize, global_txps: Vec<u32>, pug_stats: &mut PugResolutionStatistics| {
//        let mut global_genes: Vec<u32> = if gene_level_eq_map {
//            global_txps.clone()
//        } else {
//            let mut gs: Vec<u32> = global_txps.iter().map(|t| tid_to_gid[*t as usize]).collect();
//            gs.sort_unstable();
//            gs.dedup();
//            gs
//        };
//        global_genes.sort_unstable();
//
//        pug_stats.total_mccs += 1;
//        if global_genes.len() > 1 {
//            pug_stats.ambiguous_mccs += 1;
//        }
//        if num_reads == 1 {
//            pug_stats.trivial_mccs += 1;
//        }
//
//        let eq_label_len = global_genes.len();
//        let payload = gene_eqclass_hash
//            .entry(global_genes)
//            .or_insert(P::new(eq_label_len));
//        payload.inc();
//
//        if P::HAS_PROBS {
//            let n_t = global_txps.len();
//            let prob = if n_t == 1 { vec![1.0] } else { vec![1.0 / n_t as f64; n_t] };
//            payload.add_probs(&prob);
//        }
//    };
//
//    for (_comp_label, comp_verts) in comps.iter() {
//        let n = comp_verts.len();
//
//        let mut local_idx: HashMap<u32, usize> = HashMap::with_capacity(n);
//        for (i, v) in comp_verts.iter().enumerate() {
//            local_idx.insert(*v, i);
//        }
//
//        // ── split each vertex into virtual per-true-UMI entries, weighted
//        // by actual read count. A node with 7 reads of UMI A and 3 of UMI B
//        // becomes two entries, so read-level UMI identity is always exact
//        // -- no top-choice guessing, no repair swaps needed. ──
//        struct VNode {
//            eq_id: u32,
//            umi: Option<Vec<u8>>,
//            count: u32,
//        }
//        let mut vnodes: Vec<VNode> = Vec::new();
//        let mut vnodes_by_vertex: Vec<Vec<usize>> = vec![Vec::new(); n];
//
//        for (i, v) in comp_verts.iter().enumerate() {
//            let (eq_id, umi_idx) = g.from_index(*v as usize);
//            let qnames: &[String] = eqmap.eqc_info[eq_id as usize]
//                .qnames
//                .get(umi_idx as usize)
//                .map(|x| x.as_slice())
//                .unwrap_or(&[]);
//
//            let mut counts: HashMap<Vec<u8>, u32> = HashMap::new();
//            let mut unmatched = 0u32;
//            for qn in qnames {
//                match read_to_true_umi.get(qn) {
//                    Some(tu) => {
//                        *counts.entry(tu.clone()).or_insert(0) += 1;
//                    }
//                    None => unmatched += 1,
//                }
//            }
//
//            for (umi, count) in counts {
//                vnodes_by_vertex[i].push(vnodes.len());
//                vnodes.push(VNode { eq_id, umi: Some(umi), count });
//            }
//            if unmatched > 0 || vnodes_by_vertex[i].is_empty() {
//                // reads with no ground truth still get emitted, just never
//                // merged with anything via UMI equality.
//                vnodes_by_vertex[i].push(vnodes.len());
//                vnodes.push(VNode { eq_id, umi: None, count: unmatched.max(1) });
//            }
//        }
//
//        let num_v = vnodes.len();
//        let mut uf = UnionFind::new(num_v);
//
//        // union two virtual entries iff their parent vertices are edge-
//        // connected in g AND their true UMIs are equal. This is the ONLY
//        // criterion for grouping -- transcript compatibility plays no part
//        // in deciding molecule membership, since true UMI is ground truth.
//        for (i, v) in comp_verts.iter().enumerate() {
//            for nv in g.neighbors_directed(g.from_index(*v as usize), Outgoing) {
//                let n_idx = g.to_index(nv) as u32;
//                if let Some(&j) = local_idx.get(&n_idx) {
//                    for &vi in &vnodes_by_vertex[i] {
//                        for &vj in &vnodes_by_vertex[j] {
//                            if let (Some(a), Some(b)) = (&vnodes[vi].umi, &vnodes[vj].umi) {
//                                if a == b {
//                                    uf.union(vi, vj);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//
//        let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
//        for vi in 0..num_v {
//            let root = uf.find(vi);
//            groups.entry(root).or_default().push(vi);
//        }
//
//        // ── one molecule per true-UMI group, ALWAYS -- never split by
//        // transcript. Gene assignment uses best-covering-transcript
//        // (majority vote across the group's reads, weighted by count),
//        // with ties kept as gene-ambiguous. ──
//        for (_root, members) in groups.iter() {
//            let num_reads: usize = members.iter().map(|&vi| vnodes[vi].count as usize).sum();
//
//            let mut support: HashMap<u32, u32> = HashMap::new();
//            for &vi in members.iter() {
//                let weight = vnodes[vi].count;
//                for t in eqmap.refs_for_eqc(vnodes[vi].eq_id).iter() {
//                    *support.entry(*t).or_insert(0) += weight;
//                }
//            }
//            let max_support = support.values().max().cloned().unwrap_or(0);
//            let mut global_txps: Vec<u32> = support
//                .into_iter()
//                .filter(|(_, c)| *c == max_support)
//                .map(|(t, _)| t)
//                .collect();
//            global_txps.sort_unstable();
//
//            emit(num_reads, global_txps, &mut pug_stats);
//        }
//    }
//
//    pug_stats
//}

//fn extract_blaze_umi(qname: &str) -> Option<&[u8]> {
//    let (_, after_barcode) = qname.split_once('_')?;
//    let (umi, _) = after_barcode.split_once('#')?;
//
//    if umi.is_empty() {
//        None
//    } else {
//        Some(umi.as_bytes())
//    }
//}
//
//fn pack_umi(umi: &[u8]) -> Option<u64> {
//    if umi.len() > 32 {
//        return None;
//    }
//
//    let mut packed = 0u64;
//
//    for &base in umi {
//        let bits = match base.to_ascii_uppercase() {
//            b'A' => 0u64,
//            b'C' => 1u64,
//            b'G' => 2u64,
//            b'T' => 3u64,
//            _ => return None,
//        };
//
//        packed = (packed << 2) | bits;
//    }
//
//    Some(packed)
//}
//
//fn true_umi_diag_enabled() -> bool {
//    std::env::var("TRUE_UMI_DIAG").map(|v| v == "1").unwrap_or(false)
//}
//
//pub fn get_num_molecules_true_umi_oracle_virtual_node<P: EqClassPayload>(
//    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
//    eqmap: &EqMap,
//    tid_to_gid: &[u32],
//    ref_names: &[String],
//    read_to_true_umi: &HashMap<String, Vec<u8>>,
//    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
//    log: &slog::Logger,
//) -> PugResolutionStatistics {
//    let gene_level_eq_map = match eqmap.map_type {
//        EqMapType::GeneLevel => true,
//        EqMapType::TranscriptLevel => false,
//    };
//
//    let do_diag = true_umi_diag_enabled();
//    let comps = weakly_connected_components(g);
//
//    let mut pug_stats = PugResolutionStatistics {
//        used_alternative_strategy: false,
//        total_mccs: 0u64,
//        ambiguous_mccs: 0u64,
//        trivial_mccs: 0u64,
//    };
//
//    let mut emit = |num_reads: usize, global_txps: Vec<u32>, pug_stats: &mut PugResolutionStatistics| {
//        let mut global_genes: Vec<u32> = if gene_level_eq_map {
//            global_txps.clone()
//        } else {
//            let mut gs: Vec<u32> = global_txps.iter().map(|t| tid_to_gid[*t as usize]).collect();
//            gs.sort_unstable();
//            gs.dedup();
//            gs
//        };
//        global_genes.sort_unstable();
//
//        pug_stats.total_mccs += 1;
//        if global_genes.len() > 1 {
//            pug_stats.ambiguous_mccs += 1;
//        }
//        if num_reads == 1 {
//            pug_stats.trivial_mccs += 1;
//        }
//
//        let eq_label_len = global_genes.len();
//        let payload = gene_eqclass_hash
//            .entry(global_genes)
//            .or_insert(P::new(eq_label_len));
//        payload.inc();
//
//        if P::HAS_PROBS {
//            let n_t = global_txps.len();
//            let prob = if n_t == 1 { vec![1.0] } else { vec![1.0 / n_t as f64; n_t] };
//            payload.add_probs(&prob);
//        }
//    };
//
//    for (_comp_label, comp_verts) in comps.iter() {
//        let n = comp_verts.len();
//
//        let mut local_idx: HashMap<u32, usize> = HashMap::with_capacity(n);
//        for (i, v) in comp_verts.iter().enumerate() {
//            local_idx.insert(*v, i);
//        }
//
//        // VNode now also carries the read names backing it, so diagnostics
//        // can print exactly which reads landed where.
//        struct VNode {
//            eq_id: u32,
//            umi_idx: u32,
//            umi: Option<Vec<u8>>,
//            count: u32,
//            qnames: Vec<String>,
//        }
//        let mut vnodes: Vec<VNode> = Vec::new();
//        let mut vnodes_by_vertex: Vec<Vec<usize>> = vec![Vec::new(); n];
//
//        for (i, v) in comp_verts.iter().enumerate() {
//            let (eq_id, umi_idx) = g.from_index(*v as usize);
//            let qnames: &[String] = eqmap.eqc_info[eq_id as usize]
//                .qnames
//                .get(umi_idx as usize)
//                .map(|x| x.as_slice())
//                .unwrap_or(&[]);
//
//            let mut counts: HashMap<Vec<u8>, u32> = HashMap::new();
//            let mut qnames_by_umi: HashMap<Vec<u8>, Vec<String>> = HashMap::new();
//            let mut unmatched_qnames: Vec<String> = Vec::new();
//
//            for qn in qnames {
//                match read_to_true_umi.get(qn) {
//                    Some(tu) => {
//                        *counts.entry(tu.clone()).or_insert(0) += 1;
//                        qnames_by_umi.entry(tu.clone()).or_default().push(qn.clone());
//                    }
//                    None => unmatched_qnames.push(qn.clone()),
//                }
//            }
//
//            for (umi, count) in counts.iter() {
//                vnodes_by_vertex[i].push(vnodes.len());
//                vnodes.push(VNode {
//                    eq_id,
//                    umi_idx,
//                    umi: Some(umi.clone()),
//                    count: *count,
//                    qnames: qnames_by_umi.get(umi).cloned().unwrap_or_default(),
//                });
//            }
//            if !unmatched_qnames.is_empty() || vnodes_by_vertex[i].is_empty() {
//                let cnt = unmatched_qnames.len().max(1) as u32;
//                vnodes_by_vertex[i].push(vnodes.len());
//                vnodes.push(VNode {
//                    eq_id,
//                    umi_idx,
//                    umi: None,
//                    count: cnt,
//                    qnames: unmatched_qnames,
//                });
//            }
//
//            // ── CASE 1 DIAGNOSTIC: this vertex/node had >1 distinct true UMI ──
//            if do_diag && vnodes_by_vertex[i].len() > 1 {
//                eprintln!(
//                    "[CASE1] impure node: eq_id={} umi_idx={} split into {} distinct true UMIs",
//                    eq_id, umi_idx, vnodes_by_vertex[i].len()
//                );
//
//                for &vidx in &vnodes_by_vertex[i] {
//                    let vn = &vnodes[vidx];
//
//                    // Use only the first read associated with this true UMI.
//                    let first_read = vn.qnames.first();
//
//                    let blaze_umi = first_read
//                        .and_then(|qname| extract_blaze_umi(qname));
//
//                    let true_umi = vn.umi.as_deref();
//
//                    let edit_distance: Option<usize> = true_umi
//                        .zip(blaze_umi)
//                        .and_then(|(true_umi, blaze_umi)| {
//                            if true_umi.len() != blaze_umi.len() {
//                                return None;
//                            }
//                        
//                            let true_umi_packed = pack_umi(true_umi)?;
//                            let blaze_umi_packed = pack_umi(blaze_umi)?;
//                            let umi_len = u8::try_from(true_umi.len()).ok()?;
//                        
//                            Some(umi_edit_distance_from_packed_shifted(
//                                blaze_umi_packed,
//                                true_umi_packed,
//                                umi_len,
//                            ))
//                        });
//                    
//                    let true_umi_str = true_umi
//                        .map(|u| String::from_utf8_lossy(u).into_owned())
//                        .unwrap_or_else(|| "<no ground truth>".to_string());
//                    
//                    let blaze_umi_str = blaze_umi
//                        .map(|u| String::from_utf8_lossy(u).into_owned())
//                        .unwrap_or_else(|| "<could not parse>".to_string());
//                    
//                    let edit_distance_str = edit_distance
//                        .map(|d| d.to_string())
//                        .unwrap_or_else(|| "NA".to_string());
//                    
//                    eprintln!(
//                        "  [CASE1] true_umi={} blaze_umi={} edit_distance={} \
//                         count={} reads={:?}",
//                        true_umi_str,
//                        blaze_umi_str,
//                        edit_distance_str,
//                        vn.count,
//                        vn.qnames
//                    );
//                }
//            }
//        }
//
//        let num_v = vnodes.len();
//        let mut uf = UnionFind::new(num_v);
//
//        for (i, v) in comp_verts.iter().enumerate() {
//            for nv in g.neighbors_directed(g.from_index(*v as usize), Outgoing) {
//                let n_idx = g.to_index(nv) as u32;
//                if let Some(&j) = local_idx.get(&n_idx) {
//                    for &vi in &vnodes_by_vertex[i] {
//                        for &vj in &vnodes_by_vertex[j] {
//                            if let (Some(a), Some(b)) = (&vnodes[vi].umi, &vnodes[vj].umi) {
//                                if a == b {
//                                    uf.union(vi, vj);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//
//        let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
//        for vi in 0..num_v {
//            let root = uf.find(vi);
//            groups.entry(root).or_default().push(vi);
//        }
//
//        for (_root, members) in groups.iter() {
//            let num_reads: usize = members.iter().map(|&vi| vnodes[vi].count as usize).sum();
//
//            // ── CASE 2 DIAGNOSTIC: same-true-UMI group with no common transcript ──
//            if do_diag && members.len() > 1 {
//                let mut acc: Option<HashSet<u32>> = None;
//                for &vi in members.iter() {
//                    let label: HashSet<u32> = eqmap.refs_for_eqc(vnodes[vi].eq_id).iter().cloned().collect();
//                    acc = Some(match acc {
//                        None => label,
//                        Some(mut a) => { a.retain(|t| label.contains(t)); a }
//                    });
//                }
//                let empty_intersection = acc.map(|s| s.is_empty()).unwrap_or(true);
//                if empty_intersection {
//                    let umi_str = members.first()
//                        .and_then(|&vi| vnodes[vi].umi.as_ref())
//                        .map(|u| String::from_utf8_lossy(u).to_string())
//                        .unwrap_or_else(|| "<mixed/none>".to_string());
//                    eprintln!(
//                        "[CASE2] true-UMI group ({} vertices, true_umi={}) has NO common transcript",
//                        members.len(), umi_str
//                    );
//                    for &vi in members.iter() {
//                        let vn = &vnodes[vi];
//                        let labels = eqmap.refs_for_eqc(vn.eq_id);
//
//                        let transcript_names: Vec<&str> = labels
//                            .iter()
//                            .map(|&tid| {
//                                ref_names
//                                    .get(tid as usize)
//                                    .map(String::as_str)
//                                    .unwrap_or("<invalid transcript ID>")
//                            })
//                            .collect();
//
//                        eprintln!(
//                            "  [CASE2]   eq_id={} umi_idx={} count={} labels={:?} transcript_names={:?} reads={:?}",
//                            vn.eq_id, vn.umi_idx, vn.count, labels, transcript_names, vn.qnames
//                        );
//                    }
//                }
//            }
//
//            let mut support: HashMap<u32, u32> = HashMap::new();
//            for &vi in members.iter() {
//                let weight = vnodes[vi].count;
//                for t in eqmap.refs_for_eqc(vnodes[vi].eq_id).iter() {
//                    *support.entry(*t).or_insert(0) += weight;
//                }
//            }
//            let max_support = support.values().max().cloned().unwrap_or(0);
//            let mut global_txps: Vec<u32> = support
//                .into_iter()
//                .filter(|(_, c)| *c == max_support)
//                .map(|(t, _)| t)
//                .collect();
//            global_txps.sort_unstable();
//
//            emit(num_reads, global_txps, &mut pug_stats);
//        }
//    }
//
//    pug_stats
//}





/// Print the run-wide oracle totals. Call once from do_quantify after all
/// worker threads have joined.
pub fn report_true_umi_oracle_totals() {
    let l = |c: &AtomicU64| c.load(AtomicOrdering::Relaxed);
    eprintln!("=== TRUE_UMI_ORACLE totals ===");
    eprintln!("  barcodes processed  = {}", l(&ORACLE_CELLS));
    eprintln!("  graph vertices      = {}", l(&ORACLE_VERTICES));
    eprintln!("  virtual sub-nodes   = {}", l(&ORACLE_VNODES));
    eprintln!("  molecules emitted   = {}", l(&ORACLE_MOLECULES));
    eprintln!("  reads in molecules  = {}", l(&ORACLE_READS));
    eprintln!("  CASE1 impure nodes  = {}", l(&ORACLE_CASE1_NODES));
    eprintln!("  CASE2 groups        = {}", l(&ORACLE_CASE2_GROUPS));
    eprintln!("  molecules w/ case1  = {}", l(&ORACLE_MOL_CASE1));
    eprintln!("  molecules w/ case2  = {}", l(&ORACLE_MOL_CASE2));
    eprintln!("  molecules w/ both   = {}", l(&ORACLE_MOL_BOTH));
    eprintln!("  molecules w/ either = {}", l(&ORACLE_MOL_EITHER));
    eprintln!("  molecules w/ neither= {}", l(&ORACLE_MOL_NEITHER));
    eprintln!("  nodes in case1 mols = {}", l(&ORACLE_NODES_CASE1));
    eprintln!("  nodes in case2 mols = {}", l(&ORACLE_NODES_CASE2));
    eprintln!("  nodes in both mols  = {}", l(&ORACLE_NODES_BOTH));
    eprintln!("  nodes in either mols= {}", l(&ORACLE_NODES_EITHER));
    eprintln!("==============================");
}

fn extract_blaze_umi(qname: &str) -> Option<&[u8]> {
    let (_, after_barcode) = qname.split_once('_')?;
    let (umi, _) = after_barcode.split_once('#')?;
    if umi.is_empty() {
        None
    } else {
        Some(umi.as_bytes())
    }
}

fn pack_umi(umi: &[u8]) -> Option<u64> {
    if umi.len() > 32 {
        return None;
    }
    let mut packed = 0u64;
    for &base in umi {
        let bits = match base.to_ascii_uppercase() {
            b'A' => 0u64,
            b'C' => 1u64,
            b'G' => 2u64,
            b'T' => 3u64,
            _ => return None,
        };
        packed = (packed << 2) | bits;
    }
    Some(packed)
}

fn true_umi_diag_enabled() -> bool {
    std::env::var("TRUE_UMI_DIAG").map(|v| v == "1").unwrap_or(false)
}

// ─────────────────────────────────────────────────────────────────────────
// DirectCollapse oracle (virtual nodes).
//
// Each graph vertex is split into one virtual sub-node per distinct true UMI
// among its reads, so read-level UMI identity is exact — no top-choice
// guessing. Grouping is decided ONLY by true-UMI equality along existing
// graph edges; transcript compatibility plays no part in membership.
//
// Transcript set per molecule:
//   non-CASE2 → every transcript common to all members (the intersection)
//   CASE2     → the full union of every member's labels
//
// If ground-truth transcripts are supplied, the highest-read-count true
// transcript replaces the set entirely (probability 1.0). Otherwise the
// probability vector follows prob_model.
// ─────────────────────────────────────────────────────────────────────────

pub fn get_num_molecules_true_umi_oracle_virtual_node<P: EqClassPayload>(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
    tid_to_gid: &[u32],
    ref_names: &[String],
    bc: u64,
    read_to_true_umi: &HashMap<String, Vec<u8>>,
    read_to_true_txp: Option<&Arc<HashMap<String, String>>>,
    prob_model: EMMode,
    gene_eqclass_hash: &mut HashMap<Vec<u32>, P, ahash::RandomState>,
    _log: &slog::Logger,
) -> PugResolutionStatistics {
    let gene_level_eq_map = match eqmap.map_type {
        EqMapType::GeneLevel => true,
        EqMapType::TranscriptLevel => false,
    };

    let do_diag = true_umi_diag_enabled();

    // transcript name -> tid, only built when true transcripts are supplied
    let txp_name_to_id: HashMap<&str, u32> = if read_to_true_txp.is_some() {
        ref_names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i as u32))
            .collect()
    } else {
        HashMap::new()
    };

    let comps = weakly_connected_components(g);

    let mut pug_stats = PugResolutionStatistics {
        used_alternative_strategy: false,
        total_mccs: 0u64,
        ambiguous_mccs: 0u64,
        trivial_mccs: 0u64,
    };

    // per-barcode counters
    let mut cell_vertices = 0u64;
    let mut cell_vnodes = 0u64;
    let mut cell_molecules = 0u64;
    let mut cell_reads = 0u64;
    let mut cell_case1 = 0u64;
    let mut cell_case2 = 0u64;
    let mut cell_mol_case1 = 0u64;
    let mut cell_mol_case2 = 0u64;
    let mut cell_mol_both = 0u64;
    let mut cell_mol_either = 0u64;
    let mut cell_mol_neither = 0u64;
    let mut cell_nodes_case1 = 0u64;
    let mut cell_nodes_case2 = 0u64;
    let mut cell_nodes_both = 0u64;
    let mut cell_nodes_either = 0u64;

    let mut emit = |num_reads: usize,
                    global_txps: Vec<u32>,
                    global_txp_prob: Vec<f64>,
                    pug_stats: &mut PugResolutionStatistics| {
        let mut global_genes: Vec<u32> = if gene_level_eq_map {
            global_txps.clone()
        } else {
            let mut gs: Vec<u32> = global_txps.iter().map(|t| tid_to_gid[*t as usize]).collect();
            gs.sort_unstable();
            gs.dedup();
            gs
        };
        global_genes.sort_unstable();

        pug_stats.total_mccs += 1;
        if global_genes.len() > 1 {
            pug_stats.ambiguous_mccs += 1;
        }
        if num_reads == 1 {
            pug_stats.trivial_mccs += 1;
        }

        let eq_label_len = global_genes.len();
        let payload = gene_eqclass_hash
            .entry(global_genes)
            .or_insert(P::new(eq_label_len));
        payload.inc();

        if P::HAS_PROBS {
            assert_eq!(global_txps.len(), global_txp_prob.len());
            payload.add_probs(&global_txp_prob);
        }
    };

    for (_comp_label, comp_verts) in comps.iter() {
        let n = comp_verts.len();
        cell_vertices += n as u64;

        let mut local_idx: HashMap<u32, usize> = HashMap::with_capacity(n);
        for (i, v) in comp_verts.iter().enumerate() {
            local_idx.insert(*v, i);
        }

        struct VNode {
            eq_id: u32,
            umi_idx: u32,
            umi: Option<Vec<u8>>,
            count: u32,
            qnames: Vec<String>,
        }
        let mut vnodes: Vec<VNode> = Vec::new();
        let mut vnodes_by_vertex: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut vertex_is_impure: Vec<bool> = vec![false; n];

        for (i, v) in comp_verts.iter().enumerate() {
            let (eq_id, umi_idx) = g.from_index(*v as usize);
            let qnames: &[String] = eqmap.eqc_info[eq_id as usize]
                .qnames
                .get(umi_idx as usize)
                .map(|x| x.as_slice())
                .unwrap_or(&[]);

            let mut counts: HashMap<Vec<u8>, u32> = HashMap::new();
            let mut qnames_by_umi: HashMap<Vec<u8>, Vec<String>> = HashMap::new();
            let mut unmatched_qnames: Vec<String> = Vec::new();

            for qn in qnames {
                match read_to_true_umi.get(qn) {
                    Some(tu) => {
                        *counts.entry(tu.clone()).or_insert(0) += 1;
                        qnames_by_umi.entry(tu.clone()).or_default().push(qn.clone());
                    }
                    None => unmatched_qnames.push(qn.clone()),
                }
            }

            for (umi, count) in counts.iter() {
                vnodes_by_vertex[i].push(vnodes.len());
                vnodes.push(VNode {
                    eq_id,
                    umi_idx,
                    umi: Some(umi.clone()),
                    count: *count,
                    qnames: qnames_by_umi.get(umi).cloned().unwrap_or_default(),
                });
            }
            if !unmatched_qnames.is_empty() || vnodes_by_vertex[i].is_empty() {
                let cnt = unmatched_qnames.len().max(1) as u32;
                vnodes_by_vertex[i].push(vnodes.len());
                vnodes.push(VNode {
                    eq_id,
                    umi_idx,
                    umi: None,
                    count: cnt,
                    qnames: unmatched_qnames,
                });
            }

            // ── CASE 1: this vertex/node had >1 distinct true UMI ──
            if vnodes_by_vertex[i].len() > 1 {
                cell_case1 += 1;
                vertex_is_impure[i] = true;

                if do_diag {
                    eprintln!(
                        "[CASE1] bc={} impure node: eq_id={} umi_idx={} split into {} distinct true UMIs",
                        bc, eq_id, umi_idx, vnodes_by_vertex[i].len()
                    );
                    for &vidx in &vnodes_by_vertex[i] {
                        let vn = &vnodes[vidx];
                        let first_read = vn.qnames.first();
                        let blaze_umi = first_read.and_then(|qname| extract_blaze_umi(qname));
                        let true_umi = vn.umi.as_deref();
                        let edit_distance: Option<usize> =
                            true_umi.zip(blaze_umi).and_then(|(true_umi, blaze_umi)| {
                                if true_umi.len() != blaze_umi.len() {
                                    return None;
                                }
                                let true_umi_packed = pack_umi(true_umi)?;
                                let blaze_umi_packed = pack_umi(blaze_umi)?;
                                let umi_len = u8::try_from(true_umi.len()).ok()?;
                                Some(umi_edit_distance_from_packed_shifted(
                                    blaze_umi_packed,
                                    true_umi_packed,
                                    umi_len,
                                ))
                            });

                        let true_umi_str = true_umi
                            .map(|u| String::from_utf8_lossy(u).into_owned())
                            .unwrap_or_else(|| "<no ground truth>".to_string());
                        let blaze_umi_str = blaze_umi
                            .map(|u| String::from_utf8_lossy(u).into_owned())
                            .unwrap_or_else(|| "<could not parse>".to_string());
                        let edit_distance_str = edit_distance
                            .map(|d| d.to_string())
                            .unwrap_or_else(|| "NA".to_string());

                        eprintln!(
                            "  [CASE1] bc={} true_umi={} blaze_umi={} edit_distance={} count={} reads={:?}",
                            bc, true_umi_str, blaze_umi_str, edit_distance_str, vn.count, vn.qnames
                        );
                    }
                }
            }
        }

        let num_v = vnodes.len();
        cell_vnodes += num_v as u64;

        // virtual sub-node -> its parent vertex (local index)
        let mut vnode_parent: Vec<usize> = vec![0usize; num_v];
        for (i, vidxs) in vnodes_by_vertex.iter().enumerate() {
            for &vidx in vidxs {
                vnode_parent[vidx] = i;
            }
        }

        // union two virtual entries iff their parent vertices are edge-connected
        // in g AND their true UMIs are equal
        let mut uf = UnionFind::new(num_v);
        for (i, v) in comp_verts.iter().enumerate() {
            for nv in g.neighbors_directed(g.from_index(*v as usize), Outgoing) {
                let n_idx = g.to_index(nv) as u32;
                if let Some(&j) = local_idx.get(&n_idx) {
                    for &vi in &vnodes_by_vertex[i] {
                        for &vj in &vnodes_by_vertex[j] {
                            if let (Some(a), Some(b)) = (&vnodes[vi].umi, &vnodes[vj].umi) {
                                if a == b {
                                    uf.union(vi, vj);
                                }
                            }
                        }
                    }
                }
            }
        }

        let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
        for vi in 0..num_v {
            let root = uf.find(vi);
            groups.entry(root).or_default().push(vi);
        }

        for (_root, members) in groups.iter() {
            let num_reads: usize = members.iter().map(|&vi| vnodes[vi].count as usize).sum();

            // ── intersection across all members: drives resolution,
            //    not just the diagnostic ──
            let mut acc: Option<HashSet<u32>> = None;
            for &vi in members.iter() {
                let label: HashSet<u32> =
                    eqmap.refs_for_eqc(vnodes[vi].eq_id).iter().cloned().collect();
                acc = Some(match acc {
                    None => label,
                    Some(mut a) => {
                        a.retain(|t| label.contains(t));
                        a
                    }
                });
            }
            let mut common_txps: Vec<u32> =
                acc.map(|s| s.into_iter().collect()).unwrap_or_default();
            common_txps.sort_unstable();
            let is_case2 = common_txps.is_empty();

            // ── CASE 2: same-true-UMI group with no common transcript ──
            if members.len() > 1 && is_case2 {
                cell_case2 += 1;

                if do_diag {
                    let umi_str = members
                        .first()
                        .and_then(|&vi| vnodes[vi].umi.as_ref())
                        .map(|u| String::from_utf8_lossy(u).to_string())
                        .unwrap_or_else(|| "<mixed/none>".to_string());
                    eprintln!(
                        "[CASE2] bc={} true-UMI group ({} vertices, true_umi={}) has NO common transcript",
                        bc, members.len(), umi_str
                    );
                    for &vi in members.iter() {
                        let vn = &vnodes[vi];
                        let labels = eqmap.refs_for_eqc(vn.eq_id);
                        let transcript_names: Vec<&str> = labels
                            .iter()
                            .map(|&tid| {
                                ref_names
                                    .get(tid as usize)
                                    .map(String::as_str)
                                    .unwrap_or("<invalid transcript ID>")
                            })
                            .collect();
                        eprintln!(
                            "  [CASE2] bc={} eq_id={} umi_idx={} count={} labels={:?} transcript_names={:?} reads={:?}",
                            bc, vn.eq_id, vn.umi_idx, vn.count, labels, transcript_names, vn.qnames
                        );
                    }
                }
            }

            // ── transcript set ──
            // non-CASE2: every transcript common to all members
            // CASE2:     the full union of every member's labels
            let mut global_txps: Vec<u32> = if !is_case2 {
                common_txps
            } else {
                let mut v: Vec<u32> = members
                    .iter()
                    .flat_map(|&vi| eqmap.refs_for_eqc(vnodes[vi].eq_id).iter().cloned())
                    .collect();
                v.sort_unstable();
                v.dedup();
                v
            };

            // ── ground-truth transcript override: highest read count wins ──
            if let Some(lookup) = read_to_true_txp {
                let mut txp_votes: HashMap<u32, u32> = HashMap::new();
                for &vi in members.iter() {
                    for qn in vnodes[vi].qnames.iter() {
                        if let Some(name) = lookup.get(qn) {
                            if let Some(&tid) = txp_name_to_id.get(name.as_str()) {
                                *txp_votes.entry(tid).or_insert(0) += 1;
                            }
                        }
                    }
                }
                if do_diag && txp_votes.len() > 1 {
                    eprintln!(
                        "[TRUE_TXP_SPLIT] bc={} {} distinct true transcripts in one true-UMI group: {:?}",
                        bc,
                        txp_votes.len(),
                        txp_votes
                            .iter()
                            .map(|(t, c)| (
                                ref_names.get(*t as usize).map(String::as_str).unwrap_or("?"),
                                *c
                            ))
                            .collect::<Vec<_>>()
                    );
                }
                // tie-break on the lowest transcript id so the result is
                // deterministic across runs (HashMap iteration order is not)
                if let Some((&tid, _)) = txp_votes
                    .iter()
                    .max_by_key(|(t, c)| (**c, std::cmp::Reverse(**t)))
                {
                    if do_diag && !global_txps.contains(&tid) {
                        eprintln!(
                            "[TRUE_TXP] bc={} molecule reassigned: set={:?} -> true={}",
                            bc, global_txps, tid
                        );
                    }
                    global_txps = vec![tid];
                }
            }

            // ── probability vector ──
            let global_txp_prob: Vec<f64> = if !P::HAS_PROBS {
                Vec::new()
            } else if global_txps.len() == 1 {
                vec![1.0]
            } else {
                match prob_model {
                    EMMode::Uniform => {
                        vec![1.0 / global_txps.len() as f64; global_txps.len()]
                    }
                    EMMode::Likelihood => {
                        // Per-transcript log-likelihood, averaged over the members
                        // that carry that transcript. Computed directly rather than
                        // via component_tx_log_scores, because on a CASE2 union no
                        // transcript is present in every member. Averaging (not
                        // summing) keeps transcripts with different member coverage
                        // comparable -- a summed score penalises broad coverage,
                        // since log-likelihoods are negative.
                        let mut tx_scores: Vec<(u32, f64)> =
                            Vec::with_capacity(global_txps.len());
                        for &t in global_txps.iter() {
                            let mut s = 0.0f64;
                            let mut k = 0u32;
                            for &vi in members.iter() {
                                let labels = eqmap.refs_for_eqc(vnodes[vi].eq_id);
                                if let Ok(ti) = labels.binary_search(&t) {
                                    s += vertex_loglik_for_tx(
                                        vnodes[vi].eq_id,
                                        vnodes[vi].umi_idx,
                                        ti,
                                        eqmap,
                                    );
                                    k += 1;
                                }
                            }
                            tx_scores.push((
                                t,
                                if k > 0 { s / k as f64 } else { f64::NEG_INFINITY },
                            ));
                        }

                        let tx_weights = softmax_log_scores(&tx_scores);
                        let mut kept: Vec<(u32, f64)> = tx_weights.to_vec();
                        kept.sort_unstable_by_key(|(t, _)| *t);

                        let sum_p: f64 = kept.iter().map(|(_, p)| *p).sum();
                        if kept.len() == global_txps.len() && sum_p.is_finite() && sum_p > 0.0 {
                            kept.iter().map(|(_, p)| *p / sum_p).collect()
                        } else {
                            if do_diag {
                                eprintln!(
                                    "[PROB_FALLBACK] bc={} case2={} n_txps={} scored={} -> uniform",
                                    bc, is_case2, global_txps.len(), kept.len()
                                );
                            }
                            vec![1.0 / global_txps.len() as f64; global_txps.len()]
                        }
                    }
                }
            };

            // ── per-molecule accounting ──
            let mut parents: Vec<usize> = members.iter().map(|&vi| vnode_parent[vi]).collect();
            parents.sort_unstable();
            parents.dedup();
            let n_parents = parents.len() as u64;

            let has_case1 = parents.iter().any(|&p| vertex_is_impure[p]);
            let has_case2 = is_case2;

            if has_case1 {
                cell_mol_case1 += 1;
                cell_nodes_case1 += n_parents;
            }
            if has_case2 {
                cell_mol_case2 += 1;
                cell_nodes_case2 += n_parents;
            }
            if has_case1 && has_case2 {
                cell_mol_both += 1;
                cell_nodes_both += n_parents;
            }
            if has_case1 || has_case2 {
                cell_mol_either += 1;
                cell_nodes_either += n_parents;
            } else {
                cell_mol_neither += 1;
            }

            cell_molecules += 1;
            cell_reads += num_reads as u64;

            emit(num_reads, global_txps, global_txp_prob, &mut pug_stats);
        }
    }

    ORACLE_CELLS.fetch_add(1, AtomicOrdering::Relaxed);
    ORACLE_VERTICES.fetch_add(cell_vertices, AtomicOrdering::Relaxed);
    ORACLE_VNODES.fetch_add(cell_vnodes, AtomicOrdering::Relaxed);
    ORACLE_MOLECULES.fetch_add(cell_molecules, AtomicOrdering::Relaxed);
    ORACLE_READS.fetch_add(cell_reads, AtomicOrdering::Relaxed);
    ORACLE_CASE1_NODES.fetch_add(cell_case1, AtomicOrdering::Relaxed);
    ORACLE_CASE2_GROUPS.fetch_add(cell_case2, AtomicOrdering::Relaxed);
    ORACLE_MOL_CASE1.fetch_add(cell_mol_case1, AtomicOrdering::Relaxed);
    ORACLE_MOL_CASE2.fetch_add(cell_mol_case2, AtomicOrdering::Relaxed);
    ORACLE_MOL_BOTH.fetch_add(cell_mol_both, AtomicOrdering::Relaxed);
    ORACLE_MOL_EITHER.fetch_add(cell_mol_either, AtomicOrdering::Relaxed);
    ORACLE_MOL_NEITHER.fetch_add(cell_mol_neither, AtomicOrdering::Relaxed);
    ORACLE_NODES_CASE1.fetch_add(cell_nodes_case1, AtomicOrdering::Relaxed);
    ORACLE_NODES_CASE2.fetch_add(cell_nodes_case2, AtomicOrdering::Relaxed);
    ORACLE_NODES_BOTH.fetch_add(cell_nodes_both, AtomicOrdering::Relaxed);
    ORACLE_NODES_EITHER.fetch_add(cell_nodes_either, AtomicOrdering::Relaxed);

    eprintln!(
        "[ORACLE_CELL] bc={} vertices={} vnodes={} molecules={} reads={} case1_nodes={} case2_groups={}",
        bc, cell_vertices, cell_vnodes, cell_molecules, cell_reads, cell_case1, cell_case2
    );
    eprintln!(
        "[ORACLE_CELL_AFFECTED] bc={} mol_case1={} mol_case2={} mol_both={} mol_either={} \
         mol_neither={} nodes_case1={} nodes_case2={} nodes_both={} nodes_either={}",
        bc, cell_mol_case1, cell_mol_case2, cell_mol_both, cell_mol_either,
        cell_mol_neither, cell_nodes_case1, cell_nodes_case2, cell_nodes_both, cell_nodes_either
    );

    pug_stats
}