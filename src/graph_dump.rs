// graph_dump.rs
//
// Dumps the raw PUG graph for a specific cell to TSV files.
// No model parameters needed — all thresholding done in Python.
//
// Usage:
//   DUMP_CELL_NUM=647 target/release/alevin-fry quant ...
//
// Output files:
//   ref_names.tsv          written once: integer tx id -> name
//   cell_647_nodes.tsv     one row per UMI vertex
//   cell_647_edges.tsv     one row per directed edge
//   cell_647_comps.tsv     one row per weakly connected component
//
// Python — filter for any transcript and colour-code the graph:
//
//   import pandas as pd, networkx as nx, matplotlib.pyplot as plt
//
//   txmap = pd.read_csv("ref_names.tsv",        sep="\t")
//   nodes = pd.read_csv("cell_647_nodes.tsv",   sep="\t")
//   edges = pd.read_csv("cell_647_edges.tsv",   sep="\t")
//
//   tx_id = txmap.loc[txmap["tx_name"] == "NM_006013", "tx_id"].iloc[0]
//
//   # vertices that contain the target transcript
//   has_tx = nodes["transcripts"].str.split(",").apply(lambda ts: str(tx_id) in ts)
//
//   # vertices whose best_tx IS the target transcript (chose it)
//   chose_tx = nodes["best_tx"] == tx_id
//
//   # build graph and colour-code:
//   #   red   = contains AND chose target tx
//   #   orange= contains but didn't choose target tx
//   #   grey  = doesn't contain target tx

use std::collections::HashMap;
use std::io::Write;
use petgraph::prelude::*;
use petgraph::visit::NodeIndexable;

use crate::eq_class::EqMap;
use crate::pugutils::{weakly_connected_components, vertex_loglik_for_tx};

// ── public helpers called from quant.rs ───────────────────────────────────

/// Returns the cell number to dump, or None if DUMP_CELL_NUM is not set.
pub fn dump_cell_target() -> Option<usize> {
    std::env::var("DUMP_CELL_NUM")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
}

/// Write ref_names.tsv once from do_quantify before workers start.
/// Maps integer transcript id -> name so Python can resolve any transcript.
/// Only written when DUMP_CELL_NUM is set.
pub fn dump_ref_names(ref_names: &[String]) {
    let path = "ref_names.tsv";
    let mut f = std::fs::File::create(path)
        .expect("could not create ref_names.tsv");
    writeln!(f, "tx_id\ttx_name").unwrap();
    for (i, name) in ref_names.iter().enumerate() {
        writeln!(f, "{}\t{}", i, name).unwrap();
    }
    eprintln!("Wrote {} ({} transcripts)", path, ref_names.len());
}

// ── internal helpers ───────────────────────────────────────────────────────

/// log BF for absorbing vertex (eq_id, umi_idx) into transcript txp.
///   log_BF = LL(txp) - LL(best_alt)
/// NEG_INFINITY when txp absent from eq class.
/// INFINITY     when txp is the only transcript.
fn compute_log_bf(eq_id: u32, umi_idx: u32, txp: u32, eqmap: &EqMap) -> f64 {
    let labels = eqmap.refs_for_eqc(eq_id);
    let tx_index = match labels.binary_search(&txp) {
        Ok(i)  => i,
        Err(_) => return f64::NEG_INFINITY,
    };
    let target_ll = vertex_loglik_for_tx(eq_id, umi_idx, tx_index, eqmap);
    let mut best_alt = f64::NEG_INFINITY;
    for (i, alt) in labels.iter().enumerate() {
        if *alt == txp { continue; }
        let s = vertex_loglik_for_tx(eq_id, umi_idx, i, eqmap);
        if s > best_alt { best_alt = s; }
    }
    if best_alt == f64::NEG_INFINITY { f64::INFINITY } else { target_ll - best_alt }
}

/// Format ±inf as finite sentinels so pandas/R parse cleanly.
fn fmt_f(x: f64) -> String {
    if x == f64::NEG_INFINITY  { "-999.0".to_owned() }
    else if x == f64::INFINITY {  "999.0".to_owned() }
    else                       { format!("{:.6}", x) }
}

/// vertex-index → comp_id lookup.
fn build_v2c(comps: &HashMap<u32, Vec<u32>, ahash::RandomState>) -> HashMap<u32, u32> {
    let mut v2c = HashMap::new();
    for (comp_id, verts) in comps.iter() {
        for v in verts.iter() { v2c.insert(*v, *comp_id); }
    }
    v2c
}

// ── main entry point ───────────────────────────────────────────────────────

/// Dump the raw PUG for `cell_num` to three TSV files.
/// Call from run_worker_thread when cell_num == dump_cell_target().
/// No model parameters needed.
pub fn dump_cell_graph(
    cell_num: usize,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    eqmap: &EqMap,
) {
    let prefix    = format!("cell_{}", cell_num);
    let num_nodes = g.node_count();
    let comps     = weakly_connected_components(g);
    let v2c       = build_v2c(&comps);

    // ── nodes ──────────────────────────────────────────────────────────────
    // Columns:
    //   vertex_id     integer graph index
    //   comp_id       weakly connected component
    //   eq_id         equivalence class id
    //   umi_idx       UMI index within eq class
    //   n_reads       reads at this UMI
    //   n_transcripts transcripts in this eq class
    //   transcripts   comma-separated tx ids sorted by LL descending
    //   best_tx       tx id with highest LL  ← the one this vertex "chose"
    //   best_ll       its log-likelihood
    //   second_ll     second-best LL  (-999 sentinel if only one tx)
    //   delta         best_ll - second_ll  (999 sentinel if only one tx)
    //   mean_ll_best  mean per-read LL for best transcript
    //
    // Python colour coding:
    //   has_tx  = transcripts column contains target tx id
    //   chose   = best_tx == target tx id
    //   colour:   red   = chose   (node strongly assigned to target)
    //             orange= has_tx but not chose  (ambiguous, target competing)
    //             grey  = no target tx at all

    let node_path = format!("{}_nodes.tsv", prefix);
    let mut nf = std::fs::File::create(&node_path).expect("could not create node file");
    writeln!(nf,
        "vertex_id\tcomp_id\teq_id\tumi_idx\tn_reads\tn_transcripts\t\
         transcripts\tbest_tx\tbest_ll\tsecond_ll\tdelta\tmean_ll_best"
    ).unwrap();

    for v in 0..num_nodes {
        let (eq_id, umi_idx) = g.from_index(v);
        let labels = eqmap.refs_for_eqc(eq_id);

        let mut tx_lls: Vec<(u32, f64)> = labels.iter()
            .enumerate()
            .map(|(i, t)| (*t, vertex_loglik_for_tx(eq_id, umi_idx, i, eqmap)))
            .collect();
        tx_lls.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        let best_tx   = tx_lls.first().map(|(t, _)| *t).unwrap_or(0);
        let best_ll   = tx_lls.first().map(|(_, l)| *l).unwrap_or(f64::NEG_INFINITY);
        let second_ll = tx_lls.get(1).map(|(_, l)| *l).unwrap_or(f64::NEG_INFINITY);
        let delta     = if second_ll == f64::NEG_INFINITY { f64::INFINITY }
                        else { best_ll - second_ll };

        let best_tx_index = labels.binary_search(&best_tx).unwrap_or(0);
        let (n_reads, mean_ll_best) =
            match eqmap.probs_for_eq_umi_tx(eq_id, umi_idx, best_tx_index) {
                Some(pv) => {
                    let n = pv.len();
                    let m = if n > 0 {
                        pv.iter().map(|&p| p.max(1e-12_f64).ln()).sum::<f64>() / n as f64
                    } else { 0.0 };
                    (n, m)
                }
                None => (0, 0.0),
            };

        let tx_str  = tx_lls.iter().map(|(t, _)| t.to_string()).collect::<Vec<_>>().join(",");
        let comp_id = v2c.get(&(v as u32)).cloned().unwrap_or(u32::MAX);

        writeln!(nf,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{:.6}",
            v, comp_id, eq_id, umi_idx, n_reads, labels.len(),
            tx_str, best_tx, best_ll,
            fmt_f(second_ll), fmt_f(delta), mean_ll_best
        ).unwrap();
    }
    eprintln!("Wrote {}", node_path);

    // ── edges ──────────────────────────────────────────────────────────────
    // Columns:
    //   from_vertex   source vertex id
    //   to_vertex     target vertex id
    //   same_comp     1 if same component
    //   share_tx      1 if they share at least one transcript
    //   bf_from_to    log BF for absorbing to_vertex using from_vertex's best_tx
    //   bf_to_from    log BF for absorbing from_vertex using to_vertex's best_tx
    //
    // Python: apply any threshold — e.g. edges[edges["bf_from_to"] >= 0]

    let edge_path = format!("{}_edges.tsv", prefix);
    let mut ef = std::fs::File::create(&edge_path).expect("could not create edge file");
    writeln!(ef,
        "from_vertex\tto_vertex\tsame_comp\tshare_tx\tbf_from_to\tbf_to_from"
    ).unwrap();

    for v in 0..num_nodes {
        let (eq_id_v, umi_idx_v) = g.from_index(v);
        let labels_v = eqmap.refs_for_eqc(eq_id_v);
        let comp_v   = v2c.get(&(v as u32)).cloned().unwrap_or(u32::MAX);

        // best transcript of this source vertex
        let best_tx_v = labels_v.iter()
            .enumerate()
            .map(|(i, t)| (*t, vertex_loglik_for_tx(eq_id_v, umi_idx_v, i, eqmap)))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(t, _)| t)
            .unwrap_or(0);

        for nv in g.neighbors_directed(g.from_index(v), Outgoing) {
            let n                    = g.to_index(nv) as usize;
            let (eq_id_n, umi_idx_n) = nv;
            let labels_n             = eqmap.refs_for_eqc(eq_id_n);
            let comp_n               = v2c.get(&(n as u32)).cloned().unwrap_or(u32::MAX);

            let same_comp = (comp_v == comp_n) as u8;
            let share_tx  = labels_v.iter().any(|t| labels_n.binary_search(t).is_ok()) as u8;

            let bf_from_to = compute_log_bf(eq_id_n, umi_idx_n, best_tx_v, eqmap);

            let best_tx_n = labels_n.iter()
                .enumerate()
                .map(|(i, t)| (*t, vertex_loglik_for_tx(eq_id_n, umi_idx_n, i, eqmap)))
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(t, _)| t)
                .unwrap_or(0);

            let bf_to_from = compute_log_bf(eq_id_v, umi_idx_v, best_tx_n, eqmap);

            writeln!(ef, "{}\t{}\t{}\t{}\t{}\t{}",
                v, n, same_comp, share_tx,
                fmt_f(bf_from_to), fmt_f(bf_to_from)
            ).unwrap();
        }
    }
    eprintln!("Wrote {}", edge_path);

    // ── components ─────────────────────────────────────────────────────────
    // Columns:
    //   comp_id       component label
    //   size          number of vertices
    //   vertices      comma-separated vertex ids
    //   n_shared_tx   transcripts shared by ALL vertices in component
    //   ambiguous     1 if vertices disagree on best_tx

    let comp_path = format!("{}_comps.tsv", prefix);
    let mut cf = std::fs::File::create(&comp_path).expect("could not create comp file");
    writeln!(cf, "comp_id\tsize\tvertices\tn_shared_tx\tambiguous").unwrap();

    for (comp_id, verts) in comps.iter() {
        let vert_str = verts.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",");
        let mut shared: Option<std::collections::HashSet<u32>> = None;
        let mut best_txs: Vec<u32> = Vec::new();

        for &v in verts.iter() {
            let (eq_id, umi_idx) = g.from_index(v as usize);
            let labels = eqmap.refs_for_eqc(eq_id);
            let tx_set: std::collections::HashSet<u32> = labels.iter().cloned().collect();
            shared = Some(match shared {
                None       => tx_set,
                Some(prev) => prev.intersection(&tx_set).cloned().collect(),
            });
            let best = labels.iter()
                .enumerate()
                .map(|(i, t)| (*t, vertex_loglik_for_tx(eq_id, umi_idx, i, eqmap)))
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(t, _)| t)
                .unwrap_or(0);
            best_txs.push(best);
        }

        let n_shared  = shared.map(|s| s.len()).unwrap_or(0);
        let ambiguous = {
            let first = best_txs.first().cloned().unwrap_or(0);
            best_txs.iter().any(|t| *t != first) as u8
        };

        writeln!(cf, "{}\t{}\t{}\t{}\t{}",
            comp_id, verts.len(), vert_str, n_shared, ambiguous
        ).unwrap();
    }
    eprintln!("Wrote {}", comp_path);

    eprintln!("=== dump done cell={} vertices={} edges={} comps={} ===",
        cell_num, num_nodes, g.edge_count(), comps.len());
}