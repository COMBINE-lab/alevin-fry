"""Step 1: remove the diagnostic / instrumentation code from the alevin-fry copy.

Removed: the molecule dump (MOLDUMP_*), the UMI-graph dump (DUMP_CELL_NUM), the coverage dump
(COVERAGE_DIAG_*), the ground-truth UMI oracle modes (--true-umi-mode edge/bfs/direct,
--true-txp-mode, --true-umi-txp-file, the TRUE_UMI_ORACLE totals), the per-cell read-name map and
the read names kept per equivalence class, and the cell 647 debug breakpoint.
Every replacement asserts the text it replaces occurs exactly as expected.
"""
import re, os, sys

ROOT = "/scratch1/zahra/source-code/alevin-fry/src"


def edit(path, pairs):
    p = os.path.join(ROOT, path)
    s = open(p).read()
    for old, new, count in pairs:
        n = s.count(old)
        assert n == count, f"{path}: expected {count} occurrence(s) of {old[:70]!r}, found {n}"
        s = s.replace(old, new)
    open(p, "w").write(s)


def cut(path, start, end, replacement="", keep_end=True):
    """Remove the text from `start` (inclusive) to `end` (exclusive when keep_end) — exactly one match."""
    p = os.path.join(ROOT, path)
    s = open(p).read()
    i = s.find(start); assert i >= 0 and s.find(start, i + 1) < 0, f"{path}: start marker {start[:60]!r} not unique"
    j = s.find(end, i); assert j >= 0, f"{path}: end marker {end[:60]!r} not found after start"
    if not keep_end:
        j += len(end)
    s = s[:i] + replacement + s[j:]
    open(p, "w").write(s)


# ---------------------------------------------------------------- quant.rs
edit("quant.rs", [
    ("use crate::coverage_diag;\nuse crate::graph_dump;\nuse crate::pugutils_dp_new::{\n    DpConfig, ResolvedMolecule, get_num_molecules_log_dp_new, get_num_molecules_true_umi_oracle,\n    get_num_molecules_true_umi_oracle_virtual_node,\n};\nuse crate::truncation;\nuse crate::truncation_groups;\nuse crate::molecule_dump;\n",
     "use crate::pugutils_dp_new::{DpConfig, ResolvedMolecule, get_num_molecules_log_dp_new};\nuse crate::truncation;\nuse crate::truncation_groups;\n", 1),
    ("LongReadEqClassPayload, MaybeQName,", "LongReadEqClassPayload,", 1),
    ("    true_umi_mode: TrueUmiOracleMode,\n    true_txp_mode: Option<bool>,\n", "", 1),
    ("    true_umi_lookup: Option<Arc<HashMap<String, Vec<u8>>>>,\n    true_txp_lookup: Option<Arc<HashMap<String, String>>>,\n", "", 1),
    ("        + MaybeQName\n", "", 4),
    ("                let mut cell_num_ok = false;\n                if cell_num == 647 {\n                    info!(log, \"debug breakpoint cell_num = 647\");\n                    cell_num_ok = true;\n                }\n\n", "", 1),
    ("                                        &log,\n                                        cell_num_ok,\n                                    );", "                                        &log,\n                                    );", 1),
    ("                let dump_cell = molecule_dump::enabled_for(cell_num);\n                let mut resolved_molecules = (config.coverage_model\n                    || config.coverage_dedup_model\n                    || config.truncation_model\n                    || dump_cell)\n",
     "                let mut resolved_molecules = (config.coverage_model\n                    || config.coverage_dedup_model\n                    || config.truncation_model)\n", 1),
    ("                                //eprintln!(\"before the init from chunk\");\n\n                                // step 1: collect umi -> qnames BEFORE init_from_chunk destroys ordering\n                                let mut refs_umi_to_qnames: HashMap<(Vec<u32>, u64), Vec<String>> =\n                                    HashMap::new();\n                                for r in c.reads.iter() {\n                                    if let Some(qname) = r.maybe_qname() {\n                                        refs_umi_to_qnames\n                                            .entry((r.refs().to_vec(), r.umi()))\n                                            .or_default()\n                                            .push(qname.to_string());\n                                    }\n                                }\n\n", "", 1),
    ("                            if let Some(target) = graph_dump::dump_cell_target() {\n                                if cell_num == target {\n                                    graph_dump::dump_cell_graph(cell_num, &g, &eq_map);\n                                }\n                            }\n\n", "", 1),
    ("                                        if let Some(diag) = coverage_diag::get() {\n                                            diag.dump_cell(\n                                                cell_num,\n                                                bc.into(),\n                                                config.barcode_len,\n                                                molecules,\n                                                coverage,\n                                                &shared.tid_to_gid,\n                                            );\n                                        }\n", "", 1),
    ("                            if let (Some(priors), Some(molecules), Some(d)) = (\n                                dump_priors.as_ref(),\n                                resolved_molecules.as_ref(),\n                                molecule_dump::get(),\n                            ) {\n                                d.dump_molecules(cell_num, bc.into(), config.barcode_len, molecules, priors, dump_logpos.as_deref());\n                            }\n\n", "", 1),
    ("    let true_umi_mode = quant_opts.true_umi_mode;\n    let true_txp_mode = quant_opts.true_txp_mode;\n    let true_umi_txp_file = quant_opts.true_umi_txp_file.clone();\n", "", 1),
    ("    if molecule_dump::get().is_some() {\n        molecule_dump::write_ref_names(&hdr.ref_names);\n    }\n    if graph_dump::dump_cell_target().is_some() {\n        graph_dump::dump_ref_names(&hdr.ref_names);\n    }\n\n", "", 1),
    ("        let true_umi_lookup = true_umi_lookup.clone();\n        let true_txp_lookup = true_txp_lookup.clone();\n", "", 1),
    ("            true_umi_mode,\n            true_txp_mode,\n", "", 1),
    ("            true_umi_lookup,\n            true_txp_lookup,\n", "", 1),
    ("    crate::pugutils_dp_new::report_true_umi_oracle_totals();\n\n", "", 1),
])
# the oracle enum, its loaders, the per-cell oracle tables and the dump priors
cut("quant.rs", "#[derive(PartialEq, Eq, Debug, Default, Clone, Copy, Serialize)]\npub enum TrueUmiOracleMode {",
    "#[derive(PartialEq, Eq, Debug, Default, Clone, Copy, Serialize)]\npub enum EditDistanceMode {")
cut("quant.rs", "/// Loads read_name -> true_umi (as raw bytes) from the ground-truth TSV.\n", "// TODO: see if we'd rather pass an structure\n")
cut("quant.rs", "                                //eprintln!(\"after the init from chunk\");\n", "                                //let coverage_bin_width: u32 = 100;\n")
cut("quant.rs", "                            let true_umi_table = shared\n", "                            // for the PUG resolution algorithm, set the hasher\n",
    replacement="""                            let (g, pair_posteriors): (_, Option<pugutils::PairPosteriors>) = match config.dedup_edges {
                                DedupEdges::Bayes | DedupEdges::BayesTx => {
                                    let (g, st, post) = pugutils::extract_graph_bayes(
                                        &eq_map,
                                        config.ed_distance_model,
                                        config.ed_far_window,
                                        if config.dedup_count == DedupCount::Expected { config.dedup_soft_floor } else { 0.5 },
                                        config.dedup_edges == DedupEdges::BayesTx,
                                        &log,
                                    );
                                    if cell_num % 250 == 0 {
                                        info!(log, "cell_num = {}: UMI edge model {:?}", cell_num, st);
                                    }
                                    (g, Some(post))
                                }
                                DedupEdges::Rule => (
                                    pugutils::extract_graph(
                                        &eq_map,
                                        config.pug_exact_umi,
                                        config.ed_distance_model,
                                        config.ed_upper_bound,
                                        far_edges(config.ed_far_bound, config.ed_upper_bound, config.ed_far_window),
                                        &log,
                                    ),
                                    None,
                                ),
                            };

""")
cut("quant.rs", "                            let pug_stats = match config.true_umi_mode {\n",
    "                            alt_resolution = pug_stats.used_alternative_strategy; // alt_res;\n",
    replacement="""                            let dp_cfg = DpConfig {
                                log_prior_odds: config.lambda_size,
                                tau_delta: config.tau_delta,
                            };
                            let pug_stats = match config.dedup_model {
                                DedupModel::Bf => get_num_molecules_log_dp_new(
                                    &g,
                                    &eq_map,
                                    &shared.tid_to_gid,
                                    &mut gene_eqc,
                                    &s,
                                    config.large_graph_thresh,
                                    dp_cfg,
                                    resolved_molecules.as_mut(),
                                    &log,
                                ),
                                DedupModel::Parsimony | DedupModel::ParsimonyWeighted => {
                                    pugutils::get_num_molecules::<P>(
                                        &g,
                                        &eq_map,
                                        &shared.tid_to_gid,
                                        &mut gene_eqc,
                                        &s,
                                        config.large_graph_thresh,
                                        config.dedup_model == DedupModel::ParsimonyWeighted,
                                        if config.dedup_cover == DedupCover::Posterior || config.dedup_count == DedupCount::Expected { pair_posteriors.as_ref() } else { None },
                                        config.dedup_count == DedupCount::Expected,
                                        resolved_molecules.as_mut(),
                                        &log,
                                    )
                                }
                            };

                            alt_resolution = pug_stats.used_alternative_strategy;
""", keep_end=False)
cut("quant.rs", "                            let dump_priors: Option<Vec<Vec<(u32, f64)>>> = if dump_cell {\n",
    "                            // The coverage the positional models use for this cell: its own, from\n")
cut("quant.rs", "    let true_umi_lookup: Option<Arc<HashMap<String, Vec<u8>>>> = match quant_opts.true_umi_mode {\n",
    "    // create our output directory\n")

# ---------------------------------------------------------------- main.rs
edit("main.rs", [
    ("    DedupCount, DedupCover, DedupEdges, DedupModel, EMMode, EditDistanceMode, EndDistribution, EndModel, ResolutionStrategy, SplicedAmbiguityModel,\n    TrueUmiOracleMode,\n};",
     "    DedupCount, DedupCover, DedupEdges, DedupModel, EMMode, EditDistanceMode, EndDistribution, EndModel, ResolutionStrategy, SplicedAmbiguityModel,\n};", 1),
    ("""    .arg(arg!(--"true-umi-mode" <MODE> "how (if at all) to use ground-truth UMIs to compute an oracle upper bound")
        .ignore_case(true)
        .value_parser(value_parser!(TrueUmiOracleMode))
        .default_value("off"))
    .arg(arg!(--"true-txp-mode" <TxpMODE> "if using the true transcript in the oracle mode")
        .ignore_case(true)
        .value_parser(value_parser!(bool)))
    .arg(arg!(--"true-umi-txp-file" <TRUEUMIFILE> "TSV file mapping read names to true UMI/gene ground truth; required unless --true-umi-mode=off")
        .value_parser(pathbuf_file_exists_validator))
""",
     """    .arg(arg!(--"true-umi-mode" <MODE> "removed; accepted for compatibility with older command lines and must be `off`")
        .hide(true)
        .default_value("off"))
""", 1),
    ("""        let true_umi_mode = *t.get_one::<TrueUmiOracleMode>("true-umi-mode").unwrap();
        let true_txp_mode: Option<bool> = t.get_one("true-txp-mode").cloned();
        let true_umi_txp_file: Option<&PathBuf> = t.get_one("true-umi-txp-file");
""",
     """        if t.get_one::<String>("true-umi-mode").is_some_and(|m| !m.eq_ignore_ascii_case("off")) {
            bail!("--true-umi-mode: the ground-truth oracle modes were removed; only `off` is accepted");
        }
""", 1),
    ("        if junction_model && true_umi_mode == TrueUmiOracleMode::DirectCollapse {\n            bail!(\"--junction-model is not supported with --true-umi-mode=direct\");\n        }\n", "", 1),
    ("        if truncation_model && true_umi_mode == TrueUmiOracleMode::DirectCollapse {\n            bail!(\"--model-truncation is not supported with --true-umi-mode=direct\");\n        }\n", "", 1),
    ("        if coverage_enabled && true_umi_mode == TrueUmiOracleMode::DirectCollapse {\n            bail!(\"coverage modes are not supported with --true-umi-mode=direct\");\n        }\n", "", 1),
    ("""        if (true_umi_mode != TrueUmiOracleMode::Off || true_txp_mode == Some(true))
            && true_umi_txp_file.is_none()
        {
            crit!(
                log,
                "--true-umi-mode was set to {:?} but no --true-umi-txp-file was provided",
                true_umi_mode
            );
            bail!("Invalid command line option");
        }
""", "", 1),
    ("            .true_umi_mode(true_umi_mode)\n            .true_txp_mode(true_txp_mode)\n            .true_umi_txp_file(true_umi_txp_file.cloned())\n", "", 1),
])

# ---------------------------------------------------------------- prog_opts.rs
edit("prog_opts.rs", [
    ("    DedupCount, DedupCover, DedupEdges, DedupModel, EMMode, EditDistanceMode, EndDistribution, EndModel, ResolutionStrategy, SplicedAmbiguityModel,\n    TrueUmiOracleMode,\n};",
     "    DedupCount, DedupCover, DedupEdges, DedupModel, EMMode, EditDistanceMode, EndDistribution, EndModel, ResolutionStrategy, SplicedAmbiguityModel,\n};", 1),
    ("    pub true_umi_mode: TrueUmiOracleMode,\n    pub true_txp_mode: Option<bool>,\n    pub true_umi_txp_file: Option<PathBuf>,\n", "", 1),
])

# ---------------------------------------------------------------- pugutils_dp_new.rs: everything after the BF resolver is oracle code
p = os.path.join(ROOT, "pugutils_dp_new.rs"); s = open(p).read()
marker = "\n\n\n\n\n\npub fn get_num_molecules_true_umi_oracle<P: EqClassPayload>("
i = s.find(marker); assert i > 0, "oracle start not found"
s = s[:i] + "\n"
s = s.replace("use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};\n", "", 1)
s2 = re.sub(r"static ORACLE_[A-Z0-9_]+: AtomicU64 = AtomicU64::new\(0\);\n", "", s)
assert s2.count("ORACLE_") == 0, "oracle statics remain"
open(p, "w").write(s2)

# ---------------------------------------------------------------- eq_class.rs
edit("eq_class.rs", [
    ("    pub qnames: Vec<Vec<String>>,\n", "", 1),
    ("                        qnames: vec![],\n", "", 4),
    ("""    pub fn qnames_for_eq_umi(&self, eqid: u32, umi_rank: usize) -> &[String] {
        self.eqc_info
            .get(eqid as usize)
            .and_then(|e| e.qnames.get(umi_rank))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

""", "", 1),
])

# ---------------------------------------------------------------- utils.rs: the read-name trait
cut("utils.rs", "// Helper trait to optionally extract qname\n", "/// Trait that allows us to paramaterize the payload stored along with\n")

# ---------------------------------------------------------------- lib.rs and the files
edit("lib.rs", [("pub mod graph_dump;\npub mod coverage_diag;\n", "", 1), ("pub mod molecule_dump;\n", "", 1)])
for f in ["graph_dump.rs", "coverage_diag.rs", "molecule_dump.rs"]:
    os.remove(os.path.join(ROOT, f))
print("step 1 edits applied")
