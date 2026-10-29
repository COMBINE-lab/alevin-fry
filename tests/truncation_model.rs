use std::collections::HashMap;

use alevin_fry::eq_class::AlignmentInterval;
use alevin_fry::pugutils_dp_new::ResolvedMolecule;
use alevin_fry::truncation::{
    cell_molecules, fit_global, fit_molecules, fit_molecules_with, FitExtras, FreezeReason,
    TruncationAbundance, TruncationEndKernels, TruncationIntervalMode, TruncationModelSpec,
    TruncationPositionalLaw,
};

/// The oarfish model as before: per-cell parameters, mixed end kernels, full positional power.
fn spec() -> TruncationModelSpec {
    TruncationModelSpec::default()
}

fn molecule(target_intervals: &[(u32, f64, u32, u32, u32)]) -> ResolvedMolecule {
    let mut transcript_probs = Vec::new();
    let mut intervals = HashMap::new();
    for &(target, probability, start, end, len) in target_intervals {
        if !intervals.contains_key(&target) {
            transcript_probs.push((target, probability));
        }
        intervals
            .entry(target)
            .or_insert_with(Vec::new)
            .push(AlignmentInterval {
                transcript_id: target,
                start,
                end,
                transcript_len: len,
            });
    }
    ResolvedMolecule {
        transcript_probs,
        intervals,
        vertices: Vec::new(),
    }
}

#[test]
fn unambiguous_molecules_keep_their_counts() {
    let mut molecules = (0..20)
        .map(|_| molecule(&[(0, 1.0, 0, 1000, 1000)]))
        .chain((0..10).map(|_| molecule(&[(1, 1.0, 0, 500, 500)])))
        .collect::<Vec<_>>();

    let fit = fit_molecules(&mut molecules, TruncationIntervalMode::All, spec());

    assert!(fit.converged);
    assert!(fit.iterations < 1000);
    assert!(fit.has_positional_evidence);
    assert_eq!(fit.molecule_count, 30);
    assert_eq!(molecules.len(), 30);
    assert!(molecules[..20]
        .iter()
        .all(|m| m.transcript_probs == [(0, 1.0)]));
    assert!(molecules[20..]
        .iter()
        .all(|m| m.transcript_probs == [(1, 1.0)]));
}

#[test]
fn empty_input_reports_no_convergence_or_iterations() {
    let mut molecules = Vec::new();
    let fit = fit_molecules(&mut molecules, TruncationIntervalMode::All, spec());

    assert!(!fit.converged);
    assert_eq!(fit.iterations, 0);
    assert_eq!(fit.molecule_count, 0);
    assert!(!fit.has_positional_evidence);
}

#[test]
fn stable_fit_reports_convergence_without_positional_evidence() {
    let mut molecules = (0..5)
        .map(|_| ResolvedMolecule {
            transcript_probs: vec![(0, 1.0)],
            intervals: HashMap::new(),
            vertices: Vec::new(),
        })
        .collect::<Vec<_>>();
    let fit = fit_molecules(&mut molecules, TruncationIntervalMode::All, spec());

    assert!(fit.converged);
    assert!(fit.iterations <= 3);
    assert_eq!(fit.molecule_count, 5);
    assert!(!fit.has_positional_evidence);
}

#[test]
fn full_length_evidence_shifts_shared_molecules_to_the_short_isoform() {
    let mut molecules = Vec::new();
    molecules.extend((0..80).map(|_| molecule(&[(0, 1.0, 0, 2000, 2000)])));
    molecules.extend((0..120).map(|i| molecule(&[(0, 1.0, 400 + i, 2000, 2000)])));
    molecules
        .extend((0..300).map(|_| molecule(&[(0, 0.5, 1000, 2000, 2000), (1, 0.5, 0, 1000, 1000)])));

    fit_molecules(&mut molecules, TruncationIntervalMode::All, spec());

    let short_isoform_mass: f64 = molecules[200..]
        .iter()
        .map(|m| {
            m.transcript_probs
                .iter()
                .find_map(|(target, probability)| (*target == 1).then_some(*probability))
                .unwrap_or(0.0)
        })
        .sum();
    assert!(
        short_isoform_mass > 150.0,
        "expected positional evidence to favor short isoform; mass={short_isoform_mass}"
    );
}

#[test]
fn union_mode_uses_one_outer_interval_per_molecule_and_transcript() {
    let make_molecules = || {
        let mut molecules = (0..80)
            .map(|_| molecule(&[(0, 1.0, 0, 2000, 2000)]))
            .collect::<Vec<_>>();
        molecules.extend((0..80).map(|i| molecule(&[(0, 1.0, 500 + i, 2000, 2000)])));
        molecules.extend((0..80).map(|_| molecule(&[(1, 1.0, 0, 1000, 1000)])));
        molecules.extend((0..40).map(|_| {
            molecule(&[
                (0, 0.5, 0, 400, 2000),
                (0, 0.5, 1500, 2000, 2000),
                (1, 0.5, 0, 1000, 1000),
            ])
        }));
        molecules
    };

    let mut all_intervals = make_molecules();
    let mut union_interval = make_molecules();

    fit_molecules(&mut all_intervals, TruncationIntervalMode::All, spec());
    fit_molecules(&mut union_interval, TruncationIntervalMode::Union, spec());

    let ambiguous_index = 240;
    let all_probability = all_intervals[ambiguous_index]
        .transcript_probs
        .iter()
        .find_map(|(target, probability)| (*target == 0).then_some(*probability))
        .unwrap();
    let union_probability = union_interval[ambiguous_index]
        .transcript_probs
        .iter()
        .find_map(|(target, probability)| (*target == 0).then_some(*probability))
        .unwrap();
    assert!(
        (all_probability - union_probability).abs() > 1e-4,
        "the interval modes should use different positional evidence"
    );
}

#[test]
fn three_prime_anchored_reads_are_learned_per_cell() {
    let mut molecules = (0..300)
        .map(|i| molecule(&[(0, 1.0, 200 + (i % 7) * 100, 2000, 2000)]))
        .collect::<Vec<_>>();
    let fit = fit_molecules(&mut molecules, TruncationIntervalMode::All, spec());

    assert!(fit.has_positional_evidence);
    if let Some(frozen) = fit.frozen_at {
        assert!(frozen >= 3 && frozen <= fit.iterations);
    }
    assert_eq!(fit.frozen_at.is_some(), fit.freeze_reason.is_some());
    assert!(fit.params.p3 > fit.params.p5);
    assert!(fit.params.k3_gap[0] > 0.5, "{:?}", fit.params);
    let total = fit.params.p3 + fit.params.p5 + fit.params.pf + fit.params.pu;
    assert!((total - 1.0).abs() < 1e-9);
}

#[test]
fn parameters_freeze_when_abundances_keep_changing() {
    // Nearly all molecules are ambiguous between two identical-looking transcripts, so the
    // abundances converge slowly.
    let mut molecules = Vec::new();
    molecules.extend((0..30).map(|i| molecule(&[(0, 1.0, i % 5, 1000, 1000)])));
    molecules.extend((0..10).map(|i| molecule(&[(1, 1.0, i % 5, 1000, 1000)])));
    molecules.extend(
        (0..2000).map(|i| molecule(&[(0, 0.5, i % 5, 1000, 1000), (1, 0.5, i % 5, 1000, 1000)])),
    );
    let fit = fit_molecules(&mut molecules, TruncationIntervalMode::All, spec());

    let frozen = fit.frozen_at.expect("parameters should freeze");
    assert!(frozen >= 4 && frozen <= fit.iterations);
    assert!(fit.freeze_reason.is_some());
}

/// Cells of 3'-anchored molecules truncated at varying lengths, some ambiguous between a long
/// isoform and the short one that is its 3' end; cell `c` has its own isoform proportions.
fn anchored_cells(n: u32) -> Vec<Vec<ResolvedMolecule>> {
    (0..n)
        .map(|c| {
            let mut molecules: Vec<_> = (0..200 + 20 * c)
                .map(|i| molecule(&[(0, 1.0, 200 + (i % 7) * 100, 2000, 2000)]))
                .collect();
            molecules.extend(
                (0..50 + 10 * c)
                    .map(|i| molecule(&[(0, 0.5, 1500 + i % 5, 2000, 2000), (1, 0.5, 500 + i % 5, 1000, 1000)])),
            );
            molecules
        })
        .collect()
}

#[test]
fn global_parameters_are_learned_from_all_cells() {
    let mut cells: Vec<_> = anchored_cells(4)
        .iter()
        .map(|m| cell_molecules(m, TruncationIntervalMode::All, 0, None))
        .collect();
    let fit = fit_global(&mut cells, 0, 1.0, 2, spec(), TruncationAbundance::Cell, 2, 1, FitExtras::default());

    assert_eq!(fit.cells, 4);
    assert_eq!(fit.molecules, cells.iter().map(|c| c.molecule_count()).sum::<usize>());
    let frozen = fit.frozen_at.expect("parameters should freeze");
    assert!(frozen >= 3 && frozen == fit.iterations);
    assert!(fit.freeze_reason.is_some());
    assert!(fit.params.p3 > fit.params.p5, "{:?}", fit.params);
    assert!(fit.params.k3_gap[0] > 0.5, "{:?}", fit.params);
}

#[test]
fn global_fit_does_not_depend_on_the_number_of_threads() {
    let mut cells: Vec<_> = anchored_cells(5)
        .iter()
        .map(|m| cell_molecules(m, TruncationIntervalMode::All, 0, None))
        .collect();
    let (a, b) = (
        fit_global(&mut cells, 0, 1.0, 1, spec(), TruncationAbundance::Cell, 2, 1, FitExtras::default()),
        fit_global(&mut cells, 0, 1.0, 3, spec(), TruncationAbundance::Cell, 2, 1, FitExtras::default()),
    );
    assert_eq!(a.iterations, b.iterations);
    for (x, y) in [
        (a.params.p3, b.params.p3),
        (a.params.pf, b.params.pf),
        (a.params.len3_mean, b.params.len3_mean),
        (a.params.k3_gap[0], b.params.k3_gap[0]),
    ] {
        assert!((x - y).abs() <= 1e-9 * x.abs().max(1.0), "{:?} {:?}", a.params, b.params);
    }
}

#[test]
fn cells_use_the_global_parameters_unchanged() {
    let mut cells = anchored_cells(3);
    let global = fit_global(
        &mut cells
            .iter()
            .map(|m| cell_molecules(m, TruncationIntervalMode::All, 0, None))
            .collect::<Vec<_>>(),
        0,
        1.0,
        1,
        spec(),
        TruncationAbundance::Cell,
        2,
        1,
        FitExtras::default(),
    );
    let molecules = &mut cells[1];
    let fit = fit_molecules_with(molecules, TruncationIntervalMode::All, Some(&global), spec(), None, 0, None);

    assert_eq!(fit.frozen_at, Some(0));
    assert_eq!(fit.freeze_reason, Some(FreezeReason::Global));
    assert!(fit.converged);
    assert!(fit.has_positional_evidence);
    assert_eq!(fit.params.p3, global.params.p3);
    assert_eq!(fit.params.k3_gap, global.params.k3_gap);
    // Unambiguous molecules keep their transcript; every molecule's probabilities sum to one.
    assert!(molecules[..220].iter().all(|m| m.transcript_probs == [(0, 1.0)]));
    for m in molecules.iter() {
        let total: f64 = m.transcript_probs.iter().map(|(_, p)| p).sum();
        assert!((total - 1.0).abs() < 1e-9);
    }
}

#[test]
fn a_sample_freezes_no_earlier_than_the_same_cells_as_all_the_data() {
    let mut cells: Vec<_> = anchored_cells(3)
        .iter()
        .map(|m| cell_molecules(m, TruncationIntervalMode::All, 0, None))
        .collect();
    let (all, sample) = (
        fit_global(&mut cells, 0, 1.0, 1, spec(), TruncationAbundance::Cell, 2, 1, FitExtras::default()),
        fit_global(&mut cells, 0, 0.01, 1, spec(), TruncationAbundance::Cell, 2, 1, FitExtras::default()),
    );
    assert!(sample.iterations >= all.iterations);
}

#[test]
fn zero_probability_candidates_are_dropped_not_fatal() {
    // A softmax over very different log-likelihoods underflows to exactly zero.
    let mut molecules = (0..20)
        .map(|_| molecule(&[(0, 1.0, 0, 1000, 1000), (1, 0.0, 0, 500, 500)]))
        .collect::<Vec<_>>();
    let fit = fit_molecules(&mut molecules, TruncationIntervalMode::All, spec());

    assert_eq!(fit.molecule_count, 20);
    assert!(molecules.iter().all(|m| m.transcript_probs == [(0, 1.0)]));
}

/// Model v2: pooled positional parameters, per-cell abundances shrunk toward the pooled isoform
/// shares, absolute end kernels, full positional power, and each molecule drawn whole.
fn v2_spec(sample_below: f64, min_share: f64) -> TruncationModelSpec {
    TruncationModelSpec {
        end_kernels: TruncationEndKernels::Absolute,
        sample_below,
        min_share,
        ..TruncationModelSpec::default()
    }
}

/// `anchored_cells` with a distinct (equivalence class, UMI) vertex per molecule, as deduplication
/// gives them, so that each molecule's draw is its own.
fn anchored_cells_with_vertices(n: u32) -> Vec<Vec<ResolvedMolecule>> {
    let mut umi = 0u64;
    anchored_cells(n)
        .into_iter()
        .map(|mut cell| {
            for m in cell.iter_mut() {
                umi += 1;
                m.vertices = vec![((umi % 97) as u32, umi * 2_654_435_761)];
            }
            cell
        })
        .collect()
}

/// Fits model v2 to `cells` and assigns every cell's molecules; returns the cells.
fn fit_v2(mut cells: Vec<Vec<ResolvedMolecule>>, s: TruncationModelSpec) -> Vec<Vec<ResolvedMolecule>> {
    let mut input: Vec<_> = cells
        .iter()
        .map(|m| cell_molecules(m, TruncationIntervalMode::All, 0, None))
        .collect();
    let global = fit_global(&mut input, 0, 1.0, 2, s, TruncationAbundance::Shrink, 2, 1, FitExtras::default());
    for molecules in cells.iter_mut() {
        fit_molecules_with(molecules, TruncationIntervalMode::All, Some(&global), s, Some(1.0), 0, None);
    }
    cells
}

fn share(m: &ResolvedMolecule, target: u32) -> f64 {
    m.transcript_probs
        .iter()
        .find_map(|(t, p)| (*t == target).then_some(*p))
        .unwrap_or(0.0)
}

#[test]
fn drawing_whole_molecules_keeps_every_molecule_once() {
    let cells = fit_v2(anchored_cells_with_vertices(3), v2_spec(1.0, 0.0));
    for (cell, original) in cells.iter().zip(anchored_cells(3)) {
        assert_eq!(cell.len(), original.len());
        for m in cell {
            assert_eq!(m.transcript_probs.len(), 1, "{:?}", m.transcript_probs);
            assert_eq!(m.transcript_probs[0].1, 1.0);
        }
    }
}

/// Cells whose two isoforms have the same length and the same read positions, so positions cannot
/// separate them: 30 molecules unique to isoform 0, 10 unique to isoform 1 and 200 shared, each with
/// its own vertex.
fn indistinguishable_cells(n: u32) -> Vec<Vec<ResolvedMolecule>> {
    let mut umi = 0u64;
    (0..n)
        .map(|_| {
            let mut cell: Vec<_> = (0..30)
                .map(|i| molecule(&[(0, 1.0, i % 5, 1000, 1000)]))
                .chain((0..10).map(|i| molecule(&[(1, 1.0, i % 5, 1000, 1000)])))
                .chain((0..200).map(|i| molecule(&[(0, 0.5, i % 5, 1000, 1000), (1, 0.5, i % 5, 1000, 1000)])))
                .collect();
            for m in cell.iter_mut() {
                umi += 1;
                m.vertices = vec![((umi % 97) as u32, umi * 2_654_435_761)];
            }
            cell
        })
        .collect()
}

#[test]
fn drawn_counts_match_the_soft_counts_in_expectation() {
    let soft = fit_v2(indistinguishable_cells(4), v2_spec(0.0, 0.0));
    let drawn = fit_v2(indistinguishable_cells(4), v2_spec(1.0, 0.0));
    let (mut expected, mut observed, mut variance) = (0.0, 0.0, 0.0);
    for (s, d) in soft.iter().zip(&drawn) {
        for (ms, md) in s[40..].iter().zip(&d[40..]) {
            let p = share(ms, 1);
            expected += p;
            variance += p * (1.0 - p);
            observed += share(md, 1);
        }
    }
    assert!(variance > 100.0, "the shared molecules should be ambiguous (variance {variance})");
    assert!(
        (observed - expected).abs() < 4.0 * variance.sqrt(),
        "drawn {observed} vs expected {expected} (sd {})",
        variance.sqrt()
    );
}

#[test]
fn drawn_counts_do_not_change_between_runs() {
    let a = fit_v2(anchored_cells_with_vertices(3), v2_spec(1.0, 0.0));
    let b = fit_v2(anchored_cells_with_vertices(3), v2_spec(1.0, 0.0));
    for (x, y) in a.iter().flatten().zip(b.iter().flatten()) {
        assert_eq!(x.transcript_probs, y.transcript_probs);
    }
}

#[test]
fn small_shares_are_dropped_and_the_rest_renormalized() {
    let cells = fit_v2(anchored_cells_with_vertices(3), v2_spec(0.0, 0.45));
    for m in cells.iter().flatten() {
        let total: f64 = m.transcript_probs.iter().map(|(_, p)| p).sum();
        assert!((total - 1.0).abs() < 1e-9);
        if m.transcript_probs.len() > 1 {
            assert!(m.transcript_probs.iter().all(|(_, p)| *p >= 0.45), "{:?}", m.transcript_probs);
        }
    }
}

/// Cells of two types. Each cell has 20 one-candidate molecules on one isoform and 2 on the other,
/// plus 20 molecules ambiguous between the two; with `differ` the second type has the one-candidate
/// counts the other way round, without it both types are the same. Positions are the same on both
/// isoforms, so the ambiguous molecules are decided by the abundances alone.
fn typed_cells(n_per_type: u32, differ: bool) -> Vec<Vec<ResolvedMolecule>> {
    (0..2 * n_per_type)
        .map(|c| {
            let (major, minor) = if differ && c >= n_per_type { (1, 0) } else { (0, 1) };
            (0..20)
                .map(|i| molecule(&[(major, 1.0, i % 5, 1000, 1000)]))
                .chain((0..2).map(|i| molecule(&[(minor, 1.0, i % 5, 1000, 1000)])))
                .chain((0..20).map(|i| molecule(&[(0, 0.5, i % 5, 1000, 1000), (1, 0.5, i % 5, 1000, 1000)])))
                .collect()
        })
        .collect()
}

fn auto_spec() -> TruncationModelSpec {
    TruncationModelSpec {
        end_kernels: TruncationEndKernels::Absolute,
        shrink_auto: true,
        ..TruncationModelSpec::default()
    }
}

/// Fits `abundance` to `typed_cells(10, differ)` and returns the share of isoform 0 that the first
/// cell's first ambiguous molecule gets, with the global fit.
fn ambiguous_share(differ: bool, abundance: TruncationAbundance) -> (f64, alevin_fry::truncation::GlobalFit) {
    let mut cells = typed_cells(10, differ);
    let mut input: Vec<_> = cells.iter().map(|m| cell_molecules(m, TruncationIntervalMode::All, 0, None)).collect();
    let global = fit_global(&mut input, 0, 1.0, 2, auto_spec(), abundance, 2, 1, FitExtras::default());
    let mut first = typed_cells(10, differ).swap_remove(0);
    let shrink = (abundance == TruncationAbundance::Shrink).then_some(1.0);
    fit_molecules_with(&mut first, TruncationIntervalMode::All, Some(&global), auto_spec(), shrink, 0, None);
    (share(&first[22], 0), global)
}

#[test]
fn auto_shrink_pools_fully_when_cells_share_usage() {
    let (auto, global) = ambiguous_share(false, TruncationAbundance::Shrink);
    let prior = global.prior.as_ref().expect("auto shrink estimates the prior strength");
    assert!(prior.alpha.is_infinite(), "identical cells give an infinite alpha, got {}", prior.alpha);
    assert_eq!(prior.loci, 1);
    assert_eq!(prior.cells, 20);
    let (pooled, _) = ambiguous_share(false, TruncationAbundance::Global);
    assert!((auto - pooled).abs() < 1e-9, "auto {auto} should equal the pooled assignment {pooled}");
}

#[test]
fn auto_shrink_lets_cells_keep_their_own_usage_when_they_differ() {
    let (auto, global) = ambiguous_share(true, TruncationAbundance::Shrink);
    let alpha = global.prior.as_ref().expect("auto shrink estimates the prior strength").alpha;
    assert!(alpha.is_finite() && alpha < 2.0, "two cell types with opposite usage give a weak prior, got alpha {alpha}");
    let (pooled, _) = ambiguous_share(true, TruncationAbundance::Global);
    assert!((pooled - 0.5).abs() < 0.02, "pooled shares split the ambiguous molecules evenly, got {pooled}");
    assert!(auto > 0.85, "a cell with 20 one-candidate molecules on isoform 0 keeps its ambiguous ones mostly on 0, got {auto}");
}

#[test]
fn gene_map_matches_names_with_and_without_versions() {
    use alevin_fry::truncation::GeneMap;
    let path = std::env::temp_dir().join(format!("truncation_genes_{}.tsv", std::process::id()));
    std::fs::write(&path, "NM_1.2\tG1\nNM_2\tG1\nNM_3.1\tG2\textra\n").unwrap();
    let names: Vec<String> = ["NM_1", "NM_2.5", "NM_3.1", "NM_9"].iter().map(|s| s.to_string()).collect();
    let g = GeneMap::load(&path, &names).unwrap();
    std::fs::remove_file(&path).ok();
    assert_eq!(g.genes, 2);
    assert_eq!(g.matched, 3);
    assert_eq!(g.gene_of[0], g.gene_of[1], "NM_1 and NM_2.5 are both G1");
    assert_ne!(g.gene_of[0], g.gene_of[2], "NM_3.1 is G2");
    assert!(g.gene_of[3] >= 2, "a transcript not in the file gets an index of its own");
}

#[test]
fn genes_split_the_loci_the_prior_acts_within() {
    let mut cells = typed_cells(10, true);
    let mut input: Vec<_> = cells.iter().map(|m| cell_molecules(m, TruncationIntervalMode::All, 0, None)).collect();
    // both transcripts in one gene: as without genes
    let same = fit_global(&mut input, 0, 1.0, 2, auto_spec(), TruncationAbundance::Shrink, 2, 1, FitExtras { genes: Some(&[0, 0]), ..Default::default() });
    let none = fit_global(&mut input, 0, 1.0, 2, auto_spec(), TruncationAbundance::Shrink, 2, 1, FitExtras::default());
    assert_eq!(same.prior.as_ref().unwrap().loci, 1);
    assert!((same.prior.as_ref().unwrap().alpha - none.prior.as_ref().unwrap().alpha).abs() < 1e-9);
    assert!(same.prior.as_ref().unwrap().within_genes && !none.prior.as_ref().unwrap().within_genes);
    // the two transcripts in different genes: no gene has two isoforms with one-candidate molecules,
    // so there is nothing to estimate from, and the prior stays the pooled shares
    let split = fit_global(&mut input, 0, 1.0, 2, auto_spec(), TruncationAbundance::Shrink, 2, 1, FitExtras { genes: Some(&[0, 1]), ..Default::default() });
    assert_eq!(split.prior.as_ref().unwrap().loci, 0);
    assert!(split.prior.as_ref().unwrap().alpha.is_infinite());
}

fn four_param_spec() -> TruncationModelSpec {
    TruncationModelSpec {
        positional_law: TruncationPositionalLaw::FourParam,
        ..TruncationModelSpec::default()
    }
}

#[test]
fn the_four_param_law_learns_three_prime_anchoring_per_cell() {
    // every read of `anchored_cells` ends at the transcript's 3' end and starts 200-800 nt in
    let mut cell = anchored_cells(1).swap_remove(0);
    let fit = fit_molecules(&mut cell, TruncationIntervalMode::All, four_param_spec());
    let four = fit.params.four.expect("the four-parameter law reports its parameters");
    assert!(four.iter().all(|v| v.is_finite()));
    assert!(four[2] > 0.9, "3'-anchored reads give a high 3'-anchored fraction, got {}", four[2]);
    assert!(four[0] < 0.1, "no read reaches the 5' end, got a full-length fraction of {}", four[0]);
}

#[test]
fn the_four_param_law_runs_with_pooled_abundances() {
    let cells = anchored_cells(4);
    let mut input: Vec<_> = cells.iter().map(|m| cell_molecules(m, TruncationIntervalMode::All, 0, None)).collect();
    let global = fit_global(&mut input, 0, 1.0, 2, four_param_spec(), TruncationAbundance::Global, 2, 1, FitExtras::default());
    assert!(global.params.four.is_some());
    let mut first = anchored_cells(4).swap_remove(0);
    let fit = fit_molecules_with(&mut first, TruncationIntervalMode::All, Some(&global), four_param_spec(), None, 0, None);
    assert_eq!(fit.params.four, global.params.four, "cells report the global parameters");
    let total: f64 = first.iter().map(|m| m.transcript_probs.iter().map(|(_, p)| p).sum::<f64>()).sum();
    assert!((total - first.len() as f64).abs() < 1e-6, "every molecule keeps one count, got {total}");
}

#[test]
fn a_coverage_model_reweights_the_candidates_before_the_positional_model() {
    use alevin_fry::eq_class::CellTxpCoverage;
    // two transcripts of 1,000 nt; the second's coverage probability is low everywhere
    let mut coverage = vec![CellTxpCoverage::new(1000, 100), CellTxpCoverage::new(1000, 100)];
    coverage[0].coverage_prob = vec![1.0; 10];
    coverage[1].coverage_prob = vec![0.1; 10];
    let make = || {
        let mut cell: Vec<_> = (0..10).map(|_| molecule(&[(0, 1.0, 100, 1000, 1000)])).collect();
        cell.extend((0..10).map(|_| molecule(&[(1, 1.0, 100, 1000, 1000)])));
        cell.push(molecule(&[(0, 0.5, 100, 1000, 1000), (1, 0.5, 100, 1000, 1000)]));
        cell
    };
    let mut plain = make();
    fit_molecules_with(&mut plain, TruncationIntervalMode::Ignore, None, spec(), None, 0, None);
    let mut weighted = make();
    fit_molecules_with(&mut weighted, TruncationIntervalMode::Ignore, None, spec(), None, 0, Some(&coverage));
    assert!((share(&plain[20], 0) - 0.5).abs() < 0.05, "without coverage the ambiguous molecule splits evenly, got {}", share(&plain[20], 0));
    assert!(share(&weighted[20], 0) > 0.85, "with coverage it goes to the well-covered transcript, got {}", share(&weighted[20], 0));
}

#[test]
fn the_junction_law_prefers_starts_at_junctions_when_the_data_do() {
    // transcript 0 has a junction at 500; transcript 1 has none; both are 1,000 nt. Record starts are
    // 1-based, so a read whose 0-based start is the junction has start 501.
    let junctions = vec![vec![500u32], vec![]];
    let make = || -> Vec<Vec<ResolvedMolecule>> {
        (0..4)
            .map(|_| {
                let mut cell: Vec<_> = (0..20).map(|_| molecule(&[(0, 1.0, 501, 1000, 1000)])).collect();
                cell.extend((0..20).map(|_| molecule(&[(1, 1.0, 501, 1000, 1000)])));
                cell.push(molecule(&[(0, 0.5, 501, 1000, 1000), (1, 0.5, 501, 1000, 1000)]));
                cell
            })
            .collect()
    };
    let cells = make();
    let mut input: Vec<_> = cells.iter().map(|m| cell_molecules(m, TruncationIntervalMode::All, 0, None)).collect();
    let global = fit_global(&mut input, 0, 1.0, 2, spec(), TruncationAbundance::Global, 2, 1, FitExtras { junctions: Some(&junctions), ..Default::default() });
    let law = global.junction.as_ref().expect("a junction law is learned");
    assert!(law.molecules > 0);
    assert!(law.factor(0, 500) > 10.0 * law.factor(1, 500), "a start at a junction is far more likely than one 500 nt into an exon");
    let mut first = make().swap_remove(0);
    fit_molecules_with(&mut first, TruncationIntervalMode::All, Some(&global), spec(), None, 0, None);
    assert!(share(&first[40], 0) > 0.6, "the ambiguous molecule favours the isoform where its start sits at a junction, got {}", share(&first[40], 0));
}

#[test]
fn auto_gene_estimates_the_prior_per_gene() {
    let gene_spec = TruncationModelSpec {
        shrink_per_gene: true,
        ..auto_spec()
    };
    let fit = |differ: bool| {
        let cells = typed_cells(10, differ);
        let mut input: Vec<_> = cells.iter().map(|m| cell_molecules(m, TruncationIntervalMode::All, 0, None)).collect();
        fit_global(&mut input, 0, 1.0, 2, gene_spec, TruncationAbundance::Shrink, 2, 1, FitExtras { genes: Some(&[0, 0]), ..Default::default() })
    };
    let same = fit(false);
    let p = same.prior.as_ref().unwrap();
    assert!(p.per_locus.as_ref().unwrap().iter().all(|a| a.is_infinite()), "identical cells: every gene pools fully");
    let differ = fit(true);
    let p = differ.prior.as_ref().unwrap();
    assert_eq!(p.loci_finite, 1);
    assert!(p.alpha_median.is_finite() && p.alpha_median < 5.0, "two types with opposite usage give the gene a weak prior, got {}", p.alpha_median);
    let mut first = typed_cells(10, true).swap_remove(0);
    fit_molecules_with(&mut first, TruncationIntervalMode::All, Some(&differ), gene_spec, Some(1.0), 0, None);
    assert!(share(&first[22], 0) > 0.85, "the cell keeps its own usage, got {}", share(&first[22], 0));
}

fn group_spec(auto: bool, k: f64) -> TruncationModelSpec {
    TruncationModelSpec {
        end_kernels: TruncationEndKernels::Absolute,
        group_shrink: k,
        group_shrink_auto: auto,
        ..TruncationModelSpec::default()
    }
}

/// `typed_cells(10, differ)` in two groups of 10 (the second group is the one whose usage differs),
/// pooled per group under `spec`; returns the share of isoform 0 the first cell's first ambiguous
/// molecule gets in group 0 and in group 1, and the fit.
fn grouped_share(differ: bool, spec: TruncationModelSpec) -> (f64, f64, alevin_fry::truncation::GlobalFit) {
    let cells = typed_cells(10, differ);
    let mut input: Vec<_> = cells
        .iter()
        .enumerate()
        .map(|(i, m)| cell_molecules(m, TruncationIntervalMode::All, (i / 10) as u32, None))
        .collect();
    let global = fit_global(&mut input, 0, 1.0, 2, spec, TruncationAbundance::Global, 2, 2, FitExtras::default());
    let mut shares = [0.0; 2];
    for (g, s) in shares.iter_mut().enumerate() {
        let mut first = typed_cells(10, differ).swap_remove(0);
        fit_molecules_with(&mut first, TruncationIntervalMode::All, Some(&global), spec, None, g as u32, None);
        *s = share(&first[22], 0);
    }
    (shares[0], shares[1], global)
}

#[test]
fn a_group_prior_pools_identical_groups_fully() {
    let (s0, s1, fit) = grouped_share(false, group_spec(true, 0.0));
    let p = fit.group_prior.as_ref().expect("the group prior strength is estimated");
    assert!(p.alpha.is_infinite(), "groups with the same usage give an infinite beta, got {}", p.alpha);
    assert_eq!(p.cells, 2, "two groups were compared");
    assert!((s0 - s1).abs() < 1e-6, "both groups then use the shares of all cells: {s0} vs {s1}");
    let (n0, _, none) = grouped_share(false, group_spec(false, 0.0));
    assert!(none.group_prior.is_none());
    assert!((s0 - n0).abs() < 0.02, "which are the independent groups' own shares here: {s0} vs {n0}");
}

#[test]
fn a_group_prior_keeps_groups_apart_only_when_their_usage_differs() {
    // group 0 has isoform 0 as its major isoform, group 1 isoform 1
    let (s0, s1, fit) = grouped_share(true, group_spec(true, 0.0));
    let beta = fit.group_prior.as_ref().expect("estimated").alpha;
    assert!(beta.is_finite() && beta < 20.0, "groups with opposite usage give a weak prior, got beta {beta}");
    assert!(s0 > 0.8 && s1 < 0.2, "each group keeps its own usage: {s0} and {s1}");
    // a very strong given prior pulls both groups to the shares of all cells
    let (t0, t1, _) = grouped_share(true, group_spec(false, 1e7));
    assert!((t0 - 0.5).abs() < 0.02 && (t1 - 0.5).abs() < 0.02, "shares of all cells: {t0} and {t1}");
    // a moderate given prior lies in between
    let (m0, _, _) = grouped_share(true, group_spec(false, 220.0));
    assert!(m0 > 0.55 && m0 < s0, "220 pseudo-molecules against 220 molecules halve the distance: {m0}");
}
