use std::collections::HashMap;

use alevin_fry::eq_class::AlignmentInterval;
use alevin_fry::pugutils_dp_new::ResolvedMolecule;
use alevin_fry::truncation::{cell_molecules, TruncationGroups, TruncationIntervalMode, TruncationPreset};
use alevin_fry::truncation_groups::{discover_groups, DiscoverySettings};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Binomial, Distribution};

/// Cells of `types.len()` types (so many cells each) over `n_genes` genes with two isoforms each,
/// about eight one-candidate molecules per gene and cell; `swapped[g]` genes use the minor isoform
/// as the major one in the second type. Returns the counts as `discover_groups` takes them, the
/// genes of the transcripts, and each cell's type.
fn synthetic(types: &[usize], n_genes: usize, swapped: usize, seed: u64) -> (Vec<Vec<(u32, u32)>>, Vec<u32>, Vec<usize>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let gene_of: Vec<u32> = (0..2 * n_genes).map(|t| (t / 2) as u32).collect();
    let mut counts = Vec::new();
    let mut labels = Vec::new();
    for (ty, &n) in types.iter().enumerate() {
        for _ in 0..n {
            let mut cell = Vec::new();
            for g in 0..n_genes {
                let n_mol = Binomial::new(16, 0.5).unwrap().sample(&mut rng);
                if n_mol == 0 {
                    continue;
                }
                let major_share = if ty == 1 && g < swapped { 0.2 } else { 0.8 };
                let a = Binomial::new(n_mol, major_share).unwrap().sample(&mut rng) as u32;
                let b = n_mol as u32 - a;
                if a > 0 {
                    cell.push((2 * g as u32, a));
                }
                if b > 0 {
                    cell.push((2 * g as u32 + 1, b));
                }
            }
            counts.push(cell);
            labels.push(ty);
        }
    }
    // shuffle the cells so that the types are not contiguous
    for i in (1..counts.len()).rev() {
        let j = rng.random_range(0..=i);
        counts.swap(i, j);
        labels.swap(i, j);
    }
    (counts, gene_of, labels)
}

fn settings(threads: usize) -> DiscoverySettings {
    DiscoverySettings {
        threads,
        ..DiscoverySettings::default()
    }
}

/// The share of cells whose group agrees with their type under the best one-to-one matching of
/// two groups to two types.
fn agreement(labels: &[u32], truth: &[usize]) -> f64 {
    let same = labels.iter().zip(truth).filter(|&(&l, &t)| l as usize == t).count();
    let flipped = labels.iter().zip(truth).filter(|&(&l, &t)| l as usize == 1 - t).count();
    same.max(flipped) as f64 / labels.len() as f64
}

#[test]
fn cells_that_share_usage_form_one_group() {
    let (counts, gene_of, _) = synthetic(&[120, 120], 400, 0, 7);
    let d = discover_groups(&counts, &gene_of, &settings(2));
    assert_eq!(d.n_groups, 1, "splits: {:?}", d.splits);
    assert!(d.labels.iter().all(|&l| l == 0));
    assert!(d.splits.iter().all(|s| !s.accepted));
    assert_eq!(d.genes, 400);
    assert_eq!(d.transcripts, 800);
}

#[test]
fn two_cell_types_are_recovered() {
    let (counts, gene_of, truth) = synthetic(&[120, 120], 400, 120, 11);
    let d = discover_groups(&counts, &gene_of, &settings(2));
    assert_eq!(d.n_groups, 2, "splits: {:?}", d.splits);
    assert!(d.splits.iter().any(|s| s.accepted && s.depth == 0));
    assert!(
        agreement(&d.labels, &truth) > 0.97,
        "agreement {} with sizes {:?}",
        agreement(&d.labels, &truth),
        d.sizes
    );
    assert_eq!(d.sizes.iter().sum::<usize>(), 240);
}

#[test]
fn discovery_does_not_depend_on_the_number_of_threads() {
    let (counts, gene_of, _) = synthetic(&[60, 60], 300, 80, 3);
    let a = discover_groups(&counts, &gene_of, &settings(1));
    let b = discover_groups(&counts, &gene_of, &settings(4));
    assert_eq!(a.labels, b.labels);
    assert_eq!(a.splits.len(), b.splits.len());
    for (x, y) in a.splits.iter().zip(&b.splits) {
        assert_eq!(x.to_string(), y.to_string());
    }
}

#[test]
fn a_small_group_is_not_split_off_below_min_cells() {
    // the second type has fewer cells than min_cells, so no side can hold it
    let (counts, gene_of, _) = synthetic(&[100, 5], 300, 100, 5);
    let d = discover_groups(&counts, &gene_of, &settings(2));
    assert_eq!(d.n_groups, 1, "splits: {:?}", d.splits);
}

#[test]
fn empty_input_gives_one_group() {
    let d = discover_groups(&[], &[], &settings(1));
    assert_eq!(d.n_groups, 1);
    assert!(d.labels.is_empty());
    let d = discover_groups(&[vec![], vec![(0, 3)]], &[0, 0], &settings(1));
    assert_eq!(d.n_groups, 1);
    assert_eq!(d.labels, vec![0, 0]);
}

#[test]
fn groups_from_labels_name_and_map_the_cells() {
    let barcodes: Vec<String> = ["AAAA", "CCCC", "GGGG"].iter().map(|s| s.to_string()).collect();
    let g = TruncationGroups::from_labels(&barcodes, &[2, 0, 2]);
    assert_eq!(g.names(), &["L0".to_string(), "L1".to_string(), "L2".to_string()]);
    assert_eq!(g.n_groups(), 4);
    assert_eq!(g.n_barcodes(), 3);
    assert_eq!(g.group_of("AAAA"), 2);
    assert_eq!(g.group_of("CCCC"), 0);
    assert_eq!(g.group_of("TTTT"), 3);
}

#[test]
fn the_preset_parses_and_prints() {
    assert_eq!("recommended".parse::<TruncationPreset>(), Ok(TruncationPreset::Recommended));
    assert_eq!("Recommended".parse::<TruncationPreset>(), Ok(TruncationPreset::Recommended));
    assert!("best".parse::<TruncationPreset>().is_err());
    assert_eq!(TruncationPreset::Recommended.to_string(), "recommended");
}

fn molecule(target_intervals: &[(u32, f64, u32, u32, u32)]) -> ResolvedMolecule {
    let mut transcript_probs = Vec::new();
    let mut intervals = HashMap::new();
    for &(target, probability, start, end, len) in target_intervals {
        if !intervals.contains_key(&target) {
            transcript_probs.push((target, probability));
        }
        intervals.entry(target).or_insert_with(Vec::new).push(AlignmentInterval {
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
fn one_candidate_counts_count_the_unambiguous_molecules_per_transcript() {
    let molecules: Vec<_> = (0..3)
        .map(|_| molecule(&[(7, 1.0, 0, 900, 1000)]))
        .chain((0..2).map(|_| molecule(&[(3, 1.0, 0, 500, 1000)])))
        .chain((0..4).map(|_| molecule(&[(7, 0.5, 0, 900, 1000), (3, 0.5, 0, 900, 1000)])))
        // two intervals on the same transcript are still one candidate
        .chain(std::iter::once(molecule(&[(3, 1.0, 0, 200, 1000), (3, 1.0, 300, 500, 1000)])))
        .collect();
    let cell = cell_molecules(&molecules, TruncationIntervalMode::All, 0, None);
    assert_eq!(cell.one_candidate_counts(), vec![(3, 3), (7, 3)]);
}
