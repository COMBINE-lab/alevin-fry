/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Deterministic barcode-correction policy and compilation.
//!
//! This module decides corrections. Consumers such as collation receive the
//! resulting sorted table and never enumerate neighbours or choose between
//! competing targets themselves.

use ahash::AHashMap;
use anyhow::{anyhow, bail};
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use std::collections::BTreeMap;
use std::fmt;
use std::hash::Hash;
use std::str::FromStr;

/// Which fixed-length, one-error neighbourhood to search.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum BarcodeNeighborhood {
    /// One nucleotide substitution (Hamming distance one).
    HammingOne,
    /// One substitution or a fixed-length one-base shift.
    ///
    /// The shift operations preserve the barcode length by dropping a base at
    /// one end and admitting any base at the other. This is the historical
    /// alevin-fry `edit-1` behaviour; it is not general Levenshtein distance.
    SubstitutionOrShiftOne,
}

impl fmt::Display for BarcodeNeighborhood {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::HammingOne => f.write_str("hamming-1"),
            Self::SubstitutionOrShiftOne => f.write_str("substitution-or-shift-1"),
        }
    }
}

impl FromStr for BarcodeNeighborhood {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "hamming-1" => Ok(Self::HammingOne),
            "substitution-or-shift-1" | "edit-1" => Ok(Self::SubstitutionOrShiftOne),
            _ => bail!(
                "unknown barcode neighbourhood '{value}'; expected hamming-1 or substitution-or-shift-1"
            ),
        }
    }
}

/// An exact rational confidence threshold.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Confidence {
    numerator: u64,
    denominator: u64,
}

impl Confidence {
    pub const RNA: Self = Self {
        numerator: 39,
        denominator: 40,
    };
    pub const ATAC: Self = Self {
        numerator: 9,
        denominator: 10,
    };

    pub fn new(numerator: u64, denominator: u64) -> anyhow::Result<Self> {
        if denominator == 0 || numerator > denominator {
            bail!(
                "barcode-correction confidence must be between zero and one (got {numerator}/{denominator})"
            );
        }
        let divisor = gcd(numerator, denominator);
        Ok(Self {
            numerator: numerator / divisor,
            denominator: denominator / divisor,
        })
    }

    pub fn numerator(self) -> u64 {
        self.numerator
    }

    pub fn denominator(self) -> u64 {
        self.denominator
    }

    fn accepts(self, winner: u128, total: u128) -> bool {
        compare_fractions(
            winner,
            total,
            u128::from(self.numerator),
            u128::from(self.denominator),
        ) != std::cmp::Ordering::Less
    }
}

/// Compare positive rational values without cross-multiplication overflow.
fn compare_fractions(
    mut lhs_num: u128,
    mut lhs_den: u128,
    mut rhs_num: u128,
    mut rhs_den: u128,
) -> std::cmp::Ordering {
    debug_assert!(lhs_den != 0 && rhs_den != 0);
    let mut reversed = false;
    loop {
        let lhs_integer = lhs_num / lhs_den;
        let rhs_integer = rhs_num / rhs_den;
        if lhs_integer != rhs_integer {
            let ordering = lhs_integer.cmp(&rhs_integer);
            return if reversed {
                ordering.reverse()
            } else {
                ordering
            };
        }

        let lhs_remainder = lhs_num % lhs_den;
        let rhs_remainder = rhs_num % rhs_den;
        match (lhs_remainder == 0, rhs_remainder == 0) {
            (true, true) => return std::cmp::Ordering::Equal,
            (true, false) => {
                return if reversed {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Less
                };
            }
            (false, true) => {
                return if reversed {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Greater
                };
            }
            (false, false) => {
                lhs_num = lhs_den;
                lhs_den = lhs_remainder;
                rhs_num = rhs_den;
                rhs_den = rhs_remainder;
                reversed = !reversed;
            }
        }
    }
}

impl fmt::Display for Confidence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}", self.numerator, self.denominator)
    }
}

impl FromStr for Confidence {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let value = value.trim();
        if let Some((numerator, denominator)) = value.split_once('/') {
            return Self::new(numerator.parse()?, denominator.parse()?);
        }

        let (whole, fractional) = value.split_once('.').unwrap_or((value, ""));
        if whole.is_empty() || !whole.bytes().all(|b| b.is_ascii_digit()) {
            bail!("invalid barcode-correction confidence '{value}'");
        }
        if !fractional.bytes().all(|b| b.is_ascii_digit()) || fractional.len() > 18 {
            bail!("invalid barcode-correction confidence '{value}'");
        }
        let denominator = 10u64
            .checked_pow(fractional.len() as u32)
            .ok_or_else(|| anyhow!("barcode-correction confidence is too precise"))?;
        let whole: u64 = whole.parse()?;
        let fractional: u64 = if fractional.is_empty() {
            0
        } else {
            fractional.parse()?
        };
        let numerator = whole
            .checked_mul(denominator)
            .and_then(|v| v.checked_add(fractional))
            .ok_or_else(|| anyhow!("barcode-correction confidence is too large"))?;
        Self::new(numerator, denominator)
    }
}

const fn gcd(mut lhs: u64, mut rhs: u64) -> u64 {
    while rhs != 0 {
        let remainder = lhs % rhs;
        lhs = rhs;
        rhs = remainder;
    }
    if lhs == 0 { 1 } else { lhs }
}

/// How collisions between distinct canonical targets are resolved.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum BarcodeResolution {
    /// Accept only one distinct canonical target.
    Unique,
    /// Use frozen exact-source counts as priors.
    Frequency {
        confidence: Confidence,
        pseudocount: u64,
    },
}

/// Complete, resolved correction policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct CorrectionSpec {
    pub barcode_len: u8,
    pub neighborhood: BarcodeNeighborhood,
    pub resolution: BarcodeResolution,
}

impl CorrectionSpec {
    pub fn validate(self) -> anyhow::Result<Self> {
        if !(1..=32).contains(&self.barcode_len) {
            bail!(
                "barcode length must be between 1 and 32 (got {})",
                self.barcode_len
            );
        }
        if let BarcodeResolution::Frequency { pseudocount, .. } = self.resolution
            && pseudocount == 0
        {
            bail!("frequency correction requires a non-zero pseudocount");
        }
        Ok(self)
    }
}

/// A permitted source barcode, its canonical target, and its frozen raw count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RetainedSource<T> {
    pub source: u64,
    pub target: T,
    pub exact_count: u64,
}

#[derive(Debug, Clone, Copy)]
struct IndexedSource<T> {
    target: T,
    exact_count: u64,
}

/// The deterministic result for one observed barcode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorrectionDecision<T> {
    Exact(T),
    Corrected(T),
    Ambiguous,
    NotFound,
}

impl<T> CorrectionDecision<T> {
    pub fn target(self) -> Option<T> {
        match self {
            Self::Exact(target) | Self::Corrected(target) => Some(target),
            Self::Ambiguous | Self::NotFound => None,
        }
    }
}

/// Counts for both distinct observed barcodes and the reads they represent.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct CorrectionStats {
    pub exact_distinct: u64,
    pub exact_reads: u64,
    pub corrected_distinct: u64,
    pub corrected_reads: u64,
    pub ambiguous_distinct: u64,
    pub ambiguous_reads: u64,
    pub not_found_distinct: u64,
    pub not_found_reads: u64,
}

/// A sorted, accepted observed-to-canonical correction table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorrectionTable<T> {
    entries: Vec<(u64, T)>,
    stats: CorrectionStats,
}

/// A compiled, policy-free view of the complete one-error topology.
///
/// Accepted entries have exactly one distinct canonical target. Barcodes in
/// `ambiguous` have more than one structurally possible target. Both vectors
/// are sorted, which lets callers choose an execution lookup without exposing
/// one from the correction engine.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructuralCorrectionTable<T> {
    entries: Vec<(u64, T)>,
    ambiguous: Vec<u64>,
}

impl<T> StructuralCorrectionTable<T> {
    pub fn entries(&self) -> &[(u64, T)] {
        &self.entries
    }

    pub fn ambiguous(&self) -> &[u64] {
        &self.ambiguous
    }
}

impl<T> CorrectionTable<T> {
    pub fn entries(&self) -> &[(u64, T)] {
        &self.entries
    }

    pub fn into_entries(self) -> Vec<(u64, T)> {
        self.entries
    }

    pub fn stats(&self) -> CorrectionStats {
        self.stats
    }
}

/// An index of retained sources used to resolve and compile corrections.
///
/// Its concrete hash table and hasher are private so callers depend only on
/// deterministic correction semantics, not an implementation detail.
#[derive(Debug, Clone)]
pub struct CorrectionIndex<T> {
    spec: CorrectionSpec,
    /// Retained sources sorted by packed barcode. The prefix offsets keep
    /// one-substitution candidate lookups in small contiguous slices without
    /// exposing or depending on a particular hasher.
    sources: Vec<(u64, IndexedSource<T>)>,
    source_offsets: Vec<usize>,
    source_suffix_bits: u8,
}

impl<T> CorrectionIndex<T>
where
    T: Copy + Eq + Ord + Hash,
{
    pub fn new<I>(spec: CorrectionSpec, sources: I) -> anyhow::Result<Self>
    where
        I: IntoIterator<Item = RetainedSource<T>>,
    {
        let spec = spec.validate()?;
        let mut retained_sources: Vec<_> = sources.into_iter().collect();
        retained_sources.sort_unstable_by_key(|retained| retained.source);
        let mut index: Vec<(u64, IndexedSource<T>)> = Vec::with_capacity(retained_sources.len());
        for retained in retained_sources {
            validate_barcode(retained.source, spec.barcode_len)?;
            match index.last() {
                Some((source, previous))
                    if *source == retained.source && previous.target != retained.target =>
                {
                    bail!(
                        "retained source barcode {} was assigned conflicting canonical targets",
                        retained.source
                    )
                }
                Some((source, previous))
                    if *source == retained.source
                        && previous.exact_count != retained.exact_count =>
                {
                    bail!(
                        "retained source barcode {} was assigned conflicting exact counts",
                        retained.source
                    )
                }
                Some((source, _)) if *source == retained.source => {}
                None => {
                    index.push((
                        retained.source,
                        IndexedSource {
                            target: retained.target,
                            exact_count: retained.exact_count,
                        },
                    ));
                }
                Some(_) => index.push((
                    retained.source,
                    IndexedSource {
                        target: retained.target,
                        exact_count: retained.exact_count,
                    },
                )),
            }
        }
        // Cap the prefix at eight bases (65,536 slots). This preserves the
        // narrow lookup ranges used for ordinary 16-base cell barcodes while
        // remaining safe for any supported barcode length through 32.
        let source_prefix_len = spec.barcode_len.div_ceil(2).min(8);
        let source_suffix_bits = 2 * (spec.barcode_len - source_prefix_len);
        let prefix_count = 1usize << (2 * source_prefix_len);
        let mut source_offsets = vec![0; prefix_count + 1];
        let mut cursor = 0usize;
        for (prefix, offset) in source_offsets.iter_mut().enumerate() {
            while cursor < index.len() && (index[cursor].0 >> source_suffix_bits) < prefix as u64 {
                cursor += 1;
            }
            *offset = cursor;
        }
        Ok(Self {
            spec,
            sources: index,
            source_offsets,
            source_suffix_bits,
        })
    }

    pub fn spec(&self) -> CorrectionSpec {
        self.spec
    }

    /// Resolve one barcode according to the configured policy.
    pub fn resolve(&self, observed: u64) -> CorrectionDecision<T> {
        if let Some(source) = self.source_for(observed) {
            return CorrectionDecision::Exact(source.target);
        }

        if self.spec.resolution == BarcodeResolution::Unique {
            return self.unique_neighbor_target(observed);
        }

        let candidates = self.candidate_sources(observed);
        if candidates.is_empty() {
            return CorrectionDecision::NotFound;
        }

        let mut targets: SmallVec<[(T, u128); 8]> = SmallVec::new();
        for source_barcode in candidates {
            let source = self
                .source_for(source_barcode)
                .expect("candidate source is retained");
            if let Some((_, weight)) = targets
                .iter_mut()
                .find(|(target, _)| *target == source.target)
            {
                *weight += u128::from(source.exact_count) + u128::from(self.pseudocount());
            } else {
                targets.push((
                    source.target,
                    u128::from(source.exact_count) + u128::from(self.pseudocount()),
                ));
            }
        }

        let BarcodeResolution::Frequency { confidence, .. } = self.spec.resolution else {
            unreachable!("Unique correction returned through its allocation-free path")
        };
        let total = targets.iter().map(|(_, weight)| *weight).sum();
        // Use the greatest target as the deterministic final tie-break,
        // preserving the tuple-max behaviour of the precursor policy.
        let (winner, winner_weight) = targets
            .into_iter()
            .max_by(|lhs, rhs| lhs.1.cmp(&rhs.1).then_with(|| lhs.0.cmp(&rhs.0)))
            .expect("candidate target list is non-empty");
        if confidence.accepts(winner_weight, total) {
            CorrectionDecision::Corrected(winner)
        } else {
            CorrectionDecision::Ambiguous
        }
    }

    /// Resolve only the candidate topology, without applying an abundance
    /// policy. This is used by sample-Frequency routing to defer only pairs
    /// having more than one structurally possible canonical sample.
    pub fn structural_resolution(&self, observed: u64) -> CorrectionDecision<T> {
        if let Some(source) = self.source_for(observed) {
            return CorrectionDecision::Exact(source.target);
        }
        self.unique_neighbor_target(observed)
    }

    /// Compile decisions for observed barcode counts. Retained identities are
    /// always included, even when their observed count is zero.
    pub fn compile_observed<I>(&self, observations: I) -> CorrectionTable<T>
    where
        I: IntoIterator<Item = (u64, u64)>,
    {
        let mut observed_counts = BTreeMap::<u64, u64>::new();
        for (barcode, count) in observations {
            let entry = observed_counts.entry(barcode).or_default();
            *entry = entry
                .checked_add(count)
                .expect("observed barcode count overflow");
        }
        for &(source, _) in &self.sources {
            observed_counts.entry(source).or_default();
        }

        let mut entries = Vec::with_capacity(observed_counts.len());
        let mut stats = CorrectionStats::default();
        for (barcode, count) in observed_counts {
            match self.resolve(barcode) {
                CorrectionDecision::Exact(target) => {
                    entries.push((barcode, target));
                    stats.exact_distinct += 1;
                    stats.exact_reads += count;
                }
                CorrectionDecision::Corrected(target) => {
                    entries.push((barcode, target));
                    stats.corrected_distinct += 1;
                    stats.corrected_reads += count;
                }
                CorrectionDecision::Ambiguous => {
                    stats.ambiguous_distinct += 1;
                    stats.ambiguous_reads += count;
                }
                CorrectionDecision::NotFound => {
                    stats.not_found_distinct += 1;
                    stats.not_found_reads += count;
                }
            }
        }
        CorrectionTable { entries, stats }
    }

    /// Compile a histogram whose observed barcodes are already distinct while
    /// aggregating accepted read counts by canonical target in the same pass.
    ///
    /// GPL's raw histograms satisfy this stronger contract. Keeping this path
    /// separate from [`Self::compile_observed`] avoids rebuilding a large,
    /// already-distinct hash histogram as a tree and avoids a second accepted
    /// lookup table solely for frequency aggregation. The returned target
    /// counts are sorted, so the correction engine still does not expose its
    /// private hash-table implementation.
    pub(crate) fn compile_distinct_observed_with_target_counts<I>(
        &self,
        observations: I,
    ) -> (CorrectionTable<T>, Vec<(T, u64)>)
    where
        I: IntoIterator<Item = (u64, u64)>,
    {
        let observations = observations.into_iter();
        let mut entries = Vec::with_capacity(observations.size_hint().0);
        let mut target_counts = AHashMap::<T, u64>::new();
        let mut stats = CorrectionStats::default();

        for (barcode, count) in observations {
            let target = match self.resolve(barcode) {
                CorrectionDecision::Exact(target) => {
                    stats.exact_distinct += 1;
                    stats.exact_reads += count;
                    Some(target)
                }
                CorrectionDecision::Corrected(target) => {
                    stats.corrected_distinct += 1;
                    stats.corrected_reads += count;
                    Some(target)
                }
                CorrectionDecision::Ambiguous => {
                    stats.ambiguous_distinct += 1;
                    stats.ambiguous_reads += count;
                    None
                }
                CorrectionDecision::NotFound => {
                    stats.not_found_distinct += 1;
                    stats.not_found_reads += count;
                    None
                }
            };
            if let Some(target) = target {
                entries.push((barcode, target));
                let target_count = target_counts.entry(target).or_default();
                *target_count = target_count
                    .checked_add(count)
                    .expect("corrected target count overflow");
            }
        }

        entries.sort_unstable_by_key(|&(barcode, _)| barcode);
        let mut missing_sources = Vec::new();
        for &(source, indexed) in &self.sources {
            if entries
                .binary_search_by_key(&source, |&(barcode, _)| barcode)
                .is_err()
            {
                missing_sources.push((source, indexed.target));
                stats.exact_distinct += 1;
            }
        }
        if !missing_sources.is_empty() {
            entries.extend(missing_sources);
            entries.sort_unstable_by_key(|&(barcode, _)| barcode);
        }
        debug_assert!(entries.windows(2).all(|pair| pair[0].0 != pair[1].0));

        let mut target_counts: Vec<_> = target_counts.into_iter().collect();
        target_counts.sort_unstable_by_key(|&(target, _)| target);
        (CorrectionTable { entries, stats }, target_counts)
    }

    /// Compile the complete theoretical one-error map used by historical
    /// `permit_map.bin` consumers.
    pub fn compile_full_neighborhood(&self) -> CorrectionTable<T> {
        let observed = self.theoretical_observations();
        self.compile_observed(observed.into_iter().map(|barcode| (barcode, 0)))
    }

    /// Compile all structurally reachable observations before applying a
    /// collision policy. This is useful when a one-pass caller must route
    /// exact/unique observations immediately and defer only ambiguous ones.
    pub fn compile_structural_neighborhood(&self) -> StructuralCorrectionTable<T> {
        let observed = self.theoretical_observations();
        let mut entries = Vec::with_capacity(observed.len());
        let mut ambiguous = Vec::new();
        for barcode in observed {
            match self.structural_resolution(barcode) {
                CorrectionDecision::Exact(target) | CorrectionDecision::Corrected(target) => {
                    entries.push((barcode, target));
                }
                CorrectionDecision::Ambiguous => ambiguous.push(barcode),
                CorrectionDecision::NotFound => {}
            }
        }
        StructuralCorrectionTable { entries, ambiguous }
    }

    fn theoretical_observations(&self) -> Vec<u64> {
        let candidates_per_source = 1 + 3 * usize::from(self.spec.barcode_len) + 8;
        let mut observed = Vec::with_capacity(self.sources.len() * candidates_per_source);
        for &(source, _) in &self.sources {
            observed.push(source);
            for_each_neighbor(
                source,
                self.spec.barcode_len,
                self.spec.neighborhood,
                |neighbor| observed.push(neighbor),
            );
        }
        observed.sort_unstable();
        observed.dedup();
        observed
    }

    fn pseudocount(&self) -> u64 {
        match self.spec.resolution {
            BarcodeResolution::Unique => 0,
            BarcodeResolution::Frequency { pseudocount, .. } => pseudocount,
        }
    }

    /// Find whether all neighbouring retained sources agree on one canonical
    /// target. Repeated constructions and aliases of that same target are
    /// harmless, so Unique and topology-only routing need no candidate vector,
    /// sort, or deduplication.
    fn unique_neighbor_target(&self, observed: u64) -> CorrectionDecision<T> {
        let mut target = None;
        let mut ambiguous = false;
        for_each_substitution(observed, self.spec.barcode_len, |source_barcode| {
            if !ambiguous && let Some(source) = self.source_for(source_barcode) {
                ambiguous = merge_unique_target(&mut target, source.target);
            }
        });
        if !ambiguous && self.spec.neighborhood == BarcodeNeighborhood::SubstitutionOrShiftOne {
            for_each_inverse_shift_candidate(observed, self.spec.barcode_len, |source_barcode| {
                if !ambiguous && let Some(source) = self.source_for(source_barcode) {
                    ambiguous = merge_unique_target(&mut target, source.target);
                }
            });
        }

        if ambiguous {
            CorrectionDecision::Ambiguous
        } else if let Some(target) = target {
            CorrectionDecision::Corrected(target)
        } else {
            CorrectionDecision::NotFound
        }
    }

    fn candidate_sources(&self, observed: u64) -> SmallVec<[u64; 128]> {
        let mut candidates = SmallVec::<[u64; 128]>::new();
        for_each_substitution(observed, self.spec.barcode_len, |candidate| {
            if self.source_for(candidate).is_some() {
                candidates.push(candidate);
            }
        });
        if self.spec.neighborhood == BarcodeNeighborhood::SubstitutionOrShiftOne {
            for_each_inverse_shift_candidate(observed, self.spec.barcode_len, |source| {
                if self.source_for(source).is_some() {
                    candidates.push(source);
                }
            });
        }
        // Shift candidates can repeat, especially for homopolymers. A source
        // prior must contribute once regardless of how many constructions
        // lead to it.
        candidates.sort_unstable();
        candidates.dedup();
        candidates
    }

    #[inline]
    fn source_for(&self, barcode: u64) -> Option<IndexedSource<T>> {
        let prefix = usize::try_from(barcode >> self.source_suffix_bits).ok()?;
        if prefix + 1 >= self.source_offsets.len() {
            return None;
        }
        let start = self.source_offsets[prefix];
        let end = self.source_offsets[prefix + 1];
        self.sources[start..end]
            .binary_search_by_key(&barcode, |&(source, _)| source)
            .ok()
            .map(|index| self.sources[start + index].1)
    }
}

fn merge_unique_target<T: Copy + Eq>(target: &mut Option<T>, candidate: T) -> bool {
    match *target {
        Some(previous) => previous != candidate,
        None => {
            *target = Some(candidate);
            false
        }
    }
}

fn validate_barcode(barcode: u64, barcode_len: u8) -> anyhow::Result<()> {
    if barcode_len < 32 && barcode >= (1u64 << (2 * barcode_len)) {
        bail!("packed barcode {barcode} does not fit declared length {barcode_len}");
    }
    Ok(())
}

fn for_each_neighbor(
    barcode: u64,
    barcode_len: u8,
    neighborhood: BarcodeNeighborhood,
    mut visit: impl FnMut(u64),
) {
    for_each_substitution(barcode, barcode_len, &mut visit);
    if neighborhood == BarcodeNeighborhood::SubstitutionOrShiftOne {
        for_each_shift_neighbor(barcode, barcode_len, visit);
    }
}

fn for_each_substitution(barcode: u64, barcode_len: u8, mut visit: impl FnMut(u64)) {
    for position in 0..usize::from(barcode_len) {
        let shift = 2 * position;
        let mask = 3u64 << shift;
        let observed_base = (barcode & mask) >> shift;
        let cleared = barcode & !mask;
        for replacement in 0..4u64 {
            if replacement != observed_base {
                visit(cleared | (replacement << shift));
            }
        }
    }
}

fn for_each_shift_neighbor(barcode: u64, barcode_len: u8, mut visit: impl FnMut(u64)) {
    // Historical fixed-length insertion/deletion constructions. Each pair is
    // a shift with an arbitrary admitted terminal base.
    for boundary in 1..usize::from(barcode_len) {
        let lower_mask = (1u64 << (2 * boundary)) - 1;
        let upper = barcode & !lower_mask;
        let lower = barcode & lower_mask;
        for admitted in 0..4u64 {
            let insertion = upper | (admitted << (2 * (boundary - 1))) | (lower >> 2);
            let deletion = upper | admitted | ((lower & !(3u64 << (2 * boundary))) << 2);
            if insertion != barcode {
                visit(insertion);
            }
            if deletion != barcode {
                visit(deletion);
            }
        }
    }
}

/// Enumerate retained sources which would generate `observed` under the
/// historical directed shift construction above. Computing the inverse per
/// observation avoids materializing O(retained * barcode_len) reverse pairs.
fn for_each_inverse_shift_candidate(observed: u64, barcode_len: u8, mut visit: impl FnMut(u64)) {
    for boundary in 1..usize::from(barcode_len) {
        let lower_mask = packed_base_mask(boundary);
        let lower_without_terminal = packed_base_mask(boundary - 1);

        // Inverse of `insertion`: the observed low (boundary - 1) bases
        // determine source bases 1..boundary, while source base 0 is free.
        let insertion_fixed = (observed & !lower_mask) | ((observed & lower_without_terminal) << 2);
        for terminal in 0..4u64 {
            visit(insertion_fixed | terminal);
        }

        // Inverse of historical `deletion`. Its OR at the boundary means the
        // two contributing source bases must be enumerated explicitly.
        let above_boundary_mask = !packed_base_mask(boundary + 1);
        let deletion_fixed =
            (observed & above_boundary_mask) | ((observed >> 2) & lower_without_terminal);
        let observed_boundary = (observed >> (2 * boundary)) & 3;
        for lower in 0..4u64 {
            for upper in 0..4u64 {
                if lower | upper == observed_boundary {
                    visit(
                        deletion_fixed
                            | (lower << (2 * (boundary - 1)))
                            | (upper << (2 * boundary)),
                    );
                }
            }
        }
    }
}

fn packed_base_mask(num_bases: usize) -> u64 {
    if num_bases == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * num_bases)) - 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn spec(neighborhood: BarcodeNeighborhood, resolution: BarcodeResolution) -> CorrectionSpec {
        CorrectionSpec {
            barcode_len: 2,
            neighborhood,
            resolution,
        }
    }

    #[test]
    fn confidence_parsing_is_exact_and_reduced() {
        assert_eq!("0.975".parse::<Confidence>().unwrap(), Confidence::RNA);
        assert_eq!("90/100".parse::<Confidence>().unwrap(), Confidence::ATAC);
        assert!("1.01".parse::<Confidence>().is_err());
        assert!("NaN".parse::<Confidence>().is_err());
        assert_eq!(
            compare_fractions(u128::MAX - 1, u128::MAX, 39, 40),
            std::cmp::Ordering::Greater
        );
    }

    #[test]
    fn exact_sources_win_over_neighbor_collisions() {
        let index = CorrectionIndex::new(
            spec(BarcodeNeighborhood::HammingOne, BarcodeResolution::Unique),
            [
                RetainedSource {
                    source: 0,
                    target: 10,
                    exact_count: 1,
                },
                RetainedSource {
                    source: 1,
                    target: 11,
                    exact_count: 1,
                },
            ],
        )
        .unwrap();
        assert_eq!(index.resolve(0), CorrectionDecision::Exact(10));
    }

    #[test]
    fn unique_counts_distinct_canonical_targets() {
        let index = CorrectionIndex::new(
            spec(BarcodeNeighborhood::HammingOne, BarcodeResolution::Unique),
            [
                RetainedSource {
                    source: 1,
                    target: 7,
                    exact_count: 4,
                },
                RetainedSource {
                    source: 2,
                    target: 7,
                    exact_count: 5,
                },
            ],
        )
        .unwrap();
        assert_eq!(index.resolve(0), CorrectionDecision::Corrected(7));
    }

    #[test]
    fn unique_rejects_multiple_targets_independent_of_input_order() {
        let retained = [
            RetainedSource {
                source: 1,
                target: 10,
                exact_count: 1,
            },
            RetainedSource {
                source: 2,
                target: 20,
                exact_count: 1,
            },
        ];
        let forward = CorrectionIndex::new(
            spec(BarcodeNeighborhood::HammingOne, BarcodeResolution::Unique),
            retained,
        )
        .unwrap();
        let reverse = CorrectionIndex::new(
            spec(BarcodeNeighborhood::HammingOne, BarcodeResolution::Unique),
            retained.into_iter().rev(),
        )
        .unwrap();
        assert_eq!(forward.resolve(0), CorrectionDecision::Ambiguous);
        assert_eq!(forward.resolve(0), reverse.resolve(0));
    }

    #[test]
    fn structural_compilation_separates_unique_and_ambiguous_routes() {
        let index = CorrectionIndex::new(
            CorrectionSpec {
                barcode_len: 1,
                neighborhood: BarcodeNeighborhood::HammingOne,
                resolution: BarcodeResolution::Frequency {
                    confidence: Confidence::RNA,
                    pseudocount: 1,
                },
            },
            [
                RetainedSource {
                    source: 0,
                    target: 10,
                    exact_count: 100,
                },
                RetainedSource {
                    source: 2,
                    target: 20,
                    exact_count: 1,
                },
            ],
        )
        .unwrap();
        let topology = index.compile_structural_neighborhood();
        assert_eq!(topology.entries(), &[(0, 10), (2, 20)]);
        assert_eq!(topology.ambiguous(), &[1, 3]);
        // The structural table deliberately remains policy-free even though
        // Frequency can accept the high-prior target for an ambiguous source.
        assert_eq!(index.resolve(1), CorrectionDecision::Corrected(10));
    }

    #[test]
    fn frequency_sums_alias_weights_and_uses_exact_threshold() {
        let resolution = BarcodeResolution::Frequency {
            confidence: Confidence::RNA,
            pseudocount: 1,
        };
        let accepted = CorrectionIndex::new(
            spec(BarcodeNeighborhood::HammingOne, resolution),
            [
                RetainedSource {
                    source: 1,
                    target: 10,
                    exact_count: 19,
                },
                RetainedSource {
                    source: 2,
                    target: 10,
                    exact_count: 18,
                },
                RetainedSource {
                    source: 3,
                    target: 20,
                    exact_count: 0,
                },
            ],
        )
        .unwrap();
        assert_eq!(accepted.resolve(0), CorrectionDecision::Corrected(10));

        let rejected = CorrectionIndex::new(
            spec(BarcodeNeighborhood::HammingOne, resolution),
            [
                RetainedSource {
                    source: 1,
                    target: 10,
                    exact_count: 37,
                },
                RetainedSource {
                    source: 2,
                    target: 20,
                    exact_count: 1,
                },
            ],
        )
        .unwrap();
        assert_eq!(rejected.resolve(0), CorrectionDecision::Ambiguous);
    }

    #[test]
    fn frequency_tie_break_is_deterministic_but_still_subject_to_confidence() {
        let index = CorrectionIndex::new(
            spec(
                BarcodeNeighborhood::HammingOne,
                BarcodeResolution::Frequency {
                    confidence: Confidence::new(1, 2).unwrap(),
                    pseudocount: 1,
                },
            ),
            [
                RetainedSource {
                    source: 1,
                    target: 10,
                    exact_count: 0,
                },
                RetainedSource {
                    source: 2,
                    target: 20,
                    exact_count: 0,
                },
            ],
        )
        .unwrap();
        assert_eq!(index.resolve(0), CorrectionDecision::Corrected(20));
    }

    #[test]
    fn shift_candidates_are_deduplicated_before_scoring() {
        // Both shift construction directions can produce the same source for
        // repetitive barcodes; it must contribute its prior only once.
        let index = CorrectionIndex::new(
            spec(
                BarcodeNeighborhood::SubstitutionOrShiftOne,
                BarcodeResolution::Unique,
            ),
            [RetainedSource {
                source: 5,
                target: 9,
                exact_count: 0,
            }],
        )
        .unwrap();
        assert_eq!(index.resolve(1), CorrectionDecision::Corrected(9));
    }

    #[test]
    fn observed_compilation_is_sorted_and_aggregates_stats() {
        let index = CorrectionIndex::new(
            spec(BarcodeNeighborhood::HammingOne, BarcodeResolution::Unique),
            [RetainedSource {
                source: 1,
                target: 1,
                exact_count: 4,
            }],
        )
        .unwrap();
        let table = index.compile_observed([(8, 3), (0, 2), (0, 4), (1, 4)]);
        assert_eq!(table.entries(), &[(0, 1), (1, 1)]);
        assert_eq!(table.stats().exact_reads, 4);
        assert_eq!(table.stats().corrected_reads, 6);
        assert_eq!(table.stats().not_found_reads, 3);
    }

    #[test]
    fn distinct_compilation_adds_unobserved_retained_identities() {
        let index = CorrectionIndex::new(
            spec(BarcodeNeighborhood::HammingOne, BarcodeResolution::Unique),
            [RetainedSource {
                source: 1,
                target: 7,
                exact_count: 0,
            }],
        )
        .unwrap();
        let (table, target_counts) = index.compile_distinct_observed_with_target_counts([(8, 3)]);
        assert_eq!(table.entries(), &[(1, 7)]);
        assert_eq!(table.stats().exact_distinct, 1);
        assert_eq!(table.stats().exact_reads, 0);
        assert_eq!(table.stats().not_found_reads, 3);
        assert!(target_counts.is_empty());
    }

    #[test]
    fn all_supported_lengths_validate_without_shift_overflow() {
        for barcode_len in 1..=32 {
            let index = CorrectionIndex::new(
                CorrectionSpec {
                    barcode_len,
                    neighborhood: BarcodeNeighborhood::SubstitutionOrShiftOne,
                    resolution: BarcodeResolution::Unique,
                },
                [RetainedSource {
                    source: u64::MAX >> (64 - 2 * usize::from(barcode_len)),
                    target: 1,
                    exact_count: 0,
                }],
            )
            .unwrap();
            let _ = index.resolve(0);
        }
    }

    #[test]
    fn shift_neighborhood_matches_the_historical_construction() {
        for barcode_len in [1usize, 2, 8, 16, 31] {
            let mask = (1u64 << (2 * barcode_len)) - 1;
            for barcode in [0, mask, mask / 3, mask / 5] {
                let expected: HashSet<_> = crate::utils::get_all_snps(barcode, barcode_len)
                    .into_iter()
                    .chain(crate::utils::get_all_indels(barcode, barcode_len))
                    .collect();
                let mut actual = HashSet::new();
                for_each_neighbor(
                    barcode,
                    barcode_len as u8,
                    BarcodeNeighborhood::SubstitutionOrShiftOne,
                    |neighbor| {
                        actual.insert(neighbor);
                    },
                );
                assert_eq!(actual, expected, "length {barcode_len}, barcode {barcode}");
            }
        }
    }

    #[test]
    fn inverse_shift_candidates_match_the_directed_forward_relation() {
        for barcode_len in 1u8..=4 {
            let universe = 1u64 << (2 * barcode_len);
            for observed in 0..universe {
                let mut expected = HashSet::new();
                for source in 0..universe {
                    for_each_shift_neighbor(source, barcode_len, |generated| {
                        if generated == observed {
                            expected.insert(source);
                        }
                    });
                }
                let mut actual = HashSet::new();
                for_each_inverse_shift_candidate(observed, barcode_len, |source| {
                    if source != observed {
                        actual.insert(source);
                    }
                });
                assert_eq!(
                    actual, expected,
                    "length {barcode_len}, observed {observed}"
                );
            }
        }
    }

    #[test]
    fn directed_shift_lookup_preserves_historical_source_to_observation_semantics() {
        let spec = CorrectionSpec {
            barcode_len: 3,
            neighborhood: BarcodeNeighborhood::SubstitutionOrShiftOne,
            resolution: BarcodeResolution::Unique,
        };
        let forward = CorrectionIndex::new(
            spec,
            [RetainedSource {
                source: 1,
                target: 1,
                exact_count: 0,
            }],
        )
        .unwrap();
        assert_eq!(forward.resolve(8), CorrectionDecision::Corrected(1));

        let reverse = CorrectionIndex::new(
            spec,
            [RetainedSource {
                source: 8,
                target: 8,
                exact_count: 0,
            }],
        )
        .unwrap();
        assert_eq!(reverse.resolve(1), CorrectionDecision::NotFound);
    }
}
