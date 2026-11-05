"""Step 7: (a) per-molecule posteriors over the molecule's own candidates instead of a dense vector over
every transcript of the cell; (b) a composition lower bound in front of the UMI distance filter.

(a) The write-back at the end of a cell's truncation fit allocated and zeroed one f64 per transcript of
the cell for every molecule and then scanned it all (tens of thousands of entries per molecule, tens of
thousands of molecules per cell). The shares are now accumulated over the molecule's own alignments,
in ascending transcript order with repeated transcripts summed in alignment order, which gives
exactly the same values and the same ordering as before.
(b) Every base of the produced UMI is matched exactly, substituted or inserted, so the distance is at
least the number of bases that cannot be matched: len - sum over bases of min(count_x, count_y).
"""
import os
P = "/scratch1/zahra/source-code/alevin-fry/src/truncation.rs"
s = open(P).read()

def rep(old, new, count=1):
    global s
    n = s.count(old); assert n == count, f"expected {count} of {old[:70]!r}, found {n}"
    s = s.replace(old, new)

rep("""    let exps = Exps::new(&par, max_len);
    for (m, molecule) in molecules.iter_mut().enumerate() {
        let alignments = cell.molecule(m);
        let mut weights: Vec<f64> = match &frozen {
            Some(fx) => (fx.off[m]..fx.off[m + 1])
                .map(|k| prev[fx.t[k]] * fx.w[k])
                .collect(),
            None => alignments
                .iter()
                .map(|a| a.weight(&prev, &par, &exps, cell.positional).0)
                .collect(),
        };
        let mut denominator: f64 = weights.iter().sum();
        if !(denominator.is_finite() && denominator > EM_DENOM_THRESH) {
            // The molecule is already deduplicated and must keep its counts, so it falls
            // back to its prior probabilities.
            weights = alignments.iter().map(|a| a.prior).collect();
            denominator = weights.iter().sum();
        }
        assert!(
            denominator.is_finite() && denominator > 0.0,
            "truncation model produced invalid final molecule probabilities"
        );
        let mut probabilities = vec![0.0; n_targets];
        for (alignment, weight) in alignments.iter().zip(weights) {
            probabilities[alignment.target] += weight / denominator;
        }
        prune_shares(&mut probabilities, spec.min_share);
        sample_small_shares(&mut probabilities, spec.sample_below, molecule_uniform(&molecule.vertices));
        molecule.transcript_probs = probabilities
            .into_iter()
            .enumerate()
            .filter_map(|(target, probability)| {
                (probability > 0.0).then_some((target_ids[target], probability))
            })
            .collect();
    }
""",
"""    let exps = Exps::new(&par, max_len);
    // per-molecule buffers, reused: the weights of its alignments, and its shares over its own
    // candidate transcripts in ascending target order (a transcript with several alignments has
    // them summed in alignment order, as the dense vector did)
    let mut weights: Vec<f64> = Vec::new();
    let mut shares: Vec<(usize, f64)> = Vec::new();
    let mut targets: Vec<usize> = Vec::new();
    let mut probabilities: Vec<f64> = Vec::new();
    for (m, molecule) in molecules.iter_mut().enumerate() {
        let alignments = cell.molecule(m);
        weights.clear();
        match &frozen {
            Some(fx) => weights.extend((fx.off[m]..fx.off[m + 1]).map(|k| prev[fx.t[k]] * fx.w[k])),
            None => weights.extend(
                alignments
                    .iter()
                    .map(|a| a.weight(&prev, &par, &exps, cell.positional).0),
            ),
        };
        let mut denominator: f64 = weights.iter().sum();
        if !(denominator.is_finite() && denominator > EM_DENOM_THRESH) {
            // The molecule is already deduplicated and must keep its counts, so it falls
            // back to its prior probabilities.
            weights.clear();
            weights.extend(alignments.iter().map(|a| a.prior));
            denominator = weights.iter().sum();
        }
        assert!(
            denominator.is_finite() && denominator > 0.0,
            "truncation model produced invalid final molecule probabilities"
        );
        shares.clear();
        shares.extend(
            alignments
                .iter()
                .zip(&weights)
                .map(|(alignment, weight)| (alignment.target, weight / denominator)),
        );
        // stable: repeated targets keep their alignment order
        shares.sort_by_key(|&(target, _)| target);
        targets.clear();
        probabilities.clear();
        for &(target, share) in &shares {
            if targets.last() == Some(&target) {
                *probabilities.last_mut().expect("a share") += share;
            } else {
                targets.push(target);
                probabilities.push(share);
            }
        }
        prune_shares(&mut probabilities, spec.min_share);
        sample_small_shares(&mut probabilities, spec.sample_below, molecule_uniform(&molecule.vertices));
        molecule.transcript_probs = targets
            .iter()
            .zip(&probabilities)
            .filter(|(_, &probability)| probability > 0.0)
            .map(|(&target, &probability)| (target_ids[target], probability))
            .collect();
    }
""")
open(P, "w").write(s)
print("step 7a: write-back over the molecule's candidates")

# ---- the composition bound in the distance module
P2 = "/scratch1/zahra/source-code/alevin-fry/src/umi_general_distance.rs"
s = open(P2).read()
rep("""/// The pigeonhole filter for the symmetrized distance:""",
"""/// A lower bound on the distance from base composition alone: every base of the
/// produced UMI is matched exactly, substituted or inserted, and an exact match
/// pairs a base of one UMI with an equal base of the other, so at most
/// `sum over bases of min(count_x, count_y)` bases are matched and the rest cost
/// one each. Holds in both directions (free deletions only reduce what can be
/// matched), so it bounds the symmetrized distance.
#[inline]
fn composition_lower_bound(x: u64, y: u64, len: usize) -> u32 {
    let mask: u64 = if len >= 32 { u64::MAX } else { (1u64 << (2 * len)) - 1 };
    let (x, y) = (x & mask, y & mask);
    const LOW: u64 = 0x5555_5555_5555_5555;
    let (xh, xl) = ((x >> 1) & LOW, x & LOW);
    let (yh, yl) = ((y >> 1) & LOW, y & LOW);
    // per-base counts: A = 00, C = 01, G = 10, T = 11
    let count = |h: u64, l: u64| -> [u32; 4] {
        let t = (h & l).count_ones();
        let g = (h & !l).count_ones();
        let c = (!h & l & LOW).count_ones();
        [len as u32 - t - g - c, c, g, t]
    };
    let (cx, cy) = (count(xh, xl), count(yh, yl));
    let matched: u32 = (0..4).map(|b| cx[b].min(cy[b])).sum();
    len as u32 - matched
}

/// The pigeonhole filter for the symmetrized distance:""")
rep("""    if x == y {
        return Some(0);
    }
    if !may_be_within_bound(x, y, len, bound as usize) {
        return None;
    }""",
"""    if x == y {
        return Some(0);
    }
    if composition_lower_bound(x, y, len) > bound {
        return None;
    }
    if !may_be_within_bound(x, y, len, bound as usize) {
        return None;
    }""")
open(P2, "w").write(s)
print("step 7b: composition bound in the distance filter")
