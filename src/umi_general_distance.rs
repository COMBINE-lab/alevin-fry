
//! Shift-tolerant edit distance for 2-bit-packed UMIs.
//! ================================================================
//!
//! Rules: deletion is free only in the true leading run (nothing of the
//! other UMI produced yet) or true trailing run (all of it already
//! produced); deletion elsewhere, insertion, and substitution all cost 1.
//! The distance is symmetrized by taking the minimum of both directions,
//! so argument order never matters.
//!
//! Operates directly on the 2-bit-packed `u64` representation -- bases are
//! extracted with a shift-and-mask rather than decoded to ASCII first, and
//! both DP rows are fixed-size stack arrays, so a single distance call
//! performs zero heap allocation.
//!
//! Bit layout: matches needletail's `BitKmer` convention -- the *first*
//! base of the sequence occupies the *highest* 2 bits of the `len`-base
//! window, and each subsequent base sits progressively closer to the
//! lowest 2 bits (needletail builds this up via `(kmer << 2) + new_base`
//! as each base is appended, so the earliest-appended base ends up
//! shifted furthest to the high end). `base_at` accounts for this by
//! reading from bit position `2 * (len - 1 - idx)` for a given `idx`, so
//! `idx = 0` correctly means "the first base," matching the true sequence
//! order rather than the raw bit order.
//!
//! This is a shift-tolerant distance (appropriate for long-read /
//! indel-prone UMI errors). It intentionally does *not* behave like a
//! plain Hamming distance: a UMI that's simply shifted by one or two
//! bases will score much closer here than it would under Hamming
//! distance, which has implications for any threshold applied to the
//! result (see caller).
//!
//! Bounded evaluation
//! ------------------
//! The deduplication only needs to know the distance when it is small
//! (at most [`BOUND`]); every larger distance is treated alike. So the
//! callers use [`umi_edit_distance_bounded`], which returns `None` for
//! any pair farther than `BOUND` apart, and which gets its speed from two
//! exact filters in front of the full dynamic program:
//!
//! 1. a *q-gram pigeonhole* test: cutting `y` into `BOUND + 1` pieces, a
//!    pair within distance `BOUND` must have at least one piece appear
//!    verbatim in `x` within `BOUND` positions of where it sits in `y`,
//!    because every unit of offset costs an insertion (or a non-free
//!    deletion) and every edit inside a piece breaks it;
//! 2. a *banded* dynamic program (cells with `|i - j| > BOUND` are never
//!    on a path of cost at most `BOUND`, as a shift of `s` costs at least
//!    `s` insertions) with an early exit as soon as every cell of a row
//!    exceeds the bound.
//!
//! Both filters are necessary conditions only, so the bounded distance is
//! exactly `min(d, BOUND + 1)` of the unbounded one: the deduplication's
//! decisions are unchanged.

/// Sequences up to this length use fixed-size stack arrays -- zero heap
/// allocation. 32 comfortably covers a 12nt UMI with headroom.
pub const MAX_UMI_LEN: usize = 32;

/// The largest distance the deduplication distinguishes; beyond it every
/// pair is "different".
pub const BOUND: u32 = 3;

/// Extract the 2-bit code for base `idx` (0-indexed, 0 = first base of the
/// sequence) out of a packed UMI of `len` bases, matching needletail's
/// `BitKmer` layout: the first base sits in the highest 2 bits of the
/// `len`-base window, the last base in the lowest 2 bits.
#[inline]
fn base_at(packed: u64, idx: usize, len: usize) -> u8 {
    ((packed >> (2 * (len - 1 - idx))) & 0b11) as u8
}

#[inline]
fn del_cost(j: usize, n: usize) -> u32 {
    if j == 0 || j == n {
        0
    } else {
        1
    }
}

/// Directional distance on packed UMIs: deletion from `x` is free only in
/// the true leading/trailing runs, costs 1 in the middle; insertion and
/// substitution always cost 1.
fn umi_distance_directional_packed(x: u64, y: u64, len: usize) -> u32 {
    let mut prev = [0u32; MAX_UMI_LEN + 1];
    let mut curr = [0u32; MAX_UMI_LEN + 1];

    for j in 0..=len {
        prev[j] = j as u32;
    }

    for i in 1..=len {
        curr[0] = 0;
        let xi = base_at(x, i - 1, len);
        for j in 1..=len {
            let yj = base_at(y, j - 1, len);
            let sub_cost = if xi == yj { 0 } else { 1 };
            let diag = prev[j - 1] + sub_cost;
            let up = prev[j] + del_cost(j, len);
            let left = curr[j - 1] + 1;
            curr[j] = diag.min(up).min(left);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev[len]
}

/// The directional distance, computed only within the band `|i - j| <= bound`
/// and only until it is known to exceed `bound`: returns `None` in that case,
/// otherwise the exact distance (which is then at most `bound`).
///
/// Why the band is exact: a path through cell `(i, j)` has produced `j` bases
/// of `y` from `i` bases of `x`. If `i > j`, at least `i - j` bases of `x` were
/// deleted; a free deletion is only possible in the leading run (then `j = 0`
/// and the remaining `len - i < len - j` bases of `x` can produce at most
/// `len - i` of `y`'s `len` bases, so at least `i` insertions follow) or the
/// trailing run (then everything of `y` is already produced). Either way the
/// path costs at least `i - j`. If `j > i`, at least `j - i` insertions
/// happened. So no cell outside the band lies on a path of cost at most
/// `bound`, and treating those cells as "more than bound" changes nothing.
#[inline]
fn umi_distance_directional_bounded(x: u64, y: u64, len: usize, bound: u32) -> Option<u32> {
    let over = bound + 1;
    let b = bound as usize;
    let mut prev = [0u32; MAX_UMI_LEN + 1];
    let mut curr = [0u32; MAX_UMI_LEN + 1];

    for j in 0..=len {
        prev[j] = if j <= b { j as u32 } else { over };
    }

    for i in 1..=len {
        let lo = i.saturating_sub(b).max(1);
        let hi = (i + b).min(len);
        // cells left of the band: more than the bound (a shift of j < i - b costs > b)
        if lo > 1 {
            curr[lo - 1] = over;
        } else {
            curr[0] = 0;
        }
        let xi = base_at(x, i - 1, len);
        let mut row_min = over;
        for j in lo..=hi {
            let yj = base_at(y, j - 1, len);
            let sub_cost = if xi == yj { 0 } else { 1 };
            let diag = prev[j - 1] + sub_cost;
            let up = if j < i + b { prev[j] + del_cost(j, len) } else { over };
            let left = curr[j - 1] + 1;
            let v = diag.min(up).min(left).min(over);
            curr[j] = v;
            row_min = row_min.min(v);
        }
        if hi < len {
            curr[hi + 1] = over;
        }
        // a free trailing deletion keeps a finished row alive: cell (i, len) is
        // reachable from (i - 1, len) at no cost, and it is inside the band only
        // when len <= i + b, which the loop above already covers
        if row_min > bound {
            return None;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    let d = prev[len];
    if d > bound { None } else { Some(d) }
}

/// Whether some piece of `y` (cut into `bound + 1` pieces) appears verbatim
/// in `x` within `bound` positions of its place in `y`.
///
/// In the direction where `x`'s leading and trailing bases may be deleted
/// for free, a path of cost at most `bound` spoils at most `bound` pieces of
/// `y`: every substitution or insertion lands on a base of `y`, every paid
/// deletion lands between two bases of `y`, and the free deletions lie
/// outside `y`'s span. So one piece is matched base for base to consecutive
/// bases of `x`, and its offset is bounded by the cost: `s` free leading
/// deletions and `d` paid deletions before the piece leave `len - s - d`
/// bases of `x` to produce the `len` bases of `y`, which takes at least
/// `s + d` insertions.
#[inline]
fn some_piece_matches(x: u64, y: u64, len: usize, bound: usize) -> bool {
    let pieces = bound + 1;
    if len < pieces {
        return true;
    }
    let piece_len = len / pieces;
    let mask: u64 = (1u64 << (2 * piece_len)) - 1;
    // piece p of y covers bases [p * piece_len, (p + 1) * piece_len)
    for p in 0..pieces {
        let start = p * piece_len;
        // the window ends at base start + piece_len - 1, whose low bit sits at 2 * (len - 1 - (start + piece_len - 1))
        let shift_y = 2 * (len - start - piece_len);
        let piece = (y >> shift_y) & mask;
        let from = start.saturating_sub(bound);
        let to = (start + bound).min(len - piece_len);
        for s in from..=to {
            let shift_x = 2 * (len - s - piece_len);
            if (x >> shift_x) & mask == piece {
                return true;
            }
        }
    }
    false
}

/// A lower bound on the distance from base composition alone: every base of the
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

/// The pigeonhole filter for the symmetrized distance: `true` when the pair
/// may be within the bound (the full test must follow), `false` when it
/// certainly is not. The symmetrized distance is the smaller of the two
/// directional distances; in the direction with free deletions on `x`, a
/// small distance leaves a piece of `y` in `x` ([`some_piece_matches`]), and
/// in the other direction a piece of `x` in `y`. The free deletions of the
/// string that is cut into pieces can spoil pieces at no cost, which is why
/// the pieces of the *other* string are what each direction guarantees, and
/// why both directions are tested.
#[inline]
fn may_be_within_bound(x: u64, y: u64, len: usize, bound: usize) -> bool {
    some_piece_matches(x, y, len, bound) || some_piece_matches(y, x, len, bound)
}

/// Shift-tolerant UMI distance on packed representations. Symmetrized by
/// taking the minimum of both directions, so argument order never
/// matters. Drop-in replacement for `afutils::count_diff_2_bit_packed`.
///
/// `umi_len` is the number of bases packed into `x` and `y` (both must be
/// packed with the same length, using needletail's `BitKmer` convention).
#[inline]
pub fn umi_edit_distance_from_packed_shifted(x: u64, y: u64, umi_len: u8) -> usize {
    let len = umi_len as usize;
    debug_assert!(
        len <= MAX_UMI_LEN,
        "umi_len exceeds MAX_UMI_LEN; increase the constant"
    );
    let fwd = umi_distance_directional_packed(x, y, len);
    let bwd = umi_distance_directional_packed(y, x, len);
    fwd.min(bwd) as usize
}

/// The shift-tolerant distance when it is at most `bound`, `None` otherwise.
/// Equal to `umi_edit_distance_from_packed_shifted` wherever that is at most
/// `bound`, at a fraction of its cost: a pigeonhole test rejects most random
/// pairs without any dynamic programming, and the rest run a banded program
/// that stops as soon as the bound is exceeded.
#[inline]
pub fn umi_edit_distance_within(x: u64, y: u64, umi_len: u8, bound: u32) -> Option<u32> {
    let len = umi_len as usize;
    debug_assert!(
        len <= MAX_UMI_LEN,
        "umi_len exceeds MAX_UMI_LEN; increase the constant"
    );
    if x == y {
        return Some(0);
    }
    if composition_lower_bound(x, y, len) > bound {
        return None;
    }
    if !may_be_within_bound(x, y, len, bound as usize) {
        return None;
    }
    let fwd = umi_distance_directional_bounded(x, y, len, bound);
    let bwd = umi_distance_directional_bounded(y, x, len, bound);
    match (fwd, bwd) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) | (None, Some(a)) => Some(a),
        (None, None) => None,
    }
}

/// [`umi_edit_distance_within`] at the deduplication's bound [`BOUND`].
#[inline]
pub fn umi_edit_distance_bounded(x: u64, y: u64, umi_len: u8) -> Option<u32> {
    umi_edit_distance_within(x, y, umi_len, BOUND)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test-only helper: pack an ACGT string into a u64 using needletail's
    /// actual convention -- `(packed << 2) | code` per base, appended in
    /// sequence order, so the first base ends up in the high bits and the
    /// last base in the low bits (matching `extend_kmer` in
    /// needletail::bitkmer).
    fn pack(seq: &str) -> u64 {
        let mut v: u64 = 0;
        for c in seq.bytes() {
            let code: u64 = match c {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("bad base"),
            };
            v = (v << 2) | code;
        }
        v
    }

    #[test]
    fn pack_matches_needletail_test_vectors() {
        // Sanity check against needletail's own test cases for
        // bytes_to_bitmer, to confirm this test helper matches the real
        // library's convention.
        assert_eq!(pack("C"), 1);
        assert_eq!(pack("TTA"), 60);
        assert_eq!(pack("AAA"), 0);
    }

    #[test]
    fn base_at_reads_first_base_as_idx_zero() {
        let packed = pack("ACGT");
        assert_eq!(base_at(packed, 0, 4), 0); // 'A'
        assert_eq!(base_at(packed, 1, 4), 1); // 'C'
        assert_eq!(base_at(packed, 2, 4), 2); // 'G'
        assert_eq!(base_at(packed, 3, 4), 3); // 'T'
    }

    #[test]
    fn identical_sequences_are_zero() {
        let x = pack("ACGTACGTACGT");
        assert_eq!(umi_edit_distance_from_packed_shifted(x, x, 12), 0);
    }

    #[test]
    fn shift_right_by_one_costs_one() {
        let x = pack("ACGTACGTACGT");
        let y = pack("CGTACGTACGTA");
        assert_eq!(umi_edit_distance_from_packed_shifted(x, y, 12), 1);
    }

    #[test]
    fn shift_right_by_two_costs_two() {
        let x = pack("ACGTACGTACGT");
        let y = pack("GTACGTACGTAC");
        assert_eq!(umi_edit_distance_from_packed_shifted(x, y, 12), 2);
    }

    #[test]
    fn middle_substitution_costs_one() {
        let x = pack("ACGTACGTACGT");
        let y = pack("ACGTAAGTACGT");
        assert_eq!(umi_edit_distance_from_packed_shifted(x, y, 12), 1);
    }

    #[test]
    fn disputed_pair_is_two() {
        let x = pack("CCTGTCACGTGT");
        let y = pack("CGTCACGTGTAC");
        assert_eq!(umi_edit_distance_from_packed_shifted(x, y, 12), 2);
    }

    #[test]
    fn distance_is_symmetric_regardless_of_argument_order() {
        let pairs = [
            ("ACGTACGTACGT", "ACGTACGTACGT"),
            ("ACGTACGTACGT", "CGTACGTACGTA"),
            ("ACGTACGTACGT", "GTACGTACGTAC"),
            ("ACGTACGTACGT", "ACGTAAGTACGT"),
            ("CCTGTCACGTGT", "CGTCACGTGTAC"),
        ];
        for (a, b) in pairs {
            let x = pack(a);
            let y = pack(b);
            assert_eq!(
                umi_edit_distance_from_packed_shifted(x, y, 12),
                umi_edit_distance_from_packed_shifted(y, x, 12)
            );
        }
    }

    /// A small deterministic generator for exhaustive-ish cross-checks.
    fn xorshift(state: &mut u64) -> u64 {
        let mut x = *state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *state = x;
        x
    }

    #[test]
    fn bounded_distance_agrees_with_the_full_one() {
        // random pairs and near pairs (a few edits of a random UMI), for UMI lengths 8 to 16
        let mut st = 0x9E37_79B9_7F4A_7C15u64;
        for len in [6usize, 8, 10, 12, 16] {
            let mask = if len == 32 { u64::MAX } else { (1u64 << (2 * len)) - 1 };
            for _ in 0..20_000 {
                let x = xorshift(&mut st) & mask;
                let y = match xorshift(&mut st) % 4 {
                    0 => xorshift(&mut st) & mask,
                    // substitutions
                    1 => {
                        let mut y = x;
                        for _ in 0..(1 + xorshift(&mut st) % 3) {
                            let pos = (xorshift(&mut st) % len as u64) as usize;
                            let base = xorshift(&mut st) & 3;
                            y = (y & !(3u64 << (2 * (len - 1 - pos)))) | (base << (2 * (len - 1 - pos)));
                        }
                        y
                    }
                    // shifts with new bases entering
                    2 => {
                        let s = 1 + (xorshift(&mut st) % 3) as usize;
                        ((x << (2 * s)) | (xorshift(&mut st) & ((1u64 << (2 * s)) - 1))) & mask
                    }
                    _ => {
                        let s = 1 + (xorshift(&mut st) % 3) as usize;
                        ((x >> (2 * s)) | ((xorshift(&mut st) & ((1u64 << (2 * s)) - 1)) << (2 * (len - s)))) & mask
                    }
                };
                let full = umi_edit_distance_from_packed_shifted(x, y, len as u8) as u32;
                for bound in 0..=5u32 {
                    let bounded = umi_edit_distance_within(x, y, len as u8, bound);
                    let expected = if full <= bound { Some(full) } else { None };
                    assert_eq!(bounded, expected, "len {len} bound {bound} x {x:#x} y {y:#x}: full {full}, bounded {bounded:?}");
                }
            }
        }
    }

    #[test]
    fn bounded_distance_exhaustive_short_umis() {
        // every pair of 5-mers: 1024 x 1024
        let len = 5usize;
        for x in 0..(1u64 << (2 * len)) {
            for y in 0..(1u64 << (2 * len)) {
                let full = umi_edit_distance_from_packed_shifted(x, y, len as u8) as u32;
                let bounded = umi_edit_distance_bounded(x, y, len as u8);
                let expected = if full <= BOUND { Some(full) } else { None };
                assert_eq!(bounded, expected, "x {x:#x} y {y:#x}");
            }
        }
    }
}
