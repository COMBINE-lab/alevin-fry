
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

/// Sequences up to this length use fixed-size stack arrays -- zero heap
/// allocation. 32 comfortably covers a 12nt UMI with headroom.
pub const MAX_UMI_LEN: usize = 32;

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
        assert_eq!(umi_edit_distance_from_packed(x, x, 12), 0);
    }

    #[test]
    fn shift_right_by_one_costs_one() {
        let x = pack("ACGTACGTACGT");
        let y = pack("CGTACGTACGTA");
        assert_eq!(umi_edit_distance_from_packed(x, y, 12), 1);
    }

    #[test]
    fn shift_right_by_two_costs_two() {
        let x = pack("ACGTACGTACGT");
        let y = pack("GTACGTACGTAC");
        assert_eq!(umi_edit_distance_from_packed(x, y, 12), 2);
    }

    #[test]
    fn middle_substitution_costs_one() {
        let x = pack("ACGTACGTACGT");
        let y = pack("ACGTAAGTACGT");
        assert_eq!(umi_edit_distance_from_packed(x, y, 12), 1);
    }

    #[test]
    fn disputed_pair_is_two() {
        let x = pack("CCTGTCACGTGT");
        let y = pack("CGTCACGTGTAC");
        assert_eq!(umi_edit_distance_from_packed(x, y, 12), 2);
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
                umi_edit_distance_from_packed(x, y, 12),
                umi_edit_distance_from_packed(y, x, 12)
            );
        }
    }
}