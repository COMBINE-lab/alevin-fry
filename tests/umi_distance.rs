//! The general (shift-tolerant) UMI distance against an independent string implementation of its
//! documented rules, and against the Hamming and Levenshtein distances it must never exceed.
use alevin_fry::umi_general_distance::umi_edit_distance_from_packed_shifted;
use alevin_fry::utils::count_diff_2_bit_packed;

fn pack(seq: &[u8]) -> u64 {
    seq.iter().fold(0u64, |v, c| (v << 2) | match c { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => panic!() })
}

/// The documented rules on strings: deleting from `x` is free before any of `y` is produced and
/// after all of it is; every other deletion, every insertion and every substitution costs 1.
fn directional(x: &[u8], y: &[u8]) -> u32 {
    let n = y.len();
    let mut prev: Vec<u32> = (0..=n as u32).collect();
    for &xi in x {
        let mut curr = vec![0u32; n + 1];
        for j in 1..=n {
            let sub = if xi == y[j - 1] { 0 } else { 1 };
            let del = if j == n { 0 } else { 1 };
            curr[j] = (prev[j - 1] + sub).min(prev[j] + del).min(curr[j - 1] + 1);
        }
        prev = curr;
    }
    prev[n]
}

fn general(x: &[u8], y: &[u8]) -> u32 {
    directional(x, y).min(directional(y, x))
}

fn levenshtein(x: &[u8], y: &[u8]) -> u32 {
    let n = y.len();
    let mut prev: Vec<u32> = (0..=n as u32).collect();
    for (i, &xi) in x.iter().enumerate() {
        let mut curr = vec![i as u32 + 1; n + 1];
        for j in 1..=n {
            let sub = if xi == y[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + sub).min(prev[j] + 1).min(curr[j - 1] + 1);
        }
        prev = curr;
    }
    prev[n]
}

fn hamming(x: &[u8], y: &[u8]) -> u32 {
    x.iter().zip(y).filter(|(a, b)| a != b).count() as u32
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13; self.0 ^= self.0 >> 7; self.0 ^= self.0 << 17; self.0
    }
    fn base(&mut self) -> u8 { b"ACGT"[(self.next() % 4) as usize] }
    fn umi(&mut self, n: usize) -> Vec<u8> { (0..n).map(|_| self.base()).collect() }
}

#[test]
fn packed_general_distance_matches_the_string_rules_and_never_exceeds_hamming_or_levenshtein() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let n = 12;
    let mut checked = 0;
    for k in 0..200_000 {
        let x = rng.umi(n);
        // random, or derived from x by shifts, substitutions and indels so that small distances are covered
        let y: Vec<u8> = match k % 6 {
            0 => rng.umi(n),
            1 => { let mut y = x[1..].to_vec(); y.push(rng.base()); y }                      // shift left
            2 => { let mut y = vec![rng.base()]; y.extend_from_slice(&x[..n - 1]); y }       // shift right
            3 => { let mut y = x.clone(); let p = (rng.next() % n as u64) as usize; y[p] = rng.base(); y }
            4 => { let mut y = x.clone(); let p = (rng.next() % n as u64) as usize; y.remove(p); y.push(rng.base()); y } // deletion inside
            _ => { let mut y = x.clone(); let p = (rng.next() % n as u64) as usize; y.insert(p, rng.base()); y.truncate(n); y } // insertion inside
        };
        let (px, py) = (pack(&x), pack(&y));
        let d = umi_edit_distance_from_packed_shifted(px, py, n as u8) as u32;
        assert_eq!(d, general(&x, &y), "{} {}", String::from_utf8_lossy(&x), String::from_utf8_lossy(&y));
        assert_eq!(d, umi_edit_distance_from_packed_shifted(py, px, n as u8) as u32, "symmetry");
        let h = count_diff_2_bit_packed(px, py) as u32;
        assert_eq!(h, hamming(&x, &y), "hamming on packed UMIs");
        assert!(d <= h, "general {} > hamming {} for {} {}", d, h, String::from_utf8_lossy(&x), String::from_utf8_lossy(&y));
        assert!(d <= levenshtein(&x, &y), "general above levenshtein");
        checked += 1;
    }
    assert_eq!(checked, 200_000);
}

#[test]
fn shifts_cost_one_per_base_and_identical_umis_cost_nothing() {
    let x = b"ACGTTGCAAGTC";
    assert_eq!(umi_edit_distance_from_packed_shifted(pack(x), pack(x), 12), 0);
    assert_eq!(umi_edit_distance_from_packed_shifted(pack(x), pack(b"CGTTGCAAGTCA"), 12), 1);
    assert_eq!(umi_edit_distance_from_packed_shifted(pack(x), pack(b"GTTGCAAGTCAA"), 12), 2);
    assert_eq!(umi_edit_distance_from_packed_shifted(pack(x), pack(b"TACGTTGCAAGT"), 12), 1);
    assert_eq!(count_diff_2_bit_packed(pack(x), pack(b"CGTTGCAAGTCA")), 10);
}
