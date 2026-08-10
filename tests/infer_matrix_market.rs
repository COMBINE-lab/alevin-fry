//! The equivalence-class count matrix `quant` writes must be readable by
//! `infer`.
//!
//! These two subcommands are meant to chain — `quant --dump-eq` produces
//! `geqc_counts.mtx` and `infer` consumes it — but they disagreed on the
//! MatrixMarket scalar type for long enough that `infer` simply could not read
//! its own pipeline's output (issue #150):
//!
//!   * `quant` builds the matrix as `sprs::TriMatI::<f32, u32>`, so the header
//!     it writes says `real`.
//!   * `infer` asked for `i32`, so the read failed with
//!     `MismatchedMatrixMarketRead(Integer, Float)`.
//!
//! sprs is strict in *both* directions, which is the part that makes this easy
//! to get wrong: widening an `integer` file to a float type fails too, so
//! "just read floats" would have broken matrices written by other tools. The
//! tests below pin both directions.

use std::path::Path;

/// Build a matrix exactly the way `quant` does (see `src/quant.rs`, which
/// accumulates UMI counts into an `f32` triplet matrix) and write it out.
fn write_quant_style_matrix(path: &Path) {
    let mut m = sprs::TriMatI::<f32, u32>::with_capacity((4, 6), 4);
    m.add_triplet(0, 0, 5f32);
    m.add_triplet(1, 2, 7f32);
    m.add_triplet(2, 5, 1f32);
    m.add_triplet(3, 1, 12f32);
    sprs::io::write_matrix_market(path, &m).expect("could not write matrix");
}

fn write_integer_matrix(path: &Path) {
    let mut m = sprs::TriMatI::<i32, u32>::with_capacity((4, 6), 3);
    m.add_triplet(0, 0, 5i32);
    m.add_triplet(1, 2, 7i32);
    m.add_triplet(3, 1, 12i32);
    sprs::io::write_matrix_market(path, &m).expect("could not write matrix");
}

/// Mirrors the read in `infer::infer`: try `real`, fall back to `integer`.
fn read_like_infer(path: &Path) -> Option<sprs::CsMatI<f32, u32>> {
    match sprs::io::read_matrix_market::<f32, u32, &Path>(path) {
        Ok(t) => Some(t.to_csr()),
        Err(_) => sprs::io::read_matrix_market::<i32, u32, &Path>(path)
            .ok()
            .map(|t| t.to_csr().map(|&v| v as f32)),
    }
}

/// The header `quant` writes today. If this ever changes, the fallback in
/// `infer` needs revisiting rather than silently taking the slow path.
#[test]
fn quant_writes_a_real_valued_matrix() {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("geqc_counts.mtx");
    write_quant_style_matrix(&p);

    let header = std::fs::read_to_string(&p).unwrap();
    let first = header.lines().next().unwrap();
    assert!(
        first.contains("real"),
        "expected quant's matrix to be real-valued, got header: {first}"
    );

    // The regression itself: reading that file as an integer type fails.
    assert!(
        sprs::io::read_matrix_market::<i32, u32, &Path>(p.as_path()).is_err(),
        "a real-valued matrix should not be readable as integer; if sprs has \
         relaxed this, the fallback in infer can be simplified"
    );
}

#[test]
fn infer_reads_the_matrix_quant_writes() {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("geqc_counts.mtx");
    write_quant_style_matrix(&p);

    let m = read_like_infer(&p).expect("infer could not read the matrix quant writes");
    assert_eq!(m.rows(), 4);
    assert_eq!(m.cols(), 6);
    assert_eq!(m.nnz(), 4);

    // Counts survive the round trip, and round rather than truncate.
    let mut got: Vec<u32> = m.iter().map(|(v, _)| v.round() as u32).collect();
    got.sort_unstable();
    assert_eq!(got, vec![1, 5, 7, 12]);
}

/// A matrix written by an older release or another tool may be integer-typed.
/// The fallback has to keep those readable.
#[test]
fn infer_still_reads_an_integer_matrix() {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("geqc_counts.mtx");
    write_integer_matrix(&p);

    // Guard the premise: this file genuinely is not readable as a float.
    assert!(
        sprs::io::read_matrix_market::<f32, u32, &Path>(p.as_path()).is_err(),
        "integer matrices used to be readable as floats; the fallback may be unnecessary"
    );

    let m = read_like_infer(&p).expect("infer could not read an integer-typed matrix");
    assert_eq!(m.nnz(), 3);
    let mut got: Vec<u32> = m.iter().map(|(v, _)| v.round() as u32).collect();
    got.sort_unstable();
    assert_eq!(got, vec![5, 7, 12]);
}
