//! `infer` over row names as `quant` writes them.
//!
//! A multi-sample `quant` names each row `<sample>_<barcode>`. `infer` used to
//! read every row name as a bare barcode k-mer, so it panicked on these. It also
//! stopped its workers after as many cells as the `--quant-subset` list had
//! entries, so a listed barcode that no row had left it waiting forever.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Output, Stdio};
use std::time::{Duration, Instant};

const BC: &str = "ACGTACGTACGTACGT";
const OTHER_BC: &str = "TTTTACGTACGTACGT";

/// The input directory of `infer`: two genes and three equivalence classes,
/// {G0}, {G1} and {G0, G1}, over the given rows, with these counts per row.
fn write_input(dir: &Path, rows: &[&str], counts: &[&[(usize, u32)]]) {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(dir.join("quants_mat_cols.txt"), "G0\nG1\n").unwrap();
    std::fs::write(dir.join("quants_mat_rows.txt"), rows.join("\n") + "\n").unwrap();
    let mut eqc = flate2::write::GzEncoder::new(
        std::fs::File::create(dir.join("gene_eqclass.txt.gz")).unwrap(),
        flate2::Compression::default(),
    );
    eqc.write_all(b"2\n3\n0\t0\n1\t1\n0\t1\t2\n").unwrap();
    eqc.finish().unwrap();
    let entries: Vec<String> = counts
        .iter()
        .enumerate()
        .flat_map(|(row, cell)| {
            cell.iter()
                .map(move |(eq, n)| format!("{} {} {}", row + 1, eq + 1, n))
        })
        .collect();
    std::fs::write(
        dir.join("geqc_counts.mtx"),
        format!(
            "%%MatrixMarket matrix coordinate real general\n{} 3 {}\n{}\n",
            rows.len(),
            entries.len(),
            entries.join("\n")
        ),
    )
    .unwrap();
}

/// Run `infer`, killing it if it has not finished within a minute.
fn infer(input: &Path, output: &Path, subset: Option<&Path>) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_alevin-fry"));
    cmd.arg("infer")
        .arg("-c")
        .arg(input.join("geqc_counts.mtx"))
        .arg("-e")
        .arg(input.join("gene_eqclass.txt.gz"))
        .arg("-o")
        .arg(output)
        .args(["-t", "2"])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    if let Some(subset) = subset {
        cmd.arg("--quant-subset").arg(subset);
    }
    let mut child = cmd.spawn().unwrap();
    let start = Instant::now();
    while child.try_wait().unwrap().is_none() {
        if start.elapsed() > Duration::from_secs(60) {
            child.kill().unwrap();
            panic!("infer did not finish within 60 s");
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    child.wait_with_output().unwrap()
}

/// The output as row name -> gene counts.
fn read_output(dir: &Path) -> HashMap<String, [f64; 2]> {
    let rows: Vec<String> = std::fs::read_to_string(dir.join("quants_mat_rows.txt"))
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect();
    let mut out: HashMap<String, [f64; 2]> = rows.iter().map(|r| (r.clone(), [0.0; 2])).collect();
    let text = std::fs::read_to_string(dir.join("quants_mat.mtx")).unwrap();
    for line in text.lines().filter(|l| !l.starts_with('%')).skip(1) {
        let v: Vec<&str> = line.split_whitespace().collect();
        let row = &rows[v[0].parse::<usize>().unwrap() - 1];
        out.get_mut(row).unwrap()[v[1].parse::<usize>().unwrap() - 1] = v[2].parse().unwrap();
    }
    out
}

fn multi_sample_input(dir: &Path) -> [String; 3] {
    let rows = [
        format!("A-A01_{BC}"),
        format!("A-A02_{BC}"),
        format!("A-A02_{OTHER_BC}"),
    ];
    write_input(
        dir,
        &rows.iter().map(String::as_str).collect::<Vec<_>>(),
        &[&[(0, 3)], &[(1, 2)], &[(0, 1), (1, 1)]],
    );
    rows
}

#[test]
fn infer_keeps_multi_sample_row_names() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("in");
    let rows = multi_sample_input(&input);

    let output = infer(&input, &tmp.path().join("out"), None);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let counts = read_output(&tmp.path().join("out"));
    assert_eq!(counts.len(), 3);
    assert_eq!(counts[&rows[0]], [3.0, 0.0]);
    assert_eq!(counts[&rows[1]], [0.0, 2.0]);
    assert_eq!(counts[&rows[2]], [1.0, 1.0]);
}

#[test]
fn infer_subset_keeps_a_barcode_in_every_sample_and_skips_absent_ones() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("in");
    let rows = multi_sample_input(&input);
    // one barcode in two samples, and one that no row has
    let subset = tmp.path().join("subset.txt");
    std::fs::write(&subset, format!("{BC}\nGGGGGGGGGGGGGGGG\n")).unwrap();

    let output = infer(&input, &tmp.path().join("out"), Some(&subset));
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let counts = read_output(&tmp.path().join("out"));
    assert_eq!(counts.len(), 2);
    assert_eq!(counts[&rows[0]], [3.0, 0.0]);
    assert_eq!(counts[&rows[1]], [0.0, 2.0]);
}

#[test]
fn infer_subset_keeps_a_named_row_only() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("in");
    let rows = multi_sample_input(&input);
    let subset = tmp.path().join("subset.txt");
    std::fs::write(&subset, format!("{}\n", rows[1])).unwrap();

    let output = infer(&input, &tmp.path().join("out"), Some(&subset));
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let counts = read_output(&tmp.path().join("out"));
    assert_eq!(counts.len(), 1);
    assert_eq!(counts[&rows[1]], [0.0, 2.0]);
}

#[test]
fn infer_subset_on_single_sample_rows() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("in");
    write_input(&input, &[BC, OTHER_BC], &[&[(0, 3)], &[(2, 4)]]);
    let subset = tmp.path().join("subset.txt");
    // with a barcode that no row has, which used to leave the workers waiting
    // (Cell Ranger's `-1` suffix is ignored, as before)
    std::fs::write(&subset, format!("{OTHER_BC}-1\nGGGGGGGGGGGGGGGG\n")).unwrap();

    let output = infer(&input, &tmp.path().join("out"), Some(&subset));
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let counts = read_output(&tmp.path().join("out"));
    assert_eq!(counts.len(), 1);
    // an ambiguous class with no unique evidence for either gene is split
    let [g0, g1] = counts[OTHER_BC];
    assert!((g0 + g1 - 4.0).abs() < 1e-3, "{g0} + {g1}");
}
