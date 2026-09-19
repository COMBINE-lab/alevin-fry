//! Exercise the actual quantification writer, including its patched size line.

use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File};
use std::path::{Path, PathBuf};

use alevin_fry::prog_opts::QuantOpts;
use alevin_fry::quant::{ResolutionStrategy, SplicedAmbiguityModel, quantify};
use libradicl::chunk::Chunk;
use libradicl::collation::{CollationManifest, SampleGroup};
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{
    RadIntId, RadType, TagDesc, TagMap, TagSection, TagSectionLabel, TagValue,
};
use libradicl::record::{MultiBarcodeReadRecord, MultiBarcodeRecordContext, RecordContext};
use libradicl::writers::RadFileWriter;
use needletail::bitkmer::bitmer_to_bytes;
use smallvec::smallvec;

const BC_LEN: u16 = 16;
const NUM_GENES: usize = 3;
type Entries = BTreeMap<(String, String), f32>;

struct Fixture {
    input: PathBuf,
    tg_map: PathBuf,
    rows: HashSet<String>,
    entries: Entries,
}

fn barcode(cell: usize) -> u64 {
    (cell as u64 + 1).wrapping_mul(2_654_435_761) & u32::MAX as u64
}

fn barcode_label(cell: usize) -> String {
    String::from_utf8(bitmer_to_bytes((barcode(cell), BC_LEN as u8))).unwrap()
}

/// Each chunk already has the sample index in b0 and one cell barcode in b1,
/// just as collate writes it. Samples deliberately reuse their cell barcodes.
fn fixture(root: &Path, cells_per_sample: usize, ambiguous_only: bool) -> Fixture {
    let input = root.join("input");
    fs::create_dir_all(&input).unwrap();
    fs::write(input.join("collate.json"), r#"{"compressed_output":false}"#).unwrap();
    let tg_map = root.join("tg_map.tsv");
    fs::write(&tg_map, "tx0\tgene0\ntx1\tgene1\ntx2\tgene2\n").unwrap();

    let mut file_tags = TagSection::new_with_label(TagSectionLabel::FileTags);
    for name in ["num_barcodes", "b0len", "b1len", "ulen"] {
        file_tags.add_tag_desc(TagDesc::new(name, RadType::Int(RadIntId::U16)));
    }
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    for name in ["b0", "b1", "u"] {
        read_tags.add_tag_desc(TagDesc::new(name, RadType::Int(RadIntId::U32)));
    }
    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    aln_tags.add_tag_desc(TagDesc::new(
        "compressed_ori_refid",
        RadType::Int(RadIntId::U32),
    ));
    let prelude = RadPrelude {
        hdr: RadHeader {
            version: libradicl::header::SpecVersion::Legacy,
            is_paired: 0,
            ref_count: NUM_GENES as u64,
            ref_names: (0..NUM_GENES).map(|i| format!("tx{i}")).collect(),
            num_chunks: 0,
        },
        file_tags,
        read_tags,
        aln_tags,
    };
    let mut tags = TagMap::with_keyset(&prelude.file_tags.tags);
    for value in [2, 8, BC_LEN, 12] {
        tags.add(TagValue::U16(value));
    }
    let context = MultiBarcodeRecordContext::get_context_from_tag_section(
        &prelude.file_tags,
        &prelude.read_tags,
        &prelude.aln_tags,
    )
    .unwrap();
    let mut writer = RadFileWriter::new(
        File::create(input.join("map.collated.rad")).unwrap(),
        &prelude,
        &tags,
    )
    .unwrap();
    let mut manifest = CollationManifest::new(vec!["sample".into(), "cell".into()]);
    let mut rows = HashSet::new();
    let mut entries = Entries::new();
    for sample in 0..2 {
        let mut num_records = 0;
        for cell in 0..cells_per_sample {
            let label = format!("sample{sample}_{}", barcode_label(cell));
            rows.insert(label.clone());
            let mut reads = Vec::new();
            if ambiguous_only || cell % 17 == 0 {
                // CR-like drops these gene-ambiguous molecules, retaining an
                // empty matrix row. Uniform EM instead assigns half to each.
                reads.push(MultiBarcodeReadRecord {
                    barcodes: smallvec![sample as u64, barcode(cell)],
                    umi: 0,
                    dirs: vec![true, true],
                    refs: vec![0, 1],
                });
            } else {
                for (gene, count) in [(0, 1 + cell % 31), (1, 1 + (cell + 7 * sample) % 29)] {
                    entries.insert((label.clone(), format!("gene{gene}")), count as f32);
                    for umi in 0..count {
                        reads.push(MultiBarcodeReadRecord {
                            barcodes: smallvec![sample as u64, barcode(cell)],
                            umi: (gene * 100 + umi) as u64,
                            dirs: vec![true],
                            refs: vec![gene as u32],
                        });
                    }
                }
            }
            num_records += reads.len() as u64;
            writer
                .write_chunk(
                    &Chunk {
                        nbytes: 0,
                        nrec: reads.len() as u32,
                        reads,
                    },
                    &context,
                )
                .unwrap();
        }
        manifest.add_sample_group(SampleGroup {
            key: sample as u64,
            name: Some(format!("sample{sample}")),
            chunk_start: (sample * cells_per_sample) as u64,
            num_chunks: cells_per_sample as u64,
            num_records,
        });
    }
    writer.finalize().unwrap();
    manifest
        .write_to_file(&input.join("collation_manifest.bin"))
        .unwrap();
    Fixture {
        input,
        tg_map,
        rows,
        entries,
    }
}

fn run_quant(
    fixture: &Fixture,
    output: &PathBuf,
    threads: u32,
    fractional: bool,
    filter: Option<&PathBuf>,
) -> anyhow::Result<()> {
    run_quant_with_bootstraps(fixture, output, threads, fractional, filter, 0, false)
}

fn run_quant_with_bootstraps(
    fixture: &Fixture,
    output: &PathBuf,
    threads: u32,
    fractional: bool,
    filter: Option<&PathBuf>,
    num_bootstraps: u32,
    summary_stat: bool,
) -> anyhow::Result<()> {
    let log = slog::Logger::root(slog::Discard, slog::o!());
    quantify(
        QuantOpts::builder()
            .input_dir(&fixture.input)
            .tg_map(&fixture.tg_map)
            .output_dir(output)
            .num_threads(threads)
            .num_bootstraps(num_bootstraps)
            .init_uniform(true)
            .summary_stat(summary_stat)
            .dump_eq(false)
            .resolution(if fractional {
                ResolutionStrategy::CellRangerLikeEm
            } else {
                ResolutionStrategy::CellRangerLike
            })
            .pug_exact_umi(false)
            .sa_model(SplicedAmbiguityModel::WinnerTakeAll)
            .small_thresh(if fractional { 0 } else { 100 })
            .large_graph_thresh(0)
            .filter_list(filter)
            .cmdline("streaming-matrix-market-test")
            .version(env!("CARGO_PKG_VERSION"))
            .log(&log)
            .build(),
    )
}

fn read_output(output: &Path) -> (HashSet<String>, Entries) {
    read_matrix_output(output, "quants_mat.mtx", &["gene0", "gene1", "gene2"])
}

fn read_matrix_output(
    output: &Path,
    filename: &str,
    expected_cols: &[&str],
) -> (HashSet<String>, Entries) {
    let dir = output.join("alevin");
    let rows_text = fs::read_to_string(dir.join("quants_mat_rows.txt")).unwrap();
    let cols_text = fs::read_to_string(dir.join("quants_mat_cols.txt")).unwrap();
    let rows: Vec<_> = rows_text.lines().collect();
    let cols: Vec<_> = cols_text.lines().collect();
    let matrix = sprs::io::read_matrix_market::<f32, u32, _>(dir.join(filename))
        .expect("quant output must be accepted by a MatrixMarket reader");
    assert_eq!(matrix.shape(), (rows.len(), cols.len()));
    assert_eq!(cols, expected_cols);
    let unique_rows: HashSet<_> = rows.iter().map(|s| s.to_string()).collect();
    assert_eq!(unique_rows.len(), rows.len());

    let mut entries = Entries::new();
    for (&value, (row, col)) in matrix.triplet_iter() {
        assert!(value > 0.0);
        let key = (
            rows[row as usize].to_string(),
            cols[col as usize].to_string(),
        );
        assert!(entries.insert(key, value).is_none(), "duplicate coordinate");
    }
    assert_eq!(matrix.nnz(), entries.len());
    // sprs permits extra body lines; independently check that nnz is exact.
    let text = fs::read_to_string(dir.join(filename)).unwrap();
    let body: Vec<_> = text.lines().filter(|line| !line.starts_with('%')).collect();
    assert_eq!(body.len(), matrix.nnz() + 1);
    assert_eq!(body[0].split_whitespace().count(), 3);
    (unique_rows, entries)
}

fn assert_no_bootstrap_files(output: &Path) {
    for filename in ["bootstraps_mean.mtx", "bootstraps_var.mtx"] {
        assert!(
            !output.join("alevin").join(filename).exists(),
            "unexpected {filename}"
        );
    }
}

#[test]
fn streamed_matrix_matches_every_barcode_and_gene_with_multiple_workers() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = fixture(tmp.path(), 2048, false);
    // More than several 524208-byte RAD meta chunks: the multiworker case can
    // actually distribute work instead of handing a tiny fixture to one worker.
    assert!(
        fs::metadata(fixture.input.join("map.collated.rad"))
            .unwrap()
            .len()
            > 2_000_000
    );
    for threads in [1, 5] {
        let output = tmp.path().join(format!("quant-{threads}"));
        run_quant(&fixture, &output, threads, false, None).unwrap();
        let (rows, entries) = read_output(&output);
        assert_eq!(rows, fixture.rows);
        assert_eq!(entries, fixture.entries);
        assert_no_bootstrap_files(&output);
    }
}

#[test]
fn streamed_matrix_round_trips_zero_nonzeros_and_fractional_values() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = fixture(tmp.path(), 4, true);
    for fractional in [false, true] {
        let output = tmp.path().join(format!("quant-{fractional}"));
        run_quant(&fixture, &output, 5, fractional, None).unwrap();
        let (rows, entries) = read_output(&output);
        assert_eq!(rows, fixture.rows);
        let expected: Entries = if fractional {
            fixture
                .rows
                .iter()
                .flat_map(|row| {
                    [(row.clone(), "gene0".into()), (row.clone(), "gene1".into())]
                        .into_iter()
                        .map(|key| (key, 0.5))
                })
                .collect()
        } else {
            Entries::new()
        };
        assert_eq!(entries, expected);
        assert_no_bootstrap_files(&output);
    }
}

#[test]
fn shared_barcode_filter_cannot_succeed_with_out_of_bounds_matrix_rows() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = fixture(tmp.path(), 2, false);
    let filter = tmp.path().join("filter.txt");
    fs::write(&filter, format!("{}\n", barcode_label(1))).unwrap();
    let output = tmp.path().join("filtered");
    // The filter has one barcode, but it selects a distinct cell in each of
    // two samples. Returning an error is acceptable until filters support this
    // case; reporting success with a one-row matrix containing row 2 is not.
    if run_quant(&fixture, &output, 1, false, Some(&filter)).is_ok() {
        let (rows, entries) = read_output(&output);
        let expected_rows = fixture
            .rows
            .iter()
            .filter(|row| row.ends_with(&barcode_label(1)))
            .cloned()
            .collect();
        assert_eq!(rows, expected_rows);
        assert_eq!(entries, fixture.entries);
    }
}

#[test]
fn streamed_bootstrap_matrices_preserve_rows_and_statistics_in_both_modes() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = fixture(tmp.path(), 2048, false);
    // With all references assigned to one gene, every cell has one bootstrap
    // category. Resampling cannot change its molecule count, so the expected
    // mean and zero variance are deterministic without controlling an RNG.
    fs::write(&fixture.tg_map, "tx0\tgene0\ntx1\tgene0\ntx2\tgene0\n").unwrap();
    let expected: Entries = fixture
        .rows
        .iter()
        .map(|row| {
            let total: f32 = fixture
                .entries
                .range((row.clone(), String::new())..)
                .take_while(|((barcode, _), _)| barcode == row)
                .map(|(_, value)| *value)
                .sum();
            // A previously gene-ambiguous row now has one unique molecule.
            ((row.clone(), "gene0".into()), total.max(1.0))
        })
        .collect();
    for summary_stat in [false, true] {
        let output = tmp.path().join(format!("bootstrap-{summary_stat}"));
        run_quant_with_bootstraps(&fixture, &output, 5, true, None, 3, summary_stat).unwrap();
        let (rows, counts) = read_matrix_output(&output, "quants_mat.mtx", &["gene0"]);
        assert_eq!(rows, fixture.rows);
        assert_eq!(counts, expected);
        let (mean_rows, means) = read_matrix_output(&output, "bootstraps_mean.mtx", &["gene0"]);
        assert_eq!(mean_rows, fixture.rows);
        assert_eq!(means, expected);
        let (variance_rows, variances) =
            read_matrix_output(&output, "bootstraps_var.mtx", &["gene0"]);
        assert_eq!(variance_rows, fixture.rows);
        assert!(variances.is_empty());
    }
}

#[test]
fn no_bootstrap_files_are_kept_when_no_cells_produce_bootstrap_means() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = fixture(tmp.path(), 0, false);
    for summary_stat in [false, true] {
        let output = tmp.path().join(format!("empty-{summary_stat}"));
        run_quant_with_bootstraps(&fixture, &output, 5, true, None, 3, summary_stat).unwrap();
        let (rows, entries) = read_output(&output);
        assert!(rows.is_empty());
        assert!(entries.is_empty());
        assert_no_bootstrap_files(&output);
    }
}

#[test]
fn streamed_usa_matrix_preserves_spliced_unspliced_and_ambiguous_columns() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = fixture(tmp.path(), 2, false);
    fs::write(
        &fixture.tg_map,
        "tx0\tgene0\tS\ntx1\tgene0\tU\ntx2\tgene1\tS\n",
    )
    .unwrap();
    let output = tmp.path().join("usa");
    run_quant(&fixture, &output, 5, false, None).unwrap();
    let cols = ["gene0", "gene1", "gene0-U", "gene1-U", "gene0-A", "gene1-A"];
    let (rows, entries) = read_matrix_output(&output, "quants_mat.mtx", &cols);
    assert_eq!(rows, fixture.rows);
    let mut expected = Entries::new();
    for sample in 0..2 {
        expected.insert(
            (
                format!("sample{sample}_{}", barcode_label(0)),
                "gene0-A".into(),
            ),
            1.0,
        );
        let row = format!("sample{sample}_{}", barcode_label(1));
        expected.insert((row.clone(), "gene0".into()), 2.0);
        expected.insert((row, "gene0-U".into()), (2 + 7 * sample) as f32);
    }
    assert_eq!(entries, expected);
    assert_no_bootstrap_files(&output);
}

// Run failures in a subprocess: a queue-draining regression must fail with a
// deadline rather than leave the test runner blocked behind a full queue.
fn wait_for_child(child: &mut std::process::Child, description: &str) -> std::process::ExitStatus {
    use std::time::{Duration, Instant};

    let deadline = Instant::now() + Duration::from_secs(30);
    loop {
        if let Some(status) = child.try_wait().unwrap() {
            return status;
        }
        if Instant::now() >= deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("{description} blocked quantification for 30 seconds");
        }
        std::thread::sleep(Duration::from_millis(25));
    }
}

#[cfg(target_os = "linux")]
#[test]
fn output_write_errors_return_without_panicking_or_blocking_the_producer() {
    use std::os::unix::fs::symlink;
    use std::process::{Command, Stdio};

    let tmp = tempfile::tempdir().unwrap();
    let fixture = fixture(tmp.path(), 20_000, false);
    fs::write(&fixture.tg_map, "tx0\tgene0\ntx1\tgene0\ntx2\tgene0\n").unwrap();
    // A single worker has a queue of four 524208-byte meta chunks. This leaves
    // many chunks unconsumed when any ordinary output buffer first fills.
    assert!(
        fs::metadata(fixture.input.join("map.collated.rad"))
            .unwrap()
            .len()
            > 20_000_000
    );
    for (name, failing_path, bootstraps) in [
        ("matrix", "alevin/quants_mat.mtx", false),
        ("barcodes", "alevin/quants_mat_rows.txt", false),
        ("features", "featureDump.txt", false),
        ("bootstrap-mean", "alevin/bootstraps_mean.mtx", true),
        ("bootstrap-variance", "alevin/bootstraps_var.mtx", true),
    ] {
        let output = tmp.path().join(name);
        fs::create_dir_all(output.join("alevin")).unwrap();
        symlink("/dev/full", output.join(failing_path)).unwrap();
        let log_path = output.join("child.log");
        let log = File::create(&log_path).unwrap();
        let mut child = Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "output_write_failure_child",
                "--ignored",
                "--nocapture",
            ])
            .env("AF_WRITE_FAILURE_FIXTURE", tmp.path())
            .env("AF_WRITE_FAILURE_OUTPUT", &output)
            .env(
                "AF_WRITE_FAILURE_BOOTSTRAPS",
                if bootstraps { "3" } else { "0" },
            )
            .stdout(Stdio::from(log.try_clone().unwrap()))
            .stderr(Stdio::from(log))
            .spawn()
            .unwrap();
        let status = wait_for_child(&mut child, &format!("{name} output failure"));
        assert!(
            status.success(),
            "{name} output failure did not return Err normally:\n{}",
            fs::read_to_string(&log_path).unwrap()
        );
        assert!(!output.join("quant.json").exists());
    }
}

#[cfg(target_os = "linux")]
#[test]
#[ignore = "subprocess helper for output_write_errors_return_without_panicking_or_blocking_the_producer"]
fn output_write_failure_child() {
    let root = PathBuf::from(std::env::var_os("AF_WRITE_FAILURE_FIXTURE").unwrap());
    let output = PathBuf::from(std::env::var_os("AF_WRITE_FAILURE_OUTPUT").unwrap());
    let num_bootstraps = std::env::var("AF_WRITE_FAILURE_BOOTSTRAPS")
        .unwrap()
        .parse()
        .unwrap();
    let fixture = Fixture {
        input: root.join("input"),
        tg_map: root.join("tg_map.tsv"),
        rows: HashSet::new(),
        entries: Entries::new(),
    };
    // One worker forces an erroring consumer to drain its queue correctly;
    // an unwinding panic fails this child test, and a deadlock hits the parent's
    // deadline. No catch_unwind hides either failure mode.
    let error = run_quant_with_bootstraps(
        &fixture,
        &output,
        1,
        num_bootstraps > 0,
        None,
        num_bootstraps,
        true,
    )
    .expect_err("writing to /dev/full must fail");
    eprintln!("quantification returned its output error: {error:#}");
}

#[cfg(target_os = "linux")]
#[test]
fn quant_cli_reports_output_error_with_context_and_exit_code_one() {
    use std::os::unix::fs::symlink;
    use std::process::{Command, Stdio};

    let tmp = tempfile::tempdir().unwrap();
    let fixture = fixture(tmp.path(), 2, false);
    fs::write(
        fixture.input.join("generate_permit_list.json"),
        r#"{"velo_mode":false}"#,
    )
    .unwrap();
    let output = tmp.path().join("cli-error");
    fs::create_dir_all(output.join("alevin")).unwrap();
    symlink("/dev/full", output.join("alevin/quants_mat.mtx")).unwrap();
    let log_path = output.join("cli.log");
    let log = File::create(&log_path).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_alevin-fry"))
        .arg("quant")
        .arg("-i")
        .arg(&fixture.input)
        .arg("-m")
        .arg(&fixture.tg_map)
        .arg("-o")
        .arg(&output)
        .args(["-t", "2", "-r", "cr-like"])
        .stdout(Stdio::from(log.try_clone().unwrap()))
        .stderr(Stdio::from(log))
        .spawn()
        .unwrap();
    let status = wait_for_child(&mut child, "CLI matrix output error");
    let diagnostic = fs::read_to_string(log_path).unwrap();
    assert_eq!(status.code(), Some(1), "{diagnostic}");
    assert!(
        diagnostic.contains("could not finalize count matrix"),
        "{diagnostic}"
    );
    assert!(diagnostic.contains("os error 28"), "{diagnostic}");
    assert!(!diagnostic.contains("panicked at"), "{diagnostic}");
    assert!(!output.join("quant.json").exists());
}

#[test]
fn truncated_collated_input_returns_error_and_joins_waiting_workers() {
    use std::process::{Command, Stdio};

    let tmp = tempfile::tempdir().unwrap();
    let fixture = fixture(tmp.path(), 32, false);
    let path = fixture.input.join("map.collated.rad");
    let length = fs::metadata(&path).unwrap().len();
    // Preserve the valid prelude and chunk headers, but interrupt the final
    // record body. Workers waiting on the producer must still be released.
    fs::OpenOptions::new()
        .write(true)
        .open(path)
        .unwrap()
        .set_len(length - 3)
        .unwrap();
    let output = tmp.path().join("truncated");
    fs::create_dir_all(&output).unwrap();
    let log_path = output.join("child.log");
    let log = File::create(&log_path).unwrap();
    let mut child = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "truncated_input_failure_child",
            "--ignored",
            "--nocapture",
        ])
        .env("AF_TRUNCATED_INPUT_FIXTURE", tmp.path())
        .env("AF_TRUNCATED_INPUT_OUTPUT", &output)
        .stdout(Stdio::from(log.try_clone().unwrap()))
        .stderr(Stdio::from(log))
        .spawn()
        .unwrap();
    let status = wait_for_child(&mut child, "truncated input");
    assert!(
        status.success(),
        "{}",
        fs::read_to_string(log_path).unwrap()
    );
    assert!(!output.join("quant.json").exists());
}

#[test]
#[ignore = "subprocess helper for truncated_collated_input_returns_error_and_joins_waiting_workers"]
fn truncated_input_failure_child() {
    let root = PathBuf::from(std::env::var_os("AF_TRUNCATED_INPUT_FIXTURE").unwrap());
    let output = PathBuf::from(std::env::var_os("AF_TRUNCATED_INPUT_OUTPUT").unwrap());
    let fixture = Fixture {
        input: root.join("input"),
        tg_map: root.join("tg_map.tsv"),
        rows: HashSet::new(),
        entries: Entries::new(),
    };
    let error = run_quant(&fixture, &output, 5, false, None)
        .expect_err("truncated RAD body must not produce successful output");
    assert!(
        format!("{error:#}").contains("could not read collated RAD input"),
        "{error:#}"
    );
}
