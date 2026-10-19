/*
 * Copyright (c) 2020-2026 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! The per-molecule table written by `quant --dump-molecules`.
//!
//! Every molecule a resolution strategy produces for a cell becomes one row of
//! `alevin/molecules.parquet`, before molecules are summed into the count
//! matrix. A row is keyed by a *representative* UMI, which is not always a
//! single observed UMI:
//!
//! * `cr-like` / `cr-like-em`: the UMI itself, or with `--umi-edit-dist 1` the
//!   corrected UMI that Hamming-1 neighbours were merged into.
//! * `parsimony*`: the root of the arborescence that covers this molecule in
//!   the parsimonious UMI graph; the covered UMIs are collapsed into it.
//! * `trivial`: the UMI itself.
//!
//! `n_umis` is the number of distinct observed UMIs collapsed into the row and
//! `reads` the number of reads supporting it (see [`schema`]).
//!
//! The file is LZ4-compressed (Parquet `LZ4_RAW`), which pyarrow >= 10,
//! polars, DuckDB and recent R `arrow` read.
//!
//! Writing is parallel: each quant worker buffers rows for its own cells,
//! encodes and compresses a whole row group without holding any lock, and only
//! the splice of the finished column chunks into the shared file is serialized.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::sync::{Arc, Mutex};

use anyhow::Context;
use arrow_array::builder::{ArrayBuilder, Int32Builder, ListBuilder, StringBuilder, UInt32Builder};
use arrow_array::{ArrayRef, Int32DictionaryArray, RecordBatch, StringArray};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use parquet::arrow::ArrowWriter;
use parquet::arrow::arrow_writer::{ArrowRowGroupWriterFactory, compute_leaves};
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use parquet::file::writer::SerializedFileWriter;
use parquet::schema::types::ColumnPath;

use crate::utils as afutils;

/// File name of the molecule table inside the `alevin` output directory.
pub const MOLECULE_TABLE_NAME: &str = "molecules.parquet";

/// Rows a worker buffers before encoding them as one Parquet row group.
///
/// Every worker holds one partly built row group, so this sets the table's
/// memory per quant thread. On a 33M-molecule 10x v3 run, 128Ki rows cost
/// about 20 MiB of peak RSS per thread and a 4% larger file than 256Ki rows,
/// which cost about 30 MiB per thread; 32Ki-64Ki rows saved little memory
/// but grew the file by 5-25%.
const ROW_GROUP_ROWS: usize = 128 * 1024;

/// Row status values, as written to the `status` column.
pub const STATUS_ASSIGNED: &str = "assigned";
pub const STATUS_MULTI_GENE: &str = "multi_gene";
pub const STATUS_LOW_SUPPORT: &str = "low_support";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct MolRow {
    umi: u64,
    reads: u32,
    n_umis: u32,
    label_start: u32,
    label_len: u32,
    low_support: bool,
}

/// The molecules a resolver produced for one cell, reused across cells.
///
/// Labels are the resolver's own gene ids: plain gene ids, or in USA mode the
/// interleaved spliced (even) / unspliced (odd) ids. They are mapped to output
/// features only when the cell is written, by [`FeatureLabeler`].
#[derive(Default, Debug)]
pub struct CellMolecules {
    rows: Vec<MolRow>,
    labels: Vec<u32>,
}

impl CellMolecules {
    pub fn clear(&mut self) {
        self.rows.clear();
        self.labels.clear();
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Record a resolved molecule whose compatible gene labels are `labels`.
    pub fn push(&mut self, umi: u64, labels: &[u32], reads: u32, n_umis: u32) {
        self.push_row(umi, labels, reads, n_umis, false);
    }

    /// Record a UMI that the cr-like low-support filter removed entirely, with
    /// its best-supported gene label(s).
    pub fn push_low_support(&mut self, umi: u64, labels: &[u32], reads: u32, n_umis: u32) {
        self.push_row(umi, labels, reads, n_umis, true);
    }

    fn push_row(&mut self, umi: u64, labels: &[u32], reads: u32, n_umis: u32, low_support: bool) {
        self.rows.push(MolRow {
            umi,
            reads,
            n_umis,
            label_start: self.labels.len() as u32,
            label_len: labels.len() as u32,
            low_support,
        });
        self.labels.extend_from_slice(labels);
    }

    /// `(umi, labels, reads, n_umis, low_support)` for each row, in push order.
    pub fn iter(&self) -> impl Iterator<Item = (u64, &[u32], u32, u32, bool)> + '_ {
        self.rows.iter().map(|r| {
            let s = r.label_start as usize;
            (
                r.umi,
                &self.labels[s..s + r.label_len as usize],
                r.reads,
                r.n_umis,
                r.low_support,
            )
        })
    }
}

/// How a molecule's labels map onto the columns of the count matrix.
#[derive(Debug, PartialEq, Eq)]
pub enum Assignment {
    /// Counted as one whole molecule of this matrix column.
    Feature(usize),
    /// Gene-ambiguous: discarded by the non-EM strategies, apportioned by the
    /// EM strategies.
    MultiGene,
}

/// Maps resolver labels to matrix columns and feature names, exactly as the
/// count matrix does, so that the table's `assigned` rows reproduce it.
pub struct FeatureLabeler {
    /// Column names of the count matrix, i.e. `quants_mat_cols.txt`.
    names: Arc<Vec<String>>,
    usa_mode: bool,
    unspliced_offset: usize,
    ambig_offset: usize,
}

impl FeatureLabeler {
    /// `names` are the matrix column names; in USA mode `num_rows / 3` of each
    /// of spliced, unspliced and ambiguous.
    pub fn new(names: Arc<Vec<String>>, usa_mode: bool) -> Self {
        let unspliced_offset = if usa_mode { names.len() / 3 } else { 0 };
        Self {
            names,
            usa_mode,
            unspliced_offset,
            ambig_offset: 2 * unspliced_offset,
        }
    }

    /// The matrix column of a single label.
    fn label_column(&self, label: u32) -> usize {
        if !self.usa_mode {
            label as usize
        } else if afutils::is_spliced(label) {
            (label >> 1) as usize
        } else {
            self.unspliced_offset + (label >> 1) as usize
        }
    }

    /// The matrix name of a single label.
    pub fn label_name(&self, label: u32) -> &str {
        &self.names[self.label_column(label)]
    }

    /// How the count matrix treats a molecule with these (sorted) labels.
    ///
    /// `em` is whether the cell's counts came from an EM strategy. There a
    /// USA-mode label set only counts as one feature when it is splicing, not
    /// gene, ambiguous; the non-EM strategies additionally apply the
    /// prefer-spliced rule of [`afutils::usa_label_column`].
    pub fn assign(&self, labels: &[u32], em: bool) -> Assignment {
        if !self.usa_mode {
            return match labels {
                [g] => Assignment::Feature(*g as usize),
                _ => Assignment::MultiGene,
            };
        }
        let single_gene = match labels {
            [_] => true,
            [g1, g2] => afutils::same_gene(*g1, *g2, true),
            _ => false,
        };
        if em && !single_gene {
            return Assignment::MultiGene;
        }
        match afutils::usa_label_column(labels, self.unspliced_offset, self.ambig_offset) {
            Some(col) => Assignment::Feature(col),
            None => Assignment::MultiGene,
        }
    }
}

/// Arrow schema of the molecule table. `sample` is present only for
/// multi-sample (e.g. 10x Flex) input.
///
/// * `cell_barcode` — corrected cell barcode, as in `quants_mat_rows.txt`
///   (without the sample prefix; that is in `sample`).
/// * `rep_umi` — the representative UMI of the molecule.
/// * `n_umis` — distinct observed UMIs collapsed into this molecule.
/// * `feature` — the count-matrix column this molecule was counted in (for a
///   `low_support` row, the single gene it was dropped from); null when the
///   molecule is compatible with several genes.
/// * `candidate_features` — when `feature` is null: every feature the
///   molecule is compatible with (for `low_support`, the tied best-supported
///   genes); null otherwise.
/// * `reads` — reads supporting the molecule's assignment. For cr-like, the
///   reads of the (corrected) UMI compatible with the winning gene(s); for
///   parsimony, the reads of every UMI in the covering arborescence.
/// * `status` — `assigned`, `multi_gene` or `low_support`.
///
/// The low-cardinality string columns are Arrow dictionary columns, so they
/// load as categoricals and a row carries a 4-byte key rather than a string.
pub fn schema(multi_sample: bool) -> SchemaRef {
    let categorical = || DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8));
    let mut fields = vec![Field::new("cell_barcode", categorical(), false)];
    if multi_sample {
        // null when the collation manifest has no name for the cell's sample
        fields.push(Field::new("sample", categorical(), true));
    }
    fields.extend([
        Field::new("rep_umi", DataType::Utf8, false),
        Field::new("n_umis", DataType::UInt32, false),
        Field::new("feature", categorical(), true),
        Field::new(
            "candidate_features",
            DataType::List(Arc::new(Field::new_list_field(DataType::Utf8, true))),
            true,
        ),
        Field::new("reads", DataType::UInt32, false),
        Field::new("status", categorical(), false),
    ]);
    Arc::new(Schema::new(fields))
}

/// Dictionary keys of the `status` column.
const STATUS_VALUES: [&str; 3] = [STATUS_ASSIGNED, STATUS_MULTI_GENE, STATUS_LOW_SUPPORT];
const STATUS_KEY_ASSIGNED: i32 = 0;
const STATUS_KEY_MULTI_GENE: i32 = 1;
const STATUS_KEY_LOW_SUPPORT: i32 = 2;

/// The shared file. Workers hold only an `Arc` to it.
pub struct MoleculeTableWriter {
    file: Mutex<Option<SerializedFileWriter<BufWriter<File>>>>,
    factory: ArrowRowGroupWriterFactory,
    schema: SchemaRef,
    labeler: FeatureLabeler,
    /// Dictionary values shared (not copied) by every batch: the matrix
    /// column names, the sample names, and the status values.
    feature_values: ArrayRef,
    sample_values: Option<ArrayRef>,
    status_values: ArrayRef,
    umi_len: usize,
}

impl MoleculeTableWriter {
    /// `sample_names`, indexed by sample index, for multi-sample input.
    pub fn create(
        path: &Path,
        sample_names: Option<&[String]>,
        labeler: FeatureLabeler,
        umi_len: usize,
    ) -> anyhow::Result<Self> {
        anyhow::ensure!(
            (1..=32).contains(&umi_len),
            "UMI length {umi_len} is outside the supported range 1..=32"
        );
        let schema = schema(sample_names.is_some());
        let file =
            File::create(path).with_context(|| format!("could not create {}", path.display()))?;
        let props = WriterProperties::builder()
            // LZ4 (Parquet's LZ4_RAW), the codec of collated RAD chunks: on 10x
            // PBMC data it wrote ~2% faster than Snappy and loaded as fast,
            // at equal size (standard reference) or ~10% larger (USA mode).
            .set_compression(Compression::LZ4_RAW)
            // UMIs are close to unique within a row group, so a dictionary only
            // costs memory and space before Parquet falls back to plain encoding.
            .set_column_dictionary_enabled(ColumnPath::from("rep_umi"), false)
            .set_created_by(format!("alevin-fry {}", env!("CARGO_PKG_VERSION")))
            .build();
        let writer = ArrowWriter::try_new(
            BufWriter::with_capacity(4 << 20, file),
            schema.clone(),
            Some(props),
        )?;
        let (file, factory) = writer.into_serialized_writer()?;
        let strings = |v: &[String]| -> ArrayRef { Arc::new(StringArray::from_iter_values(v)) };
        Ok(Self {
            file: Mutex::new(Some(file)),
            factory,
            schema,
            feature_values: strings(&labeler.names),
            sample_values: sample_names.map(strings),
            status_values: Arc::new(StringArray::from_iter_values(STATUS_VALUES)),
            labeler,
            umi_len,
        })
    }

    /// Encode and compress `batch` as one row group, then append it. Only the
    /// append holds the lock.
    fn write_row_group(&self, batch: &RecordBatch) -> anyhow::Result<()> {
        let mut writers = self.factory.create_column_writers(0)?;
        let mut writer_iter = writers.iter_mut();
        for (array, field) in batch.columns().iter().zip(self.schema.fields()) {
            for leaves in compute_leaves(field, array)? {
                writer_iter
                    .next()
                    .context("molecule table column writers exhausted")?
                    .write(&leaves)?;
            }
        }
        let chunks = writers
            .into_iter()
            .map(|w| w.close())
            .collect::<parquet::errors::Result<Vec<_>>>()?;

        let mut guard = self
            .file
            .lock()
            .map_err(|_| anyhow::anyhow!("molecule table lock was poisoned"))?;
        let file = guard
            .as_mut()
            .context("molecule table was already finalized")?;
        let mut row_group = file.next_row_group()?;
        for chunk in chunks {
            chunk.append_to_row_group(&mut row_group)?;
        }
        row_group.close()?;
        Ok(())
    }

    /// Write the footer. Every worker must have flushed its [`MoleculeBatch`].
    /// Returns the number of rows in the table.
    pub fn finish(&self) -> anyhow::Result<u64> {
        let file = self
            .file
            .lock()
            .map_err(|_| anyhow::anyhow!("molecule table lock was poisoned"))?
            .take()
            .context("molecule table was already finalized")?;
        let metadata = file.close().context("could not finalize molecule table")?;
        Ok(metadata.file_metadata().num_rows() as u64)
    }
}

/// A worker's buffered rows; becomes one row group when full.
pub struct MoleculeBatch {
    out: Arc<MoleculeTableWriter>,
    /// one key per row into `cell_barcodes`, which holds each cell once
    cell_key: Int32Builder,
    cell_barcodes: StringBuilder,
    sample_key: Option<Int32Builder>,
    rep_umi: StringBuilder,
    n_umis: UInt32Builder,
    feature_key: Int32Builder,
    candidate_features: ListBuilder<StringBuilder>,
    reads: UInt32Builder,
    status_key: Int32Builder,
    rows: usize,
    umi_buf: Vec<u8>,
}

impl MoleculeBatch {
    pub fn new(out: Arc<MoleculeTableWriter>) -> Self {
        Self {
            cell_key: Int32Builder::new(),
            cell_barcodes: StringBuilder::new(),
            sample_key: out.sample_values.as_ref().map(|_| Int32Builder::new()),
            rep_umi: StringBuilder::new(),
            n_umis: UInt32Builder::new(),
            feature_key: Int32Builder::new(),
            candidate_features: ListBuilder::new(StringBuilder::new()),
            reads: UInt32Builder::new(),
            status_key: Int32Builder::new(),
            rows: 0,
            umi_buf: vec![0; out.umi_len],
            out,
        }
    }

    /// Append one cell's molecules. `sample` is the cell's sample index for
    /// multi-sample input; `em` is whether the cell's counts came from an EM
    /// strategy (see [`FeatureLabeler::assign`]).
    pub fn add_cell(
        &mut self,
        cell_barcode: &str,
        sample: Option<usize>,
        mols: &CellMolecules,
        em: bool,
    ) -> anyhow::Result<()> {
        if mols.is_empty() {
            return Ok(());
        }
        let cell_key = self.cell_barcodes.len() as i32;
        // Like quant's barcode output, a sample index without a name is left
        // unlabelled rather than failing the run.
        let num_samples = self.out.sample_values.as_ref().map_or(0, |v| v.len());
        let sample_key = sample.filter(|&s| s < num_samples).map(|s| s as i32);
        self.cell_barcodes.append_value(cell_barcode);
        let labeler = &self.out.labeler;
        for (umi, labels, reads, n_umis, low_support) in mols.iter() {
            self.cell_key.append_value(cell_key);
            if let Some(k) = self.sample_key.as_mut() {
                k.append_option(sample_key);
            }
            decode_umi(umi, &mut self.umi_buf);
            // decode_umi writes only ACGT, so this is valid UTF-8.
            self.rep_umi
                .append_value(std::str::from_utf8(&self.umi_buf).expect("ACGT is UTF-8"));
            self.n_umis.append_value(n_umis);
            self.reads.append_value(reads);
            let assignment = if low_support {
                match labels {
                    [l] => Assignment::Feature(labeler.label_column(*l)),
                    _ => Assignment::MultiGene,
                }
            } else {
                labeler.assign(labels, em)
            };
            self.status_key
                .append_value(match (low_support, &assignment) {
                    (true, _) => STATUS_KEY_LOW_SUPPORT,
                    (false, Assignment::Feature(_)) => STATUS_KEY_ASSIGNED,
                    (false, Assignment::MultiGene) => STATUS_KEY_MULTI_GENE,
                });
            match assignment {
                Assignment::Feature(col) => {
                    self.feature_key.append_value(col as i32);
                    self.candidate_features.append_null();
                }
                Assignment::MultiGene => {
                    self.feature_key.append_null();
                    for &l in labels {
                        self.candidate_features
                            .values()
                            .append_value(labeler.label_name(l));
                    }
                    self.candidate_features.append(true);
                }
            }
        }
        self.rows += mols.len();
        if self.rows >= ROW_GROUP_ROWS {
            self.flush()?;
        }
        Ok(())
    }

    /// Write any buffered rows as a (possibly short) row group.
    pub fn flush(&mut self) -> anyhow::Result<()> {
        if self.rows == 0 {
            return Ok(());
        }
        let dict = |keys: &mut Int32Builder, values: ArrayRef| -> anyhow::Result<ArrayRef> {
            Ok(Arc::new(Int32DictionaryArray::try_new(
                keys.finish(),
                values,
            )?))
        };
        let out = self.out.clone();
        let cell_values: ArrayRef = Arc::new(self.cell_barcodes.finish());
        let mut columns: Vec<ArrayRef> = vec![dict(&mut self.cell_key, cell_values)?];
        if let (Some(k), Some(v)) = (self.sample_key.as_mut(), out.sample_values.as_ref()) {
            columns.push(dict(k, v.clone())?);
        }
        columns.extend([
            Arc::new(self.rep_umi.finish()) as ArrayRef,
            Arc::new(self.n_umis.finish()),
            dict(&mut self.feature_key, out.feature_values.clone())?,
            Arc::new(self.candidate_features.finish()),
            Arc::new(self.reads.finish()),
            dict(&mut self.status_key, out.status_values.clone())?,
        ]);
        self.rows = 0;
        let batch = RecordBatch::try_new(out.schema.clone(), columns)?;
        out.write_row_group(&batch)
    }
}

/// Decode a 2-bit packed UMI (A=0, C=1, G=2, T=3; first base in the most
/// significant bits) into `out`, whose length is the UMI length.
fn decode_umi(umi: u64, out: &mut [u8]) {
    const NT: [u8; 4] = *b"ACGT";
    let n = out.len();
    for (i, b) in out.iter_mut().enumerate() {
        *b = NT[((umi >> (2 * (n - 1 - i))) & 3) as usize];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(v: &[&str]) -> Arc<Vec<String>> {
        Arc::new(v.iter().map(|s| s.to_string()).collect())
    }

    #[test]
    fn decodes_packed_umis() {
        let mut out = [0u8; 4];
        // A C G T
        decode_umi(0b00_01_10_11, &mut out);
        assert_eq!(&out, b"ACGT");
        let mut out = [0u8; 3];
        decode_umi(0b11_11_00, &mut out);
        assert_eq!(&out, b"TTA");
    }

    #[test]
    fn cell_molecules_round_trip() {
        let mut m = CellMolecules::default();
        m.push(7, &[1, 3], 5, 2);
        m.push_low_support(9, &[4], 1, 1);
        let rows: Vec<_> = m.iter().collect();
        assert_eq!(rows[0], (7, &[1u32, 3][..], 5, 2, false));
        assert_eq!(rows[1], (9, &[4u32][..], 1, 1, true));
        m.clear();
        assert!(m.is_empty());
    }

    #[test]
    fn standard_mode_assignment() {
        let l = FeatureLabeler::new(names(&["g0", "g1", "g2"]), false);
        assert_eq!(l.assign(&[2], false), Assignment::Feature(2));
        assert_eq!(l.assign(&[0, 1], false), Assignment::MultiGene);
        assert_eq!(l.assign(&[0, 1], true), Assignment::MultiGene);
        assert_eq!(l.label_name(1), "g1");
    }

    /// Write two cells (one per batch, so two row groups) and read them back.
    #[test]
    fn table_round_trips_through_parquet() {
        use arrow_array::Array;
        use arrow_array::cast::AsArray;
        use arrow_array::types::{Int32Type, UInt32Type};
        use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join(MOLECULE_TABLE_NAME);
        let samples = vec!["s0".to_string(), "s1".to_string()];
        let out = Arc::new(
            MoleculeTableWriter::create(
                &path,
                Some(&samples),
                FeatureLabeler::new(names(&["g0", "g1", "g2"]), false),
                4,
            )
            .unwrap(),
        );
        let mut mols = CellMolecules::default();
        let mut a = MoleculeBatch::new(out.clone());
        mols.push(0b00_01_10_11, &[2], 7, 1); // ACGT
        mols.push(0b11_11_11_11, &[0, 1], 3, 2); // TTTT
        a.add_cell("AAAC", Some(1), &mols, false).unwrap();
        let mut b = MoleculeBatch::new(out.clone());
        mols.clear();
        mols.push_low_support(0, &[1], 2, 1); // AAAA
        b.add_cell("CCCC", Some(9), &mols, false).unwrap(); // unnamed sample
        a.flush().unwrap();
        b.flush().unwrap();
        assert_eq!(out.finish().unwrap(), 3);

        let reader = ParquetRecordBatchReaderBuilder::try_new(File::open(&path).unwrap())
            .unwrap()
            .build()
            .unwrap();
        let batches: Vec<RecordBatch> = reader.map(Result::unwrap).collect();
        let batch = arrow_select::concat::concat_batches(&batches[0].schema(), &batches).unwrap();
        assert_eq!(batch.schema(), schema(true));
        let dict_strings = |name: &str| -> Vec<Option<String>> {
            let d = batch
                .column_by_name(name)
                .unwrap()
                .as_dictionary::<Int32Type>();
            let v = d.values().as_string::<i32>();
            d.keys()
                .iter()
                .map(|k| k.map(|k| v.value(k as usize).to_string()))
                .collect()
        };
        let s = |x: &str| Some(x.to_string());
        assert_eq!(
            dict_strings("cell_barcode"),
            [s("AAAC"), s("AAAC"), s("CCCC")]
        );
        assert_eq!(dict_strings("sample"), [s("s1"), s("s1"), None]);
        assert_eq!(dict_strings("feature"), [s("g2"), None, s("g1")]);
        assert_eq!(
            dict_strings("status"),
            [s("assigned"), s("multi_gene"), s("low_support")]
        );
        let umis: Vec<_> = batch
            .column_by_name("rep_umi")
            .unwrap()
            .as_string::<i32>()
            .iter()
            .collect();
        assert_eq!(umis, [Some("ACGT"), Some("TTTT"), Some("AAAA")]);
        let u32s = |name: &str| -> Vec<u32> {
            batch
                .column_by_name(name)
                .unwrap()
                .as_primitive::<UInt32Type>()
                .values()
                .to_vec()
        };
        assert_eq!(u32s("reads"), [7, 3, 2]);
        assert_eq!(u32s("n_umis"), [1, 2, 1]);
        let cands = batch
            .column_by_name("candidate_features")
            .unwrap()
            .as_list::<i32>();
        assert!(cands.is_null(0) && cands.is_null(2));
        let c1 = cands.value(1);
        let c1: Vec<_> = c1.as_string::<i32>().iter().flatten().collect();
        assert_eq!(c1, ["g0", "g1"]);
    }

    #[test]
    fn usa_mode_assignment_matches_the_matrix_rules() {
        // two genes: columns g0 g1 | g0-U g1-U | g0-A g1-A
        let l = FeatureLabeler::new(names(&["g0", "g1", "g0-U", "g1-U", "g0-A", "g1-A"]), true);
        // labels: g0 spliced = 0, g0 unspliced = 1, g1 spliced = 2, g1 unspliced = 3
        assert_eq!(l.assign(&[2], false), Assignment::Feature(1));
        assert_eq!(l.assign(&[3], false), Assignment::Feature(3));
        // splicing-ambiguous within one gene -> ambiguous column, EM or not
        assert_eq!(l.assign(&[0, 1], false), Assignment::Feature(4));
        assert_eq!(l.assign(&[0, 1], true), Assignment::Feature(4));
        // cross-gene, exactly one spliced: prefer-spliced for non-EM only
        assert_eq!(l.assign(&[1, 2], false), Assignment::Feature(1));
        assert_eq!(l.assign(&[1, 2], true), Assignment::MultiGene);
        // two spliced genes -> gene-ambiguous
        assert_eq!(l.assign(&[0, 2], false), Assignment::MultiGene);
        assert_eq!(l.label_name(3), "g1-U");
    }
}
