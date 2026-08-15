/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Stable handoff of GPL correction decisions to collation and ATAC sorting.

use crate::barcode_correction::{CorrectionSpec, CorrectionTable};
use anyhow::{Context, bail};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

pub const CORRECTION_PLAN_FILENAME: &str = "correction_plan.bin";
const MAGIC: &[u8; 8] = b"AFCORR\0\0";
const FORMAT_VERSION: u16 = 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct BarcodeCorrection {
    pub observed: u64,
    pub corrected: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CellCorrectionScope {
    /// `None` for an ordinary single-barcode experiment; otherwise the
    /// canonical sample barcode.
    pub sample_barcode: Option<u64>,
    pub spec: CorrectionSpec,
    pub corrections: Vec<BarcodeCorrection>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CorrectionPlan {
    pub sample_barcode_len: Option<u8>,
    pub cell_barcode_len: u8,
    pub sample_spec: Option<CorrectionSpec>,
    pub sample_corrections: Vec<BarcodeCorrection>,
    pub cell_scopes: Vec<CellCorrectionScope>,
}

impl CorrectionPlan {
    pub fn single(spec: CorrectionSpec, table: &CorrectionTable<u64>) -> Self {
        Self {
            sample_barcode_len: None,
            cell_barcode_len: spec.barcode_len,
            sample_spec: None,
            sample_corrections: Vec::new(),
            cell_scopes: vec![CellCorrectionScope {
                sample_barcode: None,
                spec,
                corrections: table
                    .entries()
                    .iter()
                    .map(|&(observed, corrected)| BarcodeCorrection {
                        observed,
                        corrected,
                    })
                    .collect(),
            }],
        }
    }

    pub fn validate(&self) -> anyhow::Result<()> {
        validate_len(self.cell_barcode_len, "cell")?;
        if let Some(len) = self.sample_barcode_len {
            validate_len(len, "sample")?;
        }
        match self.sample_barcode_len {
            None => {
                if self.sample_spec.is_some() || !self.sample_corrections.is_empty() {
                    bail!("single-barcode correction plan contains sample-correction information");
                }
                if self.cell_scopes.len() != 1 || self.cell_scopes[0].sample_barcode.is_some() {
                    bail!(
                        "single-barcode correction plan must contain exactly one global cell scope"
                    );
                }
            }
            Some(_) => {
                if self.cell_scopes.is_empty()
                    || self
                        .cell_scopes
                        .iter()
                        .any(|scope| scope.sample_barcode.is_none())
                {
                    bail!(
                        "multi-barcode correction plan must contain sample-scoped cell corrections"
                    );
                }
            }
        }
        if let Some(spec) = self.sample_spec {
            spec.validate()?;
            if Some(spec.barcode_len) != self.sample_barcode_len {
                bail!("correction plan sample specification has the wrong barcode length");
            }
        }
        validate_corrections(
            &self.sample_corrections,
            self.sample_barcode_len.unwrap_or(1),
            "sample",
        )?;

        let mut previous_scope = None;
        for scope in &self.cell_scopes {
            scope.spec.validate()?;
            if scope.spec.barcode_len != self.cell_barcode_len {
                bail!("correction plan cell scope has the wrong barcode length");
            }
            if let Some(previous) = previous_scope
                && previous >= scope.sample_barcode
            {
                bail!("correction plan cell scopes are not strictly ordered");
            }
            if let (Some(sample_barcode), Some(sample_len)) =
                (scope.sample_barcode, self.sample_barcode_len)
            {
                validate_barcode(sample_barcode, sample_len, "sample scope")?;
            }
            previous_scope = Some(scope.sample_barcode);
            validate_corrections(&scope.corrections, self.cell_barcode_len, "cell")?;
        }
        Ok(())
    }

    pub fn write_to(&mut self, output_dir: &Path) -> anyhow::Result<PathBuf> {
        self.sample_corrections
            .sort_unstable_by_key(|entry| entry.observed);
        self.cell_scopes
            .sort_unstable_by_key(|scope| scope.sample_barcode);
        for scope in &mut self.cell_scopes {
            scope
                .corrections
                .sort_unstable_by_key(|entry| entry.observed);
        }
        self.validate()?;

        std::fs::create_dir_all(output_dir).with_context(|| {
            format!(
                "could not create correction plan directory {}",
                output_dir.display()
            )
        })?;
        let path = output_dir.join(CORRECTION_PLAN_FILENAME);
        let temporary = output_dir.join(format!(".{CORRECTION_PLAN_FILENAME}.tmp"));
        {
            let file = File::create(&temporary).with_context(|| {
                format!("could not create correction plan {}", temporary.display())
            })?;
            let mut writer = BufWriter::new(file);
            writer.write_all(MAGIC)?;
            writer.write_all(&FORMAT_VERSION.to_le_bytes())?;
            bincode::serialize_into(&mut writer, self)
                .context("could not serialize correction plan")?;
            writer.flush()?;
        }
        std::fs::rename(&temporary, &path).with_context(|| {
            format!(
                "could not install correction plan {} as {}",
                temporary.display(),
                path.display()
            )
        })?;
        Ok(path)
    }

    pub fn read_from(path: &Path) -> anyhow::Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("could not open correction plan {}", path.display()))?;
        let mut reader = BufReader::new(file);
        let mut magic = [0u8; MAGIC.len()];
        reader.read_exact(&mut magic).with_context(|| {
            format!("correction plan {} has a truncated header", path.display())
        })?;
        if &magic != MAGIC {
            bail!("correction plan {} has invalid magic", path.display());
        }
        let mut version = [0u8; 2];
        reader.read_exact(&mut version)?;
        let version = u16::from_le_bytes(version);
        if version != FORMAT_VERSION {
            bail!(
                "correction plan {} uses unsupported format version {} (expected {})",
                path.display(),
                version,
                FORMAT_VERSION
            );
        }
        let plan: Self = bincode::deserialize_from(&mut reader)
            .context("could not deserialize correction plan")?;
        let mut trailing = [0_u8; 1];
        if reader.read(&mut trailing)? != 0 {
            bail!("correction plan {} contains trailing data", path.display());
        }
        plan.validate()?;
        Ok(plan)
    }
}

fn validate_len(len: u8, kind: &str) -> anyhow::Result<()> {
    if !(1..=32).contains(&len) {
        bail!("correction plan declares invalid {kind} barcode length {len}");
    }
    Ok(())
}

fn validate_corrections(
    corrections: &[BarcodeCorrection],
    barcode_len: u8,
    kind: &str,
) -> anyhow::Result<()> {
    let limit = if barcode_len == 32 {
        None
    } else {
        Some(1u64 << (2 * barcode_len))
    };
    let mut previous = None;
    for entry in corrections {
        if let Some(limit) = limit
            && (entry.observed >= limit || entry.corrected >= limit)
        {
            bail!("correction plan contains a {kind} barcode outside its declared length");
        }
        if previous.is_some_and(|value| value >= entry.observed) {
            bail!("correction plan {kind} corrections are not strictly ordered");
        }
        previous = Some(entry.observed);
    }
    Ok(())
}

fn validate_barcode(barcode: u64, barcode_len: u8, kind: &str) -> anyhow::Result<()> {
    if barcode_len < 32 && barcode >= (1_u64 << (2 * barcode_len)) {
        bail!("correction plan contains a {kind} barcode outside its declared length");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::barcode_correction::{
        BarcodeNeighborhood, BarcodeResolution, CorrectionIndex, RetainedSource,
    };

    #[test]
    fn deterministic_round_trip() {
        let spec = CorrectionSpec {
            barcode_len: 2,
            neighborhood: BarcodeNeighborhood::HammingOne,
            resolution: BarcodeResolution::Unique,
        };
        let index = CorrectionIndex::new(
            spec,
            [RetainedSource {
                source: 1,
                target: 1,
                exact_count: 3,
            }],
        )
        .unwrap();
        let table = index.compile_observed([(0, 2), (1, 3)]);
        let mut plan = CorrectionPlan::single(spec, &table);
        plan.cell_scopes[0].corrections.reverse();

        let dir = tempfile::tempdir().unwrap();
        let path = plan.write_to(dir.path()).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        let loaded = CorrectionPlan::read_from(&path).unwrap();
        assert_eq!(loaded, plan);

        loaded.clone().write_to(dir.path()).unwrap();
        assert_eq!(bytes, std::fs::read(path).unwrap());
    }

    #[test]
    fn present_malformed_artifact_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join(CORRECTION_PLAN_FILENAME);
        std::fs::write(&path, b"not a correction plan").unwrap();
        assert!(CorrectionPlan::read_from(&path).is_err());
    }

    #[test]
    fn mixed_scope_shapes_and_trailing_data_are_errors() {
        let spec = CorrectionSpec {
            barcode_len: 2,
            neighborhood: BarcodeNeighborhood::HammingOne,
            resolution: BarcodeResolution::Unique,
        };
        let malformed = CorrectionPlan {
            sample_barcode_len: Some(2),
            cell_barcode_len: 2,
            sample_spec: None,
            sample_corrections: Vec::new(),
            cell_scopes: vec![CellCorrectionScope {
                sample_barcode: None,
                spec,
                corrections: Vec::new(),
            }],
        };
        assert!(malformed.validate().is_err());

        let index = CorrectionIndex::new(
            spec,
            [RetainedSource {
                source: 1,
                target: 1,
                exact_count: 1,
            }],
        )
        .unwrap();
        let mut valid = CorrectionPlan::single(spec, &index.compile_observed([(1, 1)]));
        let dir = tempfile::tempdir().unwrap();
        let path = valid.write_to(dir.path()).unwrap();
        let mut bytes = std::fs::read(&path).unwrap();
        bytes.push(0);
        std::fs::write(&path, bytes).unwrap();
        assert!(CorrectionPlan::read_from(&path).is_err());
    }
}
