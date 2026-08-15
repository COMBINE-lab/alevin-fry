/*
 * Copyright (c) 2020-2026 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Development utility for checking the compact correction handoff against
//! the historical single-barcode map.

use alevin_fry::correction_plan::{CORRECTION_PLAN_FILENAME, CorrectionPlan};
use anyhow::Context;
use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    let directory = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .context("usage: compare_correction_artifacts GPL_DIRECTORY")?;
    let legacy: HashMap<u64, u64> =
        bincode::deserialize_from(File::open(directory.join("permit_map.bin"))?)?;
    let plan = CorrectionPlan::read_from(&directory.join(CORRECTION_PLAN_FILENAME))?;
    let scope = plan
        .cell_scopes
        .iter()
        .find(|scope| scope.sample_barcode.is_none())
        .context("compact plan has no global cell scope")?;
    let compact: HashMap<_, _> = scope
        .corrections
        .iter()
        .map(|entry| (entry.observed, entry.corrected))
        .collect();

    let legacy_only = legacy
        .iter()
        .filter(|(observed, corrected)| compact.get(observed) != Some(corrected))
        .count();
    let compact_only = compact
        .iter()
        .filter(|(observed, corrected)| legacy.get(observed) != Some(corrected))
        .count();
    println!(
        "legacy={} compact={} legacy_only_or_changed={} compact_only_or_changed={}",
        legacy.len(),
        compact.len(),
        legacy_only,
        compact_only
    );
    Ok(())
}
