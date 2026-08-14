/*
 * Copyright (c) 2020-2026 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Square-bin aggregation for Visium HD.
//!
//! Visium HD reads out a grid of 2 µm squares, which is finer than most
//! analyses want. Space Ranger reports the same data at several bin sizes by
//! summing squares into larger squares, and names each bin
//! `s_{size_um:03}um_{row:05}_{col:05}`. Binning is integer division of the row
//! and column by the ratio of the two bin sizes
//! (`lib/rust/barcode/src/binned.rs` in spaceranger 4.0.1).
//!
//! This does that aggregation on a `quant` count matrix. What it cannot do is
//! invent the grid: alevin-fry sees nucleotide barcodes, and the map from a
//! barcode to its square is a property of the slide. Either the matrix rows are
//! already named as square bins, or a barcode-to-position file has to say where
//! each barcode sits.

use anyhow::{Context, bail};
use num_format::{Locale, ToFormattedString};
use serde_json::json;
use slog::info;
use std::collections::HashMap;
use std::fmt;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

use crate::prog_opts::BinSpatialOpts;

/// Prefix Space Ranger gives square-bin barcodes.
const SQUARE_BIN_PREFIX: &str = "s";
/// Suffix on the bin size field.
const MICROMETER: &str = "um";

/// A square bin on the Visium HD grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SquareBinIndex {
    /// Row of the bin, counting from the top.
    pub row: usize,
    /// Column of the bin, counting from the left.
    pub col: usize,
    /// Side of the bin in micrometres.
    pub size_um: u32,
}

impl SquareBinIndex {
    /// The bin containing this one, `bin_scale` times as wide on a side.
    ///
    /// This is Space Ranger's `binned`: integer division of the coordinates,
    /// which is what makes the bins tile without gaps or overlap.
    pub fn binned(self, bin_scale: u32) -> SquareBinIndex {
        SquareBinIndex {
            row: self.row / bin_scale as usize,
            col: self.col / bin_scale as usize,
            size_um: self.size_um * bin_scale,
        }
    }
}

impl fmt::Display for SquareBinIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{SQUARE_BIN_PREFIX}_{:03}{MICROMETER}_{:05}_{:05}",
            self.size_um, self.row, self.col
        )
    }
}

impl FromStr for SquareBinIndex {
    type Err = anyhow::Error;

    /// Parses `s_002um_00123_00456`, with or without a trailing `-<gem group>`.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let body = match s.split_once('-') {
            Some((barcode, _gem_group)) => barcode,
            None => s,
        };
        let mut parts = body.split('_');

        let prefix = parts.next().unwrap_or_default();
        if prefix != SQUARE_BIN_PREFIX {
            bail!(
                "square-bin barcodes start with `{SQUARE_BIN_PREFIX}_`; `{s}` starts with `{prefix}`"
            );
        }
        let size_um: u32 = parts
            .next()
            .and_then(|p| p.strip_suffix(MICROMETER))
            .context("could not read the bin size")
            .and_then(|p| p.parse().context("could not parse the bin size"))
            .with_context(|| format!("while parsing `{s}` as a square-bin barcode"))?;
        let row: usize = parts
            .next()
            .context("could not read the row")
            .and_then(|p| p.parse().context("could not parse the row"))
            .with_context(|| format!("while parsing `{s}` as a square-bin barcode"))?;
        let col: usize = parts
            .next()
            .context("could not read the column")
            .and_then(|p| p.parse().context("could not parse the column"))
            .with_context(|| format!("while parsing `{s}` as a square-bin barcode"))?;

        Ok(SquareBinIndex { row, col, size_um })
    }
}

/// Read a barcode-to-position file: `barcode<TAB>row<TAB>col`, one per line,
/// `#` comments and blank lines skipped. Extra columns are ignored, so a slide
/// layout carrying pixel coordinates as well can be passed through unchanged.
fn read_positions(path: &Path) -> anyhow::Result<HashMap<String, (usize, usize)>> {
    let f =
        std::fs::File::open(path).with_context(|| format!("could not open {}", path.display()))?;
    let mut out = HashMap::new();
    for (lineno, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split(['\t', ',']);
        let (Some(bc), Some(row), Some(col)) = (it.next(), it.next(), it.next()) else {
            bail!(
                "{}:{}: expected at least three fields (barcode, row, column), found `{line}`",
                path.display(),
                lineno + 1
            );
        };
        let row: usize = row
            .trim()
            .parse()
            .with_context(|| format!("{}:{}: bad row `{row}`", path.display(), lineno + 1))?;
        let col: usize = col
            .trim()
            .parse()
            .with_context(|| format!("{}:{}: bad column `{col}`", path.display(), lineno + 1))?;
        out.insert(bc.trim().to_owned(), (row, col));
    }
    Ok(out)
}

/// Sum a `quant` matrix into square bins `bin_size_um` on a side.
pub fn bin_spatial(opts: BinSpatialOpts) -> anyhow::Result<usize> {
    let log = opts.log;
    let alevin_dir = opts.input_dir.join("alevin");
    let rows_path = alevin_dir.join("quants_mat_rows.txt");
    let cols_path = alevin_dir.join("quants_mat_cols.txt");
    let mtx_path = alevin_dir.join("quants_mat.mtx");

    let barcodes: Vec<String> = BufReader::new(
        std::fs::File::open(&rows_path)
            .with_context(|| format!("could not open {}", rows_path.display()))?,
    )
    .lines()
    .collect::<Result<Vec<_>, _>>()?
    .into_iter()
    .map(|l| l.trim().to_owned())
    .filter(|l| !l.is_empty())
    .collect();

    // Where each barcode sits on the grid. Either the row names are already
    // square-bin barcodes, or a positions file says.
    let positions = match &opts.barcode_positions {
        Some(p) => Some(read_positions(p)?),
        None => None,
    };

    let source_bins: Vec<SquareBinIndex> = match &positions {
        None => {
            let parsed: Result<Vec<_>, _> = barcodes
                .iter()
                .map(|b| b.parse::<SquareBinIndex>())
                .collect();
            parsed.context(
                "the matrix rows are not square-bin barcodes, so --barcode-positions is required \
                 to say where each barcode sits on the slide",
            )?
        }
        Some(pos) => {
            let mut v = Vec::with_capacity(barcodes.len());
            let mut missing = 0usize;
            let mut first_missing = String::new();
            for b in &barcodes {
                match pos.get(b) {
                    Some(&(row, col)) => v.push(SquareBinIndex {
                        row,
                        col,
                        size_um: opts.source_size_um,
                    }),
                    None => {
                        if missing == 0 {
                            first_missing = b.clone();
                        }
                        missing += 1;
                        // placed out of grid; filtered out below
                        v.push(SquareBinIndex {
                            row: usize::MAX,
                            col: usize::MAX,
                            size_um: opts.source_size_um,
                        });
                    }
                }
            }
            if missing > 0 {
                info!(
                    log,
                    "{} of {} barcodes had no position and were dropped (first: {})",
                    missing.to_formatted_string(&Locale::en),
                    barcodes.len().to_formatted_string(&Locale::en),
                    first_missing
                );
            }
            v
        }
    };

    // Bin size ratio has to be a whole number: bins tile by integer division.
    let source_size_um = source_bins
        .iter()
        .find(|b| b.row != usize::MAX)
        .map(|b| b.size_um)
        .unwrap_or(opts.source_size_um);
    if !opts.bin_size_um.is_multiple_of(source_size_um) {
        bail!(
            "--bin-size-um {} is not a whole multiple of the source bin size {} um",
            opts.bin_size_um,
            source_size_um
        );
    }
    let bin_scale = opts.bin_size_um / source_size_um;
    if bin_scale == 0 {
        bail!(
            "--bin-size-um {} is smaller than the source bin size {} um; binning only aggregates",
            opts.bin_size_um,
            source_size_um
        );
    }

    // Assign each source row to its destination bin.
    let mut bin_of_row: Vec<Option<usize>> = vec![None; barcodes.len()];
    let mut bin_ids: HashMap<SquareBinIndex, usize> = HashMap::new();
    let mut bins: Vec<SquareBinIndex> = Vec::new();
    for (i, sb) in source_bins.iter().enumerate() {
        if sb.row == usize::MAX {
            continue;
        }
        let dest = sb.binned(bin_scale);
        let id = *bin_ids.entry(dest).or_insert_with(|| {
            bins.push(dest);
            bins.len() - 1
        });
        bin_of_row[i] = Some(id);
    }

    // Bins are emitted in grid order, which does not depend on how the input
    // happened to be ordered.
    let mut order: Vec<usize> = (0..bins.len()).collect();
    order.sort_unstable_by_key(|&i| bins[i]);
    let mut new_id = vec![0usize; bins.len()];
    for (rank, &i) in order.iter().enumerate() {
        new_id[i] = rank;
    }
    let sorted_bins: Vec<SquareBinIndex> = order.iter().map(|&i| bins[i]).collect();

    info!(
        log,
        "binning {} source squares of {} um into {} squares of {} um (scale {})",
        barcodes.len().to_formatted_string(&Locale::en),
        source_size_um,
        sorted_bins.len().to_formatted_string(&Locale::en),
        opts.bin_size_um,
        bin_scale
    );

    // Sum the matrix.
    let f = std::fs::File::open(&mtx_path)
        .with_context(|| format!("could not open {}", mtx_path.display()))?;
    let mut lines = BufReader::new(f).lines();
    let mut header = lines
        .next()
        .transpose()?
        .context("matrix market file was empty")?;
    while header.starts_with('%') {
        header = lines
            .next()
            .transpose()?
            .context("matrix market file had no dimension line")?;
    }
    let dims: Vec<usize> = header
        .split_whitespace()
        .map(str::parse)
        .collect::<Result<Vec<_>, _>>()
        .context("could not parse the matrix market dimension line")?;
    if dims.len() != 3 {
        bail!("expected 3 fields on the matrix market dimension line, found {dims:?}");
    }
    let (num_rows_in, num_features) = (dims[0], dims[1]);
    if num_rows_in != barcodes.len() {
        bail!(
            "the matrix has {} rows but quants_mat_rows.txt lists {} barcodes",
            num_rows_in,
            barcodes.len()
        );
    }

    let mut acc: HashMap<(usize, u32), f64> = HashMap::new();
    let mut dropped = 0usize;
    for line in lines {
        let line = line?;
        let mut it = line.split_whitespace();
        let (Some(i), Some(j), Some(v)) = (it.next(), it.next(), it.next()) else {
            continue;
        };
        let i: usize = i.parse()?;
        let j: u32 = j.parse()?;
        let v: f64 = v.parse()?;
        match bin_of_row[i - 1] {
            Some(b) => *acc.entry((new_id[b], j)).or_insert(0.0) += v,
            None => dropped += 1,
        }
    }
    if dropped > 0 {
        info!(
            log,
            "dropped {} matrix entries belonging to barcodes with no position",
            dropped.to_formatted_string(&Locale::en)
        );
    }

    // Write the binned matrix.
    let out_alevin = opts.output_dir.join("alevin");
    std::fs::create_dir_all(&out_alevin)
        .with_context(|| format!("could not create {}", out_alevin.display()))?;

    let mut entries: Vec<((usize, u32), f64)> = acc.into_iter().collect();
    entries.sort_unstable_by_key(|&((r, c), _)| (r, c));

    let mut trimat =
        sprs::TriMatI::<f32, u32>::with_capacity((sorted_bins.len(), num_features), entries.len());
    let mut total = 0.0f64;
    for ((row, col), val) in &entries {
        trimat.add_triplet(*row, (*col - 1) as usize, *val as f32);
        total += *val;
    }
    sprs::io::write_matrix_market(out_alevin.join("quants_mat.mtx"), &trimat)?;

    {
        let mut w = BufWriter::new(std::fs::File::create(
            out_alevin.join("quants_mat_rows.txt"),
        )?);
        for b in &sorted_bins {
            writeln!(w, "{b}")?;
        }
    }
    if cols_path.is_file() {
        std::fs::copy(&cols_path, out_alevin.join("quants_mat_cols.txt"))
            .with_context(|| format!("could not copy {}", cols_path.display()))?;
    }

    let meta = json!({
        "cmd": opts.cmdline,
        "version_str": opts.version,
        "source_bins": barcodes.len(),
        "source_size_um": source_size_um,
        "bin_size_um": opts.bin_size_um,
        "bin_scale": bin_scale,
        "num_bins": sorted_bins.len(),
        "num_features": num_features,
        "nonzeros": entries.len(),
        "total_counts": total,
        "dropped_entries": dropped,
    });
    std::fs::write(
        opts.output_dir.join("bin_spatial.json"),
        serde_json::to_string_pretty(&meta)?,
    )?;

    info!(
        log,
        "wrote {} bins x {} features, {} nonzeros, {} total counts",
        sorted_bins.len().to_formatted_string(&Locale::en),
        num_features.to_formatted_string(&Locale::en),
        entries.len().to_formatted_string(&Locale::en),
        (total.round() as u64).to_formatted_string(&Locale::en)
    );
    Ok(sorted_bins.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_and_formats_like_space_ranger() {
        // the assertions from spaceranger's own `test_square_bin`
        let b = SquareBinIndex {
            row: 1,
            col: 2,
            size_um: 16,
        };
        assert_eq!(b.to_string(), "s_016um_00001_00002");
        assert_eq!(b, "s_016um_00001_00002".parse().unwrap());
        assert_eq!(b, "s_016um_00001_00002-1".parse().unwrap());
    }

    #[test]
    fn rejects_non_spatial_barcodes() {
        assert!("ACGTACGTACGTACGT".parse::<SquareBinIndex>().is_err());
        assert!("s_016um_00001".parse::<SquareBinIndex>().is_err());
        assert!("s_16_00001_00002".parse::<SquareBinIndex>().is_err());
    }

    #[test]
    fn binning_is_integer_division() {
        let b = SquareBinIndex {
            row: 7,
            col: 9,
            size_um: 2,
        };
        // 2 um squares into 8 um squares: scale 4
        let big = b.binned(4);
        assert_eq!(
            big,
            SquareBinIndex {
                row: 1,
                col: 2,
                size_um: 8
            }
        );
        // the whole 4x4 block lands in the same bin
        for r in 4..8 {
            for c in 8..12 {
                let cell = SquareBinIndex {
                    row: r,
                    col: c,
                    size_um: 2,
                };
                assert_eq!(cell.binned(4), big);
            }
        }
        // and the neighbouring block does not
        assert_ne!(
            SquareBinIndex {
                row: 8,
                col: 9,
                size_um: 2
            }
            .binned(4),
            big
        );
    }

    #[test]
    fn bins_sort_in_grid_order() {
        let mut v = [
            SquareBinIndex {
                row: 1,
                col: 5,
                size_um: 8,
            },
            SquareBinIndex {
                row: 0,
                col: 9,
                size_um: 8,
            },
            SquareBinIndex {
                row: 1,
                col: 2,
                size_um: 8,
            },
        ];
        v.sort_unstable();
        assert_eq!(v[0].row, 0);
        assert_eq!((v[1].row, v[1].col), (1, 2));
        assert_eq!((v[2].row, v[2].col), (1, 5));
    }
}
