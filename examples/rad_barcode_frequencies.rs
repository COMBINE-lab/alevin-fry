/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Materialize raw, orientation-compatible cell-barcode frequencies from RAD.
//!
//! This is a developer measurement utility, not a supported user command.

use alevin_fry::utils::write_permit_list_freq;
use anyhow::{Context, bail};
use bio_types::strand::Strand;
use dashmap::DashMap;
use libradicl::header::RadPrelude;
use libradicl::record::{
    AlevinFryReadRecord, CollatableMappedRecord, ConvertiblePrimitiveInteger, KnownSize,
    MappedRecord, RecordContext,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Arc;

fn count_barcodes<R, B>(
    mut reader: BufReader<File>,
    prelude: &RadPrelude,
    expected: Strand,
    threads: usize,
) -> anyhow::Result<HashMap<u64, u64, ahash::RandomState>>
where
    B: ConvertiblePrimitiveInteger,
    u64: From<B>,
    R: MappedRecord + CollatableMappedRecord<B> + KnownSize,
    <R as MappedRecord>::ParsingContext: RecordContext + Clone + Send,
{
    let counts = Arc::new(DashMap::<u64, u64, ahash::RandomState>::default());
    let mut chunk_reader = libradicl::readers::ParallelChunkReader::<R>::new(
        prelude,
        NonZeroUsize::new(threads).unwrap(),
    );
    std::thread::scope(|scope| -> anyhow::Result<()> {
        let mut handles = Vec::with_capacity(threads);
        for _ in 0..threads {
            let chunks = chunk_reader.chunk_iter();
            let counts = counts.clone();
            handles.push(scope.spawn(move || {
                for meta_chunk in chunks {
                    for chunk in meta_chunk.iter() {
                        for record in &chunk.reads {
                            if record.has_alignment_on_strand(expected) {
                                let barcode = u64::from(record.collate_key());
                                *counts.entry(barcode).or_insert(0) += 1;
                            }
                        }
                    }
                }
            }));
        }
        chunk_reader.start(&mut reader, None::<fn(u64, u64)>)?;
        for handle in handles {
            handle.join().expect("frequency worker panicked");
        }
        Ok(())
    })?;
    Ok(Arc::into_inner(counts)
        .context("frequency count map still has multiple owners")?
        .into_iter()
        .collect())
}

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args_os().skip(1);
    let rad = PathBuf::from(args.next().context("missing map.rad argument")?);
    let output = PathBuf::from(args.next().context("missing output frequency file")?);
    let expected = match args
        .next()
        .context("missing orientation (fw, rc, or both)")?
        .to_string_lossy()
        .as_ref()
    {
        "fw" => Strand::Forward,
        "rc" => Strand::Reverse,
        "both" | "either" => Strand::Unknown,
        other => bail!("unknown orientation '{other}'; expected fw, rc, or both"),
    };
    let threads = args
        .next()
        .map(|value| value.to_string_lossy().parse::<usize>())
        .transpose()?
        .unwrap_or(16)
        .max(1);
    if args.next().is_some() {
        bail!("usage: rad_barcode_frequencies <map.rad> <output.bin> <fw|rc|both> [threads]");
    }

    let file = File::open(&rad).with_context(|| format!("could not open {}", rad.display()))?;
    let mut reader = BufReader::new(file);
    let prelude = RadPrelude::from_bytes(&mut reader)?;
    let file_tags = prelude.file_tags.parse_tags_from_bytes(&mut reader)?;
    if let Some(value) = file_tags.get("num_barcodes") {
        let count: u16 = value
            .try_into()
            .context("RAD num_barcodes tag does not fit in u16")?;
        if count > 1 {
            bail!("use multi_correction_plan_report for multi-barcode RAD");
        }
    }
    if prelude.aln_tags.has_tag("pos")
        || (prelude.aln_tags.has_tag("as")
            && prelude.aln_tags.has_tag("start")
            && prelude.aln_tags.has_tag("end"))
    {
        bail!("this report utility currently expects the standard short-read RNA RAD record");
    }
    let barcode_len: u16 = file_tags
        .get("cblen")
        .context("RAD file has no cblen tag")?
        .try_into()
        .context("RAD cblen tag does not fit in u16")?;
    let counts = count_barcodes::<AlevinFryReadRecord, u64>(reader, &prelude, expected, threads)?;
    write_permit_list_freq(&output, barcode_len, &counts)
        .map_err(|error| anyhow::anyhow!("could not write {}: {error}", output.display()))?;
    eprintln!(
        "wrote {} distinct barcode counts to {}",
        counts.len(),
        output.display()
    );
    Ok(())
}
