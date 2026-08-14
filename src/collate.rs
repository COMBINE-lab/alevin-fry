/*
 * Copyright (c) 2020-2026 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use ahash::AHashMap;
use anyhow::{Context, anyhow};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use slog::{crit, debug, info};
//use anyhow::{anyhow, Result};
use crate::constants as afconst;
use crate::utils::InternalVersionInfo;
use bio_types::strand::{Strand, StrandError};
use crossbeam_queue::ArrayQueue;
// use dashmap::DashMap;

use libradicl::chunk;
use libradicl::collation::{CollationManifest, SampleGroup};
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::multi_collation::{
    MultiBarcodeCollationOptions, MultiBarcodeCollationPlan, MultiBarcodeSampleCorrection,
    collate_multi_barcode,
};
use libradicl::rad_types::{self, RadIntId};
use libradicl::record::{
    AlevinFryReadRecordWithPositionT, AlevinFryRecordContext, CollatableMappedRecord,
    ConvertiblePrimitiveInteger, KnownSize, MappedRecord, MultiBarcodeRecordContext,
    ScLongReadRecordContext, ScLongReadRecordT,
};
use libradicl::schema::TempCellInfo;
use libradicl::single_collation::{
    SingleBarcodeCollationOptions, SingleBarcodeCollationPlan, collate_single_barcode,
};

use num_format::{Locale, ToFormattedString};
use scroll::{Pread, Pwrite};
use serde_json::json;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::{BufWriter, Cursor, Read, Seek, Write};
use std::iter::FromIterator;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use crate::utils as afutils;
use crate::utils::KnownRecordType;

#[allow(clippy::too_many_arguments)]
pub fn collate<P1, P2>(
    input_dir: P1,
    rad_dir: P2,
    num_threads: u32,
    max_records: u32,
    compress_out: bool,
    cmdline: &str,
    version_str: &str,
    //expected_ori: Strand,
    log: &slog::Logger,
) -> anyhow::Result<()>
where
    P1: Into<PathBuf>,
    P2: AsRef<Path>,
{
    collate_with_memory_limit(
        input_dir,
        rad_dir,
        num_threads,
        max_records,
        None,
        compress_out,
        cmdline,
        version_str,
        log,
    )
}

/// Collate with an optional byte-based working-memory budget.
///
/// When `memory_limit` is `None`, the historical `max_records` value is
/// retained and translated to a byte budget for the libradicl engine. When a
/// byte limit is supplied, it becomes authoritative and `max_records` is
/// derived for compatibility with collation paths that still use that knob.
#[allow(clippy::too_many_arguments)]
pub fn collate_with_memory_limit<P1, P2>(
    input_dir: P1,
    rad_dir: P2,
    num_threads: u32,
    max_records: u32,
    memory_limit: Option<u64>,
    compress_out: bool,
    cmdline: &str,
    version_str: &str,
    log: &slog::Logger,
) -> anyhow::Result<()>
where
    P1: Into<PathBuf>,
    P2: AsRef<Path>,
{
    const BYTES_PER_LEGACY_RECORD: u64 = 72;
    const MIN_MEMORY_LIMIT: u64 = 256 * 1024 * 1024;
    let memory_budget_bytes = memory_limit
        .unwrap_or_else(|| u64::from(max_records).saturating_mul(BYTES_PER_LEGACY_RECORD))
        .max(MIN_MEMORY_LIMIT);
    let max_records = if memory_limit.is_some() {
        (memory_budget_bytes.min(2 * 1024 * 1024 * 1024) / BYTES_PER_LEGACY_RECORD)
            .clamp(1, u64::from(u32::MAX)) as u32
    } else {
        max_records
    };

    let input_dir = input_dir.into();
    let parent = std::path::Path::new(input_dir.as_path());

    // open the metadata file and read the json
    let gpl_path = parent.join("generate_permit_list.json");
    let meta_data_file = File::open(&gpl_path)
        .with_context(|| format!("Could not open the file {:?}.", gpl_path.display()))?;
    let mdata: serde_json::Value = serde_json::from_reader(meta_data_file)?;

    let calling_version = InternalVersionInfo::from_str(version_str)?;
    let vd: InternalVersionInfo;
    match mdata.get("version_str") {
        Some(vs) => match vs.as_str() {
            Some(s) => {
                vd = InternalVersionInfo::from_str(s)?;
            }
            None => {
                return Err(anyhow!("The version_str field must be a string"));
            }
        },
        None => {
            return Err(anyhow!(
                "The generate_permit_list.json file does not contain a version_str field. Please re-run the generate-permit-list step with a newer version of alevin-fry"
            ));
        }
    };

    if let Err(es) = calling_version.is_compatible_with(&vd) {
        return Err(anyhow!(es));
    }

    // Check if this is multi-barcode mode (no global permit_freq.bin,
    // but sample_info.json exists)
    let is_multi_barcode = mdata
        .get("multi_barcode")
        .and_then(|v| v.as_bool())
        .unwrap_or(false);

    if is_multi_barcode {
        // Multi-barcode: skip permit_freq.bin loading, go directly to collate_with_temp
        // which will detect the multi-barcode record type and use per-sample permit maps.
        // total_to_collate is computed from sample_info.json
        let sample_info_file = File::open(parent.join("sample_info.json"))
            .context("couldn't open sample_info.json for multi-barcode collation")?;
        let sample_info: serde_json::Value = serde_json::from_reader(sample_info_file)?;
        let total_to_collate = sample_info["matched_reads"].as_u64().unwrap_or(0);

        info!(
            log,
            "Multi-barcode mode: {} total reads to collate", total_to_collate
        );

        // tsv_map is empty — do_collate_multi_bc will build its own from per-sample data
        return collate_with_temp_memory(
            input_dir,
            rad_dir,
            num_threads,
            max_records,
            memory_budget_bytes,
            Vec::new(), // tsv_map not used for multi-barcode
            total_to_collate,
            compress_out,
            cmdline,
            version_str,
            log,
        );
    }

    // if only an *old* version of the permit_freq is present, then complain and exit
    if parent.join("permit_freq.tsv").exists() && !parent.join("permit_freq.bin").exists() {
        crit!(
            log,
            "The file permit_freq.bin doesn't exist, please rerun alevin-fry generate-permit-list command."
        );
        // std::process::exit(1);
        return Err(anyhow!("execution terminated unexpectedly"));
    }

    // open file
    let freq_file =
        std::fs::File::open(parent.join("permit_freq.bin")).context("couldn't open file")?;

    // header buffer
    let mut rbuf = [0u8; 8];

    // read header
    let mut rdr = BufReader::new(&freq_file);
    rdr.read_exact(&mut rbuf)
        .context("couldn't read freq file header")?;
    let freq_file_version = rbuf
        .pread::<u64>(0)
        .context("couldn't read freq file version")?;
    // make sure versions match
    if freq_file_version > afconst::PERMIT_FILE_VER {
        crit!(
            log,
            "The permit_freq.bin file had version {}, but this version of alevin-fry requires version {}",
            freq_file_version,
            afconst::PERMIT_FILE_VER
        );
        return Err(anyhow!("execution terminated unexpectedly"));
    }

    // read the barcode length
    rdr.read_exact(&mut rbuf)
        .context("couldn't read freq file buffer")?;
    let _bc_len = rbuf
        .pread::<u64>(0)
        .context("couldn't read freq file barcode length")?;

    // read the barcode -> frequency hashmap
    let freq_hm: HashMap<u64, u64> =
        bincode::deserialize_from(rdr).context("couldn't deserialize barcode to frequency map.")?;
    let total_to_collate = freq_hm.values().sum();
    let mut tsv_map = Vec::from_iter(freq_hm);

    // sort this so that we deal with largest cells (by # of reads) first
    // sort in _descending_ order by count.
    tsv_map.sort_unstable_by_key(|&a: &(u64, u64)| std::cmp::Reverse(a.1));

    collate_with_temp_memory(
        input_dir,
        rad_dir,
        num_threads,
        max_records,
        memory_budget_bytes,
        tsv_map,
        total_to_collate,
        compress_out,
        cmdline,
        version_str,
        log,
    )
}

fn get_orientation(mdata: &serde_json::Value) -> Result<Strand, StrandError> {
    // next line is ugly — should be a better way.  We need a char to
    // get the strand, so we get the correct field as a `str` then
    // use the chars iterator and get the first char.
    let ori_str: char = mdata["expected_ori"]
        .as_str()
        .unwrap()
        .chars()
        .next()
        .unwrap();
    Strand::from_char(&ori_str)
}

#[derive(Debug)]
enum FilterType {
    Filtered,
    Unfiltered,
}

fn get_filter_type(mdata: &serde_json::Value, log: &slog::Logger) -> FilterType {
    if let Some(fts) = mdata.get("permit-list-type") {
        match fts.as_str() {
            Some("unfiltered") => FilterType::Unfiltered,
            Some("filtered") => FilterType::Filtered,
            _ => FilterType::Filtered,
        }
    } else {
        info!(
            log,
            "permit-list-type key not present in JSON file; assuming list is filtered."
        );
        FilterType::Filtered
    }
}

fn get_most_ambiguous_record(mdata: &serde_json::Value, log: &slog::Logger) -> usize {
    if let Some(mar) = mdata.get("max-ambig-record") {
        match mar.as_u64() {
            Some(mv) => mv as usize,
            _ => 2500_usize,
        }
    } else {
        info!(
            log,
            "max-ambig-record key not present in JSON file; using default of 2,500. Please consider upgrading alevin-fry."
        );
        2500_usize
    }
}

/// Correct unmapped barcode counts for multi-barcode (Flex) data.
///
/// Reads the raw unmapped_bc_count.bin (self-describing format with per-field barcodes),
/// corrects both sample and cell barcodes, and writes unmapped_bc_count_collated.bin
/// with corrected per-field barcodes preserved for per-sample accuracy.
fn correct_unmapped_counts_multi_bc(
    unmapped_file: &std::path::Path,
    sample_permit_map: &HashMap<u64, u64>,
    sample_bc_to_idx: &HashMap<u64, usize>,
    samples: &[MultiBarcodeSampleCorrection],
    parent: &std::path::Path,
) {
    use libradicl::unmapped::{CollatedUnmappedCounts, UnmappedBcFormat, UnmappedBcRecordReader};

    // Multi-barcode: keyed by (corrected_sample, corrected_cell) for per-sample accuracy
    let mut collated = CollatedUnmappedCounts::new_multi(vec![
        libradicl::rad_types::RadIntId::U32,
        libradicl::rad_types::RadIntId::U32,
    ]);

    if let Ok(i_file) = File::open(unmapped_file) {
        let mut br = BufReader::new(i_file);

        // Try to read the self-describing header
        let format = match UnmappedBcFormat::read_header(&mut br) {
            Ok(Some(fmt)) => fmt,
            Ok(None) => {
                let s_path = parent.join("unmapped_bc_count_collated.bin");
                collated
                    .write_to_file(&s_path)
                    .expect("could not write collated file.");
                return;
            }
            Err(_) => {
                let s_path = parent.join("unmapped_bc_count_collated.bin");
                collated
                    .write_to_file(&s_path)
                    .expect("could not write collated file.");
                return;
            }
        };

        let mut reader = UnmappedBcRecordReader::new(format);

        while let Ok(Some((barcodes, count))) = reader.read_record(&mut br) {
            if barcodes.len() < 2 {
                // Single-barcode record — shouldn't happen in multi-BC path but handle gracefully
                let cell_bc = barcodes[0];
                collated.insert_single(cell_bc, count);
                continue;
            }

            // barcodes[0] = sample BC, barcodes[last] = cell BC
            let sample_bc = barcodes[0];
            let cell_bc = *barcodes.last().unwrap();

            // Correct sample BC
            let corrected_sample = match sample_permit_map.get(&sample_bc) {
                Some(&cs) => cs,
                None => continue,
            };
            let sample_idx = match sample_bc_to_idx.get(&corrected_sample) {
                Some(&idx) => idx,
                None => continue,
            };

            // Apply the same identity-first, unique-neighbor correction used
            // for mapped records during scatter.
            let Some(corrected_cell) = samples
                .get(sample_idx)
                .and_then(|sample| sample.correct_barcode(cell_bc))
            else {
                continue;
            };

            // Per-sample accuracy: key by (corrected_sample, corrected_cell)
            collated.insert_multi(corrected_sample, corrected_cell, count);
        }
    }

    let s_path = parent.join("unmapped_bc_count_collated.bin");
    collated
        .write_to_file(&s_path)
        .expect("could not write collated unmapped bc count.");
}

/// Correct unmapped barcode counts for single-barcode data.
///
/// Reads unmapped_bc_count.bin (self-describing format), applies cell BC
/// correction, writes unmapped_bc_count_collated.bin.
fn correct_unmapped_counts(
    correct_map: &HashMap<u64, u64>,
    unmapped_file: &std::path::Path,
    parent: &std::path::Path,
) {
    use libradicl::unmapped::{CollatedUnmappedCounts, UnmappedBcFormat, UnmappedBcRecordReader};

    let mut collated = CollatedUnmappedCounts::new_single(libradicl::rad_types::RadIntId::U32);

    if let Ok(i_file) = File::open(unmapped_file) {
        let mut br = BufReader::new(i_file);

        // Try new self-describing format first
        match UnmappedBcFormat::read_header(&mut br) {
            Ok(Some(fmt)) => {
                // New format: read structured records
                let mut reader = UnmappedBcRecordReader::new(fmt);
                while let Ok(Some((bcs, count))) = reader.read_record(&mut br) {
                    let raw_bc = bcs[0];
                    if let Some(&corrected) = correct_map.get(&raw_bc) {
                        collated.insert_single(corrected, count);
                    }
                }
            }
            Ok(None) => {
                // Empty file — nothing to do
            }
            Err(_) => {
                // Couldn't read header — try legacy format (raw u64+u32 pairs)
                // Reopen the file since we consumed the first byte
                if let Ok(i_file) = File::open(unmapped_file) {
                    let mut br = BufReader::new(i_file);
                    let mut rbuf = [0u8; std::mem::size_of::<u64>() + std::mem::size_of::<u32>()];
                    while br.read_exact(&mut rbuf[..]).is_ok() {
                        let k = rbuf.pread::<u64>(0).unwrap();
                        let v = rbuf.pread::<u32>(std::mem::size_of::<u64>()).unwrap();
                        if let Some(&ck) = correct_map.get(&k) {
                            collated.insert_single(ck, v);
                        }
                    }
                }
            }
        }
    }

    let s_path = parent.join("unmapped_bc_count_collated.bin");
    collated
        .write_to_file(&s_path)
        .expect("could not write collated unmapped bc count.");
}

#[allow(clippy::too_many_arguments)]
fn do_collate_single_barcode<P1, P2, A>(
    input_dir: P1,
    rad_dir: P2,
    rec_context: AlevinFryRecordContext,
    prelude: RadPrelude,
    mut br: BufReader<A>,
    end_header_pos: u64,
    num_threads: u32,
    max_records: u32,
    memory_budget_bytes: u64,
    tsv_map: Vec<(u64, u64)>,
    total_to_collate: u64,
    compress_out: bool,
    cmdline: &str,
    version: &str,
    log: &slog::Logger,
) -> anyhow::Result<()>
where
    P1: Into<PathBuf>,
    P2: AsRef<Path>,
    A: Read + Seek,
{
    let input_dir = input_dir.into();
    let parent = input_dir.as_path();
    let rad_parent = rad_dir.as_ref();
    let input_rad_path = rad_parent.join("map.rad");
    let expected_output_chunks = tsv_map.len() as u64;

    let meta_data_file = File::open(parent.join("generate_permit_list.json"))
        .context("could not open generate_permit_list.json")?;
    let metadata: serde_json::Value = serde_json::from_reader(&meta_data_file)?;
    let expected_orientation =
        get_orientation(&metadata).map_err(|error| anyhow!("could not read strand: {error}"))?;
    let velo_mode = metadata["velo_mode"]
        .as_bool()
        .context("could not read velo_mode from generate-permit-list metadata")?;

    let output_name = if velo_mode {
        "velo.map.collated.rad"
    } else if compress_out {
        "map.collated.rad.sz"
    } else {
        "map.collated.rad"
    };
    let output_path = parent.join(output_name);
    if output_path.exists() {
        std::fs::remove_file(&output_path)?;
    }
    let output = Arc::new(Mutex::new(BufWriter::with_capacity(
        1024 * 1024,
        File::create(&output_path)?,
    )));

    let correction_file =
        File::open(parent.join("permit_map.bin")).context("could not open permit_map.bin")?;
    let correction_map: HashMap<u64, u64> = bincode::deserialize_from(correction_file)
        .context("could not deserialize permit_map.bin")?;
    correct_unmapped_counts(
        &correction_map,
        &rad_parent.join("unmapped_bc_count.bin"),
        parent,
    );

    let header_end = br.get_mut().stream_position()? - br.buffer().len() as u64;
    let mut header = vec![0_u8; header_end as usize];
    File::open(&input_rad_path)?.read_exact(&mut header)?;
    let chunk_count_offset = (end_header_pos - std::mem::size_of::<u64>() as u64) as usize;
    header[chunk_count_offset..chunk_count_offset + 8]
        .copy_from_slice(&expected_output_chunks.to_le_bytes());
    if compress_out {
        let mut encoder = snap::write::FrameEncoder::new(Vec::with_capacity(header.len()));
        encoder.write_all(&header)?;
        header = encoder.into_inner()?;
    }
    output
        .lock()
        .map_err(|_| anyhow!("collated RAD output mutex was poisoned"))?
        .write_all(&header)?;

    let num_workers = (num_threads as usize).saturating_sub(1).max(1);
    let max_records_per_bucket = u64::from(max_records / num_workers as u32 + 1);
    let mut group_to_bucket = AHashMap::with_capacity(tsv_map.len());
    let mut bucket_records = Vec::<u64>::new();
    let mut current_records = 0_u64;
    for (barcode, frequency) in tsv_map {
        if current_records >= max_records_per_bucket {
            bucket_records.push(current_records);
            current_records = 0;
        }
        group_to_bucket.insert(barcode, bucket_records.len() as u32);
        current_records += frequency;
    }
    if current_records > 0 {
        bucket_records.push(current_records);
    }
    if bucket_records.is_empty() {
        return Err(anyhow!(
            "single-barcode permit list contains no output cells"
        ));
    }

    let correction_map = correction_map.into_iter().collect::<AHashMap<_, _>>();
    let plan = Arc::new(SingleBarcodeCollationPlan::new(
        correction_map,
        group_to_bucket,
        bucket_records.len(),
    )?);

    {
        let collate_metadata = json!({
            "cmd": cmdline,
            "version_str": version,
            "compressed_output": compress_out,
            "collation_mode": "optimized",
            "memory_budget_bytes": memory_budget_bytes,
        });
        let mut metadata_file = File::create(parent.join("collate.json"))?;
        serde_json::to_writer_pretty(&mut metadata_file, &collate_metadata)?;
    }

    let stats = collate_single_barcode(
        &mut br,
        prelude.hdr.num_chunks,
        rec_context,
        plan,
        output.clone(),
        parent,
        expected_orientation.into(),
        SingleBarcodeCollationOptions {
            num_threads: num_threads as usize,
            memory_budget_bytes,
            compress_output: compress_out,
        },
    )?;

    if stats.records_scattered != total_to_collate {
        return Err(anyhow!(
            "expected to collate {total_to_collate} records but retained {}",
            stats.records_scattered
        ));
    }
    if stats.output_chunks != expected_output_chunks {
        return Err(anyhow!(
            "expected {expected_output_chunks} output chunks but wrote {}",
            stats.output_chunks
        ));
    }
    output
        .lock()
        .map_err(|_| anyhow!("collated RAD output mutex was poisoned"))?
        .flush()?;
    info!(
        log,
        "Optimized single-barcode collation: {} records in {:.2}s scatter, {} chunks in {:.2}s gather ({} spool files, {} gather workers, {} KiB flush threshold)",
        stats.records_scattered.to_formatted_string(&Locale::en),
        stats.scatter_duration.as_secs_f64(),
        stats.output_chunks.to_formatted_string(&Locale::en),
        stats.gather_duration.as_secs_f64(),
        stats.num_scatter_workers,
        stats.num_gather_workers,
        stats.spool_flush_limit / 1024,
    );
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::manual_clamp)]
pub fn do_collate_with_temp<
    P1,
    P2,
    A: Read + std::io::Seek,
    B: ConvertiblePrimitiveInteger + std::convert::From<u64>,
    R: MappedRecord + KnownSize + CollatableMappedRecord<B>,
>(
    input_dir: P1,
    rad_dir: P2,
    rec_context: <R as MappedRecord>::ParsingContext,
    prelude: RadPrelude,
    mut br: BufReader<A>,
    end_header_pos: u64,
    num_threads: u32,
    max_records: u32,
    tsv_map: Vec<(u64, u64)>,
    total_to_collate: u64,
    compress_out: bool,
    cmdline: &str,
    version: &str,
    log: &slog::Logger,
) -> anyhow::Result<()>
where
    P1: Into<PathBuf>,
    P2: AsRef<Path>,
    u64: From<B>,
    // note; the 'static below simply means that the parsing context doesn't borrow anything so it
    // can be used in the closure.
    <R as MappedRecord>::ParsingContext:
        std::marker::Sync + Send + std::clone::Clone + 'static + std::fmt::Debug,
{
    let i_dir = std::path::Path::new(rad_dir.as_ref());
    let input_rad_path = i_dir.join("map.rad");

    // the number of corrected cells we'll write
    let expected_output_chunks = tsv_map.len() as u64;
    // the parent input directory
    let input_dir = input_dir.into();
    let parent = std::path::Path::new(input_dir.as_path());

    let n_workers = if num_threads > 1 {
        (num_threads - 1) as usize
    } else {
        1
    };

    // open the metadata file and read the json
    let meta_data_file = File::open(parent.join("generate_permit_list.json"))
        .context("could not open the generate_permit_list.json file.")?;
    let mdata: serde_json::Value = serde_json::from_reader(&meta_data_file)?;

    // velo_mode
    let velo_mode = mdata["velo_mode"]
        .as_bool()
        .context("couldn't read velo_mode from meta data")?;
    let expected_ori: Strand = match get_orientation(&mdata) {
        Ok(o) => o,
        Err(e) => {
            crit!(
                log,
                "Error reading strand info from {:#?} :: {}",
                &meta_data_file,
                e
            );
            return Err(anyhow!(e));
        }
    };

    let hdr = &prelude.hdr;
    info!(
        log,
        "paired : {:?}, ref_count : {}, num_chunks : {}, expected_ori : {:?}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en),
        expected_ori
    );

    let filter_type = get_filter_type(&mdata, log);
    let most_ambig_record = get_most_ambiguous_record(&mdata, log);

    // log the filter type
    info!(log, "filter_type = {:?}", filter_type);
    info!(
        log,
        "collated rad file {} be compressed",
        if compress_out { "will" } else { "will not" }
    );
    // because :
    // https://superuser.com/questions/865710/write-to-newfile-vs-overwriting-performance-issue
    let cfname = if velo_mode {
        "velo.map.collated.rad"
    } else if compress_out {
        "map.collated.rad.sz"
    } else {
        "map.collated.rad"
    };

    // writing the collate metadata
    {
        let collate_meta = json!({
            "cmd" : cmdline,
            "version_str" : version,
            "compressed_output" : compress_out,
        });

        let cm_path = parent.join("collate.json");
        let mut cm_file =
            std::fs::File::create(cm_path).context("could not create metadata file.")?;

        let cm_info_string =
            serde_json::to_string_pretty(&collate_meta).context("could not format json.")?;
        cm_file
            .write_all(cm_info_string.as_bytes())
            .context("cannot write to collate.json file")?;
    }

    let oname = parent.join(cfname);
    if oname.exists() {
        std::fs::remove_file(&oname)
            .with_context(|| format!("could not remove {}", oname.display()))?;
    }

    let ofile = File::create(parent.join(cfname))
        .with_context(|| format!("couldn't create directory {}", cfname))?;
    let owriter = Arc::new(Mutex::new(BufWriter::with_capacity(1048576, ofile)));

    // get the correction map
    let cmfile = std::fs::File::open(parent.join("permit_map.bin"))
        .context("couldn't open output permit_map.bin file")?;
    let correct_map: Arc<HashMap<u64, u64>> = Arc::new(bincode::deserialize_from(&cmfile).unwrap());

    // NOTE: the assumption of where the unmapped file will be
    // should be robustified
    let unmapped_file = i_dir.join("unmapped_bc_count.bin");
    correct_unmapped_counts(&correct_map, &unmapped_file, parent);

    info!(
        log,
        "deserialized correction map of length : {}",
        correct_map.len().to_formatted_string(&Locale::en)
    );

    // the exact position at the end of the header + file tags
    let pos = br.get_mut().stream_position().unwrap() - (br.buffer().len() as u64);

    // copy the header
    {
        // we want to copy up to the end of the header
        // minus the num chunks (sizeof u64), and then
        // write the actual number of chunks we expect.
        let chunk_bytes = std::mem::size_of::<u64>() as u64;
        let take_pos = end_header_pos - chunk_bytes;

        // This temporary file pointer and buffer will be dropped
        // at the end of this block (scope).
        let mut rfile = File::open(&input_rad_path).context("Couldn't open input RAD file")?;
        let mut hdr_buf = Cursor::new(vec![0u8; pos as usize]);

        rfile
            .read_exact(hdr_buf.get_mut())
            .context("couldn't read input file header")?;
        hdr_buf.set_position(take_pos);
        hdr_buf
            .write_all(&expected_output_chunks.to_le_bytes())
            .context("couldn't write num_chunks")?;
        hdr_buf.set_position(0);

        // compress the header buffer to a compressed buffer
        if compress_out {
            let mut compressed_buf =
                snap::write::FrameEncoder::new(Cursor::new(Vec::<u8>::with_capacity(pos as usize)));
            compressed_buf
                .write_all(hdr_buf.get_ref())
                .context("could not compress the output header.")?;
            hdr_buf = compressed_buf
                .into_inner()
                .context("couldn't unwrap the FrameEncoder.")?;
            hdr_buf.set_position(0);
        }

        if let Ok(mut oput) = owriter.lock() {
            oput.write_all(hdr_buf.get_ref())
                .context("could not write the output header.")?;
        }
    }

    // TODO: see if we can do this without the Arc
    let mut output_cache = Arc::new(HashMap::<u64, Arc<libradicl::TempBucket>>::new());

    // max_records is the max size of each intermediate file
    let mut total_allocated_records = 0;
    let mut allocated_records = 0;
    let mut temp_buckets = vec![(
        0,
        0,
        Arc::new(libradicl::TempBucket::from_id_and_parent(0, parent)),
    )];

    let max_records_per_thread = (max_records / n_workers as u32) + 1;
    // The tsv_map tells us, for each "true" barcode
    // how many records belong to it.  We can scan this information
    // to determine what true barcodes we will keep in memory.
    let mut num_bucket_chunks = 0u32;
    {
        let moutput_cache = Arc::make_mut(&mut output_cache);
        for rec in tsv_map.iter() {
            // corrected barcode points to the bucket
            // file.
            moutput_cache.insert(rec.0, temp_buckets.last().unwrap().2.clone());
            allocated_records += rec.1;
            num_bucket_chunks += 1;
            if allocated_records >= (max_records_per_thread as u64) {
                temp_buckets.last_mut().unwrap().0 = num_bucket_chunks;
                temp_buckets.last_mut().unwrap().1 = allocated_records as u32;
                let tn = temp_buckets.len() as u32;
                temp_buckets.push((
                    0,
                    0,
                    Arc::new(libradicl::TempBucket::from_id_and_parent(tn, parent)),
                ));
                total_allocated_records += allocated_records;
                allocated_records = 0;
                num_bucket_chunks = 0;
            }
        }
    }
    if num_bucket_chunks > 0 {
        temp_buckets.last_mut().unwrap().0 = num_bucket_chunks;
        temp_buckets.last_mut().unwrap().1 = allocated_records as u32;
    }
    total_allocated_records += allocated_records;
    info!(log, "Generated {} temporary buckets.", temp_buckets.len());

    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
        )
        .expect("ProgressStyle template was invalid")
        .progress_chars("╢▌▌░╟");

    let pbar_inner = ProgressBar::with_draw_target(
        Some(hdr.num_chunks),
        ProgressDrawTarget::stderr_with_hz(5u8), // update at most 5 times/sec.
    );

    pbar_inner.set_style(sty.clone());
    pbar_inner.tick();

    // create a thread-safe queue based on the number of worker threads
    let q = Arc::new(ArrayQueue::<(usize, Vec<u8>)>::new(4 * n_workers));

    // the number of cells left to process
    let chunks_to_process = Arc::new(AtomicUsize::new(hdr.num_chunks as usize));

    let mut thread_handles: Vec<thread::JoinHandle<u64>> = Vec::with_capacity(n_workers);

    let min_rec_len = 24usize; // smallest size an individual record can be loaded in memory
    let max_rec = max_records as usize;
    let num_buckets = temp_buckets.len();
    let num_threads = n_workers;
    let loc_buffer_size = (min_rec_len + (most_ambig_record * 4_usize) - 4_usize).max(
        (1000_usize.max((min_rec_len * max_rec) / (num_buckets * num_threads))).min(262_144_usize),
    ); //131072_usize);

    // for each worker, spawn off a thread
    for _worker in 0..n_workers {
        // each thread will need to access the work queue
        let in_q = q.clone();
        // the output cache and correction map
        let oc = output_cache.clone();
        let correct_map = correct_map.clone();
        // the number of chunks remaining to be processed
        let chunks_remaining = chunks_to_process.clone();
        // and knowledge of the UMI and BC types
        let nbuckets = temp_buckets.len();
        let loc_temp_buckets = temp_buckets.clone();
        //let owrite = owriter.clone();
        // now, make the worker thread
        let rec_context = rec_context.clone();
        let handle = std::thread::spawn(move || {
            // old code
            //let mut local_buffers = vec![Cursor::new(vec![0u8; loc_buffer_size]); nbuckets];

            // new approach (how much does this extra complexity matter?)
            // to avoid having a vector of cursors, where each cursor points to
            // a completely different vector (thus scattering memory of threads
            // and incurring the extra overhead for the capacity of the inner
            // vectors), we will have one backing chunk of memory.
            // NOTE: once stabilized, maybe using as_chunks_mut here
            // will be simpler (https://doc.rust-lang.org/std/primitive.slice.html#method.as_chunks_mut)

            // the memory that will back our temporary buffers
            let mut local_buffer_backing = vec![0u8; loc_buffer_size * nbuckets];
            // the vector of cursors we will use to write into our temporary buffers
            let mut local_buffers: Vec<Cursor<&mut [u8]>> = Vec::with_capacity(nbuckets);
            // The below is a bit tricky in rust but we basically break off each mutable slice
            // piece by piece.  Since `as_mut_slice(n)` returns the slices [0,n), [n,end) we
            // expect to chop off a first part of size `loc_buffer_size` a total of `nbuckets`
            // times.
            let mut tslice = local_buffer_backing.as_mut_slice();
            for _ in 0..nbuckets {
                let (first, rest) = tslice.split_at_mut(loc_buffer_size);
                //let brange = (bn*loc_buffer_size..(bn+1)*loc_buffer_size);
                local_buffers.push(Cursor::new(first));
                tslice = rest;
            }

            // pop from the work queue until everything is
            // processed
            while chunks_remaining.load(Ordering::SeqCst) > 0 {
                if let Some((_chunk_num, buf)) = in_q.pop() {
                    chunks_remaining.fetch_sub(1, Ordering::SeqCst);
                    let mut nbr = BufReader::new(&buf[..]);
                    libradicl::dump_corrected_cb_chunk_to_temp_file_generic::<B, _, R>(
                        &mut nbr,
                        &rec_context,
                        &correct_map,
                        &expected_ori,
                        &oc,
                        &mut local_buffers,
                        loc_buffer_size,
                    );
                }
            }

            // empty any remaining local buffers
            for (bucket_id, lb) in local_buffers.iter().enumerate() {
                let len = lb.position() as usize;
                if len > 0 {
                    let mut filebuf = loc_temp_buckets[bucket_id].2.bucket_writer.lock().unwrap();
                    filebuf.write_all(&lb.get_ref()[0..len]).unwrap();
                }
            }
            // return something more meaningful
            0
        });

        thread_handles.push(handle);
    } // for each worker

    // read each chunk
    pbar_inner.reset();
    let pb_msg = format!(
        "processing {} / {} total records",
        total_allocated_records, total_to_collate
    );
    pbar_inner.set_message(pb_msg);

    // read chunks from the input file and pass them to the
    // worker threads.
    let mut buf = vec![0u8; 65536];
    for cell_num in 0..(hdr.num_chunks as usize) {
        let (nbytes_chunk, nrec_chunk) = chunk::Chunk::<R>::read_header(&mut br);
        buf.resize(nbytes_chunk as usize, 0);
        buf.pwrite::<u32>(nbytes_chunk, 0)?;
        buf.pwrite::<u32>(nrec_chunk, 4)?;
        br.read_exact(&mut buf[8..]).unwrap();

        let mut bclone = (cell_num, buf.clone());
        // keep trying until we can push this payload
        while let Err(t) = q.push(bclone) {
            bclone = t;
            // no point trying to push if the queue is full
            while q.is_full() {}
        }
        pbar_inner.inc(1);
    }
    pbar_inner.finish();

    // wait for the worker threads to finish
    for h in thread_handles.drain(0..) {
        match h.join() {
            Ok(_) => {}
            Err(_e) => {
                info!(log, "thread panicked");
            }
        }
    }
    pbar_inner.finish_with_message("partitioned records into temporary files.");
    drop(q);

    // At this point, we are done with the "scatter"
    // phase of writing the records to the corresponding
    // intermediate files.  Now, we'll begin the gather
    // phase of collating the temporary files and merging
    // them into the final output file.

    for (i, temp_bucket) in temp_buckets.iter().enumerate() {
        // make sure we flush each temp bucket
        temp_bucket
            .2
            .bucket_writer
            .lock()
            .unwrap()
            .flush()
            .context("could not flush temporary output file!")?;
        // a sanity check that we have the correct number of records
        // and the expected number of bytes in each file
        let expected = temp_bucket.1;
        let observed = temp_bucket.2.num_records_written.load(Ordering::SeqCst);
        assert_eq!(expected, observed);

        let md = std::fs::metadata(parent.join(format!("bucket_{}.tmp", i)))?;
        let expected_bytes = temp_bucket.2.num_bytes_written.load(Ordering::SeqCst);
        let observed_bytes = md.len();
        assert_eq!(expected_bytes, observed_bytes);
    }

    //std::process::exit(1);

    // to hold the temp buckets threads will process
    let slack = (n_workers / 2).max(1_usize);
    let temp_bucket_queue_size = slack + n_workers;
    let fq = Arc::new(ArrayQueue::<(
        u32,
        u32,
        std::sync::Arc<libradicl::TempBucket>,
    )>::new(temp_bucket_queue_size));
    // the number of cells left to process
    let buckets_to_process = Arc::new(AtomicUsize::new(temp_buckets.len()));

    let pbar_gather = ProgressBar::new(temp_buckets.len() as u64);
    pbar_gather.set_style(sty);
    pbar_gather.tick();

    // for each worker, spawn off a thread
    for _worker in 0..n_workers {
        // each thread will need to access the work queue
        let in_q = fq.clone();
        // the output cache and correction map
        let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
        let mut cmap = HashMap::<u64, TempCellInfo, ahash::RandomState>::with_hasher(s);
        // alternative strategy
        // let mut cmap = HashMap::<u64, libradicl::CorrectedCbChunk, ahash::RandomState>::with_hasher(s);

        // the number of chunks remaining to be processed
        let buckets_remaining = buckets_to_process.clone();
        // have access to the input directory
        let input_dir: PathBuf = input_dir.clone();
        // the output file
        let owriter = owriter.clone();
        // and the progress bar
        let pbar_gather = pbar_gather.clone();
        let rec_context = rec_context.clone();
        // now, make the worker threads
        let handle = std::thread::spawn(move || {
            let ctx = rec_context;
            let mut local_chunks = 0u64;
            let parent = std::path::Path::new(&input_dir);
            // pop from the work queue until everything is
            // processed
            while buckets_remaining.load(Ordering::SeqCst) > 0 {
                if let Some(temp_bucket) = in_q.pop() {
                    buckets_remaining.fetch_sub(1, Ordering::SeqCst);
                    cmap.clear();

                    let fname = parent.join(format!("bucket_{}.tmp", temp_bucket.2.bucket_id));
                    // create a new handle for reading
                    let tfile = std::fs::File::open(&fname).expect("couldn't open temporary file.");
                    let mut treader = BufReader::new(tfile);

                    local_chunks += libradicl::collate_temporary_bucket_twopass_generic::<B, _, _, R>(
                        &mut treader,
                        &ctx,
                        temp_bucket.1,
                        &owriter,
                        compress_out,
                        &mut cmap,
                    ) as u64;

                    // we don't need the file or reader anymore
                    drop(treader);
                    std::fs::remove_file(fname).expect("could not delete temporary file.");

                    pbar_gather.inc(1);
                }
            }
            local_chunks
        });
        thread_handles.push(handle);
    } // for each worker

    // push the temporary buckets onto the work queue to be dispatched
    // by the worker threads.
    for temp_bucket in temp_buckets {
        let mut bclone = temp_bucket.clone();
        // keep trying until we can push this payload
        while let Err(t) = fq.push(bclone) {
            bclone = t;
            // no point trying to push if the queue is full
            while fq.is_full() {}
        }
        let expected = temp_bucket.1;
        let observed = temp_bucket.2.num_records_written.load(Ordering::SeqCst);
        assert_eq!(expected, observed);
    }

    // wait for all of the workers to finish
    let mut num_output_chunks = 0u64;
    for h in thread_handles.drain(0..) {
        match h.join() {
            Ok(c) => {
                num_output_chunks += c;
            }
            Err(_e) => {
                info!(log, "thread panicked");
            }
        }
    }
    pbar_gather.finish_with_message("gathered all temp files.");

    // make sure we wrote the same number of records that our
    // file suggested we should.
    assert_eq!(total_allocated_records, total_to_collate);

    info!(
        log,
        "writing num output chunks ({}) to header",
        num_output_chunks.to_formatted_string(&Locale::en)
    );

    info!(
        log,
        "expected number of output chunks {}",
        expected_output_chunks.to_formatted_string(&Locale::en)
    );

    assert_eq!(
        expected_output_chunks,
        num_output_chunks,
        "expected to write {} chunks but wrote {}",
        expected_output_chunks.to_formatted_string(&Locale::en),
        num_output_chunks.to_formatted_string(&Locale::en),
    );

    owriter.lock().unwrap().flush()?;
    info!(
        log,
        "finished collating input rad file {:?}.",
        i_dir.join("map.rad")
    );
    Ok(())
}

/// Historical record-count-based collation entry point.
#[allow(clippy::too_many_arguments)]
pub fn collate_with_temp<P1, P2>(
    input_dir: P1,
    rad_dir: P2,
    num_threads: u32,
    max_records: u32,
    tsv_map: Vec<(u64, u64)>,
    total_to_collate: u64,
    compress_out: bool,
    cmdline: &str,
    version: &str,
    log: &slog::Logger,
) -> anyhow::Result<()>
where
    P1: Into<PathBuf>,
    P2: AsRef<Path>,
{
    const BYTES_PER_LEGACY_RECORD: u64 = 72;
    const MIN_MEMORY_LIMIT: u64 = 256 * 1024 * 1024;
    let memory_budget_bytes = u64::from(max_records)
        .saturating_mul(BYTES_PER_LEGACY_RECORD)
        .max(MIN_MEMORY_LIMIT);
    collate_with_temp_memory(
        input_dir,
        rad_dir,
        num_threads,
        max_records,
        memory_budget_bytes,
        tsv_map,
        total_to_collate,
        compress_out,
        cmdline,
        version,
        log,
    )
}

#[allow(clippy::too_many_arguments, clippy::manual_clamp)]
fn collate_with_temp_memory<P1, P2>(
    input_dir: P1,
    rad_dir: P2,
    num_threads: u32,
    max_records: u32,
    memory_budget_bytes: u64,
    tsv_map: Vec<(u64, u64)>,
    total_to_collate: u64,
    compress_out: bool,
    cmdline: &str,
    version: &str,
    log: &slog::Logger,
) -> anyhow::Result<()>
where
    P1: Into<PathBuf>,
    P2: AsRef<Path>,
{
    let i_dir = std::path::Path::new(rad_dir.as_ref());

    if !i_dir.exists() {
        crit!(
            log,
            "the input RAD path {:?} does not exist",
            rad_dir.as_ref()
        );
        return Err(anyhow!("invalid input"));
    }

    let input_rad_path = i_dir.join("map.rad");
    let i_file = File::open(&input_rad_path).context("couldn't open input RAD file")?;
    let mut br = BufReader::new(i_file);

    let hdr = RadHeader::from_bytes(&mut br)?;

    // the exact position at the end of the header,
    // precisely sizeof(u64) bytes beyond the num_chunks field.
    let end_header_pos = br.get_ref().stream_position().unwrap() - (br.buffer().len() as u64);

    info!(
        log,
        "paired : {:?}, ref_count : {}, num_chunks : {}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en),
    );

    // file-level
    let fl_tags = rad_types::TagSection::from_bytes(&mut br)?;
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = rad_types::TagSection::from_bytes(&mut br)?;
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = rad_types::TagSection::from_bytes(&mut br)?;
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    // create the prelude and rebind the variables we need
    let prelude = RadPrelude::from_header_and_tag_sections(hdr, fl_tags, rl_tags, al_tags);

    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    info!(log, "File-level tag values {:?}", file_tag_map);

    let rec_type = afutils::get_record_type_from_prelude(&prelude, &file_tag_map);

    match rec_type {
        KnownRecordType::RnaLong(_bc_len) => {
            info!(log, "record type is long read single-cell RNA-seq");
            // long-read single cell
            info!(log, "long read single-cell");
            let parsing_context = prelude.get_record_context::<ScLongReadRecordContext>()?;
            do_collate_with_temp::<_, _, _, u64, ScLongReadRecordT<u64>>(
                input_dir,
                &rad_dir,
                parsing_context,
                prelude,
                br,
                end_header_pos,
                num_threads,
                max_records,
                tsv_map.clone(),
                total_to_collate,
                compress_out,
                cmdline,
                version,
                log,
            )
        }
        KnownRecordType::AtacSeq(_bc_len) => {
            info!(log, "record type is short read single-cell ATAC-seq");
            anyhow::bail!("To process atac-seq data, you should use the \"atac\" sub-command");
        }
        KnownRecordType::RnaShortPos(_bc_len) => {
            // alevin-fry with positions
            info!(log, "short read single-cell with position");
            let parsing_context = prelude.get_record_context::<AlevinFryRecordContext>()?;
            match parsing_context.bct {
                RadIntId::U64 | RadIntId::U32 | RadIntId::U16 | RadIntId::U8 => {
                    do_collate_with_temp::<_, _, _, u64, AlevinFryReadRecordWithPositionT<u64>>(
                        input_dir,
                        &rad_dir,
                        parsing_context,
                        prelude,
                        br,
                        end_header_pos,
                        num_threads,
                        max_records,
                        tsv_map.clone(),
                        total_to_collate,
                        compress_out,
                        cmdline,
                        version,
                        log,
                    )
                }
                RadIntId::U128 => {
                    unimplemented!()
                }
                _ => {
                    unimplemented!()
                }
            }
        }
        KnownRecordType::RnaShort(_bc_len) => {
            info!(log, "short read single-cell without poisition");
            let parsing_context = prelude.get_record_context::<AlevinFryRecordContext>()?;
            match parsing_context.bct {
                RadIntId::U64 | RadIntId::U32 | RadIntId::U16 | RadIntId::U8 => {
                    do_collate_single_barcode(
                        input_dir,
                        &rad_dir,
                        parsing_context,
                        prelude,
                        br,
                        end_header_pos,
                        num_threads,
                        max_records,
                        memory_budget_bytes,
                        tsv_map.clone(),
                        total_to_collate,
                        compress_out,
                        cmdline,
                        version,
                        log,
                    )
                }
                RadIntId::U128 => {
                    unimplemented!()
                }
                _ => {
                    unimplemented!()
                }
            }
        }
        KnownRecordType::RnaShortMultiBC(cell_bc_len, num_bc) => {
            info!(
                log,
                "record type is multi-barcode single-cell RNA-seq ({} barcode levels, cell BC len = {})",
                num_bc,
                cell_bc_len,
            );
            let parsing_context = prelude.get_record_context::<MultiBarcodeRecordContext>()?;
            info!(log, "Using the optimized libradicl collator");
            do_collate_multi_bc_fast(
                input_dir,
                &rad_dir,
                parsing_context,
                cell_bc_len as u32,
                prelude,
                br,
                end_header_pos,
                num_threads,
                memory_budget_bytes,
                total_to_collate,
                compress_out,
                cmdline,
                version,
                log,
            )
        }
    }
}

/// Multi-barcode hierarchical collation — optimized single-pass mode.
///
/// This function handles RAD files with multiple barcodes per read (e.g., 10x Flex).
/// It performs a single-pass scatter where records are bucketed by corrected
/// sample/cell identity, followed by a parallel gather.
///
/// The output is hierarchically grouped: all chunks for a sample are contiguous,
/// and all records for each corrected cell are collated into one chunk. Sample
/// ranges follow manifest order; cell order within a sample and record order
/// within a cell are unspecified.
#[allow(clippy::too_many_arguments)]
fn do_collate_multi_bc_fast<P1, P2, A: Read + Seek>(
    input_dir: P1,
    rad_dir: P2,
    rec_context: MultiBarcodeRecordContext,
    cell_bc_len: u32,
    prelude: RadPrelude,
    mut br: BufReader<A>,
    end_header_pos: u64,
    num_threads: u32,
    memory_budget_bytes: u64,
    total_to_collate: u64,
    compress_out: bool,
    cmdline: &str,
    version: &str,
    log: &slog::Logger,
) -> anyhow::Result<()>
where
    P1: Into<PathBuf>,
    P2: AsRef<Path>,
{
    let i_dir = std::path::Path::new(rad_dir.as_ref());
    let input_rad_path = i_dir.join("map.rad");
    let input_dir = input_dir.into();
    let parent = std::path::Path::new(input_dir.as_path());

    let n_workers = if num_threads > 1 {
        (num_threads - 1) as usize
    } else {
        1
    };

    // Compute the bit-shift and mask for building composite keys from
    // (sample_idx, corrected_cell_bc).  The cell barcode occupies the low
    // bits; the sample index is shifted above it.
    let cell_bc_bits = (cell_bc_len * 2) as u64; // 2 bits per nucleotide
    let cell_bc_mask = (1u64 << cell_bc_bits) - 1;
    let make_composite_key = |sample_idx: u64, cell_bc: u64| -> u64 {
        (sample_idx << cell_bc_bits) | (cell_bc & cell_bc_mask)
    };

    // Read metadata
    let meta_data_file = File::open(parent.join("generate_permit_list.json"))
        .context("could not open the generate_permit_list.json file.")?;
    let mdata: serde_json::Value = serde_json::from_reader(&meta_data_file)?;
    let expected_ori: Strand =
        get_orientation(&mdata).map_err(|e| anyhow!("Error reading strand info: {}", e))?;

    // Load sample_info.json to get sample metadata
    let sample_info_file = File::open(parent.join("sample_info.json")).context(
        "could not open sample_info.json — was generate-permit-list run with --sample-bc-list?",
    )?;
    let sample_info: serde_json::Value = serde_json::from_reader(&sample_info_file)?;

    let num_samples = sample_info["num_samples"]
        .as_u64()
        .context("couldn't read num_samples from sample_info.json")? as usize;
    let sample_entries = sample_info["samples"]
        .as_array()
        .context("couldn't read samples array from sample_info.json")?;

    info!(log, "Loading permit maps for {} samples", num_samples);

    // Verify the composite key (sample_idx << cell_bc_bits) | cell_bc fits
    // in u64.  We need ceil(log2(num_samples)) + cell_bc_bits <= 64.
    {
        let sample_id_bits = if num_samples <= 1 {
            0u64
        } else {
            64 - (num_samples as u64 - 1).leading_zeros() as u64
        };
        if sample_id_bits + cell_bc_bits > 64 {
            return Err(anyhow!(
                "Cannot collate: {} samples requires {} bits for the sample index, \
                 plus {} bits for {}bp cell barcodes = {} bits total, which exceeds \
                 the 64-bit composite key capacity.",
                num_samples,
                sample_id_bits,
                cell_bc_bits,
                cell_bc_len,
                sample_id_bits + cell_bc_bits
            ));
        }
    }

    // Load sample_permit_map.bin
    let sample_map_file = File::open(parent.join("sample_permit_map.bin"))
        .context("couldn't open sample_permit_map.bin")?;
    let sample_permit_map: HashMap<u64, u64> =
        bincode::deserialize_from(BufReader::new(sample_map_file))
            .map_err(|e| anyhow!("couldn't deserialize sample_permit_map.bin: {}", e))?;

    // Build sample_bc -> sample_idx mapping
    let mut sample_bc_to_idx: HashMap<u64, usize> = HashMap::new();
    let mut sample_names: Vec<String> = Vec::new();
    for (idx, entry) in sample_entries.iter().enumerate() {
        let bc_str = entry["barcode"].as_str().unwrap_or("0x0");
        let bc = u64::from_str_radix(bc_str.trim_start_matches("0x"), 16).unwrap_or(0);
        sample_bc_to_idx.insert(bc, idx);
        sample_names.push(
            entry["name"]
                .as_str()
                .unwrap_or(&format!("{:x}", bc))
                .to_string(),
        );
    }

    // Load per-sample valid barcodes from permit_freq.bin (NOT the huge permit_map.bin).
    // The valid barcodes are the identity-correction tier (self-correcting).
    // For non-identity correction, we build a BarcodeLookupMap per sample from
    // the valid barcodes and do on-the-fly 1-edit neighbor lookup during scatter.
    let mut per_sample_valid_bcs: Vec<Vec<u64>> = Vec::new();
    let mut all_cell_freqs: Vec<(u64, u64, usize)> = Vec::new(); // (corrected_cell_bc, freq, sample_idx)

    // cell_bc_len is passed in from the caller (extracted from the RAD file tags)

    for (sample_idx, _entry) in sample_entries.iter().enumerate() {
        let sample_name = &sample_names[sample_idx];
        let sample_dir = parent.join(format!("sample_{}", sample_name));

        let freq_path = sample_dir.join("permit_freq.bin");
        if freq_path.exists() {
            let freq_file = File::open(&freq_path)?;
            let mut freq_reader = BufReader::new(freq_file);
            // Read header: version (u64) + bc_len (u64)
            let mut hdr_buf = [0u8; 16];
            freq_reader.read_exact(&mut hdr_buf)?;
            let _bc_len_from_file = hdr_buf.pread::<u64>(8).unwrap_or(16) as u32;
            let freq_map: HashMap<u64, u64> = bincode::deserialize_from(freq_reader)
                .map_err(|e| anyhow!("couldn't deserialize {}: {}", freq_path.display(), e))?;

            let valid_bcs: Vec<u64> = freq_map.keys().copied().collect();
            if valid_bcs.is_empty() {
                debug!(log, "Sample '{}': no retained cell barcodes", sample_name);
            } else {
                info!(
                    log,
                    "Sample '{}': {} valid cell barcodes from permit_freq.bin",
                    sample_name,
                    valid_bcs.len().to_formatted_string(&Locale::en),
                );
            }

            for (bc, freq) in &freq_map {
                all_cell_freqs.push((*bc, *freq, sample_idx));
            }

            per_sample_valid_bcs.push(valid_bcs);
        } else {
            debug!(log, "Sample '{}': no permit freq (no reads)", sample_name);
            per_sample_valid_bcs.push(Vec::new());
        }
    }

    // Sort cell frequencies: group by sample_idx first, then by frequency descending
    all_cell_freqs.sort_by(|a, b| a.2.cmp(&b.2).then(b.1.cmp(&a.1)));

    // Build a sparse-plate-sidx -> manifest-ordinal map.  `sample_idx` so far
    // is a position in `sample_entries` (the full chemistry plate, e.g. 384
    // wells for 10x Flex v2), so it is sparse: only the wells actually used
    // in this run end up in `all_cell_freqs`.  The collation manifest and
    // quant's `sample_names` Vec are densely packed (one entry per present
    // sample, sorted by sidx), and quant.rs indexes that Vec directly with
    // the value stored in barcodes[0].  We must therefore write the
    // *manifest ordinal*, not the sparse plate sidx, into barcodes[0].
    // See: https://github.com/COMBINE-lab/simpleaf/issues/195
    //
    // Implementation: direct-addressed Vec<usize> indexed by plate sidx.
    // The key space is dense, small, and known (0..num_samples), so a Vec
    // is strictly better than a HashMap (~1.5 KB, single L1-resident load,
    // no hash per record across 10^8-10^9 scatter records).  Sentinel
    // `usize::MAX` flags unused plate positions; the scatter loop only ever
    // looks up sidxs that appeared in `all_cell_freqs`, so the sentinel
    // exists for defensive debugging only.
    let mut sidx_to_ord: Vec<usize> = vec![usize::MAX; num_samples];
    {
        let mut present: Vec<usize> = all_cell_freqs.iter().map(|(_, _, s)| *s).collect();
        present.sort_unstable();
        present.dedup();
        for (ord, sidx) in present.into_iter().enumerate() {
            sidx_to_ord[sidx] = ord;
        }
    }
    let sidx_to_ord = Arc::new(sidx_to_ord);

    // Aggregate the manifest metadata while `all_cell_freqs` is resident. The
    // same sorted vector is consumed below when constructing the bucket map,
    // avoiding a second multi-million-entry `tsv_map` allocation.
    let mut sample_chunk_counts: HashMap<usize, (u64, u64)> = HashMap::new();
    for &(_, freq, sample_index) in &all_cell_freqs {
        let entry = sample_chunk_counts.entry(sample_index).or_insert((0, 0));
        entry.0 += 1;
        entry.1 += freq;
    }

    let expected_output_chunks = all_cell_freqs.len() as u64;
    info!(
        log,
        "Total cells across all samples: {}, total records to collate: {}",
        expected_output_chunks.to_formatted_string(&Locale::en),
        total_to_collate.to_formatted_string(&Locale::en),
    );

    // Create output file
    let cfname = if compress_out {
        "map.collated.rad.sz"
    } else {
        "map.collated.rad"
    };

    // Write collate metadata
    {
        let collate_meta = json!({
            "cmd": cmdline,
            "version_str": version,
            "compressed_output": compress_out,
            "multi_barcode": true,
            "num_samples": num_samples,
            "collation_mode": "optimized",
            "memory_budget_bytes": memory_budget_bytes,
        });
        let cm_path = parent.join("collate.json");
        let mut cm_file = File::create(cm_path)?;
        let cm_str = serde_json::to_string_pretty(&collate_meta)?;
        cm_file.write_all(cm_str.as_bytes())?;
    }

    let oname = parent.join(cfname);
    if oname.exists() {
        std::fs::remove_file(&oname)?;
    }
    let ofile = File::create(&oname)?;
    let owriter = Arc::new(Mutex::new(BufWriter::with_capacity(1048576, ofile)));

    // Copy header with updated num_chunks
    let pos = br.get_mut().stream_position().unwrap() - (br.buffer().len() as u64);
    {
        let chunk_bytes = std::mem::size_of::<u64>() as u64;
        let take_pos = end_header_pos - chunk_bytes;
        let mut rfile = File::open(&input_rad_path)?;
        let mut hdr_buf = Cursor::new(vec![0u8; pos as usize]);
        rfile.read_exact(hdr_buf.get_mut())?;
        hdr_buf.set_position(take_pos);
        hdr_buf.write_all(&expected_output_chunks.to_le_bytes())?;
        hdr_buf.set_position(0);

        if compress_out {
            let mut compressed_buf =
                snap::write::FrameEncoder::new(Cursor::new(Vec::<u8>::with_capacity(pos as usize)));
            compressed_buf.write_all(hdr_buf.get_ref())?;
            hdr_buf = compressed_buf.into_inner()?;
            hdr_buf.set_position(0);
        }

        if let Ok(mut oput) = owriter.lock() {
            oput.write_all(hdr_buf.get_ref())?;
        }
    }

    // Partition corrected sample/cell groups into logical gather buckets.
    // Physical temporary storage is owned by libradicl and is bounded by the
    // number of scatter workers rather than the number of logical buckets.
    let num_cell_groups = all_cell_freqs.len();
    let mut group_to_bucket = AHashMap::<u64, u32>::with_capacity(num_cell_groups);
    let mut allocated_records: u64 = 0;
    let mut bucket_summaries = Vec::new(); // (groups, expected records)
    // Keep logical bucket sizing independent of the byte budget. Smaller
    // budgets reduce spool buffering and gather concurrency in libradicl;
    // creating more logical buckets would instead increase extent metadata
    // and can perversely raise RSS.
    let max_records_per_thread = (30_000_000 / n_workers as u32) + 1;
    let mut num_bucket_chunks = 0u32;
    let mut current_sample = None;
    for (cell_barcode, frequency, sample_index) in all_cell_freqs {
        if num_bucket_chunks > 0
            && (current_sample != Some(sample_index)
                || allocated_records >= max_records_per_thread as u64)
        {
            bucket_summaries.push((num_bucket_chunks, allocated_records));
            allocated_records = 0;
            num_bucket_chunks = 0;
        }
        current_sample = Some(sample_index);
        let bucket_id = bucket_summaries.len() as u32;
        let composite_key = make_composite_key(sample_index as u64, cell_barcode);
        group_to_bucket.insert(composite_key, bucket_id);
        allocated_records += frequency;
        num_bucket_chunks += 1;
    }
    if num_bucket_chunks > 0 {
        bucket_summaries.push((num_bucket_chunks, allocated_records));
    }
    info!(
        log,
        "Generated {} logical buckets for multi-barcode collation.",
        bucket_summaries.len()
    );

    let observed_sample_to_index: AHashMap<u64, usize> = sample_permit_map
        .iter()
        .filter_map(|(&observed, corrected)| {
            sample_bc_to_idx
                .get(corrected)
                .map(|&sample_index| (observed, sample_index))
        })
        .collect();
    let total_valid_cells: usize = per_sample_valid_bcs.iter().map(Vec::len).sum();
    let mut sample_corrections = Vec::with_capacity(per_sample_valid_bcs.len());
    for (sample_index, valid_barcodes) in per_sample_valid_bcs.into_iter().enumerate() {
        let ordinal = sidx_to_ord[sample_index];
        let correction = MultiBarcodeSampleCorrection::new(
            if ordinal == usize::MAX {
                0
            } else {
                ordinal as u64
            },
            valid_barcodes,
            cell_bc_len,
        );
        sample_corrections.push(correction);
    }
    info!(
        log,
        "Cell correction: {} valid barcodes across {} sample positions (bc_len={})",
        total_valid_cells.to_formatted_string(&Locale::en),
        num_samples,
        cell_bc_len,
    );

    let plan = Arc::new(MultiBarcodeCollationPlan::new(
        observed_sample_to_index,
        sample_corrections,
        group_to_bucket,
        cell_bc_bits as u32,
        bucket_summaries.len(),
    )?);
    info!(
        log,
        "Starting libradicl multi-barcode collation with {} workers and a {:.2} GiB memory budget",
        n_workers,
        memory_budget_bytes as f64 / (1024.0 * 1024.0 * 1024.0),
    );
    let engine_stats = collate_multi_barcode(
        &mut br,
        prelude.hdr.num_chunks,
        rec_context.clone(),
        plan.clone(),
        owriter.clone(),
        parent,
        expected_ori.into(),
        MultiBarcodeCollationOptions {
            num_threads: num_threads as usize,
            memory_budget_bytes,
            compress_output: compress_out,
        },
    )?;
    let total_output_chunks = engine_stats.output_chunks;
    info!(
        log,
        "Scatter phase complete: {} records in {:.2}s; gather produced {} chunks in {:.2}s ({} spool files, {} gather workers, {} KiB flush threshold)",
        engine_stats
            .records_scattered
            .to_formatted_string(&Locale::en),
        engine_stats.scatter_duration.as_secs_f64(),
        total_output_chunks.to_formatted_string(&Locale::en),
        engine_stats.gather_duration.as_secs_f64(),
        engine_stats.num_scatter_workers,
        engine_stats.num_gather_workers,
        engine_stats.spool_flush_limit / 1024,
    );

    // Correct unmapped barcode counts for multi-barcode data
    let unmapped_file = i_dir.join("unmapped_bc_count.bin");
    correct_unmapped_counts_multi_bc(
        &unmapped_file,
        &sample_permit_map,
        &sample_bc_to_idx,
        plan.samples(),
        parent,
    );

    // Build the manifest from the aggregate counts computed before collation.
    let mut manifest = CollationManifest::new(vec!["sample".to_string(), "cell".to_string()]);

    let mut chunk_offset: u64 = 0;
    let mut sample_indices: Vec<usize> = sample_chunk_counts.keys().copied().collect();
    sample_indices.sort();
    for sidx in sample_indices {
        let (num_cells, num_records) = sample_chunk_counts[&sidx];
        if sidx < sample_entries.len() {
            let bc_str = &sample_entries[sidx]["barcode"];
            let bc =
                u64::from_str_radix(bc_str.as_str().unwrap_or("0").trim_start_matches("0x"), 16)
                    .unwrap_or(0);
            manifest.add_sample_group(SampleGroup {
                key: bc,
                name: Some(sample_names[sidx].clone()),
                chunk_start: chunk_offset,
                num_chunks: num_cells,
                num_records,
            });
        }
        chunk_offset += num_cells;
    }

    // Write collation manifest
    let manifest_path = parent.join("collation_manifest.bin");
    manifest.write_to_file(&manifest_path)?;
    info!(
        log,
        "Wrote collation manifest: {} samples, {} total chunks",
        manifest.sample_groups.len(),
        total_output_chunks,
    );

    // Flush output
    if let Ok(mut oput) = owriter.lock() {
        oput.flush()?;
    }
    drop(owriter);

    // Backpatch num_chunks in the output file header
    if !compress_out {
        let chunk_bytes = std::mem::size_of::<u64>() as u64;
        let nc_pos = end_header_pos - chunk_bytes;
        let mut ofile = std::fs::OpenOptions::new().write(true).open(&oname)?;
        ofile.seek(std::io::SeekFrom::Start(nc_pos))?;
        ofile.write_all(&total_output_chunks.to_le_bytes())?;
        info!(
            log,
            "Backpatched num_chunks to {} in output file", total_output_chunks
        );
    }

    info!(
        log,
        "Multi-barcode collation complete: {} output chunks",
        total_output_chunks.to_formatted_string(&Locale::en),
    );

    Ok(())
}
/*
#[allow(clippy::too_many_arguments, clippy::manual_clamp)]
pub fn collate_with_temp<P1, P2>(
    input_dir: P1,
    rad_dir: P2,
    num_threads: u32,
    max_records: u32,
    tsv_map: Vec<(u64, u64)>,
    total_to_collate: u64,
    compress_out: bool,
    cmdline: &str,
    version: &str,
    log: &slog::Logger,
) -> anyhow::Result<()>
where
    P1: Into<PathBuf>,
    P2: AsRef<Path>,
{
    // the number of corrected cells we'll write
    let expected_output_chunks = tsv_map.len() as u64;
    // the parent input directory
    let input_dir = input_dir.into();
    let parent = std::path::Path::new(input_dir.as_path());

    let n_workers = if num_threads > 1 {
        (num_threads - 1) as usize
    } else {
        1
    };

    // open the metadata file and read the json
    let meta_data_file = File::open(parent.join("generate_permit_list.json"))
        .context("could not open the generate_permit_list.json file.")?;
    let mdata: serde_json::Value = serde_json::from_reader(&meta_data_file)?;

    // velo_mode
    let velo_mode = mdata["velo_mode"]
        .as_bool()
        .context("couldn't read velo_mode from meta data")?;
    let expected_ori: Strand = match get_orientation(&mdata) {
        Ok(o) => o,
        Err(e) => {
            crit!(
                log,
                "Error reading strand info from {:#?} :: {}",
                &meta_data_file,
                e
            );
            return Err(anyhow!(e));
        }
    };

    let filter_type = get_filter_type(&mdata, log);
    let most_ambig_record = get_most_ambiguous_record(&mdata, log);

    // log the filter type
    info!(log, "filter_type = {:?}", filter_type);
    info!(
        log,
        "collated rad file {} be compressed",
        if compress_out { "will" } else { "will not" }
    );
    // because :
    // https://superuser.com/questions/865710/write-to-newfile-vs-overwriting-performance-issue
    let cfname = if velo_mode {
        "velo.map.collated.rad"
    } else if compress_out {
        "map.collated.rad.sz"
    } else {
        "map.collated.rad"
    };

    // writing the collate metadata
    {
        let collate_meta = json!({
            "cmd" : cmdline,
            "version_str" : version,
            "compressed_output" : compress_out,
        });

        let cm_path = parent.join("collate.json");
        let mut cm_file =
            std::fs::File::create(cm_path).context("could not create metadata file.")?;

        let cm_info_string =
            serde_json::to_string_pretty(&collate_meta).context("could not format json.")?;
        cm_file
            .write_all(cm_info_string.as_bytes())
            .context("cannot write to collate.json file")?;
    }

    let oname = parent.join(cfname);
    if oname.exists() {
        std::fs::remove_file(&oname)
            .with_context(|| format!("could not remove {}", oname.display()))?;
    }

    let ofile = File::create(parent.join(cfname))
        .with_context(|| format!("couldn't create directory {}", cfname))?;
    let owriter = Arc::new(Mutex::new(BufWriter::with_capacity(1048576, ofile)));

    let i_dir = std::path::Path::new(rad_dir.as_ref());

    if !i_dir.exists() {
        crit!(
            log,
            "the input RAD path {:?} does not exist",
            rad_dir.as_ref()
        );
        return Err(anyhow!("invalid input"));
    }

    let input_rad_path = i_dir.join("map.rad");
    let i_file = File::open(&input_rad_path).context("couldn't open input RAD file")?;
    let mut br = BufReader::new(i_file);

    let hdr = RadHeader::from_bytes(&mut br)?;

    // the exact position at the end of the header,
    // precisely sizeof(u64) bytes beyond the num_chunks field.
    let end_header_pos = br.get_ref().stream_position().unwrap() - (br.buffer().len() as u64);

    info!(
        log,
        "paired : {:?}, ref_count : {}, num_chunks : {}, expected_ori : {:?}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en),
        expected_ori
    );

    // file-level
    let fl_tags = rad_types::TagSection::from_bytes(&mut br)?;
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = rad_types::TagSection::from_bytes(&mut br)?;
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = rad_types::TagSection::from_bytes(&mut br)?;
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    // create the prelude and rebind the variables we need
    let prelude = RadPrelude::from_header_and_tag_sections(hdr, fl_tags, rl_tags, al_tags);
    let hdr = &prelude.hdr;
    let rl_tags = &prelude.read_tags;

    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br);
    info!(log, "File-level tag values {:?}", file_tag_map);

    // the exact position at the end of the header + file tags
    let pos = br.get_ref().stream_position().unwrap() - (br.buffer().len() as u64);

    // copy the header
    {
        // we want to copy up to the end of the header
        // minus the num chunks (sizeof u64), and then
        // write the actual number of chunks we expect.
        let chunk_bytes = std::mem::size_of::<u64>() as u64;
        let take_pos = end_header_pos - chunk_bytes;

        // This temporary file pointer and buffer will be dropped
        // at the end of this block (scope).
        let mut rfile = File::open(&input_rad_path).context("Couldn't open input RAD file")?;
        let mut hdr_buf = Cursor::new(vec![0u8; pos as usize]);

        rfile
            .read_exact(hdr_buf.get_mut())
            .context("couldn't read input file header")?;
        hdr_buf.set_position(take_pos);
        hdr_buf
            .write_all(&expected_output_chunks.to_le_bytes())
            .context("couldn't write num_chunks")?;
        hdr_buf.set_position(0);

        // compress the header buffer to a compressed buffer
        if compress_out {
            let mut compressed_buf =
                snap::write::FrameEncoder::new(Cursor::new(Vec::<u8>::with_capacity(pos as usize)));
            compressed_buf
                .write_all(hdr_buf.get_ref())
                .context("could not compress the output header.")?;
            hdr_buf = compressed_buf
                .into_inner()
                .context("couldn't unwrap the FrameEncoder.")?;
            hdr_buf.set_position(0);
        }

        if let Ok(mut oput) = owriter.lock() {
            oput.write_all(hdr_buf.get_ref())
                .context("could not write the output header.")?;
        }
    }

    // get the correction map
    let cmfile = std::fs::File::open(parent.join("permit_map.bin"))
        .context("couldn't open output permit_map.bin file")?;
    let correct_map: Arc<HashMap<u64, u64>> = Arc::new(bincode::deserialize_from(&cmfile).unwrap());

    // NOTE: the assumption of where the unmapped file will be
    // should be robustified
    let unmapped_file = i_dir.join("unmapped_bc_count.bin");
    correct_unmapped_counts(&correct_map, &unmapped_file, parent);

    info!(
        log,
        "deserialized correction map of length : {}",
        correct_map.len().to_formatted_string(&Locale::en)
    );

    let bct = rl_tags.tags[0].typeid;
    let umit = rl_tags.tags[1].typeid;

    let cc = chunk::AlevinFryChunkContext {
        num_chunks: hdr.num_chunks,
        bc_type: rad_types::encode_type_tag(bct).expect("valid barcode tag type"),
        umi_type: rad_types::encode_type_tag(umit).expect("valid umi tag type"),
    };

    // TODO: see if we can do this without the Arc
    let mut output_cache = Arc::new(HashMap::<u64, Arc<libradicl::TempBucket>>::new());

    // max_records is the max size of each intermediate file
    let mut total_allocated_records = 0;
    let mut allocated_records = 0;
    let mut temp_buckets = vec![(
        0,
        0,
        Arc::new(libradicl::TempBucket::from_id_and_parent(0, parent)),
    )];

    let max_records_per_thread = (max_records / n_workers as u32) + 1;
    // The tsv_map tells us, for each "true" barcode
    // how many records belong to it.  We can scan this information
    // to determine what true barcodes we will keep in memory.
    let mut num_bucket_chunks = 0u32;
    {
        let moutput_cache = Arc::make_mut(&mut output_cache);
        for rec in tsv_map.iter() {
            // corrected barcode points to the bucket
            // file.
            moutput_cache.insert(rec.0, temp_buckets.last().unwrap().2.clone());
            allocated_records += rec.1;
            num_bucket_chunks += 1;
            if allocated_records >= (max_records_per_thread as u64) {
                temp_buckets.last_mut().unwrap().0 = num_bucket_chunks;
                temp_buckets.last_mut().unwrap().1 = allocated_records as u32;
                let tn = temp_buckets.len() as u32;
                temp_buckets.push((
                    0,
                    0,
                    Arc::new(libradicl::TempBucket::from_id_and_parent(tn, parent)),
                ));
                total_allocated_records += allocated_records;
                allocated_records = 0;
                num_bucket_chunks = 0;
            }
        }
    }
    if num_bucket_chunks > 0 {
        temp_buckets.last_mut().unwrap().0 = num_bucket_chunks;
        temp_buckets.last_mut().unwrap().1 = allocated_records as u32;
    }
    total_allocated_records += allocated_records;
    info!(log, "Generated {} temporary buckets.", temp_buckets.len());

    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
        )
        .expect("ProgressStyle template was invalid")
        .progress_chars("╢▌▌░╟");

    let pbar_inner = ProgressBar::with_draw_target(
        Some(cc.num_chunks),
        ProgressDrawTarget::stderr_with_hz(5u8), // update at most 5 times/sec.
    );

    pbar_inner.set_style(sty.clone());
    pbar_inner.tick();

    // create a thread-safe queue based on the number of worker threads
    let q = Arc::new(ArrayQueue::<(usize, Vec<u8>)>::new(4 * n_workers));

    // the number of cells left to process
    let chunks_to_process = Arc::new(AtomicUsize::new(cc.num_chunks as usize));

    let mut thread_handles: Vec<thread::JoinHandle<u64>> = Vec::with_capacity(n_workers);

    let min_rec_len = 24usize; // smallest size an individual record can be loaded in memory
    let max_rec = max_records as usize;
    let num_buckets = temp_buckets.len();
    let num_threads = n_workers;
    let loc_buffer_size = (min_rec_len + (most_ambig_record * 4_usize) - 4_usize).max(
        (1000_usize.max((min_rec_len * max_rec) / (num_buckets * num_threads))).min(262_144_usize),
    ); //131072_usize);

    // for each worker, spawn off a thread
    for _worker in 0..n_workers {
        // each thread will need to access the work queue
        let in_q = q.clone();
        // the output cache and correction map
        let oc = output_cache.clone();
        let correct_map = correct_map.clone();
        // the number of chunks remaining to be processed
        let chunks_remaining = chunks_to_process.clone();
        // and knowledge of the UMI and BC types
        let bc_type = rad_types::decode_int_type_tag(cc.bc_type).expect("unknown barcode type id.");
        let umi_type =
            rad_types::decode_int_type_tag(cc.umi_type).expect("unknown barcode type id.");
        let nbuckets = temp_buckets.len();
        let loc_temp_buckets = temp_buckets.clone();
        //let owrite = owriter.clone();
        // now, make the worker thread
        let handle = std::thread::spawn(move || {
            // old code
            //let mut local_buffers = vec![Cursor::new(vec![0u8; loc_buffer_size]); nbuckets];

            // new approach (how much does this extra complexity matter?)
            // to avoid having a vector of cursors, where each cursor points to
            // a completely different vector (thus scattering memory of threads
            // and incurring the extra overhead for the capacity of the inner
            // vectors), we will have one backing chunk of memory.
            // NOTE: once stabilized, maybe using as_chunks_mut here
            // will be simpler (https://doc.rust-lang.org/std/primitive.slice.html#method.as_chunks_mut)

            // the memory that will back our temporary buffers
            let mut local_buffer_backing = vec![0u8; loc_buffer_size * nbuckets];
            // the vector of cursors we will use to write into our temporary buffers
            let mut local_buffers: Vec<Cursor<&mut [u8]>> = Vec::with_capacity(nbuckets);
            // The below is a bit tricky in rust but we basically break off each mutable slice
            // piece by piece.  Since `as_mut_slice(n)` returns the slices [0,n), [n,end) we
            // expect to chop off a first part of size `loc_buffer_size` a total of `nbuckets`
            // times.
            let mut tslice = local_buffer_backing.as_mut_slice();
            for _ in 0..nbuckets {
                let (first, rest) = tslice.split_at_mut(loc_buffer_size);
                //let brange = (bn*loc_buffer_size..(bn+1)*loc_buffer_size);
                local_buffers.push(Cursor::new(first));
                tslice = rest;
            }

            // pop from the work queue until everything is
            // processed
            while chunks_remaining.load(Ordering::SeqCst) > 0 {
                if let Some((_chunk_num, buf)) = in_q.pop() {
                    chunks_remaining.fetch_sub(1, Ordering::SeqCst);
                    let mut nbr = BufReader::new(&buf[..]);
                    libradicl::dump_corrected_cb_chunk_to_temp_file(
                        &mut nbr,
                        &bc_type,
                        &umi_type,
                        &correct_map,
                        &expected_ori,
                        &oc,
                        &mut local_buffers,
                        loc_buffer_size,
                    );
                }
            }

            // empty any remaining local buffers
            for (bucket_id, lb) in local_buffers.iter().enumerate() {
                let len = lb.position() as usize;
                if len > 0 {
                    let mut filebuf = loc_temp_buckets[bucket_id].2.bucket_writer.lock().unwrap();
                    filebuf.write_all(&lb.get_ref()[0..len]).unwrap();
                }
            }
            // return something more meaningful
            0
        });

        thread_handles.push(handle);
    } // for each worker

    // read each chunk
    pbar_inner.reset();
    let pb_msg = format!(
        "processing {} / {} total records",
        total_allocated_records, total_to_collate
    );
    pbar_inner.set_message(pb_msg);

    // read chunks from the input file and pass them to the
    // worker threads.
    let mut buf = vec![0u8; 65536];
    for cell_num in 0..(cc.num_chunks as usize) {
        let (nbytes_chunk, nrec_chunk) = chunk::Chunk::<AlevinFryReadRecord>::read_header(&mut br);
        buf.resize(nbytes_chunk as usize, 0);
        buf.pwrite::<u32>(nbytes_chunk, 0)?;
        buf.pwrite::<u32>(nrec_chunk, 4)?;
        br.read_exact(&mut buf[8..]).unwrap();

        let mut bclone = (cell_num, buf.clone());
        // keep trying until we can push this payload
        while let Err(t) = q.push(bclone) {
            bclone = t;
            // no point trying to push if the queue is full
            while q.is_full() {}
        }
        pbar_inner.inc(1);
    }
    pbar_inner.finish();

    // wait for the worker threads to finish
    for h in thread_handles.drain(0..) {
        match h.join() {
            Ok(_) => {}
            Err(_e) => {
                info!(log, "thread panicked");
            }
        }
    }
    pbar_inner.finish_with_message("partitioned records into temporary files.");
    drop(q);

    // At this point, we are done with the "scatter"
    // phase of writing the records to the corresponding
    // intermediate files.  Now, we'll begin the gather
    // phase of collating the temporary files and merging
    // them into the final output file.

    for (i, temp_bucket) in temp_buckets.iter().enumerate() {
        // make sure we flush each temp bucket
        temp_bucket
            .2
            .bucket_writer
            .lock()
            .unwrap()
            .flush()
            .context("could not flush temporary output file!")?;
        // a sanity check that we have the correct number of records
        // and the expected number of bytes in each file
        let expected = temp_bucket.1;
        let observed = temp_bucket.2.num_records_written.load(Ordering::SeqCst);
        assert_eq!(expected, observed);

        let md = std::fs::metadata(parent.join(format!("bucket_{}.tmp", i)))?;
        let expected_bytes = temp_bucket.2.num_bytes_written.load(Ordering::SeqCst);
        let observed_bytes = md.len();
        assert_eq!(expected_bytes, observed_bytes);
    }

    //std::process::exit(1);

    // to hold the temp buckets threads will process
    let slack = (n_workers / 2).max(1_usize);
    let temp_bucket_queue_size = slack + n_workers;
    let fq = Arc::new(ArrayQueue::<(
        u32,
        u32,
        std::sync::Arc<libradicl::TempBucket>,
    )>::new(temp_bucket_queue_size));
    // the number of cells left to process
    let buckets_to_process = Arc::new(AtomicUsize::new(temp_buckets.len()));

    let pbar_gather = ProgressBar::new(temp_buckets.len() as u64);
    pbar_gather.set_style(sty);
    pbar_gather.tick();

    // for each worker, spawn off a thread
    for _worker in 0..n_workers {
        // each thread will need to access the work queue
        let in_q = fq.clone();
        // the output cache and correction map
        let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
        let mut cmap = HashMap::<u64, TempCellInfo, ahash::RandomState>::with_hasher(s);
        // alternative strategy
        // let mut cmap = HashMap::<u64, libradicl::CorrectedCbChunk, ahash::RandomState>::with_hasher(s);

        // the number of chunks remaining to be processed
        let buckets_remaining = buckets_to_process.clone();
        // and knowledge of the UMI and BC types
        let bc_type =
            rad_types::decode_int_type_tag(cc.bc_type).context("unknown barcode type id.")?;
        let umi_type =
            rad_types::decode_int_type_tag(cc.umi_type).context("unknown umi type id.")?;
        // have access to the input directory
        let input_dir: PathBuf = input_dir.clone();
        // the output file
        let owriter = owriter.clone();
        // and the progress bar
        let pbar_gather = pbar_gather.clone();

        // now, make the worker threads
        let handle = std::thread::spawn(move || {
            let mut local_chunks = 0u64;
            let parent = std::path::Path::new(&input_dir);
            // pop from the work queue until everything is
            // processed
            while buckets_remaining.load(Ordering::SeqCst) > 0 {
                if let Some(temp_bucket) = in_q.pop() {
                    buckets_remaining.fetch_sub(1, Ordering::SeqCst);
                    cmap.clear();

                    let fname = parent.join(format!("bucket_{}.tmp", temp_bucket.2.bucket_id));
                    // create a new handle for reading
                    let tfile = std::fs::File::open(&fname).expect("couldn't open temporary file.");
                    let mut treader = BufReader::new(tfile);

                    local_chunks += libradicl::collate_temporary_bucket_twopass(
                        &mut treader,
                        &bc_type,
                        &umi_type,
                        temp_bucket.1,
                        &owriter,
                        compress_out,
                        &mut cmap,
                    ) as u64;

                    // we don't need the file or reader anymore
                    drop(treader);
                    std::fs::remove_file(fname).expect("could not delete temporary file.");

                    pbar_gather.inc(1);
                }
            }
            local_chunks
        });
        thread_handles.push(handle);
    } // for each worker

    // push the temporary buckets onto the work queue to be dispatched
    // by the worker threads.
    for temp_bucket in temp_buckets {
        let mut bclone = temp_bucket.clone();
        // keep trying until we can push this payload
        while let Err(t) = fq.push(bclone) {
            bclone = t;
            // no point trying to push if the queue is full
            while fq.is_full() {}
        }
        let expected = temp_bucket.1;
        let observed = temp_bucket.2.num_records_written.load(Ordering::SeqCst);
        assert_eq!(expected, observed);
    }

    // wait for all of the workers to finish
    let mut num_output_chunks = 0u64;
    for h in thread_handles.drain(0..) {
        match h.join() {
            Ok(c) => {
                num_output_chunks += c;
            }
            Err(_e) => {
                info!(log, "thread panicked");
            }
        }
    }
    pbar_gather.finish_with_message("gathered all temp files.");

    // make sure we wrote the same number of records that our
    // file suggested we should.
    assert_eq!(total_allocated_records, total_to_collate);

    info!(
        log,
        "writing num output chunks ({}) to header",
        num_output_chunks.to_formatted_string(&Locale::en)
    );

    info!(
        log,
        "expected number of output chunks {}",
        expected_output_chunks.to_formatted_string(&Locale::en)
    );

    assert_eq!(
        expected_output_chunks,
        num_output_chunks,
        "expected to write {} chunks but wrote {}",
        expected_output_chunks.to_formatted_string(&Locale::en),
        num_output_chunks.to_formatted_string(&Locale::en),
    );

    owriter.lock().unwrap().flush()?;
    info!(
        log,
        "finished collating input rad file {:?}.",
        i_dir.join("map.rad")
    );
    Ok(())
}
*/
