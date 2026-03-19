/*
 * Copyright (c) 2020-2026 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use anyhow::{Context, anyhow};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use slog::{crit, info};
//use anyhow::{anyhow, Result};
use crate::constants as afconst;
use crate::utils::InternalVersionInfo;
use bio_types::strand::{Strand, StrandError};
use crossbeam_queue::ArrayQueue;
// use dashmap::DashMap;

use libradicl::chunk;
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{self, RadIntId};
use libradicl::record::{
    AlevinFryReadRecordT, AlevinFryReadRecordWithPositionT, AlevinFryRecordContext,
    CollatableMappedRecord, CollatableRecordHeader, ConvertiblePrimitiveInteger,
    HierarchicallyCollatable, KnownSize, MappedRecord, MultiBarcodeReadRecord,
    MultiBarcodeRecordContext, RecordHeader, ScLongReadRecordContext, ScLongReadRecordT,
};
use libradicl::schema::TempCellInfo;
use libradicl::collation::{CollationManifest, SampleGroup};

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
    let is_multi_barcode = mdata.get("multi_barcode")
        .and_then(|v| v.as_bool())
        .unwrap_or(false);

    if is_multi_barcode {
        // Multi-barcode: skip permit_freq.bin loading, go directly to collate_with_temp
        // which will detect the multi-barcode record type and use per-sample permit maps.
        // total_to_collate is computed from sample_info.json
        let sample_info_file = File::open(parent.join("sample_info.json"))
            .context("couldn't open sample_info.json for multi-barcode collation")?;
        let sample_info: serde_json::Value = serde_json::from_reader(sample_info_file)?;
        let total_to_collate = sample_info["matched_reads"]
            .as_u64()
            .unwrap_or(0);

        info!(log, "Multi-barcode mode: {} total reads to collate", total_to_collate);

        // tsv_map is empty — do_collate_multi_bc will build its own from per-sample data
        return collate_with_temp(
            input_dir,
            rad_dir,
            num_threads,
            max_records,
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

    collate_with_temp(
        input_dir,
        rad_dir,
        num_threads,
        max_records,
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
    per_sample_identity: &[std::collections::HashSet<u64, ahash::RandomState>],
    per_sample_lookup: &[libradicl::BarcodeLookupMap],
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
                collated.write_to_file(&s_path).expect("could not write collated file.");
                return;
            }
            Err(_) => {
                let s_path = parent.join("unmapped_bc_count_collated.bin");
                collated.write_to_file(&s_path).expect("could not write collated file.");
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

            // Correct cell BC (tiered: identity fast path, then per-sample 1-edit lookup)
            let corrected_cell = if sample_idx < per_sample_identity.len()
                && per_sample_identity[sample_idx].contains(&cell_bc)
            {
                cell_bc
            } else if sample_idx < per_sample_lookup.len() {
                match per_sample_lookup[sample_idx].find_neighbors(cell_bc, false) {
                    (Some(idx), 1) => per_sample_lookup[sample_idx].barcodes[idx],
                    _ => continue,
                }
            } else {
                continue
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
    correct_map: &Arc<HashMap<u64, u64>>,
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
                    do_collate_with_temp::<_, _, _, u64, AlevinFryReadRecordT<u64>>(
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
        KnownRecordType::RnaShortMultiBC(cell_bc_len, num_bc) => {
            info!(
                log,
                "record type is multi-barcode single-cell RNA-seq ({} barcode levels, cell BC len = {})",
                num_bc,
                cell_bc_len,
            );
            let parsing_context = prelude.get_record_context::<MultiBarcodeRecordContext>()?;

            // Choose collation mode based on metadata or CLI flag.
            // The two-round mode is more modular and generalizable to N levels;
            // the fast mode is a single-pass optimization for 2-level protocols.
            // Default to fast mode for now.
            //
            // TODO: read --collation-mode from CLI and pass through
            let use_fast_path = true;

            if use_fast_path {
                info!(log, "Using fast single-pass collation mode");
                do_collate_multi_bc_fast(
                    input_dir,
                    &rad_dir,
                    parsing_context,
                    prelude,
                    br,
                    end_header_pos,
                    num_threads,
                    max_records,
                    total_to_collate,
                    compress_out,
                    cmdline,
                    version,
                    log,
                )
            } else {
                info!(log, "Using two-round collation mode");
                do_collate_multi_bc_two_round(
                    input_dir,
                    &rad_dir,
                    parsing_context,
                    prelude,
                    br,
                    end_header_pos,
                    num_threads,
                    max_records,
                    total_to_collate,
                    compress_out,
                    cmdline,
                    version,
                    log,
                )
            }
        }
    }
}

/// Multi-barcode hierarchical collation — fast single-pass mode.
///
/// This function handles RAD files with multiple barcodes per read (e.g., 10x Flex).
/// It performs a single-pass collation where records are bucketed by corrected cell
/// barcode (grouped by sample), producing a hierarchically-ordered output.
///
/// The output is a single collated RAD file with chunks ordered as:
///   [sample_0/cell_0, sample_0/cell_1, ..., sample_1/cell_0, ...]
/// plus a `collation_manifest.bin` sidecar describing sample boundaries.
#[allow(clippy::too_many_arguments)]
fn do_collate_multi_bc_fast<P1, P2, A: Read + Seek>(
    input_dir: P1,
    rad_dir: P2,
    rec_context: MultiBarcodeRecordContext,
    prelude: RadPrelude,
    mut br: BufReader<A>,
    end_header_pos: u64,
    num_threads: u32,
    max_records: u32,
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

    // Read metadata
    let meta_data_file = File::open(parent.join("generate_permit_list.json"))
        .context("could not open the generate_permit_list.json file.")?;
    let mdata: serde_json::Value = serde_json::from_reader(&meta_data_file)?;
    let expected_ori: Strand = get_orientation(&mdata)
        .map_err(|e| anyhow!("Error reading strand info: {}", e))?;

    // Load sample_info.json to get sample metadata
    let sample_info_file = File::open(parent.join("sample_info.json"))
        .context("could not open sample_info.json — was generate-permit-list run with --sample-bc-list?")?;
    let sample_info: serde_json::Value = serde_json::from_reader(&sample_info_file)?;

    let num_samples = sample_info["num_samples"]
        .as_u64()
        .context("couldn't read num_samples from sample_info.json")? as usize;
    let sample_entries = sample_info["samples"]
        .as_array()
        .context("couldn't read samples array from sample_info.json")?;

    info!(log, "Loading permit maps for {} samples", num_samples);

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

    // Get cell barcode length from file tags
    let _cell_bc_tag = format!("b{}len", prelude.read_tags.tags.len().saturating_sub(2));
    // Try to read from the file_tag_map if available, otherwise use a default
    let _cell_bc_len: u32 = {
        let _ftm = prelude.file_tags.parse_tags_from_bytes(&mut std::io::Cursor::new(Vec::<u8>::new())).ok();
        // Use 16 as default cell BC length (standard for 10x)
        16u32
    };

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
            let freq_map: HashMap<u64, u64> =
                bincode::deserialize_from(freq_reader)
                    .map_err(|e| anyhow!("couldn't deserialize {}: {}", freq_path.display(), e))?;

            let valid_bcs: Vec<u64> = freq_map.keys().copied().collect();
            info!(
                log,
                "Sample '{}': {} valid cell barcodes from permit_freq.bin",
                sample_name,
                valid_bcs.len().to_formatted_string(&Locale::en),
            );

            for (bc, freq) in &freq_map {
                all_cell_freqs.push((*bc, *freq, sample_idx));
            }

            per_sample_valid_bcs.push(valid_bcs);
        } else {
            info!(log, "Sample '{}': no permit freq (no reads)", sample_name);
            per_sample_valid_bcs.push(Vec::new());
        }
    }

    // Sort cell frequencies: group by sample_idx first, then by frequency descending
    all_cell_freqs.sort_by(|a, b| {
        a.2.cmp(&b.2)
            .then(b.1.cmp(&a.1))
    });

    // Build the tsv_map equivalent: (composite_key, freq) pairs
    // composite_key encodes both sample and cell for ordering
    // We use (sample_idx << 48) | cell_bc as composite key for bucket ordering
    let tsv_map: Vec<(u64, u64)> = all_cell_freqs
        .iter()
        .map(|(bc, freq, sidx)| {
            let composite = ((*sidx as u64) << 48) | (*bc & 0x0000_FFFF_FFFF_FFFF);
            (composite, *freq)
        })
        .collect();

    let expected_output_chunks = tsv_map.len() as u64;
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

    // Build output_cache: composite_key -> TempBucket
    let mut output_cache = Arc::new(HashMap::<u64, Arc<libradicl::TempBucket>>::new());
    let _total_allocated_records: u64 = 0;
    let mut allocated_records: u64 = 0;
    let mut temp_buckets = vec![(
        0u32,
        0u32,
        Arc::new(libradicl::TempBucket::from_id_and_parent(0, parent)),
    )];
    let max_records_per_thread = (max_records / n_workers as u32) + 1;
    let mut num_bucket_chunks = 0u32;
    {
        let moutput_cache = Arc::make_mut(&mut output_cache);
        for rec in tsv_map.iter() {
            moutput_cache.insert(rec.0, temp_buckets.last().unwrap().2.clone());
            allocated_records += rec.1;
            num_bucket_chunks += 1;
            if allocated_records >= max_records_per_thread as u64 {
                temp_buckets.last_mut().unwrap().0 = num_bucket_chunks;
                temp_buckets.last_mut().unwrap().1 = allocated_records as u32;
                let tn = temp_buckets.len() as u32;
                temp_buckets.push((
                    0,
                    0,
                    Arc::new(libradicl::TempBucket::from_id_and_parent(tn, parent)),
                ));
                allocated_records = 0;
                num_bucket_chunks = 0;
            }
        }
    }
    if num_bucket_chunks > 0 {
        temp_buckets.last_mut().unwrap().0 = num_bucket_chunks;
        temp_buckets.last_mut().unwrap().1 = allocated_records as u32;
    }
    info!(log, "Generated {} temporary buckets for multi-barcode collation.", temp_buckets.len());

    let nbuckets = temp_buckets.len();
    let hdr = &prelude.hdr;

    // Progress bar
    let sty = ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}")
        .expect("invalid template")
        .progress_chars("##-");
    let pbar_inner = ProgressBar::with_draw_target(
        Some(hdr.num_chunks),
        ProgressDrawTarget::stderr_with_hz(5u8),
    );
    pbar_inner.set_style(sty.clone());
    pbar_inner.tick();

    // Work queue for scatter phase
    let _q = Arc::new(ArrayQueue::<(usize, Vec<u8>)>::new(4 * n_workers));
    let _chunks_to_process = Arc::new(AtomicUsize::new(hdr.num_chunks as usize));

    // SCATTER PHASE: Read chunks, correct both barcodes, write to temp buckets
    let _thread_handles: Vec<thread::JoinHandle<u64>> = Vec::with_capacity(n_workers);

    let min_rec_len = MultiBarcodeReadRecord::nbytes(1, &rec_context);
    let _loc_buffer_size =
        (min_rec_len * max_records as usize / (nbuckets * n_workers)).clamp(1000, 262_144usize);

    let sample_permit_map = Arc::new(sample_permit_map);
    let sample_bc_to_idx = Arc::new(sample_bc_to_idx);
    let rec_context = Arc::new(rec_context);

    // Build tiered cell correction structure:
    // Tier 1 (per-sample identity): HashSet of valid barcodes for each sample.
    //   A valid barcode corrects to itself within that sample. This is the fast path.
    //   The valid set differs per sample — the same barcode may be valid in one sample
    //   but not another (different cell populations per sample).
    // Tier 2 (per-sample BarcodeLookupMap): for barcodes NOT in the identity set,
    //   do on-the-fly 1-edit neighbor lookup. Much cheaper than pre-materializing
    //   the full correction HashMap (which can be 100M+ entries).
    let cell_correction = {
        // Get cell BC length from the first non-empty sample's freq file
        let cell_bclen = {
            let freq_path = sample_names.iter().find_map(|name| {
                let p = parent.join(format!("sample_{}/permit_freq.bin", name));
                if p.exists() { Some(p) } else { None }
            });
            if let Some(fp) = freq_path {
                let mut f = File::open(fp)?;
                let mut hdr = [0u8; 16];
                f.read_exact(&mut hdr)?;
                hdr.pread::<u64>(8).unwrap_or(16) as u32
            } else {
                16u32
            }
        };

        // Per-sample identity sets: valid barcodes that self-correct within each sample
        let per_sample_identity: Vec<std::collections::HashSet<u64, ahash::RandomState>> =
            per_sample_valid_bcs
                .iter()
                .map(|valid_bcs| valid_bcs.iter().copied().collect())
                .collect();

        // Per-sample BarcodeLookupMap for 1-edit correction
        let per_sample_lookup: Vec<libradicl::BarcodeLookupMap> = per_sample_valid_bcs
            .iter()
            .map(|valid_bcs| {
                libradicl::BarcodeLookupMap::new(valid_bcs.clone(), cell_bclen)
            })
            .collect();

        let total_identity: usize = per_sample_identity.iter().map(|s| s.len()).sum();
        let total_lookup: usize = per_sample_lookup.iter().map(|m| m.barcodes.len()).sum();
        info!(
            log,
            "Cell correction: {} total identity barcodes across {} samples (fast path), {} total in lookup maps (bc_len={})",
            total_identity.to_formatted_string(&Locale::en),
            num_samples,
            total_lookup.to_formatted_string(&Locale::en),
            cell_bclen,
        );
        Arc::new((per_sample_identity, per_sample_lookup))
    };

    // SCATTER PHASE: Multi-threaded scatter using header-peek pattern
    info!(log, "Starting scatter phase ({} worker threads)...", n_workers);

    let q = Arc::new(ArrayQueue::<(usize, Vec<u8>)>::new(4 * n_workers));
    let chunks_to_process = Arc::new(AtomicUsize::new(hdr.num_chunks as usize));
    let mut thread_handles: Vec<thread::JoinHandle<u64>> = Vec::with_capacity(n_workers);

    let min_rec_len = MultiBarcodeReadRecord::nbytes(1, &rec_context);
    let loc_buffer_size =
        (min_rec_len * max_records as usize / (nbuckets * n_workers)).clamp(1000, 262_144usize);

    let sample_permit_map = Arc::new(sample_permit_map);
    let sample_bc_to_idx = Arc::new(sample_bc_to_idx);
    let rec_context = Arc::new(rec_context);

    for _worker in 0..n_workers {
        let in_q = q.clone();
        let chunks_remaining = chunks_to_process.clone();
        let oc = output_cache.clone();
        let sample_map = sample_permit_map.clone();
        let sample_idx_map = sample_bc_to_idx.clone();
        let cell_corr = cell_correction.clone();
        let ctx = rec_context.clone();
        let loc_temp_buckets = temp_buckets.clone();
        let expected_ori_mfo: libradicl::rad_types::MappedFragmentOrientation =
            expected_ori.into();

        let handle = thread::spawn(move || {
            // Thread-local buffer backing: one contiguous allocation, split into per-bucket cursors
            let mut local_buffer_backing = vec![0u8; loc_buffer_size * nbuckets];
            let mut local_buffers: Vec<Cursor<&mut [u8]>> = Vec::with_capacity(nbuckets);
            let mut tslice = local_buffer_backing.as_mut_slice();
            for _ in 0..nbuckets {
                let (first, rest) = tslice.split_at_mut(loc_buffer_size);
                local_buffers.push(Cursor::new(first));
                tslice = rest;
            }

            let mut tbuf = vec![0u8; 4096]; // skip buffer for unmatched records
            let mut processed: u64 = 0;

            while chunks_remaining.load(Ordering::SeqCst) > 0 {
                if let Some((_chunk_num, buf)) = in_q.pop() {
                    chunks_remaining.fetch_sub(1, Ordering::SeqCst);

                    // Parse chunk using the same pattern as dump_corrected_cb_chunk_to_temp_file_generic
                    let mut reader = BufReader::new(&buf[..]);

                    // Read chunk header
                    let mut hdr_buf = [0u8; 8];
                    reader.read_exact(&mut hdr_buf).unwrap();
                    let _nbytes = hdr_buf.pread::<u32>(0).unwrap();
                    let nrec = hdr_buf.pread::<u32>(4).unwrap();

                    for _ in 0..nrec {
                        // Two-pass read: header first (peek), then body
                        let mut tup = MultiBarcodeReadRecord::from_bytes_collatable_header(
                            &mut reader, &ctx,
                        ).expect("could read multi-BC header");

                        // Extract sample BC (level 0) and look up correction
                        let sample_bc: u64 = tup.barcodes[0];
                        let corrected_sample = match sample_map.get(&sample_bc) {
                            Some(&cs) => cs,
                            None => {
                                // Skip alignment data for unmatched sample BC
                                let req_len = MultiBarcodeReadRecord::nbytes_aln(&ctx)
                                    * tup.naln() as usize;
                                if req_len > tbuf.len() { tbuf.resize(req_len, 0); }
                                reader.read_exact(&mut tbuf[..req_len]).unwrap();
                                if tbuf.len() > 4096 { tbuf.resize(4096, 0); tbuf.shrink_to_fit(); }
                                continue;
                            }
                        };
                        let sample_idx = match sample_idx_map.get(&corrected_sample) {
                            Some(&idx) => idx,
                            None => {
                                let req_len = MultiBarcodeReadRecord::nbytes_aln(&ctx)
                                    * tup.naln() as usize;
                                if req_len > tbuf.len() { tbuf.resize(req_len, 0); }
                                reader.read_exact(&mut tbuf[..req_len]).unwrap();
                                if tbuf.len() > 4096 { tbuf.resize(4096, 0); tbuf.shrink_to_fit(); }
                                continue;
                            }
                        };

                        // Look up cell BC correction (tiered: per-sample identity fast path, then 1-edit lookup)
                        let cell_bc: u64 = tup.collate_key();
                        let (ref per_sample_identity, ref per_sample_lookup) = *cell_corr;
                        let corrected_cell = if per_sample_identity[sample_idx].contains(&cell_bc) {
                            // Fast path: valid barcode in this sample, corrects to itself
                            cell_bc
                        } else {
                            // Slow path: 1-edit neighbor lookup in per-sample BarcodeLookupMap
                            match per_sample_lookup[sample_idx].find_neighbors(cell_bc, false) {
                                (Some(idx), 1) => {
                                    // Unique 1-edit neighbor found — correct to it
                                    per_sample_lookup[sample_idx].barcodes[idx]
                                }
                                _ => {
                                    // No correction available — skip alignment data
                                    let req_len = MultiBarcodeReadRecord::nbytes_aln(&ctx)
                                        * tup.naln() as usize;
                                    if req_len > tbuf.len() { tbuf.resize(req_len, 0); }
                                    reader.read_exact(&mut tbuf[..req_len]).unwrap();
                                    if tbuf.len() > 4096 { tbuf.resize(4096, 0); tbuf.shrink_to_fit(); }
                                    continue;
                                }
                            }
                        };

                        // Compute composite key for bucket lookup
                        let composite_key =
                            ((sample_idx as u64) << 48) | (corrected_cell & 0x0000_FFFF_FFFF_FFFF);

                        if let Some(bucket) = oc.get(&composite_key) {
                            // Read full record (alignments) with orientation filtering
                            let mut rr = MultiBarcodeReadRecord::from_bytes_with_header_retain_ori(
                                &mut reader, &mut tup, &ctx, &expected_ori_mfo,
                            );

                            if rr.is_empty() {
                                continue;
                            }

                            // Set corrected barcodes
                            rr.set_collation_key_at_level(0, corrected_sample);
                            rr.set_collate_key(corrected_cell);

                            let na = tup.naln() as usize;
                            let nb = MultiBarcodeReadRecord::nbytes(na as u32, &ctx);

                            let buffidx = bucket.bucket_id as usize;
                            let bcursor = &mut local_buffers[buffidx];
                            let len = bcursor.position() as usize;

                            // Flush if buffer would overflow
                            if len + nb >= loc_buffer_size {
                                let mut filebuf = bucket.bucket_writer.lock().unwrap();
                                filebuf.write_all(&bcursor.get_ref()[..len]).unwrap();
                                bcursor.set_position(0);
                            }

                            // Write corrected record to local buffer
                            rr.write(bcursor, &ctx).expect("can write record");
                            bucket.num_records_written.fetch_add(1, Ordering::SeqCst);
                            bucket.num_bytes_written.fetch_add(nb as u64, Ordering::SeqCst);
                            processed += 1;
                        } else {
                            // No bucket for this composite key — skip alignment data
                            let req_len = MultiBarcodeReadRecord::nbytes_aln(&ctx)
                                * tup.naln() as usize;
                            if req_len > tbuf.len() { tbuf.resize(req_len, 0); }
                            reader.read_exact(&mut tbuf[..req_len]).unwrap();
                            if tbuf.len() > 4096 { tbuf.resize(4096, 0); tbuf.shrink_to_fit(); }
                        }
                    }
                }
            }

            // Flush remaining local buffers
            for (bucket_id, lb) in local_buffers.iter().enumerate() {
                let len = lb.position() as usize;
                if len > 0 {
                    let mut filebuf = loc_temp_buckets[bucket_id].2.bucket_writer.lock().unwrap();
                    filebuf.write_all(&lb.get_ref()[..len]).unwrap();
                }
            }
            processed
        });
        thread_handles.push(handle);
    }

    // Main thread: read chunks and push to work queue
    let mut buf = vec![0u8; 65536];
    for cell_num in 0..(hdr.num_chunks as usize) {
        let (nbytes_chunk, nrec_chunk) =
            chunk::Chunk::<MultiBarcodeReadRecord>::read_header(&mut br);
        buf.resize(nbytes_chunk as usize, 0);
        buf.pwrite::<u32>(nbytes_chunk, 0).unwrap();
        buf.pwrite::<u32>(nrec_chunk, 4).unwrap();
        br.read_exact(&mut buf[8..]).unwrap();

        let mut bclone = (cell_num, buf.clone());
        while let Err(t) = q.push(bclone) {
            bclone = t;
            while q.is_full() {}
        }
        pbar_inner.inc(1);
    }
    pbar_inner.finish_and_clear();

    // Wait for workers
    let mut total_scattered: u64 = 0;
    for handle in thread_handles {
        total_scattered += handle.join().unwrap();
    }
    info!(log, "Scatter phase complete: {} records scattered", total_scattered.to_formatted_string(&Locale::en));

    // Flush all temp bucket writers
    for (_, _, bucket) in &temp_buckets {
        let bw = &bucket.bucket_writer;
        let mut writer = bw.lock().unwrap();
        writer.flush().unwrap();
    }

    // Correct unmapped barcode counts for multi-barcode data
    let unmapped_file = i_dir.join("unmapped_bc_count.bin");
    {
        let (ref per_sample_id, ref per_sample_lu) = *cell_correction;
        correct_unmapped_counts_multi_bc(
            &unmapped_file,
            &sample_permit_map,
            &sample_bc_to_idx,
            per_sample_id,
            per_sample_lu,
            parent,
        );
    }

    // GATHER PHASE: Merge temp buckets into final collated file
    info!(log, "Starting gather phase ({} worker threads)...", n_workers);

    let pbar_gather = ProgressBar::with_draw_target(
        Some(temp_buckets.len() as u64),
        ProgressDrawTarget::stderr_with_hz(5u8),
    );
    pbar_gather.set_style(sty);
    pbar_gather.tick();

    // Multi-threaded gather: workers pop temp buckets from a queue,
    // call collate_temporary_bucket_twopass_generic, delete temp files.
    let fq = Arc::new(ArrayQueue::<(u32, u32, Arc<libradicl::TempBucket>)>::new(
        (n_workers / 2).max(1) + n_workers,
    ));
    let buckets_to_process = Arc::new(AtomicUsize::new(temp_buckets.len()));

    let mut gather_handles: Vec<thread::JoinHandle<u64>> = Vec::with_capacity(n_workers);

    for _worker in 0..n_workers {
        let in_q = fq.clone();
        let buckets_remaining = buckets_to_process.clone();
        let input_dir_clone: PathBuf = parent.to_path_buf();
        let owriter_clone = owriter.clone();
        let pbar_clone = pbar_gather.clone();
        let ctx = rec_context.clone();

        let handle = thread::spawn(move || {
            let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
            let mut cmap = HashMap::<u64, TempCellInfo, ahash::RandomState>::with_hasher(s);
            let mut local_chunks = 0u64;
            let p = std::path::Path::new(&input_dir_clone);

            while buckets_remaining.load(Ordering::SeqCst) > 0 {
                if let Some(temp_bucket) = in_q.pop() {
                    buckets_remaining.fetch_sub(1, Ordering::SeqCst);
                    cmap.clear();

                    let nrec = temp_bucket.2.num_records_written.load(Ordering::SeqCst);
                    if nrec == 0 {
                        pbar_clone.inc(1);
                        continue;
                    }

                    let fname = p.join(format!("bucket_{}.tmp", temp_bucket.2.bucket_id));
                    let tfile = File::open(&fname).expect("couldn't open temporary file");
                    let mut treader = BufReader::new(tfile);

                    local_chunks += libradicl::collate_temporary_bucket_twopass_generic::<
                        u64, _, _, MultiBarcodeReadRecord,
                    >(
                        &mut treader, &ctx, nrec, &owriter_clone, compress_out, &mut cmap,
                    ) as u64;

                    drop(treader);
                    std::fs::remove_file(&fname).ok();
                    pbar_clone.inc(1);
                }
            }
            local_chunks
        });
        gather_handles.push(handle);
    }

    // Push temp buckets to the gather queue
    for tb in temp_buckets.iter() {
        let mut bclone = (tb.0, tb.1, tb.2.clone());
        while let Err(t) = fq.push(bclone) {
            bclone = t;
            while fq.is_full() {}
        }
    }

    // Wait for gather workers
    let mut total_output_chunks: u64 = 0;
    for handle in gather_handles {
        total_output_chunks += handle.join().unwrap();
    }
    pbar_gather.finish_and_clear();

    // Build manifest from the tsv_map (which encodes sample→cell mappings)
    let mut manifest = CollationManifest::new(vec![
        "sample".to_string(),
        "cell".to_string(),
    ]);

    // Group cells by sample from the tsv_map composite keys
    let mut sample_chunk_counts: HashMap<usize, (u64, u64)> = HashMap::new(); // sample_idx -> (num_cells, num_records)
    for (composite_key, freq) in &tsv_map {
        let sidx = (*composite_key >> 48) as usize;
        let entry = sample_chunk_counts.entry(sidx).or_insert((0, 0));
        entry.0 += 1;
        entry.1 += freq;
    }

    let mut chunk_offset: u64 = 0;
    let mut sample_indices: Vec<usize> = sample_chunk_counts.keys().copied().collect();
    sample_indices.sort();
    for sidx in sample_indices {
        let (num_cells, num_records) = sample_chunk_counts[&sidx];
        if sidx < sample_entries.len() {
            let bc_str = &sample_entries[sidx]["barcode"];
            let bc = u64::from_str_radix(
                bc_str.as_str().unwrap_or("0").trim_start_matches("0x"),
                16,
            ).unwrap_or(0);
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
        info!(log, "Backpatched num_chunks to {} in output file", total_output_chunks);
    }

    info!(
        log,
        "Multi-barcode collation complete: {} output chunks",
        total_output_chunks.to_formatted_string(&Locale::en),
    );

    Ok(())
}

/// Multi-barcode hierarchical collation — two-round mode.
///
/// Round 1: Reads the multi-barcode RAD file, corrects sample BCs, and writes
///          per-sample intermediate RAD files (same format, filtered to one sample).
/// Round 2: For each sample, runs standard cell-level collation using the existing
///          do_collate_with_temp machinery, then appends the output chunks to
///          the final collated file and builds the collation manifest.
///
/// This mode is more modular and generalizes to N barcode levels by adding more rounds.
#[allow(clippy::too_many_arguments)]
fn do_collate_multi_bc_two_round<P1, P2, A: Read + Seek>(
    input_dir: P1,
    rad_dir: P2,
    rec_context: MultiBarcodeRecordContext,
    prelude: RadPrelude,
    mut br: BufReader<A>,
    end_header_pos: u64,
    num_threads: u32,
    max_records: u32,
    _total_to_collate: u64,
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
    let input_dir: PathBuf = input_dir.into();
    let parent = std::path::Path::new(input_dir.as_path());

    // Read metadata
    let meta_data_file = File::open(parent.join("generate_permit_list.json"))?;
    let mdata: serde_json::Value = serde_json::from_reader(&meta_data_file)?;
    let expected_ori: Strand = get_orientation(&mdata)
        .map_err(|e| anyhow!("Error reading strand info: {}", e))?;

    // Load sample info
    let sample_info_file = File::open(parent.join("sample_info.json"))
        .context("could not open sample_info.json")?;
    let sample_info: serde_json::Value = serde_json::from_reader(&sample_info_file)?;
    let num_samples = sample_info["num_samples"].as_u64().unwrap() as usize;
    let sample_entries = sample_info["samples"].as_array().unwrap();

    // Load sample permit map
    let sample_map_file = File::open(parent.join("sample_permit_map.bin"))?;
    let sample_permit_map: HashMap<u64, u64> =
        bincode::deserialize_from(BufReader::new(sample_map_file))
            .map_err(|e| anyhow!("couldn't deserialize sample_permit_map.bin: {}", e))?;

    // Build sample_bc -> index mapping and sample names
    let mut sample_bc_to_idx: HashMap<u64, usize> = HashMap::new();
    let mut sample_names: Vec<String> = Vec::new();
    for (idx, entry) in sample_entries.iter().enumerate() {
        let bc_str = entry["barcode"].as_str().unwrap_or("0x0");
        let bc = u64::from_str_radix(bc_str.trim_start_matches("0x"), 16).unwrap_or(0);
        sample_bc_to_idx.insert(bc, idx);
        sample_names.push(entry["name"].as_str().unwrap_or(&format!("{:x}", bc)).to_string());
    }

    info!(log, "Two-round collation: Round 1 — scatter by sample...");

    // === ROUND 1: Scatter records into per-sample intermediate RAD files ===
    let round1_dir = parent.join("_collate_round1");
    std::fs::create_dir_all(&round1_dir)?;

    // Create per-sample RAD writers. Each intermediate file has the same header
    // as the input but will contain only that sample's records.
    let mut per_sample_writers: Vec<Option<BufWriter<File>>> = Vec::new();
    let mut per_sample_nchunks: Vec<u64> = vec![0; num_samples];
    let mut per_sample_nrecs: Vec<u64> = vec![0; num_samples];

    for name in sample_names.iter() {
        let sample_rad_path = round1_dir.join(format!("sample_{}.rad", name));
        let f = File::create(&sample_rad_path)?;
        let mut writer = BufWriter::with_capacity(1048576, f);

        // Copy the header from the input RAD (will backpatch num_chunks later)
        let mut rfile = File::open(&input_rad_path)?;
        let pos = br.get_mut().stream_position().unwrap() - (br.buffer().len() as u64);
        let mut hdr_bytes = vec![0u8; pos as usize];
        rfile.read_exact(&mut hdr_bytes)?;
        // Set num_chunks to 0 (will backpatch)
        let chunk_bytes = std::mem::size_of::<u64>() as u64;
        let nc_pos = (end_header_pos - chunk_bytes) as usize;
        hdr_bytes[nc_pos..nc_pos + 8].copy_from_slice(&0u64.to_le_bytes());
        writer.write_all(&hdr_bytes)?;

        per_sample_writers.push(Some(writer));
    }

    let _expected_ori_mfo: libradicl::rad_types::MappedFragmentOrientation = expected_ori.into();
    let hdr = &prelude.hdr;

    // Read chunks and scatter to per-sample files
    let pbar = ProgressBar::with_draw_target(
        Some(hdr.num_chunks),
        ProgressDrawTarget::stderr_with_hz(5u8),
    );
    pbar.set_style(
        ProgressStyle::with_template("[{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} chunks")
            .unwrap()
            .progress_chars("##-"),
    );

    for _ in 0..hdr.num_chunks {
        let c = chunk::Chunk::<MultiBarcodeReadRecord>::from_bytes(&mut br, &rec_context);
        // Group records by sample
        let mut per_sample_records: Vec<Vec<&MultiBarcodeReadRecord>> =
            (0..num_samples).map(|_| Vec::new()).collect();

        for read in &c.reads {
            let sample_bc: u64 = read.barcodes[0];
            if let Some(&corrected) = sample_permit_map.get(&sample_bc)
                && let Some(&idx) = sample_bc_to_idx.get(&corrected)
            {
                per_sample_records[idx].push(read);
            }
        }

        // Write per-sample chunks
        for (idx, records) in per_sample_records.iter().enumerate() {
            if records.is_empty() {
                continue;
            }
            if let Some(ref mut writer) = per_sample_writers[idx] {
                // Write chunk: nbytes (u32) + nrec (u32) + records
                let mut chunk_buf = Vec::new();
                // Placeholder for nbytes and nrec
                chunk_buf.extend_from_slice(&0u32.to_le_bytes());
                chunk_buf.extend_from_slice(&(records.len() as u32).to_le_bytes());
                for rec in records {
                    rec.write(&mut chunk_buf, &rec_context)?;
                }
                // Backpatch nbytes
                let nbytes = chunk_buf.len() as u32;
                chunk_buf[0..4].copy_from_slice(&nbytes.to_le_bytes());
                writer.write_all(&chunk_buf)?;
                per_sample_nchunks[idx] += 1;
                per_sample_nrecs[idx] += records.len() as u64;
            }
        }
        pbar.inc(1);
    }
    pbar.finish_and_clear();

    // Flush and backpatch num_chunks in per-sample files
    for (idx, writer_opt) in per_sample_writers.iter_mut().enumerate() {
        if let Some(writer) = writer_opt.take() {
            drop(writer); // flush
            let name = &sample_names[idx];
            let sample_rad_path = round1_dir.join(format!("sample_{}.rad", name));
            let mut f = std::fs::OpenOptions::new()
                .write(true)
                .open(&sample_rad_path)?;
            let _pos = br.get_mut().stream_position().unwrap() - (br.buffer().len() as u64);
            let chunk_bytes = std::mem::size_of::<u64>() as u64;
            let nc_pos = (end_header_pos - chunk_bytes) as usize;
            f.seek(std::io::SeekFrom::Start(nc_pos as u64))?;
            f.write_all(&per_sample_nchunks[idx].to_le_bytes())?;
            info!(
                log,
                "Round 1: sample '{}' — {} chunks, {} records",
                name, per_sample_nchunks[idx], per_sample_nrecs[idx],
            );
        }
    }

    info!(log, "Two-round collation: Round 2 — per-sample cell collation...");

    // === ROUND 2: Per-sample cell-level collation ===
    // For each sample, run the existing collation machinery on the intermediate file,
    // using the cell barcode (collate_key()) as the collation key.

    let cfname = if compress_out { "map.collated.rad.sz" } else { "map.collated.rad" };
    let oname = parent.join(cfname);
    if oname.exists() {
        std::fs::remove_file(&oname)?;
    }
    let ofile = File::create(&oname)?;
    let mut final_writer = BufWriter::with_capacity(1048576, ofile);

    // Write the header once (from the first sample's intermediate file, with total num_chunks)
    let _total_cells: u64 = {
        // We'll count cells during per-sample collation
        // For now, write header with placeholder
        let mut rfile = File::open(&input_rad_path)?;
        let pos_val = br.get_mut().stream_position().unwrap() - (br.buffer().len() as u64);
        let mut hdr_bytes = vec![0u8; pos_val as usize];
        rfile.read_exact(&mut hdr_bytes)?;
        // Will backpatch num_chunks after all samples are processed
        let nc_pos = (end_header_pos - std::mem::size_of::<u64>() as u64) as usize;
        hdr_bytes[nc_pos..nc_pos + 8].copy_from_slice(&0u64.to_le_bytes());

        if compress_out {
            let mut compressed =
                snap::write::FrameEncoder::new(Cursor::new(Vec::<u8>::with_capacity(hdr_bytes.len())));
            compressed.write_all(&hdr_bytes)?;
            let cbuf = compressed.into_inner()?;
            final_writer.write_all(cbuf.get_ref())?;
        } else {
            final_writer.write_all(&hdr_bytes)?;
        }
        0u64
    };
    drop(final_writer);

    let mut manifest = CollationManifest::new(vec!["sample".to_string(), "cell".to_string()]);
    let mut total_output_chunks: u64 = 0;

    let owriter = Arc::new(Mutex::new(BufWriter::with_capacity(
        1048576,
        std::fs::OpenOptions::new().append(true).open(&oname)?,
    )));

    for (sample_idx, name) in sample_names.iter().enumerate() {
        let sample_rad_path = round1_dir.join(format!("sample_{}.rad", name));
        if per_sample_nrecs[sample_idx] == 0 {
            info!(log, "Round 2: skipping sample '{}' (no records)", name);
            continue;
        }

        // Load this sample's cell permit maps
        let sample_dir = parent.join(format!("sample_{}", name));
        let freq_path = sample_dir.join("permit_freq.bin");
        let map_path = sample_dir.join("permit_map.bin");

        if !freq_path.exists() || !map_path.exists() {
            info!(log, "Round 2: skipping sample '{}' (no permit maps)", name);
            continue;
        }

        // Read permit_freq.bin to get cell barcode frequency list
        let freq_file = File::open(&freq_path)?;
        let mut freq_reader = BufReader::new(freq_file);
        let mut freq_hdr = [0u8; 16];
        freq_reader.read_exact(&mut freq_hdr)?;
        let cell_freq_map: HashMap<u64, u64> =
            bincode::deserialize_from(freq_reader)
                .map_err(|e| anyhow!("couldn't deserialize {}: {}", freq_path.display(), e))?;

        // Sort by frequency descending (same as standard collate)
        let mut tsv_map: Vec<(u64, u64)> = cell_freq_map.into_iter().collect();
        tsv_map.sort_by(|a, b| b.1.cmp(&a.1));

        let sample_total_cells = tsv_map.len() as u64;
        let sample_total_records: u64 = tsv_map.iter().map(|x| x.1).sum();

        // Open the sample's intermediate RAD file
        let sample_file = File::open(&sample_rad_path)?;
        let mut sample_br = BufReader::new(sample_file);
        let sample_prelude = RadPrelude::from_bytes(&mut sample_br)?;
        let _sample_ftm = sample_prelude.file_tags.parse_tags_from_bytes(&mut sample_br)?;
        let sample_ctx = sample_prelude.get_record_context::<MultiBarcodeRecordContext>()?;

        // Compute end_header_pos for the sample file
        let _sample_end_hdr_pos = sample_br.get_mut().stream_position().unwrap()
            - (sample_br.buffer().len() as u64);

        info!(
            log,
            "Round 2: collating sample '{}' — {} cells, {} records",
            name, sample_total_cells, sample_total_records,
        );

        // Use the existing collation machinery for cell-level collation.
        // We use collate_temporary_bucket_twopass_generic indirectly by
        // building temp buckets and scattering, same as do_collate_with_temp.
        let chunk_start = total_output_chunks;

        // Load cell permit map
        let cell_map_file = File::open(&map_path)?;
        let cell_correct_map: Arc<HashMap<u64, u64>> =
            Arc::new(bincode::deserialize_from(BufReader::new(cell_map_file))
                .map_err(|e| anyhow!("couldn't deserialize {}: {}", map_path.display(), e))?);

        // Build output_cache for this sample's cells
        let n_workers = if num_threads > 1 { (num_threads - 1) as usize } else { 1 };
        let max_records_per_thread = (max_records / n_workers as u32) + 1;

        let mut sample_output_cache = HashMap::<u64, Arc<libradicl::TempBucket>>::new();
        let mut sample_temp_buckets = vec![(
            0u32, 0u32,
            Arc::new(libradicl::TempBucket::from_id_and_parent(0, &round1_dir)),
        )];
        let mut alloc = 0u64;
        let mut nbchunks = 0u32;

        for (bc, freq) in &tsv_map {
            sample_output_cache.insert(*bc, sample_temp_buckets.last().unwrap().2.clone());
            alloc += freq;
            nbchunks += 1;
            if alloc >= max_records_per_thread as u64 {
                sample_temp_buckets.last_mut().unwrap().0 = nbchunks;
                sample_temp_buckets.last_mut().unwrap().1 = alloc as u32;
                let tn = sample_temp_buckets.len() as u32;
                sample_temp_buckets.push((
                    0, 0,
                    Arc::new(libradicl::TempBucket::from_id_and_parent(tn, &round1_dir)),
                ));
                alloc = 0;
                nbchunks = 0;
            }
        }
        if nbchunks > 0 {
            sample_temp_buckets.last_mut().unwrap().0 = nbchunks;
            sample_temp_buckets.last_mut().unwrap().1 = alloc as u32;
        }

        let sample_output_cache = Arc::new(sample_output_cache);
        let nbuckets = sample_temp_buckets.len();

        // Scatter phase for this sample
        let _expected_ori_mfo2: libradicl::rad_types::MappedFragmentOrientation =
            expected_ori.into();
        let min_rec_len = MultiBarcodeReadRecord::nbytes(1, &sample_ctx);
        let _total_sample_recs = per_sample_nrecs[sample_idx];
        let loc_buffer_size = (min_rec_len * max_records as usize / (nbuckets * n_workers.max(1)))
            .clamp(1000, 262_144usize);

        // Single-threaded scatter for simplicity (per-sample files are smaller)
        for _ in 0..per_sample_nchunks[sample_idx] {
            let mut hdr_buf = [0u8; 8];
            sample_br.read_exact(&mut hdr_buf).unwrap();
            let nbytes_chunk = u32::from_le_bytes([hdr_buf[0], hdr_buf[1], hdr_buf[2], hdr_buf[3]]);
            let _nrec_chunk = u32::from_le_bytes([hdr_buf[4], hdr_buf[5], hdr_buf[6], hdr_buf[7]]);

            // Read the chunk data and use libradicl's dump function
            let mut buf = vec![0u8; nbytes_chunk as usize];
            buf[0..8].copy_from_slice(&hdr_buf);
            sample_br.read_exact(&mut buf[8..]).unwrap();

            let mut nbr = BufReader::new(Cursor::new(&buf[..]));
            let mut backing = vec![0u8; loc_buffer_size * nbuckets];
            let mut local_buffers: Vec<Cursor<&mut [u8]>> = Vec::with_capacity(nbuckets);
            let mut tslice = backing.as_mut_slice();
            for _ in 0..nbuckets {
                let (first, rest) = tslice.split_at_mut(loc_buffer_size);
                local_buffers.push(Cursor::new(first));
                tslice = rest;
            }

            libradicl::dump_corrected_cb_chunk_to_temp_file_generic::<u64, _, MultiBarcodeReadRecord>(
                &mut nbr,
                &sample_ctx,
                &cell_correct_map,
                &expected_ori,
                &sample_output_cache,
                &mut local_buffers,
                loc_buffer_size,
            );

            // Flush local buffers
            for (_composite_key, bucket) in sample_output_cache.iter() {
                let bid = bucket.bucket_id as usize;
                if bid < local_buffers.len() {
                    let pos = local_buffers[bid].position() as usize;
                    if pos > 0 {
                        let bw = &bucket.bucket_writer;
                        let mut w = bw.lock().unwrap();
                        w.write_all(&local_buffers[bid].get_ref()[..pos]).unwrap();
                        local_buffers[bid].set_position(0);
                    }
                }
            }
        }

        // Flush all sample temp buckets
        for (_, _, bucket) in &sample_temp_buckets {
            let bw = &bucket.bucket_writer;
            let mut w = bw.lock().unwrap();
            w.flush().unwrap();
        }

        // Gather phase for this sample
        let mut sample_chunks: u64 = 0;
        for (_, _, bucket) in &sample_temp_buckets {
            let nrec = bucket.num_records_written.load(Ordering::SeqCst);
            if nrec == 0 { continue; }

            let bid = bucket.bucket_id;
            let tpath = round1_dir.join(format!("bucket_{}.tmp", bid));
            if !tpath.exists() { continue; }
            let tfile = File::open(&tpath)?;
            let mut treader = BufReader::new(tfile);

            let mut cmap: HashMap<u64, TempCellInfo, ahash::RandomState> = HashMap::default();
            let nc = libradicl::collate_temporary_bucket_twopass_generic::<
                u64, _, _, MultiBarcodeReadRecord,
            >(&mut treader, &sample_ctx, nrec, &owriter, compress_out, &mut cmap);

            sample_chunks += nc as u64;
            std::fs::remove_file(&tpath).ok();
        }

        total_output_chunks += sample_chunks;

        // Add to manifest
        let bc_str = sample_entries[sample_idx]["barcode"].as_str().unwrap_or("0x0");
        let bc = u64::from_str_radix(bc_str.trim_start_matches("0x"), 16).unwrap_or(0);
        manifest.add_sample_group(SampleGroup {
            key: bc,
            name: Some(name.clone()),
            chunk_start,
            num_chunks: sample_chunks,
            num_records: per_sample_nrecs[sample_idx],
        });

        info!(log, "Round 2: sample '{}' — {} output chunks", name, sample_chunks);
    }

    // Flush final output
    if let Ok(mut oput) = owriter.lock() {
        oput.flush()?;
    }

    // Backpatch total num_chunks in the output file header
    {
        let mut f = std::fs::OpenOptions::new().write(true).open(&oname)?;
        let nc_pos = (end_header_pos - std::mem::size_of::<u64>() as u64) as usize;
        if !compress_out {
            f.seek(std::io::SeekFrom::Start(nc_pos as u64))?;
            f.write_all(&total_output_chunks.to_le_bytes())?;
        }
        // For compressed output, backpatching is not straightforward; skip for now.
    }

    // Write collation manifest
    let manifest_path = parent.join("collation_manifest.bin");
    manifest.write_to_file(&manifest_path)?;

    // Write collate metadata
    {
        let collate_meta = json!({
            "cmd": cmdline,
            "version_str": version,
            "compressed_output": compress_out,
            "multi_barcode": true,
            "num_samples": num_samples,
            "collation_mode": "two-round",
        });
        let cm_path = parent.join("collate.json");
        let mut cm_file = File::create(cm_path)?;
        serde_json::to_writer_pretty(&mut cm_file, &collate_meta)?;
    }

    // Cleanup round 1 intermediate files
    for name in &sample_names {
        let p = round1_dir.join(format!("sample_{}.rad", name));
        std::fs::remove_file(&p).ok();
    }
    std::fs::remove_dir(&round1_dir).ok();

    info!(
        log,
        "Two-round collation complete: {} samples, {} total output chunks",
        manifest.sample_groups.len(),
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
