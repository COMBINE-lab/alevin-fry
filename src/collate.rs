/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use anyhow::{anyhow, Context};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use slog::{crit, info};
//use anyhow::{anyhow, Result};
use crate::constants as afconst;
use crate::utils::InternalVersionInfo;
use crossbeam_queue::ArrayQueue;
// use dashmap::DashMap;

use libradicl::chunk;
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types;
use libradicl::record::AtacSeqReadRecord;
use libradicl::schema::TempCellInfo;

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
            return Err(anyhow!("The generate_permit_list.json file does not contain a version_str field. Please re-run the generate-permit-list step with a newer version of alevin-fry"));
        }
    };

    if let Err(es) = calling_version.is_compatible_with(&vd) {
        return Err(anyhow!(es));
    }

    // if only an *old* version of the permit_freq is present, then complain and exit
    if parent.join("permit_freq.tsv").exists() && !parent.join("permit_freq.bin").exists() {
        crit!(log, "The file permit_freq.bin doesn't exist, please rerun alevin-fry generate-permit-list command.");
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
        crit!(log,
               "The permit_freq.bin file had version {}, but this version of alevin-fry requires version {}",
               freq_file_version, afconst::PERMIT_FILE_VER
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

    /*
    let est_num_rounds = (total_to_collate as f64 / max_records as f64).ceil() as u64;
    info!(
    log,
    "estimated that collation would require {} passes over input.", est_num_rounds
    );
    // if est_num_rounds > 2 {
    info!(log, "executing temporary file scatter-gather strategy.");
    */

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

    /*} else {
    info!(log, "executing multi-pass strategy.");
    collate_in_memory_multipass(
        input_dir,
        rad_dir,
        num_threads,
        max_records,
        tsv_map,
        total_to_collate,
        log,
    )
    }*/
}

#[derive(Debug)]
enum FilterType {
    Filtered,
    Unfiltered,
}

fn get_filter_type(mdata: &serde_json::Value, log: &slog::Logger) -> FilterType {
    if let Some(fts) = mdata.get("permit-list-type") {
        let ft = match fts.as_str() {
            Some("unfiltered") => FilterType::Unfiltered,
            Some("filtered") => FilterType::Filtered,
            _ => FilterType::Filtered,
        };
        ft
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

fn get_num_chunks(mdata: &serde_json::Value, log: &slog::Logger) -> anyhow::Result<u64> {
    if let Some(mar) = mdata.get("num-chunks") {
        match mar.as_u64() {
            Some(mv) => Ok(mv),
            _ => Err(anyhow!("Error parsing num-chunks")),
        }
    } else {
        info!(log, "num-chunks key not present in JSON file;");
        Err(anyhow!("num-chunks key not present"))
    }
}

fn correct_unmapped_counts(
    correct_map: &Arc<HashMap<u64, u64>>,
    unmapped_file: &std::path::Path,
    parent: &std::path::Path,
) {
    let i_file = File::open(unmapped_file).unwrap();
    let mut br = BufReader::new(i_file);

    // enough to hold a key value pair (a u64 key and u32 value)
    let mut rbuf = [0u8; std::mem::size_of::<u64>() + std::mem::size_of::<u32>()];

    let mut unmapped_count: HashMap<u64, u32> = HashMap::new();

    // pre-populate the output map with all valid keys
    // keys (corrected barcodes) with no unmapped reads
    // will simply have a value of 0.
    //for (&_ubc, &cbc) in correct_map.iter() {
    //    unmapped_count.entry(cbc).or_insert(0);
    //}

    // collect all of the information from the existing
    // serialized map (that may contain repeats)
    while br.read_exact(&mut rbuf[..]).is_ok() {
        let k = rbuf.pread::<u64>(0).unwrap();
        let v = rbuf.pread::<u32>(std::mem::size_of::<u64>()).unwrap();
        // get the corrected key for the raw key
        if let Some((&_rk, &ck)) = correct_map.get_key_value(&k) {
            *unmapped_count.entry(ck).or_insert(0) += v;
        }
    }

    let s_path = parent.join("unmapped_bc_count_collated.bin");
    let s_file = std::fs::File::create(s_path).expect("could not create serialization file.");
    let mut s_writer = BufWriter::new(&s_file);
    bincode::serialize_into(&mut s_writer, &unmapped_count)
        .expect("couldn't serialize corrected unmapped bc count.");
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

    let filter_type = get_filter_type(&mdata, log);
    let most_ambig_record = get_most_ambiguous_record(&mdata, log);
    let num_chunks = get_num_chunks(&mdata, log)?;

    // log the filter type
    info!(log, "filter_type = {:?}", filter_type);
    info!(
        log,
        "collated rad file {} be compressed",
        if compress_out { "will" } else { "will not" }
    );
    // because :
    // https://superuser.com/questions/865710/write-to-newfile-vs-overwriting-performance-issue
    let cfname = if compress_out {
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
        "paired : {:?}, ref_count : {}, num_chunks : {}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        num_chunks.to_formatted_string(&Locale::en)
    );

    // file-level
    let fl_tags = rad_types::TagSection::from_bytes(&mut br)?;
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = rad_types::TagSection::from_bytes(&mut br)?;
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = rad_types::TagSection::from_bytes(&mut br)?;
    info!(log, "read {:?} alignment-level tags", al_tags.tags.len());

    // create the prelude and rebind the variables we need
    let prelude = RadPrelude::from_header_and_tag_sections(hdr, fl_tags, rl_tags, al_tags);
    let rl_tags = &prelude.read_tags;

    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br);
    info!(log, "File-level tag values {:?}", file_tag_map);

    let bct = rl_tags.tags[0].typeid;

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

    let cc = chunk::ChunkConfigAtac {
        num_chunks,
        bc_type: libradicl::rad_types::encode_type_tag(bct).expect("valid barcode tag type"),
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
                    libradicl::dump_corrected_cb_chunk_to_temp_file_atac(
                        &mut nbr,
                        &bc_type,
                        &correct_map,
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
        let (nbytes_chunk, nrec_chunk) = chunk::Chunk::<AtacSeqReadRecord>::read_header(&mut br);
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
        // and knowledge of the BC types
        let bc_type =
            rad_types::decode_int_type_tag(cc.bc_type).context("unknown barcode type id.")?;

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

                    local_chunks += libradicl::collate_temporary_bucket_twopass_atac(
                        &mut treader,
                        &bc_type,
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
