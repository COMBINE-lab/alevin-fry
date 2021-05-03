/*
 * Copyright (c) 2020-2021 Rob Patro, Avi Srivastava, Hirak Sarkar, Dongze He, Mohsen Zakeri.
 *
 * This file is part of alevin-fry
 * (see https://github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

extern crate bincode;
extern crate indicatif;
extern crate slog;
use crate as libradicl;

use self::indicatif::{ProgressBar, ProgressStyle};
use self::slog::{crit, info};
//use anyhow::{anyhow, Result};
use crate::utils::InternalVersionInfo;
use bio_types::strand::{Strand, StrandError};
use crossbeam_queue::ArrayQueue;
// use dashmap::DashMap;
use self::libradicl::schema::TempCellInfo;
use num_format::{Locale, ToFormattedString};
use scroll::{Pread, Pwrite};
use serde_json::json;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::{BufWriter, Cursor, Read, Seek, SeekFrom, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

pub fn collate(
    input_dir: String,
    rad_dir: String,
    num_threads: u32,
    max_records: u32,
    compress_out: bool,
    version_str: &str,
    //expected_ori: Strand,
    log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {
    let parent = std::path::Path::new(&input_dir);

    // open the metadata file and read the json
    let gpl_path = parent.join("generate_permit_list.json");
    let meta_data_file =
        File::open(&gpl_path).expect(&format!("Could not open the file {:?}.", gpl_path)[..]);
    let mdata: serde_json::Value = serde_json::from_reader(meta_data_file)?;

    let calling_version = InternalVersionInfo::from_str(version_str);
    let vd: InternalVersionInfo;
    match mdata.get("version_str") {
        Some(vs) => match vs.as_str() {
            Some(s) => {
                vd = InternalVersionInfo::from_str(s);
            }
            None => {
                return Err("The version_str field must be a string".into());
            }
        },
        None => {
            return Err("The generate_permit_list.json file does not contain a version_str field. Please re-run the generate-permit-list step with a newer version of alevin-fry".into());
        }
    };

    if let Err(es) = calling_version.is_compatible_with(&vd) {
        return Err(es.into());
    }

    type TsvRec = (u64, u64);
    let mut tsv_map = Vec::<TsvRec>::new(); //HashMap::<u64, u64>::new();

    let freq_file =
        std::fs::File::open(parent.join("permit_freq.tsv")).expect("couldn't open file");
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(freq_file);

    let mut total_to_collate = 0;
    for result in rdr.deserialize() {
        let record: TsvRec = result?;
        tsv_map.push(record);
        total_to_collate += record.1;
    }
    // sort this so that we deal with largest cells (by # of reads) first
    // sort in _descending_ order by count.
    quickersort::sort_by_key(&mut tsv_map[..], |&a: &(u64, u64)| std::cmp::Reverse(a.1));

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

fn correct_unmapped_counts(
    correct_map: &Arc<HashMap<u64, u64>>,
    unmapped_file: &std::path::Path,
    parent: &std::path::Path,
) {
    let i_file = File::open(&unmapped_file).unwrap();
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
    let s_file = std::fs::File::create(&s_path).expect("could not create serialization file.");
    let mut s_writer = BufWriter::new(&s_file);
    bincode::serialize_into(&mut s_writer, &unmapped_count)
        .expect("couldn't serialize corrected unmapped bc count.");
}

#[allow(clippy::too_many_arguments)]
pub fn collate_with_temp(
    input_dir: String,
    rad_dir: String,
    num_threads: u32,
    max_records: u32,
    tsv_map: Vec<(u64, u64)>,
    total_to_collate: u64,
    compress_out: bool,
    log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {
    // the number of corrected cells we'll write
    let expected_output_chunks = tsv_map.len() as u64;
    // the parent input directory
    let parent = std::path::Path::new(&input_dir);

    let n_workers = if num_threads > 1 {
        (num_threads - 1) as usize
    } else {
        1
    };

    // open the metadata file and read the json
    let meta_data_file = File::open(parent.join("generate_permit_list.json"))
        .expect("could not open the generate_permit_list.json file.");
    let mdata: serde_json::Value = serde_json::from_reader(&meta_data_file)?;

    // velo_mode
    let velo_mode = mdata["velo_mode"].as_bool().unwrap();
    let expected_ori: Strand;
    match get_orientation(&mdata) {
        Ok(o) => {
            expected_ori = o;
        }
        Err(e) => {
            crit!(
                log,
                "Error reading strand info from {:#?} :: {}",
                &meta_data_file,
                e
            );
            return Err(e.into());
        }
    }

    let filter_type = get_filter_type(&mdata, log);

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
            "compressed_output" : compress_out,
        });

        let cm_path = parent.join("collate.json");
        let mut cm_file = std::fs::File::create(&cm_path).expect("could not create metadata file.");

        let cm_info_string =
            serde_json::to_string_pretty(&collate_meta).expect("could not format json.");
        cm_file
            .write_all(cm_info_string.as_bytes())
            .expect("cannot write to collate.json file");
    }

    let oname = parent.join(cfname);
    if oname.exists() {
        std::fs::remove_file(oname)?;
    }

    let ofile = File::create(parent.join(cfname)).unwrap();
    let owriter = Arc::new(Mutex::new(BufWriter::with_capacity(1048576, ofile)));

    let i_dir = std::path::Path::new(&rad_dir);

    if !i_dir.exists() {
        crit!(log, "the input RAD path {} does not exist", rad_dir);
        return Err("invalid input".into());
    }

    let input_rad_path = i_dir.join("map.rad");
    let i_file = File::open(&input_rad_path).unwrap();
    let mut br = BufReader::new(i_file);

    let hdr = libradicl::RadHeader::from_bytes(&mut br);

    // the exact position at the end of the header,
    // precisely sizeof(u64) bytes beyond the num_chunks field.
    let end_header_pos =
        br.get_ref().seek(SeekFrom::Current(0)).unwrap() - (br.buffer().len() as u64);

    info!(
        log,
        "paired : {:?}, ref_count : {}, num_chunks : {}, expected_ori : {:?}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en),
        expected_ori
    );

    // file-level
    let fl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    let ft_vals = libradicl::FileTags::from_bytes(&mut br);
    info!(log, "File-level tag values {:?}", ft_vals);

    let bct = rl_tags.tags[0].typeid;
    let umit = rl_tags.tags[1].typeid;

    // the exact position at the end of the header + file tags
    let pos = br.get_ref().seek(SeekFrom::Current(0)).unwrap() - (br.buffer().len() as u64);

    // copy the header
    {
        // we want to copy up to the end of the header
        // minus the num chunks (sizeof u64), and then
        // write the actual number of chunks we expect.
        let chunk_bytes = std::mem::size_of::<u64>() as u64;
        let take_pos = end_header_pos - chunk_bytes;

        // This temporary file pointer and buffer will be dropped
        // at the end of this block (scope).
        let mut rfile = File::open(&input_rad_path).unwrap();
        let mut hdr_buf = Cursor::new(vec![0u8; pos as usize]);

        rfile
            .read_exact(&mut hdr_buf.get_mut())
            .expect("couldn't read input file header");
        hdr_buf.set_position(take_pos);
        hdr_buf
            .write_all(&expected_output_chunks.to_le_bytes())
            .expect("couldn't write num_chunks");
        hdr_buf.set_position(0);

        // compress the header buffer to a compressed buffer
        if compress_out {
            let mut compressed_buf =
                snap::write::FrameEncoder::new(Cursor::new(Vec::<u8>::with_capacity(pos as usize)));
            compressed_buf
                .write_all(hdr_buf.get_ref())
                .expect("could not compress the output header.");
            hdr_buf = compressed_buf
                .into_inner()
                .expect("couldn't unwrap the FrameEncoder.");
            hdr_buf.set_position(0);
        }

        if let Ok(mut oput) = owriter.lock() {
            oput.write_all(hdr_buf.get_ref())
                .expect("could not write the output header.");
        }
    }

    // get the correction map
    let cmfile = std::fs::File::open(parent.join("permit_map.bin")).unwrap();
    let correct_map: Arc<HashMap<u64, u64>> = Arc::new(bincode::deserialize_from(&cmfile).unwrap());

    // NOTE: the assumption of where the unmapped file will be
    // should be robustified
    let unmapped_file = i_dir.join("unmapped_bc_count.bin");
    correct_unmapped_counts(&correct_map, &unmapped_file, &parent);

    info!(
        log,
        "deserialized correction map of length : {}",
        correct_map.len().to_formatted_string(&Locale::en)
    );

    let cc = libradicl::ChunkConfig {
        num_chunks: hdr.num_chunks,
        bc_type: bct,
        umi_type: umit,
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
        .progress_chars("╢▌▌░╟");

    let pbar_inner = ProgressBar::new(cc.num_chunks);
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
    let num_threads = n_workers as usize;
    let loc_buffer_size =
        (1000_usize.max((min_rec_len * max_rec) / (num_buckets * num_threads))).min(131072_usize);

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
        let bc_type = libradicl::decode_int_type_tag(cc.bc_type).expect("unknown barcode type id.");
        let umi_type =
            libradicl::decode_int_type_tag(cc.umi_type).expect("unknown barcode type id.");
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
    pbar_inner.set_message(&format!(
        "processing {} / {} total records",
        total_allocated_records, total_to_collate
    ));

    // read chunks from the input file and pass them to the
    // worker threads.
    let mut buf = vec![0u8; 65536];
    for cell_num in 0..(cc.num_chunks as usize) {
        let (nbytes_chunk, nrec_chunk) = libradicl::Chunk::read_header(&mut br);
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
            .expect("could not flush temporary output file!");
        // a sanity check that we have the correct number of records
        // and the expected number of bytes in each file
        let expected = temp_bucket.1;
        let observed = temp_bucket.2.num_records_written.load(Ordering::SeqCst);
        assert!(expected == observed);

        let md = std::fs::metadata(parent.join(&format!("bucket_{}.tmp", i)))?;
        let expected_bytes = temp_bucket.2.num_bytes_written.load(Ordering::SeqCst);
        let observed_bytes = md.len();
        assert!(expected_bytes == observed_bytes);
    }

    //std::process::exit(1);

    // to hold the temp buckets threads will process
    let slack = ((n_workers / 2) as usize).max(1_usize);
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
        let bc_type = libradicl::decode_int_type_tag(cc.bc_type).expect("unknown barcode type id.");
        let umi_type =
            libradicl::decode_int_type_tag(cc.umi_type).expect("unknown barcode type id.");
        // have access to the input directory
        let input_dir = input_dir.clone();
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

                    let fname = parent.join(&format!("bucket_{}.tmp", temp_bucket.2.bucket_id));
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
        assert!(expected == observed);
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
    assert!(total_allocated_records == total_to_collate);

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

/*
pub fn collate_in_memory_multipass(
    input_dir: String,
    rad_dir: String,
    num_threads: u32,
    max_records: u32,
    tsv_map: Vec<(u64, u64)>,
    total_to_collate: u64,
    log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {
    let parent = std::path::Path::new(&input_dir);

    // open the metadata file and read the json
    let meta_data_file = File::open(parent.join("generate_permit_list.json"))
        .expect("could not open the generate_permit_list.json file.");
    let mdata: serde_json::Value = serde_json::from_reader(meta_data_file)?;

    let velo_mode = mdata["velo_mode"].as_bool().unwrap();
    // info!(log, "velo_mode = {:?}", velo_mode);

    let expected_ori = get_orientation(&mdata, log);
    info!(log, "expected_ori = {:?}", expected_ori);

    let filter_type = get_filter_type(&mdata, log);
    info!(log, "filter_type = {:?}", filter_type);
    match filter_type {
        FilterType::Unfiltered => {
            unimplemented!();
        }
        FilterType::Filtered => {}
    };

    // because :
    // https://superuser.com/questions/865710/write-to-newfile-vs-overwriting-performance-issue
    let cfname = if velo_mode {
        "velo.map.collated.rad"
    } else {
        "map.collated.rad"
    };

    let oname = parent.join(cfname);
    if oname.exists() {
        std::fs::remove_file(oname)?;
    }

    let mut ofile = File::create(parent.join(cfname)).unwrap();
    let i_dir = std::path::Path::new(&rad_dir);

    if !i_dir.exists() {
        crit!(log, "the input RAD path {} does not exist", rad_dir);
        std::process::exit(1);
    }

    let i_file = File::open(i_dir.join("map.rad")).unwrap();
    let mut br = BufReader::new(i_file);

    let hdr = libradicl::RadHeader::from_bytes(&mut br);

    let end_header_pos =
        br.get_ref().seek(SeekFrom::Current(0)).unwrap() - (br.buffer().len() as u64);

    info!(
        log,
        "paired : {:?}, ref_count : {:?}, num_chunks : {:?}, expected_ori : {:?}",
        hdr.is_paired != 0,
        hdr.ref_count,
        hdr.num_chunks,
        expected_ori
    );
    // file-level
    let fl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = libradicl::TagSection::from_bytes(&mut br);
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    let ft_vals = libradicl::FileTags::from_bytes(&mut br);
    info!(log, "File-level tag values {:?}", ft_vals);

    let bct = rl_tags.tags[0].typeid;
    let umit = rl_tags.tags[1].typeid;

    let pos = br.get_ref().seek(SeekFrom::Current(0)).unwrap() - (br.buffer().len() as u64);

    // copy the header
    {
        br.get_mut()
            .seek(SeekFrom::Start(0))
            .expect("could not get read pointer.");
        let mut br2 = BufReader::new(br.get_ref());
        std::io::copy(&mut br2.by_ref().take(pos), &mut ofile).expect("couldn't copy header.");
    }

    // make sure that the buffer is empty
    // and that br starts reading from exactly
    // where we expect.
    if !br.buffer().is_empty() {
        br.consume(br.buffer().len());
    }

    br.get_mut()
        .seek(SeekFrom::Start(pos))
        .expect("could not get read pointer.");

    // get the correction map
    let cmfile = std::fs::File::open(parent.join("permit_map.bin")).unwrap();
    let correct_map: Arc<HashMap<u64, u64>> = Arc::new(bincode::deserialize_from(&cmfile).unwrap());

    let unmapped_file = i_dir.join("unmapped_bc_count.bin");
    correct_unmapped_counts(&correct_map, &unmapped_file, &parent);

    info!(
        log,
        "deserialized correction map of length : {:?}",
        correct_map.len()
    );

    let cc = libradicl::ChunkConfig {
        num_chunks: hdr.num_chunks,
        bc_type: bct,
        umi_type: umit,
    };

    let owriter = Arc::new(Mutex::new(BufWriter::with_capacity(1048576, ofile)));
    let output_cache = Arc::new(DashMap::<u64, libradicl::CorrectedCbChunk>::new());
    let mut allocated_records;
    let mut total_allocated_records = 0;
    let mut last_idx = 0;
    let mut num_output_chunks = 0;

    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
        )
        .progress_chars("╢▌▌░╟");

    let pbar_inner = ProgressBar::new(cc.num_chunks);
    pbar_inner.set_style(sty);
    pbar_inner.tick();

    // create a thread-safe queue based on the number of worker threads
    let n_workers = if num_threads > 1 {
        (num_threads - 1) as usize
    } else {
        1
    };
    let q = Arc::new(ArrayQueue::<(usize, Vec<u8>)>::new(4 * n_workers));

    while last_idx < tsv_map.len() {
        allocated_records = 0;
        output_cache.clear();
        // the number of cells left to process
        let chunks_to_process = Arc::new(AtomicUsize::new(cc.num_chunks as usize));

        // The tsv_map tells us, for each "true" barcode
        // how many records belong to it.  We can scan this information
        // to determine what true barcodes we will keep in memory.
        let init_offset = last_idx;
        for (i, rec) in tsv_map[init_offset..].iter().enumerate() {
            output_cache.insert(
                rec.0,
                libradicl::CorrectedCbChunk::from_label_and_counter(rec.0, rec.1 as u32),
            );
            allocated_records += rec.1;
            last_idx = i + 1;
            if allocated_records >= (max_records as u64) {
                break;
            }
        }

        num_output_chunks += output_cache.len();
        total_allocated_records += allocated_records;
        last_idx += init_offset;

        let mut thread_handles: Vec<thread::JoinHandle<_>> = Vec::with_capacity(n_workers);
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
            let bc_type =
                libradicl::decode_int_type_tag(cc.bc_type).expect("unknown barcode type id.");
            let umi_type =
                libradicl::decode_int_type_tag(cc.umi_type).expect("unknown barcode type id.");
            let owrite = owriter.clone();
            // now, make the worker thread
            let handle = std::thread::spawn(move || {
                // pop from the work queue until everything is
                // processed
                while chunks_remaining.load(Ordering::SeqCst) > 0 {
                    if let Some((_chunk_num, buf)) = in_q.pop() {
                        chunks_remaining.fetch_sub(1, Ordering::SeqCst);
                        let mut nbr = BufReader::new(&buf[..]);
                        libradicl::process_corrected_cb_chunk(
                            &mut nbr,
                            &bc_type,
                            &umi_type,
                            &correct_map,
                            &expected_ori,
                            &oc,
                            &owrite,
                        );
                    }
                }
            });

            thread_handles.push(handle);
        } // for each worker

        // read each chunk
        pbar_inner.reset();
        pbar_inner.set_message(&format!(
            "processing {} / {} total records",
            total_allocated_records, total_to_collate
        ));
        let mut buf = vec![0u8; 65536];
        for cell_num in 0..(cc.num_chunks as usize) {
            let (nbytes_chunk, nrec_chunk) = libradicl::Chunk::read_header(&mut br);
            buf.resize(nbytes_chunk as usize, 0);
            buf.pwrite::<u32>(nbytes_chunk, 0)?;
            buf.pwrite::<u32>(nrec_chunk, 4)?;
            br.read_exact(&mut buf[8..]).unwrap();
            loop {
                if !q.is_full() {
                    let r = q.push((cell_num, buf.clone()));
                    if r.is_ok() {
                        pbar_inner.inc(1);
                        break;
                    }
                }
            }
        }
        pbar_inner.finish();

        for h in thread_handles {
            match h.join() {
                Ok(_) => {}
                Err(_e) => {
                    info!(log, "thread panicked");
                }
            }
        }

        // reset the reader to start of the chunks
        if total_allocated_records < total_to_collate {
            br.get_ref().seek(SeekFrom::Start(pos)).unwrap();
        }
    }

    // make sure we wrote the same number of records that our
    // file suggested we should.
    info!(
        log,
        "\n\ntotal_allocated_records = {}, total_to_collate = {}\n\n",
        total_allocated_records,
        total_to_collate
    );
    assert!(total_allocated_records == total_to_collate);

    pbar_inner.finish_with_message("collated all records.");

    info!(
        log,
        "writing num output chunks ({:?}) to header", num_output_chunks
    );

    owriter.lock().unwrap().flush()?;
    owriter
        .lock()
        .unwrap()
        .get_ref()
        .seek(SeekFrom::Start(
            end_header_pos - (std::mem::size_of::<u64>() as u64),
        ))
        .expect("couldn't seek in output file");
    owriter
        .lock()
        .unwrap()
        .write_all(&num_output_chunks.to_le_bytes())
        .expect("couldn't write to output file.");

    info!(
        log,
        "finished collating input rad file {:?}.",
        i_dir.join("map.rad")
    );
    Ok(())
}
*/

// alternative collate strategy
/*
let handle = std::thread::spawn(move || {
    let mut local_chunks = 0u64;
    let parent = std::path::Path::new(&input_dir);
    // pop from the work queue until everything is
    // processed
    while buckets_remaining.load(Ordering::SeqCst) > 0 {
        if let Some(temp_bucket) = in_q.pop() {
            buckets_remaining.fetch_sub(1, Ordering::SeqCst);
            let new_cap = temp_bucket.0 as usize;
            if cmap.capacity() < new_cap {
                cmap.reserve(temp_bucket.0 as usize);
            } else {
                // cmap.truncate(temp_bucket.0 as usize);
                // cmap.shrink_to(temp_bucket.0 as usize);
                cmap.shrink_to_fit();
            }
            cmap.clear();

            let fname = parent.join(&format!("bucket_{}.tmp", temp_bucket.2.bucket_id));
            // create a new handle for reading
            let tfile = std::fs::File::open(&fname).expect("couldn't open temporary file.");
            let mut treader = BufReader::new(tfile);

            libradicl::collate_temporary_bucket(
                &mut treader,
                &bc_type,
                &umi_type,
                temp_bucket.0,
                temp_bucket.1,
                &mut cmap,
            );

            // we don't need the file or reader anymore
            drop(treader);
            std::fs::remove_file(fname).expect("could not delete temporary file.");

            // go through, add a header to each chunk
            // and flush the chunk to the global output
            // file
            for v in cmap.values_mut() {
                libradicl::dump_chunk(v, &owriter);
            }
            local_chunks += cmap.len() as u64;
            pbar_gather.inc(1);
        }
    }
    local_chunks
});
*/
