use crate::constants as afconst;
use crate::utils as afutils;
use afutils::InternalVersionInfo;
use anyhow::{anyhow, Context};
use crossbeam_queue::ArrayQueue;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use slog::{crit, info};

use libradicl::chunk;
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types;
use libradicl::record::AtacSeqReadRecord;
use libradicl::schema::CollateKey;

use crate::atac::collate::{
    correct_unmapped_counts, get_filter_type, get_most_ambiguous_record, get_num_chunks,
};
use crate::atac::utils as atac_utils;
use itertools::Itertools;
use num_format::{Locale, ToFormattedString};
use scroll::{Pread, Pwrite};
use serde_json::json;
use std::cmp;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufReader, Cursor, Read, Seek, Write};
use std::iter::FromIterator;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct HitInfo {
    pub chr: u32,
    pub start: u32,
    pub frag_len: u16,
    pub barcode: u64,
    pub count: u16,
    // rec_id: u64,
}

impl Ord for HitInfo {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        if self.chr != other.chr {
            self.chr.cmp(&other.chr)
        } else if self.start != other.start {
            self.start.cmp(&other.start)
        } else if self.frag_len != other.frag_len {
            self.frag_len.cmp(&other.frag_len)
        } else {
            self.barcode.cmp(&other.barcode)
        }
    }
}
impl PartialOrd for HitInfo {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

pub fn write_bed_string<W: ?Sized + Write>(
    writer: &mut W,
    hit_info_vec: &[HitInfo],
    ref_names: &[String],
    bc_len: u16,
    rev: bool,
) -> anyhow::Result<()> {
    let bc_len: u8 = bc_len.try_into().unwrap();
    for (count, hinfo) in hit_info_vec.iter().dedup_by_with_count(|x, y| x == y) {
        if hinfo.frag_len < afconst::MAX_ATAC_FRAG_LEN {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                &ref_names[hinfo.chr as usize],
                hinfo.start,
                (hinfo.start + hinfo.frag_len as u32),
                atac_utils::get_bc_string(&hinfo.barcode, rev, bc_len),
                count
            )?;
        }
    }
    Ok(())
}

#[allow(dead_code)]
pub fn get_bed_string(
    hit_info_vec: &[HitInfo],
    ref_names: &[String],
    bc_len: u16,
    rev: bool,
) -> anyhow::Result<String> {
    let mut s = Vec::<u8>::new();
    write_bed_string(&mut s, hit_info_vec, ref_names, bc_len, rev)?;
    Ok(String::from_utf8_lossy(&s).into_owned())
}

#[allow(clippy::too_many_arguments, clippy::manual_clamp)]
pub fn sort_temp_bucket<T: Read + Seek>(
    reader: &mut BufReader<T>,
    bct: &rad_types::RadIntId,
    barcode_len: u16,
    rc: bool,
    ref_names: &[String],
    nrec: u32,
    parent: &Path,
    buck_id: u32,
    compress: bool,
) -> anyhow::Result<()> {
    let mut hit_info_vec: Vec<HitInfo> = Vec::with_capacity(nrec as usize);
    for _ in 0..(nrec as usize) {
        // read the header of the record
        // we don't bother reading the whole thing here
        // because we will just copy later as need be
        let tup = AtacSeqReadRecord::from_bytes_record_header(reader, bct);
        // read the alignment records from the input file
        if tup.1 > 1 {
            continue;
        }
        let rr = AtacSeqReadRecord::from_bytes_with_header(reader, tup.0, tup.1);
        hit_info_vec.push(HitInfo {
            chr: rr.refs[0],
            start: rr.start_pos[0],
            frag_len: rr.frag_lengths[0],
            barcode: rr.bc,
            count: 0,
        })
    }
    hit_info_vec.sort_unstable();

    if compress {
        let bname = parent.join(format!("{}.bed.gz", buck_id));
        afutils::remove_file_if_exists(&bname)?;
        let bd = File::create(&bname)
            .with_context(|| format!("could not create temporary bed file {}", bname.display()))?;
        let mut bd = std::io::BufWriter::new(bd);
        let mut compressor = libdeflater::Compressor::new(libdeflater::CompressionLvl::default());
        let mut bed_vec = Vec::<u8>::new();
        write_bed_string(&mut bed_vec, &hit_info_vec, ref_names, barcode_len, rc)?;
        let buff_sz = compressor.gzip_compress_bound(bed_vec.len());
        let mut compressed_vec = vec![0; buff_sz];
        let compress_res = compressor.gzip_compress(&bed_vec, &mut compressed_vec);
        match compress_res {
            Ok(nbytes) => {
                bd.write_all(&compressed_vec[0..nbytes])?;
            }
            Err(e) => {
                anyhow::bail!(
                    "failed to gzip compress temporary bucket with libdeflater : {:#}",
                    e
                );
            }
        }
    } else {
        let bname = parent.join(format!("{}.bed", buck_id));
        afutils::remove_file_if_exists(&bname)?;
        let bd = File::create(&bname)
            .with_context(|| format!("could not create temporary bed file {}", bname.display()))?;
        let mut bw = std::io::BufWriter::new(bd);
        write_bed_string(&mut bw, &hit_info_vec, ref_names, barcode_len, rc)?;
    }
    Ok(())
    // write_bed(&mut bd, &h_updated, &ref_names, barcode_len, rc);
}

#[allow(clippy::too_many_arguments)]
pub fn sort<P1, P2>(
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

    let rc: bool = mdata["gpl_options"]["rc"].as_bool().unwrap();

    if let Err(es) = calling_version.is_compatible_with(&vd) {
        return Err(anyhow!(es));
    }

    if !parent.join("bin_recs.bin").exists() || !parent.join("bin_lens.bin").exists() {
        crit!(log, "bin file containing records does not exist");
        // std::process::exit(1);
        return Err(anyhow!("execution terminated unexpectedly"));
    }

    let bin_count_file =
        std::fs::File::open(parent.join("bin_recs.bin")).context("couldn't open file")?;
    let bin_rec_counts: Vec<u64> =
        bincode::deserialize_from(bin_count_file).context("couldn't open bin counts file.")?;

    let bin_len_file =
        std::fs::File::open(parent.join("bin_lens.bin")).context("couldn't open file")?;
    let bin_lens: Vec<u64> =
        bincode::deserialize_from(bin_len_file).context("couldn't open bin length file.")?;

    let freq_file =
        std::fs::File::open(parent.join("permit_freq.bin")).context("couldn't open file")?;

    // header buffer
    let mut rbuf = [0u8; 8];

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
    let total_to_collate: u64 = freq_hm.values().sum();
    let tsv_map = Vec::from_iter(freq_hm);
    sort_with_temp(
        input_dir,
        rad_dir,
        num_threads,
        max_records,
        bin_rec_counts,
        bin_lens,
        tsv_map,
        total_to_collate,
        compress_out,
        cmdline,
        version_str,
        rc,
        log,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::manual_clamp)]
pub fn sort_with_temp<P1, P2>(
    input_dir: P1,
    rad_dir: P2,
    num_threads: u32,
    max_records: u32,
    bin_recs: Vec<u64>,
    bin_lens: Vec<u64>,
    tsv_map: Vec<(u64, u64)>,
    _total_to_collate: u64,
    compress_out: bool,
    cmdline: &str,
    version: &str,
    rc: bool,
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
        "sorted bed file {} be compressed",
        if compress_out { "will" } else { "will not" }
    );
    // because :
    // https://superuser.com/questions/865710/write-to-newfile-vs-overwriting-performance-issue

    // writing the collate metadata
    {
        let collate_meta = json!({
            "cmd" : cmdline,
            "version_str" : version,
            "compressed_output" : compress_out,
        });
        let cm_path = parent.join("sort.json");

        let mut cm_file =
            std::fs::File::create(cm_path).context("could not create metadata file.")?;

        let cm_info_string =
            serde_json::to_string_pretty(&collate_meta).context("could not format json.")?;
        cm_file
            .write_all(cm_info_string.as_bytes())
            .context("cannot write to sort.json file")?;
    }

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
    let ref_names = hdr.ref_names.clone();
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

    let binding = file_tag_map?;
    let barcode_tag = binding.get("cblen").expect("tag map must contain cblen");
    let barcode_len: u16 = barcode_tag.try_into()?;

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
    // let mut total_allocated_records = 0;
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
        for (i, nrec) in bin_recs.iter().enumerate() {
            // corrected barcode points to the bucket
            // file.
            moutput_cache.insert(i as u64, temp_buckets.last().unwrap().2.clone());
            allocated_records += nrec;
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
                // total_allocated_records += allocated_records;
                allocated_records = 0;
                num_bucket_chunks = 0;
            }
        }
    }
    if num_bucket_chunks > 0 {
        temp_buckets.last_mut().unwrap().0 = num_bucket_chunks;
        temp_buckets.last_mut().unwrap().1 = allocated_records as u32;
    }
    // total_allocated_records += allocated_records;
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

    // // the number of cells left to process
    let chunks_to_process = Arc::new(AtomicUsize::new(cc.num_chunks as usize));

    let mut thread_handles: Vec<thread::JoinHandle<u64>> = Vec::with_capacity(n_workers);

    let min_rec_len = 24usize; // smallest size an individual record can be loaded in memory
    let max_rec = max_records as usize;
    let num_buckets = temp_buckets.len();
    let num_threads = n_workers;
    let loc_buffer_size = (min_rec_len + (most_ambig_record * 4_usize) - 4_usize).max(
        (1000_usize.max((min_rec_len * max_rec) / (num_buckets * num_threads))).min(262_144_usize),
    ); //131072_usize);

    // // for each worker, spawn off a thread
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
        let loc_bin_lens = bin_lens.clone();

        //let owrite = owriter.clone();
        // now, make the worker thread
        let handle = std::thread::spawn(move || {
            // old code
            // let mut local_buffers = vec![Cursor::new(vec![0u8; loc_buffer_size]); nbuckets];

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

            let size_range = 100000;
            let closure_get_bin_id = |pos: u32, ref_id: usize| {
                atac_utils::get_bin_id(pos, ref_id, size_range, &loc_bin_lens.clone())
            };

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
                        CollateKey::Pos(Box::new(closure_get_bin_id)),
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
    // let pb_msg = format!(
    //     "processing {} / {} total records",
    //     total_allocated_records, total_to_collate
    // );
    // pbar_inner.set_message(pb_msg);

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
    // let mut num_output_chunks = 0u64;
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

    for temp_bucket in temp_buckets.iter() {
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
        // let expected = temp_bucket.1;
        // let observed = temp_bucket.2.num_records_written.load(Ordering::SeqCst);
        // assert_eq!(expected, observed);

        // let md = std::fs::metadata(parent.join(format!("bucket_{}.tmp", i)))?;
        // let expected_bytes = temp_bucket.2.num_bytes_written.load(Ordering::SeqCst);
        // let observed_bytes = md.len();
        // assert_eq!(expected_bytes, observed_bytes);
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

    // // for each worker, spawn off a thread
    for _worker in 0..n_workers {
        // each thread will need to access the work queue
        let in_q = fq.clone();
        // let barcode_len = barcode_len.clone();
        // the output cache and correction map
        // let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
        // let mut cmap = HashMap::<u64, TempCellInfo, ahash::RandomState>::with_hasher(s);
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

        let r_names = ref_names.clone();
        // and the progress bar
        let pbar_gather = pbar_gather.clone();

        // now, make the worker threads
        let handle = std::thread::spawn(move || {
            let local_chunks = 0u64;
            let parent = std::path::Path::new(&input_dir);
            // pop from the work queue until everything is
            // processed
            while buckets_remaining.load(Ordering::SeqCst) > 0 {
                if let Some(temp_bucket) = in_q.pop() {
                    buckets_remaining.fetch_sub(1, Ordering::SeqCst);
                    // cmap.clear();

                    let fname = parent.join(format!("bucket_{}.tmp", temp_bucket.2.bucket_id));
                    // create a new handle for reading
                    let tfile = std::fs::File::open(&fname).expect("couldn't open temporary file.");
                    let mut treader = BufReader::new(tfile);

                    let _ = sort_temp_bucket(
                        &mut treader,
                        &bc_type,
                        barcode_len,
                        rc,
                        &r_names,
                        temp_bucket.2.num_records_written.load(Ordering::SeqCst),
                        parent,
                        temp_bucket.2.bucket_id,
                        compress_out,
                    );

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
    for temp_bucket in &temp_buckets {
        let mut bclone = temp_bucket.clone();
        // keep trying until we can push this payload
        while let Err(t) = fq.push(bclone) {
            bclone = t;
            // no point trying to push if the queue is full
            while fq.is_full() {}
        }
    }

    // wait for all of the workers to finish
    // let mut num_output_chunks = 0u64;
    for h in thread_handles.drain(0..) {
        match h.join() {
            Ok(_c) => {
                // num_output_chunks += c;
            }
            Err(_e) => {
                info!(log, "thread panicked");
            }
        }
    }
    pbar_gather.finish_with_message("gathered all temp files.");

    // make sure we wrote the same number of records that our
    // file suggested we should.
    // assert_eq!(total_allocated_records, total_to_collate);

    // info!(
    //     log,
    //     "writing num output chunks ({}) to header",
    //     num_output_chunks.to_formatted_string(&Locale::en)
    // );

    info!(
        log,
        "expected number of output chunks {}",
        expected_output_chunks.to_formatted_string(&Locale::en)
    );

    let bed_sfx = if compress_out { ".bed.gz" } else { ".bed" };

    let bedname = parent.join(format!("map{}", bed_sfx));
    /*
    if compress_out {
        let bedname = parent.join("map.bed.gz");
        afutils::remove_file_if_exists(&bedname)?;
        let out_bed_file = File::create(&bedname).with_context(|| {
            format!(
                "could not create target output bed file {}",
                bedname.display()
            )
        })?;
        let mut encoder = GzEncoder::new(out_bed_file, Compression::default());
        for i in 0..temp_buckets.len() {
            let temp_bed_name = parent.join(format!("{}.bed.gz", i));
            let mut input = File::open(&temp_bed_name)?;
            io::copy(&mut input, &mut encoder)?;
            std::fs::remove_file(&temp_bed_name)?;
        }
        encoder.finish()?;
    } else {
        let bedname = parent.join("map.bed");
    */
    afutils::remove_file_if_exists(&bedname)?;
    let mut out_bed_file = File::create(&bedname).with_context(|| {
        format!(
            "could not create target output bed file {}",
            bedname.display()
        )
    })?;
    for i in 0..temp_buckets.len() {
        let temp_bed_name = parent.join(format!("{}{}", i, bed_sfx));
        let mut input = File::open(&temp_bed_name)?;
        io::copy(&mut input, &mut out_bed_file)?;
        std::fs::remove_file(&temp_bed_name)?;
        // }
    }
    // assert_eq!(
    //     expected_output_chunks,
    //     num_output_chunks,
    //     "expected to write {} chunks but wrote {}",
    //     expected_output_chunks.to_formatted_string(&Locale::en),
    //     num_output_chunks.to_formatted_string(&Locale::en),
    // );

    info!(log, "merging temp files");
    Ok(())
}
