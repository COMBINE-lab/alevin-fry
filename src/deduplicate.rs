use crate::prog_opts::DeduplicateOpts;
use anyhow::Context;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use libradicl::{
    readers::ParallelRadReader,
    record::AtacSeqReadRecord,
};
use num_format::ToFormattedString;
use slog::info;
use std::cmp;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;
use std::num::NonZeroUsize;
use std::sync::{
    atomic:: Ordering,
    Arc, Mutex,
    atomic::AtomicU32
};
use std::thread;
pub type MetaChunk = (usize, usize, u32, u32, Vec<u8>);
use crate::utils as af_utils;
use itertools::Itertools;
use num_format::Locale;

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct HitInfo {
    chr: u32,
    start: u32,
    frag_len: u16,
    barcode: u64,
    count: u16,
    // rec_id: u64,
}

impl Ord for HitInfo {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        if self.chr != other.chr {
            self.chr.cmp(&other.chr)
        } else if self.start != other.start {
            self.start.cmp(&other.start)
        } else {
            self.frag_len.cmp(&other.frag_len)
        }
    }
}
impl PartialOrd for HitInfo {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

pub fn write_bed(
    bd_writer_lock: &Arc<Mutex<File>>,
    hit_info_vec: &[HitInfo],
    ref_names: &[String],
    rev: bool,
    bc_len: u8,
    count_frag: &Arc<AtomicU32>
) {
    let mut s = "".to_string();
    for i in 0..hit_info_vec.len() {
        if hit_info_vec[i].frag_len < 2000 {
            let s2 = [
                ref_names[hit_info_vec[i].chr as usize].clone(),
                hit_info_vec[i].start.to_string(),
                (hit_info_vec[i].start + hit_info_vec[i].frag_len as u32).to_string(),
                af_utils::get_bc_string(&hit_info_vec[i].barcode, rev, bc_len),
                hit_info_vec[i].count.to_string(),
            ]
            .join("\t");
            s.push_str(&s2);
            s.push('\n');
        }
        else {
            count_frag.fetch_add(1, Ordering::SeqCst);
        }
    }

    let bd_lock = bd_writer_lock.lock();
    let writer = &mut *bd_lock.unwrap();
    writer.write_all(s.as_bytes()).unwrap();
}

pub fn deduplicate(dedup_opts: DeduplicateOpts) -> anyhow::Result<()> {
    let parent = std::path::Path::new(dedup_opts.input_dir);
    let log = dedup_opts.log;
    let collate_md_file =
        File::open(parent.join("collate.json")).context("could not open the collate.json file.")?;
    let collate_md: serde_json::Value = serde_json::from_reader(&collate_md_file)?;

    // is the collated RAD file compressed?
    let compressed_input = collate_md["compressed_output"]
        .as_bool()
        .context("could not read compressed_output field from collate metadata.")?;

    if compressed_input {
        let i_file =
            File::open(parent.join("map.collated.rad.sz")).context("run collate before quant")?;
        let br = snap::read::FrameDecoder::new(BufReader::new(&i_file));

        info!(
            log,
            "quantifying from compressed, collated RAD file {:?}", i_file
        );
        Ok(())
        // do_deduplicate(br, dedup_opts)
    } else {
        let i_file =
            File::open(parent.join("map.collated.rad")).context("run collate before quant")?;
        let metadata = i_file.metadata()?;
        let file_len = metadata.len();
        let br = BufReader::new(i_file);

        info!(
            log,
            "quantifying from uncompressed, collated RAD file {:?}",
            parent.join("map.collated.rad")
        );

        do_deduplicate(br, dedup_opts, file_len)
    }
}

pub fn do_deduplicate(
    br: BufReader<File>,
    dedup_opts: DeduplicateOpts,
    file_len: u64,
) -> anyhow::Result<()> {
    let num_threads = dedup_opts.num_threads;

    let n_workers = if num_threads > 1 {
        (num_threads - 1) as usize
    } else {
        1
    };
    
    let mut rad_reader = ParallelRadReader::<AtacSeqReadRecord, BufReader<File>>::new(
        br,
        NonZeroUsize::new(n_workers).unwrap(),
    );
    let refs = &rad_reader.prelude.hdr.ref_names;
    let log = dedup_opts.log;
    let prelude = &rad_reader.prelude;
    let hdr = &prelude.hdr;

    let num_multimappings = Arc::new(AtomicU32::new(0 as u32));
    let num_dedup = Arc::new(AtomicU32::new(0 as u32));
    let num_frag_counts = Arc::new(AtomicU32::new(0 as u32)); // fragments larger than 2000
    let num_non_mapped_pair = Arc::new(AtomicU32::new(0 as u32)); // fragments larger than 2000

    if let Ok(summary) = prelude.summary(None) {
        println!("{}", summary);
    }
    info!(
        log,
        "paired : {:?}, ref_count : {}, num_chunks : {}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en)
    );

    let parent = std::path::Path::new(dedup_opts.input_dir);
    let bed_path = parent.join("map.bed");
    let num_chunks = hdr.num_chunks;

    let bc_unmapped_file = File::open(parent.join("unmapped_bc_count_collated.bin")).unwrap();
    let bc_unmapped_map: Arc<HashMap<u64, u32>> =
        Arc::new(bincode::deserialize_from(&bc_unmapped_file).unwrap());

    let fl_tags = &prelude.file_tags;
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = &prelude.read_tags;
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = &prelude.aln_tags;
    info!(log, "read {:?} alignment-level tags", al_tags.tags.len());

    let file_tag_map = &rad_reader.file_tag_map;
    info!(log, "File-level tag values {:?}", file_tag_map);

    let barcode_tag = file_tag_map
        .get("cblen")
        .expect("tag map must contain cblen");
    let barcode_len: u16 = barcode_tag.try_into()?;

    let pbar = ProgressBar::with_draw_target(
        Some(num_chunks),
        ProgressDrawTarget::stderr_with_hz(5u8), // update at most 5 times/sec.
    );
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .expect("ProgressStyle template was invalid.")
            .progress_chars("╢▌▌░╟"),
    );

    let bed_writer = Arc::new(Mutex::new(File::create(bed_path).unwrap()));
    let mut thread_handles: Vec<thread::JoinHandle<usize>> = Vec::with_capacity(n_workers);

    for _worker in 0..n_workers {
        let rd = rad_reader.is_done();
        let q = rad_reader.get_queue();
        let bd = bed_writer.clone();
        let refs = refs.clone();
        let num_multimappings = num_multimappings.clone();
        let num_dedup = num_dedup.clone();
        let num_frag_counts = num_frag_counts.clone();
        let num_non_mapped_pair = num_non_mapped_pair.clone();

        let unmapped_count = bc_unmapped_map.clone();

        let handle = std::thread::spawn(move || {
            let mut nrec_processed = 0_usize;
            while !rd.load(Ordering::SeqCst) {
                while let Some(meta_chunk) = q.pop() {
                    for c in meta_chunk.iter() {
                        nrec_processed += c.nrec as usize;
                        let mut hit_info_vec: Vec<HitInfo> = Vec::with_capacity(c.nrec as usize);
                        // println!("Chunk :: nbytes: {}, nrecs: {}", c.nbytes, c.nrec);
                        assert_eq!(c.nrec as usize, c.reads.len());
                        for r in c.reads.iter() {
                            let na = r.refs.len();
                            // add a field tracking such counts
                            if na == 1 && r.map_type[0] == 4 {
                                hit_info_vec.push(HitInfo {
                                    chr: r.refs[0],
                                    start: r.start_pos[0],
                                    frag_len: r.frag_lengths[0],
                                    barcode: r.bc,
                                    count: 0,
                                })
                            }
                            else if na > 1 {
                                num_multimappings.fetch_add(1, Ordering::SeqCst);
                                continue;
                            }
                            else {
                                num_non_mapped_pair.fetch_add(1, Ordering::SeqCst);
                            }
                        }
                        hit_info_vec.sort_unstable();
                        let mut h_updated: Vec<HitInfo> = Vec::with_capacity(hit_info_vec.len());
                        for (count, hv) in hit_info_vec.iter_mut().dedup_with_count() {
                            hv.count = count as u16;
                            h_updated.push(*hv);
                            if count > 1 {
                                num_dedup.fetch_add(1, Ordering::SeqCst);
                            }
                        }
                        drop(hit_info_vec);
                        write_bed(&bd, &h_updated, &refs, dedup_opts.rev, barcode_len as u8, &num_frag_counts);
                    }
                }
            }
            nrec_processed
        });
        thread_handles.push(handle);
    }
    let header_offset = rad_reader.get_byte_offset();
    let pbar = ProgressBar::new(file_len - header_offset);
    pbar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-"),
    );
    pbar.set_draw_target(ProgressDrawTarget::stderr_with_hz(5));
    let cb = |new_bytes: u64, _new_rec: u64| {
        pbar.inc(new_bytes);
    };
    let _ = rad_reader.start_chunk_parsing(Some(cb)); //libradicl::readers::EMPTY_METACHUNK_CALLBACK);
    let mut total_processed = 0;
    for handle in thread_handles {
        total_processed += handle.join().expect("The parsing thread panicked");
    }
    pbar.finish_with_message(format!(
        "finished parsing RAD file; processed {} total records\n",
        total_processed
    ));
    info!(
        log,
        "Number of records with greater than 1 mapping {}", num_multimappings.load(Ordering::SeqCst)
        );
    info!(
        log,
        "Number of records that are deduplicated {}", num_dedup.load(Ordering::SeqCst)
        );
    info!(
        log,
        "Number of records that are not mapped pairs {}", num_non_mapped_pair.load(Ordering::SeqCst)
        );
    info!(
        log,
        "Number of records that have frag length > 2000 {}", num_frag_counts.load(Ordering::SeqCst)
        );
    Ok(())
}
