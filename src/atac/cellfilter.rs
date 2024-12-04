use crate::atac::prog_opts::GenPermitListOpts;
use crate::utils as afutils;
use anyhow::{anyhow, Context};
use bstr::io::BufReadExt;
// use indexmap::map::IndexMap;
use itertools::Itertools;
use libradicl::exit_codes;
use libradicl::rad_types;
use libradicl::rad_types::TagValue;
use libradicl::utils::has_data_left;
use libradicl::BarcodeLookupMap;
use libradicl::{
    chunk,
    header::RadPrelude,
    record::{AtacSeqReadRecord, AtacSeqRecordContext},
};
use num_format::{Locale, ToFormattedString};
use serde::Serialize;
use serde_json::json;
use slog::crit;
use slog::info;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

/// Initialize the index map with key being references and position
/// Take the largest reference length (from chromosome)
/// For each chromosome divide into ranges of size_range uptil max_size
/// If we had ref length corresponding to each chromosome, we would divide the reference based on its length only
pub fn initialize_rec_list(
    b_lens: &mut [u64],
    ref_lens: &[u64],
    size_range: u64,
) -> anyhow::Result<u64> {
    let mut tot_zeros = 0;
    let num_refs = ref_lens.len();

    for r in 0..num_refs {
        let nrange = (ref_lens[r] as f32 / size_range as f32).ceil() as u64;
        let cum = b_lens[r] + nrange;
        b_lens[r + 1] = cum;
        tot_zeros += nrange;
    }
    Ok(tot_zeros)
}

#[derive(Clone, Debug, Serialize)]
pub enum CellFilterMethod {
    // use the distance method to
    // automatically find the knee
    // in the curve
    KneeFinding,
    // barcodes will be provided in the
    // form of an *unfiltered* external
    // permit list
    UnfilteredExternalList(PathBuf, usize),
}

pub fn update_barcode_hist_unfiltered(
    hist: &mut HashMap<u64, u64, ahash::RandomState>,
    unmatched_bc: &mut Vec<u64>,
    max_ambiguity_read: &mut usize,
    chunk: &chunk::Chunk<AtacSeqReadRecord>,
    bins: &mut [u64],
    blens: &[u64],
    size_range: u64,
) -> usize {
    let mut num_strand_compat_reads = 0usize;
    for r in &chunk.reads {
        num_strand_compat_reads += 1;
        *max_ambiguity_read = r.refs.len().max(*max_ambiguity_read);
        // lookup the barcode in the map of unfiltered known
        // barcodes
        match hist.get_mut(&r.bc) {
            // if we find a match, increment the count
            Some(c) => *c += 1,
            // otherwise, push this into the unmatched list
            None => {
                unmatched_bc.push(r.bc);
            }
        }

        // for (i,_j) in r.start_pos.iter().enumerate() {
        if r.start_pos.len() == 1 {
            let i = 0;
            let ref_id = r.refs[i];
            let sp = r.start_pos[i];
            let bid = sp as u64 / size_range;
            let ind: usize = (blens[ref_id as usize] + bid) as usize;
            bins[ind] += 1;
        }
    }
    num_strand_compat_reads
}

fn populate_unfiltered_barcode_map<T: Read>(
    br: BufReader<T>,
    first_bclen: &mut usize,
    rev_bc: bool,
) -> HashMap<u64, u64, ahash::RandomState> {
    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut hm = HashMap::with_hasher(s);

    // read through the external unfiltered barcode list
    // and generate a vector of encoded barcodes
    // let mut kv = Vec::<u64>::new();
    for l in br.byte_lines().flatten() {
        if *first_bclen == 0 {
            *first_bclen = l.len();
        } else {
            assert_eq!(
                *first_bclen,
                l.len(),
                "found barcodes of different lengths {} and {}",
                *first_bclen,
                l.len()
            );
        }
        if let Some((_, km, _)) =
            needletail::bitkmer::BitNuclKmer::new(&l[..], l.len() as u8, false).next()
        {
            if rev_bc {
                let km_rev = needletail::bitkmer::reverse_complement(km);
                hm.insert(km_rev.0, 0);
            } else {
                hm.insert(km.0, 0);
            }
        }
    }
    hm
}

#[allow(clippy::unnecessary_unwrap, clippy::too_many_arguments)]
fn process_unfiltered(
    mut hm: HashMap<u64, u64, ahash::RandomState>,
    mut unmatched_bc: Vec<u64>,
    file_tag_map: &rad_types::TagMap,
    filter_meth: &CellFilterMethod,
    output_dir: &PathBuf,
    version: &str,
    max_ambiguity_read: usize,
    num_chunks: u32,
    cmdline: &str,
    log: &slog::Logger,
    bmax: u64,
    gpl_opts: &GenPermitListOpts,
) -> anyhow::Result<u64> {
    // the smallest number of reads we'll allow per barcode
    let min_freq = match filter_meth {
        CellFilterMethod::UnfilteredExternalList(_, min_reads) => {
            info!(log, "minimum num reads for barcode pass = {}", *min_reads);
            *min_reads as u64
        }
        _ => {
            unimplemented!();
        }
    };

    // the set of barcodes we'll keep
    let mut kept_bc = Vec::<u64>::new();

    // iterate over the count map
    for (k, v) in hm.iter_mut() {
        // if this satisfies our requirement for the minimum count
        // then keep this barcode
        if *v >= min_freq {
            kept_bc.push(*k);
        } else {
            // otherwise, we have to add this barcode's
            // counts to our unmatched list
            for _ in 0..*v {
                unmatched_bc.push(*k);
            }
            // and then reset the counter for this barcode to 0
            *v = 0u64;
        }
    }

    // drop the absent barcodes from hm
    hm.retain(|_, &mut v| v > 0);

    // how many we will keep
    let num_passing = kept_bc.len();
    info!(
        log,
        "num_passing = {}",
        num_passing.to_formatted_string(&Locale::en)
    );
    let barcode_tag = file_tag_map
        .get("cblen")
        .expect("tag map must contain cblen");
    let barcode_len: u16 = barcode_tag.try_into()?;

    // let ref_lens = file_tag_map
    //     .get("ref_lengths")
    //     .expect("tag map must contain ref_lengths");
    // now, we create a second barcode map with just the barcodes
    // for cells we will keep / rescue.
    let bcmap2 = BarcodeLookupMap::new(kept_bc, barcode_len as u32);
    info!(
        log,
        "found {} cells with non-trivial number of reads by exact barcode match",
        bcmap2.barcodes.len().to_formatted_string(&Locale::en)
    );

    // finally, we'll go through the set of unmatched barcodes
    // and try to rescue those that have a *unique* neighbor in the
    // list of retained barcodes.

    //let mut found_exact = 0usize;
    let mut found_approx = 0usize;
    let mut ambig_approx = 0usize;
    let mut not_found = 0usize;

    let start_unmatched_time = Instant::now();

    unmatched_bc.sort_unstable();

    let mut distinct_unmatched_bc = 0usize;
    let mut distinct_recoverable_bc = 0usize;

    // mapping the uncorrected barcode to what it corrects to
    let mut corrected_list = Vec::<(u64, u64)>::with_capacity(1_000_000);

    for (count, ubc) in unmatched_bc.iter().dedup_with_count() {
        // try to find the unmatched barcode, but
        // look up to 1 edit away
        match bcmap2.find_neighbors(*ubc, false) {
            // if we have a match
            (Some(x), n) => {
                let cbc = bcmap2.barcodes[x];
                // if the uncorrected barcode had a
                // single, unique retained neighbor
                if cbc != *ubc && n == 1 {
                    // then increment the count of this
                    // barcode by 1 (because we'll correct to it)
                    if let Some(c) = hm.get_mut(&cbc) {
                        *c += count as u64;
                        corrected_list.push((*ubc, cbc));
                    }
                    // this counts as an approximate find
                    found_approx += count;
                    distinct_recoverable_bc += 1;
                }
                // if we had > 1 single-mismatch neighbor
                // then don't keep the barcode, but remember
                // the count of such events
                if n > 1 {
                    ambig_approx += count;
                }
            }
            // if we had no single-mismatch neighbor
            // then this barcode is not_found and gets
            // dropped.
            (None, _) => {
                not_found += count;
            }
        }
        distinct_unmatched_bc += 1;
    }
    let unmatched_duration = start_unmatched_time.elapsed();
    let num_corrected = distinct_recoverable_bc as u64;

    info!(
        log,
        "There were {} distinct unmatched barcodes, and {} that can be recovered",
        distinct_unmatched_bc,
        distinct_recoverable_bc
    );
    info!(
        log,
        "Matching unmatched barcodes to retained barcodes took {:?}", unmatched_duration
    );
    info!(log, "Of the unmatched barcodes\n============");
    info!(
        log,
        "\t{} had exactly 1 single-edit neighbor in the retained list",
        found_approx.to_formatted_string(&Locale::en)
    );
    info!(
        log,
        "\t{} had >1 single-edit neighbor in the retained list",
        ambig_approx.to_formatted_string(&Locale::en)
    );
    info!(
        log,
        "\t{} had no neighbor in the retained list",
        not_found.to_formatted_string(&Locale::en)
    );

    let parent = std::path::Path::new(output_dir);
    let o_path = parent.join("permit_freq.bin");

    match afutils::write_permit_list_freq(&o_path, barcode_len, &hm) {
        Ok(_) => {}
        Err(error) => {
            panic!("Error: {}", error);
        }
    };

    /*
    // don't need this right now
    let s_path = parent.join("bcmap.bin");
    let s_file = std::fs::File::create(&s_path).expect("could not create serialization file.");
    let mut s_writer = BufWriter::new(&s_file);
    bincode::serialize_into(&mut s_writer, &bcmap2).expect("couldn't serialize barcode list.");
    */

    // now that we are done with using hm to count, we can repurpose it as
    // the correction map.
    for (k, v) in hm.iter_mut() {
        // each present barcode corrects to itself
        *v = *k;
    }
    for (uncorrected, corrected) in corrected_list.iter() {
        hm.insert(*uncorrected, *corrected);
    }

    let pm_path = parent.join("permit_map.bin");
    let pm_file = std::fs::File::create(pm_path).context("could not create serialization file.")?;
    let mut pm_writer = BufWriter::new(&pm_file);
    bincode::serialize_into(&mut pm_writer, &hm)
        .context("couldn't serialize permit list mapping.")?;

    let meta_info = json!({
    "version_str" : version,
    "max-ambig-record" : max_ambiguity_read,
    "num-chunks" : num_chunks,
    "cmd" : cmdline,
    "permit-list-type" : "unfiltered",
    "gpl_options" : &gpl_opts,
    "max-rec-in-bin": bmax
    });

    let m_path = parent.join("generate_permit_list.json");
    let mut m_file = std::fs::File::create(m_path).context("could not create metadata file.")?;

    let meta_info_string =
        serde_json::to_string_pretty(&meta_info).context("could not format json.")?;
    m_file
        .write_all(meta_info_string.as_bytes())
        .context("cannot write to generate_permit_list.json file")?;

    info!(
        log,
        "total number of distinct corrected barcodes : {}",
        num_corrected.to_formatted_string(&Locale::en)
    );

    Ok(num_corrected)
}

pub fn generate_permit_list(gpl_opts: GenPermitListOpts) -> anyhow::Result<u64> {
    let rad_dir = gpl_opts.input_dir;
    let output_dir = gpl_opts.output_dir;

    let filter_meth = gpl_opts.fmeth.clone();
    let version = gpl_opts.version;
    let cmdline = gpl_opts.cmdline;
    let log = gpl_opts.log;
    let rc = gpl_opts.rc;
    let mut num_chunks = 0;
    let size_range: u64 = 100000;
    let i_dir = std::path::Path::new(&rad_dir);
    let o_dir_path = std::path::Path::new(&output_dir);
    // let ref_lens = read_ref_lengths(i_dir)?;

    // should we assume this condition was already checked
    // during parsing?
    if !i_dir.exists() {
        crit!(
            log,
            "the input rad path {} does not exist",
            rad_dir.display()
        );
        // std::process::exit(1);
        return Err(anyhow!("execution terminated unexpectedly"));
    }

    let mut first_bclen = 0usize;
    let mut unfiltered_bc_counts = None;

    let p_exi = o_dir_path
        .try_exists()
        .expect("Can't check existence of path");
    if !p_exi {
        std::fs::create_dir_all(o_dir_path).with_context(|| {
            format!(
                "couldn't create path to output directory {}",
                o_dir_path.display()
            )
        })?;
    }

    if let CellFilterMethod::UnfilteredExternalList(fname, _) = &filter_meth {
        let i_file = File::open(fname).context("could not open input file")?;
        let br = BufReader::new(i_file);
        unfiltered_bc_counts = Some(populate_unfiltered_barcode_map(br, &mut first_bclen, rc));
        info!(
            log,
            "number of unfiltered bcs read = {}",
            unfiltered_bc_counts
                .as_ref()
                .unwrap()
                .len()
                .to_formatted_string(&Locale::en)
        );
    }

    let i_file = File::open(i_dir.join("map.rad")).context("could not open input rad file")?;
    let mut br = BufReader::new(i_file);
    let prelude = RadPrelude::from_bytes(&mut br)?;
    let hdr = &prelude.hdr;

    info!(
        log,
        "paired : {:?}, ref_count : {}, num_chunks : {}",
        hdr.is_paired != 0,
        hdr.ref_count.to_formatted_string(&Locale::en),
        hdr.num_chunks.to_formatted_string(&Locale::en),
    );

    // file-level
    let fl_tags = &prelude.file_tags;
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = &prelude.read_tags;
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());

    const BNAME: &str = "barcode";
    // let mut bct: Option<RadType> = None;

    for rt in &rl_tags.tags {
        // if this is one of our tags
        if rt.name == BNAME && !rt.typeid.is_int_type() {
            crit!(
                log,
                "currently only RAD types 1--4 are supported for 'b' tags."
            );
            std::process::exit(exit_codes::EXIT_UNSUPPORTED_TAG_TYPE);
            // if rt.name == BNAME {
            //     bct = Some(rt.typeid);
            // }
        }
    }

    // alignment-level
    let al_tags = &prelude.aln_tags;
    info!(log, "read {:?} alignment-level tags", al_tags.tags.len());

    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    info!(log, "File-level tag values {:?}", file_tag_map);

    let ref_tag = file_tag_map
        .get("ref_lengths")
        .expect("tag must contain ref lengths");
    let TagValue::ArrayU64(ref_lens) = ref_tag else {
        todo!()
    };

    let record_context = prelude.get_record_context::<AtacSeqRecordContext>()?;
    let mut num_reads: usize = 0;

    let mut blens: Vec<u64> = vec![0; ref_lens.len() + 1];
    let tot_bins = initialize_rec_list(&mut blens, ref_lens, size_range);
    let mut bins: Vec<u64> = vec![0; tot_bins.unwrap() as usize];

    // if dealing with the unfiltered type
    // the set of barcodes that are not an exact match for any known barcodes
    let mut unmatched_bc: Vec<u64>;
    // let mut num_orientation_compat_reads = 0usize;
    let mut max_ambiguity_read = 0usize;
    let mut num_orientation_compat_reads = 0usize;
    // Tracking if a unique or a multihit

    match filter_meth {
        CellFilterMethod::UnfilteredExternalList(_, _min_reads) => {
            unmatched_bc = Vec::with_capacity(10000000);
            // the unfiltered_bc_count map must be valid in this branch

            if let Some(mut hmu) = unfiltered_bc_counts {
                while has_data_left(&mut br).expect("encountered error reading input file") {
                    let c = chunk::Chunk::<AtacSeqReadRecord>::from_bytes(&mut br, &record_context);
                    num_orientation_compat_reads += update_barcode_hist_unfiltered(
                        &mut hmu,
                        &mut unmatched_bc,
                        &mut max_ambiguity_read,
                        &c,
                        &mut bins,
                        &blens,
                        size_range,
                    );
                    num_chunks += 1;
                    num_reads += c.reads.len();
                }

                let bin_recs_path = output_dir.join("bin_recs.bin");
                let br_file = std::fs::File::create(bin_recs_path)
                    .expect("could not create serialization file.");
                let mut br_writer = BufWriter::new(&br_file);
                bincode::serialize_into(&mut br_writer, &bins)
                    .expect("couldn't serialize bins recs.");

                let bin_lens_path = output_dir.join("bin_lens.bin");
                let bl_file = std::fs::File::create(bin_lens_path)
                    .expect("could not create serialization file.");
                let mut bl_writer = BufWriter::new(&bl_file);
                bincode::serialize_into(&mut bl_writer, &blens)
                    .expect("couldn't serialize bins lengths.");

                info!(
                    log,
                    "observed {} reads ({} orientation consistent) in {} chunks --- max ambiguity read occurs in {} refs",
                    num_reads.to_formatted_string(&Locale::en),
                    num_orientation_compat_reads.to_formatted_string(&Locale::en),
                    num_chunks.to_formatted_string(&Locale::en),
                    max_ambiguity_read.to_formatted_string(&Locale::en)
                );
                process_unfiltered(
                    hmu,
                    unmatched_bc,
                    &file_tag_map,
                    &filter_meth,
                    output_dir,
                    version,
                    max_ambiguity_read,
                    num_chunks,
                    cmdline,
                    log,
                    *bins.iter().max().unwrap(),
                    &gpl_opts,
                )
            } else {
                Ok(0)
            }
        }
        _ => Ok(0),
    }
}
