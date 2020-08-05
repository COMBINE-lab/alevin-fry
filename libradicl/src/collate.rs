// Copyright 2020 Rob Patro, Avi Srivastava. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

extern crate bincode;
extern crate indicatif;
extern crate serde;
extern crate slog;

use self::indicatif::{ProgressBar, ProgressStyle};
use self::slog::info;
use bio_types::strand::Strand;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};

use crate as libradicl;

pub fn collate(
    input_dir: String,
    rad_file: String,
    max_records: u32,
    expected_ori: Strand,
    log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {
    let parent = std::path::Path::new(&input_dir);
    let mut ofile = File::create(parent.join("map.collated.rad")).unwrap();
    let i_file = File::open(&rad_file).unwrap();

    let mut br = BufReader::new(i_file);

    let hdr = libradicl::RADHeader::from_bytes(&mut br);

    let end_header_pos =
        br.get_ref().seek(SeekFrom::Current(0)).unwrap() - (br.buffer().len() as u64);

    info!(
        log,
        "paired : {:?}, ref_count : {:?}, num_chunks : {:?}",
        hdr.is_paired != 0,
        hdr.ref_count,
        hdr.num_chunks
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

    let mut owriter = BufWriter::new(ofile);

    type TSVRec = (u64, u64);
    let mut tsv_map = Vec::<(u64, u64)>::new(); //HashMap::<u64, u64>::new();

    let freq_file =
        std::fs::File::open(parent.join("permit_freq.tsv")).expect("couldn't open file");
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(freq_file);

    let mut total_to_collate = 0;
    for result in rdr.deserialize() {
        let record: TSVRec = result?;
        tsv_map.push(record);
        total_to_collate += record.1;
    }

    // sort this so that we deal with largest cells (by # of reads) first
    // sort in _descending_ order by count.
    quickersort::sort_by_key(&mut tsv_map[..], |&a: &(u64, u64)| std::cmp::Reverse(a.1));
    //println!("{:?}", tsv_map);

    // get the correction map
    let cmfile = std::fs::File::open(parent.join("permit_map.bin")).unwrap();
    let correct_map: HashMap<u64, u64> = bincode::deserialize_from(&cmfile).unwrap();

    info!(
        log,
        "deserialized correction map of length : {:?}",
        correct_map.len()
    );

    let mut pass_num = 1;
    let cc = libradicl::ChunkConfig {
        num_chunks: hdr.num_chunks,
        bc_type: bct,
        umi_type: umit,
    };

    let mut output_cache = HashMap::<u64, libradicl::CorrectedCBChunk>::new();
    let mut allocated_records;
    let mut total_allocated_records = 0;
    let mut last_idx = 0;
    let mut num_output_chunks = 0;

    let pbar = ProgressBar::new(total_to_collate);
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );

    while last_idx < tsv_map.len() {
        allocated_records = 0;
        output_cache.clear();

        // The tsv_map tells us, for each "true" barcode
        // how many records belong to it.  We can scan this information
        // to determine what true barcodes we will keep in memory.
        let init_offset = last_idx;
        for (i, rec) in tsv_map[init_offset..].iter().enumerate() {
            output_cache.insert(
                rec.0,
                libradicl::CorrectedCBChunk::from_label_and_counter(rec.0, rec.1),
            );
            allocated_records += rec.1;
            last_idx = i + 1;
            if allocated_records >= (max_records as u64) {
                //info!(log, "pass {:?} will collate {:?} records", pass_num, allocated_records);
                break;
            }
        }

        num_output_chunks += output_cache.len();
        total_allocated_records += allocated_records;
        last_idx += init_offset;

        // collect the output for the current barcode set
        libradicl::collect_records(&mut br, &cc, &correct_map, &expected_ori, &mut output_cache);

        // dump the output we have
        libradicl::dump_output_cache(&mut owriter, &output_cache, &cc);

        // reset the reader to start of the chunks
        br.get_ref().seek(SeekFrom::Start(pos)).unwrap();
        pass_num += 1;
        //info!(log, "total collated {:?} / {:?}", total_allocated_records, total_to_collate);
        //info!(log, "last index processed {:?} / {:?}", last_idx, tsv_map.len());
        pbar.inc(allocated_records);
    }

    // make sure we wrote the same number of records that our
    // file suggested we should.
    assert!(total_allocated_records == total_to_collate);

    info!(
        log,
        "writing num output chunks ({:?}) to header", num_output_chunks
    );
    owriter.flush()?;
    owriter
        .get_ref()
        .seek(SeekFrom::Start(
            end_header_pos - (std::mem::size_of::<u64>() as u64),
        ))
        .expect("couldn't seek in output file");
    owriter
        .write_all(&num_output_chunks.to_le_bytes())
        .expect("couldn't write to output file.");
    let pb_msg = format!("finished collating in {} passes", pass_num);
    pbar.finish_with_message(&pb_msg);

    info!(log, "finished collating input rad file {:?}.", rad_file);
    Ok(())
}
