extern crate indicatif;
extern crate needletail;
extern crate rust_htslib;
extern crate scroll;
extern crate slog;

use self::indicatif::{ProgressBar, ProgressStyle};
use self::slog::info;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Cursor, Seek, SeekFrom, Write};
// use std::sync::{Arc, Mutex};
// use needletail::bitkmer::*;
use rand::Rng;
use rust_htslib::bam::HeaderView;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use std::str;

use crate as libradicl;

#[allow(dead_code)]
fn get_random_nucl() -> &'static str {
    let nucl = vec!["A", "T", "G", "C"];
    let mut rng = rand::thread_rng();
    let idx = rng.gen_range(0, 4);
    return nucl[idx];
    // match idx {
    //     0 => {return "A";},
    //     1 => {return "A";},
    //     2 => {return "A";},
    //     3 => {return "A";},
    // }
}

#[allow(dead_code)]
// replace first occurance of 'N'
// if there are more than one 'N',
// ignore the string
// non-random replacement to avoid
// stochasticity
// https://github.com/COMBINE-lab/salmon/blob/master/src/AlevinUtils.cpp#L789
pub fn cb_string_to_u64(cb_str: &[u8]) -> Result<u64, Box<dyn Error>> {
    let mut cb_id: u64 = 0;
    for (idx, nt) in cb_str.iter().rev().enumerate() {
        let offset = idx * 2;
        match nt {
            65 | 78 => (),              // A | N 00
            67 => cb_id |= 1 << offset, // C 01
            71 => cb_id |= 2 << offset, // G 10
            84 => cb_id |= 3 << offset, // T 11
            _ => panic!("unknown nucleotide {}", nt),
        };
    }

    Ok(cb_id)
}

#[allow(dead_code)]
pub fn tid_2_contig(h: &HeaderView) -> HashMap<u32, String> {
    let mut dict: HashMap<u32, String> = HashMap::with_capacity(46);
    for (i, t) in h
        .target_names()
        .iter()
        .map(|a| str::from_utf8(a).unwrap())
        .enumerate()
    {
        dict.insert(i as u32, t.to_owned());
    }
    dict
}

pub fn bam2rad(input_file: String, rad_file: String, num_threads: u32, log: &slog::Logger) {
    let oname = Path::new(&rad_file);
    let parent = oname.parent().unwrap();
    std::fs::create_dir_all(&parent).unwrap();

    if oname.exists() {
        std::fs::remove_file(oname).expect("could not be deleted");
    }
    let ofile = File::create(&rad_file).unwrap();

    let mut bam = bam::Reader::from_path(&input_file).unwrap();
    let bam_bytes = fs::metadata(&input_file).unwrap().len();
    info! {
        log,
        "Bam file size in bytes {:?}",
        bam_bytes
    };

    if num_threads > 1 {
        bam.set_threads((num_threads as usize) - 1).unwrap();
    } else {
        bam.set_threads(1).unwrap();
    }

    let hdrv = bam.header().to_owned();
    // let tid_lookup: HashMap<u32, String>  = tid_2_contig(&hdrv);
    let mut data = Cursor::new(vec![]);
    // initialize the header (do we need this ?)
    // let mut hdr = libradicl::RADHeader::from_bam_header(&hdrv);
    // number of chunks would be decided when we know
    // let end_header_pos2 = hdr.get_size();
    // info!(
    //     log,
    //     "end header pos {:?}",
    //     end_header_pos2,
    // );

    // file writer
    // let owriter = Arc::new(Mutex::new(BufWriter::with_capacity(1048576, ofile)));
    let mut owriter = BufWriter::with_capacity(1048576, ofile);
    // intermediate buffer

    // write the header
    {
        let is_paired = 0u8;
        data.write_all(&is_paired.to_le_bytes())
            .expect("couldn't write to output file");
        let ref_count = hdrv.target_count() as u64;
        data.write_all(&ref_count.to_le_bytes())
            .expect("couldn't write to output file");
        // create longest buffer
        for t in hdrv.target_names().iter() {
            let name_size = t.len() as u16;
            data.write_all(&name_size.to_le_bytes())
                .expect("coudn't write to output file");
            data.write_all(&t).expect("coudn't write to output file");
        }
        let initial_num_chunks = 0u64;
        data.write_all(&initial_num_chunks.to_le_bytes())
            .expect("coudn't write to output file");
    }

    // test the header
    {
        info!(log, "ref count: {:?} ", hdrv.target_count(),);
    }

    // keep a pointer to header pos
    let end_header_pos =
        data.seek(SeekFrom::Current(0)).unwrap() - std::mem::size_of::<u64>() as u64;

    // check header position
    info!(log, "end header pos: {:?}", end_header_pos,);

    // ### start of tags

    // Tags we will have
    // write the tag meta-information section
    {
        // TODO: get the first record for creating flags

        // file-level
        let mut num_tags = 2u16;
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");
        // type-id
        let mut typeid = 2u8;
        let mut cb_tag_str = "cblen";
        let mut umi_tag_str = "ulen";

        // str - type
        libradicl::write_str_bin(&cb_tag_str, &libradicl::RADIntID::U16, &mut data);
        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");

        // str - type
        libradicl::write_str_bin(&umi_tag_str, &libradicl::RADIntID::U16, &mut data);
        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");

        // read-level
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");
        cb_tag_str = "b";
        umi_tag_str = "u";
        // TODO: make it conditional
        typeid = 3u8;

        libradicl::write_str_bin(&cb_tag_str, &libradicl::RADIntID::U16, &mut data);
        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");
        libradicl::write_str_bin(&umi_tag_str, &libradicl::RADIntID::U16, &mut data);

        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");

        // alignment-level
        num_tags = 1u16;
        data.write_all(&num_tags.to_le_bytes())
            .expect("couldn't write to output file");

        // reference id
        let refid_str = "compressed_ori_refid";
        typeid = 3u8;
        libradicl::write_str_bin(&refid_str, &libradicl::RADIntID::U16, &mut data);
        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");

        let bclen = 16u16;
        let umilen = 10u16;
        data.write_all(&bclen.to_le_bytes())
            .expect("coudn't write to output file");
        data.write_all(&umilen.to_le_bytes())
            .expect("coudn't write to output file");
    }

    // owriter.lock().unwrap().write_all(data.get_ref()).unwrap();
    owriter.write_all(data.get_ref()).unwrap();

    let mut num_output_chunks = 0u64;
    let mut local_nrec = 0u32;
    // let initial_cond : bool = false ;

    // allocate data
    let buf_limit = 10000u32;
    data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 24) as usize));
    data.write_all(&local_nrec.to_le_bytes()).unwrap();
    data.write_all(&local_nrec.to_le_bytes()).unwrap();

    // calculate number of records
    // let mut total_number_of_records = 0u64;
    // for r in bam.records(){
    //     total_number_of_records += 1;
    // }
    // info!(log, "total number of records in bam {:?}", total_number_of_records);

    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
        )
        .progress_chars("╢▌▌░╟");

    let expected_bar_length = bam_bytes / ((buf_limit as u64) * 24);
    // let expected_bar_length = 50u64 ;// bam_bytes / ((buf_limit as u64) * 24);

    let pbar_inner = ProgressBar::new(expected_bar_length as u64);
    pbar_inner.set_style(sty);
    pbar_inner.tick();

    let mut rec = bam::Record::new();
    // history for records that is
    // first seen
    let mut old_qname = String::from("");
    let mut bc = 0u64;
    let mut umi = 0u64;
    let mut tid_list = Vec::<u32>::new();
    //for r in bam.records(){
    while bam.read(&mut rec).unwrap() {
        // let rec = r.unwrap();
        let is_reverse = rec.is_reverse();
        let qname_str = str::from_utf8(rec.qname()).unwrap().to_owned();
        let qname = qname_str;
        let mut tid = rec.tid() as u32;
        if qname == old_qname {
            if !is_reverse {
                tid |= 0x80000000;
            }
            tid_list.push(tid);
            // local_nrec += 1;
            continue;
        }
        // if this is new read and we need to write info
        // for the last read, _unless_ this is the very
        // first read, in which case we shall continue
        if !tid_list.is_empty() {
            assert!(!tid_list.is_empty(), "Trying to write empty tid_list");
            let na = tid_list.len();
            data.write_all(&(na as u32).to_le_bytes()).unwrap();
            //bc
            data.write_all(&(bc as u32).to_le_bytes()).unwrap();
            //umi
            data.write_all(&(umi as u32).to_le_bytes()).unwrap();
            //write tid list
            for t in tid_list.iter() {
                data.write_all(&t.to_le_bytes()).unwrap();
            }
        }

        // dump if we reach the buf_limit
        if local_nrec > buf_limit {
            data.set_position(0);
            let nbytes = (data.get_ref().len()) as u32;
            let nrec = local_nrec;
            // info!(log,"local nrec {:?}-{:?}", local_nrec, tid_list.len());
            data.write_all(&nbytes.to_le_bytes()).unwrap();
            data.write_all(&nrec.to_le_bytes()).unwrap();
            //owriter.lock().unwrap().write_all(data.get_ref()).unwrap();
            owriter.write_all(data.get_ref()).unwrap();
            pbar_inner.inc(1);
            // if num_output_chunks%100 == 0 {
            //    print!("Processed {} chunks\r", num_output_chunks);
            // }

            num_output_chunks += 1;
            local_nrec = 0;
            data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 24) as usize));
            data.write_all(&local_nrec.to_le_bytes()).unwrap();
            data.write_all(&local_nrec.to_le_bytes()).unwrap();

            // for debugging
            // if num_output_chunks > expected_bar_length-1 {
            //     break;
            // }
        }
        // let tname = tid_lookup.get(&(rec.tid() as u32)).unwrap();
        // let qname_string = str::from_utf8(rec.qname()).unwrap();

        // if this is a new read update the old variables
        {
            let bc_string_in = str::from_utf8(rec.aux(b"CR").unwrap().string()).unwrap();
            let umi_string_in = str::from_utf8(rec.aux(b"UR").unwrap().string()).unwrap();

            let bc_string = bc_string_in.replacen('N', "A", 1);
            let umi_string = umi_string_in.replacen('N', "A", 1);
            if let Some(_pos) = bc_string.find('N') {
                continue;
            }
            if let Some(_pos) = umi_string.find('N') {
                continue;
            }

            // convert to u64 following
            // https://github.com/k3yavi/flash/blob/master/src-rs/src/fragments.rs#L162-L176
            bc = cb_string_to_u64(bc_string.as_bytes()).unwrap();
            umi = cb_string_to_u64(umi_string.as_bytes()).unwrap();
            old_qname = qname.clone();
            tid_list.clear();
            if !is_reverse {
                tid |= 0x80000000;
            }
            tid_list.push(tid);
            local_nrec += 1;
        }
        // println!("{:?}\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}",
        //      qname_string,
        //      tname,
        //      bc_string,
        //      umi_string,
        //      num_output_chunks,
        //      local_nrec,
        // );
        //na TODO: make is independed of 16-10 length criterion
        // for debugging
        // if num_output_chunks > expected_bar_length-1 {
        //     break;
        // }
    }

    if local_nrec > 0 {
        // println!("In the residual writing part");
        // first fill the buffer with the last remaining read
        if !tid_list.is_empty() {
            assert!(!tid_list.is_empty(), "Trying to write empty tid_list");
            let na = tid_list.len();
            data.write_all(&(na as u32).to_le_bytes()).unwrap();
            //bc
            data.write_all(&(bc as u32).to_le_bytes()).unwrap();
            //umi
            data.write_all(&(umi as u32).to_le_bytes()).unwrap();
            //write tid list
            for t in tid_list.iter() {
                data.write_all(&t.to_le_bytes()).unwrap();
            }
        }

        data.set_position(0);
        let nbytes = (data.get_ref().len()) as u32;
        let nrec = local_nrec;
        data.write_all(&nbytes.to_le_bytes()).unwrap();
        data.write_all(&nrec.to_le_bytes()).unwrap();
        // owriter.lock().unwrap().write_all(data.get_ref()).unwrap();
        owriter.write_all(data.get_ref()).unwrap();
        num_output_chunks += 1;
    }
    pbar_inner.finish_with_message("wrote all records.");

    // update chunk size
    println!("");
    info!(log, "{:?} chunks written", num_output_chunks,);

    // owriter.lock().unwrap().flush();
    // owriter
    //     .lock()
    //     .unwrap()
    //     .get_ref()
    //     .seek(SeekFrom::Start(
    //         end_header_pos,
    //     ))
    //     .expect("couldn't seek in output file");
    // owriter
    //     .lock()
    //     .unwrap()
    //     .write_all(&num_output_chunks.to_le_bytes())
    //     .expect("couldn't write to output file.");
    owriter.flush().expect("File buffer could not be flushed");
    owriter
        .seek(SeekFrom::Start(end_header_pos))
        .expect("couldn't seek in output file");
    owriter
        .write_all(&num_output_chunks.to_le_bytes())
        .expect("couldn't write to output file.");

    info!(log, "finished writing to {:?}.", rad_file);
}
