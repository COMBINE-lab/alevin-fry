extern crate rust_htslib;
extern crate slog;
extern crate needletail;
extern crate indicatif;
extern crate scroll;

use self::slog::{info};
use self::indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::{BufWriter, Cursor, Seek, SeekFrom, Write};
use std::sync::{Arc, Mutex};
use needletail::bitkmer::*;
use rust_htslib::{bam, bam::Read};
use rust_htslib::bam::HeaderView;
use std::collections::HashMap;
use std::str;
use rand::Rng;

use crate as libradicl;

fn get_random_nucl() -> &'static str {
    let nucl = vec!["A","T","G","C"];
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

pub fn tid_2_contig(h: &HeaderView) -> HashMap<u32, String> {
	let mut dict: HashMap<u32, String> = HashMap::with_capacity(46);
	for (i,t) in h.target_names()
				  .iter().map(|a| str::from_utf8(a).unwrap())
				  .enumerate() {
		dict.insert(i as u32, t.to_owned());
	}
	dict
}

pub fn bam2rad(
    input_file: String,
    rad_file: String,
    log: &slog::Logger,

){
    let ofile = File::create(&rad_file).unwrap();

    let mut bam  = bam::Reader::from_path(&input_file).unwrap();
    let hdrv = bam.header().to_owned();
    let tid_lookup: HashMap<u32, String>  = tid_2_contig(&hdrv);
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
    let owriter = Arc::new(Mutex::new(BufWriter::with_capacity(1048576, ofile)));
    // intermediate buffer

    // write the header
    {
        let is_paired = 0u8 ;
        data
            .write_all(&is_paired.to_le_bytes())
            .expect("couldn't write to output file");
        let ref_count = hdrv.target_count() as u64 ;
        data    
            .write_all(&ref_count.to_le_bytes())
            .expect("couldn't write to output file");
        // create longest buffer
        for (i, t) in hdrv.target_names().iter().enumerate() {
            let name_size = t.len() as u16 ;
            data    
                .write_all(&name_size.to_le_bytes())
                .expect("coudn't write to output file");
            data
                .write_all(&t)
                .expect("coudn't write to output file");
        }
        let initial_num_chunks = 0u64;
        data    
            .write_all(&initial_num_chunks.to_le_bytes())
            .expect("coudn't write to output file");
    }
    // keep a pointer to header pos
    let mut end_header_pos = 
        data.seek(SeekFrom::Current(0)).unwrap()
        - std::mem::size_of::<u64>() as u64;
    
    // check header position
    info!(
        log,
        "end header pos: {:?}",
        end_header_pos,
    );

    // ### start of tags 

    // Tags we will have
    // write the tag meta-information section

    {
        // file-level
        let mut num_tags = 2u16;
        data
            .write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");
        // type-id
        let mut typeid = 2u8;
        let mut cb_tag_str = "cblen";
        let mut umi_tag_str = "ulen";

        // str - type
        libradicl::write_str_bin(
            &cb_tag_str,
            &libradicl::RADIntID::U16,
            &mut data,
        );
        data    
            .write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");
        
        // str - type
        libradicl::write_str_bin(
            &umi_tag_str,
            &libradicl::RADIntID::U16,
            &mut data,
        );
        data    
            .write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");


        // read-level
        data    
            .write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");
        cb_tag_str = "b";
        umi_tag_str = "u";
        // TODO: make it conditional
        typeid = 3u8;
        
        libradicl::write_str_bin(
            &cb_tag_str,
            &libradicl::RADIntID::U16,
            &mut data,
        );
        data    
            .write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");
        libradicl::write_str_bin(
            &umi_tag_str,
            &libradicl::RADIntID::U16,
            &mut data,
        );

        data    
            .write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");
        
        // alignment-level
        num_tags = 1u16 ;
        data.write_all(&num_tags.to_le_bytes())
            .expect("couldn't write to output file");
        
        // reference id
        let refid_str = "compressed_ori_refid";
        typeid = 3u8;
        libradicl::write_str_bin(
            &refid_str,
            &libradicl::RADIntID::U16,
            &mut data,
        );
        data    
            .write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");
        
        let bclen = 16u16;
        let umilen = 10u16;
        data    
            .write_all(&bclen.to_le_bytes())
            .expect("coudn't write to output file");
        data    
            .write_all(&umilen.to_le_bytes())
            .expect("coudn't write to output file");
    }

    owriter.lock().unwrap().write_all(data.get_ref()).unwrap();

    let mut num_output_chunks = 0u64;

    let mut local_nrec = 0u32;
    let initial_cond : bool = false ;
    // reset data
    data = Cursor::new(Vec::<u8>::with_capacity((5000 * 24) as usize));
    data.write_all(&local_nrec.to_le_bytes()).unwrap();
    data.write_all(&local_nrec.to_le_bytes()).unwrap();
    
    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
        )
        .progress_chars("╢▌▌░╟");
    
    let mut expected_bar_length = 1000000u64;
    if let (_, Some(expected_size)) = bam.records().size_hint() {
        expected_bar_length = (expected_size as u64) / 5000;
    }
    let pbar_inner = ProgressBar::new(expected_bar_length as u64);
    pbar_inner.set_style(sty);
    pbar_inner.tick();

    for r in bam.records(){
        let rec = r.unwrap();
        let is_reverse = rec.is_reverse();
        let tid = rec.tid() as u32;
        // let tname = tid_lookup.get(&(rec.tid() as u32)).unwrap();
        // let qname_string = str::from_utf8(rec.qname()).unwrap();
        let bc_string_in = str::from_utf8( rec.aux(b"CR").unwrap().string() ).unwrap();
        let umi_string_in = str::from_utf8( rec.aux(b"UR").unwrap().string() ).unwrap();
        let bclen = bc_string_in.len();
        let umilen = umi_string_in.len();

        //println!("{:?}",qname_string);

        let bc_string = bc_string_in.replacen('N', get_random_nucl(), bclen);
        let umi_string = umi_string_in.replacen('N', get_random_nucl(), umilen);
        let mut bc_bytes = BitNuclKmer::new(bc_string.as_bytes(), bclen as u8, false);
        let (_, bc, _) = bc_bytes.next().expect("can't extract barcode");
        
        let mut umi_bytes = BitNuclKmer::new(umi_string.as_bytes(), umilen as u8, false);
        let (_, umi, _) = umi_bytes.next().expect("can't extract umi");

        // println!("{:?}\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}",
        //      qname_string,
        //      tname,
        //      bc_string,
        //      umi_string,
        //      num_output_chunks,
        //      local_nrec,
        // );
        //na
        data.write_all(&(1 as u32).to_le_bytes()).unwrap();
        //bc
        data.write_all(&(bc.0 as u32).to_le_bytes()).unwrap();
        //umi        
        data.write_all(&(umi.0 as u32).to_le_bytes()).unwrap();
        let mut tid_dir = tid | 0x80000000; 
        if is_reverse{
            tid_dir = tid | 0x00000000 ;
        }
        data.write_all(&tid_dir.to_le_bytes()).unwrap();
        local_nrec += 1;


        if local_nrec > 5000{
            data.set_position(0);
            let nbytes = (data.get_ref().len()) as u32;
            let nrec = local_nrec;
            data.write_all(&nbytes.to_le_bytes()).unwrap();
            data.write_all(&nrec.to_le_bytes()).unwrap();
            owriter.lock().unwrap().write_all(data.get_ref()).unwrap();
            pbar_inner.inc(1);

            num_output_chunks += 1;
            local_nrec = 0;
            data = Cursor::new(Vec::<u8>::with_capacity((5000 * 24) as usize));
            data.write_all(&local_nrec.to_le_bytes()).unwrap();
            data.write_all(&local_nrec.to_le_bytes()).unwrap();
        }
    }

    if local_nrec > 0 {
        data.set_position(0);
        let nbytes = (data.get_ref().len()) as u32;
        let nrec = local_nrec;
        data.write_all(&nbytes.to_le_bytes()).unwrap();
        data.write_all(&nrec.to_le_bytes()).unwrap();
        owriter.lock().unwrap().write_all(data.get_ref()).unwrap();
        num_output_chunks += 1;
    }
    pbar_inner.finish_with_message("wrote all records.");

    // update chunk size
    owriter
        .lock()
        .unwrap()
        .get_ref()
        .seek(SeekFrom::Start(
            end_header_pos,
        ))
        .expect("couldn't seek in output file");
    owriter
        .lock()
        .unwrap()
        .write_all(&num_output_chunks.to_le_bytes())
        .expect("couldn't write to output file.");

    info!(log, "finished collating input rad file {:?}.", rad_file);
    
}