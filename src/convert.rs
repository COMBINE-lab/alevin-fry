/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use anyhow::bail;
use indicatif::{ProgressBar, ProgressStyle};
use slog::{crit, info};

//use num_format::{Locale};
use std::fs;
use std::fs::File;
use std::io::{stdout, BufReader, BufWriter, Cursor, Seek, SeekFrom, Write};
// use std::sync::{Arc, Mutex};
//

use noodles::bam as nbam;
use noodles::{bgzf, sam};
use sam::alignment::record::data::field::{tag::Tag as SamTag, value::Value as SamTagValue};

use libradicl::rad_types::{self, RadIntId, RadType, TagDesc, TagSection, TagSectionLabel};
use libradicl::utils::MASK_LOWER_31_U32;
use libradicl::{
    chunk,
    header::{RadHeader, RadPrelude},
    record::{AlevinFryReadRecord, AlevinFryRecordContext},
};

use needletail::bitkmer::*;
use rand::Rng;
use std::error::Error;
use std::path::Path;
use std::str;

// pub fn reset_signal_pipe_handler() -> Result<()> {
//     #[cfg(target_family = "unix")]
//     {
//         use nix::sys::signal;

//         unsafe {
//             signal::signal(signal::Signal::SIGPIPE, signal::SigHandler::SigDfl)
//                 .map_err(|e| Error::Other(e.to_string()))?;
//         }
//     }

//     Ok(())
// }

#[allow(dead_code)]
fn get_random_nucl() -> &'static str {
    let nucl = ["A", "T", "G", "C"];
    let mut rng = rand::thread_rng();
    let idx = rng.gen_range(0..4);
    nucl[idx]
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
// https://github.com/k3yavi/flash/blob/master/src-rs/src/fragments.rs#L162-L176
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

#[inline]
fn write_barcode<W: Write>(barcode_t: &RadType, barcode: u64, w: &mut W) -> anyhow::Result<()> {
    match barcode_t {
        RadType::Int(int_t) => match int_t {
            RadIntId::U8 => {
                let v = barcode as u8;
                w.write_all(&v.to_le_bytes())?;
            }
            RadIntId::U16 => {
                let v = barcode as u16;
                w.write_all(&v.to_le_bytes())?;
            }
            RadIntId::U32 => {
                let v = barcode as u32;
                w.write_all(&v.to_le_bytes())?;
            }
            RadIntId::U64 => {
                let v = barcode;
                w.write_all(&v.to_le_bytes())?;
            }
            RadIntId::U128 => {
                todo!("support for u128-encoded barcodes is not yet available");
            }
        },
        _ => bail!("invalid type!"),
    }
    Ok(())
}

#[inline]
fn write_list<W: Write>(
    tid_list: &[u32],
    bct: &RadType,
    bc: u64,
    umit: &RadType,
    umi: u64,
    w: &mut W,
) -> anyhow::Result<()> {
    assert!(!tid_list.is_empty(), "Trying to write empty tid_list");
    let na = tid_list.len() as u32;
    w.write_all(&na.to_le_bytes()).unwrap();
    //bc
    write_barcode(bct, bc, w)?;
    //umi
    write_barcode(umit, umi, w)?;
    //write tid list
    for t in tid_list.iter() {
        w.write_all(&t.to_le_bytes()).unwrap();
    }
    Ok(())
}

pub fn bam2rad<P1, P2>(
    input_file: P1,
    rad_file: P2,
    num_threads: u32,
    log: &slog::Logger,
) -> anyhow::Result<()>
where
    P1: AsRef<Path>,
    P2: AsRef<Path>,
{
    const CR: SamTag = SamTag::new(b'C', b'R');
    const UR: SamTag = SamTag::new(b'U', b'R');

    let oname = Path::new(rad_file.as_ref());
    let parent = oname.parent().unwrap();
    std::fs::create_dir_all(parent).unwrap();

    if oname.exists() {
        std::fs::remove_file(oname).expect("could not be deleted");
    }
    let ofile = File::create(rad_file.as_ref()).unwrap();

    // number of bytes in the input BAM file
    let bam_bytes = fs::metadata(&input_file).unwrap().len();
    info! {
    log,
    "Bam file size in bytes {:?}",
    bam_bytes
    };

    // reading input BAM using Noodles
    // example from https://github.com/zaeleus/noodles/issues/227
    let file = File::open(&input_file)?;

    let mut reader: Box<dyn sam::alignment::io::Read<_>> =
        match input_file.as_ref().extension().and_then(|ext| ext.to_str()) {
            Some("bam") | Some("BAM") => {
                let decomp_threads = std::num::NonZeroUsize::new(if num_threads > 1 {
                    (num_threads - 1) as usize
                } else {
                    1_usize
                })
                .expect("invalid nonzero usize");

                println!(
                    "parsing BAM file using {:?} decompression threads",
                    decomp_threads
                );

                let decoder: Box<dyn std::io::BufRead> = Box::new(
                    bgzf::MultithreadedReader::with_worker_count(decomp_threads, file),
                );
                Box::new(nbam::io::Reader::from(decoder))
            }
            Some("sam") | Some("SAM") => {
                let inner: Box<dyn std::io::BufRead> = Box::new(BufReader::new(file));
                Box::new(sam::io::Reader::from(inner))
            }
            _ => {
                bail!("unsupported input file format, must end with bam/BAM or sam/SAM");
            }
        };

    let hdrv = reader.read_alignment_header()?;
    let mut data = Cursor::new(vec![]);

    // intermediate buffer
    let mut owriter = BufWriter::with_capacity(1_048_576, ofile);

    // write the header
    {
        // NOTE: The is_paired flag is not present in the
        // SAM file and so isn't meaningful as written here.
        // Also, the num_chunks will be 0 in this header
        // currently.
        let rad_header = RadHeader::from_bam_header(&hdrv);
        rad_header.write(&mut data)?;
    }

    // test the header
    {
        info!(log, "ref count: {:?} ", hdrv.reference_sequences().len(),);
    }

    // keep a pointer to header pos we need this to fill in the correct
    // number of chunks later on.
    let end_header_pos = data.stream_position().unwrap() - std::mem::size_of::<u64>() as u64;

    // check header position
    info!(log, "end header pos: {:?}", end_header_pos,);

    let mut record_it = reader.alignment_records(&hdrv).peekable();

    // ### start of tags
    // get the first record for creating flags
    let rec_res = record_it.peek();
    let first_record_exists = rec_res.is_some();
    if !first_record_exists {
        crit!(log, "bam file had no records!");
        std::process::exit(1);
    }
    let rec = match rec_res.unwrap() {
        Ok(x) => x.as_ref(),
        Err(e) => bail!("{}", e),
    };

    let bc_typeid: RadType;
    let umi_typeid: RadType;

    use bstr::BStr;
    // Tags we will have
    // write the tag meta-information section
    {
        // file-level
        let mut file_tags = TagSection::new_with_label(TagSectionLabel::FileTags);
        file_tags.add_tag_desc(TagDesc {
            name: "cblen".to_owned(),
            typeid: RadType::Int(RadIntId::U16),
        });
        file_tags.add_tag_desc(TagDesc {
            name: "ulen".to_owned(),
            typeid: RadType::Int(RadIntId::U16),
        });

        file_tags.write(&mut data)?;

        // read-level
        let flag_data = rec.data();
        let bc_string_in: &str = if let Some(Ok(bcs)) = flag_data.get(&CR) {
            match bcs {
                SamTagValue::String(bstr) => str::from_utf8(<BStr as AsRef<[u8]>>::as_ref(bstr))?,
                _ => {
                    bail!("cannot convert non-string (Z) tag into barcode string.");
                }
            }
        } else {
            panic!("Input record missing CR tag!")
        };

        let umi_string_in: &str = if let Some(Ok(umis)) = flag_data.get(&UR) {
            match umis {
                SamTagValue::String(bstr) => str::from_utf8(<BStr as AsRef<[u8]>>::as_ref(bstr))?,
                _ => {
                    bail!("cannot convert non-string (Z) tag into umi string.");
                }
            }
        } else {
            panic!("Input record missing UR tag!")
        };
        let bclen = bc_string_in.len() as u16;
        let umilen = umi_string_in.len() as u16;

        // type is conditional on barcode and umi length
        bc_typeid = match bclen {
            1..=4 => RadType::Int(rad_types::RadIntId::U8),
            5..=8 => RadType::Int(rad_types::RadIntId::U16),
            9..=16 => RadType::Int(rad_types::RadIntId::U32),
            17..=32 => RadType::Int(rad_types::RadIntId::U64),
            l => {
                crit!(log, "cannot encode barcode of length {} > 32", l);
                std::process::exit(1);
            }
        };

        umi_typeid = match umilen {
            1..=4 => RadType::Int(rad_types::RadIntId::U8),
            5..=8 => RadType::Int(rad_types::RadIntId::U16),
            9..=16 => RadType::Int(rad_types::RadIntId::U32),
            17..=32 => RadType::Int(rad_types::RadIntId::U64),
            l => {
                crit!(log, "cannot encode umi of length {} > 32", l);
                std::process::exit(1);
            }
        };

        let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
        read_tags.add_tag_desc(TagDesc {
            name: "b".to_owned(),
            typeid: bc_typeid,
        });
        read_tags.add_tag_desc(TagDesc {
            name: "u".to_owned(),
            typeid: umi_typeid,
        });
        read_tags.write(&mut data)?;

        // alignment-level
        let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
        aln_tags.add_tag_desc(TagDesc {
            name: "compressed_ori_refid".to_owned(),
            typeid: RadType::Int(RadIntId::U32),
        });
        aln_tags.write(&mut data)?;

        // done with tag descriptions
        // now write the values associated with the file-level tags
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
    let buf_limit = 10_000u32;
    data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 24) as usize));
    data.write_all(&local_nrec.to_le_bytes()).unwrap();
    data.write_all(&local_nrec.to_le_bytes()).unwrap();

    // empiricaly derived factor of size of bam vs rad
    // encoding for records.
    let approx_bam_to_rad_factor = 6.258_f64;

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
        .expect("ProgressStyle template was invalid.")
        .progress_chars("╢▌▌░╟");

    let expected_bar_length =
        bam_bytes / (((buf_limit as f64) * 24_f64) * approx_bam_to_rad_factor).round() as u64;

    let pbar_inner = ProgressBar::new(expected_bar_length);
    pbar_inner.set_style(sty);
    pbar_inner.tick();

    // history for records that is
    // first seen
    let mut old_qname = String::from("");
    let mut bc = 0u64;
    let mut umi = 0u64;
    let mut tid_list = Vec::<u32>::new();
    //for r in bam.records(){
    loop {
        let rec_res = record_it.next();
        if rec_res.is_none() {
            break;
        }
        // the iterator returns a result, so
        // make sure it's an Ok variant here.
        let rec = match rec_res.unwrap() {
            Ok(x) => x,
            Err(e) => bail!("{}", e),
        };

        let flags = rec.flags()?;

        let is_reverse = flags.is_reverse_complemented();
        let qname_str = rec.name().expect("valid name").to_string().to_owned();
        let qname = qname_str;
        let mut tid = rec.reference_sequence_id(&hdrv).unwrap().unwrap() as u32;
        if qname == old_qname {
            if !is_reverse {
                tid |= MASK_LOWER_31_U32;
            }
            tid_list.push(tid);
            // local_nrec += 1;
            continue;
        }

        // if this is new read and we need to write info
        // for the last read, _unless_ this is the very
        // first read, in which case we shall continue
        if !tid_list.is_empty() {
            write_list(&tid_list, &bc_typeid, bc, &umi_typeid, umi, &mut data)?;
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
        }

        // if this is a new read update the old variables
        {
            let flag_data = rec.data();
            let bc_string_in: &str = if let Some(Ok(bcs)) = flag_data.get(&CR) {
                match bcs {
                    SamTagValue::String(bstr) => {
                        str::from_utf8(<BStr as AsRef<[u8]>>::as_ref(bstr))?
                    }
                    _ => {
                        bail!("cannot convert non-string (Z) tag into umi string.");
                    }
                }
            } else {
                panic!("Input record missing CR tag!")
            };

            let umi_string_in: &str = if let Some(Ok(umis)) = flag_data.get(&UR) {
                match umis {
                    SamTagValue::String(bstr) => {
                        str::from_utf8(<BStr as AsRef<[u8]>>::as_ref(bstr))?
                    }
                    _ => {
                        bail!("cannot convert non-string (Z) tag into umi string.");
                    }
                }
            } else {
                panic!("Input record missing UR tag!")
            };

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
            old_qname.clone_from(&qname);
            tid_list.clear();
            if !is_reverse {
                tid |= MASK_LOWER_31_U32;
            }
            tid_list.push(tid);
            local_nrec += 1;
        }
    }

    if local_nrec > 0 {
        // println!("In the residual writing part");
        // first fill the buffer with the last remaining read
        if !tid_list.is_empty() {
            write_list(&tid_list, &bc_typeid, bc, &umi_typeid, umi, &mut data)?;
        }

        data.set_position(0);
        let nbytes = (data.get_ref().len()) as u32;
        let nrec = local_nrec;
        data.write_all(&nbytes.to_le_bytes()).unwrap();
        data.write_all(&nrec.to_le_bytes()).unwrap();

        owriter.write_all(data.get_ref()).unwrap();
        num_output_chunks += 1;
    }
    pbar_inner.finish_with_message("wrote all records.");

    // update chunk size
    println!();
    info!(log, "{:?} chunks written", num_output_chunks,);

    owriter.flush().expect("File buffer could not be flushed");
    owriter
        .seek(SeekFrom::Start(end_header_pos))
        .expect("couldn't seek in output file");
    owriter
        .write_all(&num_output_chunks.to_le_bytes())
        .expect("couldn't write to output file.");

    info!(log, "finished writing to {:?}.", rad_file.as_ref());
    Ok(())
}

pub fn view<P>(rad_file: P, print_header: bool, log: &slog::Logger)
where
    P: AsRef<Path>,
{
    let _read_num = view2(rad_file, print_header, log).unwrap();
}
pub fn view2<P>(rad_file: P, print_header: bool, log: &slog::Logger) -> anyhow::Result<u64>
where
    P: AsRef<Path>,
{
    let i_file = File::open(rad_file).unwrap();
    let mut br = BufReader::new(i_file);
    let prelude = RadPrelude::from_bytes(&mut br)?;
    let hdr = &prelude.hdr;
    let rl_tags = &prelude.read_tags;

    // right now, we only handle BC and UMI types of U8—U64, so validate that
    const BNAME: &str = "b";
    const UNAME: &str = "u";

    let mut bct: Option<RadType> = None;
    let mut umit: Option<RadType> = None;

    for rt in &rl_tags.tags {
        // if this is one of our tags
        if rt.name == BNAME || rt.name == UNAME {
            if !rt.typeid.is_int_type() {
                crit!(
                    log,
                    "currently only RAD types 1--4 are supported for 'b' and 'u' tags."
                );
                std::process::exit(libradicl::exit_codes::EXIT_UNSUPPORTED_TAG_TYPE);
            }

            if rt.name == BNAME {
                bct = Some(rt.typeid);
            }
            if rt.name == UNAME {
                umit = Some(rt.typeid);
            }
        }
    }
    assert!(bct.is_some(), "barcode type tag was missing!");
    assert!(umit.is_some(), "umi type tag was missing!");

    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    info!(log, "File-level tag map {:?}", file_tag_map);

    let barcode_tag = file_tag_map
        .get("cblen")
        .expect("tag map must contain cblen");
    let barcode_len: u16 = barcode_tag.try_into()?;

    let umi_tag = file_tag_map.get("ulen").expect("tag map must contain ulen");
    let umi_len: u16 = umi_tag.try_into()?;

    let mut num_reads: u64 = 0;
    let record_context = prelude.get_record_context::<AlevinFryRecordContext>()?;

    let stdout = stdout(); // get the global stdout entity
    let stdout_l = stdout.lock();
    let mut handle = BufWriter::new(stdout_l); // optional: wrap that handle in a buffer

    if print_header {
        for i in 0usize..hdr.ref_names.len() {
            match writeln!(handle, "{}:{}", i, hdr.ref_names[i]) {
                Ok(_) => {}
                Err(_) => {
                    return Ok(i as u64);
                }
            };
        }
    }

    let mut id = 0usize;
    for _ in 0..(hdr.num_chunks as usize) {
        let c = chunk::Chunk::<AlevinFryReadRecord>::from_bytes(&mut br, &record_context);
        for read in c.reads.iter() {
            let bc_mer: BitKmer = (read.bc, barcode_len as u8);
            let umi_mer: BitKmer = (read.umi, umi_len as u8);

            // let umi = str::from_utf8(&umi_).unwrap();
            let num_entries = read.refs.len();
            for i in 0usize..num_entries {
                let tid = &hdr.ref_names[read.refs[i] as usize];
                match writeln!(
                    handle,
                    "ID:{}\tHI:{}\tNH:{}\tCB:{}\tUMI:{}\tDIR:{:?}\t{}",
                    id,
                    i + 1,
                    num_entries,
                    unsafe { std::str::from_utf8_unchecked(&bitmer_to_bytes(bc_mer)[..]) },
                    unsafe { std::str::from_utf8_unchecked(&bitmer_to_bytes(umi_mer)[..]) },
                    read.dirs[i],
                    tid,
                ) {
                    Ok(_) => {
                        num_reads += 1;
                    }
                    Err(_) => {
                        // head broken pipe
                        // https://github.com/rust-lang/rust/issues/46016#issuecomment-605624865
                        return Ok(num_reads);
                    }
                };
            }
            id += 1;
        }
    }

    Ok(num_reads)
}
