extern crate needletail;
extern crate num;
extern crate scroll;
extern crate slog;

use self::slog::info;
use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};

use crate as libradicl;


pub fn read_bulkRAD(
    rad_file: String,
    num_threads: u32,
    log: &slog::Logger,
) -> Result<(), Box<dyn std::error::Error>> {

    let i_file = File::open(&rad_file).unwrap();
    let mut br = BufReader::new(i_file);

    let hdr = libradicl::RADBulkHeader::from_bytes(&mut br);

    // let end_header_pos =
        // br.get_ref().seek(SeekFrom::Current(0)).unwrap() - (br.buffer().len() as u64);

    info!(
        log,
        "paired : {:?}, ref_count : {:?}, num_chunks : {:?}",
        hdr.is_paired != 0,
        hdr.ref_count,
        hdr.num_chunks
    );
    // file-level
    let fl_tags = libradicl::TagSection::from_bulk_bytes(&mut br);
    info!(log, "read {:?} file-level tags", fl_tags.tags.len());
    // read-level
    let rl_tags = libradicl::TagSection::from_bulk_bytes(&mut br);
    info!(log, "read {:?} read-level tags", rl_tags.tags.len());
    // alignment-level
    let al_tags = libradicl::TagSection::from_bulk_bytes(&mut br);
    info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

    // let ft_vals = libradicl::FileTags::from_bytes(&mut br);
    // info!(log, "File-level tag values {:?}", ft_vals);

    let mut num_reads = 0;
    // for _ in 0..(hdr.num_chunks as usize) {
        
    // while !br.eo
    let mut cntr = 0;
    loop {
        match libradicl::BulkChunk::from_bytes(&mut br, & al_tags, & rl_tags, &mut cntr) {
            Ok(v) => num_reads += v.reads.len(),
            Err(_) => break,
        }
    }

    // let pos = br.get_ref().seek(SeekFrom::Current(0)).unwrap() - (br.buffer().len() as u64);

    info!(log, "finished collating input rad file {:?}. Observed {:?} reads", rad_file, num_reads);
    Ok(())
}