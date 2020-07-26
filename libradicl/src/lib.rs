// scroll now, explore nom later
extern crate fasthash;
extern crate scroll;

use scroll::Pread;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::io::{BufWriter, Write};
use std::vec::Vec;

pub mod collate;
pub mod quant;
pub mod schema;
pub mod utils;
pub mod pugutils;
pub mod em;

// Name of the program, to be used in diagnostic messages.
static LIB_NAME: &str = "libradicl";

pub fn lib_name() -> &'static str {
    LIB_NAME
}

pub struct RADHeader {
    pub is_paired: u8,
    pub ref_count: u64,
    pub ref_names: Vec<String>,
    pub num_chunks: u64,
}

pub struct TagDesc {
    pub name: String,
    pub typeid: u8,
}

pub struct TagSection {
    pub tags: Vec<TagDesc>,
}

// The below are currently hard-coded
// until we decide how to solve this
// generally
#[derive(Debug)]
pub struct FileTags {
    pub bclen: u16,
    pub umilen: u16,
}
#[derive(Debug)]
pub struct ReadRecord {
    pub bc: u64,
    pub umi: u64,
    pub dirs: Vec<bool>,
    pub refs: Vec<u32>,
}
#[derive(Debug)]
pub struct Chunk {
    pub nbytes: u32,
    pub nrec: u32,
    pub reads: Vec<ReadRecord>,
}

#[derive(Debug)]
pub struct CorrectedCBChunk {
    remaining_records: u64,
    corrected_bc: u64,
    umis: Vec<u64>,
    ref_offsets: Vec<u32>,
    ref_ids: Vec<u32>,
}

impl CorrectedCBChunk {
    pub fn from_counter(num_remain: u64) -> CorrectedCBChunk {
        CorrectedCBChunk {
            remaining_records: num_remain,
            corrected_bc: 0,
            umis: Vec::<u64>::with_capacity(num_remain as usize),
            ref_offsets: Vec::<u32>::with_capacity(num_remain as usize),
            ref_ids: Vec::<u32>::with_capacity(3 * num_remain as usize),
        }
    }
}

pub enum RADIntID {
    U8,
    U16,
    U32,
    U64,
}

pub struct ChunkConfig {
    pub num_chunks: u64,
    pub bc_type: u8,
    pub umi_type: u8,
}

pub fn collect_records<T: Read>(
    reader: &mut BufReader<T>,
    chunk_config: &ChunkConfig,
    correct_map: &HashMap<u64, u64>,
    output_cache: &mut HashMap<u64, CorrectedCBChunk>,
) {
    // NOTE: since the chunks are independent, this part could be multithreaded
    for _ in 0..(chunk_config.num_chunks as usize) {
        match (chunk_config.bc_type, chunk_config.umi_type) {
            (3, 3) => {
                process_corrected_cb_chunk(
                    reader,
                    RADIntID::U32,
                    RADIntID::U32,
                    correct_map,
                    output_cache,
                );
            }
            (3, 4) => {
                process_corrected_cb_chunk(
                    reader,
                    RADIntID::U32,
                    RADIntID::U64,
                    correct_map,
                    output_cache,
                );
            }
            (4, 3) => {
                process_corrected_cb_chunk(
                    reader,
                    RADIntID::U64,
                    RADIntID::U32,
                    correct_map,
                    output_cache,
                );
            }
            (4, 4) => {
                process_corrected_cb_chunk(
                    reader,
                    RADIntID::U64,
                    RADIntID::U32,
                    correct_map,
                    output_cache,
                );
            }
            (_, _) => println!("types not supported"),
        }
    }
}

fn read_into_u64<T: Read>(reader: &mut BufReader<T>, rt: &RADIntID) -> u64 {
    let mut rbuf = [0u8; 8];
    let v: u64;
    match rt {
        RADIntID::U8 => {
            reader.read_exact(&mut rbuf[0..1]).unwrap();
            v = rbuf.pread::<u8>(0).unwrap() as u64;
        }
        RADIntID::U16 => {
            reader.read_exact(&mut rbuf[0..2]).unwrap();
            v = rbuf.pread::<u16>(0).unwrap() as u64;
        }
        RADIntID::U32 => {
            reader.read_exact(&mut rbuf[0..4]).unwrap();
            v = rbuf.pread::<u32>(0).unwrap() as u64;
        }
        RADIntID::U64 => {
            reader.read_exact(&mut rbuf[0..8]).unwrap();
            v = rbuf.pread::<u64>(0).unwrap();
        }
    }
    v
}

impl ReadRecord {
    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>, bct: &RADIntID, umit: &RADIntID) -> Self {
        let mut rbuf = [0u8; 255];

        reader.read_exact(&mut rbuf[0..4]).unwrap();
        let na = rbuf.pread::<u32>(0).unwrap();
        let bc = read_into_u64(reader, bct);
        let umi = read_into_u64(reader, umit);

        let mut rec = Self {
            bc,
            umi,
            dirs: Vec::with_capacity(na as usize),
            refs: Vec::with_capacity(na as usize),
        };

        for _ in 0..(na as usize) {
            reader.read_exact(&mut rbuf[0..4]).unwrap();
            let v = rbuf.pread::<u32>(0).unwrap();
            let dir = (v & 0x80000000) != 0;
            rec.dirs.push(dir);
            rec.refs.push(v & 0x7FFFFFFF);
        }

        rec
    }

    pub fn from_bytes_keep_ori<T: Read>(
        reader: &mut BufReader<T>,
        bct: &RADIntID,
        umit: &RADIntID,
    ) -> Self {
        let mut rbuf = [0u8; 255];

        reader.read_exact(&mut rbuf[0..4]).unwrap();
        let na = rbuf.pread::<u32>(0).unwrap();

        let bc = read_into_u64(reader, bct);
        let umi = read_into_u64(reader, umit);

        let mut rec = Self {
            bc,
            umi,
            dirs: Vec::with_capacity(na as usize),
            refs: Vec::with_capacity(na as usize),
        };

        for _ in 0..(na as usize) {
            reader.read_exact(&mut rbuf[0..4]).unwrap();
            let v = rbuf.pread::<u32>(0).unwrap();
            //let dir = (v & 0x80000000) != 0;
            //rec.dirs.push( dir );
            rec.refs.push(v);
        }

        rec
    }
}

pub fn process_corrected_cb_chunk<T: Read>(
    reader: &mut BufReader<T>,
    bct: RADIntID,
    umit: RADIntID,
    correct_map: &HashMap<u64, u64>,
    output_cache: &mut HashMap<u64, CorrectedCBChunk>,
) {
    /*
    let (s, r) = unbounded();

    thread::spawn(
      |correct_map, output_cache| {
        for c in r.iter() {
             // update the corresponding corrected chunk entry
             v.umis.push(rr.umi);
             let ref_offset = v.ref_ids.len() as u32;
             v.ref_offsets.push(ref_offset);
             v.ref_ids.extend(&rr.refs);
             v.remaining_records -= 1;
        }
      }
    )
    */
    let mut buf = [0u8; 8];

    // get the number of bytes and records for
    // the next chunk
    reader.read_exact(&mut buf).unwrap();
    let _nbytes = buf.pread::<u32>(0).unwrap();
    let nrec = buf.pread::<u32>(4).unwrap();

    // for each record, read it
    for _ in 0..(nrec as usize) {
        let rr = ReadRecord::from_bytes_keep_ori(reader, &bct, &umit);

        // if this record had a correct or correctable barcode
        if let Some(corrected_id) = correct_map.get(&rr.bc) {
            if let Some(v) = output_cache.get_mut(corrected_id) {
                // update the corresponding corrected chunk entry
                v.umis.push(rr.umi);
                let ref_offset = v.ref_ids.len() as u32;
                v.ref_offsets.push(ref_offset);
                v.ref_ids.extend(&rr.refs);
                v.remaining_records -= 1;
            }
        }
    }
}

fn as_u8_slice(v: &[u32]) -> &[u8] {
    unsafe {
        std::slice::from_raw_parts(
            v.as_ptr() as *const u8,
            v.len() * std::mem::size_of::<u32>(),
        )
    }
}

pub fn dump_output_cache(
    owriter: &mut BufWriter<File>,
    output_cache: &HashMap<u64, CorrectedCBChunk>,
) {
    for (_bc, chunk) in output_cache.iter() {
        // number of bytes
        let mut nbytes: u32 = 0;
        nbytes += (chunk.ref_ids.len() * 4) as u32;
        // umis
        nbytes += (chunk.umis.len() * 4) as u32;
        // barcodes
        nbytes += (chunk.umis.len() * 4) as u32;
        // num alignment fields
        nbytes += (chunk.umis.len() * 4) as u32;

        let nrec = chunk.umis.len() as u32;

        owriter
            .write_all(&nbytes.to_le_bytes())
            .expect("couldn't write output.");
        owriter
            .write_all(&nrec.to_le_bytes())
            .expect("couldn't write output.");

        for i in 0..chunk.umis.len() {
            let s = chunk.ref_offsets[i];
            let e = if i == chunk.umis.len() - 1 {
                chunk.ref_ids.len() as u32
            } else {
                chunk.ref_offsets[i + 1]
            };

            // num alignments
            let num_aln = (e - s) as u32;
            // cb
            let cb = chunk.corrected_bc as u32;
            // umi
            let umi = chunk.umis[i] as u32;

            owriter
                .write_all(&num_aln.to_le_bytes())
                .expect("couldn't write output.");
            owriter
                .write_all(&cb.to_le_bytes())
                .expect("couldn't write output.");
            owriter
                .write_all(&umi.to_le_bytes())
                .expect("couldn't write output.");
            owriter
                .write_all(as_u8_slice(&chunk.ref_ids[(s as usize)..(e as usize)]))
                .expect("couldn't write output.");
        }
    }
}

impl Chunk {
    pub fn read_header<T: Read>(reader: &mut BufReader<T>) -> (u32, u32) {
        let mut buf = [0u8; 8];

        reader.read_exact(&mut buf).unwrap();
        let nbytes = buf.pread::<u32>(0).unwrap();
        let nrec = buf.pread::<u32>(4).unwrap();
        (nbytes, nrec)
    }

    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>, bct: RADIntID, umit: RADIntID) -> Self {
        let mut buf = [0u8; 8];

        reader.read_exact(&mut buf).unwrap();
        let nbytes = buf.pread::<u32>(0).unwrap();
        let nrec = buf.pread::<u32>(4).unwrap();
        let mut c = Self {
            nbytes,
            nrec,
            reads: Vec::with_capacity(nrec as usize),
        };

        for _ in 0..(nrec as usize) {
            c.reads.push(ReadRecord::from_bytes(reader, &bct, &umit));
        }

        c
    }
}

impl FileTags {
    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>) -> Self {
        let mut buf = [0u8; 4];
        reader.read_exact(&mut buf).unwrap();

        Self {
            bclen: buf.pread::<u16>(0).unwrap(),
            umilen: buf.pread::<u16>(2).unwrap(),
        }
    }
}

impl TagDesc {
    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>) -> TagDesc {
        // space for the string length (1 byte)
        // the longest string possible (255 char)
        // and the typeid
        let mut buf = [0u8; 257];
        reader.read_exact(&mut buf[0..2]).unwrap();
        let str_len = buf.pread::<u16>(0).unwrap() as usize;

        // read str_len + 1 to get the type id that follows the string
        reader.read_exact(&mut buf[0..str_len + 1]).unwrap();
        TagDesc {
            name: std::str::from_utf8(&buf[0..str_len]).unwrap().to_string(),
            typeid: buf.pread(str_len).unwrap(),
        }
    }
}

impl TagSection {
    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>) -> TagSection {
        let mut buf = [0u8; 2];
        reader.read_exact(&mut buf).unwrap();
        let num_tags = buf.pread::<u16>(0).unwrap() as usize;

        let mut ts = TagSection {
            tags: Vec::with_capacity(num_tags),
        };

        for _ in 0..num_tags {
            ts.tags.push(TagDesc::from_bytes(reader));
        }

        ts
    }
}

impl RADHeader {
    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>) -> RADHeader {
        let mut rh = RADHeader {
            is_paired: 0,
            ref_count: 0,
            ref_names: vec![],
            num_chunks: 0,
        };

        // size of the longest allowable string.
        let mut buf = [0u8; 65536];
        reader.read_exact(&mut buf[0..9]).unwrap();
        rh.is_paired = buf.pread(0).unwrap();
        rh.ref_count = buf.pread::<u64>(1).unwrap();

        // we know how many names we will read in.
        rh.ref_names.reserve_exact(rh.ref_count as usize);

        let mut num_read = 0u64;
        while num_read < rh.ref_count {
            reader.read_exact(&mut buf[0..2]).unwrap();
            let l: usize = buf.pread::<u16>(0).unwrap() as usize;
            reader.read_exact(&mut buf[0..l]).unwrap();
            rh.ref_names
                .push(std::str::from_utf8(&buf[0..l]).unwrap().to_string());
            num_read += 1;
        }

        reader.read_exact(&mut buf[0..8]).unwrap();
        rh.num_chunks = buf.pread::<u64>(0).unwrap();
        rh
    }
}

pub fn update_barcode_hist(
    hist: &mut HashMap<u64, u64, fasthash::RandomState<fasthash::sea::Hash64>>,
    chunk: &Chunk,
) {
    for r in &chunk.reads {
        *hist.entry(r.bc).or_insert(0) += 1;
    }
}

pub fn permit_list_from_threshold(
    hist: &HashMap<u64, u64, fasthash::RandomState<fasthash::sea::Hash64>>,
    min_freq: u64,
) -> Vec<u64> {
    let valid_bc: Vec<u64> = hist
        .iter()
        .filter_map(|(k, v)| if v >= &min_freq { Some(*k) } else { None })
        .collect();
    valid_bc
}
