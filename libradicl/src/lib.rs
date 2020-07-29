// scroll now, explore nom later
extern crate fasthash;
extern crate needletail;
extern crate num;
extern crate quickersort;
extern crate scroll;

use bio_types::strand::*;
use needletail::bitkmer::*;
use num::cast::AsPrimitive;
use scroll::Pread;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::io::{BufWriter, Write};
use std::vec::Vec;

pub mod cellfilter;
pub mod collate;
pub mod em;
pub mod pugutils;
pub mod quant;
pub mod schema;
pub mod utils;

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
    pub fn from_label_and_counter(corrected_bc_in: u64, num_remain: u64) -> CorrectedCBChunk {
        CorrectedCBChunk {
            remaining_records: num_remain,
            corrected_bc: corrected_bc_in,
            umis: Vec::<u64>::with_capacity(num_remain as usize),
            ref_offsets: Vec::<u32>::with_capacity(num_remain as usize),
            ref_ids: Vec::<u32>::with_capacity(3 * num_remain as usize),
        }
    }
}

#[derive(Copy, Clone)]
pub enum RADIntID {
    U8,
    U16,
    U32,
    U64,
}

impl RADIntID {
    pub fn bytes_for_type(&self) -> usize {
        match self {
            Self::U8 => std::mem::size_of::<u8>(),
            Self::U16 => std::mem::size_of::<u16>(),
            Self::U32 => std::mem::size_of::<u32>(),
            Self::U64 => std::mem::size_of::<u64>(),
        }
    }

    pub fn write_to<T>(&self, v: T, owriter: &mut BufWriter<File>) -> std::io::Result<()>
    where
        T: AsPrimitive<u8>
            + AsPrimitive<u16>
            + AsPrimitive<u32>
            + AsPrimitive<u64>
            + AsPrimitive<usize>
            + AsPrimitive<i8>
            + AsPrimitive<i16>
            + AsPrimitive<i32>
            + AsPrimitive<i64>
            + AsPrimitive<isize>
    {
        match self {
            Self::U8 => {
                let vo: u8 = v.as_();
                owriter.write_all(&vo.to_le_bytes())
            }
            Self::U16 => {
                let vo: u16 = v.as_();
                owriter.write_all(&vo.to_le_bytes())
            }
            Self::U32 => {
                let vo: u32 = v.as_();
                owriter.write_all(&vo.to_le_bytes())
            }
            Self::U64 => {
                let vo: u64 = v.as_();
                owriter.write_all(&vo.to_le_bytes())
            }
        }
    }
}

pub struct ChunkConfig {
    pub num_chunks: u64,
    pub bc_type: u8,
    pub umi_type: u8,
}

#[derive(Copy, Clone)]
pub enum RADType {
    BOOL,
    U8,
    U16,
    U32,
    U64,
    F32,
    F64,
}

pub fn decode_type_tag(type_id: u8) -> Option<RADType> {
    match type_id {
        0 => Some(RADType::BOOL),
        1 => Some(RADType::U8),
        2 => Some(RADType::U16),
        3 => Some(RADType::U32),
        4 => Some(RADType::U64),
        5 => Some(RADType::F32),
        6 => Some(RADType::F64),
        _ => None,
    }
}

pub fn decode_int_type_tag(type_id: u8) -> Option<RADIntID> {
    match type_id {
        1 => Some(RADIntID::U8),
        2 => Some(RADIntID::U16),
        3 => Some(RADIntID::U32),
        4 => Some(RADIntID::U64),
        _ => None,
    }
}

pub fn collect_records<T: Read>(
    reader: &mut BufReader<T>,
    chunk_config: &ChunkConfig,
    correct_map: &HashMap<u64, u64>,
    expected_ori: &Strand,
    output_cache: &mut HashMap<u64, CorrectedCBChunk>,
) {
    // NOTE: since the chunks are independent, this part could be multithreaded
    let bc_type = decode_int_type_tag(chunk_config.bc_type).expect("unknown barcode type id.");
    let umi_type = decode_int_type_tag(chunk_config.umi_type).expect("unknown barcode type id.");

    for _ in 0..(chunk_config.num_chunks as usize) {
        process_corrected_cb_chunk(
            reader,
            &bc_type,
            &umi_type,
            correct_map,
            expected_ori,
            output_cache,
        );
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
    pub fn is_empty(&self) -> bool {
        self.refs.is_empty()
    }

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
        expected_ori: &Strand,
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

            // fw if the leftmost bit is 1, otherwise rc
            let strand = if (v & utils::MASK_LOWER_31_U32) > 0 {
                Strand::Forward
            } else {
                Strand::Reverse
            };

            if expected_ori.same(&strand) {
                rec.refs.push(v & utils::MASK_TOP_BIT_U32);
            }
        }
        // make sure these are sorted in this step.
        quickersort::sort(&mut rec.refs[..]);
        rec
    }
}

pub fn process_corrected_cb_chunk<T: Read>(
    reader: &mut BufReader<T>,
    bct: &RADIntID,
    umit: &RADIntID,
    correct_map: &HashMap<u64, u64>,
    expected_ori: &Strand,
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
        let rr = ReadRecord::from_bytes_keep_ori(reader, &bct, &umit, expected_ori);
        // if this record had a correct or correctable barcode
        if let Some(corrected_id) = correct_map.get(&rr.bc) {
            if let Some(v) = output_cache.get_mut(corrected_id) {
                // update the corresponding corrected chunk entry
                v.remaining_records -= 1;
                // if there are no alignments in the record
                // (potentially b/c of orientation filtering)
                // then don't push info on to the vector.
                if rr.is_empty() {
                    continue;
                }
                v.umis.push(rr.umi);
                let ref_offset = v.ref_ids.len() as u32;
                v.ref_offsets.push(ref_offset);
                v.ref_ids.extend(&rr.refs);
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
    mut owriter: &mut BufWriter<File>,
    output_cache: &HashMap<u64, CorrectedCBChunk>,
    chunk_config: &ChunkConfig,
) {
    // NOTE: since the chunks are independent, this part could be multithreaded
    let bc_type = decode_int_type_tag(chunk_config.bc_type).expect("unknown barcode type id.");
    let umi_type = decode_int_type_tag(chunk_config.umi_type).expect("unknown barcode type id.");

    for (_bc, chunk) in output_cache.iter() {
        // number of bytes
        let mut nbytes: u32 = 0;
        let bytes_for_u32 = std::mem::size_of::<u32>();
        let bytes_for_bc = bc_type.bytes_for_type();
        let bytes_for_umi = umi_type.bytes_for_type();

        // for reference IDs in this chunk
        nbytes += (chunk.ref_ids.len() * bytes_for_u32) as u32;
        // umis
        nbytes += (chunk.umis.len() * bytes_for_umi) as u32;
        // barcodes
        nbytes += (chunk.umis.len() * bytes_for_bc) as u32;
        // num alignment fields
        nbytes += (chunk.umis.len() * bytes_for_u32) as u32;

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
            owriter
                .write_all(&num_aln.to_le_bytes())
                .expect("couldn't write output.");

            bc_type
                .write_to(*_bc, &mut owriter)
                .expect("couldn't write output.");
            umi_type
                .write_to(chunk.umis[i], &mut owriter)
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

    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>, bct: &RADIntID, umit: &RADIntID) -> Self {
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

pub fn permit_list_from_file(ifile: String, bclen: u16) -> Vec<u64> {
    let f = File::open(ifile).expect("couldn't open input barcode file.");
    let br = BufReader::new(f);
    let mut bc = Vec::<u64>::with_capacity(10_000);

    for l in br.lines() {
        let line = l.expect("couldn't read line from barcode file.");
        let mut bnk = BitNuclKmer::new(line.as_bytes(), bclen as u8, false);
        let (_, k, _) = bnk.next().expect("can't extract kmer");
        bc.push(k.0);
    }
    bc
}
