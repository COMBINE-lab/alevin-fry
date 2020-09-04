// Copyright 2020 Rob Patro, Avi Srivastava. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// scroll now, explore nom later
extern crate fasthash;
extern crate needletail;
extern crate num;
extern crate quickersort;
extern crate sce;
extern crate scroll;

use bio_types::strand::*;
use dashmap::DashMap;
use needletail::bitkmer::*;
use num::cast::AsPrimitive;
use scroll::Pread;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read};
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::vec::Vec;

pub mod cellfilter;
pub mod collate;
pub mod em;
pub mod exit_codes;
pub mod pugutils;
pub mod quant;
pub mod schema;
pub mod utils;
pub mod bulkrad;

// Name of the program, to be used in diagnostic messages.
static LIB_NAME: &str = "libradicl";

pub fn lib_name() -> &'static str {
    LIB_NAME
}

pub struct RADBulkHeader {
    pub is_paired: u8,
    pub ref_count: u64,
    pub ref_names: Vec<String>,
    pub ref_lens: Vec<u32>,
    pub num_chunks: u64,
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
    nrec: u32,
    data: Cursor<Vec<u8>>, /*,
                           umis: Vec<u64>,
                           ref_offsets: Vec<u32>,
                           ref_ids: Vec<u32>,
                           */
}

impl CorrectedCBChunk {
    pub fn from_label_and_counter(corrected_bc_in: u64, num_remain: u64) -> CorrectedCBChunk {
        let mut cc = CorrectedCBChunk {
            remaining_records: num_remain,
            corrected_bc: corrected_bc_in,
            nrec: 0u32,
            data: Cursor::new(Vec::<u8>::with_capacity((num_remain * 24) as usize))
            //umis: Vec::<u64>::with_capacity(num_remain as usize),
            //ref_offsets: Vec::<u32>::with_capacity(num_remain as usize),
            //ref_ids: Vec::<u32>::with_capacity(3 * num_remain as usize),
        };
        let dummy = 0u32;
        cc.data.write_all(&dummy.to_le_bytes()).unwrap();
        cc.data.write_all(&dummy.to_le_bytes()).unwrap();
        cc
    }
}

#[derive(Debug)]
pub enum AlignmentTag {
    RefId,
    Pos1, Pos2,
    Ori1, Ori2,
    Score1, Score2,
    Cigar1, Cigar2,
    NotListed,
}

impl std::str::FromStr for AlignmentTag {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match &*s.to_lowercase() {
            "refid"  => Ok(AlignmentTag::RefId),
            "leftpos" => Ok(AlignmentTag::Pos1),            
            "leftori" => Ok(AlignmentTag::Ori1),
            "leftscore" => Ok(AlignmentTag::Score1),
            "leftcigar" => Ok(AlignmentTag::Cigar1),
            "rightpos" => Ok(AlignmentTag::Pos2),
            "rightori" => Ok(AlignmentTag::Ori2),
            "rightscore" => Ok(AlignmentTag::Score2),
            "rightcigar" => Ok(AlignmentTag::Cigar2),

            _ => Err(format!("'{}' is not a bulk alignment standard type", s)),
        }
    }
}

#[derive(Debug)]
pub struct PEAlignment {
    pub ref_id: u32,
    pub pos1: u32,
    pub pos2: u32,
    pub ori1: bool,
    pub ori2: bool,
    pub score1: f32,
    pub score2: f32,
    pub cigar1: String,
    pub cigar2: String,
    // pub extra_tags: TagSection
}

#[derive(Debug)]
pub struct SEAlignment {
    pub ref_id: u32,
    pub pos: u32,
    pub dir: bool,
    pub score: f32,
    pub cigar: String,
    // pub extra_tags: TagSection
}

#[derive(Debug)]
pub struct BulkReadRecord {
    pub read_name: String,
    pub read_len: u32,
    pub alignments: Vec<PEAlignment>,
}

#[derive(Debug)]
pub struct BulkChunk {
    pub nbytes: u32,
    pub nrec: u32,
    pub reads: Vec<BulkReadRecord>,
}

#[derive(Copy, Clone)]
pub enum RADIntID {
    U8,
    U16,
    U32,
    U64,
}

pub trait PrimitiveInteger:
    AsPrimitive<u8>
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
}

impl<
        T: AsPrimitive<u8>
            + AsPrimitive<u16>
            + AsPrimitive<u32>
            + AsPrimitive<u64>
            + AsPrimitive<usize>
            + AsPrimitive<i8>
            + AsPrimitive<i16>
            + AsPrimitive<i32>
            + AsPrimitive<i64>
            + AsPrimitive<isize>,
    > PrimitiveInteger for T
{
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

    /// Based on the variant of the current enum, write the value `v`  
    /// out using `owrite`.  Here, `v` is bound to be some primitive
    /// integer type.  It is the responsibility of the caller to ensure
    /// that, if `v` is wider than the enum type on which this function
    /// is called, no important information is lost by discarding the higher
    /// order bits.
    pub fn write_to<T: PrimitiveInteger, U: Write>(
        &self,
        v: T,
        owriter: &mut U,
    ) -> std::io::Result<()> {
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
    STR,
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
        7 => Some(RADType::STR),
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

/*
pub fn collect_records<T: Read>(
    reader: &mut BufReader<T>,
    chunk_config: &ChunkConfig,
    correct_map: &HashMap<u64, u64>,
    expected_ori: &Strand,
    output_cache: &DashMap<u64, CorrectedCBChunk>,
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
*/

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

    pub fn from_bytes_record_header<T: Read>(
        reader: &mut BufReader<T>,
        bct: &RADIntID,
        umit: &RADIntID,
    ) -> (u64, u64, u32) {
        let mut rbuf = [0u8; 4];
        reader.read_exact(&mut rbuf).unwrap();
        let na = u32::from_le_bytes(rbuf); //.pread::<u32>(0).unwrap();
        let bc = read_into_u64(reader, bct);
        let umi = read_into_u64(reader, umit);
        (bc, umi, na)
    }

    pub fn from_bytes_with_header_keep_ori<T: Read>(
        reader: &mut BufReader<T>,
        bc: u64,
        umi: u64,
        na: u32,
        expected_ori: &Strand,
    ) -> Self {
        let mut rbuf = [0u8; 255];
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

impl BulkReadRecord {
    pub fn is_empty(&self) -> bool {
        self.alignments.is_empty()
    }

    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>, readTags: & TagSection, alignment_tags: & TagSection) -> Self {
        let mut rbuf = [0u8; 255];

        reader.read_exact(&mut rbuf[0..4]).unwrap();
        let na = rbuf.pread::<u32>(0).unwrap();

        let mut rec = Self {
            read_name: String::from(""),
            read_len: 0,
            alignments: Vec::with_capacity(na as usize),
        };
        
        for _ in 0..(na as usize) {
            rec.alignments.push(PEAlignment::from_bytes(reader, alignment_tags));
        }
        // for _ in 0..(na as usize) {
        //     reader.read_exact(&mut rbuf[0..4]).unwrap();
        //     let v = rbuf.pread::<u32>(0).unwrap();
        //     let dir = (v & 0x80000000) != 0;
        //     rec.dirs.push(dir);
        //     rec.refs.push(v & 0x7FFFFFFF);
        // }

        rec
    }
}

impl PEAlignment {
    pub fn is_empty(&self) -> bool {
        self.score1 == 0.0
    }

    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>, alignment_tags: & TagSection) -> Self {
        let mut rbuf = [0u8; 255];
        let mut rec = Self {
            ref_id: 0,
            pos1: 0, pos2: 0,
            ori1: true, ori2: true,
            score1: 0.0, score2: 0.0,
            cigar1: String::from(""), cigar2: String::from(""),
            // extra_tags: TagSection { tags: Vec::with_capacity(0)}
        };
        for idx in 0..alignment_tags.tags.len() {
            match alignment_tags.tags[idx].name.parse().unwrap_or(AlignmentTag::NotListed) {
                AlignmentTag::RefId => {
                    reader.read_exact(&mut rbuf[0..4]).unwrap();
                    rec.ref_id = rbuf.pread::<u32>(0).unwrap();
                },
                AlignmentTag::Pos1 => {
                    reader.read_exact(&mut rbuf[0..4]).unwrap();
                    rec.pos1 = rbuf.pread::<u32>(0).unwrap();
                },
                AlignmentTag::Pos2 => {
                    reader.read_exact(&mut rbuf[0..4]).unwrap();
                    rec.pos2 = rbuf.pread::<u32>(0).unwrap();
                },
                AlignmentTag::Ori1 => {
                    reader.read_exact(&mut rbuf[0..1]).unwrap();
                    rec.ori1 = rbuf.pread::<u8>(0).unwrap() == 1;
                },
                AlignmentTag::Ori2 => {
                    reader.read_exact(&mut rbuf[0..1]).unwrap();
                    rec.ori2 = rbuf.pread::<u8>(0).unwrap() == 1;
                },
                AlignmentTag::Score1 => {
                    reader.read_exact(&mut rbuf[0..4]).unwrap();
                    rec.score1 = rbuf.pread::<f32>(0).unwrap();
                },
                AlignmentTag::Score2 => {
                    reader.read_exact(&mut rbuf[0..4]).unwrap();
                    rec.score2 = rbuf.pread::<f32>(0).unwrap();
                },
                AlignmentTag::Cigar1 => {
                    reader.read_exact(&mut rbuf[0..1]).unwrap();
                    let string_len = rbuf.pread::<u8>(0).unwrap() as usize;
                    let mut strbuf = vec![0u8; string_len]; 
                    reader.read_exact(&mut strbuf).unwrap();
                    rec.cigar1 = String::from_utf8(strbuf).unwrap();
                },
                AlignmentTag::Cigar2 => {
                    reader.read_exact(&mut rbuf[0..1]).unwrap();
                    let string_len = rbuf.pread::<u8>(0).unwrap() as usize;
                    let mut strbuf = vec![0u8; string_len]; 
                    reader.read_exact(&mut strbuf).unwrap();
                    rec.cigar2 = String::from_utf8(strbuf).unwrap();
                },
                AlignmentTag::NotListed => {
                    // rec.extra_tags.tags.push(alignmentTags.tags[idx])
                }
            }
            // For any of the tags that are bulk alignment standard tags:alignmentTags
            
        }
        // for _ in 0..(na as usize) {
        //     reader.read_exact(&mut rbuf[0..4]).unwrap();
        //     let v = rbuf.pread::<u32>(0).unwrap();
        //     let dir = (v & 0x80000000) != 0;
        //     rec.dirs.push(dir);
        //     rec.refs.push(v & 0x7FFFFFFF);
        // }

        rec
    }
}

#[inline]
pub fn dump_chunk(v: &mut CorrectedCBChunk, owriter: &Mutex<BufWriter<File>>) {
    v.data.set_position(0);
    let nbytes = (v.data.get_ref().len()) as u32;
    let nrec = v.nrec;
    v.data.write_all(&nbytes.to_le_bytes()).unwrap();
    v.data.write_all(&nrec.to_le_bytes()).unwrap();
    owriter.lock().unwrap().write_all(v.data.get_ref()).unwrap();
}

pub fn collate_temporary_bucket<T: Read>(
    reader: &mut BufReader<T>,
    bct: &RADIntID,
    umit: &RADIntID,
    nchunks: u32,
    nrec: u32,
    output_cache: &mut HashMap<u64, CorrectedCBChunk>,
) {
    let mut tbuf = [0u8; 65536];
    // estimated average number of records per barcode
    // this is just for trying to pre-allocate buffers
    // right; should not affect correctness
    let est_num_rec = (nrec / nchunks) + 1;

    // for each record, read it
    for _ in 0..(nrec as usize) {
        // read the header of the record
        // we don't bother reading the whole thing here
        // because we will just copy later as need be
        let tup = ReadRecord::from_bytes_record_header(reader, &bct, &umit);

        // get the entry for this chunk, or create a new one
        let v = output_cache
            .entry(tup.0)
            .or_insert_with(|| CorrectedCBChunk::from_label_and_counter(tup.0, est_num_rec as u64));

        // keep track of the number of records we're writing
        (*v).nrec += 1;
        // write the num align
        let na = tup.2;
        (*v).data.write_all(&na.to_le_bytes()).unwrap();
        // write the corrected barcode
        bct.write_to(tup.0, &mut (*v).data).unwrap();
        umit.write_to(tup.1, &mut (*v).data).unwrap();
        // read the alignment records
        reader.read_exact(&mut tbuf[0..(4 * na as usize)]).unwrap();
        // write them
        (*v).data.write_all(&tbuf[..(4 * na as usize)]).unwrap();
    }
}

pub fn process_corrected_cb_chunk<T: Read>(
    reader: &mut BufReader<T>,
    bct: &RADIntID,
    umit: &RADIntID,
    correct_map: &HashMap<u64, u64>,
    expected_ori: &Strand,
    output_cache: &DashMap<u64, CorrectedCBChunk>,
    owriter: &Mutex<BufWriter<File>>,
) {
    let mut buf = [0u8; 8];
    let mut tbuf = [0u8; 65536];

    // get the number of bytes and records for
    // the next chunk
    reader.read_exact(&mut buf).unwrap();
    let _nbytes = buf.pread::<u32>(0).unwrap();
    let nrec = buf.pread::<u32>(4).unwrap();
    // for each record, read it
    for _ in 0..(nrec as usize) {
        let tup = ReadRecord::from_bytes_record_header(reader, &bct, &umit);
        //let rr = ReadRecord::from_bytes_keep_ori(reader, &bct, &umit, expected_ori);
        // if this record had a correct or correctable barcode
        if let Some(corrected_id) = correct_map.get(&tup.0) {
            let rr = ReadRecord::from_bytes_with_header_keep_ori(
                reader,
                tup.0,
                tup.1,
                tup.2,
                expected_ori,
            );

            if let Some(mut v) = output_cache.get_mut(corrected_id) {
                // update the corresponding corrected chunk entry
                v.remaining_records -= 1;
                let last_record = v.remaining_records == 0;
                // if there are no alignments in the record
                // (potentially b/c of orientation filtering)
                // then don't push info on to the vector.
                if rr.is_empty() {
                    if last_record {
                        dump_chunk(&mut v, owriter);
                    }
                    continue;
                }
                v.nrec += 1;
                let na = rr.refs.len() as u32;
                v.data.write_all(&na.to_le_bytes()).unwrap();
                bct.write_to(*corrected_id, &mut v.data).unwrap();
                umit.write_to(rr.umi, &mut v.data).unwrap();
                v.data.write_all(as_u8_slice(&rr.refs[..])).unwrap();
                if last_record {
                    dump_chunk(&mut v, owriter);
                }
            }
        } else {
            reader
                .read_exact(&mut tbuf[0..(4 * (tup.2 as usize))])
                .unwrap();
        }
    }
}

pub struct TempBucket {
    pub bucket_id: u32,
    pub bucket_writer: Arc<Mutex<BufWriter<File>>>,
    pub num_chunks: u32,
    pub num_records: u32,
    pub num_records_written: AtomicU32,
    pub num_bytes_written: AtomicU64,
}

impl TempBucket {
    pub fn from_id_and_parent(bucket_id: u32, parent: &std::path::Path) -> Self {
        TempBucket {
            bucket_id,
            bucket_writer: Arc::new(Mutex::new(BufWriter::new(
                File::create(parent.join(&format!("bucket_{}.tmp", bucket_id))).unwrap(),
            ))),
            num_chunks: 0u32,
            num_records: 0u32,
            num_records_written: AtomicU32::new(0u32),
            num_bytes_written: AtomicU64::new(0u64),
        }
    }
}

pub fn dump_corrected_cb_chunk_to_temp_file<T: Read>(
    reader: &mut BufReader<T>,
    bct: &RADIntID,
    umit: &RADIntID,
    correct_map: &HashMap<u64, u64>,
    expected_ori: &Strand,
    output_cache: &HashMap<u64, Arc<TempBucket>>,
    local_buffers: &mut [Cursor<Vec<u8>>],
) {
    let mut buf = [0u8; 8];
    let mut tbuf = vec![0u8; 65536];
    //let mut tcursor = Cursor::new(tbuf);
    //tcursor.set_position(0);

    // get the number of bytes and records for
    // the next chunk
    reader.read_exact(&mut buf).unwrap();
    let _nbytes = buf.pread::<u32>(0).unwrap();
    let nrec = buf.pread::<u32>(4).unwrap();

    let bc_bytes = bct.bytes_for_type();
    let umi_bytes = bct.bytes_for_type();
    let na_bytes = std::mem::size_of::<u32>();
    let target_id_bytes = std::mem::size_of::<u32>();

    // for each record, read it
    for _ in 0..(nrec as usize) {
        let tup = ReadRecord::from_bytes_record_header(reader, &bct, &umit);
        //let rr = ReadRecord::from_bytes_keep_ori(reader, &bct, &umit, expected_ori);
        // if this record had a correct or correctable barcode
        if let Some(corrected_id) = correct_map.get(&tup.0) {
            let rr = ReadRecord::from_bytes_with_header_keep_ori(
                reader,
                tup.0,
                tup.1,
                tup.2,
                expected_ori,
            );

            if rr.is_empty() {
                continue;
            }
            if let Some(v) = output_cache.get(corrected_id) {
                // if this is a valid barcode, then
                // write the corresponding entry to the
                // thread-local buffer for this bucket

                // update number of written records
                v.num_records_written.fetch_add(1, Ordering::SeqCst);
                let nb = (rr.refs.len() * target_id_bytes + na_bytes + bc_bytes + umi_bytes) as u64;

                v.num_bytes_written.fetch_add(nb, Ordering::SeqCst);
                let buffidx = v.bucket_id as usize;
                let bcursor = &mut local_buffers[buffidx];
                let na = rr.refs.len() as u32;
                bcursor.write_all(&na.to_le_bytes()).unwrap();
                bct.write_to(*corrected_id, bcursor).unwrap();
                umit.write_to(rr.umi, bcursor).unwrap();
                bcursor.write_all(as_u8_slice(&rr.refs[..])).unwrap();
                let len = bcursor.position() as usize;

                // if the thread-local buffer for this bucket is
                // greater than the flush size, then flush to file
                if len > 65536 {
                    let mut filebuf = v.bucket_writer.lock().unwrap();
                    filebuf
                        .write_all(&bcursor.get_ref()[0..len as usize])
                        .unwrap();
                    // and reset the local buffer cursor
                    bcursor.set_position(0);
                }
            }
        } else {
            reader
                .read_exact(&mut tbuf[0..(target_id_bytes * (tup.2 as usize))])
                .unwrap();
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

//pub fn dump_output_cache(
//    mut owriter: &mut BufWriter<File>,
//    output_cache: &DashMap<u64, CorrectedCBChunk>,
//    chunk_config: &ChunkConfig,
//) {
//    // NOTE: since the chunks are independent, this part could be multithreaded
//    let bc_type = decode_int_type_tag(chunk_config.bc_type).expect("unknown barcode type id.");
//    let umi_type = decode_int_type_tag(chunk_config.umi_type).expect("unknown barcode type id.");
//
//    for entry_ref in output_cache.iter() {
//        let _bc = entry_ref.key();
//        let chunk = entry_ref.value();
//        // number of bytes
//        let mut nbytes: u32 = 0;
//        let bytes_for_u32 = std::mem::size_of::<u32>();
//
//        let bytes_for_bc = bc_type.bytes_for_type();
//        let bytes_for_umi = umi_type.bytes_for_type();
//
//        // new
//        /*nbytes += chunk.data.get_ref().len() as u32;
//        owriter
//            .write_all(&nbytes.to_le_bytes())
//            .expect("couldn't write output.");
//        let nrec = chunk.nrec;
//        owriter
//            .write_all(&nrec.to_le_bytes())
//            .expect("couldn't write output.");
//        */
//        owriter
//            .write_all(&chunk.data.get_ref())
//            .expect("couldn't write output.");
//        // end new
//        /*
//        // for reference IDs in this chunk
//        nbytes += (chunk.ref_ids.len() * bytes_for_u32) as u32;
//        // umis
//        nbytes += (chunk.umis.len() * bytes_for_umi) as u32;
//        // barcodes
//        nbytes += (chunk.umis.len() * bytes_for_bc) as u32;
//        // num alignment fields
//        nbytes += (chunk.umis.len() * bytes_for_u32) as u32;
//
//        let nrec = chunk.umis.len() as u32;
//
//        owriter
//            .write_all(&nbytes.to_le_bytes())
//            .expect("couldn't write output.");
//        owriter
//            .write_all(&nrec.to_le_bytes())
//            .expect("couldn't write output.");
//
//        for i in 0..chunk.umis.len() {
//            let s = chunk.ref_offsets[i];
//            let e = if i == chunk.umis.len() - 1 {
//                chunk.ref_ids.len() as u32
//            } else {
//                chunk.ref_offsets[i + 1]
//            };
//
//            // num alignments
//            let num_aln = (e - s) as u32;
//            owriter
//                .write_all(&num_aln.to_le_bytes())
//                .expect("couldn't write output.");
//
//            bc_type
//                .write_to(*_bc, &mut owriter)
//                .expect("couldn't write output.");
//            umi_type
//                .write_to(chunk.umis[i], &mut owriter)
//                .expect("couldn't write output.");
//            owriter
//                .write_all(as_u8_slice(&chunk.ref_ids[(s as usize)..(e as usize)]))
//                .expect("couldn't write output.");
//        }
//        */
//    }
//}

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

impl BulkChunk {
    pub fn read_header<T: Read>(reader: &mut BufReader<T>) -> (u32, u32) {
        let mut buf = [0u8; 8];

        reader.read_exact(&mut buf).unwrap();
        let nbytes = buf.pread::<u32>(0).unwrap();
        let nrec = buf.pread::<u32>(4).unwrap();
        (nbytes, nrec)
    }

    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>, read_tags: & TagSection, alignment_tags: & TagSection) -> Result<Self, String> {
        let mut buf = [0u8; 8];

        if let Err(e) = reader.read_exact(&mut buf) {
            return Err(format!("EOF"));
        }
        
        let nbytes = buf.pread::<u32>(0).unwrap();
        let nrec = buf.pread::<u32>(4).unwrap();
        let mut c = Self {
            nbytes,
            nrec,
            reads: Vec::with_capacity(nrec as usize),
        };

        for _ in 0..(nrec as usize) {
            c.reads.push(BulkReadRecord::from_bytes(reader, read_tags, alignment_tags));
        }

        Ok(c)
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

impl RADBulkHeader {
    pub fn from_bytes<T: Read>(reader: &mut BufReader<T>) -> RADBulkHeader {
        let mut rh = RADBulkHeader {
            is_paired: 0,
            ref_count: 0,
            ref_names: vec![],
            ref_lens: vec![],
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
            reader.read_exact(&mut buf[0..1]).unwrap();
            let l: usize = buf.pread::<u8>(0).unwrap() as usize;
            reader.read_exact(&mut buf[0..l]).unwrap();
            rh.ref_names
                .push(std::str::from_utf8(&buf[0..l]).unwrap().to_string());
            reader.read_exact(&mut buf[0..4]).unwrap();
            rh.ref_lens.push(buf.pread::<u32>(0).unwrap());
            num_read += 1;
        }

        reader.read_exact(&mut buf[0..8]).unwrap();
        rh.num_chunks = buf.pread::<u64>(0).unwrap();
        rh
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
    expected_ori: &Strand,
) {
    match expected_ori {
        Strand::Unknown => {
            for r in &chunk.reads {
                *hist.entry(r.bc).or_insert(0) += 1;
            }
        }
        Strand::Forward => {
            for r in &chunk.reads {
                if r.dirs.iter().any(|&x| x) {
                    *hist.entry(r.bc).or_insert(0) += 1;
                }
            }
        }
        Strand::Reverse => {
            for r in &chunk.reads {
                if r.dirs.iter().any(|&x| !x) {
                    *hist.entry(r.bc).or_insert(0) += 1;
                }
            }
        }
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
