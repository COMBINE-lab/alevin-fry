/*
 * Copyright (c) 2020-2021 Rob Patro, Avi Srivastava, Hirak Sarkar, Dongze He, Mohsen Zakeri.
 *
 * This file is part of alevin-fry
 * (see https://github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

// scroll now, explore nom later

use crate as libradicl;

use self::libradicl::rad_types::{CorrectedCbChunk, RadIntId, ReadRecord};
use self::libradicl::schema::TempCellInfo;
#[allow(unused_imports)]
use ahash::{AHasher, RandomState};
use bio_types::strand::*;
use dashmap::DashMap;
use scroll::Pread;
use serde::{Deserialize, Serialize};

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom};
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::vec::Vec;

pub mod exit_codes;
pub mod rad_types;
pub mod schema;
pub mod utils;

// Name of the program, to be used in diagnostic messages.
static LIB_NAME: &str = "libradicl";

pub fn lib_name() -> &'static str {
    LIB_NAME
}

#[derive(Serialize, Deserialize, Debug)]
#[allow(dead_code)]
pub struct BarcodeLookupMap {
    pub barcodes: Vec<u64>,
    //pub counts: Vec<usize>,
    offsets: Vec<usize>,
    bclen: u32,
    prefix_len: u32,
    suffix_len: u32,
}

impl BarcodeLookupMap {
    pub fn new(mut kv: Vec<u64>, bclen: u32) -> BarcodeLookupMap {
        let prefix_len = ((bclen + 1) / 2) as u64;
        let suffix_len = bclen - prefix_len as u32;

        let _prefix_bits = 2 * prefix_len;
        let suffix_bits = 2 * suffix_len;

        kv.sort_unstable();

        let pref_mask = ((4usize.pow(prefix_len as u32) - 1) as u64) << (suffix_bits);
        let mut offsets = vec![0; 4usize.pow(prefix_len as u32) + 1];
        let mut prev_ind = 0xFFFF;

        for (n, &v) in kv.iter().enumerate() {
            let ind = ((v & pref_mask) >> (suffix_bits)) as usize;
            if ind != prev_ind {
                for item in offsets.iter_mut().take(ind).skip(prev_ind + 1) {
                    *item = n;
                }
                offsets[ind] = n;
                prev_ind = ind;
            }
        }
        for item in offsets.iter_mut().skip(prev_ind + 1) {
            *item = kv.len();
        }

        //let nbc = kv.len();
        BarcodeLookupMap {
            barcodes: kv,
            //counts: vec![0usize; nbc],
            offsets,
            bclen,
            prefix_len: prefix_len as u32,
            suffix_len,
        }
    }

    #[allow(dead_code)]
    pub fn barcode_for_idx(&self, idx: usize) -> u64 {
        self.barcodes[idx]
    }

    pub fn find_exact(&self, query: u64) -> Option<usize> {
        let mut ret: Option<usize> = None;

        // extract the prefix we will use to search
        let suffix_bits = 2 * self.suffix_len;
        let query_pref = query >> suffix_bits;

        // the range of entries having query_pref as their prefix
        let qrange = std::ops::Range {
            start: self.offsets[query_pref as usize],
            end: self.offsets[(query_pref + 1) as usize],
        };

        let qs = qrange.start as usize;

        // if we can, then we return the found barcode and that there was 1 best hit
        if let Ok(res) = self.barcodes[qrange].binary_search(&query) {
            ret = Some(qs + res);
        }
        ret
    }

    /// The find function searches for the barcode `query` in the
    /// BarcodeLookupMap.  It returns a tuple `(Option<usize>, usize)` where
    /// the first element is either Some(usize) or None.  If
    /// Some(usize) is returned, this is the *index* of a matching/neighboring barcode
    /// if None is returned, then no match was found.  The second element is either
    /// 0, 1 or 2.  If 0, no match was found; if 1 a unique match was found, if 2
    /// then 2 or more equally good matches were found.
    ///
    /// The parameter `try_exact` controls whether a an exact search is performed
    /// or not.  If this parameter is true, an exact search is performed before
    /// a neighbor search.  Otherwise, the exact search is skipped.
    pub fn find_neighbors(&self, query: u64, try_exact: bool) -> (Option<usize>, usize) {
        let mut ret: Option<usize> = None;

        // extract the prefix we will use to search
        let pref_bits = 2 * self.prefix_len;
        let suffix_bits = 2 * self.suffix_len;
        let mut query_pref = query >> suffix_bits;
        let mut num_neighbors = 0usize;

        // the range of entries having query_pref as their prefix
        let qrange = std::ops::Range {
            start: self.offsets[query_pref as usize],
            end: self.offsets[(query_pref + 1) as usize],
        };

        let qs = qrange.start as usize;

        if try_exact {
            // first, we try to find exactly.
            // if we can, then we return the found barcode and that there was 1 best hit
            if let Ok(res) = self.barcodes[qrange.clone()].binary_search(&query) {
                ret = Some(qs + res);
                num_neighbors += 1;
                return (ret, num_neighbors);
            }
        }

        // othwerwise, we fall back to the 1 mismatch search
        // NOTE: We stop here as soon as we find at most 2 neighbors
        // for the query.  Thus, we only distinguish between the
        // the cases where the query has 1 neighbor, or 2 or more neighbors.

        // if we match the prefix exactly, we will look for possible matches
        // that are 1 mismatch off in the suffix.
        if !(std::ops::Range::<usize>::is_empty(&qrange)) {
            // the initial offset of suffixes for this prefix
            let qs = qrange.start as usize;

            // for each position in the suffix
            for i in (0..suffix_bits).step_by(2) {
                let bit_mask = 3 << (i);

                // for each nucleotide
                for nmod in 1..4 {
                    let nucl = 0x3 & ((query >> i) + nmod);
                    let nquery = (query & (!bit_mask)) | (nucl << i);

                    if let Ok(res) = self.barcodes[qrange.clone()].binary_search(&nquery) {
                        ret = Some(qs + res);
                        num_neighbors += 1;
                        if num_neighbors >= 2 {
                            return (ret, num_neighbors);
                        }
                    }
                }
            }
        }

        {
            // if we get here we've had either 0 or 1 matches holding the prefix fixed
            // so we will now hold the suffix fixed and consider possible mutations of the prefix.

            // for each position in the prefix
            for i in (suffix_bits..(suffix_bits + pref_bits)).step_by(2) {
                let bit_mask = 3 << i;

                // for each nucleotide
                for nmod in 1..4 {
                    let nucl = 0x3 & ((query >> i) + nmod);
                    let nquery = (query & (!bit_mask)) | (nucl << i);

                    query_pref = nquery >> suffix_bits;

                    let qrange = std::ops::Range {
                        start: self.offsets[query_pref as usize],
                        end: self.offsets[(query_pref + 1) as usize],
                    };
                    let qs = qrange.start as usize;
                    if let Ok(res) = self.barcodes[qrange].binary_search(&nquery) {
                        ret = Some(qs + res);
                        num_neighbors += 1;
                        if num_neighbors >= 2 {
                            return (ret, num_neighbors);
                        }
                    }
                }
            }
        }

        (ret, num_neighbors)
    }
}

impl CorrectedCbChunk {
    pub fn from_label_and_counter(corrected_bc_in: u64, num_remain: u32) -> CorrectedCbChunk {
        let mut cc = CorrectedCbChunk {
            remaining_records: num_remain,
            corrected_bc: corrected_bc_in,
            nrec: 0u32,
            data: Cursor::new(Vec::<u8>::with_capacity((num_remain * 24) as usize)), //umis: Vec::<u64>::with_capacity(num_remain as usize),
                                                                                     //ref_offsets: Vec::<u32>::with_capacity(num_remain as usize),
                                                                                     //ref_ids: Vec::<u32>::with_capacity(3 * num_remain as usize),
        };
        let dummy = 0u32;
        cc.data.write_all(&dummy.to_le_bytes()).unwrap();
        cc.data.write_all(&dummy.to_le_bytes()).unwrap();
        cc
    }
}

#[derive(Debug, Clone)]
pub struct GlobalEqCellList {
    cell_ids: Vec<usize>,
    count: u32,
}

impl GlobalEqCellList {
    pub fn from_umi_and_count(bc_mer: usize, count: u32) -> GlobalEqCellList {
        let mut cc = GlobalEqCellList {
            cell_ids: Vec::new(),
            count: 0,
        };
        cc.cell_ids.push(bc_mer);
        cc.count += count;
        cc
    }
    pub fn add_element(&mut self, bc_mer: usize, count: u32) {
        self.cell_ids.push(bc_mer);
        self.count += count;
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

#[inline]
pub fn dump_chunk(v: &mut CorrectedCbChunk, owriter: &Mutex<BufWriter<File>>) {
    v.data.set_position(0);
    let nbytes = (v.data.get_ref().len()) as u32;
    let nrec = v.nrec;
    v.data.write_all(&nbytes.to_le_bytes()).unwrap();
    v.data.write_all(&nrec.to_le_bytes()).unwrap();
    owriter.lock().unwrap().write_all(v.data.get_ref()).unwrap();
}

/// Given a `BufReader<T>` from which to read a set of records that
/// should reside in the same collated bucket, this function will
/// collate the records by cell barcode, filling them into a chunk of
/// memory exactly as they will reside on disk.  If `compress` is true
/// the collated chunk will be compressed, and then the result will be
/// written to the output guarded by `owriter`.
pub fn collate_temporary_bucket_twopass<T: Read + Seek, U: Write>(
    reader: &mut BufReader<T>,
    bct: &RadIntId,
    umit: &RadIntId,
    nrec: u32,
    owriter: &Mutex<U>,
    compress: bool,
    cb_byte_map: &mut HashMap<u64, TempCellInfo, ahash::RandomState>,
) -> usize {
    let mut tbuf = vec![0u8; 65536];
    let mut total_bytes = 0usize;
    let header_size = 2 * std::mem::size_of::<u32>() as u64;
    let size_of_u32 = std::mem::size_of::<u32>();
    let size_of_bc = bct.bytes_for_type();
    let size_of_umi = umit.bytes_for_type();

    let calc_record_bytes = |num_aln: usize| -> usize {
        size_of_u32 + size_of_bc + size_of_umi + (size_of_u32 * num_aln)
    };

    // read each record
    for _ in 0..(nrec as usize) {
        // read the header of the record
        // we don't bother reading the whole thing here
        // because we will just copy later as need be
        let tup = ReadRecord::from_bytes_record_header(reader, bct, umit);

        // get the entry for this chunk, or create a new one
        let v = cb_byte_map.entry(tup.0).or_insert(TempCellInfo {
            offset: header_size,
            nbytes: header_size as u32,
            nrec: 0_u32,
        });

        // read the alignment records from the input file
        let na = tup.2 as usize;
        let req_size = size_of_u32 * na;
        if tbuf.len() < req_size {
            tbuf.resize(req_size, 0);
        }
        reader.read_exact(&mut tbuf[0..(size_of_u32 * na)]).unwrap();
        // compute the total number of bytes this record requires
        let nbytes = calc_record_bytes(na);
        (*v).offset += nbytes as u64;
        (*v).nbytes += nbytes as u32;
        (*v).nrec += 1;
        total_bytes += nbytes as usize;
    }

    // each cell will have a header (8 bytes each)
    total_bytes += cb_byte_map.len() * header_size as usize;
    let mut output_buffer = Cursor::new(vec![0u8; total_bytes]);

    // loop over all distinct cell barcodes, write their
    // corresponding chunk header, and compute what the
    // offset in `output_buffer` is where the corresponding
    // records should start.
    let mut next_offset = 0u64;
    for (_, v) in cb_byte_map.iter_mut() {
        // jump to the position where this chunk should start
        // and write the header
        output_buffer.set_position(next_offset);
        let cell_bytes = (*v).nbytes as u32;
        let cell_rec = (*v).nrec as u32;
        output_buffer.write_all(&cell_bytes.to_le_bytes()).unwrap();
        output_buffer.write_all(&cell_rec.to_le_bytes()).unwrap();
        // where we will start writing records for this cell
        (*v).offset = output_buffer.position();
        // the number of bytes allocated to this chunk
        let nbytes = (*v).nbytes as u64;
        // the next record will start after this one
        next_offset += nbytes;
    }

    // now each key points to where we should write the next record for the CB
    // reset the input pointer
    reader
        .get_mut()
        .seek(SeekFrom::Start(0))
        .expect("could not get read pointer.");

    // for each record, read it
    for _ in 0..(nrec as usize) {
        // read the header of the record
        // we don't bother reading the whole thing here
        // because we will just copy later as need be
        let tup = ReadRecord::from_bytes_record_header(reader, bct, umit);

        // get the entry for this chunk, or create a new one
        if let Some(v) = cb_byte_map.get_mut(&tup.0) {
            output_buffer.set_position(v.offset);

            // write the num align
            let na = tup.2 as usize;
            let nau32 = na as u32;
            output_buffer.write_all(&nau32.to_le_bytes()).unwrap();

            // write the corrected barcode
            bct.write_to(tup.0, &mut output_buffer).unwrap();
            umit.write_to(tup.1, &mut output_buffer).unwrap();

            // read the alignment records
            reader
                .read_exact(&mut tbuf[0..(size_of_u32 as usize * na)])
                .unwrap();
            // write them
            output_buffer
                .write_all(&tbuf[..(size_of_u32 as usize * na)])
                .unwrap();

            (*v).offset = output_buffer.position();
        } else {
            panic!("should not have any barcodes we can't find");
        }
    }

    output_buffer.set_position(0);

    if compress {
        // compress the contents of output_buffer to compressed_output
        let mut compressed_output =
            snap::write::FrameEncoder::new(Cursor::new(Vec::<u8>::with_capacity(total_bytes)));
        compressed_output
            .write_all(output_buffer.get_ref())
            .expect("could not compress the output chunk.");

        output_buffer = compressed_output
            .into_inner()
            .expect("couldn't unwrap the FrameEncoder.");
        output_buffer.set_position(0);
    }

    owriter
        .lock()
        .unwrap()
        .write_all(output_buffer.get_ref())
        .unwrap();

    cb_byte_map.len()
}

pub fn collate_temporary_bucket<T: Read>(
    reader: &mut T,
    bct: &RadIntId,
    umit: &RadIntId,
    _nchunks: u32,
    nrec: u32,
    output_cache: &mut HashMap<u64, CorrectedCbChunk, ahash::RandomState>,
) {
    let mut tbuf = [0u8; 65536];
    // estimated average number of records per barcode
    // this is just for trying to pre-allocate buffers
    // right; should not affect correctness
    let est_num_rec = 1; //(nrec / nchunks) + 1;

    // for each record, read it
    for _ in 0..(nrec as usize) {
        // read the header of the record
        // we don't bother reading the whole thing here
        // because we will just copy later as need be
        let tup = ReadRecord::from_bytes_record_header(reader, bct, umit);

        // get the entry for this chunk, or create a new one
        let v = output_cache
            .entry(tup.0)
            .or_insert_with(|| CorrectedCbChunk::from_label_and_counter(tup.0, est_num_rec));

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
    reader: &mut T,
    bct: &RadIntId,
    umit: &RadIntId,
    correct_map: &HashMap<u64, u64>,
    expected_ori: &Strand,
    output_cache: &DashMap<u64, CorrectedCbChunk>,
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
        let tup = ReadRecord::from_bytes_record_header(reader, bct, umit);
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

/// Represents a temporary bucket of barcodes whose records will
/// be written together and then collated later in memory.
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
            bucket_writer: Arc::new(Mutex::new(BufWriter::with_capacity(
                4096_usize,
                File::create(parent.join(&format!("bucket_{}.tmp", bucket_id))).unwrap(),
            ))),
            num_chunks: 0u32,
            num_records: 0u32,
            num_records_written: AtomicU32::new(0u32),
            num_bytes_written: AtomicU64::new(0u64),
        }
    }
}

/// Read an input chunk from `reader` and write the
/// resulting records to the corresponding in-memory
/// buffers `local_buffers`.  As soon as any buffer
/// reaches `flush_limit`, flush the buffer by writing
/// it to the `output_cache`.
#[allow(clippy::too_many_arguments)]
pub fn dump_corrected_cb_chunk_to_temp_file<T: Read>(
    reader: &mut BufReader<T>,
    bct: &RadIntId,
    umit: &RadIntId,
    correct_map: &HashMap<u64, u64>,
    expected_ori: &Strand,
    output_cache: &HashMap<u64, Arc<TempBucket>>,
    local_buffers: &mut [Cursor<&mut [u8]>],
    flush_limit: usize,
) {
    let mut buf = [0u8; 8];
    let mut tbuf = vec![0u8; 4096];
    //let mut tcursor = Cursor::new(tbuf);
    //tcursor.set_position(0);

    // get the number of bytes and records for
    // the next chunk
    reader.read_exact(&mut buf).unwrap();
    let _nbytes = buf.pread::<u32>(0).unwrap();
    let nrec = buf.pread::<u32>(4).unwrap();

    let bc_bytes = bct.bytes_for_type();
    let umi_bytes = umit.bytes_for_type();
    let na_bytes = std::mem::size_of::<u32>();
    let target_id_bytes = std::mem::size_of::<u32>();

    // for each record, read it
    for _ in 0..(nrec as usize) {
        let tup = ReadRecord::from_bytes_record_header(reader, bct, umit);

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

                // the total number of bytes this record will take
                let nb = (rr.refs.len() * target_id_bytes + na_bytes + bc_bytes + umi_bytes) as u64;

                // the buffer index for this corrected barcode
                let buffidx = v.bucket_id as usize;
                // the current cursor for this buffer
                let bcursor = &mut local_buffers[buffidx];
                // the current position of the cursor
                let len = bcursor.position() as usize;

                // if writing the next record (nb bytes) will put us over
                // the flush size for the thread-local buffer for this bucket
                // then first flush the buffer to file.
                if len + nb as usize >= flush_limit {
                    let mut filebuf = v.bucket_writer.lock().unwrap();
                    filebuf
                        .write_all(&bcursor.get_ref()[0..len as usize])
                        .unwrap();
                    // and reset the local buffer cursor
                    bcursor.set_position(0);
                }

                // now, write the record to the buffer
                let na = rr.refs.len() as u32;
                bcursor.write_all(&na.to_le_bytes()).unwrap();
                bct.write_to(*corrected_id, bcursor).unwrap();
                umit.write_to(rr.umi, bcursor).unwrap();
                bcursor.write_all(as_u8_slice(&rr.refs[..])).unwrap();

                // update number of written records
                v.num_records_written.fetch_add(1, Ordering::SeqCst);
                // update number of written bytes
                v.num_bytes_written.fetch_add(nb, Ordering::SeqCst);
            }
        } else {
            // in this branch, we don't have access to a correct barcode for
            // what we observed, so we need to discard the remaining part of
            // the record.
            let req_len = target_id_bytes * (tup.2 as usize);
            let do_resize = req_len > tbuf.len();

            if do_resize {
                tbuf.resize(req_len, 0);
            }

            reader
                .read_exact(&mut tbuf[0..(target_id_bytes * (tup.2 as usize))])
                .unwrap();

            if do_resize {
                tbuf.resize(4096, 0);
                tbuf.shrink_to_fit();
            }
        }
    }
}

pub fn as_u8_slice(v: &[u32]) -> &[u8] {
    unsafe {
        std::slice::from_raw_parts(
            v.as_ptr() as *const u8,
            v.len() * std::mem::size_of::<u32>(),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::BarcodeLookupMap;

    #[test]
    fn test_barcode_lookup_map() {
        let barcode_sv_even = vec![
            b"AACC", b"AAGG", b"CAGT", b"CATT", b"GACC", b"GATA", b"TCAG", b"TCGT",
        ];
        let barcode_sv_odd = vec![
            b"AACCA", b"AAGGC", b"CAGTA", b"CATTG", b"GACCG", b"GATAC", b"TCAGA", b"TCGTG",
        ];

        let mut barcode_even = Vec::with_capacity(barcode_sv_even.len());
        for b in barcode_sv_even.clone() {
            if let Some((_, km, _)) =
                needletail::bitkmer::BitNuclKmer::new(&b[..], b.len() as u8, false).next()
            {
                barcode_even.push(km.0);
            }
        }
        let mut barcode_odd = Vec::with_capacity(barcode_sv_odd.len());
        for b in barcode_sv_odd.clone() {
            if let Some((_, km, _)) =
                needletail::bitkmer::BitNuclKmer::new(&b[..], b.len() as u8, false).next()
            {
                barcode_odd.push(km.0);
            }
        }

        let me = BarcodeLookupMap::new(barcode_even, 4);
        let mo = BarcodeLookupMap::new(barcode_odd, 5);

        let x = b"CAGA";
        if let Some((_, et, _)) =
            needletail::bitkmer::BitNuclKmer::new(&x[..], x.len() as u8, false).next()
        {
            assert_eq!((Some(2), 1), me.find_neighbors(et.0, false));
        }
        let x = b"CAATG";
        if let Some((_, et, _)) =
            needletail::bitkmer::BitNuclKmer::new(&x[..], x.len() as u8, false).next()
        {
            assert_eq!((Some(3), 1), mo.find_neighbors(et.0, false));
        }
    }
}
