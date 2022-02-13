/// some (hopefully) generally useful I/O related utilities
use crossbeam_queue::ArrayQueue;
use indicatif::ProgressBar;
use scroll::Pwrite;

use libradicl::rad_types;
use std::collections::HashSet;
use std::io::Read;
use std::sync::Arc;

pub(crate) type MetaChunk = (usize, usize, u32, u32, Vec<u8>);

pub(crate) fn fill_work_queue<T: Read>(
    q: Arc<ArrayQueue<MetaChunk>>,
    mut br: T,
    num_chunks: usize,
    pbar: &ProgressBar,
) -> Result<(), Box<dyn std::error::Error>> {
    const BUFSIZE: usize = 524208;
    // the buffer that will hold our records
    let mut buf = vec![0u8; BUFSIZE];
    // the number of bytes currently packed into the chunk
    let mut cbytes = 0u32;
    // the number of records currently packed into the chunk
    let mut crec = 0u32;
    // the number of cells in the current chunk
    let mut cells_in_chunk = 0usize;
    // the offset of the first cell in this chunk
    let mut first_cell = 0usize;
    // if we had to expand the buffer already and should
    // forcibly push the current buffer onto the queue
    let mut force_push = false;
    // the number of bytes and records in the next chunk header
    let mut nbytes_chunk = 0u32;
    let mut nrec_chunk = 0u32;

    // we include the endpoint here because we will not actually
    // copy a chunk in the first iteration (since we have not yet
    // read the header, which comes at the end of the loop).
    for chunk_num in 0..=num_chunks {
        // in the first iteration we've not read a header yet
        // so we can't fill a chunk, otherwise we read the header
        // at the bottom of the previous iteration of this loop, and
        // we will fill in the buffer appropriately here.
        if chunk_num > 0 {
            // if the current cell (the cell whose header we read in the last iteration of
            // the loop) alone is too big for the buffer, then resize the buffer to be big enough
            if nbytes_chunk as usize > buf.len() {
                // if we had to resize the buffer to fit this cell, then make sure we push
                // immediately in the next round
                force_push = true;
                let chunk_resize = nbytes_chunk as usize + cbytes as usize;
                buf.resize(chunk_resize, 0);
            }

            // copy the data for the current chunk into the buffer
            let boffset = cbytes as usize;
            buf.pwrite::<u32>(nbytes_chunk, boffset)?;
            buf.pwrite::<u32>(nrec_chunk, boffset + 4)?;
            br.read_exact(&mut buf[(boffset + 8)..(boffset + nbytes_chunk as usize)])
                .unwrap();
            cells_in_chunk += 1;
            cbytes += nbytes_chunk;
            crec += nrec_chunk;
        }

        // in the last iteration of the loop, we will have read num_chunks headers already
        // and we are just filling up the buffer with the last cell, and there will be no more
        // headers left to read, so skip this
        if chunk_num < num_chunks {
            let (nc, nr) = rad_types::Chunk::read_header(&mut br);
            nbytes_chunk = nc;
            nrec_chunk = nr;
        }

        // determine if we should dump the current buffer to the work queue
        if force_push  // if we were told to push this chunk
	    || // or if adding the next cell to this chunk would exceed the buffer size
	    ((cbytes + nbytes_chunk) as usize > buf.len() && cells_in_chunk > 0)
	    || // of if this was the last chunk
	    chunk_num == num_chunks
        {
            // launch off these cells on the queue
            let mut bclone = (first_cell, cells_in_chunk, cbytes, crec, buf.clone());
            // keep trying until we can push this payload
            while let Err(t) = q.push(bclone) {
                bclone = t;
                // no point trying to push if the queue is full
                while q.is_full() {}
            }
            pbar.inc(cells_in_chunk as u64);

            // offset of the first cell in the next chunk
            first_cell += cells_in_chunk;
            // reset the counters
            cells_in_chunk = 0;
            cbytes = 0;
            crec = 0;
            buf.resize(BUFSIZE, 0);
            force_push = false;
        }
    }
    Ok(())
}

/// This function is the same as `fill_work_queue`, except that
/// when parsing the input file, it ignores (i.e. does not enqueue)
/// any cell whose barcode is not in `keep_set`.
pub(crate) fn fill_work_queue_filtered<T: Read>(
    keep_set: HashSet<u64, ahash::RandomState>,
    rl_tags: &rad_types::TagSection,
    q: Arc<ArrayQueue<MetaChunk>>,
    mut br: T,
    num_chunks: usize,
    pbar: &ProgressBar,
) -> Result<(), Box<dyn std::error::Error>> {
    let bct = rl_tags.tags[0].typeid;
    let umit = rl_tags.tags[1].typeid;
    let bc_type = rad_types::decode_int_type_tag(bct).expect("unsupported barcode type id.");
    let umi_type = rad_types::decode_int_type_tag(umit).expect("unsupported umi type id.");

    const BUFSIZE: usize = 524208;
    // the buffer that will hold our records
    let mut buf = vec![0u8; BUFSIZE];
    // the number of bytes currently packed into the chunk
    let mut cbytes = 0u32;
    // the number of records currently packed into the chunk
    let mut crec = 0u32;
    // the number of cells in the current chunk
    let mut cells_in_chunk = 0usize;
    // the offset of the first cell in this chunk
    let mut first_cell = 0usize;
    // if we had to expand the buffer already and should
    // forcibly push the current buffer onto the queue
    let mut force_push = false;
    // the number of bytes and records in the next chunk header
    let mut nbytes_chunk = 0u32;
    let mut nrec_chunk = 0u32;

    // we include the endpoint here because we will not actually
    // copy a chunk in the first iteration (since we have not yet
    // read the header, which comes at the end of the loop).
    for chunk_num in 0..=num_chunks {
        // in the first iteration we've not read a header yet
        // so we can't fill a chunk, otherwise we read the header
        // at the bottom of the previous iteration of this loop, and
        // we will fill in the buffer appropriately here.
        if chunk_num > 0 {
            // if the currenc cell (the cell whose header we read in the last iteration of
            // the loop) alone is too big for the buffer, than resize the buffer to be big enough
            if nbytes_chunk as usize > buf.len() {
                // if we had to resize the buffer to fit this cell, then make sure we push
                // immediately in the next round, unless we are skipping it's barcode
                force_push = true;
                let chunk_resize = nbytes_chunk as usize + cbytes as usize;
                buf.resize(chunk_resize, 0);
            }

            // copy the data for the current chunk into the buffer
            let boffset = cbytes as usize;
            buf.pwrite::<u32>(nbytes_chunk, boffset)?;
            buf.pwrite::<u32>(nrec_chunk, boffset + 4)?;
            br.read_exact(&mut buf[(boffset + 8)..(boffset + nbytes_chunk as usize)])
                .unwrap();
            // get the barcode for this chunk
            let (bc, _umi) =
                rad_types::Chunk::peek_record(&buf[boffset + 8..], &bc_type, &umi_type);
            if keep_set.contains(&bc) {
                cells_in_chunk += 1;
                cbytes += nbytes_chunk;
                crec += nrec_chunk;
            } else {
                // if we are skipping this cell, and it
                // triggered a force_push, then undo that
                force_push = false;
            }
        }

        // in the last iteration of the loop, we will have read num_chunks headers already
        // and we are just filling up the buffer with the last cell, and there will be no more
        // headers left to read, so skip this
        if chunk_num < num_chunks {
            let (nc, nr) = rad_types::Chunk::read_header(&mut br);
            nbytes_chunk = nc;
            nrec_chunk = nr;
        }

        // determine if we should dump the current buffer to the work queue
        if force_push  // if we were told to push this chunk
	    || // or if adding the next cell to this chunk would exceed the buffer size
	    ((cbytes + nbytes_chunk) as usize > buf.len() && cells_in_chunk > 0)
	    || // of if this was the last chunk
	    chunk_num == num_chunks
        {
            // launch off these cells on the queue
            let mut bclone = (first_cell, cells_in_chunk, cbytes, crec, buf.clone());
            // keep trying until we can push this payload
            while let Err(t) = q.push(bclone) {
                bclone = t;
                // no point trying to push if the queue is full
                while q.is_full() {}
            }
            pbar.inc(cells_in_chunk as u64);

            // offset of the first cell in the next chunk
            first_cell += cells_in_chunk;
            // reset the counters
            cells_in_chunk = 0;
            cbytes = 0;
            crec = 0;
            buf.resize(BUFSIZE, 0);
            force_push = false;
        }
    }
    Ok(())
}
