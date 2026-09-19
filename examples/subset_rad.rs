// Carve a small RAD holding the first K chunks of an input RAD, rewriting the
// header's num_chunks. Chunks are copied verbatim (each on disk is exactly its
// declared `nbytes` including the 8-byte [nbytes][nrec] header). For building
// small fixtures from large real RADs.
//
// Usage: subset_rad <in.rad> <out.rad> <num_chunks>
use libradicl::header::RadPrelude;
use std::io::{BufReader, Read, Seek, SeekFrom, Write};

fn main() -> anyhow::Result<()> {
    let mut a = std::env::args().skip(1);
    let inp = a.next().expect("in.rad");
    let outp = a.next().expect("out.rad");
    let k: u64 = a.next().expect("num_chunks").parse()?;

    let mut br = BufReader::new(std::fs::File::open(&inp)?);
    let mut prelude = RadPrelude::from_bytes(&mut br)?;
    let _ = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    let first_chunk = br.stream_position()?;
    let orig_chunks = prelude.hdr.num_chunks;
    anyhow::ensure!(k <= orig_chunks, "requested {k} > available {orig_chunks}");
    prelude.hdr.num_chunks = k;

    // prelude + file-tag values are the bytes [0, first_chunk); copy them, but
    // with the rewritten num_chunks. Easiest: write prelude (descriptors) then
    // copy the file-tag-value bytes [desc_end, first_chunk). Find desc_end.
    let mut br2 = BufReader::new(std::fs::File::open(&inp)?);
    let _ = RadPrelude::from_bytes(&mut br2)?;
    let desc_end = br2.stream_position()?;

    let mut out = std::io::BufWriter::new(std::fs::File::create(&outp)?);
    prelude.write(&mut out)?;
    // copy file-tag values verbatim
    let mut src = std::fs::File::open(&inp)?;
    src.seek(SeekFrom::Start(desc_end))?;
    let mut ftv = vec![0u8; (first_chunk - desc_end) as usize];
    src.read_exact(&mut ftv)?;
    out.write_all(&ftv)?;

    // copy K chunks
    src.seek(SeekFrom::Start(first_chunk))?;
    for _ in 0..k {
        let mut hb = [0u8; 8];
        src.read_exact(&mut hb)?;
        let nbytes = u32::from_le_bytes(hb[0..4].try_into().unwrap()) as usize;
        out.write_all(&hb)?;
        let mut payload = vec![0u8; nbytes - 8];
        src.read_exact(&mut payload)?;
        out.write_all(&payload)?;
    }
    out.flush()?;
    eprintln!("wrote {k} chunks (of {orig_chunks}) -> {outp}");
    Ok(())
}
