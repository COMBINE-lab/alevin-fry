// Standalone chunk-offset index builder for a collated RAD file.
//
// PROTOTYPE (measurement-quality, not upstream). One-time cost; in production
// this would be written at collate time. NOT counted in the timed quant
// comparison.
//
// Scans `map.collated.rad` once (seek-skip, not a full read), recording the
// absolute byte offset of every chunk header. Writes a sidecar
// `map.collated.rad.chunkidx`:
//     [num_chunks: u64 LE][offset_0 .. offset_num_chunks : u64 LE]  (num_chunks+1 offsets)
// offset_0 is the first chunk (right after prelude + file-tag map); the final
// offset equals the file length. The parallel reader in `quant` uses these to
// split the file into disjoint byte ranges on chunk boundaries.
//
// Usage: build_chunkidx <af_quant_dir | path/to/map.collated.rad>

use anyhow::{Context, ensure};
use libradicl::header::RadPrelude;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn main() -> anyhow::Result<()> {
    let arg = std::env::args()
        .nth(1)
        .context("usage: build_chunkidx <af_quant_dir | map.collated.rad>")?;
    let p = PathBuf::from(&arg);
    let rad_path: PathBuf = if p.is_dir() {
        p.join("map.collated.rad")
    } else {
        p
    };
    ensure!(
        rad_path.exists(),
        "collated RAD not found: {}",
        rad_path.display()
    );

    let t0 = Instant::now();

    // Parse the prelude + file-tag map exactly as do_quantify_dispatch does, so
    // the position we record for chunk 0 matches where the quant reader begins.
    let f = File::open(&rad_path)?;
    let file_len = f.metadata()?.len();
    let mut reader = BufReader::new(f);
    let prelude = RadPrelude::from_bytes(&mut reader)
        .context("could not parse RAD prelude; is this a collated RAD?")?;
    let _file_tag_map = prelude
        .file_tags
        .parse_tags_from_bytes(&mut reader)
        .map_err(|e| anyhow::anyhow!("could not parse file-tag map: {e}"))?;
    let num_chunks = prelude.hdr.num_chunks as usize;
    let offset0 = reader.stream_position()?;
    drop(reader);

    eprintln!(
        "num_chunks = {num_chunks}, first-chunk offset = {offset0}, file_len = {file_len}"
    );

    // Scan with a fresh handle, tracking the offset manually (offset += nbytes)
    // to avoid a stream_position syscall per chunk. Read only the 8-byte header
    // of each chunk, then seek past its payload.
    let mut f = File::open(&rad_path)?;
    f.seek(SeekFrom::Start(offset0))?;
    let mut offsets: Vec<u64> = Vec::with_capacity(num_chunks + 1);
    offsets.push(offset0);
    let mut off = offset0;
    let mut hdr = [0u8; 8];
    for c in 0..num_chunks {
        f.read_exact(&mut hdr)
            .with_context(|| format!("failed reading header of chunk {c}; RAD truncated?"))?;
        let nbytes = u32::from_le_bytes([hdr[0], hdr[1], hdr[2], hdr[3]]) as u64;
        ensure!(
            nbytes >= 8,
            "chunk {c} declares nbytes={nbytes} (< 8); corrupt RAD"
        );
        off += nbytes;
        offsets.push(off);
        // seek past the remaining payload (we already consumed the 8-byte header)
        if c + 1 < num_chunks {
            f.seek(SeekFrom::Current(nbytes as i64 - 8))?;
        }
    }

    // The load-bearing sanity check: "nbytes includes the 8-byte header" must be
    // exactly right, or every downstream byte range is wrong. If it is right,
    // the running offset lands precisely on EOF.
    ensure!(
        offsets.len() == num_chunks + 1,
        "expected {} offsets, got {}",
        num_chunks + 1,
        offsets.len()
    );
    ensure!(
        *offsets.last().unwrap() == file_len,
        "final offset {} != file length {} — chunk framing assumption wrong",
        offsets.last().unwrap(),
        file_len
    );

    let idx_path = sidecar_path(&rad_path);
    let mut out = File::create(&idx_path)?;
    out.write_all(&(num_chunks as u64).to_le_bytes())?;
    for o in &offsets {
        out.write_all(&o.to_le_bytes())?;
    }
    out.flush()?;

    eprintln!(
        "wrote {} ({} bytes) in {:.2}s; total chunk bytes = {}",
        idx_path.display(),
        (num_chunks + 2) * 8,
        t0.elapsed().as_secs_f64(),
        file_len - offset0
    );
    Ok(())
}

fn sidecar_path(rad: &Path) -> PathBuf {
    let mut s = rad.as_os_str().to_owned();
    s.push(".chunkidx");
    PathBuf::from(s)
}
