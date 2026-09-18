// Standalone exercise of libradicl::chunk::validate_first_chunk_layout on a raw
// map.rad: parse the prelude, then run the field-completeness self-check on the
// first chunk. Confirms the check does not false-positive on valid formats.
use anyhow::Context;
use libradicl::header::RadPrelude;
use std::io::{BufReader, Seek};

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .context("usage: check_layout <map.rad>")?;
    let mut br = BufReader::new(std::fs::File::open(&path)?);
    let prelude = RadPrelude::from_bytes(&mut br)?;
    // Consume file-tag values so the reader sits at the first chunk.
    let _ = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    let pos = br.stream_position()?;
    println!(
        "{path}: major={} minor={} num_chunks={} read_tags={} aln_tags={} first_chunk@{pos}",
        prelude.hdr.version.major(),
        prelude.hdr.version.minor(),
        prelude.hdr.num_chunks,
        prelude.read_tags.tags.len(),
        prelude.aln_tags.tags.len(),
    );
    match libradicl::chunk::validate_first_chunk_layout(
        &mut br,
        &prelude.read_tags,
        &prelude.aln_tags,
    ) {
        Ok(()) => println!("  field-completeness: OK"),
        Err(e) => {
            println!("  field-completeness: FAILED -> {e:#}");
            std::process::exit(1);
        }
    }
    Ok(())
}
