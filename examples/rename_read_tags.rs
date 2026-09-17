// Rewrite a RAD's prelude renaming read-level tags by an old=new map, copying
// chunk data verbatim. Roles/types are preserved. For validation only (e.g. to
// let a name-bridge consumer read a role-renamed collated RAD).
//
// Usage: rename_read_tags <in.rad> <out.rad> old1=new1 [old2=new2 ...]
use libradicl::header::RadPrelude;
use std::collections::HashMap;
use std::io::{BufReader, Seek, SeekFrom, Write};

fn main() -> anyhow::Result<()> {
    let mut a = std::env::args().skip(1);
    let inp = a.next().expect("in.rad");
    let outp = a.next().expect("out.rad");
    let map: HashMap<String, String> = a
        .map(|p| {
            let (o, n) = p.split_once('=').expect("expected old=new");
            (o.to_string(), n.to_string())
        })
        .collect();

    let mut br = BufReader::new(std::fs::File::open(&inp)?);
    let mut prelude = RadPrelude::from_bytes(&mut br)?;
    let desc_end = br.stream_position()?;

    for t in &mut prelude.read_tags.tags {
        if let Some(n) = map.get(&t.name) {
            t.name = n.clone();
        }
    }

    let mut out = std::io::BufWriter::new(std::fs::File::create(&outp)?);
    prelude.write(&mut out)?;
    let mut rest = std::fs::File::open(&inp)?;
    rest.seek(SeekFrom::Start(desc_end))?;
    std::io::copy(&mut rest, &mut out)?;
    out.flush()?;
    eprintln!("renamed {} read tag(s) -> {outp}", map.len());
    Ok(())
}
