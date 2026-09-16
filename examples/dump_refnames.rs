// Print the RAD header's reference names, one per line.
use libradicl::header::RadHeader;
use std::io::BufReader;

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: dump_refnames <rad>");
    let mut br = BufReader::new(std::fs::File::open(&path)?);
    let hdr = RadHeader::from_bytes(&mut br)?;
    for name in &hdr.ref_names {
        println!("{name}");
    }
    eprintln!("refs={}", hdr.ref_names.len());
    Ok(())
}
