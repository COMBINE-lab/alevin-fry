// byteorder now, explore nom later
extern crate byteorder;

use byteorder::{ByteOrder, LittleEndian};
use std::vec::{Vec};
use std::io::{Read, BufReader};
use std::fs::File;

// Name of the program, to be used in diagnostic messages.
static LIB_NAME: &str = "libradicl";

pub fn lib_name() -> &'static str {
  LIB_NAME
}

pub struct RADHeader {
  pub is_paired: u8,
  pub ref_count: u64,
  pub ref_names: Vec<String>,
  pub num_chunks: u64
}

impl RADHeader {
  pub fn from_bytes(reader: &mut BufReader<File>) -> RADHeader {
    //let mut buffer = reader.fill_buf().unwrap();
    let mut rh = RADHeader{ is_paired : 0, ref_count : 0, ref_names : vec![], num_chunks : 0 };

    let mut buf = [0u8; 65536];
    reader.read_exact(&mut buf[0..1]).unwrap();
    rh.is_paired = buf[0];
    reader.read_exact(&mut buf[0..8]).unwrap();
    rh.ref_count = LittleEndian::read_u64(&buf);
    
    let mut num_read = 0u64;
    while num_read < rh.ref_count {
      reader.read_exact(&mut buf[0..2]).unwrap();
      let l : usize = LittleEndian::read_u16(&buf) as usize;
      reader.read_exact(&mut buf[0..l]).unwrap();
      rh.ref_names.push( std::str::from_utf8(&buf[0..l]).unwrap().to_string() );
      println!("{:?}",rh.ref_names[num_read as usize]);
      num_read += 1;
    }

    reader.read_exact(&mut buf[0..8]).unwrap();
    rh.num_chunks = LittleEndian::read_u64(&buf);
    rh
  }
}