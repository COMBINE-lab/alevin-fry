// scroll now, explore nom later
extern crate scroll;

use scroll::{ctx, Pread, LE, Endian};

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

pub struct TagDesc {
  pub name: String,
  pub typeid: u8
}


pub struct TagSection {
  pub tags: Vec<TagDesc>
}


impl TagDesc {
  pub fn from_bytes(reader: &mut BufReader<File>) -> TagDesc {
    // space for the string length (1 byte)
    // the longest string possible (255 char)
    // and the typeid
    let mut buf = [0u8; 257];
    reader.read_exact(&mut buf[0..2]).unwrap();
    let str_len = buf.pread::<u16>(0).unwrap() as usize;
    println!("STR LEN : {:?}", str_len);
    // read str_len + 1 to get the type id that follows the string
    reader.read_exact(&mut buf[0..str_len+1]).unwrap();
    TagDesc {
      name : std::str::from_utf8(&buf[0..str_len]).unwrap().to_string(),
      typeid : buf.pread(str_len).unwrap()
    }
  }
}


impl TagSection {
  pub fn from_bytes(reader: &mut BufReader<File>) -> TagSection {
    let mut buf = [0u8; 2];
    reader.read_exact(&mut buf).unwrap();
    let num_tags = buf.pread::<u16>(0).unwrap() as usize;
    println!("Will read {:?}", num_tags);

    let mut ts = TagSection{ tags: Vec::with_capacity(num_tags) };

    for _ in 0..num_tags {
      ts.tags.push( TagDesc::from_bytes(reader) );
      println!("TAG: {:?}", ts.tags.last().unwrap().name );
    }

    ts
  }
}

impl RADHeader {
  pub fn from_bytes(reader: &mut BufReader<File>) -> RADHeader {
    
    let mut rh = RADHeader{ is_paired : 0, ref_count : 0, ref_names : vec![], num_chunks : 0 };

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
      let l : usize = buf.pread::<u16>(0).unwrap() as usize;
      reader.read_exact(&mut buf[0..l]).unwrap();
      rh.ref_names.push( std::str::from_utf8(&buf[0..l]).unwrap().to_string() );
      //println!("{:?}",rh.ref_names[num_read as usize]);
      num_read += 1;
    }

    reader.read_exact(&mut buf[0..8]).unwrap();
    rh.num_chunks = buf.pread::<u64>(0).unwrap();
    rh
  }
}