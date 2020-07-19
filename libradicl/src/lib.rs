// scroll now, explore nom later
extern crate scroll;
extern crate fasthash;

use scroll::{Pread}; 
use std::vec::{Vec};
use std::collections::HashMap;
use std::io::{Read, BufReader};
use std::fs::File;

mod utils;

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

// The below are currently hard-coded
// until we decide how to solve this 
// generally
#[derive(Debug)]
pub struct FileTags {
  pub bclen : u16,
  pub umilen : u16
}
#[derive(Debug)]
pub struct ReadRecord {
  pub bc : u64,
  pub umi: u64,
  pub dirs: Vec<bool>,
  pub refs: Vec<u32>
}
#[derive(Debug)]
pub struct Chunk {
  pub nbytes: u32,
  pub nrec: u32,
  pub reads: Vec<ReadRecord>
}

pub enum RADIntID {
  U8,
  U16,
  U32,
  U64
}

fn read_into_u64(reader: &mut BufReader<File>, rt : &RADIntID ) -> u64 {
  let mut rbuf = [0u8; 8];
  let v : u64;
  match rt {
    RADIntID::U8 => {
      reader.read_exact(&mut rbuf[0..1]).unwrap();
      v = rbuf.pread::<u8>(0).unwrap() as u64;
    },
    RADIntID::U16 => {
      reader.read_exact(&mut rbuf[0..2]).unwrap();
      v = rbuf.pread::<u16>(0).unwrap() as u64;
    },
    RADIntID::U32 => {
      reader.read_exact(&mut rbuf[0..4]).unwrap();
      v = rbuf.pread::<u32>(0).unwrap() as u64;
    },
    RADIntID::U64 => {
      reader.read_exact(&mut rbuf[0..8]).unwrap();
      v = rbuf.pread::<u64>(0).unwrap();
    }
  }
  v
}

impl ReadRecord {
  pub fn from_bytes(reader: &mut BufReader<File>, bct : &RADIntID, umit : &RADIntID) -> Self {
    let mut rbuf = [0u8; 255];
    
    reader.read_exact(&mut rbuf[0..4]).unwrap();
    let na = rbuf.pread::<u32>(0).unwrap();

    let bc = read_into_u64(reader, bct); 
    let umi = read_into_u64(reader, umit);

    let mut rec = Self {
      bc,
      umi,
      dirs : Vec::with_capacity(na as usize),
      refs : Vec::with_capacity(na as usize)
    };

    for _ in 0..(na as usize) {
      reader.read_exact(&mut rbuf[0..4]).unwrap();
      let v = rbuf.pread::<u32>(0).unwrap();
      let dir = (v & 0x80000000) != 0; 
      rec.dirs.push( dir );
      rec.refs.push( v & 0x7FFFFFFF );
    }

    rec
  }
}

impl Chunk {
  pub fn from_bytes(reader: &mut BufReader<File>, bct : RADIntID, umit : RADIntID) -> Self {
    let mut buf = [0u8;8];

    reader.read_exact(&mut buf).unwrap();
    let nbytes = buf.pread::<u32>(0).unwrap();
    let nrec = buf.pread::<u32>(4).unwrap();

    let mut c = Self {
      nbytes,
      nrec,
      reads : Vec::with_capacity(nrec as usize)
    };

    for _ in 0..(nrec as usize) {
      c.reads.push(ReadRecord::from_bytes(reader, &bct, &umit));
    }

    c
  }
}

impl FileTags {
  pub fn from_bytes(reader: &mut BufReader<File>) -> Self {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf).unwrap();

    Self {
      bclen : buf.pread::<u16>(0).unwrap(),
      umilen : buf.pread::<u16>(2).unwrap()
    }

  }
}

impl TagDesc {
  pub fn from_bytes(reader: &mut BufReader<File>) -> TagDesc {
    // space for the string length (1 byte)
    // the longest string possible (255 char)
    // and the typeid
    let mut buf = [0u8; 257];
    reader.read_exact(&mut buf[0..2]).unwrap();
    let str_len = buf.pread::<u16>(0).unwrap() as usize;
    
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

    let mut ts = TagSection{ tags: Vec::with_capacity(num_tags) };

    for _ in 0..num_tags {
      ts.tags.push( TagDesc::from_bytes(reader) );
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
      num_read += 1;
    }

    reader.read_exact(&mut buf[0..8]).unwrap();
    rh.num_chunks = buf.pread::<u64>(0).unwrap();
    rh
  }
}


pub fn update_barcode_hist(hist: &mut HashMap<u64, u64, fasthash::RandomState<fasthash::sea::Hash64>>, chunk: &Chunk) -> () {

  for r in &chunk.reads {
     *hist.entry(r.bc).or_insert(0) += 1;
  }

}