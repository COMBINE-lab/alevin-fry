"""Step 2: drop read names from the long-read RAD records.

The record keeps a one-byte read-name length after the UMI so that files written with names (older
`convert` runs) still parse: readers skip the announced bytes without allocating, writers always
write a length of 0, and the per-record size accounts for exactly one byte. The read-name String and
its per-record allocation are gone.
"""
import os
P = "/scratch1/zahra/source-code/libradicl/libradicl/src/record.rs"
s = open(P).read()

def rep(old, new, count=1):
    global s
    n = s.count(old); assert n == count, f"expected {count} of {old[:60]!r}, found {n}"
    s = s.replace(old, new)

# helper
rep("""/// Header information for an [AtacSeqReadRecord]
pub struct AtacSeqReadRecordHeader {""",
"""/// Skips the read-name field of a long-read record: a one-byte length followed by that many bytes.
/// Files written by this code store a length of 0; older files store the name, which is not needed.
#[inline]
pub fn skip_read_name<T: Read>(reader: &mut T) {
    let mut len_buf = [0u8; 1];
    reader.read_exact(&mut len_buf).expect("could not read the read-name length");
    let n = len_buf[0] as usize;
    if n > 0 {
        let mut skip = [0u8; 255];
        reader.read_exact(&mut skip[..n]).expect("could not skip the read name");
    }
}

/// Header information for an [AtacSeqReadRecord]
pub struct AtacSeqReadRecordHeader {""")

# collatable header: the name byte(s) belong to the header
rep("""        let bc = rad_io::read_into::<T, B>(reader, &context.bct);
        // NOTE: We likely will want to make the UMI generic as well
        let umi = rad_io::read_into_u64(reader, &context.umit);
        Ok(Self::CollatableRecordHeader { naln: na, bc, umi })
    }

    fn peek_collatable_header(
        buf: &[u8],
        ctx: &<Self as MappedRecord>::ParsingContext,
    ) -> anyhow::Result<Self::CollatableRecordHeader> {
        let na_size = mem::size_of::<u32>();
        let bc_size = ctx.bct.bytes_for_type();

        let na = buf.pread::<u32>(0).unwrap();

        let bc: B = match ctx.bct {
            RadIntId::U8 => NewU8(buf.pread::<u8>(na_size).unwrap()).into(),
            RadIntId::U16 => NewU16(buf.pread::<u16>(na_size).unwrap()).into(),
            RadIntId::U32 => NewU32(buf.pread::<u32>(na_size).unwrap()).into(),
            RadIntId::U64 => NewU64(buf.pread::<u64>(na_size).unwrap()).into(),
            RadIntId::U128 => NewU128(buf.pread::<u128>(na_size).unwrap()).into(),
            _ => panic!("signed barcode integer encodings are not supported"),
        };
        let umi = match ctx.umit {
            RadIntId::U8 => buf.pread::<u8>(na_size + bc_size).unwrap() as u64,
            RadIntId::U16 => buf.pread::<u16>(na_size + bc_size).unwrap() as u64,
            RadIntId::U32 => buf.pread::<u32>(na_size + bc_size).unwrap() as u64,
            RadIntId::U64 => buf.pread::<u64>(na_size + bc_size).unwrap(),
            RadIntId::U128 => panic!("u128 is currently not supported as a umi type"),
            _ => panic!("signed umi integer encodings are not supported"),
        };
        Ok(Self::CollatableRecordHeader { naln: na, bc, umi })
    }

    fn set_collate_key(&mut self, k: B) {
        self.bc = k;
    }
    fn collate_key(&self) -> B {
        self.bc
    }

    #[inline]
    fn from_bytes_with_header_retain_ori<T: Read>(
        reader: &mut T,
        hdr: &mut Self::CollatableRecordHeader,
        _ctx: &<Self as MappedRecord>::ParsingContext,
        expected_ori: &MappedFragmentOrientation,
    ) -> Self {
        let na = hdr.naln;
        let bc = hdr.bc;
        let umi = hdr.umi;
        let mut rbuf = [0u8; 255];

        let mut rec = Self {
            bc,
            umi,
            dirs: Vec::with_capacity(na as usize),
            refs: Vec::with_capacity(na as usize),
            as_scores: Vec::with_capacity(na as usize),
            starts: Vec::with_capacity(na as usize),
            ends: Vec::with_capacity(na as usize),
            tlens: Vec::with_capacity(na as usize),
        };
""",
"""        let bc = rad_io::read_into::<T, B>(reader, &context.bct);
        // NOTE: We likely will want to make the UMI generic as well
        let umi = rad_io::read_into_u64(reader, &context.umit);
        // the read-name field follows the UMI and is part of the header
        skip_read_name(reader);
        Ok(Self::CollatableRecordHeader { naln: na, bc, umi })
    }

    fn peek_collatable_header(
        buf: &[u8],
        ctx: &<Self as MappedRecord>::ParsingContext,
    ) -> anyhow::Result<Self::CollatableRecordHeader> {
        let na_size = mem::size_of::<u32>();
        let bc_size = ctx.bct.bytes_for_type();

        let na = buf.pread::<u32>(0).unwrap();

        let bc: B = match ctx.bct {
            RadIntId::U8 => NewU8(buf.pread::<u8>(na_size).unwrap()).into(),
            RadIntId::U16 => NewU16(buf.pread::<u16>(na_size).unwrap()).into(),
            RadIntId::U32 => NewU32(buf.pread::<u32>(na_size).unwrap()).into(),
            RadIntId::U64 => NewU64(buf.pread::<u64>(na_size).unwrap()).into(),
            RadIntId::U128 => NewU128(buf.pread::<u128>(na_size).unwrap()).into(),
            _ => panic!("signed barcode integer encodings are not supported"),
        };
        let umi = match ctx.umit {
            RadIntId::U8 => buf.pread::<u8>(na_size + bc_size).unwrap() as u64,
            RadIntId::U16 => buf.pread::<u16>(na_size + bc_size).unwrap() as u64,
            RadIntId::U32 => buf.pread::<u32>(na_size + bc_size).unwrap() as u64,
            RadIntId::U64 => buf.pread::<u64>(na_size + bc_size).unwrap(),
            RadIntId::U128 => panic!("u128 is currently not supported as a umi type"),
            _ => panic!("signed umi integer encodings are not supported"),
        };
        Ok(Self::CollatableRecordHeader { naln: na, bc, umi })
    }

    fn set_collate_key(&mut self, k: B) {
        self.bc = k;
    }
    fn collate_key(&self) -> B {
        self.bc
    }

    #[inline]
    fn from_bytes_with_header_retain_ori<T: Read>(
        reader: &mut T,
        hdr: &mut Self::CollatableRecordHeader,
        _ctx: &<Self as MappedRecord>::ParsingContext,
        expected_ori: &MappedFragmentOrientation,
    ) -> Self {
        // the header (including the read-name field) was consumed by from_bytes_collatable_header
        let na = hdr.naln;
        let bc = hdr.bc;
        let umi = hdr.umi;
        let mut rbuf = [0u8; 255];

        let mut rec = Self {
            bc,
            umi,
            dirs: Vec::with_capacity(na as usize),
            refs: Vec::with_capacity(na as usize),
            as_scores: Vec::with_capacity(na as usize),
            starts: Vec::with_capacity(na as usize),
            ends: Vec::with_capacity(na as usize),
            tlens: Vec::with_capacity(na as usize),
        };
""")

# header writer: the zero-length name
rep("""        ctx.umit
            .write_to(self.umi, writer)
            .context("couldn't write umi field for record")?;
        Ok(())
    }
}

// ==== ATAC seq read""",
"""        ctx.umit
            .write_to(self.umi, writer)
            .context("couldn't write umi field for record")?;
        // the read-name field: always empty
        writer
            .write_all(&[0u8])
            .context("couldn't write the read-name length for record")?;
        Ok(())
    }
}

// ==== ATAC seq read""")

# the full-record parser and writer
rep("""        let (bc, umi, na) = Self::from_bytes_record_header(reader, &ctx.bct, &ctx.umit);
        let mut rec = Self {
            bc,
            umi,
            dirs: Vec::with_capacity(na as usize),
            refs: Vec::with_capacity(na as usize),
            as_scores: Vec::with_capacity(na as usize),
            starts: Vec::with_capacity(na as usize),
            ends: Vec::with_capacity(na as usize),
            tlens: Vec::with_capacity(na as usize),
        };
""",
"""        let (bc, umi, na) = Self::from_bytes_record_header(reader, &ctx.bct, &ctx.umit);
        skip_read_name(reader);
        let mut rec = Self {
            bc,
            umi,
            dirs: Vec::with_capacity(na as usize),
            refs: Vec::with_capacity(na as usize),
            as_scores: Vec::with_capacity(na as usize),
            starts: Vec::with_capacity(na as usize),
            ends: Vec::with_capacity(na as usize),
            tlens: Vec::with_capacity(na as usize),
        };
""")
rep("""        ctx.umit
            .write_to(self.umi, writer)
            .context("couldn't write umi field for record")?;

        for i in 0..(na as usize) {
            let ref_idx = self.refs[i];
            let dir = self.dirs[i];
            let as_i32 = self.as_scores[i];""",
"""        ctx.umit
            .write_to(self.umi, writer)
            .context("couldn't write umi field for record")?;
        // the read-name field: always empty
        writer
            .write_all(&[0u8])
            .context("couldn't write the read-name length for record")?;

        for i in 0..(na as usize) {
            let ref_idx = self.refs[i];
            let dir = self.dirs[i];
            let as_i32 = self.as_scores[i];""")

# size accounting
rep("""impl<B: ConvertiblePrimitiveInteger> KnownSize for ScLongReadRecordT<B> {
    fn nbytes(na: u32, ctx: &<Self as MappedRecord>::ParsingContext) -> usize {
        // for na field
        std::mem::size_of::<u32>() +
        // for barcode type
        ctx.bct.bytes_for_type() +
        // for the umi
        ctx.umit.bytes_for_type() +
        // for each alignment a 
        (na as usize * Self::nbytes_aln(ctx))
    }""",
"""impl<B: ConvertiblePrimitiveInteger> KnownSize for ScLongReadRecordT<B> {
    fn nbytes(na: u32, ctx: &<Self as MappedRecord>::ParsingContext) -> usize {
        // for na field
        std::mem::size_of::<u32>() +
        // for barcode type
        ctx.bct.bytes_for_type() +
        // for the umi
        ctx.umit.bytes_for_type() +
        // the read-name length byte (the name itself is never written)
        1 +
        // for each alignment a 
        (na as usize * Self::nbytes_aln(ctx))
    }""")
open(P, "w").write(s)
print("step 2: libradicl patched")
