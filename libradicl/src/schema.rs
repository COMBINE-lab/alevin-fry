/*
 * Copyright (c) 2020-2021 Rob Patro, Avi Srivastava, Hirak Sarkar, Dongze He, Mohsen Zakeri.
 *
 * This file is part of alevin-fry
 * (see https://github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use bio_types;

#[allow(unused_imports)]
use bio_types::strand::Strand;

pub struct TempCellInfo {
    pub offset: u64,
    pub nbytes: u32,
    pub nrec: u32,
}

#[derive(Debug)]
pub struct ProtocolInfo {
    // TODO: only makes sense
    // for single-strand protocols
    // right now.  Expand to be generic.
    pub expected_ori: Strand,
}
