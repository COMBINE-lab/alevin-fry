/*
 * Copyright (c) 2020-2021 Rob Patro, Avi Srivastava, Hirak Sarkar, Dongze He, Mohsen Zakeri.
 *
 * This file is part of alevin-fry
 * (see https://github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

use std::fs::File;

pub const MASK_TOP_BIT_U32: u32 = 0x7FFFFFFF;
pub const MASK_LOWER_31_U32: u32 = 0x80000000;
pub const SPLICE_MASK_U32: u32 = 0xFFFFFFFE;

#[allow(dead_code)]
#[derive(Debug)]
pub struct InternalVersionInfo {
    major: u32,
    minor: u32,
    patch: u32,
}

impl InternalVersionInfo {
    pub fn from_string(vs: &str) -> Self {
        let versions: Vec<u32> = vs.split('.').map(|s| s.parse::<u32>().unwrap()).collect();
        assert_eq!(
            versions.len(),
            3,
            "The version string should be of the format x.y.z; it was {}",
            vs
        );
        Self {
            major: versions[0],
            minor: versions[1],
            patch: versions[2],
        }
    }

    pub fn is_compatible_with(&self, other: &InternalVersionInfo) -> Result<(), String> {
        if self.major == other.major && self.minor == other.minor {
            Ok(())
        } else {
            let s = format!(
                "version {:?} is incompatible with version {:?}",
                self, other
            );
            Err(s)
        }
    }
}

pub fn is_velo_mode(input_dir: String) -> bool {
    let parent = std::path::Path::new(&input_dir);
    // open the metadata file and read the json
    let meta_data_file = File::open(parent.join("generate_permit_list.json"))
        .expect("could not open the generate_permit_list.json file.");
    let mdata: serde_json::Value = serde_json::from_reader(meta_data_file)
        .expect("could not deseralize generate_permit_list.json");
    let vm = mdata.get("velo_mode");
    match vm {
        Some(v) => v.as_bool().unwrap_or(false),
        None => false,
    }
}

#[cfg(test)]
mod tests {
    use self::libradicl::utils::*;
    use crate as libradicl;

    #[test]
    fn test_version_info() {
        let vi = InternalVersionInfo::from_string("1.2.3");
        assert_eq!(
            vi,
            InternalVersionInfo{1, 2, 3}
        );
    }
}
