use crate::quant::{ResolutionStrategy, SplicedAmbiguityModel};
use clap;
use std::path::{Path, PathBuf};

impl clap::ValueEnum for ResolutionStrategy {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            Self::Trivial,
            Self::CellRangerLike,
            Self::CellRangerLikeEm,
            Self::ParsimonyEm,
            Self::Parsimony,
            Self::ParsimonyGeneEm,
            Self::ParsimonyGene,
        ]
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        match self {
            Self::Trivial => Some(clap::builder::PossibleValue::new("trivial")),
            Self::CellRangerLike => Some(clap::builder::PossibleValue::new("cr-like")),
            Self::CellRangerLikeEm => Some(clap::builder::PossibleValue::new("cr-like-em")),
            Self::ParsimonyEm => Some(clap::builder::PossibleValue::new("parsimony-em")),
            Self::Parsimony => Some(clap::builder::PossibleValue::new("parsimony")),
            Self::ParsimonyGeneEm => Some(clap::builder::PossibleValue::new("parsimony-gene-em")),
            Self::ParsimonyGene => Some(clap::builder::PossibleValue::new("parsimony-gene")),
        }
    }
}

impl clap::ValueEnum for SplicedAmbiguityModel {
    fn value_variants<'a>() -> &'a [Self] {
        &[Self::PreferAmbiguity, Self::WinnerTakeAll]
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        match self {
            Self::PreferAmbiguity => Some(clap::builder::PossibleValue::new("prefer-ambig")),
            Self::WinnerTakeAll => Some(clap::builder::PossibleValue::new("winner-take-all")),
        }
    }
}

/// Checks if the path pointed to by v exists.  It can be
/// any valid entity (e.g. disk file, FIFO, directory, etc.).
/// If there is any issue with permissions or failure to properly
/// resolve symlinks, or if the path is wrong, it returns
/// an Err(String), else Ok(PathBuf).
pub fn pathbuf_file_exists_validator(v: &str) -> Result<PathBuf, String> {
    // NOTE: we explicitly *do not* check `is_file()` here
    // since we want to return true even if the path is to
    // a FIFO/named pipe.
    if !Path::new(v).exists() {
        Err(String::from("No valid file was found at this path."))
    } else {
        Ok(PathBuf::from(v))
    }
}

/// Checks if the path pointed to by v exists and is
/// a valid directory on disk.  If there is any issue
/// with permissions or failure to properly
/// resolve symlinks, or if the path is wrong, it returns
/// an Err(String), else Ok(PathBuf).
pub fn pathbuf_directory_exists_validator(v: &str) -> Result<PathBuf, String> {
    if !Path::new(v).is_dir() {
        Err(String::from("No valid directory was found at this path."))
    } else {
        Ok(PathBuf::from(v))
    }
}
