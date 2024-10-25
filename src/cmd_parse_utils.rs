use std::path::{Path, PathBuf};

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
