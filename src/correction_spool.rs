/*
 * Copyright (c) 2020-2024 COMBINE-lab.
 *
 * This file is part of alevin-fry
 * (see https://www.github.com/COMBINE-lab/alevin-fry).
 *
 * License: 3-clause BSD, see https://opensource.org/licenses/BSD-3-Clause
 */

//! Private bounded spill/replay for ambiguous sample/cell pairs.

use anyhow::{Context, bail};
use snap::read::FrameDecoder;
use snap::write::FrameEncoder;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, ErrorKind, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_SPOOL: AtomicU64 = AtomicU64::new(0);

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SpoolStats {
    pub runs: u64,
    pub pairs: u64,
    pub encoded_pairs: u64,
    pub uncompressed_bytes: u64,
    pub compressed_bytes: u64,
}

/// A worker-local sequence of independently sorted, run-length encoded runs.
///
/// The file is removed on drop, including error unwinds. Runs do not need a
/// global merge because replay is purely additive.
pub(crate) struct DeferredPairSpool {
    path: PathBuf,
    writer: Option<FrameEncoder<BufWriter<File>>>,
    buffer: Vec<(u64, u64)>,
    max_pairs: usize,
    stats: SpoolStats,
    finished: bool,
}

impl DeferredPairSpool {
    pub fn new(directory: &Path, worker: usize, max_pairs: usize) -> anyhow::Result<Self> {
        std::fs::create_dir_all(directory).with_context(|| {
            format!(
                "could not create correction spool directory {}",
                directory.display()
            )
        })?;
        let (path, file) = create_unique_file(directory, worker)?;
        Ok(Self {
            path,
            writer: Some(FrameEncoder::new(BufWriter::new(file))),
            buffer: Vec::with_capacity(max_pairs.max(1)),
            max_pairs: max_pairs.max(1),
            stats: SpoolStats::default(),
            finished: false,
        })
    }

    pub fn push(&mut self, sample: u64, cell: u64) -> anyhow::Result<()> {
        if self.finished {
            bail!("attempted to append to a finished correction spool");
        }
        self.buffer.push((sample, cell));
        self.stats.pairs += 1;
        if self.buffer.len() >= self.max_pairs {
            self.flush_run()?;
        }
        Ok(())
    }

    pub fn finish(&mut self) -> anyhow::Result<()> {
        if self.finished {
            return Ok(());
        }
        self.flush_run()?;
        let encoder = self
            .writer
            .take()
            .expect("unfinished correction spool has a writer");
        let mut writer = encoder.into_inner().map_err(|error| error.into_error())?;
        writer.flush()?;
        drop(writer);
        self.stats.compressed_bytes = std::fs::metadata(&self.path)?.len();
        self.finished = true;
        Ok(())
    }

    pub fn replay(
        &mut self,
        mut consume: impl FnMut(u64, u64, u64) -> anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        self.finish()?;
        let file = File::open(&self.path)
            .with_context(|| format!("could not open correction spool {}", self.path.display()))?;
        let mut decoder = FrameDecoder::new(BufReader::new(file));
        loop {
            let Some(entries) = read_optional_u64(&mut decoder)? else {
                break;
            };
            for _ in 0..entries {
                let sample = read_u64(&mut decoder)?;
                let cell = read_u64(&mut decoder)?;
                let count = read_u64(&mut decoder)?;
                if count == 0 {
                    bail!("correction spool contains a zero-length run");
                }
                consume(sample, cell, count)?;
            }
        }
        Ok(())
    }

    pub fn stats(&self) -> SpoolStats {
        self.stats
    }

    #[cfg(test)]
    fn path(&self) -> &Path {
        &self.path
    }

    fn flush_run(&mut self) -> anyhow::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        self.buffer.sort_unstable();
        let encoded_pairs = 1 + self
            .buffer
            .windows(2)
            .filter(|pair| pair[0] != pair[1])
            .count();

        let writer = self
            .writer
            .as_mut()
            .expect("unfinished correction spool has a writer");
        writer.write_all(&(encoded_pairs as u64).to_le_bytes())?;
        let mut start = 0;
        while start < self.buffer.len() {
            let pair = self.buffer[start];
            let mut end = start + 1;
            while end < self.buffer.len() && self.buffer[end] == pair {
                end += 1;
            }
            writer.write_all(&pair.0.to_le_bytes())?;
            writer.write_all(&pair.1.to_le_bytes())?;
            writer.write_all(&((end - start) as u64).to_le_bytes())?;
            start = end;
        }
        self.stats.runs += 1;
        self.stats.encoded_pairs += encoded_pairs as u64;
        self.stats.uncompressed_bytes += 8 + 24 * encoded_pairs as u64;
        self.buffer.clear();
        Ok(())
    }
}

impl Drop for DeferredPairSpool {
    fn drop(&mut self) {
        self.writer.take();
        let _ = std::fs::remove_file(&self.path);
    }
}

fn create_unique_file(directory: &Path, worker: usize) -> anyhow::Result<(PathBuf, File)> {
    loop {
        let ordinal = NEXT_SPOOL.fetch_add(1, Ordering::Relaxed);
        let path = directory.join(format!(
            ".af-correction-{}-{worker}-{ordinal}.sz",
            std::process::id()
        ));
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(file) => return Ok((path, file)),
            Err(error) if error.kind() == ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(error).with_context(|| {
                    format!("could not create correction spool {}", path.display())
                });
            }
        }
    }
}

fn read_optional_u64(reader: &mut impl Read) -> anyhow::Result<Option<u64>> {
    let mut bytes = [0u8; 8];
    let mut filled = 0;
    while filled < bytes.len() {
        match reader.read(&mut bytes[filled..]) {
            Ok(0) if filled == 0 => return Ok(None),
            Ok(0) => bail!("correction spool is truncated"),
            Ok(read) => filled += read,
            Err(error) if error.kind() == ErrorKind::Interrupted => {}
            Err(error) => return Err(error.into()),
        }
    }
    Ok(Some(u64::from_le_bytes(bytes)))
}

fn read_u64(reader: &mut impl Read) -> anyhow::Result<u64> {
    read_optional_u64(reader)?.context("correction spool is truncated")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_one_and_many_runs_replay_identically_and_cleanup() {
        for (max_pairs, pairs, expected_runs) in [
            (4, Vec::new(), 0),
            (4, vec![(1, 2), (1, 2), (3, 4)], 1),
            (2, vec![(1, 2), (1, 2), (3, 4), (1, 2), (3, 4)], 3),
        ] {
            let dir = tempfile::tempdir().unwrap();
            let path;
            {
                let mut spool = DeferredPairSpool::new(dir.path(), 0, max_pairs).unwrap();
                path = spool.path().to_path_buf();
                for &(sample, cell) in &pairs {
                    spool.push(sample, cell).unwrap();
                }
                let mut replayed = std::collections::BTreeMap::new();
                spool
                    .replay(|sample, cell, count| {
                        *replayed.entry((sample, cell)).or_insert(0) += count;
                        Ok(())
                    })
                    .unwrap();
                let mut expected = std::collections::BTreeMap::new();
                for pair in pairs {
                    *expected.entry(pair).or_insert(0) += 1;
                }
                assert_eq!(replayed, expected);
                assert_eq!(spool.stats().runs, expected_runs);
                assert!(path.exists());
            }
            assert!(!path.exists());
        }
    }

    #[test]
    fn corrupt_spool_is_reported_and_still_cleaned_up() {
        let dir = tempfile::tempdir().unwrap();
        let path;
        {
            let mut spool = DeferredPairSpool::new(dir.path(), 0, 2).unwrap();
            spool.push(1, 2).unwrap();
            spool.finish().unwrap();
            path = spool.path().to_path_buf();
            std::fs::write(&path, b"not a snappy frame").unwrap();
            assert!(spool.replay(|_, _, _| Ok(())).is_err());
            assert!(path.exists());
        }
        assert!(!path.exists());
    }
}
