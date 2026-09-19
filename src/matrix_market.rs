//! Bounded, buffered output for MatrixMarket coordinate matrices.

use std::fmt::Write as FmtWrite;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::sync::Mutex;

pub(crate) const MATRIX_BUFFER_CAPACITY: usize = 256 * 1024;
const NNZ_FIELD_WIDTH: usize = 20;
// A decimal usize needs at most usize::BITS characters; 64 more cover any
// displayed f32 and the separators. Reserve this slack once so the entry
// crossing the batch threshold does not double the String's allocation.
const ENTRY_CAPACITY: usize = 2 * usize::BITS as usize + 64;
const BANNER: &str = "%%MatrixMarket matrix coordinate real general\n% written by sprs\n";

enum WriterState {
    Open,
    Finished,
    Failed,
}

/// A seekable matrix stream whose nonzero count is filled in on completion.
///
/// Coordinates can arrive in any order. Each batch is checked against the
/// declared dimensions before its bytes are appended. An I/O failure makes
/// the stream unusable: retrying a partially written batch could duplicate
/// coordinates or overwrite data after a failed seek.
pub(crate) struct MatrixMarketWriter<W: Write + Seek> {
    writer: BufWriter<W>,
    rows: usize,
    cols: usize,
    nnz: u64,
    nnz_offset: u64,
    state: WriterState,
}

impl<W: Write + Seek> MatrixMarketWriter<W> {
    pub(crate) fn new(inner: W, rows: usize, cols: usize) -> io::Result<Self> {
        let mut writer = BufWriter::with_capacity(MATRIX_BUFFER_CAPACITY, inner);
        writer.write_all(BANNER.as_bytes())?;
        let size_prefix = format!("{rows} {cols} ");
        writer.write_all(size_prefix.as_bytes())?;
        let nnz_offset = (BANNER.len() + size_prefix.len()) as u64;
        writer.write_all(&[b' '; NNZ_FIELD_WIDTH])?;
        writer.write_all(b"\n")?;
        Ok(Self {
            writer,
            rows,
            cols,
            nnz: 0,
            nnz_offset,
            state: WriterState::Open,
        })
    }

    pub(crate) fn append(&mut self, batch: &MatrixBatch) -> io::Result<()> {
        self.require_open()?;
        if batch.nnz == 0 {
            return Ok(());
        }
        if batch.max_row >= self.rows || batch.max_col >= self.cols {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "matrix batch contains zero-based row {} or column {} outside dimensions {} x {}",
                    batch.max_row, batch.max_col, self.rows, self.cols
                ),
            ));
        }
        let nnz = self.nnz.checked_add(batch.nnz).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "matrix entry count overflow")
        })?;
        if let Err(error) = self.writer.write_all(batch.text.as_bytes()) {
            self.state = WriterState::Failed;
            return Err(error);
        }
        self.nnz = nnz;
        Ok(())
    }

    /// Flush all entries, patch the reserved count field, and flush that patch.
    /// Successful finalization is idempotent; appends are subsequently rejected.
    pub(crate) fn finish(&mut self) -> io::Result<u64> {
        if matches!(self.state, WriterState::Finished) {
            return Ok(self.nnz);
        }
        self.require_open()?;
        let result = (|| {
            // BufWriter::seek flushes the body before moving the underlying
            // stream. Keep formatting inside the buffer so each padding space
            // in the fixed-width field does not become a filesystem write.
            self.writer.seek(SeekFrom::Start(self.nnz_offset))?;
            write!(self.writer, "{:<width$}", self.nnz, width = NNZ_FIELD_WIDTH)?;
            self.writer.flush()
        })();
        match result {
            Ok(()) => {
                self.state = WriterState::Finished;
                Ok(self.nnz)
            }
            Err(error) => {
                self.state = WriterState::Failed;
                Err(error)
            }
        }
    }

    fn require_open(&self) -> io::Result<()> {
        match self.state {
            WriterState::Open => Ok(()),
            WriterState::Finished => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "matrix stream has already been finalized",
            )),
            WriterState::Failed => Err(io::Error::other(
                "matrix stream cannot continue after an I/O failure",
            )),
        }
    }
}

/// Worker-local text staging, flushed while traversing entries rather than at
/// cell boundaries. Even a single very large cell cannot grow this buffer
/// beyond the threshold plus one formatted entry. Empty batches allocate no
/// staging storage. Coordinates passed to `push` are zero-based.
#[derive(Default)]
pub(crate) struct MatrixBatch {
    text: String,
    nnz: u64,
    max_row: usize,
    max_col: usize,
}

impl MatrixBatch {
    pub(crate) fn push<W: Write + Seek>(
        &mut self,
        row: usize,
        col: usize,
        val: f32,
        output: &Mutex<MatrixMarketWriter<W>>,
    ) -> io::Result<()> {
        let one_based_row = row.checked_add(1).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "matrix row index overflow")
        })?;
        let one_based_col = col.checked_add(1).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "matrix column index overflow")
        })?;
        let nnz = self.nnz.checked_add(1).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "matrix batch count overflow")
        })?;
        if self.text.capacity() == 0 {
            self.text
                .try_reserve_exact(MATRIX_BUFFER_CAPACITY + ENTRY_CAPACITY)
                .map_err(io::Error::other)?;
        }
        writeln!(self.text, "{one_based_row} {one_based_col} {val}").map_err(io::Error::other)?;
        self.nnz = nnz;
        self.max_row = self.max_row.max(row);
        self.max_col = self.max_col.max(col);
        if self.text.len() >= MATRIX_BUFFER_CAPACITY {
            self.flush(output)?;
        }
        Ok(())
    }

    pub(crate) fn flush<W: Write + Seek>(
        &mut self,
        output: &Mutex<MatrixMarketWriter<W>>,
    ) -> io::Result<()> {
        if self.nnz == 0 {
            return Ok(());
        }
        output
            .lock()
            .map_err(|_| io::Error::other("matrix output lock is poisoned"))?
            .append(self)?;
        // Preserve the batch on failure, and reuse its allocation on success.
        self.text.clear();
        self.nnz = 0;
        self.max_row = 0;
        self.max_col = 0;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[derive(Clone, Copy)]
    enum Fault {
        Write,
        ShortWrite,
        Seek,
        Flush,
    }

    #[derive(Default)]
    struct TestStream {
        data: Cursor<Vec<u8>>,
        write_lengths: Vec<usize>,
        fault: Option<Fault>,
    }

    impl Write for TestStream {
        fn write(&mut self, data: &[u8]) -> io::Result<usize> {
            if matches!(self.fault, Some(Fault::Write)) {
                return Err(io::Error::other("injected write failure"));
            }
            let length = if matches!(self.fault, Some(Fault::ShortWrite)) {
                self.fault = Some(Fault::Write);
                data.len().min(7)
            } else {
                data.len()
            };
            self.write_lengths.push(length);
            self.data.write(&data[..length])
        }

        fn flush(&mut self) -> io::Result<()> {
            if matches!(self.fault, Some(Fault::Flush)) {
                Err(io::Error::other("injected flush failure"))
            } else {
                Ok(())
            }
        }
    }

    impl Seek for TestStream {
        fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
            if matches!(self.fault, Some(Fault::Seek)) {
                Err(io::Error::other("injected seek failure"))
            } else {
                self.data.seek(pos)
            }
        }
    }

    // (nrows, ncols, nnz) header triple plus the parsed (row, col, value) entries.
    type ParsedMatrix = ((usize, usize, usize), Vec<(usize, usize, f32)>);

    fn parse_matrix(bytes: &[u8]) -> ParsedMatrix {
        let text = std::str::from_utf8(bytes).unwrap();
        assert!(text.starts_with(BANNER));
        let mut lines = text.lines().filter(|line| !line.starts_with('%'));
        let dims: Vec<usize> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|field| field.parse().unwrap())
            .collect();
        assert_eq!(dims.len(), 3);
        let entries: Vec<_> = lines
            .map(|line| {
                let fields: Vec<_> = line.split_whitespace().collect();
                assert_eq!(fields.len(), 3);
                let row = fields[0].parse::<usize>().unwrap();
                let col = fields[1].parse::<usize>().unwrap();
                assert!((1..=dims[0]).contains(&row));
                assert!((1..=dims[1]).contains(&col));
                (row, col, fields[2].parse::<f32>().unwrap())
            })
            .collect();
        assert_eq!(entries.len(), dims[2]);
        ((dims[0], dims[1], dims[2]), entries)
    }

    #[test]
    fn round_trips_empty_and_fractional_matrices() {
        let mut empty = MatrixMarketWriter::new(Cursor::new(Vec::new()), 0, 5).unwrap();
        assert_eq!(empty.finish().unwrap(), 0);
        assert_eq!(empty.finish().unwrap(), 0);
        assert_eq!(parse_matrix(empty.writer.get_ref().get_ref()).0, (0, 5, 0));
        assert!(empty.append(&MatrixBatch::default()).is_err());

        let output = Mutex::new(MatrixMarketWriter::new(Cursor::new(Vec::new()), 3, 5).unwrap());
        let mut batch = MatrixBatch::default();
        batch.push(2, 4, 0.25, &output).unwrap();
        batch.push(0, 0, 0.0, &output).unwrap();
        batch.flush(&output).unwrap();
        batch.push(1, 2, -3.5, &output).unwrap();
        batch.flush(&output).unwrap();
        let mut writer = output.lock().unwrap();
        assert_eq!(writer.finish().unwrap(), 3);
        let (dims, entries) = parse_matrix(writer.writer.get_ref().get_ref());
        assert_eq!(dims, (3, 5, 3));
        assert_eq!(entries, vec![(3, 5, 0.25), (1, 1, 0.0), (2, 3, -3.5)]);
    }

    #[test]
    fn one_large_cell_uses_multiple_bounded_batches() {
        const ENTRIES: usize = 50_000;
        let output =
            Mutex::new(MatrixMarketWriter::new(Cursor::new(Vec::new()), 1, ENTRIES).unwrap());
        let values = [
            0.25,
            f32::MAX,
            f32::MIN_POSITIVE,
            f32::from_bits(1),
            -f32::from_bits(1),
        ];
        let mut batch = MatrixBatch::default();
        let mut allocated_capacity = None;
        for col in 0..ENTRIES {
            batch
                .push(0, col, values[col % values.len()], &output)
                .unwrap();
            assert!(batch.text.len() < MATRIX_BUFFER_CAPACITY);
            let initial_capacity = *allocated_capacity.get_or_insert(batch.text.capacity());
            assert!(initial_capacity >= MATRIX_BUFFER_CAPACITY + ENTRY_CAPACITY);
            assert_eq!(batch.text.capacity(), initial_capacity);
        }
        assert!(
            output.lock().unwrap().nnz > 0,
            "the cell must flush before its last entry"
        );
        batch.flush(&output).unwrap();
        assert_eq!(batch.nnz, 0);
        assert!(batch.text.is_empty());
        let mut writer = output.lock().unwrap();
        assert_eq!(writer.finish().unwrap(), ENTRIES as u64);
        let (_, entries) = parse_matrix(writer.writer.get_ref().get_ref());
        for (col, &(row, actual_col, value)) in entries.iter().enumerate() {
            assert_eq!(
                (row, actual_col, value),
                (1, col + 1, values[col % values.len()])
            );
        }
    }

    #[test]
    fn coalesces_small_batches_and_header_padding() {
        let output = Mutex::new(MatrixMarketWriter::new(TestStream::default(), 100, 10).unwrap());
        let mut batch = MatrixBatch::default();
        for row in 0..100 {
            for col in 0..10 {
                batch.push(row, col, 1.25, &output).unwrap();
            }
            batch.flush(&output).unwrap();
        }
        let mut writer = output.lock().unwrap();
        assert!(writer.writer.get_ref().write_lengths.is_empty());
        assert_eq!(writer.finish().unwrap(), 1000);
        let stream = writer.writer.get_ref();
        assert_eq!(
            stream.write_lengths.len(),
            2,
            "one body write plus one patch write"
        );
        assert_eq!(stream.write_lengths[1], NNZ_FIELD_WIDTH);
        assert_eq!(parse_matrix(stream.data.get_ref()).0, (100, 10, 1000));
    }

    #[test]
    fn rejects_bounds_before_writing_and_keeps_failed_batch() {
        for (rows, cols, row, col) in [(2, 1, 2, 0), (2, 1, 0, 1), (0, 1, 0, 0)] {
            let output =
                Mutex::new(MatrixMarketWriter::new(TestStream::default(), rows, cols).unwrap());
            let mut batch = MatrixBatch::default();
            batch.push(row, col, 1.0, &output).unwrap();
            assert_eq!(
                batch.flush(&output).unwrap_err().kind(),
                io::ErrorKind::InvalidInput
            );
            assert_eq!(batch.nnz, 1);
            assert!(!batch.text.is_empty());
            let writer = output.lock().unwrap();
            assert_eq!(writer.nnz, 0);
            assert!(writer.writer.get_ref().write_lengths.is_empty());
        }
    }

    #[test]
    fn rejects_coordinate_and_count_overflow() {
        let output = Mutex::new(MatrixMarketWriter::new(TestStream::default(), 1, 1).unwrap());
        let mut batch = MatrixBatch::default();
        assert!(batch.push(usize::MAX, 0, 1.0, &output).is_err());
        assert!(batch.push(0, usize::MAX, 1.0, &output).is_err());
        assert!(batch.text.is_empty());
        batch.nnz = u64::MAX;
        assert!(batch.push(0, 0, 1.0, &output).is_err());
        assert!(batch.text.is_empty());

        let mut batch = MatrixBatch::default();
        batch.push(0, 0, 1.0, &output).unwrap();
        output.lock().unwrap().nnz = u64::MAX;
        assert!(batch.flush(&output).is_err());
        assert!(
            output
                .lock()
                .unwrap()
                .writer
                .get_ref()
                .write_lengths
                .is_empty()
        );
    }

    #[test]
    fn propagates_finalization_write_seek_and_flush_failures() {
        for fault in [Fault::Write, Fault::ShortWrite, Fault::Seek, Fault::Flush] {
            let mut writer = MatrixMarketWriter::new(TestStream::default(), 0, 0).unwrap();
            writer.writer.get_mut().fault = Some(fault);
            assert!(writer.finish().is_err());
            assert!(
                writer.finish().is_err(),
                "a partially completed stream must not be retried"
            );
            assert!(writer.append(&MatrixBatch::default()).is_err());
        }
    }

    #[test]
    fn preserves_batch_on_automatic_flush_failure() {
        let output = Mutex::new(MatrixMarketWriter::new(TestStream::default(), 1, 50_000).unwrap());
        output.lock().unwrap().writer.get_mut().fault = Some(Fault::ShortWrite);
        let mut batch = MatrixBatch::default();
        let mut failed = false;
        for col in 0..50_000 {
            if batch.push(0, col, 1.25, &output).is_err() {
                failed = true;
                break;
            }
        }
        assert!(failed);
        assert!(batch.nnz > 0);
        assert!(batch.text.len() >= MATRIX_BUFFER_CAPACITY);
        let mut writer = output.lock().unwrap();
        assert_eq!(writer.nnz, 0);
        assert!(writer.finish().is_err());
    }

    #[test]
    fn reports_poisoned_output_as_an_error() {
        let output = std::sync::Arc::new(Mutex::new(
            MatrixMarketWriter::new(Cursor::new(Vec::new()), 1, 1).unwrap(),
        ));
        let other = output.clone();
        assert!(
            std::thread::spawn(move || {
                let _guard = other.lock().unwrap();
                panic!("inject lock poisoning");
            })
            .join()
            .is_err()
        );
        let mut batch = MatrixBatch::default();
        batch.push(0, 0, 1.0, &output).unwrap();
        assert_eq!(
            batch.flush(&output).unwrap_err().kind(),
            io::ErrorKind::Other
        );
        assert_eq!(batch.nnz, 1);
    }
}
