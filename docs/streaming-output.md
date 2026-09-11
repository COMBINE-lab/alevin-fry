# Streaming quantification output

Quantification writes the count matrix and optional bootstrap mean/variance
matrices without retaining triplets for the whole experiment. Each worker
formats entries into a reusable 256 KiB batch outside the matrix lock. A batch
is appended under the lock when it reaches the threshold, including midway
through a large cell; each worker flushes its final partial batch before exiting.
The file writer also has a 256 KiB buffer.

Matrix serialization memory therefore scales with workers and active output
matrices rather than total nonzeros. Each active worker batch reserves 256 KiB
plus room for one formatted entry, and each active matrix has one 256 KiB file
buffer. Bootstrap inference can still retain the current cell's replicate
vectors; the change removes accumulation across cells.

Rows are assigned while appending paired barcode and feature records under
their existing lock. Both text records are formatted before acquiring that
lock. Matrix coordinates retain the assigned row index, so batches from
different workers may interleave safely. Barcode/feature files use explicit
256 KiB buffers too.

The MatrixMarket header reserves 20 characters for its final nonzero count.
Finalization seeks through `BufWriter`, which flushes the body, writes the
padded count through the buffer, and explicitly flushes the patch. Trailing
whitespace on the size line is intentional. Batches are checked against the
declared row and column dimensions before writing. A barcode filter that
selects more sample/cell pairs than its declared row count returns an error;
support for counting such filtered pairs is a separate concern.

Output I/O failures propagate through the worker result. Other workers stop
quantifying and drain queued input, allowing the producer to finish without
blocking on its bounded queue. All workers are joined before returning an
error. Reader errors and final buffer flush errors also propagate, and failure
paths do not write new success metadata. The CLI preserves the error context
and exits unsuccessfully instead of replacing these errors with a panic.
This lifecycle handles returned I/O errors; it is not general recovery from
arbitrary panics in quantification or record parsing.

Bootstrap mean and variance calculations preserve their existing semantics.
`--summary-stat` retains the population-variance convention used by bootstrap
inference; summarizing full replicate vectors retains the sample-variance
convention. No bootstrap files are produced when bootstrapping is disabled,
and files are removed on successful completion when no bootstrap means were
produced, matching the previous output file set.

Unit tests cover bounded batches spanning a large cell, numeric roundtrips,
coordinate/count overflow, bounds validation, write coalescing, poisoned locks,
and injected write/seek/flush failures. Integration tests compare every
barcode/gene count with one and multiple workers, cover empty and fractional
output and USA axes, exercise deterministic bootstrap statistics in both modes,
and use subprocess deadlines to test output failures and truncated input.

For performance measurements, compare release builds with the same dependency
lockfile, target CPU flags, thread count, input, and quantification options.
Use separate Cargo target directories for different source checkouts. Measure
wall time and peak RSS without tracing, and collect filesystem-write counts in
separate traced runs. Align complete barcode and gene axes when comparing
matrices; aggregate row/column sums are insufficient to establish equivalence.
