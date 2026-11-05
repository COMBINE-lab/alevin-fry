# Optimized alevin-fry (long-read, per-cell truncation prior): what changed, step by step

Source: `/scratch1/zahra/source-code/alevin-fry` and `/scratch1/zahra/source-code/libradicl` (the
path dependency), copied on 2026-10-30 from `/nfshomes/zzare/ScLongRead/develeop-refactor/alevin-fry`
(commit 758965d, "all truncation models") and from the `libradicl` working tree next to it. Nothing
is committed: every change is in the working tree (`git diff` in both repositories shows it all),
and this directory holds the exact scripts that made each step, so the steps can be replayed or
reviewed one at a time.

Build (fat LTO, about 75 s):

```
export CMAKE=/nfshomes/zzare/miniconda3/bin/cmake CARGO_TARGET_DIR=/scratch5/zahra/tmp/optim/target
cd /scratch1/zahra/source-code/alevin-fry && cargo build --release
```

The configuration everything was measured and verified on is the recommended one of the per-cell
truncation report, with the likelihood deduplication engine and the cell-wide calibrated edge model
(not the transcript-length one):

```
alevin-fry quant -i <collated dir> -m <transcript map> -t 8 -o <out> -r parsimony-em \
  --lambda-size=0.0 --tau-delta=0.0 \
  --ed-distance-model general --ed-upper-bound 2 --dedup-edges bayes --ed-far-window 10 \
  --model-truncation --truncation-interval-mode all --truncation-end-kernels absolute \
  --truncation-cell-prior 0.25 --truncation-positional-power 0.25
```

`--dedup-model` defaults to `bf`, the likelihood engine. `--true-umi-mode off` is still accepted (and
ignored) so that existing command lines keep working.

## Verification on every data set of the per-cell report (`verify_all.sh`, 2026-11-01)

Original scratch binary against the optimized one (step 10), 8 threads each, in both deduplication
configurations: the requested one (`--dedup-edges bayes --ed-far-window 10`) and the report's PC13
configuration (rule edges, general ed2). All twelve comparisons are **bit-identical** (same cells,
same nonzero entries, maximum absolute difference 0; `/scratch5/zahra/tmp/optim/verify/compare.log`).
Wall times below were measured with all six data sets running at once on the same machine, so they
are comparable within a row, not across rows.

| data set | config | original | optimized | peak RSS original → optimized |
|---|---|---|---|---|
| new ONT sim (1,467 cells) | bayes | 10:12 | 5:23 | 3.05 → 2.81 GB |
| new ONT sim | rule | 10:39 | 4:36 | 3.06 → 2.80 GB |
| new PacBio sim | bayes | 3:59 | 2:36 | 1.37 → 1.34 GB |
| new PacBio sim | rule | 3:54 | 2:37 | 1.37 → 1.34 GB |
| old ONT sim | bayes | 7:20 | 3:28 | 2.60 → 2.33 GB |
| old ONT sim | rule | 8:37 | 2:36 | 2.51 → 2.40 GB |
| old PacBio sim | bayes | 2:42 | 1:48 | 1.09 → 1.10 GB |
| old PacBio sim | rule | 2:30 | 1:38 | 1.10 → 1.08 GB |
| cell lines ONT (8 lines) | bayes | 4:26 | 1:49 | 1.07 → 1.01 GB |
| cell lines ONT | rule | 4:16 | 1:14 | 0.97 → 1.09 GB |
| cell lines PacBio (8 lines) | bayes | 11:45 | 4:28 | 1.71 → 1.49 GB |
| cell lines PacBio | rule | 12:18 | 2:33 | 1.57 → 1.54 GB |

## Verification during development

Each step was checked on the 96-cell TKSM pilot (ONT, UMI 12, threshold 0.95, 5.6 M reads) and the
final binary also on the full new-ONT simulation (1,467 cells, 21 GB collated RAD): the
quantification matrix is **bit-identical** to the original scratch binary's (`cmp_quants.py`: same
cells, same nonzeros, maximum absolute difference 0). The truncation-model and UMI-distance test
suites pass (`cargo test --release`), and the bounded distance is cross-checked against the full one
on 100,000 random pairs for UMI lengths 6 to 16 and bounds 0 to 5, plus every pair of 5-mers.

## Timings (`bench.sh`: pilot, 8 threads, best of two runs, each run compared with the baseline)

| state | wall | CPU | peak RSS |
|---|---|---|---|
| original scratch binary | 30.8 s | 211 s | 2.39 GB |
| after steps 1–2 (diagnostics, read names) | 30.2 s | 207 s | 2.19 GB |
| after steps 3–5 (bounded distance, reused tables, no spinning) | 21.2 s | 141 s | 2.32 GB |
| after step 6 (fat LTO) | 20.9 s | 139 s | 2.29 GB |
| after step 7 (write-back over own candidates, composition bound) | 13.5 s | 88 s | 2.21 GB |
| after step 8 (ahash lookup maps; no measurable change) | 13.9 s | 90 s | 2.26 GB |
| final (step 5 amended: sleep instead of yield) | 13.5 s | 87 s | 2.21 GB |
| after step 9 (queue holds only the chunk's bytes) | 13.7 s | 88 s | 2.22 GB |

Full simulation (new ONT, 1,467 cells, 21 GB collated RAD), 8 threads: original 8:40 wall,
4,068 s user + 35 s system CPU, 3.17 GB peak; final binary 4:31 wall, 1,845 s user + 18 s system,
2.97 GB peak (the step-8 binary with the yield loop: 4:32, 1,886 s + 252 s system, 2.93 GB); the
matrix over 1,467 cells and 193,169 transcripts is identical in both runs.

Where the time went before (perf, pilot): 35% in the UMI edit distance, 21% in a `Vec::from_iter`
that was the fit's write-back (step 7), 8% in the RAD reader, 8% in the truncation EM. After: the
bounded distance and its filters 23%, the truncation EM's per-alignment terms about 22%, the rest
spread over the deduplication, the equivalence-class map and the reader.

## The steps

### Step 0: the code that produced the reported results, into the copy
The copy starts from your repository and takes the sources of the scratch build that produced the
reported results (`/scratch5/zahra/tmp/percell/alevin-fry`): the per-cell prior
(`--truncation-cell-prior`, `--truncation-cell-iterations`), the likelihood deduplication engine
(`--dedup-model bf`), the calibrated edge model (`--dedup-edges bayes` and `bayes-tx`,
`--ed-far-bound`, `--ed-far-window`), the expected-count and posterior-cover options, and the
UMI-length fix (the UMI length comes from the RAD file's `ulen` tag instead of a hard-coded 12). Your
repository's sources are a strict subset of these (the diff is additions only), so nothing of yours
was dropped. Verified identical to the scratch binary.

### Step 1: diagnostics and instrumentation removed (`step1_remove_diagnostics.py`)
Removed: the molecule dump (`MOLDUMP_*` environment variables), the UMI-graph dump (`DUMP_CELL_NUM`),
the coverage dump (`COVERAGE_DIAG_*`), the ground-truth UMI oracle (`--true-umi-mode edge|bfs|direct`,
`--true-txp-mode`, `--true-umi-txp-file`, the `TRUE_UMI_ORACLE` totals and `[ORACLE_CELL]` prints,
about 1,300 lines of `pugutils_dp_new.rs`), the cell-647 debug breakpoint and the EM debug prints it
switched on, the per-cell map from (candidate set, UMI) to read names and the read names kept per
equivalence class. The last two ran on every cell of every run: a HashMap keyed by cloned
`Vec<u32>` candidate lists holding `String` read names, built from all reads of the cell, then one
`Vec` clone per equivalence class per UMI to fill a field nothing read. Also removed the imports they
left unused.

### Step 2: read names out of the RAD records (`step2_read_names.py`, in `libradicl`)
`libradicl` was reverted to its committed state (your uncommitted change there was exactly the
read-name support; it is kept in `libradicl_qname_user.patch`) and the long-read record now treats the
read-name field as one length byte that is always 0: `convert` writes a 0, the record writer writes a
0, the parsers skip the announced number of bytes on a stack buffer without allocating, and the
per-record size counts one byte. Files written earlier (with names) still read, since the parser
skips the name. Files written now are 35% smaller (pilot: `map.rad` 715 MB instead of 1.10 GB,
`map.collated.rad` 699 MB instead of 1.08 GB), and `convert`, `collate` and the chunk reader no longer
build a `String` per read. This also fixes a latent bug: the old size estimate assumed every read
name was 68 characters (true for BLAZE-style names with 16-base barcodes and 12-base UMIs, false for
other UMI lengths or other tools), and `collate` positions records by that estimate. The short-read
record type goes back to upstream (no name byte), so standard piscem RAD files parse again. The
`view` command prints `-` in the read-name column.

### Step 3: bounded UMI edit distance (`step3_edit_distance.py`, `umi_general_distance.rs`)
The deduplication only needs to know whether two UMIs are within distance 3 (calibrated edge model)
or within the far bound (rule model); every larger distance means "different". The distance function
now takes that bound and returns `None` beyond it, with exact filters in front of the dynamic
program: a composition lower bound (the bases that cannot be matched, `len - sum over bases of
min(count_x, count_y)`; added in step 7), a pigeonhole test (one of `bound + 1` pieces of either UMI
must appear verbatim in the other within `bound` positions; the pieces of *either* string, because
the free end-deletions of the shift-tolerant distance can spoil pieces of the string being cut at no
cost, which the first version missed and the random cross-check caught), and a banded dynamic
program (cells with |i − j| > bound cannot lie on a path of cost at most bound, as a shift of s costs
at least s insertions) that stops as soon as a whole row exceeds the bound. Every filter is a
necessary condition, so the decisions are unchanged. On random 12-mers 27% of pairs reach the
dynamic program; a 2-gram count bound would bring that to 19%, not enough to be worth its tables.

### Step 4: the truncation model's lookup tables reused (`step4_exps_reuse.py`)
Every EM iteration of a cell's fit rebuilt from scratch two tables with one entry per base up to the
longest transcript (the end-kernel density and its running sum) and the length-law tables:
megabytes allocated, zeroed and filled per iteration, about 50 iterations per cell. They are now
refreshed in place into buffers that live for the whole fit, with the same arithmetic in the same
order (bit-identical values); the same change in the two global fits.

### Step 5: the chunk reader waits instead of spinning (`step5_reader_wait.py`, `libradicl/src/readers.rs`)
When the work queue was full the reader thread spun in an empty loop at 100% CPU, counted as user
time. It first became `yield_now` (which moved the waste into 252 s of system time on the full
simulation), then a 100 µs sleep: a chunk takes milliseconds to process, so the workers never wait
for the reader. Wall time is unchanged, CPU time stops being burnt.

### Step 6: build profile (`Cargo.toml`)
`lto = "fat"`, `codegen-units = 1`: cross-crate inlining of the libradicl parsers and the hashing and
vector code in the hot loops. Worth about 2%; the release build takes 75 s instead of 35 s.

### Step 7: posteriors over the molecule's own candidates (`step7_writeback_and_bound.py`)
The largest single win. At the end of a cell's fit, every molecule allocated and zeroed one f64 per
transcript of the cell (tens of thousands) and then scanned all of them to keep the nonzero shares:
with 30,000 molecules per cell that is gigabytes touched per cell, which perf showed as 24% in
`Vec::from_iter` plus 4% in `memset`. The shares are now accumulated over the molecule's own
alignments, in ascending transcript order with repeated transcripts summed in alignment order,
which gives exactly the same values in the same order; the per-molecule buffers are reused. The
composition bound of step 3 was added in the same step. The pooled assignment (`assign_pooled`,
`--truncation-abundance global|shrink`) still uses the dense vector; it is not on the per-cell path.

### Step 8: fixed-seed ahash for the per-cell lookup maps (`step8_hashers.py`)
The interval map per (equivalence class, UMI), the vertex position ranges and the pair posteriors are
only ever looked up by key, or iterated to build another lookup map, so the hasher cannot change a
result; std's SipHash was replaced by ahash with fixed seeds. No measurable change on the pilot; kept
because it cannot hurt.

### Step 9: the work queue holds only each meta-chunk's bytes (`step9_chunk_blob.py`, `libradicl/src/readers.rs`)
The reader's buffer starts at 512 KB and grows to the largest cell seen, never shrinks, and every
meta-chunk pushed to the queue (four slots per worker) was a clone of the whole buffer; the iterator
only parses the first `num_bytes`, so the clone is now that prefix. Bytes seen by the workers are
unchanged (pilot identical). No measurable change in peak memory or time on the pilot, so the queue
was not what set the peak; kept because it removes needless copying on data sets with big cells.

### Step 10: no second copy of the reference names (`step10_ref_names_copy.py`)
The reference names are needed (see Memory below), but the fork also cloned all of them into an
`Arc` for the worker threads, whose only consumer was the removed oracle. The group, gene and oracle
file loaders run on the main thread and now take a slice of the header's names; the worker field is
gone. About 24 MB less per run on a 507 k-transcript reference; pilot identical (13.6 s, 2.18 GB).

## Memory

The reference names themselves have to stay: the RAD header stores one name per reference, and
quantification needs them to translate the transcript-to-gene file (which is written in names) into
the file's transcript ids and to write the names into `quants_mat_cols.txt`; the original alevin-fry
does exactly this (`hdr.ref_names` → name-to-id map → `parse_tg_map`). They are read once (about
12 MB of strings for 507 k transcripts) and now kept once.

Peak memory is almost entirely per worker thread: the pilot peaks at 0.74 GB with 2 threads and
2.2 GB with 8, so `-t` sets the footprint. The steps above take the pilot from 2.39 GB to 2.21 GB and
the full simulation from 3.17 GB to 2.97 GB (removed read names and per-cell name maps, no dense
per-molecule vectors, reused kernel tables, trimmed queue clones).

Where it goes (massif, system allocator, one worker, the largest pilot cell with 67,000 molecules;
`massif_largest_cell_peak.txt`), peak 550 MB:

| what | MB | per cell or fixed |
|---|---|---|
| transcript-to-gene maps (`parse_tg_map`: two `HashMap<String, u32>` over 507 k transcripts and their genes, plus a rehash) | ~100 | fixed |
| reference names (`RadPrelude`) and their copy in `do_quantify` | ~36 | fixed |
| the work queue (four meta-chunks per worker) | ~54 | fixed per worker |
| the per-molecule interval maps (`intervals_for_resolved_vertices`: one `HashMap<u32, Vec<AlignmentInterval>>` per resolved molecule) | ~44 | per cell |
| the resolved molecules themselves (transcript probabilities, vertices) | ~44 | per cell |
| the truncation fit: the alignment array (26) and the end-kernel tables (one f64 per base up to the longest candidate transcript, four tables, ~33) | ~60 | per cell |
| the equivalence-class map with its per-vertex intervals | ~24 | per cell |
| the parsed reads of the cell (six vectors per read) | ~24 | per cell |

So a big cell costs a worker about 200 MB while it is being processed, and the peak of a run is
the base plus the sum over workers of the cells they happen to hold at the same moment. The
remaining per-cell items are the model's data structures, not waste; the next reductions would be a
flat interval array per molecule instead of a map (about 20 MB per big cell, and it would also make
the coverage model's summation order deterministic) and kernel tables indexed by bin instead of by
base (the running sum would then be computed differently, so not bit-identical). Neither was done.

## What was not changed, and two things to know

* Every model option is still there (pooled and shrunk abundances, groups and their discovery, the
  coverage, junction and four-parameter models, the end models, bootstraps, USA mode);
  `--truncation-dump` stays, it is an analysis output of the groups workflow, not instrumentation.
* **Results depend on read order.** Re-running `convert`, `generate-permit-list`, `collate` and
  `quant` with the *original* binary on the pilot BAM gives a matrix that differs from the earlier run
  of the same pipeline (sum of counts 2,520,333 against 2,520,426; 2.8 M entries differ by small
  amounts): parallel `convert` and `collate` do not preserve the order of reads, and the calibrated
  edge model samples every 20th far pair in enumeration order, which depends on it. The optimized
  binary behaves the same way, and its full-pipeline run differs from the baseline by the same amount
  (`pipeline_test.sh`). A canonical ordering of UMIs within equivalence classes would make
  quantification independent of read order; it changes which pairs are sampled, so it was not done.
* The next real speedups are structural: enumerating only the near UMI pairs through an index (needs
  an order-independent sample of far pairs for the edge model's calibration, which changes results
  slightly), and caching per-alignment invariants in the truncation E-step (its reciprocals change
  rounding unless done carefully).

## Tools in this directory

* `bench.sh <binary> <tag> [runs] [threads]`: times the pilot and compares with the baseline matrix in
  `/scratch5/zahra/tmp/optim/baseline_pilot`.
* `cmp_quants.py <dirA> <dirB>`: aligns two quantifications by barcode and reports the differences.
* `pipeline_test.sh <binary> <out dir>`: convert, permit list, collate and quant on the pilot BAM.
* `umi_general_distance_v2.rs`: the distance module as installed in `src/`.
* `libradicl_qname_user.patch`: your read-name change to libradicl, as it was before step 2.
