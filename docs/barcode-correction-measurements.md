# Deterministic barcode-correction measurements

Measurements were collected on 2026-08-15 against alevin-fry `7fee462`
(the repaired PR #171 baseline) and libradicl PR #56. Times are wall-clock
seconds and RSS values come from GNU `time -v`. The benchmark inputs were on a
warm local filesystem; the three-run comparisons used the same commands,
thread counts, and output settings.

## Default Flex performance gate

The full Flex v2 input was
`/scratch1/rob/flex_bench/inv/flex-t32/af_map/map.rad`: 39,902,036,069 bytes
and 1,994,961,121 RAD records. GPL used 32 threads, Exact sample correction,
Hamming-1 Unique cell correction, the local 384-sample list, the 737,280-cell
whitelist, and a minimum exact count of 10. Collation used 32 threads and a
2 GiB memory budget.

| Stage | Implementation | Wall times (s) | Median (s) | Peak RSS runs (KiB) | Median RSS (KiB) | Median change |
|---|---|---:|---:|---:|---:|---:|
| GPL | immediate baseline | 32.64, 32.21, 32.15 | 32.21 | 2,578,804; 2,640,924; 2,485,752 | 2,578,804 | — |
| GPL | compiled correction | 31.77, 30.94, 31.32 | 31.32 | 1,753,088; 1,792,952; 1,780,888 | 1,780,888 | wall −2.8%; RSS −30.9% |
| Collate | immediate baseline | 24.14, 24.50, 24.66 | 24.50 | 1,302,004; 1,291,912; 1,313,240 | 1,302,004 | — |
| Collate | fused compiled lookup | 18.57, 21.92, 22.82 | 21.92 | 1,244,672; 1,231,896; 1,255,520 | 1,244,672 | wall −10.5%; RSS −4.4% |

The default correction plan was 150,659,022 bytes. Cell-correction compilation
took about 1.94–1.99 seconds and the plan write took about 66 ms. Three GPL
runs produced byte-identical plans despite independent worker scheduling.

The baseline and feature collated RADs had the same exact byte length
(37,402,724,053), identical collation manifests, and identical sample names,
barcodes, read totals, and cell totals. Record order inside a cell is not part
of the RAD contract and differed between parallel runs, so whole-file hashes
are not expected to match.

### 0.18.0 registry-dependency release gate

Immediately before release, the full Flex GPL and collation measurements were
repeated to compare current `master` (using its pinned libradicl Git revision)
with the 0.18.0 release candidate (using `libradicl = "0.18.0"` from
crates.io). Both builds used their locked release dependency graph.

| Stage | Build | Wall times (s) | Median (s) | Peak RSS runs (KiB) | Median RSS (KiB) | RC change |
|---|---|---:|---:|---:|---:|---:|
| GPL | current master | 10.24, 10.27, 10.23 | 10.24 | 1,575,844; 1,598,820; 1,584,600 | 1,584,600 | — |
| GPL | 0.18.0 RC | 9.95, 10.16, 10.37 | 10.16 | 1,590,256; 1,582,072; 1,619,352 | 1,590,256 | wall −0.8%; RSS +0.4% |
| Collate | current master | 20.64, 20.99, 20.85 | 20.85 | 1,224,836; 1,255,840; 1,230,088 | 1,230,088 | — |
| Collate | 0.18.0 RC | 21.20, 20.04, 21.14 | 21.14 | 1,237,232; 1,238,248; 1,252,016 | 1,238,248 | wall +1.4%; RSS +0.7% |

All six GPL runs produced the same correction-plan hash. All six collations
produced the same manifest hash, 2,423,889 chunks, 1,869,082,132 written
records, and a 37,402,724,053-byte output. The registry transition therefore
cleared the 3% wall-time and peak-RSS release gates without changing the
output contract.

The 2.6 GB unfiltered single-barcode PBMC GPL was also measured three times at
8 threads. Baseline wall times were 2.25, 2.27, and 2.26 seconds (median 2.26)
with median RSS 357,124 KiB. Feature times were 2.24, 2.25, and 2.35 seconds
(median 2.25) with median RSS 328,104 KiB: wall −0.4% and RSS −8.1%. Both
implementations reported exactly 73,151 exact targets, 76,734 corrected
distinct observations / 983,308 corrected reads, 2,942 ambiguous observations /
62,254 ambiguous reads, and 484,771 observations / 1,979,827 reads without a
candidate. The three feature plans were byte-identical.

Collation of the same single-barcode input used 8 threads and a 2 GiB memory
budget. Baseline wall times were 1.77, 1.74, and 1.78 seconds (median 1.77)
with median RSS 906,500 KiB. Compiled-plan times were 1.64, 1.66, and 1.65
seconds (median 1.65) with median RSS 896,504 KiB: wall −6.8% and RSS
−1.1%. During this measurement, retaining a second decoded correction map
was found to amplify allocator high-water use; the final handoff instead
borrows the one map for unmapped-read accounting and then moves it into
libradicl's fused lookup.

## Final hash and count aggregation optimization

A subsequent release-candidate pass replaced private multi-sample routing and
cell-count maps with AHash-based execution structures, made the immutable cell
whitelist an `AHashSet`, and aggregated counts in bounded worker-local maps.
The correction engine and persisted artifacts remained ordered and
deterministic; the hash implementation is private and does not choose targets.

The full Flex v2 input above was measured with 32 threads and three runs per
main comparison:

| Implementation | Median wall (s) | Median peak RSS (MiB) |
|---|---:|---:|
| Original resolved-correction GPL | 31.39 | 1,735.5 |
| Fused route, SipHash | 28.80 | 1,697.9 |
| Fused route, AHash | 19.43 | 1,741.5 |
| Immutable AHash whitelist | 17.02 | 1,703.0 |
| Bounded worker-local counts (final) | **9.35** | **1,550.2** |

The final implementation reduced median wall time by 70.2% and peak RSS by
10.7% relative to the original resolved-correction implementation.  The
forced-cell path improved from 15.65 seconds / 1,402.3 MiB to 9.84 seconds /
1,352.4 MiB.  Cell-Frequency and sample-Frequency completed in 8.03 and 8.14
seconds and reproduced their pre-change correction artifacts byte-for-byte.

A six-run FoldHash comparison at the fused-routing stage measured 19.22 seconds
versus 19.43 seconds for AHash.  The roughly 1% and inconsistent difference did
not justify another direct dependency or a different fixed-seed tradeoff.

On a 28.4-million-record, 3.8 GB ATAC RAD, 13 paired warm-cache sort runs
improved from 0.76 to 0.72 seconds median (5.3%).  The BED outputs were
byte-identical.

## Neighbourhood effects

The following tables compare Hamming-1 with the historical
substitution-or-shift-1 neighbourhood. “Extra rejection” means that Hamming
accepted an observation but an added shift candidate made Unique ambiguous or
made Frequency fall below confidence. “Changed target” means both policies
accepted but chose different canonical targets. Each entry is
`distinct observations / represented reads`.

### Filtered 5k PBMC

This used the existing `all_freq.bin` from the 11 GB filtered 5k PBMC run:
2,770,385 distinct observations, 309,547,179 reads, and 4,784 retained targets.

| Resolution | Identical target | Shift-only rescue | Extra rejection | Changed target |
|---|---:|---:|---:|---:|
| Unique | 230,854 / 286,528,454 | 86,295 / 1,649,338 | 236 / 4,427 | 0 / 0 |
| Frequency, 97.5% | 230,854 / 286,528,454 | 86,295 / 1,649,338 | 236 / 4,427 | 0 / 0 |

For both policies, 4,784 target frequencies changed. The total absolute target
read delta was 1,645,101 and the largest single-target delta was 8,055 reads.

### Unfiltered 1k PBMC

This used the 2.6 GB unfiltered PBMC RAD and its whitelist: 637,598 distinct
orientation-compatible observations, 54,030,084 reads, and 73,151 retained
targets.

| Resolution | Identical target | Shift-only rescue | Extra rejection | Changed target |
|---|---:|---:|---:|---:|
| Unique | 146,775 / 51,930,928 | 32,877 / 295,696 | 3,110 / 57,075 | 0 / 0 |
| Frequency, 97.5% | 150,782 / 52,038,609 | 33,812 / 306,331 | 632 / 1,097 | 610 / 8,653 |

Unique changed 15,036 target frequencies (268,797 total absolute read delta;
maximum 1,975). Frequency changed 15,034 target frequencies (322,740 total
absolute read delta; maximum 2,210).

### Full Flex v2

The sample correction policy was held at Exact, so the comparison covers cell
correction within the same canonical sample. Reads rejected at the sample
level are reported separately and are not included in the cell outcomes.

| Resolution | Identical target | Shift-only rescue | Extra rejection | Changed target |
|---|---:|---:|---:|---:|
| Unique | 9,001,481 / 1,867,649,181 | 4,224,798 / 28,846,196 | 413,698 / 1,432,951 | 0 / 0 |
| Frequency, 97.5% | 10,021,860 / 1,872,552,240 | 4,512,244 / 31,011,577 | 100,878 / 242,241 | 55,188 / 196,478 |

Unique changed 974,071 corrected target frequencies (28,014,137 total
absolute read delta; maximum 4,147). Frequency changed 975,401 target
frequencies (31,208,398 total absolute read delta; maximum 4,200).

The Hamming Unique and Frequency correction phases took 1.99 and 2.06 seconds;
the shift-aware phases took 4.54 and 4.78 seconds. Their correction-plan sizes
were respectively 150.7 MB, 162.9 MB, 211.6 MB, and 233.5 MB. Whole GPL wall
time was 31.50 seconds for Hamming Unique and 34.07 seconds for shift-aware
Unique.

These results support retaining shift-aware correction as explicit
compatibility behaviour, but not changing the Hamming default for unfiltered,
multi-sample, or ATAC paths in this campaign.

## Sample Frequency and spill/replay

On the full Flex input with 256 MiB of deferred-buffer budget, sample Frequency
routed:

- 1,937,727,066 exact reads;
- 29,783,475 structurally unique non-exact reads;
- 27,450,580 reads with no structural target; and
- zero structurally ambiguous reads.

Because the local Flex sample list has no ambiguous Hamming-1 observations, the
correct production behaviour was zero spill runs and zero spool bytes. The run
finished in 32.74 seconds with 1,915,412 KiB peak RSS and produced a
152,236,015-byte plan. Its plan was byte-identical before and after the hot-path
optimization that replaced shared per-read prior updates and repeated
neighbour discovery with worker-local prior counters and a compiled direct
ambiguity index. That optimization reduced this run from 131.93 seconds to
32.74 seconds.

The synthetic end-to-end spill test deliberately creates one ambiguous and
one impossible sample barcode. It observes 101 exact reads, 10 deferred reads,
and 5 immediately rejected reads; forces a spill run; resolves all 10 deferred
reads to the high-prior sample at 97.5%; replays them into the cell histogram;
collates 110 records for that sample plus the other sample's one exact record;
and verifies that the RAII spool file is gone. Unit tests separately force
zero, one, and multiple runs, duplicate run-length aggregation, replay
equality, corruption reporting, and cleanup.

## Reproducible reports

The developer examples `correction_policy_report`,
`multi_correction_plan_report`, and `rad_barcode_frequencies` generate the
policy reports above. `compare_correction_artifacts` checks that a compact
single-barcode plan has exactly the same decisions as its legacy map. They
live under `examples/` and are deliberately not installed as supported
alevin-fry commands.
