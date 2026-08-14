# Benchmarks

Three [divan](https://docs.rs/divan) benchmark targets covering the routines a
dependency or algorithm change is most likely to move:

| target       | covers                                                                                  |
|--------------|-----------------------------------------------------------------------------------------|
| `em`         | `em_optimize_subset` (the per-cell EM driver behind the `*-em` resolutions), and loading the dumped equivalence classes back into the flat CSR layout |
| `barcode`    | `generate_permitlist_map`, `generate_whitelist_set`, `get_all_one_edit_neighbors`         |
| `collate_io` | Snappy encode/decode of a temporary-bucket-sized buffer, and bincode round-trip of the barcode-correction map |

```sh
cargo bench                      # all three
cargo bench --bench em           # one
cargo bench --bench em -- batch  # filter by name
```

## Running on real data

The benchmarks default to a deterministic synthetic generator so that
`cargo bench` works in a fresh checkout, and they print a note on stderr when
they do. Synthetic numbers are fine for comparing two branches against each
other; they are not representative of real input, and should not be quoted as
if they were.

Point `AF_BENCH_DATA` at a directory holding real artifacts to switch every
bench over:

```text
$AF_BENCH_DATA/
  gene_eqclass.txt.gz   # quant -r cr-like-em --dump-eqclasses
  geqc_counts.mtx       # ditto
  permit_list.txt       # one barcode per line, e.g. 10x_v3_permit.txt
  permit_map.bin        # generate-permit-list output
  map.collated.rad      # collate
```

Nothing is committed to the repository: the smallest useful fixture set is
around 1.2 GB. The set above is what the project's own CI toy dataset produces:

```sh
wget https://umd.box.com/shared/static/uctm7oh0f2ef32aa3dmvk6zvbb8764up -O toy_data.tar.gz
tar xzf toy_data.tar.gz
fry=target/release/alevin-fry
$fry generate-permit-list -u toy_data/10x_v3_permit.txt -d fw -i toy_data/alevin_map -o gpl
$fry collate -i gpl -r toy_data/alevin_map -t 8
$fry quant -r cr-like-em --small-thresh 0 --dump-eqclasses --use-mtx \
     -m toy_data/t2g_3col.tsv -i gpl -o quant -t 8

mkdir fixtures && cd fixtures
ln -s ../quant/alevin/gene_eqclass.txt.gz ../quant/alevin/geqc_counts.mtx .
ln -s ../gpl/map.collated.rad .
ln -s ../gpl/permit_map.bin .
ln -s ../toy_data/10x_v3_permit.txt permit_list.txt
```

## Comparing two branches

divan reports a distribution per branch but does not pair samples across
builds, so a difference is only believable when it is larger than the spread of
the baseline. The protocol used for the performance changes in this repository:

```sh
export AF_BENCH_DATA=/path/to/fixtures
git checkout master && cargo bench --bench em 2>&1 | tee "${TMPDIR:?}/before.txt"
git checkout my-branch && cargo bench --bench em 2>&1 | tee "${TMPDIR:?}/after.txt"
```

Compare the **median** columns, not the fastest, and re-run both sides on an
otherwise idle machine. Deltas under roughly 3% on these targets are inside the
run-to-run noise observed on an M-series laptop and should not be reported as
wins.

A microbenchmark delta is also not a pipeline delta: `collate` is I/O bound and
`quant` is threaded, so any claim about end-to-end time has to come from timing
the actual subcommands on a real dataset.
