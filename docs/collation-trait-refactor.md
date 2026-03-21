# Multi-barcode collation: trait-based grouping key refactor

## The problem

The generic function `collate_temporary_bucket_twopass_generic` in libradicl groups records within each temp bucket by `tup.collate_key().into()` (line 338 of `libradicl/src/lib.rs`). For `MultiBarcodeReadRecordT`, `collate_key()` returns `self.barcodes.last()` — the cell barcode only.

In multiplexed (multi-barcode) experiments like 10x Flex, different samples can share the same cell barcode. During collation, cells from different samples are packed into temp buckets sequentially. When the bucket boundary falls mid-sample, a single bucket can contain cells from multiple samples. If two cells from different samples share the same cell barcode and land in the same bucket, the gather phase merges their records into a single output chunk — because they have the same `collate_key()`.

This caused two observable bugs:

1. **Fewer output chunks than expected:** The collated RAD file had 594,830 chunks instead of the expected 595,200 (a 4-plex Flex experiment). The ~370 missing chunks were cross-sample merges.

2. **Wrong sample assignment:** The quant step determined each cell's sample by mapping its chunk index to a sample via the collation manifest. Since the multi-threaded gather phase writes chunks in non-deterministic order, chunk positions didn't match the manifest's expected ordering. This caused ~2,000 cells to be assigned to the wrong sample, producing duplicate `sample_cell` barcode strings in the output.

## Current fix (in alevin-fry)

We fixed this with four commits on the alevin-fry `master` branch:

### 1. Multi-barcode-aware gather (`collate_multi_bc_bucket`)

Added a local function `collate_multi_bc_bucket` in `src/collate.rs` that replaces `collate_temporary_bucket_twopass_generic` for the multi-barcode gather phase. It computes a composite grouping key `(sample_idx << cell_bc_bits) | (cell_bc & cell_bc_mask)` from all barcodes in the record header, ensuring cells from different samples always produce separate chunks.

The shift is computed dynamically from the actual cell barcode length (`cell_bc_len * 2` bits per nucleotide) rather than being hardcoded, so it works for any chemistry (including Flex v2 with 10bp sample barcodes).

This composite key is computed consistently across four locations:
- `tsv_map` construction (bucket assignment planning)
- Scatter phase (bucket routing)
- Gather phase (chunk grouping within each bucket)
- Manifest construction (extracting sample boundaries from composite keys)

### 2. Sample identity from records, not chunk position

The quant step now reads the sample index directly from `barcodes[0]` of each chunk's first record, instead of using a positional `cell_sample_idx[chunk_number]` lookup via the manifest. This makes sample assignment authoritative per-record, independent of chunk ordering.

The sample index extractor is passed as an `Option<Arc<dyn Fn(&R) -> usize + Send + Sync>>` closure into `do_quantify`. For `MultiBarcodeReadRecord` it reads `rec.barcodes[0] as usize`; for single-barcode types it's `None`.

### 3. Sample index (not barcode) in collated records

The scatter phase writes the integer sample index into `barcodes[0]` instead of the canonical sample barcode sequence. This means quant reads `barcodes[0]` as a direct index into a `Vec<String>` of sample names — no HashMap needed. Since sample_idx, canonical barcode, and sample name are all 1-to-1 after correction, the compact integer index is sufficient and anything else can be recovered from `sample_info.json`.

### 4. Overflow check

Added a check that `ceil(log2(num_samples)) + cell_bc_bits <= 64` before collation begins. Fails with a clear error if the composite key would overflow u64.

## The plan: upstream trait-based collation key to libradicl

The current fix works but is architecturally unsatisfying: `collate_multi_bc_bucket` in alevin-fry is essentially a copy of `collate_temporary_bucket_twopass_generic` from libradicl with a different grouping key. The proper fix is to make the grouping key trait-based so the generic function handles both record types.

### Changes to libradicl

**1. Add `collation_group_key()` to `CollatableMappedRecord`:**

```rust
// In record.rs, inside the CollatableMappedRecord trait:

/// Returns the key used to group records into output chunks during collation.
/// For single-barcode records, this is the cell barcode (same as collate_key()).
/// For multi-barcode records, this is a composite key encoding all barcode
/// levels (e.g., sample + cell) to prevent cross-sample merging.
fn collation_group_key(&self) -> u64 {
    // Default: same as collate_key() — correct for single-barcode types
    self.collate_key().into()
}
```

**2. Override for `MultiBarcodeReadRecordT`:**

The multi-barcode override computes `(barcodes[0] << cell_bc_bits) | barcodes[last]`. This requires `cell_bc_bits` to be available. Options:

- **Option A:** Store `cell_bc_bits` on `MultiBarcodeRecordContext` (the parsing context) and pass it through. The context already stores `bc_types` which encodes the barcode widths.
- **Option B:** Compute it from the barcode type: for `bc_types[last]`, the bit width is `bc_type.bytes_for_type() * 8` (since 2-bit encoding for nucleotides means 4 nucleotides per byte, but the stored type is the integer type, so U32 = 32 bits).
- **Option C:** Add `collation_group_key()` to `CollatableRecordHeader` instead of the record itself, since the header is what's available during the gather's two-pass read. The header already has `barcodes: SmallVec<[B; N]>`.

Option C is probably cleanest since `collate_temporary_bucket_twopass_generic` calls `from_bytes_collatable_header` and uses the header's `collate_key()`. We'd add `collation_group_key()` to `CollatableRecordHeader` (the trait bound, not a concrete type).

Actually, looking more carefully: `CollatableRecordHeader` is an associated type, not a trait. The trait bound on it is just `RecordHeader`. We'd need to add the method there:

```rust
// In record.rs:
pub trait RecordHeader {
    fn naln(&self) -> u32;
    fn collation_group_key(&self) -> u64;  // NEW
    fn write_fields<W: Write>(&self, writer: &mut W, ctx: &...) -> ...;
}
```

With `MultiBarcodeReadRecordHeader` implementing it using the composite key, and other header types defaulting to `collate_key().into()`.

**3. Update `collate_temporary_bucket_twopass_generic`:**

Change line 338 from:
```rust
.entry(tup.collate_key().into())
```
to:
```rust
.entry(tup.collation_group_key())
```

And similarly at line 402 (second pass).

### Changes to alevin-fry after the upstream fix

- Delete `collate_multi_bc_bucket` entirely
- Revert the gather call back to `collate_temporary_bucket_twopass_generic`
- The scatter, quant, and overflow check changes remain as-is

### What stays in alevin-fry regardless

- Scatter-phase logic: writing `sample_idx` into `barcodes[0]`, correcting cell BCs
- Quant-phase logic: reading sample index from records via the extractor closure
- Composite key computation for `tsv_map` and `output_cache` (bucket assignment)
- Manifest construction from composite keys
- Overflow check

## Test to validate the fix

The key property to verify: **no duplicate `sample_cell` barcode strings in the output**, even when cells from different samples share the same cell barcode.

### Automated test

```
1. Use the 4-plex Flex v1 dataset:
   - FASTQs: flex_data_3/4plex_human_colorectal_kidney_scFFPE_multiplex_fastqs/
   - Probe set: flex_data/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv

2. Run the full pipeline (simpleaf multiplex-quant or individual steps):
   - generate-permit-list with --sample-bc-list and --sample-correction-mode exact
   - collate
   - quant with cr-like resolution

3. Check quants_mat_rows.txt:
   - Total rows must equal the sum of num_cells across all samples in sample_info.json
     (595,200 for this dataset with min_reads=10)
   - Number of unique rows must equal total rows (zero duplicates)

4. Compare against CellRanger per-sample output:
   - Per-gene Pearson correlation > 0.9999 for all 4 samples
   - Per-cell Pearson correlation > 0.9999 for all 4 samples
   - These per-cell correlations were ~0.67-0.96 before the fix due to
     cross-sample contamination from merged chunks

5. Verify sample assignment:
   - For each row in quants_mat_rows.txt, the sample prefix must match
     the sample that the manifest assigns to that chunk's position
   - Before the fix, ~2,000 rows had mismatched assignments
```

### Specific regression test for the trait refactor

After upstreaming to libradicl:

```
1. Run the same pipeline as above
2. Verify identical output to the alevin-fry-only fix:
   - Same quants_mat_rows.txt (zero duplicates)
   - Same quants_mat.mtx (identical nnz count and element values)
3. Verify single-barcode data is unaffected:
   - Run a standard 10x Chromium v3 dataset through collate + quant
   - Output must be identical to the pre-refactor version
```
