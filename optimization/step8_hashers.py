"""Step 8: fixed-seed ahash for the per-cell lookup maps of the edge model.

`alignment_intervals` (one entry per (equivalence class, UMI) vertex), the vertex position ranges
and the pair posteriors are only ever looked up by key, or iterated to build another lookup map,
so the hasher cannot change any result; SipHash (std's default) is several times slower than ahash
on these small keys and these maps are probed millions of times per cell.
"""
import os
ROOT = "/scratch1/zahra/source-code/alevin-fry/src"


def edit(path, pairs):
    p = os.path.join(ROOT, path)
    s = open(p).read()
    for old, new, count in pairs:
        n = s.count(old)
        assert n == count, f"{path}: expected {count} of {old[:70]!r}, found {n}"
        s = s.replace(old, new)
    open(p, "w").write(s)


edit("pugutils.rs", [
    ("pub type PairPosteriors = HashMap<((u32, u32), (u32, u32)), f32>;",
     "pub type PairPosteriors = HashMap<((u32, u32), (u32, u32)), f32, ahash::RandomState>;\n\n/// The hasher of the per-cell lookup maps: fixed seeds, so a run is reproducible.\nfn lookup_hasher() -> ahash::RandomState {\n    ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64)\n}", 1),
    ("type VertexRanges = HashMap<(u32, u64), SmallVec<[(u32, u32, u32, u32, u32); 4]>>;",
     "type VertexRanges = HashMap<(u32, u64), SmallVec<[(u32, u32, u32, u32, u32); 4]>, ahash::RandomState>;", 1),
    ("    let mut out: VertexRanges = HashMap::new();", "    let mut out: VertexRanges = HashMap::with_hasher(lookup_hasher());", 1),
    ("    let mut posteriors: PairPosteriors = HashMap::new();", "    let mut posteriors: PairPosteriors = HashMap::with_hasher(lookup_hasher());", 1),
])
edit("eq_class.rs", [
    ("    pub alignment_intervals: Option<HashMap<(u32, u64), Vec<AlignmentInterval>>>,",
     "    pub alignment_intervals: Option<HashMap<(u32, u64), Vec<AlignmentInterval>, ahash::RandomState>>,", 1),
    ("                .get_or_insert_with(HashMap::new)",
     "                .get_or_insert_with(|| HashMap::with_hasher(ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64)))", 1),
])
print("step 8: ahash lookup maps")
