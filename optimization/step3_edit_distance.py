"""Step 3: bounded UMI edit distance (pigeonhole filter + banded DP with early exit) in both UMI-graph builders.

The deduplication never needs a distance above its bound (3 in the calibrated edge model; the far
bound in the rule model), so the distance is evaluated only up to that bound. Decisions are unchanged.
"""
import os, shutil

ROOT = "/scratch1/zahra/source-code/alevin-fry/src"
shutil.copy("/scratch5/zahra/tmp/optim/umi_general_distance_v2.rs", os.path.join(ROOT, "umi_general_distance.rs"))


def edit(path, pairs):
    p = os.path.join(ROOT, path)
    s = open(p).read()
    for old, new, count in pairs:
        n = s.count(old)
        assert n == count, f"{path}: expected {count} of {old[:70]!r}, found {n}"
        s = s.replace(old, new)
    open(p, "w").write(s)


edit("pugutils.rs", [
    ("use crate::umi_general_distance::umi_edit_distance_from_packed_shifted;",
     "use crate::umi_general_distance::umi_edit_distance_within;", 1),
    # --- the rule model (extract_graph): any distance of far_bound or more is no edge
    ("""    let umi_len: u8 = umi_len();
    let mut has_edge = |ex: u32, x: &(u64, u32), ey: u32, y: &(u64, u32)| -> PugEdgeType {
        let hdist = if pug_exact_umi {
            if x.0 == y.0 { 0 } else { usize::MAX }
        } else {
            if ed_distance_model == EditDistanceMode::General {
                umi_edit_distance_from_packed_shifted(x.0, y.0, umi_len)
            } else if ed_distance_model == EditDistanceMode::Levenshtein {""",
     """    let umi_len: u8 = umi_len();
    // distances of `far_bound` or more all mean "no edge", so the general distance is only evaluated up to far_bound - 1
    let general_bound: u32 = far_bound.saturating_sub(1) as u32;
    let mut has_edge = |ex: u32, x: &(u64, u32), ey: u32, y: &(u64, u32)| -> PugEdgeType {
        let hdist = if pug_exact_umi {
            if x.0 == y.0 { 0 } else { usize::MAX }
        } else {
            if ed_distance_model == EditDistanceMode::General {
                umi_edit_distance_within(x.0, y.0, umi_len, general_bound).map_or(usize::MAX, |d| d as usize)
            } else if ed_distance_model == EditDistanceMode::Levenshtein {""", 1),
    # --- the calibrated edge model (extract_graph_bayes): distances above MAX_D are "different molecules"
    ("""    let dist = |x: u64, y: u64| -> usize {
        match ed_distance_model {
            EditDistanceMode::General => umi_edit_distance_from_packed_shifted(x, y, umi_len),
            EditDistanceMode::Levenshtein => umi_edit_distance_from_packed(x, y, umi_len),
            EditDistanceMode::Hamming => afutils::count_diff_2_bit_packed(x, y),
        }
    };""",
     """    // `Some(d)` for a distance of at most MAX_D, `None` beyond it (the exact value is not needed there)
    let dist = |x: u64, y: u64| -> Option<usize> {
        match ed_distance_model {
            EditDistanceMode::General => umi_edit_distance_within(x, y, umi_len, MAX_D as u32).map(|d| d as usize),
            EditDistanceMode::Levenshtein => Some(umi_edit_distance_from_packed(x, y, umi_len)).filter(|d| *d <= MAX_D),
            EditDistanceMode::Hamming => Some(afutils::count_diff_2_bit_packed(x, y)).filter(|d| *d <= MAX_D),
        }
    };""", 1),
    ("""        ($ex:expr, $xi:expr, $x:expr, $ey:expr, $yi:expr, $y:expr) => {{
            let d = dist($x.0, $y.0);
            if d == 0 {
                graph.add_edge(($ex, $xi), ($ey, $yi), ());
                graph.add_edge(($ey, $yi), ($ex, $xi), ());
            } else if d <= MAX_D {
                let a = agree($ex, $x.0, $ey, $y.0);
                let crowd = stratum_of(eqmap.eqc_info[$ex as usize].umis.len().max(eqmap.eqc_info[$ey as usize].umis.len()));
                cands.push(($ex, $xi, $ey, $yi, d as u8, a, $x.1, $y.1, crowd));
            } else {""",
     """        ($ex:expr, $xi:expr, $x:expr, $ey:expr, $yi:expr, $y:expr) => {{
            let d = dist($x.0, $y.0);
            if d == Some(0) {
                graph.add_edge(($ex, $xi), ($ey, $yi), ());
                graph.add_edge(($ey, $yi), ($ex, $xi), ());
            } else if let Some(d) = d {
                let a = agree($ex, $x.0, $ey, $y.0);
                let crowd = stratum_of(eqmap.eqc_info[$ex as usize].umis.len().max(eqmap.eqc_info[$ey as usize].umis.len()));
                cands.push(($ex, $xi, $ey, $yi, d as u8, a, $x.1, $y.1, crowd));
            } else {""", 1),
])
print("step 3 edits applied")
