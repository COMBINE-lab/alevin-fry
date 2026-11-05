"""Compare two alevin-fry quant outputs (quants_mat.mtx + rows) after aligning cells by barcode.
usage: cmp_quants.py <dirA> <dirB>   (each holding alevin/quants_mat.mtx, quants_mat_rows.txt, quants_mat_cols.txt)
"""
import sys
import numpy as np, scipy.io, scipy.sparse as sp


def load(d):
    A = scipy.io.mmread(f"{d}/alevin/quants_mat.mtx").tocsr()
    rows = [l.strip() for l in open(f"{d}/alevin/quants_mat_rows.txt")]
    cols = [l.strip() for l in open(f"{d}/alevin/quants_mat_cols.txt")]
    return A, rows, cols


A, ra, ca = load(sys.argv[1]); B, rb, cb = load(sys.argv[2])
assert ca == cb, "column sets differ"
ia = {b: i for i, b in enumerate(ra)}
common = [b for b in rb if b in ia]
print(f"cells: {len(ra)} vs {len(rb)}, common {len(common)}; columns {len(ca)}")
A2 = A[[ia[b] for b in common]]; B2 = B[[rb.index(b) for b in common]]
D = (A2 - B2).tocsr()
if D.nnz:
    mx = np.abs(D.data).max(); n = (np.abs(D.data) > 1e-6).sum()
else:
    mx = 0.0; n = 0
print(f"nonzeros: {A2.nnz} vs {B2.nnz}; sums {A2.sum():.3f} vs {B2.sum():.3f}")
print(f"max |diff| = {mx:.3g}; entries differing by > 1e-6: {n}")
print("IDENTICAL" if n == 0 and A2.nnz == B2.nnz else "DIFFERENT")
