"""Groups of cells that share isoform usage, found from the one-candidate molecules of a truncation dump
(`--truncation-dump`): bisecting clustering of the cells' usage residuals, each split kept only when it
improves the held-out likelihood of molecules not used to make or to choose it.

Per gene with several isoforms, a cell's one-candidate counts x_cg (over the gene's transcripts) are
multinomial around the shares of its group. Each leaf's molecules are thinned at random into three
parts: the first gives the residual features (x - n phi) / sqrt(n phi (1 - phi)) and the candidate
splits (k-means with several k on the leading singular vectors; each cluster against the rest), the
second ranks the candidates by their held-out likelihood, the third confirms the best one. A side's
shares are its counts shrunk toward the leaf's shares with a strength per gene estimated from how much
the two sides differ beyond sampling (as --truncation-group-shrink auto-gene does), so that a small
side is penalized only on the genes where it does differ. A split is accepted when the confirmation
gain exceeds `min_gain` nats in total and `min_gain_per_cell` per cell. Leaves below `min_cells` are
not split. Without any accepted split every cell is in one group: the pooled model.
`alevin-fry quant --truncation-discover-groups [t2g]` runs this procedure inside the quantification (the Rust
implementation in src/truncation_groups.rs); this script explores the groups from a dump, with a `--truth` check.
usage: truncation_groups.py <dump dir> <t2g> <out groups.tsv> [--truth barcode<TAB>label] ...
  <dump dir>: the directory given to `alevin-fry quant --truncation-dump`; <t2g>: transcript<TAB>gene;
  the output is a barcode<TAB>group file for `--truncation-groups`. Needs numpy, pandas, scipy, scikit-learn.
"""
import sys, argparse, time
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.sparse.linalg import svds
from sklearn.cluster import KMeans

ap = argparse.ArgumentParser()
ap.add_argument("dump"); ap.add_argument("t2g"); ap.add_argument("out")
ap.add_argument("--truth", default=None, help="barcode<TAB>label file to evaluate against")
ap.add_argument("--min-cells", type=int, default=8)
ap.add_argument("--min-gain", type=float, default=20.0, help="nats by which the confirmation gain must beat the best random side of the same size")
ap.add_argument("--null-sides", type=int, default=5, help="random sides of the same size the confirmation gain is compared with")
ap.add_argument("--min-gain-per-cell", type=float, default=0.02)
ap.add_argument("--components", type=int, default=20)
ap.add_argument("--min-gene-molecules", type=int, default=50)
ap.add_argument("--gate-z", type=float, default=4.0, help="null standard deviations a gene's between-side deviation must exceed to get a strength of its own")
ap.add_argument("--gene-prior-weight", type=float, default=250.0, help="pairs of pseudo-molecules pulling a gene's strength toward the common one")
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--max-depth", type=int, default=6)
a = ap.parse_args()
t0 = time.time()
rng = np.random.default_rng(a.seed)
strip = lambda s: s.split(".")[0]

# ---- the dump: cells x transcripts one-candidate counts, restricted to genes with >= 2 isoforms seen
d = pd.read_csv(f"{a.dump}/unique_counts.tsv.gz", sep="\t", dtype={"barcode": str, "transcript": str, "unique": np.int32, "candidates": np.int32})
d = d[d["unique"] > 0]
t2g = pd.read_csv(a.t2g, sep="\t", header=None, usecols=[0, 1], names=["t", "g"], dtype=str)
t2g["t"] = t2g["t"].map(strip); t2g = t2g.drop_duplicates("t").set_index("t")["g"]
d["gene"] = t2g.reindex(d["transcript"].map(strip).values).values
d = d.dropna(subset=["gene"])
cells = np.array(sorted(d["barcode"].unique()))
ci = pd.Index(cells).get_indexer(d["barcode"])
per_t = d.groupby(["gene", "transcript"])["unique"].sum().reset_index()
per_g = per_t.groupby("gene").agg(n_iso=("transcript", "size"), mols=("unique", "sum"))
genes = per_g[(per_g["n_iso"] >= 2) & (per_g["mols"] >= a.min_gene_molecules)].index
per_t = per_t[per_t["gene"].isin(genes)]
txps = per_t["transcript"].values; gene_of = per_t["gene"].values
ti = pd.Index(txps).get_indexer(d["transcript"]); keep = ti >= 0
X = sp.csr_matrix((d["unique"].values[keep].astype(np.float64), (ci[keep], ti[keep])), shape=(len(cells), len(txps)))
gi = pd.Index(pd.unique(gene_of)).get_indexer(gene_of)           # transcript -> gene index
n_genes = gi.max() + 1
G = sp.csr_matrix((np.ones(len(txps)), (np.arange(len(txps)), gi)), shape=(len(txps), n_genes))  # transcript x gene
print(f"cells {len(cells):,}; genes with >= 2 isoforms and >= {a.min_gene_molecules} one-candidate molecules: {len(genes):,}; transcripts {len(txps):,}; molecules {X.sum():,.0f} ({time.time() - t0:.0f}s)", flush=True)


def thin(X, parts):
    """Multinomial thinning of a count matrix into `parts` equal parts."""
    C = X.tocoo()
    draws = rng.multinomial(1, np.full(parts, 1.0 / parts), size=int(C.data.sum())) if C.data.sum() < 5e7 else None
    out = []
    remaining = C.data.astype(np.int64).copy()
    for p in range(parts):
        if p == parts - 1:
            x = remaining
        else:
            x = rng.binomial(remaining, 1.0 / (parts - p))
            remaining = remaining - x
        M = sp.csr_matrix((x.astype(np.float64), (C.row, C.col)), shape=X.shape); M.eliminate_zeros(); out.append(M)
    return out


def colsum(Xs):
    return np.asarray(Xs.sum(0)).ravel()


def plain_shares(x):
    """Each transcript's share of its gene from summed counts x (uniform where the gene has none)."""
    n = G.T @ x
    uni = 1.0 / np.asarray(G.sum(0)).ravel()[gi]
    return np.where(n[gi] > 0, x / np.maximum(n[gi], 1e-300), uni)


def eb_shares(x_side, x_rest, phi):
    """Shares of each side shrunk toward phi with a strength per gene from how much the sides differ
    beyond sampling: method of moments over the two sides with the pooled share plugged in (so that
    a unit's deviation counts with weight 1 - 2 n/N + S2/N^2), shrunk toward the common strength
    with `gene_prior_weight` pairs; returns (shares_side, shares_rest, beta per gene)."""
    n_s = (G.T @ x_side)[gi]; n_r = (G.T @ x_rest)[gi]
    v = phi * (1 - phi)
    N = n_s + n_r
    S2_all = n_s ** 2 + n_r ** 2
    dev = np.zeros(len(phi)); sn = np.zeros(len(phi)); sn2 = np.zeros(len(phi)); sn3 = np.zeros(len(phi))
    for x, n in ((x_side, n_s), (x_rest, n_r)):
        ok = n >= 2
        dev += np.where(ok, (x - n * phi) ** 2, 0.0)
        sn += np.where(ok, n, 0.0); sn2 += np.where(ok, n * n, 0.0); sn3 += np.where(ok, n ** 3, 0.0)
    Npos = np.maximum(N, 1e-300)
    e = sn - sn2 / Npos
    with np.errstate(invalid="ignore", divide="ignore"):
        D = np.where(N > 0, (sn2 - sn) - 2 * (sn3 - sn2) / Npos + sn2 * (S2_all - N) / Npos ** 2, 0.0)
    exc = dev - v * e
    good = (D > 0) & (N > 0)
    num_g = np.bincount(gi, weights=np.where(good, D * v, 0.0), minlength=n_genes)
    den_g = np.bincount(gi, weights=np.where(good, exc, 0.0), minlength=n_genes)
    # the gate: a gene's own strength only where its deviation is beyond sampling (null variance
    # about 2 D sum_t v_t^2); the common strength is that of the gated genes
    null_var = np.bincount(gi, weights=np.where(good, 2 * D * v * v, 0.0), minlength=n_genes)
    use = (num_g > 0) & (den_g > a.gate_z * np.sqrt(null_var))
    common = den_g[use].sum() / max(num_g[use].sum(), 1e-300) if use.any() else 0.0
    r = (den_g + a.gene_prior_weight * common) / (num_g + a.gene_prior_weight)
    beta = np.where(r > 0, np.maximum(1.0 / np.maximum(r, 1e-300) - 1.0, 1.0), np.inf)   # never below one pseudo-molecule
    beta = np.where(use, beta, np.inf)
    b = np.minimum(beta[gi], 1e9)
    s_side = (x_side + b * phi) / np.maximum(n_s + b, 1e-300)
    s_rest = (x_rest + b * phi) / np.maximum(n_r + b, 1e-300)
    return s_side, s_rest, beta


def loglik_rows(Xs, phi):
    return np.asarray(Xs @ np.log(np.maximum(phi, 1e-300)).reshape(-1, 1)).ravel()


def variable_genes(Xs, phi):
    """Genes whose isoform shares vary across the rows of Xs beyond sampling: the pairs-weighted
    excess of sum_c (x - n phi)^2 over its multinomial expectation, as a z-score against the
    spread of the per-cell terms; returns a boolean mask over transcripts (True: keep)."""
    n = (Xs @ G).tocsr()[:, gi]                                     # cells x transcripts, n_cg
    E = sp.csr_matrix(n.multiply(phi.reshape(1, -1)))
    V = E.multiply((1 - phi).reshape(1, -1))
    term = ((Xs - E).power(2) - V).tocsr()                          # per cell and transcript, mean 0 under no variation
    exc_t = np.asarray(term.sum(0)).ravel()
    sd_t = np.sqrt(np.asarray(term.power(2).sum(0)).ravel() + 1e-12)
    exc_g = np.bincount(gi, weights=exc_t, minlength=n_genes); sd_g = np.sqrt(np.bincount(gi, weights=sd_t ** 2, minlength=n_genes))
    z = exc_g / np.maximum(sd_g, 1e-12)
    keep_g = z > 2.0
    if keep_g.sum() < 10:
        keep_g = np.ones(n_genes, bool)
    return keep_g[gi], int(keep_g.sum())


def residuals(Xs, phi, keep=None):
    """Pearson residuals of the one-candidate counts around shares phi (sparse where n_cg = 0),
    restricted to the transcripts in `keep`."""
    n = (Xs @ G).tocsr()
    E = sp.csr_matrix(n[:, gi].multiply(phi.reshape(1, -1)))
    V = E.multiply((1 - phi).reshape(1, -1)).tocsr()
    R = (Xs - E).tocsr()
    R.data /= np.sqrt(np.asarray(V[R.nonzero()]).ravel() + 0.25)
    R.data = np.clip(R.data, -10, 10)
    if keep is not None:
        R = R[:, np.where(keep)[0]]
    return R


def gain_of(side, Xa, Xb):
    """Held-out gain of splitting the leaf into side / rest: shares from Xa (EB-shrunk toward the
    leaf's), likelihood on Xb, against the leaf's shares from Xa."""
    xa = colsum(Xa); phi = plain_shares(xa)
    s_in, s_out, beta = eb_shares(colsum(Xa[side]), colsum(Xa[~side]), phi)
    ll_leaf = loglik_rows(Xb, phi)
    ll_split = np.where(side, loglik_rows(Xb, s_in), loglik_rows(Xb, s_out))
    return float((ll_split - ll_leaf).sum()), int(np.isfinite(beta).sum())


def candidates(Z, n, min_cells):
    """Candidate sides from an embedding Z: k-means clusters for several k, agglomerative (Ward)
    clusters for several k, and the tails of each component; each as a boolean mask (cluster vs rest)."""
    from sklearn.cluster import AgglomerativeClustering
    seen = set(); out = []
    def add(mask, how):
        if min(mask.sum(), (~mask).sum()) < min_cells:
            return
        key = mask.tobytes()
        if key in seen:
            return
        seen.add(key); out.append((mask, how))
    for kk in [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]:
        if kk > n // min_cells:
            break
        lab = KMeans(n_clusters=kk, n_init=10, random_state=a.seed).fit(Z).labels_
        for c in range(kk):
            add(lab == c, f"k-means k={kk}")
        if n <= 6000:
            lab = AgglomerativeClustering(n_clusters=kk, linkage="ward").fit(Z).labels_
            for c in range(kk):
                add(lab == c, f"ward k={kk}")
    for j in range(Z.shape[1]):
        z = (Z[:, j] - np.median(Z[:, j])) / (1.4826 * np.median(np.abs(Z[:, j] - np.median(Z[:, j]))) + 1e-12)
        for sign in (1, -1):
            for cut in (3.0, 5.0):
                add(sign * z > cut, f"component {j} tail {sign * cut:+.0f}")
    return out


def split(rows, depth):
    """Try to split the cells `rows`; returns (left, right, gain) or None."""
    if len(rows) < 2 * a.min_cells or depth >= a.max_depth:
        return None
    X1, X2, X3, X4 = thin(X[rows], 4)
    Xf = X1 + X2                                  # features and the shares the selection scores with
    phi_f = plain_shares(colsum(Xf))
    keep, n_var = variable_genes(Xf, phi_f)
    R = residuals(Xf, phi_f, keep)
    k = min(a.components, min(R.shape) - 1)
    if k < 2 or R.nnz == 0:
        return None
    U, S, _ = svds(R, k=k, random_state=a.seed)
    Z = U * S
    best = None
    for side, how in candidates(Z, len(rows), a.min_cells):
        g, _ = gain_of(side, Xf, X3)
        if best is None or g > best[1]:
            best = (side, g, how)
    if best is None:
        return None
    side = best[0]
    # confirmation on molecules that neither made nor chose the split: shares from parts 1-3, likelihood on part 4
    gain, genes_finite = gain_of(side, Xf + X3, X4)
    # the null: random sides of the same size, scored the same way (diffuse variation between
    # cells, or overfitting, raises their gains too; a real group must beat them clearly)
    null = []
    for _ in range(a.null_sides):
        fake = np.zeros(len(rows), bool); fake[rng.choice(len(rows), int(side.sum()), replace=False)] = True
        null.append(gain_of(fake, Xf + X3, X4)[0])
    null_max = max(null)
    per_cell = (gain - null_max) / len(rows)
    ok = gain > null_max + a.min_gain and per_cell > a.min_gain_per_cell
    print(f"  depth {depth}: {len(rows)} cells, {n_var} variable genes -> {side.sum()} + {(~side).sum()} ({best[2]}; selection gain {best[1]:.1f}), "
          f"confirmation gain {gain:.1f} nats against random sides {null_max:.1f} at most ({per_cell:.3f}/cell above), {genes_finite} genes with a finite strength: {'ACCEPT' if ok else 'reject'}", flush=True)
    return (rows[side], rows[~side], gain) if ok else None


leaves = [(np.arange(len(cells)), 0)]
final = []
while leaves:
    rows, depth = leaves.pop()
    s = split(rows, depth)
    if s is None:
        final.append(rows)
    else:
        leaves.extend([(s[0], depth + 1), (s[1], depth + 1)])
labels = np.zeros(len(cells), int)
for k, rows in enumerate(final):
    labels[rows] = k
K = len(final)
if K > 1:
    # refinement: each cell to the group whose (shrunk) shares fit all its molecules best
    x_all = colsum(X); phi_all = plain_shares(x_all)
    for it in range(15):
        LL = np.zeros((len(cells), K))
        for k in range(K):
            s_in, _, _ = eb_shares(colsum(X[labels == k]), colsum(X[labels != k]), phi_all)
            LL[:, k] = loglik_rows(X, s_in)
        new = LL.argmax(1)
        for k in range(K):
            if (new == k).sum() < a.min_cells:
                new[labels == k] = k
        moved = int((new != labels).sum())
        labels = new
        print(f"  refinement {it}: {moved} cells moved", flush=True)
        if moved == 0:
            break
sizes = np.bincount(labels, minlength=K)
print(f"groups: {K} ({', '.join(str(s) for s in sorted(sizes, reverse=True))}) ({time.time() - t0:.0f}s)")
pd.DataFrame({"barcode": cells, "group": [f"L{k}" for k in labels]}).to_csv(a.out, sep="\t", index=False, header=False)
if a.truth:
    tr = pd.read_csv(a.truth, sep="\t", header=None, names=["bc", "lab"], dtype=str).drop_duplicates("bc").set_index("bc")["lab"]
    lab = tr.reindex(cells).fillna("?").values
    from sklearn.metrics import adjusted_rand_score
    print(f"adjusted Rand index vs truth: {adjusted_rand_score(lab, labels):.3f}")
    print(pd.crosstab(pd.Series(lab, name="truth"), pd.Series([f"L{k}" for k in labels], name="found")))
