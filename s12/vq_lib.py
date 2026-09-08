"""s12 QUANTUM-ROLE agent -- the combinatorial instance family and its Hamiltonian.

THE PROBLEM (VQ-A): ANCHORED MULTI-SEGMENT FRAGMENT ASSEMBLY
============================================================
A target of length n is cut into k contiguous segments.  For each segment s a small set of
m structurally distinct library pieces is retrieved (BLOSUM62 on the sub-sequence, then
agglomerative clustering to m medoids so the m choices are genuinely different shapes).
Each piece is PLACED by superposing it onto the corresponding residues of a deployable
ANCHOR structure (the shipped synthesis `fit_ca` -- no native anywhere).  An assembly is a
choice of one piece per segment; the emitted structure is the concatenation of the placed
pieces.

The objective is the SHIPPED leave-fold-out distogram Bayes risk of the assembled CA trace,
bit-identical to `I.shipped_score`:

    F(x) = (1/npairs) sum_{p=(i,j)} risk_p( || X_i(x) - X_j(x) || )

Because a placed piece's coordinates depend only on that segment's choice, a pair (i,j)
with i in segment s and j in segment t depends on x_s and x_t ONLY.  Therefore

    npairs * F(x) = sum_s h_s(x_s) + sum_{s<t} J_st(x_s, x_t)                 (EXACT)

with no approximation of any kind: h_s is the sum of risk over the pairs internal to s and
J_st the sum over the pairs that straddle s and t.  `verify_decomposition` asserts this to
machine precision against a from-scratch recomputation of `I.shipped_score`.

ENCODINGS
---------
one-hot   qubit (s,a) = 1 iff segment s takes piece a;  k*m qubits;  a genuine 2-local
          Ising Hamiltonian once the one-hot penalty  P * sum_s (sum_a x_sa - 1)^2  is
          expanded.  Infeasible basis states exist (2^(km) states, m^k feasible).
index     ceil(log2 m) qubits per segment address the piece directly; no constraint, the
          whole register is feasible, but the energy is NOT 2-local in the qubits.

Both are exercised: the one-hot form is the faithful QUBO a VQE is built for, the index
form is the one that uses the register efficiently.
"""
from __future__ import annotations
import os, sys, math, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I
from s12 import asm_lib as A

CACHE = os.path.join(ROOT, "s12", "cache", "vq")
os.makedirs(CACHE, exist_ok=True)


# --------------------------------------------------------------------------- geometry
def kabsch_place(P, T):
    """Superpose each (L,3) in P (m,L,3) onto T (L,3).  Returns the moved copies."""
    return I.superpose_batch(np.asarray(P, float), np.asarray(T, float))


def segments(n, k, minlen=4):
    """k contiguous segments of [0,n), as equal as possible, each >= minlen."""
    if k * minlen > n:
        return None
    base = n // k
    rem = n - base * k
    out, s = [], 0
    for q in range(k):
        L = base + (1 if q < rem else 0)
        out.append((s, s + L))
        s += L
    assert out[-1][1] == n and min(b - a for a, b in out) >= minlen
    return out


# --------------------------------------------------------------------------- risk
def pair_risk(dg, D):
    """Per-pair Bayes risk for pair-distance rows D (b, npairs) -- the shipped lookup,
    element for element (`I.shipped_score` is this .mean(1))."""
    grid = dg["grid"]; risk = dg["risk"]
    g = np.clip(((np.asarray(D, float) - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    return risk[np.arange(risk.shape[0])[None, :], g]


# --------------------------------------------------------------------------- instance
def retrieve_pieces(fold, L, sub_codes, n_ret=200, m=8, seed=0):
    """m structurally distinct library pieces of length L for a sub-sequence.

    Top-`n_ret` by BLOSUM62 against `sub_codes` (deployable retrieval), then agglomerative
    average-linkage clustering of the n_ret x n_ret CA-RMSD matrix into m clusters; the
    medoid of each cluster is kept.  This makes the m choices span the retrieved set's
    structural diversity instead of duplicating its top entry -- the same device the
    aggregation agent used to build its 16 hypotheses.
    """
    b = A.bank(int(fold), int(L))
    B = A.B62()
    sim = B[np.asarray(sub_codes, int)[None, :], b["S"].astype(int)].sum(1)
    n_ret = min(n_ret, len(sim))
    top = np.argsort(-sim, kind="stable")[:n_ret]
    W = b["W"][top]                                   # (n_ret, L, 3)
    if n_ret <= m:
        return W, top
    Wc = W - W.mean(1, keepdims=True)
    P = A.rmsd_many(Wc, (Wc ** 2).sum((1, 2)), W)     # (n_ret, n_ret), analytic Kabsch
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    Z = linkage(squareform((P + P.T) / 2.0, checks=False), method="average")
    lab = fcluster(Z, t=m, criterion="maxclust")
    keep = []
    for c in np.unique(lab):
        idx = np.where(lab == c)[0]
        keep.append(idx[int(np.argmin(P[np.ix_(idx, idx)].mean(1)))])
    keep = np.array(sorted(keep)[:m], int)
    while len(keep) < m:                              # fcluster can return < m clusters
        extra = [t for t in range(n_ret) if t not in set(keep.tolist())]
        keep = np.array(sorted(keep.tolist() + extra[: m - len(keep)]), int)
    return W[keep], top[keep]


def make_instance(t, k=2, m=8, n_ret=200, minlen=4, anchor="fit_ca"):
    """Build one VQ-A instance.  Returns a dict with the exact 2-local tables.

    Nothing here touches the native: retrieval is BLOSUM on the target sequence, the
    anchor is the shipped emitted structure, the objective is the shipped distogram.
    """
    pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
    segs = segments(n, k, minlen)
    if segs is None:
        return None
    rec = I.shipped_record(pdb)
    T = np.asarray(rec[anchor], float)                # (n,3) deployable anchor
    dg = I.distogram(pdb, seq, fold)
    pi, pj = I.pair_index(n)
    npairs = len(pi)
    codes = np.array([I.ALPHABET.index(c) for c in seq], int)

    placed, src = [], []
    for (a, b_) in segs:
        W, top = retrieve_pieces(fold, b_ - a, codes[a:b_], n_ret=n_ret, m=m)
        placed.append(kabsch_place(W, T[a:b_]))       # (m, L, 3)
        src.append(top)
    seg_of = np.empty(n, int)
    for s, (a, b_) in enumerate(segs):
        seg_of[a:b_] = s

    # ---- exact 2-local tables -------------------------------------------------------
    # a pair is internal to s (both endpoints in s) or straddles (s,t)
    ss, tt = seg_of[pi], seg_of[pj]
    h = [np.zeros(m) for _ in range(k)]
    J = {}
    for s in range(k):
        msk = (ss == s) & (tt == s)
        if msk.any():
            a0 = segs[s][0]
            d = np.linalg.norm(placed[s][:, pi[msk] - a0, :] - placed[s][:, pj[msk] - a0, :], axis=-1)
            grid = dg["grid"]; risk = np.asarray(dg["risk"][msk], np.float64)
            g = np.clip(((d - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
            h[s] = risk[np.arange(msk.sum())[None, :], g].sum(1)
    for s in range(k):
        for u in range(s + 1, k):
            msk = (ss == s) & (tt == u)
            if not msk.any():
                J[(s, u)] = np.zeros((m, m)); continue
            a0, b0 = segs[s][0], segs[u][0]
            X = placed[s][:, pi[msk] - a0, :]         # (m, q, 3)
            Y = placed[u][:, pj[msk] - b0, :]         # (m, q, 3)
            d = np.linalg.norm(X[:, None, :, :] - Y[None, :, :, :], axis=-1)   # (m,m,q)
            grid = dg["grid"]; risk = np.asarray(dg["risk"][msk], np.float64)
            g = np.clip(((d - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
            J[(s, u)] = risk[np.arange(msk.sum())[None, None, :], g].sum(2)
    return dict(pdb=pdb, n=n, fold=fold, seq=seq, k=k, m=m, segs=segs, npairs=npairs,
                placed=placed, h=h, J=J, anchor=T, dg=dg, pi=pi, pj=pj, src=src,
                seg_of=seg_of)


# --------------------------------------------------------------------------- objective
def assemble(inst, cfg):
    """Emitted CA point cloud for configuration cfg (length-k piece indices)."""
    X = np.empty((inst["n"], 3))
    for s, (a, b_) in enumerate(inst["segs"]):
        X[a:b_] = inst["placed"][s][int(cfg[s])]
    return X


def assemble_batch(inst, cfgs):
    """(b,k) configurations -> (b,n,3) assembled point clouds, vectorised."""
    cfgs = np.atleast_2d(np.asarray(cfgs, int))
    X = np.empty((len(cfgs), inst["n"], 3))
    for s, (a, b_) in enumerate(inst["segs"]):
        X[:, a:b_] = inst["placed"][s][cfgs[:, s]]
    return X


def objective_true(inst, cfgs):
    """Ground truth: rebuild each assembly and score it with `I.shipped_score`."""
    cfgs = np.atleast_2d(np.asarray(cfgs, int))
    X = np.stack([assemble(inst, c) for c in cfgs])
    D = I.pair_dists(X, inst["pi"], inst["pj"])
    return I.shipped_score(inst["dg"], D)


def objective_tables(inst, cfgs):
    """The 2-local Hamiltonian evaluated on cfgs -- must equal `objective_true`."""
    cfgs = np.atleast_2d(np.asarray(cfgs, int))
    out = np.zeros(len(cfgs))
    for s in range(inst["k"]):
        out += inst["h"][s][cfgs[:, s]]
    for (s, u), Jm in inst["J"].items():
        out += Jm[cfgs[:, s], cfgs[:, u]]
    return out / inst["npairs"]


def enumerate_all(inst):
    """All m^k configurations and their exact energies (small instances only)."""
    k, m = inst["k"], inst["m"]
    grids = np.meshgrid(*[np.arange(m)] * k, indexing="ij")
    cfgs = np.stack([g.ravel() for g in grids], 1)
    return cfgs, objective_tables(inst, cfgs)


def objective_true64(inst, cfgs):
    """`objective_true` with the risk table promoted to float64 before the mean.

    `I.shipped_score` accumulates the 78-760 element risk lookup in float32; the table
    Hamiltonian accumulates block sums in float64.  In exact arithmetic they are the same
    number, so the honest verification compares like with like AND separately asserts that
    the two routes address IDENTICAL risk-grid bins (`bins_identical`).
    """
    cfgs = np.atleast_2d(np.asarray(cfgs, int))
    X = np.stack([assemble(inst, c) for c in cfgs])
    D = I.pair_dists(X, inst["pi"], inst["pj"])
    dg = inst["dg"]; grid = dg["grid"]
    g = np.clip(((D - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    risk = np.asarray(dg["risk"], np.float64)
    return risk[np.arange(risk.shape[0])[None, :], g].mean(1)


def bins_identical(inst, cfgs):
    """True iff the risk-grid bin of every pair is the same under both routes.

    The table route stores, per block, the bins of the distances between PLACED pieces;
    the truth route recomputes them from the assembled cloud.  Identical bins + exact
    block partition of the pair list is the whole content of the decomposition claim.
    """
    cfgs = np.atleast_2d(np.asarray(cfgs, int))
    X = np.stack([assemble(inst, c) for c in cfgs])
    D = I.pair_dists(X, inst["pi"], inst["pj"])
    grid = inst["dg"]["grid"]
    g_true = np.clip(((D - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    # rebuild the same bins from the per-block placed coordinates
    g_blk = np.empty_like(g_true)
    ss, tt = inst["seg_of"][inst["pi"]], inst["seg_of"][inst["pj"]]
    for s in range(inst["k"]):
        for u in range(s, inst["k"]):
            msk = (ss == s) & (tt == u)
            if not msk.any():
                continue
            a0, b0 = inst["segs"][s][0], inst["segs"][u][0]
            P = inst["placed"][s][cfgs[:, s]][:, inst["pi"][msk] - a0, :]
            Qq = inst["placed"][u][cfgs[:, u]][:, inst["pj"][msk] - b0, :]
            d = np.linalg.norm(P - Qq, axis=-1)
            g_blk[:, msk] = np.clip(((d - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    return bool(np.array_equal(g_true, g_blk))


def verify_decomposition(inst, n_rand=64, seed=0):
    """|tables - truth| over random configurations, plus the exact-bin assertion."""
    rng = np.random.default_rng(seed)
    cfgs = rng.integers(0, inst["m"], size=(n_rand, inst["k"]))
    a = objective_tables(inst, cfgs)
    b32 = objective_true(inst, cfgs)
    b64 = objective_true64(inst, cfgs)
    return dict(max_abs_vs_shipped_f32=float(np.abs(a - b32).max()),
                max_abs_vs_f64=float(np.abs(a - b64).max()),
                max_rel_vs_f64=float((np.abs(a - b64) / np.abs(b64)).max()),
                bins_identical=bins_identical(inst, cfgs))


# --------------------------------------------------------------------------- encodings
def onehot_ising(inst, penalty=None):
    """One-hot QUBO/Ising:  Q (N,N) upper-triangular + linear q (N,), N = k*m.

    E(x) = sum_i q_i x_i + sum_{i<j} Q_ij x_i x_j  with the one-hot penalty folded in.
    Returned in x in {0,1}; `ising_from_qubo` converts to +-1 spins.
    """
    k, m = inst["k"], inst["m"]
    N = k * m
    q = np.zeros(N); Q = np.zeros((N, N))
    idx = lambda s, a: s * m + a
    for s in range(k):
        for a in range(m):
            q[idx(s, a)] += inst["h"][s][a] / inst["npairs"]
    for (s, u), Jm in inst["J"].items():
        for a in range(m):
            for b_ in range(m):
                i, j = idx(s, a), idx(u, b_)
                lo, hi = (i, j) if i < j else (j, i)
                Q[lo, hi] += Jm[a, b_] / inst["npairs"]
    if penalty is None:
        if m ** k <= 200000:
            _, E = enumerate_all(inst)
            span = float(E.max() - E.min())
        else:
            span = 1.0
        penalty = 2.0 * max(span, 1e-6)
    # P * sum_s (sum_a x_sa - 1)^2  = P * sum_s [ -sum_a x_sa + 2 sum_{a<b} x_sa x_sb ] + P*k
    for s in range(k):
        for a in range(m):
            q[idx(s, a)] -= penalty
            for b_ in range(a + 1, m):
                Q[idx(s, a), idx(s, b_)] += 2.0 * penalty
    return dict(q=q, Q=Q, offset=penalty * k, penalty=penalty, N=N, k=k, m=m)


def qubo_energy(H, X):
    """E for bit rows X (b,N)."""
    X = np.atleast_2d(np.asarray(X, float))
    return X @ H["q"] + np.einsum("bi,ij,bj->b", X, H["Q"], X) + H["offset"]


def ising_from_qubo(H):
    """x = (1 - z)/2  ->  E(z) = c + sum h_i z_i + sum_{i<j} Jz_ij z_i z_j."""
    q, Q, N = H["q"], H["Q"], H["N"]
    Am = Q + Q.T                             # symmetric, zero diagonal; A_ij = Q_ij (i<j)
    r = Am.sum(1)
    hz = -q / 2.0 - r / 4.0
    Jz = np.triu(Am, 1) / 4.0
    c = H["offset"] + q.sum() / 2.0 + Am.sum() / 8.0
    return dict(h=hz, J=Jz, const=c, N=N)


def ising_energy(Iz, Z):
    Z = np.atleast_2d(np.asarray(Z, float))
    return Iz["const"] + Z @ Iz["h"] + np.einsum("bi,ij,bj->b", Z, np.triu(Iz["J"], 1), Z)


def energy_vector_onehot(inst, H):
    """Full 2^N energy vector for the one-hot register (N = k*m <= ~16)."""
    N = H["N"]
    B = ((np.arange(1 << N)[:, None] >> np.arange(N - 1, -1, -1)[None, :]) & 1).astype(float)
    return qubo_energy(H, B), B


def energy_vector_index(inst):
    """Full 2^(k*bits) energy vector for the INDEX encoding (bits = ceil(log2 m)).

    Piece indices >= m are wrapped modulo m so every basis state is feasible.
    """
    k, m = inst["k"], inst["m"]
    bits = int(math.ceil(math.log2(m)))
    N = k * bits
    codes = np.arange(1 << N)
    cfg = np.empty((1 << N, k), int)
    for s in range(k):
        sh = (k - 1 - s) * bits
        cfg[:, s] = ((codes >> sh) & ((1 << bits) - 1)) % m
    return objective_tables(inst, cfg), cfg, N


# --------------------------------------------------------------------------- caching
def instance_path(pdb, k, m, n_ret, minlen, anchor):
    return os.path.join(CACHE, f"inst_{pdb}_k{k}_m{m}_r{n_ret}_L{minlen}_{anchor}.npz")


def cached_instance(t, k=2, m=8, n_ret=200, minlen=4, anchor="fit_ca"):
    """`make_instance` with the expensive parts (retrieval, clustering, tables) on disk."""
    p = instance_path(t["pdb"], k, m, n_ret, minlen, anchor)
    if os.path.exists(p):
        z = np.load(p, allow_pickle=False)
        kk = int(z["k"]); mm = int(z["m"]); n = int(z["n"])
        segs = [tuple(int(x) for x in r) for r in z["segs"]]
        inst = dict(pdb=str(t["pdb"]), n=n, fold=int(t["fold"]), seq=t["seq"], k=kk, m=mm,
                    segs=segs, npairs=int(z["npairs"]),
                    placed=[z[f"placed{s}"] for s in range(kk)],
                    h=[z[f"h{s}"] for s in range(kk)],
                    J={(int(a), int(b)): z[f"J{a}_{b}"] for a, b in z["Jkeys"]},
                    anchor=z["anchor"], src=[z[f"src{s}"] for s in range(kk)],
                    seg_of=z["seg_of"], pi=z["pi"], pj=z["pj"])
        inst["dg"] = I.distogram(t["pdb"], t["seq"], t["fold"])
        return inst
    inst = make_instance(t, k=k, m=m, n_ret=n_ret, minlen=minlen, anchor=anchor)
    if inst is None:
        return None
    out = dict(k=inst["k"], m=inst["m"], n=inst["n"], npairs=inst["npairs"],
               segs=np.array(inst["segs"], int), anchor=inst["anchor"],
               seg_of=inst["seg_of"], pi=inst["pi"], pj=inst["pj"],
               Jkeys=np.array(sorted(inst["J"].keys()), int).reshape(-1, 2))
    for s in range(inst["k"]):
        out[f"placed{s}"] = inst["placed"][s]
        out[f"h{s}"] = inst["h"][s]
        out[f"src{s}"] = inst["src"][s]
    for (a, b_), Jm in inst["J"].items():
        out[f"J{a}_{b_}"] = Jm
    np.savez_compressed(p, **out)
    return inst
