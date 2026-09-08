"""s12 QUANTUM-ROLE -- classical solvers for the VQ-A assembly QUBO and the VQ-B subset QUBO.

Every solver here works on the SAME tables the quantum arm sees, counts its objective
evaluations, and reports wall time, so the comparison in `vq_FINDINGS.md` is at matched
budget rather than matched effort-in-words.

VQ-A is a categorical problem (one choice of m per segment); the natural neighbourhood is
"change one segment's piece".  VQ-B is a binary subset problem; the natural neighbourhood is
flip / swap.  Both families get: exhaustive (ground truth), greedy, greedy + local search,
simulated annealing, beam search, and -- for VQ-A -- an exact integer program via
`scipy.optimize.milp` on the standard Glover linearisation of the quadratic assignment.
"""
from __future__ import annotations
import time
import numpy as np


class Counter:
    def __init__(self):
        self.n = 0
        self.t0 = time.perf_counter()

    def add(self, k=1):
        self.n += int(k)

    def done(self, cfg, val):
        return dict(cfg=np.asarray(cfg).tolist(), value=float(val), evals=int(self.n),
                    wall=float(time.perf_counter() - self.t0))


# ===================================================================== VQ-A (categorical)
def a_energy(h, J, cfgs, npairs):
    """Vectorised table energy for (b,k) integer configurations."""
    cfgs = np.atleast_2d(np.asarray(cfgs, int))
    out = np.zeros(len(cfgs))
    for s in range(len(h)):
        out += h[s][cfgs[:, s]]
    for (s, u), Jm in J.items():
        out += Jm[cfgs[:, s], cfgs[:, u]]
    return out / npairs


def a_exhaustive(h, J, npairs):
    k, m = len(h), len(h[0])
    c = Counter()
    grids = np.meshgrid(*[np.arange(m)] * k, indexing="ij")
    cfgs = np.stack([g.ravel() for g in grids], 1)
    E = a_energy(h, J, cfgs, npairs)
    c.add(len(cfgs))
    b = int(np.argmin(E))
    r = c.done(cfgs[b], E[b])
    r["all_E"] = E
    r["all_cfg"] = cfgs
    return r


def a_partial(h, J, cfg, assigned, npairs):
    """Energy of the terms whose every endpoint is in `assigned` (an exact prefix)."""
    cfg = np.atleast_2d(np.asarray(cfg, int))
    aset = set(int(x) for x in assigned)
    out = np.zeros(len(cfg))
    for s in aset:
        out += h[s][cfg[:, s]]
    for (s, u), Jm in J.items():
        if s in aset and u in aset:
            out += Jm[cfg[:, s], cfg[:, u]]
    return out / npairs


def a_greedy(h, J, npairs, order=None):
    """Forward greedy: fix segments in `order`, each to the piece minimising the energy of
    the terms among the already-assigned segments (a true prefix objective)."""
    k, m = len(h), len(h[0])
    c = Counter()
    cfg = np.zeros(k, int)
    order = list(range(k)) if order is None else list(order)
    assigned = []
    for s in order:
        cand = np.tile(cfg, (m, 1))
        cand[:, s] = np.arange(m)
        E = a_partial(h, J, cand, assigned + [s], npairs)
        c.add(m)
        cfg[s] = int(np.argmin(E))
        assigned.append(s)
    return c.done(cfg, a_energy(h, J, cfg, npairs)[0])


def a_local_search(h, J, npairs, cfg0, c=None):
    """Steepest-descent 1-opt (change one segment) to convergence."""
    k, m = len(h), len(h[0])
    c = Counter() if c is None else c
    cfg = np.asarray(cfg0, int).copy()
    cur = a_energy(h, J, cfg, npairs)[0]
    while True:
        best, bcfg = cur, None
        for s in range(k):
            cand = np.tile(cfg, (m, 1))
            cand[:, s] = np.arange(m)
            E = a_energy(h, J, cand, npairs)
            c.add(m)
            b = int(np.argmin(E))
            if E[b] < best - 1e-12:
                best, bcfg = E[b], cand[b]
        if bcfg is None:
            break
        cfg, cur = bcfg, best
    return c.done(cfg, cur)


def a_greedy_ls(h, J, npairs):
    g = a_greedy(h, J, npairs)
    c = Counter()
    c.n = g["evals"]
    return a_local_search(h, J, npairs, g["cfg"], c)


def a_anneal(h, J, npairs, steps=3000, seed=0, T0=None, T1=None, restarts=1):
    k, m = len(h), len(h[0])
    c = Counter()
    rng = np.random.default_rng(seed)
    bestc, beste = None, np.inf
    for r in range(restarts):
        cfg = rng.integers(0, m, k)
        cur = a_energy(h, J, cfg, npairs)[0]
        c.add(1)
        if T0 is None:
            probe = a_energy(h, J, rng.integers(0, m, (64, k)), npairs)
            c.add(64)
            t0, t1 = float(probe.std()) + 1e-9, (float(probe.std()) + 1e-9) * 1e-3
        else:
            t0, t1 = T0, T1
        for it in range(steps):
            T = t0 * (t1 / t0) ** (it / max(steps - 1, 1))
            s = int(rng.integers(0, k))
            new = cfg.copy()
            new[s] = int(rng.integers(0, m))
            e = a_energy(h, J, new, npairs)[0]
            c.add(1)
            if e < cur or rng.random() < np.exp(-(e - cur) / max(T, 1e-12)):
                cfg, cur = new, e
            if cur < beste:
                beste, bestc = cur, cfg.copy()
    return c.done(bestc, beste)


def a_beam(h, J, npairs, width=8):
    """Beam search over segments in index order, keeping `width` partial assemblies.

    Partial energy uses only the blocks among the already-assigned segments, which is an
    exact prefix of the total (no bound needed, the terms are simply not yet present).
    """
    k, m = len(h), len(h[0])
    c = Counter()
    beam = [np.zeros(0, int)]
    for s in range(k):
        cand = []
        for pref in beam:
            for a in range(m):
                cand.append(np.concatenate([pref, [a]]))
        cand = np.stack(cand)
        E = np.zeros(len(cand))
        for u in range(s + 1):
            E += h[u][cand[:, u]]
        for (u, v), Jm in J.items():
            if u <= s and v <= s:
                E += Jm[cand[:, u], cand[:, v]]
        c.add(len(cand))
        keep = np.argsort(E)[:width]
        beam = [cand[i] for i in keep]
    E = a_energy(h, J, np.stack(beam), npairs)
    b = int(np.argmin(E))
    return c.done(beam[b], E[b])


def a_milp(h, J, npairs, time_limit=30.0):
    """Exact MILP via the Glover linearisation of the one-hot quadratic assignment.

    Binary x_{s,a} with sum_a x_{s,a} = 1 per segment, and y_{sa,tb} >= x_{sa} + x_{tb} - 1,
    y >= 0.  With NON-NEGATIVE quadratic coefficients (risk is a Bayes risk, so J >= 0)
    the minimisation drives y down to exactly x_{sa} x_{tb}, so the relaxation of y to a
    continuous variable is exact and no upper linking constraint is needed.
    """
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import lil_matrix
    k, m = len(h), len(h[0])
    c = Counter()
    xs = [(s, a) for s in range(k) for a in range(m)]
    xi = {p: i for i, p in enumerate(xs)}
    ys = [(s, a, u, b) for (s, u) in J for a in range(m) for b in range(m)]
    yi = {p: len(xs) + i for i, p in enumerate(ys)}
    nv = len(xs) + len(ys)
    cvec = np.zeros(nv)
    for (s, a) in xs:
        cvec[xi[(s, a)]] = h[s][a] / npairs
    for (s, a, u, b) in ys:
        cvec[yi[(s, a, u, b)]] = J[(s, u)][a, b] / npairs
    rows, lo, hi = [], [], []
    A = lil_matrix((k + len(ys), nv))
    r = 0
    for s in range(k):
        for a in range(m):
            A[r, xi[(s, a)]] = 1.0
        lo.append(1.0); hi.append(1.0); r += 1
    for (s, a, u, b) in ys:
        A[r, yi[(s, a, u, b)]] = 1.0
        A[r, xi[(s, a)]] = -1.0
        A[r, xi[(u, b)]] = -1.0
        lo.append(-1.0); hi.append(np.inf); r += 1
    integrality = np.zeros(nv)
    integrality[: len(xs)] = 1
    res = milp(c=cvec, constraints=LinearConstraint(A.tocsr(), lo, hi),
               integrality=integrality, bounds=Bounds(0, 1),
               options=dict(time_limit=time_limit))
    if res.x is None:
        return c.done(np.zeros(k, int), np.inf)
    X = res.x[: len(xs)].reshape(k, m)
    cfg = X.argmax(1)
    return c.done(cfg, a_energy(h, J, cfg, npairs)[0])


def a_random(h, J, npairs, n=None, seed=0):
    """Uniform random sampling at a matched evaluation budget -- the null the VQE must beat."""
    k, m = len(h), len(h[0])
    c = Counter()
    rng = np.random.default_rng(seed)
    n = int(n if n is not None else m ** k)
    cfgs = rng.integers(0, m, (n, k))
    E = a_energy(h, J, cfgs, npairs)
    c.add(n)
    b = int(np.argmin(E))
    return c.done(cfgs[b], E[b])


# ======================================================================== VQ-B (subsets)
def b_energy(M, b, const, X):
    """agg's H(x) = x'Mx/k^2 - 2 b'x/k + const, k = sum(x); +inf for k < 2."""
    X = np.atleast_2d(np.asarray(X, float))
    kk = X.sum(1)
    out = np.full(len(X), np.inf)
    ok = kk >= 2
    if ok.any():
        Xo = X[ok]; ko = kk[ok]
        out[ok] = (np.einsum("bi,ij,bj->b", Xo, M, Xo) / ko ** 2
                   - 2.0 * (Xo @ b) / ko + const)
    return out


def b_exhaustive(M, bv, const):
    N = len(bv)
    c = Counter()
    B = ((np.arange(1 << N)[:, None] >> np.arange(N - 1, -1, -1)[None, :]) & 1).astype(float)
    E = b_energy(M, bv, const, B)
    c.add(len(B))
    i = int(np.argmin(E))
    r = c.done(B[i], E[i])
    r["all_E"] = E
    return r


def b_greedy(M, bv, const):
    N = len(bv)
    c = Counter()
    cur = np.zeros(N)
    best, bx = np.inf, None
    for _ in range(N):
        off = np.where(cur == 0)[0]
        if not len(off):
            break
        C = np.tile(cur, (len(off), 1))
        C[np.arange(len(off)), off] = 1
        E = b_energy(M, bv, const, C)
        c.add(len(off))
        i = int(np.argmin(E))
        cur = C[i]
        if E[i] < best:
            best, bx = E[i], cur.copy()
    if bx is None:
        bx, best = cur, b_energy(M, bv, const, cur)[0]
    return c.done(bx, best)


def b_local_search(M, bv, const, x0, c=None):
    """1-flip steepest descent to convergence."""
    N = len(bv)
    c = Counter() if c is None else c
    x = np.asarray(x0, float).copy()
    cur = b_energy(M, bv, const, x)[0]
    while True:
        C = np.tile(x, (N, 1))
        C[np.arange(N), np.arange(N)] = 1 - C[np.arange(N), np.arange(N)]
        E = b_energy(M, bv, const, C)
        c.add(N)
        i = int(np.argmin(E))
        if E[i] < cur - 1e-12:
            x, cur = C[i], E[i]
        else:
            break
    return c.done(x, cur)


def b_greedy_ls(M, bv, const):
    g = b_greedy(M, bv, const)
    c = Counter()
    c.n = g["evals"]
    return b_local_search(M, bv, const, g["cfg"], c)


def b_anneal(M, bv, const, steps=3000, seed=0, restarts=1):
    N = len(bv)
    c = Counter()
    rng = np.random.default_rng(seed)
    beste, bx = np.inf, None
    for _ in range(restarts):
        x = (rng.random(N) < 0.5).astype(float)
        if x.sum() < 2:
            x[rng.choice(N, 2, replace=False)] = 1
        cur = b_energy(M, bv, const, x)[0]
        c.add(1)
        probe = b_energy(M, bv, const, (rng.random((64, N)) < 0.5).astype(float))
        c.add(64)
        pf = probe[np.isfinite(probe)]
        t0 = float(pf.std()) + 1e-9 if len(pf) else 1.0
        t1 = t0 * 1e-3
        for it in range(steps):
            T = t0 * (t1 / t0) ** (it / max(steps - 1, 1))
            y = x.copy()
            i = int(rng.integers(0, N))
            y[i] = 1 - y[i]
            e = b_energy(M, bv, const, y)[0]
            c.add(1)
            if np.isfinite(e) and (e < cur or rng.random() < np.exp(-(e - cur) / max(T, 1e-12))):
                x, cur = y, e
            if cur < beste:
                beste, bx = cur, x.copy()
    return c.done(bx, beste)


def b_random(M, bv, const, n, seed=0):
    N = len(bv)
    c = Counter()
    rng = np.random.default_rng(seed)
    X = (rng.random((n, N)) < 0.5).astype(float)
    E = b_energy(M, bv, const, X)
    c.add(n)
    i = int(np.nanargmin(np.where(np.isfinite(E), E, np.inf)))
    return c.done(X[i], E[i])
