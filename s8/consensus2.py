"""How far does filter-then-consensus go, and what stops it?

The one thing that worked
S8-8 (`s8/inband.py`) established the only intervention in three sprints with a confidence
interval excluding zero: filter the shipped BLOSUM top-500 to the shipped score's own top
75, then return the CONSENSUS MEDOID of that subset.

    3.282 A vs 3.454 baseline,  -0.172 [-0.316, -0.027],  74W/47L,
    chosen leave-fold-out over 20 arms, negative in all five folds,
    bootstrap P(gain) 0.991, sign test p 0.0177, frac <2 A 0.214 -> 0.278.

The principle: the shipped score is a good coarse FILTER and a bad fine RANKER; consensus
is the reverse.  The score's own ordering INSIDE its top 75 is worth 0.10 A; consensus over
that identical set is worth 0.27 A.

This module does not re-derive that.  It pushes it until it stops improving and reports
where the ceiling is:

  1. THE CEILING.  For every filtered set, the ORACLE best inside it alongside what
     consensus returns -- the headroom consensus leaves on the table.
  2. THE JOINT SWEEP.  8 filter sizes x 7 filter keys x 16 consensus operators, because
     the optimum is joint and m was previously chosen over a narrow arm set.
  3. BETTER OPERATORS.  Weighted / trimmed / iterative / cluster-restricted medoids, the
     coordinate average, and the TORSION-SPACE circular mean (which needs no rebuild-free
     construction: bond geometry is exact by construction).  Circular means are computed
     as atan2 of summed unit vectors -- the naive mean of 170 and -170 degrees is 0, and
     that bug has appeared in this repo before.
  4. MULTI-HYPOTHESIS.  S8-2 measured ~2-3 populated basins in the near-native band, so a
     single medoid may sit between them.  Consensus per cluster, best-of-k as a labelled
     DIAGNOSTIC ceiling, and what a native-free rule picks among them.
  5. DOES REFINEMENT STACK.  A bounded Legacy relaxation of the consensus output.

INSTRUMENT
The same 126-target tuning instrument, the same K=500 BLOSUM pools, rebuilt through
`s8.generate`'s cached universes.  `stage_build` ASSERTS the baseline reproduces the
shipped score's 3.454 A and pool best 1.711 A before anything else runs.

DEPLOYABLE vs ORACLE -- the split is structural
Every filter key and every consensus operator is a function of the candidate coordinates,
their torsions, the target sequence and the fold's own out-of-fold distogram.  `rr` (each
candidate's CA-RMSD to the native) and `nat_ca` enter ONLY as reporting labels and in rows
explicitly named `oracle`.  `stage_leak` NaN-poisons both and asserts every filter order
and every operator output is bit-identical.

TIES
`sel_of` is imported from `s8.inband` rather than rewritten.  `np.argmin` returns the FIRST
tied index and the candidate order is informative (BLOSUM rank), so a signal with many ties
reads out the sort order and can fake a spectacular win; `sel_of` averages the true RMSD
over the whole tied argmin set, which is the expectation under random tie-breaking.

    python -m s8.consensus2 build     # per-target cache: pools, scores, pairwise RMSD
    python -m s8.consensus2 sweep     # the joint filter x size x operator sweep
    python -m s8.consensus2 multi     # multi-hypothesis output, best-of-k, native-free rules
    python -m s8.consensus2 refine    # does a bounded Legacy relaxation stack on top?
    python -m s8.consensus2 dev       # ONE pre-registered dev-24 pass
    python -m s8.consensus2 leak      # NaN-poison assertion for every deployable arm
    python -m s8.consensus2 report    # print everything already computed
"""
import json
import math
import os
import sys
import time
import zlib

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np                                                          # noqa: E402

os.environ.setdefault("NT", "2")
try:
    import torch
    torch.set_num_threads(2)
except Exception:                                                           # pragma: no cover
    torch = None

import protein_geometry as geo                                              # noqa: E402
from s7 import audit, debias, poolsize                                      # noqa: E402
from s8 import generate as gen                                              # noqa: E402
from s8 import inband                                                       # noqa: E402
from s8.inband import sel_of, _paired                                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "consensus2_cache")
SWEEP = os.path.join(CACHE, "sweep")
DCACHE = os.path.join(HERE, "inband_devpool")      # read-only: built by s8/inband.py

K = 500
SEED = 20260907
CA_CA = 3.80          #: the virtual CA-CA bond length ideal geometry builds with

#: Filter sizes.  The mandate's dense sweep; 500 = no filter (the whole pool).
MS = (10, 25, 50, 75, 100, 150, 250, 500)

#: The pre-registered incumbent, for reference rows.
INCUMBENT = ("sc", 75, "medoid")

#: THE single arm carried to dev-24, pinned here BEFORE `stage_dev` was ever run and
#: asserted by `test_consensus2.py`.
#:
#: The leave-fold-out sweep over all 1120 arms picks the SYNTHESIS family in every fold and
#: the cell `sc|75` in three of five (`sc|100`/`sc|150` in the other two, within 0.008 A).
#: The arm taken forward is that cell's synthesis made PHYSICALLY VALID -- `fit`, the
#: torsion projection -- rather than the raw coordinate average, because the average's mean
#: CA-CA bond is 2.961 A against a native 3.812 and it is therefore not a peptide.  That
#: choice costs 0.156 A of tuning-set RMSD (3.048 -> 3.204) and is not negotiable.
DEV_ARM = ("sc", 75, "fit")


# ============================================================ small numerics
def superpose_batch(P, ref):
    """Batched Kabsch: rotate every ``(n,3)`` in ``P`` onto ``ref``. Reflections forbidden.

    Numerically identical to `protein_geometry.kabsch_superpose` applied one at a time --
    `t_superpose_matches_geo` holds it to 1e-9 -- but written batched because the
    coordinate-average operators superpose up to 500 structures per cell per target.
    """
    P = np.asarray(P, float)
    ref = np.asarray(ref, float)
    if P.ndim == 2:
        P = P[None]
    Pm = P.mean(1, keepdims=True)
    Rm = ref.mean(0, keepdims=True)
    Pc = P - Pm
    Rc = ref - Rm
    H = np.einsum("bni,nj->bij", Pc, Rc)
    U, S, Vt = np.linalg.svd(H)
    det = np.linalg.det(np.matmul(np.transpose(Vt, (0, 2, 1)),
                                  np.transpose(U, (0, 2, 1))))
    Dm = np.zeros((len(P), 3, 3))
    Dm[:, 0, 0] = 1.0
    Dm[:, 1, 1] = 1.0
    Dm[:, 2, 2] = np.sign(det)
    R = np.matmul(np.matmul(np.transpose(Vt, (0, 2, 1)), Dm),
                  np.transpose(U, (0, 2, 1)))
    return np.einsum("bij,bnj->bni", R, Pc) + Rm


def circmean(A, w=None, axis=0):
    """Circular mean of angles in radians: atan2 of the summed unit vectors.

    The naive arithmetic mean of +170 and -170 degrees is 0 degrees, which is the exact
    opposite of the right answer (180).  This bug has appeared in this repository before,
    so `t_circmean_wraparound` pins the wrap case explicitly.
    """
    A = np.asarray(A, float)
    s = np.average(np.sin(A), axis=axis, weights=w)
    c = np.average(np.cos(A), axis=axis, weights=w)
    return np.arctan2(s, c)


def _rg(X):
    """Radius of gyration of a CA trace."""
    X = np.asarray(X, float)
    return float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))


def _rank(v):
    from scipy.stats import rankdata
    return rankdata(np.asarray(v, float))


# ============================================================ stage: build
def _cpath(pdbid):
    return os.path.join(CACHE, f"{pdbid}.npz")


def build_one(pdbid):
    """Everything one target's sweep needs.  `rr`/`nat_ca` are REPORTING LABELS only."""
    u = gen.load_univ(pdbid)
    pred = poolsize.load_pred(pdbid)
    n = int(u["n"])
    idx = np.asarray(u["order"], int)[:K]
    W = np.asarray(u["W"], float)[idx]
    PHI = np.asarray(u["PHI"], float)[idx]
    PSI = np.asarray(u["PSI"], float)[idx]
    D = gen.pair_D(W, n)
    sc = gen.shipped_score(pred, D)
    typ = np.abs(D - D.mean(0)[None, :]).mean(1)
    B = len(W)
    P = np.zeros((B, B), np.float32)
    for a in range(B):
        P[a] = audit.kabsch_rmsd_batch(W, W[a])
    # the ideal-geometry rebuild, so the rebuild displacement is measured not assumed
    RB = geo.build_backbone_batch(PHI, PSI)["CA"]
    reb = np.array([geo.rmsd(geo.kabsch_superpose(RB[b], W[b]), W[b]) for b in range(B)])
    nat = np.asarray(u["nat_ca"], float)
    return {"pdb": pdbid, "n": n, "fold": int(u["fold"]), "seq": str(u["seq"]),
            "W": W.astype(np.float32), "PHI": PHI.astype(np.float32),
            "PSI": PSI.astype(np.float32), "D": D.astype(np.float32),
            "sc": sc.astype(np.float32), "typ": typ.astype(np.float32),
            "P": P, "reb": reb.astype(np.float32),
            "rr_reb": audit.kabsch_rmsd_batch(RB, nat).astype(np.float32),
            "rr": np.asarray(u["rr"], float)[idx].astype(np.float32),
            "nat_ca": nat.astype(np.float32)}


def stage_build(verbose=True):
    os.makedirs(CACHE, exist_ok=True)
    tgs = gen.cached_targets()
    todo = [p for p in tgs if not os.path.exists(_cpath(p))]
    if verbose:
        print(f"{len(todo)}/{len(tgs)} targets to build", flush=True)
    for k, pdbid in enumerate(todo):
        t0 = time.time()
        blob = build_one(pdbid)
        tmp = _cpath(pdbid) + ".tmp.npz"
        np.savez_compressed(tmp, **blob)
        os.replace(tmp, _cpath(pdbid))
        if verbose:
            print(f"[{k+1:3d}/{len(todo)}] {pdbid} n={blob['n']} "
                  f"({time.time()-t0:.1f}s)", flush=True)
    # instrument validity: the shipped baseline must reproduce sprint 7/8 exactly
    v = validity()
    if verbose:
        print(f"VALIDITY  shipped score {v['score']:.3f} (expect 3.454)  "
              f"pool best {v['pool_best']:.3f} (expect 1.711)  "
              f"medoid75 {v['medoid75']:.3f} (expect 3.282)", flush=True)
    assert abs(v["score"] - 3.454) < 5e-3, v
    assert abs(v["pool_best"] - 1.711) < 5e-3, v
    assert abs(v["medoid75"] - 3.282) < 5e-3, v
    return v


def load_cache(pdbid):
    z = np.load(_cpath(pdbid), allow_pickle=True)
    out = {"pdb": str(z["pdb"]), "n": int(z["n"]), "fold": int(z["fold"]),
           "seq": str(z["seq"])}
    for k in ("W", "PHI", "PSI", "D", "sc", "typ", "P", "reb", "rr_reb", "rr", "nat_ca"):
        out[k] = z[k].astype(float)
    return out


def cached():
    if not os.path.isdir(CACHE):
        return []
    return sorted(f[:-4] for f in os.listdir(CACHE) if f.endswith(".npz"))


def validity():
    """Reproduce the three numbers this whole study stands on."""
    sc, pb, md = [], [], []
    for pdbid in cached():
        c = load_cache(pdbid)
        rr = c["rr"]
        sc.append(sel_of(c["sc"], rr))
        pb.append(float(rr.min()))
        sub = np.argsort(c["sc"], kind="stable")[:75]
        Ps = c["P"][np.ix_(sub, sub)]
        md.append(sel_of(Ps.sum(1) / max(len(sub) - 1, 1), rr[sub]))
    return {"n": len(sc), "score": float(np.mean(sc)),
            "pool_best": float(np.mean(pb)), "medoid75": float(np.mean(md))}


# ============================================================ DEPLOYABLE: filters
#: Filter keys.  Each returns an ORDER over the pool -- lowest key first -- and the filter
#: keeps the first `m`.  Nothing here reads a native.
FILTERS = ("sc", "ty", "scty", "scty2", "ref1", "refmed", "rand")


def filter_orders(c):
    """{name: (K,) int order} for every filter key. DEPLOYABLE."""
    sc, typ, P = c["sc"], c["typ"], c["P"]
    out = {}
    out["sc"] = np.argsort(sc, kind="stable")
    out["ty"] = np.argsort(typ, kind="stable")
    rs, rt = _rank(sc), _rank(typ)
    out["scty"] = np.argsort(rs + rt, kind="stable")
    out["scty2"] = np.argsort(rs + 0.5 * rt, kind="stable")
    # S8-3's form: predict a STRUCTURE, rank by CA-RMSD to it.  Two native-free
    # references: the score's own top-1, and the medoid of the score's top-75.
    i0 = int(out["sc"][0])
    out["ref1"] = np.argsort(P[:, i0], kind="stable")
    s75 = out["sc"][:75]
    md = P[np.ix_(s75, s75)].sum(1)
    j0 = int(s75[int(np.argmin(md))])
    out["refmed"] = np.argsort(P[:, j0], kind="stable")
    rng = np.random.default_rng(SEED + zlib.crc32(c["pdb"].encode()))
    out["rand"] = rng.permutation(len(sc))
    return out


# ============================================================ DEPLOYABLE: operators
#: Consensus operators.  `SEL_OPS` return a score vector over the subset (lower = chosen,
#: a REAL pool candidate); `CON_OPS` return a CONSTRUCTED (n,3) trace whose RMSD to the
#: native is measured directly.
SEL_OPS = ("medoid", "medoid_sq", "typic", "wmed_ty", "wmed_sc", "trim25", "trim50",
           "iter2", "clmed15", "clmed3", "snapavg", "score")
CON_OPS = ("avg", "avg_trim", "avg_iter", "avg_scale", "avg_wty",
           "tors", "tors_trim", "reb_med")
OPS = SEL_OPS + CON_OPS


def _medoid(Ps):
    """Mean CA-RMSD of each member to every OTHER member -- the incumbent's statistic."""
    return Ps.sum(1) / max(Ps.shape[0] - 1, 1)


def _trim_keep(Ps, frac):
    """Indices (into the subset) of the `frac` fraction most central members."""
    md = _medoid(Ps)
    k = max(2, int(round(frac * len(md))))
    return np.argsort(md, kind="stable")[:k]


def _largest_cluster(Ps, thr=None, kmax=None):
    """Members of the largest average-linkage cluster. `thr` in A, or `kmax` clusters."""
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    B = Ps.shape[0]
    if B < 3:
        return np.arange(B)
    d = squareform(np.maximum((Ps + Ps.T) / 2.0, 0.0), checks=False)
    Z = linkage(d, method="average")
    lab = (fcluster(Z, t=thr, criterion="distance") if thr is not None
           else fcluster(Z, t=min(kmax, B), criterion="maxclust"))
    vals, cnt = np.unique(lab, return_counts=True)
    return np.where(lab == vals[int(np.argmax(cnt))])[0]


def op_value(op, c, sub):
    """Run one consensus operator on one filtered subset. DEPLOYABLE.

    Returns ``("sel", v)`` -- a score vector over `sub`, lower is chosen -- or
    ``("con", C)`` -- a constructed (n, 3) CA trace.
    """
    P, W, D = c["P"], c["W"], c["D"]
    PHI, PSI = c["PHI"], c["PSI"]
    Ps = P[np.ix_(sub, sub)]
    B = len(sub)
    md = _medoid(Ps)

    if op == "medoid":
        return "sel", md
    if op == "medoid_sq":
        return "sel", (Ps ** 2).sum(1) / max(B - 1, 1)
    if op == "typic":
        Ds = D[sub]
        return "sel", np.abs(Ds - Ds.mean(0)[None, :]).mean(1)
    if op == "score":
        return "sel", c["sc"][sub]
    if op in ("wmed_ty", "wmed_sc"):
        v = c["typ"][sub] if op == "wmed_ty" else c["sc"][sub]
        z = (v - v.mean()) / max(v.std(), 1e-12)
        w = np.exp(-z)
        w = w / w.sum()
        return "sel", Ps @ w
    if op in ("trim25", "trim50"):
        keep = _trim_keep(Ps, 0.75 if op == "trim25" else 0.50)
        v = np.full(B, np.inf)
        v[keep] = Ps[np.ix_(keep, keep)].sum(1) / max(len(keep) - 1, 1)
        return "sel", v
    if op == "iter2":
        cur = np.arange(B)
        for _ in range(2):
            m0 = cur[int(np.argmin(Ps[np.ix_(cur, cur)].sum(1)))]
            r = Ps[m0, cur]
            keep = cur[r <= np.median(r)]
            if len(keep) < 3:
                break
            cur = keep
        v = np.full(B, np.inf)
        v[cur] = Ps[np.ix_(cur, cur)].sum(1) / max(len(cur) - 1, 1)
        return "sel", v
    if op in ("clmed15", "clmed3"):
        cl = (_largest_cluster(Ps, thr=1.5) if op == "clmed15"
              else _largest_cluster(Ps, kmax=3))
        v = np.full(B, np.inf)
        v[cl] = Ps[np.ix_(cl, cl)].sum(1) / max(len(cl) - 1, 1)
        return "sel", v
    if op == "snapavg":
        C = _avg_struct(W[sub], Ps)
        return "sel", audit.kabsch_rmsd_batch(W[sub], C)
    if op == "avg":
        return "con", _avg_struct(W[sub], Ps)
    if op == "avg_trim":
        keep = _trim_keep(Ps, 0.75)
        return "con", _avg_struct(W[sub][keep], Ps[np.ix_(keep, keep)])
    if op == "avg_iter":
        # The Frechet / Karcher mean: re-superpose every member onto the RUNNING average
        # and re-average until it stops moving.  One-shot `avg` superposes onto the medoid
        # once, which is the first iterate of exactly this.
        C = _avg_struct(W[sub], Ps)
        for _ in range(5):
            C2 = superpose_batch(W[sub], C).mean(0)
            moved = float(np.abs(C2 - C).max())
            C = C2
            if moved < 1e-9:
                break
        return "con", C
    if op == "avg_scale":
        # The shrinkage control.  Averaging superposed structures CONTRACTS the chain --
        # the mean of points scattered about a curve lies inside it -- so the average is
        # not a physical backbone and a shrunken object can flatter CA-RMSD.  Rescale it
        # about its centroid to restore the 3.80 A CA-CA step.  If the gain survives,
        # contraction is not what buys it.
        C = _avg_struct(W[sub], Ps)
        st = np.linalg.norm(C[1:] - C[:-1], axis=1).mean()
        return "con", (C - C.mean(0)) * (CA_CA / max(st, 1e-9)) + C.mean(0)
    if op == "avg_wty":
        v = c["typ"][sub]
        z = (v - v.mean()) / max(v.std(), 1e-12)
        w = np.exp(-z)
        b = int(np.argmin(md))
        return "con", np.average(superpose_batch(W[sub], W[sub][b]), axis=0, weights=w)
    if op == "tors":
        return "con", _tors_struct(PHI[sub], PSI[sub])
    if op == "tors_trim":
        keep = _trim_keep(Ps, 0.75)
        return "con", _tors_struct(PHI[sub][keep], PSI[sub][keep])
    if op == "reb_med":
        # CONTROL that isolates the rebuild cost: take the medoid, throw away its real
        # coordinates and rebuild it from its own torsions.  Any torsion-space operator
        # pays this, so it is the right comparator for `tors`.
        b = int(np.argmin(md))
        return "con", _tors_struct(PHI[sub][b][None], PSI[sub][b][None])
    raise KeyError(op)


def _avg_struct(Wsub, Ps):
    """Coordinate average after superposing every member onto the subset's medoid."""
    b = int(np.argmin(_medoid(Ps)))
    return superpose_batch(Wsub, Wsub[b]).mean(0)


def _tors_struct(PHI, PSI):
    """CIRCULAR mean of the torsions, rebuilt with ideal bond geometry.

    Bond lengths and angles are exact by construction, so this construction has no
    'is it a physical chain' defect -- what it pays instead is the rebuild displacement,
    which `reb_med` measures on its own.
    """
    ph = circmean(np.asarray(PHI, float), axis=0)
    ps = circmean(np.asarray(PSI, float), axis=0)
    return geo.build_backbone_batch(ph[None], ps[None])["CA"][0]


# ============================================================ stage: sweep
def _spath(pdbid):
    return os.path.join(SWEEP, f"{pdbid}.json")


def sweep_one(c):
    """Every (filter, m, operator) cell for one target, plus the oracle ceilings."""
    rr, nat = c["rr"], c["nat_ca"]
    orders = filter_orders(c)
    cells = {}
    for f in FILTERS:
        o = orders[f]
        for m in MS:
            sub = o[:m]
            rs = rr[sub]
            cell = {"oracle": float(rs.min()), "submean": float(rs.mean()),
                    "ops": {}, "disp": {}, "geom": {}}
            for op in OPS:
                kind, v = op_value(op, c, sub)
                if kind == "sel":
                    cell["ops"][op] = sel_of(v, rs)
                else:
                    cell["ops"][op] = float(
                        audit.kabsch_rmsd_batch(v[None], nat)[0])
                    # how far the CONSTRUCTED structure sits from any real candidate,
                    # and whether it is still a physical chain at all
                    cell["disp"][op] = float(
                        audit.kabsch_rmsd_batch(c["W"][sub], v).min())
                    st = np.linalg.norm(v[1:] - v[:-1], axis=1)
                    cell["geom"][op] = [float(st.mean()), float(st.std()),
                                        float(_rg(v))]
            cells[f"{f}|{m}"] = cell
    s75 = orders["sc"][:75]
    natst = np.linalg.norm(nat[1:] - nat[:-1], axis=1)
    return {"pdb": c["pdb"], "fold": c["fold"], "n": c["n"],
            "nat_rg": _rg(nat), "nat_step": float(natst.mean()),
            "med_rg": _rg(c["W"][int(s75[int(np.argmin(
                _medoid(c["P"][np.ix_(s75, s75)])))])]),
            "pool_best": float(rr.min()), "pool_mean": float(rr.mean()),
            "score": sel_of(c["sc"], rr),
            "reb": float(c["reb"].mean()),
            "rr_reb_gap": float((c["rr_reb"] - rr).mean()),
            "cells": cells}


def stage_sweep(verbose=True):
    os.makedirs(SWEEP, exist_ok=True)
    tgs = cached()
    per = []
    for k, pdbid in enumerate(tgs):
        p = _spath(pdbid)
        if os.path.exists(p):
            with open(p) as fh:
                per.append(json.load(fh))
            continue
        t0 = time.time()
        row = sweep_one(load_cache(pdbid))
        tmp = p + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(row, fh, default=float)
        os.replace(tmp, p)
        per.append(row)
        if verbose:
            print(f"[{k+1:3d}/{len(tgs)}] {pdbid} ({time.time()-t0:.1f}s)", flush=True)
    return _analyse_sweep(per)


def _arm_vec(per, f, m, op):
    return np.array([r["cells"][f"{f}|{m}"]["ops"][op] for r in per], float)


def _analyse_sweep(per):
    base = np.array([r["score"] for r in per], float)
    inc = _arm_vec(per, *INCUMBENT)
    folds = np.array([r["fold"] for r in per], int)
    out = {"n_targets": len(per), "K": K, "MS": list(MS), "FILTERS": list(FILTERS),
           "OPS": list(OPS), "incumbent": list(INCUMBENT),
           "base_score": float(base.mean()),
           "base_se": float(base.std(ddof=1) / math.sqrt(len(base))),
           "incumbent_sel": float(inc.mean()),
           "pool_best": float(np.mean([r["pool_best"] for r in per])),
           "pool_mean": float(np.mean([r["pool_mean"] for r in per])),
           "rebuild": {
               "displacement": float(np.mean([r["reb"] for r in per])),
               "rr_penalty": float(np.mean([r["rr_reb_gap"] for r in per]))},
           "nat_rg": float(np.mean([r["nat_rg"] for r in per])),
           "nat_step": float(np.mean([r["nat_step"] for r in per])),
           "med_rg_ratio": float(np.mean([r["med_rg"] / r["nat_rg"] for r in per])),
           "cells": {}, "arms": {}, "lfo": {}, "per_target": per}

    def summarise(v):
        return {"sel": float(v.mean()),
                "se": float(v.std(ddof=1) / math.sqrt(len(v))),
                "median": float(np.median(v)), "sd": float(v.std(ddof=1)),
                "worst": float(v.max()), "best": float(v.min()),
                "f2": float((v < 2.0).mean()), "f15": float((v < 1.5).mean()),
                "vs_score": _paired(v, base), "vs_incumbent": _paired(v, inc)}

    # the ceiling table: oracle-in-filter vs what consensus returns
    for f in FILTERS:
        for m in MS:
            key = f"{f}|{m}"
            orc = np.array([r["cells"][key]["oracle"] for r in per], float)
            sm = np.array([r["cells"][key]["submean"] for r in per], float)
            best_op = min(OPS, key=lambda o: _arm_vec(per, f, m, o).mean())
            bv = _arm_vec(per, f, m, best_op)
            out["cells"][key] = {
                "oracle_in_filter": float(orc.mean()),
                "submean": float(sm.mean()),
                "medoid": float(_arm_vec(per, f, m, "medoid").mean()),
                "best_op": best_op, "best_op_sel": float(bv.mean()),
                "headroom_medoid": float((_arm_vec(per, f, m, "medoid") - orc).mean()),
                "headroom_best": float((bv - orc).mean())}
    for f in FILTERS:
        for m in MS:
            for op in OPS:
                v = _arm_vec(per, f, m, op)
                d = summarise(v)
                dsp = [r["cells"][f"{f}|{m}"]["disp"].get(op) for r in per]
                if dsp[0] is not None:
                    d["disp"] = float(np.mean(dsp))
                    G = np.array([r["cells"][f"{f}|{m}"]["geom"][op] for r in per], float)
                    d["step"] = float(G[:, 0].mean())
                    d["step_sd"] = float(G[:, 1].mean())
                    d["rg_ratio"] = float(np.mean(
                        G[:, 2] / np.array([r["nat_rg"] for r in per], float)))
                out["arms"][f"{f}|{m}:{op}"] = d

    # ------------------------------------------------- leave-fold-out arm selection
    # No target contributes to choosing the arm it is then scored under.  `v` MUST stay
    # aligned with `per` (and therefore with `base`): mis-alignment here inflated a CI by
    # 2x in S8-8 on an identical mean.
    def lfo(cands, name):
        v = np.full(len(per), np.nan)
        picked = []
        for fd in range(5):
            tr = folds != fd
            if not tr.any() or not (folds == fd).any():
                continue
            best = min(cands, key=lambda a: float(
                _arm_vec(per, *a)[tr].mean()))
            picked.append(f"fold{fd}:{best[0]}|{best[1]}:{best[2]}")
            hv = _arm_vec(per, *best)
            v[folds == fd] = hv[folds == fd]
        assert np.isfinite(v).all(), "leave-fold-out coverage is incomplete"
        d = summarise(v)
        d["picked_per_fold"] = picked
        out["lfo"][name] = d
        return d

    allc = [(f, m, op) for f in FILTERS for m in MS for op in OPS]
    lfo(allc, "all")
    lfo([("sc", m, "medoid") for m in MS], "sc_medoid_m")
    lfo([(f, m, "medoid") for f in FILTERS for m in MS], "medoid_anyfilter")
    lfo([("sc", m, op) for m in MS for op in SEL_OPS], "sc_selops")
    lfo([("sc", m, op) for m in MS for op in OPS], "sc_anyop")
    _write("consensus2_sweep.json", out)
    _print_sweep(out)
    return out


def _write(name, obj):
    p = os.path.join(HERE, name)
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
    os.replace(tmp, p)


def _read(name):
    with open(os.path.join(HERE, name)) as fh:
        return json.load(fh)


def _print_sweep(out):
    print(f"\n=== JOINT FILTER x SIZE x OPERATOR SWEEP — {out['n_targets']} targets, "
          f"K={out['K']} ===")
    print(f"  pool best {out['pool_best']:.3f}  pool mean {out['pool_mean']:.3f}  "
          f"SHIPPED score {out['base_score']:.3f} +- {out['base_se']:.3f}  "
          f"incumbent(sc75:medoid) {out['incumbent_sel']:.3f}")
    print(f"  rebuild displacement {out['rebuild']['displacement']:.3f} A; "
          f"rebuilding a candidate costs {out['rebuild']['rr_penalty']:+.3f} A of "
          f"CA-RMSD to the native")

    print("\n--- CEILING: oracle best inside the filtered set vs what consensus returns")
    print(f"{'cell':12} {'oracle':>7} {'submean':>8} {'medoid':>7} {'headroom':>9} "
          f"{'best op':>10} {'best sel':>9}")
    for f in out["FILTERS"]:
        for m in out["MS"]:
            d = out["cells"][f"{f}|{m}"]
            print(f"{f+chr(124)+str(m):12} {d['oracle_in_filter']:7.3f} {d['submean']:8.3f} "
                  f"{d['medoid']:7.3f} {d['headroom_medoid']:9.3f} "
                  f"{d['best_op']:>10} {d['best_op_sel']:9.3f}")

    print("\n--- TOP 25 ARMS by mean selected CA-RMSD (descriptive; the LFO row is the "
          "honest number)")
    _arm_table(out, sorted(out["arms"].items(), key=lambda kv: kv[1]["sel"])[:25])

    print("\n--- OPERATOR at the incumbent filter (sc, m=75)")
    _arm_table(out, [(f"sc|75:{op}", out["arms"][f"sc|75:{op}"]) for op in out["OPS"]])

    print("\n--- SIZE sweep, sc filter, medoid operator")
    _arm_table(out, [(f"sc|{m}:medoid", out["arms"][f"sc|{m}:medoid"]) for m in out["MS"]])

    print("\n--- IS THE CONSTRUCTED AVERAGE STILL A CHAIN?  (native CA-CA "
          f"{out['nat_step']:.3f} A; the medoid, a real window, has rg/rg_nat "
          f"{out['med_rg_ratio']:.3f})")
    print(f"{'arm':22} {'sel':>7} {'CA-CA':>7} {'sd':>6} {'rg/rg_nat':>10} "
          f"{'to nearest real':>16}")
    for op in out["OPS"]:
        d = out["arms"].get(f"sc|75:{op}", {})
        if "step" not in d:
            continue
        print(f"{'sc|75:' + op:22} {d['sel']:7.3f} {d['step']:7.3f} "
              f"{d['step_sd']:6.3f} {d['rg_ratio']:10.3f} {d['disp']:16.3f}")

    print("\n--- LEAVE-FOLD-OUT arm selection (the honest numbers)")
    _arm_table(out, list(out["lfo"].items()))
    for k, d in out["lfo"].items():
        print(f"   {k:18} picks {d['picked_per_fold']}")


def _arm_table(out, rows):
    print(f"{'arm':22} {'sel':>7} {'+-se':>6} {'med':>7} {'<2A':>6} {'<1.5A':>6} "
          f"{'d vs score':>11} {'95% CI':>19} {'W/L':>9} {'d vs inc':>9}")
    for a, d in rows:
        p = d["vs_score"]
        q = d["vs_incumbent"]
        ci = f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}]"
        print(f"{a:22} {d['sel']:7.3f} {d['se']:6.3f} {d['median']:7.3f} "
              f"{d['f2']:6.3f} {d['f15']:6.3f} {p['mean_diff']:+11.3f} {ci:>19} "
              f"{p['n_better']:3d}/{p['n_worse']:<4d} {q['mean_diff']:+9.3f}")


# ============================================================ stage: fit
#: Cells the torsion-projection study runs on.  Kept small: it is an optimisation per
#: target per cell, not a table lookup.
FIT_CELLS = (("sc", 75), ("ref1", 100))
#: The cell the multistart (torsion-free) form is measured on, because it costs 4x.
FITMS_CELL = ("sc", 75)
#: Generic starting conformations for the torsion projection, in degrees: extended,
#: alpha-helix, beta-strand, polyproline II.  Using these instead of the medoid's own
#: (phi, psi) makes the projection runnable on ANY pool of bare CA traces -- the dev pools
#: store coordinates and no torsions -- and on 24 tuning targets the two agree to
#: -0.010 A of selected CA-RMSD, so nothing is lost by not having the torsions.
FIT_STARTS = ((-120.0, 130.0), (-57.0, -47.0), (-139.0, 135.0), (-75.0, 145.0))


def fit_torsions(C, phi0, psi0, maxiter=300):
    """The ideal-geometry chain CLOSEST to a constructed consensus. DEPLOYABLE.

    The coordinate average is not a physical backbone -- averaging superposed structures
    contracts the chain, because the mean of points scattered about a curve lies inside
    it.  Rather than rescale it (`avg_scale`, a crude global fix), project it onto the
    manifold of ideal-geometry chains: minimise CA-RMSD(build(phi, psi), C) over the
    torsions, started from the medoid's own.  The result has EXACT bond geometry and is
    the nearest thing to the consensus that a real peptide could be.

    Reads only `C` and a starting torsion pair, both native-free.
    """
    from scipy.optimize import minimize
    C = np.asarray(C, float)
    n = len(C)
    eps = 1e-5
    E = np.eye(2 * n) * eps

    def build(X):
        return geo.build_backbone_batch(X[:, :n], X[:, n:])["CA"]

    def fg(x):
        # The finite-difference gradient is the whole cost, so all 2n perturbations are
        # built in ONE batched call.  Handed to scipy as `jac` it is 27x fewer builds per
        # iteration than letting L-BFGS-B difference the scalar objective itself (39 s ->
        # ~1 s per target), and the numbers agree.
        X = np.vstack([x[None], x[None] + E])
        r = audit.kabsch_rmsd_batch(build(X), C)
        return float(r[0]), (r[1:] - r[0]) / eps

    x0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
    f0 = fg(x0)[0]
    r = minimize(fg, x0, jac=True, method="L-BFGS-B", options={"maxiter": maxiter})
    return build(r.x[None])[0], r.x[:n].copy(), r.x[n:].copy(), float(r.fun), float(f0)


def fit_consensus(C, extra=None, maxiter=300):
    """`fit_torsions` from several starts, keeping whichever lands nearest `C`.

    DEPLOYABLE and torsion-free: the starts are fixed generic conformations, so this needs
    only the constructed consensus trace.  `extra` optionally adds one more start (the
    medoid's own torsions, when the pool carries them).
    """
    n = len(C)
    best = None
    starts = [(np.full(n, math.radians(a)), np.full(n, math.radians(b)))
              for a, b in FIT_STARTS]
    if extra is not None:
        starts.append((np.asarray(extra[0], float), np.asarray(extra[1], float)))
    for ph, ps in starts:
        got = fit_torsions(C, ph, ps, maxiter=maxiter)
        if best is None or got[3] < best[3]:
            best = got
    return best


def stage_fit(cells=FIT_CELLS, verbose=True):
    """Does projecting the coordinate consensus back onto a real chain keep the gain?"""
    tgs = cached()
    per = []
    for k, pdbid in enumerate(tgs):
        c = load_cache(pdbid)
        orders = filter_orders(c)
        nat = c["nat_ca"]
        row = {"pdb": pdbid, "fold": c["fold"], "score": sel_of(c["sc"], c["rr"]),
               "cells": {}}
        for f, m in cells:
            sub = orders[f][:m]
            Ps = c["P"][np.ix_(sub, sub)]
            b = int(np.argmin(_medoid(Ps)))
            gi = int(sub[b])
            C = _avg_struct(c["W"][sub], Ps)
            F, _fp, _fs, d1, d0 = fit_torsions(C, c["PHI"][gi], c["PSI"][gi])
            st = np.linalg.norm(F[1:] - F[:-1], axis=1)
            cell = {"medoid": float(c["rr"][gi]),
                    "avg": float(audit.kabsch_rmsd_batch(C[None], nat)[0]),
                    "fit": float(audit.kabsch_rmsd_batch(F[None], nat)[0]),
                    "fit_to_avg": d1, "start_to_avg": d0,
                    "step": float(st.mean()), "step_sd": float(st.std())}
            if (f, m) == FITMS_CELL:
                G, _gp, _gs, d2, _ = fit_consensus(C)
                cell["fitms"] = float(audit.kabsch_rmsd_batch(G[None], nat)[0])
                cell["fitms_to_avg"] = d2
            row["cells"][f"{f}|{m}"] = cell
        per.append(row)
        if verbose and (k + 1) % 10 == 0:
            print(f"  {k+1}/{len(tgs)}", flush=True)
    base = np.array([r["score"] for r in per], float)
    out = {"n_targets": len(per), "cells": [list(x) for x in cells],
           "base_score": float(base.mean()), "per_target": per, "summary": {}}
    for f, m in cells:
        key = f"{f}|{m}"
        d = {}
        med = np.array([r["cells"][key]["medoid"] for r in per], float)
        for nm in ("medoid", "avg", "fit", "fitms"):
            if nm not in per[0]["cells"][key]:
                continue
            v = np.array([r["cells"][key][nm] for r in per], float)
            d[nm] = {"sel": float(v.mean()),
                     "se": float(v.std(ddof=1) / math.sqrt(len(v))),
                     "median": float(np.median(v)),
                     "f2": float((v < 2.0).mean()), "f15": float((v < 1.5).mean()),
                     "vs_score": _paired(v, base), "vs_medoid": _paired(v, med)}
        for nm in ("fit_to_avg", "start_to_avg", "step", "step_sd", "fitms_to_avg"):
            if nm in per[0]["cells"][key]:
                d[nm] = float(np.mean([r["cells"][key][nm] for r in per]))
        out["summary"][key] = d
    _write("consensus2_fit.json", out)
    _print_fit(out)
    return out


def _print_fit(out):
    print(f"\n=== PROJECTING THE CONSENSUS BACK ONTO A REAL CHAIN — "
          f"{out['n_targets']} targets ===")
    print(f"  SHIPPED score {out['base_score']:.3f}")
    print(f"{'cell:form':16} {'medoid':>7} {'avg':>7} {'fit':>7} {'d(fit,avg)':>11} "
          f"{'CA-CA':>7} {'sd':>8} {'fit vs score':>13} {'95% CI':>19} {'W/L':>9}")
    for key, d in out["summary"].items():
        for nm in ("fit", "fitms"):
            if nm not in d:
                continue
            p = d[nm]["vs_score"]
            q = d[nm]["vs_medoid"]
            print(f"{key + ':' + nm:16} {d['medoid']['sel']:7.3f} {d['avg']['sel']:7.3f} "
                  f"{d[nm]['sel']:7.3f} {d['fit_to_avg']:11.3f} {d['step']:7.3f} "
                  f"{d['step_sd']:8.1e} {p['mean_diff']:+13.3f} "
                  f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}]".rjust(20) +
                  f" {p['n_better']:3d}/{p['n_worse']:<4d}"
                  f"  vs medoid {q['mean_diff']:+.3f} "
                  f"[{q['ci95'][0]:+.3f},{q['ci95'][1]:+.3f}] "
                  f"{q['n_better']}/{q['n_worse']}")


# ============================================================ stage: native
#: Cells the native-percentile diagnostic runs on.
NATIVE_CELLS = (("sc", 75), ("sc", 500), ("ref1", 100), ("rand", 75))
#: Criteria whose native percentile is measured.  Every one is native-free AS A CRITERION;
#: dropping the native into the set it ranks is the ORACLE part.
NAT_CRIT = ("medoid", "medoid_sq", "typic", "d_to_avg", "d_to_fit", "d_to_avgtrim",
            "score")


def native_percentiles(c, f, m):
    """ORACLE DIAGNOSTIC.  Where does the NATIVE fall under each consensus criterion?

    The coordinator's integration agent measured that the native sits at the 82.8th
    percentile (median 97.3) of the medoid criterion inside the score's top 75 -- i.e. the
    medoid is NOT finding the native-like candidate, it is returning the pool's MODE.  That
    bounds every operator that returns a pool MEMBER: if the criterion ranks the native at
    the 83rd percentile, no better filtering or better medoid rule can reach it.

    This reproduces that measurement on this instrument and asks whether ANY variant of the
    criterion -- squared, distance-matrix typicality, or distance to the SYNTHESISED
    consensus rather than to the other members -- places the native lower.  A variant that
    did would carry genuine nativeness information.

    The augmented pairwise matrix needs no new Kabsch work: the native's row is exactly
    `rr[sub]`, each candidate's CA-RMSD to the native, which the cache already holds.
    """
    n = c["n"]
    sub = filter_orders(c)[f][:m]
    W, rr, nat = c["W"][sub], c["rr"][sub], c["nat_ca"]
    Ps = c["P"][np.ix_(sub, sub)]
    B = len(sub)
    A = np.zeros((B + 1, B + 1))
    A[:B, :B] = Ps
    A[:B, B] = rr
    A[B, :B] = rr
    out = {}

    def pct(v):
        """Percentile of the native (last element) under a lower-is-better criterion."""
        v = np.asarray(v, float)
        return 100.0 * float((v[:-1] < v[-1]).mean()), bool(v[-1] <= v[:-1].min())

    out["medoid"] = pct(A.sum(1) / B)
    out["medoid_sq"] = pct((A ** 2).sum(1) / B)
    Dn = gen.pair_D(nat[None], n)
    Da = np.vstack([c["D"][sub], Dn])
    out["typic"] = pct(np.abs(Da - Da.mean(0)[None, :]).mean(1))
    # distance to the SYNTHESISED consensus, built from the candidates only
    C = _avg_struct(W, Ps)
    out["d_to_avg"] = pct(np.concatenate(
        [audit.kabsch_rmsd_batch(W, C), audit.kabsch_rmsd_batch(nat[None], C)]))
    keep = _trim_keep(Ps, 0.75)
    Ct = _avg_struct(W[keep], Ps[np.ix_(keep, keep)])
    out["d_to_avgtrim"] = pct(np.concatenate(
        [audit.kabsch_rmsd_batch(W, Ct), audit.kabsch_rmsd_batch(nat[None], Ct)]))
    # one generic start, not four: this is a percentile diagnostic, not a deliverable,
    # and the four-start form was measured to agree to 0.004 A on the full instrument
    F = fit_torsions(C, np.full(n, math.radians(-120.0)),
                     np.full(n, math.radians(130.0)), maxiter=150)[0]
    out["d_to_fit"] = pct(np.concatenate(
        [audit.kabsch_rmsd_batch(W, F), audit.kabsch_rmsd_batch(nat[None], F)]))
    pred = poolsize.load_pred(c["pdb"])
    out["score"] = pct(np.concatenate(
        [gen.shipped_score(pred, c["D"][sub]), gen.shipped_score(pred, Dn)]))
    # where synthesis sits relative to the MODE and the native
    b = int(np.argmin(_medoid(Ps)))
    geo_row = {
        "medoid_rr": float(rr[b]),
        "oracle_member": float(rr.min()),
        "avg_rr": float(audit.kabsch_rmsd_batch(C[None], nat)[0]),
        "fit_rr": float(audit.kabsch_rmsd_batch(F[None], nat)[0]),
        "avg_to_mode": float(audit.kabsch_rmsd_batch(C[None], W[b])[0]),
        "fit_to_mode": float(audit.kabsch_rmsd_batch(F[None], W[b])[0]),
        "avg_to_nearest": float(audit.kabsch_rmsd_batch(W, C).min()),
        "spread": float(Ps[np.triu_indices(B, 1)].mean())}
    return out, geo_row


def stage_native(cells=NATIVE_CELLS, verbose=True):
    tgs = cached()
    per = []
    for k, pdbid in enumerate(tgs):
        c = load_cache(pdbid)
        row = {"pdb": pdbid, "fold": c["fold"], "cells": {}}
        for f, m in cells:
            pcts, g = native_percentiles(c, f, m)
            row["cells"][f"{f}|{m}"] = {"pct": {a: v[0] for a, v in pcts.items()},
                                        "argmin": {a: v[1] for a, v in pcts.items()},
                                        "geo": g}
        per.append(row)
        if verbose and (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tgs)}", flush=True)
    out = {"n_targets": len(per), "cells": [list(x) for x in cells],
           "crits": list(NAT_CRIT), "per_target": per, "summary": {}}
    for f, m in cells:
        key = f"{f}|{m}"
        d = {"crit": {}}
        for a in NAT_CRIT:
            v = np.array([r["cells"][key]["pct"][a] for r in per], float)
            d["crit"][a] = {
                "mean_pct": float(v.mean()), "median_pct": float(np.median(v)),
                "se": float(v.std(ddof=1) / math.sqrt(len(v))),
                "frac_below_10": float((v < 10).mean()),
                "n_argmin": int(sum(r["cells"][key]["argmin"][a] for r in per))}
        for g in ("medoid_rr", "oracle_member", "avg_rr", "fit_rr", "avg_to_mode",
                  "fit_to_mode", "avg_to_nearest", "spread"):
            d[g] = float(np.mean([r["cells"][key]["geo"][g] for r in per]))
        out["summary"][key] = d
    _write("consensus2_native.json", out)
    _print_native(out)
    return out


def _print_native(out):
    print(f"\n=== WHERE DOES THE NATIVE FALL UNDER THE CONSENSUS CRITERION? — "
          f"{out['n_targets']} targets (ORACLE DIAGNOSTIC) ===")
    for key, d in out["summary"].items():
        print(f"  --- {key}")
        print(f"{'criterion':14} {'mean pct':>9} {'median':>8} {'+-se':>6} "
              f"{'<10th pct':>10} {'argmin':>7}")
        for a in out["crits"]:
            v = d["crit"][a]
            print(f"{a:14} {v['mean_pct']:9.1f} {v['median_pct']:8.1f} {v['se']:6.1f} "
                  f"{v['frac_below_10']:10.3f} {v['n_argmin']:4d}/{out['n_targets']}")
        print(f"    THE BOUND:  mode (medoid) {d['medoid_rr']:.3f} | "
              f"oracle best MEMBER {d['oracle_member']:.3f} | "
              f"synthesis avg {d['avg_rr']:.3f}, fit {d['fit_rr']:.3f}")
        print(f"    synthesis sits {d['avg_to_mode']:.3f} A from the mode "
              f"(fit {d['fit_to_mode']:.3f}), {d['avg_to_nearest']:.3f} A from the "
              f"nearest real member; subpool spread {d['spread']:.3f} A")


# ============================================================ stage: multi
#: The cell the multi-hypothesis study runs on.  Fixed to the incumbent filter so the
#: comparison is against a number already on the record.
MULTI_CELL = ("sc", 75)
MULTI_K = (1, 2, 3, 5)


def multi_one(c, f=MULTI_CELL[0], m=MULTI_CELL[1], ks=MULTI_K):
    """Consensus per cluster: best-of-k (ORACLE diagnostic) and native-free picks."""
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    rr = c["rr"]
    sub = filter_orders(c)[f][:m]
    Ps = c["P"][np.ix_(sub, sub)]
    rs = rr[sub]
    d = squareform(np.maximum((Ps + Ps.T) / 2.0, 0.0), checks=False)
    Z = linkage(d, method="average")
    nat = c["nat_ca"]
    row = {"pdb": c["pdb"], "fold": c["fold"], "score": sel_of(c["sc"], rr), "k": {}}
    for k in ks:
        lab = fcluster(Z, t=min(k, len(sub)), criterion="maxclust")
        hyp = {"medoid": [], "avg": []}
        size, mscore, tscore = [], [], []
        for lv in np.unique(lab):
            mem = np.where(lab == lv)[0]
            v = np.full(len(sub), np.inf)
            v[mem] = Ps[np.ix_(mem, mem)].sum(1) / max(len(mem) - 1, 1)
            hyp["medoid"].append(sel_of(v, rs))
            C = _avg_struct(c["W"][sub][mem], Ps[np.ix_(mem, mem)])
            hyp["avg"].append(float(audit.kabsch_rmsd_batch(C[None], nat)[0]))
            size.append(int(len(mem)))
            mscore.append(float(c["sc"][sub][mem].mean()))
            tscore.append(float(np.min(v[mem])))
        cell = {"n_clusters": len(size), "sizes": size}
        for op, h in hyp.items():
            h = np.array(h, float)
            sfx = "" if op == "medoid" else "_avg"
            cell[f"best_of_k{sfx}"] = float(h.min())              # ORACLE DIAGNOSTIC
            cell[f"worst_of_k{sfx}"] = float(h.max())
            cell[f"mean_of_k{sfx}"] = float(h.mean())
            cell[f"pick_largest{sfx}"] = float(h[int(np.argmax(size))])   # DEPLOYABLE
            cell[f"pick_score{sfx}"] = float(h[int(np.argmin(mscore))])   # DEPLOYABLE
            cell[f"pick_tight{sfx}"] = float(h[int(np.argmin(tscore))])   # DEPLOYABLE
        row["k"][str(k)] = cell
    return row


def stage_multi(verbose=True):
    tgs = cached()
    per = []
    for k, pdbid in enumerate(tgs):
        per.append(multi_one(load_cache(pdbid)))
        if verbose and (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tgs)}", flush=True)
    base = np.array([r["score"] for r in per], float)
    sw = _read("consensus2_sweep.json")
    inc = np.array([r["cells"][f"{INCUMBENT[0]}|{INCUMBENT[1]}"]["ops"][INCUMBENT[2]]
                    for r in sw["per_target"]], float)
    out = {"n_targets": len(per), "cell": list(MULTI_CELL), "ks": list(MULTI_K),
           "base_score": float(base.mean()), "incumbent_sel": float(inc.mean()),
           "per_target": per, "summary": {}}
    for k in MULTI_K:
        s = {}
        for key in ("best_of_k", "worst_of_k", "mean_of_k", "pick_largest",
                    "pick_score", "pick_tight",
                    "best_of_k_avg", "worst_of_k_avg", "mean_of_k_avg",
                    "pick_largest_avg", "pick_score_avg", "pick_tight_avg"):
            v = np.array([r["k"][str(k)][key] for r in per], float)
            s[key] = {"sel": float(v.mean()),
                      "se": float(v.std(ddof=1) / math.sqrt(len(v))),
                      "f2": float((v < 2.0).mean()),
                      "vs_score": _paired(v, base), "vs_incumbent": _paired(v, inc)}
        s["mean_n_clusters"] = float(np.mean([r["k"][str(k)]["n_clusters"] for r in per]))
        s["mean_largest_frac"] = float(np.mean(
            [max(r["k"][str(k)]["sizes"]) / sum(r["k"][str(k)]["sizes"]) for r in per]))
        out["summary"][str(k)] = s
    _write("consensus2_multi.json", out)
    _print_multi(out)
    return out


def _print_multi(out):
    print(f"\n=== MULTI-HYPOTHESIS — {out['n_targets']} targets, filter "
          f"{out['cell'][0]}{out['cell'][1]} ===")
    print(f"  SHIPPED score {out['base_score']:.3f}   single medoid "
          f"{out['incumbent_sel']:.3f}")
    for sfx, lbl in (("", "MEDOID per cluster"), ("_avg", "AVERAGE per cluster")):
        print(f"  --- {lbl}")
        print(f"{'k':>3} {'ncl':>5} {'lgfrac':>7} {'best-of-k*':>11} {'mean-of-k':>10} "
              f"{'worst':>8} {'largest':>9} {'byscore':>9} {'tightest':>9}")
        for k in out["ks"]:
            s = out["summary"][str(k)]
            print(f"{k:3d} {s['mean_n_clusters']:5.2f} {s['mean_largest_frac']:7.3f} "
                  f"{s['best_of_k' + sfx]['sel']:11.3f} "
                  f"{s['mean_of_k' + sfx]['sel']:10.3f} "
                  f"{s['worst_of_k' + sfx]['sel']:8.3f} "
                  f"{s['pick_largest' + sfx]['sel']:9.3f} "
                  f"{s['pick_score' + sfx]['sel']:9.3f} "
                  f"{s['pick_tight' + sfx]['sel']:9.3f}")
    print("  * best-of-k is an ORACLE DIAGNOSTIC: it reads the native to choose among the "
          "k hypotheses.")
    for k in out["ks"]:
        if k == 1:
            continue
        for key in ("pick_largest", "pick_score", "pick_tight",
                    "pick_largest_avg", "pick_score_avg", "pick_tight_avg"):
            p = out["summary"][str(k)][key]["vs_incumbent"]
            print(f"   k={k} {key:14} vs single medoid {p['mean_diff']:+.3f} "
                  f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] "
                  f"{p['n_better']}/{p['n_worse']}")


# ============================================================ stage: refine
def stage_refine(box_deg=(5.0, 15.0), n_targets=None, verbose=True):
    """Does a bounded Legacy relaxation of the CONSENSUS OUTPUT stack on top of it?

    `legacy_refine.refine_angles` is a bounded pattern search in torsion space, so it needs
    torsions: the consensus medoid is a real window carrying its parent's (phi, psi), and
    the refined structure is a REBUILD.  `reb_med` in the sweep is therefore the right
    zero-step comparator -- it isolates what the rebuild alone costs from what the
    relaxation adds.

    Sprint-6 work (`work/refine_study.py`, quoted in `legacy_refine`'s docstring) already
    prices this from oracle and native starts (+0.27 to +1.84 A); this measures it from the
    start we actually ship.  Amber relaxation is being measured concurrently by another
    agent in `s8/relax.py`; its committed `relax_sweep.json` shows k=100/10 restrained
    minimisation moves a candidate 0.0003-0.011 A of CA-RMSD, so it is not duplicated here.
    """
    import legacy_refine as lr
    tgs = cached()
    if n_targets:
        tgs = tgs[:n_targets]
    f, m, _ = INCUMBENT
    per = []
    for k, pdbid in enumerate(tgs):
        c = load_cache(pdbid)
        sub = filter_orders(c)[f][:m]
        Ps = c["P"][np.ix_(sub, sub)]
        gi = int(sub[int(np.argmin(_medoid(Ps)))])
        nat = c["nat_ca"]
        C = _avg_struct(c["W"][sub], Ps)
        F, fph, fps, _, _ = fit_consensus(C)
        starts = {"medoid": (c["PHI"][gi], c["PSI"][gi]), "fit": (fph, fps)}
        row = {"pdb": pdbid, "fold": c["fold"],
               "medoid": float(c["rr"][gi]),
               "avg": float(audit.kabsch_rmsd_batch(C[None], nat)[0]),
               "medoid_reb": float(audit.kabsch_rmsd_batch(
                   geo.build_backbone_batch(c["PHI"][gi][None],
                                            c["PSI"][gi][None])["CA"], nat)[0]),
               "fit": float(audit.kabsch_rmsd_batch(F[None], nat)[0])}
        for nm, (ph, ps) in starts.items():
            for bd in box_deg:
                r = lr.refine_angles(c["seq"], ph, ps, box_deg=bd, steps=40, pop=48,
                                     seed=0, polish=True)
                CA = np.asarray(r["CA"], float)
                row[f"{nm}_box{bd:g}"] = float(
                    audit.kabsch_rmsd_batch(CA[None], nat)[0])
                row[f"{nm}_dE{bd:g}"] = float(r["energy"] - r["energy0"])
        per.append(row)
        if verbose and (k + 1) % 10 == 0:
            print(f"  {k+1}/{len(tgs)}", flush=True)
    out = {"n_targets": len(per), "arm": f"{f}|{m}", "boxes": list(box_deg),
           "per_target": per, "starts": {}}
    for nm, base_key in (("medoid", "medoid_reb"), ("fit", "fit")):
        base = np.array([r[base_key] for r in per], float)
        real = np.array([r["medoid" if nm == "medoid" else "fit"] for r in per], float)
        d = {"start": float(base.mean()), "reference": float(real.mean()), "box": {}}
        for bd in box_deg:
            v = np.array([r[f"{nm}_box{bd:g}"] for r in per], float)
            d["box"][f"{bd:g}"] = {
                "sel": float(v.mean()),
                "se": float(v.std(ddof=1) / math.sqrt(len(v))),
                "f2": float((v < 2.0).mean()),
                "vs_start": _paired(v, base), "vs_reference": _paired(v, real),
                "dE": float(np.mean([r[f"{nm}_dE{bd:g}"] for r in per]))}
        out["starts"][nm] = d
    med = np.array([r["medoid"] for r in per], float)
    out["medoid"] = float(med.mean())
    out["avg"] = float(np.mean([r["avg"] for r in per]))
    out["fit"] = float(np.mean([r["fit"] for r in per]))
    out["rebuild_cost"] = _paired(
        np.array([r["medoid_reb"] for r in per], float), med)
    _write("consensus2_refine.json", out)
    _print_refine(out)
    return out


def _print_refine(out):
    print(f"\n=== DOES REFINEMENT STACK ON THE CONSENSUS OUTPUT? — {out['n_targets']} "
          f"targets, filter {out['arm']} ===")
    print(f"  consensus medoid (a real window)          {out['medoid']:.3f}")
    p = out["rebuild_cost"]
    print(f"  coordinate average (NOT a chain)          {out['avg']:.3f}")
    print(f"  torsion projection of that average        {out['fit']:.3f}")
    print(f"  the medoid's own ideal-geometry rebuild   "
          f"{out['starts']['medoid']['start']:.3f}  "
          f"({p['mean_diff']:+.3f} [{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] "
          f"{p['n_better']}/{p['n_worse']}) -- THE REBUILD COST")
    for nm, d in out["starts"].items():
        print(f"  --- Legacy relaxation started from the {nm.upper()} "
              f"(start {d['start']:.3f})")
        for bd, b in d["box"].items():
            a = b["vs_start"]
            print(f"      box {bd:>4} deg   {b['sel']:.3f} +- {b['se']:.3f}   "
                  f"{a['mean_diff']:+.3f} [{a['ci95'][0]:+.3f},{a['ci95'][1]:+.3f}]  "
                  f"{a['n_better']}/{a['n_worse']}   dE {b['dE']:+.1f} kcal   "
                  f"<2 A {b['f2']:.3f}")


# ============================================================ stage: dev
def stage_dev():
    """The single dev pass. One pre-registered arm (`DEV_ARM`), no iteration.

    The dev pools are `s8/inband_devpool`, built and asserted against `s7/audit_cache` by
    `s8.inband.stage_devpool`.  This module does not rebuild or alter them.
    """
    inband.stage_devpool(verbose=False)
    f, m, op = DEV_ARM
    #: the dev natives, used ONLY as reporting labels -- `inband_devpool` stores `rr` but
    #: not coordinates, and a CONSTRUCTED consensus needs a native to be scored against.
    nats = {q.pdb: np.asarray(q.ca, float) for q in debias.dev_targets()}
    rows = []
    for pdbid in sorted(x[:-4] for x in os.listdir(DCACHE) if x.endswith(".npz")):
        z = np.load(os.path.join(DCACHE, f"{pdbid}.npz"), allow_pickle=True)
        W = z["W"].astype(float)
        rr = z["rr"].astype(float)
        n = int(z["n"])
        pred = poolsize.load_pred(pdbid)
        D = gen.pair_D(W, n)
        sc = gen.shipped_score(pred, D)
        P = np.zeros((len(W), len(W)), float)
        for a in range(len(W)):
            P[a] = audit.kabsch_rmsd_batch(W, W[a])
        c = {"pdb": pdbid, "n": n, "fold": int(z["fold"]), "seq": str(z["seq"]),
             "W": W, "PHI": None, "PSI": None, "D": D, "sc": sc,
             "typ": np.abs(D - D.mean(0)[None, :]).mean(1), "P": P}
        sub = filter_orders(c)[f][:m]
        row = {"pdb": pdbid, "n": n, "fold": int(z["fold"]),
               "best": float(rr.min()), "mean": float(rr.mean()),
               "score": sel_of(sc, rr)}
        for nm in (op, "medoid"):
            if nm == "fit":
                Ps = c["P"][np.ix_(sub, sub)]
                F = fit_consensus(_avg_struct(W[sub], Ps))[0]
                row[nm] = float(audit.kabsch_rmsd_batch(F[None], nats[pdbid])[0])
                continue
            kind, v = op_value(nm, c, sub)
            row[nm] = (sel_of(v, rr[sub]) if kind == "sel" else
                       float(audit.kabsch_rmsd_batch(v[None], nats[pdbid])[0]))
        row["arm"] = row[op]
        rows.append(row)
    base = np.array([r["score"] for r in rows], float)
    v = np.array([r["arm"] for r in rows], float)
    med = np.array([r["medoid"] for r in rows], float)
    out = {"n_targets": len(rows), "arm": f"{f}|{m}:{op}", "per_target": rows,
           "medoid75": {"sel": float(med.mean()),
                        "se": float(med.std(ddof=1) / math.sqrt(len(med))),
                        "vs_score": _paired(med, base)},
           "pool_best": float(np.mean([r["best"] for r in rows])),
           "score": {"sel": float(base.mean()),
                     "se": float(base.std(ddof=1) / math.sqrt(len(base))),
                     "median": float(np.median(base)),
                     "f2": float((base < 2.0).mean())},
           "arm_stats": {"sel": float(v.mean()),
                         "se": float(v.std(ddof=1) / math.sqrt(len(v))),
                         "median": float(np.median(v)),
                         "f2": float((v < 2.0).mean()),
                         "vs_score": _paired(v, base)}}
    _write("consensus2_dev.json", out)
    print(f"\n=== DEV-24, ONE PRE-REGISTERED PASS — arm {out['arm']} ===")
    print(f"  SHIPPED score {out['score']['sel']:.3f} +- {out['score']['se']:.3f}")
    p = out["arm_stats"]["vs_score"]
    q = out["medoid75"]["vs_score"]
    print(f"  {'sc|75:medoid (S8-8)':22} {out['medoid75']['sel']:.3f} +- "
          f"{out['medoid75']['se']:.3f}   d {q['mean_diff']:+.3f} "
          f"[{q['ci95'][0]:+.3f},{q['ci95'][1]:+.3f}]  {q['n_better']}/{q['n_worse']}")
    print(f"  {out['arm']:22} {out['arm_stats']['sel']:.3f} +- "
          f"{out['arm_stats']['se']:.3f}   d {p['mean_diff']:+.3f} "
          f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}]  {p['n_better']}/{p['n_worse']}")
    return out


# ============================================================ stage: leak
def stage_leak(n_targets=12):
    """Every deployable arm must be bit-identical when the native is replaced by NaN.

    The filtered SUBSET is recomputed from the poisoned cache too -- if a filter key had
    read a native the subset itself would change, which is a stronger test than fixing the
    index set first.
    """
    bad, worst = [], 0.0
    for pdbid in cached()[:n_targets]:
        c = load_cache(pdbid)
        p = dict(c)
        p["rr"] = np.full_like(c["rr"], np.nan)
        p["nat_ca"] = np.full_like(c["nat_ca"], np.nan)
        oa, ob = filter_orders(c), filter_orders(p)
        for f in FILTERS:
            if not np.array_equal(oa[f], ob[f]):
                bad.append((pdbid, f, "filter order differs"))
            for m in MS:
                sa, sb = oa[f][:m], ob[f][:m]
                for op in OPS:
                    ka, va = op_value(op, c, sa)
                    kb, vb = op_value(op, p, sb)
                    va, vb = np.asarray(va, float), np.asarray(vb, float)
                    if not np.isfinite(vb[np.isfinite(va)]).all():
                        bad.append((pdbid, f"{f}{m}:{op}", "nan"))
                    fin = np.isfinite(va) & np.isfinite(vb)
                    d = float(np.abs(va[fin] - vb[fin]).max()) if fin.any() else 0.0
                    if d > 0:
                        bad.append((pdbid, f"{f}{m}:{op}", d))
                    worst = max(worst, d)
    out = {"n_targets": min(n_targets, len(cached())),
           "n_arms": len(FILTERS) * len(MS) * len(OPS),
           "worst_abs_diff": worst, "violations": bad[:50], "clean": not bad}
    _write("consensus2_leak.json", out)
    print(f"leakage poison test: {out['n_targets']} targets x {out['n_arms']} arms, "
          f"worst |diff| {worst:.3e}, {'CLEAN' if out['clean'] else 'VIOLATIONS'}")
    return out


# ============================================================ report
def stage_report():
    for name, fn in (("consensus2_sweep.json", _print_sweep),
                     ("consensus2_fit.json", _print_fit),
                     ("consensus2_native.json", _print_native),
                     ("consensus2_multi.json", _print_multi),
                     ("consensus2_refine.json", _print_refine)):
        p = os.path.join(HERE, name)
        if os.path.exists(p):
            fn(_read(name))


STAGES = {"build": stage_build, "sweep": stage_sweep, "fit": stage_fit,
          "native": stage_native, "multi": stage_multi, "refine": stage_refine,
          "dev": stage_dev, "leak": stage_leak, "report": stage_report}

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in STAGES:
        print(f"usage: python -m s8.consensus2 [{'|'.join(STAGES)}]")
        raise SystemExit(2)
    STAGES[sys.argv[1]]()
