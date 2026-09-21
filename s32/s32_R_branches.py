#!/usr/bin/env python
"""S32 LANE R -- the projection-branch census, its ORACLE ceiling, and the native-free
criteria that might pick a branch.

PREREG: s32/PREREG_S32_R.md, committed at 02754f5a before the first number existed.

BASIS.  Built-chain CA-RMSD, tuning126, n=126.  The chain of a cloud C is
`fit_prior(... lam=0) -> fit_prior(... pen, lam=0.3)`, the same two rungs production walks,
emitted by `core.project`'s bit-exact builder.  PRODUCTION (`PROD`) is
`s12.instrument.project(C, seq, fold)["ca"]`, i.e. `lam_path(multi=True)` verbatim, run IN
THIS JOB so that every contrast is paired inside one process (contract rule 3).

THE CLOUD.  `s29/results/s29_O_structs/<pdb>.npz["prod"]` is the object the contract's
3.2105 / 3.0483 are defined on.  It differs from `bench_results/cache/<PROD_KEY>/<pdb>.json
["avg_ca"]` by ~1e-14 on 126/126 targets -- a recomputation, not a copy -- and that is NOT
nothing on this stage: the module's own note prices the amplification at ~1e13.  Both
clouds are projected here (`PROD` and `PROD_CACHECLOUD`) so the floor is measured, not
assumed.

WHAT A "BRANCH" IS.  One converged (phi, psi) solution of the projection, reached from one
start.  Start families (all native-free to construct):
    GEN4    the 4 production generic constant-torsion starts       -> this is production
    GEN4D   the same 4 taken DIRECTLY to lam=0.3 (production's `fit_multi` at the top rung,
            so production's whole candidate set is a subset of this job's branch set)
    MEM75   the phi/psi of the production top-75 filtered pool members
    RANDk   75 windows drawn uniformly from the same target's universe, 3 pinned draws

ORACLE.  `rmsd_nat` per branch is a LABEL.  Nothing in the start construction, the branch
set, or any criterion reads it.  Every arm that selects on it is labelled
ORACLE / NOT DEPLOYABLE in the artefact key itself.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                                          # noqa: E402
from core import project as pj                                           # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")
STRUCTS = os.path.join(RESULTS, "s32_R_structs")
S29_STRUCTS = os.path.join(ROOT, "s29", "results", "s29_O_structs")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(STRUCTS, exist_ok=True)

LAM = 0.3
MAXITER = 300
GRAD = "exact"
N_RAND = 75
RAND_DRAWS = 3
RAND_SEED0 = 32_000_017          # pinned here, in the committed source


# ------------------------------------------------------------------ inputs
def prod_cloud(pdb):
    with np.load(os.path.join(S29_STRUCTS, f"{pdb}.npz"), allow_pickle=True) as z:
        return (np.asarray(z["prod"], float), str(z["seq"]), int(z["fold"]), int(z["n"]))


def top75(u, pdb):
    """(indices into the universe of the production top-75, medoid position within it).

    Reproduced, not trusted: the 75 members are re-averaged with the production operator
    and the result is checked against the persisted cloud.
    """
    pool = I.pool_idx(u)
    sub = np.asarray(I.shipped_record(pdb)["sub"], int)
    idx = pool[sub]
    W = np.asarray(u["W"], float)[idx].astype(np.float32).astype(float)   # reference f32
    P = I.pairwise_rmsd(W).astype(np.float32).astype(float)
    b = int(np.argmin(P.sum(1) / max(P.shape[0] - 1, 1)))
    C = I.superpose_batch(W, W[b]).mean(0)
    return idx, b, C, W


# ------------------------------------------------------------------ branches
def starts_for(u, idx, n, rng_seed):
    """Every start, as (family, tag, phi0, psi0).  No native anywhere."""
    out = []
    for a, b in pj.STARTS:
        out.append(("GEN4", f"gen_{a:.0f}_{b:.0f}",
                    np.full(n, math.radians(a)), np.full(n, math.radians(b))))
    PHI = np.asarray(u["PHI"], float)
    PSI = np.asarray(u["PSI"], float)
    for r, k in enumerate(idx):
        out.append(("MEM75", f"mem_{r}", PHI[k].copy(), PSI[k].copy()))
    nw = len(PHI)
    for d in range(RAND_DRAWS):
        rng = np.random.default_rng(rng_seed + 7919 * d)
        pick = rng.choice(nw, size=min(N_RAND, nw), replace=False)
        for r, k in enumerate(pick):
            out.append((f"RAND{d}", f"rand{d}_{r}", PHI[k].copy(), PSI[k].copy()))
    return out


def run_branches(C, seq, fold, n, starts):
    """Two rungs per start -- lam=0 then continue to lam=0.3 -- plus the DIRECT lam=0.3 fit
    for the four generic starts, which is the other half of production's candidate set."""
    pen = pj.make_penalty("ramah", seq, int(fold))
    rows, PH, PS = [], [], []
    for fam, tag, p0, s0 in starts:
        g0 = pj.fit_prior(C, p0, s0, pen=None, lam=0.0, maxiter=MAXITER, grad=GRAD)
        g1 = pj.fit_prior(C, g0[1], g0[2], pen=pen, lam=LAM, maxiter=MAXITER, grad=GRAD)
        rows.append(dict(fam=fam, tag=tag, obj0=float(g0[3]), obj1=float(g1[3]),
                         d_to_C=float(g1[5]), d0_to_C=float(g0[5])))
        PH.append(np.asarray(g1[1], float))
        PS.append(np.asarray(g1[2], float))
    for a, b in pj.STARTS:                                        # direct top-rung fits
        g = pj.fit_prior(C, np.full(n, math.radians(a)), np.full(n, math.radians(b)),
                         pen=pen, lam=LAM, maxiter=MAXITER, grad=GRAD)
        rows.append(dict(fam="GEN4D", tag=f"gend_{a:.0f}_{b:.0f}", obj0=float("nan"),
                         obj1=float(g[3]), d_to_C=float(g[5]), d0_to_C=float("nan")))
        PH.append(np.asarray(g[1], float))
        PS.append(np.asarray(g[2], float))
    return rows, np.array(PH), np.array(PS)


# ------------------------------------------------------------------ native-free criteria
def criteria(seq, fold, pdb, PH, PS, CA, C, Wmem):
    """Every native-free per-branch score.  Fixed in the PREREG before any number existed."""
    n = len(seq)
    out = {}
    rama = pj.RamaPenalty(seq, int(fold), kind="rama")
    rama20 = pj.RamaPenalty(seq, int(fold), kind="rama20")
    ramah = pj.RamaHingePenalty(seq, int(fold))
    out["rama_nlp"] = np.asarray(rama(PH, PS), float)
    out["rama20_nlp"] = np.asarray(rama20(PH, PS), float)
    out["ramah"] = np.asarray(ramah(PH, PS), float)
    gly = np.array([c == "G" for c in seq])
    m = pj._interior_mask(n) & ~gly
    if not m.any():
        m = pj._interior_mask(n)
    out["posphi_frac"] = (PH[:, m] > 0).mean(1).astype(float)

    dg = I.distogram(pdb, seq, int(fold))
    i, j = I.pair_index(n)
    D = I.pair_dists(CA, i, j)
    out["disto_risk"] = np.asarray(I.shipped_score(dg, D), float)
    out["disto_mae"] = np.abs(D - np.asarray(dg["expected"], float)[None, :]).mean(1)

    #: Legacy's fitted-weight combination, from the branch's own torsions
    import core
    lf = core.backend("legacy")
    from core import geometry as geo
    BB = geo.build_backbone_batch(PH, PS)

    class _Rep:
        def __init__(s, A, B):
            s._phi = np.asarray(A, float)[:1].T.copy()
            s._psi = np.asarray(B, float)[:1].T.copy()
            s.n_states = 1
            s.n_residues = s._phi.shape[0]
    bl = lf.BatchLegacy(seq, _Rep(PH, PS))
    T = bl.terms_from_coords(BB, phi=PH, psi=PS)
    wv = np.array([lf.FITTED_WEIGHTS.get(t, 0.0) for t in lf.TERMS], float)
    out["legacy"] = np.asarray(T, float) @ wv

    #: typicality -- mean CA-RMSD of the branch to the production top-75 members
    out["typicality"] = np.array([I.kabsch_rmsd_batch(Wmem, ca).mean() for ca in CA])
    #: geometry secondaries, carried from row one (contract rule 15)
    vb = np.linalg.norm(CA[:, 1:] - CA[:, :-1], axis=-1)
    out["vbond_mean"] = vb.mean(1)
    out["vbond_sd"] = vb.std(1)
    cen = CA - CA.mean(1, keepdims=True)
    out["rg"] = np.sqrt((cen ** 2).sum(-1).mean(1))
    out["d_to_cloud"] = I.kabsch_rmsd_batch(CA, C)
    return out


# ------------------------------------------------------------------ one target
def one(pdb, verbose=True):
    t0 = time.time()
    C, seq, fold, n = prod_cloud(pdb)
    u = I.load_univ(pdb)
    assert str(u["seq"]) == seq and int(u["fold"]) == fold and int(u["n"]) == n
    nat = np.asarray(u["nat_ca"], float)                         # ORACLE: label only
    idx, medoid_pos, Crec, Wmem = top75(u, pdb)
    cloud_recon_err = float(np.abs(Crec - C).max())

    #: the cache's own cloud, 1e-14 away from the one the contract is defined on
    Ccache = np.asarray(I.shipped_record(pdb)["avg_ca"], float)
    cloud_cache_delta = float(np.abs(Ccache - C).max())

    #: PRODUCTION, run here, both clouds
    pr = I.project(C, seq, fold)
    pr_cache = I.project(Ccache, seq, fold)

    starts = starts_for(u, idx, n, RAND_SEED0 + (abs(hash(pdb)) % 1000))
    rows, PH, PS = run_branches(C, seq, fold, n, starts)
    CA = np.asarray(pj.build_ca_exact(PH, PS), float)
    crit = criteria(seq, fold, pdb, PH, PS, CA, C, Wmem)
    rmsd = I.kabsch_rmsd_batch(CA, nat)                          # ORACLE LABEL
    d_prod = I.kabsch_rmsd_batch(CA, np.asarray(pr["ca"], float))

    #: distinct branches at the threshold `core.project.degeneracy` calls `n_same_point`
    B = len(CA)
    P = np.zeros((B, B))
    for a in range(B):
        P[a] = I.kabsch_rmsd_batch(CA, CA[a])
    #: CONNECTED COMPONENTS of the `P < 1e-3` adjacency, by union-find.
    #: BUG FIXED 2026-09-21 (found by lane V's adversary pass, before any n_distinct was
    #: quoted).  The previous line was `lab[P[a] < 1e-3] = nc`, which OVERWRITES labels that
    #: an earlier seed had already assigned -- it counted greedy seeds, not components, and
    #: it stole members from earlier clusters.  `n_distinct` underwrites this lane's claim
    #: that the branches are structurally distinct, so it has to be the real quantity.
    par = np.arange(B)

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for a in range(B):
        for b in np.where(P[a] < 1e-3)[0]:
            ra, rb = find(a), find(int(b))
            if ra != rb:
                par[ra] = rb
    roots = np.array([find(a) for a in range(B)])
    _, lab = np.unique(roots, return_inverse=True)
    nc = int(lab.max()) + 1

    out = dict(
        pdb=pdb, n=n, fold=fold, seq=seq,
        cloud_rmsd=float(I.ca_rmsd(C, nat)),
        cloud_cache_rmsd=float(I.ca_rmsd(Ccache, nat)),
        cloud_cache_delta=cloud_cache_delta, cloud_recon_err=cloud_recon_err,
        cloud_vbond_mean=float(np.linalg.norm(C[1:] - C[:-1], axis=-1).mean()),
        cloud_vbond_sd=float(np.linalg.norm(C[1:] - C[:-1], axis=-1).std()),
        cloud_rg=float(np.sqrt(((C - C.mean(0)) ** 2).sum(-1).mean())),
        prod_chain=float(I.ca_rmsd(pr["ca"], nat)),
        prod_fit_chain=float(I.ca_rmsd(pr["fit_ca"], nat)),
        prod_cachecloud_chain=float(I.ca_rmsd(pr_cache["ca"], nat)),
        prod_d_to_cloud=float(I.ca_rmsd(pr["ca"], C)),
        medoid_pos=medoid_pos,
        n_branches=B, n_distinct=int(nc),
        fam=[r["fam"] for r in rows], tag=[r["tag"] for r in rows],
        cluster=lab.tolist(),
        obj0=[r["obj0"] for r in rows], obj1=[r["obj1"] for r in rows],
        d_to_C=[r["d_to_C"] for r in rows], d0_to_C=[r["d0_to_C"] for r in rows],
        rmsd_nat=rmsd.tolist(), d_to_prod=d_prod.tolist(),
        secs=round(time.time() - t0, 2),
    )
    for k, v in crit.items():
        out[k] = np.asarray(v, float).tolist()
    np.savez_compressed(os.path.join(STRUCTS, f"{pdb}.npz"),
                        phi=PH.astype(np.float32), psi=PS.astype(np.float32),
                        fam=np.array([r["fam"] for r in rows]),
                        cluster=lab, rmsd_nat=rmsd.astype(np.float32),
                        prod_ca=np.asarray(pr["ca"], float), seq=seq, fold=fold, n=n)
    if verbose:
        print("%-6s n=%2d B=%3d distinct=%3d prod=%.4f cloud=%.4f "
              "branch[min %.4f med %.4f max %.4f] %.1fs"
              % (pdb, n, B, nc, out["prod_chain"], out["cloud_rmsd"],
                 rmsd.min(), np.median(rmsd), rmsd.max(), out["secs"]), flush=True)
    return out


def rows_path(shard):
    return os.path.join(RESULTS, f"s32_R_branches_shard{int(shard)}.jsonl")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)
    pdbs = [t["pdb"] for t in I.targets()]
    mine = [p for k, p in enumerate(pdbs) if k % a.nshard == a.shard]
    if a.limit:
        mine = mine[:a.limit]
    path = rows_path(a.shard)
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["pdb"])
                except Exception:                                        # noqa: BLE001
                    pass
    for p in mine:
        if p in done:
            continue
        r = one(p)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(r) + "\n")
    print("shard %d done: %d targets -> %s" % (a.shard, len(mine), path), flush=True)


if __name__ == "__main__":
    main()
