#!/usr/bin/env python
"""s27/s28_C_readout.py -- S28 lane C, Part 2: readouts that consume a ranker's in-pool
information without the uniform average's bottleneck (`s27/PREREG_S28_C.md` section 2 and
addendum 2).

Every readout keeps the PRODUCTION retained set (the shipped score's top-75, ties by the stable
per-target key) and changes only how the set is consumed.  All three are convex combinations
of the same 75 members (stated in addendum 2):

  (a) MEDNB(k)        the top-75 member with the best ranker value plus its k nearest
                      neighbours by CA-RMSD inside the top-75; uniform average of the k+1
                      (`readout_uniform`, their own medoid frame).  k = 74 is production.
  (b) TRIM(q)         drop the ceil(q*75) members with the worst ranker value; uniform average
                      of the rest.  q = 0 is production.
  (c) DIVW(beta,gam)  all 75 superposed on the production medoid; weighted mean with
                      w_i ~ exp(-beta z_i) * (sum_j exp(-P_ij^2 / 2 sigma^2))^(-gamma),
                      z_i = zrank of the ranker within the 75, sigma = median off-diagonal
                      P.  (0, 0) is production to floating point.

Controls, matched in the operator's space: the same operator with the ranker rank-permuted
within the 75 by a stable per-target key.  NATIVE-FREE: no function here reads `nat_ca` or
`oracle_rr`; the ORACLE endpoint is evaluated only in `oracle_rows`.

    python s27/s28_C_readout.py selftest
    python s27/s28_C_readout.py cloud  [--limit N]          # every cell, point cloud
    python s27/s28_C_readout.py chain  [--limit N] [--arms ...]   # built chain, primary cells
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s27 import ham_lib as HL              # noqa: E402
from s27 import run_pool as RP             # noqa: E402

RESULTS = os.path.join(HERE, "results")
CLOUD_ROWS = os.path.join(RESULTS, "s28_C_readout_cloud_rows.jsonl")
CHAIN_ROWS = os.path.join(RESULTS, "s28_C_readout_chain_rows.jsonl")
M = 75
KS = (10, 20, 40)
QS = (0.05, 0.10, 0.20)
BETAS = (0.5, 1.0, 2.0)
#: (beta, gamma) cells: the identity, the two one-factor references, the beta ladder at gamma=1
DIVW_CELLS = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.5, 1.0), (1.0, 1.0), (2.0, 1.0))
PRIMARY_CHAIN = ["PROD", "MEDNB[CONS,k=20]", "MEDNB[CONS~perm,k=20]",
                 "TRIM[CONS,q=0.1]", "TRIM[CONS~perm,q=0.1]", "TRIM[DISTPOT,q=0.1]", "TRIM[DISTPOT~perm,q=0.1]",
                 "DIVW[CONS,b=1,g=1]", "DIVW[CONS~perm,b=1,g=1]", "DIVW[CONS,b=1,g=0]", "DIVW[CONS~perm,b=1,g=0]",
                 "DIVW[CONS,b=0,g=1]"]


# ============================================================================ the readouts
def production_set(ch, pdb):
    k = len(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(k)
    return RP.topm(ch["DIS"], M, key), key


def _fmt(x):
    return ("%g" % x)


def readout_mednb(W75, P, r, k):
    """(a): the member with the lowest ranker value r, plus its k nearest by P; uniform average
    of those k+1 in their own medoid frame.  Returns (C, member indices)."""
    m = len(W75)
    k = int(min(k, m - 1))
    seed = int(np.argmin(r))               # r carries no exact ties by construction (see below)
    order = np.argsort(P[seed], kind="stable")
    nb = [int(a) for a in order if a != seed][:k]
    idx = np.array([seed] + nb, int)
    C, _ = I.coordinate_average(W75[idx], P[np.ix_(idx, idx)])
    return C, idx


def readout_trim(W75, P, r, q):
    """(b): drop the ceil(q*m) members with the highest ranker value; uniform average of the
    rest in their own medoid frame."""
    m = len(W75)
    nd = int(math.ceil(q * m - 1e-9)) if q > 0 else 0
    keep = np.argsort(-r, kind="stable")[nd:]
    keep = np.sort(keep)
    C, _ = I.coordinate_average(W75[keep], P[np.ix_(keep, keep)])
    return C, keep


def divw_weights(P, r, beta, gamma):
    """(c): w_i ~ exp(-beta z_i) * dens_i^(-gamma), z = zrank(r) within the set, dens_i =
    sum_j exp(-P_ij^2 / (2 sigma^2)), sigma = median off-diagonal P.  Normalised to sum 1."""
    m = len(r)
    z = HL.zrank(r) if beta != 0 else np.zeros(m)
    off = P[np.triu_indices(m, 1)]
    sigma = float(np.median(off)) if len(off) else 1.0
    sigma = max(sigma, 1e-6)
    dens = np.exp(-(P ** 2) / (2.0 * sigma ** 2)).sum(1)
    logw = -beta * z - gamma * np.log(dens)
    logw = logw - logw.max()
    w = np.exp(logw)
    return w / w.sum()


def readout_divw(W75, P, r, beta, gamma):
    """(c): all members superposed on the production medoid (the medoid of P), weighted mean."""
    b = I.medoid(P)
    Wm = I.superpose_batch(W75, W75[b])
    w = divw_weights(P, r, beta, gamma)
    return (w[:, None, None] * Wm).sum(0), w


def ranker_within(ch, top, name, pdb, perm=False):
    """The ranker's values over the retained set; `perm` rank-permutes them within the set by a
    stable per-target key (the matched control).  Exact ties are broken by a stable key so no
    argmin ever reads array order."""
    r = np.asarray(ch[name], float)[top]
    if perm:
        r = r[RP.rng_for(pdb, f"s28C-perm-{name}").permutation(len(r))]
    key = RP.rng_for(pdb, f"s28C-tie-{name}").random(len(r))
    # a strictly ordered surrogate: rank with ties broken by the key (monotone in r)
    order = np.lexsort((key, r))
    rr = np.empty(len(r))
    rr[order] = np.arange(len(r), dtype=float)
    return rr


def all_arms():
    arms = [("PROD", None)]
    for k in KS:
        for perm in (False, True):
            arms.append((f"MEDNB[CONS{'~perm' if perm else ''},k={k}]", ("mednb", "CONS", perm, k)))
    for rk in ("CONS", "DISTPOT"):
        for q in QS:
            for perm in (False, True):
                arms.append((f"TRIM[{rk}{'~perm' if perm else ''},q={_fmt(q)}]", ("trim", rk, perm, q)))
    for beta, gamma in DIVW_CELLS:
        if beta == 0 and gamma == 0:
            arms.append((f"DIVW[CONS,b=0,g=0]", ("divw", "CONS", False, (beta, gamma))))
            continue
        arms.append((f"DIVW[CONS,b={_fmt(beta)},g={_fmt(gamma)}]", ("divw", "CONS", False, (beta, gamma))))
        if beta > 0:
            arms.append((f"DIVW[CONS~perm,b={_fmt(beta)},g={_fmt(gamma)}]", ("divw", "CONS", True, (beta, gamma))))
    return arms


def emit(cand, ch, arm_spec, top=None, P=None):
    """NATIVE-FREE emission of one arm's point cloud.  Returns (C, info)."""
    if top is None:
        top, _ = production_set(ch, cand.pdb)
    W75 = np.asarray(cand.W, float)[top]
    if P is None:
        P = I.pairwise_rmsd(W75)
    if arm_spec is None:
        C, b = I.coordinate_average(W75, P)
        return C, dict(n_members=M, medoid=int(top[b]))
    kind, rk, perm, par = arm_spec
    r = ranker_within(ch, top, rk, cand.pdb, perm)
    if kind == "mednb":
        C, idx = readout_mednb(W75, P, r, par)
        return C, dict(n_members=int(len(idx)), seed_member=int(top[idx[0]]))
    if kind == "trim":
        C, keep = readout_trim(W75, P, r, par)
        return C, dict(n_members=int(len(keep)))
    if kind == "divw":
        C, w = readout_divw(W75, P, r, par[0], par[1])
        ess = float(1.0 / (w ** 2).sum())
        return C, dict(n_members=M, ess=ess, w_max=float(w.max()))
    raise KeyError(kind)


# ============================================================================ ORACLE evaluation
def oracle_rows(pdb, arms, mode):
    """ORACLE: `nat_ca` is read here and nowhere else.  Point cloud for every arm; the built
    chain (production projection) when mode == 'chain'."""
    cand, ch, _ = RP.channels_for(pdb)
    top, _ = production_set(ch, pdb)
    W75 = np.asarray(cand.W, float)[top]
    P = I.pairwise_rmsd(W75)
    rows = []
    for name, spec in arms:
        t1 = time.time()
        C, info = emit(cand, ch, spec, top, P)
        row = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), arm=name, basis_cloud="point_cloud",
                   rmsd_cloud=float(I.ca_rmsd(C, cand.nat_ca)), **info)
        if mode == "chain":
            ca = H.readout_projected(cand, C)
            row["rmsd_chain"] = float(I.ca_rmsd(ca, cand.nat_ca))
            row["basis_chain"] = "built_chain"
        row["secs"] = float(time.time() - t1)
        rows.append(row)
    return rows


def run(mode, limit=0, arms=None):
    from s25 import phys_lib as P
    path = CLOUD_ROWS if mode == "cloud" else CHAIN_ROWS
    all_ = all_arms()
    if arms:
        want = set(arms)
        all_ = [a for a in all_ if a[0] in want]
        missing = want - {a[0] for a in all_}
        assert not missing, f"unknown arms {missing}"
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                    done.add((r["pdb"], r["arm"]))
                except Exception:
                    pass
    pdbs = P.targets()[:limit] if limit else P.targets()
    t0 = time.time()
    for n, pdb in enumerate(pdbs):
        todo = [a for a in all_ if (pdb, a[0]) not in done]
        if not todo:
            continue
        rows = oracle_rows(pdb, todo, mode)
        with open(path, "a", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        print(f"  [{mode} {n+1}/{len(pdbs)}] {pdb} {len(rows)} arms {sum(r['secs'] for r in rows):.1f}s "
              f"(elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", path)


# ============================================================================ analysis
def _load_rows(path):
    by = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            by.setdefault(r["arm"], {})[r["pdb"]] = r
    return by


def analyse(mode, print_all=True):
    """Every arm against production (paired, `ST.compare`), against its matched-random control,
    the FAIL18 / non-FAIL18 strata, and `ST.best_of_k_within` on each family's grid."""
    from s24 import stats_lib as ST
    path = CLOUD_ROWS if mode == "cloud" else CHAIN_ROWS
    key = "rmsd_cloud" if mode == "cloud" else "rmsd_chain"
    basis = "POINT CLOUD" if mode == "cloud" else "BUILT CHAIN"
    by = _load_rows(path)
    pdbs = sorted(by["PROD"])
    complete = [a for a in by if all(p in by[a] for p in pdbs)]
    folds = np.array([int(by["PROD"][p]["fold"]) for p in pdbs])
    fail = np.array([p in I.FAIL18 for p in pdbs])                  # ORACLE stratum label
    prod = np.array([by["PROD"][p][key] for p in pdbs])
    out = {"mode": mode, "basis": basis, "n": len(pdbs), "prod_mean": float(prod.mean()), "arms": {}}
    if mode == "cloud":
        rand = {}
        with open(os.path.join(RESULTS, "pool_rows.jsonl"), encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                if r["config"] == "DIS":
                    rand[r["pdb"]] = r["rand_mean"]
        out["random75_mean"] = float(np.mean([rand[p] for p in pdbs]))
    for a in complete:
        if a == "PROD":
            continue
        v = np.array([by[a][p][key] for p in pdbs])
        r = ST.compare(v, prod, folds, names=pdbs, label=f"{a} vs PROD ({basis})")
        d = v - prod
        rec = {"mean": float(v.mean()), "vs_prod": {k: x for k, x in r.items() if k != "concentration"},
               "concentration": r["concentration"],
               "fail18_effect": float(d[fail].mean()), "fail18_se": float(d[fail].std(ddof=1) / math.sqrt(fail.sum())),
               "non_fail18_effect": float(d[~fail].mean()), "non_fail18_se": float(d[~fail].std(ddof=1) / math.sqrt((~fail).sum()))}
        ctrl = a.replace("[CONS,", "[CONS~perm,").replace("[DISTPOT,", "[DISTPOT~perm,")
        if ctrl != a and ctrl in complete:
            vc = np.array([by[ctrl][p][key] for p in pdbs])
            rc = ST.compare(v, vc, folds, names=pdbs, label=f"{a} vs its permuted-ranker control ({basis})")
            rec["vs_control"] = {k: x for k, x in rc.items() if k != "concentration"}
            rec["control_mean"] = float(vc.mean())
        if "ess" in by[a][pdbs[0]]:
            rec["ess_mean"] = float(np.mean([by[a][p]["ess"] for p in pdbs]))
        if "n_members" in by[a][pdbs[0]]:
            rec["n_members"] = int(by[a][pdbs[0]]["n_members"])
        out["arms"][a] = rec
        if print_all:
            print(ST.fmt(r))
            print(f"    FAIL18 {rec['fail18_effect']:+.4f} (SE {rec['fail18_se']:.4f})   non-FAIL18 {rec['non_fail18_effect']:+.4f} (SE {rec['non_fail18_se']:.4f})"
                  + (f"   vs control: effect {rec['vs_control']['effect']:+.4f} MDE {rec['vs_control']['mde']:.4f} fold CI [{rec['vs_control']['ci95_fold'][0]:+.4f}, {rec['vs_control']['ci95_fold'][1]:+.4f}]" if "vs_control" in rec else ""))
    # grid pricing per family (real-ranker cells only; the identity is in every grid)
    fams = {"MEDNB": [f"MEDNB[CONS,k={k}]" for k in KS],
            "TRIM_CONS": [f"TRIM[CONS,q={_fmt(q)}]" for q in QS],
            "TRIM_DISTPOT": [f"TRIM[DISTPOT,q={_fmt(q)}]" for q in QS],
            "DIVW": [f"DIVW[CONS,b={_fmt(b)},g={_fmt(g)}]" for b, g in DIVW_CELLS if not (b == 0 and g == 0)]}
    out["grids"] = {}
    for fam, cells in fams.items():
        cells = [c for c in cells if c in complete]
        if len(cells) < 2:
            continue
        Mx = np.column_stack([prod] + [np.array([by[c][p][key] for p in pdbs]) for c in cells])
        w = ST.best_of_k_within(Mx)
        out["grids"][fam] = {"cells": ["PROD"] + cells, **{k: v for k, v in w.items()}}
        if print_all:
            print(f"  grid {fam}: oracle-min gain {w['observed_gain']:+.4f}  split-half {w['split_half']:+.4f} ({100*w['split_half_frac']:.0f}%)  k_eff {w['k_eff']:.1f}  -> {w['verdict']}")
    ST.save_atomic(os.path.join(RESULTS, f"s28_C_readout_{mode}_summary.json"), out, module_file=__file__)
    return out


def selftest():
    """Synthetic pool: identities reproduce the uniform average; the permuted control is a
    different operator; weights are a convex combination."""
    rng = np.random.default_rng(1)
    m, n = 75, 12
    base = np.cumsum(rng.normal(size=(n, 3)) * 2.0, 0)
    W = base[None] + rng.normal(scale=0.8, size=(m, n, 3))
    P = I.pairwise_rmsd(W)
    r = rng.random(m)
    C0, _ = I.coordinate_average(W, P)
    Ca, idx = readout_mednb(W, P, r, 74)
    assert len(idx) == 75 and np.allclose(Ca, C0, atol=1e-9)
    Cb, keep = readout_trim(W, P, r, 0.0)
    assert len(keep) == 75 and np.allclose(Cb, C0, atol=1e-9)
    Cc, w = readout_divw(W, P, r, 0.0, 0.0)
    assert np.allclose(w, 1.0 / m) and np.allclose(Cc, C0, atol=1e-9)
    Cd, wd = readout_divw(W, P, r, 1.0, 1.0)
    assert abs(wd.sum() - 1) < 1e-12 and (wd > 0).all() and not np.allclose(Cd, C0)
    Ce, ke = readout_trim(W, P, r, 0.1)
    assert len(ke) == 67
    Cf, idf = readout_mednb(W, P, r, 20)
    assert len(idf) == 21 and idf[0] == int(np.argmin(r))
    print("  s28_C_readout selftest OK")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["selftest", "cloud", "chain", "analyse_cloud", "analyse_chain"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--arms", default="", help="';'-separated arm names (names contain commas)")
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    if a.mode == "selftest":
        selftest()
    elif a.mode.startswith("analyse_"):
        analyse(a.mode.split("_", 1)[1])
    else:
        arms = [x for x in a.arms.split(";") if x] if a.arms else (PRIMARY_CHAIN if a.mode == "chain" else None)
        run(a.mode, a.limit, arms)


if __name__ == "__main__":
    main()
