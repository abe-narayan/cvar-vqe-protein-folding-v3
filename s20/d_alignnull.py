"""SPRINT 20, AGENT D, BLOCK C -- the null for the shared-bias alignment statistic.

Pre-registered in `s20/PREREG_D.md` Block C, falsifier F-D1 fixed before this ran.

`s19/a_source.py` reports cross-family alignment of the COHERENT error component as
`corrC(A,B) = corr(d(X_A) - d_nat, d(X_B) - d_nat)`, and `s19/CLAIMS.md` S2 reads
sepprior 0.575 / helix 0.598 against a same-seed ceiling of 0.898 as "two thirds of the
harmful component is reproduced by ZERO-INFORMATION references".

That statistic has a null nobody subtracted.  Both arguments are scored against the SAME
native, so every family's coherent residual contains the identical term -(d_nat - d_typ).
Two structures that share nothing but the target will correlate.  This module measures how
much, three ways:

    N1 pool x pool    two independent windows from this target's own K=500 retrieval pool
    N2 fit x fit      the SAME s15/align_lib.fit, from the SAME cached start, with the SAME
                      1/sd^2 weights, driven to the distance field of a randomly drawn pool
                      window -- two independent draws.  MATCHED IN THE OPERATOR'S SPACE.
    N3 cross-target   a window from a DIFFERENT target of the same length

ORACLE: `nat` is read to score.  No inference decision reads a native quantity.
"""
from __future__ import annotations
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I          # noqa: E402
from s15 import distcal as C             # noqa: E402
from s15 import align_lib as A           # noqa: E402
from s15 import seed as SD               # noqa: E402
from s19 import a_fit as F               # noqa: E402

OUT = os.path.join(HERE, "results", "D_C_ALIGNNULL")
os.makedirs(OUT, exist_ok=True)
SALT = "s20_D_alignnull"


def corr(a, b):
    a = np.asarray(a, float) - np.mean(a)
    b = np.asarray(b, float) - np.mean(b)
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else np.nan


def pdists(X, i, j):
    X = np.asarray(X, float)
    rv = X[..., i, :] - X[..., j, :]
    return np.sqrt((rv * rv).sum(-1))


def run(nmax=None):
    tg = I.targets()
    if nmax:
        tg = tg[:nmax]
    data = C.gather(tg)
    by_n = {}
    for t in I.targets():
        by_n.setdefault(int(t["n"]), []).append(t["pdb"])

    rows = []
    t0 = time.time()
    path = os.path.join(OUT, "alignnull.json")
    for c, t in enumerate(tg):
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        d = data[pdb]
        i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
        rng = SD.stable_rng(pdb, "alignnull", salt=SALT)
        phi0, psi0, _avg = F.start(pdb, seq, fold)

        u = I.load_univ(pdb)
        pool = I.pool_idx(u)                       # the shipped K=500 BLOSUM pool
        W = np.asarray(u["W"], float)[pool]
        d_nat = pdists(nat, i, j)

        # R independent DRAWS per target, averaged, so the reported per-target value is not
        # one lucky pair of windows.  N1 and N2 use the SAME draws, so their difference
        # isolates the fit operator and nothing else.
        R = 8
        n1s, n2s, n3s, fitq = [], [], [], []
        others = [p for p in by_n.get(n, []) if p != pdb]
        Wo_cache = {}
        for _r in range(R):
            k1, k2 = rng.choice(len(W), size=2, replace=False)
            # ---- N1: two independent pool windows, no fit at all
            n1s.append(corr(pdists(W[k1], i, j) - d_nat, pdists(W[k2], i, j) - d_nat))
            # ---- N2: the SAME fit, SAME start, SAME weights, driven to those SAME two
            #          windows' distance fields.  Matched in the operator's space.
            fits = []
            for k in (k1, k2):
                tgt_d = np.maximum(pdists(W[k], i, j), 2.0)
                ph, ps, _f = A.fit(tgt_d, sd, i, j, phi0, psi0)
                Xk = I.build_ca(ph, ps)
                fitq.append(float(I.ca_rmsd(Xk, W[k])))
                fits.append(pdists(Xk, i, j) - d_nat)
            n2s.append(corr(fits[0], fits[1]))
            # ---- N3: windows of the same length from DIFFERENT targets
            if len(others) >= 2:
                o1, o2 = rng.choice(len(others), size=2, replace=False)
                X = []
                for o in (others[o1], others[o2]):
                    if o not in Wo_cache:
                        uo = I.load_univ(o)
                        Wo_cache[o] = np.asarray(uo["W"], float)[I.pool_idx(uo)]
                    Wo = Wo_cache[o]
                    X.append(pdists(Wo[rng.integers(0, len(Wo))], i, j) - d_nat)
                n3s.append(corr(X[0], X[1]))
        n1 = float(np.nanmean(n1s))
        n2 = float(np.nanmean(n2s))
        n3 = float(np.nanmean(n3s)) if n3s else np.nan

        # ---- diagnostics that price the mechanism
        # the native's own deviation from a sequence-blind separation profile
        sep = np.asarray(d["sep"], float)
        pool_d = pdists(W, i, j)                      # (500, npairs)
        d_typ = pool_d.mean(0)                        # this target's pool-mean field
        rows.append({"pdb": pdb, "n": n, "fold": fold, "npairs": int(len(i)),
                     "N1_poolxpool": n1, "N2_fitxfit": n2, "N3_crosstarget": float(n3),
                     "var_nat_dev": float(np.var(d_nat - d_typ)),
                     "var_pool_dev": float(np.mean(np.var(pool_d - d_typ, axis=0))),
                     "fit_reach_rmsd": float(np.mean(fitq)),
                     "sep_mean": float(sep.mean())})
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            json.dump({"rows": rows, "complete": False}, open(path, "w"))
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)

    json.dump({"rows": rows, "complete": len(rows) == len(tg), "n_expected": len(tg)},
              open(path, "w"))
    if len(rows) == len(tg) and nmax is None:
        open(os.path.join(OUT, "COMPLETE"), "w").write(f"n={len(rows)}\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(OUT, "alignnull.json")))["rows"]
    src = json.load(open(os.path.join(ROOT, "s19", "results", "a_source.json")))["rows"]
    src = {r["pdb"]: r for r in src}
    pdbs = [r["pdb"] for r in rows]
    folds = np.array([r["fold"] for r in rows])
    rng = np.random.default_rng(20200907)

    def boot(x, B=4000):
        x = np.asarray(x, float)
        x = x[np.isfinite(x)]
        k = len(x)
        m = x[rng.integers(0, k, size=(B, k))].mean(1)
        return float(x.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)), k

    g = lambda k: np.array([r[k] for r in rows], float)                       # noqa: E731
    s = lambda k: np.array([src[p][k] if p in src else np.nan for p in pdbs], float)  # noqa: E731

    print(f"\nn = {len(rows)}.  corrC(A,B) = corr(d(X_A) - d_nat, d(X_B) - d_nat).\n")
    print(f"  {'quantity':<44}{'mean':>8}{'95% CI':>20}{'n':>5}")
    nulls = {}
    for k, lab in (("N1_poolxpool", "NULL N1  two pool windows, no fit"),
                   ("N2_fitxfit", "NULL N2  same fit/start/weights, MATCHED"),
                   ("N3_crosstarget", "NULL N3  windows from other targets")):
        m, lo, hi, kk = boot(g(k))
        nulls[k] = (m, lo, hi)
        print(f"  {lab:<44}{m:>8.3f}   [{lo:+.3f},{hi:+.3f}]{kk:>5}")
    print()
    fam = {}
    for k, lab in (("corrC|d20|d20_s1", "same arch, different seed -- the CEILING"),
                   ("corrC|deployed|d20", "same arch, different regularisation"),
                   ("corrC|deployed|pairnet", "DIFFERENT architecture"),
                   ("corrC|deployed|poolmean", "retrieval pool mean"),
                   ("corrC|deployed|sepprior", "ZERO-INFO sequence-blind separation prior"),
                   ("corrC|deployed|helix", "ZERO-INFO constant alpha-helix")):
        m, lo, hi, kk = boot(s(k))
        fam[k] = m
        print(f"  {lab:<44}{m:>8.3f}   [{lo:+.3f},{hi:+.3f}]{kk:>5}")

    print(f"\n  EXCESS OVER THE MATCHED NULL N2 (paired, per target)\n")
    print(f"  {'pair':<44}{'raw':>8}{'excess':>9}{'95% CI':>20}{'W/L':>10}{'folds':>7}")
    exc = {}
    n2 = g("N2_fitxfit")
    for k, lab in (("corrC|d20|d20_s1", "same arch, different seed -- the CEILING"),
                   ("corrC|deployed|d20", "same arch, different regularisation"),
                   ("corrC|deployed|pairnet", "DIFFERENT architecture"),
                   ("corrC|deployed|poolmean", "retrieval pool mean"),
                   ("corrC|deployed|sepprior", "ZERO-INFO separation prior"),
                   ("corrC|deployed|helix", "ZERO-INFO constant alpha-helix")):
        dif = s(k) - n2
        m, lo, hi, kk = boot(dif)
        good = np.isfinite(dif)
        sg = sum(1 for f in np.unique(folds)
                 if np.sign(np.nanmean(dif[folds == f])) == np.sign(m))
        exc[k] = {"raw": fam[k], "excess": m, "ci": [lo, hi],
                  "W": int((dif[good] > 0).sum()), "L": int((dif[good] < 0).sum()),
                  "folds": int(sg)}
        print(f"  {lab:<44}{fam[k]:>8.3f}{m:>+9.3f}   [{lo:+.3f},{hi:+.3f}]"
              f"{int((dif[good]>0).sum()):>6}/{int((dif[good]<0).sum()):<4}{sg:>5}/5")

    ceil = exc["corrC|d20|d20_s1"]["excess"]
    print(f"\n  SHARE OF THE CEILING, before and after subtracting the null:")
    print(f"  {'pair':<44}{'as published':>14}{'null-subtracted':>18}")
    for k, lab in (("corrC|deployed|d20", "same arch, different regularisation"),
                   ("corrC|deployed|pairnet", "DIFFERENT architecture"),
                   ("corrC|deployed|poolmean", "retrieval pool mean"),
                   ("corrC|deployed|sepprior", "ZERO-INFO separation prior"),
                   ("corrC|deployed|helix", "ZERO-INFO constant alpha-helix")):
        pub = fam[k] / fam["corrC|d20|d20_s1"]
        adj = exc[k]["excess"] / ceil if ceil != 0 else np.nan
        print(f"  {lab:<44}{pub:>14.2f}{adj:>18.2f}")

    m2, lo2, hi2, _ = boot(g("N2_fitxfit"))
    print(f"\n  F-D1: null N2 = {m2:.3f} [{lo2:.3f},{hi2:.3f}]")
    if hi2 < 0.40:
        print("  F-D1 DOES NOT FIRE: null below 0.40, the attack FAILS, S2 survives as stated.")
    elif lo2 <= 0.575 <= hi2 or m2 >= 0.575:
        print("  F-D1 FIRES at the strong clause: the zero-information references are AT "
              "the null.")
    else:
        print("  F-D1 fires at the weak clause: the null is substantial and must be "
              "subtracted; the references are not AT it.")

    json.dump({"nulls": nulls, "families": fam, "excess": exc, "n": len(rows)},
              open(os.path.join(OUT, "alignnull_report.json"), "w"), indent=1)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.report:
        report()
    else:
        run(a.n)
