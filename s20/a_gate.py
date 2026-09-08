"""SPRINT 20, AGENT A -- is it the DISTOGRAM, or any score?  And is the agreement more than
"both are typical peptides"?

`a_sel.py` found that swapping the CORPUS leaves the emitted error aligned at 0.951-0.968 with the
incumbent while swapping the SELECTOR drops it to 0.79.  Two follow-ups, both cheap, both POST-HOC
(neither was pre-registered):

  1. THIRD SELECTOR.  Rank the 500 candidates by Ca-RMSD to a constant ideal alpha-helix -- a
     ZERO-INFORMATION score that knows nothing about the target beyond its length (brief section 6
     rule 4; Sprint 19 C4 measured this as the WORST gate on the board).  If a helix gate also
     pulls the emitted error onto the incumbent's, the carrier is "any score ordering", not the
     distogram.  If only the distogram does, the selector is a specific lever.

  2. PARTIAL CORRELATION ON THE ZERO-INFORMATION MODE.  Every arm here emits a typical peptide of
     the target's length, and the constant alpha-helix alone aligns with the incumbent at 0.814.
     So a raw rho of 0.96 may be mostly that generic mode.  `rho_partial(e_A, e_B | e_helix)` asks
     how much the two agree BEYOND what a structure with no information agrees with them about.
     This is the honest version of the alignment statistic and it is reported beside the raw one.

BASIS: POINT CLOUD throughout (`I.coordinate_average`), stated because it is contracted.

Run:  python -m s20.a_gate
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402
from s20 import a_src as S                   # noqa: E402

OUT = os.path.join(S.RESULTS, "a_gate.json")
CORP = ["blosum", "univ"]
SEL = ["disto", "helixgate", "rand"]
ARMS = [f"{c}.{s}" for c in CORP for s in SEL]
K, M = 500, 75


def _partial(a, b, z):
    """Pearson(a, b) with z partialled out of both."""
    A = np.stack([a, b, z]).astype(float)
    A = A - A.mean(1, keepdims=True)
    ra = A[0] - (A[0] @ A[2]) / max(A[2] @ A[2], 1e-12) * A[2]
    rb = A[1] - (A[1] @ A[2]) / max(A[2] @ A[2], 1e-12) * A[2]
    d = np.sqrt((ra @ ra) * (rb @ rb))
    return float(ra @ rb / d) if d > 0 else np.nan


def run(tg=None, limit=None):
    tg = tg if tg is not None else I.targets()
    if limit:
        tg = tg[:limit]
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        u = I.load_univ(pdb)
        nat = u["nat_ca"]
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        dtrue = np.linalg.norm(nat[i] - nat[j], axis=1)
        order = np.asarray(u["order"], int)
        W_all = np.asarray(u["W"], float)
        rr = np.asarray(u["rr"], float)
        nw = len(W_all)
        hel = I.build_ca(np.full(n, np.deg2rad(-63.0)), np.full(n, np.deg2rad(-42.0)))
        e_hel = np.linalg.norm(hel[i] - hel[j], axis=1) - dtrue

        idx = {"blosum": order[:K],
               "univ": SD.stable_rng(pdb, "s20A_univ").choice(nw, size=min(K, nw),
                                                              replace=False)}
        e = {"pdb": pdb, "n": n, "fold": fold}
        X = {}
        for cp in CORP:
            ix = idx[cp]
            W = W_all[ix]
            D = I.pair_dists(W, i, j)
            sc = I.shipped_score(dg, D.astype(np.float32).astype(float))
            hg = I.kabsch_rmsd_batch(W, hel)      # ZERO-INFORMATION gate: closeness to a helix
            pick = {"disto": np.argsort(sc, kind="stable")[:M],
                    "helixgate": np.argsort(hg, kind="stable")[:M],
                    "rand": SD.stable_rng(pdb, cp, "s20A_selrand").choice(
                        len(ix), size=min(M, len(ix)), replace=False)}
            for sl in SEL:
                a = f"{cp}.{sl}"
                p = pick[sl]
                X[a], _b = I.coordinate_average(W[p])
                e[f"rmsd|{a}"] = float(I.ca_rmsd(X[a], nat))
                e[f"member_rmsd|{a}"] = float(rr[ix][p].mean())
                e[f"sel_ceiling|{a}"] = float(rr[ix][p].min())
                e[f"bond|{a}"] = float(np.linalg.norm(X[a][1:] - X[a][:-1], axis=1).mean())
        E = {a: np.linalg.norm(X[a][i] - X[a][j], axis=1) - dtrue for a in ARMS}
        E["helix"] = e_hel
        names = ARMS + ["helix"]
        Z = S._std(np.stack([E[a] for a in names]))
        C = Z @ Z.T
        for a in range(len(names)):
            for b in range(a + 1, len(names)):
                e[f"rho|{names[a]}|{names[b]}"] = float(C[a, b])
                if "helix" not in (names[a], names[b]):
                    e[f"rhoP|{names[a]}|{names[b]}"] = _partial(E[names[a]], E[names[b]], e_hel)
                    e[f"xrmsd|{names[a]}|{names[b]}"] = float(I.ca_rmsd(X[names[a]], X[names[b]]))
        rows.append(e)
        if (c + 1) % 25 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"),
                      default=float)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(S.RESULTS, "a_gate.COMPLETE"), "w").write("ok\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    folds = np.array([r["fold"] for r in rows], int)
    pdbs = [r["pdb"] for r in rows]
    inc = "blosum.disto"
    print(f"\n=== SELECTOR IDENTITY   n = {len(rows)}   BASIS: POINT CLOUD ===")
    print(f"  {'arm':<18}{'realis':>8}{'member':>8}{'selceil':>9}{'bond':>7}"
          f"{'rho vs inc':>11}{'rhoP vs inc':>12}")
    for a in ARMS + ["helix"]:
        if a == "helix":
            k = f"rho|{inc}|helix"
            print(f"  {'helix (zero-info)':<18}{'-':>8}{'-':>8}{'-':>9}{3.804:>7.3f}"
                  f"{np.nanmean(g(k)):>11.3f}{'-':>12}")
            continue
        k = f"rho|{inc}|{a}"
        kp = f"rhoP|{inc}|{a}"
        rv = f"{np.nanmean(g(k)):>11.3f}" if k in rows[0] else f"{1.0:>11.3f}"
        rp = f"{np.nanmean(g(kp)):>12.3f}" if kp in rows[0] else f"{1.0:>12.3f}"
        print(f"  {a:<18}{g('rmsd|'+a).mean():>8.3f}{g('member_rmsd|'+a).mean():>8.3f}"
              f"{g('sel_ceiling|'+a).mean():>9.3f}{g('bond|'+a).mean():>7.3f}{rv}{rp}")
    print("\n  --- raw rho vs rho PARTIALLED on the zero-information helix mode ---")
    for a, b in [("blosum.disto", "univ.disto"), ("blosum.disto", "blosum.helixgate"),
                 ("blosum.disto", "univ.helixgate"), ("blosum.disto", "blosum.rand"),
                 ("blosum.disto", "univ.rand"), ("blosum.helixgate", "univ.helixgate"),
                 ("blosum.rand", "univ.rand")]:
        k = f"rho|{a}|{b}" if f"rho|{a}|{b}" in rows[0] else f"rho|{b}|{a}"
        kp = f"rhoP|{a}|{b}" if f"rhoP|{a}|{b}" in rows[0] else f"rhoP|{b}|{a}"
        xk = f"xrmsd|{a}|{b}" if f"xrmsd|{a}|{b}" in rows[0] else f"xrmsd|{b}|{a}"
        print(f"  {a+' ~ '+b:<40}rho {np.nanmean(g(k)):>6.3f}   partial "
              f"{np.nanmean(g(kp)):>6.3f}   RMSD {g(xk).mean():>5.3f}")
    print("\n  --- realised, paired vs the incumbent (POINT CLOUD basis) ---")
    base = g("rmsd|" + inc)
    for a in ARMS:
        if a == inc:
            continue
        st = I.paired(g("rmsd|" + a), base, folds=folds, names=pdbs)
        print(f"  {a:<18}{g('rmsd|'+a).mean():>7.3f}  {st['mean_diff']:+.3f} "
              f"[{st['ci95'][0]:+.3f},{st['ci95'][1]:+.3f}]  {st['n_better']}W/{st['n_worse']}L")
    return rows


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run(limit=int(sys.argv[1]) if len(sys.argv) > 1 else None)
