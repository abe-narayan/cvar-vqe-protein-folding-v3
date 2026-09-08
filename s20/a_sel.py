"""SPRINT 20, AGENT A -- CORPUS or SELECTOR?  A 4x2 factorial on where the shared bias enters.

The coordinator's constraint ("any Q2 test whose candidate space is fitted to the retrieval pool
is asking a foregone question") exposes a confound in `a_src.py` that I did not pre-register:
every arm there swaps the CORPUS but keeps the SELECTOR, and the selector is the shipped
leave-fold-out distogram, which is trained on the same library the corpus comes from.  So a null
in `a_src` cannot distinguish "the corpus is not the carrier" from "the corpus is not the carrier
BECAUSE the selector is".

This module crosses them.  Nothing here was pre-registered; it is reported as a POST-HOC
DECOMPOSITION, not as a test of F-A1.

    corpus    blosum  order[:500]            fitted to: the target sequence (BLOSUM62), library
              pep     org == True,  top 500  fitted to: the target sequence; PEPTIDE corpus only
              prot    org == False, top 500  fitted to: the target sequence; PROTEIN fragments only
              univ    500 uniform at random  fitted to: NOTHING (length only)

    selector  disto   shipped Bayes-risk distogram score, lowest 75   fitted to: the library
              rand    75 uniformly at random                          fitted to: NOTHING

`blosum x disto` is the incumbent and reproduces the shipped top-75 exactly (gate in `a_src`).
Every emitted object is a POINT CLOUD (contracted); the basis is stated on every table.

Run:  python -m s20.a_sel
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

OUT = os.path.join(S.RESULTS, "a_sel.json")
CORP = ["blosum", "pep", "prot", "univ"]
SEL = ["disto", "rand"]
ARMS = [f"{c}.{s}" for c in CORP for s in SEL]
K, M = 500, 75


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
        org = np.asarray(u["org"], bool)
        W_all = np.asarray(u["W"], float)
        rr = np.asarray(u["rr"], float)
        nw = len(org)

        idx = {"blosum": order[:K], "pep": order[org[order]][:K],
               "prot": order[~org[order]][:K],
               "univ": SD.stable_rng(pdb, "s20A_univ").choice(nw, size=min(K, nw),
                                                              replace=False)}
        e = {"pdb": pdb, "n": n, "fold": fold}
        X, EM = {}, {}
        for cp in CORP:
            ix = idx[cp]
            W = W_all[ix]
            D = I.pair_dists(W, i, j)
            sc = I.shipped_score(dg, D.astype(np.float32).astype(float))
            pick = {"disto": np.argsort(sc, kind="stable")[:M],
                    "rand": SD.stable_rng(pdb, cp, "s20A_selrand").choice(
                        len(ix), size=min(M, len(ix)), replace=False)}
            for sl in SEL:
                a = f"{cp}.{sl}"
                p = pick[sl]
                Wt = W[p]
                X[a], _b = I.coordinate_average(Wt)
                EM[a] = D[p] - dtrue[None]
                P = I.pairwise_rmsd(Wt)
                iu = np.triu_indices(len(P), 1)
                bl = np.linalg.norm(X[a][1:] - X[a][:-1], axis=1)
                e[f"rmsd|{a}"] = float(I.ca_rmsd(X[a], nat))
                e[f"member_rmsd|{a}"] = float(rr[ix][p].mean())
                e[f"gen_ceiling|{a}"] = float(rr[ix].min())
                e[f"sel_ceiling|{a}"] = float(rr[ix][p].min())
                e[f"cov2|{a}"] = float((rr[ix][p] < 2.0).mean())
                e[f"div_geo|{a}"] = float(P[iu].mean())
                e[f"bond|{a}"] = float(bl.mean())
        Es = S._std(np.stack([np.linalg.norm(X[a][i] - X[a][j], axis=1) - dtrue for a in ARMS]))
        C = Es @ Es.T
        for a in range(len(ARMS)):
            for b in range(a + 1, len(ARMS)):
                e[f"rho|{ARMS[a]}|{ARMS[b]}"] = float(C[a, b])
                e[f"xrmsd|{ARMS[a]}|{ARMS[b]}"] = float(I.ca_rmsd(X[ARMS[a]], X[ARMS[b]]))
        for a in ARMS:
            Z = S._std(EM[a])
            iu = np.triu_indices(len(Z), 1)
            e[f"memrho_in|{a}"] = float((Z @ Z.T)[iu].mean())
            e[f"memrho_vs_inc|{a}"] = float((Z @ S._std(EM["blosum.disto"]).T).mean())
        rows.append(e)
        if (c + 1) % 20 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"),
                      default=float)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(S.RESULTS, "a_sel.COMPLETE"), "w").write("ok\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    folds = np.array([r["fold"] for r in rows], int)
    pdbs = [r["pdb"] for r in rows]
    print(f"\n=== CORPUS x SELECTOR   n = {len(rows)}   BASIS: POINT CLOUD ===")
    print(f"  {'arm':<16}{'realis':>8}{'member':>8}{'genceil':>9}{'selceil':>9}"
          f"{'cov<2A':>8}{'divGEO':>8}{'bond':>7}{'rho vs inc':>11}{'memrho':>9}")
    for a in ARMS:
        k = (f"rho|blosum.disto|{a}" if f"rho|blosum.disto|{a}" in rows[0] else None)
        rv = np.nanmean(g(k)) if k else 1.0
        print(f"  {a:<16}{g('rmsd|'+a).mean():>8.3f}{g('member_rmsd|'+a).mean():>8.3f}"
              f"{g('gen_ceiling|'+a).mean():>9.3f}{g('sel_ceiling|'+a).mean():>9.3f}"
              f"{g('cov2|'+a).mean():>8.3f}{g('div_geo|'+a).mean():>8.3f}"
              f"{g('bond|'+a).mean():>7.3f}{rv:>11.3f}{g('memrho_vs_inc|'+a).mean():>9.3f}")
    print("\n  --- what moves the emitted error, corpus or selector? (rho, POINT CLOUD) ---")
    pairs = [("blosum.disto", "pep.disto", "CORPUS   swap, selector held (disto)"),
             ("pep.disto", "prot.disto", "CORPUS   partition, selector held (disto)"),
             ("blosum.disto", "blosum.rand", "SELECTOR swap, corpus held (blosum)"),
             ("pep.disto", "pep.rand", "SELECTOR swap, corpus held (pep)"),
             ("blosum.disto", "univ.disto", "corpus -> uninformative, selector held"),
             ("blosum.disto", "univ.rand", "BOTH removed"),
             ("blosum.rand", "univ.rand", "corpus swap with NO selector anywhere"),
             ("pep.rand", "prot.rand", "CORPUS partition with NO selector anywhere")]
    for a, b, lab in pairs:
        k = f"rho|{a}|{b}" if f"rho|{a}|{b}" in rows[0] else f"rho|{b}|{a}"
        xk = f"xrmsd|{a}|{b}" if f"xrmsd|{a}|{b}" in rows[0] else f"xrmsd|{b}|{a}"
        print(f"  {a+' ~ '+b:<32}rho {np.nanmean(g(k)):>6.3f}   "
              f"RMSD(X_A,X_B) {g(xk).mean():>5.3f}   {lab}")
    print("\n  --- realised, paired vs the incumbent arm (POINT CLOUD basis) ---")
    base = g("rmsd|blosum.disto")
    for a in ARMS[1:]:
        st = I.paired(g("rmsd|" + a), base, folds=folds, names=pdbs)
        print(f"  {a:<16}{st['mean_diff']:+.3f} [{st['ci95'][0]:+.3f},{st['ci95'][1]:+.3f}]"
              f"  {st['n_better']}W/{st['n_worse']}L")
    return rows


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run(limit=int(sys.argv[1]) if len(sys.argv) > 1 else None)
