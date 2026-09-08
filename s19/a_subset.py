"""SPRINT 19, AGENT A, Q1b -- is the HARMFUL error in a NATIVE-FREE-DETECTABLE subset?

The pre-registered question: repair 25% of the pairs (set `dhat := dtrue` on them) and ask
whether any NATIVE-FREE criterion for choosing the 25% beats a size-matched RANDOM choice.

Two families of arm, and both are needed:

  raw       repair the subset, leave everything else alone.  Different criteria remove
            different amounts of squared error, so a criterion can win purely by picking the
            big residuals.  Each arm reports its achieved residual RMS.
  matched   the same subset repaired, then the SURVIVING residual rescaled so that the
            per-target residual RMS equals `real_rms * sqrt(0.75)` for EVERY arm -- the RMS
            a random 25% repair achieves in expectation.  This is the Sprint-18 G2e design:
            identical residual magnitude, different spatial pattern.

CRITERIA (all native-free except the labelled ORACLE ceiling):
  nf_model    leave-fold-out linear model of |r| from native-free features (a_topo, R2 0.28)
  hi_sd       least confident quartile          lo_sd    most confident quartile
  long_sep    largest separation                short_sep smallest separation
  terminal    closest to a chain terminus       multimode highest predictive entropy
  random      the mandatory size-matched control
  ORACLE_absr the true largest |r| -- the ceiling of perfect detection, never a result

Run:  python -m s19.a_subset
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

from s19 import a_lib as L                   # noqa: E402
from s19 import a_fit as F                   # noqa: E402
from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

OUT = os.path.join(L.RESULTS, "a_subset.json")
Q = 0.25
CRIT = ["nf_model", "hi_sd", "lo_sd", "long_sep", "short_sep", "terminal", "multimode",
        "random", "ORACLE_absr"]
FEATS = ["sep", "term", "sd", "ent", "maxp", "nmode", "near1", "dhat", "deg_m",
         "pred_contact", "hyd_m", "chg_prod", "vol_m", "hel_m", "she_m", "ss_same",
         "ss_H", "ss_E"]


def nf_error_models(data, pdbs, folds):
    """Leave-fold-out linear predictors of |r| from native-free features."""
    out = {}
    for f in sorted(set(np.asarray(folds).tolist())):
        tr = [p for p in pdbs if data[p]["fold"] != f]
        X = np.column_stack([np.concatenate([data[p][k] for p in tr]) for k in FEATS])
        y = np.concatenate([data[p]["absr"] for p in tr])          # ORACLE label, in-fold only
        mu, sg = X.mean(0), np.maximum(X.std(0), 1e-9)
        Z = np.column_stack([(X - mu) / sg, np.ones(len(X))])
        coef, *_ = np.linalg.lstsq(Z, y, rcond=None)
        out[f] = (mu, sg, coef)
    return out


def run(tg=None):
    tg = tg if tg is not None else I.targets()
    data, pdbs, folds = L.gather_all(tg)
    nfm = nf_error_models(data, pdbs, folds)
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        d = data[pdb]
        i, j, sd, nat, n = d["i"], d["j"], d["sd"], d["nat"], d["n"]
        dhat, dtrue = d["dhat"], d["dtrue"]
        r = dhat - dtrue
        P = len(r)
        k = max(1, int(round(Q * P)))
        rms = float(np.sqrt((r ** 2).mean()))
        target_rms = rms * np.sqrt(1.0 - Q)
        rng = SD.stable_rng(pdb, "s19A_subset")
        phi0, psi0, avg = F.start(pdb, d["seq"], d["fold"])
        mu, sg, coef = nfm[d["fold"]]
        Z = np.column_stack([(np.column_stack([d[f] for f in FEATS]) - mu) / sg, np.ones(P)])
        pred = Z @ coef

        scores = {
            "nf_model": pred, "hi_sd": sd, "lo_sd": -sd, "long_sep": d["sep"],
            "short_sep": -d["sep"], "terminal": -d["term"], "multimode": d["ent"],
            "random": rng.random(P), "ORACLE_absr": np.abs(r),
        }
        e = {"pdb": pdb, "n": n, "fold": d["fold"], "npairs": P, "k": k,
             "avg": float(I.ca_rmsd(avg, nat)), "resid_rms": rms}
        rm, _ = F.fit_rmsd(dhat, sd, i, j, phi0, psi0, nat)
        e["real"] = rm
        for cr in CRIT:
            S = np.argsort(-np.asarray(scores[cr], float))[:k]
            rr = r.copy()
            rr[S] = 0.0
            e[cr] = F.fit_rmsd(dtrue + rr, sd, i, j, phi0, psi0, nat)[0]
            e[cr + "_rms"] = float(np.sqrt((rr ** 2).mean()))
            s = target_rms / max(float(np.sqrt((rr ** 2).mean())), 1e-9)
            e[cr + "_m"] = F.fit_rmsd(dtrue + rr * s, sd, i, j, phi0, psi0, nat)[0]
        rows.append(e)
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"),
                      default=float)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(L.RESULTS, "a_subset.COMPLETE"), "w").write("ok\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rng = SD.stable_rng("s19A", "subset", "report")
    g = lambda k: np.array([r[k] for r in rows])                     # noqa: E731
    real = g("real")
    rnd = g("random")
    rndm = g("random_m")
    print(f"\n=== Q1b  IS THE HARM IN A DETECTABLE SUBSET?   n = {len(rows)}, q = {Q} ===")
    print(f"  reproduction gate: real = {real.mean():.3f}   resid RMS "
          f"{g('resid_rms').mean():.3f}")
    print(f"\n  {'criterion':<14}{'RMSD raw':>10}{'residRMS':>10}{'vs random':>26}"
          f"   |  {'RMSD matched':>13}{'vs random_m':>26}")
    tab = {}
    for cr in CRIT:
        v, vm = g(cr), g(cr + "_m")
        m1, l1, h1 = L.boot(v - rnd, rng)
        m2, l2, h2 = L.boot(vm - rndm, rng)
        d1, d2 = v - rnd, vm - rndm
        print(f"  {cr:<14}{v.mean():>10.3f}{g(cr+'_rms').mean():>10.3f}"
              f"   {m1:+.3f} [{l1:+.3f},{h1:+.3f}] {int((d1<0).sum()):>3}W/{int((d1>0).sum())}L"
              f"   |  {vm.mean():>13.3f}"
              f"   {m2:+.3f} [{l2:+.3f},{h2:+.3f}] {int((d2<0).sum()):>3}W/{int((d2>0).sum())}L")
        tab[cr] = {"raw": float(v.mean()), "matched": float(vm.mean()),
                   "resid_rms": float(g(cr + "_rms").mean()),
                   "vs_random_raw": {"diff": m1, "ci": [l1, h1]},
                   "vs_random_matched": {"diff": m2, "ci": [l2, h2]}}
    print("\n  ORACLE_absr is the ceiling of perfect detection and is NEVER a result.")
    json.dump(tab, open(os.path.join(L.RESULTS, "a_subset_report.json"), "w"), indent=1,
              default=float)
    return tab


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
