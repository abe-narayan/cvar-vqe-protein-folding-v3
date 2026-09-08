"""AGENT D, Sprint 19 -- Block P5: THE HEADROOM TEST for the coordinator's multimodality arm.

Pre-registered in `s19/PREREG_D.md` block P5.  Exactly `s17/refine.py`'s `refine_full`
instrument -- same start (ideal-geometry projection of the coordinate average of the shipped
top-75), same objective `sum ((d - target)/sd)^2`, same optimiser `s15/align_lib.fit`, same
leave-fold-out separation debias fitted on the FULL `I.targets()`.  **Only the target vector
changes between arms; `sd` is bit-identical across every arm** (the Sprint-18 `objceil.py`
line-163 defect was exactly a control that moved the weights too).

    raw               the deployed dhat            -- validity gate: must reproduce s17 refine_full
    med               the Bayes-L1 median of the model's own distribution (native-free)
    mode1             the highest-density mode                          (native-free)
    ORACLE_bestmode   the mode nearest the truth  -- the CEILING of all mode selection
    minofN_null       same |offsets| from the mean, random signs, minimised the same way
                                                  -- the min-of-N-matched null for the ceiling
    rand_flip         mean +- |bestmode - mean| with a random sign -- matched-magnitude null
    ORACLE_true       the true distances -- the parameterisation floor, for scale only

Run:  python -m s19.d_headroom [n_targets]
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I            # noqa: E402
from s14 import avgspace as AV             # noqa: E402
from s15 import align_lib as A             # noqa: E402
from s15 import distcal as C               # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s19 import d_modality as M            # noqa: E402

#: an immutable result directory per run.  The first launch of this module was detached with
#: `nohup ... &`, killed when its parent shell exited, and a relaunch then CONTENDED for the same
#: JSON -- two writers, one file, and a reader could see either one's partial state.  The run tag
#: makes the directory unique so that can never happen again; `D_P5_headroom` is kept as the name
#: of the canonical run.
OUT = os.path.join(HERE, "results", os.environ.get("D_P5_TAG", "D_P5_headroom"))
os.makedirs(OUT, exist_ok=True)
CENTRES, WIDTHS = M.CENTRES, M.WIDTHS

ARMS = ["raw", "med", "mode1", "ORACLE_bestmode", "minofN_null", "rand_flip", "ORACLE_true"]


def aim_points(pdb, exp, dtrue, rng):
    """All alternative aim points, from the cached distogram.  Returns dict arm -> vector."""
    z = np.load(os.path.join(ROOT, "s12", "cache", f"disto_{pdb}.npz"))
    p = np.asarray(z["prob"], float)
    risk, grid = np.asarray(z["risk"], float), np.asarray(z["grid"], float)
    med = grid[np.argmin(risk, axis=1)]                 # argmin unaffected by the w scaling
    Mask = M.dens_modes(p)
    has = Mask.any(1)
    dens = np.where(Mask, p / WIDTHS[None], -np.inf)
    mode1 = np.where(has, CENTRES[np.argmax(dens, 1)], exp)
    pos = np.where(Mask, CENTRES[None], np.nan)
    dd = np.abs(pos - dtrue[:, None])
    best = np.where(has, CENTRES[np.nanargmin(np.where(np.isnan(dd), np.inf, dd), 1)], exp)
    off = np.abs(pos - exp[:, None])
    sgn = rng.choice([-1.0, 1.0], size=off.shape)
    cand = exp[:, None] + sgn * off
    cd = np.where(np.isnan(cand), np.inf, np.abs(cand - dtrue[:, None]))
    nullbest = np.where(has, cand[np.arange(len(exp)), np.argmin(cd, 1)], exp)
    delta = np.abs(best - exp)
    flip = exp + rng.choice([-1.0, 1.0], size=len(exp)) * delta
    return {"med": med, "mode1": mode1, "ORACLE_bestmode": best,
            "minofN_null": nullbest, "rand_flip": flip}


def run(nmax=None):
    tg = I.targets()
    data = C.gather(tg)                    # full list -- the subset trap, BRIEF section 7
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    use = tg if nmax is None else tg[:nmax]
    rows, t0 = [], time.time()
    path = os.path.join(OUT, "headroom.json")
    for c, t in enumerate(use):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        d = data[pdb]
        i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
        b = deb[fold](d["sep"])                       # the additive debias, identical for all arms
        rng = SD.stable_rng(pdb, "s19_D_headroom")

        W = np.asarray(AV.top75_windows(pdb)[0], float)
        P = I.pairwise_rmsd(W)
        avg, _bm = I.coordinate_average(W, P)
        pr = I.project(np.asarray(avg, float), seq, fold)
        phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)

        tgt = {"raw": d["dhat"], "ORACLE_true": d["dtrue"] + b}
        tgt.update(aim_points(pdb, d["dhat"], d["dtrue"], rng))
        e = {"pdb": pdb, "n": n, "fold": fold,
             "avg": float(I.ca_rmsd(np.asarray(avg, float), nat)),
             "proj": float(I.ca_rmsd(np.asarray(pr["ca"], float), nat))}
        for a in ARMS:
            v = np.maximum(tgt[a] - b, 2.0)
            pf, qf, fv = A.fit(v, sd, i, j, phi0, psi0)     # sd IDENTICAL in every arm
            e[a] = float(I.ca_rmsd(I.build_ca(pf, qf), nat))
            e["obj_" + a] = float(fv)
            e["shift_" + a] = float(np.abs(v - np.maximum(d["dhat"] - b, 2.0)).mean())
        rows.append(e)
        if (c + 1) % 5 == 0 or c + 1 == len(use):
            json.dump({"rows": rows, "complete": False}, open(path, "w"))
            print(f"  {c+1}/{len(use)}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"rows": rows, "complete": len(rows) == len(use), "n_expected": len(use)},
              open(path, "w"))
    if len(rows) == len(use) and nmax is None:
        open(os.path.join(OUT, "COMPLETE"), "w").write(f"n={len(rows)}\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(OUT, "headroom.json")))["rows"]
    rng = SD.stable_rng("D", "headreport")
    g = lambda k: np.asarray([r[k] for r in rows], float)     # noqa: E731

    def boot(dif, B=4000):
        k = len(dif)
        m = dif[rng.integers(0, k, size=(B, k))].mean(1)
        return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))

    raw = g("raw")
    print(f"\nn = {len(rows)}.  Start = projection of the coordinate average; arms differ ONLY "
          f"in the target vector.\n")
    print(f"  {'arm':<18}{'RMSD':>8}{'median':>9}{'vs raw [95% CI]':>26}{'W/L':>9}"
          f"{'objective':>12}{'|shift|':>9}")
    print(f"  {'avg (start)':<18}{g('avg').mean():>8.3f}{np.median(g('avg')):>9.3f}")
    print(f"  {'proj':<18}{g('proj').mean():>8.3f}{np.median(g('proj')):>9.3f}")
    for a in ARMS:
        v = g(a)
        mu, lo, hi = boot(v - raw)
        print(f"  {a:<18}{v.mean():>8.3f}{np.median(v):>9.3f}   "
              f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]"
              f"{int((v < raw).sum()):>5}/{int((v > raw).sum())}"
              f"{g('obj_' + a).mean():>12.1f}{g('shift_' + a).mean():>9.3f}")
    # the decisive contrast
    mu, lo, hi = boot(g("ORACLE_bestmode") - g("minofN_null"))
    print(f"\n  ORACLE_bestmode - minofN_null (the mode INFORMATION): "
          f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]")
    # validity gate
    try:
        s17 = {r["pdb"]: r["refine_full"] for r in
               json.load(open(os.path.join(ROOT, "s17", "results", "refine.json")))["rows"]}
        dd = np.abs(np.asarray([r["raw"] - s17[r["pdb"]] for r in rows if r["pdb"] in s17]))
        print(f"  VALIDITY GATE vs s17 refine_full: max|d| = {dd.max():.4f}, "
              f"n within 0.01 = {int((dd < 0.01).sum())}/{len(dd)}")
    except Exception as ex:                                   # noqa: BLE001
        print("  validity gate unavailable:", ex)


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else None)
