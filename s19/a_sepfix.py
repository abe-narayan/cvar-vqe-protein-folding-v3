"""SPRINT 19, AGENT A -- IS THE HARMFUL COHERENT MODE ESTIMABLE NATIVE-FREE?

THE LAST LIVE QUESTION IN THIS LANE.  `a_struct` (n=126) decomposed the residual's coherent
harm into named modes and only ONE owns it, magnitude-matched, against the 1.000 A Sprint-18 gap:

    per-target OFFSET        -0.174 A   (removing it HURTS)
    per-target STRETCH       -0.069 A   (removing it HURTS)
    per-target SEP PROFILE   +0.525 A   <-- 52.5% of the gap
    per-residue ADDITIVE     +0.073 A   (NOT MEASURED, CI spans zero)

So the harmful coherent mode is a **per-target separation profile**: the model gets this
target's distance-versus-separation curve wrong in a target-specific way, five numbers per
target.  Removing it with ORACLE knowledge takes 3.610 -> 3.085 at IDENTICAL residual RMS.

Five numbers per target is small enough that a native-free estimator might exist.  This module
tests whether one does.  For each separation bin b, the ORACLE correction is
`mean_b(dhat) - mean_b(dtrue)`; a native-free arm replaces `mean_b(dtrue)` with an estimate:

  pool_sep      the shipped top-75 retrieval windows' own mean distance in bin b
  poolmed_sep   the pool medoid's
  helix_sep     a constant ideal alpha-helix's        <- ZERO-INFORMATION (brief section 8)
  global_sep    the leave-fold-out sequence-blind separation prior's  <- ZERO-INFORMATION
  rand_sep      per-bin shifts drawn to match `pool_sep`'s shift magnitudes  <- MATCHED-RANDOM

THE TRAP THIS ARM MUST CLEAR, stated before the run.  Moving `dhat` toward the retrieval pool
moves it toward the object whose coordinate average IS the 3.048 A incumbent, so an arm can
"improve" purely by regressing to what the pipeline already does.  Two references are therefore
reported beside every arm: the coordinate average (3.048) and the FULL replacement of `dhat` by
the pool's mean distances (`poolmean` in `a_source`, 3.487).  **A five-parameter profile match
only means something if it beats BOTH, and beats its own zero-information and matched-random
controls.**  Anything between 3.487 and 3.610 is interpolation and is reported as such.

Every arm is NATIVE-FREE except `ORACLE_no_sep`, which is the ceiling and is never a result.

Run:  python -m s19.a_sepfix
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
from s14 import avgspace as AV               # noqa: E402
from s15 import seed as SD                   # noqa: E402

OUT = os.path.join(L.RESULTS, "a_sepfix.json")
NB = 5
ARMS = ["real", "ORACLE_no_sep", "pool_sep", "poolmed_sep", "helix_sep", "global_sep",
        "rand_sep", "pool_sep_half", "poolfull"]


def binof(sep):
    return np.clip(np.searchsorted([2, 4, 6, 8, 11, 999], sep, side="right") - 1, 0, NB - 1)


def bin_means(v, b, w):
    """Weighted mean of `v` in each separation bin (empty bins -> nan)."""
    out = np.full(NB, np.nan)
    for k in range(NB):
        m = b == k
        if m.any():
            out[k] = float((w[m] * v[m]).sum() / w[m].sum())
    return out


def shift_to(dhat, b, w, ref):
    """Move each bin of `dhat` so its weighted mean equals `ref`'s. Empty bins untouched."""
    cur = bin_means(dhat, b, w)
    delta = np.where(np.isfinite(cur) & np.isfinite(ref), cur - ref, 0.0)
    return dhat - delta[b], delta


def run(tg=None):
    tg = tg if tg is not None else I.targets()
    data, pdbs, folds = L.gather_all(tg)

    # the leave-fold-out sequence-blind separation prior (ZERO-INFORMATION)
    sp = {}
    for f in sorted(set(folds.tolist())):
        tr = [p for p in pdbs if data[p]["fold"] != f]
        S = np.concatenate([data[p]["sep"] for p in tr])
        T = np.concatenate([data[p]["dtrue"] for p in tr])
        bb = binof(S)
        sp[f] = np.array([T[bb == k].mean() if (bb == k).any() else np.nan
                          for k in range(NB)])

    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        d = data[pdb]
        i, j, sep, sd, nat, n = d["i"], d["j"], d["sep"], d["sd"], d["nat"], d["n"]
        dhat, dtrue = d["dhat"], d["dtrue"]
        w = 1.0 / sd ** 2
        b = binof(sep)
        rng = SD.stable_rng(pdb, "s19A_sepfix")
        phi0, psi0, avg = F.start(pdb, d["seq"], d["fold"])
        rms = float(np.sqrt(((dhat - dtrue) ** 2).mean()))

        W = np.asarray(AV.top75_windows(pdb)[0], float)
        DW = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=-1)
        d_pool = DW.mean(0)
        d_med = DW[I.medoid(I.pairwise_rmsd(W))]
        ca_h = I.build_ca(np.full(n, np.deg2rad(-63.0)), np.full(n, np.deg2rad(-42.0)))
        d_hel = np.linalg.norm(ca_h[i] - ca_h[j], axis=1)

        refs = {"pool_sep": bin_means(d_pool, b, w),
                "poolmed_sep": bin_means(d_med, b, w),
                "helix_sep": bin_means(d_hel, b, w),
                "global_sep": sp[d["fold"]],
                "ORACLE_no_sep": bin_means(dtrue, b, w)}
        fields = {"real": dhat, "poolfull": d_pool}
        deltas = {}
        for k, ref in refs.items():
            fields[k], deltas[k] = shift_to(dhat, b, w, ref)
        # MATCHED-RANDOM: the same shift magnitudes, signs and bin assignment randomised
        dl = deltas["pool_sep"]
        rp = rng.permutation(np.abs(dl)) * rng.choice([-1.0, 1.0], size=NB)
        fields["rand_sep"] = dhat - rp[b]
        fields["pool_sep_half"] = dhat - 0.5 * deltas["pool_sep"][b]
        # the ORACLE arm is reported magnitude-matched too, as in a_struct
        v = fields["ORACLE_no_sep"] - dtrue
        fields["ORACLE_no_sep"] = dtrue + v * (rms / max(float(np.sqrt((v ** 2).mean())), 1e-9))

        e = {"pdb": pdb, "n": n, "fold": d["fold"],
             "avg": float(I.ca_rmsd(avg, nat)), "resid_rms": rms,
             "pool_shift_rms": float(np.sqrt((dl ** 2).mean()))}
        for a in ARMS:
            fld = np.maximum(fields[a], 2.0)
            e[a] = F.fit_rmsd(fld, sd, i, j, phi0, psi0, nat)[0]
            e[a + "_rms"] = float(np.sqrt(((fld - dtrue) ** 2).mean()))
            e[a + "_mae"] = float(np.abs(fld - dtrue).mean())
        rows.append(e)
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"),
                      default=float)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(L.RESULTS, "a_sepfix.COMPLETE"), "w").write("ok\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rng = SD.stable_rng("s19A", "sepfix", "report")
    g = lambda k: np.array([r[k] for r in rows], float)               # noqa: E731
    real = g("real")
    folds = g("fold")
    print(f"\n=== IS THE HARMFUL MODE ESTIMABLE NATIVE-FREE?   n = {len(rows)} ===")
    print(f"  coordinate average (the incumbent)     {g('avg').mean():.3f}")
    print(f"  reproduction gate: real                {real.mean():.3f}   (objceil a0.0 3.610)")
    print(f"  mean per-bin shift asked for by pool   {g('pool_shift_rms').mean():.3f} A\n")
    print(f"  {'arm':<16}{'RMSD':>8}{'med':>8}{'MAE*':>7}{'residRMS':>10}"
          f"   vs real                              folds")
    tab = {}
    for a in ARMS:
        v = g(a)
        d = v - real
        m, lo, hi = L.boot(d, rng)
        sg = sum(1 for f in np.unique(folds)
                 if m != 0 and np.sign(d[folds == f].mean()) == np.sign(m))
        print(f"  {a:<16}{v.mean():>8.3f}{np.median(v):>8.3f}{g(a+'_mae').mean():>7.3f}"
              f"{g(a+'_rms').mean():>10.3f}   {m:+.3f} [{lo:+.3f},{hi:+.3f}]"
              f"  {int((d<0).sum()):>3}W/{int((d>0).sum())}L  {sg}/5")
        tab[a] = {"rmsd": float(v.mean()), "median": float(np.median(v)),
                  "mae": float(g(a + "_mae").mean()),
                  "resid_rms": float(g(a + "_rms").mean()),
                  "vs_real": {"diff": m, "ci": [lo, hi], "W": int((d < 0).sum()),
                              "L": int((d > 0).sum()), "folds": sg}}
    print("  * MAE is an ORACLE DIAGNOSTIC and is not an outcome.")
    print("\n  --- the controls that decide it ---")
    for a, b in (("pool_sep", "rand_sep"), ("pool_sep", "helix_sep"),
                 ("pool_sep", "global_sep"), ("pool_sep", "poolfull")):
        s, rec = L.report_pair(f"{a} - {b}", g(a), g(b), rng, folds)
        print("  " + s)
        tab[f"_{a}_vs_{b}"] = rec
    json.dump(tab, open(os.path.join(L.RESULTS, "a_sepfix_report.json"), "w"), indent=1,
              default=float)
    return tab


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
