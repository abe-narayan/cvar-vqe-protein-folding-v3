"""s19/agentC_kv2.py -- Q1, THE DECLARED EXTENSION: is the diversity channel CAUSAL?

**This module was written AFTER seeing `agentC_kv`'s pre-registered primary, and it says so.**
It is not a rescue: `PREREG_C.md` section 3.4's primary stands as written and as reported, and
every arm here can only make the mechanism claim HARDER to sustain, never easier.

WHAT THE PRIMARY LEFT UNSETTLED.  Section 3.4 measured that a Legacy gate collapses the ensemble's
ambiguity `D` relative to a matched-random gate, and that the diversity channel accounts for more
than all of the damage.  But `D` and "being a physics score" are confounded in those arms: every
arm that lost diversity was also an arm that ordered on an energy.  Two arms break the confound,
and both are NATIVE-FREE (`D` needs no native):

    rand_lowD    of R = 200 matched-random subsets of the SAME size, the one with the LOWEST D
    rand_highD   ... the one with the HIGHEST D

**Neither contains one bit of physics, of sequence, or of any score.**  If the damage is carried by
diversity collapse, `rand_lowD` must hurt like a physics gate and `rand_highD` must not -- and the
gap between them must be produced by nothing but `D`.  If instead the damage needs a *score* to
appear, the two random-extreme arms will be indistinguishable and the mechanism claim fails.

Two further score-free arms bracket the axis:

    fps          pure greedy farthest-point on Ca-RMSD: the MAXIMUM-diversity gate, no score
    medoid       the m candidates closest to the pool medoid: the MINIMUM-diversity gate, no score

    python -m s19.agentC_kv2
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                                   # noqa: E402
from s14.avgspace import top75_windows                            # noqa: E402
from s15 import seed as SD                                        # noqa: E402
from s18 import phys_lib as PL                                    # noqa: E402
from s19 import agentC_lib as CL                                  # noqa: E402
from s19 import agentC_kv as KV                                   # noqa: E402

R_SCREEN = 200            # native-free draws screened on D


def run_target(t):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    nat = np.asarray(u["nat_ca"], float)
    K = len(W)
    m = max(2, int(round(K * (1.0 - KV.F_PRIMARY))))
    d_win = np.asarray(I.kabsch_rmsd_batch(W, nat), float)
    P = np.asarray(I.pairwise_rmsd(W), float)
    rng = SD.stable_rng(pdb, "s19C_kv2")

    #: the screen: D of R random subsets.  NATIVE-FREE -- D is a property of the candidates only.
    draws, Ds = [], []
    for _ in range(R_SCREEN):
        idx = CL.gate_random(K, m, rng)
        Y, _Yb, _C, _dv = CL.common_frame_members(W[idx], PHI[idx], PSI[idx])
        Ds.append(float(np.sum((Y - Y.mean(0)) ** 2) / (m * n)))
        draws.append(idx)
    Ds = np.asarray(Ds, float)
    lo, hi = int(np.argmin(Ds)), int(np.argmax(Ds))

    idxs = {"rand_lowD": draws[lo], "rand_highD": draws[hi],
            "fps": CL._farthest_point(P, m, int(np.argmin(P.mean(1)))),
            "medoid": np.sort(np.argsort(P[int(np.argmin(P.mean(1)))])[:m])}

    arms = {}
    for nm, idx in idxs.items():
        Y, _Yb, _C, _dv = CL.common_frame_members(W[idx], PHI[idx], PSI[idx])
        k = CL.kv_of(Y, nat, d_win=d_win[idx])
        k["div_s18"] = PL.diversity(W[idx])
        arms[nm] = k
    return {"pdb": pdb, "n": n, "fold": int(t["fold"]), "K": int(K), "m": int(m),
            "D2_screen_min": float(Ds.min()), "D2_screen_max": float(Ds.max()),
            "D2_screen_mean": float(Ds.mean()), "arms": arms}


def run(targets=None, out="agentC_kv2.json", verbose=True):
    import json
    tg = targets if targets is not None else I.targets()
    path = os.path.join(CL.RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    cfg = {"R_SCREEN": R_SCREEN, "F_PRIMARY": KV.F_PRIMARY,
           "note": "declared extension, written after the section-3.4 primary was read"}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        PL.mem_ok(0.8)
        rows.append(run_target(t))
        if verbose and len(rows) % 10 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
        KV._write(out, rows, cfg, len(tg))
    KV._write(out, rows, cfg, len(tg))
    return rows


def report():
    o = PL.read_complete(os.path.join(CL.RESULTS, "agentC_kv2.json"), need=126)
    o1 = PL.read_complete(os.path.join(CL.RESULTS, "agentC_kv.json"), need=126)
    r2 = {r["pdb"]: r for r in o["rows"]}
    rows = o1["rows"]
    pdbs = [r["pdb"] for r in rows]
    folds = np.array([r["fold"] for r in rows])

    def c1(arm, key):
        return np.array([r["arms"][arm].get(key, np.nan) for r in rows], float)

    def c2(arm, key):
        return np.array([r2[p]["arms"][arm].get(key, np.nan) for p in pdbs], float)

    print("=" * 100)
    print("Q1 EXTENSION -- SCORE-FREE DIVERSITY MANIPULATION.  Declared extension, n = %d."
          % len(rows))
    print("   Every arm below contains NO score, NO physics and NO sequence information.")
    print("=" * 100)
    print(f"\n   {'arm':<14}{'readout':>9}{'E_mem':>9}{'D':>9}{'err_cos':>9}{'shared':>8}"
          f"      vs matched-random of the same count")
    for a, col in (("rand_lowD", c2), ("rand_highD", c2), ("medoid", c2), ("fps", c2),
                   ("legacy", c1), ("legacy_clust", c1), ("helix", c1),
                   ("rand", c1), ("none", c1)):
        p = PL.paired(col(a, "readout"), c1("rand", "readout"), folds=folds, names=pdbs)
        print(f"   {a:<14}{np.nanmean(col(a,'readout')):>9.3f}{np.nanmean(col(a,'E_mem')):>9.3f}"
              f"{np.nanmean(col(a,'D')):>9.3f}{np.nanmean(col(a,'err_cos')):>9.3f}"
              f"{np.nanmean(col(a,'shared_frac')):>8.3f}      "
              + ("--" if a == "rand" else
                 f"{p['mean']:+.4f} [{p.get('ci_fold',p['ci'])[0]:+.4f},"
                 f"{p.get('ci_fold',p['ci'])[1]:+.4f}] {p['W']}W/{p['L']}L"))
    print("\n   THE DECISIVE CONTRAST (both arms are random subsets; only D differs):")
    p = PL.paired(c2("rand_lowD", "readout"), c2("rand_highD", "readout"),
                  folds=folds, names=pdbs)
    ci = p.get("ci_fold", p["ci"])
    print(f"      rand_lowD - rand_highD = {p['mean']:+.4f} [{ci[0]:+.4f}, {ci[1]:+.4f}]  "
          f"med {p['median']:+.4f}  {p['W']}W/{p['L']}L")
    print(f"      D:  lowD {np.nanmean(c2('rand_lowD','D')):.3f}   "
          f"highD {np.nanmean(c2('rand_highD','D')):.3f}   "
          f"legacy {np.nanmean(c1('legacy','D')):.3f}   rand {np.nanmean(c1('rand','D')):.3f}")
    print(f"      E_mem: lowD {np.nanmean(c2('rand_lowD','E_mem')):.3f}   "
          f"highD {np.nanmean(c2('rand_highD','E_mem')):.3f}     "
          f"(if these differ, the contrast is NOT diversity alone -- reported, not hidden)")

    print("\n   THE ONE-PARAMETER LAW, pooled over every arm and target "
          "(readout^2 = E_mem^2 - D^2 is EXACT; what is measured is how much E_mem^2 MOVES):")
    A1 = [a for a in rows[0]["arms"]]
    A2 = [a for a in o["rows"][0]["arms"]]
    dr, dd, de = [], [], []
    for a in A1:
        if a == "rand":
            continue
        dr.append(c1(a, "readout2") - c1("rand", "readout2"))
        dd.append(c1("rand", "D2") - c1(a, "D2"))
        de.append(c1(a, "E_mem2") - c1("rand", "E_mem2"))
    for a in A2:
        dr.append(c2(a, "readout2") - c1("rand", "readout2"))
        dd.append(c1("rand", "D2") - c2(a, "D2"))
        de.append(c2(a, "E_mem2") - c1("rand", "E_mem2"))
    dr = np.concatenate(dr); dd = np.concatenate(dd); de = np.concatenate(de)
    print(f"      n = {len(dr)} (arm, target) pairs over {len(A1)+len(A2)-1} arms")
    print(f"      corr(Delta readout^2,  -Delta D^2 )   = {np.corrcoef(dr, dd)[0,1]:+.4f}")
    print(f"      corr(Delta readout^2, Delta E_mem^2)  = {np.corrcoef(dr, de)[0,1]:+.4f}")
    print(f"      sd(-Delta D^2) {dd.std():.4f}   sd(Delta E_mem^2) {de.std():.4f}   "
          f"sd(Delta readout^2) {dr.std():.4f}")
    return o


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    if a.report:
        report()
    else:
        run(I.targets()[:5] if a.smoke else None,
            out="_SMOKE_agentC_kv2.json" if a.smoke else "agentC_kv2.json")
