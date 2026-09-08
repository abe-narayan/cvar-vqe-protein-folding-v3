"""SPRINT 15, ALIGN, TASK 3 -- THE HONEST TEST: is alignment a LEVER or a DESCRIPTION?

The alignment story makes one checkable prediction:

    an intervention that improves ALIGNMENT should improve emitted RMSD
    WITHOUT improving raw distance or torsion accuracy.

If every RMSD gain is accompanied by a raw-accuracy gain, alignment is a re-description of
accuracy and must be reported as such.  This module reads `s15/results/align_fit.json` and
runs four tests, all paired, all with bootstrap CIs:

  T1  per arm: the paired change in RMSD, in torsion RMS error, in emitted-distance MAE, and
      in alignment, side by side.  An arm that wins on RMSD while LOSING or TYING on both raw
      axes is direct evidence for the lever.
  T2  pooled over every (target, arm) pair: regress the RMSD change on the raw-accuracy change
      and the alignment change together, standardised, and report the partial correlation of
      the RMSD change with the alignment change CONTROLLING for raw accuracy.
  T3  robustness vs alignment: is the combination additive?  If robustness and alignment are
      the same lever the combination is sub-additive.  The interaction is
      `delta(combo) - delta(robust) - delta(align)`, with a CI.
  T4  do the two levers help the SAME targets?  Correlation of the per-target gains.

    python -m s15.align_axes
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I              # noqa: E402
from s15 import info_lib as L                # noqa: E402

PATH = os.path.join(ROOT, "s15", "results", "align_fit.json")
KEYS = ("rmsd", "tors_rms_deg", "dist_mae", "resid_medz", "align", "align_self")


def load():
    with open(PATH) as fh:
        d = json.load(fh)
    return d["rows"], d


def mat(rows, arm, key):
    return np.asarray([r["arms"][arm][key] for r in rows], float)


def main():
    rows, doc = load()
    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([r["fold"] for r in rows], int)
    arms = [a for a in sorted(rows[0]["arms"]) if a != "squared"]
    base = {k: mat(rows, "squared", k) for k in KEYS}

    print(f"TASK 3 -- BOTH AXES FOR EVERY ARM (n = {len(rows)}), paired against `squared`\n")
    print(f"{'arm':<26}{'dRMSD':>9}{'  CI95':<20}{'dTORSION(deg)':>15}{'  CI95':<20}"
          f"{'dDISTMAE':>10}{'dALIGN':>9}{'verdict':>12}")
    out = {"per_arm": {}}
    for a in arms:
        v = {k: mat(rows, a, k) for k in KEYS}
        pr = {k: I.paired(v[k], base[k], folds=folds, names=pdbs) for k in KEYS}
        dr, dt, dd, da = (pr["rmsd"], pr["tors_rms_deg"], pr["dist_mae"], pr["align"])
        win_rmsd = dr["ci95"][1] < 0
        raw_better = dt["ci95"][1] < 0 or dd["ci95"][1] < 0
        if win_rmsd and not raw_better:
            verdict = "LEVER"
        elif win_rmsd and raw_better:
            verdict = "confounded"
        elif dr["ci95"][0] > 0:
            verdict = "HURTS"
        else:
            verdict = "null"
        out["per_arm"][a] = {k: pr[k] for k in KEYS}
        out["per_arm"][a]["verdict"] = verdict
        print(f"{a:<26}{dr['mean_diff']:>+9.3f}"
              f"  [{dr['ci95'][0]:+.3f},{dr['ci95'][1]:+.3f}]  "
              f"{dt['mean_diff']:>+13.2f}"
              f"  [{dt['ci95'][0]:+.2f},{dt['ci95'][1]:+.2f}]  "
              f"{dd['mean_diff']:>+8.3f}{da['mean_diff']:>+9.3f}{verdict:>12}")

    # ------------------------------------------------------------------ T2 pooled regression
    print("\nT2 -- POOLED over every (target, arm): does the alignment change explain the "
          "RMSD change\n     once raw accuracy is controlled for?")
    keep = [a for a in arms if not a.startswith("ORACLE") and not a.endswith("_lfo")
            and "INFOLD" not in a]
    dR = np.concatenate([mat(rows, a, "rmsd") - base["rmsd"] for a in keep])
    dT = np.concatenate([mat(rows, a, "tors_rms_deg") - base["tors_rms_deg"] for a in keep])
    dD = np.concatenate([mat(rows, a, "dist_mae") - base["dist_mae"] for a in keep])
    dA = np.concatenate([mat(rows, a, "align") - base["align"] for a in keep])
    dAs = np.concatenate([mat(rows, a, "align_self") - base["align_self"] for a in keep])
    N = len(dR)

    def z(x):
        s = x.std()
        return (x - x.mean()) / s if s > 0 else x * 0

    X = np.column_stack([np.ones(N), z(dT), z(dD), z(dA)])
    beta, *_ = np.linalg.lstsq(X, z(dR), rcond=None)
    r_par = L.pearson(L.residualise(dR, dT, dD), L.residualise(dA, dT, dD))
    print(f"  n = {N} (target, arm) pairs over {len(keep)} native-free arms")
    print(f"  marginal corr(dRMSD, dTORSION) {L.pearson(dR, dT):+.3f}   "
          f"corr(dRMSD, dDISTMAE) {L.pearson(dR, dD):+.3f}   "
          f"corr(dRMSD, dALIGN) {L.pearson(dR, dA):+.3f}")
    print(f"  standardised betas: torsion {beta[1]:+.3f}  distMAE {beta[2]:+.3f}  "
          f"align {beta[3]:+.3f}")
    print(f"  PARTIAL corr(dRMSD, dALIGN | dTORSION, dDISTMAE) = {r_par:+.3f}   "
          f"(native-free-J version {L.pearson(L.residualise(dR, dT, dD), L.residualise(dAs, dT, dD)):+.3f})")
    out["pooled"] = {"n": int(N), "arms": keep,
                     "corr_dRMSD_dTORS": L.pearson(dR, dT),
                     "corr_dRMSD_dDIST": L.pearson(dR, dD),
                     "corr_dRMSD_dALIGN": L.pearson(dR, dA),
                     "beta_tors": float(beta[1]), "beta_dist": float(beta[2]),
                     "beta_align": float(beta[3]), "partial_align": float(r_par)}

    # ------------------------------------------------------------------ T3 additivity
    print("\nT3 -- ARE ROBUSTNESS AND ALIGNMENT THE SAME LEVER?  (interaction = combo minus "
          "the two parts)")
    out["additivity"] = {}
    for combo, ra, aa in (("cauchy1_x_term", "cauchy_s1.0", "term2_t0.25"),
                          ("cauchy1_x_damp", "cauchy_s1.0", "damp_k2.0"),
                          ("term_x_damp", "term2_t0.25", "damp_k2.0")):
        if combo not in rows[0]["arms"]:
            continue
        c = mat(rows, combo, "rmsd") - base["rmsd"]
        p1 = mat(rows, ra, "rmsd") - base["rmsd"]
        p2 = mat(rows, aa, "rmsd") - base["rmsd"]
        inter = c - (p1 + p2)
        pr = I.paired(c, p1 + p2, folds=folds, names=pdbs)
        out["additivity"][combo] = {"parts": [ra, aa], "d_combo": float(c.mean()),
                                    "d_part1": float(p1.mean()), "d_part2": float(p2.mean()),
                                    "interaction": pr}
        tag = ("SUB-additive (one lever)" if pr["ci95"][0] > 0 else
               "SUPER-additive" if pr["ci95"][1] < 0 else "additive (two levers)")
        print(f"  {combo:<18} combo {c.mean():+.3f} = {p1.mean():+.3f} ({ra}) + "
              f"{p2.mean():+.3f} ({aa})   interaction {pr['mean_diff']:+.3f} "
              f"[{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}]  -> {tag}")

    # ------------------------------------------------------------------ T4 same targets?
    print("\nT4 -- do the two families help the SAME targets?  (per-target gain correlation)")
    fam = {"robust(cauchy1)": "cauchy_s1.0", "terminal(t0.25)": "term2_t0.25",
           "damage(k2)": "damp_k2.0", "rank(f0.5)": "rank_f0.5"}
    fam = {k: v for k, v in fam.items() if v in rows[0]["arms"]}
    ks = list(fam)
    G = {k: mat(rows, fam[k], "rmsd") - base["rmsd"] for k in ks}
    print("      " + "".join(f"{k:>18}" for k in ks))
    out["gain_corr"] = {}
    for k in ks:
        line = f"{k:<18}"
        for k2 in ks:
            r = L.pearson(G[k], G[k2])
            out["gain_corr"][f"{k}|{k2}"] = r
            line += f"{r:>18.3f}"
        print(line)

    L.jwrite("align_axes", out)
    return out


if __name__ == "__main__":
    main()
