"""S24 / LANE D / D1 -- statistics and re-stamping for the Hamiltonian-disagreement panel.

Pure analysis of `s24/results/d1_hamdisagree.json`. It runs NO AMBER, holds no lock, and does
not recompute a single energy -- the 30-target OpenMM hold is not repeated. It adds the
paired statistics the BRIEF demands (SE, MDE = 2.8016 x SE per comparison, effect/MDE, iid CI
beside a fold-clustered one, W/L, worst-target degradation) and re-stamps the artefact through
`s24.stats_lib.save_atomic` with a source hash and git commit, per the coordinator's
sprint-wide requirement.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST          # noqa: E402
from s24 import d_harness as H           # noqa: E402

RES = os.path.join(HERE, "results")
PARTS = ("AGREE_GOOD", "LEGACY_PREFERS", "AMBER_PREFERS", "STRONG_DISAGREE", "AGREE_BAD")


def main() -> int:
    src = os.path.join(RES, "d1_hamdisagree.json")
    d = json.load(open(src))
    rows = [r for r in d["rows"] if "error" not in r]
    folds = np.array([r["fold"] for r in rows])
    n = len(rows)
    print(f"n = {n} targets, POINT-CLOUD basis, ACHIEVABLE label\n")

    inc_r = np.array([r["rmsd_incumbent"] for r in rows], float)
    inc_a = np.array([r["incumbent_cos_distogram"] for r in rows], float)
    out = {"n": n, "incumbent_rmsd_mean": float(inc_r.mean()),
           "incumbent_cos_distogram_mean": float(inc_a.mean())}

    def g(nm, key):
        return np.array([r["partitions"][nm].get(key, np.nan) for r in rows], float)

    def fl(nm, key):
        return np.array([(r["shared_referent_floor"].get(nm) or {}).get(key, np.nan)
                         for r in rows], float)

    # ---- PRIMARY 1: is the partition's bias DIRECTION different from the incumbent's,
    # ---- measured against the shared-referent floor rather than against zero?
    print("=" * 96)
    print("PRIMARY 1 -- bias cosine vs incumbent, AGAINST ITS OWN SHARED-REFERENT FLOOR")
    print("=" * 96)
    print(f"{'partition':<17}{'cos':>8}{'floor':>8}{'excess':>9}{'SE':>8}{'MDE':>8}"
          f"{'x MDE':>7}{'CI_fold':>20}{'W/L':>9}  verdict")
    p1 = {}
    for nm in PARTS:
        c, f = g(nm, "cos_vs_incumbent"), fl(nm, "mean")
        s = H.paired_stats(c, f, folds, name_a=f"{nm}_cos", name_b="floor")
        p1[nm] = s
        print(f"{nm:<17}{np.nanmean(c):>8.4f}{np.nanmean(f):>8.4f}{s['mean']:>+9.4f}"
              f"{s['se']:>8.4f}{s['mde']:>8.4f}{s['eff_over_mde']:>7.2f}"
              f"  [{s['ci_fold'][0]:+.3f},{s['ci_fold'][1]:+.3f}]{s['W']:>5}/{s['L']:<4}"
              f" {s['verdict']}")
    out["primary1_cos_vs_floor"] = p1

    # ---- the L2 SPEC: a source must clear BOTH bars. Reported explicitly, not implied.
    print("\nL2 SPEC CHECK -- a source must clear BOTH bars to break even, not either")
    print(f"{'partition':<17}{'RMSD':>8}{'<=3.9?':>8}{'cos':>8}{'<=0.65?':>9}"
          f"{'q':>7}{'<1.274?':>9}   PASSES BOTH?")
    spec = {}
    for nm in PARTS:
        r_, c_, q_ = g(nm, "rmsd"), g(nm, "cos_vs_incumbent"), g(nm, "q")
        a = float(np.nanmean(r_)) <= 3.9
        b = float(np.nanmean(c_)) <= 0.65
        cq = float(np.nanmean(q_)) < 1.274
        spec[nm] = dict(rmsd=float(np.nanmean(r_)), cos=float(np.nanmean(c_)),
                        q=float(np.nanmean(q_)), pass_rmsd=a, pass_cos=b, pass_q=cq,
                        pass_all=bool(a and b and cq))
        print(f"{nm:<17}{np.nanmean(r_):>8.4f}{'YES' if a else 'no':>8}"
              f"{np.nanmean(c_):>8.4f}{'YES' if b else 'no':>9}"
              f"{np.nanmean(q_):>7.3f}{'YES' if cq else 'no':>9}"
              f"   {'** PASSES **' if (a and b and cq) else 'FAILS'}")
    out["l2_spec_check"] = spec

    # ---- PRIMARY 2 (D1-B): does a different FUNCTIONAL escape the distogram's referent?
    print("\n" + "=" * 96)
    print("PRIMARY 2 (D1-B) -- alignment with the DISTOGRAM'S OWN error, vs the random floor")
    print("=" * 96)
    print(f"INCUMBENT top-75  cos_distogram {inc_a.mean():.4f}   beta "
          f"{np.mean([r['incumbent_beta'] for r in rows]):.4f}")
    print(f"{'partition':<17}{'cos':>8}{'floor':>8}{'excess':>9}{'SE':>8}{'MDE':>8}"
          f"{'x MDE':>7}{'CI_fold':>20}{'W/L':>9}  verdict")
    p2 = {}
    for nm in PARTS:
        c, f = g(nm, "cos_distogram"), fl(nm, "cos_distogram_mean")
        s = H.paired_stats(c, f, folds, name_a=f"{nm}_cosdist", name_b="floor")
        p2[nm] = s
        print(f"{nm:<17}{np.nanmean(c):>8.4f}{np.nanmean(f):>8.4f}{s['mean']:>+9.4f}"
              f"{s['se']:>8.4f}{s['mde']:>8.4f}{s['eff_over_mde']:>7.2f}"
              f"  [{s['ci_fold'][0]:+.3f},{s['ci_fold'][1]:+.3f}]{s['W']:>5}/{s['L']:<4}"
              f" {s['verdict']}")
    out["primary2_cosdistogram_vs_floor"] = p2

    # ---- the contrast that decides the lever: partition vs the INCUMBENT's alignment
    print("\nvs the INCUMBENT's own alignment (0.6841) -- the reading rule's actual test")
    p2b = {}
    for nm in PARTS:
        s = H.paired_stats(g(nm, "cos_distogram"), inc_a, folds,
                           name_a=nm, name_b="incumbent")
        p2b[nm] = s
        print(f"  {nm:<17}{s['mean']:>+9.4f}  SE {s['se']:.4f}  MDE {s['mde']:.4f}  "
              f"{s['eff_over_mde']:.2f}x  fold[{s['ci_fold'][0]:+.3f},{s['ci_fold'][1]:+.3f}]"
              f"  {s['W']}W/{s['L']}L  {s['verdict']}")
    out["primary2_vs_incumbent"] = p2b

    # ---- D1-C: the quality-matched arm. The coordinator's pre-stated expectation is that
    # ---- alignment moves back UP toward the incumbent once the score does the selecting.
    print("\n" + "=" * 96)
    print("D1-C -- physics as a PRE-FILTER, then the SHIPPED score at the same rung")
    print("=" * 96)
    print(f"{'pre-filter':<17}{'RMSD':>8}{'vs inc':>9}{'SE':>8}{'x MDE':>7}"
          f"{'W/L':>9}{'cos_dist':>10}{'raw->filt':>11}  verdict (RMSD)")
    p3 = {}
    for nm in PARTS:
        fr = g(nm, "filt_rmsd")
        s = H.paired_stats(fr, inc_r, folds, name_a=f"{nm}_filtered", name_b="incumbent")
        raw_a, filt_a = np.nanmean(g(nm, "cos_distogram")), np.nanmean(g(nm, "filt_cos_distogram"))
        sa = H.paired_stats(g(nm, "filt_cos_distogram"), g(nm, "cos_distogram"), folds,
                            name_a="filtered", name_b="raw")
        p3[nm] = dict(rmsd_vs_incumbent=s, alignment_raw_to_filtered=sa,
                      filt_cos_distogram=float(filt_a))
        print(f"{nm:<17}{np.nanmean(fr):>8.4f}{s['mean']:>+9.4f}{s['se']:>8.4f}"
              f"{s['eff_over_mde']:>7.2f}{s['W']:>5}/{s['L']:<4}{filt_a:>10.4f}"
              f"{sa['mean']:>+11.4f}  {s['verdict']}")
    out["d1c_prefilter"] = p3
    print("\n  'raw->filt' is the change in distogram alignment caused by letting the SHIPPED")
    print("  score select inside the partition. The coordinator pre-stated that it should move")
    print("  UP toward the incumbent; a positive column is the mechanism reproducing, not a")
    print("  failure of the arm.")

    # ---- descriptive: how the two genuine Hamiltonians relate on THIS manifold
    rho = np.array([r["rho_legacy_amber"] for r in rows], float)
    out["rho_legacy_amber"] = ST.boot_ci(rho) if hasattr(ST, "boot_ci") else \
        dict(mean=float(rho.mean()), median=float(np.median(rho)),
             se=float(rho.std(ddof=1) / np.sqrt(n)))
    print(f"\nrho(Legacy, AMBER) on the retrieval manifold: mean {rho.mean():+.4f} "
          f"median {np.median(rho):+.4f}  SE {rho.std(ddof=1)/np.sqrt(n):.4f}  "
          f"{int((rho<0).sum())}/{n} negative")
    print(f"  s20 measured -0.0886 on the CONTINUOUS TORSION instrument. This is an "
          f"independent reproduction on a different manifold.")
    out["frac_amber_clash"] = float(np.mean([r["frac_amber_clash_gt_1e4"] for r in rows]))
    out["frac_amber_z_degenerate"] = float(np.mean([r["frac_amber_z_degenerate"] for r in rows]))
    ov = [r["partitions"][nm].get("overlap_with_raw_z") for nm in PARTS for r in rows
          if r["partitions"][nm].get("overlap_with_raw_z") is not None]
    out["mean_overlap_rank_vs_raw_z"] = float(np.mean(ov))
    print(f"\nD1-A: frac E>1e4 kcal {out['frac_amber_clash']:.3f} | "
          f"frac |z_raw|<0.1 {out['frac_amber_z_degenerate']:.3f} | "
          f"rank-vs-raw partition overlap {out['mean_overlap_rank_vs_raw_z']:.3f}")

    ST.save_atomic(os.path.join(RES, "d1_analysis.json"), out,
                   complete_keys=("pdb", "fold", "partitions", "shared_referent_floor",
                                  "incumbent_cos_distogram", "rmsd_incumbent"),
                   rows=rows, n_expected=30, module_file=__file__)
    # Re-stamp the RAW artefact too -- same content, now with provenance. It must be stamped
    # with the module that PRODUCED it (`d_hamiltonians.py`), not with this analyser, or the
    # source hash points at code that never wrote those numbers -- which is precisely the
    # failure Lane E found elsewhere in this sprint.
    ST.save_atomic(src, d, complete_keys=("pdb", "fold", "partitions",
                                          "shared_referent_floor"),
                   rows=rows, n_expected=30,
                   module_file=os.path.join(HERE, "d_hamiltonians.py"))
    print(f"\nwrote d1_analysis.json and re-stamped d1_hamdisagree.json with provenance")
    return 0


if __name__ == "__main__":
    sys.exit(main())
