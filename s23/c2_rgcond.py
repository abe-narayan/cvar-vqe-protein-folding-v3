"""s23/c2_rgcond.py -- H_C2: rg_z-CONDITIONED AGGREGATION WIDTH.

Pre-registered in `s23/PREREG_C.md` before this file read a single RMSD. See that file for the
full six-axis fork list and the falsifier. Summary of the design actually run:

  * classical m-ladder from `s21/results/poolgap.json` (avg_500/150/75/20/5/1), bit-verified
    avg_75 mean = 3.0483 = the incumbent.
  * rg_z (and rg_gap, secondary) from `s23/qc_lib.build_rg_table` -- the full 126-target
    extension of `s21/rgsign.py`'s recipe, verified bit-identical on the 75-target overlap.
  * ONE fitted scalar per fold (a threshold `tau` on rg_z), TWO pre-declared rungs never fit
    (m_lo=75 for rg_z<tau, m_hi=150 for rg_z>=tau), nested nested 5-fold CV (leave-one-pinned-
    fold-out; tau fit on the other 4, applied once to the held-out fold).
  * secondary/robustness, reported regardless of sign: rg_gap in place of rg_z; and the
    REVERSED direction (m_hi for LOW rg_z) as an adversarial cherry-picking control.
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

from s23 import qc_lib as QL   # noqa: E402

POOLGAP = os.path.join(ROOT, "s21", "results", "poolgap.json")
M_LO, M_HI = 75, 150
RUNG_KEY = {75: "avg_75", 150: "avg_150", 500: "avg_500", 20: "avg_20", 5: "avg_5", 1: "avg_1"}


def load_frame():
    pg = json.load(open(POOLGAP))
    assert pg["complete"] and pg["n_expected"] == 126
    ladder = {r["pdb"]: r for r in pg["rows"]}
    rg = {r["pdb"]: r for r in QL.build_rg_table()}
    common = sorted(set(ladder) & set(rg))
    assert len(common) == 126, f"expected 126 overlapping targets, got {len(common)}"
    pdbs = np.array(common)
    fold = np.array([ladder[p]["fold"] for p in common], int)
    fixed75 = np.array([ladder[p][RUNG_KEY[75]] for p in common], float)
    rung150 = np.array([ladder[p][RUNG_KEY[150]] for p in common], float)
    rgz = np.array([rg[p]["rg_z"] for p in common], float)
    rggap = np.array([rg[p]["rg_gap"] for p in common], float)
    # sanity: fixed75 must reproduce the incumbent mean exactly
    assert abs(float(fixed75.mean()) - 3.0483380938795324) < 1e-9
    return dict(pdb=pdbs, fold=fold, m_lo_rmsd=fixed75, m_hi_rmsd=rung150,
                rg_z=rgz, rg_gap=rggap)


def _route(feat_train, m_lo_train, m_hi_train, feat_eval, m_lo_eval, m_hi_eval,
          reverse: bool = False):
    """Grid-search tau on TRAIN, minimising train mean RMSD; apply on EVAL. Returns
    (routed_eval, tau, train_mean_at_tau)."""
    cands = np.unique(feat_train)
    # midpoints between consecutive sorted unique values, plus the extremes (all-lo / all-hi)
    if len(cands) > 1:
        mids = (cands[:-1] + cands[1:]) / 2.0
    else:
        mids = np.array([])
    taus = np.concatenate([[cands.min() - 1.0], mids, [cands.max() + 1.0]])
    best_tau, best_mean = None, np.inf
    for tau in taus:
        if not reverse:
            lo_mask = feat_train < tau     # low feature -> m_lo (declared direction)
        else:
            lo_mask = feat_train >= tau    # reversed adversarial control
        pred = np.where(lo_mask, m_lo_train, m_hi_train)
        m = float(pred.mean())
        if m < best_mean:
            best_mean, best_tau = m, float(tau)
    if not reverse:
        lo_eval = feat_eval < best_tau
    else:
        lo_eval = feat_eval >= best_tau
    routed_eval = np.where(lo_eval, m_lo_eval, m_hi_eval)
    return routed_eval, best_tau, best_mean


def nested_cv(feat: np.ndarray, fold: np.ndarray, m_lo_rmsd: np.ndarray,
             m_hi_rmsd: np.ndarray, reverse: bool = False) -> dict:
    n = len(feat)
    routed = np.empty(n, float)
    taus = {}
    for f in np.unique(fold):
        train = fold != f
        held = fold == f
        r_eval, tau, train_mean = _route(feat[train], m_lo_rmsd[train], m_hi_rmsd[train],
                                         feat[held], m_lo_rmsd[held], m_hi_rmsd[held],
                                         reverse=reverse)
        routed[held] = r_eval
        taus[int(f)] = tau
    return {"routed": routed, "taus": taus}


def run():
    fr = load_frame()
    out = {"n": 126, "m_lo": M_LO, "m_hi": M_HI, "primary": {}, "secondary": {}}

    # ---------------- PRIMARY: rg_z, declared direction ----------------
    cv = nested_cv(fr["rg_z"], fr["fold"], fr["m_lo_rmsd"], fr["m_hi_rmsd"], reverse=False)
    ps = QL.paired_stats(cv["routed"], fr["m_lo_rmsd"], fr["fold"],
                         label="routed(rg_z) - fixed(m=75)")
    out["primary"] = {
        "routed_mean": float(cv["routed"].mean()), "fixed_mean": float(fr["m_lo_rmsd"].mean()),
        "taus_per_fold": cv["taus"], "stats": ps, "verdict": QL.verdict(ps),
        "n_routed_to_hi": int((cv["routed"] == fr["m_hi_rmsd"]).sum()),
    }

    # ---------------- SECONDARY 1: rg_gap in place of rg_z ----------------
    cvg = nested_cv(fr["rg_gap"], fr["fold"], fr["m_lo_rmsd"], fr["m_hi_rmsd"], reverse=False)
    psg = QL.paired_stats(cvg["routed"], fr["m_lo_rmsd"], fr["fold"],
                          label="routed(rg_gap) - fixed(m=75)")
    out["secondary"]["rg_gap"] = {
        "routed_mean": float(cvg["routed"].mean()), "taus_per_fold": cvg["taus"],
        "stats": psg, "verdict": QL.verdict(psg),
    }

    # ---------------- SECONDARY 2: REVERSED direction, adversarial control ----------------
    cvr = nested_cv(fr["rg_z"], fr["fold"], fr["m_lo_rmsd"], fr["m_hi_rmsd"], reverse=True)
    psr = QL.paired_stats(cvr["routed"], fr["m_lo_rmsd"], fr["fold"],
                          label="routed(rg_z,REVERSED) - fixed(m=75)")
    out["secondary"]["reversed_direction"] = {
        "routed_mean": float(cvr["routed"].mean()), "taus_per_fold": cvr["taus"],
        "stats": psr, "verdict": QL.verdict(psr),
    }

    # ---------------- context: the m=150 fixed arm and the oracle-in-sample rung choice ------
    ps150 = QL.paired_stats(fr["m_hi_rmsd"], fr["m_lo_rmsd"], fr["fold"],
                            label="fixed(m=150) - fixed(m=75)")
    out["context"] = {"fixed_m150_vs_m75": ps150}

    out["complete"] = True
    QL.save_json("c2_rgcond.json", out)
    return out


def report(out=None):
    if out is None:
        out = json.load(open(os.path.join(QL.RESULTS, "c2_rgcond.json")))
    print("=== H_C2 -- rg_z-conditioned aggregation width ===")
    print(f"fixed m=75 (incumbent):   {out['primary']['fixed_mean']:.4f}")
    print(f"fixed m=150 (context):    {out['context']['fixed_m150_vs_m75']['mean_diff'] + out['primary']['fixed_mean']:.4f}")
    p = out["primary"]
    s = p["stats"]
    print(f"\nPRIMARY routed(rg_z):     {p['routed_mean']:.4f}  "
          f"({p['n_routed_to_hi']}/126 routed to m=150)")
    print(f"  diff vs fixed: {s['mean_diff']:+.4f}  SE {s['se']:.4f}  MDE {s['mde']:.4f}  "
          f"|eff|/MDE {s['abs_over_mde']:.2f}")
    print(f"  CI iid  {s['ci95_iid']}   CI fold {s['ci95_fold']}   W/L {s['n_better']}/{s['n_worse']}")
    print(f"  VERDICT: {p['verdict']}")
    print(f"  per-fold tau: {p['taus_per_fold']}")
    for name, d in out["secondary"].items():
        s2 = d["stats"]
        print(f"\nSECONDARY [{name}]: routed_mean={d['routed_mean']:.4f}  "
              f"diff={s2['mean_diff']:+.4f}  |eff|/MDE={s2['abs_over_mde']:.2f}  "
              f"CI_fold={s2['ci95_fold']}  VERDICT={d['verdict']}")


if __name__ == "__main__":
    o = run()
    report(o)
