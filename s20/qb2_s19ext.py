"""SPRINT 20 / WORKSTREAM B -- RECOVERY of the SECOND unanalysed Sprint-19 quantum artefact.

    python -m s20.qb2_s19ext

`s19/results/qb_ext.json` (40/40 targets, COMPLETE on disk) holds the Sprint-19 lane's
"go deeper than an alpha sweep" panel -- warm starts, adaptive alpha, CVaR+uniform and
CVaR+annealing mixtures, a deeper ansatz, alpha = 0.02 and 0.50, and **four independent ansatz
seeds of both the trained and the untrained circuit**.  It was never analysed.  Workstream D's
closure (`s20/results/D_QB_CLOSE`) analysed `qb_main.json`; this module analyses `qb_ext.json`,
which D did not read, and it is scored against `qb_main`'s own arms ON THE SAME 40 TARGETS so
that nothing is compared across different target sets.

Nothing here is re-run.  Pure analysis of a prior lane's completed artefact.
"""
from __future__ import annotations

import json
import os

import numpy as np

from s20 import qb2_lib as L

EXT = os.path.join(L.ROOT, "s19", "results", "qb_ext.json")
MAIN = os.path.join(L.ROOT, "s19", "results", "qb_main.json")
XARMS = ["x_warm", "x_adapt", "x_mixunif50", "x_mixanneal", "x_mps3f", "x_a0.02", "x_a0.50"]
REF = ["pool500", "c_marg", "c_metroH", "c_helix", "q_a1.00", "q_a0.25", "q_untrained"]


def report():
    with open(EXT) as fh:
        e = json.load(fh)
    with open(MAIN) as fh:
        m = json.load(fh)
    ks = sorted(k for k in e if isinstance(e[k], dict) and "arms" in e[k])
    ks = [k for k in ks if k in m]
    folds = [e[p]["fold"] for p in ks]
    out = {"n": len(ks), "targets": ks}

    def g(src, a, key="realised_ORACLE"):
        return np.array([src[p]["arms"][a][key] for p in ks])

    S = np.array([[e[p]["arms"][f"s_q0.25_{s}"]["realised_ORACLE"] for s in range(4)]
                  for p in ks])
    U = np.array([[e[p]["arms"][f"s_untrained_{s}"]["realised_ORACLE"] for s in range(4)]
                  for p in ks])
    print(f"=== RECOVERED: s19/results/qb_ext.json, n = {len(ks)} targets "
          f"(a declared prefix of the 126), realised = BUILT chain ===")
    print("\n--- SEED SENSITIVITY, the number this lane most needed and never had ---")
    out["seed"] = {
        "trained_mean": float(S.mean()), "trained_within_target_sd": float(S.std(1).mean()),
        "trained_best_of_4": float(S.min(1).mean()),
        "untrained_mean": float(U.mean()), "untrained_within_target_sd": float(U.std(1).mean()),
        "untrained_best_of_4": float(U.min(1).mean())}
    print(f"  trained   4 seeds: mean {S.mean():.4f}  WITHIN-TARGET sd over seeds "
          f"{S.std(1).mean():.4f}  best-of-4 {S.min(1).mean():.4f}")
    print(f"  untrained 4 seeds: mean {U.mean():.4f}  WITHIN-TARGET sd over seeds "
          f"{U.std(1).mean():.4f}  best-of-4 {U.min(1).mean():.4f}")
    print(f"  >>> the 126-target instrument's MDE is 0.084 A; the ANSATZ SEED alone moves a "
          f"target by {S.std(1).mean():.3f} A.")
    t = L.paired_ci(S.mean(1), U.mean(1), folds=folds)
    out["trained_vs_untrained_4seed"] = t
    print(f"  trained - untrained (4-seed means, the MANDATORY control): {t['mean']:+.4f} "
          f"[{t['ci95'][0]:+.4f}, {t['ci95'][1]:+.4f}] W/L {t['W']}/{t['L']} "
          f"folds {t.get('folds_same_sign')}/5 sig={int(t['sig'])}")

    print("\n--- the extension arms, and the reference arms ON THE SAME 40 TARGETS ---")
    print(f"  {'arm':<20}{'realised':>10}{'cloud':>10}{'gen':>10}{'sel':>10}")
    out["means"] = {}
    rows = [("q0.25 (4-seed mean)", S.mean(1), None)]
    for a in XARMS:
        rows.append((a, g(e, a), ("ext", a)))
    for a in REF:
        rows.append((a, g(m, a), ("main", a)))
    for lab, v, src in rows:
        if src is None:
            print(f"  {lab:<20}{v.mean():>10.4f}{'-':>10}{'-':>10}{'-':>10}")
            out["means"][lab] = {"realised": float(v.mean())}
            continue
        s = e if src[0] == "ext" else m
        out["means"][lab] = {
            "realised": float(v.mean()),
            "cloud": float(g(s, src[1], "avg_ORACLE").mean()),
            "gen": float(g(s, src[1], "gen_best_ORACLE").mean()),
            "sel": float(g(s, src[1], "sel_best_ORACLE").mean())}
        q = out["means"][lab]
        print(f"  {lab:<20}{q['realised']:>10.4f}{q['cloud']:>10.4f}{q['gen']:>10.4f}"
              f"{q['sel']:>10.4f}")

    print("\n--- every extension arm vs the 4-seed trained baseline (negative = the arm wins) ---")
    out["vs_base"] = {}
    for a in XARMS:
        t = L.paired_ci(g(e, a), S.mean(1), folds=folds)
        out["vs_base"][a] = t
        print(f"  {a:<14}{t['mean']:+8.4f} [{t['ci95'][0]:+.4f},{t['ci95'][1]:+.4f}] "
              f"W/L {t['W']}/{t['L']} folds {t.get('folds_same_sign')}/5 sig={int(t['sig'])}")

    print("\n--- and against the two arms that actually matter (SAME 40 targets) ---")
    out["vs_ref"] = {}
    for a in XARMS:
        for b in ("pool500", "c_marg"):
            t = L.paired_ci(g(e, a), g(m, b), folds=folds)
            out["vs_ref"].setdefault(a, {})[b] = t
            print(f"  {a:<14} vs {b:<9}{t['mean']:+8.4f} [{t['ci95'][0]:+.4f},"
                  f"{t['ci95'][1]:+.4f}] W/L {t['W']}/{t['L']} "
                  f"folds {t.get('folds_same_sign')}/5 sig={int(t['sig'])}")
    L.write("qb2_s19ext", out, complete=True)
    return out


if __name__ == "__main__":
    report()
