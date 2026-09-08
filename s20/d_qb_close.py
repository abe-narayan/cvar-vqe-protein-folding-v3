"""SPRINT 20, AGENT D -- close Sprint 19 section 5 (CVaR-VQE), NOT DELIVERED at freeze.

Pure analysis of `s19/results/qb_main.json`, which is COMPLETE on disk (126/126 targets,
17 arms, no missing cells) but was never reported.  Executes `s19/PREREG_B.md` sections 6-9
exactly as written: P1 realised, P2 eps-dominance in the (M,D) plane, P3 generation ceiling,
plus every kill rule in section 8.  No structure is rebuilt; every number here comes out of
Agent B's own artefact.

TARGET is the unit.  Every RMSD column is ORACLE post-hoc scoring.
"""
from __future__ import annotations
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I     # noqa: E402

OUT = os.path.join(HERE, "results", "D_QB_CLOSE")
os.makedirs(OUT, exist_ok=True)
MDE = 0.084

Q = ["q_a0.05", "q_a0.25", "q_a1.00", "q_anneal"]
CLASSICAL = ["c_prod0.05", "c_prod0.25", "c_chain0.05", "c_chain0.25", "c_cem0.05",
             "c_metroL", "c_metroH", "c_anneal", "c_lbfgs"]
ZEROINFO = ["c_marg", "c_helix"]


def boot(d, rng, B=4000):
    d = np.asarray(d, float)
    k = len(d)
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def fold_boot(d, folds, rng, B=2000):
    """Fold-CLUSTERED bootstrap: resample folds, then targets within each drawn fold."""
    d = np.asarray(d, float)
    folds = np.asarray(folds)
    uf = np.unique(folds)
    idx = {f: np.flatnonzero(folds == f) for f in uf}
    out = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, len(uf), len(uf))
        sel = np.concatenate([idx[uf[p]][rng.integers(0, len(idx[uf[p]]), len(idx[uf[p]]))]
                              for p in pick])
        out[b] = d[sel].mean()
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def cmp_arm(name, a, b, folds, rng):
    dif = np.asarray(a, float) - np.asarray(b, float)
    m, lo, hi = boot(dif, rng)
    flo, fhi = fold_boot(dif, folds, rng)
    sgn = sum(1 for f in np.unique(folds) if np.sign(dif[folds == f].mean()) == np.sign(m))
    return {"name": name, "mean_a": float(np.mean(a)), "mean_b": float(np.mean(b)),
            "diff": m, "ci": [lo, hi], "fold_ci": [flo, fhi],
            "median": float(np.median(dif)),
            "W": int((dif < 0).sum()), "L": int((dif > 0).sum()),
            "folds_same_sign": int(sgn), "n_folds": int(len(np.unique(folds)))}


def line(r):
    return (f"  {r['name']:<36}{r['mean_a']:>7.3f} vs {r['mean_b']:>7.3f}  "
            f"{r['diff']:+.3f} [{r['ci'][0]:+.3f},{r['ci'][1]:+.3f}] "
            f"fold[{r['fold_ci'][0]:+.3f},{r['fold_ci'][1]:+.3f}] "
            f"med {r['median']:+.3f} {r['W']}W/{r['L']}L f{r['folds_same_sign']}/{r['n_folds']}")


def main():
    d = json.load(open(os.path.join(ROOT, "s19", "results", "qb_main.json")))
    pdbs = sorted([k for k in d if not k.startswith("_")])
    tg = {t["pdb"]: t for t in I.targets()}
    assert len(pdbs) == 126, len(pdbs)
    folds = np.array([int(tg[p]["fold"]) for p in pdbs])
    arms = list(d[pdbs[0]]["arms"].keys())
    for p in pdbs:
        assert set(d[p]["arms"]) == set(arms), p

    def G(a, k):
        return np.array([d[p]["arms"][a][k] for p in pdbs], float)

    rng = np.random.default_rng(20200907)
    rep = {"n": len(pdbs), "arms": arms, "mde": MDE,
           "source": "s19/results/qb_main.json", "written": d.get("_written")}

    idres = np.array([max(abs(d[p]["arms"][a]["identity_resid"]) for a in arms) for p in pdbs])
    rep["md_identity_max_resid"] = float(idres.max())
    print(f"(M,D) identity  readout^2 = M^2 - D^2 : max |resid| = {idres.max():.3e} A^2  "
          f"({'PASS' if idres.max() < 1e-9 else 'FAIL'})\n")

    print(f"{'arm':<14}{'realised':>10}{'sel_best':>10}{'gen_best':>10}{'gen@500':>9}"
          f"{'M':>8}{'D':>8}{'cov<2.5':>9}{'evals':>8}")
    tab = {}
    for a in arms:
        row = {k: float(G(a, k).mean()) for k in
               ("realised_ORACLE", "sel_best_ORACLE", "gen_best_ORACLE",
                "gen_best_sub500_ORACLE", "M", "D", "cov_2.5_ORACLE", "evals",
                "avg_ORACLE", "obj_best")}
        tab[a] = row
        print(f"{a:<14}{row['realised_ORACLE']:>10.3f}{row['sel_best_ORACLE']:>10.3f}"
              f"{row['gen_best_ORACLE']:>10.3f}{row['gen_best_sub500_ORACLE']:>9.3f}"
              f"{row['M']:>8.3f}{row['D']:>8.3f}{row['cov_2.5_ORACLE']:>9.1f}{row['evals']:>8.0f}")
    rep["table"] = tab

    bq = min(Q, key=lambda a: tab[a]["realised_ORACLE"])
    bc = min(CLASSICAL, key=lambda a: tab[a]["realised_ORACLE"])
    print(f"\n=== P1 (PRIMARY, realised) === best quantum={bq}  best classical={bc}\n")
    p1 = [cmp_arm(f"P1 {bq} - {bc}", G(bq, "realised_ORACLE"), G(bc, "realised_ORACLE"), folds, rng),
          cmp_arm(f"{bq} - q_untrained", G(bq, "realised_ORACLE"),
                  G("q_untrained", "realised_ORACLE"), folds, rng)]
    for z in ZEROINFO + ["pool500"]:
        p1.append(cmp_arm(f"{bq} - {z}", G(bq, "realised_ORACLE"),
                          G(z, "realised_ORACLE"), folds, rng))
    for c in ("c_prod0.25", "c_chain0.25", "c_metroL", "c_anneal", "c_lbfgs"):
        p1.append(cmp_arm(f"{bq} - {c}", G(bq, "realised_ORACLE"),
                          G(c, "realised_ORACLE"), folds, rng))
    for r in p1:
        print(line(r))
    rep["P1"] = p1
    rep["best_quantum"] = bq
    rep["best_classical"] = bc

    print("\n=== P3 (generation ceiling) ===\n")
    bqg = min(Q, key=lambda a: tab[a]["gen_best_ORACLE"])
    bcg = min(CLASSICAL, key=lambda a: tab[a]["gen_best_ORACLE"])
    p3 = [cmp_arm(f"P3 full 8192: {bqg} - {bcg}", G(bqg, "gen_best_ORACLE"),
                  G(bcg, "gen_best_ORACLE"), folds, rng),
          cmp_arm(f"P3 @500: {bqg} - {bcg}", G(bqg, "gen_best_sub500_ORACLE"),
                  G(bcg, "gen_best_sub500_ORACLE"), folds, rng),
          cmp_arm(f"{bqg} - q_untrained (gen)", G(bqg, "gen_best_ORACLE"),
                  G("q_untrained", "gen_best_ORACLE"), folds, rng),
          cmp_arm(f"{bqg}@500 - pool500 (matched count)", G(bqg, "gen_best_sub500_ORACLE"),
                  G("pool500", "gen_best_ORACLE"), folds, rng),
          cmp_arm(f"{bqg} - c_marg (gen)", G(bqg, "gen_best_ORACLE"),
                  G("c_marg", "gen_best_ORACLE"), folds, rng),
          cmp_arm(f"{bqg} - c_helix (gen)", G(bqg, "gen_best_ORACLE"),
                  G("c_helix", "gen_best_ORACLE"), folds, rng)]
    for r in p3:
        print(line(r))
    rep["P3"] = p3

    def eps(px, py, others):
        return float(min(max(mc - px, py - dc) for mc, dc in others))

    print("\n=== P2 (eps-dominance, (M,D) plane) ===\n")
    Mc = {a: G(a, "M") for a in arms}
    Dc = {a: G(a, "D") for a in arms}
    cls = CLASSICAL + ZEROINFO + ["pool500"]
    p2 = {}
    for a in Q + ["q_untrained"]:
        e = np.array([eps(Mc[a][t], Dc[a][t], [(Mc[c][t], Dc[c][t]) for c in cls])
                      for t in range(len(pdbs))])
        m, lo, hi = boot(e, rng)
        p2[a] = {"mean": m, "ci": [lo, hi], "frac_pos": float((e > 0).mean())}
    for c in cls:
        sib = [x for x in cls if x != c]
        e = np.array([eps(Mc[c][t], Dc[c][t], [(Mc[s][t], Dc[s][t]) for s in sib])
                      for t in range(len(pdbs))])
        m, lo, hi = boot(e, rng)
        p2[c] = {"mean": m, "ci": [lo, hi], "frac_pos": float((e > 0).mean()), "loo": True}
    for k, v in sorted(p2.items(), key=lambda kv: -kv[1]["mean"]):
        print(f"  {k:<16}{'(LOO)' if v.get('loo') else '     '} eps = {v['mean']:+.4f} "
              f"[{v['ci'][0]:+.4f},{v['ci'][1]:+.4f}]  frac>0 {v['frac_pos']:.2f}")
    rep["P2"] = p2

    json.dump(rep, open(os.path.join(OUT, "qb_close.json"), "w"), indent=1)
    open(os.path.join(OUT, "COMPLETE"), "w").write(f"n={len(pdbs)}\n")
    print(f"\nwritten {OUT}")


if __name__ == "__main__":
    main()
