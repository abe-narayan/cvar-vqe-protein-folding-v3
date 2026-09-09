"""s23/d_scale_placebo2.py -- WORKSTREAM D: CORRECTED PLACEBO, self-caught defect in placebo.py.

d_scale_placebo.py's within-target atom-permutation placebo is CONFOUNDED: permuting nat's row
order destroys coherent chain geometry entirely, producing a degenerate rotation-fit problem (mean
RMSD 7.5 A vs the real 3.05 A) that is not a fair-difficulty null.  Its "387% mechanical share"
number is an artefact of comparing improvement fractions across two problems of very different
difficulty and is WITHDRAWN, reported rather than silently dropped (project convention, c.f. S22
M7).

CORRECTED CONSTRUCTION.  Pair each target's OWN averaged pool structure C_t with the NATIVE
structure of a DIFFERENT, randomly chosen target of the SAME chain length n (a real, coherently
folded chain -- same difficulty class as the true pairing, just the wrong identity).  This is the
project's standard cross-target placebo pattern (c.f. S22 M1's within-length-decile shuffle).  20
random cross-target draws per target, averaged.
"""
from __future__ import annotations
import json, os, sys
from collections import defaultdict
import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402

TOPM = 75
NDRAW = 20


def kabsch_quad(C, T):
    Cc0 = C - C.mean(0); Tc0 = T - T.mean(0)
    H = Cc0.T @ Tc0
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.eye(3); D[2, 2] = d
    R = Vt.T @ D @ U.T
    Crot = Cc0 @ R.T
    n = Crot.shape[0]
    Cc = float((Crot ** 2).sum()); Ct = float((Crot * Tc0).sum()); Tt = float((Tc0 ** 2).sum())
    return Cc, Ct, Tt, n


def run():
    tg = I.targets()
    by_len = defaultdict(list)
    for t in tg:
        by_len[t["n"]].append(t["pdb"])

    nat_cache = {}
    C_cache = {}
    for t in tg:
        pdb = t["pdb"]
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]
        nat_cache[pdb] = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        o = np.argsort(sc, kind="stable")
        C, _b = I.coordinate_average(W[o[:TOPM]])
        C_cache[pdb] = np.asarray(C, float)

    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]; n = t["n"]
        C = C_cache[pdb]; nat = nat_cache[pdb]
        Cc, Ct, Tt, nn = kabsch_quad(C, nat)
        rmsd1 = float(np.sqrt(max(Cc - 2 * Ct + Tt, 0.0) / nn))
        rmsd_star_real = float(np.sqrt(max(Tt - Ct * Ct / Cc, 0.0) / nn))

        others = [p for p in by_len[n] if p != pdb]
        rng = SD.stable_rng(pdb, "s23dscaleplacebo2")
        if not others:
            rows.append({"pdb": pdb, "n": int(n), "fold": int(t["fold"]), "no_peer": True,
                         "real_rmsd1": rmsd1, "real_star": rmsd_star_real})
            continue
        draws = rng.choice(others, size=min(NDRAW, len(others) * 5), replace=True)
        p1, pstar = [], []
        for op in draws:
            Ccp, Ctp, Ttp, npp = kabsch_quad(C, nat_cache[op])
            p1.append(float(np.sqrt(max(Ccp - 2 * Ctp + Ttp, 0.0) / npp)))
            pstar.append(float(np.sqrt(max(Ttp - Ctp * Ctp / Ccp, 0.0) / npp)))
        rows.append({
            "pdb": pdb, "n": int(n), "fold": int(t["fold"]), "no_peer": False,
            "real_rmsd1": rmsd1, "real_star": rmsd_star_real, "real_delta": rmsd_star_real - rmsd1,
            "placebo_rmsd1_mean": float(np.mean(p1)), "placebo_star_mean": float(np.mean(pstar)),
            "placebo_delta_mean": float(np.mean(np.array(pstar) - np.array(p1))),
        })
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)

    with open(os.path.join(RES, "d_scale_placebo2.json"), "w") as fh:
        json.dump({"rows": rows, "ndraw": NDRAW}, fh)
    report(rows)


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "d_scale_placebo2.json")))["rows"]
    rows = [r for r in rows if not r.get("no_peer")]
    g = lambda k: np.array([r[k] for r in rows], float)  # noqa: E731
    real1 = g("real_rmsd1"); realstar = g("real_star")
    plac1 = g("placebo_rmsd1_mean"); placstar = g("placebo_star_mean")

    print("\nn = %d (targets with a same-length peer).\n" % len(rows))
    print("  REAL      : rmsd1 %.4f -> oracle %.4f   delta %+.4f  (frac %.1f%%)"
          % (real1.mean(), realstar.mean(), (realstar - real1).mean(),
             100 * (real1 - realstar).mean() / real1.mean()))
    print("  PLACEBO   : rmsd1 %.4f -> oracle %.4f   delta %+.4f  (frac %.1f%%)"
          % (plac1.mean(), placstar.mean(), (placstar - plac1).mean(),
             100 * (plac1 - placstar).mean() / plac1.mean()))
    frac_abs = 100.0 * (placstar - plac1).mean() / (realstar - real1).mean()
    print("\n  mechanical-artefact share of the REAL in-sample magnitude (absolute A): %.1f%%" % frac_abs)
    print("  mechanical-artefact share, RELATIVE/fractional basis                  : %.1f%%"
          % (100 * (100 * (plac1 - placstar).mean() / plac1.mean())
             / (100 * (real1 - realstar).mean() / real1.mean())))


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
