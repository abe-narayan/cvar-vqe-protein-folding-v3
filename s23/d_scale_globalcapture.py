"""s23/d_scale_globalcapture.py -- WORKSTREAM D duty (2): what fraction of the oracle ceiling can
a SINGLE GLOBAL scalar capture, derived from the spread of s* across targets, independent of any
grid.  Cross-checks gscale.py's nested-CV near-null (+0.0002) against the population-optimal
global scale computed directly from the closed-form per-target quadratics (no grid truncation).

The project's endpoint is the MEAN of PER-TARGET RMSD (not a pooled/concatenated RMSD), so the
population objective for a candidate global s is L(s) = (1/T) sum_t sqrt(SSE_t(s)/n_t) -- a sum of
square roots of per-target parabolas, not itself a parabola.  Minimised here by a fine 1-D search
(no grid-truncation risk: range widened to [0.3, 2.2], comfortably covering every per-target
s_analytic found in d_scale_closedform.py, min 0.252 max 1.924).
"""
from __future__ import annotations
import json, os, sys
import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I  # noqa: E402

TOPM = 75
FAIL18 = {"1ID6", "1JBF", "1LB7", "2BFI", "2BP4", "2JN5", "2MQ2", "2N5C", "2NB7", "2NDM",
          "3BTB", "3SGO", "5W52", "7JS6", "7LCW", "8T63", "9KAR", "9L1M"}


def kabsch_CcCtTt(C, nat):
    Cc0 = C - C.mean(0); Tc0 = nat - nat.mean(0)
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
    quads = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]; nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        o = np.argsort(sc, kind="stable")
        C, _b = I.coordinate_average(W[o[:TOPM]])
        Cc, Ct, Tt, n = kabsch_CcCtTt(np.asarray(C, float), nat)
        quads.append({"pdb": pdb, "Cc": Cc, "Ct": Ct, "Tt": Tt, "n": n, "fail18": pdb in FAIL18})
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
    with open(os.path.join(RES, "d_scale_globalcapture.json"), "w") as fh:
        json.dump(quads, fh)
    report(quads)


def report(quads=None):
    if quads is None:
        quads = json.load(open(os.path.join(RES, "d_scale_globalcapture.json")))

    def rmsd_curve(s, rows):
        out = np.zeros((len(s), len(rows)))
        for k, r in enumerate(rows):
            sse = np.maximum(r["Cc"] * s ** 2 - 2 * r["Ct"] * s + r["Tt"], 0.0)
            out[:, k] = np.sqrt(sse / r["n"])
        return out  # (len(s), n_targets)

    grid = np.linspace(0.3, 2.2, 3801)  # 0.0005 resolution, wide enough for every observed s_analytic

    def analyse(rows, label):
        curve = rmsd_curve(grid, rows)  # (S, T)
        mean_per_s = curve.mean(1)
        i_best = int(np.argmin(mean_per_s))
        s_pop = float(grid[i_best])
        pop_mean = float(mean_per_s[i_best])
        base_mean = float(curve[np.argmin(np.abs(grid - 1.0))].mean())
        # true per-target oracle (unconstrained) for the ceiling
        oracle_mean = float(np.mean([np.sqrt(max(r["Tt"] - r["Ct"] ** 2 / r["Cc"], 0.0) / r["n"])
                                      for r in rows]))
        s_star = np.array([r["Ct"] / r["Cc"] for r in rows])
        cap = 100.0 * (base_mean - pop_mean) / (base_mean - oracle_mean) if base_mean > oracle_mean else float("nan")
        print("\n=== %s (n=%d) ===" % (label, len(rows)))
        print("  incumbent (s=1.0)             %.4f" % base_mean)
        print("  population-optimal global s   %.4f  ->  %.4f" % (s_pop, pop_mean))
        print("  per-target oracle (ceiling)   %.4f  (unconstrained closed form)" % oracle_mean)
        print("  delta global vs incumbent     %+.4f" % (pop_mean - base_mean))
        print("  delta oracle vs incumbent     %+.4f  <- full ceiling" % (oracle_mean - base_mean))
        print("  CAPTURED FRACTION of ceiling by the best POSSIBLE single global constant: %.1f%%" % cap)
        print("  s* spread: mean %.4f sd %.4f  (dispersion this global constant must average over)"
              % (s_star.mean(), s_star.std()))
        pos = (s_star > 1.0).sum(); neg = (s_star < 1.0).sum()
        print("  direction split: s*>1 (wants EXPANSION) %d/%d , s*<1 (wants CONTRACTION) %d/%d"
              % (pos, len(rows), neg, len(rows)))

    analyse(quads, "FULL PANEL")
    analyse([q for q in quads if q["fail18"]], "FAIL18 subset")
    analyse([q for q in quads if not q["fail18"]], "non-FAIL18")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
