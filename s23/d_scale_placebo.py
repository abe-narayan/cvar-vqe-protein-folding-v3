"""s23/d_scale_placebo.py -- WORKSTREAM D: HOW BIG IS THE MECHANICAL SELECTION-ARTEFACT ALONE?

d_scale_closedform.py proves SSE(1)-SSE(s*) = (Cc-Ct)^2/Cc >= 0 ALWAYS for the in-sample fit
(equality iff s*=1 exactly), so a near-100% win rate is GUARANTEED by construction and is not, by
itself, evidence of real per-target structure.  d_scale_transfer.py already answers the decisive
question (does it generalise -- yes, ~99%).  This script answers a narrower, complementary one:
how much of the in-sample MAGNITUDE (-0.34 A) would a pure mechanical artefact produce with ZERO
real atom correspondence, holding each target's own point-cloud GEOMETRY fixed?

METHOD.  For each target, randomly permute the row order of `nat` (breaks the true residue-to-
residue correspondence between the averaged pool structure C and the native structure T, while
leaving both point clouds' own shapes -- Cc, Tt, Rg, bond-length distribution -- untouched).
Recompute the FULL closed-form fit (rotation + scale) on the permuted pairing.  Average over 20
random permutations per target.  This is the project's own placebo convention (c.f. S22 M1's
within-length-decile shuffle).
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

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402

TOPM = 75
NPERM = 20


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
    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]; nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        o = np.argsort(sc, kind="stable")
        C, _b = I.coordinate_average(W[o[:TOPM]])
        C = np.asarray(C, float)

        Cc, Ct, Tt, n = kabsch_quad(C, nat)
        rmsd1 = float(np.sqrt(max(Cc - 2 * Ct + Tt, 0.0) / n))
        rmsd_star_real = float(np.sqrt(max(Tt - Ct * Ct / Cc, 0.0) / n))

        rng = SD.stable_rng(pdb, "s23dscaleplacebo")
        placebo_rmsd1, placebo_star = [], []
        for _p in range(NPERM):
            perm = rng.permutation(n)
            Ccp, Ctp, Ttp, np_ = kabsch_quad(C, nat[perm])
            placebo_rmsd1.append(float(np.sqrt(max(Ccp - 2 * Ctp + Ttp, 0.0) / np_)))
            placebo_star.append(float(np.sqrt(max(Ttp - Ctp * Ctp / Ccp, 0.0) / np_)))

        rows.append({
            "pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]),
            "real_rmsd1": rmsd1, "real_star": rmsd_star_real, "real_delta": rmsd_star_real - rmsd1,
            "placebo_rmsd1_mean": float(np.mean(placebo_rmsd1)),
            "placebo_star_mean": float(np.mean(placebo_star)),
            "placebo_delta_mean": float(np.mean(np.array(placebo_star) - np.array(placebo_rmsd1))),
            "placebo_win_rate": float(np.mean(np.array(placebo_star) <= np.array(placebo_rmsd1) + 1e-9)),
        })
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)

    with open(os.path.join(RES, "d_scale_placebo.json"), "w") as fh:
        json.dump({"rows": rows, "nperm": NPERM}, fh)
    report(rows)


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "d_scale_placebo.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)  # noqa: E731
    print("\nn = %d.  %d random atom-order permutations of `nat` per target (geometry held fixed).\n"
          % (len(rows), NPERM))
    print("  REAL in-sample delta (oracle - incumbent)     : %+.4f" % g("real_delta").mean())
    print("  PLACEBO in-sample delta (no true correspondence): %+.4f" % g("placebo_delta_mean").mean())
    print("  placebo win rate (placebo star <= placebo rmsd1): %.1f%%" % (100 * g("placebo_win_rate").mean()))
    frac = 100.0 * g("placebo_delta_mean").mean() / g("real_delta").mean()
    print("\n  MECHANICAL ARTEFACT SHARE OF THE IN-SAMPLE MAGNITUDE: %.1f%%" % frac)
    print("  (i.e. even with ZERO true atom correspondence, a 1-parameter in-sample scale fit")
    print("   recovers this fraction of the observed -0.34 A in-sample number by construction alone)")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
