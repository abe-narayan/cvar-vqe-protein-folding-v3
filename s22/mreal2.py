"""s22/mreal.py -- IS THE PER-TARGET m HEADROOM REAL, OR IS IT NOISE-FITTING?

The routing ceiling (L1) is 0.482 A over 13 arms, and **73% of it (0.350 A) sits in the m-ladder
alone** -- six arms of ONE operator on ONE source.  That looks like the most accessible routing
problem in the project.

    BUT SPRINT 21 L23 ESTABLISHED THE LADDER IS FLAT: m=150 is +0.023 and m=20 is +0.042
    against m=75, both CIs spanning zero, median tied-set size 2.0 of 6 rungs.

**A per-target argmin over six largely-flat, strongly-correlated arms captures per-target NOISE.**
Its ORACLE value is then guaranteed to look like headroom while containing nothing transferable --
the same min-of-K structure that made D9's top-M ladder look like ordering skill when it was
min-of-M.

THE DECISIVE TEST IS A SPLIT-HALF, AND IT NEEDS NO NEW MACHINERY.  Each target's K=500 pool splits
into two disjoint halves.  Build the whole m-ladder on each half independently.  Then:

    SELECT m on half A            ->   EVALUATE that m on half B
    versus the FIXED incumbent m  ->   EVALUATED on the same half B

If per-target m is a real property of the target, the selection transfers and beats fixed m.
If it is noise, it does not transfer, and the 0.350 A is an artefact of selecting on the outcome.

OPERATOR FORKS, per BRIEF SS4 rule 0.  I hold a DIRECTIONAL hypothesis (I expect no transfer), so
every fork is named together with the direction it would push.

    functional     DECLARED the shipped Bayes-risk score for the filter, as deployed.
                   NOT TAKEN the squared functional (S21 C2: interchangeable at the argmin here).
    basis          DECLARED the pool's own window coordinates, matched on BOTH sides of every
                   comparison.  NOT TAKEN the rebuild.
    readout        DECLARED top-m coordinate average -- the shipped operator, and the object whose
                   m is in question.  NOT TAKEN argmin or medoid; those are different arms, priced
                   separately in L1.
    normalisation  none; RMSD in Angstroms.  NOT TAKEN per-target z-scoring, which would hide that
                   the mean is set by the targets an arm FAILS on.
    null           FIXED m=75 evaluated on the SAME half.  This asks whether SELECTING m beats NOT
                   selecting it, on matched data.  NOT TAKEN the half-A oracle, which is the
                   quantity under suspicion and would beg the question.
    THE LABEL      DECLARED the continuous half-B RMSD.  NOT TAKEN "did the selected m win", a
                   binarised label whose threshold would depend on the ladder's flatness -- exactly
                   the covariate-dependent-threshold trap of S21 L25.

  Hypothesis   per-target m does NOT transfer; the m-ladder headroom is largely noise-fitting.
  Falsifier    if selected-m beats fixed-m on held-out halves past that comparison's own MDE with
               a CI excluding zero, per-target m IS real and routing on it is live.
  Null         fixed m=75 on the same half.
  Budget       126 targets x 8 repeats x 2 halves x 6 rungs.
  Promotion    none.  This prices a ceiling; it does not ship an arm.

TIES.  An argmin on a flat ladder reads the index order -- the trap that once manufactured a
1.386 A winner in this project.  The selected set is the FULL TIED SET at the half's minimum and the
outcome is averaged over it, which is the expectation under random tie-breaking.
"""
from __future__ import annotations

import json
import os
import sys

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

MS = (500, 150, 75, 20, 5, 1)
REPS = 8
FIXED = 75


def _save(o, name="mreal.json"):
    p = os.path.join(RES, name)
    t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, reps: %d" % (len(tg), REPS), flush=True)
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]
        nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"])
        i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        rng = SD.stable_rng(pdb, "s22mreal")
        rep = []
        for _r in range(REPS):
            perm = rng.permutation(len(W))
            halves = (perm[:len(W) // 2], perm[len(W) // 2:])
            lad = []
            for h in halves:
                o = h[np.argsort(sc[h], kind="stable")]
                row = []
                for m in MS:
                    sel = o[:min(m, len(o))]
                    C, _b = I.coordinate_average(W[sel])
                    row.append(float(I.ca_rmsd(np.asarray(C, float), nat)))
                lad.append(row)
            A, B = np.array(lad[0]), np.array(lad[1])
            tiedA = np.flatnonzero(A <= A.min() + 1e-12)
            tiedB = np.flatnonzero(B <= B.min() + 1e-12)
            rep.append({"selA_onB": float(B[tiedA].mean()), "selB_onA": float(A[tiedB].mean()),
                        "fixed_onB": float(B[MS.index(FIXED)]), "fixed_onA": float(A[MS.index(FIXED)]),
                        "oracleB": float(B.min()), "oracleA": float(A.min()),
                        "n_tiedA": int(len(tiedA)),
                        #: the two controls that separate PER-TARGET selection from a GLOBAL m fix
                        "ladA": A.tolist(), "ladB": B.tolist(),
                        "randm_onB": float(B[rng.integers(0, len(MS))])})
        agg = lambda k: float(np.mean([x[k] for x in rep]))            # noqa: E731
        rows.append({"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]),
                     "ladA_mean": np.mean([x["ladA"] for x in rep], 0).tolist(),
                     "ladB_mean": np.mean([x["ladB"] for x in rep], 0).tolist(),
                     "randm_heldout": agg("randm_onB"),
                     "sel_heldout": 0.5 * (agg("selA_onB") + agg("selB_onA")),
                     "fixed_heldout": 0.5 * (agg("fixed_onB") + agg("fixed_onA")),
                     "oracle_insample": 0.5 * (agg("oracleA") + agg("oracleB")),
                     "n_tied": agg("n_tiedA")})
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)}, "mreal2.json")

    need = ("sel_heldout", "fixed_heldout", "oracle_insample", "n_tied")
    ok = len(rows) == len(tg) and all(all(np.isfinite(r[k]) for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "reps": REPS, "ms": list(MS), "fixed": FIXED}, "mreal2.json")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "mreal2.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)                # noqa: E731
    rng = SD.stable_rng("mreal", "rep")

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        b = x[rng.integers(0, k, size=(4000, k))].mean(1)
        return (x.mean(), se, 2.8016 * se,
                float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)), int((x < 0).sum()))

    print("\nn = %d targets, %d split-half repeats, half-pool of 250.\n" % (len(rows), REPS))
    print("  %-36s%8s" % ("arm", "RMSD"))
    for k, lab in (("fixed_heldout", "FIXED m=75, held out"),
                   ("sel_heldout", "SELECTED m (half A) -> half B"),
                   ("oracle_insample", "IN-SAMPLE oracle m (the suspect)")):
        print("  %-36s%8.3f" % (lab, g(k).mean()))

    m, se, mde, lo, hi, w = st(g("sel_heldout") - g("fixed_heldout"))
    verdict = "TRANSFERS" if (hi < 0 and abs(m) > mde) else "DOES NOT TRANSFER"
    print("\n  PRIMARY  selected - fixed, HELD OUT : %+.3f SE %.3f MDE %.3f [%+.3f,%+.3f] %dW/%dL"
          % (m, se, mde, lo, hi, w, len(rows) - w))
    print("           -> %s" % verdict)
    m2, se2, mde2, lo2, hi2, _w2 = st(g("oracle_insample") - g("fixed_heldout"))
    print("  in-sample oracle - fixed            : %+.3f [%+.3f,%+.3f]   <- the APPARENT headroom"
          % (m2, lo2, hi2))
    if m2 < 0:
        frac = 100.0 * min(max(m / m2, 0.0), 1.0)
        print("\n  FRACTION OF THE APPARENT HEADROOM THAT ACTUALLY TRANSFERS: %.0f%%" % frac)
    print("  median tied rungs at the half-ladder minimum: %.2f of %d" % (np.median(g("n_tied")), len(MS)))
    print("\n  Falsifier was: selected beats fixed past its own MDE with a CI excluding zero.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
