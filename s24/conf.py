"""s24/conf.py -- IF THE OUTPUT INHERITS THE PRIOR'S ERROR, WEIGHT THE PRIOR BY ITS OWN CONFIDENCE.
ONE GLOBAL EXPONENT, NESTED CV.

WHY, AND WHY THIS IS NOT ANOTHER WAY OF RE-CONSUMING THE PRIOR.  L5 measured three things at n=126:
the emitted cloud is CLOSER to the distogram's prediction than the native is (RMS 2.2358 vs 3.2973,
111W/15L, 3.0x MDE); the share of the prior's own error reappearing in the output is beta = 0.520
against a mismatched-distogram placebo floor of 0.352; and Workstream C found that applying the
selector raises alignment with the prior's error by ~0.20 in every arm, including a zero-information
one.  **Selection installs the prior's error into the output.**

The shipped score is `risk[pair, g(D_pair)].mean(1)` -- a UNIFORM mean over pairs.  Every pair
contributes equally regardless of how well the distogram actually knows that pair.  The distogram
publishes a per-pair `sd`.  So the single most direct consequence of L5 is: **weight each pair by the
prior's own confidence in it, and the transferred error should fall.**

    w_p  proportional to  sd_p^(-k),  normalised to sum to 1
    score_k(candidate) = sum_p w_p * risk[p, g(D_p)]

**k = 0 recovers the shipped score EXACTLY**, which is what makes this a clean null rather than a
reparameterisation: the arm contains its own incumbent as an interior point of the grid, and the
identity is asserted numerically before anything else is reported.

Project memory records that eleven ways of re-consuming the prior failed to survive dev
(`score-axis-does-not-transfer`).  Every one of those re-consumed the prior's POINT ESTIMATE or its
MAE -- a second read of the same number.  This does not touch the point estimate at all; it changes
the RELATIVE WEIGHT OF PAIRS inside the risk functional, using a quantity (`sd`) the shipped score
currently ignores entirely.  It is a different object, and it is the one L5 names.  Stated here so
the distinction is on the record and can be attacked rather than assumed.

ONE GLOBAL SCALAR, BY CONSTRAINT.  The standing finite-sample bound from this project's own
sigma ~ 0.41 A gives even a one-global-threshold router a 0.39 A generalisation gap at n ~ 100 per
fold, and richer classes need 280-1045 targets.  A single scalar fitted on training folds is the only
adaptive object this instrument supports, and five per-target routers have already failed exactly as
the bound predicts.  k is that scalar.  Nothing here is fitted per target.

OPERATOR FORKS, per rule 0.  NOTE: enumerated by the coordinator, who has a stake; Lane E has been
asked to re-enumerate.  Recorded, not hidden.

    functional     DECLARED the shipped Bayes-risk table `risk`, unmodified, re-weighted across pairs
                   by the distogram's own published `sd`.  NOT TAKEN modifying the risk table itself,
                   NOT TAKEN the posterior mean `expected`, and NOT TAKEN a weight fitted freely per
                   pair (which the sample size forbids).
    basis          DECLARED point cloud, medoid frame, on both sides -- the arm and its own null are
                   the same operator with k=0 vs k fitted.  NOT TAKEN a built chain.
    readout        DECLARED the shipped uniform top-75 coordinate average, unchanged, so the ONLY
                   thing that varies is which 75 candidates are selected.  NOT TAKEN re-optimising m,
                   NOT TAKEN any weighting of members (s23 L5: uniform is optimal).
    normalisation  DECLARED weights normalised to sum to 1 so that k=0 is bit-identical to the
                   shipped score, and sd floored at a small epsilon so a zero-sd pair cannot take
                   infinite weight.  NOT TAKEN unnormalised weights (which would change the score's
                   scale and nothing else) and NOT TAKEN z-scoring the risk per pair.
    null           DECLARED k = 0, the shipped score, on identical data -- does confidence weighting
                   beat NOT weighting.  NOT TAKEN the per-target oracle k, which is the ceiling under
                   test and is reported separately and labelled ORACLE.
    THE LABEL      DECLARED continuous Ca-RMSD, and SEPARATELY beta, the prior-inheritance
                   coefficient from s24/referent.py.  Both are reported because they answer different
                   questions: RMSD says whether it helps, beta says whether it helped BY THE
                   MECHANISM L5 predicts.  An RMSD gain with beta unmoved would mean the mechanism is
                   wrong even if the arm wins.  NOT TAKEN RMSD alone, and NOT TAKEN any binarised
                   win rate.

  H   weighting pairs by the prior's own confidence reduces the share of the prior's error that
      selection transfers, and therefore lowers RMSD.
  Falsifier  the nested-CV arm failing to beat k=0 past its own MDE with a fold-clustered CI
      excluding zero.  Secondary falsifier: an RMSD gain unaccompanied by a fall in beta, which would
      refute the stated mechanism regardless of the endpoint.
  Nesting  k is fitted on the 4 training folds and applied to the held-out fold.  The in-sample fit is
      reported separately and labelled as the optimism, never as the result.

  Natives are read for EVALUATION ONLY.  k is never fitted on a held-out fold.
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

TOPM = 75
KGRID = np.array([-1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
EPS = 1e-3


def _save(o, name="conf.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _avg(members):
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    return I.superpose_batch(members, members[b]).mean(0)


def _score_k(dg, D, k):
    """The shipped Bayes-risk score with pairs re-weighted by the prior's own confidence.

    k = 0 reproduces `I.shipped_score` exactly (asserted in run()).
    """
    grid = np.asarray(dg["grid"], float); risk = np.asarray(dg["risk"], float)
    g = np.clip(((np.asarray(D, float) - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    r = risk[np.arange(risk.shape[0])[None, :], g]                # (b, npairs)
    if k == 0.0:
        return r.mean(1)
    w = np.maximum(np.asarray(dg["sd"], float), EPS) ** (-float(k))
    w = w / w.sum()
    return (r * w[None, :]).sum(1)


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, m=%d, k grid %s" % (len(tg), TOPM, KGRID.tolist()), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        Wall = np.asarray(u["W"], float); nat = np.asarray(u["nat_ca"], float)
        n_res = int(u["n"]); i, j = I.pair_index(n_res)
        dg = I.distogram(pdb, u["seq"], u["fold"])
        Wpool = Wall[I.pool_idx(u)]
        D = I.pair_dists(Wpool, i, j)

        #: the k=0 identity gate -- if this ever fails the whole arm is void
        s0 = _score_k(dg, D, 0.0); ship = np.asarray(I.shipped_score(dg, D), float)
        assert np.allclose(s0, ship, atol=1e-12), (pdb, float(np.abs(s0 - ship).max()))

        Dt = I.pair_dists(nat[None], i, j)[0]
        Dh = np.asarray(dg["grid"], float)[np.argmin(np.asarray(dg["risk"], float), axis=1)]
        eP = Dh - Dt; den = float((eP * eP).sum())

        rm, be = [], []
        for k in KGRID:
            sk = _score_k(dg, D, float(k))
            sel = Wpool[np.argsort(sk, kind="stable")[:TOPM]]
            C = _avg(sel)
            rm.append(float(I.ca_rmsd(C, nat)))
            Dc = I.pair_dists(C[None], i, j)[0]
            be.append(float(((Dc - Dt) * eP).sum() / den) if den > 0 else float("nan"))

        rows.append({"pdb": pdb, "n": n_res, "fold": int(u["fold"]),
                     "rmsd": rm, "beta": be})
        if (c_i + 1) % 20 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    ok = len(rows) == len(tg) and all(len(r["rmsd"]) == len(KGRID) and len(r["beta"]) == len(KGRID)
                                      for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "kgrid": KGRID.tolist(), "topm": TOPM, "eps": EPS})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "conf.json")))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    R = np.array([r["rmsd"] for r in rows], float)          # (T, K)
    B = np.array([r["beta"] for r in rows], float)
    k0 = int(np.argmin(np.abs(KGRID)))
    assert KGRID[k0] == 0.0
    base = R[:, k0]; base_b = B[:, k0]
    rng = SD.stable_rng("conf", "rep")

    def st(x):
        x = np.asarray(x, float); n = len(x)
        se = x.std(ddof=1) / np.sqrt(n)
        F = sorted(set(fold))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se, float(np.percentile(fs, 2.5)),
                float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    print("\nn = %d.  Point cloud, shipped top-75 uniform average, only the SELECTION varies.\n"
          % len(rows))
    print("  k = 0 reproduces the shipped score exactly.  Incumbent %.4f, beta %.4f\n"
          % (base.mean(), base_b.mean()))
    print("  %8s%10s%10s%12s" % ("k", "RMSD", "beta", "vs k=0"))
    for a, k in enumerate(KGRID):
        print("  %8.2f%10.4f%10.4f%+12.4f" % (k, R[:, a].mean(), B[:, a].mean(),
                                              (R[:, a] - base).mean()))

    #: --- NESTED CV: fit k on the 4 training folds, apply to the held-out fold
    held = np.empty(len(rows)); held_b = np.empty(len(rows)); kused = np.empty(len(rows))
    for f in sorted(set(fold)):
        tr = fold != f; te = fold == f
        a = int(np.argmin(R[tr].mean(0)))
        held[te] = R[te][:, a]; held_b[te] = B[te][:, a]; kused[te] = KGRID[a]
    a_all = int(np.argmin(R.mean(0)))

    print("\n  PRIMARY (paired vs k = 0, the shipped score):")
    for lab, x in (("CONFIDENCE WEIGHT, NESTED CV  <-- the deployable arm", held - base),
                   ("in-sample best k (leaky, the optimism)", R[:, a_all] - base),
                   ("per-target ORACLE k (CEILING)", R.min(1) - base)):
        m, se, mde, flo, fhi, w = st(x)
        v = ("BEATS" if (fhi < 0 and abs(m) > mde) else
             "worse" if (flo > 0 and abs(m) > mde) else "not measured")
        tm = " [TYPE-M ZONE]" if mde > 0 and 0.7 <= abs(m) / mde <= 1.3 else ""
        print("    %-52s%+.4f SE %.4f MDE %.4f fold[%+.4f,%+.4f] %3dW/%3dL %s%s"
              % (lab, m, se, mde, flo, fhi, w, len(rows) - w, v, tm))
    print("    k selected by nested CV: mean %.3f, values %s"
          % (kused.mean(), sorted(set(kused.tolist()))))
    print("    in-sample best k = %.2f" % KGRID[a_all])

    print("\n  MECHANISM CHECK (the secondary falsifier):")
    m, se, mde, flo, fhi, w = st(held_b - base_b)
    print("    beta, nested arm minus k=0:  %+.4f SE %.4f MDE %.4f fold[%+.4f,%+.4f]  %3dW/%3dL"
          % (m, se, mde, flo, fhi, w, len(rows) - w))
    print("    An RMSD gain WITHOUT a fall in beta refutes the stated mechanism even if the arm wins.")
    print("\n  Falsifier: the nested arm failing to beat k=0 past its own MDE with a fold CI"
          " excluding zero.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
