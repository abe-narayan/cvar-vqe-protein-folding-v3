"""s23/diverse.py -- THE ONE THING L5 AND L9 ACTUALLY PREDICT: DELIBERATELY DIVERSIFY THE
RETAINED SET.  A DIRECTIONAL PREDICTION WITH AN INTERNAL GRADIENT, AND A MATCHED NULL.

WHY THIS ARM EXISTS, AND WHY IT IS NOT ANOTHER CONSENSUS TWEAK.  Two independent results this
sprint say the same thing about the shipped operator:

  L5  removing 4 members at RANDOM from the top-75 costs nothing; removing the 4 most
      geometrically DEVIANT members costs +0.142 A, CI excluding zero.  The outliers carry error
      that CANCELS.
  L9  the exact decomposition mean_k|e_k|^2 = |ebar|^2 + mean_k|d_k|^2 puts 68% of the pool's
      squared error in the COMMON term and 32% in the idiosyncratic one.  Averaging removes the
      32% and cannot touch the 68%.

Every arm this project has run on the retained set TIGHTENS it -- medoid, clustering, sharp
weights, small m, amplitude weighting -- and every one of them failed.  The mechanism says why:
they attack the cancellation.  **The mechanism's own prediction is the opposite move, and it has
never been tested: choose the 75 members to be as MUTUALLY DIVERSE as possible, at matched
quality band and matched set size.**  If diversity is what makes averaging work, buying more of it
should pay, and buying less of it should cost.

MATCHED QUALITY IS THE WHOLE DESIGN.  Simply widening m mixes diversity with quality and has
already been measured (global best m = 75).  Here every arm draws exactly 75 members from the SAME
top-150 score band, so the arms differ in DIVERSITY and in nothing else.

  A0  incumbent          top-75 by the shipped score              (the deployment reference)
  A1  DIVERSE            greedy farthest-point 75 of the top-150  (the mechanism's prediction)
  A2  RANDOM             75 drawn uniformly from the top-150      (the MATCHED NULL for A1)
  A3  TIGHT              greedy closest-point 75 of the top-150   (reversed-direction control)

**The prediction is an ORDER, not a single contrast: A1 < A2 < A3.**  A directional gradient with
a reversed arm is much harder to fake than one winning comparison, and A3 is the adversarial
control the BRIEF requires: if A3 does not lose, the mechanism is wrong even if A1 happens to win.

OPERATOR FORKS, per BRIEF SS4 rule 0.  Each names the alternative not taken.
NOTE ON RULE 0: enumerated by the coordinator, who has a stake; no independent party enumerated
them.  Recorded, not hidden.

    functional     DECLARED the shipped Bayes-risk distogram score, as deployed, for the band.
                   NOT TAKEN the squared functional, or any re-scoring.
    basis          DECLARED point cloud, each arm in its OWN medoid frame, exactly as
                   `I.coordinate_average` does for the incumbent.  NOT TAKEN a rebuilt chain, and
                   NOT TAKEN forcing all arms into one shared frame (which would make the arms
                   differ by frame as well as by membership).
    readout        DECLARED the uniform mean of the 75 selected members -- uniform because L5
                   showed uniform is optimal.  NOT TAKEN any weighting.
    normalisation  DECLARED diversity as pairwise Ca-RMSD, the same metric the medoid already
                   minimises, so "diverse" means diverse in the metric the operator lives in.
                   NOT TAKEN a distance-matrix or torsion-space diversity.
    null           DECLARED A2, uniform random 75 from the same top-150 band, averaged over 8
                   seeded draws.  It matches band, size and readout and carries NO diversity
                   information -- the only difference from A1 is the selection rule.
                   NOT TAKEN A0 as the null: A0 is a different band and is reported separately as
                   the deployment reference.
    THE LABEL      DECLARED continuous Ca-RMSD.  NOT TAKEN any binarised win rate.

  Band K = 150 = 2 x the shipped m, fixed a priori, not tuned.  m = 75, the shipped width.
  H  diversity is the active ingredient in coordinate averaging, so A1 < A2 < A3 in RMSD.
  Falsifier  A1 failing to beat A2 past its own MDE with a fold-clustered CI excluding zero.
     Secondary falsifier: A3 failing to LOSE to A2 on the same standard, which would refute the
     mechanism regardless of what A1 does.
  Native structures are read for EVALUATION ONLY and never enter a selection rule.
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
BAND = 150
NDRAW = 8


def _save(o, name="diverse.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _greedy(P, m, far=True):
    """Greedy farthest-point (far=True) or closest-point (far=False) selection of m rows.

    Seeded with row 0 -- the BEST-SCORING member of the band -- so both directions start from the
    same, score-determined, native-free point and differ only in the direction of the greedy step.
    """
    k = len(P)
    chosen = [0]
    d = P[0].copy()
    d[0] = np.inf if not far else -np.inf
    for _ in range(m - 1):
        c = int(np.argmax(d)) if far else int(np.argmin(d))
        chosen.append(c)
        d = np.minimum(d, P[c])
        d[chosen] = -np.inf if far else np.inf
    return np.asarray(chosen, int)


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, band %d, m %d, %d random draws" % (len(tg), BAND, TOPM, NDRAW), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]; nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        o = np.argsort(sc, kind="stable")
        B = W[o[:BAND]]
        P = I.pairwise_rmsd(B)

        def arm(sel):
            C, _ = I.coordinate_average(B[sel])
            return float(I.ca_rmsd(C, nat)), float(P[np.ix_(sel, sel)].sum() / (len(sel) ** 2 - len(sel)))

        s_div = _greedy(P, TOPM, far=True)
        s_tgt = _greedy(P, TOPM, far=False)
        r0, sp0 = arm(np.arange(TOPM))                    # A0 incumbent, top-75 by score
        r1, sp1 = arm(s_div)                              # A1 diverse
        r3, sp3 = arm(s_tgt)                              # A3 tight
        rng = SD.stable_rng("diverse", pdb)
        rr, ss = [], []
        for _ in range(NDRAW):                            # A2 matched random null
            s = rng.choice(BAND, TOPM, replace=False)
            a, b = arm(s); rr.append(a); ss.append(b)

        rows.append({
            "pdb": pdb, "n": int(u["n"]), "fold": int(u["fold"]),
            "A0_incumbent": r0, "A1_diverse": r1, "A3_tight": r3,
            "A2_random": float(np.mean(rr)), "A2_sd": float(np.std(rr, ddof=1)),
            "A2_draws": [float(x) for x in rr],
            "spread_A0": sp0, "spread_A1": sp1, "spread_A2": float(np.mean(ss)), "spread_A3": sp3,
        })
        if (c_i + 1) % 20 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("A0_incumbent", "A1_diverse", "A2_random", "A3_tight", "spread_A1")
    ok = len(rows) == len(tg) and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "band": BAND, "topm": TOPM, "ndraw": NDRAW})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "diverse.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    fold = g("fold"); rng = SD.stable_rng("diverse", "rep")

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        F = sorted(set(fold.astype(int)))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se,
                float(np.percentile(fs, 2.5)), float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    print("\nn = %d.  All arms: %d members from the top-%d score band, uniform coordinate average,"
          " point cloud.\n" % (len(rows), TOPM, BAND))
    print("  %-34s%9s   %s" % ("arm", "RMSD", "mean pairwise spread of the retained set"))
    for lab, k, sk in (("A0 incumbent  top-75 by score", "A0_incumbent", "spread_A0"),
                       ("A1 DIVERSE    farthest-point", "A1_diverse", "spread_A1"),
                       ("A2 RANDOM     matched null", "A2_random", "spread_A2"),
                       ("A3 TIGHT      closest-point", "A3_tight", "spread_A3")):
        print("  %-34s%9.4f   %.3f A" % (lab, g(k).mean(), g(sk).mean()))

    print("\n  PRIMARY and controls (paired):")
    for lab, x in (("A1 DIVERSE  -  A2 RANDOM   <-- PRIMARY", g("A1_diverse") - g("A2_random")),
                   ("A3 TIGHT    -  A2 RANDOM   <-- reversed control",
                    g("A3_tight") - g("A2_random")),
                   ("A1 DIVERSE  -  A0 incumbent", g("A1_diverse") - g("A0_incumbent")),
                   ("A2 RANDOM   -  A0 incumbent", g("A2_random") - g("A0_incumbent")),
                   ("A3 TIGHT    -  A0 incumbent", g("A3_tight") - g("A0_incumbent"))):
        m, se, mde, flo, fhi, w = st(x)
        v = "BEATS" if (fhi < 0 and abs(m) > mde) else ("WORSE" if (flo > 0 and abs(m) > mde)
                                                        else "not measured")
        print("    %-42s%+.4f SE %.4f MDE %.3f fold[%+.4f,%+.4f] %3dW/%3dL  %s"
              % (lab, m, se, mde, flo, fhi, w, len(rows) - w, v))

    d = g("A2_sd").mean()
    print("\n  sd of A2 across its %d draws, per target: %.4f A -- the null's own draw noise" % (NDRAW, d))
    print("  Falsifier: A1 failing to beat A2 past its MDE with a fold CI excluding zero.")
    print("  Secondary falsifier: A3 failing to LOSE to A2 on the same standard.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
