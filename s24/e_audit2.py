"""s24/e_audit2.py -- LANE E, ROUND 2.  (C3) THE UNIFORM-SUBSET LADDER FOR qmatch/L3, AND
(B7) THE ARM-DISJOINT VERSION OF THE 0.693 RATIO.

(C3)  The coordinator's L3 argues that in `qmatch.py` the undeclared score-prefix fork runs TOWARD
      a bow, none appeared, and therefore "A and B' are parallel" is STRENGTHENED.  That is the
      coordinator reasoning about the coordinator's own result, so Lane E measures it.
      `qmatch.py:137` builds `concatenate([A[:a], B[:b]])` where BOTH `A` and `B'` are held in
      shipped-score order, so every interior arm keeps the best-scoring members of BOTH parents
      while the endpoints keep all 75.  In `biasalign` only one half was affected.
      NULL: the identical ladder with a uniform a-subset of A and a uniform b-subset of B'.
      The coordinator's claim is FALSIFIED if removing the fork moves the ladder DOWN -- i.e. if
      the fork was suppressing a bow rather than creating one.

(B7)  The BRIEF quotes `cos(A,C)/cos(C1,C2) = 0.693`.  Numerator and denominator SHARE the arm
      `e_C1` as well as the native referent, so it is not a ratio of two independent estimates.
      The arm-disjoint version is `cos(A,C2)/cos(C1,C2)`, whose numerator and denominator share
      only `e_C2`, plus the fully-disjoint-numerator average of the two.  Reported beside the
      quoted figure so the BRIEF can be corrected or confirmed.

OPERATOR FORKS (Lane E, no stake in the outcome).
    functional     DECLARED the shipped Bayes-risk score, unchanged, for A and B', reproducing
                   qmatch's own RNG stream bit-for-bit.  NOT TAKEN any re-scoring.
    basis          DECLARED point cloud, bias vectors in the native frame, RMSD in each set's own
                   medoid frame -- identical to the files under audit.  NOT TAKEN a built chain.
    readout        DECLARED uniform mean of exactly 75 members at every ladder point.
                   NOT TAKEN unmatched sizes.
    normalisation  DECLARED the bow in absolute A AND as a fraction of that ladder's S(0).
                   NOT TAKEN absolute A alone.
    null           DECLARED the uniform-subset ladder (both halves) as the fork's null, and the
                   analytic no-operator ladder as the shared-referent floor for c_eff.
                   NOT TAKEN "no bow", which assumes the conclusion.
    THE LABEL      DECLARED continuous Ca-RMSD and the raw cosines.  NOT TAKEN a binarised call.

  Natives are read for EVALUATION ONLY, after every candidate set is fixed.
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
from s24 import stats_lib as ST           # noqa: E402
from s24.e_bowaudit import _avg, _bias, _cos, _ladder, fit_c, MIX, TOPM   # noqa: E402

DRAW = 2000            # qmatch's B' candidate-set size; read from qmatch.py, not re-chosen
NRND = 4
OUT = os.path.join(RES, "e_audit2.json")
NEED = ("Lq_obs", "Lq_rnd", "Lq_ana", "cos_A_B", "cos_A_C1", "cos_A_C2", "cos_C1_C2")


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, m=%d, DRAW=%d, %d subset perms" % (len(tg), TOPM, DRAW, NRND), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        Wall = np.asarray(u["W"], float); nat = np.asarray(u["nat_ca"], float)
        N = len(Wall); n_res = int(u["n"])
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(n_res)
        idx = I.pool_idx(u); Wpool = Wall[idx]
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(Wpool, i, j)), float)
        A = Wpool[np.argsort(sc, kind="stable")[:TOPM]]

        #: --- qmatch's B', reproducing its RNG stream exactly
        rq = SD.stable_rng("qmatch", pdb)
        mask = np.ones(N, bool); mask[idx] = False
        pool_free = np.flatnonzero(mask)
        pick = rq.choice(pool_free, min(DRAW, len(pool_free)), replace=False)
        Wf = Wall[pick]
        scf = np.asarray(I.shipped_score(dg, I.pair_dists(Wf, i, j)), float)
        B = Wf[np.argsort(scf, kind="stable")[:TOPM]]

        #: --- biasalign's C1/C2, reproducing ITS stream exactly (for B7)
        rb = SD.stable_rng("biasalign", pdb)
        C1 = Wall[rb.choice(N, TOPM, replace=False)]
        C2 = Wall[rb.choice(N, TOPM, replace=False)]

        eA = _bias(_avg(A), nat); eB = _bias(_avg(B), nat)
        eC1 = _bias(_avg(C1), nat); eC2 = _bias(_avg(C2), nat)

        row = {"pdb": pdb, "n": n_res, "fold": int(u["fold"]),
               "cos_A_B": _cos(eA, eB), "cos_A_C1": _cos(eA, eC1), "cos_A_C2": _cos(eA, eC2),
               "cos_C1_C2": _cos(eC1, eC2),
               "norm_A": float(np.linalg.norm(eA)), "norm_B": float(np.linalg.norm(eB))}

        #: --- (C3) the qmatch ladder, as shipped and with the fork removed on BOTH halves
        row["Lq_obs"] = _ladder(A, B, nat)
        ra = SD.stable_rng("e_audit2", pdb)
        rl = [_ladder(A[ra.permutation(TOPM)], B[ra.permutation(TOPM)], nat) for _ in range(NRND)]
        row["Lq_rnd"] = {k: float(np.mean([d[k] for d in rl])) for k in rl[0]}
        row["Lq_rnd_sd"] = {k: float(np.std([d[k] for d in rl], ddof=1)) for k in rl[0]}
        row["Lq_ana"] = {"m%d_%d" % (a, b):
                         float(np.linalg.norm((1 - b / float(TOPM)) * eA + (b / float(TOPM)) * eB)
                               / np.sqrt(n_res)) for a, b in MIX}
        rows.append(row)
        if (c_i + 1) % 10 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            ST.save_atomic(OUT, {"rows": rows}, complete_keys=NEED, rows=rows,
                           n_expected=len(tg), module_file=__file__)

    ST.save_atomic(OUT, {"rows": rows, "topm": TOPM, "draw": DRAW, "mix": MIX, "nrnd": NRND},
                   complete_keys=NEED, rows=rows, n_expected=len(tg), module_file=__file__)
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    keys = ["m%d_%d" % (a, b) for a, b in MIX]
    lams = np.array([b / float(TOPM) for _, b in MIX])
    fold = np.array([r["fold"] for r in rows], int)
    pdbs = [r["pdb"] for r in rows]
    g = lambda k: np.array([r[k] for r in rows], float)     # noqa: E731

    def bows(tag):
        M = np.array([[r[tag][k] for k in keys] for r in rows], float)
        line = (1 - lams)[None, :] * M[:, :1] + lams[None, :] * M[:, -1:]
        return M, M - line

    print("\n" + "=" * 96)
    print("LANE E ROUND 2.  n = %d." % len(rows))
    print("=" * 96)

    print("\n  (C3) THE qmatch LADDER, WITH AND WITHOUT THE SCORE-PREFIX FORK ON BOTH HALVES")
    for tag, lab in (("Lq_obs", "A[:a] + B'[:b]              AS SHIPPED (score prefix, BOTH halves)"),
                     ("Lq_rnd", "A[perm][:a] + B'[perm][:b]  FORK REMOVED (uniform subsets)"),
                     ("Lq_ana", "analytic |(1-l)eA + l eB'|  NO OPERATOR")):
        M, B = bows(tag)
        print("\n    %s" % lab)
        print("      %-10s%9s%12s%12s%10s" % ("mixture", "RMSD", "vs 75/0", "vs line", "bow/S(0)"))
        for c, k in enumerate(keys):
            print("      %-10s%9.4f%+12.4f%+12.4f%+10.4f"
                  % (k, M[:, c].mean(), (M[:, c] - M[:, 0]).mean(), B[:, c].mean(),
                     (B[:, c] / M[:, 0]).mean()))

    print("\n    PAIRED, per target: SHIPPED ladder MINUS FORK-REMOVED ladder, at each point.")
    Mo, Bo = bows("Lq_obs"); Mr, Br = bows("Lq_rnd")
    for c, k in enumerate(keys):
        if c in (0, len(keys) - 1):
            continue
        r = ST.compare(Mo[:, c], Mr[:, c], fold, names=pdbs, label="%s  shipped vs fork-removed" % k)
        print("      %-10s %+.4f  med %+.4f  SE %.4f  MDE %.4f  e/MDE %+.2f  fold[%+.4f,%+.4f] "
              "%3dW/%3dL  %s" % (k, r["effect"], r["median_effect"], r["se"], r["mde"],
                                 r["effect_over_mde"], r["ci95_fold"][0], r["ci95_fold"][1],
                                 r["n_better"], r["n_worse"], r["verdict"]))
    print("      (NEGATIVE = the shipped ladder sits LOWER, i.e. the prefix CREATES apparent bow.")
    print("       POSITIVE = the prefix was SUPPRESSING a bow, and the coordinator's L3 reading")
    print("       -- 'the fork runs toward a bow and none appeared, so parallel is strengthened' --")
    print("       is WRONG in direction.)")

    print("\n    FITTED c_eff PER LADDER vs THE DIRECT COSINE cos(A,B')")
    d = g("cos_A_B")
    for tag, lab in (("Lq_obs", "qmatch as shipped"), ("Lq_rnd", "qmatch fork removed"),
                     ("Lq_ana", "analytic, no operator  <- SHARED-REFERENT FLOOR")):
        ce = np.array([fit_c([r[tag][k] for k in keys], lams)[0] for r in rows])
        f = np.isfinite(ce) & np.isfinite(d)
        print("      %-42s c_eff mn %+.4f  md %+.4f   direct mn %+.4f md %+.4f   rho %.3f"
              % (lab, ce[f].mean(), np.median(ce[f]), d[f].mean(), np.median(d[f]),
                 np.corrcoef(ce[f], d[f])[0, 1]))

    print("\n    THE ENDPOINT CONTRAST L3 ACTUALLY RESTS ON:")
    for tag, lab in (("Lq_obs", "shipped"), ("Lq_rnd", "fork removed")):
        M, _ = bows(tag)
        best = int(np.argmin(M.mean(0)))
        r = ST.compare(M[:, best], M[:, 0], fold, names=pdbs,
                       label="%s: best interior %s vs m75_0" % (lab, keys[best]))
        print("      %-14s best point %-8s %s" % (lab, keys[best], r["verdict"]))
        print(ST.fmt(r))

    print("\n  (B7) THE RATIO, WITH AND WITHOUT A SHARED ARM")
    ctl = g("cos_C1_C2")
    for k, lab in (("cos_A_C1", "cos(A,C1)   SHARES e_C1 with the control"),
                   ("cos_A_C2", "cos(A,C2)   SHARES e_C2 with the control"),
                   ("cos_A_B", "cos(A,B')   qmatch's quality-matched source")):
        v = g(k)
        print("    %-44s mean %+.4f  median %+.4f  sd %.4f" % (lab, v.mean(), np.median(v), v.std()))
    print("    %-44s mean %+.4f  median %+.4f  sd %.4f"
          % ("cos(C1,C2)  WITHIN-SOURCE CONTROL", ctl.mean(), np.median(ctl), ctl.std()))
    a1, a2 = g("cos_A_C1"), g("cos_A_C2")
    print("\n    ratio of means, as the BRIEF quotes it:")
    print("      cos(A,C1)/cos(C1,C2) = %.4f   <- the quoted 0.693, numerator SHARES e_C1"
          % (a1.mean() / ctl.mean()))
    print("      cos(A,C2)/cos(C1,C2) = %.4f   <- ARM-DISJOINT numerator" % (a2.mean() / ctl.mean()))
    print("      mean of the two      = %.4f   <- symmetric, uses both draws once"
          % (((a1 + a2) / 2).mean() / ctl.mean()))
    r = ST.compare(a2, a1, fold, names=pdbs, label="cos(A,C2) - cos(A,C1): is the shared arm biasing?")
    print()
    print(ST.fmt(r))
    print("\n    A ratio built from a shared arm is biased only if that arm's realisation is")
    print("    correlated with both quantities.  The contrast above measures exactly that, and")
    print("    C1 and C2 are exchangeable by construction, so this is the whole test.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
