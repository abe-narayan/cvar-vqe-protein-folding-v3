"""s23/errdecomp.py -- WHY THE PER-TARGET SCALE IS UNREACHABLE: A DECOMPOSITION, AND THE ONE
PARAMETER-FREE ESTIMATOR THE THEORY ACTUALLY IMPLIES.

WHY THIS EXISTS.  Three lanes converged on the same statement this sprint by three routes:
coordinate averaging is a VARIANCE-REDUCTION operator (L5: removing 4 geometric outliers costs
+0.142 A while removing 4 at random costs nothing), the per-target scale s* is real at 126W/0L
(L2, ceiling -0.3403 in closed form, L6a), and s* is a property of the (pool, reference) PAIR
rather than of the fold (L6d: a mismatched same-length native produces an equal or larger
apparent gain).  L6d is an EMPIRICAL placebo.  This file asks whether the same conclusion follows
from a DECOMPOSITION -- i.e. whether it is a theorem about what a pool can see, not a fact about
this panel.

THE ALGEBRA, STATED BEFORE THE MEASUREMENT.  Write the m retained members, in the frame the
average is actually taken in, as w_k = t + e_k with t the native.  The emitted point cloud is
c = t + ebar, ebar = mean_k e_k.  Kabsch fits rotation and translation only, so with c and t both
centred the optimal scalar is exactly

    s* = <c,t> / |c|^2                                        (D's closed form, no grid)
       = (|c|^2 - <c,ebar>) / |c|^2  =  1 - <c,ebar>/|c|^2

so s* is set by the COMMON component ebar and by nothing else.  Decompose each member's error
about the ensemble mean, d_k = w_k - c = e_k - ebar.  Then

    mean_k |e_k|^2  =  |ebar|^2  +  mean_k |d_k|^2
                       ^^^^^^^^     ^^^^^^^^^^^^^^
                       COMMON       IDIOSYNCRATIC
                       invisible    the ensemble's own dispersion, fully observable

**The pool can measure the second term exactly and the first term not at all.**  A bias shared by
every member is, by construction, invisible from inside the ensemble: it moves c and every w_k
together and leaves every within-pool statistic unchanged.  If -- and only if -- the errors were
i.i.d. across members would the invisible term be recoverable from the visible one, via
|ebar|^2 = mean_k|d_k|^2 / (m-1).

THAT GIVES A REAL ARM, NOT ONLY A DIAGNOSIS.  Under the i.i.d. model the estimator

    s_hat = 1 - mean_k|d_k|^2 / ((m-1) * |c|^2)

is native-free, parameter-free, fitted to nothing, and derived rather than searched.  It is the
ONLY scale estimator this project has proposed that comes out of a model instead of a grid.  It
is therefore also a clean test OF the model: if the i.i.d. assumption held it would capture the
ceiling; the size of the gap between s_hat and s* IS the size of the common-mode error.

OPERATOR FORKS, per BRIEF SS4 rule 0.  Directional hypothesis, so each names its alternative.
NOTE ON RULE 0: these were enumerated by the coordinator, who has a stake in the outcome; no
independent party enumerated them.  That is a weakening of the rule and is recorded, not hidden.

    functional     DECLARED the shipped Bayes-risk distogram score and the shipped top-75 gate,
                   identical to L2's arm.  NOT TAKEN any re-scored or re-sized pool.
    basis          DECLARED point cloud in the medoid frame -- the frame the average is ACTUALLY
                   taken in, so d_k is the dispersion the shipped operator really averages over.
                   NOT TAKEN individually superposing each member onto the native (an oracle
                   frame, which would shrink d_k by fitting m extra rotations).
    readout        DECLARED a single global multiplier about the centroid, as in L2.
                   NOT TAKEN scaling about the medoid.
    normalisation  DECLARED sums of squares over all n residues x 3 coordinates, so |ebar|^2/n is
                   exactly RMSD_c^2 and the identity is checkable.  NOT TAKEN per-residue means.
    null           DECLARED s = 1.0, the incumbent, on identical data -- does the derived
                   estimator beat NOT scaling.  NOT TAKEN the oracle s*, the ceiling under test.
    THE LABEL      DECLARED continuous Ca-RMSD, and separately the common-mode FRACTION
                   f = |ebar|^2 / (|ebar|^2 + mean_k|d_k|^2), the mechanism quantity.
                   NOT TAKEN any binarised "did the estimator win".

  H  the errors of retrieved pool members are substantially i.i.d., so the observable dispersion
     predicts the invisible common mode, s_hat approximates s*, and a derived native-free scale
     captures a real share of the -0.3403 ceiling.
  Falsifier  s_hat failing to beat s = 1.0 past its own MDE with a fold-clustered CI excluding
     zero.  If it fires, the i.i.d. model is refuted at the same time, and the measured f says by
     how much: f >> 1/m is exactly the statement that pool errors are dominated by a shared bias.
  Prediction under the null model  f = 1/m = 1/75 = 0.0133.
  Oracle arms are labelled ORACLE and are EVALUATION ONLY; nothing here is chosen with them.
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


def _save(o, name="errdecomp.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _kabsch_R(P, Q):
    """Rotation R with (P - Pbar) @ R.T ~ (Q - Qbar).  Rotation only; scale is NOT fitted."""
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    return Vt.T @ D @ U.T


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, m = %d" % (len(tg), TOPM), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]; nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        o = np.argsort(sc, kind="stable")
        top = W[o[:TOPM]]
        m = len(top)

        #: the frame the shipped average is ACTUALLY taken in
        P = I.pairwise_rmsd(top); b = I.medoid(P)
        Wm = I.superpose_batch(top, top[b])
        C = Wm.mean(0)

        #: ONE rigid transform for the whole ensemble: the Kabsch fit of C onto the native,
        #: applied identically to every member, so within-ensemble geometry is untouched.
        R = _kabsch_R(C, nat)
        cbar = C.mean(0); tbar = nat.mean(0)
        c_al = (C - cbar) @ R.T
        W_al = (Wm - cbar) @ R.T
        t_al = nat - tbar
        n_res = len(nat)

        ebar = c_al - t_al                       # COMMON error (oracle, eval only)
        d = W_al - c_al[None]                    # IDIOSYNCRATIC error (fully observable)
        S_common = float((ebar ** 2).sum())
        S_idio = float((d ** 2).sum(axis=(1, 2)).mean())
        Cc = float((c_al ** 2).sum())             # |c|^2
        Ct = float((c_al * t_al).sum())           # <c,t>

        s_star = Ct / Cc                          # ORACLE, closed form
        s_hat = 1.0 - S_idio / ((m - 1) * Cc)     # DERIVED, native-free, parameter-free

        rows.append({
            "pdb": pdb, "n": int(n_res), "fold": int(u["fold"]), "m": int(m),
            "rmsd_1": float(I.ca_rmsd(C, nat)),
            "rmsd_shat": float(I.ca_rmsd(c_al * s_hat, t_al)),
            "rmsd_star": float(I.ca_rmsd(c_al * s_star, t_al)),
            "s_star": float(s_star), "s_hat": float(s_hat),
            "S_common": S_common, "S_idio": S_idio, "Cc": Cc, "Ct": Ct,
            "rmsd_c_check": float(np.sqrt(S_common / n_res)),
            "f_common": float(S_common / (S_common + S_idio)),
            "iid_pred_common": float(S_idio / (m - 1)),
        })
        if (c_i + 1) % 20 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("s_star", "s_hat", "S_common", "S_idio", "f_common", "rmsd_shat")
    ok = len(rows) == len(tg) and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "topm": TOPM})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "errdecomp.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    fold = g("fold"); base = g("rmsd_1"); rng = SD.stable_rng("errdecomp", "rep")

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        bb = x[rng.integers(0, k, size=(4000, k))].mean(1)
        F = sorted(set(fold.astype(int)))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se,
                float(np.percentile(bb, 2.5)), float(np.percentile(bb, 97.5)),
                float(np.percentile(fs, 2.5)), float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    m = float(g("m").mean())
    print("\nn = %d.  All arms POINT CLOUDS, medoid frame, m = %.0f.\n" % (len(rows), m))

    print("  IDENTITY CHECK (the decomposition is exact, not approximate):")
    print("    max |RMSD(c,t) - sqrt(S_common/n)| = %.3e A"
          % np.abs(base - g("rmsd_c_check")).max())

    print("\n  ERROR DECOMPOSITION, per target, summed over all n x 3 coordinates:")
    print("    COMMON        |ebar|^2         mean %10.2f   (INVISIBLE from inside the pool)"
          % g("S_common").mean())
    print("    IDIOSYNCRATIC mean_k |d_k|^2   mean %10.2f   (fully observable)"
          % g("S_idio").mean())
    f = g("f_common")
    print("    common-mode FRACTION f          mean %.4f   median %.4f   10th %.4f  90th %.4f"
          % (f.mean(), np.median(f), np.percentile(f, 10), np.percentile(f, 90)))
    print("    i.i.d. model predicts f = 1/m = %.4f   -> OBSERVED / PREDICTED = %.1fx"
          % (1.0 / m, f.mean() * m))
    print("    targets with f > 1/m:  %d / %d" % (int((f > 1.0 / m).sum()), len(rows)))
    print("    the i.i.d. estimate of the common term is short by a factor of %.1f on average"
          % (g("S_common") / g("iid_pred_common")).mean())

    print("\n  SCALE ESTIMATORS:")
    print("    s* (ORACLE, closed form)     mean %.4f  sd %.4f   [%.3f, %.3f]"
          % (g("s_star").mean(), g("s_star").std(), g("s_star").min(), g("s_star").max()))
    print("    s_hat (DERIVED, native-free) mean %.4f  sd %.4f   [%.3f, %.3f]"
          % (g("s_hat").mean(), g("s_hat").std(), g("s_hat").min(), g("s_hat").max()))
    sp = np.corrcoef(np.argsort(np.argsort(g("s_hat"))), np.argsort(np.argsort(g("s_star"))))[0, 1]
    print("    Spearman rho(s_hat, s*) = %+.3f     Pearson %+.3f"
          % (sp, np.corrcoef(g("s_hat"), g("s_star"))[0, 1]))

    print("\n  PRIMARY (paired vs the incumbent, s = 1.0, mean %.4f):" % base.mean())
    for lab, x in (("s_hat  DERIVED native-free  <-- the deployable arm", g("rmsd_shat") - base),
                   ("s*     ORACLE closed form (CEILING)", g("rmsd_star") - base)):
        mm, se, mde, lo, hi, flo, fhi, w = st(x)
        verdict = "BEATS" if (fhi < 0 and abs(mm) > mde) else "not measured"
        print("    %-46s%+.4f SE %.4f MDE %.3f fold[%+.4f,%+.4f] %3dW/%3dL  %s"
              % (lab, mm, se, mde, flo, fhi, w, len(rows) - w, verdict))

    cap = (g("rmsd_shat") - base).mean() / (g("rmsd_star") - base).mean()
    print("\n  fraction of the closed-form ceiling captured by the derived estimator: %.1f%%"
          % (100 * cap))
    print("  Falsifier: s_hat failing to beat s=1.0 past its own MDE with a fold CI excluding zero.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
