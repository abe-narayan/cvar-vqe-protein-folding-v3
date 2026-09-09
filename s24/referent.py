"""s24/referent.py -- IS THE DISTOGRAM THE SHARED REFERENT?  DOES THE PIPELINE INHERIT ITS PRIOR'S
ERROR?  THE MEASUREMENT L3 DEMANDS.

WHY.  L3 measured that a retrieval-free, quality-matched candidate source selected by the SAME
shipped distogram score has bias cosine +0.9432 against the incumbent -- above the +0.9330
within-source control -- and a mixture curve that is flat against linear interpolation.  Selection by
the score makes candidate sets parallel; provenance buys nothing.  The natural explanation is the
project's own `shared-referent-floor`: two quantities measured against a COMMON reference correlate
by construction, and both sources were selected against the same distogram prediction.

If that is right, the 68% common-mode error of s23 L9 is not a property of the candidates at all.
It is the prior's error, transmitted through the score into whatever candidates the score is given.

THIS IS MEASURABLE IN DISTANCE SPACE, WITH NO STRUCTURE OPERATIONS AND NO MODEL.  For each target,
over the same pair set the score uses:

    Dhat  the per-pair Bayes-optimal distance -- literally what the shipped score wants, taken as
          grid[argmin(risk)] rather than the posterior mean, because argmin(risk) is the point the
          score actually pulls candidates toward.  `expected` is reported as a declared secondary.
    Dt    the native pair distances
    Dc    the pair distances of the emitted cloud (shipped top-75 uniform coordinate average)

  PRIMARY 1, and it needs no correlation machinery to interpret:
        RMS|Dc - Dhat|   against   RMS|Dt - Dhat|
      **Is the emitted structure closer to the prediction than the TRUTH is?**  If it is, the
      pipeline is over-fitting the prior and the prior's error is being written into the output.

  PRIMARY 2, the inheritance coefficient:
        beta = <Dc - Dt, Dhat - Dt> / |Dhat - Dt|^2
      the fraction of the prior's own error that reappears in the emitted structure.  beta = 1 is
      total inheritance; beta = 0 is none.

THE CONTROL, WHICH IS THE WHOLE DIFFICULTY.  Both `Dc - Dt` and `Dhat - Dt` are measured against the
same referent Dt, and subtracting a common term from two independent quantities induces positive
covariance by construction (Cov(X-Z, Y-Z) = Var(Z) for independent X, Y, Z).  Project memory records
this exact trap turning a "2/3 sequence-independent" claim into "1/5".  So every number above is
reported beside a PLACEBO computed identically with a MISMATCHED same-length target's distogram:
same construction, same referent, no true relationship.  **The placebo is the floor, and only the
excess over the placebo is evidence.**

OPERATOR FORKS, per rule 0.  NOTE: enumerated by the coordinator, who has a stake; Lane E has been
asked to re-enumerate.  Recorded, not hidden.

    functional     DECLARED Dhat = grid[argmin(risk)], the score's own per-pair optimum, because the
                   question is what the SCORE pulls toward.  NOT TAKEN the posterior mean `expected`
                   (reported as a secondary, and any disagreement between them reported).
    basis          DECLARED distance space over the score's own pair set (min_sep=2), where the prior
                   and the output are directly commensurable.  NOT TAKEN Cartesian space, which would
                   require a common frame the prior does not have.
    readout        DECLARED the shipped top-75 uniform coordinate average, as deployed.
                   NOT TAKEN the medoid or any other aggregation.
    normalisation  DECLARED RMS over pairs for the distances and an unnormalised projection for beta.
                   NOT TAKEN per-pair z-scoring by the distogram's own sd, which would let a
                   confident-and-wrong prior hide inside its own confidence.
    null           DECLARED the mismatched same-length distogram placebo, computed through the
                   identical pipeline.  NOT TAKEN zero, which ignores the shared-referent floor.
    THE LABEL      DECLARED both the closeness comparison and beta, reported together, because
                   either alone can be read the wrong way.  NOT TAKEN a correlation coefficient,
                   which hides the magnitude.

  H  the emitted structure has inherited a substantial share of the prior's error: beta well above
     its placebo, and the cloud closer to Dhat than the native is.
  Falsifier  beta at or below its placebo, and the cloud no closer to Dhat than the native is.  Then
     the prior is NOT the shared referent and L3's alignment needs a different explanation.

  Natives are read for EVALUATION ONLY.  Nothing here selects anything.
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


def _save(o, name="referent.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _avg(members):
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    return I.superpose_batch(members, members[b]).mean(0)


def _dhat(dg):
    """The score's own per-pair optimum: the distance minimising the shipped Bayes risk."""
    return np.asarray(dg["grid"], float)[np.argmin(np.asarray(dg["risk"], float), axis=1)]


def _stats(Dc, Dt, Dh):
    eP = Dh - Dt                     # the prior's own error
    eC = Dc - Dt                     # the emitted structure's error
    den = float((eP * eP).sum())
    return {
        "rms_cloud_to_hat": float(np.sqrt(((Dc - Dh) ** 2).mean())),
        "rms_nat_to_hat": float(np.sqrt(((Dt - Dh) ** 2).mean())),
        "rms_cloud_to_nat": float(np.sqrt((eC ** 2).mean())),
        "beta": float((eC * eP).sum() / den) if den > 0 else float("nan"),
        "cos": float((eC * eP).sum() / np.sqrt((eC * eC).sum() * den)) if den > 0 else float("nan"),
    }


def run():
    tg = I.targets()
    by_len = {}
    for t in tg:
        by_len.setdefault(int(t["n"]), []).append(t["pdb"])
    rows = []
    print("targets: %d" % len(tg), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        Wall = np.asarray(u["W"], float); nat = np.asarray(u["nat_ca"], float)
        n_res = int(u["n"]); i, j = I.pair_index(n_res)
        dg = I.distogram(pdb, u["seq"], u["fold"])
        idx = I.pool_idx(u); Wpool = Wall[idx]
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(Wpool, i, j)), float)
        A = Wpool[np.argsort(sc, kind="stable")[:TOPM]]
        C = _avg(A)

        Dt = I.pair_dists(nat[None], i, j)[0]
        Dc = I.pair_dists(C[None], i, j)[0]
        Dh = _dhat(dg)
        Dexp = np.asarray(dg["expected"], float)

        row = {"pdb": pdb, "n": n_res, "fold": int(u["fold"]),
               "rmsd": float(I.ca_rmsd(C, nat))}
        row["real"] = _stats(Dc, Dt, Dh)
        row["real_expected"] = _stats(Dc, Dt, Dexp)          # declared secondary functional

        #: --- PLACEBO: a mismatched same-length target's distogram, identical construction
        peers = [p for p in by_len[n_res] if p != pdb]
        if peers:
            rng = SD.stable_rng("referent", pdb)
            other = peers[int(rng.integers(0, len(peers)))]
            ou = I.load_univ(other)
            odg = I.distogram(other, ou["seq"], ou["fold"])
            row["placebo_pdb"] = other
            row["placebo"] = _stats(Dc, Dt, _dhat(odg))
        rows.append(row)
        if (c_i + 1) % 20 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("real", "real_expected", "placebo", "rmsd")
    ok = len(rows) == len(tg) and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "topm": TOPM})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "referent.json")))["rows"]
    rows = [r for r in rows if "placebo" in r]
    fold = np.array([r["fold"] for r in rows], int)
    rng = SD.stable_rng("referent", "rep")
    G = lambda a, k: np.array([r[a][k] for r in rows], float)      # noqa: E731

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        F = sorted(set(fold))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se, float(np.percentile(fs, 2.5)),
                float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    print("\nn = %d (targets with a same-length peer for the placebo).  Distance space, min_sep=2.\n"
          % len(rows))
    print("  PRIMARY 1 -- is the emitted cloud CLOSER to the prediction than the NATIVE is?")
    print("    RMS |Dc - Dhat|   (cloud to prediction)   %.4f A" % G("real", "rms_cloud_to_hat").mean())
    print("    RMS |Dt - Dhat|   (native to prediction)  %.4f A" % G("real", "rms_nat_to_hat").mean())
    m, se, mde, flo, fhi, w = st(G("real", "rms_cloud_to_hat") - G("real", "rms_nat_to_hat"))
    v = ("CLOSER than the native" if (fhi < 0 and abs(m) > mde) else
         "further" if (flo > 0 and abs(m) > mde) else "not measured")
    print("    difference  %+.4f  SE %.4f  MDE %.4f  fold[%+.4f,%+.4f]  %3dW/%3dL  -> %s"
          % (m, se, mde, flo, fhi, w, len(rows) - w, v))
    print("    RMS |Dc - Dt|     (cloud to native)       %.4f A" % G("real", "rms_cloud_to_nat").mean())

    print("\n  PRIMARY 2 -- beta, the share of the PRIOR's error that reappears in the OUTPUT")
    for a, lab in (("real", "REAL   distogram"), ("placebo", "PLACEBO mismatched same-length"),
                   ("real_expected", "real, posterior-mean functional (secondary)")):
        b = G(a, "beta"); c = G(a, "cos")
        print("    %-44s beta %+.4f  median %+.4f   cos %+.4f"
              % (lab, b.mean(), np.median(b), c.mean()))
    m, se, mde, flo, fhi, w = st(G("real", "beta") - G("placebo", "beta"))
    v = ("EXCESS OVER THE FLOOR" if (flo > 0 and abs(m) > mde) else
         "below the floor" if (fhi < 0 and abs(m) > mde) else "not measured")
    print("    real - placebo  %+.4f  SE %.4f  MDE %.4f  fold[%+.4f,%+.4f]  %3dW/%3dL  -> %s"
          % (m, se, mde, flo, fhi, len(rows) - w, w, v))
    print("\n  The placebo is the shared-referent FLOOR.  Only the excess over it is evidence.")
    print("  Falsifier: beta at or below its placebo AND the cloud no closer to Dhat than the native.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
