"""s24/biasalign.py -- ARE TWO CANDIDATE SOURCES WRONG IN THE SAME DIRECTION?  AND DOES MIXING THEM
HELP?  THE ZERO-COST PROTOTYPE OF THE DIRECTIVE'S CENTRAL CAUSAL EXPERIMENT.

WHY THIS RUNS BEFORE A MODEL IS TRAINED.  The Sprint 24 directive's whole plan rests on one premise:
that a second candidate source can supply structural information whose ERROR IS LESS CORRELATED with
the retrieval pool's.  `s24/srcdecomp.py` has already established the sources' marginal statistics.
It cannot answer the premise, because the premise is about the ANGLE BETWEEN their error vectors,
not about either one's magnitude.

**That angle is measurable today, with no neural network, using a source that already exists.**  A
uniform draw from the same universe of same-length windows is a genuine second source: real protein
geometry, zero retrieval information, and -- critically -- drawn from the SAME corpus a
library-trained generator would be trained on.  So the angle between its bias and retrieval's bias
is an upper bound on how decorrelated a corpus-trained generator can be expected to be.  If two
sources that share only the corpus are already nearly parallel, then the corpus is the bias, and a
model trained on it inherits the bias no matter how it is parameterised.

THE CONTROL THAT MAKES THE ANGLE MEAN ANYTHING.  Two averages of real protein windows will point in
similar directions for trivial reasons -- both are contracted, both are protein-like.  So a raw
cosine is uninterpretable on its own.  The comparison that matters is against the WITHIN-SOURCE
ceiling: the cosine between two INDEPENDENT draws from the SAME source, which is what perfect bias
sharing looks like on this instrument.

    cos(e_A, e_B)  /  cos(e_C1, e_C2)   ->  1  means "as aligned as a source is with itself"
                                        ->  0  means genuinely independent bias

WHAT THIS FILE MEASURES

  (1) bias angles, in the native frame, between
        A = top-75 by the shipped distogram score            (the incumbent)
        B = top-75 by BLOSUM `order` alone                   (retrieval, no score)
        C = 75 uniform draws from the universe               (corpus, no retrieval)
        C' = a second, independent uniform draw               (the WITHIN-SOURCE control)

  (2) matched-size MIXTURES.  Every arm emits exactly 75 members through the identical readout, so
      the arms differ in composition and in nothing else:
          75/0, 60/15, 50/25, 38/37, 25/50, 0/75   (score-ranked / library-uniform)
      Under the error model, if the two sources' biases were independent, a mixture must beat both
      parents somewhere in the interior.  If every interior point is a straight-line interpolation of
      the endpoints, the biases are parallel and mixing buys nothing -- which is the same statement
      as (1), arrived at through the endpoint that actually matters, RMSD.

  (3) the directive's SS17 union, run honestly: pool the shipped K=500 with 500 uniform library
      windows, score ALL 1000 with the SAME shipped functional, take the SAME top-75.  This asks
      whether the existing scorer, offered genuinely new candidates, picks them up -- and the share
      of the merged top-75 drawn from the non-retrieved half prices what BLOSUM retrieval is
      contributing over the distogram score alone.

OPERATOR FORKS, per the standing rule 0.  Each names the alternative not taken.
NOTE: enumerated by the coordinator, who has a stake.  Recorded, not hidden.

    functional     DECLARED the shipped Bayes-risk distogram score, unchanged, for A and for the SS17
                   union.  NOT TAKEN a re-tuned or re-normalised score for the enlarged pool.
    basis          DECLARED point cloud throughout.  Bias vectors are compared in the NATIVE frame --
                   the only frame in which two different sources' errors are commensurable -- and
                   each arm's RMSD is still computed in its own medoid frame, as deployed.
                   NOT TAKEN comparing bias vectors in either source's own medoid frame, which would
                   make the angle a function of the frame rather than of the error.
    readout        DECLARED uniform mean of exactly 75 members for every mixture arm, so size is
                   matched across the whole ladder.  NOT TAKEN letting the union arm keep more members.
    normalisation  DECLARED cosine of the raw bias vectors, and separately the ratio to the
                   within-source control.  NOT TAKEN correlation of coordinates, which would remove
                   the mean displacement that IS the bias.
    null           DECLARED the within-source cosine cos(e_C1, e_C2) -- two independent draws from one
                   source, the value that "shares all of its bias" takes on this instrument.
                   NOT TAKEN zero, which is the wrong null for structures that are all protein-like.
    THE LABEL      DECLARED continuous Ca-RMSD for the mixtures and the raw cosine for the angles.
                   NOT TAKEN any binarised "is it decorrelated".

  H  a corpus-drawn second source has bias substantially unaligned with retrieval's, so some interior
     mixture beats the pure incumbent.
  Falsifier  every interior mixture lying at or above the straight line between the endpoints, AND
     cos(e_A, e_C) reaching the within-source control.  That combination says the corpus carries the
     bias, and it predicts -- before any training run -- that a generator trained on this corpus will
     not reduce common-mode error.

  Native structures are read for EVALUATION ONLY, and only after every candidate set is fixed.
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
MIX = [(75, 0), (60, 15), (50, 25), (38, 37), (25, 50), (0, 75)]


def _save(o, name="biasalign.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _kabsch_R(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def _avg(members):
    """The shipped operator: superpose on the medoid, uniform mean.  Returns the point cloud."""
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    return I.superpose_batch(members, members[b]).mean(0)


def _bias(C, nat):
    """The error vector of an emitted cloud, expressed in the NATIVE frame.

    Rotation and translation only -- exactly what Kabsch fits and exactly what RMSD quotients out --
    so what remains is the part of the error that RMSD actually charges for.
    """
    R = _kabsch_R(C, nat)
    c_al = (C - C.mean(0)) @ R.T
    t_al = nat - nat.mean(0)
    return c_al - t_al


def _cos(a, b):
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    return float((a * b).sum() / (na * nb)) if na > 0 and nb > 0 else float("nan")


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, m = %d, mixtures %s" % (len(tg), TOPM, MIX), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        Wall = np.asarray(u["W"], float); nat = np.asarray(u["nat_ca"], float)
        order = np.asarray(u["order"], int); N = len(Wall)
        n_res = int(u["n"])
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(n_res)

        idx = I.pool_idx(u); Wpool = Wall[idx]
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(Wpool, i, j)), float)
        o = np.argsort(sc, kind="stable")

        rng = SD.stable_rng("biasalign", pdb)
        c1 = rng.choice(N, TOPM, replace=False)
        c2 = rng.choice(N, TOPM, replace=False)          # independent draw: the WITHIN-SOURCE control

        A = Wpool[o[:TOPM]]; B = Wall[order[:TOPM]]
        C1 = Wall[c1]; C2 = Wall[c2]
        eA = _bias(_avg(A), nat); eB = _bias(_avg(B), nat)
        eC1 = _bias(_avg(C1), nat); eC2 = _bias(_avg(C2), nat)

        #: --- matched-size mixtures: same readout, same count, only composition changes
        mixes = {}
        for a, b in MIX:
            sel = np.concatenate([A[:a], C1[:b]], axis=0) if b else A[:a]
            if a == 0:
                sel = C1[:b]
            mixes["m%d_%d" % (a, b)] = float(I.ca_rmsd(_avg(sel), nat))

        #: --- the directive's SS17 union: merge the candidate POOLS, score with the SAME functional
        lib500 = rng.choice(N, min(500, N), replace=False)
        Wu = np.concatenate([Wpool, Wall[lib500]], axis=0)
        src = np.concatenate([np.zeros(len(Wpool), int), np.ones(len(lib500), int)])
        scu = np.asarray(I.shipped_score(dg, I.pair_dists(Wu, i, j)), float)
        ou = np.argsort(scu, kind="stable")[:TOPM]
        union_rmsd = float(I.ca_rmsd(_avg(Wu[ou]), nat))

        rows.append({
            "pdb": pdb, "n": n_res, "fold": int(u["fold"]), "n_universe": int(N),
            "cos_A_B": _cos(eA, eB), "cos_A_C": _cos(eA, eC1), "cos_B_C": _cos(eB, eC1),
            "cos_C1_C2": _cos(eC1, eC2),
            "norm_A": float(np.linalg.norm(eA)), "norm_C": float(np.linalg.norm(eC1)),
            "mix": mixes,
            "union_rmsd": union_rmsd,
            "union_lib_frac": float(src[ou].mean()),
            "rmsd_A": float(I.ca_rmsd(_avg(A), nat)),
        })
        if (c_i + 1) % 10 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("cos_A_C", "cos_C1_C2", "mix", "union_rmsd", "union_lib_frac")
    ok = len(rows) == len(tg) and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "topm": TOPM, "mix": MIX})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "biasalign.json")))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    rng = SD.stable_rng("biasalign", "rep")
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        F = sorted(set(fold))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se,
                float(np.percentile(fs, 2.5)), float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    print("\nn = %d.  Bias vectors in the NATIVE frame; RMSD in each arm's own medoid frame.\n"
          % len(rows))
    print("  BIAS-DIRECTION COSINES")
    ctl = g("cos_C1_C2")
    for k, lab in (("cos_A_B", "score-75  vs  BLOSUM-75      (both retrieved)"),
                   ("cos_A_C", "score-75  vs  library-75     <-- THE QUESTION"),
                   ("cos_B_C", "BLOSUM-75 vs  library-75"),
                   ("cos_C1_C2", "library-75 vs library-75'    <-- WITHIN-SOURCE CONTROL")):
        v = g(k)
        print("    %-46s mean %+.4f  median %+.4f  sd %.4f  [%+.3f,%+.3f]"
              % (lab, v.mean(), np.median(v), v.std(), np.percentile(v, 10), np.percentile(v, 90)))
    print("\n    ratio cos(A,C) / cos(C1,C2) = %.3f"
          % (g("cos_A_C").mean() / ctl.mean()))
    print("    -> 1.0 means the retrieved pool and a blind corpus draw are as aligned as the corpus")
    print("       is with ITSELF, i.e. the bias is a property of the CORPUS, not of the retrieval rule.")

    print("\n  MATCHED-SIZE MIXTURES (score-ranked / library-uniform, always 75 members)")
    keys = [k for k in rows[0]["mix"]]
    M = {k: np.array([r["mix"][k] for r in rows], float) for k in keys}
    e0, e1 = M[keys[0]], M[keys[-1]]
    print("    %-14s%9s%12s%12s" % ("mixture", "RMSD", "vs 75/0", "vs the line"))
    for idx_k, k in enumerate(keys):
        frac = float(k.split("_")[1]) / TOPM
        line = (1 - frac) * e0 + frac * e1          # the parallel-bias prediction
        d = M[k] - e0
        m, se, mde, flo, fhi, w = st(d)
        dl = (M[k] - line).mean()
        v = ("BEATS" if (fhi < 0 and abs(m) > mde) else
             "worse" if (flo > 0 and abs(m) > mde) else "ns")
        print("    %-14s%9.4f%+12.4f%+12.4f   SE %.4f MDE %.4f %3dW/%3dL %s"
              % (k, M[k].mean(), m, dl, se, mde, w, len(rows) - w, v))
    print("    'vs the line' is the departure from a straight interpolation of the two endpoints.")
    print("    Independent biases MUST bow below the line; parallel biases sit on it.")

    print("\n  DIRECTIVE SS17 UNION: shipped 500 + 500 uniform library, SAME score, SAME top-75")
    m, se, mde, flo, fhi, w = st(g("union_rmsd") - g("rmsd_A"))
    v = ("BEATS" if (fhi < 0 and abs(m) > mde) else
         "worse" if (flo > 0 and abs(m) > mde) else "not measured")
    print("    union vs incumbent   %+.4f  SE %.4f  MDE %.4f  fold[%+.4f,%+.4f]  %3dW/%3dL  %s"
          % (m, se, mde, flo, fhi, w, len(rows) - w, v))
    lf = g("union_lib_frac")
    print("    share of the merged top-75 drawn from the NON-RETRIEVED half: mean %.3f  median %.3f"
          % (lf.mean(), np.median(lf)))
    print("    (0.0 = the score never prefers a blind library window over a retrieved one;")
    print("     0.5 = retrieval contributes nothing the score cannot recover on its own.)")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
