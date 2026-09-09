"""s24/srcdecomp.py -- WHERE DOES THE 68% COMMON-MODE ERROR COME FROM?  LIBRARY, RETRIEVAL RULE,
OR SCORE?  THE TIER-1 SCREEN THAT DECIDES WHETHER A LEARNED GENERATOR CAN HELP AT ALL.

WHY THIS RUNS BEFORE ANY MODEL IS BUILT.  Sprint 23 (L9) measured, as an exact identity, that 67.6%
of the shipped top-75's squared error is a bias its members hold IN COMMON, against 1.3% under an
i.i.d. model.  The Sprint 24 directive's central hypothesis is that a learned generator can produce
candidates whose errors are LESS correlated.  **But that hypothesis is not one hypothesis, it is
three, and they have opposite consequences for what to build:**

  (H-lib)    the shared bias is a property of the FRAGMENT LIBRARY -- the corpus of real protein
             windows simply does not contain peptide-like geometry, so every window drawn from it is
             wrong in the same direction.  A generator TRAINED ON THAT LIBRARY inherits the bias and
             cannot remove it.  The fix would have to be a different corpus or a genuinely
             extrapolating model, not a resampling of the same manifold.

  (H-retr)   the shared bias is introduced by the BLOSUM RETRIEVAL RULE -- sequence matching selects
             a correlated subset of an otherwise diverse library.  A generator trained on the library
             at large CAN remove it, because the information is present in the corpus and only the
             selection rule is discarding it.  This is the case in which the directive's plan works.

  (H-score)  the shared bias is introduced by the DISTOGRAM SCORE gate on top of retrieval.  Then the
             cheapest fix is upstream of any neural network at all.

These are separated by ONE measurement, and it costs one pass over the instrument.  Same readout,
same m, same basis, same frame -- only the SOURCE of the 75 members changes:

    P0  top-75 by the shipped distogram score          = the incumbent            (score + retrieval)
    P1  top-75 by BLOSUM `order` alone                                            (retrieval, no score)
    P2  75 windows drawn UNIFORMLY AT RANDOM from the whole universe              (library, no retrieval)
    P3  the full K=500 shipped pool, averaged                                     (reference)
    P4  all universe windows, averaged                                            (the library's own mean)

`u["W"]` is the FULL universe of same-length windows (e.g. 13,247 for a 14-mer) and `u["order"]` is
the BLOSUM ranking over it, so P2 is a genuine library sample that no retrieval rule has touched.

WHAT EACH OUTCOME MEANS, WRITTEN DOWN BEFORE THE NUMBERS EXIST.

    |e|^2(P2) >= |e|^2(P0) and f(P2) >= f(P0)   -> H-lib.  Retrieval is HELPING; the residual bias is
        the corpus itself.  A generator trained on this corpus is predicted NOT to reduce common-mode
        error, and the directive's Source-B arm should be expected to fail for a reason that is
        knowable NOW rather than after training a model.  The productive variant becomes the RESIDUAL
        generator (directive SS21/SS39), which does not need the corpus to contain the answer.

    |e|^2(P2) <  |e|^2(P0)                      -> H-retr.  The library carries less common error than
        the retrieved subset, i.e. retrieval concentrates bias.  A learned generator over the library
        has real headroom and the directive's main plan is on solid ground.

    f(P1) ~ f(P0) but f(P2) much lower          -> the bias enters at retrieval, not at scoring.
    f(P1) much lower than f(P0)                 -> the bias enters at the SCORE, which is cheap to fix.

OPERATOR FORKS, per the standing rule 0.  Each names the alternative not taken.
NOTE: enumerated by the coordinator, who has a stake.  Recorded, not hidden.

    functional     DECLARED the shipped Bayes-risk distogram score for P0, exactly as deployed, and
                   the shipped BLOSUM `order` for P1.  NOT TAKEN any re-scoring or re-ranking.
    basis          DECLARED point cloud, each arm in its OWN medoid frame -- the frame its own average
                   is taken in.  NOT TAKEN a shared frame across arms, and NOT TAKEN a built chain.
    readout        DECLARED the uniform mean of exactly m=75 members for P0-P2, so the arms differ in
                   SOURCE and in nothing else.  NOT TAKEN matching on anything but count.
    normalisation  DECLARED sums of squares over all n x 3 coordinates, so |ebar|^2/n is exactly
                   RMSD^2 and the identity is checkable per arm.  NOT TAKEN per-residue means.
    null           DECLARED P2, the uniform library draw -- the zero-retrieval-information source.
                   It is a PLAUSIBLE control (real protein windows of the right length), not a
                   degenerate one, per the project's standing rule that a zero-information control
                   must be plausible rather than uniform-on-the-torus.
                   NOT TAKEN a synthetic or randomised-geometry control.
    THE LABEL      DECLARED the common-mode fraction f AND the absolute common term |ebar|^2, because
                   f alone can fall while the absolute error rises.  Both are reported for every arm.
                   NOT TAKEN f alone, and NOT TAKEN any binarised "did the source win".

  Falsifier for the sprint's main plan: if P2 shows f and |ebar|^2 at or above P0's, then the corpus,
  not the retrieval rule, carries the shared bias, and a generator trained on the corpus is predicted
  to inherit it.  That prediction is then testable directly and cheaply.

  Native structures are read for EVALUATION ONLY.  No arm uses native information to select anything.
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
NDRAW = 8          # P2 is stochastic; average the statistic over this many independent draws


def _save(o, name="srcdecomp.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _kabsch_R(P, Q):
    """Rotation R with (P - Pbar) @ R.T ~ (Q - Qbar).  Rotation only; scale is NOT fitted."""
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def _decomp(members, nat):
    """The L9 decomposition for one candidate set, in the frame its own average is taken in.

    Returns rmsd, S_common (|ebar|^2), S_idio (mean_k |d_k|^2), f, s_star, |c|^2.
    """
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    Wm = I.superpose_batch(members, members[b])
    C = Wm.mean(0)
    R = _kabsch_R(C, nat)
    c_al = (C - C.mean(0)) @ R.T
    W_al = (Wm - C.mean(0)) @ R.T
    t_al = nat - nat.mean(0)
    ebar = c_al - t_al
    d = W_al - c_al[None]
    S_common = float((ebar ** 2).sum())
    S_idio = float((d ** 2).sum(axis=(1, 2)).mean())
    Cc = float((c_al ** 2).sum()); Ct = float((c_al * t_al).sum())
    return {
        "rmsd": float(I.ca_rmsd(C, nat)),
        "S_common": S_common, "S_idio": S_idio,
        "f_common": float(S_common / (S_common + S_idio)),
        "member_rmsd": float(np.sqrt((S_common + S_idio) / len(nat))),
        "s_star": float(Ct / Cc), "Cc": Cc,
        "spread": float(P.sum() / (len(P) ** 2 - len(P))),
    }


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, m = %d, %d random-library draws" % (len(tg), TOPM, NDRAW), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        Wall = np.asarray(u["W"], float)
        nat = np.asarray(u["nat_ca"], float)
        order = np.asarray(u["order"], int)
        N = len(Wall)

        idx = I.pool_idx(u)                       # the shipped K=500 BLOSUM pool
        Wpool = Wall[idx]
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(Wpool, i, j)), float)
        o = np.argsort(sc, kind="stable")

        row = {"pdb": pdb, "n": int(u["n"]), "fold": int(u["fold"]), "n_universe": int(N)}
        row["P0_score75"] = _decomp(Wpool[o[:TOPM]], nat)          # incumbent
        row["P1_blosum75"] = _decomp(Wall[order[:TOPM]], nat)      # retrieval, no score
        rng = SD.stable_rng("srcdecomp", pdb)
        draws = [_decomp(Wall[rng.choice(N, TOPM, replace=False)], nat) for _ in range(NDRAW)]
        row["P2_library75"] = {k: float(np.mean([d[k] for d in draws])) for k in draws[0]}
        row["P2_sd"] = {k: float(np.std([d[k] for d in draws], ddof=1)) for k in draws[0]}
        row["P3_pool500"] = _decomp(Wpool, nat)                    # the whole shipped pool
        #: P4 -- the library's own mean.  Subsampled when the universe is large: the average of a
        #: uniform sample IS an unbiased estimate of the library mean, and 2000 is far past the
        #: point where the sample mean has converged for this purpose.
        sub = rng.choice(N, min(N, 300), replace=False)
        row["P4_library_all"] = _decomp(Wall[sub], nat)
        rows.append(row)

        if (c_i + 1) % 10 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("P0_score75", "P1_blosum75", "P2_library75", "P3_pool500", "P4_library_all")
    ok = len(rows) == len(tg) and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "topm": TOPM, "ndraw": NDRAW})
    report(rows)
    return rows


ARMS = [("P0_score75", "P0  top-75 by shipped score   [INCUMBENT]"),
        ("P1_blosum75", "P1  top-75 by BLOSUM order    [retrieval, no score]"),
        ("P2_library75", "P2  75 UNIFORM from universe  [library, no retrieval]"),
        ("P3_pool500", "P3  full K=500 shipped pool   [reference]"),
        ("P4_library_all", "P4  library mean (300 draw)   [reference]")]


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "srcdecomp.json")))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    rng = SD.stable_rng("srcdecomp", "rep")

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        F = sorted(set(fold))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se,
                float(np.percentile(fs, 2.5)), float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    G = lambda a, k: np.array([r[a][k] for r in rows], float)      # noqa: E731
    print("\nn = %d.  Every arm: uniform coordinate average, point cloud, its own medoid frame.\n"
          % len(rows))
    print("  mean universe size: %.0f windows per target\n"
          % np.mean([r["n_universe"] for r in rows]))
    print("  %-46s%9s%9s%12s%12s%9s" % ("source", "RMSD", "member", "|ebar|^2", "idio", "f"))
    for a, lab in ARMS:
        print("  %-46s%9.4f%9.4f%12.2f%12.2f%9.4f"
              % (lab, G(a, "rmsd").mean(), G(a, "member_rmsd").mean(),
                 G(a, "S_common").mean(), G(a, "S_idio").mean(), G(a, "f_common").mean()))

    print("\n  PAIRED CONTRASTS -- the decision this file exists to make:")
    base_f = G("P0_score75", "f_common"); base_c = G("P0_score75", "S_common")
    base_r = G("P0_score75", "rmsd")
    for a, lab in ARMS[1:]:
        for k, bl, unit in (("f_common", base_f, "f"), ("S_common", base_c, "|ebar|^2"),
                            ("rmsd", base_r, "RMSD")):
            m, se, mde, flo, fhi, w = st(G(a, k) - bl)
            v = ("LOWER" if (fhi < 0 and abs(m) > mde) else
                 "HIGHER" if (flo > 0 and abs(m) > mde) else "not measured")
            print("    %-30s %-9s vs P0  %+10.4f SE %8.4f MDE %8.4f fold[%+.4f,%+.4f] %3dW/%3dL  %s"
                  % (lab.split("[")[0].strip(), unit, m, se, mde, flo, fhi, w, len(rows) - w, v))

    print("\n  READING THE VERDICT (declared before the run, see the module docstring):")
    dc = (G("P2_library75", "S_common") - base_c).mean()
    df = (G("P2_library75", "f_common") - base_f).mean()
    if dc >= 0 and df >= 0:
        print("    P2 carries MORE common error than P0 in both absolute and fractional terms.")
        print("    -> H-lib.  Retrieval is HELPING.  The residual shared bias is the CORPUS.")
        print("    -> A generator trained on this corpus is predicted to INHERIT the bias.")
        print("    -> The residual-generator route (directive SS21/SS39) is the one with headroom.")
    elif dc < 0:
        print("    P2 carries LESS absolute common error than P0.")
        print("    -> H-retr.  Retrieval CONCENTRATES bias; a library-trained generator has headroom.")
    else:
        print("    Mixed: fraction and absolute term disagree.  Report both; do not collapse them.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
