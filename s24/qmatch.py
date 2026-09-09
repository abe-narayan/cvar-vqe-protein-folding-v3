"""s24/qmatch.py -- A QUALITY-MATCHED SECOND SOURCE, BUILT FOR FREE, TO SHARPEN THE GENERATOR SPEC.

WHY.  `s24/biasalign.py` (L2) established two things that point opposite ways: the blind-library
source's bias is genuinely non-parallel to retrieval's (cosine 0.647 against a within-source control
of 0.933, and the mixture curve bows 0.226 A below its own endpoints), yet every mixture LOSES
because that source is 0.76 A worse standalone.  The derived spec says a second source must come in
under ~1.27x the incumbent's error to break even.

**That spec has never been tested against an actual source that meets it, because no such source
exists yet -- Lanes B and C are building one.**  Waiting to find out is the expensive way.  A
quality-matched, retrieval-free source can be constructed today at zero training cost:

    B'  =  top-75 by the SHIPPED distogram score, drawn from 2000 uniform library windows
           with the shipped BLOSUM pool EXCLUDED

B' has retrieval-free provenance -- no BLOSUM sequence matching touched it -- but it is score-
selected, so its quality is far above a blind draw.  It is the closest zero-cost stand-in for
"a generator that produces candidates the scorer likes", which is exactly what Lanes B and C will
deliver if they succeed.  If a mixture with B' wins, the spec is achievable and the sprint's main
plan has a demonstrated path.  If a mixture with B' still loses despite meeting the spec, then the
spec derived in L2 is too generous and both generator lanes need a harder target -- which is far
better learned now than after two training runs.

OPERATOR FORKS, per rule 0.  NOTE: enumerated by the coordinator, who has a stake.  Lane E has been
asked to re-enumerate these independently; that is the mitigation, and it is recorded, not hidden.

    functional     DECLARED the shipped Bayes-risk distogram score, unchanged, for BOTH A and B', so
                   the two sources differ in PROVENANCE and not in how they were scored.
                   NOT TAKEN a different or re-tuned functional for B'.
    basis          DECLARED point cloud.  Bias vectors in the NATIVE frame (the only frame in which
                   two sources' errors are commensurable); each arm's RMSD in its own medoid frame,
                   as deployed.  NOT TAKEN a shared frame for RMSD, NOT TAKEN a built chain.
    readout        DECLARED uniform mean of exactly 75 members for every arm, size matched across
                   the whole ladder.  NOT TAKEN letting any arm keep more members.
    normalisation  DECLARED a 2000-window draw for B''s candidate set, fixed a priori at 4x the
                   shipped m and comparable to the shipped K=500's own selectivity.
                   NOT TAKEN sweeping the draw size, which would tune B''s quality against the label.
    null           DECLARED the incumbent A at 75/0 -- does adding B' beat NOT adding it.  The
                   library-uniform source of L2 remains the zero-information reference and is
                   re-reported here so the two ladders are readable side by side.
                   NOT TAKEN an oracle-selected B', which would beg the question.
    THE LABEL      DECLARED continuous Ca-RMSD, plus the bias cosine as a separate reported quantity.
                   NOT TAKEN any binarised win rate, and NOT TAKEN fcommon, which L1 showed selects
                   backwards.

  H  a retrieval-free source at matched quality has non-parallel bias, so some interior mixture BEATS
     the incumbent.
  Falsifier  every interior mixture at or above the incumbent past its own MDE with a fold-clustered
     CI excluding zero.  If it fires, provenance alone does not buy decorrelation at usable quality,
     and the generator lanes need a target harder than "meets the L2 spec".

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

TOPM = 75
DRAW = 2000
MIX = [(75, 0), (60, 15), (50, 25), (38, 37), (25, 50), (0, 75)]


def _save(o, name="qmatch.json"):
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
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    return I.superpose_batch(members, members[b]).mean(0)


def _bias(C, nat):
    R = _kabsch_R(C, nat)
    return (C - C.mean(0)) @ R.T - (nat - nat.mean(0))


def _cos(a, b):
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    return float((a * b).sum() / (na * nb)) if na > 0 and nb > 0 else float("nan")


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, m=%d, B' drawn from %d library windows" % (len(tg), TOPM, DRAW), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        Wall = np.asarray(u["W"], float); nat = np.asarray(u["nat_ca"], float)
        N = len(Wall); n_res = int(u["n"])
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(n_res)

        idx = I.pool_idx(u); Wpool = Wall[idx]
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(Wpool, i, j)), float)
        A = Wpool[np.argsort(sc, kind="stable")[:TOPM]]

        #: B' -- retrieval-free provenance, score-selected quality.  The shipped pool is EXCLUDED so
        #: no BLOSUM-retrieved window can enter B' by accident.
        rng = SD.stable_rng("qmatch", pdb)
        mask = np.ones(N, bool); mask[idx] = False
        pool_free = np.flatnonzero(mask)
        pick = rng.choice(pool_free, min(DRAW, len(pool_free)), replace=False)
        Wf = Wall[pick]
        scf = np.asarray(I.shipped_score(dg, I.pair_dists(Wf, i, j)), float)
        B = Wf[np.argsort(scf, kind="stable")[:TOPM]]

        eA = _bias(_avg(A), nat); eB = _bias(_avg(B), nat)
        mixes = {}
        for a, b in MIX:
            sel = B[:b] if a == 0 else (np.concatenate([A[:a], B[:b]], axis=0) if b else A[:a])
            mixes["m%d_%d" % (a, b)] = float(I.ca_rmsd(_avg(sel), nat))

        rows.append({
            "pdb": pdb, "n": n_res, "fold": int(u["fold"]),
            "rmsd_A": float(I.ca_rmsd(_avg(A), nat)),
            "rmsd_B": float(I.ca_rmsd(_avg(B), nat)),
            "cos_A_B": _cos(eA, eB),
            "q": float(np.linalg.norm(eB) / np.linalg.norm(eA)),
            "mix": mixes,
        })
        if (c_i + 1) % 10 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("rmsd_A", "rmsd_B", "cos_A_B", "q", "mix")
    ok = len(rows) == len(tg) and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "topm": TOPM, "draw": DRAW, "mix": MIX})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "qmatch.json")))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    rng = SD.stable_rng("qmatch", "rep")
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        F = sorted(set(fold))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se, float(np.percentile(fs, 2.5)),
                float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    print("\nn = %d.  A = shipped top-75.  B' = top-75 by the SAME score from %d retrieval-free"
          " library windows.\n" % (len(rows), DRAW))
    print("  A  (incumbent)            %.4f" % g("rmsd_A").mean())
    print("  B' (retrieval-free)       %.4f      q = |e_B|/|e_A| mean %.3f  median %.3f"
          % (g("rmsd_B").mean(), g("q").mean(), np.median(g("q"))))
    print("  cos(bias_A, bias_B')      %+.4f  median %+.4f     [L2 control: within-source +0.9330]"
          % (g("cos_A_B").mean(), np.median(g("cos_A_B"))))
    print("  SPEC CHECK: B' meets the L2 quality bar (q < 1.274)? %s   bias bar (cos <= 0.65)? %s"
          % ("YES" if np.median(g("q")) < 1.274 else "NO",
             "YES" if g("cos_A_B").mean() <= 0.65 else "NO"))

    keys = list(rows[0]["mix"])
    M = {k: np.array([r["mix"][k] for r in rows], float) for k in keys}
    e0, e1 = M[keys[0]], M[keys[-1]]
    print("\n  MATCHED-SIZE MIXTURES  (A / B', always 75 members)")
    print("    %-14s%9s%12s%12s" % ("mixture", "RMSD", "vs 75/0", "vs the line"))
    best = None
    for k in keys:
        frac = float(k.split("_")[1]) / TOPM
        line = (1 - frac) * e0 + frac * e1
        d = M[k] - e0
        m, se, mde, flo, fhi, w = st(d)
        v = ("BEATS" if (fhi < 0 and abs(m) > mde) else
             "worse" if (flo > 0 and abs(m) > mde) else "ns")
        print("    %-14s%9.4f%+12.4f%+12.4f   SE %.4f MDE %.4f %3dW/%3dL %s"
              % (k, M[k].mean(), m, (M[k] - line).mean(), se, mde, w, len(rows) - w, v))
        if k != keys[0] and (best is None or M[k].mean() < M[best].mean()):
            best = k
    print("\n  best interior/endpoint arm: %s at %.4f  (incumbent %.4f)"
          % (best, M[best].mean(), e0.mean()))
    print("  Falsifier: every interior mixture at or above the incumbent past its own MDE.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
