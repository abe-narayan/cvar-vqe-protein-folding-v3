#!/usr/bin/env python
"""s32/s32_P_placebo.py -- S32 LANE P: the deployable direction hiding inside the placebo, priced
against its SHARED-REFERENT FLOOR before anyone calls it a lead.

`s32_P_sign.py` found that a MISMATCHED same-length native reproduces 87.1% of pc1's alignment
with the true common mode.  A mismatched deposited structure is available at INFERENCE -- it is
library data, not this target's native -- so `mu_hat = c - t'` is a genuinely DEPLOYABLE estimate
of `mu = c - t`.  Before that is called a lead it has to clear the project's own trap:

    `shared-referent-floor`: two quantities measured against a COMMON reference correlate BY
    CONSTRUCTION.  mu and mu_hat share `c`.  With t = m + xi and t' = m + xi' independent,
    cos(mu, mu_hat) -> |c-m|^2 / (|c-m|^2 + sigma^2) with NO information about t whatsoever.

So this file measures the cos AND the floor, in the same script, and reports the excess.
FLOOR, measured not assumed: the same statistic with t' replaced by a structure drawn from the
same same-length population and, separately, the value predicted by the identity above using the
measured population mean `m`.  ORACLE / NOT DEPLOYABLE (every cos reads the native to score it).

    python s32/s32_P_placebo.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I                                        # noqa: E402
from s24 import stats_lib as ST                                       # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s32_P_placebo.json")
DIM, SEED = 128, 32_0_6613
AGRID = np.linspace(0.0, 3.0, 121)


def main():
    ts = I.targets()
    by_n = {}
    for t in ts:
        by_n.setdefault(int(t["n"]), []).append(t["pdb"])
    rows = []
    for t in ts:
        pdb = t["pdb"]; n = int(t["n"])
        u = I.load_univ(pdb)
        DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(float)
        pool = np.asarray(u["order"], int)[:500]
        W = np.asarray(u["W"], float)[pool]
        top = np.argsort(DIS, kind="stable")[:DIM]
        Wt = W[top]
        Pt = I.pairwise_rmsd(Wt).astype(np.float32).astype(float)
        b = I.medoid(Pt)
        X = I.superpose_batch(Wt, Wt[b]).reshape(DIM, -1)
        c = X.mean(0)
        tf = I.superpose_batch(np.asarray(u["nat_ca"], float)[None], Wt[b])[0].reshape(-1)
        mu = c - tf; nmu = float(np.linalg.norm(mu))

        others = [p for p in by_n[n] if p != pdb]
        cs, mm = [], []
        for p in others[:6]:
            tp = I.superpose_batch(np.asarray(I.load_univ(p)["nat_ca"], float)[None],
                                   Wt[b])[0].reshape(-1)
            mp = c - tp
            cs.append(float(mp @ mu / max(np.linalg.norm(mp) * nmu, 1e-12)))
            mm.append(tp)
        #: the POPULATION-MEAN direction: c minus the mean of the same-length deposited
        #: structures.  This is the floor's own carrier -- if the placebo's cos is just the
        #: shared `c`, this arm reproduces it while containing nothing about t' either.
        cos_pop = float("nan")
        if mm:
            tbar = np.mean(mm, 0)
            mb = c - tbar
            cos_pop = float(mb @ mu / max(np.linalg.norm(mb) * nmu, 1e-12))
        #: and the pure shrink-toward-the-centroid direction, which is what a correction along
        #: any of these actually does (lane Q: hull projection is entirely shrinkage)
        Cn = c.reshape(n, 3)
        sh = (Cn - Cn.mean(0)).reshape(-1)   # the SCALE direction: contract c about its own centre
        cos_shrink = float(sh @ mu / max(np.linalg.norm(sh) * nmu, 1e-12))
        rows.append(dict(pdb=pdb, n=n, fold=int(t["fold"]), n_other=len(cs),
                         mu_rmsd=float(nmu / np.sqrt(n)),
                         cos_placebo=float(np.mean(cs)) if cs else float("nan"),
                         cos_placebo_all=[float(v) for v in cs],
                         cos_popmean=cos_pop, cos_shrink=cos_shrink))

    R = [r for r in rows if r["n_other"] > 0]
    folds = np.array([r["fold"] for r in R])
    gl = sorted(set(folds.tolist()))
    g = lambda k: np.array([r[k] for r in R], float)                     # noqa: E731
    mur = g("mu_rmsd")

    def fold_se(x):
        m = np.array([np.nanmean(x[folds == f]) for f in gl])
        return float(m.std(ddof=1) / np.sqrt(len(m)))

    def payoff(cs):
        """One global step fitted leave-fold-out; a wrong sign costs."""
        al = np.empty(len(R))
        for f in gl:
            tr = folds != f
            j = int(np.argmin([np.sqrt(np.maximum(mur ** 2 - 2 * a * mur * cs + a * a, 0))[tr].mean()
                               for a in AGRID]))
            al[folds == f] = AGRID[j]
        e = np.sqrt(np.maximum(mur ** 2 - 2 * al * mur * cs + al ** 2, 0.0))
        return float(mur.mean() - e.mean()), float(al.mean()), e

    out = dict(prereg="s32/PREREG_S32_P.md @ 33dfe0d3 (H-P3 extension, EXPLORATORY)",
               basis="CA POINT CLOUD, score-top-128 medoid frame",
               ORACLE="ORACLE / NOT DEPLOYABLE -- every cos and payoff reads the native to SCORE "
                      "a direction that is itself native-free",
               n=len(R), n_excluded_no_same_length_partner=len(rows) - len(R), arms={})
    base = None
    for k, lab in (("cos_placebo", "c - t' for a MISMATCHED same-length deposited native"),
                   ("cos_popmean", "c - mean(same-length deposited natives)"),
                   ("cos_shrink", "c itself: pure shrink toward the frame origin")):
        cs = g(k)
        pay, al, e = payoff(cs)
        o = ST.compare(e, mur, folds=folds, names=[r["pdb"] for r in R],
                       label="corrected along %s - uncorrected (CLOUD, mu-RMSD)" % k,
                       seed_parts=("s32Ppl", str(SEED)))
        out["arms"][k] = dict(what=lab, cos_mean=float(np.nanmean(cs)), cos_se=fold_se(cs),
                              cos_median=float(np.nanmedian(cs)),
                              frac_positive=float(np.nanmean(cs > 0)),
                              payoff_A=pay, step=al,
                              vs_uncorrected={kk: o[kk] for kk in
                                              ("effect", "se", "mde", "effect_over_mde",
                                               "ci95_fold", "folds_same_sign", "n_better",
                                               "n_worse", "verdict")})
        if base is None:
            base = out["arms"][k]
    p, s = out["arms"]["cos_placebo"], out["arms"]["cos_shrink"]
    out["VERDICT"] = dict(
        shared_referent_floor=s["cos_mean"],
        placebo_excess_over_floor=float(p["cos_mean"] - s["cos_mean"]),
        reading="the placebo direction's cos is reported ONLY beside the floor; a cos that does "
                "not exceed the shrink-toward-origin arm is the shared `c` and nothing else")
    out["baseline_mu_rmsd"] = float(mur.mean())
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
