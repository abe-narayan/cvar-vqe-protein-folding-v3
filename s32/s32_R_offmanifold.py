#!/usr/bin/env python
"""S32 LANE R -- PRODUCER for `s32/results/s32_R_offmanifold_source.json`.

    python s32/s32_R_offmanifold.py

THE QUESTION (S32-L(R2), EXPLORATORY R-13).  The production average sits `d` off the
valid-chain manifold and pays a quadrature tax for it.  **Where does `d` come from?**  The
hypothesis is that it is the AVERAGING ARTEFACT -- the Jensen contraction you get from
meaning structures that disagree -- and therefore predictable from the pool's own spread with
no native anywhere.

    d        = Ca-RMSD(production chain, production cloud)        NATIVE-FREE
    spread75 = mean pairwise Ca-RMSD of the 75 members averaged   NATIVE-FREE

THE SHARED-REFERENT FLOOR IS MEASURED FIRST, not asserted.  Both quantities are RMSDs over
the same target and so share its scale; two quantities measured against a common reference
correlate by construction.  The null permutes `spread75` WITHIN chain-length strata, which
destroys the relation while preserving every shared referent, and the observed value must
clear that null's maximum -- not its 95th percentile.

WHY IT EXISTS AS A FILE.  This computation was first run as an inline shell heredoc and the
artefact had **no committed producer**.  Lane V's AUDIT 11 caught two others of mine on the
dilation arms; lane R's own producer audit caught this one and the sparse-control one.  Same
defect, three instances, one lane.

BASIS: `d`, `spread75` and `cos` are per-target diagnostics; `e` (cloud RMSD) and `cos` read
the native for EVALUATION only.  No endpoint claim is made here.
"""
from __future__ import annotations

import glob
import json
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST                                          # noqa: E402
from s12 import instrument as I                                          # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")
N_PERM = 4000
PERM_SEED = 7                  # pinned here, in the committed source


def partial(a, b, ctrl):
    """Spearman(a, b | ctrl), via residuals of the ranks on the control's rank."""
    from scipy.stats import rankdata
    ra, rb, rc = (rankdata(x).astype(float) for x in (a, b, ctrl))
    X = np.c_[np.ones(len(rc)), rc]
    res = lambda y: y - X @ np.linalg.lstsq(X, y, rcond=None)[0]         # noqa: E731
    return float(np.corrcoef(res(ra), res(rb))[0, 1])


def main():
    from scipy.stats import spearmanr
    lad = {}
    for f in sorted(glob.glob(os.path.join(RESULTS, "s32_R_laddernull_shard*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    lad[r["pdb"]] = r

    rows = []
    for p in [t["pdb"] for t in I.targets() if t["pdb"] in lad]:
        u = I.load_univ(p)
        pool = I.pool_idx(u)
        sub = np.asarray(I.shipped_record(p)["sub"], int)
        W = np.asarray(u["W"], float)[pool[sub]].astype(np.float32).astype(float)
        P = I.pairwise_rmsd(W)
        spread = float(P[np.triu_indices(len(W), 1)].mean())
        g = lad[p]["rungs"]["prod"]
        e, d, ch = g["cloud_rmsd"], g["d_to_cloud"], g["chain_rmsd"]
        rows.append(dict(pdb=p, fold=lad[p]["fold"], n=lad[p]["n"], spread75=spread,
                         d=d, e=e, chain=ch,
                         cos=(e ** 2 + d ** 2 - ch ** 2) / max(2 * e * d, 1e-12)))

    D = np.array([r["d"] for r in rows])
    S = np.array([r["spread75"] for r in rows])
    E = np.array([r["e"] for r in rows])
    N = np.array([r["n"] for r in rows], float)
    CO = np.array([r["cos"] for r in rows])
    F = np.array([r["fold"] for r in rows], int)

    rho = float(spearmanr(D, S).statistic)
    rng = np.random.default_rng(PERM_SEED)
    nulls = np.empty(N_PERM)
    for t in range(N_PERM):
        Sp = S.copy()
        for nn in np.unique(N):
            m = N == nn
            if m.sum() > 1:
                Sp[m] = rng.permutation(Sp[m])
        nulls[t] = spearmanr(D, Sp).statistic

    out = dict(
        n=len(D), INCOMPLETE=len(D) < 126, perm_seed=PERM_SEED,
        spearman_d_spread=rho, pearson=float(np.corrcoef(D, S)[0, 1]),
        r2_linear=float(np.corrcoef(D, S)[0, 1] ** 2),
        ols_slope=float(np.polyfit(S, D, 1)[0]),
        ols_intercept=float(np.polyfit(S, D, 1)[1]),
        partial_given_n=partial(D, S, N), partial_given_e=partial(D, S, E),
        within_n_permutation_null=dict(
            mean=float(nulls.mean()), p95=float(np.percentile(nulls, 95)),
            p999=float(np.percentile(nulls, 99.9)), max=float(nulls.max()),
            n_draws=N_PERM,
            what="spread75 permuted WITHIN chain-length strata: destroys the relation, "
                 "preserves every shared referent"),
        per_fold_spearman={int(q): float(spearmanr(D[F == q], S[F == q]).statistic)
                           for q in sorted(set(F.tolist()))},
        ratio_d_over_spread=dict(mean=float((D / S).mean()), sd=float((D / S).std()),
                                 cv=float((D / S).std() / (D / S).mean()),
                                 what="cv this large means the relation is MONOTONE, not a "
                                      "proportionality -- quote rho, never a coefficient"),
        cos_prod=dict(mean=float(CO.mean()),
                      se=float(CO.std(ddof=1) / np.sqrt(len(CO))),
                      median=float(np.median(CO)), n_positive=int((CO > 0).sum())),
        spearman_d_e=float(spearmanr(D, E).statistic),
        spearman_spread_e=float(spearmanr(S, E).statistic),
        d_mean=float(D.mean()), spread_mean=float(S.mean()),
        rows=rows,
        note=("EXPLORATORY R-13. d and spread75 are NATIVE-FREE on both sides; e, chain and "
              "cos read the native for EVALUATION only."))
    ST.save_atomic(os.path.join(RESULTS, "s32_R_offmanifold_source.json"), out,
                   module_file=__file__)

    print("n = %d%s" % (out["n"], "  (INCOMPLETE)" if out["INCOMPLETE"] else ""))
    print("d(prod) mean %.4f   top-75 member spread mean %.4f" % (D.mean(), S.mean()))
    print("spearman(d, spread)        %+.4f   pearson %+.4f   R2 %.4f"
          % (rho, out["pearson"], out["r2_linear"]))
    print("  partial | chain length    %+.4f" % out["partial_given_n"])
    print("  partial | native error    %+.4f" % out["partial_given_e"])
    print("  per fold                  %s"
          % {k: round(v, 3) for k, v in out["per_fold_spearman"].items()})
    nl = out["within_n_permutation_null"]
    print("  within-n permutation null: mean %+.4f  p95 %+.4f  p99.9 %+.4f  MAX %+.4f (%d draws)"
          % (nl["mean"], nl["p95"], nl["p999"], nl["max"], nl["n_draws"]))
    print("  ratio d/spread mean %.4f sd %.4f cv %.3f -> MONOTONE, not a proportionality"
          % (out["ratio_d_over_spread"]["mean"], out["ratio_d_over_spread"]["sd"],
             out["ratio_d_over_spread"]["cv"]))
    print("\ncos(prod) mean %+.4f (SE %.4f) median %+.4f  %d/%d positive"
          % (out["cos_prod"]["mean"], out["cos_prod"]["se"], out["cos_prod"]["median"],
             out["cos_prod"]["n_positive"], out["n"]))
    print("wrote", os.path.join(RESULTS, "s32_R_offmanifold_source.json"))


if __name__ == "__main__":
    main()
