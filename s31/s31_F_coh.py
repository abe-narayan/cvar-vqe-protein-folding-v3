#!/usr/bin/env python
"""s31/s31_F_coh.py -- S31 LANE F: does a terminal operator LEAVE THE POOL'S AFFINE HULL?

Lane B's theorem (`s31/s31_B2_inpool.py:256-268`): for any readout `C = sum_m a_m x_m` with
`sum_m a_m = 1`, the pool's common mode passes through with coefficient EXACTLY one, so
`coh = corr(pair error, common mode)` is a function of the CONCENTRATION of `a` and not of the
ranker.  Only a readout that leaves the affine hull can move it.

**S30's 0.6931 admission bar is STRUCK from this file (2026-09-21 00:55) -- it grades a
CORRECTOR (the distogram's own prediction error) and not a READOUT (the emitted structure's
error).  See the note at the constants.  This file grades nothing; it reports LEVELS, a
PAIRED CONTRAST, and the affine-hull residual, which is the certificate that actually bears.**

`coh` is **ORACLE** -- both of its arguments need the native.  It is a diagnostic of WHY an
operator works and is never a deployable gate.  Every number here carries that label.

Arms, all on the shipped top-75 of the 500-member pool, CA point cloud:

    uniform_mean_pairspace   sanity: must be 1.0000 exactly (lane B's derivation)
    AVG                      the shipped uniform coordinate average          -- affine
    MED                      the consensus medoid MEMBER                     -- affine (a delta)
    AVG_RG                   AVG rescaled about its centroid to the members' mean Rg
    AVG_SEP                  AVG rescaled PER SEPARATION, re-embedded by MDS
                             (registered as 'NOT affine' -- MEASURED to be affine in
                              practice: hull residual 0.045 A RMS per coordinate)
    ORACLE_best              the best member of the 75                       -- affine (a delta)

    python s31/s31_F_coh.py
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

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s31.s31_F_terminal import sepprof, sepprof_batch, mds   # noqa: E402

RESULTS = os.path.join(HERE, "results")
OUT = os.path.join(RESULTS, "s31_F_coh.json")
ROWS = os.path.join(RESULTS, "s31_F_coh_rows.jsonl")
POOL_K, M, MIN_SEP, SEED = 500, 75, 2, 31006

#: THE 0.6931 ADMISSION BAR IS STRUCK (2026-09-21 00:55).  Lane V read S30 at source: S30's
#: `coh` grades a CORRECTOR -- its first argument is `(expected - d_nat) - correction`, the
#: DISTOGRAM'S OWN PREDICTION ERROR, an INPUT to scoring, and the 0.6931 row is the
#: `correction = 0` baseline of exactly that object.  THIS file's `coh` grades a READOUT --
#: the EMITTED STRUCTURE's error, the OUTPUT.  They share only `mu`.  The arithmetic settles
#: it without interpretation: the same pipeline, the same 126 targets and the same `mu` give
#: 0.6931 in the corrector table and 0.9780 in the readout table, 0.285 apart, because they
#: are two different errors.  Contract rule 8, sixth instance.  NO READOUT-SPACE BAR HAS
#: BEEN ESTABLISHED, so this file reports LEVELS and PAIRED CONTRASTS and grades nothing.


def one(t):
    pdb = t["pdb"]
    u = I.load_univ(pdb); rec = I.shipped_record(pdb); dg = I.distogram(pdb)
    n = int(u["n"])
    pool = np.asarray(u["order"], int)[:POOL_K]
    W = np.asarray(u["W"], float)[pool]
    rr = np.asarray(u["rr"], float)[pool]
    nat = np.asarray(u["nat_ca"], float)
    ii, jj = I.pair_index(n, MIN_SEP)
    D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
    order = np.argsort(I.shipped_score(dg, D), kind="stable")
    top = order[:M]
    if set(top.tolist()) != set(np.asarray(rec["sub"], int).tolist()):
        raise RuntimeError("REPRODUCTION GATE FAILED on %s" % pdb)

    pi, pj = I.pair_index(n, MIN_SEP)
    Dp = I.pair_dists(W[top], pi, pj)                     # (75, P)
    dn = I.pair_dists(nat[None], pi, pj)[0]               # ORACLE
    mu = (Dp - dn[None, :]).mean(0)                       # ORACLE common mode
    sdmu = float(mu.std())

    def coh_of(X):
        e = I.pair_dists(np.asarray(X, float)[None], pi, pj)[0] - dn
        if sdmu < 1e-12 or float(np.std(e)) < 1e-12:
            return float("nan")
        return float(np.corrcoef(e, mu)[0, 1])

    Pt = I.pairwise_rmsd(W[top])
    b = int(I.medoid(Pt))
    C = I.superpose_batch(W[top], W[top[b]]).mean(0)
    rg_members = np.sqrt(((W[top] - W[top].mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    lam = float(rg_members.mean()) / max(float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean())), 1e-9)
    Crg = (C - C.mean(0)) * lam + C.mean(0)
    p_mem = np.array(sepprof_batch(W[top], n), float)
    p_avg = np.array(sepprof(C, n), float)
    f = p_mem / np.maximum(p_avg, 1e-9)
    Dc = np.linalg.norm(C[:, None, :] - C[None, :, :], axis=2)
    sep = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    F = np.ones((n, n)); nz = sep > 0; F[nz] = f[sep[nz] - 1]
    Csep = mds(Dc * F)

    row = dict(pdb=pdb, n=n, fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18), sd_mu=sdmu)
    e_unif = Dp.mean(0) - dn
    row["uniform_mean_pairspace"] = (float(np.corrcoef(e_unif, mu)[0, 1])
                                     if sdmu > 1e-12 else float("nan"))
    row["AVG"] = coh_of(C)
    row["MED"] = coh_of(W[top[b]])
    row["AVG_RG"] = coh_of(Crg)
    row["AVG_SEP"] = coh_of(Csep)
    row["ORACLE_best"] = coh_of(W[top[int(np.argmin(rr[top]))]])
    row["cloud_AVG"] = I.ca_rmsd(C, nat)
    row["cloud_AVG_SEP"] = I.ca_rmsd(Csep, nat)
    # how far from the affine hull is AVG_SEP?  least-squares projection of Csep onto the
    # affine hull of the superposed 75, in COORDINATE space, after superposing on the medoid.
    Sup = I.superpose_batch(W[top], W[top[b]])
    Z = I.superpose_batch(np.asarray(Csep, float)[None], C)[0]
    A = Sup.reshape(M, -1)
    y = Z.reshape(-1)
    A1 = np.vstack([A, np.ones((1, A.shape[1]))])          # enforce sum a = 1 via a big penalty
    y1 = np.concatenate([y, [0.0]])
    Aa = np.vstack([A.T, np.ones((1, M)) * 1e3])
    yy = np.concatenate([y, [1e3]])
    a, *_ = np.linalg.lstsq(Aa, yy, rcond=None)
    resid = float(np.sqrt(((A.T @ a - y) ** 2).mean()))
    row["hull_residual_rms_per_coord"] = resid
    row["hull_residual_rms_per_coord_AVG"] = 0.0           # AVG is IN the hull by construction
    del A1, y1
    return row


def main():
    done = set()
    if os.path.exists(ROWS):
        for ln in open(ROWS):
            if ln.strip():
                done.add(json.loads(ln)["pdb"])
    ts = I.targets()
    for t in ts:
        if t["pdb"] in done:
            continue
        r = one(t)
        with open(ROWS, "a") as fh:
            fh.write(json.dumps(r) + "\n")
    R = [json.loads(ln) for ln in open(ROWS) if ln.strip()]
    R.sort(key=lambda r: r["pdb"])
    assert len(R) == 126, len(R)
    folds = np.array([r["fold"] for r in R]); names = [r["pdb"] for r in R]
    out = dict(seed=SEED, n=126,
               STRUCK="S30's 0.6931 admission bar is NOT applicable here -- it grades the "
                      "distogram's own prediction error (a CORRECTOR), not an emitted "
                      "structure's error (a READOUT). Same pipeline, same 126, same mu, "
                      "0.6931 vs 0.9780 -- two different errors. The claim is carried by "
                      "the PAIRED CONTRAST and by the affine-hull residual, not by a bar.",
               ORACLE="ORACLE / NOT DEPLOYABLE -- coh needs the native in BOTH arguments; it is a "
                      "diagnostic of why an operator works, never a deployable gate",
               basis="CA POINT CLOUD, pair-distance space, min_sep 2, shipped top-75",
               theorem="lane B, s31_B2_inpool.py:256-268 -- sum(a) = 1 passes the common mode "
                       "through with coefficient exactly one, so coh depends on the CONCENTRATION "
                       "of a and not on the ranker; only a readout OUTSIDE the affine hull moves it")
    arms = {}
    for k in ("uniform_mean_pairspace", "AVG", "MED", "AVG_RG", "AVG_SEP", "ORACLE_best"):
        v = np.array([r[k] for r in R], float)
        ok = np.isfinite(v)
        arms[k] = dict(mean=float(v[ok].mean()), median=float(np.median(v[ok])),
                       sd=float(v[ok].std(ddof=1)), n=int(ok.sum()),
                       NO_BAR="0.6931 is a CORRECTOR-space quantity and is struck here; "
                              "no readout-space bar exists, so this level is not graded")
    out["coh"] = arms
    a_sep = np.array([r["AVG_SEP"] for r in R], float)
    a_avg = np.array([r["AVG"] for r in R], float)
    o = ST.compare(a_sep, a_avg, folds=folds, names=names,
                   label="coh(AVG_SEP) - coh(AVG) (ORACLE diagnostic)", seed_parts=("s31Fcoh",))
    out["coh_AVG_SEP_minus_AVG"] = dict(
        effect=o["effect"], se=o["se"], mde=o["mde"], effect_over_mde=o["effect_over_mde"],
        ci95_fold=o["ci95_fold"], folds_same_sign=o["folds_same_sign"],
        W=o["n_better"], L=o["n_worse"])
    hr = np.array([r["hull_residual_rms_per_coord"] for r in R], float)
    out["affine_hull_departure"] = dict(
        AVG_SEP_residual_rms_A=dict(mean=float(hr.mean()), median=float(np.median(hr)),
                                    min=float(hr.min()), max=float(hr.max())),
        AVG_residual_rms_A=0.0,
        note="RMS per-coordinate residual of AVG_SEP against the best affine (sum a = 1) "
             "combination of the superposed 75. AVG is exactly in the hull, so its residual is 0 "
             "by construction and is printed as the reference, not as a measurement.")
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
