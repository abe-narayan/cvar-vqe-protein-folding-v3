#!/usr/bin/env python
"""s29/s29_M_F2_supply.py -- F2 STAGE 0: the projection-free supply audit.

Pre-registered in `s29/PREREG_S29_M_F2.md` (+ addenda 1, 2, 3). Read that first: the stopping
rule, the shrink twin's justification, the bits currency and the three over-reading caveats all
live there, not here.

THE QUESTION. The shell profile's ORACLE version is worth 0.75 A on the point cloud
(`s12/obj_FINDINGS.md:188`, 2.299 vs 3.048) -- a displacement cosine of about 0.66, ~4.7x the
0.140 that lane T's bound (S29-L23) assumes no native-free field can exceed (assumption B2).
Five native-free supplies of that profile already exist and all lose (`:306-331`). So: how much
of that 0.66 survives leave-fold-out, in the bound's own currency?

WHY THIS COSTS ALMOST NOTHING. The cosine that decides B2 needs NO projection: emitting a point
cloud is a score, a tie-safe top-75 and a 75-member coordinate average. The built chain is spent
only on an arm that has already cleared (prereg addendum 1, A1.3).

THE OBJECTIVE, S12's own swap convention (`s12/obj_profile.py:19`) so the comparison is
like-for-like:  target_p = q'_{sep(p)} + ( E[d]_p - q_disto_{sep(p)} )   -- only the profile moves.

ARMS
    PROD          q' = q_disto            the shipped distogram profile (the incumbent)
    POOL          q' = q_pool             pure typicality, r = 1, the zero-information arm
    RATIO         q' = q_pool * r_hat     LFO ridge on native-free features            [primary]
    RSHRINK       q' = q_pool * (1+s*z)   RATIO's realised shrink, random direction    [the twin]
    ORACLE_RATIO  q' = q_pool * r_true    the ceiling of this parameterisation   ORACLE DIAGNOSTIC
    ORACLE_PROF   q' = q_true             s12's 2.402 / 2.299 row                ORACLE DIAGNOSTIC

REPORTED PER ARM: point-cloud RMSD (ORACLE), the displacement cosine against PROD (ORACLE), the
native's percentile in its own pool (the meter axis that SURVIVED -- contract addendum 5; rule 20's
cosine justification is withdrawn and is not cited here), the realised shrink s, the emitted bond
and Rg, and bits_delivered = 7 - log2(rank of the ORACLE-best member of the FIXED production
top-128 under the arm's own score).

    python s26/jobrun.py --agent S29M --tag CPU --name m_f2_supply --est-ram 1.0 -- \
        python s29/s29_M_F2_supply.py --full
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import run_pool as RP             # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
CELLS = os.path.join(RESULTS, "s29_M_F2_cells.npz")
ROWS = os.path.join(RESULTS, "s29_M_F2_supply_rows.jsonl")

PROBE = ["1A13", "1I6Y", "1M02", "2BFI", "2LWS", "2MP9",
         "2P5H", "5Z5W", "6MBM", "7JGX", "8HVS", "9KAR"]
ARMS = ("PROD", "POOL", "RATIO", "RSHRINK", "ORACLE_RATIO", "ORACLE_PROF")
NBITS = 7                      # log2(128): the production top-128 register
TOP128 = 128


# ============================================================ per-target primitives
def target_data(pdb):
    """Everything one target needs, all cached. `nat`/`q_true` are ORACLE and used only for
    the diagnostics that say so."""
    cand, ch, _ = RP.channels_for(pdb)
    dg = I.distogram(pdb, cand.seq, cand.fold)
    i, j = I.pair_index(cand.n)
    D = I.pair_dists(cand.W, i, j).astype(np.float32).astype(float)     # (500, npairs)
    sep = (j - i).astype(int)
    us = np.unique(sep)
    exp = np.asarray(dg["expected"], float)                              # (npairs,)
    sd = np.asarray(dg["sd"], float)
    w = 1.0 / np.maximum(sd + 0.5, 1e-6)
    w = w / max(w.mean(), 1e-12)
    nat = np.asarray(cand.nat_ca, float)                                 # ORACLE
    dnat = np.linalg.norm(nat[i] - nat[j], axis=1)                       # ORACLE
    q_true = np.array([dnat[sep == s].mean() for s in us])               # ORACLE
    q_disto = np.array([exp[sep == s].mean() for s in us])
    q_pool = np.array([D[:, sep == s].mean() for s in us])
    return dict(pdb=pdb, cand=cand, n=int(cand.n), fold=int(cand.fold), seq=cand.seq,
                i=i, j=j, D=D, sep=sep, us=us, exp=exp, w=w, nat=nat,
                q_true=q_true, q_disto=q_disto, q_pool=q_pool)


def features_of(d):
    """(n_shells, F) NATIVE-FREE features for the profile ratio, one row per shell.

    Named blocks, per prereg section 3: the distogram's own ratio (the thing to beat), the
    compactness block (lane L S29-L19: the in-band axis is compactness), the pool's own shape,
    length, shell index, and composition/properties.
    """
    from core.predict import _PROP_TABLE
    us, D, sep = d["us"], d["D"], d["sep"]
    r_disto = d["q_disto"] / np.maximum(d["q_pool"], 1e-9)
    # compactness block, all native-free
    W = d["cand"].W
    rg_mem = np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(2).mean(1))     # (500,)
    rg_pool_mean, rg_pool_sd = float(rg_mem.mean()), float(rg_mem.std())
    rg_pool_skew = float(((rg_mem - rg_mem.mean()) ** 3).mean() / max(rg_pool_sd ** 3, 1e-9))
    n = d["n"]
    # the distogram's implied Rg: Rg^2 = (1/(2 n^2)) sum_ij d_ij^2, with |i-j|=1 at 3.81 A
    s2 = 2.0 * (d["exp"] ** 2).sum() + 2.0 * (n - 1) * 3.81 ** 2
    rg_pred = float(np.sqrt(max(s2 / (2.0 * n * n), 1e-9)))
    rg_ratio = rg_pred / max(rg_pool_mean, 1e-9)         # prediction/pool disagreement (memory)
    A = I.ALPHABET
    comp = np.array([d["seq"].count(a) for a in A], float) / max(n, 1)
    idx = np.array([A.index(a) if a in A else 0 for a in d["seq"]])
    props = _PROP_TABLE[idx]
    base = np.concatenate([[rg_pool_mean, rg_pool_sd, rg_pool_skew, rg_pred, rg_ratio,
                            float(n)], comp, props.mean(0), props.std(0)])
    rows = []
    for k, s in enumerate(us):
        pool_spread = float(D[:, sep == s].std())
        rows.append(np.concatenate([[r_disto[k], float(s), float(s) / n, pool_spread,
                                     float(d["q_pool"][k])], base]))
    return np.asarray(rows, float)


def score_of(d, qprime):
    """The shipped weighted-L1 score against the swapped-profile target (S12's convention)."""
    tgt = qprime[np.searchsorted(d["us"], d["sep"])] + (d["exp"] - d["q_disto"][np.searchsorted(d["us"], d["sep"])])
    return (np.abs(d["D"] - tgt[None, :]) * d["w"][None, :]).mean(1)


def native_score(d, qprime):
    """The same score evaluated on the NATIVE structure -- for the percentile meter (ORACLE)."""
    dnat = np.linalg.norm(d["nat"][d["i"]] - d["nat"][d["j"]], axis=1)
    k = np.searchsorted(d["us"], d["sep"])
    tgt = qprime[k] + (d["exp"] - d["q_disto"][k])
    return float((np.abs(dnat - tgt) * d["w"]).mean())


def emit(d, sc, key):
    """Score -> tie-safe top-75 -> uniform coordinate average. NO projection."""
    top = np.lexsort((key, RP.zr(sc)))[:RP.M]
    C, _ = H.readout_uniform(d["cand"], top)
    return C, top


def bits_delivered(d, sc, prod_top128):
    """7 - log2(rank of the ORACLE-best member of the FIXED production top-128 under `sc`).

    ORACLE DIAGNOSTIC. The candidate SET is held fixed at production's own prefix, so only the
    RANKING varies. Caveat (prereg addendum 2): bits about the BEST member are not bits the
    mean-consuming readout can spend.
    """
    rr = np.asarray(d["cand"].oracle_rr, float)[prod_top128]      # ORACLE
    best_local = int(np.argmin(rr))
    order = np.argsort(sc[prod_top128], kind="stable")
    rank = int(np.flatnonzero(order == best_local)[0]) + 1
    return float(NBITS - math.log2(rank)), rank


# ============================================================ the LFO ridge
def ridge_fit(X, y, lam):
    Xb = np.hstack([X, np.ones((len(X), 1))])
    A = Xb.T @ Xb + lam * np.eye(Xb.shape[1])
    A[-1, -1] -= lam                                   # do not penalise the intercept
    return np.linalg.solve(A, Xb.T @ y)


def ridge_pred(X, beta):
    return np.hstack([X, np.ones((len(X), 1))]) @ beta


def fit_lfo(cells, folds, lams=(1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3)):
    """Leave-fold-out r_hat, with the ridge penalty chosen by nested CV INSIDE the training
    folds only. Standardisation is fitted on the training folds too."""
    X = np.vstack([c["X"] for c in cells])
    y = np.concatenate([c["r_true"] for c in cells])
    g = np.concatenate([np.full(len(c["r_true"]), c["fold"]) for c in cells])
    out = np.zeros_like(y)
    chosen = {}
    for f in sorted(set(folds)):
        tr, te = g != f, g == f
        # nested CV over the training folds
        inner = sorted(set(g[tr].tolist()))
        best, best_lam = np.inf, lams[0]
        for lam in lams:
            err = []
            for h in inner:
                itr, ite = tr & (g != h), tr & (g == h)
                mu, sd = X[itr].mean(0), X[itr].std(0) + 1e-9
                b = ridge_fit((X[itr] - mu) / sd, y[itr], lam)
                err.append(np.mean((ridge_pred((X[ite] - mu) / sd, b) - y[ite]) ** 2))
            m = float(np.mean(err))
            if m < best:
                best, best_lam = m, lam
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        b = ridge_fit((X[tr] - mu) / sd, y[tr], best_lam)
        out[te] = ridge_pred((X[te] - mu) / sd, b)
        chosen[int(f)] = float(best_lam)
    # split back per target
    k = 0
    for c in cells:
        m = len(c["r_true"])
        c["r_hat"] = out[k:k + m]
        k += m
    return chosen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    pdbs = P.targets() if a.full else PROBE
    t0 = time.time()

    # ---- pass 1: build every target's cells (cheap; no readout yet)
    cells, data = [], {}
    for k, pdb in enumerate(pdbs):
        d = target_data(pdb)
        data[pdb] = d
        cells.append(dict(pdb=pdb, fold=d["fold"], X=features_of(d),
                          r_true=d["q_true"] / np.maximum(d["q_pool"], 1e-9)))
        if (k + 1) % 20 == 0:
            print("  built %d/%d (%.1f min)" % (k + 1, len(pdbs), (time.time() - t0) / 60), flush=True)
    folds = sorted({c["fold"] for c in cells})
    chosen = fit_lfo(cells, folds)
    print("ridge lambda per fold:", chosen, flush=True)

    # ---- the supply gap, before any readout
    r_true = np.concatenate([c["r_true"] for c in cells])
    r_hat = np.concatenate([c["r_hat"] for c in cells])
    r_disto = np.concatenate([c["X"][:, 0] for c in cells])
    supply = dict(
        n_cells=int(len(r_true)),
        corr_rhat_rtrue=float(np.corrcoef(r_hat, r_true)[0, 1]),
        corr_rdisto_rtrue=float(np.corrcoef(r_disto, r_true)[0, 1]),
        sd_rtrue=float(np.std(r_true - 1.0)), sd_rhat=float(np.std(r_hat - 1.0)),
        sd_rdisto=float(np.std(r_disto - 1.0)),
        shrink_rhat=float(np.std(r_hat - 1.0) / max(np.std(r_true - 1.0), 1e-9)),
        shrink_rdisto=float(np.std(r_disto - 1.0) / max(np.std(r_true - 1.0), 1e-9)),
        var_explained_rhat=float(1.0 - np.var(r_hat - r_true) / max(np.var(r_true), 1e-12)),
        var_explained_rdisto=float(1.0 - np.var(r_disto - r_true) / max(np.var(r_true), 1e-12)),
        ridge_lambda_per_fold=chosen)
    print("SUPPLY GAP:", json.dumps({k: v for k, v in supply.items() if k != "ridge_lambda_per_fold"},
                                    indent=1), flush=True)

    # ---- pass 2: emit each arm's cloud (NO projection) and measure
    done = set()
    if os.path.exists(ROWS):
        for line in open(ROWS, encoding="utf-8"):
            try:
                r = json.loads(line)
                done.add((r["pdb"], r["arm"]))
            except Exception:
                pass
    for c in cells:
        pdb = c["pdb"]
        d = data[pdb]
        if all((pdb, arm) in done for arm in ARMS):
            continue
        key = RP.rng_for(pdb, "tiekey").random(d["cand"].k)
        s_real = float(np.std(c["r_hat"] - 1.0) / max(np.std(c["r_true"] - 1.0), 1e-9))
        z = RP.rng_for(pdb, "F2shrink").normal(0.0, 1.0, len(c["r_true"]))
        z = z / max(z.std(), 1e-9) * np.std(c["r_true"] - 1.0)
        qp = dict(PROD=d["q_disto"], POOL=d["q_pool"],
                  RATIO=d["q_pool"] * c["r_hat"],
                  RSHRINK=d["q_pool"] * (1.0 + s_real * z),
                  ORACLE_RATIO=d["q_pool"] * c["r_true"],
                  ORACLE_PROF=d["q_true"])
        sc_prod = score_of(d, qp["PROD"])
        prod_top128 = np.lexsort((key, RP.zr(sc_prod)))[:TOP128]
        C_prod, _ = emit(d, sc_prod, key)
        e_prod = (I.superpose_batch(C_prod[None], d["nat"])[0] - d["nat"]).ravel()   # ORACLE
        rows = []
        for arm in ARMS:
            sc = score_of(d, qp[arm])
            C, top = emit(d, sc, key)
            bits, rank = bits_delivered(d, sc, prod_top128)
            # the displacement cosine (ORACLE, post hoc), the bound's own currency
            Cs = I.superpose_batch(C[None], d["nat"])[0]
            Ps = I.superpose_batch(C_prod[None], d["nat"])[0]
            u = (Cs - Ps).ravel()
            cos = float(np.dot(u, -e_prod) / max(np.linalg.norm(u) * np.linalg.norm(e_prod), 1e-12))
            pct = float((sc < native_score(d, qp[arm])).mean())      # native's percentile, ORACLE
            bond = float(np.linalg.norm(np.diff(C, axis=0), axis=1).mean())
            rg = float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean()))
            rows.append(dict(pdb=pdb, n=d["n"], fold=d["fold"], arm=arm,
                             rmsd_cloud=float(I.ca_rmsd(C, d["nat"])),
                             cos_vs_prod=cos, bits=bits, best_rank=rank,
                             native_pct=pct, shrink=s_real if arm == "RATIO" else
                             (s_real if arm == "RSHRINK" else float("nan")),
                             bond=bond, rg=rg,
                             overlap_prod=float(len(set(top.tolist())
                                                    & set(np.lexsort((key, RP.zr(sc_prod)))[:RP.M].tolist())) / RP.M)))
        with open(ROWS, "a", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
    np.savez_compressed(CELLS, **{c["pdb"]: np.vstack([c["r_true"], c["r_hat"]]) for c in cells})
    ST.save_atomic(os.path.join(RESULTS, "s29_M_F2_supply.json"),
                   dict(supply=supply, n_targets=len(pdbs), probe=not a.full,
                        secs=time.time() - t0), module_file=__file__)
    print("done: %s (%.1f min)" % (ROWS, (time.time() - t0) / 60))


if __name__ == "__main__":
    sys.exit(main())
