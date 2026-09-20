#!/usr/bin/env python
"""s29/s29_M_F1.py -- F1: swap ONLY the selection functional in the shipped pipeline.

Pre-registered in `s29/PREREG_S29_M_F1.md` (read it first; the falsifier, the rule-10 declaration
and the reason the coordinator's monotone control is replaced are all there).

Arms, all through the IDENTICAL pool / top-75 / uniform average / production projection:
    PROD    the shipped L1 Bayes risk                     (reproduces chain_rows.jsonl :: DIS)
    LOG     sum_p -log(prob[p, bin(d_p)] + 1e-4)          lane D's metered cost_nll  [PRIMARY]
    LOGW    the same, weighted by the shipped 1/(sd+0.5)  (LOG - LOGW prices the weight)
    L2RISK  mean_p w_p sum_c prob[p,c] (d_p - C_c)^2      matched functional control
    LOGPERM LOG scoring pair p's distance against pair perm(p)'s posterior  information control
    SWAPCTL PROD's top-75 with k_swap members exchanged at random for candidates from ranks
            76..500, where k_swap is the count PROD keeps and LOG drops -- the ZERO-INFORMATION
            re-ordering at LOG's own exchange rate, 4 draws (PREREG addendum 1)

Per target it also records the MECHANISM (contract rule 18): the emitted cloud's mean virtual
CA-CA bond and radius of gyration against the native's, the projection price, the top-75 set
overlap with PROD, the tie fraction at the cut, the posterior-floor hit rate, and the ORACLE
error-direction cosine against PROD (S24 L2/L3 parallel-bias check).

    python s26/jobrun.py --agent S29M --tag CPU --name m_f1_probe --est-ram 1.0 -- \
        python s29/s29_M_F1.py --probe
    python s26/jobrun.py --agent S29M --tag CPU --name m_f1_full --est-ram 1.0 -- \
        python s29/s29_M_F1.py --full

Rows (resumable, one line per (pdb, arm)): s29/results/s29_M_F1_rows.jsonl
"""
from __future__ import annotations

import argparse
import json
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
from core import predict as PRD            # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
ROWS = os.path.join(RESULTS, "s29_M_F1_rows.jsonl")

#: inherited from lane X's metered cost, NOT chosen here (`s29/s29_X_config.py:107`)
EPS_P = 1e-4
#: the pinned probe set (`s29/s29_X_config.py:120`, = P.targets()[::11][:12])
PROBE = ["1A13", "1I6Y", "1M02", "2BFI", "2LWS", "2MP9",
         "2P5H", "5Z5W", "6MBM", "7JGX", "8HVS", "9KAR"]
ARMS = ("PROD", "LOG", "LOGW", "L2RISK", "LOGPERM", "SWAPCTL")
#: SWAPCTL draws (PREREG addendum 1: the zero-information re-ordering at LOG's own
#: exchange rate, matched in the operator's space)
N_SWAP_DRAWS = 4
NATIVE_BOND = 3.8122          # ARCHITECTURE.md section 0, the natives' mean virtual CA-CA bond


# ============================================================ the five scores
def scores_for(dg, D, w):
    """Every arm's score vector for one target. D (k, npairs) A, w (npairs,) the shipped weight.

    PROD is NOT recomputed here -- it comes from the cached DIS channel, so the comparator is
    bit-identical to the production anchor (audit check 5). The other four are built from the
    SAME posterior `dg`.
    """
    prob = np.clip(np.asarray(dg["prob"], float), 0.0, 1.0)      # (npairs, 17)
    cen = np.asarray(dg["centres"], float)                        # (17,)
    lp = -np.log(prob + EPS_P)                                    # (npairs, 17) nats
    b = np.digitize(D, PRD.BIN_EDGES)                             # (k, npairs) bin of each distance
    rows = np.arange(lp.shape[0])[None, :]
    log_pp = lp[rows, b]                                          # (k, npairs)
    out = {}
    out["LOG"] = log_pp.sum(1)
    out["LOGW"] = (log_pp * w[None, :]).mean(1)
    risk2 = (prob[None, :, :] * (D[:, :, None] - cen[None, None, :]) ** 2).sum(2)  # (k, npairs)
    out["L2RISK"] = (risk2 * w[None, :]).mean(1)
    # the floor-hit rate: pairs whose occupied bin carries essentially zero posterior mass
    floor = float((prob[rows, b] <= EPS_P).mean())
    return out, log_pp, floor, lp, b, rows


def geometry_of(C, nat):
    """The mechanism quantities of an emitted cloud."""
    C = np.asarray(C, float)
    bond = float(np.linalg.norm(np.diff(C, axis=0), axis=1).mean())
    rg = float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean()))
    rg_nat = float(np.sqrt(((nat - nat.mean(0)) ** 2).sum(1).mean()))
    return bond, rg, rg_nat


def err_vector(X, nat):
    """The emitted structure's error vector in the native's frame (rigid body removed)."""
    Xs = I.superpose_batch(np.asarray(X, float)[None], np.asarray(nat, float))[0]
    return (Xs - np.asarray(nat, float)).ravel()


def swap_control(cand, prod_order, k_swap, pdb, prod_err, floor):
    """The zero-information re-ordering at LOG's own exchange rate (PREREG addendum 1).

    Exchange `k_swap` members of PROD's top-75, chosen uniformly at random, for `k_swap`
    candidates drawn uniformly from PROD's ranks 76..500. `N_SWAP_DRAWS` independent draws
    through the IDENTICAL average and projection; the arm's value is the mean over draws.
    """
    inside = np.asarray(prod_order[:RP.M], int)
    outside = np.asarray(prod_order[RP.M:], int)
    ch_v, cl_v, bd_v, ov_v, cos_v = [], [], [], [], []
    t0 = time.time()
    for d in range(N_SWAP_DRAWS):
        rng = RP.rng_for(pdb, "F1swap%d" % d)
        if k_swap > 0:
            drop = rng.choice(len(inside), size=k_swap, replace=False)
            add = rng.choice(len(outside), size=k_swap, replace=False)
            top = np.concatenate([np.delete(inside, drop), outside[add]])
        else:
            top = inside.copy()
        C, _ = H.readout_uniform(cand, top)
        ca = H.readout_projected(cand, C)
        bond, rg, rg_nat = geometry_of(C, cand.nat_ca)
        cl_v.append(float(I.ca_rmsd(C, cand.nat_ca)))
        ch_v.append(float(I.ca_rmsd(ca, cand.nat_ca)))
        bd_v.append(bond)
        ov_v.append(float(len(set(top.tolist()) & set(inside.tolist())) / RP.M))
        e = err_vector(ca, cand.nat_ca)
        cos_v.append(float(np.dot(e, prod_err) / max(np.linalg.norm(e) * np.linalg.norm(prod_err), 1e-12))
                     if prod_err is not None else float("nan"))
    return dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), arm="SWAPCTL",
                rmsd_cloud=float(np.mean(cl_v)), rmsd_chain=float(np.mean(ch_v)),
                rmsd_chain_sd=float(np.std(ch_v, ddof=1)) if len(ch_v) > 1 else 0.0,
                rmsd_chain_draws=[float(v) for v in ch_v],
                bond_cloud=float(np.mean(bd_v)), bond_chain=float("nan"),
                rg_cloud=float("nan"), rg_native=float(np.sqrt(((cand.nat_ca - cand.nat_ca.mean(0)) ** 2).sum(1).mean())),
                overlap_prod=float(np.mean(ov_v)), tie_frac_at_cut=float("nan"),
                floor_rate=floor, medoid=-1, k_swap=int(k_swap), n_draws=int(N_SWAP_DRAWS),
                cos_err_vs_prod=float(np.nanmean(cos_v)), err_norm=float("nan"),
                secs=float(time.time() - t0))


def run_target(pdb, done):
    todo = [a for a in ARMS if (pdb, a) not in done]
    if not todo:
        return []
    cand, ch, _ = RP.channels_for(pdb)
    dg = I.distogram(pdb, cand.seq, cand.fold)
    i, j = I.pair_index(cand.n)
    D = I.pair_dists(cand.W, i, j).astype(np.float32).astype(float)
    # the shipped per-pair weight, exactly as core/predict.py:416-417 builds it
    sd = np.asarray(dg["sd"], float)
    w = 1.0 / np.maximum(sd + 0.5, 1e-6)
    w = w / max(w.mean(), 1e-12)
    sc, log_pp, floor, lp, b, rows_ix = scores_for(dg, D, w)
    sc["PROD"] = np.asarray(ch["DIS"], float)
    # DEFECT CAUGHT ON THE 12-TARGET PROBE, 2026-09-20, and fixed here rather than quietly:
    # the first version permuted the per-pair values AFTER the gather (`log_pp[:, perm].sum(1)`),
    # which for a SUM aggregator is the IDENTITY -- the control reproduced LOG to every digit.
    # The information control must break the correspondence BETWEEN a pair's distance and a
    # pair's posterior, i.e. permute the POSTERIOR's rows and gather with the TRUE bins.
    perm = RP.rng_for(pdb, "F1perm").permutation(lp.shape[0])
    sc["LOGPERM"] = lp[perm][rows_ix, b].sum(1)

    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    prod_order = np.lexsort((key, RP.zr(sc["PROD"])))
    prod_top = set(prod_order[:RP.M].tolist())
    log_top = set(np.lexsort((key, RP.zr(sc["LOG"])))[:RP.M].tolist())
    k_swap = int(len(prod_top - log_top))
    prod_err = None
    rows = []
    for arm in ARMS:
        if arm == "SWAPCTL":
            if arm in todo:
                rows.append(swap_control(cand, prod_order, k_swap, pdb, prod_err, floor))
            continue
        if arm != "PROD" and arm not in todo:
            continue          # already on disk; PROD is always recomputed for prod_err
        t0 = time.time()
        E = RP.zr(sc[arm])
        order = np.lexsort((key, E))
        top = order[:RP.M]
        C, medoid = H.readout_uniform(cand, top)
        ca = H.readout_projected(cand, C)
        bond, rg, rg_nat = geometry_of(C, cand.nat_ca)
        cb, _, _ = geometry_of(ca, cand.nat_ca)
        e = err_vector(ca, cand.nat_ca)
        if arm == "PROD":
            prod_err = e
        cos = float(np.dot(e, prod_err) / max(np.linalg.norm(e) * np.linalg.norm(prod_err), 1e-12)) \
            if prod_err is not None else float("nan")
        cut = float(E[top].max())
        if arm not in todo:
            continue
        rows.append(dict(
            pdb=pdb, n=int(cand.n), fold=int(cand.fold), arm=arm,
            rmsd_cloud=float(I.ca_rmsd(C, cand.nat_ca)),
            rmsd_chain=float(I.ca_rmsd(ca, cand.nat_ca)),
            bond_cloud=bond, bond_chain=cb, rg_cloud=rg, rg_native=rg_nat,
            overlap_prod=float(len(set(top.tolist()) & prod_top) / RP.M),
            tie_frac_at_cut=float((np.abs(E - cut) <= 1e-12).sum() / RP.M),
            floor_rate=floor, medoid=int(medoid),
            cos_err_vs_prod=cos, err_norm=float(np.linalg.norm(e)),
            secs=float(time.time() - t0)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    pdbs = PROBE if a.probe else (P.targets() if a.full else PROBE)
    done = set()
    if os.path.exists(ROWS):
        with open(ROWS, encoding="utf-8") as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                    done.add((r["pdb"], r["arm"]))
                except Exception:
                    pass
    t0 = time.time()
    for k, pdb in enumerate(pdbs):
        rows = run_target(pdb, done)
        if rows:
            with open(ROWS, "a", encoding="utf-8") as fh:
                for r in rows:
                    fh.write(json.dumps(r) + "\n")
            for r in rows:
                done.add((r["pdb"], r["arm"]))
        print("[%3d/%d] %s %d rows (%.1f min elapsed)" % (k + 1, len(pdbs), pdb, len(rows),
                                                          (time.time() - t0) / 60.0), flush=True)
    print("done:", ROWS)


if __name__ == "__main__":
    sys.exit(main())
