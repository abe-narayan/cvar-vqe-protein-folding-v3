#!/usr/bin/env python
"""s32/s32_Q0_invariance.py -- LANE Q, question Q0.

FALSIFY S31 section 5.1's claim that the deployed CVaR-VQE Hamiltonian
`E = _zrank(pool["sc"][o])` is a target-independent constant to within 0.0407 max-norm.

Registered in `s32/PREREG_S32_Q.md` at commit a8f9d6a7.  Three hypotheses:

  Q0-H1  the residual is EXACTLY the tie pattern and nothing else
  Q0-H2  no other target-dependent quantity enters `min_p F(p)`
  Q0-H3  the residual is worth < 0.7x MDE at the endpoint

The self-test that can fail (contract rule 5).  S31's own verification was run on random
floats, WHICH NEVER TIE, so it was structurally incapable of testing its own caveat.  This
one runs on the real 126 pools AND on two synthetic vectors built to break the claim:
a deliberately tied vector (must NOT give the ramp) and a deliberately unsorted vector
(must NOT give the ramp).  If those two do not break it, this test is decoration.

    python s32/s32_Q0_invariance.py [--limit N] [--draws 8]
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

from scipy.stats import rankdata                                        # noqa: E402

import core                                                            # noqa: E402
from core.pipeline import VQE_LFO, _zrank, consensus_medoid            # noqa: E402
from s8 import consensus2 as cc                                        # noqa: E402
from s31 import s31_P_substitute as P                                  # noqa: E402

CACHE = os.path.join(ROOT, "s31", "results", "s31_C_cache")
OUT = os.path.join(HERE, "results")
NQ = 128                                                               # 2**cfg.vqe_qubits


# ------------------------------------------------------------------ the self-test
def self_test():
    """Two inputs CONSTRUCTED so the claim must break on them, plus one that must hold."""
    ramp = _zrank(np.arange(NQ, dtype=float))

    tied = np.arange(NQ, dtype=float)
    tied[10:20] = tied[10]                                             # one block of 10 ties
    e_tied = _zrank(tied)
    broke_tied = float(np.abs(e_tied - ramp).max())

    rng = np.random.default_rng(0)
    uns = rng.permutation(NQ).astype(float)                            # deliberately UNSORTED
    e_uns = _zrank(uns)
    broke_uns = float(np.abs(e_uns - ramp).max())

    clean = np.sort(rng.normal(size=NQ))                               # sorted, no ties
    ok_clean = float(np.abs(_zrank(clean) - ramp).max())

    out = {"tied_maxdev": broke_tied, "unsorted_maxdev": broke_uns,
           "sorted_distinct_maxdev": ok_clean,
           "tied_breaks_ramp": broke_tied > 1e-6,
           "unsorted_breaks_ramp": broke_uns > 1e-6,
           "sorted_distinct_is_ramp": ok_clean < 1e-12}
    assert out["tied_breaks_ramp"], "self-test cannot fail: a tied vector gave the ramp"
    assert out["unsorted_breaks_ramp"], "self-test cannot fail: an unsorted vector gave the ramp"
    assert out["sorted_distinct_is_ramp"], "a sorted distinct vector did NOT give the ramp"
    return out


# ------------------------------------------------------------------ tie blocks
def blocks_of(v):
    """Run-lengths of EXACTLY equal consecutive values in `v` (v assumed non-decreasing)."""
    v = np.asarray(v, float)
    cuts = np.flatnonzero(v[1:] != v[:-1]) + 1
    edges = np.concatenate(([0], cuts, [len(v)]))
    return np.diff(edges)


def E_from_blocks(sizes):
    """The standardised rank vector implied by a block-size sequence (midranks)."""
    r = np.empty(int(np.sum(sizes)), float)
    p = 0
    for b in np.asarray(sizes, int):
        r[p:p + b] = p + (b + 1.0) / 2.0                                # 1-indexed midrank
        p += b
    return (r - r.mean()) / max(r.std(), 1e-12)


# ------------------------------------------------------------------ readouts
def readouts(W, Pt, p, nat):
    """(convex cloud RMSD, selection cloud RMSD, selected index) for one weight vector."""
    geo = core.backend("geometry")
    b = consensus_medoid(Pt, p)
    Sup = cc.superpose_batch(W, W[b])
    ww = np.asarray(p, float)
    ww = ww / max(float(ww.sum()), 1e-12)
    C = np.tensordot(ww, Sup, axes=(0, 0))
    return (float(geo.ca_rmsd(C, nat)), float(geo.ca_rmsd(W[b], nat)), int(b), C)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--draws", type=int, default=8)
    args = ap.parse_args()

    st = self_test()
    print("[self-test]", json.dumps(st), flush=True)

    geo = core.backend("geometry")
    names = sorted(f[:-4] for f in os.listdir(CACHE) if f.endswith(".npz"))
    if args.limit:
        names = names[:args.limit]

    rows = []
    t0 = time.time()
    for idx, pdb in enumerate(names):
        z = np.load(os.path.join(CACHE, "%s.npz" % pdb))
        order = np.asarray(z["order"], int)
        dis = np.asarray(z["dis"], float)
        W = np.asarray(z["W"], float)
        nat = np.asarray(z["nat"], float)
        fold = int(z["fold"])
        o = order[:NQ]
        sc = dis[o]

        # ---- Q0-H1: structure
        mono = bool(np.all(np.diff(sc) >= 0.0))
        E_tie = _zrank(sc)
        E_ramp = _zrank(np.arange(NQ, dtype=float))
        sizes = blocks_of(sc)
        E_blk = E_from_blocks(sizes)
        n_tied_pos = int(NQ - len(sizes))                              # elements in excess of blocks
        maxdev = float(np.abs(E_tie - E_ramp).max())
        blk_err = float(np.abs(E_tie - E_blk).max())                   # must be ~0 if H1 holds

        # ---- Q0-H2: the other paths
        alpha, T = VQE_LFO[fold % len(VQE_LFO)]
        len_top = int(len(order))

        # ---- Q0-H3: pricing
        Wo = W[o]
        Pt = geo.pairwise_ca_rmsd(Wo)
        p_tie, _s, _f = P.closed_form(E_tie, alpha, T)
        p_ramp, _s, _f = P.closed_form(E_ramp, alpha, T)
        cv_t, sel_t, b_t, C_t = readouts(Wo, Pt, p_tie, nat)
        cv_r, sel_r, b_r, C_r = readouts(Wo, Pt, p_ramp, nat)

        # matched-space control: the SAME block sizes, relocated along the ramp
        rng = np.random.default_rng(abs(hash(("s32Q0", pdb))) % (2 ** 32))
        perm_cv, perm_sel, perm_dev = [], [], []
        for _d in range(args.draws):
            szp = rng.permutation(sizes)
            E_p = E_from_blocks(szp)
            p_p, _s, _f = P.closed_form(E_p, alpha, T)
            cvp, selp, _bp, _Cp = readouts(Wo, Pt, p_p, nat)
            perm_cv.append(cvp)
            perm_sel.append(selp)
            perm_dev.append(float(np.abs(E_p - E_ramp).max()))

        rows.append({
            "pdb": pdb, "fold": fold, "alpha": float(alpha), "T": float(T),
            "n": int(z["n"]), "k": int(z["k"]), "len_order": len_top,
            "monotone": mono, "n_blocks": int(len(sizes)), "n_tied_excess": n_tied_pos,
            "maxdev_ramp": maxdev, "blockmodel_err": blk_err,
            "tv_p": float(0.5 * np.abs(p_tie - p_ramp).sum()),
            "cloud_tie": cv_t, "cloud_ramp": cv_r,
            "sel_tie": sel_t, "sel_ramp": sel_r,
            "selidx_tie": b_t, "selidx_ramp": b_r,
            "cloud_move": float(np.sqrt(((C_t - C_r) ** 2).sum() / len(C_t))),
            "perm_cloud_mean": float(np.mean(perm_cv)),
            "perm_cloud_sd": float(np.std(perm_cv, ddof=1)) if len(perm_cv) > 1 else 0.0,
            "perm_cloud_draws": [float(x) for x in perm_cv],
            "perm_sel_mean": float(np.mean(perm_sel)),
            "perm_maxdev_mean": float(np.mean(perm_dev)),
        })
        if idx % 20 == 0:
            print("  %3d/%d %s  %.1fs" % (idx + 1, len(names), pdb, time.time() - t0), flush=True)

    os.makedirs(OUT, exist_ok=True)
    tmp = os.path.join(OUT, "s32_Q0_rows.jsonl.%d.tmp" % os.getpid())
    with open(tmp, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    os.replace(tmp, os.path.join(OUT, "s32_Q0_rows.jsonl"))

    # ---------------------------------------------------------- verdicts
    sys.path.insert(0, os.path.join(ROOT, "s24"))
    import stats_lib as S

    nm = [r["pdb"] for r in rows]
    fo = [r["fold"] for r in rows]
    cloud_t = np.array([r["cloud_tie"] for r in rows])
    cloud_r = np.array([r["cloud_ramp"] for r in rows])
    sel_t = np.array([r["sel_tie"] for r in rows])
    sel_r = np.array([r["sel_ramp"] for r in rows])
    perm_m = np.array([r["perm_cloud_mean"] for r in rows])

    cmp_cloud = S.compare(cloud_t, cloud_r, folds=fo, names=nm, label="Q0_TIE-Q0_RAMP_cloud")
    cmp_sel = S.compare(sel_t, sel_r, folds=fo, names=nm, label="Q0_TIE_SEL-Q0_RAMP_SEL_cloud")
    cmp_perm = S.compare(perm_m, cloud_r, folds=fo, names=nm, label="Q0_PERM-Q0_RAMP_cloud")

    summary = {
        "prereg_commit": "a8f9d6a7",
        "self_test": st,
        "n": len(rows),
        "H1": {
            "monotone_all": bool(all(r["monotone"] for r in rows)),
            "n_with_ties": int(sum(1 for r in rows if r["n_tied_excess"] > 0)),
            "n_exactly_ramp": int(sum(1 for r in rows if r["maxdev_ramp"] < 1e-12)),
            "maxdev_max": float(max(r["maxdev_ramp"] for r in rows)),
            "maxdev_mean": float(np.mean([r["maxdev_ramp"] for r in rows])),
            "blockmodel_err_max": float(max(r["blockmodel_err"] for r in rows)),
            "tied_excess_mean": float(np.mean([r["n_tied_excess"] for r in rows])),
            "tied_excess_max": int(max(r["n_tied_excess"] for r in rows)),
            "H1_holds": bool(all(r["monotone"] for r in rows)
                             and max(r["blockmodel_err"] for r in rows) < 1e-9
                             and all((r["maxdev_ramp"] > 1e-12) == (r["n_tied_excess"] > 0)
                                     for r in rows)),
        },
        "H2": {
            "len_order_min": int(min(r["len_order"] for r in rows)),
            "distinct_alphaT": sorted({(r["alpha"], r["T"]) for r in rows}),
            "dim_always_128": bool(all(r["len_order"] >= NQ for r in rows)),
        },
        "H3": {
            "tv_p_mean": float(np.mean([r["tv_p"] for r in rows])),
            "tv_p_max": float(max(r["tv_p"] for r in rows)),
            "sel_same": int(sum(1 for r in rows if r["selidx_tie"] == r["selidx_ramp"])),
            "cloud_move_mean": float(np.mean([r["cloud_move"] for r in rows])),
            "cloud_move_max": float(max(r["cloud_move"] for r in rows)),
            "cmp_cloud": cmp_cloud, "cmp_sel": cmp_sel, "cmp_perm": cmp_perm,
            "perm_draw_sd_mean": float(np.mean([r["perm_cloud_sd"] for r in rows])),
        },
        "basis": "CA point cloud (diagnostic); built chain NOT run -- see verdict",
        "elapsed_s": time.time() - t0,
    }
    tmp = os.path.join(OUT, "s32_Q0_invariance.json.%d.tmp" % os.getpid())
    with open(tmp, "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    os.replace(tmp, os.path.join(OUT, "s32_Q0_invariance.json"))

    print(json.dumps({k: v for k, v in summary.items() if k != "H3"}, indent=2, default=float))
    for key in ("cmp_cloud", "cmp_sel", "cmp_perm"):
        c = summary["H3"][key]
        print("%-34s eff %+.4f  SE %.4f  MDE %.4f  %.2fx  %dW/%dL/%dT  median %+.4f"
              % (c["label"], c["effect"], c["se"], c["mde"], c["effect_over_mde"],
                 c["n_better"], c["n_worse"], c["n_tied"], c["median_effect"]))
    print("sel same on %d/%d; TV(p) mean %.3e max %.3e"
          % (summary["H3"]["sel_same"], len(rows),
             summary["H3"]["tv_p_mean"], summary["H3"]["tv_p_max"]))


if __name__ == "__main__":
    main()
