#!/usr/bin/env python
"""s30/s30_Q_analyse.py -- lane Q, analysis of the sparse-readout ceiling-vs-bits measurement.

Reads `s30/results/s30_Q_sparse_rows.jsonl` (126 targets) and answers F1-F4 of
`s30/PREREG_S30_Q_sparse.md`. POINT CLOUD basis throughout; every printed figure says so.
The ORACLE/ORACLE corner is cross-checked against S29 lane O's already-chained rows.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s30.s30_Q_sparse import log2C          # noqa: E402

ROWS = os.path.join(HERE, "results", "s30_Q_sparse_rows.jsonl")
OUT = os.path.join(HERE, "results", "s30_Q_sparse.json")
O_CHAIN = os.path.join(ROOT, "s29", "results")


def load_rows():
    out = {}
    for L in open(ROWS):
        L = L.strip()
        if L:
            r = json.loads(L)
            out[r["pdb"]] = r
    return out


def load_o_chain():
    """S29 lane O's paired (cloud, chain) arm values -- the calibration set and the anchor."""
    pairs = {}
    for f in glob.glob(os.path.join(O_CHAIN, "s29_O_chain_rows*.jsonl")):
        for L in open(f):
            L = L.strip()
            if not L:
                continue
            r = json.loads(L)
            pairs[(r["pdb"], r["item"])] = (r["rmsd_cloud"], r["rmsd_chain"])
    return pairs


def chain_price_fit(pairs):
    """The projection price as a function of cloud accuracy (S29 D5). Returns (a, b) for
    chain - cloud = a + b * cloud, fitted over every paired arm, and the fit's R^2."""
    c = np.array([v[0] for v in pairs.values()])
    ch = np.array([v[1] for v in pairs.values()])
    d = ch - c
    A = np.vstack([np.ones_like(c), c]).T
    coef, *_ = np.linalg.lstsq(A, d, rcond=None)
    pred = A @ coef
    ss = 1.0 - ((d - pred) ** 2).sum() / ((d - d.mean()) ** 2).sum()
    return float(coef[0]), float(coef[1]), float(ss), len(c)


def main():
    rows = load_rows()
    pdbs = [t["pdb"] for t in I.targets() if t["pdb"] in rows]
    n = len(pdbs)
    folds = ST.pinned_folds(pdbs)
    fm = np.array([p in I.FAIL18 for p in pdbs])
    print("== lane Q sparse-readout ceiling. n = %d targets. POINT CLOUD basis." % n)
    print("   EVERY arm whose name is not D_* or H_* reads the native and is ORACLE.\n")

    names = sorted({a for p in pdbs for a in rows[p]["arms"]})
    V, BITS = {}, {}
    for a in names:
        v = np.array([rows[p]["arms"][a]["rmsd"] for p in pdbs])
        V[a] = v
        sb = rows[pdbs[0]]["arms"][a].get("sup_bits", 0.0)
        wb = rows[pdbs[0]]["arms"][a].get("w_bits", 0.0)
        BITS[a] = (float(sb), float(wb))
    prod = V["prod"]
    ref = np.array([[rows[p]["argmin_ref"][str(B)] for B in range(10)] for p in pdbs])  # (n,10)

    def se(x):
        return float(np.std(x, ddof=1) / math.sqrt(len(x)))

    def line(a):
        v = V[a]; sb, wb = BITS[a]
        tb = sb + wb
        return dict(arm=a, mean=float(v.mean()), se=se(v), sup_bits=sb, w_bits=wb,
                    total_bits=(None if not np.isfinite(tb) else float(tb)),
                    fail18=float(v[fm].mean()), other108=float(v[~fm].mean()))

    out = dict(n=n, arms={a: line(a) for a in names}, provenance=ST.provenance(__file__))

    # ---------------- the argmin reference curve (the budget this family must beat)
    print("-- REFERENCE: the deployed readout's own curve. argmin over the top-2**B of the score")
    print("   order costs EXACTLY B bits (point cloud, ORACLE readout).")
    out["argmin_ref"] = {}
    for B in range(10):
        m = float(ref[:, B].mean())
        out["argmin_ref"][B] = dict(mean=m, se=se(ref[:, B]), bits=float(B))
        print("     B = %d bits   top-%-3d   %.4f A" % (B, 2 ** B, m))
    print("     (production, the m = 75 uniform average: %.4f A)\n" % prod.mean())

    # ---------------- F1: does the SUPPORT carry the gain?
    print("-- F1: ORACLE support + UNIFORM weights, as a share of the ORACLE/ORACLE gain")
    out["F1"] = {}
    for S in ("pool", "top75"):
        for s in (2, 3, 5, 10, 20):
            A, Bv = V.get("A_%s_s%d" % (S, s)), V.get("B_%s_s%d" % (S, s))
            if A is None:
                continue
            gain = prod.mean() - A.mean()
            got = prod.mean() - Bv.mean()
            frac = got / gain if gain > 0 else float("nan")
            out["F1"]["%s_s%d" % (S, s)] = dict(A=float(A.mean()), B=float(Bv.mean()),
                                                gain=float(gain), recovered=float(got),
                                                share=float(frac))
            print("     %-6s s=%-2d  A(oracle w) %.4f   B(uniform w) %.4f   recovers %5.1f%% of %.4f A"
                  % (S, s, A.mean(), Bv.mean(), 100 * frac, gain))
    f1 = out["F1"]["pool_s2"]["share"]
    out["F1_verdict"] = ("CONFIRMED" if f1 >= 0.50 else "REFUTED")
    print("   F1 (s=2 pool, threshold 50%%): %.1f%% -> %s\n" % (100 * f1, out["F1_verdict"]))

    # ---------------- F2: how many bits do the WEIGHTS cost?
    print("-- F2: ORACLE support, weights on the 1/L simplex lattice (L levels)")
    out["F2"] = {}
    for S in ("pool",):
        for s in (2, 3, 5, 10, 20):
            A = V.get("A_%s_s%d" % (S, s))
            if A is None:
                continue
            rec = {}
            for L in (1, 2, 3, 4, 8, 16, 32):
                q = V.get("Q_%s_s%d_L%d" % (S, s, L))
                if q is None:
                    continue
                rec[L] = dict(mean=float(q.mean()), gap=float(q.mean() - A.mean()),
                              w_bits=BITS["Q_%s_s%d_L%d" % (S, s, L)][1])
            out["F2"]["%s_s%d" % (S, s)] = dict(A=float(A.mean()), grid=rec)
            good = [L for L, d in sorted(rec.items()) if d["gap"] <= 0.05]
            print("     %-6s s=%-2d  continuous %.4f | %s | first L within 0.05 A: %s"
                  % (S, s, A.mean(),
                     "  ".join("L%d %.3f(+%.3f,%.1fb)" % (L, d["mean"], d["gap"], d["w_bits"])
                               for L, d in sorted(rec.items())),
                     (good[0] if good else "none")))
    g2 = out["F2"]["pool_s2"]["grid"]
    first = [L for L in sorted(g2) if g2[L]["gap"] <= 0.05]
    out["F2_verdict"] = ("CONFIRMED" if (first and first[0] <= 4) else "REFUTED")
    print("   F2 (s=2 pool, within 0.05 A at L<=4): first L = %s -> %s\n"
          % (first[0] if first else "none", out["F2_verdict"]))

    # ---------------- F3: is there a NATIVE-FREE support rule with skill?
    print("-- F3: how the support is chosen, with ORACLE weights held on (point cloud)")
    out["F3"] = {}
    for S in ("pool",):
        for s in (2, 3, 5, 10, 20):
            A = V.get("A_%s_s%d" % (S, s))
            if A is None:
                continue
            C, G, E = V["C_%s_s%d" % (S, s)], V["G_%s_s%d" % (S, s)], V["E_%s_s%d" % (S, s)]
            best_nf = min(C.mean(), G.mean())
            out["F3"]["%s_s%d" % (S, s)] = dict(
                A=float(A.mean()), C_scoreprefix=float(C.mean()), G_diversity=float(G.mean()),
                E_random_mean=float(E.mean()), best_native_free=float(best_nf),
                nf_minus_random=float(best_nf - E.mean()), oracle_minus_nf=float(best_nf - A.mean()))
            print("     %-6s s=%-2d  ORACLE sup %.4f | score-prefix %.4f  diversity %.4f  "
                  "RANDOM(mean of 8) %.4f | best native-free beats random by %+.4f"
                  % (S, s, A.mean(), C.mean(), G.mean(), E.mean(), E.mean() - best_nf))
    d3 = out["F3"]["pool_s2"]["nf_minus_random"]
    out["F3_verdict"] = ("CONFIRMED" if (-d3) >= 0.20 else "REFUTED")
    print("   F3 (s=2 pool, native-free beats random by >=0.20 A): %+.4f -> %s\n"
          % (-d3, out["F3_verdict"]))

    # ---------------- the deployable (0-bit, native-free-construction) rows
    print("-- DEPLOYABLE rows: construction reads NO native (only the RMSD readout is ORACLE)")
    out["deployable"] = {}
    for S in ("pool",):
        for s in (2, 3, 5, 10, 20):
            D, H = V.get("D_%s_s%d" % (S, s)), V.get("H_%s_s%d" % (S, s))
            if D is None:
                continue
            out["deployable"]["%s_s%d" % (S, s)] = dict(score_prefix_uniform=float(D.mean()),
                                                        diversity_uniform=float(H.mean()))
            print("     %-6s s=%-2d  top-s average %.4f   diversity-s average %.4f   (production %.4f)"
                  % (S, s, D.mean(), H.mean(), prod.mean()))
    print()

    # ---------------- F4: the frontier. ceiling vs TOTAL bits, against the argmin at equal bits
    print("-- F4 (DECISIVE): fully-priced sparse arms (support bits + weight bits, no continuous")
    print("   weights) against the argmin readout at the SAME budget. Point cloud.")
    priced = [(a, BITS[a][0] + BITS[a][1], V[a].mean()) for a in names
              if np.isfinite(BITS[a][0] + BITS[a][1]) and not a.startswith(("E_", "F_", "prod"))]
    out["F4"] = {}
    wins = []
    for B in range(0, 21):
        cands = [(m, a, b) for (a, b, m) in priced if b <= B + 1e-9]
        if not cands:
            continue
        m, a, b = min(cands)
        r = float(ref[:, B].mean()) if B < 10 else None
        rec = dict(budget=B, best_arm=a, arm_bits=float(b), arm_mean=float(m), argmin_mean=r,
                   sparse_minus_argmin=(None if r is None else float(m - r)))
        out["F4"][B] = rec
        if r is not None:
            flag = "SPARSE WINS" if m < r - 1e-9 else ""
            if m < r - 1e-9:
                wins.append(B)
            print("     B = %2d bits | best priced sparse %-16s (%.1f b) %.4f | argmin %.4f | %+.4f %s"
                  % (B, a, b, m, r, m - r, flag))
        else:
            print("     B = %2d bits | best priced sparse %-16s (%.1f b) %.4f | argmin  n/a (pool is 2**8.97)"
                  % (B, a, b, m))
    # the sprint's governing question: does the sparse family behave differently on the TAIL?
    print("\n   The same frontier on the FAIL18 tail and the other 108 (S30's governing split):")
    out["F4_strata"] = {}
    for B in (5, 7, 9):
        cands = [(m, a, b) for (a, b, m) in priced if b <= B + 1e-9]
        if not cands:
            continue
        _, a, b = min(cands)
        sa, so = float(V[a][fm].mean()), float(V[a][~fm].mean())
        ra, ro = float(ref[fm, B].mean()), float(ref[~fm, B].mean())
        out["F4_strata"][B] = dict(arm=a, sparse_fail18=sa, argmin_fail18=ra,
                                   sparse_108=so, argmin_108=ro,
                                   diff_fail18=sa - ra, diff_108=so - ro)
        print("     B = %d | FAIL18: sparse %.4f vs argmin %.4f (%+.4f) | 108: sparse %.4f vs "
              "argmin %.4f (%+.4f)" % (B, sa, ra, sa - ra, so, ro, so - ro))
    print()

    out["F4_wins_at"] = wins
    out["F4_verdict"] = ("CONFIRMED" if any(4 <= b <= 12 for b in wins) else "REFUTED")
    print("   F4 (a priced sparse arm beats the argmin at equal bits, B in 4..12): %s -> %s\n"
          % (wins if wins else "never", out["F4_verdict"]))

    # ---------------- F4 ON THE ENDPOINT: lane O's ALREADY-CHAINED arms, no new compute
    pairs = load_o_chain()
    print("-- F4 ON THE BUILT CHAIN (the endpoint, contract rule 1). These are S29 lane O's")
    print("   already-chained arms; nothing here is a projection. Bits by the same ledger.")
    chain_arms = [
        ("best1_top75",     "argmin over top-75",            math.log2(75),          None),
        ("best1_top128",    "argmin over top-128",           7.0,                    None),
        ("best1_pool",      "argmin over the 500",           math.log2(500),         None),
        ("sparse_top75_s2", "2 of top-75 + CONTINUOUS w",    log2C(75, 2),           float("inf")),
        ("sparse_top75_s5", "5 of top-75 + CONTINUOUS w",    log2C(75, 5),           float("inf")),
        ("sparse_pool_s2",  "2 of 500 + CONTINUOUS w",       log2C(500, 2),          float("inf")),
        ("sparse_pool_s5",  "5 of 500 + CONTINUOUS w",       log2C(500, 5),          float("inf")),
        ("sparse_pool_s20", "20 of 500 + CONTINUOUS w",      log2C(500, 20),         float("inf")),
        ("bestm128",        "prefix size m in top-128",      7.0,                    None),
        ("prod",            "production, m = 75 uniform",    0.0,                    None),
    ]
    ch = {}
    for (p, it), (c, k) in pairs.items():
        ch.setdefault(it, []).append(k)
    out["F4_chain"] = {}
    print("     %-17s  %-30s  %8s  %s" % ("arm", "what it is", "bits", "BUILT CHAIN (A)"))
    for nm, desc, sb, wb in chain_arms:
        if nm not in ch or len(ch[nm]) < 100:
            continue
        m = float(np.mean(ch[nm]))
        tb = sb if wb is None else float("inf")
        out["F4_chain"][nm] = dict(desc=desc, sup_bits=float(sb), w_bits=(None if wb is None else "inf"),
                                   chain=m, n=len(ch[nm]))
        print("     %-17s  %-30s  %8s  %.4f" % (nm, desc, ("%.1f" % sb) if wb is None else "%.1f+inf" % sb, m))
    if "sparse_top75_s2" in ch and "best1_top128" in ch:
        d = float(np.mean(ch["sparse_top75_s2"])) - float(np.mean(ch["best1_top128"]))
        out["F4_chain_headline"] = dict(
            sparse_top75_s2=float(np.mean(ch["sparse_top75_s2"])), sparse_bits=log2C(75, 2),
            best1_top128=float(np.mean(ch["best1_top128"])), argmin_bits=7.0, diff=d)
        print("\n     HEADLINE: 2 of top-75 with FREE CONTINUOUS weights (%.1f bits + unbounded) is"
              % log2C(75, 2))
        print("     %.4f A, while the argmin over the top-128 (7.0 bits, no weights) is %.4f A."
              % (np.mean(ch["sparse_top75_s2"]), np.mean(ch["best1_top128"])))
        print("     The argmin WINS by %.4f A using %.1f FEWER bits and no weight channel at all.\n"
              % (d, log2C(75, 2) - 7.0))

    # ---------------- order-statistic pricing of the random-support arm
    print("-- ORDER-STATISTIC CONTROL (contract rule 8): the RANDOM-support arm is reported as the")
    print("   MEAN of 8 draws; its per-target MIN is priced here and never used as an arm.")
    out["order_stat"] = {}
    for S in ("pool",):
        for s in (2, 5):
            M = np.array([rows[p]["arms"]["E_%s_s%d" % (S, s)]["draws"] for p in pdbs])
            bk = ST.best_of_k_within(M, seed_parts=("s30Q", "rand", S, str(s)))
            out["order_stat"]["%s_s%d" % (S, s)] = bk
            print("     %-6s s=%-2d  mean-of-8 %.4f  min-of-8 %.4f  observed gain %+.4f, "
                  "across-target null %+.4f (%.0f%% accounted), split-half %+.4f (%.0f%%)"
                  % (S, s, M.mean(), M.min(1).mean(), bk["observed_gain"], bk["null_across_targets"],
                     100 * bk["share_accounted"], bk["split_half"], 100 * bk["split_half_frac"]))
    print()

    # ---------------- greedy-gap control
    ex = V.get("EX64_s2_cont"); gr = V.get("T64_s2_cont")
    if ex is not None and gr is not None:
        out["greedy_gap"] = dict(exhaustive=float(ex.mean()), greedy=float(gr.mean()),
                                 gap=float(gr.mean() - ex.mean()))
        print("-- GREEDY-GAP CONTROL at T = 64, s = 2 (C(64,2) = 2016 pairs enumerated):")
        print("     exhaustive %.4f   greedy %.4f   gap %+.4f A  -> greedy arms are %s\n"
              % (ex.mean(), gr.mean(), gr.mean() - ex.mean(),
                 "faithful" if abs(gr.mean() - ex.mean()) <= 0.05 else "LOWER BOUNDS on the ceiling"))

    # ---------------- chain calibration (the endpoint is the built chain)
    a0, b0, r2, npair = chain_price_fit(pairs)
    out["chain_fit"] = dict(intercept=a0, slope=b0, r2=r2, n=npair)
    print("-- CHAIN PROJECTION (S29 D5). Fitted on %d paired lane-O arm values:" % npair)
    print("     chain - cloud = %+.4f %+.4f * cloud   (R^2 %.3f)" % (a0, b0, r2))
    print("     ESTIMATES ONLY -- the built chain is the endpoint and the decisive arms are chained.")

    def to_chain(m):
        return m + a0 + b0 * m
    key_arms = ["prod", "A_pool_s2", "B_pool_s2", "Q_pool_s2_L4", "C_pool_s2", "D_pool_s2",
                "G_pool_s2", "H_pool_s2", "E_pool_s2"]
    out["chain_est"] = {}
    for a in key_arms:
        if a in V:
            out["chain_est"][a] = dict(cloud=float(V[a].mean()), chain_est=float(to_chain(V[a].mean())))
            print("     %-16s cloud %.4f -> chain est %.4f" % (a, V[a].mean(), to_chain(V[a].mean())))

    ST.save_atomic(OUT, out, module_file=__file__)
    print("\nwrote", OUT)
    return out


if __name__ == "__main__":
    main()
