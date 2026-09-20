#!/usr/bin/env python
"""s29/s29_D_bcheck.py -- lane D: the adversary check of lane B's S29-L25 (the exhaustive
pair search over the deployed pool, and its claim that the f-optimal subset is not a prefix).

Written from the ENTRY's description, not from `s29/s29_B_tta.py`'s flow: my pair scan is a
single vectorised pass with no chunking, my tie handling averages the OUTCOME over the tied
argmin set (`ST.argmin_tied`'s rule), and I evaluate f under BOTH readings of the shipped risk
table -- the piecewise-linear `Surrogate` lane B used and the EXACT shipped lookup
(`s12.instrument.shipped_score`) -- because the two differ by about 0.006 to 0.009 in absolute
value and three of lane B's twelve gaps are smaller than that.

Checks, in the coordinator's order:
  A. is the search exhaustive?             C(N, 2) pairs enumerated and counted against the entry
  B. tie handling                          the tied argmin SET's size, and whether the prefix is
                                           inside it; lane B takes `tie[0]`, the first in ARRAY
                                           order, which on a DIS-sorted pool is the best-ranked
                                           member and therefore biased TOWARD the prefix
  C. the frame                             `Frame(W, DIS top-75)` -- the deployed readout's own
                                           medoid frame -- asserted against `readout_uniform`
  D. f vs the shipped score                the gap recomputed under the exact lookup
  E. the 7JGX zero gap                     is the argmin genuinely the prefix, or an identity?

    python s29/s29_D_bcheck.py pairs [--targets 7JGX,5Z5W,...]
    python s29/s29_D_bcheck.py power

Results `s29/results/s29_D_bcheck_pairs.json`, `s29_D_bcheck_power.json`.
ORACLE only in the RMSD columns, which are attached after every selection, as lane B did.
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
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_A_amp as AMP           # noqa: E402

RESULTS = os.path.join(HERE, "results")
B_ROWS = os.path.join(RESULTS, "s29_B_tta_subset_rows.jsonl")
M_PROD = 75
#: lane B's 12 targets (S27 T11's trainability set), and its reported best pairs
B_TARGETS = ["1A13", "1I6Y", "1M02", "2BFI", "2LWS", "2MP9", "2P5H", "5Z5W", "6MBM", "7JGX", "8HVS", "9KAR"]


def surrogate_values(sur, Cb):
    """The piecewise-linear risk reading, batched: identical arithmetic to Surrogate.value_grad."""
    Cb = np.asarray(Cb, float)
    out = np.empty(len(Cb))
    ar = np.arange(sur.npairs)
    for a in range(0, len(Cb), 4096):
        b = Cb[a:a + 4096]
        d = np.linalg.norm(b[:, sur.i, :] - b[:, sur.j, :], axis=2)
        u = (d - sur.g0) / sur.dg
        i0 = np.clip(np.floor(u).astype(int), 0, sur.G - 2)
        f = u - i0
        r0 = sur.risk[ar[None, :], i0]
        sl = sur.slope[ar[None, :], i0]
        out[a:a + 4096] = (r0 + sl * f * sur.dg).mean(1)
    return out


def shipped_values(dg, i, j, Cb):
    """The EXACT shipped lookup (`I.shipped_score`), the production functional."""
    Cb = np.asarray(Cb, float)
    out = np.empty(len(Cb))
    for a in range(0, len(Cb), 4096):
        b = Cb[a:a + 4096]
        D = I.pair_dists(b, i, j)
        out[a:a + 4096] = np.asarray(I.shipped_score(dg, D.astype(np.float32).astype(float)), float)
    return out


def pair_check(pdb, n_pool=500):
    from s24 import d_harness as H
    from s27 import run_pool as RP
    t0 = time.time()
    cand, ch, _ = RP.channels_for(pdb)
    E_real = RP.zr(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    order = RP.topm(E_real, cand.k, key)
    top = np.sort(RP.topm(E_real, M_PROD, key)).astype(int)
    frame = AMP.Frame(cand.W, top)
    dg = I.distogram(pdb, cand.seq, cand.fold)
    sur = AMP.Surrogate(dg, int(cand.n))
    idx = np.asarray(order[:min(n_pool, cand.k)], int)
    Wp = frame.Wp[idx]
    N = len(idx)
    # C: the frame is the deployed readout's own
    Cd, _ = H.readout_uniform(cand, top)
    Cframe = frame.Wp[top].mean(0)
    frame_dev = float(np.abs(np.asarray(Cd, float) - Cframe).max())
    # A: exhaustive enumeration, ONE pass, no chunk-local argmin
    ii, jj = np.triu_indices(N, 1)
    n_pairs = int(len(ii))
    mid = 0.5 * (Wp[ii] + Wp[jj])
    v_sur = surrogate_values(sur, mid)
    v_shp = shipped_values(dg, sur.i, sur.j, mid)
    out = dict(pdb=pdb, n=int(cand.n), N=N, n_pairs=n_pairs,
               n_pairs_expected=int(N * (N - 1) // 2), exhaustive=bool(n_pairs == N * (N - 1) // 2),
               frame_max_dev_vs_deployed_average=frame_dev, secs_scan=time.time() - t0)
    for tag, v in (("surrogate", v_sur), ("shipped", v_shp)):
        mn = float(v.min())
        tie = np.flatnonzero(v <= mn + 1e-15)
        prefix_v = float(v[(ii == 0) & (jj == 1)][0])
        best_ranks = sorted((int(ii[tie[0]]), int(jj[tie[0]])))
        # the OUTCOME averaged over the tied argmin set (ST.argmin_tied's rule)
        rr = np.asarray(cand.oracle_rr, float)                       # ORACLE, attached after
        Cb_tie = 0.5 * (Wp[ii[tie]] + Wp[jj[tie]])
        rmsd_tie = np.array([I.ca_rmsd(c, cand.nat_ca) for c in Cb_tie])   # ORACLE
        out[tag] = dict(
            best_value=mn, prefix_value=prefix_v, gap=float(prefix_v - mn),
            n_tied=int(len(tie)), best_ranks_array_order=best_ranks,
            tie_contains_prefix=bool(any((ii[t] == 0 and jj[t] == 1) for t in tie)),
            argmin_is_prefix=bool(best_ranks == [0, 1]),
            tie_rank_min=int(min(ii[tie].min(), jj[tie].min())),
            tie_rank_max=int(max(ii[tie].max(), jj[tie].max())),
            rmsd_over_tie_mean=float(rmsd_tie.mean()), rmsd_over_tie_sd=float(rmsd_tie.std()),
            rmsd_prefix=float(I.ca_rmsd(0.5 * (Wp[0] + Wp[1]), cand.nat_ca)))              # ORACLE
    # B: what lane B's chunk-local rule would have produced (its own arithmetic, reproduced)
    best_v, best_pair, CH = np.inf, (-1, -1), 20000
    for a in range(0, n_pairs, CH):
        sl = slice(a, a + CH)
        vals = v_sur[sl]
        t = int(np.argmin(vals))
        tie = np.flatnonzero(vals <= vals[t] + 1e-15)
        if vals[t] < best_v:
            best_v = float(vals[t]); best_pair = (int(ii[sl][tie[0]]), int(jj[sl][tie[0]]))
    out["laneB_chunk_rule"] = dict(value=best_v, pair=sorted(best_pair),
                                   agrees_with_single_pass=bool(abs(best_v - out["surrogate"]["best_value"]) < 1e-15))
    out["secs"] = time.time() - t0
    return out


def power_check():
    """Claim 2's price: what n would be needed to measure the direction of +0.4278 A at 0.59x MDE."""
    rows = [json.loads(l) for l in open(B_ROWS, encoding="utf-8") if l.strip()] if os.path.exists(B_ROWS) else []
    out = dict(check="lane B's claim 2 (the f-optimal m=5 subset's average vs production), power arithmetic",
               n_observed=12, effect=0.4278, over_mde=0.59)
    # from effect / MDE = 0.59 at n = 12: MDE = 2.8016 * sd / sqrt(n)  =>  sd
    mde = out["effect"] / out["over_mde"]
    sd = mde * np.sqrt(out["n_observed"]) / ST.MDE_K
    out["implied_sd"] = float(sd)
    out["implied_se_at_12"] = float(sd / np.sqrt(12))
    for n in (12, 24, 40, 60, 90, 126):
        se = sd / np.sqrt(n)
        rd = ST.retrodesign(out["effect"], se)
        out[f"n_{n}"] = dict(se=float(se), mde=float(ST.MDE_K * se), over_mde=float(out["effect"] / (ST.MDE_K * se)),
                             power=rd["power"], type_m=rd["type_m"], type_s=rd["type_s"])
    # the n at which the observed effect would clear its own MDE
    n_needed = int(np.ceil((ST.MDE_K * sd / out["effect"]) ** 2))
    out["n_to_clear_its_own_MDE"] = n_needed
    out["n_to_reach_power_0.8"] = int(np.ceil((2.8016 * sd / out["effect"]) ** 2))
    if rows:
        d = []
        for r in rows:
            for k in ("rmsd_greedy", "rmsd_sub5", "rmsd_m5"):
                if k in r:
                    d.append(float(r[k])); break
        out["rows_found"] = len(rows)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["pairs", "power"])
    ap.add_argument("--targets", default="7JGX,5Z5W,2BFI,1A13,2MP9,9KAR")
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    if a.mode == "power":
        res = power_check()
        print(json.dumps(res, indent=1))
        ST.save_atomic(os.path.join(RESULTS, "s29_D_bcheck_power.json"), res, module_file=__file__)
        return
    rows = []
    for pdb in a.targets.split(","):
        r = pair_check(pdb)
        rows.append(r)
        s, sh = r["surrogate"], r["shipped"]
        print(f"  {pdb}: {r['n_pairs']} pairs (exhaustive {r['exhaustive']}), frame dev {r['frame_max_dev_vs_deployed_average']:.2e}")
        print(f"     SURROGATE best {s['best_ranks_array_order']} value {s['best_value']:.4f} prefix {s['prefix_value']:.4f} "
              f"gap {s['gap']:+.4f}  ties {s['n_tied']}  argmin is prefix: {s['argmin_is_prefix']}")
        print(f"     SHIPPED   best {sh['best_ranks_array_order']} value {sh['best_value']:.4f} prefix {sh['prefix_value']:.4f} "
              f"gap {sh['gap']:+.4f}  ties {sh['n_tied']}  argmin is prefix: {sh['argmin_is_prefix']}")
        print(f"     lane B's chunk rule: pair {r['laneB_chunk_rule']['pair']} value {r['laneB_chunk_rule']['value']:.4f} "
              f"(agrees with a single pass: {r['laneB_chunk_rule']['agrees_with_single_pass']})   {r['secs']:.1f}s", flush=True)
    res = dict(check="adversary check of lane B's S29-L25 pair search", rows=rows,
               n_targets=len(rows),
               all_exhaustive=bool(all(r["exhaustive"] for r in rows)),
               frame_max_dev=float(max(r["frame_max_dev_vs_deployed_average"] for r in rows)),
               n_argmin_is_prefix_surrogate=int(sum(r["surrogate"]["argmin_is_prefix"] for r in rows)),
               n_argmin_is_prefix_shipped=int(sum(r["shipped"]["argmin_is_prefix"] for r in rows)),
               n_sign_flips_between_readings=int(sum(r["surrogate"]["argmin_is_prefix"] != r["shipped"]["argmin_is_prefix"]
                                                     for r in rows)))
    ST.save_atomic(os.path.join(RESULTS, "s29_D_bcheck_pairs.json"), res, module_file=__file__)
    print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
