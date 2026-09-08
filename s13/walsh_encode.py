"""SPRINT 13 PAULI-SPECTRUM, EXPERIMENT W4 -- ENCODING DEPENDENCE.

The energy is a function on the k-ary configuration space.  Its RESIDUE-order spectrum D_d
is a property of that function and no index encoding can change it.  Its PAULI-WEIGHT
spectrum V_w is a property of the function COMPOSED WITH THE ENCODING, and different
encodings of the same physics give different weights.  This experiment measures how much.

Encodings compared, all on identical energy tables:

  binary       state j -> the b-bit binary numeral j
  gray         state j -> gray(j) = j ^ (j>>1)
  prior_gray   states first sorted by the leakage-safe empirical prior occupancy, then gray
  energy_gray  states first sorted by the residue's MARGINAL mean energy, then gray
  perm_rand    200 draws of an independent uniform permutation per residue (the null)
  perm_best/worst  the best/worst of an exhaustive search over the k! single permutation
               applied to every residue (k=4 -> 24 candidates)
  onehot       L*k qubits, one bit per (residue, state); the physical energy is defined ONLY
               on the feasible subspace, so TWO extensions are reported --
               "mean" (infeasible -> mean feasible energy) and "penalty" (mean + lambda *
               sum_r (popcount_r - 1)^2).  The dependence on that choice is the point.

    python -m s13.walsh_encode
"""
from __future__ import annotations

import itertools
import sys
import time

import numpy as np

from s13 import walsh_lib as W
from s13 import qarch_lib as Q

TARGETS = ["1A13", "1A1P", "2BFI", "1DEP", "1ID6", "1CEK"]
CELLS = [(6, 4), (5, 4), (4, 8)]
ONEHOT_CELLS = [(4, 4), (5, 4)]


def relabel_spectrum(T, perms, m, b):
    return W.spectrum(W.apply_relabel(T, perms).ravel(), m, bits_per_res=b)


def cell(pdb_id, L, k, rng):
    sp = W.space_for(pdb_id, L, k)
    b, m = sp.bits_per_res, L * sp.bits_per_res
    S = W.enumerate_states(L, k)
    _, tot = W.legacy_tables(sp, S)
    T = np.asarray(tot, float).reshape((k,) * L)

    P = Q.empirical_prior(sp)                       # (L, k), leakage-safe, no native
    marg = np.array([[T.take(j, axis=r).mean() for j in range(k)] for r in range(L)])

    ident = [np.arange(k) for _ in range(L)]
    enc = {}
    enc["binary"] = ident
    enc["gray"] = [W.perm_gray(k) for _ in range(L)]
    # sort states by prior occupancy (descending) then lay a gray code over them
    g = W.perm_gray(k)
    enc["prior_gray"] = [g[np.argsort(np.argsort(-P[r]))] for r in range(L)]
    enc["energy_gray"] = [g[np.argsort(np.argsort(marg[r]))] for r in range(L)]

    row = {"pdb": pdb_id, "L": L, "k": k, "m_bits": m, "bits_per_res": b, "enc": {}}
    for name, perms in enc.items():
        s = relabel_spectrum(T, perms, m, b)
        row["enc"][name] = {kk: s[kk] for kk in
                            ("W_mean", "W_eff", "V1", "V2", "V3", "tail_gt3",
                             "D_mean", "D1", "D_tail_gt2", "spread", "var",
                             "parseval_rel_err", "order_weight_bound_holds")}
        row["enc"][name]["V"] = s["V"]
        row["enc"][name]["D"] = s["D"]

    # ---- the null: independent random relabellings
    Ws, V1s = [], []
    for _ in range(200):
        perms = [rng.permutation(k) for _ in range(L)]
        s = relabel_spectrum(T, perms, m, b)
        Ws.append(s["W_mean"]); V1s.append(s["V1"])
    row["perm_rand"] = {"n": 200, "W_mean_mean": float(np.mean(Ws)),
                        "W_mean_sd": float(np.std(Ws, ddof=1)),
                        "W_mean_min": float(np.min(Ws)), "W_mean_max": float(np.max(Ws)),
                        "V1_mean": float(np.mean(V1s)), "V1_max": float(np.max(V1s))}

    # ---- exhaustive over the single permutation applied to every residue
    best = None; worst = None
    for p in itertools.permutations(range(k)):
        s = relabel_spectrum(T, [np.array(p)] * L, m, b)
        c = {"perm": list(p), "W_mean": s["W_mean"], "V1": s["V1"],
             "tail_gt3": s["tail_gt3"]}
        if best is None or c["W_mean"] < best["W_mean"]:
            best = c
        if worst is None or c["W_mean"] > worst["W_mean"]:
            worst = c
    row["perm_uniform_best"] = best
    row["perm_uniform_worst"] = worst
    # invariance check: D must be identical across every encoding
    Ds = np.array([row["enc"][n]["D"] for n in row["enc"]])
    row["D_invariance_max_abs_dev"] = float(np.abs(Ds - Ds[0]).max())
    print(f"  {pdb_id} L={L} k={k} m={m}: "
          + "  ".join(f"{n} {row['enc'][n]['W_mean']:.3f}" for n in row["enc"])
          + f"  rand {row['perm_rand']['W_mean_mean']:.3f}"
          f"+-{row['perm_rand']['W_mean_sd']:.3f}"
          f"  [{best['W_mean']:.3f},{worst['W_mean']:.3f}]"
          f"  Dinv {row['D_invariance_max_abs_dev']:.2e}", flush=True)
    return row


def onehot_cell(pdb_id, L, k):
    sp = W.space_for(pdb_id, L, k)
    b = sp.bits_per_res
    S = W.enumerate_states(L, k)
    _, tot = W.legacy_tables(sp, S)
    T = np.asarray(tot, float).reshape((k,) * L)
    rng_ = float(np.percentile(T, 99) - np.percentile(T, 1))
    out = {"pdb": pdb_id, "L": L, "k": k,
           "m_binary": L * b, "m_onehot": L * k,
           "binary": {kk: v for kk, v in
                      W.spectrum(T.ravel(), L * b, bits_per_res=b).items()
                      if kk in ("W_mean", "V1", "V2", "tail_gt3", "D_mean", "var")}}
    for tag, kwargs in [("onehot_mean", dict(fill="mean")),
                        ("onehot_pen_1x", dict(fill="penalty", penalty=rng_)),
                        ("onehot_pen_10x", dict(fill="penalty", penalty=10 * rng_))]:
        E = W.onehot_table(T, k, **kwargs)
        s = W.spectrum(E, L * k)
        out[tag] = {kk: s[kk] for kk in ("W_mean", "W_eff", "V1", "V2", "V3", "tail_gt3",
                                         "var", "parseval_rel_err")}
        out[tag]["V"] = s["V"]
    print(f"  ONEHOT {pdb_id} L={L} k={k}: binary W {out['binary']['W_mean']:.3f} "
          f"(m={L*b})  onehot-mean W {out['onehot_mean']['W_mean']:.3f} "
          f"pen1x {out['onehot_pen_1x']['W_mean']:.3f} "
          f"pen10x {out['onehot_pen_10x']['W_mean']:.3f} (m={L*k})", flush=True)
    return out


def main(argv):
    rng = np.random.default_rng(7)
    out = {"rows": [], "onehot": []}
    for (L, k) in CELLS:
        for p in TARGETS:
            Q.wait_for_memory(0.7, tag="walsh_encode")
            out["rows"].append(cell(p, L, k, rng))
            W.write("walsh_encode", out)
    for (L, k) in ONEHOT_CELLS:
        for p in TARGETS[:3]:
            Q.wait_for_memory(0.7, tag="walsh_encode")
            out["onehot"].append(onehot_cell(p, L, k))
            W.write("walsh_encode", out)
    print("wrote", W.write("walsh_encode", out))


if __name__ == "__main__":
    main(sys.argv[1:])
