"""SPRINT 13 PAULI-SPECTRUM, EXPERIMENT W0 -- VERIFY THE TRANSFORM BEFORE TRUSTING IT.

Nothing downstream is worth reading unless these pass:

 1. the FWHT agrees with the naive O(4^m) Walsh matrix on random data;
 2. Parseval holds and the inverse transform rebuilds a REAL Legacy energy table to
    machine precision;
 3. the 1-local torsion prior -- which is 1-local BY CONSTRUCTION in the residue variables
    -- comes out with residue-order share D_1 = 1 exactly and Pauli weight supported on
    w <= log2(k).  This is the correctness check the brief asks for, and it already shows
    the encoding point: a "1-local" prior is NOT weight-1 in qubits;
 4. the unbiased Krawtchouk pair estimator reproduces the exact spectrum;
 5. the order/weight bound d <= |S| <= d*b holds on every Walsh subset.

    python -m s13.walsh_verify
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from s13 import walsh_lib as W
from s13 import qarch_lib as Q

OUT = {}


def naive_walsh(a):
    m = int(np.log2(len(a)))
    idx = np.arange(len(a), dtype=np.uint64)
    H = (-1.0) ** (W.popcount(idx[:, None] & idx[None, :]) % 2)
    return H @ a / len(a)


def main():
    t0 = time.time()
    rng = np.random.default_rng(0)

    # ---- 1. transform correctness against the naive matrix ---------------------
    a = rng.normal(size=1 << 10)
    c_fast = W.walsh_coeffs(a)
    c_slow = naive_walsh(a)
    OUT["fwht_vs_naive_max_abs"] = float(np.abs(c_fast - c_slow).max())
    OUT["fwht_vs_naive_rel"] = float(np.abs(c_fast - c_slow).max() / np.abs(c_slow).max())

    # ---- 5. order/weight bound -------------------------------------------------
    OUT["order_weight_bound"] = {f"m={m},b={b}": W.check_order_weight_bound(m, b)
                                 for m, b in [(8, 2), (12, 2), (12, 3), (10, 1), (16, 4)]}

    # ---- 2,3. real Legacy table + the 1-local prior ----------------------------
    checks = []
    for pdb_id, L, k in [("1A13", 5, 4), ("1A13", 8, 2), ("2BFI", 5, 4), ("1A1P", 4, 8)]:
        sp = W.space_for(pdb_id, L, k)
        b = sp.bits_per_res
        m = L * b
        S = W.enumerate_states(L, k)
        comp, tot = W.legacy_tables(sp, S)
        spec = W.spectrum(tot, m, bits_per_res=b)
        P = Q.empirical_prior(sp)
        pri = Q.prior_energy(P, S)
        pspec = W.spectrum(pri, m, bits_per_res=b)
        row = {
            "pdb": pdb_id, "L": L, "k": k, "m_bits": m,
            "legacy_parseval_rel_err": spec["parseval_rel_err"],
            "legacy_reconstruct_max_rel_err": spec["reconstruct_max_rel_err"],
            "legacy_var": spec["var"], "legacy_W_mean": spec["W_mean"],
            "legacy_D_mean": spec["D_mean"],
            "prior_parseval_rel_err": pspec["parseval_rel_err"],
            "prior_reconstruct_max_rel_err": pspec["reconstruct_max_rel_err"],
            "prior_D1": pspec["D1"],
            "prior_D1_err": abs(pspec["D1"] - 1.0),
            "prior_V": pspec["V"],
            "prior_mass_above_b": float(sum(pspec["V"][b + 1:])),
            "prior_W_mean": pspec["W_mean"],
            "order_weight_bound_holds": spec["order_weight_bound_holds"],
        }
        checks.append(row)
        print(f"  {pdb_id} L={L} k={k} m={m}: parseval {spec['parseval_rel_err']:.2e} "
              f"recon {spec['reconstruct_max_rel_err']:.2e} priorD1 {pspec['D1']:.15f} "
              f"prior mass>b {row['prior_mass_above_b']:.2e} priorW {pspec['W_mean']:.4f}",
              flush=True)
    OUT["exact_checks"] = checks

    # ---- 4. sampling estimator against the exact spectrum ---------------------
    sp = W.space_for("1A13", 6, 4)
    m = 12
    S = W.enumerate_states(6, 4)
    _, tot = W.legacy_tables(sp, S)
    exact = W.spectrum(tot, m, bits_per_res=2)
    samp = {}
    for N in (1024, 4096, 16384):
        j = rng.integers(0, 1 << m, size=N)
        est = W.sampled_spectrum(j.astype(np.uint64), tot[j], m, blocks=8, seed=1)
        samp[str(N)] = {
            "V1_exact": exact["V1"], "V1_est": est["V1"], "V1_sd": est["V_block_sd"][1],
            "V2_exact": exact["V2"], "V2_est": est["V2"], "V2_sd": est["V_block_sd"][2],
            "V3_exact": exact["V3"], "V3_est": est["V3"], "V3_sd": est["V_block_sd"][3],
            "W_exact": exact["W_mean"], "W_est": est["W_mean"], "W_sd": est["W_mean_sd"],
            "var_exact": exact["var"], "var_est": est["var_est"],
            "L1_over_weights": float(np.abs(np.array(est["V_norm"])
                                            - np.array(exact["V"])).sum()),
        }
        print(f"  sampled N={N}: W {est['W_mean']:.3f}+-{est['W_mean_sd']:.3f} "
              f"(exact {exact['W_mean']:.3f})  L1 {samp[str(N)]['L1_over_weights']:.4f}",
              flush=True)
    OUT["sampling_check"] = samp
    OUT["seconds"] = time.time() - t0
    p = W.write("walsh_verify", OUT)
    print("wrote", p, f"in {OUT['seconds']:.1f}s")


if __name__ == "__main__":
    main()
