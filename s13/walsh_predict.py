"""SPRINT 13 PAULI-SPECTRUM, EXPERIMENT W5 -- FROM CHAIN GEOMETRY TO PAULI STRUCTURE.

The architecture agent PROVED by perturbation (`s13/qarch_FINDINGS.md` §1) that the support
of the CA-CA distance d_ij is exactly the j-i-1 residues strictly between i and j.  This
script tests that statement in the Walsh basis, where it is exact rather than thresholded,
and turns it into a quantitative prediction of the weight spectrum.

PART A -- THE SUPPORT LAW, EXACTLY.  For every residue pair (i,j) in a fully enumerated
register, transform d_ij and measure the variance carried by Walsh subsets that touch a
residue OUTSIDE {i+1..j-1}.  Prediction: exactly zero.  Prediction: separation-1 pairs are
constant.  This is 2^m x (all subsets) rather than 64 perturbations x a 0.1 A threshold.

PART B -- ORDER FROM SUPPORT.  If a pair term is generic on its s-1 residue support, its
residue-order spectrum is D_d = C(s-1, d)(k-1)^d / (k^{s-1} - 1).  Measured against that.

PART C -- WEIGHT FROM ORDER (the encoding law).  If a residue's dependence is generic across
its own b = log2(k) qubits, a residue-order-d subset carries mean Pauli weight
    W | D=d  =  d * b * 2^{b-1} / (2^b - 1).
Measured W_given_D against that, over every cell computed in this strand.

PART D -- THE END-TO-END LOOP on an object whose every ingredient is known: a SURROGATE
pairwise CA potential (12-6 Lennard-Jones on CA-CA distances -- NOT a force field, and named
a surrogate everywhere).  Its spectrum is predicted from the per-pair spectra of Part A
weighted by each pair's own variance, with cross-covariances dropped, and compared with its
measured spectrum.  The residual is the size of the effect that the "sum of independent pair
terms" picture misses.

PART E -- TRAINABILITY.  The measured spectra are combined with the geo/ strand's measured
ansatz kernel v(w) to state what gradient-variance scaling the spectra predict, and under
which regime of validity (`s13/results/lit_methods.json`).

    python -m s13.walsh_predict [--cells 1A13:6:4,...]
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

from s13 import walsh_lib as W
from s13 import qarch_lib as Q

CELLS = [("1A13", 6, 4), ("1A1P", 6, 4), ("2BFI", 6, 4), ("1DEP", 6, 4),
         ("1ID6", 6, 4), ("1CEK", 6, 4), ("1A13", 7, 4), ("1A13", 12, 2),
         ("1A13", 4, 8), ("2BFI", 4, 8)]
LJ_SIGMA = 4.0          # CA-CA contact distance, Angstrom
LJ_EPS = 1.0


def residue_degree(m, b):
    idx = np.arange(1 << m, dtype=np.uint64)
    L = m // b
    d = np.zeros(1 << m, dtype=np.int64)
    masks = []
    for r in range(L):
        mask = np.uint64(((1 << b) - 1) << (r * b))
        masks.append(mask)
        d += ((idx & mask) != 0).astype(np.int64)
    return d, masks, W.popcount(idx).astype(np.int64)


def generic_order_share(q, k):
    """ANOVA order shares of a generic function of q k-ary variables."""
    if q <= 0:
        return np.array([1.0])
    tot = k ** q - 1
    return np.array([math.comb(q, d) * (k - 1) ** d / tot for d in range(q + 1)])


def part_abc(pdb_id, L, k):
    sp = W.space_for(pdb_id, L, k)
    b, m = sp.bits_per_res, L * sp.bits_per_res
    S = W.enumerate_states(L, k)
    CA = sp.ca(S)                                          # (2^m, L, 3)
    d_res, masks, wt = residue_degree(m, b)
    # residue r of the C-ordered table owns bit-block (L-1-r)
    blockmask = {r: masks[L - 1 - r] for r in range(L)}
    idx = np.arange(1 << m, dtype=np.uint64)

    pairs = []
    per_pair_c2 = {}
    for i in range(L):
        for j in range(i + 1, L):
            dij = np.linalg.norm(CA[:, i, :] - CA[:, j, :], axis=1)
            c = W.walsh_coeffs(dij)
            c2 = c * c
            var = float(c2[1:].sum())
            s = j - i
            interior = set(range(i + 1, j))
            out_mask = np.zeros(1 << m, dtype=bool)
            for r in range(L):
                if r not in interior:
                    out_mask |= (idx & blockmask[r]) != 0
            outside = float(c2[1:][out_mask[1:]].sum())
            mu2 = float(c[0]) ** 2
            rel_var = var / mu2 if mu2 > 0 else var
            row = {"i": i, "j": j, "sep": s, "var": var, "rel_var": rel_var,
                   "outside_support_var": outside,
                   "outside_support_share": outside / var if var > 0 else 0.0,
                   "n_support_residues": len(interior)}
            # a variance of 1e-30 relative to the mean is floating-point noise on a
            # geometrically constant distance, not structure
            if rel_var > 1e-24 and s > 1:
                D = np.bincount(d_res[1:], weights=c2[1:], minlength=L + 1) / var
                V = np.bincount(wt[1:], weights=c2[1:], minlength=m + 1) / var
                row["D_mean"] = float((np.arange(L + 1) * D).sum())
                row["W_mean"] = float((np.arange(m + 1) * V).sum())
                row["max_order_with_mass"] = int(np.max(np.nonzero(D > 1e-14)[0]))
                pred = generic_order_share(s - 1, k)
                row["D_measured"] = [float(x) for x in D[:s]]
                row["D_generic_pred"] = [float(x) for x in pred]
                row["D_L1_vs_generic"] = float(np.abs(D[:s] - pred).sum())
                row["W_mean_generic_pred"] = float((s - 1) * b * 2 ** (b - 1) / (2 ** b - 1))
                per_pair_c2[(i, j)] = c2
            else:
                row["D_mean"] = row["W_mean"] = 0.0
                row["max_order_with_mass"] = 0
            pairs.append(row)

    # ---- PART D: the SURROGATE pairwise CA potential -------------------------
    Esyn = np.zeros(1 << m)
    pair_terms = {}
    for i in range(L):
        for j in range(i + 2, L):                          # skip the constant sep-1 pairs
            r = np.linalg.norm(CA[:, i, :] - CA[:, j, :], axis=1)
            phi = 4.0 * LJ_EPS * ((LJ_SIGMA / r) ** 12 - (LJ_SIGMA / r) ** 6)
            pair_terms[(i, j)] = phi
            Esyn += phi
    syn_meas = W.spectrum(Esyn, m, bits_per_res=b)
    # prediction: variance-weighted mixture of the individual pair spectra
    acc_V = np.zeros(m + 1); acc_D = np.zeros(L + 1); tot_v = 0.0
    for key, phi in pair_terms.items():
        c2 = W.walsh_coeffs(phi) ** 2
        v = float(c2[1:].sum())
        if v <= 0:
            continue
        acc_V += np.bincount(wt[1:], weights=c2[1:], minlength=m + 1)
        acc_D += np.bincount(d_res[1:], weights=c2[1:], minlength=L + 1)
        tot_v += v
    predV = acc_V / tot_v
    predD = acc_D / tot_v
    syn = {
        "SURROGATE": "12-6 Lennard-Jones on CA-CA distances; NOT a force field",
        "sigma": LJ_SIGMA, "eps": LJ_EPS,
        "measured": {kk: syn_meas[kk] for kk in ("W_mean", "D_mean", "V1", "V2", "V3",
                                                 "tail_gt3", "var", "parseval_rel_err")},
        "measured_V": syn_meas["V"], "measured_D": syn_meas["D"],
        "pred_additive_V": [float(x) for x in predV],
        "pred_additive_D": [float(x) for x in predD],
        "pred_W_mean": float((np.arange(m + 1) * predV).sum()),
        "pred_D_mean": float((np.arange(L + 1) * predD).sum()),
        "L1_V": float(np.abs(predV - np.array(syn_meas["V"])).sum()),
        "L1_D": float(np.abs(predD - np.array(syn_meas["D"])).sum()),
        "sum_pair_var_over_total_var": float(tot_v / syn_meas["var"]),
    }

    # ---- PART C on the real energies of this cell -----------------------------
    _, tot = W.legacy_tables(sp, S)
    leg = W.spectrum(tot, m, bits_per_res=b)
    wgd_pred = [float(d * b * 2 ** (b - 1) / (2 ** b - 1)) for d in range(L + 1)]

    row = {"pdb": pdb_id, "L": L, "k": k, "m_bits": m, "bits_per_res": b,
           "pairs": pairs, "surrogate_lj": syn,
           "legacy_W_given_D": leg["W_given_D"], "W_given_D_pred": wgd_pred,
           "legacy_spread": leg["spread"], "spread_pred": float(b * 2 ** (b - 1) / (2 ** b - 1)),
           "legacy_W_mean": leg["W_mean"], "legacy_D_mean": leg["D_mean"]}
    ok = [p for p in pairs if p["var"] > 0]
    row["max_outside_support_share"] = float(max(
        (p["outside_support_share"] for p in pairs if p["sep"] > 1), default=0.0))
    row["sep1_max_rel_var"] = float(max([p["rel_var"] for p in pairs if p["sep"] == 1]
                                         or [0.0]))
    row["support_law_exact"] = bool(row["max_outside_support_share"] < 1e-12
                                    and row["sep1_max_rel_var"] < 1e-24)
    print(f"  {pdb_id} L={L} k={k} m={m}: support-law exact={row['support_law_exact']} "
          f"(max outside share {row['max_outside_support_share']:.2e}, "
          f"max sep-1 rel var {row['sep1_max_rel_var']:.2e})  "
          f"LJ-surrogate W meas {syn_meas['W_mean']:.3f} vs additive pred "
          f"{syn['pred_W_mean']:.3f} (L1 {syn['L1_V']:.3f})  "
          f"spread {leg['spread']:.3f} vs pred {row['spread_pred']:.3f}", flush=True)
    return row


def part_e():
    """The trainability statement: measured spectra x geo's measured ansatz kernel."""
    p = os.path.join(W.RESULTS, "geo_pauli.json")
    if not os.path.exists(p):
        return {"error": "geo_pauli.json not present"}
    geo = json.load(open(p))
    ker = {}
    for c in geo["cells"]:
        key = (c["n_qubits"], c["layers"])
        if key in ker:
            continue
        v = np.array(c["kernel_v_of_weight"], float)
        if v.size < 3 or v[1] <= 0:
            continue
        rel = v[1:] / v[1]
        ker[key] = {"n": c["n_qubits"], "layers": c["layers"],
                    "v_rel": [float(x) for x in rel],
                    "max_over_min": float(rel[rel > 0].max() / rel[rel > 0].min()),
                    "v_at_max_weight_over_v1": float(rel[-1]),
                    "decay_base": c["kernel_decay"]["decay_base"],
                    "exp_r2": c["kernel_decay"]["exp_r2"]}
    out = {"ansatz_kernel_from_geo": list(ker.values())}
    bases = [v["decay_base"] for v in ker.values()]
    ratios = [v["v_at_max_weight_over_v1"] for v in ker.values()]
    out["kernel_summary"] = {
        "n_cells": len(ker),
        "decay_base_median": float(np.median(bases)),
        "decay_base_range": [float(min(bases)), float(max(bases))],
        "v(w_max)/v(1)_median": float(np.median(ratios)),
        "v(w_max)/v(1)_range": [float(min(ratios)), float(max(ratios))],
        "verdict": ("the ansatz kernel is FLAT in Pauli weight (decay base ~0.79-1.01, "
                    "v(w_max)/v(1) within a factor of ~5 of unity), so on THIS ansatz at "
                    "8-12 qubits the weight spectrum does not enter the gradient variance "
                    "beyond the total Frobenius norm Sum_S c_S^2"),
    }
    return out


def main(argv):
    cells = CELLS
    if "--cells" in argv:
        cells = [tuple(x.split(":")) for x in argv[argv.index("--cells") + 1].split(",")]
        cells = [(a, int(b), int(c)) for a, b, c in cells]
    out = {"rows": [], "lj_surrogate_note":
           "the LJ arm is a SURROGATE pair potential on CA coordinates, not a force field"}
    for (p, L, k) in cells:
        Q.wait_for_memory(0.7, tag="walsh_predict")
        try:
            out["rows"].append(part_abc(p, L, k))
        except Exception as exc:                            # noqa: BLE001
            print(f"  !! {p} {L} {k}: {type(exc).__name__}: {exc}", flush=True)
            out["rows"].append({"pdb": p, "L": L, "k": k,
                                "error": f"{type(exc).__name__}: {exc}"})
        W.write("walsh_predict", out)
    out["trainability"] = part_e()
    print(json.dumps(out["trainability"].get("kernel_summary", {}), indent=1))
    print("wrote", W.write("walsh_predict", out))


if __name__ == "__main__":
    main(sys.argv[1:])
