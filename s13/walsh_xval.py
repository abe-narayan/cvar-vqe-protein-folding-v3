"""SPRINT 13 PAULI-SPECTRUM -- INDEPENDENT CROSS-VALIDATION AGAINST THE geo/ STRAND.

Two agents computed Pauli-weight spectra of the same energy models by separate code paths.
Nobody has measured this quantity for a molecular force field before, so there is no
literature to check it against; two independent implementations agreeing is the only
available credibility check.  This script performs it at two levels.

LEVEL 1 -- THE TRANSFORM.  Take geo's OWN cached energy tables (`s13/cache/geo_E_*.npy`),
run MY transform (`walsh_lib.spectrum`) over them, and compare against the numbers in
`s13/results/geo_pauli.json` coefficient by coefficient.  This isolates the transform from
the energy table.

LEVEL 2 -- THE ENERGY TABLE.  Rebuild the same cells from scratch through my own code path
under geo's library convention and compare the tables elementwise.

A NOTE ON THE CONVENTION, because it explains a discrepancy that would otherwise look like a
bug: geo builds the library as `library_for(seq[:L], k, exclude_seq=FULL target sequence)`;
my own sweeps use `library_for(seq[:L], k, exclude_seq=seq[:L])`.  Both are leakage-safe,
but they hold out different pools and therefore give DIFFERENT libraries and different
energies.  Cross-validation is run under geo's convention.

    python -m s13.walsh_xval
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys

import numpy as np

from s13 import walsh_lib as W
from s13 import qarch_lib as Q
import torsion_lib2 as tl2


def geo_space(pdb_id, L, k):
    """geo's convention: prefix sequence, FULL target sequence held out."""
    t = {x["pdb"]: x for x in Q.I.targets()}[pdb_id]
    full = t["seq"]
    sub = full[:L]
    tab = tl2.library_for(sub, k, full)
    sp = Q.Space.__new__(Q.Space)
    sp.pdb, sp.seq, sp.n, sp.fold, sp.k = pdb_id, sub, L, int(t["fold"]), k
    sp.rep = tl2.PerResidueTorsion(sub, tab, chi_bits=False)
    sp.PHI = np.ascontiguousarray(sp.rep._phi, float)
    sp.PSI = np.ascontiguousarray(sp.rep._psi, float)
    sp.bits_per_res = int(np.log2(k))
    sp.n_qubits_binary = L * sp.bits_per_res
    sp._rows = np.arange(L)
    return sp


def main():
    geo = json.load(open(os.path.join(W.RESULTS, "geo_pauli.json")))
    ref = {}
    for c in geo["cells"]:
        for model in ("legacy", "amber", "amber_soft"):
            if model in c and isinstance(c[model], dict) and "spectrum" in c[model]:
                ref.setdefault((c["pdb"], c["L"], c["k"], model), c[model]["spectrum"])

    out = {"level1_transform": [], "level2_energy_table": [], "convention_note": (
        "geo: library_for(seq[:L], k, exclude_seq=FULL seq); "
        "walsh sweeps: library_for(seq[:L], k, exclude_seq=seq[:L]). Different libraries.")}

    # ---------------- LEVEL 1: my transform on geo's tables --------------------
    for path in sorted(glob.glob(os.path.join(W.CACHE, "geo_E_*.npy"))):
        mm = re.match(r"geo_E_(\w+)_(\d+)_(\d+)_(\w+)\.npy", os.path.basename(path))
        if not mm:
            continue
        pdb_id, L, k, model = mm.group(1), int(mm.group(2)), int(mm.group(3)), mm.group(4)
        E = np.load(path)
        m = int(round(np.log2(E.size)))
        if E.size != 1 << m:
            continue
        s = W.spectrum(E, m, bits_per_res=int(np.log2(k)))
        r = ref.get((pdb_id, L, k, model))
        row = {"pdb": pdb_id, "L": L, "k": k, "model": model, "m_bits": m,
               "mine_mean_weight": s["W_mean"], "mine_var": s["var"],
               "mine_parseval_rel_err": s["parseval_rel_err"]}
        if r:
            row["geo_mean_weight"] = r["mean_weight"]
            row["geo_var"] = r["parseval_var"]
            row["abs_diff_mean_weight"] = abs(s["W_mean"] - r["mean_weight"])
            row["rel_diff_var"] = abs(s["var"] - r["parseval_var"]) / max(abs(r["parseval_var"]), 1e-300)
            gv = np.array(r["weight_share"])
            mv = np.array(s["V"][:len(gv)])
            row["max_abs_diff_weight_share"] = float(np.abs(gv - mv).max())
            row["L1_weight_share"] = float(np.abs(gv - mv).sum())
        out["level1_transform"].append(row)
        if r:
            print(f"  L1 {pdb_id} L={L} k={k} {model:10s} m={m}: "
                  f"dW {row['abs_diff_mean_weight']:.3e}  dVar/Var {row['rel_diff_var']:.3e} "
                  f" maxdV_w {row['max_abs_diff_weight_share']:.3e}", flush=True)

    # ---------------- LEVEL 2: rebuild the energy table independently ----------
    for pdb_id, L, k in [("1A13", 5, 4), ("2BFI", 5, 4), ("1A1P", 4, 8), ("1A13", 10, 2)]:
        p = os.path.join(W.CACHE, f"geo_E_{pdb_id}_{L}_{k}_legacy.npy")
        if not os.path.exists(p):
            continue
        sp = geo_space(pdb_id, L, k)
        S = W.enumerate_states(L, k)
        _, mine = W.legacy_tables(sp, S)
        theirs = np.load(p)
        d = np.abs(np.asarray(mine, float) - theirs)
        row = {"pdb": pdb_id, "L": L, "k": k, "model": "legacy",
               "max_abs_diff": float(d.max()),
               "max_rel_diff": float(d.max() / max(np.abs(theirs).max(), 1e-300)),
               "n": int(len(theirs))}
        out["level2_energy_table"].append(row)
        print(f"  L2 {pdb_id} L={L} k={k} legacy: max abs diff {row['max_abs_diff']:.3e} "
              f"(rel {row['max_rel_diff']:.3e}) over {row['n']} configurations", flush=True)

    # ---------------- the concentration diagnostic -----------------------------
    conc = []
    for path in sorted(glob.glob(os.path.join(W.CACHE, "geo_E_*.npy"))
                       + glob.glob(os.path.join(W.CACHE, "walsh_E_*.npy"))):
        base = os.path.basename(path)
        mm = re.match(r"\w+_E_(\w+)_(\d+)_(\d+)_(\w+)\.npy", base)
        if not mm:
            continue
        pdb_id, L, k, model = mm.group(1), int(mm.group(2)), int(mm.group(3)), mm.group(4)
        E = np.load(path)
        m = int(round(np.log2(E.size)))
        if E.size != 1 << m or m < 6:
            continue
        conc.append(dict(concentration(E, m), pdb=pdb_id, L=L, k=k, model=model,
                         m_bits=m, source=base.split("_")[0]))
    out["concentration"] = conc
    print("\n  CONCENTRATION (share of the Walsh variance carried by the single most "
          "extreme configuration, and L1 distance of the weight spectrum from the "
          "pure-delta binomial):")
    for r in sorted(conc, key=lambda r: (r["model"], r["m_bits"])):
        print(f"    {r['model']:11s} {r['pdb']} L={r['L']} k={r['k']} m={r['m_bits']:2d}: "
              f"top1 {r['var_share_top1']:.4f}  top10 {r['var_share_top10']:.4f}  "
              f"W {r['W_mean']:.3f} (delta ref {r['binomial_W_mean']:.3f})  "
              f"L1_from_delta {r['L1_from_binomial']:.4f}", flush=True)
    print("wrote", W.write("walsh_xval", out))


def concentration(E, m):
    """How much of the Walsh variance is carried by the few most extreme configurations,
    and how close the weight spectrum is to that of a single delta spike (binomial)."""
    from math import comb
    E = np.asarray(E, float)
    n = 1 << m
    s = W.spectrum(E, m)
    dev = (E - E.mean()) ** 2 / n
    o = np.argsort(-dev)
    binom = np.array([comb(m, w) for w in range(m + 1)], float)
    binom /= (2 ** m - 1)
    V = np.array(s["V"])
    return {"var": s["var"], "W_mean": s["W_mean"],
            "binomial_W_mean": float((np.arange(m + 1) * binom).sum()),
            "L1_from_binomial": float(np.abs(V - binom).sum()),
            "var_share_top1": float(dev[o[0]] / s["var"]),
            "var_share_top10": float(dev[o[:10]].sum() / s["var"]),
            "var_share_top1pct": float(dev[o[:max(1, n // 100)]].sum() / s["var"]),
            "E_min": float(E.min()), "E_max": float(E.max())}


if __name__ == "__main__":
    main()
