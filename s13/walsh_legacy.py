"""SPRINT 13 PAULI-SPECTRUM, EXPERIMENT W1 -- THE EXACT SPECTRUM OF THE LEGACY POTENTIAL.

Exact enumeration of the whole register, exact Walsh-Hadamard transform, exact Pauli-weight
spectrum -- for the 11-term Legacy total, for EVERY term separately, and for the 1-local
torsion prior, over three crossing ladders so that qubit count, peptide length and
states-per-residue are not confounded:

    k=4 : L = 4..9    -> m =  8,10,12,14,16,18     (the sprint's working encoding)
    k=2 : L = 8..16   -> m =  8..16                (twice the peptide at equal m)
    k=8 : L = 4,5     -> m = 12,15                 (same peptide, more states)

Targets: the six s12 tuning targets geo/ uses, one per leave-fold-out fold, two of them
FAIL18 members.  benchmark60 and dev24 are untouched.

    python -m s13.walsh_legacy [--ladder k4|k2|k8] [--targets A,B] [--maxm 18]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

from s13 import walsh_lib as W
from s13 import qarch_lib as Q
from core import energy as et

TARGETS = ["1A13", "1A1P", "2BFI", "1DEP", "1ID6", "1CEK"]
LADDERS = {
    "k4": [(L, 4) for L in (4, 5, 6, 7, 8, 9)],
    "k2": [(L, 2) for L in (8, 10, 12, 14, 16)],
    "k8": [(L, 8) for L in (4, 5)],
}


def cell(pdb_id, L, k, save_table=True):
    t0 = time.time()
    sp = W.space_for(pdb_id, L, k)
    b = sp.bits_per_res
    m = L * b
    S = W.enumerate_states(L, k)
    comp, tot = W.legacy_tables(sp, S, chunk=4096)
    dt_energy = time.time() - t0

    spec = W.spectrum(tot, m, bits_per_res=b, want_top=8)
    row = {"pdb": pdb_id, "L": L, "k": k, "m_bits": m, "bits_per_res": b,
           "n_configs": int(len(S)), "seconds_energy": dt_energy,
           "ms_per_config": 1000.0 * dt_energy / len(S),
           "total": spec}

    # ---- per-term spectra, and each term's weighted share of the total variance
    wts = et.DEFAULT_WEIGHTS
    vt = spec["var"]
    terms = {}
    for name in Q.LEGACY_TERMS:
        e = np.asarray(comp[name], float)
        s = W.spectrum(e, m, bits_per_res=b)
        wt = float(wts.get(name, 0.0))
        s["weight"] = wt
        s["weighted_var"] = wt * wt * s["var"]
        s["weighted_var_share"] = s["weighted_var"] / vt if vt > 0 else float("nan")
        # covariance-aware share: contribution of this term to Var(total)
        s["cov_share"] = float(wt * np.cov(e, np.asarray(tot, float))[0, 1] / vt) if vt > 0 \
            else float("nan")
        for key in ("top_terms",):
            s.pop(key, None)
        terms[name] = s
    row["terms"] = terms

    # ---- the 1-local prior on the same register (correctness anchor)
    P = Q.empirical_prior(sp)
    row["prior"] = W.spectrum(Q.prior_energy(P, S), m, bits_per_res=b)

    if save_table and m >= 12:
        np.save(os.path.join(W.CACHE, f"walsh_E_{pdb_id}_{L}_{k}_legacy.npy"),
                np.asarray(tot, np.float64))
    print(f"  {pdb_id} L={L} k={k} m={m}: {dt_energy:.1f}s "
          f"({row['ms_per_config']:.3f} ms/cfg)  W_mean {spec['W_mean']:.3f}  "
          f"D_mean {spec['D_mean']:.3f}  V1 {spec['V1']:.4f} tail>3 {spec['tail_gt3']:.4f}  "
          f"parseval {spec['parseval_rel_err']:.1e}", flush=True)
    return row


def main(argv):
    ladders = ["k4", "k2", "k8"]
    targets = TARGETS
    maxm = 18
    i = 0
    while i < len(argv):
        if argv[i] == "--ladder":
            ladders = argv[i + 1].split(","); i += 2
        elif argv[i] == "--targets":
            targets = argv[i + 1].split(","); i += 2
        elif argv[i] == "--maxm":
            maxm = int(argv[i + 1]); i += 2
        else:
            i += 1
    out = {"targets": targets, "ladders": ladders, "rows": []}
    fn = f"walsh_legacy_{'_'.join(ladders)}_m{maxm}_{len(targets)}t"
    for lad in ladders:
        for (L, k) in LADDERS[lad]:
            if L * int(np.log2(k)) > maxm:
                continue
            for pdb_id in targets:
                # this process's own peak RSS is ~150 MB at m=18 (measured), far below the
                # 1.2 GB cap; the wait threshold is set accordingly and logged.
                Q.wait_for_memory(0.7, tag="walsh_legacy")
                try:
                    out["rows"].append(dict(cell(pdb_id, L, k), ladder=lad))
                except Exception as exc:                       # noqa: BLE001
                    print(f"  !! {pdb_id} L={L} k={k}: {type(exc).__name__}: {exc}",
                          flush=True)
                    out["rows"].append({"pdb": pdb_id, "L": L, "k": k, "ladder": lad,
                                        "error": f"{type(exc).__name__}: {exc}"})
                W.write(fn, out)
    print("wrote", W.write(fn, out))


if __name__ == "__main__":
    main(sys.argv[1:])
