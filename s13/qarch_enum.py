"""SPRINT 13 QUANTUM-ARCHITECTURE, EXPERIMENT 1 -- small-instance ground truth.

The nine n=9 tuning targets at k=4 give 4^9 = 262,144 configurations: the whole discrete
torsion state space, enumerable.  For every configuration we cache

    rmsd      ORACLE post-hoc CA-RMSD of the ideal-geometry CA trace to the native
    legacy    all 11 `core.energy` terms (and their DEFAULT_WEIGHTS total)
    prior     the 1-local torsion-prior energy under the leakage-safe empirical prior

and, on a labelled 3,000-configuration subsample (uniform / prior-sampled / oracle
near-native band), the genuine AMBER ff14SB/GBn2 five-term single point.

That gives certified optima for every candidate Hamiltonian on a real instance family,
which is what the trainability agent and the classical controls need to score against.

Output: `s13/results/qarch_enum_<PDB>.npz` (the arrays) + `qarch_enum.json` (the summary).
The npz files are the cache; nothing downstream re-enumerates.

    python -m s13.qarch_enum [PDB ...]
"""
from __future__ import annotations

import itertools
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import qarch_lib as Q          # noqa: E402

K = 4
N_AMBER_UNIFORM = 1200
N_AMBER_PRIOR = 1200
N_AMBER_BAND = 600
SMALL = ["1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5"]


def all_configs(n, k):
    """(k^n, n) int8 enumeration in odometer order; index = sum s_i k^(n-1-i)."""
    return np.array(list(itertools.product(range(k), repeat=n)), dtype=np.int8)


def enumerate_target(pdb_id, seed=0, amber=True, skip_existing=True):
    t0 = time.time()
    done = os.path.join(Q.RESULTS, f"qarch_enum_{pdb_id}.npz")
    if skip_existing and os.path.exists(done):
        z = np.load(done)
        if "amber_total" in z.files and np.isfinite(np.asarray(z["amber_total"])).any():
            print(f"{pdb_id}: cached, skipping", flush=True)
            r = np.asarray(z["rmsd"], float); a = np.asarray(z["amber_total"], float)
            ai = np.asarray(z["amber_idx"], int)
            return dict(pdb=pdb_id, n=int(z["n"]), k=int(z["k"]),
                        qubits=int(int(z["n"]) * np.log2(int(z["k"]))),
                        configs=int(len(r)), rmsd_min=float(r.min()),
                        rmsd_max=float(r.max()), rmsd_mean=float(r.mean()),
                        rmsd_median=float(np.median(r)),
                        snap_rmsd=float(r[int(z["snap_index"])]),
                        legacy_argmin_rmsd=float(r[int(np.argmin(z["legacy"]))]),
                        prior_argmin_rmsd=float(r[int(np.argmin(z["prior"]))]),
                        amber_n=int(len(ai)),
                        amber_argmin_rmsd=float(r[ai[int(np.nanargmin(a))]]),
                        wall_s=0.0, path=done, cached=True)
    sp = Q.Space(pdb_id, K)
    n = sp.n
    S = all_configs(n, K)
    B = len(S)
    print(f"{pdb_id}: n={n} k={K} configs={B} qubits={sp.n_qubits_binary}", flush=True)

    rmsd = sp.rmsd(S).astype(np.float32)
    print(f"  rmsd done [{time.time()-t0:.0f}s]  min={rmsd.min():.3f}", flush=True)

    comp = Q.legacy_components(sp, S, chunk=8192)
    legacy = Q.legacy_total(comp).astype(np.float32)
    print(f"  legacy done [{time.time()-t0:.0f}s]", flush=True)

    Pemp = Q.empirical_prior(sp)
    pri = Q.prior_energy(Pemp, S).astype(np.float32)
    snap = sp.ORACLE_snap()
    snap_idx = int(np.ravel_multi_index(tuple(snap), (K,) * n))

    # ---- AMBER subsample (labelled populations) --------------------------
    rng = np.random.default_rng(seed)
    idx_u = rng.choice(B, N_AMBER_UNIFORM, replace=False)
    Sp = sp.sample_prior(Pemp, N_AMBER_PRIOR, rng)
    idx_p = np.ravel_multi_index(tuple(Sp.T), (K,) * n)
    band = np.argsort(rmsd)[:max(2000, B // 100)]
    idx_b = rng.choice(band, N_AMBER_BAND, replace=False)
    a_idx = np.unique(np.concatenate([idx_u, idx_p, idx_b, [snap_idx]]))
    a_kind = np.zeros(len(a_idx), np.int8)
    a_kind[np.isin(a_idx, idx_p)] = 1
    a_kind[np.isin(a_idx, idx_b)] = 2

    amb_tot = np.full(len(a_idx), np.nan, np.float64)
    amb_comp = {}
    if amber:
        Q.wait_for_memory(1.5, pdb_id)
        ta = time.time()
        c, tt = Q.amber_energies(sp, S[a_idx], components=True, progress=400)
        amb_tot = tt
        amb_comp = {k: v.astype(np.float32) for k, v in c.items()}
        print(f"  amber {len(a_idx)} in {time.time()-ta:.0f}s "
              f"({(time.time()-ta)/len(a_idx)*1000:.1f} ms/cfg)", flush=True)

    out = dict(pdb=pdb_id, n=n, k=K, seq=sp.seq, fold=sp.fold,
               rmsd=rmsd, legacy=legacy, prior=pri.astype(np.float32),
               snap_index=snap_idx, snap_states=snap.astype(np.int8),
               prior_table=Pemp.astype(np.float32),
               PHI=sp.PHI.astype(np.float32), PSI=sp.PSI.astype(np.float32),
               amber_idx=a_idx.astype(np.int64), amber_kind=a_kind,
               amber_total=amb_tot.astype(np.float32))
    for t in Q.LEGACY_TERMS:
        out["leg_" + t] = comp[t].astype(np.float32)
    for t, v in amb_comp.items():
        out["amb_" + t] = v
    path = os.path.join(Q.RESULTS, f"qarch_enum_{pdb_id}.npz")
    np.savez_compressed(path, **out)

    summ = dict(pdb=pdb_id, n=n, k=K, qubits=sp.n_qubits_binary, configs=int(B),
                rmsd_min=float(rmsd.min()), rmsd_max=float(rmsd.max()),
                rmsd_mean=float(rmsd.mean()), rmsd_median=float(np.median(rmsd)),
                snap_rmsd=float(rmsd[snap_idx]),
                legacy_argmin_rmsd=float(rmsd[int(np.argmin(legacy))]),
                prior_argmin_rmsd=float(rmsd[int(np.argmin(pri))]),
                amber_n=int(len(a_idx)),
                amber_argmin_rmsd=(float(rmsd[a_idx[int(np.nanargmin(amb_tot))]])
                                   if amber else None),
                wall_s=round(time.time() - t0, 1), path=path)
    print("  " + json.dumps(summ), flush=True)
    return summ


def main(pdbs=None, amber=True):
    pdbs = pdbs or SMALL
    rows = []
    for p in pdbs:
        rows.append(enumerate_target(p, amber=amber))
        Q.write("qarch_enum", {"what": "full k=4 enumeration of the n=9 targets",
                               "n_expected": len(pdbs), "per_target": rows})
    return rows


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    main(args or None, amber="--no-amber" not in sys.argv)
