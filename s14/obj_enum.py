"""SPRINT 14, OBJ -- extend the fully-enumerated k=4 ground-truth set.

s13 enumerated the nine n=9 targets (4^9 = 262,144 configs each).  Nine targets is not
enough for leave-fold-out training with held-out targets in every fold, so this module
enumerates every tuning126 target with n <= 11 at k = 4:

    n=9   262,144 configs   9 targets   (already cached in s13/results)
    n=10  1,048,576         11 targets
    n=11  4,194,304         13 targets

giving 33 targets, >= 6 per pinned fold, and 6.9e7 exactly-labelled structures.

For every configuration we cache, exactly as s13 did:
    rmsd      ORACLE post-hoc CA-RMSD of the ideal-geometry CA trace   (LABEL ONLY)
    legacy    DEFAULT_WEIGHTS total, plus all 11 core.energy terms
    prior     1-local torsion-prior energy under the leakage-safe empirical prior
plus a labelled AMBER subsample (uniform / prior-sampled / near-native band).

Nothing here is used at inference.  `rmsd` is a training label and an evaluation axis.

    python -m s14.obj_enum 10          # every n=10 target
    python -m s14.obj_enum 11          # every n=11 target
    python -m s14.obj_enum 1N9U 1TOR   # named targets

Output: s14/cache/obj_enum_<PDB>.npz  (same schema as s13/results/qarch_enum_<PDB>.npz)
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s13 import qarch_lib as Q          # noqa: E402
from s12 import instrument as I         # noqa: E402

K = 4
CACHE = os.path.join(ROOT, "s14", "cache")
RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

N_AMBER_UNIFORM = 700
N_AMBER_PRIOR = 700
N_AMBER_BAND = 400


def enum_configs(n, k, lo, hi):
    """Rows [lo, hi) of the odometer enumeration of {0..k-1}^n, as (hi-lo, n) int8.

    Row index r decodes as s_i = (r // k^(n-1-i)) % k, identical to
    np.array(list(itertools.product(range(k), repeat=n))) but vectorised.
    """
    r = np.arange(lo, hi, dtype=np.int64)
    out = np.empty((hi - lo, n), dtype=np.int8)
    for i in range(n):
        out[:, i] = (r // (k ** (n - 1 - i))) % k
    return out


def enum_path(pdb_id):
    p = os.path.join(CACHE, f"obj_enum_{pdb_id}.npz")
    if os.path.exists(p):
        return p
    q = os.path.join(ROOT, "s13", "results", f"qarch_enum_{pdb_id}.npz")
    return q if os.path.exists(q) else p


def have(pdb_id):
    return os.path.exists(enum_path(pdb_id))


def load_enum(pdb_id):
    return np.load(enum_path(pdb_id))


def enumerate_target(pdb_id, seed=0, amber=True, chunk=131072, skip_existing=True):
    t0 = time.time()
    dest = os.path.join(CACHE, f"obj_enum_{pdb_id}.npz")
    if skip_existing and have(pdb_id):
        print(f"{pdb_id}: cached at {enum_path(pdb_id)}", flush=True)
        return dict(pdb=pdb_id, cached=True, path=enum_path(pdb_id))

    sp = Q.Space(pdb_id, K)
    n = sp.n
    B = K ** n
    # 14 float32 columns over B configs, plus decode scratch of `chunk` rows.
    need_gb = (14 * 4 * B) / 1e9 + 0.6
    Q.wait_for_memory(max(1.2, need_gb), pdb_id)
    print(f"{pdb_id}: n={n} k={K} configs={B} qubits={sp.n_qubits_binary} "
          f"(reserving ~{need_gb:.2f} GB)", flush=True)

    rmsd = np.empty(B, np.float32)
    legacy = np.empty(B, np.float32)
    comps = {t: np.empty(B, np.float32) for t in Q.LEGACY_TERMS}
    for a in range(0, B, chunk):
        b = min(a + chunk, B)
        S = enum_configs(n, K, a, b)
        rmsd[a:b] = sp.rmsd(S)
        c = Q.legacy_components(sp, S, chunk=8192)
        for t in Q.LEGACY_TERMS:
            comps[t][a:b] = c[t]
        legacy[a:b] = Q.legacy_total(c)
        del S, c
        if (a // chunk) % 8 == 0:
            print(f"  {b}/{B} [{time.time()-t0:.0f}s]", flush=True)
    print(f"  rmsd+legacy done [{time.time()-t0:.0f}s] rmsd_min={rmsd.min():.3f}",
          flush=True)

    Pemp = Q.empirical_prior(sp)
    Lp = -np.log(Pemp).astype(np.float32)                       # (n, k)
    prior = np.zeros(B, np.float32)
    for a in range(0, B, chunk):
        b = min(a + chunk, B)
        S = enum_configs(n, K, a, b)
        prior[a:b] = Lp[np.arange(n)[None, :], S].sum(1)
        del S
    snap = sp.ORACLE_snap()
    snap_idx = int(np.ravel_multi_index(tuple(snap), (K,) * n))

    rng = np.random.default_rng(seed)
    idx_u = rng.choice(B, N_AMBER_UNIFORM, replace=False)
    Sp = sp.sample_prior(Pemp, N_AMBER_PRIOR, rng)
    idx_p = np.ravel_multi_index(tuple(Sp.T), (K,) * n)
    band = np.argpartition(rmsd, max(2000, B // 100))[:max(2000, B // 100)]
    idx_b = rng.choice(band, N_AMBER_BAND, replace=False)
    a_idx = np.unique(np.concatenate([idx_u, idx_p, idx_b, [snap_idx]]))
    a_kind = np.zeros(len(a_idx), np.int8)
    a_kind[np.isin(a_idx, idx_p)] = 1
    a_kind[np.isin(a_idx, idx_b)] = 2

    amb_tot = np.full(len(a_idx), np.nan, np.float64)
    amb_comp = {}
    if amber:
        Q.wait_for_memory(1.2, pdb_id + "/amber")
        ta = time.time()
        Sa = np.stack([enum_configs(n, K, int(j), int(j) + 1)[0] for j in a_idx])
        c, tt = Q.amber_energies(sp, Sa, components=True, progress=400)
        amb_tot = tt
        amb_comp = {kk: v.astype(np.float32) for kk, v in c.items()}
        print(f"  amber {len(a_idx)} in {time.time()-ta:.0f}s "
              f"({(time.time()-ta)/len(a_idx)*1000:.1f} ms/cfg)", flush=True)

    out = dict(pdb=pdb_id, n=n, k=K, seq=sp.seq, fold=sp.fold,
               rmsd=rmsd, legacy=legacy, prior=prior,
               snap_index=snap_idx, snap_states=snap.astype(np.int8),
               prior_table=Pemp.astype(np.float32),
               PHI=sp.PHI.astype(np.float32), PSI=sp.PSI.astype(np.float32),
               amber_idx=a_idx.astype(np.int64), amber_kind=a_kind,
               amber_total=amb_tot.astype(np.float32))
    for t in Q.LEGACY_TERMS:
        out["leg_" + t] = comps[t]
    for t, v in amb_comp.items():
        out["amb_" + t] = v
    tmp = dest + ".tmp.npz"
    np.savez(tmp, **out)          # uncompressed: 3-5x faster, disk is not the constraint
    os.replace(tmp, dest)

    summ = dict(pdb=pdb_id, n=n, k=K, qubits=sp.n_qubits_binary, configs=int(B),
                fold=int(sp.fold), seq=sp.seq,
                rmsd_min=float(rmsd.min()), rmsd_max=float(rmsd.max()),
                rmsd_mean=float(rmsd.mean()), rmsd_median=float(np.median(rmsd)),
                snap_rmsd=float(rmsd[snap_idx]),
                legacy_argmin_rmsd=float(rmsd[int(np.argmin(legacy))]),
                prior_argmin_rmsd=float(rmsd[int(np.argmin(prior))]),
                amber_n=int(len(a_idx)),
                amber_argmin_rmsd=(float(rmsd[a_idx[int(np.nanargmin(amb_tot))]])
                                   if amber else None),
                wall_s=round(time.time() - t0, 1), path=dest)
    print("  " + json.dumps(summ), flush=True)
    del rmsd, legacy, comps, prior
    return summ


def targets_with_n(nmax):
    return [t["pdb"] for t in I.targets() if t["n"] <= nmax]


def targets_exactly(n):
    return [t["pdb"] for t in I.targets() if t["n"] == n]


def main(args):
    amber = "--no-amber" not in args
    args = [a for a in args if not a.startswith("-")]
    if len(args) == 1 and args[0].isdigit():
        pdbs = targets_exactly(int(args[0]))
    elif args:
        pdbs = args
    else:
        pdbs = targets_with_n(11)
    rows = []
    for p in pdbs:
        rows.append(enumerate_target(p, amber=amber))
        with open(os.path.join(RESULTS, "obj_enum_log.json"), "w") as fh:
            json.dump({"n_expected": len(pdbs), "done": len(rows),
                       "per_target": rows}, fh, indent=1, default=str)
    return rows


if __name__ == "__main__":
    main(sys.argv[1:])
