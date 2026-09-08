"""SPRINT 13 PAULI-SPECTRUM, EXPERIMENT W2/W3 -- THE EXACT SPECTRUM OF GENUINE AMBER.

`core.amber.single_point` is ff14SB + GBn2 with NO minimisation, on the structure built from
the discrete torsion configuration.  It is a function of the bitstring, so it too is exactly
a weighted sum of Pauli-Z strings.  At ~13-28 ms per DISTINCT configuration the whole
register is enumerable up to m = 14 (16,384 configurations, ~4-8 min); above that the
unbiased Krawtchouk pair estimator of `walsh_lib.sampled_spectrum` is used and its
sampling error reported.

Per-term spectra come from `single_point(..., components=True)`:
    bond, angle, torsion, nonbonded, solvation
and the brief's "AMBER with its steric term removed" is `total - nonbonded` (OpenMM's
NonbondedForce carries the Lennard-Jones steric repulsion together with the Coulomb term;
there is no way to remove only the r^-12 wall through this API, and the arm is named
accordingly).  A CAPPED variant, E_cap = min(E, percentile_99(E)), is also reported: it is a
DIFFERENT observable and is labelled as one, but an unbounded observable has unbounded Pauli
norm and no VQE can be run on one, so its spectrum is the operationally relevant number.

AMBER's 92% memory guard is polled, not caught.

    python -m s13.walsh_amber [--time] [--exact] [--sampled]
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

from s13 import walsh_lib as W
from s13 import qarch_lib as Q

TARGETS = ["1A13", "1A1P", "2BFI", "1DEP", "1ID6", "1CEK"]
# Ordered by priority: the per-term decomposition (coordinator's instruction) first, the
# length-scaling ladder afterwards, the expensive m=16 exact cells last.
EXACT_CELLS = [((6, 4), TARGETS), ((7, 4), TARGETS), ((5, 4), TARGETS),
               ((4, 4), TARGETS), ((12, 2), TARGETS[:3]), ((8, 2), TARGETS[:3]),
               ((4, 8), TARGETS[:3]), ((8, 4), TARGETS[:2])]
SAMPLED_CELLS = [((8, 4), TARGETS[:2]), ((9, 4), TARGETS[:2])]
MEM_CEIL = 90.5                      # AMBER raises above 92; keep a margin


def guard(tag=""):
    from core import amber as am
    while True:
        p = am.memory_percent()
        if p <= MEM_CEIL:
            return p
        print(f"  [{tag}] memory {p:.0f}% > {MEM_CEIL}%; waiting", flush=True)
        time.sleep(20)


def amber_table(sp, S, components=True, tag=""):
    """Genuine ff14SB/GBn2 single points over the given configurations."""
    from core import amber as am
    n = len(S)
    tot = np.empty(n)
    comps = {t: np.empty(n) for t in am.AMBER_TERMS} if components else None
    t0 = time.time()
    for b in range(n):
        if b % 512 == 0:
            guard(tag)
        r = am.single_point(sp.seq, sp.rep, [int(x) for x in S[b]],
                            components=components, threads=1)
        tot[b] = float(r["energy"])
        if components:
            for t in am.AMBER_TERMS:
                comps[t][b] = float(r["components"][t])
        if b and b % 2000 == 0:
            el = time.time() - t0
            print(f"    [{tag}] {b}/{n}  {1000*el/b:.1f} ms/cfg  "
                  f"eta {(n-b)*el/b/60:.1f} min", flush=True)
    return comps, tot, time.time() - t0


def summarise(tot, comps, m, b, tag):
    out = {}
    spec = W.spectrum(tot, m, bits_per_res=b, want_top=8)
    out["total"] = spec
    cap = float(np.percentile(tot, 99))
    out["capped_p99"] = dict(W.spectrum(np.minimum(tot, cap), m, bits_per_res=b),
                             cap_value=cap, MODIFIED_OBSERVABLE=True)
    cap95 = float(np.percentile(tot, 95))
    out["capped_p95"] = dict(W.spectrum(np.minimum(tot, cap95), m, bits_per_res=b),
                             cap_value=cap95, MODIFIED_OBSERVABLE=True)
    out["dynamic_range_orders"] = float(np.log10(max(tot.max() - tot.min(), 1e-300))
                                        - np.log10(max(abs(np.median(tot)), 1e-300)))
    if comps is not None:
        vt = spec["var"]
        out["terms"] = {}
        for name, e in comps.items():
            s = W.spectrum(np.asarray(e, float), m, bits_per_res=b)
            s["var_share"] = s["var"] / vt if vt > 0 else float("nan")
            s["cov_share"] = float(np.cov(np.asarray(e, float),
                                          np.asarray(tot, float))[0, 1] / vt) \
                if vt > 0 else float("nan")
            s.pop("top_terms", None)
            out["terms"][name] = s
        nb = np.asarray(comps["nonbonded"], float)
        minus = np.asarray(tot, float) - nb
        out["minus_nonbonded"] = W.spectrum(minus, m, bits_per_res=b)
        out["sum_of_components_max_abs_err"] = float(
            np.abs(sum(np.asarray(v, float) for v in comps.values())
                   - np.asarray(tot, float)).max())
    return out


def exact_cell(pdb_id, L, k, components=True):
    sp = W.space_for(pdb_id, L, k)
    b, m = sp.bits_per_res, L * sp.bits_per_res
    S = W.enumerate_states(L, k)
    tag = f"{pdb_id}/{L}/{k}"
    comps, tot, dt = amber_table(sp, S, components=components, tag=tag)
    np.save(os.path.join(W.CACHE, f"walsh_E_{pdb_id}_{L}_{k}_amber.npy"), tot)
    if comps is not None:
        np.savez(os.path.join(W.CACHE, f"walsh_C_{pdb_id}_{L}_{k}_amber.npz"), **comps)
    row = {"pdb": pdb_id, "L": L, "k": k, "m_bits": m, "bits_per_res": b,
           "n_configs": int(len(S)), "seconds": dt,
           "ms_per_config": 1000.0 * dt / len(S), "exact": True}
    row.update(summarise(tot, comps, m, b, tag))
    s = row["total"]
    print(f"  AMBER {tag} m={m}: {dt/60:.1f} min  W_mean {s['W_mean']:.3f} "
          f"D_mean {s['D_mean']:.3f} V1 {s['V1']:.4f} tail>3 {s['tail_gt3']:.4f} "
          f"| capped99 W {row['capped_p99']['W_mean']:.3f} "
          f"| parseval {s['parseval_rel_err']:.1e}", flush=True)
    return row


def sampled_cell(pdb_id, L, k, N=8192, seed=0):
    sp = W.space_for(pdb_id, L, k)
    b, m = sp.bits_per_res, L * sp.bits_per_res
    rng = np.random.default_rng(seed)
    j = rng.integers(0, 1 << m, size=N).astype(np.uint64)
    # decode bit index -> state vector (binary encoding, residue r owns bits [r*b, r*b+b))
    S = np.empty((N, L), dtype=np.int64)
    for r in range(L):
        S[:, r] = (j.astype(np.int64) >> ((L - 1 - r) * b)) & ((1 << b) - 1)
    tag = f"{pdb_id}/{L}/{k}/sampled"
    comps, tot, dt = amber_table(sp, S, components=True, tag=tag)
    _, leg = W.legacy_tables(sp, S)
    est = W.sampled_spectrum(j, tot, m, blocks=8, seed=1)
    cap = float(np.percentile(tot, 99))
    est_cap = W.sampled_spectrum(j, np.minimum(tot, cap), m, blocks=8, seed=1)
    est_leg = W.sampled_spectrum(j, np.asarray(leg, float), m, blocks=8, seed=1)
    row = {"pdb": pdb_id, "L": L, "k": k, "m_bits": m, "N": N, "exact": False,
           "seconds": dt, "ms_per_config": 1000.0 * dt / N,
           "amber": est, "amber_capped_p99": dict(est_cap, cap_value=cap),
           "legacy_same_samples": est_leg}
    for name, e in comps.items():
        row[f"term_{name}"] = W.sampled_spectrum(j, np.asarray(e, float), m, blocks=8, seed=1)
    print(f"  AMBER-SAMPLED {tag} m={m} N={N}: {dt/60:.1f} min  "
          f"W {est['W_mean']:.3f}+-{est['W_mean_sd']:.3f}  "
          f"capped W {est_cap['W_mean']:.3f}+-{est_cap['W_mean_sd']:.3f}  "
          f"legacy W {est_leg['W_mean']:.3f}+-{est_leg['W_mean_sd']:.3f}", flush=True)
    return row


def timing():
    from core import amber as am
    guard("time")
    sp = W.space_for("1A13", 6, 4)
    rng = np.random.default_rng(0)
    S = rng.integers(0, 4, size=(120, 6))
    out = {}
    for comp in (False, True):
        am.single_point(sp.seq, sp.rep, [int(x) for x in S[0]], components=comp, threads=1)
        t0 = time.time()
        for b in range(1, 101):
            am.single_point(sp.seq, sp.rep, [int(x) for x in S[b]],
                            components=comp, threads=1)
        out[f"components={comp}"] = 1000.0 * (time.time() - t0) / 100
    out["memory_percent"] = am.memory_percent()
    print("  AMBER ms/distinct-config:", out, flush=True)
    return out


def main(argv):
    do_time = "--time" in argv or len(argv) == 0
    do_exact = "--exact" in argv or len(argv) == 0
    do_samp = "--sampled" in argv or len(argv) == 0
    out = {}
    if do_time:
        out["timing"] = timing()
        W.write("walsh_amber", out)
    if do_exact:
        out["exact"] = []
        for (L, k), tg in EXACT_CELLS:
            for p in tg:
                guard("cell")
                try:
                    out["exact"].append(exact_cell(p, L, k))
                except Exception as exc:                      # noqa: BLE001
                    print(f"  !! {p} L={L} k={k}: {type(exc).__name__}: {exc}", flush=True)
                    out["exact"].append({"pdb": p, "L": L, "k": k,
                                         "error": f"{type(exc).__name__}: {exc}"})
                W.write("walsh_amber", out)
    if do_samp:
        out["sampled"] = []
        for (L, k), tg in SAMPLED_CELLS:
            for p in tg:
                guard("cell")
                try:
                    out["sampled"].append(sampled_cell(p, L, k))
                except Exception as exc:                      # noqa: BLE001
                    print(f"  !! {p} L={L} k={k}: {type(exc).__name__}: {exc}", flush=True)
                W.write("walsh_amber", out)
    print("wrote", W.write("walsh_amber", out))


if __name__ == "__main__":
    main(sys.argv[1:])
