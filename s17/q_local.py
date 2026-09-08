"""SPRINT 17 / QUANTUM -- T4, THE LOCAL VQE REFORMULATION (sprint section 23).

PRE-REGISTERED IN `s17/PREREG_quantum.md` (P3).

THE PROPOSAL.  Instead of asking a VQE to discover global structure, run it on a SMALL LOCAL
CONFORMATIONAL NEIGHBOURHOOD around a retrieval candidate:

    retrieval candidate -> local neighbourhood -> small encoded subproblem
                        -> CVaR-VQE proposal distribution -> new nearby candidates

HYPOTHESIS (P3).  The local proposal distribution reaches an (M, D) point outside what
classical local search reaches on the SAME neighbourhood at the same cost.

EXPECTED OUTCOME, AND WHY THIS IS A TRAP.  A neighbourhood small enough for a local VQE to be
tractable is small enough to EXHAUSTIVELY ENUMERATE.  A window of `w` residues at k = 4 has
`4^w` configurations: 256 at w = 4, 1,024 at w = 5, 4,096 at w = 6.  Exhaustive enumeration
delivers the CERTIFIED optimum and the EXACT Boltzmann law of the subproblem at any
temperature, for `4^w` objective evaluations and no hyperparameters.  The reformulation
removes its own justification, and this module measures exactly that rather than asserting
it: every arm is run at a budget BELOW `4^w` so the comparison is not rigged, and the
certified enumeration is printed as the price of skipping the whole question.

THE CANDIDATE IS NATIVE-FREE.  The base configuration is the mode of the retrieval torsion
prior -- `prior` is exactly additive (verified), so its Boltzmann mode is the per-residue
argmin, a genuine retrieval candidate costing ZERO reads of the deployed objective.  The
zero-information control required by the brief is the base candidate itself (do nothing), and
the matched-random control is a random neighbourhood draw of the same count.

READOUT.  The same (member error, diversity) plane as everything else in this workstream,
plus the set best (ORACLE ceiling of what the proposal put on the table) and the realised
argmin under the objective.

RUN:  python -m s17.q_local run
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from s12 import instrument as I
from s14 import vqe_lib as V
from s14 import vqe_run as R
from s15 import seed as SD
from s17 import q_lib as L

W = 6                     # window width in residues -> 4^6 = 4,096 configurations
SUB_BUDGET = 1024         # a QUARTER of the exhaustive cost, so the VQE is not rigged
SUB_SHOTS = 64
ALPHAS = (0.10, 0.25, 1.00)
SEEDS = (0, 1, 2)
MSET = 40                 # the proposal set a consumer would take


def base_candidate(ins):
    """NATIVE-FREE. The retrieval torsion prior's mode: per-residue argmin of an exactly
    additive potential.  Zero reads of the deployed objective."""
    f, _, _ = L.prior_factors(ins)
    return np.asarray(f.argmin(1), np.int64)


def window_index_map(ins, base_states, a, w):
    """Full-register indices of every configuration in the `w`-residue window at offset `a`,
    with all other residues held at `base_states`.  Index algebra only -- no energy read."""
    n, k = ins.n, ins.k
    pw = k ** np.arange(n - 1, -1, -1)
    base_i = int(base_states @ pw)
    sub = np.arange(k ** w, dtype=np.int64)
    dig = (sub[:, None] // (k ** np.arange(w - 1, -1, -1))[None, :]) % k
    off = ((dig - base_states[a:a + w][None, :]) * pw[a:a + w][None, :]).sum(1)
    return base_i + off


def md_of(ins, full_idx, rng, m=MSET, rep=6):
    acc = {kk: [] for kk in ("M", "D", "readout", "set_best")}
    full_idx = np.asarray(full_idx, np.int64)
    for _ in range(rep):
        pick = rng.choice(full_idx, size=m, replace=full_idx.size < m)
        r = L.md_plane(ins.ca(pick), ins.nat)
        for kk in acc:
            acc[kk].append(r[kk])
    return {kk: (float(np.sqrt(np.mean(np.square(v))))
                 if kk in ("M", "D", "readout") else float(np.mean(v)))
            for kk, v in acc.items()}


def run_cell(pdb, seed, w=W, budget=SUB_BUDGET):
    ins = L.inst(pdb)
    E = ins.hamil()
    rmsd = ins.rmsd                                   # ORACLE
    bs = base_candidate(ins)
    pw = ins.k ** np.arange(ins.n - 1, -1, -1)
    base_i = int(bs @ pw)
    rng = SD.stable_rng(pdb, seed, "local", salt=L.SALT)
    nwin = ins.n - w + 1
    out = {"_meta": {"pdb": pdb, "seed": seed, "n": ins.n, "w": w,
                     "n_windows": int(nwin), "sub_N": int(ins.k ** w),
                     "budget": int(budget),
                     "exhaustive_cost": int(ins.k ** w),
                     "base_rmsd_ORACLE": float(rmsd[base_i]),
                     "base_E": float(E[base_i]),
                     "space_best_ORACLE": float(rmsd.min())}}
    rows = {}
    for a in range(nwin):
        fi = window_index_map(ins, bs, a, w)
        Es = E[fi]                                    # materialised; arms are CHARGED per read
        nq = 2 * w
        cell = {}

        def add(tag, sub_seen, used):
            u = np.unique(np.asarray(sub_seen, np.int64))
            es = Es[np.asarray(sub_seen, np.int64)]
            mn = es.min()
            tied = np.unique(np.asarray(sub_seen)[es == mn])
            d = {"n_distinct": int(u.size), "used": int(used),
                 "best_e": float(mn),
                 "gap_to_certified": float(mn - Es.min()),
                 "reaches_certified": bool(mn <= Es.min() + 1e-12),
                 "argmin_rmsd_ORACLE": float(rmsd[fi[tied]].mean()),
                 "set_best_ORACLE": float(rmsd[fi[u]].min())}
            d.update(md_of(ins, fi[np.asarray(sub_seen, np.int64)], rng))
            cell[tag] = d

        # ---- quantum proposal distributions --------------------------------
        for al in ALPHAS:
            v = R.run(Es, nq, al, budget, shots=SUB_SHOTS, ansatz=L.ANSATZ, seed=seed,
                      rmsd=rmsd[fi], bits_per_res=2, exact_dist=False, keep_seen=True,
                      lr=0.15, baseline="const")
            add(f"vqe_a{al}", v["seen"], v["evals"])
        # ---- classical arms on the SAME neighbourhood, same budget ---------
        for tag, fn in (("greedy", V.search_greedy), ("anneal", V.search_anneal),
                        ("random", V.search_random)):
            c = fn(Es, w, ins.k, budget, SD.stable_rng(pdb, seed, f"{tag}{a}", salt=L.SALT))
            add(tag, c.all_seen(), c.used)
        for t in (0.1, 0.5, 2.0):
            c = L.search_metropolis(Es, w, ins.k, budget, t * float(Es.std()),
                                    SD.stable_rng(pdb, seed, f"m{t}{a}", salt=L.SALT))
            add(f"metro_T{t}", c.all_seen(), c.used)
        # ---- THE CONTROL THAT CLOSES IT: exhaustive enumeration ------------
        allsub = np.arange(Es.size)
        add("exhaustive", allsub, int(Es.size))
        # the exact Boltzmann law of the subproblem, at the VQE's own entropy, drawn
        # `budget` times -- available for free once the enumeration is paid for
        for t in (0.2, 1.0, 5.0):
            q = np.exp(-(Es - Es.min()) / max(t * float(Es.std()), 1e-12))
            q /= q.sum()
            s = SD.stable_rng(pdb, seed, f"eb{t}{a}", salt=L.SALT).choice(
                Es.size, size=budget, p=q)
            add(f"exact_boltz_T{t}", s, int(Es.size))
        # ---- zero-information control: the base candidate, unchanged -------
        cell["base_do_nothing"] = {"argmin_rmsd_ORACLE": float(rmsd[base_i]),
                                   "set_best_ORACLE": float(rmsd[base_i]),
                                   "used": 0}
        rows[str(a)] = cell
    out["windows"] = rows
    return out


def run(targets=L.TARGETS9, seeds=SEEDS):
    done = L.ck_load("local")
    for pdb in targets:
        for sd in seeds:
            key = f"{pdb}_{sd}"
            if key in done:
                print(f"  skip {key}", flush=True)
                continue
            L.gate(key)
            t0 = time.time()
            r = run_cell(pdb, sd)
            L.ck("local", key, r)
            done[key] = r
            w0 = r["windows"]["0"]
            print(f"  {key} {time.time()-t0:.0f}s  win0 vqe0.25 M/D "
                  f"{w0['vqe_a0.25']['M']:.3f}/{w0['vqe_a0.25']['D']:.3f} "
                  f"cert {w0['vqe_a0.25']['reaches_certified']}  "
                  f"greedy cert {w0['greedy']['reaches_certified']}", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "run":
        tg = sys.argv[2].split(",") if len(sys.argv) > 2 else L.TARGETS9
        run(tg)
    else:
        raise SystemExit(mode)
