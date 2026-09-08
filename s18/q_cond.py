"""SPRINT 18 / QUANTUM-ADVERSARIAL -- A4: the pre-registered quantum conditional, measured.

PRE-REGISTERED IN `s18/PREREG_quantum.md` BEFORE THIS FILE EXISTED.

    If the new objective contains useful higher-order correlations, then
    (a) greedy 1-opt should CEASE to certify its optimum at tiny budget, and
    (b) the CNOT-free product ansatz should become measurably WORSE than the entangled VQE at
        matched budget, with a target-level interval excluding zero.

`budget` mode measures (a) on all three objectives -- `full`, `RA`, `W1` -- so the question
"is the new objective harder to search than the old one" is answered by a paired comparison
rather than by an assertion.  `vqe` mode measures (b) on `RA`, the object the brief's section
4 defines and the only one with a continuous-torsion analogue.

The zero-information control (uniform random over the register) and the matched-random control
(an equal-count uniform draw) are on every arm.

RUN:  python -m s18.q_cond budget
      python -m s18.q_cond vqe
"""
from __future__ import annotations

import gc
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s14 import vqe_lib as V                       # noqa: E402
from s14 import vqe_run as R                       # noqa: E402
from s15 import seed as SD                         # noqa: E402
from s16 import qphase_lib as QP                   # noqa: E402
from s17 import q_lib as L                         # noqa: E402
from s17 import q_pareto as PA                     # noqa: E402
from s18 import q_anova as A                       # noqa: E402

RESULTS = os.path.join(ROOT, "s18", "results")
os.makedirs(RESULTS, exist_ok=True)
SALT = "s18quantum"
BUDGETS = (20, 36, 64, 128, 256, 512, 1024, 2048, 8192)
SEEDS = (0, 1, 2, 3)
VQE_ALPHAS = (0.05, 0.25, 1.00)


def free_gb():
    import subprocess
    try:
        r = subprocess.run(
            ["powershell", "-NoProfile", "-Command",
             "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"],
            capture_output=True, text=True, timeout=60)
        return float(r.stdout.strip()) / 1048576.0
    except Exception:
        return float("nan")


def build_objects(ins):
    """`full`, `RA`, `W1` for one target, all rank-uniformised as in Sprint 17."""
    E = V.uniformise(ins.hamil().astype(np.float64))
    n, k = ins.n, ins.k
    pw = k ** np.arange(n - 1, -1, -1)
    ra = A.residue_additive(E, n, k, pw)
    # W1 without a full FWHT: the weight-1 coefficients are functions of the marginals alone
    nq = ins.n_qubits
    bits = ins.bits_per_res
    idx = np.arange(ins.N, dtype=np.int64)
    w1 = np.full(ins.N, float(E.mean()))
    codes = np.arange(k)
    for i in range(n):
        s = (idx // pw[i]) % k
        m = np.bincount(s, weights=E, minlength=k) / np.bincount(s, minlength=k)
        for b in range(bits):
            bitpos = bits - 1 - b
            chi = 1.0 - 2.0 * ((codes >> bitpos) & 1)
            c = float((m * chi).mean())
            w1 += c * chi[s]
        del s
    return {"full": E, "RA": ra, "W1": w1}


# ============================================================ (a) the budget floor
def budget_cell(pdb, budgets=BUDGETS, seeds=SEEDS):
    ins = QP.inst(pdb)
    objs = build_objects(ins)
    out = {"pdb": pdb, "n": ins.n, "k": ins.k, "N": int(ins.N), "fold": ins.fold,
           "budgets": list(budgets), "seeds": list(seeds), "budget": {}}
    for name, E in objs.items():
        opt = float(E.min())
        n_opt = int((E <= opt + 1e-12).sum())
        row = {"one_pass_cost": int(ins.n * (ins.k - 1) + 1),
               "closed_form_cost": (int(ins.n * ins.k) if name == "RA" else
                                    (int(2 * ins.n_qubits) if name == "W1" else None)),
               "n_optimal_configs": n_opt,
               "certified_argmin_rmsd_ORACLE": float(ins.rmsd[E <= opt + 1e-12].mean()),
               "greedy_hit": [], "metro_hit": [], "rand_hit": [],
               "greedy_rmsd_ORACLE": [], "metro_rmsd_ORACLE": []}
        for b in budgets:
            gh, mh, rh, gr, mr = [], [], [], [], []
            for s in seeds:
                c = V.search_greedy(E, ins.n, ins.k, b,
                                    SD.stable_rng(pdb, s, f"g{name}{b}", salt=SALT))
                gh.append(bool(c.best_e <= opt + 1e-12)); gr.append(float(ins.rmsd[c.best_i]))
                c = L.search_metropolis(E, ins.n, ins.k, b, 0.1 * float(E.std()),
                                        SD.stable_rng(pdb, s, f"m{name}{b}", salt=SALT))
                mh.append(bool(c.best_e <= opt + 1e-12)); mr.append(float(ins.rmsd[c.best_i]))
                rr = SD.stable_rng(pdb, s, f"r{name}{b}", salt=SALT)
                pick = rr.integers(0, ins.N, size=b)
                rh.append(bool(E[pick].min() <= opt + 1e-12))
            row["greedy_hit"].append(float(np.mean(gh)))
            row["metro_hit"].append(float(np.mean(mh)))
            row["rand_hit"].append(float(np.mean(rh)))
            row["greedy_rmsd_ORACLE"].append(float(np.mean(gr)))
            row["metro_rmsd_ORACLE"].append(float(np.mean(mr)))
        out["budget"][name] = row
    del objs
    gc.collect()
    return out


def run_budget(targets=None):
    targets = list(targets or QP.TARGETS19)
    have = A.load("qcond")
    for p in targets:
        if p in have:
            print(f"[skip] {p}", flush=True)
            continue
        fg = free_gb()
        if fg < 0.9:
            print(f"[STOP] free RAM {fg:.2f} GB < 0.90; stopping per the brief's rule",
                  flush=True)
            return
        t0 = time.time()
        r = budget_cell(p)
        A.ck("qcond", p, r)
        bi = list(BUDGETS).index(1024)
        print(f"[{p}] {time.time()-t0:5.1f}s free {fg:.2f}GB  greedy@1024: " +
              "  ".join(f"{o} {100*r['budget'][o]['greedy_hit'][bi]:.0f}%"
                        for o in ("full", "RA", "W1")), flush=True)
    A.ck("qcond", "_complete", len(A.load("qcond")) >= len(QP.TARGETS19))


# ============================================================ (b) entangled vs CNOT-free
def vqe_cell(pdb, seed=0, obj="RA", alphas=VQE_ALPHAS):
    ins = QP.inst(pdb)
    E = build_objects(ins)[obj]
    rng = SD.stable_rng(pdb, seed, "readout", salt=SALT)
    res = {"_meta": {"pdb": pdb, "seed": seed, "fold": ins.fold, "n": ins.n, "obj": obj,
                     "budget": L.BUDGET, "shots": L.SHOTS, "lr": 0.15,
                     "baseline": "const"}}
    for ans in ("mps2f", "mps2fn"):
        for al in alphas:
            c = L.Cost()
            v = R.run(E, ins.n_qubits, al, L.BUDGET, shots=L.SHOTS, ansatz=ans, seed=seed,
                      rmsd=ins.rmsd, bits_per_res=ins.bits_per_res, exact_dist=False,
                      keep_seen=True, lr=0.15, baseline="const")
            seen = np.asarray(v["seen"], np.int64)
            r, _ = PA.arm_readout(ins, E, seen, rng,
                                  c.stop(v["evals"], circuit_samples=v["evals"],
                                         grad_passes=v["iters"],
                                         distinct=int(np.unique(seen).size)))
            r["iters"] = int(v["iters"])
            res[f"{ans}_a{al}"] = r
    # controls
    for nm, pick in (("random", SD.stable_rng(pdb, seed, "rnd", salt=SALT)
                      .integers(0, ins.N, size=L.BUDGET)),):
        c = L.Cost()
        r, _ = PA.arm_readout(ins, E, np.asarray(pick, np.int64), rng,
                              c.stop(L.BUDGET, circuit_samples=0, grad_passes=0,
                                     distinct=int(np.unique(pick).size)))
        res[nm] = r
    c = L.Cost()
    cg = V.search_greedy(E, ins.n, ins.k, L.BUDGET,
                         SD.stable_rng(pdb, seed, "gr", salt=SALT))
    sg = np.asarray(cg.all_seen(), np.int64)
    r, _ = PA.arm_readout(ins, E, sg, rng,
                          c.stop(sg.size, circuit_samples=0, grad_passes=0,
                                 distinct=int(np.unique(sg).size)))
    res["greedy"] = r
    del E
    gc.collect()
    return res


def run_vqe(targets=None, seed=0, obj="RA"):
    targets = list(targets or QP.TARGETS19)
    tag = f"qvqe_{obj}"
    have = A.load(tag)
    for p in targets:
        key = f"{p}_{seed}"
        if key in have:
            print(f"[skip] {key}", flush=True)
            continue
        fg = free_gb()
        if fg < 0.9:
            print(f"[STOP] free RAM {fg:.2f} GB < 0.90; stopping per the brief's rule",
                  flush=True)
            return
        t0 = time.time()
        r = vqe_cell(p, seed, obj)
        A.ck(tag, key, r)
        print(f"[{key}] {time.time()-t0:5.0f}s free {fg:.2f}GB  " +
              "  ".join(f"{a}:{r[f'mps2f_a{a}']['M75']:.3f}/{r[f'mps2fn_a{a}']['M75']:.3f}"
                        for a in VQE_ALPHAS), flush=True)
    A.ck(tag, "_complete", len(A.load(tag)) >= len(QP.TARGETS19))


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else "budget"
    rest = sys.argv[2:]
    if m == "budget":
        run_budget(rest or None)
    elif m == "vqe":
        obj = rest[0] if rest and rest[0] in ("RA", "W1", "full") else "RA"
        tg = [x for x in rest if x not in ("RA", "W1", "full")]
        run_vqe(tg or None, obj=obj)
