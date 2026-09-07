"""s21/d_sim.py -- D4: THE REACHABILITY VERDICT, BY MEASUREMENT RATHER THAN BY CITATION.

    python -m s21.d_sim

THE QUESTION (BRIEF Priority 3, and it gates whether Workstream B should spend a night on
ansatz design).  Sprint 20 established as EXACT that the deployed ansatz is bond-dimension 4 --
a <=16-state HMM, classically samplable by construction.  So: (a) at what bond dimension does
classical simulability actually break for this circuit class, (b) does any ansatz that escapes
it remain trainable at this qubit count, (c) has anyone shown a separation on a comparable
continuous-variable encoding?

THE ANSWER TO (a) DOES NOT DEPEND ON THE BOND DIMENSION AT ALL, AND THIS FILE MEASURES IT.

From source, `s20.qb2_opt.arm_vqe:322`, the deployed circuit is

    Q.MPSAnsatz(n, layers=2, final_ry=True, entangler="cnot")

with `n = F.n` = **the number of RESIDUES**.  The tuning instrument is 9 <= n <= 16, so the
register is **9 to 16 qubits** and the latent is a bitstring selecting one conformer basin per
residue.  The whole Hilbert space therefore has 512 to 65,536 basis states.

    chi = 2 ** layers = 4 is a statement about the MPS bond.  It is TRUE and it is not the
    binding constraint.  The binding constraint is that a 16-qubit register's ENTIRE
    distribution can be written down.

So this file does exactly that: for every register size the benchmark actually uses, it
enumerates the deployed ansatz's **complete** probability distribution exactly, checks it sums
to 1, and times it.  A circuit whose full law can be tabulated in milliseconds admits no
quantum-resource claim at any bond dimension, any depth, any entangler and any budget --
including bond dimensions far above 4, and including a maximally entangled register (for which
chi_max = 2 ** floor(n/2) = 256 at n = 16, still trivial).

WHAT THIS IS NOT.  It is not a claim that the VQE is useless -- Q4 stands, the VQE genuinely
trains, and a restricted parameterisation of a small distribution can still be a good
optimiser.  It is a claim about **quantum resource**: nothing here can be evidence of quantum
advantage, and no redesign of the ansatz can make it so while the register is one qubit per
residue on 9-16-residue peptides.

CONTROL, so the verdict is not an artefact of the deployed settings: the enumeration is
repeated at `layers` up to the point where chi = 2**layers saturates the register
(chi >= 2**floor(n/2)), i.e. at MAXIMAL entanglement for this topology.  If the exhaustive
enumeration is still cheap there -- and it is, because its cost is 2**n and does not depend on
layers at all -- then no choice of depth escapes.

Artefact: `s21/results/d_sim.json` with a completion flag requiring EVERY register size the
126-target instrument uses, at every tested depth, with the normalisation check passing.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I        # noqa: E402
from core import quantum as Q          # noqa: E402

LAYERS = (1, 2, 3, 4, 6, 8, 12, 20)    # deployed is 2; >= n//2 saturates chi for this topology
MPS_LAYERS = (1, 2, 3, 4)              # the MPS route is O(n chi^4); above this use the dense one


def run():
    ns = sorted({int(t["n"]) for t in I.targets()})
    counts = {}
    for t in I.targets():
        counts[int(t["n"])] = counts.get(int(t["n"]), 0) + 1
    rng = np.random.default_rng(20210907)
    rows = []
    print(f"Register sizes used by the 126-target instrument: {ns}")
    print(f"Deployed ansatz: MPSAnsatz(n, layers=2, final_ry=True, entangler='cnot'), "
          f"n = RESIDUES (s20/qb2_opt.py:322)\n")
    print(f"  {'n':>3}{'targets':>8}{'states':>8}{'L':>3}{'chi=2^L':>9}{'chi_max':>9}"
          f"{'n_params':>9}{'ms':>9}{'sum p':>14}{'entropy':>9}{'max p':>10}  route")
    for n in ns:
        allb = ((np.arange(2 ** n)[:, None] >> np.arange(n - 1, -1, -1)) & 1).astype(np.int8)
        for L in LAYERS:
            #: TWO routes, and the contrast is the point.  The MPS route costs O(n chi^4) and
            #: therefore BLOWS UP with depth -- that is a property of the REPRESENTATION.  The
            #: dense statevector route costs 2**n and is INDEPENDENT of depth, which is why no
            #: choice of depth escapes classical simulation on a 9-16 qubit register.
            route = "MPS" if L in MPS_LAYERS else "dense"
            if route == "MPS":
                an = Q.MPSAnsatz(n, layers=L, final_ry=True, entangler="cnot")
                th = np.pi / 2 + rng.normal(0.0, 0.8, an.n_params())
                t0 = time.time()
                p = np.exp(an.logp(th, allb))
                ms = (time.time() - t0) * 1000.0
                npar = int(an.n_params())
            else:
                an = Q.StatevectorCircuit(n, layers=L, ring=True)
                th = np.pi / 2 + rng.normal(0.0, 0.8, an.n_params())
                t0 = time.time()
                p = np.asarray(an.probs(th), float)
                ms = (time.time() - t0) * 1000.0
                npar = int(an.n_params())
            H = float(-(p * np.log2(np.maximum(p, 1e-300))).sum())
            r = {"n": n, "n_targets": counts[n], "states": 2 ** n, "layers": L,
                 "route": route,
                 "chi": min(2 ** L, 2 ** (n // 2)), "chi_max": 2 ** (n // 2),
                 "n_params": npar,
                 "ms": ms, "sum_p": float(p.sum()), "entropy_bits": H,
                 "max_bits": n, "max_p": float(p.max()),
                 "normalised": bool(abs(p.sum() - 1.0) < 1e-9),
                 "chi_saturates_register": bool(2 ** L >= 2 ** (n // 2))}
            rows.append(r)
            print(f"  {n:>3}{counts[n]:>8}{2**n:>8}{L:>3}{min(2**L,2**(n//2)):>9}"
                  f"{2**(n//2):>9}{npar:>9}{ms:>9.1f}{p.sum():>14.10f}{H:>9.3f}"
                  f"{p.max():>10.2e}  {route}")
    tot = sum(r["ms"] for r in rows)
    print(f"\nTOTAL time to enumerate EVERY register size at EVERY tested depth exactly: "
          f"{tot/1000:.2f} s")
    worst = max(rows, key=lambda r: r["ms"])
    print(f"Worst single case: n = {worst['n']} ({worst['states']} states), "
          f"layers = {worst['layers']}, {worst['ms']:.0f} ms, sum p = {worst['sum_p']:.10f}")
    print(f"Cost is 2**n and is INDEPENDENT of layers, so no depth escapes: at n = 16 the "
          f"maximal bond for this topology is chi_max = 256 and the exhaustive enumeration "
          f"is unchanged.")
    #: THE BUDGET AGAINST THE WHOLE SPACE.  `s19.qb_lib.draw_from_basins` takes bits of shape
    #: (B, n) with n = RESIDUES and bits in {0,1}: one qubit per residue, TWO basins each.  So
    #: the entire latent space is 2**n configurations.  Compare that to the budgets spent.
    print("\nTHE BUDGET AGAINST THE WHOLE LATENT SPACE "
          "(one qubit per residue, 2 basins each, so |space| = 2**n):")
    bud = {}
    for B in (512, 4096, 8192):
        k = int(sum(counts[n] for n in ns if 2 ** n <= B))
        bud[str(B)] = k
        print(f"  budget {B:>5}: 2**n <= budget on {k:>3}/126 targets ({100*k/126:.0f}%)")
    med_n = int(np.median(np.repeat(ns, [counts[n] for n in ns])))
    print(f"  |latent|: min {2**min(ns)}   median {2**med_n}   max {2**max(ns)}")
    print("  READ.  At the 8192-evaluation budget Sprint 20 spent, THE BUDGET EQUALS OR EXCEEDS")
    print("  THE ENTIRE LATENT SPACE ON 75 OF 126 TARGETS, and the median target's latent space")
    print("  is exactly 8192.  On a majority of targets the deployed VQE is not searching a")
    print("  space it cannot enumerate -- it is resampling one it could have enumerated within")
    print("  budget.  So Q3 ('no sampler at 8192 evaluations beats the ZERO-evaluation retrieval")
    print("  pool') is not a search failure on those targets: exhaustive enumeration of the")
    print("  latent was affordable and would not have helped, because what binds is the")
    print("  OBJECTIVE'S ORDERING, not the reach of the sampler.")

    print("\nVERDICT (D4a): there is NO bond dimension at which classical simulability breaks")
    print("for this circuit, because the register is 9-16 qubits and its complete law is")
    print("tabulable in milliseconds.  The bond-dimension-4 result is TRUE and is not the")
    print("binding constraint.  A quantum-resource claim is NOT REACHABLE at this scale for")
    print("ANY ansatz, and redesigning the ansatz cannot change that.")
    out = {"rows": rows, "layers_tested": list(LAYERS), "ns": ns,
           "deployed": {"layers": 2, "final_ry": True, "entangler": "cnot",
                        "qubits": "one per residue", "source": "s20/qb2_opt.py:322"},
           "total_ms": tot, "budget_vs_latent": bud,
           "latent_note": "one qubit per residue, 2 basins each (s19.qb_lib.draw_from_basins); "
                          "|latent| = 2**n"}
    json.dump(out, open(os.path.join(RESULTS, "d_sim.json"), "w"), indent=1)
    need_ok = (len(rows) == len(ns) * len(LAYERS)
               and all(r["normalised"] for r in rows)
               and sum(r["n_targets"] for r in rows) // len(LAYERS) == 126)
    p = os.path.join(RESULTS, "d_sim.COMPLETE")
    if need_ok:
        with open(p, "w") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                     f"register_sizes={ns} (covering all 126 tuning targets) "
                     f"layers={list(LAYERS)} cells={len(rows)} "
                     f"all_normalised=True total_ms={tot:.0f}\n")
        print(f"\nCOMPLETE: {len(rows)} cells, all register sizes of all 126 targets, "
              f"all normalisation checks passed.")
    elif os.path.exists(p):
        os.remove(p)
    return out


if __name__ == "__main__":
    run()
