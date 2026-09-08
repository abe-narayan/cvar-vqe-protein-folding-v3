"""s12 QUANTUM-ROLE -- END-TO-END: what the component does to the EMITTED structure.

The assembled point cloud is not an ideal-geometry chain, exactly as the production
coordinate average is not; the production path therefore projects it onto the
ideal-geometry manifold (`I.project`, lam=0 arm = `fit_ca`).  This runs that same path for
every arm, so the numbers are comparable with the shipped 3.2041.

Arms (family A2_64, 12 qubits, 4,096 assemblies, all 126 targets):
  anchor        the shipped `fit_ca` itself (the incumbent, 3.2041)
  exact         the certified exhaustive optimum of the Hamiltonian
  greedy_ls     greedy + 1-opt local search
  vqe           `run_global_cvar_vqe` under a 3,200-unique-evaluation shared budget
  vqe_pweight   the VQE's DISTRIBUTION readout: sampled-frequency-weighted coordinate
                average of the assemblies it visited (the "return a distribution, not a
                point" claim, in the form the pipeline could actually consume)

Usage: python -m s12.vq_e2e [family] [budget]
"""
from __future__ import annotations
import os, sys, json, time, math, warnings
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import vq_lib as V
from s12 import vq_classical as C
from s12 import vq_quantum as Qm
from s12 import vq_qrun as R
from s12.vq_build import FAMILIES
from s12.vq_run import native


def run(fam="A2_64", budget=3200, layers=3, alpha=0.25, shots=16, restarts=1,
        seed=0, out=None):
    from core import quantum as Q
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    cfgf = FAMILIES[fam]
    tg = [t for t in I.targets() if V.segments(t["n"], cfgf["k"], cfgf["minlen"]) is not None]
    rows = []
    t00 = time.time()
    for q, t in enumerate(tg):
        while I.free_gb() < 1.2:
            print("  waiting on memory", flush=True); time.sleep(20)
        inst = V.cached_instance(t, **cfgf)
        E, N, cfgmap, H = R.energy_of(fam, inst)
        nat = native(t["pdb"])
        exh = C.a_exhaustive(inst["h"], inst["J"], inst["npairs"])
        Eall, cfgs = exh.pop("all_E"), exh.pop("all_cfg")
        Xall = V.assemble_batch(inst, cfgs)
        gls = C.a_greedy_ls(inst["h"], inst["J"], inst["npairs"])

        ham = Qm.qubo_hamiltonian(E, N, eval_budget=budget)
        res = Q.run_global_cvar_vqe(ham, layers=layers, alpha=alpha, shots=shots,
                                    restarts=restarts, seed=seed, optimizer="SPSA",
                                    device="lightning.qubit", final_shots=max(256, shots * 4))
        bs = res["best_seen_bitstring"]
        c_vqe = cfgmap[int(bs, 2)] if cfgmap is not None else R.state_to_cfg(fam, inst, int(bs, 2), cfgmap, H)
        # distribution readout: the states the VQE actually charged for, weighted by the
        # Boltzmann-free CVaR tail it optimised is not available post hoc, so use the
        # final circuit's own sampled frequencies over the visited (cached) states.
        seen = sorted(ham._cache.keys())
        Xw = None
        if seen:
            idx = np.array([int(s, 2) for s in seen], int)
            en = np.array([ham._cache[s] for s in seen], float)
            # CVaR tail of what the search saw, uniformly averaged -- the honest
            # "distribution rather than a point" readout of a budgeted sampler
            kk = max(1, int(round(alpha * len(en))))
            tail = idx[np.argsort(en)[:kk]]
            cf = np.stack([cfgmap[i] if cfgmap is not None else R.state_to_cfg(fam, inst, int(i), cfgmap, H)
                           for i in tail if (cfgmap is not None or R.state_to_cfg(fam, inst, int(i), cfgmap, H) is not None)])
            Xs = V.assemble_batch(inst, cf)
            Xw = I.superpose_batch(Xs, Xs[0]).mean(0)

        arms = {"anchor": inst["anchor"],
                "exact": V.assemble(inst, exh["cfg"]),
                "greedy_ls": V.assemble(inst, gls["cfg"]),
                "vqe": V.assemble(inst, c_vqe)}
        if Xw is not None:
            arms["vqe_cvar_tail_avg"] = Xw
        row = dict(pdb=t["pdb"], n=t["n"], fold=t["fold"], raw={}, proj={},
                   vqe_best_energy=float(res["best_seen_energy"]),
                   opt=float(Eall.min()),
                   vqe_entropy_bits=float(res["distribution_entropy_bits"]))
        for nm, X in arms.items():
            row["raw"][nm] = float(I.ca_rmsd(X, nat))
            if nm == "anchor":
                row["proj"][nm] = row["raw"][nm]      # already on the manifold
            else:
                pr = I.project(X, t["seq"], t["fold"])
                row["proj"][nm] = float(I.ca_rmsd(pr["fit_ca"], nat))
        rows.append(row)
        if q % 5 == 0:
            print(f"  [e2e] {q+1}/{len(tg)} {t['pdb']} {time.time()-t00:.0f}s "
                  f"free={I.free_gb():.2f}", flush=True)
            I.write(out or f"vq_e2e_{fam}", dict(family=fam, budget=budget, rows=rows))
    I.write(out or f"vq_e2e_{fam}", dict(family=fam, budget=budget, rows=rows))
    return rows


if __name__ == "__main__":
    fam = sys.argv[1] if len(sys.argv) > 1 else "A2_64"
    bud = int(sys.argv[2]) if len(sys.argv) > 2 else 3200
    run(fam, bud)
