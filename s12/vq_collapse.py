"""s12 QUANTUM-ROLE -- does CVaR prevent collapse, and is the surviving entropy meaningful?

The one honest argument on the record for CVaR specifically is that it stops the state
collapsing onto the argmin and returns a DISTRIBUTION.  Three things are needed to settle
it, and only the first is usually shown:

1. MEASURED entropy as a function of alpha at T = 0 (no entropy term to do the work).
2. An ANALYTIC bound.  CVaR_alpha is minimised by ANY distribution that puts at least
   `alpha` mass on the argmin -- the upper `1 - alpha` of the distribution is completely
   unconstrained.  So the maximum-entropy exact minimiser has

       H_max(alpha) = -alpha log2 alpha - (1-alpha) log2((1-alpha)/(dim-1))

   CVaR does not *seek* diversity; it is *indifferent* above its tail.  Comparing measured
   entropy with this bound says whether what survives is sought or merely unpenalised.
3. A CONVERGENCE control.  If the entropy at small alpha falls when the optimiser is given
   5-10x the iterations, the breadth was under-optimisation, not a property of CVaR.

Usage: python -m s12.vq_collapse [family] [n_targets]
"""
from __future__ import annotations
import os, sys, json, math, time
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import vq_lib as V
from s12 import vq_quantum as Qm
from s12 import vq_qrun as R
from s12.vq_build import FAMILIES


def h_max_bits(alpha, dim):
    a = float(alpha)
    if a >= 1.0:
        return 0.0
    return float(-a * math.log2(a) - (1 - a) * math.log2((1 - a) / (dim - 1)))


def run(fam="A2_64", ntarg=4, alphas=(0.05, 0.1, 0.25, 0.5, 1.0),
        iters_list=(50, 200, 600), layers=3, seed=0, out=None):
    cfgf = FAMILIES[fam]
    tg = R.stratified(I.targets(), ntarg, cfgf["k"], cfgf["minlen"])
    rows = []
    t00 = time.time()
    for t in tg:
        while I.free_gb() < 1.2:
            print("  waiting on memory", flush=True); time.sleep(20)
        inst = V.cached_instance(t, **cfgf)
        E, N, cfgmap, H = R.energy_of(fam, inst)
        dim = len(E)
        row = dict(pdb=t["pdb"], n_qubits=N, dim=dim, opt=float(E.min()), arms=[])
        for a in alphas:
            for it in iters_list:
                r = Qm.run(E, alpha=a, T=0.0, n=N, layers=layers, iters=it, seed=seed)
                ro = Qm.readouts(E, r["p"])
                p = np.asarray(r["p"], float)
                row["arms"].append(dict(
                    alpha=a, iters=it, entropy_bits=r["entropy_bits"],
                    h_max_bits=h_max_bits(a, dim), max_bits=float(N),
                    cvar=r["cvar"], p_top=ro["p_top"],
                    mass_on_argmin=float(p[int(np.argmin(E))]),
                    gap_argmax=float(ro["argmax_energy"] - float(E.min())),
                    wall=r["wall"]))
        rows.append(row)
        print(f"  [collapse] {t['pdb']} {time.time()-t00:.0f}s free={I.free_gb():.2f}",
              flush=True)
        I.write(out or f"vq_collapse_{fam}", dict(family=fam, rows=rows))
    # pooled
    combos = {}
    for r in rows:
        for a in r["arms"]:
            combos.setdefault((a["alpha"], a["iters"]), []).append(a)
    summ = []
    for (a, it), lst in sorted(combos.items()):
        summ.append(dict(alpha=a, iters=it,
                         entropy_bits=float(np.mean([x["entropy_bits"] for x in lst])),
                         h_max_bits=float(np.mean([x["h_max_bits"] for x in lst])),
                         mass_on_argmin=float(np.mean([x["mass_on_argmin"] for x in lst])),
                         p_top=float(np.mean([x["p_top"] for x in lst])),
                         gap_argmax=float(np.mean([x["gap_argmax"] for x in lst])),
                         wall=float(np.mean([x["wall"] for x in lst]))))
    I.write(out or f"vq_collapse_{fam}", dict(family=fam, summary=summ, rows=rows))
    print(json.dumps(summ, indent=1))
    return summ


if __name__ == "__main__":
    fam = sys.argv[1] if len(sys.argv) > 1 else "A2_64"
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    run(fam, nt)
