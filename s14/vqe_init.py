"""SPRINT 14 / VQE -- PART D, the initialisation question, run standalone.

The brief is emphatic: an initialisation that hands VQE the answer is not a VQE result, and
the only way to see that is to print the initialisation's own RMSD next to the outcome.
This is the reduced, affordable form of `s14/vqe_ansatz.py`'s D4/D5 -- the full module was
stopped when the shared box fell to ~0.07 cores per process and would not have finished.

THE HONEST CONTROL, and it is the point of the module. "After VQE" is not compared against
the initialisation's MEAN.  It is compared against `init_best`: the best of the SAME NUMBER
OF SHOTS drawn from the INITIAL distribution with no optimisation at all.  That is the arm a
VQE has to beat to have contributed anything, because best-of-N sampling is free and is what
the shipped pipeline already does.  DELTA = after - init_best.  Only a negative DELTA is a
contribution from the optimisation.

    python -m s14.vqe_init
"""
from __future__ import annotations

import numpy as np

from s12 import instrument as I
from s14 import vqe_lib as V
from s14 import vqe_run as R


def main(targets=("1CS9", "2MK7", "2P5H", "6EY3"), signals=(0.1, 0.3),
         budget=10240, shots=512, seeds=(0, 1, 2)):
    from s13.qarch_lib import Space, empirical_prior, ORACLE_prior
    print("=" * 104)
    print("PART D5. INITIALISATION: acceleration, or handing over the answer?")
    print("  Every row prints the RMSD of its OWN initialisation before the outcome.")
    print("  DELTA = after VQE - best of the same shots from the INITIAL distribution.")
    print("=" * 104)
    out = {}
    for sig in signals:
        print(f"\n--- signal {sig}, budget {budget} evaluations, "
              f"{len(targets)} targets x {len(seeds)} seeds ---")
        print(f"  {'initialisation':44s} {'INIT mean':>10s} {'init best':>10s} "
              f"{'after VQE':>10s} {'DELTA':>8s} {'certified':>10s}")
        rowacc = {}
        for pdb in targets:
            e = V.Enum(pdb)
            E = V.blend_objective(e.legacy, e.rmsd, sig)
            sp = Space(pdb, k=4)
            arms = {"random (native-free)": ("random", None),
                    "uniform / zero (native-free)": ("zero", None),
                    "empirical prior warm start (NATIVE-FREE)":
                        ("prior", empirical_prior(sp))}
            for q in (0.5, 0.8, 0.95):
                arms[f"ORACLE prior q={q} warm start (ORACLE DIAGNOSTIC)"] = \
                    ("prior", ORACLE_prior(sp, q))
            for name, (mode, P) in arms.items():
                rs = [R.run(E, e.n_qubits, 0.25, budget=budget, shots=shots,
                            ansatz="mps2f", seed=s, rmsd=e.rmsd, init=mode, prior=P,
                            bits_per_res=2, exact_dist=False) for s in seeds]
                rowacc.setdefault(name, {"i0": [], "ib": [], "af": [], "ex": []})
                rowacc[name]["i0"] += [r["init_rmsd_mean"] for r in rs]
                rowacc[name]["ib"] += [r["init_rmsd_best"] for r in rs]
                rowacc[name]["af"] += [r["rmsd_returned"] for r in rs]
                rowacc[name]["ex"].append(float(e.rmsd[int(np.argmin(E))]))
        for name, a in rowacc.items():
            i0, ib, af = np.mean(a["i0"]), np.mean(a["ib"]), np.mean(a["af"])
            out.setdefault(str(sig), {})[name] = {
                "init_mean": float(i0), "init_best": float(ib), "after": float(af),
                "delta": float(af - ib), "certified": float(np.mean(a["ex"])),
                "per_run_after": a["af"], "per_run_init_best": a["ib"]}
            print(f"  {name:44s} {i0:10.3f} {ib:10.3f} {af:10.3f} {af-ib:+8.3f} "
                  f"{np.mean(a['ex']):10.3f}")
        # the paired test that decides it, over runs
        print()
        for name, a in rowacc.items():
            p = I.paired(np.array(a["af"]), np.array(a["ib"]))
            print(f"    {name[:42]:44s} after vs init-best: {p['mean_diff']:+.3f} "
                  f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] W/L "
                  f"{p['n_better']}/{p['n_worse']}")
            out[str(sig)][name]["paired_vs_init_best"] = p
    print("\n  A NEGATIVE mean difference is the only thing that can be called a")
    print("  contribution from the optimisation.  A warm start that lowers 'init best'")
    print("  and leaves DELTA at zero has handed over the answer, not accelerated a search.")
    V.write("vqe_init", out)
    print("\nwritten -> s14/results/vqe_init.json")


if __name__ == "__main__":
    main()
