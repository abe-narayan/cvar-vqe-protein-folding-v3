"""SPRINT 20, AGENT D, BLOCK E -- is `H_AMBER` an energy, or an energy composed with a relaxation?

Pre-registered in `s20/PREREG_D.md` Block E, falsifier F-D2 fixed before this ran.

`s20/BRIEF.md` section 4 Q1 says "change only H_Legacy <-> H_AMBER" with everything else held.
`amber_hamiltonian.AmberHamiltonian.energy` runs
    openmm.LocalEnergyMinimizer.minimize(ctx, tol=2.0, maxIterations=50)
under a k=100 kcal/mol/A^2 positional restraint, then sets k_rest=0 and reports the UNRESTRAINED
energy at the MINIMISED coordinates.  `core.quantum.FoldingHamiltonian.energy` evaluates in place.
So H_AMBER = E o Relax_50 and the swap is not a swap of two energies on one landscape.

This measures the operator, on the SAME bitstrings, changing ONLY the iteration cap:

    (a) the energy shift E o Relax_50 - E o Relax_1
    (b) the SPEARMAN rank correlation between the two orderings -- the ordering is what CVaR eats
    (c) how often the cap actually BINDS (Z6/Z6b: report the calls the bound BOUND)
    (d) n_collapsed, the +inf guard's firing count

OpenMM's maxIterations=0 means UNBOUNDED, so 1 is the smallest honest setting; the two arms then
differ ONLY in the cap, which is the matched control.

Nothing here reads a native quantity.  No structure is emitted.
"""
from __future__ import annotations
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402

OUT = os.path.join(HERE, "results", "D_E_AMBEROP")
os.makedirs(OUT, exist_ok=True)
SALT = "s20_D_amberop"
NBITS_DRAW = 64


def run(pdbs=None, nbits=NBITS_DRAW):
    import torsion_lib2 as tl2
    import amber_hamiltonian as AH
    from core import quantum as Q
    from scipy.stats import spearmanr

    tg = {t["pdb"]: t for t in I.targets()}
    pdbs = pdbs or ["1CS9", "1A1P", "1A13"]          # 9, 13, 14 residues -- pinned, shortest first
    rows = []
    t0 = time.time()
    for pdb in pdbs:
        seq = tg[pdb]["seq"]
        tab = tl2.library_for(seq, 8, seq)            # production torsion_window = 8
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        rng = SD.stable_rng(pdb, "amberop", salt=SALT)
        bits = ["".join(str(b) for b in rng.integers(0, 2, rep.n_bits)) for _ in range(nbits)]

        rec = {"pdb": pdb, "n": int(tg[pdb]["n"]), "n_bits": int(rep.n_bits),
               "n_states": int(rep.n_states), "bits_per_residue": int(rep.bits_per_residue)}

        # --- the two AMBER arms: identical object, ONLY the iteration cap differs
        E = {}
        for cap in (1, 50):
            h = AH.AmberHamiltonian(seq, rep, minimization_steps=cap)
            e = np.array([h.energy(b) for b in bits], float)
            E[cap] = e
            rec[f"n_collapsed_{cap}"] = int(h.n_collapsed)
            rec[f"t_minimize_{cap}"] = float(h.t_minimize)
            del h

        # --- did the cap BIND?  A call is bound if the 50-step arm moved measurably past
        #     where a 51st step would have taken it; the practical instrument is whether
        #     doubling the cap still changes the energy.
        h100 = AH.AmberHamiltonian(seq, rep, minimization_steps=100)
        e100 = np.array([h100.energy(b) for b in bits], float)
        rec["n_collapsed_100"] = int(h100.n_collapsed)
        del h100

        # Per-arm definedness: the relaxation is not cosmetic if it is what makes the
        # objective FINITE.  Reported per arm before any joint mask is applied.
        rec["n_finite_relax1"] = int(np.isfinite(E[1]).sum())
        rec["n_finite_relax50"] = int(np.isfinite(E[50]).sum())
        rec["n_finite_relax100"] = int(np.isfinite(e100).sum())
        fin = np.isfinite(E[1]) & np.isfinite(E[50]) & np.isfinite(e100)
        rec["n_finite"] = int(fin.sum())
        rec["n_draw"] = int(len(bits))
        d = E[50][fin] - E[1][fin]
        rec["shift_mean"] = float(d.mean())
        rec["shift_med"] = float(np.median(d))
        rec["shift_sd"] = float(d.std())
        rec["spearman_50_vs_1"] = float(spearmanr(E[50][fin], E[1][fin]).correlation)
        rec["spearman_100_vs_50"] = float(spearmanr(e100[fin], E[50][fin]).correlation)
        # the cap binds on a call if 100 steps still moves it by more than 1e-3 kcal/mol
        moved = np.abs(e100[fin] - E[50][fin]) > 1e-3
        rec["frac_cap_bound"] = float(moved.mean())
        rec["bound_shift_med"] = float(np.median(np.abs(e100[fin] - E[50][fin])[moved])) \
            if moved.any() else 0.0
        # --- the Legacy arm on the SAME bitstrings, for the coupling comparison
        hl = Q.FoldingHamiltonian(seq, rep)
        el = np.array([hl.energy(b) for b in bits], float)
        finl = fin & np.isfinite(el)
        rec["spearman_legacy_vs_amber50"] = float(
            spearmanr(el[finl], E[50][finl]).correlation)
        rec["spearman_legacy_vs_amber1"] = float(
            spearmanr(el[finl], E[1][finl]).correlation)
        rows.append(rec)
        print(f"  {pdb} n={rec['n']} bits={rec['n_bits']} finite={rec['n_finite']}/{rec['n_draw']} "
              f"shift {rec['shift_mean']:+.1f} rho50v1 {rec['spearman_50_vs_1']:+.3f} "
              f"capbound {rec['frac_cap_bound']:.2f} ({time.time()-t0:.0f}s)", flush=True)
        json.dump({"rows": rows, "complete": False}, open(os.path.join(OUT, "amberop.json"), "w"))

    full = (len(pdbs) == 3 and nbits == NBITS_DRAW)
    json.dump({"rows": rows, "complete": bool(full), "n_expected": len(pdbs),
               "n_draw": int(nbits), "is_full_preregistered_config": bool(full)},
              open(os.path.join(OUT, "amberop.json"), "w"))
    # A COMPLETE flag is written ONLY for the full pre-registered configuration.  Writing one
    # after a smoke is precisely the hazard this lane flagged in three other artefacts, and
    # this module did it once before the guard was added.
    if full:
        open(os.path.join(OUT, "COMPLETE"), "w").write(f"targets={len(rows)} draws={nbits}\n")
    else:
        p = os.path.join(OUT, "COMPLETE")
        if os.path.exists(p):
            os.remove(p)
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(OUT, "amberop.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    print("\n=== F-D2: is H_AMBER an energy, or E o Relax_50? ===\n")
    print(f"  {'pdb':<7}{'n':>3}{'bits':>6}{'finite':>8}{'E50-E1 mean':>13}{'rho(50,1)':>11}"
          f"{'rho(100,50)':>13}{'cap bound':>11}{'collapsed@50':>14}")
    for r in rows:
        print(f"  {r['pdb']:<7}{r['n']:>3}{r['n_bits']:>6}{r['n_finite']:>4}/{r['n_draw']:<3}"
              f"{r['shift_mean']:>13.1f}{r['spearman_50_vs_1']:>11.3f}"
              f"{r['spearman_100_vs_50']:>13.3f}{r['frac_cap_bound']:>11.2f}"
              f"{r['n_collapsed_50']:>14}")
    print(f"\n  mean Spearman(E o Relax_50, E o Relax_1) = {g('spearman_50_vs_1').mean():.3f}")
    print(f"  mean fraction of calls on which the 50-iteration cap BINDS = "
          f"{g('frac_cap_bound').mean():.3f}")
    print(f"  mean Spearman(E_legacy, E_amber o Relax_50) = "
          f"{g('spearman_legacy_vs_amber50').mean():+.3f}")
    print(f"  mean Spearman(E_legacy, E_amber o Relax_1)  = "
          f"{g('spearman_legacy_vs_amber1').mean():+.3f}")
    tot = g("n_draw").sum()
    print(f"\n  DEFINEDNESS -- the relaxation is what makes the objective finite:")
    print(f"    finite at Relax_1   {int(g('n_finite_relax1').sum())}/{int(tot)}"
          f"   ({g('n_finite_relax1').sum()/tot:.0%})")
    print(f"    finite at Relax_50  {int(g('n_finite_relax50').sum())}/{int(tot)}"
          f"   ({g('n_finite_relax50').sum()/tot:.0%})")
    print(f"    finite at Relax_100 {int(g('n_finite_relax100').sum())}/{int(tot)}"
          f"   ({g('n_finite_relax100').sum()/tot:.0%})")
    rho = g("spearman_50_vs_1").mean()
    cap = g("frac_cap_bound").mean()
    print()
    if rho >= 0.95 and cap < 0.10:
        print("  F-D2 DOES NOT FIRE. The relaxation is a near-monotone reparameterisation and the")
        print("  cap is rarely reached: Q1 is executable roughly as briefed. Objection WITHDRAWN.")
    else:
        why = []
        if rho < 0.95:
            why.append(f"Spearman {rho:.3f} < 0.95 (the ORDERING changes, and CVaR eats the ordering)")
        if cap >= 0.10:
            why.append(f"the cap binds on {cap:.0%} of calls (>= 10%)")
        print("  F-D2 FIRES: " + "; ".join(why) + ".")
        print("  Q1 as briefed compares E_legacy against E_amber o Relax_50. Any lane running it")
        print("  must add the arm E_amber o Relax_1, or label every landscape conclusion")
        print("  RELAXATION-CONFOUNDED.")
    json.dump({"rows": rows, "mean_spearman_50_vs_1": float(rho),
               "mean_frac_cap_bound": float(cap)},
              open(os.path.join(OUT, "amberop_report.json"), "w"), indent=1)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbs", type=str, default=None)
    ap.add_argument("--nbits", type=int, default=NBITS_DRAW)
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.report:
        report()
    else:
        run(a.pdbs.split(",") if a.pdbs else None, a.nbits)
