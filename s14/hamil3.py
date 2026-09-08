"""SPRINT 14, coordinator -- testing the ENER agent's strongest open lead at scale.

THE LEAD.  The ENER agent's sweep of every individual Legacy component left exactly one
survivor as a HYPOTHESIS: `leg_torsion` selects 2.954 A, -0.835 [-1.607,-0.105] against its
comparator, 6W/3L, and it clears the constant-alpha-helix baseline and survives the helix
controls.  If that held it would be the best native-free selection result in the project.

But it was measured on the **nine fully enumerated targets only**, and **two of those nine
carry 68% of the effect**.  With n=9 and that concentration the null cannot discriminate,
which is why they filed it as a hypothesis rather than a finding.  Their own note says the
helix channel is dead on that ensemble (rho +0.003), so the usual confound is absent -- but
absence of the usual confound on nine targets is not evidence.

This tests it on all 126 targets, which is fourteen times the sample, through the same
native-free machinery as `s14/hamil.py` so the numbers compose with the rest of the sprint.
Three questions, in order:

  1. Does `leg_torsion` alone still select well at n=126, and does the concentration survive?
  2. Does adding it to the structural Hamiltonian (retrieval prior + distogram) help, or is
     it redundant with the prior?  Both are local backbone-conformation terms, so redundancy
     is the null hypothesis and complementarity would be the surprise.
  3. Under a STRICTLY UNIFORM proposal as well as the mixed one, because the C7 correction
     showed that a term agreeing with the proposal inflates its own in-decile correlation.

Reported on all four axes the corrected brief now requires -- global rho, in-decile rho,
argmin RMSD, decile-mean RMSD -- because in-decile alone is non-monotone in objective
quality and ambiguous between objectives differing fourfold in real skill.

NATIVE-FREE generation; the native is read only to score.

Run:
    python -m s14.hamil3
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402
from s14 import hamil as H                 # noqa: E402
from s14 import ladder as L                # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(RESULTS, exist_ok=True)

N_SAMPLE = 3000

# (name, w_prior, w_disto, w_torsion) on standardised terms
ARMS = [
    ("prior_only",            1.0, 0.0, 0.0),
    ("disto_only",            0.0, 1.0, 0.0),
    ("legtorsion_only",       0.0, 0.0, 1.0),
    ("prior+disto",           0.75, 0.25, 0.0),
    ("legtorsion+disto",      0.0, 0.25, 0.75),
    ("all_three",             0.50, 0.25, 0.25),
    ("legtorsion+prior",      0.50, 0.0, 0.50),
]


def analyse(t, proposal, n_sample=N_SAMPLE, seed=0):
    from s13.qarch_lib import legacy_components
    T = H.Terms(t["pdb"], t["seq"], t["n"], t["fold"])
    rng = np.random.default_rng(hash((t["pdb"], proposal, seed)) % (2 ** 32))

    if proposal == "uniform":
        S = rng.integers(0, T.k, size=(n_sample, T.n))
    else:
        S1 = H._sample(T.P, n_sample // 2, rng)
        S2 = rng.integers(0, T.k, size=(n_sample - n_sample // 2, T.n))
        S = np.vstack([S1, S2])

    ep = H._z(T.e_prior(S))
    ed = H._z(T.e_disto(S))
    comp = legacy_components(T.sp, S)
    et_ = H._z(np.asarray(comp["torsion"], float))
    rm = T.ORACLE_rmsd(S)                                   # ORACLE, post hoc only

    out = {"pdb": t["pdb"], "fold": int(t["fold"]),
           "oracle_best_sampled": float(rm.min()), "arms": {}}
    for name, wp, wd, wt in ARMS:
        E = wp * ep + wd * ed + wt * et_
        order = np.argsort(E)
        dec = order[:max(1, len(E) // 10)]
        out["arms"][name] = {
            "rho_all": H._spearman(E, rm),
            "rho_low_decile": H._spearman(E[dec], rm[dec]),
            "argmin_rmsd": float(rm[int(order[0])]),
            "decile_mean_rmsd": float(rm[dec].mean()),
        }
    return out


def run(proposals=("mixed", "uniform"), targets=None):
    tg = targets if targets is not None else I.targets()
    inc = L.incumbent_rmsd()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    fail = np.isin(pdbs, I.FAIL18)
    ref = np.asarray([inc[p] for p in pdbs], float)

    out = {"n": len(tg), "n_sample": N_SAMPLE, "incumbent": float(ref.mean()),
           "proposals": {}}
    for prop in proposals:
        rows = [analyse(t, prop) for t in tg]
        block = {"per_target": {r["pdb"]: r for r in rows},
                 "oracle_best_sampled": float(np.mean([r["oracle_best_sampled"]
                                                       for r in rows])),
                 "arms": {}}
        for name, *_ in ARMS:
            arg = np.asarray([r["arms"][name]["argmin_rmsd"] for r in rows], float)
            block["arms"][name] = {
                "rho_all": float(np.mean([r["arms"][name]["rho_all"] for r in rows])),
                "rho_low_decile": float(np.mean([r["arms"][name]["rho_low_decile"]
                                                 for r in rows])),
                "argmin_mean": float(arg.mean()),
                "argmin_FAIL18": float(arg[fail].mean()),
                "decile_mean": float(np.mean([r["arms"][name]["decile_mean_rmsd"]
                                              for r in rows])),
                "vs_incumbent": I.paired(arg, ref, folds=folds, names=pdbs)}
        out["proposals"][prop] = block

    with open(os.path.join(RESULTS, "hamil3.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_hamil3", out, n_expected=len(tg))

    print(f"incumbent {ref.mean():.3f}   (ENER's lead: leg_torsion selected 2.954 A "
          f"on 9 enumerated targets, 2 of 9 carrying 68%)\n")
    for prop, block in out["proposals"].items():
        print(f"--- proposal: {prop}   ORACLE best sampled "
              f"{block['oracle_best_sampled']:.3f}")
        print(f"{'arm':<20}{'rho all':>9}{'rho dec':>9}{'argmin':>9}{'decile':>9}"
              f"{'FAIL18':>9}{'vs incumbent':>26}{'top10 share':>13}")
        for name, *_ in ARMS:
            a = block["arms"][name]
            v = a["vs_incumbent"]
            ts = v.get("top10_share")
            tss = f"{ts:.2f}" if ts is not None else "  -  "
            print(f"{name:<20}{a['rho_all']:>9.3f}{a['rho_low_decile']:>9.3f}"
                  f"{a['argmin_mean']:>9.3f}{a['decile_mean']:>9.3f}"
                  f"{a['argmin_FAIL18']:>9.3f}"
                  f"   {v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]"
                  f"{tss:>13}")
        print()
    return out


if __name__ == "__main__":
    run()
