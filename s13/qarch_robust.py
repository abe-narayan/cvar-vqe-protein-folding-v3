"""SPRINT 13 QUANTUM-ARCHITECTURE, EXPERIMENT 3b -- AMBER's scale is not standardisable.

`qarch_validity` fits per-residue normalisation constants mu = n*a, sd = sqrt(n)*b on a
native-free reference population.  For Legacy and the prior that works.  For AMBER it does
NOT: the ff14SB/GBn2 single point of an UNMINIMISED discrete-torsion structure is dominated
by van der Waals singularities, so the sample mean over a random population came out at
5.2e14 kcal/mol with sd 2.8e16.  A variance-standardised combination built on those
constants is arithmetically meaningless, and the measured consequence is that
`z_prior + z_amber` emits +0.545 A WORSE than a random draw.

So this file does the honest version: rank-based and median/MAD-based normalisations, which
are scale-free and immune to the singularity, applied to the SAME configurations.  It also
reports the raw scale statistics that force the choice, because "AMBER cannot be combined
by variance standardisation on unminimised torsion configurations" is itself a result.

Output: `s13/results/qarch_robust.json`.

    python -m s13.qarch_robust
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import qarch_lib as Q          # noqa: E402
from s13.qarch_validity import _metrics  # noqa: E402


def robust_z(x):
    x = np.asarray(x, float)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return (x - med) / (1.4826 * mad + 1e-12)


def rank_z(x):
    r = Q._rank(np.asarray(x, float))
    return (r - r.mean()) / (r.std() + 1e-12)


def main():
    pdbs = sorted(f[len("qarch_enum_"):-4] for f in os.listdir(Q.RESULTS)
                  if f.startswith("qarch_enum_") and f.endswith(".npz"))
    rows = []
    for p in pdbs:
        z = np.load(os.path.join(Q.RESULTS, f"qarch_enum_{p}.npz"))
        rmsd = np.asarray(z["rmsd"], float)
        legacy = np.asarray(z["legacy"], float)
        prior = np.asarray(z["prior"], float)
        a_idx = np.asarray(z["amber_idx"], int)
        a_kind = np.asarray(z["amber_kind"], int)
        a_tot = np.asarray(z["amber_total"], float)
        m = a_kind == 0                                    # uniform subsample only
        A, R = a_tot[m], rmsd[a_idx[m]]
        L, P = legacy[a_idx[m]], prior[a_idx[m]]
        rec = {"pdb": p, "n": int(z["n"]),
               "amber_scale": {"mean": float(A.mean()), "sd": float(A.std()),
                               "median": float(np.median(A)),
                               "mad": float(np.median(np.abs(A - np.median(A)))),
                               "p99": float(np.percentile(A, 99)),
                               "max": float(A.max()),
                               "frac_above_1e4_kcal": float((A > 1e4).mean())},
               "legacy_scale": {"mean": float(L.mean()), "sd": float(L.std()),
                                "median": float(np.median(L)),
                                "mad": float(np.median(np.abs(L - np.median(L))))},
               "objectives": {}}
        cand = {
            "amber_rank": rank_z(A), "amber_robustz": robust_z(A),
            "legacy_rank": rank_z(L), "prior_rank": rank_z(P),
            "rank_prior_plus_amber": rank_z(P) + rank_z(A),
            "rank_prior_plus_legacy": rank_z(P) + rank_z(L),
            "rank_prior_plus_legacy_plus_amber": rank_z(P) + rank_z(L) + rank_z(A),
            "robustz_prior_plus_amber": robust_z(P) + robust_z(A),
        }
        for nm, E in cand.items():
            rec["objectives"][nm] = _metrics(E, R, None, "amber_uniform_subsample")
        rows.append(rec)
        print(f"  {p}: amber median {rec['amber_scale']['median']:.3g} "
              f"mean {rec['amber_scale']['mean']:.3g} max {rec['amber_scale']['max']:.3g} "
              f"({rec['amber_scale']['frac_above_1e4_kcal']*100:.0f}% above 1e4 kcal/mol)",
              flush=True)
    Q.write("qarch_robust", {"what": "scale-free normalisations of AMBER and the "
                                     "combinations built on them", "rows": rows})
    # aggregate
    names = list(rows[0]["objectives"])
    print(f"\n{'objective':<36s} {'rho':>7s} {'top1':>7s} {'pool':>7s} {'delta':>8s} {'W/L':>6s}")
    for nm in names:
        rho = np.array([r["objectives"][nm]["spearman_E_rmsd"] for r in rows])
        t1 = np.array([r["objectives"][nm]["top1_rmsd"] for r in rows])
        pm = np.array([r["objectives"][nm]["rmsd_mean"] for r in rows])
        d = t1 - pm
        print(f"{nm:<36s} {np.nanmean(rho):+7.3f} {t1.mean():7.3f} {pm.mean():7.3f} "
              f"{d.mean():+8.3f} {int((d<0).sum())}/{int((d>0).sum())}")
    return rows


if __name__ == "__main__":
    main()
