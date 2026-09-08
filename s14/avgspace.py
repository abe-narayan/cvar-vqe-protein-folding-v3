"""SPRINT 14, coordinator -- WHERE you average decides what you get.

Finding C2 left a mechanism as an inference: the retrieval pool's per-residue torsion
distribution is a BETTER torsion predictor than any trained sequence model in this project
(phi 33.6 deg, psi 59.2 deg, against the learned predictor's 36.1 / 62.4) and yet emits a
WORSE structure (4.072 A) than the pipeline consuming the identical windows (3.204 A).

The proposed mechanism is that torsion errors compound along the chain while coordinate
errors cancel in place.  Averaging in torsion space commits every residue's error into the
integration; averaging in coordinate space lets independent errors cancel where they are.

That inference is confounded, because the incumbent also runs a multi-start projection onto
the ideal-geometry manifold and the torsion arm does not.  This isolates the two:

    A  coordinate average -> project        the incumbent route
    B  coordinate average, no projection    averaging space alone
    C  torsion circular mean -> build       torsion space alone
    D  torsion circular mean -> build -> project    torsion space, projection restored

A vs C differs in both.  B vs C isolates the AVERAGING SPACE with no projection anywhere.
C vs D isolates the PROJECTION with the averaging space held fixed.

All four consume exactly the same top-75 windows the shipped filter selects, so the
candidate set is held constant and only the aggregation operator moves.

NATIVE-FREE generation; the native is read only to score.

Run:
    python -m s14.avgspace
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
from s14 import retprior as R              # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(RESULTS, exist_ok=True)


def top75_windows(pdb):
    """The exact window set the shipped distogram filter keeps, and their torsions."""
    u = I.load_univ(pdb)
    p = I.pool_idx(u)
    rec = I.shipped_record(pdb)
    sub = np.asarray(rec["sub"], int)
    idx = p[sub]
    return (np.asarray(u["W"], float)[idx], np.asarray(u["PHI"], float)[idx],
            np.asarray(u["PSI"], float)[idx], u)


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    arms = {k: [] for k in ["A_coord_avg_project", "B_coord_avg_raw",
                            "C_torsion_mean_build", "D_torsion_mean_project"]}
    pdbs, folds = [], []

    for t in tg:
        pdb, seq, fold = t["pdb"], t["seq"], int(t["fold"])
        W, PHI, PSI, u = top75_windows(pdb)
        nat = u["nat_ca"]

        # --- coordinate space
        C, _ = I.coordinate_average(W)
        arms["B_coord_avg_raw"].append(I.ca_rmsd(C, nat))
        pr = I.project(C, seq, fold)
        arms["A_coord_avg_project"].append(I.ca_rmsd(pr["fit_ca"], nat))

        # --- torsion space, identical windows
        phi = R.circ_mean(PHI, axis=0)
        psi = R.circ_mean(PSI, axis=0)
        ca = I.build_ca(phi, psi)
        arms["C_torsion_mean_build"].append(I.ca_rmsd(ca, nat))
        pr2 = I.project(ca, seq, fold)
        arms["D_torsion_mean_project"].append(I.ca_rmsd(pr2["fit_ca"], nat))

        pdbs.append(pdb); folds.append(fold)

    folds = np.asarray(folds, int)
    arms = {k: np.asarray(v, float) for k, v in arms.items()}
    fail = np.isin(pdbs, I.FAIL18)

    out = {"n": len(pdbs), "per_target": {k: dict(zip(pdbs, map(float, v)))
                                          for k, v in arms.items()},
           "summary": {k: {**I.summary(v), "FAIL18": float(v[fail].mean())}
                       for k, v in arms.items()},
           "contrasts": {}}

    def con(a, b, label, why):
        out["contrasts"][label] = {"why": why,
                                   **I.paired(arms[a], arms[b], folds=folds, names=pdbs)}

    con("C_torsion_mean_build", "A_coord_avg_project", "torsion_vs_incumbent",
        "both effects together")
    con("C_torsion_mean_build", "B_coord_avg_raw", "averaging_space_alone",
        "same windows, no projection on either side: isolates WHERE the average is taken")
    con("D_torsion_mean_project", "C_torsion_mean_build", "projection_alone",
        "torsion averaging held fixed: isolates what the projection step is worth")
    con("A_coord_avg_project", "B_coord_avg_raw", "projection_on_coordinate_average",
        "the projection's value in the incumbent route, for comparison")
    con("D_torsion_mean_project", "A_coord_avg_project", "residual_after_matching_projection",
        "both projected: whatever remains is the averaging space alone")

    with open(os.path.join(RESULTS, "avgspace.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_avgspace", out, n_expected=len(pdbs))

    print(f"{'arm':<28}{'mean':>8}{'median':>8}{'<2A':>7}{'FAIL18':>9}")
    for k, v in out["summary"].items():
        print(f"{k:<28}{v['mean']:>8.3f}{v['median']:>8.3f}"
              f"{v['frac_under_2.0']:>7.2f}{v['FAIL18']:>9.3f}")
    print()
    for k, v in out["contrasts"].items():
        print(f"{k:<38}{v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]"
              f"  {v['n_better']}W/{v['n_worse']}L  drop10 {v['drop_top10_mean_diff']:+.3f}")
    return out


if __name__ == "__main__":
    run()
