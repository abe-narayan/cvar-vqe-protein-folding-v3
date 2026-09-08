"""SPRINT 14, coordinator -- the ensemble bridge, and the classical limit of a VQE state.

THE ARGUMENT.  Two coordinator findings combine into an architecture.

  C2  The retrieval pool's position-specific torsion distribution is the best torsion
      channel in this project (phi 33.6 deg, psi 59.2 deg, beating the trained
      leave-fold-out sequence predictor's 36.1 / 62.4) and still emits a worse structure
      (4.072 A) than the pipeline consuming the identical windows (3.204 A).

  C3  Mean absolute torsion error does not order emitted accuracy.  What a backbone builder
      integrates is the CUMULATIVE error, so committing to a single torsion vector pays a
      compounding cost that averaging in coordinate space does not.

If the mechanism is right then a torsion-space method is penalised for COMMITTING TO ONE
CONFIGURATION, not for being in torsion space.  And a variational quantum state does not
commit to one configuration: it *is* a distribution over torsion configurations.  Measure it
B times, build B structures, and average those in COORDINATE space, and the compounding cost
is paid independently by each sample and then cancels -- while the quantum search still
chooses the distribution.

That is the sprint brief's section 19 (set-level objectives, posterior mean, structural
consensus) arriving as a consequence of a measurement rather than as a menu option, and it
is the natural marriage of the project's two architectures.

WHAT THIS FILE MEASURES.  The CLASSICAL LIMIT of that architecture.  A VQE prepared to
exactly reproduce the leakage-safe prior, with no optimisation at all, is the same thing as
sampling the prior.  So:

    sample B configurations from a native-free torsion prior
        -> build B CA traces
        -> coordinate-average (superpose on the medoid, mean)
        -> project onto the ideal-geometry manifold
        -> score

This sets the bar the VQE must clear.  It is the honest control for any later claim that
VQE/CVaR contributes: if the quantum state merely reproduces the prior, it must land here,
and anything it adds has to be measured ON TOP of this number, not on top of the single
torsion vector at 4.072 A.  Getting this wrong is exactly the failure mode the brief's
section 23 warns about -- letting classical preprocessing steal the credit.

NATIVE-FREE generation throughout; the native is read only to score.

Run:
    python -m s14.ensemble
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
from s14 import ladder as L                # noqa: E402
from s14 import retprior as R              # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(RESULTS, exist_ok=True)

SIZES = [1, 3, 5, 10, 20, 50, 75, 200]


def sample_states(P, B, rng):
    """B configurations from a per-residue (n, k) distribution."""
    cum = np.cumsum(P, axis=1)
    u = rng.random((B, P.shape[0]))
    return (u[:, :, None] > cum[None, :, :-1]).sum(2)


def consensus(CA, project_it, seq, fold):
    """Coordinate-average an ensemble of CA traces, optionally project to ideal geometry."""
    C, _ = I.coordinate_average(CA)
    if not project_it:
        return C
    return I.project(C, seq, fold)["fit_ca"]


def run(sizes=SIZES, n_seed=3, arm="top75", targets=None):
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    fail = np.isin(pdbs, I.FAIL18)

    # the reference columns
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)

    curves = {"discrete_k4_project": {}, "discrete_k4_raw": {},
              "continuous_windows_project": {}}
    single = {"torsion_circmean_build": [], "torsion_circmean_project": []}

    for t in tg:
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        u = I.load_univ(pdb)
        nat = u["nat_ca"]
        sp = L.space(pdb, 4)

        # --- the native-free priors
        P = R.state_prior(pdb, arm=arm, k=4)["P"]           # (n, 4) retrieval-conditioned
        PHI_w, PSI_w, _ = R.windows(pdb, arm)               # the raw window torsions

        # --- the single-vector arms, for reference
        phi_m = R.circ_mean(PHI_w, axis=0); psi_m = R.circ_mean(PSI_w, axis=0)
        ca_m = I.build_ca(phi_m, psi_m)
        single["torsion_circmean_build"].append(I.ca_rmsd(ca_m, nat))
        single["torsion_circmean_project"].append(
            I.ca_rmsd(I.project(ca_m, seq, fold)["fit_ca"], nat))

        rows = np.arange(n)
        for B in sizes:
            v_p, v_r, v_c = [], [], []
            for s in range(n_seed):
                rng = np.random.default_rng(hash((pdb, B, s)) % (2 ** 32))

                # discrete: sample the k=4 state prior the quantum encoding uses
                S = sample_states(P, B, rng)
                CA = I.build_ca(sp.PHI[rows[None, :], S], sp.PSI[rows[None, :], S])
                CA = np.atleast_3d(CA).reshape(B, n, 3)
                v_r.append(I.ca_rmsd(consensus(CA, False, seq, fold), nat))
                v_p.append(I.ca_rmsd(consensus(CA, True, seq, fold), nat))

                # continuous: draw whole retrieved windows' torsions, no discretisation
                pick = rng.integers(0, PHI_w.shape[0], size=B)
                CAc = I.build_ca(PHI_w[pick], PSI_w[pick]).reshape(B, n, 3)
                v_c.append(I.ca_rmsd(consensus(CAc, True, seq, fold), nat))

            curves["discrete_k4_project"].setdefault(B, []).append(float(np.mean(v_p)))
            curves["discrete_k4_raw"].setdefault(B, []).append(float(np.mean(v_r)))
            curves["continuous_windows_project"].setdefault(B, []).append(float(np.mean(v_c)))

    out = {"n": len(pdbs), "arm": arm, "n_seed": n_seed, "sizes": sizes,
           "reference": {"incumbent_synthesis": {**I.summary(ref),
                                                 "FAIL18": float(ref[fail].mean())}},
           "single_vector": {}, "curves": {}, "contrasts": {}}

    for k, v in single.items():
        a = np.asarray(v, float)
        out["single_vector"][k] = {**I.summary(a), "FAIL18": float(a[fail].mean()),
                                   "vs_incumbent": I.paired(a, ref, folds=folds, names=pdbs)}

    for cname, cur in curves.items():
        out["curves"][cname] = {}
        for B, vals in cur.items():
            a = np.asarray(vals, float)
            out["curves"][cname][str(B)] = {
                **I.summary(a), "FAIL18": float(a[fail].mean()),
                "vs_incumbent": I.paired(a, ref, folds=folds, names=pdbs)}

    # the contrast the architecture rests on
    b_best = max(sizes)
    a_ens = np.asarray(curves["discrete_k4_project"][b_best], float)
    a_one = np.asarray(single["torsion_circmean_build"], float)
    out["contrasts"]["ensemble_vs_single_vector"] = {
        "why": ("same native-free torsion information; the only change is whether the "
                "method commits to one configuration or averages many in coordinate space"),
        **I.paired(a_ens, a_one, folds=folds, names=pdbs)}

    with open(os.path.join(RESULTS, "ensemble.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_ensemble", out, n_expected=len(pdbs))

    print(f"incumbent synthesis pipeline            {ref.mean():6.3f}   "
          f"FAIL18 {ref[fail].mean():.3f}")
    for k, v in out["single_vector"].items():
        d = v["vs_incumbent"]
        print(f"{k:<40}{v['mean']:6.3f}   FAIL18 {v['FAIL18']:.3f}   "
              f"vs inc {d['mean_diff']:+.3f}")
    print()
    print(f"{'B':>5}" + "".join(f"{c:>30}" for c in curves))
    for B in sizes:
        line = f"{B:>5}"
        for cname in curves:
            r = out["curves"][cname][str(B)]
            line += f"{r['mean']:>12.3f} ({r['vs_incumbent']['mean_diff']:+.3f})   "
        print(line)
    e = out["contrasts"]["ensemble_vs_single_vector"]
    print(f"\nensemble(B={b_best}) vs single torsion vector: {e['mean_diff']:+.3f} "
          f"[{e['ci95'][0]:+.3f},{e['ci95'][1]:+.3f}]  {e['n_better']}W/{e['n_worse']}L"
          f"  drop10 {e['drop_top10_mean_diff']:+.3f}")
    return out


if __name__ == "__main__":
    run()
