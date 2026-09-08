"""SPRINT 14, coordinator -- the proposal control the structural Hamiltonian owes.

THE CONFOUND.  `s14/hamil.py` draws half its configurations from the retrieval-conditioned
torsion prior and then scores them with an objective whose w=0 arm IS that prior.  A high
in-decile rank correlation could therefore partly measure the agreement between the proposal
and the scorer rather than any structural skill: configurations the prior likes are
over-represented, so the low-energy decile of the prior term is also the dense region of the
sample, and dense regions of a good prior are near-native for reasons that have nothing to
do with the objective ordering them.

This runs the identical analysis under proposals that break that coupling:

  mixed     half prior, half uniform      the arm used in `hamil.py`
  uniform   uniform over all k^n states   the prior cannot agree with the proposal at all
  prior     entirely from the prior       the confound at its maximum, as an upper bound

If the in-decile correlation survives the uniform proposal, it is a property of the
objective.  If it collapses, it was the proposal, and the `hamil.py` numbers must be
re-reported against the uniform arm.  Either way the answer is reportable; that is the point
of running it before quoting anything.

A second control is included: **a sequence-blind twin**, in which the prior term is replaced
by the class back-off prior (which carries generic Ramachandran but almost no
sequence-specific information).  If the structural Hamiltonian's skill is unchanged, the
skill is generic geometry rather than target-specific information -- the same decomposition
that showed 88% of the Sprint 13 library's value was generic Ramachandran.

Run:
    python -m s14.hamil_control
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

PROPOSALS = ["mixed", "uniform", "prior"]
WEIGHTS = [0.0, 0.5, 1.0]


def propose(T, kind, B, rng):
    if kind == "uniform":
        return rng.integers(0, T.k, size=(B, T.n))
    if kind == "prior":
        return H._sample(T.P, B, rng)
    S1 = H._sample(T.P, B // 2, rng)
    S2 = rng.integers(0, T.k, size=(B - B // 2, T.n))
    return np.vstack([S1, S2])


def analyse(t, kind, blind, n_sample, weights, seed=0):
    T = H.Terms(t["pdb"], t["seq"], t["n"], t["fold"])
    if blind:
        # sequence-blind twin: the class back-off prior, generic Ramachandran only
        P = L.prior(t["pdb"], 4)
        T.P = P
        T.logP = np.log(np.clip(P, 1e-12, None))
    rng = np.random.default_rng(hash((t["pdb"], kind, blind, seed)) % (2 ** 32))
    S = propose(T, kind, n_sample, rng)
    ep, ed = H._z(T.e_prior(S)), H._z(T.e_disto(S))
    rm = T.ORACLE_rmsd(S)
    out = {}
    for w in weights:
        E = (1.0 - w) * ep + w * ed
        order = np.argsort(E)
        dec = order[:max(1, len(E) // 10)]
        out[str(w)] = {"rho_all": H._spearman(E, rm),
                       "rho_low_decile": H._spearman(E[dec], rm[dec]),
                       "argmin_rmsd": float(rm[int(order[0])]),
                       "decile_mean_rmsd": float(rm[dec].mean())}
    out["_sample_mean_rmsd"] = float(rm.mean())
    out["_sample_best_rmsd"] = float(rm.min())
    return out


def run(n_sample=2000, targets=None, weights=WEIGHTS):
    tg = targets if targets is not None else I.targets()
    inc = L.incumbent_rmsd()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    ref = np.asarray([inc[p] for p in pdbs], float)

    res = {}
    for kind in PROPOSALS:
        for blind in (False, True):
            if blind and kind != "mixed":
                continue                     # the blind twin is only run on the used arm
            tag = kind + ("_SEQBLIND" if blind else "")
            rows = [analyse(t, kind, blind, n_sample, weights) for t in tg]
            res[tag] = {"per_target": dict(zip(pdbs, rows)),
                        "sample_mean_rmsd": float(np.mean([r["_sample_mean_rmsd"]
                                                           for r in rows])),
                        "sample_best_rmsd": float(np.mean([r["_sample_best_rmsd"]
                                                           for r in rows])),
                        "arms": {}}
            for w in weights:
                k = str(w)
                arg = np.asarray([r[k]["argmin_rmsd"] for r in rows], float)
                res[tag]["arms"][k] = {
                    "rho_all": float(np.mean([r[k]["rho_all"] for r in rows])),
                    "rho_low_decile": float(np.mean([r[k]["rho_low_decile"]
                                                     for r in rows])),
                    "argmin_mean": float(arg.mean()),
                    "decile_mean": float(np.mean([r[k]["decile_mean_rmsd"] for r in rows])),
                    "vs_incumbent": I.paired(arg, ref, folds=folds, names=pdbs)}

    out = {"n": len(tg), "n_sample": n_sample, "results": res,
           "incumbent": float(ref.mean())}
    with open(os.path.join(RESULTS, "hamil_control.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_hamil_control", out, n_expected=len(tg))

    print(f"incumbent {ref.mean():.3f}\n")
    print(f"{'proposal':<20}{'w':>5}{'rho all':>10}{'rho decile':>12}{'argmin':>9}"
          f"{'decile':>9}{'sampled best':>14}")
    for tag, r in res.items():
        for w in weights:
            a = r["arms"][str(w)]
            print(f"{tag:<20}{w:>5.2f}{a['rho_all']:>10.3f}{a['rho_low_decile']:>12.3f}"
                  f"{a['argmin_mean']:>9.3f}{a['decile_mean']:>9.3f}"
                  f"{r['sample_best_rmsd']:>14.3f}")
    return out


if __name__ == "__main__":
    run()
