"""SPRINT 14, coordinator -- the native-free STRUCTURAL Hamiltonian.

WHY THIS AND NOT ANOTHER ENERGY.  Sprint 13 closed the physical energies as search
objectives: Legacy's rank correlation with accuracy inside its own low-energy decile is
+0.043, raw AMBER's -0.088, and Legacy's certified global optimum on nine fully enumerated
targets is +0.139 A WORSE than random sampling.  Optimising either one harder makes
structures worse.  So a VQE pointed at Legacy or AMBER is optimising a landscape whose
minimum is known to be bad, and no ansatz, encoding or CVaR schedule repairs that.

But "the objective must be a molecular energy" was never a requirement -- it was an
inherited assumption.  The project already owns two native-free channels that are
STRUCTURAL rather than energetic, and they carry complementary information:

  * LOCAL.  The retrieval pool's position-specific torsion distribution (`s14/retprior.py`),
    the best torsion channel measured anywhere in this project: phi 33.6 deg, psi 59.2 deg,
    beating the trained leave-fold-out sequence predictor's 36.1 / 62.4.  It is 1-local in
    the torsion encoding, so it is exactly diagonal and cheap.
  * NON-LOCAL.  The shipped leave-fold-out distogram's Bayes-risk score over CA-CA pair
    distances -- the same object the production filter uses, consumed here as an objective
    over CONFIGURATIONS rather than as a filter over retrieved windows.

Neither reaches 2.0 A alone.  The question this file asks is whether their SUM has a
low-energy region that is structurally accurate -- which is the precondition for any
variational search to be worth running at all.

By the Sprint 13 locality theorem the distogram term is genuinely many-body: `d_ij` depends
on exactly the `j-i-1` residues strictly between i and j, so a separation-8 pair is a
14-qubit interaction at k=4 and no 2-local Ising form of it exists.  That is a real cost and
it is stated here rather than hidden.

NORMALISATION.  The brief requires terms be put on comparable scales before any combination,
and warns that a term dominating numerically is a normalisation artefact and not evidence of
superiority.  Both terms are standardised per target by their own mean and sd over the same
sampled configuration set, so the weight `w` is interpretable and the arms differ only in w.

EVERYTHING HERE IS NATIVE-FREE.  The native is read only by `ORACLE_*` scoring functions,
post hoc, to measure whether the objective's low-energy region is any good.

Run:
    python -m s14.hamil
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

WEIGHTS = [0.0, 0.25, 0.5, 0.75, 1.0]      # w = share on the distogram term
N_SAMPLE = 4000


class Terms:
    """Cached per-target machinery for the native-free structural Hamiltonian."""

    def __init__(self, pdb, seq, n, fold, k=4, arm="top75"):
        self.pdb, self.seq, self.n, self.fold, self.k = pdb, seq, int(n), int(fold), k
        self.sp = L.space(pdb, k)
        self.P = R.state_prior(pdb, arm=arm, k=k)["P"]           # (n, k) native-free
        self.logP = np.log(np.clip(self.P, 1e-12, None))
        self.dg = I.distogram(pdb, seq, fold)                    # leave-fold-out, native-free
        self.i, self.j = I.pair_index(self.n)
        self._rows = np.arange(self.n)

    # ---- term 1: local torsion prior.  1-local, diagonal, exact.
    def e_prior(self, S):
        S = np.atleast_2d(np.asarray(S, int))
        return -self.logP[self._rows[None, :], S].sum(1)

    # ---- term 2: distogram Bayes risk.  Many-body by the locality theorem.
    def e_disto(self, S, chunk=2000):
        S = np.atleast_2d(np.asarray(S, int))
        out = np.empty(len(S))
        for a in range(0, len(S), chunk):
            b = S[a:a + chunk]
            CA = I.build_ca(self.sp.PHI[self._rows[None, :], b],
                            self.sp.PSI[self._rows[None, :], b])
            CA = np.asarray(CA, float).reshape(len(b), self.n, 3)
            D = I.pair_dists(CA, self.i, self.j)
            out[a:a + chunk] = I.shipped_score(self.dg, D)
        return out

    def ORACLE_rmsd(self, S):
        return self.sp.rmsd(S)


def _z(x):
    x = np.asarray(x, float)
    s = x.std()
    return (x - x.mean()) / (s if s > 1e-12 else 1.0)


def _spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / d) if d > 0 else 0.0


def analyse_target(t, n_sample=N_SAMPLE, weights=WEIGHTS, seed=0):
    """Sample the prior, score every arm, and report BOTH axes separately."""
    T = Terms(t["pdb"], t["seq"], t["n"], t["fold"])
    rng = np.random.default_rng(hash((t["pdb"], seed)) % (2 ** 32))

    # native-free proposal distribution: half from the prior, half uniform, so the
    # low-energy decile is not an artefact of the proposal agreeing with one term
    B = n_sample
    S1 = _sample(T.P, B // 2, rng)
    S2 = rng.integers(0, T.k, size=(B - B // 2, T.n))
    S = np.vstack([S1, S2])

    ep, ed = _z(T.e_prior(S)), _z(T.e_disto(S))
    rm = T.ORACLE_rmsd(S)                                  # ORACLE, post hoc only

    res = {"pdb": t["pdb"], "n": T.n, "fold": int(t["fold"]),
           "rmsd_best_sampled": float(rm.min()), "rmsd_mean_sampled": float(rm.mean()),
           "arms": {}}
    for w in weights:
        E = (1.0 - w) * ep + w * ed
        order = np.argsort(E)
        dec = order[:max(1, len(E) // 10)]
        res["arms"][str(w)] = {
            "rho_all": _spearman(E, rm),
            "rho_low_decile": _spearman(E[dec], rm[dec]),
            "argmin_rmsd": float(rm[int(order[0])]),
            "top10_mean_rmsd": float(rm[order[:10]].mean()),
            "decile_mean_rmsd": float(rm[dec].mean()),
            "decile_best_rmsd": float(rm[dec].min()),
        }
    return res


def _sample(P, B, rng):
    cum = np.cumsum(P, axis=1)
    u = rng.random((B, P.shape[0]))
    return (u[:, :, None] > cum[None, :, :-1]).sum(2)


def run(weights=WEIGHTS, n_sample=N_SAMPLE, targets=None):
    tg = targets if targets is not None else I.targets()
    rows = [analyse_target(t, n_sample, weights) for t in tg]

    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([r["fold"] for r in rows], int)
    fail = np.isin(pdbs, I.FAIL18)
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)

    out = {"n": len(rows), "n_sample": n_sample, "weights": weights,
           "per_target": rows, "summary": {},
           "reference": {"incumbent": float(ref.mean()),
                         "best_sampled_ORACLE": float(np.mean([r["rmsd_best_sampled"]
                                                               for r in rows])),
                         "mean_sampled": float(np.mean([r["rmsd_mean_sampled"]
                                                        for r in rows]))}}
    for w in weights:
        k = str(w)
        arg = np.asarray([r["arms"][k]["argmin_rmsd"] for r in rows], float)
        out["summary"][k] = {
            "rho_all": float(np.mean([r["arms"][k]["rho_all"] for r in rows])),
            "rho_low_decile": float(np.mean([r["arms"][k]["rho_low_decile"] for r in rows])),
            "argmin": {**I.summary(arg), "FAIL18": float(arg[fail].mean())},
            "top10_mean": float(np.mean([r["arms"][k]["top10_mean_rmsd"] for r in rows])),
            "decile_mean": float(np.mean([r["arms"][k]["decile_mean_rmsd"] for r in rows])),
            "vs_incumbent": I.paired(arg, ref, folds=folds, names=pdbs),
        }

    with open(os.path.join(RESULTS, "hamil.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_hamil", out, n_expected=len(rows))

    r = out["reference"]
    print(f"incumbent {r['incumbent']:.3f}   ORACLE best of the {n_sample} sampled "
          f"{r['best_sampled_ORACLE']:.3f}   mean sampled {r['mean_sampled']:.3f}\n")
    print(f"{'w(disto)':>9}{'rho all':>10}{'rho decile':>12}{'argmin':>9}{'top10':>9}"
          f"{'decile':>9}{'FAIL18':>9}{'vs incumbent':>26}")
    for w in weights:
        s = out["summary"][str(w)]
        v = s["vs_incumbent"]
        print(f"{w:>9.2f}{s['rho_all']:>10.3f}{s['rho_low_decile']:>12.3f}"
              f"{s['argmin']['mean']:>9.3f}{s['top10_mean']:>9.3f}{s['decile_mean']:>9.3f}"
              f"{s['argmin']['FAIL18']:>9.3f}"
              f"   {v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]")
    return out


if __name__ == "__main__":
    run()
