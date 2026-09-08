"""SPRINT 15, coordinator -- fit the FULL predicted distance distribution, not its mean.

WHAT `s15/distgeo.py` LEFT ON THE TABLE.

The weighted least-squares fit uses `dhat` and `sd` — the distogram's mean and spread. But
`I.distogram` returns `prob`, a full **17-bin distribution over distance for every pair**.
Compressing that to a mean and a standard deviation throws away everything the distribution
knows about **asymmetry and multimodality**, and a distance distribution over a flexible
peptide is very often neither symmetric nor unimodal: a pair is frequently either in contact
or not, with little mass between.

That compression is also exactly where the measured defects live. `s15/distacc.py` found the
mean is biased by **+0.113 A at separation 2-3 rising to +1.492 A at 11-15**, and the sd is
**2.633x over-confident**. A mean-and-sd summary of a bimodal distribution is biased toward
the gap between the modes and reports a spread that describes neither. **The full
distribution may be well calibrated even where its first two moments are not** — that is the
hypothesis this module tests, and it is cheap.

THE OBJECTIVE. Maximum likelihood under the distogram's own predicted distribution:

    f(phi, psi)  =  - sum_{|i-j|>=2}  log P_ij( d_ij(phi, psi) )

This is strictly more informative than least squares, which is its Gaussian special case. It
is also the *natural* objective: the distogram is a probabilistic model, and this asks for the
conformation it considers most probable, rather than the conformation closest to its means.

IMPLEMENTATION. `prob` is (npairs, 17) on `centres`. A cubic-smoothed log-probability is
evaluated by linear interpolation in d with an analytic derivative, floored to keep the
objective finite outside the support. The gradient chains through
`d(log P)/dd * dd/dCA * dCA/d(phi,psi)`, reusing the exact O(n) torsion gradient.

ARMS, all native-free except where marked:

    ls_mean_sd          the s15/distgeo predictive arm, for reference
    ml_full             maximum likelihood under the full 17-bin distribution
    ml_full_debias      the same, on separation-debiased bin centres (leave-fold-out)
    ml_plus_prior       ML restraints + the retrieval torsion prior as a regulariser
                        -- this is the Family A Hamiltonian in its classical limit
    ORACLE_ml_true      ML against a distribution sharply peaked on the TRUE distances,
                        a consistency check that the ML machinery can reach the ceiling

The `ml_plus_prior` arm matters beyond its number: it is the first objective in this project
that combines a **generative distance restraint** with the **best measured torsion channel**
(retrieval prior, phi 33.6 / psi 59.2 deg), and it is exactly the diagonal Hamiltonian a
Family A VQE would encode. Its classical limit must be measured before any quantum arm is
built, or the quantum arm has no control.

Run:
    python -m s15.distml
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

from scipy.optimize import minimize          # noqa: E402

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD             # noqa: E402
from core import project as pj               # noqa: E402
from s15 import distgeo as D                 # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)

LOGP_FLOOR = -12.0          # keeps the objective finite far outside the predicted support


class LogPTable:
    """Piecewise-linear log P(d) per pair on a NON-UNIFORM grid, with an analytic derivative.

    `prob` (P, B) on `centres` (B,).  Outside the support the log-probability is clamped to
    LOGP_FLOOR with zero gradient, so a fit that starts far away is pulled by the restraints
    that ARE informative rather than by numerical noise from the ones that are not.

    CORRECTED 2026-09-05.  The first version of this class computed the bin index as
    `floor((d - c[0]) / (c[1] - c[0]))`, i.e. it assumed a UNIFORM grid.  The distogram's
    centres are not uniform -- they run 0.5 A apart at short range and 4.0 A apart at long
    range -- so every query above about 5 A was read from the wrong bin, badly (a query at
    8.50 A read the bin centred on 7.25 A; a query at 15.00 A read the bin centred on
    17.50 A).  Every `ml_*` number produced before this date is void; see
    `s15/coord_FINDINGS.md` K4-RETRACTED.  The lookup is now a `searchsorted` on the real
    centres with per-interval widths, which is correct for any monotone grid.
    """

    def __init__(self, prob, centres, floor=LOGP_FLOOR):
        p = np.asarray(prob, float)
        p = np.maximum(p, 1e-12)
        p = p / p.sum(1, keepdims=True)
        self.lp = np.maximum(np.log(p), floor)          # (P, B)
        self.c = np.asarray(centres, float)             # (B,) strictly increasing
        assert np.all(np.diff(self.c) > 0), "centres must be strictly increasing"
        self.B = len(self.c)
        self.wid = np.diff(self.c)                      # (B-1,) per-interval width
        # forward differences of the piecewise-linear interpolant, per interval
        self.dlp = np.zeros_like(self.lp)
        self.dlp[:, :-1] = (self.lp[:, 1:] - self.lp[:, :-1]) / self.wid[None, :]

    def _locate(self, d):
        """Interval index and fractional position, correct for a non-uniform grid."""
        k = np.clip(np.searchsorted(self.c, d, side="right") - 1, 0, self.B - 2)
        t = (d - self.c[k]) / self.wid[k]
        return k, t

    def __call__(self, d):
        """Returns (logP, dlogP/dd).  `d` may be (P,) or (..., P)."""
        d = np.asarray(d, float)
        k, t = self._locate(d)
        cols = np.arange(d.shape[-1])
        lo = self.lp[cols, k]
        g = self.dlp[cols, k]
        val = lo + g * (t * self.wid[k])
        inside = (d >= self.c[0]) & (d <= self.c[-1])
        return np.where(inside, val, LOGP_FLOOR), np.where(inside, g, 0.0)

    def rowsum(self, D):
        """(-log-likelihood is summed by the caller) sum_p logP_p(D[..., p]) over pairs.

        Vectorised over leading axes, for scoring a whole enumerated configuration space.
        """
        val, _ = self(np.asarray(D, float))
        return val.sum(-1)


class ShiftedLogP:
    """`P_corrected(d) = P(d + bias)` -- a per-pair additive correction of the SUPPORT.

    The distogram over-predicts distance by a separation-dependent bias, so shifting the
    query is exactly the debiasing.  The derivative is unchanged by a constant shift.
    """

    def __init__(self, table, bias):
        self.t = table
        self.bias = np.asarray(bias, float)

    def __call__(self, d):
        return self.t(d + self.bias)

    def rowsum(self, D):
        val, _ = self(np.asarray(D, float))
        return val.sum(-1)


def fit_ml(table, i, j, phi0, psi0, prior=None, lam=0.0, maxiter=400):
    """Maximise sum log P_ij(d_ij), optionally with a per-residue torsion log-prior."""
    n = len(phi0)

    def fg(x):
        phi = x[:n]; psi = x[n:]
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        r = CA[i] - CA[j]
        d = np.maximum(np.sqrt((r * r).sum(1)), 1e-9)
        lp, dlp = table(d)
        f = -float(lp.sum())
        coef = (-dlp / d)[:, None] * r
        gCA = np.zeros_like(CA)
        np.add.at(gCA, i, coef)
        np.add.at(gCA, j, -coef)
        gt = pj._torsion_grad(G, CA, gCA)
        if prior is not None and lam > 0.0:
            fp, gp = prior(phi, psi)
            f += lam * fp
            gt = gt + lam * gp
        return f, gt

    x0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
    res = minimize(fg, x0, jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    return res.x[:n], res.x[n:], float(res.fun)


def vonmises_prior(mu_phi, mu_psi, kappa):
    """Smooth circular pull toward the retrieval-conditioned torsion means.

    -sum kappa*cos(theta - mu) is the von Mises log-density up to a constant, so this is a
    genuine probabilistic prior rather than an ad-hoc penalty, and its gradient is exact.
    """
    mp = np.asarray(mu_phi, float); ms = np.asarray(mu_psi, float)
    kp = np.asarray(kappa, float)

    def f(phi, psi):
        val = -float((kp * (np.cos(phi - mp) + np.cos(psi - ms))).sum())
        g = np.concatenate([kp * np.sin(phi - mp), kp * np.sin(psi - ms)])
        return val, g
    return f


def run(targets=None, n_start=6, lam=1.0, kappa=2.0):
    from s14 import retprior as R
    from s14 import ladder as L
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    fail = np.isin(pdbs, I.FAIL18)

    # leave-fold-out separation bias, reused from s15.distcal's definition
    from s15 import distcal as C
    data = C.gather(tg)
    corr = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        corr[f], _ = C.fit_correction(data, train, "sep")

    arms = ["ls_mean_sd", "ml_full", "ml_full_debias", "ml_plus_prior", "ORACLE_ml_true"]
    res = {a: [] for a in arms}
    path = os.path.join(RESULTS, "distml.json")

    for c, t in enumerate(tg):
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        d = data[pdb]
        nat = d["nat"]; i, j = d["i"], d["j"]
        dg = I.distogram(pdb, seq, fold)
        prob = np.asarray(dg["prob"], float)
        centres = np.asarray(dg["centres"], float)
        rng = SD.stable_rng(pdb)
        S = D.starts(pdb, seq, n, fold, n_start, rng)

        # retrieval-conditioned circular means, for the prior arm
        PHI, PSI, _ = R.windows(pdb, "top75")
        mu_phi = R.circ_mean(PHI, axis=0); mu_psi = R.circ_mean(PSI, axis=0)
        pri = vonmises_prior(mu_phi, mu_psi, np.full(n, kappa))

        # a distribution sharply peaked on the true distances (ORACLE consistency check)
        tp = np.exp(-0.5 * ((centres[None, :] - d["dtrue"][:, None]) / 0.6) ** 2)

        #: the distogram over-predicts by bias(sep), so the corrected law is
        #: P_corrected(d) = P(d + bias).  `corr` returns (x - bias), so evaluating it at
        #: zero yields -bias, and querying the raw table at (d + bias) is the correction.
        neg_bias = corr[fold](np.zeros_like(d["sep"]), d["sep"])      # = -bias(sep)
        base_tab = LogPTable(prob, centres)
        tables = {
            "ml_full": base_tab,
            "ml_full_debias": ShiftedLogP(base_tab, -neg_bias),
            "ml_plus_prior": base_tab,
            "ORACLE_ml_true": LogPTable(tp, centres),
        }

        for a in arms:
            best = None
            for phi0, psi0, _tag in S:
                if a == "ls_mean_sd":
                    phi, psi, f, _ = D.fit_distances(
                        d["dhat"], 1.0 / (d["sd"] ** 2), i, j, phi0, psi0)
                else:
                    phi, psi, f = fit_ml(
                        tables[a], i, j, phi0, psi0,
                        prior=(pri if a == "ml_plus_prior" else None),
                        lam=(lam if a == "ml_plus_prior" else 0.0))
                if best is None or f < best[0]:
                    best = (f, float(I.ca_rmsd(I.build_ca(phi, psi), nat)))
            res[a].append(best[1])

        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": {a: res[a] for a in arms}, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    base = np.asarray(res["ls_mean_sd"], float)
    out = {"n": len(pdbs), "incumbent": float(ref.mean()), "lam": lam, "kappa": kappa,
           "arms": {}, "per_target": {a: dict(zip(pdbs, map(float, res[a]))) for a in arms}}
    for a in arms:
        v = np.asarray(res[a], float)
        out["arms"][a] = {**I.summary(v), "median": float(np.median(v)),
                          "FAIL18": float(v[fail].mean()),
                          "frac_under_2_5": float((v < 2.5).mean()),
                          "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs),
                          "vs_ls": I.paired(v, base, folds=folds, names=pdbs)}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_distml", out, n_expected=len(pdbs))

    print(f"\nincumbent {ref.mean():.3f}\n")
    print(f"{'arm':<20}{'mean':>8}{'median':>8}{'<2A':>7}{'<2.5A':>8}{'FAIL18':>9}"
          f"{'vs least-squares':>24}{'vs incumbent':>24}")
    for a in arms:
        s = out["arms"][a]
        vl, vi = s["vs_ls"], s["vs_incumbent"]
        print(f"{a:<20}{s['mean']:>8.3f}{s['median']:>8.3f}{s['frac_under_2.0']:>7.2f}"
              f"{s['frac_under_2_5']:>8.2f}{s['FAIL18']:>9.3f}"
              f"  {vl['mean_diff']:+.3f} [{vl['ci95'][0]:+.3f},{vl['ci95'][1]:+.3f}]"
              f"  {vi['mean_diff']:+.3f} [{vi['ci95'][0]:+.3f},{vi['ci95'][1]:+.3f}]")
    return out


if __name__ == "__main__":
    run()
