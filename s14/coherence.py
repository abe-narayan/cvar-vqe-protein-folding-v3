"""SPRINT 14, coordinator -- torsion-error COHERENCE, and why sigma is not the axis.

THE ANOMALY THAT STARTED THIS.  In `s14/retprior.py`, two emitters built from the same
retrieved windows came out backwards:

    L2b top75 circular mean      phi 33.6 deg  psi 59.2 deg   ->  4.072 A
    L2d top75 similarity-weighted phi 38.5 deg  psi 62.4 deg   ->  3.514 A

The arm with 4.9 deg WORSE phi error emits a 0.56 A BETTER structure.  Mean absolute
torsion error therefore does not order emitted accuracy, and the whole restraint surface
-- sigma x coverage -> RMSD, the sprint's most attractive route -- is indexed on sigma.

THE HYPOTHESIS.  Torsion errors compound along a chain.  What a backbone builder actually
integrates is the CUMULATIVE error, so two error fields with identical marginal magnitude
give very different structures depending on whether successive errors reinforce or cancel:

  * i.i.d. errors partially cancel; the CA trace random-walks around the native.
  * COHERENT errors (a constant bias, or positively autocorrelated along the chain) turn
    into systematic curvature and the chain leaves the native monotonically.
  * ANTI-coherent errors cancel almost exactly and cost far less than sigma predicts.

This is the Sprint 12 error-coherence law reappearing in torsion space: at identical
accuracy, coherent mistakes and i.i.d. mistakes have opposite consequences.

WHY IT MATTERS TO THIS SPRINT.  Sprint 12's restraint surface corrupted native torsions
with I.I.D. Gaussian noise.  If real predictors -- retrieval priors, TALOS-class shift
methods, sequence models -- have COHERENT errors, that surface is priced on the wrong
error model, and every sigma threshold derived from it (including the "sigma 12 deg
reaches 1.486 A" figure and the "incumbent equals sigma ~29 deg" equivalence) needs a
coherence coordinate before it can be trusted.

ORACLE DIAGNOSTIC.  Everything here reads native torsions.  It prices an error model; it
is not a method and emits no prediction.

Run:
    python -m s14.coherence
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


def wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


def native_torsions(pdb):
    import peptide_db as pdbm
    p = pdbm.by_pdb(pdb)
    return np.asarray(p.phi, float), np.asarray(p.psi, float)


# --------------------------------------------------------------- part 1: measure real errors
def error_field(emit, pdb, seq, n, fold):
    """Signed torsion error of an emitter against the native, in radians."""
    rng = np.random.default_rng(0)
    phi, psi = emit(pdb, seq, int(n), int(fold), rng)
    nphi, npsi = native_torsions(pdb)
    return wrap(phi - nphi), wrap(psi - npsi)


def coherence_stats(e):
    """Bias, magnitude and along-chain autocorrelation of one error field."""
    e = np.asarray(e, float)
    if len(e) < 3:
        return None
    mae = float(np.rad2deg(np.abs(e)).mean())
    bias = float(np.rad2deg(e.mean()))
    x = e - e.mean()
    denom = float((x * x).sum())
    r1 = float((x[:-1] * x[1:]).sum() / denom) if denom > 1e-12 else 0.0
    # the quantity a builder actually integrates
    cum = float(np.rad2deg(np.abs(np.cumsum(e))).mean())
    return {"mae_deg": mae, "bias_deg": bias, "lag1_autocorr": r1,
            "mean_abs_cumulative_deg": cum,
            "coherence_ratio": cum / max(mae, 1e-9)}


def measure_real_emitters():
    """Coherence of every native-free emitter built so far."""
    emitters = dict(L.EMITTERS)
    emitters.update(R.EMITTERS)
    tg = I.targets()
    rows = {}
    for name, (emit, _) in emitters.items():
        ph, ps, out = [], [], []
        for t in tg:
            ephi, epsi = error_field(emit, t["pdb"], t["seq"], t["n"], t["fold"])
            a = coherence_stats(ephi[1:])
            b = coherence_stats(epsi[:-1])
            if a and b:
                ph.append(a); ps.append(b)
            u = I.load_univ(t["pdb"])
            rng = np.random.default_rng(0)
            p2, s2 = emit(t["pdb"], t["seq"], int(t["n"]), int(t["fold"]), rng)
            out.append(I.ca_rmsd(I.build_ca(p2, s2), u["nat_ca"]))
        rows[name] = {
            "rmsd": float(np.mean(out)),
            "phi": {k: float(np.mean([x[k] for x in ph])) for k in ph[0]},
            "psi": {k: float(np.mean([x[k] for x in ps])) for k in ps[0]},
        }
    return rows


# ------------------------------------------------- part 2: the controlled coherence surface
def synth_error(n, sigma_deg, mode, rng):
    """An error field of prescribed marginal sigma and prescribed coherence.

    Every mode is rescaled to the SAME marginal standard deviation, so the arms differ
    only in the along-chain correlation structure and never in magnitude.
    """
    s = np.deg2rad(sigma_deg)
    if mode == "iid":
        e = rng.normal(size=n)
    elif mode == "bias":                       # perfectly coherent: one constant offset
        e = np.full(n, rng.normal())
    elif mode == "ar_pos":                     # positively autocorrelated, rho = 0.7
        e = np.empty(n); e[0] = rng.normal()
        for i in range(1, n):
            e[i] = 0.7 * e[i - 1] + rng.normal() * np.sqrt(1 - 0.49)
    elif mode == "ar_neg":                     # anti-correlated, rho = -0.7
        e = np.empty(n); e[0] = rng.normal()
        for i in range(1, n):
            e[i] = -0.7 * e[i - 1] + rng.normal() * np.sqrt(1 - 0.49)
    elif mode == "alternating":                # maximal cancellation
        e = rng.normal(size=n) * np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
        e = e - e.mean()
    else:
        raise ValueError(mode)
    sd = e.std()
    if sd > 1e-12:
        e = e / sd
    return e * s


MODES = ["iid", "bias", "ar_pos", "ar_neg", "alternating"]
SIGMAS = [5, 10, 12, 15, 20, 25, 30, 40]


def coherence_surface(sigmas=SIGMAS, modes=MODES, n_rep=12, targets=None):
    """ORACLE.  Corrupt native torsions at matched sigma but different coherence, rebuild,
    and score.  This is the missing coordinate of the Sprint 12 restraint surface."""
    tg = targets if targets is not None else I.targets()
    cells = {}
    for mode in modes:
        for sg in sigmas:
            vals = []
            for t in tg:
                u = I.load_univ(t["pdb"])
                nphi, npsi = native_torsions(t["pdb"])
                n = len(nphi)
                rr = []
                for rep in range(n_rep):
                    rng = np.random.default_rng(hash((t["pdb"], mode, sg, rep)) % (2 ** 32))
                    phi = wrap(nphi + synth_error(n, sg, mode, rng))
                    psi = wrap(npsi + synth_error(n, sg, mode, rng))
                    rr.append(I.ca_rmsd(I.build_ca(phi, psi), u["nat_ca"]))
                vals.append(float(np.mean(rr)))
            cells[f"{mode}_s{sg}"] = {"mode": mode, "sigma": sg,
                                      "mean": float(np.mean(vals)),
                                      "median": float(np.median(vals)),
                                      "frac_under_2": float((np.asarray(vals) < 2.0).mean()),
                                      "per_target": vals}
            print(f"  {mode:<12} sigma {sg:>3} deg -> {np.mean(vals):6.3f} A"
                  f"   (<2 A {np.mean(np.asarray(vals) < 2.0):.2f})")
    return cells


def sigma_equivalent(cells, mode, rmsd_target):
    """Interpolate the sigma at which a given coherence mode reaches a target RMSD."""
    xs = sorted({c["sigma"] for c in cells.values() if c["mode"] == mode})
    ys = [next(c["mean"] for c in cells.values() if c["mode"] == mode and c["sigma"] == x)
          for x in xs]
    if rmsd_target <= ys[0]:
        return float(xs[0])
    if rmsd_target >= ys[-1]:
        return float(xs[-1])
    return float(np.interp(rmsd_target, ys, xs))


def run():
    print("PART 1 -- coherence of the real native-free emitters")
    real = measure_real_emitters()
    print(f"{'emitter':<26}{'RMSD':>7}{'phi MAE':>9}{'phi r1':>8}{'phi bias':>10}"
          f"{'psi MAE':>9}{'psi r1':>8}{'cum/MAE':>9}")
    for k, v in sorted(real.items(), key=lambda kv: kv[1]["rmsd"]):
        print(f"{k:<26}{v['rmsd']:>7.3f}{v['phi']['mae_deg']:>9.1f}"
              f"{v['phi']['lag1_autocorr']:>8.3f}{v['phi']['bias_deg']:>10.1f}"
              f"{v['psi']['mae_deg']:>9.1f}{v['psi']['lag1_autocorr']:>8.3f}"
              f"{v['phi']['coherence_ratio']:>9.2f}")

    print("\nPART 2 -- ORACLE coherence surface, matched sigma, 126 targets")
    cells = coherence_surface()

    inc = 3.2040761603809194
    equiv = {m: sigma_equivalent(cells, m, inc) for m in MODES}
    two = {m: sigma_equivalent(cells, m, 2.0) for m in MODES}
    print("\nsigma equivalent to the incumbent 3.204 A, by coherence mode:")
    for m in MODES:
        print(f"  {m:<12} {equiv[m]:6.1f} deg      sigma to reach 2.0 A: {two[m]:6.1f} deg")

    out = {"real_emitters": real, "surface": cells,
           "sigma_equivalent_to_incumbent": equiv, "sigma_to_reach_2A": two,
           "note": ("Sprint 12's restraint surface used iid corruption only. These arms are "
                    "matched in marginal sigma and differ only in along-chain correlation.")}
    with open(os.path.join(RESULTS, "coherence.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_coherence", out, n_expected=len(MODES) * len(SIGMAS))
    return out


if __name__ == "__main__":
    run()
