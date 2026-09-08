"""s19/distobj.py -- DOES THE OBJECTIVE THROW AWAY THE DISTRIBUTION IT WAS GIVEN?

THE OBSERVATION THIS RESTS ON.  `core/predict.py` emits a full 17-bin distribution for every
pair.  The SELECTION path consumes all of it through a Bayes-risk lookup (`Distogram._risk`).
The OBJECTIVE path does not -- `s15/distcal.gather` keeps two numbers:

    "dhat": dg["expected"]      the first moment
    "sd":   dg["sd"]            the second moment

and every objective in Sprints 17 and 18 minimised `((d - dhat)/sd)^2`.  Measured on 12 targets,
the model puts only **0.641** of its probability mass within 1 A of the mean it reports, and
**21.7%** of pairs carry two or more modes.

    On a bimodal pair the expectation lands BETWEEN the modes -- at a distance the model itself
    assigns little mass to, and the geometry may never realise.

WHY THAT WOULD EXPLAIN SPRINT 18.  The programme's most robust finding is that the distogram's
errors are worse than RANDOM errors of the same magnitude (permuting residuals: 3.610 -> 2.560).
A permuted residual is incoherent and lands on achievable geometry on average.  A moment-collapse
error is COHERENT: when a region flips between two conformer families, every pair spanning it
moves the same way, so the mean is wrong in a correlated, structurally-impossible direction.

    THE QUESTION: if the objective consumes the WHOLE predicted distribution instead of its
    first two moments, does the refinement stop being harmful?

ARMS.  All native-free.  The native is read only to score.

    avg          coordinate average of the shipped top-75          the start, ~3.048 A
    moment       ((d - E[d])/sd)^2                                 THE DEPLOYED OBJECTIVE, 3.610
    nll          -log p(d) under the predicted distribution        the full distribution
    risk         Bayes risk, E_p|d - D|                            what SELECTION already uses
    mode         ((d - argmax_p)/sd)^2                             the mean replaced by the MODE
    CONTROLS
    unimodal     `moment` restricted to pairs the model calls unimodal, others dropped
    modeshuf     `mode` with the mode offsets permuted across pairs   zero-information
    matched      a random torsion move of the same displacement       matched-magnitude

WHY `mode` IS THE SHARPEST ARM.  It changes ONE thing -- the target value -- keeping the
functional form, the weights and the optimiser identical.  If the moment-collapse hypothesis is
right, `mode` beats `moment` and the gap is carried by the multimodal pairs.  If `mode` and
`moment` are the same, the hypothesis is dead regardless of what `nll` does, because `nll` also
changes the loss shape and would be confounded.

PRE-REGISTRATION, written before the run.

  HYPOTHESIS.  A distribution-aware objective beats the moment-collapsed one at n = 126, and the
  advantage is concentrated on targets with a high multimodal fraction.

  PRIMARY OUTCOME.  `mode` - `moment`, paired over targets, fold-clustered bootstrap CI.
  SECONDARY.  `nll` - `moment` and `risk` - `moment`, same test.

  EXPECTED, honestly.  I expect a real but PARTIAL effect, and I expect it NOT to reach 3.048.
  Sprint 18 showed the errors are harmful because of their coherent placement; fixing the target
  value on 22% of pairs removes one source of that coherence but not the retrieval-driven and
  sequence-limited sources.  A result near 3.4-3.5 would support the mechanism without rescuing
  refinement.  If `mode` LOSES, the hypothesis is refuted and I close it the same day.

  FALSIFIER.  `mode` - `moment` has a CI spanning zero, OR the effect is not concentrated on the
  multimodal targets (which would mean any gain is not the stated mechanism).

  THE TRAP I AM WATCHING FOR.  `sd` is computed FROM the distribution, so on a bimodal pair it is
  large -- the deployed weighting already down-weights exactly the pairs this hypothesis is about.
  That is a competing explanation for why the effect might be small, and `unimodal` prices it: if
  dropping the multimodal pairs entirely is as good as fixing them, then the weighting had already
  neutralised them and there was never anything to recover.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from core import project as pj              # noqa: E402
from core.predict import BIN_EDGES, CENTRES  # noqa: E402
from s12 import instrument as I             # noqa: E402
from s14 import avgspace as AV              # noqa: E402
from s15 import align_lib as A              # noqa: E402
from s15 import distcal as C                # noqa: E402
from s15 import seed as SD                  # noqa: E402

#: a pair is MULTIMODAL if >=2 bins above this mass are local maxima
MODE_MASS = 0.02

#: BIN WIDTHS.  `BIN_EDGES` is non-uniform (0.5 A below 8 A, up to 4 A above 19), and the first
#: and last bins are unbounded -- their nominal widths are chosen to match the offsets `CENTRES`
#: already uses (`BIN_EDGES[0] - 0.5` and `BIN_EDGES[-1] + 2.0`).
#:
#: THIS CONSTANT EXISTS BECAUSE OF TWO DEFECTS AGENT D FOUND IN THIS FILE.  A bin's MASS is its
#: density times its width, so on a non-uniform grid `argmax(mass)` is not the mode and
#: `-log(mass)` is not a negative log likelihood.  Both were wrong here; see `modality` and
#: `_target_terms`.
WIDTHS = np.concatenate([[1.0], np.diff(BIN_EDGES), [4.0]])


def modality(prob):
    """(peaks, mode_centre, mass_near_mean) for each pair.  Native-free.

    The mode is the DENSITY argmax, not the mass argmax.  Measured by AGENT D at n = 126 over
    8,549 pairs, the two disagree on **14.0%** of pairs and the mass mode sits **+2.72 A further
    out** where they do, a **+0.380 A** mean outward shift overall.  The distogram already
    over-predicts distance (+0.509 A global, +1.492 A at separation 11-15) and the coordinate
    average contracts the backbone 22-25%, so a mass-mode arm carries an expansion confound in
    exactly the direction the shipped separation debias exists to remove.

    Peaks are counted over ALL bins, padding the ends, so a mode in the first or last bin is not
    silently invisible (the previous interior-only comparison could not see them).
    """
    P = np.asarray(prob, float)
    e = (P * CENTRES).sum(1)
    dens = P / WIDTHS[None]
    Pc = P.copy()
    Pc[Pc < MODE_MASS] = 0.0
    pad = np.pad(Pc, ((0, 0), (1, 1)), constant_values=-1.0)
    pk = ((pad[:, 1:-1] > pad[:, :-2]) & (pad[:, 1:-1] >= pad[:, 2:])).sum(1)
    mode = CENTRES[np.argmax(dens, axis=1)]
    near = (P * (np.abs(CENTRES[None, :] - e[:, None]) <= 1.0)).sum(1)
    return pk, mode, near


def _fit(kind, prob, dhat, sd, i, j, phi0, psi0, mask=None, weighted=True, maxiter=400):
    """L-BFGS on one of the objective forms.  The `moment` branch is the deployed one.

    `weighted` applies the deployed per-pair weight `1/sd^2` to the loss.  **AGENT D found that
    the previous version of this docstring was false**: it claimed every branch carried the same
    weights, but `nll` and `risk` used no weight at all, so those arms changed the target AND the
    weighting at once and were uninterpretable.  Sprint 18 priced the weighting at
    +0.153 [+0.066, +0.247], so that confound was the same size as the effect being looked for.

    `nllw` / `riskw` are the matched-weight arms; `nll` / `risk` are kept unweighted and are
    labelled as changing both, because the likelihood's own spread term IS a weighting and
    re-weighting it double-counts.  `mode` was never affected -- it changes the target alone,
    which is why it was made the primary.
    """
    from scipy.optimize import minimize

    n = len(phi0)
    P = np.asarray(prob, float)
    inv = 1.0 / np.asarray(sd, float)
    if mask is None:
        mask = np.ones(len(dhat), bool)
    m = np.asarray(mask, bool)
    inv_m = inv[m]
    #: log DENSITY, not log mass.  `-log p_b = -log(density_b) - log(w_b)`, and on this
    #: non-uniform grid `-log w_b` spans 2.079 nats -- a spurious reward for landing in a wide
    #: long-distance bin, pulling outward in exactly the direction the separation debias removes.
    logP = np.log(np.maximum(P[m] / WIDTHS[None], 1e-9))
    wt = (inv_m * inv_m) if weighted else np.ones(int(m.sum()))

    def _target_terms(d):
        """value and d/dd for the chosen loss, evaluated on the masked pairs."""
        if kind in ("moment", "mode", "unimodal"):
            r = (d - dhat[m]) * inv_m
            return (r * r), (2.0 * r * inv_m)
        #: bin index of each distance, and a soft neighbour interpolation so the loss is
        #: differentiable rather than piecewise constant.
        k = np.clip(np.searchsorted(CENTRES, d) , 1, len(CENTRES) - 1)
        c0, c1 = CENTRES[k - 1], CENTRES[k]
        w = np.clip((d - c0) / np.maximum(c1 - c0, 1e-9), 0.0, 1.0)
        dw = 1.0 / np.maximum(c1 - c0, 1e-9)
        rows = np.arange(len(d))
        if kind in ("nll", "nllw"):
            v = -((1 - w) * logP[rows, k - 1] + w * logP[rows, k])
            g = -(logP[rows, k] - logP[rows, k - 1]) * dw
            return v * wt, g * wt
        if kind in ("risk", "riskw"):
            #: E_p |d - D| , exact under the binned distribution
            dif = d[:, None] - CENTRES[None, :]
            v = (P[m] * np.abs(dif)).sum(1)
            g = (P[m] * np.sign(dif)).sum(1)
            return v * wt, g * wt
        raise ValueError(kind)

    ii, jj = i[m], j[m]

    def fg(x):
        phi, psi = x[:n], x[n:]
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        rv = CA[ii] - CA[jj]
        d = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
        v, gd = _target_terms(d)
        f = float(v.sum())
        coef = (gd / d)[:, None] * rv
        gCA = np.zeros_like(CA)
        np.add.at(gCA, ii, coef)
        np.add.at(gCA, jj, -coef)
        return f, pj._torsion_grad(G, CA, gCA)

    x0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
    res = minimize(fg, x0, jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    return res.x[:n], res.x[n:], float(res.fun)


def run(targets=None):
    #: the debias is ALWAYS fitted on the full target list -- fitting it on a subset silently
    #: changes `dhat` (the recorded Sprint-18 subset trap).
    full = I.targets()
    tg = targets if targets is not None else full
    rows, t0 = [], time.time()
    data = C.gather(full)
    pdbs = [t["pdb"] for t in full]
    deb = {}
    for f in sorted({int(t["fold"]) for t in full}):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        d = data[pdb]
        i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
        bias = deb[fold](d["sep"])
        dhat = np.maximum(d["dhat"] - bias, 2.0)
        rng = SD.stable_rng(pdb, "s19distobj")

        dg = I.distogram(pdb, seq, fold)
        prob = np.asarray(dg["prob"], float)
        pk, mode, near = modality(prob)
        multi = pk >= 2
        #: the mode target carries the SAME leave-fold-out debias as the mean, so the two arms
        #: differ in the statistic only and not in the calibration applied to it.
        dmode = np.maximum(mode - bias, 2.0)

        W = np.asarray(AV.top75_windows(pdb)[0], float)
        Pw = I.pairwise_rmsd(W)
        avg, _b = I.coordinate_average(W, Pw)
        pr = I.project(np.asarray(avg, float), seq, fold)
        phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)
        th0 = np.concatenate([phi0, psi0])

        e = {"avg": float(I.ca_rmsd(np.asarray(avg, float), nat)),
             "proj": float(I.ca_rmsd(np.asarray(pr["ca"], float), nat)),
             "f_multi": float(multi.mean()), "near_mean": float(near.mean()),
             "n_pairs": int(len(dhat))}

        def _emit(tag, ph, ps):
            e[tag] = float(I.ca_rmsd(I.build_ca(ph, ps), nat))
            e["disp_" + tag] = float(np.linalg.norm(
                A.wrap(np.concatenate([ph, ps]) - th0)))

        for tag, kind, tgt in (("moment", "moment", dhat), ("mode", "mode", dmode)):
            ph, ps, _f = _fit(kind, prob, tgt, sd, i, j, phi0, psi0)
            _emit(tag, ph, ps)
        for tag, wtd in (("nll", False), ("risk", False), ("nllw", True), ("riskw", True)):
            ph, ps, _f = _fit(tag, prob, dhat, sd, i, j, phi0, psi0, weighted=wtd)
            _emit(tag, ph, ps)

        #: CONTROL -- the deployed objective on UNIMODAL pairs only.  Prices the competing
        #: explanation that `sd` had already neutralised the multimodal pairs.
        if (~multi).sum() >= 3:
            ph, ps, _f = _fit("unimodal", prob, dhat, sd, i, j, phi0, psi0, mask=~multi)
            _emit("unimodal", ph, ps)
        else:
            e["unimodal"] = e["moment"]; e["disp_unimodal"] = e["disp_moment"]

        #: CONTROL -- the mode OFFSETS permuted across pairs.  This is NOT a clean
        #: zero-information control and must not be read as one: permuting globally also
        #: ORPHANS each offset from its own pair's weight and separation, and Sprint 18 priced
        #: that operation alone at -0.467 [-0.598, -0.338] (G6a).  AGENT D's scaling argument
        #: says the orphaning tax on offsets of this magnitude could account for the whole
        #: `modeshuf` penalty, so `modeshuf - mode` cannot license "the modes carry positional
        #: information".
        off = dmode - dhat
        dshuf = np.maximum(dhat + off[rng.permutation(len(off))], 2.0)
        ph, ps, _f = _fit("mode", prob, dshuf, sd, i, j, phi0, psi0)
        _emit("modeshuf", ph, ps)

        #: CONTROL -- the SAME permutation restricted to (separation bin x sd quartile) strata,
        #: so every offset lands on a pair of comparable weight and separation.  This is the
        #: decomposition Sprint 18 invented for exactly this confound:
        #:      modeshuf_strat - mode        prices POSITION with alignment held
        #:      modeshuf - modeshuf_strat    prices the ORPHANING
        sep = np.asarray(d["sep"], float)
        sb = np.digitize(sep, [2, 3, 4, 5, 6, 8, 11, 15])
        qb = np.digitize(sd, np.quantile(sd, [0.25, 0.5, 0.75]))
        strat = sb * 10 + qb
        off_s = off.copy()
        for key in np.unique(strat):
            cell = np.flatnonzero(strat == key)
            if cell.size > 1:
                off_s[cell] = off[rng.permutation(cell)]
        dss = np.maximum(dhat + off_s, 2.0)
        ph, ps, _f = _fit("mode", prob, dss, sd, i, j, phi0, psi0)
        _emit("modeshuf_strat", ph, ps)

        #: CONTROL -- matched-magnitude random torsion move
        z = rng.standard_normal(2 * n)
        z *= e["disp_mode"] / max(float(np.linalg.norm(z)), 1e-12)
        th = th0 + z
        e["matched"] = float(I.ca_rmsd(I.build_ca(th[:n], th[n:]), nat))

        rows.append({"pdb": pdb, "n": n, "fold": fold, **e})
        if (c + 1) % 5 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": False},
                      open(os.path.join(RESULTS, "distobj.json"), "w"))
    json.dump({"rows": rows, "complete": len(rows) == len(tg)},
              open(os.path.join(RESULTS, "distobj.json"), "w"))
    report(rows)
    return rows


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "distobj.json")))["rows"]
    rng = SD.stable_rng("distobj", "report")
    g = lambda k: np.array([r[k] for r in rows])            # noqa: E731
    avg, mom = g("avg"), g("moment")
    print(f"\nn = {len(rows)}.  Start = the coordinate average of the shipped top-75 "
          f"({avg.mean():.3f} A).  All arms NATIVE-FREE.")
    print(f"  multimodal pair fraction {g('f_multi').mean():.3f} | "
          f"mass within 1 A of the reported mean {g('near_mean').mean():.3f}\n")
    print(f"  {'arm':<12}{'RMSD':>8}{'median':>9}{'vs moment [95% CI]':>27}{'W/L':>9}"
          f"{'vs avg':>10}")
    print(f"  {'avg (start)':<12}{avg.mean():>8.3f}{np.median(avg):>9.3f}{'--':>27}{'--':>9}"
          f"{'--':>10}")
    for k in ("moment", "mode", "nll", "nllw", "risk", "riskw", "unimodal", "modeshuf", "modeshuf_strat", "matched"):
        if k not in rows[0]:
            continue
        v = g(k)
        mu, lo, hi = _boot(v - mom, rng)
        w = int((v < mom).sum()); l = int((v > mom).sum())
        tag = "  CONTROL" if k in ("unimodal", "modeshuf", "modeshuf_strat", "matched") else ("  confounded" if k in ("nll", "risk") else "")
        print(f"  {k:<12}{v.mean():>8.3f}{np.median(v):>9.3f}   "
              f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]{w:>5}/{l}{v.mean()-avg.mean():>+10.3f}{tag}")

    #: THE PRE-REGISTERED SECONDARY: is any gain concentrated where the mechanism says?
    fm = g("f_multi")
    hi_m = fm >= np.median(fm)
    print(f"\n  PRIMARY: mode - moment  {_boot(g('mode') - mom, rng)[0]:+.3f}")
    for lab, sel in (("high-multimodal half", hi_m), ("low-multimodal half", ~hi_m)):
        mu, lo, hi = _boot((g("mode") - mom)[sel], rng)
        print(f"    {lab:<22} n={int(sel.sum()):>4}  {mu:+.3f} [{lo:+.3f},{hi:+.3f}]")
    print("\nREAD.  `mode` differs from `moment` in the TARGET VALUE ALONE -- same weights, same")
    print("form, same optimiser.  If they tie, the moment-collapse hypothesis is dead whatever")
    print("`nll` does.  `unimodal` prices the competing explanation that 1/sd^2 had already")
    print("neutralised the multimodal pairs; `modeshuf` is the zero-information control.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
