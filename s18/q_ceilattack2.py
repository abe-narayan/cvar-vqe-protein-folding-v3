"""SPRINT 18 / QUANTUM-ADVERSARIAL -- ATTACK 5 on `objceil`: `shuf_paired` changes TWO things,
and only one of them is the error's assignment to pairs.

THE DEFECT.  `s18/objceil.py` builds its strongest control as

    pi = rng.permutation(len(r))
    pw = dtrue + r[pi] * random_sign
    A.fit(pw, sd[pi], i, j, phi0, psi0)                 <-- NOTE: sd[pi], not sd

so the arm permutes the residual AND passes a PERMUTED WEIGHT VECTOR into the fit.  The
deployed functional is `sum_p (d_p - dhat_p)^2 / sd_p^2`; with `sd[pi]` the arm is no longer
that functional -- it is a DIFFERENT objective, weighted by someone else's confidences.

`shuffled` (global permutation of r, weights NOT permuted) reaches 2.609 A.
`shuf_paired` (the same permutation of r, weights permuted too) reaches 2.072 A.
The entire extra 0.537 A therefore comes from permuting the WEIGHTS, not from anything about
where the errors land.  The coordinator's reading -- "the distogram's errors are harmful
because of their assignment to pairs" -- is confounded with "the distogram's 1/sd^2 confidences
are anti-informative", which is a different claim with a different prescription.

THE ISOLATING CONTROLS, at n = 126, on the identical start and machinery:

    wperm_only    the REAL dhat, with a permuted weight vector.  Residuals untouched, their
                  assignment to pairs untouched.  If this alone beats alpha = 0 by most of
                  the 1.537 A, the effect is a WEIGHTING result and `shuf_paired` should not
                  be read as an error-assignment result.
    wflat         the REAL dhat with UNIFORM weights (every 1/sd^2 replaced by its mean).  The
                  zero-information version of the same question: does the distogram's
                  confidence carry ANY usable weighting information?
    shufr_wkeep   r permuted, weights kept  -- reproduces `shuffled` in this module, so the
                  two modules' numbers are comparable on the same start.
    shufr_wperm   r permuted, weights permuted -- reproduces `shuf_paired`.
    a1_wperm      PERFECT distances with a permuted weight vector.  If the weight permutation
                  is a real gain it should be neutral here (nothing to mis-weight), and if it
                  still helps, the weights are simply bad.

Every arm is an ORACLE DIAGNOSTIC except `wperm_only`/`wflat`, which read no native distance
at all and are therefore NATIVE-FREE -- and that matters: if `wflat` alone improves the
refinement, it is a deployable change, not a ceiling.

RUN:  python -m s18.q_ceilattack2 run
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s14 import avgspace as AV               # noqa: E402
from s15 import align_lib as AL              # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import seed as SD                   # noqa: E402

SALT = "s18qceil2"
OUT = os.path.join(RESULTS, "q_ceilattack2.json")


def _rmsd_of(phi, psi, nat):
    return float(I.ca_rmsd(I.build_ca(phi, psi), nat))


def cell(t, data, deb):
    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    d = data[pdb]
    i, j, sd, nat = d["i"], d["j"], np.asarray(d["sd"], float), np.asarray(d["nat"], float)
    dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
    dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
    r = dhat - dtrue
    rng = SD.stable_rng(pdb, "ceilattack2", salt=SALT)

    W = np.asarray(AV.top75_windows(pdb)[0], float)
    P = I.pairwise_rmsd(W)
    avg, _b = I.coordinate_average(W, P)
    pr = I.project(np.asarray(avg, float), seq, fold)
    phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)

    pi = rng.permutation(len(r))
    sgn = rng.choice([-1.0, 1.0], size=len(r))
    sd_perm = sd[pi]
    sd_flat = np.full_like(sd, float(np.sqrt(1.0 / np.mean(1.0 / sd ** 2))))

    arms = {
        # NATIVE-FREE -- these read no native distance
        "a0_wkeep": (dhat, sd),
        "wperm_only": (dhat, sd_perm),
        "wflat": (dhat, sd_flat),
        # ORACLE DIAGNOSTICS
        "shufr_wkeep": (np.maximum(dtrue + r[pi] * sgn, 2.0), sd),
        "shufr_wperm": (np.maximum(dtrue + r[pi] * sgn, 2.0), sd_perm),
        "shufr_wflat": (np.maximum(dtrue + r[pi] * sgn, 2.0), sd_flat),
        "a1_wkeep": (dtrue, sd),
        "a1_wperm": (dtrue, sd_perm),
    }
    e = {"pdb": pdb, "n": n, "fold": fold,
         "avg": float(I.ca_rmsd(np.asarray(avg, float), nat))}
    for nm, (dd, ww) in arms.items():
        p_, q_, _f = AL.fit(dd, ww, i, j, phi0, psi0)
        e[nm] = _rmsd_of(p_, q_, nat)
    return e


def run(limit=None):
    tg = I.targets()
    if limit:
        tg = tg[:int(limit)]
    data = C.gather(tg)
    pdbs = [t["pdb"] for t in tg]
    deb = {}
    for f in sorted({int(t["fold"]) for t in tg}):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))
    rows = []
    if os.path.exists(OUT):
        try:
            rows = json.load(open(OUT))["rows"]
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        rows.append(cell(t, data, deb))
        if len(rows) % 10 == 0:
            json.dump({"rows": rows, "complete": False}, open(OUT, "w"), default=float)
            print(f"  {len(rows)}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    report(rows)


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rng = SD.stable_rng("ceilattack2", "report", salt=SALT)
    g = lambda k: np.array([r[k] for r in rows])                 # noqa: E731
    P = print
    a0 = g("a0_wkeep")
    P("=" * 100)
    P(f"ATTACK 5 -- `shuf_paired` CHANGES TWO THINGS.   n = {len(rows)}")
    P("=" * 100)
    P(f"  coordinate average (start)   {g('avg').mean():.3f}")
    P(f"  alpha = 0, real weights      {a0.mean():.3f}   (objceil's a0.0 = 3.610)")
    P("")
    P(f"{'arm':<16s}{'reads native?':>15s}{'RMSD':>9s}{'median':>9s}{'vs a0':>26s}{'W/L':>9s}")
    lab = {"a0_wkeep": "no", "wperm_only": "no", "wflat": "no",
           "shufr_wkeep": "ORACLE", "shufr_wperm": "ORACLE", "shufr_wflat": "ORACLE",
           "a1_wkeep": "ORACLE", "a1_wperm": "ORACLE"}
    for k in ("a0_wkeep", "wperm_only", "wflat", "shufr_wkeep", "shufr_wperm",
              "shufr_wflat", "a1_wkeep", "a1_wperm"):
        v = g(k)
        m, lo, hi = _boot(v - a0, rng)
        w = int((v < a0 - 1e-9).sum()); l = int((v > a0 + 1e-9).sum())
        P(f"{k:<16s}{lab[k]:>15s}{v.mean():>9.3f}{np.median(v):>9.3f}"
          f"{f'{m:+.3f} [{lo:+.3f},{hi:+.3f}]':>26s}{f'{w}/{l}':>9s}")
    P("")
    d1 = g("shufr_wperm") - g("shufr_wkeep")
    m, lo, hi = _boot(d1, rng)
    P(f"  the WEIGHT permutation alone, holding the residual permutation fixed:")
    P(f"     shufr_wperm - shufr_wkeep = {m:+.3f} [{lo:+.3f}, {hi:+.3f}]")
    d2 = g("wperm_only") - a0
    m2, lo2, hi2 = _boot(d2, rng)
    P(f"  the WEIGHT permutation alone, with the REAL distogram and NO residual permutation:")
    P(f"     wperm_only - a0           = {m2:+.3f} [{lo2:+.3f}, {hi2:+.3f}]")
    d3 = g("wflat") - a0
    m3, lo3, hi3 = _boot(d3, rng)
    P(f"  UNIFORM weights, real distogram (NATIVE-FREE, therefore deployable if it helps):")
    P(f"     wflat - a0                = {m3:+.3f} [{lo3:+.3f}, {hi3:+.3f}]")
    P("")
    P("  READ.  If `wperm_only` and `wflat` carry most of `shuf_paired`'s advantage, the")
    P("  coordinator's 'errors are harmful because of their ASSIGNMENT to pairs' is confounded")
    P("  with 'the distogram's 1/sd^2 confidences are anti-informative' -- a different claim")
    P("  with a different, and deployable, prescription.")
    P("")


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else "run"
    if m == "run":
        run(sys.argv[2] if len(sys.argv) > 2 else None)
    else:
        report()
