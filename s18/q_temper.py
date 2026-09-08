"""SPRINT 18 / QUANTUM-ADVERSARIAL -- ATTACK 5b: TEMPERING THE DISTOGRAM'S OWN WEIGHTS.

WHY THIS ARM EXISTS, AND WHAT IT DISTINGUISHES.

Attack 4 (`q_ceilattack.py`) found that at IDENTICAL residual RMS, an error model that fixes
the pairs the distogram was already CONFIDENT about reaches 2.145 A while one that fixes its
LOW-confidence pairs reaches 2.697 A -- so the confidence ORDERING carries real information
about where accuracy pays.

Attack 5 (`q_ceilattack2.py`) asks the opposite-looking question: is the `1/sd^2` WEIGHTING of
those same confidences, inside the fit, helping or hurting?

Both can be true at once, and they are different claims with different designs:

  READING 1 -- the confidences are informative but `1/sd^2` OVER-weights them.  Then a tempered
               weight `sd^-p` with 0 < p < 2 should beat BOTH uniform (p = 0) and the deployed
               (p = 2), and `p` is a deployable parameter.
  READING 2 -- the confidences are miscalibrated AS MAGNITUDES even though their ordering is
               informative.  Then the response should be monotone with no interior optimum:
               flat wins, and the ordering has to be consumed some other way.

The sweep distinguishes them directly.  The fit's loss is `sum_p ((d_p - dhat_p)/s_p)^2`, so
passing `s = sd^(p/2)` realises the weight `sd^-p` exactly:

    p = 0    uniform            (identical to `wflat`; a global weight scale cannot move an argmin)
    p = 0.5  strongly tempered
    p = 1.0  tempered
    p = 2.0  the deployed objective
    p = 3, 4 SHARPER than deployed -- included because attack 5's `wflat` came out WORSE than
             the deployed weighting, so the response may be monotone INCREASING in p and the
             deployed value may not be the optimum in the direction anyone expected.

**NATIVE-FREE.**  Every arm here reads the real `dhat` and the real `sd`.  Nothing reads a
native distance.  If an interior `p` wins, it is DEPLOYABLE, not a ceiling -- which is why it is
run at n = 126 and reported with a paired fold-clustered interval rather than as a diagnostic.

MANDATORY CONTROLS.  Zero-information: `p = 0`, which deletes the confidence information
entirely.  Matched-random: `wperm` -- the deployed `1/sd^2` weights PERMUTED across pairs, so
the weight multiset and its entropy are identical and only the assignment is destroyed.

RUN:  python -m s18.q_temper run
      python -m s18.q_temper report
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

SALT = "s18qtemper"
OUT = os.path.join(RESULTS, "q_temper.json")
PS = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0)


def _rmsd_of(phi, psi, nat):
    return float(I.ca_rmsd(I.build_ca(phi, psi), nat))


def cell(t, data, deb):
    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    d = data[pdb]
    i, j, sd, nat = d["i"], d["j"], np.asarray(d["sd"], float), np.asarray(d["nat"], float)
    dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
    rng = SD.stable_rng(pdb, "temper", salt=SALT)

    W = np.asarray(AV.top75_windows(pdb)[0], float)
    P = I.pairwise_rmsd(W)
    avg, _b = I.coordinate_average(W, P)
    pr = I.project(np.asarray(avg, float), seq, fold)
    phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)

    e = {"pdb": pdb, "n": n, "fold": fold,
         "avg": float(I.ca_rmsd(np.asarray(avg, float), nat)),
         "sd_min": float(sd.min()), "sd_max": float(sd.max()),
         "sd_cv": float(sd.std() / max(sd.mean(), 1e-12))}
    for p in PS:
        s_eff = sd ** (p / 2.0)                     # weight = sd^-p
        p_, q_, _f = AL.fit(dhat, s_eff, i, j, phi0, psi0)
        e[f"p{p}"] = _rmsd_of(p_, q_, nat)
    # MATCHED-RANDOM: the deployed weights permuted (identical multiset and entropy)
    pi = rng.permutation(len(sd))
    p_, q_, _f = AL.fit(dhat, sd[pi], i, j, phi0, psi0)
    e["wperm"] = _rmsd_of(p_, q_, nat)
    return e


def extend(ps):
    """Add missing `p` rungs to an existing complete sweep, in place, without re-running the
    rungs already measured.  Used because the first launch of this module ran PS = (0..2) and
    the sweep past the deployed exponent -- the only direction in which the monotone curve
    could still have an optimum -- had to be added afterwards."""
    d = json.load(open(OUT))
    rows = d["rows"]
    tg = [t for t in I.targets() if t["pdb"] in {r["pdb"] for r in rows}]
    data = C.gather(tg)
    pdbs = [t["pdb"] for t in tg]
    deb = {}
    for f in sorted({int(t["fold"]) for t in tg}):
        train = [q for q in pdbs if data[q]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))
    by = {r["pdb"]: r for r in rows}
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb, seq, fold = t["pdb"], t["seq"], int(t["fold"])
        r = by[pdb]
        if all(f"p{x}" in r for x in ps):
            continue
        dd = data[pdb]
        i, j, sd = dd["i"], dd["j"], np.asarray(dd["sd"], float)
        nat = np.asarray(dd["nat"], float)
        dhat = np.maximum(dd["dhat"] - deb[fold](dd["sep"]), 2.0)
        W = np.asarray(AV.top75_windows(pdb)[0], float)
        P = I.pairwise_rmsd(W)
        avg, _b = I.coordinate_average(W, P)
        pr = I.project(np.asarray(avg, float), seq, fold)
        phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)
        for x in ps:
            if f"p{x}" in r:
                continue
            p_, q_, _f = AL.fit(dhat, sd ** (x / 2.0), i, j, phi0, psi0)
            r[f"p{x}"] = _rmsd_of(p_, q_, nat)
        if (c + 1) % 10 == 0:
            json.dump({"rows": rows, "complete": False}, open(OUT, "w"), default=float)
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"rows": rows, "complete": True}, open(OUT, "w"), default=float)
    report(rows)


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


def _bootf(dif, folds, rng, B=4000):
    """Fold-CLUSTERED paired bootstrap."""
    dif = np.asarray(dif, float); folds = np.asarray(folds)
    uf = np.unique(folds); idx = {f: np.flatnonzero(folds == f) for f in uf}
    out = np.empty(B)
    for b in range(B):
        take = []
        for f in rng.choice(uf, size=uf.size, replace=True):
            m = idx[f]
            take.append(rng.choice(m, size=m.size, replace=True))
        out[b] = dif[np.concatenate(take)].mean()
    return float(dif.mean()), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rng = SD.stable_rng("temper", "report", salt=SALT)
    g = lambda k: np.array([r[k] for r in rows])              # noqa: E731
    folds = [r["fold"] for r in rows]
    P = print
    dep = g("p2.0")
    P("=" * 100)
    P(f"ATTACK 5b -- TEMPERING THE DISTOGRAM'S WEIGHTS.  weight = sd^-p.   n = {len(rows)}")
    P("            EVERY ARM IS NATIVE-FREE: real dhat, real sd, no native distance read.")
    P("=" * 100)
    P(f"  coordinate average (start)        {g('avg').mean():.3f}")
    P(f"  sd dispersion: CV {np.mean([r['sd_cv'] for r in rows]):.3f}, "
      f"min {np.mean([r['sd_min'] for r in rows]):.3f}, "
      f"max {np.mean([r['sd_max'] for r in rows]):.3f}")
    P("")
    P(f"{'arm':<10s}{'weight':<12s}{'RMSD':>9s}{'median':>9s}"
      f"{'vs deployed p=2':>28s}{'W/L':>9s}{'folds':>7s}")
    arms = [(f"p{p}", f"sd^-{p}") for p in PS] + [("wperm", "1/sd^2 permuted")]
    for k, lab in arms:
        v = g(k)
        m, lo, hi = _bootf(v - dep, folds, rng)
        w = int((v < dep - 1e-9).sum()); l = int((v > dep + 1e-9).sum())
        uf = np.unique(folds)
        fm = np.array([(v - dep)[np.asarray(folds) == f].mean() for f in uf])
        same = int((np.sign(fm) == np.sign(m)).sum()) if abs(m) > 1e-12 else uf.size
        P(f"{k:<10s}{lab:<12s}{v.mean():>9.3f}{np.median(v):>9.3f}"
          f"{f'{m:+.3f} [{lo:+.3f},{hi:+.3f}]':>28s}{f'{w}/{l}':>9s}"
          f"{f'{same}/{uf.size}':>7s}")
    P("")
    flat = g("p0.0")
    best = min(PS, key=lambda p: g(f"p{p}").mean())
    P(f"  best p on the mean: p = {best}  ({g(f'p{best}').mean():.3f} A)")
    m, lo, hi = _bootf(flat - g("wperm"), folds, rng)
    P(f"  uniform (p=0) - permuted 1/sd^2 : {m:+.3f} [{lo:+.3f}, {hi:+.3f}]")
    P("")
    P("  READING 1 (confidences informative, 1/sd^2 OVER-weights them) is supported if an")
    P("  INTERIOR p beats BOTH p = 0 and p = 2.  READING 2 (miscalibrated as magnitudes) is")
    P("  supported if the response is monotone with the optimum at p = 0.")
    P("")
    P("  A CORRECTION TO THIS MODULE'S OWN PRE-REGISTERED READING, PRESERVED IN PLACE.  The")
    P("  pre-registration said 'if uniform merely MATCHES the permuted weights, sd carries no")
    P("  usable weighting information either way'.  THAT INFERENCE IS INVALID and the")
    P("  coordinator caught it.  Uniform and permuted are TWO WAYS OF DESTROYING THE SAME")
    P("  INFORMATION, so their agreeing is a consistency check on the null, not evidence of")
    P("  absence.  What decides the question is that BOTH are significantly WORSE than the")
    P("  deployed weighting.  The correct reading: sd carries real usable weighting")
    P("  information, destroying it costs ~0.15 A, and the information is in WHICH PAIR gets")
    P("  which weight rather than in the spread of the weight values.")
    P("")


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else "run"
    if m == "run":
        run(sys.argv[2] if len(sys.argv) > 2 else None)
    elif m == "extend":
        extend([float(x) for x in (sys.argv[2:] or ["3.0", "4.0"])])
    else:
        report()
