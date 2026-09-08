"""SPRINT 15, coordinator -- CAN THE COHERENT ERROR COMPONENT BE REMOVED NATIVE-FREE?

WHAT MAKES THIS THE RIGHT NEXT EXPERIMENT.

`s15/errstruct.py` destroyed one property of the distogram's real error vector at a time and refitted.
On a 6-target read, one arm dominated everything else:

    surrogate               mean    vs real
    real                   2.832    +0.000
    ORACLE_shuffle_signs   1.825    **-1.007 [-1.314, -0.664]**
    ORACLE_gauss_matched   2.501    -0.332 [-0.868, +0.195]
    ORACLE_shuffle_pairs   2.577    -0.256 [-0.884, +0.431]
    ORACLE_winsor_z3       2.711    -0.121 [-0.322, +0.057]
    debias_sep             3.012    +0.179 [-0.361, +0.744]   (native-free, and WORSE)

**Keeping every error magnitude exactly and randomising only the SIGNS is worth about an ångström**,
with a confidence interval excluding zero at n = 6. Outliers are worth little, separation structure
is worth little, and matching an i.i.d. Gaussian of the same RMS is worth little. The expensive
property is that the errors **agree in direction**: the distogram is not noisy about the answer so
much as consistently wrong in one direction per target, and a least-squares fit follows it there.

The measured mean bias on those targets is **+0.865 Å**, so the structure is over-predicted in size,
which is the same pathology K2 found by separation and the opposite of the 25.8% contraction that
coordinate averaging produces.

**And the obvious fix does not work.** `debias_sep` — the leave-fold-out separation-bias profile,
which is the correction this project already had — makes it slightly WORSE (+0.179). A profile fitted
across targets removes the *average* coherent component, and the coherent component is **per target**:
subtracting a population mean from a target whose own bias is larger or smaller than the mean moves
it the wrong way about as often as the right way. That is why K2's correction bought so little and it
is a result worth stating on its own.

THE QUESTION THIS MODULE ASKS. `ORACLE_shuffle_signs` is an ORACLE diagnostic — it needs the true
error to destroy it, so it is not a method. It says the coherent per-target component is worth ~1 Å.
**Can that component be estimated without the native?**

THE CHANNEL THAT COULD DO IT. The K = 500 retrieval pool supplies a complete CA-CA distance matrix
per member and therefore an independent estimate of the target's overall SIZE, from 500
sequence-matched fragments, with no training. K3 measured it as slightly more accurate than the
distogram (MAE 2.249 vs 2.386) with the **opposite bias sign** at short and medium separation. K6
then measured it as a poor *objective* (the native sits at its 66.5th percentile). Those are two
different uses. **This is a third**: not as an objective and not as a restraint, but as a
native-free reference for a single scalar per target.

ARMS. `dhat` is corrected as `d' = a * dhat + b`, with `(a, b)` estimated per target from a
native-free reference, then fitted exactly as `s15/distgeo.py` fits the raw arm.

    real              no correction                                          NATIVE-FREE, the control
    debias_sep        the existing leave-fold-out separation profile          NATIVE-FREE
    pool_scale        a from the pool's distances, b = 0                      NATIVE-FREE
    pool_shift        a = 1, b from the pool's distances                      NATIVE-FREE
    pool_affine       both, least squares against the pool                    NATIVE-FREE
    pool_scale_rg     a from the ratio of radii of gyration                   NATIVE-FREE
    ORACLE_scale      a chosen against the TRUE distances                     ORACLE ceiling
    ORACLE_affine     both chosen against the TRUE distances                  ORACLE ceiling
    ORACLE_perfect    dhat replaced by dtrue                                  ORACLE, 0.611 A

The two ORACLE correction arms are the informative ceilings. `ORACLE_scale` says how much of the
1 Å is reachable by ANY per-target isotropic rescaling; if it is small, the coherent component is not
a scale error and the pool arms cannot possibly work, which would be worth learning in one run rather
than three. `ORACLE_perfect` is the floor of the whole channel.

WHAT WOULD FALSIFY THE IDEA. If `ORACLE_scale` recovers most of the ångström but the pool arms
recover none of it, the component is real and estimable in principle but the pool is the wrong
reference. If `ORACLE_scale` itself recovers little, the coherent component is not a scale at all and
this whole direction closes — a clean negative, obtained cheaply.

Run:
    python -m s15.scale
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

from s12 import instrument as I              # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import distgeo as D                 # noqa: E402
from s15 import pooldist as P                # noqa: E402
from s15 import seed as SD                   # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)

ARMS = ("real", "debias_sep", "pool_scale", "pool_shift", "pool_affine", "pool_scale_rg",
        "ORACLE_scale", "ORACLE_affine", "ORACLE_perfect")


def _scale_to(dhat, ref, w):
    """Least-squares multiplicative factor taking `dhat` onto `ref`."""
    return float((w * dhat * ref).sum() / max((w * dhat * dhat).sum(), 1e-12))


def _affine_to(dhat, ref, w):
    """Least-squares (a, b) with ref ~ a*dhat + b, weighted."""
    X = np.column_stack([dhat, np.ones_like(dhat)])
    W = np.diag(w) if False else w[:, None]
    A = X.T @ (W * X)
    y = X.T @ (w * ref)
    try:
        a, b = np.linalg.solve(A, y)
    except np.linalg.LinAlgError:
        return 1.0, 0.0
    return float(a), float(b)


def rg_from_dist(d, n):
    """Radius of gyration from a complete CA-CA distance set: Rg^2 = (1/n^2) sum_{i<j} d_ij^2.

    Only the |i-j|>=2 pairs are available here, so this is a consistent proxy rather than the
    exact Rg -- it is used only as a RATIO between two distance sets, where the missing near
    pairs largely cancel.
    """
    return float(np.sqrt((d * d).sum() / max(n * n, 1)))


def run(targets=None, n_start=4):
    from s14 import ladder as L
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    fail = np.isin(pdbs, I.FAIL18)

    data = C.gather(tg)
    bias = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        bias[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    res = {a: [] for a in ARMS}
    coefs = []
    path = os.path.join(RESULTS, "scale.json")

    for c, t in enumerate(tg):
        p = t["pdb"]; d = data[p]
        i, j, sep, sd = d["i"], d["j"], d["sep"], d["sd"]
        dtrue, dhat, n = d["dtrue"], d["dhat"], d["n"]
        w = 1.0 / sd ** 2

        Dm, _sim, _i, _j = P.pool_distances(p, n)
        d_pool = np.median(Dm, axis=0)                       # NATIVE-FREE reference

        a_pool = _scale_to(dhat, d_pool, w)
        b_pool = float((w * (d_pool - dhat)).sum() / max(w.sum(), 1e-12))
        aa_pool, bb_pool = _affine_to(dhat, d_pool, w)
        a_rg = rg_from_dist(d_pool, n) / max(rg_from_dist(dhat, n), 1e-9)
        a_true = _scale_to(dhat, dtrue, w)
        aa_true, bb_true = _affine_to(dhat, dtrue, w)
        coefs.append({"pdb": p, "a_pool": a_pool, "a_rg": a_rg, "a_true": a_true,
                      "b_pool": b_pool, "bias": float((dhat - dtrue).mean())})

        tgt = {
            "real": dhat,
            "debias_sep": dhat - bias[d["fold"]](sep),
            "pool_scale": a_pool * dhat,
            "pool_shift": dhat + b_pool,
            "pool_affine": aa_pool * dhat + bb_pool,
            "pool_scale_rg": a_rg * dhat,
            "ORACLE_scale": a_true * dhat,
            "ORACLE_affine": aa_true * dhat + bb_true,
            "ORACLE_perfect": dtrue,
        }
        starts = D.starts(p, d["seq"], n, d["fold"], n_start, SD.stable_rng(p, "scale"))
        for a in ARMS:
            best = None
            for phi0, psi0, _tag in starts:
                phi, psi, f, _ = D.fit_distances(np.maximum(tgt[a], 2.0), w, i, j, phi0, psi0)
                if best is None or f < best[0]:
                    best = (f, float(I.ca_rmsd(I.build_ca(phi, psi), d["nat"])))
            res[a].append(best[1])
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": res, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    base = np.asarray(res["real"], float)
    out = {"n": len(pdbs), "incumbent": float(ref.mean()), "coefficients": coefs, "arms": {},
           "per_target": {a: dict(zip(pdbs, map(float, res[a]))) for a in ARMS}}
    for a in ARMS:
        v = np.asarray(res[a], float)
        out["arms"][a] = {**I.summary(v), "median": float(np.median(v)),
                          "FAIL18": float(v[fail].mean()),
                          "frac_under_2_5": float((v < 2.5).mean()),
                          "vs_real": I.paired(v, base, folds=folds, names=pdbs),
                          "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs)}
    ap = np.asarray([c["a_pool"] for c in coefs]); at = np.asarray([c["a_true"] for c in coefs])
    ar = np.asarray([c["a_rg"] for c in coefs])
    out["scale_agreement"] = {
        "a_true_mean": float(at.mean()), "a_pool_mean": float(ap.mean()),
        "a_rg_mean": float(ar.mean()),
        "corr_pool_true": float(np.corrcoef(ap, at)[0, 1]),
        "corr_rg_true": float(np.corrcoef(ar, at)[0, 1]),
        "mae_pool_true": float(np.abs(ap - at).mean()),
        "mae_rg_true": float(np.abs(ar - at).mean())}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_scale", out, n_expected=len(tg))

    sa = out["scale_agreement"]
    print(f"\nn = {len(pdbs)}   incumbent {ref.mean():.3f}")
    print(f"the per-target scale the TRUE distances want: mean {sa['a_true_mean']:.4f}")
    print(f"  the pool estimates it at {sa['a_pool_mean']:.4f} "
          f"(corr {sa['corr_pool_true']:+.3f}, MAE {sa['mae_pool_true']:.4f})")
    print(f"  the Rg ratio estimates it at {sa['a_rg_mean']:.4f} "
          f"(corr {sa['corr_rg_true']:+.3f}, MAE {sa['mae_rg_true']:.4f})\n")
    print(f"{'arm':<18}{'mean':>8}{'median':>8}{'<2A':>7}{'<2.5A':>8}{'FAIL18':>9}"
          f"{'vs real':>24}{'vs incumbent':>24}")
    for a in ARMS:
        s = out["arms"][a]
        vr, vi = s["vs_real"], s["vs_incumbent"]
        print(f"{a:<18}{s['mean']:>8.3f}{s['median']:>8.3f}{s['frac_under_2.0']:>7.2f}"
              f"{s['frac_under_2_5']:>8.2f}{s['FAIL18']:>9.3f}"
              f"  {vr['mean_diff']:+.3f} [{vr['ci95'][0]:+.3f},{vr['ci95'][1]:+.3f}]"
              f"  {vi['mean_diff']:+.3f} [{vi['ci95'][0]:+.3f},{vi['ci95'][1]:+.3f}]")
    print("\n(`ORACLE_scale` is the ceiling of ANY per-target isotropic rescaling. If it is small,"
          "\n the coherent component is not a scale error and this direction closes. If it is large"
          "\n and the pool arms recover none of it, the component is estimable in principle and the"
          "\n pool is the wrong reference.)")
    return out


if __name__ == "__main__":
    run()
