"""SPRINT 14, coordinator -- can anything NATIVE-FREE supply the per-target sign?

THE ONE DIRECTION WITH LEVERAGE, AND THIS IS THE CHEAPEST TEST OF IT.

C19 established the sprint's deepest result: near-native ordering is learnable to 0.986 WITHIN
a target with an overfitting gap of -0.0005, and collapses to 0.600 ACROSS targets. Capacity is
saturated by a linear model. The one direction left with any leverage is therefore not a better
objective but a **conditioning signal**: anything that supplies, at inference, the per-target
SIGN of the in-band ordering axis.

The flip diagnostic named that axis precisely. Per-target skill correlates **+0.909 with the
native's z-scored radius of gyration** and +0.951 with rho(contacts, RMSD). So the sign is a
compactness sign -- whether, for this particular target, more-compact structures are better or
worse than the population expectation.

**But the radius of gyration used there is a NATIVE quantity.** The whole question is whether
any native-free predictor of it exists. And one does, in principle, already in the pipeline:
the shipped leave-fold-out distogram predicts every CA-CA pair distance, and a radius of
gyration follows from pair distances in closed form:

    Rg^2  =  (1 / (2 N^2))  *  sum_ij  d_ij^2

So the distogram's `expected` column IS a native-free Rg predictor, for free, and nobody has
ever read it that way. A second, independent native-free proxy is the mean Rg of the retrieved
window pool.

WHAT THIS DECIDES.
  * If a native-free proxy tracks the native Rg well, **the per-target sign is predictable**,
    C19's bottleneck has a key, and the next sprint has a concrete target.
  * If it does not, then the one remaining direction is closed too, and the sprint's negative
    is complete in the strongest possible sense -- not "we did not find a conditioning signal"
    but "the obvious one does not exist".

Either answer is decisive, which is why it is worth the twenty minutes.

NATIVE-FREE except for the ORACLE_ scoring column, which is the thing being predicted.

Run:
    python -m s14.signpred
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

RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(RESULTS, exist_ok=True)


def rg_from_coords(X):
    """Radius of gyration of an (n, 3) trace."""
    X = np.asarray(X, float)
    return float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))


def rg_from_pair_distances(n, i, j, d, fill_short=True):
    """Rg from pair distances, exact:  Rg^2 = (1/(2 N^2)) * sum_ij d_ij^2.

    The distogram covers only |i-j| >= 2, so the short-separation terms are filled with
    ideal-geometry constants (3.80 A adjacent, 6.0 A for i,i+2 -- the latter varies with
    conformation but contributes a small, nearly constant share and is applied identically
    to every arm, so it cannot create a spurious correlation ACROSS targets).
    """
    d = np.asarray(d, float)
    tot = 2.0 * float((d ** 2).sum())          # each (i,j) counted once -> double it
    if fill_short:
        tot += 2.0 * (n - 1) * (3.80 ** 2)     # |i-j| = 1
        if n > 2:
            tot += 2.0 * (n - 2) * (6.00 ** 2)  # |i-j| = 2 approximation
    return float(np.sqrt(tot / (2.0 * n * n)))


def spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    den = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / den) if den > 0 else 0.0


def pearson(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a - a.mean(); b = b - b.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else 0.0


def boot_ci(a, b, fn=pearson, n_boot=4000, seed=0):
    a = np.asarray(a, float); b = np.asarray(b, float)
    rng = np.random.default_rng(seed)
    n = len(a)
    vals = [fn(a[k], b[k]) for k in (rng.integers(0, n, n) for _ in range(n_boot))]
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows = []
    for t in tg:
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        u = I.load_univ(pdb)

        # ---- ORACLE: the quantity to be predicted
        rg_native = rg_from_coords(u["nat_ca"])

        # ---- native-free proxy 1: the shipped leave-fold-out distogram's expected distances
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        rg_disto = rg_from_pair_distances(n, i, j, dg["expected"])

        # ---- native-free proxy 2: mean Rg of the retrieved window pool
        p = I.pool_idx(u)
        W = np.asarray(u["W"], float)[p]
        rg_pool = float(np.mean([rg_from_coords(w) for w in W]))

        # ---- native-free proxy 3: the incumbent's own emitted structure
        rec = I.shipped_record(pdb)
        rg_emit = rg_from_coords(np.asarray(rec["fit_ca"], float))

        # ---- native-free control: chain length alone (compactness scales with n)
        rows.append({"pdb": pdb, "n": n, "fold": fold,
                     "rg_native_ORACLE": rg_native, "rg_disto": rg_disto,
                     "rg_pool": rg_pool, "rg_emit": rg_emit})

    def col(k):
        return np.asarray([r[k] for r in rows], float)

    nat = col("rg_native_ORACLE")
    proxies = {"distogram-predicted Rg": col("rg_disto"),
               "retrieval pool mean Rg": col("rg_pool"),
               "incumbent emitted Rg": col("rg_emit"),
               "chain length alone (control)": col("n")}

    # The sign is a per-target property RELATIVE to what length alone predicts, so the
    # decisive quantity is the LENGTH-RESIDUALISED correlation: does the proxy know anything
    # about this target's compactness beyond how long it is?
    def residualise(x, n_):
        A = np.column_stack([np.ones_like(n_), n_])
        beta, *_ = np.linalg.lstsq(A, x, rcond=None)
        return x - A @ beta

    nvec = col("n")
    nat_r = residualise(nat, nvec)

    out = {"n_targets": len(rows), "per_target": rows, "raw": {}, "length_residualised": {}}
    for name, v in proxies.items():
        out["raw"][name] = {"pearson": pearson(v, nat), "spearman": spearman(v, nat),
                            "ci95": boot_ci(v, nat)}
        if name.startswith("chain length"):
            continue
        vr = residualise(v, nvec)
        out["length_residualised"][name] = {
            "pearson": pearson(vr, nat_r), "spearman": spearman(vr, nat_r),
            "ci95": boot_ci(vr, nat_r)}

    with open(os.path.join(RESULTS, "signpred.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_signpred", out, n_expected=len(rows))

    print(f"Can a NATIVE-FREE proxy predict the native radius of gyration?  "
          f"n = {len(rows)} targets\n")
    print(f"{'proxy':<32}{'pearson':>10}{'spearman':>10}{'CI95':>24}")
    for name, s in out["raw"].items():
        print(f"{name:<32}{s['pearson']:>10.3f}{s['spearman']:>10.3f}"
              f"   [{s['ci95'][0]:+.3f},{s['ci95'][1]:+.3f}]")
    print(f"\n-- LENGTH-RESIDUALISED (the decisive test: does the proxy know this target's "
          f"compactness\n   beyond how long it is?) --")
    print(f"{'proxy':<32}{'pearson':>10}{'spearman':>10}{'CI95':>24}")
    for name, s in out["length_residualised"].items():
        print(f"{name:<32}{s['pearson']:>10.3f}{s['spearman']:>10.3f}"
              f"   [{s['ci95'][0]:+.3f},{s['ci95'][1]:+.3f}]")
    print("\nReference: per-target in-band skill correlates +0.909 with the native Rg "
          "z-score (ORACLE).")
    return out


if __name__ == "__main__":
    run()
