"""SPRINT 16, ENERGY -- RESOLVING LEGACY'S ELEVEN WEIGHTS HONESTLY.

THE PROBLEM.  `core.energy.DEFAULT_WEIGHTS` is documented in the repository as "a
variance-balanced starting point, not a fit"; `WEIGHT_ORIGIN` marks seven of the eleven
as `empirical ... MUST be calibrated`, and `energy_quality.calibrate_weights` -- the
function that was supposed to set them -- was never run.  Comparing a NEVER-FITTED
statistical potential against a force field whose parameters were fitted for decades
compares TRAINING STATUS, not physics.

THE PATH TAKEN: **fit them, leave-fold-out, and report the ablation.**  The alternative
(declare them fixed and limit every Legacy claim) is the fallback if the fit does not
generalise, and this module measures which of the two the evidence supports.

THE INSTRUMENT.  The 19 exhaustively enumerated targets -- 9 at n = 9 (262,144 configs)
and 10 at n = 10 (1,048,576) -- with all eleven Legacy components and an exact ORACLE
CA-RMSD stored for EVERY configuration.  Complete spaces, so an argmin is a CERTIFIED
optimum, not a search result.

LEAKAGE.  The fit consumes ORACLE CA-RMSD as a training LABEL.  That is permitted only
leave-fold-out and only with the label declared, which it is here: `fit_fold(f)` trains
on the targets whose pinned fold != f and is scored ONLY on the targets whose fold == f.
No held-out target's native touches its own weight vector.  Every number produced by the
fitted arm carries the FITTED tag; every number from `DEFAULT_WEIGHTS` carries UNFITTED.
**A fitted and an unfitted system are never presented as equally trained.**

PRE-REGISTERED (written before the first fold returned):
  * Fit = ordinary least squares of the eleven PER-TARGET-STANDARDISED components onto
    the PER-TARGET-STANDARDISED ORACLE RMSD, pooled over the training targets.  Per-target
    standardisation is what stops one target's energy scale from owning the fit; it is the
    same device `LegacyField._standardise` already uses.
  * Ridge lambda is fixed at 1e-3 x trace(X'X)/p, one value, not tuned.
  * Success = the FITTED weights beat UNFITTED on held-out folds on BOTH of
    (a) mean Spearman rho(E, ORACLE RMSD), and (b) the certified-argmin RMSD, with a
    paired CI over the 19 targets excluding zero on (b).
  * If the fit wins on rho and NOT on argmin, that is reported as a rho-only win and the
    ranking claim stays dead -- improving a proxy while the ranking does not move is this
    project's recorded "better matrix, worse ranking" failure mode.

    python -m s16.energy_weights
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s15.seed import stable_rng            # noqa: E402
from s16 import energy_lib as L            # noqa: E402

NSUB = 40000            # configurations sampled per target for the FIT (pre-registered)
RIDGE = 1e-3


def _load(path):
    z = np.load(path)
    return dict(pdb=str(z["pdb"]), n=int(z["n"]), fold=int(z["fold"]),
                rmsd=np.asarray(z["rmsd"], np.float32),
                comp=np.column_stack([np.asarray(z["leg_" + t], np.float32)
                                      for t in L.LEG_TERMS]),
                snap_index=int(z["snap_index"]))


def _std(x):
    x = np.asarray(x, np.float64)
    s = x.std(0)
    s = np.where(s > 0, s, 1.0)
    return (x - x.mean(0)) / s


def fit_weights(train, ridge=RIDGE):
    """OLS/ridge of standardised components onto standardised ORACLE RMSD.

    Returns an eleven-vector on the ORIGINAL component scale, so it drops straight into
    `legacy_total_from`.  Sign convention: lower energy should mean lower RMSD, so the
    fitted coefficient on a component that RISES with RMSD is POSITIVE.
    """
    Xs, ys, scale = [], [], []
    for d in train:
        rng = stable_rng("s16energy", "wfit", d["pdb"])
        ix = rng.permutation(len(d["rmsd"]))[:NSUB]
        C = d["comp"][ix].astype(np.float64)
        sd = C.std(0); sd = np.where(sd > 0, sd, 1.0)
        Xs.append((C - C.mean(0)) / sd)
        ys.append(_std(d["rmsd"][ix].astype(np.float64)[:, None])[:, 0])
        scale.append(sd)
    X = np.vstack(Xs); y = np.concatenate(ys)
    A = X.T @ X
    lam = ridge * np.trace(A) / A.shape[0]
    beta = np.linalg.solve(A + lam * np.eye(A.shape[0]), X.T @ y)
    # back to the raw component scale, using the mean per-target sd
    sd_mean = np.mean(np.vstack(scale), 0)
    return beta / sd_mean, beta


def evaluate(d, w, tag):
    """Score one target's FULL enumeration under weight vector `w`.  ORACLE RMSD is the
    evaluation label only."""
    E = (d["comp"].astype(np.float64) @ np.asarray(w, float))
    r = d["rmsd"].astype(np.float64)
    amin, ntie = L.argmin_tied(E, r)
    pct = float((E < E[d["snap_index"]]).mean())
    return {"arm": tag, "pdb": d["pdb"], "fold": d["fold"], "n": d["n"],
            "rho": L.spearman(E, r), "argmin_rmsd": amin, "n_tied": ntie,
            "mean_rmsd": float(r.mean()), "best_rmsd": float(r.min()),
            "argmin_minus_mean": amin - float(r.mean()),
            "ORACLE_native_pctile": pct,
            "top100_mean_rmsd": float(r[np.argsort(E)[:100]].mean())}


def run():
    data = [_load(p) for p in L.ENUM_FILES]
    print(f"loaded {len(data)} enumerated targets, "
          f"{sum(len(d['rmsd']) for d in data):,} labelled configurations", flush=True)
    w0 = L.legacy_weight_vector()                       # UNFITTED
    folds = sorted({d["fold"] for d in data})
    rows, fitted_w = [], {}
    for f in folds:
        tr = [d for d in data if d["fold"] != f]
        te = [d for d in data if d["fold"] == f]
        if not te:
            continue
        w, beta = fit_weights(tr)
        fitted_w[int(f)] = {"w": w.tolist(), "beta_std": beta.tolist(),
                            "n_train_targets": len(tr)}
        print(f"fold {f}: trained on {len(tr)} targets, scoring {len(te)}", flush=True)
        for d in te:
            rows.append(evaluate(d, w0, "UNFITTED"))
            rows.append(evaluate(d, w, "FITTED_LOFO"))
    return rows, fitted_w, w0


def report(rows, fitted_w, w0):
    import collections
    by = collections.defaultdict(dict)
    for r in rows:
        by[r["pdb"]][r["arm"]] = r
    pdbs = sorted(by)
    fold = np.array([by[p]["UNFITTED"]["fold"] for p in pdbs])
    out = {"n_targets": len(pdbs), "per_target": rows,
           "fitted_weights_by_fold": fitted_w,
           "unfitted_weights": dict(zip(L.LEG_TERMS, w0.tolist()))}
    for key in ("rho", "argmin_rmsd", "argmin_minus_mean", "ORACLE_native_pctile",
                "top100_mean_rmsd"):
        a = np.array([by[p]["FITTED_LOFO"][key] for p in pdbs], float)
        b = np.array([by[p]["UNFITTED"][key] for p in pdbs], float)
        out[key] = {"FITTED_LOFO_mean": float(a.mean()), "UNFITTED_mean": float(b.mean()),
                    "paired_fitted_minus_unfitted": I.paired(a, b, folds=fold, names=pdbs)}
    # the "random member" reference: mean RMSD of the enumeration
    out["mean_rmsd_of_space"] = float(np.mean(
        [by[p]["UNFITTED"]["mean_rmsd"] for p in pdbs]))
    out["best_rmsd_of_space"] = float(np.mean(
        [by[p]["UNFITTED"]["best_rmsd"] for p in pdbs]))
    return out


if __name__ == "__main__":
    rows, fw, w0 = run()
    rep = report(rows, fw, w0)
    p = L.write("energy_weights", rep)
    print(f"\nwrote {p}")
    for key in ("rho", "argmin_rmsd", "argmin_minus_mean", "ORACLE_native_pctile"):
        d = rep[key]; pr = d["paired_fitted_minus_unfitted"]
        print(f"{key:24s} FITTED {d['FITTED_LOFO_mean']:+.4f}  UNFITTED "
              f"{d['UNFITTED_mean']:+.4f}  diff {pr['mean_diff']:+.4f} "
              f"[{pr['ci95'][0]:+.4f},{pr['ci95'][1]:+.4f}] W/L {pr['n_better']}/{pr['n_worse']}")
    print(f"mean RMSD of the enumerated space (the random control): "
          f"{rep['mean_rmsd_of_space']:.4f}; best {rep['best_rmsd_of_space']:.4f}")
