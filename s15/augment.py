"""SPRINT 15, coordinator -- AUGMENT the retrieval pool with restraint-SOLVED conformers.

THE ARITHMETIC THAT MOTIVATES THIS, WHICH IS ALREADY MEASURED.

Sprint 12 established the terminal operator's behaviour as a regression with R^2 = 0.89:

    d_out  =  1.16 * d_set_mean  +  0.04 * d_set_best

The emitted structure tracks the **mean** of the aggregated set, and is almost blind to its
**best** member. A perfect rank-1 selection is worth -1.74 A through an argmin and only -0.03 A
through the m=75 average. This is why fourteen sprints of work on ranking bought so little: the
pipeline's last stage does not consume ranking.

But it consumes the set MEAN, and that gives a lever nobody has pulled. On the 126 targets:

    the retrieval circular-mean start the fits are handed      4.072 A
    a restraint-fitted structure (raw predicted distogram)     3.644 A

**A restraint-fitted conformer is 0.43 A better than the set it would be averaged into.** By the
regression above, mixing structures better than the set mean into the set must move the set mean
down, and the emitted structure with it -- and the fits are *independent* of the retrieval errors,
because they are solved from a learned distance predictor rather than drawn from a sequence-matched
library. Averaging correlated errors cancels nothing; averaging independent ones does.

That is the whole hypothesis, and it is worth stating how modest it is: this does not need the
generative arm to BEAT the incumbent. K1-CORRECTED showed it does not (3.644 vs 3.204, +0.440
[+0.290, +0.592]). It only needs the fits to beat the **set mean** of the pool they join, which is
a much weaker requirement and is already measured to hold by 0.43 A.

WHY THIS IS NOT DOUBLE-COUNTING. The distogram already appears in the incumbent, as the FILTER
that selects 75 windows from 500. Here it appears as a **generator** of new geometry that was never
in the library. The two uses are different in kind: filtering can only reorder structures the
library happens to contain, whereas solving can produce a conformation the library does not contain
at all. The `distogram_filter_only` arm is included precisely so this distinction is measured
rather than argued.

ARMS, all native-free:

    incumbent_avg      coordinate average of the production top-75          the control
    fits_only          coordinate average of the restraint-fitted ensemble
    augment_w          average of (top-75 + n_fit fits), fits weighted w    the arm under test
    augment_replace    replace the worst-objective k of the 75 with fits
    ORACLE_bestfit     the average with the single best-RMSD fit added      ORACLE, ceiling only

`w` and `k` are chosen LEAVE-FOLD-OUT, on the four training folds, by RMSD on those folds only --
which is legitimate because the folds are pinned and disjoint, and is exactly how every other
learned component in this project is fitted. The fold being reported never contributes to the
choice.

WHAT WOULD FALSIFY IT. If `augment_w` does not beat `incumbent_avg` at any w, then either the fits'
errors are not independent of the retrieval errors after all -- they share the distogram, so this is
a live possibility and is the most likely failure mode -- or the set-mean regression does not
extrapolate to sets containing structures of a different kind. Both are worth knowing, and the
second would be a genuine limit on Sprint 12's most-used law.

Run:
    python -m s15.augment
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
from s15 import robust as Rb                 # noqa: E402
from s15 import seed as SD                   # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)

WEIGHTS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0)      # weight per fit, relative to a pool member
REPLACE = (0, 10, 25, 50)


def top75(pdb, n):
    """The production top-75 window set, exactly as the incumbent aggregates it."""
    u = I.load_univ(pdb)
    p = I.pool_idx(u, 500)
    rec = I.shipped_record(pdb)
    sub = np.asarray(rec["sub"], int)
    return np.asarray(u["W"], float)[p][sub], u


def fit_ensemble(pdb, seq, n, fold, dhat, sd, i, j, n_fit=8, loss="squared", scale=1.0):
    """A diverse ensemble of restraint-solved conformers. NATIVE-FREE throughout.

    Returns the structures and their objective values; the objective is kept because the
    `augment_replace` arm needs an ordering and it must not be an RMSD ordering.
    """
    rng = SD.stable_rng(pdb, "augment")
    starts = D.starts(pdb, seq, n, fold, n_fit, rng)
    ens, objs = [], []
    for phi0, psi0, _tag in starts:
        phi, psi, f = Rb.fit_robust(dhat, sd, i, j, phi0, psi0, loss, scale)
        ens.append(np.asarray(I.build_ca(phi, psi), float))
        objs.append(f)
    return np.asarray(ens), np.asarray(objs)


def _avg(sets, weights):
    """Weighted version of the production aggregation operator.

    `I.coordinate_average(W, P)` takes a PAIRWISE-RMSD matrix as its second argument, not
    weights -- it superposes on the medoid and takes an UNWEIGHTED mean. To weight the members
    we reproduce the operator explicitly: same medoid rule, same superposition, weighted mean.
    With uniform weights this is bit-identical to `I.coordinate_average`, which `_avg_selfcheck`
    asserts.
    """
    W = np.concatenate(sets, axis=0)
    w = np.concatenate([np.full(len(s), float(wi)) for s, wi in zip(sets, weights)])
    keep = w > 0
    W = W[keep]; w = w[keep]
    b = I.medoid(I.pairwise_rmsd(W))
    S = I.superpose_batch(W, W[b])
    return (S * w[:, None, None]).sum(0) / w.sum()


def _avg_selfcheck(W):
    """Uniform weights must reproduce the production operator exactly."""
    a = _avg([np.asarray(W, float)], [1.0])
    b, _ = I.coordinate_average(np.asarray(W, float))
    return float(np.abs(a - b).max())


def run_target(t, deb, n_fit=8, weights=WEIGHTS, replace=REPLACE):
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    T75, u = top75(pdb, n)
    nat = np.asarray(u["nat_ca"], float)
    i, j = I.pair_index(n)
    sep = (j - i).astype(float)
    dg = I.distogram(pdb, seq, fold)
    dhat = np.maximum(np.asarray(dg["expected"], float) - deb(sep), 2.0)
    sd = np.maximum(np.asarray(dg["sd"], float), 1e-3)

    ens, objs = fit_ensemble(pdb, seq, n, fold, dhat, sd, i, j, n_fit=n_fit)
    #: the incumbent's own scores, used only to decide which pool members are worst -- the
    #: production filter score, never an RMSD.
    rec = I.shipped_record(pdb)
    sc = rec.get("score", rec.get("scores"))
    pool_score = (np.asarray(sc, float)[:len(T75)] if sc is not None
                  else np.zeros(len(T75)))

    #: DIVERSITY CHECK, mandatory. The set-mean law (d_out = 1.16 d_set_mean + 0.04 d_set_best)
    #: was fitted on sets with real diversity and is OUT OF DOMAIN on a collapsed set: the QGEOM
    #: workstream measured a CVaR-VQE returning 2.7 distinct structures out of 2,048 draws whose
    #: "set mean" looked 0.85 A better while its actual coordinate average was worse. This module
    #: rests on that law, so the diversity of the augmenting set is reported beside every RMSD.
    P_fit = I.pairwise_rmsd(ens)
    P_pool = I.pairwise_rmsd(T75)
    iu = np.triu_indices(len(ens), 1)
    row = {"pdb": pdb, "n": n, "fold": fold,
           "set_mean_top75": float(I.kabsch_rmsd_batch(T75, nat).mean()),
           "fit_mean": float(np.mean([I.ca_rmsd(c, nat) for c in ens])),
           "fit_best_ORACLE": float(np.min([I.ca_rmsd(c, nat) for c in ens])),
           "fit_spread": float(P_fit[iu].mean()),
           "fit_frac_pairs_distinct": float((P_fit[iu] > 0.25).mean()),
           "pool_spread": float(P_pool[np.triu_indices(len(T75), 1)].mean()),
           "arms": {}}

    def emit(Cavg):
        return float(I.ca_rmsd(I.project(Cavg, seq, fold)["fit_ca"], nat))

    row["arms"]["incumbent_avg"] = emit(_avg([T75], [1.0]))
    row["arms"]["fits_only"] = emit(_avg([ens], [1.0]))
    for w in weights:
        row["arms"][f"augment_w{w}"] = emit(_avg([T75, ens], [1.0, w]))
    order = np.argsort(-pool_score) if pool_score.any() else np.arange(len(T75))
    for k in replace:
        keep = T75 if k == 0 else T75[order[:max(1, len(T75) - k)]]
        row["arms"][f"augment_replace{k}"] = emit(_avg([keep, ens], [1.0, 1.0]))
    #: ORACLE ceiling -- the single best fit added at high weight. Diagnostic only.
    b = int(np.argmin([I.ca_rmsd(c, nat) for c in ens]))
    row["arms"]["ORACLE_bestfit"] = emit(_avg([T75, ens[b][None]], [1.0, 8.0]))
    return row


def run(targets=None, n_fit=8):
    from s14 import ladder as L
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    fail = np.isin(pdbs, I.FAIL18)

    data = C.gather(tg)
    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sep, _fn=fn: -_fn(np.zeros_like(sep), sep))

    rows = []
    path = os.path.join(RESULTS, "augment.json")
    for c, t in enumerate(tg):
        rows.append(run_target(t, deb[int(t["fold"])], n_fit=n_fit))
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": rows, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    arms = list(rows[0]["arms"])
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    base = np.asarray([r["arms"]["incumbent_avg"] for r in rows], float)

    # ---- leave-fold-out choice of w and k, on RMSD over TRAINING folds only
    def lfo(prefix, opts):
        picked, vals = {}, []
        for f in sorted(set(folds)):
            tr = folds != f
            best = min(opts, key=lambda o: np.mean(
                [rows[q]["arms"][f"{prefix}{o}"] for q in range(len(rows)) if tr[q]]))
            picked[int(f)] = best
        for q, r in enumerate(rows):
            vals.append(r["arms"][f"{prefix}{picked[int(folds[q])]}"])
        return np.asarray(vals, float), picked

    v_w, pick_w = lfo("augment_w", WEIGHTS)
    v_k, pick_k = lfo("augment_replace", REPLACE)

    out = {"n": len(rows), "incumbent": float(ref.mean()), "rows": rows,
           "lfo_weight": pick_w, "lfo_replace": pick_k, "arms": {}}
    named = {a: np.asarray([r["arms"][a] for r in rows], float) for a in arms}
    named["augment_LFO"] = v_w
    named["augment_replace_LFO"] = v_k
    for a, v in named.items():
        out["arms"][a] = {**I.summary(v), "median": float(np.median(v)),
                          "FAIL18": float(v[fail].mean()),
                          "frac_under_2_5": float((v < 2.5).mean()),
                          "vs_incumbent_avg": I.paired(v, base, folds=folds, names=pdbs),
                          "vs_shipped": I.paired(v, ref, folds=folds, names=pdbs)}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_augment", out, n_expected=len(tg))

    sm = float(np.mean([r["set_mean_top75"] for r in rows]))
    fm = float(np.mean([r["fit_mean"] for r in rows]))
    fs = float(np.mean([r["fit_spread"] for r in rows]))
    ps = float(np.mean([r["pool_spread"] for r in rows]))
    fd = float(np.mean([r["fit_frac_pairs_distinct"] for r in rows]))
    out["fit_spread_mean"] = fs
    out["pool_spread_mean"] = ps
    out["fit_frac_pairs_distinct_mean"] = fd
    print(f"\nn = {len(rows)}   incumbent {ref.mean():.3f}   "
          f"top-75 set mean {sm:.3f}   fit mean {fm:.3f}   "
          f"(the hypothesis needs fit mean < set mean: {'YES' if fm < sm else 'NO'})")
    print("DIVERSITY CHECK -- the set-mean law is out of domain on a collapsed set:")
    print(f"  fit ensemble mean pairwise RMSD {fs:.3f} A, retrieval pool {ps:.3f} A, "
          f"fraction of fit pairs separated by >0.25 A: {fd:.3f}")
    if fd < 0.5:
        print("  ** WARNING: the fit ensemble is largely collapsed, so the set-mean "
              "argument that motivates this module does NOT apply, and the arms "
              "below must not be read through it. **")
    print(f"leave-fold-out weight per fold: {pick_w}    replace-k per fold: {pick_k}\n")
    print(f"{'arm':<24}{'mean':>8}{'median':>8}{'<2A':>7}{'<2.5A':>8}{'FAIL18':>9}"
          f"{'vs incumbent avg':>26}")
    for a in ["incumbent_avg", "fits_only"] + [f"augment_w{w}" for w in WEIGHTS] + \
             [f"augment_replace{k}" for k in REPLACE] + \
             ["augment_LFO", "augment_replace_LFO", "ORACLE_bestfit"]:
        s = out["arms"][a]
        v = s["vs_incumbent_avg"]
        print(f"{a:<24}{s['mean']:>8.3f}{s['median']:>8.3f}{s['frac_under_2.0']:>7.2f}"
              f"{s['frac_under_2_5']:>8.2f}{s['FAIL18']:>9.3f}"
              f"   {v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]")
    print("\n(only `augment_LFO` and `augment_replace_LFO` are honest headline numbers; the"
          "\n individual w and k rows are the curve they were chosen from, and ORACLE_bestfit"
          "\n is a ceiling that can never be reported as a prediction.)")
    return out


if __name__ == "__main__":
    run()
