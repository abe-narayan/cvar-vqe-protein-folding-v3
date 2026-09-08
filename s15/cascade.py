"""SPRINT 15, coordinator -- the full cascade: GENERATION -> SELECTION -> AGGREGATION -> FINAL.

WHY THIS MODULE EXISTS. The brief requires every architecture to report four distinct numbers
and the three gaps between them, because a single final RMSD hides which stage is actually
failing:

    G  what quality exists somewhere in the generated ensemble   (ORACLE readout)
    S  what the method can identify without the native            (objective argmin)
    A  what an estimator recovers by combining the ensemble       (coordinate consensus)
    F  what is finally reported, after physics                    (validity + refinement)

Sprint 14 measured that **aggregation is worth 3.4x the objective** (coordinate averaging
buys 0.587 A where the whole contribution of any objective is 0.171 A) and that a *perfect*
ranker emits 2.474 A through a decile but **1.219 A through a top-100** -- a 1.26 A swing from
the terminal operator alone. So `A` is where the leverage is, and `S - A` is the number that
decides whether selection matters at all.

THE ARCHITECTURE BEING MEASURED. Restraint-driven **generation**, not candidate selection:

  1. build a per-pair distance likelihood from one or more native-free channels
     (the learned distogram, the retrieval pool's empirical histogram, or their product),
     optionally separation-debiased leave-fold-out;
  2. fit continuous torsions to it by L-BFGS from many native-free starts, using the exact
     analytic torsion gradient -- **each start gives an ensemble member**;
  3. read G, S and A off that ensemble;
  4. apply physics: Legacy as a clash gate (its only measured-honest role) and a single
     AMBER relaxation, reported separately so neither can silently inflate the headline.

This is the classical limit of Families A and C. **It must be measured before any quantum arm
is built**, because Sprint 14's decisive negative was that running a variational optimiser is
worse than best-of-N from an untrained circuit -- so best-of-N from *this* ensemble is the
control every quantum claim will have to beat.

STRICT DISCIPLINE. Selection uses the OBJECTIVE ONLY, never the RMSD. `G` reads the native and
is labelled ORACLE everywhere it appears. Physical validity is checked, and any structure
excluded is counted and reported, per the brief's requirement that exclusions never silently
improve a mean.

Run:
    python -m s15.cascade
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
from s15 import seed as SD             # noqa: E402
from s15 import distgeo as D                 # noqa: E402
from s15 import distml as M                  # noqa: E402
from s15 import pooldist as P                # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)


# ------------------------------------------------------------------ native-free start bank
def start_bank(pdb, seq, n, fold, n_start, rng):
    """A diverse bank of native-free initialisations.

    Diversity matters here in a way it did not for a single fit: the ensemble IS the product,
    so starts that converge to the same basin are wasted members.
    """
    from s14 import retprior as R
    PHI, PSI, _ = R.windows(pdb, "top75")
    out = [(R.circ_mean(PHI, axis=0), R.circ_mean(PSI, axis=0), "retrieval_circmean"),
           (np.full(n, np.deg2rad(-63.0)), np.full(n, np.deg2rad(-42.0)), "helix"),
           (np.full(n, np.deg2rad(-120.0)), np.full(n, np.deg2rad(130.0)), "extended")]
    m = PHI.shape[0]
    idx = rng.permutation(m)
    for s in range(max(0, n_start - len(out))):
        k = int(idx[s % m])
        out.append((PHI[k].copy(), PSI[k].copy(), f"pool_{k}"))
    return out[:n_start]


def build_likelihood(pdb, seq, n, fold, channel, debias=None):
    """Per-pair distance log-likelihood from native-free channels."""
    dg = I.distogram(pdb, seq, fold)
    centres = np.asarray(dg["centres"], float)
    i, j = I.pair_index(n)
    sep = (j - i).astype(float)
    tab_dg = M.LogPTable(np.asarray(dg["prob"], float), centres)
    if debias is not None:
        tab_dg = M.ShiftedLogP(tab_dg, debias(sep))
    if channel == "distogram":
        return tab_dg, i, j
    Dm, sim, _i, _j = P.pool_distances(pdb, n)
    tab_pool = P.hist_table(Dm, centres)
    if channel == "pool":
        return tab_pool, i, j
    if channel == "combined":
        return P.SumLogP(tab_dg, tab_pool), i, j
    raise ValueError(channel)


# --------------------------------------------------------------------------- validity gate
def clash_fraction(CA, cutoff=3.6):
    """Fraction of |i-j|>=2 CA pairs closer than `cutoff`. A minimal steric gate."""
    n = CA.shape[0]
    i, j = I.pair_index(n)
    d = np.sqrt(((CA[i] - CA[j]) ** 2).sum(1))
    return float((d < cutoff).mean())


def run_target(t, channel="combined", n_start=16, debias=None, gate=0.02):
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)
    tab, i, j = build_likelihood(pdb, seq, n, fold, channel, debias)
    rng = SD.stable_rng(pdb, channel)
    S = start_bank(pdb, seq, n, fold, n_start, rng)

    ens, objs, clashes = [], [], []
    for phi0, psi0, _tag in S:
        phi, psi, f = M.fit_ml(tab, i, j, phi0, psi0)
        ca = np.asarray(I.build_ca(phi, psi), float)
        ens.append(ca); objs.append(f); clashes.append(clash_fraction(ca))
    ens = np.asarray(ens); objs = np.asarray(objs); clashes = np.asarray(clashes)
    rms = np.asarray([I.ca_rmsd(c, nat) for c in ens])          # ORACLE, post hoc only

    keep = clashes <= gate
    n_excluded = int((~keep).sum())
    if not keep.any():                                          # never let a gate empty a target
        keep = np.ones(len(ens), bool)
        n_excluded = 0

    # ---- the four numbers
    G = float(rms.min())                                        # ORACLE
    S_idx = int(np.argmin(np.where(keep, objs, np.inf)))
    S_val = float(rms[S_idx])
    Ck, _ = I.coordinate_average(ens[keep])
    A_val = float(I.ca_rmsd(Ck, nat))
    Cp = I.project(Ck, seq, fold)["fit_ca"]
    F_val = float(I.ca_rmsd(Cp, nat))

    # aggregation over the objective-best half, the Sprint 12 "set mean not set best" law
    order = np.argsort(np.where(keep, objs, np.inf))
    half = max(2, int(keep.sum() // 2))
    Ch, _ = I.coordinate_average(ens[order[:half]])
    A_half = float(I.ca_rmsd(Ch, nat))

    return {"pdb": pdb, "n": n, "fold": fold, "G_ORACLE": G, "S": S_val, "A": A_val,
            "A_besthalf": A_half, "F": F_val, "n_excluded": n_excluded,
            "clash_mean": float(clashes.mean()), "ens_mean_rmsd": float(rms.mean()),
            "ens_spread": float(rms.std())}


def run(targets=None, channel="combined", n_start=16, use_debias=True):
    from s14 import ladder as L
    from s15 import distcal as C
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    fail = np.isin(pdbs, I.FAIL18)

    debias_by_fold = None
    if use_debias:
        data = C.gather(tg)
        debias_by_fold = {}
        for f in sorted(set(folds)):
            train = [p for p in pdbs if data[p]["fold"] != f]
            fn, _ = C.fit_correction(data, train, "sep")
            #: fn(x, sep) = x - bias(sep); we need +bias to shift the SUPPORT
            debias_by_fold[f] = (lambda sep, _fn=fn: -_fn(np.zeros_like(sep), sep))

    rows = []
    path = os.path.join(RESULTS, f"cascade_{channel}.json")
    for c, t in enumerate(tg):
        db = debias_by_fold[int(t["fold"])] if debias_by_fold else None
        rows.append(run_target(t, channel=channel, n_start=n_start, debias=db))
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": rows, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    keys = ["G_ORACLE", "S", "A", "A_besthalf", "F"]
    out = {"n": len(rows), "channel": channel, "n_start": n_start,
           "use_debias": use_debias, "incumbent": float(ref.mean()),
           "rows": rows, "stages": {}}
    for k in keys:
        v = np.asarray([r[k] for r in rows], float)
        out["stages"][k] = {**I.summary(v), "median": float(np.median(v)),
                            "FAIL18": float(v[fail].mean()),
                            "frac_under_2_5": float((v < 2.5).mean()),
                            "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs)}
    g = {k: np.asarray([r[k] for r in rows], float) for k in keys}
    out["gaps"] = {"G_to_S": float((g["S"] - g["G_ORACLE"]).mean()),
                   "S_to_A": float((g["A"] - g["S"]).mean()),
                   "A_to_F": float((g["F"] - g["A"]).mean())}
    out["excluded_total"] = int(sum(r["n_excluded"] for r in rows))
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write(f"s15_cascade_{channel}", out, n_expected=len(tg))

    print(f"\nchannel={channel}  starts={n_start}  debias={use_debias}   "
          f"incumbent {ref.mean():.3f}   structures excluded by the clash gate: "
          f"{out['excluded_total']}\n")
    print(f"{'stage':<14}{'mean':>8}{'median':>8}{'<2A':>7}{'<2.5A':>8}{'FAIL18':>9}"
          f"{'vs incumbent':>26}")
    for k in keys:
        s = out["stages"][k]
        v = s["vs_incumbent"]
        print(f"{k:<14}{s['mean']:>8.3f}{s['median']:>8.3f}{s['frac_under_2.0']:>7.2f}"
              f"{s['frac_under_2_5']:>8.2f}{s['FAIL18']:>9.3f}"
              f"   {v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]")
    print(f"\ngaps:  G->S {out['gaps']['G_to_S']:+.3f}   S->A {out['gaps']['S_to_A']:+.3f}"
          f"   A->F {out['gaps']['A_to_F']:+.3f}      (the largest is the next target)")
    return out


if __name__ == "__main__":
    run()
