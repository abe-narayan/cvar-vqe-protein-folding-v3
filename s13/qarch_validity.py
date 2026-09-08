"""SPRINT 13 QUANTUM-ARCHITECTURE, EXPERIMENT 3 -- OBJECTIVE VALIDITY BEFORE OPTIMISATION.

Sprint 12's hardest lesson: the certified exact optimum of an objective that does not rank
the native emits a WORSE structure than doing nothing.  So before any VQE runs, every
candidate Hamiltonian is measured on a real population of configurations drawn from the
discrete torsion space:

    * the percentile of the native-snapped configuration under it,
    * the percentile of the space's BEST configuration (the thing an optimiser should find),
    * Spearman rho(energy, CA-RMSD), globally and IN-BAND,
    * whether the argmin beats a random draw, and by how much,
    * the mean RMSD of the m lowest-energy configurations, m = 1, 10, 100,
    * the binned energy-vs-RMSD scatter.

On the n=9 targets the whole 4^9 space is enumerated, so every number here is EXACT for
Legacy and the prior -- no sampling error at all.  AMBER is measured on the labelled
2,955-configuration subsample cached by `qarch_enum` (uniform / prior-sampled / oracle
near-native band), and every AMBER row says which population it is on.

NORMALISATION.  Legacy and AMBER differ in units, scale, dynamic range and length
dependence, so a combination needs a principled normalisation, not a tuned weight.  Each
channel is standardised on a NATIVE-FREE reference population (configurations drawn from
the leakage-safe empirical prior) whose mean and sd are fitted ON TRAINING FOLDS as
per-residue constants, mu_X = n*a_X and sd_X = sqrt(n)*b_X, and the combination is the
unit-weight sum of the resulting z-scores.  No weight anywhere is chosen by looking at a
reported number; the lambda sweep that is reported is an explicitly labelled diagnostic.

    python -m s13.qarch_validity
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import qarch_lib as Q          # noqa: E402

BAND_FRAC = 0.01          # ORACLE in-band definition: best 1% by true CA-RMSD
TOPM = (1, 10, 100)


def _metrics(E, rmsd, snap_idx, tag, band_idx=None, rng=None):
    """Every validity number for one (energy, rmsd) pair on one population."""
    E = np.asarray(E, float).copy(); rmsd = np.asarray(rmsd, float)
    ok = np.isfinite(E)
    if not ok.all():                 # keep indices aligned; a blown-up energy is the worst
        E[~ok] = (E[ok].max() if ok.any() else 0.0) + 1.0
    B = len(E)
    best_i = int(np.argmin(rmsd))
    out = {"population": tag, "n": B,
           "rmsd_mean": float(rmsd.mean()), "rmsd_best": float(rmsd.min()),
           "spearman_E_rmsd": Q.spearman(E, rmsd),
           "pct_of_best_config": Q.percentile_of(E[best_i], E)}
    if snap_idx is not None:
        out["pct_of_native_snap"] = Q.percentile_of(E[snap_idx], E)
    order = np.argsort(E, kind="mergesort")
    for m in TOPM:
        if m <= B:
            # average over the tied-argmin set is the correct estimator (s12 trap:
            # np.argmin on a tied signal reads the enumeration order, not the energy)
            if m == 1:
                tie = np.flatnonzero(E == E.min())
                out["top1_rmsd"] = float(rmsd[tie].mean())
                out["top1_n_tied"] = int(len(tie))
            else:
                out[f"top{m}_rmsd"] = float(rmsd[order[:m]].mean())
    out["random_draw_rmsd"] = float(rmsd.mean())
    out["argmin_minus_random"] = out["top1_rmsd"] - out["random_draw_rmsd"]
    out["argmin_rmsd_percentile"] = Q.percentile_of(out["top1_rmsd"], rmsd)
    if band_idx is not None and len(band_idx) > 8:
        out["inband_n"] = int(len(band_idx))
        out["inband_spearman"] = Q.spearman(E[band_idx], rmsd[band_idx])
        out["inband_top1_rmsd"] = float(rmsd[band_idx][np.argsort(E[band_idx])[:1]].mean())
        out["inband_mean_rmsd"] = float(rmsd[band_idx].mean())
    # binned scatter (10 energy deciles -> mean rmsd)
    q = np.quantile(E, np.linspace(0, 1, 11))
    b = np.clip(np.searchsorted(q[1:-1], E), 0, 9)
    out["scatter_decile_rmsd"] = [float(rmsd[b == d].mean()) if (b == d).any() else None
                                  for d in range(10)]
    return out


def _norm_stats(rows, key):
    """Per-residue normalisation constants fitted across targets: mu = n*a, sd = sqrt(n)*b."""
    n = np.array([r["n"] for r in rows], float)
    mu = np.array([r[key + "_mu"] for r in rows], float)
    sd = np.array([r[key + "_sd"] for r in rows], float)
    return float((mu / n).mean()), float((sd / np.sqrt(n)).mean())


def per_target(pdb_id, rng):
    path = os.path.join(Q.RESULTS, f"qarch_enum_{pdb_id}.npz")
    z = np.load(path)
    n = int(z["n"]); k = int(z["k"])
    rmsd = np.asarray(z["rmsd"], float)
    legacy = np.asarray(z["legacy"], float)
    prior = np.asarray(z["prior"], float)
    snap_idx = int(z["snap_index"])
    B = len(rmsd)
    band = np.argsort(rmsd)[:max(1, int(BAND_FRAC * B))]

    sp = Q.Space(pdb_id, k, seq=str(z["seq"]), n=n, fold=int(z["fold"]))
    Pemp = np.asarray(z["prior_table"], float)
    # a native-free reference population for the normalisation constants
    Sref = sp.sample_prior(Pemp, 4000, rng)
    ref_idx = np.ravel_multi_index(tuple(Sref.T), (k,) * n)

    rec = {"pdb": pdb_id, "n": n, "k": k, "fold": int(z["fold"]), "configs": B,
           "legacy_mu": float(legacy[ref_idx].mean()), "legacy_sd": float(legacy[ref_idx].std()),
           "prior_mu": float(prior[ref_idx].mean()), "prior_sd": float(prior[ref_idx].std()),
           "snap_rmsd": float(rmsd[snap_idx]), "best_rmsd": float(rmsd.min()),
           "mean_rmsd": float(rmsd.mean()), "arrays": path}

    # ---- single-channel objectives, EXACT over the whole space -----------
    objs = {"prior_empirical": prior, "legacy_total": legacy}
    for t in Q.LEGACY_TERMS:
        objs["leg_" + t] = np.asarray(z["leg_" + t], float)
    # ORACLE prior quality sweep -- labelled diagnostic
    S_all = np.array(np.unravel_index(np.arange(B), (k,) * n)).T.astype(np.int8)
    for q in (0.4, 0.6, 0.8, 1.0):
        objs[f"ORACLE_prior_q{q}"] = Q.prior_energy(Q.ORACLE_prior(sp, q), S_all)

    rec["full_space"] = {name: _metrics(E, rmsd, snap_idx, "uniform_full_space", band)
                         for name, E in objs.items()}
    rec["_S_all_shape"] = list(S_all.shape)

    # ---- prior-sampled population (the population a VQE actually visits) --
    Sp = sp.sample_prior(Pemp, 20000, rng)
    pidx = np.ravel_multi_index(tuple(Sp.T), (k,) * n)
    pband = np.argsort(rmsd[pidx])[:max(8, len(pidx) // 100)]
    rec["prior_population"] = {
        name: _metrics(E[pidx], rmsd[pidx], None, "prior_sampled", pband)
        for name, E in objs.items() if not name.startswith("leg_")}

    # ---- AMBER on its labelled subsample ---------------------------------
    a_idx = np.asarray(z["amber_idx"], int)
    a_kind = np.asarray(z["amber_kind"], int)
    a_tot = np.asarray(z["amber_total"], float)
    if np.isfinite(a_tot).any():
        rec["amber_mu"] = float(np.nanmean(a_tot[a_kind == 1]))
        rec["amber_sd"] = float(np.nanstd(a_tot[a_kind == 1]))
        sub = {}
        for kind, tag in ((0, "amber_uniform_subsample"), (1, "amber_prior_subsample"),
                          (2, "ORACLE_amber_nearnative_band")):
            m = a_kind == kind
            if m.sum() < 20:
                continue
            r = rmsd[a_idx[m]]
            si = int(np.flatnonzero(a_idx[m] == int(z["snap_index"]))[0]) \
                if (a_idx[m] == int(z["snap_index"])).any() else None
            sub[tag] = {
                "amber_total": _metrics(a_tot[m], r, si, tag),
                "legacy_total": _metrics(legacy[a_idx[m]], r, si, tag),
                "prior_empirical": _metrics(prior[a_idx[m]], r, si, tag)}
            for t in ("bond", "angle", "torsion", "nonbonded", "solvation"):
                key = "amb_" + t
                if key in z.files:
                    sub[tag]["amb_" + t] = _metrics(np.asarray(z[key], float)[m], r, si, tag)
        rec["amber"] = sub
    return rec


def combinations(rows, rng):
    """Variance-standardised combinations, weights derived not tuned."""
    a_leg, b_leg = _norm_stats(rows, "legacy")
    a_pri, b_pri = _norm_stats(rows, "prior")
    have_amb = all("amber_mu" in r for r in rows)
    a_amb, b_amb = _norm_stats(rows, "amber") if have_amb else (None, None)
    norm = {"legacy_per_res_mu": a_leg, "legacy_per_sqrt_res_sd": b_leg,
            "prior_per_res_mu": a_pri, "prior_per_sqrt_res_sd": b_pri,
            "amber_per_res_mu": a_amb, "amber_per_sqrt_res_sd": b_amb}
    out = []
    for r in rows:
        z = np.load(r["arrays"])
        n, k = int(z["n"]), int(z["k"])
        rmsd = np.asarray(z["rmsd"], float)
        B = len(rmsd)
        band = np.argsort(rmsd)[:max(1, int(BAND_FRAC * B))]
        zl = (np.asarray(z["legacy"], float) - n * a_leg) / (np.sqrt(n) * b_leg)
        zp = (np.asarray(z["prior"], float) - n * a_pri) / (np.sqrt(n) * b_pri)
        snap_idx = int(z["snap_index"])
        rec = {"pdb": r["pdb"], "n": n, "full_space": {}}
        rec["full_space"]["z_prior_plus_legacy"] = _metrics(zp + zl, rmsd, snap_idx,
                                                            "uniform_full_space", band)
        # DIAGNOSTIC lambda sweep -- reported, never used to pick a weight
        rec["lambda_sweep_DIAGNOSTIC"] = {
            str(lam): _metrics(zp + lam * zl, rmsd, snap_idx, "uniform_full_space", band)
                      ["spearman_E_rmsd"]
            for lam in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 1e6)}
        if have_amb:
            a_idx = np.asarray(z["amber_idx"], int); a_kind = np.asarray(z["amber_kind"], int)
            a_tot = np.asarray(z["amber_total"], float)
            m = a_kind == 0
            za = (a_tot[m] - n * a_amb) / (np.sqrt(n) * b_amb)
            zpa = zp[a_idx[m]]
            rec["amber_uniform_subsample"] = {
                "z_prior_plus_amber": _metrics(zpa + za, rmsd[a_idx[m]], None,
                                               "amber_uniform_subsample"),
                "z_prior_plus_legacy": _metrics(zpa + zl[a_idx[m]], rmsd[a_idx[m]], None,
                                                "amber_uniform_subsample")}
        out.append(rec)
    return norm, out


def main():
    rng = np.random.default_rng(0)
    pdbs = sorted(f[len("qarch_enum_"):-4] for f in os.listdir(Q.RESULTS)
                  if f.startswith("qarch_enum_") and f.endswith(".npz"))
    rows = []
    for p in pdbs:
        rows.append(per_target(p, rng))
        print(f"  {p}: legacy rho={rows[-1]['full_space']['legacy_total']['spearman_E_rmsd']:+.3f} "
              f"top1={rows[-1]['full_space']['legacy_total']['top1_rmsd']:.2f} "
              f"(mean {rows[-1]['mean_rmsd']:.2f}, best {rows[-1]['best_rmsd']:.2f}); "
              f"prior rho={rows[-1]['full_space']['prior_empirical']['spearman_E_rmsd']:+.3f}",
              flush=True)
    norm, comb = combinations(rows, rng)
    Q.write("qarch_validity", {"what": "objective validity of every candidate Hamiltonian "
                                       "over the exhaustively enumerated k=4 torsion space",
                               "band_frac": BAND_FRAC, "normalisation": norm,
                               "per_target": rows, "combinations": comb})
    return rows


if __name__ == "__main__":
    main()
