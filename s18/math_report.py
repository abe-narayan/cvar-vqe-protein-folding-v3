"""SPRINT 18 / MATH -- aggregation and the two mandatory controls for the lattice arm.

Reads the sealed `s18/results/math_lattice.json` (and `math_mu.json` when present) and emits
`s18/results/math_report.json` plus the tables quoted in `s18/math_FINDINGS.md`.

CONTROLS, both required by the brief on every arm:

  ZERO-INFORMATION REFERENCE -- the mean RMSD of a uniformly drawn configuration of the same
  enumerated space.  Any objective's argmin must be read against it.

  MATCHED-RANDOM OPERATION -- a random Walsh subspace projection retaining the SAME variance
  fraction as the strict weight-<=1 truncation, with the constant kept.  This prices "throw
  away 39% of the objective's variance at random" and is the correct null for "throwing away
  the high-degree part helped".  If the matched-random projection also improves the argmin,
  the degree-1 story is a story about DAMAGE, not about degree.

RUN:  python -m s18.math_report            (aggregate + controls)
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from s18 import math_lib as M
from s15 import seed as SD                       # noqa: E402


def _load(tag):
    p = os.path.join(M.RESULTS, f"math_{tag}.json")
    if not os.path.exists(p):
        return None
    d = json.load(open(p))
    if not d.get("_complete"):
        print(f"  WARNING: {tag} artefact is NOT complete -- refusing to aggregate")
        return None
    return d


def matched_random_control(n_draw=24):
    """Random Walsh subspace projections at the strict truncation's retained variance."""
    from s16 import qphase_lib as QP
    from s14 import vqe_lib as V
    out = {}
    for pdb in QP.TARGETS19:
        ins = QP.inst(pdb)
        nq, N = ins.n_qubits, ins.N
        E = V.uniformise(ins.hamil().astype(np.float64))
        rmsd = ins.rmsd
        c = M.fwht(E)
        deg = M.popcount(np.arange(N, dtype=np.int64))
        p = c ** 2
        tot = float(p.sum() - p[0])
        keep_le1 = deg <= 1
        frac = float(p[keep_le1 & (deg > 0)].sum() / tot)
        vals = []
        for t in range(n_draw):
            rng = SD.stable_rng(pdb, t, "randproj", salt=M.SALT)
            # greedily add randomly ordered non-constant coefficients until the retained
            # variance fraction first reaches `frac` -- matched in magnitude, random in place
            order = rng.permutation(np.flatnonzero(deg > 0))
            cs = np.cumsum(p[order])
            m = int(np.searchsorted(cs, frac * tot)) + 1
            cc = np.zeros(N)
            cc[0] = c[0]
            cc[order[:m]] = c[order[:m]]
            g = M.ifwht(cc)
            vals.append(float(rmsd[g == g.min()].mean()))
        out[pdb] = {"retained_var_frac": frac, "rand_proj_argmin_ORACLE": vals,
                    "rand_proj_mean": float(np.mean(vals)),
                    "zero_info_reference_mean_rmsd_ORACLE": float(rmsd.mean()),
                    "space_best_ORACLE": float(rmsd.min())}
        print(f"  {pdb} frac {frac:.3f}  randproj {np.mean(vals):.3f} "
              f"(sd {np.std(vals):.3f})  zero-info {rmsd.mean():.3f}", flush=True)
        M.ck("control", f"ctrl_{pdb}", out[pdb])
    if len(out) == len(QP.TARGETS19):
        M.seal("control", {"module": "s18/math_report.py", "n_draw": n_draw,
                           "salt": M.SALT})
    return out


def continuous_varfrac(n_eval=4096):
    """The continuous analogue of the lattice retained-variance table.

    A FRESH sample from the same reference measure (a different stable_rng stream from the
    frozen quadrature, so this is an out-of-sample check of the projection, not a refit),
    then `Var_mu(E_trunc) / Var_mu(E_full)` and the orthogonality residual
    `|Var(E) - Var(E_le1) - Var(E_ge2)| / Var(E)`.  Native-free.
    """
    from s18 import math_anova as A
    from s18 import math_iface as MI
    from s12 import instrument as I
    from s15 import seed as SD
    tg = I.targets()
    done = M.ck_load("varfrac")
    for t in tg:
        pdb = t["pdb"]
        if f"vf_{pdb}" in done:
            continue
        p = os.path.join(M.CACHE, f"anova_{pdb}_{A.MU_DEFAULT}_{A.NSAMP}_{A.GRID}.npz")
        if not os.path.exists(p):
            continue
        o = MI.build(pdb, target=t)
        tt = o.t
        rng = SD.stable_rng(pdb, "varfrac", n_eval, salt=M.SALT)
        PHI, PSI = A.mu_samples(pdb, tt.seq, tt.fold, tt.n, tt.mu, n_eval, rng)
        ef = np.array([o.E_full(PHI[k], PSI[k])[0] for k in range(n_eval)])
        er = np.array([o.E_res(PHI[k], PSI[k])[0] for k in range(n_eval)])
        ea = np.array([o.E_ang(PHI[k], PSI[k])[0] for k in range(n_eval)])
        v = ef.var()
        row = {"pdb": pdb, "n": tt.n, "fold": tt.fold,
               "var_frac_residue": float(er.var() / max(v, 1e-300)),
               "var_frac_angle": float(ea.var() / max(v, 1e-300)),
               "orth_resid_residue": float(abs(v - er.var() - (ef - er).var()) / max(v, 1e-300)),
               "orth_resid_angle": float(abs(v - ea.var() - (ef - ea).var()) / max(v, 1e-300)),
               "corr_res_full": float(np.corrcoef(er, ef)[0, 1]),
               "corr_ang_full": float(np.corrcoef(ea, ef)[0, 1]),
               "spearman_res_full": M.spearman(er, ef),
               "spearman_ang_full": M.spearman(ea, ef)}
        M.ck("varfrac", f"vf_{pdb}", row)
        print(f"  {pdb} n{tt.n}  varfrac res {row['var_frac_residue']:.3f}  "
              f"ang {row['var_frac_angle']:.3f}  orth {row['orth_resid_residue']:.2e}  "
              f"rho {row['spearman_res_full']:.3f}", flush=True)
    d = M.ck_load("varfrac")
    ks = [k for k in d if k.startswith("vf_")]
    if ks:
        print(f"\n  n = {len(ks)}   residue-additive retains "
              f"{np.mean([d[k]['var_frac_residue'] for k in ks]):.3f} of Var_mu(E); "
              f"angle-additive {np.mean([d[k]['var_frac_angle'] for k in ks]):.3f}")
        print(f"  Spearman with the full objective on mu draws: residue "
              f"{np.mean([d[k]['spearman_res_full'] for k in ks]):.3f}, angle "
              f"{np.mean([d[k]['spearman_ang_full'] for k in ks]):.3f}")
    if len(ks) == len(tg):
        M.seal("varfrac", {"module": "s18/math_report.py", "n_eval": n_eval})
    return d


def aggregate():
    d = _load("lattice")
    if d is None:
        raise SystemExit("lattice artefact missing or incomplete")
    ks = sorted(k for k in d if k.startswith("cell_"))
    pdbs = [d[k]["pdb"] for k in ks]
    folds = np.asarray([d[k]["fold"] for k in ks], int)
    g = lambda tag, f: np.array([d[k][tag][f] for k in ks])     # noqa: E731
    rep = {"n_targets": len(ks), "pdbs": pdbs}

    print("=" * 92)
    print("L1/L2  THE BRIDGE, VERIFIED EXACTLY (max |ANOVA - Walsh| relative to max|E|)")
    print("=" * 92)
    for tag in ("uniformised", "raw"):
        a = g(tag, "L1_rel").max()
        b = g(tag, "L2_rel").max()
        o1 = g(tag, "L3_max_orth_residual_qubit").max()
        o2 = g(tag, "L3_max_orth_residual_residue").max()
        v1 = g(tag, "L3_var_budget_resid_qubit").max()
        v2 = g(tag, "L3_var_budget_resid_residue").max()
        print(f"  {tag:<13} qubit-ANOVA == Walsh w<=1     worst {a:.2e}")
        print(f"  {'':<13} residue-ANOVA == w<=1+intra   worst {b:.2e}")
        print(f"  {'':<13} orthogonality residual        worst {max(o1, o2):.2e}")
        print(f"  {'':<13} Var(E)=Var(P E)+Var(E-P E)    worst {max(v1, v2):.2e}")
        rep[f"bridge_{tag}"] = {"L1_worst_rel": float(a), "L2_worst_rel": float(b),
                                "orth_worst": float(max(o1, o2)),
                                "var_budget_worst": float(max(v1, v2))}

    print()
    print("=" * 92)
    print("L3/L4  THE TWO OBJECTS ARE DIFFERENT, AND ONLY ONE OF THEM IS THE 2.411")
    print("=" * 92)
    for tag in ("uniformised", "raw"):
        full = g(tag, "argmin_rmsd_full_ORACLE")
        w1 = g(tag, "argmin_rmsd_walsh_le1_ORACLE")
        res = g(tag, "argmin_rmsd_residue_add_ORACLE")
        vq = g(tag, "L3_var_frac_qubit_le1")
        vr = g(tag, "L3_var_frac_residue_le1")
        print(f"\n  --- {tag} ---   retained variance: strict w<=1 {vq.mean():.3f}   "
              f"residue-additive {vr.mean():.3f}   "
              f"(the strict object DISCARDS {100*(vr.mean()-vq.mean())/vr.mean():.0f}% "
              f"of the per-residue field)")
        print(f"  {'object':<26}{'argmin RMSD':>12}{'median':>9}"
              f"{'vs full [95% CI, paired]':>32}{'W/L/T':>10}{'fold sign':>11}{'rho w RMSD':>12}")
        rows = [("full objective", full, g(tag, "rho_full_rmsd_ORACLE")),
                ("strict Walsh weight<=1", w1, g(tag, "rho_walsh_le1_rmsd_ORACLE")),
                ("residue-additive ANOVA", res, g(tag, "rho_residue_add_rmsd_ORACLE"))]
        for name, v, rho in rows:
            if name == "full objective":
                print(f"  {name:<26}{v.mean():>12.3f}{np.median(v):>9.3f}"
                      f"{'--':>32}{'--':>10}{'--':>11}{rho.mean():>12.3f}")
                continue
            s = M.paired(v, full, folds)
            print(f"  {name:<26}{v.mean():>12.3f}{np.median(v):>9.3f}"
                  f"   {s['mean_diff']:+.3f} [{s['ci95'][0]:+.3f},{s['ci95'][1]:+.3f}]"
                  f"{s['n_better']:>5}/{s['n_worse']}/{s['n_tie']}"
                  f"{s['folds_same_sign']:>11}{rho.mean():>12.3f}")
            rep[f"{tag}_{name}"] = s
        rep[f"{tag}_var_frac"] = {"strict_w1": float(vq.mean()),
                                  "residue_additive": float(vr.mean())}
        rep[f"{tag}_means"] = {"full": float(full.mean()), "walsh_le1": float(w1.mean()),
                               "residue_add": float(res.mean())}

    print()
    print("=" * 92)
    print("L5  ENCODING INVARIANCE -- is the strict object a property of the PHYSICS?")
    print("=" * 92)
    W = np.array([d[k]["L5"]["walsh_le1_argmin_rmsd_ORACLE"] for k in ks])   # (T, P)
    R = np.array([d[k]["L5"]["residue_add_argmin_rmsd_ORACLE"] for k in ks])
    inv = all(d[k]["L5"]["residue_add_invariant"] for k in ks)
    base_w = np.array([d[k]["L5"]["unpermuted_walsh_le1"] for k in ks])
    full_u = g("uniformised", "argmin_rmsd_full_ORACLE")
    m = W.mean(0)
    print(f"  residue-additive projection is EXACTLY invariant on all targets : {inv}")
    print(f"  strict weight<=1, mean over 19 targets, per relabelling         : "
          f"{m.mean():.3f} +- {m.std():.3f}   [{m.min():.3f}, {m.max():.3f}]")
    print(f"  the SHIPPED labelling                                           : "
          f"{base_w.mean():.3f}   <-- the published 2.411")
    print(f"  the FULL objective                                              : "
          f"{full_u.mean():.3f}")
    print(f"  relabellings whose 19-target mean is WORSE than the full object : "
          f"{np.mean(m > full_u.mean()):.3f}")
    print(f"  percentile of the shipped labelling inside the relabelling null : "
          f"{np.mean(m <= base_w.mean()):.3f}")
    print(f"  per-target sd of the strict argmin across relabellings, mean    : "
          f"{W.std(1).mean():.3f} A   (residue-additive: {R.std(1).mean():.3e})")
    rep["L5"] = {"residue_add_invariant": bool(inv),
                 "walsh_null_mean": float(m.mean()), "walsh_null_sd": float(m.std()),
                 "walsh_null_min": float(m.min()), "walsh_null_max": float(m.max()),
                 "shipped_labelling": float(base_w.mean()),
                 "full": float(full_u.mean()),
                 "frac_relabellings_worse_than_full": float(np.mean(m > full_u.mean())),
                 "shipped_percentile_in_null": float(np.mean(m <= base_w.mean())),
                 "per_target_sd_mean": float(W.std(1).mean())}

    c = _load("control")
    if c is not None:
        print()
        print("=" * 92)
        print("CONTROLS")
        print("=" * 92)
        rp = np.array([c[f"ctrl_{p}"]["rand_proj_mean"] for p in pdbs])
        zi = np.array([c[f"ctrl_{p}"]["zero_info_reference_mean_rmsd_ORACLE"]
                       for p in pdbs])
        sb = np.array([c[f"ctrl_{p}"]["space_best_ORACLE"] for p in pdbs])
        w1 = g("uniformised", "argmin_rmsd_walsh_le1_ORACLE")
        s = M.paired(w1, rp, folds)
        print(f"  zero-information reference (uniform draw)   {zi.mean():.3f}")
        print(f"  space best (ORACLE ceiling)                 {sb.mean():.3f}")
        print(f"  matched-random projection, same var frac    {rp.mean():.3f}")
        print(f"  strict weight<=1 vs matched-random          {s['mean_diff']:+.3f} "
              f"[{s['ci95'][0]:+.3f},{s['ci95'][1]:+.3f}]  "
              f"{s['n_better']}/{s['n_worse']}/{s['n_tie']}")
        rep["controls"] = {"zero_info": float(zi.mean()), "space_best": float(sb.mean()),
                           "rand_proj": float(rp.mean()), "w1_vs_randproj": s}

    mu = _load("mu")
    if mu is not None:
        print()
        print("=" * 92)
        print("L6  REFERENCE-MEASURE SENSITIVITY (residue-additive ANOVA under three mu)")
        print("=" * 92)
        full_u = g("uniformised", "argmin_rmsd_full_ORACLE")
        rep["mu"] = {}
        print(f"  {'mu':<10}{'argmin RMSD':>12}{'median':>9}"
              f"{'vs full [95% CI, paired]':>32}{'W/L/T':>10}{'fold sign':>11}{'var frac':>10}{'qubit-fact':>12}")
        for name in ("uniform", "pool", "rama"):
            try:
                v = np.array([mu[f"mu_{p}"][name]["argmin_rmsd_residue_add_ORACLE"]
                              for p in pdbs])
                vf = np.array([mu[f"mu_{p}"][name]["var_frac_residue_le1_mu_weighted"]
                               for p in pdbs])
                qf = all(mu[f"mu_{p}"][name]["qubit_factorises"] for p in pdbs)
            except KeyError:
                print(f"  {name:<10} unavailable")
                continue
            s = M.paired(v, full_u, folds)
            print(f"  {name:<10}{v.mean():>12.3f}{np.median(v):>9.3f}"
                  f"   {s['mean_diff']:+.3f} [{s['ci95'][0]:+.3f},{s['ci95'][1]:+.3f}]"
                  f"{s['n_better']:>5}/{s['n_worse']}/{s['n_tie']}"
                  f"{s['folds_same_sign']:>11}{vf.mean():>10.3f}{str(qf):>12}")
            rep["mu"][name] = {"mean": float(v.mean()), "paired_vs_full": s,
                               "var_frac": float(vf.mean()), "qubit_factorises": bool(qf)}

    sc = _load("anova_selfcheck")
    if sc is not None:
        rep["anova_selfcheck"] = sc["checks"]

    M.ck("report", "aggregate", rep)
    M.seal("report", {"module": "s18/math_report.py",
                      "source_hashes": {"lattice": d.get("_config_hash")}})
    return rep


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("all", "control"):
        matched_random_control()
    if mode in ("all", "varfrac"):
        continuous_varfrac()
    if mode in ("all", "agg"):
        aggregate()
