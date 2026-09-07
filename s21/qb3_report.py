"""SPRINT 21 / WORKSTREAM B -- analysis.

    python -m s21.qb3_report a      BLOCK A  -- ansatz vs bond dimension vs trainability vs RMSD
    python -m s21.qb3_report q      BLOCK Q  -- optimiser battery incl. Fisher natural gradient
    python -m s21.qb3_report e      BLOCK E  -- the encoding gauge sweep
    python -m s21.qb3_report en     BLOCK E  -- the Nelder-Mead gauge sweep

RULES APPLIED THROUGHOUT (BRIEF sections 7, 9):
  * TARGET is the unit.  Seeds are averaged WITHIN a target before any test, so the 4 seeds buy
    precision and never inflate n.
  * Mean, MEDIAN, W/L, paired target-level bootstrap CI and fold behaviour, in that order.
  * Every arm is priced against `best_of_N` at matched budget -- never an initialisation mean.
  * A |difference| below the 0.084 A MDE with a zero-spanning CI is NOT MEASURED, never "matched".
  * A comparison-only arm (Nelder-Mead, best_of_N) is EXACTLY invariant under a strictly monotone
    transform of the objective, so an AMB-vs-AMBc row for such an arm is an identity and is
    excluded from any conditioning table rather than counted as a null.
"""
from __future__ import annotations

import sys

import numpy as np

from s21 import qb3_lib as K

MDE = 0.084


def cell(rows, pdbs, key_fn, seeds=(0, 1, 2, 3), field="rmsd_ORACLE"):
    """Per-target value, averaged over seeds.  NaN where a target is missing the cell."""
    out = []
    for p in pdbs:
        v = []
        for s in seeds:
            r = rows.get(p, {}).get("rows", {}).get(key_fn(s))
            if r is not None and np.isfinite(r.get(field, np.nan)):
                v.append(float(r[field]))
        out.append(float(np.mean(v)) if v else float("nan"))
    return np.asarray(out, float)


def folds(rows, pdbs):
    return np.asarray([rows[p]["fold"] if p in rows else -1 for p in pdbs], int)


def line(name, a, b, fl, width=30):
    d = K.paired_ci(a, b, folds=fl)
    if "ci95" not in d:
        return f"{name:<{width}} n<3"
    sig = "*" if d.get("sig") else " "
    mde = "" if abs(d["mean"]) >= MDE else "  [<MDE]"
    fs = d.get("folds_same_sign", "-")
    return (f"{name:<{width}} {d['mean']:+7.3f}{sig} med {d['median']:+7.3f} "
            f"[{d['ci95'][0]:+.3f},{d['ci95'][1]:+.3f}] {d['W']:3d}W/{d['L']:<3d}L "
            f"folds {fs}/5  a={d['mean_a']:.3f} b={d['mean_b']:.3f}{mde}")


def _pdbs(rows):
    return sorted(rows)


# ============================================================ BLOCK A
def report_a():
    rows = K.read21("qb3_a")
    pdbs = _pdbs(rows)
    fl = folds(rows, pdbs)
    cfg = K.read21("qb3_config")
    print(f"BLOCK A -- ansatz zoo.  n_targets = {len(pdbs)}, seeds = {cfg['SEEDS']}, "
          f"budget = {cfg['B']}, shots = {cfg['A_SHOTS']}, alpha = {cfg['A_ALPHA']}")
    print("Readout: argmin over everything evaluated, BUILT CHAIN basis.\n")
    arms = [a for a in cfg["A_ANSATZ"]] + ["untrained_L2", "best_of_N"]
    for kind in cfg["kinds_A"]:
        print(f"=== {kind} ===")
        print("  -- geometry and trainability (measured, not quoted) --")
        print(f"  {'ansatz':<16} {'P':>4} {'chi_an':>6} {'chi_amp':>7} {'chi_p':>6} "
              f"{'H0':>6} {'H1':>6} {'|g|':>10} {'|g|sd':>10} {'cov':>4}")
        for a in arms:
            g = {}
            for f in ("n_params", "a_chi_analytic", "g1_chi_amp_max", "g1_chi_prob_max",
                      "g0_entropy_bits", "g1_entropy_bits", "gnorm_mean", "gnorm_sd"):
                v = [rows[p]["rows"][f"{kind}|{a}|{s}"].get(f)
                     for p in pdbs for s in cfg["SEEDS"]
                     if f"{kind}|{a}|{s}" in rows[p]["rows"]]
                v = [x for x in v if isinstance(x, (int, float)) and np.isfinite(x)]
                g[f] = float(np.mean(v)) if v else float("nan")
            cov = len({p for p in pdbs if f"{kind}|{a}|0" in rows[p]["rows"]})
            print(f"  {a:<16} {g['n_params']:>4.0f} {g['a_chi_analytic']:>6.1f} "
                  f"{g['g1_chi_amp_max']:>7.1f} {g['g1_chi_prob_max']:>6.1f} "
                  f"{g['g0_entropy_bits']:>6.2f} {g['g1_entropy_bits']:>6.2f} "
                  f"{g['gnorm_mean']:>10.3g} {g['gnorm_sd']:>10.3g} {cov:>4d}")
        print("\n  -- Ca-RMSD, paired vs the DEPLOYED mps_L2 (negative = the arm is better) --")
        base = cell(rows, pdbs, lambda s, a="mps_L2", k=kind: f"{k}|{a}|{s}", cfg["SEEDS"])
        bon = cell(rows, pdbs, lambda s, k=kind: f"{k}|best_of_N|{s}", cfg["SEEDS"])
        for a in arms:
            if a == "mps_L2":
                continue
            x = cell(rows, pdbs, lambda s, a=a, k=kind: f"{k}|{a}|{s}", cfg["SEEDS"])
            print("   " + line(f"{a} - mps_L2", x, base, fl))
        print("\n  -- vs the MANDATORY best_of_N at matched budget --")
        for a in arms:
            if a == "best_of_N":
                continue
            x = cell(rows, pdbs, lambda s, a=a, k=kind: f"{k}|{a}|{s}", cfg["SEEDS"])
            print("   " + line(f"{a} - best_of_N", x, bon, fl))
        print("\n  -- absolute means (built chain) --")
        for a in arms:
            x = cell(rows, pdbs, lambda s, a=a, k=kind: f"{k}|{a}|{s}", cfg["SEEDS"])
            print(f"   {a:<20} {np.nanmean(x):6.3f}  med {np.nanmedian(x):6.3f}")
        pool = np.asarray([rows[p]["ref"]["pool_mean_ORACLE"] for p in pdbs])
        pbest = np.asarray([rows[p]["ref"]["pool_best_ORACLE"] for p in pdbs])
        print(f"   {'[pool mean ORACLE]':<20} {pool.mean():6.3f}")
        print(f"   {'[pool best ORACLE]':<20} {pbest.mean():6.3f}")
        print()
    # objective-vs-structure decoupling, within block
    print("=== objective gained vs RMSD gained (F-A3 / BRIEF section 3) ===")
    for kind in cfg["kinds_A"]:
        bo = cell(rows, pdbs, lambda s, k=kind: f"{k}|best_of_N|{s}", cfg["SEEDS"], "best_std")
        br = cell(rows, pdbs, lambda s, k=kind: f"{k}|best_of_N|{s}", cfg["SEEDS"])
        for a in cfg["A_ANSATZ"] + ["untrained_L2"]:
            o = cell(rows, pdbs, lambda s, a=a, k=kind: f"{k}|{a}|{s}", cfg["SEEDS"], "best_std")
            r = cell(rows, pdbs, lambda s, a=a, k=kind: f"{k}|{a}|{s}", cfg["SEEDS"])
            m = np.isfinite(o) & np.isfinite(r) & np.isfinite(bo) & np.isfinite(br)
            if m.sum() < 5:
                continue
            rho = K.L.spearman((bo - o)[m], (br - r)[m])
            print(f"  {kind:5s} {a:<16} rho(obj gained, RMSD gained) = {rho:+.3f}  n={m.sum()}")


# ============================================================ BLOCK Q
def report_q():
    rows = K.read21("qb3_q")
    pdbs = _pdbs(rows)
    fl = folds(rows, pdbs)
    cfg = K.read21("qb3_config")
    arms = [f"opt_{o}" for o in cfg["Q_OPTS"]] + ["untrained_L2", "best_of_N"]
    print(f"BLOCK Q -- optimiser battery.  n_targets = {len(pdbs)}, seeds = {cfg['SEEDS']}, "
          f"budget = {cfg['B']}, shots = {cfg['Q_SHOTS']}, alpha = {cfg['Q_ALPHA']}")
    print("Preconditioner is the CLASSICAL FISHER of p_theta(b), not the Fubini-Study metric.\n")
    for kind in cfg["kinds_Q"]:
        print(f"=== {kind} ===")
        base = cell(rows, pdbs, lambda s, k=kind: f"{k}|opt_adam|{s}", cfg["SEEDS"])
        bon = cell(rows, pdbs, lambda s, k=kind: f"{k}|best_of_N|{s}", cfg["SEEDS"])
        print("  -- Ca-RMSD vs the deployed adam --")
        for a in arms:
            if a == "opt_adam":
                continue
            x = cell(rows, pdbs, lambda s, a=a, k=kind: f"{k}|{a}|{s}", cfg["SEEDS"])
            print("   " + line(f"{a} - adam", x, base, fl))
        print("  -- Ca-RMSD vs the MANDATORY best_of_N --")
        for a in arms:
            if a == "best_of_N":
                continue
            x = cell(rows, pdbs, lambda s, a=a, k=kind: f"{k}|{a}|{s}", cfg["SEEDS"])
            print("   " + line(f"{a} - best_of_N", x, bon, fl))
        # `best_std` is the LOWEST SINGLE STANDARDISED ENERGY the arm evaluated -- the search
        # quality column.  It is NOT the CVaR training loss; that is `cvar_first -> cvar_last`,
        # printed separately below.  Conflating the two would be exactly the kind of label drift
        # this sprint has corrected twice.
        print("  -- BEST EVALUATED ENERGY (search quality; NOT the CVaR loss) --")
        for a in arms:
            if a in ("opt_adam", "best_of_N"):
                continue
            x = cell(rows, pdbs, lambda s, a=a, k=kind: f"{k}|{a}|{s}", cfg["SEEDS"], "best_std")
            y = cell(rows, pdbs, lambda s, k=kind: f"{k}|opt_adam|{s}", cfg["SEEDS"], "best_std")
            print("   " + line(f"{a} - adam [OBJECTIVE]", x, y, fl))
        print("  -- CVaR trajectory and conditioning --")
        for a in arms:
            g = {}
            for f in ("cvar_first", "cvar_last", "gnorm_mean", "fisher_cond_med",
                      "param_disp", "iters"):
                v = [rows[p]["rows"][f"{kind}|{a}|{s}"].get(f)
                     for p in pdbs for s in cfg["SEEDS"]
                     if f"{kind}|{a}|{s}" in rows[p]["rows"]]
                v = [x for x in v if isinstance(x, (int, float)) and np.isfinite(x)]
                g[f] = float(np.mean(v)) if v else float("nan")
            print(f"   {a:<16} cvar {g['cvar_first']:+9.3f} -> {g['cvar_last']:+9.3f}  "
                  f"|g| {g['gnorm_mean']:9.3g}  F.cond {g['fisher_cond_med']:9.3g}  "
                  f"|dtheta| {g['param_disp']:7.3f}  it {g['iters']:.1f}")
        print()


# ============================================================ BLOCK E
def report_e(tag="e"):
    rows = K.read21(f"qb3_{tag}")
    pdbs = _pdbs(rows)
    fl = folds(rows, pdbs)
    cfg = K.read21("qb3_config")
    arms = cfg["E_ARMS"] if tag == "e" else ["nelder"]
    charts = cfg["E_CHARTS"]
    print(f"BLOCK E ({tag}) -- the encoding as a GAUGE question.  n_targets = {len(pdbs)}, "
          f"seeds = {cfg['SEEDS']}, budget = {cfg['B']}")
    print("Every chart below is PHYSICALLY IDENTICAL; `s` and `r` are exact gauge parameters.\n")
    for kind in cfg["kinds_E"]:
        for arm in arms:
            print(f"=== {kind} / {arm} ===")
            th = cell(rows, pdbs, lambda s, k=kind, a=arm: f"{k}|{a}|th|{s}", cfg["SEEDS"])
            print("  -- Ca-RMSD vs the `th` baseline (negative = the chart is better) --")
            for c in charts:
                if c == "th":
                    continue
                x = cell(rows, pdbs, lambda s, c=c, k=kind, a=arm: f"{k}|{a}|{c}|{s}",
                         cfg["SEEDS"])
                print("   " + line(f"{c} - th", x, th, fl))
            print("  -- WITHIN-GAUGE-ORBIT contrasts (F-E1: these must be zero) --")
            for pair in (("emb_r05", "emb_r1"), ("emb_r2", "emb_r1"),
                         ("th_s05", "th"), ("th_s2", "th")):
                a1 = cell(rows, pdbs, lambda s, c=pair[0], k=kind, a=arm: f"{k}|{a}|{c}|{s}",
                          cfg["SEEDS"])
                a2 = cell(rows, pdbs, lambda s, c=pair[1], k=kind, a=arm: f"{k}|{a}|{c}|{s}",
                          cfg["SEEDS"])
                print("   " + line(f"{pair[0]} - {pair[1]}", a1, a2, fl))
            print("  -- mechanism columns --")
            print(f"   {'chart':<12} {'RMSD':>7} {'obj':>10} {'iters':>7} {'disp_ang':>9} "
                  f"{'rad_last':>9} {'wrapfire':>9}")
            for c in charts:
                r = cell(rows, pdbs, lambda s, c=c, k=kind, a=arm: f"{k}|{a}|{c}|{s}",
                         cfg["SEEDS"])
                o = cell(rows, pdbs, lambda s, c=c, k=kind, a=arm: f"{k}|{a}|{c}|{s}",
                         cfg["SEEDS"], "best_std")
                it = cell(rows, pdbs, lambda s, c=c, k=kind, a=arm: f"{k}|{a}|{c}|{s}",
                          cfg["SEEDS"], "iters")
                dp = cell(rows, pdbs, lambda s, c=c, k=kind, a=arm: f"{k}|{a}|{c}|{s}",
                          cfg["SEEDS"], "disp_ang")
                rl = cell(rows, pdbs, lambda s, c=c, k=kind, a=arm: f"{k}|{a}|{c}|{s}",
                          cfg["SEEDS"], "radius_last_med")
                wf = cell(rows, pdbs, lambda s, c=c, k=kind, a=arm: f"{k}|{a}|{c}|{s}",
                          cfg["SEEDS"], "n_wrap_fired")
                print(f"   {c:<12} {np.nanmean(r):7.3f} {np.nanmean(o):10.3f} "
                      f"{np.nanmean(it):7.1f} {np.nanmean(dp):9.3f} "
                      f"{np.nanmean(rl):9.3f} {np.nanmean(wf):9.1f}")
            print()


# ============================================================ BLOCK S: shots x iterations
def report_s():
    """PREREG addendum §6. At fixed budget `shots * iterations = B`, so the shot count IS an 8x
    sweep of optimisation effort with the encoding, chart, dimension, probe size, ansatz, starts
    and budget all held exactly fixed. The direct test of "fewer steps helps".
    """
    rows = K.read21("qb3_s")
    pdbs = _pdbs(rows)
    fl = folds(rows, pdbs)
    cfg = K.read21("qb3_config")
    shots = cfg["S_SHOTS"]
    B = cfg["B"]
    print(f"BLOCK S -- shots x iterations at fixed budget B={B}.  n_targets={len(pdbs)}, "
          f"seeds={cfg['SEEDS']}, ansatz mps_L2, adam, alpha={cfg['S_ALPHA']}")
    print("  more shots = FEWER gradient steps.  The step count is the ONLY thing that changes.\n")
    for kind in cfg["kinds_S"]:
        if not any(f"{kind}|sh{shots[0]}|0" in rows[p]["rows"] for p in pdbs):
            continue
        cov = len([p for p in pdbs if f"{kind}|sh{shots[0]}|0" in rows[p]["rows"]])
        print(f"=== {kind} (n={cov}) ===")
        print(f"  {'shots':>6} {'steps':>6} {'CVaR last':>10} {'RMSD mean':>10} {'RMSD med':>9} "
              f"{'|dtheta|':>9}")
        for sh in shots:
            c = cell(rows, pdbs, lambda s, sh=sh, k=kind: f"{k}|sh{sh}|{s}",
                     cfg["SEEDS"], "cvar_last")
            r = cell(rows, pdbs, lambda s, sh=sh, k=kind: f"{k}|sh{sh}|{s}", cfg["SEEDS"])
            d = cell(rows, pdbs, lambda s, sh=sh, k=kind: f"{k}|sh{sh}|{s}",
                     cfg["SEEDS"], "param_disp")
            print(f"  {sh:>6} {B // sh:>6} {np.nanmean(c):>10.3f} {np.nanmean(r):>10.3f} "
                  f"{np.nanmedian(r):>9.3f} {np.nanmean(d):>9.3f}")
        b = cell(rows, pdbs, lambda s, k=kind: f"{k}|best_of_N|{s}", cfg["SEEDS"])
        print(f"  {'bestN':>6} {'0':>6} {'-':>10} {np.nanmean(b):>10.3f} {np.nanmedian(b):>9.3f}")
        base = cell(rows, pdbs, lambda s, sh=shots[-1], k=kind: f"{k}|sh{sh}|{s}", cfg["SEEDS"])
        cbase = cell(rows, pdbs, lambda s, sh=shots[-1], k=kind: f"{k}|sh{sh}|{s}",
                     cfg["SEEDS"], "cvar_last")
        print(f"\n  -- vs the FEWEST-STEPS arm (shots={shots[-1]}, {B // shots[-1]} steps) --")
        for sh in shots[:-1]:
            r = cell(rows, pdbs, lambda s, sh=sh, k=kind: f"{k}|sh{sh}|{s}", cfg["SEEDS"])
            c = cell(rows, pdbs, lambda s, sh=sh, k=kind: f"{k}|sh{sh}|{s}",
                     cfg["SEEDS"], "cvar_last")
            print("   " + line(f"sh{sh} ({B // sh} steps) [RMSD]", r, base, fl, 28))
            print("   " + line(f"sh{sh} ({B // sh} steps) [CVaR]", c, cbase, fl, 28))
        print("\n  -- every arm vs best_of_N [RMSD] --")
        for sh in shots:
            r = cell(rows, pdbs, lambda s, sh=sh, k=kind: f"{k}|sh{sh}|{s}", cfg["SEEDS"])
            print("   " + line(f"sh{sh} ({B // sh} steps)", r, b, fl, 28))
        print()


# ============================================================ latent concentration
def report_conc(kind="DIST"):
    """Does a MORE CONCENTRATED latent build a better structure?

    L4 is binding here: on this instrument a RAW rank correlation against Ca-RMSD is largely a
    TARGET-DIFFICULTY measurement, so every number is reported both ways and the partial is the
    one that counts.  The control is `pool_mean_ORACLE` -- an ORACLE difficulty proxy, labelled.

    The second panel removes the difficulty confound by construction: WITHIN each target, rank
    the eight ansatz rungs by trained latent entropy and by Ca-RMSD and correlate.  Target
    difficulty is constant inside a target, so it cannot drive that number.
    """
    rows = K.read21("qb3_a")
    cfg = K.read21("qb3_config")
    pdbs = _pdbs(rows)
    diff = np.array([rows[x]["ref"]["pool_mean_ORACLE"] for x in pdbs])
    print(f"LATENT CONCENTRATION vs STRUCTURAL OUTCOME -- {kind}, n={len(pdbs)}, "
          f"seeds={cfg['SEEDS']}")
    print(f"  {'ansatz':<16} {'raw rho(H,RMSD)':>16} {'partial':>9} "
          f"{'raw rho(maxp,RMSD)':>19} {'partial':>9}")
    out = {}
    for a in list(cfg["A_ANSATZ"]) + ["untrained_L2"]:
        r = cell(rows, pdbs, lambda s, a=a: f"{kind}|{a}|{s}", cfg["SEEDS"])
        h = cell(rows, pdbs, lambda s, a=a: f"{kind}|{a}|{s}", cfg["SEEDS"], "g1_entropy_bits")
        m = cell(rows, pdbs, lambda s, a=a: f"{kind}|{a}|{s}", cfg["SEEDS"], "g1_max_prob")
        ok = np.isfinite(r) & np.isfinite(h)
        if ok.sum() < 8:
            continue
        out[a] = {"raw_H": K.L.spearman(h[ok], r[ok]),
                  "partial_H": K.L.partial_spearman(h[ok], r[ok], diff[ok]),
                  "raw_maxp": K.L.spearman(m[ok], r[ok]),
                  "partial_maxp": K.L.partial_spearman(m[ok], r[ok], diff[ok]), "n": int(ok.sum())}
        print(f"  {a:<16} {out[a]['raw_H']:>16.3f} {out[a]['partial_H']:>9.3f} "
              f"{out[a]['raw_maxp']:>19.3f} {out[a]['partial_maxp']:>9.3f}")
    rs = []
    for x in pdbs:
        hh, rr = [], []
        for a in cfg["A_ANSATZ"]:
            v = [rows[x]["rows"][f"{kind}|{a}|{s}"] for s in cfg["SEEDS"]
                 if f"{kind}|{a}|{s}" in rows[x]["rows"]]
            if not v:
                continue
            hh.append(np.mean([q["g1_entropy_bits"] for q in v]))
            rr.append(np.mean([q["rmsd_ORACLE"] for q in v]))
        if len(hh) >= 5:
            rs.append(K.L.spearman(np.array(hh), np.array(rr)))
    b = K.boot_mean_ci(np.array(rs))
    print(f"\n  WITHIN-TARGET across the ansatz rungs (difficulty constant by construction):"
          f"\n    mean rho = {b['mean']:+.3f} [{b['ci95'][0]:+.3f}, {b['ci95'][1]:+.3f}], "
          f"n={b['n']} targets  (positive = more entropy goes with WORSE RMSD)")
    K.write21("qb3_concentration", {"kind": kind, "per_ansatz": out,
                                    "within_target_rho": b,
                                    "difficulty_control": "pool_mean_ORACLE (ORACLE)"},
              complete=True)


# ============================================================ train-H vs readout-H
def report_x(tag="q"):
    """BRIEF section 2, open question (b): a READOUT Hamiltonian different from the TRAINING one.

    CVaR is the TRAINING objective; the readout is an argmin.  Nothing forces the argmin to be
    taken under the same H the sampler trained on.  Every arm already evaluated a set of
    configurations, so re-taking the argmin under a different H costs NO budget and NO extra
    structure.  The diagonal (`xread_K` where K is the training kind) is an IDENTITY CHECK and
    must reproduce the arm's own `rmsd_ORACLE` exactly -- reported, not assumed.
    """
    rows = K.read21(f"qb3_{tag}")
    pdbs = _pdbs(rows)
    fl = folds(rows, pdbs)
    cfg = K.read21("qb3_config")
    kinds = cfg["kinds_Q"] if tag == "q" else cfg["kinds_A"]
    arms = ([f"opt_{o}" for o in cfg["Q_OPTS"]] if tag == "q" else list(cfg["A_ANSATZ"])) + \
        ["untrained_L2", "best_of_N"]
    reads = ("DIST", "LEG")
    print(f"TRAIN-H x READOUT-H, block {tag}.  n_targets = {len(pdbs)}, seeds = {cfg['SEEDS']}")
    print("Both axes are the SAME evaluated set; only which H takes the argmin changes.\n")
    # identity check first
    bad, ncmp = 0, 0
    for p in pdbs:
        for kk, r in rows[p]["rows"].items():
            tk = kk.split("|")[0]
            if tk in reads and f"xread_{tk}_rmsd" in r:
                ncmp += 1
                if abs(r[f"xread_{tk}_rmsd"] - r["rmsd_ORACLE"]) > 1e-12:
                    bad += 1
    print(f"IDENTITY GATE  xread_<train H> == the arm's own readout: "
          f"{ncmp - bad}/{ncmp} exact\n")
    print(f"  {'train H':>6} {'arm':<16} " + " ".join(f"{'read ' + r:>12}" for r in reads)
          + f" {'own':>8}")
    for tk in kinds:
        for a in arms:
            own = cell(rows, pdbs, lambda s, a=a, k=tk: f"{k}|{a}|{s}", cfg["SEEDS"])
            vs = [cell(rows, pdbs, lambda s, a=a, k=tk: f"{k}|{a}|{s}", cfg["SEEDS"],
                       f"xread_{r}_rmsd") for r in reads]
            if not np.isfinite(own).any():
                continue
            print(f"  {tk:>6} {a:<16} " + " ".join(f"{np.nanmean(v):12.3f}" for v in vs)
                  + f" {np.nanmean(own):8.3f}")
    print("\n  -- the question: does a PHYSICS-trained set, read out by the distogram, beat a "
          "distogram-trained set read out by the distogram? --")
    for tk in kinds:
        if tk == "DIST":
            continue
        for a in arms:
            x = cell(rows, pdbs, lambda s, a=a, k=tk: f"{k}|{a}|{s}", cfg["SEEDS"],
                     "xread_DIST_rmsd")
            # the DIST-trained run's DIST readout IS its own readout, by the identity the
            # gate above verifies -- fall back to it where the field was not persisted.
            y = cell(rows, pdbs, lambda s, a=a: f"DIST|{a}|{s}", cfg["SEEDS"],
                     "xread_DIST_rmsd")
            y2 = cell(rows, pdbs, lambda s, a=a: f"DIST|{a}|{s}", cfg["SEEDS"])
            y = np.where(np.isfinite(y), y, y2)
            if np.isfinite(x).sum() < 5 or np.isfinite(y).sum() < 5:
                continue
            print("   " + line(f"train {tk} - train DIST (both read DIST) [{a}]", x, y, fl, 46))
    print("\n  -- and: on a physics-trained set, is the DISTOGRAM readout better than the "
          "physics readout it trained on? --")
    for tk in kinds:
        if tk == "DIST":
            continue
        for a in arms:
            x = cell(rows, pdbs, lambda s, a=a, k=tk: f"{k}|{a}|{s}", cfg["SEEDS"],
                     "xread_DIST_rmsd")
            y = cell(rows, pdbs, lambda s, a=a, k=tk: f"{k}|{a}|{s}", cfg["SEEDS"])
            if np.isfinite(x).sum() < 5:
                continue
            print("   " + line(f"read DIST - read {tk} (train {tk}) [{a}]", x, y, fl, 46))


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "a"
    {"a": report_a, "q": report_q, "e": lambda: report_e("e"),
     "en": lambda: report_e("en"), "xa": lambda: report_x("a"),
     "xq": lambda: report_x("q"), "conc": report_conc, "s": report_s}[which]()


if __name__ == "__main__":
    main()
