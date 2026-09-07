"""s21/a_report.py -- the mandatory matrix's report, in the order BRIEF section 4 requires.

    final Ca-RMSD (mean AND median) -> target-level paired difference -> CI -> W/L -> folds
    ... and ONLY THEN objective convergence, gradient geometry, CVaR behaviour, diversity.

RMSD is the endpoint; everything else is explanation.  Every averaged arm is labelled AVERAGED.
Every ORACLE quantity is labelled ORACLE and is never presented as achieved.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s18 import phys_lib as PL          # noqa: E402
from s21 import a_matrix as A           # noqa: E402

MDE = 0.084


def _load(path=None):
    p = path or os.path.join(A.RESULTS, "a_matrix.json")
    o = json.load(open(p))
    return o


def _cell(rows, h, arm, seed, key="rmsd_ORACLE"):
    out = []
    for r in rows:
        c = r["cells"].get(f"{h}|{arm}|{seed}")
        out.append(float(c[key]) if c and key in c else np.nan)
    return np.array(out, float)


def _seedmean(rows, h, arm, key="rmsd_ORACLE"):
    """The per-target mean over this arm's seeds.  Sprint 20 measured within-target ansatz-seed
    sd at 0.200 A = 2.4x MDE, so a single seed on a VARIATIONAL arm is NOT MEASURED and the
    seed mean is the unit.  Non-variational controls carry no ansatz seed (see A.seeds_for)."""
    M = np.array([_cell(rows, h, arm, s, key) for s in A.seeds_for(arm)], float)
    return np.nanmean(M, axis=0), M


def _mde(p, a, b):
    """MDE IS PER COMPARISON: 2.8016 x SE at 80% power.  The project's pooled 0.084 A constant
    is wrong by up to 84x in BOTH directions on an individual contrast, so every line carries
    its own SE and the MDE its own paired sd implies."""
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d[np.isfinite(d)]
    p["se"] = float(d.std(ddof=1) / max(np.sqrt(len(d)), 1e-12)) if len(d) > 1 else float("nan")
    p["mde"] = 2.8016 * p["se"]
    return p


def _P(a, b, folds=None, names=None):
    """Paired contrast with its OWN standard error and its OWN minimum detectable effect."""
    return _mde(PL.paired(np.asarray(a, float), np.asarray(b, float),
                          folds=folds, names=names), a, b)


def _line(p):
    ci, cf = p["ci"], p.get("ci_fold", p["ci"])
    se = p.get("se", float("nan"))
    md = p.get("mde", float("nan"))
    return (f"{p['mean']:+.3f} med {p['median']:+.3f} se {se:.3f} (MDE {md:.3f})  "
            f"iid[{ci[0]:+.3f},{ci[1]:+.3f}]  "
            f"fold[{cf[0]:+.3f},{cf[1]:+.3f}]  {p['W']}W/{p['L']}L")


def _verdict(p):
    lo, hi = p["ci"]
    if not np.isfinite(p["mean"]) or not np.isfinite(lo):
        return "n<3, no statistic"
    md = p.get("mde", MDE)
    if not np.isfinite(md):
        md = MDE
    if lo < 0 < hi:
        return "NOT MEASURED" if (hi - lo) > 2 * md else "matched"
    if abs(p["mean"]) < md:
        return "below its own MDE"
    return "BETTER" if p["mean"] < 0 else "WORSE"


def main(path=None):
    o = _load(path)
    rows = o["rows"]
    folds = np.array([r["fold"] for r in rows])
    names = [r["pdb"] for r in rows]
    n = len(rows)

    print(f"ARTEFACT complete={o.get('complete')}  n_rows={o.get('n_rows')}"
          f"  n_expected={o.get('n_expected')}  n_complete_rows={o.get('n_complete_rows')}"
          f"  cfg={o.get('cfg_hash')}")
    print(f"CONFIG budget={A.BUDGET} candidate evaluations/arm  shots={A.SHOTS}  "
          f"alpha={A.ALPHA}  seeds={list(A.SEEDS)}  norm={A.NORM}")
    print(f"n = {n} targets.  BASIS: BUILT CHAIN (single structure from continuous torsions "
          f"via build_ca_exact).")
    print("MDE IS PER COMPARISON: every line prints its own SE and its own 2.8016 x SE minimum")
    print("detectable effect.  The project's pooled 0.084 A constant is NOT used as a threshold")
    print("here -- it is wrong by up to 84x in both directions on an individual contrast.")
    if n < 40:
        print(f"POWER WARNING: n={n} is a COMPUTE cut (52 ms/AMBER point under contention vs a")
        print("6 ms nominal).  Per-cell configuration is intact; only power is reduced.  Read")
        print("every zero-spanning interval here as NOT MEASURED, never as 'matched'.")

    d0 = np.array([r["ORACLE_pool_mean"] for r in rows], float)
    d1 = np.array([r["ORACLE_pool_best"] for r in rows], float)
    print(f"\nORACLE references (never achieved): retrieval pool mean {d0.mean():.3f}, "
          f"pool best member {d1.mean():.3f}")

    # =================================================================== 1. RMSD
    print(f"\n{'='*100}")
    print("1.  FINAL Ca-RMSD -- THE ENDPOINT.  Mean and median over targets, per Hamiltonian "
          "x arm.")
    print("    Readout = argmin of the TRAINING Hamiltonian over everything seen "
          "(single structure, the pillar's own operator).")
    print(f"{'='*100}")
    hdr = f"  {'Hamiltonian':<14}" + "".join(f"{a:>17}" for a in A.ARMS)
    print(hdr)
    tab = {}
    for h in A.HAMS:
        line = f"  {h:<14}"
        for arm in A.ARMS:
            v, _M = _seedmean(rows, h, arm)
            tab[(h, arm)] = v
            line += f"{np.nanmean(v):>9.3f}/{np.nanmedian(v):<7.3f}"
        print(line)
    print(f"  (mean / median.  Variational arms {A.VARIATIONAL} are the mean over "
          f"{len(A.SEEDS)} ansatz seeds; non-variational controls over {len(A.CTRL_SEEDS)} "
          f"-- they carry no ansatz seed, a DECLARED reduction.)")

    # =================================================================== 2. paired
    print(f"\n{'='*100}")
    print("2.  TARGET-LEVEL PAIRED DIFFERENCES, CI, W/L, FOLDS -- in that order, before "
          "any diagnostic.")
    print(f"{'='*100}")
    ref = "Distance"
    print(f"\n  2a.  Every Hamiltonian's VQE arm minus the STRUCTURAL-ONLY reference "
          f"({ref}, VQE).  Negative = physics helps.")
    base = tab[(ref, "vqe")]
    for h in A.HAMS:
        if h == ref:
            continue
        p = _P(tab[(h, "vqe")], base, folds, names)
        print(f"    {h:<14}{_line(p)}   -> {_verdict(p)}")

    print(f"\n  2b.  MANDATORY CONTROL -- VQE minus BEST-OF-N FROM THE UNTRAINED CIRCUIT, "
          f"matched shots and budget.")
    print("       (theta never stepped.  This is NOT an initialisation mean.)")
    for h in A.HAMS:
        p = _P(tab[(h, "vqe")], tab[(h, "vqe_untrained")], folds, names)
        print(f"    {h:<14}{_line(p)}   -> {_verdict(p)}")

    print(f"\n  2c.  MATCHED CLASSICAL CONTROLS at identical budget, per Hamiltonian.")
    for ctl in ("best_of_N", "metro"):
        print(f"    VQE minus {ctl}:")
        for h in A.HAMS:
            p = _P(tab[(h, "vqe")], tab[(h, ctl)], folds, names)
            print(f"      {h:<14}{_line(p)}   -> {_verdict(p)}")

    print(f"\n  2d.  ZERO-INFORMATION reference (constant alpha-helix + matched torsion noise, "
          f"same Hamiltonian, same budget).")
    print("       As a REFERENCE this is a null; as a GATE it is the worst arm measured "
          "(BRIEF section 7 rule 4).")
    for h in A.HAMS:
        p = _P(tab[(h, "vqe")], tab[(h, "helix")], folds, names)
        print(f"    {h:<14}{_line(p)}   -> {_verdict(p)}")

    print(f"\n  2e.  PER-FOLD means of the primary contrast (each Hamiltonian's VQE vs "
          f"{ref} VQE):")
    for h in A.HAMS:
        if h == ref:
            continue
        p = _P(tab[(h, "vqe")], base, folds=folds)
        pf = p.get("per_fold", {})
        print(f"    {h:<14}" + "  ".join(f"f{k}:{v:+.3f}" for k, v in sorted(pf.items())))

    # =================================================================== 3. seeds
    print(f"\n{'='*100}")
    print("3.  ANSATZ-SEED SENSITIVITY -- the reason 4 seeds is the minimum.")
    print(f"{'='*100}")
    for h in A.HAMS:
        _v, M = _seedmean(rows, h, "vqe")
        _vu, MU = _seedmean(rows, h, "vqe_untrained")
        print(f"  {h:<14}within-target seed sd {np.nanmean(np.nanstd(M, axis=0)):.3f} A "
              f"(untrained {np.nanmean(np.nanstd(MU, axis=0)):.3f}); "
              f"best-of-4-seeds {np.nanmean(np.nanmin(M, axis=0)):.3f} vs "
              f"mean {np.nanmean(M):.3f}")

    # =================================================================== 4. readouts
    print(f"\n{'='*100}")
    print("4.  THE READOUT COLUMN.  *** A HAMILTONIAN'S SIGN DEPENDS ON THE READOUT. ***")
    print("    Coordinator, n=126 pool-restricted: Legacy BEATS its matched control at")
    print("    tail_member (-0.41 [-0.56,-0.26]) and LOSES to it at tail_medoid (+0.30) and")
    print("    tail_avg (+0.33).  Both are true.  So no row is quoted without its readout, and")
    print("    'which Hamiltonian is best' is reported PER READOUT, never unqualified.")
    print("    Mechanism is ERROR COHERENCE, not diversity: averaging cancels i.i.d. error and")
    print("    PRESERVES systematic error, and Legacy's tail carries a coherent compactness")
    print("    bias (Sprint 20 L2c).  The diversity account was refuted by its own sign control.")
    print(f"{'='*100}")
    for a in A.TAIL_ALPHAS:
        print(f"\n  alpha = {a:g}   (mean Ca-RMSD, A; every cell is the {len(A.SEEDS)}-seed mean "
              f"of the VQE arm)")
        print(f"    {'Hamiltonian':<14}" +
              "".join(f"{r:>12}" for r in ("tail_min", "tail_member", "tail_medoid",
                                           "tail_avg")) +
              "   |" + "".join(f"{'rand_'+r:>13}" for r in ("member", "medoid", "avg")))
        for h in A.HAMS:
            line = f"    {h:<14}"
            for pre in ("tailmin", "tailmember", "tailmedoid", "tailavg"):
                v, _ = _seedmean(rows, h, "vqe", f"{pre}{a:g}")
                line += f"{np.nanmean(v):>12.3f}"
            line += "   |"
            for pre in ("randmember", "randmedoid", "randavg"):
                v, _ = _seedmean(rows, h, "vqe", f"{pre}{a:g}")
                line += f"{np.nanmean(v):>13.3f}"
            print(line)
        print(f"\n    PAIRED vs the MATCHED-COUNT RANDOM control of the SAME readout "
              f"(negative = the Hamiltonian's tail carries information over chance):")
        for pre, rnd, lab in (("tailmember", "randmember", "member (single structure)"),
                              ("tailmedoid", "randmedoid", "medoid (single structure)"),
                              ("tailavg", "randavg", "AVERAGED (point cloud)")):
            print(f"      readout = {lab}")
            for h in A.HAMS:
                x, _ = _seedmean(rows, h, "vqe", f"{pre}{a:g}")
                y, _ = _seedmean(rows, h, "vqe", f"{rnd}{a:g}")
                p = _P(x, y, folds, names)
                print(f"        {h:<14}{_line(p)}   -> {_verdict(p)}")

    print(f"\n  IDENTITY CHECK (a check, never a discovery): tail_min == argmin over everything")
    print("  seen, because the alpha-tail of the TRAINING H contains that H's own minimum.")
    idv = []
    for h in A.HAMS:
        for s in A.seeds_for("vqe"):
            idv.append(_cell(rows, h, "vqe", s, "identity_tailmin_vs_argmin"))
    idv = np.concatenate(idv)
    print(f"    max |tail_min - argmin| = {np.nanmax(idv):.3e} over "
          f"{int(np.sum(np.isfinite(idv)))} cells")

    # =================================================================== 5. cross
    print(f"\n{'='*100}")
    print(f"5.  CROSS-READOUT, all {len(A.SEEDS)} variational seeds.  rows = TRAINING H "
          "(defines the candidate set), cols = READOUT H (argmin inside it).  Single structure.")
    print("    COORDINATOR PRIORITY: the pool-restricted, same-H, order-based case is closed")
    print("    analytically and verified at tail_min == argmin, max difference exactly 0, so the")
    print("    only cells that can be informative are those breaking that scope -- readout H")
    print("    DIFFERENT from training H.  The diagonal is the argmin of the training H.")
    print(f"{'='*100}")
    print(f"  {'train|read':<14}" + "".join(f"{h[:11]:>13}" for h in A.HAMS) + f"{'argmin(tr)':>13}")
    X = {}
    for h in A.HAMS:
        line = f"  {h:<14}"
        for h2 in A.HAMS:
            v, _ = _seedmean(rows, h, "vqe", f"x_global|{h2}")
            X[(h, h2)] = v
            line += f"{np.nanmean(v):>13.3f}"
        print(line + f"{np.nanmean(tab[(h,'vqe')]):>13.3f}")
    dg = max(abs(np.nanmean(X[(h, h)]) - np.nanmean(tab[(h, "vqe")])) for h in A.HAMS)
    print(f"  IDENTITY CHECK (a check, not a discovery): max |diag - argmin(train)| = {dg:.3e}")
    print("\n  BEST OFF-DIAGONAL vs BEST DIAGONAL, paired over targets.  Both cells are chosen")
    print("  ON this instrument -- a max over 42 and 7 cells -- so this is a CEILING, not a method.")
    bd = min(A.HAMS, key=lambda h: np.nanmean(X[(h, h)]))
    offs = [(h, h2) for h in A.HAMS for h2 in A.HAMS if h != h2]
    bo = min(offs, key=lambda q: np.nanmean(X[q]))
    p = _P(X[bo], X[(bd, bd)], folds, names)
    print(f"    {bo[0]}>{bo[1]} vs {bd}>{bd}: {_line(p)}   -> {_verdict(p)}")
    sb, _ = _seedmean(rows, "Distance", "vqe", "ORACLE_seen_best")
    sm, _ = _seedmean(rows, "Distance", "vqe", "ORACLE_seen_mean")
    print(f"  ORACLE ceilings on the Distance VQE's OWN generated set (never achieved): "
          f"best {np.nanmean(sb):.3f}, mean {np.nanmean(sm):.3f}")

    print("\n  ORACLE MARGINAL rank skill of each component against true RMSD, on the VQE's own")
    print("  generated candidates.  POSITIVE = lower energy goes with lower RMSD = correct order.")
    print("  A MARGINAL IS NOT A CONTRIBUTION: Legacy is strongly rank-correlated with the")
    print("  deployed distogram, and its PARTIAL given that distogram is -0.008 fold[-0.062,")
    print("  +0.068], 25/42 negative (findings S1.3).  Read the Legacy row as a marginal ONLY.")
    for k in A.COMPONENTS:
        v = _cell(rows, "Distance", "vqe", A.SEEDS[0], f"ORACLE_rho|{k}")
        print(f"    {k:<6}{np.nanmean(v):+.4f} (median {np.nanmedian(v):+.4f}, "
              f"{int(np.nansum(v>0))}/{np.sum(np.isfinite(v))} positive)")

    # =================================================================== 6. CVaR / opt
    print(f"\n{'='*100}")
    print("6.  EXPLANATION ONLY -- objective convergence, gradient geometry, CVaR behaviour.")
    print(f"{'='*100}")
    print(f"  {'Hamiltonian':<14}{'cvar_first':>12}{'cvar_last':>12}{'d_cvar':>10}"
          f"{'|grad|':>11}{'grad sd':>10}{'ESS/shots':>11}{'cos(g,g_a1)':>13}{'|dtheta|':>10}")
    for h in A.HAMS:
        vals = {}
        for k in ("cvar_first", "cvar_last", "gnorm_mean", "gnorm_sd", "ess_frac",
                  "cos_grad_vs_alpha1", "param_disp"):
            v, _ = _seedmean(rows, h, "vqe", k)
            vals[k] = np.nanmean(v)
        print(f"  {h:<14}{vals['cvar_first']:>12.3f}{vals['cvar_last']:>12.3f}"
              f"{vals['cvar_last']-vals['cvar_first']:>10.3f}{vals['gnorm_mean']:>11.4f}"
              f"{vals['gnorm_sd']:>10.4f}{vals['ess_frac']:>11.3f}"
              f"{vals['cos_grad_vs_alpha1']:>13.3f}{vals['param_disp']:>10.3f}")
    print("  d_cvar < 0 means the objective FELL.  Read it against column 1 of section 1:")
    print("  'optimised the objective' and 'got a better structure' are different claims.")

    print("\n  OBJECTIVE-vs-STRUCTURE: paired (objective gained) vs (RMSD gained), "
          "trained minus untrained:")
    for h in A.HAMS:
        cf, _ = _seedmean(rows, h, "vqe", "cvar_first")
        cl, _ = _seedmean(rows, h, "vqe", "cvar_last")
        dq = cl - cf
        dr = tab[(h, "vqe")] - tab[(h, "vqe_untrained")]
        m = np.isfinite(dq) & np.isfinite(dr)
        rho = PL.spearman(dq[m], dr[m]) if m.sum() > 3 else float("nan")
        print(f"    {h:<14}rho(objective gained, RMSD gained) = {rho:+.3f}  "
              f"(n={int(m.sum())})")

    print("\n  NON-FINITE / SENTINEL COUNTS -- a guard that never fires is not evidence:")
    for h in A.HAMS:
        v = np.nansum([_cell(rows, h, "vqe", s, "n_nonfinite") for s in A.seeds_for("vqe")])
        u = np.nansum([_cell(rows, h, "vqe", s, "used") for s in A.seeds_for("vqe")])
        print(f"    {h:<14}{int(v)} sentinel substitutions over {int(u)} evaluations "
              f"({100.0*v/max(u,1):.2f}%)")

    # =================================================================== 7. audits
    print(f"\n{'='*100}")
    print("7.  NORMALISATION AUDIT -- the two forms that were NOT declared, on the HYBRID "
          "cells (seed %d)." % A.SEEDS[0])
    print("    Single-component cells are EXACTLY invariant (gate G-A3), so only hybrids "
          "can move.")
    print(f"{'='*100}")
    #: DEFECT FOUND IN MY OWN TABLE (coordinator fork review, 2026-09-07).  The audit arms are
    #: run at SEED 0 only, but the first published version of this table put the DECLARED form's
    #: FOUR-SEED MEAN beside them -- an unmatched comparison across the seed dimension, on an
    #: instrument whose seed sd is 0.51-1.06 A.  The declared column is now the SEED-0 value, so
    #: all three forms are read on the identical circuit.  The 4-seed mean is printed beside it,
    #: labelled, and is NOT the comparator.
    print(f"  {'Hamiltonian':<14}{'declared@s0':>13}" +
          "".join(f"{a + '@s0':>13}" for a in A.AUDIT_NORMS) +
          f"{'WINNER':>14}{'| declared 4-seed':>19}")
    best = []
    for h in ("Leg+Amb", "Dist+Leg", "Dist+Amb", "Dist+Leg+Amb"):
        vals = {"declared": float(np.nanmean(_cell(rows, h, "vqe", A.SEEDS[0])))}
        for nm in A.AUDIT_NORMS:
            v = np.array([float(r["cells"].get(f"AUD{nm}|{h}|vqe|{A.SEEDS[0]}", {})
                                .get("rmsd_ORACLE", np.nan)) for r in rows])
            vals[nm] = float(np.nanmean(v))
        w = min(vals, key=vals.get)
        best.append(vals[w])
        line = f"  {h:<14}" + "".join(f"{vals[k]:>13.3f}"
                                      for k in ("declared",) + tuple(A.AUDIT_NORMS))
        print(line + f"{w:>14}" + f"{np.nanmean(tab[(h,'vqe')]):>19.3f}")
    print("  If a non-declared form wins, the declared choice is recorded as WRONG on that cell.")
    d0 = float(np.nanmean(_cell(rows, "Distance", "vqe", A.SEEDS[0])))
    print(f"\n  THE FORK IS NOT LOAD-BEARING, and this is the strong form of the answer:")
    print(f"    best normalisation PER CELL = {[round(b, 3) for b in best]}, min {min(best):.3f}")
    print(f"    Distance-only at the same seed  = {d0:.3f}  (it carries NO normalisation at all)")
    print(f"    margin {min(best) - d0:+.3f} A -- Distance still wins, but by far less than an")
    print(f"    unmatched table would have suggested, and the margin is deep inside the")
    print(f"    0.51-1.06 A seed noise, so it is NOT MEASURED in either direction.")

    # =================================================================== 8. gates
    print(f"\n{'='*100}")
    print("8.  GATES -- how many times each FIRED (a gate that never fires passed vacuously).")
    print(f"{'='*100}")
    rel = np.array([r["G_A1_rel"] for r in rows], float)
    nk = int(np.sum([r["G_A1_n"] for r in rows]))
    print(f"  G-A1 AmberSP bit-exact vs core.amber.refine_coords(k=0, steps=-1): "
          f"max rel {np.nanmax(rel):.3e} over {nk} comparisons (FIRED {nk}x)")
    gp = os.path.join(A.RESULTS, "a_matrix_gate.json")
    if os.path.exists(gp):
        g = json.load(open(gp))
        print(f"  G-A2 conditioning order-preserving: {g['G_A2_max_order_breaks']} breaks "
              f"over {g['n_compare_G_A2']} ranks (FIRED {g['n_compare_G_A2']}x)")
        print(f"  G-A3 single-component argmin invariant to normalisation: "
              f"{g['G_A3_max_argmin_breaks']} breaks over {g['n_compare_G_A3']} checks")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
