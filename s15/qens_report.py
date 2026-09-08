"""SPRINT 15 / QENS -- the analysis of the replication.  Reads only checkpoints; no compute.

Everything the brief demands, emitted together and never as fields a reader can select from:
paired mean difference, bootstrap 95% CI, median, mean/sd, W/L, and one PASS/FAIL
null-calibrated concentration verdict.

TWO UNITS OF ANALYSIS, both printed, the second primary.  The original result treated
(target, seed) cells as independent units, which they are not -- seeds within a target are
repeated measurements of the same instance.  `cell` is the original convention (n = targets x
seeds).  `target` averages seeds first (n = targets) and is the honest unit for a claim about
targets.

ORACLE.  Every RMSD column is post-hoc from native coordinates and is ORACLE.

    python -m s15.qens_report
"""
from __future__ import annotations

import json
import os

import numpy as np

from s15 import qens_lib as Q

READOUTS = ("rand5_coordavg_rmsd", "rand20_coordavg_rmsd", "rand75_coordavg_rmsd",
            "setrand_coordavg_rmsd", "set_coordavg_rmsd", "set_mean_rmsd",
            "set_best_rmsd")

#: THE EMPIRICAL FALSE-POSITIVE FLOOR, measured by the RESTRAINT workstream on a comparison
#: that is ZERO BY CONSTRUCTION (maximum likelihood against a constant-width Gaussian *is*
#: least squares): it returned +0.081 [+0.014, +0.169] on one start draw and -0.003 on
#: another, and absolute arm means vary with sd 0.132 A across four draws.  A 95% interval
#: excluded zero on a null effect from the start draw alone.  Any paired mean at or below
#: this magnitude is UNRESOLVED in this machinery regardless of its interval.
FLOOR = 0.08
MDS = ("md5_coordavg_rmsd", "md20_coordavg_rmsd", "md75_coordavg_rmsd")
DIV = ("n_distinct", "draw_entropy_bits", "mean_pairwise_rmsd")


def load(tag="repl"):
    p = os.path.join(Q.RESULTS, f"qens_{tag}.json")
    with open(p) as fh:
        return json.load(fh)


def split(key):
    pdb, a, s = key.split("|")
    return pdb, float(a[1:]), int(s[1:])


def collect(cells, alpha, arm_a, arm_b, ro):
    """-> (per-cell diffs, per-target mean diffs, per-target labels, raw a, raw b)."""
    va, vb, tgt = [], [], []
    for k, c in cells.items():
        pdb, aa, s = split(k)
        if abs(aa - alpha) > 1e-12:
            continue
        if arm_a not in c or arm_b not in c:
            continue
        if ro not in c[arm_a] or ro not in c[arm_b]:
            continue
        va.append(float(c[arm_a][ro]))
        vb.append(float(c[arm_b][ro]))
        tgt.append(pdb)
    va, vb, tgt = np.array(va), np.array(vb), np.array(tgt)
    order = sorted(set(tgt.tolist()))
    ta = np.array([va[tgt == t].mean() for t in order])
    tb = np.array([vb[tgt == t].mean() for t in order])
    return va, vb, ta, tb, order


def row(name, va, vb, seed=0):
    r = Q.verdict_line(name, va, vb, seed=seed)
    return r


def table(cells, alphas, arm_a, arm_b, readouts=READOUTS, title="", subset=None):
    print()
    print("=" * 128)
    print(title or f"{arm_a} - {arm_b}   (NEGATIVE = {arm_a} better)")
    print(f"  `<=FLOOR` marks |mean| <= {FLOOR} A, the empirical false-positive floor of this "
          f"machinery (a null effect returned +0.081 [+0.014,+0.169] on one start draw).")
    print("=" * 128)
    res = {}
    if subset is not None:
        cells = {k: v for k, v in cells.items() if split(k)[0] in subset}
    for ro in readouts:
        print(f"\n  --- {ro} ---   (ORACLE post-hoc)")
        print(f"  {'unit':>7s} {'alpha':>6s} {'n':>4s} {'A':>8s} {'B':>8s} {'diff':>9s} "
              f"{'CI95':>21s} {'median':>8s} {'W/L':>9s} {'m/sd':>6s} {'verdict':>11s} "
              f"{'concentration':>20s}")
        for a in alphas:
            va, vb, ta, tb, order = collect(cells, a, arm_a, arm_b, ro)
            if va.size == 0:
                continue
            for unit, x, y in (("cell", va, vb), ("target", ta, tb)):
                r = row(f"{ro}|a{a}|{unit}", x, y)
                r["below_floor"] = bool(abs(r["mean"]) <= FLOOR)
                res[f"{ro}|a{a}|{unit}"] = r
                print(f"  {unit:>7s} {a:6.3f} {r['n']:4d} {np.mean(x):8.3f} "
                      f"{np.mean(y):8.3f} {r['mean']:+9.4f} "
                      f"[{r['ci_lo']:+8.4f},{r['ci_hi']:+8.4f}] {r['median']:+8.4f} "
                      f"{r['win']:4d}/{r['loss']:<4d} {r['mean_over_sd']:+6.2f} "
                      f"{r['verdict']:>11s} {r.get('conc_verdict',''):>20s}"
                      + ("  <=FLOOR" if r["below_floor"] else ""))
    return res


def diversity_table(cells, alphas, arms):
    print()
    print("=" * 128)
    print("DIVERSITY, beside every RMSD.  A set with one distinct member has a coordinate")
    print("average equal to that member; the set-mean law is OUT OF DOMAIN there.")
    print("=" * 128)
    out = {}
    print(f"  {'alpha':>6s} {'arm':>12s} {'n_distinct':>11s} {'median':>7s} "
          f"{'entropy(bits)':>14s} {'mean pairwise RMSD':>19s} {'collapsed cells':>16s}")
    for a in alphas:
        for arm in arms:
            nd, en, mp = [], [], []
            for k, c in cells.items():
                if abs(split(k)[1] - a) > 1e-12 or arm not in c:
                    continue
                nd.append(c[arm]["n_distinct"])
                en.append(c[arm]["draw_entropy_bits"])
                mp.append(c[arm]["mean_pairwise_rmsd"])
            if not nd:
                continue
            nd = np.array(nd, float)
            coll = float((nd <= 5).mean())
            out[f"a{a}|{arm}"] = {"n_distinct_mean": float(nd.mean()),
                                  "n_distinct_median": float(np.median(nd)),
                                  "entropy_bits": float(np.mean(en)),
                                  "mean_pairwise_rmsd": float(np.mean(mp)),
                                  "frac_cells_collapsed_le5": coll}
            print(f"  {a:6.3f} {arm:>12s} {nd.mean():11.1f} {np.median(nd):7.0f} "
                  f"{np.mean(en):14.3f} {np.mean(mp):19.3f} {coll:16.2f}"
                  + ("   <-- COLLAPSED" if coll > 0.5 else ""))
    return out


def per_target(cells, alpha, arm_a, arm_b, ro):
    va, vb, ta, tb, order = collect(cells, alpha, arm_a, arm_b, ro)
    print()
    print(f"  PER-TARGET, {ro}, alpha={alpha}, {arm_a} - {arm_b}  (ORACLE post-hoc)")
    print(f"  {'target':>8s} {'n':>4s} {arm_a:>9s} {arm_b:>9s} {'diff':>9s} {'set':>7s}")
    rows = {}
    for i, t in enumerate(order):
        m = np.array([split(k)[0] == t for k in cells if abs(split(k)[1] - alpha) < 1e-12])
        cnt = int(m.sum())
        d = ta[i] - tb[i]
        grp = "orig9" if t in Q.TARGETS9 else "new10"
        rows[t] = {"a": float(ta[i]), "b": float(tb[i]), "diff": float(d), "n": cnt,
                   "set": grp}
        print(f"  {t:>8s} {cnt:4d} {ta[i]:9.3f} {tb[i]:9.3f} {d:+9.4f} {grp:>7s}")
    return rows


def per_seed(cells, alpha, arm_a, arm_b, ro):
    ds = {}
    for k, c in cells.items():
        pdb, aa, s = split(k)
        if abs(aa - alpha) > 1e-12 or ro not in c.get(arm_a, {}):
            continue
        ds.setdefault(s, []).append(c[arm_a][ro] - c[arm_b][ro])
    print()
    print(f"  PER-SEED, {ro}, alpha={alpha}, {arm_a} - {arm_b}  (ORACLE post-hoc)")
    print(f"  {'seed':>5s} {'n':>4s} {'mean diff':>10s}")
    out = {}
    for s in sorted(ds):
        out[s] = {"n": len(ds[s]), "mean": float(np.mean(ds[s]))}
        print(f"  {s:5d} {len(ds[s]):4d} {np.mean(ds[s]):+10.4f}")
    return out


def alpha_map(cells3, readouts=("rand5_coordavg_rmsd", "rand20_coordavg_rmsd",
                                "rand75_coordavg_rmsd", "setrand_coordavg_rmsd",
                                "md5_coordavg_rmsd", "set_mean_rmsd")):
    print()
    print("=" * 128)
    print("R3. THE ALPHA MAP -- is the reported non-monotonicity real?")
    print("   (NEGATIVE = VQE better than the untrained circuit.  Target-level unit, n=19.)")
    print("=" * 128)
    alphas = sorted({split(k)[1] for k in cells3}, reverse=True)
    out = {}
    for ro in readouts:
        print(f"\n  --- {ro} ---  (ORACLE post-hoc)")
        print(f"  {'alpha':>6s} {'n_t':>4s} {'VQE':>8s} {'untr':>8s} {'diff':>9s} "
              f"{'CI95':>21s} {'W/L':>8s} {'m/sd':>6s} {'verdict':>11s} "
              f"{'distinct':>9s} {'H bits':>7s}")
        for a in alphas:
            va, vb, ta, tb, order = collect(cells3, a, "vqe", "untrained", ro)
            if not va.size:
                continue
            r = Q.verdict_line("", ta, tb)
            nd = np.mean([c["vqe"]["n_distinct"] for k, c in cells3.items()
                          if abs(split(k)[1] - a) < 1e-12])
            hb = np.mean([c["entropy_bits"] for k, c in cells3.items()
                          if abs(split(k)[1] - a) < 1e-12])
            out[f"{ro}|a{a}"] = {**r, "n_distinct": float(nd), "entropy_bits": float(hb)}
            print(f"  {a:6.3f} {r['n']:4d} {np.mean(ta):8.3f} {np.mean(tb):8.3f} "
                  f"{r['mean']:+9.4f} [{r['ci_lo']:+8.4f},{r['ci_hi']:+8.4f}] "
                  f"{r['win']:3d}/{r['loss']:<4d} {r['mean_over_sd']:+6.2f} "
                  f"{r['verdict']:>11s} {nd:9.1f} {hb:7.2f}")
    return out


def operator_law(cells, alphas):
    """Does the coordinate-average difference track the SET-MEAN difference?

    The recorded law is `d_out = 1.16 d_set_mean + 0.04 d_set_best`.  If the ensemble effect
    is a LOCATION effect it must show up as a set-mean difference of the right size; if the
    coordinate average moves without the set mean moving, the effect is shape or diversity.
    """
    print()
    print("=" * 128)
    print("DECOMPOSITION: is the coordinate-average effect explained by the SET MEAN?")
    print("  OLS of d(rand5_coordavg) on d(set_mean) and d(set_best), cell level.")
    print("=" * 128)
    out = {}
    for a in alphas:
        x1 = collect(cells, a, "vqe", "untrained", "set_mean_rmsd")
        x2 = collect(cells, a, "vqe", "untrained", "set_best_rmsd")
        y = collect(cells, a, "vqe", "untrained", "rand5_coordavg_rmsd")
        if not y[0].size:
            continue
        d1, d2, dy = x1[0] - x1[1], x2[0] - x2[1], y[0] - y[1]
        A = np.column_stack([np.ones_like(d1), d1, d2])
        beta, *_ = np.linalg.lstsq(A, dy, rcond=None)
        pred = A @ beta
        ss = 1.0 - ((dy - pred) ** 2).sum() / max(((dy - dy.mean()) ** 2).sum(), 1e-30)
        out[str(a)] = {"n": int(dy.size), "const": float(beta[0]),
                       "b_set_mean": float(beta[1]), "b_set_best": float(beta[2]),
                       "R2": float(ss),
                       "d_set_mean": float(d1.mean()), "d_set_best": float(d2.mean()),
                       "d_coordavg": float(dy.mean())}
        print(f"  alpha={a:5.3f} n={dy.size:4d}  d_coordavg={dy.mean():+.4f} = "
              f"{beta[0]:+.4f} + {beta[1]:+.4f}*d_setmean({d1.mean():+.4f}) "
              f"+ {beta[2]:+.4f}*d_setbest({d2.mean():+.4f})   R2={ss:.3f}")
    return out


def original_table():
    """The ORIGINAL E2 result recomputed from `qgeom_ens.json`, at both units of analysis.

    This is the left-hand column of the power comparison.  Its arm names are `vqe` and
    `control`; `control` is the untrained circuit, i.e. this workstream's `untrained`.
    """
    p = os.path.join(Q.RESULTS, "qgeom_ens.json")
    with open(p) as fh:
        ens = json.load(fh)["E2_generator"]
    cells = {k: {"vqe": v["vqe"], "untrained": v["control"]} for k, v in ens.items()}
    print()
    print("=" * 128)
    print("THE ORIGINAL RESULT, recomputed from qgeom_ens.json (9 targets x 3 seeds)")
    print("=" * 128)
    alphas = sorted({split(k)[1] for k in cells}, reverse=True)
    return table(cells, alphas, "vqe", "untrained",
                 readouts=("rand5_coordavg_rmsd", "rand20_coordavg_rmsd",
                           "rand75_coordavg_rmsd", "set_coordavg_rmsd",
                           "set_mean_rmsd", "set_best_rmsd"),
                 title="ORIGINAL n=27: vqe - untrained  (NEGATIVE = VQE better)")


def reference_levels(pdbs=None):
    """Where the arms sit against fixed reference points on the SAME sub-registers.

    ORACLE DIAGNOSTIC row-by-row: `certified_argmin` is the RMSD of the objective's exact
    global minimiser over all 4,096 configurations (computable only because the space is
    enumerated), `space_best` and `space_mean` are the sub-register's own best and mean.
    These say whether a "win" is a property of the OPTIMISER or of the OBJECTIVE.
    """
    pdbs = pdbs or Q.TARGETS19
    print()
    print("=" * 128)
    print("REFERENCE LEVELS on the same 12-qubit sub-registers (ORACLE DIAGNOSTIC)")
    print("=" * 128)
    rows = {}
    print(f"  {'target':>8s} {'space best':>11s} {'space mean':>11s} "
          f"{'certified argmin of the objective':>34s} {'argmin pctile':>14s}")
    for p in pdbs:
        st = Q.Struct(p)
        E = Q.hamil_sub(st)
        i = int(np.argmin(E))
        pct = float((st.rmsd < st.rmsd[i]).mean())
        rows[p] = {"space_best": float(st.rmsd.min()), "space_mean": float(st.rmsd.mean()),
                   "certified_argmin_rmsd": float(st.rmsd[i]), "argmin_percentile": pct}
        print(f"  {p:>8s} {st.rmsd.min():11.3f} {st.rmsd.mean():11.3f} "
              f"{st.rmsd[i]:34.3f} {pct:14.3f}")
    m = {k: float(np.mean([r[k] for r in rows.values()])) for k in
         ("space_best", "space_mean", "certified_argmin_rmsd", "argmin_percentile")}
    print(f"  {'MEAN':>8s} {m['space_best']:11.3f} {m['space_mean']:11.3f} "
          f"{m['certified_argmin_rmsd']:34.3f} {m['argmin_percentile']:14.3f}")
    return {"per_target": rows, "mean": m}


def location_table(cells, alphas):
    """LOCATION: the distribution-level expected RMSD, free of sampling noise. ORACLE."""
    print()
    print("=" * 128)
    print("LOCATION: E_p[RMSD] under each distribution (ORACLE post-hoc, no sampling noise)")
    print("=" * 128)
    out = {}
    print(f"  {'alpha':>6s} {'n':>4s} {'E_p[RMSD] VQE':>14s} {'E_p[RMSD] untrained':>20s} "
          f"{'diff':>9s} {'CI95':>21s} {'W/L':>9s} {'verdict':>11s}")
    for a in alphas:
        v = [c["meta"]["Ep_rmsd_vqe"] for k, c in cells.items()
             if abs(split(k)[1] - a) < 1e-12]
        u = [c["meta"]["Ep_rmsd_init"] for k, c in cells.items()
             if abs(split(k)[1] - a) < 1e-12]
        if not v:
            continue
        r = Q.verdict_line("", v, u)
        out[str(a)] = r
        print(f"  {a:6.3f} {r['n']:4d} {np.mean(v):14.3f} {np.mean(u):20.3f} "
              f"{r['mean']:+9.4f} [{r['ci_lo']:+8.4f},{r['ci_hi']:+8.4f}] "
              f"{r['win']:4d}/{r['loss']:<4d} {r['verdict']:>11s}")
    return out


def objective_axis(cells4, alphas):
    """R4: THE OBJECTIVE AXIS, reported separately from the structural one.

    The objective is rank-uniformised over the 4,096-configuration sub-register, so a value is
    directly the mean rank-percentile of what a distribution samples: **0 is the certified
    global optimum and 0.5 is a uniform random draw.**  `p_argmin` is the probability the
    final distribution places on the exact global minimiser.
    """
    print()
    print("=" * 128)
    print("R4. THE OBJECTIVE AXIS (rank-uniformised: 0 = certified optimum, 0.5 = random draw)")
    print("=" * 128)
    out = {}
    print(f"  {'alpha':>6s} {'n':>5s} {'E_p[E] VQE':>11s} {'E_p[E] untr':>12s} "
          f"{'E_p[E] anneal':>14s} {'min E anneal':>13s} {'p(argmin)':>10s} "
          f"{'mode E':>8s} {'mode RMSD':>10s}")
    for a in alphas:
        r = [c for k, c in cells4.items() if abs(split(k)[1] - a) < 1e-12]
        if not r:
            continue
        f = lambda k: float(np.mean([x[k] for x in r]))
        out[str(a)] = {k: f(k) for k in ("Ep_E_vqe", "Ep_E_init", "Ep_E_anneal_tail",
                                         "min_E_anneal", "p_argmin", "mode_E",
                                         "Ep_rmsd_vqe", "mode_rmsd")}
        out[str(a)]["n"] = len(r)
        out[str(a)]["frac_p_argmin_zero"] = float(np.mean([x["p_argmin"] < 1e-12
                                                           for x in r]))
        pr = Q.verdict_line("", [x["Ep_E_vqe"] for x in r],
                            [x["Ep_E_anneal_tail"] for x in r])
        out[str(a)]["vqe_minus_anneal_objective"] = pr
        print(f"  {a:6.3f} {len(r):5d} {f('Ep_E_vqe'):11.5f} {f('Ep_E_init'):12.5f} "
              f"{f('Ep_E_anneal_tail'):14.5f} {f('min_E_anneal'):13.6f} "
              f"{f('p_argmin'):10.5f} {f('mode_E'):8.5f} {f('mode_rmsd'):10.3f}")
    print("\n  A VQE that reached the certified optimum would show E_p[E] -> 0 and")
    print("  p(argmin) -> 1.  The classical annealer's `min E` is its best visited value.")
    return out


def objective_quality_law(cells, alphas, ref, ro="rand5_coordavg_rmsd"):
    """THE MECHANISM: the per-target win is a function of the OBJECTIVE's quality, not the circuit.

    `argmin_percentile` is where the objective's own certified global minimiser sits in the
    sub-register's ORACLE RMSD distribution -- 0 means the objective's optimum IS the best
    structure in the space, 1 means it is the worst.  It is a property of the OBJECTIVE on
    that target and has nothing to do with any optimiser.  If the paired VQE-minus-untrained
    difference tracks it, the "quantum win" is the objective's low-energy region being good,
    reached by a method that concentrates on it.
    """
    print()
    print("=" * 128)
    print("MECHANISM: per-target win vs the OBJECTIVE's own certified-argmin quality "
          "(ORACLE DIAGNOSTIC)")
    print("=" * 128)
    per = ref["per_target"]
    tg = sorted(per)
    P = np.array([per[t]["argmin_percentile"] for t in tg])
    A = np.array([per[t]["certified_argmin_rmsd"] for t in tg])
    out = {}
    print(f"  {'alpha':>6s} {'rho(diff, argmin percentile)':>30s} "
          f"{'rho(diff, argmin RMSD)':>24s}")
    for a in alphas:
        va, vb, ta, tb, order = collect(cells, a, "vqe", "untrained", ro)
        if not va.size:
            continue
        idx = {t: i for i, t in enumerate(order)}
        D = np.array([ta[idx[t]] - tb[idx[t]] for t in tg if t in idx])
        Pp = np.array([P[i] for i, t in enumerate(tg) if t in idx])
        Aa = np.array([A[i] for i, t in enumerate(tg) if t in idx])
        r1, r2 = float(Q.G.spearman(D, Pp)), float(Q.G.spearman(D, Aa))
        out[str(a)] = {"n": int(D.size), "rho_vs_argmin_percentile": r1,
                       "rho_vs_argmin_rmsd": r2}
        print(f"  {a:6.3f} {r1:+30.3f} {r2:+24.3f}")
    print("\n  A POSITIVE rho means the VQE wins on the targets where the objective's own")
    print("  optimum is GOOD and loses where it is bad -- i.e. the effect is a property of")
    print("  the OBJECTIVE, reached by a method that concentrates on it.")
    return out


def concentration_law(cells3, cells1=None, ro="rand5_coordavg_rmsd"):
    """THE COORDINATOR'S PREDICTION: the win tracks DISTINCT-CONFIGURATION COUNT, not alpha.

    QRESTRAINT found that the sign of its VQE-versus-untrained comparison depends on the
    evaluation budget, and proposed that alpha's effect here is the same phenomenon on a
    different axis -- alpha controls concentration and concentration controls how many
    distinct configurations a fixed budget yields, so a non-monotone alpha response would be
    a monotone response in DIVERSITY seen through a non-monotone map.

    Tested three ways, on the pooled alpha map:
      (a) the marginal correlation of the paired difference with log2(distinct) and log(alpha)
      (b) a two-predictor OLS, to see which survives the other
      (c) the WITHIN-ALPHA correlation with distinct count, which alpha cannot explain at all
    """
    print()
    print("=" * 128)
    print("CONCENTRATION LAW: does the win track DISTINCT-CONFIGURATION COUNT rather than alpha?")
    print("=" * 128)
    rows = []
    for k, c in (cells3 or {}).items():
        if ro not in c.get("vqe", {}):
            continue
        pdb, a, s = split(k)
        rows.append((pdb, a, s, c["vqe"][ro] - c["untrained"][ro],
                     c["vqe"]["n_distinct"], c["untrained"]["n_distinct"],
                     c.get("entropy_bits", np.nan)))
    for k, c in (cells1 or {}).items():
        if ro not in c.get("vqe", {}):
            continue
        pdb, a, s = split(k)
        rows.append((pdb, a, s, c["vqe"][ro] - c["untrained"][ro],
                     c["vqe"]["n_distinct"], c["untrained"]["n_distinct"],
                     c["meta"].get("entropy_bits", np.nan)))
    if not rows:
        print("  no cells")
        return {}
    d = np.array([r[3] for r in rows], float)
    nd = np.array([r[4] for r in rows], float)
    al = np.array([r[1] for r in rows], float)
    x1 = np.log2(np.maximum(nd, 1.0))
    x2 = np.log(al)
    sp = Q.G.spearman
    out = {"n": int(d.size),
           "rho_diff_vs_log2distinct": float(sp(d, x1)),
           "rho_diff_vs_logalpha": float(sp(d, x2)),
           "rho_log2distinct_vs_logalpha": float(sp(x1, x2))}
    A = np.column_stack([np.ones_like(d), x1, x2])
    beta, *_ = np.linalg.lstsq(A, d, rcond=None)
    pred = A @ beta
    out["ols"] = {"const": float(beta[0]), "b_log2distinct": float(beta[1]),
                  "b_logalpha": float(beta[2]),
                  "R2": float(1 - ((d - pred) ** 2).sum()
                              / max(((d - d.mean()) ** 2).sum(), 1e-30))}
    print(f"  n = {d.size} pooled cells, readout {ro}")
    print(f"  rho(diff, log2 distinct)   = {out['rho_diff_vs_log2distinct']:+.3f}")
    print(f"  rho(diff, log alpha)       = {out['rho_diff_vs_logalpha']:+.3f}")
    print(f"  rho(log2 distinct, log alpha) = "
          f"{out['rho_log2distinct_vs_logalpha']:+.3f}   [the two predictors' own overlap]")
    print(f"  OLS diff = {beta[0]:+.4f} {beta[1]:+.4f}*log2(distinct) "
          f"{beta[2]:+.4f}*log(alpha)   R2 = {out['ols']['R2']:.3f}")
    print()
    print(f"  WITHIN each alpha (which alpha cannot explain):")
    print(f"  {'alpha':>6s} {'n':>4s} {'mean distinct':>14s} {'rho(diff, distinct)':>21s} "
          f"{'mean diff':>10s}")
    within = {}
    for a in sorted(set(al.tolist()), reverse=True):
        m = al == a
        if m.sum() < 8:
            continue
        r = float(sp(d[m], nd[m]))
        within[str(a)] = {"n": int(m.sum()), "rho": r,
                          "mean_distinct": float(nd[m].mean()),
                          "mean_diff": float(d[m].mean())}
        print(f"  {a:6.3f} {int(m.sum()):4d} {nd[m].mean():14.1f} {r:+21.3f} "
              f"{d[m].mean():+10.4f}")
    out["within_alpha"] = within
    out["mean_within_alpha_rho"] = float(np.mean([v["rho"] for v in within.values()]))
    print(f"\n  mean within-alpha rho(diff, distinct) = "
          f"{out['mean_within_alpha_rho']:+.3f}")
    return out


def main():
    d = load("repl")
    cells = d.get("R1_cells", {})
    cells3 = d.get("R3_cells", {})
    alphas = sorted({split(k)[1] for k in cells}, reverse=True)
    print("=" * 128)
    print(f"QENS REPLICATION REPORT.  R1 cells: {len(cells)}   R3 cells: {len(cells3)}")
    if "R0_reproduce" in d:
        print(f"  R0 bit-level reproduction of the original 108 cells: max abs diff "
              f"{d['R0_reproduce']['max_abs_diff']:.3e} over "
              f"{d['R0_reproduce']['n_cells']} cells")
    print("=" * 128)
    out = {}
    out["original"] = original_table()
    if cells:
        out["primary"] = table(cells, alphas, "vqe", "untrained",
                               title="R1 PRIMARY: vqe - untrained  (NEGATIVE = VQE better)")
        out["matched_div"] = table(cells, alphas, "vqe", "untrained", readouts=MDS,
                                   title="R2 MATCHED DIVERSITY: vqe - untrained at equal "
                                         "distinct-configuration count")
        out["diversity"] = diversity_table(cells, alphas,
                                           ("vqe", "untrained", "uniform", "anneal",
                                            "anneal_cost", "tilt"))
        for b in ("tilt", "anneal", "anneal_cost", "uniform"):
            out[f"vs_{b}"] = table(cells, alphas, "vqe", b,
                                   readouts=("rand5_coordavg_rmsd",
                                             "rand75_coordavg_rmsd",
                                             "setrand_coordavg_rmsd", "set_mean_rmsd"),
                                   title=f"CLASSICAL CONTROL: vqe - {b}  "
                                         f"(NEGATIVE = VQE better)")
        out["subset_orig9"] = table(cells, alphas, "vqe", "untrained",
                                    readouts=("rand5_coordavg_rmsd",),
                                    title="SUBSET: the ORIGINAL nine targets only",
                                    subset=set(Q.TARGETS9))
        out["subset_new10"] = table(cells, alphas, "vqe", "untrained",
                                    readouts=("rand5_coordavg_rmsd",),
                                    title="SUBSET: the TEN NEW targets only",
                                    subset=set(Q.TARGETS10))
        out["per_target"] = {str(a): per_target(cells, a, "vqe", "untrained",
                                                "rand5_coordavg_rmsd") for a in alphas}
        out["per_seed"] = {str(a): per_seed(cells, a, "vqe", "untrained",
                                            "rand5_coordavg_rmsd") for a in alphas}
        out["location"] = location_table(cells, alphas)
        out["operator_law"] = operator_law(cells, alphas)
        out["reference"] = reference_levels()
        out["objective_quality_law"] = objective_quality_law(cells, alphas,
                                                             out["reference"])
    if cells3:
        out["alpha_map"] = alpha_map(cells3)
    if cells3 or cells:
        out["concentration_law"] = concentration_law(cells3, cells)
    cells4 = d.get("R4_cells", {})
    if cells4:
        out["objective_axis"] = objective_axis(cells4,
                                               sorted({split(k)[1] for k in cells4},
                                                      reverse=True))
    Q.ck("report", "report", out)
    print("\nwritten -> s15/results/qens_report.json")
    return out


if __name__ == "__main__":
    main()
