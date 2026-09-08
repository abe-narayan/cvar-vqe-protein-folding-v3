"""s19/agentC_report.py -- every AGENT C table, from complete artefacts only.

Nothing in this module computes a scientific quantity; it reads `s19/results/*.json` (refusing a
partial through `s18.phys_lib.read_complete`) and prints the pre-registered comparisons with the
statistics `s18.phys_lib.paired` produces, unchanged, so Sprint-18 and Sprint-19 numbers come out
of literally the same estimator.

    python -m s19.agentC_report kv        # Q1, the ambiguity decomposition
    python -m s19.agentC_report pareto    # Q2, the AMBER frontier
    python -m s19.agentC_report reject    # Q3, hard rejection
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

from s18 import phys_lib as PL                                    # noqa: E402
from s19 import agentC_lib as CL                                  # noqa: E402

RES = CL.RESULTS


def _load(name, need=126):
    return PL.read_complete(os.path.join(RES, name), need=need)


def _fmt(p, unit="A"):
    ci = p.get("ci_fold", p["ci"])
    return (f"{p['mean']:+.4f} [{ci[0]:+.4f}, {ci[1]:+.4f}]  med {p['median']:+.4f}  "
            f"{p['W']}W/{p['L']}L  n={p['n']}")


# ======================================================================= Q1
def kv(fname="agentC_kv.json"):
    o = _load(fname)
    rows = o["rows"]
    pdbs = [r["pdb"] for r in rows]
    folds = np.array([r["fold"] for r in rows])
    arms = list(rows[0]["arms"].keys())

    def col(arm, key):
        return np.array([r["arms"][arm].get(key, np.nan) for r in rows], float)

    print("=" * 100)
    print("Q1  THE AMBIGUITY DECOMPOSITION OF GATE DAMAGE     n = %d, m = %d of K = %d"
          % (len(rows), rows[0]["m"], rows[0]["K"]))
    print("    pre-registration: s19/PREREG_C.md section 3;  identity verified in gate GC0")
    print("=" * 100)

    print("\n-- A.  THE IDENTITY HOLDS ON EVERY ARM (EXACT, not a finding)")
    worst = max(float(np.nanmax(np.abs(col(a, "kv_resid")))) for a in arms)
    print(f"   max |readout^2 - (E_mem^2 - D^2)| over {len(arms)} arms x {len(rows)} targets "
          f"= {worst:.3e} A^2")

    print("\n-- B.  THE LEVELS.  readout = the emitted coordinate average's Ca-RMSD (the metric).")
    hdr = f"   {'arm':<16}{'readout':>9}{'E_mem':>9}{'D':>9}{'FRAME':>9}{'<d_reb>':>9}" \
          f"{'d_win_mean':>12}{'d_win_best':>11}{'err_cos':>9}{'shared':>8}"
    print(hdr)
    for a in arms:
        print(f"   {a:<16}{np.nanmean(col(a,'readout')):>9.3f}{np.nanmean(col(a,'E_mem')):>9.3f}"
              f"{np.nanmean(col(a,'D')):>9.3f}{np.nanmean(col(a,'FRAME')):>9.3f}"
              f"{np.nanmean(col(a,'d_reb_mean')):>9.3f}{np.nanmean(col(a,'d_mean')):>12.3f}"
              f"{np.nanmean(col(a,'d_best')):>11.3f}{np.nanmean(col(a,'err_cos')):>9.3f}"
              f"{np.nanmean(col(a,'shared_frac')):>8.3f}")

    print("\n-- C.  THE SPRINT-18 RESULT REPRODUCES (readout, gate vs MATCHED-RANDOM, same count)")
    for a in arms:
        if a in ("none", "rand"):
            continue
        p = PL.paired(col(a, "readout"), col("rand", "readout"), folds=folds, names=pdbs)
        print(f"   {a:<16} vs rand   {_fmt(p)}")
    p = PL.paired(col("rand", "readout"), col("none", "readout"), folds=folds, names=pdbs)
    print(f"   {'rand':<16} vs none   {_fmt(p)}      <- the cost of halving alone")

    print("\n-- D.  THE PRE-REGISTERED PRIMARY: the three channels, gate vs matched-random.")
    print("       Delta(readout^2) = Delta(E_mem^2) - Delta(D^2)   is EXACT, term by term.")
    print(f"   {'arm':<16}{'D(readout^2)':>14}{'D(E_mem^2)':>13}{'-D(D^2)':>11}"
          f"{'D(<d_reb^2>)':>14}{'D(FRAME^2)':>13}{'S_D':>8}{'S_F':>8}{'S_q':>8}")
    prim = {}
    for a in arms:
        if a in ("none", "rand"):
            continue
        dr = col(a, "readout2") - col("rand", "readout2")
        de = col(a, "E_mem2") - col("rand", "E_mem2")
        dd = col("rand", "D2") - col(a, "D2")            # -Delta(D^2): positive = diversity LOST
        dq = col(a, "d_reb2_mean") - col("rand", "d_reb2_mean")
        df = col(a, "FRAME2") - col("rand", "FRAME2")
        tot = float(np.mean(dr))
        prim[a] = dict(dr=dr, de=de, dd=dd, dq=dq, df=df)
        print(f"   {a:<16}{np.mean(dr):>14.4f}{np.mean(de):>13.4f}{np.mean(dd):>11.4f}"
              f"{np.mean(dq):>14.4f}{np.mean(df):>13.4f}"
              f"{np.mean(dd)/tot:>8.2f}{np.mean(df)/tot:>8.2f}{np.mean(dq)/tot:>8.2f}"
              if abs(tot) > 1e-12 else "")
    print("       S_D = share of the damage carried by DIVERSITY LOSS, S_F by FRAME, "
          "S_q by member quality.")

    print("\n-- E.  EACH CHANNEL WITH ITS OWN CI (paired, fold-clustered).  ZERO is the null.")
    z = np.zeros(len(rows))
    for a in arms:
        if a in ("none", "rand"):
            continue
        P = prim[a]
        print(f"   {a}")
        for lab, v in (("Delta(readout^2)", P["dr"]), ("-Delta(D^2)   [diversity lost]", P["dd"]),
                       ("Delta(E_mem^2)", P["de"]), ("Delta(FRAME^2)", P["df"]),
                       ("Delta(<d_reb^2>)", P["dq"])):
            print(f"      {lab:<32}{_fmt(PL.paired(v, z, folds=folds, names=pdbs))}")

    print("\n-- F.  THE OPERATOR LAW, RECOMPUTED ON MY OWN ARMS (not imported).")
    print("       law: d_out = 1.16 * d_set_mean + 0.04 * d_set_best   (the ledger's fit)")
    pred_none = 1.16 * col("none", "d_mean") + 0.04 * col("none", "d_best")
    for a in arms:
        if a == "none":
            continue
        pred = 1.16 * col(a, "d_mean") + 0.04 * col(a, "d_best")
        obs = col(a, "readout") - col("none", "readout")
        prd = pred - pred_none
        miss = obs - prd
        print(f"   {a:<16} predicts {np.mean(prd):+.4f}   observes {np.mean(obs):+.4f}   "
              f"miss {_fmt(PL.paired(miss, z, folds=folds, names=pdbs))}")

    print("\n-- G.  SECTION 3.5: THE DIVERSITY-PRESERVING GATES (H-C2).")
    for a in ("legacy_spread", "legacy_clust"):
        if a not in arms:
            continue
        for ref in ("legacy", "rand", "none"):
            p = PL.paired(col(a, "readout"), col(ref, "readout"), folds=folds, names=pdbs)
            print(f"   {a:<16} vs {ref:<8} {_fmt(p)}")
        print(f"      D: {np.nanmean(col(a,'D')):.3f}  vs legacy {np.nanmean(col('legacy','D')):.3f}"
              f"  vs rand {np.nanmean(col('rand','D')):.3f}"
              f"   |  E_mem: {np.nanmean(col(a,'E_mem')):.3f} vs "
              f"{np.nanmean(col('legacy','E_mem')):.3f} / {np.nanmean(col('rand','E_mem')):.3f}")
    return o


# ======================================================================= Q3
def reject(fname="agentC_reject.json"):
    o = _load(fname)
    rows = o["rows"]
    pdbs = [r["pdb"] for r in rows]
    folds = np.array([r["fold"] for r in rows])
    arms = list(rows[0]["arms"].keys())

    def col(arm, key="readout"):
        return np.array([r["arms"][arm].get(key, np.nan) for r in rows], float)

    print("=" * 100)
    print("Q3  HARD REJECTION -- can Legacy reject a class of IMPOSSIBLE candidate?")
    print("    pre-registration: s19/PREREG_C.md section 5.  PRIMARY r = 10.")
    print("=" * 100)
    rs = sorted({int(a.split("@")[1]) for a in arms if "@" in a})
    for r in rs:
        tag = "  <-- PRE-REGISTERED PRIMARY" if r == 10 else ""
        print(f"\n-- r = {r}{tag}")
        for a in [x for x in arms if x.endswith(f"@{r}")]:
            base = a.split("@")[0]
            if base in ("rand", "divkeep"):
                continue
            for ref in ("rand", "divkeep", "none"):
                rk = f"{ref}@{r}" if ref != "none" else "none"
                if rk not in arms:
                    continue
                p = PL.paired(col(a), col(rk), folds=folds, names=pdbs)
                print(f"   {base:<12} vs {ref:<9} {_fmt(p)}")
        for ctrl in ("divkeep", "rand"):
            k = f"{ctrl}@{r}"
            if k in arms:
                p = PL.paired(col(k), col("none"), folds=folds, names=pdbs)
                print(f"   {ctrl:<12} vs {'none':<9} {_fmt(p)}   [CONTROL, score-free]")
        print(f"   {'levels':<12} " + "  ".join(
            f"{a}={np.nanmean(col(a)):.3f}" for a in arms if a.endswith(f"@{r}") or a == "none"))
    print("\n-- DOSE-RESPONSE: the cost of rejecting r, by rule (mean Ca-RMSD, n = 126)")
    bases = sorted({a.split("@")[0] for a in arms if "@" in a})
    print(f"   {'rule':<12}" + "".join(f"{('r='+str(r)):>10}" for r in rs))
    print(f"   {'none':<12}" + "".join(f"{np.nanmean(col('none')):>10.4f}" for _ in rs))
    for b in bases:
        print(f"   {b:<12}" + "".join(
            f"{np.nanmean(col(f'{b}@{r}')):>10.4f}" if f"{b}@{r}" in arms else f"{'--':>10}"
            for r in rs))
    print("\n-- HOW MANY CANDIDATES ARE ACTUALLY IMPOSSIBLE?  (min heavy-atom distance, "
          "ideal-geometry rebuild)")
    for k, lab in (("n_impossible_2p0", "< 2.0 A"), ("n_impossible_2p6", "< 2.6 A")):
        v = np.array([r[k] for r in rows], float)
        print(f"   candidates per 75 with min_heavy {lab}: mean {v.mean():.2f}  "
              f"median {np.median(v):.1f}  max {v.max():.0f}  "
              f"targets with none {(v == 0).sum()}/{len(v)}")
    mh = np.array([r["min_heavy_pool_min"] for r in rows], float)
    print(f"   worst min_heavy in the pool: mean {mh.mean():.3f} A, min {mh.min():.3f} A")
    return o


# ======================================================================= Q2
def pareto(fname="agentC_pareto.json"):
    o = _load(fname)
    rows = o["rows"]
    pdbs = [r["pdb"] for r in rows]
    folds = np.array([r["fold"] for r in rows])
    arms = [a for a in rows[0]["arms"].keys()]
    VAL = ("n_clash_2A", "n_clash_2p6A", "min_heavy", "bond_strain", "angle_strain",
           "rama_favoured", "rama_outlier", "cis_frac", "chirality_L_frac")

    def col(arm, key):
        return np.array([r["arms"][arm].get(key, np.nan) for r in rows], float)

    def vcol(arm, key):
        return np.array([r["arms"][arm].get("valid", {}).get(key, np.nan) for r in rows], float)

    print("=" * 100)
    print("Q2  AMBER AS A CONSTRAINED TERMINAL OPERATOR -- TWO AXES, NEVER ONE SCALAR")
    print("    pre-registration: s19/PREREG_C.md section 4.  n = %d" % len(rows))
    print("=" * 100)

    ok = {a: (np.isfinite(col(a, "energy")) & (col(a, "energy") <= 1000.0)
              & np.array([bool(r["arms"][a].get("converged", False)) for r in rows]))
          for a in arms}
    allok = np.ones(len(rows), bool)
    for a in arms:
        allok &= ok[a]
    print(f"\n-- CONVERGENCE GATE (declared in PREREG_C.md section 2, before use): "
          f"finite E <= 1000 kcal/mol, restraint off.")
    for a in arms:
        ex = [pdbs[i] for i in np.flatnonzero(~ok[a])]
        print(f"   {a:<18} excluded {len(ex):>2}  {' '.join(ex)}")
    exall = [pdbs[i] for i in np.flatnonzero(~allok)]
    print(f"   {'ALL-ARM GATE':<18} excluded {len(exall):>2}  {' '.join(exall)}")

    inp = col(arms[0], "input_rmsd")
    print(f"\n-- INPUT (identical for every arm): mean Ca-RMSD {np.nanmean(inp):.4f} A, "
          f"gated {np.nanmean(inp[allok]):.4f} A")
    ivalid = {k: np.array([rows[i]["input_valid"][k] for i in range(len(rows))], float)
              for k in VAL}
    print("   input validity: " + "  ".join(f"{k}={np.nanmean(ivalid[k][allok]):.3f}" for k in VAL))

    print("\n-- AXIS 1: ACCURACY.  Every arm against ITS OWN GATED INPUT, all-arm gate.")
    print(f"   {'arm':<18}{'Ca-RMSD':>9}{'delta vs input (fold CI)':>44}{'ca_disp':>9}")
    for a in arms:
        g = allok & ok[a]
        p = PL.paired(col(a, "rmsd")[g], inp[g], folds=folds[g],
                      names=[pdbs[i] for i in np.flatnonzero(g)])
        print(f"   {a:<18}{np.nanmean(col(a,'rmsd')[g]):>9.3f}   {_fmt(p)}"
              f"   {np.nanmean(col(a,'ca_disp')[g]):>6.3f}")

    print("\n-- AXIS 2: PHYSICAL VALIDITY.  Same structures, same gate.  Never fused into a scalar.")
    print(f"   {'arm':<18}" + "".join(f"{k[:11]:>12}" for k in VAL))
    print(f"   {'INPUT':<18}" + "".join(f"{np.nanmean(ivalid[k][allok]):>12.3f}" for k in VAL))
    for a in arms:
        g = allok & ok[a]
        print(f"   {a:<18}" + "".join(f"{np.nanmean(vcol(a,k)[g]):>12.3f}" for k in VAL))

    steps_cfg = int(o["config"].get("STEPS", 0))
    print("\n-- THE ITERATION BOUND: how many times it actually FIRED (the vacuity rule).")
    print("   An inertness certificate means nothing unless the bound was exercised; a capped arm")
    print("   is a DIFFERENT OPERATOR ('N iterations of restrained relaxation') and is labelled.")
    if steps_cfg == 0:
        print("   `steps = 0`: NO BOUND IS CONFIGURED.  This run is the deployed unbounded protocol,")
        print("   so there is nothing to fire and nothing to certify.  The bound was WITHDRAWN after")
        print("   the defect that motivated it was refuted -- see agentC_FINDINGS.md section 4.8.")
    tot = 0
    for a in ([] if steps_cfg == 0 else arms):
        h = np.array([bool(r["arms"][a].get("hit_cap", False)) for r in rows])
        tot += int(h.sum())
        if h.sum():
            print(f"   {a:<18} hit the bound on {int(h.sum()):>3} / {len(rows)} targets"
                  f"   ({' '.join([pdbs[i] for i in np.flatnonzero(h)][:8])}"
                  f"{' ...' if h.sum() > 8 else ''})")
    if steps_cfg != 0:
        print(f"   {'TOTAL':<18} {tot} arm-target calls hit the bound"
              f"{'  -- EXERCISED, so the inertness certificate is not vacuous' if tot else '  -- ZERO firings: the certificate is VACUOUS for this run'}")

    print("\n-- THE PARETO TEST (pre-registered): vs the incumbent k30, on BOTH axes.")
    print("   PAIRWISE gating: each arm against k30 on the targets where BOTH converge, which is")
    print("   the brief's 'a gated arm against its own gated input'.  The all-arm gate above costs")
    print(f"   {int((~allok).sum())} of {len(rows)} targets to ONE divergent arm and shifts the input mean from "
          f"{np.nanmean(inp):.3f} to {np.nanmean(inp[allok]):.3f} A, so it is reported but not used here.")
    ref = "k30"
    VAL2 = ("n_clash_2p6A", "bond_strain", "rama_favoured", "rama_outlier", "cis_frac")
    for a in arms:
        if a == ref:
            continue
        g = ok[a] & ok[ref]
        nm = [pdbs[i] for i in np.flatnonzero(g)]
        pr = PL.paired(col(a, "rmsd")[g], col(ref, "rmsd")[g], folds=folds[g], names=nm)
        print(f"   {a:<18} Ca {_fmt(pr)}")
        line = []
        for k in VAL2:
            q = PL.paired(vcol(a, k)[g], vcol(ref, k)[g], folds=folds[g], names=nm)
            ci = q.get("ci_fold", q["ci"])
            star = "*" if (ci[0] > 0 or ci[1] < 0) else " "
            line.append(f"{k[:9]} {q['mean']:+.4f}[{ci[0]:+.4f},{ci[1]:+.4f}]{star}")
        print(f"   {'':<18} " + "  ".join(line))
    print("   * = fold-clustered CI excludes zero.  For n_clash/bond_strain/rama_outlier/cis a")
    print("     POSITIVE value is WORSE than k30; for rama_favoured a positive value is BETTER.")
    if "frame_null" in o:
        fn = o["frame_null"]
        print(f"\n-- THE ROTATED-FRAME NULL (exactly zero by construction), reported with its "
              f"MAXIMUM, per BRIEF section 7:")
        print(f"   n = {fn['n']}  MAX |delta Ca-RMSD| = {fn['max_abs']:.6f} A   "
              f"mean {fn['mean_abs']:.6f} A   max |delta energy| {fn['max_abs_energy']:.4f} kcal/mol")
    return o


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "kv"
    {"kv": kv, "pareto": pareto, "reject": reject}[which]()
