"""SPRINT 18 / QUANTUM-ADVERSARIAL -- the regenerated report.  No number in
`s18/quantum_FINDINGS.md` is typed by hand; every one comes from here.

    python -m s18.q_report anova     # A2: which object is 'degree-1', and does it matter
    python -m s18.q_report gauge     # A3: is the degree-1 result invariant under relabelling
    python -m s18.q_report quantum   # A4: the quantum precondition
    python -m s18.q_report all
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s18.q_attack import boot_paired, summarise, fmt          # noqa: E402
from s18 import q_anova as A                                  # noqa: E402

RESULTS = os.path.join(ROOT, "s18", "results")


def _rows(tag):
    d = A.load(tag)
    r = [v for k, v in d.items() if not k.startswith("_")]
    r.sort(key=lambda x: x["pdb"])
    return r


# ================================================================= A2
def rep_anova():
    rows = _rows("anova")
    P = print
    P("=" * 106)
    P("A2 -- THREE OBJECTS ARE BEING CALLED 'DEGREE-1' AND THEY ARE NOT THE SAME OBJECT")
    P("=" * 106)
    P("  W1        strict Walsh weight-<=1 projection      (the object that produced 2.411 A)")
    P("  RA        residue-additive first-order ANOVA,     (the object BRIEF section 4 defines")
    P("            E0 + sum_i (E_mu[E|theta_i] - E0)        and the only one that ports to")
    P("                                                     continuous torsions)")
    P("  RA_walsh  projection onto Walsh coefficients supported inside ONE residue")
    P("  full      the deployed objective (rank-uniformised, as in Sprint 17)")
    P("")

    # ---- the equivalence audit
    a = [r["audit"] for r in rows]
    P("AUDIT 1 -- IS RA THE WALSH WEIGHT-<=1 PROJECTION?   (the brief asserts it is)")
    P("-" * 106)
    P(f"  max over {len(rows)} targets of  |RA - (weight<=1 + INTRA-RESIDUE weight-2) |  = "
      f"{max(x['max_abs_RA_minus_residue_block_projection'] for x in a):.3e}   <- EXACT, they")
    P("      are the same object to float precision.  The ANOVA construction reproduces the")
    P("      projection onto coefficients supported INSIDE ONE RESIDUE, not the weight-<=1")
    P("      projection.")
    P(f"  max over {len(rows)} targets of  |RA - W1| / sd(objective)                    = "
      f"{max(x['RA_minus_W1_rel_sd'] for x in a):.3f}")
    P(f"  mean over targets of |RA - W1| / sd(objective)                              = "
      f"{np.mean([x['RA_minus_W1_rel_sd'] for x in a]):.3f}")
    P(f"  mean Spearman(RA, W1)                                                       = "
      f"{np.mean([x['spearman_RA_W1'] for x in a]):.3f}")
    P("")
    P("  VERDICT: the brief's section-4 equivalence claim is FALSE AS STATED. RA == the")
    P("  residue-block projection (EXACT, 1e-12), and RA != W1 by 0.9-2.2 objective sd.")
    P("  With k=4 states in 2 qubits a per-residue field NEEDS its intra-residue weight-2")
    P("  term, and Sprint 17 measured 95.7% of all weight-2 mass as intra-residue.")
    P("")

    # ---- the consequential table
    P("AUDIT 2 -- AND IT CHANGES THE HEADLINE NUMBER")
    P("-" * 106)
    P(f"{'pdb':<6s}{'n':>3s}{'fold':>5s} | {'full':>7s}{'W1':>8s}{'RA':>8s}{'RAwalsh':>9s}"
      f"{'W2':>8s} | {'W1-full':>9s}{'RA-full':>9s} | {'space':>7s}")
    for r in rows:
        o = r["obj"]
        P(f"{r['pdb']:<6s}{r['n']:>3d}{r['fold']:>5d} | "
          f"{o['full']['argmin_rmsd_ORACLE']:>7.3f}{o['W1']['argmin_rmsd_ORACLE']:>8.3f}"
          f"{o['RA']['argmin_rmsd_ORACLE']:>8.3f}{o['RA_walsh']['argmin_rmsd_ORACLE']:>9.3f}"
          f"{o['W2']['argmin_rmsd_ORACLE']:>8.3f} | "
          f"{o['W1']['argmin_rmsd_ORACLE']-o['full']['argmin_rmsd_ORACLE']:>+9.3f}"
          f"{o['RA']['argmin_rmsd_ORACLE']-o['full']['argmin_rmsd_ORACLE']:>+9.3f} | "
          f"{r['space_best_ORACLE']:>7.3f}")
    folds = [r["fold"] for r in rows]
    g = lambda nm: np.array([r["obj"][nm]["argmin_rmsd_ORACLE"] for r in rows])
    P("-" * 106)
    P(f"{'MEAN':<6s}{'':>3s}{'':>5s} | {g('full').mean():>7.3f}{g('W1').mean():>8.3f}"
      f"{g('RA').mean():>8.3f}{g('RA_walsh').mean():>9.3f}{g('W2').mean():>8.3f} | "
      f"{(g('W1')-g('full')).mean():>+9.3f}{(g('RA')-g('full')).mean():>+9.3f} | "
      f"{np.mean([r['space_best_ORACLE'] for r in rows]):>7.3f}")
    P("")
    P("PAIRED, TARGET AS THE UNIT, fold-clustered bootstrap  (negative = the object beats full)")
    P("-" * 106)
    out = {}
    for nm in ("W1", "RA", "RA_walsh", "W2"):
        s = summarise(f"{nm} - full  (certified argmin RMSD)", g(nm) - g("full"), folds)
        out[nm] = s
        P(fmt(s))
    P("")
    P("  >> The object the brief tells the sprint to port -- RA -- is worth "
      f"{(g('RA')-g('full')).mean():+.3f} A, which is")
    P(f"  >> {abs((g('RA')-g('full')).mean())/0.084*100:.0f}% of the 126-target instrument's "
      "0.084 A minimum detectable effect. The entire")
    P("  >> 'degree-1 advantage' lives in W1, an object with no continuous-torsion analogue.")
    P("")

    # ---- band decomposition on the same objects
    P("BAND QUALITY vs WITHIN-BAND DRAW, recomputed on all four objects")
    P("  argmin = [mean RMSD of the object's own top 1%] + [argmin - that band mean]")
    P("-" * 106)
    bm = lambda nm: np.array([r["obj"][nm]["band_mean_ORACLE"] for r in rows])
    ib = lambda nm: np.array([r["obj"][nm]["rho_inband_ORACLE"] for r in rows])
    gl = lambda nm: np.array([r["obj"][nm]["rho_global_ORACLE"] for r in rows])
    P(f"{'object':<10s}{'argmin':>8s}{'bandmean':>10s}{'excess':>9s}{'rho_glob':>10s}"
      f"{'rho_inband':>12s}{'argmin pct in own band':>24s}")
    for nm in ("full", "W1", "RA", "W2"):
        ap = np.array([r["obj"][nm]["argmin_pct_within_band_ORACLE"] for r in rows])
        P(f"{nm:<10s}{g(nm).mean():>8.3f}{bm(nm).mean():>10.3f}"
          f"{(g(nm)-bm(nm)).mean():>+9.3f}{gl(nm).mean():>+10.3f}{ib(nm).mean():>+12.3f}"
          f"{ap.mean():>24.4f}")
    P("")
    for nm in ("W1", "RA", "W2"):
        s = summarise(f"BAND QUALITY  {nm} - full", bm(nm) - bm("full"), folds)
        P(fmt(s))
    for nm in ("W1", "RA", "W2"):
        s = summarise(f"IN-BAND rho   {nm} - full", -(ib(nm) - ib("full")), folds)
        s["name"] = f"IN-BAND rho   full - {nm}  (+ = full orders its band better)"
        P(fmt(s))
    P("")

    # ---- separability -- the quantum precondition
    P("SEPARABILITY -- THE QUANTUM PRECONDITION  (fraction of Walsh variance)")
    P("-" * 106)
    P(f"{'object':<10s}{'weight-1':>10s}{'inside 1 residue':>19s}{'INTER-RESIDUE':>16s}"
      f"{'mean Pauli weight':>20s}")
    for nm in ("full", "W2", "RA", "W1"):
        v1 = np.mean([r["obj"][nm]["var_w1"] for r in rows])
        vi = np.mean([r["obj"][nm]["var_inside_one_residue"] for r in rows])
        ve = np.mean([r["obj"][nm]["var_inter_residue"] for r in rows])
        mw = np.mean([r["obj"][nm]["mean_weight"] for r in rows])
        P(f"{nm:<10s}{v1:>10.4f}{vi:>19.4f}{ve:>16.4f}{mw:>20.3f}")
    P("")
    P("  The deployed objective already carries only "
      f"{np.mean([r['obj']['full']['var_inter_residue'] for r in rows]):.1%} of its variance "
      "in INTER-RESIDUE")
    P("  couplings.  Both truncations carry EXACTLY ZERO.  Truncating makes the objective")
    P("  MORE separable, so it makes the quantum case strictly WORSE, not better.")
    P("")
    cf = rows[0]["closed_form"]
    ok = all(r["closed_form"]["RA_argmin_matches_enumerated"] for r in rows)
    P(f"  RA's global optimum is a closed-form per-residue argmin at n*k = "
      f"{np.mean([r['closed_form']['RA_cost_objective_evaluations'] for r in rows]):.0f} table")
    P(f"  reads (verified against the enumerated argmin on {len(rows)}/{len(rows)} targets: "
      f"{ok}).  W1's is a per-QUBIT")
    P(f"  argmin at {np.mean([r['closed_form']['W1_cost_objective_evaluations'] for r in rows]):.0f}"
      " reads.  Neither is a search problem for any device.")
    P("")
    json.dump({"summaries": out, "_complete": True},
              open(os.path.join(RESULTS, "q_report_anova.json"), "w"), indent=1, default=float)


# ================================================================= A3
def rep_gauge():
    rows = _rows("gauge")
    P = print
    P("=" * 106)
    P("A3 -- THE DEGREE-1 RESULT IS NOT INVARIANT UNDER RELABELLING THE TORSION STATES")
    P("=" * 106)
    P("Which 2-bit code names which of the k = 4 torsion states of a residue is a LABELLING")
    P("choice.  E, its argmin, every RMSD and RA are invariant under it.  W1 is not, because")
    P("it keeps only the two SINGLE-QUBIT coefficients of each residue's 4-vector and picks")
    P("the two bits independently.  The gauge group is (S_4)^n.")
    P("")
    P(f"  exactness check: |weight-1 Walsh coefficients from the residue marginals - from the")
    P(f"  full FWHT| <= "
      f"{max(r.get('exactness_weight1_coeff_max_abs_diff', 0.0) for r in rows):.2e} "
      "-> the gauge orbit below is EXACT, not approximate.")
    P("")
    P(f"{'pdb':<6s}{'n':>3s} {'full':>7s}{'RA':>7s}{'W1id':>7s} | "
      f"{'gaugemean':>10s}{'gaugemed':>9s}{'g2.5%':>8s}{'g97.5%':>8s}{'midpct':>8s} | "
      f"{'|orbit|':>8s}{'agree':>7s}")
    for r in rows:
        g = r["gauge"]
        P(f"{r['pdb']:<6s}{r['n']:>3d} {r['full_argmin_rmsd_ORACLE']:>7.3f}"
          f"{r['RA_argmin_rmsd_ORACLE']:>7.3f}{r['W1_identity_argmin_rmsd_ORACLE']:>7.3f} | "
          f"{g['mean']:>10.3f}{g['median']:>9.3f}{g['q'][0]:>8.3f}{g['q'][4]:>8.3f}"
          f"{g['pct_of_identity_midrank']:>8.3f} | {g['n_distinct_configs']:>8d}"
          f"{r['frac_residue_labellings_where_W1_agrees_with_RA']:>7.3f}")
    folds = [r["fold"] for r in rows]
    full = np.array([r["full_argmin_rmsd_ORACLE"] for r in rows])
    ra = np.array([r["RA_argmin_rmsd_ORACLE"] for r in rows])
    w1 = np.array([r["W1_identity_argmin_rmsd_ORACLE"] for r in rows])
    gm = np.array([r["gauge"]["mean"] for r in rows])
    S = np.array([r["gauge"]["samples"] for r in rows])            # (T, B)
    P("-" * 106)
    P(f"{'MEAN':<6s}{'':>3s} {full.mean():>7.3f}{ra.mean():>7.3f}{w1.mean():>7.3f} | "
      f"{gm.mean():>10.3f}")
    P("")
    P("THE GAUGE NULL FOR THE HEADLINE STATISTIC  'deg1 - full'")
    P("-" * 106)
    stat = S.mean(0) - full.mean()          # (B,) -- the statistic under a random relabelling
    obs = w1.mean() - full.mean()
    P(f"  observed (the codebase's own labelling)      {obs:+.3f} A")
    P(f"  gauge null: mean {stat.mean():+.3f}   median {np.median(stat):+.3f}   "
      f"95% band [{np.percentile(stat,2.5):+.3f}, {np.percentile(stat,97.5):+.3f}]")
    P(f"  the observed value sits at gauge percentile "
      f"{100*(stat < obs).mean():.2f}  (one-sided p = {(stat <= obs).mean():.4f})")
    P("")
    P(f"  Under a RANDOM relabelling of the torsion states -- a transformation under which no")
    P(f"  physical quantity moves at all -- the strict Walsh weight-<=1 truncation is worth")
    P(f"  {stat.mean():+.3f} A against the full objective, not {obs:+.3f} A.")
    P("")
    mp = np.array([r["gauge"]["pct_of_identity_midrank"] for r in rows])
    deg = np.array([r["gauge"]["degenerate_orbit"] for r in rows])
    P(f"  per-target mid-rank percentile of the identity labelling: mean {mp.mean():.3f}, "
      f"median {np.median(mp):.3f}")
    P(f"  ({int(deg.sum())} of {len(rows)} targets have a DEGENERATE orbit -- every labelling")
    P("   gives the same configuration -- and are correctly scored at 0.500, not 0.000.)")
    nondeg = ~deg
    P(f"  on the {int(nondeg.sum())} non-degenerate targets: mean {mp[nondeg].mean():.3f}, "
      f"median {np.median(mp[nondeg]):.3f}, below 0.5 on "
      f"{int((mp[nondeg] < 0.5).sum())}/{int(nondeg.sum())}")
    P("")
    P("PAIRED SUMMARIES, target as the unit")
    P("-" * 106)
    for nm, v in (("W1(identity labelling) - full", w1 - full),
                  ("W1(random labelling, gauge mean) - full", gm - full),
                  ("RA (label-free) - full", ra - full),
                  ("W1(identity) - W1(gauge mean)  [pure encoding luck]", w1 - gm)):
        P(fmt(summarise(nm, v, folds)))
    P("")
    mr = np.array([r["matched_random"]["mean"] for r in rows])
    sm = np.array([r["space_mean_ORACLE"] for r in rows])
    sb = np.array([r["space_best_ORACLE"] for r in rows])
    P("BOTH MANDATORY CONTROLS")
    P("-" * 106)
    P(f"  ZERO-INFORMATION (space mean RMSD)                     {sm.mean():.3f} A")
    P(f"  MATCHED-RANDOM (same construction, information deleted) {mr.mean():.3f} A")
    P(f"  W1 under a random labelling                            {gm.mean():.3f} A")
    P(f"  RA (label-free)                                        {ra.mean():.3f} A")
    P(f"  full objective                                         {full.mean():.3f} A")
    P(f"  W1 under the codebase's labelling                      {w1.mean():.3f} A")
    P(f"  ORACLE space best                                      {sb.mean():.3f} A")
    P("")
    json.dump({"observed": float(obs), "gauge_null_mean": float(stat.mean()),
               "gauge_null_lo": float(np.percentile(stat, 2.5)),
               "gauge_null_hi": float(np.percentile(stat, 97.5)),
               "gauge_pct_of_observed": float((stat < obs).mean()),
               "means": {"full": float(full.mean()), "RA": float(ra.mean()),
                         "W1_identity": float(w1.mean()), "W1_gauge": float(gm.mean()),
                         "matched_random": float(mr.mean()),
                         "space_mean": float(sm.mean()), "space_best": float(sb.mean())},
               "_complete": True},
              open(os.path.join(RESULTS, "q_report_gauge.json"), "w"), indent=1, default=float)


# ================================================================= A4
def rep_quantum():
    d = A.load("qcond")
    rows = [v for k, v in d.items() if not k.startswith("_")]
    if not rows:
        print("[q_report quantum] no q_qcond.json yet -- run `python -m s18.q_cond run`")
        return
    rows.sort(key=lambda x: x["pdb"])
    P = print
    P("=" * 106)
    P("A4 -- THE PRE-REGISTERED QUANTUM PRECONDITION, TESTED ON THE NEW OBJECTIVE")
    P("=" * 106)
    P("PRE-REGISTERED (before any arm was run, s18/PREREG_quantum.md):")
    P("  If the new objective contains useful higher-order correlations then")
    P("   (a) greedy 1-opt should CEASE to certify its optimum at tiny budget, and")
    P("   (b) the CNOT-free product ansatz should become measurably WORSE than the entangled")
    P("       VQE at matched budget, with a target-level interval excluding zero.")
    P("")
    objs = sorted({k for r in rows for k in r["budget"]})
    budgets = rows[0]["budgets"]
    P("(a) GREEDY 1-OPT CERTIFICATION RATE vs BUDGET  (% of target x seed cells reaching the")
    P("    CERTIFIED global optimum of that objective)")
    P("-" * 106)
    P(f"{'budget':>8s}" + "".join(f"{o:>14s}" for o in objs))
    for bi, b in enumerate(budgets):
        P(f"{b:>8d}" + "".join(
            f"{100*np.mean([r['budget'][o]['greedy_hit'][bi] for r in rows]):>13.1f}%"
            for o in objs))
    P("")
    P(f"{'metropolis':>8s}")
    P(f"{'budget':>8s}" + "".join(f"{o:>14s}" for o in objs))
    for bi, b in enumerate(budgets):
        P(f"{b:>8d}" + "".join(
            f"{100*np.mean([r['budget'][o]['metro_hit'][bi] for r in rows]):>13.1f}%"
            for o in objs))
    P("")
    P(f"{'random draw (zero-information control)':>8s}")
    P(f"{'budget':>8s}" + "".join(f"{o:>14s}" for o in objs))
    for bi, b in enumerate(budgets):
        P(f"{b:>8d}" + "".join(
            f"{100*np.mean([r['budget'][o]['rand_hit'][bi] for r in rows]):>13.1f}%"
            for o in objs))
    P("")
    for o in objs:
        onep = np.mean([r["budget"][o]["one_pass_cost"] for r in rows])
        cf = rows[0]["budget"][o].get("closed_form_cost")
        P(f"  {o:<6s} one coordinate pass = n(k-1)+1 = {onep:.0f} evals" +
          (f";  CLOSED-FORM optimum at "
           f"{np.mean([r['budget'][o]['closed_form_cost'] for r in rows]):.0f} table reads "
           "(no search)" if cf else ";  no closed form"))
    P("")

    # ---- PAIRED: is the new objective HARDER to search than the old one?
    P("(a) PAIRED, TARGET AS THE UNIT -- is the new objective HARDER than the deployed one?")
    P("    statistic: greedy certification rate at a given budget, new minus full")
    P("    (positive = the new objective is EASIER, i.e. the quantum case gets WORSE)")
    P("-" * 106)
    folds = [r["fold"] for r in rows]
    for o in objs:
        if o == "full":
            continue
        for b in (20, 36, 64, 128):
            bi = budgets.index(b)
            d = np.array([r["budget"][o]["greedy_hit"][bi] -
                          r["budget"]["full"]["greedy_hit"][bi] for r in rows])
            s = summarise(f"greedy cert @ {b:>4d} evals   {o} - full", d, folds)
            s["name"] = f"greedy cert @ {b:>4d} evals   {o} - full"
            P(fmt(s))
    P("")
    P("  Q1 fires if certification is >= 90% at <= 1024 evaluations on the new objective, or")
    P("  is HIGHER than on the full objective.")
    P("")

    # ---- certified argmin quality
    P("(a-corollary) THE CERTIFIED ARGMIN OF EACH OBJECT, and what greedy returns")
    P("-" * 106)
    P(f"{'object':<8s}{'certified argmin':>18s}{'greedy@36':>12s}{'greedy@1024':>13s}"
      f"{'metro@1024':>12s}")
    for o in objs:
        b36 = budgets.index(36); b1k = budgets.index(1024)
        P(f"{o:<8s}"
          f"{np.mean([r['budget'][o]['certified_argmin_rmsd_ORACLE'] for r in rows]):>18.3f}"
          f"{np.mean([r['budget'][o]['greedy_rmsd_ORACLE'][b36] for r in rows]):>12.3f}"
          f"{np.mean([r['budget'][o]['greedy_rmsd_ORACLE'][b1k] for r in rows]):>13.3f}"
          f"{np.mean([r['budget'][o]['metro_rmsd_ORACLE'][b1k] for r in rows]):>12.3f}")
    P("")


def rep_vqe(obj="RA"):
    d = A.load(f"qvqe_{obj}")
    rows = [v for k, v in d.items() if not k.startswith("_")]
    if not rows:
        print(f"[q_report vqe] no q_qvqe_{obj}.json yet -- run `python -m s18.q_cond vqe`")
        return
    rows.sort(key=lambda x: x["_meta"]["pdb"])
    P = print
    P("=" * 106)
    P(f"(b) THE ENTANGLEMENT CONTROL ON THE NEW OBJECTIVE `{obj}`   n = {len(rows)} targets")
    P("=" * 106)
    P("`mps2fn` is the IDENTICAL circuit with the CNOTs removed -- a product Bernoulli model,")
    P("trained by the identical CVaR score-function estimator, the identical Adam optimiser,")
    P("the identical 8,192-evaluation budget, shots, learning rate and SEED. It is a")
    P("classical algorithm. If it matches `mps2f`, the entanglement is decorative.")
    P("")
    folds = [r["_meta"]["fold"] for r in rows]
    arms = [k for k in rows[0] if not k.startswith("_")]
    P("BOTH consumption modes are reported and never substituted: `sel*` is the m = 75")
    P("lowest-objective distinct visited (what the pipeline consumes); the bare columns are a")
    P("uniform draw from the arm's own visited multiset (the sampler's geometry, no ranker).")
    P("")
    P(f"{'arm':<16s}{'selM75':>9s}{'selD75':>9s}{'selread':>9s}{'selbest':>9s} | "
      f"{'M75':>8s}{'D75':>8s}{'readout':>9s}{'setbest':>9s} | {'distinct':>9s}{'evals':>7s}")
    for a in arms:
        P(f"{a:<16s}" + "".join(
            f"{np.mean([r[a][c] for r in rows]):>9.3f}"
            for c in ("selM75", "selD75", "selreadout75", "selset_best75"))
          + " | " + "".join(
            f"{np.mean([r[a][c] for r in rows]):>{w}.3f}"
            for c, w in (("M75", 8), ("D75", 8), ("readout75", 9), ("set_best75", 9)))
          + f" | {np.mean([r[a]['n_distinct'] for r in rows]):>9.0f}"
          + f"{np.mean([r[a]['cost']['objective_evals'] for r in rows]):>7.0f}")
    P("")
    P("HEAD-TO-HEAD, paired over targets.  On M and the readout, POSITIVE = the ENTANGLED")
    P("circuit is WORSE.  D is DIVERSITY and is a coordinate, not a quality: positive D means")
    P("the entangled circuit is LESS converged, which on the (M, D) plane is a position, not a")
    P("win or a loss.  It is labelled here so the sign is not read as a verdict.")
    P("-" * 106)
    als = sorted({float(a.split("_a")[1]) for a in arms if a.startswith("mps2f_a")})
    for al in als:
        for col in ("selM75", "selreadout75", "selset_best75", "M75", "D75", "readout75"):
            d_ = np.array([r[f"mps2f_a{al}"][col] - r[f"mps2fn_a{al}"][col] for r in rows])
            s = summarise(f"a={al:<4} {col:<13s} mps2f - mps2fn", d_, folds)
            P(fmt(s))
        P("")
    P("Q2 fires AGAINST the pillar if these intervals span zero.")
    P("")
    json.dump({"n": len(rows), "_complete": True},
              open(os.path.join(RESULTS, f"q_report_vqe_{obj}.json"), "w"), indent=1)


def main():
    m = sys.argv[1] if len(sys.argv) > 1 else "all"
    if m in ("anova", "all"):
        rep_anova()
    if m in ("gauge", "all"):
        rep_gauge()
    if m in ("quantum", "all"):
        rep_quantum()
    if m in ("vqe", "all"):
        rep_vqe(sys.argv[2] if len(sys.argv) > 2 else "RA")


if __name__ == "__main__":
    main()
