"""SPRINT 17 / QUANTUM -- every table in `s17/quantum_FINDINGS.md`, regenerated from the
checkpoints so no number in the findings is typed by hand.

    python -m s17.q_report pareto      # P1, the Pareto test, and P2's cost table
    python -m s17.q_report novelty     # P5, mixtures
    python -m s17.q_report repr        # P4, ensemble reweighting + alpha-as-temperature
    python -m s17.q_report theory      # P6, the Walsh spectrum and the representability proof
    python -m s17.q_report local       # P3, the local VQE
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from s12 import instrument as I
from s17 import q_lib as L
from s17 import q_pareto as P


def _cells(tag):
    d = L.ck_load(tag)
    return {k: v for k, v in d.items() if not k.startswith("_")}


def _by_target(cells):
    """Seeds averaged WITHIN a target before anything else. The unit is the target."""
    out = {}
    for k, v in cells.items():
        pdb = k.rsplit("_", 1)[0]
        out.setdefault(pdb, []).append(v)
    return out


def _agg_point(runs, arm, m=75, pre=""):
    """The arm's (M, D) point for one target: quadratic mean over seeds, because M and D
    enter the identity quadratically and an arithmetic mean of RMSDs is the recorded trap."""
    Ms = [r[arm][f"{pre}M{m}"] for r in runs if arm in r and f"{pre}M{m}" in r[arm]]
    Ds = [r[arm][f"{pre}D{m}"] for r in runs if arm in r and f"{pre}D{m}" in r[arm]]
    Rs = [r[arm][f"{pre}readout{m}"] for r in runs if arm in r and f"{pre}M{m}" in r[arm]]
    Bs = [r[arm]["set_best_rmsd"] for r in runs if arm in r]
    Ns = [r[arm]["n_distinct"] for r in runs if arm in r]
    if not Ms:
        return None
    return {"M": float(np.sqrt(np.mean(np.square(Ms)))),
            "D": float(np.sqrt(np.mean(np.square(Ds)))),
            "readout": float(np.sqrt(np.mean(np.square(Rs)))),
            "set_best": float(np.mean(Bs)),
            "n_distinct": float(np.mean(Ns))}


def pareto(m=75, pre=""):
    cells = _cells("pareto")
    bt = _by_target(cells)
    pdbs = sorted(bt)
    folds = np.array([L.inst(p).fold for p in pdbs])
    arms = sorted({a for r in cells.values() for a in r if not a.startswith("_")})
    Q = [a for a in arms if a.startswith("vqe")]
    C = [a for a in arms if a in P.CLASSICAL_MATCHED]
    X = [a for a in arms if a in P.CIRCUIT_DERIVED]

    pts = {p: {a: _agg_point(bt[p], a, m, pre) for a in arms} for p in pdbs}

    mode = ("UNIFORM draw from the arm's visited multiset (no ranker)" if pre == ""
            else "TOP-m BY THE OBJECTIVE among distinct visited (the pipeline's own cut)")
    print(f"\n=== P1  THE PARETO TEST, set size m = {m}, consumption = {mode}, "
          f"n = {len(pdbs)} targets, seeds averaged within target ===")
    print("     lower M is better, higher D is better.  tau = "
          f"{L.TAU} A.  floor = {L.FLOOR} A.\n")

    # ---------------- the two curves, means over targets --------------------
    print(f"{'arm':22s} {'M':>7s} {'D':>7s} {'readout':>8s} {'set_best':>9s} "
          f"{'distinct':>9s} {'obj_evals':>10s} {'wall_s':>7s}")
    order = Q + ["|"] + C + ["|"] + X
    for a in order:
        if a == "|":
            print("-" * 84)
            continue
        v = [pts[p][a] for p in pdbs if pts[p].get(a)]
        if not v:
            continue
        cost = np.mean([r[a]["cost"]["objective_evals"] for r in
                        sum(bt.values(), []) if a in r])
        wall = np.mean([r[a]["cost"]["wall_s"] for r in sum(bt.values(), []) if a in r])
        print(f"{a:22s} {np.mean([x['M'] for x in v]):7.3f} "
              f"{np.mean([x['D'] for x in v]):7.3f} "
              f"{np.mean([x['readout'] for x in v]):8.3f} "
              f"{np.mean([x['set_best'] for x in v]):9.3f} "
              f"{np.mean([x['n_distinct'] for x in v]):9.0f} "
              f"{cost:10.0f} {wall:7.1f}")

    # ---------------- THE PRIMARY STATISTIC: epsilon-dominance ---------------
    print(f"\n--- P1 PRIMARY: additive epsilon-dominance of each QUANTUM arm against the "
          f"{len(C)}-point classical attainable set")
    print("    eps = min_c max( M_c - M_vqe , D_vqe - D_c ),  in Angstroms.")
    print("    eps > 0  =>  NO classical point dominates the VQE point  =>  a QUANTUM WIN.")
    print("    eps <= 0 =>  some classical point weakly dominates it.")
    print("    Defined on EVERY target, including those where the VQE out-diversifies every")
    print("    classical arm -- the case the conditional M-margin below leaves undefined.\n")
    print(f"{'arm':22s} {'mean eps':>9s} {'median':>8s} {'CI95':>20s} "
          f"{'eps>0':>8s} {'foldsign':>9s}")
    rows = {}
    for a in Q:
        eps = np.asarray([L.eps_dominance(pts[p][a],
                                          [pts[p][c] for c in C if pts[p].get(c)])
                          for p in pdbs], float)
        v = L.verdict(eps, folds=folds, label=a)
        rows[a] = {"eps": eps.tolist(), **v}
        print(f"{a:22s} {v['mean']:+9.3f} {v['median']:+8.3f} "
              f"[{v['ci95'][0]:+7.3f},{v['ci95'][1]:+7.3f}] "
              f"{int((eps > 0).sum()):3d}/{len(pdbs):<4d} "
              f"{v['folds_same_sign']}/{v['n_folds']}")
    print("\n    THE NULL THIS MUST BE READ AGAINST: the same statistic for each CLASSICAL "
          "arm, leave-one-out against the rest of the classical set.")
    for a in C:
        eps = np.asarray([L.eps_dominance(pts[p][a],
                                          [pts[p][c] for c in C
                                           if c != a and pts[p].get(c)])
                          for p in pdbs], float)
        if not np.isfinite(eps).all():
            continue
        print(f"{a:22s} {eps.mean():+9.3f} {np.median(eps):+8.3f} "
              f"{'':>20s} {int((eps > 0).sum()):3d}/{len(pdbs):<4d}")

    # ---------------- the conditional M-margin, secondary --------------------
    print(f"\n--- SECONDARY: M-margin = M_vqe - min{{ M_c : D_c >= D_vqe - tau }}. "
          "NEGATIVE = quantum win.")
    print("    UNDEFINED where the VQE out-diversifies every classical arm; that count is "
          "printed and those cells are NOT dropped from the primary above.\n")
    print(f"{'arm':22s} {'mean':>8s} {'median':>8s} {'CI95':>20s} {'W/L':>7s} "
          f"{'undef':>8s}")
    for a in Q:
        marg = np.asarray([L.dominance_margin(pts[p][a],
                                              [pts[p][c] for c in C if pts[p].get(c)])
                           for p in pdbs], float)
        rows[a]["M_margin"] = marg.tolist()
        nd = int((~np.isfinite(marg)).sum())
        mm = marg[np.isfinite(marg)]
        if mm.size < 3:
            print(f"{a:22s} {'--':>8s} {'--':>8s} {'--':>20s} {'--':>7s} "
                  f"{nd:3d}/{len(pdbs):<4d}")
            continue
        v = L.verdict(mm, label=a)
        print(f"{a:22s} {v['mean']:+8.3f} {v['median']:+8.3f} "
              f"[{v['ci95'][0]:+7.3f},{v['ci95'][1]:+7.3f}] {v['W']:3d}/{v['L']:<3d} "
              f"{nd:3d}/{len(pdbs):<4d}")

    # the falsifier names a NARROWER classical set -- "the classical thermostat plus 1-opt".
    # Report it separately so the pre-registered rule is scored on its own terms.
    Cn = [a for a in C if a.startswith("metro") or a == "greedy"]
    print(f"\n    THE FALSIFIER'S LITERAL SET -- thermostat ladder + 1-opt only "
          f"({len(Cn)} points), epsilon-dominance; eps > 0 = quantum win:")
    for a in Q:
        marg = np.asarray([L.eps_dominance(pts[p][a],
                                           [pts[p][c] for c in Cn if pts[p].get(c)])
                           for p in pdbs], float)
        v = L.verdict(marg, folds=folds, label=a)
        rows.setdefault(a, {})["eps_narrow"] = marg.tolist()
        rows[a]["verdict_narrow"] = {k: v[k] for k in ("mean", "median", "ci95", "W", "L")}
        print(f"{a:22s} {v['mean']:+8.3f} {v['median']:+8.3f} "
              f"[{v['ci95'][0]:+7.3f},{v['ci95'][1]:+7.3f}] "
              f"{int((marg > 0).sum()):3d}/{len(pdbs):<4d}")

    # the same for the circuit-derived arms, as a reference
    print("\n    (reference, not classical: circuit-derived arms through the same rule)")
    for a in X:
        marg = np.asarray([L.eps_dominance(pts[p][a],
                                           [pts[p][c] for c in C if pts[p].get(c)])
                           for p in pdbs], float)
        v = L.verdict(marg, folds=folds, label=a)
        print(f"{a:22s} {v['mean']:+8.3f} {v['median']:+8.3f} "
              f"[{v['ci95'][0]:+7.3f},{v['ci95'][1]:+7.3f}] "
              f"{int((marg > 0).sum()):3d}/{len(pdbs):<4d}")

    # ---------------- joint Pareto front membership -------------------------
    print("\n--- per-target JOINT Pareto front (quantum + classical together): "
          "how often is each arm non-dominated?")
    memb = {a: 0 for a in Q + C}
    for p in pdbs:
        names = [a for a in Q + C if pts[p].get(a)]
        pp = [pts[p][a] for a in names]
        for i in L.pareto_front(pp):
            memb[names[i]] += 1
    for a in sorted(memb, key=lambda z: -memb[z]):
        if memb[a]:
            fam = "QUANTUM " if a.startswith("vqe") else "classical"
            print(f"   {fam} {a:22s} {memb[a]:3d}/{len(pdbs)}")
    print("   quantum arms on the joint front, total cells: "
          f"{sum(memb[a] for a in Q)} ; classical: {sum(memb[a] for a in C)}")

    # ---------------- the reverse diversity preference ----------------------
    print("\n--- the OPPOSITE preference (a consumer that wants LOW diversity, e.g. an "
          "argmin readout): margin against classical points with D <= D_vqe + tau")
    for a in Q:
        marg = []
        for p in pdbs:
            cl = [pts[p][c] for c in C if pts[p].get(c)
                  and pts[p][c]["D"] <= pts[p][a]["D"] + L.TAU]
            marg.append(pts[p][a]["M"] - min(c["M"] for c in cl) if cl else np.nan)
        marg = np.asarray(marg, float)
        if np.isfinite(marg).all():
            v = L.verdict(marg, folds=folds, label=a)
            print(f"{a:22s} {v['mean']:+8.3f} [{v['ci95'][0]:+7.3f},{v['ci95'][1]:+7.3f}] "
                  f"{v['W']:3d}/{v['L']:<3d}")

    # ---------------- the identity, aggregated ------------------------------
    resid = []
    for p in pdbs:
        for a in arms:
            q = pts[p].get(a)
            if q:
                resid.append(q["readout"] ** 2 - (q["M"] ** 2 - q["D"] ** 2))
    print(f"\n--- identity check across every aggregated point: max |readout^2 - (M^2 - D^2)| "
          f"= {np.abs(resid).max():.3e} A^2  (n = {len(resid)})")
    out = {"m": m, "n_targets": len(pdbs), "targets": pdbs,
           "points": pts, "margins": rows,
           "front_membership": memb}
    with open(os.path.join(L.RESULTS,
                           f"quantum_pareto_report_{pre or 'unif'}m{m}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


def novelty():
    cells = _cells("pareto")
    bt = _by_target(cells)
    pdbs = sorted(bt)
    print(f"\n=== P5  MIXTURES: does the VQE contribute UNIQUE high-quality structures? "
          f"n = {len(pdbs)} targets ===")
    print("   'elite' = the ORACLE true top 0.1% of the register. Support sets are the "
          "arms' visited configurations.\n")
    keys = ("n_elite", "elite_found_by_quantum_union", "elite_found_by_classical_union",
            "elite_only_quantum", "elite_only_classical", "n_quantum_union",
            "n_classical_union", "frac_quantum_support_novel")
    agg = {k: [] for k in keys}
    for p in pdbs:
        for r in bt[p]:
            for k in keys:
                agg[k].append(r["_novelty"][k])
    for k in keys:
        print(f"   {k:36s} mean {np.mean(agg[k]):10.3f}   median "
              f"{np.median(agg[k]):10.3f}   total {np.sum(agg[k]):10.0f}")
    print("\n   LEAVE-ONE-OUT novelty (symmetric across arms): elite configurations found "
          "by an arm and by NO other arm")
    arms = sorted({a for r in cells.values() for a in r if not a.startswith("_")})
    loo = {a: np.mean([r[a]["elite_unique_loo"] for r in cells.values() if a in r])
           for a in arms}
    hit = {a: np.mean([r[a]["elite_hits"] for r in cells.values() if a in r])
           for a in arms}
    for a in sorted(arms, key=lambda z: -hit[z]):
        fam = "QUANTUM " if a.startswith("vqe") else ("circuit  " if a in P.CIRCUIT_DERIVED
                                                      else "classical")
        print(f"   {fam} {a:22s} elite hits {hit[a]:8.2f}   unique(LOO) {loo[a]:6.3f}")


def theory():
    d = L.ck_load("theory")
    rep = d.get("representability")
    if rep:
        print("\n=== P6.4  REPRESENTABILITY (EXACT) ===")
        print(f"   mps2f reaches a computational basis state exactly: min max-prob = "
              f"{rep['min_max_prob']:.12f}; distinct basis states reached in 64 tries = "
              f"{rep['distinct_basis_states_reached_in_64_tries']}/64; PASS = "
              f"{rep['PASS_delta_representable']}")
        print("   => the unbudgeted variational optimum of CVaR_alpha at ANY alpha is a "
              "DELTA on the objective's certified argmin (D = 0).")
    specs = {k[5:]: v for k, v in d.items() if k.startswith("spec_")}
    if not specs:
        return
    print(f"\n=== P6.1  WALSH (PAULI-Z) SPECTRUM, n = {len(specs)} targets ===")
    for field in ("hamil_uniformised", "hamil_raw", "prior_uniformised",
                  "legacy_uniformised", "TRUTH_uniformised_ORACLE"):
        rows = [v[field] for v in specs.values() if field in v]
        if not rows:
            continue
        mw = np.mean([r["mean_weight"] for r in rows])
        cum = np.mean([r["cum_var_frac"][:9] for r in rows], axis=0)
        print(f"   {field:28s} mean Pauli weight {mw:5.2f}   cumulative variance fraction "
              "by weight <= 1..8: " + " ".join(f"{x:.3f}" for x in cum[1:9]))
    print("\n=== P6.2  DEGREE TRUNCATION of the uniformised deployed objective ===")
    print("   (what a weight-<=d approximation retains; ORACLE columns are post-hoc)")
    tr = [v["hamil_uniformised"]["truncation"] for v in specs.values()
          if "truncation" in v.get("hamil_uniformised", {})]
    if tr:
        degs = sorted(tr[0], key=int)
        print(f"   {'d':>3s} {'rho w/ full':>12s} {'rho in-band':>12s} "
              f"{'rho w/ RMSD':>12s} {'argmin RMSD':>12s} {'argmin pct':>11s}")
        for dd in degs:
            print(f"   {dd:>3s} {np.mean([t[dd]['rho_with_full'] for t in tr]):12.3f} "
                  f"{np.mean([t[dd]['rho_inband_with_full'] for t in tr]):12.3f} "
                  f"{np.mean([t[dd]['rho_with_rmsd_ORACLE'] for t in tr]):12.3f} "
                  f"{np.mean([t[dd]['argmin_rmsd_ORACLE'] for t in tr]):12.3f} "
                  f"{np.mean([t[dd]['argmin_pct_ORACLE'] for t in tr]):11.3f}")
        print(f"   {'full':>3s} {1.0:12.3f} {1.0:12.3f} "
              f"{np.mean([v.get('rho_full', np.nan) for v in specs.values()]):12s}"
              if False else
              f"   full  certified argmin RMSD "
              f"{np.mean([v['certified_argmin_rmsd_ORACLE'] for v in specs.values()]):.3f}"
              f"   space best "
              f"{np.mean([v['rmsd_best_in_space_ORACLE'] for v in specs.values()]):.3f}")
    rng_ = {k[6:]: v for k, v in d.items() if k.startswith("range_")}
    if rng_:
        print(f"\n=== P6.1b  WHERE THE WEIGHT-2 WALSH MASS SITS (qubit separation), "
              f"n = {len(rng_)} targets ===")
        print("   2 qubits encode one residue, so separation 1 is INTRA-RESIDUE.")
        for f in ("hamil_uniformised", "TRUTH_uniformised_ORACLE"):
            ms = np.mean([v[f]["mean_separation"] for v in rng_.values()])
            w2 = np.mean([v[f]["frac_within_2_qubits"] for v in rng_.values()])
            w1 = np.mean([v[f]["frac_within_1_residue"] for v in rng_.values()])
            print(f"   {f:28s} mean separation {ms:.2f} qubits   "
                  f"frac at separation <= 1 (intra-residue) {w1:.3f}   <= 2 {w2:.3f}")
    bud = {k[7:]: v for k, v in d.items() if k.startswith("budget_")}
    if bud:
        print(f"\n=== P6.3  HOW CHEAP IS THIS SEARCH PROBLEM? n = {len(bud)} targets, "
              "8 seeds, fraction reaching the CERTIFIED GLOBAL OPTIMUM ===")
        bs = [b for b in sorted(bud[list(bud)[0]], key=lambda z: (z.isdigit(), int(z)
                                                                 if z.isdigit() else 0))
              if b.isdigit()]
        print(f"   {'budget':>8s} {'greedy %cert':>13s} {'greedy RMSD':>12s} "
              f"{'metro %cert':>12s} {'targets at 100%':>16s}")
        for b in bs:
            gc = np.mean([v[b]["greedy_frac_certified"] for v in bud.values()])
            gr = np.mean([v[b]["greedy_rmsd_ORACLE"] for v in bud.values()])
            mc = np.mean([v[b]["metro_frac_certified"] for v in bud.values()])
            n100 = sum(1 for v in bud.values() if v[b]["greedy_frac_certified"] >= 1.0)
            print(f"   {b:>8s} {100*gc:13.1f} {gr:12.3f} {100*mc:12.1f} "
                  f"{n100:6d}/{len(bud):<9d}")
        print(f"   one-pass coordinate-descent cost n(k-1)+1 = "
              f"{np.mean([v['one_pass_cost_n_times_k'] for v in bud.values()]):.0f} "
              f"evaluations; the CVaR-VQE budget used throughout Sprints 15-17 is 8,192.")
    g = d.get("gradient")
    if g:
        print("\n=== P6  GRADIENT ESTIMATOR SIGNAL-TO-NOISE ===")
        for nq in sorted({r["n_qubits"] for r in g}):
            rr = [r for r in g if r["n_qubits"] == nq]
            print(f"   n_qubits {nq:3d}  |g_exact| {np.mean([r['norm_exact'] for r in rr]):.4f}"
                  f"  |g_sampled| {np.mean([r['norm_sampled'] for r in rr]):.4f}"
                  f"  cos(sampled, exact) {np.mean([r['cos_sampled_vs_exact'] for r in rr]):+.4f}"
                  f"  n_params {rr[0]['n_params']}")


def repr_report():
    cells = _cells("repr")
    if not cells:
        print("no repr cells yet")
        return
    pdbs = sorted({k.rsplit("_", 1)[0] for k in cells})
    folds = np.array([L.inst(p).fold for p in pdbs])
    alphas = sorted({k for v in cells.values() for k in v if k.startswith("a")},
                    key=lambda z: float(z[1:]))

    def per_target(alpha, field):
        out = []
        for p in pdbs:
            vals = [v[alpha][field] for k, v in cells.items()
                    if k.rsplit("_", 1)[0] == p and alpha in v]
            out.append(np.mean(vals) if vals else np.nan)
        return np.asarray(out, float)

    print(f"\n=== ALPHA AS A TEMPERATURE, EXACT (n = {len(pdbs)} targets) ===")
    print("   KL(q_theta || family) in BITS at the best-fit temperature. Small = q_theta IS "
          "that law.")
    print(f"   {'alpha':>6s} {'H(q) bits':>10s} {'T_eff/Boltz':>12s} {'KL bits':>9s} "
          f"{'T_eff/p0tilt':>13s} {'KL bits':>9s} {'KL to p0':>9s}")
    for a in alphas:
        print(f"   {a[1:]:>6s} {np.nanmean(per_target(a,'H_q_bits')):10.2f} "
              f"{np.nanmean(per_target(a,'T_eff_boltzmann')):12.4f} "
              f"{np.nanmean(per_target(a,'KL_bits_boltzmann')):9.3f} "
              f"{np.nanmean(per_target(a,'T_eff_p0tilt')):13.4f} "
              f"{np.nanmean(per_target(a,'KL_bits_p0tilt')):9.3f} "
              f"{np.nanmean(per_target(a,'KL_bits_to_p0')):9.3f}")

    print(f"\n=== P4  IS log q_theta A BETTER RANKER THAN THE OBJECTIVE IT WAS TRAINED ON? ===")
    print("   in-band = the objective's own best 1%, the only band a selector consumes.")
    print(f"   {'alpha':>6s} {'rho_ib -E':>10s} {'rho_ib q_th':>12s} {'rho_ib mf':>10s} "
          f"{'q-E: mean':>10s} {'CI95':>18s} {'W/L':>7s} | {'mf-E mean':>10s}")
    for a in alphas:
        e = per_target(a, "rho_inband_negE")
        q = per_target(a, "rho_inband_logq_theta")
        m = per_target(a, "rho_inband_logq_meanfield")
        # higher rho = better ranking; the paired difference is q - E, POSITIVE = better
        v = L.verdict(-(q - e), folds=folds, label=a)   # verdict is signed 'negative better'
        vm = L.verdict(-(m - e), folds=folds, label=a)
        print(f"   {a[1:]:>6s} {np.nanmean(e):10.3f} {np.nanmean(q):12.3f} "
              f"{np.nanmean(m):10.3f} {-v['mean']:+10.3f} "
              f"[{-v['ci95'][1]:+7.3f},{-v['ci95'][0]:+7.3f}] {v['W']:3d}/{v['L']:<3d} | "
              f"{-vm['mean']:+10.3f}")
    print("\n   THE LOAD-BEARING PAIRED TEST for P4: log q_theta MINUS the mean-field "
          "control, in-band, per target (positive = the circuit adds ranking skill the "
          "zero-correlation classical model does not):")
    for a in alphas:
        q = per_target(a, "rho_inband_logq_theta")
        m = per_target(a, "rho_inband_logq_meanfield")
        v = L.verdict(-(q - m), folds=folds, label=a)
        print(f"   {a[1:]:>6s} q_theta - meanfield {-v['mean']:+.3f} "
              f"[{-v['ci95'][1]:+.3f},{-v['ci95'][0]:+.3f}]  median {-v['median']:+.3f}  "
              f"W/L {v['W']}/{v['L']}  folds same sign {v['folds_same_sign']}/{v['n_folds']}")

    print("\n   global (not in-band), same three:")
    for a in alphas:
        print(f"   {a[1:]:>6s} -E {np.nanmean(per_target(a,'rho_global_negE')):+.3f}  "
              f"q_theta {np.nanmean(per_target(a,'rho_global_logq_theta')):+.3f}  "
              f"meanfield {np.nanmean(per_target(a,'rho_global_logq_meanfield')):+.3f}  "
              f"rho(log q, -E) {np.nanmean(per_target(a,'rho_logq_vs_negE')):+.3f}  "
              f"rho(log q_mf, -E) {np.nanmean(per_target(a,'rho_logqmf_vs_negE')):+.3f}")

    print(f"\n=== section 24  THE REWEIGHTING READOUT on a fixed 500-candidate pool "
          f"(picks of 75) ===")
    print(f"   {'alpha':>6s} {'arm':>20s} {'M':>7s} {'D':>7s} {'readout':>8s} "
          f"{'set_best':>9s} {'argmax RMSD':>12s} {'H(w) bits':>10s}")
    uni = [v["reweight_uniform"] for v in cells.values()]
    print(f"   {'-':>6s} {'uniform(control)':>20s} "
          f"{np.mean([u['M'] for u in uni]):7.3f} {np.mean([u['D'] for u in uni]):7.3f} "
          f"{np.mean([u['readout'] for u in uni]):8.3f} "
          f"{np.mean([u['set_best'] for u in uni]):9.3f} "
          f"{np.mean([u['argmax_rmsd'] for u in uni]):12.3f} "
          f"{np.mean([u['H_weights_bits'] for u in uni]):10.2f}")
    for a in alphas:
        for arm in ("vqe", "boltzmann", "meanfield", "p0tilt", "random_matched",
                    "vqe_priorpool", "boltzmann_priorpool", "meanfield_priorpool",
                    "random_matched_priorpool"):
            rr = [v[a]["reweight"][arm] for v in cells.values()
                  if a in v and arm in v[a]["reweight"]]
            if not rr:
                continue
            print(f"   {a[1:]:>6s} {arm:>20s} {np.mean([u['M'] for u in rr]):7.3f} "
                  f"{np.mean([u['D'] for u in rr]):7.3f} "
                  f"{np.mean([u['readout'] for u in rr]):8.3f} "
                  f"{np.mean([u['set_best'] for u in rr]):9.3f} "
                  f"{np.mean([u['argmax_rmsd'] for u in rr]):12.3f} "
                  f"{np.mean([u['H_weights_bits'] for u in rr]):10.2f}")


def local():
    cells = _cells("local")
    if not cells:
        print("no local cells yet")
        return
    pdbs = sorted({k.rsplit("_", 1)[0] for k in cells})
    folds = np.array([L.inst(p).fold for p in pdbs])
    arms = sorted({a for v in cells.values() for w in v["windows"].values() for a in w})
    meta = list(cells.values())[0]["_meta"]
    print(f"\n=== P3  THE LOCAL VQE, window w = {meta['w']} residues "
          f"({meta['sub_N']} configurations), budget {meta['budget']} "
          f"= {meta['budget']/meta['sub_N']:.2f}x the exhaustive cost ===")
    print(f"   n = {len(pdbs)} targets, windows averaged within target, seeds averaged "
          "within target.\n")

    def per_target(arm, field):
        out = []
        for p in pdbs:
            vals = []
            for k, v in cells.items():
                if k.rsplit("_", 1)[0] != p:
                    continue
                for w in v["windows"].values():
                    if arm in w and field in w[arm]:
                        vals.append(w[arm][field])
            out.append(np.mean(vals) if vals else np.nan)
        return np.asarray(out, float)

    print(f"{'arm':20s} {'M':>7s} {'D':>7s} {'readout':>8s} {'set_best':>9s} "
          f"{'argminRMSD':>11s} {'%certified':>11s} {'used':>7s}")
    for a in ["base_do_nothing"] + [x for x in arms if x != "base_do_nothing"]:
        M = per_target(a, "M"); D = per_target(a, "D")
        print(f"{a:20s} {np.nanmean(M):7.3f} {np.nanmean(D):7.3f} "
              f"{np.nanmean(per_target(a,'readout')):8.3f} "
              f"{np.nanmean(per_target(a,'set_best_ORACLE')):9.3f} "
              f"{np.nanmean(per_target(a,'argmin_rmsd_ORACLE')):11.3f} "
              f"{100*np.nanmean(per_target(a,'reaches_certified')):11.1f} "
              f"{np.nanmean(per_target(a,'used')):7.0f}")
    print("\n--- dominance margin of each local-VQE arm against the local classical set")
    Q = [a for a in arms if a.startswith("vqe")]
    C = [a for a in arms if a.startswith(("greedy", "anneal", "random", "metro",
                                          "exhaustive", "exact_boltz"))]
    marg = {}
    for a in Q:
        d = []
        for p in pdbs:
            pa = {"M": per_target(a, "M")[pdbs.index(p)],
                  "D": per_target(a, "D")[pdbs.index(p)]}
            cl = [{"M": per_target(c, "M")[pdbs.index(p)],
                   "D": per_target(c, "D")[pdbs.index(p)]} for c in C]
            d.append(L.eps_dominance(pa, cl))
        d = np.asarray(d, float)
        v = L.verdict(d, folds=folds, label=a)
        marg[a] = v
        print(f"{a:20s} mean {v['mean']:+7.3f} median {v['median']:+7.3f} "
              f"CI [{v['ci95'][0]:+7.3f},{v['ci95'][1]:+7.3f}] W/L {v['W']}/{v['L']}  "
              f"eps>0 {int((d > 0).sum())}/{len(pdbs)}")


def sched(m=75, pre=""):
    """T6: shots-vs-iterations, restarts, and THE ENTANGLEMENT CONTROL.

    The product-ansatz arm (`mps2fn`, the identical circuit with the CNOTs removed) is a
    classical variational sampler: sampling and score-function gradients of a fully
    factorised Bernoulli model are trivially classical.  It is therefore merged into the
    CLASSICAL attainable set and the mps2f arms are re-scored against it.
    """
    sc = _cells("sched")
    pc = _cells("pareto")
    if not sc:
        print("no sched cells yet")
        return
    keys = sorted(set(sc) & set(pc))
    pdbs = sorted({k.rsplit("_", 1)[0] for k in keys})
    folds = np.array([L.inst(p).fold for p in pdbs])
    bts = {p: [sc[k] for k in keys if k.rsplit("_", 1)[0] == p] for p in pdbs}
    btp = {p: [pc[k] for k in keys if k.rsplit("_", 1)[0] == p] for p in pdbs}
    sarms = sorted({a for v in sc.values() for a in v if not a.startswith("_")})
    print(f"\n=== T6  SHOTS / RESTARTS / THE ENTANGLEMENT CONTROL, m = {m}, "
          f"n = {len(pdbs)} targets, seeds {sorted({k.rsplit('_',1)[1] for k in keys})} ===\n")
    print(f"{'arm':22s} {'M':>7s} {'D':>7s} {'readout':>8s} {'set_best':>9s} "
          f"{'distinct':>9s} {'iters':>6s}")
    spts = {}
    for a in sarms:
        pts = {p: _agg_point(bts[p], a, m, pre) for p in pdbs}
        spts[a] = pts
        v = [pts[p] for p in pdbs if pts[p]]
        it = np.mean([r[a].get("iters", np.nan) for r in sum(bts.values(), []) if a in r])
        print(f"{a:22s} {np.mean([x['M'] for x in v]):7.3f} "
              f"{np.mean([x['D'] for x in v]):7.3f} "
              f"{np.mean([x['readout'] for x in v]):8.3f} "
              f"{np.mean([x['set_best'] for x in v]):9.3f} "
              f"{np.mean([x['n_distinct'] for x in v]):9.0f} {it:6.0f}")
    # merged Pareto: the quantum arms of q_pareto against classical + product-ansatz VI
    C = [a for a in P.CLASSICAL_MATCHED]
    Q = list(P.QUANTUM_ARMS)
    prod = [a for a in sarms if a.startswith("product")]
    print(f"\n--- dominance margin of the mps2f arms against classical ({len(C)}) PLUS the "
          f"product-ansatz classical variational sampler ({len(prod)})")
    print(f"{'arm':22s} {'mean':>8s} {'median':>8s} {'CI95':>20s} {'W/L':>7s} {'outside':>8s}")
    for a in Q:
        marg = []
        for p in pdbs:
            qp = _agg_point(btp[p], a, m, pre)
            cl = [x for x in ([_agg_point(btp[p], c, m, pre) for c in C]
                              + [spts[z][p] for z in prod]) if x]
            marg.append(L.eps_dominance(qp, cl) if qp else np.nan)
        marg = np.asarray(marg, float)
        if not np.isfinite(marg).all():
            continue
        v = L.verdict(marg, folds=folds, label=a)
        print(f"{a:22s} {v['mean']:+8.3f} {v['median']:+8.3f} "
              f"[{v['ci95'][0]:+7.3f},{v['ci95'][1]:+7.3f}] "
              f"{int((marg > 0).sum()):3d}/{len(pdbs):<4d}")
    print("\n--- head-to-head: entangled mps2f vs the CNOT-free product ansatz at matched "
          "alpha, budget, estimator, optimiser and seed (paired over targets)")
    for al, pa in (("vqe_a0.05", "product_a0.05"), ("vqe_a0.25", "product_a0.25"),
                   ("vqe_a1.00", "product_a1.0")):
        if pa not in spts:
            continue
        dM, dD, dR = [], [], []
        for p in pdbs:
            q = _agg_point(btp[p], al, m, pre)
            r = spts[pa][p]
            if q and r:
                dM.append(q["M"] - r["M"]); dD.append(q["D"] - r["D"])
                dR.append(q["readout"] - r["readout"])
        vM = L.verdict(np.asarray(dM), folds=folds, label=al)
        vD = L.verdict(np.asarray(dD), folds=folds, label=al)
        vR = L.verdict(np.asarray(dR), folds=folds, label=al)
        print(f"   {al} - {pa}:  dM {vM['mean']:+.3f} [{vM['ci95'][0]:+.3f},"
              f"{vM['ci95'][1]:+.3f}] W/L {vM['W']}/{vM['L']}   dD {vD['mean']:+.3f} "
              f"[{vD['ci95'][0]:+.3f},{vD['ci95'][1]:+.3f}]   dreadout {vR['mean']:+.3f} "
              f"[{vR['ci95'][0]:+.3f},{vR['ci95'][1]:+.3f}] W/L {vR['W']}/{vR['L']}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "pareto"
    if mode == "pareto":
        pareto(int(sys.argv[2]) if len(sys.argv) > 2 else 75,
               sys.argv[3] if len(sys.argv) > 3 else "")
    elif mode == "novelty":
        novelty()
    elif mode == "theory":
        theory()
    elif mode == "repr":
        repr_report()
    elif mode == "local":
        local()
    elif mode == "sched":
        sched(int(sys.argv[2]) if len(sys.argv) > 2 else 75,
              sys.argv[3] if len(sys.argv) > 3 else "")
    else:
        raise SystemExit(mode)
