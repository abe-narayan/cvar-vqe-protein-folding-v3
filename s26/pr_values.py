"""PR lane (Presentation), Sprint 26: the registry of every number that appears on a slide
or in the speaker notes of `vqe_research_overview.pptx`.

Each entry is read from the artefact it names at build time, so the deck is regenerable and
no number is typed by hand into a slide, a note or a figure. The registry is written out as
`s26/pr_values.json` by `s26/pr_build_deck.py` so the record of every value and its path is a
file, not prose.

Rule 1 of the lane contract: nothing here opens `s9/final_report.json`,
`results/benchmark_manifest.json` or any benchmark record. The one benchmark number the deck
carries (claim C06 of `s26/EXAMINATION.md`) is typed from the claim ledger and labelled
NAMED_NOT_OPENED; its two component means are asserted by `tests/test_pipeline.py`.

Bases are named on every RMSD entry: `built_chain` (production, `rmsd_arm`), `point_cloud`
(`rmsd_avg`, an intermediate), `relaxed` (`rmsd_full`), `single_window` (a pool member),
`s8_selection` (the S25 quantum instrument: one medoid out of 128).
"""
from __future__ import annotations

import json
import math
import os
from collections import OrderedDict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROD_CACHE = os.path.join(ROOT, "bench_results", "cache", "1fc9f2dcf489e2fb")


def _j(rel):
    with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
        return json.load(fh)


def _slope(rows):
    """Least-squares slope of log2(var_g0) against n over the given rows (the S25 fit)."""
    n = np.array([r["n"] for r in rows], float)
    y = np.log2(np.array([r["var_g0"] for r in rows], float))
    return float(np.polyfit(n, y, 1)[0])


def load_values():
    """Return an OrderedDict token -> dict(value, path, basis, status, note)."""
    V = OrderedDict()

    def put(token, value, path, status="SOURCED", basis="", note=""):
        V[token] = dict(value=value, path=path, status=status, basis=basis, note=note)

    # ------------------------------------------------------------------ headline (C01..C03)
    base = _j("bench_results/baseline_tuning126.json")["science"]
    put("ARM_MEAN", base["rmsd_arm"]["mean"], "bench_results/baseline_tuning126.json :: science/rmsd_arm/mean",
        basis="built_chain", note="claim C01; also s26/results/e_reproduce.json :: summary/rmsd_arm/mean_fresh")
    put("ARM_MEDIAN", base["rmsd_arm"]["median"], "bench_results/baseline_tuning126.json :: science/rmsd_arm/median", basis="built_chain")
    put("ARM_MIN", base["rmsd_arm"]["min"], "bench_results/baseline_tuning126.json :: science/rmsd_arm/min", basis="built_chain")
    put("ARM_MAX", base["rmsd_arm"]["max"], "bench_results/baseline_tuning126.json :: science/rmsd_arm/max", basis="built_chain")
    put("ARM_FRAC2", base["rmsd_arm"]["frac_under_2.0"], "bench_results/baseline_tuning126.json :: science/rmsd_arm/frac_under_2.0", basis="built_chain")
    put("AVG_MEAN", base["rmsd_avg"]["mean"], "bench_results/baseline_tuning126.json :: science/rmsd_avg/mean",
        basis="point_cloud", note="claim C02; an intermediate, not a structure")
    put("FULL_MEAN", base["rmsd_full"]["mean"], "bench_results/baseline_tuning126.json :: science/rmsd_full/mean",
        basis="relaxed", note="claim C03; 9KAR ends above the 1000 kcal/mol gate and is inside this mean")
    put("POOL_BEST", base["pool_best"]["mean"], "bench_results/baseline_tuning126.json :: science/pool_best/mean",
        basis="single_window", note="claim C07, ORACLE: best member of the K=500 pool")
    put("TOPM_BEST", base["top_m_best"]["mean"], "bench_results/baseline_tuning126.json :: science/top_m_best/mean",
        basis="single_window", note="ORACLE: best member of the shipped top-75")
    put("SHIPPED_ARGMIN", base["shipped"]["mean"], "bench_results/baseline_tuning126.json :: science/shipped/mean",
        basis="single_window", note="the shipped distogram argmin (one window)")
    put("N_TARGETS", int(base["n"]), "bench_results/baseline_tuning126.json :: science/n")

    rep = _j("s26/results/e_reproduce.json")
    rows = rep["rows"]
    arm = np.array([r["rmsd_arm"] for r in rows], float)
    put("ARM_FRAC3", float(np.mean(arm < 3.0)), "s26/results/e_reproduce.json :: rows[*]/rmsd_arm (fraction under 3.0 A)",
        status="DERIVED", basis="built_chain")
    worst = rows[int(np.argmax(arm))]
    put("ARM_WORST_PDB", worst["pdb"], "s26/results/e_reproduce.json :: rows[argmax rmsd_arm]/pdb", status="DERIVED", basis="built_chain")
    best = rows[int(np.argmin(arm))]
    put("ARM_BEST_PDB", best["pdb"], "s26/results/e_reproduce.json :: rows[argmin rmsd_arm]/pdb", status="DERIVED", basis="built_chain")
    put("REPRO_MAX_DISAGREE", rep["summary"]["rmsd_arm"]["max_abs_disagreement"],
        "s26/results/e_reproduce.json :: summary/rmsd_arm/max_abs_disagreement", note="fresh instrument vs stored record, 126 targets")
    put("GAP_ARM_MINUS_AVG", rep["summary"]["projection_gap_arm_minus_avg"]["mean"],
        "s26/results/e_reproduce.json :: summary/projection_gap_arm_minus_avg/mean", basis="built_chain minus point_cloud")
    put("BOND_AVG_MEAN", rep["summary"]["virtual_bond"]["avg_ca_mean_bond"], "s26/results/e_reproduce.json :: summary/virtual_bond/avg_ca_mean_bond",
        basis="point_cloud", note="mean virtual CA-CA bond of the raw average")
    put("BOND_AVG_WORST", rep["summary"]["virtual_bond"]["avg_ca_worst_bond"], "s26/results/e_reproduce.json :: summary/virtual_bond/avg_ca_worst_bond",
        basis="point_cloud")

    # ------------------------------------------------------------------ the two overlay targets
    for pdb, tag in (("1S9Z", "T030"), ("9KAR", "HARD")):
        rec = _j(f"bench_results/cache/1fc9f2dcf489e2fb/{pdb}.json")
        put(f"{tag}_PDB", pdb, f"bench_results/cache/1fc9f2dcf489e2fb/{pdb}.json :: pdb")
        put(f"{tag}_N", int(rec["n"]), f"bench_results/cache/1fc9f2dcf489e2fb/{pdb}.json :: n")
        put(f"{tag}_ARM", rec["rmsd_arm"], f"bench_results/cache/1fc9f2dcf489e2fb/{pdb}.json :: rmsd_arm", basis="built_chain")
        put(f"{tag}_POOL_BEST", rec["pool_best"], f"bench_results/cache/1fc9f2dcf489e2fb/{pdb}.json :: pool_best", basis="single_window", note="ORACLE")
        put(f"{tag}_TOPM_BEST", rec["top_m_best"], f"bench_results/cache/1fc9f2dcf489e2fb/{pdb}.json :: top_m_best", basis="single_window", note="ORACLE")
        put(f"{tag}_AVG", rec["rmsd_avg"], f"bench_results/cache/1fc9f2dcf489e2fb/{pdb}.json :: rmsd_avg", basis="point_cloud")
    trace = _j("s26/results/e_trace_1S9Z.json")
    put("T030_N_WINDOWS", int(_j("bench_results/cache/1fc9f2dcf489e2fb/1S9Z.json")["n_windows"]),
        "bench_results/cache/1fc9f2dcf489e2fb/1S9Z.json :: n_windows", note="windows in the leakage-safe universe of T030")

    # ------------------------------------------------------------------ benchmark (C04..C06): typed from the claim ledger
    put("BENCH_DELTA", 0.0103, "s9/final_report.json (NOT OPENED, Rule 1); claim C06 of s26/EXAMINATION.md; docs/FINDINGS.md:4505; README.md:14",
        status="NAMED_NOT_OPENED", basis="benchmark synthesis minus shipped baseline, paired, n=60")
    put("BENCH_CI_LO", -0.1596, "same as BENCH_DELTA", status="NAMED_NOT_OPENED")
    put("BENCH_CI_HI", 0.1803, "same as BENCH_DELTA", status="NAMED_NOT_OPENED")
    put("BENCH_WL", "31W/29L", "same as BENCH_DELTA", status="NAMED_NOT_OPENED")
    put("BENCH_FULL", 2.9610, "s9/final_report.json :: dist/full/mean (NOT OPENED); asserted within 5e-4 by tests/test_pipeline.py:338-347 (passed, s26/TEST_RUN.md)",
        status="SOURCED_BY_TEST", basis="benchmark, 60 targets")
    put("BENCH_SHIPPED", 2.9507, "s9/final_report.json :: dist/shipped/mean (NOT OPENED); same tests", status="SOURCED_BY_TEST", basis="benchmark, 60 targets")

    # ------------------------------------------------------------------ controls and ORACLE ladder (C08..C17)
    ps = _j("s25/results/phys_suite.json")
    assert ps["basis"] == "point_cloud"
    put("RANDOM75", ps["random_null_rank"]["mean"], "s25/results/phys_suite.json :: random_null_rank/mean", basis="point_cloud",
        note="claim C08; the artefact declares basis = point_cloud; its production counterpart is AVG_MEAN")
    put("SUITE_INCUMBENT", ps["incumbent"], "s25/results/phys_suite.json :: incumbent", basis="point_cloud")
    cfg = {ps["configs_rank"][k]["name"]: ps["configs_rank"][k]["mean"] for k in ps["configs_rank"]}
    for name, tok in (("Distogram", "SUITE_DIST"), ("AMBER+Distogram", "SUITE_AMBER_DIST"), ("Legacy+Distogram", "SUITE_LEGACY_DIST"),
                      ("Legacy+AMBER+Distogram", "SUITE_LAD"), ("Legacy+AMBER", "SUITE_LEGACY_AMBER"), ("Legacy", "SUITE_LEGACY"), ("AMBER", "SUITE_AMBER")):
        put(tok, cfg[name], f"s25/results/phys_suite.json :: configs_rank/*/mean (name == {name})", basis="point_cloud", note="claim C29")
    c1 = ps["controls_rank"]["1"]["vs_random"]
    c2 = ps["controls_rank"]["2"]["vs_random"]
    put("LEGACY_VS_RANDOM", c1["effect"], "s25/results/phys_suite.json :: controls_rank/1/vs_random/effect", basis="point_cloud", note="claim C30")
    put("LEGACY_VS_RANDOM_CI", c1["ci95_fold"], "s25/results/phys_suite.json :: controls_rank/1/vs_random/ci95_fold", basis="point_cloud")
    put("LEGACY_VS_RANDOM_FOLDS", int(c1["folds_same_sign"]), "s25/results/phys_suite.json :: controls_rank/1/vs_random/folds_same_sign")
    put("LEGACY_VS_RANDOM_MDE", c1["mde"], "s25/results/phys_suite.json :: controls_rank/1/vs_random/mde")
    put("AMBER_VS_RANDOM", c2["effect"], "s25/results/phys_suite.json :: controls_rank/2/vs_random/effect", basis="point_cloud", note="claim C30")
    put("AMBER_VS_RANDOM_CI", c2["ci95_fold"], "s25/results/phys_suite.json :: controls_rank/2/vs_random/ci95_fold", basis="point_cloud")
    put("AMBER_VS_RANDOM_FOLDS", int(c2["folds_same_sign"]), "s25/results/phys_suite.json :: controls_rank/2/vs_random/folds_same_sign")
    put("AMBER_VS_RANDOM_MDE", c2["mde"], "s25/results/phys_suite.json :: controls_rank/2/vs_random/mde")

    lb = _j("results/summary/leaderboard.json")
    assert lb["basis"] == "built_chain_bb"
    lbrows = {r["configuration"]: r for r in lb["rows"]}
    put("LB_PROD_MEAN", lbrows["production"]["mean"], "results/summary/leaderboard.json :: rows[production]/mean", basis="built_chain",
        note="the results-lab rebuild (re-projected, PDB round trip); production cache is ARM_MEAN")
    for cfgname in ("distogram", "amber_distogram", "legacy_distogram", "legacy_amber_distogram", "legacy_amber", "legacy", "amber"):
        if cfgname in lbrows:
            put("LB_" + cfgname.upper(), lbrows[cfgname]["mean"], f"results/summary/leaderboard.json :: rows[{cfgname}]/mean", basis="built_chain",
                note="built-chain mean of the same seven-configuration row (L28 item 3)")

    s14 = _j("s12/results/s14_ladder.json")["rows"]["L0_constant_helix"]
    put("HELIX", s14["mean"], "s12/results/s14_ladder.json :: rows/L0_constant_helix/mean", basis="built_chain",
        note="claim C09; zero-information control, constant alpha-helix built at ideal geometry")
    put("HELIX_MINUS_ARM", s14["mean"] - base["rmsd_arm"]["mean"], "s12/results/s14_ladder.json :: rows/L0_constant_helix/mean minus bench_results/baseline_tuning126.json :: science/rmsd_arm/mean",
        status="DERIVED", basis="built_chain", note="ledger L14 quotes the same paired value +0.8500")
    te = _j("s13/results/tors_eval.json")["reports"]["a_pepPos"]["build"]
    put("TORS", te["mean"], "s13/results/tors_eval.json :: reports/a_pepPos/build/mean", basis="built_chain", note="claim C10; sequence-only torsion predictor")
    put("TORS_MINUS_ARM", te["mean"] - base["rmsd_arm"]["mean"], "s13/results/tors_eval.json :: reports/a_pepPos/build/mean minus ARM_MEAN",
        status="DERIVED", basis="built_chain", note="ledger L14: +0.5557, SE 0.1318, MDE 0.3694, 37W/89L")

    pl = _j("s24/results/priorladder.json")["rows"]
    m0 = float(np.mean([r["MASS0.0"] for r in pl])); m01 = float(np.mean([r["MASS0.1"] for r in pl])); m1 = float(np.mean([r["MASS1.0"] for r in pl]))
    put("PRIOR_GAMMA0", m0, "s24/results/priorladder.json :: rows[*]/MASS0.0 (mean over 126)", status="DERIVED", basis="point_cloud", note="equals AVG_MEAN")
    put("PRIOR_PERFECT", m1, "s24/results/priorladder.json :: rows[*]/MASS1.0 (mean over 126)", status="DERIVED", basis="point_cloud", note="claim C11, ORACLE")
    put("PRIOR_SLOPE", (m01 - m0) / 0.1, "s24/results/priorladder.json :: (mean MASS0.1 - mean MASS0.0) / 0.1", status="DERIVED", basis="point_cloud", note="claim C15")
    put("PRIOR_GAMMA_FOR_3A", (m0 - 3.0) / (m0 - m01) * 0.1, "s24/results/priorladder.json :: linear interpolation on the first rung", status="DERIVED", basis="point_cloud", note="claim C16")
    ed = _j("s23/results/errdecomp.json")["rows"]
    put("COMMON_MODE", float(np.mean([r["f_common"] for r in ed])), "s23/results/errdecomp.json :: rows[*]/f_common (mean over 126)", status="DERIVED", note="claim C17")
    put("UNIVERSE_BEST", _j("bench_results/recon_library_saturation.json")["universe_best"]["small"],
        "bench_results/recon_library_saturation.json :: universe_best/small", basis="single_window", note="claim C12, ORACLE: best window in the whole library")
    put("TORSION_CEILING", _j("s13/results/ceiling_report.json")["k"]["4"]["descent"]["mean"],
        "s13/results/ceiling_report.json :: k/4/descent/mean", basis="built_chain (torsion-space rebuild)", note="claim C13, ORACLE")
    put("DISTGEO_TRUE", _j("s15/results/distgeo.json")["arms"]["ORACLE_true_distances"]["mean"],
        "s15/results/distgeo.json :: arms/ORACLE_true_distances/mean", basis="distance-geometry fit", note="claim C14, ORACLE")

    # ------------------------------------------------------------------ phi (C26, derived by the Adversary, L32)
    phi = _j("s26/results/a_c26_phi_mae.json")["summary"]
    put("PHI_SEQ", phi["arms"]["p_grid"]["pooled_mae_phi_deg"], "s26/results/a_c26_phi_mae.json :: summary/arms/p_grid/pooled_mae_phi_deg",
        note="claim C26 re-derived (L32); full sequence context")
    put("PHI_BLIND", phi["arms"]["n_marg"]["pooled_mae_phi_deg"], "s26/results/a_c26_phi_mae.json :: summary/arms/n_marg/pooled_mae_phi_deg",
        note="sequence-blind corpus marginal")
    put("PHI_N_RES", int(phi["n_residues_phi"]), "s26/results/a_c26_phi_mae.json :: summary/n_residues_phi")

    # ------------------------------------------------------------------ the quantum component (C18..C25, C31..C33)
    qv = _j("s25/results/q_verify.json")["results"]
    put("Q_QUBITS", int(qv["core/pipeline Config.vqe_qubits"]), "s25/results/q_verify.json :: results/'core/pipeline Config.vqe_qubits'")
    put("Q_LAYERS", int(qv["core/pipeline Config.vqe_layers"]), "s25/results/q_verify.json :: results/'core/pipeline Config.vqe_layers'")
    put("Q_PARAMS", int(qv["parameters = layers * qubits"]), "s25/results/q_verify.json :: results/'parameters = layers * qubits'")
    put("Q_ITERS", int(qv["core/pipeline Config.vqe_iters"]), "s25/results/q_verify.json :: results/'core/pipeline Config.vqe_iters'")
    put("Q_DIM", int(qv["hypothesis-set size 2**vqe_qubits"]), "s25/results/q_verify.json :: results/'hypothesis-set size 2**vqe_qubits'")
    put("Q_SV_ERR", qv["max |p_statevector - p_dense|"], "s25/results/q_verify.json :: results/'max |p_statevector - p_dense|'", note="claim C21 (stored as a string)")
    put("Q_PS_COS", qv["cos(param-shift, exact FD), min over 12"], "s25/results/q_verify.json :: results/'cos(param-shift, exact FD), min over 12'", note="claim C22")
    put("Q_PS_REL", qv["rel |g_ps - g_fd| / |g_fd|, max over 12"], "s25/results/q_verify.json :: results/'rel |g_ps - g_fd| / |g_fd|, max over 12'", note="claim C23")
    put("Q_FE_REL", qv["rel err of free-energy grad, max over 4"], "s25/results/q_verify.json :: results/'rel err of free-energy grad, max over 4'", note="claim C23")
    put("Q_SAMPLING", qv["sampling calls inside run_cvar_vqe/free_energy"], "s25/results/q_verify.json :: results/'sampling calls inside run_cvar_vqe/free_energy'")
    put("Q_RNG", qv["rng use inside run_cvar_vqe"], "s25/results/q_verify.json :: results/'rng use inside run_cvar_vqe'")
    put("Q_CELLS", int(qv["n_cells"]), "s25/results/q_verify.json :: results/n_cells", note="claim C32")
    put("Q_VIOLATIONS", int(qv["SUBSET-hood violations  <- THE THEOREM"]), "s25/results/q_verify.json :: results/'SUBSET-hood violations  <- THE THEOREM'", note="claim C32")
    put("Q_EQUALITY", qv["full-support cells with EXACT value equality"], "s25/results/q_verify.json :: results/'full-support cells with EXACT value equality'", note="claim C32")
    put("Q_TAIL_COS", qv["cos(exact-expectation, TAIL baseline)"], "s25/results/q_verify.json :: results/'cos(exact-expectation, TAIL baseline)'", note="claim C24 (the sourced member)")
    put("Q_LFO", qv["VQE_LFO (alpha, T) per fold"], "s25/results/q_verify.json :: results/'VQE_LFO (alpha, T) per fold'")
    put("Q_ALPHA1_FOLDS", qv["folds running alpha == 1.0 (NO tail constraint)"], "s25/results/q_verify.json :: results/'folds running alpha == 1.0 (NO tail constraint)'")
    put("Q_MPS_ERR", qv["max |p_MPS - p_dense| over 24 configs"], "s25/results/q_verify.json :: results/'max |p_MPS - p_dense| over 24 configs'", note="generation lane; not the selector")

    qa = _j("s25/results/q_alpha.json")["results"]
    put("Q_SHARE_NO_TAIL", qa["share_of_targets_with_no_tail_constraint"], "s25/results/q_alpha.json :: results/share_of_targets_with_no_tail_constraint")
    va = qa["vs_no_circuit"]["VQE_LFO - argmin (shipped)"]
    put("Q_VS_ARGMIN", va["mean"], "s25/results/q_alpha.json :: results/vs_no_circuit/'VQE_LFO - argmin (shipped)'/mean", basis="s8_selection")
    put("Q_VS_ARGMIN_SE", va["se"], "s25/results/q_alpha.json :: .../'VQE_LFO - argmin (shipped)'/se")
    put("Q_VS_ARGMIN_MDE", va["mde"], "s25/results/q_alpha.json :: .../'VQE_LFO - argmin (shipped)'/mde")
    put("Q_VS_ARGMIN_X", va["eff_over_mde"], "s25/results/q_alpha.json :: .../'VQE_LFO - argmin (shipped)'/eff_over_mde")
    put("Q_VS_ARGMIN_W", int(va["W"]), "s25/results/q_alpha.json :: .../'VQE_LFO - argmin (shipped)'/W")
    put("Q_VS_ARGMIN_L", int(va["L"]), "s25/results/q_alpha.json :: .../'VQE_LFO - argmin (shipped)'/L")
    put("Q_VS_ARGMIN_T", int(va["ties"]), "s25/results/q_alpha.json :: .../'VQE_LFO - argmin (shipped)'/ties")
    put("Q_VS_ARGMIN_FOLDS", int(va["folds_same_sign"]), "s25/results/q_alpha.json :: .../'VQE_LFO - argmin (shipped)'/folds_same_sign")
    put("Q_VS_ARGMIN_VERDICT", va["verdict"], "s25/results/q_alpha.json :: .../'VQE_LFO - argmin (shipped)'/verdict")
    vb = qa["vs_no_circuit"]["VQE_LFO - Boltzmann T=0.3"]
    put("Q_VS_BOLTZ", vb["mean"], "s25/results/q_alpha.json :: results/vs_no_circuit/'VQE_LFO - Boltzmann T=0.3'/mean", basis="s8_selection")
    ex = qa["expressivity_vs_boltzmann"]["circuit - exact Boltzmann @T=0.3"]
    put("Q_CIRC_VS_GIBBS", ex["mean"], "s25/results/q_alpha.json :: results/expressivity_vs_boltzmann/'circuit - exact Boltzmann @T=0.3'/mean", basis="s8_selection", note="claim C19")
    put("Q_CIRC_VS_GIBBS_X", ex["eff_over_mde"], "s25/results/q_alpha.json :: .../'circuit - exact Boltzmann @T=0.3'/eff_over_mde", note="claim C19")
    put("Q_RHO_H", qa["cell_correlations"]["H"], "s25/results/q_alpha.json :: results/cell_correlations/H", note="claim C31")
    put("Q_OFF_CURVE", qa["circuit_off_curve"]["mean_residual_vqe"], "s25/results/q_alpha.json :: results/circuit_off_curve/mean_residual_vqe")

    qg = _j("s25/results/q_gibbs.json")["results"]
    lad = [d for d in qg["divergence_ladder"] if d["T"] == 0.3][0]["ladder"]
    dep = lad["50 Adam steps (DEPLOYED)"]
    put("Q_KL_NATS", dep["kl_nats"], "s25/results/q_gibbs.json :: results/divergence_ladder/[T=0.3]/ladder/'50 Adam steps (DEPLOYED)'/kl_nats", note="claim C18")
    put("Q_TV", dep["tv"], "s25/results/q_gibbs.json :: .../'50 Adam steps (DEPLOYED)'/tv")
    put("Q_H_TRAINED", dep["entropy_bits"], "s25/results/q_gibbs.json :: .../'50 Adam steps (DEPLOYED)'/entropy_bits")
    put("Q_H_GIBBS", lad["the Gibbs optimum itself"]["entropy_bits"], "s25/results/q_gibbs.json :: .../'the Gibbs optimum itself'/entropy_bits")
    tc = {round(d["T"], 2): d for d in qg["training_control"]}
    put("Q_GAP_MIN", min(d["frac_gap_closed"] for d in tc.values()), "s25/results/q_gibbs.json :: results/training_control/*/frac_gap_closed (min over T)", status="DERIVED", note="claim C33")
    put("Q_GAP_MAX", max(d["frac_gap_closed"] for d in tc.values()), "s25/results/q_gibbs.json :: results/training_control/*/frac_gap_closed (max over T)", status="DERIVED", note="claim C33")
    put("Q_GAP_T03", tc[0.3]["frac_gap_closed"], "s25/results/q_gibbs.json :: results/training_control/[T=0.3]/frac_gap_closed")
    put("Q_BEATS_BEST_OF_N", all(d["beats_best_of_N"] for d in tc.values()), "s25/results/q_gibbs.json :: results/training_control/*/beats_best_of_N", status="DERIVED")
    put("Q_BEST_OF_N", int(tc[0.3]["N"]), "s25/results/q_gibbs.json :: results/training_control/[T=0.3]/N")
    sp = qg["spectrum_target_independence"]
    put("Q_LADDER_WORST_FRAC", sp["worst_frac_of_range"], "s25/results/q_gibbs.json :: results/spectrum_target_independence/worst_frac_of_range", note="claim C20")
    put("Q_LADDER_N_CHECKED", int(sp["n_checked"]), "s25/results/q_gibbs.json :: results/spectrum_target_independence/n_checked")

    qp = _j("s25/results/q_plateau.json")["results"]
    for cell, tok in (("linear_alpha1_T0", "SLOPE_LIN"), ("cvar_alpha025_T0", "SLOPE_A025_T0"), ("cvar_alpha01_T0", "SLOPE_A01_T0"),
                      ("deployed_a1_T03", "SLOPE_A1_T03"), ("deployed_a025_T03", "SLOPE_A025_T03")):
        put(tok, qp[cell]["log2_slope_per_qubit"], f"s25/results/q_plateau.json :: results/{cell}/log2_slope_per_qubit", note="claim C25; n = 4..13, depth 3, theta ~ N(0, 0.6^2)")
    ns = sorted({r["n"] for r in qp["linear_alpha1_T0"]["rows"]})
    put("SWEEP_N_MIN", int(min(ns)), "s25/results/q_plateau.json :: results/linear_alpha1_T0/rows[*]/n (min)")
    put("SWEEP_N_MAX", int(max(ns)), "s25/results/q_plateau.json :: results/linear_alpha1_T0/rows[*]/n (max)")
    draws = [r["n_theta"] for r in qp["linear_alpha1_T0"]["rows"]]
    put("SWEEP_DRAWS_MIN", int(min(draws)), "s25/results/q_plateau.json :: results/linear_alpha1_T0/rows[*]/n_theta (min)")
    put("SWEEP_DRAWS_MAX", int(max(draws)), "s25/results/q_plateau.json :: results/linear_alpha1_T0/rows[*]/n_theta (max)")
    ds = {r["layers"]: r["var_g0"] for r in qp["depth_sweep_n7"]}
    put("DEPTH_L1", ds[1], "s25/results/q_plateau.json :: results/depth_sweep_n7/[layers=1]/var_g0")
    put("DEPTH_L3", ds[3], "s25/results/q_plateau.json :: results/depth_sweep_n7/[layers=3]/var_g0")
    put("DEPTH_L4", ds[4], "s25/results/q_plateau.json :: results/depth_sweep_n7/[layers=4]/var_g0")
    put("DEPTH_L12", ds[12], "s25/results/q_plateau.json :: results/depth_sweep_n7/[layers=12]/var_g0")
    rt = qp["nonlinearity_variance_ratio"]
    put("RATIO_N7_A025", rt["7"][0], "s25/results/q_plateau.json :: results/nonlinearity_variance_ratio/7[0]", note="Var[grad CVaR_0.25] / Var[grad mean] at n = 7")
    put("RATIO_N13_A025", rt["13"][0], "s25/results/q_plateau.json :: results/nonlinearity_variance_ratio/13[0]")
    gk = _j("s13/results/geo_kernel.json")["scaling"]["lay8_v_w1"]
    put("GEO_BASE_DEPTH8", gk["decay_base"], "s13/results/geo_kernel.json :: scaling/lay8_v_w1/decay_base", note="the S13 depth-8 decay base, the 2-design-like rate")
    put("GEO_SLOPE_DEPTH8", math.log2(gk["decay_base"]), "s13/results/geo_kernel.json :: log2(scaling/lay8_v_w1/decay_base)", status="DERIVED")
    gp = _j("s13/results/geo_pauli.json")["cells"]
    ratio = [c["legacy"]["ratio_meas_over_pred_exact"] for c in gp]
    put("PAULI_RATIO_LEGACY", float(np.median(ratio)), "s13/results/geo_pauli.json :: cells[*]/legacy/ratio_meas_over_pred_exact (median)", status="DERIVED",
        note="measured over predicted gradient variance, per-string spectrum prediction, %d cells" % len(ratio))
    put("PAULI_N_CELLS", len(ratio), "s13/results/geo_pauli.json :: len(cells)")

    # ------------------------------------------------------------------ A2 and A4 (L27, L35)
    qd = _j("s26/results/q_dla.json")["results"]
    put("DLA_N7_L1", int(qd["fixed"]["n7_L1"]["dim"]), "s26/results/q_dla.json :: results/fixed/n7_L1/dim")
    put("DLA_N7_L2", int(qd["fixed"]["n7_L2"]["dim"]), "s26/results/q_dla.json :: results/fixed/n7_L2/dim")
    put("DLA_N7_L3", int(qd["fixed"]["n7_L3"]["dim"]), "s26/results/q_dla.json :: results/fixed/n7_L3/dim")
    put("DLA_SO128", int(qd["fixed"]["n7_L3"]["dim_so"]), "s26/results/q_dla.json :: results/fixed/n7_L3/dim_so")
    put("DLA_POOL_V_N7", int(qd["pools"]["V_n7"]["dim"]), "s26/results/q_dla.json :: results/pools/V_n7/dim")
    put("DLA_POOL_L2_N7", int(qd["pools"]["L2_n7"]["dim"]), "s26/results/q_dla.json :: results/pools/L2_n7/dim")
    agree = [v["agree"] for k, v in qd["numeric"].items() if "agree" in v]
    put("DLA_NUMERIC_AGREE", f"{sum(agree)}/{len(agree)}", "s26/results/q_dla.json :: results/numeric/*/agree", status="DERIVED")
    ad = qd["adapt_sets"]
    a1 = [v for k, v in ad.items() if v["alpha"] == 1.0]
    put("ADAPT_A1_DIM", max(v["ladder"][-1]["dim"] for v in a1), "s26/results/q_dla.json :: results/adapt_sets/[alpha=1]/ladder[-1]/dim (max over pools and optimisers)", status="DERIVED")
    a025L2 = [v for k, v in ad.items() if v["alpha"] == 0.25 and v["pool"] == "L2"]
    put("ADAPT_A025_L2_DIM", max(v["ladder"][-1]["dim"] for v in a025L2), "s26/results/q_dla.json :: results/adapt_sets/[alpha=0.25, pool L2]/ladder[-1]/dim (max)", status="DERIVED")

    qvar = _j("s26/results/q_var.json")["results"]
    grown = qvar["grown"]
    def matched(key):
        return [r for r in grown[key]["rows"] if r.get("P_matched")]
    put("A4_V_LIN", _slope(matched("V:adam_best:linear_alpha1_T0")), "s26/results/q_var.json :: results/grown/'V:adam_best:linear_alpha1_T0'/rows[P_matched] (log2 slope)", status="DERIVED", note="ledger L35: -0.079")
    put("A4_L2_LIN", _slope(matched("L2:adam_best:linear_alpha1_T0")), "s26/results/q_var.json :: results/grown/'L2:adam_best:linear_alpha1_T0'/rows[P_matched]", status="DERIVED", note="L35: +0.006")
    put("A4_V_A1_T03", _slope(matched("V:adam_best:deployed_a1_T03")), "s26/results/q_var.json :: results/grown/'V:adam_best:deployed_a1_T03'/rows[P_matched]", status="DERIVED", note="L35: +0.035")
    put("A4_L2_A1_T03", _slope(matched("L2:adam_best:deployed_a1_T03")), "s26/results/q_var.json :: results/grown/'L2:adam_best:deployed_a1_T03'/rows[P_matched]", status="DERIVED", note="L35: -0.008")
    put("A4_V_A025_T03", _slope(matched("V:adam_best:deployed_a025_T03")), "s26/results/q_var.json :: results/grown/'V:adam_best:deployed_a025_T03'/rows[P_matched]", status="DERIVED", note="L35: -0.246")
    put("A4_L2_A025_T03", _slope(matched("L2:adam_best:deployed_a025_T03")), "s26/results/q_var.json :: results/grown/'L2:adam_best:deployed_a025_T03'/rows[P_matched]", status="DERIVED", note="L35: -0.302")
    put("A4_FIXED_A025_T03", qvar["fixed"]["deployed_a025_T03"]["log2_slope_per_qubit"], "s26/results/q_var.json :: results/fixed/deployed_a025_T03/log2_slope_per_qubit")
    nd = [r["n_distinct_ops"] for key in ("V:adam_best:linear_alpha1_T0", "L2:adam_best:linear_alpha1_T0", "V:adam_best:deployed_a1_T03", "L2:adam_best:deployed_a1_T03") for r in matched(key)]
    put("A4_A1_DISTINCT_MIN", int(min(nd)), "s26/results/q_var.json :: results/grown/[alpha=1 cells]/rows[P_matched]/n_distinct_ops (min)", status="DERIVED")
    put("A4_A1_DISTINCT_MAX", int(max(nd)), "s26/results/q_var.json :: results/grown/[alpha=1 cells]/rows[P_matched]/n_distinct_ops (max)", status="DERIVED")
    nd2 = [r["n_distinct_ops"] for r in matched("L2:adam_best:deployed_a025_T03")]
    put("A4_L2_A025_DISTINCT_MIN", int(min(nd2)), "s26/results/q_var.json :: results/grown/'L2:adam_best:deployed_a025_T03'/rows[P_matched]/n_distinct_ops (min)", status="DERIVED")
    put("A4_L2_A025_DISTINCT_MAX", int(max(nd2)), "s26/results/q_var.json :: results/grown/'L2:adam_best:deployed_a025_T03'/rows[P_matched]/n_distinct_ops (max)", status="DERIVED")
    fixed_by_cell = {cell: {r["n"]: r["var_g0"] for r in qvar["fixed"][cell]["rows"]} for cell in qvar["fixed"]}
    ratios = []
    for key in grown:
        cell = key.split(":")[-1]
        for r in matched(key):
            ratios.append(r["var_g0"] / fixed_by_cell[cell][r["n"]])
    put("A4_RATIO_MIN", float(min(ratios)), "s26/results/q_var.json :: grown var_g0 / fixed var_g0 at matched (cell, n), min over the matched rows", status="DERIVED", note="L35: 1.63")
    put("A4_RATIO_MAX", float(max(ratios)), "s26/results/q_var.json :: the same, max", status="DERIVED", note="L35: 370")
    put("A4_N_MATCHED", len(ratios), "s26/results/q_var.json :: number of matched rows", status="DERIVED")
    put("A4_REPRO_WORST_REL", qvar["reproduction"]["worst_rel"], "s26/results/q_var.json :: results/reproduction/worst_rel", note="S25 n = 7 rows reproduced")

    probe = _j("s26/results/probe/1A13.json")
    put("PROD_KL_IDEAL", probe["product_diagnostics"]["ideal_ladder"]["kl_gibbs_to_product"], "s26/results/probe/1A13.json :: product_diagnostics/ideal_ladder/kl_gibbs_to_product")
    put("PROD_KL_1A13", probe["product_diagnostics"]["zrank"]["kl_gibbs_to_product"], "s26/results/probe/1A13.json :: product_diagnostics/zrank/kl_gibbs_to_product")
    put("PROD_KL_RY7", probe["arms"]["adaptV_lbfgs_zrank_P7"]["kl_to_gibbs"], "s26/results/probe/1A13.json :: arms/adaptV_lbfgs_zrank_P7/kl_to_gibbs", note="a 7-parameter RY layer, L-BFGS-B")
    put("PROD_KL_FIXED", probe["arms"]["fixed_zrank_it50"]["kl_to_gibbs"], "s26/results/probe/1A13.json :: arms/fixed_zrank_it50/kl_to_gibbs", note="the deployed 21-parameter circuit, 50 Adam steps")

    # ------------------------------------------------------------------ B1 (L13) and C3 (L24, L39)
    b1 = _j("s26/results/b1_feasibility.json")
    ram = b1["measured_2026-09-13"]["c_size_and_ram"]
    put("B1_RESIDENT_GB", ram["ram_as_fair_esm_loads_it_GB"]["total_resident_minimum"], "s26/results/b1_feasibility.json :: measured_2026-09-13/c_size_and_ram/ram_as_fair_esm_loads_it_GB/total_resident_minimum")
    put("B1_FP16_GB", ram["ram_fp16_GB_parameters_only"]["total"], "s26/results/b1_feasibility.json :: .../ram_fp16_GB_parameters_only/total")
    put("B1_DOWNLOAD_GB", ram["download_bytes_http_head"]["total_GB"], "s26/results/b1_feasibility.json :: .../download_bytes_http_head/total_GB")
    put("B1_VERDICT", b1["verdict"].split(".")[0], "s26/results/b1_feasibility.json :: verdict")
    put("B1_HEADROOM_GB", 4.4, "s26/LEDGER.md L13 (governor_state 00:19: 64.8% of 16.75 GB used, ceiling 93%)", status="LEDGER", note="campaign headroom at the time of the B1 measurement")

    c3 = _j("s26/results/ph_c3_stage1.json")["summary"]
    def c3put(tok, key, field, **kw):
        put(tok, c3[key][field], f"s26/results/ph_c3_stage1.json :: summary/'{key}'/{field}", **kw)
    c3put("C3_AMBER_VS_NONE", "stage1 AMBER (relaxed) minus do-nothing (built chain)", "effect", basis="relaxed minus built_chain")
    c3put("C3_AMBER_VS_NONE_CI", "stage1 AMBER (relaxed) minus do-nothing (built chain)", "ci95_fold")
    c3put("C3_AMBER_VS_NONE_X", "stage1 AMBER (relaxed) minus do-nothing (built chain)", "effect_over_mde")
    c3put("C3_AMBER_VS_RANDOM", "stage1 AMBER minus matched-magnitude RANDOM (16 draws)", "effect", basis="relaxed minus displaced built_chain")
    c3put("C3_AMBER_VS_RANDOM_CI", "stage1 AMBER minus matched-magnitude RANDOM (16 draws)", "ci95_fold")
    c3put("C3_AMBER_VS_RANDOM_X", "stage1 AMBER minus matched-magnitude RANDOM (16 draws)", "effect_over_mde")
    c3put("C3_AMBER_VS_RANDOM_MDE", "stage1 AMBER minus matched-magnitude RANDOM (16 draws)", "mde")
    c3put("C3_AMBER_VS_RANDOM_FOLDS", "stage1 AMBER minus matched-magnitude RANDOM (16 draws)", "folds_same_sign")
    c3put("C3_AMBER_VS_MEMBER", "stage1 AMBER minus matched-magnitude TOWARD-MEMBER (16 draws)", "effect", basis="relaxed minus displaced built_chain")
    c3put("C3_AMBER_VS_MEMBER_CI", "stage1 AMBER minus matched-magnitude TOWARD-MEMBER (16 draws)", "ci95_fold")
    c3put("C3_RANDOM_VS_NONE", "stage1 RANDOM minus do-nothing", "effect", basis="displaced built_chain minus built_chain")
    c3put("C3_RANDOM_VS_NONE_CI", "stage1 RANDOM minus do-nothing", "ci95_fold")
    c3put("C3_MEMBER_VS_NONE", "stage1 TOWARD-MEMBER minus do-nothing", "effect", basis="displaced built_chain minus built_chain", note="provisional until the registered replication lands (L39)")
    put("C3_MAG", _j("s26/results/ph_c3_nativefree.json")["summary"]["mag_sup"]["mean"], "s26/results/ph_c3_nativefree.json :: summary/mag_sup/mean", note="per-atom RMS displacement of the CA trace by the relaxation, after superposition")
    rows3 = _j("s26/results/ph_c3_stage1.json")["rows"]
    put("C3_COS", float(np.mean([r["cos_amber"] for r in rows3])), "s26/results/ph_c3_stage1.json :: rows[*]/cos_amber (mean)", status="DERIVED", note="ORACLE cosine between the relaxation displacement and the true residual")
    put("C3_COS_FRAC_POS", c3["frac_cos_amber_positive"], "s26/results/ph_c3_stage1.json :: summary/frac_cos_amber_positive")
    nf = _j("s26/results/ph_c3_nativefree.json")["summary"]
    put("C3_FRAC_E0_OVER_1E4", nf["frac_e0_over_1e4"], "s26/results/ph_c3_nativefree.json :: summary/frac_e0_over_1e4")
    put("C3_N_CONVERGED", int(nf["n_converged_le_1000"]), "s26/results/ph_c3_nativefree.json :: summary/n_converged_le_1000")
    put("C3_E1_MAX", nf["e1_max"], "s26/results/ph_c3_nativefree.json :: summary/e1_max", note="9KAR")
    put("C3_BOND_AFTER", nf["amb_bond_mean"]["mean"], "s26/results/ph_c3_nativefree.json :: summary/amb_bond_mean/mean")
    put("C3_BOND_BEFORE", nf["ca_bond_mean"]["mean"], "s26/results/ph_c3_nativefree.json :: summary/ca_bond_mean/mean")

    # ------------------------------------------------------------------ L43 (steric reject), L44/L58 (the 2/60 bound), L53 (strain)
    rr = _j("s26/results/ph_reject_report.json")["report"]["point_cloud"]["1e4"]
    put("REJECT_R_VS_ANCHOR", rr["R_vs_anchor_all"]["effect"], "s26/results/ph_reject_report.json :: report/point_cloud/1e4/R_vs_anchor_all/effect",
        basis="point_cloud", note="L43; Type-M zone (1.17x MDE), tail-carried (L54)")
    put("REJECT_R_VS_ANCHOR_CI", rr["R_vs_anchor_all"]["ci95_fold"], "s26/results/ph_reject_report.json :: report/point_cloud/1e4/R_vs_anchor_all/ci95_fold", basis="point_cloud")
    put("REJECT_R_VS_ANCHOR_MEDIAN", rr["R_vs_anchor_all"]["median_effect"], "s26/results/ph_reject_report.json :: report/point_cloud/1e4/R_vs_anchor_all/median_effect", basis="point_cloud")
    put("REJECT_R_VS_RANDR", rr["R_vs_RANDR_all"]["effect"], "s26/results/ph_reject_report.json :: report/point_cloud/1e4/R_vs_RANDR_all/effect",
        basis="point_cloud", note="L43; Type-M zone (1.26x MDE)")
    wb = _j("s26/results/w_selfcopy_bound.json")
    er = wb["verdict"]["arm"]["envelope_readings_A"]
    put("LEAK_BOUND_ARM_MEANCI", er["mean_ci_limit"], "s26/results/w_selfcopy_bound.json :: verdict/arm/envelope_readings_A/mean_ci_limit", basis="built_chain", note="L44/L58; (2/60) x the own-native envelope's fold-CI limit")
    put("LEAK_BOUND_ARM_WORST", er["worst_target"], "s26/results/w_selfcopy_bound.json :: verdict/arm/envelope_readings_A/worst_target", basis="built_chain", note="L58; worst single target " + str(er["worst_target_pdb"]))
    put("LEAK_CLASS", wb["verdict"]["arm"]["class_under_every_reading"], "s26/results/w_selfcopy_bound.json :: verdict/arm/class_under_every_reading")
    put("LEAK_BENCH_CI_HALF", wb["bench_ci_half"], "s26/results/w_selfcopy_bound.json :: bench_ci_half", note="half-width of the benchmark CI, the materiality scale")
    we = _j("s26/results/w_selfcopy_endpoint.json")["A"]
    put("LEAK_DEV_PRICE_FIT", -we["ge06_mean_delta_all126"]["fit"], "s26/results/w_selfcopy_endpoint.json :: -A/ge06_mean_delta_all126/fit (clean minus production, lam = 0 chain)",
        status="DERIVED", basis="fit (lam = 0 chain)", note="L44: S10-4's +0.0004 re-derived; fold CI [-0.0001, +0.0010], MDE 0.0012")
    put("LEAK_DEV_PRICE_ARM", -we["ge06_mean_delta_all126"]["arm"], "s26/results/w_selfcopy_endpoint.json :: -A/ge06_mean_delta_all126/arm (clean minus production)",
        status="DERIVED", basis="built_chain", note="L44: +0.0018 [-0.0003, +0.0039]")
    put("LEAK_DEV4_BOUND_ARM", wb["signed_bounds_gated"]["A_real4"]["arm"],
        "s26/results/w_selfcopy_bound.json :: signed_bounds_gated/A_real4/arm (dev-proxy price, (2/60) x the signed max over the four dev self-copies)",
        basis="built_chain", note="L44 Part D: 0.0023, IMMATERIAL")
    stt = _j("s26/results/ph_strain.json")
    sm = stt["summary"]["partial_n_rg"]["moved"]
    put("STRAIN_RHO_MOVED", sm["rho"], "s26/results/ph_strain.json :: summary/partial_n_rg/moved/rho", note="L53; Spearman with rmsd_arm, partial on n and Rg; replicated (ph_strain_rep.json)")
    put("STRAIN_RHO_MOVED_CI", sm["ci95_fold"], "s26/results/ph_strain.json :: summary/partial_n_rg/moved/ci95_fold")
    srows = stt["rows"]
    mv = np.array([r["moved"] for r in srows]); ar = np.array([r["rmsd_arm"] for r in srows])
    order = np.argsort(mv, kind="stable")
    qm = [float(ar[order[0:32]].mean()), float(ar[order[32:64]].mean()), float(ar[order[64:96]].mean()), float(ar[order[96:]].mean())]
    put("STRAIN_QUARTILE_MEANS", qm, "s26/results/ph_strain.json :: rows[*]/{moved, rmsd_arm}: mean rmsd_arm by quartile of moved (32/32/32/30)",
        status="DERIVED", basis="built_chain", note="L53: 2.286 / 2.936 / 3.758 / 3.923; a calibration flag, never a gain")

    # ------------------------------------------------------------------ Proposal C: C1 closures (L26), C2 anchor (L56/L57)
    put("C1_S12_REAL_N8", _j("s12/results/agg_dec_v2_n8.json")["mean_raw"], "s12/results/agg_dec_v2_n8.json :: mean_raw", basis="point_cloud (weighted average)",
        note="L26: set-transformer ranker, real features, 8 training targets")
    put("C1_S12_REAL_FULL", _j("s12/results/agg_dec_v2.json")["mean_raw"], "s12/results/agg_dec_v2.json :: mean_raw", basis="point_cloud (weighted average)",
        note="L26: the same, full training set (ntrain None = all)")
    put("C1_S12_LEAK_N8", _j("s12/results/agg_dec_v2_oracle_n8.json")["mean_raw"], "s12/results/agg_dec_v2_oracle_n8.json :: mean_raw", basis="point_cloud (weighted average)",
        note="L26: the same model with a leaked (ORACLE) label, 8 training targets")
    put("C1_S12_LEAK_FULL", _j("s12/results/agg_dec_v2_oracle.json")["mean_raw"], "s12/results/agg_dec_v2_oracle.json :: mean_raw", basis="point_cloud (weighted average)",
        note="ORACLE label, full training set")
    ib = _j("s17/results/inband.json")["rows"]
    put("C1_S17_BAND_BEST", float(np.mean([r["cells"]["25"]["band_best"] for r in ib])), "s17/results/inband.json :: rows[*]/cells/25/band_best (mean over 126)",
        status="DERIVED", basis="single_window", note="L26: a perfect ranker inside the shipped top-25 (ORACLE)")
    put("C1_S17_RAND", float(np.mean([r["cells"]["25"]["rand"] for r in ib])), "s17/results/inband.json :: rows[*]/cells/25/rand (mean over 126)",
        status="DERIVED", basis="single_window", note="a random member of the top-25")
    put("C1_S17_DIST_PICK", float(np.mean([r["cells"]["25"]["dist_pick"] for r in ib])), "s17/results/inband.json :: rows[*]/cells/25/dist_pick (mean over 126)",
        status="DERIVED", basis="single_window", note="the distogram argmin")
    pa = _j("s26/results/p_ladder_report_shipped_s0.json")
    put("C2_ANCHOR_SEL_MAXABS", pa["sel"]["maxabs_vs_cache"], "s26/results/p_ladder_report_shipped_s0.json :: sel/maxabs_vs_cache", note="L56: the anchor vs the production cache, selection")
    put("C2_ANCHOR_CLOUD_MAXABS", pa["cloud"]["maxabs_vs_cache"], "s26/results/p_ladder_report_shipped_s0.json :: cloud/maxabs_vs_cache", note="L56: point cloud")
    put("C2_ANCHOR_ARM", pa["arm"]["mean"], "s26/results/p_ladder_report_shipped_s0.json :: arm/mean", basis="built_chain (leaderboard-rebuild basis, L57)")
    put("C2_ANCHOR_ARM_MAXABS", pa["arm"]["maxabs_vs_cache"], "s26/results/p_ladder_report_shipped_s0.json :: arm/maxabs_vs_cache", note="L57: the projection lands in a different local optimum on 120/126")
    mdir = os.path.join(ROOT, "s26", "models", "p_ladder")
    present = set(os.listdir(mdir)) if os.path.isdir(mdir) else set()
    cands = {f.split("_fold")[0] for f in present if f.endswith("_s0.pt") and "probe" not in f}
    rungs_trained = sorted(r for r in cands if all(f"{r}_fold{k}_s0.pt" in present for k in range(5)))
    put("C2_RUNGS_TRAINED", len(rungs_trained), "s26/models/p_ladder/<rung>_fold{0..4}_s0.pt (rungs with all five fold checkpoints; L41 lists noesm, conly, pca32, wide, pca32f, pca128, esm8m)",
        status="DERIVED", note=", ".join(rungs_trained))
    put("C2_N_RUNGS", 11, "s26/PREREG_C2.md / s26/PROPOSAL_C.md: rungs shipped, noesm, conly, pca32, pca32f, pca128, raw, esm8m, wide, pairnet, mix", status="LEDGER")

    # ------------------------------------------------------------------ tests, corpus, leak counts
    tr = _j("s26/results/test_run.json")["combined"]
    put("TESTS_TOTAL", int(tr["total"]), "s26/results/test_run.json :: combined/total")
    put("TESTS_PASSED", int(tr["passed"]), "s26/results/test_run.json :: combined/passed")
    put("TESTS_SKIPPED", int(tr["skipped"]), "s26/results/test_run.json :: combined/skipped")
    put("TESTS_FAILED", int(tr["failed"]), "s26/results/test_run.json :: combined/failed")
    pe = _j("s26/results/p_probe_esm.json")
    put("CORPUS_PEPTIDES", 787, "s26/results/p_probe_esm.json (ledger L11: 787 peptides + 6,003 fragments = 6,790 unique training sequences)", status="LEDGER")
    put("CORPUS_FRAGMENTS", 6003, "same", status="LEDGER")
    ia = _j("s26/results/i_identity_audit.json")
    put("SELFCOPY_DEV", "4/126", "s26/results/i_identity_audit.json (ledger L15); tests/test_data.py:72-74", status="SOURCED", note="claim C28")
    put("SELFCOPY_BENCH", "2/60", "s26/results/i_identity_audit.json (count only, ledger L15/L18)", status="SOURCED", note="claim C28, count only")
    put("MDE_FACTOR", 2.8016, "s26/LANE_CONTRACT.md section 2; s25/results/q_verify.json :: results/'MDE == 2.8016*SE exactly' (True)", status="SOURCED")
    put("SCIENCE_DELTA_MAX", max(abs(v) for v in _j("bench_results/compare_tuning126.json")["science_delta"].values()),
        "bench_results/compare_tuning126.json :: science_delta (max |value|)", status="DERIVED", note="1-worker baseline vs 8-worker optimised run")
    return V


def fmt(V, token, spec=None):
    """Format a registry value for a slide: numbers by spec, strings and lists verbatim."""
    v = V[token]["value"]
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, (int, np.integer)) and spec is None:
        return f"{int(v):,}"
    if isinstance(v, (float, np.floating)):
        return format(float(v), spec or ".4f")
    if isinstance(v, list):
        if all(isinstance(x, (int, np.integer)) and not isinstance(x, bool) for x in v):
            return ", ".join(str(int(x)) for x in v)
        return "[" + ", ".join(format(float(x), spec or "+.4f") for x in v) + "]"
    return str(v)


if __name__ == "__main__":
    V = load_values()
    for k, d in V.items():
        print(f"{k:28s} {str(d['value'])[:40]:40s} {d['status']:16s} {d['basis']:34s} {d['path']}")
