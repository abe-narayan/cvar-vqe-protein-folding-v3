#!/usr/bin/env python
"""s31/s31_F_analyse.py -- S31 LANE F analysis: F1 (the terminal operator + the gate) and
F2 (what distinguishes the tail).  Consumes `s31/results/s31_F_terminal_rows.jsonl` only.

Pre-registered in `s31/PREREG_S31_F.md`.  Endpoint basis is stated on every number.

    python s31/s31_F_analyse.py
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s24 import stats_lib as ST            # noqa: E402
from s12 import instrument as I            # noqa: E402

RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "s31_F_terminal_rows.jsonl")
OUT = os.path.join(RESULTS, "s31_F_analyse.json")
SEED = 31006
NCOMP = [0]


def load():
    R = [json.loads(ln) for ln in open(ROWS) if ln.strip()]
    R.sort(key=lambda r: r["pdb"])
    if len(R) != 126:
        raise RuntimeError("expected 126 rows, have %d -- the run is not finished" % len(R))
    return R


def col(R, k):
    return np.array([r[k] for r in R], float)


def cmp2(a, b, folds, names, label):
    NCOMP[0] += 1
    o = ST.compare(a, b, folds=folds, names=names, label=label, seed_parts=("s31F", str(SEED)))
    return o


def brief(o):
    return dict(label=o["label"], mean_a=o["mean_a"], mean_b=o["mean_b"], effect=o["effect"],
                median_effect=o["median_effect"], se=o["se"], mde=o["mde"],
                effect_over_mde=o["effect_over_mde"], ci95_fold=o["ci95_fold"],
                folds_same_sign=o["folds_same_sign"], n_folds=o["n_folds"],
                W=o["n_better"], L=o["n_worse"], T=o["n_tied"], n=o["n"])


def verdict(o):
    """The project's fixed rule, quoted from the gate and never from `.verdict` (rule 8)."""
    r = abs(o["effect_over_mde"])
    ci = o["ci95_fold"]
    ex = bool(ci is not None and (ci[0] > 0 or ci[1] < 0))
    if r >= 1.0 and ex:
        return "MEASURED (|effect| >= 1.0x MDE, fold CI excludes zero)"
    if r >= 1.0 and not ex:
        return "NOT MEASURED (>=1.0x MDE but fold CI spans zero)"
    if r >= 0.7:
        return "NOT MEASURED (0.7-1.0x MDE)"
    return "NOT A RESULT (< 0.7x MDE); may support a NULL, never a presence"


def spearman(x, y):
    from scipy.stats import rankdata
    x = np.asarray(x, float); y = np.asarray(y, float)
    a, b = rankdata(x), rankdata(y)
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 0 else float("nan")


def two_group(v, folds, mask, label, n_boot=4000):
    """Unpaired hi-minus-lo contrast with a fold-clustered bootstrap (s30_G_disp's construction)."""
    NCOMP[0] += 1
    v = np.asarray(v, float); folds = np.asarray(folds); mask = np.asarray(mask, bool)
    a, b = v[mask], v[~mask]
    fa, fb = folds[mask], folds[~mask]
    rng = np.random.default_rng(SEED)
    eff = float(a.mean() - b.mean())
    uf = sorted(set(folds.tolist()))
    draws = []
    for _ in range(n_boot):
        pick = rng.choice(uf, size=len(uf), replace=True)
        va = np.concatenate([a[fa == f] for f in pick])
        vb = np.concatenate([b[fb == f] for f in pick])
        if va.size < 2 or vb.size < 2:
            continue
        draws.append(va.mean() - vb.mean())
    draws = np.asarray(draws, float)
    lo, hi = float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))
    se_iid = float(np.sqrt(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size))
    se = max(se_iid, float(draws.std(ddof=1)))
    mde = 2.8016 * se
    per = []
    for f in uf:
        aa, bb = a[fa == f], b[fb == f]
        per.append(float(aa.mean() - bb.mean()) if aa.size and bb.size else float("nan"))
    ss = [s for s in per if np.isfinite(s)]
    return dict(label=label, effect=eff, mean_hi=float(a.mean()), mean_lo=float(b.mean()),
                n_hi=int(a.size), n_lo=int(b.size), se=se, mde=float(mde),
                effect_over_mde=float(eff / mde) if mde else float("nan"),
                ci95_fold=[lo, hi], excludes_zero=bool(lo > 0 or hi < 0),
                folds_same_sign="%d/%d" % (sum(1 for s in ss if np.sign(s) == np.sign(eff)), len(ss)),
                per_fold=per)


def main():
    R = load()
    names = [r["pdb"] for r in R]
    folds = np.array([r["fold"] for r in R])
    out = dict(seed=SEED, n=len(R), rows=ROWS,
               prereg="s31/PREREG_S31_F.md @ 8a14edea + 49ee7c92")

    cl_A, cl_M, cl_R = col(R, "cloud_AVG"), col(R, "cloud_MED"), col(R, "cloud_AVG_RG")
    cl_S = col(R, "cloud_AVG_SEP")
    ch_A, ch_M, ch_R = col(R, "chain_AVG"), col(R, "chain_MED"), col(R, "chain_AVG_RG")
    ch_S = col(R, "chain_AVG_SEP")
    P_A, P_M, P_R = col(R, "P_AVG"), col(R, "P_MED"), col(R, "P_AVG_RG")
    P_S = col(R, "P_AVG_SEP")
    MOVE_A, MOVE_M = col(R, "move_AVG"), col(R, "move_MED")
    DISP, S, s_b = col(R, "DISP"), col(R, "S"), col(R, "s_b")

    # ---------------------------------------------------------------- reproduction
    rec_avg, rec_fit = col(R, "rec_avg"), col(R, "rec_fit")
    out["reproduction"] = dict(
        basis="CA point cloud / built chain, n=126",
        cloud_AVG_mean=float(cl_A.mean()), record_rmsd_avg_mean=float(rec_avg.mean()),
        cloud_max_abs_dev=float(np.abs(cl_A - rec_avg).max()),
        chain_AVG_mean=float(ch_A.mean()), record_rmsd_fit_mean=float(rec_fit.mean()),
        chain_max_abs_dev=float(np.abs(ch_A - rec_fit).max()),
        note="chain_AVG is production RE-PROJECTED IN THIS PROCESS and is the comparator for "
             "every arm below (contract rule 31: the projection seed is not pinned, and the "
             "cross-instrument spread is 0.0107 A).")

    # ---------------------------------------------------------------- SS2: the derivation
    pred = np.sqrt(cl_A ** 2 + s_b ** 2)
    out["S2_theorem"] = dict(
        claim="d_MED_cloud^2 ~ B^2 + s_b^2 with B = d_AVG_cloud -- derived in the prereg before "
              "any number existed",
        basis="CA point cloud",
        pred_mean=float(pred.mean()), obs_mean=float(cl_M.mean()),
        mean_signed_err=float((cl_M - pred).mean()),
        mean_abs_err=float(np.abs(cl_M - pred).mean()),
        rho_pred_obs=spearman(pred, cl_M),
        frac_MED_worse_than_AVG_on_cloud=float((cl_M > cl_A).mean()),
        rho_s_b_vs_cloud_gap=spearman(s_b, cl_M - cl_A),
        rho_DISP_vs_cloud_gap=spearman(DISP, cl_M - cl_A),
        rho_DISP_vs_P_AVG=spearman(DISP, P_A),
        note="the derivation predicted BOTH sides of the inequality grow with the spread")

    # the registered pre-check (rule 22): is the cloud gap already bigger than the budget?
    g_cloud = cmp2(cl_M, cl_A, folds, names, "F.precheck.cloud MED-AVG (CA POINT CLOUD)")
    out["S2_precheck"] = dict(cmp=brief(g_cloud), verdict=verdict(g_cloud),
                              budget_P_AVG_mean=float(P_A.mean()),
                              dead_on_arrival=bool(g_cloud["effect"] > P_A.mean()))

    # ---------------------------------------------------------------- F1: the operators
    ops = {}
    for nm, v in (("MED", ch_M), ("AVG_RG", ch_R), ("AVG_SEP", ch_S)):
        o = cmp2(v, ch_A, folds, names, "F1.%s - AVG (BUILT CHAIN)" % nm)
        ops[nm] = dict(cmp=brief(o), verdict=verdict(o))
        oc = cmp2({"MED": cl_M, "AVG_RG": cl_R, "AVG_SEP": cl_S}[nm], cl_A, folds, names,
                  "F1.%s - AVG (CA POINT CLOUD)" % nm)
        ops[nm]["cloud_cmp"] = brief(oc)
    ops["MED"]["record_cross_check"] = dict(
        source="s12/agg_FINDINGS.md:43-70 (predates this lane; found only AFTER the registration "
               "was committed -- declared as lane F's own defect in PREREG s10.4)",
        record_medoid75_cloud=3.2822, record_avg75_cloud=3.0483, record_delta=0.2339,
        record_ci=[0.162, 0.305], record_WL="34/92",
        this_run_cloud_MED=float(cl_M.mean()), this_run_cloud_AVG=float(cl_A.mean()),
        this_run_delta=float((cl_M - cl_A).mean()),
        prereg_predicted_chain_delta=0.0717)
    ops["AVG_RG"]["record_cross_check"] = dict(
        source="s23/LEDGER.md:143-176 via s31/PREREG_S31_C.md A1.1 -- the native-free `pool` scale "
               "arm, CLOSED at +0.095 with the CI excluding zero on the bad side",
        note="AVG_RG is a REGISTERED NEGATIVE CONTROL here, not a candidate")
    out["F1_operators"] = dict(basis="BUILT CHAIN, n=126, paired, fold-clustered CI on pinned folds",
                               chain_AVG_mean=float(ch_A.mean()),
                               chain_MED_mean=float(ch_M.mean()),
                               chain_AVG_RG_mean=float(ch_R.mean()),
                               chain_AVG_SEP_mean=float(ch_S.mean()), arms=ops)

    # ---------------------------------------------------------------- F1: the gates
    med_disp = float(np.median(DISP))
    hi = DISP > med_disp
    gates = {}

    # G0: median split, registered direction high-DISP -> MED
    g0 = np.where(hi, ch_M, ch_A)
    o = cmp2(g0, ch_A, folds, names, "F1.G0 median-DISP gate (hi->MED) - AVG (BUILT CHAIN)")
    gates["G0"] = dict(rule="DISP > median(DISP) -> MED else AVG; 0 free parameters; NATIVE-FREE",
                       threshold=med_disp, n_to_MED=int(hi.sum()),
                       cmp=brief(o), verdict=verdict(o))
    # G0 reversed -- reported because the registered direction can fire against me
    g0r = np.where(hi, ch_A, ch_M)
    o = cmp2(g0r, ch_A, folds, names, "F1.G0rev median-DISP gate (LO->MED) - AVG (BUILT CHAIN)")
    gates["G0_reversed"] = dict(rule="DISP <= median(DISP) -> MED else AVG; NATIVE-FREE; this is "
                                     "the OPPOSITE of the registered direction and is reported "
                                     "only because the registered one may fire against the lane",
                                cmp=brief(o), verdict=verdict(o))

    # G1: leave-fold-out threshold on native chain RMSD -- ORACLE-ADJACENT
    grid = np.quantile(DISP, np.linspace(0.05, 0.95, 19))
    g1 = np.empty(len(R))
    tau_by_fold = {}
    for f in sorted(set(folds.tolist())):
        tr = folds != f
        best, bt = None, None
        for tau in grid:
            v = np.where(DISP[tr] > tau, ch_M[tr], ch_A[tr]).mean()
            if best is None or v < best:
                best, bt = v, float(tau)
        tau_by_fold[int(f)] = bt
        g1[folds == f] = np.where(DISP[folds == f] > bt, ch_M[folds == f], ch_A[folds == f])
    o = cmp2(g1, ch_A, folds, names, "F1.G1 leave-fold-out DISP threshold - AVG (BUILT CHAIN)")
    gates["G1"] = dict(rule="DISP > tau, tau fitted LEAVE-FOLD-OUT on native chain RMSD",
                       ORACLE="ORACLE-ADJACENT / the threshold touches the native -- a diagnostic, "
                              "NOT a clean deployable",
                       tau_by_fold=tau_by_fold, cmp=brief(o), verdict=verdict(o))

    # G2: per-target min -- ORACLE, with a split-half transfer arm
    g2 = np.minimum(ch_M, ch_A)
    o = cmp2(g2, ch_A, folds, names, "F1.G2 per-target min(AVG,MED) - AVG (BUILT CHAIN)")
    Mmat = np.stack([ch_A, ch_M], 1)
    sh = ST.split_half_transfer(Mmat, seed_parts=("s31F", str(SEED)))
    gates["G2"] = dict(rule="per-target min(AVG, MED)",
                       ORACLE="ORACLE / NOT DEPLOYABLE -- this is best-of-2, not skill",
                       cmp=brief(o), verdict=verdict(o), split_half_transfer=sh)

    # G3: the MOVE gate -- EXPLORATORY (registered mid-run, prereg s9)
    med_move = float(np.median(MOVE_A))
    him = MOVE_A > med_move
    g3 = np.where(him, ch_M, ch_A)
    o = cmp2(g3, ch_A, folds, names, "F1.G3 median-MOVE gate (hi->MED) - AVG (BUILT CHAIN)")
    Mmat3 = np.stack([ch_A, ch_M], 1)
    gates["G3"] = dict(rule="MOVE(AVG) > median(MOVE(AVG)) -> MED else AVG; 0 free parameters; "
                            "NATIVE-FREE (both arguments of MOVE exist at inference)",
                       EXPLORATORY="registered mid-run (PREREG s9) -- if this is the only arm that "
                                   "clears it is exploratory and requires confirmation, never the "
                                   "lane's confirmed result",
                       threshold=med_move, n_to_MED=int(him.sum()),
                       cmp=brief(o), verdict=verdict(o),
                       split_half_transfer=ST.split_half_transfer(
                           Mmat3, seed_parts=("s31F", str(SEED), "g3")))
    g3r = np.where(him, ch_A, ch_M)
    o = cmp2(g3r, ch_A, folds, names, "F1.G3rev median-MOVE gate (LO->MED) - AVG (BUILT CHAIN)")
    gates["G3_reversed"] = dict(rule="the OPPOSITE direction, reported because the registered one "
                                     "may fire against the lane", cmp=brief(o), verdict=verdict(o))
    # G3tau: leave-fold-out MOVE threshold -- ORACLE-ADJACENT
    gridm = np.quantile(MOVE_A, np.linspace(0.05, 0.95, 19))
    g3t = np.empty(len(R)); tau_m = {}
    for f in sorted(set(folds.tolist())):
        tr = folds != f
        best, bt = None, None
        for tau in gridm:
            v = np.where(MOVE_A[tr] > tau, ch_M[tr], ch_A[tr]).mean()
            if best is None or v < best:
                best, bt = v, float(tau)
        tau_m[int(f)] = bt
        g3t[folds == f] = np.where(MOVE_A[folds == f] > bt, ch_M[folds == f], ch_A[folds == f])
    o = cmp2(g3t, ch_A, folds, names, "F1.G3tau leave-fold-out MOVE threshold - AVG (BUILT CHAIN)")
    gates["G3tau"] = dict(rule="MOVE(AVG) > tau, tau fitted LEAVE-FOLD-OUT on native chain RMSD",
                          ORACLE="ORACLE-ADJACENT / the threshold touches the native -- a "
                                 "diagnostic, NOT a clean deployable",
                          tau_by_fold=tau_m, cmp=brief(o), verdict=verdict(o))
    out["F1_gates"] = gates
    out["F1_gate_variables"] = dict(
        rho_DISP_vs_chain_MED_minus_AVG=spearman(DISP, ch_M - ch_A),
        rho_MOVE_vs_chain_MED_minus_AVG=spearman(MOVE_A, ch_M - ch_A),
        rho_DISP_vs_MOVE=spearman(DISP, MOVE_A),
        rho_MOVE_vs_P_AVG=spearman(MOVE_A, P_A),
        MOVE_AVG_mean=float(MOVE_A.mean()), MOVE_MED_mean=float(MOVE_M.mean()),
        note="MOVE(MED) ~ 0 is the check that the medoid is already a chain; a nonzero MOVE(MED) "
             "would mean the pool members are not ideal-geometry backbones")

    # ---------------------------------------------------------------- F1: the mechanism contrast
    d_MA = ch_M - ch_A
    out["F1_mechanism"] = dict(
        basis="BUILT CHAIN; hi-minus-lo contrast of chain(MED) - chain(AVG) across the DISP median",
        registered_expectation="NEGATIVE (the medoid should help most where the pool diverges)",
        contrast=two_group(d_MA, folds, hi, "F1.mech hi-minus-lo DISP of chain(MED)-chain(AVG)"),
        rho_DISP_vs_dMA=spearman(DISP, d_MA),
        mean_dMA_hi=float(d_MA[hi].mean()), mean_dMA_lo=float(d_MA[~hi].mean()),
        by_DISP_tertile=[dict(t=int(q),
                              disp_lo=float(np.quantile(DISP, q / 3.0)),
                              disp_hi=float(np.quantile(DISP, (q + 1) / 3.0)),
                              mean_dMA=float(d_MA[(DISP >= np.quantile(DISP, q / 3.0)) &
                                                  (DISP <= np.quantile(DISP, (q + 1) / 3.0))].mean()),
                              n=int(((DISP >= np.quantile(DISP, q / 3.0)) &
                                     (DISP <= np.quantile(DISP, (q + 1) / 3.0))).sum()))
                         for q in range(3)])

    # ---------------------------------------------------------------- AVG_SEP: mechanism + risks
    d_SA = ch_S - ch_A
    def worst18_(v):
        m = np.zeros(len(R), bool); m[np.argsort(-v)[:18]] = True; return m
    pool_mean0, pool_best0 = col(R, "pool_mean"), col(R, "pool_best")
    tails0 = {"T_POOL_worst18_by_pool_mean": worst18_(pool_mean0),
              "T_BEST_worst18_by_pool_best": worst18_(pool_best0),
              "T_CHAIN_worst18_by_production_chain": worst18_(ch_A),
              "FAIL18_DIAGNOSTIC_ONLY": np.array([r["fail18"] for r in R], bool)}
    out["AVG_SEP_mechanism"] = dict(
        basis="BUILT CHAIN, paired, n=126",
        registered_prediction="if the separation-band story is right the gain must CONCENTRATE on "
                              "high-dispersion targets, because averaging's distortion scales with "
                              "spread (rho(DISP, contraction) = +0.947). A FLAT profile in "
                              "dispersion refutes the mechanism even if the number is good.",
        full_paired_distribution=dict(
            mean=float(d_SA.mean()), median=float(np.median(d_SA)), sd=float(d_SA.std(ddof=1)),
            W=int((d_SA < 0).sum()), L=int((d_SA > 0).sum()),
            worst_degradation=float(d_SA.max()),
            worst_target=names[int(np.argmax(d_SA))],
            best_improvement=float(d_SA.min()),
            best_target=names[int(np.argmin(d_SA))],
            p10=float(np.percentile(d_SA, 10)), p90=float(np.percentile(d_SA, 90))),
        by_DISP_half=two_group(d_SA, folds, hi, "AVG_SEP.hi-minus-lo DISP of chain(SEP)-chain(AVG)"),
        mean_hi=float(d_SA[hi].mean()), mean_lo=float(d_SA[~hi].mean()),
        rho_DISP_vs_dSA=spearman(DISP, d_SA),
        rho_MOVE_vs_dSA=spearman(MOVE_A, d_SA),
        by_tail={k: dict(tail_mean=float(d_SA[m].mean()), rest_mean=float(d_SA[~m].mean()),
                         n_tail=int(m.sum()),
                         CAVEAT=("FAIL18 is defined by the filter's own recall "
                                 "(s12/instrument.py:271-278) and cannot measure it; diagnostic "
                                 "cross-check only" if k.startswith("FAIL18") else None))
                 for k, m in tails0.items()},
        gates=dict(
            REGISTERED="PREREG s12 (fifth amendment), appended at 93/126 rows with NO AVG_SEP "
                       "aggregate computed; both G4 and G5 are EXPLORATORY",
            G4=dict(rule="MOVE(AVG) > median -> AVG_SEP else AVG; native-free, 0 parameters",
                    cmp=brief(cmp2(np.where(MOVE_A > np.median(MOVE_A), ch_S, ch_A), ch_A,
                                   folds, names, "G4 MOVE gate -> AVG_SEP"))),
            G5=dict(rule="DISP > median -> AVG_SEP else AVG; native-free, 0 parameters",
                    cmp=brief(cmp2(np.where(hi, ch_S, ch_A), ch_A, folds, names,
                                   "G5 DISP gate -> AVG_SEP"))),
            G6=dict(rule="per-target min(AVG, AVG_SEP)",
                    ORACLE="ORACLE / NOT DEPLOYABLE -- best-of-2, not skill",
                    cmp=brief(cmp2(np.minimum(ch_S, ch_A), ch_A, folds, names,
                                   "G6 per-target min(AVG, AVG_SEP)")),
                    split_half_transfer=ST.split_half_transfer(
                        np.stack([ch_A, ch_S], 1), seed_parts=("s31F", str(SEED), "g6")))),
        coherence_cross_check=dict(
            source="s31/results/s31_F_coh.json",
            note="ORACLE / NOT DEPLOYABLE: coh(AVG_SEP) = 0.9689 against AVG's 0.9780 and the "
                 "0.6931 bar, 0/126 targets under the bar, and the direct affine-hull residual is "
                 "0.045 A RMS per coordinate. AVG_SEP does NOT leave the pool's affine hull, so "
                 "any endpoint gain it shows is NOT explained by lane B's non-affine mechanism."))
    out["physical_validity"] = dict(
        virtual_bond_A_from_the_separation_profile=dict(
            native=float(np.mean([r["prof_nat"][0] for r in R])),
            members=float(np.mean([r["prof_members"][0] for r in R])),
            cloud_AVG=float(np.mean([r["prof_cloud_AVG"][0] for r in R])),
            cloud_MED=float(np.mean([r["prof_cloud_MED"][0] for r in R])),
            cloud_AVG_SEP=float(np.mean([r["prof_cloud_AVG_SEP"][0] for r in R])),
            chain_AVG=float(np.mean([r["prof_chain_AVG"][0] for r in R])),
            chain_MED=float(np.mean([r["prof_chain_MED"][0] for r in R])),
            chain_AVG_SEP=float(np.mean([r["prof_chain_AVG_SEP"][0] for r in R]))),
        note="every CHAIN arm is the output of stage 3b, which is parameterised by (phi, psi) on "
             "IDEAL peptide geometry -- so peptide bond geometry, chirality and continuity are "
             "guaranteed BY CONSTRUCTION for all arms, and the virtual-bond column is the check "
             "that the builder was actually used. The CLOUD arms are NOT structures and their "
             "virtual bond says so; that is the point, not a defect.")

    # ---------------------------------------------------------------- F1: geometry vs accuracy
    rg_nat = col(R, "rg_nat")
    bA = np.array([r["bond_cloud_AVG"][0] for r in R]); bAs = np.array([r["bond_cloud_AVG"][1] for r in R])
    bM = np.array([r["bond_cloud_MED"][0] for r in R]); bMs = np.array([r["bond_cloud_MED"][1] for r in R])
    bcA = np.array([r["bond_chain_AVG"][0] for r in R]); bcAs = np.array([r["bond_chain_AVG"][1] for r in R])
    bcM = np.array([r["bond_chain_MED"][0] for r in R]); bcMs = np.array([r["bond_chain_MED"][1] for r in R])
    out["F1_geometry_vs_accuracy"] = dict(
        note="(a) accuracy = built-chain CA RMSD; (b) geometry = projection penalty P and "
             "Rg contraction. Reported separately, and the accuracy half governs the verdict.",
        contraction_cloud=dict(
            AVG_rg_over_nat=float((col(R, "rg_cloud_AVG") / rg_nat).mean()),
            MED_rg_over_nat=float((col(R, "rg_cloud_MED") / rg_nat).mean()),
            AVG_rg_over_members=float((col(R, "rg_cloud_AVG") / col(R, "rg_members_mean")).mean()),
            AVG_contraction_pct=float((1 - col(R, "rg_cloud_AVG") / col(R, "rg_members_mean")).mean() * 100),
            hi_DISP=float((1 - col(R, "rg_cloud_AVG")[hi] / col(R, "rg_members_mean")[hi]).mean() * 100),
            lo_DISP=float((1 - col(R, "rg_cloud_AVG")[~hi] / col(R, "rg_members_mean")[~hi]).mean() * 100)),
        bond_CA_CA=dict(cloud_AVG_mean=float(bA.mean()), cloud_AVG_sd=float(bAs.mean()),
                        cloud_MED_mean=float(bM.mean()), cloud_MED_sd=float(bMs.mean()),
                        chain_AVG_mean=float(bcA.mean()), chain_AVG_sd=float(bcAs.mean()),
                        chain_MED_mean=float(bcM.mean()), chain_MED_sd=float(bcMs.mean()),
                        note="a real backbone is ~3.80 A with a small sd; the AVG cloud is not a "
                             "structure and its bond statistics say so"),
        projection_penalty=dict(P_AVG_mean=float(P_A.mean()), P_MED_mean=float(P_M.mean()),
                                P_AVG_RG_mean=float(P_R.mean()),
                                P_AVG_hi=float(P_A[hi].mean()), P_AVG_lo=float(P_A[~hi].mean()),
                                contrast=two_group(P_A, folds, hi, "F1.geom hi-minus-lo DISP of P(AVG)")))

    # ---------------------------------------------------------------- the separation band
    # The CORRECTED averaging mechanism (prereg s10): a separation-dependent SHAPE distortion,
    # short contracted / long expanded, crossing unity near |i-j| = 8 -- not a uniform scale.
    SMAX = 15
    def ratio_tab(key, ref="prof_nat"):
        """mean over targets of prof_op[s]/prof_ref[s], per separation s (ragged n handled)."""
        acc = [[] for _ in range(SMAX)]
        for r in R:
            a = np.asarray(r[key], float); b = np.asarray(r[ref], float)
            for s in range(len(a)):
                acc[s].append(a[s] / b[s])
        return [float(np.mean(v)) if v else float("nan") for v in acc], \
               [len(v) for v in acc]

    prof_tab, prof_n = {}, None
    for key in ("prof_members", "prof_cloud_AVG", "prof_cloud_MED", "prof_cloud_AVG_SEP",
                "prof_chain_AVG", "prof_chain_MED", "prof_chain_AVG_SEP"):
        prof_tab[key], prof_n = ratio_tab(key)
    prof_nf = {}
    for key in ("prof_cloud_AVG", "prof_cloud_MED", "prof_chain_AVG", "prof_chain_MED"):
        prof_nf[key], _ = ratio_tab(key, ref="prof_members")

    def band_dev(r, key, lo, hi, ref="prof_nat"):
        """mean |ratio - 1| over separations in [lo, hi] for ONE target."""
        a = np.asarray(r[key], float); b = np.asarray(r[ref], float)
        s = np.arange(1, len(a) + 1)
        m = (s >= lo) & (s <= hi)
        return float(np.abs(a[m] / b[m] - 1.0).mean()) if m.any() else float("nan")

    bands = {}
    for key in ("prof_cloud_AVG", "prof_cloud_MED", "prof_chain_AVG", "prof_chain_MED",
                "prof_chain_AVG_SEP"):
        sh_ = np.array([band_dev(r, key, 1, 6) for r in R])
        lg_ = np.array([band_dev(r, key, 7, 99) for r in R])
        bands[key] = dict(short_1_6=float(np.nanmean(sh_)), long_7plus=float(np.nanmean(lg_)),
                          n_with_long=int(np.isfinite(lg_).sum()))
    dev_long_A = np.array([band_dev(r, "prof_chain_AVG", 7, 99) for r in R])
    dev_long_M = np.array([band_dev(r, "prof_chain_MED", 7, 99) for r in R])
    dev_short_A = np.array([band_dev(r, "prof_chain_AVG", 1, 6) for r in R])
    dev_short_M = np.array([band_dev(r, "prof_chain_MED", 1, 6) for r in R])
    okl = np.isfinite(dev_long_A) & np.isfinite(dev_long_M)
    NCOMP[0] += 2
    out["MECH_separation_band"] = dict(
        basis="ORACLE diagnostic (ratios are against the NATIVE's own per-separation distances); "
              "no parameter is tuned here",
        per_separation_vs_native=dict(s=list(range(1, SMAX + 1)), n_targets=prof_n, **prof_tab),
        per_separation_vs_members_NATIVE_FREE=dict(s=list(range(1, SMAX + 1)), **prof_nf),
        band_abs_deviation_from_native=bands,
        registered_falsifier="the medoid's profile must be FLATTER than the average's in the "
                             "s >= 7 band; if it is not, the separation-band framing is refuted "
                             "for this operator pair",
        long_band_MED_minus_AVG=float((dev_long_M - dev_long_A)[okl].mean()),
        short_band_MED_minus_AVG=float((dev_short_M - dev_short_A).mean()),
        cmp_long=brief(cmp2(dev_long_M[okl], dev_long_A[okl], folds[okl],
                            [names[i] for i in np.where(okl)[0]],
                            "MECH.long-band |ratio-1| MED - AVG (BUILT CHAIN, ORACLE)")),
        cmp_short=brief(cmp2(dev_short_M, dev_short_A, folds, names,
                             "MECH.short-band |ratio-1| MED - AVG (BUILT CHAIN, ORACLE)")),
        rho_long_band_distortion_vs_dMA=spearman(dev_long_A[okl], (ch_M - ch_A)[okl]),
        rho_short_band_distortion_vs_dMA=spearman(dev_short_A, ch_M - ch_A),
        withdrawn="the 25.8%% backbone contraction is WITHDRAWN (s15/coord_FINDINGS.md:914-921); "
                  "the corrected figure is 3.5%% against true distances and the distortion is "
                  "separation-dependent, not uniform")

    # ---------------------------------------------------------------- F2: the tail
    pool_mean, pool_best = col(R, "pool_mean"), col(R, "pool_best")
    set_mean, set_best = col(R, "set_mean"), col(R, "set_best")
    fail18 = np.array([r["fail18"] for r in R], bool)

    def worst18(v):
        m = np.zeros(len(R), bool)
        m[np.argsort(-v)[:18]] = True
        return m

    tails = {"T_POOL": worst18(pool_mean), "T_BEST": worst18(pool_best),
             "T_CHAIN": worst18(ch_A), "FAIL18": fail18}
    stats = dict(chain_AVG=ch_A, chain_MED=ch_M, chain_AVG_SEP=ch_S,
                 cloud_AVG=cl_A, MOVE_AVG=MOVE_A,
                 long_band_dev_AVG=dev_long_A, short_band_dev_AVG=dev_short_A,
                 pool_mean=pool_mean,
                 pool_best=pool_best, set_mean=set_mean, set_best=set_best,
                 DISP=DISP, S=S, s_b=s_b, P_AVG=P_A, P_MED=P_M,
                 rank_best_in_pool=col(R, "rank_best_in_pool").astype(float),
                 frac_top75_better_than_avg=col(R, "frac_top75_better_than_avg"),
                 n_top75_under3=col(R, "n_top75_under3").astype(float),
                 n_distinct=col(R, "n_distinct").astype(float),
                 d_MED_minus_AVG_chain=d_MA,
                 recall_gap_set_best_minus_pool_best=set_best - pool_best,
                 selection_gap_chain_minus_set_best=ch_A - set_best,
                 S_over_B=S / np.maximum(cl_A, 1e-9))
    F2 = {}
    for tn, tm in tails.items():
        ent = dict(members=sorted([names[i] for i in np.where(tm)[0]]),
                   overlap_with_T_CHAIN=int((tm & tails["T_CHAIN"]).sum()),
                   overlap_with_FAIL18=int((tm & fail18).sum()))
        if tn == "FAIL18":
            ent["CAVEAT"] = ("FAIL18 is defined by the FILTER's own recall (s12/instrument.py:"
                             "271-278) and therefore cannot measure filter recall; it appears "
                             "here as a diagnostic cross-check only")
        if tn == "T_BEST":
            ent["CAVEAT"] = "ORACLE stratum (defined on pool_best), but filter-independent"
        for sn, sv in stats.items():
            ent[sn] = dict(tail=float(sv[tm].mean()), rest=float(sv[~tm].mean()),
                           delta=float(sv[tm].mean() - sv[~tm].mean()))
        F2[tn] = ent
    out["F2_tails"] = F2
    out["F2_note"] = ("every statistic derived from rr / nat_ca is an ORACLE diagnostic; "
                      "nothing here tunes a parameter")

    # F2: the two gaps, decomposed on the built chain
    out["F2_decomposition"] = dict(
        basis="mixed: pool_best/set_best are CA POINT CLOUD ORACLE member RMSDs, chain_AVG is the "
              "BUILT CHAIN -- the two are NOT the same object and the gap is an upper bound",
        whole=dict(pool_best=float(pool_best.mean()), set_best=float(set_best.mean()),
                   set_mean=float(set_mean.mean()), cloud_AVG=float(cl_A.mean()),
                   chain_AVG=float(ch_A.mean())),
        T_CHAIN=dict(pool_best=float(pool_best[tails["T_CHAIN"]].mean()),
                     set_best=float(set_best[tails["T_CHAIN"]].mean()),
                     set_mean=float(set_mean[tails["T_CHAIN"]].mean()),
                     cloud_AVG=float(cl_A[tails["T_CHAIN"]].mean()),
                     chain_AVG=float(ch_A[tails["T_CHAIN"]].mean())),
        T_POOL=dict(pool_best=float(pool_best[tails["T_POOL"]].mean()),
                    set_best=float(set_best[tails["T_POOL"]].mean()),
                    set_mean=float(set_mean[tails["T_POOL"]].mean()),
                    cloud_AVG=float(cl_A[tails["T_POOL"]].mean()),
                    chain_AVG=float(ch_A[tails["T_POOL"]].mean())))

    # ---- F2: is the tail MORE filter-limited than readout-limited?  MEASURED, not tabulated.
    # Single basis throughout: pool_best, set_best and cloud_AVG are all CA POINT CLOUD.
    fl = set_best - pool_best                       # what the FILTER throws away
    rl = cl_A - set_best                            # what the READOUT loses on what it kept
    out["F2_filter_vs_readout_MEASURED"] = dict(
        basis="CA POINT CLOUD throughout -- pool_best, set_best and cloud_AVG are the same object, "
              "so nothing here is quoted across bases",
        ORACLE="ORACLE diagnostic: pool_best and set_best need the native",
        whole=dict(filter_loss=float(fl.mean()), readout_loss=float(rl.mean())),
        strata={tn: dict(filter_loss_tail=float(fl[tm].mean()),
                         filter_loss_rest=float(fl[~tm].mean()),
                         filter_growth=float(fl[tm].mean() / fl.mean()),
                         readout_loss_tail=float(rl[tm].mean()),
                         readout_loss_rest=float(rl[~tm].mean()),
                         readout_growth=float(rl[tm].mean() / rl.mean()),
                         excess_filter_over_readout=two_group(
                             fl - rl, folds, tm,
                             "F2.%s tail-minus-rest of (filter loss - readout loss)" % tn),
                         CAVEAT=("FAIL18 is the filter's own zero-recall set and cannot measure "
                                 "filter recall -- cross-check only"
                                 if tn.startswith("FAIL18") else
                                 ("defined by the OUTCOME, so partly downstream of the filter"
                                  if tn == "T_CHAIN" else "filter-independent stratum")))
                for tn, tm in tails.items()})

    # ---- the COMPLETE readout-choice ceiling over all four operators (UNREGISTERED arm)
    Mall = np.stack([ch_A, ch_M, ch_R, ch_S], 1)
    o4 = cmp2(Mall.min(1), ch_A, folds, names,
              "ORACLE min over all FOUR operators - AVG (BUILT CHAIN)")
    cnt = np.bincount(Mall.argmin(1), minlength=4)
    from scipy.stats import chisquare
    cs = chisquare(cnt)
    g4 = Mall.min(1) - ch_A
    out["ORACLE_readout_choice_ceiling"] = dict(
        UNREGISTERED="this arm was not pre-registered and is declared as such",
        ORACLE="ORACLE / NOT DEPLOYABLE -- a per-target minimum over K = 4, the construction "
               "S31-L11 audited; read with that entry's discipline",
        basis="BUILT CHAIN", arms=["AVG", "MED", "AVG_RG", "AVG_SEP"],
        cmp=brief(o4), mean=float(Mall.min(1).mean()),
        argmin_counts=cnt.tolist(), chi2=float(cs.statistic), p_vs_uniform=float(cs.pvalue),
        frac_targets_where_winner_is_not_AVG=float(1.0 - cnt[0] / len(R)),
        per_target=dict(mean=float(g4.mean()), median=float(np.median(g4)),
                        p90=float(np.percentile(g4, 90)), min=float(g4.min()),
                        n_tied_with_AVG=int((g4 == 0).sum())),
        split_half_transfer=dict(
            value=ST.split_half_transfer(Mall, seed_parts=("s31F", str(SEED), "of4")),
            DO_NOT_QUOTE="degenerate at small K when one column dominates globally -- the CI "
                         "collapses to a point; recorded so that nobody quotes it later"))

    out["multiplicity"] = dict(registered=44, emitted=int(NCOMP[0]))
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({k: out[k] for k in ("reproduction", "S2_theorem", "S2_precheck",
                                          "F1_operators", "F1_gates", "F1_gate_variables",
                                          "F1_mechanism", "AVG_SEP_mechanism",
                                          "physical_validity", "F1_geometry_vs_accuracy",
                                          "MECH_separation_band", "multiplicity")},
                     indent=1, default=float))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
