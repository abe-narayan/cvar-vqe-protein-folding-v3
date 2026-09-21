#!/usr/bin/env python
"""s31/s31_B_verify.py -- asserts every number lane B quotes in ledger S31-L10.

Reads only the committed artefacts. Any drift between the prose and the JSON fails loudly, so
the ledger entry cannot outlive the data it rests on.

    python s31/s31_B_verify.py        # exits non-zero on any mismatch
"""
import io, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "results")


def L(f):
    return json.load(io.open(os.path.join(R, f + ".json"), encoding="utf-8"))


CHECKS = []


def chk(lab, got, want, tol=5e-4):
    CHECKS.append((lab, got, want, tol, got is not None and abs(got - want) <= tol))


def main():
    a, s, p, b = (L("s31_B1_achirality"), L("s31_B2_sweep_parity"),
                  L("s31_B2_inpool"), L("s31_B3_graph"))
    Rw = {r["cost"]: r for r in s["rows"]}

    # --- B1: Lemma B1 ------------------------------------------------------------------
    chk("F2 every ff14SB torsion phase is 0 or pi", a["F2"]["max_dist_to_0_or_pi_rad"], 0.0, 1e-12)
    chk("F2 n_torsions", a["F2"]["n_torsions"], 1340, 0)
    chk("F2 no CMAP", float(a["F2"]["has_cmap"]), 0.0, 0)
    chk("F1 reflection residual", a["F1"]["max_rel"], 3.673e-06, 1e-8)
    chk("F1 ROTATION CONTROL is LARGER", a["F1"]["max_rel_rotation_control"], 5.377e-06, 1e-8)
    chk("F1 reflection/rotation", a["F1"]["reflection_over_rotation"], 0.68, 0.01)
    chk("F1 PeriodicTorsionForce exact", a["F1"]["max_rel_per_group"]["PeriodicTorsionForce"],
        1.98e-16, 1e-17)
    chk("F3 Legacy torsion odd share", a["F3"]["per_term"]["torsion"]["mean_odd_share"], 0.3006, 1e-3)
    chk("F3 Legacy steric odd share", a["F3"]["per_term"]["steric"]["mean_odd_share"], 13.016, 0.01)
    for t in ("hbond_local", "hbond_longrange", "coop_helix", "coop_sheet", "compactness"):
        chk("F3 %s is EXACTLY achiral" % t, a["F3"]["per_term"][t]["mean_odd_share"], 0.0, 1e-12)

    # --- the sweep reproduces the PUBLISHED S30 values ----------------------------------
    for c, v in (("LEG_steric", 0.25148809523809523), ("CAGEO", 0.23412698412698413),
                 ("LEG_torsion", 0.22123015873015872), ("DIS", 0.03571428571428571),
                 ("DSSPHB", -0.18353174603174602), ("RAMA", 0.1498015873015873)):
        chk("sweep %s == published S30" % c, Rw[c]["tot"]["effect"], v, 1e-9)
    chk("LEG_torsion odd contrast", Rw["LEG_torsion"]["odd"]["effect"], 0.2698, 1e-3)
    chk("LEG_torsion odd xMDE", Rw["LEG_torsion"]["odd"]["x_mde"], 2.78, 0.01)
    chk("LEG_torsion pref(circ_best vs PROD)", Rw["LEG_torsion"]["terms"]["tot"]["pref_best"], 0.4048, 1e-3)
    chk("LEG_steric pref(circ_best vs PROD)", Rw["LEG_steric"]["terms"]["tot"]["pref_best"], 0.5079, 1e-3)
    chk("LEG_torsion contrast vs GAUSS_MATCHED", Rw["LEG_torsion"]["terms"]["tot"]["contrast_vs_GAUSS_MATCHED"], 0.2222, 1e-3)
    chk("DIS odd part EXACTLY zero", Rw["DIS"]["odd"]["effect"], 0.0, 1e-15)
    chk("CONTACT odd part EXACTLY zero", Rw["CONTACT"]["odd"]["effect"], 0.0, 1e-15)

    # --- B2 in-pool ---------------------------------------------------------------------
    chk("M1 LEG_torsion odd variance share", p["M1_parity"]["leg_torsion_odd_share_mean"], 0.3001, 1e-3)
    chk("M2 RAND_SIGNED z in the pool", p["M2_saturation"]["z_RAND_SIGNED_mean"], 2.363, 1e-2)
    chk("M2 circ_best z in the pool", p["M2_saturation"]["z_circ_best_mean"], 1.129, 1e-2)
    chk("M0 mirror gap / pool sd", p["M0_mirror_gap"]["mean_gap_over_pool_sd"], 0.402, 1e-3)
    chk("M3 even half, 500 band", p["M3_rho"]["LEG_tors_even.p500"]["mean"], 0.2428, 1e-3)
    chk("M3 odd half, 500 band", p["M3_rho"]["LEG_tors_odd.p500"]["mean"], 0.0183, 1e-3)
    chk("M3 odd half xMDE, 500 band", p["M3_rho"]["LEG_tors_odd.p500"]["x_mde"], 0.28, 0.01)
    chk("G2 LEG_torsion in band xMDE", p["M3_rho"]["LEG_tors.sub75"]["x_mde"], 0.65, 0.01)
    chk("DIS in band xMDE (NOT MEASURED)", p["M3_rho"]["DIS.sub75"]["x_mde"], 0.83, 0.01)
    chk("HELIX_CONST rho, 500 band", p["M3_rho"]["HELIX_CONST.p500"]["mean"], 0.3413, 1e-3)
    chk("M6 XTWIST_4 rho, 500 band", p["M3_rho"]["XTWIST_4.p500"]["mean"], -0.0104, 1e-3)
    chk("M6 XTWABS_4 beats it, 500 band", p["M3_rho"]["XTWABS_4.p500"]["mean"], 0.1124, 1e-3)
    g = p["G5_vs_helix_control"]["LEG_tors.sub75"]
    chk("G5 in band effect", g["effect_channel_minus_control"], -0.0901, 1e-3)
    chk("G5 in band xMDE", g["x_mde"], -1.81, 0.01)
    chk("G5 in band fold CI upper < 0", g["ci95_fold"][1], -0.0639, 1e-3)

    # --- M5 coherence -------------------------------------------------------------------
    m = p["M5_coherence"]
    chk("M5 uniform mean coh EXACTLY 1", m["uniform_mean_pairspace"]["mean"], 1.0, 1e-12)
    chk("M5 argmin by the shipped DIS", m["argmin_DIS"]["mean"], 0.8288, 1e-3)
    chk("M5 argmin by a RANDOM member", m["argmin_RANDOM"]["mean"], 0.8244, 1e-3)
    chk("M5 ORACLE best member", m["ORACLE_best"]["mean"], 0.6708, 1e-3)
    chk("M5 coordinate average", m["coordinate_average"]["mean"], 0.9780, 1e-3)

    # --- B3 graph -----------------------------------------------------------------------
    chk("C1 geodesic vs direct", b["C1_collapse"]["geodesic"]["mean"], 0.8671, 1e-3)
    chk("C1 commute vs direct", b["C1_collapse"]["commute"]["mean"], 0.7916, 1e-3)
    chk("C2 medoid criterion | DIS", b["node_partial_on_medoid"]["medoid_crit"]["mean"], 0.2094, 1e-3)
    chk("C2 commute centrality | medoid", b["node_partial_on_medoid"]["commute_cent"]["mean"], 0.0121, 1e-3)
    chk("C2 geodesic centrality | medoid", b["node_partial_on_medoid"]["geo_cent"]["mean"], 0.0377, 1e-3)

    # --- provenance ---------------------------------------------------------------------
    chk("in-pool n = 126", p["n_targets"], 126, 0)
    chk("sweep n = 126", s["n"], 126, 0)
    chk("graph n = 126", b["n"], 126, 0)
    tot = (s["comparisons_emitted"] + p["comparisons_emitted"] + b["comparisons_emitted"])
    chk("lane B comparisons = 146", tot, 146, 0)

    bad = [c for c in CHECKS if not c[4]]
    for lab, got, want, tol, good in CHECKS:
        print("  %s  %-46s got %-14s want %s" % ("ok " if good else "FAIL", lab, got, want))
    print("\n%d/%d verified" % (len(CHECKS) - len(bad), len(CHECKS)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
