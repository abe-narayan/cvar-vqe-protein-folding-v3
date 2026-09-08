"""Aggregate analysis over the A1/A2/A3 artefacts, for `s22/agentA_FINDINGS.md`.

Run after all three artefacts exist. Prints everything FINDINGS needs; nothing here re-runs any
training -- it is pure post-hoc arithmetic over the persisted JSON.
"""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")


def load(name):
    with open(os.path.join(RES, name)) as fh:
        return json.load(fh)


def paired_ci(a, b, n_boot=4000, seed=0):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    d = a - b
    n = len(d)
    if n < 3:
        return {"n": int(n), "mean": float(d.mean()) if n else float("nan")}
    rng = np.random.default_rng(seed)
    bs = np.array([d[rng.integers(0, n, n)].mean() for _ in range(n_boot)])
    se = float(d.std(ddof=1) / math.sqrt(n))
    return {"n": int(n), "mean": float(d.mean()), "median": float(np.median(d)), "se": se,
           "mde": float(2.8016 * se), "ci95": [float(np.percentile(bs, 2.5)),
                                                float(np.percentile(bs, 97.5))],
           "W": int((d < 0).sum()), "L": int((d > 0).sum())}


def section_a1():
    print("=" * 70, "\nA1\n", "=" * 70)
    d = load("a1_gauge.json")
    rows = d["rows"]
    g1_fired = sum(r["gate1_fired"] for r in rows)
    g1_checks = sum(r["gate1_checks"] for r in rows)
    print(f"GATE 1 (soundness, corrected): {g1_fired}/{g1_checks}")
    g2 = {"argmin_rr": 0, "cvar": 0, "tailavg": 0}
    g2c = 0
    for r in rows:
        for k in g2:
            g2[k] += r["gate2_fired"][k]
        g2c += r["gate2_checks"]
    print(f"GATE 2 (achieved gauge, corrected), fired/{g2c}:", g2)

    diffs_tailavg, diffs_cvar, diffs_argmin = [], [], []
    for r in rows:
        cm = r["canon_mean"]
        for p in r["perm_rows"]:
            diffs_tailavg.append(p["r_tailavg"] - cm["tailavg"])
            diffs_cvar.append(p["final_cvar"] - cm["cvar"])
            diffs_argmin.append(p["argmin_rr"] - cm["argmin_rr"])
    diffs_tailavg = np.array(diffs_tailavg); diffs_cvar = np.array(diffs_cvar)
    diffs_argmin = np.array(diffs_argmin)
    print(f"n perm-cells = {len(diffs_tailavg)}")
    print(f"tailavg gauge deviation: mean {diffs_tailavg.mean():.4f} sd {diffs_tailavg.std():.4f} "
         f"mean|d| {np.abs(diffs_tailavg).mean():.4f} median|d| {np.median(np.abs(diffs_tailavg)):.4f}")
    for thr in (0.05, 0.2, 0.5):
        print(f"  pct |d|>{thr}A: {(np.abs(diffs_tailavg) > thr).mean():.3f}")
    print(f"cvar gauge deviation: mean|d| {np.abs(diffs_cvar).mean():.5f}")
    print(f"argmin_rr gauge deviation: mean|d| {np.abs(diffs_argmin).mean():.5f} "
         f"(expect ~0, argmin set is gauge-invariant post-fix)")

    # RMSD decomposition
    R_pool = np.array([r["R_pool"] for r in rows])
    R_score_argmin = np.array([r["R_score_argmin"] for r in rows])
    R_score_tailavg = np.array([r["R_score_tailavg"] for r in rows])
    R_vqe_argmin = np.array([r["canon_mean"]["argmin_rr"] for r in rows])
    R_vqe_tailavg = np.array([r["canon_mean"]["tailavg"] for r in rows])
    R_rand_tailavg = np.array([r["random_control"]["r_tailavg_mean"] for r in rows])
    R_rand_argmin = np.array([r["random_control"]["r_argmin_mean"] for r in rows])
    R_untrained_argmin = np.array([r["untrained_argmin_rr"] for r in rows])
    R_untrained_tailavg = np.array([r["untrained_tailavg_rr"] for r in rows])
    print("\nRMSD decomposition, mean over 16 targets:")
    for name, arr in [("R_pool (ORACLE ceiling)", R_pool),
                      ("R_score_argmin (classical exact)", R_score_argmin),
                      ("R_VQE_argmin (trained, canonical)", R_vqe_argmin),
                      ("R_score_tailavg (classical exact top-a avg)", R_score_tailavg),
                      ("R_VQE_tailavg (trained, canonical, T=0)", R_vqe_tailavg),
                      ("R_rand_argmin (matched-count random)", R_rand_argmin),
                      ("R_rand_tailavg (matched-count random)", R_rand_tailavg),
                      ("R_untrained_argmin (best_of_N, untrained circuit)", R_untrained_argmin),
                      ("R_untrained_tailavg (best_of_N, untrained circuit)", R_untrained_tailavg)]:
        print(f"  {name:52s} {np.mean(arr):.3f}")

    print("\nVQE_argmin vs classical exact argmin (should be ~0, theorem):",
         paired_ci(R_vqe_argmin, R_score_argmin))
    print("VQE_tailavg vs classical exact tailavg:", paired_ci(R_vqe_tailavg, R_score_tailavg))
    print("VQE_tailavg vs matched random control:", paired_ci(R_vqe_tailavg, R_rand_tailavg))
    print("classical exact tailavg vs matched random control:",
         paired_ci(R_score_tailavg, R_rand_tailavg))

    sa_found = np.mean([r["sa_control"]["found_global"] for r in rows])
    print(f"\nSA control: found global optimum in {sa_found*100:.0f}% of targets "
         f"(budget = K = 500 evals)")

    # legacy secondary
    n_leg = sum(1 for r in rows if r.get("legacy"))
    print(f"\nH_Legacy secondary pillar: ran on {n_leg} targets")
    return d


def section_a2a():
    print("\n" + "=" * 70, "\nA2a (collapse dynamics)\n", "=" * 70)
    d = load("a2a_collapse.json")
    rows = d["rows"]
    cks = d["checkpoints"]
    alphas = d["alphas"]
    # ESS and tail_size vs iteration, pooled over targets/seeds, per alpha
    for alpha in alphas:
        print(f"\nalpha={alpha}")
        for ck in cks:
            ess_vals, ts_vals, rmsd_vals = [], [], []
            for r in rows:
                for cell in r["cells"]:
                    if cell["alpha"] != alpha:
                        continue
                    snap = cell["checkpoints"].get(str(ck)) or cell["checkpoints"].get(ck)
                    if snap is None:
                        continue
                    ess_vals.append(snap["ess"])
                    ts_vals.append(snap["tail_size"])
                    rmsd_vals.append(snap["r_tailavg"])
            if ess_vals:
                print(f"  it={ck:4d}  ess={np.mean(ess_vals):7.2f}  "
                     f"tail_size={np.mean(ts_vals):6.1f}  "
                     f"tailavg_rmsd={np.nanmean(rmsd_vals):.3f}")
    return d


def section_a2b():
    print("\n" + "=" * 70, "\nA2b (entropy-regularised degeneracy)\n", "=" * 70)
    d = load("a2b_entropy.json")
    rows = d["rows"]
    Ts = d["Ts"]

    # Rebuild each target's pool once, for the classical size-matched top-m control (the bar
    # the coordinator specified: NOT random, since s21 priced "bigger tail averages better"
    # at -0.28 to -0.31 A independent of the energy -- that is the null, not the finding).
    import s22.qcand_lib as QC
    pools = {r["pdb"]: QC.build_pool(r["pdb"]) for r in rows}

    for entangler in ("cnot", "none"):
        print(f"\nentangler={entangler}")
        for T in Ts:
            ess_vals, ts_vals, rmsd_vals, rand_rmsd, classical_rmsd = [], [], [], [], []
            target_level = {}   # pdb -> list of (vqe_rmsd, classical_rmsd) across seeds
            for r in rows:
                pool = pools[r["pdb"]]
                for cell in r["cells"]:
                    if cell["entangler"] != entangler or cell["T"] != T:
                        continue
                    m = max(1, cell["n_tail_cands"])
                    cls = QC.classical_exact_sort(pool, pool["score_dist"], m / len(pool["score_dist"]))
                    ess_vals.append(cell["ess"])
                    ts_vals.append(cell["tail_size"])
                    rmsd_vals.append(cell["r_tailavg"])
                    rand_rmsd.append(cell["random_control_matched_m"]["r_tailavg_mean"])
                    classical_rmsd.append(cls["r_tailavg"])
                    target_level.setdefault(r["pdb"], []).append(
                        (cell["r_tailavg"], cls["r_tailavg"]))
            if not ess_vals:
                continue
            ci_rand = paired_ci(rmsd_vals, rand_rmsd)
            ci_cls = paired_ci(rmsd_vals, classical_rmsd)
            # target-level (mean over seeds per target) for the PRIMARY comparison, per PREREG
            tv, cv = [], []
            for pdb, lst in target_level.items():
                tv.append(np.mean([x[0] for x in lst]))
                cv.append(np.mean([x[1] for x in lst]))
            ci_cls_target = paired_ci(tv, cv)
            print(f"  T={T:5.2f} ess={np.mean(ess_vals):7.2f} tail_size={np.mean(ts_vals):6.1f} "
                 f"VQEtailavg={np.nanmean(rmsd_vals):.3f} rand={np.nanmean(rand_rmsd):.3f} "
                 f"classical_topm={np.nanmean(classical_rmsd):.3f}")
            print(f"      vs RANDOM (matched m):   mean={ci_rand.get('mean'):+.4f} "
                 f"MDE={ci_rand.get('mde'):.4f} CI={ci_rand.get('ci95')} "
                 f"W/L={ci_rand.get('W')}/{ci_rand.get('L')}  [cell-level, n={ci_rand.get('n')}]")
            print(f"      vs CLASSICAL top-m (matched m, PRIMARY BAR): "
                 f"mean={ci_cls_target.get('mean'):+.4f} MDE={ci_cls_target.get('mde'):.4f} "
                 f"CI={ci_cls_target.get('ci95')} W/L={ci_cls_target.get('W')}/{ci_cls_target.get('L')} "
                 f"[target-level, n={ci_cls_target.get('n')}]")
    return d


def section_a3():
    print("\n" + "=" * 70, "\nA3 (multi-stage)\n", "=" * 70)
    d = load("a3_multistage.json")
    print("PRIMARY staged - single (tail-average RMSD):", d["primary_staged_minus_single_tailavg"])
    print("SECONDARY staged - single (argmin RMSD):", d["secondary_staged_minus_single_argmin"])
    print("argmin tied to exact optimum:", d["argmin_tied_exact_optimum"])
    return d


if __name__ == "__main__":
    section_a1()
    section_a2a()
    section_a2b()
    section_a3()
