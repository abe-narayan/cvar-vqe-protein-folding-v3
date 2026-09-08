"""SPRINT 20 / WORKSTREAM B -- Q2: does a genuine CVaR-VQE's coherent error DECORRELATE
from the retrieval pool's?

    python -m s20.qb2_q2 run  [n]
    python -m s20.qb2_q2 report

THE DECOMPOSITION, AND AN IDENTITY STATED SO IT IS NOT MISTAKEN FOR A FINDING.
Sprint 19 (`s19/a_source.py`) splits a *predictor's* error exactly:

    r_F = dhat_F - d_true = (d(X_F) - d_true) + (dhat_F - d(X_F)) = r_coh + r_inc

For a *generator* there is no `dhat`: the family emits a structure, and a structure realises its
own distances exactly.  So `r_inc == 0` and the generator's whole error IS its coherent
component.  **EXACT, an identity.**  The live question is therefore only whether
`rho(r_coh_VQE, r_coh_pool)` is low.

THE CEILING.  "Shared" is never judged against zero.  Two independently seeded runs of the SAME
generator set the ceiling for sharing-because-same-method, and every cross-family number is read
against THAT (the rule `s19/a_source` established).

THE FALSIFIER, adopted verbatim from `s20/BRIEF.md`: if the VQE candidates align with the pool at
the cross-architecture level (~0.75+), quantum sampling is another view of the same bias and the
generator role closes.

Error diversity is reported SEPARATELY from geometric diversity (BRIEF section 4): two
populations can be geometrically different and make the same structural mistake.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

from s20 import qb2_lib as L
from s20 import qb2_opt as OP
from s20 import qb2_run as R
from s12 import instrument as I
from s15 import seed as SD
from s19 import qb_lib as QB
from core import quantum as Q

BUDGET = 4096
M_SEL = 75
ALPHA = 0.25


def _terminal(Zset, tgt, m=M_SEL):
    """The FROZEN native-free terminal operator: shipped Bayes-risk top-m, coordinate average,
    ideal-geometry projection.  `avg` is the POINT CLOUD, `ca` is the BUILT chain -- both are
    returned and every downstream number states which basis it is on (the L20 rule)."""
    from core import project as pj
    phi, psi = L.unpack(Zset, tgt["n"])
    CA = np.asarray(pj.build_ca_exact(phi, psi), float)
    D = I.pair_dists(CA, tgt["i"], tgt["j"])
    sc = I.shipped_score(tgt["dg"], D)
    sel = np.argsort(sc, kind="stable")[:min(m, len(CA))]
    W = CA[sel]
    avg, _b = I.coordinate_average(W)
    pr = I.project(np.asarray(avg, float), tgt["seq"], tgt["fold"])
    return {"avg": np.asarray(avg, float), "ca": np.asarray(pr["ca"], float),
            "sel": W}


def _err(ca, tgt):
    """The realisable (coherent) distance error of a STRUCTURE, over CA pairs |i-j| >= 2."""
    d = np.linalg.norm(ca[tgt["i"]] - ca[tgt["j"]], axis=-1)
    dt = np.linalg.norm(tgt["nat"][tgt["i"]] - tgt["nat"][tgt["j"]], axis=-1)
    return d - dt


def _emit_vqe(tgt, seed, budget=BUDGET, alpha=ALPHA, entangler="cnot", train=True):
    h = L.Ham("DIST", tgt, budget=budget, keep=True)
    h.med_ref, h.iqr_ref = 0.0, 1.0
    F = OP.Field(h, tgt); F.set_pen(1e12)
    rng = SD.stable_rng(tgt["pdb"], "q2vqe", seed, alpha, entangler, train, salt=L.SALT)
    OP.arm_vqe(F, None, rng, alpha=alpha, shots=64, train=train, entangler=entangler)
    phi, psi, _e = h.seen()
    return L.pack(phi, psi)


def _emit_bestofn(tgt, seed, budget=BUDGET):
    rng = SD.stable_rng(tgt["pdb"], "q2bon", seed, salt=L.SALT)
    return L.basin_starts(tgt, budget, rng)


def _emit_metro(tgt, seed, budget=BUDGET):
    h = L.Ham("DIST", tgt, budget=budget, keep=True)
    h.med_ref, h.iqr_ref = 0.0, 1.0
    F = OP.Field(h, tgt); F.set_pen(1e12)
    rng = SD.stable_rng(tgt["pdb"], "q2metro", seed, salt=L.SALT)
    OP.arm_metro(F, None, rng)
    phi, psi, _e = h.seen()
    return L.pack(phi, psi)


def _emit_helix(tgt, seed, budget=500):
    """ZERO-INFORMATION, plausible: a constant ideal alpha-helix plus a matched isotropic
    torsion perturbation.  Uniform-on-the-torus is NOT used (BRIEF section 6 rule 4)."""
    rng = SD.stable_rng(tgt["pdb"], "q2helix", seed, salt=L.SALT)
    n = tgt["n"]
    phi = np.full((budget, n), np.deg2rad(-57.0)) + rng.normal(0, 0.35, (budget, n))
    psi = np.full((budget, n), np.deg2rad(-47.0)) + rng.normal(0, 0.35, (budget, n))
    return L.pack(L.wrap(phi), L.wrap(psi))


FAMILIES = {
    "pool500": lambda t, s: L.pack(t["PHI"], t["PSI"]),
    "vqe": lambda t, s: _emit_vqe(t, s),
    "vqe_prod": lambda t, s: _emit_vqe(t, s, entangler="none"),
    "vqe_untrained": lambda t, s: _emit_vqe(t, s, train=False),
    "bestofN": lambda t, s: _emit_bestofn(t, s),
    "metro": lambda t, s: _emit_metro(t, s),
    "helix": lambda t, s: _emit_helix(t, s),
}
SEEDS = (0, 1)


def run(pdbs):
    tag = "q2"
    out = R.ck_load(tag)
    for ti, pdb in enumerate(pdbs):
        if pdb in out:
            continue
        t0 = time.time()
        tgt = L.target(pdb)
        rec = {"n": tgt["n"], "fold": tgt["fold"], "fam": {}}
        for fam, fn in FAMILIES.items():
            for s in SEEDS:
                Z = fn(tgt, s)
                if len(Z) == 0:
                    continue
                T = _terminal(Z, tgt)
                e = _err(T["ca"], tgt)
                # geometric diversity of the SELECTED set, and its member quality (ORACLE)
                rr = I.kabsch_rmsd_batch(T["sel"], tgt["nat"])
                P = I.pairwise_rmsd(T["sel"][:40])
                rec["fam"][f"{fam}|{s}"] = {
                    "err": e.tolist(),
                    "rmsd_built_ORACLE": float(I.ca_rmsd(T["ca"], tgt["nat"])),
                    "rmsd_cloud_ORACLE": float(I.ca_rmsd(T["avg"], tgt["nat"])),
                    "member_mean_ORACLE": float(rr.mean()),
                    "member_best_ORACLE": float(rr.min()),
                    "geom_diversity": float(P[np.triu_indices(len(P), 1)].mean()),
                    "n_emitted": int(len(Z)),
                }
        out[pdb] = rec
        R.ck_save(tag, out)
        R.log(tag, f"[{ti+1}/{len(pdbs)}] {pdb} n={tgt['n']} {time.time()-t0:.0f}s "
              + " ".join(f"{f}:{rec['fam'][f+'|0']['rmsd_built_ORACLE']:.2f}"
                         for f in FAMILIES if f + "|0" in rec["fam"]))
    with open(os.path.join(L.RESULTS, "qb2_q2_COMPLETE"), "w") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S"))
    return out


def _corr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 4:
        return np.nan
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt((a ** 2).sum() * (b ** 2).sum())
    return float((a * b).sum() / d) if d > 0 else np.nan


def report():
    d = R.ck_load("q2")
    if not d:
        print("no qb2_q2.json"); return {}
    pdbs = sorted(d)
    folds = [d[p]["fold"] for p in pdbs]
    print(f"\n=== Q2 DECORRELATION, n = {len(pdbs)} targets ===")
    print("rho of the REALISABLE (coherent) distance error of the BUILT chain, per target,")
    print("against the retrieval pool's own.  Judged against the same-generator ceiling.")
    out = {"n": len(pdbs), "targets": pdbs}
    # the ceiling: the SAME generator, two seeds
    print("\n--- the CEILING: same generator, two independent seeds ---")
    out["ceiling"] = {}
    for fam in FAMILIES:
        v = [_corr(d[p]["fam"][f"{fam}|0"]["err"], d[p]["fam"][f"{fam}|1"]["err"])
             for p in pdbs if f"{fam}|0" in d[p]["fam"] and f"{fam}|1" in d[p]["fam"]]
        if v:
            ci = L.boot_mean_ci(np.array(v))
            out["ceiling"][fam] = ci
            print(f"  {fam:<15} {ci['mean']:+.3f}  CI [{ci['ci95'][0]:+.3f}, {ci['ci95'][1]:+.3f}]")
    print("\n--- vs the RETRIEVAL POOL (the falsifier line is 0.75) ---")
    out["vs_pool"] = {}
    for fam in FAMILIES:
        if fam == "pool500":
            continue
        v = [_corr(d[p]["fam"][f"{fam}|0"]["err"], d[p]["fam"]["pool500|0"]["err"])
             for p in pdbs if f"{fam}|0" in d[p]["fam"]]
        if v:
            ci = L.boot_mean_ci(np.array(v))
            out["vs_pool"][fam] = ci
            fired = "FALSIFIER FIRES" if ci["ci95"][0] >= 0.75 else (
                "below 0.75" if ci["ci95"][1] < 0.75 else "spans 0.75")
            print(f"  {fam:<15} {ci['mean']:+.3f}  CI [{ci['ci95'][0]:+.3f}, {ci['ci95'][1]:+.3f}]"
                  f"   {fired}")
    print("\n--- quality and diversity, reported SEPARATELY from error correlation ---")
    print(f"  {'family':<15}{'built RMSD':>12}{'cloud RMSD':>12}{'member mean':>13}"
          f"{'member best':>13}{'geom div':>10}")
    out["quality"] = {}
    for fam in FAMILIES:
        k = f"{fam}|0"
        rows = [d[p]["fam"][k] for p in pdbs if k in d[p]["fam"]]
        if not rows:
            continue
        q = {h: L.boot_mean_ci(np.array([r[h] for r in rows]))
             for h in ("rmsd_built_ORACLE", "rmsd_cloud_ORACLE", "member_mean_ORACLE",
                       "member_best_ORACLE", "geom_diversity")}
        out["quality"][fam] = q
        print(f"  {fam:<15}{q['rmsd_built_ORACLE']['mean']:>12.3f}"
              f"{q['rmsd_cloud_ORACLE']['mean']:>12.3f}"
              f"{q['member_mean_ORACLE']['mean']:>13.3f}"
              f"{q['member_best_ORACLE']['mean']:>13.3f}"
              f"{q['geom_diversity']['mean']:>10.3f}")
    # paired: every family against the pool on BUILT RMSD
    print("\n--- built-chain RMSD vs pool500, paired (negative = the family wins) ---")
    out["rmsd_vs_pool"] = {}
    base = np.array([d[p]["fam"]["pool500|0"]["rmsd_built_ORACLE"] for p in pdbs])
    for fam in FAMILIES:
        if fam == "pool500":
            continue
        k = f"{fam}|0"
        if not all(k in d[p]["fam"] for p in pdbs):
            continue
        v = np.array([d[p]["fam"][k]["rmsd_built_ORACLE"] for p in pdbs])
        t = L.paired_ci(v, base, folds=folds)
        out["rmsd_vs_pool"][fam] = t
        print(f"  {fam:<15} {t['mean']:+.3f} [{t['ci95'][0]:+.3f},{t['ci95'][1]:+.3f}] "
              f"W/L {t['W']}/{t['L']} sig={int(t['sig'])}")
    # ---- THE NULL-REFERENCED TEST, and a pre-registration mis-specification.
    # PREREG_B section 5 adopted the brief's absolute falsifier line of 0.75.  That line was
    # calibrated on cross-architecture PREDICTOR families.  Measured here, the ZERO-INFORMATION
    # constant helix already correlates with the pool at ~0.79 -- i.e. the threshold sits BELOW
    # its own zero-information null for GENERATORS, so the absolute test cannot discriminate.
    # The pre-registration is kept unedited; the valid test is against the measured null and
    # against each generator's OWN same-seed ceiling, and both are reported here.
    print("")
    print("--- NULL-REFERENCED (the valid test): rho(F, pool) minus a reference's ---")
    out["null_referenced"] = {}
    v = np.array([_corr(d[p]["fam"]["vqe|0"]["err"], d[p]["fam"]["pool500|0"]["err"])
                  for p in pdbs])
    ceil = np.array([_corr(d[p]["fam"]["vqe|0"]["err"], d[p]["fam"]["vqe|1"]["err"])
                     for p in pdbs])
    for lab, ref in (("helix (ZERO-INFORMATION)", "helix"), ("metro (classical thermostat)", "metro"),
                     ("bestofN (no optimisation)", "bestofN"), ("vqe_prod (CNOTs deleted)", "vqe_prod")):
        r = np.array([_corr(d[p]["fam"][f"{ref}|0"]["err"], d[p]["fam"]["pool500|0"]["err"])
                      for p in pdbs])
        t = L.paired_ci(v, r, folds=folds)
        out["null_referenced"][ref] = t
        print(f"  vqe-vs-pool  minus  {lab:<32} {t['mean']:+.4f} "
              f"[{t['ci95'][0]:+.4f},{t['ci95'][1]:+.4f}] W/L {t['W']}/{t['L']} sig={int(t['sig'])}")
    t = L.paired_ci(v, ceil, folds=folds)
    out["null_referenced"]["own_ceiling"] = t
    print(f"  vqe-vs-pool  minus  its OWN same-seed ceiling      {t['mean']:+.4f} "
          f"[{t['ci95'][0]:+.4f},{t['ci95'][1]:+.4f}] W/L {t['W']}/{t['L']} sig={int(t['sig'])}")
    out["shared_fraction"] = {
        "vqe_cross_over_own_ceiling": float(v.mean() / ceil.mean()),
        "vqe_per_target_ratio_mean": float(np.mean(v / ceil))}
    print("")
    print(f"  >>> share of the VQE's OWN reproducible error structure that is shared with the "
          f"retrieval pool: {v.mean()/ceil.mean():.4f} "
          f"(per-target ratio mean {np.mean(v/ceil):.4f})")
    L.write("qb2_report_q2", out, complete=True)
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    n_t = int(sys.argv[2]) if len(sys.argv) > 2 else L.N_SUBSET
    if mode == "run":
        run([t["pdb"] for t in L.subset(L.N_SUBSET)[:n_t]])
    report()


if __name__ == "__main__":
    main()
