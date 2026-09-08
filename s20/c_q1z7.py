"""s20/c_q1z7.py -- DOES COMPACTNESS ALIGN WITH THE POOL'S OWN ERROR?  Testing Sprint 19's Z7.

DECLARED EXTENSION, written after `s20/results/c_q1_report.txt` section 4 was read and labelled
as such.  It was written after the Legacy-only/AMBER-only partition separated on RADIUS OF
GYRATION at -0.448 A [-0.526, -0.371], 124W/2L -- i.e. after Q1 established that **Legacy is a
compactness model** -- and it tests the standing claim that this is exactly the property that
makes a score harmful.

THE CLAIM UNDER TEST.  `s19/CLAIMS.md` Z7, label **PLAUSIBLE -- testable**:

    "Any score that prefers compact, well-formed, pool-typical geometry is selecting TOWARD the
     pool's own systematic error."

It was inferred from two facts (score gates raise ALIGN; the zero-information constant-alpha-helix
gate is the WORST arm measured) and was never measured directly.  It can be, on exactly the
instrument Q1 already built.

THE MEASUREMENT.  In the deployed averaging operator's own common frame, with `e_i = Y_i - T`,
`b_pool = mean_i e_i` and

    align_i = <e_i, b_pool> / (||b_pool|| sqrt(n))          [ORACLE, angstrom, mean_i = readout]

a gate that prefers HIGH-align candidates displaces the emitted mean along the direction the pool
is already wrong -- which is s19's C1/C2 mechanism, exactly.  So Z7 is the statement

    corr( compactness-preference , align )  >  0   across the candidates of a target,

and it is a per-candidate, per-target Spearman.  Scores tested, all native-free:

    rg            radius of gyration                 (LOW rg = compact -> expect rho(rg, align) < 0)
    z_leg         genuine Legacy total
    leg_compact   Legacy's `compactness` component alone
    z_amb         genuine AMBER single point
    e_disto       the shipped leave-fold-out distogram score
    helix_d       Ca-RMSD to a constant ideal alpha-helix  -- ZERO INFORMATION, and s19's WORST
                  gate.  If Z7 is right this should be the STRONGEST predictor of align of all.

FALSIFIER, registered here before the run.  Z7 is **REFUTED** if `rho(rg, align)` and
`rho(helix_d, align)` do not have the sign Z7 requires with fold-aware CIs excluding zero, or if
the zero-information `helix_d` is NOT among the strongest -- because Z7's whole force comes from
the zero-information gate being the worst one.  A matched control is unnecessary for a
correlation between two per-candidate quantities on a fixed candidate set, but a
LABEL-PERMUTATION null (shuffle `align` within the target) is reported beside every row.

    python -m s20.c_q1z7
"""
from __future__ import annotations

import os
import sys
import json
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

RESULTS = os.path.join(HERE, "results")
from s12 import instrument as I                     # noqa: E402
from s14.avgspace import top75_windows              # noqa: E402
from s15 import seed as SD                          # noqa: E402
from s16 import energy_lib as EL                    # noqa: E402
from s17 import phys_lib as P17                     # noqa: E402
from s18 import phys_lib as PL                      # noqa: E402
from s19 import agentC_lib as CL                    # noqa: E402
from s20 import c_q1 as Q1                          # noqa: E402

SCORES = ("rg", "leg_compact", "z_leg", "z_amb", "e_disto", "helix_d", "min_heavy")


# ---------------------------------------------------------------------------------------
# A SELF-CAUGHT LABELLING DEFECT.  `s12.instrument.paired`'s `ci95` is a PLAIN i.i.d.
# target-level bootstrap; its `folds` argument only adds a per-fold mean breakdown and does
# NOT cluster the resample.  Sprint 19's numbers were made with `s18.phys_lib.paired`, which
# returns BOTH an i.i.d. `ci` and a FOLD-CLUSTERED `ci_fold`, and "fold-aware" in this
# programme means the latter (BRIEF section 8).  Every CI below is therefore produced by
# `PL.paired` and the FOLD-CLUSTERED interval is the one quoted, with the i.i.d. one printed
# beside it wherever both matter.  Caught by reading `s12/instrument.py:188` rather than
# trusting the parameter name.
def PP(a, b, folds=None):
    st = PL.paired(np.asarray(a, float), np.asarray(b, float), folds=folds)
    return {"mean_diff": st["mean"], "ci95": st.get("ci_fold", st["ci"]),
            "ci_iid": st["ci"], "median_diff": st["median"],
            "n_better": st["W"], "n_worse": st["L"],
            "mean_a": st["mean_a"], "mean_b": st["mean_b"], "n": st["n"]}


def run(out="c_q1_z7.json"):
    from core import geometry as geo
    down = Q1.load_down()
    rows = []
    t0 = time.time()
    for t in I.targets():
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        W, PHI, PSI, u = top75_windows(pdb)
        W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
        nat = np.asarray(u["nat_ca"], float)
        K = len(W)
        Y, Ybar, _C, _d = CL.common_frame_members(W, PHI, PSI)
        T = I.superpose_batch(nat[None], Ybar)[0]
        e = Y - T
        b = e.mean(0)
        nb = float(np.sqrt((b ** 2).sum()))
        align = ((e.reshape(K, -1) @ b.reshape(-1)) / nb) / np.sqrt(n)   # ORACLE
        d_reb = np.asarray(I.kabsch_rmsd_batch(Y, nat), float)           # ORACLE

        dr = down[pdb]
        comp = {k: np.asarray(v, float) for k, v in dr["leg_terms"].items()}
        e_leg = np.asarray(EL.legacy_total_from(comp), float)
        e_amb = np.asarray(dr["s_amber_sp"], float)
        bb = geo.build_backbone_batch(PHI, PSI)
        CAb = np.asarray(bb["CA"], float)
        rg = np.sqrt(((CAb - CAb.mean(1, keepdims=True)) ** 2).sum(-1).mean(-1))
        hb = P17.helix_backbone(n)
        helix_d = np.asarray(I.kabsch_rmsd_batch(CAb, np.asarray(hb["CA"], float)), float)
        dg = I.distogram(pdb, seq, fold)
        ii, jj = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
        D = np.linalg.norm(CAb[:, ii, :] - CAb[:, jj, :], axis=-1)
        e_dis = np.asarray(I.shipped_score(dg, D), float)
        mh = np.array([EL.panel({a: np.asarray(bb[a][x], float) for a in bb}, seq)["min_heavy"]
                       for x in range(K)], float)

        sc = {"rg": rg, "leg_compact": comp["compactness"], "z_leg": e_leg, "z_amb": e_amb,
              "e_disto": e_dis, "helix_d": helix_d, "min_heavy": mh}
        rng = SD.stable_rng(pdb, "s20C_z7")
        rec = {"pdb": pdb, "n": n, "fold": fold, "readout": float(I.ca_rmsd(Ybar, nat)),
               "rho_align": {}, "rho_align_null_sd": {}, "rho_dreb": {},
               "rho_align_dreb": Q1._spearman(align, d_reb)}
        #: THE EXACT SET-LEVEL QUANTITY, which is what s19's ALIGN actually is.
        #: With b_S the survivors' mean error, Delta_S = b_S - b_pool and u = b_pool/||b_pool||,
        #: `s19.agentC_kv3.stats` computes  align = <Delta_S, u>/sqrt(n)  -- and since
        #: align_i := <e_i, u>/sqrt(n) here, that is EXACTLY
        #:
        #:      ALIGN_S = mean_S align_i  -  mean_K align_i          (identity, no free factor)
        #:
        #: A DEFECT I INTRODUCED AND CAUGHT BY READING s19's SOURCE RATHER THAN ITS PROSE
        #: (BRIEF section 10).  My first draft wrote `sqrt(n) * (mean_S align - readout)`,
        #: which is s19's quantity times sqrt(n) ~ 3.7 -- and it produced a helix ALIGN of
        #: +0.188 against s19's +0.0476, a 4x disagreement I would otherwise have had to
        #: explain away.  With the factor removed the two sprints agree; see the report.
        #: A rank correlation is only a proxy for this; the CUT is the operator.
        rec["align_excess"] = {}
        rec["align_excess_null"] = {}
        m = 38
        pool_mean = float(align.mean())
        for k, v in sc.items():
            rec["rho_align"][k] = Q1._spearman(v, align)
            rec["rho_dreb"][k] = Q1._spearman(v, d_reb)
            nul = [Q1._spearman(v, align[rng.permutation(K)]) for _ in range(50)]
            rec["rho_align_null_sd"][k] = float(np.std(nul))
            idx = CL.keep_lowest(v, m)
            rec["align_excess"][k] = float(align[idx].mean() - pool_mean)
        nulls = [float(align[rng.permutation(K)[:m]].mean() - pool_mean)
                 for _ in range(200)]
        rec["align_excess_null"] = {"mean": float(np.mean(nulls)), "sd": float(np.std(nulls))}
        rows.append(rec)
        if len(rows) % 20 == 0:
            print(f"  {len(rows)}/126 ({time.time()-t0:.0f}s)", flush=True)
    cfg = {"SCORES": list(SCORES), "n_perm": 50,
           "note": "declared extension; tests s19 CLAIMS Z7"}
    obj = {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg),
           "n_rows": len(rows), "n_expected": 126, "complete": len(rows) >= 126}
    with open(os.path.join(RESULTS, out), "w") as fh:
        json.dump(obj, fh)
    report(obj)
    return obj


def report(o=None, out="c_q1_z7.json"):
    if o is None:
        o = json.load(open(os.path.join(RESULTS, out)))
    rows = o["rows"]
    folds = np.array([r["fold"] for r in rows], int)
    print("=" * 106)
    print(f"Q1c  TESTING SPRINT 19's Z7 DIRECTLY:  does a compactness-preferring score select "
          f"TOWARD the pool's own error?   n = {len(rows)}")
    print("     align_i = <e_i, b_pool>/(||b_pool|| sqrt(n))  in the deployed operator's common "
          "frame.  ORACLE.")
    print("     A score whose LOW end (the end a gate keeps) has HIGH align is harmful.")
    print("=" * 106)
    print(f"\n{'score':<14s}{'rho(score, align)':>20s}{'CI95':>24s}{'median':>9s}"
          f"{'W/L':>10s}{'perm sd':>9s}{'rho(score, ORACLE d)':>22s}")
    for k in o["config"]["SCORES"]:
        v = np.array([r["rho_align"][k] for r in rows], float)
        d = np.array([r["rho_dreb"][k] for r in rows], float)
        ns = np.array([r["rho_align_null_sd"][k] for r in rows], float)
        st = PP(v, np.zeros_like(v), folds=folds)
        star = "  *" if (st["ci95"][0] > 0 or st["ci95"][1] < 0) else ""
        print(f"{k:<14s}{v.mean():+20.4f}  [{st['ci95'][0]:+9.4f},{st['ci95'][1]:+9.4f}]"
              f"{np.median(v):+9.4f}{int((v<0).sum()):>5d}/{int((v>0).sum()):<4d}"
              f"{ns.mean():9.4f}{d.mean():+22.4f}{star}")
    print(f"\n  THE EXACT SET-LEVEL QUANTITY -- s19's ALIGN itself, not a rank proxy.")
    print(f"  ALIGN_S = <b_S - b_pool, u>/sqrt(n) = mean_S align - mean_K align   (s19 kv3 units), "
          f"m = 38 of 75.")
    print(f"  POSITIVE = the gate pushed the emitted mean ALONG the pool's own error.")
    nl = np.array([r["align_excess_null"]["mean"] for r in rows], float)
    print(f"\n{'score':<14s}{'ALIGN_S (A)':>14s}{'excess over matched-random':>30s}"
          f"{'median':>9s}{'W/L':>10s}")
    for k in o["config"]["SCORES"]:
        v = np.array([r["align_excess"][k] for r in rows], float)
        st = PP(v, nl, folds=folds)
        star = "  *" if (st["ci95"][0] > 0 or st["ci95"][1] < 0) else ""
        print(f"{k:<14s}{v.mean():+14.4f}{st['mean_diff']:+14.4f}"
              f"  [{st['ci95'][0]:+7.4f},{st['ci95'][1]:+7.4f}]{st['median_diff']:+9.4f}"
              f"{st['n_better']:>5d}/{st['n_worse']:<4d}{star}")
    print(f"{'matched-random':<14s}{nl.mean():+14.4f}   (200 draws per target)")
    ad = np.array([r["rho_align_dreb"] for r in rows], float)
    st = PP(ad, np.zeros_like(ad), folds=folds)
    print(f"\n    reference: rho(align, ORACLE d_reb) = {ad.mean():+.4f} "
          f"[{st['ci95'][0]:+.4f}, {st['ci95'][1]:+.4f}] -- how much 'being aligned with the "
          f"pool's error' is just 'being a bad candidate'")
    print("\n    READ THE SIGNS.  A gate keeps the LOWEST scores.  So a score is harmful when a")
    print("    LOW score goes with a HIGH align, i.e. when rho(score, align) is NEGATIVE.")
    print("=" * 106)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--report":
        report()
    else:
        run()
