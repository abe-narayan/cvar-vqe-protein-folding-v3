"""SPRINT 19 / AGENT B -- the report.  Every number in `agentB_FINDINGS.md` is printed here.

    python -m s19.qb_report main      # the arm table and the pre-registered P1/P2/P3 tests
    python -m s19.qb_report theory    # T1 (CVaR indifference) and T2 (the restricted form)
    python -m s19.qb_report gauge     # the Z2^n gauge orbit
    python -m s19.qb_report law       # does the terminal operator consume the mean here too?
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I           # noqa: E402
from s17 import q_lib as QL               # noqa: E402
from s19 import qb_lib as L               # noqa: E402
from s19 import qb_arms as AR             # noqa: E402

MDE = 0.084


def rows(tag="main"):
    d = L.load(tag)
    r = [v for k, v in d.items() if not k.startswith("_")]
    return sorted(r, key=lambda x: x["pdb"])


def col(rs, arm, key):
    return np.array([r["arms"][arm][key] for r in rs], float)


def verdict(a, b, folds, label, rng_seed=0):
    """Paired a - b (negative = a better), with everything brief section 9 demands."""
    d = np.asarray(a, float) - np.asarray(b, float)
    p = I.paired(np.asarray(a, float), np.asarray(b, float), seed=rng_seed)
    cl = QL.cluster_paired(d, folds, seed=rng_seed)
    sig = (p["ci95"][0] > 0) or (p["ci95"][1] < 0)
    return {"label": label, "n": p["n"], "mean_a": p["mean_a"], "mean_b": p["mean_b"],
            "mean": p["mean_diff"], "median": p["median_diff"], "ci95": p["ci95"],
            "W": p["n_better"], "L": p["n_worse"], "sig": bool(sig),
            "cluster_ci95": cl["cluster_ci95"], "folds_same_sign": cl["folds_same_sign"],
            "above_MDE": bool(abs(p["mean_diff"]) >= MDE)}


def _line(v):
    return (f"  {v['label']:<44}{v['mean']:+8.3f}  med {v['median']:+7.3f}  "
            f"[{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]  "
            f"W/L {v['W']:3d}/{v['L']:3d}  folds {v['folds_same_sign']}/5"
            f"{'  SIG' if v['sig'] else ''}{'' if v['above_MDE'] else '  <MDE'}")


def main(tag="main"):
    rs = rows(tag)
    folds = np.array([r["fold"] for r in rs], int)
    names = list(AR.arms().keys())
    print(f"\nn = {len(rs)} targets.  Budget = {L.BUDGET} objective evaluations for every "
          f"budgeted arm.\nCONTINUOUS TORSION SPACE.  No lattice, no torsion discretisation.")
    print("Every RMSD column is ORACLE (post-hoc scoring of a native-free decision).\n")
    hdr = (f"  {'arm':<14}{'evals':>7}{'gen best':>10}{'sel best':>10}{'realised':>10}"
           f"{'avg':>8}{'M':>8}{'D':>8}{'cov<2.5':>9}{'obj best':>11}")
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    tbl = {}
    for nm in names:
        g = lambda k: col(rs, nm, k)                                    # noqa: E731
        tbl[nm] = {"gen": g("gen_best_ORACLE"), "sel": g("sel_best_ORACLE"),
                   "real": g("realised_ORACLE"), "avg": g("avg_ORACLE"),
                   "M": g("M"), "D": g("D"), "cov": g("cov_2.5_ORACLE"),
                   "obj": g("obj_best"), "ev": g("evals")}
        t = tbl[nm]
        print(f"  {nm:<14}{t['ev'].mean():7.0f}{t['gen'].mean():10.3f}"
              f"{t['sel'].mean():10.3f}{t['real'].mean():10.3f}{t['avg'].mean():8.3f}"
              f"{t['M'].mean():8.3f}{t['D'].mean():8.3f}{t['cov'].mean():9.1f}"
              f"{t['obj'].mean():11.1f}")

    qb = min(AR.QUANTUM, key=lambda n: tbl[n]["real"].mean())
    cb = min(AR.CLASSICAL, key=lambda n: tbl[n]["real"].mean())
    qg = min(AR.QUANTUM, key=lambda n: tbl[n]["gen"].mean())
    cg = min(AR.CLASSICAL, key=lambda n: tbl[n]["gen"].mean())
    print(f"\n  best quantum arm by realised: {qb}     best classical: {cb}")
    print(f"  best quantum arm by generation: {qg}    best classical: {cg}")

    out = {"n": len(rs), "best_quantum_realised": qb, "best_classical_realised": cb,
           "best_quantum_gen": qg, "best_classical_gen": cg,
           "means": {nm: {k: float(np.mean(v)) for k, v in tbl[nm].items()}
                     for nm in names}, "tests": {}}

    print("\n=== P1 PRIMARY -- REALISED (the frozen terminal operator).  "
          "Negative = quantum better. ===")
    p1 = []
    for other, lab in ((cb, f"vs best classical {cb}"),
                       ("q_untrained", "vs q_untrained (MANDATORY best-of-N)"),
                       ("c_marg", "vs c_marg (zero-information marginals)"),
                       ("c_helix", "vs c_helix (zero-information helix)"),
                       ("c_prod0.25", "vs c_prod0.25 (CNOTs deleted)"),
                       ("c_chain0.25", "vs c_chain (classical learned proposal)"),
                       ("pool500", "vs pool500 (the incumbent generator)")):
        v = verdict(tbl[qb]["real"], tbl[other]["real"], folds, f"{qb} {lab}")
        p1.append(v); print(_line(v))
    out["tests"]["P1"] = p1

    print("\n=== P3 -- GENERATION CEILING (pool best).  Negative = quantum better. ===")
    p3 = []
    for other, lab in ((cg, f"vs best classical {cg}"),
                       ("q_untrained", "vs q_untrained"),
                       ("c_marg", "vs c_marg (zero-information)"),
                       ("c_helix", "vs c_helix (zero-information)"),
                       ("c_prod0.25", "vs c_prod0.25 (CNOTs deleted)"),
                       ("c_chain0.25", "vs c_chain"),
                       ("pool500", "vs pool500 (the incumbent generator)")):
        v = verdict(tbl[qg]["gen"], tbl[other]["gen"], folds, f"{qg} {lab}")
        p3.append(v); print(_line(v))
    out["tests"]["P3"] = p3

    print("\n=== P3b -- SELECTION CEILING (best of the frozen native-free top-75) ===")
    p3b = []
    for other in (cg, "q_untrained", "c_marg", "pool500"):
        v = verdict(tbl[qg]["sel"], tbl[other]["sel"], folds, f"{qg} sel vs {other}")
        p3b.append(v); print(_line(v))
    out["tests"]["P3b"] = p3b

    print("\n=== P2 -- the (M, D) PARETO plane: additive epsilon-dominance, per target. ===")
    print("    eps > 0 = NO classical point dominates.  The leave-one-out row for each")
    print("    CLASSICAL arm is the null Sprint 17 records: a classical arm scored against")
    print("    its own siblings looked exactly like the VQE.")
    cls = list(AR.CLASSICAL) + list(AR.NULLS) + ["pool500"]
    eps = {}
    for nm in names:
        e = []
        for k, r in enumerate(rs):
            p = {"M": r["arms"][nm]["M"], "D": r["arms"][nm]["D"]}
            cset = [{"M": r["arms"][c]["M"], "D": r["arms"][c]["D"]}
                    for c in cls if c != nm]
            e.append(QL.eps_dominance(p, cset))
        eps[nm] = np.asarray(e, float)
    print(f"\n  {'arm':<14}{'mean eps':>10}{'median':>9}{'CI95':>22}{'eps>0':>8}")
    for nm in names:
        p = I.paired(eps[nm], np.zeros_like(eps[nm]))
        print(f"  {nm:<14}{eps[nm].mean():10.3f}{np.median(eps[nm]):9.3f}"
              f"   [{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}]"
              f"{int((eps[nm] > 0).sum()):8d}/{len(rs)}")
    out["eps"] = {nm: float(eps[nm].mean()) for nm in names}
    out["eps_pos"] = {nm: int((eps[nm] > 0).sum()) for nm in names}

    print("\n=== The budget-trap control, stated on its own: does TRAINING help? ===")
    for nm in AR.QUANTUM:
        for key, lab in (("real", "realised"), ("gen", "generation"), ("D", "diversity")):
            v = verdict(tbl[nm][key], tbl["q_untrained"][key], folds,
                        f"{nm} - q_untrained  ({lab})")
            print(_line(v))
    print("\n=== What the CNOTs buy: latent nearest-neighbour mutual information ===")
    for nm in ("q_a0.25", "c_prod0.25", "q_untrained"):
        mi = np.array([r["arms"][nm].get("nn_mi_bits_mean", np.nan) for r in rs], float)
        ent = np.array([r["arms"][nm].get("entropy_bits", np.nan) for r in rs], float)
        ess = np.array([r["arms"][nm].get("ess", np.nan) for r in rs], float)
        print(f"  {nm:<14} nn_MI {np.nanmean(mi):+.4f} bits   entropy "
              f"{np.nanmean(ent):6.3f} bits   ESS {np.nanmean(ess):8.1f}")

    with open(os.path.join(L.RESULTS, f"qb_report_{tag}.json"), "w") as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return out


def law(tag="main"):
    """Does the LEDGER's terminal-operator law hold in this lane's own data?

        d_out = a * d_set_mean + b * d_set_best + c
    """
    rs = rows(tag)
    names = list(AR.arms().keys())
    y, xm, xb = [], [], []
    for r in rs:
        for nm in names:
            a = r["arms"][nm]
            y.append(a["realised_ORACLE"]); xm.append(a["sel_mean_ORACLE"])
            xb.append(a["sel_best_ORACLE"])
    X = np.column_stack([xm, xb, np.ones(len(y))])
    co, *_ = np.linalg.lstsq(X, np.asarray(y, float), rcond=None)
    pred = X @ co
    r2 = 1.0 - ((np.asarray(y) - pred) ** 2).sum() / ((np.asarray(y) - np.mean(y)) ** 2).sum()
    print(f"\n  d_out = {co[0]:.3f}*d_set_mean + {co[1]:.3f}*d_set_best + {co[2]:.3f}"
          f"   R2 = {r2:.3f}   n = {len(y)} (arm, target) cells")
    print(f"  LEDGER (a different instrument): 1.16*mean + 0.04*best, R2 0.89.")
    print(f"  => one Angstrom of GENERATION ceiling is worth {co[1]:.3f} A of realised RMSD.")
    return {"coef": co.tolist(), "r2": float(r2), "n": len(y)}


def theory():
    p1 = os.path.join(L.RESULTS, "qb_t1.json")
    if os.path.exists(p1):
        d = json.load(open(p1))
        print(f"\n=== T1 -- CVaR INDIFFERENCE, verified.  n = {len(d['rows'])} targets, "
              f"M = {d['rows'][0]['M']} continuous configurations. ===")
        print("   Three laws, all EXACTLY CVaR-optimal at the same alpha, and what they emit.\n")
        for a in ("0.05", "0.25", "0.5"):
            sp = np.array([r["alphas"][a]["cvar_spread"] for r in d["rows"]])
            print(f"  alpha = {a}   worst-case CVaR spread max(A,B,C) - min(A,B,C) over targets "
                  f"= {sp.max():.3e}")
            for nm in ("A_delta", "B_uniform_tail", "C_worst_tail"):
                g = lambda k: np.array([r["alphas"][a][nm][k] for r in d["rows"]])  # noqa
                print(f"     {nm:<16} D {g('D').mean():6.3f}   set best "
                      f"{g('set_best_ORACLE').mean():6.3f}   set mean "
                      f"{g('set_mean_ORACLE').mean():6.3f}   coord-avg "
                      f"{g('avg_ORACLE').mean():6.3f}")
            print()
    p2 = os.path.join(L.RESULTS, "qb_t2.json")
    if os.path.exists(p2):
        d = json.load(open(p2))
        print(f"=== T1b -- the RESTRICTED-family form on the actual mps2f circuit.  "
              f"n = {len(d['rows'])} targets, K = {d['rows'][0]['K']} random parameters ===\n")
        for a in ("0.05", "0.25", "0.5"):
            g = lambda k: np.array([r["alphas"][a][k] for r in d["rows"]], float)  # noqa
            bd = np.array([r["alphas"][a]["band_avg_range_ORACLE"] for r in d["rows"]])
            bc = np.array([r["alphas"][a]["band_rel_cvar_spread"] for r in d["rows"]])
            print(f"  alpha = {a}:  rho(CVaR, D) = {g('rho_cvar_D').mean():+.3f}   "
                  f"rho(CVaR, coord-avg RMSD) = {g('rho_cvar_avg_ORACLE').mean():+.3f}   "
                  f"rho(CVaR, gen best) = {g('rho_cvar_gen_ORACLE').mean():+.3f}")
            print(f"            inside the best-12.5% CVaR band: CVaR varies "
                  f"{100*bc.mean():.1f}% while the emitted coordinate average spans "
                  f"{(bd[:, 1]-bd[:, 0]).mean():.3f} A")


def gauge(tag="gauge_PARTIAL_n20"):
    d = L.load(tag)
    rs = sorted([v for k, v in d.items() if not k.startswith("_")],
                key=lambda x: x["pdb"])
    if not rs:
        print("no gauge artefact"); return
    pcs, ids, orb = [], [], []
    for r in rs:
        v = np.array([o["realised"] for o in r["orbit"]], float)
        x = float(r["identity"]["realised"])
        lt = float((v < x).sum()); eq = float((v == x).sum())
        pcs.append((lt + 0.5 * eq) / len(v))            # MID-RANK (the Sprint-18 tie trap)
        ids.append(x); orb.append(v.mean())
    pcs = np.asarray(pcs); ids = np.asarray(ids); orb = np.asarray(orb)
    print(f"\n=== Z2^n GAUGE ORBIT, n = {len(rs)} targets, "
          f"{len(rs[0]['orbit'])} random relabellings each ===")
    print(f"  identity labelling realised   {ids.mean():.3f}")
    print(f"  gauge-orbit mean realised     {orb.mean():.3f}")
    p = I.paired(ids, orb)
    print(f"  identity - orbit mean         {p['mean_diff']:+.3f}  med "
          f"{p['median_diff']:+.3f}  [{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}]  "
          f"W/L {p['n_better']}/{p['n_worse']}")
    print(f"  identity's MID-RANK percentile in its own orbit: mean {pcs.mean():.3f}  "
          f"median {np.median(pcs):.3f}  below 0.05 on {int((pcs < 0.05).sum())}/{len(pcs)}")
    print("  FIRES (result is encoding-dependent) iff the mean percentile < 0.05.")
    return {"identity": float(ids.mean()), "orbit": float(orb.mean()),
            "pct_mean": float(pcs.mean())}


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else "main"
    {"main": main, "theory": theory, "gauge": gauge, "law": law}[m]()
