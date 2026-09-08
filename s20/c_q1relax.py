"""s20/c_q1relax.py -- Q1, RE-SPECIFIED: the two potentials compared through a SHARED operator.

WHY THIS MODULE EXISTS.  Mid-sprint the coordinator's `s20/LEDGER.md` L6 audit established that
on the production k=8 register the AMBER Hamiltonian is not `E_amber` at all -- it is
`E_amber o Relax_50`, and **the relaxation is constitutive**: the as-built ideal-geometry
structures are in hard steric overlap, unrelaxed AMBER is not finite on 42% of that register, and
the 50-iteration cap binds on 192/192 calls.  So every prior "Legacy and AMBER rank differently"
statement, MINE INCLUDED (`s20/results/c_q1_report.txt` section 1), compared `E_legacy` against
`E_amber o Relax_50` -- two functions that differ by an operator as well as by a potential.

**IT REPRODUCES ON MY DOMAIN, IN A WORSE FORM, AND I MEASURED IT BEFORE WRITING THIS.**  On the
9450 shipped top-75 ideal-geometry rebuilds the unrelaxed AMBER single point is numerically
finite everywhere -- and

    median 16062 kcal/mol,  p99 3.1e13,  max 5.5e23,
    53.5% above 1e4 kcal/mol,  31.1% above 1e6

against a physical range of about -1170 .. -500 for a relaxed peptide of this size.  The collapse
is therefore **not a lattice artefact**; it is a property of the ideal-geometry rebuild, and on
this domain it is *finite-but-meaningless* rather than +inf, which is more dangerous because
nothing raises an exception.  This is the project's own
`pauli-spectrum-delta-spike-artefact` in a new place: **an unconditioned AMBER energy measures its
worst clash.**

THE RE-SPECIFIED COMPARISON.  Two arms, both pre-registered here before the run:

  (b) PRIMARY -- give Legacy the SAME preprocessing.  Score BOTH potentials at the
      AMBER-relaxed coordinates, so the relaxation is a shared operator and CANCELS:

          rho( E_legacy o Relax_s ,  E_amber o Relax_s )        s in {1, 50}

      This is the only genuine "change only H" available, and it asks the question that matters:
      what do the two POTENTIALS disagree about when both see the same geometry?

  (a) DIAGNOSTIC -- restrict to the subset where the UNRELAXED AMBER energy is physical, and run
      `E_legacy`, `E_amber o Relax_1`, `E_amber o Relax_50` on that common domain.  The
      `Relax_50 - Relax_1` contrast then isolates **the relaxation's own effect, separated from
      the physics**.  The restriction is part of the operator and is reported as such, with
      `n_collapsed` per arm.

COLLAPSE THRESHOLD, DECLARED BEFORE THE RUN.  `collapsed := E_amber_unrelaxed > 1e4 kcal/mol` --
three orders of magnitude above the physical range, so it cannot be a borderline call.  The census
is also printed at 1e3 and 1e6 so the choice is visible rather than load-bearing.

THE RELAXATION.  `s19.agentC_pareto._run_ref`, which gate GC2 certifies is bit-identical to
`core.amber._run` and to the deployed `refine_coords` (0.00e+00 A / 0.00e+00 kcal/mol), at the
deployed restraint constant `core.amber.K_MODERATE = 10.0` -- the same constant
`core.pipeline.Config.amber_k` uses -- with `steps` in {1, 50}.  `steps` is a CAP here, not a
convergence criterion, so **every arm reports how many calls actually hit it** (s19 Z6).

    python -m s20.c_q1relax --smoke     # 2 targets
    python -m s20.c_q1relax             # the declared subset
    python -m s20.c_q1relax --report
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
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                     # noqa: E402
from s14.avgspace import top75_windows              # noqa: E402
from s15 import seed as SD                          # noqa: E402
from s16 import energy_lib as EL                    # noqa: E402
from s18 import phys_lib as PL                      # noqa: E402
from s19 import agentC_pareto as CP                 # noqa: E402
from s20 import c_q1 as Q1                          # noqa: E402

STEPS = (1, 50)
K_RELAX = 10.0                 # core.amber.K_MODERATE == core.pipeline.Config.amber_k
TOL = 1.0
COLLAPSE = 1e4                 # kcal/mol, declared before the run
N_TARGETS = 40
ATOMS = ("N", "CA", "C", "O", "CB")


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


def subset(n=N_TARGETS):
    tg = I.targets()
    rng = SD.stable_rng("s20C_q1relax", "subset")
    idx = np.sort(rng.permutation(len(tg))[:n])
    return [tg[int(i)] for i in idx]


def legacy_at(seq, bb_batch, w):
    """Genuine Legacy at ARBITRARY coordinates -- torsions extracted from those same
    coordinates, so nothing about the structure is taken from anywhere else."""
    from core import energy as et
    from core import geometry as geo
    PHI, PSI = geo.extract_torsions_batch(bb_batch["N"], bb_batch["CA"], bb_batch["C"])
    comp = et.components_batch(seq, bb_batch, PHI, PSI)
    M = np.column_stack([np.asarray(comp[t], float) for t in EL.LEG_TERMS])
    return M @ w, PHI, PSI


def run_target(t):
    from core import geometry as geo
    import torsion_lib2 as tl2
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    W, PHI, PSI, u = top75_windows(pdb)
    PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    nat = np.asarray(u["nat_ca"], float)
    K = len(PHI)
    w = EL.legacy_weight_vector()

    #: the as-built candidates -- the identical ideal-geometry rebuilds Q1 scored.
    bb0 = geo.build_backbone_batch(PHI, PSI)
    e_leg0, _p, _q = legacy_at(seq, bb0, w)
    d_reb0 = np.asarray(I.kabsch_rmsd_batch(np.asarray(bb0["CA"], float), nat), float)

    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    H = CP.builder(seq, rep, CP.SETS["bb"])

    out = {"pdb": pdb, "n": n, "fold": fold, "K": int(K)}
    rel = {}
    for s in STEPS:
        Eam = np.empty(K); walls = np.empty(K)
        bbr = None
        for b in range(K):
            c = {a: np.asarray(bb0[a][b], float) for a in bb0}
            pos = H._assemble(H._heavy_positions(c, chi1=None))
            r = CP._run_ref(H, pos, pos, K_RELAX, steps=int(s), tolerance=TOL)
            Eam[b] = float(r["energy"])
            walls[b] = float(r["wall"])
            #: `core.amber._split` returns whichever backbone atoms the topology carries;
            #: glycine has no CB, so the key set is discovered rather than assumed.  Legacy's
            #: CB-pair terms are simply absent for a glycine, exactly as in the scalar path.
            if bbr is None:
                bbr = {a: np.empty((K, n, 3)) for a in r["backbone"]}
            for a in bbr:
                bbr[a][b] = np.asarray(r["backbone"][a], float)
        e_leg_r, _p, _q = legacy_at(seq, bbr, w)
        rel[s] = {"e_amber": Eam, "e_legacy": e_leg_r,
                  "d_reb": np.asarray(I.kabsch_rmsd_batch(bbr["CA"], nat), float),
                  "wall": walls}
    CP.drop_builders()

    #: the UNRELAXED AMBER single point (Relax_0) -- read from the s18 artefact, whose Legacy
    #: components gate GC20a already certified are these candidates' (0.00e+00).
    e_amb0 = np.asarray(Q1.load_down()[pdb]["s_amber_sp"], float) if _DOWN is None \
        else np.asarray(_DOWN[pdb]["s_amber_sp"], float)

    out["collapse"] = {f"gt_{th:g}": int((e_amb0 > th).sum()) for th in (1e3, COLLAPSE, 1e6)}
    out["n_collapsed"] = int((e_amb0 > COLLAPSE).sum())
    ok = e_amb0 <= COLLAPSE                     # the PHYSICAL subset, arm (a)'s domain

    sp = Q1._spearman
    out["rho"] = {
        #: what every previous statement actually measured
        "leg0_vs_amb0": sp(e_leg0, e_amb0),
        #: (b) THE PRIMARY -- shared operator, it cancels
        "leg50_vs_amb50": sp(rel[50]["e_legacy"], rel[50]["e_amber"]),
        "leg1_vs_amb1": sp(rel[1]["e_legacy"], rel[1]["e_amber"]),
        #: the mismatched comparisons, kept so the reader can see the operator's contribution
        "leg0_vs_amb50": sp(e_leg0, rel[50]["e_amber"]),
        "leg0_vs_amb1": sp(e_leg0, rel[1]["e_amber"]),
        #: how much the relaxation MOVES each potential's own ordering
        "amb0_vs_amb50": sp(e_amb0, rel[50]["e_amber"]),
        "amb1_vs_amb50": sp(rel[1]["e_amber"], rel[50]["e_amber"]),
        "leg0_vs_leg50": sp(e_leg0, rel[50]["e_legacy"]),
        #: ORACLE, the only thing that matters for a ranker
        "amb50_vs_d": sp(rel[50]["e_amber"], rel[50]["d_reb"]),
        "leg50_vs_d": sp(rel[50]["e_legacy"], rel[50]["d_reb"]),
        "amb0_vs_d": sp(e_amb0, d_reb0),
        "leg0_vs_d": sp(e_leg0, d_reb0),
    }
    #: (a) THE DIAGNOSTIC -- the same quantities on the PHYSICAL subset only.
    out["n_ok"] = int(ok.sum())
    if ok.sum() >= 8:
        out["rho_ok"] = {
            "leg0_vs_amb0": sp(e_leg0[ok], e_amb0[ok]),
            "leg0_vs_amb1": sp(e_leg0[ok], rel[1]["e_amber"][ok]),
            "leg0_vs_amb50": sp(e_leg0[ok], rel[50]["e_amber"][ok]),
            "leg50_vs_amb50": sp(rel[50]["e_legacy"][ok], rel[50]["e_amber"][ok]),
            "amb1_vs_amb50": sp(rel[1]["e_amber"][ok], rel[50]["e_amber"][ok]),
        }
    else:
        out["rho_ok"] = None

    #: levels and displacement, so the relaxation's size is visible and not only its rank effect
    out["levels"] = {
        "e_amb0_median": float(np.median(e_amb0)),
        "e_amb1_median": float(np.median(rel[1]["e_amber"])),
        "e_amb50_median": float(np.median(rel[50]["e_amber"])),
        "e_leg0_mean": float(e_leg0.mean()),
        "e_leg50_mean": float(rel[50]["e_legacy"].mean()),
        "d_reb0_mean": float(d_reb0.mean()),
        "d_reb1_mean": float(rel[1]["d_reb"].mean()),
        "d_reb50_mean": float(rel[50]["d_reb"].mean()),
        "wall50_mean": float(rel[50]["wall"].mean()),
        "wall50_max": float(rel[50]["wall"].max()),
        "wall1_mean": float(rel[1]["wall"].mean()),
    }
    #: THE CAP FIRING COUNT (s19 Z6).  `steps` is a cap, not a convergence criterion: a call
    #: whose energy still moves between 50 and 100 iterations HIT the cap.  Measured on a
    #: deterministic 8-candidate probe per target rather than assumed.
    probe = np.linspace(0, K - 1, min(8, K)).astype(int)
    hit = 0
    H = CP.builder(seq, rep, CP.SETS["bb"])
    for b in probe:
        c = {a: np.asarray(bb0[a][b], float) for a in bb0}
        pos = H._assemble(H._heavy_positions(c, chi1=None))
        r50 = CP._run_ref(H, pos, pos, K_RELAX, steps=50, tolerance=TOL)
        r100 = CP._run_ref(H, pos, pos, K_RELAX, steps=100, tolerance=TOL)
        if abs(float(r50["energy"]) - float(r100["energy"])) > 1e-3:
            hit += 1
    CP.drop_builders()
    out["cap_fired_50_of_probe"] = int(hit)
    out["cap_probe_n"] = int(len(probe))
    return out


_DOWN = None


def run(n=N_TARGETS, out="c_q1relax.json", verbose=True):
    global _DOWN
    _DOWN = Q1.load_down()
    tg = subset(n)
    path = os.path.join(RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    cfg = {"STEPS": list(STEPS), "K_RELAX": K_RELAX, "TOL": TOL, "COLLAPSE": COLLAPSE,
           "subset": [t["pdb"] for t in tg],
           "note": "coordinator LEDGER L6 re-specification: shared Relax operator"}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        PL.mem_ok(0.55)
        rows.append(run_target(t))
        if verbose:
            r = rows[-1]
            print(f"  {len(rows)}/{len(tg)} {t['pdb']} collapsed={r['n_collapsed']}/{r['K']} "
                  f"rho(L0,A0)={r['rho']['leg0_vs_amb0']:+.3f} "
                  f"rho(L50,A50)={r['rho']['leg50_vs_amb50']:+.3f} "
                  f"rho(A0,A50)={r['rho']['amb0_vs_amb50']:+.3f} "
                  f"cap {r['cap_fired_50_of_probe']}/{r['cap_probe_n']} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        _write(out, rows, cfg, len(tg))
    _write(out, rows, cfg, len(tg))
    print(f"DONE {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)


def _write(out, rows, cfg, n_expected):
    obj = {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg),
           "n_rows": len(rows), "n_expected": int(n_expected),
           "complete": bool(len(rows) >= int(n_expected))}
    p = os.path.join(RESULTS, out)
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, p)


def report(out="c_q1relax.json"):
    o = json.load(open(os.path.join(RESULTS, out)))
    rows = o["rows"]
    n = len(rows)
    tag = "" if o.get("complete") else f"  *** PARTIAL n={n}/{o.get('n_expected')} ***"
    folds = np.array([r["fold"] for r in rows], int)
    print("=" * 104)
    print(f"Q1d  THE TWO POTENTIALS THROUGH A SHARED OPERATOR   n = {n} targets x 75 "
          f"candidates{tag}")
    print("     coordinator LEDGER L6 re-specification.  Relax = core.amber at k = 10.0 "
          "(the DEPLOYED constant),")
    print("     steps in {1, 50}; Legacy re-scored at the SAME relaxed coordinates so the "
          "operator cancels.")
    print("=" * 104)

    print("\n0.  THE COLLAPSE CENSUS ON THIS DOMAIN  (unrelaxed AMBER single point, kcal/mol)")
    K = sum(r["K"] for r in rows)
    for th in ("gt_1000", "gt_10000", "gt_1e+06"):
        c = sum(r["collapse"].get(th, 0) for r in rows)
        print(f"    candidates above {th.replace('gt_',''):>8s} kcal/mol: {c:6d} / {K} "
              f"({100*c/K:5.1f}%)")
    nc = np.array([r["n_collapsed"] for r in rows], float)
    print(f"    DECLARED collapsed (> {o['config']['COLLAPSE']:g}): mean {nc.mean():.1f} of 75 "
          f"per target; {int((nc == 0).sum())}/{n} targets with none; "
          f"{int((nc >= 60).sum())}/{n} targets with 60+")
    cf = np.array([r["cap_fired_50_of_probe"] for r in rows], float)
    cp = np.array([r["cap_probe_n"] for r in rows], float)
    print(f"    THE 50-ITERATION CAP FIRED on {int(cf.sum())}/{int(cp.sum())} probe calls "
          f"(|dE| > 1e-3 kcal/mol between 50 and 100 iterations) -- s19 Z6: a cap that never "
          f"fires certifies nothing, and this one fires.")

    print("\n1.  LEVELS -- what the relaxation actually does (medians over candidates, "
          "then over targets)")
    for k in ("e_amb0_median", "e_amb1_median", "e_amb50_median", "e_leg0_mean",
              "e_leg50_mean", "d_reb0_mean", "d_reb1_mean", "d_reb50_mean",
              "wall1_mean", "wall50_mean", "wall50_max"):
        v = np.array([r["levels"][k] for r in rows], float)
        print(f"    {k:<18s} median over targets {np.median(v):16.4f}   mean {v.mean():16.4f}")

    print("\n2.  *** THE RE-SPECIFIED PRIMARY ***  per-target Spearman over the 75 candidates")
    print(f"{'comparison':<26s}{'mean':>9s}{'CI95':>24s}{'median':>9s}{'>0':>7s}"
          f"{'|rho|>0.8':>11s}")
    order = ["leg0_vs_amb0", "leg0_vs_amb1", "leg0_vs_amb50",
             "leg1_vs_amb1", "leg50_vs_amb50",
             "amb0_vs_amb50", "amb1_vs_amb50", "leg0_vs_leg50",
             "leg0_vs_d", "leg50_vs_d", "amb0_vs_d", "amb50_vs_d"]
    lab = {"leg0_vs_amb0": "E_leg  vs E_amb  (raw)",
           "leg0_vs_amb1": "E_leg  vs A o R1   MIXED",
           "leg0_vs_amb50": "E_leg  vs A o R50  MIXED",
           "leg1_vs_amb1": "L o R1  vs A o R1  SHARED",
           "leg50_vs_amb50": "L o R50 vs A o R50 SHARED",
           "amb0_vs_amb50": "AMBER: raw vs R50",
           "amb1_vs_amb50": "AMBER: R1 vs R50",
           "leg0_vs_leg50": "Legacy: raw vs R50",
           "leg0_vs_d": "E_leg  vs ORACLE d",
           "leg50_vs_d": "L o R50 vs ORACLE d",
           "amb0_vs_d": "E_amb  vs ORACLE d",
           "amb50_vs_d": "A o R50 vs ORACLE d"}
    for k in order:
        v = np.array([r["rho"][k] for r in rows], float)
        st = PP(v, np.zeros_like(v), folds=folds)
        print(f"{lab[k]:<26s}{v.mean():+9.4f}  [{st['ci95'][0]:+9.4f},{st['ci95'][1]:+9.4f}]"
              f"{np.median(v):+9.4f}{(v>0).mean():7.2f}{(np.abs(v)>0.8).mean():11.2f}")

    a = np.array([r["rho"]["leg50_vs_amb50"] for r in rows], float)
    b = np.array([r["rho"]["leg0_vs_amb0"] for r in rows], float)
    c = np.array([r["rho"]["leg0_vs_amb50"] for r in rows], float)
    st = PP(a, b, folds=folds)
    print(f"\n    SHARED(R50) minus RAW:   {st['mean_diff']:+.4f} "
          f"[{st['ci95'][0]:+.4f}, {st['ci95'][1]:+.4f}]  med {st['median_diff']:+.4f}  "
          f"{st['n_better']}W/{st['n_worse']}L")
    st = PP(a, c, folds=folds)
    print(f"    SHARED(R50) minus MIXED: {st['mean_diff']:+.4f} "
          f"[{st['ci95'][0]:+.4f}, {st['ci95'][1]:+.4f}]  med {st['median_diff']:+.4f}  "
          f"{st['n_better']}W/{st['n_worse']}L"
          f"    <- how much of the 'disagreement' was the OPERATOR, not the potential")

    ok = [r for r in rows if r.get("rho_ok")]
    if ok:
        fo = np.array([r["fold"] for r in ok], int)
        nk = np.array([r["n_ok"] for r in ok], float)
        print(f"\n3.  DIAGNOSTIC (a) -- restricted to the PHYSICAL subset "
              f"(unrelaxed AMBER <= {o['config']['COLLAPSE']:g} kcal/mol).")
        print(f"    THE RESTRICTION IS PART OF THE OPERATOR: mean {nk.mean():.1f} of 75 "
              f"candidates survive it, on {len(ok)}/{n} targets with >= 8 survivors.")
        print(f"{'comparison':<26s}{'mean':>9s}{'CI95':>24s}{'median':>9s}")
        for k in ("leg0_vs_amb0", "leg0_vs_amb1", "leg0_vs_amb50", "leg50_vs_amb50",
                  "amb1_vs_amb50"):
            v = np.array([r["rho_ok"][k] for r in ok], float)
            st = PP(v, np.zeros_like(v), folds=fo)
            print(f"{lab[k]:<26s}{v.mean():+9.4f}  [{st['ci95'][0]:+9.4f},"
                  f"{st['ci95'][1]:+9.4f}]{np.median(v):+9.4f}")
        v1 = np.array([r["rho_ok"]["leg0_vs_amb50"] for r in ok], float)
        v0 = np.array([r["rho_ok"]["leg0_vs_amb1"] for r in ok], float)
        st = PP(v1, v0, folds=fo)
        print(f"\n    THE RELAXATION'S OWN EFFECT, physics held fixed "
              f"(R50 - R1 on the common domain): {st['mean_diff']:+.4f} "
              f"[{st['ci95'][0]:+.4f}, {st['ci95'][1]:+.4f}]  med {st['median_diff']:+.4f}  "
              f"{st['n_better']}W/{st['n_worse']}L")
    print("=" * 104)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--n", type=int, default=N_TARGETS)
    a = ap.parse_args()
    if a.report:
        report()
    elif a.smoke:
        run(n=2, out="_SMOKE_c_q1relax.json")
    else:
        run(n=a.n)
