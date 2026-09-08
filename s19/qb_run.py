"""SPRINT 19 / AGENT B -- the sweep.  One heavy process, BLAS capped at 1.

    python -m s19.qb_run run              # all 126 targets, every arm, seed 0
    python -m s19.qb_run run 30           # a declared prefix
    python -m s19.qb_run gauge            # the pre-registered Z2^n gauge orbit
    python -m s19.qb_run verify           # P0 checks, run before any claim
"""
from __future__ import annotations

import gc
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I           # noqa: E402
from s19 import qb_lib as L               # noqa: E402
from s19 import qb_arms as AR             # noqa: E402


def run(limit=None, seed=0, tag="main"):
    tg = I.targets()
    if limit:
        tg = tg[:int(limit)]
    L.gather_full()
    have = L.load(tag)
    A = AR.arms()
    t00 = time.time()
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        if pdb in have:
            continue
        # MEMORY DISCIPLINE.  This process's own footprint is ~0.3 GB; the box is shared
        # with sibling workstreams, so a hard stop at the first dip hands the sweep to
        # whoever launched last.  Wait for the box instead, and stop only if it never
        # recovers.  Every wall-clock number in this lane is CONTENDED and is not used
        # for any cost comparison.
        fg = I.free_gb()
        for _w in range(30):
            if fg >= 0.60:
                break
            print(f"[wait] free RAM {fg:.2f} GB < 0.60, sleeping 60s", flush=True)
            time.sleep(60)
            fg = I.free_gb()
        if fg < 0.60:
            print(f"[STOP] free RAM {fg:.2f} GB < 0.60 for 30 min", flush=True)
            return
        t0 = time.time()
        tgt = L.target(pdb)
        row = {"pdb": pdb, "n": tgt["n"], "fold": tgt["fold"], "arms": {}}
        for nm, (k, kw) in A.items():
            row["arms"][nm] = AR.run_arm(tgt, nm, k, kw, seed=seed)
        L.ck(tag, pdb, row)
        have[pdb] = row
        del tgt
        gc.collect()
        r = row["arms"]
        print(f"[{c+1}/{len(tg)}] {pdb} n={row['n']} {time.time()-t0:5.1f}s "
              f"free={fg:.1f}GB  q0.25 {r['q_a0.25']['realised_ORACLE']:.3f}  "
              f"untr {r['q_untrained']['realised_ORACLE']:.3f}  "
              f"chain {r['c_chain0.25']['realised_ORACLE']:.3f}  "
              f"marg {r['c_marg']['realised_ORACLE']:.3f}  "
              f"lbfgs {r['c_lbfgs']['realised_ORACLE']:.3f}", flush=True)
    L.ck(tag, "_complete", int(len([k for k in L.load(tag) if not k.startswith("_")])))
    print(f"done {time.time()-t00:.0f}s", flush=True)


def gauge(limit=30, seed=0, n_draw=24, tag="gauge_PARTIAL_n20"):
    """THE PRE-REGISTERED Z2^n GAUGE TEST.

    `b_i <-> 1-b_i` (which basin is called 0) is the ENTIRE labelling freedom of this
    latent -- unlike Sprint 18's (S4)^n, the per-residue orbit has two elements and the
    group is Z2^n.  The objective, the family of representable continuous laws and every
    RMSD are invariant; the ANSATZ's inductive bias need not be.  The identity labelling's
    mid-rank percentile in its own orbit is the statistic (mid-rank because a strict `<`
    on a degenerate orbit manufactured a result in Sprint 18).
    """
    from s15 import seed as SD
    tg = I.targets()[:int(limit)]
    L.gather_full()
    have = L.load(tag)
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        if pdb in have:
            continue
        t0 = time.time()
        tgt = L.target(pdb)
        rng = SD.stable_rng(pdb, "gaugeflips", salt=L.SALT)
        base = AR.run_arm(tgt, "q_a0.25", "circuit",
                          dict(alpha=0.25, ansatz="mps2f", train=True), seed=seed)
        vals = []
        mu0, kap0 = tgt["mu"].copy(), tgt["kap"].copy()
        w0 = tgt["wmarg"].copy()
        for g in range(int(n_draw)):
            f = rng.integers(0, 2, tgt["n"]).astype(bool)
            tgt["mu"] = mu0.copy(); tgt["kap"] = kap0.copy(); tgt["wmarg"] = w0.copy()
            tgt["mu"][f] = mu0[f][:, ::-1]
            tgt["kap"][f] = kap0[f][:, ::-1]
            tgt["wmarg"][f] = w0[f][:, ::-1]
            r = AR.run_arm(tgt, "q_a0.25", "circuit",
                           dict(alpha=0.25, ansatz="mps2f", train=True), seed=seed)
            vals.append({"realised": r["realised_ORACLE"], "gen": r["gen_best_ORACLE"],
                         "M": r["M"], "D": r["D"]})
        tgt["mu"], tgt["kap"], tgt["wmarg"] = mu0, kap0, w0
        row = {"pdb": pdb, "n": tgt["n"], "fold": tgt["fold"],
               "identity": {"realised": base["realised_ORACLE"],
                            "gen": base["gen_best_ORACLE"],
                            "M": base["M"], "D": base["D"]},
               "orbit": vals}
        L.ck(tag, pdb, row)
        have[pdb] = row
        print(f"[{c+1}/{len(tg)}] {pdb} {time.time()-t0:.0f}s  identity "
              f"{base['realised_ORACLE']:.3f}  orbit mean "
              f"{np.mean([v['realised'] for v in vals]):.3f}", flush=True)
        del tgt
        gc.collect()
    L.ck(tag, "_complete", int(len([k for k in L.load(tag) if not k.startswith("_")])))


def verify():
    """P0 -- the checks that must pass before any number in this lane is read."""
    from s15 import align_lib as A
    from s17 import q_lib as QL
    out = {}
    print("== P0.1  the (M, D) identity  readout^2 = M^2 - D^2 ==")
    L.gather_full()
    tgt = L.target(I.targets()[0]["pdb"])
    rng = L.SD.stable_rng("verify", salt=L.SALT)
    worst = 0.0
    from core import project as pj
    for m in (3, 17, 75, 200):
        b = rng.integers(0, 2, (m, tgt["n"]))
        phi, psi = L.draw_from_basins(b, tgt["mu"], tgt["kap"], rng)
        W = np.asarray(pj.build_ca_exact(phi, psi), float)
        r = QL.md_plane(W, tgt["nat"])
        worst = max(worst, abs(r["identity_resid"]))
    out["md_identity_max_abs_A2"] = worst
    out["md_identity_PASS"] = bool(worst < 1e-9)
    print(f"   max |resid| = {worst:.3e} A^2   PASS={out['md_identity_PASS']}")

    print("== P0.2  qb_lib.Obj reproduces s15/align_lib.fit's functional EXACTLY ==")
    o = L.new_obj(tgt)
    b = rng.integers(0, 2, (4, tgt["n"]))
    phi, psi = L.draw_from_basins(b, tgt["mu"], tgt["kap"], rng)
    # `A.fit` returns (phi, psi, f) at ITS optimum; evaluate MY functional at THAT point.
    # Equality there is an exact identity check of the functional, not of the optimiser.
    ph, ps, theirs = [], [], []
    for q in range(4):
        p_, q_, f = A.fit(tgt["dhat"], tgt["sd"], tgt["i"], tgt["j"], phi[q], psi[q])
        ph.append(p_); ps.append(q_); theirs.append(f)
    mine, _ = o.raw(np.array(ph), np.array(ps))
    d = float(np.max(np.abs(mine - np.array(theirs))))
    out["obj_vs_alignlib_max_abs"] = d
    out["obj_PASS"] = bool(d < 1e-8 * max(1.0, float(np.max(np.abs(mine)))))
    print(f"   max |E_mine - align_lib.fit(.).f| at align_lib's own optimum = {d:.3e}   PASS={out['obj_PASS']}")

    print("== P0.3  the sampler's law is CONTINUOUS and FULL-SUPPORT ==")
    b = rng.integers(0, 2, (4096, tgt["n"]))
    phi, psi = L.draw_from_basins(b, tgt["mu"], tgt["kap"], rng)
    u = np.unique(np.round(phi, 9))
    out["distinct_phi_values"] = int(u.size)
    out["kappa_max"] = float(tgt["kap"].max())
    out["kappa_min"] = float(tgt["kap"].min())
    out["full_support_PASS"] = bool(np.isfinite(tgt["kap"]).all()
                                    and tgt["kap"].max() <= 40.0 + 1e-9)
    print(f"   distinct phi values in 4096x{tgt['n']} draws = {u.size} "
          f"(a lattice would give <= {2*tgt['n']}); kappa in "
          f"[{tgt['kap'].min():.2f}, {tgt['kap'].max():.2f}]  "
          f"PASS={out['full_support_PASS']}")

    print("== P0.4  the budget is a HARD cap and every arm gets the same one ==")
    A2 = AR.arms()
    ev = {}
    for nm, (k, kw) in A2.items():
        oo = L.new_obj(tgt)
        if k == "pool":
            ev[nm] = 0
            continue
        getattr(AR, {"circuit": "run_circuit", "chain": "run_chain", "cem": "run_cem",
                     "metro": "run_metro", "lbfgs": "run_lbfgs", "marg": "run_marg",
                     "helix": "run_helix"}[k])(tgt, oo, seed=0, **kw)
        ev[nm] = int(oo.used)
    out["evals"] = ev
    bad = {k: v for k, v in ev.items() if v not in (0,) and not (L.BUDGET <= v <= L.BUDGET + 8)}
    out["budget_PASS"] = not bad
    print(f"   {ev}\n   PASS={out['budget_PASS']}  offenders={bad}")

    print("== P0.5  no hash(), no native quantity in any decision ==")
    import subprocess
    r = subprocess.run(["grep", "-n", "hash(", "s19/qb_lib.py", "s19/qb_arms.py",
                        "s19/qb_run.py", "s19/qb_theory.py", "s19/qb_report.py"],
                       capture_output=True, text=True, cwd=ROOT)
    out["hash_hits"] = r.stdout.strip()
    print(f"   grep 'hash(' -> {out['hash_hits'] or 'NONE'}")
    import json
    with open(os.path.join(L.RESULTS, "qb_verify.json"), "w") as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o2: o2.tolist() if hasattr(o2, "tolist") else str(o2))
    return out


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else "run"
    rest = sys.argv[2:]
    if m == "run":
        run(limit=int(rest[0]) if rest else None)
    elif m == "gauge":
        gauge(limit=int(rest[0]) if rest else 30)
    elif m == "verify":
        verify()
