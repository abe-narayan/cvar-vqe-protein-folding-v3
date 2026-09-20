#!/usr/bin/env python
"""s29/s29_D_m6.py -- lane D: THE FIXED-PROFILE CONTROL (pre-registered in `s29/PREREG_S29_D_m6.md`).

Lane T's S29-L15 (M6): the deployed CVaR is exactly flat on 437 of 511 simplex directions, the
entropy term makes the optimum uniform above the VaR and exponentially enhanced below it, and the
energies are the same standardised rank ladder on every target -- so the whole quantum stage is a
TARGET-INDEPENDENT weight profile over RANKS, and at the endpoint it reduces to one number: where
the prefix cuts. This script builds that control and runs it against the deployed CVaR-VQE arm.

  profile   compute p*(alpha, T) ONCE on the standardised rank ladder, with no target, no circuit
            and no optimiser (Sion's minimax form of the project's own free energy, verified
            against `core.quantum.free_energy` at the same alpha, T), and report its prefix m*.
  cloud     the point-cloud contrast, from `s27/results/vqe_rows.jsonl` alone (zero compute):
            per-target |RMSD(top-m* prefix) - RMSD(deployed arm)|, the agreement counts at the
            built-chain input floor and at 1e-6, and `ST.compare` for every m* in the sweep.
  chain     the same on the BUILT CHAIN (the reporting basis), both sides projected in this job
            from clouds built by the same code path.

    python s29/s29_D_m6.py profile
    python s29/s29_D_m6.py cloud [--seed 0]
    python s29/s29_D_m6.py chain [--limit N] [--seed 0]
    python s29/s29_D_m6.py analyse_chain

Every RMSD is the endpoint (ORACLE by definition); no native quantity chooses any parameter here,
and the control contains no per-target computation at all.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s29 import s29_D_cost_audit as M      # noqa: E402

RESULTS = os.path.join(HERE, "results")
VQE_ROWS = os.path.join(ROOT, "s27", "results", "vqe_rows.jsonl")
ALPHA, TEMP, LAYERS, ITERS = 0.18, 0.5, 3, 80          # the S27 deployed cell
K, NQ = 500, 9
FLOOR = 0.006          # the built chain's mean input floor (S28-L18 / RETRACTIONS_S28 R2)
FLOOR_TAIL = 0.02      # its per-target tail (12 of 126 above this)
SWEEP = (70, 71, 74, 75, 80)


# ============================================================ the profile, computed once
def rank_ladder(k=K, n_qubits=NQ):
    """The deployed energies: zrank of the score over the k candidates, padded to 2^n.

    `s24.d_harness.arm_vqe` encodes `E = zrank(score)` through `s22.qcand_lib.Encoding`; the
    SCORE's values do not enter, only its ORDER, so the ladder is the same on every target up to
    ties (S25 L17: 1.18% of range).  Padding follows `Encoding` (the pad states sit at the top).
    """
    from s22 import qcand_lib as QC
    from s27 import run_pool as RP
    E = QC.Encoding(RP.zr(np.arange(k, dtype=float))).E
    return np.asarray(E, float)


def fixed_profile(E, alpha=ALPHA, T=TEMP):
    """p*(alpha, T): the exact minimiser of F(p) = CVaR_alpha(E; p) - T H(p) over the simplex.

    Sion's minimax form (lane T, S29-L15): p*(x) proportional to exp((t* - E_x)_+ / (alpha T)),
    with t* the maximiser of t - T log sum_x exp((t - E_x)_+/(alpha T)).  Computed here by a
    one-dimensional golden-section search on t and verified against the project's own
    `core.quantum.cvar_exact` free energy, so the form is checked rather than trusted.
    """
    E = np.asarray(E, float)

    def p_of(t):
        z = np.maximum(t - E, 0.0) / (alpha * T)
        z = z - z.max()
        p = np.exp(z)
        return p / p.sum()

    def dual(t):
        z = np.maximum(t - E, 0.0) / (alpha * T)
        mx = z.max()
        return float(t - T * (mx + np.log(np.exp(z - mx).sum())))
    lo, hi = float(E.min() - 5.0), float(E.max() + 5.0)
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - gr * (b - a), a + gr * (b - a)
    for _ in range(200):
        if dual(c) < dual(d):
            a = c
        else:
            b = d
        c, d = b - gr * (b - a), a + gr * (b - a)
    t_star = 0.5 * (a + b)
    p = p_of(t_star)
    return p, t_star


def profile_report(alpha=ALPHA, T=TEMP):
    from core import quantum as Q
    from s22 import qcand_lib as QC
    E = rank_ladder()
    p, t_star = fixed_profile(E, alpha, T)
    F_star = float(Q.cvar_from_probs(E, p, alpha)[0] - T * (-(p[p > 0] * np.log(p[p > 0])).sum()))
    face = QC.exact_face(E, p, alpha)
    m_star = int(face["tail_size"])
    # controls: the uniform state, and a numerical descent from uniform (no circuit) as a check
    u = np.ones_like(E) / len(E)
    F_unif = float(Q.cvar_from_probs(E, u, alpha)[0] - T * np.log(len(E)))
    m_unif = int(QC.exact_face(E, u, alpha)["tail_size"])
    # a mirror-descent refinement of p*, to confirm the closed form is the minimiser
    q = u.copy()
    for _ in range(4000):
        v, quant, dp = Q.cvar_exact(E, q, alpha)
        g = dp / alpha + T * (np.log(np.maximum(q, 1e-300)) + 1.0)
        q = q * np.exp(-0.05 * (g - g.mean()))
        q = np.maximum(q, 1e-300); q = q / q.sum()
    F_md = float(Q.cvar_from_probs(E, q, alpha)[0] - T * (-(q[q > 0] * np.log(q[q > 0])).sum()))
    out = dict(check="the fixed profile p*(alpha, T), computed once with no target",
               alpha=alpha, T=T, k=K, n_qubits=NQ, D=int(len(E)),
               t_star=float(t_star), F_star=F_star, m_star=m_star,
               p_max=float(p.max()), p_min=float(p.min()),
               entropy_bits=float(-(p[p > 0] * np.log2(p[p > 0])).sum()),
               participation_ratio=float(1.0 / (p ** 2).sum()),
               enhancement_across_prefix=float(p[np.argsort(E)][0] / p[np.argsort(E)][max(m_star - 1, 0)]),
               F_uniform=F_unif, m_uniform=m_unif,
               F_mirror_descent=F_md, m_mirror_descent=int(QC.exact_face(E, q, alpha)["tail_size"]),
               closed_form_minus_descent=float(F_star - F_md),
               note=("the closed form is verified against a 4,000-step mirror descent on the project's own "
                     "free energy; F_star below F_uniform confirms the optimum is not the uniform state"))
    print(json.dumps({k: v for k, v in out.items() if k != "note"}, indent=1))
    return out


# ============================================================ the point-cloud contrast (zero compute)
def vqe_rows(seed=0, config="DIS"):
    rows = [json.loads(l) for l in open(VQE_ROWS, encoding="utf-8") if l.strip()]
    r = [x for x in rows if x["config"] == config and int(x["seed"]) == int(seed)]
    r.sort(key=lambda x: x["pdb"])
    return r


def cloud(seed=0, m_star=None):
    """The point-cloud contrast from the stored rows: the deployed arm against the fixed prefix.

    `rmsd_topm` in the rows is the classical top-m prefix at the arm's OWN realised m (the
    set-equality theorem, S28-L21); `rmsd_top75` is the fixed 75-prefix, which is production.
    For any other m* the prefix average has to be built, which `--recompute` does; the registered
    primary and the 75 cell need no computation at all.
    """
    rows = vqe_rows(seed)
    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    m = np.array([int(r["m"]) for r in rows])
    rv = np.array([float(r["rmsd_vqe"]) for r in rows])
    rt = np.array([float(r["rmsd_topm"]) for r in rows])
    r75 = np.array([float(r["rmsd_top75"]) for r in rows])
    out = dict(check="fixed-profile control vs the deployed CVaR-VQE arm (POINT CLOUD)", seed=int(seed),
               n=len(rows), m_star_profile=m_star, alpha=ALPHA, T=TEMP,
               m_stats=dict(mean=float(m.mean()), sd=float(m.std(ddof=1)), min=int(m.min()), max=int(m.max()),
                            eq_counts={str(v): int((m == v).sum()) for v in sorted(set(SWEEP) | ({m_star} if m_star else set()))}),
               set_equality=dict(note="the tail set IS the top-m prefix (S28-L21, 1.1e-13); so the sets are "
                                      "equal exactly when m(t) == m*, and |m(t) - m*| is the set distance",
                                 mean_abs_dm={}, exact_equal={}, jaccard_mean={}),
               vqe_mean=float(rv.mean()), topm_mean=float(rt.mean()), top75_mean=float(r75.mean()),
               vqe_vs_topm_max_abs=float(np.abs(rv - rt).max()),
               cells={})
    for ms in sorted(set(SWEEP) | ({int(m_star)} if m_star else set())):
        if ms == 75:
            arm = r75
        else:
            arm = None
        out["set_equality"]["mean_abs_dm"][str(ms)] = float(np.abs(m - ms).mean())
        out["set_equality"]["exact_equal"][str(ms)] = int((m == ms).sum())
        out["set_equality"]["jaccard_mean"][str(ms)] = float(np.mean(np.minimum(m, ms) / np.maximum(m, ms)))
        if arm is None:
            continue
        d = np.abs(arm - rv)
        c = ST.compare(arm, rv, folds, names=pdbs, label=f"fixed prefix m*={ms} minus the deployed VQE arm (point cloud, seed {seed})",
                       seed_parts=("s29Dm6",))
        out["cells"][str(ms)] = dict(mean_control=float(arm.mean()), mean_deployed=float(rv.mean()),
                                     n_within_floor=int((d <= FLOOR).sum()), n_within_floor_tail=int((d <= FLOOR_TAIL).sum()),
                                     n_within_1e6=int((d <= 1e-6).sum()), max_abs_diff=float(d.max()),
                                     median_abs_diff=float(np.median(d)),
                                     contrast={k: v for k, v in c.items() if k != "concentration"}, fmt=ST.fmt(c))
    return out


def cloud_recompute(seed=0, m_list=SWEEP, m_star=None):
    """Build the top-m* prefix average for every m* in the sweep, per target (a few seconds each)."""
    from s24 import d_harness as H
    from s27 import run_pool as RP
    from s25 import phys_lib as P
    rows = {r["pdb"]: r for r in vqe_rows(seed)}
    ms = sorted(set(m_list) | ({int(m_star)} if m_star else set()))
    out = {}
    for q, pdb in enumerate(P.targets()):
        cand, ch, _ = RP.channels_for(pdb)
        key = RP.rng_for(pdb, "tiekey").random(cand.k)
        order = np.lexsort((key, np.asarray(ch["DIS"], float)))
        rec = {}
        for mm in ms:
            C, _ = H.readout_uniform(cand, np.sort(order[:mm]))
            rec[str(mm)] = float(I.ca_rmsd(C, cand.nat_ca))                  # ORACLE (the endpoint)
        mdep = int(rows[pdb]["m"])
        Cd, _ = H.readout_uniform(cand, np.sort(order[:mdep]))
        rec["deployed"] = float(I.ca_rmsd(Cd, cand.nat_ca))
        rec["m_deployed"] = mdep
        rec["stored_vqe"] = float(rows[pdb]["rmsd_vqe"])
        out[pdb] = rec
        if (q + 1) % 30 == 0:
            print(f"  [{q+1}/126]", flush=True)
    return out


# ============================================================ the built chain
def chain(seed=0, m_list=SWEEP, m_star=None, limit=0):
    """Both sides projected in THIS job from clouds built by the same code path (S28-L43)."""
    from s24 import d_harness as H
    from s27 import run_pool as RP
    from s25 import phys_lib as P
    rows = {r["pdb"]: r for r in vqe_rows(seed)}
    ms = sorted(set(m_list) | ({int(m_star)} if m_star else set()))
    path = os.path.join(RESULTS, f"s29_D_m6_chain_rows_seed{seed}.jsonl")
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    done.add(json.loads(line)["pdb"])
    pdbs = P.targets()[:limit] if limit else P.targets()
    t0 = time.time()
    for q, pdb in enumerate(pdbs):
        if pdb in done:
            continue
        t1 = time.time()
        cand, ch, _ = RP.channels_for(pdb)
        key = RP.rng_for(pdb, "tiekey").random(cand.k)
        order = np.lexsort((key, np.asarray(ch["DIS"], float)))
        rec = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), seed=int(seed),
                   m_deployed=int(rows[pdb]["m"]), arms={})
        for tag, mm in [("deployed", int(rows[pdb]["m"]))] + [(f"m{mm}", mm) for mm in ms]:
            C, _ = H.readout_uniform(cand, np.sort(order[:mm]))
            pr = I.project(np.asarray(C, float), cand.seq, cand.fold)
            ca = np.asarray(pr["ca"], float)
            rec["arms"][tag] = dict(m=int(mm), rmsd_cloud=float(I.ca_rmsd(C, cand.nat_ca)),
                                    rmsd_chain=float(I.ca_rmsd(ca, cand.nat_ca)))      # ORACLE
        rec["secs"] = time.time() - t1
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
        a = rec["arms"]
        print(f"  [{q+1}/{len(pdbs)}] {pdb} m_dep {rec['m_deployed']:3d}  chain deployed {a['deployed']['rmsd_chain']:.4f} "
              f"m75 {a['m75']['rmsd_chain']:.4f}  diff {a['m75']['rmsd_chain']-a['deployed']['rmsd_chain']:+.4f}  "
              f"{rec['secs']:.1f}s (elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", path)
    return path


def analyse_chain(seed=0):
    path = os.path.join(RESULTS, f"s29_D_m6_chain_rows_seed{seed}.jsonl")
    rows = [json.loads(l) for l in open(path, encoding="utf-8") if l.strip()]
    rows.sort(key=lambda r: r["pdb"])
    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    dep = np.array([r["arms"]["deployed"]["rmsd_chain"] for r in rows])
    dep_c = np.array([r["arms"]["deployed"]["rmsd_cloud"] for r in rows])
    out = dict(check="fixed-profile control vs the deployed CVaR-VQE arm (BUILT CHAIN, the reporting basis)",
               seed=int(seed), n=len(rows), deployed_mean_chain=float(dep.mean()), deployed_mean_cloud=float(dep_c.mean()),
               m_deployed_mean=float(np.mean([r["m_deployed"] for r in rows])), cells={})
    tags = [t for t in rows[0]["arms"] if t != "deployed"]
    for tag in tags:
        a = np.array([r["arms"][tag]["rmsd_chain"] for r in rows])
        ac = np.array([r["arms"][tag]["rmsd_cloud"] for r in rows])
        d = np.abs(a - dep)
        c = ST.compare(a, dep, folds, names=pdbs, label=f"fixed prefix {tag} minus the deployed VQE arm (BUILT CHAIN, seed {seed})",
                       seed_parts=("s29Dm6",))
        cc = ST.compare(ac, dep_c, folds, names=pdbs, label=f"fixed prefix {tag} minus the deployed VQE arm (point cloud, seed {seed})",
                        seed_parts=("s29Dm6",))
        out["cells"][tag] = dict(m=int(rows[0]["arms"][tag]["m"]), mean_chain=float(a.mean()), mean_cloud=float(ac.mean()),
                                 n_within_floor=int((d <= FLOOR).sum()), n_within_floor_tail=int((d <= FLOOR_TAIL).sum()),
                                 n_within_1e6=int((d <= 1e-6).sum()), max_abs_diff=float(d.max()), median_abs_diff=float(np.median(d)),
                                 contrast_chain={k: v for k, v in c.items() if k != "concentration"},
                                 contrast_cloud={k: v for k, v in cc.items() if k != "concentration"},
                                 fmt_chain=ST.fmt(c), fmt_cloud=ST.fmt(cc))
        print(f"  {tag:8s} chain mean {a.mean():.4f} vs deployed {dep.mean():.4f}   effect {c['effect']:+.4f} "
              f"{c['effect_over_mde']:+.2f}x fold [{c['ci95_fold'][0]:+.4f},{c['ci95_fold'][1]:+.4f}]  "
              f"within floor {int((d <= FLOOR).sum())}/{len(rows)}  within 0.02 {int((d <= FLOOR_TAIL).sum())}/{len(rows)}  "
              f"max |diff| {d.max():.4f}")
    # the m-ladder slope on the chain: how much RMSD one unit of m is worth
    mm = np.array([out["cells"][t]["m"] for t in tags], float)
    yy = np.array([out["cells"][t]["mean_chain"] for t in tags], float)
    if len(mm) >= 2:
        A = np.column_stack([mm, np.ones_like(mm)])
        sl = float(np.linalg.lstsq(A, yy, rcond=None)[0][0])
        out["m_ladder_slope_A_per_unit_m"] = sl
        out["m_ladder_range"] = [float(mm.min()), float(mm.max())]
        print(f"  m-ladder slope on the chain: {sl:+.5f} A per unit of m over m in [{mm.min():.0f}, {mm.max():.0f}] "
              f"(the ONE scalar the quantum stage still chooses)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["profile", "cloud", "chain", "analyse_chain"])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--recompute", action="store_true")
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    if a.mode == "profile":
        res = profile_report()
        ST.save_atomic(os.path.join(RESULTS, "s29_D_m6_profile.json"), res, module_file=__file__)
    elif a.mode == "cloud":
        prof = profile_report()
        res = cloud(a.seed, m_star=prof["m_star"])
        if a.recompute:
            res["recomputed"] = cloud_recompute(a.seed, m_star=prof["m_star"])
        res["profile"] = prof
        for tag, c in res["cells"].items():
            print(c["fmt"])
        ST.save_atomic(os.path.join(RESULTS, f"s29_D_m6_cloud_seed{a.seed}.json"), res, module_file=__file__)
    elif a.mode == "chain":
        prof = profile_report()
        chain(a.seed, m_star=prof["m_star"], limit=a.limit)
    else:
        res = analyse_chain(a.seed)
        ST.save_atomic(os.path.join(RESULTS, f"s29_D_m6_chain_seed{a.seed}.json"), res, module_file=__file__)
    print("ok")


if __name__ == "__main__":
    main()
