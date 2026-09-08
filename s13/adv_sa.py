"""SPRINT 13 -- ADVERSARIAL AUDIT 4: is the annealer optimising, or is the arm broken?

C3(a) claims "optimising Legacy is worse than random sampling" (+0.102 A, 58W/67L).  That
conclusion is only about the OBJECTIVE if the OPTIMISER works.  `coord_search.anneal` runs
32 INDEPENDENT chains for (5000-32)//32 = 155 single-residue moves each -- so each chain
sees 155 proposals over ~13 residues while the random control sees 5,000 draws.  That is a
plausible way to measure a broken annealer rather than a broken objective.

Arms, all native-free, all charged the SAME 5,000 Legacy evaluations:

    random        5,000 uniform draws, argmin Legacy                       (the control)
    sa32          coord_search.anneal verbatim, batch=32                   (the sprint's arm)
    sa8 / sa1     the same annealer at batch 8 and batch 1 (155 -> 4,968 sequential moves)
    cd_legacy     multi-restart coordinate descent ON LEGACY, native-free  (strong optimiser)
    deep          sa1 at a 100,000-evaluation budget                       (20x the budget)

The discriminating quantity is not RMSD, it is **the Legacy energy each arm reaches**.  If
no arm beats `random`'s minimum energy, the arm measures the annealer.  If the arms DO reach
much lower energy and the structures get WORSE, C3(a) is strengthened and the mechanism is
in plain view.

    python -m s13.adv_sa [n_targets] [budget]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s13 import coord_search as CS         # noqa: E402
import torsion_lib2 as tl2                 # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
K = 4
SEED = 20260905


def cd_legacy(n, k, f, budget, rng, restarts=8):
    """Native-free multi-restart coordinate descent, ALL n*k neighbours per call.

    The whole 1-flip neighbourhood is scored in one batch (best-improvement / Jacobi
    order), which is the same algorithm with far fewer python-level calls.  Every
    candidate is still charged to the budget.
    """
    bestS, bestE = None, np.inf
    idx = np.arange(n)
    while budget.used < budget.cap:
        s = rng.integers(0, k, n)
        if not budget.take(1):
            break
        E = float(f(s[None])[0])
        while budget.used < budget.cap:
            cand = np.repeat(s[None, :], n * k, axis=0)
            cand[np.arange(n * k), np.repeat(idx, k)] = np.tile(np.arange(k), n)
            if not budget.take(n * k):
                break
            e = f(cand)
            b = int(np.nanargmin(e))
            if e[b] >= E - 1e-12:
                break
            E = float(e[b]); s = cand[b].copy()
        if E < bestE:
            bestE, bestS = E, s.copy()
    return bestS, bestE


def run_target(t, budget_cap=5000, deep_cap=50000):
    rng = np.random.default_rng(SEED + hash(t["pdb"]) % 99991)
    seq, n = t["seq"], t["n"]
    u = I.load_univ(t["pdb"]); nat = u["nat_ca"]
    tab = tl2.library_for(seq, K, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    PHI = np.asarray(rep._phi, float); PSI = np.asarray(rep._psi, float)
    idx = np.arange(n)

    def rmsd_of(S):
        S = np.atleast_2d(np.asarray(S, int))
        W = I.build_ca(PHI[idx[None, :], S], PSI[idx[None, :], S])
        return I.kabsch_rmsd_batch(W, nat)

    def f(S):
        return CS.legacy_batch(seq, PHI, PSI, np.atleast_2d(np.asarray(S, int)))

    out = {"pdb": t["pdb"], "n": n, "fold": t["fold"]}

    # --- the control, exactly as coord_search runs it
    Sr = rng.integers(0, K, size=(budget_cap, n))
    Er = f(Sr)
    j = int(np.nanargmin(Er))
    out["random_E"] = float(Er[j]); out["random_rmsd"] = float(rmsd_of(Sr[j])[0])
    out["random_meanE"] = float(np.nanmean(Er)); out["random_sdE"] = float(np.nanstd(Er))
    out["random_best_possible_rmsd"] = float(np.nanmin(rmsd_of(Sr)))

    for nm, batch in (("sa32", 32), ("sa8", 8)):
        b = CS.Budget(budget_cap)
        s, e = CS.anneal(n, K, f, b, np.random.default_rng(rng.integers(1 << 30)), batch=batch)
        out[nm + "_E"] = float(e); out[nm + "_rmsd"] = float(rmsd_of(s)[0]); out[nm + "_evals"] = b.used

    b = CS.Budget(budget_cap)
    s, e = cd_legacy(n, K, f, b, np.random.default_rng(rng.integers(1 << 30)))
    out["cd_legacy_E"] = float(e); out["cd_legacy_rmsd"] = float(rmsd_of(s)[0])
    out["cd_legacy_evals"] = b.used

    if deep_cap:
        # 20x the budget through the FASTEST strong optimiser available: batched SA.
        b = CS.Budget(deep_cap)
        s, e = CS.anneal(n, K, f, b, np.random.default_rng(rng.integers(1 << 30)), batch=128)
        out["deep_E"] = float(e); out["deep_rmsd"] = float(rmsd_of(s)[0]); out["deep_evals"] = b.used
        # and a matched random-sampling control at the SAME 20x budget, chunked so the
        # pairwise term does not allocate a multi-GB array.
        bestE2, bestS2 = np.inf, None
        for _ in range(deep_cap // budget_cap):
            Sc = rng.integers(0, K, size=(budget_cap, n)); Ec = f(Sc)
            jc = int(np.nanargmin(Ec))
            if Ec[jc] < bestE2:
                bestE2, bestS2 = float(Ec[jc]), Sc[jc].copy()
        out["random20x_E"] = float(bestE2); out["random20x_rmsd"] = float(rmsd_of(bestS2)[0])

    # ORACLE labels for interpretation only
    out["ORACLE_native_snapE"] = None
    return out


def report(rows):
    arms = ["random", "sa32", "sa8", "cd_legacy", "deep", "random20x"]
    g = lambda a, f: np.array([r.get(a + "_" + f, np.nan) for r in rows], float)   # noqa: E731
    base_E = g("random", "E"); base_R = g("random", "rmsd")
    print(f"\n{len(rows)} targets, k={K}, Legacy objective.  "
          f"'E' is the LOWEST Legacy energy the arm reached (lower = better optimisation).\n")
    hdr = f"{'arm':12s} {'mean E':>10s} {'dE vs rand':>11s} {'mean RMSD':>10s} {'dRMSD':>8s} {'W/L rmsd':>10s}"
    print(hdr); print("-" * len(hdr))
    out = {"n": len(rows), "arms": {}}
    for a in arms:
        E = g(a, "E"); R = g(a, "rmsd")
        if not np.isfinite(E).any():
            continue
        d = I.paired(R[np.isfinite(R)], base_R[np.isfinite(R)])
        out["arms"][a] = {"mean_E": float(np.nanmean(E)), "mean_dE_vs_random": float(np.nanmean(E - base_E)),
                          "n_E_lower_than_random": int(np.nansum(E < base_E - 1e-9)),
                          "rmsd": I.summary(R[np.isfinite(R)]), "paired_rmsd_vs_random": d}
        print(f"{a:12s} {np.nanmean(E):10.3f} {np.nanmean(E-base_E):11.3f} {np.nanmean(R):10.3f} "
              f"{d['mean_diff']:+8.3f} {d['n_better']:4d}/{d['n_worse']:<4d}")
    # the mechanism: does deeper optimisation buy worse structure?
    dE = g("deep", "E") - base_E; dR = g("deep", "rmsd") - base_R
    ok = np.isfinite(dE) & np.isfinite(dR)
    if ok.sum() > 5:
        from scipy import stats
        s = stats.spearmanr(dE[ok], dR[ok])
        out["depth_vs_damage"] = {"spearman_dE_vs_dRMSD": [float(s.statistic), float(s.pvalue)],
                                  "n": int(ok.sum())}
        print(f"\nrho(energy gained below random, RMSD gained above random) over "
              f"{ok.sum()} targets = {s.statistic:+.3f}  p={s.pvalue:.2g}")
        print("   (negative = the deeper the optimisation, the WORSE the structure)")
    return out


def main(n_targets=126, budget=5000, deep=50000):
    tg = I.targets()
    sel = tg if n_targets >= len(tg) else [tg[i] for i in range(0, len(tg), max(1, len(tg) // n_targets))][:n_targets]
    path = os.path.join(RESULTS, "adv_sa.json")
    rows = json.load(open(path))["per_target"] if os.path.exists(path) else []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for i, t in enumerate(sel):
        if t["pdb"] in done:
            continue
        rows.append(run_target(t, budget, deep))
        if i % 5 == 0 or i == len(sel) - 1:
            json.dump({"what": "is the annealer optimising", "budget": budget,
                       "deep_budget": deep, "per_target": rows}, open(path, "w"), indent=1)
            r = rows[-1]
            print(f"  {i+1}/{len(sel)} {t['pdb']} randE={r['random_E']:.2f} sa32E={r['sa32_E']:.2f} "
                  f"cdE={r['cd_legacy_E']:.2f} deepE={r['deep_E']:.2f} | "
                  f"randR={r['random_rmsd']:.2f} deepR={r['deep_rmsd']:.2f} [{time.time()-t0:.0f}s]",
                  flush=True)
    json.dump({"what": "is the annealer optimising", "budget": budget, "deep_budget": deep,
               "per_target": rows}, open(path, "w"), indent=1)
    rep = report(rows)
    json.dump(rep, open(os.path.join(RESULTS, "adv_sa_report.json"), "w"), indent=1)


if __name__ == "__main__":
    nt = int(sys.argv[1]) if len(sys.argv) > 1 else 126
    bg = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
    main(nt, bg)
