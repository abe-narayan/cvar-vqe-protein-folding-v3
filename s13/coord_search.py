"""SPRINT 13, COORDINATOR EXPERIMENT 3 -- what does a NATIVE-FREE search in torsion space emit?

Everything else in this sprint is upstream of this number.  `s13/ceiling.py` showed the
discrete torsion space CONTAINS a 1.594 A answer at k=4 in ~26 qubits.  This asks the only
question that matters next: **can a native-free objective, optimised properly, find it?**

It is also the sprint's mandatory classical control.  The literature survey turned up
Boulebnane et al. (npj QI 2023), who found that on comparable problems "the performance of
QAOA can be matched by random sampling up to a small overhead", casting "serious doubt" on
the approach even in an extremely simplified noiseless setting.  So **uniform random sampling
at matched objective-evaluation budget is an arm here, not an afterthought.**  Any quantum
result later in this sprint has to beat this table, not just the retrieval baseline.

ARMS, all at a MATCHED objective-evaluation budget, all native-free at inference:

    random        uniform sampling of the configuration space; return the argmin.  THE CONTROL.
    sa_legacy     simulated annealing on the genuine 11-term Legacy energy.
    sa_prior      SA on the 1-local torsion prior alone (class-state occupancy; no learned
                  prediction, so this is what the representation gives you for free).
    sa_both       SA on prior + Legacy, each standardised on the same sample so the weights
                  are not tuned.
    sa_amber_nc   SA on AMBER electrostatics + solvation ("nonclash"), a MODIFIED potential,
                  named as one, on a subset only -- 28 ms per distinct call.
    ORACLE_desc   coordinate descent on true CA-RMSD.  The representation ceiling. DIAGNOSTIC.

Reported: emitted CA-RMSD of each arm's returned configuration (the built chain is already
ideal geometry, so no projection is needed and none is applied), against the 3.213 A shipped
retrieval pipeline and the 1.594 A representation ceiling.

Native quantities are labels only.  `ORACLE_desc` is the sole arm that reads one.

    python -m s13.coord_search [budget] [n_targets] [--amber]
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
from s13 import ceiling as CE              # noqa: E402
from s13.coord_objval import prior_logp    # noqa: E402
import peptide_db as pdb                   # noqa: E402
import torsion_lib2 as tl2                 # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
os.makedirs(RESULTS, exist_ok=True)
K = 4
SEED = 20260905


class Budget:
    """Exact objective-evaluation accounting. Budget parity is a correctness condition."""

    def __init__(self, cap):
        self.cap = int(cap); self.used = 0

    def take(self, m=1):
        if self.used + m > self.cap:
            return False
        self.used += m
        return True


def legacy_batch(seq, PHI, PSI, S):
    import core.energy as EN
    from core import geometry as geo
    n = S.shape[1]; idx = np.arange(n)
    phi = PHI[idx[None, :], S]; psi = PSI[idx[None, :], S]
    BB = geo.build_backbone_batch(phi, psi)
    return np.asarray(EN.totals_batch(EN.components_batch(seq, BB, phi=phi, psi=psi)), float)


def amber_nc_batch(seq, rep, S):
    """AMBER electrostatics + implicit solvation only. A MODIFIED potential, named as one."""
    import core.amber as AM
    out = np.full(len(S), np.nan)
    for _ in range(120):
        if AM.memory_percent() < AM.MEMORY_LIMIT_PERCENT - 1.0:
            break
        time.sleep(10.0)
    for b, s in enumerate(S):
        try:
            r = AM.single_point(seq, rep, s, components=True)
            c = r.get("components", {}) or {}
            out[b] = float(c.get("nonbonded_elec", c.get("electrostatic", 0.0))
                           + c.get("gb", c.get("solvation", 0.0)))
        except Exception:                                              # noqa: BLE001
            pass
    return out


def anneal(n, k, energy_of, budget, rng, batch=32, t0=None, t1=None):
    """Batched simulated annealing over per-residue states. Charges every evaluation."""
    S = rng.integers(0, k, size=(batch, n))
    if not budget.take(batch):
        return S[0], np.inf
    E = energy_of(S)
    fin = np.isfinite(E)
    scale = float(np.nanstd(E[fin])) if fin.sum() > 2 else 1.0
    t0 = t0 if t0 is not None else max(scale, 1e-6)
    t1 = t1 if t1 is not None else max(scale * 1e-3, 1e-9)
    best_i = int(np.nanargmin(np.where(fin, E, np.inf)))
    best_S, best_E = S[best_i].copy(), float(E[best_i])
    steps = max(1, (budget.cap - budget.used) // batch)
    for it in range(steps):
        T = t0 * (t1 / t0) ** (it / max(steps - 1, 1))
        P = S.copy()
        r = rng.integers(0, n, size=batch)
        P[np.arange(batch), r] = rng.integers(0, k, size=batch)
        if not budget.take(batch):
            break
        Ep = energy_of(P)
        d = Ep - E
        acc = (d < 0) | (rng.random(batch) < np.exp(-np.clip(d, 0, 50) / max(T, 1e-12)))
        acc &= np.isfinite(Ep)
        S[acc] = P[acc]; E[acc] = Ep[acc]
        j = int(np.nanargmin(np.where(np.isfinite(E), E, np.inf)))
        if np.isfinite(E[j]) and E[j] < best_E:
            best_E, best_S = float(E[j]), S[j].copy()
    return best_S, best_E


def run_target(t, budget_cap=5000, want_amber=False):
    rng = np.random.default_rng(SEED + hash(t["pdb"]) % 99991)
    seq, n = t["seq"], t["n"]
    u = I.load_univ(t["pdb"]); nat = u["nat_ca"]
    tab = tl2.library_for(seq, K, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    PHI = np.asarray(rep._phi, float); PSI = np.asarray(rep._psi, float)
    idx = np.arange(n)

    def rmsd_of(S):
        S = np.atleast_2d(S)
        W = I.build_ca(PHI[idx[None, :], S], PSI[idx[None, :], S])
        return I.kabsch_rmsd_batch(W, nat)

    lp = prior_logp(seq, K, seq)

    def f_legacy(S):
        return legacy_batch(seq, PHI, PSI, S)

    def f_prior(S):
        return -lp[idx[None, :], S].sum(1)

    #: standardisation for the combined objective is computed on ONE shared random sample,
    #: not tuned against any reported number.
    S0 = rng.integers(0, K, size=(256, n))
    l0 = f_legacy(S0); p0 = f_prior(S0)
    ls, lm = float(np.nanstd(l0)) or 1.0, float(np.nanmean(l0))
    ps, pm = float(np.nanstd(p0)) or 1.0, float(np.nanmean(p0))

    def f_both(S):
        return (f_legacy(S) - lm) / ls + (f_prior(S) - pm) / ps

    out = {"pdb": t["pdb"], "n": n, "fold": t["fold"], "budget": int(budget_cap),
           "qubits": int(n * np.log2(K))}

    # --- THE CONTROL: uniform random sampling at the same budget
    b = Budget(budget_cap)
    Sr = rng.integers(0, K, size=(budget_cap, n)); b.take(budget_cap)
    Er = f_legacy(Sr)
    out["random_legacy"] = float(rmsd_of(Sr[int(np.nanargmin(Er))])[0])
    out["random_best_possible"] = float(np.nanmin(rmsd_of(Sr)))          # ORACLE diagnostic

    for name, fn in (("sa_legacy", f_legacy), ("sa_prior", f_prior), ("sa_both", f_both)):
        bb = Budget(budget_cap)
        s, e = anneal(n, K, fn, bb, np.random.default_rng(rng.integers(1 << 30)))
        out[name] = float(rmsd_of(s)[0]); out[name + "_evals"] = bb.used
        out[name + "_E"] = float(e)

    if want_amber:
        bb = Budget(min(budget_cap, 1500))
        s, e = anneal(n, K, lambda S: amber_nc_batch(seq, rep, S), bb,
                      np.random.default_rng(rng.integers(1 << 30)), batch=16)
        out["sa_amber_nc"] = float(rmsd_of(s)[0]); out["sa_amber_nc_evals"] = bb.used

    p = pdb.by_pdb(t["pdb"])
    s_snap = CE.snap(PHI, PSI, np.asarray(p.phi, float), np.asarray(p.psi, float))
    s_desc, r_desc = CE.descent(PHI, PSI, nat, s_snap)
    out["ORACLE_descent"] = float(r_desc)                                # ceiling, DIAGNOSTIC
    out["ORACLE_snap"] = float(rmsd_of(s_snap)[0])
    return out


def main(budget_cap=5000, n_targets=126, want_amber=False):
    tg = I.targets()
    if n_targets < len(tg):
        step = max(1, len(tg) // n_targets)
        tg = [tg[i] for i in range(0, len(tg), step)][:n_targets]
    tag = f"b{budget_cap}" + ("_amber" if want_amber else "")
    path = os.path.join(RESULTS, f"coord_search_{tag}.json")
    rows = json.load(open(path))["per_target"] if os.path.exists(path) else []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for i, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        rows.append(run_target(t, budget_cap, want_amber))
        if i % 5 == 0 or i == len(tg) - 1:
            json.dump({"what": "native-free search in discrete torsion space", "k": K,
                       "budget": budget_cap, "per_target": rows}, open(path, "w"), indent=1)
            r = rows[-1]
            print(f"  {i+1}/{len(tg)} {t['pdb']} rand={r['random_legacy']:.2f} "
                  f"sa_leg={r['sa_legacy']:.2f} sa_pri={r['sa_prior']:.2f} "
                  f"sa_both={r['sa_both']:.2f} ORACLE={r['ORACLE_descent']:.2f} "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    json.dump({"what": "native-free search in discrete torsion space", "k": K,
               "budget": budget_cap, "per_target": rows}, open(path, "w"), indent=1)
    report(rows, tag)


def report(rows, tag=""):
    f18 = set(I.FAIL18); isf = np.array([r["pdb"] in f18 for r in rows])
    folds = [r["fold"] for r in rows]; names = [r["pdb"] for r in rows]
    arms = [a for a in ("random_legacy", "sa_prior", "sa_legacy", "sa_both", "sa_amber_nc",
                        "ORACLE_snap", "ORACLE_descent") if a in rows[0]]
    g = lambda a: np.array([r.get(a, np.nan) for r in rows], float)      # noqa: E731
    out = {"n": len(rows), "arms": {}}
    print(f"\n{len(rows)} targets, k={K}, budget={rows[0]['budget']} evaluations, "
          f"mean {rows[0]['qubits']}+ qubits.  Shipped retrieval pipeline = 3.213 A.\n")
    hdr = (f"{'arm':16s} {'mean':>7s} {'median':>7s} {'<2A':>6s} {'<1.5A':>6s} "
           f"{'FAIL18':>7s} {'other108':>8s}")
    print(hdr); print("-" * len(hdr))
    for a in arms:
        v = g(a); ok = np.isfinite(v)
        out["arms"][a] = {"summary": I.summary(v[ok]),
                          "FAIL18": float(np.nanmean(v[isf])),
                          "other108": float(np.nanmean(v[~isf]))}
        print(f"{a:16s} {np.nanmean(v):7.3f} {np.nanmedian(v):7.3f} "
              f"{np.nanmean(v < 2.0):6.2f} {np.nanmean(v < 1.5):6.2f} "
              f"{np.nanmean(v[isf]):7.3f} {np.nanmean(v[~isf]):8.3f}")
    base = g("random_legacy")
    out["paired_vs_random"] = {}
    print("\npaired against the RANDOM-SAMPLING control at matched budget "
          "(negative = the search beats random):")
    for a in arms:
        if a.startswith("ORACLE") or a == "random_legacy":
            continue
        v = g(a); ok = np.isfinite(v) & np.isfinite(base)
        d = I.paired(v[ok], base[ok], folds=list(np.array(folds)[ok]),
                     names=list(np.array(names)[ok]))
        out["paired_vs_random"][a] = d
        print(f"  {a:16s} {d['mean_diff']:+.4f} [{d['ci95'][0]:+.4f},{d['ci95'][1]:+.4f}] "
              f"{d['n_better']}W/{d['n_worse']}L drop10 {d['drop_top10_mean_diff']:+.4f}")
    with open(os.path.join(RESULTS, f"coord_search_{tag}_report.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    bc = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 126
    am = "--amber" in sys.argv
    main(bc, nt, am)
