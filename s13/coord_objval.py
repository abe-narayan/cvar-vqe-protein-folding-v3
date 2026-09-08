"""SPRINT 13, COORDINATOR EXPERIMENT 1 -- does ANY energy rank the native in torsion space?

This is the number the sprint lives or dies on, and it must be measured before anything is
built on top of it.  `s13/ceiling.py` established that the discrete torsion space is not the
barrier: at k=4 (26 qubits) it contains a 1.594 A answer, at k=8 (39 qubits) a 1.184 A one.
So the architecture now depends entirely on whether a NATIVE-FREE objective prefers those
configurations over the rest of the space.

Sprint 12's hardest lesson was that the certified exact optimum of its assembly Hamiltonian
emitted 0.169 A WORSE than its own anchor, because the objective was anti-correlated with
truth.  Optimising an objective that does not rank the native harder makes the structure
worse.  So: measure the objective first, optimise second.

METHOD.  Per target, draw N configurations uniformly from the k-state torsion space, plus
three reference configurations that are NOT random: the native snapped to its nearest library
states, the coordinate-descent optimum from `s13.ceiling`, and the all-zeros state.  For every
configuration compute, in one pass, the true CA-RMSD (ORACLE label, evaluation only) and every
candidate objective:

    legacy        the genuine 11-term knowledge-based energy (`core.energy.components_batch`)
    amber         genuine ff14SB/GBn2 single-point (`core.amber.single_point`, 6 ms warm)
    amber_min     the same after restrained minimisation, on a small subset only (5.9 s each)
    prior         the library's own state log-frequency -- a 1-local torsion prior with NO
                  sequence-specific prediction in it, so it is the null a real predictor
                  must beat

Then report, per objective: Spearman rho(E, RMSD), the native-snap's percentile under E, the
RMSD of the objective's argmin over the sample, and how that compares with a random draw and
with the sample's own best.  A useful objective has rho well above 0, puts the native near the
bottom, and has an argmin better than random.

EVERYTHING NATIVE HERE IS A LABEL.  `rr` / `nat_ca` / native torsions are used to SCORE and to
place the native reference; no objective reads them.  All arms are DIAGNOSTIC.

    python -m s13.coord_objval [n_samples] [k]
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
import peptide_db as pdb                   # noqa: E402
import torsion_lib2 as tl2                 # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
os.makedirs(RESULTS, exist_ok=True)
SEED = 20260905


def prior_logp(seq, k, exclude_seq=None):
    """1-local torsion prior: the OCCUPANCY of each library state, per residue class.

    `library_for` defaults to `mode="class"` -- k states per residue CLASS (GENERAL / GLY /
    PRO / PRE_PRO), clustered over the database with `exclude_seq` held out.  The states
    themselves carry no occupancy, so this recomputes it the only honest way: assign every
    held-out observation of that class to its nearest state and take the frequencies.

    This is a genuine, deployable, leakage-safe torsion prior with NO learned prediction in
    it -- the class of a residue is known from the sequence alone.  It is therefore the NULL
    that any trained per-residue torsion predictor has to beat, and quoting a predictor's
    value without it would be meaningless.
    """
    import representations as reps
    exclude_seq = seq if exclude_seq is None else exclude_seq
    tab = np.asarray(tl2.library_for(seq, k, exclude_seq), float)     # (n, k, 2)
    entries = pdb.holdout(exclude_seq) if exclude_seq else pdb.load()
    lib = tl2.ContextLibrary(entries)
    cls = reps.residue_classes(seq, len(seq))
    occ_by_class = {}
    for c in set(cls):
        pool = np.asarray(lib.pools.get((c,)) or lib.pools[(reps.CLASS_GENERAL,)], float)
        occ_by_class[c] = pool
    out = np.zeros((len(seq), k))
    for i, c in enumerate(cls):
        pool = occ_by_class[c]
        st = tab[i]                                                   # (k, 2)
        d = (np.abs(_wrap_ang(pool[:, None, 0] - st[None, :, 0]))
             + np.abs(_wrap_ang(pool[:, None, 1] - st[None, :, 1])))
        a = np.argmin(d, axis=1)
        cnt = np.bincount(a, minlength=k).astype(float) + 0.5         # Laplace
        out[i] = cnt / cnt.sum()
    return np.log(np.maximum(out, 1e-9))


def _wrap_ang(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


def energies_for(seq, rep, PHI, PSI, S, want_amber=True):
    """(legacy, amber) for a batch of state assignments S (B, n)."""
    import core.energy as EN
    import core.amber as AM
    from core import geometry as geo
    B, n = S.shape
    idx = np.arange(n)
    phi = PHI[idx[None, :], S]; psi = PSI[idx[None, :], S]
    BB = geo.build_backbone_batch(phi, psi)
    comp = EN.components_batch(seq, BB, phi=phi, psi=psi)
    lg = EN.totals_batch(comp)
    am = np.full(B, np.nan)
    if want_amber:
        #: `core.amber` refuses to start an OpenMM context above 92% physical memory, which
        #: a busy agent fleet trips regularly.  Wait for headroom rather than silently
        #: returning NaN, because a NaN column would quietly become "AMBER has no signal".
        for _ in range(60):
            if AM.memory_percent() < AM.MEMORY_LIMIT_PERCENT - 1.0:
                break
            time.sleep(10.0)
        n_fail = 0
        for b in range(B):
            try:
                am[b] = AM.single_point(seq, rep, S[b])["energy"]
            except Exception:                                          # noqa: BLE001
                n_fail += 1
                am[b] = np.nan
        if n_fail:
            print(f"    [amber: {n_fail}/{B} calls failed]", flush=True)
    return np.asarray(lg, float), am, BB


def run_target(t, n_samples=160, k=4, want_amber=True, rng=None):
    rng = rng or np.random.default_rng(SEED + hash(t["pdb"]) % 99991)
    seq, n = t["seq"], t["n"]
    u = I.load_univ(t["pdb"]); nat = u["nat_ca"]
    tab = tl2.library_for(seq, k, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    PHI = np.asarray(rep._phi, float); PSI = np.asarray(rep._psi, float)

    p = pdb.by_pdb(t["pdb"])
    phi0 = np.asarray(p.phi, float); psi0 = np.asarray(p.psi, float)     # ORACLE
    s_snap = CE.snap(PHI, PSI, phi0, psi0)                               # ORACLE reference
    s_desc, _ = CE.descent(PHI, PSI, nat, s_snap)                        # ORACLE reference

    S = rng.integers(0, k, size=(n_samples, n))
    S = np.vstack([S, s_snap[None, :], s_desc[None, :], np.zeros((1, n), int)])
    tag = ["rand"] * n_samples + ["native_snap", "descent", "zeros"]

    lg, am, _ = energies_for(seq, rep, PHI, PSI, S, want_amber)
    W = I.build_ca(PHI[np.arange(n)[None, :], S], PSI[np.arange(n)[None, :], S])
    rr = I.kabsch_rmsd_batch(W, nat)                                     # ORACLE label

    lp = prior_logp(seq, k, seq)
    pr = -lp[np.arange(n)[None, :], S].sum(1)                            # lower = more typical

    out = {"pdb": t["pdb"], "n": n, "fold": t["fold"], "k": k,
           "n_samples": int(n_samples), "qubits": int(n * np.log2(k)),
           "rmsd_rand_mean": float(rr[:n_samples].mean()),
           "rmsd_rand_best": float(rr[:n_samples].min()),
           "rmsd_native_snap": float(rr[n_samples]), "rmsd_descent": float(rr[n_samples + 1]),
           "objectives": {}}
    from scipy.stats import spearmanr
    for name, E in (("legacy", lg), ("amber", am), ("prior", pr)):
        E = np.asarray(E, float)
        ok = np.isfinite(E)
        if ok.sum() < 10:
            out["objectives"][name] = None
            continue
        r = spearmanr(E[ok][:n_samples], rr[ok][:n_samples])
        # where does the native-snapped configuration sit under this objective?
        pct = float(np.mean(E[:n_samples] < E[n_samples])) if np.isfinite(E[n_samples]) else None
        am_idx = int(np.nanargmin(np.where(ok[:n_samples], E[:n_samples], np.inf)))
        out["objectives"][name] = {
            "rho_E_vs_rmsd": float(r.statistic), "p": float(r.pvalue),
            "native_percentile": pct,
            "argmin_rmsd": float(rr[am_idx]),
            "E_native": float(E[n_samples]) if np.isfinite(E[n_samples]) else None,
            "E_descent": float(E[n_samples + 1]) if np.isfinite(E[n_samples + 1]) else None,
            "E_rand_mean": float(np.nanmean(E[:n_samples])),
            "descent_percentile": (float(np.mean(E[:n_samples] < E[n_samples + 1]))
                                   if np.isfinite(E[n_samples + 1]) else None)}
    out["tag"] = tag[:3]
    return out


def main(n_samples=160, k=4, want_amber=True):
    tg = I.targets()
    path = os.path.join(RESULTS, f"coord_objval_k{k}.json")
    rows = json.load(open(path))["per_target"] if os.path.exists(path) else []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for i, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        rows.append(run_target(t, n_samples, k, want_amber))
        if i % 10 == 0 or i == len(tg) - 1:
            json.dump({"what": "objective validity in the discrete torsion space",
                       "n_samples": n_samples, "k": k, "per_target": rows},
                      open(path, "w"), indent=1)
            o = rows[-1]["objectives"]
            print(f"  {i+1}/{len(tg)} {t['pdb']} rand={rows[-1]['rmsd_rand_mean']:.2f} "
                  f"snap={rows[-1]['rmsd_native_snap']:.2f} | "
                  + " ".join(f"{nm}:rho={o[nm]['rho_E_vs_rmsd']:+.2f},pct={o[nm]['native_percentile']:.2f}"
                             for nm in ("legacy", "amber", "prior") if o.get(nm))
                  + f"  [{time.time()-t0:.0f}s]", flush=True)
    json.dump({"what": "objective validity in the discrete torsion space",
               "n_samples": n_samples, "k": k, "per_target": rows}, open(path, "w"), indent=1)
    report(rows, k)


def report(rows, k=4):
    f18 = set(I.FAIL18); isf = np.array([r["pdb"] in f18 for r in rows])
    out = {"n": len(rows), "k": k, "objectives": {}}
    print(f"\nk={k}  random-config mean RMSD {np.mean([r['rmsd_rand_mean'] for r in rows]):.3f}"
          f"   native-snap {np.mean([r['rmsd_native_snap'] for r in rows]):.3f}"
          f"   descent {np.mean([r['rmsd_descent'] for r in rows]):.3f}\n")
    hdr = (f"{'objective':10s} {'rho':>7s} {'nat pct':>8s} {'desc pct':>9s} {'argmin RMSD':>12s} "
           f"{'rand best':>10s} {'F18 rho':>8s}")
    print(hdr); print("-" * len(hdr))
    for nm in ("legacy", "amber", "prior"):
        v = [r["objectives"].get(nm) for r in rows]
        ok = np.array([x is not None for x in v])
        if ok.sum() == 0:
            continue
        rho = np.array([x["rho_E_vs_rmsd"] if x else np.nan for x in v], float)
        pct = np.array([x["native_percentile"] if x and x["native_percentile"] is not None else np.nan for x in v], float)
        dpc = np.array([x["descent_percentile"] if x and x["descent_percentile"] is not None else np.nan for x in v], float)
        arg = np.array([x["argmin_rmsd"] if x else np.nan for x in v], float)
        rb = np.array([r["rmsd_rand_best"] for r in rows], float)
        out["objectives"][nm] = {
            "rho_mean": float(np.nanmean(rho)), "rho_median": float(np.nanmedian(rho)),
            "rho_frac_positive": float(np.nanmean(rho > 0)),
            "native_percentile_mean": float(np.nanmean(pct)),
            "descent_percentile_mean": float(np.nanmean(dpc)),
            "argmin_rmsd_mean": float(np.nanmean(arg)),
            "rand_best_mean": float(np.nanmean(rb)),
            "FAIL18_rho": float(np.nanmean(rho[isf])), "other108_rho": float(np.nanmean(rho[~isf]))}
        o = out["objectives"][nm]
        print(f"{nm:10s} {o['rho_mean']:+7.3f} {o['native_percentile_mean']:8.3f} "
              f"{o['descent_percentile_mean']:9.3f} {o['argmin_rmsd_mean']:12.3f} "
              f"{o['rand_best_mean']:10.3f} {o['FAIL18_rho']:+8.3f}")
    print("\nREAD: rho > 0 means lower energy goes with lower RMSD (good). native/descent "
          "percentile is the FRACTION OF RANDOM CONFIGS THAT SCORE BETTER than the good "
          "structure -- near 0 is what a usable objective looks like, 0.5 is chance.")
    with open(os.path.join(RESULTS, f"coord_objval_k{k}_report.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    ns = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    kk = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    main(ns, kk)
