"""SPRINT 15, coordinator -- the retrieval pool is a SECOND, independent distance channel.

THE OBSERVATION. `s15/distacc.py` showed the whole gap between the oracle distance fit
(0.697 A) and the predicted one (2.79 A) is distogram error, and that the distogram is
systematically biased (+0.113 A at separation 2-3 rising to **+1.492 A at 11-15**) with
uncertainties **2.633x over-confident**. So the restraint channel is the bottleneck, and the
obvious question is whether a better one exists.

One does, and it has been sitting in the universe cache for fourteen sprints. The K=500
BLOSUM-retrieved window pool contains **coordinates**, so every pool member supplies a
complete CA-CA distance matrix. The pool therefore gives a **per-pair empirical distance
distribution** for the target — not a learned prediction, but a direct structural sample from
sequence-matched fragments.

Sprint 14 used this pool for its **torsions** and found the resulting per-residue prior beats
the trained leave-fold-out sequence model (phi 33.6 vs 36.1, psi 59.2 vs 62.4) with no
training at all. **Nobody has ever used its distances.**

WHY IT MIGHT BE BETTER, AND WHY IT MIGHT BE WORSE. Better: it is non-parametric, it needs no
training, it carries genuine multimodality (a pair is often either in contact or not), and its
errors are structurally generated rather than model-generated, so they should be **independent
of the distogram's**. Worse: it inherits the pool's composition, which Sprint 12 measured as
about 80% protein fragments whose sequence-structure channel is far weaker than that of
isolated peptides, and it has no mechanism to know which pool members are relevant.

**Two channels with independent errors is the interesting case**, because Sprint 14 measured
that Legacy and AMBER are genuinely decorrelated (+0.096 truth-partialled) and yet unfusable —
fusion gain goes as the SQUARE of the weaker channel's skill. That arithmetic is the thing to
beat here, and it must be checked before any fusion is claimed.

ARMS
    distogram          the learned leave-fold-out predictor (reference)
    pool_mean          pool circular/arithmetic mean distance, inverse-variance weighted
    pool_hist          maximum likelihood under the pool's empirical distance histogram
    pool_sim_weighted  the same, weighted by BLOSUM similarity to the target
    combined           product of the two likelihoods -- the fusion arm
    ORACLE_true        fit to the true distance matrix, the ceiling

Every arm is native-free except the ORACLE. The pool is built by BLOSUM retrieval from a
leakage-safe library with the target and everything above 0.6 identity removed; the `rr`
column (oracle RMSD) is never read.

Run:
    python -m s15.pooldist
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
    os.environ.setdefault(_v, "2")

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD             # noqa: E402
from s15 import distgeo as D                 # noqa: E402
from s15 import distml as M                  # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)

SEPS = [(2, 3), (4, 5), (6, 7), (8, 10), (11, 15)]


def pool_distances(pdb, n, arm="pool", k=500):
    """(m, npairs) CA-CA distances of every retrieved window, plus BLOSUM similarities.

    NATIVE-FREE: `order` and `sim` are retrieval outputs; `rr` and `nat_ca` are not read.
    """
    u = I.load_univ(pdb)
    p = I.pool_idx(u, k)
    if arm == "top75":
        rec = I.shipped_record(pdb)
        p = p[np.asarray(rec["sub"], int)]
    W = np.asarray(u["W"], float)[p]
    sim = np.asarray(u["sim"], float)[p]
    i, j = I.pair_index(n)
    Dm = np.sqrt(((W[:, i] - W[:, j]) ** 2).sum(-1))          # (m, npairs)
    return Dm, sim, i, j


def hist_table(Dm, centres, w=None, smooth=1.0):
    """Empirical per-pair distance histogram on the distogram's own bin centres."""
    m, P = Dm.shape
    edges = np.concatenate([[centres[0] - (centres[1] - centres[0]) / 2],
                            (centres[:-1] + centres[1:]) / 2,
                            [centres[-1] + (centres[-1] - centres[-2]) / 2]])
    ww = np.ones(m) if w is None else np.asarray(w, float)
    prob = np.empty((P, len(centres)))
    for q in range(P):
        h, _ = np.histogram(Dm[:, q], bins=edges, weights=ww)
        prob[q] = h
    prob = prob + smooth                                        # Laplace, keeps support
    prob = prob / prob.sum(1, keepdims=True)
    return M.LogPTable(prob, centres)


class SumLogP:
    """Product of independent likelihoods = sum of log-densities. The fusion arm."""

    def __init__(self, *tables):
        self.tables = tables

    def __call__(self, d):
        vals = [t(d) for t in self.tables]
        return (sum(v[0] for v in vals), sum(v[1] for v in vals))

    def rowsum(self, D):
        """Per-configuration log-likelihood sum, vectorised over leading axes."""
        return sum(t.rowsum(D) for t in self.tables)


def error_profile(targets=None, k=500):
    """How accurate is the POOL's distance estimate, against the distogram's? ORACLE scoring."""
    tg = targets if targets is not None else I.targets()
    ep, dp, sp = [], [], []
    for t in tg:
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        Dm, sim, i, j = pool_distances(pdb, n, k=k)
        dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
        dg = I.distogram(pdb, seq, fold)
        ep.append(np.median(Dm, axis=0) - dtrue)                # pool median error
        dp.append(np.asarray(dg["expected"], float) - dtrue)    # distogram error
        sp.append((j - i).astype(float))
    e_pool = np.concatenate(ep); e_dg = np.concatenate(dp); s = np.concatenate(sp)
    rows = []
    for lo, hi in SEPS:
        m = (s >= lo) & (s <= hi)
        if not m.any():
            continue
        rows.append({"sep": f"{lo}-{hi}", "n": int(m.sum()),
                     "pool_mae": float(np.abs(e_pool[m]).mean()),
                     "pool_bias": float(e_pool[m].mean()),
                     "disto_mae": float(np.abs(e_dg[m]).mean()),
                     "disto_bias": float(e_dg[m].mean()),
                     "err_corr": float(np.corrcoef(e_pool[m], e_dg[m])[0, 1])})
    glob = {"pool_mae": float(np.abs(e_pool).mean()), "pool_bias": float(e_pool.mean()),
            "disto_mae": float(np.abs(e_dg).mean()), "disto_bias": float(e_dg.mean()),
            "err_corr": float(np.corrcoef(e_pool, e_dg)[0, 1]), "n": int(len(e_pool))}
    return {"by_separation": rows, "global": glob}


def run(targets=None, n_start=6, k=500):
    from s14 import ladder as L
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    fail = np.isin(pdbs, I.FAIL18)

    arms = ["distogram", "pool_mean", "pool_hist", "pool_sim_weighted", "combined",
            "ORACLE_true"]
    res = {a: [] for a in arms}
    path = os.path.join(RESULTS, "pooldist.json")

    for c, t in enumerate(tg):
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        Dm, sim, i, j = pool_distances(pdb, n, k=k)
        dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
        dg = I.distogram(pdb, seq, fold)
        centres = np.asarray(dg["centres"], float)
        tab_dg = M.LogPTable(np.asarray(dg["prob"], float), centres)
        tab_pool = hist_table(Dm, centres)
        wsim = np.exp((sim - sim.max()) / 2.0)
        tab_sim = hist_table(Dm, centres, w=wsim)
        tab_comb = SumLogP(tab_dg, tab_pool)

        d_pool = np.median(Dm, axis=0)
        sd_pool = np.maximum(Dm.std(axis=0), 0.25)
        rng = SD.stable_rng(pdb)
        S = D.starts(pdb, seq, n, fold, n_start, rng)

        for a in arms:
            best = None
            for phi0, psi0, _tag in S:
                if a == "pool_mean":
                    phi, psi, f, _ = D.fit_distances(d_pool, 1.0 / sd_pool ** 2,
                                                     i, j, phi0, psi0)
                elif a == "ORACLE_true":
                    phi, psi, f, _ = D.fit_distances(dtrue, np.ones_like(dtrue),
                                                     i, j, phi0, psi0)
                else:
                    tab = {"distogram": tab_dg, "pool_hist": tab_pool,
                           "pool_sim_weighted": tab_sim, "combined": tab_comb}[a]
                    phi, psi, f = M.fit_ml(tab, i, j, phi0, psi0)
                if best is None or f < best[0]:
                    best = (f, float(I.ca_rmsd(I.build_ca(phi, psi), nat)))
            res[a].append(best[1])

        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": res, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    base = np.asarray(res["distogram"], float)
    out = {"n": len(pdbs), "incumbent": float(ref.mean()), "arms": {},
           "per_target": {a: dict(zip(pdbs, map(float, res[a]))) for a in arms}}
    for a in arms:
        v = np.asarray(res[a], float)
        out["arms"][a] = {**I.summary(v), "median": float(np.median(v)),
                          "FAIL18": float(v[fail].mean()),
                          "frac_under_2_5": float((v < 2.5).mean()),
                          "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs),
                          "vs_distogram": I.paired(v, base, folds=folds, names=pdbs)}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_pooldist", out, n_expected=len(pdbs))

    print(f"\nincumbent {ref.mean():.3f}\n")
    print(f"{'arm':<20}{'mean':>8}{'median':>8}{'<2A':>7}{'<2.5A':>8}{'FAIL18':>9}"
          f"{'vs distogram':>24}{'vs incumbent':>24}")
    for a in arms:
        s = out["arms"][a]
        vd, vi = s["vs_distogram"], s["vs_incumbent"]
        print(f"{a:<20}{s['mean']:>8.3f}{s['median']:>8.3f}{s['frac_under_2.0']:>7.2f}"
              f"{s['frac_under_2_5']:>8.2f}{s['FAIL18']:>9.3f}"
              f"  {vd['mean_diff']:+.3f} [{vd['ci95'][0]:+.3f},{vd['ci95'][1]:+.3f}]"
              f"  {vi['mean_diff']:+.3f} [{vi['ci95'][0]:+.3f},{vi['ci95'][1]:+.3f}]")
    return out


if __name__ == "__main__":
    run()
