"""s12 QUANTUM-ROLE -- null controls for the one attractive correlation in this work.

Section 3c found `corr(F, CA-RMSD) = +0.449` over the 4,096 assemblies of a target and the
exact optimum at the 24.8th RMSD percentile -- i.e. the shipped objective looks like it
ranks inside this anchored family, which is NOT what the record says it does elsewhere.
The sprint requires a null control for any attractive correlation.  Two are run:

NULL 1 -- shuffled objective.  The per-pair risk rows of the shipped distogram are
permuted across pairs (same 17-bin marginal shapes, same everything, wrong pair
assignment) and the whole 2-local Hamiltonian is rebuilt from the same placed pieces.  If
a shuffled distogram ranks the assemblies as well, the correlation is a property of the
geometry family, not of the prediction.

NULL 2 -- radius-of-gyration control.  `rg` of an assembly drives both its risk (a
too-compact chain hits low-distance bins) and its CA-RMSD.  The partial correlation
`corr(F, RMSD | rg)` says how much of the +0.449 survives once compactness is removed.
This is the control that overturned "physics ranks real geometry" in the record.

Usage: python -m s12.vq_null <family> [n_targets]
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import vq_lib as V
from s12.vq_build import FAMILIES
from s12.vq_run import native


def tables_from_risk(inst, risk, dg_grid):
    """Rebuild (h, J) from the instance's PLACED pieces against an arbitrary risk table."""
    k, m = inst["k"], inst["m"]
    pi, pj, seg_of = inst["pi"], inst["pj"], inst["seg_of"]
    ss, tt = seg_of[pi], seg_of[pj]
    h = [np.zeros(m) for _ in range(k)]
    J = {}
    for s in range(k):
        msk = (ss == s) & (tt == s)
        if msk.any():
            a0 = inst["segs"][s][0]
            d = np.linalg.norm(inst["placed"][s][:, pi[msk] - a0, :]
                               - inst["placed"][s][:, pj[msk] - a0, :], axis=-1)
            g = np.clip(((d - dg_grid[0]) / 0.05).astype(np.int32), 0, len(dg_grid) - 1)
            R = np.asarray(risk[msk], np.float64)
            h[s] = R[np.arange(msk.sum())[None, :], g].sum(1)
    for s in range(k):
        for u in range(s + 1, k):
            msk = (ss == s) & (tt == u)
            if not msk.any():
                J[(s, u)] = np.zeros((m, m)); continue
            a0, b0 = inst["segs"][s][0], inst["segs"][u][0]
            X = inst["placed"][s][:, pi[msk] - a0, :]
            Y = inst["placed"][u][:, pj[msk] - b0, :]
            d = np.linalg.norm(X[:, None, :, :] - Y[None, :, :, :], axis=-1)
            g = np.clip(((d - dg_grid[0]) / 0.05).astype(np.int32), 0, len(dg_grid) - 1)
            R = np.asarray(risk[msk], np.float64)
            J[(s, u)] = R[np.arange(msk.sum())[None, None, :], g].sum(2)
    return h, J


def partial_corr(a, b, c):
    """corr(a, b | c) by residualising both on c."""
    a, b, c = [np.asarray(x, float) for x in (a, b, c)]
    A = np.c_[np.ones(len(c)), c]
    ra = a - A @ np.linalg.lstsq(A, a, rcond=None)[0]
    rb = b - A @ np.linalg.lstsq(A, b, rcond=None)[0]
    if ra.std() < 1e-12 or rb.std() < 1e-12:
        return np.nan
    return float(np.corrcoef(ra, rb)[0, 1])


def run(fam="A2_64", ntarg=None, n_shuffle=3, seed=0, out=None):
    cfgf = FAMILIES[fam]
    tg = [t for t in I.targets() if V.segments(t["n"], cfgf["k"], cfgf["minlen"]) is not None]
    if ntarg:
        tg = tg[:: max(1, len(tg) // ntarg)][:ntarg]
    rows = []
    t00 = time.time()
    for q, t in enumerate(tg):
        inst = V.cached_instance(t, **cfgf)
        nat = native(t["pdb"])
        cfgs, E = V.enumerate_all(inst)
        Xall = V.assemble_batch(inst, cfgs)
        rr = I.kabsch_rmsd_batch(Xall, nat)
        rg = np.sqrt(((Xall - Xall.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        rg_nat = float(np.sqrt(((nat - nat.mean(0)) ** 2).sum(1).mean()))
        drg = np.abs(rg - rg_nat)
        row = dict(pdb=t["pdb"], n=t["n"], fold=t["fold"],
                   corr_E_rmsd=float(np.corrcoef(E, rr)[0, 1]),
                   corr_E_rg=float(np.corrcoef(E, rg)[0, 1]),
                   corr_rmsd_rg=float(np.corrcoef(rr, rg)[0, 1]),
                   partial_corr_E_rmsd_given_rg=partial_corr(E, rr, rg),
                   pct_of_opt=float((rr < rr[int(np.argmin(E))]).mean()),
                   pct_of_rg_opt=float((rr < rr[int(np.argmin(drg))]).mean()),
                   rmsd_of_opt=float(rr[int(np.argmin(E))]),
                   rmsd_of_rg_opt=float(rr[int(np.argmin(drg))]),
                   rmsd_best=float(rr.min()), rmsd_mean=float(rr.mean()),
                   shuffles=[])
        dg = inst["dg"]
        rng = np.random.default_rng(seed)
        for s in range(n_shuffle):
            perm = rng.permutation(dg["risk"].shape[0])
            hN, JN = tables_from_risk(inst, dg["risk"][perm], dg["grid"])
            iN = dict(inst); iN["h"] = hN; iN["J"] = JN
            EN = V.objective_tables(iN, cfgs)
            row["shuffles"].append(dict(
                corr_E_rmsd=float(np.corrcoef(EN, rr)[0, 1]),
                pct_of_opt=float((rr < rr[int(np.argmin(EN))]).mean()),
                rmsd_of_opt=float(rr[int(np.argmin(EN))]),
                partial_corr_E_rmsd_given_rg=partial_corr(EN, rr, rg)))
        rows.append(row)
        if q % 20 == 0:
            print(f"  [null] {q+1}/{len(tg)} {time.time()-t00:.0f}s free={I.free_gb():.2f}",
                  flush=True)
    I.write(out or f"vq_null_{fam}", dict(family=fam, n_shuffle=n_shuffle, rows=rows))
    # pooled
    f18 = np.array([r["pdb"] in I.FAIL18 for r in rows])
    def col(k):
        return np.array([r[k] for r in rows], float)
    sh = lambda k: np.array([np.mean([s[k] for s in r["shuffles"]]) for r in rows], float)
    summ = {
        "n": len(rows),
        "corr_E_rmsd": [float(col("corr_E_rmsd").mean()), float(col("corr_E_rmsd")[f18].mean()),
                        float(col("corr_E_rmsd")[~f18].mean())],
        "corr_E_rmsd_SHUFFLED": [float(sh("corr_E_rmsd").mean()), float(sh("corr_E_rmsd")[f18].mean()),
                                 float(sh("corr_E_rmsd")[~f18].mean())],
        "partial_given_rg": [float(np.nanmean(col("partial_corr_E_rmsd_given_rg"))),
                             float(np.nanmean(col("partial_corr_E_rmsd_given_rg")[f18])),
                             float(np.nanmean(col("partial_corr_E_rmsd_given_rg")[~f18]))],
        "partial_given_rg_SHUFFLED": [float(np.nanmean(sh("partial_corr_E_rmsd_given_rg"))),
                                      float(np.nanmean(sh("partial_corr_E_rmsd_given_rg")[f18])),
                                      float(np.nanmean(sh("partial_corr_E_rmsd_given_rg")[~f18]))],
        "corr_E_rg": float(col("corr_E_rg").mean()),
        "corr_rmsd_rg": float(col("corr_rmsd_rg").mean()),
        "pct_of_opt": float(col("pct_of_opt").mean()),
        "pct_of_opt_SHUFFLED": float(sh("pct_of_opt").mean()),
        "pct_of_rg_opt": float(col("pct_of_rg_opt").mean()),
        "rmsd_of_opt": float(col("rmsd_of_opt").mean()),
        "rmsd_of_opt_SHUFFLED": float(sh("rmsd_of_opt").mean()),
        "rmsd_of_rg_opt": float(col("rmsd_of_rg_opt").mean()),
        "rmsd_best": float(col("rmsd_best").mean()),
        "rmsd_mean": float(col("rmsd_mean").mean()),
        "paired_opt_vs_shuffled_rmsd": I.paired(col("rmsd_of_opt"), sh("rmsd_of_opt"),
                                                folds=np.array([r["fold"] for r in rows])),
    }
    I.write(out or f"vq_null_{fam}", dict(family=fam, n_shuffle=n_shuffle,
                                          summary=summ, rows=rows))
    print(json.dumps(summ, indent=1))
    return summ


if __name__ == "__main__":
    fam = sys.argv[1] if len(sys.argv) > 1 else "A2_64"
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else None
    run(fam, nt)
