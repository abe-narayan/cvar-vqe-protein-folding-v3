"""Coordinator item 3: is there a DEPLOYABLE per-target proxy for objective quality that
lets the pipeline choose the averaging cardinality m per target?

Labels (leave-fold-out training labels only): the per-target best m on the SHIPPED-objective
slice of the operator x quality surface.  Features: deployable per-target statistics only.
Controls: an oracle router (ceiling), a constant router (floor = always m=75), a
label-permutation null, and a per-fold breakdown.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_common as A
from s12 import agg_surface as S
from s12 import agg_features as AF
from s12 import agg_meta

GRIDS = {"g3": [20, 75, 300], "g4": [5, 25 if 25 in S.MS else 20, 75, 220], "g2": [1, 75],
         "g5": [3, 12, 50, 150, 500]}


def features():
    """Deployable per-target features (nothing from rr / nat)."""
    tg = I.targets()
    X, names = [], None
    for t in tg:
        d = A.load(t["pdb"], mem=False)
        o = AF.load(t["pdb"])
        mt = agg_meta.load(t["pdb"])
        sc = d["sc"].astype(float)
        sd = d["sd"].astype(float)
        e = d["exp"].astype(float)
        P = o["P"].astype(float)
        srt = np.sort(sc)
        f = {
            "n": t["n"],
            "sc_mean": sc.mean(), "sc_sd": sc.std(),
            "sc_gap_1_75": srt[74] - srt[0], "sc_gap_75_500": srt[-1] - srt[74],
            "sc_skew": float(((sc - sc.mean()) ** 3).mean() / (sc.std() ** 3 + 1e-9)),
            "dg_sd_mean": sd.mean(), "dg_sd_max": sd.max(),
            "dg_exp_mean": e.mean(), "dg_exp_sd": e.std(),
            "cons_mean": P.mean(), "cons_sd": P.std(),
            "cons_med": float(np.median(P)),
            "sim_mean": float(mt["sim"][o["idx"]].mean()), "sim_sd": float(mt["sim"][o["idx"]].std()),
            "org_frac": float(mt["org"][o["idx"]].mean()),
            "rg_sd": float(np.std(np.sqrt(((o["Asup"].astype(float) -
                     o["Asup"].astype(float).mean(1, keepdims=True)) ** 2).sum(-1).mean(-1)))),
            "sim_full_mean": float(mt["sim"].mean()),
        }
        if names is None:
            names = sorted(f)
        X.append([f[k] for k in names])
    return np.array(X, float), names


def run(out="agg_router"):
    tg = I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.array([t["fold"] for t in tg])
    r = json.load(open(os.path.join(I.RESULTS, "agg_surface.json")))
    X, names = features()
    res = {"features": names, "grids": {}}
    rng = np.random.default_rng(0)
    for gname, grid in GRIDS.items():
        grid = [m for m in grid if m in S.MS]
        V = np.array([[r[p][f"shipped|m{m}"] for m in grid] for p in pdbs])   # (126, |grid|)
        y = V.argmin(1)
        base = V[:, grid.index(75)] if 75 in grid else V[:, -1]
        arms = {"const_m75": base, "ORACLE_router": V.min(1)}
        for tag, perm in (("lfo", False), ("null", True)):
            pred = np.zeros(len(pdbs), int)
            for f in sorted(set(folds)):
                tr = folds != f; te = folds == f
                Xtr = X[tr]; mu = Xtr.mean(0); sg = Xtr.std(0) + 1e-9
                Ztr = (Xtr - mu) / sg; Zte = (X[te] - mu) / sg
                # multi-output ridge on the RMSD of each cardinality (regression, then argmin)
                Ytr = V[tr].copy()
                if perm:
                    Ytr = Ytr[rng.permutation(len(Ytr))]
                Aa = Ztr.T @ Ztr + 30.0 * np.eye(Ztr.shape[1])
                Bb = Ztr.T @ (Ytr - Ytr.mean(0))
                Wc = np.linalg.solve(Aa, Bb)
                pr = Zte @ Wc + Ytr.mean(0)
                pred[te] = pr.argmin(1)
            arms[f"ridge_{tag}"] = V[np.arange(len(pdbs)), pred]
            arms[f"ridge_{tag}_choices"] = None
            res.setdefault("choice", {})[f"{gname}_{tag}"] = [int(grid[p]) for p in pred]
        st = {}
        for k, v in arms.items():
            if v is None:
                continue
            s = I.paired(np.asarray(v, float), base, folds=folds, names=pdbs)
            st[k] = {"mean": float(np.mean(v)), "d": s["mean_diff"], "ci": s["ci95"],
                     "WL": f"{s['n_better']}/{s['n_worse']}", "drop10": s["drop_top10_mean_diff"],
                     "per_fold": s["per_fold"]}
        res["grids"][gname] = {"grid": grid, "stats": st,
                               "label_hist": np.bincount(y, minlength=len(grid)).tolist()}
        print(gname, json.dumps({k: round(v["mean"], 4) for k, v in st.items()}), flush=True)
    I.write(out, res)


if __name__ == "__main__":
    run()
