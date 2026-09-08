"""COORDINATOR EXPERIMENT 1 -- how much does the architecture know about the SEQUENCE?

Every number in this project is quoted against the shipped distogram argmin (3.454).  No
arm has ever been quoted against the one control that prices the whole architecture:

    what does the pipeline emit when it is told NOTHING about the target sequence?

The pipeline has exactly two sequence-conditioned stages -- BLOSUM62 retrieval (which 500
windows enter the pool) and the distogram filter (which 75 of them survive).  Everything
else (superpose, average, project, relax) is sequence-blind geometry.  So a 2x2 ablation
separates them, and its lower-left cell is the "typical peptide-like conformation of this
length" baseline that the entire architecture must beat to have learned anything at all.

    arm            retrieval        filter          what it isolates
    ---------------------------------------------------------------------------
    shipped        BLOSUM top-500   distogram 75    the incumbent (3.204 emitted)
    no_filter      BLOSUM top-500   random 75       what the distogram adds
    no_retrieval   random 500       distogram 75    what BLOSUM adds
    blind          random 500       random 75       NOTHING about the sequence

Random arms are averaged over R seeds so the comparison is against the EXPECTATION of the
blind pipeline, not one draw of it.  Every arm goes through the identical downstream path
(medoid superposition -> coordinate average -> multi-start projection at lam=0.3), so the
only thing that differs between cells is which windows the operator sees.

This is a DIAGNOSTIC.  Nothing here is deployable and nothing here is tuned.

    python -m s12.coord_null [n_seeds]
"""
import os
import sys
import json
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402

K = 500
M = 75
SEED = 20260905


def one_target(t, seeds):
    u = I.load_univ(t["pdb"])
    n, seq, fold = t["n"], t["seq"], t["fold"]
    nat = u["nat_ca"]
    W_all = u["W"]
    nw = len(W_all)
    dg = I.distogram(t["pdb"], seq, fold)
    i, j = I.pair_index(n)
    blos = I.pool_idx(u, K)

    def emit(idx):
        """The pipeline's terminal operator on a candidate index set -> emitted CA-RMSD."""
        Wq = W_all[idx].astype(np.float32).astype(float)      # reference float32 round trip
        C, _ = I.coordinate_average(Wq)
        out = I.project(C, seq, fold, lam=0.3)
        return (I.ca_rmsd(out["ca"], nat), I.ca_rmsd(out["fit_ca"], nat))

    def top_m(pool):
        D = I.pair_dists(W_all[pool].astype(np.float32).astype(float), i, j)
        sc = I.shipped_score(dg, D)
        return pool[np.argsort(sc, kind="stable")[:M]]

    res = {"pdb": t["pdb"], "n": n, "fold": fold, "n_windows": int(nw),
           "pool_best_blosum": float(u["rr"][blos].min()),
           "pool_mean_blosum": float(u["rr"][blos].mean())}

    a, f = emit(top_m(blos));       res["shipped_arm"], res["shipped_fit"] = a, f
    res["shipped_sel"] = float(u["rr"][blos][int(np.argmin(
        I.shipped_score(dg, I.pair_dists(W_all[blos].astype(np.float32).astype(float), i, j))))])

    for name, get in (("no_filter", "bf"), ("no_retrieval", "rd"), ("blind", "rr")):
        arms, fits, pbest, pmean = [], [], [], []
        for s in seeds:
            rng = np.random.default_rng(SEED + 1000 * s + hash(t["pdb"]) % 997)
            if get == "bf":
                idx = blos[rng.permutation(len(blos))[:M]]
                pool = blos
            else:
                pool = rng.permutation(nw)[:min(K, nw)]
                idx = top_m(pool) if get == "rd" else pool[rng.permutation(len(pool))[:M]]
            a, f = emit(idx)
            arms.append(a); fits.append(f)
            pbest.append(float(u["rr"][pool].min())); pmean.append(float(u["rr"][pool].mean()))
        res[name + "_arm"] = float(np.mean(arms)); res[name + "_arm_sd"] = float(np.std(arms))
        res[name + "_fit"] = float(np.mean(fits))
        res[name + "_pool_best"] = float(np.mean(pbest)); res[name + "_pool_mean"] = float(np.mean(pmean))
    return res


def main(n_seeds=3):
    seeds = list(range(n_seeds))
    tg = I.targets()
    out_path = os.path.join(I.RESULTS, "coord_null.json")
    rows = []
    if os.path.exists(out_path):
        rows = json.load(open(out_path)).get("per_target", [])
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        rows.append(one_target(t, seeds))
        if k % 5 == 0 or k == len(tg) - 1:
            _dump(rows, seeds, out_path)
            print(f"  {k + 1}/{len(tg)} {t['pdb']} shipped={rows[-1]['shipped_arm']:.3f} "
                  f"blind={rows[-1]['blind_arm']:.3f}  [{time.time() - t0:.0f}s]", flush=True)
    _dump(rows, seeds, out_path)
    _report(rows)


def _dump(rows, seeds, path):
    with open(path, "w") as fh:
        json.dump({"what": "sequence-information ablation of the whole architecture",
                   "n_seeds": len(seeds), "per_target": rows}, fh, indent=1)


def _report(rows):
    names = [r["pdb"] for r in rows]; folds = [r["fold"] for r in rows]
    f18 = set(I.FAIL18); isf = np.array([r["pdb"] in f18 for r in rows])
    arms = {k: np.array([r[k] for r in rows], float)
            for k in ("shipped_sel", "shipped_arm", "no_filter_arm", "no_retrieval_arm", "blind_arm")}
    tab = {k: I.summary(v) for k, v in arms.items()}
    cmp = {
        "distogram_value  (blind -> no_retrieval)": I.paired(arms["no_retrieval_arm"], arms["blind_arm"], folds=folds, names=names),
        "blosum_value     (blind -> no_filter)":    I.paired(arms["no_filter_arm"], arms["blind_arm"], folds=folds, names=names),
        "both             (blind -> shipped)":      I.paired(arms["shipped_arm"], arms["blind_arm"], folds=folds, names=names),
        "distogram_given_blosum (no_filter -> shipped)": I.paired(arms["shipped_arm"], arms["no_filter_arm"], folds=folds, names=names),
        "blosum_given_distogram (no_retrieval -> shipped)": I.paired(arms["shipped_arm"], arms["no_retrieval_arm"], folds=folds, names=names),
    }
    sub = {"FAIL18": {k: float(v[isf].mean()) for k, v in arms.items()},
           "other108": {k: float(v[~isf].mean()) for k, v in arms.items()}}
    out = {"what": "How much of the emitted structure is sequence information?",
           "arms": tab, "paired": cmp, "subgroups": sub}
    I.write("coord_null_report", out)
    print("\n=== emitted CA-RMSD (lam=0.3 arm, no AMBER) ===")
    for k, v in tab.items():
        print(f"  {k:22s} {v['mean']:.4f}   median {v['median']:.3f}  <2A {v['frac_under_2.0']:.3f}")
    print("\n=== paired (negative = first arm better) ===")
    for k, v in cmp.items():
        print(f"  {k:48s} {v['mean_diff']:+.4f} [{v['ci95'][0]:+.4f},{v['ci95'][1]:+.4f}] "
              f"{v['n_better']}W/{v['n_worse']}L  drop10 {v['drop_top10_mean_diff']:+.4f}")
    print("\n=== subgroups ===")
    for g, d in sub.items():
        print(f"  {g:9s} " + "  ".join(f"{k.replace('_arm',''):13s}{v:.3f}" for k, v in d.items()))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
