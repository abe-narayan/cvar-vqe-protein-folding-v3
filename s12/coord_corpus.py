"""COORDINATOR EXPERIMENT 4 -- the retrieval library is 80% drawn from a corpus in which
sequence carries almost no structural information.

`s12/coord_infochannel.py` measured the sequence->structure channel model-free, by pairing
every same-length member of each corpus against every other and correlating sequence
identity with CA-RMSD:

    corpus                            rho(identity, RMSD), n = 9..16
    peptide database (787 chains)     -0.17 .. -0.29     (all p < 1e-5)
    protein fragments (6,003)         -0.030 .. -0.040

The peptide corpus carries roughly SEVEN TIMES the sequence->structure correlation of the
protein-fragment corpus.  And inside a target's own retrieval universe -- which is ~80%
fragments -- the channel measures rho(BLOSUM sim, true RMSD) = **-0.066**, exactly the
mixture those two numbers predict, and **+0.001 on the FAIL18**, i.e. nothing at all.

The mechanism is not mysterious.  A 12-residue stretch cut out of a folded protein has the
conformation its tertiary context imposes; its own sequence did not choose it.  An isolated
peptide's conformation is all its sequence has.  `fragment_db`'s own docstring says as
much and the record already knows the training-distribution version of this (S7-2: adding
protein-fragment pairs made the distogram monotonically WORSE).  What has not been asked is
the RETRIEVAL version: fragments are 80% of the pool the filter has to rank.

The forensics agent adds the piece that makes this actionable: on the FAIL18 the shipped
top-75 UNDER-selects peptide-derived windows (0.177 of the set) while the near-native band
is ENRICHED in them (0.389).  The filter is steering toward the uninformative corpus
exactly where the informative one holds the answer.

WHY THIS IS NOT A REPEAT OF A CLOSED EXPERIMENT.  `s8/generate.py`'s `gen_pep` arm measured
a peptide-only pool and found it null: selected 3.4451 against base 3.4540, mean_diff
-0.0089 [-0.098, +0.080].  **That was measured on the ARGMIN terminal operator.**  The same
file records the peptide-only pool's MEAN at 4.321 against base 4.453 -- a better pool
distribution -- and S7-6/S8-11 establish that the synthesis operator responds to the pool
MEAN while argmin responds to the pool BEST.  So the arm was measured against the one
terminal operator least able to use what it improved.  This module re-measures it on the
operator the system actually ships.

ARMS (all K=500, then shipped-score top-75, coordinate average, projection at lambda=0.3):
    base        the shipped BLOSUM pool                                     (incumbent)
    pep         BLOSUM top-500 restricted to peptide-derived windows
    frag        BLOSUM top-500 restricted to fragment-derived windows       (the other half)
    boost       BLOSUM + a bonus on peptide provenance, bonus swept

CONTROL.  `pep` and `base` differ in SELECTIVITY as well as provenance: the peptide sub-
universe is ~20% the size, so a top-500 out of it is a much less selective cut, and less
selective cuts have better means for reasons that have nothing to do with provenance.  The
`selctl` arm fixes that: a RANDOM sub-universe of the same size as the peptide one, cut to
top-500 by the same BLOSUM key.  A `pep` gain that `selctl` reproduces is a selectivity
artefact.

    python -m s12.coord_corpus
"""
import os
import sys
import json
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402

BOOSTS = (2.0, 5.0, 10.0)


def run_target(t, seeds=(0, 1)):
    u = I.load_univ(t["pdb"])
    n, seq, fold, nat = t["n"], t["seq"], t["fold"], u["nat_ca"]
    W_all = u["W"]; rr = u["rr"]; sim = np.asarray(u["sim"], float)
    org = np.asarray(u["org"], bool)
    dg = I.distogram(t["pdb"], seq, fold)
    i, j = I.pair_index(n)

    def pool_from(order):
        return np.asarray(order[:I.K], int)

    def chain(pool):
        Wp = W_all[pool].astype(np.float32).astype(float)
        sc = I.shipped_score(dg, I.pair_dists(Wp, i, j))
        sub = np.argsort(sc, kind="stable")[:I.M]
        C, _ = I.coordinate_average(Wp[sub])
        o = I.project(C, seq, fold, lam=0.3)
        return {"arm": I.ca_rmsd(o["ca"], nat), "fit": I.ca_rmsd(o["fit_ca"], nat),
                "argmin": float(rr[pool][int(np.argmin(sc))]),
                "pool_best": float(rr[pool].min()), "pool_mean": float(rr[pool].mean()),
                "top75_best": float(rr[pool][sub].min()),
                "top75_mean": float(rr[pool][sub].mean()),
                "frac_pep_pool": float(org[pool].mean()),
                "frac_pep_top75": float(org[pool][sub].mean())}

    res = {"pdb": t["pdb"], "n": n, "fold": fold, "n_windows": int(len(rr)),
           "frac_pep_universe": float(org.mean())}

    # incumbent: the pinned stable argsort of -sim over the whole universe
    res["base"] = chain(pool_from(np.asarray(u["order"], int)))

    for name, mask in (("pep", org), ("frag", ~org)):
        idx = np.where(mask)[0]
        if len(idx) < I.K:
            res[name] = None
            continue
        sub_order = idx[np.argsort(-sim[idx], kind="stable")]
        res[name] = chain(pool_from(sub_order))
        res[name]["sub_universe"] = int(len(idx))

    # provenance BOOST: one key, one parameter, base at bonus 0
    for b in BOOSTS:
        key = sim + b * org.astype(float)
        res[f"boost{b}"] = chain(pool_from(np.argsort(-key, kind="stable")))

    # SELECTIVITY control: a random sub-universe the size of the peptide one
    npep = int(org.sum())
    if npep >= I.K:
        arms = []
        for s in seeds:
            rng = np.random.default_rng(7717 * s + hash(t["pdb"]) % 9973)
            idx = rng.permutation(len(rr))[:npep]
            sub_order = idx[np.argsort(-sim[idx], kind="stable")]
            arms.append(chain(pool_from(sub_order)))
        res["selctl"] = {k: float(np.mean([a[k] for a in arms])) for k in arms[0]}
    else:
        res["selctl"] = None
    return res


def main():
    tg = I.targets()
    path = os.path.join(I.RESULTS, "coord_corpus.json")
    rows = json.load(open(path))["per_target"] if os.path.exists(path) else []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        rows.append(run_target(t))
        if k % 5 == 0 or k == len(tg) - 1:
            json.dump({"what": "corpus provenance in retrieval", "per_target": rows},
                      open(path, "w"), indent=1)
            r = rows[-1]
            pep = r["pep"]["arm"] if r["pep"] else float("nan")
            print(f"  {k+1}/{len(tg)} {t['pdb']} base={r['base']['arm']:.2f} pep={pep:.2f} "
                  f"boost5={r['boost5.0']['arm']:.2f} pepfrac={r['frac_pep_universe']:.2f} "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    json.dump({"what": "corpus provenance in retrieval", "per_target": rows},
              open(path, "w"), indent=1)
    report(rows)


def report(rows):
    arms = ["base", "pep", "frag", "selctl"] + [f"boost{b}" for b in BOOSTS]
    ok = [r for r in rows if all(r.get(a) for a in arms)]
    print(f"\n{len(ok)}/{len(rows)} targets have every arm")
    names = [r["pdb"] for r in ok]; folds = [r["fold"] for r in ok]
    f18 = set(I.FAIL18); isf = np.array([r["pdb"] in f18 for r in ok])
    g = lambda a, k: np.array([r[a][k] for r in ok], float)                 # noqa: E731
    out = {"n": len(ok), "arms": {}, "paired_vs_base": {}, "subgroups": {}}
    for a in arms:
        out["arms"][a] = {k: float(g(a, k).mean()) for k in
                          ("arm", "fit", "argmin", "pool_best", "pool_mean", "top75_best",
                           "top75_mean", "frac_pep_pool", "frac_pep_top75")}
        out["arms"][a]["arm_summary"] = I.summary(g(a, "arm"))
        if a != "base":
            out["paired_vs_base"][a] = I.paired(g(a, "arm"), g("base", "arm"),
                                                folds=folds, names=names)
        out["subgroups"][a] = {"FAIL18": float(g(a, "arm")[isf].mean()),
                               "other108": float(g(a, "arm")[~isf].mean())}
    out["paired_pep_vs_selctl"] = I.paired(g("pep", "arm"), g("selctl", "arm"),
                                           folds=folds, names=names)
    I.write("coord_corpus_report", out)
    print(f"\n{'arm':10s} {'emitted':>8s} {'argmin':>8s} {'poolbest':>9s} {'poolmean':>9s} "
          f"{'t75best':>8s} {'pep%pool':>9s} {'pep%75':>7s}")
    for a in arms:
        v = out["arms"][a]
        print(f"{a:10s} {v['arm']:8.4f} {v['argmin']:8.4f} {v['pool_best']:9.4f} "
              f"{v['pool_mean']:9.4f} {v['top75_best']:8.4f} {v['frac_pep_pool']:9.3f} "
              f"{v['frac_pep_top75']:7.3f}")
    print("\npaired vs base (negative = better):")
    for a, v in out["paired_vs_base"].items():
        print(f"  {a:10s} {v['mean_diff']:+.4f} [{v['ci95'][0]:+.4f},{v['ci95'][1]:+.4f}] "
              f"{v['n_better']}W/{v['n_worse']}L  drop10 {v['drop_top10_mean_diff']:+.4f}  "
              f"folds {['%+.2f' % x for x in v['per_fold'].values()]}")
    v = out["paired_pep_vs_selctl"]
    print(f"\nSELECTIVITY CONTROL  pep vs selctl: {v['mean_diff']:+.4f} "
          f"[{v['ci95'][0]:+.4f},{v['ci95'][1]:+.4f}] {v['n_better']}W/{v['n_worse']}L")
    print("\nsubgroups:")
    for a in arms:
        s = out["subgroups"][a]
        print(f"  {a:10s} FAIL18 {s['FAIL18']:.4f}   other108 {s['other108']:.4f}")


if __name__ == "__main__":
    main()
