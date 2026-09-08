"""ADVERSARIAL AUDIT 2 -- is the FAIL18 a stable object, or a threshold artefact?

FAIL18 is defined by ZERO RECALL: none of the K=500 pool members within
(pool_best + BAND=1.5 A) survives the top-M=75 filter.  Three arbitrary constants.  Every
"FAIL18 vs other-108" contrast in the sprint inherits whatever instability they carry.

  (a) BAND sweep      0.75 ... 3.0 A, and an ABSOLUTE band (rr <= 2.5 A) which does not
                      reference pool_best at all
  (b) POOL SIZE K     250 / 500 / 1000  (changes both the pool and pool_best)
  (c) FILTER SIZE M   50 / 75 / 100
  (d) CIRCULARITY     zero recall <=> band_score_pct_min > M/K = 0.15 BY DEFINITION.  So
                      the reported 0.391 vs 0.009 contrast is partly definitional.  The
                      honest control is the LOW-RECALL group (1-3 band members kept), which
                      is NOT in FAIL18 but is subject to the same near-definitional force.
  (e) SHRINKAGE       recompute the headline subgroup contrasts (emitted RMSD gap, rho gap,
                      fibril/lasso Fisher) under every perturbed membership.

There is no random seed anywhere in this chain -- retrieval is a deterministic stable
argsort and the filter is a deterministic gather -- so there is nothing to reseed; the only
instability available is in the three constants.  That is stated, not assumed: verified by
recomputing the whole selfcheck twice.
"""
from __future__ import annotations
import os, sys, json, re, itertools
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import instrument as I

BANDS = [0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0]
KS = [250, 500, 1000]
MS = [50, 75, 100]


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(len(a | b), 1)


def main():
    tg = I.targets()
    per = {}
    for t in tg:
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        dg = I.distogram(pdb, t["seq"], t["fold"])
        i, j = I.pair_index(t["n"])
        rec = I.shipped_record(pdb)
        d = dict(n=t["n"], fold=t["fold"], emitted=float(rec["rmsd_fit"]))
        for k in KS:
            p = I.pool_idx(u, k)
            W = u["W"][p]
            rr = u["rr"][p]
            sc = I.shipped_score(dg, I.pair_dists(W, i, j))
            order = np.argsort(sc, kind="stable")
            d[f"K{k}"] = dict(rr=rr, order=order, sc=sc)
        per[pdb] = d
        print(pdb, flush=True)

    names = [t["pdb"] for t in tg]
    # ---------- membership under every (band, K, M)
    sets = {}
    for band, k, m in itertools.product(BANDS, KS, MS):
        z = []
        for pdb in names:
            g = per[pdb][f"K{k}"]
            rr = g["rr"]; sub = g["order"][:m]
            b = np.where(rr <= rr.min() + band)[0]
            if not np.isin(b, sub).any():
                z.append(pdb)
        sets[f"band{band}_K{k}_M{m}"] = z
    # absolute band, no reference to pool_best
    for thr, k, m in itertools.product([2.0, 2.5, 3.0], KS, MS):
        z = []
        for pdb in names:
            g = per[pdb][f"K{k}"]
            rr = g["rr"]; sub = g["order"][:m]
            b = np.where(rr <= thr)[0]
            if len(b) and not np.isin(b, sub).any():
                z.append(pdb)
        sets[f"abs{thr}_K{k}_M{m}"] = z

    base = set(I.FAIL18)
    stab = {kk: dict(n=len(v), jaccard=jaccard(v, base),
                     kept=len(set(v) & base), added=len(set(v) - base),
                     lost=sorted(base - set(v)), extra=sorted(set(v) - base))
            for kk, v in sets.items()}

    # per-target membership frequency across the whole grid
    freq = {p: float(np.mean([p in v for v in sets.values()])) for p in names}

    # ---------- (e) shrinkage of the headline contrasts under each membership
    emitted = np.array([per[p]["emitted"] for p in names])
    rho = np.array([spearmanr(per[p]["K500"]["sc"], per[p]["K500"]["rr"]).statistic for p in names])
    with open(os.path.join(I.RESULTS, "adv_meta.json")) as fh:
        am = json.load(fh)
    # rebuild the fibril|lasso flag from my own audit (adv_meta stored counts only), so
    # re-derive it here from fail_headers.json which I already cross-checked
    with open(os.path.join(I.RESULTS, "fail_headers.json")) as fh:
        hd = json.load(fh)["per_target"]
    lf = np.array([bool(hd[p]["fibril"]) or bool(hd[p]["lasso"]) for p in names])

    def contrasts(members):
        f = np.array([p in set(members) for p in names], bool)
        if f.sum() < 3 or (~f).sum() < 3:
            return None
        a = int((lf & f).sum()); b = int((lf & ~f).sum())
        return dict(n_fail=int(f.sum()),
                    emitted_fail=float(emitted[f].mean()), emitted_other=float(emitted[~f].mean()),
                    emitted_gap=float(emitted[f].mean() - emitted[~f].mean()),
                    rho_fail=float(rho[f].mean()), rho_other=float(rho[~f].mean()),
                    rho_gap=float(rho[f].mean() - rho[~f].mean()),
                    lf_fail=a, lf_other=b,
                    lf_p=float(fisher_exact([[a, int(f.sum()) - a],
                                             [b, int((~f).sum()) - b]])[1]))

    shrink = {kk: contrasts(v) for kk, v in sets.items()}

    # ---------- (d) circularity control: the LOW-RECALL group
    recall = {}
    for p in names:
        g = per[p]["K500"]; rr = g["rr"]; sub = g["order"][:75]
        b = np.where(rr <= rr.min() + 1.5)[0]
        recall[p] = int(np.isin(b, sub).sum())
    low = [p for p in names if 1 <= recall[p] <= 3]
    hi = [p for p in names if recall[p] > 3]
    zero = [p for p in names if recall[p] == 0]

    def bandpct(p, m=75):
        g = per[p]["K500"]; rr = g["rr"]; sc = g["sc"]
        b = np.where(rr <= rr.min() + 1.5)[0]
        pct = np.array([(sc < sc[x]).mean() for x in b])
        return float(pct.mean()), float(pct.min())

    bp = {p: bandpct(p) for p in names}
    grp = {}
    for nm, gg in (("zero(FAIL18)", zero), ("low_recall_1_3", low), ("high_recall_4+", hi)):
        grp[nm] = dict(n=len(gg),
                       band_pct_mean=float(np.mean([bp[p][0] for p in gg])),
                       band_pct_min=float(np.mean([bp[p][1] for p in gg])),
                       emitted=float(np.mean([per[p]["emitted"] for p in gg])),
                       rho=float(np.mean([rho[names.index(p)] for p in gg])),
                       lf_frac=float(np.mean([lf[names.index(p)] for p in gg])))

    out = dict(base=sorted(base), stability=stab, freq=freq, shrinkage=shrink,
               recall_groups=grp, recall=recall,
               summary=dict(
                   n_grid=len(sets),
                   jaccard_mean=float(np.mean([v["jaccard"] for v in stab.values()])),
                   jaccard_min=float(min(v["jaccard"] for v in stab.values())),
                   always_in=[p for p in names if freq[p] == 1.0],
                   never_in=[p for p in names if freq[p] == 0.0],
                   core_18_always=[p for p in sorted(base) if freq[p] == 1.0],
                   base_freq={p: freq[p] for p in sorted(base)},
               ))
    print(json.dumps(out["summary"], indent=1))
    print(json.dumps(grp, indent=1))
    I.write("adv_fail18", out)


if __name__ == "__main__":
    main()
