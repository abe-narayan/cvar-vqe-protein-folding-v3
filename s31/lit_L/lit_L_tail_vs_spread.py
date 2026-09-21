"""Does the deposited-ensemble spread explain the FAIL18 tail?

ORACLE DIAGNOSTIC of the BENCHMARK, not a predictor.
"""
import csv, json, os, sys, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = r"C:\Users\abena\Protein-Folding-Algorithm"

rows = list(csv.DictReader(open(os.path.join(ROOT, "results/summary/results.csv"))))
cfgs = collections.Counter(r["configuration"] for r in rows)
print("configurations:", dict(cfgs))

spread = {r["pdb"]: r for r in json.load(open(os.path.join(HERE, "ens_spread.json")))}


def arm(cfg):
    d = {}
    for r in rows:
        if r["configuration"] != cfg or r["status"] != "present":
            continue
        try:
            d[r["pdb_id"].upper()] = float(r["rmsd"])
        except (ValueError, TypeError):
            pass
    return d


for cfg in cfgs:
    d = arm(cfg)
    if len(d) < 100:
        print(f"  skip {cfg}: n={len(d)}")
        continue
    v = np.array(list(d.values()))
    print(f"\n=== {cfg}: n={len(d)} mean={v.mean():.4f} median={np.median(v):.4f}")

    pairs = [(p, r, spread[p]["spread"], spread[p]["m1"])
             for p, r in d.items()
             if p in spread and spread[p]["spread"] == spread[p]["spread"]]
    if len(pairs) < 50:
        print("   too few matched"); continue
    P_ = np.array([x[1] for x in pairs])
    S_ = np.array([x[2] for x in pairs])
    M_ = np.array([x[3] for x in pairs])
    n = len(pairs)

    def r_and_ci(a, b):
        r = float(np.corrcoef(a, b)[0, 1])
        z = np.arctanh(r); se = 1 / np.sqrt(len(a) - 3)
        return r, float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))

    r1 = r_and_ci(P_, S_); r2 = r_and_ci(P_, M_)
    print(f"   matched n={n}")
    print(f"   corr(production RMSD, ensemble spread)   = {r1[0]:+.4f}  95% CI [{r1[1]:+.4f},{r1[2]:+.4f}]")
    print(f"   corr(production RMSD, model1->medoid)    = {r2[0]:+.4f}  95% CI [{r2[1]:+.4f},{r2[2]:+.4f}]")

    o = np.argsort(-P_)
    k = 18
    worst, rest = o[:k], o[k:]
    print(f"   WORST {k} by production RMSD : mean RMSD {P_[worst].mean():.4f}"
          f"   mean ensemble spread {S_[worst].mean():.4f}   m1->medoid {M_[worst].mean():.4f}")
    print(f"   the other {n-k}              : mean RMSD {P_[rest].mean():.4f}"
          f"   mean ensemble spread {S_[rest].mean():.4f}   m1->medoid {M_[rest].mean():.4f}")
    ds = S_[worst].mean() - S_[rest].mean()
    se = np.sqrt(S_[worst].var(ddof=1)/k + S_[rest].var(ddof=1)/len(rest))
    print(f"   spread difference tail-minus-rest = {ds:+.4f}  SE {se:.4f}  "
          f"MDE {2.8016*se:.4f}  -> {abs(ds)/(2.8016*se):.2f}x MDE")

    hi = S_ > 2.0
    print(f"   targets with deposited spread > 2.0 A: n={hi.sum()}  "
          f"mean production RMSD {P_[hi].mean():.4f}  vs {P_[~hi].mean():.4f} for the rest")
