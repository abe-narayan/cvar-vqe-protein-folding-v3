"""s22/results/d_rg_route.py -- THE ANGSTROM CONVERSION TEST.
Pre-registered in s22/PREREG_D.md (D6) BEFORE these numbers were read.
Does rg_z / rg_gap (L27/L28's compactness-disagreement channel), used as a NATIVE-FREE
per-target router feature, let an ACHIEVABLE (held-out-fold) selector beat the fixed
incumbent (avg_75, 3.048) on the FULL n=126 pool?  Candidate router targets: choose between
avg_75 (trust the score) and each of the wider/narrower ladder rungs + medoid, since L19/L27's
mechanism (extended/disagreement targets are a CONCENTRATED WRONG REGION for the score) predicts
that on large-|rg_z| targets a selection-lighter arm (wider average, or medoid) should do
relatively better than the score's top-75.
"""
import json
import numpy as np

rg = {r["pdb"]: r for r in json.load(open("s22/results/d_rg_full126.json"))["rows"]}
pg = {r["pdb"]: r for r in json.load(open("s21/results/poolgap.json"))["rows"]}
common = sorted(set(rg) & set(pg))
assert len(common) == 126

ARMS = ["avg_500", "avg_150", "avg_75", "avg_20", "avg_5", "avg_1", "medoid_all", "medoid_75"]
X = np.array([[pg[p][a] for a in ARMS] for p in common])
inc_col = ARMS.index("avg_75")
folds = np.array([pg[p]["fold"] for p in common])
rgz = np.array([rg[p]["rg_z"] for p in common])
rggap = np.array([rg[p]["rg_gap"] for p in common])
n = np.array([pg[p]["n"] for p in common])

def stat(x):
    x = np.asarray(x, float); k = len(x)
    se = x.std(ddof=1)/np.sqrt(k)
    return x.mean(), se, 2.8016*se

print(f"n = {len(common)} targets.  Oracle over the 8-arm ladder = {X.min(1).mean():.4f}  "
      f"(fixed avg_75 = {X[:,inc_col].mean():.4f})")
print(f"corr(rg_z, n) = {np.corrcoef(rgz, n)[0,1]:+.4f}   corr(rg_gap, n) = {np.corrcoef(rggap,n)[0,1]:+.4f}")

# in-sample oracle check: is there ANY monotone relationship between |rg_z| and which arm wins?
winner = X.argmin(1)
print("\nmean |rg_z| by winning arm (in-sample, diagnostic only, NOT a router):")
for k, a in enumerate(ARMS):
    sel = winner == k
    if sel.sum() == 0: continue
    print(f"  {a:12s} n={sel.sum():3d}  mean rg_z={rgz[sel].mean():+.3f}  mean rg_gap={rggap[sel].mean():+.3f}")

def held_out_router(feature, cand_arms, k_bins=3):
    """For each fold: on training folds, bin `feature` into k_bins quantile bins; within each
    bin pick whichever of cand_arms has the lowest TRAINING mean; apply to held-out fold."""
    achieved = np.empty(len(common))
    for f in np.unique(folds):
        te = np.where(folds == f)[0]
        tr = np.where(folds != f)[0]
        edges = np.quantile(feature[tr], np.linspace(0, 1, k_bins+1))
        edges[0] -= 1e-9; edges[-1] += 1e-9
        bin_tr = np.clip(np.digitize(feature[tr], edges)-1, 0, k_bins-1)
        best = {}
        for b in range(k_bins):
            idx = tr[bin_tr == b]
            if len(idx) == 0: continue
            means = X[idx][:, [ARMS.index(a) for a in cand_arms]].mean(0)
            best[b] = cand_arms[int(np.argmin(means))]
        bin_te = np.clip(np.digitize(feature[te], edges)-1, 0, k_bins-1)
        for t, b in zip(te, bin_te):
            arm = best.get(b, "avg_75")
            achieved[t] = X[t, ARMS.index(arm)]
    return achieved

fixed = X[:, inc_col]
for feat_name, feat in (("rg_z", rgz), ("rg_gap", rggap), ("|rg_z|", np.abs(rgz))):
    for cand in (["avg_75", "avg_500"], ["avg_75", "avg_150"], ["avg_75", "medoid_75"],
                 ARMS):
        ach = held_out_router(feat, cand, k_bins=3)
        d = ach - fixed
        m, se, mde = stat(d)
        print(f"  router[{feat_name:7s} -> {str(cand):45s}]  achieved {ach.mean():.4f}  "
              f"vs fixed {fixed.mean():.4f}  d={m:+.4f} SE {se:.4f} MDE {mde:.4f}  "
              f"W/L {(d<0).sum()}/{(d>0).sum()}")
