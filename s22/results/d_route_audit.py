"""s22/results/d_route_audit.py -- WORKSTREAM D independent audit/reconstruction of the
coordinator's s22 LEDGER L1 (the 0.482 A routing headroom over 13 native-free arms).

Run: python -m s22.results.d_route_audit   (or `python s22/results/d_route_audit.py` from repo root)
"""
import json, itertools
import numpy as np

POOLGAP = json.load(open("s21/results/poolgap.json"))["rows"]
LATSEL  = json.load(open("s21/results/latentsel.json"))["rows"]

pg = {r["pdb"]: r for r in POOLGAP}
ls = {r["pdb"]: r for r in LATSEL}
common = sorted(set(pg) & set(ls))
print(f"poolgap targets: {len(pg)}  latentsel targets: {len(ls)}  common: {len(common)}")
assert len(common) == 126, "expected the full 126-target panel in both files"

# ---- 0. Sanity: identity/duplicate columns inside poolgap.json --------------------------
diffs_a1_pa = np.array([pg[p]["avg_1"] - pg[p]["pool_argmin"] for p in common])
diffs_a75_m = np.array([pg[p]["avg_75"] - pg[p]["avg75_medoid75"] for p in common])
print(f"\n[duplicate check] avg_1 vs pool_argmin (poolgap):      max|diff| = {np.abs(diffs_a1_pa).max():.2e}")
print(f"[duplicate check] avg_75 vs avg75_medoid75 (poolgap):   max|diff| = {np.abs(diffs_a75_m).max():.2e}")
print("  -> both are IDENTITY columns inside poolgap.json, not independent arms.")

# ---- 1. Basis price: window (poolgap) vs rebuilt-onto-latent-manifold (latentsel) --------
d_argmin_basis = np.array([ls[p]["pool_argmin"] - pg[p]["avg_1"] for p in common])
d_avg75_basis  = np.array([ls[p]["pool_avg75"]  - pg[p]["avg_75"] for p in common])
d_ship_check   = np.array([ls[p]["ship_avg75"]  - pg[p]["avg_75"] for p in common])
def stat(x):
    x = np.asarray(x, float); n = len(x)
    se = x.std(ddof=1) / np.sqrt(n)
    return x.mean(), se, 2.8016*se

for name, d in (("pool_argmin: rebuilt - window", d_argmin_basis),
                ("pool_avg75:  rebuilt - window", d_avg75_basis),
                ("ship_avg75 vs poolgap avg_75 (both window, should be ~0)", d_ship_check)):
    m, se, mde = stat(d)
    print(f"\n[basis price] {name}")
    print(f"    mean {m:+.4f}  SE {se:.4f}  MDE {mde:.4f}  max|.|={np.abs(d).max():.4f}  "
          f"n_nonzero={int((np.abs(d)>1e-9).sum())}/{len(d)}")


# ---- 2. Candidate 13-arm reconstructions -------------------------------------------------
POOLGAP_LADDER = ["avg_500", "avg_150", "avg_75", "avg_20", "avg_5", "avg_1"]   # 6, window basis
POOLGAP_MEDOID = ["medoid_all", "medoid_75"]                                    # 2, window basis
LAT_ARMS       = ["lat_rand1", "lat_argmin", "lat_avg75", "lat_medoid"]         # 4, built-chain/latent

def arm_value(name, p):
    if name in pg[p]:
        return pg[p][name]
    return ls[p][name]

def build_table(names):
    X = np.array([[arm_value(nm, p) for nm in names] for p in common])
    return X

def route_report(names, label):
    X = build_table(names)                      # (126, K)
    means = X.mean(0)
    best_i = int(np.argmin(means))
    best_fixed_name, best_fixed = names[best_i], means[best_i]
    per_target_oracle = X.min(1)
    oracle_mean = per_target_oracle.mean()
    headroom = best_fixed - oracle_mean
    winner = X.argmin(1)
    counts = np.bincount(winner, minlength=len(names))
    incumbent_idx = names.index("avg_75") if "avg_75" in names else None
    inc_rmsd = np.array([pg[p]["avg_75"] for p in common])
    per_target_headroom = np.array([pg[p]["avg_75"] for p in common]) - per_target_oracle
    order = np.argsort(-per_target_headroom)
    cum = np.cumsum(per_target_headroom[order])
    total = cum[-1]
    top21 = cum[20] / total if total > 0 else float("nan")
    top50 = cum[49] / total if total > 0 else float("nan")
    corr = np.corrcoef(per_target_headroom, inc_rmsd)[0, 1]
    n_inc_best = int((per_target_headroom <= 1e-9).sum())
    print(f"\n=== {label}  (K={len(names)} arms: {names}) ===")
    print(f"  best fixed arm: {best_fixed_name} = {best_fixed:.4f}")
    print(f"  oracle per-target routing (over these K arms) = {oracle_mean:.4f}")
    print(f"  HEADROOM (best_fixed - oracle) = {best_fixed - oracle_mean:.4f}")
    print(f"  headroom vs avg_75 as fixed baseline           = {inc_rmsd.mean() - oracle_mean:.4f}")
    print(f"  winner counts: {dict(zip(names, counts.tolist()))}")
    print(f"  concentration (baseline=avg_75): top21/126={top21*100:.1f}%  top50/126={top50*100:.1f}%  "
          f"median per-target headroom={np.median(per_target_headroom):.3f}  "
          f"corr(headroom, incumbent RMSD)={corr:+.3f}  n targets where incumbent already ORACLE-best={n_inc_best}")
    return dict(best_fixed=best_fixed, oracle=oracle_mean, headroom=best_fixed-oracle_mean)

# reconstruction A: literal reading -- 6 + 2 + 4 = 12 arms (as coordinator's prose parses)
A = POOLGAP_LADDER + POOLGAP_MEDOID + LAT_ARMS
route_report(A, "RECON A: m-ladder(6) + medoids(2) + 4 latent arms = 12")

# reconstruction B: A + pool_argmin[rebuilt onto latent manifold] (13, keeps window pool_argmin AND rebuild)
B = A + ["pool_argmin"]
route_report(B, "RECON B: A + latentsel.pool_argmin (rebuilt) = 13")

# reconstruction C: A + pool_avg75[rebuilt] instead
C = A + ["pool_avg75"]
route_report(C, "RECON C: A + latentsel.pool_avg75 (rebuilt) = 13")

# reconstruction D: A + BOTH rebuilt pool arms = 14 (one too many, shown for bracketing)
D = A + ["pool_argmin", "pool_avg75"]
route_report(D, "RECON D: A + BOTH rebuilt pool arms = 14 (bracket)")

# reconstruction E: drop the window avg_1/pool_argmin duplicate pair issue -- use ladder without
# avg_1 (since avg_1 IS the window pool_argmin) and add rebuilt pool_argmin + rebuilt pool_avg75
# + all 4 latent arms = (avg_500,150,75,20,5)=5 + medoids(2) + rebuilt pool_argmin,pool_avg75(2) + lat(4) = 13
E = ["avg_500","avg_150","avg_75","avg_20","avg_5"] + POOLGAP_MEDOID + ["pool_argmin","pool_avg75"] + LAT_ARMS
route_report(E, "RECON E: ladder-without-avg_1(5) + medoids(2) + rebuilt pool argmin+avg(2) + lat(4) = 13")


# ---- 3. THE ACHIEVABLE-ROUTER TEST: length-bucket router, proper 5-fold held-out CV -------
# Zero-cost feature "n" (peptide length) is already present in both source jsons -- no new
# heavy computation (no distogram/AMBER call) is needed for this minimal, honest first pass.
print("\n" + "="*90)
print("THE ACHIEVABLE ROUTER: can ANY router beat the best fixed arm (3.048) under proper")
print("held-out-fold scoring?  Router = 'best arm on TRAINING folds, within this target's own")
print("length-quartile (quartile boundaries fit on TRAINING folds only)'.")
print("="*90)

NAMES = A  # the 12-arm reconstruction; result is insensitive to 12 vs 13 vs 14 per above
X = build_table(NAMES)                       # (126, K)
ns = np.array([pg[p]["n"] for p in common])
folds = np.array([pg[p]["fold"] for p in common])
inc_col = NAMES.index("avg_75")

def route_predict(train_idx, test_idx, X, ns, k_bins=4):
    """Fit length-quantile bin edges + per-bin best arm on TRAIN, apply to TEST."""
    edges = np.quantile(ns[train_idx], np.linspace(0, 1, k_bins + 1))
    edges[0] -= 1e-6; edges[-1] += 1e-6
    bin_of_train = np.digitize(ns[train_idx], edges) - 1
    best_arm_per_bin = {}
    for b in range(k_bins):
        sel = train_idx[bin_of_train == b]
        if len(sel) == 0:
            continue
        best_arm_per_bin[b] = int(np.argmin(X[sel].mean(0)))
    bin_of_test = np.clip(np.digitize(ns[test_idx], edges) - 1, 0, k_bins - 1)
    pred = np.array([X[t, best_arm_per_bin.get(b, inc_col)] for t, b in zip(test_idx, bin_of_test)])
    return pred

achieved = np.empty(len(common))
for f in np.unique(folds):
    test_idx = np.where(folds == f)[0]
    train_idx = np.where(folds != f)[0]
    achieved[test_idx] = route_predict(train_idx, test_idx, X, ns, k_bins=4)

fixed_baseline = X[:, inc_col]
d = achieved - fixed_baseline
m, se, mde = stat(d)
lo = m - 1.96*se; hi = m + 1.96*se
print(f"\nlength-quartile router (5-fold held-out) mean RMSD = {achieved.mean():.4f}")
print(f"best FIXED arm (avg_75, in-sample)                  = {fixed_baseline.mean():.4f}")
print(f"achieved - fixed  = {m:+.4f}  SE {se:.4f}  MDE {mde:.4f}  approx 95% CI [{lo:+.4f},{hi:+.4f}]")
print(f"W/L (achieved beats fixed per target): {(d<0).sum()}/{(d>0).sum()}")
print(f"\nfor reference, the ORACLE ceiling on these same 12 arms = {X.min(1).mean():.4f}  "
      f"(fraction of the {fixed_baseline.mean()-X.min(1).mean():.4f} A ceiling captured: "
      f"{100*(fixed_baseline.mean()-achieved.mean())/(fixed_baseline.mean()-X.min(1).mean()):.1f}%)")
