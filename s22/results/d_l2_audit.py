"""s22/results/d_l2_audit.py -- audit of s22 LEDGER L2 (the retrieval/latent difficulty
interaction) for a regression-to-the-mean artefact, per the coordinator's own objection."""
import json
import numpy as np

LATSEL = json.load(open("s21/results/latentsel.json"))["rows"]
rows = LATSEL
inc = np.array([r["ship_avg75"] for r in rows])
lat = np.array([r["lat_avg75"] for r in rows])
n   = np.array([r["n"] for r in rows])
fold = np.array([r["fold"] for r in rows])
d = lat - inc
N = len(rows)
rng = np.random.default_rng(220907)

def quartile_stats(strat, d, label):
    qs = np.quantile(strat, [0, .25, .5, .75, 1.0])
    qs[0] -= 1e-9
    bins = np.digitize(strat, qs) - 1
    bins = np.clip(bins, 0, 3)
    print(f"\n  stratifier: {label}")
    for b in range(4):
        sel = bins == b
        x = d[sel]
        se = x.std(ddof=1)/np.sqrt(sel.sum())
        W = int((x < 0).sum()); L = int((x > 0).sum())
        print(f"    Q{b+1} (n={sel.sum():3d}, mean strat={strat[sel].mean():7.3f}): "
              f"d mean {x.mean():+.4f}  SE {se:.4f}  95%CI [{x.mean()-1.96*se:+.4f},{x.mean()+1.96*se:+.4f}]  "
              f"W/L {W}/{L}")
    return bins

print("="*90)
print("L2 REPLICATION on the REAL stratifier used (incumbent's own RMSD = ship_avg75)")
print("="*90)
quartile_stats(inc, d, "ship_avg75 (ORACLE-consuming, the L2 stratifier)")

print("\n" + "="*90)
print("INDEPENDENT CHECK 1: stratify by peptide length n (native-free, weak difficulty proxy,")
print("NOT mechanically coupled to (lat_avg75 - ship_avg75) the way ship_avg75 itself is)")
print("="*90)
quartile_stats(n, d, "n (length)")

print("\n" + "="*90)
print("INDEPENDENT CHECK 2: PLACEBO -- within each n-decile, randomly re-pair lat_avg75 with")
print("a DIFFERENT target's ship_avg75 (breaking the true per-target pairing while preserving")
print("each variable's own marginal AND its association with n).  If regression-to-the-mean")
print("alone can produce L2's monotone quartile trend, it will show up here too, on FAKE pairs.")
print("="*90)
n_deciles = np.digitize(n, np.quantile(n, np.linspace(0,1,6))[1:-1])
placebo_trends = []
B = 500
for b in range(B):
    lat_shuf = lat.copy()
    for dec in np.unique(n_deciles):
        idx = np.where(n_deciles == dec)[0]
        lat_shuf[idx] = lat[rng.permutation(idx)]
    dp = lat_shuf - inc
    qs = np.quantile(inc, [0,.25,.5,.75,1.0]); qs[0]-=1e-9
    bins = np.clip(np.digitize(inc, qs)-1, 0, 3)
    means = [dp[bins==k].mean() for k in range(4)]
    placebo_trends.append(means)
placebo_trends = np.array(placebo_trends)
print("\n  placebo (n=%d resamples), mean +/- sd of quartile means under FAKE pairing:" % B)
for k in range(4):
    print(f"    Q{k+1}: {placebo_trends[:,k].mean():+.4f} +/- {placebo_trends[:,k].std():.4f}   "
          f"(2.5-97.5 pct [{np.percentile(placebo_trends[:,k],2.5):+.4f},"
          f"{np.percentile(placebo_trends[:,k],97.5):+.4f}])")
print("\n  REAL quartile means (from the top block, ORACLE stratifier) for comparison:")
qs = np.quantile(inc, [0,.25,.5,.75,1.0]); qs[0]-=1e-9
bins = np.clip(np.digitize(inc, qs)-1, 0, 3)
for k in range(4):
    print(f"    Q{k+1}: {d[bins==k].mean():+.4f}")

print("\n" + "="*90)
print("DECILE ROBUSTNESS + LEAVE-ONE-TARGET-OUT DRIVER CHECK on the REAL data")
print("="*90)
deciles_inc = np.digitize(inc, np.quantile(inc, np.linspace(0,1,11))[1:-1])
print("  decile means of d = lat_avg75 - ship_avg75, by incumbent-RMSD decile (low=easy):")
for k in range(10):
    sel = deciles_inc == k
    if sel.sum() == 0: continue
    print(f"    D{k+1:2d} (n={sel.sum():2d}, incumbent mean {inc[sel].mean():.2f}): d = {d[sel].mean():+.4f}")

# jackknife the hardest-quartile mean to see if a few targets drive it
qs = np.quantile(inc, [0,.25,.5,.75,1.0]); qs[0]-=1e-9
bins = np.clip(np.digitize(inc, qs)-1, 0, 3)
hardest = d[bins==3]
jk = np.array([np.delete(hardest, i).mean() for i in range(len(hardest))])
print(f"\n  hardest-quartile jackknife: full mean {hardest.mean():+.4f}, "
      f"leave-one-out range [{jk.min():+.4f},{jk.max():+.4f}] over {len(hardest)} targets")

print("\n" + "="*90)
print("THE CLEAN VERSION: is the trend MORE than a straight-line regression of lat on inc predicts?")
print("For bivariate data, E[lat | inc=x] is LINEAR in x under no special interaction; d = lat-inc")
print("is therefore mechanically decreasing in x with slope (beta-1) where beta is the OLS slope of")
print("lat on inc.  A genuine 'crosses zero and reverses which SOURCE wins' finding needs curvature")
print("or a slope steeper than what a single global beta implies -- not just 'd decreases in x'.")
print("="*90)
r = np.corrcoef(inc, lat)[0,1]
beta, alpha0 = np.polyfit(inc, lat, 1)
resid = lat - (alpha0 + beta*inc)
print(f"\n  corr(ship_avg75, lat_avg75) = {r:+.4f}   OLS slope beta = {beta:.4f}  intercept {alpha0:.4f}")
print(f"  mechanical prediction: d = lat - inc = {alpha0:.4f} + ({beta:.4f}-1)*inc + resid")
print(f"  so d's slope on inc is forced to be (beta-1) = {beta-1:.4f} by the LINEAR fit ALONE,")
print(f"  with NO 'routing' interpretation required -- it is what 'lat correlates less than 1:1")
print(f"  with inc' mechanically implies for d=lat-inc.")
qs = np.quantile(inc, [0,.25,.5,.75,1.0]); qs[0]-=1e-9
bins = np.clip(np.digitize(inc, qs)-1, 0, 3)
print("\n  quartile means: REAL d vs the LINEAR-FIT-IMPLIED d (using each target's own inc, no")
print("  per-target noise) vs the RESIDUAL (real d minus the linear prediction) --")
print("  residual is what's left AFTER removing the mechanical linear component:")
for k in range(4):
    sel = bins==k
    d_real = d[sel].mean()
    d_lin  = (alpha0 + beta*inc[sel] - inc[sel]).mean()
    d_res  = resid[sel].mean()
    print(f"    Q{k+1}: real d {d_real:+.4f}   linear-mechanical component {d_lin:+.4f}   "
          f"residual (genuine nonlinearity) {d_res:+.4f}")
print(f"\n  R^2 of linear fit = {r**2:.4f}  (i.e. {100*r**2:.1f}% of lat's variance is 'explained' by inc)")

print("\n" + "="*90)
print("Bootstrap CI on the RESIDUAL quartile means (fold-clustered, target-level)")
print("="*90)
B = 20000
rng2 = np.random.default_rng(9)
for k in range(4):
    sel = np.where(bins==k)[0]
    x = resid[sel]
    boot = x[rng2.integers(0, len(x), size=(B, len(x)))].mean(1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    print(f"    Q{k+1} residual: mean {x.mean():+.4f}  95% CI [{lo:+.4f},{hi:+.4f}]  n={len(x)}")
