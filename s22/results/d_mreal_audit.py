"""s22/results/d_mreal_audit.py -- independent reproduction + attack on s22/mreal.py's
'65% of the m-ladder headroom transfers' finding (the campaign's new load-bearing positive)."""
import json
import numpy as np

d = json.load(open("s22/results/mreal.json"))
rows = d["rows"]
print(f"n={len(rows)}  reps={d['reps']}  ms={d['ms']}  fixed={d['fixed']}  complete={d['complete']}")

sel = np.array([r["sel_heldout"] for r in rows])
fix = np.array([r["fixed_heldout"] for r in rows])
orc = np.array([r["oracle_insample"] for r in rows])
ntied = np.array([r["n_tied"] for r in rows])
folds = None
try:
    folds = np.array([r["fold"] for r in rows])
except KeyError:
    pass

rng = np.random.default_rng(9)

def stat_iid(x, B=20000):
    x = np.asarray(x, float); k = len(x)
    se = x.std(ddof=1) / np.sqrt(k)
    b = x[rng.integers(0, k, size=(B, k))].mean(1)
    return x.mean(), se, 2.8016*se, np.percentile(b, 2.5), np.percentile(b, 97.5)

def stat_fold(x, folds, B=20000):
    fs = np.unique(folds)
    groups = [x[folds == f] for f in fs]
    out = np.empty(B)
    for b in range(B):
        pick = [groups[i] for i in rng.integers(0, len(groups), len(groups))]
        out[b] = np.concatenate(pick).mean()
    return np.percentile(out, 2.5), np.percentile(out, 97.5)

d1 = sel - fix
m, se, mde, lo, hi = stat_iid(d1)
print(f"\n[REPRODUCE] selected - fixed, held out: {m:+.4f} SE {se:.4f} MDE {mde:.4f} "
      f"iid CI [{lo:+.4f},{hi:+.4f}]  W/L {int((d1<0).sum())}/{int((d1>0).sum())}")
if folds is not None:
    flo, fhi = stat_fold(d1, folds)
    print(f"            fold-clustered 95% CI: [{flo:+.4f},{fhi:+.4f}]  "
          f"(mreal.py's own report() does NOT compute this -- BRIEF's own standing practice, e.g. "
          f"d_cvarop.py, reports both beside every primary)")

d2 = orc - fix
m2, se2, mde2, lo2, hi2 = stat_iid(d2)
print(f"\n[REPRODUCE] in-sample oracle - fixed (apparent headroom): {m2:+.4f} [{lo2:+.4f},{hi2:+.4f}]")
frac = 100.0*min(max(m/m2, 0.0), 1.0)
print(f"[REPRODUCE] fraction transferring: {frac:.1f}%")
print(f"[REPRODUCE] median tied rungs: {np.median(ntied):.2f} of 6")

# ---- attack 1: driven by a handful of targets? jackknife the primary mean ----
jk = np.array([np.delete(d1, i).mean() for i in range(len(d1))])
print(f"\n[ATTACK] jackknife leave-one-out range of the primary mean: "
      f"[{jk.min():+.4f},{jk.max():+.4f}]  (full mean {d1.mean():+.4f})")

# ---- attack 2: correlate the per-target effect with n_tied and with target difficulty ----
print(f"\n[ATTACK] corr(d1, n_tied) = {np.corrcoef(d1, ntied)[0,1]:+.3f}  "
      f"(a target with MORE tied rungs has LESS distinct m-information; if the effect concentrates"
      f" on low-n_tied targets that supports 'real distinction', if on high-n_tied it is suspicious)")
lo_tied = ntied <= np.median(ntied)
print(f"    d1 mean on n_tied<=median ({lo_tied.sum()} targets): {d1[lo_tied].mean():+.4f}")
print(f"    d1 mean on n_tied> median ({(~lo_tied).sum()} targets): {d1[~lo_tied].mean():+.4f}")

# ---- attack 3: does the fixed-m=75 half-pool baseline match what a FAIR half-pool comparison needs?
# i.e. is m=75 even the best FIXED choice on half-pools, or does half-pool truncation favour a
# different m as the fair baseline?
print(f"\n[ATTACK] oracle_insample mean {orc.mean():.4f} vs fixed(m=75) mean {fix.mean():.4f} vs "
      f"full-pool incumbent 3.048 (from s21 poolgap.json) -- half-pool fixed-m75 is "
      f"{'WORSE' if fix.mean() > 3.048 else 'BETTER'} than the full pool's own m=75, "
      f"by {fix.mean()-3.048:+.4f} (expected: smaller pool -> worse average, sanity check)")

# ---- attack 4: per-repeat variance -- was aggregation-then-CI (n=126) done correctly, i.e. is the
# reported SE about 126 independent targets, not (126 x 8 reps) pseudo-replicated?  Confirm from
# the source: rows already average over 8 reps x 2 directions BEFORE this file's stat() runs, so
# no pseudoreplication at the population level -- checked by construction (see mreal.py:122-127).
print("\n[CONFIRMED BY READING SOURCE] rows in mreal.json are already averaged over the 8 reps and "
      "both A/B directions per target (mreal.py lines 122-127) BEFORE the population CI is computed, "
      "so there is no repeat-level pseudoreplication in the reported SE -- this was checked, not assumed.")

print("\n" + "="*90)
print("ATTACK 5 (the one the coordinator flagged): the half-pool m=500 rung is DEGENERATE with the")
print("whole half (250 members), since min(500,250)=250 in source. Does the transfer finding depend")
print("on that redundant rung being IN the selectable set? mreal.json does not store the full 6-value")
print("per-half ladder, only the already-collapsed sel/fixed/oracle summaries, so this cannot be")
print("re-derived WITHOUT rerunning against source data with the m=500 rung excluded. Flagging as")
print("OPEN rather than fabricating a number: the raw per-(target,rep,half) 6-vectors were not saved.")
print("="*90)
