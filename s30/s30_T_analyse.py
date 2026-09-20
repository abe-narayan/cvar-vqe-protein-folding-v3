"""S30 lane T -- analysis of s30_T_bits.json against the bars in s30/PREREG_S30_T.md."""
import json
import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "s24"))
from s24 import stats_lib as ST  # noqa: E402

d = json.load(open(os.path.join(ROOT, "s30", "results", "s30_T_bits.json")))
rows = d["rows"]
NS = d["NS"]
pdbs = [r["pdb"] for r in rows]
n = len(rows)

# pinned folds
folds = ST.pinned_folds(pdbs) if hasattr(ST, "pinned_folds") else None

L = np.array([r["ladder"] for r in rows])          # (126, 10) best-of-N in the deployed order
R = np.log2(np.array(NS, dtype=float))
mean_lad = L.mean(0)
se_lad = L.std(0, ddof=1) / np.sqrt(n)

print("=" * 78)
print("M1  THE ORACLE ORDER-STATISTIC LADDER INSIDE THE DEPLOYED TOP-500 (n=126)")
print("=" * 78)
print("   N      R=log2N     mean best-of-N (A)     SE      marginal A/bit")
prev = None
for k, N in enumerate(NS):
    marg = ""
    if prev is not None:
        dr = R[k] - R[k - 1]
        marg = f"{(prev - mean_lad[k]) / dr:8.4f}"
    print(f"{N:6d}   {R[k]:7.3f}   {mean_lad[k]:14.4f}   {se_lad[k]:7.4f}   {marg}")
    prev = mean_lad[k]

# ---- P1a: fit D(R) = a + c 2^{-R/gamma} on the first 9 rungs (N=1..256, exact powers of 2)
from scipy.optimize import curve_fit  # noqa: E402


def model(r, a, c, g):
    return a + c * np.power(2.0, -r / g)


Rp = R[:9]
Yp = mean_lad[:9]
p0 = [1.0, 3.0, 3.0]
popt, _ = curve_fit(model, Rp, Yp, p0=p0, maxfev=200000)
a, c, g = popt
pred = model(Rp, *popt)
ss_res = float(((Yp - pred) ** 2).sum())
ss_tot = float(((Yp - Yp.mean()) ** 2).sum())
r2 = 1 - ss_res / ss_tot
print(f"\nP1a FIT  D(R) = a + c*2^(-R/gamma):  a={a:.4f} A   c={c:.4f}   gamma={g:.4f}   R2={r2:.6f}")
print(f"   bar: R2>=0.99 and gamma in [2,6]; falsified if R2<0.97 or gamma outside [1.5,8]")
print(f"   VERDICT P1a: {'HELD' if (r2>=0.99 and 2<=g<=6) else ('FALSIFIED' if (r2<0.97 or g<1.5 or g>8) else 'PARTIAL')}")

marg7 = c * np.log(2) / g * 2 ** (-7.0 / g)
print(f"\nP1b  marginal -dD/dR at R=7 (analytic from the fit): {marg7:.4f} A/bit")
emp7 = (mean_lad[6] - mean_lad[7])  # N=64 -> 128, one bit
print(f"     empirical one-bit step 64->128: {emp7:.4f} A/bit")
print(f"   bar: [0.05,0.14]; falsified outside [0.03,0.20]")
print(f"   VERDICT P1b: {'HELD' if 0.05<=marg7<=0.14 else ('FALSIFIED' if (marg7<0.03 or marg7>0.20) else 'PARTIAL')}")

print(f"\nP1c  fitted floor a = {a:.4f} A ; D(7) - a = {mean_lad[7]-a:.4f} A")
print(f"     true pool floor (best of all 500)          = {mean_lad[9]:.4f} A")
print(f"     true UNIVERSE floor (best of all nw)       = {np.mean([r['rr_univ_min'] for r in rows]):.4f} A")
print(f"   bar: a<=1.6 and D(7)-a>=0.4; falsified if D(7)-a<0.2 or a>2.0")
ok1c = (a <= 1.6) and (mean_lad[7] - a >= 0.4)
bad1c = (mean_lad[7] - a < 0.2) or (a > 2.0)
print(f"   VERDICT P1c: {'HELD' if ok1c else ('FALSIFIED' if bad1c else 'PARTIAL')}")

# ---- M2
print()
print("=" * 78)
print("M2  WHAT THE 3,252 BITS OF RETRIEVAL CHOICE BOUGHT")
print("=" * 78)
blos500 = L[:, 9]
rand500 = np.array([r["rand500"] for r in rows])
print(f"   BLOSUM top-500 ORACLE best : {blos500.mean():.4f} A  (SE {blos500.std(ddof=1)/np.sqrt(n):.4f})")
print(f"   random  500  ORACLE best   : {rand500.mean():.4f} A  (SE {rand500.std(ddof=1)/np.sqrt(n):.4f})")
try:
    cmp = ST.compare(blos500, rand500, folds)
    ST.render(cmp, "M2  BLOSUM-500 best - RANDOM-500 best (NEGATIVE = the key helps)")
except Exception as e:
    dd = blos500 - rand500
    se = dd.std(ddof=1) / np.sqrt(n)
    print(f"   effect {dd.mean():+.4f}  median {np.median(dd):+.4f}  SE {se:.4f}  MDE {2.8016*se:.4f}"
          f"  effect/MDE {dd.mean()/(2.8016*se):+.2f}   W/L {(dd<0).sum()}/{(dd>0).sum()}  [{e}]")

# search-equivalent bits: where on the RANDOM ladder does BLOSUM-500 land?
beq = []
for r in rows:
    rl = {int(k): v for k, v in r["rand_ladder"].items()}
    Ns = np.array(sorted(rl))
    Vs = np.array([rl[k] for k in Ns])
    target = r["ladder"][9]
    # random ladder is decreasing in N; find N* with V(N*) = target by log-linear interpolation
    if target >= Vs[0]:
        nstar = 1.0
    elif target <= Vs[-1]:
        nstar = float(Ns[-1])         # saturated: lower bound only
    else:
        k = int(np.searchsorted(-Vs, -target))
        k = min(max(k, 1), len(Ns) - 1)
        x0, x1 = np.log2(Ns[k - 1]), np.log2(Ns[k])
        y0, y1 = Vs[k - 1], Vs[k]
        nstar = 2 ** (x0 + (y0 - target) * (x1 - x0) / (y0 - y1)) if y0 != y1 else float(Ns[k])
    beq.append(np.log2(nstar / 500.0))
beq = np.array(beq)
sat = np.array([r["ladder"][9] <= min(float(v) for v in r["rand_ladder"].values()) for r in rows])
print(f"\nP2a  search-equivalent bits  b_ret = log2(N*/500):  mean {beq.mean():+.4f}  median {np.median(beq):+.4f}"
      f"  SE {beq.std(ddof=1)/np.sqrt(n):.4f}")
print(f"     saturated (N* beyond the measured random ladder, lower bound only): {int(sat.sum())}/{n}")
print(f"     positive on {int((beq>0).sum())}/{n} targets")
print(f"   bar: b_ret <= 3.0 bits; falsified if > 4.0")
print(f"   VERDICT P2a: {'HELD' if beq.mean()<=3.0 else ('FALSIFIED' if beq.mean()>4.0 else 'PARTIAL')}")
dd = blos500 - rand500
print(f"\nP2b  BLOSUM-500 minus random-500 = {dd.mean():+.4f} A  (registered [+0.10,+0.30] as a GAIN,"
      f" i.e. dd in [-0.30,-0.10]; falsified outside dd in [-0.45, 0.0])")
print(f"   VERDICT P2b: {'HELD' if -0.30<=dd.mean()<=-0.10 else 'FALSIFIED / SIGN REVERSED' if dd.mean()>0 else 'PARTIAL'}")

# ---- M3
print()
print("=" * 78)
print("M3  RANK COLLAPSE -- THE PRE-CHECK BEFORE ANY LIFTED OBJECTIVE")
print("=" * 78)
rsc = np.array([r["rs_coord"] for r in rows])
rsd = np.array([r["rs_dist"] for r in rows])
tc = np.array([r["top_coord"] for r in rows])
td = np.array([r["top_dist"] for r in rows])
dco = np.array([r["dcoord"] for r in rows])
npr = np.array([r["npairs"] for r in rows])
for nm, v, dim in [("A_coord (500 x 3n)", rsc, dco), ("A_dist (500 x npairs)", rsd, npr)]:
    print(f"   r_stable {nm}: mean {v.mean():.3f}  median {np.median(v):.3f}"
          f"  min {v.min():.3f}  max {v.max():.3f}   nominal dim mean {dim.mean():.1f}")
print(f"   lambda1/sum  coord {tc.mean():.3f}   dist {td.mean():.3f}")
print(f"   bar P3a: rs_coord in [1.5,6], rs_dist in [1.2,4]; the clause that matters is rs_dist>10 -> FALSIFIED")
ok3 = (1.5 <= rsc.mean() <= 6) and (1.2 <= rsd.mean() <= 4)
print(f"   VERDICT P3a: {'HELD' if ok3 else ('FALSIFIED' if rsd.mean()>10 else 'PARTIAL')}")
print(f"   P3b decision rule: rs_dist mean {rsd.mean():.3f} -> "
      f"{'CLOSED at the encoding level (rs<2.0)' if rsd.mean()<2.0 else f'a ceil(rs)+1 = {int(np.ceil(rsd.mean()))+1}-parameter family'}")

# reachable-set cap in bits (Sauer-Shelah on halfspaces of VC dim d+1)
from math import lgamma, log2  # noqa: E402


def logC(N, k):
    if k < 0 or k > N:
        return -np.inf
    return (lgamma(N + 1) - lgamma(k + 1) - lgamma(N - k + 1)) / np.log(2)


for dim, nm in [(int(np.ceil(rsd.mean())) + 1, "effective (stable rank of A_dist)"),
                (int(round(npr.mean())), "nominal npairs"),
                (int(round(dco.mean())) - 6, "nominal 3n-6")]:
    cap = log2(sum(2 ** logC(500, i) for i in range(0, min(dim + 1, 500) + 1)))
    print(f"   reachable-tail cap, halfspace class of VC dim {dim+1:3d} ({nm:34s}): {cap:8.1f} bits"
          f"   vs log2 C(500,75) = {logC(500,75):.1f} bits of set choice")
