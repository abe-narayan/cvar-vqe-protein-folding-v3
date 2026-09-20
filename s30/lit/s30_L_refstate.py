"""Closed-form arithmetic on the DOPE-style finite-sphere reference state.
NOT a measurement. Pure arithmetic on a published closed form, evaluated at the
chain lengths of this instrument (9-16 aa) against a typical small protein.

Reference state (Shen & Sali 2006, DOPE): non-interacting points uniformly
distributed in a ball of radius a = sqrt(5/3) * Rg.
Pair-distance density of two uniform points in a ball of radius a:
    f(d;a) = (3 d^2 / a^3) * (1 - (3/4)(d/a) + (1/16)(d/a)^3),  0 <= d <= 2a
(verified: integrates to 1 over [0,2a]).

The infinite / ideal-gas reference is f_inf(d) ~ d^2.
The ratio  R(d;a) = f(d;a) / (d^2 normalised)  = 1 - (3/4)t + (1/16)t^3,  t = d/a
is the FINITE-SIZE CORRECTION FACTOR.  It is exactly 1 at d=0 and exactly 0 at d=2a.
"""
import numpy as np

def a_of_n(n, c=2.2, p=0.38):
    """a = sqrt(5/3) Rg, Rg from the folded-protein scaling law used by s27 RG_LAW."""
    rg = c * n ** p
    return np.sqrt(5.0 / 3.0) * rg, rg

def corr(d, a):
    t = np.asarray(d, float) / a
    out = 1.0 - 0.75 * t + t ** 3 / 16.0
    return np.where(t <= 2.0, out, 0.0)

print("=== ball radius a = sqrt(5/3)*Rg, Rg = 2.2 n^0.38 (s27 RG_LAW) ===")
print(f"{'n':>5} {'Rg':>7} {'a':>7} {'2a (support)':>13}")
for n in (9, 13, 16, 30, 60, 150):
    a, rg = a_of_n(n)
    print(f"{n:>5} {rg:7.2f} {a:7.2f} {2*a:13.2f}")

print()
print("=== finite-size correction factor R(d;a) = f_ball/f_ideal, by chain length ===")
ds = [4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 15.0]
print(f"{'n':>5} " + " ".join(f"{d:>7.0f}A" for d in ds))
for n in (9, 13, 16, 30, 60, 150):
    a, _ = a_of_n(n)
    row = corr(np.array(ds), a)
    print(f"{n:>5} " + " ".join(f"{v:8.3f}" for v in row))

print()
print("=== the SAME correction expressed as an energy, kT units: -ln R ===")
print("    (this is what a fixed-reference potential ADDS to a candidate's score")
print("     relative to a size-matched one, per pair, at that distance)")
print(f"{'n':>5} " + " ".join(f"{d:>7.0f}A" for d in ds))
for n in (9, 13, 16, 150):
    a, _ = a_of_n(n)
    r = corr(np.array(ds), a)
    row = np.where(r > 0, -np.log(np.maximum(r, 1e-12)), np.inf)
    print(f"{n:>5} " + " ".join((f"{v:8.3f}" if np.isfinite(v) else "     inf") for v in row))

print()
print("=== scale sensitivity: d/dlnRg of the total reference term, n=13 ===")
print("How much does the reference term move if the candidate is 10% more compact?")
n = 13
a0, rg0 = a_of_n(n)
# a uniform 10% contraction of the structure: all d -> 0.9 d, Rg -> 0.9 Rg
# FIXED reference (a stays a0):   sum_pairs -ln R(0.9 d; a0) - (-ln R(d; a0))
# MATCHED reference (a -> 0.9a0): R(0.9d; 0.9a0) == R(d; a0) EXACTLY (scale invariant)
rng = np.random.default_rng(0)
# a plausible Ca-Ca distance set for a compact 13-mer: use the ideal-ball distribution itself
u = rng.random(200000)
# inverse-cdf by rejection on f(d;a0)
cand = rng.random(400000) * 2 * a0
acc = rng.random(400000) * (3 * (2 * a0 / np.sqrt(3)) ** 2 / a0 ** 3)
fv = (3 * cand ** 2 / a0 ** 3) * corr(cand, a0)
d_s = cand[acc < fv][:20000]
for lam in (0.90, 0.95, 1.05, 1.10):
    fixed = -np.log(np.maximum(corr(lam * d_s, a0), 1e-12)).sum() / len(d_s)
    base = -np.log(np.maximum(corr(d_s, a0), 1e-12)).sum() / len(d_s)
    print(f"  scale {lam:4.2f}: fixed-reference shift = {fixed - base:+7.4f} kT per pair"
          f"   | size-matched reference shift = {0.0:+7.4f} (exact, by scale invariance)")
