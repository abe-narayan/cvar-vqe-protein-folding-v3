"""Finite-shot bias of the empirical CVaR, at the DEPLOYED cell's parameters.

Deployed cell (s29 ledger): n = 9 qubits (512 basis states), alpha = 0.18,
shots = 2048 in the global driver (core/quantum.py:179).  Tail = ceil(0.18*2048) = 369 shots.

The empirical lower-tail CVaR is the mean of the k = ceil(alpha*S) SMALLEST of S draws
from the circuit's distribution over basis states, each carrying its (diagonal) energy.

CLAIM UNDER TEST: E[empirical CVaR] != true CVaR, the gap depends on the STATE (not just on
its true CVaR), and at equal true CVaR the gap is smaller for a more CONCENTRATED state.
If so, the finite-shot objective is not a noisy copy of the exact objective: it has a
different argmin, and the difference systematically rewards concentration.
"""
import numpy as np

rng = np.random.default_rng(7)
N, ALPHA, SHOTS = 512, 0.18, 2048
K = int(np.ceil(ALPHA * SHOTS))

# a plausible diagonal spectrum: energies of 512 candidate-index states, spread ~ 1 unit
E = np.sort(rng.normal(0.0, 1.0, N))

def true_cvar(p):
    """exact lower-tail CVaR of the discrete distribution (p over the sorted energies)."""
    c = np.cumsum(p)
    j = np.searchsorted(c, ALPHA)
    j = min(j, N - 1)
    w = p.copy()
    w[j] = ALPHA - (c[j - 1] if j > 0 else 0.0)
    w[j + 1:] = 0.0
    return float((w[:j + 1] * E[:j + 1]).sum() / ALPHA)

def emp_cvar(p, reps=4000):
    idx = rng.choice(N, size=(reps, SHOTS), p=p)
    x = np.sort(E[idx], axis=1)[:, :K]
    return x.mean(axis=1)

def softmax_state(beta):
    z = np.exp(-beta * (E - E.min()))
    return z / z.sum()

print(f"n=9 -> {N} states, alpha={ALPHA}, shots={SHOTS}, tail k={K}")
print()
print("A. bias as a function of concentration (Boltzmann-like states, beta = inverse temp)")
print(f"{'beta':>6} {'PR (eff #states)':>17} {'true CVaR':>11} {'E[emp CVaR]':>12} {'BIAS':>9} {'sd':>7}")
rows = []
for beta in (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0):
    p = softmax_state(beta)
    pr = 1.0 / (p ** 2).sum()
    tc = true_cvar(p)
    ec = emp_cvar(p)
    rows.append((beta, pr, tc, ec.mean(), ec.mean() - tc, ec.std()))
    print(f"{beta:6.1f} {pr:17.1f} {tc:11.4f} {ec.mean():12.4f} {ec.mean()-tc:+9.4f} {ec.std():7.4f}")

print()
print("B. THE DECISIVE TEST -- two states MATCHED on true CVaR, differing in concentration.")
print("   If the finite-shot objective were a noisy copy of the exact one, the two would tie.")
# build a 2-parameter family: mass q spread uniformly over the lowest L states, rest Boltzmann-ish
def two_block(L, q):
    p = np.zeros(N)
    p[:L] = q / L
    tailmass = 1.0 - q
    w = np.exp(-0.5 * (E[L:] - E[L]))
    p[L:] = tailmass * w / w.sum()
    return p

# search for pairs with (near-)equal true CVaR at different L
cands = []
for L in (2, 4, 8, 16, 32, 64, 128):
    for q in np.linspace(0.02, 0.95, 260):
        p = two_block(L, q)
        cands.append((L, q, true_cvar(p)))
target = np.median([c[2] for c in cands])
best = {}
for L, q, tc in cands:
    if L not in best or abs(tc - target) < abs(best[L][2] - target):
        best[L] = (L, q, tc)
print(f"   matched at true CVaR ~ {target:.4f}")
print(f"{'L (support)':>12} {'q':>6} {'true CVaR':>11} {'E[emp CVaR]':>12} {'BIAS':>9} {'PR':>8}")
for L in sorted(best):
    Lv, q, tc = best[L]
    if abs(tc - target) > 0.02:
        continue
    p = two_block(Lv, q)
    ec = emp_cvar(p)
    pr = 1.0 / (p ** 2).sum()
    print(f"{Lv:12d} {q:6.3f} {tc:11.4f} {ec.mean():12.4f} {ec.mean()-tc:+9.4f} {pr:8.1f}")

print()
print("C. does the bias change the ARGMIN? sweep beta, compare exact vs finite-shot minimiser")
betas = np.linspace(0.0, 40.0, 81)
tcs, ecs = [], []
for b in betas:
    p = softmax_state(b)
    tcs.append(true_cvar(p))
    ecs.append(emp_cvar(p, reps=3000).mean())
tcs, ecs = np.array(tcs), np.array(ecs)
print(f"   exact-CVaR minimising beta      = {betas[tcs.argmin()]:.1f}  (CVaR {tcs.min():.4f})")
print(f"   finite-shot minimising beta     = {betas[ecs.argmin()]:.1f}  (E[emp] {ecs.min():.4f})")
print(f"   exact CVaR at the finite-shot minimiser = {tcs[ecs.argmin()]:.4f}"
      f"   (penalty {tcs[ecs.argmin()] - tcs.min():+.4f})")
print()
print("D. the flat global-minimiser set (Barkoutsos Prop 5.1) under finite shots")
print("   states with EQUAL exact CVaR but different concentration:")
lo = E[:int(ALPHA * N)]
for L in (int(ALPHA*N), int(ALPHA*N)//2, 8, 4, 2, 1):
    p = np.zeros(N)
    p[:L] = ALPHA / L
    p[L:] = (1 - ALPHA) / (N - L)
    tc = true_cvar(p)
    ec = emp_cvar(p)
    print(f"   support {L:4d}: exact {tc:8.4f}   E[emp] {ec.mean():8.4f}   bias {ec.mean()-tc:+8.4f}   sd {ec.std():.4f}")
