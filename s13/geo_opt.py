"""SPRINT 13 GEO -- EXPERIMENT 7: a geometry-aware optimiser at MATCHED evaluation budget.

ARMS (all device-realistic: energies come from SAMPLED bitstrings, never from the enumerated
table, and every unique bitstring is charged to a shared `budget.BudgetedEnergyModel`)

  ``qng``     quantum natural gradient: ``theta <- theta - eta (g + lambda I)^-1 grad``,
              with ``g`` the Fubini-Study metric.  Regularisation and the metric's
              conditioning are reported, since that is where such results usually die.
  ``adam``    Adam on the same sampled score-function gradient (baseline="const").
  ``gd``      plain gradient descent, same gradient.
  ``spsa``    SPSA on the sampled CVaR value -- gradient-free, two evaluations per step.
  ``random``  **MANDATORY CONTROL.**  Uniform random bitstrings.  Boulebnane et al. (npj QI
              2023) found variational performance on a 20-qubit noiseless instance "can be
              matched by random sampling up to a small overhead".  If the variational arms
              do not beat it at matched budget, that is the headline.

BUDGET ACCOUNTING.  One evaluation per UNIQUE bitstring, cache hits free -- the project's
own contract.  The circuit simulation, the probability Jacobian and the Fubini-Study metric
cost ZERO energy evaluations, which is honest for a simulator and is stated as a limitation:
on hardware the metric would need O(P^2) extra circuit evaluations (no extra *energy*
evaluations, but real quantum cost).  That asymmetry is reported, not hidden.

ACCURACY IS REPORTED SEPARATELY FROM THE OBJECTIVE.  Every arm returns its best-seen
bitstring; the ORACLE CA-RMSD is rebuilt FROM THAT BITSTRING (the record's
`rmsd-reporting-basis-mismatch` defect: never score one structure's RMSD against another's
energy).  An arm that optimises the objective better and returns a worse structure is
reported as exactly that.

    python -m s13.geo_opt
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                        # noqa: E402
from budget import BudgetExhausted                     # noqa: E402
from core.quantum import Adam                          # noqa: E402

ALPHA = 0.25
LAYERS = 3
SHOTS = 256
#: Budgets are kept a small fraction of the register (6 % and 25 % at L=6, 1.6 % and 6 % at
#: L=7, 0.4 % and 1.6 % at L=8).  At L=5 a 512-evaluation budget covers half the space and
#: every arm trivially finds the optimum, which measures nothing.
BUDGETS = (256, 1024)
SEEDS = (0, 1, 2, 3, 4)
CELLS = [(p, L, 4) for p in ("1A13", "1A1P", "2BFI") for L in (6, 7, 8)]


class Runner:
    """Shared sampling + budgeted evaluation + best-seen tracking for every arm."""

    def __init__(self, E, circ, budget, seed):
        self.model = G.TableModel(E, circ.n, eval_budget=budget)
        self.circ = circ
        self.rng = np.random.default_rng(seed)
        self.best_e = float("inf")
        self.best_j = None
        self.n_calls = 0
        # Gradient- and step-size units. Legacy energies live on ~1e1 and AMBER
        # single-point energies on ~1e2-1e18, so an unscaled learning rate would compare
        # units rather than optimisers. The scale is estimated ONCE from the FIRST batch
        # of sampled energies the arm pays for -- deployable information, no peeking at
        # the enumerated table -- and every arm of every model uses the same rule.
        self.scale = None
        # TERMINATION. The budget only charges for UNIQUE bitstrings, so once the
        # variational state concentrates every sampled bitstring is a cache hit and the
        # optimiser spends nothing -- a `while True` on BudgetExhausted then never
        # terminates. That is not a budget violation (the arm genuinely stops buying new
        # information) but it is a non-terminating loop, so `stalled` counts consecutive
        # iterations that bought nothing and the arms stop after `STALL_LIMIT` of them.
        # The stall count is REPORTED: an arm that converges early and stops paying is a
        # real result about that arm, not a defect to hide.
        self.stalled = 0
        self._last_spent = 0

    def note_iteration(self) -> bool:
        """True if this iteration bought at least one new energy evaluation."""
        spent = self.model.n_energy_evaluations
        new = spent > self._last_spent
        self._last_spent = spent
        self.stalled = 0 if new else self.stalled + 1
        return new

    def evaluate(self, js):
        out = np.empty(len(js))
        for a, j in enumerate(js):
            out[a] = self.model.energy_index(int(j))
            if out[a] < self.best_e:
                self.best_e, self.best_j = float(out[a]), int(j)
        self.n_calls += len(js)
        if self.scale is None and len(out) >= 8:
            iqr = float(np.quantile(out, 0.75) - np.quantile(out, 0.25))
            self.scale = iqr if iqr > 0 else 1.0
        return out

    def sample(self, theta, shots):
        p = self.circ.probs(theta)
        return self.rng.choice(len(p), size=shots, p=p), p

    def cvar_and_grad(self, theta, alpha, shots, J=None):
        """Sampled CVaR and its score-function gradient, baseline='const' (never 'tail')."""
        x, p = self.sample(theta, shots)
        e = self.evaluate(x)
        m = max(1, int(round(alpha * shots)))
        q = np.partition(e, m - 1)[m - 1]
        val = float(np.sort(e)[:m].mean())
        w = np.where(e <= q, (e - q) / alpha, 0.0) / (self.scale or 1.0)
        w = w - w.mean()                             # constant in x -> unbiased
        if J is None:
            PR = self.circ.probs_batch(self.circ._shift_grid(theta, np.pi / 2))
            J = (PR[0::2] - PR[1::2]) / 2.0
        Gl = J[:, x] / np.maximum(p[x], 1e-15)[None, :]
        return val, (Gl * w[None, :]).mean(1)

    def cvar_only(self, theta, alpha, shots):
        x, _ = self.sample(theta, shots)
        e = self.evaluate(x)
        m = max(1, int(round(alpha * shots)))
        return float(np.sort(e)[:m].mean()) / (self.scale or 1.0)


#: Consecutive zero-charge iterations after which an arm is declared converged and stopped.
STALL_LIMIT = 40


def arm_random(E, circ, budget, seed):
    r = Runner(E, circ, budget, seed)
    try:
        while r.stalled < STALL_LIMIT:
            r.evaluate(r.rng.integers(0, 1 << circ.n, size=64))
            r.note_iteration()
    except BudgetExhausted:
        pass
    return r, {"stalled_out": r.stalled >= STALL_LIMIT}


def arm_grad(E, circ, budget, seed, kind, alpha=ALPHA, shots=SHOTS, lr=0.10,
             lam=1e-3):
    r = Runner(E, circ, budget, seed)
    th = G.init_theta(circ.n_params(), seed)
    opt = Adam(circ.n_params(), lr=lr) if kind == "adam" else None
    conds, lams, iters = [], [], 0
    try:
        while r.stalled < STALL_LIMIT:
            if kind == "qng":
                g = G.fs_metric(circ, th)
                w = np.linalg.eigvalsh((g + g.T) / 2)
                conds.append(float(w.max() / max(w.min(), 1e-300)))
                lams.append(lam)
                _val, grad = r.cvar_and_grad(th, alpha, shots)
                step = np.linalg.solve(g + lam * np.eye(len(g)), grad)
                th = th - lr * step
            else:
                _val, grad = r.cvar_and_grad(th, alpha, shots)
                th = opt.step(th, grad) if opt is not None else th - lr * grad
            iters += 1
            r.note_iteration()
    except BudgetExhausted:
        pass
    return r, {"iters": iters, "stalled_out": r.stalled >= STALL_LIMIT,
               "metric_cond_mean": float(np.mean(conds)) if conds else None,
               "metric_cond_max": float(np.max(conds)) if conds else None,
               "qng_lambda": lam if kind == "qng" else None}


def arm_spsa(E, circ, budget, seed, alpha=ALPHA, shots=SHOTS, a=0.60, c=0.20):
    r = Runner(E, circ, budget, seed)
    th = G.init_theta(circ.n_params(), seed)
    rng = np.random.default_rng(seed + 991)
    it = 0
    try:
        while r.stalled < STALL_LIMIT:
            it += 1
            ak = a / (it + 10) ** 0.602
            ck = c / it ** 0.101
            d = rng.choice([-1.0, 1.0], size=th.size)
            fp = r.cvar_only(th + ck * d, alpha, shots)
            fm = r.cvar_only(th - ck * d, alpha, shots)
            th = th - ak * (fp - fm) / (2 * ck) * d
            r.note_iteration()
    except BudgetExhausted:
        pass
    return r, {"iters": it, "stalled_out": r.stalled >= STALL_LIMIT}


def run_cell(pdb, L, k, budget):
    n = int(round(np.log2(k))) * L
    circ = G.circuit(n, layers=LAYERS)
    R = G.rmsd_table(pdb, L, k)
    rows = []
    for v in G.VARIANTS:
        E = G.variant_table(pdb, L, k, v)
        emin, emax = float(E.min()), float(E.max())
        emed = float(np.median(E))
        for seed in SEEDS:
            for arm, fn in (("random", arm_random),
                            ("gd", lambda *A: arm_grad(*A, kind="gd")),
                            ("adam", lambda *A: arm_grad(*A, kind="adam")),
                            ("qng", lambda *A: arm_grad(*A, kind="qng")),
                            ("spsa", arm_spsa)):
                t0 = time.time()
                r, extra = fn(E, circ, budget, seed)
                rows.append({
                    "pdb": pdb, "L": L, "k": k, "n_qubits": n, "model": v, "arm": arm,
                    "budget": budget, "seed": seed, "alpha": ALPHA, "shots": SHOTS,
                    "layers": LAYERS,
                    "unique_evals_spent": int(r.model.n_energy_evaluations),
                    "best_energy": r.best_e,
                    "best_energy_gap_frac": float((r.best_e - emin)
                                                  / max(emax - emin, 1e-300)),
                    "best_energy_gap_robust": float(
                        (r.best_e - emin) / max(emed - emin, 1e-300)),
                    "best_energy_percentile": float((E < r.best_e).mean() * 100),
                    "scale_used": r.scale,
                    "ORACLE_rmsd_of_best_seen": float(R[r.best_j])
                    if r.best_j is not None else None,
                    "ORACLE_rmsd_pool_best": float(R.min()),
                    "ORACLE_rmsd_pool_mean": float(R.mean()),
                    "wall_s": time.time() - t0, **extra})
        print("  %-5s L=%-2d k=%d %-11s budget=%d done" % (pdb, L, k, v, budget),
              flush=True)
    return rows


def main():
    rows = []
    for pdb, L, k in CELLS:
        if not os.path.exists(os.path.join(G.CACHE, G._tag(pdb, L, k, "amber") + ".npy")):
            continue
        for b in BUDGETS:
            G.gate(tries=1)
            rows += run_cell(pdb, L, k, b)
            G.write("geo_opt", {"what": "QNG / Adam / GD / SPSA / uniform-random at matched "
                                        "unique-evaluation budget", "rows": rows})
    G.write("geo_opt", {"what": "QNG / Adam / GD / SPSA / uniform-random at matched "
                                "unique-evaluation budget", "rows": rows})


if __name__ == "__main__":
    main()
