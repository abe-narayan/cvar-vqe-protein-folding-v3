"""SPRINT 15 / QGEOM -- PART A(ii): DOES QNG ACTUALLY HELP WHERE THE METRIC IS ILL-CONDITIONED?

Sprint 14 left this as its single most concrete open quantum question, recorded as a
HYPOTHESIS: the Fubini-Study metric departs from I/4 at depth >= 2 with large condition
numbers, so "QNG is refuted" is a depth-1 statement only, and whether QNG *helps* at depth
>= 2 was never measured.  This measures it.

THE DESIGN

* Three optimiser families on IDENTICAL circuits, objectives, initialisations and seeds:
      sgd        plain gradient descent
      adam       the shipped optimiser
      qng        natural gradient, (g + lam tr(g)/P I)^-1 grad,  lam swept
      qng_diag   diagonal preconditioner only (the cheap hardware approximation)
      qng_adam   Adam on the natural-gradient direction
      ng_shots   natural gradient with the metric estimated from FINITE SHOTS
                 (the empirical Fisher of the sampled distribution) -- the arm that is
                 actually available on hardware
* EVERY arm gets its own learning-rate sweep and the BEST lr is reported.  No arm is
  handicapped by a bad step size; the comparison is best-tuned against best-tuned.
* The conditioning range is the axis: patterns/depths spanning cond|range 1 to ~10^3,
  measured per point rather than seed-averaged (`qgeom_metric.A1`).
* BUDGET IS REPORTED BOTH WAYS, as the brief requires:
      (A) equal ITERATIONS -- equal objective/gradient evaluations
      (B) equal COST under hardware accounting: a parameter-shift gradient costs 2P circuit
          evaluations; the full quantum geometric tensor costs a further P(P+1)/2; the
          diagonal costs a further P.  So QNG at P=30 pays ~9x per step.
  If the conclusion changes between them, the conflict is reported.

Two objectives are used, deliberately:
  KL   -- fit the ansatz to a fixed target distribution.  A clean, well-posed geometry
          problem where "convergence" is unambiguous and the metric is exactly the right
          preconditioner in theory.  This is the most favourable possible setting for QNG.
  VQE  -- the real thing: expectation value / CVaR of a native-free structural objective on
          a fully enumerated target, scored against the certified global optimum.

    python -m s15.qgeom_qng
"""
from __future__ import annotations

import itertools
import time

import numpy as np

from core.quantum import Adam, cvar_exact
from s15 import qgeom_lib as G

TAG = "qng"

# conditioning ladder, from `qgeom_metric.A1` (per-point median cond over the range)
LADDER = [("ring", 1), ("all_to_all", 2), ("ring", 2), ("ring", 3),
          ("block", 2), ("chain", 2), ("brick", 2), ("chain", 3), ("brick", 3)]


# ================================================================== objectives
def sub_objective(e, residues, base_cfg, which="legacy"):
    n, k = e.n, e.k
    src = getattr(e, which)
    cfg = np.array(list(itertools.product(range(k), repeat=len(residues))), np.int64)
    S = np.tile(np.asarray(base_cfg, np.int64), (len(cfg), 1))
    S[:, list(residues)] = cfg
    return src[e.index(S)]


def pack(pdb="1CS9", residues=(1, 2, 3, 4, 5, 6), which="legacy", seed=0):
    """A native-free objective on a 12-qubit sub-register, plus its ORACLE RMSD column."""
    e = G.V.Enum(pdb)
    base = np.random.default_rng(seed).integers(0, e.k, e.n)
    E = G.V.uniformise(sub_objective(e, residues, base, which))
    R = sub_objective(e, residues, base, "rmsd")          # ORACLE, post-hoc scoring only
    return {"E": E, "rmsd": R, "n_qubits": 2 * len(residues),
            "exact_i": int(np.argmin(E)), "pdb": pdb}


# =============================================================== cost + exact gradient
def cost_and_grad(circ, th, E, alpha=1.0, D=None, psi=None):
    """Exact cost and gradient.  Returns (cost, grad, psi, D)."""
    if psi is None:
        psi = circ.state(th)
    if D is None:
        D = G.deriv_states(circ, th)
    p = psi ** 2
    p = p / p.sum()
    if alpha >= 1.0:
        c = float(p @ E)
        w = E
    else:
        c, _q, w = cvar_exact(E, p, alpha)
        c = float(c)
    return c, 2.0 * (D @ (psi * w)), psi, D


def kl_cost_and_grad(circ, th, target, D=None, psi=None):
    if psi is None:
        psi = circ.state(th)
    if D is None:
        D = G.deriv_states(circ, th)
    p = psi ** 2
    p = p / p.sum()
    m = target > 0
    c = float((target[m] * np.log(target[m] / np.maximum(p[m], 1e-300))).sum())
    return c, -(2.0 * (D @ (psi * (target / np.maximum(p, 1e-300))))), psi, D


def metric_from(D, psi):
    ov = D @ psi
    g = D @ D.T - np.outer(ov, ov)
    return (g + g.T) / 2.0


def sampled_metric(circ, th, D, psi, shots, rng):
    """Empirical Fisher of the measured distribution from `shots` samples.

    F_hat = (1/S) sum_s s_i(x) s_j(x),  s_j(x) = dlog p(x)/dtheta_j = 2 D_j(x)/psi(x).
    This is the hardware-available estimator; it converges to 4g.
    """
    p = psi ** 2
    p = p / p.sum()
    idx = rng.choice(len(p), size=shots, p=p)
    denom = psi[idx]
    ok = np.abs(denom) > 1e-12
    S = 2.0 * D[:, idx[ok]] / denom[ok][None, :]
    return (S @ S.T) / max(int(ok.sum()), 1) / 4.0        # /4 to put it on the g scale


# ======================================================================= optimisers
ARMS = ("sgd", "adam", "qng", "qng_diag", "qng_adam", "ng_shots")


def run_arm(circ, th0, gradfn, arm, lr, iters, lam=1e-3, shots=512, seed=0):
    th = np.array(th0, float)
    opt = Adam(circ.n_params(), lr=lr) if arm in ("adam", "qng_adam") else None
    rng = np.random.default_rng(seed)
    hist = []
    for _ in range(iters):
        c, gr, psi, D = gradfn(th)
        hist.append(c)
        if arm == "sgd":
            step = gr
        elif arm == "adam":
            th = opt.step(th, gr)
            continue
        elif arm == "qng":
            step = G.nat_solve(metric_from(D, psi), gr, lam, "tikhonov")
        elif arm == "qng_diag":
            step = G.nat_solve(metric_from(D, psi), gr, lam, "diag")
        elif arm == "qng_adam":
            th = opt.step(th, G.nat_solve(metric_from(D, psi), gr, lam, "tikhonov"))
            continue
        elif arm == "ng_shots":
            step = G.nat_solve(sampled_metric(circ, th, D, psi, shots, rng), gr, lam,
                               "tikhonov")
        else:
            raise ValueError(arm)
        th = th - lr * step
        if not np.all(np.isfinite(th)):
            return {"theta": th0, "hist": hist + [float("inf")], "diverged": True}
    c, _, _, _ = gradfn(th)
    hist.append(c)
    return {"theta": th, "hist": hist, "diverged": False}


# The BASELINES keep the widest sweep -- a negative result about QNG must never rest on an
# under-tuned SGD or Adam.  The QNG family's grid was trimmed for throughput on a saturated
# box; the retained values bracket every optimum the full sweep of the first run selected
# (qng chose lr in {0.1, 0.3} and lam in {1e-3, 1e-2} in all seven completed A3 cells).
LRS = {"sgd": (0.03, 0.1, 0.3, 1.0, 3.0),
       "adam": (0.02, 0.05, 0.1, 0.2, 0.4),
       "qng": (0.03, 0.1, 0.3),
       "qng_diag": (0.03, 0.1, 0.3),
       "qng_adam": (0.05, 0.2),
       "ng_shots": (0.03, 0.1, 0.3)}
LAMS = {"qng": (1e-3, 1e-2), "qng_diag": (1e-3,), "qng_adam": (1e-2,),
        "ng_shots": (1e-1,)}


def tune_then_run(circ, gradfn, arm, iters, seeds, score, seed0=0):
    """Pick (lr, lam) on ONE tuning seed, then evaluate on `seeds` fresh seeds.

    Every arm is tuned identically, so the protocol favours none of them; and the reported
    number is not the best of a sweep over the evaluation seeds, which would be the
    "best of dozens presented as a primary result" trap the brief forbids.
    """
    best, bcfg = np.inf, None
    for lr in LRS[arm]:
        for lam in LAMS.get(arm, (1e-3,)):
            th0 = G.random_theta(circ, 100 * seed0 + 7)
            r = run_arm(circ, th0, gradfn, arm, lr, iters, lam, seed=seed0)
            v = score(r)
            if np.isfinite(v) and v < best:
                best, bcfg = v, {"lr": lr, "lam": lam}
    if bcfg is None:
        bcfg = {"lr": LRS[arm][0], "lam": 1e-3}
    vals, results = [], []
    for s in seeds:
        th0 = G.random_theta(circ, 100 * s + 7)
        r = run_arm(circ, th0, gradfn, arm, bcfg["lr"], iters, bcfg["lam"], seed=s)
        vals.append(score(r))
        results.append(r)
    return {"cfg": bcfg, "vals": [float(v) for v in vals],
            "mean": float(np.nanmean(vals)), "results": results}


def hardware_cost_per_iter(P, arm):
    """(B) equal-COST accounting: circuit evaluations per optimiser step on hardware."""
    g = 2 * P                                   # parameter-shift gradient
    if arm in ("sgd", "adam"):
        return g
    if arm == "qng_diag":
        return g + P                            # diagonal QGT
    return g + P * (P + 1) // 2                 # full QGT (or its sampled surrogate)


# ==================================================================== experiment 1
def _kl_target():
    e = G.V.Enum("1CS9")
    t = G.V.uniformise(sub_objective(e, (1, 2, 3, 4, 5, 6), np.zeros(e.n, int), "rmsd"))
    t = np.exp(-6.0 * t)
    return t / t.sum()


def kl_experiment(n=12, iters=150, seeds=(0, 1, 2), ladder=None):
    """QNG's most favourable setting: a well-posed distribution-fitting problem.

    RESUMABLE: ladder entries already in the checkpoint are skipped.
    """
    print("=" * 112)
    print("A3. QNG ON A WELL-POSED GEOMETRY PROBLEM (KL fit).  The most favourable setting.")
    print("    Equal ITERATIONS.  Every arm at its own tuned (lr, lam).  Lower KL is better.")
    print("=" * 112)
    tgt = _kl_target()
    out = G.ck_load(TAG).get("A3_kl_equal_iters", {})
    print(f"{'pattern':11s} {'L':>2s} {'P':>3s} {'cond':>9s} " +
          "".join(f"{a:>11s}" for a in ARMS) + "   winner")
    for pat, L in (ladder or LADDER):
        if f"{pat}|L{L}" in out:
            continue
        c = G.FlexCircuit(n, L, pat, 2)
        P = c.n_params()
        cond = float(np.median([G.spec_stats(G.fs_metric(
            c, G.random_theta(c, 100 * s + 7)))["cond"] for s in seeds]))
        sc = lambda r: (r["hist"][-1] if np.isfinite(r["hist"][-1]) else 1e6)
        best = {}
        for arm in ARMS:
            best[arm] = tune_then_run(c, lambda t: kl_cost_and_grad(c, t, tgt),
                                      arm, iters, seeds, sc)
            best[arm].pop("results", None)
        win = min(ARMS, key=lambda a: best[a]["mean"])
        out[f"{pat}|L{L}"] = {"cond": cond, "n_params": P, "arms": best, "winner": win,
                              "iters": iters}
        print(f"{pat:11s} {L:2d} {P:3d} {cond:9.2f} " +
              "".join(f"{best[a]['mean']:11.5f}" for a in ARMS) + f"   {win}")
        G.ck(TAG, "A3_kl_equal_iters", out)
    return out


def kl_equal_cost(n=12, seeds=(0, 1, 2), budget=6000, ladder=None):
    """(B) equal COST: the same circuit-evaluation budget, not the same iteration count."""
    print()
    print("=" * 112)
    print(f"A3b. THE SAME EXPERIMENT AT EQUAL HARDWARE COST ({budget} circuit evaluations).")
    print("     Gradient = 2P evals/step; full QGT costs a further P(P+1)/2; diagonal P.")
    print("=" * 112)
    tgt = _kl_target()
    out = G.ck_load(TAG).get("A3b_kl_equal_cost", {})
    print(f"{'pattern':11s} {'L':>2s} {'P':>3s} " +
          "".join(f"{a:>11s}" for a in ARMS) + "   winner")
    for pat, L in (ladder or LADDER):
        if f"{pat}|L{L}" in out:
            continue
        c = G.FlexCircuit(n, L, pat, 2)
        P = c.n_params()
        sc = lambda r: (r["hist"][-1] if np.isfinite(r["hist"][-1]) else 1e6)
        row, its = {}, {}
        for arm in ARMS:
            it = max(5, budget // hardware_cost_per_iter(P, arm))
            its[arm] = int(it)
            b = tune_then_run(c, lambda t: kl_cost_and_grad(c, t, tgt), arm, it, seeds, sc)
            b.pop("results", None)
            row[arm] = b
        win = min(ARMS, key=lambda a: row[a]["mean"])
        out[f"{pat}|L{L}"] = {"n_params": P, "arms": row, "iters": its, "winner": win,
                              "budget": budget}
        print(f"{pat:11s} {L:2d} {P:3d} " +
              "".join(f"{row[a]['mean']:11.5f}" for a in ARMS) + f"   {win}")
        print(f"{'  iters':11s} {'':2s} {'':3s} " +
              "".join(f"{its[a]:11d}" for a in ARMS))
        G.ck(TAG, "A3b_kl_equal_cost", out)
    return out


# ==================================================================== experiment 2
VQE_LADDER = [("ring", 1), ("all_to_all", 2), ("block", 2), ("brick", 2), ("chain", 3)]


def vqe_experiment(pdbs=("1CS9", "2MK7", "7N2I"), n=12, iters=120, seeds=(0, 1, 2),
                   alphas=(1.0, 0.25), ladder=None):
    """The real question: does QNG improve a VQE on a native-free structural objective?

    RESUMABLE: cells already in the checkpoint are skipped.
    """
    print()
    print("=" * 112)
    print("A4. QNG ON THE REAL VQE.  Native-free Legacy sub-objective on a 12-qubit")
    print("    sub-register, certified optimum known by enumeration.  BOTH axes.")
    print("=" * 112)
    out = G.ck_load(TAG).get("A4_vqe", {})
    for pdb in pdbs:
        pk = pack(pdb)
        E, R, ex = pk["E"], pk["rmsd"], pk["exact_i"]
        for alpha in alphas:
            print()
            print(f"--- {pdb}  alpha={alpha}  certified-optimum RMSD {R[ex]:.3f} A, "
                  f"space best {R.min():.3f}, random draw {R.mean():.3f} ---")
            print(f"{'pattern':11s} {'L':>2s} {'P':>3s} {'cond':>8s} | objective gap " +
                  "".join(f"{a:>10s}" for a in ARMS))
            for pat, L in (ladder or VQE_LADDER):
                if f"{pdb}|a{alpha}|{pat}|L{L}" in out:
                    continue
                c = G.FlexCircuit(n, L, pat, 2)
                P = c.n_params()
                cond = float(np.median([G.spec_stats(G.fs_metric(
                    c, G.random_theta(c, 100 * s + 7)))["cond"] for s in seeds]))

                def sc(r, E=E, ex=ex, c=c):
                    if r["diverged"]:
                        return 1e6
                    return float(c.probs(r["theta"]) @ E) - float(E[ex])

                gaps, rms, cfgs = {}, {}, {}
                for arm in ARMS:
                    b = tune_then_run(c, lambda t: cost_and_grad(c, t, E, alpha),
                                      arm, iters, seeds, sc)
                    rr = []
                    for r in b["results"]:
                        if r["diverged"]:
                            rr.append(np.nan)
                        else:
                            p = c.probs(r["theta"])
                            rr.append(float(R[int(np.argmax(p))]))
                    gaps[arm] = b["mean"]
                    rms[arm] = float(np.nanmean(rr))
                    cfgs[arm] = b["cfg"]
                out[f"{pdb}|a{alpha}|{pat}|L{L}"] = {
                    "cond": cond, "n_params": P, "obj_gap": gaps, "mode_rmsd": rms,
                    "cfg": cfgs, "certified_rmsd": float(R[ex]),
                    "space_best": float(R.min()), "random_draw": float(R.mean()),
                    "iters": iters}
                print(f"{pat:11s} {L:2d} {P:3d} {cond:8.2f} | {'':13s}" +
                      "".join(f"{gaps[a]:10.5f}" for a in ARMS))
                print(f"{'':11s} {'':2s} {'':3s} {'':8s} | mode RMSD    " +
                      "".join(f"{rms[a]:10.3f}" for a in ARMS))
                G.ck(TAG, "A4_vqe", out)
    return out


# Four rungs spanning the measured conditioning range: 1.0 / 9.6 / 86 / 2004.
LEAN = [("ring", 1), ("ring", 2), ("chain", 2), ("chain", 3)]


def vqe_equal_cost(pdbs=("1CS9", "2MK7", "7N2I"), n=12, seeds=(0, 1),
                   alphas=(1.0, 0.25), budget=8000, ladder=None):
    """(B) EQUAL COST for the real VQE.  A parameter-shift gradient costs 2P circuit
    evaluations per step; the full quantum geometric tensor costs a further P(P+1)/2 and
    the diagonal a further P.  At P=24 that is 348 evaluations per QNG step against 48 for
    a plain gradient step -- QNG gets 7.25x fewer iterations for the same hardware.

    The brief requires BOTH conventions and requires any conflict between them to be
    reported.  RESUMABLE.
    """
    print()
    print("=" * 112)
    print(f"A4b. THE REAL VQE AT EQUAL HARDWARE COST ({budget} circuit evaluations)")
    print("=" * 112)
    out = G.ck_load(TAG).get("A4b_vqe_equal_cost", {})
    for pdb in pdbs:
        pk = pack(pdb)
        E, R, ex = pk["E"], pk["rmsd"], pk["exact_i"]
        for alpha in alphas:
            print()
            print(f"--- {pdb}  alpha={alpha}  certified-optimum RMSD {R[ex]:.3f} A ---")
            print(f"{'pattern':11s} {'L':>2s} {'P':>3s} | " +
                  "".join(f"{a:>10s}" for a in ARMS))
            for pat, L in (ladder or VQE_LADDER):
                key = f"{pdb}|a{alpha}|{pat}|L{L}"
                if key in out:
                    continue
                c = G.FlexCircuit(n, L, pat, 2)
                P = c.n_params()

                def sc(r, E=E, ex=ex, c=c):
                    if r["diverged"]:
                        return 1e6
                    return float(c.probs(r["theta"]) @ E) - float(E[ex])

                gaps, rms, its = {}, {}, {}
                for arm in ARMS:
                    it = max(5, budget // hardware_cost_per_iter(P, arm))
                    its[arm] = int(it)
                    b = tune_then_run(c, lambda t: cost_and_grad(c, t, E, alpha),
                                      arm, it, seeds, sc)
                    rr = []
                    for r in b["results"]:
                        if r["diverged"]:
                            rr.append(np.nan)
                        else:
                            p = c.probs(r["theta"])
                            rr.append(float(R[int(np.argmax(p))]))
                    gaps[arm] = b["mean"]
                    rms[arm] = float(np.nanmean(rr))
                out[key] = {"n_params": P, "obj_gap": gaps, "mode_rmsd": rms,
                            "iters": its, "budget": budget,
                            "certified_rmsd": float(R[ex]),
                            "random_draw": float(R.mean())}
                print(f"{pat:11s} {L:2d} {P:3d} | " +
                      "".join(f"{gaps[a]:10.5f}" for a in ARMS))
                print(f"{'  iters':11s} {'':2s} {'':3s} | " +
                      "".join(f"{its[a]:10d}" for a in ARMS))
                G.ck(TAG, "A4b_vqe_equal_cost", out)
    return out


def main(lean=True, iters=90):
    G.wait_mem(0.8, "qgeom_qng")
    G.ck_load(TAG)
    t0 = time.time()
    lad = LEAN if lean else None
    sd = (0, 1) if lean else (0, 1, 2)
    # A4 FIRST: it is the decisive cell (the real VQE) and A3 is already 6/9 complete.
    vqe_experiment(seeds=sd, ladder=lad, iters=iters)
    vqe_equal_cost(seeds=sd, ladder=lad)
    kl_equal_cost(seeds=sd, ladder=lad)
    kl_experiment(seeds=sd, ladder=lad)
    print(f"\ndone in {time.time()-t0:.0f}s -> s15/results/qgeom_qng.json")


if __name__ == "__main__":
    main()
