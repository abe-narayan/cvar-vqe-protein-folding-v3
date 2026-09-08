"""SPRINT 15 / QGEOM -- PART D: CVaR, EXTENDED NOT REPEATED.

Sprint 14 verified the CVaR value estimator to 1.7e-15 and recorded three defects (the
`baseline="tail"` gradient bias with closed form `-c grad P(E<q)`; a sampled-CVaR upward bias
at non-integer `alpha*N`; and `dCVaR/dp == 0` iff `p(argmin E) >= alpha`).  **None of that is
re-verified here.**  It also established that `alpha` is not a learning rate, that small alpha
collapses CVaR to an argmin-finder, and that plain expectation-value VQE (alpha=1) returns the
best structure.

THE OPEN QUESTION THIS MODULE ATTACKS.  All of that priced CVaR as a **single-structure**
optimiser.  The recorded law that the terminal operator consumes the **set MEAN**
(`d_out = 1.16 d_set_mean + 0.04 d_set_best`, R2 0.89) says the readout that actually matters
is a property of a SET, and CVaR at level `alpha` is exactly a device for shaping the mass in
a tail -- i.e. for shaping a set.  So:

    Is there any regime in which tail-shaping is the right move, and what decides it?

Four experiments, each with the MANDATORY control the brief specifies -- **best-of-N from the
same initial distribution at the same budget**, the arm Sprint 14 found running VQE loses to
0/12 by +0.65 to +1.32 A.

    D1  TRAJECTORY.  entropy, effective sample size, tail mass, mode concentration,
        distribution movement, expected structural quality, and top-m SET quality including
        the coordinate average, tracked through the optimisation at every alpha.
    D2  THE CONTROL TEST at the set level.  Matched budget, three readouts:
        argmin-by-objective (Sprint 14's), top-m set MEAN (the recorded terminal operator),
        and the coordinate average of the top-m (the real operator).
    D3  THE REGIME SWEEP.  objective quality x alpha x readout, looking for ANY cell where
        tail-shaping wins.  Pre-registered prediction: none, because concentration is the
        wrong move when discrimination binds and the set readouts reward BREADTH.
    D4  CVaR AS A CONSTRAINT, not as the objective.  The one use of a risk measure this
        project has never tried: minimise the expected structural term subject to the TAIL of
        a clash term.  This is what CVaR was invented for and it is not what Sprint 14 tested.

GENEROSITY TO VQE IS DELIBERATE.  Every VQE arm gets the EXACT gradient from the statevector
(no shot noise, no SPSA), and the evaluation budget is charged only for the samples drawn for
the readout.  The control gets the same number of samples.  If VQE still loses, it loses from
its best case.

NATIVE-FREE.  Every objective is Legacy / retrieval prior / distogram.  `Enum.rmsd` and the
native CA trace are read only for post-hoc scoring.

    python -m s15.qgeom_cvar
"""
from __future__ import annotations

import time

import numpy as np

from core.quantum import Adam, cvar_exact
from s12 import instrument as I
from s15 import qgeom_lib as G

TAG = "cvar"

TARGETS = ("1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5")
ALPHAS = (1.0, 0.25, 0.05, 0.01)


# ============================================================ structural readouts
class Struct:
    """CA traces and set-level readouts for one enumerated target.

    `residues=None` uses the full n=9 register (18 qubits).  A tuple restricts to a
    sub-register with the other residues frozen at a fixed native-free base configuration;
    every structure is still a real CA trace and every RMSD a real Kabsch CA-RMSD, so a
    sub-register experiment is a smaller instance of the same problem, not a proxy.
    """

    def __init__(self, pdb, residues=None):
        import itertools
        self.e = G.V.Enum(pdb)
        u = I.load_univ(pdb)
        self.nat = np.asarray(u["nat_ca"], float)
        self.pdb = pdb
        e = self.e
        self.res = tuple(range(e.n)) if residues is None else tuple(residues)
        self.n_qubits = 2 * len(self.res)
        base = np.random.default_rng(0).integers(0, e.k, e.n)
        cfg = np.array(list(itertools.product(range(e.k), repeat=len(self.res))),
                       np.int64)
        S = np.tile(base, (len(cfg), 1))
        S[:, list(self.res)] = cfg
        self.full = e.index(S)                        # sub index -> full config index
        self.rmsd = e.rmsd[self.full]                 # ORACLE, post-hoc scoring only

    def build(self, idx):
        S = self.e.states(self.full[np.asarray(idx, np.int64)])
        r = np.arange(self.e.n)
        phi = self.e.PHI[r[None, :], S]
        psi = self.e.PSI[r[None, :], S]
        return np.array([I.build_ca(phi[i], psi[i]) for i in range(len(S))])

    def coord_avg_rmsd(self, idx):
        """RMSD of the COORDINATE AVERAGE of a set -- the recorded terminal operator.

        No projection: the brief records restoring valid bonds costs 0.157 A by any route,
        a constant offset that cannot reorder arms.  Reported as `coordavg_raw`.
        """
        W = self.build(idx)
        if len(W) == 1:
            return float(I.kabsch_rmsd_batch(W, self.nat)[0])
        A = I.coordinate_average(W)
        A = A[0] if isinstance(A, tuple) else A
        return float(I.kabsch_rmsd_batch(np.asarray(A, float)[None], self.nat)[0])

    def set_readouts(self, idx, E, ms=(1, 5, 20, 75)):
        """Everything a downstream operator could consume from a drawn SET."""
        idx = np.asarray(idx, np.int64)
        if idx.size == 0:
            return {}
        order = idx[np.argsort(E[idx], kind="mergesort")]
        R = self.rmsd                                     # ORACLE, post-hoc scoring only
        out = {"argmin_rmsd": float(R[order[0]]),
               "set_best_rmsd": float(R[idx].min()),
               "set_mean_rmsd": float(R[idx].mean()),
               "n_drawn": int(idx.size), "n_distinct": int(np.unique(idx).size)}
        for m in ms:
            top = order[:m]
            out[f"top{m}_mean_rmsd"] = float(R[top].mean())
            out[f"top{m}_best_rmsd"] = float(R[top].min())
            out[f"top{m}_coordavg_rmsd"] = self.coord_avg_rmsd(np.unique(top))
        return out


# ==================================================================== the VQE
def make_circ(nq, layers=2, pattern="ring"):
    return G.FlexCircuit(nq, layers, pattern, 2)


def cvar_cost_grad(circ, th, E, alpha, lowmem=True):
    """Exact cost and gradient via ADJOINT differentiation -- O(gates), not O(P*gates).

    Verified against the O(P) derivative-state route to 2.4e-15 (`qgeom_lib.verify_adjoint`).
    """
    psi = circ.state(th)
    p = psi ** 2
    p = p / p.sum()
    if alpha >= 1.0:
        c, w = float(p @ E), E
    else:
        c, _q, w = cvar_exact(E, p, alpha)
        c = float(c)
    return c, G.grad_adjoint(circ, th, 2.0 * psi * w), p


def dist_stats(p, E, alpha, p0=None):
    q = np.quantile(E, alpha) if alpha < 1.0 else E.max()
    return {"entropy_bits": float(-(p[p > 0] * np.log2(p[p > 0])).sum()),
            "ess_frac": float(1.0 / np.sum(p ** 2) / len(p)),
            "max_p": float(p.max()),
            "tail_mass": float(p[E < q].sum()),
            "tv_from_init": float(0.5 * np.abs(p - p0).sum()) if p0 is not None else 0.0}


def run_vqe(circ, E, alpha, iters=200, lr=0.10, seed=0, lowmem=False, track_every=0,
            struct=None, rmsd=None):
    th = G.random_theta(circ, seed)
    opt = Adam(circ.n_params(), lr=lr)
    p0 = circ.probs(th).copy()
    traj = []
    for it in range(iters):
        c, gr, p = cvar_cost_grad(circ, th, E, alpha, lowmem)
        if track_every and (it % track_every == 0 or it == iters - 1):
            row = {"iter": it, "cost": c, **dist_stats(p, E, alpha, p0)}
            if rmsd is not None:
                row["E_p_rmsd"] = float(p @ rmsd)
                row["mode_rmsd"] = float(rmsd[int(np.argmax(p))])
                for m in (5, 20, 75):
                    top = np.argsort(p)[::-1][:m]
                    row[f"topP{m}_mean_rmsd"] = float(rmsd[top].mean())
                    if struct is not None and it % (track_every * 4) == 0:
                        row[f"topP{m}_coordavg_rmsd"] = struct.coord_avg_rmsd(top)
            traj.append(row)
        th = opt.step(th, gr)
    p = circ.probs(th)
    return {"theta": th, "p_final": p, "p_init": p0, "traj": traj,
            "final": dist_stats(p, E, alpha, p0)}


# =========================================================================== D1
def trajectory(pdbs=("1CS9", "7N2I"), alphas=ALPHAS, iters=200, seed=0, w=0.25,
               residues=None):
    print("=" * 112)
    print("D1. THE CVaR TRAJECTORY -- what tail-shaping does to the DISTRIBUTION")
    print("=" * 112)
    from s14.vqe_hamil import tabulate, combine
    out = {}
    for pdb in pdbs:
        st = Struct(pdb, residues)
        pr, ds = tabulate(pdb)
        E = G.V.uniformise(combine(pr, ds, w)[st.full])
        circ = make_circ(st.n_qubits)
        for a in alphas:
            r = run_vqe(circ, E, a, iters, seed=seed, lowmem=True, track_every=20,
                        struct=st, rmsd=st.rmsd)
            out[f"{pdb}|a{a}"] = {"traj": r["traj"], "final": r["final"]}
            print(f"\n--- {pdb} alpha={a} ---")
            print(f"{'iter':>5s} {'cost':>9s} {'H bits':>7s} {'ESSfrac':>8s} "
                  f"{'max p':>7s} {'tailmass':>9s} {'TV':>6s} {'E_p[RMSD]':>10s} "
                  f"{'mode':>7s} {'topP20mean':>11s} {'topP20avg':>10s}")
            for row in r["traj"]:
                print(f"{row['iter']:5d} {row['cost']:9.5f} {row['entropy_bits']:7.2f} "
                      f"{row['ess_frac']:8.5f} {row['max_p']:7.4f} "
                      f"{row['tail_mass']:9.4f} {row['tv_from_init']:6.3f} "
                      f"{row.get('E_p_rmsd', np.nan):10.3f} "
                      f"{row.get('mode_rmsd', np.nan):7.3f} "
                      f"{row.get('topP20_mean_rmsd', np.nan):11.3f} "
                      f"{row.get('topP20_coordavg_rmsd', np.nan):10.3f}")
            G.ck(TAG, "D1_trajectory", out)
    return out


# =========================================================================== D2
def control_test(pdbs=TARGETS, alphas=ALPHAS, iters=200, seeds=(0, 1, 2), budget=4096,
                 w=0.25, residues=None, tag="D2_control"):
    """THE MANDATORY CONTROL: best-of-N from the SAME initial distribution, same budget."""
    print()
    print("=" * 112)
    print(f"D2. MATCHED-BUDGET CONTROL TEST ({budget} samples for every arm).")
    print("    Control = best-of-N from the arm's OWN untrained initial distribution.")
    print("    Readouts: argmin-by-objective, top-m set MEAN, top-m COORDINATE AVERAGE.")
    print("=" * 112)
    from s14.vqe_hamil import tabulate, combine
    prev = G.ck_load(TAG)
    out = prev.get(tag, {})
    for pdb in pdbs:
        st = Struct(pdb, residues)
        pr, ds = tabulate(pdb)
        E = G.V.uniformise(combine(pr, ds, w)[st.full])
        circ = make_circ(st.n_qubits)
        for a in alphas:
            for s in seeds:
                key = f"{pdb}|a{a}|s{s}"
                if key in out:
                    continue
                t0 = time.time()
                rng = np.random.default_rng(1000 * s + 7)
                r = run_vqe(circ, E, a, iters, seed=s, lowmem=True)
                ctrl_idx = rng.choice(len(E), size=budget, p=r["p_init"])
                vqe_idx = rng.choice(len(E), size=budget, p=r["p_final"])
                out[key] = {"control": st.set_readouts(ctrl_idx, E),
                            "vqe": st.set_readouts(vqe_idx, E),
                            "final": r["final"],
                            "secs": time.time() - t0}
                G.ck(TAG, tag, out)
                c, v = out[key]["control"], out[key]["vqe"]
                print(f"{key:20s} argmin {c['argmin_rmsd']:6.3f}->{v['argmin_rmsd']:6.3f}  "
                      f"top20mean {c['top20_mean_rmsd']:6.3f}->{v['top20_mean_rmsd']:6.3f}  "
                      f"top20avg {c['top20_coordavg_rmsd']:6.3f}->"
                      f"{v['top20_coordavg_rmsd']:6.3f}  ({out[key]['secs']:.0f}s)")
    summarise_control(out, alphas)
    return out


def summarise_control(out, alphas=ALPHAS):
    print()
    print("=" * 112)
    print("D2b. PAIRED, per readout and per alpha.  NEGATIVE = VQE better than the control.")
    print("=" * 112)
    keys = ("argmin_rmsd", "top5_mean_rmsd", "top20_mean_rmsd", "top75_mean_rmsd",
            "top5_coordavg_rmsd", "top20_coordavg_rmsd", "top75_coordavg_rmsd",
            "set_mean_rmsd", "set_best_rmsd")
    res = {}
    for ro in keys:
        print(f"\n  --- {ro} ---")
        print(f"  {'alpha':>6s} {'n':>4s} {'VQE - control':>14s} {'CI95':>24s} "
              f"{'W/L':>7s} {'mean/sd':>8s} {'verdict':>12s}")
        for a in alphas:
            va, ca = [], []
            for k, v in out.items():
                if f"|a{a}|" not in k or ro not in v["vqe"]:
                    continue
                va.append(v["vqe"][ro]); ca.append(v["control"][ro])
            pr = G.paired(va, ca)
            res[f"{ro}|a{a}"] = pr
            if pr.get("n", 0) == 0:
                continue
            print(f"  {a:6.2f} {pr['n']:4d} {pr['mean']:+14.4f} "
                  f"{f'[{pr['ci_lo']:+.4f}, {pr['ci_hi']:+.4f}]':>24s} "
                  f"{f'{pr['win']}/{pr['loss']}':>7s} {pr['mean_over_sd']:8.3f} "
                  f"{pr['verdict']:>12s}")
    G.ck(TAG, "D2b_paired", res)
    return res


# =========================================================================== D3
def regime_sweep(pdbs=("1CS9", "2MK7", "2P5H", "6EY3", "7N2I", "9UV5"),
                 signals=(0.0, 0.1, 0.3, 0.7), alphas=ALPHAS, iters=150,
                 seeds=(0, 1), budget=2048, residues=(1, 2, 3, 4, 5, 6)):
    """Objective quality x alpha x readout: is there ANY regime where tail-shaping wins?"""
    print()
    print("=" * 112)
    print("D3. THE REGIME SWEEP.  ORACLE DIAGNOSTIC objective-quality knob (blend), so the")
    print("    question 'at what objective quality does tail-shaping pay?' is askable.")
    print("=" * 112)
    prev = G.ck_load(TAG)
    out = prev.get("D3_regime", {})
    for pdb in pdbs:
        st = Struct(pdb, residues)
        circ = make_circ(st.n_qubits)
        for sig in signals:
            E = G.V.blend_objective(st.e.legacy[st.full], st.rmsd, sig)
            for a in alphas:
                for s in seeds:
                    key = f"{pdb}|sig{sig}|a{a}|s{s}"
                    if key in out:
                        continue
                    rng = np.random.default_rng(1000 * s + 7)
                    r = run_vqe(circ, E, a, iters, seed=s, lowmem=True)
                    ci = rng.choice(len(E), size=budget, p=r["p_init"])
                    vi = rng.choice(len(E), size=budget, p=r["p_final"])
                    out[key] = {"control": st.set_readouts(ci, E, ms=(1, 20, 75)),
                                "vqe": st.set_readouts(vi, E, ms=(1, 20, 75)),
                                "entropy_bits": r["final"]["entropy_bits"]}
                    G.ck(TAG, "D3_regime", out)
            print(f"  {pdb} signal {sig} done", flush=True)
    summarise_regime(out, signals, alphas)
    return out


def summarise_regime(out, signals, alphas):
    print()
    print("=" * 112)
    print("D3b. WHERE, IF ANYWHERE, DOES TAIL-SHAPING PAY?")
    print("     Two questions per cell: (i) does alpha<1 beat alpha=1?  (ii) does it beat")
    print("     the best-of-N control?  A regime for CVaR needs BOTH.")
    print("=" * 112)
    res = {}
    for ro in ("argmin_rmsd", "top20_mean_rmsd", "top20_coordavg_rmsd",
               "top75_mean_rmsd", "top75_coordavg_rmsd"):
        print(f"\n  --- {ro} ---")
        print(f"  {'signal':>7s} " + "".join(f"{f'a={a}':>22s}" for a in alphas))
        for sig in signals:
            cells = []
            for a in alphas:
                v = [x["vqe"][ro] for k, x in out.items()
                     if f"|sig{sig}|a{a}|" in k and ro in x["vqe"]]
                c = [x["control"][ro] for k, x in out.items()
                     if f"|sig{sig}|a{a}|" in k and ro in x["control"]]
                if not v:
                    cells.append("--")
                    continue
                pr = G.paired(v, c)
                res[f"{ro}|sig{sig}|a{a}"] = {"vqe_mean": float(np.mean(v)),
                                              "ctrl_mean": float(np.mean(c)), **pr}
                cells.append(f"{np.mean(v):.3f}/{np.mean(c):.3f}({pr['win']}/{pr['loss']})")
            print(f"  {sig:7.2f} " + "".join(f"{x:>22s}" for x in cells))
    print("\n  Each cell is  VQE / control (W/L).  W = VQE better.")
    G.ck(TAG, "D3b_regime_paired", res)
    return res


# =========================================================================== D4
def constraint_role(pdbs=("1CS9", "2MK7", "2P5H", "6EY3", "7N2I", "9UV5"),
                    iters=200, seeds=(0, 1, 2), budget=4096, mus=(0.5, 2.0),
                    residues=(1, 2, 3, 4, 5, 6)):
    """CVaR as a RISK MEASURE on a constraint term -- the use it was invented for.

    C = E_p[disto]  +  mu * RISK(steric),   RISK in {E_p, CVaR_alpha at the UPPER tail}.
    The upper tail of a clash term is what a constraint cares about; the expectation is not.
    Implemented by feeding `-steric` to the lower-tail CVaR and negating, so the same
    verified estimator is used.
    """
    print()
    print("=" * 112)
    print("D4. CVaR AS A CONSTRAINT, NOT AS THE OBJECTIVE.")
    print("    minimise  E[distogram]  +  mu * RISK(steric clash)")
    print("    RISK = expectation (the incumbent) vs CVaR of the UPPER tail (a risk measure)")
    print("=" * 112)
    from s14.vqe_hamil import tabulate
    prev = G.ck_load(TAG)
    out = prev.get("D4_constraint", {})
    for pdb in pdbs:
        st = Struct(pdb, residues)
        _pr, ds = tabulate(pdb)
        D_ = G.V.uniformise(ds[st.full])
        Cl = G.V.uniformise(st.e.leg["steric"][st.full])
        thr = float(np.quantile(Cl, 0.90))            # "violating" = worst decile of clash
        circ = make_circ(st.n_qubits)
        for mu in mus:
            for risk in ("expectation", "cvar0.10", "cvar0.02"):
                for s in seeds:
                    key = f"{pdb}|mu{mu}|{risk}|s{s}"
                    if key in out:
                        continue
                    th = G.random_theta(circ, s)
                    opt = Adam(circ.n_params(), lr=0.10)
                    p0 = circ.probs(th).copy()
                    for _ in range(iters):
                        psi = circ.state(th)
                        p = psi ** 2
                        p = p / p.sum()
                        if risk == "expectation":
                            wr = Cl
                        else:
                            a = float(risk.split("cvar")[1])
                            wr = -cvar_exact(-Cl, p, a)[2]
                        wtot = D_ + mu * wr
                        th = opt.step(th, G.grad_adjoint(circ, th, 2.0 * psi * wtot))
                    p = circ.probs(th)
                    rng = np.random.default_rng(1000 * s + 7)
                    idx = rng.choice(len(D_), size=budget, p=p)
                    ci = rng.choice(len(D_), size=budget, p=p0)
                    sel = D_ + mu * Cl
                    ro = st.set_readouts(idx, sel, ms=(1, 20, 75))
                    ro_c = st.set_readouts(ci, sel, ms=(1, 20, 75))
                    out[key] = {"vqe": ro, "control": ro_c,
                                "violate_mass": float(p[Cl > thr].sum()),
                                "violate_mass_init": float(p0[Cl > thr].sum()),
                                "entropy_bits": float(-(p[p > 0]
                                                        * np.log2(p[p > 0])).sum())}
                    G.ck(TAG, "D4_constraint", out)
        print(f"  {pdb} done", flush=True)
    summarise_constraint(out, mus)
    return out


def summarise_constraint(out, mus=(0.5, 2.0)):
    print()
    print("=" * 112)
    print("D4b. DOES A RISK MEASURE ON THE CONSTRAINT BEAT AN EXPECTATION ON IT?")
    print("=" * 112)
    res = {}
    for mu in mus:
        print(f"\n  --- mu = {mu} ---")
        print(f"  {'risk':>12s} {'violate mass':>13s} {'H bits':>8s} {'argmin':>8s} "
              f"{'top20mean':>10s} {'top20avg':>9s} {'vs ctrl top20avg':>17s}")
        for risk in ("expectation", "cvar0.10", "cvar0.02"):
            rows = [v for k, v in out.items() if f"|mu{mu}|{risk}|" in k]
            if not rows:
                continue
            f = lambda kk, src="vqe": float(np.mean([r[src][kk] for r in rows]))
            pr = G.paired([r["vqe"]["top20_coordavg_rmsd"] for r in rows],
                          [r["control"]["top20_coordavg_rmsd"] for r in rows])
            res[f"mu{mu}|{risk}"] = {"violate_mass": float(np.mean(
                [r["violate_mass"] for r in rows])),
                "entropy_bits": float(np.mean([r["entropy_bits"] for r in rows])),
                "argmin": f("argmin_rmsd"), "top20_mean": f("top20_mean_rmsd"),
                "top20_coordavg": f("top20_coordavg_rmsd"), "vs_control": pr}
            print(f"  {risk:>12s} {np.mean([r['violate_mass'] for r in rows]):13.4f} "
                  f"{np.mean([r['entropy_bits'] for r in rows]):8.2f} "
                  f"{f('argmin_rmsd'):8.3f} {f('top20_mean_rmsd'):10.3f} "
                  f"{f('top20_coordavg_rmsd'):9.3f} "
                  f"{pr['mean']:+8.3f} {pr['verdict']:>8s}")
    G.ck(TAG, "D4b_constraint_summary", res)
    return res


SUB = (1, 2, 3, 4, 5, 6)


def main():
    G.wait_mem(1.0, "qgeom_cvar")
    G.ck_load(TAG)
    t0 = time.time()
    trajectory(residues=SUB)
    control_test(residues=SUB, tag="D2_control_sub12")
    regime_sweep(residues=SUB)
    constraint_role(residues=SUB)
    control_test(pdbs=TARGETS, alphas=(1.0, 0.25, 0.05), seeds=(0,), iters=120,
                 residues=None, tag="D2_control_full18")
    print(f"\ndone in {time.time()-t0:.0f}s -> s15/results/qgeom_cvar.json")


if __name__ == "__main__":
    main()
