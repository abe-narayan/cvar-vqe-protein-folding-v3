"""SPRINT 14 / VQE -- PART D: the ansatz as a research variable, and initialisation.

Families compared, all simulated EXACTLY (real-amplitude RY + CNOT, complex only for the
QAOA arm):

    none          RY layers, no entangler.  The product-distribution control that says
                  whether the entanglement does anything at all.
    chain         RY + nearest-neighbour CNOT chain.  The shipped hardware-efficient ansatz.
    ring          chain plus the wrap-around CNOT.  The shipped default.
    brick         even pairs then odd pairs -- the same connectivity at half the depth.
    block         CNOTs only INSIDE each residue's log2(k)-qubit block.  Torsion-aware and
                  deliberately reducible: the distribution factorises over residues.
    block_chain   inside-block CNOTs plus ONE link between consecutive residue blocks.
                  This is the "problem-inspired / geometry-aware" arm: it entangles along
                  the covalent chain, which is the axis the locality theorem says matters.
    all_to_all    every ordered pair.  The expressivity upper bound at this width.

MEASURED: Fubini-Study metric spectrum and condition number, parameter redundancy (metric
rank), gradient variance vs width (barren plateau), expressivity against a deliberately
multimodal target with SEEDS (s13 left this at one seed per cell and classified it OPEN),
convergence, and end-to-end structural accuracy at matched budget.

INITIALISATION IS TREATED AS A CONFOUND, NOT A FEATURE.  Every arm reports the RMSD of its
OWN initialisation next to the RMSD after optimisation.  An initialisation already at the
answer is not a VQE result and the table is built so that cannot be hidden.

    python -m s14.vqe_ansatz
"""
from __future__ import annotations

import numpy as np

from core import quantum as Q
from s14 import vqe_lib as V
from s14 import vqe_run as R


# =================================================== a configurable exact circuit
class FlexCircuit:
    """``layers x (RY on every wire, a fixed CNOT pattern)``, exact real statevector.

    Identical arithmetic to `core.quantum.StatevectorCircuit` -- the entangler is a fixed
    permutation of basis indices composed once -- but the pattern is a parameter.
    """

    def __init__(self, n, layers, pattern="ring", bits_per_res=2, final_ry=False):
        self.n, self.layers = int(n), int(layers)
        self.pattern, self.b = pattern, int(bits_per_res)
        self.final_ry = bool(final_ry)
        self.dim = 1 << self.n
        self.pairs = self._pairs()
        self._perm = self._permutation()

    def _pairs(self):
        n, b = self.n, self.b
        p = self.pattern
        if p == "none":
            return []
        if p == "chain":
            return [(q, q + 1) for q in range(n - 1)]
        if p == "ring":
            return [(q, q + 1) for q in range(n - 1)] + ([(n - 1, 0)] if n > 2 else [])
        if p == "brick":
            return ([(q, q + 1) for q in range(0, n - 1, 2)]
                    + [(q, q + 1) for q in range(1, n - 1, 2)])
        if p == "block":
            out = []
            for s in range(0, n - b + 1, b):
                out += [(s + i, s + i + 1) for i in range(b - 1)]
            return out
        if p == "block_chain":
            out = []
            for s in range(0, n - b + 1, b):
                out += [(s + i, s + i + 1) for i in range(b - 1)]
            out += [(s + b - 1, s + b) for s in range(0, n - b, b)]
            return out
        if p == "all_to_all":
            return [(i, j) for i in range(n) for j in range(n) if i != j]
        raise ValueError(p)

    def _permutation(self):
        n, dim = self.n, self.dim
        cur = np.arange(dim, dtype=np.int64)
        idx = np.arange(dim, dtype=np.int64)
        for c, t in self.pairs:
            cb = (idx >> (n - 1 - c)) & 1
            cur = cur[np.where(cb == 1, idx ^ (1 << (n - 1 - t)), idx)]
        return cur

    def n_params(self):
        return self.n * (self.layers + int(self.final_ry))

    def states_batch(self, thetas):
        TH = np.atleast_2d(np.asarray(thetas, float))
        B = TH.shape[0]
        nb = self.layers + int(self.final_ry)
        TH = TH.reshape(B, nb, self.n)
        CO, SI = np.cos(TH / 2.0), np.sin(TH / 2.0)
        psi = np.zeros((B, self.dim))
        psi[:, 0] = 1.0
        for L in range(nb):
            for q in range(self.n):
                c = CO[:, L, q][:, None, None]
                s = SI[:, L, q][:, None, None]
                v = psi.reshape(B, 1 << q, 2, 1 << (self.n - q - 1))
                a = v[:, :, 0, :].copy()
                bb = v[:, :, 1, :]
                v[:, :, 0, :] = c * a - s * bb
                v[:, :, 1, :] = s * a + c * bb
                psi = v.reshape(B, self.dim)
            if L < self.layers and self.pairs:
                psi = psi[:, self._perm]
        return psi

    def state(self, th):
        return self.states_batch(np.asarray(th, float)[None, :])[0]

    def probs(self, th):
        p = self.state(th) ** 2
        return p / p.sum()

    def probs_batch(self, TH):
        p = self.states_batch(TH) ** 2
        return p / p.sum(1, keepdims=True)

    def _shift_grid(self, th, h):
        th = np.asarray(th, float)
        P = th.size
        TH = np.repeat(th[None, :], 2 * P, axis=0)
        r = np.arange(P)
        TH[2 * r, r] += h
        TH[2 * r + 1, r] -= h
        return TH


PATTERNS = ("none", "chain", "ring", "brick", "block", "block_chain", "all_to_all")


# ---------------------------------------------------------------- Fubini-Study
def fubini_study(circ, theta, h=1e-5):
    """``g_ij = Re<d_i psi|d_j psi> - <d_i psi|psi><psi|d_j psi>`` by central differences.

    Real amplitudes, so the imaginary (Berry) part vanishes identically and the metric is
    the real symmetric part.  Independent of `core.quantum`'s machinery on purpose.
    """
    th = np.asarray(theta, float)
    P = th.size
    psi = circ.state(th)
    D = np.empty((P, circ.dim))
    for k in range(P):
        tp = th.copy(); tp[k] += h
        tm = th.copy(); tm[k] -= h
        D[k] = (circ.state(tp) - circ.state(tm)) / (2 * h)
    ov = D @ psi
    return D @ D.T - np.outer(ov, ov)


def metric_study(n=10, bits_per_res=2, seeds=6):
    print("=" * 100)
    print("D1. FUBINI-STUDY METRIC per entangler pattern and depth")
    print("    The brief records the metric is EXACTLY I/4 at depth 1 for the shipped")
    print("    ansatz and contains no Hamiltonian.  Re-tested here for TORSION-AWARE")
    print("    patterns: does entangling along the residue blocks change the geometry?")
    print("=" * 100)
    print(f"{'pattern':13s} {'L':>2s} {'P':>3s} {'mean g_ii':>10s} {'max|g_ij| off':>14s} "
          f"{'cond number':>12s} {'rank':>5s} {'||g - I/4||':>12s}")
    out = {}
    for pat in PATTERNS:
        for L in (1, 2, 3):
            gs = []
            for s in range(seeds):
                r = np.random.default_rng(100 * s + 7)
                c = FlexCircuit(n, L, pat, bits_per_res)
                gs.append(fubini_study(c, r.normal(0, 0.9, c.n_params())))
            g = np.mean(gs, axis=0)
            P = g.shape[0]
            off = g - np.diag(np.diag(g))
            ev = np.linalg.eigvalsh(np.mean(gs, axis=0))
            cond = float(ev.max() / max(ev.min(), 1e-15)) if ev.min() > 1e-12 else np.inf
            rank = int(np.linalg.matrix_rank(g, tol=1e-9))
            dI = float(np.abs(g - np.eye(P) / 4).max())
            out.setdefault(pat, {})[L] = {
                "mean_gii": float(np.mean(np.diag(g))),
                "max_offdiag": float(np.abs(off).max()),
                "cond": cond, "rank": rank, "n_params": P,
                "max_dev_from_I_over_4": dI}
            print(f"{pat:13s} {L:2d} {P:3d} {np.mean(np.diag(g)):10.6f} "
                  f"{np.abs(off).max():14.6f} {cond:12.3f} {rank:5d} {dI:12.2e}")
    print("\n  A metric equal to I/4 means the natural gradient IS the ordinary gradient")
    print("  and QNG can buy nothing.  Deviation from I/4 is what would make QNG useful.")
    return out


def metric_contains_no_hamiltonian(n=10):
    print()
    print("=" * 100)
    print("D1b. DOES THE METRIC DEPEND ON THE ENERGY MODEL?  (re-test at the block ansatz)")
    print("=" * 100)
    e = V.Enum("1CS9")
    r = np.random.default_rng(0)
    for pat in ("ring", "block_chain"):
        c = FlexCircuit(n, 2, pat, 2)
        th = r.normal(0, 0.9, c.n_params())
        g = fubini_study(c, th)
        print(f"  {pat:13s} metric is a function of theta only; max |g(E1) - g(E2)| = "
              f"{0.0:.3e} by construction (no E enters `fubini_study`).")
    print("  CONFIRMED STRUCTURALLY, not merely numerically: the Fubini-Study metric is")
    print("  built from the STATE and its derivatives alone.  No energy model can appear.")
    print("  The energy selects WHERE on the manifold you go, not the manifold's shape.")


# ------------------------------------------------------------- barren plateaus
def barren(seeds=80):
    print()
    print("=" * 100)
    print("D2. GRADIENT VARIANCE vs WIDTH per pattern (barren-plateau behaviour)")
    print("    Var of one parameter's gradient of the EXPECTATION value (alpha=1),")
    print("    conditioned objective, random angles uniform on [0, 2pi).")
    print("=" * 100)
    print(f"{'pattern':13s} " + "".join(f"{f'n={n}':>12s}" for n in (6, 8, 10, 12))
          + f"{'decay/qubit':>13s}")
    out = {}
    for pat in PATTERNS:
        vs = []
        for n in (6, 8, 10, 12):
            r0 = np.random.default_rng(n)
            E = V.uniformise(r0.standard_normal(1 << n))
            c = FlexCircuit(n, 3, pat, 2)
            g1 = []
            for s in range(seeds):
                r = np.random.default_rng(1000 * s + n)
                th = r.uniform(0, 2 * np.pi, c.n_params())
                dp = Q.cvar_exact(E, c.probs(th), 1.0)[2]
                PR = c.probs_batch(c._shift_grid(th, np.pi / 2))
                g1.append(float(((PR[0::2] - PR[1::2]) @ dp / 2.0)[0]))
            vs.append(float(np.var(g1)))
        sl = np.polyfit([6, 8, 10, 12], np.log2(vs), 1)[0]
        out[pat] = {"vars": vs, "log2_slope_per_qubit": float(sl)}
        print(f"{pat:13s} " + "".join(f"{v:12.3e}" for v in vs) + f"{sl:13.3f}")
    print("\n  A slope of -1.0 per qubit is the textbook barren plateau (Var ~ 2^-n).")
    return out


# ---------------------------------------------------------------- expressivity
def expressivity(pdb="1CS9", n_res=6, seeds=6, iters=300):
    print()
    print("=" * 100)
    print(f"D3. EXPRESSIVITY -- fit a deliberately MULTIMODAL target ({seeds} seeds)")
    print("    Target: uniform over the 8 lowest-CA-RMSD configurations of a real")
    print("    12-qubit sub-register.  ORACLE DIAGNOSTIC in its TARGET only; the measured")
    print("    quantity is a property of the ansatz.  s13 ran ONE seed and classified this")
    print("    OPEN -- this closes it.")
    print("=" * 100)
    e = V.Enum(pdb)
    residues = tuple(range(1, 1 + n_res))
    base = np.random.default_rng(0).integers(0, e.k, e.n)
    from s14.vqe_encoding import sub_objective
    r = sub_objective(e, residues, base, "rmsd")
    modes = np.argsort(r)[:8]
    tgt = np.zeros(len(r)); tgt[modes] = 1 / 8
    nq = n_res * 2
    print(f"{'pattern':13s} " + "".join(f"{f'L={L}':>10s}" for L in (1, 2, 3, 4, 6))
          + f"{'  best modes captured':>22s}")
    out = {}
    for pat in PATTERNS:
        row, bestj = [], 0
        for L in (1, 2, 3, 4, 6):
            kls = []
            for s in range(seeds):
                rg = np.random.default_rng(7 * s + L)
                c = FlexCircuit(nq, L, pat, 2)
                th = rg.normal(0, 0.8, c.n_params())
                opt = Q.Adam(c.n_params(), lr=0.1)
                for _ in range(iters):
                    p = np.maximum(c.probs(th), 1e-15)
                    # d KL(tgt||p) / dp = -tgt/p
                    PR = c.probs_batch(c._shift_grid(th, np.pi / 2))
                    g = (PR[0::2] - PR[1::2]) @ (-tgt / p) / 2.0
                    th = opt.step(th, g)
                p = np.maximum(c.probs(th), 1e-300)
                kls.append(float((tgt[modes] * np.log(tgt[modes] / p[modes])).sum()))
            row.append(float(np.median(kls)))
            bestj = max(bestj, int(round(8 * np.exp(-min(kls)))))
        out[pat] = {"kl_median": row, "layers": [1, 2, 3, 4, 6]}
        print(f"{pat:13s} " + "".join(f"{v:10.4f}" for v in row)
              + f"{min(8, bestj):22d}")
    print("\n  KL = ln(8/j) exactly when the circuit has put its mass on j of the 8 modes:")
    print("  2.079 = 1 mode, 1.386 = 2, 0.693 = 4, 0.000 = all 8.  MEDIAN over seeds.")
    return out


# ------------------------------------------------------ end to end + initialisation
def end_to_end(pdb="1CS9", signals=(0.0, 0.1, 0.3, 1.0), budget=20480, seeds=(0, 1, 2)):
    print()
    print("=" * 100)
    print(f"D4. END TO END on {pdb} at matched budget ({budget} objective evaluations)")
    print("    ansatz family x signal level, and INITIALISATION RMSD next to the outcome.")
    print("=" * 100)
    e = V.Enum(pdb)
    ans = ["onelayer", "mps1f", "mps2f", "mps3f", "mps2", "mps2n"]
    out = {}
    for sig in signals:
        E = V.blend_objective(e.legacy, e.rmsd, sig)
        ex = int(np.argmin(E))
        print(f"\n--- signal {sig} (rho {V.spearman(E, e.rmsd):+.3f}, certified optimum "
              f"{e.rmsd[ex]:.3f} A, space best {e.rmsd.min():.3f}) ---")
        print(f"{'ansatz':10s} {'P':>4s} {'init RMSD':>10s} {'RMSD ret':>9s} {'sd':>6s} "
              f"{'mode RMSD':>10s} {'objgap':>8s} {'H bits':>7s} {'sel gap':>8s}")
        for a in ans:
            rs = [R.run(E, e.n_qubits, 0.25, budget=budget, shots=512, ansatz=a,
                        seed=s, rmsd=e.rmsd) for s in seeds]
            row = {k: float(np.mean([r[k] for r in rs]))
                   for k in ("init_rmsd_mean", "rmsd_returned", "mode_rmsd",
                             "entropy_bits", "rmsd_best_seen", "best_e")}
            row["rmsd_sd"] = float(np.std([r["rmsd_returned"] for r in rs]))
            row["objective_gap"] = row["best_e"] - E[ex]
            row["selection_gap"] = row["rmsd_returned"] - row["rmsd_best_seen"]
            row["n_params"] = rs[0]["n_params"]
            out.setdefault(str(sig), {})[a] = row
            print(f"{a:10s} {row['n_params']:4d} {row['init_rmsd_mean']:10.3f} "
                  f"{row['rmsd_returned']:9.3f} {row['rmsd_sd']:6.3f} "
                  f"{row['mode_rmsd']:10.3f} {row['objective_gap']:8.4f} "
                  f"{row['entropy_bits']:7.2f} {row['selection_gap']:8.3f}")
    return out


def initialisation(pdb="1CS9", signal=0.1, budget=20480, seeds=(0, 1, 2, 3, 4)):
    print()
    print("=" * 100)
    print("D5. INITIALISATION: acceleration, or handing over the answer?")
    print("    THE RULE: an initialisation already at the answer is not a VQE result.")
    print("    So every row prints the RMSD of its OWN initialisation first.")
    print("=" * 100)
    from s13.qarch_lib import Space, empirical_prior, ORACLE_prior
    e = V.Enum(pdb)
    E = V.blend_objective(e.legacy, e.rmsd, signal)
    ex = int(np.argmin(E))
    sp = Space(pdb, k=4)
    Pemp = empirical_prior(sp)                       # leakage-safe, native-free
    arms = {"random (native-free)": ("random", None),
            "uniform/zero (native-free)": ("zero", None),
            "empirical prior warm start (native-free)": ("prior", Pemp)}
    for q in (0.5, 0.8, 0.95):
        arms[f"ORACLE prior q={q} warm start (ORACLE)"] = ("prior", ORACLE_prior(sp, q))
    print(f"signal {signal}: rho {V.spearman(E, e.rmsd):+.3f}, certified optimum "
          f"{e.rmsd[ex]:.3f} A, space best {e.rmsd.min():.3f}, random draw "
          f"{e.rmsd.mean():.3f}")
    print()
    print(f"{'initialisation':44s} {'INIT RMSD':>10s} {'init best':>10s} "
          f"{'after VQE':>10s} {'mode':>8s} {'DELTA':>8s}")
    out = {}
    for name, (mode, P) in arms.items():
        rs = [R.run(E, e.n_qubits, 0.25, budget=budget, shots=512, ansatz="mps2f",
                    seed=s, rmsd=e.rmsd, init=mode, prior=P, bits_per_res=2)
              for s in seeds]
        i0 = float(np.mean([r["init_rmsd_mean"] for r in rs]))
        ib = float(np.mean([r["init_rmsd_best"] for r in rs]))
        af = float(np.mean([r["rmsd_returned"] for r in rs]))
        md = float(np.mean([r["mode_rmsd"] for r in rs]))
        out[name] = {"init_rmsd_mean": i0, "init_rmsd_best": ib, "after": af,
                     "mode": md, "delta": af - ib}
        print(f"{name:44s} {i0:10.3f} {ib:10.3f} {af:10.3f} {md:8.3f} {af-ib:+8.3f}")
    print("\n  'init best' is the best of the same number of shots drawn from the")
    print("  INITIAL distribution with no optimisation at all -- the honest control.")
    print("  DELTA = after VQE minus that control.  A negative DELTA is the only thing")
    print("  that can be called a contribution from the optimisation.")
    return out


def main():
    V.wait_for_memory(1.0, "vqe_ansatz")
    out = {}
    out["metric"] = metric_study()
    metric_contains_no_hamiltonian()
    out["barren"] = barren()
    out["expressivity"] = expressivity()
    out["end_to_end"] = end_to_end()
    out["initialisation"] = initialisation()
    V.write("vqe_ansatz", out)
    print("\nwritten -> s14/results/vqe_ansatz.json")


if __name__ == "__main__":
    main()
