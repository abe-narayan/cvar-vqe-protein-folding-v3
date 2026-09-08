"""SPRINT 15 / QGEOM -- PART A(i): the metric, verified, and the Sprint 14 table corrected.

WHAT THIS MODULE ESTABLISHES

V.  Verification.  Exact analytic derivatives vs Sprint 14's central differences; the
    classical Fisher information of the measured distribution vs the Fubini-Study metric;
    the uniform reference point.

A1. THE SPRINT 14 CONDITION NUMBERS ARE SEED-AVERAGED AND ARE NOT CONDITION NUMBERS OF ANY
    METRIC.  `s14/vqe_ansatz.metric_study` averages g over 6 seeds and then takes the
    condition number of the AVERAGE.  g depends on theta; the average of six metrics at six
    different points is not the metric at any point.  Recomputed per point, the spread is
    enormous and the recorded numbers are wrong in BOTH directions.

A2. THE GRADIENT LIES EXACTLY IN THE RANGE OF THE METRIC -- a theorem, then measured.
    For any cost that depends on theta only through the state,
        grad_j C = <d_j psi | dC/dpsi>  =  (D u)_j,     u = dC/dpsi
    and, since <d_j psi|psi> = 0 identically for a normalised real state,
        g = D D^T,   range(g) = range(D) ,   so   grad C in range(g)  ALWAYS.
    A rank-deficient metric therefore has an INFINITE condition number that costs nothing:
    the flat directions are exactly the directions the gradient never points along.
    This turns Sprint 14's HYPOTHESIS ("the badly-conditioned directions are the reducible
    ones") into a proof for the exactly-null directions, and leaves the small-but-nonzero
    eigendirections as the empirical question, measured here as the gradient's spectral
    profile in the eigenbasis of g.

    python -m s15.qgeom_metric
"""
from __future__ import annotations

import numpy as np

from s15 import qgeom_lib as G

TAG = "metric"
PATTERNS = ("none", "chain", "ring", "brick", "block", "block_chain", "all_to_all")


# ============================================================================ V
def verify(n=8, seeds=4):
    print("=" * 96)
    print("V. VERIFICATION -- three independent checks before any claim is made")
    print("=" * 96)
    out = {}
    d_fd, d_fim, d_uni = [], [], []
    for pat in PATTERNS:
        for L in (1, 2, 3):
            c = G.FlexCircuit(n, L, pat, 2)
            for s in range(seeds):
                th = G.random_theta(c, 10 * s + L)
                g = G.fs_metric(c, th)
                d_fd.append(np.abs(g - G.fs_metric_fd(c, th)).max())
                d_fim.append(np.abs(G.classical_fim(c, th) - 4.0 * g).max())
            d_uni.append(np.abs(c.probs(G.uniform_theta(c)) - 1.0 / c.dim).max())
    out["max_abs_exact_minus_fd"] = float(np.max(d_fd))
    out["max_abs_classicalFIM_minus_4g"] = float(np.max(d_fim))
    out["max_abs_uniform_theta_dev"] = float(np.max(d_uni))
    print(f"  exact analytic metric vs s14 central differences : max abs diff "
          f"{out['max_abs_exact_minus_fd']:.3e}  (h=1e-5 truncation)")
    print(f"  classical Fisher information vs 4 * Fubini-Study : max abs diff "
          f"{out['max_abs_classicalFIM_minus_4g']:.3e}  <-- EXACT IDENTITY")
    print(f"  uniform_theta gives the uniform distribution     : max abs dev "
          f"{out['max_abs_uniform_theta_dev']:.3e}")
    print()
    print("  F = 4g is a THEOREM here, not a coincidence.  Amplitudes are real, so")
    print("  p = psi^2, dp_j = 2 psi d_j psi, and F_ij = sum_x dp_i dp_j / p = 4 D D^T,")
    print("  while <d_j psi|psi> = 0 for a normalised state makes g = D D^T.  Therefore")
    print("  the metric a SAMPLING optimiser of a DIAGONAL cost lives on is the quantum")
    print("  geometric tensor up to a factor 4, absorbed into the learning rate:")
    print("  QNG and classical natural gradient are THE SAME ALGORITHM in this setting.")
    G.ck(TAG, "verify", out)
    return out


# =========================================================================== A1
def metric_table(n=10, seeds=12):
    """Per-point metric statistics, with the seed-averaged artefact printed alongside."""
    print()
    print("=" * 96)
    print("A1. THE METRIC PER ENTANGLER AND DEPTH -- PER POINT, not seed-averaged")
    print("    s14 averaged g over seeds then took cond(mean g).  That is not the")
    print("    condition number of the metric at any point.  Both are printed.")
    print("=" * 96)
    print(f"{'pattern':12s} {'L':>2s} {'P':>3s} {'rank':>5s} {'g_ii':>8s} "
          f"{'cond|range median':>18s} {'[min,max]':>22s} {'s14-style cond(mean g)':>23s}")
    out = {}
    for pat in PATTERNS:
        for L in (1, 2, 3):
            c = G.FlexCircuit(n, L, pat, 2)
            gs, cs, rk, gii = [], [], [], []
            for s in range(seeds):
                th = np.random.default_rng(100 * s + 7).normal(0, 0.9, c.n_params())
                g = G.fs_metric(c, th)
                st = G.spec_stats(g)
                gs.append(g)
                cs.append(st["cond"])
                rk.append(st["rank"])
                gii.append(st["mean_gii"])
            gm = np.mean(gs, axis=0)
            s14c = G.spec_stats(gm)["cond_full"]
            row = {"n_params": c.n_params(), "rank_median": int(np.median(rk)),
                   "rank_min": int(np.min(rk)), "rank_max": int(np.max(rk)),
                   "mean_gii": float(np.mean(gii)),
                   "cond_range_median": float(np.median(cs)),
                   "cond_range_min": float(np.min(cs)), "cond_range_max": float(np.max(cs)),
                   "cond_of_seed_averaged_metric_S14_STYLE": float(s14c),
                   "seeds": seeds}
            out[f"{pat}|L{L}"] = row
            print(f"{pat:12s} {L:2d} {c.n_params():3d} {int(np.median(rk)):5d} "
                  f"{np.mean(gii):8.6f} {np.median(cs):18.2f} "
                  f"{f'[{np.min(cs):.1f}, {np.max(cs):.1f}]':>22s} {s14c:23.3f}")
            G.ck(TAG, "A1_table", out)
    print()
    print("  `cond|range` is the condition number over the NON-NULL eigendirections;")
    print("  the null directions are exactly gradient-free (A2) so they cannot hurt.")
    print("  `rank` < P is exact parameter redundancy.")
    return out


def rank_law(n=12, seeds=6):
    """The block ansatz's rank is set by the block Hilbert dimension, not by P."""
    print()
    print("=" * 96)
    print("A1b. THE RANK OF A BLOCK-FACTORISED ANSATZ IS CAPPED BY THE BLOCK, NOT BY P")
    print("=" * 96)
    print(f"{'pattern':12s} {'b':>2s} {'L':>2s} {'P':>3s} {'rank':>5s} "
          f"{'predicted':>10s} {'match':>6s}")
    out = {}
    for b in (2, 3):
        nb = n // b * b
        for L in (1, 2, 3, 4):
            for pat in ("none", "block"):
                c = G.FlexCircuit(nb, L, pat, b)
                rks = [G.spec_stats(G.fs_metric(c, G.random_theta(c, s)))["rank"]
                       for s in range(seeds)]
                nblocks = nb // b if pat == "block" else nb
                dof = (2 ** b - 1) if pat == "block" else 1
                pred = min(c.n_params(), nblocks * dof)
                ok = int(np.median(rks)) == pred
                out[f"{pat}|b{b}|L{L}"] = {"P": c.n_params(), "rank": int(np.median(rks)),
                                           "predicted": int(pred), "match": bool(ok)}
                print(f"{pat:12s} {b:2d} {L:2d} {c.n_params():3d} "
                      f"{int(np.median(rks)):5d} {pred:10d} {str(ok):>6s}")
    print()
    print("  A `block` ansatz's state is a PRODUCT over blocks, so its reachable manifold")
    print("  has at most (2^b - 1) real dimensions per block.  Parameters beyond that are")
    print("  EXACTLY redundant, the metric is EXACTLY singular, and the infinite condition")
    print("  number s14 reported for `none` is the same phenomenon one level up.")
    G.ck(TAG, "A1b_rank_law", out)
    return out


# =========================================================================== A2
def gradient_in_range(n=10, seeds=8):
    """THE THEOREM, measured: grad C lies in range(g); and where in the spectrum it sits."""
    print()
    print("=" * 96)
    print("A2. THE GRADIENT LIES IN THE RANGE OF THE METRIC (theorem), AND ITS SPECTRAL")
    print("    PROFILE (measurement).  If the gradient avoided the small eigendirections,")
    print("    QNG could buy nothing.  If it lives there, QNG has something to fix.")
    print("=" * 96)
    e = G.V.Enum("1CS9")
    E = G.V.uniformise(e.legacy)                      # native-free, conditioned
    print(f"{'pattern':12s} {'L':>2s} {'rank':>5s} {'|g_null|/|g|':>13s} "
          f"{'top-decile eig share':>21s} {'bottom-decile share':>20s} {'cond|range':>11s}")
    out = {}
    for pat in PATTERNS:
        for L in (1, 2, 3):
            c = G.FlexCircuit(n, L, pat, 2)
            # a matched sub-objective on the first n qubits of the enumerated register
            Esub = E[: c.dim] if E.size >= c.dim else np.resize(E, c.dim)
            Esub = G.V.uniformise(Esub)
            nulls, tops, bots, conds, rks = [], [], [], [], []
            for s in range(seeds):
                th = G.random_theta(c, 100 * s + 7)
                psi = c.state(th)
                D = G.deriv_states(c, th)
                g = D @ D.T - np.outer(D @ psi, psi @ D.T if D.ndim == 2 else 0)
                g = (g + g.T) / 2
                # exact gradient of the EXPECTATION value of Esub
                grad = 2.0 * (D @ (psi * Esub))
                w, Vv = np.linalg.eigh(g)
                w = np.clip(w, 0, None)
                comp = (Vv.T @ grad) ** 2
                tot = comp.sum()
                tol = 1e-9 * max(w.max(), 1e-30)
                nullmask = w <= tol
                nulls.append(float(comp[nullmask].sum() / max(tot, 1e-300)))
                pos = np.where(~nullmask)[0]
                order = pos[np.argsort(w[pos])]          # ascending eigenvalue
                nd = max(1, len(order) // 10)
                bots.append(float(comp[order[:nd]].sum() / max(tot, 1e-300)))
                tops.append(float(comp[order[-nd:]].sum() / max(tot, 1e-300)))
                conds.append(float(w[pos].max() / w[pos].min()) if len(pos) else np.nan)
                rks.append(int(len(pos)))
            row = {"rank": int(np.median(rks)),
                   "grad_share_in_null_space": float(np.mean(nulls)),
                   "grad_share_top_eig_decile": float(np.mean(tops)),
                   "grad_share_bottom_eig_decile": float(np.mean(bots)),
                   "cond_range_median": float(np.nanmedian(conds))}
            out[f"{pat}|L{L}"] = row
            print(f"{pat:12s} {L:2d} {row['rank']:5d} "
                  f"{row['grad_share_in_null_space']:13.3e} "
                  f"{row['grad_share_top_eig_decile']:21.4f} "
                  f"{row['grad_share_bottom_eig_decile']:20.4f} "
                  f"{row['cond_range_median']:11.2f}")
            G.ck(TAG, "A2_gradient_in_range", out)
    print()
    print("  `|g_null|/|g|` is the fraction of the gradient's squared norm lying in the")
    print("  metric's EXACT null space.  The theorem says 0; the measurement is the check.")
    print("  A uniform gradient would put 0.10 in each decile.")
    return out


def main():
    G.wait_mem(0.8, "qgeom_metric")
    G.ck_load(TAG)
    verify()
    metric_table()
    rank_law()
    gradient_in_range()
    print("\nwritten -> s15/results/qgeom_metric.json")


if __name__ == "__main__":
    main()
