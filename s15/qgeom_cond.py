"""SPRINT 15 / QGEOM -- PART B: DOES TARGET-SPECIFIC CONDITIONING RESHAPE THE MANIFOLD?

Sprint 14 established that the Fubini-Study metric is bit-identical across energy models at
matched parameters (0.000e+00), and `qgeom_metric` re-establishes it STRUCTURALLY: g is a
function of the state and its derivatives alone, so no energy array can enter it.  The energy
model selects WHICH REGION of the manifold is visited; it does not reshape the manifold.

THE OPEN HALF, and it is this module's whole subject: **conditioning**.  A target-conditioned
VQE does not merely score differently -- it is PREPARED differently.  If target-specific
information drives the state into a geometrically distinguished region (different eigenvalue
spectrum, different effective rank, different metric volume, different trainability), then
the chain

    target-specific information -> state manifold -> metric/curvature -> trainability

is real.  If it does not, that is a clean negative and closes the question.

PRE-REGISTERED, WRITTEN BEFORE THE RUN
--------------------------------------
H0 (null)          the conditioned point is geometrically indistinguishable from a random
                   point at the same depth.
H1 (concentration) the geometry tracks the ENTROPY of whatever distribution the state was
                   prepared for, and NOT its target-specificity.  Conditioned and
                   entropy-matched-scrambled points then agree with each other and both
                   differ from uniform.
H2 (thesis)        target-specific conditioning is geometrically special: the conditioned
                   point differs from an entropy-matched scramble of the SAME prior.
My prior is on H1.  H2 is the paper's centrepiece if it survives; H1 is a clean mechanism;
H0 is a clean negative.  The discriminating control is the ENTROPY-MATCHED SCRAMBLE, which
holds every per-residue distribution fixed and destroys only the assignment to residues /
states.

THE ARMS, all fitted with the identical optimiser, budget and seeds:
    uniform       theta giving the exactly uniform distribution (an analytic point)
    random        a random theta at the shipped init scale
    COND          fitted to the target's own retrieval-conditioned torsion prior
                  (`s14.retprior.state_prior`, the best channel in the project)
    scram_state   fitted to the SAME prior with the k state labels permuted per residue
                  -- identical per-residue entropies, target information destroyed
    scram_res     fitted to the SAME prior with the residues permuted
                  -- identical multiset of distributions, wrong positions
    cross         fitted to a DIFFERENT target's real prior
    dirichlet     fitted to a synthetic product prior with matched total entropy

NATIVE-FREE.  The retrieval prior is a BLOSUM retrieval output.  `Enum.rmsd` is read only
for post-hoc scoring in the trainability table and is labelled there.

    python -m s15.qgeom_cond
"""
from __future__ import annotations

import time

import numpy as np

from s14 import retprior as RP
from s15 import qgeom_lib as G

TAG = "cond"

TARGETS = ("1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5")
ARMS = ("uniform", "random", "COND", "scram_state", "scram_res", "cross", "dirichlet")


# ------------------------------------------------------------------- entropy-matched nulls
def dirichlet_matched(P, seed=0, tries=300):
    """A synthetic product prior whose TOTAL entropy matches P's, with no target content."""
    tgt = G.prior_entropy_bits(P)
    r = np.random.default_rng(seed)
    best, bd = None, np.inf
    lo, hi = 0.02, 20.0
    for _ in range(tries):
        a = np.sqrt(lo * hi)
        Q = r.dirichlet(np.full(P.shape[1], a), size=P.shape[0])
        h = G.prior_entropy_bits(Q)
        if abs(h - tgt) < bd:
            best, bd = Q, abs(h - tgt)
        if h < tgt:
            lo = a
        else:
            hi = a
        lo, hi = max(lo, 0.01), min(hi, 50.0)
        if hi / lo < 1.01:
            lo, hi = 0.02, 20.0                     # restart the bracket, keep sampling
    return best


def arm_targets(pdb, P_by_pdb, seed=0):
    """The seven conditioning targets for one PDB, as full product distributions."""
    P = P_by_pdb[pdb]
    others = [q for q in P_by_pdb if q != pdb]
    Pc = P_by_pdb[others[hash(pdb) % len(others)]]
    if Pc.shape != P.shape:
        Pc = P_by_pdb[others[0]]
    return {
        "COND": P,
        "scram_state": G.scramble_prior(P, seed, "state"),
        "scram_res": G.scramble_prior(P, seed, "residue"),
        "cross": Pc,
        "dirichlet": dirichlet_matched(P, seed),
    }


# ----------------------------------------------------------------------- the measurement
def geometry_at(circ, th, E=None, lowmem=True):
    """Everything the brief asks for at one point of the manifold."""
    g = G.fs_metric(circ, th, lowmem)
    st = G.spec_stats(g)
    psi = circ.state(th)
    p = psi ** 2
    p = p / p.sum()
    st["entropy_bits"] = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
    st["max_p"] = float(p.max())
    if E is not None:
        D = G.deriv_states_lowmem(circ, th) if lowmem else G.deriv_states(circ, th)
        gr = 2.0 * (D @ (psi * E))
        st["grad_norm"] = float(np.linalg.norm(gr))
        st["grad_var"] = float(np.var(gr))
        w, Vv = np.linalg.eigh(g)
        w = np.clip(w, 0, None)
        comp = (Vv.T @ gr) ** 2
        tot = max(comp.sum(), 1e-300)
        tol = 1e-9 * max(w.max(), 1e-30)
        pos = np.where(w > tol)[0]
        order = pos[np.argsort(w[pos])]
        nd = max(1, len(order) // 10)
        st["grad_share_null"] = float(comp[w <= tol].sum() / tot)
        st["grad_share_bottom_decile"] = float(comp[order[:nd]].sum() / tot)
        st["grad_share_top_decile"] = float(comp[order[-nd:]].sum() / tot)
    return st


KEYS = ("cond", "eff_rank_frac", "log_pseudo_det", "mean_gii", "max_dev_from_I4",
        "max_offdiag", "rank", "eig_max", "eig_min_pos", "entropy_bits", "max_p",
        "grad_norm", "grad_var", "grad_share_bottom_decile")


def run_target(pdb, P_by_pdb, pattern="block", L=2, iters=300, seed=0,
               snaps=(0, 25, 75, 150, 299), residues=(1, 2, 3, 4, 5, 6)):
    """`residues=None` uses the FULL register; a tuple uses that sub-register.

    A sub-register is still genuine target-specific conditioning: the prior rows are that
    target's own retrieval posterior for those residues.  It is 64x cheaper than n=18, so
    the statistics are over 9 targets x 7 arms rather than over two anecdotes.
    """
    from s15.qgeom_qng import sub_objective
    e = G.V.Enum(pdb)
    res = tuple(range(e.n)) if residues is None else tuple(residues)
    nq = 2 * len(res)
    circ = G.FlexCircuit(nq, L, pattern, 2)
    base = np.random.default_rng(0).integers(0, e.k, e.n)
    E = G.V.uniformise(sub_objective(e, res, base, "legacy"))     # native-free
    lowmem = nq >= 16
    row = {"n_qubits": nq, "n_params": circ.n_params(), "pattern": pattern, "layers": L,
           "residues": list(res),
           "prior_entropy_bits": G.prior_entropy_bits(P_by_pdb[pdb][list(res)])}
    P_by_pdb = {p: v[list(res)] for p, v in P_by_pdb.items()}
    tg = arm_targets(pdb, P_by_pdb, seed)
    row["arm_entropy_bits"] = {k: G.prior_entropy_bits(v) for k, v in tg.items()}

    row["arms"] = {}
    th_u = G.uniform_theta(circ)
    row["arms"]["uniform"] = {"geom": geometry_at(circ, th_u, E, lowmem), "kl": 0.0}
    th_r = G.random_theta(circ, seed)
    row["arms"]["random"] = {"geom": geometry_at(circ, th_r, E, lowmem),
                             "kl": float("nan")}

    for name, Pk in tg.items():
        t = G.product_target(Pk)
        fit = G.fit_kl(circ, t, theta0=G.random_theta(circ, seed), iters=iters, lr=0.08,
                       lowmem=lowmem, track=snaps)
        row["arms"][name] = {"geom": geometry_at(circ, fit["theta"], E, lowmem),
                             "kl": fit["kl"], "tv": fit["tv"],
                             "kl_history_tail": fit["kl_history"][-1:],
                             "snaps": {str(k): {kk: v[kk] for kk in
                                                ("cond", "eff_rank_frac",
                                                 "log_pseudo_det", "kl")
                                                if kk in v}
                                       for k, v in fit["snaps"].items()}}
    return row


def main(pattern="block", L=2, iters=300, targets=TARGETS, residues=(1, 2, 3, 4, 5, 6),
         tag_suffix="", seeds=(0, 1, 2)):
    G.wait_mem(0.8, "qgeom_cond")
    prev = G.ck_load(TAG)
    key = f"B1_{pattern}L{L}{tag_suffix}"
    done = prev.get(key, {})
    P_by_pdb = {p: RP.state_prior(p, "top75", 4)["P"] for p in targets}
    print("=" * 112)
    print(f"B1. GEOMETRY AT THE CONDITIONED POINT -- ansatz {pattern} L={L}, "
          f"{len(targets)} enumerated targets x {len(seeds)} seeds")
    print("=" * 112)
    print(f"{'target':10s} {'arm':12s} {'KL':>8s} {'H(p) bits':>10s} {'rank':>5s} "
          f"{'cond':>10s} {'effrank':>8s} {'log vol':>9s} {'g_ii':>8s} "
          f"{'|g-I/4|':>8s} {'grad var':>10s} {'g_bot10':>8s}")
    for pdb in targets:
        for sd in seeds:
            cell = f"{pdb}|s{sd}"
            if cell in done:
                continue
            t0 = time.time()
            row = run_target(pdb, P_by_pdb, pattern, L, iters, seed=sd,
                             residues=residues)
            done[cell] = row
            for a in ARMS:
                g = row["arms"][a]["geom"]
                kl = row["arms"][a]["kl"]
                print(f"{cell:10s} {a:12s} {kl:8.4f} {g['entropy_bits']:10.3f} "
                      f"{g['rank']:5d} {g['cond']:10.2f} {g['eff_rank_frac']:8.3f} "
                      f"{g['log_pseudo_det']:9.2f} {g['mean_gii']:8.5f} "
                      f"{g['max_dev_from_I4']:8.4f} {g['grad_var']:10.3e} "
                      f"{g.get('grad_share_bottom_decile', float('nan')):8.4f}")
            print(f"           ({time.time()-t0:.0f}s)")
            G.ck(TAG, key, done)
    summarise(done, pattern, L, tag_suffix)
    return done


def summarise(done, pattern, L, tag_suffix=""):
    done = {k: v for k, v in done.items() if "|s" in k}   # drop pre-seed-format cells
    print()
    print("=" * 112)
    print(f"B2. PAIRED COMPARISONS ACROSS {len(done)} TARGETS -- ansatz {pattern} L={L}")
    print("    The decisive row is COND vs scram_state: identical per-residue entropies,")
    print("    target information destroyed.  A null there refutes H2.")
    print("=" * 112)
    stats = ("cond", "eff_rank_frac", "log_pseudo_det", "max_dev_from_I4", "grad_var",
             "entropy_bits", "grad_share_bottom_decile")
    pairs = [("COND", "uniform"), ("COND", "random"), ("COND", "scram_state"),
             ("COND", "scram_res"), ("COND", "cross"), ("COND", "dirichlet"),
             ("scram_state", "uniform")]
    out = {}
    for st in stats:
        print(f"\n  --- {st} ---")
        print(f"  {'comparison':26s} {'mean diff':>11s} {'CI95':>24s} "
              f"{'W/L':>7s} {'med':>10s} {'verdict':>12s}")
        for a, b in pairs:
            xa = [done[p]["arms"][a]["geom"].get(st, np.nan) for p in done]
            xb = [done[p]["arms"][b]["geom"].get(st, np.nan) for p in done]
            if st == "cond":
                xa = [np.log10(v) if v > 0 else np.nan for v in xa]
                xb = [np.log10(v) if v > 0 else np.nan for v in xb]
            pr = G.paired(xa, xb)
            out[f"{st}|{a}-{b}"] = pr
            if pr.get("n", 0) == 0:
                continue
            print(f"  {a+' - '+b:26s} {pr['mean']:+11.4f} "
                  f"{f'[{pr['ci_lo']:+.4f}, {pr['ci_hi']:+.4f}]':>24s} "
                  f"{f'{pr['win']}/{pr['loss']}':>7s} {pr['median']:+10.4f} "
                  f"{pr['verdict']:>12s}")
    print("\n  (`cond` is compared on log10 so the ratio is the effect size.)")
    print("  W = first arm LOWER.  Fit quality is checked separately: if COND fits worse")
    print("  than its scramble the comparison is confounded and is reported as such.")
    kl_c = [done[p]["arms"]["COND"]["kl"] for p in done]
    kl_s = [done[p]["arms"]["scram_state"]["kl"] for p in done]
    prk = G.paired(kl_c, kl_s)
    out["FITQUALITY|COND-scram_state"] = prk
    print(f"\n  FIT-QUALITY CONTROL  KL(COND) - KL(scram_state) = {prk['mean']:+.5f} "
          f"[{prk['ci_lo']:+.5f}, {prk['ci_hi']:+.5f}]  W/L {prk['win']}/{prk['loss']} "
          f"-> {prk['verdict']}")
    G.ck(TAG, f"B2_paired_{pattern}L{L}{tag_suffix}", out)
    return out


if __name__ == "__main__":
    main()
