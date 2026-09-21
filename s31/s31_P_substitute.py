#!/usr/bin/env python
"""s31/s31_P_substitute.py -- S31 lane P: THE DECISIVE SUBSTITUTION.

Pre-registration: `s31/PREREG_S31_P.md`.  Read it first; it was written before any number
here existed and it names the falsifier in Angstrom on the built chain.

WHAT THIS DOES.  Lane L (S31-L4 / `s31/LIT_L.md` L1.1) proved the deployed CVaR free energy
`F(p) = CVaR_alpha(E;p) - T H(p)` (`core/quantum.py:993`) is CONVEX in `p` with a closed-form
global minimiser pinned by one scalar:

    s*   = argmax_s { s - T log sum_i exp( (s - E_i)_+ / (alpha T) ) }      1-D, concave
    p*_i propto exp( (s* - E_i)_+ / (alpha T) )                             hinged Gibbs

and that the shipped `run_cvar_vqe` is strictly worse in 12/12 synthetic cells by an
EXPRESSIVITY gap.  This script substitutes `p*` for the circuit's `p` and scores the
endpoint on tuning126, at the BUILT CHAIN, in six arms:

    A  quantum OFF                       stored canonical prod cloud       -> project
    B  ON, p_vqe, SELECTION readout       o[consensus_medoid(block, p_vqe)] -> project
    C  ON, p*,    SELECTION readout       o[consensus_medoid(block, p*)]    -> project
    D  ON, p_vqe, CONVEX readout          average_weighted(Wo, block, p_vqe)-> project
    E  ON, p*,    CONVEX readout          average_weighted(Wo, block, p*)   -> project
    F  ON, uniform, CONVEX readout        average_weighted(Wo, block, 1)    -> project

Arm F is the SHIPPED zero-information ablation (`core/pipeline.py:1116`), used as built.
Arms B/D reproduce the current quantum arm.  Arm A is the deployable baseline, 3.2105 A.

THE PROJECTION IS DETERMINISTIC BUT CHAOTIC (S31-L4, defect D-B): amplification ~1e13 from
input bits to branch choice.  THEREFORE ALL SIX ARMS OF A TARGET ARE PROJECTED IN THE SAME
PROCESS, IN ONE CALL, FROM THE SAME STORED POOL CACHE, and arm A is re-projected from the
stored canonical cloud rather than quoted.  No built-chain claim below 0.0107 A.

NATIVE-FREE.  `p*` is a function of `(E, alpha, T)` only.  `E = zrank(dis[o])` is the shipped
distogram Bayes-risk score of the pool windows under the HELD-OUT-FOLD model; `(alpha, T)`
come from the leave-fold-out table `core.pipeline.VQE_LFO`, used as shipped.  The ONLY oracle
read in this file is `nat_ca`, used to score RMSD after every structure is final, and it is
read in exactly one function, `_oracle_rmsd`.

    python s31/s31_P_substitute.py check
    python s31/s31_P_substitute.py run [--shard i --n-shards N] [--limit N] [--pdbs A,B]
    python s31/s31_P_substitute.py analyse
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import core                                          # noqa: E402
from core import pipeline as PL                      # noqa: E402
from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from s31 import s31_C_cache as CA                    # noqa: E402

RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "s31_P_rows.jsonl")
OUT = os.path.join(RESULTS, "s31_P_substitute.json")
STRUCTS = os.path.join(ROOT, "s29", "results", "s29_O_structs")
CANON_ROWS = os.path.join(ROOT, "s29", "results", "s29_O_chain_rows.jsonl")

#: the deployed quantum-stage settings, read from the shipped Config rather than retyped.
NQ = PL.PROD.vqe_qubits            # 7
LAYERS = PL.PROD.vqe_layers        # 3
ITERS = PL.PROD.vqe_iters          # 50
VSEED = PL.PROD.vqe_seed           # 0
DIM = 1 << NQ                      # 128
M75 = PL.PROD.m                    # 75
ARMS = ("A", "B", "C", "D", "E", "F")
#: S31-L4 defect D-B: the instrument's own reprojection noise floor on the built chain.
CHAIN_FLOOR = 0.0107
CANON_PROD_CHAIN = 3.2105


# ============================================================ the closed form (lane L L1.1)
def _phi(s, E, alpha, T):
    """The 1-D concave dual `s - T log sum_i exp((s - E_i)_+ / (alpha T))`."""
    z = np.maximum(s - E, 0.0) / (alpha * T)
    m = float(z.max())
    return float(s) - T * (m + math.log(float(np.exp(z - m).sum())))


def _p_of_s(s, E, alpha, T):
    z = np.maximum(s - E, 0.0) / (alpha * T)
    z = z - z.max()
    w = np.exp(z)
    return w / w.sum()


def closed_form(E, alpha, T, iters=400):
    """Golden-section maximisation of the concave dual; returns (p_star, s_star, phi_star).

    For alpha = 1 the dual is FLAT for s >= max(E) and `p_of_s` is constant there (it reduces
    to the plain Gibbs distribution exp(-E/(alpha T))), so landing anywhere in the plateau is
    the same answer -- asserted in `check`.
    """
    E = np.asarray(E, float)
    lo, hi = float(E.min()) - 1.0, float(E.max()) + 1.0
    g = (math.sqrt(5.0) - 1.0) / 2.0
    c, d = hi - g * (hi - lo), lo + g * (hi - lo)
    fc, fd = _phi(c, E, alpha, T), _phi(d, E, alpha, T)
    for _ in range(iters):
        if fc > fd:
            hi, d, fd = d, c, fc
            c = hi - g * (hi - lo)
            fc = _phi(c, E, alpha, T)
        else:
            lo, c, fc = c, d, fd
            d = lo + g * (hi - lo)
            fd = _phi(d, E, alpha, T)
        if hi - lo < 1e-15:
            break
    s = 0.5 * (lo + hi)
    return _p_of_s(s, E, alpha, T), float(s), _phi(s, E, alpha, T)


def free_energy_of(qm, E, p, alpha, T):
    """`F = CVaR_alpha(E;p) - T H(p)` using the SHIPPED `cvar_exact`.  Returns (F, CVaR, H, q)."""
    v, q, _dp = qm.cvar_exact(E, p, alpha)
    H = float(-(p * np.log(np.maximum(p, 1e-15))).sum())
    return float(v) - T * H, float(v), H, float(q)


def grad_F(qm, E, p, alpha, T):
    """dF/dp, exactly as `core.quantum.free_energy` forms it (`d = dp - T*dH`)."""
    _v, _q, dp = qm.cvar_exact(E, p, alpha)
    lp = np.log(np.maximum(p, 1e-300))
    return dp + T * (lp + 1.0)


def grad_F_at_s(E, p, alpha, T, s):
    """The same subgradient, but with the R-U variable held at `s` instead of at
    `cvar_exact`'s discrete quantile.

    WHY BOTH EXIST, and it is not a redundancy.  At the optimum `p*` the inner max over `s`
    is attained on a FLAT PLATEAU `[E_(k), E_(k+1)]` whenever the cumulative mass reaches
    `alpha` exactly at a data point -- which is precisely the stationarity condition.  Every
    `s` in that plateau is a valid maximiser and gives the SAME CVaR value but a DIFFERENT
    element of the subdifferential.  `cvar_exact` returns the element at the left endpoint;
    the dual returns the one at `s*`.  So `grad_F` is not constant at `p*` and that is a
    property of the kink, NOT a violation of optimality -- the rigorous certificate is the
    duality gap, which is one-sided and needs no differentiability at all.  This function
    checks the algebra: at `s*` the subgradient IS exactly constant, analytically.
    """
    dp = np.where(E < s, (E - s) / alpha, 0.0)
    return dp + T * (np.log(np.maximum(p, 1e-300)) + 1.0)


def mirror_descent(qm, E, alpha, T, iters=3000, eta0=0.3):
    """An INDEPENDENT optimiser over the full simplex -- entropic mirror descent.

    Its job is to try to beat `p*`.  If it ever does by more than 1e-6 the closed form is
    not the global optimum and the lane aborts.  Landing ABOVE `p*` is uninformative (it
    only means this optimiser did not converge), so the abort condition is one-sided.
    """
    N = len(E)
    p = np.full(N, 1.0 / N)
    best = free_energy_of(qm, E, p, alpha, T)[0]
    for t in range(1, iters + 1):
        g = grad_F(qm, E, p, alpha, T)
        l = np.log(np.maximum(p, 1e-300)) - (eta0 / math.sqrt(t)) * g
        l -= l.max()
        p = np.exp(l)
        p /= p.sum()
        f = free_energy_of(qm, E, p, alpha, T)[0]
        if f < best:
            best = f
    return float(best)


def tv(a, b):
    return float(0.5 * np.abs(np.asarray(a, float) - np.asarray(b, float)).sum())


def ess(w):
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / max(float((w * w).sum()), 1e-300))


def bits(H):
    return float(H) / math.log(2.0)


# ============================================================ the one oracle read
def _oracle_rmsd(X, nat):
    """ORACLE.  The ONLY native read in this file, and it runs after a structure is final."""
    return float(I.ca_rmsd(np.asarray(X, float), np.asarray(nat, float)))


# ============================================================ one target, all six arms
def seq_fold_map():
    return {t["pdb"]: (t["seq"], int(t["fold"]), int(t["n"])) for t in I.targets()}


_SFM = None


def target_row(pdb, verbose=True):
    global _SFM
    if _SFM is None:
        _SFM = seq_fold_map()
    t0 = time.time()
    seq, fold, n = _SFM[pdb]
    z = CA.load(pdb)
    W = np.asarray(z["W"], float)                       # NATIVE-FREE
    order = np.asarray(z["order"], int)                 # NATIVE-FREE (DIS score order)
    dis = np.asarray(z["dis"], float)                   # NATIVE-FREE (shipped distogram score)
    top75 = np.asarray(z["top75"], int)                 # NATIVE-FREE
    assert int(z["fold"]) == fold and int(z["n"]) == n
    nat = np.asarray(I.load_univ(pdb)["nat_ca"], float)   # ORACLE, scoring only

    # ---- the hypothesis set and the energy, exactly as `core/pipeline.py:860-866` ----
    o = order[:DIM]
    assert len(o) == DIM, "%s: pool smaller than 2^%d" % (pdb, NQ)
    E = PL._zrank(dis[o])
    alpha, T = PL.VQE_LFO[int(fold) % len(PL.VQE_LFO)]
    block = I.pairwise_rmsd(W[o])
    Wo = W[o]

    qm = core.backend("quantum")

    # ---- the two distributions ----
    tq0 = time.time()
    p_vqe, cvar_v, H_v, _circ = qm.run_cvar_vqe(E, alpha, T, n=NQ, layers=LAYERS,
                                                iters=ITERS, seed=VSEED)
    t_vqe = time.time() - tq0
    p_vqe = np.asarray(p_vqe, float)
    tq1 = time.time()
    p_star, s_star, phi_star = closed_form(E, alpha, T)
    t_star = time.time() - tq1

    F_star, cv_star, H_star, q_star = free_energy_of(qm, E, p_star, alpha, T)
    F_vqe, cv_vqe, H_vqe, _q_vqe = free_energy_of(qm, E, p_vqe, alpha, T)
    gs = grad_F_at_s(E, p_star, alpha, T, s_star)
    kkt_star = float(np.abs(gs - gs.mean()).max())       # KKT at s*: analytically exact
    gq = grad_F(qm, E, p_star, alpha, T)
    kkt_shipped = float(np.abs(gq - gq.mean()).max())    # at cvar_exact's plateau endpoint
    F_mirror = mirror_descent(qm, E, alpha, T)

    # ---- the three readouts ----
    loc_B = PL.consensus_medoid(block, p_vqe)
    loc_C = PL.consensus_medoid(block, p_star)
    clk = PL.Clock()
    clouds = {}
    f28 = os.path.join(STRUCTS, "%s.npz" % pdb)
    with np.load(f28, allow_pickle=True) as zs:
        clouds["A"] = np.asarray(zs["prod"], float)      # the STORED canonical cloud
    C_chk, _b = I.coordinate_average(W[top75])
    prod_dev = float(np.abs(C_chk - clouds["A"]).max())
    clouds["B"] = Wo[loc_B]
    clouds["C"] = Wo[loc_C]
    clouds["D"] = PL.average_weighted(Wo, block, p_vqe, clk)
    clouds["E"] = PL.average_weighted(Wo, block, p_star, clk)
    clouds["F"] = PL.average_weighted(Wo, block, np.ones(DIM), clk)

    # ---- ONE projection pass, same process, same stored clouds ----
    row = dict(pdb=pdb, n=n, fold=fold, alpha=float(alpha), T=float(T),
               prod_cloud_dev=prod_dev,
               s_star=s_star, phi_star=phi_star, q_star=q_star,
               F_star=F_star, F_vqe=F_vqe, F_mirror=F_mirror,
               dual_gap=float(F_star - phi_star), kkt_star_max=kkt_star,
               kkt_shipped_max=kkt_shipped, s_minus_q=float(s_star - q_star),
               F_gap_vqe_minus_star=float(F_vqe - F_star),
               mirror_minus_star=float(F_mirror - F_star),
               cvar_star=cv_star, cvar_vqe=cv_vqe,
               H_star_bits=bits(H_star), H_vqe_bits=bits(H_vqe),
               ess_star=ess(p_star), ess_vqe=ess(p_vqe),
               tv_star_vqe=tv(p_star, p_vqe),
               tv_star_unif=tv(p_star, np.full(DIM, 1.0 / DIM)),
               tv_vqe_unif=tv(p_vqe, np.full(DIM, 1.0 / DIM)),
               sel_local_B=int(loc_B), sel_local_C=int(loc_C),
               sel_B=int(o[loc_B]), sel_C=int(o[loc_C]),
               disagree=bool(int(loc_B) != int(loc_C)),
               secs_vqe=t_vqe, secs_star=t_star)
    for a in ARMS:
        X = clouds[a]
        row["cloud_%s" % a] = _oracle_rmsd(X, nat)       # ORACLE
        pr = I.project(X, seq, fold)
        row["chain_%s" % a] = _oracle_rmsd(pr["ca"], nat)  # ORACLE
    row["secs"] = time.time() - t0
    if verbose:
        print("  %s n=%2d f%d a=%.2f T=%.2f | Fgap %+.5f TV %.3f H* %.2f/%.2f b | "
              "medoid %s | chain A %.3f B %.3f C %.3f D %.3f E %.3f F %.3f (%.1fs)"
              % (pdb, n, fold, alpha, T, row["F_gap_vqe_minus_star"], row["tv_star_vqe"],
                 row["H_star_bits"], row["H_vqe_bits"],
                 "DIFF" if row["disagree"] else "same",
                 row["chain_A"], row["chain_B"], row["chain_C"],
                 row["chain_D"], row["chain_E"], row["chain_F"], row["secs"]), flush=True)
    return row


# ============================================================ rows io
A2_ROWS = os.path.join(RESULTS, "s31_P_a2_rows.jsonl")


def a2_row(pdb, verbose=True):
    """ARM A2 -- THE CODE-PATH CONTROL, and it is not optional.

    Arm A is `I.coordinate_average(W[top75])` (the stored canonical cloud).  Arms D/E/F are
    `core.pipeline.average_weighted(...)`.  With UNIFORM weights on the SAME top-75 those two
    are the same operator -- superpose on the medoid, take the mean -- but they are DIFFERENT
    IMPLEMENTATIONS (`s12.instrument.superpose_batch` vs `s8.consensus2.superpose_batch`) and
    they agree only to ~1e-14.  S31-L4 measured that a 1e-14 relative cloud perturbation moved
    one target's built chain by 0.511 A, amplification ~1e13.

    So `F - A` and `D - A` contain an unknown amount of pure implementation noise, and without
    this control I could not tell that from a result.  A2 is `average_weighted` with uniform
    weights on the top-75: same set as A, same weights as A, other implementation.  `A2 - A`
    IS the code-path noise on the built chain, measured rather than assumed, and `F - A2` is
    then the clean top-75 -> top-128 widening effect with both sides on one implementation.
    """
    global _SFM
    if _SFM is None:
        _SFM = seq_fold_map()
    t0 = time.time()
    seq, fold, n = _SFM[pdb]
    z = CA.load(pdb)
    W = np.asarray(z["W"], float)
    top = np.sort(np.asarray(z["top75"], int))
    nat = np.asarray(I.load_univ(pdb)["nat_ca"], float)   # ORACLE, scoring only
    P75 = I.pairwise_rmsd(W[top])
    C1, _b = I.coordinate_average(W[top], P75)
    C2 = PL.average_weighted(W[top], P75, np.ones(len(top)), PL.Clock())
    with np.load(os.path.join(STRUCTS, "%s.npz" % pdb), allow_pickle=True) as zs:
        C0 = np.asarray(zs["prod"], float)
    pr = I.project(C2, seq, fold)
    row = dict(pdb=pdb, n=n, fold=fold,
               cloud_dev_A2_vs_A=float(np.abs(C2 - C0).max()),
               cloud_dev_recomputed_vs_stored=float(np.abs(C1 - C0).max()),
               cloud_A2=_oracle_rmsd(C2, nat),                 # ORACLE
               chain_A2=_oracle_rmsd(pr["ca"], nat),           # ORACLE
               secs=time.time() - t0)
    if verbose:
        print("  %s A2 chain %.4f  cloud dev vs A %.3e (%.1fs)"
              % (pdb, row["chain_A2"], row["cloud_dev_A2_vs_A"], row["secs"]), flush=True)
    return row


def rows_path(shard=None):
    return ROWS if shard is None else ROWS.replace(".jsonl", ".s%d.jsonl" % int(shard))


def e_constancy(pdbs):
    """MEASURE the zero-information claim at the HAMILTONIAN, not only at the state.

    `E = _zrank(dis[order[:128]])`.  `rankdata` of 128 distinct values is 1..128 and
    standardising that is a constant vector, so `E` is the SAME diagonal on every target up
    to the pool's tie structure.  This reports the largest per-element deviation of any
    target's `E` from the tie-free reference `zrank(0..127)`, and how many targets carry a
    tie in their top-128 at all -- the only mechanism by which `E` can vary.

    It also reports how many targets' top-128 SET differs between the DIS order used here
    (production tie key) and `np.argsort(sc, kind="stable")` as `core/pipeline.py:757` forms
    it, because an exact tie at the 128 boundary can admit a different candidate.
    """
    ref = PL._zrank(np.arange(DIM, dtype=float))
    worst, n_tied, n_setdiff, setdiff = 0.0, 0, 0, []
    for p in pdbs:
        z = CA.load(p)
        dis = np.asarray(z["dis"], float)
        order = np.asarray(z["order"], int)
        o = order[:DIM]
        E = PL._zrank(dis[o])
        worst = max(worst, float(np.abs(np.sort(E) - ref).max()))
        if len(np.unique(dis[o])) < DIM:
            n_tied += 1
        st = np.argsort(dis, kind=PL.PROD.tie_break)[:DIM]
        if set(o.tolist()) != set(st.tolist()):
            n_setdiff += 1
            setdiff.append(p)
    return dict(max_abs_dev_from_tiefree_reference=worst, n_targets_with_ties_in_top128=n_tied,
                n_targets_top128_set_differs_from_stable_argsort=n_setdiff,
                top128_set_differs_pdbs=setdiff, reference_first5=ref[:5].tolist(),
                reference_last3=ref[-3:].tolist())


def read_a2():
    out = {}
    for f in sorted(glob.glob(A2_ROWS.replace(".jsonl", "*.jsonl"))):
        with open(f) as fh:
            for line in fh:
                if line.strip():
                    r = json.loads(line)
                    out[r["pdb"]] = r
    return out


def phase_a2(pdbs, shard=None, n_shards=1):
    done = read_a2()
    if shard is not None:
        pdbs = [p for k, p in enumerate(pdbs) if k % int(n_shards) == int(shard)]
    todo = [p for p in pdbs if p not in done]
    path = A2_ROWS if shard is None else A2_ROWS.replace(".jsonl", ".s%d.jsonl" % int(shard))
    print("a2 shard %s/%s: %d to do -> %s" % (shard, n_shards, len(todo),
                                              os.path.basename(path)), flush=True)
    for p in todo:
        append_row(path, a2_row(p))
    print("a2 rows:", path)


def read_rows():
    out = {}
    for f in sorted(glob.glob(ROWS.replace(".jsonl", "*.jsonl"))):
        if "_a2_" in f:
            continue
        with open(f) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                prev = out.get(r["pdb"])
                if prev is not None and abs(prev["chain_E"] - r["chain_E"]) > 1e-9:
                    raise RuntimeError("row disagreement for %s: chain_E %.9f vs %.9f"
                                       % (r["pdb"], prev["chain_E"], r["chain_E"]))
                out[r["pdb"]] = r
    return out


def append_row(path, row):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a") as fh:
        fh.write(json.dumps(row) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


# ============================================================ phases
def phase_check():
    """Plumbing + the closed form's own certificates, before any endpoint number."""
    print("backends:", core.backend_report())
    print("deployed quantum settings: n=%d layers=%d iters=%d seed=%d dim=%d m=%d"
          % (NQ, LAYERS, ITERS, VSEED, DIM, M75))
    print("VQE_LFO:", PL.VQE_LFO)
    qm = core.backend("quantum")

    # 1. the closed form against an independent optimiser, on lane L's synthetic ladder
    print("\n-- closed form vs mirror descent vs the shipped circuit (synthetic z-rank ladder) --")
    rng = np.random.default_rng(0)
    for n in (7, 9):
        for alpha in (0.1, 0.25, 1.0):
            for T in (0.1, 0.3):
                E = PL._zrank(rng.permutation(1 << n).astype(float))
                p_s, s_s, phi_s = closed_form(E, alpha, T)
                F_s = free_energy_of(qm, E, p_s, alpha, T)[0]
                F_m = mirror_descent(qm, E, alpha, T)
                p_v, _c, _h, _ = qm.run_cvar_vqe(E, alpha, T, n=n, layers=LAYERS,
                                                 iters=ITERS, seed=VSEED)
                F_v = free_energy_of(qm, E, np.asarray(p_v, float), alpha, T)[0]
                gs = grad_F_at_s(E, p_s, alpha, T, s_s)
                gq = grad_F(qm, E, p_s, alpha, T)
                print("  n=%d a=%.2f T=%.2f | F* %.9f dual %.2e kkt(s*) %.2e kkt(q) %.2e | "
                      "mirror %+.2e | vqe %+.6f | TV %.3f"
                      % (n, alpha, T, F_s, abs(F_s - phi_s),
                         float(np.abs(gs - gs.mean()).max()),
                         float(np.abs(gq - gq.mean()).max()), F_m - F_s, F_v - F_s,
                         tv(p_s, p_v)))
                assert abs(F_s - phi_s) < 1e-6, "strong duality failed"
                assert F_m - F_s > -1e-6, "mirror descent beat the closed form"
                assert float(np.abs(gs - gs.mean()).max()) < 1e-9, "KKT at s* not constant"
                if alpha >= 1.0:
                    # the alpha = 1 plateau: p* must be the plain Gibbs distribution
                    gb = np.exp(-(E - E.min()) / (alpha * T))
                    gb /= gb.sum()
                    assert tv(p_s, gb) < 1e-10, "alpha=1 closed form is not Gibbs"

    # 2. the cache reproduces the stored canonical clouds bit-for-bit
    print("\n-- cache -> canonical prod cloud --")
    dev = []
    for pdb in [t["pdb"] for t in I.targets()][:8]:
        z = CA.load(pdb)
        C, _b = I.coordinate_average(np.asarray(z["W"], float),
                                     None) if False else I.coordinate_average(
            np.asarray(z["W"], float)[np.asarray(z["top75"], int)])
        with np.load(os.path.join(STRUCTS, "%s.npz" % pdb), allow_pickle=True) as zs:
            dev.append(float(np.abs(C - np.asarray(zs["prod"], float)).max()))
    print("  max |recomputed - stored| over 8 targets: %.3e (0.0 = bit-identical)" % max(dev))
    assert max(dev) == 0.0, "cache does not reproduce the stored canonical clouds"

    # 3. the native-free trace of p*
    print("\n-- p* native-free trace --")
    print("  p* = f(E, alpha, T);  E = zrank(dis[order[:128]]);  dis = shipped distogram")
    print("  Bayes-risk score under the HELD-OUT-fold model;  (alpha,T) = VQE_LFO[fold].")
    print("  The only oracle read in this module is `_oracle_rmsd`, after structures are final.")
    src = open(os.path.abspath(__file__)).read()
    #: the needle is assembled at run time so that this assertion does not count ITSELF.
    needle = '["' + "nat" + '_ca"]'
    assert src.count(needle) == 1, "more than one native read site"
    assert ("oracle" + "_rr") not in src and ('z["' + 'rr"]') not in src, "an oracle rr read crept in"
    print("  asserted: exactly one `nat_ca` read site, no oracle `rr` read.")
    print("\ncheck OK")


def phase_run(pdbs, shard=None, n_shards=1):
    done = read_rows()
    if shard is not None:
        pdbs = [p for k, p in enumerate(pdbs) if k % int(n_shards) == int(shard)]
    todo = [p for p in pdbs if p not in done]
    path = rows_path(shard)
    print("run shard %s/%s: %d targets, %d done, %d to do -> %s"
          % (shard, n_shards, len(pdbs), len(pdbs) - len(todo), len(todo),
             os.path.basename(path)), flush=True)
    t0 = time.time()
    for k, pdb in enumerate(todo):
        append_row(path, target_row(pdb))
        if (k + 1) % 10 == 0 or k == len(todo) - 1:
            print("  [%d/%d] elapsed %.1f min" % (k + 1, len(todo), (time.time() - t0) / 60),
                  flush=True)
    print("rows:", path)


def _canon_prod():
    v = {}
    with open(CANON_ROWS) as fh:
        for line in fh:
            r = json.loads(line)
            if r["item"] == "prod":
                v[r["pdb"]] = float(r["rmsd_chain"])
    return v


def phase_analyse(write=True):
    rows = read_rows()
    pdbs = [t["pdb"] for t in I.targets()]
    have = [p for p in pdbs if p in rows]
    if len(have) != len(pdbs):
        print("WARNING: %d/%d targets present -- analysis is PARTIAL" % (len(have), len(pdbs)))
    pdbs = have
    folds = ST.pinned_folds(pdbs)
    G = lambda k: np.array([rows[p][k] for p in pdbs], float)          # noqa: E731

    chain = {a: G("chain_%s" % a) for a in ARMS}
    cloud = {a: G("cloud_%s" % a) for a in ARMS}

    # ---- certificates ----
    cert = dict(
        n=len(pdbs),
        max_dual_gap=float(np.abs(G("dual_gap")).max()),
        max_kkt_star=float(G("kkt_star_max").max()),
        max_kkt_shipped=float(G("kkt_shipped_max").max()),
        max_abs_s_minus_q=float(np.abs(G("s_minus_q")).max()),
        min_mirror_minus_star=float(G("mirror_minus_star").min()),
        n_vqe_worse=int((G("F_gap_vqe_minus_star") > 0).sum()),
        n_vqe_better=int((G("F_gap_vqe_minus_star") < 0).sum()),
        F_gap_mean=float(G("F_gap_vqe_minus_star").mean()),
        F_gap_min=float(G("F_gap_vqe_minus_star").min()),
        F_gap_max=float(G("F_gap_vqe_minus_star").max()),
        max_prod_cloud_dev=float(G("prod_cloud_dev").max()),
        tv_star_vqe_mean=float(G("tv_star_vqe").mean()),
        tv_star_vqe_min=float(G("tv_star_vqe").min()),
        tv_star_vqe_max=float(G("tv_star_vqe").max()),
        tv_star_unif_mean=float(G("tv_star_unif").mean()),
        tv_vqe_unif_mean=float(G("tv_vqe_unif").mean()),
        H_star_bits_mean=float(G("H_star_bits").mean()),
        H_vqe_bits_mean=float(G("H_vqe_bits").mean()),
        ess_star_mean=float(G("ess_star").mean()),
        ess_vqe_mean=float(G("ess_vqe").mean()),
        n_disagree=int(sum(bool(rows[p]["disagree"]) for p in pdbs)),
        disagree_pdbs=[p for p in pdbs if rows[p]["disagree"]],
        secs_vqe_mean=float(G("secs_vqe").mean()),
        secs_star_mean=float(G("secs_star").mean()),
    )
    #: per-alpha breakdown.  Lane L's L1.1 verified the closed form at T = 0.1 and T = 0.05;
    #: the DEPLOYED T is 0.3 for every fold (`core/pipeline.py:113`), so none of its 12 cells
    #: is at the deployed temperature.  This block is here because the entropy ordering is
    #: NOT the same at T = 0.3 as in L1.1's table, and the difference is alpha-dependent.
    #: TARGET-INDEPENDENCE.  `E = zrank(dis[o])` is the standardised rank of 128 values, so up
    #: to the pool's tie structure it is THE SAME VECTOR on every target.  `p*` is a function
    #: of `(E, alpha, T)` alone and `run_cvar_vqe` is seeded at 0, so BOTH distributions are
    #: (to within the ties) one fixed vector per `alpha`.  Measured, not asserted: the
    #: within-alpha spread of every scalar summary of them.
    cert["E_constancy"] = e_constancy(pdbs)
    al = G("alpha")
    cert["target_independence"] = {}
    for a_val in sorted(set(al.tolist())):
        m = al == a_val
        blk = {}
        for k in ("H_star_bits", "H_vqe_bits", "ess_star", "ess_vqe", "F_star", "F_vqe",
                  "tv_star_vqe"):
            v = G(k)[m]
            blk[k] = dict(mean=float(v.mean()), sd=float(v.std(ddof=1)),
                          range=float(v.max() - v.min()))
        cert["target_independence"]["%.2f" % a_val] = blk

    cert["by_alpha"] = {}
    for a_val in sorted(set(al.tolist())):
        m = al == a_val
        cert["by_alpha"]["%.2f" % a_val] = dict(
            n=int(m.sum()),
            folds=sorted(set(G("fold")[m].astype(int).tolist())),
            H_star_bits=float(G("H_star_bits")[m].mean()),
            H_vqe_bits=float(G("H_vqe_bits")[m].mean()),
            ess_star=float(G("ess_star")[m].mean()), ess_vqe=float(G("ess_vqe")[m].mean()),
            tv_star_vqe=float(G("tv_star_vqe")[m].mean()),
            tv_star_unif=float(G("tv_star_unif")[m].mean()),
            tv_vqe_unif=float(G("tv_vqe_unif")[m].mean()),
            F_gap=float(G("F_gap_vqe_minus_star")[m].mean()),
            n_disagree=int(sum(bool(rows[p]["disagree"]) for p, k in zip(pdbs, m) if k)),
            chain_A=float(chain["A"][m].mean()), chain_D=float(chain["D"][m].mean()),
            chain_E=float(chain["E"][m].mean()), chain_F=float(chain["F"][m].mean()))

    # arm A reproduction of the canonical endpoint
    can = _canon_prod()
    a_can = np.array([can[p] for p in pdbs], float)
    cert["armA_mean"] = float(chain["A"].mean())
    cert["canonical_prod_mean_same_targets"] = float(a_can.mean())
    cert["armA_minus_canonical_mean"] = float(chain["A"].mean() - a_can.mean())
    cert["armA_max_abs_dev_per_target"] = float(np.abs(chain["A"] - a_can).max())
    cert["armA_n_moved_over_floor"] = int((np.abs(chain["A"] - a_can) > CHAIN_FLOOR).sum())

    # ---- the endpoint table ----
    table = {}
    for a in ARMS:
        table[a] = dict(chain_mean=float(chain[a].mean()),
                        chain_se=float(chain[a].std(ddof=1) / math.sqrt(len(pdbs))),
                        chain_median=float(np.median(chain[a])),
                        cloud_mean=float(cloud[a].mean()))

    # ---- comparisons (family registered in s31/MULTIPLICITY.md, k = 9) ----
    CMP = [("P1_E_minus_D", "E", "D", "PRIMARY  convex readout: p* vs VQE p (BUILT CHAIN)"),
           ("P2_C_minus_B", "C", "B", "PRIMARY  selection readout: p* vs VQE p (BUILT CHAIN)"),
           ("S1_B_minus_A", "B", "A", "VQE selection vs quantum-OFF production (BUILT CHAIN)"),
           ("S2_C_minus_A", "C", "A", "p* selection vs quantum-OFF production (BUILT CHAIN)"),
           ("S3_D_minus_A", "D", "A", "VQE convex vs quantum-OFF production (BUILT CHAIN)"),
           ("S4_E_minus_A", "E", "A", "p* convex vs quantum-OFF production (BUILT CHAIN)"),
           ("S5_F_minus_A", "F", "A", "uniform convex vs quantum-OFF production (BUILT CHAIN)"),
           ("S6_E_minus_F", "E", "F", "p* convex vs the SHIPPED uniform ablation (BUILT CHAIN)"),
           ("S7_D_minus_F", "D", "F", "VQE convex vs the SHIPPED uniform ablation (BUILT CHAIN)")]
    def spread(d):
        """The per-target |delta| distribution -- AMENDMENT 1's reporting obligation.

        Lane A's structure (per-target 0.102 A cancelling to -0.0044 A in the mean) is why a
        paired mean alone is not a sufficient report of this experiment: a null mean with a
        large spread means the substitution CHANGES THE ANSWER AND THE CHANGES CANCEL, which
        is a different sentence from `p*` makes no difference.
        """
        ad = np.abs(np.asarray(d, float))
        return dict(abs_mean=float(ad.mean()), abs_median=float(np.median(ad)),
                    abs_p90=float(np.percentile(ad, 90)), abs_max=float(ad.max()),
                    n_zero=int((ad == 0.0).sum()),
                    n_above_chain_floor=int((ad > CHAIN_FLOOR).sum()),
                    rms=float(np.sqrt((ad ** 2).mean())))

    #: ARM A2, the code-path control.  Present only after `phase_a2` has run; if it is
    #: missing the lane says so rather than silently dropping the control.
    a2 = read_a2()
    have_a2 = all(p in a2 for p in pdbs)
    if have_a2:
        chain["A2"] = np.array([a2[p]["chain_A2"] for p in pdbs], float)
        cloud["A2"] = np.array([a2[p]["cloud_A2"] for p in pdbs], float)
        table["A2"] = dict(chain_mean=float(chain["A2"].mean()),
                           chain_se=float(chain["A2"].std(ddof=1) / math.sqrt(len(pdbs))),
                           chain_median=float(np.median(chain["A2"])),
                           cloud_mean=float(cloud["A2"].mean()))
        cert["a2_max_cloud_dev_vs_A"] = float(max(a2[p]["cloud_dev_A2_vs_A"] for p in pdbs))
        CMP = CMP + [("X1_A2_minus_A", "A2", "A",
                      "CODE-PATH CONTROL: average_weighted(uniform, top-75) vs "
                      "coordinate_average(top-75) -- the SAME operator, two implementations "
                      "agreeing to ~1e-14, through a projection that amplifies ~1e13 "
                      "(BUILT CHAIN)"),
                     ("X2_F_minus_A2", "F", "A2",
                      "top-75 -> top-128 widening, both sides on the SAME implementation "
                      "(BUILT CHAIN)")]
    else:
        cert["a2_present"] = False

    comps = {}
    for key, a, b, lab in CMP:
        c = ST.compare(chain[a], chain[b], folds, names=pdbs, label=lab,
                       seed_parts=("s31P", key))
        c["below_chain_floor"] = bool(abs(c["effect"]) < CHAIN_FLOOR)
        c["per_target_spread"] = spread(chain[a] - chain[b])
        comps[key] = c
    # the cloud basis, carried as a secondary basis and labelled as one
    comps_cloud = {}
    for key, a, b, lab in CMP:
        cc = ST.compare(cloud[a], cloud[b], folds, names=pdbs,
                        label=lab.replace("BUILT CHAIN", "CA CLOUD"),
                        seed_parts=("s31P", key, "cloud"))
        cc["per_target_spread"] = spread(cloud[a] - cloud[b])
        comps_cloud[key] = cc

    out = dict(provenance=ST.provenance(__file__), n=len(pdbs), pdbs=pdbs,
               settings=dict(n_qubits=NQ, layers=LAYERS, iters=ITERS, seed=VSEED, dim=DIM,
                             m=M75, VQE_LFO={str(k): list(v) for k, v in PL.VQE_LFO.items()},
                             backends=core.backend_report()),
               certificates=cert, endpoint=table, comparisons=comps,
               comparisons_cloud_basis=comps_cloud,
               canonical_prod_chain=CANON_PROD_CHAIN, chain_floor=CHAIN_FLOOR)
    if write:
        ST.save_atomic(OUT, out, rows=None, n_expected=None, module_file=__file__)
    render(out)
    return out


def render(out):
    c, tb = out["certificates"], out["endpoint"]
    n = out["n"]
    L = []
    L.append("== S31 lane P -- the closed-form substitution, n = %d, BUILT CHAIN ==" % n)
    L.append("")
    L.append("CERTIFICATES (the closed form IS the global optimum of the deployed objective)")
    L.append("  strong duality   max |F(p*) - phi(s*)|        %.3e" % c["max_dual_gap"])
    L.append("  KKT at s*        max |grad F(p*) - mean|      %.3e  (analytic identity)"
             % c["max_kkt_star"])
    L.append("  KKT at cvar_exact's quantile                  %.3e  (the plateau, |s*-q| <= %.2e;"
             % (c["max_kkt_shipped"], c["max_abs_s_minus_q"]))
    L.append("                   a different valid subgradient element at the kink, NOT a defect)")
    L.append("  mirror descent   min (F_mirror - F*)          %+.3e  (>= 0 means it never beat p*)"
             % c["min_mirror_minus_star"])
    L.append("  circuit worse    %d / %d targets, F gap mean %+.5f  [%.5f, %.5f]"
             % (c["n_vqe_worse"], n, c["F_gap_mean"], c["F_gap_min"], c["F_gap_max"]))
    L.append("  TV(p*, p_vqe)    mean %.3f  [%.3f, %.3f]      H  p* %.2f bits vs vqe %.2f bits"
             % (c["tv_star_vqe_mean"], c["tv_star_vqe_min"], c["tv_star_vqe_max"],
                c["H_star_bits_mean"], c["H_vqe_bits_mean"]))
    L.append("  ESS              p* %.1f vs vqe %.1f of %d     TV to uniform: p* %.3f vqe %.3f"
             % (c["ess_star_mean"], c["ess_vqe_mean"], DIM,
                c["tv_star_unif_mean"], c["tv_vqe_unif_mean"]))
    L.append("  cost             p* %.4f s vs run_cvar_vqe %.4f s per target"
             % (c["secs_star_mean"], c["secs_vqe_mean"]))
    L.append("")
    L.append("ARM A reproduces the canonical endpoint")
    L.append("  arm A %.4f vs canonical s29 prod %.4f on the same targets (delta %+.4f);"
             % (c["armA_mean"], c["canonical_prod_mean_same_targets"],
                c["armA_minus_canonical_mean"]))
    L.append("  worst per-target reprojection deviation %.4f A, %d targets above the %.4f A floor"
             % (c["armA_max_abs_dev_per_target"], c["armA_n_moved_over_floor"], CHAIN_FLOOR))
    L.append("")
    ec = c["E_constancy"]
    L.append("THE HAMILTONIAN DIAGONAL IS A TARGET-INDEPENDENT CONSTANT")
    L.append("  E = zrank(dis[order[:128]]); rankdata of 128 distinct values is 1..128, so E is")
    L.append("  the fixed vector %s ... %s on EVERY target."
             % (np.round(ec["reference_first5"], 6).tolist(),
                np.round(ec["reference_last3"], 6).tolist()))
    L.append("  max |sort(E_target) - tie-free reference| over %d targets: %.3e"
             % (n, ec["max_abs_dev_from_tiefree_reference"]))
    L.append("  targets with any tie inside the top-128: %d/%d   top-128 SET differs from"
             % (ec["n_targets_with_ties_in_top128"], n))
    L.append("  np.argsort(sc, kind='%s') on %d/%d targets %s"
             % (PL.PROD.tie_break, ec["n_targets_top128_set_differs_from_stable_argsort"], n,
                ec["top128_set_differs_pdbs"]))
    L.append("")
    L.append("TARGET-INDEPENDENCE OF BOTH DISTRIBUTIONS (E = zrank of 128 values is the SAME")
    L.append("vector on every target up to ties, so p* and p_vqe are ONE FIXED RANK-WEIGHTING")
    L.append("per alpha -- they carry no target-specific information at all)")
    for a_val, v in sorted(c["target_independence"].items()):
        L.append("  alpha=%s  H*  %.4f bits (sd %.1e, range %.1e)   ESS* %.2f of %d (sd %.1e)"
                 % (a_val, v["H_star_bits"]["mean"], v["H_star_bits"]["sd"],
                    v["H_star_bits"]["range"], v["ess_star"]["mean"], DIM, v["ess_star"]["sd"]))
        L.append("            H_vqe %.4f bits (sd %.1e, range %.1e)  ESS_vqe %.2f (sd %.1e)"
                 % (v["H_vqe_bits"]["mean"], v["H_vqe_bits"]["sd"], v["H_vqe_bits"]["range"],
                    v["ess_vqe"]["mean"], v["ess_vqe"]["sd"]))
    L.append("")
    L.append("BY ALPHA (the DEPLOYED T is 0.3 on every fold; L1.1's 12 cells were at T=0.1/0.05,")
    L.append("so none of them is at the deployed temperature -- and the entropy ORDERING differs)")
    for a_val, v in sorted(c["by_alpha"].items()):
        L.append("  alpha=%s folds %s n=%d | H* %.2f vs H_vqe %.2f bits | ESS %.1f vs %.1f | "
                 "TV(p*,vqe) %.3f | TV to unif: p* %.3f vqe %.3f | Fgap %+.4f | medoid DIFF %d/%d"
                 % (a_val, v["folds"], v["n"], v["H_star_bits"], v["H_vqe_bits"],
                    v["ess_star"], v["ess_vqe"], v["tv_star_vqe"], v["tv_star_unif"],
                    v["tv_vqe_unif"], v["F_gap"], v["n_disagree"], v["n"]))
        L.append("            chain A %.4f  D %.4f  E %.4f  F %.4f"
                 % (v["chain_A"], v["chain_D"], v["chain_E"], v["chain_F"]))
    L.append("")
    L.append("THE DISAGREEMENT COUNT (lane L's caveat)")
    L.append("  consensus_medoid(block, p*) != consensus_medoid(block, p_vqe) on %d / %d targets"
             % (c["n_disagree"], n))
    L.append("")
    L.append("ENDPOINT, BUILT CHAIN (mean +- SE; CA cloud in brackets)")
    names = {"A": "A  quantum OFF (production)        ",
             "B": "B  ON, VQE p,   SELECTION readout  ",
             "C": "C  ON, p*,      SELECTION readout  ",
             "D": "D  ON, VQE p,   CONVEX readout     ",
             "E": "E  ON, p*,      CONVEX readout     ",
             "F": "F  ON, uniform, CONVEX readout     ",
             "A2": "A2 code-path control (unif top-75)"}
    for a in (ARMS + ("A2",) if "A2" in tb else ARMS):
        L.append("  %s %.4f +- %.4f   median %.4f   [cloud %.4f]"
                 % (names[a], tb[a]["chain_mean"], tb[a]["chain_se"], tb[a]["chain_median"],
                    tb[a]["cloud_mean"]))
    L.append("")
    L.append("COMPARISONS (paired n = %d, fold-clustered, MDE = 2.8016 x SE)" % n)
    for k, v in out["comparisons"].items():
        gate = ("NULL" if abs(v["effect_over_mde"]) < 0.7 else
                "NOT MEASURED" if abs(v["effect_over_mde"]) < 1.0 else "MEASURED")
        if v["below_chain_floor"]:
            gate += " / BELOW CHAIN FLOOR"
        ci = v.get("ci95_fold") or [float("nan")] * 2
        L.append("  %-14s %+.4f A  SE %.4f  MDE %.4f  %.2fx  %dW/%dL/%dT  fold CI [%+.4f, %+.4f]  %s"
                 % (k, v["effect"], v["se"], v["mde"], abs(v["effect_over_mde"]),
                    v["n_better"], v["n_worse"], v["n_tied"], ci[0], ci[1], gate))
        L.append("                 median %+.4f   %d/%d folds same sign   verdict: %s"
                 % (v["median_effect"], v.get("n_folds", 0) and v["folds_same_sign"],
                    v.get("n_folds", 0), v["verdict"]))
    L.append("")
    L.append("PER-TARGET |DELTA| DISTRIBUTION (AMENDMENT 1: a null mean with a wide spread is")
    L.append("NOT 'no effect' -- it is 'the answer changes on most targets and cancels')")
    for k, v in out["comparisons"].items():
        s = v["per_target_spread"]
        L.append("  %-14s |d| mean %.4f  median %.4f  p90 %.4f  max %.4f  rms %.4f  "
                 "exactly 0 on %d/%d  above %.4f A floor on %d/%d"
                 % (k, s["abs_mean"], s["abs_median"], s["abs_p90"], s["abs_max"], s["rms"],
                    s["n_zero"], n, CHAIN_FLOOR, s["n_above_chain_floor"], n))
    L.append("")
    L.append("SAME, CA CLOUD basis (a DIFFERENT object from the built chain)")
    for k, v in out["comparisons_cloud_basis"].items():
        s = v["per_target_spread"]
        gate = ("NULL" if abs(v["effect_over_mde"]) < 0.7 else
                "NOT MEASURED" if abs(v["effect_over_mde"]) < 1.0 else "MEASURED")
        L.append("  %-14s %+.4f A  SE %.4f  %.2fx  %s   |d| mean %.4f  p90 %.4f  max %.4f"
                 % (k, v["effect"], v["se"], abs(v["effect_over_mde"]), gate,
                    s["abs_mean"], s["abs_p90"], s["abs_max"]))
    print("\n".join(L))
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["check", "run", "a2", "analyse"])
    ap.add_argument("--shard", type=int, default=None)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--pdbs", type=str, default="")
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    if a.phase == "check":
        phase_check()
        return
    if a.phase == "analyse":
        phase_analyse()
        return
    pdbs = [t["pdb"] for t in I.targets()]
    if a.pdbs:
        want = [s.strip() for s in a.pdbs.split(",") if s.strip()]
        pdbs = [p for p in pdbs if p in want]
    if a.limit:
        pdbs = pdbs[:a.limit]
    if a.phase == "a2":
        phase_a2(pdbs, a.shard, a.n_shards)
    else:
        phase_run(pdbs, a.shard, a.n_shards)


if __name__ == "__main__":
    main()
