#!/usr/bin/env python
"""s29/s29_B_tta.py -- S29 LANE B, MEASUREMENT 5: "TAIL THEN AGGREGATE".

Pre-registered in `s29/PREREG_S29_B.md` ADDENDUM 2 (read it first).  Derived by lane T in
S29-L15 (Q1: the deployed CVaR is exactly constant on 437 of 511 simplex directions, the entropy
term fills them uniformly, and the whole quantum stage reduces at the endpoint to one number m)
and S29-L17 section 4(b).

THE OBJECT
    F(theta) = CVaR_alpha(E; p) - T H(p) + lam f(R_alpha(p))
    R_alpha(p) = (1/alpha) [ sum_{x in S(p)} p_x W_x + (alpha - mass(S)) W_{x_q} ]

R_alpha is the CVaR tail's own coordinate average AS A FUNCTION OF p, so the objective sees which
candidates populate the tail and with what weight -- the freedom S29-L15 proves the deployed
objective is indifferent to.  By the envelope argument dR/dp_y = (W_y - W_{x_q})/alpha on the
strict tail and 0 above the VaR, so

    dF/dp_y = (E_y - q)/alpha - T dH/dp_y + lam <grad f(R), W_y - W_{x_q}>/alpha   (strict tail)
            =            0    - T dH/dp_y                                          (above the VaR)

and the parameter-shift chain is the deployed one: g = (PR[0::2] - PR[1::2]) @ dF/dp / 2.

THE FLATNESS MEASUREMENT (prereg B2.3, gate G5), in the basis lane T counts in: delta_y = e_y -
e_{y0} with y0 a fixed state ABOVE the VaR, y != y0, so there are D - 1 = 511 directions and the
D - m - 1 of them with y also above the VaR are the ones the CVaR term cannot see.  Reported per
term, for the EMITTED STRUCTURE as well as for the objective, and as the OVERLAP, which is the
quantity the gate reads:

    flat_obj(term)   fraction of the 511 directions on which that term is EXACTLY constant
    flat_readout     the same for the emitted structure C(p)
    overlap          among the directions on which C moves, the fraction on which the
                     INFORMATION-BEARING part of the objective (CVaR + lam f, entropy excluded
                     because it carries no information) is exactly constant

NATIVE-FREE / ORACLE SPLIT: nothing in the objective, the readout or the flatness measurement
reads a native.  Natives are read only in `oracle_rmsd`.

USAGE
    python s29/s29_B_tta.py --flat [--limit N]      measurement 5a, the mechanism (gate G5)
    python s29/s29_B_tta.py --analyse-flat
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import quantum as Q                       # noqa: E402
from s12 import instrument as I                     # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
FLAT_ROWS = os.path.join(RESULTS, "s29_B_tta_flat_rows.jsonl")
FLAT_OUT = os.path.join(RESULTS, "s29_B_tta_flat.json")
SUBSET_ROWS = os.path.join(RESULTS, "s29_B_tta_subset_rows.jsonl")
END_ROWS = os.path.join(RESULTS, "s29_B_tta_end_rows.jsonl")
CHAIN_ROWS = os.path.join(RESULTS, "s29_B_tta_chain_rows.jsonl")
N_UNTRAINED = 16
FAIL18 = {"1ID6", "1JBF", "1LB7", "2BFI", "2BP4", "2JN5", "2MQ2", "2N5C", "2NB7",
          "2NDM", "3BTB", "3SGO", "5W52", "7JS6", "7LCW", "8T63", "9KAR", "9L1M"}

ALPHA, TEMP, LAYERS, ITERS, LR, M_PROD = 0.18, 0.5, 3, 80, 0.15, 75
LAM_GRID = (0.0, 0.1, 0.3, 1.0, 3.0)                # prereg B2.6, fixed
ZERO_TOL = 0.0                                      # "exactly constant" means exactly


# ============================================================================ the tail operator
def tail_lambda(E: np.ndarray, p: np.ndarray, alpha: float):
    """The CVaR tail's weights lambda_x(p), its strict support and its boundary state.

    lambda_x = p_x / alpha on the strict tail (E_x < q), (alpha - mass)/alpha on the boundary
    state x_q, 0 above.  sum lambda = 1 by construction.  Same ordering convention as
    `core.quantum.cvar_exact` (stable sort on E), so the tail here IS the deployed tail.
    """
    E, p = np.asarray(E, float), np.asarray(p, float)
    o = np.argsort(E, kind="stable")
    c = np.cumsum(p[o])
    k = int(np.searchsorted(c, alpha, side="left"))
    k = min(k, len(E) - 1)
    q = float(E[o[k]])
    strict = o[:k]
    mass = float(c[k - 1]) if k > 0 else 0.0
    lam = np.zeros_like(p)
    lam[strict] = p[strict] / alpha
    x_q = int(o[k])
    lam[x_q] += (alpha - mass) / alpha
    return lam, strict, x_q, q, mass


def r_alpha(E, p, alpha, Wf):
    """R_alpha(p) flattened (3 n_res,). `Wf` is (D, 3 n_res) with ZERO rows on the padding."""
    lam, strict, x_q, q, mass = tail_lambda(E, p, alpha)
    return lam @ np.asarray(Wf, float), lam, strict, x_q, q


def dR_dp_rows(Wf, strict, x_q, alpha):
    """(len(strict), 3 n_res): dR/dp_y for y on the strict tail. Zero for every other y."""
    Wf = np.asarray(Wf, float)
    return (Wf[strict] - Wf[x_q][None, :]) / alpha


# ================================================================================ the objective
def dF_dp(E, p, alpha, T, lam_f, Wf, sur, n_res):
    """dF/dp for F = CVaR - T H + lam f(R_alpha), returned per TERM so flatness is auditable."""
    v, q, d_cvar = Q.cvar_exact(E, p, alpha)
    lp = np.log(np.maximum(p, 1e-15))
    H = float(-(p * lp).sum())
    d_ent = -T * (-(lp + 1.0))                       # the gradient of (-T H)
    R, lam_w, strict, x_q, _ = r_alpha(E, p, alpha, Wf)
    d_f = np.zeros_like(p)
    s_val = float("nan")
    if sur is not None:
        s_val, dS_dC = sur.value_grad(R.reshape(n_res, 3))
        gC = dS_dC.ravel()
        d_f[strict] = dR_dp_rows(Wf, strict, x_q, alpha) @ gC
    return dict(cvar=d_cvar, entropy=d_ent, f=d_f,
                total=d_cvar + d_ent + lam_f * d_f,
                info_bearing=d_cvar + lam_f * d_f,
                value=float(v - T * H + (lam_f * s_val if np.isfinite(s_val) else 0.0)),
                cvar_value=float(v), H=H, s_val=s_val, R=R, lam_w=lam_w,
                strict=strict, x_q=x_q)


def objective_theta(circ, theta, E, alpha, T, lam_f, Wf, sur, n_res):
    """F and its EXACT gradient by the deployed parameter-shift chain."""
    p = circ.probs(theta)
    d = dF_dp(E, p, alpha, T, lam_f, Wf, sur, n_res)
    PR = circ.probs_batch(circ._shift_grid(np.asarray(theta, float), np.pi / 2))
    g = (PR[0::2] - PR[1::2]) @ d["total"] / 2.0
    return d["value"], g, p, d


def run_tta_vqe(E, alpha, T, lam_f, Wf, sur, n_res, n, layers=LAYERS, iters=ITERS,
                seed=0, lr=LR):
    """`core.quantum.run_cvar_vqe` with the tail-average term. At lam = 0 the loop is the
    deployed one line for line (same RNG draw, same Adam, same order of operations)."""
    circ = Q.StatevectorCircuit(n, layers)
    rng = np.random.default_rng(seed)
    th = rng.normal(0.0, 0.6, circ.n_params())
    m = np.zeros_like(th)
    v = np.zeros_like(th)
    for t in range(1, iters + 1):
        if lam_f == 0.0:
            _, g, _, _, _ = Q.free_energy(circ, th, E, alpha, T)
        else:
            _, g, _, _ = objective_theta(circ, th, E, alpha, T, lam_f, Wf, sur, n_res)
        m = 0.9 * m + 0.1 * g
        v = 0.999 * v + 0.001 * g * g
        th = th - lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-8)
    f, g, p, d = objective_theta(circ, th, E, alpha, T, max(lam_f, 0.0), Wf, sur, n_res)
    return p, th, circ, d


# ============================================================== the flatness measurement (G5)
def flat_report(E, p, alpha, T, lam_f, Wf, sur, n_res) -> Dict:
    """Gate G5. Directions delta_y = e_y - e_{y0}, y0 the HIGHEST-E state (above the VaR), which
    is lane T's basis and gives D - 1 = 511 directions of which D - m - 1 sit above the VaR."""
    d = dF_dp(E, p, alpha, T, lam_f, Wf, sur, n_res)
    D = len(p)
    y0 = int(np.argmax(np.asarray(E, float)))
    ys = np.array([y for y in range(D) if y != y0])
    out = dict(D=D, n_directions=int(len(ys)), y0=y0,
               m_strict=int(len(d["strict"])), x_q=int(d["x_q"]))
    for term in ("cvar", "entropy", "f", "total", "info_bearing"):
        dd = d[term]
        out[f"flat_{term}"] = float(np.mean(np.abs(dd[ys] - dd[y0]) <= ZERO_TOL))
        out[f"nflat_{term}"] = int(np.sum(np.abs(dd[ys] - dd[y0]) <= ZERO_TOL))
    # --- the emitted structure's own flat set, for BOTH readouts
    # TTA readout: dC/dp_y = (W_y - W_xq)/alpha on the strict tail, 0 above the VaR
    dC = np.zeros(D)
    rows = dR_dp_rows(Wf, d["strict"], d["x_q"], alpha)
    dC[d["strict"]] = np.linalg.norm(rows, axis=1)
    out["flat_readout_tta"] = float(np.mean(np.abs(dC[ys] - dC[y0]) <= ZERO_TOL))
    out["nmove_readout_tta"] = int(np.sum(np.abs(dC[ys] - dC[y0]) > ZERO_TOL))
    # deployed readout: a uniform average over the tail SET, piecewise constant in p
    out["flat_readout_deployed"] = 1.0
    out["nmove_readout_deployed"] = 0
    # --- the OVERLAP, which is what gate G5 reads
    moving = ys[np.abs(dC[ys] - dC[y0]) > ZERO_TOL]
    dd = d["info_bearing"]
    out["overlap_tta"] = (float(np.mean(np.abs(dd[moving] - dd[y0]) <= ZERO_TOL))
                          if len(moving) else float("nan"))
    out["overlap_deployed"] = float("nan")            # 0 moving directions: undefined, by design
    out["dC_norm_tail_mean"] = float(np.mean(np.linalg.norm(rows, axis=1))) if len(rows) else 0.0
    out["s_val"] = d["s_val"]
    out["cvar_value"] = d["cvar_value"]
    out["H_nats"] = d["H"]
    return out


# ===================================================================== per target, native-free
def flat_target(pdb: str) -> List[Dict]:
    from s22 import qcand_lib as QC
    from s27 import run_pool as RP
    from s27 import s28_A_amp as AMP
    cand, ch, _ = RP.channels_for(pdb)
    E_real = RP.zr(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    top = np.sort(RP.topm(E_real, M_PROD, key)).astype(int)
    frame = AMP.Frame(cand.W, top)
    dg = I.distogram(pdb, cand.seq, cand.fold)
    sur = AMP.Surrogate(dg, int(cand.n))
    enc = QC.Encoding(E_real)
    E, dim, k = enc.E, enc.dim, cand.k
    Wf = np.zeros((dim, 3 * int(cand.n)))
    Wf[:k] = frame.Wf                                  # padding states carry NO structure
    rows = []
    # the production point: the DEPLOYED optimum (lam = 0), seed 0, exactly as shipped
    p0, th0, circ, _ = run_tta_vqe(E, ALPHA, TEMP, 0.0, Wf, sur, int(cand.n), enc.n_qubits, seed=0)
    lam_w0, strict0, x_q0, _, _ = tail_lambda(E, p0, ALPHA)
    tail_idx = np.sort(np.union1d(strict0, [x_q0])).astype(int)
    tail_idx = tail_idx[tail_idx < k]                  # real candidates only (padding never in tail)
    # the DEPLOYED readout operator, in the retained set's OWN medoid frame, exactly as shipped
    C_uni = I.coordinate_average(cand.W[tail_idx])[0]
    C_prod = I.coordinate_average(cand.W[top])[0]
    # frame.Wf[top].mean(0) is the SAME object: lane A's Frame superposes every window onto the
    # medoid of `top`, which is the frame `I.coordinate_average` averages in, so the two arrays
    # are bit-identical.  Asserted on the raw arrays, which is strictly stronger than an RMSD
    # check: `I.ca_rmsd` has a Kabsch/SVD floor of ~1.3e-7 A on IDENTICAL input (measured here
    # on 1A13: max|diff| exactly 0.0, ca_rmsd 1.274e-07), so any "identical to 1e-8" claim made
    # through `ca_rmsd` cannot be met and the first version of this gate mis-fired on it.
    assert np.max(np.abs(C_prod - frame.Wf[top].mean(0).reshape(int(cand.n), 3))) == 0.0
    R0 = (lam_w0 @ Wf).reshape(int(cand.n), 3)
    pad_mass = float(p0[k:].sum())
    for lam_f in LAM_GRID:
        fr = flat_report(E, p0, ALPHA, TEMP, lam_f, Wf, sur, int(cand.n))
        fr.update(pdb=pdb, lam=float(lam_f), n_res=int(cand.n), fold=int(cand.fold), k=int(k),
                  pad_mass=pad_mass, at="deployed_optimum_seed0",
                  label="NATIVE-FREE PROPERTY MEASUREMENT")
        rows.append(fr)
    # the zero-lam readout anchor (prereg B2.4): R_alpha vs the deployed uniform tail average
    # vs production, ORACLE-scored, at the SAME deployed optimum.
    def orc(Cc):
        if cand.nat_ca is None or not np.isfinite(np.asarray(cand.nat_ca, float)).all():
            return float("nan")
        return float(I.ca_rmsd(np.asarray(Cc, float), cand.nat_ca))
    rows.append(dict(pdb=pdb, lam=None, kind="readout_anchor", n_res=int(cand.n),
                     fold=int(cand.fold), k=int(k), at="deployed_optimum_seed0",
                     label="ORACLE-SCORED ACHIEVABLE READOUTS",
                     rmsd_prod=orc(C_prod), rmsd_tail_uniform=orc(C_uni), rmsd_R_alpha=orc(R0),
                     m_tail=int(len(tail_idx)), pad_mass=pad_mass,
                     tail_is_prefix=bool(np.array_equal(
                         np.sort(tail_idx),
                         np.sort(np.asarray(RP.topm(E_real, len(tail_idx), key), int)))),
                     lam_w_max=float(lam_w0.max()), lam_w_ratio=float(
                         lam_w0.max() / max(lam_w0[lam_w0 > 0].min(), 1e-300)),
                     rg_prod=float(np.sqrt(((C_prod - C_prod.mean(0)) ** 2).sum(1).mean())),
                     rg_R_alpha=float(np.sqrt(((R0 - R0.mean(0)) ** 2).sum(1).mean())),
                     bond_prod=float(np.linalg.norm(np.diff(C_prod, axis=0), axis=1).mean()),
                     bond_R_alpha=float(np.linalg.norm(np.diff(R0, axis=0), axis=1).mean())))
    return rows


# ======================================== the set-equality counterexample on a real instrument
def argmin_untied(vals, rng, atol=1e-12, rtol=1e-9):
    """The argmin of vals, with EXACT ties broken by a stable RNG draw, never by array order.

    The project's named failure mode (memory tie-breaking-leaks-the-pool-order, and lane D's
    S29-L38 check of S29-L25): on a pool held in DIS-sorted order, np.argmin on a tied score
    reads the POOL ORDER, which here is biased toward the energy prefix -- precisely the
    hypothesis under test.  s24.stats_lib.argmin_tied averages an OUTCOME over the tie set,
    which is the right tool when an outcome exists; here the object is an INDEX, so the tie is
    broken by an independent uniform draw and the tie-set size is returned and recorded.
    """
    vals = np.asarray(vals, float)
    m = float(np.nanmin(vals))
    tie = np.flatnonzero(np.isclose(vals, m, atol=atol, rtol=rtol))
    if len(tie) == 1:
        return int(tie[0]), 1
    return int(tie[rng.integers(len(tie))]), int(len(tie))


def sur_values(sur, Cb: np.ndarray, chunk: int = 4096) -> np.ndarray:
    """The surrogate's VALUE for a batch of clouds `Cb` (B, n_res, 3). Same arithmetic as
    `s28_A_amp.Surrogate.value_grad`, batched, value only."""
    Cb = np.asarray(Cb, float)
    B = len(Cb)
    ar = np.arange(sur.npairs)
    out = np.empty(B)
    for a in range(0, B, chunk):
        b = Cb[a:a + chunk]
        v = b[:, sur.i] - b[:, sur.j]
        d = np.linalg.norm(v, axis=2)
        u = (d - sur.g0) / sur.dg
        i0 = np.clip(np.floor(u).astype(int), 0, sur.G - 2)
        f = u - i0
        r0 = sur.risk[ar[None, :], i0]
        sl = sur.slope[ar[None, :], i0]
        out[a:a + chunk] = (r0 + sl * f * sur.dg).mean(1)
    return out


def subset_target(pdb: str, n_pool: int = 500, greedy_m: int = 5) -> List[Dict]:
    """Is the f-optimal SMALL SUBSET a prefix of the energy order? (S29-L17 section 4.5.)

    Exhaustive over all C(n_pool, 2) pairs of the DIS-ordered pool, then a greedy + local-search
    subset of size `greedy_m`. Native-free: f is the shipped distogram risk of the SUBSET AVERAGE
    and the ordering is E. The ORACLE RMSDs are attached afterwards, for reporting only.
    """
    from s27 import run_pool as RP
    from s27 import s28_A_amp as AMP
    cand, ch, _ = RP.channels_for(pdb)
    E_real = RP.zr(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    order = RP.topm(E_real, cand.k, key)
    top = np.sort(RP.topm(E_real, M_PROD, key)).astype(int)
    frame = AMP.Frame(cand.W, top)
    sur = AMP.Surrogate(I.distogram(pdb, cand.seq, cand.fold), int(cand.n))
    n_res = int(cand.n)
    idx = np.asarray(order[:min(n_pool, cand.k)], int)
    Wp = frame.Wp[idx]                                     # (N, n_res, 3), the readout's frame
    N = len(idx)
    rank_of = {int(c): r for r, c in enumerate(idx)}

    def orc(Cc):
        if cand.nat_ca is None or not np.isfinite(np.asarray(cand.nat_ca, float)).all():
            return float("nan")
        return float(I.ca_rmsd(np.asarray(Cc, float), cand.nat_ca))

    f_single = sur_values(sur, Wp)                         # the per-state criterion
    ii, jj = np.triu_indices(N, 1)
    rng_tie = RP.rng_for(pdb, "s29B_tie")
    best_v, best_pair, n_tied = np.inf, (-1, -1), 1
    CH = 20000
    all_vals = np.empty(len(ii))
    for a in range(0, len(ii), CH):
        i2, j2 = ii[a:a + CH], jj[a:a + CH]
        all_vals[a:a + CH] = sur_values(sur, 0.5 * (Wp[i2] + Wp[j2]))
    # ONE global argmin over all C(N,2) pairs, ties broken by a stable RNG draw, never by array
    # order (the chunked version could also lose a global tie that straddled two chunks)
    t, n_tied = argmin_untied(all_vals, rng_tie)
    best_v, best_pair = float(all_vals[t]), (int(ii[t]), int(jj[t]))
    prefix_pair_v = float(sur_values(sur, 0.5 * (Wp[0] + Wp[1])[None])[0])
    # the best pair by the PER-STATE criterion (the two lowest f_single)
    ps = np.argsort(f_single, kind="stable")[:2]
    perstate_pair_v = float(sur_values(sur, 0.5 * (Wp[ps[0]] + Wp[ps[1]])[None])[0])
    # greedy + one pass of local search at size greedy_m
    t0i, tie0 = argmin_untied(f_single, rng_tie)
    sel, tie_max = [int(t0i)], int(tie0)
    for _ in range(greedy_m - 1):
        cand_ids = np.array([c for c in range(N) if c not in sel])
        Cb = (Wp[sel].sum(0)[None] + Wp[cand_ids]) / (len(sel) + 1.0)
        v = sur_values(sur, Cb)
        bi, nt = argmin_untied(v, rng_tie)
        tie_max = max(tie_max, nt)
        sel.append(int(cand_ids[bi]))
    improved, n_swaps = True, 0
    while improved:
        improved = False
        cur = float(sur_values(sur, Wp[sel].mean(0)[None])[0])
        for pos in range(len(sel)):
            rest = [s for q, s in enumerate(sel) if q != pos]
            cand_ids = np.array([c for c in range(N) if c not in rest])
            Cb = (Wp[rest].sum(0)[None] + Wp[cand_ids]) / (len(rest) + 1.0)
            v = sur_values(sur, Cb)
            b, nt = argmin_untied(v, rng_tie)
            tie_max = max(tie_max, nt)
            if v[b] < cur - 1e-12:
                sel = rest + [int(cand_ids[b])]
                cur, improved, n_swaps = float(v[b]), True, n_swaps + 1
    greedy_v = float(sur_values(sur, Wp[sel].mean(0)[None])[0])
    prefix_m_v = float(sur_values(sur, Wp[:greedy_m].mean(0)[None])[0])
    return [dict(pdb=pdb, n_res=n_res, fold=int(cand.fold), N=N, kind="subset_counterexample",
                 label="NATIVE-FREE SELECTION; the RMSD columns are ORACLE",
                 best_pair_ranks=[int(best_pair[0]), int(best_pair[1])],
                 best_pair_f=best_v, n_tied_in_best_chunk=int(n_tied),
                 prefix_pair_f=prefix_pair_v, perstate_pair_ranks=[int(ps[0]), int(ps[1])],
                 perstate_pair_f=perstate_pair_v,
                 pair_is_prefix=bool(sorted(best_pair) == [0, 1]),
                 pair_gap_vs_prefix=float(prefix_pair_v - best_v),
                 greedy_m=greedy_m, greedy_ranks=sorted(int(s) for s in sel),
                 greedy_f=greedy_v, prefix_m_f=prefix_m_v,
                 greedy_is_prefix=bool(sorted(sel) == list(range(greedy_m))),
                 greedy_gap_vs_prefix=float(prefix_m_v - greedy_v),
                 rmsd_best_pair=orc(0.5 * (Wp[best_pair[0]] + Wp[best_pair[1]])),
                 rmsd_prefix_pair=orc(0.5 * (Wp[0] + Wp[1])),
                 rmsd_greedy=orc(Wp[sel].mean(0)),
                 rmsd_prefix_m=orc(Wp[:greedy_m].mean(0)),
                 rmsd_prod=orc(frame.Wf[top].mean(0).reshape(n_res, 3)),
                 rank_of_min_f_single=int(t0i), max_tie_set=int(tie_max),
                 n_swaps=int(n_swaps))]


# ================================================= MEASUREMENT 5b: THE ENDPOINT (prereg addendum 3)
def fixed_profile(E, alpha, T, n_grid: int = 4001):
    """M6 (S29-L15): the exact minimax optimum of the deployed objective over the rank ladder,
    p*(x) proportional to exp((t* - E_x)_+ / (alpha T)) with
    t* = argmax_t [ t - T log sum_x exp((t - E_x)_+ / (alpha T)) ].

    TARGET-INDEPENDENT: it depends on (alpha, T) and the ladder only, so applying it to a target's
    own DIS order involves no circuit, no optimiser and no per-target computation.  That is what
    makes it strictly stronger than "a classical equivalent" (contract rule 15).
    """
    E = np.asarray(E, float)
    ts = np.linspace(E.min() - 1.0, E.max() + 1.0, int(n_grid))
    best_t, best_g = None, -np.inf
    for t in ts:
        z = np.maximum(t - E, 0.0) / (alpha * T)
        mx = float(z.max())
        g = t - T * (mx + math.log(float(np.exp(z - mx).sum())))
        if g > best_g:
            best_g, best_t = g, float(t)
    z = np.maximum(best_t - E, 0.0) / (alpha * T)
    z = z - z.max()
    p = np.exp(z)
    return p / p.sum(), best_t, best_g


def _readouts(E, p, alpha, Wf, k, W, n_res):
    """R_alpha (the quantity the objective steers) and the DEPLOYED uniform average over the same
    tail SET (which it does not).  Returns (C_Ralpha, C_tailset, tail_idx, lam_w)."""
    lam_w, strict, x_q, _, _ = tail_lambda(E, p, alpha)
    R = (lam_w @ Wf).reshape(n_res, 3)
    tail = np.sort(np.union1d(strict, [x_q])).astype(int)
    tail = tail[tail < k]
    C_set = I.coordinate_average(np.asarray(W, float)[tail])[0]
    return R, C_set, tail, lam_w


def endpoint_target(pdb):
    """Every arm and control of measurement 5b for one target. Selection is native-free; the
    RMSD columns are ORACLE evaluations of achievable selections."""
    from s22 import qcand_lib as QC
    from s27 import run_pool as RP
    from s27 import s28_A_amp as AMP
    cand, ch, _ = RP.channels_for(pdb)
    E_real = RP.zr(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    top = np.sort(RP.topm(E_real, M_PROD, key)).astype(int)
    top_lex = np.sort(np.lexsort((key, np.asarray(ch["DIS"], float)))[:M_PROD]).astype(int)
    assert np.array_equal(top, top_lex), pdb + ": top-75 differs from lane P's lexsort path"
    frame = AMP.Frame(cand.W, top)
    sur = AMP.Surrogate(I.distogram(pdb, cand.seq, cand.fold), int(cand.n))
    enc = QC.Encoding(E_real)
    E, dim, k, n_res = enc.E, enc.dim, int(cand.k), int(cand.n)
    Wf = np.zeros((dim, 3 * n_res))
    Wf[:k] = frame.Wf
    C_prod = I.coordinate_average(cand.W[top])[0]

    def orc(Cc):
        if cand.nat_ca is None or not np.isfinite(np.asarray(cand.nat_ca, float)).all():
            return float("nan")
        Cc = np.asarray(Cc, float)
        return float(I.ca_rmsd(Cc, cand.nat_ca)) if np.isfinite(Cc).all() else float("nan")

    base = dict(pdb=pdb, n=n_res, fold=int(cand.fold), k=k, basis="point_cloud",
                label="ACHIEVABLE SELECTION, ORACLE-SCORED", rmsd_prod=orc(C_prod),
                fail18=bool(pdb in FAIL18))
    rows = []

    def emit(arm, p, extra=None):
        R, C_set, tail, lam_w = _readouts(E, p, ALPHA, Wf, k, cand.W, n_res)
        pr_real = p[:k] / max(float(p[:k].sum()), 1e-300)
        d = dict(base)
        d.update(arm=arm, m_tail=int(len(tail)),
                 pr=float(1.0 / np.sum(pr_real ** 2)), pad_mass=float(p[k:].sum()),
                 jac75=float(len(set(tail.tolist()) & set(top.tolist()))
                             / max(1, len(set(tail.tolist()) | set(top.tolist())))),
                 lam_w_ratio=float(lam_w.max() / max(float(lam_w[lam_w > 0].min()), 1e-300)),
                 rmsd_Ralpha=orc(R), rmsd_tailset=orc(C_set),
                 rg_Ralpha=float(np.sqrt(((R - R.mean(0)) ** 2).sum(1).mean())),
                 bond_Ralpha=float(np.linalg.norm(np.diff(R, axis=0), axis=1).mean()),
                 tail_is_prefix=bool(np.array_equal(
                     np.sort(tail),
                     np.sort(np.asarray(RP.topm(E_real, len(tail), key), int)))),
                 C_Ralpha=R.tolist(), C_tailset=np.asarray(C_set, float).tolist())
        if extra:
            d.update(extra)
        rows.append(d)

    for lam_f in LAM_GRID:
        for seed in (0, 1):
            t0 = time.time()
            p, th, circ, dd = run_tta_vqe(E, ALPHA, TEMP, lam_f, Wf, sur, n_res,
                                          enc.n_qubits, seed=seed)
            emit("vqe|lam%g|s%d" % (lam_f, seed), p,
                 dict(secs=float(time.time() - t0), F=float(dd["value"]),
                      cvar=float(dd["cvar_value"]), H_nats=float(dd["H"]),
                      s_val=float(dd["s_val"]), seed=seed, lam=float(lam_f), source="vqe"))
        circ = Q.StatevectorCircuit(enc.n_qubits, LAYERS)
        rng = np.random.default_rng(10000 + int(lam_f * 1000))
        best = None
        for _ in range(N_UNTRAINED):
            th = rng.normal(0.0, 0.6, circ.n_params())
            pu = circ.probs(th)
            du = dF_dp(E, pu, ALPHA, TEMP, lam_f, Wf, sur, n_res)
            if best is None or du["value"] < best[0]:
                best = (du["value"], pu)
        emit("untr16|lam%g" % lam_f, best[1],
             dict(F=float(best[0]), lam=float(lam_f), source="untrained16", seed=-1))

    p_fix, t_star, g_star = fixed_profile(E, ALPHA, TEMP)
    emit("fixed_profile|M6", p_fix,
         dict(t_star=float(t_star), g_star=float(g_star), source="fixed_profile", seed=-2,
              lam=None))

    prod_row = dict(base)
    prod_row.update(arm="production", source="production", seed=-3, lam=None,
                    m_tail=int(len(top)), pr=float(M_PROD), pad_mass=0.0, jac75=1.0,
                    lam_w_ratio=1.0, rmsd_Ralpha=orc(C_prod), rmsd_tailset=orc(C_prod),
                    rg_Ralpha=float(np.sqrt(((C_prod - C_prod.mean(0)) ** 2).sum(1).mean())),
                    bond_Ralpha=float(np.linalg.norm(np.diff(C_prod, axis=0), axis=1).mean()),
                    tail_is_prefix=True, C_Ralpha=np.asarray(C_prod, float).tolist(),
                    C_tailset=np.asarray(C_prod, float).tolist())
    rows.append(prod_row)
    return rows


def chain_main(arms, readout="C_Ralpha", shard=""):
    """Built chain for named arms, resumable per (arm, pdb, readout)."""
    from s24 import d_harness as H
    from s29 import s29_B_compat as C
    rows = C.load_all(END_ROWS)
    by = {(r["arm"], r["pdb"]): r for r in rows}
    path = C.shard_path(CHAIN_ROWS, shard)
    done = set()
    for r in C.load_all(CHAIN_ROWS):
        done.add((r["arm"], r["pdb"], r["readout"]))
    pdbs = C._shard(sorted(set(r["pdb"] for r in rows)), shard)
    t0 = time.time()
    for i, pdb in enumerate(pdbs):
        todo = [a for a in arms if (a, pdb, readout) not in done and (a, pdb) in by]
        if not todo:
            continue
        cand = H.Candidates.from_universe(pdb, k=500)
        nat_ok = cand.nat_ca is not None and np.isfinite(np.asarray(cand.nat_ca, float)).all()
        out = []
        for a in todo:
            r = by[(a, pdb)]
            Cc = np.asarray(r[readout], float)
            t1 = time.time()
            ca = H.readout_projected(cand, Cc)
            out.append(dict(arm=a, pdb=pdb, readout=readout, n=int(cand.n), fold=int(cand.fold),
                            basis="built_chain", label="ACHIEVABLE, ORACLE-SCORED",
                            rmsd_cloud=float(I.ca_rmsd(Cc, cand.nat_ca)) if nat_ok else float("nan"),
                            rmsd_chain=float(I.ca_rmsd(ca, cand.nat_ca)) if nat_ok else float("nan"),
                            fail18=bool(pdb in FAIL18), secs=float(time.time() - t1)))
        with open(path, "a", encoding="utf-8") as fh:
            for r in out:
                fh.write(json.dumps(r) + "\n")
        print("  [chain %d/%d] %s %d arms %.1fs (elapsed %.1f min)"
              % (i + 1, len(pdbs), pdb, len(out), sum(r["secs"] for r in out),
                 (time.time() - t0) / 60.0), flush=True)
    print("done:", path, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--flat", action="store_true")
    ap.add_argument("--subset", action="store_true")
    ap.add_argument("--endpoint", action="store_true")
    ap.add_argument("--chain", type=str, default="")
    ap.add_argument("--readout", type=str, default="C_Ralpha")
    ap.add_argument("--targets", type=str, default="12")
    ap.add_argument("--shard", type=str, default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--analyse-flat", action="store_true")
    a = ap.parse_args(argv)
    if a.flat:
        from s29 import s29_B_compat as C
        pdbs = C.picks_12()
        if a.limit:
            pdbs = pdbs[:a.limit]
        C.run_phase("flat", pdbs, FLAT_ROWS, flat_target)
    if a.subset:
        from s25 import phys_lib as P
        from s29 import s29_B_compat as C
        pdbs = C.picks_12() if a.targets == "12" else list(P.targets())
        if a.limit:
            pdbs = pdbs[:a.limit]
        pdbs = C._shard(pdbs, a.shard)
        C.run_phase("subset", pdbs, C.shard_path(SUBSET_ROWS, a.shard), subset_target)
    if a.endpoint:
        from s25 import phys_lib as P
        from s29 import s29_B_compat as C
        pdbs = C.picks_12() if a.targets == "12" else list(P.targets())
        if a.limit:
            pdbs = pdbs[:a.limit]
        pdbs = C._shard(pdbs, a.shard)
        C.run_phase("end", pdbs, C.shard_path(END_ROWS, a.shard), endpoint_target)
    if a.chain:
        chain_main([x for x in a.chain.split(",") if x], a.readout, a.shard)
    if a.analyse_flat:
        analyse_flat()
    return 0


def analyse_flat() -> Dict:
    from s24 import stats_lib as ST
    from s29 import s29_B_compat as C
    rows = C.load_all(FLAT_ROWS)
    flat = [r for r in rows if r.get("kind") != "readout_anchor"]
    anch = [r for r in rows if r.get("kind") == "readout_anchor"]
    pdbs = sorted({r["pdb"] for r in flat})
    out = dict(kind="s29_B_tta_flat", lane="B", sprint=29, n=len(pdbs), pdbs=pdbs,
               lam_grid=list(LAM_GRID), by_lam={},
               note="native-free property measurement; the readout anchor is ORACLE-scored")
    for lam in LAM_GRID:
        rs = [r for r in flat if r["lam"] == lam]
        if not rs:
            continue
        cell = dict(n=len(rs))
        for kk in ("flat_cvar", "flat_entropy", "flat_f", "flat_total", "flat_info_bearing",
                   "flat_readout_tta", "flat_readout_deployed", "overlap_tta",
                   "m_strict", "nmove_readout_tta", "dC_norm_tail_mean"):
            v = np.array([r[kk] for r in rs], float)
            cell[kk] = float(np.nanmedian(v))
            cell[kk + "_min"] = float(np.nanmin(v))
            cell[kk + "_max"] = float(np.nanmax(v))
        out["by_lam"][str(lam)] = cell
    if anch:
        pd = sorted({r["pdb"] for r in anch})
        folds = ST.pinned_folds(pd)
        v_prod = np.array([[r for r in anch if r["pdb"] == p][0]["rmsd_prod"] for p in pd])
        for nm in ("rmsd_tail_uniform", "rmsd_R_alpha"):
            v = np.array([[r for r in anch if r["pdb"] == p][0][nm] for p in pd])
            o = ST.compare(v, v_prod, folds=folds, names=pd,
                           label=f"{nm} - production (point cloud, ORACLE)")
            out.setdefault("readout_anchor", {})[nm] = dict(
                mean=float(v.mean()), effect=o["effect"], se=o["se"], mde=o["mde"],
                x_mde=o["effect_over_mde"], ci_fold=o["ci95_fold"], verdict=o["verdict"],
                fmt=ST.fmt(o))
        out["readout_anchor"]["prod_mean"] = float(v_prod.mean())
        out["readout_anchor"]["lam_w_ratio_median"] = float(
            np.median([r["lam_w_ratio"] for r in anch]))
        out["readout_anchor"]["rg_ratio_median"] = float(
            np.median([r["rg_R_alpha"] / r["rg_prod"] for r in anch]))
        out["readout_anchor"]["bond_ratio_median"] = float(
            np.median([r["bond_R_alpha"] / r["bond_prod"] for r in anch]))
    ST.save_atomic(FLAT_OUT, out)
    print(json.dumps(out, indent=1)[:5000])
    return out


if __name__ == "__main__":
    raise SystemExit(main())
