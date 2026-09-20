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
    # frame.Wf[top].mean(0) is the same object in the top-75 medoid frame; assert they agree
    assert abs(float(I.ca_rmsd(C_prod, frame.Wf[top].mean(0).reshape(int(cand.n), 3)))) < 1e-8
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
    best_v, best_pair = np.inf, (-1, -1)
    CH = 20000
    for a in range(0, len(ii), CH):
        i2, j2 = ii[a:a + CH], jj[a:a + CH]
        vals = sur_values(sur, 0.5 * (Wp[i2] + Wp[j2]))
        t = int(np.argmin(vals))
        # ties: average over the tied argmin set rather than reading array order
        tie = np.flatnonzero(vals <= vals[t] + 1e-15)
        if vals[t] < best_v:
            best_v = float(vals[t])
            best_pair = (int(i2[tie[0]]), int(j2[tie[0]]))
            n_tied = int(len(tie))
    prefix_pair_v = float(sur_values(sur, 0.5 * (Wp[0] + Wp[1])[None])[0])
    # the best pair by the PER-STATE criterion (the two lowest f_single)
    ps = np.argsort(f_single, kind="stable")[:2]
    perstate_pair_v = float(sur_values(sur, 0.5 * (Wp[ps[0]] + Wp[ps[1]])[None])[0])
    # greedy + one pass of local search at size greedy_m
    sel = [int(np.argmin(f_single))]
    for _ in range(greedy_m - 1):
        cand_ids = np.array([c for c in range(N) if c not in sel])
        Cb = (Wp[sel].sum(0)[None] + Wp[cand_ids]) / (len(sel) + 1.0)
        v = sur_values(sur, Cb)
        sel.append(int(cand_ids[int(np.argmin(v))]))
    improved = True
    while improved:
        improved = False
        cur = float(sur_values(sur, Wp[sel].mean(0)[None])[0])
        for pos in range(len(sel)):
            rest = [s for q, s in enumerate(sel) if q != pos]
            cand_ids = np.array([c for c in range(N) if c not in rest])
            Cb = (Wp[rest].sum(0)[None] + Wp[cand_ids]) / (len(rest) + 1.0)
            v = sur_values(sur, Cb)
            b = int(np.argmin(v))
            if v[b] < cur - 1e-12:
                sel = rest + [int(cand_ids[b])]
                cur, improved = float(v[b]), True
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
                 rank_of_min_f_single=int(np.argmin(f_single)))]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--flat", action="store_true")
    ap.add_argument("--subset", action="store_true")
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
        from s29 import s29_B_compat as C
        pdbs = C.picks_12()
        if a.limit:
            pdbs = pdbs[:a.limit]
        C.run_phase("subset", pdbs, SUBSET_ROWS, subset_target)
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
