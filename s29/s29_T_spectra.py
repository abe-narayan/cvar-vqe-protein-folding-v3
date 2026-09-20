#!/usr/bin/env python
"""s29/s29_T_spectra.py -- LANE T, section 3c: the spectra of three pool compatibility
matrices, their STABLE RANKS as a function of register size, and the hop-only gradient
variance of the CENTERED observables at n = 4..9.

PROPERTY MEASUREMENT.  Reads no native, no RMSD, no oracle field.  Nothing here is
deployable and nothing is selected by any of it; it exists to test one derivation
(`s29/THEORY.md` section 3):

    Var_theta[ dF/dtheta_k ]  ~=  ||A - (tr A / D) I||_F^2 / D^2        (2-design leading term)
                              =   r_stable(A) / D^2   when ||A||_2 = 1,

so the decay of a unit-spectral-norm observable's gradient variance is set by its STABLE
RANK r_stable = ||A||_F^2 / ||A||_2^2 and by nothing else about it (not its rank-one-ness,
not its off-diagonality).  The prediction that matters for S29: centering the graph (removing
the Perron/typicality mode) leaves the stable rank BOUNDED, so the -2 per qubit decay
survives centering.

THE THREE MATRICES (all native-free, all on the same sub-pool of the 2^n best DIS candidates)
    A     the S28 Gaussian similarity graph, A_ij = exp(-d_ij^2 / 2 sigma^2), zero diagonal,
          sigma = median off-diagonal d, normalised to unit spectral norm (S28 lane B).
    A_c   its double centering, H A H with H = I - 11^T/N, unit spectral norm: the graph with
          the uniform (typicality) mode projected out.
    G     the SIGNED AGREEMENT matrix, G = Delta Delta^T / n_res with Delta the members'
          deviations from the pool mean in the medoid frame: inner products of deviations,
          positive for members that deviate the same way.  Rank <= 3 n_res, by construction.

USAGE
    python s29/s29_T_spectra.py --spectra [--limit N]     phase 1 (spectra, seconds/target)
    python s29/s29_T_spectra.py --grad    [--limit N]     phase 2 (gradient variance)
    python s29/s29_T_spectra.py --analyse                 -> s29/results/s29_T_spectra.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
SPEC_ROWS = os.path.join(RESULTS, "s29_T_spectra_rows.jsonl")
GRAD_ROWS = os.path.join(RESULTS, "s29_T_grad_rows.jsonl")
OUT = os.path.join(RESULTS, "s29_T_spectra.json")

NS = [4, 5, 6, 7, 8, 9]
N_THETA = 120          # S28 lane B's draw count, same seed law, so the rows are comparable
SEED = 1009
INIT_SD = 0.6


# ------------------------------------------------------------------ the 12 targets (S27 T11)
def picks():
    from s25 import phys_lib as P
    return P.targets()[::11][:12]


def _done(path, key="pdb"):
    s = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    s.add(json.loads(line)[key])
                except Exception:
                    pass
    return s


def _append(path, rows):
    with open(path, "a", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")


# ------------------------------------------------------------------------ the three matrices
def spectrum_stats(M, name):
    """Eigen-statistics of a real symmetric matrix, in the convention the theory uses."""
    M = np.asarray(M, float)
    N = len(M)
    w = np.linalg.eigvalsh(M)
    o = np.argsort(-np.abs(w))
    wa = w[o]                                      # by decreasing |lambda|
    lam1 = float(wa[0])
    fro2 = float((w ** 2).sum())
    tr = float(np.trace(M))
    Mt = M - (tr / N) * np.eye(N)                  # the traceless part (the identity has no gradient)
    fro2_t = float((Mt ** 2).sum())
    return dict(name=name, N=N, lam_max=float(w[-1]), lam_min=float(w[0]),
                lam1_abs=lam1, lam2_over_lam1=float(wa[1] / wa[0]),
                lam3_over_lam1=float(wa[2] / wa[0]),
                trace=tr, fro2=fro2,
                r_stable=float(fro2 / lam1 ** 2),
                r_stable_traceless=float(fro2_t / lam1 ** 2),
                # participation of the eigenvalue-squared distribution: an effective rank
                eff_rank_sq=float(fro2 ** 2 / (w ** 4).sum()),
                top_share=float(wa[0] ** 2 / fro2))


def perron_stats(M):
    w, v = np.linalg.eigh(np.asarray(M, float))
    v1 = v[:, -1]
    N = len(M)
    u = np.ones(N) / np.sqrt(N)
    return dict(perron_uniform_overlap2=float((v1 @ u) ** 2),
                perron_pr_over_dim=float(1.0 / (N * (v1 ** 4).sum())))


def double_center(M):
    M = np.asarray(M, float)
    N = len(M)
    r = M.mean(1, keepdims=True)
    c = M.mean(0, keepdims=True)
    return M - r - c + M.mean()


def deviation_gram(W, D):
    """G = Delta Delta^T / n_res with Delta the deviations from the pool mean in the MEDOID
    frame (every member superposed onto the medoid; the frame the readout averages in)."""
    from s12 import instrument as I
    W = np.asarray(W, float)
    b = int(np.argmin(D.mean(1)))                 # medoid; float ties are not expected, no key needed
    Sup = I.superpose_batch(W, W[b])              # (k, n_res, 3), all in the medoid frame
    c = Sup.mean(0)
    Delta = (Sup - c).reshape(len(Sup), -1)       # (k, 3 n_res)
    n_res = W.shape[1]
    return (Delta @ Delta.T) / float(n_res), b, Delta


def band_stats(D):
    N = len(D)
    iu = np.triu_indices(N, 1)
    d = D[iu]
    return dict(d_min=float(d.min()), d_max=float(d.max()), d_mean=float(d.mean()),
                d_median=float(np.median(d)), d_sd=float(d.std()),
                d2_mean=float((d ** 2).mean()), cv=float(d.std() / d.mean()))


# ------------------------------------------------------------------------------- phase 1
def spectra_target(pdb):
    from s27 import run_pool as RP
    from s27 import s28_B_hop as B
    cand, ch, _ = RP.channels_for(pdb)
    E_full = RP.zr(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    order = RP.topm(E_full, cand.k, key)
    Dfull = B.pairwise_rmsd_matrix(cand.W)
    rows = []
    for n in NS:
        dim = 1 << n
        sub = order[:min(dim, cand.k)]
        Ds = Dfull[np.ix_(sub, sub)]
        g = B.kernel_graph(Ds)
        A = g["A"]                                   # unit spectral norm, zero diagonal
        Ac = double_center(A)
        Ac = Ac / float(np.abs(np.linalg.eigvalsh(Ac)).max())
        G, medoid, Delta = deviation_gram(cand.W[sub], Ds)
        Gn = G / float(np.abs(np.linalg.eigvalsh(G)).max())
        # the pool's own shape spectrum (the deviation covariance), for the 3a prediction
        gw = np.linalg.eigvalsh(G)
        gw = np.sort(gw)[::-1]
        f1 = float(gw[0] / gw.sum())
        bs = band_stats(Ds)
        row = dict(pdb=pdb, n=n, dim=dim, k_used=int(len(sub)), n_res=int(cand.n),
                   sigma=float(g["sigma"]), lam_max_raw=float(g["lam_max"]),
                   band=bs,
                   shape_f1=f1,
                   shape_eff_dim=float(gw.sum() ** 2 / (gw ** 2).sum()),
                   pred_lam2_over_lam1=float(f1 * bs["d2_mean"] / (2.0 * g["sigma"] ** 2)),
                   A=spectrum_stats(A, "A"), A_c=spectrum_stats(Ac, "A_c"),
                   G=spectrum_stats(Gn, "G"))
        row["A"].update(perron_stats(A))
        row["G"].update(perron_stats(Gn))
        row["A_c"].update(perron_stats(Ac))
        rows.append(row)
    return rows


# ------------------------------------------------------------------------------- phase 2
def grad_target(pdb):
    """hop-only gradient variance for A, A_c and G at n = 4..9, S28's draw law exactly."""
    from s27 import run_pool as RP
    from s27 import s28_B_hop as B
    cand, ch, _ = RP.channels_for(pdb)
    E_full = RP.zr(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    order = RP.topm(E_full, cand.k, key)
    Dfull = B.pairwise_rmsd_matrix(cand.W)
    rows = []
    for n in NS:
        dim = 1 << n
        sub = order[:min(dim, cand.k)]
        Ds = Dfull[np.ix_(sub, sub)]
        g = B.kernel_graph(Ds)
        A = g["A"]
        Ac = double_center(A)
        Ac = Ac / float(np.abs(np.linalg.eigvalsh(Ac)).max())
        G, _, _ = deviation_gram(cand.W[sub], Ds)
        Gn = G / float(np.abs(np.linalg.eigvalsh(G)).max())
        E = B.deployed_E(dim)
        for nm, M in (("A", A), ("A_c", Ac), ("G", Gn)):
            Mp = M if len(M) == dim else B.pad_graph(M, dim)
            r = B.measure_hop(n, B.LAYERS, B.ALPHA, B.TEMP, 1.0, N_THETA, SEED, INIT_SD,
                              E, Mp, "hop_only")
            r.update(pdb=pdb, matrix=nm, fro2=float((Mp ** 2).sum()),
                     r_stable=float((Mp ** 2).sum() / max(float(np.abs(np.linalg.eigvalsh(Mp)).max()) ** 2, 1e-30)),
                     pred_var=float((Mp ** 2).sum() / dim ** 2))
            rows.append(r)
    return rows


# ------------------------------------------------------------------------------- analysis
def analyse():
    import numpy as np
    srows = [json.loads(l) for l in open(SPEC_ROWS, encoding="utf-8")] if os.path.exists(SPEC_ROWS) else []
    grows = [json.loads(l) for l in open(GRAD_ROWS, encoding="utf-8")] if os.path.exists(GRAD_ROWS) else []
    out = dict(kind="s29_T_spectra", lane="T", sprint=29,
               note="property measurement; no native, no RMSD", ns=NS)
    # --- spectra, median over targets per (n, matrix)
    spec = {}
    for n in NS:
        cell = {}
        for mat in ("A", "A_c", "G"):
            rs = [r[mat] for r in srows if r["n"] == n]
            if not rs:
                continue
            cell[mat] = {k: float(np.median([r[k] for r in rs]))
                         for k in ("lam2_over_lam1", "lam3_over_lam1", "r_stable",
                                   "r_stable_traceless", "eff_rank_sq", "top_share",
                                   "perron_uniform_overlap2", "perron_pr_over_dim")}
            cell[mat]["n_targets"] = len(rs)
        base = [r for r in srows if r["n"] == n]
        if base:
            cell["pool"] = dict(
                sigma=float(np.median([r["sigma"] for r in base])),
                d_min=float(np.median([r["band"]["d_min"] for r in base])),
                d_max=float(np.median([r["band"]["d_max"] for r in base])),
                d_mean=float(np.median([r["band"]["d_mean"] for r in base])),
                d_sd=float(np.median([r["band"]["d_sd"] for r in base])),
                cv=float(np.median([r["band"]["cv"] for r in base])),
                d_max_over_sigma=float(np.median([r["band"]["d_max"] / r["sigma"] for r in base])),
                shape_f1=float(np.median([r["shape_f1"] for r in base])),
                shape_eff_dim=float(np.median([r["shape_eff_dim"] for r in base])),
                pred_lam2_over_lam1=float(np.median([r["pred_lam2_over_lam1"] for r in base])),
                obs_lam2_over_lam1=float(np.median([r["A"]["lam2_over_lam1"] for r in base])),
                pred_over_obs=float(np.median([r["pred_lam2_over_lam1"] / r["A"]["lam2_over_lam1"]
                                               for r in base])))
        spec[str(n)] = cell
    out["spectra"] = spec
    # --- gradient variance: measured vs predicted, and the fitted slope per matrix
    gr = {}
    for mat in ("A", "A_c", "G"):
        per_n, meas, pred = {}, [], []
        for n in NS:
            rs = [r for r in grows if r["matrix"] == mat and r["n"] == n]
            if not rs:
                continue
            m = float(np.median([r["var_g0"] for r in rs]))
            p = float(np.median([r["pred_var"] for r in rs]))
            per_n[str(n)] = dict(var_g0_median=m, pred_var_median=p,
                                 ratio_meas_over_pred=float(m / p) if p > 0 else None,
                                 r_stable_median=float(np.median([r["r_stable"] for r in rs])),
                                 n_targets=len(rs))
            meas.append((n, m))
            pred.append((n, p))
        if len(meas) > 2:
            ns_ = np.array([x[0] for x in meas], float)
            v = np.array([x[1] for x in meas], float)
            pv = np.array([x[1] for x in pred], float)
            gr[mat] = dict(per_n=per_n,
                           log2_slope_measured=float(np.polyfit(ns_, np.log2(v), 1)[0]),
                           log2_slope_predicted=float(np.polyfit(ns_, np.log2(pv), 1)[0]),
                           log2_slope_measured_48=float(np.polyfit(ns_[ns_ <= 8], np.log2(v[ns_ <= 8]), 1)[0]))
        else:
            gr[mat] = dict(per_n=per_n)
    out["gradient"] = gr
    from s24 import stats_lib as ST
    ST.save_atomic(OUT, out)
    print(json.dumps({"spectra_n9": spec.get("9"), "gradient": {k: {kk: vv for kk, vv in v.items() if kk != "per_n"}
                                                                for k, v in gr.items()}}, indent=1))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--spectra", action="store_true")
    ap.add_argument("--grad", action="store_true")
    ap.add_argument("--analyse", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)
    if a.spectra or a.grad:
        pdbs = picks()
        if a.limit:
            pdbs = pdbs[:a.limit]
        for phase, path, fn in (("spectra", SPEC_ROWS, spectra_target),
                                ("grad", GRAD_ROWS, grad_target)):
            if not getattr(a, phase):
                continue
            done = _done(path)
            for pdb in pdbs:
                if pdb in done:
                    print(f"[{phase}] {pdb} done", flush=True)
                    continue
                t0 = time.time()
                _append(path, fn(pdb))
                print(f"[{phase}] {pdb} {time.time() - t0:.1f}s", flush=True)
    if a.analyse:
        analyse()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
