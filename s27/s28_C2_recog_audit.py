#!/usr/bin/env python
"""s27/s28_C2_recog_audit.py -- S28 lane C2: the recognition audit (`s27/PREREG_S28_C2.md`).

ORACLE DIAGNOSTIC THROUGHOUT.  Which native-free scorer in the S27 library, if any, prefers the
ORACLE structures the amplitude family can express (S28-L1b: 0.29 A best-of-5, 0.36 to 0.42
single start; the random 27-subspace at 0.61) to the production average?  Every structure
except PROD and the DIRECTIONS of the two controls is chosen against the native.

Structures per target (all (n, 3) CA clouds in lane A's frame):
  PROD          the production average (`s27/results/s28_A_structs/<pdb>.npz :: prod`)
  ORACLE circ_best   `:: oracle_circ` (best of 5 starts)
  ORACLE circ_s0     regenerated single-start seed-0 optimum (lane A's `oracle_circuit_ceiling`)
  ORACLE sub0        regenerated random-27-subspace least squares, subspace 0
  NATIVE(aff500)     `:: oracle_aff500` (the native to ~1e-7 by dimension counting)
  RAND_SIGNED[d]     control (i): random signed affine combination at circ_best's distance from PROD
  GAUSS_0.3[d]       control (ii): PROD + Gaussian noise at circ_best's RMSD to the native
  GAUSS_MATCHED[d]   control (ii-b): PROD + Gaussian noise at circ_best's distance from PROD

Scorers: the CA-level S27 channels through a one-Context adapter (`ca_scores`), three
pool-relative adapters (CONS_POOL, DMAP_CONS_POOL, POOLGO_POOL), and, through the production
projection, the backbone channels (RAMA, DSSPHB, ELEC, TORS_CONS_POOL, the 11 LEG terms and
their total).  Lower is better everywhere.

    python s27/s28_C2_recog_audit.py selftest
    python s27/s28_C2_recog_audit.py ca    [--limit N]     # CA-level audit, rows per target
    python s27/s28_C2_recog_audit.py chain [--limit N]     # built-chain audit (8 projections/target)
    python s27/s28_C2_recog_audit.py analyse_ca | analyse_chain
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import ham_lib as HL              # noqa: E402
from s27 import run_pool as RP             # noqa: E402

RESULTS = os.path.join(HERE, "results")
STRUCTS = os.path.join(RESULTS, "s28_A_structs")
A_ROWS = os.path.join(RESULTS, "s28_A_oracle_rows.jsonl")
CA_ROWS = os.path.join(RESULTS, "s28_C2_ca_rows.jsonl")
CHAIN_ROWS = os.path.join(RESULTS, "s28_C2_chain_rows.jsonl")
SALT = "s28C2"
N_DRAWS = 4
LADDER = ["PROD", "sub0", "circ_s0", "circ_best", "NATIVE"]
CONTROLS = ["RAND_SIGNED", "GAUSS_0.3", "GAUSS_MATCHED"]
CA_SCORERS = ["DIS", "DIS_MEAN", "CONTACT_LL", "DISTPOT", "CONTACT", "ENV", "HP", "RG_LAW", "RG_UNIV",
              "EXVOL", "CAGEO", "SS_MATCH", "CONS_POOL", "DMAP_CONS_POOL", "POOLGO_POOL"]
CHAIN_SCORERS = ["RAMA", "DSSPHB", "ELEC", "TORS_CONS_POOL", "LEG"] + ["LEG_" + t for t in HL.LEG_TERMS]
CA_FUNS = {"DIS_MEAN": HL.h_dis_mean, "CONTACT_LL": HL.h_contact_ll, "DISTPOT": HL.h_distpot,
           "CONTACT": HL.h_contact, "ENV": HL.h_env, "HP": HL.h_hp, "RG_LAW": HL.h_rg_law,
           "RG_UNIV": HL.h_rg_univ, "EXVOL": HL.h_exvol, "CAGEO": HL.h_cageo, "SS_MATCH": HL.h_ss_match}
RAMA = np.load(os.path.join(ROOT, "s8", "generate_rama.npz"))["cnt"]


# ============================================================================ adapters
def make_context(cand, W, universe, dg, PHI=None, PSI=None):
    """A `ham_lib.Context` whose 'pool' is the stack of structures to score.  NATIVE-FREE: the
    universe fits and the distogram are the target's own; nothing reads `nat_ca`."""
    W = np.asarray(W, float)
    k = len(W)
    c = SimpleNamespace(pdb=cand.pdb, seq=cand.seq, n=int(cand.n), fold=int(cand.fold), k=k, W=W,
                        PHI=np.zeros((k, cand.n)) if PHI is None else np.asarray(PHI, float),
                        PSI=np.zeros((k, cand.n)) if PSI is None else np.asarray(PSI, float))
    return HL.Context(c, universe, dg, RAMA[int(cand.fold)])


def ca_scores(cand, W, universe, dg):
    """Every CA-level scorer of `W` (m, n, 3): dict name -> (m,) values, lower is better."""
    cx = make_context(cand, W, universe, dg)
    out = {}
    i, j = I.pair_index(cand.n)
    D = I.pair_dists(np.asarray(W, float), i, j)
    out["DIS"] = np.asarray(I.shipped_score(dg, D.astype(np.float32).astype(float)), float)
    for nm, fn in CA_FUNS.items():
        out[nm] = np.asarray(fn(cx), float)
    pool = np.asarray(cand.W, float)
    out["CONS_POOL"] = np.array([float(I.kabsch_rmsd_batch(pool, w).mean()) for w in np.asarray(W, float)])
    i2, j2 = I.pair_index(cand.n, 2)
    Dp = I.pair_dists(pool, i2, j2)
    med = np.median(Dp, axis=0)
    D2 = I.pair_dists(np.asarray(W, float), i2, j2)
    out["DMAP_CONS_POOL"] = np.abs(D2 - med[None, :]).mean(1)
    sep = j2 - i2
    m3 = sep >= 3
    f = (Dp[:, m3] < 8.0).mean(0)
    out["POOLGO_POOL"] = -((D2[:, m3] < 8.0) * f[None, :]).sum(1)
    return out


def chain_scores(cand, CA, PHI, PSI, universe, dg):
    """The backbone scorers of projected chains: RAMA, DSSPHB, ELEC, TORS_CONS_POOL, LEG terms."""
    cx = make_context(cand, CA, universe, dg, PHI, PSI)
    out = {"RAMA": np.asarray(HL.h_rama(cx), float), "DSSPHB": np.asarray(HL.h_dssphb(cx), float),
           "ELEC": np.asarray(HL.h_elec(cx), float)}
    # TORS_CONS against the POOL's circular means (the S27 channel uses the pool it scores)
    def cmean(x):
        return np.arctan2(np.sin(x).mean(0), np.cos(x).mean(0))
    mphi, mpsi = cmean(np.asarray(cand.PHI, float)), cmean(np.asarray(cand.PSI, float))
    dphi = np.abs(np.angle(np.exp(1j * (np.asarray(PHI, float) - mphi[None, :]))))
    dpsi = np.abs(np.angle(np.exp(1j * (np.asarray(PSI, float) - mpsi[None, :]))))
    out["TORS_CONS_POOL"] = (dphi + dpsi).sum(1)
    lt = HL.legacy_terms(cx)
    tot = np.zeros(len(CA))
    for k, v in lt.items():
        out[k] = np.asarray(v, float)
        tot = tot + out[k]
    out["LEG"] = tot
    return out


# ============================================================================ structures
def load_A(pdb):
    z = np.load(os.path.join(STRUCTS, f"{pdb}.npz"))
    S = {"PROD": np.asarray(z["prod"], float), "circ_best": np.asarray(z["oracle_circ"], float),
         "NATIVE": np.asarray(z["oracle_aff500"], float)}
    row = None
    with open(A_ROWS, encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            if r["pdb"] == pdb:
                row = r
                break
    return S, row


def geom_of(C):
    """Mean virtual CA-CA bond and radius of gyration (native-free descriptors)."""
    C = np.asarray(C, float)
    return dict(bond=float(np.linalg.norm(np.diff(C, axis=0), axis=1).mean()),
                rg=float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean())))


def rmsd_between(A, B):
    return float(I.ca_rmsd(np.asarray(A, float), np.asarray(B, float)))


def scale_to_distance(base, target_struct, dist, lo=0.0, hi=64.0, tol=1e-7):
    """C(s) = (1 - s) base + s target; find s >= 0 with CA-RMSD(C(s), base) == dist (bisection;
    the RMSD after superposition is monotone in s for a fixed direction)."""
    base = np.asarray(base, float); target_struct = np.asarray(target_struct, float)

    def f(s):
        return rmsd_between((1 - s) * base + s * target_struct, base) - dist
    while f(hi) < 0 and hi < 1e6:
        hi *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    s = 0.5 * (lo + hi)
    return (1 - s) * base + s * target_struct, s


def rand_signed(frame, prod, dist, pdb, draw, cols=27):
    """Control (i): a random signed affine combination of the posed windows, moved to `dist`
    from PROD along its own displacement.  NATIVE-FREE direction; ORACLE scale."""
    rng = SD.stable_rng(pdb, "randsigned", int(draw), salt=SALT)
    B = rng.normal(0.0, 1.0, (frame.k, cols))
    z = rng.normal(0.0, 1.0, cols)
    w = B @ z
    w = w / w.sum() if abs(w.sum()) > 1e-9 else (w - w.mean() + 1.0 / frame.k)
    Cw = (w @ frame.Wf).reshape(frame.n, 3)
    C, s = scale_to_distance(prod, Cw, dist)
    return C, dict(s=float(s), frac_neg=float((w < 0).mean()))


def gauss_perturb(prod, dist, pdb, draw, tag):
    rng = SD.stable_rng(pdb, "gauss" + tag, int(draw), salt=SALT)
    noise = rng.normal(0.0, 1.0, np.asarray(prod).shape)
    C, s = scale_to_distance(prod, np.asarray(prod, float) + noise, dist)
    return C, dict(s=float(s))


def build_structures(pdb, draw_offset=0):
    """All structures for one target and their ORACLE RMSDs.  Reads the native (ORACLE) for the
    labels and for the two regenerated ORACLE structures.  `draw_offset` selects the control
    draws (0 to 3 registered; 4 to 7 the second seed)."""
    from s27 import s28_A_amp as A
    from core import quantum as Q
    cand, dis, top, dg = A.load_pool(pdb)
    frame = A.Frame(cand.W, top)
    S, arow = load_A(pdb)
    from s24 import d_harness as H
    Cd, _ = H.readout_uniform(cand, top)
    assert np.abs(S["PROD"] - Cd).max() < 1e-10, f"{pdb}: PROD differs from the deployed average"
    circ = Q.StatevectorCircuit(A.N_QUBITS, A.LAYERS)
    oc = A.oracle_circuit_ceiling(circ, frame, cand.nat_ca, starts=1)         # ORACLE
    S["circ_s0"] = np.asarray(oc["C"], float)
    r_ref = float(arow["arms"]["oracle_circ"]["per_start"][0])
    assert abs(oc["rmsd_cloud"] - r_ref) < 1e-6, f"{pdb}: circ_s0 {oc['rmsd_cloud']} != A's per_start[0] {r_ref}"
    rs, Cs, ws = A.oracle_subspace_ls(frame, cand.nat_ca, A.subspace_matrix(pdb, 0, frame.k))  # ORACLE
    S["sub0"] = np.asarray(Cs, float)
    assert abs(rs - float(arow["arms"]["oracle_sub"]["per_sub"][0])) < 1e-6, f"{pdb}: sub0 mismatch"
    d_native = rmsd_between(S["circ_best"], cand.nat_ca)                        # ORACLE scale
    d_prod = rmsd_between(S["circ_best"], S["PROD"])                            # ORACLE scale
    meta = {"d_circ_best_native": d_native, "d_circ_best_prod": d_prod}
    for d in range(draw_offset, draw_offset + N_DRAWS):
        S[f"RAND_SIGNED[{d}]"], m = rand_signed(frame, S["PROD"], d_prod, pdb, d)
        meta[f"rand_signed_{d}"] = m
        S[f"GAUSS_0.3[{d}]"], _ = gauss_perturb(S["PROD"], d_native, pdb, d, "03")
        S[f"GAUSS_MATCHED[{d}]"], _ = gauss_perturb(S["PROD"], d_prod, pdb, d, "m")
    oracle_rmsd = {k: rmsd_between(v, cand.nat_ca) for k, v in S.items()}       # ORACLE labels
    return cand, dg, frame, S, meta, oracle_rmsd


def names_in_order(S):
    return [k for k in LADDER if k in S] + [k for k in S if k not in LADDER]


# ============================================================================ per-target rows
def ca_row(pdb, draw_offset=0):
    t0 = time.time()
    cand, dg, frame, S, meta, orr = build_structures(pdb, draw_offset)
    u = I.load_univ(pdb)
    names = names_in_order(S)
    W = np.stack([S[k] for k in names])
    sc = ca_scores(cand, W, u, dg)
    del u
    row = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), names=names, meta=meta,
               oracle_rmsd={k: float(orr[k]) for k in names},
               geom={k: geom_of(S[k]) for k in names},
               scores={nm: [float(x) for x in v] for nm, v in sc.items()},
               d_prod={k: rmsd_between(S[k], S["PROD"]) for k in names}, secs=float(time.time() - t0))
    return row, S


def chain_row(pdb, S=None, cand=None, dg=None):
    from s24 import d_harness as H
    t0 = time.time()
    if S is None:
        cand, dg, frame, S, meta, orr = build_structures(pdb)
    names = [k for k in LADDER if k in S] + ["RAND_SIGNED[0]", "GAUSS_0.3[0]", "GAUSS_MATCHED[0]"]
    CA, PHI, PSI = [], [], []
    for k in names:
        pr = I.project(np.asarray(S[k], float), cand.seq, cand.fold)
        CA.append(np.asarray(pr["ca"], float)); PHI.append(np.asarray(pr["phi"], float)); PSI.append(np.asarray(pr["psi"], float))
    CA, PHI, PSI = np.stack(CA), np.stack(PHI), np.stack(PSI)
    u = I.load_univ(pdb)
    sc = chain_scores(cand, CA, PHI, PSI, u, dg)
    sc_ca = ca_scores(cand, CA, u, dg)                       # the CA scorers on the projected chains
    del u
    row = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), names=names, basis="built_chain",
               oracle_rmsd_chain={k: rmsd_between(CA[a], cand.nat_ca) for a, k in enumerate(names)},   # ORACLE
               oracle_rmsd_cloud={k: rmsd_between(S[k], cand.nat_ca) for k in names},                  # ORACLE
               scores={nm: [float(x) for x in v] for nm, v in sc.items()},
               scores_ca_on_chain={nm: [float(x) for x in v] for nm, v in sc_ca.items()},
               geom_chain={k: geom_of(CA[a]) for a, k in enumerate(names)},
               geom_cloud={k: geom_of(S[k]) for k in names},
               secs=float(time.time() - t0))
    return row


def run(mode, limit=0, draw_offset=0):
    from s25 import phys_lib as P
    path = CA_ROWS if mode == "ca" else CHAIN_ROWS
    if draw_offset:
        path = path.replace(".jsonl", f"_seed{draw_offset // N_DRAWS + 1}.jsonl")
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["pdb"])
                except Exception:
                    pass
    pdbs = P.targets()[:limit] if limit else P.targets()
    t0 = time.time()
    for n, pdb in enumerate(pdbs):
        if pdb in done:
            continue
        row = ca_row(pdb, draw_offset)[0] if mode == "ca" else chain_row(pdb)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(row) + "\n")
        orr = row["oracle_rmsd"] if mode == "ca" else row["oracle_rmsd_chain"]
        print(f"  [{mode} {n+1}/{len(pdbs)}] {pdb} n={row['n']} ORACLE rmsd PROD {orr['PROD']:.3f} sub0 {orr['sub0']:.3f} "
              f"circ_s0 {orr['circ_s0']:.3f} circ_best {orr['circ_best']:.3f} native {orr['NATIVE']:.3f} "
              f"{row['secs']:.1f}s (elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", path)


# ============================================================================ statistics
def wilson(k, n, z=1.959964):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [c - h, c + h]


def fold_ci(ind, folds, n_boot=4000, seed_parts=("s28C2",)):
    ind = np.asarray(ind, float); folds = np.asarray(folds)
    F = sorted(set(folds.tolist()))
    rng = SD.stable_rng(*seed_parts, salt=SALT)
    bs = np.array([np.concatenate([ind[folds == q] for q in rng.choice(F, len(F), replace=True)]).mean()
                   for _ in range(n_boot)])
    return [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]


def pairwise_logistic_fit(Delta, alpha, iters=50):
    """No intercept, L2 on scaled columns (training RMS; never centred), sign-augmented pairs."""
    X = np.asarray(Delta, float)
    mu = np.zeros(X.shape[1]); sd = np.sqrt((X ** 2).mean(0)) + 1e-9     # SCALE only: centring would erase the preference
    Z = (X - mu) / sd
    A = np.vstack([Z, -Z]); y = np.r_[np.zeros(len(Z)), np.ones(len(Z))]
    w = np.zeros(Z.shape[1])

    def loss(b):
        eta = A @ b
        return float((np.logaddexp(0.0, eta) - y * eta).sum() + 0.5 * alpha * (b @ b))
    cur = loss(w)
    for _ in range(iters):
        pr = 1.0 / (1.0 + np.exp(-np.clip(A @ w, -500, 500)))
        g = A.T @ (pr - y) + alpha * w
        Hm = (A * (pr * (1 - pr))[:, None]).T @ A + alpha * np.eye(len(w))
        step = np.linalg.solve(Hm + 1e-10 * np.eye(len(w)), g)
        t = 1.0; nxt = w - t * step; nl = loss(nxt)
        while nl > cur + 1e-12 and t > 1e-4:
            t *= 0.5; nxt = w - t * step; nl = loss(nxt)
        w, prev, cur = nxt, cur, nl
        if abs(prev - cur) < 1e-10 * max(1.0, abs(cur)):
            break
    return {"w": w, "mu": mu, "sd": sd}


def pairwise_decision(m, Delta):
    return ((np.asarray(Delta, float) - m["mu"]) / m["sd"]) @ m["w"]


def nested_pairwise(Delta, folds, alphas=np.logspace(-2, 4, 13), return_models=False):
    """Held-out decision values; negative = the ORACLE structure is preferred.  alpha chosen by
    inner leave-one-fold-out accuracy (ties -> larger alpha)."""
    Delta = np.asarray(Delta, float); folds = np.asarray(folds)
    F = sorted(set(folds.tolist()))
    dec = np.zeros(len(Delta)); chosen = {}; models = {}
    for f in F:
        tr = folds != f
        inner = [g for g in F if g != f]
        best, best_a = -np.inf, alphas[-1]
        for a in alphas:
            acc = []
            for g in inner:
                itr = tr & (folds != g); ite = folds == g
                m = pairwise_logistic_fit(Delta[itr], a)
                acc.append(float((pairwise_decision(m, Delta[ite]) < 0).mean()))
            sc = float(np.mean(acc))
            if sc >= best - 1e-12:
                best, best_a = sc, a
        m = pairwise_logistic_fit(Delta[tr], best_a)
        dec[folds == f] = pairwise_decision(m, Delta[folds == f])
        chosen[int(f)] = float(best_a)
        models[int(f)] = m
    if return_models:
        return dec, chosen, models
    return dec, chosen


def analyse(mode, n_null=200, n_max_null=500, print_all=True, tag=""):
    path = CA_ROWS if mode == "ca" else CHAIN_ROWS
    if tag:
        path = path.replace(".jsonl", f"_{tag}.jsonl")
    rows = [json.loads(l) for l in open(path, encoding="utf-8")]
    rows.sort(key=lambda r: r["pdb"])
    if mode == "chain":
        # S28-L34(e): one multiplicity null over ALL scorers seen on the projected chains, the 16
        # backbone scorers plus the 15 CA scorers re-evaluated on the chains (suffix "@chain")
        for r in rows:
            merged = dict(r["scores"])
            for nm, v in r.get("scores_ca_on_chain", {}).items():
                merged[nm + "@chain"] = v
            r["scores"] = merged
    pdbs = [r["pdb"] for r in rows]
    folds = np.array([int(r["fold"]) for r in rows])
    fail = np.array([p in I.FAIL18 for p in pdbs])
    n = len(rows)
    scorers = list(rows[0]["scores"].keys())
    key_r = "oracle_rmsd" if mode == "ca" else "oracle_rmsd_chain"
    out = {"mode": mode, "n": n, "scorers": scorers, "ladder_oracle_rmsd_mean": {}, "pref": {}, "ladder_rho": {},
           "contrast": {}, "strata": {}}
    # ladder RMSDs (ORACLE)
    for k in LADDER:
        out["ladder_oracle_rmsd_mean"][k] = float(np.mean([r[key_r][k] for r in rows]))
    for c in CONTROLS:
        vals = [np.mean([r[key_r][nm] for nm in r["names"] if nm.startswith(c + "[")]) for r in rows]
        out["ladder_oracle_rmsd_mean"][c] = float(np.mean(vals))
    # per scorer
    def S(r, nm, k):
        return r["scores"][nm][r["names"].index(k)]

    def below(a, b):
        """1 if a < b, 0.5 on an exact tie (S28-L34(a)), 0 otherwise."""
        return 1.0 if a < b else (0.5 if a == b else 0.0)
    # geometry beside every preference (S28-L34(b))
    gkey = "geom" if mode == "ca" else "geom_chain"
    if gkey in rows[0]:
        out["geom"] = {}
        for k in LADDER + CONTROLS:
            ks = [k] if k in LADDER else [nm2 for nm2 in rows[0]["names"] if nm2.startswith(k + "[")]
            out["geom"][k] = dict(bond=float(np.mean([r[gkey][kk]["bond"] for r in rows for kk in ks])),
                                  rg=float(np.mean([r[gkey][kk]["rg"] for r in rows for kk in ks])))
    pref_ind = {}
    for nm in scorers:
        rec = {}
        for k in LADDER[1:]:
            ind = np.array([below(S(r, nm, k), S(r, nm, "PROD")) for r in rows])
            ties = int(sum(1 for r in rows if S(r, nm, k) == S(r, nm, "PROD")))
            rec[k] = dict(pref=float(ind.mean()), wilson=wilson(float(ind.sum()), n), fold_ci=fold_ci(ind, folds, seed_parts=("pref", nm, k)),
                          fail18=float(ind[fail].mean()), fail18_k=float(ind[fail].sum()), other=float(ind[~fail].mean()), ties=ties)
            pref_ind[(nm, k)] = ind
        for c in CONTROLS:
            draws = [nm2 for nm2 in rows[0]["names"] if nm2.startswith(c + "[")]
            ind = np.array([np.mean([below(S(r, nm, d), S(r, nm, "PROD")) for d in draws if d in r["names"]]) for r in rows])
            ties = int(sum(1 for r in rows for d in draws if d in r["names"] and S(r, nm, d) == S(r, nm, "PROD")))
            rec[c] = dict(pref=float(ind.mean()), fold_ci=fold_ci(ind, folds, seed_parts=("pref", nm, c)),
                          fail18=float(ind[fail].mean()), other=float(ind[~fail].mean()), n_draws=len(draws), ties=ties)
            pref_ind[(nm, c)] = ind
        # HEAD-TO-HEAD (addendum 1): the ORACLE structure against the control itself, production absent
        for c in ("RAND_SIGNED", "GAUSS_MATCHED"):
            draws = [nm2 for nm2 in rows[0]["names"] if nm2.startswith(c + "[")]
            ind = np.array([np.mean([below(S(r, nm, "circ_best"), S(r, nm, d)) for d in draws if d in r["names"]]) for r in rows])
            rec["h2h_" + c] = dict(pref=float(ind.mean()), fold_ci=fold_ci(ind, folds, seed_parts=("h2h", nm, c)),
                                   fail18=float(ind[fail].mean()), other=float(ind[~fail].mean()))
        out["pref"][nm] = rec
        # ladder rho (per target, over the 5 ladder points)
        from scipy.stats import spearmanr
        rhos = []
        for r in rows:
            v = [S(r, nm, k) for k in LADDER]; d = [r[key_r][k] for k in LADDER]
            rho = spearmanr(v, d).correlation
            rhos.append(0.0 if not np.isfinite(rho) else float(rho))
        rhos = np.array(rhos)
        out["ladder_rho"][nm] = dict(mean=float(rhos.mean()), fold_ci=fold_ci(rhos, folds, seed_parts=("rho", nm)),
                                     frac_pos=float((rhos > 0).mean()), fail18=float(rhos[fail].mean()), other=float(rhos[~fail].mean()))
        # the recognition contrasts (paired indicators)
        cc = {}
        for c in CONTROLS:
            r_ = ST.compare(pref_ind[(nm, "circ_best")], pref_ind[(nm, c)], folds, names=pdbs,
                            label=f"{nm}: pref(ORACLE circ_best) - pref({c})")
            cc[c] = {k: v for k, v in r_.items() if k != "concentration"}
            r2_ = ST.compare(pref_ind[(nm, "circ_s0")], pref_ind[(nm, c)], folds, names=pdbs,
                             label=f"{nm}: pref(ORACLE circ_s0) - pref({c})")
            cc[c + "|circ_s0"] = {k: v for k, v in r2_.items() if k != "concentration"}
        out["contrast"][nm] = cc
    # falsifier check per scorer
    verdicts = {}
    for nm in scorers:
        p = out["pref"][nm]["circ_best"]; c = out["contrast"][nm]["RAND_SIGNED"]
        p0 = out["pref"][nm]["circ_s0"]; c0 = out["contrast"][nm]["RAND_SIGNED|circ_s0"]
        clause1 = p["fold_ci"][0] > 0.5
        clause2 = c["ci95_fold"][0] > 0.0
        anti = p["fold_ci"][1] < 0.5
        uninformative = p["ties"] > n / 2
        verdicts[nm] = dict(clause1_above_half=bool(clause1), clause2_beats_rand_signed=bool(clause2), anti_recognition=bool(anti),
                            uninformative_ties=bool(uninformative), fires=bool(clause1 and clause2 and not uninformative),
                            fires_on_circ_s0=bool(p0["fold_ci"][0] > 0.5 and c0["ci95_fold"][0] > 0.0 and not uninformative))
    out["verdicts"] = verdicts
    # max-over-scorers sign-flip null for pref(circ_best)
    rng = SD.stable_rng("maxnull", mode, salt=SALT)
    obs = {nm: out["pref"][nm]["circ_best"]["pref"] for nm in scorers}
    best = max(obs, key=obs.get)
    D = np.column_stack([pref_ind[(nm, "circ_best")] for nm in scorers])   # (n, S) indicators
    mx = np.empty(n_max_null)
    for t in range(n_max_null):
        eps = rng.random(n) < 0.5
        Dp = np.where(eps[:, None], D, 1 - D)                              # flip which structure is 'ORACLE'
        mx[t] = Dp.mean(0).max()
    out["max_null"] = dict(best_scorer=best, best_pref=obs[best], null_mean=float(mx.mean()), null_p95=float(np.percentile(mx, 95)),
                           p_max=float((mx >= obs[best]).mean()), n_scorers=len(scorers))
    # nested linear combination (pairwise logistic) with a sign-flip null
    Delta = np.column_stack([[S(r, nm, "circ_best") - S(r, nm, "PROD") for r in rows] for nm in scorers])
    ok = np.isfinite(Delta).all(1)
    dec, chosen, models = nested_pairwise(Delta[ok], folds[ok], return_models=True)
    acc = float((dec < 0).mean())
    # the same held-out rule applied to the controls' differences (control (i) for the combination)
    lc_ctrl = {}
    fk = folds[ok]
    for c in CONTROLS:
        draws = [nm2 for nm2 in rows[0]["names"] if nm2.startswith(c + "[")]
        fr = []
        for d in draws:
            Dc = np.column_stack([[S(r, nm, d) - S(r, nm, "PROD") for r in rows] for nm in scorers])[ok]
            dc = np.zeros(len(Dc))
            for f, m in models.items():
                dc[fk == f] = pairwise_decision(m, Dc[fk == f])
            fr.append((dc < 0).astype(float))
        ind_c = np.mean(fr, axis=0)
        r_ = ST.compare((dec < 0).astype(float), ind_c, fk, label=f"linear combination: pref(circ_best) - pref({c})")
        # head-to-head for the combination: the rule applied to score(circ_best) - score(control)
        h2h = []
        for d in draws:
            Dh = np.column_stack([[S(r, nm, "circ_best") - S(r, nm, d) for r in rows] for nm in scorers])[ok]
            dh = np.zeros(len(Dh))
            for f, m in models.items():
                dh[fk == f] = pairwise_decision(m, Dh[fk == f])
            h2h.append((dh < 0).astype(float))
        ind_h = np.mean(h2h, axis=0)
        lc_ctrl[c] = dict(pref=float(ind_c.mean()), fold_ci=fold_ci(ind_c, fk, seed_parts=("lc", c)),
                          contrast_effect=float(r_["effect"]), contrast_fold_ci=r_["ci95_fold"],
                          h2h_pref=float(ind_h.mean()), h2h_fold_ci=fold_ci(ind_h, fk, seed_parts=("lch2h", c)))
    null = np.empty(n_null)
    rng2 = SD.stable_rng("lcnull", mode, salt=SALT)
    for t in range(n_null):
        eps = np.where(rng2.random(int(ok.sum())) < 0.5, 1.0, -1.0)
        dp, _ = nested_pairwise(Delta[ok] * eps[:, None], folds[ok])
        null[t] = float((dp < 0).mean())
    ind_lc = np.zeros(n); ind_lc[ok] = (dec < 0).astype(float)
    out["linear_combination"] = dict(held_out_sign_acc=acc, n_used=int(ok.sum()), alphas=chosen, controls=lc_ctrl, null_mean=float(null.mean()),
                                     null_p95=float(np.percentile(null, 95)), p=float((null >= acc).mean()),
                                     fold_ci=fold_ci(ind_lc[ok], folds[ok], seed_parts=("lc",)),
                                     fail18=float(ind_lc[fail & ok].mean()), other=float(ind_lc[~fail & ok].mean()))
    pool_member_control(mode, out, rows, folds, fail, print_all=print_all)
    ST.save_atomic(os.path.join(RESULTS, f"s28_C2_{mode}{('_' + tag) if tag else ''}_summary.json"), out, module_file=__file__)
    if print_all:
        print(f"  ORACLE ladder mean RMSD ({mode}): " + "  ".join(f"{k} {v:.3f}" for k, v in out["ladder_oracle_rmsd_mean"].items()))
        if "geom" in out:
            print("  geometry (mean virtual bond A / Rg A): " + "  ".join(f"{k} {v['bond']:.2f}/{v['rg']:.2f}" for k, v in out["geom"].items()))
        print("  scorer            pref(circ_best) [fold CI]        pref(circ_s0) pref(sub0) pref(NATIVE) | RAND_SIGNED GAUSS_0.3 GAUSS_M | rho_ladder | vs RAND_SIGNED (fold CI) | FAIL18/108 | verdict")
        for nm in scorers:
            p = out["pref"][nm]; c = out["contrast"][nm]["RAND_SIGNED"]; v = verdicts[nm]
            print(f"  {nm:16s} {p['circ_best']['pref']:.3f} [{p['circ_best']['fold_ci'][0]:.3f},{p['circ_best']['fold_ci'][1]:.3f}]   "
                  f"{p['circ_s0']['pref']:.3f}  {p['sub0']['pref']:.3f}  {p['NATIVE']['pref']:.3f}   | "
                  f"{p['RAND_SIGNED']['pref']:.3f}  {p['GAUSS_0.3']['pref']:.3f}  {p['GAUSS_MATCHED']['pref']:.3f} | "
                  f"{out['ladder_rho'][nm]['mean']:+.3f} | {c['effect']:+.3f} [{c['ci95_fold'][0]:+.3f},{c['ci95_fold'][1]:+.3f}] | "
                  f"{p['circ_best']['fail18_k']:.0f}/18,{p['circ_best']['other']:.2f} | ties {p['circ_best']['ties']:3d} | "
                  f"{'FIRES' if v['fires'] else ('anti' if v['anti_recognition'] else 'no')}{' (s0 too)' if v['fires_on_circ_s0'] else ''}")
        m = out["max_null"]; lc = out["linear_combination"]
        print(f"  best single {m['best_scorer']} pref {m['best_pref']:.3f}; max-over-{m['n_scorers']} sign-flip null mean {m['null_mean']:.3f} p95 {m['null_p95']:.3f} p_max {m['p_max']:.3f}")
        print(f"  linear combination: held-out sign accuracy {lc['held_out_sign_acc']:.3f} [fold CI {lc['fold_ci'][0]:.3f},{lc['fold_ci'][1]:.3f}]  null mean {lc['null_mean']:.3f} p95 {lc['null_p95']:.3f} p {lc['p']:.3f}  FAIL18 {lc['fail18']:.2f} other {lc['other']:.2f}")
        print("    the same held-out rule on the controls: " + "  ".join(f"{c} {v['pref']:.3f} [{v['fold_ci'][0]:.3f},{v['fold_ci'][1]:.3f}] (contrast {v['contrast_effect']:+.3f} [{v['contrast_fold_ci'][0]:+.3f},{v['contrast_fold_ci'][1]:+.3f}]; head-to-head {v['h2h_pref']:.3f} [{v['h2h_fold_ci'][0]:.3f},{v['h2h_fold_ci'][1]:.3f}])" for c, v in lc["controls"].items()))
        print("  head-to-head (ORACLE circ_best scored below the control itself; production absent):")
        for nm in scorers:
            p = out["pref"][nm]
            print(f"    {nm:16s} vs RAND_SIGNED {p['h2h_RAND_SIGNED']['pref']:.3f} [{p['h2h_RAND_SIGNED']['fold_ci'][0]:.3f},{p['h2h_RAND_SIGNED']['fold_ci'][1]:.3f}]   vs GAUSS_MATCHED {p['h2h_GAUSS_MATCHED']['pref']:.3f} [{p['h2h_GAUSS_MATCHED']['fold_ci'][0]:.3f},{p['h2h_GAUSS_MATCHED']['fold_ci'][1]:.3f}]")
    return out


# ============================================================================ pool-member control (S28-L36)
POOL_CHANNEL = {"DIS": "DIS", "DIS_MEAN": "DIS_MEAN", "CONTACT_LL": "CONTACT_LL", "DISTPOT": "DISTPOT",
                "CONTACT": "CONTACT", "ENV": "ENV", "HP": "HP", "RG_LAW": "RG_LAW", "RG_UNIV": "RG_UNIV",
                "EXVOL": "EXVOL", "CAGEO": "CAGEO", "SS_MATCH": "SS_MATCH", "RAMA": "RAMA", "DSSPHB": "DSSPHB",
                "ELEC": "ELEC", "LEG": "LEG"}
POOL_CHANNEL.update({"LEG_" + t: "LEG_" + t for t in HL.LEG_TERMS})


def pool_member_control(mode, out, rows, folds, fail, print_all=True):
    """Lane D's S28-L36 control: a real pool member (a genuine protein trace with no information
    about the native).  For every scorer with a pool channel of the same name (the three
    pool-relative adapters are excluded), pct(X) = the fraction of the 500 pool members that score
    BETTER than X plus half the ties (X's percentile in its own pool; pool values from
    `s27/cache/<pdb>.npz`, the S27 pool rows the adapters reproduce).  pref(pool member vs PROD)
    = pct(PROD); the contrast pref(circ_best vs PROD) - pct(PROD) (paired, fold CI); head-to-head
    vs a pool member = 1 - pct(circ_best).  On the chain the CA scorers evaluated on the
    projected chains ("@chain") are compared with the pool members' CA-channel values."""
    from s27 import run_pool as RP
    cache = {}
    for r in rows:
        z = np.load(os.path.join(RP.CACHE, f"{r['pdb']}.npz"))
        cache[r["pdb"]] = {k: np.asarray(z[k], float) for k in z.files if k != "cost_ms"}
    scorers = out["scorers"]
    res = {}
    for nm in scorers:
        base = nm.replace("@chain", "")
        if base not in POOL_CHANNEL:
            continue
        ch = POOL_CHANNEL[base]
        pct = {k: [] for k in LADDER + CONTROLS}
        for r in rows:
            pool = cache[r["pdb"]][ch]
            for k in LADDER + CONTROLS:
                ks = [k] if k in LADDER else [nm2 for nm2 in r["names"] if nm2.startswith(k + "[")]
                vals = []
                for kk in ks:
                    x = r["scores"][nm][r["names"].index(kk)]
                    vals.append(float((pool < x).mean() + 0.5 * (pool == x).mean()))
                pct[k].append(float(np.mean(vals)))
        pct = {k: np.array(v) for k, v in pct.items()}
        ind_cb = np.array([out_pref_ind(out, rows, nm, "circ_best")])[0]
        ind_s0 = np.array([out_pref_ind(out, rows, nm, "circ_s0")])[0]
        c1 = ST.compare(ind_cb, pct["PROD"], folds, label=f"{nm}: pref(circ_best vs PROD) - pref(pool member vs PROD)")
        c2 = ST.compare(ind_s0, pct["PROD"], folds, label=f"{nm}: pref(circ_s0 vs PROD) - pref(pool member vs PROD)")
        res[nm] = dict(pool_channel=ch,
                       pct_median={k: float(np.median(v)) for k, v in pct.items()},
                       pct_mean={k: float(v.mean()) for k, v in pct.items()},
                       pref_pool_member_vs_prod=float(pct["PROD"].mean()),
                       contrast_circ_best={k: v for k, v in c1.items() if k != "concentration"},
                       contrast_circ_s0={k: v for k, v in c2.items() if k != "concentration"},
                       h2h_circ_best_beats_pool_member=float((1 - pct["circ_best"]).mean()),
                       h2h_native_beats_pool_member=float((1 - pct["NATIVE"]).mean()),
                       h2h_rand_signed_beats_pool_member=float((1 - pct["RAND_SIGNED"]).mean()),
                       fail18_contrast=float((ind_cb - pct["PROD"])[fail].mean()))
    out["pool_member_control"] = res
    if print_all:
        print("  pool-member control (S28-L36): pref(cb vs PROD) - pct(PROD) [fold CI], x MDE | h2h cb beats a pool member, NATIVE, RAND_SIGNED | pct(PROD) med, pct(cb), pct(NATIVE)")
        for nm, v in res.items():
            c = v["contrast_circ_best"]; c0 = v["contrast_circ_s0"]
            print(f"    {nm:20s} {c['effect']:+.3f} [{c['ci95_fold'][0]:+.3f},{c['ci95_fold'][1]:+.3f}] {c['effect_over_mde']:+.2f}x  (circ_s0 {c0['effect']:+.3f} {c0['effect_over_mde']:+.2f}x) | "
                  f"{v['h2h_circ_best_beats_pool_member']:.3f}  {v['h2h_native_beats_pool_member']:.3f}  {v['h2h_rand_signed_beats_pool_member']:.3f} | "
                  f"{v['pct_median']['PROD']:.2f}  {v['pct_mean']['circ_best']:.3f}  {v['pct_mean']['NATIVE']:.3f}")
    return res


def out_pref_ind(out, rows, nm, k):
    """Rebuild the (0, 0.5, 1) preference indicator of structure k over PROD for scorer nm."""
    def S(r, nm, k):
        return r["scores"][nm][r["names"].index(k)]
    return np.array([1.0 if S(r, nm, k) < S(r, nm, "PROD") else (0.5 if S(r, nm, k) == S(r, nm, "PROD") else 0.0) for r in rows])


# ============================================================================ selftest
def selftest():
    rng = np.random.default_rng(0)
    n, m = 12, 6
    base = np.cumsum(rng.normal(size=(n, 3)) * 2.5, 0)
    prod = base + rng.normal(scale=0.5, size=(n, 3))
    tgt = base + rng.normal(scale=3.0, size=(n, 3))
    C, s = scale_to_distance(prod, tgt, 1.25)
    assert abs(rmsd_between(C, prod) - 1.25) < 1e-6
    Cg, _ = gauss_perturb(prod, 0.3, "TEST", 0, "03")
    assert abs(rmsd_between(Cg, prod) - 0.3) < 1e-6
    # pairwise logistic: a planted preference is found, a random one is not
    folds = np.repeat(np.arange(5), 26)[:126]
    D = rng.normal(size=(126, 5)); D[:, 0] -= 1.5       # scorer 0 prefers ORACLE
    dec, _ = nested_pairwise(D, folds)
    acc = float((dec < 0).mean())
    Dr = rng.normal(size=(126, 5))
    dr, _ = nested_pairwise(Dr, folds)
    accr = float((dr < 0).mean())
    print(f"  pairwise logistic: planted {acc:.3f}, random {accr:.3f}")
    assert acc > 0.85 and accr < 0.7
    print("  s28_C2 selftest OK")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["selftest", "ca", "chain", "analyse_ca", "analyse_chain"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--draw-offset", type=int, default=0, help="0 = registered draws 0..3; 4 = the second control seed")
    ap.add_argument("--tag", default="", help="analyse: rows file suffix (e.g. seed2)")
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    if a.mode == "selftest":
        selftest()
    elif a.mode.startswith("analyse_"):
        analyse(a.mode.split("_", 1)[1], tag=a.tag)
    else:
        run(a.mode, a.limit, a.draw_offset)


if __name__ == "__main__":
    main()
