#!/usr/bin/env python
"""s29/s29_D_cost_audit.py -- lane D: THE COST-RMSD METER (S29 contract rule 19; brief S29D duty 2).

EVERY NUMBER THIS METER PRINTS IS ORACLE (it reads the native to score the ladder, to point the
gradient and to place the native in the pool).  The meter never tunes anything: it takes a cost
function f and reports, on the 126 dev targets with fold-clustered CIs, the four numbers the
charter asks of every objective before an endpoint run:

  (a) LADDER   Spearman rho(f, RMSD) per target over a ladder of structures, mean over targets.
               A cost whose reduction means RMSD reduction has rho -> +1 (f rises with RMSD).
               Three ladders are reported: S28 (C2's five rungs PROD, sub0, circ_s0, circ_best,
               NATIVE: the reproduction anchor, -0.402 on the chain / -0.182 at CA level for the
               shipped cost), CHARTER (the six rungs the S29 brief names: the native-free circuit
               optimum, a random signed combination, production, a Gaussian perturbation, the
               ORACLE circuit structure, the native) and FULL (all nine rungs).
  (b) COSINE   cos(-grad f at production, direction to the native), rigid-body components
               removed from both (S28-L23b machinery, `s27/s28_A2_local.py`); analytic for the
               shipped cost, central finite differences (batched, step --fd-h) for anything else;
               a random-direction reference (16 draws per target) beside it.
  (c) PCTILE   the native's percentile in its own pool under f: the fraction of the 500 pool
               members that score BETTER (lower) than the native posed in the frame
               (`s27/s28_A_objdiag.py`; 0.369 for the shipped surrogate S~).
  (d) PREF     pref(ORACLE circ_best vs PROD) = share of targets on which f scores the 0.29 A
               ORACLE structure below production (ties 0.5), beside the POOL-MEMBER CONTROL
               (S28-L36, `s27/s28_D_c2_poolmember.py`): pref(random pool member vs PROD) =
               PROD's percentile in the pool, and their paired contrast through `ST.compare`.
               A cost that recognises nativeness must prefer the ORACLE structure to production
               MORE often than it prefers an arbitrary real trace to production.

The cost function f is either a NAME from the S27 channel library (DIS = the shipped cost,
DIS_SURR = its linear-interpolation surrogate S~, or any of DIS_MEAN, CONTACT_LL, DISTPOT,
CONTACT, ENV, HP, RG_LAW, RG_UNIV, EXVOL, CAGEO, SS_MATCH, CONS_POOL, DMAP_CONS_POOL,
POOLGO_POOL at CA level; RAMA, DSSPHB, ELEC, TORS_CONS_POOL, LEG, LEG_<term> on the chain) or a
python CALLABLE given as `module:function`, with the contract

    f(W, ctx) -> (m,) float, lower is better, for W an (m, n, 3) stack of CA clouds and ctx a
    per-target context (pdb, seq, n, fold, k, the posed pool W (k, n, 3), PHI, PSI, the
    distogram dg, the DIS top-75 index `top`, the pool's DIS scores `dis`, the fold's `rama`
    counts, `universe()` (lazy), `frame`).  ctx.nat_ca and ctx.oracle_rr are NaN BY
    CONSTRUCTION: a cost that reads them emits NaN and the meter refuses it (the NaN-poison
    test is built in).  A callable that returns one float for a single (n, 3) cloud also works
    (the meter loops), but the batched form is what makes finite differences cheap.

Structures per target (CA level, lane A's frame, `s27/results/s28_A_structs/`): PROD (`prod`,
asserted equal to the deployed average), circ_best (`oracle_circ`, ORACLE best-of-5), NATIVE
(`oracle_aff500`, the native to 1e-7), circ_opt (`<pdb>_recog.npz :: circ_l1_i80`, lane A's
NATIVE-FREE recognition optimum at lam = 1, S28-L18b); regenerated exactly as C2 did and
CACHED in `s29/results/s29_D_ladder_structs/<pdb>.npz` with their assertions against the S28
artefacts: sub0 and circ_s0 (ORACLE), RAND_SIGNED[0], GAUSS_0.3[0], GAUSS_MATCHED[0] (controls,
native-free directions, ORACLE scale).  `build-cache --chain` also stores every rung's
production projection (CA, phi, psi) so any f can be metered on the BUILT CHAIN, the
reporting basis.  `--basis chain-s28rows` reads C2's stored chain scores instead (the shipped
cost's -0.402 / 0.071 come from there; only for the 31 S28 scorer names).

    python s29/s29_D_cost_audit.py build-cache [--limit N] [--chain] [--rebuild]
    python s29/s29_D_cost_audit.py meter --f DIS [--basis ca|chain|chain-s28rows] [--limit N]
                                         [--fd-h 1e-3] [--tag NAME] [--out PATH]
    python s29/s29_D_cost_audit.py meter --f s29.s29_X_cost:my_cost --basis ca
    python s29/s29_D_cost_audit.py selftest

Tests: `tests/test_s29_D.py` (the S28 values for the shipped cost reproduce to the third
decimal through this code; the callable path; the NaN poison; finite differences agree with the
analytic gradient).  Results: `s29/results/s29_D_cost_audit_<tag>_<basis>.json`.
"""
from __future__ import annotations

import argparse
import importlib
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
from s27 import s28_A_amp as A             # noqa: E402
from s27 import s28_A2_local as A2         # noqa: E402
from s27 import s28_C2_recog_audit as C2   # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(RESULTS, "s29_D_ladder_structs")
S28_RESULTS = os.path.join(ROOT, "s27", "results")
S27_CACHE = os.path.join(ROOT, "s27", "cache")
SALT = "s29D"

#: rung names, in the order stored; ORACLE rungs are chosen against the native
RUNGS = ["circ_opt", "RAND_SIGNED[0]", "PROD", "GAUSS_MATCHED[0]", "GAUSS_0.3[0]",
         "sub0", "circ_s0", "circ_best", "NATIVE"]
ORACLE_RUNGS = {"sub0", "circ_s0", "circ_best", "NATIVE"}
LADDERS = {
    "S28": ["PROD", "sub0", "circ_s0", "circ_best", "NATIVE"],                        # C2's LADDER
    "CHARTER": ["circ_opt", "RAND_SIGNED[0]", "PROD", "GAUSS_MATCHED[0]", "circ_best", "NATIVE"],
    "FULL": list(RUNGS),
}
CIRC_OPT_KEY = "circ_l1_i80"          # lane A's native-free recognition optimum, lam = 1, 80 iters
N_COS_REF = 16
CA_NAMES = ["DIS", "DIS_SURR"] + C2.CA_SCORERS[1:]          # DIS listed once
CHAIN_NAMES = list(C2.CHAIN_SCORERS)
SHRINK_E = 0.3            # the descent displacement at which the shrink signature is read (contract addendum 20)
FD_H_DEFAULT = 1e-3       # h = 1e-3 A: cosine within 3e-4 of the analytic on the piecewise-linear S~; step-function costs need --fd-h 0.5 (A2's smoothed FD)


def _safe(k):
    return k.replace("[", "_").replace("]", "")


# ============================================================================ per-target context
class Ctx:
    """What a cost may read.  NATIVE-FREE: nat_ca / oracle_rr are NaN by construction."""

    def __init__(self, pdb, cand, frame, dg, top, dis):
        self.pdb, self.seq, self.n, self.fold, self.k = pdb, cand.seq, int(cand.n), int(cand.fold), int(cand.k)
        self.W = np.asarray(frame.Wp, float)               # the pool, posed in the frame
        self.PHI = np.asarray(cand.PHI, float); self.PSI = np.asarray(cand.PSI, float)
        self.dg = dg; self.top = np.asarray(top, int); self.dis = np.asarray(dis, float)
        self.frame = frame
        self.rama = C2.RAMA[int(cand.fold)]
        self.nat_ca = np.full((self.n, 3), np.nan)          # POISONED
        self.oracle_rr = np.full(self.k, np.nan)            # POISONED
        self._u = None
        self._cand = cand

    def universe(self):
        if self._u is None:
            self._u = I.load_univ(self.pdb)
        return self._u

    def cand_poisoned(self):
        """A pool object with the natives NaN, for the S27 adapters (which never read them)."""
        c = SimpleNamespace(pdb=self.pdb, seq=self.seq, n=self.n, fold=self.fold, k=self.k,
                            W=np.asarray(self._cand.W, float), PHI=self.PHI, PSI=self.PSI,
                            nat_ca=self.nat_ca, oracle_rr=self.oracle_rr)
        return c


def load_target(pdb):
    cand, dis, top, dg = A.load_pool(pdb)
    frame = A.Frame(cand.W, top)
    sur = A.Surrogate(dg, cand.n)
    return cand, dis, top, dg, frame, sur


# ============================================================================ the cost wrapper
class Cost:
    """A named S27 channel or a callable, with values on a stack and a gradient at one cloud."""

    def __init__(self, spec, fd_h=FD_H_DEFAULT):
        self.spec = str(spec); self.fd_h = float(fd_h)
        self.fn = None
        if ":" in self.spec:
            mod, fn = self.spec.split(":", 1)
            self.fn = getattr(importlib.import_module(mod), fn)
            self.kind = "callable"; self.name = self.spec.replace(":", ".")
        elif self.spec in CA_NAMES:
            self.kind = "ca"; self.name = self.spec
        elif self.spec in CHAIN_NAMES:
            self.kind = "chain"; self.name = self.spec
        else:
            raise KeyError("unknown cost %r; names: %s | %s, or module:function" % (spec, CA_NAMES, CHAIN_NAMES))
        self.analytic = self.spec in ("DIS", "DIS_SURR", "RG_LAW", "EXVOL")
        self._mode = None            # callables: "batched" or "loop", decided on the first call
        self.errors = []             # every exception a callable raised, kept for the refusal message

    # ---- values -------------------------------------------------------------------------
    def values(self, W, ctx, sur=None, PHI=None, PSI=None):
        """(m,) values of f on the stack W (m, n, 3); PHI/PSI (m, n) for the chain scorers."""
        W = np.asarray(W, float)
        if W.ndim == 2:
            W = W[None]
        if self.kind == "callable":
            return self._call(W, ctx)
        if self.spec == "DIS":
            i, j = I.pair_index(ctx.n)
            D = I.pair_dists(W, i, j)
            return np.asarray(I.shipped_score(ctx.dg, D.astype(np.float32).astype(float)), float)
        if self.spec == "DIS_SURR":
            return np.array([sur.value_grad(w)[0] for w in W])
        if self.kind == "ca":
            sc = C2.ca_scores(ctx.cand_poisoned(), W, ctx.universe(), ctx.dg)
            return np.asarray(sc[self.spec], float)
        if PHI is None:
            raise ValueError("%s is a backbone scorer: it needs projected chains (--basis chain)" % self.spec)
        sc = C2.chain_scores(ctx.cand_poisoned(), W, PHI, PSI, ctx.universe(), ctx.dg)
        return np.asarray(sc[self.spec], float)

    def _call(self, W, ctx):
        """Call the callable, batched if it takes a stack and per-structure if it does not.

        A cost that reads ctx.nat_ca (NaN by construction) either returns NaN or RAISES: both
        are a refusal, so every exception is recorded in `self.errors` and returned as NaN,
        and `run_meter` refuses the cost naming the failure.  Nothing is silently swallowed.
        """
        errs = []
        for mode in (("batched", "loop") if self._mode is None else (self._mode,)):
            try:
                if mode == "batched":
                    v = np.asarray(self.fn(W, ctx), float).ravel()
                    if v.shape != (len(W),):
                        raise ValueError("batched call returned shape %s for a stack of %d" % (v.shape, len(W)))
                else:
                    v = np.array([float(np.asarray(self.fn(w, ctx), float).ravel()[0]) for w in W])
                self._mode = mode
                return v
            except Exception as e:                       # noqa: BLE001 -- recorded, never hidden
                errs.append("%s call: %s: %s" % (mode, type(e).__name__, e))
        self.errors.extend(errs)
        return np.full(len(W), np.nan)

    # ---- gradient at one cloud ----------------------------------------------------------
    def grad(self, C0, ctx, sur=None):
        """(value, dF/dC (n, 3), method).  Analytic where S28 had it, else central differences."""
        C0 = np.asarray(C0, float)
        if self.spec in ("DIS", "DIS_SURR"):
            v, g = sur.value_grad(C0)
            return float(v), g, "analytic (surrogate S~ gradient, S28-L23b)"
        if self.spec == "RG_LAW":
            v, g = A2.grad_rg_law(C0, ctx.n)
            return float(v), g, "analytic"
        if self.spec == "EXVOL":
            v, g = A2.grad_exvol(C0)
            return float(v), g, "analytic"
        n = len(C0); h = self.fd_h
        stack = [C0]
        for k in range(3 * n):
            e = np.zeros(3 * n); e[k] = h
            stack.append(C0 + e.reshape(n, 3)); stack.append(C0 - e.reshape(n, 3))
        vals = self.values(np.stack(stack), ctx, sur=sur)
        g = ((vals[1::2] - vals[2::2]) / (2 * h)).reshape(n, 3)
        return float(vals[0]), g, "central finite differences, h = %g A" % h


# ============================================================================ the ladder cache
def cache_path(pdb):
    return os.path.join(CACHE, f"{pdb}.npz")


def build_rungs(pdb, cand, frame, dg, top, with_chain=False, check=True):
    """Regenerate the rungs exactly as C2's `build_structures` (ORACLE for sub0 / circ_s0 and for
    every scale), assert them against the S28 artefacts, and return the dict of (n, 3) clouds
    plus their ORACLE RMSDs (and the projected chains if asked)."""
    from core import quantum as Q
    from s24 import d_harness as H
    z = np.load(os.path.join(A.STRUCTS, f"{pdb}.npz"))
    S = {"PROD": np.asarray(z["prod"], float), "circ_best": np.asarray(z["oracle_circ"], float),
         "NATIVE": np.asarray(z["oracle_aff500"], float)}
    zr = np.load(os.path.join(A.STRUCTS, f"{pdb}_recog.npz"))
    S["circ_opt"] = np.asarray(zr[CIRC_OPT_KEY], float)
    Cd, _ = H.readout_uniform(cand, top)
    assert np.abs(S["PROD"] - Cd).max() < 1e-10, f"{pdb}: PROD differs from the deployed average"
    arow = C2.load_A(pdb)[1]
    circ = Q.StatevectorCircuit(A.N_QUBITS, A.LAYERS)
    oc = A.oracle_circuit_ceiling(circ, frame, cand.nat_ca, starts=1)                 # ORACLE
    S["circ_s0"] = np.asarray(oc["C"], float)
    dev_s0 = abs(oc["rmsd_cloud"] - float(arow["arms"]["oracle_circ"]["per_start"][0]))
    rs, Cs, _ = A.oracle_subspace_ls(frame, cand.nat_ca, A.subspace_matrix(pdb, 0, frame.k))  # ORACLE
    S["sub0"] = np.asarray(Cs, float)
    dev_sub = abs(rs - float(arow["arms"]["oracle_sub"]["per_sub"][0]))
    if check:
        assert dev_s0 < 1e-6, f"{pdb}: circ_s0 {oc['rmsd_cloud']} != lane A's per_start[0] (dev {dev_s0})"
        assert dev_sub < 1e-6, f"{pdb}: sub0 deviates from lane A's per_sub[0] by {dev_sub}"
    d_native = C2.rmsd_between(S["circ_best"], cand.nat_ca)                                # ORACLE scale
    d_prod = C2.rmsd_between(S["circ_best"], S["PROD"])
    S["RAND_SIGNED[0]"], meta_rs = C2.rand_signed(frame, S["PROD"], d_prod, pdb, 0)
    S["GAUSS_0.3[0]"], _ = C2.gauss_perturb(S["PROD"], d_native, pdb, 0, "03")
    S["GAUSS_MATCHED[0]"], _ = C2.gauss_perturb(S["PROD"], d_prod, pdb, 0, "m")
    orr = {k: C2.rmsd_between(S[k], cand.nat_ca) for k in RUNGS}                             # ORACLE labels
    out = {"S": S, "oracle_rmsd": orr, "meta": dict(d_circ_best_native=d_native, d_circ_best_prod=d_prod,
                                                    dev_circ_s0=dev_s0, dev_sub0=dev_sub, **{"rand_signed_" + k: v for k, v in meta_rs.items()})}
    # the regenerated controls against C2's stored shipped-DIS scores (a regeneration check)
    crow = _c2_row(pdb, "ca")
    if crow is not None:
        i, j = I.pair_index(cand.n)
        devs = {}
        for k in ("RAND_SIGNED[0]", "GAUSS_0.3[0]", "GAUSS_MATCHED[0]", "sub0", "circ_s0"):
            D = I.pair_dists(S[k][None], i, j)
            mine = float(I.shipped_score(dg, D.astype(np.float32).astype(float))[0])
            devs[k] = abs(mine - float(crow["scores"]["DIS"][crow["names"].index(k)]))
        out["meta"]["dev_vs_c2_dis"] = devs
        if check:
            assert max(devs.values()) < 1e-5, f"{pdb}: regenerated rungs disagree with C2's stored DIS: {devs}"
    if with_chain:
        ch = {}
        for k in RUNGS:
            pr = I.project(S[k], cand.seq, cand.fold)
            ch[k] = dict(ca=np.asarray(pr["ca"], float), phi=np.asarray(pr["phi"], float), psi=np.asarray(pr["psi"], float))
        out["chain"] = ch
        out["oracle_rmsd_chain"] = {k: C2.rmsd_between(ch[k]["ca"], cand.nat_ca) for k in RUNGS}      # ORACLE
    return out


def save_rungs(pdb, cand, b):
    os.makedirs(CACHE, exist_ok=True)
    f = cache_path(pdb)
    old = {}
    if os.path.exists(f):
        with np.load(f, allow_pickle=False) as zz:
            old = {k: np.array(zz[k]) for k in zz.files}
    d = dict(old)
    for k in RUNGS:
        d[_safe(k)] = b["S"][k]
        d["rmsd_" + _safe(k)] = np.array(float(b["oracle_rmsd"][k]))
    if "chain" in b:
        for k in RUNGS:
            for q in ("ca", "phi", "psi"):
                d[f"chain_{q}_{_safe(k)}"] = b["chain"][k][q]
            d["rmsdchain_" + _safe(k)] = np.array(float(b["oracle_rmsd_chain"][k]))
    d["meta"] = np.array(json.dumps(b["meta"]))
    d["seq"] = np.array(cand.seq); d["fold"] = np.array(int(cand.fold)); d["n"] = np.array(int(cand.n))
    tmp = f + f".tmp{os.getpid()}.npz"
    np.savez_compressed(tmp, **d)
    A.replace_retry(tmp, f)


def load_rungs(pdb, need_chain=False):
    f = cache_path(pdb)
    if not os.path.exists(f):
        return None
    with np.load(f, allow_pickle=False) as z:
        keys = set(z.files)
        if any(_safe(k) not in keys for k in RUNGS):
            return None
        if need_chain and any(f"chain_ca_{_safe(k)}" not in keys for k in RUNGS):
            return None
        out = {"S": {k: np.array(z[_safe(k)]) for k in RUNGS},
               "oracle_rmsd": {k: float(z["rmsd_" + _safe(k)]) for k in RUNGS},
               "meta": json.loads(str(z["meta"]))}
        if f"chain_ca_{_safe(RUNGS[0])}" in keys:
            out["chain"] = {k: {q: np.array(z[f"chain_{q}_{_safe(k)}"]) for q in ("ca", "phi", "psi")} for k in RUNGS}
            out["oracle_rmsd_chain"] = {k: float(z["rmsdchain_" + _safe(k)]) for k in RUNGS}
    return out


_C2_ROWS = {}


def _c2_row(pdb, mode):
    if mode not in _C2_ROWS:
        p = C2.CA_ROWS if mode == "ca" else C2.CHAIN_ROWS
        _C2_ROWS[mode] = {}
        if os.path.exists(p):
            with open(p, encoding="utf-8") as fh:
                for line in fh:
                    if line.strip():
                        r = json.loads(line); _C2_ROWS[mode][r["pdb"]] = r
    return _C2_ROWS[mode].get(pdb)


def build_cache(pdbs, with_chain=False, rebuild=False):
    t0 = time.time()
    for q, pdb in enumerate(pdbs):
        have = load_rungs(pdb, need_chain=with_chain)
        if have is not None and not rebuild:
            continue
        t1 = time.time()
        cand, dis, top, dg, frame, sur = load_target(pdb)
        b = build_rungs(pdb, cand, frame, dg, top, with_chain=with_chain)
        save_rungs(pdb, cand, b)
        o = b["oracle_rmsd"]
        msg = " ".join(f"{k} {o[k]:.3f}" for k in RUNGS)
        if with_chain:
            msg += " | chain " + " ".join(f"{k} {b['oracle_rmsd_chain'][k]:.3f}" for k in RUNGS)
        print(f"  [cache {q+1}/{len(pdbs)}] {pdb} n={cand.n} ORACLE rmsd: {msg}  {time.time()-t1:.1f}s "
              f"(elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("cache:", CACHE)


# ============================================================================ per-target meter
def spearman(x, y):
    """Spearman rho; NaN when either input is constant (C2 mapped that to 0.0; reported both ways)."""
    import warnings
    from scipy.stats import spearmanr
    x = np.asarray(x, float); y = np.asarray(y, float)
    if np.ptp(x) == 0 or np.ptp(y) == 0 or not (np.isfinite(x).all() and np.isfinite(y).all()):
        return float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = spearmanr(x, y).correlation
    return float(r) if np.isfinite(r) else float("nan")


def pct_in_pool(pool, v):
    """Percentile of value v in the pool: share scoring BETTER (lower) plus half the ties."""
    pool = np.asarray(pool, float)
    return float((pool < v).mean() + 0.5 * (pool == v).mean())


def meter_target(pdb, cost: Cost, basis="ca", need_grad=True):
    """Everything the four numbers need for one target.  ORACLE throughout (labelled)."""
    t0 = time.time()
    cand, dis, top, dg, frame, sur = load_target(pdb)
    ctx = Ctx(pdb, cand, frame, dg, top, dis)
    rungs = load_rungs(pdb, need_chain=(basis == "chain"))
    if rungs is None:
        raise FileNotFoundError(f"{pdb}: no ladder cache for basis {basis}; run `build-cache`"
                                + (" --chain" if basis == "chain" else ""))
    row = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), fail18=bool(pdb in I.FAIL18), basis=basis,
               cost=cost.name, oracle=True)
    S = rungs["S"]
    # ---- f on the rungs
    if basis == "ca":
        W = np.stack([S[k] for k in RUNGS])
        fv = cost.values(W, ctx, sur=sur)
        orr = rungs["oracle_rmsd"]
    else:
        W = np.stack([rungs["chain"][k]["ca"] for k in RUNGS])
        PHI = np.stack([rungs["chain"][k]["phi"] for k in RUNGS]); PSI = np.stack([rungs["chain"][k]["psi"] for k in RUNGS])
        fv = cost.values(W, ctx, sur=sur, PHI=PHI, PSI=PSI)
        orr = rungs["oracle_rmsd_chain"]
    if not np.isfinite(fv).all():
        row["nan_on_rungs"] = [k for k, v in zip(RUNGS, fv) if not np.isfinite(v)]
    row["f"] = {k: float(v) for k, v in zip(RUNGS, fv)}
    row["rmsd"] = {k: float(orr[k]) for k in RUNGS}                                    # ORACLE
    # ---- (a) ladders
    row["rho"] = {}
    for nm, L in LADDERS.items():
        row["rho"][nm] = spearman([row["f"][k] for k in L], [row["rmsd"][k] for k in L])
    # ---- f on the pool (real traces; CA level) and the native posed in the frame
    if cost.kind in ("ca", "callable") or cost.spec in ("DIS", "DIS_SURR"):
        if cost.kind == "ca" and cost.spec in C2.POOL_CHANNEL and os.path.exists(os.path.join(S27_CACHE, f"{pdb}.npz")):
            with np.load(os.path.join(S27_CACHE, f"{pdb}.npz")) as zc:
                pool = np.asarray(zc[C2.POOL_CHANNEL[cost.spec]], float)   # the S27 pool rows (exact)
            row["pool_source"] = "s27/cache"
        else:
            pool = cost.values(frame.Wp, ctx, sur=sur)
            row["pool_source"] = "evaluated on the posed pool"
        natp = I.superpose_batch(np.asarray(cand.nat_ca, float)[None], frame.ref)[0]      # ORACLE
        f_nat = float(cost.values(natp[None], ctx, sur=sur)[0])
        row["f_native_posed"] = f_nat
        ok = np.isfinite(pool)
        row["pool_mean"] = float(pool[ok].mean()) if ok.any() else float("nan")
        row["pool_min"] = float(pool[ok].min()) if ok.any() else float("nan")
        row["native_pctile_strict"] = float((pool < f_nat).mean())                # objdiag's definition
        row["native_pctile"] = pct_in_pool(pool, f_nat)
        row["pct"] = {k: pct_in_pool(pool, row["f"][k]) for k in RUNGS}
        row["pool_finite"] = int(np.isfinite(pool).sum())
    else:
        row["pool_source"] = "none (backbone scorer: no CA pool channel)"
        row["pct"] = {}
    # ---- (d) preferences vs PROD (ties 0.5)
    fp = row["f"]["PROD"]
    row["pref"] = {k: (1.0 if row["f"][k] < fp else (0.5 if row["f"][k] == fp else 0.0)) for k in RUNGS if k != "PROD"}
    # ---- (b) the gradient cosine at production (CA level only; the projection is not differentiable)
    if basis == "ca" and need_grad:
        C0 = S["PROD"]
        val, g, method = cost.grad(C0, ctx, sur=sur)
        gp = A2.remove_rigid(g, C0)
        u = A2.remove_rigid(A2.oracle_direction(C0, cand.nat_ca), C0)                  # ORACLE
        row["cos"] = A2.cosine(-gp, u)
        row["grad_rms"] = A2.rms(gp); row["u_rms"] = A2.rms(u); row["grad_method"] = method
        row["grad_zero_frac"] = float((np.abs(g) < 1e-300).mean())
        #: CONTRACT ADDENDUM 20 (the shrink rule, lane T's theorem 2, S29-L7): a positive cosine
        #: is purchasable with no information by shrinking the target map toward typicality,
        #: which CONTRACTS the emitted structure.  The descent step's geometry is the signature:
        #: descend 0.3 A along -grad f and report what happened to the virtual bond and Rg.
        if A2.rms(gp) > 1e-15:
            Cs = C0 - SHRINK_E * gp / A2.rms(gp)
            gp0, gps = A.struct_diag(C0), A.struct_diag(Cs)
            row["shrink"] = dict(e=SHRINK_E, bond_prod=gp0["bond"], rg_prod=gp0["rg"],
                                 bond_step=gps["bond"], rg_step=gps["rg"],
                                 bond_ratio=gps["bond"] / max(gp0["bond"], 1e-12),
                                 rg_ratio=gps["rg"] / max(gp0["rg"], 1e-12))
        rng = SD.stable_rng(pdb, "s28A2_cosref", salt=A2.SALT)                      # A2's own reference draws
        row["cos_random_ref"] = [A2.cosine(A2.remove_rigid(rng.normal(size=C0.shape), C0), u) for _ in range(N_COS_REF)]
    row["secs"] = time.time() - t0
    return row


def meter_target_from_s28_rows(pdb, cost: Cost):
    """The shipped cost's chain numbers straight from C2's stored chain rows (S28-L48 artefact)."""
    r = _c2_row(pdb, "chain")
    if r is None:
        raise FileNotFoundError("no C2 chain row for " + pdb)
    base = cost.spec
    key = "scores_ca_on_chain" if base in C2.CA_SCORERS else "scores"
    if base not in r[key]:
        raise KeyError(f"{base} is not among C2's stored chain scorers")
    names = r["names"]
    f = {k: float(r[key][base][names.index(k)]) for k in names}
    row = dict(pdb=pdb, n=int(r["n"]), fold=int(r["fold"]), fail18=bool(pdb in I.FAIL18), basis="chain-s28rows",
               cost=cost.name, oracle=True, f=f, rmsd={k: float(r["oracle_rmsd_chain"][k]) for k in names})
    row["rho"] = {}
    for nm, L in LADDERS.items():
        L2 = [k for k in L if k in names]
        row["rho"][nm] = spearman([f[k] for k in L2], [row["rmsd"][k] for k in L2]) if len(L2) >= 3 else float("nan")
        if len(L2) < len(L):
            row["rho"][nm + "_note"] = "rungs absent in C2's chain rows: " + ",".join(k for k in L if k not in names)
    fp = f["PROD"]
    row["pref"] = {k: (1.0 if f[k] < fp else (0.5 if f[k] == fp else 0.0)) for k in names if k != "PROD"}
    if base in C2.POOL_CHANNEL and os.path.exists(os.path.join(S27_CACHE, f"{pdb}.npz")):
        with np.load(os.path.join(S27_CACHE, f"{pdb}.npz")) as zc:
            pool = np.asarray(zc[C2.POOL_CHANNEL[base]], float)
        row["pct"] = {k: pct_in_pool(pool, f[k]) for k in names}
        row["pool_source"] = "s27/cache (real pool traces) vs PROJECTED structures: cross-basis, as S28-L48"
    else:
        row["pct"] = {}
    return row


# ============================================================================ the summary
def _summ(x, folds, fail, label):
    """Mean, median, SE and both CIs of a per-target quantity, NaNs dropped and counted.

    With fewer than 3 finite targets (a cost whose gradient is zero everywhere gives NONE) every
    key is still present and the CIs are None, so a caller can print the row without a crash --
    the defect `tests/test_s29_D.py::test_all_nan_axis_does_not_crash_the_render` pins.
    """
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    out = dict(label=label, n=int(ok.sum()), n_nan=int((~ok).sum()))
    if ok.sum() < 3:
        out.update(mean=float("nan"), median=float("nan"), se=float("nan"), ci95_iid=None, ci95_fold=None,
                   per_fold=None, folds_same_sign=0, n_pos=0, n_neg=0,
                   fail18_mean=float("nan"), other108_mean=float("nan"))
        return out
    r = ST.compare(x[ok], np.zeros(int(ok.sum())), folds[ok], label=label, seed_parts=("s29D",))
    out.update(mean=r["effect"], median=r["median_effect"], se=r["se"], ci95_iid=r["ci95_iid"], ci95_fold=r["ci95_fold"],
               per_fold=r["per_fold"], folds_same_sign=r["folds_same_sign"], n_pos=int((x[ok] > 0).sum()),
               n_neg=int((x[ok] < 0).sum()), fail18_mean=float(np.nanmean(x[fail])) if fail.any() else float("nan"),
               other108_mean=float(np.nanmean(x[~fail])) if (~fail).any() else float("nan"))
    return out


def summarise(rows, cost: Cost, basis):
    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    fail = np.array([r["fail18"] for r in rows])
    names = [k for k in RUNGS if all(k in r["f"] for r in rows)]
    out = dict(cost=cost.name, spec=cost.spec, kind=cost.kind, basis=basis, n=len(rows), oracle=True,
               rungs=names, ladders={k: [x for x in v if x in names] for k, v in LADDERS.items()},
               mean_f={k: float(np.nanmean([r["f"][k] for r in rows])) for k in names},
               mean_rmsd={k: float(np.mean([r["rmsd"][k] for r in rows])) for k in names})
    # (a)
    out["ladder_rho"] = {nm: _summ([r["rho"].get(nm, np.nan) for r in rows], folds, fail, f"ladder rho ({nm}) vs 0")
                         for nm in LADDERS}
    for nm in LADDERS:
        out["ladder_rho"][nm]["n_undefined_set_to_0_as_C2"] = int(sum(1 for r in rows if not np.isfinite(r["rho"].get(nm, np.nan))))
        out["ladder_rho"][nm]["mean_C2_convention"] = float(np.mean([0.0 if not np.isfinite(r["rho"].get(nm, np.nan)) else r["rho"][nm] for r in rows]))
    # (b)
    if all("cos" in r for r in rows):
        out["cosine"] = _summ([r["cos"] for r in rows], folds, fail, "cos(-grad f, native direction) vs 0")
        out["cosine"]["n_defined"] = int(np.isfinite([r["cos"] for r in rows]).sum())
        sh = [r["shrink"] for r in rows if "shrink" in r]
        out["cosine"]["shrink_signature"] = (dict(
            n=len(sh), e=SHRINK_E,
            bond_ratio_mean=float(np.mean([x["bond_ratio"] for x in sh])),
            rg_ratio_mean=float(np.mean([x["rg_ratio"] for x in sh])),
            n_contracting=int(sum(1 for x in sh if x["rg_ratio"] < 1.0)),
            bond_prod_mean=float(np.mean([x["bond_prod"] for x in sh])),
            rg_prod_mean=float(np.mean([x["rg_prod"] for x in sh])),
            note="contract addendum 20 / S29-L7: a cosine gain bought by shrinking the target map CONTRACTS "
                 "the emitted structure; read this beside the cosine and beside the native percentile, which "
                 "the shrink moves the wrong way. Production's own trace is already 22% contracted (S23 L1).")
            if sh else None)
        out["cosine"]["random_ref_mean_abs"] = float(np.mean([np.mean(np.abs(r["cos_random_ref"])) for r in rows]))
        out["cosine"]["random_ref_p95_abs"] = float(np.percentile(np.concatenate([np.abs(r["cos_random_ref"]) for r in rows]), 95))
        out["cosine"]["grad_method"] = rows[0]["grad_method"]
        out["cosine"]["grad_rms_mean"] = float(np.mean([r["grad_rms"] for r in rows]))
        out["cosine"]["grad_zero_frac_mean"] = float(np.mean([r["grad_zero_frac"] for r in rows]))
        cs = np.array([r["cos"] for r in rows], float); pr = np.array([r["rmsd"]["PROD"] for r in rows], float)
        ok = np.isfinite(cs)
        out["cosine"]["spearman_cos_vs_prod_rmsd"] = spearman(cs[ok], pr[ok]) if ok.sum() >= 3 else float("nan")
    else:
        out["cosine"] = None
    # (c)
    if all("native_pctile" in r for r in rows):
        out["native_pctile"] = _summ([r["native_pctile"] for r in rows], folds, fail, "native percentile in pool")
        out["native_pctile"]["mean_strict_objdiag"] = float(np.mean([r["native_pctile_strict"] for r in rows]))
        out["native_pctile"]["n_prod_below_native"] = int(sum(1 for r in rows if r["f"]["PROD"] < r["f_native_posed"]))
        out["native_pctile"]["n_prod_below_pool_min"] = int(sum(1 for r in rows if r["f"]["PROD"] < r["pool_min"]))
        out["native_pctile"]["mean_f_native_posed"] = float(np.mean([r["f_native_posed"] for r in rows]))
        out["native_pctile"]["mean_pool_min"] = float(np.mean([r["pool_min"] for r in rows]))
        out["native_pctile"]["mean_pool_mean"] = float(np.mean([r["pool_mean"] for r in rows]))
    else:
        out["native_pctile"] = None
    # (d)
    pref = {}
    for k in names:
        if k == "PROD":
            continue
        ind = np.array([r["pref"][k] for r in rows])
        s = _summ(ind, folds, fail, f"pref({k} vs PROD)")
        s["ties"] = int((ind == 0.5).sum())
        pref[k] = s
    out["pref"] = pref
    have_pct = all(r.get("pct") for r in rows)
    if have_pct:
        pct_prod = np.array([r["pct"]["PROD"] for r in rows])
        ind_cb = np.array([r["pref"]["circ_best"] for r in rows])
        c = ST.compare(ind_cb, pct_prod, folds, names=pdbs, label=f"{cost.name}: pref(ORACLE circ_best vs PROD) - pref(pool member vs PROD)",
                       seed_parts=("s29D",))
        out["pool_member_control"] = dict(
            pref_circ_best_vs_prod=float(ind_cb.mean()), pref_pool_member_vs_prod=float(pct_prod.mean()),
            pct_prod_median=float(np.median(pct_prod)),
            contrast={k: v for k, v in c.items() if k != "concentration"}, fmt=ST.fmt(c),
            h2h_beats_pool_member={k: float(np.mean([1 - r["pct"][k] for r in rows])) for k in names},
            pct_mean={k: float(np.mean([r["pct"][k] for r in rows])) for k in names},
            pool_source=rows[0].get("pool_source"))
        if "RAND_SIGNED[0]" in names:
            ind_rs = np.array([r["pref"]["RAND_SIGNED[0]"] for r in rows])
            c2 = ST.compare(ind_cb, ind_rs, folds, names=pdbs, label=f"{cost.name}: pref(ORACLE circ_best) - pref(RAND_SIGNED[0])",
                            seed_parts=("s29D",))
            out["pool_member_control"]["contrast_vs_rand_signed"] = {k: v for k, v in c2.items() if k != "concentration"}
            out["pool_member_control"]["fmt_vs_rand_signed"] = ST.fmt(c2)
    else:
        out["pool_member_control"] = None
    out["nan_rows"] = [r["pdb"] for r in rows if r.get("nan_on_rungs")]
    out["text"] = render(out)
    return out


def render(o):
    L = []
    L.append(f"COST-RMSD METER (ORACLE DIAGNOSTIC; every number reads the native): cost {o['cost']} [{o['kind']}], basis {o['basis']}, n={o['n']}")
    L.append("  rung mean ORACLE RMSD / mean f: " + "  ".join(f"{k} {o['mean_rmsd'][k]:.3f}/{o['mean_f'][k]:.4g}" for k in o["rungs"]))
    L.append("  (a) LADDER Spearman rho(f, RMSD) per target, mean [fold CI] (+1 = f rises with RMSD, the good sign):")
    for nm, s in o["ladder_rho"].items():
        if s.get("ci95_fold"):
            L.append(f"      {nm:8s} {s['mean']:+.4f} [{s['ci95_fold'][0]:+.3f}, {s['ci95_fold'][1]:+.3f}]  SE {s['se']:.4f}  median {s['median']:+.3f}  "
                     f"pos {s['n_pos']}/{s['n']}  FAIL18 {s['fail18_mean']:+.3f} / 108 {s['other108_mean']:+.3f}  "
                     f"(undefined->0 as C2: {s['n_undefined_set_to_0_as_C2']}, mean {s['mean_C2_convention']:+.4f})  rungs {o['ladders'][nm]}")
    c = o.get("cosine")
    if c and c["ci95_fold"] is None:
        L.append(f"  (b) COSINE: UNDEFINED on {c['n_nan']} of {c['n'] + c['n_nan']} targets [{c['grad_method']}] -- "
                 f"{100*c['grad_zero_frac_mean']:.1f}% of the gradient components are exactly zero. This cost is a step "
                 f"function at this step size: re-run with --fd-h 0.5 for A2's SMOOTHED finite difference (S28-L23b), "
                 f"and never quote a cosine computed on the surviving targets -- they are the ones nearest a discontinuity.")
    elif c:
        L.append(f"  (b) COSINE cos(-grad f at production, native direction), rigid body removed [{c['grad_method']}]"
                 + (f" -- UNDEFINED on {c['n_nan']} of {c['n'] + c['n_nan']} targets (zero gradient); the mean below is on the survivors and is a SELECTED SUBSAMPLE, not a measurement" if c["n_nan"] else "") + ":")
        L.append(f"      mean {c['mean']:+.4f}  SE {c['se']:.4f}  median {c['median']:+.4f}  fold CI [{c['ci95_fold'][0]:+.3f}, {c['ci95_fold'][1]:+.3f}]  "
                 f"pos {c['n_pos']}/{c['n']}  FAIL18 {c['fail18_mean']:+.3f} / 108 {c['other108_mean']:+.3f}  "
                 f"| random-direction |cos| mean {c['random_ref_mean_abs']:.3f} p95 {c['random_ref_p95_abs']:.3f}  "
                 f"| grad RMS {c['grad_rms_mean']:.3g}, zero components {100*c['grad_zero_frac_mean']:.0f}%  | Spearman(cos, prod RMSD) {c['spearman_cos_vs_prod_rmsd']:+.3f}")
        sh = c.get("shrink_signature")
        if sh:
            L.append(f"      SHRINK SIGNATURE (contract addendum 20, S29-L7): descending {sh['e']} A along -grad f multiplies the mean virtual bond by "
                     f"{sh['bond_ratio_mean']:.4f} and Rg by {sh['rg_ratio_mean']:.4f} ({sh['n_contracting']}/{sh['n']} targets contract; production's own bond "
                     f"{sh['bond_prod_mean']:.2f} A, Rg {sh['rg_prod_mean']:.2f} A, already 22% contracted). A cosine gain with a ratio below 1 is a shrink, "
                     f"not information: read it beside (c), which the shrink moves the wrong way.")
    else:
        L.append("  (b) COSINE: not defined on this basis (the projection is not differentiable)")
    p = o.get("native_pctile")
    if p and p["ci95_fold"] is None:
        L.append(f"  (c) PCTILE: UNDEFINED on {p['n_nan']} of {p['n'] + p['n_nan']} targets (NaN cost on the native or the pool)")
        p = None
    if p:
        L.append(f"  (c) PCTILE native's percentile in its own pool under f (share of members scoring better, half ties): "
                 f"mean {p['mean']:.4f} [{p['ci95_fold'][0]:.3f}, {p['ci95_fold'][1]:.3f}]  (strict, objdiag: {p['mean_strict_objdiag']:.4f})  "
                 f"median {p['median']:.3f}  FAIL18 {p['fail18_mean']:.3f} / 108 {p['other108_mean']:.3f}  | PROD below native on {p['n_prod_below_native']}/{o['n']}, "
                 f"below the pool's best on {p['n_prod_below_pool_min']}/{o['n']}  | mean f: native {p['mean_f_native_posed']:.4g} pool-min {p['mean_pool_min']:.4g} pool-mean {p['mean_pool_mean']:.4g}")
    else:
        L.append("  (c) PCTILE: not available (no CA pool channel for this cost)")
    L.append("  (d) PREF share of targets on which f scores X below PROD (ties 0.5), mean [fold CI]:")
    for k, s in o["pref"].items():
        ci = ("[%.3f, %.3f]" % tuple(s["ci95_fold"])) if s["ci95_fold"] else "[CI undefined: %d finite]" % s["n"]
        L.append(f"      {k:16s} {s['mean']:.3f} {ci}  ties {s['ties']}  FAIL18 {s['fail18_mean']:.2f} / 108 {s['other108_mean']:.2f}")
    m = o.get("pool_member_control")
    if m:
        cc = m["contrast"]
        L.append(f"      POOL-MEMBER CONTROL (S28-L36): pref(circ_best vs PROD) {m['pref_circ_best_vs_prod']:.3f}  vs  pref(random pool member vs PROD) = pct(PROD) {m['pref_pool_member_vs_prod']:.3f} (median {m['pct_prod_median']:.2f})")
        L.append(f"        contrast {cc['effect']:+.4f}  SE {cc['se']:.4f}  MDE {cc['mde']:.4f}  {cc['effect_over_mde']:+.2f}x  fold CI [{cc['ci95_fold'][0]:+.3f}, {cc['ci95_fold'][1]:+.3f}]  folds same sign {cc['folds_same_sign']}/{cc['n_folds']}   "
                 f"[pool: {m['pool_source']}]")
        L.append("        head-to-head, X scored below a random pool member: " + "  ".join(f"{k} {v:.3f}" for k, v in m["h2h_beats_pool_member"].items()))
        if "contrast_vs_rand_signed" in m:
            c2 = m["contrast_vs_rand_signed"]
            L.append(f"        vs RAND_SIGNED[0] (C2 clause 2): pref(circ_best) - pref(RAND_SIGNED) {c2['effect']:+.4f}  {c2['effect_over_mde']:+.2f}x  fold CI [{c2['ci95_fold'][0]:+.3f}, {c2['ci95_fold'][1]:+.3f}]")
    if o.get("nan_rows"):
        L.append(f"  NaN on some rung for {len(o['nan_rows'])} targets: {o['nan_rows'][:10]}")
    L.append("  The meter tunes nothing; a cost passes to an endpoint run only on its OWN pre-registered falsifier, and every number above is ORACLE.")
    return "\n".join(L)


# ============================================================================ driver
def run_meter(spec, basis="ca", limit=0, fd_h=FD_H_DEFAULT, tag=None, out=None, pdbs=None, quiet=False, save=True):
    from s25 import phys_lib as P
    cost = Cost(spec, fd_h=fd_h)
    if pdbs is None:
        pdbs = P.targets()[:limit] if limit else P.targets()
    rows = []
    t0 = time.time()
    for q, pdb in enumerate(pdbs):
        if basis == "chain-s28rows":
            row = meter_target_from_s28_rows(pdb, cost)
        else:
            row = meter_target(pdb, cost, basis=basis)
        if row.get("nan_on_rungs") and cost.kind == "callable":
            raise ValueError(f"{pdb}: the cost gave NaN on {row['nan_on_rungs']} under the POISONED context "
                             "(ctx.nat_ca / ctx.oracle_rr are NaN): it reads a native quantity or is undefined; "
                             "the meter refuses it. Errors: " + "; ".join(cost.errors[:4]))
        rows.append(row)
        if not quiet and ((q + 1) % 10 == 0 or q + 1 == len(pdbs)):
            print(f"  [{q+1}/{len(pdbs)}] {pdb} rho S28 {row['rho']['S28']:+.3f} CHARTER {row['rho']['CHARTER']:+.3f}"
                  + (f" cos {row['cos']:+.3f}" if "cos" in row else "")
                  + (f" pct {row['native_pctile']:.3f}" if "native_pctile" in row else "")
                  + f"  ({(time.time()-t0)/60:.1f} min)", flush=True)
    summ = summarise(rows, cost, basis)
    summ["rows"] = rows
    summ["fd_h"] = fd_h
    if save:
        os.makedirs(RESULTS, exist_ok=True)
        tag = tag or cost.name.replace("/", "_")
        path = out or os.path.join(RESULTS, f"s29_D_cost_audit_{tag}_{basis}.json")
        ST.save_atomic(path, summ, rows=rows, n_expected=len(pdbs), complete_keys=("pdb", "f", "rho", "pref"), module_file=__file__)
        summ["path"] = path
    if not quiet:
        print(summ["text"])
        if save:
            print("wrote", summ["path"])
    return summ


# ============================================================================ selftest
def selftest():
    """Synthetic: a cost equal to the ORACLE RMSD must meter at rho +1, cos +1, pctile 0, pref 1;
    a cost that ignores structure meters at rho NaN (undefined), pref 0.5 everywhere."""
    rng = np.random.default_rng(0)
    n = 8
    nat = np.cumsum(rng.normal(size=(n, 3)) * 3.0, 0)
    C0 = nat + rng.normal(scale=1.5, size=(n, 3))
    u = A2.remove_rigid(A2.oracle_direction(C0, nat), C0)
    r, g = A.oracle_rmsd_grad(C0, nat)
    cos = A2.cosine(-A2.remove_rigid(g, C0), u)
    assert cos > 0.95, cos
    assert abs(spearman([1, 2, 3, 4], [1, 2, 3, 4]) - 1.0) < 1e-12
    assert math.isnan(spearman([1, 1, 1], [1, 2, 3]))
    assert pct_in_pool(np.array([1.0, 2.0, 3.0]), 2.0) == 0.5
    print("  s29_D_cost_audit selftest OK (the RMSD-gradient cosine at a random point is %.3f)" % cos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["build-cache", "meter", "selftest"])
    ap.add_argument("--f", default="DIS", help="cost name (DIS, DIS_SURR, an S27 channel) or module:function")
    ap.add_argument("--basis", default="ca", choices=["ca", "chain", "chain-s28rows"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--fd-h", type=float, default=FD_H_DEFAULT)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--chain", action="store_true", help="build-cache: also store every rung's production projection")
    ap.add_argument("--rebuild", action="store_true")
    a = ap.parse_args()
    if a.mode == "selftest":
        selftest()
        return
    from s25 import phys_lib as P
    pdbs = P.targets()[:a.limit] if a.limit else P.targets()
    if a.mode == "build-cache":
        build_cache(pdbs, with_chain=a.chain, rebuild=a.rebuild)
    else:
        run_meter(a.f, basis=a.basis, limit=a.limit, fd_h=a.fd_h, tag=a.tag, out=a.out, pdbs=pdbs)


if __name__ == "__main__":
    main()
