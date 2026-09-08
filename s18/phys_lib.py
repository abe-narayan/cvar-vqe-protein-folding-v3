"""s18/phys_lib.py -- the PHYSICS workstream's shared machinery for Sprint 18.

Pre-registration: `s18/PREREG_phys.md`, written before any number here existed.

WHAT IS IN HERE

  * `ContactTerm`   -- `core.energy.contact_term` (genuine Legacy, MJ-corrected, weight 1.0,
                       NEVER fitted) lifted into continuous torsion space, with a batched
                       finite-difference torsion gradient.  Gate G0 checks it bit-for-bit
                       against `s16.energy_lib.legacy_components_of_windows`; gate G1 checks
                       the gradient machinery against `s15.align_lib`'s exact analytic
                       gradient on the pure-distance objective.
  * `pool_scale`    -- the pre-registered NATIVE-FREE normalisation:
                       s = sd_pool(E_le1_dist) / sd_pool(E_contact) over the ranker-neutral
                       K = 500 BLOSUM62 retrieval pool.  `lam_c = 1` then means one pool
                       standard deviation against one pool standard deviation.
  * `shuffled_mj`   -- the zero-information control for a CONTACT term.  Shuffling the
                       distogram controls the DISTANCE term and says nothing about this one;
                       the right null here is a residue-label permutation of the MJ matrix,
                       which keeps the functional form and the magnitude and destroys the
                       sequence information.
  * statistics      -- paired fold-clustered bootstrap with median and W/L (Sprint 17's
                       `phys_lib.paired`, carried forward verbatim in behaviour), the gated
                       comparison that refuses to be separated from its exclusion count, and
                       `write` with an explicit `complete` flag.

WHY FINITE DIFFERENCES AND NOT AN ANALYTIC GRADIENT.  The contact energy is a function of CB,
and CB is a function of (N, CA, C) of the same residue; `core.project._torsion_grad` chains a
gradient defined on CA only.  Writing the N/C chain rule by hand is the kind of thing that
produces a wrong number that looks right.  With n <= 16 the whole central-difference stencil is
4n+1 <= 65 structures and `core.geometry.build_backbone_batch` builds all of them in one
vectorised pass, so the exact-looking route buys nothing and risks a silent defect.  The cosine
switch is C1, so the differences are accurate; G1 measures that rather than assuming it.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from typing import Dict, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                      # noqa: E402
from s15 import seed as SD                           # noqa: E402

ATOMS = ("N", "CA", "C", "O", "CB")

#: The pre-registered ladder.  SIX arms.  Not extended, whatever the result.
LAM_LADDER = (-1.0, -0.5, 0.0, 0.5, 1.0, 2.0)
LAM_PRIMARY = 1.0

#: `core.energy.contact_term`'s own constants, restated here only so a reader can see them.
CONTACT_D0, CONTACT_DC = 4.5, 8.5
CONTACT_MIN_SEP = 3

FD_H = 1e-5          # radians, central differences


# ==========================================================================
# 1.  THE CONTACT TERM IN CONTINUOUS TORSION SPACE
# ==========================================================================
class ContactTerm:
    """`leg_contact` as a differentiable function of (phi, psi).

    Genuine Legacy: the MJ-corrected table from `core.energy`, the cosine switch at
    (4.5, 8.5) A over CB pairs at |i-j| >= 3, at `DEFAULT_WEIGHTS["contact"] = 1.0`.
    Nothing here is fitted and nothing is a surrogate.
    """

    def __init__(self, seq: str, mj: Optional[np.ndarray] = None):
        from core import energy as et
        self.seq = seq
        self.n = len(seq)
        di, dj, sep = et.pair_index(self.n)
        m3 = sep >= CONTACT_MIN_SEP
        self.di = np.asarray(di[m3], int)
        self.dj = np.asarray(dj[m3], int)
        if mj is None:
            _b, _q, MJ = et.sequence_arrays(seq, True)
        else:
            MJ = np.asarray(mj, float)
        self.mj_pair = np.asarray(MJ[self.di, self.dj], float)
        self.weight = float(et.DEFAULT_WEIGHTS["contact"])

    # ---------------------------------------------------------------- energy
    def of_cb(self, CB: np.ndarray) -> np.ndarray:
        """Energy of one (n,3) or a batch (B,n,3) of CB sets.  Unweighted, as Legacy reports it."""
        from core import energy as et
        CB = np.asarray(CB, float)
        single = CB.ndim == 2
        if single:
            CB = CB[None]
        d = np.linalg.norm(CB[:, self.di] - CB[:, self.dj], axis=2)
        v = et.switch(d, CONTACT_D0, CONTACT_DC) @ self.mj_pair
        return float(v[0]) if single else v

    def of_torsions(self, phi, psi):
        """Energy of one (n,) or a batch (B,n) torsion set."""
        from core import geometry as geo
        phi = np.asarray(phi, float)
        psi = np.asarray(psi, float)
        single = phi.ndim == 1
        PH = phi[None] if single else phi
        PS = psi[None] if single else psi
        CB = geo.build_backbone_batch(PH, PS)["CB"]
        out = self.of_cb(CB)
        return float(np.atleast_1d(out)[0]) if single else np.asarray(out, float)

    # -------------------------------------------------------------- gradient
    def fg(self, phi, psi, h: float = FD_H):
        """(value, gradient over the 2n torsions) by BATCHED central differences.

        One `build_backbone_batch` call of 4n+1 structures per evaluation.
        """
        phi = np.asarray(phi, float)
        psi = np.asarray(psi, float)
        n = self.n
        x = np.concatenate([phi, psi])
        X = np.repeat(x[None], 4 * n + 1, axis=0)
        for k in range(2 * n):
            X[1 + 2 * k, k] += h
            X[2 + 2 * k, k] -= h
        E = self.of_torsions(X[:, :n], X[:, n:])
        g = (E[1::2] - E[2::2]) / (2.0 * h)
        return float(E[0]), g


def shuffled_mj(seq: str, rng) -> np.ndarray:
    """ZERO-INFORMATION control for the contact term.

    A permutation of the residue labels feeding the MJ table.  The functional form, the pair
    set, the switch and the magnitude distribution are all preserved; the map from *which*
    residue sits at position i to its contact preference is destroyed.  This is the control a
    contact term needs; a shuffled distogram is the control the DISTANCE term needs and is
    silent about this one.
    """
    from core import energy as et
    perm = rng.permutation(len(seq))
    shuffled = "".join(seq[k] for k in perm)
    _b, _q, MJ = et.sequence_arrays(shuffled, True)
    return np.asarray(MJ, float)


# ==========================================================================
# 2.  THE PRE-REGISTERED NATIVE-FREE NORMALISATION
# ==========================================================================
def pool_scale(e_dist_pool: np.ndarray, e_contact_pool: np.ndarray) -> Dict[str, float]:
    """s = sd_pool(E_dist) / sd_pool(E_contact), both over the SAME candidate pool.

    NATIVE-FREE by construction: only candidate energies enter.  Returns the two standard
    deviations beside the ratio so a reader can see what set the scale, and a `degenerate`
    flag rather than a silent divide-by-zero.
    """
    a = np.asarray(e_dist_pool, float)
    b = np.asarray(e_contact_pool, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    sa = float(a.std()) if len(a) else float("nan")
    sb = float(b.std()) if len(b) else float("nan")
    deg = (not np.isfinite(sa)) or (not np.isfinite(sb)) or sb < 1e-9
    return {"sd_dist": sa, "sd_contact": sb,
            "s": float("nan") if deg else sa / sb, "degenerate": bool(deg),
            "n_pool": int(min(len(a), len(b)))}


# ==========================================================================
# 3.  GATES -- these run before any scientific number is quoted
# ==========================================================================
def gate_G0(n_targets: int = 8, k: int = 64, verbose: bool = True) -> Dict[str, object]:
    """G0 -- the torsion-space contact term IS genuine Legacy's `contact` component.

    Compared against `s16.energy_lib.legacy_components_of_windows`, which is the call Sprint 17
    used for every Legacy number on the record.
    """
    from s16 import energy_lib as EL
    tg = I.targets()[:n_targets]
    worst = 0.0
    rows = []
    for t in tg:
        pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
        u = I.load_univ(pdb)
        order = np.asarray(u["order"], int)[:k]
        PHI = np.asarray(u["PHI"], float)[order]
        PSI = np.asarray(u["PSI"], float)[order]
        ref = np.asarray(EL.legacy_components_of_windows(seq, PHI, PSI)["contact"], float)
        mine = ContactTerm(seq).of_torsions(PHI, PSI)
        e = float(np.abs(ref - mine).max())
        worst = max(worst, e)
        rows.append({"pdb": pdb, "n": n, "k": int(len(order)), "max_abs_err": e})
        if verbose:
            print(f"  G0 {pdb:>6}  max|dE| = {e:.3e}", flush=True)
    out = {"gate": "G0", "max_abs_err": worst, "tol": 1e-9,
           "passed": bool(worst < 1e-9), "rows": rows}
    if verbose:
        print(f"  G0 {'PASS' if out['passed'] else 'FAIL'}: max|dE| = {worst:.3e} over "
              f"{len(rows)} targets x {k} candidates")
    return out


def gate_G1(n_targets: int = 6, verbose: bool = True) -> Dict[str, object]:
    """G1 -- the finite-difference torsion gradient machinery is correct.

    Measured on the PURE DISTANCE objective, where `s15.align_lib` supplies an exact analytic
    gradient, so the check is of the differencing and not of a second hand-written chain rule.
    The same stencil is then reused for the contact term.
    """
    from s15 import distcal as C
    from core import project as pj
    from s15.robust import rho

    tg = I.targets()[:n_targets]
    data = C.gather(tg)
    worst = 0.0
    rows = []
    for t in tg:
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        d = data[pdb]
        i, j, sd = d["i"], d["j"], d["sd"]
        dhat = np.maximum(np.asarray(d["dhat"], float), 2.0)
        rng = SD.stable_rng(pdb, "s18G1")
        phi = rng.uniform(-np.pi, np.pi, n)
        psi = rng.uniform(-np.pi, np.pi, n)
        inv = 1.0 / np.asarray(sd, float)

        def f_only(PH, PS):
            CA = np.asarray(pj.build_ca_exact(np.atleast_2d(PH), np.atleast_2d(PS)), float)
            rv = CA[:, i] - CA[:, j]
            dd = np.maximum(np.sqrt((rv * rv).sum(2)), 1e-9)
            v, _g = rho((dd - dhat) * inv, "squared", 1.0)
            return v.sum(1)

        # exact analytic gradient, the one `align_lib.fit` optimises with
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        rv = CA[i] - CA[j]
        dd = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
        r = (dd - dhat) * inv
        _v, gr = rho(r, "squared", 1.0)
        coef = ((gr * inv) / dd)[:, None] * rv
        gCA = np.zeros_like(CA)
        np.add.at(gCA, i, coef)
        np.add.at(gCA, j, -coef)
        g_exact = np.asarray(pj._torsion_grad(G, CA, gCA), float)

        # the SAME stencil `ContactTerm.fg` uses
        x = np.concatenate([phi, psi])
        X = np.repeat(x[None], 4 * n + 1, axis=0)
        for kk in range(2 * n):
            X[1 + 2 * kk, kk] += FD_H
            X[2 + 2 * kk, kk] -= FD_H
        E = f_only(X[:, :n], X[:, n:])
        g_fd = (E[1::2] - E[2::2]) / (2.0 * FD_H)

        rel = float(np.abs(g_fd - g_exact).max() / max(np.abs(g_exact).max(), 1e-12))
        worst = max(worst, rel)
        rows.append({"pdb": pdb, "n": n, "rel_err": rel})
        if verbose:
            print(f"  G1 {pdb:>6}  rel|dg| = {rel:.3e}", flush=True)
    out = {"gate": "G1", "max_rel_err": worst, "tol": 1e-5,
           "passed": bool(worst < 1e-5), "rows": rows}
    if verbose:
        print(f"  G1 {'PASS' if out['passed'] else 'FAIL'}: max rel err = {worst:.3e}")
    return out


# ==========================================================================
# 4.  STATISTICS -- target as the unit, always with median and W/L
# ==========================================================================
def paired(a, b, folds=None, names=None, seed=0, n_boot=4000):
    """Paired a - b (negative = `a` better).  i.i.d. bootstrap CI AND a fold-clustered twin.

    Carried forward from `s17/phys_lib.paired` unchanged in behaviour so Sprint 17 and Sprint 18
    numbers are directly comparable.  Non-finite pairs are dropped and the surviving `n` is
    reported, so a shrinking denominator is visible rather than silent.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if not m.all():
        folds = None if folds is None else np.asarray(folds)[m]
        names = None if names is None else [names[i] for i in np.flatnonzero(m)]
        a, b = a[m], b[m]
    d = a - b
    n = len(d)
    if n < 3:
        return {"n": int(n), "mean": float("nan"), "median": float("nan"),
                "ci": [float("nan")] * 2, "W": 0, "L": 0,
                "mean_a": float("nan"), "mean_b": float("nan")}
    rng = np.random.default_rng(seed)
    bs = d[rng.integers(0, n, size=(n_boot, n))].mean(1)
    out = {"n": int(n), "mean_a": float(a.mean()), "mean_b": float(b.mean()),
           "median_a": float(np.median(a)), "median_b": float(np.median(b)),
           "mean": float(d.mean()), "median": float(np.median(d)),
           "ci": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
           "W": int((d < 0).sum()), "L": int((d > 0).sum())}
    if folds is not None:
        f = np.asarray(folds)
        uf = np.unique(f)
        idx = {u: np.flatnonzero(f == u) for u in uf}
        cb = np.empty(n_boot)
        for b_ in range(n_boot):
            pick = rng.integers(0, len(uf), len(uf))
            cb[b_] = d[np.concatenate([idx[uf[p]] for p in pick])].mean()
        out["ci_fold"] = [float(np.percentile(cb, 2.5)), float(np.percentile(cb, 97.5))]
        out["per_fold"] = {int(u): float(d[idx[u]].mean()) for u in uf}
    if names is not None:
        o = np.argsort(d)
        out["top10"] = [(str(names[k]), float(d[k])) for k in o[:10]]
    return out


def gated(a, b, energies, converged=None, **kw):
    """A gated comparison AND its ungated twin AND the exclusion count -- inseparable.

    THE CONVERGENCE GATE, declared in `PREREG_phys.md` before use: final energy finite and
    <= `core.amber.CONVERGE_MAX_KCAL` (1000 kcal/mol).  `a` and `b` must ALREADY be the arm and
    **its own gated input** -- this function does not know which is which, and comparing a gated
    arm to an ungated instrument-wide baseline is a caller error the brief names explicitly.
    """
    e = np.asarray(energies, float)
    ok = np.isfinite(e) & (e <= 1000.0)
    if converged is not None:
        ok &= np.asarray(converged, bool)
    fl = kw.get("folds")
    nm = kw.get("names")
    out = {"gated": paired(np.asarray(a, float)[ok], np.asarray(b, float)[ok],
                           folds=(None if fl is None else np.asarray(fl)[ok]),
                           names=(None if nm is None else [nm[i] for i in np.flatnonzero(ok)]),
                           seed=kw.get("seed", 0)),
           "ungated": paired(a, b, folds=fl, names=nm, seed=kw.get("seed", 0)),
           "n_excluded": int((~ok).sum())}
    if nm is not None:
        out["excluded"] = [str(nm[i]) for i in np.flatnonzero(~ok)]
    return out


def spearman(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return float("nan")
    from scipy.stats import rankdata
    a = rankdata(x[m])
    b = rankdata(y[m])
    if a.std() < 1e-12 or b.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def argmin_tied(score, value, rng=None):
    """Mean `value` over the FULL tied argmin set.

    The recorded trap: `np.argmin` on a tied score reads the cache's sort order, which on this
    project is sometimes the ORACLE order, and once invented a 1.386 A winner.  Averaging over
    the tie set is the fix the ledger prescribes.
    """
    s = np.asarray(score, float)
    v = np.asarray(value, float)
    m = np.isfinite(s)
    if not m.any():
        return float("nan")
    lo = np.nanmin(s[m])
    tie = np.flatnonzero(m & (s <= lo + 1e-12))
    return float(v[tie].mean())


def diversity(W: np.ndarray) -> float:
    """Mean pairwise Ca-RMSD inside an ensemble of (B, n, 3) structures."""
    W = np.asarray(W, float)
    if len(W) < 2:
        return float("nan")
    P = I.pairwise_rmsd(W)
    iu = np.triu_indices(len(W), 1)
    return float(np.asarray(P)[iu].mean())


# ==========================================================================
# 5.  ARTEFACTS
# ==========================================================================
def cfg_hash(obj) -> str:
    return hashlib.blake2b(json.dumps(obj, sort_keys=True, default=str).encode(),
                           digest_size=8).hexdigest()


def write(name: str, obj: dict, n_expected: Optional[int] = None, subdir: str = ""):
    """Write with an EXPLICIT `complete` flag and row count.

    Sprint 17 lost a finished 126-row artefact to a 9-row partial written by a racing process
    (ledger L29).  Nothing in this workstream writes a result file without the flag.
    """
    d = os.path.join(RESULTS, subdir) if subdir else RESULTS
    os.makedirs(d, exist_ok=True)
    rows = obj.get("rows") or obj.get("per_target") or []
    obj = dict(obj, n_rows=len(rows))
    if n_expected is not None:
        obj["n_expected"] = int(n_expected)
        obj["complete"] = bool(len(rows) >= int(n_expected))
    else:
        obj.setdefault("complete", True)
    p = os.path.join(d, name if name.endswith(".json") else name + ".json")
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, p)
    return p


def read_complete(path: str, need: Optional[int] = None):
    """Refuse to read a partial artefact as a complete one."""
    with open(path) as fh:
        o = json.load(fh)
    rows = o.get("rows") or o.get("per_target") or []
    if not o.get("complete", False):
        raise RuntimeError(f"{path}: artefact is NOT complete ({len(rows)} rows)")
    if need is not None and len(rows) < need:
        raise RuntimeError(f"{path}: {len(rows)} rows < {need} expected")
    return o


def mem_ok(need_gb: float = 1.6) -> float:
    """Yield until `need_gb` of physical memory is free.  Three agents share ~11 GB."""
    import time
    for _ in range(600):
        f = I.free_gb()
        if f >= need_gb:
            return f
        time.sleep(2.0)
    return I.free_gb()


if __name__ == "__main__":
    print("PHYSICS gates -- these run before any scientific number is quoted.\n")
    g0 = gate_G0()
    g1 = gate_G1()
    write("gates_phys", {"G0": g0, "G1": g1,
                         "passed": bool(g0["passed"] and g1["passed"])})
    print(f"\nBOTH GATES {'PASS' if g0['passed'] and g1['passed'] else 'FAIL'}")
