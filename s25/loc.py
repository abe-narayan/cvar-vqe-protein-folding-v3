"""s25/loc.py -- THE LOCATION DIAGNOSTIC.  WHERE DOES THE SCORE'S TARGET SIT, AND CAN IT BE MOVED?

NO RMSD IS COMPUTED IN THIS FILE.  Natives are read for DIAGNOSIS ONLY -- nothing is fitted, nothing
is selected, no production object is modified.  Same standing as `s25/calib.py`.  Pre-registered in
`s25/PREREG_RMSD.md`; the fork list went to the coordinator before this ran.

WHY.  The shipped score is `risk(t) = sum_c p_c |t - C_c|`, an L1 functional whose minimiser is the
posterior MEDIAN, so it has a per-pair TARGET and that target is an object one can name and measure.
This file measures where that target sits and how far any alternative reading of the same posterior
moves it.

**CORRECTION, ISSUED AFTER THIS FILE RAN AND KEPT HERE RATHER THAN QUIETLY EDITED.**  This module was
written under s25 L2's framing -- "width is inert, LOCATION is not".  **s25 L7 RETRACTS that**: the
unconfounded experiment (mean-preserving widening versus pure translation, n=126) finds location and
width EQUALLY FLAT, all eight cells at 0.38-0.58x their own MDE, and the response to location
SYMMETRIC about a posterior that is biased -0.41 A against the natives.  See `s25/LEDGER.md` L7.
**Nothing measured in this file depends on that framing** -- it computes where the target is and how
far each functional moves it, which is true either way -- but the motivating sentence was wrong and
the reader should know the prior on every family here is lower than it was when they were chosen.

THE INCUMBENT'S TARGET IS AN OBJECT WE CAN NAME AND MEASURE.  It is `grid[argmin_t risk_p[t]]` -- the
posterior median, on the deployed 0.05 A grid.  Five one-parameter families move that statistic, each
containing the incumbent as a BIT-EXACT interior point:

    A1 POWER  rho(u) = |u|^q               q->0 mode,  q=1 MEDIAN (identity),  q=2 mean
    A2 TRUNC  rho(u) = min(|u|, delta)     delta->0 mode-seeking,  delta=40 >= max|grid-C| = identity
    B  QUANT  rho(u) = 2[tau u+ + (1-tau) u-]   tau = 0.5 is identity; tau<0.5 targets SHORTER
    C  SHELL  location shifted by alpha * off(sep)     alpha = 0 is identity
    D  METRIC location moved toward the rank-3 EDM projection of the whole field;  eta = 0 identity

A1/A2 are the SHAPE axis (multimodality, and "a different risk functional").  B/C are the SHIFT axis
in two deliberately different normalisations -- B self-scaled by the posterior's own spread, C in
absolute Angstroms off a per-separation profile.  D is the CONSTRAINT axis and is the only arm in
this lane that uses information across pairs rather than reshaping each pair independently.

WHAT IS MEASURED.  For every candidate location `L`, against the move to truth, per pair:

    dL = L - L_med        dT = D_native - L_med
    gam_eff = <dL,dT>/|dT|^2      the fraction of the way to truth actually travelled
    cos     = <dL,dT>/(|dL||dT|)  how much of the move is aligned at all
    |dL|/|dT|                     how far it moved

`gam_eff` is the PRIOR LADDER'S OWN CURRENCY: the ladder's TILT arm is location-only and tracks MASS
to 0.024 A at gamma=0.1, so gamma is a legitimate unit for a pure location move.  **The ladder's
gamma has cos = 1 by construction.  Any real functional has cos < 1 and its orthogonal component adds
error the ladder never carried.  Therefore -2.1496 * gam_eff is an UPPER BOUND on the endpoint gain,
never a prediction, and it is labelled as a bound at every appearance.**

THREE FURTHER MEASUREMENTS, each able to kill a family before it costs an end-to-end run:

  D1 THE CANCELLATION TEST -- the coordinator's warning made falsifiable.  s24 Workstream C measured
     the distogram's signed offset at +0.4062 and the whole legal universe of REAL protein windows at
     +0.4419 against these same natives, so the offset is a corpus/native scale mismatch the
     CANDIDATES SHARE.  PREDICTION STATED BEFORE THE RUN: pool_signed[shell] matches prior_signed
     [shell] to within ~0.1 A at every shell, and families B and C therefore FAIL end-to-end despite
     a positive gam_eff, because shifting only the prior moves the target away from where the
     candidates live and breaks a cancellation that currently helps.
     ADDITION 4 (coordinator, before run 1): the offset is reported for the full K=500 pool AND for
     the shipped-score TOP-75, because selection could move it and the candidates that actually
     reach the readout are the ones whose cancellation matters.
  D2 THE MULTIMODAL SUB-POPULATION.  On pairs with >= 2 peaks above 0.02 (24.1% of all pairs, s25
     L1): is the native closer to the NEAREST MODE than to the median?  If not, the mode-seeking
     mechanism has no headroom and A2 is dead before it runs.
     ADDITION 2 (coordinator, PRE-REGISTERED BEFORE RUN 1, not requested after the numbers existed):
     EVERY family's gam_eff and cosine are reported STRATIFIED by that criterion, and the end-to-end
     module will report its delta the same way.  Prediction 3 says family A is open precisely because
     the multimodal pairs live there, so the stratification IS the mechanism test -- but only if it
     is declared in advance.  s24 Lane C correctly declined to promote a stratification requested
     after its numbers existed; declaring this one now costs nothing and makes it usable.
  D4 GRID QUANTISATION (ADDITION 3, coordinator).  The risk table is evaluated on 760 grid points at
     0.05 A but the posterior lives on 17 IRREGULAR centres, 0.5 A apart at the short end and 3.5 A
     at the long end.  As q -> 0 the minimiser approaches the mode of a BINNED posterior, so small q
     may be measuring bin quantisation rather than mode-seeking, and the artefact would be strongly
     separation-dependent because the bins are.  Every arm reports the realised minimiser's distance
     to the nearest bin centre and the fraction sitting exactly on one; if it collapses onto centres
     the small-q end is bounded and reported as quantisation, not read as a result.
  D3 ZERO-INFORMATION CONTROLS MATCHED IN THE OPERATOR'S SPACE (`control-must-match-the-operators-
     space`, the project's most repeated error).  (i) PERM: the arm's own per-pair move vector
     permuted within the target -- identical magnitude distribution, no pair correspondence.
     (ii) XSHELL: the per-shell offset profile of a different, randomly assigned target.

RULE 0 -- SIX OPERATOR FORKS, each naming the alternative NOT taken.  Enumerated by me, who has a
stake in the answer; the audit lane should re-enumerate.

    functional     DECLARED every location computed as the argmin over the DEPLOYED 0.05 A grid of a
                   risk row built from the shipped `prob` with only the loss KERNEL replaced, so the
                   incumbent's own location is produced by the same code path as every alternative.
                   Ties are averaged over the argmin set (`consensus-is-the-only-in-band-
                   discriminator`: np.argmin on a tied signal reads the array's order).
                   NOT TAKEN reading `dg["expected"]` as "the location", which is the posterior MEAN
                   and is NOT what the shipped L1 risk minimises -- that substitution is exactly the
                   estimator mismatch s24 L9-A had to reconcile between two lanes.
    basis          DECLARED distance space over the score's own pair set, min_sep=2 -- the space the
                   functional consumes.  NOT TAKEN Cartesian space, which the prior does not have,
                   and NOT TAKEN RMSD, which this file does not compute at all.
    readout        DECLARED per-pair statistics aggregated per TARGET and then averaged over targets,
                   so long targets do not dominate.  NOT TAKEN pooling all pairs of all targets,
                   which weights by n^2 (reported beside it as `pooled` so the difference is visible).
    normalisation  DECLARED gam_eff UNWEIGHTED over pairs as primary, with the w-weighted version
                   beside it.  NOT TAKEN weighting by w silently, which would let the shipped
                   confidence weighting -- a channel s24 L6 and s25 L2 both found inert -- decide
                   which family looks best.  NOTE THE STANDING DEPENDENCY: the "width is inert"
                   framing rests on s25 L2, and the audit lane is checking whether w = shell/(sd+0.5)^g
                   makes width non-inert through the WEIGHT even when it is inert through the risk
                   curve's SHAPE.  Reporting both normalisations is what makes this lane survive
                   either answer.
    null           DECLARED the identity parameter of each family, which must reproduce the shipped
                   risk table BIT-EXACTLY through this module's own kernel path (asserted per target,
                   np.array_equal), plus the two matched zero-information moves of D3.  NOT TAKEN
                   assuming the identity holds because the algebra says so.
    THE LABEL      DECLARED gam_eff AND cos AND |dL| AND signed error AND MAE, reported together,
                   because a move can travel toward truth on average while being mostly orthogonal
                   noise, and the two have opposite consequences.  NOT TAKEN MAE alone -- `prior-mae-
                   prices-selected-rmsd` prices MAE's correlation with selected RMSD at r = 0.19, and
                   `error-shape-not-mae-decides-ranking` shows the best-MAE arm ranked 16th of 28.

  Falsifier for the whole location direction: every family showing gam_eff <= 0, or gam_eff > 0 with
  cos so low that the orthogonal component dominates, would say the location cannot be moved toward
  truth by any reshaping of the consuming functional, and Phase I's open question is answered NO on
  the functional side without spending an end-to-end run.
"""
from __future__ import annotations

import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I           # noqa: E402
from s24 import stats_lib as ST           # noqa: E402
from core import predict as dgm           # noqa: E402

CEN = np.asarray(dgm.CENTRES, float)
SHELLS = ((2, 2), (3, 3), (4, 5), (6, 8), (9, 1000))
SHELL_NAMES = ("2-2", "3-3", "4-5", "6-8", "9-15")

QGRID = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
DGRID = (0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 40.0)
TGRID = (0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60)
EGRID = (0.0, 0.25, 0.5, 0.75, 1.0)
OUT = os.path.join(RES, "loc.json")

#: the deployed grid, taken from a genuine Distogram at run time -- never hard-coded here
_GRID = None


# ----------------------------------------------------------------- kernels and the argmin
def _kern_power(U, q):
    return np.abs(U) ** q


def _kern_trunc(U, d):
    return np.minimum(np.abs(U), d)


def _kern_quant(U, tau):
    return 2.0 * (tau * np.maximum(U, 0.0) + (1.0 - tau) * np.maximum(-U, 0.0))


def _argmin_grid(R, grid):
    """Location = mean grid value over the ARGMIN SET, not the first index numpy returns.

    `consensus-is-the-only-in-band-discriminator` (the methodological trap): np.argmin on a tied
    signal reads the array's own order and can manufacture a winner.  An L1 risk over a coarse
    17-bin posterior on a 0.05 A grid ties routinely.
    """
    m = R.min(1, keepdims=True)
    tie = R <= m + 1e-12
    return (tie * grid[None, :]).sum(1) / tie.sum(1)


# ----------------------------------------------------------------- metric projection (family D)
def _edm_project(M, w2, chain=3.81, iters=300):
    """Nearest realisable Ca distance field: classical MDS to rank 3, then deterministic SMACOF.

    M is the full (n,n) location matrix with the |i-j|=1 entries set to the ideal Ca-Ca 3.81 A (the
    score does not predict them; they are geometry, not prediction).  w2 is the (n,n) weight.
    No randomness anywhere: the init is the classical-MDS solution, so the arm is reproducible.
    """
    n = len(M)
    D2 = M ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    G = -0.5 * J @ D2 @ J
    ev, V = np.linalg.eigh((G + G.T) / 2.0)
    k = np.argsort(ev)[::-1][:3]
    X = V[:, k] * np.sqrt(np.maximum(ev[k], 0.0))[None, :]
    W = w2.copy()
    np.fill_diagonal(W, 0.0)
    s = W.sum(1)
    Vm = np.diag(s) - W
    Vp = np.linalg.pinv(Vm)
    for _ in range(iters):
        Dx = np.linalg.norm(X[:, None] - X[None], axis=-1)
        with np.errstate(divide="ignore", invalid="ignore"):
            B = np.where(Dx > 1e-9, -W * M / np.maximum(Dx, 1e-9), 0.0)
        np.fill_diagonal(B, 0.0)
        np.fill_diagonal(B, -B.sum(1))
        Xn = Vp @ (B @ X)
        if np.abs(Xn - X).max() < 1e-9:
            X = Xn
            break
        X = Xn
    return np.linalg.norm(X[:, None] - X[None], axis=-1)


# ----------------------------------------------------------------- per-target work
def _shell_of(sep):
    s = np.zeros(len(sep), int)
    for k, (a, b) in enumerate(SHELLS):
        s[(sep >= a) & (sep <= b)] = k
    return s


def _agg(dL, dT, w):
    """<dL,dT>, |dT|^2, |dL|^2 -- unweighted and w-weighted.  gam_eff/cos are formed in report()."""
    return dict(ip=float((dL * dT).sum()), tt=float((dT * dT).sum()), ll=float((dL * dL).sum()),
                ipw=float((w * dL * dT).sum()), ttw=float((w * dT * dT).sum()),
                llw=float((w * dL * dL).sum()),
                signed=float(dL.mean()), mae=float(np.abs(dT - dL).mean()))


def _sub(dL, dT, m):
    """The three sums on a sub-population, so gam_eff AND cosine can be formed per stratum."""
    if not m.any():
        return None
    return [float((dL[m] * dT[m]).sum()), float((dT[m] * dT[m]).sum()), float((dL[m] * dL[m]).sum()),
            int(m.sum())]


def _quant(L):
    """ADDITION 3: how far the realised minimiser sits from the nearest of the 17 bin CENTRES.

    Near zero means the arm has collapsed onto the binning and is measuring quantisation rather
    than the location statistic it claims to move.
    """
    d = np.abs(L[:, None] - CEN[None, :]).min(1)
    return [float(d.mean()), float((d < 0.026).mean())]


def run():
    global _GRID
    tg = I.targets()
    rows = []
    rng = np.random.default_rng(20250908)
    perm_target = rng.permutation(len(tg))          # D3 (ii): whose profile each target borrows
    print("s25/loc.py  LOCATION DIAGNOSTIC.  n=%d.  No RMSD is computed here." % len(tg), flush=True)

    #: PASS 1 -- per-pair native/prior/pool signed errors, needed for the family-C profile and D1
    cache = []
    for t in tg:
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        dg = I.distogram(pdb)
        prob = np.asarray(dg["prob"], float)
        pi, pj = np.asarray(dg["i"]), np.asarray(dg["j"])
        d0 = dgm.Distogram(u["seq"], prob, pi, pj)
        #: GATE -- the shipped artefact must be reproduced by this module's own kernel path, bit-exact
        if _GRID is None:
            _GRID = np.asarray(d0.grid, float)
        U = _GRID[:, None] - CEN[None, :]
        mine = np.ascontiguousarray((prob[:, None, :] * _kern_power(U, 1.0)[None]).sum(2)
                                    * np.asarray(d0.w, float)[:, None], dtype=np.float32)
        if not np.array_equal(mine, np.asarray(d0._risk)):
            raise SystemExit("KERNEL PATH IS NOT BIT-EXACT on %s (max %g). Stopping rather than "
                             "approximating the shipped loss silently."
                             % (pdb, float(np.abs(mine - np.asarray(d0._risk)).max())))
        idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]
        Dt = I.pair_dists(np.asarray(u["nat_ca"], float)[None], pi, pj)[0]
        Dp = I.pair_dists(W, pi, pj)                                     # (500, npairs)
        #: ADDITION 4 -- the offset of the candidates that actually REACH THE READOUT
        top = np.argsort(I.shipped_score(dg, Dp), kind="stable")[:75]
        cache.append(dict(pdb=pdb, n=int(u["n"]), fold=int(u["fold"]), seq=u["seq"],
                          prob=prob, i=pi, j=pj, w=np.asarray(d0.w, float), Dt=Dt,
                          pool_mean=Dp.mean(0), top_mean=Dp[top].mean(0),
                          sd=np.asarray(d0.sd, float),
                          expected=np.asarray(d0.expected, float)))
        del u, W, Dp
    print("  pass 1 complete: kernel path bit-exact on all %d targets." % len(cache), flush=True)

    #: the family-C per-shell offset profile.  ORACLE and IN-SAMPLE here -- this is a diagnostic; the
    #: end-to-end module fits it on TRAINING FOLDS ONLY.
    allsep = np.concatenate([c["j"] - c["i"] for c in cache])
    allsh = _shell_of(allsep)
    prof = np.zeros(len(SHELLS))
    for k in range(len(SHELLS)):
        pass  # filled below once L_med exists

    #: PASS 2 -- locations
    Lmed_all = []
    for c in cache:
        prob = c["prob"]
        U = _GRID[:, None] - CEN[None, :]
        R = prob @ _kern_power(U, 1.0).T
        Lmed_all.append(_argmin_grid(R, _GRID))
    for k in range(len(SHELLS)):
        m = allsh == k
        prof[k] = float((np.concatenate(Lmed_all)[m] - np.concatenate([c["Dt"] for c in cache])[m]).mean())
    print("  family-C per-shell median offset profile (ORACLE, in-sample): %s"
          % np.round(prof, 4).tolist(), flush=True)

    for ci, c in enumerate(cache):
        prob = c["prob"]; Dt = c["Dt"]; w = c["w"]
        sep = (c["j"] - c["i"]).astype(int)
        sh = _shell_of(sep)
        U = _GRID[:, None] - CEN[None, :]
        Lmed = Lmed_all[ci]
        dT = Dt - Lmed
        row = {"pdb": c["pdb"], "n": c["n"], "fold": c["fold"], "npairs": int(len(Dt))}
        row["arms"] = {}

        #: the multimodal criterion, computed BEFORE the arms so every arm can be stratified by it
        p = prob
        pk = np.zeros(p.shape, bool)
        pk[:, 1:-1] = (p[:, 1:-1] > p[:, :-2]) & (p[:, 1:-1] > p[:, 2:]) & (p[:, 1:-1] > 0.02)
        pk[:, 0] = (p[:, 0] > p[:, 1]) & (p[:, 0] > 0.02)
        pk[:, -1] = (p[:, -1] > p[:, -2]) & (p[:, -1] > 0.02)
        mm = pk.sum(1) >= 2

        def add(tag, L):
            dL = L - Lmed
            z = _agg(dL, dT, w)
            z["shell"] = {SHELL_NAMES[k]: _sub(dL, dT, sh == k) for k in range(len(SHELLS))}
            z["mm"] = _sub(dL, dT, mm)
            z["uni"] = _sub(dL, dT, ~mm)
            z["quant"] = _quant(L)
            row["arms"][tag] = z

        add("IDENT", Lmed)
        add("MEAN", c["expected"])
        add("MODE", CEN[prob.argmax(1)])
        for q in QGRID:
            add("POWER%.2f" % q, _argmin_grid(prob @ _kern_power(U, q).T, _GRID))
        for dd in DGRID:
            add("TRUNC%.1f" % dd, _argmin_grid(prob @ _kern_trunc(U, dd).T, _GRID))
        for tau in TGRID:
            add("QUANT%.2f" % tau, _argmin_grid(prob @ _kern_quant(U, tau).T, _GRID))
        #: family C -- shift the location by alpha * (its shell's offset)
        for a in (0.0, 0.5, 1.0):
            add("SHELL%.2f" % a, Lmed - a * prof[sh])
        #: family D -- metric projection of the whole location field
        n = c["n"]
        M = np.zeros((n, n)); Wm = np.zeros((n, n))
        M[c["i"], c["j"]] = Lmed; M[c["j"], c["i"]] = Lmed
        Wm[c["i"], c["j"]] = w; Wm[c["j"], c["i"]] = w
        for k in range(n - 1):
            M[k, k + 1] = M[k + 1, k] = 3.81
            Wm[k, k + 1] = Wm[k + 1, k] = float(w.mean())
        P = _edm_project(M, Wm)
        Lproj = P[c["i"], c["j"]]
        for e in EGRID:
            add("METRIC%.2f" % e, (1.0 - e) * Lmed + e * Lproj)
        #: D3 controls, matched in the operator's space
        best = Lmed - 1.0 * prof[sh]
        add("CTRL_PERM", Lmed + rng.permutation(best - Lmed))
        oth = cache[perm_target[ci]]
        othsh = _shell_of((oth["j"] - oth["i"]).astype(int))
        othprof = np.zeros(len(SHELLS))
        for k in range(len(SHELLS)):
            m = othsh == k
            othprof[k] = float((Lmed_all[perm_target[ci]][m] - oth["Dt"][m]).mean()) if m.any() else 0.0
        add("CTRL_XSHELL", Lmed - othprof[sh])

        #: D1 -- the cancellation test.  Signed offset of the PRIOR and of the POOL, per shell.
        row["shell"] = {}
        for k in range(len(SHELLS)):
            m = sh == k
            if not m.any():
                continue
            row["shell"][SHELL_NAMES[k]] = dict(
                npairs=int(m.sum()),
                prior_med=float((Lmed[m] - Dt[m]).mean()),
                prior_mean=float((c["expected"][m] - Dt[m]).mean()),
                pool=float((c["pool_mean"][m] - Dt[m]).mean()),
                top75=float((c["top_mean"][m] - Dt[m]).mean()),
                sd=float(c["sd"][m].mean()))

        #: D2 -- the multimodal sub-population
        row["multimodal"] = dict(frac=float(mm.mean()), npairs=int(mm.sum()))
        if mm.any():
            nearmode = CEN[np.argmin(np.where(pk[mm], np.abs(CEN[None, :] - Dt[mm][:, None]),
                                              np.inf), 1)]
            row["multimodal"].update(
                err_med=float(np.abs(Lmed[mm] - Dt[mm]).mean()),
                err_nearest_mode=float(np.abs(nearmode - Dt[mm]).mean()),
                err_top_mode=float(np.abs(CEN[p[mm].argmax(1)] - Dt[mm]).mean()),
                unimodal_err_med=float(np.abs(Lmed[~mm] - Dt[~mm]).mean()) if (~mm).any() else None)
        rows.append(row)
        if (ci + 1) % 25 == 0:
            print("  %d/%d" % (ci + 1, len(cache)), flush=True)

    ST.save_atomic(OUT, {"rows": rows, "shell_profile": prof.tolist(),
                         "qgrid": list(QGRID), "dgrid": list(DGRID), "tgrid": list(TGRID),
                         "egrid": list(EGRID), "shells": list(SHELL_NAMES)},
                   complete_keys=("arms", "shell", "multimodal"), rows=rows,
                   n_expected=len(tg), module_file=__file__)
    report(rows)
    return rows


# ----------------------------------------------------------------- reporting
def report(rows=None):
    import json
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    names = sorted(rows[0]["arms"].keys())
    order = (["IDENT", "MEAN", "MODE"]
             + ["POWER%.2f" % q for q in QGRID] + ["TRUNC%.1f" % d for d in DGRID]
             + ["QUANT%.2f" % t for t in TGRID] + ["SHELL%.2f" % a for a in (0.0, 0.5, 1.0)]
             + ["METRIC%.2f" % e for e in EGRID] + ["CTRL_PERM", "CTRL_XSHELL"])
    order = [o for o in order if o in names]

    print("\n" + "=" * 100)
    print("s25/loc.py  LOCATION DIAGNOSTIC.  n = %d targets.  DISTANCE SPACE, min_sep=2." % len(rows))
    print("ORACLE DIAGNOSTIC -- natives read for diagnosis only.  NO RMSD IS COMPUTED IN THIS FILE.")
    print("=" * 100)
    print("\n  gam_eff = <dL,dT>/|dT|^2 with dL = L - L_med, dT = D_native - L_med.")
    print("  The prior ladder's gamma has cos = 1 BY CONSTRUCTION; these moves do not, so")
    print("  -2.1496 * gam_eff is an UPPER BOUND on the endpoint gain, never a prediction.\n")
    def _gc(a, key=None):
        """(gam_eff, cos, |dL|/|dT|) averaged PER TARGET.  key selects a stratum."""
        ge, co, rl = [], [], []
        for r in rows:
            z = r["arms"][a]
            v = ([z["ip"], z["tt"], z["ll"]] if key is None
                 else (z["shell"][key][:3] if key in SHELL_NAMES and z["shell"].get(key)
                       else (z[key][:3] if z.get(key) else None)))
            if v is None or v[1] <= 1e-12:
                continue
            ge.append(v[0] / v[1])
            co.append(v[0] / max(np.sqrt(v[2] * v[1]), 1e-12) if v[2] > 1e-18 else 0.0)
            rl.append(np.sqrt(v[2] / v[1]))
        if not ge:
            return None
        return np.mean(ge), np.mean(co), np.mean(rl)

    print("  %-13s %9s %9s %9s %9s %9s %10s %8s %7s" %
          ("arm", "gam_eff", "cos", "|dL|/|dT|", "gam_w", "MAE", "bound A", "d2centre", "on_ctr"))
    for a in order:
        g, co, rl = _gc(a)
        gw = np.mean([r["arms"][a]["ipw"] / max(r["arms"][a]["ttw"], 1e-12) for r in rows])
        mae = np.mean([r["arms"][a]["mae"] for r in rows])
        qd = np.mean([r["arms"][a]["quant"][0] for r in rows])
        qf = np.mean([r["arms"][a]["quant"][1] for r in rows])
        print("  %-13s %+9.4f %+9.4f %9.4f %+9.4f %9.4f %10s %8.3f %7.3f"
              % (a, g, co, rl, gw, mae, "%+.4f" % (-2.1496 * g), qd, qf))
    print("\n  d2centre / on_ctr are ADDITION 3, the quantisation diagnostic: mean distance from the")
    print("  realised minimiser to the nearest of the 17 IRREGULAR bin centres, and the fraction")
    print("  sitting on one.  An arm whose minimiser has collapsed onto the centres is measuring")
    print("  BIN QUANTISATION, not the location statistic it claims to move, and the effect is")
    print("  separation-dependent because the bins are (0.5 A apart at 4 A, 3.5 A at 21 A).")

    print("\n  --- gam_eff / cos PER SHELL (ADDITION 1).  cos is the honest discount on the bound.")
    hdr = "  %-13s" % "arm" + "".join("%18s" % s for s in SHELL_NAMES)
    print(hdr)
    for a in order:
        cells = []
        for s in SHELL_NAMES:
            v = _gc(a, s)
            cells.append("       --       " if v is None else "%+8.4f/%+7.3f" % (v[0], v[1]))
        print("  %-13s" % a + "".join("%18s" % c for c in cells))

    print("\n  --- gam_eff / cos STRATIFIED BY MULTIMODALITY (ADDITION 2, pre-registered before the")
    print("      run).  MM = posterior carries >= 2 peaks above 0.02, calib.py's exact criterion.")
    print("  %-13s %20s %20s" % ("arm", "MULTIMODAL", "UNIMODAL"))
    for a in order:
        vm, vu = _gc(a, "mm"), _gc(a, "uni")
        f = lambda v: "      --      " if v is None else "%+9.4f/%+7.3f" % (v[0], v[1])  # noqa: E731
        print("  %-13s %20s %20s" % (a, f(vm), f(vu)))

    print("\n  --- D1  THE CANCELLATION TEST.  Signed offset (estimate - native), per shell.")
    print("  PRE-REGISTERED PREDICTION: pool matches prior to within ~0.1 A at every shell, which")
    print("  would mean the offset is shared by the candidates and shifting only the prior breaks")
    print("  a cancellation.  s24 Workstream C: distogram +0.4062, real protein windows +0.4419.\n")
    print("  %-8s %8s %11s %11s %13s %11s %11s %11s" %
          ("shell", "pairs", "prior_med", "prior_mean", "POOL(K=500)", "TOP-75", "pri-pool",
           "pri-top75"))
    for s in SHELL_NAMES:
        v = [r["shell"][s] for r in rows if s in r["shell"]]
        if not v:
            continue
        pm = np.mean([x["prior_med"] for x in v]); po = np.mean([x["pool"] for x in v])
        tp = np.mean([x["top75"] for x in v])
        print("  %-8s %8d %+11.4f %+11.4f %+13.4f %+11.4f %+11.4f %+11.4f"
              % (s, sum(x["npairs"] for x in v), pm,
                 np.mean([x["prior_mean"] for x in v]), po, tp, pm - po, pm - tp))
    am = lambda k: np.mean([np.mean([x[k] for x in r["shell"].values()]) for r in rows])  # noqa: E731
    print("  %-8s %8s %+11.4f %11s %+13.4f %+11.4f %+11.4f %+11.4f"
          % ("ALL", "", am("prior_med"), "", am("pool"), am("top75"),
             am("prior_med") - am("pool"), am("prior_med") - am("top75")))

    print("\n  --- D2  THE MULTIMODAL SUB-POPULATION (>= 2 peaks above 0.02).")
    f = np.mean([r["multimodal"]["frac"] for r in rows])
    em = np.mean([r["multimodal"]["err_med"] for r in rows if "err_med" in r["multimodal"]])
    en = np.mean([r["multimodal"]["err_nearest_mode"] for r in rows if "err_med" in r["multimodal"]])
    et = np.mean([r["multimodal"]["err_top_mode"] for r in rows if "err_med" in r["multimodal"]])
    eu = np.mean([r["multimodal"]["unimodal_err_med"] for r in rows
                  if r["multimodal"].get("unimodal_err_med") is not None])
    print("    multimodal fraction %.3f   (s25 L1 published 0.241)" % f)
    print("    on those pairs:  |median - native| %.4f    |NEAREST mode - native| %.4f  (ORACLE"
          " choice of mode)" % (em, en))
    print("                     |TOP mode - native| %.4f   <- the achievable mode-seeking arm" % et)
    print("    on unimodal pairs: |median - native| %.4f" % eu)

    print("\n  --- D3  ZERO-INFORMATION CONTROLS, matched in the operator's space.")
    print("    CTRL_PERM   = the SHELL1.00 move permuted within the target (same magnitudes, no")
    print("                  pair correspondence).")
    print("    CTRL_XSHELL = another target's per-shell profile (same shape, no target"
          " correspondence).")
    print("    Read them against SHELL1.00 in the table above.\n")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
