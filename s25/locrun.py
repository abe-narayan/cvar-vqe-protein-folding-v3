"""s25/locrun.py -- CAN THE SCORE'S LOCATION BE MOVED, AND DOES THE ENDPOINT RESPOND?  END TO END.

**PRIMARY ENDPOINT: mean full-chain Ca-RMSD, 126 cluster-disjoint dev targets, POINT CLOUD basis,
shipped uniform top-75 coordinate average in the medoid frame.  Incumbent 3.0483 A.**
The pool, the readout and the per-pair weight `w` are IDENTICAL in every arm.  The ONLY thing that
varies is the risk table, and therefore WHICH 75 candidates are selected.

Pre-registered in `s25/PREREG_RMSD.md`; run 1 (`s25/loc.py`) is the diagnostic that chose these
families and its results were sent to the coordinator before this ran.

------------------------------------------------------------------------------------------------
WHAT RUN 1 ESTABLISHED, AND WHAT IT CHANGED

gam_eff = <dL,dT>/|dT|^2, the fraction of the way to truth a location move travels, in the prior
ladder's own currency.  -2.1496*gam_eff is an UPPER BOUND on the endpoint gain, never a prediction,
because the ladder's gamma has cos = 1 by construction and none of these moves do.

    arm            gam_eff      cos     MAE     on_bin_centre
    IDENT          +0.0000   +0.000   2.3968   1.000
    MODE / TRUNC   -0.0189   -0.068   2.4906   1.000     mode-seeking is WRONG-SIGNED
    SHELL1.00      +0.0136   +0.094   2.3993   0.189     its own permuted control gets +0.0108
    MEAN / POWER2  +0.0582   +0.253   2.3386   0.063
    METRIC1.00     +0.2497   +0.495   2.1756   0.067     largest, and predicted to convert worst

**THE INCUMBENT'S TARGET IS QUANTISED.**  The L1 minimiser of a discrete distribution sits on an
atom, so the shipped score's effective per-pair target takes only **17 distinct values** (4.0, 4.75,
... 21.0, 25.0), spaced 3.5 A apart at the top of the range -- up to +/-1.75 A of pure location
error at long separation from binning alone.  Across 30 arms in run 1 the SIGN of gam_eff tracks
whether the arm de-quantises: every arm with gam_eff > 0 has on_ctr <= 0.19, every arm with
gam_eff <= 0 has on_ctr = 1.000.  Hence family F below, which is the primary.

------------------------------------------------------------------------------------------------
THE SIX FAMILIES.  ONE GLOBAL PARAMETER EACH, INCUMBENT AN EXACT INTERIOR POINT.

Two implementation styles, both of which reproduce the shipped table BIT-EXACTLY at the identity
parameter (asserted per target with `np.array_equal`, and additionally cross-checked per target
against `s25/results/temper.json`'s independently computed f=1 column):

  (a) KERNEL REPLACEMENT -- the loss shape changes, the posterior does not.
      A POWER(q)    rho = |u|^q             q=1 IDENTITY.  q<1 mode-seeking, q=2 = L2 = mean.
      B TRUNC(d)    rho = min(|u|, d)       d=40 >= max|grid-C| = 35.95 IDENTITY.  d->0 mode-seeking.
      C QUANT(tau)  DROPPED BEFORE RUNNING.  A pure location shift, self-scaled by the posterior's
                    own spread.  Closed for free: run 1's D1 predicted it would fail end-to-end
                    despite a positive gam_eff, and the audit lane's A3/A8 measured exactly that --
                    a global translation is SYMMETRIC and null (-0.30 A costs +0.0206, +0.30 A costs
                    +0.0191).  No compute spent.
      F DEQUANT(b)  each bin's point mass spread UNIFORMLY over ITS OWN support, half-width
                    b*HALF[c], CENTRE UNMOVED.  b=0 IDENTITY; b=1 is the histogram taken literally.
                    rho_c(u) = |u| if |u|>=h_c else (u^2+h_c^2)/(2 h_c).
                    The risk derivative becomes sum_c p_c clip((t-C_c)/h_c, -1, 1), continuous and
                    strictly increasing, so **the minimiser is continuous instead of quantised**.
                    **NO MASS CROSSES A BIN AND NO CENTRE MOVES**, which is what separates it from
                    the widening arms: those move mass BETWEEN bins.  b > 1 does cross boundaries
                    and is a smoothing rather than the histogram's own resolution -- b = 1 is the
                    principled point and b = 1.5, 2.0 are carried as an extension, labelled.
      F'DEQUANTI(b) IDENTICAL except the two CATCH-ALL bins are left as point masses.  Mass in bin 0
                    means "< 4.5 A" and in bin 16 means ">= 23 A"; neither is uniform over a window,
                    and they sit exactly where L6 says the quantisation error is largest.  DEQUANTI
                    assumes NOTHING about them and is **the conservative PRIMARY**; DEQUANT is
                    reported beside it so the reader can see what the assumption is worth.

  (b) PER-PAIR TRANSLATION -- the location moves, the posterior's SHAPE is preserved exactly.
      D SHELL(a)    DROPPED BEFORE RUNNING, and CLOSED BY TWO INDEPENDENT ROUTES.  (i) Run 1: at
                    a=1 with an ORACLE in-sample profile it travels gam_eff 0.0136, its own
                    magnitude-matched permuted control gets 0.0108 (79% of it), and ANOTHER
                    TARGET'S profile scores higher than the target's own.  (ii) The audit lane's
                    A8: a separation-graded shift correlating r=+0.977 with L1's own signed-error
                    profile is null on the endpoint (+0.0270, 0.41x MDE) -- and it was already
                    running, unlabelled, inside temper.py's SD arm.  No compute spent.
      E METRIC(e)   target moved to (1-e)*L_med + e*L_proj, L_proj = the deterministic rank-3 EDM
                    projection of the whole location field (classical MDS init, weighted SMACOF, no
                    randomness).  e=0 IDENTITY.  **RETEST OF A CLOSED DIRECTION**, labelled as such:
                    `core/predict.py`'s own `realize` docstring records "projecting the matrix and
                    scoring against the projection measurably hurts", and `better-matrix-worse-
                    ranking` item 1 prices it at MAE 2.13->1.96 with in-band rho 0.379->0.334.  Both
                    were measured on the decoy bank the project has since RETRACTED (`decoy-bank-
                    not-a-pool-proxy`) and neither had an interpolation parameter, which is what
                    makes e=0 as an exact interior point a different test.

**`w` IS NEVER TOUCHED BY ANY ARM, AND THAT MATTERS MORE THAN I THOUGHT.**  The audit lane has since
established that `_score_weights()` returns `(ones, 1.0)` because `score_weights.json` is
deliberately absent, so `w = shell/(sd+0.5)^g` reduces to `w = 1/(sd+0.5)`, normalised to mean 1 --
**the only per-pair weight in the shipped score is a pure function of the posterior's WIDTH.**  No
arm here modifies `prob` or `CENTRES`, so `sd` and therefore `w` are BIT-IDENTICAL on both sides of
every contrast in this file.  Whatever the audit concludes about the width channel cannot confound
anything measured here.  NOTE THE ALTERNATIVE NOT TAKEN: a genuinely de-quantised representation has
a larger effective sd (by h^2/3 within each bin) and would earn a different `w`; recomputing it is
a pair-reweighting-by-confidence intervention, which is closed (s24 L6), so `w` stays frozen.

------------------------------------------------------------------------------------------------
THE FRAME THIS MODULE WAS DESIGNED UNDER HAS BEEN RETRACTED.  WHAT SURVIVES, AND WHY THESE ARMS DO.

s25 L7 retracts L2's "the ranking responds to LOCATION and not to WIDTH".  The unconfounded
experiment -- mean-preserving widening against pure translation, n=126 -- is FLAT ON BOTH AXES:

    pure WIDTH     f=1.3/1.5/2.0/3.0     +0.0139/+0.0208/+0.0203/+0.0253   0.51/0.56/0.43/0.38x MDE
    pure LOCATION  d=-0.9/-0.6/-0.3/+0.3 +0.0617/+0.0322/+0.0206/+0.0191   0.58/0.40/0.40/0.44x MDE

All eight NOT MEASURED, and the response to location is SYMMETRIC about the shipped posterior even
though that posterior is biased -0.41 A against the natives -- moving the prior TOWARD truth is
worth the same as moving it away.  **So the prior for everything below is much lower than it was
when the families were chosen, and I am saying so before the numbers exist.**

THE ONE DISTINCTION KEEPING THESE ARMS ALIVE, STATED PLAINLY BECAUSE IT IS DOING A LOT OF WORK:
**the audit's location arm is a RIGID TRANSLATION -- one constant added to every pair.  Every arm
here is an ADAPTIVE PER-PAIR change whose size and sign are set by each pair's own posterior shape.**
A rigid shift being flat does not imply an adaptive one is.  That is the whole of the remaining
case, and if these arms are also flat the two results together close the functional side.

AND THE SCOPING LINE I AM ASKED TO PLACE EACH ARM AGAINST, BEFORE SEEING ANY ENDPOINT.  Reaching
gamma needs genuinely BETTER location -- new information -- not a RE-READING of the location already
in the posterior:

    POWER, TRUNC          RE-READINGS.  Median -> mean/mode/quantile of the SAME posterior.  No new
                          information of any kind.  On the audit's scoping these should not move the
                          endpoint, and I expect them not to.
    QUANT, SHELL          DROPPED, see above.  SHELL is DOUBLY CLOSED --
                          dead against its own permuted control in run 1, and independently
                          reproduced as the unlabelled per-separation confound inside temper.py's
                          SD arm.  Carried at a=1 only, as a two-route confirmation, not a test.
    METRIC                NEW information -- the metric constraint is not in the posterior at all.
                          But see P1: I predict the ranking cannot see it.
    DEQUANT / DEQUANTI    **On the re-reading side of the line, strictly speaking**: the bin support
                          is not information about THIS target.  What it removes is a REPRESENTATION
                          defect -- a 17-atom discretisation of a continuous quantity -- rather than
                          an estimation error.  That is a real distinction but it is not new
                          information, so my honest prior is LOW, and lower than it was before L7.

THE SCALE CONFOUND IN THE SHAPE FAMILIES, NAMED AND CONTROLLED.  `|u|^q`, `min(|u|,d)` and the
DEQUANT kernel change the location statistic AND the relative weight of one pair against another
(a pair with large typical |u| gains weight at q=2 and loses it at q<1).  Every kernel family is
therefore run TWICE: raw, and SCALE-MATCHED -- each pair's risk row rescaled so its mean over the
grid equals the identity row's, which is exactly 1.0 at the identity parameter and isolates SHAPE
from WEIGHT.  Same device as `priorladder`'s MASSFIXW and `temper`'s SDFIXW.

------------------------------------------------------------------------------------------------
TWO PREDICTIONS REGISTERED WITH THE COORDINATOR BEFORE THIS RAN

P1  **METRIC will show the LARGEST gam_eff and the SMALLEST endpoint movement per unit gam_eff of
    any family.**  Every candidate in the pool is itself realisable, so the component the projection
    removes lies orthogonal to the manifold the candidates live on and adds a near-CONSTANT offset
    to every candidate's score.  gam_eff is measured against the native, which is a point IN the
    set, so it counts that component; the ranking cannot see it.  A shared-referent floor applies
    too: orthogonal projection onto a set containing the native gives gam_eff = (|dL|/|dT|)^2 =
    0.202 by identity alone, so METRIC's genuine excess is 0.048, not 0.250.
P2  `amp = (dL' Sigma dL) / (|dL|^2 * tr(Sigma)/npairs)`, with Sigma the POOL's own distance
    covariance, measures how much of a move points along directions the candidates actually vary in.
    **amp << 1 for METRIC; amp at or above 1 for MEAN/POWER2 and DEQUANT.**  If P1 holds and P2
    explains it, `gam_eff` is established as insufficient on its own and `amp` is the missing
    discount -- which is a transferable methodological result independent of any endpoint move.

RULE 0 -- THE SIX APPROVED FORKS, PLUS A SEVENTH THE COORDINATOR REQUIRED BEFORE DEQUANT COULD RUN.

    BIN SUPPORT    DECLARED the bin edges read VERBATIM from `core.predict.BIN_EDGES` -- the same
    (fork 7)       array `_labels` digitises the training distances against, so they are published
                   in the source and are not assumed -- with `CENTRES[1:-1]` ASSERTED equal to their
                   exact midpoints at import time.  `HALF[c]` is the largest symmetric window about
                   the published centre that does not cross a published boundary: 0.25-2.0 A on the
                   interior, and 0.5 / 2.0 on the two catch-alls, where the centre is NOT a midpoint
                   (CENTRES[0] = 4.5-0.5, CENTRES[16] = 23+2.0).
                   NOT TAKEN inferring the edges as midpoints between adjacent CENTRES, which is
                   what a reader of the artefact alone would have to do -- the artefact publishes
                   `centres` only -- and which gets both catch-alls wrong by construction.
                   NOT TAKEN assigning the catch-alls a support at all in the PRIMARY arm: DEQUANTI
                   leaves them as point masses, and DEQUANT (which does assign one) is reported
                   beside it so the assumption's value is visible rather than baked in.

    functional     DECLARED the shipped Bayes-risk form with `w` and `grid` taken from a genuine
                   `core.predict.Distogram` and ONLY the loss kernel (a) or the evaluation point (b)
                   replaced; the identity parameter asserted `np.array_equal` to the object's own
                   `_risk` PER TARGET, and the identity arm's per-target RMSD asserted equal to
                   `temper.json`'s independently computed f=1 column.  NOT TAKEN re-implementing the
                   risk table or `w`; NOT TAKEN modifying `w`.
    basis          DECLARED point cloud, medoid frame, on BOTH sides of every contrast.  NOT TAKEN
                   the built-chain basis (0.156 A of pure operator choice).
    readout        DECLARED the shipped uniform top-75 coordinate average, unchanged in every arm.
                   NOT TAKEN re-optimising m, weighting members, or clustering (s23 L5: uniform is
                   the in-sample optimum and removing geometric outliers COSTS +0.142 A).
    normalisation  DECLARED every kernel family run raw AND scale-matched, side by side.  NOT TAKEN
                   choosing one and reporting it as "the" result.
    null           DECLARED each family's identity parameter, reproducing the incumbent bit-for-bit
                   through this module's own path, so the null is a self-check on the scoring path
                   rather than a comparison against a pinned constant.  NOT TAKEN comparing to
                   3.0483 alone.
    THE LABEL      DECLARED held-out Ca-RMSD AND gam_eff AND amp AND held-out z_sd, reported
                   together, plus the endpoint delta stratified by the target's multimodal fraction
                   (pre-registered by the coordinator before run 1).  NOT TAKEN RMSD alone; NOT
                   TAKEN any binarised win rate, which cannot diagnose a best-of-K arm.

  H  the score's per-pair location can be moved toward truth using no new information -- by
     de-quantising the binned posterior, by consuming it with a different risk functional, or by
     imposing metric realisability -- and the endpoint responds.
  Falsifier  every family's nested-CV arm failing to beat its own identity parameter past its own
     MDE with a fold-clustered CI excluding zero.  That answers Phase I's open question NO on the
     functional side, which is a result and will be reported as one.
"""
from __future__ import annotations

import json
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
from s24 import residlib as RL            # noqa: E402
from s25 import loc as LOC                # noqa: E402
from core import predict as dgm           # noqa: E402

CEN = np.asarray(dgm.CENTRES, float)
BE = np.asarray(dgm.BIN_EDGES, float)

#: FORK 7 -- THE BIN SUPPORT CONVENTION.  See the docstring.  The edges are NOT assumed: they are
#: read verbatim from `core.predict.BIN_EDGES`, the same array `features`/`_labels` digitises
#: against, and the interior centres are ASSERTED to be their exact midpoints.  HALF[c] is the
#: largest symmetric window about the published centre that does not cross a published boundary.
HALF = np.empty(len(CEN))
HALF[0] = BE[0] - CEN[0]                                  # 0.5, the catch-all "< 4.5"
HALF[1:-1] = np.minimum(CEN[1:-1] - BE[:-1], BE[1:] - CEN[1:-1])
HALF[-1] = CEN[-1] - BE[-1]                               # 2.0, the catch-all ">= 23"
assert np.array_equal(CEN[1:-1], 0.5 * (BE[1:] + BE[:-1])) and (HALF > 0).all()
INTERIOR = np.zeros(len(CEN), bool)
INTERIOR[1:-1] = True                                     # the bins whose support is EXACT
TOPM = 75
OUT = os.path.join(RES, "locrun.json")
LOCK = os.path.join(RES, "LOCK_TRAIN")

QGRID = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)                  # identity q = 1.0   (index 2)
DGRID = (1.0, 2.0, 3.0, 5.0, 8.0, 40.0)                   # identity d = 40.0  (index 5)
TGRID = (0.40, 0.45, 0.50, 0.55, 0.60)                    # identity tau = 0.5 (index 2)
BGRID = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0)             # identity b = 0.0   (index 0)
AGRID = (0.0, 1.0)                                        # identity a = 0.0   (index 0)
EGRID = (0.0, 0.25, 0.5, 0.75, 1.0)                       # identity e = 0.0   (index 0)
IDENT = {"POWER": 2, "TRUNC": 5, "DEQUANT": 0, "DEQUANTI": 0, "METRIC": 0}
GRIDS = {"POWER": QGRID, "TRUNC": DGRID, "DEQUANT": BGRID, "DEQUANTI": BGRID,
         "METRIC": EGRID}
#: REPRIORITISED after s25 L7.  METRIC first (the only arm that ADDS information across pairs),
#: DEQUANTI/DEQUANT second (pair-specific, removes a real discretisation error), POWER third and
#: honestly a RE-READING.  QUANT and SHELL are DROPPED -- closed for free by the D1 prediction and
#: the audit lane's A3/A8, and no compute is spent on them.
FAMS = ("METRIC", "DEQUANTI", "DEQUANT", "POWER", "TRUNC")
KFAMS = ("DEQUANTI", "DEQUANT", "POWER", "TRUNC")
SCALED = ("DEQUANTI", "DEQUANT", "POWER", "TRUNC")        # the kernel families


# ---------------------------------------------------------------- kernels, K is (760, 17)
def _k_power(G, q):
    return np.abs(G[:, None] - CEN[None, :]) ** q


def _k_trunc(G, d):
    return np.minimum(np.abs(G[:, None] - CEN[None, :]), d)


def _k_quant(G, tau):
    U = G[:, None] - CEN[None, :]
    return 2.0 * (tau * np.maximum(U, 0.0) + (1.0 - tau) * np.maximum(-U, 0.0))


def _k_dequant(G, b, interior_only=False):
    """Each bin's mass spread uniformly over its own support; centre unmoved.  b=0 -> |u| exactly.

    `interior_only` leaves the two CATCH-ALL bins (mass there means "<= 4.5" and ">= 23", not
    "uniform over a window") as point masses.  That arm assumes NOTHING about their support and is
    the conservative primary; the full arm is reported beside it.
    """
    U = G[:, None] - CEN[None, :]
    A = np.abs(U)
    if b <= 0.0:
        return A
    h = b * HALF[None, :].copy()
    if interior_only:
        h = np.where(INTERIOR[None, :], h, 0.0)
    return np.where(A >= h, A, (U * U + h * h) / (2.0 * np.maximum(h, 1e-300)))


KERN = {"POWER": _k_power, "TRUNC": _k_trunc, "QUANT": _k_quant, "DEQUANT": _k_dequant,
        "DEQUANTI": lambda G, b: _k_dequant(G, b, interior_only=True)}


def _lock():
    fd = os.open(LOCK, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, ("%d locrun.py\n" % os.getpid()).encode())
    os.close(fd)


def _unlock():
    try:
        os.remove(LOCK)
    except OSError:
        pass


def run():
    tg = I.targets()
    tmp = json.load(open(os.path.join(RES, "temper.json")))["rows"]
    ref = {r["pdb"]: float(r["SD"][0]) for r in tmp}          # the independent f=1 column
    lrows = json.load(open(os.path.join(RES, "loc.json")))["rows"]
    lmm = {r["pdb"]: float(r["multimodal"]["frac"]) for r in lrows}
    #: LEAVE-ONE-FOLD-OUT per-shell offset profiles for family D, built from run 1's per-target
    #: per-shell offsets.  A target's own fold NEVER contributes to the profile applied to it.
    folds = sorted({r["fold"] for r in lrows})
    prof = {}
    for q in folds:
        acc = {}
        for r in lrows:
            if r["fold"] == q:
                continue
            for s, v in r["shell"].items():
                a, b = acc.setdefault(s, [0.0, 0])
                acc[s] = [a + v["prior_med"] * v["npairs"], b + v["npairs"]]
        prof[q] = {s: (a / b) for s, (a, b) in acc.items()}
    print("LOFO shell profiles (family D):")
    for q in folds:
        print("  fold %d held out: %s" % (q, {k: round(v, 4) for k, v in prof[q].items()}))

    rows = []
    for ci, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        seq = u["seq"]; nat = np.asarray(u["nat_ca"], float); n = int(u["n"])
        pi, pj = I.pair_index(n)
        dg = I.distogram(pdb, seq, u["fold"])
        prob = np.asarray(dg["prob"], float)
        d0 = dgm.Distogram(seq, prob, np.asarray(dg["i"]), np.asarray(dg["j"]))
        G = np.asarray(d0.grid, float)
        w = np.asarray(d0.w, float)
        shipped = np.asarray(d0._risk)

        #: GATE 1 -- this module's own kernel path must reproduce the shipped table BIT-EXACTLY
        base = np.einsum("pc,tc->pt", prob, _k_power(G, 1.0))
        if not np.array_equal(np.ascontiguousarray(base * w[:, None], dtype=np.float32), shipped):
            raise SystemExit("KERNEL PATH NOT BIT-EXACT on %s. Stopping." % pdb)
        base_mean = base.mean(1)                             # for the scale-matched variants

        Wp = np.asarray(u["W"], float)[I.pool_idx(u)]
        D = I.pair_dists(Wp, pi, pj)
        Dt = I.pair_dists(nat[None], pi, pj)[0]
        X = D - D.mean(0, keepdims=True)                     # for the amp diagnostic
        mvar = float((X * X).sum() / (len(D) - 1) / D.shape[1])
        sep = (pj - pi).astype(int)
        sh = LOC._shell_of(sep)
        Lmed = LOC._argmin_grid(base, G)
        dT = Dt - Lmed

        #: family E -- the deterministic rank-3 EDM projection of the location field
        M = np.zeros((n, n)); Wm = np.zeros((n, n))
        M[pi, pj] = Lmed; M[pj, pi] = Lmed
        Wm[pi, pj] = w; Wm[pj, pi] = w
        for k in range(n - 1):
            M[k, k + 1] = M[k + 1, k] = 3.81
            Wm[k, k + 1] = Wm[k + 1, k] = float(w.mean())
        Lproj = LOC._edm_project(M, Wm)[pi, pj]

        row = {"pdb": pdb, "n": n, "fold": int(u["fold"]), "mmfrac": lmm[pdb],
               "Emean": float(np.asarray(d0.expected, float).mean()),
               "Lmed": float(Lmed.mean()), "Lproj": float(Lproj.mean())}

        def emit(rk, tag):
            """Score, select top-75, shipped readout, Ca-RMSD.  Also stores gam_eff and amp.

            The table is rounded to float32 and read back as float64 -- exactly what
            `core.predict.Distogram` stores and what `temper.py`/`priorladder.py` consume -- so the
            identity arms are bit-comparable to those modules rather than merely close.
            """
            sc = np.asarray(I.shipped_score(
                {"grid": G, "risk": np.asarray(rk.astype(np.float32), float)}, D), float)
            sel = Wp[np.argsort(sc, kind="stable")[:TOPM]]
            row[tag] = float(I.ca_rmsd(RL.readout(sel), nat))
            L = LOC._argmin_grid(rk / np.maximum(w[:, None], 1e-12), G)
            dL = L - Lmed
            ll = float((dL * dL).sum())
            row["g_" + tag] = float((dL * dT).sum() / max((dT * dT).sum(), 1e-12))
            row["a_" + tag] = (float((X @ dL) @ (X @ dL) / (len(D) - 1) / max(ll, 1e-12) / mvar)
                               if ll > 1e-12 else 0.0)
            row["c_" + tag] = float((dL * dT).sum() / max(np.sqrt(ll * float((dT * dT).sum())), 1e-12))
            row["L_" + tag] = float(L.mean())
            row["onc_" + tag] = float((np.abs(L[:, None] - CEN[None, :]).min(1) < 0.026).mean())

        var0 = (prob * (CEN[None] - (prob * CEN[None]).sum(1)[:, None]) ** 2).sum(1)
        for fam in KFAMS:
            for a, v in enumerate(GRIDS[fam]):
                R = np.einsum("pc,tc->pt", prob, KERN[fam](G, v))
                emit(R * w[:, None], "%s_%d" % (fam, a))
                if fam in SCALED:
                    sc = base_mean / np.maximum(R.mean(1), 1e-300)
                    emit(R * sc[:, None] * w[:, None], "%sS_%d" % (fam, a))
            #: ST.achieved inputs -- what the operator DELIVERS, not what its grid label claims.
            #: The smeared distribution's MEAN must be unchanged (every window is symmetric about
            #: its own published centre) and its SD must inflate by exactly sqrt(1 + E[h^2/3]/var).
            if fam in ("DEQUANT", "DEQUANTI"):
                for a, v in enumerate(GRIDS[fam]):
                    h = v * HALF * (INTERIOR if fam == "DEQUANTI" else 1.0)
                    E0 = (prob * CEN[None]).sum(1)
                    row["dqmean_%s_%d" % (fam, a)] = float(E0.mean())
                    row["dqsmear_%s_%d" % (fam, a)] = float((prob * CEN[None]).sum(1).mean())
                    row["dqsd_%s_%d" % (fam, a)] = float(
                        (np.sqrt(var0 + (prob * (h[None] ** 2) / 3.0).sum(1))
                         / np.sqrt(np.maximum(var0, 1e-12))).mean())
        #: pair-level stratification of the intervention: modify ONLY the multimodal pairs
        mmask = np.zeros(len(prob), bool)
        p = prob
        pk = np.zeros(p.shape, bool)
        pk[:, 1:-1] = (p[:, 1:-1] > p[:, :-2]) & (p[:, 1:-1] > p[:, 2:]) & (p[:, 1:-1] > 0.02)
        pk[:, 0] = (p[:, 0] > p[:, 1]) & (p[:, 0] > 0.02)
        pk[:, -1] = (p[:, -1] > p[:, -2]) & (p[:, -1] > 0.02)
        mmask[:] = pk.sum(1) >= 2
        for fam, v in (("POWER", 2.0), ("TRUNC", 1.0), ("DEQUANTI", 1.0)):
            R = np.einsum("pc,tc->pt", prob, KERN[fam](G, v))
            R = np.where(mmask[:, None], R, base)
            emit(R * w[:, None], "MM_" + fam)

        for fam, S in (("METRIC", EGRID),):
            for a, v in enumerate(S):
                s = v * (Lmed - Lproj)
                A = G[None, :, None] + s[:, None, None] - CEN[None, None, :]
                row["int_%s_%d" % (fam, a)] = float(((1 - v) * Lmed + v * Lproj).mean())
                R = np.einsum("pc,ptc->pt", prob, np.abs(A))
                if v == 0.0 and not np.array_equal(R, base):
                    raise SystemExit("TRANSLATION PATH NOT IDENTICAL AT s=0 on %s" % pdb)
                emit(R * w[:, None], "%s_%d" % (fam, a))

        #: GATE 2 -- every family's identity arm must equal the INDEPENDENTLY computed f=1 column
        for fam in FAMS:
            for tag in ("%s_%d" % (fam, IDENT[fam]),
                        *(("%sS_%d" % (fam, IDENT[fam]),) if fam in SCALED else ())):
                if abs(row[tag] - ref[pdb]) > 1e-9:
                    raise SystemExit("IDENTITY ARM %s = %.6f but temper.json f=1 = %.6f on %s"
                                     % (tag, row[tag], ref[pdb], pdb))
        rows.append(row)
        if (ci + 1) % 10 == 0:
            print("  %d/%d  (identity %.4f)" % (ci + 1, len(tg),
                                                np.mean([r["POWER_2"] for r in rows])), flush=True)
            ST.save_atomic(OUT, {"rows": rows}, rows=rows, n_expected=len(tg), module_file=__file__)
        del u, Wp, D, X

    need = [k for k in rows[0] if k not in ("pdb", "n", "fold", "mmfrac")]
    ST.save_atomic(OUT, {"rows": rows, "grids": {k: list(v) for k, v in GRIDS.items()},
                         "ident": IDENT, "topm": TOPM},
                   complete_keys=need, rows=rows, n_expected=len(tg), module_file=__file__)
    report(rows)
    return rows


# ---------------------------------------------------------------- reporting
def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    pdbs = [r["pdb"] for r in rows]
    mmf = np.array([r["mmfrac"] for r in rows], float)
    gv = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731

    print("\n" + "=" * 104)
    print("s25/locrun.py  n = %d.  POINT CLOUD, shipped uniform top-75 coordinate average." % len(rows))
    print("Pool, readout and per-pair weight w are IDENTICAL in every arm; only the risk table")
    print("varies, so the only thing that changes is WHICH 75 candidates are selected.")
    print("=" * 104)
    base = gv("POWER_2")
    print("\n  identity arm through this module's own path: %.4f  (pinned incumbent 3.0483,"
          % base.mean())
    print("  and asserted equal per target to temper.json's independent f=1 column)\n")

    print("  %-12s %8s %9s %9s %9s %9s %9s" %
          ("arm", "param", "RMSD", "vs ident", "gam_eff", "cos", "amp"))
    for fam in FAMS:
        for suff in ("", "S"):
            if suff and fam not in SCALED:
                continue
            print("  %s" % ("-" * 70))
            for a, v in enumerate(GRIDS[fam]):
                tag = "%s%s_%d" % (fam, suff, a)
                x = gv(tag)
                print("  %-12s %8.2f %9.4f %+9.4f %+9.4f %+9.4f %9.3f"
                      % (fam + suff, v, x.mean(), (x - base).mean(),
                         gv("g_" + tag).mean(), gv("c_" + tag).mean(), gv("a_" + tag).mean()))
    print("  %s" % ("-" * 70))
    for tag in ("MM_POWER", "MM_TRUNC", "MM_DEQUANTI"):
        x = gv(tag)
        print("  %-12s %8s %9.4f %+9.4f %+9.4f %+9.4f %9.3f"
              % (tag, "mm-only", x.mean(), (x - base).mean(), gv("g_" + tag).mean(),
                 gv("c_" + tag).mean(), gv("a_" + tag).mean()))

    print("WHAT THE OPERATORS ACTUALLY DELIVERED (ST.achieved -- an arm's name is a claim and it")
    print("  is checkable.  s25's SD(f) delivered 1.43x at a nominal 1.5 while ALSO shifting the")
    print("  posterior mean by -0.375 A, and only the nominal was ever printed).")
    a1 = ST.achieved([gv("Emean").mean()], [gv("L_POWER_5").mean()],
                     "POWER q=2 must realise the posterior MEAN")
    print("    %-50s nominal %.4f achieved %.4f ratio %.4f%s"
          % (a1["label"], a1["nominal"][0], a1["achieved"][0], a1["ratio"][0],
             "  FLAG" if a1["flag"] else ""))
    for a, v in enumerate(GRIDS["METRIC"]):
        if v == 0.0:
            continue
        n_ = gv("int_METRIC_%d" % a).mean(); m_ = gv("L_METRIC_%d" % a).mean()
        r = ST.achieved([n_], [m_], "METRIC e=%.2f target = (1-e)Lmed + e Lproj" % v)
        print("    %-50s nominal %.4f achieved %.4f ratio %.4f%s"
              % (r["label"], n_, m_, r["ratio"][0], "  FLAG" if r["flag"] else ""))
    for fam in ("DEQUANTI", "DEQUANT"):
        for a, v in enumerate(GRIDS[fam]):
            if v == 0.0:
                continue
            r = ST.achieved([gv("dqmean_%s_%d" % (fam, a)).mean()],
                            [gv("dqsmear_%s_%d" % (fam, a)).mean()],
                            "%s b=%.2f must PRESERVE the posterior mean" % (fam, v))
            print("    %-50s ratio %.6f  sd inflation %.4f  minimiser OFF-centre %.3f%s"
                  % (r["label"], r["ratio"][0], gv("dqsd_%s_%d" % (fam, a)).mean(),
                     1.0 - gv("onc_%s_%d" % (fam, a)).mean(), "  FLAG" if r["flag"] else ""))
    print("")
    print("  CAVEAT CARRIED ON EVERY FIGURE FROM HERE: 4 of the 126 dev targets carry a VERBATIM")
    print("  self-sequence window at BLOSUM rank #1 (core/data.py's longer-normalised identity).")
    print("")
    print("\n  P1/P2, the predictions registered before this run.  CONVERSION = endpoint delta per")
    print("  unit gam_eff.  P1 says METRIC has the largest gam_eff and the worst conversion; P2")
    print("  says amp (how much of the move points along directions the POOL varies in) explains it.")
    print("  %-12s %10s %10s %12s %9s" % ("family", "gam_eff", "delta", "delta/gam", "amp"))
    for fam, a in (("METRIC", 4), ("DEQUANTI", 4), ("DEQUANT", 4), ("POWER", 5)):
        tag = "%s_%d" % (fam, a)
        g = gv("g_" + tag).mean(); d = (gv(tag) - base).mean()
        print("  %-12s %10.4f %+10.4f %12s %9.3f"
              % (tag, g, d, ("%+.3f" % (d / g)) if abs(g) > 1e-6 else "--", gv("a_" + tag).mean()))

    print("\n  NESTED CV over the 5 pinned folds -- one global parameter per family, chosen on the")
    print("  TRAINING folds by mean RMSD, applied to the held-out fold.  Paired against the family's")
    print("  own identity arm, which is the incumbent bit-for-bit.")
    for fam in FAMS:
        for suff in ("", "S"):
            if suff and fam not in SCALED:
                continue
            M = np.stack([gv("%s%s_%d" % (fam, suff, a)) for a in range(len(GRIDS[fam]))], 1)
            held = np.empty(len(rows)); pick = []
            for q in sorted(set(fold.tolist())):
                tr = fold != q; te = fold == q
                a = int(np.argmin(M[tr].mean(0)))
                held[te] = M[te][:, a]; pick.append(GRIDS[fam][a])
            r = ST.compare(held, base, fold, names=pdbs,
                           label="%s%s nested-CV vs identity" % (fam, suff))
            print("\n" + ST.fmt(r))
            print("    parameter chosen per fold: %s" % pick)
            #: the ceiling: the single best global parameter with FULL LEAKAGE on all 126
            a = int(np.argmin(M.mean(0)))
            print("    FULL-LEAKAGE best single parameter %s -> %.4f = %+.4f  (the family's ceiling)"
                  % (GRIDS[fam][a], M[:, a].mean(), M[:, a].mean() - base.mean()))
            #: the per-target oracle, scored with the SHARED library's within-grid null.  The
            #: pooled `best_of_k_null` is right for a candidate pool and WRONG for a K-column grid;
            #: temper.py's home-rolled own-row null was proven tautological (83% on pure noise).
            b = ST.best_of_k_within(M)
            print("    per-target ORACLE %+.4f;  VALID across-target null %+.4f -> %.0f%% accounted,"
                  " residual %+.4f, k_eff %.1f of %d"
                  % (M.min(1).mean() - base.mean(), b["null_across_targets"],
                     100 * b["share_accounted"], b["residual"], b["k_eff"], b["k"]))
            print("    SPLIT-HALF TRANSFER (self-nulling, the number to quote) %+.4f = %.0f%% of the"
                  " oracle -> %s" % (b["split_half"], 100 * b["split_half_frac"], b["verdict"]))
            print("    [the INVALID own-row null would have said %.0f%% -- printed so the trap stays"
                  " visible]" % (100 * b["share_accounted_INVALID"]))

    print("\n  STRATIFIED BY THE TARGET'S MULTIMODAL FRACTION (pre-registered by the coordinator")
    print("  BEFORE run 1).  A pair-level split is not expressible in a target-level endpoint, so")
    print("  the target's multimodal fraction is the stratifier and the MM_* arms above are the")
    print("  pair-level version -- the intervention applied to multimodal pairs ONLY.")
    hi = mmf >= np.median(mmf)
    print("  %-12s %8s %12s %12s" % ("arm", "param", "high-MM(%d)" % hi.sum(),
                                     "low-MM(%d)" % (~hi).sum()))
    for fam, a in (("METRIC", 4), ("DEQUANTI", 4), ("DEQUANT", 4), ("POWER", 5), ("TRUNC", 0)):
        tag = "%s_%d" % (fam, a)
        d = gv(tag) - base
        print("  %-12s %8.2f %+12.4f %+12.4f" % (fam, GRIDS[fam][a], d[hi].mean(), d[~hi].mean()))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        _lock()
        try:
            run()
        finally:
            _unlock()
