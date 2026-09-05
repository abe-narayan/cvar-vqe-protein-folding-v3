"""Does constraining the projection to plausible TORSION SPACE buy accuracy, validity, or both?

THE INCUMBENT AND ITS DEFECT
----------------------------
S8-11's winning arm (`s8/consensus2.py`, `sc|75:fit`) coordinate-averages the score-filtered
top-75 and then PROJECTS that average back onto the manifold of ideal-geometry chains:

    minimise  CA-RMSD( build(phi, psi), C )   over (phi, psi)

3.204 A against a 3.454 baseline, replicated on dev-24 at 3.221.  `s8/audit8.py`
independently reproduced it and then inspected the chains it returns.  The CA TRACE is
clean -- bonds exact to 1.7e-15, no pseudo-angle outside 75-150 deg on any of 126 targets,
fewer sub-4 A non-local contacts than the deposited natives.  The TORSIONS are not:

    12.1% positive phi against 4.1% in real library windows, and a pseudo-angle mean of
    93.5 deg against the natives' 104.0 -- the projection has no Ramachandran term and
    returns a somewhat over-helical conformation.

The objective above is indifferent to where in (phi, psi) it lands.  Real backbones occupy
a small, strongly non-uniform region of that space; the projection is free to leave it and
does.  This module adds the missing term:

    minimise  CA-RMSD( build(phi, psi), C )  +  lambda * penalty(phi, psi)

and sweeps `lambda` from 0 -- which must reproduce `fit` EXACTLY, and is asserted to --
upward, reporting the whole trade-off between torsion plausibility and selected CA-RMSD.

THE FOUR PENALTIES
------------------
`rama`   Per-residue -log P(phi, psi | residue class) under a wrapped-Gaussian-smoothed
         histogram estimated from TRAINING FOLDS ONLY.  Four classes -- GLY, PRO, PRE-PRO,
         GENERAL -- because glycine and proline have genuinely different allowed regions
         (glycine's positive-phi rate in this corpus is 42.1%, everything else's is 5.66%)
         and a single pooled table would either forbid glycine's real basins or license
         everyone else's excursions into them.  Bilinear interpolation of the LOG density
         on the wrapped grid, so the penalty is continuous and its finite-difference
         gradient is meaningful.
`rama20` The same construction conditioned on all 20 amino acids rather than 4 classes.
         The control for "is class conditioning enough".
`vm`     The SEQUENCE-CONDITIONED mixture of von Mises from `s5/torsion.py`, reused rather
         than reinvented: a dilated-convolution net emitting a per-residue K-component
         mixture over (phi, psi), trained leave-one-fold-out on the same corpus.  Strictly
         more informative than `rama` if the sequence signal is real, and the arm that
         tests whether it is.
`phip`   A smooth one-sided barrier on positive phi for non-glycine and nothing else.  The
         CRUDE control: it knows the symptom and not the distribution.  S7 measured that
         driving the positive-phi rate to exactly zero OVERSHOOTS real protein, so this arm
         is expected to fix the number and not the conformation, and it is here to be
         beaten rather than to win.

Plus one modifier, `w` (`--wdata`): weight the data term per residue by the local agreement
of the candidates being averaged, so the prior dominates exactly where the coordinate
average is least determined.  This needs a WEIGHTED Kabsch, which is implemented here and
tested against the unweighted one at uniform weights.

GEOMETRICALLY INERT PARAMETERS -- a correction to the audit's headline number
-----------------------------------------------------------------------------
`protein_geometry.build_backbone_batch` consumes `phi[1:]` and `psi[:-1]`.  `phi[0]` and
`psi[n-1]` are NEVER READ: they are free parameters of the optimisation that cannot move
the structure, so at the end of an unconstrained fit they hold whatever the optimiser left
there, which is noise.  The audit's 12.1% averages over all n residues including `phi[0]`,
so roughly one residue in n of that figure is a coin flip.  The same convention is what
puts `phi[0] = -60 deg` and `psi[n-1] = -45 deg` as literal constants in `peptide_db`.

This module therefore reports positive-phi on the CONSTRAINED set (`phi[1:]`) as the
honest number, alongside the all-residue number for like-for-like comparison with the
audit, and it excludes the inert parameters from the penalty.

DEPLOYABLE vs REPORTING
-----------------------
Every penalty is a function of (sequence, fold, torsions).  The Ramachandran tables are
accumulated from `peptide_db.load()` restricted to `folds[seq] != fold` plus
`distogram._fold_fragments(fold, 5)` -- the production discipline, and the identical corpus
`s8/generate.py` uses.  The von Mises net is `s5.torsion.train_fold(fold)`, same rule.
`rr` and `nat_ca` enter only as reporting labels.  `stage_leak` NaN-poisons both and
asserts every returned structure is bit-identical.

WHAT IT BOUGHT: VALIDITY, NOT ACCURACY
--------------------------------------
Stated plainly, because the mandate's bar was to beat 3.204 and this does not.

**No accuracy.**  The chosen arm is 3.202 against the incumbent's 3.204 at m=75, and the
filter-size table says that is a cell and not an effect: d(arm, fit) is +0.001, +0.002,
-0.002, -0.000 at m = 25, 50, 75, 150.  The concentration analysis closes it -- the
arm-vs-fit delta has a top-10 share of 1.034 and a mean excluding those ten of exactly
+0.000, which is the arithmetic statement that there is nothing there.  It is a uniform
null: -0.002 in both achievability strata and -0.005 to -0.000 across all five folds.
The tuning bound on the effect is |d| < 0.005 A.

**Validity, and for free.**  Positive phi on the constrained non-glycine set falls from
17.5% to 5.5%, landing on the 5.40% measured over the library windows the pools are drawn
from and inside the 2.90-5.66% band real backbones occupy.  The plateau holds -- 6.5%,
5.7%, 5.5%, 4.6% across the four filter sizes.  Nothing that already worked is disturbed:
bond geometry stays exact to 5.3e-15 A, clashes, rg and the pseudo-angle move in the fourth
decimal.

**The mechanism, and why the crude penalty won.**  A CA trace does not determine (phi, psi)
uniquely -- there are generally two ideal-geometry solutions at nearly the same distance
from the average, one plausible and one not -- and the unconstrained projection takes
whichever the data term reached.  So the defect costs the objective nothing, and removing
it costs the objective nothing.  The barrier beat both full density priors precisely
BECAUSE it is crude: its gradient is zero wherever phi is already negative, so it only
touches the residues that are wrong, while a mean-log-density penalty pulls plausible
residues toward their basin modes too and pays +0.008 to +0.054 A for the privilege.
Sparse support, not sophistication.

**What the prior does NOT fix, and this is the part to carry forward.**  Fixing the sign of
phi is not fixing the conformation.  Residues in no canonical basin fall from 22.8% to
14.3% against a real window's 1.6%; symmetric KL against the training corpus falls from
9.87 to 8.80 against a real window's 1.84.  And the residues it rescues become MORE
helical, not less -- alphaR 63.9% -> 70.1%, beta 7.4% -> 11.9% against a window's 22.4%.
The pseudo-angle stays at 93.2 deg against the natives' 104.0.  No torsion prior tested
here moves that: `vm`, the only sequence-conditioned one, moves it furthest (+3.0 deg) and
charges +0.054 A.  A marginal Ramachandran prior structurally CANNOT fix an over-helical
bias, because alpha is its mode -- it reinforces the defect it was brought in to correct.

**Dev-24, one pre-registered pass** (endpoints fixed before it ran, and it is the third
pass on this family, declared).  `fit` reproduces S8-11 exactly at 3.221, -0.255
[-0.501, -0.008], 16/8.  The arm is 3.220, d vs fit -0.001 [-0.008, +0.006], and positive
phi goes 17.66% -> 4.50% -- the defect and its correction both transfer to within a
percentage point.  Non-inferiority holds an order of magnitude inside its 0.05 A margin.
No accuracy gain is claimed and none was testable: dev-24's SE on a mean is 0.349 A.

**THE ARM TO CARRY FORWARD IS `ramah@0.3`, NOT THE ONE THE RULE PICKED.**  The
pre-declared rule ranked by RMSD subject to positive phi landing in the real band, and
positive phi is one bit of a two-dimensional distribution that both candidates satisfy
(5.45% and 5.62%).  On the distribution itself the hinge dominates: residues in NO
canonical basin fall to 2.90% against the barrier's 14.34%, a real library window's 1.59%
and the training corpus's 3.30%, and symmetric KL to 6.78 against 8.80.  The barrier only
sees the sign of phi; a residue at phi=-100, psi=-170 is invisible to it and forbidden by
the hinge.  The cost is +0.004 A [-0.005, +0.013], not distinguishable from zero.  This
choice is POST HOC on a criterion the original rule did not contain, and is labelled as
such; the dev pass stands on the arm that was pinned.

**The convergence hypothesis, tested and both true and beside the point.**  It was put to
this module that the over-helical bias is the projection faithfully absorbing the
compaction of its own target (mean CA-CA bond 2.961 A against a native 3.812), so a less
converged fit would be a more plausible one and no prior would be needed.  Measured within
targets over 13 iteration counts with no prior: rho(distance-to-average, positive-phi) is
-0.507 mean and negative on 98% of 109 targets, and the rate runs 1.1% -> 17.5% as the fit
tightens.  The mechanism is real.  The BETWEEN-target correlation is +0.559 -- the opposite
sign -- because hard targets both converge less and carry more positive phi, which is why
it had to be measured within.  But early stopping cannot exploit it: reaching the real rate
needs maxiter 3, which costs 3.204 -> 3.395, giving back three quarters of the S8-11 gain.
The prior reaches the same rate at full convergence for -0.002.  At matched plausibility
the prior is 0.19 A better than early stopping, and that is what the machinery is for.

    python -m s8.project prior     # build + inspect the fold-disciplined torsion priors
    python -m s8.project sweep     # the lambda x penalty trade-off curve at m=75
    python -m s8.project ms        # the filter-size plateau for the chosen arm
    python -m s8.project init      # torsion-seeded vs generic starts, WITH a prior
    python -m s8.project valid     # the full torsion + geometry validity table
    python -m s8.project conc      # per-target concentration, the audit's two fragilities
    python -m s8.project leak      # NaN-poison assertion
    python -m s8.project dev       # ONE pre-registered dev-24 pass (declared arm only)
    python -m s8.project report    # print everything already computed
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np                                                          # noqa: E402

os.environ.setdefault("NT", "2")
try:
    import torch
    torch.set_num_threads(2)
except Exception:                                                           # pragma: no cover
    torch = None

import protein_geometry as geo                                              # noqa: E402
from s7 import audit, debias, poolsize                                      # noqa: E402
from s8 import consensus2 as cc                                             # noqa: E402
from s8 import generate as gen                                              # noqa: E402
from s8 import inband                                                       # noqa: E402
from s8.inband import sel_of, _paired                                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

#: Ramachandran grid.  36 bins = 10 deg, matching `s8.generate.RB` so the two tables are
#: directly comparable, and coarse enough that a 787-chain + 6k-fragment corpus fills it.
RB = 36
#: Wrapped-Gaussian smoothing width, in BINS, applied to the counts before taking logs.
#: 1.5 bins = 15 deg.  Ramachandran basins are ~40-60 deg across, so this smooths sampling
#: noise without merging alpha and beta.
SIGMA_BINS = 1.5
#: Dirichlet pseudo-count per bin, added AFTER smoothing, so an unvisited corner has a
#: finite log-density and the penalty cannot be infinite.
PSEUDO = 0.5

#: Residue classes for the conditional table.  Order is the priority order.
CLASSES = ("gen", "gly", "pro", "prepro")
AA20 = "ARNDCQEGHILKMFPSTWYV"

#: The penalties this module can put on the projection.
PENALTIES = ("rama", "rama20", "vm", "phip", "ramah")

#: `ramah`'s hinge threshold, as a PERCENTILE of the log-density of real residues under
#: the same table.  5.0 means "penalise a residue only if it is less plausible than the
#: least plausible 5% of real residues" -- the mandate's calibration rule, applied to the
#: whole density rather than to the sign of phi.  It is measured per fold and per class
#: from training-fold residues, never chosen by looking at an RMSD.
HINGE_PCT = 5.0

#: The lambda ladder.  Geometric, because the interesting scale is unknown to within an
#: order of magnitude and a linear ladder would spend most of its cells in one regime.
#: The top of the ladder was moved DOWN from 3.0 to 1.0 and a 0.001 rung added after a
#: single-target pilot showed every penalty saturating by lam ~ 0.3 while the interesting
#: transition -- where the positive-phi rate crosses the 5.4% library reference -- happens
#: between 0.001 and 0.01.  The pilot ran before any aggregate was computed.
LAMS = (0.0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0)

#: Filter sizes for the plateau check -- the mandate's set.
MS = (25, 50, 75, 150)

#: The tuning cell everything is developed on.  `sc` = the shipped score's own order.
CELL_F = "sc"
CELL_M = 75

#: `phip` barrier softness, in units of sin(phi).  0.15 puts the barrier at 0.97 by
#: phi = +30 deg and 0.03 by phi = -30 deg: sharp enough to be a barrier, smooth enough
#: that a finite-difference gradient at eps = 1e-5 rad is meaningful.
PHIP_TAU = 0.15


def _p(name):
    return os.path.join(HERE, name)


def _write(name, obj):
    tmp = _p(name) + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, _p(name))


def _read(name, default=None):
    try:
        with open(_p(name)) as fh:
            return json.load(fh)
    except Exception:
        return default


def _wrap(a):
    return (np.asarray(a, float) + math.pi) % (2 * math.pi) - math.pi


# ============================================================ the torsion priors
def res_classes(seq):
    """``(n,)`` int class index per residue.  GLY / PRO / PRE-PRO / GENERAL.

    Priority matters and is stated rather than implied: a glycine FOLLOWED by proline is
    scored as glycine, because glycine's own backbone freedom dominates the pre-proline
    restriction.  Proline itself is scored as proline for the same reason.
    """
    n = len(seq)
    out = np.zeros(n, int)
    for i, c in enumerate(seq):
        if c == "G":
            out[i] = 1
        elif c == "P":
            out[i] = 2
        elif i + 1 < n and seq[i + 1] == "P":
            out[i] = 3
    return out


def _bin(a):
    """Wrapped bin index on the RB grid for angles in radians."""
    return np.mod(((np.asarray(a, float) + math.pi) / (2 * math.pi) * RB).astype(int), RB)


def _accum(counts, seq, phi, psi, aa_counts=None):
    """Accumulate ONE chain into the class table (and optionally the 20-aa table).

    `phi[0]` and `psi[-1]` are the geometrically inert parameters `peptide_db` fills with
    the constants -60 / -45 deg.  Including them would put a spike of pure convention at
    those two coordinates, so the interior mask drops them: only residues 1..n-2 carry a
    real (phi, psi) PAIR.
    """
    n = len(seq)
    if n < 3:
        return
    cls = res_classes(seq)
    bp, bs = _bin(phi), _bin(psi)
    sl = slice(1, n - 1)
    np.add.at(counts, (cls[sl], bp[sl], bs[sl]), 1.0)
    if aa_counts is not None:
        a = np.array([AA20.find(c) for c in seq], int)
        ok = a[sl] >= 0
        np.add.at(aa_counts, (a[sl][ok], bp[sl][ok], bs[sl][ok]), 1.0)


def _corpus(fold, n_folds=5):
    """Every chain a model for `fold` is allowed to see.  The production discipline."""
    import distogram as dgm
    import peptide_db as db
    folds = db.folds(n_folds)
    peps = [q for q in db.load() if folds[q.seq] != fold]
    return peps + list(dgm._fold_fragments(fold, n_folds))


def build_tables(verbose=True):
    """Per-fold (4, RB, RB) class counts and (20, RB, RB) amino-acid counts. Cached."""
    got = _read("project_prior.json")
    if got is not None:
        return got
    cls = np.zeros((5, len(CLASSES), RB, RB))
    aa = np.zeros((5, 20, RB, RB))
    ref = []
    for f in range(5):
        t0 = time.time()
        ent = _corpus(f)
        for e in ent:
            _accum(cls[f], e.seq, np.asarray(e.phi, float), np.asarray(e.psi, float), aa[f])
        # reference statistics on the SAME corpus, interior residues only
        pp = np.concatenate([_wrap(np.asarray(e.phi, float))[1:len(e.seq) - 1]
                             for e in ent if len(e.seq) >= 3])
        sq = "".join(e.seq[1:len(e.seq) - 1] for e in ent if len(e.seq) >= 3)
        g = np.array([c == "G" for c in sq])
        ref.append({"fold": f, "chains": len(ent), "residues": int(cls[f].sum()),
                    "posphi_all": float((pp > 0).mean()),
                    "posphi_nongly": float((pp[~g] > 0).mean()),
                    "posphi_gly": float((pp[g] > 0).mean()),
                    "gly_frac": float(g.mean())})
        if verbose:
            print(f"  fold {f}: {len(ent)} chains, {int(cls[f].sum())} interior residues, "
                  f"pos-phi non-Gly {ref[-1]['posphi_nongly']:.4f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    out = {"RB": RB, "classes": list(CLASSES),
           "cls": cls.astype(int).tolist(), "aa": aa.astype(int).tolist(),
           "corpus_ref": ref}
    _write("project_prior.json", out)
    return out


def _wrap_smooth(C, sigma_bins=SIGMA_BINS):
    """Wrapped-Gaussian smoothing of a (..., RB, RB) count table, by FFT convolution."""
    k = np.arange(RB)
    d = np.minimum(k, RB - k)
    g = np.exp(-0.5 * (d / sigma_bins) ** 2)
    g = g / g.sum()
    G = np.fft.rfft(g)
    A = np.fft.irfft(np.fft.rfft(C, axis=-1) * G, n=RB, axis=-1)
    A = np.fft.irfft(np.fft.rfft(A, axis=-2) * G[:, None], n=RB, axis=-2)
    return np.maximum(A, 0.0)


_TAB = {}


def logp_tables(kind="rama"):
    """``(5, C, RB, RB)`` log-densities: smoothed, pseudo-counted, normalised per class."""
    if kind not in _TAB:
        raw = build_tables(verbose=False)
        C = np.asarray(raw["cls" if kind == "rama" else "aa"], float)
        S = _wrap_smooth(C) + PSEUDO
        _TAB[kind] = np.log(S / S.sum((2, 3), keepdims=True))
    return _TAB[kind]


def _bilinear(L, phi, psi):
    """Bilinear interpolation of a (K, RB, RB) LOG-density at per-residue (phi, psi).

    `phi`/`psi` are ``(B, n)``; `L` is indexed by a ``(n,)`` per-residue key supplied by
    the caller through `L` already gathered to ``(n, RB, RB)``.  Wrapping is exact on both
    axes.  Interpolating the LOG rather than the density keeps the penalty finite
    everywhere and makes its gradient the natural one.
    """
    u = (np.asarray(phi, float) + math.pi) / (2 * math.pi) * RB - 0.5
    v = (np.asarray(psi, float) + math.pi) / (2 * math.pi) * RB - 0.5
    i0 = np.floor(u).astype(int)
    j0 = np.floor(v).astype(int)
    fu, fv = u - i0, v - j0
    i0 %= RB
    j0 %= RB
    i1, j1 = (i0 + 1) % RB, (j0 + 1) % RB
    r = np.arange(L.shape[0])[None, :]
    return ((1 - fu) * (1 - fv) * L[r, i0, j0] + fu * (1 - fv) * L[r, i1, j0]
            + (1 - fu) * fv * L[r, i0, j1] + fu * fv * L[r, i1, j1])


class RamaPenalty:
    """-mean log P(phi, psi) over the CONSTRAINED residues, from training folds only."""

    def __init__(self, seq, fold, kind="rama"):
        self.seq, self.fold, self.kind = seq, fold, kind
        n = len(seq)
        T = logp_tables(kind)[fold]
        key = (res_classes(seq) if kind == "rama"
               else np.array([max(AA20.find(c), 0) for c in seq], int))
        self.L = T[key]                                     # (n, RB, RB)
        #: residues whose (phi, psi) PAIR both drive geometry: phi[i] for i>=1, psi[i]
        #: for i<=n-2.  Index 0 and n-1 are half-inert, so the pair density is not the
        #: right object there and they are excluded from the penalty entirely.
        self.mask = np.zeros(n, bool)
        self.mask[1:n - 1] = True
        if not self.mask.any():                             # n < 3, degenerate
            self.mask[:] = True

    def __call__(self, phi, psi):
        lp = _bilinear(self.L, np.atleast_2d(phi), np.atleast_2d(psi))
        return -lp[:, self.mask].mean(1)


class VMPenalty:
    """Mean NLL under `s5.torsion`'s sequence-conditioned von Mises mixture.

    `s5.torsion.for_target` selects the leave-fold-out model by the target's own identity
    fold, which is the same discipline the tables use.  The mixture is evaluated with
    `torch.no_grad`; nothing here trains.
    """

    def __init__(self, seq, fold):
        from s5 import torsion as t5
        self.t5 = t5
        self.prior = t5.for_target(seq)
        n = len(seq)
        self.mask = np.zeros(n, bool)
        self.mask[1:n - 1] = True
        if not self.mask.any():
            self.mask[:] = True

    def __call__(self, phi, psi):
        P = torch.as_tensor(np.atleast_2d(phi), dtype=torch.float32)
        S = torch.as_tensor(np.atleast_2d(psi), dtype=torch.float32)
        B = P.shape[0]
        p = {k: v[None].expand(B, *v.shape) for k, v in self.prior.p.items()}
        with torch.no_grad():
            ll = self.t5.mixture_logpdf(p, P, S)
        return (-ll.numpy()[:, self.mask]).mean(1)


class PhiPosPenalty:
    """Smooth one-sided barrier on positive phi, non-glycine only.  The CRUDE control.

    ``sigmoid(sin(phi) / tau)`` rather than ``sigmoid(phi / tau)``, because phi lives on a
    circle.  A bare logistic in phi is discontinuous at the wrap: it would read +179 deg as
    maximally forbidden and -179 deg -- the same conformation to within 2 deg -- as
    maximally allowed, and the optimiser would escape the barrier by walking phi past
    +180.  The sine form is exactly 2*pi-periodic, is 0.5 at BOTH transitions (0 and
    +-180 deg, where the sign of phi is genuinely undefined), and satisfies
    ``s(phi) + s(-phi) = 1``.
    """

    def __init__(self, seq, fold, tau=PHIP_TAU):
        self.tau = tau
        n = len(seq)
        m = np.zeros(n, bool)
        m[1:n - 1] = True
        if not m.any():
            m[:] = True
        self.mask = m & np.array([c != "G" for c in seq])
        if not self.mask.any():                             # all-glycine: nothing to hold
            self.mask = m

    def __call__(self, phi, psi):
        s = 1.0 / (1.0 + np.exp(-np.sin(np.atleast_2d(phi)) / self.tau))
        return s[:, self.mask].mean(1)


_HINGE = {}


def hinge_thresholds(pct=HINGE_PCT):
    """``(5, C)`` per-fold, per-class log-density below which a residue counts as
    implausible: the `pct`-th percentile of the log-density of REAL residues of that class
    under that fold's own table.

    Computed by weighting each grid cell by how many real residues fall in it, so the
    threshold is a percentile over residues and not over grid area -- the two differ by a
    lot, because most of Ramachandran space is empty and would otherwise set the bar.
    """
    key = round(pct, 6)
    if key not in _HINGE:
        raw = build_tables(verbose=False)
        cnt = np.asarray(raw["cls"], float)                     # (5, C, RB, RB) counts
        L = logp_tables("rama")                                 # (5, C, RB, RB) log p
        out = np.zeros(L.shape[:2])
        for f in range(L.shape[0]):
            for k in range(L.shape[1]):
                v = L[f, k].ravel()
                w = cnt[f, k].ravel()
                o = np.argsort(v)
                cw = np.cumsum(w[o])
                tot = cw[-1]
                if tot <= 0:                                    # no data: never fires
                    out[f, k] = v.min()
                    continue
                out[f, k] = v[o[int(np.searchsorted(cw, tot * pct / 100.0))]]
        _HINGE[key] = out
    return _HINGE[key]


class RamaHingePenalty(RamaPenalty):
    """The lesson of the survey sweep, applied to the whole density instead of to phi.

    The one-sided phi barrier beat both full density priors on the statistic all three
    target, and it was the only one that cost nothing.  The reason is SUPPORT: the
    barrier's gradient is exactly zero wherever phi is already negative, so it acts only on
    the residues that are wrong, while a mean-log-density penalty has gradient everywhere
    and drags plausible residues toward their basin modes as well.  That drag is what buys
    the +0.008 to +0.054 A the density priors cost.

    So: hinge the log-density at a threshold calibrated to real residues.

        pen = mean_i  max(0, t_class(i) - log P(phi_i, psi_i | class(i)))

    A residue already as plausible as 95% of real residues contributes nothing and feels no
    force.  A residue in a region real protein does not occupy is pulled until it is merely
    unusual, and then released.  This is the mandate's "calibrate to the measured rate
    rather than to zero" applied to the density rather than to the sign of phi, and unlike
    the phi barrier it can also see psi.
    """

    def __init__(self, seq, fold, kind="rama", pct=HINGE_PCT):
        super().__init__(seq, fold, kind=kind)
        t = hinge_thresholds(pct)[fold]
        self.t = t[res_classes(seq)]                            # (n,)

    def __call__(self, phi, psi):
        lp = _bilinear(self.L, np.atleast_2d(phi), np.atleast_2d(psi))
        return np.maximum(self.t[None, :] - lp, 0.0)[:, self.mask].mean(1)


def make_penalty(kind, seq, fold):
    if kind == "ramah":
        return RamaHingePenalty(seq, fold)
    if kind in ("rama", "rama20"):
        return RamaPenalty(seq, fold, kind=kind)
    if kind == "vm":
        return VMPenalty(seq, fold)
    if kind == "phip":
        return PhiPosPenalty(seq, fold)
    raise KeyError(kind)


# ============================================================ weighted superposition
def wkabsch_rmsd_batch(P, ref, w):
    """Weighted CA-RMSD of every ``(n,3)`` in `P` against `ref`.  Reflections forbidden.

    The same construction as `s7.audit.kabsch_rmsd_batch` with a per-atom weight: weighted
    centroids, weighted cross-covariance, and the residual read off the singular values
    with the smallest sign-flipped when det(V U^T) < 0.  At uniform weights it is
    numerically identical to the unweighted routine, which `test_project.py` asserts.
    """
    P = np.asarray(P, float)
    ref = np.asarray(ref, float)
    if P.ndim == 2:
        P = P[None]
    w = np.asarray(w, float)
    w = w / w.sum()
    Pc = P - (w[None, :, None] * P).sum(1, keepdims=True)
    Rc = ref - (w[:, None] * ref).sum(0, keepdims=True)
    Wp = Pc * w[None, :, None]
    H = np.einsum("bni,nj->bij", Wp, Rc)
    U, S, Vt = np.linalg.svd(H)
    det = np.linalg.det(np.einsum("bji,bkj->bik", Vt, U))
    S = S.copy()
    S[:, -1] *= np.sign(det)
    resid = (Wp * Pc).sum((1, 2)) + (w[:, None] * Rc * Rc).sum() - 2.0 * S.sum(1)
    return np.sqrt(np.maximum(resid, 0.0))


def agreement_weights(Wsub, C, floor=0.25):
    """Per-residue data weights from the LOCAL AGREEMENT of the averaged candidates.

    After superposing the subset onto its medoid, residue `i`'s scatter `s_i` is the RMS
    displacement of the candidates from the average at that position.  Weight is
    ``1 / (s_i + s0)`` with ``s0`` the median scatter, so the profile is bounded (no single
    tightly-agreed residue can take the whole weight) and the prior is left to dominate
    exactly where the average is least determined.  `floor` keeps the flattest possible
    profile from collapsing to a delta.
    """
    d = np.sqrt(((Wsub - C[None]) ** 2).sum(-1))            # (B, n)
    s = np.sqrt((d ** 2).mean(0))
    s0 = max(float(np.median(s)), 1e-6)
    w = 1.0 / (s + s0)
    w = w / w.mean()
    return np.maximum(w, floor)


# ============================================================ the constrained projection
#: The four generic starting conformations `s8.consensus2.fit_consensus` uses.  Imported
#: by value rather than reimplemented so the lambda = 0 reproduction is exact.
STARTS = cc.FIT_STARTS


def fit_prior(C, phi0, psi0, pen=None, lam=0.0, w=None, maxiter=300):
    """The ideal-geometry chain nearest `C` that is ALSO plausible in torsion space.

        minimise   RMSD( build(phi, psi), C )  +  lam * pen(phi, psi)

    At ``lam = 0`` and ``w = None`` this is `s8.consensus2.fit_torsions` line for line --
    same objective, same finite-difference gradient, same L-BFGS-B call -- and
    `test_project.py` asserts the two agree to machine precision.

    The gradient is the whole cost, so all 2n perturbations are built in ONE batched call
    and the penalty is evaluated on the SAME batch.  Handed to scipy as `jac` this is ~27x
    fewer builds per iteration than letting L-BFGS-B difference the scalar objective.
    """
    from scipy.optimize import minimize
    C = np.asarray(C, float)
    n = len(C)
    eps = 1e-5
    E = np.eye(2 * n) * eps

    def build(X):
        return geo.build_backbone_batch(X[:, :n], X[:, n:])["CA"]

    def fg(x):
        X = np.vstack([x[None], x[None] + E])
        B = build(X)
        r = (audit.kabsch_rmsd_batch(B, C) if w is None
             else wkabsch_rmsd_batch(B, C, w))
        if lam and pen is not None:
            r = r + lam * pen(X[:, :n], X[:, n:])
        return float(r[0]), (r[1:] - r[0]) / eps

    x0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
    f0 = fg(x0)[0]
    r = minimize(fg, x0, jac=True, method="L-BFGS-B", options={"maxiter": maxiter})
    ph, ps = r.x[:n].copy(), r.x[n:].copy()
    F = build(r.x[None])[0]
    d = float(audit.kabsch_rmsd_batch(F[None], C)[0])
    return F, ph, ps, float(r.fun), f0, d


def fit_multi(C, pen=None, lam=0.0, w=None, extra=None, maxiter=300):
    """`fit_prior` from every generic start plus an optional extra, keeping the best."""
    n = len(C)
    starts = [(np.full(n, math.radians(a)), np.full(n, math.radians(b)))
              for a, b in STARTS]
    if extra is not None:
        starts.append((np.asarray(extra[0], float), np.asarray(extra[1], float)))
    best = None
    for ph, ps in starts:
        got = fit_prior(C, ph, ps, pen=pen, lam=lam, w=w, maxiter=maxiter)
        if best is None or got[3] < best[3]:
            best = got
    return best


def lam_path(C, pen, lams, w=None, extra=None, maxiter=300, multi=False):
    """The path in lambda.  Solve at lam=0 from the generic starts, then walk up the ladder.

    ``multi=False`` is pure CONTINUATION: each solve is warm-started from the one below it
    and nothing else.  It is the cheap form and it is what the first full sweep used.

    ``multi=True`` adds the four generic starts back at EVERY rung and keeps whichever
    lands the lowest objective.  This is not a refinement of taste -- the continuation path
    demonstrably under-optimises.  The CA trace does not determine (phi, psi) uniquely: for
    a given pseudo-angle and CA dihedral there are generally two ideal-geometry torsion
    solutions, one Ramachandran-plausible and one not, and they differ in RMSD-to-C by
    almost nothing.  A path started at the lam=0 optimum sits in whichever of those the
    data term happened to reach and L-BFGS cannot cross the barrier between them: on 1A1P,
    1CS9 and 1I6Y the multi-start finds a STRICTLY LOWER objective than the path at the
    same lambda, with a positive-phi rate half as large.  The prior's real job is choosing
    among near-degenerate solutions, and it can only do that if the optimiser can see them.

    Returns ``{lam: (F, phi, psi, fval, f0, d_to_C)}`` with the lam=0 entry EXACTLY the
    incumbent `fit`.
    """
    out = {}
    cur = fit_multi(C, pen=None, lam=0.0, w=w, extra=extra, maxiter=maxiter)
    for lam in lams:
        if lam == 0.0:
            out[lam] = cur if w is None else fit_multi(C, pen=None, lam=0.0, w=w,
                                                       extra=extra, maxiter=maxiter)
            continue
        got = fit_prior(C, cur[1], cur[2], pen=pen, lam=lam, w=w, maxiter=maxiter)
        if multi:
            alt = fit_multi(C, pen=pen, lam=lam, w=w, extra=extra, maxiter=maxiter)
            if alt[3] < got[3]:
                got = alt
        out[lam] = got
        cur = got
    return out


# ============================================================ instrument
def targets():
    return cc.cached()


def cell_of(c, f=CELL_F, m=CELL_M):
    """The filtered subset, its medoid index, and the coordinate average.  DEPLOYABLE."""
    sub = cc.filter_orders(c)[f][:m]
    Ps = c["P"][np.ix_(sub, sub)]
    b = int(np.argmin(cc._medoid(Ps)))
    C = cc._avg_struct(c["W"][sub], Ps)
    return sub, b, C


# ============================================================ torsion + geometry stats
def torsion_report(seq, phi, psi):
    """Positive-phi rates and the (phi, psi) occupancy of ONE torsion set.

    `all` is the audit's convention -- every residue, including the geometrically inert
    `phi[0]`.  `con` is the constrained set `phi[1:]`, which is the honest number.
    """
    n = len(seq)
    a = _wrap(phi)
    g = np.array([ch == "G" for ch in seq])
    con = np.zeros(n, bool)
    con[1:] = True
    ng = con & ~g
    out = {"posphi_all": float((a > 0).mean()),
           "posphi_all_nongly": float((a[~g] > 0).mean()) if (~g).any() else float("nan"),
           "posphi_con": float((a[con] > 0).mean()) if con.any() else float("nan"),
           "posphi_con_nongly": float((a[ng] > 0).mean()) if ng.any() else float("nan"),
           "n_con_nongly": int(ng.sum()), "n_pos_con_nongly": int((a[ng] > 0).sum())}
    return out


def geometry_report(P):
    """Bonds, CA-CA-CA pseudo-angles, non-local clashes, CA-trace chirality, rg.

    Identical in definition to `s8.audit8.geometry_report`, so the numbers here are
    directly comparable to the audit's table.
    """
    P = np.asarray(P, float)
    n = len(P)
    step = np.linalg.norm(P[1:] - P[:-1], axis=1)
    u = P[1:] - P[:-1]
    u = u / np.linalg.norm(u, axis=1, keepdims=True)
    ang = np.degrees(np.arccos(np.clip((u[:-1] * u[1:]).sum(1), -1, 1)))
    pang = 180.0 - ang
    if n >= 4:
        b0, b1, b2 = P[1:-2] - P[:-3], P[2:-1] - P[1:-2], P[3:] - P[2:-1]
        n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
        mm = np.cross(n1, b1 / np.linalg.norm(b1, axis=1, keepdims=True))
        dih = np.degrees(np.arctan2((mm * n2).sum(1), (n1 * n2).sum(1)))
    else:
        dih = np.array([np.nan])
    i, j = np.nonzero((np.arange(n)[None, :] - np.arange(n)[:, None]) >= 3)
    nl = np.linalg.norm(P[i] - P[j], axis=1) if len(i) else np.array([np.inf])
    d2 = ((P[:, None, :] - P[None, :, :]) ** 2).sum(-1)
    return {"step_mean": float(step.mean()), "step_sd": float(step.std()),
            "pseudoangle_mean": float(pang.mean()),
            "frac_pseudoangle_out_of_75_150": float(((pang < 75) | (pang > 150)).mean()),
            "min_nonlocal_CA": float(nl.min()),
            "frac_nonlocal_under_4A": float((nl < 4.0).mean()),
            "frac_dihedral_positive": float(np.mean(dih > 0)),
            "rg": float(np.sqrt(d2.sum() / (2.0 * n * n)))}


def _agg(rows):
    ks = rows[0].keys()
    return {k: float(np.mean([r[k] for r in rows])) for k in ks}


# ============================================================ stage: prior
def stage_prior(verbose=True):
    """Build the fold-disciplined tables and state what the corpus says the target is."""
    raw = build_tables(verbose=verbose)
    L = logp_tables("rama")
    out = {"corpus_ref": raw["corpus_ref"], "classes": list(CLASSES), "per_class": {}}
    cls = np.asarray(raw["cls"], float).sum(0)                      # pooled over folds
    ii, jj = np.meshgrid(np.arange(RB), np.arange(RB), indexing="ij")
    phic = (ii + 0.5) / RB * 2 * math.pi - math.pi
    for k, nm in enumerate(CLASSES):
        c = cls[k]
        out["per_class"][nm] = {"residues": float(c.sum()),
                                "posphi": float(c[phic > 0].sum() / max(c.sum(), 1))}
    # what the smoothed table itself implies, per fold, as a self-consistency check
    P = np.exp(L)
    out["table_posphi"] = {nm: float(P[:, k][:, phic > 0].sum() / 5.0)
                           for k, nm in enumerate(CLASSES)}
    # the reference rates the emitted structures are judged against
    nat, win = [], []
    for t in targets():
        c = cc.load_cache(t)
        # library windows carry real torsions; natives do not (CA only) -- so the
        # native reference is the pseudo-angle, and the window reference is pos-phi.
        seq = c["seq"]
        for b in range(min(50, len(c["PHI"]))):
            win.append(torsion_report(seq, c["PHI"][b], c["PSI"][b])["posphi_con_nongly"])
        nat.append(geometry_report(c["nat_ca"])["pseudoangle_mean"])
    out["ref_window_posphi_con_nongly"] = float(np.nanmean(win))
    out["ref_native_pseudoangle"] = float(np.mean(nat))
    _write("project_prior_report.json", out)
    _print_prior(out)
    return out


def _print_prior(out):
    print("\n=== FOLD-DISCIPLINED TORSION PRIOR ===")
    for r in out["corpus_ref"]:
        print(f"  fold {r['fold']}: {r['chains']:5d} chains  {r['residues']:7d} interior "
              f"residues   pos-phi all {r['posphi_all']:.4f}  "
              f"non-Gly {r['posphi_nongly']:.4f}  Gly {r['posphi_gly']:.4f}")
    print("  class occupancy (pooled counts / smoothed table):")
    for nm in out["classes"]:
        print(f"    {nm:8} {out['per_class'][nm]['residues']:9.0f} residues   "
              f"pos-phi {out['per_class'][nm]['posphi']:.4f} / "
              f"{out['table_posphi'][nm]:.4f}")
    print(f"  REFERENCE  library-window pos-phi (constrained, non-Gly) "
          f"{out['ref_window_posphi_con_nongly']:.4f}")
    print(f"  REFERENCE  native CA-CA-CA pseudo-angle mean "
          f"{out['ref_native_pseudoangle']:.2f} deg")


# ============================================================ stage: sweep
def _row_stats(F, seq, nat):
    r = {"rmsd": float(audit.kabsch_rmsd_batch(F[None], nat)[0])}
    return r


#: The narrower ladder the multi-start sweep runs, chosen from the continuation sweep.
LAMS_X = (0.0, 0.01, 0.03, 0.1, 0.3)
#: The ladder for the SPARSE-SUPPORT penalties.  `phip` and `ramah` are zero on the
#: plausible region, so their effective strength per ACTIVE residue is much larger than
#: their lambda suggests and the interesting range sits higher.
LAMS_H = (0.0, 0.03, 0.1, 0.3, 1.0)


def stage_sweep(pens=PENALTIES, lams=LAMS, m=CELL_M, resume=True, verbose=True,
                multi=False, tag=""):
    """THE TRADE-OFF CURVE.  lambda x penalty at the tuning cell, 126 targets, resumable."""
    name = f"project_sweep{tag}_m{m}.json"
    out = _read(name, {"m": m, "lams": list(lams), "pens": list(pens), "multi": multi,
                       "per": {}}) \
        if resume else {"m": m, "lams": list(lams), "pens": list(pens), "multi": multi,
                        "per": {}}
    per = out["per"]
    tgs = targets()
    todo = [t for t in tgs if t not in per]
    if verbose:
        print(f"{len(todo)}/{len(tgs)} targets to sweep (m={m})", flush=True)
    for k, pdbid in enumerate(todo):
        t0 = time.time()
        c = cc.load_cache(pdbid)
        seq, fold, nat = c["seq"], int(c["fold"]), c["nat_ca"]
        sub, b, C = cell_of(c, CELL_F, m)
        gi = int(sub[b])
        row = {"fold": fold, "n": int(c["n"]), "seq": seq,
               "base": sel_of(c["sc"], c["rr"]),
               "medoid": float(c["rr"][gi]),
               "pool_best": float(c["rr"].min()),
               "pool_best_sub": float(c["rr"][sub].min()),
               "avg": float(audit.kabsch_rmsd_batch(C[None], nat)[0]),
               "arms": {}}
        for pk in pens:
            try:
                pen = make_penalty(pk, seq, fold)
            except Exception as exc:                        # pragma: no cover
                row["arms"][pk] = {"error": repr(exc)}
                continue
            path = lam_path(C, pen, lams, multi=multi)
            row["arms"][pk] = {
                f"{lam:g}": {"rmsd": float(audit.kabsch_rmsd_batch(v[0][None], nat)[0]),
                             "d_to_C": v[5], "obj": v[3],
                             "pen": float(pen(v[1], v[2])[0]),
                             "tor": torsion_report(seq, v[1], v[2]),
                             "geo": geometry_report(v[0])}
                for lam, v in path.items()}
        per[pdbid] = row
        _write(name, out)
        if verbose:
            f0 = row["arms"][pens[0]]["0"]["rmsd"]
            print(f"[{k+1:3d}/{len(todo)}] {pdbid} n={row['n']} base={row['base']:.3f} "
                  f"fit={f0:.3f} ({time.time()-t0:.1f}s)", flush=True)
    if len(per) == len(tgs):
        out["summary"] = _summarise_sweep(out)
        _write(name, out)
        _print_sweep(out)
    return out


def _summarise_sweep(out):
    per = out["per"]
    tgs = sorted(per)
    base = np.array([per[t]["base"] for t in tgs], float)
    med = np.array([per[t]["medoid"] for t in tgs], float)
    fit = np.array([per[t]["arms"][out["pens"][0]]["0"]["rmsd"] for t in tgs], float)
    S = {"n": len(tgs), "base": float(base.mean()), "medoid": float(med.mean()),
         "fit": float(fit.mean()), "cells": {}}
    for pk in out["pens"]:
        for lam in out["lams"]:
            key = f"{lam:g}"
            if key not in per[tgs[0]]["arms"].get(pk, {}):
                continue
            v = np.array([per[t]["arms"][pk][key]["rmsd"] for t in tgs], float)
            tp = np.array([per[t]["arms"][pk][key]["tor"]["posphi_con_nongly"]
                           for t in tgs], float)
            ta = np.array([per[t]["arms"][pk][key]["tor"]["posphi_all"] for t in tgs], float)
            pa = np.array([per[t]["arms"][pk][key]["geo"]["pseudoangle_mean"]
                           for t in tgs], float)
            dc = np.array([per[t]["arms"][pk][key]["d_to_C"] for t in tgs], float)
            S["cells"][f"{pk}|{key}"] = {
                "sel": float(v.mean()),
                "se": float(v.std(ddof=1) / math.sqrt(len(v))),
                "median": float(np.median(v)),
                "f2": float((v < 2.0).mean()), "f15": float((v < 1.5).mean()),
                "worst": float(v.max()), "best": float(v.min()),
                "vs_base": _paired(v, base), "vs_fit": _paired(v, fit),
                "vs_medoid": _paired(v, med),
                "posphi_con_nongly": float(np.nanmean(tp)),
                "posphi_all": float(ta.mean()),
                "pseudoangle": float(pa.mean()), "d_to_C": float(dc.mean())}
    return S


def _print_sweep(out):
    S = out["summary"]
    print(f"\n=== TORSION-CONSTRAINED PROJECTION — trade-off curve, m={out['m']}, "
          f"{S['n']} targets, "
          f"{'MULTI-START' if out.get('multi') else 'continuation'} ===")
    print(f"  baseline (shipped score) {S['base']:.3f}   medoid {S['medoid']:.3f}   "
          f"fit (lam=0) {S['fit']:.3f}")
    print(f"{'arm':14} {'sel':>7} {'d vs 3.454':>11} {'95% CI':>18} {'W/L':>9} "
          f"{'d vs fit':>9} {'95% CI':>18} {'W/L':>9} {'phi+ nG':>8} {'phi+ all':>9} "
          f"{'pang':>7} {'f<2':>6}")
    for key, d in S["cells"].items():
        a, b_ = d["vs_base"], d["vs_fit"]
        print(f"{key:14} {d['sel']:7.3f} {a['mean_diff']:+11.3f} "
              f"[{a['ci95'][0]:+.3f},{a['ci95'][1]:+.3f}]".rjust(19) +
              f" {a['n_better']:3d}/{a['n_worse']:<4d} {b_['mean_diff']:+9.3f} " +
              f"[{b_['ci95'][0]:+.3f},{b_['ci95'][1]:+.3f}]".rjust(19) +
              f" {b_['n_better']:3d}/{b_['n_worse']:<4d} "
              f"{d['posphi_con_nongly']:8.4f} {d['posphi_all']:9.4f} "
              f"{d['pseudoangle']:7.2f} {d['f2']:6.3f}")


# ============================================================ stage: conv
#: Iteration counts the convergence trace is sampled at.  Geometric, and it reaches 300 --
#: `fit_torsions`'s own `maxiter` -- so the last snapshot IS the incumbent.
SNAPS = (1, 2, 3, 5, 8, 12, 20, 30, 50, 80, 120, 200, 300)


def fit_trace(C, phi0, psi0, snaps=SNAPS, maxiter=300):
    """`fit_prior` at lam=0, returning the iterate at each of `snaps` iterations.

    One solve gives the whole convergence curve, so a maxiter ladder costs no more than a
    single fit.  Snapshots past termination repeat the final iterate, which is the correct
    reading of "run L-BFGS-B with this maxiter" when it converged sooner.
    """
    from scipy.optimize import minimize
    C = np.asarray(C, float)
    n = len(C)
    eps = 1e-5
    E = np.eye(2 * n) * eps

    def build(X):
        return geo.build_backbone_batch(X[:, :n], X[:, n:])["CA"]

    def fg(x):
        X = np.vstack([x[None], x[None] + E])
        r = audit.kabsch_rmsd_batch(build(X), C)
        return float(r[0]), (r[1:] - r[0]) / eps

    x0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
    seen = []

    def cb(xk):
        seen.append(np.asarray(xk, float).copy())

    r = minimize(fg, x0, jac=True, method="L-BFGS-B", callback=cb,
                 options={"maxiter": maxiter})
    seen.append(np.asarray(r.x, float).copy())
    out = {}
    for k in snaps:
        x = seen[min(k, len(seen)) - 1]
        out[k] = (x, fg(x)[0])
    return out


def stage_conv(m=CELL_M, resume=True, verbose=True):
    """IS THE OVER-HELICAL BIAS A CONVERGENCE ARTEFACT?  The coordinator's hypothesis.

    A second implementation of the same projection reports 5.5% positive phi where this one
    reports 12.1% (17.5% on the constrained set), and it stops 0.92 A from the coordinate
    average where this one stops 0.785.  The proposed mechanism: the average is a compacted
    non-physical object -- mean CA-CA bond 2.961 A against a native 3.812 -- so the closer
    the projection fits it, the more of that compaction it must absorb as helical torsion,
    and a less converged fit is therefore a more plausible one by accident.

    That is a WITHIN-target claim about the optimisation path, so it needs a within-target
    measurement.  This runs the four generic starts to convergence with NO prior and reads
    the iterate out at 13 iteration counts, giving the whole positive-phi-versus-convergence
    curve for the price of one fit.  At each snapshot the start with the lowest objective
    SO FAR is taken, which is what "multi-start with this maxiter" actually means.

    The comparison that settles the practical question is then direct: early stopping and a
    prior both trade something for plausibility, and this says what each one's price is.
    """
    name = f"project_conv_m{m}.json"
    out = _read(name, {"m": m, "snaps": list(SNAPS), "per": {}}) if resume else \
        {"m": m, "snaps": list(SNAPS), "per": {}}
    per = out["per"]
    tgs = targets()
    todo = [t for t in tgs if t not in per]
    if verbose:
        print(f"{len(todo)}/{len(tgs)} targets for the convergence trace", flush=True)
    for k, pdbid in enumerate(todo):
        t0 = time.time()
        c = cc.load_cache(pdbid)
        seq, nat = c["seq"], c["nat_ca"]
        sub, b, C = cell_of(c, CELL_F, m)
        n = len(C)
        traces = [fit_trace(C, np.full(n, math.radians(a)), np.full(n, math.radians(bb)))
                  for a, bb in STARTS]
        row = {"fold": int(c["fold"]), "n": int(c["n"]),
               "base": sel_of(c["sc"], c["rr"]), "snaps": {}}
        for s in SNAPS:
            j = int(np.argmin([tr[s][1] for tr in traces]))
            x = traces[j][s][0]
            ph, ps = x[:n], x[n:]
            F = geo.build_backbone_batch(ph[None], ps[None])["CA"][0]
            row["snaps"][str(s)] = {
                "d_to_C": float(traces[j][s][1]),
                "rmsd": float(audit.kabsch_rmsd_batch(F[None], nat)[0]),
                "tor": torsion_report(seq, ph, ps),
                "pang": geometry_report(F)["pseudoangle_mean"]}
        per[pdbid] = row
        _write(name, out)
        if verbose and (k + 1) % 20 == 0:
            print(f"  {k+1}/{len(todo)} ({time.time()-t0:.1f}s)", flush=True)
    if len(per) == len(tgs):
        out["summary"] = _summarise_conv(out)
        _write(name, out)
        _print_conv(out)
    return out


def _summarise_conv(out):
    per = out["per"]
    tgs = sorted(per)
    base = np.array([per[t]["base"] for t in tgs], float)
    full = np.array([per[t]["snaps"][str(SNAPS[-1])]["rmsd"] for t in tgs], float)
    S = {"n": len(tgs), "base": float(base.mean()), "snaps": {}}
    for s in out["snaps"]:
        k = str(s)
        v = np.array([per[t]["snaps"][k]["rmsd"] for t in tgs], float)
        d = np.array([per[t]["snaps"][k]["d_to_C"] for t in tgs], float)
        p = np.array([per[t]["snaps"][k]["tor"]["posphi_con_nongly"] for t in tgs], float)
        pa = np.array([per[t]["snaps"][k]["tor"]["posphi_all"] for t in tgs], float)
        g = np.array([per[t]["snaps"][k]["pang"] for t in tgs], float)
        S["snaps"][k] = {"sel": float(v.mean()), "d_to_C": float(d.mean()),
                         "posphi": float(np.nanmean(p)), "posphi_all": float(pa.mean()),
                         "pang": float(g.mean()),
                         "vs_base": _paired(v, base), "vs_full": _paired(v, full)}
    # the within-target question, asked directly: does tightening RAISE positive phi?
    a = np.array([[per[t]["snaps"][str(s)]["d_to_C"] for s in out["snaps"]] for t in tgs])
    b = np.array([[per[t]["snaps"][str(s)]["tor"]["posphi_con_nongly"]
                   for s in out["snaps"]] for t in tgs])
    rho = []
    for i in range(len(tgs)):
        u, v = a[i], b[i]
        ok = np.isfinite(u) & np.isfinite(v)
        if ok.sum() > 3 and np.std(u[ok]) > 1e-12 and np.std(v[ok]) > 1e-12:
            rho.append(float(np.corrcoef(u[ok], v[ok])[0, 1]))
    S["within_rho_d_posphi"] = {"mean": float(np.mean(rho)), "n": len(rho),
                                "median": float(np.median(rho)),
                                "frac_negative": float(np.mean(np.array(rho) < 0))}
    return S


def _print_conv(out):
    S = out["summary"]
    print(f"\n=== CONVERGENCE vs VALIDITY, NO PRIOR — m={out['m']}, {S['n']} targets ===")
    print(f"  baseline {S['base']:.3f}")
    print(f"{'maxiter':>8} {'d to avg':>9} {'selected':>9} {'d vs base':>10} "
          f"{'95% CI':>18} {'phi+ nG':>8} {'phi+ all':>9} {'pang':>7}")
    for k, d in S["snaps"].items():
        p = d["vs_base"]
        print(f"{k:>8} {d['d_to_C']:9.4f} {d['sel']:9.3f} {p['mean_diff']:+10.3f} " +
              f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}]".rjust(19) +
              f" {d['posphi']:8.4f} {d['posphi_all']:9.4f} {d['pang']:7.2f}")
    r = S["within_rho_d_posphi"]
    print(f"  WITHIN-target rho(distance-to-average, positive-phi) over the trace: "
          f"mean {r['mean']:+.3f}, median {r['median']:+.3f}, "
          f"{r['frac_negative']:.2f} of {r['n']} targets negative")
    print("  (the convergence hypothesis predicts a NEGATIVE rho: tighter fit, more "
          "positive phi)")


# ============================================================ stage: ms (the plateau)
def stage_ms(pen_kind, lam, ms=MS, resume=True, verbose=True, multi=True):
    """Is the chosen (penalty, lambda) a plateau in filter size or a lucky cell?"""
    name = f"project_ms_{pen_kind}.json"
    out = _read(name, {"pen": pen_kind, "lam": lam, "ms": list(ms), "multi": multi,
                       "per": {}}) \
        if resume else {"pen": pen_kind, "lam": lam, "ms": list(ms), "multi": multi,
                        "per": {}}
    per = out["per"]
    tgs = targets()
    todo = [t for t in tgs if t not in per]
    if verbose:
        print(f"{len(todo)}/{len(tgs)} targets for the m-plateau "
              f"({pen_kind}, lam={lam})", flush=True)
    for k, pdbid in enumerate(todo):
        t0 = time.time()
        c = cc.load_cache(pdbid)
        seq, fold, nat = c["seq"], int(c["fold"]), c["nat_ca"]
        pen = make_penalty(pen_kind, seq, fold)
        row = {"fold": fold, "n": int(c["n"]), "base": sel_of(c["sc"], c["rr"]),
               "cells": {}}
        for m in ms:
            sub, b, C = cell_of(c, CELL_F, m)
            path = lam_path(C, pen, (0.0, lam), multi=multi)
            row["cells"][str(m)] = {
                "medoid": float(c["rr"][int(sub[b])]),
                "avg": float(audit.kabsch_rmsd_batch(C[None], nat)[0]),
                "fit": float(audit.kabsch_rmsd_batch(path[0.0][0][None], nat)[0]),
                "arm": float(audit.kabsch_rmsd_batch(path[lam][0][None], nat)[0]),
                "posphi": torsion_report(seq, path[lam][1],
                                         path[lam][2])["posphi_con_nongly"],
                "pang": geometry_report(path[lam][0])["pseudoangle_mean"]}
        per[pdbid] = row
        _write(name, out)
        if verbose:
            print(f"[{k+1:3d}/{len(todo)}] {pdbid} ({time.time()-t0:.1f}s)", flush=True)
    if len(per) == len(tgs):
        out["summary"] = _summarise_ms(out)
        _write(name, out)
        _print_ms(out)
    return out


def _summarise_ms(out):
    per = out["per"]
    tgs = sorted(per)
    base = np.array([per[t]["base"] for t in tgs], float)
    S = {"n": len(tgs), "base": float(base.mean()), "cells": {}}
    for m in out["ms"]:
        k = str(m)
        fit = np.array([per[t]["cells"][k]["fit"] for t in tgs], float)
        arm = np.array([per[t]["cells"][k]["arm"] for t in tgs], float)
        S["cells"][k] = {
            "fit": float(fit.mean()), "arm": float(arm.mean()),
            "arm_se": float(arm.std(ddof=1) / math.sqrt(len(arm))),
            "arm_median": float(np.median(arm)),
            "f2": float((arm < 2.0).mean()), "f15": float((arm < 1.5).mean()),
            "fit_vs_base": _paired(fit, base),
            "arm_vs_base": _paired(arm, base), "arm_vs_fit": _paired(arm, fit),
            "posphi": float(np.nanmean([per[t]["cells"][k]["posphi"] for t in tgs])),
            "pang": float(np.mean([per[t]["cells"][k]["pang"] for t in tgs]))}
    return S


def _print_ms(out):
    S = out["summary"]
    print(f"\n=== FILTER-SIZE PLATEAU — {out['pen']}, lam={out['lam']}, "
          f"{S['n']} targets ===")
    print(f"{'m':>5} {'fit':>7} {'arm':>7} {'d(fit,base)':>12} {'d(arm,base)':>12} "
          f"{'95% CI':>18} {'d(arm,fit)':>11} {'95% CI':>18} {'W/L':>9} "
          f"{'phi+':>7} {'pang':>7}")
    for k, d in S["cells"].items():
        a, b_ = d["arm_vs_base"], d["arm_vs_fit"]
        print(f"{k:>5} {d['fit']:7.3f} {d['arm']:7.3f} "
              f"{d['fit_vs_base']['mean_diff']:+12.3f} {a['mean_diff']:+12.3f} " +
              f"[{a['ci95'][0]:+.3f},{a['ci95'][1]:+.3f}]".rjust(19) +
              f" {b_['mean_diff']:+11.3f} " +
              f"[{b_['ci95'][0]:+.3f},{b_['ci95'][1]:+.3f}]".rjust(19) +
              f" {b_['n_better']:3d}/{b_['n_worse']:<4d} {d['posphi']:7.4f} "
              f"{d['pang']:7.2f}")


# ============================================================ stage: init
def stage_init(pen_kind, lam, m=CELL_M, resume=True, verbose=True):
    """Does the STARTING POINT matter once there is a prior?

    The audit found 4 generic starts reproduce the torsion-seeded projection to 0.004 A
    WITHOUT a prior, which says the unconstrained basin is not the issue.  With a prior the
    landscape is different, so this measures four forms on identical inputs:

      `path`  the continuation path this module uses -- generic starts at lam=0, then walk
      `cold`  a cold multi-start AT lam, from the four generic conformations
      `seed`  a single cold solve started from the MEDOID'S OWN torsions
      `coldx` cold multi-start including the medoid's torsions as a fifth start
    """
    name = f"project_init_{pen_kind}.json"
    out = _read(name, {"pen": pen_kind, "lam": lam, "m": m, "per": {}}) \
        if resume else {"pen": pen_kind, "lam": lam, "m": m, "per": {}}
    per = out["per"]
    tgs = targets()
    todo = [t for t in tgs if t not in per]
    if verbose:
        print(f"{len(todo)}/{len(tgs)} targets for the init comparison", flush=True)
    for k, pdbid in enumerate(todo):
        t0 = time.time()
        c = cc.load_cache(pdbid)
        seq, fold, nat = c["seq"], int(c["fold"]), c["nat_ca"]
        pen = make_penalty(pen_kind, seq, fold)
        sub, b, C = cell_of(c, CELL_F, m)
        gi = int(sub[b])
        ex = (c["PHI"][gi], c["PSI"][gi])
        forms = {}
        pth = lam_path(C, pen, (0.0, lam))
        forms["path"] = pth[lam]
        forms["pathx"] = lam_path(C, pen, (0.0, lam), multi=True)[lam]
        forms["cold"] = fit_multi(C, pen=pen, lam=lam)
        forms["seed"] = fit_prior(C, ex[0], ex[1], pen=pen, lam=lam)
        forms["coldx"] = fit_multi(C, pen=pen, lam=lam, extra=ex)
        # the modifier, on the continuation path
        w = agreement_weights(c["W"][sub], C)
        wp = lam_path(C, pen, (0.0, lam), w=w, multi=True)
        forms["wdata"] = wp[lam]                  # weighted data term + the prior
        forms["wdata0"] = wp[0.0]                 # weighted data term ALONE, the control
        row = {"fold": fold, "n": int(c["n"]), "base": sel_of(c["sc"], c["rr"]),
               "fit": float(audit.kabsch_rmsd_batch(pth[0.0][0][None], nat)[0]),
               "forms": {nm: {"rmsd": float(audit.kabsch_rmsd_batch(v[0][None], nat)[0]),
                              "obj": v[3], "d_to_C": v[5],
                              "posphi": torsion_report(seq, v[1],
                                                       v[2])["posphi_con_nongly"]}
                         for nm, v in forms.items()}}
        per[pdbid] = row
        _write(name, out)
        if verbose:
            print(f"[{k+1:3d}/{len(todo)}] {pdbid} ({time.time()-t0:.1f}s)", flush=True)
    if len(per) == len(tgs):
        out["summary"] = _summarise_init(out)
        _write(name, out)
        _print_init(out)
    return out


def _summarise_init(out):
    per = out["per"]
    tgs = sorted(per)
    base = np.array([per[t]["base"] for t in tgs], float)
    fit = np.array([per[t]["fit"] for t in tgs], float)
    S = {"n": len(tgs), "base": float(base.mean()), "fit": float(fit.mean()), "forms": {}}
    for nm in per[tgs[0]]["forms"]:
        v = np.array([per[t]["forms"][nm]["rmsd"] for t in tgs], float)
        o = np.array([per[t]["forms"][nm]["obj"] for t in tgs], float)
        S["forms"][nm] = {"sel": float(v.mean()),
                          "se": float(v.std(ddof=1) / math.sqrt(len(v))),
                          "obj": float(o.mean()),
                          "f2": float((v < 2.0).mean()),
                          "vs_base": _paired(v, base), "vs_fit": _paired(v, fit),
                          "posphi": float(np.nanmean(
                              [per[t]["forms"][nm]["posphi"] for t in tgs]))}
    return S


def _print_init(out):
    S = out["summary"]
    print(f"\n=== INITIALISATION, WITH A PRIOR — {out['pen']}, lam={out['lam']}, "
          f"{S['n']} targets ===")
    print(f"  base {S['base']:.3f}   fit (lam=0) {S['fit']:.3f}")
    print(f"{'form':8} {'sel':>7} {'obj':>8} {'d vs fit':>9} {'95% CI':>18} {'W/L':>9} "
          f"{'phi+':>7}")
    for nm, d in S["forms"].items():
        p = d["vs_fit"]
        print(f"{nm:8} {d['sel']:7.3f} {d['obj']:8.4f} {p['mean_diff']:+9.3f} " +
              f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}]".rjust(19) +
              f" {p['n_better']:3d}/{p['n_worse']:<4d} {d['posphi']:7.4f}")


# ============================================================ stage: valid
#: Ramachandran basins, in degrees, as (phi_lo, phi_hi, psi_lo, psi_hi).  Conventional
#: boxes; `other` is whatever falls in none of them.  They exist for READABILITY -- the
#: quantitative distribution comparison is the symmetric KL below, which uses no boxes.
BASINS = {"alphaR": (-160, -20, -120, 50),
          "beta":   (-180, -20, 50, 180),
          "alphaL": (20, 125, -30, 100)}


def basin_of(phi, psi):
    """``(n,)`` basin label index per residue; ``len(BASINS)`` means none of them."""
    a = np.degrees(_wrap(phi))
    b = np.degrees(_wrap(psi))
    out = np.full(a.shape, len(BASINS), int)
    for k, (p0, p1, s0, s1) in enumerate(BASINS.values()):
        hit = (a >= p0) & (a <= p1) & (b >= s0) & (b <= s1)
        out = np.where(hit & (out == len(BASINS)), k, out)
    return out


def _hist2d(phi, psi, w=None):
    H = np.zeros((RB, RB))
    np.add.at(H, (_bin(phi), _bin(psi)), 1.0 if w is None else w)
    return H


def _symkl(P, Q, eps=1e-9):
    """Symmetric KL between two (RB, RB) histograms, normalised, in nats."""
    P = P / max(P.sum(), eps) + eps
    Q = Q / max(Q.sum(), eps) + eps
    P, Q = P / P.sum(), Q / Q.sum()
    return float((P * np.log(P / Q)).sum() + (Q * np.log(Q / P)).sum())


def stage_valid(pen_kind, lam, m=CELL_M, verbose=True, tag="x"):
    """The full geometry + torsion validity table: has the prior broken what worked?

    Unlike the sweep stages this RECOMPUTES the chosen arm rather than reading summary
    statistics back, because the distribution comparison needs the emitted torsions
    themselves and the sweep only persisted their summaries.

    Five columns throughout, and the comparators are the point: `arm` and `fit` are what
    the pipeline emits with and without the prior; `avg` is the coordinate average they are
    both projections of; `win` is the medoid LIBRARY WINDOW -- a real fragment of a real
    protein, and the fairest reference for a construction of this length; `nat` is the
    deposited native.

    What is exact by construction and therefore not measured: every bond length
    (BOND_N_CA / BOND_CA_C / BOND_C_N), every bond angle (ANGLE_N_CA_C / ANGLE_CA_C_N /
    ANGLE_C_N_CA) and omega, which `build_backbone_batch` fixes at OMEGA_TRANS for every
    residue of every structure it builds.  The virtual CA-CA bond is asserted below to
    1e-9 rather than reported as a distribution, because a distribution of a constant is
    not informative -- what IS informative is that the prior has not perturbed it.
    """
    tgs = targets()
    rows = {"fit": [], "arm": [], "nat": [], "win": [], "avg": []}
    tor = {"fit": [], "arm": [], "win": []}
    H = {k: np.zeros((RB, RB)) for k in ("fit", "arm", "win")}
    Bs = {k: np.zeros(len(BASINS) + 1) for k in ("fit", "arm", "win")}
    worst_bond = 0.0
    per = {}
    for k, t in enumerate(tgs):
        c = cc.load_cache(t)
        seq = c["seq"]
        n = len(seq)
        sub, b, C = cell_of(c, CELL_F, m)
        gi = int(sub[b])
        pen = make_penalty(pen_kind, seq, int(c["fold"]))
        path = lam_path(C, pen, (0.0, lam), multi=True)
        got = {"fit": path[0.0], "arm": path[lam]}
        rows["nat"].append(geometry_report(c["nat_ca"]))
        rows["win"].append(geometry_report(c["W"][gi]))
        rows["avg"].append(geometry_report(C))
        #: the residues whose (phi, psi) PAIR is geometrically live -- 1..n-2.  The
        #: distribution comparison must use the same set the penalty does, or it scores
        #: two free parameters that no builder ever reads.
        sl = slice(1, n - 1)
        for nm, v in got.items():
            rows[nm].append(geometry_report(v[0]))
            tor[nm].append(torsion_report(seq, v[1], v[2]))
            H[nm] += _hist2d(v[1][sl], v[2][sl])
            np.add.at(Bs[nm], basin_of(v[1][sl], v[2][sl]), 1.0)
            st = np.linalg.norm(v[0][1:] - v[0][:-1], axis=1)
            worst_bond = max(worst_bond, float(np.abs(st - st.mean()).max()))
        tor["win"].append(torsion_report(seq, c["PHI"][gi], c["PSI"][gi]))
        H["win"] += _hist2d(c["PHI"][gi][sl], c["PSI"][gi][sl])
        np.add.at(Bs["win"], basin_of(c["PHI"][gi][sl], c["PSI"][gi][sl]), 1.0)
        per[t] = {"arm_rmsd": float(audit.kabsch_rmsd_batch(got["arm"][0][None],
                                                            c["nat_ca"])[0]),
                  "fit_rmsd": float(audit.kabsch_rmsd_batch(got["fit"][0][None],
                                                            c["nat_ca"])[0])}
        if verbose and (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tgs)}", flush=True)
    #: the REAL reference distribution: the fold-disjoint training corpus, pooled.
    Href = np.asarray(build_tables(verbose=False)["cls"], float).sum((0, 1))
    Bref = np.zeros(len(BASINS) + 1)
    ii, jj = np.meshgrid((np.arange(RB) + 0.5) / RB * 360 - 180,
                         (np.arange(RB) + 0.5) / RB * 360 - 180, indexing="ij")
    np.add.at(Bref, basin_of(np.radians(ii), np.radians(jj)).ravel(), Href.ravel())
    names = list(BASINS) + ["other"]
    out = {"pen": pen_kind, "lam": lam, "m": m, "n": len(tgs),
           "worst_bond_deviation": worst_bond,
           "geom": {k: _agg(v) for k, v in rows.items()},
           "tors": {k: {kk: float(np.nanmean([r[kk] for r in v])) for kk in v[0]}
                    for k, v in tor.items()},
           "symkl_vs_corpus": {k: _symkl(H[k], Href) for k in H},
           "basins": {k: dict(zip(names, (Bs[k] / max(Bs[k].sum(), 1)).tolist()))
                      for k in Bs},
           "basins_ref": dict(zip(names, (Bref / Bref.sum()).tolist())),
           "per_target": per}
    assert worst_bond < 1e-9, worst_bond
    _write(f"project_valid_{pen_kind}.json", out)
    _print_valid(out)
    return out


def _print_valid(out):
    print(f"\n=== VALIDITY — {out['pen']}, lam={out['lam']}, {out['n']} targets ===")
    ks = ["step_mean", "step_sd", "pseudoangle_mean", "frac_pseudoangle_out_of_75_150",
          "min_nonlocal_CA", "frac_nonlocal_under_4A", "frac_dihedral_positive", "rg"]
    print(f"{'quantity':34} {'arm':>10} {'fit':>10} {'avg':>10} {'window':>10} "
          f"{'native':>10}")
    for k in ks:
        g = out["geom"]
        print(f"{k:34} {g['arm'][k]:10.4f} {g['fit'][k]:10.4f} {g['avg'][k]:10.4f} "
              f"{g['win'][k]:10.4f} {g['nat'][k]:10.4f}")
    print(f"{'--- torsions ---':34}")
    for k in ("posphi_all", "posphi_all_nongly", "posphi_con", "posphi_con_nongly"):
        t = out["tors"]
        print(f"{k:34} {t['arm'][k]:10.4f} {t['fit'][k]:10.4f} {'':>10} "
              f"{t['win'][k]:10.4f} {'':>10}")
    print(f"{'--- (phi,psi) distribution ---':34}")
    k = out["symkl_vs_corpus"]
    print(f"{'sym KL vs training corpus':34} {k['arm']:10.4f} {k['fit']:10.4f} "
          f"{'':>10} {k['win']:10.4f} {'':>10}")
    for nm in out["basins_ref"]:
        b = out["basins"]
        print(f"{'occupancy ' + nm:34} {b['arm'][nm]:10.4f} {b['fit'][nm]:10.4f} "
              f"{'':>10} {b['win'][nm]:10.4f} {out['basins_ref'][nm]:10.4f}")
    print(f"{'  (last column = training corpus)':34}")
    print(f"  bonds/angles/omega are EXACT by construction; worst CA-CA deviation "
          f"within a chain {out['worst_bond_deviation']:.2e} A")


# ============================================================ stage: conc
def stage_conc(pen_kind, lam, m=CELL_M, verbose=True, tag="x"):
    """The audit's two fragility findings, applied to THIS arm.

    (1) Is the mean carried by ten targets?  Sorted per-target gains, the median gain, and
        the mean recomputed with the top-10 contributors removed.
    (2) Is the gain confined to targets whose POOL contains a sub-2 A candidate?  The audit
        found `medoid75` gains only on those 73 and LOSES 0.096 A on the other 53.
    """
    sw = _read(f"project_sweep{tag}_m{m}.json")
    assert sw and sw.get("summary"), "run `sweep` to completion first"
    per = sw["per"]
    tgs = sorted(per)
    key = f"{lam:g}"
    base = np.array([per[t]["base"] for t in tgs], float)
    fit = np.array([per[t]["arms"][pen_kind]["0"]["rmsd"] for t in tgs], float)
    arm = np.array([per[t]["arms"][pen_kind][key]["rmsd"] for t in tgs], float)
    pb = np.array([per[t]["pool_best"] for t in tgs], float)

    def conc(v, ref, nm):
        d = v - ref
        o = np.argsort(d)
        return {"name": nm, "mean": float(d.mean()), "median": float(np.median(d)),
                "top10_share": float(d[o[:10]].sum() / d.sum()) if d.sum() != 0 else
                float("nan"),
                "mean_ex_top10": float(d[o[10:]].mean()),
                "n_better": int((d < 0).sum()), "n_worse": int((d > 0).sum()),
                "worst_targets": [[tgs[i], float(d[i])] for i in o[-5:][::-1]],
                "best_targets": [[tgs[i], float(d[i])] for i in o[:5]]}

    near = pb < 2.0
    out = {"pen": pen_kind, "lam": lam, "m": m, "n": len(tgs),
           "arm_vs_base": conc(arm, base, "arm vs base"),
           "arm_vs_fit": conc(arm, fit, "arm vs fit"),
           "fit_vs_base": conc(fit, base, "fit vs base"),
           "near": {"n": int(near.sum()),
                    "arm_vs_base": _paired(arm[near], base[near]),
                    "arm_vs_fit": _paired(arm[near], fit[near]),
                    "fit_vs_base": _paired(fit[near], base[near])},
           "far": {"n": int((~near).sum()),
                   "arm_vs_base": _paired(arm[~near], base[~near]),
                   "arm_vs_fit": _paired(arm[~near], fit[~near]),
                   "fit_vs_base": _paired(fit[~near], base[~near])}}
    # per-fold, the leave-fold-out discipline check
    folds = np.array([per[t]["fold"] for t in tgs], int)
    out["per_fold"] = {int(f): {"n": int((folds == f).sum()),
                                "arm_vs_base": float((arm - base)[folds == f].mean()),
                                "arm_vs_fit": float((arm - fit)[folds == f].mean())}
                       for f in sorted(set(folds.tolist()))}
    _write(f"project_conc_{pen_kind}.json", out)
    _print_conc(out)
    return out


def _print_conc(out):
    print(f"\n=== CONCENTRATION — {out['pen']}, lam={out['lam']}, {out['n']} targets ===")
    for k in ("fit_vs_base", "arm_vs_base", "arm_vs_fit"):
        d = out[k]
        print(f"  {d['name']:14} mean {d['mean']:+.3f}  median {d['median']:+.3f}  "
              f"top-10 share {d['top10_share']:+.3f}  mean ex-top10 "
              f"{d['mean_ex_top10']:+.3f}  {d['n_better']}/{d['n_worse']}")
    for grp in ("near", "far"):
        g = out[grp]
        lbl = "pool best < 2 A" if grp == "near" else "pool best >= 2 A"
        print(f"  {lbl:16} n={g['n']:3d}  "
              f"fit {g['fit_vs_base']['mean_diff']:+.3f}  "
              f"arm {g['arm_vs_base']['mean_diff']:+.3f} "
              f"[{g['arm_vs_base']['ci95'][0]:+.3f},{g['arm_vs_base']['ci95'][1]:+.3f}]  "
              f"arm vs fit {g['arm_vs_fit']['mean_diff']:+.3f} "
              f"[{g['arm_vs_fit']['ci95'][0]:+.3f},{g['arm_vs_fit']['ci95'][1]:+.3f}]")
    print("  per fold (arm vs base / arm vs fit): " + "  ".join(
        f"{f}:{d['arm_vs_base']:+.3f}/{d['arm_vs_fit']:+.3f}"
        for f, d in out["per_fold"].items()))


# ============================================================ stage: leak
def stage_leak(pen_kind="rama", lam=0.03, n_targets=12, m=CELL_M):
    """NaN-poison the native and assert every returned structure is BIT-IDENTICAL.

    The subset is recomputed from the poisoned cache too, so a filter that read a native
    would change the subset and be caught, not just a penalty that did.
    """
    worst, bad = 0.0, []
    for pdbid in targets()[:n_targets]:
        c = cc.load_cache(pdbid)
        p = dict(c)
        p["rr"] = np.full_like(c["rr"], np.nan)
        p["nat_ca"] = np.full_like(c["nat_ca"], np.nan)
        for src in (c, p):
            pass
        outs = []
        for src in (c, p):
            sub, b, C = cell_of(src, CELL_F, m)
            pen = make_penalty(pen_kind, src["seq"], int(src["fold"]))
            path = lam_path(C, pen, (0.0, lam), multi=True)
            outs.append((sub.copy(), path[0.0][0], path[lam][0], path[lam][1]))
        if not np.array_equal(outs[0][0], outs[1][0]):
            bad.append((pdbid, "subset"))
            continue
        for a, b2 in zip(outs[0][1:], outs[1][1:]):
            d = float(np.abs(np.asarray(a) - np.asarray(b2)).max())
            worst = max(worst, d)
            if d != 0.0:
                bad.append((pdbid, d))
    out = {"n_targets": n_targets, "pen": pen_kind, "lam": lam,
           "worst_abs_diff": worst, "bad": bad}
    _write("project_leak.json", out)
    print(f"\n=== LEAKAGE AUDIT === {n_targets} targets, penalty {pen_kind}, lam {lam}: "
          f"worst |diff| {worst:.3e}, {len(bad)} failures")
    assert not bad, bad
    return out


# ============================================================ stage: dev
#: THE arm carried to dev-24, pinned here BEFORE `stage_dev` was ever run and asserted by
#: `test_project.py`.  Filled by `pin_dev_arm` from the tuning result, then frozen.
DEV_ARM = _read("project_devarm.json")

#: The band of positive-phi rates (constrained set, non-glycine) that counts as REAL.
#: Its ends are the two extreme measurements this module made of real backbones -- 2.90%
#: over the whole fold-disjoint corpus and 5.66% over `peptide_db` chains alone, with the
#: library windows the pools actually draw from at 5.40% in between.  S7 established that
#: driving the rate to exactly ZERO overshoots real protein, so the target is this band and
#: not the origin, and an arm that lands below 2.90% is as mis-calibrated as one above
#: 5.66%.
REAL_POSPHI = (0.029, 0.0566)


def choose_arm(tag="h", m=CELL_M, band=REAL_POSPHI):
    """THE PRE-DECLARED SELECTION RULE, written down before the sweep it reads finished.

    Among cells whose positive-phi rate lands INSIDE `REAL_POSPHI`, take the lowest
    selected CA-RMSD on the 126-target tuning instrument; break ties toward the smaller
    lambda, because a smaller lambda perturbs the incumbent objective less.  Cells outside
    the band are not eligible at any RMSD -- an arm that fixes the number by overshooting
    it has not fixed anything.

    Returns ``(penalty, lam, sel, posphi)`` or None if no cell qualifies.
    """
    o = _read(f"project_sweep{tag}_m{m}.json")
    if not o or not o.get("summary"):
        return None
    best = None
    for key, d in o["summary"]["cells"].items():
        pk, lam = key.split("|")
        if not (band[0] <= d["posphi_con_nongly"] <= band[1]):
            continue
        cand = (round(d["sel"], 6), float(lam), pk, d["posphi_con_nongly"])
        if best is None or cand[:2] < best[:2]:
            best = cand
    return None if best is None else (best[2], best[1], best[0], best[3])


#: The non-inferiority margin on the dev pass's SECONDARY endpoint, in Angstroms.
#: 0.05 A is an order of magnitude above what the tuning set bounds the effect at
#: (|d| < 0.005 A over four filter sizes and 126 paired targets) and an order of magnitude
#: below the gain the arm must not damage (S8-11's 0.250 A).
DEV_MARGIN = 0.05


def pin_dev_arm(pen_kind, lam, m, note=""):
    """Write the pre-registration file.  Refuses to overwrite an existing pin.

    THIS IS THE THIRD DEV-24 PASS ON THE `fit` FAMILY (S8-8, S8-11, this), and it is
    declared as such.  What makes it defensible is that it is not a third attempt at the
    same question: S8-8 and S8-11 spent dev on ACCURACY claims, and this pass does not make
    one.  The tuning instrument has already bounded the accuracy effect at |d| < 0.005 A
    across four filter sizes on 126 paired targets, and dev-24 -- whose SE on a mean is
    ~0.25 A -- cannot resolve that and is not being asked to.

    Endpoints, fixed here before the pass runs:

      PRIMARY, and the only thing the pass can actually decide: does the positive-phi
        correction TRANSFER?  17.5% on tuning without the prior, 5.5% with it.  This is a
        deterministic property of the emitted structure, so it is the endpoint dev has
        power for.
      SECONDARY, non-inferiority: is the paired RMSD difference against `fit`, computed in
        the SAME pass on the SAME targets from the SAME pools, inside +-`DEV_MARGIN`?
      NOT AN ENDPOINT: any accuracy gain.  No claim of one will be made from this pass
        whatever it returns, and the tuning result already says there is none to find.
    """
    if os.path.exists(_p("project_devarm.json")):
        raise RuntimeError("dev arm already pinned; it is not re-negotiable")
    _write("project_devarm.json", {
        "pen": pen_kind, "lam": lam, "m": m, "multi": True,
        "pinned_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "pass_number_on_family": 3,
        "primary_endpoint": "positive-phi (constrained, non-Gly) transfers "
                            "from 0.175 toward the 0.029-0.0566 real band",
        "secondary_endpoint": f"paired d(arm, fit) within +-{DEV_MARGIN} A "
                              f"(non-inferiority, not superiority)",
        "not_an_endpoint": "any accuracy gain; tuning bounds it at |d| < 0.005 A",
        "note": note})
    return _read("project_devarm.json")


def stage_dev():
    """THE SINGLE DEV PASS.  One pre-registered arm, no iteration.

    Dev pools are `s8/inband_devpool`, built and asserted by `s8.inband.stage_devpool`.
    They carry no torsions, so the projection runs from the generic starts only -- which is
    the deployable form anyway.
    """
    arm = _read("project_devarm.json")
    assert arm, "pin the arm with `pin_dev_arm` before running dev"
    pen_kind, lam, m = arm["pen"], arm["lam"], arm["m"]
    multi = arm.get("multi", True)
    inband.stage_devpool(verbose=False)
    nats = {q.pdb: np.asarray(q.ca, float) for q in debias.dev_targets()}
    rows = []
    for pdbid in sorted(x[:-4] for x in os.listdir(cc.DCACHE) if x.endswith(".npz")):
        z = np.load(os.path.join(cc.DCACHE, f"{pdbid}.npz"), allow_pickle=True)
        W = z["W"].astype(float)
        rr = z["rr"].astype(float)
        n = int(z["n"])
        seq = str(z["seq"])
        fold = int(z["fold"])
        pred = poolsize.load_pred(pdbid)
        D = gen.pair_D(W, n)
        sc = gen.shipped_score(pred, D)
        P = np.zeros((len(W), len(W)), float)
        for a in range(len(W)):
            P[a] = audit.kabsch_rmsd_batch(W, W[a])
        c = {"pdb": pdbid, "n": n, "fold": fold, "seq": seq, "W": W, "D": D, "sc": sc,
             "typ": np.abs(D - D.mean(0)[None, :]).mean(1), "P": P}
        sub, b, C = cell_of(c, CELL_F, m)
        pen = make_penalty(pen_kind, seq, fold)
        path = lam_path(C, pen, (0.0, lam), multi=multi)
        nat = nats[pdbid]
        rows.append({"pdb": pdbid, "n": n, "fold": fold,
                     "best": float(rr.min()), "score": sel_of(sc, rr),
                     "medoid": float(rr[int(sub[b])]),
                     "fit": float(audit.kabsch_rmsd_batch(path[0.0][0][None], nat)[0]),
                     "arm": float(audit.kabsch_rmsd_batch(path[lam][0][None], nat)[0]),
                     "posphi": torsion_report(seq, path[lam][1],
                                              path[lam][2])["posphi_con_nongly"],
                     "posphi_fit": torsion_report(seq, path[0.0][1],
                                                  path[0.0][2])["posphi_con_nongly"]})
    base = np.array([r["score"] for r in rows], float)
    fit = np.array([r["fit"] for r in rows], float)
    v = np.array([r["arm"] for r in rows], float)
    med = np.array([r["medoid"] for r in rows], float)
    out = {"n_targets": len(rows), "arm": f"{CELL_F}|{m}:fit+{pen_kind}@{lam}",
           "per_target": rows,
           "pool_best": float(np.mean([r["best"] for r in rows])),
           "score": float(base.mean()), "medoid": float(med.mean()),
           "fit": float(fit.mean()),
           "fit_vs_base": _paired(fit, base),
           "arm_stats": {"sel": float(v.mean()),
                         "se": float(v.std(ddof=1) / math.sqrt(len(v))),
                         "median": float(np.median(v)),
                         "f2": float((v < 2.0).mean()),
                         "vs_base": _paired(v, base), "vs_fit": _paired(v, fit)},
           "posphi_arm": float(np.nanmean([r["posphi"] for r in rows])),
           "posphi_fit": float(np.nanmean([r["posphi_fit"] for r in rows]))}
    _write("project_dev.json", out)
    _print_dev(out)
    return out


def _print_dev(out):
    print(f"\n=== DEV-24, ONE PRE-REGISTERED PASS — {out['arm']} ===")
    print(f"  shipped score {out['score']:.3f}   pool best {out['pool_best']:.3f}   "
          f"medoid {out['medoid']:.3f}")
    p, q = out["fit_vs_base"], out["arm_stats"]["vs_base"]
    r = out["arm_stats"]["vs_fit"]
    print(f"  fit (lam=0)  {out['fit']:.3f}   d {p['mean_diff']:+.3f} "
          f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}]  {p['n_better']}/{p['n_worse']}")
    print(f"  ARM          {out['arm_stats']['sel']:.3f} +- {out['arm_stats']['se']:.3f}"
          f"   d vs base {q['mean_diff']:+.3f} "
          f"[{q['ci95'][0]:+.3f},{q['ci95'][1]:+.3f}]  {q['n_better']}/{q['n_worse']}"
          f"   d vs fit {r['mean_diff']:+.3f} "
          f"[{r['ci95'][0]:+.3f},{r['ci95'][1]:+.3f}]  {r['n_better']}/{r['n_worse']}")
    print(f"  pos-phi (constrained, non-Gly)  fit {out['posphi_fit']:.4f} -> "
          f"arm {out['posphi_arm']:.4f}")


# ============================================================ report
def stage_report():
    for nm, fn in (("project_prior_report.json", _print_prior),
                   (f"project_sweep_m{CELL_M}.json", lambda o: _print_sweep(o)),
                   (f"project_sweepx_m{CELL_M}.json", lambda o: _print_sweep(o)),
                   (f"project_sweeph_m{CELL_M}.json", lambda o: _print_sweep(o)),
                   ("project_leak.json", None)):
        o = _read(nm)
        if o is None:
            continue
        if fn is None:
            print(f"\n=== LEAKAGE === worst |diff| {o['worst_abs_diff']:.3e}, "
                  f"{len(o['bad'])} failures")
        elif nm.startswith("project_sweep") and not o.get("summary"):
            print(f"\n(sweep incomplete: {len(o['per'])} targets)")
        else:
            fn(o)
    for pat, fn in (("project_ms_", _print_ms), ("project_init_", _print_init),
                    ("project_valid_", _print_valid), ("project_conc_", _print_conc)):
        for f in sorted(os.listdir(HERE)):
            if f.startswith(pat) and f.endswith(".json"):
                o = _read(f)
                if o and (o.get("summary") or pat.startswith("project_valid")
                          or pat.startswith("project_conc")):
                    fn(o)
    o = _read("project_dev.json")
    if o:
        _print_dev(o)


STAGES = {"prior": stage_prior, "sweep": stage_sweep, "sweepx": stage_sweep,
          "sweeph": stage_sweep, "conv": stage_conv, "ms": stage_ms,
          "init": stage_init, "valid": stage_valid, "conc": stage_conc,
          "leak": stage_leak, "dev": stage_dev, "report": stage_report}


def main(argv):
    if not argv or argv[0] not in STAGES:
        print(__doc__)
        return 1
    cmd, rest = argv[0], argv[1:]
    fn = STAGES[cmd]
    if cmd in ("ms", "init", "valid", "conc"):
        fn(rest[0], float(rest[1]))
    elif cmd == "sweep" and rest:
        stage_sweep(pens=tuple(rest[0].split(",")), m=int(rest[1]) if len(rest) > 1
                    else CELL_M)
    elif cmd == "sweepx":
        stage_sweep(pens=tuple(rest[0].split(",")) if rest else PENALTIES,
                    lams=LAMS_X, multi=True, tag="x")
    elif cmd == "sweeph":
        stage_sweep(pens=tuple(rest[0].split(",")) if rest else ("phip", "ramah"),
                    lams=LAMS_H, multi=True, tag="h")
    else:
        fn()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
