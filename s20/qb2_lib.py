"""SPRINT 20 / WORKSTREAM B -- machinery for the LEGACY-vs-AMBER LANDSCAPE question (Q1).

Pre-registered in `s20/PREREG_B.md` before this file existed.

WHAT IS HELD CONSTANT AND WHAT CHANGES
--------------------------------------
Sequence, torsion representation (CONTINUOUS `z = (phi, psi) in T^{2n}`), candidate set, ansatz,
qubit count, initialisation, optimiser, evaluation budget, CVaR rule, measurement budget and
convergence threshold are all held fixed.  **Only the Hamiltonian changes.**

    LEG   core.energy.components_batch . DEFAULT_WEIGHTS      (genuine Legacy, never fitted)
    AMB   ff14SB / GBn2 SINGLE POINT, no minimisation         (genuine AMBER)
    AMBc  sign(E) * log1p(|E|) of the SAME AMB call           (a MONOTONE reparameterisation)
    DIST  the deployed distogram functional                   (reference axis only)

`AMBc` is the load-bearing control.  A strictly monotone transform leaves the argmin, the full
ranking of configurations and every level set UNCHANGED; it changes only gradient magnitudes,
curvature and dynamic range.  So any AMB/AMBc difference is EXACTLY not a property of the
optimisation problem's solution set.  That is an identity, labelled EXACT, not a finding.

THE SCALE CONTROL (PREREG section 3, BRIEF section 6 rule 1).  Every landscape metric is computed
on the STANDARDISED field  Ehat = (E - mu_ref)/sd_ref,  with mu_ref/sd_ref measured on the
target's own K=500 BLOSUM retrieval pool.  Without it, comparing a gradient norm of an energy in
units of 1e8 with one in units of 1e1 is a units comparison, not a physics one.

THE LEAN AMBER CALL.  `core.amber._run` does a restraint-parameter loop, a getPositions, a
_split and a restraint-RMSD on every call; for a pure single point none of that is needed.
`AmberSP` reuses the same `AmberHamiltonian` context and is asserted BIT-EXACT (0.00e+00
relative difference) against `core.amber.refine_coords(k_restraint=0, steps=-1)` on every target
before any number is read (`AmberSP.verify`).  It is ~9x faster and is arithmetically the same
call.

BUDGET.  One budget unit = one evaluation of the objective at one continuous configuration,
finite-difference probes INCLUDED.  Identical for every arm.
"""
from __future__ import annotations

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

from s12 import instrument as I                 # noqa: E402
from s13 import qarch_lib as QA                 # noqa: E402
from s15 import seed as SD                      # noqa: E402
from s16 import energy_lib as EL                # noqa: E402
from s19 import qb_lib as QB                    # noqa: E402
from core import geometry as geo                # noqa: E402

SALT = "s20qb"
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

FD_H = 0.02          # radians, central-difference step (PREREG section 3)
N_SUBSET = 20


def wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


# ============================================================ target subset
def subset(n=N_SUBSET):
    """Deterministic fold-balanced, length-spread subset of the 126 tuning targets.

    Drawn once with `stable_rng` and persisted in every artefact.  benchmark60 and dev24 are
    never read.
    """
    tg = I.targets()
    folds = sorted({int(t["fold"]) for t in tg})
    per = max(1, n // len(folds))
    out = []
    for f in folds:
        cand = sorted([t for t in tg if int(t["fold"]) == f], key=lambda t: (t["n"], t["pdb"]))
        # length-spread: take `per` evenly spaced ranks in the length ordering
        idx = np.linspace(0, len(cand) - 1, per).round().astype(int)
        out += [cand[k] for k in dict.fromkeys(idx)]
    # top up deterministically if rounding left us short
    if len(out) < n:
        rng = SD.stable_rng("subset", n, salt=SALT)
        rest = [t for t in tg if t["pdb"] not in {o["pdb"] for o in out}]
        rest = sorted(rest, key=lambda t: t["pdb"])
        out += [rest[k] for k in rng.permutation(len(rest))[:n - len(out)]]
    return sorted(out, key=lambda t: t["pdb"])[:n]


# ============================================================ the AMBER single point
class AmberSP:
    """Genuine ff14SB/GBn2 single point from CONTINUOUS torsions, no minimisation.

    Bit-exact against `core.amber.refine_coords(k_restraint=0.0, steps=-1)` -- asserted by
    `verify()`, which is called on every target before any number is read (BRIEF section 10).
    """

    def __init__(self, seq, rep):
        from core import amber as am
        self.am = am
        self.seq = seq
        self.rep = rep
        self.H = am.builder_for(seq, rep, "CPU", 1)
        self.n_calls = 0

    def _e(self, coords):
        from openmm import unit
        H = self.H
        pos = H._assemble(H._heavy_positions(coords, chi1=None))
        ctx = H.context
        ctx.setPositions(pos * unit.nanometer)
        ctx.setParameter("k_rest", 0.0)
        return float(ctx.getState(getEnergy=True).getPotentialEnergy()
                     .value_in_unit(unit.kilocalorie_per_mole))

    def batch(self, phi, psi):
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        c = geo.build_backbone_batch(phi, psi)
        out = np.empty(len(phi))
        for b in range(len(phi)):
            out[b] = self._e({k: v[b] for k, v in c.items()})
        self.n_calls += len(phi)
        return out

    def verify(self, phi, psi, m=4):
        """BIT-EXACTNESS GATE.  Returns the max relative difference actually measured, and
        the number of comparisons made -- a gate that fires zero times is not evidence
        (BRIEF section 6 rule 3)."""
        phi = np.atleast_2d(np.asarray(phi, float))[:m]
        psi = np.atleast_2d(np.asarray(psi, float))[:m]
        c = geo.build_backbone_batch(phi, psi)
        worst, k = 0.0, 0
        for b in range(len(phi)):
            cc = {kk: v[b] for kk, v in c.items()}
            ref = float(self.am.refine_coords(self.seq, self.rep, cc, k_restraint=0.0,
                                              steps=-1, tolerance=1e9, threads=1,
                                              memo=False)["energy"])
            got = self._e(cc)
            worst = max(worst, abs(ref - got) / max(1.0, abs(ref)))
            k += 1
        return worst, k


# ============================================================ the four Hamiltonians
class Ham:
    """One Hamiltonian with a hard, identically-counted evaluation budget.

    `kind` in {"LEG", "AMB", "AMBc", "DIST"}.  `raw()` is unbudgeted and unrecorded (post-hoc
    scoring and the free reference-ensemble calibration only); `__call__` is budgeted and
    records the emitted set.
    """

    def __init__(self, kind, tgt, budget=10 ** 9, keep=False, sp=None):
        self.kind = kind
        self.t = tgt
        self.budget = int(budget)
        self.used = 0
        self.keep = keep
        self.sp = sp
        self._phi, self._psi, self._e = [], [], []
        self.mu_ref = 0.0
        self.sd_ref = 1.0
        self.med_ref = 0.0
        self.iqr_ref = 1.0

    @property
    def left(self):
        return max(0, self.budget - self.used)

    # -- the genuine energies -------------------------------------------------
    def raw(self, phi, psi):
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        k = self.kind
        if k == "LEG":
            comp = EL.legacy_components_of_windows(self.t["seq"], phi, psi)
            return np.asarray(EL.legacy_total_from(comp), float)
        if k in ("AMB", "AMBc"):
            e = self.sp.batch(phi, psi)
            if k == "AMBc":
                e = np.sign(e) * np.log1p(np.abs(e))
            return e
        if k == "DIST":
            from core import project as pj
            CA = np.asarray(pj.build_ca_exact(phi, psi), float)
            rv = CA[:, self.t["i"], :] - CA[:, self.t["j"], :]
            d = np.sqrt((rv * rv).sum(-1))
            return ((d - self.t["dhat"]) ** 2 * self.t["inv2"]).sum(-1)
        raise ValueError(k)

    def __call__(self, phi, psi):
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        k = min(len(phi), self.left)
        if k <= 0:
            return np.zeros(0)
        phi, psi = phi[:k], psi[:k]
        e = self.raw(phi, psi)
        self.used += k
        if self.keep:
            self._phi.append(phi.copy()); self._psi.append(psi.copy()); self._e.append(e)
        return e

    def seen(self):
        if not self._phi:
            return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros(0)
        return (np.concatenate(self._phi), np.concatenate(self._psi),
                np.concatenate(self._e))

    # -- the scale control ----------------------------------------------------
    def calibrate(self, PHI, PSI, e=None):
        """FREE, unbudgeted.  The scale control, on the target's own K=500 retrieval pool.

        PRE-REGISTRATION DEVIATION, DECLARED (PREREG_B section 3 specified mean/sd).  Measured
        on the first target: raw `AMB` over 500 REAL retrieval windows has
        `sd_ref = 1.0e12` and a range of **13.2 decades** -- an r^-12 clash tail.  Dividing by
        that sd sends every other structure's standardised value, gradient norm and Hessian
        eigenvalue to ~1e-9, i.e. the mean/sd control is NOT A SCALE for this variable and would
        have silently reported "AMBER's landscape is flat".  A ROBUST centre/scale
        (median, IQR/1.349) is substituted as the primary and BOTH are persisted; the mean/sd
        pair and `log10_range` are reported as the diagnostic that forced the substitution.
        The pre-registration is kept unedited and this deviation is recorded in the findings.

        Non-finite values (AMBER's collapse sentinel) are dropped and counted, because a gate
        that never fires is not evidence."""
        e = self.raw(PHI, PSI) if e is None else np.asarray(e, float)
        ok = np.isfinite(e)
        ev = e[ok]
        self.mu_ref = float(np.mean(ev)) if ok.any() else 0.0
        self.sd_ref = float(np.std(ev)) if ok.any() else 1.0
        if not np.isfinite(self.sd_ref) or self.sd_ref <= 0:
            self.sd_ref = 1.0
        if ok.sum() > 3:
            q1, q2, q3 = np.percentile(ev, [25, 50, 75])
            self.med_ref = float(q2)
            self.iqr_ref = float(max((q3 - q1) / 1.349, 1e-12))
        else:
            self.med_ref, self.iqr_ref = 0.0, 1.0
        return {"mu_ref": self.mu_ref, "sd_ref": self.sd_ref,
                "med_ref": self.med_ref, "iqr_ref": self.iqr_ref,
                "n_ref": int(len(e)), "n_nonfinite_ref": int((~ok).sum()),
                "log10_range": float(np.log10(np.ptp(ev) + 1e-30)) if ok.sum() > 1 else float("nan"),
                "sd_over_iqr": float(self.sd_ref / self.iqr_ref),
                "skew_p99_over_p50": float(np.percentile(np.abs(ev), 99)
                                           / max(1e-30, np.percentile(np.abs(ev), 50)))
                if ok.sum() > 3 else float("nan"),
                "ref_std_max": float(np.max(self.std(ev))) if ok.any() else 6.0,
                "ref_std_p99": float(np.percentile(self.std(ev), 99)) if ok.sum() > 3 else 6.0}

    def std(self, e):
        """The PRIMARY scale control: robust standardisation.  Scale-invariant landscape
        metrics (frac_neg, cond, aniso, n_localmin, acorr_len, tv/range) do not depend on this
        choice at all; only the scale-carrying ones (gradient norm, lmax, barrier, range) do."""
        return (np.asarray(e, float) - self.med_ref) / self.iqr_ref

    def zstd(self, e):
        """The pre-registered mean/sd standardisation, kept for the record."""
        return (np.asarray(e, float) - self.mu_ref) / self.sd_ref


# ============================================================ z <-> (phi, psi)
def unpack(z, n):
    z = np.atleast_2d(np.asarray(z, float))
    return z[:, :n], z[:, n:]


def pack(phi, psi):
    return np.concatenate([np.atleast_2d(phi), np.atleast_2d(psi)], axis=1)


def Ez(ham, z, n, budgeted=True):
    phi, psi = unpack(z, n)
    return ham(phi, psi) if budgeted else ham.raw(phi, psi)


# ============================================================ landscape panel
def fd_grad(ham, z, n, h=FD_H, budgeted=True):
    """Central finite-difference gradient.  Costs 2*2n budget units."""
    d = 2 * n
    Z = np.repeat(np.atleast_2d(z), 2 * d, axis=0)
    for k in range(d):
        Z[2 * k, k] += h
        Z[2 * k + 1, k] -= h
    e = Ez(ham, Z, n, budgeted)
    if len(e) < 2 * d:
        return None
    return (e[0::2] - e[1::2]) / (2 * h)


def fd_hess(ham, z, n, h=FD_H, budgeted=False):
    """Central finite-difference Hessian.  Costs 2d + 2d(d-1) budget units (unbudgeted by
    default; the Hessian is a diagnostic, not an optimiser step)."""
    d = 2 * n
    z = np.atleast_1d(np.asarray(z, float))
    e0 = float(Ez(ham, z[None], n, budgeted)[0])
    Zp = np.repeat(z[None], d, 0); Zm = np.repeat(z[None], d, 0)
    for k in range(d):
        Zp[k, k] += h; Zm[k, k] -= h
    ep = Ez(ham, Zp, n, budgeted); em = Ez(ham, Zm, n, budgeted)
    H = np.zeros((d, d))
    np.fill_diagonal(H, (ep + em - 2 * e0) / h ** 2)
    P, M = [], []
    idx = []
    for a in range(d):
        for b in range(a + 1, d):
            zp = z.copy(); zp[a] += h; zp[b] += h; P.append(zp)
            zm = z.copy(); zm[a] -= h; zm[b] -= h; M.append(zm)
            idx.append((a, b))
    if idx:
        Epp = Ez(ham, np.array(P), n, budgeted)
        Emm = Ez(ham, np.array(M), n, budgeted)
        for k, (a, b) in enumerate(idx):
            v = (Epp[k] + Emm[k] - ep[a] - em[a] - ep[b] - em[b] + 2 * e0) / (2 * h ** 2)
            H[a, b] = H[b, a] = v
    return H


def spectrum(H, scale=1.0):
    """The declared Hessian panel, on the STANDARDISED field.  LOCAL by construction --
    never used to infer global topology (PREREG section 3)."""
    Hs = np.asarray(H, float) / scale
    if not np.isfinite(Hs).all():
        return {k: float("nan") for k in
                ("frac_neg", "cond", "gap", "frac_nearzero", "aniso", "lmax", "lmin_abs",
                 "trace")}
    w = np.linalg.eigvalsh((Hs + Hs.T) / 2.0)
    a = np.abs(w)
    amax = a.max() if a.size else 0.0
    keep = a >= 1e-8 * max(amax, 1e-300)
    return {"frac_neg": float((w < 0).mean()),
            "cond": float(a[keep].max() / a[keep].min()) if keep.any() else float("inf"),
            "gap": float((np.sort(w)[1] - np.sort(w)[0]) / (abs(np.sort(w)[0]) + 1e-12))
            if len(w) > 1 else float("nan"),
            "frac_nearzero": float((a < 1e-3 * max(amax, 1e-300)).mean()),
            "aniso": float(amax / (a.mean() + 1e-300)),
            "lmax": float(w.max()), "lmin_abs": float(a.min()), "trace": float(w.sum())}


def line_scan(ham, z, n, rng, npts=65, budgeted=False):
    """1-D ruggedness along a random unit direction, +/- pi.  Scale-free indices only."""
    d = 2 * n
    u = rng.normal(size=d); u /= np.linalg.norm(u)
    ts = np.linspace(-np.pi, np.pi, npts)
    Z = np.atleast_2d(z) + ts[:, None] * u[None, :]
    e = Ez(ham, Z, n, budgeted)
    if len(e) < npts or not np.isfinite(e).all():
        return None
    es = ham.std(e)
    # local minima per 2pi
    lm = int(((es[1:-1] < es[:-2]) & (es[1:-1] < es[2:])).sum())
    # detrended autocorrelation length, in radians
    x = es - es.mean()
    if x.std() <= 0:
        acl = float("nan")
    else:
        ac = np.correlate(x, x, "full")[len(x) - 1:]
        ac = ac / ac[0]
        below = np.where(ac < 1.0 / np.e)[0]
        acl = float(below[0] * (ts[1] - ts[0])) if len(below) else float(ts[-1] - ts[0])
    rng_e = float(es.max() - es.min())
    tv = float(np.abs(np.diff(es)).sum())
    return {"n_localmin": lm, "acorr_len_rad": acl,
            "tv_over_range": tv / rng_e if rng_e > 0 else float("nan"),
            "range_std": rng_e}


def barrier(ham, za, zb, n, npts=33, budgeted=False):
    ts = np.linspace(0.0, 1.0, npts)
    Z = (1 - ts)[:, None] * np.atleast_2d(za) + ts[:, None] * np.atleast_2d(zb)
    e = Ez(ham, Z, n, budgeted)
    if len(e) < npts or not np.isfinite(e).all():
        return float("nan")
    es = ham.std(e)
    return float(es.max() - max(es[0], es[-1]))


# ============================================================ starts
def basin_starts(tgt, R, rng):
    """R independent starts drawn from the SAME per-residue von Mises basin mixtures every
    arm is initialised from.  The mixtures are NATIVE-FREE (fitted on the K=500 pool)."""
    n = tgt["n"]
    b = (rng.random((R, n)) < tgt["wmarg"][:, 1][None, :]).astype(np.int64)
    phi, psi = QB.draw_from_basins(b, tgt["mu"], tgt["kap"], rng)
    return pack(phi, psi)


def cached_start(pdb):
    """`s19/cache/start_<pdb>.npz` -- the deterministic projected start, identical for every
    arm and lane (BRIEF section 5)."""
    p = os.path.join(ROOT, "s19", "cache", f"start_{pdb}.npz")
    if not os.path.exists(p):
        return None
    z = np.load(p)
    return pack(np.asarray(z["phi"], float), np.asarray(z["psi"], float))[0]


# ============================================================ ORACLE scoring
def rmsd_of(z, tgt):
    """ORACLE, post-hoc only.  Full-chain Ca-RMSD, all residues, proper rotations."""
    from core import project as pj
    phi, psi = unpack(z, tgt["n"])
    CA = np.asarray(pj.build_ca_exact(phi, psi), float)
    return I.kabsch_rmsd_batch(CA, tgt["nat"])


# ============================================================ target bundle
def target(pdb):
    t = QB.target(pdb)
    t["inv2"] = 1.0 / np.asarray(t["sd"], float) ** 2
    sp = QA.Space(pdb, 4)
    t["rep"] = sp.rep
    return t


def hams(tgt, kinds=("LEG", "AMB", "AMBc", "DIST"), budget=10 ** 9, keep=False):
    sp = None
    if "AMB" in kinds or "AMBc" in kinds:
        sp = AmberSP(tgt["seq"], tgt["rep"])
    return {k: Ham(k, tgt, budget=budget, keep=keep, sp=sp) for k in kinds}, sp


# ============================================================ io
def write(name, obj, complete=True):
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    tmp = path + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for _ in range(30):
        try:
            os.replace(tmp, path)
            break
        except PermissionError:
            time.sleep(0.4)
    if complete:
        with open(path.replace(".json", "") + "_COMPLETE", "w") as fh:
            fh.write(time.strftime("%Y-%m-%d %H:%M:%S"))
    return path


def read(name):
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    with open(path) as fh:
        return json.load(fh)


def spearman(a, b):
    return QA.spearman(a, b)


def _ranks(x):
    return QA._rank(np.asarray(x, float))


def partial_spearman(x, y, z):
    """Rank-based partial correlation of x and y given z.

    BINDING REQUIREMENT (coordinator, s20 LEDGER L4): on this instrument a RAW rho against
    Ca-RMSD is largely a TARGET-DIFFICULTY measurement -- every circuit-side distributional
    metric collapsed from |rho| ~ 0.25 to ~0.05 once pool difficulty was partialled out.  So
    every landscape metric here is reported BOTH ways, and the partial is the one that counts.

    Computed on ranks: regress rank(x) and rank(y) each on rank(z) (with an intercept) and
    correlate the residuals.  With one control this equals the usual partial-correlation
    formula; it is an algebraic identity, not a finding.
    """
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if ok.sum() < 5:
        return float("nan")
    rx, ry, rz = _ranks(x[ok]), _ranks(y[ok]), _ranks(z[ok])
    A = np.column_stack([np.ones(ok.sum()), rz])
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    d = np.sqrt((ex ** 2).sum() * (ey ** 2).sum())
    return float((ex * ey).sum() / d) if d > 0 else float("nan")


def boot_mean_ci(x, n_boot=4000, seed=0):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return {"n": int(len(x)), "mean": float(np.mean(x)) if len(x) else float("nan"),
                "ci95": [float("nan"), float("nan")]}
    rng = np.random.default_rng(seed)
    bs = np.array([x[rng.integers(0, len(x), len(x))].mean() for _ in range(n_boot)])
    return {"n": int(len(x)), "mean": float(x.mean()), "median": float(np.median(x)),
            "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]}


def paired_ci(a, b, folds=None, n_boot=4000, seed=0):
    """Paired target-level bootstrap.  Negative = a better.  Median and W/L beside the mean;
    a mean is never reported alone (BRIEF section 8)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    d = a - b
    if len(d) < 3:
        return {"n": int(len(d)), "mean": float(np.mean(d)) if len(d) else float("nan")}
    rng = np.random.default_rng(seed)
    bs = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(n_boot)])
    out = {"n": int(len(d)), "mean_a": float(a.mean()), "mean_b": float(b.mean()),
           "mean": float(d.mean()), "median": float(np.median(d)),
           "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
           "W": int((d < 0).sum()), "L": int((d > 0).sum())}
    out["sig"] = bool(out["ci95"][0] * out["ci95"][1] > 0)
    if folds is not None:
        f = np.asarray(folds)[ok]
        per = {int(k): float(d[f == k].mean()) for k in np.unique(f)}
        out["per_fold"] = per
        s = np.sign(d.mean())
        out["folds_same_sign"] = int(sum(1 for v in per.values() if np.sign(v) == s))
    return out


# ============================================================ the RELAXATION arms
class AmberRelax:
    """`E_amber o Relax_k` -- the AMBER objective the DEPLOYED VQE Hamiltonian actually uses.

    THE CONFOUND THIS EXISTS TO MEASURE.  `amber_hamiltonian.AmberHamiltonian.energy` is not
    an energy function of the configuration: it builds coordinates, applies a
    `restraint_k = 100 kcal/mol/A^2` positional restraint on every heavy atom, runs
    `openmm.LocalEnergyMinimizer.minimize(tolerance, minimization_steps=50)`, then reports the
    UNRESTRAINED energy at the MINIMISED coordinates.  So `H_AMBER = E o Relax`, and every
    landscape metric measured on it is a property of the composition, not of AMBER physics.
    This class makes the relaxation an explicit, ablatable arm (steps = 1 and steps = 50)
    instead of a hidden term.

    `steps = 1` is the smallest honest setting: OpenMM reads `maxIterations = 0` as
    UNBOUNDED, so 0 is not the floor.

    THE COLLAPSE SENTINEL, reported rather than assumed.  `core.amber._run` -- the entry point
    `AmberSP` matches bit-exactly -- **never calls `_is_collapsed`** (verified in source: the
    only call site is `AmberHamiltonian._evaluate`, line 1168).  So the bare single-point arm
    has NO infinite-valued region by construction, not merely an unfired guard.  For the
    relaxed arms the deployed predicate (`collapse_mode="geometry"`: any heavy-atom pair closer
    than `collapse_min_contact = 1.05 A`) is evaluated explicitly and COUNTED, so `n_collapsed`
    is a measured number on every arm.
    """

    def __init__(self, seq, rep, k_restraint=100.0, steps=50, tolerance=2.0,
                 collapse_min_contact=1.05):
        from core import amber as am
        self.am = am
        self.seq = seq
        self.rep = rep
        self.k = float(k_restraint)
        self.steps = int(steps)
        self.tol = float(tolerance)
        self.cmin = float(collapse_min_contact)
        self.H = am.builder_for(seq, rep, "CPU", 1)
        self.n_calls = 0
        self.n_collapsed = 0

    def _e(self, coords):
        from openmm import unit
        import openmm
        H = self.H
        pos = H._assemble(H._heavy_positions(coords, chi1=None))
        ctx = H.context
        ctx.setPositions(pos * unit.nanometer)
        for j, idx in enumerate(H._restraint_idx):
            H._rest_force.setParticleParameters(j, idx, pos[idx].tolist())
        H._rest_force.updateParametersInContext(ctx)
        ctx.setParameter("k_rest", self.k * self.am._K_SCALE)
        openmm.LocalEnergyMinimizer.minimize(ctx, self.tol, self.steps)
        ctx.setParameter("k_rest", 0.0)
        st = ctx.getState(getEnergy=True, getPositions=True)
        e = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        p = np.asarray(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer), float)
        hv = p[H._heavy_order] * 10.0
        d = np.linalg.norm(hv[:, None, :] - hv[None, :, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        if d.min() < self.cmin:
            self.n_collapsed += 1
        return float(e)

    def batch(self, phi, psi):
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        c = geo.build_backbone_batch(phi, psi)
        out = np.empty(len(phi))
        for b in range(len(phi)):
            out[b] = self._e({k: v[b] for k, v in c.items()})
        self.n_calls += len(phi)
        return out

    def verify(self, phi, psi, m=3):
        """BIT-EXACTNESS GATE against `core.amber.refine_coords` at the same settings."""
        phi = np.atleast_2d(np.asarray(phi, float))[:m]
        psi = np.atleast_2d(np.asarray(psi, float))[:m]
        c = geo.build_backbone_batch(phi, psi)
        worst, k = 0.0, 0
        for b in range(len(phi)):
            cc = {kk: v[b] for kk, v in c.items()}
            ref = float(self.am.refine_coords(self.seq, self.rep, cc, k_restraint=self.k,
                                              steps=self.steps, tolerance=self.tol,
                                              threads=1, memo=False)["energy"])
            got = self._e(cc)
            worst = max(worst, abs(ref - got) / max(1.0, abs(ref)))
            k += 1
        return worst, k
