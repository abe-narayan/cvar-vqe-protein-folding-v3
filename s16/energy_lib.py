"""SPRINT 16, ENERGY -- shared machinery for the two mandated energy pillars.

Genuine Legacy/MJ (`core.energy`, eleven components, `DEFAULT_WEIGHTS`) and genuine
all-atom AMBER ff14SB/GBn2 (`core.amber`).  Nothing here substitutes a learned surrogate
for either; the Legacy total is always the weighted sum of the eleven measured components
and the AMBER numbers always come from OpenMM.

FOUR THINGS THIS MODULE FIXES, all of which bind every statistic downstream.

1. `gate()` -- the pre-declared AMBER convergence gate.  Implemented in
   `core.amber.convergence_flags` (declared there, in the code, before it was applied);
   this is the consumer side.  Four of 126 targets end minimisation above 1000 kcal/mol.

2. `frame_null()` -- the standing rigid-invariance regression check.  ff14SB, GBn2, the
   positional restraint and every RMSD in this project are rigid-invariant, so relaxing
   the SAME structure in a rotated lab frame must return EXACTLY zero.  It does not.

3. `condition()` -- monotone conditioning.  37.6% of cached AMBER single points exceed
   1e6 kcal/mol, 3.2% exceed 1e12 and the maximum is 2.6e20, so the 99th percentile is
   itself 1e12-1e15 and winsorising there conditions nothing.  Only a MONOTONE map (rank
   or signed-log) may be applied before a spectral or moment-based claim, because a
   monotone map preserves every ordering statement and destroys the moment artefact.

4. `binding_mask()` -- the AMBER data rule.  `amber_kind == 0 AND amber_idx != snap_index`
   for every tail / in-band / argmin / top-k statistic, with per-target n printed.
   `s14/results/obj_floor.json` is RETRACTED and is never opened here.

ORACLE labelling: every function whose name contains `ORACLE`, and every field named
`rmsd`/`rr`/`nat_*`, is a native-derived evaluation label.  None of them enters a
selection, a threshold or a fitted parameter except where the docstring says
"ORACLE DIAGNOSTIC" or "leave-fold-out label".
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import glob
import json
import math
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s15.seed import stable_rng            # noqa: E402

RESULTS = os.path.join(ROOT, "s16", "results")
CACHE = os.path.join(ROOT, "s16", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

ATOMS = ("N", "CA", "C", "O", "CB")
IDEAL = dict(N_CA=1.458, CA_C=1.525, C_N=1.329, C_O=1.231,
             ang_N_CA_C=111.0, ang_CA_C_N=116.2, ang_C_N_CA=121.7, ca_ca=3.804)


# ==========================================================================
# 1.  THE CONVERGENCE GATE  (consumer side; the rule lives in core.amber)
# ==========================================================================
from core.amber import CONVERGE_MAX_KCAL, convergence_flags   # noqa: E402,F401


def gate(energies, max_kcal=CONVERGE_MAX_KCAL):
    """Boolean keep-mask over an array of FINAL AMBER energies (restraint off)."""
    e = np.asarray(energies, float)
    return np.isfinite(e) & (e <= float(max_kcal))


def gated_paired(a, b, energies, folds=None, names=None, seed=0,
                 max_kcal=CONVERGE_MAX_KCAL):
    """`I.paired(a, b)` reported BOTH ways: ungated and convergence-gated.

    Never returns one without the other.  `n_excluded` and the excluded names are always
    printed, because a gate that is not counted is a filter.
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    keep = gate(energies, max_kcal)
    out = {"ungated": I.paired(a, b, folds=folds, names=names, seed=seed),
           "n_excluded": int((~keep).sum()),
           "excluded": [str(names[i]) for i in np.where(~keep)[0]] if names is not None
                       else np.where(~keep)[0].tolist(),
           "max_kcal": float(max_kcal)}
    if keep.sum() >= 3:
        out["gated"] = I.paired(a[keep], b[keep],
                                folds=None if folds is None else np.asarray(folds)[keep],
                                names=None if names is None else
                                      [names[i] for i in np.where(keep)[0]], seed=seed)
    return out


# ==========================================================================
# 2.  THE STANDING FRAME-INVARIANCE NULL
# ==========================================================================
def random_rigid(rng):
    """Proper rotation (det = +1) and translation.  A reflection is NOT a symmetry of
    ff14SB, so O(3) is deliberately reduced to SO(3)."""
    Q, R = np.linalg.qr(rng.normal(size=(3, 3)))
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1.0
    return Q, rng.normal(scale=10.0, size=3)


#: The regression check's pre-declared PASS band, on the CONVERGED subset only.
#: |mean dRMSD| <= 0.005 A and max |dRMSD| <= 0.05 A.  Both numbers come from the
#: Sprint 15 gated measurement (-0.0005 [-0.0046, +0.0035]) and are stated here before
#: any Sprint 16 draw was taken.
FRAME_TOL_MEAN = 0.005
FRAME_TOL_MAX = 0.05


def frame_null(targets=None, draw=1, k=30.0, max_kcal=CONVERGE_MAX_KCAL, verbose=True):
    """STANDING NULL.  Relax the same averaged backbone in two lab frames.

    Zero by construction: ff14SB + GBn2 + the positional restraint + Kabsch RMSD are all
    rigid-invariant.  Anything this returns is the AMBER minimiser's own numerical floor.
    Returns per-target rows and the PASS/FAIL verdict on the CONVERGED subset.
    """
    from core import amber as am
    import torsion_lib2 as tl2
    from s14.avgspace import top75_windows
    from s15.phys_repl import averaged_backbone_from

    tg = targets if targets is not None else I.targets()
    rows = []
    for t in tg:
        pdb, seq = t["pdb"], t["seq"]
        W, PHI, PSI, u = top75_windows(pdb)
        avg, C_ca, _ = averaged_backbone_from(W, PHI, PSI)
        rng = stable_rng("s16energy", "frame", int(draw), pdb)
        Q, tv = random_rigid(rng)
        rot = {a: v @ Q.T + tv for a, v in avg.items()}
        tab = tl2.library_for(seq, 4, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        out = {}
        for tag, coords in (("A", avg), ("B", rot)):
            r = am.refine_coords(seq, rep, coords, k_restraint=float(k), steps=0,
                                 tolerance=1.0, threads=1, memo=False)
            out[tag] = dict(rmsd=float(I.ca_rmsd(np.asarray(r["ca"], float), u["nat_ca"])),
                            energy=float(r["energy"]), converged=bool(r["converged"]))
        rows.append(dict(pdb=pdb, fold=int(t["fold"]),
                         rmsd_A=out["A"]["rmsd"], rmsd_B=out["B"]["rmsd"],
                         d_rmsd=out["B"]["rmsd"] - out["A"]["rmsd"],
                         e_A=out["A"]["energy"], e_B=out["B"]["energy"],
                         converged=bool(out["A"]["converged"] and out["B"]["converged"])))
        if verbose:
            r0 = rows[-1]
            print(f"  frame {pdb} dRMSD {r0['d_rmsd']:+.5f} conv {r0['converged']}",
                  flush=True)
    d = np.array([r["d_rmsd"] for r in rows])
    ok = np.array([r["converged"] for r in rows])
    verdict = {}
    for tag, m in (("ungated", np.ones(len(d), bool)), ("gated", ok)):
        if m.sum() == 0:
            continue
        verdict[tag] = dict(n=int(m.sum()), mean=float(d[m].mean()),
                            sd=float(d[m].std(ddof=1)) if m.sum() > 1 else 0.0,
                            max_abs=float(np.abs(d[m]).max()),
                            n_nonzero=int((np.abs(d[m]) > 1e-9).sum()))
    v = verdict.get("gated", verdict.get("ungated"))
    verdict["PASS"] = bool(abs(v["mean"]) <= FRAME_TOL_MEAN and
                           v["max_abs"] <= FRAME_TOL_MAX)
    verdict["tol"] = dict(mean=FRAME_TOL_MEAN, max_abs=FRAME_TOL_MAX)
    verdict["n_excluded_by_gate"] = int((~ok).sum())
    verdict["excluded"] = [r["pdb"] for r in rows if not r["converged"]]
    return rows, verdict


# ==========================================================================
# 3.  MONOTONE CONDITIONING
# ==========================================================================
def condition(E, how="signed_log", ref=None):
    """Monotone conditioning of an energy array.  Order-preserving BY CONSTRUCTION.

    `signed_log`  x -> sign(x - ref) * log1p(|x - ref|), ref = median by default.
                  Strictly increasing, so every "A ranks below B" statement survives.
    `rank`        x -> its within-array rank scaled to [0, 1].  Also strictly monotone
                  (ties averaged), and completely tail-free.
    `winsor99`    the WRONG one, kept so the report can quote it: clipping at the 99th
                  percentile of a distribution whose 99th percentile is 1e12-1e15.

    Winsorising is NOT monotone-preserving above the clip (it is constant there), which
    is exactly why it fails: it maps the entire pathological tail onto one value while
    leaving that value 1e12 away from the bulk.
    """
    x = np.asarray(E, float)
    fin = np.isfinite(x)
    if how == "rank":
        out = np.full(x.shape, np.nan)
        r = _avg_rank(x[fin])
        out[fin] = r / max(len(r) - 1, 1)
        return out
    if how == "winsor99":
        hi = np.percentile(x[fin], 99.0)
        return np.clip(x, None, hi)
    r = float(np.median(x[fin])) if ref is None else float(ref)
    d = x - r
    return np.sign(d) * np.log1p(np.abs(d))


def _avg_rank(x):
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x), float)
    r[order] = np.arange(len(x), dtype=float)
    xs = x[order]
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def tail_share(E, top=10):
    """Fraction of the total (centred) variance carried by the `top` extreme values.

    The diagnostic that made the Pauli-spectrum artefact visible: on RAW AMBER the top
    ten configurations carry 99.6% of the Walsh variance.
    """
    x = np.asarray(E, float)
    x = x[np.isfinite(x)]
    c = x - x.mean()
    v = (c ** 2)
    idx = np.argsort(-np.abs(c))[:top]
    return float(v[idx].sum() / v.sum()) if v.sum() > 0 else float("nan")


# ==========================================================================
# 4.  THE BINDING AMBER DATA RULE
# ==========================================================================
def binding_mask(z):
    """`amber_kind == 0 AND amber_idx != snap_index`.  Returns (positions, report).

    `positions` indexes the AMBER subsample arrays (`amber_total`, `amb_*`).  Use
    `z['amber_idx'][positions]` to reach the full enumeration.
    """
    kind = np.asarray(z["amber_kind"], int)
    idx = np.asarray(z["amber_idx"], int)
    snap = int(z["snap_index"])
    m = (kind == 0) & (idx != snap)
    rep = {"n_amber": int(len(idx)), "n_kind0": int((kind == 0).sum()),
           "n_binding": int(m.sum()), "snap_in_kind0": bool((idx[kind == 0] == snap).any())}
    return np.where(m)[0], rep


ENUM_FILES = (sorted(glob.glob(os.path.join(ROOT, "s13", "results", "qarch_enum_*.npz")))
              + sorted(glob.glob(os.path.join(ROOT, "s14", "cache", "obj_enum_*.npz"))))
LEG_TERMS = ("steric", "contact", "hbond_local", "hbond_longrange", "coop_helix",
             "coop_sheet", "solvation", "electrostatic", "aromatic", "torsion",
             "compactness")


def load_enum(path, terms=True):
    z = np.load(path)
    out = {"pdb": str(z["pdb"]), "n": int(z["n"]), "k": int(z["k"]),
           "seq": str(z["seq"]), "fold": int(z["fold"]),
           "rmsd": np.asarray(z["rmsd"], float),
           "legacy": np.asarray(z["legacy"], float),
           "prior": np.asarray(z["prior"], float),
           "snap_index": int(z["snap_index"]),
           "amber_idx": np.asarray(z["amber_idx"], int),
           "amber_kind": np.asarray(z["amber_kind"], int),
           "amber_total": np.asarray(z["amber_total"], float),
           "PHI": np.asarray(z["PHI"], float), "PSI": np.asarray(z["PSI"], float),
           "snap_states": np.asarray(z["snap_states"], int)}
    if terms:
        out["comp"] = {t: np.asarray(z["leg_" + t], float) for t in LEG_TERMS}
        out["amb"] = {t: np.asarray(z["amb_" + t], float)
                      for t in ("bond", "angle", "torsion", "nonbonded", "solvation")}
    out["_z"] = z
    return out


# ==========================================================================
# 5.  LEGACY -- the genuine eleven-component potential
# ==========================================================================
def legacy_weight_vector(w=None):
    from core import energy as et
    w = et.DEFAULT_WEIGHTS if w is None else w
    return np.array([float(w[t]) for t in LEG_TERMS])


def legacy_total_from(comp, w=None):
    """The eleven-term weighted sum.  Genuine Legacy: no surrogate, no learned stand-in
    unless `w` is an explicitly-fitted vector, which every caller must label."""
    W = legacy_weight_vector() if w is None else np.asarray(w, float)
    M = np.column_stack([np.asarray(comp[t], float) for t in LEG_TERMS])
    return M @ W


def legacy_components_of_windows(seq, PHI, PSI, chunk=2048):
    """The eleven Legacy components of a batch of (B, n) torsion sets.

    Ideal-geometry rebuild, exactly as `s13.qarch_lib.legacy_components` does it, so the
    numbers are comparable with the enumerated caches.
    """
    from core import energy as et
    from core import geometry as geo
    PHI = np.atleast_2d(np.asarray(PHI, float))
    PSI = np.atleast_2d(np.asarray(PSI, float))
    B = len(PHI)
    out = {t: np.empty(B) for t in LEG_TERMS}
    for a in range(0, B, chunk):
        ph, ps = PHI[a:a + chunk], PSI[a:a + chunk]
        c = geo.build_backbone_batch(ph, ps)
        comp = et.components_batch(seq, c, ph, ps)
        for t in LEG_TERMS:
            out[t][a:a + chunk] = comp[t]
    return out


# ==========================================================================
# 6.  THE STEREOCHEMISTRY PANEL  (the repair benchmark's measuring instrument)
# ==========================================================================
def panel(c, seq):
    """Full stereochemistry panel of one backbone dict (angstrom, keys N/CA/C[/O/CB]).

    Every axis the repair benchmark scores, on ONE structure:
      rama_favoured / rama_allowed / rama_outlier   (three-way, not the binary rama_ok)
      n_clash_2A, n_clash_2p6A, min_heavy           (steric)
      bond_strain, angle_strain                     (rms relative deviation from ideal)
      cis_frac                                      (|omega| < 90 deg, non-PRO)
      chirality_L_frac, chirality_mean              (CB improper)
    """
    from s14.ener_geom import _ang, _dih
    N, CA, Cc = (np.asarray(c[k], float) for k in ("N", "CA", "C"))
    O = np.asarray(c["O"], float) if "O" in c else None
    CB = np.asarray(c["CB"], float) if "CB" in c else None
    g = dict(N_CA=np.linalg.norm(CA - N, axis=1),
             CA_C=np.linalg.norm(Cc - CA, axis=1),
             C_N=np.linalg.norm(N[1:] - Cc[:-1], axis=1),
             ca_ca=np.linalg.norm(CA[1:] - CA[:-1], axis=1),
             ang_N_CA_C=_ang(N, CA, Cc),
             ang_CA_C_N=_ang(CA[:-1], Cc[:-1], N[1:]),
             ang_C_N_CA=_ang(Cc[:-1], N[1:], CA[1:]))
    if O is not None:
        g["C_O"] = np.linalg.norm(O - Cc, axis=1)
    out = {k: float(np.mean(v)) for k, v in g.items()}
    out.update({k + "_min": float(np.min(v)) for k, v in g.items()})

    # -- bond and angle strain, separately (the published scalar fuses them)
    bkeys = [k for k in ("N_CA", "CA_C", "C_N", "C_O", "ca_ca") if k in g]
    akeys = [k for k in ("ang_N_CA_C", "ang_CA_C_N", "ang_C_N_CA") if k in g]
    out["bond_strain"] = float(np.sqrt(np.mean(np.concatenate(
        [((g[k] - IDEAL[k]) / IDEAL[k]) ** 2 for k in bkeys]))))
    out["angle_strain"] = float(np.sqrt(np.mean(np.concatenate(
        [((g[k] - IDEAL[k]) / IDEAL[k]) ** 2 for k in akeys]))))
    dev = [(out[k] - IDEAL[k]) / IDEAL[k] for k in IDEAL if k in out]
    out["geom_rms_rel_dev"] = float(np.sqrt(np.mean(np.square(dev))))
    out["bond_contraction"] = float(np.mean([
        (out[k] - IDEAL[k]) / IDEAL[k] for k in ("N_CA", "CA_C", "C_N", "ca_ca")]))

    # -- omega / cis peptide
    om = np.degrees(_dih(CA[:-1], Cc[:-1], N[1:], CA[1:]))
    out["omega"] = float(np.mean(np.abs(om)))
    non_pro = np.array([i for i in range(len(seq) - 1) if seq[i + 1] != "P"], int)
    sel = om[non_pro] if len(non_pro) else om
    out["cis_frac"] = float((np.abs(sel) < 90.0).mean()) if len(sel) else 0.0
    out["omega_dev"] = float(np.mean(np.abs(np.abs(om) - 180.0)))

    # -- chirality
    if CB is not None:
        ok = np.array([i for i, a in enumerate(seq) if a != "G"], int)
        if len(ok):
            imp = np.degrees(_dih(N[ok], CA[ok], Cc[ok], CB[ok]))
            out["chirality_mean"] = float(imp.mean())
            out["chirality_L_frac"] = float((imp < 0).mean())
            out["n_D_centres"] = int((imp >= 0).sum())

    # -- Ramachandran, three-way
    ph, ps = torsions_of(c)
    out.update(rama3(ph, ps, seq))

    # -- steric
    arrs, res = [], []
    for a in ATOMS:
        if a in c:
            v = np.asarray(c[a], float); arrs.append(v); res.append(np.arange(len(v)))
    Hh = np.vstack(arrs); Ri = np.concatenate(res)
    Dm = np.linalg.norm(Hh[:, None, :] - Hh[None, :, :], axis=-1)
    Dm = np.where(np.abs(Ri[:, None] - Ri[None, :]) >= 2, Dm, 9e9)
    out["min_heavy"] = float(Dm.min())
    out["n_clash_2A"] = int((Dm < 2.0).sum() // 2)
    out["n_clash_2p6A"] = int((Dm < 2.6).sum() // 2)
    return out


def torsions_of(c):
    from s14.ener_geom import _dih
    N, CA, Cc = (np.asarray(c[k], float) for k in ("N", "CA", "C"))
    n = len(CA)
    phi = np.zeros(n); psi = np.zeros(n)
    phi[1:] = _dih(Cc[:-1], N[1:], CA[1:], Cc[1:])
    psi[:-1] = _dih(N[:-1], CA[:-1], Cc[:-1], N[1:])
    return phi, psi


def wrap_deg(a):
    """Radians -> degrees in (-180, 180].

    NECESSARY, and its absence is a defect in the incumbent.  `core.project.lam_path`
    optimises torsions with an unconstrained L-BFGS, so it returns UNWRAPPED angles --
    +1022.72 deg where it means -57.28, -1115.93 where it means -35.93.  Every fixed
    Ramachandran window test (`s14.ener_avgrefine.rama_ok` and its callers) compares raw
    degrees against those windows and therefore scores an unwrapped favoured residue as
    an outlier.  See s16/energy_FINDINGS.md 3.4.
    """
    return (np.degrees(np.asarray(a, float)) + 180.0) % 360.0 - 180.0


def rama_ok(phi, psi):
    """The project's incumbent BINARY favoured test, with the torsions WRAPPED."""
    p = wrap_deg(phi)[1:-1]; s = wrap_deg(psi)[1:-1]
    a = (p > -160) & (p < -20) & (s > -120) & (s < 50)
    b = (p > -180) & (p < -40) & (s > 90) & (s < 180)
    l = (p > 20) & (p < 100) & (s > -20) & (s < 90)
    return float((a | b | l).mean()) if len(p) else float("nan")


def rama3(phi, psi, seq=None):
    """Three-way favoured / allowed / outlier over interior residues.

    `favoured` is the incumbent binary test (unchanged, so it stays comparable).
    `allowed` widens each of the three regions by 20 degrees on every side; anything
    outside both is an OUTLIER.  Fixed constants, declared here, no tuning.
    """
    p = wrap_deg(phi)[1:-1]; s = wrap_deg(psi)[1:-1]
    if len(p) == 0:
        return {"rama_favoured": float("nan"), "rama_allowed": float("nan"),
                "rama_outlier": float("nan"), "n_rama": 0}

    def _reg(pad):
        a = (p > -160 - pad) & (p < -20 + pad) & (s > -120 - pad) & (s < 50 + pad)
        b = (p > -180 - pad) & (p < -40 + pad) & (s > 90 - pad) & (s < 180 + pad)
        l = (p > 20 - pad) & (p < 100 + pad) & (s > -20 - pad) & (s < 90 + pad)
        return a | b | l

    fav = _reg(0.0)
    wide = _reg(20.0)
    return {"rama_favoured": float(fav.mean()),
            "rama_allowed": float((wide & ~fav).mean()),
            "rama_outlier": float((~wide).mean()),
            "n_rama": int(len(p))}


# ==========================================================================
# 7.  small helpers
# ==========================================================================
def spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3:
        return float("nan")
    ra = _avg_rank(a) - _avg_rank(a).mean()
    rb = _avg_rank(b) - _avg_rank(b).mean()
    d = math.sqrt(float((ra ** 2).sum() * (rb ** 2).sum()))
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


def auroc(score, label):
    """P(score of a positive > score of a negative), ties at 0.5.  `label` boolean."""
    s = np.asarray(score, float); y = np.asarray(label, bool)
    ok = np.isfinite(s)
    s, y = s[ok], y[ok]
    if y.all() or (~y).all() or len(s) < 3:
        return float("nan")
    r = _avg_rank(s)
    n1 = int(y.sum()); n0 = int((~y).sum())
    return float((r[y].sum() - n1 * (n1 - 1) / 2.0) / (n1 * n0))


def argmin_tied(score, value):
    """Mean `value` over the FULL tied argmin set.

    `np.argmin` on a tied signal reads the input's sort order, which on this project's
    caches is the ORACLE sort order.  Recorded trap; averaged over instead.
    """
    s = np.asarray(score, float); v = np.asarray(value, float)
    m = np.nanmin(s)
    tie = np.where(s <= m + 1e-12 * max(1.0, abs(m)))[0]
    return float(v[tie].mean()), int(len(tie))


def write(name, obj, n_expected=None):
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    if n_expected is not None and isinstance(obj, dict):
        rows = obj.get("per_target") or obj.get("rows") or []
        obj = dict(obj, n_rows=len(rows), n_expected=int(n_expected))
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return path


def read(name):
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def cpu_hold(limit=90.0, tag="", max_wait=120.0):
    """Yield to a contended box, but never deadlock on it.

    The box is shared with sibling agents.  An unbounded hold at a 96% threshold blocked
    this workstream indefinitely while three other processes sat at 100%, so the hold is
    now BOUNDED: yield for up to `max_wait` seconds, then proceed.  These jobs run
    OpenMM at `Threads=1`, ~1 core of a 6.43-core-equivalent box, so proceeding after the
    yield does not itself breach the ceiling.
    """
    import subprocess, time
    cmd = ("$c=Get-CimInstance Win32_Processor|Measure-Object -Property LoadPercentage "
           "-Average; [int]$c.Average")
    t0 = time.time()
    while True:
        try:
            v = float(subprocess.run(["powershell", "-NoProfile", "-Command", cmd],
                                     capture_output=True, text=True, timeout=60
                                     ).stdout.strip())
        except Exception:
            return -1.0
        if v <= limit or time.time() - t0 >= max_wait:
            return v
        print(f"  [{tag}] CPU {v:.0f}% > {limit}; holding", flush=True)
        time.sleep(20)


def mem_hold(min_gb=1.6, tag="", max_wait=600.0):
    """Block until physical memory frees up.  THIS is the gate that actually bites:
    `core.amber.memory_guard` raises MemoryError above a 92% ceiling and killed a shard
    of the Sprint 16 ablation mid-run."""
    import time
    t0 = time.time()
    while True:
        g = I.free_gb()
        if g >= min_gb or time.time() - t0 >= max_wait:
            return g
        print(f"  [{tag}] free {g:.2f} GB < {min_gb}; waiting", flush=True)
        time.sleep(15)
