"""s17/phys_lib.py -- Ca-preserving AMBER repair posed as CONSTRAINED OPTIMISATION.

    minimise  E_ff14SB/GBn2(x)   subject to   ||CA_i(x) - CA_i(x0)|| <= eps  for every i

Two exact implementations of the constraint, both built on `core.amber`'s genuine
ff14SB/GBn2 Hamiltonian (nothing here is a surrogate, and no arithmetic in `core/amber.py`
is touched):

  HARD   eps = 0.  The Ca particle masses are set to zero, which OpenMM's
         `LocalEnergyMinimizer` treats as frozen.  The Ca coordinates come out bit-identical
         to the input, so the Ca-RMSD cost is EXACTLY zero -- by construction, not by
         measurement.  That is declared in PREREG_phys.md and is why the eps = 0 rung is
         judged on validity alone.

  FLAT   eps > 0.  A per-atom flat-bottom external restraint
         0.5 * k_flat * max(0, |x - x0| - eps)^2 on the Ca set, with k_flat large.  This is
         the exact penalty form of the per-atom eps-ball; the residual violation is bounded
         by F/k_flat and is measured and reported, never assumed.

N, C, O, CB, sidechains and every hydrogen are FREE at every rung.  The incumbent's harmonic
N/CA/C restraint (k = 30) is carried as a comparison arm, not as a rung of the ladder.

One OpenMM Context serves the whole free-Ca ladder: the restraint forces carry global
parameters (kh_ca, kf_ca, eps_ca, kh_bb) so a rung is a `setParameter` call, not a rebuild.
The hard rung needs a second Context because particle mass is a System property.

CONVERGENCE GATE.  `core.amber.convergence_flags` -- final potential energy with the
restraint switched off, finite and <= 1000 kcal/mol.  Declared before use, reported with its
exclusion count and the excluded ids, never silently applied.
"""
from __future__ import annotations

import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                     # noqa: E402
from s15 import seed as SD                          # noqa: E402
from s16 import energy_lib as EL                    # noqa: E402

ATOMS = ("N", "CA", "C", "O", "CB")

#: k for the flat-bottom wall, kcal/mol/A^2, on the CROSS-CHECK rungs only.
#:
#: WHY THE MAIN LADDER IS HARMONIC AND NOT FLAT-BOTTOM.  The constrained problem
#:      min E(x)  s.t.  ||CA(x) - CA(x0)|| <= eps
#: has the Lagrangian family  min E(x) + (k/2)||CA(x) - CA(x0)||^2, and sweeping k traces
#: exactly the same Pareto frontier as sweeping eps -- k is the dual variable of eps.  The
#: harmonic family is SMOOTH; the flat-bottom family has a kink at |x-x0| = eps which
#: L-BFGS handles badly.  Measured on 1A13 at eps = 0.10 A: the wall needs k_f = 1e6 to hold
#: the constraint to 0.004 A and then costs 78 s per target per rung, against 7.6 s for the
#: harmonic arm; at an affordable k_f = 1e3 it OVERSHOOTS the constraint by 0.267 A, i.e. it
#: is not enforcing the constraint at all.  So the ladder is swept in k, the frontier is
#: plotted against the REALISED Ca displacement (the quantity the constraint bounds, measured
#: not assumed), and two flat-bottom rungs are carried as a cross-check that the two
#: parameterisations trace the same frontier.
K_FLAT = 100000.0

#: The ladder.  `name`, kh_ca, kf_ca, eps_ca, kh_bb, hard.
#: kh_* are harmonic kcal/mol/A^2 (Ca-only unless kh_bb > 0); kf_ca/eps_ca the flat wall.
LADDER: Tuple[Tuple[str, float, float, float, float, bool], ...] = (
    ("cafix",     0.0,    0.0,  0.00,  0.0, True),     # eps = 0 EXACTLY, hard (mass 0)
    ("ca1000", 1000.0,    0.0,  0.00,  0.0, False),    # extremely tight
    ("ca300",   300.0,    0.0,  0.00,  0.0, False),
    ("ca100",   100.0,    0.0,  0.00,  0.0, False),    # tight
    ("ca30",     30.0,    0.0,  0.00,  0.0, False),    # moderate (incumbent k, Ca only)
    ("ca10",     10.0,    0.0,  0.00,  0.0, False),
    ("ca3",       3.0,    0.0,  0.00,  0.0, False),    # loose
    ("free",      0.0,    0.0,  0.00,  0.0, False),    # unrestrained
    ("k30",      30.0,    0.0,  0.00, 30.0, False),    # THE INCUMBENT: harmonic N/CA/C
    #: flat-bottom cross-check rungs (expensive; subsample only)
    ("eps010",    0.0, K_FLAT,  0.10,  0.0, False),
    ("eps025",    0.0, K_FLAT,  0.25,  0.0, False),
)
#: The rungs run at full n = 126.
LADDER_FULL = ("cafix", "ca1000", "ca300", "ca100", "ca30", "ca10", "ca3", "free", "k30")
#: Cross-check rungs, pre-registered subsample only.
LADDER_XCHECK = ("eps010", "eps025")
LADDER_BY_NAME = {r[0]: r for r in LADDER}


# --------------------------------------------------------------------------- the box
class ConstrainedBox:
    """One target's genuine AMBER Hamiltonian plus the two constrained-minimisation contexts.

    Built lazily; `close()` frees both contexts and the systems behind them.
    """

    def __init__(self, seq: str, rep, threads: int = 1):
        from core import amber as am
        self.am = am
        self.seq = seq
        self.rep = rep
        self.H = am.builder_for(seq, rep, "CPU", threads)
        self._soft = None       # (context, ca_force, bb_force, system)
        self._hard = None       # (context, system)
        self.threads = threads

    # -- system construction ------------------------------------------------
    def _clone_system(self):
        """A fresh System with the same forces as the builder's, minus its restraint."""
        import openmm
        from openmm import app, unit
        H = self.H
        sysc = H.forcefield.createSystem(
            H.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None,
            rigidWater=False,
            removeCMMotion=False,
            implicitSolventKappa=0.0 / unit.nanometer)
        from core.amber import _GROUP_OF
        for f in sysc.getForces():
            f.setForceGroup(_GROUP_OF.get(f.__class__.__name__, 5))
        return sysc

    def _ca_idx(self) -> List[int]:
        H = self.H
        return [H._heavy_index[(i, "CA")] for i in range(len(self.seq))]

    def _bb_idx(self) -> List[int]:
        H = self.H
        return [H._heavy_index[(i, nm)] for i in range(len(self.seq)) for nm in ("N", "C")]

    def _make_soft(self):
        import openmm
        from openmm import unit
        sysc = self._clone_system()
        # Ca: harmonic + flat-bottom wall, both switchable by global parameter.
        ca = openmm.CustomExternalForce(
            "0.5*kh_ca*d2 + 0.5*kf_ca*max(0, sqrt(d2)-eps_ca)^2; "
            "d2=(x-x0)^2+(y-y0)^2+(z-z0)^2")
        for p in ("kh_ca", "kf_ca", "eps_ca"):
            ca.addGlobalParameter(p, 0.0)
        for nm in ("x0", "y0", "z0"):
            ca.addPerParticleParameter(nm)
        for idx in self._ca_idx():
            ca.addParticle(idx, [0.0, 0.0, 0.0])
        ca.setForceGroup(7)
        sysc.addForce(ca)
        # N and C: harmonic only (the incumbent arm).
        bb = openmm.CustomExternalForce(
            "0.5*kh_bb*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        bb.addGlobalParameter("kh_bb", 0.0)
        for nm in ("x0", "y0", "z0"):
            bb.addPerParticleParameter(nm)
        for idx in self._bb_idx():
            bb.addParticle(idx, [0.0, 0.0, 0.0])
        bb.setForceGroup(8)
        sysc.addForce(bb)
        ctx = openmm.Context(sysc, openmm.VerletIntegrator(0.001 * unit.picoseconds),
                             self.H.platform, self.H.platform_properties)
        self._soft = (ctx, ca, bb, sysc)

    def _make_hard(self):
        import openmm
        from openmm import unit
        sysc = self._clone_system()
        for idx in self._ca_idx():
            sysc.setParticleMass(idx, 0.0)
        ctx = openmm.Context(sysc, openmm.VerletIntegrator(0.001 * unit.picoseconds),
                             self.H.platform, self.H.platform_properties)
        self._hard = (ctx, sysc)

    def close(self):
        for slot in ("_soft", "_hard"):
            v = getattr(self, slot)
            if v is not None:
                try:
                    del v
                except Exception:
                    pass
                setattr(self, slot, None)

    # -- the run ------------------------------------------------------------
    def relax(self, coords: Dict[str, np.ndarray], rung: str,
              steps: int = 0, tolerance: float = 1.0) -> Dict[str, object]:
        """One constrained minimisation of `coords` at ladder rung `rung`.

        `coords` is a backbone dict in angstrom (N, CA, C, O, CB); sidechains are built by
        `core.amber.build_full_structure` exactly as `refine_coords` does.
        """
        from openmm import unit
        import openmm
        name, kh_ca, kf_ca, eps_ca, kh_bb, hard = LADDER_BY_NAME[rung]
        H = self.H
        t0 = time.time()
        heavy_nm = H._heavy_positions(coords, chi1=None)
        pos = H._assemble(heavy_nm)                       # nm, hydrogens in frozen frames

        if hard:
            if self._hard is None:
                self._make_hard()
            ctx = self._hard[0]
            ctx.setPositions(pos * unit.nanometer)
            e0 = _energy(ctx)
            openmm.LocalEnergyMinimizer.minimize(ctx, float(tolerance), int(steps))
            energy = _energy(ctx)
        else:
            if self._soft is None:
                self._make_soft()
            ctx, caf, bbf, _ = self._soft
            ctx.setPositions(pos * unit.nanometer)
            for j, idx in enumerate(self._ca_idx()):
                caf.setParticleParameters(j, idx, pos[idx].tolist())
            caf.updateParametersInContext(ctx)
            for j, idx in enumerate(self._bb_idx()):
                bbf.setParticleParameters(j, idx, pos[idx].tolist())
            bbf.updateParametersInContext(ctx)
            _zero(ctx)
            e0 = _energy_raw(ctx)
            ctx.setParameter("kh_ca", float(kh_ca) * _KS)
            ctx.setParameter("kf_ca", float(kf_ca) * _KS)
            ctx.setParameter("eps_ca", float(eps_ca) * 0.1)      # A -> nm
            ctx.setParameter("kh_bb", float(kh_bb) * _KS)
            openmm.LocalEnergyMinimizer.minimize(ctx, float(tolerance), int(steps))
            _zero(ctx)
            energy = _energy_raw(ctx)

        out_nm = np.asarray(
            ctx.getState(getPositions=True).getPositions(asNumpy=True)
            .value_in_unit(unit.nanometer), dtype=float)
        out = _split(H, out_nm)
        ca_in = np.asarray(coords["CA"], float)
        d = out["ca"] - ca_in
        out.update({
            "rung": name, "energy": float(energy), "energy_initial": float(e0),
            "eps": float(eps_ca), "k_flat": float(kf_ca), "kh_ca": float(kh_ca),
            "kh_bb": float(kh_bb), "hard": bool(hard),
            #: realised Ca displacement, in the INPUT frame (no superposition): this is the
            #: quantity the constraint bounds, and it is not an RMSD.
            "ca_disp_rms": float(np.sqrt((d ** 2).sum(1).mean())),
            "ca_disp_max": float(np.sqrt((d ** 2).sum(1)).max()),
            #: constraint violation: how far the flat-bottom wall was actually pushed through
            "ca_violation": float(max(0.0, np.sqrt((d ** 2).sum(1)).max() - eps_ca)),
            "wall": time.time() - t0,
        })
        out.update(self.am.convergence_flags(energy, e0))
        return out


    # -- single point ------------------------------------------------------
    def energy_point(self, coords: Dict[str, np.ndarray],
                     components: bool = False) -> Dict[str, float]:
        """Genuine ff14SB/GBn2 single point on `coords` -- NO minimisation at all.

        The structure is assembled exactly as a minimisation arm assembles it (heavy atoms
        from the backbone dict via `build_full_structure`, hydrogens in frozen local frames),
        so a single point and a relaxation are read on the same object.
        """
        from openmm import unit
        H = self.H
        if self._soft is None:
            self._make_soft()
        ctx = self._soft[0]
        pos = H._assemble(H._heavy_positions(coords, chi1=None))
        ctx.setPositions(pos * unit.nanometer)
        _zero(ctx)
        out = {"energy": _energy(ctx)}
        if components:
            from core.amber import AMBER_TERMS, _TERM_GROUP
            for t in AMBER_TERMS:
                out[t] = ctx.getState(getEnergy=True, groups={_TERM_GROUP[t]}) \
                    .getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        return out


_KS = None      # set below from core.amber._K_SCALE


def _zero(ctx):
    for p in ("kh_ca", "kf_ca", "eps_ca", "kh_bb"):
        ctx.setParameter(p, 0.0)


def _energy(ctx):
    from openmm import unit
    return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole)


def _energy_raw(ctx):
    """Potential energy with every restraint switched off (caller must have zeroed them)."""
    return _energy(ctx)


def _split(H, pos_nm):
    n = len(H.sequence)
    ang = pos_nm * 10.0
    bb = {}
    for nm in ("N", "CA", "C", "O"):
        bb[nm] = np.array([ang[H._heavy_index[(i, nm)]] for i in range(n)])
    cb = []
    for i in range(n):
        k = (i, "CB")
        cb.append(ang[H._heavy_index[k]] if k in H._heavy_index else bb["CA"][i])
    bb["CB"] = np.array(cb)
    heavy_names = list(H._heavy_names)
    heavy = np.array([ang[H._heavy_index[k]] for k in heavy_names])
    return {"ca": bb["CA"], "backbone": bb, "heavy": heavy, "heavy_names": heavy_names}


def _init_scale():
    global _KS
    if _KS is None:
        from core.amber import _K_SCALE
        _KS = _K_SCALE
    return _KS


# --------------------------------------------------------------------------- validity
def validity(bb: Dict[str, np.ndarray], seq: str) -> Dict[str, float]:
    """The full stereochemistry panel, on the SAME structure whose RMSD is quoted.

    `s16.energy_lib.panel` -- wrapped torsions (the Sprint-16 unwrapped-torsion defect is
    fixed there), three-way Ramachandran, clash counts, minimum heavy separation, bond and
    angle strain, cis fraction, omega deviation, chirality.
    """
    c = {a: np.asarray(bb[a], float) for a in ATOMS if a in bb}
    return EL.panel(c, seq)


VALIDITY_KEYS = ("rama_favoured", "rama_allowed", "rama_outlier", "n_clash_2A",
                 "n_clash_2p6A", "min_heavy", "bond_strain", "angle_strain",
                 "geom_rms_rel_dev", "cis_frac", "omega_dev", "chirality_L_frac")


# --------------------------------------------------------------------------- controls
HELIX_PHI, HELIX_PSI = -57.0, -47.0


def constant_backbone(n: int, phi_deg: float, psi_deg: float) -> Dict[str, np.ndarray]:
    """A constant-torsion ideal-geometry backbone of length n.

    UNITS TRAP, found the hard way in this workstream and recorded rather than silently
    fixed: `core.geometry.build_backbone` takes torsions in RADIANS (`psi[i] + math.pi`
    inside it).  Passing degrees builds a DIFFERENT constant conformation -- a helix
    specified as (-57, -47) degrees comes back with recomputed torsions (-25.9, -172.9)
    degrees and scores Ramachandran-favoured 0.000 rather than 1.000.  Everything that
    feeds this builder must convert.
    """
    from core import geometry as geo
    return geo.build_backbone(np.full(n, np.radians(phi_deg)),
                              np.full(n, np.radians(psi_deg)))


def helix_backbone(n: int) -> Dict[str, np.ndarray]:
    """ZERO-INFORMATION reference: a constant ideal alpha-helix of the same length.

    It scores Ramachandran 1.000 and zero clashes BY CONSTRUCTION.  No validity statistic
    is quotable without it.
    """
    return constant_backbone(n, HELIX_PHI, HELIX_PSI)


STRAND_PHI, STRAND_PSI = -139.0, 135.0


def strand_backbone(n: int) -> Dict[str, np.ndarray]:
    """ZERO-INFORMATION reference: a constant ideal beta-strand of the same length."""
    return constant_backbone(n, STRAND_PHI, STRAND_PSI)


def matched_random_disp(ca: np.ndarray, rms: float, rng) -> np.ndarray:
    """A matched-magnitude isotropic random Ca displacement (rigid part removed)."""
    ca = np.asarray(ca, float)
    v = rng.normal(size=ca.shape)
    v = v - v.mean(0, keepdims=True)                     # remove translation
    # remove the infinitesimal rotation component about the centroid
    X = ca - ca.mean(0, keepdims=True)
    for k in range(3):
        w = np.zeros(3); w[k] = 1.0
        g = np.cross(np.tile(w, (len(X), 1)), X)
        ng = (g ** 2).sum()
        if ng > 1e-12:
            v = v - (v * g).sum() / ng * g
    cur = float(np.sqrt((v ** 2).sum(1).mean()))
    return ca + v * (rms / max(cur, 1e-12))


def random_rigid(rng):
    """A proper rotation (det = +1) and translation.  Reflections excluded -- a reflection
    is not a symmetry of ff14SB."""
    return EL.random_rigid(rng)


# --------------------------------------------------------------------------- statistics
def paired(a, b, folds=None, names=None, seed=0, n_boot=4000):
    """Paired a - b with an i.i.d. bootstrap CI and a FOLD-CLUSTERED twin beside it."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    #: Drop pairs where either side is undefined (e.g. a recall with no near-native member).
    #: Reported as `n`, so a shrinking n is visible rather than silent.
    m = np.isfinite(a) & np.isfinite(b)
    if not m.all():
        folds = None if folds is None else np.asarray(folds)[m]
        names = None if names is None else [names[i] for i in np.flatnonzero(m)]
        a, b = a[m], b[m]
    d = a - b; n = len(d)
    if n < 3:
        return {"n": int(n), "mean": float("nan"), "median": float("nan"),
                "ci": [float("nan"), float("nan")], "W": 0, "L": 0,
                "mean_a": float("nan"), "mean_b": float("nan")}
    rng = np.random.default_rng(seed)
    bs = np.array([d[rng.integers(0, n, n)].mean() for _ in range(n_boot)])
    out = {"n": int(n), "mean_a": float(a.mean()), "mean_b": float(b.mean()),
           "mean": float(d.mean()), "median": float(np.median(d)),
           "ci": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
           "W": int((d < 0).sum()), "L": int((d > 0).sum())}
    if folds is not None:
        f = np.asarray(folds); uf = np.unique(f)
        idx = {u: np.flatnonzero(f == u) for u in uf}
        cb = []
        for _ in range(n_boot):
            pick = rng.integers(0, len(uf), len(uf))
            sel = np.concatenate([idx[uf[p]] for p in pick])
            cb.append(d[sel].mean())
        cb = np.array(cb)
        out["ci_fold"] = [float(np.percentile(cb, 2.5)), float(np.percentile(cb, 97.5))]
        out["per_fold"] = {int(u): float(d[idx[u]].mean()) for u in uf}
    if names is not None:
        o = np.argsort(d)
        out["top10"] = [(str(names[k]), float(d[k])) for k in o[:10]]
    return out


def gated(a, b, energies, **kw):
    """A gated comparison AND its ungated twin AND the exclusion count -- inseparable."""
    e = np.asarray(energies, float)
    ok = np.isfinite(e) & (e <= 1000.0)
    out = {"gated": paired(np.asarray(a)[ok], np.asarray(b)[ok],
                           folds=(np.asarray(kw["folds"])[ok] if kw.get("folds") is not None else None),
                           names=([kw["names"][i] for i in np.flatnonzero(ok)]
                                  if kw.get("names") is not None else None),
                           seed=kw.get("seed", 0)),
           "ungated": paired(a, b, folds=kw.get("folds"), names=kw.get("names"),
                             seed=kw.get("seed", 0)),
           "n_excluded": int((~ok).sum())}
    if kw.get("names") is not None:
        out["excluded"] = [str(kw["names"][i]) for i in np.flatnonzero(~ok)]
    return out


_init_scale()
