"""s21/a_matrix.py -- THE MANDATORY MATRIX (BRIEF section 4), with genuine CVaR-VQE.

    Legacy | AMBER | Legacy+AMBER | Distance | Distance+Legacy | Distance+AMBER | D+L+A

Same candidate problem, same representation, same budget, same seeds, every component
SEPARABLY EVALUABLE at all times: `MultiHam` keeps each genuine component's RAW value for every
candidate it ever scores, so Legacy-only, AMBER-only, hybrid and structural-only are all
readable off one run and nothing is buried in an opaque score.

WHAT IS GENUINE HERE, and where each piece comes from.

  H_Legacy   `core.energy` eleven components at DEFAULT_WEIGHTS, never fitted, via
             `s16.energy_lib.legacy_components_of_windows` on the ideal-geometry rebuild.
  H_AMBER    genuine ff14SB/GBn2 SINGLE POINT via `s20.qb2_lib.AmberSP`, which is asserted
             bit-exact against `core.amber.refine_coords(k_restraint=0, steps=-1)` by its own
             `verify()`.  NO minimisation: H_AMBER = E, not E o Relax_50.
  H_Distance the shipped leave-fold-out distogram's inverse-variance-weighted squared
             deviation -- `s20.qb2_lib.Ham("DIST")`, the deployed structural objective.
  selector   `s20.qb2_opt.arm_vqe` -- `core.quantum.MPSAnsatz` (RY + CNOT chain + final RY)
             simulated EXACTLY, sampled, with the analytic CVaR score-function gradient at
             `baseline="const"` (the CORRECTED estimator, never the recorded `tail` defect).
             The bitstring selects a per-residue conformer basin and the ANGLE IS DRAWN
             CONTINUOUSLY from that basin's von Mises component -- a continuous density with
             full support on the torus.  No lattice.

READOUT.  Workstream D verified in source what the pillar emits: `vqe_bitstring` (argmin over
the final distribution's samples), `vqe_modal_bitstring` (the mode), `best_seen_bitstring`
(argmin over everything seen).  NONE is an average of a tail.  **CVaR is the TRAINING
objective; the READOUT is an argmin.**  So the primary readout here is `Field.best_z`, the
argmin of the training Hamiltonian over every candidate the arm evaluated.  An AVERAGED tail
readout is also emitted and is LABELLED AVERAGED throughout -- it is a different (Sprint-19)
lever and must never be compared to a single-structure arm as though it were one.

NORMALISATION -- DECLARED BEFORE ANY RMSD WAS SEEN (PREREG_A stage 1, BRIEF section 6).

  Stage 1 measured the three energies over each target's own K=500 retrieval pool.  AMBER's
  per-target range is ~13 decades with a median 41% of members above 1e4 kcal/mol; Legacy and
  Distance are within one or two decades.  A raw sum is therefore an AMBER-outlier detector
  wearing a hybrid's name.

  DECLARED FORM, per component, fitted FREE and NATIVE-FREE on the target's own K=500 pool:

      z_c(e) = ( cond_c(e) - median_pool(cond_c) ) / ( IQR_pool(cond_c) / 1.349 )
      cond_AMB = signed log about the pool median;  cond_DIST = cond_LEG = identity
      H_combined = sum over components of z_c

  THE PROPERTY THAT MAKES THIS DEFENSIBLE, and it is an identity, not a result:
  `cond` is STRICTLY MONOTONE, so for a SINGLE-component Hamiltonian every argmin readout is
  EXACTLY unchanged by it -- `argmin cond(E) == argmin E`.  Conditioning cannot flip a
  single-Hamiltonian ranking; it changes only the CVaR value, the gradient and the relative
  weight a component carries inside a SUM.  The normalisation is therefore a decision about
  the HYBRID cells only, and it is declared here before any RMSD was read.

  AUDIT arms, computed alongside and reported whatever they say: `robustz` (the same robust
  standardisation with NO conditioning -- Sprint 20's primary) and `rank` (pool-quantile
  rank-to-normal, the form `s21/tailprice.py` declared).  If an audit arm wins, the declared
  choice is recorded as WRONG, not quietly swapped.

BUDGET.  Matched in CANDIDATE EVALUATIONS, which is the space the SELECTOR works in (BRIEF
section 7 rule 1).  One candidate costs one unit whether the Hamiltonian has one component or
three, so a hybrid gets no extra look at the space -- it only gets a different opinion about
what it saw.  Wall-clock is NOT matched and is reported separately.

CONTROLS, on every quantum claim (BRIEF section 7):
  vqe_untrained  best-of-N from the UNTRAINED circuit, identical shots and budget, theta never
                 stepped.  Never an initialisation mean.
  best_of_N      B i.i.d. draws from the same von Mises basin mixtures.  Zero optimisation.
  metro          continuous single-residue Metropolis in torsion space, matched move class.
  helix          ZERO-INFORMATION but plausible: constant ideal alpha-helix plus matched
                 isotropic torsion noise.  Uniform-on-the-torus is NOT a zero-information
                 control (BRIEF section 7 rule 4).

OPERATOR FORKS (BRIEF section 7, rule 0 -- added 2026-09-07).  This lane has a DIRECTIONAL
hypothesis: falsifier F-A1 expects no Legacy/AMBER Hamiltonian to beat the structural arm.  So
every fork is enumerated and THE ALTERNATIVE NOT TAKEN IS NAMED, with the direction it would
have pushed.  I designed these comparisons and I have a stake in their direction, so this list
is offered to the coordinator for independent fork review rather than treated as sufficient.

  FUNCTIONAL   H_AMBER is the BARE single point.
               NOT TAKEN: `E o Relax_50`, the deployed object (Sprint 20 L7c: the relaxation is
               CONSTITUTIVE).  Relaxation would very likely make AMBER LESS catastrophic, so
               this fork points TOWARD my expected conclusion and is the one to distrust most.
  BASIS        energies on the ideal-geometry rebuild; RMSD on `build_ca_exact` of the same
               torsions.  NOT TAKEN: scoring the distogram on the rebuild too (54/75 agreement
               with the shipped filter, vs 75/75 on the window).  Direction: unknown; measured
               instead -- the pool run carries BOTH bases and they agree to ~0.06 A.
  READOUT      four readouts emitted, NONE privileged.  NOT TAKEN: quoting `tail_avg` alone,
               which is exactly where Legacy looks worst -- that would have pushed toward my
               expected conclusion, which is why the readout column is mandatory here.
  NORMALISATION signed-log + robust z on the target's own pool.  NOT TAKEN: the raw sum and the
               pool-quantile rank transform.  Both are computed as declared audits and reported
               whatever they say.  Direction: the raw sum is an AMBER-outlier detector and would
               have pushed toward my expected conclusion.
  NULL         matched-count random subsets drawn from THE ARM'S OWN generated set.
               NOT TAKEN: uniform-on-the-torus (BRIEF section 7 rule 4 -- not zero-information)
               and the initialisation MEAN (forbidden; best-of-N from the untrained circuit is
               used instead).
  SEED COUNT   4 seeds on variational arms, 2 on the non-variational controls.
               **THIS FORK WAS OMITTED FROM MY FIRST LIST and added on coordinator review.**
               NOT TAKEN: more seeds, or quoting best-of-seeds.  It is the fork that most
               constrains the table: seed sd is 0.51-1.06 A (2-5x Sprint 20's 0.200 at 16x the
               budget) and best-of-4 beats the 4-seed mean by 0.65-1.36 A, which is LARGER than
               most row-to-row gaps.  Consequence: only the DIRECTION of the matrix is safe;
               every row-to-row MAGNITUDE sits inside seed noise.
  BUDGET       matched in CANDIDATE EVALUATIONS.  NOT TAKEN: matching WALL CLOCK, which would
               give the single-component Hamiltonians roughly 3x more candidates than the
               three-component hybrids.  **This fork points AGAINST my expected conclusion** --
               it would favour Legacy-only and AMBER-only arms -- and it is named for that
               reason.

    python -m s21.a_matrix gate            # G-A1 .. G-A3, must pass before any number
    python -m s21.a_matrix smoke           # 2 targets, writes _SMOKE_, never a result
    python -m s21.a_matrix run [n]
    python -m s21.a_matrix report
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                   # noqa: E402
from s15 import seed as SD                        # noqa: E402
from s16 import energy_lib as EL                  # noqa: E402
from s18 import phys_lib as PL                    # noqa: E402
from s20 import qb2_lib as L                      # noqa: E402
from s20 import qb2_opt as OP                     # noqa: E402

SALT = "s21A"

#: the mandatory matrix, in the brief's own order.  Values are the genuine components.
MATRIX = {
    "Legacy":     ("LEG",),
    "AMBER":      ("AMB",),
    "Leg+Amb":    ("LEG", "AMB"),
    "Distance":   ("DIST",),
    "Dist+Leg":   ("DIST", "LEG"),
    "Dist+Amb":   ("DIST", "AMB"),
    "Dist+Leg+Amb": ("DIST", "LEG", "AMB"),
}
HAMS = tuple(MATRIX)
COMPONENTS = ("DIST", "LEG", "AMB")

#: TARGET COUNT.  DECLARED, and it is a COMPUTE cut, not a scientific one.  The shared box
#: measured 52 ms per AMBER single point against a 6 ms nominal under four-way contention; at
#: BUDGET=512 one target costs ~44,000 genuine ff14SB/GBn2 single points, so n=20 did not fit the
#: window.  n is reduced and the per-cell configuration is left intact, because a thinner cell
#: would change what is measured while a smaller n only changes the power -- and the power is
#: printed as a warning on every table.  MDE at n=12 is FAR above the instrument's 0.084 A.
N_TARGETS = int(os.environ.get("S21A_N", 12))

#: DECLARED BEFORE THE RUN.  Not swept, not tuned on RMSD.  The env overrides exist ONLY so a
#: code-path smoke can run in seconds; every reported number uses the defaults, and the config
#: (with its hash) is persisted in the artefact so a smoke can never be mistaken for a result.
BUDGET = int(os.environ.get("S21A_BUDGET", 512))   # candidate evaluations, identical every arm
SHOTS = int(os.environ.get("S21A_SHOTS", 32))      # -> 16 CVaR gradient steps at BUDGET=512
ALPHA = 0.25           # CVaR level, Sprint 20's deployed value
#: 4 seeds MINIMUM on any VARIATIONAL arm -- Sprint 20 measured within-target ansatz-seed sd at
#: 0.200 A = 2.4x MDE, so fewer is NOT MEASURED.  The non-variational controls carry no ansatz
#: seed at all; they are run at CTRL_SEEDS and that reduction is DECLARED here, not discovered.
SEEDS = tuple(range(int(os.environ.get("S21A_SEEDS", 4))))
CTRL_SEEDS = tuple(range(int(os.environ.get("S21A_CTRL_SEEDS", 2))))
NORM = "signedlog_robustz"
AUDIT_NORMS = ("robustz", "rank")
TAIL_ALPHAS = (0.05, 0.15, 0.30)

VARIATIONAL = ("vqe", "vqe_untrained")
CONTROLS = ("best_of_N", "metro", "helix")
ARMS = VARIATIONAL + CONTROLS
AUDIT = os.environ.get("S21A_AUDIT", "1") == "1"


def seeds_for(arm):
    return SEEDS if arm in VARIATIONAL else CTRL_SEEDS


# ==========================================================================
# 1.  the separably-evaluable combined Hamiltonian
# ==========================================================================
class Comp:
    """ONE genuine energy component with its own free, native-free pool calibration.

    `raw` is the genuine energy and is never modified.  `cond` is a strictly monotone
    conditioning (identity except for AMBER); `z` is the declared standardisation.  The two
    audit standardisations are computed from the same pool sample, so switching normalisation
    never re-runs physics.
    """

    def __init__(self, kind, tgt, sp=None):
        self.kind = kind
        self.t = tgt
        self.sp = sp
        self.cal = {}

    # -- the genuine energy, unbudgeted ------------------------------------
    def raw(self, phi, psi):
        k = self.kind
        if k == "LEG":
            comp = EL.legacy_components_of_windows(self.t["seq"], phi, psi)
            return np.asarray(EL.legacy_total_from(comp), float)
        if k == "AMB":
            return np.asarray(self.sp.batch(phi, psi), float)
        if k == "DIST":
            from core import project as pj
            CA = np.asarray(pj.build_ca_exact(phi, psi), float)
            rv = CA[:, self.t["i"], :] - CA[:, self.t["j"], :]
            d = np.sqrt((rv * rv).sum(-1))
            return ((d - self.t["dhat"]) ** 2 * self.t["inv2"]).sum(-1)
        raise ValueError(k)

    # -- the declared conditioning ------------------------------------------
    def cond(self, e):
        """STRICTLY MONOTONE.  For a single-component Hamiltonian, argmin is EXACTLY
        invariant under this map -- an identity, stated so it is not read as a finding."""
        e = np.asarray(e, float)
        if self.kind != "AMB":
            return e
        d = e - self.cal.get("ref", 0.0)
        return np.sign(d) * np.log1p(np.abs(d))

    def calibrate(self, PHI, PSI):
        """FREE and UNBUDGETED.  Fitted on the target's own K=500 retrieval pool, which is
        native-free.  Persists all three normalisations from ONE physics evaluation."""
        e = self.raw(PHI, PSI)
        ok = np.isfinite(e)
        ev = e[ok]
        self.cal["ref"] = float(np.median(ev)) if ok.any() else 0.0
        c = self.cond(e)[ok]
        q1, q2, q3 = np.percentile(c, [25, 50, 75])
        r1, r2, r3 = np.percentile(ev, [25, 50, 75])
        self.cal.update({
            "med": float(q2), "iqr": float(max((q3 - q1) / 1.349, 1e-12)),
            "med_raw": float(r2), "iqr_raw": float(max((r3 - r1) / 1.349, 1e-12)),
            "pool_sorted": np.sort(ev),
            "n_ref": int(len(e)), "n_nonfinite_ref": int((~ok).sum()),
            "mean": float(ev.mean()), "sd": float(ev.std()),
            "log10_range": float(np.log10(np.ptp(ev) + 1e-30)) if ok.sum() > 1 else float("nan"),
            "frac_gt_1e4": float((ev > 1e4).mean()),
            "skew": float(((ev - ev.mean()) ** 3).mean() / max(ev.std(), 1e-30) ** 3),
        })
        return self.cal

    # -- the three standardisations ----------------------------------------
    def z(self, e, norm=NORM):
        e = np.asarray(e, float)
        if norm == "signedlog_robustz":
            return (self.cond(e) - self.cal["med"]) / self.cal["iqr"]
        if norm == "robustz":
            return (e - self.cal["med_raw"]) / self.cal["iqr_raw"]
        if norm == "rank":
            from scipy.special import ndtri
            P = self.cal["pool_sorted"]
            N = len(P)
            r = np.searchsorted(P, e, side="left").astype(float)
            return ndtri(np.clip((r + 0.5) / (N + 1.0), 1e-6, 1 - 1e-6))
        raise ValueError(norm)


class MultiHam:
    """A budgeted SUM of standardised genuine components that never stops being separable.

    Duck-typed to `s20.qb2_opt.Field`'s expectations (`left`, `__call__`, `std`, `raw`), so the
    entire Sprint-20 optimiser and CVaR-VQE stack runs over it unmodified.

    THE BUDGET IS IN CANDIDATES.  One candidate costs one unit whether the Hamiltonian has one
    component or three: a hybrid gets no extra look at the space.
    """

    def __init__(self, kinds, tgt, comps, budget=10 ** 9, keep=True, norm=NORM):
        self.kinds = tuple(kinds)
        self.t = tgt
        self.c = {k: comps[k] for k in self.kinds}
        self.budget = int(budget)
        self.used = 0
        self.keep = bool(keep)
        self.norm = norm
        self._phi, self._psi = [], []
        self._raw = {k: [] for k in self.kinds}

    @property
    def left(self):
        return max(0, self.budget - self.used)

    def raw(self, phi, psi):
        """UNBUDGETED.  Returns the standardised SUM; per-component raws are available from
        `raw_components`.  `Field` calls `std()` on this, which is the identity here."""
        return self.raw_components(phi, psi)[0]

    def raw_components(self, phi, psi):
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        per = {k: self.c[k].raw(phi, psi) for k in self.kinds}
        s = np.zeros(len(phi))
        for k in self.kinds:
            s = s + self.c[k].z(per[k], self.norm)
        return s, per

    def std(self, e):
        return np.asarray(e, float)          # `raw` is already standardised

    def __call__(self, phi, psi):
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        k = min(len(phi), self.left)
        if k <= 0:
            return np.zeros(0)
        phi, psi = phi[:k], psi[:k]
        s, per = self.raw_components(phi, psi)
        self.used += k
        if self.keep:
            self._phi.append(phi.copy()); self._psi.append(psi.copy())
            for q in self.kinds:
                self._raw[q].append(np.asarray(per[q], float))
        return s

    def seen(self):
        if not self._phi:
            return np.zeros((0, 0)), np.zeros((0, 0)), {}
        return (np.concatenate(self._phi), np.concatenate(self._psi),
                {k: np.concatenate(v) for k, v in self._raw.items()})


# ==========================================================================
# 2.  one target: build the components once, calibrate once
# ==========================================================================
def build(pdb, tries=40):
    """MEASURED, not assumed: `s21/tailprice.py` was KILLED at target 13/126 by
    `core.amber.memory_guard` ("physical memory at 96% exceeds the 92% ceiling") because
    sibling workstreams filled the shared box.  The guard is right; the caller must be
    resilient.  Hold for memory and retry -- never proceed by disabling the guard."""
    import time
    tgt = L.target(pdb)
    sp = None
    for a in range(tries):
        EL.mem_hold(min_gb=1.45, tag="s21A/matrix", max_wait=300.0)
        try:
            sp = L.AmberSP(tgt["seq"], tgt["rep"])
            break
        except MemoryError:
            print(f"  [s21A] memory guard fired (attempt {a+1}/{tries}); backing off", flush=True)
            time.sleep(30)
    if sp is None:
        raise MemoryError("could not obtain an OpenMM context within the retry budget")
    comps = {k: Comp(k, tgt, sp) for k in COMPONENTS}
    cal = {k: {kk: vv for kk, vv in comps[k].calibrate(tgt["PHI"], tgt["PSI"]).items()
               if kk != "pool_sorted"} for k in COMPONENTS}
    return tgt, comps, sp, cal


# ==========================================================================
# 3.  gates -- every one reports how many times it FIRED (BRIEF section 7 rule 3)
# ==========================================================================
def gate(n_targets=3, verbose=True):
    """G-A1  AmberSP is bit-exact against `core.amber.refine_coords(k=0, steps=-1)`.
       G-A2  the conditioning is EXACTLY order-preserving on real pool energies.
       G-A3  a single-component MultiHam's argmin equals the genuine energy's argmin under
             every one of the three normalisations -- the identity the write-up leans on."""
    out = {"rows": [], "n_compare_G_A1": 0, "n_compare_G_A2": 0, "n_compare_G_A3": 0}
    worst_sp, worst_ord, worst_arg = 0.0, 0, 0
    for t in L.subset(64)[:n_targets]:
        pdb = t["pdb"]
        tgt, comps, sp, _cal = build(pdb)
        PHI, PSI = tgt["PHI"][:64], tgt["PSI"][:64]
        rel, k = sp.verify(PHI, PSI, m=4)
        worst_sp = max(worst_sp, rel); out["n_compare_G_A1"] += k
        row = {"pdb": pdb, "sp_rel": rel, "sp_n": k}
        for c in COMPONENTS:
            e = comps[c].raw(PHI, PSI)
            o1 = np.argsort(e, kind="stable")
            o2 = np.argsort(comps[c].cond(e), kind="stable")
            bad = int((o1 != o2).sum())
            worst_ord = max(worst_ord, bad); out["n_compare_G_A2"] += len(e)
            row[f"ord_break_{c}"] = bad
            for nm in (NORM,) + AUDIT_NORMS:
                h = MultiHam((c,), tgt, comps, budget=10 ** 9, keep=False, norm=nm)
                d = int(np.argmin(h.raw(PHI, PSI)) != int(np.argmin(e)))
                worst_arg = max(worst_arg, d); out["n_compare_G_A3"] += 1
                row[f"argmin_break_{c}_{nm}"] = d
        out["rows"].append(row)
    out.update({"G_A1_max_rel": worst_sp, "G_A2_max_order_breaks": int(worst_ord),
                "G_A3_max_argmin_breaks": int(worst_arg),
                "passed": bool(worst_sp < 1e-9 and worst_ord == 0 and worst_arg == 0)})
    if verbose:
        print(f"G-A1 AmberSP bit-exact: max rel {worst_sp:.3e} over "
              f"{out['n_compare_G_A1']} comparisons  (FIRED {out['n_compare_G_A1']}x)")
        print(f"G-A2 conditioning order-preserving: {out['G_A2_max_order_breaks']} breaks over "
              f"{out['n_compare_G_A2']} ranks  (FIRED {out['n_compare_G_A2']}x)")
        print(f"G-A3 single-component argmin invariant to normalisation: "
              f"{out['G_A3_max_argmin_breaks']} breaks over {out['n_compare_G_A3']} checks")
        print(f"  -> {'PASS' if out['passed'] else 'FAIL'}")
    return out


# ==========================================================================
# 4.  one arm
# ==========================================================================
def _rmsd(z, tgt):
    return float(L.rmsd_of(np.atleast_2d(z), tgt)[0])


def _avg_rmsd(PHI, PSI, idx, tgt):
    """AVERAGED readout -- the deployed coordinate average of a candidate subset.  POINT CLOUD,
    not a built structure.  Labelled AVERAGED everywhere it appears."""
    from core import project as pj
    CA = np.asarray(pj.build_ca_exact(PHI[idx], PSI[idx]), float)
    if len(idx) == 1:
        return float(I.ca_rmsd(CA[0], tgt["nat"]))
    a, _b = I.coordinate_average(CA)
    return float(I.ca_rmsd(np.asarray(a, float), tgt["nat"]))


def run_arm(name, hname, tgt, comps, cal, seed, norm=NORM, budget=BUDGET,
            cross=False, sp=None):
    kinds = MATRIX[hname]
    h = MultiHam(kinds, tgt, comps, budget=budget, keep=True, norm=norm)
    F = OP.Field(h, tgt)
    F.set_pen(12.0)                       # standardised units; identical for every arm
    rng = SD.stable_rng(tgt["pdb"], hname, name, seed, norm, salt=SALT)
    t0 = time.time()
    if name == "vqe":
        info = OP.arm_vqe(F, None, rng, alpha=ALPHA, shots=SHOTS, train=True)
    elif name == "vqe_untrained":
        info = OP.arm_vqe(F, None, rng, alpha=ALPHA, shots=SHOTS, train=False)
    elif name == "helix":
        info = _arm_helix(F, rng)
    else:
        info = OP.ARMS[name](F, None, rng)
    wall = time.time() - t0

    PHI, PSI, per = h.seen()
    B = len(PHI)
    s = np.zeros(B)
    for k in kinds:
        s = s + comps[k].z(per[k], norm)

    out = {"arm": name, "ham": hname, "seed": int(seed), "norm": norm,
           "used": int(h.used), "n_seen": int(B), "wall": wall,
           "n_nonfinite": int(F.n_nonfinite),
           "best_std": float(F.best),
           #: PRIMARY READOUT -- argmin of the TRAINING Hamiltonian over everything seen.
           #: This is the operator the deployed pillar has (`best_seen_bitstring`).
           "rmsd_ORACLE": _rmsd(F.best_z, tgt) if F.best_z is not None else float("nan")}
    for kk in ("iters", "ess_frac", "gnorm_mean", "gnorm_sd", "cos_grad_vs_alpha1",
               "cvar_first", "cvar_last", "param_disp", "latent_entropy_bits",
               "accept_rate", "steps"):
        if kk in info:
            out[kk] = info[kk]

    if B:
        order = np.argsort(s, kind="stable")
        #: ORACLE, post-hoc scoring only: the true RMSD of every candidate this arm generated.
        d = np.asarray(L.rmsd_of(L.pack(PHI, PSI), tgt), float)
        out["ORACLE_seen_best"] = float(d.min())
        out["ORACLE_seen_mean"] = float(d.mean())

        #: ================= THE READOUT COLUMN =================================
        #: COORDINATOR, 2026-09-07, measured at n=126: **Legacy's sign FLIPS with the readout**
        #: -- tail_member -0.41 [-0.56,-0.26] (BEATS its control) vs tail_medoid +0.30
        #: [+0.20,+0.40] and tail_avg +0.33 [+0.21,+0.45] (both WORSE, 35W/91L and 41W/85L).
        #: Both are true.  So "which Hamiltonian is best" is meaningless unqualified: EVERY row
        #: below states its readout, and the four readouts are reported side by side, never
        #: collapsed.  The mechanism is ERROR COHERENCE, not diversity: averaging cancels i.i.d.
        #: error and PRESERVES systematic error, and Legacy's tail shares a coherent compactness
        #: bias (0.45 A more compact, 124W/2L, Sprint 20 L2c).  The diversity explanation was
        #: REFUTED by its own sign control -- maximising set diversity is dead on all three
        #: bands and MINIMISING it helps (-0.214 [-0.397,-0.040]) -- so no diversity term is
        #: built into any arm here.
        #:
        #: Each readout gets its OWN matched-count random control drawn from THIS ARM'S OWN
        #: generated set, because a control must be matched in the space the operator works in.
        rng_r = SD.stable_rng(tgt["pdb"], hname, name, seed, "readout", salt=SALT)
        for a in TAIL_ALPHAS:
            m = max(1, int(round(a * B)))
            T = order[:m]
            out[f"tailavg{a:g}"] = _avg_rmsd(PHI, PSI, T, tgt)          # AVERAGED
            out[f"tailmedoid{a:g}"] = _medoid_rmsd(PHI, PSI, T, tgt)    # single structure
            pick = rng_r.choice(m, size=min(8, m), replace=(m < 8))
            out[f"tailmember{a:g}"] = float(np.mean(d[T[pick]]))        # single structure
            out[f"tailmin{a:g}"] = float(d[T[0]])                       # == argmin, by identity
            #: matched-count random controls, one per readout
            ra, rm, rb = [], [], []
            for _ in range(8):
                S = rng_r.choice(B, m, replace=False)
                ra.append(_avg_rmsd(PHI, PSI, S, tgt))
                rm.append(_medoid_rmsd(PHI, PSI, S, tgt))
                rb.append(float(d[S[rng_r.integers(0, m)]]))
            out[f"randavg{a:g}"] = float(np.mean(ra))
            out[f"randmedoid{a:g}"] = float(np.mean(rm))
            out[f"randmember{a:g}"] = float(np.mean(rb))
        #: legacy key names kept so nothing downstream silently reads a missing field
        for a in TAIL_ALPHAS:
            out[f"AVERAGED_tail{a:g}"] = out[f"tailavg{a:g}"]
            out[f"medoid_tail{a:g}"] = out[f"tailmedoid{a:g}"]

        #: ================= CROSS-READOUT =====================================
        #: COORDINATOR PRIORITY.  The pool-restricted, same-H, order-based case is closed
        #: analytically and confirmed at tail_min == argmin in all 9 cells, max difference
        #: exactly 0 -- so the informative cells are the ones breaking those scope conditions,
        #: above all a READOUT H DIFFERENT FROM THE TRAINING H.  Now computed on ALL FOUR
        #: variational seeds, not one, because a single seed is NOT MEASURED.
        if cross:
            miss = [k for k in COMPONENTS if k not in kinds]
            for k in miss:
                per[k] = comps[k].raw(PHI, PSI)
            Z = {k: comps[k].z(per[k], norm) for k in COMPONENTS}
            for h2 in HAMS:
                s2 = np.zeros(B)
                for k in MATRIX[h2]:
                    s2 = s2 + Z[k]
                out[f"x_global|{h2}"] = float(d[int(np.argmin(s2))])
                for a in TAIL_ALPHAS:
                    m = max(1, int(round(a * B)))
                    T = order[:m]
                    out[f"x{a:g}|{h2}"] = float(d[T[int(np.argmin(s2[T]))]])
                    #: the same cross-readout at the OTHER readouts, so the readout column
                    #: and the cross-readout matrix are the same object seen twice.
                    o2 = T[np.argsort(s2[T], kind="stable")]
                    out[f"xavg{a:g}|{h2}"] = _avg_rmsd(PHI, PSI, o2[:max(1, m // 2)], tgt)
            for k in COMPONENTS:
                out[f"ORACLE_rho|{k}"] = PL.spearman(Z[k], d)
            #: IDENTITY CHECK, printed as a check and never as a discovery: the alpha-tail of
            #: the TRAINING H contains that H's own minimum, so its diagonal cross-readout is
            #: exactly the argmin over everything seen.
            out["identity_tailmin_vs_argmin"] = float(
                abs(out[f"x{TAIL_ALPHAS[0]:g}|{hname}"] - out["rmsd_ORACLE"]))
    return out


def _medoid_rmsd(PHI, PSI, idx, tgt):
    """Single-structure readout of the same tail: its medoid.  No averaging."""
    from core import project as pj
    CA = np.asarray(pj.build_ca_exact(PHI[idx], PSI[idx]), float)
    if len(idx) == 1:
        return float(I.ca_rmsd(CA[0], tgt["nat"]))
    b = I.medoid(I.pairwise_rmsd(CA))
    return float(I.ca_rmsd(CA[b], tgt["nat"]))


def _arm_helix(F, rng, batch=64):
    """ZERO-INFORMATION but PLAUSIBLE reference: a constant ideal alpha-helix plus matched
    isotropic torsion noise, scored by the same Hamiltonian at the same budget.  BRIEF
    section 7 rule 4 -- uniform on the torus is NOT a zero-information control."""
    n = F.n
    while F.left > 0:
        k = min(batch, F.left)
        phi = L.wrap(np.full((k, n), np.deg2rad(-57.0)) + rng.normal(0, 0.35, (k, n)))
        psi = L.wrap(np.full((k, n), np.deg2rad(-47.0)) + rng.normal(0, 0.35, (k, n)))
        F(L.pack(phi, psi))
    return {}


# ==========================================================================
# 5.  the run
# ==========================================================================
def _path(out):
    return os.path.join(RESULTS, out)


CFG = {"MATRIX": {k: list(v) for k, v in MATRIX.items()}, "BUDGET": BUDGET, "SHOTS": SHOTS,
       "ALPHA": ALPHA, "SEEDS": list(SEEDS), "CTRL_SEEDS": list(CTRL_SEEDS), "NORM": NORM, "AUDIT_NORMS": list(AUDIT_NORMS),
       "ARMS": list(ARMS), "TAIL_ALPHAS": list(TAIL_ALPHAS), "N_TARGETS": N_TARGETS,
       "cross_readout_on": "vqe, all 4 seeds",
       "readout": "argmin of the training H over everything seen (Field.best_z)",
       "selector": "s20.qb2_opt.arm_vqe -- genuine CVaR-VQE, continuous torsions"}


def _required(row):
    """A row is complete only if EVERY Hamiltonian x arm x seed cell it promises is present,
    plus the cross-readout on the declared seed.  The flag must require the FULL configuration,
    not the subset that happened to be called (BRIEF section 6)."""
    need = [f"{h}|{a}|{s}" for h in HAMS for a in ARMS for s in seeds_for(a)]
    if not all(k in row.get("cells", {}) for k in need):
        return False
    if not all(f"x_global|{h2}" in row["cells"][f"{h}|vqe|{s}"]
               for h in HAMS for h2 in HAMS for s in SEEDS):
        return False
    if not all(f"{p}{a:g}" in row["cells"][f"{h}|{arm}|{s}"]
               for h in HAMS for arm in ARMS for s in seeds_for(arm) for a in TAIL_ALPHAS
               for p in ("tailavg", "tailmedoid", "tailmember", "tailmin",
                         "randavg", "randmedoid", "randmember")):
        return False
    #: the AUDIT normalisations are part of the declared configuration, so when they are enabled
    #: their cells are part of COMPLETE.  A flag that ignores a configured block is the same
    #: hazard as one that counts skipped rows -- it certifies a subset as the whole.
    if AUDIT and not all(f"AUD{nm}|{h}|vqe|{SEEDS[0]}" in row["cells"]
                         for nm in AUDIT_NORMS
                         for h in ("Leg+Amb", "Dist+Leg", "Dist+Amb", "Dist+Leg+Amb")):
        return False
    return True


def _write(rows, n_expected, out):
    n_ok = sum(1 for r in rows if _required(r))
    obj = {"rows": rows, "config": CFG, "cfg_hash": PL.cfg_hash(CFG),
           "n_rows": len(rows), "n_expected": int(n_expected),
           "n_complete_rows": n_ok,
           "complete": bool(n_ok == int(n_expected) and len(rows) == int(n_expected))}
    p = _path(out)
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, p)


def run(n=None, out="a_matrix.json", audit=True):
    n = N_TARGETS if n is None else n
    tg = L.subset(n)
    rows = []
    if os.path.exists(_path(out)):
        try:
            rows = [r for r in json.load(open(_path(out))).get("rows", []) if _required(r)]
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for t in tg:
        pdb = t["pdb"]
        if pdb in done:
            continue
        tgt, comps, sp, cal = build(pdb)
        rel, nk = sp.verify(tgt["PHI"], tgt["PSI"], m=3)
        rec = {"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]),
               "cal": cal, "G_A1_rel": float(rel), "G_A1_n": int(nk), "cells": {}}
        d0 = np.asarray(L.rmsd_of(L.pack(tgt["PHI"], tgt["PSI"]), tgt), float)
        rec["ORACLE_pool_mean"] = float(d0.mean())
        rec["ORACLE_pool_best"] = float(d0.min())
        for hname in HAMS:
            for arm in ARMS:
                for s in seeds_for(arm):
                    r = run_arm(arm, hname, tgt, comps, cal, s,
                                cross=(arm == "vqe"), sp=sp)
                    rec["cells"][f"{hname}|{arm}|{s}"] = r
        #: AUDIT normalisations, on the HYBRID cells only (the identity above makes the
        #: single-component cells invariant), one seed, declared before the run.
        if audit and AUDIT:
            for nm in AUDIT_NORMS:
                for hname in ("Leg+Amb", "Dist+Leg", "Dist+Amb", "Dist+Leg+Amb"):
                    r = run_arm("vqe", hname, tgt, comps, cal, SEEDS[0], norm=nm)
                    rec["cells"][f"AUD{nm}|{hname}|vqe|{SEEDS[0]}"] = r
        rows.append(rec)
        _write(rows, len(tg), out)
        print(f"  {len(rows)}/{len(tg)}  {pdb}  ({time.time()-t0:.0f}s)", flush=True)
    _write(rows, len(tg), out)
    print(f"DONE {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)
    return rows


if __name__ == "__main__":
    a = sys.argv[1] if len(sys.argv) > 1 else "run"
    if a == "gate":
        g = gate()
        json.dump(g, open(_path("a_matrix_gate.json"), "w"),
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
        sys.exit(0 if g["passed"] else 1)
    elif a == "smoke":
        run(n=2, out="_SMOKE_a_matrix.json")
    elif a == "report":
        from s21 import a_report                                   # noqa: F401
        a_report.main()
    else:
        run(n=int(sys.argv[2]) if len(sys.argv) > 2 else 20)
