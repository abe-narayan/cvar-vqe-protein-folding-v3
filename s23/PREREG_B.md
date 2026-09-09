# PREREG_B — WORKSTREAM B, TORSION SPACE & GEOMETRY REPAIR (Sprint 23)

Never edited after data is read. Addenda are dated and appended, never in-place edits.

Git commit at write time: `a15406c82245412bc151fca51a6140401135e8a1`.

---

## PART 0 — THE AUDIT (H3), DONE FIRST, REPORTED BEFORE ANY NEW RUN

**Question:** did the "coordinate averaging beats torsion averaging by 1.024 Å" result
(`s14/avgspace.py`, arm `averaging_space_alone`, memory file `averaging-space-beats-the-objective.md`)
use a circular mean of angles or a naive arithmetic mean across the ±π branch cut?

**Source read:** `s14/avgspace.py:82-83` calls `R.circ_mean(PHI, axis=0)` /
`R.circ_mean(PSI, axis=0)`, where `R = s14.retprior`. `s14/retprior.py:53-66`:

```python
def circ_mean(A, w=None, axis=0):
    A = np.asarray(A, float)
    if w is None:
        c, s = np.cos(A).mean(axis), np.sin(A).mean(axis)
    else:
        ...weighted version...
    return np.arctan2(s, c)
```

**This is `atan2(mean sin, mean cos)` — a genuine circular mean, not a naive arithmetic mean.**
Reproduced from the artefact: `s14/results/avgspace.json`, contrast `averaging_space_alone`
(`C_torsion_mean_build` vs `B_coord_avg_raw`) = **+1.0241 [+0.6984,+1.3647], n=126, 41W/85L** —
exactly the memory figure, confirming this IS the experiment in question and it is not a
transcription of a different, buggier run.

**DISPOSITION: H3 STAYS CLOSED.** Per the BRIEF's own branching rule ("If it was circular, the
direction stays closed; report that and move all your effort to (2)"), no von-Mises-weighted or
basin-conditioned circular-mean rerun is registered here. The 1.024 Å loss is a real property of
averaging in torsion space with a correct circular mean on this pool, not a branch-cut bug.

**One qualification, stated but not chased under this budget (H3 is closed, not reopened by it):**
`circ_mean` is *unweighted* and *not basin-conditioned* — a circular mean of a genuinely bimodal
per-residue distribution lands between the modes, which is exactly the BRIEF's own warning. The
audit's job was to rule out a *bug*; it does not certify that basin-conditioning is worthless, only
that the closure was not an artefact of the branch-cut. Reopening basin-conditioned torsion
averaging is out of scope for this sprint's budget given H5 is the higher-priority, lower-risk item
and AMBER/OpenMM is the serialised bottleneck this workstream owns; named here so it is not silently
dropped.

**No live coordinator channel was available in this environment (no `ListAgents`/multi-agent
session found) to send this result ahead of the run as the BRIEF requested — it is reported here,
first, timestamped before any H5 result was read, and repeated at the top of `agentB_FINDINGS.md`.**

---

## PART 1 — H5: Cα-PRESERVING GEOMETRY REPAIR

**Mechanism read from source before any run.** `core/amber.py` — `RESTRAINED_BACKBONE = ("N", "CA",
"C")` (line 831), consumed by `AmberHamiltonian._index_topology` to build `self._restraint_idx`,
which is a **module-level global read once at Hamiltonian-construction time**, not a per-call
argument to `refine_coords`. The deployed/incumbent restraint (`core/pipeline.py:_relax_inner`,
`s16`'s canon) therefore restrains **N, C, and Cα together** — the whole rigid peptide-plane frame —
not Cα alone. The potential is `0.5*k*|r-r0|^2` per restrained atom (kcal/mol/Å²), so at fixed k an
atom pulled by internal force F settles at displacement `d = F/k` (module docstring, lines
1269-1280). `K_WEAK=1` (`d~0.1-1.0`Å), `K_MODERATE=10` (`d~0.03-0.3`Å, the pipeline default),
`K_STRONG=100` (`d~0.01-0.1`Å, the BRIEF's cited deployed repair constant). **Nobody has restrained
Cα alone and let N/C relax**, which is the mechanistic difference this experiment isolates: with
N and C free, the peptide-plane geometry (bond lengths/angles, ω-planarity) can relax to its ff14SB
optimum without being fought by a frozen N/C position, while the Cα trace itself — the thing
Cα-RMSD scores — stays pinned.

**Restraint-scope implementation.** `RESTRAINED_BACKBONE` is read at `AmberHamiltonian.__init__` →
`_build_topology` → `_index_topology`. Bypassing the shared `builder_for`/`_BUILDERS` cache (which
keys on sequence/representation/platform/threads only, not restraint scope — reusing it across
scopes would silently serve a builder frozen with the wrong restrained-atom set), this experiment
constructs its own two `AmberHamiltonian` objects per target — one with
`core.amber.RESTRAINED_BACKBONE` monkey-patched to `("CA",)`, one left at the module default
`("N","CA","C")` — restored immediately after each construction, and evaluated through
`core.amber._run_memo` (the same function `refine_coords` calls) so every returned field
(`ca`, `backbone`, `energy`, `energy_initial`, `restraint_rmsd`, `converged`, `converge_reason`,
`wall`, `components`) is byte-identical in provenance to the public API path.

**Basis.** The **point cloud** — `top75_windows(pdb)` → `I.coordinate_average` (the incumbent's own
3.048 Å arm, `s14.avgspace`'s `B_coord_avg_raw`) — rebuilt to a full backbone via
`s15.phys_repl.averaged_backbone_from(W, PHI, PSI)`, which returns `avg` (N/CA/C/O/CB dict whose CA
matches the point-cloud CA to within `dev`, reported per target) and `C_ca` (the raw point-cloud CA
itself). **This is a point cloud on both sides of every comparison it enters and is never compared
directly to a rebuilt/projected chain's RMSD** — where the ideal-geometry projection (`I.project`,
arm A in `s14/s16`'s notation) appears, it is reported as its own separately-labelled row, not
blended into the point-cloud table.

**Targets.** The pre-registered, native-free, fold×length-stratified subsample from
`s16.repair.shape_subsample(n=30)` — reused verbatim (same construction, same seed via
`stable_rng("s16repair","shape_subsample",30)`) rather than drawn fresh, both to avoid a second
best-of-K subsample choice this sprint and because it is already the project's audited convention
for a repair-family shape check. **n=30 is a shape/trade-off characterisation sample, not a
126-target claim of the mean-RMSD primary endpoint** — stated up front so a 30-target CI is never
misread as moving the BRIEF's headline number.

**Restraint ladder.** `k ∈ {0 (none), 1, 3, 10, 30, 100, 300, 1000, 3000}` kcal/mol/Å², log-spaced
across three orders of magnitude either side of `K_MODERATE`, landmarks `K_WEAK=1`,
`K_MODERATE=10`, `K_STRONG=100` included so every point on the ladder is comparable to a value this
project has already used elsewhere. `k=0` is restraint-scope-invariant (no restraint atoms move
differently at k=0) and is run once, not twice, per target. `steps=0` (minimise to convergence,
OpenMM's own convention) throughout, matching `s16`'s canon and the pipeline default
`amber_steps=0` — depth is not swept here (out of scope; `s16` pass C already covers minimisation
depth at one k and found no live interaction, and this workstream's budget is the restraint-scope
axis, not depth).

**Six-axis operator fork (Rule 0), declared before any RMSD is read:**

| axis | DECLARED | NOT TAKEN, and why |
|---|---|---|
| functional | Cartesian harmonic positional restraint, `0.5·k·\|r-r0\|²` per restrained atom, restrained SET = Cα only | (a) the incumbent's N+CA+C set — kept as the matched comparator, run at every k on the same ladder, not dropped; (b) a virtual Cα–Cα *distance* restraint (bond-length-only, no absolute position) — would target the contraction defect (3.12 vs 3.82 Å) more directly than an absolute-position tether, but changes the functional form entirely and is named as the natural next experiment if H5 shows the position tether under-repairs the bond-length axis specifically |
| basis | point cloud (`top75` coordinate average, rebuilt to a full backbone via `averaged_backbone_from`) | the ideal-geometry PROJECTED chain (arm A) as the *input* to repair — not run as an input because the projected chain already has near-ideal bonded geometry by construction, so repairing it would mostly measure noise; it IS reported as a separate reference row (already-valid-geometry ceiling) |
| readout | SEVEN validity axes reported separately (Cα-displacement, Cα-RMSD-to-native [ORACLE, eval-only], clash count, cis-peptide fraction, Ramachandran-outlier fraction, ω deviation from 180°, bond-length deviation, bond-angle deviation) — never collapsed into one score | a single composite "validity score" (e.g. weighted sum) — explicitly forbidden by the BRIEF ("REPORT VALIDITY AS SEPARATE AXES, NEVER COLLAPSED") |
| normalisation | explicit force-constant ladder (kcal/mol/Å², log-spaced) — an EXPLICIT-k sweep | an IMPLICIT sweep parameterised by target Cα-displacement (choose k per target to hit a fixed d) — would hide the functional's own d=F/k relationship inside a per-target fit and cost a root-find per target per structure; the explicit-k ladder already traces the same curve and is cheaper and more transparent |
| null | (1) `k=0` = the raw point cloud, unrepaired — prices the ACCURACY-best/VALIDITY-worst end; (2) the incumbent's own N+CA+C restraint at the SAME k values — prices whether CA-only is actually different from the deployed scope rather than just "any restraint helps"; (3) a rotational-invariance spot-check (3 targets, 2 lab frames) — position restraints are already established rotation-covariant machinery (`s16`), so this is a correctness spot-check on the CA-only variant specifically, not a full re-derivation | a size-matched random-displacement control (`s16`'s `rand` arm) — not rerun here because this experiment's question is "does restricting the restrained SET change the trade-off," which the matched N+CA+C comparator answers directly; `rand` is already priced in `s16` for the general "does AMBER move toward the truth" question |
| THE LABEL | Cα-RMSD to native is ORACLE, reported for evaluation only, NEVER used to pick a k — the seven axes above are what a deployment choice would condition on, all native-free except the RMSD row itself | picking a "winning" k by RMSD and promoting it as a selector input — forbidden by BRIEF Pillar 1 (classical methods are refinement/diagnostic, never the selector) |

**Hypothesis.** Across the k-ladder, the Cα-only restraint traces a Pareto-dominant or
Pareto-comparable frontier against the matched N+CA+C restraint on the validity axes at equal or
lower Cα-displacement/Cα-RMSD cost — i.e., for at least one k, Cα-only repairs clashes/Ramachandran/
ω-planarity/bond-angle strain by a margin exceeding a matched-SE threshold while moving the Cα trace
less (or the RMSD-to-native less) than the N+CA+C restraint does at any k.

**Falsifier.** No k on the Cα-only ladder beats the N+CA+C ladder's Pareto frontier (i.e., for every
Cα-only point, some N+CA+C point achieves equal-or-better validity at equal-or-lower Cα displacement)
— in which case restricting the restrained set buys nothing over the incumbent's own scope, and the
result is reported as such rather than reframed post hoc. **A flat Cα-RMSD curve across k while
validity improves is explicitly NOT a failure** — the BRIEF states this directly: "a repair that
holds Cα-RMSD flat while fixing clashes is a real result even if RMSD does not improve." The
falsifier is about the VALIDITY trade-off against the matched comparator, not about beating RMSD.

**Promotion.** None from this file alone; a win needs the ablation the BRIEF's Hard Rules require
("ablate anything that wins") — at minimum, re-running the winning k on a second, disjoint stratified
subsample before any claim reaches the 126-target instrument.

---

## LOGGING (BRIEF §5)

experiment id `agentB_h5_carestraint` · git commit `a15406c82245412bc151fca51a6140401135e8a1` ·
seed `stable_rng("s16repair","shape_subsample",30)` for target selection, no RNG in the repair
runs themselves (AMBER path is deterministic at Threads=1, CPU platform, DeterministicForces) ·
target ids: `s16.repair.shape_subsample(30)`, printed into the results artefact · config: k-ladder
above × 2 restraint scopes + 1 unrestrained null, `steps=0`, `tolerance=1.0`, CPU/threads=1 ·
Hamiltonian: `AmberHamiltonian` (ff14SB/GBn2), never the selector · ansatz/qubits/optimiser/α/
shots: n/a (classical repair only) · ensemble size: n/a (one point-cloud structure per target) ·
post-processing: none beyond the restrained minimisation itself · AMBER settings: as above ·
repair settings: as above · mean/median/worst Cα-RMSD, CI, per-target delta vs incumbent (k=0)
reported per k per scope in `s23/results/agentB_h5_carestraint_k{K}_scope{SCOPE}.json` and rolled
up in `s23/results/agentB_h5_summary.json`.

**AMBER/OpenMM SERIALISATION.** This workstream is the sole authorised lane this sprint. Start/
release of the AMBER path is announced in-band (no live inter-agent channel was reachable from this
session) in the accompanying message and in `agentB_FINDINGS.md`'s header.

Dated 2026-09-08, before any H5 run is executed.

---

## ADDENDUM, dated 2026-09-08, after the n=30 sweep completed (30/30, 0 errors)

**The falsifier as strictly stated does NOT fire, but not because Cα-only wins outright.** The
two restraint scopes are **Pareto-incomparable**: Cα-only Pareto-dominates N+Cα+C on bond-length/
bond-angle deviation and on convergence rate at k≳30-100 (unanimous 30W/0L on bond-length
deviation by k=1000); N+Cα+C Pareto-dominates Cα-only on ω-planarity/cis-fraction/Rama-outlier
fraction at the same k range (cis fraction reaches 64% under Cα-only at k=3000 against ~0% for
N+Cα+C). Cα displacement shows a significant CROSSOVER at k≈30 (Cα-only moves the trace MORE
than N+Cα+C below that k, LESS above it) with a clean mechanistic explanation (§ full derivation
in `agentB_FINDINGS.md`: N+Cα+C's geometric over-determination stabilises Cα even without a
direct restraint on it at low k). Cα-RMSD-to-native is a wash (mostly null, sign inconsistent
across k). Clash count never differed (both scopes always clash-free on this basis). Full results,
every paired comparison with iid+fold-clustered CI and MDE: `agentB_FINDINGS.md` Part 1.

This is reported as the pre-registered falsifier's own logic requires: a genuine result either
way, not reframed post hoc into a clean win for the hypothesis.
