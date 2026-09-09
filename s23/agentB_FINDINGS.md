# agentB_FINDINGS — WORKSTREAM B, Sprint 23: torsion space & geometry repair

Full pre-registration: `s23/PREREG_B.md` (written and timestamped before any run). Read-first
sources per the BRIEF: `s23/BRIEF.md`, `s22/LEDGER.md`, `s22/CLAIMS.md`, `s21/LEDGER.md`,
`core/project.py`, `s12/instrument.py`, `core/amber.py`, `core/geometry.py`.

**No live inter-agent channel was reachable from this session** (`ListAgents` is not an
available tool here) — the circular-mean audit result below, which the BRIEF asked to be sent
"as soon as you have it," is instead reported here, first, and was determined and written down
**before** any H5 run started.

---

## PART 0 — H3 AUDIT: THE 1.024 Å TORSION-AVERAGING LOSS IS NOT A BRANCH-CUT BUG

**THE DIRECTION STAYS CLOSED.**

Source: `s14/avgspace.py` line 82-83 —

```python
phi = R.circ_mean(PHI, axis=0)
psi = R.circ_mean(PSI, axis=0)
```

`R` is `s14.retprior`; `s14/retprior.py:53-66`:

```python
def circ_mean(A, w=None, axis=0):
    A = np.asarray(A, float)
    c, s = np.cos(A).mean(axis), np.sin(A).mean(axis)
    return np.arctan2(s, c)
```

This is `atan2(mean sin, mean cos)` — a **genuine circular mean**, not `mean(angles)` across
the ±π branch cut. The averaging-space-alone contrast reproduces exactly from the persisted
artefact (`s14/results/avgspace.json`, contrast `averaging_space_alone`, `C_torsion_mean_build`
vs `B_coord_avg_raw`): **+1.0241 [+0.6984, +1.3647], n=126, 41W/85L** — the precise figure in
project memory, confirming this is the experiment in question, not a mislabelled or superseded
run.

**Per the BRIEF's own branch rule** ("If it was circular, the direction stays closed"), no
rerun is registered. **H3 is closed.** All this workstream's compute budget went to H5.

**One qualification, named but not chased under this sprint's budget:** `circ_mean` is
unweighted and not basin-conditioned. A circular mean of a genuinely bimodal per-residue
distribution lands between the modes — exactly the failure mode the BRIEF itself warned about.
The audit rules out a *bug*; it does not certify that basin-conditioned averaging (cluster the
top-m torsions per residue into modes first, average within the dominant mode, per BRIEF §H3)
is worthless. That is a real, still-open question, left explicitly named rather than silently
dropped, and deprioritised behind H5 because H5 is the BRIEF's own "high-priority, low-risk"
item and this workstream owns the sprint's one serialised AMBER lane.

---

## PART 1 — H5: Ca-PRESERVING GEOMETRY REPAIR

**AMBER/OpenMM lane: OPENED** at the start of this experiment (`s23/agentB_h5_carestraint.py`,
`--mode sweep`), **sole authorised lane this sprint** per the BRIEF. Released when the sweep
below completes or is stopped; see the completion flag in
`s23/results/agentB_h5_carestraint.json` (`"complete": true/false`, `n_rows`/`n_expected`).

### Mechanism confirmed from source before running (see PREREG_B.md for the full derivation)

`core.amber.RESTRAINED_BACKBONE = ("N", "CA", "C")` — a **module-level global**, read once at
`AmberHamiltonian` construction (`_index_topology`), **not** a per-call argument to
`refine_coords`. The deployed/incumbent restrained minimisation (`core/pipeline.py:relax`,
default `amber_k=10` = `K_MODERATE`, `amber_steps=0`; the BRIEF's cited `E∘Relax_50` uses
`k_rest=100` = `K_STRONG`) restrains **N, Cα and C together** — the whole rigid peptide-plane
frame — never Cα alone. Potential: `0.5·k·|r-r0|²` per restrained atom, `d = F/k` at fixed
internal force F. **This experiment builds a genuinely separate Hamiltonian with
`RESTRAINED_BACKBONE` monkey-patched to `("CA",)`** (restrained-atom indices frozen at
construction, verified by an assertion that the CA-only index count equals the residue count
and the full-scope count equals 3×residue count) and sweeps the same k-ladder on both scopes.

### Basis

Point cloud throughout: `top75` coordinate average (the incumbent's own 3.048 Å arm),
rebuilt to a full backbone via `s15.phys_repl.averaged_backbone_from` so AMBER has bonded atoms
to act on. **Never compared to a rebuilt/projected chain's RMSD in this file** — where the
ideal-geometry projection appears it is its own separately-labelled row. Native Cα-RMSD is
read for **evaluation only** (ORACLE) and never used to choose a setting.

### Targets

`s16.repair.shape_subsample(30)` reused verbatim — the project's own pre-registered,
native-free, fold×length-stratified n=30 subsample, same seed. **This is a shape/trade-off
characterisation, not a 126-target claim on the BRIEF's primary endpoint.**

### Sweep

k ∈ {1, 3, 10, 30, 100, 300, 1000, 3000} kcal/mol/Å² × scope ∈ {Cα-only, N+Cα+C}, plus a
single k=0 (unrestrained) row common to both scopes. `steps=0` (minimise to convergence),
`tolerance=1.0`, CPU/Threads=1/DeterministicForces (matches `s16`'s canon and the pipeline
default depth).

### Results

**COMPLETE: 30/30 targets, 0 errors.** Artefacts: `s23/results/agentB_h5_carestraint_sweep.json`
(`"complete": true`, full per-target/per-setting rows), `s23/results/agentB_h5_carestraint_sweep_stats.json`
(every paired comparison below, iid + fold-clustered CI, MDE, W/L, worst-target delta).

**Basis note, stated once for the whole section:** every row below is a **point cloud**
(the incumbent's own `top75` coordinate average, rebuilt to a full backbone). The incumbent
row (no AMBER at all) on this n=30 subsample is **3.493 Å** — NOT the BRIEF's 126-target 3.048 Å;
this is a smaller, differently-composed sample, reported for internal comparison only, never as
a competing headline number.

#### Means, n=30, point-cloud basis throughout

| scope | k | Cα-RMSD (ORACLE) | Cα-disp | clash | cis frac | Rama outlier | ω-dev (°) | bond-len dev | bond-angle dev | converged |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| — | incumbent (no AMBER) | 3.493 | — | — | — | — | — | — | — | — |
| — | k=0 (unrestrained) | 3.940 | — | 0.00 | 0.099 | 0.161 | 22.4 | 0.0113 | 0.0246 | 30/30 |
| Cα | 1 | 3.779 | 1.247 | 0.00 | 0.100 | 0.172 | 23.0 | 0.0107 | 0.0240 | 28/30 |
| Cα | 3 | 3.774 | 1.099 | 0.00 | 0.099 | 0.160 | 23.5 | 0.0097 | 0.0228 | 28/30 |
| Cα | 10 | 3.718 | 0.954 | 0.00 | 0.125 | 0.191 | 29.8 | 0.0069 | 0.0205 | 28/30 |
| Cα | 30 | 3.641 | 0.784 | 0.00 | 0.224 | 0.180 | 49.7 | 0.0062 | 0.0188 | 28/30 |
| Cα | 100 | 3.570 | 0.507 | 0.00 | 0.514 | 0.301 | 92.0 | 0.0088 | 0.0234 | 28/30 |
| Cα | 300 | 3.535 | 0.332 | 0.00 | 0.599 | 0.353 | 99.0 | 0.0129 | 0.0334 | 28/30 |
| Cα | 1000 | 3.517 | 0.216 | 0.00 | 0.631 | 0.421 | 100.8 | 0.0161 | 0.0433 | 24/30 |
| Cα | 3000 | 3.512 | 0.163 | 0.00 | 0.637 | 0.463 | 103.6 | 0.0190 | 0.0478 | 22/30 |
| N+Cα+C | 1 | 3.759 | 1.124 | 0.00 | 0.098 | 0.161 | 23.2 | 0.0100 | 0.0231 | 28/30 |
| N+Cα+C | 3 | 3.717 | 1.022 | 0.00 | 0.095 | 0.181 | 23.4 | 0.0080 | 0.0204 | 28/30 |
| N+Cα+C | 10 | 3.671 | 0.918 | 0.00 | 0.098 | 0.196 | 24.5 | 0.0059 | 0.0168 | 28/30 |
| N+Cα+C | 30 | 3.631 | 0.827 | 0.00 | 0.103 | 0.178 | 26.4 | 0.0156 | 0.0168 | 28/30 |
| N+Cα+C | 100 | 3.585 | 0.683 | 0.00 | 0.093 | 0.174 | 26.9 | 0.0439 | 0.0366 | 28/30 |
| N+Cα+C | 300 | 3.541 | 0.507 | 0.00 | 0.049 | 0.168 | 24.2 | 0.0942 | 0.0551 | 22/30 |
| N+Cα+C | 1000 | 3.515 | 0.317 | 0.00 | 0.003 | 0.204 | 17.0 | 0.1689 | 0.0515 | 11/30 |
| N+Cα+C | 3000 | 3.507 | 0.206 | 0.00 | 0.000 | 0.218 | 10.2 | 0.2352 | 0.0728 | 10/30 |

(`N+Cα+C` at the pipeline's own defaults — `K_MODERATE=10`, `K_STRONG=100` — is the incumbent's
deployed restraint scope; `Cα` is this experiment's proposal.)

#### The finding: A REAL, MECHANISTICALLY CLEAR TRADE-OFF — NOT A CLEAN WIN EITHER WAY

**The falsifier does NOT fire in its strict form**, but not because Cα-only wins outright: it
survives because **the two scopes are Pareto-incomparable**, each dominating the other on a
disjoint subset of axes. Every number below is `paired_stats` (iid + fold-clustered CI, `n=30`,
`s23.qc_lib`), verdict per its own MDE:

**1. Bonded geometry (bond-length, bond-angle deviation): Cα-only wins, and it is not close at
high k.** Crossover at k≈10-30 (below it, full is marginally better; both are small regardless).
Above k=30, Cα-only is **MEASURED** better at every k, growing to **UNANIMOUS**:

    ca − full, bond_len_rel_dev,  k=1000:  -0.1528  SE 0.0183  MDE 0.0512  iid[-0.187,-0.117]  fold[-0.190,-0.115]  30W/0L
    ca − full, bond_len_rel_dev,  k=3000:  -0.2162  SE 0.0258  MDE 0.0723  iid[-0.264,-0.165]  fold[-0.269,-0.164]  30W/0L
    ca − full, bond_angle_rel_dev, k=100:  -0.0133  SE 0.0039  MDE 0.0108  iid[-0.021,-0.006]  fold[-0.022,-0.004]  18W/12L

At `k=1000`, restraining N and C alongside Cα leaves the incumbent's own restraint scope with
bond lengths **~17% relative deviation from ideal** (`0.1689`, e.g. peptide bonds pinned near
whatever the raw point cloud's imperfect spacing was) against Cα-only's **1.6%** — mechanistically
exactly the prediction in `PREREG_B.md`: freeing N and C lets the peptide-plane bond lengths
relax to their ff14SB optimum even while Cα stays pinned.

**2. Peptide-bond planarity (ω-deviation, cis fraction) and Ramachandran outliers: N+Cα+C wins,
and it is not close at high k either.** Full-backbone restraint DOMINATES from k≈10 upward:

    ca − full, cis_frac,       k=100:  +0.4217  SE 0.0608  MDE 0.1704  iid[+0.304,+0.536]  fold[+0.302,+0.541]   0W/20L
    ca − full, cis_frac,       k=1000: +0.6287  SE 0.0802  MDE 0.2247  iid[+0.471,+0.776]  fold[+0.468,+0.794]   0W/22L
    ca − full, omega_dev_mean, k=1000: +83.80   SE 9.83    MDE 27.54   iid[+64.4,+101.9]   fold[+62.4,+105.2]    0W/30L
    ca − full, rama_outlier,   k=1000: +0.2166  SE 0.0451  MDE 0.1262  iid[+0.135,+0.309]  fold[+0.148,+0.284]   0W/18L

At `k≥100`, Cα-only's own escape valve is exactly what the BRIEF's mechanism section predicted
but in the opposite piece of geometry from bond length: with N and C free, the minimiser can
satisfy the hard Cα position restraint (and whatever local bond-length/angle strain the input
carries) by **flipping the peptide bond toward cis** rather than translating Cα — cis fraction
climbs to **64%** of peptide bonds at `k=3000` (`0W/22L` against full's near-zero), which is not
a subtle validity cost, it is wholesale loss of one of the most basic all-atom sanity checks.

**3. Cα displacement: a genuine, significant CROSSOVER at k≈30**, not a win for either scope
throughout:

    ca − full, ca_displacement, k=1:    +0.1229  SE 0.0141  MDE 0.0396  iid[+0.096,+0.152]  fold[+0.109,+0.139]   1W/29L  (Cα-only moves MORE)
    ca − full, ca_displacement, k=10:   +0.0351  SE 0.0082  MDE 0.0230  iid[+0.020,+0.051]  fold[+0.021,+0.049]   5W/25L  (Cα-only moves MORE)
    ca − full, ca_displacement, k=30:   -0.0427  SE 0.0110  MDE 0.0307  iid[-0.065,-0.022]  fold[-0.060,-0.025]  19W/11L  (Cα-only moves LESS)
    ca − full, ca_displacement, k=100:  -0.1755  SE 0.0245  MDE 0.0687  iid[-0.221,-0.128]  fold[-0.230,-0.121]  23W/7L   (Cα-only moves LESS)

**Mechanism for the low-k regime (counter-intuitive, worth stating plainly):** at weak k, pinning
the WHOLE rigid peptide-plane frame (N, Cα, C together) leaves Cα almost immobile purely from the
local bond-length/angle rigidity of ff14SB — N and C being anchored geometrically over-determines
Cα's own position even with no restraint directly on it. Restraining Cα ALONE at the same weak k
leaves N and C free to swing on their bonded terms and drag Cα with them, so **Cα actually moves
MORE under a dedicated-but-weak Cα restraint than under a stronger-in-aggregate whole-frame
restraint at the same nominal k.** Only once k is strong enough (≈30+) that the direct positional
term dominates any bonded drag does Cα-only deliver the (obviously expected) lower displacement.

**4. Clash count: no differentiation at all.** `n_clash_heavy_2A` is exactly 0.00 for every
scope/k cell, `0W/0L` at every comparison — this basis (real retrieved fragments, coordinate-
averaged) essentially never produces a 2 Å heavy-atom clash after ANY restrained minimisation,
restrained or not. **A real null, reported as one**, not evidence either scope "fixes clashes" —
there was nothing to fix on this pool.

**5. Cα-RMSD-to-native (ORACLE, evaluation only, never used to pick k): small, mostly NOT
significant, sign inconsistent across k.** `ca − full` is MEASURED worse for Cα-only only at
`k=3` (`+0.057`) and direction-established at `k=10`; every other k is a null (CI includes zero
on the fold-clustered draw). **Accuracy is a wash between the two scopes** — the BRIEF's own
framing that a validity-only win is a real result even at flat RMSD applies directly here.

**6. Convergence rate degrades for BOTH scopes at high k, and MUCH faster for the incumbent
scope.** `Cα`: 28/30 through k=300, dropping to 24/30 (k=1000), 22/30 (k=3000). `N+Cα+C`: 28/30
through k=100, then **22/30 (k=300), 11/30 (k=1000), 10/30 (k=3000)** — at the two highest
settings, **most targets fail the project's own `CONVERGE_MAX_KCAL=1000` gate** under the
incumbent's own restrained-atom scope. This is itself a validity-relevant finding independent of
the paired comparisons above (which use whatever the minimiser returned, converged or not, per
the BRIEF's "reported, never silently applied" convergence-gate convention) — an operational
reason NOT to run the incumbent's scope at very strong k, on top of its bond-length cost.

#### Cross-check against the established project fact, on a DIFFERENT basis (stated as such)

Every restrained-minimisation arm — **both scopes, every k, including the incumbent's own
`K_MODERATE=10`/`K_STRONG=100`** — is **significantly WORSE on Cα-RMSD than doing no AMBER
repair at all** on this basis:

    full_k10_minus_incumbent    +0.1784  SE 0.0362  MDE 0.1014  iid[+0.108,+0.249]  fold[+0.101,+0.270]  4W/26L  MEASURED
    full_k100_minus_incumbent   +0.0923  SE 0.0248  MDE 0.0695  iid[+0.045,+0.139]  fold[+0.049,+0.135]  6W/24L  MEASURED
    k0_minus_incumbent          +0.4468  SE 0.1104  MDE 0.3093  iid[+0.247,+0.672]  fold[+0.319,+0.592]  6W/24L  MEASURED

This **replicates, independently and on a different basis, the established fact** (`s21/LEDGER.md`
L8: the deployed full-backbone AMBER repair on the PROJECTED/ideal-geometry chain is
`+0.0207 [SE 0.0034], 6.09 SE` worse than no repair, n=126). Here the input is the **raw point
cloud** (not the ideal-geometry-projected chain the pipeline actually feeds `relax()`), a starting
structure with far more bonded strain to begin with, so a much larger correction
(`+0.09` to `+0.45` depending on k) is expected and is **not comparable in magnitude to S21's
number** — **stated explicitly, per basis discipline, rather than implying a bigger effect on the
same quantity.** The DIRECTION replicates cleanly; the MAGNITUDE does not transfer across bases.
Unrestrained (`k=0`) minimisation is the single worst arm measured in this file, consistent with
the restraint's whole purpose (preventing the fold from drifting during minimisation) mattering
regardless of which atoms it is anchored to.

### Ablation (BRIEF Hard Rule: "ablate anything that wins")

**Nothing here is promoted as a winning arm** — this is a characterised trade-off, not a
single-number win, so there is no arm to ablate into deployment. The one directional claim that
would need an ablation before ANY promotion is narrow and explicit: *"if a downstream step needs
bonded-geometry fidelity (bond length/angle) more than peptide-bond planarity, Cα-only restraint
at k≈100-300 is the better choice of the two scopes tested."* Per `PREREG_B.md`'s promotion
clause, that claim would need re-running on a second, disjoint stratified subsample before it
reaches the 126-target instrument; not done here (budget), named as the next step.

### AMBER/OpenMM lane: RELEASED

Opened at the start of the H5 sweep; released after `s23/results/agentB_h5_carestraint_sweep.json`
reached `"complete": true` (30/30 targets, 0 errors, ~70 minutes wall clock, single-threaded CPU
platform, DeterministicForces, per `PREREG_B.md`'s logging block). No other lane's AMBER/OpenMM
activity was observed in `s23/results/` during this run. No live inter-agent channel was available
to announce this in real time (no `ListAgents`/multi-agent session reachable from this session);
recorded here and in the accompanying chat message instead.
