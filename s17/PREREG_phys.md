# SPRINT 17 — PHYSICS workstream pre-registration

Written before any Sprint-17 AMBER minimisation or Legacy statistic was drawn.
Every experiment below states hypothesis · expected outcome · strongest control · success
criterion · falsifier. `s17/phys_FINDINGS.md` says for each rule whether it fired.

Standing data rules inherited and binding here:
* **Convergence gate** — a restrained minimisation is CONVERGED iff its final potential energy,
  restraint switched off, is finite and ≤ 1000 kcal/mol (`core.amber.CONVERGE_MAX_KCAL`).
  Reported with its exclusion count and the excluded PDB ids, every time.
* **Rotated-lab-frame null** — reported with its **MAXIMUM**, not only its mean. Pre-declared band
  `|mean| ≤ 0.005 Å AND max ≤ 0.05 Å` on the converged subset (`s16/energy_lib.FRAME_TOL_*`),
  unchanged. Every converged full-depth Sprint-16 arm FAILED it (mean −0.00095, max 0.1373 Å).
* **The AMBER path's own floor** is ≈0.004 Å gated / up to 0.038 Å ungated. The 0.08 Å multi-start
  floor does **not** transfer to it (no RNG on the path).
* **No validity statistic is quoted alone.** A constant ideal α-helix scores Ramachandran 1.000 and
  zero clashes by construction and beats AMBER 0W/65L. Every validity claim is a **conjunction**
  with staying near the input, and carries the α-helix zero-information reference.
* **Legacy weights are `DEFAULT_WEIGHTS` and are never fitted.**
* No native information in any predictor. The sealed 60-target benchmark is not read.
* `stable_rng` for all randomness; target is the unit; paired fold-clustered CIs, medians and W/L.

---

## E1 — Cα-preserving AMBER repair as CONSTRAINED OPTIMISATION

**The problem posed properly.**

> minimise `E_AMBER(x)` subject to `‖CA_i(x) − CA_i(x₀)‖ ≤ ε` for every residue i,
> for a ladder of ε; and the dual, minimise Cα displacement subject to validity ≥ threshold.

Implemented as a per-atom **flat-bottom** external restraint on the Cα set,
`0.5·k_f·max(0, |x−x₀| − ε)²` with `k_f = 1000 kcal/mol/Å²`, which is the exact penalty form of
the ε-ball constraint; and at ε = 0 as a **hard** constraint (Cα particle masses set to zero, which
OpenMM's `LocalEnergyMinimizer` treats as frozen). N, C, O, CB, sidechains and hydrogens are
**free at every rung**; the incumbent's harmonic N/CA/C restraint at k = 30 is carried as the
comparison arm, not as a rung.

**THE TRAP I AM DECLARING BEFORE I RUN IT.** At ε = 0 the Cα coordinates are unchanged *exactly*,
so Cα-RMSD is *identically* that of doing nothing. The pre-registered falsifier handed to me —
"accuracy cost above 0.10 Å against doing nothing" — **cannot fire at that rung**: it is true by
construction, not by measurement. I record this rather than claim the tautology as a result. The
scientific content of the ε = 0 rung is therefore **entirely on the validity axis**, and the
accuracy question only becomes non-trivial for ε > 0. The falsifier is kept verbatim for the ε > 0
rungs, where it can fire, and the ε = 0 rung is judged on H1b below.

**H1a (accuracy, ε > 0 rungs).** Flat-bottom-ε AMBER costs < 0.10 Å of mean Cα-RMSD against
doing nothing, for every ε ≤ 0.5 Å.
*Expected*: true, and trivially so — §3.4 of `s16/repair_FINDINGS.md` shows the RMSD cost of an
operator whose displacement is orthogonal to the residual is fixed by its magnitude alone, and a
per-atom ε-ball caps that magnitude at ε.
*Control*: **do nothing** (zero-information, the operator's own input, 3.0498 Å) and a
**matched-magnitude random Cα displacement** at the same realised RMS displacement, 3 `stable_rng`
draws.
*Success*: mean |Δ| < 0.10 Å with the CI reported.
*Falsifier*: mean cost ≥ 0.10 Å at any ε ≤ 0.5 Å.

**H1b (validity at zero accuracy cost — the real question).** With Cα **hard-fixed**, restrained
ff14SB/GBn2 minimisation still reaches the k = 30 arm's steric validity: zero heavy-atom clashes
below 2.0 Å on ≥ 90% of targets, and Ramachandran-favoured ≥ 0.824 (the k = 30 arm's 0.874 less
0.05).
*Expected*: **I expect this to FAIL on the bond/angle axis and possibly on Ramachandran.** The
input is a coordinate average whose backbone is contracted (`memory: averaging contracts the
backbone 25.8%`); freezing Cα freezes the Cα–Cα distances at their contracted values, and no
choice of N/C/O/CB can then produce ideal peptide geometry. If it fails there while succeeding on
clashes, the frontier is real but partial.
*Control*: the constant ideal α-helix (which scores rama 1.000 / zero clashes by construction and
is therefore the reference every validity number must be quoted against), the do-nothing input, the
ideal-geometry projection, and the incumbent k = 30 arm.
*Success*: both criteria met.
*Falsifier*: ≥ 1 clash below 2.0 Å on more than 10% of targets, **or** Ramachandran-favoured below
the do-nothing input's 0.836.

**H1c (the frontier is worth having).** There exists an ε at which validity is statistically
indistinguishable from the unrestrained k = 0 arm's while the Cα cost is < 0.05 Å.
*Falsifier*: no ε achieves both; validity rises only as Cα cost rises, i.e. the frontier is a
straight trade with no free lunch.

**Reported at every rung**: Ramachandran (favoured / allowed / outlier, wrapped torsions),
clashes < 2.0 Å and < 2.6 Å, minimum heavy separation, bond strain, angle strain, rms relative
geometric deviation, cis-peptide fraction, ω deviation, final AMBER energy, gate exclusions,
wall clock, realised Cα RMS displacement, and Cα-RMSD to native (ORACLE).

---

## E2 — Legacy, decomposed: which components carry information, and what is *not* in the distance model

**H2a.** At least one Legacy component has in-band ranking skill over pool candidates that the
distance model does not already contain.
*Measured*: per-component AUROC / AUPRC for "this candidate is in the target's best decile", the
**false-negative rate on near-native (< 2.0 Å) candidates**, false-positive rate, calibration, and
— the actual question — **incremental value beyond distance features by nested ablation**
(distance features alone vs distance + Legacy component, leave-fold-out, target as unit).
*Control*: a random score of the same distribution; and the distance-only nested model, which is
the only baseline that can answer "what does Legacy add".
*Success*: a component whose nested-ablation ΔAUROC (or Δ in-band Spearman) has a fold-clustered CI
excluding zero.
*Falsifier*: every component's incremental value over the distance features has a CI containing
zero — Legacy is redundant given the distance model.
*Expected*: `torsion` and `steric` carry the information (they are the two with detection skill in
`s16` §3.2), and both are largely **redundant with a Ramachandran prior** rather than with the
distance model; the MJ contact term carries nothing.

## E3 — Legacy as a GATE, and the danger

**H3.** `all candidates → Legacy rejection of the worst f → distance selection` beats **matched
random rejection of the same count** on realized Cα-RMSD.
*Primary danger, stated first*: **Legacy may remove exactly the candidates that matter.** Sprint 16
already caught the Legacy filter improving the set mean 3.551 → 3.531 while destroying the set best
2.306 → 2.416.
*Measured*: realized selected RMSD, **and near-native recall through the gate** — the fraction of
targets whose sub-2.0 Å candidates all survive the gate, and the set best before/after.
*Control*: matched random rejection of the same count, 3 `stable_rng` draws; and the no-gate arm.
*Success*: realized RMSD improvement with a fold-clustered CI excluding zero **and** near-native
recall not below the matched random control's.
*Falsifier*: recall of sub-2.0 Å candidates through the gate is below the matched random control's
— then the gate is a failure even if the mean improves. Also falsified if realized RMSD vs matched
random has a CI containing zero (the Sprint-16 result, +0.0118 [−0.0242, +0.0492], predicts this).
*Expected*: **falsified.** The prior is that this reproduces Sprint 16.

## E4 — The decisive Legacy-vs-AMBER experiment on IDENTICAL candidates

On one candidate set per target, compute distance score, Legacy total and its eleven components,
AMBER single-point energy (~6 ms; shortlist only, and the shortlist's recall cost is reported),
AMBER relaxed result on a sub-shortlist, and the validity panel. Evaluate each as **selector**,
**gate**, **rank feature**, and **repair operator**.
*Success criterion*: a statement of the form "Legacy knows X that AMBER does not, and AMBER knows Y
that Legacy does not", each half supported by a contrast with a CI excluding zero against its own
matched random control.
*Falsifier*: neither energy beats its matched random control on any of the four roles — in which
case the answer is that they know the same nothing, which is itself the result.
*Shortlist honesty*: the shortlist is native-free (top-M by distance score). Its **recall cost** —
the fraction of targets whose full-universe best candidate is outside the shortlist, and the
ceiling the shortlist imposes — is reported before any AMBER number.

## E5 — Disagreement as information

**H5.** Where distance, Legacy and AMBER rank differently, the disagreement predicts when the
distance selector is unreliable.
*Features (all native-free, target-level)*: Kendall τ and Spearman between each pair of rankers on
the shortlist, top-k candidate overlap, entropy of the selectors' choices, rank of each selector's
pick under the others' scores.
*Measured*: correlation of each feature with the distance selector's realized error, and with its
**regret** (realized − pool best), leave-fold-out.
*Control*: a random target-level feature of the same distribution; and the strongest native-free
non-disagreement predictor already known (pool typicality / consensus), so that any skill is
*incremental*, not re-discovered.
*Success*: a fold-clustered CI on the leave-fold-out Spearman excluding zero, **and** incremental
over typicality.
*Falsifier*: disagreement features have no leave-fold-out skill, or none beyond typicality.
*Expected*: weak positive at best. `memory: nothing ranks within the pool`, and disagreement is a
second-order functional of rankers that individually have no in-band skill.
