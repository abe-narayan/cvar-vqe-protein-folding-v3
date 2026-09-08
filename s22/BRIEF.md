# SPRINT 22 — FINAL CAMPAIGN. BINDING CONTRACT.

## 0. THE THREE PILLARS — NON-NEGOTIABLE

1. **GENUINE CVaR-VQE IS THE SELECTOR.** Not a classical argmin, not a heuristic, not
   "quantum-inspired". Classical methods are CONTROLS.
2. **`H_Legacy` GENUINE and INDEPENDENTLY EVALUABLE.** Never absorbed into a learned score.
3. **`H_AMBER` GENUINE and INDEPENDENTLY EVALUABLE.** Legacy and AMBER separately measurable at all
   times.

**PRIMARY ENDPOINT: mean full-chain Cα-RMSD, 126-target instrument.** Lower energy / loss / CVaR /
convergence is NOT an improvement. **Incumbent 3.048 Å** (point cloud, top-75 coordinate average) /
**3.204 Å** (built chain). Targets: **<2.5 Å**, aggressively **<2.0 Å**.

**The 60-target benchmark is SEALED.** Do not inspect, probe, tune against, or derive statistics from
it. Hash-verify at start and end.

---

## 1. WHAT SPRINT 21 CLOSED — DO NOT REDO

| closed | evidence |
|---|---|
| **Search** | exhaustive enumeration of the whole `2**n` latent; exact argmin **ties** a zero-eval pool (+0.085); ORACLE over the identical set **1.72 Å better, 0/122 reversals**; 512 random draws reach the global argmin (−0.045) |
| **Ansatz / χ** | χ=1→48; on AMBER `best_of_N` **wins the ladder**; χ orders nothing |
| **Optimiser** | no optimiser beats `best_of_N` on RMSD; training the circuit +0.003 |
| **Encoding θ vs (sin,cos)** | 5,760 cells, 5 independent legs; **radius ≡ scale, the same knob** |
| **Training Hamiltonian** | +0.056 [−0.020,+0.129], null |
| **Generative latent as a SOURCE** | −0.103 Å at the ORACLE **of the union**; closed by arithmetic |
| **Consensus selector family** | 9 arms, **0%** of the in-pool gap |
| **Legacy→AMBER continuation** | staged is WORSE (+0.233); a raw λ ladder **never visits an intermediate H** |
| **AMBER preconditioning** | drives AMBER **31× deeper** into the steric singularity |
| **Hamiltonian-vs-Hamiltonian disagreement** | no per-target skill; whole ORACLE gain is min-of-k bias |

**Both physics energies are worse than a matched RANDOM tail at n=126 (42/42 CIs excluding zero).**
Distance alone beats all seven matrix cells.

---

## 2. WHAT IS OPEN, AND THE ONE MEASURED NUMBER THAT MATTERS

**THE ROUTING CEILING — the largest deployable lever identified in the project.**

    best single fixed arm (pool top-75 average)   3.048
    ORACLE per-target routing over 13 arms        2.566
    HEADROOM                                      0.482 A

**Every one of the 13 native-free arms wins on some target** (winner counts 1–17, no arm dominant).
The correct baseline for a router is **the best FIXED arm, 3.048** — a router must beat that.

> **A NULL I ALREADY GOT WRONG, RECORDED SO NO LANE REPEATS IT.** I tried a permutation null that
> shuffled targets *within* each arm. It returned **1.115**, i.e. *better* than the observed oracle.
> **It is mis-specified**: arms are strongly correlated because TARGET DIFFICULTY dominates, and
> shuffling destroys that shared factor, letting the min draw from easy targets. **Do not use it.**
> The oracle routing value is a legitimate ceiling (like `latent_oracle`); the min-of-K concern
> applies to *estimating a router's performance*, not to the ceiling itself. **Any router must be
> scored on held-out folds, never on the folds that selected it.**

**The other live lead: the compactness/disagreement channel (S21 L27/L28).** `rg_z`, `rg_gap` —
scale-free contrasts of the distogram's predicted Rg against the pool's realised Rg. Partial ρ ≈
**0.34–0.38**, clears a pooled 17-predictor permuted bar, survives length stratification (effect is
*larger* stratified), a native-free difficulty control, and the pool's own realised extension.
**DEMONSTRATED AS A SIGNAL. ÅNGSTRÖM VALUE NOT MEASURED — that is the campaign's central conversion.**

**Untested and mechanistically motivated:** error-placement weighting `H_D* = Σ s_ij w_ij (Δd_ij)²`
with native-free `s_ij` (Sprints 18–19 established equal residual magnitude has wildly different
structural consequence depending on placement). And candidate-state quantum encoding `H|i⟩ = E_i|i⟩`.

---

## 3. HARD RULES

1. **RMSD-first.** Every claim ends in Cα-RMSD or is labelled as not doing so.
2. **ORACLE vs ACHIEVED vs FLOOR** labelled at every appearance. ORACLE quantities are ceilings.
3. **MDE = 2.8016 × SE, per comparison.** Never a constant. Report SE beside every mean. Report
   effects as a **multiple of their own MDE**; 0.7–1.3× is the Type-M zone — directions safe,
   magnitudes upper bounds.
4. **BASIS.** Point-cloud and built-chain RMSD are never compared. State the basis on both sides.
5. **Completion flags demand the FULL KEY SET**, not a row count. A flag that never fires is vacuous.
6. **Atomic writes AND a per-run output path.** Atomicity protects readers from a half-written file;
   it does **nothing** against a second writer. Derive the filename from the config.
7. **`pkill -f` succeeding is not evidence a process died.** Check the thing, not the wrapper.

## 4. RULE 0 — OPERATOR FORKS (the rule Sprint 21 paid for)

For any comparison with a **directional hypothesis**, enumerate **six** axes in the module docstring
and **name the alternative not taken**: **functional, basis, readout, normalisation, null, and THE
LABEL.**

- *Why:* one primary carried three unstated forks worth **+0.215 Å against a true effect of +0.025**
  — 89% operator, all three pointing the author's way.
- **Clause 2:** where one person both designs a comparison and has a stake in its direction, the forks
  are enumerated by someone who does not. **Send the fork list BEFORE the run, not the result after.**
  A reviewer acquires a stake the moment they write a verdict — a fork review is itself an artefact.
- **Clause 3:** a best-of-K result's null is the distribution of the **maximum**. Best-of-10 |AUC| =
  0.625 is the *median* of pure noise at n=75.
- **THE LABEL:** a binarised label whose threshold depends on a covariate **manufactures a predictor
  of that covariate**, and will pass a correct permutation null.

## 5. PRE-REGISTRATION — mandatory before any large experiment

Hypothesis · primary endpoint · expected mechanism · **falsifier** · null · matched control · budget ·
promotion criterion. **Falsifiers are never moved afterward.** A pre-registration that names a
falsifier buys you the apparatus that refutes you — that is its value, not the honesty.

## 6. READ THE BODY, NOT THE INDEX LINE

Project memory is an index of summaries. Bodies carry scope limits, supersessions, the variable a
quantity is actually indexed by, and open directions. **Two lanes cost themselves a confound and a
retraction in one night by quoting index lines** — and the same unread section held the sprint's only
positive result.

## 7. CONTROLS

- A control must be matched in the space the **operator** works in.
- Zero-information controls must be **plausible** (constant α-helix), never uniform-on-the-torus.
- Never an initialisation mean; use **best-of-N from the untrained circuit**.
- When an analytic null and a measured null disagree, **the measurement is the null**.
- **`min` and `mean` readouts consume different statistics of a set.** No statement about an
  objective's usefulness is well-posed until the readout is fixed.

## 8. COMPUTE

~95% CPU, never intentionally above 97%. **AMBER/OpenMM context is SERIALISED — one lane at a time,
announce before taking it.** The box fits two heavy jobs; 3+ breach RAM.
