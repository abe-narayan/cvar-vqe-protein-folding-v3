# SPRINT 25 — ENDGAME. PHASE I: RMSD. PHASE II: FREEZE, CLEAN, RESULTS LAB.

**PRIMARY ENDPOINT: mean full-chain Cα-RMSD, 126 cluster-disjoint dev targets, point-cloud basis.
Incumbent 3.0483 Å. TARGET < 3.0 Å, and as far below as the evidence honestly allows.**
**The 60-target benchmark is SEALED.** No inspection, no RMSD, no tuning, no model selection.

**PHASE ORDER IS BINDING: research → best architecture → freeze → cleanup → result generation →
renderer → final report.** Do not build the results website while the architecture is still moving.

---

## 0. THE FOUR PILLARS — NON-NEGOTIABLE IN THE FINAL PRODUCTION PATH

1. **Genuine CVaR-VQE as selector.** Real variational state preparation, a Hamiltonian-derived
   objective, a real probability distribution, a real CVaR tail that materially participates. A
   classical ranking implementation is permitted **as a control**, never as the production selector.
2. **Genuine `H_Legacy`** — the 11-term potential at `DEFAULT_WEIGHTS`, never fitted, never proxied.
3. **Genuine `H_AMBER`** — ff14SB/GBn2 through the real OpenMM path, with bit-exactness checks. It
   stays in the final system even though it is not the strongest selector.
4. **The distogram is first-class**, not a side feature. Its role is to be determined by experiment.

Final shape: sequence → prior/representation → candidates → Legacy/AMBER evaluation → CVaR-VQE →
CVaR tail → ensemble → structure. Component ORDER may change if an experiment justifies it.

---

## 1. WHAT THE COORDINATOR HAS ALREADY MEASURED THIS SPRINT — READ BEFORE PROPOSING ANYTHING

`s25/LEDGER.md` L1 and L2, both n=126, both complete, both pre-registered.

**L1 — the posterior is over-confident by 2×.** z_sd 1.9962 against nominal 1.0; 90% coverage is
actually 0.6159. Nearly centred (z_mean −0.052). Over-confidence peaks at mid separations. 24.1% of
pairs are multimodal. Signed error grows monotonically with separation, −0.048 → −0.589 Å.

**L2 — fixing that makes RMSD WORSE, and the mechanism is now known.** A calibration-fitted widening
drove held-out z_sd to 0.9834 — essentially perfect — and cost **+0.0538 Å**. Tempering changes
calibration threefold and the endpoint by 0.003 Å.

> **RETRACTED IN FULL BY L7 — DO NOT BUILD ON THE SENTENCES BELOW.** They read as: the score is a
> location-based ranker and width is inert. **The audit lane falsified that on three independent
> grounds and every RMSD magnitude in L2 was BELOW ITS OWN MDE.** See `s25/LEDGER.md` L7 before
> choosing a direction.

**WHAT ACTUALLY STANDS.** (i) L1's over-confidence, at 1.66–2.09 across robust estimators, with the
defect living in the TAILS rather than the core. (ii) L2's NULL: recalibrating the posterior does not
buy RMSD. (iii) The width/variance/confidence family is closed — **but for the right reason (no
response) rather than the retracted one (wrong channel).**

**AND THE SCOPING THAT REPLACES IT.** Unconfounded, **location and width are equally flat**: pure
width at f=3.0 is +0.0253 and pure translation at −0.30 Å is +0.0206, all eight arms NOT MEASURED with
every fold CI including zero. The response to location is **symmetric** about the shipped posterior
even though it is biased −0.41 Å against the natives — **moving the prior toward the truth is worth
nothing and is indistinguishable from moving it away.**

> **Reaching γ requires genuinely BETTER location — new information — not a re-reading of the location
> already present in the posterior.** Re-readings of the existing posterior are now measured out at
> these amplitudes.

---

## 2. DO NOT REDO — closed with evidence

| closed | evidence |
|---|---|
| Posterior width / variance / confidence calibration | s25 L2: calibration achieved exactly, RMSD +0.054 |
| Pair reweighting by prior confidence | s24 L6: nested CV +0.0314, 33W/93L |
| Candidate generation, any provenance | s24: union +0.0022, quality-matched cosine 0.9432, 4 generated unions null |
| From-scratch torsion generator | s24: sequence channel ~0, chain correlation 7.2% of max, corpus smaller than the library |
| Residual generator | s24 L7: −0.032 Å at four stacked oracles, 0.71× MDE |
| Functional lever — filter, partition AND score | s24 L10/L16: oracle w with full leakage worth 0.0148 Å |
| Deeper VQE / larger χ / optimiser roulette / encoding | s21–s22, 5,760 cells |
| Probability-weighted CVaR readout | s23 L8: worse at 3 of 4 temperatures |
| AMBER refinement or Cα-only repair | s23 L11: 17 of 17 settings at or worse than no repair |
| Clustering / diversity routing / global scale / weighted averaging | s23 L5, L7, L10; s24 L15-A |

**Anything here needs a materially different formulation, a stated mechanism, and a pre-registered
falsifier — not a retry.**

---

## 3. THE ONE OPEN QUESTION IN PHASE I

> **Can the posterior's LOCATION be improved from information already available — and does the
> endpoint respond?**

The prior ladder says the endpoint's sensitivity is −2.15 Å per unit γ at the origin and that
**γ = 0.0225 reaches 3.0 Å**. L2 says the ranking responds to location and not to width. So the
question is narrow, well-posed, and cheap to attack. Candidate mechanisms — each needs its own
falsifier:

- **per-separation location correction.** The signed error grows −0.048 → −0.589 Å with separation.
  **But s24 Workstream C measured that this offset is NOT a distogram defect** — real protein windows
  carry the same offset (+0.4419) against these natives, so it is a corpus/native scale mismatch that
  the candidates share, and shifting only the prior may break a cancellation that currently helps.
  **State that prediction before running.**
- **length conditioning, residue-identity conditioning, empirical cross-target location priors.**
- **exploiting the 24.1% of pairs that are multimodal** — the median of a bimodal posterior sits
  between the modes, which is the same failure the project has met three times (circular means,
  damped samplers, amplitude shrinkage: *the safe arm is the absurd one*).
- **a different risk functional** whose minimiser is not the median — but this changes the deployed
  score and must be justified, nested-CV'd, and controlled.
- **pairwise consistency / covariance** — the posterior treats pairs independently and real distance
  matrices are metrically constrained.

---

## 4. HARD RULES

- **MDE = 2.8016 × SE, per comparison.** Report SE and effect/MDE. **0.7–1.3× is the Type-M zone and
  is not a result.**
- **Any oracle formed as a minimum over K variants must be scored against the distribution of the
  MINIMUM**, with `share_accounted`, residual, and `k_eff`. **W/L cannot diagnose it** — a best-of-K
  arm wins nearly everywhere by construction. Prefer a **split-half transfer**, which nulls itself.
- **Fold-clustered CIs beside iid.** iid CIs on this instrument are anticonservative: each of the 5
  distogram fold models saw ~100 of the other 125 dev natives.
- **Basis discipline.** Point-cloud and built-chain RMSD are never compared; the gap is 0.156 Å of
  pure operator choice. State the basis on both sides.
- **ORACLE / ACHIEVABLE / PRODUCTION** labelled at every appearance. Never promote an oracle.
- **Nested CV for anything fitted.** In-sample fits are reported as the optimism, never as the result.
- **Rule 0 — six operator forks** in the module docstring, each NAMING THE ALTERNATIVE NOT TAKEN:
  functional, basis, readout, normalisation, null, **and THE LABEL**. Send the list to the coordinator
  BEFORE the run. Enumeration by someone with a stake is a declared weakening.
- **Normalisation between Hamiltonians must be stated mathematically.** Raw AMBER on unrelaxed windows
  is clash-dominated — 57.5% above 1e4 kcal, worst 7.1e18 — which puts 99.4% of a pool inside
  |z|<0.1. Use rank standardisation, as the pipeline already does.
- **Provenance stamping is mandatory** for anything promoted:
  `ST.save_atomic(path, obj, complete_keys=NEED, rows=rows, n_expected=126, module_file=__file__)`.
  `complete` requires the FULL KEY SET, never a row count.
- **Never edit a module while a job launched from it is still running.** An artefact was written from
  vanished source once already; a `pkill` failed silently.

## 5. PHASE II — only after the architecture is frozen
Standardised 7-configuration comparison suite (Legacy, AMBER, Distogram, and the four combinations),
**same targets, same folds, same metric, same basis, CVaR-VQE as the matched selector throughout**;
≥50 targets with structures and renders per configuration, preferably all 126; deterministic
`T001__method.pdb` naming; `results/summary/{results.json,results.csv,leaderboard.csv,final_report.md,
professor_brief.md}`; `ARCHITECTURE.md`; a 3D renderer with native overlay and RMSD shown
prominently; production/research/archive separation; tests and build clean.

## 6. THE STANDARD
If <3.0 Å cannot be reached honestly, **do not manufacture it.** Report the closest defensible result
and state exactly what remains limiting. Never turn an oracle into an achievement, never a fluctuation
into a discovery, never hide a negative, never change the evaluation basis to improve a number.
