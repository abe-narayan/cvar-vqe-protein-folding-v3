# PREREG — LANE Q (CVaR-VQE / QUANTUM), SPRINT 25

Filed before any directional run. Point-cloud vs built-chain never compared.
MDE = 2.8016 x SE, per comparison. 0.7-1.3x is the Type-M zone and is not a result.

---

## Q-A — PROPERTY VERIFICATION (NON-DIRECTIONAL, NO FORK LIST)

Re-verification of the four standing facts I inherited, from source and by measurement.
These read no RMSD and have no outcome that can favour a hypothesis, so Rule 0's fork
enumeration does not apply (the s24 Lane-D precedent for `d_setequality_proof.py`).
**Falsifiers are registered instead.**

| # | claim under test | FALSIFIER — what would make me report it wrong |
|---|---|---|
| A1 | `MPSAnsatz` bond dimension is `2**layers` by construction; no SVD, no truncation | any state where the padded MPS's amplitudes differ from an independent dense simulation of the same gate list by more than 1e-12; or any `chi` not equal to `2**layers` |
| A2 | `run_cvar_vqe` is Adam on the EXACT parameter-shift gradient over a real statevector | parameter-shift gradient disagreeing with central finite differences on the same exact objective beyond FD truncation error; or a code path that samples |
| A3 | realised CVaR tail support is a SUBSET of an initial prefix of the energy order, equal to it under full support | one subset violation over my own adversarial families, built independently of `d_setequality_proof.py` and containing EXACT zeros (not epsilon floors — that is the s24 REV1 dead-assertion failure) |
| A4 | the harness's `paired_stats` MDE fix is in, and has no siblings | the 0.7x MDE rule absent; or any other place where a missing/unmeasured quantity is defaulted to a PASS |

## Q-C — HARNESS AUDIT (NON-DIRECTIONAL). Same status, same falsifier discipline.

---

## Q-B — DIRECTIONAL. **WHAT DOES THE CVaR TAIL ACTUALLY CONTRIBUTE AT THE DEPLOYED SETTING?**

Re-analysis of the stored artefact `s8/integrate_vqe.json` (n=126 targets, 5 pinned folds,
n_qubits=7, layers=3, iters=50, 3 alphas x 3 temperatures, full per-target arrays). No new
pipeline run; the data are fixed and already on disk.

**DISCLOSURE, MADE BEFORE THE FORKS BECAUSE IT WEAKENS THEM.** I read the artefact's stored
per-arm MARGINAL MEANS before filing this. They are printed in the file. What I have NOT
computed and am pre-registering here is the PAIRED contrast, its SE, its MDE, its
fold-clustered CI and its verdict. Enumeration by someone who has seen the marginals is a
declared weakening and I declare it.

**THE QUESTION.** The inherited statement is *"CVaR earns its keep as a REGULARISER at
readout: at T=0.1, alpha=1 collapses to 0.076 bits and degenerates to a plain argmin at
3.4540 A; alpha=0.1 holds 6.36 bits and reaches 3.3414 A — +0.113 A of genuine
contribution."* That comparison holds T fixed at 0.1. **The DEPLOYED table `VQE_LFO` in
`core/pipeline.py:118` runs T=0.3 on all five folds and alpha=1.0 on three of them.** So the
question is whether the alpha effect survives at the temperature actually shipped, or whether
alpha and T are substitutes and the +0.113 A is an artefact of the T=0.1 slice.

**PREDICTION, STATED BEFORE THE PAIRED STATS.** I expect the alpha effect to be
temperature-dependent and to be much smaller or absent at T=0.3, because both knobs act on
the same quantity — how concentrated p_theta is at readout — and the entropy term is the more
direct control of it. If so, the honest attribution is to *non-collapse*, not to *the tail*.

### RULE 0 — SIX OPERATOR FORKS, EACH NAMING THE ALTERNATIVE NOT TAKEN

1. **FUNCTIONAL.** The contrast is the per-target CA-RMSD of the arm's returned structure,
   differenced pairwise between two grid cells.
   *NOT TAKEN:* the CVaR objective value, the free energy, or the state entropy as the
   endpoint. Rejected because the objective is not the endpoint, and the whole reason this
   question is live is that an objective-level story ("the tail prevents collapse") was read
   as an endpoint-level contribution.
2. **BASIS.** Every number in Q-B is INTERNAL to `s8/integrate_vqe.json`; all its arms share
   one instrument and one basis (single-candidate selection out of a 128-candidate filtered
   set, not the top-75 coordinate average).
   *NOT TAKEN:* comparing any Q-B number to the s25 incumbent 3.0483 A. Different operator,
   different readout. **No Q-B figure will be quoted against 3.0483.**
3. **READOUT.** The stored `p_theta`-weighted consensus medoid, exactly as the artefact ran it.
   *NOT TAKEN:* re-running with the uniform coordinate-average readout, or with a
   probability-weighted CVaR readout (s23 L8 closed the latter at 4 of 4 temperatures).
   Taking it would confound an alpha effect with a readout change.
4. **NORMALISATION.** Energies as stored: `zrank` of the shipped distogram score over the
   top-2^7 filtered candidates.
   *NOT TAKEN:* a raw or moment z-score of the score. Rejected on
   `pauli-spectrum-delta-spike-artefact` / s24 D §4 — a moment z-score of an unconditioned
   energy measures its worst outlier. Rank standardisation is monotone, so it changes no
   ordering and, by the set-equality theorem, no tail membership.
5. **NULL.** The artefact's own no-circuit controls: `argmin` (the collapse limit),
   `boltz_T*` (Boltzmann weights over the same energies, no circuit at all) and
   `topfrac_*` (a flat average over the classical top fraction). These are matched in the
   OPERATOR'S space — same candidates, same energies, same readout — per
   `control-must-match-the-operators-space`.
   *NOT TAKEN:* differencing against an initialisation mean or against a uniform-random
   draw. Banned by `concentration-is-wrong-when-discrimination-binds` and
   `zero-information-control-must-be-plausible`.
6. **THE LABEL.** The label is the **paired difference in mean CA-RMSD, n=126**, with SE,
   effect/MDE, iid CI beside a fold-clustered CI, W/L, median beside mean, and the verdict
   from the FIXED `d_harness.paired_stats` rule (NULL / TYPE-M 0.7-1.3x / UNDERPOWERED <0.7x
   / MEASURED).
   *NOT TAKEN:* labelling on the marginal means (which is exactly how "+0.113 A" was
   produced), on W/L, or on a CI alone. A 5-cluster bootstrap CI excluding zero is NOT a
   licence to say MEASURED — that was s24 Lane D's own recorded defect.

### FALSIFIER FOR Q-B

**If the alpha effect at T=0.3 and T=1.0 is in the same direction and of comparable
magnitude to the T=0.1 effect, my prediction is wrong, the inherited "+0.113 A is the CVaR
tail" attribution stands unchanged, and I report it as confirmed.**

### WHAT I WILL NOT DO WITH THE RESULT

No grid cell will be re-selected on this outcome. `VQE_LFO` is a leave-fold-out table and it
stays pinned. Q-B can only change what is *said* about the component in `s25/QUANTUM.md`; it
cannot change what the component *does*.

---

## Q-D — CONDITIONAL, MATCHED-HARNESS RUNS FOR UPSTREAM RMSD-LANE CHANGES

Not yet triggered. If a change from the RMSD lane reaches me as promising, it goes through
`s24/d_harness.run_both` unchanged — genuine CVaR-VQE arm and the size-matched classical
control on identical candidates and an identical energy vector — and **its own six forks are
filed to the coordinator before that run.** No Q-D run happens on this pre-registration.

---

## STANDING CONSTRAINTS ACKNOWLEDGED

* LOCK_TRAIN and LOCK_AMBER are NOT taken by this lane without telling the coordinator first.
  Q-A/Q-B/Q-C take neither.
* No subagents.
* Artefacts provenance-stamped with `s24.stats_lib.save_atomic`.
* The 60-target benchmark is sealed and is not touched.
