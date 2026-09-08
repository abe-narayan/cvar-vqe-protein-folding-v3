# PRE-REGISTRATION — FEATURES workstream, Sprint 17

Written **before** any number in `s17/features_FINDINGS.md` was produced. Every experiment below
carries hypothesis · expected outcome · strongest control · success criterion · falsifier.

---

## 0. The problem this lane attacks, and why it is the last open one

Problem **B** (selection), specifically the **in-band** sub-problem localised by L4/L5/L10 and
SELECT §5.

The shipped distance objective has global Spearman 0.568 and **in-band Spearman 0.131**; its argmin
inside its own top-25 is **−0.014 [−0.118, +0.094], 62W/64L** against a random pick from that same
top-25. Three independent closures (58 distogram functionals; 29 native-free in-band signals;
27 target-level calibration features) all land at zero.

**Every signal closed so far is geometric or energetic** — computed from a candidate's own
coordinates or from a physical potential. **No learned sequence representation has ever been tested
in-band.** The record contains one overturned result claiming ESM buys **0.288 Å on SELECTION**
(p = 0.005, n = 126, `s7/repr_select.py`) — but that was measured **globally**, through a *retrained
distance predictor*, never in-band and never as a per-candidate feature.

Two readings are possible and they are the experiment:

* the 0.288 Å is **in-band content** → the most valuable thing in the programme;
* the 0.288 Å is **more garbage-filtering** (a better global filter, same in-band void) → the last
  open feature class closes and the sprint's selection answer is final.

---

## 1. THE BAR — stated before any measurement

The distance objective is **not** the bar. It has no in-band skill, so beating it is meaningless.
A representation-derived feature counts as having in-band skill only if, in a **native-free**
band (the shipped objective's own top-25 or top-75 at K = 500), it beats **all** of:

1. **random-in-band** (exact expectation = band mean), and
2. **`NULL_alpha`** — CA-RMSD of each candidate to a **constant ideal α-helix**, and
3. **`NULL_beta`** — CA-RMSD to a **constant ideal β-strand**,

on **selected CA-RMSD**, paired at the **target** level, with a **fold-clustered** interval
excluding zero, medians and W/L printed beside every mean.

A positive in-band Spearman whose *selected* column does not clear (1)–(3) is **not** a result.

## 2. UNITS — declared in advance (L3's failure mode)

Three axes appear and they are not interchangeable:

    pairwise ordering ACCURACY (null 0.500)  ->  Kendall tau = 2*acc - 1
    Spearman rho_S
    Gaussian-copula rho = sin(pi*tau/2)  =  2*sin(pi*rho_S/6) from Spearman

`s17/selection_theory.py` supplies `rho_from_accuracy`, `accuracy_from_rho`, `rho_from_spearman`.
Every table prints **the measured Spearman** and, beside it, **a directly-counted pairwise ordering
accuracy** (concordant pairs, not converted), so the axis is never inferred. Any requirement quoted
is quoted on both axes.

## 3. TRAP GUARDS — declared in advance

* **Min-of-N (L8, SELECT §6).** No "best of V features chosen with labels" number is quoted
  without the zero-information minimum over V random picks from the same band beside it. Every
  feature's **sign is fixed a priori** by its physical meaning, so the battery is not a hidden
  best-of-2V screen; where a sign must be chosen it is chosen **leave-fold-out**.
* **Tie trap (S8).** Selection is tie-safe via `sel_lib.sel_of`, which averages the outcome over
  the whole tied argmin set. A constant feature then scores exactly the band mean, as a
  zero-information signal must.
* **Boundary trap.** Any ladder-selected hyper-parameter landing on its boundary is reported as a
  substitution, not a parameter.
* **Leakage.** All training labels (native local geometry, native SS) are used **inside training
  folds only**; inference is native-free. The ESM caches read here are sequence-only — no structure
  of the scored target reaches them. ESM-2's **contact head is supervised on structure** (of *other*
  proteins, not the target) and is flagged wherever it is used.
* **Sealed benchmark.** Not read, not probed, not derived.

---

## F0 — INVENTORY (no hypothesis; an audit of what exists)

Establish what representation data exists, per target and per candidate window, at what cost, and
whether it is fold-honest. Reported as fact, not as a result.

---

## F1 — ESM contact head vs the candidate's realised geometry  [in-band]

**Hypothesis.** ESM-2's contact head is the only structure-supervised part of the model and is a
*different channel* from the distogram (attention-derived, not the shipped predictor). If any
learned representation discriminates in-band, this is the most likely one.

**Expected.** Small positive in-band Spearman; selected RMSD better than band-random.

**Controls.** `NULL_alpha`, `NULL_beta`, `NULL_chance`, random-in-band; plus `dist_shipped` as a
*reference* (not a bar).

**Success.** Clears §1's bar in top-25 or top-75.

**Falsifier.** In-band Spearman CI contains zero AND selected-vs-band-random CI contains zero in
both native-free bands → the ESM contact channel carries no in-band information.

## F2 — ESM contact CONFIDENCE as a per-pair WEIGHT on the shipped residual  [in-band]

**Hypothesis.** L5 closed 58 functionals of the distogram, including its *own* sd-based weights.
An ESM-contact-derived weight is an **exogenous** weighting — information from outside the
distogram deciding which pairs to trust.

**Expected.** Weak. L5's mechanism (you cannot fix an in-band problem by re-weighting a global
signal) predicts this fails, but the weight is exogenous, which is the one thing L5 did not test.

**Control.** The shipped weighting itself, and a **uniform** weighting, on the same instrument.

**Success / falsifier.** As §1.

## F3 — Sequence-conditioned pseudo-likelihood of the candidate window's OWN source sequence

**Hypothesis.** A retrieved window is a real fragment with its own sequence. ESM-2 gives a
per-position distribution over amino acids **conditioned on the target sequence**; scoring a
candidate by how likely its source sequence is under the *target's* context is a learned-sequence
measure of "does this fragment belong here" that BLOSUM (a position-independent substitution
matrix, measured null) cannot express.

**Expected.** This is the cleanest new channel. Prior is weak — `structure-and-sequence-are-decoupled`
records 12% identity between a target and its structurally nearest window.

**Controls.** (i) `blosum_sim` — the position-independent substitution score, already measured null;
(ii) a **background-frequency** pseudo-likelihood with no sequence context at all (zero-information);
(iii) §1's constant conformations.

**Success.** Clears §1's bar AND beats the BLOSUM control with a fold-clustered CI excluding zero
(otherwise it is BLOSUM re-derived, not a representation result).

**Falsifier.** Does not separate from the background-frequency control.

## F4 — Does the embedding predict THIS residue's local environment, and does agreement discriminate?

**Hypothesis.** Per-residue ESM embeddings predict local backbone environment; a candidate whose
realised local environment agrees with the prediction is more likely native-like.

**Model.** Leave-fold-out ridge from the target's per-residue ESM-PCA32 embedding (plus its
neighbours) to the residue's native local CA descriptor `[d(r,r+2), d(r,r+3), d(r,r+4)]`. Labels are
native structures **of training-fold targets only**. Candidate score = mean standardised
disagreement between prediction and the candidate's realised descriptor.

**THE CONTROL THAT DECIDES IT.** The identical model with a **one-hot** input block. This is the
exact contrast the 0.288 Å claim rests on. A third arm predicts the training-set **mean** descriptor
(no sequence at all) — the zero-information reference for this family.

**Success.** ESM arm clears §1's bar AND beats the one-hot arm with a fold-clustered CI excluding
zero.

**Falsifier.** ESM arm indistinguishable from one-hot in-band → the recorded 0.288 Å is **not**
in-band content.

## F5 — Predicted secondary structure vs the candidate's realised secondary structure

**Hypothesis.** As F4, coarsened to a 3-class per-residue SS label. Lower capacity, so more
learnable from 126 targets.

**Model / controls / criteria.** Identical to F4 (ESM vs one-hot vs marginal), with a multinomial
logistic head and a per-residue log-likelihood score. SS is assigned from **CA coordinates only**
by one fixed rule applied identically to native and candidate, so the two are commensurable.

**Falsifier.** As F4.

## F6 — Embedding-space agreement between the target sequence and the window's SOURCE sequence

**Hypothesis.** The retrieved window has its own sequence; embedding it and comparing to the
target's embedding is a learned similarity that BLOSUM cannot express.

**Cost gate, declared before running.** Per-candidate embeddings **do not exist** in the repo
(F0). Computing ESM-2-650M for the full universe (13k–27k windows × 126 targets ≈ 1.6M sequences)
is out of budget on this box. Computing it for the **band only** (top-75 of K = 500 → ≤ 9,450
sequences) is feasible. **F6 runs only if F1–F5 leave the lane open and the box has headroom**;
otherwise it is reported as NOT RUN with its cost, which is an honest outcome and not a silent gap.

**Controls.** `blosum_sim`, F3's background control, §1's constant conformations.

**MODEL LADDER, declared before running (added when the box's contention made the cost gate
bind).** `esm2_t6_8M_UR50D` (30 MB, ~80× cheaper, already on disk) is run **first, as a declared
scout**, at B = 25. `esm2_t33_650M_UR50D` — the shipped model — is run **only if the scout shows
in-band signal**. Two model sizes are **two different instruments**: their numbers are never
pooled, the model that produced each number is recorded in the artefact and quoted in every table,
and a null from the scout is reported as *a null from the 8M model*, not as a null from ESM-2.
This ladder is a compute decision, not a scientific one, and is written down before either run so
it cannot be reconstructed afterwards as a choice made on results.

## F7 — ADDED MID-LANE (pre-registered before running, after F1–F5 returned)

**Why it was added.** F1 returned the largest in-band **rank correlation** any native-free signal
has shown in this sprint (ESM contact head, ρ_S = +0.116 [+0.047, +0.184] at top-75) while its
**argmin** stayed null. A signal that orders the band a little but cannot pick its best is useless
as a ranker and may still be useful as a **shortlist constructor** — a role L12/L14 identify as
unoccupied, and which the coordinator independently flagged (the shipped score's top-M shortlist
has an ORACLE best *significantly worse* than a matched random shortlist of the same size).

**Hypothesis.** A shortlist built from the ESM contact score retains near-native members that the
distance score's top-B discards.

**Expected.** Weak. L12's mechanism is a property of *any* score-ranked shortlist, so a second
score should inherit it.

**Controls, all three mandatory.** (i) the shipped distance top-B, (ii) a **matched random**
shortlist of the same size B, (iii) the union/interleave of the two score shortlists at the same
total size — otherwise a union arm is priced against a smaller shortlist and wins for free.

**Success.** ORACLE best inside the ESM-built shortlist beats the distance-built shortlist of the
**same size** with a fold-clustered interval excluding zero.

**Falsifier.** No improvement over the distance shortlist, or no improvement over matched random
— in which case the shortlist-construction role is closed for this feature class too.

**Declared in advance:** this arm reports a **CEILING** (the shortlist's ORACLE best). A ceiling
gain is worthless unless a realized readout can consume it, so the realized argmin and medoid
inside each shortlist are reported beside it, per L14.

---

## LANE-LEVEL FALSIFIER

> **If no representation-derived feature beats a constant ideal α-helix in-band with a
> fold-clustered interval excluding zero, the last open feature class is closed and the sprint's
> selection answer is final.**

If any feature *does* clear the bar, the immediate follow-up is **its ceiling, not a model**: feed
the measured in-band ρ through `s17/selection_theory.py`'s band curve and report the RMSD it would
deliver at each shortlist size, before anyone builds a predictor around it.
