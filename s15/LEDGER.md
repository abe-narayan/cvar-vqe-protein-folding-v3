# LEDGER

Every experiment run in Sprint 15, in the order it was run, with its verdict and its artefact.
**Append-only.** A row is never edited once written; a superseding row is added below and the
original is marked, so the sequence of belief is recoverable — which matters more here than the
final state, because three entries were retracted by the person who wrote them.

Legend: **✓** result stands · **✗** refuted · **⊘** retracted (the measurement itself was invalid)
· **↻** superseded by a later, better-powered run · **⋯** in flight

---

## PHASE 0 — AUDIT AND INSTRUMENT FREEZE

| # | experiment | module | artefact | verdict |
|---|---|---|---|---|
| 0.1 | rebuild `peptide_db.npz` from 1,524 PDB files | `audit_repro` | `audit_repro_db.json` | ✓ 787/787 bit-identical, max diff 0.0 |
| 0.2 | rebuild all 126 window universes | `audit_repro` | `audit_repro_univ.json` | ✓ bit-identical, 2.35 M windows |
| 0.3 | reproduce `pool_best` from PDB files, not cache | `audit_repro` | — | ✓ 1.7108244199364904 exactly |
| 0.4 | end-to-end pipeline reproduction | `audit_pipe` | `audit_pipe_run.json` | ✓ `sub`, `ca`, `fit_ca`, `amber_ca` all max-diff 0.0 |
| 0.5 | determinism across torch thread counts | `audit_determinism` | `audit_determinism_t{1,2,8}.json` | ✓ zero argmin flips; ✗ distogram not bit-stable (`prob` 2.0e-6, `risk` 1.3e-4) |
| 0.6 | independent RMSD reimplementation (Horn quaternion) | `audit_rmsd` | `audit_rmsd.json` | ✓ 2.04e-13 Å over 63,000 structures; 15/15 edge cases including the mirror test |
| 0.7 | AMBER ff14SB/GBn2 parameter verification | `audit_amber` | `audit_amber_recompute.json` | ✓ bit-exact against an independently built ForceField |
| 0.8 | AMBER cost, memoisation defeated | `audit_amber` | `audit_amber_cost_sweep.json` | ✗ **8.3–23.3 ms, not 28 ms** — budget matching over-charged AMBER 2–3× |
| 0.9 | Legacy total vs its 11 components | `audit_energy` | `audit_energy_legacy.json` | ✓ worst 2.7e-5 over 1.28e7 structures |
| 0.10 | CVaR vs Rockafellar–Uryasev | `audit_energy` | `audit_energy_cost.json` | ✓ max error 4.2e-14 over 400 cases |
| 0.11 | load-bearing constant sweep | `audit_sensitivity` | `audit_sensitivity.json` | ✗ **`BAND = 1.5 Å` has no derivation**; `n_zero_recall` 45→2 across 0.5→3.0 |
| 0.12 | oracle-conditioning of the cached AMBER subset | `audit_amber` | `audit_amber_oracle.json` | ⊘ the brief's "40% oracle-conditioned" is **transposed**: 21.7% oracle, 40% unbiased |
| 0.13 | the `amber_kind == 0` stratum in its tails | `audit_amber` | `audit_amber_snap_inband.json` | ✗ **not unbiased in tails** — `a_idx` force-includes the ORACLE snap index |
| 0.14 | leakage disclosure | `audit_provenance` | `audit_provenance.json` | ✓ 16/126 verbatim own-fold windows, 4 self-hits; **measured impact on `pool_best` exactly zero** |
| 0.15 | `s14/results/obj_floor.json` AMBER rows | `audit_energy` | — | ⊘ **RETRACTED** — no `amber_kind` mask, yet flagged complete |

**Gate verdict: PASS**, with one retraction, one correction to the brief, one new defect, four
disclosures.

---

## THE GENERATIVE ARCHITECTURE

| # | experiment | module | artefact | verdict |
|---|---|---|---|---|
| 1.1 | analytic torsion gradient vs central differences | `distgeo` | — | ✓ **cosine 1.000000000000000**; four terminal torsions inert (`|g|` ratio 9e-17) |
| 1.2 | torsion distance geometry, 8 targets | `distgeo` | `distgeo_smoke.json` | ↻ ORACLE 0.697, predicted 2.792 — **superseded by 1.3, which reverses the sign** |
| 1.3 | torsion distance geometry, all 126 | `distgeo` | `distgeo.json` | ✓ ORACLE **0.611** (−2.593 vs incumbent); predicted **3.644** (**+0.440 [+0.290, +0.592]**, a loss) |
| 1.4 | distogram error profile and sd calibration | `distacc` | `distacc.json` | ✓ MAE 2.386, **bias +0.509** rising to +1.492 at separation 11–15, **z-sd 2.633** |
| 1.5 | leave-fold-out separation debiasing | `distcal` | `distcal.json` | ✓ legitimate and native-free; helps the fit |
| 1.6 | maximum likelihood on the full 17-bin distribution | `distml` | `distml.json` | ⊘ **VOID** — `LogPTable` used a uniform-grid lookup on a non-uniform grid |
| 1.7 | the "1.06 Å quantisation ceiling" | `distml` | — | ⊘ **RETRACTED same day** — an artefact of 1.6, not a measurement |
| 1.8 | `LogPTable` correction + derivative re-verification | `distml` | — | ✓ `searchsorted` on real centres; **cosine 1.000000000000000**, max err 5e-09 |
| 1.9 | the retrieval pool as a second distance channel | `pooldist` | — | ✓ pool MAE **2.249** vs distogram 2.386, **opposite bias sign** at short/medium separation |
| 1.10 | multi-start seeding | `seed` | — | ⊘ **`hash(pdb)` is salted per process** — no multi-start constant was reproducible. Fixed with `blake2b`; verified identical across two interpreters |

---

## THE MECHANISM EXPERIMENTS

| # | experiment | module | artefact | verdict |
|---|---|---|---|---|
| 2.1 | where does the objective rank the native? | `feasible` | `feasible.json` | ✓ native at the **34.7th percentile**, argmin on 4/126 — reproduces the standing 36.8th/3-of-126 through independent code |
| 2.2 | is the objective's argmin better than random? | `feasible` | `feasible.json` | ✓ **−0.949 Å [−1.147, −0.753]**, ρ positive on 88.1% of targets — **the rebuttal to Roget et al.** |
| 2.3 | ORACLE control on 2.1/2.2 | `feasible` | `feasible.json` | ✓ native at percentile **0.000 on 126/126**, ρ +0.862 — machinery verified |
| 2.4 | the pool channel as an objective | `feasible` | `feasible.json` | ✗ **worst arm** — native at the 66.5th percentile, *above* chance. A better estimator and a worse objective |
| 2.5 | Family B feasibility | `feasible` | `feasible.json` | ✗ **falsified as declared** — at ε = 1σ only 45.8% of pairs admit the native; median z 1.336, max z **7.237** |
| 2.6 | the full G→S→A→F cascade, 126 targets | `cascade` | `cascade_combined.json` | ✗ **does not beat the incumbent** — F **+0.117 [+0.050, +0.188]**; A at parity (CI includes zero) |
| 2.7 | cascade decomposition | `cascade` | `cascade_combined.json` | ✓ G→S **+0.751**, S→A **−0.360**, A→F **+0.170**; filtering the ensemble by the objective makes aggregation *worse* |
| 2.8 | distance-accuracy phase diagram | `distacc` | `distacc.json` | ⋯ in flight — converts the programme into a stated accuracy requirement |
| 2.9 | robust and redescending losses | `robust` | `robust.json` | ⋯ queued; 8-target smoke shows the redescending signature (welsch: better median 2.466, worse mean 2.967) |
| 2.10 | graduated non-convexity | `robust` | `robust.json` | ⋯ queued — the standard cure for 2.9's failure mode |
| 2.11 | undoing the averaging contraction before projection | `expand` | `expand.json` | ⋯ in flight — attacks the A→F gap; would improve the **production** pipeline if it works |
| 2.12 | augmenting the retrieval pool with solved conformers | `augment` | `augment.json` | ⋯ queued — needs only the fits to beat the set mean, which holds by 0.43 Å |

---

## THE AGENT WORKSTREAMS

| # | workstream | artefact | headline |
|---|---|---|---|
| 3.1 | Phase 0 audit | `audit_FINDINGS.md` | ✓ gate PASS; rows 0.1–0.15 above |
| 3.2 | literature and novelty boundary | `LITERATURE.md` (961 lines, 67 papers) | ✗ two claims taken (Roget for lattice certified optima; QuPepFold for in-loop CVaR); ✓ the untrained-circuit control has no precedent |
| 3.3 | information channel audit | `info_FINDINGS.md` | ✓ the i.i.d. surface over-prices real channels — 0.6–1.7 Å native-free, up to 2.9 Å ORACLE; alignment 0.565/0.665 against a 0.945 null, **but a constant α-helix scores 0.709**; 2.0 Å **not supported**, 2.5 Å **marginal** |
| 3.4 | | | ✗ its own headline refuted — the native-free sign proxy does not transfer to real pools |
| 3.5 | quantum geometry and CVaR | `qgeom_FINDINGS.md` (1,197 lines) | ✓ classical Fisher = **4×** Fubini–Study to 3.3e-16; ✗ Sprint 14's condition numbers are a seed-averaging artefact; ✓ the gradient avoids the metric's small directions; ✗ **the budget convention decides the QNG conclusion** (5/5 at equal iterations, 0/7 at equal cost) |
| 3.6 | representation ladder | `restraint_FINDINGS.md` (675 lines) | ✓ **75% of my "quantisation ceiling" was my bug**; the real cost is +0.381 Å [+0.156, +0.614] and it is bin **placement**, not count; ✗ ML still does not beat least squares — the mechanism is **calibration** (ML wins on exactly the 52 calibrated targets, by nothing on the other 74); ✓ **an empirical false-positive floor of ~0.08 Å**, demonstrated on a provably-null comparison |
| 3.7 | the quantum experiment | `qrestraint_FINDINGS.md` (766 lines) | ✗ **the Sprint 14 negative does not reproduce** — VQE beats the untrained control at budget 8,192 (−0.202, 13/19), but the win is **objective-independent**, **absent at 32,768**, and **never beats greedy** (0/16); ✓ **the CVaR gradient defect HELPS** (2.4× smaller gradient norm = a step-size reduction); the binding quantity is a selection gap **identical in both arms** |
| 3.8 | alignment engineering | `align_FINDINGS.md` (760 lines) | ✓ **FIVE** exactly-null torsion directions, not four — the fifth a distributed whole-chain crank, derivable as `2n−4` parameters onto a `2n−5` shape space; ✓ ORACLE rotation is worth **−1.822 Å**; ✗ **nothing native-free reaches it** (best −0.007), so **alignment is a description, not a lever — Family D stays empty**; ✗ **the terminal-gap lever does not transfer** (perfect native termini worth −0.000 Å); ✗ **robust losses refuted at n = 126**, firing `robust.py`'s own falsification condition |
| 3.9 | replication of the quantum positive | `qens_FINDINGS.md` | ✗ **DEMOTED** — the 27 "cells" are 9 targets × 3 seeds; at target level α = 1 is **null**, α = 0.05 **fails** concentration, and the control was charged 2,048 evaluations against the arm's **819,200** |
| 3.10 | adversarial consistency audit | `verify_FINDINGS.md` (744 lines) | 178 claims traced, **141 matched, 21 mismatched, 16 untraceable**; **13 blockers, all applied**; the benchmark confirmed untouched |
| 3.11 | the AMBER positive, replicated | `phys_FINDINGS.md` | **WEAKENED on accuracy, CONFIRMED on validity.** Replicates (5/5, sd 0.0076) and is **not** a multi-start artefact — the pipeline has no RNG, so **the coordinator's floor did not apply**; but **absent on the frame-reproducible 69%** (38W/46L, 3 frames), **erased by an exact rotation null** (+0.0117), and **4 targets silently non-converged** (8.9×10⁸ kcal/mol; `refine_coords` has no gate). Validity understated: **0.466→0.874 Ramachandran (116W/2L)**, **1.397→0.000 clashes (63W/0L)**, after two Sprint 14 reporting defects |

---

## RUNNING TALLY OF RETRACTIONS BY THE COORDINATOR

Kept as a visible count, because the brief requires the record to show corrections rather than a clean
narrative, and because the rate is itself information about how much the rest should be trusted.

1. **K1 on 8 targets** — the predicted arm reported at 2.792 Å; the full instrument gives 3.644 Å and
   **reverses the sign** of the comparison. *(Superseded, not retracted; the 8-target run was correct
   for 8 targets and is left in place unedited.)*
2. **K4, the quantisation ceiling** — retracted the same day. A uniform-grid bin lookup on a
   non-uniform grid. **75% of the claimed 1.06 Å was the bug**; the real cost is +0.381 Å.
3. **The 25.8% contraction figure** — reused in four places without checking; measured at **3.5%**.
4. **The multi-start seeding** — `hash(pdb)` is salted per process, so no multi-start constant was
   reproducible. Found by an agent, outside its brief.
5. **"Five times the headroom"** — in eight documents including the paper's abstract, and **not
   derivable from any number in the sprint**. The supported ratios are 3.19× and 2.07×.
6. **An ORACLE arm presented as "real errors"** in four documents, with the control that matters (a
   constant α-helix scoring 0.709) dropped in relaying.
7. **"2,048 draws"** — wrong in more than twenty places; the experiment drew 4,096.
8. **The in-tail ρ called the global ρ** — making "generative objectives have better global ordering"
   false as stated.
9. **The decorrelation ladder's "geometrically impossible" claim** — an artefact of eight targets;
   at n = 40 the uncontrolled ladder is monotone, and the requirement moves from 30% to **69%**.
10. **The claimed quantum positive** — led with in the dossier after the replication workstream had
    already demoted it; the control was charged **2,048 evaluations against the arm's 819,200**.
11. **Three sets of numbers quoted from before the seeding fix**, after writing a document forbidding
    exactly that.

**Two of the first four were caught by checking machinery against its own stated assumptions, not by
scientific intuition — and in both cases the wrong result was more plausible than the right one.
Items 5–11 were caught by an adversarial audit that loaded the artefacts and compared them to the
prose. None was caught by thinking harder about the science.**

