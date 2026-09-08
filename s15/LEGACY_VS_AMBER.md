# LEGACY VERSUS AMBER

The two energy models are mandated components of this programme, so their behaviour is measured
rather than assumed. This document states what each is, what each can and cannot do on 9–16 residue
peptides, where each has a defensible role, and — because it changes how every prior comparison must
be read — what the Phase 0 audit found wrong with the way they were previously compared.

---

## 1. WHAT THEY ARE

**Legacy** is a statistical potential with **eleven components** combined by a fixed weight vector.
The total is verified to be the `DEFAULT_WEIGHTS`-weighted sum of its components to within **2.7e-5
over 1.28e7 structures**, so the implementation is sound.

> **Its eleven weights were never fitted.** They are documented in the repository as "a
> variance-balanced starting point, not a fit." This is disclosed rather than presented as a design
> choice, and it is a strong claim to have to defend under review.

**AMBER** is ff14SB with GBn2 implicit solvent. Parameters verified against an independently
constructed ForceField at **maximum difference 0.0**; the golden energy −489.9138948277905 reproduces
**bit-exactly**.

---

## 2. FOUR CORRECTIONS FROM THE PHASE 0 AUDIT THAT CHANGE HOW COMPARISONS READ

These come first because every earlier Legacy-versus-AMBER number in this project was computed under
at least one of them.

### 2.1 AMBER costs 8–14 ms, not 28 ms — so budget matching over-charged it 2–3×

With memoisation defeated and hits counted, a single-point evaluation costs **8.3 ms at 127 atoms
rising to 23.3 ms at 237 atoms** — **8–14 ms on the n = 9–10 enumerated targets**. A repeated
structure costs **0.33 ms**, a 33× memoisation speedup. (`s13/qarch_lib.amber_energies` still
documents "~6 ms"; both figures were wrong, in opposite directions.)

**Direction of the error is safe.** Budget matching that charged 28 ms gave AMBER *fewer* evaluations
than a fair match. Every negative result about AMBER is therefore **understated**, not inflated.

### 2.2 The oracle-conditioning share was transposed

Measured across all 19 enumerated files: **21.7% oracle-conditioned** (`amber_kind == 2`), **40%
unbiased**, 39% prior-conditioned. Earlier records had this as "40% oracle-conditioned". The −0.401 Å
effect size was correct (−0.4009 on the nine, −0.4243 on nineteen).

### 2.3 The "unbiased" stratum is not unbiased in its tails

The index array **force-includes the ORACLE snap index**, which lands in the `kind == 0` stratum on
**16 of 19 files**. There it is:

- the **minimum-RMSD member on 6/19**;
- in the **RMSD-lowest 1% on 10/16**;
- **inside the < 1.5 Å in-band set on 14/19** — touching **40% / 25% / 17% / 15% / 11%** of all
  in-band pairs on five targets whose bands hold only 5, 8, 12, 13 and 19 members.

`kind == 0` is unbiased on the **mean** (+0.0003 Å) and **not in the tails**.

> **Binding rule: every tail, in-band, argmin or top-k statistic uses
> `amber_kind == 0 AND amber_idx != snap_index`, with the per-target n printed beside it.**

Direction is again safe — the headline 0.463 for AMBER could only have been flattered — but the
number is not clean.

### 2.4 One artefact is retracted outright

**Every AMBER row in `s14/results/obj_floor.json`.** It computes on the full labelled subsample with
**no `amber_kind` mask** (0.42 Å oracle-conditioned), contains no stratification string, and is
flagged `"complete": true`. Corrected values live in the Sprint 14 findings; the artefact was never
regenerated. Every other consumer is correctly masked.

---

## 3. NEITHER ENERGY RANKS THE NATIVE

The central negative, and it holds for both models.

- **The native is not the argmin.** In torsion space, Legacy's **certified optimum is worse than
  random**, and a zero-information constant α-helix beats the random control — so "beats random"
  demonstrates nothing on this instrument.
- **On matched pools with a radius-of-gyration control, AMBER is significantly *worse* than the
  control and places the native at the 51st percentile.** An earlier claim that physics ranks real
  geometry was **overturned** by that experiment.
- **All-atom AMBER reranking inside a real pool is worth +0.004 Å.** Nothing ranks within the pool.
- **Native anti-ranking is distributed, not a single-term defect** — it is a property of the total
  energy, not of one component that could be removed.
- **A polarizable force field does not fix it.** AMOEBA nearly flips one target and fails on another
  through van der Waals; the ceiling extends to polarizable physics.

### 3.1 The published context, which takes part of this claim

> **⚠ SPRINT 16 ADDITION — 2026-09-06, RETRACT workstream. A second source takes more of it, and it
> is about AMBER specifically rather than about contact potentials.**
>
> > Shao Q, Zhu W. **"Assessing AMBER force fields for protein folding in an implicit solvent."**
> > *Phys. Chem. Chem. Phys.* **20**(10):7206–7216 (2018). doi:10.1039/c7cp08010g. PMID 29480910.
>
> Enhanced-sampling MD over six AMBER force fields **including FF14SB and FF14SBonlysc, each with
> GB-Neck2** — **our exact pair** — on TC10b (20 res, α), HP35 (35, α), **1E0Q (17 res, β-hairpin)**
> and GTT (~35, β). Verbatim: *"the measured counterparts are significantly discrepant in the cases
> of larger or β-structured peptides (HP35, 1E0Q, and GTT)"*, and *"a combination ... able to
> describe all aspects of the folding transitions towards the native structures of all the considered
> peptides was not identified"*. Sprint 15 §2.4 downgraded this threat on the ground that the AMBER
> peptide literature does not test ff14SB. **That ground is false for this source, and the downgrade
> is reversed.** Abstract obtained verbatim from the Europe PMC record; full text not obtained.
>
> **Two defences that do NOT work, stated first:** the systems are not all longer than ours (1E0Q is
> a 17-mer), and folding free energy is not a weaker observable than potential-energy ranking — it is
> a **stronger** one, because it is Boltzmann-weighted. Their result is upstream of ours and
> **predicts** it.
>
> **What the source does not contain, and what we may therefore still claim:** no conformer ranking
> of a fixed pool, no enumerated space, no certified global optimum, no per-target CI, no
> filtering/ordering decomposition, and no polarizable comparison. The paper is also *mixed* — it
> endorses **FF14SBonlysc/GB-Neck2 as "a reasonably balanced combination"** and reports the α-helical
> 20-mer as working.
>
> **Corrected position for this document:** the bullets above are **confirmations of documented
> expected behaviour**, not discoveries. Cite Shao & Zhu as the expectation; claim the observable,
> the statistics, the AMOEBA extension, and — untouched by any of this literature — the
> **stereochemical-repair** result (Ramachandran 0.466 → 0.874, 116W/2L; clashes 1.397 → 0.000,
> 63W/0L). See `s15/LITERATURE.md` §2.4 and `s16/retract_FINDINGS.md`.

Roget et al. (arXiv:2606.21241) enumerate every peptide up to length 15 across 12,446 PDB structures
and report that the minimum-cost conformation has on average a **larger** RMSD than a randomly chosen
feasible one — with the mechanism that the Miyazawa–Jernigan contact potential was parameterised on
proteins above 50 residues, so ρ(cost, RMSD) is negative at 9–16 residues.

**What survives as ours** is the same phenomenon for a **certified optimum of an all-atom force field
plus a statistical potential in continuous torsion space**. Their length-parameterisation mechanism
cannot cover ff14SB. The claim must cite them as prior art for the phenomenon and be stated as an
extension to a different potential class — never as a discovery.

**And the counterpart matters.** For a *distance-restraint* objective at the same chain lengths, the
argmin over a 500-member pool is **−0.949 Å [−1.147, −0.753]** *better* than random, winning on 79%
of targets. Their pathology is a property of contact potentials, not of short-chain structure
prediction.

---

## 4. WHERE EACH IS GENUINELY DECORRELATED — AND WHY THAT DOES NOT HELP

Truth-partialled, Legacy and AMBER are genuinely independent: **+0.096** partial correlation, with
cross-channel error correlations of 0.04–0.26.

**And fusing them buys +0.004 to +0.011 Å.** The gain goes as the **square** of the weaker channel's
skill, and both channels' skills are near zero.

This sprint supplies a sharper and more useful account of why. The relevant question is not whether
two models err on the same *candidates* but whether their errors point at the same *wrong structure*.
For two estimates at mean distance `r` from the truth and `s` from each other, the midpoint sits at
~~`√(r² − (s/2)²)` — a law verified to **0.162 Å** across all 126 targets on the distance channels.~~
**⚠ CORRECTED, Sprint 16, 2026-09-06 (RETRACT):** the expression is the two-member Krogh–Vedelsby
(1995) ambiguity decomposition, `√((r₁²+r₂²)/2 − (s/2)²)`, and it is **exact** (2.65e−15 in a common
frame on all 126 targets). The published 0.162 Å was 54% an arithmetic-vs-quadratic mean error
(corrected residual **+0.075 Å**) and 46% a frame convention. `s16/retract_FINDINGS.md`.
The conclusion drawn from it is unaffected: fusion is transformative
only as `s → 2r`. Neither the Legacy/AMBER pair nor any other pair in this project comes close.

---

## 5. THE ROLES THAT SURVIVE

| | role | evidence |
|---|---|---|
| **Legacy** | **clash gate, and nothing else** | its only measured-honest role; its weights were never fitted, so a ranking role is indefensible |
| **AMBER** | **stereochemical repair**, and the objective in the constrained family | relaxation at k = 30 delivers a **large, robust validity gain**; its accuracy effect is **WEAKENED** and no longer quotable as an accuracy claim — see below |

### 5.1 The accuracy claim, WEAKENED

Replication under a pre-registered decision rule (five bootstrap redraws of the 75-window input;
threshold derived before any draw returned):

| | |
|---|---|
| per-draw effect | −0.0145 / −0.0235 / −0.0207 / −0.0226 / −0.0354 |
| mean, sd | **−0.0233 Å**, sd **0.0076 Å** — clears the pre-registered 0.0085 threshold |
| same sign | 5/5 (8/8 with the frame draws) |
| caveat, pre-registered then confirmed | r(projection cost, effect) = **−0.826**; the bootstrap axis **flatters** the effect |

**It is not a multi-start artefact** — this pipeline has **no RNG at all** and reproduces
bit-identically across interpreters (max |Δ| = 0.000e+00, n = 126), so the sprint's 0.081 Å
multi-start floor **cannot reach it**. A floor had to be built for this pipeline, and it is
**≈0.004 Å gated**.

**But it fails three of its own tests:**

1. **An exact null returns non-zero.** Relaxing the *same structure in a rotated lab frame* is zero by
   construction — ff14SB/GBn2, the restraint and every RMSD are rigid-invariant. Measured:
   **+0.0117 [−0.0034, +0.0383]**, sd 0.137, max **1.51 Å**. On that frame the effect becomes
   **−0.0102 [−0.0332, +0.0215]**, a CI **including zero**.
2. **The effect is absent on most of the instrument.** On the **84 targets (69%) where the minimiser
   is frame-reproducible** it is **−0.009 [−0.024, +0.006] with a *losing* 38W/46L**, replicated on
   three independent frames. **The whole effect sits on the 31% where the AMBER minimiser is not even
   frame-reproducible.** Sprint 14's "5/5 folds same sign" also does not fully replicate (one draw
   gives 3/5).
3. **NEW DEFECT — `refine_coords` has no convergence gate.** 1D6X, 1MF6, 2NB7 and 7BX2 end
   minimisation **above 1000 kcal/mol**, 1MF6 at **8.9 × 10⁸ kcal/mol** with 2 clashes, and are
   **silently scored into the published mean**. Those four alone move the null from −0.0005 to
   +0.0117; gated, the null collapses to −0.0005 [−0.0046, +0.0035] and the effect reproduces across
   frames to 0.0003 Å.

> **The accuracy claim may no longer be quoted as "−0.022 Å at valid geometry."**

**Cost, for completeness:** **no budget-matched comparison enters that figure** — it compares two
one-shot post-processing operators — and the single-point cost is the wrong number anyway: this arm
is a **full minimisation at 12.57 s mean**, ≈1,500× a single point and ≈2.8× the projection.

### 5.2 The VALIDITY claim, CONFIRMED — and badly understated, because of two Sprint 14 defects

Separating the two axes turned up two reporting defects in the original:

- **Arm A's Ramachandran fraction was scored on a DIFFERENT STRUCTURE from its RMSD** — the λ = 0.3
  arm's geometry quoted against the λ = 0 arm's RMSD.
- **Arm A's clash rate was a hard-coded placeholder** (`[0.0] * len`), not a measurement.

Measured correctly:

| quantity | baseline → relaxed | W/L |
|---|---|---|
| Ramachandran-favoured | **0.466 → 0.874** *(not 0.734 → 0.874)* | **116 / 2** |
| heavy-atom clashes | **1.397 → 0.000** | **63 / 0** |
| minimum heavy separation | **+0.843 Å** | 3 / 123 |

**This half is nowhere near any noise floor.** AMBER's defensible role is **stereochemical repair**,
and on that axis it is larger than the record said. What it is not is an accuracy gain.

**Legacy in generation actively hurts.** Used as a generative term rather than a gate it damages the
answer; used as a corrector of the selector it is huge on single cells and worth only −0.16 Å
(p = 0.23) aggregated.

---

## 6. A NUMERICAL CAUTION FOR ANY SPECTRAL ANALYSIS OF THESE ENERGIES

**37.6% of cached AMBER single points exceed 1e6 kcal/mol**, 3.2% exceed 1e12, and the maximum is
**2.6e20**. The 99th percentile is itself 1e12–1e15.

So **winsorising at the 99th percentile is not sufficient conditioning**, and a Pauli or Walsh
spectrum of the unconditioned energy measures its worst steric clash and nothing else — the top-10
configurations carry 99.6% of raw AMBER's Walsh variance, with the weight landing on
`Binomial(m, 1/2)` exactly. Condition monotonically before any spectral claim.

A related category error to avoid: the locality theorem (`d_ij` depends on exactly the `j − i − 1`
residues between `i` and `j`, agreement 1.0000) applies to the **Cα geometry**. Both energies are
**full-register**; neither decomposes over a bounded neighbourhood. "AMBER is less local than Legacy"
is not a true statement, and not a false one — it is a comparison of two things that are equally
non-local.

---

## 7. WHAT THIS DOCUMENT DOES NOT CLAIM

- It does not claim either energy is broken. Both implementations are verified correct (§1); what
  fails is their ability to **rank near-native peptide conformations**, which is a claim about the
  models' domain, not their code.
- It does not claim AMBER is worse than Legacy or the reverse. On the axis that matters — in-band
  ranking — both are at chance, and the budget-matched comparisons that appeared to separate them
  were computed with AMBER under-charged 2–3× and with an oracle index in its "unbiased" stratum.
- It does not claim the −0.022 Å relaxation gain scales. It was measured at k = 30 on this
  instrument.
