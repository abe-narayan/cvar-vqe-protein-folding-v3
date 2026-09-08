# SPRINT 15 — PHYS findings

**One question: is the programme's only positive physics result real, or is it below the noise
floor?** The result under test is AMBER ff14SB/GBn2 restrained relaxation at `k_restraint = 30`,
written up in `s15/LEGACY_VS_AMBER.md` §5 and `s15/ARCHITECTURE.md` as **−0.022 Å
[−0.036, −0.009]** against the incumbent ideal-geometry projection.

Tiering: **DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **HYPOTHESIS** / **REFUTED**.

Module: `s15/phys_repl.py` (one module, 5 modes). Artefacts: `s15/results/phys_repl_canon0.json`,
`phys_repl_boot{1..5}.json`, `phys_repl_frame{1,2,3}.json`, `phys_repl_startperm{1,2,3}.json`,
`phys_nullstruct_frame{1,2,3}.json`, `phys_repl_report.json`.

---

## VERDICT

> **WEAKENED.** The effect survives replication — it keeps its sign on **8 of 8** independent
> draws and its across-draw sd is **0.0076 Å**, so the RESTRAINT workstream's 0.081 Å floor
> **does not transfer to this pipeline and must not be applied to it** — but it fails three of
> its own tests: an exact null that *must* be zero returns **+0.0117 [−0.0034, +0.0383]** and
> destroys the effect's significance on that draw; **69% of the instrument contributes nothing
> (−0.009 [−0.024, +0.006], 38W/46L, replicated on three frames)** with the whole effect carried
> by the 31% of targets where the AMBER minimiser is not even frame-reproducible; and four
> targets carry a **silently non-converged minimisation** (final energy up to 8.9 × 10⁸
> kcal/mol) that no gate anywhere in the pipeline catches.
>
> **The accuracy claim may not be quoted as "−0.022 Å at valid geometry" any longer.**
> The *validity* claim, by contrast, is **CONFIRMED and was badly understated**: measured against
> the correct baseline it is **0.466 → 0.874 Ramachandran-favoured on 116 targets to 2**, and
> **1.40 → 0.00 heavy-atom clashes on 63 to 0**. That half is nowhere near any noise floor.

| claim | tier |
|---|---|
| the k = 30 pipeline is deterministic; the 0.081 Å multi-start floor cannot reach it | **DEMONSTRATED** |
| the effect's sign and magnitude replicate across 5 bootstrap input draws (sd 0.0076 Å) | **DEMONSTRATED** |
| the exact frame null is non-zero and can erase the effect's significance | **DEMONSTRATED** |
| `refine_coords` has no convergence gate; 4/126 targets are silently non-converged | **DEMONSTRATED** (new defect) |
| the effect is absent on the frame-reproducible 69% of targets | **DEMONSTRATED**, ×3 frames |
| Sprint 14 scored arm A's Ramachandran on a different structure from its RMSD | **DEMONSTRATED** (defect) |
| Sprint 14's arm-A clash rate was a hard-coded placeholder | **DEMONSTRATED** (defect) |
| the validity gain (Ramachandran, clashes, min heavy separation) | **DEMONSTRATED**, large |
| no budget-matched comparison enters the −0.022 Å figure | **DEMONSTRATED** |
| permuting the projection's four fixed starts can flip basins | **REFUTED** (my own hypothesis) |
| the effect is a lawful function of the projection cost, not noise | **HYPOTHESIS**, r = −0.826 on n = 6 draws |

---

## 0. PRE-REGISTERED DECISION RULE, recorded before the replication draws returned

Written to this file while `canon` was still running and before a single `boot` draw existed, so
it cannot be retrofitted. Reproduced verbatim, then scored.

| verdict | condition | outcome |
|---|---|---|
| **CONFIRMED** | both exact nulls sit at a floor ≪ 0.022 Å (null \|mean\| and CI half-width < 0.005 Å), **and** the effect keeps its sign on ≥ 4 of 5 replication draws with an across-draw sd small against \|mean\| | **half met.** Replication: 5/5 same sign, sd 0.0076. Nulls: startperm exactly 0; **frame null \|mean\| = 0.0117 and CI half-width 0.021 — FAILS** |
| **WEAKENED** | the sign holds but the across-draw sd is comparable to the effect (≳ \|mean\|/2), **or** a majority of replication draws carry a CI including zero | sd 0.0076 is 0.33 × \|mean\| — under the bar; 4/5 boot CIs exclude zero |
| **REFUTED** | the sign flips on ≥ 2 of 5 replication draws, **or** a null floor is of order 0.022 Å, **or** the pooled replication effect's CI comfortably includes zero | sign never flips; pooled boot is significant; **but the ungated null floor's CI reaches +0.038 Å, which IS of order 0.022** |

**The rule does not return a clean verdict, and saying so is part of the result.** The
replication axis passes cleanly and the null axis fails cleanly, so the honest reading is
**WEAKENED** with the reason stated: the effect is reproducible, and the machinery that measures
it is not.

**The quantitative threshold, derived before the draws returned.** The within-draw standard error
is `0.0789 / √126 = 0.00703 Å`, so the quoted effect is 3.1 SE from zero. Adding a between-draw
sd `s_d`, the effect stops being two standard errors from zero once **s_d > 0.0085 Å**.
**Measured: s_d = 0.00759 Å**, total SE 0.01035, effect/SE = **2.25**. It clears the bar it was
given, with almost nothing to spare. For scale, the RESTRAINT workstream's imported floor
(0.081 Å) is **ten times** the measured spread.

---

## 1. THE QUESTION AS POSED IS NOT ANSWERABLE, BECAUSE THIS PIPELINE HAS NO START DRAWS

**DEMONSTRATED.** The brief asked for the k = 30 effect to be re-measured "across at least five
independent start/seed draws", on the model of the RESTRAINT workstream's false-positive floor
(R1.5). That floor has a specific mechanism: `np.random.default_rng(hash(pdb) % 2**32)` seeded a
**random multi-start** L-BFGS, `hash()` is salted per process, so every interpreter drew a
different start set and a comparison that is zero by construction returned +0.081 [+0.014, +0.169].

**No such mechanism exists in the k = 30 path.** Every stage was audited for randomness:

| stage | randomness |
|---|---|
| `s14.avgspace.top75_windows` | none — a cache read of the shipped window filter |
| `I.coordinate_average` | none — medoid + Kabsch superposition |
| `core.geometry.build_backbone_batch` | none |
| `I.project` → `core.project.lam_path(multi=True)` | **none.** `multi=True` means "all four of the FIXED generic starts". `core.project.STARTS` is the literal 4-tuple `(-120,130), (-57,-47), (-139,135), (-75,145)` |
| `torsion_lib2.library_for` | `seed = 1`, a fixed constant, and it places **side chains only** |
| `core.amber.refine_coords(steps=0)` | `random.seed(0)` before hydrogen placement (a fixed constant); OpenMM CPU platform, `Threads=1`, `DeterministicForces=true`; no thermostat, no MD |

`core/amber.py` states this in its own docstring: *"Determinism comes from the CPU platform at
Threads=1 with DeterministicForces, and from the starting structure being a pure function of its
inputs: no random seeding, no thermostat, no MD."*

**Verified across processes, which is the check A2 says the sprint never had.** The whole
construction was re-run in a fresh interpreter, on all 126 targets, and compared against the
Sprint 14 artefact:

| arm | max \|Δ\| over 126 targets |
|---|---|
| arm A, the ideal-geometry projection | **0.000e+00 Å** |
| arm B + AMBER k = 30 | **0.000e+00 Å** |

**Bit-identical on both arms**, different interpreter, different day, `PYTHONHASHSEED` unset in
both.

> **Correction to the coordinator's record (R1.5's table).** AMBER relaxation at k = 30 is listed
> there as an effect that "must now be replicated across start draws before it is quoted".
> Applied literally that is a category error — there are no start draws in this pipeline, and
> importing the 0.081 Å figure would have priced this result against a noise source that cannot
> reach it. **This does not make the effect real.** It means the floor had to be *constructed* for
> this pipeline rather than imported, which is what §3 does — and the constructed floor turns out
> to bite for an entirely different reason.

---

## 2. CROSS-PROCESS REPRODUCTION — and the input-matching confound, which runs the safe way

**DEMONSTRATED.** Bit-identical, both arms, n = 126 (§1).

One matching detail, checked because the two arms are not fed the same object. Arm A projects
`C_ca = I.coordinate_average(W)`; the AMBER arm relaxes `avg`, the all-atom average whose CA block
differs from `C_ca` by a **mean 0.159 Å, max 0.373 Å** per-coordinate deviation (the per-window
ideal-geometry rebuild residual, most of which cancels in the mean). So AMBER's input is not
exactly arm A's input.

Measured, the induced RMSD offset is **+0.0015 Å [−0.0051, +0.0082], 59W/67L** — AMBER's starting
structure is very slightly *worse*. Adjusting for it moves the effect from −0.0221 to
**−0.0236 [−0.0360, −0.0119], 78W/48L, median −0.0137**. The confound runs *against* AMBER, so
correcting it strengthens the claim. Recorded because it was not in the original record.

---

## 3. THE NULL, BUILT FOR THIS PIPELINE

Two comparisons that are **zero by construction** in the AMBER path, run through the identical
machinery — the analogue of "maximum likelihood against a constant-width Gaussian is least
squares".

### 3.1 Exact null A — the same structure relaxed in a different LAB FRAME

Apply a random proper rotation and translation to the averaged backbone before handing it to
both arms. Every term is rigid-invariant by construction: ff14SB and GBn2 are translation- and
rotation-invariant; the positional restraint is `0.5·k·|x − x₀|²` with `x₀` set from the same
(rotated) input, so it rotates with the structure; `fit_prior` minimises a Kabsch RMSD; and
`ca_rmsd` is a Kabsch RMSD. **The frame-to-frame difference is mathematically zero.** Anything it
measures is this pipeline's own numerical floor. (Reflections are excluded — `det = +1` is
enforced — because a reflection is *not* a symmetry of ff14SB.)

**DEMONSTRATED. It is very far from zero, and it found a defect nobody was looking for.**

| exceedance of the null \|Δ RMSD\|, frame 1 | targets |
|---|---|
| > 1e−6 Å | **126 / 126** |
| > 1e−4 | 110 |
| > 1e−3 | 56 |
| > 1e−2 | 28 (22%) |
| > 5e−2 | 8 |
| > 1e−1 | 3 |
| > 1 Å | **1 (1MF6, +1.5100 Å)** |

Arm A is, as predicted, essentially frame-invariant (mean −0.00016, sd 0.0019, max 0.0142) — its
objective is a Kabsch RMSD, so a rigid motion cannot reach it. **All of the noise is in the AMBER
minimiser.**

### 3.1a NEW DEFECT: `refine_coords` has no convergence gate

**DEMONSTRATED, and it is in no prior record.** Four targets — **1D6X, 1MF6, 2NB7, 7BX2** — end
their restrained minimisation at a final energy **above 1000 kcal/mol** (1785, −439/8.9e8, 3685,
1493), and 1MF6 in the rotated frame ends at **8.9 × 10⁸ kcal/mol**, with two sub-2 Å clashes and
a 0.13 geometry deviation. `refine_coords(steps=0)` asks OpenMM to minimise until converged and
**never checks that it did**; the resulting RMSD is scored into the published mean like any other.

Those four targets alone move the exact null's mean from **−0.0005 to +0.0117 Å** and its
per-target sd from **0.024 to 0.137**.

> **Recommendation, and it is cheap.** Gate `refine_coords` on its own final energy — or apply the
> `bond + angle > 1000` strain gate that already exists in the codebase but is not used here — and
> re-emit. Nothing in the pipeline currently notices a failed minimisation. Gating changes the
> canon effect from −0.02207 to **−0.02339 (n = 123)**, so the fix is free of direction risk.

### 3.1b What the null says, across three independent frames

| stratum | n | effect, frame 0 | effect, frame 1 | **the null** (zero by construction) |
|---|---|---|---|---|
| **ALL, as published** | 126 | −0.02207 [−0.03591, −0.00857] 67W/59L | **−0.01018 [−0.03323, +0.02148]** | **+0.01173 [−0.00342, +0.03827]**, sd 0.137 |
| **converged only** | 122 | −0.02345 [−0.03713, −0.00954] 65W/57L | −0.02378 [−0.03832, −0.00953] | −0.00049 [−0.00464, +0.00349], sd 0.024 |
| converged, minimiser **stable** (\|ΔE\| ≤ 0.5) | 84 | **−0.00912 [−0.02374, +0.00553] 38W/46L** | −0.00902 [−0.02356, +0.00578] | +0.00007 [−0.00015, +0.00029], sd 0.001 |
| converged, minimiser **unstable** | 38 | −0.05514 [−0.08321, −0.02540] 27W/11L | −0.05641 [−0.08801, −0.02487] | −0.00173 [−0.01527, +0.01140], sd 0.043 |

Replicated on two further independent frames (`phys_nullstruct_frame{2,3}.json`):

| frame draw | the null, n = 126 | effect, whole instrument | **effect on the frame-STABLE stratum** |
|---|---|---|---|
| 0 (canon) | — | −0.02207 [−0.03591, −0.00857] 67W/59L | — |
| 1 | **+0.01173** [−0.00342, +0.03827], sd 0.137, max 1.510 | −0.01018 [−0.03323, **+0.02148**] | −0.00912 [−0.02374, **+0.00553**] **38W/46L** (n = 84) |
| 2 | **−0.00432** [−0.00894, +0.00015], sd 0.026, max 0.138 | −0.02639 [−0.04097, −0.01245] | −0.00453 [−0.01954, **+0.01013**] **36W/49L** (n = 85) |
| 3 | **+0.00021** [−0.00486, +0.00520], sd 0.029, max 0.172 | −0.02187 [−0.03577, −0.00800] | −0.00847 [−0.02386, **+0.00691**] **38W/47L** (n = 85) |

**Three things follow, and they do not all point the same way.**

1. **As the pipeline currently stands, the effect is inside its own floor.** The exact null's 95%
   interval reaches **+0.038 Å**, and re-measuring in a rotated coordinate frame — a change that
   is mathematically no change at all — moves the headline from −0.0221 to **−0.0102 with a
   confidence interval that includes zero**.
2. **The cause is identified and fixable, and fixing it makes the effect slightly larger.** With
   the non-converged minimisations gated out, the null collapses to
   **−0.0005 Å [−0.0046, +0.0035]** — a genuine floor of about **0.004 Å**, five times below the
   effect — and the effect reproduces across frames to **0.0003 Å** (−0.02345 versus −0.02378).
   **So this is a data-quality defect, not a demonstration that the physics is noise.**
3. **But the null also localises the effect, and that is the damaging part.** On the **84 targets
   (69%) where the minimiser is frame-reproducible to 0.5 kcal/mol, the effect is
   −0.009 [−0.024, +0.006] and loses 38W/46L** — indistinguishable from zero with a *losing*
   win/loss record, and the same on all three frames (−0.009 / −0.005 / −0.008; 38/46, 36/49,
   38/47). The entire published effect is carried by the **38 targets (31%) where the minimiser is
   frame-sensitive**; on those it is −0.055 [−0.083, −0.025] and it does reproduce
   (−0.0551 versus −0.0564).

Point 3 admits two readings and the data discriminate. *Either* the effect is noise from an
unreliable minimiser — **refuted**, because on that stratum it reproduces across an independent
frame to 0.0013 Å — *or* the relaxation genuinely helps only where it has real work to do, i.e.
where the restrained minimum is far enough from the input that the basin choice is delicate.
**The second reading is what the data support**, and it is a sharper statement than "AMBER
relaxation is worth −0.022 Å". It is also a much narrower claim.

### 3.2 Exact null B — permuting `fit_multi`'s four fixed starts

Null A perturbs floating point only; both arms stay in the same basin, so it certifies the
numerical floor and nothing more. R1.5's floor came from something stronger — a multi-start
optimiser landing in **different basins**. The only multi-start in this pipeline is `fit_multi`,
which loops over the four fixed `STARTS` and keeps the lowest objective with a **strict `<`**.
Permuting their order is an argmin over an unordered set and is therefore zero by construction —
*unless* two starts tie within float error, in which case the strict `<` breaks the tie **by
order** and the emitted structure flips basins. That is this project's own recorded trap ("an
argmin over a tied signal silently read the ORACLE sort order") turned into a null.

**Result: exactly zero. My own hypothesis is REFUTED.**

| draw | permutation | max \|Δ\|, lam = 0 arm | max \|Δ\|, lam = 0.3 arm | targets with \|Δ\| > 0 |
|---|---|---|---|---|
| startperm 1 | [0, 3, 2, 1] | **0.000e+00** | 0.000e+00 | **0 / 126** |
| startperm 2 | [0, 3, 1, 2] | **0.000e+00** | 0.000e+00 | **0 / 126** |
| startperm 3 | [0, 1, 3, 2] | **0.000e+00** | 0.000e+00 | **0 / 126** |

The four generic starts are never within float error of one another on this instrument, so the
strict `<` never has to break a tie and no basin flip is available. **Arm A contributes nothing
whatever to this pipeline's noise floor.** I proposed this null expecting it to be the strong one;
it is the empty one, and the useful consequence is that **every ångström of irreproducibility in
this comparison lives in the AMBER minimiser.**

---

## 4. REPLICATION ACROSS INDEPENDENT DRAWS

Since no start draw exists (§1), the replication axis is a **bootstrap of the 75 windows** that
enter the coordinate average — 75 draws with replacement, ~47.8 distinct. Both arms consume the
identical resampled set, so each draw stays start-matched, but the object being projected and
relaxed is a different, equally legitimate draw of the pipeline's input. Seeded with
`SD.stable_rng("phys15", "boot", draw, pdb)`.

**DEMONSTRATED**, n = 126 per draw, 5 draws, plus the 3 frame draws and 3 startperm draws as
degenerate controls.

| draw | n | arm A | arm B + AMBER | **mean** | CI 95% | median | W/L | mean/sd | conc. verdict | gated (E ≤ 1e3) |
|---|---|---|---|---|---|---|---|---|---|---|
| **canon (draw 0)** | 126 | 3.205 | 3.183 | **−0.0221** | [−0.0359, −0.0086] | −0.0069 | 67/59 | −0.280 | PASS | −0.0234 (n=123) |
| boot 1 | 126 | 3.222 | 3.208 | **−0.0145** | [−0.0318, **+0.0028**] | −0.0048 | 68/58 | −0.147 | PASS | −0.0198 (n=122) |
| boot 2 | 126 | 3.225 | 3.202 | **−0.0235** | [−0.0403, −0.0084] | −0.0027 | 66/60 | −0.255 | PASS | −0.0220 (n=124) |
| boot 3 | 126 | 3.205 | 3.184 | **−0.0207** | [−0.0379, −0.0026] | −0.0055 | 68/58 | −0.207 | PASS | −0.0234 (n=124) |
| boot 4 | 126 | 3.222 | 3.199 | **−0.0226** | [−0.0417, −0.0043] | −0.0082 | 70/56 | −0.208 | PASS | −0.0182 (n=124) |
| boot 5 | 126 | 3.227 | 3.192 | **−0.0354** | [−0.0511, −0.0210] | −0.0096 | 77/49 | −0.417 | PASS | −0.0352 (n=124) |
| frame 1 | 126 | 3.205 | 3.195 | −0.0102 | [−0.0332, **+0.0215**] | −0.0018 | 66/60 | −0.064 | PASS | — |
| frame 2 | 126 | 3.205 | 3.179 | −0.0264 | [−0.0410, −0.0125] | −0.0043 | 68/58 | −0.322 | PASS | — |
| frame 3 | 126 | 3.205 | 3.183 | −0.0219 | [−0.0358, −0.0080] | −0.0010 | 64/62 | −0.265 | PASS | — |
| startperm 1–3 | 126 | 3.205 | 3.183 | −0.0221 | [−0.0359, −0.0086] | −0.0069 | 67/59 | −0.280 | PASS | — |

**Across-draw summary:**

| axis | n draws | mean | **sd** | range | sign |
|---|---|---|---|---|---|
| **boot** (the replication axis) | 5 | **−0.02332** | **0.00759** | [−0.0354, −0.0145] | **5/5 negative** |
| frame (the exact-null axis) | 3 | −0.0195 | 0.0084 | [−0.0264, −0.0102] | 3/3 negative |
| startperm (degenerate by construction) | 3 | −0.0221 | 0.0000 | — | — |

**The sign never moves. The significance does.** 4 of 5 boot draws and 2 of 3 frame draws exclude
zero; boot 1 and frame 1 do not. Total SE including draw variance = 0.01035 Å, giving
**effect / SE = 2.25** — still nominally significant, and only just.

### 4.1 The confound in the boot axis, pre-registered and then confirmed

A bootstrap keeps only ~48 of 75 distinct windows, so the average is noisier, the backbone
contracts more, and the projection cost `A − B` rises. The effect is known to scale with that
cost (ρ = −0.51 across targets on canon), so a *larger* effect on boot draws would be mechanism,
not evidence of stability. Measured:

| draw | distinct windows | projection cost A − B | effect |
|---|---|---|---|
| canon | 75.0 | +0.1554 | −0.0221 |
| boot 1 | 47.9 | +0.1540 | −0.0145 |
| boot 2 | 47.7 | +0.1643 | −0.0235 |
| boot 3 | 47.6 | +0.1594 | −0.0207 |
| boot 4 | 47.5 | +0.1538 | −0.0226 |
| boot 5 | 48.2 | **+0.1693** | **−0.0354** |

**r(projection cost, effect) = −0.826 across the six draws.** **HYPOTHESIS** (n = 6 draws, so this
is suggestive, not established): the effect is a lawful function of how much re-expansion the
projection has to do, not a random quantity. boot 5 has both the largest projection cost and the
largest effect; boot 1 and boot 4 have the smallest costs and the smallest effects. If this holds,
the boot axis is *not* a neutral resampling — it perturbs the input in the direction that
amplifies the effect, and the replication is correspondingly flattered.

### 4.2 Per-fold and FAIL18 — Sprint 14's "5/5 folds same sign" does NOT fully replicate

FAIL18 is reported for cross-sprint comparability only; per the Phase 0 audit it is a
`BAND = 1.5` threshold artefact and no argument may rest on its membership.

| draw | fold 0 | fold 1 | fold 2 | fold 3 | fold 4 | FAIL18 | other | folds negative |
|---|---|---|---|---|---|---|---|---|
| canon | −0.0213 | −0.0138 | −0.0165 | −0.0195 | −0.0357 | −0.0113 | −0.0239 | **5/5** |
| boot 1 | −0.0279 | **+0.0327** | −0.0183 | **+0.0020** | −0.0490 | **+0.0128** | −0.0190 | **3/5** |
| boot 2 | −0.0249 | −0.0130 | −0.0260 | −0.0122 | −0.0368 | −0.0335 | −0.0218 | 5/5 |
| boot 3 | −0.0339 | **+0.0063** | −0.0148 | −0.0117 | −0.0424 | −0.0190 | −0.0210 | 4/5 |
| boot 4 | −0.0239 | −0.0108 | −0.0159 | −0.0148 | −0.0420 | −0.0083 | −0.0249 | 5/5 |
| boot 5 | −0.0506 | −0.0205 | −0.0267 | −0.0186 | −0.0541 | −0.0633 | −0.0307 | 5/5 |

Sprint 14 wrote that the per-fold table "is the genuinely informative one" and that 5/5 same sign
"is the real evidence for uniformity, and it is independent of the drop-top machinery". **Under
redraw that evidence is weaker than stated**: two of five bootstrap draws lose a fold (boot 1
loses two, and flips FAIL18 positive). Fold 4 is the only fold negative on every draw and it is
consistently the largest contributor.

---

## 5. ACCURACY VERSUS VALIDITY — and two Sprint 14 reporting defects found while separating them

### 5.1 Defect: arm A's Ramachandran fraction was scored on a DIFFERENT STRUCTURE from its RMSD

**DEMONSTRATED.** `s14/ener_avgrefine.run_target` scores

```python
rec["projected_rmsd"]    = I.ca_rmsd(pr["fit_ca"], nat)        # the lam = 0 arm
rec["projected_rama_ok"] = rama_ok(pr["phi"], pr["psi"])       # the lam = 0.3 arm
```

`I.project` returns the unpenalised `fit_ca` **and** the `ramah`-penalised `ca`/`phi`/`psi`, and
throws away the lam = 0 arm's torsions — so the Ramachandran statistic *had* to come from the
penalised arm. Those are two different structures, and the penalised one is the arm explicitly
optimised for Ramachandran plausibility. **The published validity comparison ("AMBER 0.874 against
the projection's 0.734") therefore compares AMBER against a structure whose RMSD is not the one
quoted.** Same class of defect as the project's recorded `rmsd-reporting-basis-mismatch`.
`s15/phys_repl.project_both` returns both arms' torsions so both statistics are on one structure.

### 5.2 Defect: arm A's clash rate was a hard-coded placeholder, not a measurement

**DEMONSTRATED.** `s14/ener_avgrefine.report` calls

```python
line("A projected", A, [0.0] * len(ok), [r["projected_rama_ok"] for r in ok], [0.0] * len(ok))
```

The third and fifth arguments are the geometry-deviation and clash columns. The geometry column is
legitimately zero (the ideal-geometry builder fixes bonds and angles by construction); **the clash
column is not** — an ideal-geometry chain can and does clash. Measured on all 126 targets, arm A
carries a **mean of 1.397 heavy-atom clashes below 2.0 Å**, and AMBER beats it **63 targets to 0**.
The published row "A projected … clashed 0.000" hid a real result.

Both defects run *against* the incumbent, so correcting them can only make AMBER's validity
advantage look larger. The numbers were still wrong.

### 5.3 The split

**DEMONSTRATED**, n = 126, `phys_repl_canon0.json`. Every row paired against arm A (the lam = 0
projection), scored on the **same structure** as its RMSD.

| axis | arm A | arm B + AMBER k = 30 | paired diff | CI 95% | W/L | median |
|---|---|---|---|---|---|---|
| **ACCURACY** — CA-RMSD to native *(ORACLE)* | 3.2052 | 3.1831 | **−0.0221** | [−0.0359, −0.0086] | 67/59 | −0.0069 |
| **VALIDITY** — 1 − Ramachandran-favoured | 0.5342 | 0.1264 | **−0.4077** | [−0.4529, −0.3628] | **116/2** | −0.4365 |
| *(the same, against the lam = 0.3 arm S14 used by mistake)* | 0.2665 | 0.1264 | −0.1400 | [−0.1885, −0.0901] | 59/26 | +0.0000 |
| **VALIDITY** — heavy-atom clashes < 2.0 Å | 1.397 | **0.000** | **−1.397** | [−1.754, −1.071] | **63/0** | −0.500 |
| **VALIDITY** — min heavy-atom separation (Å) | 1.963 | 2.806 | **+0.843** | [+0.766, +0.919] | 3/123 | +0.785 |
| geometry — rms rel. bond/angle deviation | 0 by construction | 0.0192 (max 0.0618) | — | — | — | — |

Ramachandran-favoured fractions, all four objects: **arm A (lam = 0) 0.466**, arm A (lam = 0.3)
0.733, raw average 0.837, **AMBER 0.874**.

**The two halves are not comparable in strength and must not be presented as one claim.**

- The **validity** half is overwhelming and was **understated by more than a factor of two**.
  Against the correctly matched baseline the Ramachandran gain is **0.466 → 0.874 on 116 targets
  to 2**, and arm A carries **1.40 sub-2.0 Å clashes per target** against AMBER's zero, on 63
  targets to 0. Neither statistic is anywhere near a noise floor — both are measured in their own
  units, but each has a confidence interval whose *near* bound is many multiples of its own
  standard error (Ramachandran 116W/2L; clashes 63W/0L; min separation 3W/123L), where the
  accuracy row is 67W/59L.
- The **accuracy** half is 0.7% of a 3.2 Å baseline, near-even on wins, and is what §3 and §4 are
  about.

> **The claim that survives without qualification is a STEREOCHEMISTRY claim, not an accuracy
> claim.** Restrained AMBER relaxation produces a materially better-formed structure than the
> ideal-geometry projection — clash-free, and Ramachandran-favoured on 87% of residues against
> the projection's 47%. That is worth having, it is robust, and it is not what the programme has
> been quoting.

---

## 6. COST ACCOUNTING — no budget-matched comparison enters this result

**DEMONSTRATED.** The brief asked whether the k = 30 result involves a budget-matched comparison
and, if so, to recompute it at the corrected AMBER cost.

**It does not.** `s14/ener_avgrefine` compares two *post-processing operators applied once to the
same structure* — the ideal-geometry projection versus one restrained relaxation. There is no
search, no evaluation count, and no budget on either side, so neither the "equal objective
evaluations" nor the "equal computational cost" convention is engaged. **The 28 ms → 8–14 ms
correction moves nothing here.**

**And the corrected single-point figure is the wrong number for this arm anyway.** The audited
8.3 ms (127 atoms) to 23.3 ms (237 atoms) is an AMBER **single point**. This arm is a full
restrained **minimisation** (`refine_coords(steps=0, tolerance=1.0)`). Measured on this machine,
n = 126, one process at `Threads=1` (contended box, so upper bounds):

| operator | mean | median | max |
|---|---|---|---|
| arm B + AMBER k = 30 | **12.57 s** | 11.84 s | 38.43 s |
| arm A, `lam_path(multi=True)` — nine L-BFGS runs | **4.46 s** | 4.50 s | 8.73 s |

So the AMBER arm is **≈1,500× an AMBER single point** and **≈2.8× the projection it is compared
against**. Should this ever be made a budget-matched comparison, AMBER is the *expensive* arm by a
factor of three, and the corrected single-point cost is not the relevant number.

*(For the reproducibility document: `s13/qarch_lib.amber_energies` still documents "~6 ms"; the
brief's `refine_coords` figure of 9.1 s is consistent with the 11.84 s median measured here under
load.)*

---

## 7. THE BINDING DATA RULE — declared not applicable, deliberately

**No cached AMBER subset is read anywhere in this workstream.** `s13/results/qarch_enum_*.npz`
and `s14/cache/` are not opened; every AMBER number here is computed fresh through
`core.amber.refine_coords` on a structure built in-process. There is therefore no `amber_kind`
stratum, no `amber_idx`, and no `snap_index` to exclude, and no tail / in-band / argmin / top-k
statistic is taken over that cache. The per-target n that *is* reported is the number of targets
in each paired comparison, printed in every table (126 throughout, 122–124 under the convergence
gate, 84–85 / 38 in the stability strata).

`s14/results/obj_floor.json` was **not read**, per the retraction. The 60-target protected
benchmark was not touched.

---

## 8. WHAT THE EFFECT IS, MECHANICALLY — it is a price of validity, not an accuracy gain

**DEMONSTRATED**, reproduced here from the Sprint 14 artefact.

| structure | CA-RMSD to native *(ORACLE)* | geometry |
|---|---|---|
| **B**, the raw all-atom coordinate average | **3.050 Å** | broken — 25.8% bond contraction, 0.248 rms rel. dev., 57/126 targets clashing |
| **A**, the ideal-geometry projection of B | 3.205 Å (+0.155) | ideal bonds/angles; **1.40 clashes/target**; Rama 0.466 |
| **B + AMBER k = 30** | 3.183 Å (+0.133) | 0.019 rms rel. dev.; **0 clashes**; Rama 0.874 |

Both operators make the structure **worse** than the unrefined average. The −0.022 Å is the
difference between two prices for restoring valid geometry — 0.133 Å versus 0.155 Å — not physics
finding a better structure. Sprint 14's own write-up says this ("the gain is a strain allowance"),
and the correlations confirm it: ρ(projection cost A − B, effect) = **−0.510** across targets, and
r = **−0.826** across draws (§4.1).

**The consequence for how it may be described.** "AMBER relaxation improves the emitted structure
by −0.022 Å" is true of the comparison as run and misleading as a summary: the operator that
produces the best RMSD here is *doing nothing at all*, and both refiners are paying for
stereochemistry. The defensible phrasing is **"restrained relaxation restores valid geometry for
about 0.02 Å less than the projection does, and produces a much better-formed structure while
doing it."**

---

## 9. THE FRAGILITY THAT WAS ALREADY VISIBLE IN THE ORIGINAL ARTEFACT

**DEMONSTRATED**, computed on `s14/results/ener_avgrefine.json`, n = 126, all reproduced.

| statistic | value | reading |
|---|---|---|
| mean | −0.02207 | the quoted effect |
| **median** | **−0.00694** | the mean is **3.2× the median** |
| per-target sd | 0.0789 | |
| **mean/sd** | **−0.280** | the concentration test has almost no power at this ratio |
| W/L | 67/59 | **sign test p = 0.533** — undetectable as a sign effect |
| Wilcoxon signed-rank | p = 0.0164 | |
| paired t-test | p = 0.0021 | |
| bootstrap CI | [−0.0359, −0.0086] | excludes zero |
| 10% / 20% trimmed mean | −0.0183 / −0.0144 | |
| drop-top-5 / drop-top-10 | −0.0144 / −0.0080 | |
| null-calibrated concentration verdict | **PASS (not concentrated)** | observed at the 62nd percentile of a uniform-effect null on both statistics; **PASS on all 11 draws in §4** |
| leave-one-target-out range | [−0.0241, −0.0203] | no single target carries it |

**This is the median-vs-mean early warning firing exactly as the ledger says it should**: a
near-even W/L with a CI excluding zero. The null-calibrated check returns PASS on every draw, and
the ledger's own caveat applies — at mean/sd = 0.28 that check *cannot settle the question either
way*. What carried the result was the CI on the mean and the 5/5 fold consistency; §4.2 shows the
fold consistency is weaker than reported, and §3.1b shows the mean is carried by 31% of targets.

**And Sprint 14's own selection caveat compounds it, unchanged.** `k = 30` was chosen on
**development-set RMSD** from the grid {0, 2, 10, 30, 100, 300}. The strain yardstick rules out
k = 100 and k = 300 but does not select k = 30 out of {2, 10, 30}; the non-circular variant of the
rule selects **k = 10, which emits +0.001 Å — no gain at all**. A pre-registered restraint constant
could have landed anywhere from +0.001 to −0.022, and that is before any noise-floor argument.

---

## 10. WHAT MUST CHANGE IN THE SHARED DOCUMENTS

Flagged, not edited — these are the coordinator's files.

1. **`s15/LEGACY_VS_AMBER.md` §5** currently reads: *"force-field relaxation at k = 30 is worth
   −0.022 Å [−0.036, −0.009] at valid geometry, with better Ramachandran statistics than the
   projection"* and *"AMBER's −0.022 Å is small and it is real"*. The Ramachandran half is right
   in direction and **wrong in magnitude** (0.466 → 0.874, not 0.734 → 0.874, and 116W/2L), the
   clash row of the Sprint 14 table is a placeholder, and "it is real" needs the §3.1b
   qualification: absent on 69% of the instrument, gated behind a convergence check the pipeline
   does not perform.
2. **`s15/ARCHITECTURE.md` line 209** carries the same −0.022 Å phrasing and needs the same
   treatment.
3. **`s15/coord_FINDINGS.md` R1.5** should record that the 0.081 Å floor's mechanism
   (a per-process-salted random multi-start) **cannot act on this pipeline**, which is
   bit-reproducible across interpreters, and that the pipeline's own constructed floor is
   **0.004 Å gated / up to 0.038 Å ungated** for a different reason.
4. **`s15/NEGATIVE_RESULTS.md`** should gain the convergence-gate defect (§3.1a) and the two
   Sprint 14 reporting defects (§5.1, §5.2).

---

## 11. WHAT I REFUTED, INCLUDING MY OWN

- **My own hypothesis that permuting `fit_multi`'s four fixed starts would flip optimiser basins
  and reproduce the 0.08 Å mechanism.** Exactly zero on 126 targets × 3 permutations × 2 lambda
  arms. The projection is not a source of irreproducibility here at all.
- **My own initial read that the frame null was too weak to matter.** A 3-target probe gave
  \|Δ\| ~ 1e−4 and I nearly cut the frame block to two draws on that basis. On the full 126 the
  same null contains a 1.51 Å outlier and a mean of +0.0117 Å. **A 3-target probe of a
  heavy-tailed quantity measures the wrong thing** — this is the project's "n ≤ 4 has reversed a
  conclusion four times" trap in a new place.
- **The brief's framing that this result needed "five independent start draws".** There are none.
  The replication had to be constructed on a different axis, and stating that was a precondition
  for the measurement being meaningful rather than theatre.

## 12. WHAT REMAINS OPEN

- **Does the effect survive a convergence gate on a pre-registered `k`?** The two corrections
  point opposite ways (gating raises it to −0.0234; the stability split says 69% of targets
  contribute nothing) and neither has been run as a single pre-registered protocol.
- **Is the frame-instability stratum a real target property or a property of the input structure?**
  If it is native-free and predictable, it is a selector: apply AMBER relaxation only where it
  helps. `phys_nullstruct_frame*.json` has the labels; nothing has been fitted.
- **Nothing here has been near the 60-target benchmark**, and the project's history is that no
  development improvement has transferred.

---

## REPRODUCTION

```
python -m s12.instrument                                   # pinned constants
python -m s15.phys_repl --mode canon --draws 0             # ~35 min, 126 targets
python -m s15.phys_repl --mode frame --draws 1             # exact null A, both arms
python -m s15.phys_repl --mode frame --draws 2,3 --noproj  # exact null A, AMBER arm only
python -m s15.phys_repl --mode startperm --draws 1,2,3 --noamber   # exact null B
python -m s15.phys_repl --mode boot --draws 1,2,3,4,5      # the replication axis, ~2.5 h
python -m s15.phys_repl --nullstruct --draws 1             # the convergence-gate analysis
python -m s15.phys_repl --report                           # every table above
```

Environment: `OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=2`, one process at a time,
checkpointed every 10 targets with resume. All seeding via `s15.seed.stable_rng`; no `hash()`.
