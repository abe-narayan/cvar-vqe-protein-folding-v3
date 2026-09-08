# SPRINT 15 — VERIFY workstream findings

Adversarial audit of the sprint's write-up against the JSON artefacts in `s15/results/`.
**Report only; nothing was edited.** Every finding below names the document, the location, the
claimed value and the artefact value.

Line numbers refer to the files as of 2026-09-06 00:47 (the mtimes at the top of this audit).

---

## 0. COVERAGE

| | |
|---|---|
| numeric tokens in the four headline documents | 2,097 (dossier 157, paper 240, results 447, coord 1,253) |
| **distinct quantitative claims traced to an artefact** | **178** |
| of which **matched exactly** (to the printed precision) | **141** |
| of which **mismatched** | **21** |
| of which **could not be traced to any JSON in `s15/results/`** | **16** |
| result JSONs opened and read | 46 of 143 files in `s15/results/` |
| `s15/*.py` modules read for control/labelling hazards | 14 |
| figures inspected (code + caption + backing data) | 16 declared, **12 rendered** |

**What "could not be traced" means, itemised** — these are not necessarily wrong, they simply
have no artefact in `s15/results/`:

1. `pool MAE 2.249` (RESULTS §3, NEGATIVE_RESULTS, LEDGER 1.9) — `LEDGER.md` itself records the
   artefact for row 1.9 as `—`. No `pooldist.json` exists.
2. `corr(fit residual, true error) = −0.513` — `coherence.py` computes no such field and
   `coherence.json` has no key for it (see BLOCKER 4).
3. The three corrected fusion controls in RESULTS §6.2 / INFORMATION_CHANNEL §4 — not in
   `coherence.json` (see CORRECTION 1; I recomputed them and they are **right**).
4. `−0.949 [−1.147, −0.753]` and `−1.032 [−1.213, −0.859]` (the positive-control CIs, RESULTS §3,
   dossier §3, PAPER §7). `feasible.json` stores the point estimates but no paired CI, and
   `feasible.log` prints none. Point estimates verified; intervals unverifiable.
5. `67 papers` (dossier §7) — `LITERATURE.md` carries 64 identifier lines and states no total.
6. `15 deliverable documents` (dossier §7) — `s15/` holds 24 `.md` files; no definition of
   "deliverable" is given.
7. RESULTS §5.2's six surrogate values, §4's welsch smoke read, and dossier §5.1's ten ladder
   values — these *have* artefacts, but the artefacts now disagree (BLOCKERS 5 and 6).

---

# BLOCKERS

## B1. The 60-target benchmark — CLEAN, with one byte-level exception that should be closed

**No document and no module in `s15/` reads the protected benchmark.** I grepped every `.py` in
`s15/` and `s12/` for `benchmark60`, `bench60`, `benchmark_manifest` and every `.md` for the same.
`s12/instrument.py` reads only `s8/generate_univ/*.npz` (the 126 tuning universes) and
`bench_results/cache/1fc9f2dcf489e2fb/` (the production run **on the tuning targets**).

**The one exception, which is not a leak but should be fixed before the freeze.**
`s15/audit_provenance.py:38` declares `results/benchmark_manifest.json` as a node with the comment
*"PINNED by BRIEF; benchmark, not touched here"*, and then `build()` (lines 85–105) calls `_hash()`
on every node's files, which **opens and reads the manifest's bytes**. The emitted record
(`results/audit_provenance.json`) carries only `total_bytes: 8002`, two mtimes and
`sample_sha256_16: a40581ad01cfd2b7` — **no target identity, no sequence, no label**. So nothing is
exposed and the freeze is not invalidated. But the node's own comment is false, and a reader
checking rather than trusting (PROTOCOL_FROZEN §5) will find a read where the code says there is
none.

*Suggested fix:* add `results/benchmark_manifest.json` to a `NO_HASH` set, or change the comment to
"presence and digest only; contents never parsed", and say so in `REPRODUCIBILITY.md`.

## B2. RESULTS §6.1 presents n = 24 numbers with no n, two paragraphs before the n = 126 versions

`RESULTS.md:236–241` (§6.1 "Correlated errors, different wrong structures"):

| claimed | artefact (`coherence.json`, n = 126) |
|---|---|
| distogram fit vs native **3.270 Å** | `fusion.r_disto` = **3.622** |
| pool fit vs native **3.469 Å** | `fusion.r_pool` = **3.659** |
| the two fits disagree by **2.980 Å** | `fusion.disagreement` = **3.138** |
| raw per-pair error correlation **+0.627** | `fusion.pool_err_corr` = **+0.602** |
| pool's unrealizable/wrong ratio **0.255** | not computed at n = 126 |

Every one is the superseded 24-target read. The table carries **no n**, and `RESULTS.md:4–6`
promises that *"Sections marked (in flight) … report the smoke read with its n stated inline."*
§6.1 is not marked and states no n. §6.2, sixteen lines later, gives 3.622 / 3.659 / 3.138 / +0.602
in a column headed **n = 126**, so the document contradicts itself on the same page.

*Suggested wording:* replace §6.1's table with the n = 126 column from §6.2 and delete the
duplicated numbers, or head the table "n = 24, superseded — see §6.2".

## B3. The same stale pair is in the dossier's and the paper's headline prose

- `FINAL_DOSSIER.md:87` (§1.5): *"Two channels whose per-pair errors correlate at **+0.627**
  nonetheless build wrong structures **2.980 Å** apart."* The law quoted four lines below is the
  n = 126 law (3.205 / 3.367 / +0.162). Correct values: **+0.602** and **3.138 Å**.
- `PAPER_DRAFT.md:270` (§6): identical sentence, identical error — and the paper's **own abstract**
  (`PAPER_DRAFT.md:47–48`) already says *"correlate at +0.602 … 3.14 Å apart"*. The abstract and §6
  of the same draft state different numbers for the same two quantities.
- `INFORMATION_CHANNEL.md:186–189` carries the same four values but **does** say "on 24 targets"
  in the sentence above, so it is honest; it should still be updated.

## B4. `−0.513` is an n = 24 number sitting inside n = 126 tables

`FINAL_DOSSIER.md:60` and `PAPER_DRAFT.md:283` both place the row
`corr(fit residual, true error) | −0.513` **inside the realizability table headed n = 126**, and
`RESULTS.md:215` asserts *"The −0.513 correlation completes it"* directly after the n = 126 table.

`coherence.py` computes no residual-vs-error correlation and `coherence.json`'s `realizability`
block has no such key. The only place −0.513 exists is `coord_FINDINGS.md:604`, in the K9 section
explicitly headed *"Measured on 24 targets"*.

*Suggested wording:* either add the correlation to `coherence.py` and re-run, or move the row out of
the table and write "−0.513 at n = 24; not yet recomputed on the full instrument."

## B5. RESULTS §5.2's surrogate table does not reproduce, and its one native-free row has flipped sign

`RESULTS.md:170–176`, labelled *(n = 6 smoke read)*, against `results/errstruct.json`
(`n_done: 10`, written 00:38 — i.e. **after** the numbers in the document):

| arm | RESULTS claims | errstruct.json n = 6 | n = 8 | n = 10 |
|---|---|---|---|---|
| `real` | 2.832 | **2.899** | 2.877 | 2.960 |
| `ORACLE_shuffle_signs` | 1.825 | **1.887** | 1.866 | 1.842 |
| `ORACLE_gauss_matched` | 2.501 | **2.444** | 2.399 | 2.401 |
| `ORACLE_shuffle_pairs` | 2.577 | **2.644** | 2.505 | 2.448 |
| `ORACLE_winsor_z3` | 2.711 | **2.764** | 2.761 | 2.848 |
| **`debias_sep`** (the only native-free row) | **3.012, i.e. +0.179 WORSE** | **2.791, i.e. −0.108 BETTER** | −0.052 better | −0.065 better |

No prefix of the artefact reproduces the table. The ORACLE rows move by ≤ 0.08 Å and the
qualitative ordering survives — but **the `debias_sep` row reverses sign at every n**. The
document's derived claim, repeated in `RESULTS.md:229` (*"separation debiasing is worth little and
sometimes negative"*), `FINAL_DOSSIER.md:78`, and `PAPER_DRAFT.md:315`, is therefore not supported
by the current fit-axis artefact. (It *is* supported on the **ranking** axis — `feasible.json` gives
`ls_debias` native pct 0.359 against `ls_pred`'s 0.347 and argmin RMSD 3.559 against 3.504 — which is
what the dossier's wording actually says. RESULTS and PAPER generalise it to the fit axis, where the
artefact now says the opposite.)

**Root cause, and why this is a blocker rather than a correction.** `s15/seed.py` replaced the
per-process-salted `hash()` seeding with `blake2b`; `errstruct.py:143` now seeds with
`SD.stable_rng(p, "errstruct")`. The numbers in the document were produced *before* that fix and are
not reproducible — which is exactly retraction #4 in `FINAL_DOSSIER.md:195` and exactly what
`PROTOCOL_FROZEN.md` §4 forbids (*"exploratory numbers produced before that fix … are not quoted as
reproducible constants"*). They are still quoted, and §4's reassurance that such numbers are
"correct but not bit-reproducible" is too strong: re-running flips a sign.

*Suggested wording:* replace the table with the current artefact at its true n, restate the
`debias_sep` claim as ranking-axis-only, and delete "sometimes negative" from the fit-axis sentences
until the full run lands.

## B6. The decorrelation ladder in FINAL_DOSSIER §5.1 does not reproduce, and its most quoted sentence is no longer true

`FINAL_DOSSIER.md:235–253` (*"(n = 8; full run queued)"*) against `results/decorr.json`
(`n_done: 10`, written 00:45 — two minutes before the dossier's mtime):

| | t = 0 | t = 0.5 | t = 1.0 |
|---|---|---|---|
| dossier, constant-magnitude (`rescaled_signs`) | 2.749 | 2.337 | **2.012** |
| **artefact, n = 8** | 2.749 ✓ | 2.345 ≈ | **1.821** |
| dossier, confound (`raw_signs`) | 2.749 | **1.709** | 2.012 |
| **artefact, n = 8** | 2.749 ✓ | **1.929** | 1.821 |

Two consequences:

1. *"Destroying the coherence is worth **0.737 Å**"* (`FINAL_DOSSIER.md:240`) becomes **0.928 Å**
   (2.749 − 1.821). The inversion *"reaching 2.5 Å requires destroying only about 30% of the
   error's coherence"* becomes **27%** — same direction, different number.
2. **The confound story collapses.** `FINAL_DOSSIER.md:251–255` says *"At t = 0.5 the uncontrolled
   ladder reads 1.709 Å, better than either endpoint — geometrically impossible as a decorrelation
   effect … Reporting it alone would have produced a spectacular and wholly spurious result. This is
   the third experiment in the sprint that inverts its own conclusion without a control."*
   In the current artefact the raw ladder is **monotone** (2.749 → 2.265 → 1.929 → 1.834 → 1.821 at
   n = 8) and never dips below its own endpoint. The "geometrically impossible" reading does not
   exist in the data on disk. At n = 10 the raw t = 0.5 beats t = 1 by 0.005 Å — noise, not 0.30 Å.

*Suggested wording:* re-state §5.1 from `decorr.json` at its true n, keep the confound **control**
(it is still the right design and the confound is still large — 0.42 Å at t = 0.5), and drop the
"spectacular and wholly spurious result" paragraph, which no longer has data behind it.

## B7. `2,048 draws` is wrong everywhere: the experiment drew 4,096

`qgeom_cvar.py:216` — `def control_test(pdbs=TARGETS, alphas=ALPHAS, iters=200, seeds=(0,1,2),
budget=4096, …)`. Every cell in `qgeom_cvar.json`'s `D2_control_sub12` records
`n_drawn: 4096` for both the VQE and the control (I checked all 108 cells; the set of values is
`{4096}`). The distinct-configuration means 2.7 / 149.4 / 314.9 / 408.0 and the control's 510.3 are
out of **4,096**, not 2,048.

Occurrences to fix:

| document | line | text |
|---|---|---|
| `RESULTS.md` | 292 | "2.7 distinct configurations out of **2,048** draws" |
| `RESULTS.md` | 328 | "3 seeds, **2,048** draws, paired over 27 cells" |
| `FINAL_DOSSIER.md` | 143 | "2.7 distinct structures out of **2,048** draws" |
| `PAPER_DRAFT.md` | 189 | "**2,048** draws, paired over 27 cells per α" |
| `PAPER_DRAFT.md` | 205 | "2.7 → 408 distinct configurations out of **2,048** draws" |
| `PAPER_DRAFT.md` | 229 | "two distinct structures out of **2,048** draws" |
| `VQE_CVaR.md` | 131, 189, 206 | same |
| `coord_FINDINGS.md` | 934, 1553, 1590 | same |
| `REVIEWER_RESPONSES.md` | 136 | same |
| `qgeom_FINDINGS.md` | 59, 437, 471, 919, 942 | same |
| `s15/figures.py` | 570 | y-axis label `"distinct configs\nout of 2,048 draws"` on **fig15** |

The separate `qgeom_ens` / E2 experiment (`qgeom_ens.py`, `ens_test(..., budget=2048)`) *does* use
2,048 — and its distinct means are 2.63 / 111.6 / 232.2 / 300.9 against a control of 374, not
2.7 / 149.4 / 314.9 / 408 against 510.3. So the two experiments have been conflated: the D2
numbers are quoted with the E2 budget.

## B8. `+0.127, z = 13.1` is the **in-tail** ρ, not the global ρ — the paper says the opposite

`PAPER_DRAFT.md:179–180`: *"The generative objectives have better **global** ρ (+0.127, z = 13.1,
against −0.042) and it does not survive into the tail."*

`qrestraint_FINDINGS.md:242–245` and my own recomputation from
`qrestraint_analysed.json['targets'][*]['_skill']`:

| objective | **global** ρ (mean over 19 targets) | **in-tail** ρ (0.01) |
|---|---|---|
| `E_combined` (generative) | **+0.523** | **+0.127** (z = 13.1) |
| `E_ml_pred` | +0.513 | — |
| `E_ls_pred` | +0.449 | — |
| `E_ls_pool` | +0.316 | — |
| `S14_disto_bayes` (the selective control actually used) | +0.485 | **+0.048** |
| `legacy` | +0.105 | **−0.042** |

So the sentence attaches the tail numbers to the word "global", and it contrasts `E_combined`
against **Legacy** while the surrounding paragraph is about the generative-vs-`S14_disto_bayes`
comparison. Mean global ρ over the four generative objectives is **+0.450**, *below*
`S14_disto_bayes`'s +0.485 — so "the generative objectives have better global ρ" is also false as
stated when the selective arm is the one in the table above it.

*Suggested wording:* "In the tail, `E_combined` reaches ρ = +0.127 (z = 13.1) against +0.048 for the
selective control and −0.042 for Legacy — the only objective with measurable in-tail ordering — and
it is not the cell where the VQE wins."

## B9. The torsion-channel headline numbers are ORACLE and are not labelled

`info_FINDINGS.md:356–361` and `:401` show that **0.565**, **1.77 Å** and **−2.941** all belong to
one arm: **`ORACLE best pool window`** — the pool member selected by its RMSD to the native. The
sign-flip cost **+1.596 [+1.385, +1.813]** is that same ORACLE arm's; `info_FINDINGS.md:381` records
that for the native-free retrieval circular mean the sign-flip term is **not significant**
(−0.150 [−0.366, +0.070]).

Unlabelled appearances:

| document | line | text |
|---|---|---|
| `FINAL_DOSSIER.md` | 48–50 | "**real errors** are 0.6–2.9 Å cheaper … (alignment 0.565 …; sign-flipping costs +1.596 Å)" |
| `PAPER_DRAFT.md` | 42 (**abstract**) | "**real errors** are 0.6–2.9 Å cheaper than i.i.d." |
| `PAPER_DRAFT.md` | 299–302 | "a **real** 64°-RMS channel emits 1.77 Å … 0.565 … +1.596 [+1.385, +1.813]" |
| `RESULTS.md` | 413–415 | same |
| `ARCHITECTURE.md` | 179–181 | same |
| `INFORMATION_CHANNEL.md` | 66–74 | "the **best fragment channel**" — closer, still not labelled ORACLE |

Two separate problems:

1. **The range 0.6–2.9 Å mixes arms.** 0.6 is the native-free retrieval pool-500 circular mean
   (−0.618); 2.9 is the ORACLE window (−2.941). The largest **native-free** discount in the table is
   the incumbent's projected torsions at **−1.685**. The honest native-free range is **0.6–1.7 Å**.
2. **The comparison null is wrong.** The headline docs compare 0.565 only to the 0.945
   random-direction null. `info_FINDINGS.md:420–422` states the control explicitly: *"a **constant
   alpha-helix** already scores 0.709. Alignment below 1 is therefore not by itself evidence of a
   good channel — it must clear the helix."* A zero-information constant helix clearing the null is
   the exact trap the project's own record warns about. Only 0.565 (ORACLE) and 0.665 (incumbent
   projection) clear the helix.

*Suggested wording:* "…real torsion errors are 0.6–1.7 Å cheaper than i.i.d. errors of the same
magnitude (native-free arms), and up to 2.9 Å cheaper for the ORACLE best-matching pool window.
Alignment 0.565 (ORACLE) / 0.665 (incumbent projection) against a 0.945 random-direction null and a
0.709 constant-α-helix control; sign-flipping the ORACLE window's error costs +1.596 Å
[+1.385, +1.813], while for the native-free retrieval mean the same operation is null."

## B10. "Five times the headroom" is not derivable from any number in the sprint

`FINAL_DOSSIER.md:29–30`, `RESULTS.md:75`, `PAPER_DRAFT.md:33` and `:107`, `ARCHITECTURE.md:34` and
`:219`, `GEOMETRY.md:88`, `REVIEWER_RESPONSES.md:58`, `coord_FINDINGS.md:1108` all assert that the
value of improving distance prediction was understated **fivefold / by roughly five times**.

The two ratios the sprint's own numbers support:

- floor ratio: 1.95 / 0.611 = **3.19×**
- headroom ratio (from the incumbent 3.204): (3.204 − 0.611) / (3.204 − 1.95) = 2.593 / 1.254 =
  **2.07×**

Neither is 5. `coord_FINDINGS.md:1104–1108` gives the derivation — *"an improvement that halved
distance error could at best recover a fraction of a 1.2 Å gap. Against a 0.61 Å floor the same
improvement has roughly 2.07× the headroom (the floor itself being 3.19× lower)"* — and the arithmetic in that sentence yields 2.07,
not 5. (The only "5" available is 3.204 / 0.611 = 5.24, which is the incumbent-to-floor ratio, a
different quantity.)

This is in the **paper's abstract**. *Suggested wording:* "…moves the floor from ≈1.95 Å to 0.611 Å,
roughly tripling it and **doubling** the headroom that distance-prediction work can address."

## B11. FINAL_DOSSIER §1.7 states a quantum positive that the QENS workstream has already demoted

`FINAL_DOSSIER.md:138–145`: *"…VQE beats best-of-N by **0.36–0.57 Å** with intervals excluding zero,
at exactly the α values where diversity survives, **concentration check PASSED** — and flagged LOW
POWER … It is under replication."*

The point estimates and the cell-level CIs are correct (`qgeom_ens.json` `E2b_paired`:
rand5 α = 1.0 −0.5671 [−1.0799, −0.0679]; rand5 α = 0.05 −0.3593 [−0.7079, −0.0325]; rand75 α = 1.0
−0.4370 [−0.8357, −0.0323]), and `qgeom_nullconc.json` does verdict `DIFFUSE (PASS)` on all of them.

But `qens_FINDINGS.md` §3, written **before** the dossier's mtime and tiered **DEMONSTRATED**,
establishes from the same data that:

- the 27 "cells" are 9 targets × 3 seeds, so the unit of analysis overstates n threefold. At the
  target level (n = 9), **α = 1 is NULL on both coordinate-average readouts** — rand5
  −0.567 [−1.234, **+0.085**], rand75 −0.437 [−1.056, **+0.160**];
- the surviving α = 0.05 cell **FAILS** the null-calibrated concentration check at the target level
  (`SIG, CONCENTRATED (FAIL)`), carried by two or three of nine targets;
- `qgeom_ens.unranked_readouts`'s `set_coordavg_rmsd` has a **readout defect** (`np.unique` returns
  sorted indices, so `u[:400]` is the 400 lowest-indexed configurations, not a random subsample) —
  it bites unequally across arms;
- the budget convention gives the VQE **819,200 objective evaluations against the control's 2,048**
  — a **400× advantage** — and that is the convention the positive result was measured under.

So "concentration check PASSED" is true only at the cell level, and "intervals excluding zero" is
true only under a convention the workstream has demonstrated is wrong. The dossier's stated "reason
to disbelieve" (non-monotone in α; 2.7 distinct at α = 1) is *not* the demonstrated one.

*Suggested wording:* "…beats best-of-N by 0.36–0.57 Å at the cell level (n = 27 = 9 targets × 3
seeds). At the target level (n = 9) the α = 1 cells are **null** and the surviving α = 0.05 cell
**fails** the concentration check; the control is also charged 2,048 evaluations against the VQE's
819,200. The QENS workstream has demonstrated all three. Under replication at 19 targets × 12 seeds
with matched diversity and matched cost. **It is not a result.**"

## B12. The `n ≥ 100` figure guard is bypassed three ways, and one bypass is a headline figure

`FINAL_DOSSIER.md:279–281`: *"11 of 16 figures rendered, the remainder blocked on runs in flight and
refused by a hard `n ≥ 100` guard that prevents a smoke read from silently becoming a figure."*

`figures.py:51–68` defines `MIN_N = 100` and `_load(name, min_n=MIN_N)`. Three escapes:

1. **`figures.py:136` — explicit override.** `fig_phase()` calls `_load("distacc.json", min_n=1)`.
   `fig02_phase_diagram.png` is the **error-shape phase diagram** — the source of the headline
   `4.12 / 1.993` vs `3.00 / 2.561` comparison in dossier §1.2, RESULTS §9.2 and PAPER §5.1 — and it
   is rendered from `n_targets_phase: 40`. *Mitigation:* the panel does print "40 targets" in a
   footnote (`figures.py:177`), so it is not silent. The dossier's word "hard" is still wrong.
2. **The guard cannot see `distacc.json` anyway.** `_load` tests `d.get("n")`; `distacc.json` has no
   top-level `n` (it uses `n_targets_phase`). So `fig05_bias_profile.png` (`figures.py:244`, default
   guard) also passes without being checked.
3. **`figures.py:494–502` — `_relayed()` bypasses `_load` entirely.** `fig14_channels` and
   **`fig15_cvar_trade`** read `relayed.json` with no guard at all. fig15's data is
   **9 targets / 27 paired cells**, and the figure states no n anywhere. It is the smallest-n figure
   in the set and the one the guard was written for.

*Suggested fix:* route `_relayed()` through a guard that reads an explicit `n` from each relayed
block, print n in every figure footnote, and change the dossier sentence to "…refused by an
`n ≥ 100` guard, with two declared exceptions (fig02, n = 40; fig15, n = 9 targets), each of which
prints its n on the panel."

## B13. Two figure captions contradict their own data

1. **`fig16_gaps.png`** (`figures.py:607–608`) is titled *"Where the ångströms go / the ensemble is
   better than the baseline; **every stage after it loses ground**."* The panel itself draws
   S → A as **−0.360 Å in the "good" colour** (`figures.py:596`, `col = C_GOOD if dd < 0`).
   Aggregation *gains* ground; `RESULTS.md:36` says so ("Aggregation recovers −0.360 Å"). The
   caption asserts the opposite of the arrow drawn beneath it.
   *Suggested title:* "…selection and projection lose ground; aggregation recovers part of it."
2. **`fig05_bias_profile.png`** (`figures.py:262`) is subtitled *"a pure additive defect that had
   never been corrected."* The sprint's central result is that this defect is **not** what makes the
   channel expensive: `RESULTS.md:404–406` — *"the distogram's separation-dependent profile is
   **not** what makes it expensive"* — and `RESULTS.md:229` / `FINAL_DOSSIER.md:78` say correcting
   it is worth little. Calling it "a pure additive defect" invites exactly the correction the sprint
   measured as worthless.
   *Suggested subtitle:* "a systematic, separation-dependent offset — and §9.2 shows correcting it
   buys almost nothing."

---

# CORRECTIONS

## C1. `coherence.json` is stale relative to `coherence.py`; the mislabelled key is still on disk

`coherence.py:170–181` carries the correction comment (*"MISLABELLING CORRECTED. `np.minimum` per
target is a PER-TARGET ORACLE SELECTION"*) and now writes three keys:
`coordavg_vs_better_channel_NATIVE_FREE`, `restraintfusion_vs_better_channel_NATIVE_FREE`,
`coordavg_vs_ORACLE_per_target_pick`.

`coherence.json` (written 00:35, two minutes **before** the fix) contains **none of them**. It still
carries only the old `coordavg_vs_best_single`, whose `mean_b` is 3.1508 — the per-target minimum.
So the corrected controls have never been computed by the module.

**The numbers in the documents are nevertheless right.** I recomputed them from the stored
per-target arrays in `coherence.json['rows']` through `s12.instrument.paired` (same seed, same 4,000
resamples):

| comparison | RESULTS §6.2 / INFORMATION_CHANNEL §4 | recomputed |
|---|---|---|
| coordavg vs the better channel globally (native-free) | −0.255 [−0.370, −0.145] | **−0.2554 [−0.3696, −0.1448]**, 82/126 ✓ |
| restraint fusion vs the better channel globally | −0.071 [−0.156, +0.014] | **−0.0710 [−0.1558, +0.0142]**, 72/126 ✓ |
| coordavg vs the per-target ORACLE pick | +0.216 [+0.133, +0.302] | matches `coordavg_vs_best_single` ✓ |
| restraint fusion vs the per-target ORACLE pick | +0.400 [+0.262, +0.546] | **+0.4003 [+0.2623, +0.5455]**, 27/126 ✓ |

*Action:* re-run `python -m s15.coherence` before the freeze so the artefact matches the paper, and
so the mislabelled key name leaves the record.

## C2. RESULTS §5.4 says `2.3×` where its own table, the dossier, the paper and fig11 say `2.4×`

`RESULTS.md:211`. `coherence.json` `ratio_real_over_null = 0.41260`, so 1/0.4126 = **2.42×**.
The 2.3× is the n = 24 value (1/0.427 = 2.34, `coord_FINDINGS.md:622`) carried into the n = 126
paragraph. `figures.py:439` renders `1/r['ratio_real_over_null']` = 2.4 on fig11, so the figure and
the text disagree.

## C3. `59%` and `56%` realizable are computed with different denominators

- `coord_FINDINGS.md:626` / `INFORMATION_CHANNEL.md:124`: **56%**, defined as
  1 − (unrealizable / wrong) = 1 − 0.879/2.009 at n = 24.
- `FINAL_DOSSIER.md:66`, `PAPER_DRAFT.md:38` and `:96`: **59%**, which is 1 − 0.4126 =
  1 − (unrealizable / **the i.i.d. null**).

On the K9 definition the n = 126 value is **1 − 0.9773/2.3386 = 58.2%**. Pick one denominator and
state it; as printed, two deliverables give different percentages for the same named quantity.

## C4. `s/r = 0.88` is the n = 24 ratio; at n = 126 it is 0.86

`FINAL_DOSSIER.md:94`, `RESULTS.md:287`, `PAPER_DRAFT.md:337`, `INFORMATION_CHANNEL.md:219`,
`coord_FINDINGS.md:741`. From `coherence.json`: r = (3.6221 + 3.6592)/2 = 3.6407, s = 3.1381,
**s/r = 0.862**. (At n = 24: 2.980/3.370 = 0.884.) The conclusion is unaffected.

## C5. The selection gap is 1.34–2.14 Å, not 1.47–2.14 Å

`FINAL_DOSSIER.md:113`, `PAPER_DRAFT.md:22` and `:182`, `VQE_CVaR.md:281`,
`REVIEWER_RESPONSES.md:121` all give **1.47–2.14 Å** and describe it as *"flat across a 16× budget
range."* From `qrestraint_analysed.json['mechanism']`, over all non-ORACLE cells at all three
budgets, `selgap_vqe`/`selgap_control` span **1.341 – 2.137**. 1.47 is the minimum of the **8,192**
slice only; at 2,048 `E_combined` gives 1.369 (control 1.358) and `E_ml_pred` 1.426. The ORACLE
collapse **0.29–0.44** is exactly right (0.2915–0.4439).

## C6. `0.0016–0.0043` for the bottom-eigendecile gradient share is a selected sub-range

`FINAL_DOSSIER.md:119`, `RESULTS.md:344–345`, `PAPER_DRAFT.md:212`, `REVIEWER_RESPONSES.md:126`.
`qgeom_metric.json['A2_gradient_in_range']` has 21 cells; `grad_share_bottom_eig_decile` spans
**0.0016 – 0.2266**:

| regime | cells | bottom-decile share |
|---|---|---|
| cond ≈ 1 (flat spectrum; the decile is not meaningful) | 7 cells | 0.0455 – **0.2266** (five of them ≥ 0.09, i.e. at or above the uniform 0.10) |
| cond > 5 (where ill-conditioning is real) | 14 cells | **0.0016 – 0.0324** |

`VQE_CVaR.md:96` scopes it correctly ("where the condition number is 688–2,677"). The four headline
documents drop the scope, so as printed the claim reads as a property of the gradient in general
when for a third of the measured cells the gradient puts *more* than uniform into the small
eigendirections.
*Suggested wording:* "where the metric is genuinely ill-conditioned (cond > 5) the bottom eigenvalue
decile carries 0.0016–0.032 of the gradient's squared norm against a uniform 0.10; at cond ≈ 1 the
decile is not a meaningful split."

## C7. `median 735, max 1.5e6` is one ansatz cell, and it is not the table's maximum

`FINAL_DOSSIER.md:200`, `RESULTS.md:352–353`, `PAPER_DRAFT.md:223`. From
`qgeom_metric.json['A1_table']`, 735.4 / 1.528e6 are **`chain|L3`**. Per-cell medians across the
21 cells span 1.0 – 735.4, and the table's true maximum is **7.475e6** (`brick|L3` and
`block_chain|L3`), 5× the quoted figure. `VQE_CVaR.md:59` scopes it correctly ("for a chain at depth
L = 3"); the headline documents do not.

## C8. `Var ~ 2^(−0.18 to −0.36 n)` understates the spread

`FINAL_DOSSIER.md:118`. `qgeom_grad.json['C1_scaling_legacy']` slopes span **−0.175 to −0.364**
(so the lower end rounds to −0.17, not −0.18); `['C3_alpha']`, the CVaR-α sweep, spans **−0.137 to
−0.321**. The load-bearing part — *"every CI excluding −1.0"* — is **true in all 15 cells**
(`excludes_textbook_minus1: true` everywhere). Suggested: `2^(−0.14 to −0.36 n)`.

## C9. Three "in flight" notes are stale; their runs have landed

| document | line | says | artefact |
|---|---|---|---|
| `RESULTS.md` | 315 | "*(pool-referenced arms in flight)*" | `expand.json` has `project_scaled_pool` (3.2993, +0.095 [+0.033, +0.158]) and `project_scaled_both` (3.2844, +0.080 [+0.035, +0.129]). `PROTOCOL_FROZEN.md` §3 already records them as **FAILED**. |
| `RESULTS.md` | 288 | "*(full-instrument run in flight)*" on the fusion law | `coherence.json` is n = 126 and §6.2 already reports it |
| `FINAL_DOSSIER.md` | 281 | "In flight: the full-instrument realizability and fusion runs (upgrading n = 6–24 results to 126)" | same — landed |
| `PROTOCOL_FROZEN.md` | §3 | "restraint-level channel fusion … in flight" | `coherence.json` answers it: −0.071 [−0.156, **+0.014**] → **FAILED** under §2.1(2) |

## C10. `11 of 16 figures rendered` is now 12

`FINAL_DOSSIER.md:279–280`. `s15/figures/` holds `fig01, 02, 03, 04, 05, 06, 09, 11, 12, 14, 15,
16` = **12** (plus `info_phase_diagram.png` from `info_figure.py`, which is not in the 16-figure
series). The four still skipped are fig07 (robust), fig08 (augment), fig10 (errstruct) and fig13
(scale) — all four correctly refused because their JSONs carry `partial`.

## C11. RESULTS §4's welsch smoke read does not reproduce

`RESULTS.md:145–148` (*"an 8-target smoke read"*): welsch median **2.466** against 2.777, "three
times as many targets under 2 Å", mean **2.967** against **2.704**. From `results/robust.json`
(`n_done: 10`):

| n | squared mean / median / <2 Å | welsch mean / median / <2 Å |
|---|---|---|
| 6 | 2.791 / 2.965 / 1 | 2.739 / 2.865 / 2 |
| **8** | **2.825 / 2.928 / 1** | **2.918 / 2.865 / 2** |
| 10 | 2.773 / 2.814 / 1 | 2.994 / 2.865 / 2 |

The redescending signature survives (better median, more targets under 2 Å, worse mean) but every
number differs and "three times as many" is **twice** as many (2 against 1). Same root cause as B5.

## C12. `PAPER_DRAFT.md:229` says "two distinct structures" where every other line says 2.7

Internal inconsistency within the same draft (§3.3 at line 205 says 2.7, §3.5 at line 229 says two).
`qgeom_cvar.json` D2 gives a mean of **2.7** distinct at α = 1 (median 1). Pick one and say which
statistic it is.

## C13. The `−0.266 / −0.145` CVaR-defect comparison is an ORACLE-objective cell, and the "corrected"
arm's own interval includes zero

`FINAL_DOSSIER.md:125–130`, `PAPER_DRAFT.md:216–219`, `RESULTS.md` §8 all state *"The biased
tail-centred baseline gives −0.266 where the corrected one gives −0.145."*
`qrestraint_FINDINGS.md:583–587`: those two figures are the **`E_ORACLE_true`** column —
`vqe_TAILDEFECT_a0.25` −0.266 [−0.541, −0.024] against `vqe_const_a0.25` −0.145 [−0.388, **+0.050**].
So the comparison is (a) on an ORACLE restraint objective and (b) between a significant cell and a
**null** one. The claim "a stated liability inverts into a mechanism" needs both facts inline.

## C14. `REVIEWER_RESPONSES.md:198` quotes n = 24 numbers in the answer that promises n inline

R3.2 says *"every small-n result carries its n in the sentence that states it"* and then cites
*"real 0.879 Å versus an i.i.d. null of 2.058 Å, ratio 0.427"* with no n. Those are the 24-target
values; the full instrument gives **0.977 / 2.369 / 0.413**.

## C15. `coord_FINDINGS.md` K10's secondary result has no forward pointer to its own reversal

`coord_FINDINGS.md:752` — *"Restraint-level fusion (3.046 Å) slightly beats coordinate-level
averaging (3.105 Å)"* — is reversed at n = 126 (3.551 against 3.367), and the file records the
reversal 300 lines earlier at `:450`. The file's stated convention ("nothing has been reordered or
edited away") makes this legitimate, but K10 is the section a reader lands on from the reading map's
"K10 a native-free law for what channel fusion is worth". A one-line `[SUPERSEDED at n = 126 — see
the K10 supersession block above]` marker would cost nothing.

---

# WEAKENINGS

## W1. "Its strongest cell is the *selective* objective" needs its criterion stated

`FINAL_DOSSIER.md:110`, `PAPER_DRAFT.md:170–176`. At budget 8,192, from
`qrestraint_analysed.json['primary']`:

| arm | mean diff | wins | CI |
|---|---|---|---|
| `E_ml_pred` (generative) | **−0.2017** | 13/19 | [−0.400, −0.035] SIG |
| `E_combined` (generative) | −0.2050 | 13/19 | [−0.424, +0.003] null |
| `E_ls_pool` (generative) | −0.1840 | 14/19 | [−0.368, −0.010] SIG |
| `S14_disto_bayes` (selective) | −0.1917 | **15/19** | [−0.320, −0.076] SIG |

By **mean difference** the strongest cell is `E_ml_pred`, not the selective control. By **wins** and
by **CI tightness** the selective control is strongest. Add: "strongest by wins (15/19) and by the
tightest interval; by point estimate `E_ml_pred` is marginally larger."

## W2. "Cells significant, of 6" has the wrong denominator at two budgets

`PAPER_DRAFT.md:176–180`. `S14_disto_bayes` was only run at 8,192, so 2,048 and 32,768 have **five**
cells, not six. The counts themselves are right (1 SIG at 2,048 — `E_ls_pool`; 3 at 8,192; 0 at
32,768). Print "1 of 5 / 3 of 6 / 0 of 5".

## W3. `2.19–2.98` is called "near-constant" in support of a "provably zero" claim

`PAPER_DRAFT.md:312–314`, `FINAL_DOSSIER.md:76`. The invariance argument ("a weighted least-squares
argmin is invariant under uniform rescaling of the weights") is exactly true — for a **uniform**
rescale. `distacc.json`'s `by_separation.z_sd` is 2.186 / 2.983 / 2.699 / 2.701 / 2.671, a 36%
spread, so the miscalibration is not uniform and a separation-dependent recalibration is **not**
covered by the theorem. Suggest: "a *uniform* recalibration is provably worth zero; the measured
miscalibration is 2.19–2.98 across separation, so a separation-dependent one is not covered by that
argument — and is measured separately (it is worth little)."

## W4. The n = 6–10 smoke reads sit at or below the multi-start noise floor

`quant_startnoise.json` (n = 126, 4 seeds) records `sd_of_mean_across_seeds = 0.132 Å`,
`mean_per_target_sd = 0.449 Å`, `max_per_target_range = 5.15 Å` from the start draw alone. At n = 6–8
the standard error of a mean from start noise alone is ≈0.45/√7 ≈ 0.17 Å. Several claims fall inside
that: RESULTS §5.2's "outliers are worth 0.12", §5.2's "separation structure 0.26", §5.3's
`ORACLE_scale` −0.090 [−0.284, +0.033] (interval already includes zero and is described as a
ceiling — acceptable), and dossier §1.4's "outlier clipping / robust losses 0.12 Å (ORACLE)". Add
the noise floor to each, or drop the second decimal.

Independent confirmation of the size of this effect: the `real` control arm — identical target,
identical `n_start=4`, identical objective — differs across modules because each seeds with its own
tag (`errstruct.py:143` `stable_rng(p,"errstruct")`, `scale.py:174` `…"scale"`, `decorr.py:124`
`…"decorr"`). On the first ten targets, `errstruct.json` gives `real[1] = 4.196` where `scale.json`
and `decorr.json` give **3.642**, and `real[9] = 4.498` against **2.406**. That is by design and is
harmless within a file, but it means **no `real` baseline is comparable across the three modules**,
and any cross-module statement built from them (e.g. reading `ORACLE_scale` from `scale.json` against
`shuffle_signs` from `errstruct.json`) is invalid. I found no such statement in the headline
documents; flagging it so none is added.

## W5. `augment.py` picks its hyperparameters on RMSD without saying so

`augment.py:220–227`: `lfo()` selects `w` and `k` by `np.mean([rows[q]["arms"][…]])` where the arm
values are **native RMSDs**. It *is* leave-fold-out, so it is permitted by `PROTOCOL_FROZEN` §2.1(5),
and `quant.py:715–740` and `align_fit.py:52–55` do the same thing **with an explicit docstring note**.
`augment.py` has no such note. Add one before the augmentation arm is quoted anywhere.

## W6. The dossier implies a benchmark run that the frozen protocol says should not happen

`FINAL_DOSSIER.md:283–286`: *"…the anticipated verdict — null — was written down before **the
confirmatory run**, which is what will make the eventual number worth anything."*
`PROTOCOL_FROZEN.md` §2.6 says the opposite: *"If nothing qualifies, the benchmark is NOT run"*, and
§3 currently shows **two FAILED, four in flight/queued, zero qualified**. Either the dossier should
say "the confirmatory run, **if any intervention qualifies**", or §2.6 should be revisited — but as
written the two documents describe different futures.

## W7. `PAPER_DRAFT.md:151` "two significant effects of opposite sign and near-equal size"

0.847 and 1.092 differ by 29%. `RESULTS.md:337` and `FINAL_DOSSIER.md:132` say only "opposite sign",
which is right. Drop "near-equal size" or say "within 30%".

## W8. `fig15_cvar_trade.png` draws null and significant points identically

`figures.py:551–560` plots all four α values of `set_mean`, `set_best` and `top75_avg` as continuous
lines with no significance marking, under the title *"two significant effects of opposite sign,
both scaled by α."* From `qgeom_cvar.json` `D2b_paired`, `set_mean|a0.01` (−0.066) and
`set_best|a0.01` (+0.087) are **NULL**, as are `top75_coordavg` at α = 1.0 and α = 0.01. The caption
generalises significance across curves on which two of twelve plotted points are not significant.
Suggest open markers for null cells and "significant at α = 1.0, 0.25, 0.05; null at α = 0.01" in the
caption.

---

# CLEAN

Things I checked in full and found correct. The absence of a finding below is informative.

**The headline cascade table (RESULTS §1, dossier §2, PAPER §7) — every cell traces exactly.**
All 25 numbers in the five-row table plus the three gaps reconcile to `cascade_combined.json` to the
printed precision: G 2.7599/2.5471/0.3492/0.4921/5.1549/−0.4442 [−0.5727, −0.3247]; S
3.5108/3.3273/0.2381/0.2857/6.2170/+0.3068 [+0.1780, +0.4392]; A
3.1508/2.9457/0.2698/0.3730/5.8384/−0.0532 [−0.1246, +0.0203]; A_besthalf
3.2058/2.9624/0.2698/0.3810/5.9573/+0.0017 [−0.0863, +0.0880]; F
3.3210/3.1561/0.2619/0.3492/6.0417/+0.1169 [+0.0496, +0.1875]; gaps +0.7509 / −0.3600 / +0.1701;
`excluded_total = 2`; `n_start = 16`; `use_debias = True`; `channel = "combined"`. The ORACLE row is
labelled ORACLE in all three documents and in `figures.py:91`.

**The collapse check (RESULTS §1).** Recomputed from `cascade_combined.json['rows']['ens_spread']`:
mean **0.5271**, median **0.4771**, **17** of 126 below 0.25. Exactly as printed.

**The ceiling table (RESULTS §2, PAPER §4, dossier §1.1) — every cell traces exactly.**
`distgeo.json`: ORACLE 0.6108 / 0.0552 / 0.8571 / −2.5933 [−2.9005, −2.2946]; predicted 1/sd²
3.6445 / 3.4660 / 0.1587 / +0.4404 [+0.2896, +0.5917]; unweighted 3.7980 / 3.6248 / 0.1190 /
+0.5940 [+0.4302, +0.7607]; start 4.0725 (+0.868). The derived quantities check too: 3.03 Å between
the arms (3.6445 − 0.6108 = 3.0337), the `sd` worth 0.154 (3.798 − 3.6445), the fit improving on its
start by 0.428 (4.0725 − 3.6445). The ORACLE arm never appears without 3.644 beside it in any of the
three documents — I checked every occurrence.

**The objective-ranking table (RESULTS §3, dossier §3, PAPER §7) — all 49 cells trace exactly** to
`feasible.json['arms']`, including the machinery verification (`ORACLE_true` at percentile 0.000 on
126/126, argmin on 126/126, ρ = +0.862 positive on 100%), the 34.7th percentile, the argmin on 4 of
126 (0.031746 × 126 = 4.0), and ρ positive on 88.1%. The ORACLE row is labelled in every table and
hatched in `fig04`. Only the two bootstrap intervals are untraceable (see §0 item 4).

**Family B feasibility (RESULTS §4, dossier §2).** `feasible.json['feasibility_mean']`: 0.248 /
0.458 / 0.712 / 0.840, median z 1.336, max z 7.237, over 8,549 pairs
(`distacc.json['error_profile']['global']['n_pairs'] = 8549`). All exact.

**The distance-accuracy phase diagram and the stated requirement (RESULTS §5.1, §9.1, §9.2, dossier
§1.2, §5, PAPER §5.1).** All eight `abs_gauss` cells reconcile to `distacc.json['cells']`
(0.554 / 0.709 / 0.939 / 1.306 / 1.450 / 1.787 / 2.075 / 2.561). The effective-RMS conversion is
correct and its moment factors are the ones in the code: `sep_scaled` 2.6008 × 1.1331 = **2.947 →
2.95**, 1.6322 × 1.1331 = **1.850 → 1.85**; `outliers` 3.0 × 1.375 = **4.125 → 4.12** with emitted
1.9931 → **1.993**. `abs_gauss` thresholds 2.8743 → **2.87** and 1.8701 → **1.87**. Distogram global
RMSE **3.7042**, MAE **2.3855**, bias **+0.5092**, z-sd **2.6334** → the "2.633× over-confident" and
"+0.509 Å expansion bias" claims are exact. The "22% reduction" is 1 − 2.874/3.704 = 22.4%.

**The realizability table (dossier §1.3, RESULTS §5.4, PAPER §5.2)** — the five residuals, the mean
|error| and the ratio all trace exactly to `coherence.json['realizability']`: 0.9773 / 0.4581 /
2.3688 / 2.0810 / 2.0700 / 2.3386 / 0.41260, paired −1.3914 [−1.6234, −1.1836] with 126/126 better.
"Within 13% of the null" checks: 2.0810/2.3688 = 0.879 and 2.0700/2.3688 = 0.874. Every surrogate row
is labelled ORACLE in all three documents. (The `2.3×` and `−0.513` defects are C2 and B4.)

**The fusion law at n = 126 (RESULTS §6.2, PAPER §6, INFORMATION_CHANNEL §4).** r_disto 3.6221,
r_pool 3.6592, r = 3.6407 → **3.641**; s 3.1381 → **3.138**; law 3.2049 → **3.205**; measured 3.3668
→ **3.367**; error 0.16189 → **+0.162**; restraint fusion 3.5511 → **3.551**; correlation +0.6021 →
**+0.602**. All exact.

> **⚠ SPRINT 16 NOTE — 2026-09-06, RETRACT workstream. This row passed, and it should not have
> reassured anyone.** Every number above traces to `coherence.json` exactly, which is what this audit
> checked. But `coherence.json` itself computed `r` as the **arithmetic** mean of the two channels'
> RMSDs where the identity requires the **quadratic** mean — so the audit verified the faithful
> transcription of a wrong number. Corrected: law **3.292**, error **+0.075**
> (`s16/retract_law.py`, n = 126, paired −0.087 [−0.119, −0.061] i.i.d.-target,
> [−0.129, −0.056] fold-clustered, W/L 126/0).
>
> **The general lesson, recorded here because this audit is where it would have been caught:**
> **a transcription audit cannot see a formula error.** Checking that a document's numbers trace to
> an artefact is orthogonal to checking that the artefact computes what its own docstring says. The
> Sprint 15 verification pass had no arm that re-derived a published identity from its definition.
> Sprint 16 added one only because the literature workstream read the formula, not the numbers. The sign-flip correction narrative in §6.2 is correct and the ORACLE
per-target pick is labelled "*a ceiling, not a control*" in every row.

**The projection gap (RESULTS §7).** `expand.json`: project_raw 3.2041; project_scaled 3.3273
(+0.1232 [+0.0653, +0.1835]); project_scaled_deb 3.2889 (+0.0849 [+0.0404, +0.1315]);
`ORACLE_scale_true` 2.9287 (−0.2754 [−0.3972, −0.1740]). The scales: truth 1.04116 → **1.0412**,
distogram 1.09968 → **1.0997**. The 25.8% withdrawal is correct and complete: measured
`contraction_vs_true_pct` **3.5046** and `contraction_vs_pred_pct` **9.5422**. I grepped all 24 `.md`
files — **every surviving occurrence of 25.8% is inside an explicit withdrawal**
(`BRIEF.md:102–106`, `coord_FINDINGS.md:49/171/769/777/781/783/927`, `FINAL_DOSSIER.md:194`,
`GEOMETRY.md:106`, `RESULTS.md:317`, `REVIEWER_RESPONSES.md:206`). Nothing stands unqualified.
`s15/figures.py:371` renders the corrected 3.5% on fig09.

**The CVaR readout table (RESULTS §8, PAPER §3.2b) — all 20 cells trace exactly** to
`qgeom_cvar.json['D2b_paired']`, including every SIG/null verdict: argmin −0.0362/+0.2507/+0.0682/
+0.0480; top20 coordavg +0.0676/+0.3723/+0.1973/+0.1772; top75 coordavg +0.1467/+0.3920/+0.2764/
+0.1640; set mean −0.8468/−0.3719/−0.1910/−0.0661; set best +1.0924/+0.3446/+0.1552/+0.0865. The
monotone-shrinkage claim holds in both series. The distinct-configuration ladder 2.7 / 149.4 / 314.9
/ 408.0 and control 510.3 are exactly right (only the *denominator* is wrong — B7).

**The classical-control claim.** *"0 wins in 16 against greedy, losing 3"* is exact: 16 `vs_greedy`
cells in `qrestraint_analysed.json['vs_classical']`, zero with a significant negative difference,
three with a significant positive one (`E_ls_pool|2048` +0.192, `E_ls_pool|8192` +0.153,
`E_ORACLE_true|2048` +0.270).

**The budget-dependence claim.** `E_ml_pred|8192` = −0.2017, 13/19, CI [−0.400, −0.035] ✓;
0 significant cells at 32,768 ✓ (all five CIs cross zero).

**"tail − bulk ordering skill is negative for every objective, including all four generative ones and
the ORACLE."** I recomputed the mean of `rho_tail_0.01 − rho_global` over all 19 targets for all nine
objectives in `qrestraint_analysed.json['targets']`: −0.374 / −0.409 / −0.272 / −0.396 / −0.419 /
−0.147 / −0.014 / −0.437 / −0.178. All nine negative. Exactly as claimed.

**Estimator and instrument verification.** `audit_cvar.json`
`max_abs_value_error_vs_Rockafellar_Uryasev = 4.1744e-14` → **4.2e-14** ✓; the second recorded defect
is inert (`median_abs_grad_error_fd_vs_envelope = 0.0`) ✓. `qgeom_metric.json['verify']`
`max_abs_classicalFIM_minus_4g = 3.3307e-16` → **3.3e-16** ✓. `grad_share_in_null_space` spans
0.0 – 1.13e-30 → **"0 to 1e-30"** ✓. The `block` ansatz is rank-deficient at depth ≥ 2
(`A2`: P = 20 → rank 15 at L2, P = 30 → rank 15 at L3) ✓. `audit_rmsd.json`: Horn vs the production
Kabsch **2.0423e-13 Å** over **63,000** structures ✓, with mirror-image and point-inversion edge
cases both returning 2.960 where a sign-unfixed SVD returns 0 ✓.

**Sensitivity constants.** `audit_sensitivity.json`: `n_zero_recall` **45 → 2** across BAND 0.5 → 3.0
with 18 at BAND 1.5 ✓ (RESULTS §1's FAIL18 footnote). The oracle-conditioning transposition
(**21.7%** oracle, **40%** unbiased) matches `audit_FINDINGS.md:499–506` ✓.

**Literature figures.** 4.89 Å mean / 4.70 median / sd 1.10 on 75 fragments (Zhang et al.,
`LITERATURE.md:204`) ✓; AlphaFold2 rank-0 at **13%** against a **20%** uniform null (Gulsevin &
Meiler, `LITERATURE.md:163/229/514`) ✓. `FINAL_DOSSIER.md:270`'s "not a like-for-like win" caveat is
present and correct.

**ORACLE labelling in code.** I read every `min()`/`np.minimum` in `s15/*.py` looking for a second
per-target-minimum mislabelled as a control. `coherence.py:170` is the one already found and
corrected. The others are all correct: `cascade.py:139` `G = float(rms.min())  # ORACLE`;
`expand.py:148` feeds `unprojected_gain_ORACLE` / `oracle_best_scale`, both ORACLE-named;
`augment.py:166` `fit_best_ORACLE`; `feasible.py:138`, `info_*.py`, `qgeom_*.py`, `qens_lib.py:335`
are all `pool_best` / `set_best` evaluation statistics, correctly named and never used as controls.
`qgeom_qng.py:249/282` and `qgeom_qngreport.py:65` take a minimum over **arms on a global mean**,
which is a winner-declaration, not a per-target oracle. **No second instance of the `coherence.py`
mislabelling exists.**

**Hyperparameters chosen by RMSD.** Three modules do it, all leave-fold-out and therefore permitted
by `PROTOCOL_FROZEN` §2.1(5): `quant.py:715–740` (documented in the artefact's own `note` field),
`align_fit.py:50–55` and `:276` (documented in the docstring, with the in-fold ceiling labelled
ORACLE), and `augment.py:220` (**not** documented — W5). No module selects a multi-start by RMSD:
`distgeo.py:158`, `errstruct.py:142`, `scale.py:174`, `decorr.py:124` and `coherence.py` all select
on the objective, and `distgeo.py`'s docstring states it.

**Figure-by-figure captions against data.** fig01 (title and every CI matches
`cascade_combined.json`; ORACLE hatched and labelled), fig03 (median z and max z interpolated live
from the artefact), fig04 (arms sorted, ORACLE hatched, chance line at 0.5), fig06 ("worth 0.611 Å,
not 1.95 Å" matches `distgeo.json`; the 1.95 reference rule is labelled "through the library"),
fig09 (3.5% — the corrected figure), fig11 (2.4× computed live from the ratio), fig12 (n and MAE
read live from `coherence.json`, n = 126), fig14 ("the best is 0.566 against a requirement of 0.638"
matches `relayed.json['channels']`) — **all correct**. Only fig05 and fig16 have caption defects
(B13) and only fig02/fig15 have guard defects (B12).

**No figure is rendered from a `partial` artefact.** `_load`'s `"partial" in d` test correctly
refuses `robust.json`, `augment.json`, `errstruct.json` and `scale.json` — the four unrendered
figures — so the four smoke-read tables that are wrong in the text (B5, B6, C11) did **not** become
figures. The guard did the job it was written for; the three escapes in B12 are all in the other
direction.
