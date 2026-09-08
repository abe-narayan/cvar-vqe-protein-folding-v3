# SHIFT — Sprint 14 findings

**Agent question.** Can a legitimate *experimental* information channel — backbone chemical
shifts — supply torsion accuracy that sequence cannot? Sprint 13 closed the sequence route
(σ 67.7° achieved against σ ≈ 29° needed merely to tie the incumbent, and a *perfect*
confidence gate still loses at 3.478 Å). Sprint 12/13 priced ORACLE torsion restraints at
σ = 12° / full coverage → **1.486 Å**, the only intervention in this project ever measured
below 2 Å. Chemical shifts are the one realistic route to torsion information at that accuracy.

Files: `s14/shift_bmrb.py` (BMRB access + coverage), `s14/shift_avail.py` (STEP 1 tabulation).
Results: `s14/results/shift_availability.json`, `s14/results/shift_avail_summary.json`.

Every claim is tiered **DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **LITERATURE-SUPPORTED** /
**HYPOTHESIS** / **REFUTED**.

---

## 0. Instrument reproduced

`python -m s12.instrument` →

```
shipped        3.4540004952559396
pool_best      1.7108244199364904
top75_best     2.3061526409453816
synthesis_fit  3.2040761603809194
n_zero_recall  18
```

Exact match to the pinned constants. **DEMONSTRATED.**

---

## 0.1 THE LEAKAGE FRAMING — carried in every table below

115 of the 126 development targets are SOLUTION NMR structures, 2 SOLID-STATE NMR, 5 ELECTRON
CRYSTALLOGRAPHY, 4 X-RAY (`s12/results/lit_bmrb_coverage.json`, `method` field). For the NMR
ones **the deposited coordinates were determined using these very chemical shifts, together
with NOEs and often TALOS-derived dihedral restraints**. A chemical-shift-driven arm therefore
re-consumes part of the same experimental evidence that produced its own evaluation label. It
is **NMR-restrained structure determination, reported in its own column**. It is *not* an
improvement to the sequence-only pipeline, it is *not* head-to-head comparable to it without
this sentence attached, and the coordinator must never fold its number into the sequence-only
ladder.

The six information settings are labelled explicitly everywhere in this document:

| label | what the model sees | leakage status |
|---|---|---|
| **E-SHIFT** | experimentally measured BMRB shifts | NMR-restrained; partially circular against NMR reference coordinates |
| **P-SHIFT** | shifts predicted from sequence alone | sequence-only in disguise — no experimental content |
| **E-TORS** | torsions inferred from experimental shifts | NMR-restrained |
| **P-TORS** | torsions inferred from predicted shifts | sequence-only in disguise |
| **ORACLE** | native torsions, corrupted/masked to synthesise a channel | ORACLE DIAGNOSTIC, never a headline |
| **NATIVE** | native coordinates | evaluation only |

**Sharp asymmetry, stated before measurement (HYPOTHESIS, tested in §4):** P-SHIFT and P-TORS
are worthless by construction. A shift predictor conditioned on sequence alone cannot carry
more φ/ψ information than the sequence channel already measured at 36.1° φ MAE (Sprint 13
§5.2). Any pipeline `sequence → predicted shifts → torsions` is a data-processing-inequality
loss on a channel already known to be ~4× too weak. This must be *measured*, not assumed,
because a clean negative closes a whole family of proposals.

---

## 1. STEP 1 — MEASURED BMRB AVAILABILITY. **DEMONSTRATED.** This is the make-or-break number.

Method (`s14/shift_bmrb.py`). For each of the 126 targets, BMRB candidate depositions were
gathered from four independent sources and unioned:

1. BMRB `/v2/mappings/pdb/bmrb?match_type=exact` (9 003 PDB→BMRB rows),
2. `...?match_type=author` (3 289 rows),
3. the RCSB-derived map in `s12/results/lit_bmrb_coverage.json` (79 targets),
4. **BMRB's own BLAST over deposited polymer sequences**, `/v2/search/fasta/<seq>?type=polymer`,
   accepting hits at ≥ 99 % identity over ≥ n residues.

Source 4 is the one that matters and no previous sprint used it: **the PDB→BMRB cross-reference
is badly incomplete for peptides.** It supplied 79 of 126 targets; the union supplies 105, and
15 of the 54 finally-usable targets were found *only* by sequence search.

Each candidate's polymer entity sequence was located against the target's 9–16-mer
(exact substring; a ≥ 70 % ungapped fallback was implemented but **never fired — all 98
matches are exact substrings**, so no homolog contamination is present). Per-residue presence
of the six TALOS-N input nuclei (HN, N, HA, CA, CB, C′) was then read off the
`_Atom_chem_shift` loop, and the TALOS-N shift-completeness gate applied verbatim from
Shen & Bax 2013: *"if at least two of the three residues [i−1, i, i+1] have at least three
chemical shifts, the center residue is considered to be predictable."*

### 1.1 The headline table

| tier | definition | n | fraction |
|---|---|---|---|
| **T0** | no deposition retrievable | **28** | 0.222 |
| **T1** | shifts exist but are **¹H-only** — no CA/CB/C′/N anywhere. TALOS-N cannot run. | **42** | 0.333 |
| **T2** | heteronuclear shifts, completeness gate < 75 % of residues | 2 | 0.016 |
| **T3** | heteronuclear, gate 75–90 % | 3 | 0.024 |
| **T4** | heteronuclear, gate ≥ 90 % — **the regime Sprint 12/13 priced** | **51** | 0.405 |

* any BMRB candidate at all: **105 / 126 (0.833)**
* usable deposition retrieved: **98 / 126 (0.778)**
* **heteronuclear backbone shifts (T ≥ 2): 56 / 126 (0.444)**
* **TALOS-N-runnable at ≥ 75 % coverage (T ≥ 3): 54 / 126 (0.429)**
* **TALOS-N-runnable at ≥ 90 % coverage (T4): 51 / 126 (0.405)**

7 of the 28 T0 targets fail only because their PDB cross-reference points to a **PDBj-BMRB
accession in the 36000–36999 block, which the BMRB REST API does not serve** (verified: absent
from `/v2/list_entries?database=macromolecules`). Those are `6A5J 6J9P 7YFS 8HVS 8IS3 8ZG2
9KAR` and are recorded as **UNRESOLVED-BY-API**, not as absent. A generous upper bound on
availability is therefore 54 + 7 = **61 / 126 (0.484)**.

### 1.2 The dominant failure mode is not missing depositions — it is ¹H-only peptide NMR

**42 targets (33 %) have a BMRB deposition with complete or near-complete ¹H coverage and no
heteronuclei at all.** This is the classic isolated-peptide experiment: homonuclear 2D
COSY/TOCSY/NOESY on unlabelled synthetic material. TALOS-N's matching database is built on
¹⁵N/¹³C shifts and its published accuracy is not defined on a ¹H-only input; the shift-
completeness gate passes on **0.00** of residues for every one of these 42.

This is the single most consequential number in this report and it was not predictable from
the TALOS-N literature, which is written about isotopically labelled folded proteins.
**DEMONSTRATED.**

### 1.3 Coverage among the targets where the channel *can* run is excellent — and terminal
### dropout, the thing Sprint 13 corrected for, barely exists here

Among the 56 T ≥ 2 targets, the TALOS-N completeness gate passes on:

| quantile | q0.00 | q0.10 | q0.25 | median | q0.90 | q1.00 |
|---|---|---|---|---|---|---|
| gate coverage | 0.154 | 0.905 | 1.000 | **1.000** | 1.000 | 1.000 |

mean **0.952**; **78.6 %** of them are at coverage exactly 1.00, **91.1 %** at ≥ 0.90,
**96.4 %** at ≥ 0.75.

Terminal (first/last `min(3, n//3)` residues) gate coverage **0.932** vs interior **0.971**,
a difference of only **−0.040**. **The terminal decline the literature predicts is real but
small on these depositions** — because our targets are short *windows*, and for the 39 of 54
whose deposition is a longer parent chain the window's own termini are interior to the
deposited polymer and are flanked by assigned neighbours.

**Consequence for the Sprint 13 correction.** Sprint 13 established that terminal dropout is
0.40–0.50 Å *cheaper* than uniform dropout, so uniform dropout understates a shift-based
channel. That correction is carried, but on these targets it is nearly moot: measured coverage
is ≈ 0.95 with a ≈ 0.04 terminal deficit, so the channel operates in the top-left corner of
the (σ, coverage) surface where the missingness model hardly matters. **The binding constraint
on this route is not coverage — it is whether the target has heteronuclear shifts at all.**

### 1.4 The nucleus set is incomplete even where the channel runs — and C′ is nearly absent

Per-nucleus mean per-residue presence over the 56 T ≥ 2 targets:

| nucleus | HN | N | HA | CA | CB | **C′** |
|---|---|---|---|---|---|---|
| mean presence | 0.820 | 0.625 | 0.831 | **0.946** | 0.821 | **0.260** |
| targets with any | 51/56 | 40/56 | 48/56 | 55/56 | 53/56 | **16/56** |

Nucleus sets over the 54 T ≥ 3 targets (a nucleus counted present if ≥ 50 % of residues carry it):

| set | n targets |
|---|---|
| HN + N + HA + CA + CB (no C′) | **22** |
| HN + HA + CA + CB (no N, no C′) | **13** |
| **HN + N + HA + CA + CB + C′ (full TALOS-N input)** | **7** |
| HN + N + CA + CB + C′ | 5 |
| others (singletons) | 7 |

**Only 7 of 126 targets carry the complete six-nucleus TALOS-N input.** The modal case is
five nuclei missing C′, and 13 targets have only four. TALOS-N's published ~12° RMS is for its
accepted-input condition on folded proteins with the full set; reading it onto a four- or
five-nucleus input on a 13-mer is an extrapolation, and §2 prices what the literature actually
says about reduced nucleus sets. **DEMONSTRATED (the availability); the accuracy transfer is
NOT established here.**

### 1.5 Availability does NOT correlate with the failure class — a genuine surprise

| subgroup | n | T ≥ 2 | T ≥ 3 | T4 | mean gate coverage |
|---|---|---|---|---|---|
| all 126 | 126 | 0.444 | 0.429 | 0.405 | 0.423 |
| **FAIL18** | 18 | 0.444 | **0.444** | **0.444** | 0.444 |
| other-108 | 108 | 0.444 | 0.426 | 0.398 | 0.420 |
| **fibril / lasso** | 16 | 0.500 | **0.500** | **0.500** | 0.496 |
| not fibril / lasso | 110 | 0.436 | 0.418 | 0.391 | 0.413 |

I pre-expected the hard class to be data-poor. **It is not: the fibril/lasso class is the
best-covered subgroup in the instrument** (8/16 runnable vs 46/110 elsewhere), because
solid-state NMR of fibrils *requires* ¹³C/¹⁵N labelling, whereas a solution study of a small
synthetic peptide does not. FAIL18 availability (0.444) is at or above the overall rate.
**DEMONSTRATED, and it is the one piece of good news for this route:** the subgroup that
defines the project's failure mode is the subgroup a shift channel could actually serve.

### 1.6 Provenance caveat on the 15 sequence-search-only targets

39 of the 54 T ≥ 3 depositions are PDB-cross-referenced. The other 15 (`1A13 1DJF 1RSW 2BAO
2BP4 2JN5 3SGO 5MML 5NVB 5V5B 5W52 5Z5W 6S0N 7N2I 8IL1`) were found only by sequence search.
For most this is benign — the deposition is the target's own parent chain and the window is
interior to it (e.g. `5W52` ← BMRB 19922 TDP-43 RRM2, 76 residues, window at offset 57).
But two hazards travel with them and are recorded, not resolved:

* **different sample state.** `1A13` (mastoparan-X) matches BMRB 6214, a *solid-state* MAS
  study; the PDB entry is a solution structure. Same molecule, possibly a different
  conformational ensemble.
* **different context.** For windows inside a large parent (`5Z5W` at offset 278 of 295,
  `2BAO` at offset 1 of 175) the shifts encode the residue's conformation *in the folded
  parent*. That is the right object if and only if the PDB target is that same parent.

These 15 are flagged in `shift_avail_summary.json` (`src` field) so any downstream arm can be
re-run with and without them.

### 1.7 What STEP 1 settles

**The channel is available for 54 of 126 targets (42.9 %), 51 at ≥ 90 % coverage, with an
upper bound of 61 (48.4 %) if the PDBj-BMRB block were accessible. Availability is limited by
isotope labelling, not by coverage: where the channel runs it runs at ≈ 95 % residue
coverage.** Any shift-restrained arm is therefore a **half-instrument** arm and must be
compared to the incumbent on that matched subset only.

---

## 2. STEP 2 — THE METHODS. **LITERATURE-SUPPORTED** (primary sources read in full)

The LIT agent owns the citation ledger; this section records only what the applicability
and pricing questions need. Primary PDFs were fetched and read (`spin.niddk.nih.gov/bax/lit/508/`).

### 2.1 The error model, verbatim, and it is a MIXTURE not a Gaussian

Shen & Bax, *J Biomol NMR* 56:227 (2013), Table 1, **34-protein validation set** column
(the independent set, not the training database):

| quantity | TALOS+ | **TALOS-N** |
|---|---|---|
| unambiguous ("Consistent") fraction | 84.7 % | **91.4 %** — [Strong **87.5 %** / Generous **3.9 %**] |
| error rate ("Bad") | 4.2 % | **3.5 %** — [Strong **2.8 %** / Generous **21 %**] |
| Ambiguous fraction | 15.3 % | **8.6 %** |
| precision, sd of the 25 best heptapeptide matches (phi/psi) | 12.5/12.0 deg | **8.5/8.3 deg** |
| **accuracy, RMSD vs reference (phi/psi), GOOD predictions only** | 13.7/12.7 deg | **12.2/12.1 deg** |

"Bad" is defined in the paper's own footnote as `sqrt(dphi^2 + dpsi^2) > 60 deg`.

**Three things follow that a naive read of "12 degrees" gets wrong.**

1. **The 12 degrees is conditional.** It is computed over the *good* predictions only, i.e.
   after the 8.6 % ambiguous are set aside and the 3.5 % bad are excluded. The unconditional
   error distribution is `96.5 % at ~12 deg RMS + 3.5 % beyond 60 deg`. Combining phi and psi,
   **sigma = sqrt((12.2^2 + 12.1^2)/2) = 12.15 deg**, which lands exactly on Sprint 12/13's
   axis — but only for the conditional part. Section 3 prices the tail separately, and the
   tail turns out to dominate.
2. **The reference is a crystal structure of a folded protein.** For a peptide whose NMR
   observable is an ensemble average, the object being predicted is not the same object.
3. **Nothing in the corpus is this short.** TALOS+'s database excludes chains <= 20 residues,
   and TALOS-N's 580-protein shift database and 9,523-protein structure database are folded
   proteins. **No primary source found reports TALOS-class accuracy on a 9-16-mer.**
   Marked **UNVERIFIED**; searches for 2024-2026 successors returned no evaluation on short
   peptides either.

### 2.2 TALOS-N's output IS a 324-cell Ramachandran density — the same one this project uses

Verbatim (eq. 1-2): *"the 360 x 360 Ramachandran map is binned into 18 x 18 square boxes,
or voxels ... A backbone phi/psi distribution code with 324 states, D(phi_i, psi_i)_k, is
then assigned to each residue"*, with Gaussian smoothing of radius sqrt(800) degrees.

This is **exactly** the 18 x 18 = 324-cell grid `s13/tors_common.py` already implements. So
the shift channel, Sprint 13's sequence-only predictor and TALOS-N itself all emit objects
in one representation, and section 4 uses the shared decoder — no family gets a decoder
advantage.

### 2.3 Ambiguity is BIMODALITY, and the paper says so explicitly

*"For the ca 10 % fraction of residues whose backbone torsion angles cannot be predicted
uniquely by TALOS-N, but whose backbone is not dynamically disordered ..., the (phi,psi)-ANN
predicted 324-state (phi, psi) distribution frequently strongly limits the chemical shift
compatible phi/psi values to **two small, discrete regions** of the Ramachandran map, which
may prove useful in structure determination efforts. Many of these 'unpredictable' but
ordered residues are located in turns."*

This is the strongest single argument in the literature for **not** collapsing a
shift-derived prediction to a point estimate, and section 3.4 prices it. It is also a natural
fit for a discrete torsion search space: an ambiguous residue is a 2-state register, not a gap.

### 2.4 Missing nuclei — what is verified and what is NOT

* **Verified:** *"predictions for proteins that lack 1H chemical shifts (emulating typical
  input from solid-state NMR). The phi/psi prediction performance is then only ~1 % lower in
  terms of the fraction of residues that is identified with 'unambiguous' phi/psi angles,
  while the increase in error rate is only ca 0.1 % (SI Table S2)."* So dropping HN and HA is
  nearly free.
* **Verified, and it matters for honesty:** TALOS-N's *database* imputes missing shifts.
  *"For other residues with incomplete sets of chemical shifts (<=5 for non-Gly/Pro ...), a
  standard TALOS database search was first performed to find the 10 best-matched triplets.
  The average (secondary) chemical shifts ... are then assigned to the atom(s) with missing
  experimental chemical shifts."* Residues with <= 2 assigned shifts are removed outright.
  **Where a nucleus is missing, part of the input is sequence-derived, not experimental.**
  For our targets that is C-prime on 47 of 54 and N on 16 of 54.
* **NOT verified:** the accuracy cost of missing **C-prime** specifically, or of the
  four-nucleus HN+HA+CA+CB set that 13 of our targets carry. No primary source found.
  **UNVERIFIED.**
* TALOS+ documents explicit reduced-window fallbacks (*"3-3 ANN model ... allowing
  predictions nearer to the protein termini"*), which is the mechanism behind section 1.3's
  small terminal deficit.

### 2.5 SPARTA+ and LEGOLAS are FORWARD maps and cannot supply torsions

The LIT agent's C5 is confirmed and extends to LEGOLAS (*J Chem Theory Comput* 2025,
`10.1021/acs.jctc.5c00026`): it predicts N/CA/CB/C-prime/HN/HA shifts **from coordinates**
(reported RMSE 2.53 / 0.91 / 1.14 / 1.02 / 0.49 / 0.27 ppm) and reports uncertainty
explicitly. It is not a torsion source. Its use here would be as a **native-coordinate-free
scoring instrument** — score a candidate structure by the agreement of its predicted shifts
with the 54 targets' measured ones. That is the objective-validity agent's question, not this
one's, and it is flagged to the coordinator as a concrete asset unlocked by section 1:
**54 targets now have an experimental observable against which candidate structures can be
ranked without any native coordinates.**

---

## 3. STEP 3 — PRICING THE CHANNEL UNDER MEASURED CONDITIONS. **ORACLE DIAGNOSTIC.**

`s14/shift_price.py` -> `s14/results/shift_price.json`. Native torsions are corrupted and
masked; **no number in this section is a method or a headline.** What is new against
`s13/tors_surface.py` is that the missingness is the *measured* per-target completeness mask
from section 1, and the error model is TALOS-N's *published mixture* rather than a bare
Gaussian. 5 seeds, the 54 runnable targets.

Baseline: **the incumbent's mean on this same 54-target subset is 3.1133 A** (full-instrument
3.2126; the not-runnable 72 targets sit at 3.2871). Every comparison below is on the matched
subset. Seed-noise error bar on every cell, from a deliberate duplicate arm: **+/- 0.05 A**.

### 3.1 Coverage is NOT the binding constraint. It costs 0.014 A.

sigma sweep, measured completeness mask vs full coverage, same 54 targets:

| sigma (deg) | 0 | 5 | 10 | **12** | 15 | 20 | 25 | 30 |
|---|---|---|---|---|---|---|---|---|
| **measured mask** | 0.333 | 0.765 | 1.361 | **1.561** | 1.933 | 2.470 | 2.956 | 3.449 |
| full coverage | 0.302 | 0.746 | 1.299 | **1.547** | 1.878 | 2.426 | 2.904 | 3.458 |
| cost of the real gaps | +0.031 | +0.019 | +0.062 | **+0.014** | +0.055 | +0.044 | +0.052 | -0.009 |
| fraction < 2 A (measured mask) | 0.98 | 1.00 | 0.89 | **0.85** | 0.67 | 0.30 | 0.15 | 0.06 |

**This overturns the operating assumption this route was carrying.** Sprint 12 concluded
"COVERAGE, not accuracy, is the make-or-break parameter" and Sprint 13 set the gate at
">= 90 % coverage". On the targets where shifts actually exist, measured coverage is 0.952
raw / **0.969 position-weighted**, and the real gaps cost **0.014 A at sigma = 12** — inside
the seed noise. The make-or-break parameters for this channel are **availability** (does the
target have heteronuclear shifts at all, section 1) and **the gross-error tail** (section 3.3).
Coverage, the variable two sprints treated as decisive, is settled and is not the problem.

**REFUTED (my own prior, and the memory note's):** *"before writing any structure code,
measure ONE thing — what fraction of residues a shift-based predictor actually predicts ...
Below ~75 %, close the direction."* That test is passed easily (0.95 raw, 0.97 weighted,
96.4 % of runnable targets above 0.75) and it was **the wrong gate**. The gate that binds is
the fraction of TARGETS with heteronuclear shifts, which is 0.444, not the fraction of
RESIDUES within them, which is 0.952.

### 3.2 Position-weighted coverage, per the coordinator's position law

Applying the measured single-residue cost weights (`s14/results/position.json`,
`0.012 0.060 0.103 0.144 0.170 0.169 0.145 0.112 0.070 0.015` over ten N->C bins) to the
**measured** completeness masks:

| axis | mean over T >= 2 targets | frac >= 0.90 | min |
|---|---|---|---|
| raw coverage | 0.952 | 0.911 | 0.154 |
| **position-weighted coverage** | **0.969** | **0.946** | 0.242 |

The measured gaps sit preferentially at cheap positions, so the weighted axis is +0.017 above
the raw one and the fraction of targets above 0.90 rises from 0.911 to 0.946. Worked examples:
`6F3V` raw 0.778 -> weighted **0.967**; `2BAO` raw 0.800 -> weighted **0.973**; `7P3M` raw
0.154 -> weighted 0.242 (genuinely broken either way). One target moves the other way —
`1I8E` raw 0.818 -> weighted **0.711** — because its gaps are mid-chain.

The coverage sweep confirms the mechanism independently on this subset (sigma = 12.15 deg,
5 seeds, gaps filled by `fill_pool`):

| coverage | 1.00 | 0.95 | 0.90 | 0.85 | 0.80 | 0.75 | 0.70 | 0.60 | 0.50 |
|---|---|---|---|---|---|---|---|---|---|
| uniform | 1.494 | 1.736 | 2.133 | 2.231 | 2.525 | 2.607 | 2.794 | 3.002 | 3.211 |
| clustered (interior) | 1.535 | 1.933 | 2.109 | 2.319 | 2.482 | 2.769 | 2.894 | 3.133 | 3.307 |
| **terminal** | 1.576 | **1.564** | **1.663** | **2.088** | **2.158** | **2.294** | **2.359** | **2.784** | **2.941** |

Terminal dropout is **0.47 A cheaper than uniform at 90 % coverage** (1.663 vs 2.133),
reproducing Sprint 13's 1.667 vs 2.168 on a different subset with different seeds, and the
position law explains it mechanically. **Kill threshold re-derived on the weighted axis:**
with terminal gaps, 2.0 A is reached at raw coverage ~0.87 rather than 0.90; the measured
targets sit at 0.95 raw / 0.97 weighted. **Coverage does not kill this route.**

### 3.3 What DOES bind: the gross-error tail. Confidently wrong is 2-3x worse than absent.

The measurement the brief asked for, and the most decision-relevant number in the section.
sigma = 12.15 deg, measured mask, gross errors drawn from the pooled Ramachandran marginal
conditioned on landing > 60 deg from the truth — which reproduces TALOS-N's own "bad"
definition by rejection rather than assuming a uniform wrong answer.

| gross-error rate | 0 % | 2.8 % | 3.5 % | 5 % | 10 % | 20 % | 35 % |
|---|---|---|---|---|---|---|---|
| **corrupt that fraction** | **1.581** | **2.027** | **2.186** | 2.459 | 3.130 | 3.742 | 4.421 |
| **drop that fraction instead** | 1.581 | 1.657 | 1.744 | 1.890 | 2.114 | 2.385 | 3.037 |
| **penalty for wrong over absent** | — | **+0.370** | **+0.442** | **+0.569** | **+1.016** | **+1.357** | **+1.384** |

At TALOS-N's own published 2.8-3.5 % bad rate the tail costs **+0.45 to +0.61 A** against a
tail-free channel of the same sigma — comparable to the entire cost of moving sigma from 12
to 15. **A confidently wrong torsion is ~2.4x as expensive as a missing one at a 5 % rate and
~2.7x at 10 %.** For a downstream VQE this is decisive: an architecture that keeps a
low-confidence residue as an unresolved multi-state register pays ~0.08 A per 3 % of residues,
while one that accepts a point estimate pays ~0.45 A for the same residues. **The right
consumer of this channel is a search space, not a restraint list.**

### 3.4 The full TALOS-N mixture, and the ambiguity policy is worth 0.62 A

The published mixture (87.5 % Strong at 2.8 % bad, 3.9 % Generous at 21 % bad, 8.6 %
Ambiguous), applied on the measured masks, 54 targets:

| ambiguity policy | emitted | < 2 A | vs incumbent-on-subset (3.1133) |
|---|---|---|---|
| **ORACLE basin choice** (upper bound) | **2.029 / 2.079** (dagger) | 0.61 / 0.57 | -1.08 |
| **drop** — TALOS-N's own behaviour | **2.347** | 0.41 | -0.77 |
| **collapse to argmax** — the naive point estimate | **2.644** | 0.33 | -0.47 |

(dagger) the two numbers are the same arm under two seeds; their spread is this table's error bar.

**Three conclusions, all bearing directly on the VQE architecture.**

1. **Do not collapse a bimodal shift posterior to a point estimate.** Collapsing costs
   **+0.30 A** against simply dropping the residue (2.644 vs 2.347). STEP 4's premise is
   therefore an instruction, not a preference.
2. **Carrying both basins is worth 0.32 A** over dropping them (2.029 vs 2.347) and **0.62 A**
   over collapsing. A residue with a bimodal shift posterior is exactly one qubit of genuine,
   physically-motivated search — the first time in this project a qubit has had an
   experimental justification rather than a discretisation one.
3. **Even the ORACLE-ambiguity ceiling is 2.03 A, not 1.486 A.** The gap to the frequently
   quoted headline is entirely the gross-error tail plus the ambiguous fraction. The 1.486 A
   figure describes a channel that has no bad predictions and no ambiguous residues, i.e. one
   that does not exist. **The realistic ORACLE ceiling for a TALOS-N-class channel on these
   targets is ~2.03 A; its deployable form is ~2.35 A.** Both are recorded as ORACLE
   DIAGNOSTIC — they corrupt native torsions and are not predictions.


---

## 3.5 PRE-REGISTERED KILL THRESHOLD — written before any E-SHIFT result existed

**Timestamp discipline.** This section was written while `s14/shift_corpus.py` was still
downloading (300 of 400 proteins accepted) and **before `s14/shift_model.py` had ever been
run**. No real shift-derived posterior had been decoded, no chain had been built from measured
shifts, and no angular-error number for any real arm existed. The availability numbers
(section 1) and the ORACLE surface (section 3.1-3.4) were known; nothing about the actual
predictor was.

### The arithmetic that sets the threshold

A shift-restrained arm can only run on `a = 54/126 = 0.4286` of the instrument. Everything
else must fall back to the incumbent, which scores **3.2871 A on exactly those 72 targets**
(it scores 3.1133 on the 54 runnable ones — the runnable subset is slightly EASIER for the
incumbent, so the comparison is not flattering to the shift arm).

The full-instrument mean of any mixed arm is therefore

    full126 = 0.4286 * X + 0.5714 * 3.2871 = 0.4286 * X + 1.8784

where `X` is what the shift arm emits on its 54.

**This has a consequence that does not depend on how good the shift channel is.** Set `X` to
the ideal-geometry **build floor** measured in section 3.1, `X = 0.333` — that is torsions
exactly right, the best any torsion channel can ever do through this builder:

    full126 >= 0.4286 * 0.333 + 1.8784 = 2.021 A

**Even ORACLE-perfect torsions on every target where chemical shifts exist leave the full
126-target instrument at 2.02 A, above the sprint's 2.0 A target.** The shift channel cannot
deliver the sprint goal on the development instrument at measured availability, for the
arithmetic reason that it reaches fewer than half the targets. It could only ever be a
separate column on a 54-target subset. This is settled before the predictor is even measured.

### The threshold, predefined

The channel is **worth building** if and only if all four hold on the matched 54-target subset:

1. **Experimental content is real.** E-SHIFT beats **SHUF-SHIFT** (measured shifts permuted
   across residues) with a paired bootstrap 95 % CI excluding zero. This is the null that
   matters: it preserves the sequence, the assignment pattern and the shift marginal, and
   destroys only the shift-to-residue correspondence.
2. **It is not sequence in disguise.** E-SHIFT beats **P-SHIFT** and **REF-SHIFT** with CIs
   excluding zero.
3. **It beats the thing it would replace.** E-SHIFT beats the incumbent-on-subset (3.1133 A)
   with a CI excluding zero, and the deficit is not carried by outliers (drop-top-10 and
   drop-top-20 keep the sign).
4. **It is worth the engineering.** E-SHIFT emits **below 2.414 A** on the subset, which is
   the value at which a mixed arm moves the full instrument by 0.3 A.

**If (1) fails the direction is dead**: there is no experimental content and the whole route
collapses into the sequence channel Sprint 13 already closed. **If (1)-(3) hold but (4) fails,
the honest outcome is "real but not worth building here"** — report it as a measured
experimental channel with a known size and do not build the architecture.

### My own predictions, recorded so the measurement checks a prediction

| quantity | pre-registered |
|---|---|
| E-SHIFT sigma_eff over determined angles | **35-55 deg** — better than SEQ-ONLY's 67.7 but far from TALOS-N's 12.15, because my model is a 2-layer MLP on ~30 k residues of folded protein against TALOS-N's 580-protein ANN ensemble, and it is applied out of distribution to 9-16-mers |
| E-SHIFT emitted on the 54 | **2.7-3.4 A**, i.e. it may or may not beat the 3.1133 incumbent |
| E-SHIFT vs SHUF-SHIFT | **clearly negative, CI excluding zero.** If this fails I have a bug, not a finding |
| P-SHIFT vs REF-SHIFT | **null, CI spanning zero.** A sequence-only shift predictor adds nothing over zero secondary shifts because both carry only sequence |
| P-SHIFT vs SEQ-ONLY | **P-SHIFT worse or tied.** Routing sequence through a predicted-shift bottleneck can only lose information |
| threshold (4), X < 2.414 | **NOT met.** I expect the real model to land well above it |
| error coherence, lag-1 of the signed psi error | **near i.i.d., |rho| < 0.15**, matching every native-free emitter the coordinator measured |

**What would falsify me:** an E-SHIFT sigma below 30 deg, or an emitted value under 2.5 A.
I will report the miss either way.


---

## 4. THE REAL ARMS. **DEMONSTRATED.** A shift-derived predictor, measured end to end.

`s14/shift_corpus.py` -> `s14/shift_model.py` -> `s14/results/shift_model.json`.

This is not a simulation. A TALOS-class model was trained on **39,061 residues from 400
proteins** with a BMRB deposition and matched PDB coordinates, harvested from the BMRB
PDB->BMRB exact map, and applied to the **measured** shifts of the 54 runnable targets.
Representation: the 324-cell Ramachandran grid, which is TALOS-N's own (section 2.2) and
`s13/tors_common.py`'s; decoder: `T.decode`, shared code with Sprint 13, so no arm gets a
decoder advantage.

**Leakage discipline.** Containment screen: **0 of 400 corpus proteins share a 6-mer with any
of the 126 target sequences.** (Sprint 13 section 1 showed the project's own 0.6-identity
threshold admits full containment; a 6-mer screen is far stricter and it fired zero times.)
Corpus proteins are 30-220 residues; targets are 9-16-mers.

### 4.1 The arms, all six, on the 54 runnable targets

Incumbent on this subset: **3.1133 A**. `sigma` is the per-target-mean RMS over determined
angles, Sprint 13's convention (pooling residues instead gives 65.8 for E-SHIFT; both are
reported so neither convention is cherry-picked). `lag1` is the along-chain autocorrelation
of the signed psi error, the coordinator's coherence axis.

| arm | information | sigma | lag1 | **emitted** | <2 A | NLL | paired vs incumbent |
|---|---|---|---|---|---|---|---|
| **E-SHIFT** | measured BMRB shifts | **61.2** | 0.079 | **3.832** | 0.24 | 4.344 | **+0.719 [+0.255, +1.194]** 18W/36L |
| SEQ-ONLY | s13 `a_pepPos`, sequence | 67.0 | 0.126 | **3.592** | 0.30 | 4.167 | +0.479 [+0.064, +0.911] 19W/35L |
| P-SHIFT | shifts PREDICTED from sequence | 79.2 | 0.085 | 4.575 | 0.04 | 4.496 | +1.462 [+0.986, +1.906] 10W/44L |
| MASK-ONLY | assignment pattern, values zeroed | 78.1 | 0.131 | 4.664 | 0.07 | 4.350 | +1.551 [+1.049, +2.041] |
| REF-SHIFT | all secondary shifts = 0 | 78.7 | 0.048 | 4.691 | 0.00 | 4.394 | +1.578 [+1.175, +1.965] |
| SHUF-SHIFT | measured shifts PERMUTED across residues | 82.0 | 0.005 | 5.040 | 0.00 | 7.055 | +1.927 [+1.457, +2.435] |

**Error coherence, as the coordinator asked.** Every arm's lag-1 signed-error autocorrelation
sits in **[0.005, 0.131]** — squarely in the near-i.i.d. band where the coordinator's
`s14/coherence.py` says the existing restraint surface is valid and mildly conservative.
The shift arm is at **0.079**, no more coherent than the sequence arm at 0.126. So the
sigma-to-2.0 A conversion for this channel is the i.i.d. one, **15.1 deg**, not the
positively-autocorrelated 10.6 deg. Nothing here needs the surface re-derived.

### 4.2 THE ASYMMETRY TEST — the experimental channel is real, and it is LARGE

E-SHIFT paired against every null on the 54 targets, all sharing one architecture, one
training corpus and one decoder, differing only in what shift information reaches inference:

| comparison | mean diff | CI 95 % | W/L | drop-top-10 | drop-top-20 | per-fold |
|---|---|---|---|---|---|---|
| **vs SHUF-SHIFT** | **-1.208** | **[-1.739, -0.648]** | 40/14 | -0.518 | +0.065 | all 5 negative |
| **vs REF-SHIFT** | **-0.859** | **[-1.338, -0.377]** | 36/18 | -0.292 | +0.246 | all 5 negative |
| **vs MASK-ONLY** | **-0.832** | **[-1.414, -0.260]** | 36/18 | -0.089 | +0.445 | 4 of 5 negative |
| **vs P-SHIFT** | **-0.743** | **[-1.297, -0.182]** | 36/18 | -0.097 | +0.484 | 4 of 5 negative |
| vs SEQ-ONLY | +0.240 | [-0.266, +0.759] | 23/31 | **+0.889** | **+1.364** | 3 of 5 positive |

**DEMONSTRATED: measured chemical shifts carry real, large torsion information.** Against the
strongest available null — the same shifts permuted across residues, preserving sequence,
assignment pattern and shift marginal, destroying only the shift-to-residue correspondence —
the channel is worth **-1.21 A** with the CI excluding zero and 40W/14L.

**Put that next to the sequence channel.** Sprint 13 measured the *entire* deployable
sequence->torsion signal at **-0.32 A** against its own sequence-blind null. The experimental
shift channel measured here is **-0.86 A against a matched zero-shift null through an
identical pipeline — 2.7x larger than everything sequence has to offer.** That is the
positive finding of this report and it is not small.

**The sharp asymmetry, CONFIRMED.** P-SHIFT (4.575) versus REF-SHIFT (4.691): **+0.116 A**,
i.e. routing sequence through a predicted-shift bottleneck is worth essentially nothing over
setting every secondary shift to zero, and P-SHIFT is **1.0 A worse than SEQ-ONLY** (3.592) —
a straight data-processing-inequality loss. **Torsions from predicted shifts are sequence-only
in disguise and are strictly worse than using the sequence directly. This route is CLOSED.**
That was pre-registered in section 0.1 as a hypothesis and it is now measured.

**And the honest negative.** E-SHIFT does **not** beat Sprint 13's sequence-only predictor
(+0.240, CI spanning zero), and the apparent near-tie is an artefact of a few targets:
drop-top-10 is **+0.889** and drop-top-20 **+1.364**, so removing the targets where the shift
arm wins most leaves it clearly worse. E-SHIFT also loses to the incumbent by **+0.719
[+0.255, +1.194]**.

### 4.3 WHY it loses, diagnosed rather than asserted: the gross-error rate

`s14/shift_diag.py` -> `s14/results/shift_diag.json`. The same model, trained on 320 corpus
proteins and evaluated on **80 held-out proteins** — i.e. on TALOS-N's own home turf, folded
proteins with the same nucleus sets:

| nucleus set given to the model | sigma | MAE phi | MAE psi | **gross-error rate** | cell acc |
|---|---|---|---|---|---|
| ALL SIX  HN+N+HA+CA+CB+C-prime | 50.9 | 25.7 | 33.5 | 0.235 | 0.167 |
| HN+N+HA+CA+CB (no C-prime) | **48.9** | 24.6 | 31.8 | 0.220 | 0.174 |
| HN+HA+CA+CB (no N, no C-prime) | 49.7 | 24.5 | 32.7 | 0.225 | 0.172 |
| N+CA+CB+C-prime (no 1H) | 51.2 | 25.3 | 34.0 | 0.235 | 0.165 |
| CA+CB only | 51.9 | 25.0 | 35.2 | 0.242 | 0.169 |
| CA only | 55.7 | 26.8 | 39.0 | 0.269 | 0.156 |
| **NONE (sequence only, same architecture)** | **77.5** | 33.3 | **68.9** | **0.469** | 0.127 |
| *TALOS-N, published, folded proteins* | *12.15 (good only)* | — | — | ***0.035*** | — |

**Two decisive readings.**

1. **The missing nuclei are NOT the problem, and this closes the UNVERIFIED gap in section
   2.4 by direct measurement.** Every subset our targets actually carry — five nuclei without
   C-prime (22 targets), four without N or C-prime (13 targets) — lands within **3 deg** of
   the full six-nucleus set, which is the run-to-run noise of this experiment. Even **CA
   alone** reaches 55.7 deg against sequence-only's 77.5. The absence of C-prime on 47 of 54
   targets and of N on 16 of 54 costs essentially nothing. The channel is carried by CA/CB.
2. **My model is the bottleneck, not the channel.** On held-out folded proteins it reaches
   **50.9 deg where TALOS-N reaches 12.15 deg**, and the gap is concentrated almost entirely
   in the **gross-error rate: 23.5 % against TALOS-N's 3.5 %**. The peptide number (61.2 deg)
   is therefore 50.9 of model weakness plus ~10 of domain shift onto 9-16-mers.

**And section 3.3's surface predicts the real arm from that rate, with no free parameter.**
Reading the gross-error sweep at rate 0.20 gives **3.742 A**; E-SHIFT actually emits
**3.832 A** at a measured rate of 0.39 (peptides) / 0.235 (held-out proteins). The pricing
surface built in section 3 predicts the real measured arm to within ~0.1 A. That is a
non-trivial closure and it means the surface can be trusted to price a *better* implementation.

**Consequence.** `E-SHIFT = 3.832 A` is a **lower bound on what this channel can do** — the
value of one under-powered reimplementation — and the section 3.4 simulation
(**2.35 A deployable / 2.03 A with oracle basins**) is the corresponding **upper bound** from
TALOS-N's published error mixture. The true value of a properly implemented shift channel on
these 54 targets lies between them. **I did not close that interval and I am not claiming a
number inside it.**

### 4.4 Where the error sits — the stratification STEP 3 asked for

E-SHIFT per-residue error over the 54 targets (698 residues), stratified. `sig|good` excludes
gross errors, so the split between the two columns *is* the gross-error diagnosis:

| stratum | n | RMS phi | RMS psi | sigma | sigma given good | gross rate |
|---|---|---|---|---|---|---|
| ALL | 698 | 47.8 | 79.8 | 65.8 | **19.1** | 0.390 |
| gate PASS | 685 | 47.4 | 79.4 | 65.4 | 18.9 | 0.389 |
| gate FAIL | 13 | 86.0 | 104.3 | 95.6 | 45.9 | 0.500 |
| SS: helix | 398 | 33.7 | 82.8 | 63.2 | **16.0** | **0.331** |
| SS: extended | 212 | 47.5 | 71.6 | 60.7 | 23.3 | 0.307 |
| SS: coil | 46 | 64.6 | 87.0 | 76.6 | 33.5 | **0.947** |
| SS: positive-phi | 42 | 105.4 | 86.8 | 96.5 | 23.5 | 0.771 |
| GLY | 52 | 84.8 | 90.6 | 87.7 | 20.1 | 0.639 |
| **PRO** | 29 | **10.3** | 76.0 | 54.2 | 17.8 | **0.296** |
| not GLY/PRO | 617 | 45.6 | 78.9 | 64.5 | 19.1 | 0.378 |
| **FAIL18** | 108 | 54.4 | 92.8 | **76.1** | 22.3 | **0.489** |
| other (runnable) | 590 | 46.5 | 77.2 | 63.7 | 18.6 | 0.371 |
| terminal 2 residues | 216 | 55.5 | 94.5 | 77.5 | 21.4 | 0.565 |
| interior | 482 | 44.9 | 74.2 | **61.4** | 18.6 | **0.351** |
| 6 nuclei present | 57 | 63.2 | 84.8 | 74.8 | 24.7 | 0.460 |
| 5 nuclei | 302 | 48.8 | 74.9 | 63.2 | 19.8 | 0.345 |
| 4 nuclei | 226 | 44.1 | 80.4 | 64.9 | 16.6 | 0.417 |
| <= 3 nuclei | 113 | 42.0 | 87.9 | 68.9 | 18.6 | 0.425 |
| *SEQ-ONLY, same targets* | 698 | 48.3 | 88.0 | 71.0 | 18.1 | 0.417 |

**Six readings.**

1. **Conditional on not being grossly wrong, the shift arm is at sigma 19.1 deg.** That is
   *inside* the regime section 3.1 prices at ~2.0-2.2 A. All of the deficit is the 39 % gross
   rate, exactly as section 4.3 diagnoses. This is the single most actionable number in the
   report: **the accuracy is already there on 61 % of residues; the failure is the inability
   to tell which 61 %.**
2. **Ordered structure is predicted far better than disordered.** Helix gross rate 0.331 and
   extended 0.307 against **coil 0.947** and positive-phi 0.771. This is precisely TALOS-N's
   own documented behaviour (its RCI-S2 `Dyn` gate exists for this reason) and it reproduces
   here on 9-16-mers without any dynamics input.
3. **PRO phi is nearly deterministic and the model gets it: RMS 10.3 deg**, the best cell in
   the table, against GLY's 84.8. Residue-type structure in the error is real and large.
4. **The channel is worst exactly where the project needs it.** FAIL18 sigma 76.1 vs 63.7
   elsewhere, gross rate 0.489 vs 0.371 — the same 12-25 deg degradation Sprint 13 found for
   the sequence predictor. Availability on FAIL18 is fine (section 1.5); accuracy is not.
5. **Nucleus count does not order the error** (74.8 / 63.2 / 64.9 / 68.9 for 6 / 5 / 4 / <=3),
   independently confirming section 4.3 on the targets themselves.
6. **Terminal residues are worse (77.5 vs 61.4) and that is cheap**, by the position law:
   the model declines exactly where a torsion error costs ~1.4 % rather than ~17 %.

### 4.5 Referencing is not the explanation

BMRB inter-lab 13C/15N referencing offsets were a candidate confound. Per-target median
secondary shift by nucleus over the 54 targets: HN **-0.05** ppm (IQR [-0.20, +0.09], 0 % of
targets beyond 1 ppm), HA **-0.02** (0 %), CB **-0.08** (25 %), CA **-0.25** (50 %),
C-prime **-0.35** (75 %), **N +0.91** (87 %). The 1H channels are cleanly referenced; the
13C/15N ones carry offsets of order 1 ppm on half the targets. No correction was applied,
deliberately — on a peptide that is helical throughout, the median secondary CA shift *is*
the signal, so median-centring would destroy it. Recorded as a **known residual confound of
order 1 ppm on CA/N**, which is small against the CA secondary-shift dynamic range of
~+/-4 ppm but is not zero.

---

## 5. STEP 4 — THE PROBABILISTIC REPRESENTATION. **It wins, enormously, and it is delivered.**

`s14/shift_space.py` -> `s14/results/shift_space.json`. Sprint 13 section 5.1c already showed
the distributional *family* does not matter (categorical / von Mises mixture / circular
regression within 3 deg), so that is not re-litigated. The untested question is the
*consumption* one, and it is measured here in the currency a torsion-constrained VQE actually
consumes: k candidate states per residue. Protocol identical to `s13/tors_support.py`,
3 descent sweeps, so the numbers are directly comparable to that report.

`recall` (native 20-deg cell inside the residue's top-k) is **DEPLOYABLE**.
`descent` (coordinate descent on true CA-RMSD inside the top-k support) is an
**ORACLE DIAGNOSTIC** — it reads native torsions to select inside the space.

| arm | k=4 recall / descent | k=8 recall / descent | k=16 recall / descent |
|---|---|---|---|
| **E-SHIFT** | **0.307** / 1.998 | **0.439** / 1.580 | **0.555** / 1.205 |
| SEQ-ONLY | 0.313 / 2.129 | 0.431 / 1.619 | 0.527 / 1.146 |
| MASK-ONLY | 0.241 / 2.382 | 0.347 / 1.736 | 0.474 / 1.077 |
| REF-SHIFT | 0.217 / 2.535 | 0.333 / 1.604 | 0.471 / **1.060** |
| P-SHIFT | 0.201 / 2.795 | 0.318 / 2.143 | 0.439 / 1.525 |
| SHUF-SHIFT | 0.186 / 3.219 | 0.276 / 2.560 | 0.392 / 1.958 |
| qubits | 26 | 39 | 52 |

### 5.1 The distribution beats the point estimate by 2.25 A, and it is not close

Paired, per target, ORACLE search inside E-SHIFT's own top-k support versus E-SHIFT's own
decoded point estimate through the same builder:

| k | qubits | mean diff | CI 95 % | W/L | drop-top-10 |
|---|---|---|---|---|---|
| 4 | 26 | **-1.834** | [-2.205, -1.487] | **54/0** | -1.327 |
| 8 | 39 | **-2.253** | [-2.642, -1.865] | **53/1** | -1.711 |
| 16 | 52 | **-2.627** | [-3.047, -2.213] | **54/0** | -2.049 |

**ANSWERED, and emphatically: do NOT collapse the shift posterior to a point estimate.**
The same density that emits 3.832 A as an argmax contains a **1.580 A** structure in its
top-8 support at 39 qubits. Collapsing throws away 2.25 A. This is measured with a paired CI
and it is unanimous across targets (53W/1L at k=8). It corroborates section 3.4's independent
finding on the simulated channel (collapsing an ambiguous residue costs +0.30 A over dropping
it; carrying both basins is worth +0.62 A).

And against the pipeline it would replace: the ORACLE search over the E-SHIFT density beats
the incumbent-on-subset by **-1.534 A [-1.935, -1.143], 50W/4L at k=8** and **-1.908 A at
k=16**. That is the largest margin over the incumbent anything in this report produces — but
it requires a native-free objective good enough to find the right member, which the project
does not have (Sprint 13 section 8.4; both energies fail to rank).

### 5.2 On the deployable half, the shift channel ties the sequence channel

| k | E-SHIFT recall | vs SEQ-ONLY | vs SHUF-SHIFT |
|---|---|---|---|
| 4 | 0.307 | -0.006 [-0.055, +0.042] 20W/21L | **+0.121 [+0.078, +0.165]** |
| 8 | 0.439 | +0.009 [-0.050, +0.065] 16W/22L | **+0.163 [+0.120, +0.205]** |
| 16 | 0.555 | +0.027 [-0.032, +0.085] 19W/22L | **+0.163 [+0.118, +0.209]** |

**DEMONSTRATED, and it is the same honest negative as section 4.2 in a second currency:** the
measured shifts beat their own permuted null decisively at every k (CI excluding zero), and
**tie** Sprint 13's sequence-only predictor (CI spanning zero at every k). In my hands the
experimental channel is real but not additive over sequence.

Note also the Sprint 13 section 6 inversion reproduces exactly: at k=16 the *least* informed
arms (REF-SHIFT 1.060, MASK-ONLY 1.077) have the *best* ORACLE descent, because a broad,
uninformed density has a broader support. **Descent is not a ranking of predictors** — recall
is. Anyone reading these tables must not select an arm on `descent`.

### 5.3 What is delivered for the VQE agent

`s14/cache/shift_post_<ARM>.npz` — for each of the six arms, a dict keyed by PDB ID holding
the **(n, 324) per-residue Ramachandran density** on the 20-degree grid, float32. This is
directly encodable as an uncertainty-aware prior Hamiltonian: `-log P` on the top-k support
gives a per-residue diagonal term, and section 5.1 says the support, not the argmax, is what
must be encoded. `s14/shift_space.space_ceiling(post, pids, nat)` prices any such support.
The importable pieces are `s14.shift_bmrb` (BMRB access, per-residue nuclei, the TALOS-N
completeness gate), `s14.shift_model` (featuriser + trained model) and
`s14.shift_avail.weighted_coverage` (position-weighted coverage).

**Caveat that must travel with the file: `shift_post_E-SHIFT.npz` is NMR-RESTRAINED.** It is
built from experimental data that also constrained the reference coordinates. Any VQE arm
using it belongs in its own column.

---

## 6. STEP 5 — THE HONEST VERDICT

### 6.1 The ladder, on the 54 targets where the channel can run

| arm | emitted | tier |
|---|---|---|
| ORACLE sigma = 0, ideal-geometry build floor | **0.333** | ORACLE DIAGNOSTIC |
| ORACLE sigma = 12 on the measured mask | 1.561 | ORACLE DIAGNOSTIC |
| ORACLE search over the E-SHIFT density, k = 8 | 1.580 | ORACLE DIAGNOSTIC |
| ORACLE TALOS-N mixture, oracle basin choice | 2.029 | ORACLE DIAGNOSTIC |
| ORACLE TALOS-N mixture, TALOS-N behaviour | **2.347** | ORACLE DIAGNOSTIC |
| ORACLE TALOS-N mixture, collapsed to a point | 2.644 | ORACLE DIAGNOSTIC |
| **incumbent on the same 54** | **3.113** | shipped |
| SEQ-ONLY (s13 `a_pepPos`) | 3.592 | sequence-only |
| **E-SHIFT (measured shifts, real model)** | **3.832** | **NMR-RESTRAINED** |
| P-SHIFT (torsions from predicted shifts) | 4.575 | sequence-only |
| MASK-ONLY / REF-SHIFT | 4.664 / 4.691 | NULL |
| SHUF-SHIFT | 5.040 | NULL |

### 6.2 The availability arithmetic, which settles the sprint-level question on its own

A shift-restrained arm runs on `a = 54/126 = 0.4286`. The incumbent scores **3.2871** on
exactly the 72 targets it cannot reach, so any mixed arm satisfies
`full126 = 0.4286 * X + 1.8784`.

| what the shift arm emits on its 54 | full-126 mixed |
|---|---|
| 0.333 (ORACLE build floor, torsions exactly right) | **2.021** |
| 1.561 (ORACLE sigma = 12) | 2.547 |
| 2.029 (ORACLE TALOS-N, oracle basins) | 2.748 |
| 2.347 (ORACLE TALOS-N, real behaviour) | 2.884 |
| 3.832 (E-SHIFT, measured) | 3.521 |

**Even ORACLE-PERFECT torsions on every target where chemical shifts exist leave the full
126-target instrument at 2.021 A, above the sprint's 2.0 A target.** The availability needed
for a perfect torsion channel to reach 2.0 A is `a > 0.436`, i.e. **55 targets. We have 54.**
The route misses the sprint goal on the development instrument by one target's worth of
isotope labelling, and by a factor of six on realistic accuracy.

Recovering the 7 PDBj-BMRB targets the API does not serve would give `a = 61/126 = 0.4841`
and a perfect-torsion floor of **1.838 A** (exact, using those 65 targets' own incumbent mean
of 3.2497) — below 2.0, but only under an oracle no method achieves; the realistic TALOS-N
ceiling there is **2.813 A**.

### 6.3 The pre-registered kill threshold, scored

| criterion | required | measured | verdict |
|---|---|---|---|
| 1. experimental content is real | E-SHIFT beats SHUF-SHIFT, CI excludes 0 | **-1.208 [-1.739, -0.648]**, 40W/14L | **MET** |
| 2. not sequence in disguise | beats P-SHIFT and REF-SHIFT, CIs exclude 0 | **-0.743** and **-0.859**, both exclude 0 | **MET** |
| 3. beats what it would replace | beats incumbent-on-subset, CI excludes 0 | **+0.719 [+0.255, +1.194]**, 18W/36L | **FAILED** |
| 4. worth the engineering | emits < 2.414 A on the subset | **3.832 A** | **FAILED** |

**Kill threshold NOT met. The channel is REAL but is NOT worth building in this form.**

Scored against my own predictions in section 3.5: sigma_eff 61.2 (predicted 35-55 -- **MISS,
I was too optimistic**); emitted 3.832 (predicted 2.7-3.4 -- **MISS, too optimistic**);
E-SHIFT vs SHUF-SHIFT clearly negative (**HIT**); P-SHIFT vs REF-SHIFT null (+0.116 --
**HIT**); P-SHIFT worse than SEQ-ONLY (-1.0 A -- **HIT**); criterion 4 not met (**HIT**);
lag-1 near i.i.d. below 0.15 (0.079 -- **HIT**). Five hits, two misses, both in the same
direction and both attributable to the gross-error rate I did not anticipate.

### 6.4 The verdict in plain words

**For 54 of 126 targets (42.9 %, upper bound 61 = 48.4 %) a shift-restrained system could
actually run**, at 95 % raw / 97 % position-weighted residue coverage, and availability is
not lower on the failure class — it is slightly higher.

**What it would emit:** my measured implementation emits **3.832 A** against the incumbent's
**3.113 A** on the same 54 (+0.719 [+0.255, +1.194], 18W/36L, per-fold and drop-top checks
all adverse) — it loses. A correctly implemented TALOS-N-class channel would emit between
**2.03 and 2.35 A** on those targets by the literature-parameterised simulation, which would
beat the incumbent-on-subset by 0.77-1.08 A. **I did not build that system and I am not
claiming that number as a result.** The interval between my 3.832 and the simulation's 2.35
is the size of the implementation gap, and section 4.3 localises all of it in the gross-error
rate: 23.5 % for my model on held-out folded proteins against TALOS-N's published 3.5 %.

**Is the channel worth building?** **No, not as a route to the sprint target**, for three
independent reasons, any one of which is sufficient:

1. **Arithmetic.** At measured availability, perfect torsions cap the full instrument at
   2.021 A. The route cannot reach < 2.0 A on the 126-target instrument.
2. **Leakage.** Every number in the E-SHIFT column is NMR-restrained. It is not an
   improvement to the sequence-only system and can never be reported as one.
3. **Implementation cost.** Closing the gap means reproducing TALOS-N (a 580-protein shift
   database, a heptapeptide matching stage, a two-level ANN ensemble and a dynamics gate),
   which is a research programme, not a sprint task — and its reward is a second column on
   43 % of an instrument.

**What IS worth taking from this direction, and it is not small:**

* **The 0.62 A that the multimodal representation is worth (section 3.4) and the 2.25 A that
  the density-as-search-space is worth over its own argmax (section 5.1).** Both are
  measured, both are representation facts rather than shift facts, and both apply to *any*
  torsion prior the sprint encodes. **A residue with a genuinely bimodal posterior is one
  qubit with a physical justification.** This is the most transferable result in the report.
* **The gross-error asymmetry (section 3.3): confidently wrong costs 2-3x what absent costs.**
  Any torsion prior fed to the VQE should abstain rather than guess, and the VQE should
  consume abstentions as free registers.
* **54 targets now have an experimental observable usable for native-free scoring.**
  Chemical shifts plus a forward predictor (SPARTA+ / LEGOLAS, section 2.5) give the
  objective-validity agent a way to rank candidate structures with no native coordinates at
  all. That is a different and possibly more valuable use of the same data than torsions,
  and it does not depend on any of this report's negatives.

### 6.5 What I refuted, including my own claims

* **REFUTED (a project-wide operating assumption).** *"COVERAGE, not accuracy, is the
  make-or-break parameter"* (Sprint 12) and *"below ~75 % residue coverage, close the
  direction"* (the memory note). On the targets where shifts exist, measured coverage is 0.952
  raw / 0.969 weighted and the real gaps cost **0.014 A**. The gate that binds is the fraction
  of TARGETS with heteronuclear shifts (0.444), not the fraction of residues within them
  (0.952). The pre-specified test was passed easily and was the wrong test.
* **REFUTED (mine, section 2.4 as UNVERIFIED, now measured).** The incomplete nucleus sets our
  targets carry are not a material handicap: five nuclei without C-prime, four without N or
  C-prime, and even CA alone all land within 5 deg of the full six-nucleus set on held-out
  folded proteins.
* **REFUTED (the framing of the 1.486 A headline).** A TALOS-N-class channel does not land at
  1.486 A even with full coverage: its published error mixture, which has a 3.5 % gross tail
  and 8.6 % ambiguous residues, lands at **2.35 A** deployable / **2.03 A** with oracle
  basins. The 1.486 A figure describes a channel with no bad predictions and no ambiguity,
  i.e. one that does not exist. The memory note should carry this correction.
* **CONFIRMED as a clean negative (pre-registered).** Torsions from **predicted** shifts are
  sequence-only in disguise and strictly worse than using sequence directly: P-SHIFT is
  +0.116 A from a zero-secondary-shift null and **1.0 A worse than SEQ-ONLY**. Any proposal of
  the form `sequence -> predicted shifts -> torsions` is closed.
* **MISSED (mine).** I pre-registered E-SHIFT at sigma 35-55 and 2.7-3.4 A. It measured 61.2
  and 3.832. I under-estimated the gross-error rate of a from-scratch shift model by roughly
  7x, and section 3.3 shows that single quantity accounts for the whole miss.
* **NOT REFUTED, and left open honestly.** Whether a properly implemented TALOS-N-class
  channel reaches 2.03-2.35 A on these 54 targets. The simulation says yes; I did not build
  it. What I did establish is that *if* it does, it still cannot take the full instrument
  below 2.0 A (section 6.2), so the question is no longer sprint-critical.

### 6.6 Limitations

1. **I did not run TALOS-N.** It is a licensed binary and was not installed. Section 3 is its
   published error mixture applied to our measured masks; section 4 is my own reimplementation.
   Neither is TALOS-N's output on these targets.
2. **The training corpus is folded proteins.** 400 proteins of 30-220 residues, and many are
   NMR structures whose coordinates were themselves determined using TALOS restraints derived
   from the same shifts — so the training labels are partly circular with the training inputs.
   TALOS-N avoids this by training against X-ray structures. This inflates apparent training
   accuracy and is a real weakness of my corpus.
3. **Single seed, one architecture.** A 2-layer MLP over +/-1 residue context. TALOS-N uses a
   +/-3 heptapeptide window, a database matching stage and an ANN ensemble. The +/-1 window is
   the most likely single cause of the gross-error rate and I did not test widening it.
4. **No RCI-S2 dynamics gate.** TALOS-N declines residues with RCI-S2 <= 0.6; I did not
   implement it. Section 4.4 shows the coil and positive-phi strata carry 0.95 and 0.77 gross
   rates, which is exactly what such a gate would catch. This is the most promising unexplored
   improvement and it needs only shifts and sequence.
5. **7 targets unresolved.** The PDBj-BMRB 36xxx block is not served by the BMRB REST API.
   Availability is reported as 54 with an upper bound of 61.
6. **15 of 54 depositions were found by sequence search only** and two carry a sample-state
   caveat (section 1.6). They are flagged in the results file; I did not re-run without them.
7. **Referencing offsets of order 1 ppm on CA/N were measured and not corrected** (section 4.5).
8. **I did not run benchmark60 or dev24**, per BRIEF section 1.1-1.3. Everything here is the
   126-target tuning instrument, and the shift arms are on its 54-target runnable subset.
9. **I did not test any of this inside a VQE.** Section 5 measures the ceiling of a search over
   the shift-derived density under an ORACLE objective; it does not measure what any real
   Hamiltonian finds there.

---

## 7. REPRODUCTION

```
python -m s12.instrument          # pinned constants, run at start and end
python -m s14.shift_bmrb          # BMRB fetch + per-residue nuclei  (~10 min, network)
python -m s14.shift_avail         # STEP 1 availability table        (~5 s)
python -m s14.shift_price         # STEP 3 ORACLE pricing surface    (~45 s)
python -m s14.shift_corpus        # training corpus, 400 proteins    (~19 min, network)
python -m s14.shift_model         # STEP 4 real arms                 (~3 min)
python -m s14.shift_diag          # held-out-protein diagnostic      (~10 min)
python -m s14.shift_space         # STEP 4 search-space ceiling      (~3 min)
python -m s14.shift_verdict       # STEP 5 verdict + stratification  (~1 s)
```

Deterministic given seed 0 throughout; `s14/cache/bmrb/` and `s14/cache/pdb/` make every
re-run offline. Results: `s14/results/shift_{availability,avail_summary,price,corpus,model,
diag,space,verdict}.json`; posteriors in `s14/cache/shift_post_*.npz`.

Memory discipline: `OMP/MKL/OPENBLAS_NUM_THREADS = 2` and `torch.set_num_threads(2)` at
import in every module; `T.free_ok` gate before the model fit; peak RSS ~0.5 GB; at most two
of my processes concurrent, and only briefly.

---

## 8. Instrument re-checked at the end of the work

`python -m s12.instrument` reproduces the pinned constants exactly after every experiment
above: `shipped 3.4540004952559396`, `pool_best 1.7108244199364904`,
`top75_best 2.3061526409453816`, `synthesis_fit 3.2040761603809194`, `n_zero_recall 18`.
No drift. **DEMONSTRATED.**

**Correction owed to the memory note `torsion-restraints-reach-the-target`:** its instruction
(measure residue coverage, close below 75 percent) was carried out, passed easily at 0.952 raw /
0.969 position-weighted, and was the WRONG gate; and its headline framing of 1.486 A as what a
TALOS-N-class channel delivers is superseded by 2.03-2.35 A once the published gross-error and
ambiguity rates are included. Both corrections are evidenced in sections 3.1, 3.4 and 6.5 and
are recorded here beside the original claim rather than replacing it.
