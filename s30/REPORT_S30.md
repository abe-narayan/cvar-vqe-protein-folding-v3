# Sprint 30 — CVaR-VQE Protein Folding: the search for the first real accuracy breakthrough

**Status: DRAFT IN PROGRESS.** Sections marked `[PENDING]` await lanes still running. The charter
requires the report only after everything is finished; this file is assembled as results land so
nothing is reconstructed from memory at the end.

Branch `s26` · benchmark: 126 dev targets, 9–16 aa · endpoint: **mean built-chain Cα RMSD**
Ledger: `s30/LEDGER.md` · Running state: `s30/STATE.md` · Contract: `s30/S30_CONTRACT.md` (28 rules)
Charter: `s30/BRIEF.md`

---

## 0. The answer, up front

[PENDING — written last, once lane P's decisive measurement lands.]

The charter asked for the first real accuracy breakthrough, with a primary target of < 3.0 Å and an
ambitious target of < 2.5 Å against production's **3.2105 Å**. It also said, explicitly, that *"a
well-evidenced ceiling argument that redirects the next three sprints is worth more than a fragile
2.98 Å."*

---

## 1. The arithmetic that governed the sprint

Computed at the open, before any experiment, from S29's stored rows:

```
production, n = 126      mean 3.2105   median 2.9661
the worst 18 targets     mean 6.2758
the other 108            mean 2.6997

cap the worst 10 at 3.00 Å  ->  mean 2.9074  (−0.3031)   beats the primary target
cap the worst 18 at 3.00 Å  ->  mean 2.7426  (−0.4680)
cap the worst 30 at 3.00 Å  ->  mean 2.5778  (−0.6328)

versus improving EVERY ONE of the 126 by 0.20 Å  ->  mean 3.0105  (−0.2000)
```

**Fixing ten targets is worth more than improving all 126 by 0.20 Å.** The endpoint is a mean over
a distribution whose median already clears 3.0, so a mechanism that works on typical targets and
leaves the tail alone is close to worthless here. This is why the sprint staffed the tail first.

The caveat carried with it throughout: capping is an ORACLE operation. It bounds the prize; it does
not deliver it.

---

## 2. What was tested

[PENDING — the eight lanes, their pre-registrations, and the multiplicity ledger.]

---

## 3. What was closed, and by what

The sprint's dominant output. Each row is a direction that is now shut, with the instrument that
shut it. Several are theorems rather than measurements.

| direction | closed by | how |
|---|---|---|
| combining the 21 native-free displacement fields | D + T | ORACLE global ρ = **0.169**, LFO **0.012** (worse than the best single field) against **0.358** needed for 3.00 Å; Gram stable rank 2.057 — the fields span ~2 directions, not 11 |
| sparse weighted readout (S29's last open ladder class) | Q | **by price** — a plain argmin dominates at every bit budget; 2-of-75 with free weights is 2.1683 Å at 11.4+ bits against argmin-over-128 at **2.1435 Å for 7.0 bits** |
| second-moment / quadric escape | Q **and** T, independently | +0.2059 Å (Q, 1.81× MDE) and +0.086 (T, 2.95×); the structured `disp2` form scores 3.3585 against production's 3.0483 |
| subset objective through an averaging readout | T | **T1**: the tail is *always* a prefix — of the order induced by ∇V at the optimum. S29 §4.3 and §4.5 are incompatible and §4.3 wins |
| generative structural spaces | X | closed **jointly with the readout**: at ρ≈0 the endpoint is near-unit-slope in the set **mean**, and width buys ceiling while costing mean |
| torsion / configuration encodings | T | arithmetically infeasible — **48 bits** for 4 basins/residue at n=12 against **7** deployed |
| common-mode correction from pool data | L | **non-identifiable**: the likelihood depends on (t, μ) only through t+μ at any K; **m_eff = 1.4** of 75 members |
| E1 — a prior on the shared bias's form | L | three independent ways |
| E2 — constraint repair | L | field-scale: removing 98% of clashes costs **+0.08 Å**; works here only because 2/n is 15.4% at n=13 |
| filter width as a free lunch | F | no knee — benefit/harm cross smoothly at ~55%/55%; the **ORACLE global argmin over k IS the shipped 75** |

---

## 4. What was established positively

### 4.1 The tail is selection-limited, not pool-limited

The ORACLE best member of the 18 worst pools is **2.2842 Å** (built chain 2.2845) — already under
the 3.00 Å cap the opening arithmetic asks for, with **13 of the 18** holding a member under 3.00
and the worst tail pool bottoming out at 3.54. **The material to fix the tail is already inside the
candidate sets the pipeline is handed.** Nothing in this claim is conditioned on the thing it
measures, which is why it survived the adversary while the entry's headline did not (§Appendix A).

Retrieval is exonerated at every stratum: BLOSUM's 500 against a random 500 of the same universe is
NOT MEASURED everywhere, most pointedly on the tail (0.01× its own MDE). **This revises what the
project had recorded** — the harm on hard targets was attributed to the retrieval corpus; it is not
retrieval.

**The strongest tail result, and it replicates on all three tail definitions:** the shipped score's
Spearman with ORACLE in-pool RMSD is **+0.6446 on the easy 108** (fold CI [+0.5905, +0.7063], 5/5
folds) and **+0.1066 on the hard 18 with the fold CI including zero**. *The score cannot order its
own pool on hard targets.* It is the distogram's own error that drives this (ρ = −0.799, −0.819
length-residualised), the chain is fully mediated, and **shape error is 83% of it — scale is
refuted as the mechanism** (partial ρ = −0.067, p = 0.46, against shape's −0.641 at p = 6.5e−16).

### 4.2 The pool is a codebook, not a channel

The charter asked where the 1.44 usable bits went and where the other 5.56 were spent. The question
is not well-posed: **the 500 deposited backbones carry the structure and the index only names it**,
so bits are not conserved across an index. Seven index bits move the ORACLE ladder 4.108 → 1.898 Å,
which through the displacement bound is ρ = 0.887 — **36.6 bits of displacement information out of
7 index bits, a 5.2× ratio**.

> The honest inversion: **the readout's 7 bits are worth five times their face value, and the
> system cannot supply even one of them.**

The **value-of-a-bit law** follows and makes allocations comparable: `D(R) = a + c·2^(−R/γ)` fits at
**R² = 0.9983** (a = 1.3312 Å, γ = 3.1636), so `−dD/dR = 0.219·(D − 1.331)` Å per bit. Candidate
indexing beats subset cardinality by **3×** — which *explains* the 3.5× S29 measured — and
**torsion/configuration encodings are arithmetically infeasible** at this width (48 bits for four
basins per residue at n = 12, against 7 deployed). The charter called the encoding the least-examined
component and the likely hidden bottleneck. **It was examined and it is not.**

### 4.3 The field library spends its rank on the wrong direction

Lane L derived that a fixed-reference statistical potential contains a **separable term that is a
pure function of scale** (`+kT·Σ ln P_ref(d_ij)`; a uniform 10% contraction of a 13-mer with zero
shape change moves it ~14 kT, while a size-matched reference moves by exactly 0.0000). It predicted,
falsifiably and with its concession pre-stated, that **one of the Gram's two effective directions
must therefore be the radial one**. Measured independently:

```
radial share of the Gram trace                        0.5798  [+0.545, +0.603]
cos(dominant principal direction, radial)             0.947 mean / 0.984 median
fraction of lambda_1 that is radial                   94.8%
stable rank with radial removed                       1.705 -> 2.642

cos(direction to the native, radial), all targets    -0.0675
cos(direction to the native, radial), FAIL18         -0.2524
```

> **The library spends the majority of its two available directions on a component that is
> orthogonal to the answer in general and *anti-aligned* on the hard targets** — and lane F had
> independently shown that shape, not scale, is 83% of the tail's error. Three lanes, one chain.

It does **not** follow that deflating scale helps: the residual 42% has no large eigenvalue, the
combination's ceiling already prices the whole span, and lane L withdrew its own deployable proposal
on exactly this point — the size-matched field *is* the deflated field, and deflation reallocates
rank without creating any.

### 4.4 Ordering exists; preference does not

On a ladder where **kind, local realism and perturbation budget are all matched** — every rung an
ideal-geometry backbone built from torsions drawn from the fold's leakage-safe Ramachandran table,
no projection, no averaging, no contraction — the verdict splits:

- **Ordering survives.** Two of 43 channels clear the bar, DIS at +0.347 with an anchor contrast of
  +0.134 (max-over-channels sign-flip null p = 0.000).
- **Preference fails on all 43.** The best `pref_near` in the library is 0.640, under the bar, and
  that channel prefers a *random pool member* to production on 0.790.
- **The largest preference effect in the library points the wrong way:** DIS prefers production to a
  0.55 Å structure on **94.2%** of targets. The project's most-cited negative, with its confound
  removed, comes out **sharper**.
- The leave-fold-out combination prefers the near-native rung to production on **93.0%** of held-out
  targets — and prefers an **arbitrary pool member on 100%** and a **3 Å rung on 100%**. Margin
  **−0.070 [−0.110, −0.028]**. *It learned "is this production?", not "is this near-native".*

**And the mechanism makes the null a theorem on this instrument.** Held-out R² at matched capacity:
local features give **ΔR² −0.089**, global features **+0.600** — and the local block is
**ORACLE-advantaged**, handed per-residue deviations from the native anchor, and still adds nothing
beyond knowing the perturbation budget. **A sum of per-residue terms cannot see a lever arm.**

Resolution, which bounds what any support rule could ever select on: DIS concordance inside the near
band is **0.520** at |Δ| = 0–0.25 Å, reaching 0.822 only above 4 Å. **Coarse triage and nothing
else.**

### 4.5 Generation is closed jointly with the readout

For a coordinate-average terminal, `set_mean² ≈ B² + S²` where **B is the endpoint itself** and S is
the set's spread. The terminal's entire value is spread extraction (`avg_gain = 0.4143·S − 0.1559`,
r = 0.885), and within target `corr(S, B) = +0.085` — **concentration is orthogonal to bias**. So a
set mean improved purely by concentration is worth **zero by algebra**, and the half that pays *is*
the endpoint, which is non-identifiable from pool data at any K.

The extreme case settles it: the most concentrated source ever built here (spread 0.568 against the
pool's 1.577) is **the worst endpoint in the record, 3.789**. **For an averaging terminal, spread is
the raw material, not a defect.**

---

## 5. What the CVaR-VQE contributed

[PENDING]

## 6. Statistical discipline and multiplicity

[PENDING]

## 7. Every hypothesis entertained and killed

[PENDING]

## 8. The twelve S29 leads: pursued, rejected, and why

[PENDING]

## 9. Literature relied on and rejected

[PENDING]

## 10. The cost/RMSD meter's baselines for every cost tested

[PENDING]

## 11. What remains open

[PENDING]

## 12. The next highest-value scientific question

[PENDING]

## 13. Architecture diagram of what was actually built

[PENDING]

---

## Appendix A — claims withdrawn during this sprint

The sprint's most distinctive feature was lanes destroying their own results, usually before anyone
asked. Recorded in full because the charter requires it and because it is the reason the surviving
results are worth anything.

### A.1 The coordinator's (mine)

| claim | how it died |
|---|---|
| **The "two qubits" synthesis** — a selecting readout over a wider register aimed at the tail | Written at 13:07, corrected at 13:08. Lane F closed filter width by ceiling (the ORACLE global argmin over k **is** the shipped 75) and showed the widening-rescues-the-tail result does not replicate on either filter-independent tail (≈0, and *reversed* at +0.73 on one) |
| **"The median target is worse than chance"** — reported to the user twice | `7 − log₂r` is right-skewed. Null mean 1.4050, null **median 0.9888**, null win rate 62.5%. "Below random on 82 of 126" is what a random ranking does to itself (null expects 78.8, z = +0.60). I read a median against a mean's baseline |
| **"LEG_torsion is at chance"** | My wording, not lane R's. Its anchor contrast is **+0.024 with the fold CI excluding zero** — it *does* order slightly above its control. It fails by being a quarter of the required margin, not by being noise |
| **The ceiling gate I handed lane X** — "measure the ceiling first, always" | An **anti-predictor**: 1 of 10 cells correct on predicting the endpoint's direction, against 5/6 for the set-mean law; mean residual 0.333 Å vs 0.0153 — **22×** |
| **"Nobody has built a typical-good generator"** | Somebody had, and **diversity is the correct design**. For an averaging terminal spread is the raw material; the most concentrated source ever built is the worst endpoint in the record |
| **"Run `selftest` and confirm it reproduces the baselines"** | `selftest` is a synthetic 8-residue check that runs in 0.4 s and **would pass on a meter whose every number had drifted**. The instruction was unsatisfiable. A real `verify` now exists |
| **The 0.3688 anchor in lane D's brief** | That is **DIS_SURR's**; the shipped cost is 0.3676. I listed it beside five `DIS` numbers |
| **"Find a source with decorrelated errors"** | Priced by lane L: a *perfectly orthogonal* channel is worth a **5.1% discount** on the requirement and must be **3.01× better** than anything owned. **Decorrelation is not the lever; skill is** |
| **Quoting stable rank 1.86 without its feature space** | 1.86 is **pair-distance** space; **coordinate space is 3.4–3.6** with k90 = 11.2 not 5.6 — and sparse readouts and second-moment constructions act on coordinates |
| **Citing S28-L48 as established** | Every rung in that ladder differs in **kind** as well as nativeness, so it measured a cross-kind preference that perception–distortion already predicts (lane R, S30-L1, posted before its own numbers existed) |

### A.2 Lane T's — four, three toward its own hypothesis

| claim | how it died |
|---|---|
| **P2b: retrieval worth +0.10–0.30 Å** | Measured **0.07**; its copula reasoning overestimated by 2–4× |
| **P4c: searched classes beat prefix by < 0.6 Å** | Measured **0.982** — missed by 64%, in the flattering direction |
| **"12 bits ≈ six coefficients — two independent routes"** | An **unregistered aside**; a coincidence, not a derivation. Withdrawn and corrected *downward* to 3.78 bits, and flagged louder precisely because it had no bar attached |
| **"HALFSPACE − PREFIX = −0.9816 Å, BETTER"** — its own headline | Lane Q named the null; lane T ran it with a common direction bank: **196% accounted**, split-half transfer 19%, and the **transferable rule lands +0.472 Å WORSE than production** (1.88× MDE, 5/5 folds). Its own diagnosis: the matched random-subset null asked whether a structured class beats an unstructured one, not whether **the winning direction is the same direction twice** |

### A.3 Lane L's

| claim | how it died |
|---|---|
| **EDM projection / triangle repair as a route** | Grepped the codebase after deriving it: S19 had already measured it and left a **named prohibition**. Refuted in the *unexpected* direction — an incoherent magnitude-matched field is worse on realisability and lands **1.24 Å better** |
| **The deployable half of its own reference-state proposal** | Its own algebra refutes it: the score decomposes additively, so the size-matched field **is** the deflated field, and deflation reallocates rank without creating any |
| **"The AMBER-relax benefit concentrates on divergent pools", 2:1** | Lowered to **roughly even** by lane L itself, on lane F's dispersion null (−0.128, CI includes zero) — adjacent evidence pointing the other way |
| **"Two independent arrivals" on confidently-wrong** | Downgraded to one literature result plus one project result on a *different object*, after lane F measured error×confidence (−0.729) as **worse** than error alone (−0.799) |

### A.4 Lane X's

| claim | how it died |
|---|---|
| **Its own S30-L10 admission condition** | Amended to be **stricter**: Δ(set mean) is two channels, not one, so `ADMIT iff Δ(bias B) + 0.298·Δ(set best) < 0` |
| **Its counting falsifier** | Volunteered as having **essentially tied its bar — 24.80% against 25%** — and therefore deciding nothing; the verdict rests on the transfer arm alone |
| **A −0.4744 correlation it could have quoted** | Flagged as **driven by one arm**; drop it and the value is +0.0851. "The strong negative is an artefact and must not be quoted" |

### A.5 Lane F's, lane D's, lane R's, lane Q's

| claim | how it died |
|---|---|
| **Lane F's F1c FAIL18 row** (+1.7674 Å) | The adversary: FAIL18 is *defined* as the targets whose top-75 retained **zero** in-band members, so **49% is forced arithmetic** and the rest is truncation-inflated. Lane F *had* checked circularity against production; the predicate is the **filter's own recall**. **"A matched control in the right space does not rescue a stratum defined by the outcome"** |
| **Lane F's own F2a and F2c** | Both refuted by its own measurements — no knee in filter width, and its widening result flagged by itself as circular |
| **Three defects in lane D's own combination arms** | An untuned ridge returning a negative cosine; a least-squares objective violating its own consistency floor (caught because ρ fell below the best single field, impossible for a true maximum); selection leakage worth 0.019 in ρ |
| **Lane R's 0.500/0.500/0.500 cell** | A **null-input artefact** — three identical values across three different questions is what caught it |
| **Lane Q's framing of my point 4** | "The framing survives; the operator does not" — identity is not only the better question but the **cheaper** one, so a sparse weighted readout is the expensive way to ask it |

### A.6 The three checklist entries this sprint earned

1. **A matched control in the right space does not rescue a stratum defined by the outcome.**
2. **When a statistic is a nonlinear transform, its mean, median and win-rate have three different
   nulls** — never read one against another's baseline.
3. **"The winning direction is the same direction twice"** is a failure mode distinct from the
   order-statistic one, and needs a common direction bank to detect.
