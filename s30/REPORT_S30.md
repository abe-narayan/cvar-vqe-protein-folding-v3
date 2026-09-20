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

The sprint's most distinctive feature was lanes destroying their own results. Recorded in full
because the charter asks for it and because it is why the surviving results are trustworthy.

[PENDING — the full table.]
