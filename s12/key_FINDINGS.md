# RETRIEVAL-KEY agent — can a sequence-predicted LOCAL BACKBONE CONFORMATION profile replace BLOSUM62?

Sprint 12. All numbers on **tuning126**. benchmark60 never touched; dev24 never touched.
Code: `s12/key_lib.py`, `key_corpus.py`, `key_oracle.py`, `key_pred.py`, `key_emit.py`,
`key_report.py`. Results: `s12/results/key_*.json`.

## 0. Instrument validation (run first)

`python -m s12.instrument` →
`shipped 3.4540004952 / pool_best 1.7108244199 / top75_best 2.3061526409 /
synthesis_fit 3.2040761604 / n_zero_recall 18`, zero-recall set == `I.FAIL18`. **Exact.**

My own `key_lib.emit()` (recompute shipped score over a pool → top-75 → coordinate average →
project) reproduces the production record to ≤ 4e-4 Å on spot checks
(1A13 2.6402 vs 2.6401, 1A1P 3.5654 vs 3.5658, 1CB3 2.9344 vs 2.9344, 1CEK 0.5252 vs 0.5251),
so the arms below are on the same code path as the shipped system.

**Metric definitions used throughout.** `pool best/mean` are over the K=500 pool. `band` is
the number of pool members within **universe-best + 1.5 Å** — deliberately referenced to the
*universe*, not to each arm's own pool best, so the number is comparable across arms (an
arm-relative band would move its own goalposts). `fit` is the emitted CA-RMSD at λ=0 (the
3.204 arm the record calls "synthesis") and `proj` at λ=0.3; they track each other to
0.01 Å throughout, so only `fit` is discussed. `argmin` is the 3.454 legacy arm.

---

## E1 — Reproducing the ORACLE torsion-bin key, and a correction to what it measures

**Hypothesis.** Retrieving K=500 by torsion-bin match instead of BLOSUM62 sum improves the pool.

**The correction.** The literature agent defined the oracle key as *the ABEGO string of the
universe's true-best window*. That key is **circular**: the window whose bin string the key
*is* necessarily scores a perfect match, so retrieval recovers the universe best essentially
by construction (`abego_pool_best == univ_best` on 17/18 FAIL18 rows of
`lit_abego_retrieval.json`, and on 126/126 of mine). No sequence predictor can ever produce
that string — a predictor predicts the **native's** local conformation. I therefore separate:

| key | definition | reachable by a predictor? |
|---|---|---|
| `bestwin` | bin string of `argmin(rr)` in the universe (what lit measured) | **no** — circular |
| `native`  | bin string of the **native structure's own** (φ,ψ) | yes, in the limit |

`s12/results/key_oracle.json`, all 126 targets (lit used 46), K=500, same stable argsort.

| key (ORACLE) | ALL best | ALL mean | ALL band | F18 best | F18 mean | F18 band | O108 best | O108 mean | O108 band |
|---|---|---|---|---|---|---|---|---|---|
| BLOSUM (incumbent) | 1.711 | 4.453 | 84.8 | 2.284 | 5.964 | 16.0 | 1.615 | 4.201 | 96.3 |
| bestwin abego4 (circular) | 1.320 | 3.640 | 202.0 | **1.640** | 5.111 | 91.6 | 1.266 | 3.395 | 220.4 |
| **native abego4** | **1.450** | **3.646** | **202.0** | **1.768** | **5.096** | **96.2** | **1.397** | 3.405 | 219.6 |
| native abego5 | 1.511 | 3.651 | 191.4 | 1.877 | 5.123 | 73.9 | 1.450 | 3.405 | 210.9 |
| native km8 | 1.522 | 3.841 | 164.3 | 1.920 | 5.324 | 67.7 | 1.455 | 3.594 | 180.4 |
| native km16 | 1.509 | 3.990 | 133.3 | 1.803 | 5.366 | 64.4 | 1.460 | 3.761 | 144.8 |
| native abego4 + BLOSUM (z-sum) | 1.545 | 3.888 | 153.7 | 2.003 | 5.614 | 50.3 | 1.468 | 3.601 | 170.9 |
| universe best (bound) | 1.313 | — | — | 1.640 | — | — | 1.259 | — | — |

(band = windows within pool_best+1.5 Å, out of 500. FAIL18 lit numbers on 18 targets:
mine 1.640 for `bestwin` matches lit's 1.640 exactly; lit's 2.284 BLOSUM matches. **Reproduced.**)

**Findings.**
1. The lit result survives de-circularisation, with a haircut: the honest oracle ceiling on
   FAIL18 pool-best is **1.768**, not 1.640 (BLOSUM 2.284). Band recall 16 → **96 of 500** —
   a 6× increase, and this is the quantity the FAIL18 defect is defined by.
2. **Alphabet resolution: coarser is better.** abego4 dominates abego5, km8 and km16 on every
   pool statistic for the *native* key. Finer alphabets are simultaneously harder to predict
   and worse as a key. The design variable resolves in favour of the easiest-to-predict
   option. All subsequent work uses **abego4**.
3. Native-vs-bestwin bin agreement is only **0.738** overall / **0.645** on FAIL18 (abego4).
   So the pre-registered "≥ 0.70 per-residue accuracy" gate cannot be read against the
   bestwin string: the *native itself* only agrees with it at 0.645 on FAIL18 and still
   delivers 1.768. The gate must be read against the **native** bins, which is what §E2 does.
4. z-summing BLOSUM into the oracle key **hurts** it (F18 best 1.768 → 2.003, band 96 → 50).
   BLOSUM actively dilutes conformational retrieval.
5. Pool **mean** improves as much as pool best (4.453 → 3.646 overall; 5.964 → 5.096 on
   FAIL18). This matters: S7-6 says emitted RMSD tracks the pool mean. Both move the right
   way, so the mechanism is not the pool-best/pool-mean trade-off that kills most pool
   interventions. §E3 tests whether that reaches the emitted structure.

---

## E1b — Does the ORACLE key reach the EMITTED structure? (the ceiling worth chasing)

Full production chain, matched K=500 → shipped distogram score → top-75 → coordinate average
→ project. All 126 targets. `s12/results/key_emit_oracle.json`, `key_report_oracle.json`.
The `blosum` arm here is my recomputation, not the shipped record, and it lands on
**3.205 / 1.711 / 2.306 / 3.454** — i.e. it *is* the shipped system.

| metric | BLOSUM | ORACLE-native abego4 | Δ |
|---|---|---|---|
| **emitted fit (λ=0)** all | **3.205** | **2.912** | **−0.293** |
| emitted proj (λ=0.3) all | 3.213 | 2.918 | −0.295 |
| emitted fit FAIL18 | 6.034 | 5.363 | −0.671 |
| emitted fit other108 | 2.734 | 2.503 | −0.231 |
| argmin arm all | 3.454 | 3.271 | −0.183 |
| pool best / mean / band all | 1.711 / 4.453 / 84.8 | 1.450 / 3.646 / 202.0 | — |
| top-75 best all | 2.306 | 2.012 | −0.294 |

Paired (candidate − incumbent, emitted fit):

| group | n | Δ | CI95 | W/L | median | drop-top10 | drop-top20 | per fold |
|---|---|---|---|---|---|---|---|---|
| all | 126 | **−0.293** | [−0.439, −0.148] | **83/43** | −0.099 | −0.129 | **−0.040** | −0.30/−0.59/−0.10/−0.19/−0.30 |
| FAIL18 | 18 | −0.670 | [−1.423, +0.074] | 15/3 | −0.498 | +0.566 | n/a | f1 −1.34, f2 −0.18, f3 **+0.87**, f4 −1.19 |
| other108 | 108 | −0.230 | [−0.346, −0.121] | 68/40 | −0.082 | −0.103 | −0.020 | all negative |

**The ceiling is real and it reaches the emitted structure**, and the sign is broad
(83W/43L, every fold negative, median −0.099). But it is **heavy-tailed**: drop-top-20
leaves only −0.040 Å. So the honest reading is *"a real, broadly-positive but
tail-dominated −0.29 Å"*, not *"a clean −0.29 Å across the board"*. The FAIL18 CI crosses
zero at n=18 and fold 3's FAIL18 targets get **worse**.

Crucially the gain is on the SYNTHESIS arm (what the system emits) and is *larger* there
than on the argmin arm (−0.293 vs −0.183), and both pool best and pool mean move the right
way — so it is not the S7-6 best-vs-mean trap. **Everything now depends on whether a
sequence predictor can reach a usable fraction of it.**

---

## E1d — The ALPHABET-RESOLUTION axis, and settling the lit-vs-forensics disagreement

The coordinator's reconciling hypothesis was that value lives at FINE resolution (3-state SS
nearly useless, finer binning sharp). **That hypothesis is refuted — the dependence runs the
other way.** All ORACLE, all 126 targets, K=500 (`key_oracle.json`).

| ORACLE key (native bins) | states | ALL best | ALL **mean** | ALL band | F18 best | F18 **mean** | F18 band |
|---|---|---|---|---|---|---|---|
| BLOSUM (incumbent) | — | 1.711 | 4.453 | 84.8 | 2.284 | 5.964 | 16.0 |
| **ss3** (H/E/C from torsions) | 3 | 1.463 | 3.700 | 180.2 | **1.690** | **5.028** | 96.2 |
| **abego4** | 4 | **1.450** | **3.646** | **202.0** | 1.768 | 5.096 | 96.2 |
| abego5 | 5 | 1.511 | 3.651 | 191.4 | 1.877 | 5.123 | 73.9 |
| km8 (Ramachandran k-means) | 8 | 1.522 | 3.841 | 164.3 | 1.920 | 5.324 | 67.7 |
| km12 | 12 | 1.501 | 3.909 | 148.9 | 1.888 | 5.430 | 59.6 |
| km16 | 16 | 1.509 | 3.990 | 133.3 | 1.803 | 5.366 | 64.4 |

Band recall and pool mean degrade **monotonically** with resolution from 4 states upward, and
3-state SS is essentially tied with (on FAIL18, slightly better than) 4-state ABEGO. The
mechanism is obvious in hindsight: a finer alphabet makes an exact-match key so specific that
almost no window matches well, so the pool fills with near-arbitrary partial matches. Coarse
is right, and coarse is also the easiest thing to predict — the design variable resolves
favourably on both axes. **Everything from here uses abego4 (and ss3 as its control).**

### The disagreement, settled

There is no real disagreement between the literature agent and the forensics agent — they
measured different mechanisms and the coordinator's summary conflated a filter with a key:

| source | arm | what it did | pool best | emitted Δ |
|---|---|---|---|---|
| S6-6 | ORACLE SS **FILTER** before the argmin | keeps windows already in the BLOSUM pool | — | −0.21 |
| forensics `class_only` | PDB structural class as key | class flags far too coarse | 1.635 | −0.027 (CI spans 0) |
| forensics **`o_shape_only`** | ORACLE native-SS agreement + size as PRIMARY key | **a key** | **1.497** | **−0.348** [−0.509,−0.198], F18 −1.116, 91W/35L |
| **this agent** `oracle_native` abego4 | ORACLE native torsion bins as PRIMARY key | **a key** | **1.450** | **−0.293** [−0.439,−0.148], F18 −0.670, 83W/43L |
| lit `abego` | ORACLE **best-window** bins as key | circular (see E1) | 1.320 | not measured |

The forensics agent's own KEY arm and mine agree to within 0.05 Å on pool best, 1.5 windows on
band recall, and 0.055 Å on the emitted effect, from two independent implementations with
different alphabets. The "SS as a key is worth almost nothing" reading comes from the FILTER
arm and the `class_only` arm, not from `o_shape_only`. **Lit and forensics both measured the
same real thing; only lit's magnitude was inflated by circularity.**

---

## E1c — Corruption ladder on the NATIVE key, and the CAPACITY NULL

`s12/key_noise.py` → `key_noise.json` / `key_noise_agg.json`. All 126 targets, pool level,
5 repeats per point. `acc` = per-residue agreement of the corrupted key with the native bins.
Three error models: `uniform` (random bin), `biased` (collapse to the string's own modal
bin — how a majority-class predictor fails), `marginal` (draw from the library's pooled
Ramachandran marginal — this is **the coordinator's capacity null** and at p=1.0 it is a
content-free key).

**all 126** (BLOSUM: best 1.711 / mean 4.453 / band 84.8)

| model, p | acc | pool best | pool **mean** | band |
|---|---|---|---|---|
| oracle p=0 | 1.000 | 1.450 | 3.646 | 202.0 |
| marginal 0.2 | 0.890 | 1.518 | 3.820 | 170.7 |
| marginal 0.3 | 0.832 | 1.576 | 3.925 | 147.7 |
| marginal 0.4 | 0.783 | 1.630 | 4.017 | 131.9 |
| marginal 0.5 | 0.730 | 1.682 | 4.114 | 120.0 |
| marginal 0.7 | 0.611 | 1.898 | 4.421 | 77.4 |
| **marginal 1.0 (CAPACITY NULL)** | **0.445** | **2.268** | **4.752** | **43.5** |
| biased 1.0 (all-one-letter key) | 0.718 | 2.029 | 4.098 | 175.9 |
| uniform 1.0 | 0.249 | 2.183 | 4.989 | 28.7 |

**FAIL18** (BLOSUM: 2.284 / 5.964 / 16.0): oracle 1.768/5.096/96.2; marginal 0.5 →
2.136/5.496/54.2; **marginal 1.0 → 2.828/5.907/17.2**; biased 1.0 → 2.640/5.731/66.7.

**The capacity null is dead at pool level.** A content-free bin key drawn i.i.d. from the
library's own Ramachandran marginal is WORSE than BLOSUM on every statistic and on every
subgroup (all: best +0.56, mean +0.30, band −41; FAIL18: best +0.54, mean −0.06, band +1.2).
The oracle key's advantage is therefore *conformational information about this target*, not
selection capacity — the failure mode that killed the assembly agent's headline. It is
re-tested end-to-end (through emission) in §E4.

**Where the gain dies (pre-registered gate, read off the `marginal` model, which is the
right yardstick — `biased` inflates `acc` because agreement with a majority-collapsed
string is not skill):**

| criterion | accuracy at which the key ties BLOSUM |
|---|---|
| pool **mean**, all 126 | **≈ 0.60** |
| pool best, all 126 | ≈ 0.72 |
| pool mean, FAIL18 | ≈ 0.55 |
| pool best, FAIL18 | ≈ 0.70 |

Chance accuracy under the marginal is 0.445, and majority-class ("all A") is 0.53. So the
predictor must beat chance by ≈ 0.15 absolute to reach the pool-mean gate.

---

## PRE-REGISTERED PREDICTION (written before any predicted-key number was computed)

The coordinator's framing (S7-6: emitted RMSD tracks pool MEAN, moves against pool BEST;
assembly's pool-best gain of 0.536 Å emitted +0.218 Å WORSE) gives a sharp prediction, and
E1b already contains its confirmation for the oracle arm:

1. The oracle key moves pool **mean** by −0.807 Å (4.453 → 3.646) *and* pool best by
   −0.261 Å. Because the MEAN moves, S7-6 predicts the emitted number moves — and it did
   (−0.293 Å). This is the mechanism that distinguishes a retrieval key from assembly: a key
   re-weights the *whole* pool, assembly adds one good member. **I therefore predict the
   emitted gain of any key arm to track its pool-MEAN gain, not its pool-BEST gain**, and I
   will regress it explicitly (§E5).
2. I predict a per-residue-accuracy-to-emitted-gain transfer of roughly
   `Δ_emitted ≈ Δ_oracle_emitted × (Δmean_arm / Δmean_oracle)`. At a plausible predictor
   accuracy of 0.65–0.75 the ladder gives Δmean ≈ −0.35 to −0.45 of the oracle's −0.807, so
   I predict **Δ_emitted ≈ −0.10 to −0.16 Å**, i.e. real but roughly half the oracle, and
   NOT enough on its own to reach the sprint's <2.0 Å objective.
3. I predict the **shuffled-label** and **marginal-random** keys land at or worse than
   BLOSUM on emitted RMSD (pool-level they already do), and the **composition-only** key
   lands between BLOSUM and the predicted key but much closer to BLOSUM — because BLOSUM
   already carries composition.
4. Risk I flag in advance: `biased 1.0` (an all-one-letter key, zero positional information,
   0.718 nominal "accuracy") gives band 175.9 and pool mean 4.098 — *better than BLOSUM's
   mean*. A predictor that collapses to the majority state would reproduce much of the
   pool-mean gain with no positional content. **The composition null and the confusion
   matrix are what separate these two, and they are load-bearing, not decoration.**

---

## E1e — The CAPACITY NULL end to end, and a result that splits the story in two

`key_emit_nulls.json` / `key_report_nulls.json`. Full chain, all 126, paired vs BLOSUM,
emitted `fit`.

| arm | group | emitted | Δ vs BLOSUM | CI95 | W/L | med | d10 | d20 |
|---|---|---|---|---|---|---|---|---|
| ORACLE native abego4 | all | 2.912 | −0.293 | [−0.439,−0.148] | 83/43 | −0.099 | −0.129 | −0.040 |
| ORACLE native **ss3** | all | 2.947 | −0.258 | [−0.400,−0.120] | 87/39 | −0.093 | −0.092 | −0.014 |
| **rand_marginal (CAPACITY NULL)** | all | 3.730 | **+0.525** | [+0.347,+0.715] | 40/86 | +0.343 | +0.674 | +0.773 |
| ORACLE native abego4 | FAIL18 | 5.363 | −0.670 | [−1.423,+0.074] | 15/3 | −0.498 | — | — |
| ORACLE native **ss3** | FAIL18 | 5.307 | **−0.726** | **[−1.438,−0.044]** | 13/5 | −0.387 | — | — |
| **rand_marginal** | **FAIL18** | 5.680 | **−0.353** | [−0.763,+0.011] | **12/6** | −0.159 | — | — |
| ORACLE native abego4 | other108 | 2.503 | −0.230 | [−0.346,−0.121] | 68/40 | −0.082 | −0.103 | −0.020 |
| ORACLE native ss3 | other108 | 2.553 | −0.180 | [−0.290,−0.070] | 74/34 | −0.079 | −0.056 | +0.008 |
| **rand_marginal** | other108 | 3.405 | **+0.672** | [+0.479,+0.872] | 28/80 | +0.428 | +0.812 | +0.920 |

**Two findings, and the second is the more important one.**

1. **Globally the capacity null is comprehensively dead.** A content-free bin key drawn
   i.i.d. from the library's Ramachandran marginal makes the system **0.525 Å worse**,
   40W/86L. So the oracle key's global −0.293 Å is conformational information about this
   target, not selection capacity. It survives the control that killed the assembly agent's
   result. Also: 3-state SS and 4-state ABEGO are **statistically indistinguishable**
   end-to-end (−0.258 vs −0.293, overlapping CIs), and on FAIL18 ss3 is the only arm whose
   CI excludes zero. The resolution axis is settled: **coarse wins.**

2. **On FAIL18 specifically, HALF the oracle gain is capacity, not information.** A random
   bin key buys −0.353 Å on FAIL18 (12W/6L) against the oracle's −0.670. The conformational
   *information* is worth only ≈ 0.32 Å beyond a content-free key there. On the other-108 the
   separation is enormous (oracle −0.230 vs random +0.672, a 0.90 Å gap), so on the easy set
   the effect is entirely informational. Mechanistically this says: **on FAIL18 the BLOSUM
   pool is not merely uninformative, it is actively harmful — almost anything that breaks
   BLOSUM's ordering helps.** That is a new and independently useful fact about the FAIL18
   defect, and it means FAIL18 gains from *any* new key must be discounted by the random-key
   baseline before being called conformational. Every FAIL18 number below is reported against
   `rand_marginal` as well as against BLOSUM.

---

## E1f — Adjudication vs the forensics agent (coordinator's items a/b/c)

**(b) Concentration test applied to every ORACLE arm I have.** This is the test the
coordinator says matters, and it is the one that hurts.

| ORACLE arm | group | Δ emitted | CI95 | W/L | median | **drop-10** | **drop-20** |
|---|---|---|---|---|---|---|---|
| native abego4 | all | −0.293 | [−0.439,−0.148] | 83/43 | −0.099 | **−0.129** | **−0.040** |
| native ss3 | all | −0.258 | [−0.400,−0.120] | 87/39 | −0.093 | **−0.092** | **−0.014** |
| native abego4 | other108 | −0.230 | [−0.346,−0.121] | 68/40 | −0.082 | −0.103 | −0.020 |
| native ss3 | other108 | −0.180 | [−0.290,−0.070] | 74/34 | −0.079 | −0.056 | +0.008 |
| *(forensics `o_shape`)* | all | −0.179 | — | — | — | *+0.004* | — |
| *(forensics `o_shape_only`)* | all | −0.348 | [−0.509,−0.198] | 91/35 | — | *−0.149* | — |

My oracle arms **do** survive drop-10 (−0.129 / −0.092), unlike the forensics `o_shape` arm
(+0.004) — but they agree closely with the forensics agent's own *pure-key* arm
`o_shape_only` (−0.149). The `+0.004` figure comes from the arm that keeps BLOSUM as the
base and adds a shape bonus; the pure-key arms from both agents agree. **They do NOT survive
drop-20**: −0.040 (abego4) and −0.014 (ss3), i.e. at or below the sprint's 0.03 Å ignore
threshold. So the correct statement about the oracle ceiling is:

> A perfect local-conformation key is worth −0.29 Å on the mean, with the right sign on
> 83/126 targets and a −0.099 Å median, but roughly **two thirds of the mean effect lives in
> 20 targets**. It is a real, broad, small-per-target effect with a heavy tail — not a
> −0.29 Å improvement you can expect on an arbitrary target.

**(c) The capacity null: done, and it splits.** See §E1e. Globally dead (+0.525 Å, 40W/86L,
so key-space gains are NOT generic selection capacity). On **FAIL18** a content-free key
buys −0.353 Å of the oracle's −0.670 — so on the FAIL18 the coordinator's suspicion is
**half right**, and every FAIL18 key gain must be discounted by that baseline.

**On item 5 ("FAIL18 is not query-limited").** Not in conflict with anything I measured, and
I think the forensics framing is right. My FAIL18 band counts are 16/500 (BLOSUM) vs 96/500
(oracle key) — a 6× change — but the *universe-level* scarcity is the binding constraint:
FAIL18 universes simply contain far fewer near-native windows, so equal in-pool fractions
still leave FAIL18 with a sixth of the band members. A key can only redistribute what the
library contains. That is consistent with the forensics decomposition (filter 1.40 Å,
retrieval 0.54 Å) and I do not dispute it: **retrieval is the smaller of the two levers.**

**(a) is the decisive item and needs my predictor's peptide-level accuracy — §E2b below.**

---

## E2a — Leakage audit of the predictor (written before its numbers)

The training corpus for fold *f* is `key_corpus.corpus(f)`, obtained by **stitching the
sliding windows of the universe file of a fold-*f* target back into their parent chains**
(consecutive windows overlap by n−1, so parents are recoverable exactly; `key_corpus.stitch`
also checks torsion continuity, not just sequence). Per fold that is ≈ 6.6 k parents /
92 k residues, ≈ 10.7 k of them peptide residues and ≈ 81 k protein-fragment residues.

| channel | status |
|---|---|
| target's own fold peptides | **excluded by construction** — the universe of a fold-*f* target contains only out-of-fold peptides, and the 126 targets are peptides |
| target's own native torsions | never in training; used ONLY to score the predictor and to build arms explicitly labelled ORACLE |
| the retrieved pool / the shipped top-75 | never touched by the predictor (this is the failure mode lit's "modal ABEGO of the top-75" arm hit: 3.890 Å) |
| labels | torsion bins of **library** parents, from their own deposited (φ,ψ) |
| the k-means discretisations | fitted unsupervised on library torsions from 12 random targets' universes. A cross-fold statistic — but abego4 and ss3, the only alphabets used for any reported arm, are **fixed analytic rules with no fitted parameter**, so the reported results have no exposure here |
| `MARGINAL` (capacity null) | a library-wide bin frequency, not a target property |
| residual exposure I do not claim to have removed | a target peptide's sequence could in principle occur inside one of the 13,751 `prots/` proteins the fragments come from. That is exactly the exposure the shipped distogram already has, and it is identical across every arm including BLOSUM, so it cannot create a *differential* effect |

---

## E2b — THE PRE-REGISTERED GATE, and it is FAILED on the FAIL18 and CLEARED on the other-108

`key_pred.py` → `key_pred_acc.json`, `key_pred_post.json`. Sliding ±7-residue context of
one-hot + Chou-Fasman/KD/charge/volume propensities → 2-layer MLP (256/128, dropout 0.2,
Adam, early stop on a parent-grouped 10 % split) → 4-way abego4 posterior. Five models, one
per pinned fold, each trained ONLY on that fold's legal corpus.

**Accuracy on held-out LIBRARY parents (the easy number):** 0.804–0.820 (`full`),
0.829–0.839 (`frag`), 0.702–0.773 (`pep`), 0.618–0.637 (`comp`), 0.521–0.531 (`shuf`).

**Accuracy on the 126 TARGETS' native bins (the number that matters).** Majority class on
target residues is **A at 0.562** (bin frequency A .562 / B .366 / G .056 / E .015):

| variant | training data | **acc all** | **acc FAIL18** | acc other108 | top-2 | mean entropy |
|---|---|---|---|---|---|---|
| **majority-class baseline** | — | **0.562** | — | — | — | — |
| `full` | peptides + fragments | 0.618 | 0.425 | 0.650 | 0.917 | 0.461 |
| **`pep`** | **peptides only** | **0.690** | **0.517** | **0.718** | 0.924 | 0.608 |
| `frag` | fragments only | 0.609 | 0.461 | 0.634 | 0.906 | 0.391 |
| `comp` (NULL) | composition only | 0.600 | **0.542** | 0.610 | 0.922 | 0.778 |
| `shuf` (NULL) | permuted labels | 0.554 | 0.375 | 0.584 | 0.927 | 0.924 |

Per-class recall, `pep`: A 0.828, B 0.624, **G 0.011, E 0.000**. `shuf` per-class recall is
exactly [1, 0, 0, 0] — it collapses to "always A", which is what a dead null should do
(0.554 ≈ the 0.562 majority).

**Verdict against the pre-registered gate (§E1c: 0.60 for the pool-mean tie, 0.72 for the
pool-best tie; FAIL18 gates 0.55 / 0.70):**

| | acc | pool-mean gate | pool-best gate | verdict |
|---|---|---|---|---|
| other-108 | **0.718** | 0.60 ✓ | 0.72 ≈ tie | **CLEARED** (mean), marginal (best) |
| **FAIL18** | **0.517** | 0.55 ✗ | 0.70 ✗ | **FAILED** |
| all 126 | 0.690 | 0.60 ✓ | 0.72 ✗ | cleared on mean only |

Stated before looking at any emitted RMSD, as pre-registered.

### Three findings inside this table that are worth more than the headline

1. **S7-2 is confirmed independently, and it is large.** Training on peptides ONLY
   (10.5 k residues) beats training on peptides + 81 k fragment residues by **+0.072**
   absolute accuracy (0.690 vs 0.618), and fragment-only is worst of the three (0.609).
   Eight times more data makes the predictor worse. The brief told me to measure this rather
   than assume it; measured, and the distribution shift is the dominant effect. **A residue's
   local conformation inside a folded protein is not what it is in an isolated peptide** —
   the same conclusion the forensics agent reached in a different model family.
2. **The composition null is alarmingly strong, and on FAIL18 it BEATS the real predictor**
   (0.542 vs 0.517). Most of what the sequence predictor knows is amino-acid propensity,
   which BLOSUM already carries. This was flagged in advance (§Prediction 4) as the most
   likely way the result dies, and it is exactly what the numbers show.
3. **ESM-2 adds nothing to this task** (`key_pred_esm_acc.json`). Adding the target's ESM-2
   PCA-32 per-residue embedding plus its context mean (64 extra features, same architecture,
   same LFO, peptide-only corpus): **0.688 / 0.508 / 0.718** against `pep`'s 0.690 / 0.517 /
   0.718. Identical to three decimal places on the other-108 and slightly *worse* on FAIL18.
   Library val-accuracy also barely moves (0.779 vs 0.773 on fold 0). So the ≈ 0.69 ceiling
   is **not** a feature-engineering failure — the strongest sequence representation this
   project has adds zero local-conformation signal over a ±7 one-hot window. (Note this does
   NOT contradict S7-11's −0.288 Å for ESM on *selection*; it says ESM's contribution there
   is not mediated by local backbone conformation.)
4. **G and E are unpredictable.** Recall 0.011 and 0.000 for `pep`. The predictor is a
   two-state A/B classifier wearing a four-state hat. Left-handed and extended-left
   conformations — precisely the states that distinguish an unusual peptide backbone from a
   generic one — carry no predicted signal at all.

### Adjudication of coordinator item (a)

The forensics agent's structural-alphabet predictor came in **below** its majority baseline
(0.375 vs 0.393). Mine comes in **above** it (0.690 vs 0.562, and 0.718 vs 0.562 on the
other-108). So item (a) is **not** a second falsification of the same kind — the two results
differ, and I believe the difference is real and explainable: (i) alphabet resolution — I
show above that coarse (3–4 state) alphabets are both more predictable and *better keys*,
and a finer structural alphabet is where their 0.393 baseline comes from; (ii) **training
corpus** — my `full` variant (peptides+fragments, the natural choice) scores 0.618, and only
dropping the 81 k fragment residues gets to 0.690. An agent that trained on the fragment-rich
corpus would land near or below baseline, which is consistent with what they report.

**But this does not rescue the hypothesis**, because the gate fails exactly where it needed
to pass: FAIL18 accuracy 0.517, below both the 0.55 pool-mean gate and the 0.542 that a
composition-only null achieves. The measured emitted numbers follow in §E3.

**Note on selecting `pep` over `full` (`key_pred_oof.json`).** Choosing the variant by its
accuracy on the 126 targets would be selection on oracle labels. I therefore did the choice
**leave-one-fold-out**: for each fold, pick the variant with the best accuracy on the other
four folds' targets. The choice is `pep` **unanimously on all five folds** (OOF accuracy
0.671–0.710 for `pep` vs 0.600–0.633 for `full` vs 0.593–0.616 for `frag`), and the assembled
OOF-selected predictor scores exactly 0.690 / 0.517. So `pep` is a deployable choice, not a
post-hoc one.

---

## E3a — Pool statistics of the DEPLOYABLE keys, and the null that kills it

`key_pool.py` → `key_pool.json` / `key_pool_agg.json`. All 126 targets, matched K=500.
`pred:X` = soft key (Σ log P), `predhard:X` = argmax key, `mix:full:w` = (1−w)·z(BLOSUM) + w·z(soft).

| key | ALL best | **ALL mean** | ALL band | F18 best | F18 mean | F18 band | O108 best | O108 mean | O108 band |
|---|---|---|---|---|---|---|---|---|---|
| BLOSUM | **1.711** | 4.453 | 84.8 | **2.284** | 5.964 | 16.0 | **1.615** | 4.201 | 96.3 |
| ORACLE native | 1.450 | **3.646** | **202.0** | 1.768 | **5.096** | **96.2** | 1.397 | **3.405** | **219.6** |
| `pred:full` | 1.983 | 4.221 | 132.4 | 3.331 | 6.336 | 6.6 | 1.758 | 3.869 | 153.4 |
| **`pred:pep`** | 2.075 | **4.082** | 147.1 | 3.778 | **5.841** | 15.5 | 1.792 | **3.789** | 169.0 |
| `pred:frag` | 1.997 | 4.263 | 130.2 | 3.170 | 6.241 | 17.7 | 1.801 | 3.934 | 148.9 |
| **`pred:shuf` (NULL)** | 2.969 | **4.099** | **145.9** | 4.634 | **5.800** | **25.7** | 2.692 | **3.816** | **165.9** |
| **`pred:comp` (NULL)** | 2.419 | 4.500 | **150.0** | 2.880 | 6.204 | **51.8** | 2.342 | 4.216 | 166.4 |
| `mix:full:0.25` | 1.703 | 4.283 | 104.9 | 2.390 | 6.004 | 14.9 | 1.588 | 3.996 | 119.9 |
| `mix:full:0.75` | 1.852 | 4.166 | 138.9 | 3.207 | 6.294 | 12.3 | 1.627 | 3.812 | 160.0 |

**Read the two bolded columns together. This is the result.**

- Every deployable predicted key makes pool **best** clearly WORSE than BLOSUM
  (1.98–2.08 vs 1.711; on FAIL18 3.2–3.8 vs 2.284). Only `mix:full:0.25`, which is 75 %
  BLOSUM, ties it.
- Every deployable predicted key makes pool **mean** BETTER than BLOSUM
  (`pred:pep` 4.082 vs 4.453, −0.371 Å) and roughly doubles band recall (147 vs 85).
- **And the shuffled-label null does exactly the same thing.** `pred:shuf` gives pool mean
  **4.099** against `pred:pep`'s 4.082 — a difference of **0.017 Å** — and band recall
  145.9 against 147.1, a difference of **1.2 windows out of 500**. On FAIL18 the null is
  *better* than the real predictor on both (5.800 vs 5.841; 25.7 vs 15.5 band).
  The composition null has the highest band recall of any deployable arm (150.0 overall,
  **51.8** on FAIL18 — more than triple BLOSUM's 16.0).

So the pool-mean and band-recall gains of the deployable key are **entirely reproduced by
keys that contain no conformational information about the target**. This is precisely the
failure mode written down in advance (§Prediction 4): a key that says "prefer windows whose
backbone is uniformly α" concentrates the pool on regular, compact, low-variance geometry —
which lowers the pool mean and raises band recall for generic reasons — whether or not the
target is α. `pred:shuf` collapses to "always A" (per-class recall [1,0,0,0]) and that alone
buys the whole effect.

The only quantity that separates the real predictor from the nulls is **pool best**, where
`pred:pep` (2.075) beats `pred:shuf` (2.969) by 0.89 Å — i.e. the real predictor does carry
information, and the information lands exactly on the statistic that the record (S7-6) and
the assembly agent both say **does not reach the emitted answer**.

Prediction, restated before the emitted measurement: `pred:pep` and `pred:shuf` will emit
within ~0.05 Å of each other; whether that is above or below BLOSUM is the open question,
but the *difference between them* is the test, and I expect it to be null.

---

## E4 — Where between PREDICTED and ORACLE does the gain appear? (`key_interp.json`)

Replace a fraction *r* of the predictor's argmax bins with the ORACLE bin; *r*=0 is the
deployable `full` key, *r*=1 the oracle key. 5 repeats per rate, all 126 targets, pool level.

| r | acc (all) | best | **mean** | band | | acc F18 | best F18 | mean F18 | band F18 |
|---|---|---|---|---|---|---|---|---|---|
| BLOSUM | — | **1.711** | 4.453 | 84.8 | | — | **2.284** | 5.964 | 16.0 |
| 0.0 (deployable) | 0.618 | 2.002 | 4.210 | 130.8 | | 0.425 | 3.425 | 6.294 | 4.4 |
| 0.2 | 0.698 | 1.813 | 4.082 | 142.0 | | 0.553 | 2.822 | 6.065 | 16.1 |
| 0.4 | 0.775 | **1.708** | 3.977 | 156.0 | | 0.649 | 2.470 | 5.826 | 25.6 |
| 0.6 | 0.848 | 1.585 | 3.856 | 172.0 | | 0.794 | **2.039** | 5.577 | 57.3 |
| 0.8 | 0.920 | 1.507 | 3.757 | 185.8 | | 0.871 | 1.816 | 5.405 | 66.2 |
| 1.0 (ORACLE) | 1.000 | 1.450 | 3.646 | 202.0 | | 1.000 | 1.768 | 5.096 | 96.2 |

The ladder is smooth and monotone, and it locates the death of the effect precisely:

- **Pool best ties BLOSUM at per-residue accuracy ≈ 0.775 overall and ≈ 0.70 on FAIL18.**
  (This independently reproduces the gate I pre-registered from the corruption ladder,
  0.72 / 0.70, by a completely different construction — corrupting an oracle vs. repairing a
  predictor. The two agree to 0.05.)
- The best predictor I can build reaches **0.690 overall and 0.517 on FAIL18**. The shortfall
  is **≈ 0.09 accuracy overall and ≈ 0.19 on FAIL18** — not a rounding error, and on FAIL18
  it is a gap of the same size as the entire distance from chance to the gate.
- Pool *mean* crosses BLOSUM immediately (already better at r=0), but §E3a showed that
  crossing is reproduced exactly by the shuffled-label null, so it is not evidence.

**This is the answer to "where does the gain die".** It dies between 0.69 and 0.78
per-residue accuracy on the whole set, and between 0.52 and 0.70 on the FAIL18 — and the
achievable accuracy sits below both thresholds. The oracle key's advantage is real and is
*not* accessible at achievable prediction accuracy, and the shortfall is largest exactly on
the subgroup that would have to move for the sprint objective to be reached.

---

## E6 — Does the ORACLE key fix the FAIL18 subclasses? (fibril / lasso), against the null

Class flags from `forensics/fail_headers.json`. Emitted `fit`, Δ vs BLOSUM in brackets.

| group | n | BLOSUM | ORACLE abego4 | ORACLE ss3 | **rand_marginal (null)** | information = oracle − null |
|---|---|---|---|---|---|---|
| all | 126 | 3.205 | 2.912 (−0.293) | 2.947 (−0.258) | 3.730 (+0.525) | **−0.818** |
| FAIL18 | 18 | 6.034 | 5.363 (−0.670) | 5.307 (−0.726) | 5.680 (−0.353) | −0.317 |
| **F18 fibril / steric zipper** | 6 | 5.809 | 4.830 (**−0.979**) | 5.006 (−0.803) | 5.525 (−0.284) | **−0.695** |
| **F18 lasso** | 4 | 6.825 | 6.188 (−0.637) | 6.129 (−0.696) | 6.847 (+0.022) | **−0.659** |
| F18 fibril **or** lasso | 10 | 6.216 | 5.373 (−0.842) | 5.455 (−0.760) | 6.054 (−0.162) | **−0.680** |
| F18 neither | 8 | 5.806 | 5.351 (−0.455) | 5.122 (−0.683) | 5.213 (**−0.592**) | **+0.137** |
| other108 | 108 | 2.734 | 2.503 (−0.230) | 2.553 (−0.180) | 3.405 (+0.672) | −0.902 |

**This is the cleanest mechanistic result in the report.** The conformational *information*
(oracle minus content-free null) is worth **−0.68 Å on the 10 fibril/lasso targets** and
**+0.14 Å (i.e. nothing) on the other 8 FAIL18 targets**, where the entire apparent gain is
the content-free "escape BLOSUM" effect. That is exactly the pattern the mechanism predicts:
fibril/steric-zipper segments are uniformly extended and lasso peptides have a distinctive
local backbone, so their *bin string* is unusual and highly discriminative; the remaining
FAIL18 targets fail for reasons a local-conformation key cannot see (partner-bound
conformations, tertiary contacts).

So *if* an accurate predictor existed, it would help precisely the 10 targets the forensics
agent identified as the structural core of the failure. It does not exist (§E2b), and on
those 10 targets the predictor is at its **worst** (FAIL18 accuracy 0.517 vs 0.718 elsewhere).

---

## E3b — EMITTED RMSD of the deployable arms. **The hypothesis is falsified.**

`key_emit_pred.json` + `key_report_all.json` / `key_report_all.log`. All 126 targets, full
production chain, matched K=500, paired vs the BLOSUM incumbent (3.205 emitted fit).

### All 126 targets

| arm | emitted | **Δ** | CI95 | **W/L** | median | drop-10 | drop-20 | per fold |
|---|---|---|---|---|---|---|---|---|
| BLOSUM (incumbent) | 3.205 | — | — | — | — | — | — | — |
| **`pred:pep`** (the hypothesis) | 3.220 | **+0.015** | [−0.114,+0.137] | **56/70** | +0.018 | +0.147 | +0.214 | −.02/+.06/+.01/+.23/−.15 |
| `mix:pep:0.25` | 3.200 | −0.005 | [−0.047,+0.037] | 71/55 | −0.008 | +0.038 | +0.058 | ~0 |
| `mix:pep:0.75` | 3.166 | −0.039 | [−0.150,+0.063] | 63/63 | +0.002 | +0.082 | +0.129 | mixed signs |
| `pred:comp` (NULL) | 3.701 | +0.496 | [+0.255,+0.734] | 41/85 | +0.073 | — | — | — |
| `pred:shuf` (NULL) | 3.984 | +0.779 | [+0.515,+1.061] | 35/91 | +0.204 | — | — | — |
| `rand_marginal` (NULL) | 3.730 | +0.525 | [+0.347,+0.715] | 40/86 | +0.343 | — | — | — |
| *ORACLE native abego4* | *2.912* | *−0.293* | *[−0.439,−0.148]* | *83/43* | *−0.099* | *−0.129* | *−0.040* | *all −* |

**No deployable arm beats BLOSUM.** The predicted key alone is +0.015 Å (worse, 56W/70L,
CI straddling zero). The best mixed arm is −0.039 Å with a CI that includes zero, **W/L
exactly 63/63**, a median of **+0.002**, and a drop-10 of **+0.082** — i.e. it is a tail
artefact with no per-target reality. Every arm is inside the brief's 0.03 Å ignore threshold
or on the wrong side of it.

### FAIL18, where the arms *do* move — and where the nulls destroy the interpretation

| arm | FAIL18 emitted | Δ | CI95 | W/L | median | drop-10 |
|---|---|---|---|---|---|---|
| BLOSUM | 6.034 | — | — | — | — | — |
| **`pred:pep`** | 5.524 | **−0.510** | **[−0.957,−0.149]** | 11/7 | −0.187 | +0.087 |
| **`pred:comp` (COMPOSITION NULL)** | **5.388** | **−0.646** | **[−1.377,−0.013]** | **12/6** | −0.157 | +0.366 |
| `rand_marginal` (CAPACITY NULL) | 5.680 | −0.353 | [−0.763,+0.011] | 12/6 | −0.159 | +0.274 |
| `pred:shuf` (SHUFFLED-LABEL NULL) | 5.888 | −0.145 | [−0.855,+0.419] | 7/11 | +0.033 | +0.750 |
| `mix:pep:0.75` | 5.592 | −0.441 | [−0.925,−0.044] | 11/7 | −0.023 | +0.168 |
| *ORACLE native* | *5.363* | *−0.670* | *[−1.423,+0.074]* | *15/3* | *−0.498* | *+0.566* |

**The FAIL18 result is the one that would have been reported as a success, and the
composition null kills it outright.** A key built from amino-acid composition alone — no
positional information whatsoever, a single posterior repeated at every residue — emits
**−0.646 Å** on the FAIL18, **more than the real conformation predictor's −0.510 Å**, with a
CI that also excludes zero and a better W/L. On this subgroup the ordering is:

> composition null (−0.646) > **real predictor (−0.510)** > capacity null (−0.353) >
> shuffled-label null (−0.145)

The real predictor sits *between two content-free controls*. The information it adds over a
composition key is **negative**. Every FAIL18 gain here is "anything but BLOSUM" (see §E1e:
on FAIL18 the BLOSUM pool is actively harmful), plus a composition effect BLOSUM already
carries in a worse form. Drop-10 is positive for every one of these arms.

### Verdict on the pre-registered predictions

| # | prediction | outcome |
|---|---|---|
| 1 | emitted gain tracks pool MEAN, not pool BEST | **partly right, and it backfired.** Every arm that improved pool mean (`pred:pep` 4.082, `pred:shuf` 4.099) emitted *worse* than BLOSUM overall. The pool-mean→emitted link that S7-6 describes did not carry a *deployable* key; only the oracle key, which improved mean AND best, moved the emitted number. |
| 2 | Δ_emitted ≈ −0.10 to −0.16 at 0.65–0.75 accuracy | **wrong, and too optimistic.** Actual: **+0.015**. The linear accuracy→gain transfer I assumed does not hold; the deployable key's error structure (A/B confusion, G/E never predicted) is worse than the corruption models used to derive the gate. |
| 3 | shuffled and marginal nulls land at or worse than BLOSUM | **right** globally (+0.779, +0.525) — but I did not predict that on FAIL18 they would land *better* than BLOSUM and that the composition null would *beat the real predictor*. |
| 4 | the composition null is the most likely way the result dies | **right, and it is exactly how it died.** |

---

## Verdict

**The retrieval-key hypothesis (H1) is CLOSED as not deployable.** Written as the record
would state it:

- The ORACLE is real: a perfect 3–4-state local-conformation key is worth **−0.293 Å**
  emitted [−0.439,−0.148], 83W/43L, and −0.68 Å of genuine information on the 10
  fibril/lasso FAIL18 targets. It is not a capacity artefact (the content-free key is
  +0.525 Å). But it is **tail-dominated** (drop-20 −0.040) and it is only about a third of
  what a perfect *filter* is worth (forensics: 1.40 Å vs 0.54 Å) — retrieval was always the
  smaller lever.
- The predictor is **as good as this data allows and not good enough**: 0.690 per-residue
  vs a 0.562 majority baseline, above the forensics agent's below-baseline result, unmoved
  by ESM-2 (0.688), and **0.517 on the FAIL18** — below the 0.55 pre-registered gate and
  below what a composition-only model achieves there (0.542).
- The gain dies between **0.69 and 0.78** per-residue accuracy (§E4), and the achievable
  accuracy is below that on the whole set and far below it on the subgroup that matters.
- End to end the deployable key is **+0.015 Å [−0.114,+0.137], 56W/70L**, and its only
  significant subgroup effect is **exceeded by a composition null**.

The binding constraint is not the alphabet, the model, or the features. It is that
**10.5 k peptide residues is the entire supply of in-distribution training data**, and the
81 k protein-fragment residues that would fix the sample size make the predictor *worse*
(0.690 → 0.618). That is S7-2 in a second model family, and it is the same wall the
forensics agent hit. More capacity does not help; more *peptide* data would, and there is
none — the record already establishes that all 204 identity clusters of 9–16-mers are spent.

**dev24 should NOT be spent on this.** There is nothing here that a confirmation split
could confirm.

---

## What I would do with the next budget (ranked)

1. **The filter, not the key.** Three threads (forensics' 1.40 Å vs 0.54 Å decomposition,
   the assembly agent's pool-best null, and this report) now agree that the terminal
   filter/aggregation is where the remaining accuracy lives, and the record's C1 correction
   says oracle top-25 *averaging* reaches 1.644 Å. The open question the record itself flags
   (C2) is **learned aggregation over the full n×n candidate-vs-objective deviation map** —
   the one class of in-band ranker that has never been tried, because every previous ranker
   saw scalar poolings. That is where I would put the next agent.
2. **Two facts from this report that are directly reusable elsewhere.**
   (a) **On FAIL18, BLOSUM is worse than a random key** (−0.353 Å just for replacing it with
   noise, 12W/6L). Nobody has exploited that: it means a *deliberately diversified* pool on
   the flagged targets is nearly free. Combined with the forensics agent's confidence signal
   (LFO ridge ρ = +0.516 with the answer, which ranks targets even though it cannot name
   catastrophes), a **selective pool-diversification rule on the top-k flagged targets** is a
   cheap, honest experiment that needs no new information channel.
   (b) **Local conformation is the FAIL18 fibril/lasso signature** (−0.68 Å of real
   information on those 10 targets, §E6). Since it cannot be *predicted*, the only way to get
   it is to *measure* it — which is what the literature agent's BMRB chemical-shift lead
   (§lit H-shift) proposes. If any external per-residue observable is admissible, this report
   quantifies exactly what it would be worth and on which targets.
3. **Do not fund**: finer structural alphabets (monotonically worse as keys, §E1d), ESM
   features for local conformation (zero, §E2b), more protein-fragment training data (makes
   it worse, §E2b), fusing keys (forensics: all keys r ≥ 0.87 with BLOSUM).

## Files

| file | contents |
|---|---|
| `s12/key_lib.py` | alphabets (ss3/abego4/abego5/km8/km12/km16), native torsions, retrieval, `emit()` |
| `s12/key_corpus.py` | window→parent stitching, per-fold LFO corpora |
| `s12/key_pred.py` | the torsion-bin predictor + all variants and nulls |
| `s12/key_oracle.py` / `key_noise.py` / `key_interp.py` | oracle alphabet sweep / corruption ladder / predicted↔oracle interpolation |
| `s12/key_emit.py` / `key_pool.py` / `key_report.py` | full-chain arms / cheap pool arms / paired statistics |
| `s12/results/key_oracle.json` | alphabet-resolution sweep, all 126 |
| `s12/results/key_noise.json`, `key_noise_agg.json` | corruption ladder + capacity null (pool) |
| `s12/results/key_pred_acc.json`, `key_pred_post.json`, `key_pred_oof.json`, `key_pred_esm_acc.json` | predictor accuracy, posteriors, OOF variant selection, ESM variant |
| `s12/results/key_pool.json`, `key_pool_agg.json` | pool statistics for every deployable key |
| `s12/results/key_emit_oracle.json`, `key_emit_nulls.json`, `key_emit_pred.json` | emitted RMSD, all arms |
| `s12/results/key_report_all.json`, `key_report_all.log` | the final paired tables |
| `s12/results/key_interp.json` | where the gain dies |
