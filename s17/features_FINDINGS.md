# FEATURES workstream — Sprint 17 findings

Pre-registration: `s17/PREREG_features.md`, written before any number below (F7's section written
before F7 ran, after F1–F5 returned; the F6 model ladder written before either F6 run).
Modules: `s17/feat_lib.py` · `feat_esm.py` · `feat_pll.py` · `feat_window.py` · `feat_shortlist.py` ·
`feat_report.py`
Artefacts: `s17/results/feat_inband.{json,log}` · `feat_shortlist.{json,log}` ·
`s17/cache/feat_pll.npz` · `s17/cache/feat_win_*.npz`

All numbers: **n = 126** tuning targets, **K = 500**, **identical candidate sets** (`sel_lib.pack`),
target as the unit, tie-safe selection, paired target-level **and** fold-clustered intervals,
medians and W/L beside every mean. The sealed 60-target benchmark was not read, probed, or derived.

---

## 0. THE ANSWER, first

**No representation-derived feature has in-band ARGMIN skill.** Nothing in this lane beats a
constant ideal α-helix *and* random-in-band on selected RMSD with a fold-clustered interval
excluding zero, in either deployable band. **The lane-level falsifier fired. The last open feature
class is closed and the sprint's selection answer is final.**

Five things the lane returned that are worth more than that sentence:

1. **ESM-2's contact head supplies the first ESM-SPECIFIC in-band content ever measured in this
   programme — on the correlation axis only.** In-band Spearman **+0.116 [+0.047, +0.184]** in the
   deployable top-75 band (counted pairwise accuracy 0.541) — **the largest in-band Spearman any
   native-free signal has reached in THAT band**, against SELECT's best-of-29 at +0.086 and the
   shipped objective's +0.065 measured here. (Band-qualified deliberately: PHYSICS measured the
   distance score at +0.126 on a *different* band, so an unqualified "largest in the programme"
   would be comparing two sets.) Against
   its **own ESM-free twin** (identical functional, uniform pair weights instead of ESM contact
   probabilities) the increment is **+0.050 [+0.005, +0.097]** at top-75 and **+0.062 [+0.030,
   +0.099]** on `oracle1.5` — CIs excluding zero. **SUPPORTED.**
   **And none of it reaches the argmin, in any band** (best selected: −0.026 [−0.202, +0.144]).

2. **That is not a paradox; it is the number the sprint needed.** Fed through
   `selection_theory`'s band curve, ρ_S = +0.116 (copula ρ 0.121) delivers **3.393 Å at B = 25**
   against a band mean of 3.502 — **0.109 Å**. The best in-band correlation measured in the
   deployable band is worth **a tenth of an Angstrom** against a 2.249 A gap, and the curve is so
   flat here that the largest in-band figure anyone in this sprint has reported (+0.131) buys only
   0.014 A more. And at **B ≤ 25 the 2.5 Å target is
   unreachable at ρ = 1**: a *perfect* in-band ranker returns 2.609 Å, because that is the band's
   best member. **The top-25 shortlist, not the ranker inside it, is the binding constraint.**
   **EXACT** — arithmetic on the band's RMSD distribution, not a discovery.

3. **The ESM-versus-one-hot contrast — the exact contrast the recorded 0.288 Å claim rests on — is
   null or NEGATIVE in-band.** At matched architecture the ESM arm is +0.036 [−0.063, +0.139]
   (top-25) and **+0.116 [−0.004, +0.243], fold [+0.019, +0.174]** (top-75) *worse* than one-hot at
   predicting local environment, and −0.019 [−0.173, +0.136] (null) on secondary structure. The
   sequence-context pseudo-likelihood does not separate from a **composition-only** control
   (+0.040 [−0.083, +0.167] at top-25; −0.002 [−0.151, +0.144] at top-75).
   **The recorded 0.288 Å is not in-band content — it is a better global FILTER.** **ESTABLISHED.**

4. **The role the signal might have had instead — shortlist construction — is closed too, and the
   closure comes with an independent reproduction of a PHYSICS result.** An ESM-ranked shortlist is
   **significantly worse than a matched random one** at every size (+0.315 to +0.449 on the
   ceiling). And the shipped score's own top-B shortlist has an ORACLE best **significantly worse
   than matched random** at every size — **+0.209 / +0.226 / +0.170** at B = 25/75/150 here against
   PHYSICS's independently written **+0.195 / +0.242 / +0.179**. **ESTABLISHED.** See §9.

5. **A candidate finding that is NOT a representation feature and needs independent confirmation.**
   The ESM-free control built only to price F1 — `cf_topd_uniform`, *the mean realised CA distance
   over the n pairs the distogram predicts closest* — is the **best-performing feature in the whole
   battery**: top-75 ρ_S +0.109 [+0.041, +0.177], selected **3.415**, **−0.135 [−0.252, −0.022]**
   against random-in-band, −0.196 against the α-helix, −0.348 against the β-strand, 72W/54L, and it
   beats the shipped objective (3.454). It is **one arm out of 25–28 features × 3 bands with no
   multiplicity correction**, its top-25 interval touches zero (−0.087 [−0.181, +0.000]), and it was built as a
   control rather than as a hypothesis. **PLAUSIBLE, explicitly not ESTABLISHED** — flagged for
   independent confirmation, not claimed. See §8.

---

## 1. Which pre-registered rules fired

| experiment | rule | fired? |
|---|---|---|
| F1 ESM contact head | falsifier: in-band ρ CI contains zero **and** selected-vs-random CI contains zero, both bands | **PARTIALLY** — ρ significant in both native-free bands and significant against its ESM-free twin; **selected null in all three bands**. Rank skill without argmin skill. |
| F2 exogenous pair weights | §1's bar | **FIRED** — `risk_wcon − risk_wone` +0.057 [−0.070, +0.188] (top-25), +0.021 [−0.107, +0.151] (top-75), −0.001 [−0.061, +0.059] (`oracle1.5`). |
| F3 window-sequence PLL | falsifier: does not separate from the background-frequency control | **FIRED** — `pll_mask − pll_bg` +0.040 / −0.002 / +0.060; and against BLOSUM it is *worse* on `oracle1.5` (+0.091 [+0.001, +0.178]). |
| F4 embedding → local env | falsifier: ESM arm indistinguishable from one-hot | **FIRED, in the wrong direction** — ESM worse, fold-clustered CI excluding zero at top-75 and `oracle1.5`. |
| F5 predicted vs realised SS | same falsifier | **FIRED** — ESM vs one-hot null in all three bands; the **zero-sequence** arm beats both. |
| F6 window embeddings | falsifier: does not separate from BLOSUM / the constant conformations | **FIRED** — residue-aligned cosine has counted ordering accuracy **0.503** (null 0.500); no arm beats the α-helix. The shipped 650M model, 2,369 window sequences embedded. |
| F7 shortlist construction | falsifier: no gain over the distance shortlist or over matched random | **FIRED, twice** — the ESM shortlist is worse than the incumbent's and **significantly worse than matched random** (+0.315 to +0.449); and the one positive arm (`mix`) is beaten by matched random at every size. |
| **LANE** | **no representation-derived feature beats a constant α-helix in-band with a fold-clustered CI excluding zero** | **FIRED** |

---

## 2. F0 — INVENTORY: what representation data exists, and what does not

Established before any modelling, because the lane's scope depends on it.

| artefact | contents | coverage | fold-honest? |
|---|---|---|---|
| `esm_small.npz` | per-residue ESM-2-650M final-layer reps, 1280-d | **360 sequences; all 126 targets present** | yes — sequence-only, a fixed pretrained artefact identical for every fold |
| `esm_pca.npz` | the shipped 1280 → 32 whitening projection | global | yes |
| `s12/cache/esm_con_targets.npz` | ESM-2 **contact-head** probabilities (n, n) | **all 126 targets** | the contact head is supervised on **other proteins'** structures, never the target — flagged at every use |
| `esm_cache.npz` | 1.5 GB, ~22k library sequences | keyed by **full parent sequence** | — |
| `s12/cache/esm_bank.npz` | 102 MB compact bank (pca32/pca128/contacts) | bank sequences | — |
| `s8/generate_univ/<pdb>.npz` `S` | **the source sequence of every candidate window**, int8 | full universe | yes |

> **PER-CANDIDATE ESM EMBEDDINGS DO NOT EXIST, and this was established before any modelling.**
> The universes hold 13k–27k length-n *windows* per target; the ESM cache is keyed by *parent*
> sequence with **no window→parent map**. Featurising the full universe is ≈1.6M short sequences —
> out of budget. Unique window sequences: **top-25 of the BLOSUM order 2,340; top-75 6,608;
> top-500 34,216.** The *band* is affordable; the universe is not; **no claim in this document
> depends on the universe.**

Two consequences, taken rather than worked around:

* **a per-target-only representation still answers a version of the question** — F1/F2/F4/F5 are
  built entirely from target-level ESM plus the candidate's own geometry, needing no new ESM
  compute at all;
* **F3 converts the per-candidate question into a per-target computation** — run ESM on the
  *target*, read the per-position distribution over amino acids, then score each candidate's **own
  source sequence** under it: **126 forward passes instead of 63,000**.

**Compute actually spent:** F1/F2/F4/F5 ≈ 35 s total (numpy). F3 ≈ 17 min (ESM-2-650M, one
unmasked pass plus n masked passes per target, 126 targets). F6 ≈ 19 min (ESM-2-650M over 2,369
unique band window sequences). F7 ≈ 2 min. One heavy process at a time, each run alone and
checkpointed. The box sat
at 100% CPU with 3 sibling agents throughout; at one point free RAM reached **0.06 GB** and the
650M load (2.6 GB) was **deferred behind a memory waiter** rather than launched into a thrashing
machine.

---

## 3. THE BAR, and why it must name its band

The distance objective is not the bar — it has no in-band skill (L10; PHYSICS reproduced this
independently at +0.0587 [−0.0395, +0.1606]), so beating it is meaningless.

| band | size | ORACLE best | random-in-band | constant α-helix | constant β-strand | dist argmin | n |
|---|---|---|---|---|---|---|---|
| **top-25** (native-free) | 25.0 | 2.609 | 3.502 | **3.486** | 3.639 | 3.454 | 126 |
| **top-75** (native-free) | 75.0 | 2.306 | 3.551 | **3.612** | 3.763 | 3.454 | 126 |
| `oracle1.5` (diagnostic) | 138.3 | 1.718 | 2.623 | **2.775** | 2.875 | 2.657 | 122 |

**The bar moves by band, and both figures in circulation reproduce here.** SELECT measured the
constant α-helix at ρ +0.053 in the shipped top-25; PHYSICS measured it at +0.291 selected with
ρ −0.135 on `pool_best + 1.5 Å`. On this instrument the helix is **−0.016 [−0.138, +0.104]**
against random in top-25 (indistinguishable from random *and* from the distance objective) and
**+0.151 [+0.096, +0.206]** against random on `oracle1.5` (significantly *worse* than random). Not
a conflict — two different sets. **Every number in this document names its band.**

`random-in-band` is the **exact** expectation of a uniform pick (= the band mean), not a few seeded
draws, so no sampling noise enters any paired difference.

**A units guard that fired before it could cost anything.** `core.geometry.build_backbone` takes
**radians**; degrees silently builds a different constant conformation. Measured, not inspected:
the radian α-helix has d(i,i+3) 5.20 Å, d(i,i+4) 6.38 Å, rise 1.58 Å/residue; the degree bug gives
d(i,i+3) 7.49 Å, which is neither helix nor strand. `feat_lib.check_constants` asserts the geometry
**and** that the CA-only SS rule labels each reference correctly, so it cannot regress.

---

## 4. GLOBAL vs IN-BAND — where each feature sits on the axis the sprint turns on

L4's axis, applied to the new features. A large global ρ with a small in-band ρ is a **garbage
filter**; that is what the shipped objective is.

| feature | GLOBAL ρ_S | in-band top-25 | in-band top-75 | in-band / global |
|---|---|---|---|---|
| `dist_shipped` | **0.568** | 0.048 | 0.065 | **0.12** |
| `risk_wone` / `risk_wconf` / `risk_wcon` | 0.555–0.566 | 0.025–0.032 | 0.049–0.068 | 0.09–0.12 |
| `cf_topd_uniform` *(ESM-free control)* | 0.449 | **0.105** | **0.109** | **0.24** |
| `esmcon_topd` | 0.409 | 0.098 | 0.109 | 0.27 |
| **`esmcon_agree`** | 0.373 | **0.106** | **0.116** | **0.31** |
| `env_esm` / `ss_const` | 0.34 | 0.028 / 0.067 | 0.024 / 0.026 | 0.07 / 0.08 |
| `cf_uniform` | 0.310 | 0.068 | 0.066 | 0.21 |
| `NULL_alpha` | 0.295 | 0.053 | 0.015 | 0.05 |
| `esmcon_corr` | 0.253 | 0.066 | 0.098 | **0.39** |
| `pll_mask` | 0.112 | 0.031 | 0.044 | **0.39** |
| `blosum_sim` | 0.054 | −0.004 | −0.013 | −0.23 |
| `NULL_chance` | −0.002 | −0.017 | −0.020 | — |
| `esmcon_bce` | −0.204 | −0.040 | −0.031 | 0.15 |

**Read.** The ESM-derived signals are a *different kind of object* from the shipped objective: they
have **one third to one half the global skill and roughly twice the in-band skill**, i.e. an
in-band/global ratio of 0.27–0.39 against the objective's 0.12. That is the one structural sense in
which this lane found something new. It is also, by §5, not enough to change an argmin.

`esmcon_bce` deserves a line of its own: the *principled* cross-entropy form, chosen a priori as
the right functional, has **negative global ρ (−0.204)** and is significantly worse than random
in-band (+0.202 [+0.047, +0.366] at top-75) and than the shipped score (+0.299 [+0.148, +0.455],
fold [+0.146, +0.447]). A principled aggregation is a hypothesis, not a guarantee.

---

## 5. F1 / F2 — the ESM contact head  [Problem B]

### 5.1 In-band correlation: real, and ESM-specific

Deployable **top-75** band, n = 126. `rhoS` measured; `acc` **counted directly** (null 0.500), not
converted; `rho_c` the copula equivalent via `selection_theory.rho_from_spearman`.

| signal | ρ_S | 95% CI | acc (counted) | ρ_c | selected | vs random-in-band | W/L |
|---|---|---|---|---|---|---|---|
| **`esmcon_agree`** | **+0.116** | **[+0.047, +0.184]** | **0.541** | 0.122 | 3.554 | +0.003 [−0.170, +0.176] | 59/67 |
| `esmcon_topd` | +0.109 | [+0.031, +0.185] | 0.538 | 0.114 | 3.610 | +0.059 [−0.112, +0.236] | 60/66 |
| `esmcon_corr` | +0.098 | [+0.029, +0.167] | 0.536 | 0.103 | 3.524 | −0.026 [−0.202, +0.144] | 64/62 |
| `esmcon_bce` | −0.031 | [−0.101, +0.038] | 0.487 | −0.033 | 3.753 | +0.202 [+0.047, +0.366] | 57/69 |
| `dist_shipped` *(reference)* | +0.065 | [+0.011, +0.121] | 0.525 | 0.068 | 3.454 | −0.097 [−0.210, +0.026] | 70/56 |
| `NULL_alpha` *(bar)* | +0.015 | [−0.078, +0.112] | 0.510 | 0.015 | 3.612 | +0.061 [−0.082, +0.201] | 60/66 |

### 5.2 THE CONTROL THAT DECIDES WHETHER F1 IS ESM AT ALL

A contact-probability-weighted count of realised contacts is dominated by **how compact the
candidate is**, and compactness is a channel SELECT already surveyed. This is BRIEF §5's failure
mode — *a quantity measured correctly and read as a different quantity* — and it is the one place
in this lane where it could have landed. So each F1 functional was given an **ESM-free twin**:
identical construction, uniform weights (or the distogram's own choice of contact-like pairs)
instead of ESM's contact probabilities. Paired on in-band ρ, **positive = ESM orders the band
better**:

| contrast | top-25 | top-75 | `oracle1.5` |
|---|---|---|---|
| **ρ(`esmcon_agree`) − ρ(`cf_uniform`)** | +0.038 [−0.016, +0.093] fold [+0.009, +0.072] | **+0.050 [+0.005, +0.097]** fold [+0.019, +0.076] | **+0.062 [+0.030, +0.099]** fold [+0.043, +0.089] |
| ρ(`esmcon_agree`) − ρ(`rg_raw`) | +0.074 [+0.008, +0.141] | +0.066 [+0.009, +0.124] | +0.075 [+0.035, +0.119] |
| **ρ(`esmcon_topd`) − ρ(`cf_topd_uniform`)** | −0.007 [−0.098, +0.080] | **+0.000 [−0.078, +0.080]** | −0.019 [−0.067, +0.030] |

**The verdict is split by functional, and both halves matter.** Weighting a soft contact map by
ESM's contact probabilities beats weighting it uniformly, by ~0.05 of in-band ρ, with intervals
excluding zero in two of three bands — **a genuine, ESM-specific in-band increment, the first
measured in this programme. SUPPORTED.** But choosing *which pairs to look at* by ESM's contact
head is worth **exactly nothing** over choosing them by the distogram's own predicted distances
(+0.000 [−0.078, +0.080]). So most of F1's in-band correlation is not the language model; the
ESM-specific part is the smaller half.

### 5.3 And none of it reaches the argmin

Not one F1 functional beats random-in-band, the α-helix, or the β-strand with an interval excluding
zero, in any band. **The signal exists on the correlation axis and vanishes on the decision axis.**

### 5.4 F2: an exogenous weight adds nothing

L5 closed the distogram's *own* weightings. An ESM-contact-derived weight is **exogenous** — a
different model deciding which pairs to trust — and is the one thing that sweep did not contain.

| contrast | top-25 | top-75 | `oracle1.5` |
|---|---|---|---|
| `risk_wcon − risk_wone` | +0.057 [−0.070, +0.188] | +0.021 [−0.107, +0.151] | −0.001 [−0.061, +0.059] |

Null everywhere. `risk_wshipcon` has the family's highest ρ_S (+0.088 [+0.031, +0.145] at top-75)
and a selected value of −0.049 [−0.167, +0.073] — F1's pattern again. **F2's falsifier fired.**
L5's mechanism extends: *you cannot fix an in-band problem by re-weighting a global signal, and it
does not matter whether the weight comes from inside or outside the distogram.*

---

## 6. F3 / F4 / F5 — the sequence channel, and the 0.288 Å claim  [Problem B]

### 6.1 F4 / F5 — the deciding controls

Matched architecture, matched training, matched labels; **only the per-residue feature block
changes**. Labels are native local CA descriptors / CA-only SS of **training-fold targets only**;
inference is native-free. **Positive = the ESM arm is worse.**

| contrast | top-25 | top-75 | `oracle1.5` |
|---|---|---|---|
| **`env_esm − env_onehot`** | +0.036 [−0.063, +0.139] | **+0.116 [−0.004, +0.243] fold [+0.019, +0.174]** | +0.061 [−0.011, +0.135] fold [+0.010, +0.109] |
| `env_esm − env_const` (vs no sequence) | +0.014 [−0.111, +0.131] | +0.108 [−0.013, +0.229] fold [+0.037, +0.172] | +0.017 [−0.062, +0.099] |
| `ss_esm − ss_onehot` | +0.118 [−0.035, +0.270] | −0.019 [−0.173, +0.136] | −0.046 [−0.120, +0.023] |

In every band the ESM block is at best indistinguishable from 20-dimensional one-hot, and at
top-75 significantly *worse* on the fold-clustered interval — also worse than an arm with **no
sequence input at all**. **F4's and F5's falsifiers both fired.**

### 6.2 F3 — the target-conditioned pseudo-likelihood of the window's own sequence

ESM-2-650M run on the *target* sequence with residue r masked, giving P(aa | context); each
candidate scored by the mean log-probability of **its own source sequence** under that
distribution. This is the learned-sequence question BLOSUM cannot express.

| signal | top-75 ρ_S | 95% CI | selected | vs random-in-band |
|---|---|---|---|---|
| `pll_mask` (masked, the honest form) | **+0.044** | **[+0.009, +0.078]** | 3.524 | −0.026 [−0.140, +0.086] |
| `pll_nat` (unmasked) | +0.014 | [−0.016, +0.046] | 3.602 | +0.051 [−0.064, +0.170] |
| `pll_bg` (**zero-information**: composition only) | +0.030 | [−0.006, +0.065] | 3.526 | −0.024 [−0.139, +0.091] |
| `blosum_sim` (substitution matrix) | −0.013 | [−0.044, +0.020] | 3.541 | −0.010 [−0.144, +0.126] |

| contrast | top-25 | top-75 | `oracle1.5` |
|---|---|---|---|
| `pll_mask − pll_bg` (context vs composition) | +0.040 [−0.083, +0.167] | −0.002 [−0.151, +0.144] | +0.060 [−0.024, +0.148] |
| `pll_mask − blosum_sim` | −0.031 [−0.167, +0.113] | −0.017 [−0.181, +0.152] | **+0.091 [+0.001, +0.178]** |

`pll_mask` has a small in-band ρ whose CI excludes zero at top-75 — but it **does not separate from
a position-independent amino-acid-composition control**, and on the ORACLE band it is significantly
*worse* than BLOSUM. **F3's falsifier fired.** The masked/unmasked gap (+0.044 vs +0.014) is real
and is the expected direction, so the instrument works; there is simply nothing for it to find.

### 6.3 Reconciling the recorded 0.288 Å

`s7/repr_select.py` measured ESM-PCA32 beating one-hot by **0.288 Å [−0.484, −0.092], p = 0.0046,
74W/45L** on selected RMSD at n = 126. That number is not challenged here — it was measured
**globally**, through a *retrained distance predictor* whose score then ranks the whole
500-candidate pool. **Global skill is garbage-filtering (L4/L10).** Measured *in-band*, as a
per-candidate feature, at matched architecture, the same contrast is null or negative in all three
bands.

> **ESM buys a better FILTER. It does not buy a DISCRIMINATOR.** The 0.288 Å is real, is global,
> and is the same kind of quantity as the distance objective's −0.999 Å against matched random —
> which L10 showed comes with **−0.014 [−0.118, +0.094]** of in-band skill. **ESTABLISHED.**

### 6.4 The zero-sequence arm wins its own family

At top-25 `ss_const` — the training-set SS **marginal**, no sequence input whatsoever — has
ρ_S +0.067 [+0.001, +0.135] and selected 3.401, **−0.101 [−0.209, +0.006]** against random-in-band
and **−0.085** against the α-helix, beating both the ESM arm (3.626) and the one-hot arm (3.508) of
its own family. The marginal is helix-dominated, so this is a soft α-helix preference wearing a
predictor's clothes — BRIEF §5's warning in miniature, and precisely why the α-helix, not the
distance objective, is the bar.

---

## 7. F6 — per-candidate window embeddings, the most direct form of the question  [Problem B]

The cost gate and the model ladder were both declared in `PREREG_features.md` **before** either
run. **The scout rung was not needed**: when the run was launched, free RAM was 5.37 GB, so
`esm2_t33_650M_UR50D` — *the shipped model* — was affordable and was used. No 8M numbers exist and
none are mixed in. `s17/feat_window.py` embedded the **2,369 unique source sequences** of the
shipped objective's own top-25 across all 126 targets (3,150 windows), 19 min of ESM-2-650M on a
contended box, cached to `s17/cache/feat_win_esm2_t33_650M_UR50D_B25.npz`.

Each candidate window has the **same length** as the target, so the target's per-residue embedding
matrix and the window's are **aligned by construction** — the similarity is a per-residue quantity,
not a pooled one. Valid only in the **top-25** band (that is what was embedded); reported for no
other band, and the instrument refuses to score a band-restricted feature on a partial band.

| signal | ρ_S | 95% CI | acc (counted) | selected | vs random-in-band | W/L |
|---|---|---|---|---|---|---|
| `emb_cos` (residue-aligned cosine) | +0.009 | [−0.035, +0.055] | 0.503 | 3.511 | +0.009 [−0.105, +0.116] | 65/61 |
| `emb_l2` (shipped whitened PCA-32 distance) | +0.019 | [−0.027, +0.067] | 0.507 | 3.509 | +0.007 [−0.107, +0.119] | 63/63 |
| `emb_pool` (mean-pooled cosine) | +0.019 | [−0.030, +0.068] | 0.506 | **3.405** | **−0.096 [−0.196, −0.002]** | 70/56 |
| `blosum_sim` (control) | −0.004 | [−0.056, +0.045] | 0.498 | 3.488 | −0.013 [−0.128, +0.095] | 65/61 |

**F6's falsifier fired.** `emb_cos` — the sharpest form of "does this window's sequence look like
the target's in the language model's space" — has **in-band ordering accuracy 0.503** against a null
of 0.500, does not separate from BLOSUM (+0.022 [−0.118, +0.158]) and does not beat the α-helix
(+0.025 [−0.126, +0.174]). Embedding-space sequence agreement carries **no in-band information at
all**. This is `structure-and-sequence-are-decoupled` measured on the decision axis.

**`emb_pool` and why it is not an exception.** It beats random-in-band and the β-strand, but the
pre-registered bar requires beating **both** constant conformations, and it does not beat the
α-helix:

    emb_pool vs random-in-band   -0.096 [-0.196,-0.002]  fold [-0.153,-0.045]  med -0.027  70W/56L  dtop +0.009
    emb_pool vs ALPHA-HELIX      -0.081 [-0.238,+0.080]  fold [-0.237,+0.115]  med +0.000  59W/59L  dtop +0.086
    emb_pool vs beta-strand      -0.234 [-0.406,-0.060]  fold [-0.343,-0.116]  med -0.133  76W/48L  dtop -0.061
    emb_pool vs BLOSUM           -0.083 [-0.241,+0.069]  fold [-0.178,+0.028]  med -0.009  65W/53L  dtop +0.080

Three independent reasons it is not a result, all of which the pre-registration anticipated:
**(i)** 59W/59L and a **median of exactly +0.000** against the α-helix — a dead heat;
**(ii)** its **drop-top-10 mean is +0.009**, i.e. the entire 0.096 Å is carried by about ten
targets, which BRIEF §4 says is not a result; **(iii)** its **in-band ρ is +0.019 and its counted
ordering accuracy is 0.506** — it has no ordering skill, so there is no mechanism for the gain, and
a selected-RMSD win with no rank correlation behind it is what luck looks like. **NULL.**

That the *pooled* form scores better than the *residue-aligned* form while having the same
(absent) ordering skill is itself the tell: pooling discards exactly the positional information a
real sequence-structure signal would need.

---

## 8. THE CEILING OF THE SIGNAL THAT DOES EXIST — this lane's most portable number

Per the brief: *if you find in-band signal, the immediate follow-up is its ceiling, not a model.*
Measured in-band ρ fed through `selection_theory`'s band curve at K = 500, n = 126:

| measured in-band ρ_S | ρ_c | acc | B = 5 | B = 10 | **B = 25** | B = 50 | B = 100 |
|---|---|---|---|---|---|---|---|
| **+0.116** `esmcon_agree` (best in programme) | 0.121 | 0.539 | 3.434 | 3.406 | **3.393** | 3.394 | 3.387 |
| +0.131 (SELECT's figure for the objective at K = 500) | 0.137 | 0.544 | 3.428 | 3.395 | 3.379 | 3.376 | 3.365 |
| +0.065 `dist_shipped` measured here | 0.068 | 0.522 | 3.456 | 3.439 | 3.441 | 3.457 | 3.465 |
| band mean (random-in-band) | 0 | 0.500 | 3.484 | 3.478 | 3.502 | 3.539 | 3.567 |
| band best (ORACLE) | 1 | 1.000 | 3.025 | 2.834 | 2.609 | 2.403 | 2.210 |

**The best in-band signal measured in the deployable band is worth ≈0.11 Å through an argmin**,
and the same curve read at SELECT's and PHYSICS's larger in-band figures moves it by hundredths. The requirement
side is worse than "hard":

| shortlist B | band ORACLE best | ρ for 2.5 Å | ρ for 2.2 Å | ρ for 2.0 Å |
|---|---|---|---|---|
| 5 | 3.025 | unreachable | unreachable | unreachable |
| 10 | 2.834 | unreachable | unreachable | unreachable |
| **25** | 2.609 | **unreachable** | unreachable | unreachable |
| 50 | 2.403 | ρ 0.95 / acc 0.899 | unreachable | unreachable |
| 100 | 2.210 | ρ 0.80 / acc 0.795 | unreachable | unreachable |

> **At B ≤ 25 the sprint's 2.5 Å target is unreachable at ρ = 1.** A perfect in-band ranker returns
> 2.609 Å because that is the band's best member. This is not a statement about ranking skill; it
> is arithmetic on the shortlist's contents. **The top-25 shortlist, not the ranker inside it, is
> the binding constraint** — which is L12/L14's conclusion arriving from a third direction.
> **EXACT.**

The measured selected values also sit **below** their own curve (`esmcon_corr` returns 3.524 where
ρ_c = 0.103 predicts ≈3.40), which by `selection_theory`'s own stated caveat means the errors are
**structured, not exchangeable** — the signal is not merely weak, it is weak in a patterned way.

### 8.1 The candidate finding, labelled honestly

`cf_topd_uniform` — mean realised CA distance over the **n pairs the distogram predicts closest** —
was built **only as F1's ESM-free control** and is the battery's best arm:

| band | ρ_S | selected | vs random-in-band | vs α-helix | vs β-strand | W/L |
|---|---|---|---|---|---|---|
| top-25 | +0.105 [+0.038, +0.170] | 3.414 | −0.087 [−0.181, **+0.000**] | −0.072 | −0.225 | 70/56 |
| **top-75** | **+0.109 [+0.041, +0.177]** | **3.415** | **−0.135 [−0.252, −0.022]** | **−0.196** | **−0.348** | 72/54 |
| `oracle1.5` | +0.069 [+0.010, +0.128] | 2.678 | +0.054 [+0.005, +0.102] *(worse)* | −0.097 | −0.198 | 50/72 |

It beats the shipped objective (3.454) and every representation feature. **It is not a
representation feature** — it is a hard top-k functional of the distogram's predicted centres, a
class L5 nominally closed. Against it: it is **one arm of 25 across 3 bands with no multiplicity
correction**, its top-25 interval **touches zero**, it **reverses sign on the ORACLE band**, and it
is a post-hoc control rather than a pre-registered hypothesis. **PLAUSIBLE, not ESTABLISHED.**
Recorded here so SELECT or AUDIT can confirm or kill it on an independent instrument; **no
architectural claim is made on it.**

### 8.2 Stratification (BRIEF §4) — the ESM contact signal is length-gated

A mean null can hide a real effect on one stratum, exactly as a mean improvement can be carried by
easy targets. `esmcon_agree`, top-25 band:

| stratum | n | ρ_S | selected vs random | W/L |
|---|---|---|---|---|
| length 9–11 | 33 | **−0.047** | **+0.203** | 13/20 |
| length 12–13 | 42 | +0.127 | +0.131 | 17/25 |
| **length 14–16** | **51** | **+0.188** | **−0.271** | 28/23 |
| pool best < 2.0 Å (easy) | 73 | +0.108 | +0.059 | 32/41 |
| pool best ≥ 2.0 Å (hard) | 53 | +0.104 | −0.112 | 26/27 |
| folds 0–4 | 23–30 | +0.066 … +0.135 | −0.143 … +0.195 | — |

**The ESM contact channel has in-band content only at length ≥ 14 and is actively harmful at
length ≤ 11** — mechanistically sensible, since a 9-mer has almost no tertiary contacts for a
contact head to predict. Length 14–16 is also where SELECT §9 puts 31.4% of the residual hard
targets. This is a **post-hoc stratification at n = 51 with a W/L of 28/23** — **PLAUSIBLE**, and
the one place in this lane where a targeted follow-up could still be worth someone's compute.

### 8.3 Min-of-N (LEDGER L8), reported because a best-of-V is otherwise unreadable

Real features, V-matched zero-information null, per band (the ladder's nearest rung is used
when the feature count is not on it, and the rung is named):

| band | ORACLE min over 25 real features | zero-info min over 25 random picks | real − zero-info |
|---|---|---|---|
| top-25 (28 real features vs a V = 30 null) | 2.695 | 2.675 | **+0.020 [−0.011, +0.053]** |
| top-75 (25 vs 25) | 2.468 | 2.530 | **−0.061 [−0.105, −0.022]** |
| `oracle1.5` (25 vs 25) | 2.012 | 1.928 | **+0.083 [+0.039, +0.129]** |

**The sign flips with the band.** SELECT measured +0.059 [+0.020, +0.100] (real loses); the
coordinator measured −0.032 [−0.067, +0.001] (touching zero); this lane gets all three signs on one
instrument by changing only the band. **No routing headroom is claimed, and the honest reading is
that this statistic is not stable enough to support a claim in either direction.** L8's retraction
stands, and gains a third refinement: *a min-of-N contrast must name its band too.*

---

## 9. F7 — is the ESM channel a SHORTLIST CONSTRUCTOR instead?  [Problems B and C]

F1 found rank skill with no argmin skill. L12/L14 identify an unoccupied role for exactly such a
signal — **choosing the shortlist** rather than choosing inside it — and the coordinator flagged the
same opening independently. `s17/feat_shortlist.py`, n = 126, K = 500. `ceiling` is the shortlist's
ORACLE best; `mix` takes B/2 from each score and is **topped up to EXACTLY B** so a union arm is
never priced against a larger shortlist.

| B | arm | **CEILING** | mean | argmin | medoid | P(sub-2 Å) |
|---|---|---|---|---|---|---|
| 25 | `dist` (incumbent) | 2.609 | 3.502 | **3.454** | 3.369 | 0.373 |
| 25 | `esmcon` | 2.849 | 4.248 | 3.776 | 3.921 | 0.341 |
| 25 | **`rand`** (matched random) | **2.400** | 4.498 | 3.527 | 3.819 | 0.389 |
| 25 | `mix` | 2.369 | 3.913 | **3.454** | 3.405 | 0.405 |
| 75 | `dist` | 2.306 | 3.551 | **3.454** | **3.282** | 0.437 |
| 75 | `esmcon` | 2.410 | 4.028 | 3.568 | 3.589 | 0.444 |
| 75 | **`rand`** | **2.080** | 4.453 | 3.500 | 3.707 | 0.444 |
| 75 | `mix` | 2.100 | 3.872 | **3.454** | 3.450 | 0.484 |
| 150 | `dist` | 2.106 | 3.606 | **3.454** | **3.383** | 0.476 |
| 150 | `esmcon` | 2.250 | 3.945 | 3.528 | 3.498 | 0.476 |
| 150 | **`rand`** | **1.935** | 4.451 | 3.543 | 3.724 | 0.484 |
| 150 | `mix` | 1.985 | 3.834 | **3.454** | 3.421 | 0.508 |

### 9.1 F7's own falsifier fired — the ESM shortlist is worse, not better

| contrast on the CEILING | B = 25 | B = 75 | B = 150 |
|---|---|---|---|
| `esmcon` vs `dist` | +0.241 [−0.013, +0.483] | +0.104 [−0.105, +0.314] | +0.145 [−0.036, +0.336] |
| **`esmcon` vs matched random** | **+0.449 [+0.172, +0.745]** | **+0.330 [+0.120, +0.558]** | **+0.315 [+0.130, +0.525]** |

An ESM-contact-ranked shortlist is no better than the incumbent's and **significantly worse than a
random shortlist of the same size**. Its realized medoid is worse than the incumbent's at every
size (+0.552 / +0.307 / +0.115). **REFUTED: the ESM contact channel has no shortlist-construction
value either.**

### 9.2 An independent reproduction of the PHYSICS shortlist result

| CEILING: `dist` vs matched random | B = 25 | B = 75 | B = 150 |
|---|---|---|---|
| this instrument | **+0.209 [+0.021, +0.402]** | **+0.226 [+0.083, +0.377]** | **+0.170 [+0.038, +0.310]** |
| PHYSICS, separately written | +0.195 | +0.242 | +0.179 |

**The shipped score's top-B shortlist has an ORACLE best significantly WORSE than a matched random
shortlist of the same size, at every size, on two independently written instruments agreeing to
0.02 Å.** This is L12's mechanism measured directly at fixed K rather than across K. **ESTABLISHED.**

### 9.3 The one positive arm, killed by its own pre-registered control

`mix` beats the incumbent on the ceiling at every size, with **both** intervals excluding zero:

    mix vs dist, CEILING:  B=25  -0.239 [-0.389,-0.108]  fold [-0.436,-0.033]
                           B=75  -0.206 [-0.334,-0.093]  fold [-0.365,-0.057]
                           B=150 -0.121 [-0.220,-0.032]  fold [-0.199,-0.037]

and P(a sub-2 Å candidate is present) rises 0.373 → 0.405, 0.437 → 0.484, 0.476 → 0.508.

**And the matched-random control says it is nothing.** `rand` reaches 2.400 / 2.080 / 1.935 against
`mix`'s 2.369 / 2.100 / 1.985 — a random shortlist of the same size does **as well or better at
every size**. The `mix` gain is not a diversification discovery; it is *partial movement toward
random*, and random is the endpoint. Without the matched-random arm this would have been written up
as a positive result. **It was pre-registered as mandatory for exactly this reason, and it fired.**

### 9.4 And the whole ceiling axis is unconsumable, as L14 predicted in advance

**The realized argmin is 3.454 for `dist`, `mix` at every B — identical to three decimals** —
because the argmin is score-driven and the score's own favourite is in every shortlist. The medoid
*degrades* whenever the shortlist ceiling improves (`mix` vs `dist`: +0.036 / +0.168 / +0.038;
`rand` has the best ceiling and the worst medoid at B = 25 and 75). Ceiling and readout move in
opposite directions because the readouts consume the shortlist **mean**, which diversification
destroys, while the ceiling is the **minimum**, which diversification improves.

> **This lane's contribution to Problem C: the shortlist-ceiling axis is now measured with a
> matched-random control at fixed K, and the control says the incumbent score is the WORST
> shortlist constructor tested while being the only one whose shortlist a readout can use.**
> L14's *"a ceiling gain is worthless unless a realized readout can use it"* is reproduced here on
> a third construction. **Problems B and C must be solved together or not at all.**

---

## 10. What I believe follows — inference from the numbers above, labelled as inference

1. **The sprint's selection answer is final, and the closure is now by measurement on four
   independent feature classes.** 58 distogram functionals (SELECT E2), 29 native-free in-band
   signals (SELECT §5), 27 target-level calibration features (SELECT E4), and — this lane — six
   classes of learned sequence representation. All land at zero on the decision axis. This is not
   "we did not find it"; it is four measured closures with matched controls attached.

2. **Stop asking for a better in-band ranker at B ≤ 25. It cannot work, and the reason is
   arithmetic, not skill.** A *perfect* in-band ranker over the shipped top-25 returns 2.609 Å.
   The sprint's 2.5 Å target is unreachable there at ρ = 1. Any further reranking work must first
   name a shortlist whose ORACLE best clears the target — which means the shortlist, not the
   ranker, is the object to change.

3. **But do not change it by diversifying, because F7 measured where that goes.** Every shortlist
   construction that improves the ceiling degrades the readout, and the endpoint of the
   ceiling-improving direction is *random*, whose realized medoid is the worst of the four arms.
   L14 said "Problems B and C must be solved together or not at all"; F7 is the third independent
   measurement saying so, and the first with a matched-random control at fixed K.

4. **The 0.288 Å ESM result should be re-labelled in the record, not retracted.** It is real and it
   is *global*. The record currently reads as though the sequence representation carries selection
   signal; the accurate statement is that it carries **filtering** signal, of the same kind and in
   the same place as the distance objective's. The practical consequence recorded with it — that
   one-hot arms are handicapped by ~0.29 Å — is unaffected. What must not survive is any inference
   from that number to *in-band* discrimination, because this lane measured that directly and it is
   null or negative.

5. **The one direction this lane leaves genuinely open is narrow, and I would not spend much on
   it.** The ESM contact channel's in-band content is **length-gated**: ρ_S +0.188 with selected
   −0.271 at length 14–16 (n = 51), against ρ_S −0.047 and +0.203 at length 9–11 (n = 33). That is
   post-hoc, uncorrected, and 28W/23L. If anyone follows it, the honest framing is a *pre-registered
   replication on length ≥ 14 only*, and §8's ceiling says the prize even if it replicates is
   ≈0.1–0.3 Å on a third of the targets — not a route to 2.0 Å.

6. **`cf_topd_uniform` (§8.1) deserves ten minutes of someone else's instrument, and nothing
   more until it survives that.** It is the best arm in the battery, it beats both zero-information
   conformations and the shipped objective at top-75 with a CI excluding zero, and it is a
   *distogram* functional — so if it replicates it belongs to SELECT's E2 lane, not this one, and
   it would mean L5's 58-variant sweep missed a hard top-k functional. If it does not replicate it
   is one uncorrected arm out of 75 tests, which is exactly what it looks like.

---

## 11. Limitations, stated rather than buried

- Every number is on the 126-target **tuning** instrument. The sealed benchmark was not touched.
- `oracle1.5` needs the native to draw and is a **diagnostic** band; it drops to n = 122 where the
  band has fewer than 8 members. Every deployable claim uses the native-free top-B bands, and the
  report filters **per band** so the target set never changes silently between tables.
- The fold-clustered bootstrap resamples **five** clusters; the intervals are correspondingly wide.
  Both intervals are printed everywhere and neither is suppressed. Where they disagree
  (`env_esm − env_onehot` at top-75: target [−0.004, +0.243], fold [+0.019, +0.174]) both are
  quoted and the weaker one governs the claim.
- **25–28 features × 3 bands (75+ tests) with no multiplicity correction.** Every individually-significant result
  in this document is flagged with that fact attached; the lane's *negative* conclusion is
  unaffected by multiplicity, and its one positive candidate (§8.1) is labelled PLAUSIBLE for
  exactly this reason.
- F4/F5's secondary structure is a **CA-only proxy** for DSSP (native traces are CA-only), applied
  by one fixed rule to native and candidate alike so the two are commensurable. It is not DSSP and
  is never called DSSP.
- F4/F5 train on ~1,300 residues from 100 training-fold targets — small. A larger corpus could
  change the *absolute* quality of those heads; it cannot change the ESM-vs-one-hot **contrast**,
  which is what the falsifier is stated against, since both arms see the same corpus.
- ESM-2's **contact head is supervised on structure** (of other proteins). It is the one component
  here that is not purely a language model, and it is flagged at every use.
- ESM contacts are an **input feature of the shipped distogram** (`core/predict.py:144,612`), so F1
  is not an independent channel from the objective — it is the same input consumed differently.
- `esmcon_agree` vs `dist_shipped` on in-band ρ is **band-dependent** (+0.051 at top-75 with the
  target-level CI touching zero; **−0.061 [−0.124, +0.000]** on `oracle1.5`, i.e. the trained
  distogram wins there). **No claim is made that the raw contact map beats the trained predictor.**
