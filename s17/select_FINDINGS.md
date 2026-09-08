# SELECT workstream — Sprint 17 findings

Pre-registration: `s17/PREREG_select.md`, written before any number below.
Modules: `s17/sel_lib.py` · `sel_obj.py` · `sel_bench.py` · `sel_inband.py` · `sel_cal.py` · `sel_hard.py`
Artefacts: `s17/results/sel_obj_K500.{json,log}` · `sel_bench.{json,log}` · `sel_inband.{json,log}` ·
`sel_cal.{json,log}`

All numbers: n = 126 tuning targets, target as the unit, identical candidate set per (target, K),
tie-safe selection, both a matched-random and a zero-information control on every arm. The sealed
60-target benchmark was not read, probed, or derived.

---

## 0. The four things that change the sprint

1. **Widening the candidate set makes the realized answer WORSE, not better.** From top-75 to the
   full universe the ORACLE ceiling improves by 0.791 Å and the distance selector's realized RMSD
   **degrades** by 0.141 Å [+0.000, +0.291], 51W/68L. The fraction of extra ceiling recovered is
   **−8.4% / −16.4% / −17.9%** at K = 500 / 2000 / full. Recall is not the blocker, and the pivot is
   stronger than "ranking is the blocker": ranking *deteriorates* as recall improves. **ESTABLISHED.**

2. **The mechanism, in one line.** As K grows the selector's global rank correlation *rises*
   (0.530 → 0.587) while its **in-band** correlation *falls* (0.200 → 0.090). The distance
   objective's skill is almost entirely garbage-versus-plausible and nearly absent where the argmin
   is decided. **ESTABLISHED.**

3. **58 functionals of the existing distogram are worth −0.009 Å [−0.128, +0.110] leave-fold-out.**
   Pair weighting (sequence separation, uncertainty, entropy), residual form (Bayes risk, NLL,
   absolute, z, log, squared), robust aggregation (Huber, Tukey, median, trimmed, CVaR),
   set-conditional standardisation and rank-normalisation, contact-vs-distance mixtures, a
   distance-objective↔typicality blend, and a spectral form — none survives fold-honest selection.
   With (2) attached this is not a disappointment but a mechanism: **you cannot fix an in-band
   problem by re-weighting a global signal.** *"The distance objective can be re-engineered into a
   better selector"* is **REFUTED**.

4. **No native-free signal available to this programme has usable in-band discrimination.**
   29 signals × 3 band definitions × n = 126: the best in-band Spearman in the deployable top-75
   band is **+0.086** (pairwise-ordering accuracy 0.529). The recorded requirement for 2.0 Å is
   accuracy **0.638**. A zero-information constant α-helix reaches +0.053 in the top-25 band,
   statistically indistinguishable from the distance objective's +0.048 there. **ESTABLISHED.**

And one methodological result that invalidated two candidate headline claims — one mine, one the
coordinator's — before either was published: **§6, the min-of-N trap.**

---

## 1. Which pre-registered rules fired

| experiment | rule | fired? |
|---|---|---|
| E1 identical-candidate instrument | falsifier: `sel-dist − sel-random` CI contains zero | **did NOT fire** — −0.999 [−1.194, −0.799], 102W/24L. The distance objective transfers. |
| E2 objective decomposition | falsifier: LFO variant's advantage over shipped has a CI containing zero | **FIRED** — −0.009 [−0.128, +0.110] |
| E3 set-conditional inference | falsifier: set conditioning adds nothing over the score alone | **FIRED** (partially, see §4) |
| E4 target-level calibration | falsifier: calibrator does not beat the **constant** baseline | **FIRED** — +0.037 to +0.083 *worse* than constant at every λ |
| E4 boundary trap | a selector pinned at its ladder boundary is a substitution | **FIRED** — best λ is the largest on the ladder (λ = 100, i.e. maximal shrinkage toward the constant) |
| E5 hard targets | falsifier: the hard set is unchanged by the full universe | **did NOT fire** — 53 → 24 |
| E6 does ceiling reach the selector | symmetric; either answer redirects | **fired in the pessimistic direction, strongly** |
| readout triangle (§5A, added on request) | not pre-registered; reported as DISCOVERED | averaging dominates both selection readouts in 39/39 cells |

### 1.1 Independent cross-check of the instrument

`s17/sel_bench.py` (this workstream) and `s17/oracle_map.py` (the coordinator's) are separately
written and share only `s12/instrument.py`. On the 120 targets both had completed, they agree to
**3.2e−7 Å** on both the ORACLE ceiling and the realized distance selector at all four widths:

| K | ORACLE (map / mine) | dist (map / mine) | max per-target deviation |
|---|---|---|---|
| 75 | 2.0859 / 2.0859 | 3.4000 / 3.4000 | 2.3e−7 |
| 500 | 1.7073 / 1.7073 | 3.4301 / 3.4301 | 2.3e−7 |
| 2000 | 1.5026 / 1.5026 | 3.5017 / 3.5017 | 2.4e−7 |
| full | 1.3141 / 1.3141 | 3.5536 / 3.5536 | 3.2e−7 |

This is a reproduction, not a transcription audit — BRIEF §5's closing rule is that a transcription
audit cannot see a formula error, and two independent implementations agreeing is the thing that
can.

---

## 2. E1 — the identical-candidate-set instrument  [Problem B]

`s17/sel_bench.py`. Every arm in a row consumes the same structures. ORACLE columns read the native.

### 2.1 Ceiling / realized / gap

| K | ORACLE | dist | GAP | legacy | consensus | random | ρ_within | ρ_band |
|---|---|---|---|---|---|---|---|---|
| 75 | 2.104 | 3.421 | +1.317 | 4.144 | 3.598 | 4.287 | 0.530 | 0.200 |
| 500 | **1.711** | 3.454 | +1.743 | 4.490 | **3.282** | 4.453 | 0.568 | 0.131 |
| 2000 | 1.504 | 3.520 | +2.016 | — | 3.344 | 4.596 | 0.581 | 0.098 |
| full | **1.313** | 3.562 | +2.249 | — | 3.461 | 4.816 | 0.587 | 0.090 |

The pinned constants reproduce exactly (K = 500 ORACLE 1.711, top-75 2.306 — see §7).

### 2.2 Every arm against its own matched control, K = 500

| arm vs control | means | Δ [95% CI] | fold-clustered CI | W/L |
|---|---|---|---|---|
| dist vs matched random | 3.454 / 4.453 | **−0.999 [−1.194, −0.799]** | [−1.118, −0.892] | 102/24 |
| consensus vs **matched random consensus** | 3.282 / 3.726 | **−0.444 [−0.672, −0.231]** | [−0.550, −0.290] | 80/46 |
| consensus vs dist | 3.282 / 3.454 | **−0.172 [−0.322, −0.031]** | [−0.282, −0.076] | 74/47 |
| legacy vs matched random | 4.490 / 4.453 | +0.036 [−0.276, +0.353] | [−0.277, +0.272] | 71/55 |
| legacy vs dist | 4.490 / 3.454 | **+1.036 [+0.694, +1.384]** | [+0.828, +1.220] | 38/88 |
| legacy_gate25 → dist vs **random gate of the same count** | 3.416 / 3.450 | −0.033 [−0.093, +0.024] | [−0.074, +0.001] | 70/56 |
| legacy_gate50 → dist vs random gate | 3.462 / 3.463 | −0.001 [−0.096, +0.096] | [−0.074, +0.078] | 61/65 |
| legacy_gate75 → dist vs random gate | 3.574 / 3.472 | +0.102 [−0.043, +0.257] | [−0.031, +0.222] | 64/62 |

**Reads.**
- The distance objective is a real selector: it beats matched random by 1 Å with a fold-clustered
  CI excluding zero, at every width. **ESTABLISHED.**
- **Legacy as a ranker is null against random** (+0.036 [−0.276, +0.353]) and **decisively worse
  than the distance objective** (+1.036). Confirms the standing result at n = 126 on an identical
  candidate set for the first time.
- **Legacy as a gate is a null, not a positive.** Against its own matched random gate the effect is
  −0.033 [−0.093, +0.024] at the 25% gate and shrinks to nothing at 50% and reverses at 75%. The
  one Legacy use the record actively supported does not survive its matched control. **NULL.**
- **Consensus is the only arm that beats both its matched control and the distance objective**, and
  it reproduces Sprint 8's −0.172 [−0.316, −0.027] to three decimals on a rebuilt instrument
  (−0.172 [−0.322, −0.031]). Independent reproduction. **SUPPORTED.**
- A control that must be read carefully: at K = 75 `consensus` and `rand_consensus` are identical by
  construction (a random 75 of 75 is the same set) and the instrument correctly reports +0.000
  [+0.000, +0.000], 0W/0L. That is a self-check, not a result.

### 2.3 An ORACLE yardstick worth having: how many random draws is the objective worth?

`bestof_m` = the best of m random draws, **picked with the native**. It is an ORACLE operation and
is a *currency*, not a competitor.

| K | dist | best-of-5 random (ORACLE) | Δ |
|---|---|---|---|
| 75 | 3.421 | 2.976 | +0.445 [+0.266, +0.629] |
| 500 | 3.454 | 3.068 | +0.386 [+0.209, +0.572] |
| full | 3.562 | 3.346 | +0.216 [+0.009, +0.435] |

**The shipped distance objective's entire discriminative content is worth less than five random
draws under an oracle.** That is the honest scale of the selection signal this programme has.

---

## 3. E6 — does the extra ceiling reach the selector?  [Problems A vs B]

| K | ceiling gained vs top-75 | realized gained | recovered |
|---|---|---|---|
| 500 | −0.393 | **+0.033** [−0.078, +0.152] | **−8.4%** |
| 2000 | −0.600 | **+0.099** [−0.022, +0.226] | **−16.4%** |
| full | −0.791 | **+0.141** [+0.000, +0.291] | **−17.9%** |

**The recovered fraction is negative at every width.** The programme's early 6-target read of
"1–13% recovered" was optimistic; at n = 126 it is negative.

**The copula gap grows with K**, monotonically:

| K | realized | copula prediction from the measured global ρ | gap | above curve |
|---|---|---|---|---|
| 75 | 3.421 | 2.982 | +0.439 | 96/126 |
| 500 | 3.454 | 2.765 | +0.689 | 106/126 |
| 2000 | 3.520 | 2.647 | +0.872 | 104/126 |
| full | 3.562 | 2.541 | **+1.021** | 115/126 |

Fed the *in-band* ρ instead, the curve overshoots the other way (K = 500: predicts 3.722 against a
realized 3.454, −0.268, 51/126 above). **The two predictions bracket the truth**, so the operative
ρ lies between 0.137 and 0.586 and much nearer the in-band end. The selector's errors are
**structured, not exchangeable**, and the structure is now identified: its skill lives outside the
band and the argmin is decided inside it.

### 3.1 WHY widening hurts — the mechanism, isolated

The obvious guess is that a wider set is simply harder to rank. It is worse than that: **widening
the candidate set removes the answer from the shortlist.**

| K | ORACLE over the whole set | **ORACLE inside the score's own top-75** | mean of that top-75 | argmin | consensus |
|---|---|---|---|---|---|
| 75 | 2.104 | 2.104 | 4.287 | 3.421 | 3.598 |
| 500 | 1.711 | 2.306 | 3.551 | 3.454 | 3.282 |
| 2000 | 1.504 | 2.406 | 3.526 | 3.520 | 3.344 |
| full | **1.313** | **2.572** | 3.557 | 3.562 | 3.461 |

    top-75 ORACLE, full vs K=75:  +0.468 [+0.304, +0.642]  fold [+0.255, +0.652]  46W/77L
    top-75 MEAN,   full vs K=75:  -0.730 [-0.906, -0.560]  fold [-0.844, -0.623]  104W/22L
    argmin,        full vs K=75:  +0.141 [+0.000, +0.291]  fold [+0.046, +0.253]  51W/68L
    consensus,     full vs K=75:  -0.138 [-0.381, +0.102]  fold [-0.324, +0.142]  70W/56L

**As the set widens, the score's shortlist gets better on average and worse at its best.** The extra
windows are more plausible-looking without being nearer the native, and they *displace* the genuinely
near-native members out of the top-75. The shortlist ceiling degrades by 0.468 Å — more than the
0.141 Å of realized degradation — so shortlist destruction accounts for the whole effect and more.

This is Sprint 16's **set-mean trap** (G4: Legacy's filter improved the set mean 3.551 → 3.531 while
destroying the set best 2.306 → 2.416) reproduced on a completely different axis: the *retrieval
width* axis rather than a filter's aggressiveness. It is the same score behaving the same way.

**Consensus is the exception**: it *improves* with K (−0.138), because it consumes the set mean,
which is exactly the quantity that improves. Argmin consumes the set best, which is the quantity
that is destroyed. **The readout and the filter must be matched to each other**, and the current
architecture pairs a best-consuming readout with a mean-improving filter.

**ESTABLISHED.** Actionable consequence: at wide K the shortlist must not be the score's top-B. A
diversity-preserving or coverage-preserving shortlist is the only way the wide-K ceiling could ever
be reached, and that is a Problem-C (ensemble construction) question, not a Problem-B one.

### 3.2 A units correction that was load-bearing for the sprint's design equation

The programme's recorded in-band ordering figures (0.986 within target, 0.600 across, 0.638
required) are **pairwise ordering accuracies with a null of 0.500**, not correlations. They must be
converted before being read onto a copula-ρ axis:

    Kendall tau = 2*accuracy - 1;   Gaussian copula rho = sin(pi*tau/2)
    from Spearman:                  rho = 2*sin(pi*rho_S/6)

Verified numerically against a 400,000-sample Gaussian copula, agreeing with the closed form to
three decimals.

| recorded accuracy | copula ρ |
|---|---|
| 0.539 (Legacy in-band) | 0.122 |
| **0.600 (cross-target in-band)** | **0.309** |
| **0.638 (requirement for 2.0 Å)** | **0.420** |
| 0.986 (within target) | 0.999 |

Reading 0.600 straight onto the ρ axis overstates the programme's position by a very large margin.
Reported to the coordinator, who retracted the affected strategic read and fixed it at source.
Every ρ in this document is a Spearman unless it says otherwise, and accuracy equivalents are
printed beside it in `sel_inband.py`.

---

## 4. E2 — the distance objective taken apart  [Problem B]

`s17/sel_obj.py`, 58 variants, identical candidate set, K = 500.

### 4.1 The headline

| | Å |
|---|---|
| ORACLE ceiling (best member) | 1.711 |
| in-sample best of 58 (`risk_rankset`) | 3.360 — **OVERFIT UPPER BOUND, not a result** |
| **leave-fold-out variant selection** | **3.445** |
| shipped objective | 3.454 |
| matched random | 4.453 |

    LFO vs shipped   -0.009 [-0.128, +0.110]  fold [-0.052, +0.022]  median +0.000  33W/40L
    LFO vs random    -1.008 [-1.216, -0.796]  fold [-1.132, -0.916]  101W/25L
    shipped vs random -0.999 [-1.194, -0.799] fold [-1.118, -0.892]  102W/24L

**E2's falsifier fired.** The space of functionals of the existing distogram is closed for
**selection** as well as for generation. This extends `score-axis-does-not-transfer` — measured on
24 dev targets in the *generation* setting — to *selection* at n = 126 with fold-honest variant
choice. **REFUTED: "the distance objective can be re-engineered into a better selector."**

The in-sample-to-LFO gap (3.360 → 3.445, 0.085 Å) is the overfitting a best-of-58 screen
manufactures, measured rather than assumed.

### 4.2 What the shipped weighting is actually worth

The shipped objective is an L1 Bayes risk weighted by 1/(sd + 0.5), with **all five sequence-
separation shell weights equal to 1.0 and γ = 1.0** (read out of `core/predict.py`, not assumed).

| ablation | Å | vs shipped | fold CI |
|---|---|---|---|
| remove uncertainty weighting (`risk_one`) | 3.528 | +0.074 [−0.031, +0.181] | **[+0.035, +0.110]** |
| sharpen it (1/(sd+0.5)²) | 3.459 | +0.005 [−0.079, +0.088] | [−0.050, +0.064] |
| soften it (1/√(sd+0.5)) | 3.429 | −0.025 [−0.099, +0.040] | [−0.066, +0.006] |
| entropy weighting instead | 3.503 | +0.049 [−0.063, +0.165] | [−0.026, +0.119] |

The uncertainty weighting **earns its place** — removing it costs 0.074 Å with a fold-clustered CI
excluding zero — and it is already close to optimal. **ESTABLISHED, narrow.**

Sequence-separation weighting, never fitted in the shipped code, is worth nothing: `risk_sep`
+0.089, `risk_isep` +0.039, long-range only +0.133, local only +0.239 — all worse.

### 4.3 The centre-blend family (the one novel construction)

`score_β` replaces the predicted centre by a blend with the candidate set's own mean profile;
β = 0 is the distance objective in residual form and β = 1 is typicality.

| β | 0.00 | 0.15 | 0.30 | 0.45 | 0.60 | 0.75 | 0.90 | 1.00 |
|---|---|---|---|---|---|---|---|---|
| RMSD (pred scale) | 3.436 | **3.423** | 3.497 | 3.572 | 3.616 | 3.648 | 3.655 | 3.806 |
| ρ_band | 0.150 | 0.136 | 0.106 | 0.062 | 0.007 | −0.033 | −0.068 | −0.087 |

Monotone degradation with β and no interior optimum worth having (β = 0.15 is −0.031 [−0.132,
+0.070], a null). **The typicality endpoint is not a usable objective**, and mixing toward it
destroys in-band ρ. **NULL.**

### 4.4 A class of terms that is identically zero, reported as a result

Triangle-inequality and chain-connectivity consistency terms **cannot rank this candidate set**:
every candidate is a real protein backbone window and satisfies both exactly by construction. Such
a term is identically zero here. It would matter only for a generator emitting distance matrices
directly. **EXACT** (a property of the candidate set, not a discovery).

---

## 5. The in-band discrimination survey — the sharpest result  [Problem B]

`s17/sel_inband.py`. 29 native-free signals, three band definitions, n = 126, K = 500.

Band definitions, both reported because they are different objects:
- **`oracle1.5`** — members within best + 1.5 Å. ORACLE-DEFINED, a **diagnostic** band; it is what
  the argmin decision spans and it is the definition the recorded 0.600/0.638 use.
- **`top25` / `top75`** — the shipped objective's own shortlist. NATIVE-FREE, what a deployed
  reranker actually sees.

### 5.1 The deployable band (top-75): band ORACLE best 2.306, band random 3.551

| signal | in-band ρ | accuracy | 95% CI on ρ | selected | vs band-random |
|---|---|---|---|---|---|
| `dist_z` | **+0.086** | 0.529 | [+0.031, +0.141] | 3.426 | −0.124 [−0.226, −0.019] |
| `dist_blend15` | +0.079 | 0.526 | [+0.026, +0.134] | 3.423 | −0.127 [−0.252, +0.001] |
| `typicality_med` | +0.071 | 0.524 | [−0.015, +0.160] | 3.597 | +0.046 |
| `dist_shipped` | +0.065 | 0.522 | [+0.011, +0.121] | 3.454 | −0.097 [−0.210, +0.026] |
| `dist_rankset` | +0.054 | 0.518 | [+0.010, +0.098] | **3.360** | **−0.191 [−0.311, −0.068]** |
| `leg_torsion` (best Legacy term) | +0.044 | 0.515 | [−0.004, +0.093] | 3.420 | −0.131 [−0.248, −0.017] |
| **`NULL_beta`** (zero-information) | −0.042 | 0.486 | [−0.143, +0.060] | 3.763 | +0.212 |
| **`NULL_randwin`** (zero-information) | +0.018 | 0.506 | [−0.072, +0.106] | 3.629 | +0.079 |
| **`NULL_alpha`** (zero-information) | +0.015 | 0.505 | [−0.078, +0.112] | 3.612 | +0.061 |
| **`NULL_chance`** | +0.003 | 0.501 | [−0.017, +0.022] | 3.581 | +0.030 |
| `retrieval_rank`, `blosum_sim` | −0.013 | 0.496 | [−0.044, +0.020] | 3.538 | −0.012 |

**Requirement for 2.0 Å: accuracy 0.638. Best available signal: 0.529.** The gap is not a tuning gap.

### 5.2 The top-25 band, where the zero-information reference is the story

| signal | in-band ρ | 95% CI |
|---|---|---|
| `typicality_med` | +0.093 | [+0.011, +0.173] |
| `dist_z` | +0.091 | [+0.035, +0.145] |
| **`NULL_alpha` — a constant ideal α-helix** | **+0.053** | [−0.035, +0.142] |
| `dist_shipped` | +0.048 | [−0.005, +0.102] |

**A zero-information constant α-helix has in-band skill statistically indistinguishable from the
shipped distance objective's inside its own top-25.** Sprint 16's most-repeated lesson, reproduced
on a new axis at essentially no cost. **ESTABLISHED.**

### 5.3 A non-reproduction that must go on the record

Sprint 8 recorded typicality as "the only signal with positive in-band skill" (tight-pool ρ +0.125,
+0.177 partialling Rg). **On this instrument it does not reproduce.** In the ORACLE-defined band its
ρ is **−0.074 [−0.127, −0.021]** — significantly *negative* — and in the top-75 band it is +0.015
[−0.060, +0.094], null. The **consensus medoid** operator still works (§2.2); the **typicality
score** does not. Those were treated as one mechanism in the record and on this instrument they are
two. The band constructions differ (S8's tight pools vs a score shortlist), which is the likely
reconciliation, but the claim as written does not transfer. **REFUTED as stated; the medoid arm is
unaffected.**

### 5.4 Legacy's eleven components, in-band

Best |ρ| in the top-75 band is `leg_torsion` at +0.044; in the ORACLE band it is `leg_torsion` at
+0.100 [+0.047, +0.154] with a selected RMSD of −0.067 [−0.129, −0.005] against band-random. Every
other component is inside the chance band. Legacy's torsion term — the one the record already
identifies as its non-tautological detector — is also its only term with in-band content, and it is
an order of magnitude short of the requirement.

---

## 5A. The readout triangle — averaging beats both selection readouts everywhere  [Problems B and C]

`s17/sel_readout.py`. Three readouts consuming the **same** shortlist, n = 126, CA-only, m = 1…300
at K = 75 / 500 / 2000, with a matched-random shortlist of the same size at every cell and
**leave-fold-out** choice of m.

The three numbers in circulation came from three different places — argmin 3.454 (shipped), medoid
3.282 (S8), average ≈3.05 (S16 all-atom) — and had never been comparable. They are now.

### K = 500 (ORACLE over the set 1.711, random-of-one 4.453, distance argmin 3.454)

| m | ORACLE@m | medoid | **average** | diversity | rand medoid | rand avg | medoid − avg [95% CI] | W/L |
|---|---|---|---|---|---|---|---|---|
| 5 | 3.025 | 3.491 | 3.231 | 1.319 | 3.974 | 3.730 | +0.260 [+0.195, +0.329] | 34/92 |
| 30 | 2.530 | 3.373 | 3.072 | 1.731 | 3.734 | 3.460 | +0.301 [+0.227, +0.375] | 29/97 |
| **75** | 2.306 | **3.282** | **3.048** | 1.902 | 3.720 | 3.430 | +0.234 [+0.161, +0.305] | 34/92 |
| 200 | 1.995 | 3.427 | 3.090 | 2.255 | 3.711 | 3.408 | +0.337 [+0.262, +0.416] | 22/104 |
| 300 | 1.886 | 3.457 | 3.171 | 2.617 | 3.712 | 3.402 | +0.285 [+0.202, +0.365] | 28/98 |

### The verdict, leave-fold-out in m

| K | argmin | medoid | **average** |
|---|---|---|---|
| 75 | 3.421 | 3.439 | **3.063** |
| 500 | 3.454 | 3.282 | **3.056** |
| 2000 | 3.520 | 3.368 | **3.064** |

**The coordinate average dominates the medoid in all 39 (K, m) cells measured, every one with a CI
excluding zero** (+0.148 to +0.337), and it beats the distance argmin at every K by 0.36–0.46 Å.
It reproduces Sprint 16's raw coordinate average (3.050) to 0.002 Å at K = 500, m = 75 — an
independent cross-instrument replication.

**Three things follow, and they are uncomfortable for the sprint's stated architecture.**

1. **On current numbers, switching from an averaging readout to a selection readout is a
   0.23–0.46 Å regression, not an improvement.** BRIEF §1's central architectural hypothesis is not
   *refuted* — the ceiling argument stands, selection's ceiling is 1.711 Å and averaging's is
   ≈3.05 Å — but the switch cannot be made on today's selectors, and §5 says the in-band
   discriminator that would justify it does not exist among 29 native-free signals.
   **The hypothesis is correct about ceilings and wrong about the immediate move.**

2. **Averaging is flat in K (3.056–3.064) while both selection readouts degrade with K.** That is
   §3.1's mechanism again: averaging consumes the set mean, which improves as K grows; both
   selection readouts consume the set best, which the score's shortlist destroys.

3. **The score's shortlist genuinely contributes to the average**: at K = 500, m = 300 the average
   over the score's top-300 is 3.171 against 3.402 over a matched random 300. So the objective is
   doing real work — as a **filter**, which is what §5 says it is good at — and the failure is
   entirely in the final pick.

**ESTABLISHED.** The honest architectural statement is that the programme currently has a good
coarse filter, a good averaging readout, and no fine ranker; and that the 1.711 Å prize requires
the last of those, which three independent closures in this document say is not available.

---

## 6. THE MIN-OF-N TRAP — the methodological result of this workstream

**Two candidate headline claims died here within an hour, one mine and one the coordinator's.**

Given V selectors and T targets, `mean_t min_v RMSD(t, v)` looks like the ceiling of per-target
calibration or routing. **It is not.** The minimum of V draws falls with V whether or not any
per-target structure exists. The matched reference is *the same statistic over V zero-information
selections from the same candidate set*.

| statistic | Å |
|---|---|
| per-target ORACLE min over the **58 real objective variants**, K = 500 | 2.482 |
| **zero-information: min over 58 random selections from the same set** | **2.148** |
| shipped | 3.454 |

    real min-of-58 vs zero-info min-of-58:  +0.335 [+0.215, +0.463]  fold [+0.151, +0.475]  43W/83L

**The apparent 0.972 Å of "per-target calibration headroom" is not merely explained by the null —
it is worse than the null.** The real variants are correlated enough with one another that their
per-target minimum does not reach the independent-draw floor.

Reproduced on the in-band survey's 25 signals inside the score's own top-25 band:

| | Å |
|---|---|
| ORACLE min over 25 real native-free signals | 2.752 |
| zero-information min over 25 random picks from the band | **2.694** |
| (min over 11 random picks, for comparison with an 11-feature portfolio) | 2.816 |
| band ORACLE best / band random | 2.609 / 3.502 |

    real vs zero-info:  +0.059 [+0.020, +0.100]  fold [+0.028, +0.102]  CI EXCLUDES ZERO

**RULE ADOPTED.** *No "per-target oracle over N selectors" may be quoted as a ceiling without the
min-over-N zero-information reference beside it.* This is BRIEF §5's shape exactly — a quantity
measured correctly and read as a different quantity — and it is now cheap to check
(`sel_cal.minN_null`, `sel_lib.topm_random`).

---

## 7. E4 — target-level calibration  [Problem B]

`s17/sel_cal.py`. Five pre-declared arms on the identical K = 500 set; 27 native-free target-level
features (length, sequence composition, hydrophobic patterning and helical moment, predicted
secondary-structure content, distogram-predicted Rg, model confidence and entropy, retrieval
similarity and dispersion, score mean/sd/skew/best-vs-median gap, pool spread, top-25 agreement).
Every feature is listed in words in `sel_cal.FEATURES` so the leakage claim is checkable without
reading code.

| arm | Å | vs dist |
|---|---|---|
| dist | 3.454 | — |
| cons25 | 3.369 | −0.085 [−0.222, +0.043] |
| **cons75** | **3.282** | **−0.172 [−0.322, −0.031]** |
| cons150 | 3.383 | −0.071 [−0.246, +0.096] |
| blend15 | 3.423 | −0.031 [−0.132, +0.070] |

**The ceiling, with its null.** Per-target ORACLE over the 5 arms **2.917**; min over 5
zero-information selections **3.056**. Here the ceiling *does* clear its null, by **0.138 Å** — a
small but real per-target routing headroom among a small, low-capacity arm set. (Contrast §6: with
58 correlated arms the same statistic is pure artefact. The headroom is real only when V is small
and the arms are genuinely different operators.)

**The three baselines, and the verdict.**

| | Å |
|---|---|
| uncalibrated (dist) | 3.454 |
| **constant** (best arm, leave-fold-out — all 5 folds chose `cons75`) | **3.282** |
| linear ridge λ = 1 | 3.365 |
| linear ridge λ = 10 | 3.334 |
| linear ridge λ = 100 | 3.319 |

    calibrated (lam=100) vs uncalibrated:  -0.135 [-0.276, -0.001]   63W/49L
    calibrated (lam=100) vs CONSTANT:      +0.037 [-0.041, +0.127]   21W/30L

**E4's falsifier fired.** The calibrator beats the uncalibrated arm but **loses to the constant
baseline at every λ**, and the loss shrinks monotonically as λ rises — i.e. the model improves by
becoming the constant. **The boundary trap also fired**: the best λ sits at the top of the ladder,
so the selector is choosing an endpoint (maximal shrinkage toward the constant policy), not a
parameter. This is precisely Sprint 15's failure mode (four native-free estimators at 0.357–0.413
against a 0.516 constant) reproduced on a new feature class.

**No feature carries the per-target margin.** Spearman of each of the 27 features against
(cons75 − dist), the margin whose sign a calibrator must supply, against a ±0.175 null band at
n = 126: largest is `f_pro` at −0.148, then `pred_helix` +0.135, `f_org` +0.114. **Nothing clears
the null band.** The constant-sign baseline is 0.587 (74 targets favour cons75, 47 favour dist).

**REFUTED at this feature class: per-target calibration of a fixed selector combination.**
The 0.138 Å of routing headroom is real and is *not* reachable by these 27 target-level features.

---

## 8. E3 — selection as probabilistic inference  [Problem B]

Partially answered by E2, which contained the decisive arms rather than needing a separate model.
Set-conditional constructions were included in the 58: per-pair standardisation across the
candidate set (`*_zset`), per-pair rank-normalisation across the set (`*_rankset`), and the
centre-blend family, which is exactly "condition the predicted centre on the candidate set".

- `risk_rankset`, the fully nonparametric set-conditional form, is the in-sample best of all 58
  (3.360, −0.094 [−0.226, +0.037] vs shipped) and is chosen by 4 of 5 folds — but its LFO advantage
  is a null.
- `risk_zset` (+0.078) and `z_zset` (+0.088) are *worse* than their unconditioned forms.

**Set conditioning helps only in its rank-normalised form and only by an amount that does not
survive fold-honest selection.** E3's falsifier fired at this feature class. A calibrated listwise
model (Plackett–Luce, conformal) was **not** built, and that is a deliberate stop: §5 shows the
in-band information available to any such model is ρ ≤ 0.086 in the deployable band, and no
likelihood over features with that content can reach the 0.638 requirement. Building it would have
been optimising the report. **OPEN, but priced: the information is not there.**

---

## 9. E5 — the 53 hard targets against the full-universe ceiling  [Problems A and B]

`s17/sel_hard.py`. All ORACLE.

| threshold | top-75 | K = 500 | **FULL universe** |
|---|---|---|---|
| pool best > 2.0 Å | 68/126 | **53/126** | **24/126** |
| pool best > 2.5 Å | 47/126 | 27/126 | **8/126** |

**29 of the 53 hard targets — 55% — were a truncation artefact.** Of the 18 zero-recall FAIL18
targets, only **5** remain above 2.0 Å over the full universe.

**What the residual 24 have in common** (bootstrap CI on the hard−easy difference; * = excludes zero):

| | hard | easy | difference |
|---|---|---|---|
| * chain length | 14.13 | 12.69 | +1.44 [+0.76, +2.15] |
| * predicted helix content | 0.272 | 0.462 | −0.191 [−0.301, −0.072] |
| * hydrophobic fraction | 0.357 | 0.440 | −0.083 [−0.154, −0.012] |
| * polar fraction | 0.260 | 0.186 | +0.074 [+0.016, +0.129] |
| * helical hydrophobic moment | 0.770 | 1.057 | −0.287 [−0.465, −0.105] |
| * mean BLOSUM similarity of the pool | 4.18 | 6.30 | −2.12 [−3.43, −0.89] |
| * universe size | 13,232 | 19,954 | −6,722 [−9,496, −4,025] |
| * pool dispersion / top-25 agreement | 2.36 / 3.08 | 2.09 / 1.96 | +0.26 / +1.13 |

Base rates, so no class statement is read without one: **length 9–11: 0/33 hard (0.0%)**;
12–13: 8/42 (19.0%); 14–16: 16/51 (31.4%). By fold: 8.7%–28.0%, no fold dominates.

**The residual hard set is: longer, less helical, less hydrophobic, more polar, with weaker
sequence matches and a smaller, more dispersed library.** That is a coherent picture — long
non-helical peptides whose library coverage is thin — and it is a **retrieval/library** problem.

**And it is a selection problem of exactly the same size as everywhere else:**

| | ORACLE (full) | realized | GAP |
|---|---|---|---|
| hard set (24) | 2.408 | 4.845 | +2.437 |
| easy set (102) | 1.056 | 3.260 | +2.204 |

The gap is essentially identical. **The hard targets are not a different selection problem; they
are the same selection problem on top of a retrieval deficit.**

**And they are now small.** The 24 residual hard targets contribute only **+0.078 Å** to the
full-universe mean ceiling of 1.313 Å. Capping every one of them at 2.0 Å moves the ceiling to
1.236 Å. **After the full universe, the hard-target problem is no longer where the angstroms are.**
E5's falsifier did not fire; the hypothesis that the hard set was largely a truncation artefact is
supported, and the direction is closed as a *selection* priority.

---

## 10. What I believe the sprint should do with this

Stated as inference from the numbers above, and labelled as such.

1. **Do not widen the candidate set for the deployed pipeline.** Widening is measured to be
   *harmful* to the realized answer (§3) and the ceiling it buys is unreachable. If AUDIT's
   redundancy null shows the wide-K ceilings are partly near-duplicate draws, that reinforces the
   conclusion; it cannot reverse it, because the realized degradation is measured directly.

2. **Do not switch the readout from averaging to selection yet.** The consensus medoid is the best
   *selection* readout — the only arm that beats both its matched control and the distance objective
   at every K ≥ 500, reproducing an S8 result to three decimals — but §5A shows the **coordinate
   average beats it in all 39 (K, m) cells measured** and beats the distance argmin by 0.36–0.46 Å.
   The selection readout's ceiling is genuinely 1.711 Å and averaging's is genuinely ≈3.05 Å, so the
   architectural argument is sound; but the switch is a regression today and stays one until an
   in-band discriminator exists.

3. **Sub-2.0 Å is not reachable by reranking on this instrument.** The requirement is in-band
   ordering accuracy 0.638; the best of 29 native-free signals reaches 0.529 in the deployable band
   and a constant α-helix reaches 0.518. Three independent closures now agree — 58 functionals of
   the distogram (§4), 29 signals in-band (§5), 27 target-level features for calibration (§7) — and
   they are closures by *measurement*, not by failure to find something.

4. **The one direction §7 leaves open is narrow and honest**: 0.138 Å of *real* per-target routing
   headroom exists among a small set of genuinely different operators (argmin vs medoid at several
   widths), it clears its zero-information null, and the 27 target-level features tested do not
   reach it. A different *kind* of per-target observable — not another summary of the same
   pipeline's outputs — would be required. The record already names one such class (a per-target
   measured observable such as chemical shifts) and it is outside this instrument.

---

## 11. Limitations, stated rather than buried

- Every number is on the 126-target **tuning** instrument. The sealed benchmark was not touched.
- ρ_band uses the ORACLE-defined band, which needs the native to draw. It is a **diagnostic**, and
  every deployable claim in §5 uses the native-free top-B bands instead.
- The fold-clustered bootstrap resamples **five** clusters. That is coarse and the intervals are
  correspondingly wide; both intervals are reported everywhere and neither is suppressed.
- `best-of-m random` is an **ORACLE** operation and is used as a currency, never as a competitor.
- The recomputed shipped score matches the instrument to 2.4e−6 (a float32 cache round-trip) and
  the **selection it induces is identical on every target tested** (max deviation exactly 0.0), so
  no variant is being priced against a subtly different baseline. Asserted in `sel_lib.selfcheck`.
- Consensus at K = 500 (3.282 Å) is the Cα of a **selected window**, before any repair stage. It is
  **not** comparable to the deployed 3.204 Å, which includes repair. No pipeline claim is made here.
- The 0.138 Å routing headroom in §7 is measured over five arms chosen by me before the run. A
  different arm set would give a different ceiling; the null-referenced procedure transfers, the
  number does not.
