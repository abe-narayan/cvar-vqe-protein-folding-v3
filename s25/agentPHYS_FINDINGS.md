# SPRINT 25 — LEGACY / AMBER PHYSICS LANE. FINDINGS.

Pre-registered in `s25/PREREG_PHYS.md`, written and on disk **before any configuration ran**.
Point-cloud basis throughout. ORACLE / ACHIEVABLE / PRODUCTION labelled at every appearance.
No benchmark inspection anywhere in this lane. `LOCK_AMBER` taken once, 22 s, announced on open
and on release, and released.

Artefacts, all `complete: true` on the full key set and provenance-stamped with `ST.save_atomic`:

    s25/results/phys_gate.json        126 rows   F0 cache-integrity gate
    s25/results/phys_suite.json      1764 cells  the 7-configuration suite
    s25/results/phys_landscape.json   126 rows   the landscape comparison
    s25/results/phys_form5.json       126 rows   FORM 5 -- the per-target audit, closed
    s25/results/suite_cells/*.json    126 files  per-target cells, resumable
    s25/results/phys_analyse.log                 the full statistics, verbatim
    s25/results/phys_form5.log                   FORM 5, verbatim

---

## 0. THE TWO RESULTS THAT CARRY THE LANE

These come before the headline table on purpose. They are what a reader should remember, and
both rest on a control matched in the operator's own space rather than on a ranking of arms.

> **1. USED AS THE SELECTOR, BOTH GENUINE PHYSICS ENERGIES ARE MEASURABLY WORSE THAN A RANDOM
> 75-CANDIDATE SUBSET OF THE SAME POOL — 5/5 FOLDS EACH. NOT NULL, WORSE.**
> Legacy +0.3302 Å (1.99× MDE), AMBER +0.4554 Å (2.42× MDE), Legacy+AMBER +0.2491 Å (1.53× MDE),
> against a null of 3.4251 Å built from 16 random 75-subsets put through the *identical*
> coordinate average. §1.4.

> **2. DESTROYING A PHYSICS CHANNEL'S CORRESPONDENCE TO CANDIDATES WHILE PRESERVING ITS MARGINAL
> EXACTLY *IMPROVES* THE ENDPOINT EVERYWHERE.**
> Legacy is worse than noise of its own amplitude by +0.327 Å; AMBER by +0.471 Å (26W/100L);
> Legacy inside a distogram-led score costs +0.124 Å against the same amount of noise. And
> **AMBER inside a distogram-led score is statistically indistinguishable from noise** —
> +0.0331 Å at **0.43× MDE**, fold CI spanning zero. That reproduces s24 L16 **with no free
> parameter at all**, which **closes the "the fit was the problem" escape route**. §1.5.

### And in four more lines

3. **The suite is anchored exactly.** Configuration 3 (Distogram) under the classical top-75
   control returns **3.048338 Å**, the incumbent, to `|diff| = 0.00e+00`. The set-equality
   theorem reads back on **1764 of 1764 cells** at `max |RMSD_VQE − RMSD_top-m| = 1.14e−13`.
   The AMBER cache was re-verified against **fresh** single points, bit-identical, max_rel 0.0.
4. **FORM 5 — the last open form of the functional lever — was run and is CLOSED.** A
   native-free per-target audit signal fails at **0.14× its own MDE** and is indistinguishable
   from a permuted-signal null (+0.0049 against +0.0053). **The functional lever is now closed
   in all five of its forms: filter, partition, fitted score, equal-weight score, per-target
   audit — and the fifth by measurement, not by omission.** §5.
5. **Every standing physics fact in the brief is verified on the full 126-target instrument**,
   with two amendments: ρ(AMBER, distogram) = **−0.0186, 95% CI [−0.0583, +0.0211] — the CI
   INCLUDES ZERO, so the correct word is ORTHOGONAL, not anti-correlated**; and ρ(Legacy,
   distogram) has **median +0.556 against mean +0.367**, so the mean is dragged by a minority
   of anti-aligned targets.
6. **My own registered falsifier F4 fired against me, and the reason is worse than the
   prediction.** Raw-moment AMBER did not merely underperform — it is **not monotone in
   float64**: on 40 of 126 targets the moment transform *changes AMBER's own ordering*, with
   exact tie blocks up to 462 of 500, and every Ångström of its apparent advantage lives on
   exactly those 40 targets. The apparent win is the **retrieval order leaking through a tie**.
   §3.

---

## 1. DELIVERABLE 1 — THE 7-CONFIGURATION COMPARISON SUITE

### 1.0 The design, stated so it can be audited rather than trusted

| held identical across all seven | value |
|---|---|
| targets | 126 dev targets, `sorted(pdb)`, same order every configuration |
| folds | the 5 **pinned** folds, read from the instrument, never recomputed |
| pool | shipped **K = 500** BLOSUM retrieval pool; `universe_idx` asserted bit-identical to `s24/cache_amber` on every target |
| basis | **point cloud** — `I.ca_rmsd(I.coordinate_average(W[S]), nat_ca)`; never compared to built-chain |
| metric | full-chain Cα-RMSD |
| **selector** | **genuine CVaR-VQE for all seven** — exact `StatevectorCircuit`, 9 qubits (512 states, 500 real + 12 padding), 3 layers, 80 Adam iters on the **exact parameter-shift gradient**, α = 0.18, T = 0.5, seed 0 |
| readout | **uniform coordinate average over the realised CVaR tail**, in the retained set's own medoid frame |
| pairing | **every contrast paired per target** |

**α = 0.18 was pinned on a native-free criterion and no other.** A 12-target probe on the
distogram channel — no RMSD read — gave median realised tail 75 at α = 0.18 (64 at 0.15, 84 at
0.20). 75 is the production rung `M_PROD`, so the VQE readout and the classical top-75 control
consume the same number of candidates. This is the only way "same readout for every
configuration" and "CVaR-VQE as the matched selector for all seven" can both hold.

**The classical arms are CONTROLS, reported as their own complete matched table for every
configuration, never as a substitute for a VQE row.**

### 1.1 NORMALISATION — the mathematics the brief demanded

For a score `x ∈ R^K` on one target's pool, with `r = rankdata(x)` (average ranks):

    zrank(x) = ( r − mean(r) ) / sd(r)

With no ties `r` is a permutation of `1..K`, so `mean(r) = (K+1)/2` and `sd(r) = sqrt((K²−1)/12)`
**exactly** — constants independent of `x`, of the target and of the functional. `zrank`
therefore maps every channel onto the *identical* empirical distribution: standardised discrete
uniform, mean 0, sd 1, support ±1.7315 at K = 500. **The combined energy of configuration S is**

    E_S = zrank( Σ_{c ∈ S} zrank(x_c) )

* *Equal weight is exact*: unit coefficients on marginals with identical mean, variance, support
  and quantile spacing. Not a hope about matched scales — a statement about matched marginals.
* *Every configuration hands the VQE an identically-scaled Hamiltonian*: without the OUTER
  `zrank`, a |S|-channel sum has `sd = sqrt(|S| + |S|(|S|−1)ρ̄)` = **1.00 / 1.41 / 1.73** at ρ̄ ≈ 0
  for |S| = 1 / 2 / 3, and CVaR trades energy against `T·H`, so **the temperature would silently
  differ between configurations** — a one-, two- and three-channel configuration would each be
  run at a different effective temperature purely because of how many channels it contains.
  **§1.6 shows this is not hypothetical: it is worth up to 0.14 Å of pure normalisation artefact
  masquerading as a physics result.** The outer `zrank` removes it exactly, and being monotone it
  cannot change any classical ordering while doing so.
* `zrank` is **idempotent**, so C1/C2/C3 are unmodified single channels, and **strictly
  monotone**, so it changes no ordering, argmin, level set or CVaR tail membership.

> **Why a bare `E_LEG + E_AMB` would have been meaningless, measured on this pool:** mean
> `|z_moment|` per channel is **AMBER 0.1127, distogram 0.7529, Legacy 0.8013**. Two energies
> with coefficient 1 each are *not* equally weighted — under raw moment AMBER enters at **1/6.7
> of the distogram's amplitude**.

### 1.2 THE GATES

    F0  cache integrity   A pool identity 126/126 PASS   B score reproduction 126/126 PASS
                          C fresh single points, 5 fold-stratified targets x 4 candidates
                            under LOCK_AMBER: BIT-IDENTICAL on all 20, max_rel 0.00e+00
                            AmberSP vs core.amber.refine_coords: 0.00e+00 over 20 comparisons
    F1  anchor            C3 classical top-75 = 3.048338 A  vs incumbent 3.048338  |d| 0.00e+00
    F2  theorem           subset-hood 1764/1764 cells; max |RMSD_VQE - RMSD_top-m| = 1.14e-13
                          equality (reported, never asserted) 1.0000; mean holes 6.59

**The cache is genuine and is used.** 63,000 ff14SB/GBn2 single points, `amber_verify_max_rel`
0.0 on all 126 targets, re-verified against *fresh* single points rather than by re-reading the
stored flag — a gate that cannot fire is not evidence, and this one recomputes.

### 1.3 THE HEADLINE — PRIMARY (rank standardisation), CVaR-VQE arm, n = 126, point cloud

    #  configuration                VQE    median      SD      Q1      Q3   best  worst    m    ent
    3  Distogram                 3.0580   2.8460  1.6510  1.8404  3.8096  0.186  8.093  74.1  8.84
    5  AMBER+Distogram           3.1317   2.9480  1.6616  1.8578  4.1879  0.310  7.305  74.0  8.83
    4  Legacy+Distogram          3.2151   2.9729  1.7078  1.8824  4.2055  0.218  7.537  75.5  8.84
    7  Legacy+AMBER+Distogram    3.2530   3.1343  1.7317  1.7603  4.3613  0.268  7.247  73.8  8.82
    6  Legacy+AMBER              3.6742   3.5787  1.9541  1.7622  4.8930  0.351  8.452  75.1  8.83
    1  Legacy                    3.7553   3.9207  1.9524  1.7935  5.0330  0.489  8.136  74.5  8.83
    2  AMBER                     3.8805   3.5791  1.6710  2.7374  4.7605  0.906  8.930  73.2  8.83
    -- random 75-subset NULL     3.4251                          [best-of-16 draws 3.1109]
    -- ORACLE best-in-pool       1.7108        ORACLE pool mean 4.4533

Against the incumbent (C3 classical top-75, 3.0483 Å), paired, `ST.compare`:

    #  configuration                effect      SE     MDE   x MDE     W/L    fold CI            verdict
    3  Distogram                   +0.0097  0.0081  0.0226   +0.43   60/64   [-0.0005,+0.0250]  NOT MEASURED
    5  AMBER+Distogram             +0.0833  0.0427  0.1196   +0.70   55/71   [+0.0444,+0.1129]  NOT MEASURED
    4  Legacy+Distogram            +0.1668  0.0651  0.1823   +0.91   47/79   [+0.0551,+0.2412]  NOT MEASURED
    7  Legacy+AMBER+Distogram      +0.2047  0.0684  0.1917   +1.07   47/79   [+0.0703,+0.3156]  WORSE (Type-M zone)
    6  Legacy+AMBER                +0.6259  0.1180  0.3305   +1.89   41/85   [+0.4056,+0.7845]  WORSE
    1  Legacy                      +0.7070  0.1180  0.3305   +2.14   38/88   [+0.4999,+0.8373]  WORSE
    2  AMBER                       +0.8322  0.1114  0.3120   +2.67   30/96   [+0.6896,+0.9684]  WORSE

**Nothing beats the incumbent. Four configurations are measurably worse and three are not
measured.** The three "not measured" rows are not near-misses in the promising direction: all
three have positive point estimates and losing W/L records.

### 1.4 THE RESULT THAT MATTERS MOST — BOTH PHYSICS ENERGIES LOSE TO A RANDOM SUBSET

Against the zero-information null **matched in the operator's space** (16 random 75-subsets of
the same pool through the *identical* coordinate average):

    #  configuration          vs random 75-subset       x MDE    W/L     folds  verdict
    2  AMBER                        +0.4554  SE 0.0672   2.42   30/96    5/5   WORSE
    1  Legacy                       +0.3302  SE 0.0593   1.99   41/85    5/5   WORSE
    6  Legacy+AMBER                 +0.2491  SE 0.0582   1.53   45/81    5/5   WORSE
    3  Distogram                    -0.3671  SE 0.0868   1.51   81/45    5/5   BETTER
    5  AMBER+Distogram              -0.2934  SE 0.0653   1.60   84/42    4/5   BETTER
    4  Legacy+Distogram             -0.2100  SE 0.0653   1.15   75/51    5/5   BETTER (Type-M zone)
    7  Legacy+AMBER+Distogram       -0.1721  SE 0.0550   1.12   74/52    5/5   BETTER (Type-M zone)

> **Used as the sole selector over the shipped retrieval pool, a genuine 11-term Legacy
> potential and a genuine ff14SB/GBn2 force field are both MEASURABLY WORSE than picking 75
> candidates at random.** Not null — worse, at 1.99× and 2.42× their own MDEs, with 5/5 folds
> agreeing in sign. This is the cleanest statement the suite produces and it is beyond reproach:
> same pool, same readout, same selector, same 126 targets, and a null matched in the operator's
> own space rather than a degenerate one.

### 1.5 THE RANK-PERMUTED CONTROL — THE PHYSICS CHANNELS ARE WORSE THAN NOISE

The control preserves each physics channel's **marginal exactly** and destroys only its
**correspondence to candidates**. Any configuration whose real arm cannot beat its permuted twin
is carrying no usable information.

Both arms are the **classical top-75** on the identical readout, which is where the permuted
control is exact and free.

    #  configuration          real(top75) permuted   effect   x MDE    W/L     folds  verdict
    1  Legacy                 3.7543     3.4269    +0.3274    1.87   44/82     5/5   WORSE
    2  AMBER                  3.8810     3.4099    +0.4712    2.34   26/100    5/5   WORSE
    6  Legacy+AMBER           3.6794     3.4264    +0.2530    1.52   41/85     5/5   WORSE
    4  Legacy+Distogram       3.2134     3.0898    +0.1236    1.00   42/84     5/5   WORSE (Type-M zone)
    7  Legacy+AMBER+Dist      3.2497     3.1201    +0.1297    1.19   46/80     5/5   WORSE (Type-M zone)
    5  AMBER+Distogram        3.1293     3.0962    +0.0331    0.43   59/67     3/5   NOT MEASURED

Read this carefully, because the direction is the whole finding. Permuting a channel does not
delete it — it converts it into noise of the **same amplitude**. So `real − permuted > 0` says
the channel is **worse than noise of its own size**, not merely uninformative.

* **Legacy and AMBER standing alone are worse than their own noise**, and their permuted arms
  land on the random null (3.4269 and 3.4099 against 3.4251) — exactly where a destroyed
  ordering should land, which is itself a soundness check on the control.
* **Inside a distogram-led score, Legacy at equal weight costs +0.124 Å against the same
  amount of noise**, and Legacy+AMBER together cost +0.130 Å.
* **AMBER inside a distogram-led score is indistinguishable from noise**: +0.0331 at **0.43×
  MDE**, fold CI [−0.029, +0.087] spanning zero, 59W/67L. This is s24 L16's central result —
  "whatever movement the real arm shows is not distinguishable from a control with all of
  AMBER's information removed" — reproduced under an entirely different parameterisation
  (**equal weight, not a fitted `w`**), on the same 126 targets, through a quantum selector
  rather than a classical sort.

### 1.6 THE m-LADDER — REPORTED, NOT HIDDEN, AND UNDER RANK NORMALISATION IT IS NEGLIGIBLE

By the s24 set-equality theorem the realised CVaR tail is the classical top-m of the
configuration's own energy order, so a configuration-to-configuration comparison could be
confounded by `m` rather than by the energy (`operator-consumes-set-mean`). The exact
decomposition, persisted per cell:

    RMSD_VQE(c) = RMSD_top75(c) + [RMSD_top-m_c(c) − RMSD_top75(c)] + eps_c

    normalisation   realised m across the seven      m-ladder term        eps
    rank            73.2 - 75.5                      +0.0011 .. +0.0097   <= 1.1e-13
    moment          35.9 - 92.8                      -0.100  .. +0.141    <= 1.1e-13

> **Under the primary normalisation the m-ladder spans 0.009 Å across all seven configurations,
> so the headline table is a comparison of ENERGIES and not of set sizes.** `eps` is at machine
> precision on every one of the 1764 cells — the theorem reading back, never a quantum effect.

**And this is the measured justification for the outer re-standardisation.** Under raw moment
the same seven functionals drive the realised tail from **35.9 to 92.8 candidates** and the
state entropy from 8.11 to 9.00 bits out of a 9-bit maximum. For Legacy the two normalisations
give a **bit-identical classical top-75** (set overlap 1.0000) and yet move the VQE's tail from
74.5 to 35.9 candidates, worth **+0.0975 Å through the m-ladder alone**. Normalisation changes
the quantum selector even where it provably cannot change the classical one.

### 1.7 THE CVaR TAIL MATERIALLY PARTICIPATES (Pillar 1)

    rank arm, all seven configurations:  entropy 8.82-8.84 bits of a 9-bit maximum
                                         ESS 59.4-61.9 of 512 basis states
                                         realised m 73.2-75.5 against a NOMINAL alpha*2^n = 93

The trained state is broad — it has not collapsed to an argmin, which is the failure mode
`core/pipeline` documents at α = 1 — and it reaches α mass in **74 states rather than the 93 a
uniform distribution would need**, i.e. it puts ~25% more mass on the low-energy prefix than
uniform. The tail is doing the selecting, and it is doing it on a trained distribution rather
than on a collapsed one.

**What the circuit contributes to the ANSWER is nonetheless nil, and that is the theorem, not a
result.** `VQE vs classical top-m` is `0.0000 ± 0.0000` on all seven configurations; `VQE vs
classical top-75` is +0.0011 to +0.0097, every one below its own MDE. Any future reader who sees
the quantum arm matching the classical control exactly should read s24 `agentD_FINDINGS.md` §1,
not celebrate.

### 1.8 FOLD-WISE AND LENGTH-STRATIFIED BREAKDOWN (rank, VQE arm)

    #  configuration              f0      f1      f2      f3      f4   |  n<=11   n<=13   n<=15    n>15
    1  Legacy                 3.6507  3.7397  3.2858  4.1283  3.9599  | 3.6461  3.2077  4.2618  4.2988
    2  AMBER                  3.8593  3.7258  3.5819  4.1549  4.0553  | 3.3744  3.4661  4.6026  4.4730
    3  Distogram              2.8095  3.0212  2.9872  3.3194  3.1519  | 2.7130  2.8450  3.7815  2.9926
    4  Legacy+Distogram       3.0436  3.2355  2.9273  3.4816  3.3780  | 3.0068  2.9699  3.7275  3.3010
    5  AMBER+Distogram        2.8648  3.0322  3.0576  3.4017  3.2851  | 2.6577  2.8461  3.8534  3.4166
    6  Legacy+AMBER           3.5763  3.5618  3.2088  4.1113  3.8948  | 3.4988  3.1481  4.2209  4.2210
    7  Legacy+AMBER+Dist      3.0734  3.1219  2.9617  3.6794  3.4190  | 2.9109  2.9566  3.8452  3.5375
    band sizes                                                          33      42      30      21

The ordering of the seven is stable across folds; **fold 2 is the one place the physics
configurations are least penalised** (C4 at 2.9273 beats C3 at 2.9872 there — the single
per-fold sign reversal in the C3-vs-C4 contrast, and the only fold on which any Legacy- or
AMBER-containing configuration beats the distogram alone). AMBER's penalty is
**strongly length-dependent**: at
n ≤ 11 it costs 0.66 Å against C3 and at 13 < n ≤ 15 it costs 0.82 Å, consistent with the
all-atom clash census growing with chain length on unrelaxed retrieval windows.

### 1.9 "THE BEST OF THE SEVEN" IS AN ORDER STATISTIC — AND THE USUAL NULL DOES NOT APPLY

    best configuration by dev mean                        C3, 3.0580 A
    mean pairwise correlation between the seven arms       0.8627   ->  k_eff = 1.13
    ORACLE per-target min over the seven                   2.7231 A   (-0.3349 vs C3)
    pooled i.i.d. best-of-7 null                           gain -1.7953, share_accounted 5.360

**`share_accounted = 5.36` must not be read as "the arm is the null". THE GENERAL FORM, which is
worth carrying beyond this lane: a best-of-K null is INAPPLICABLE when the K arms are not
independent, and `k_eff` is how you know.** It is a signal that the null itself is inapplicable: `ST.best_of_k_null` resamples seven values i.i.d. from the pooled
distribution, and these seven arms have mean pairwise correlation 0.86 — they are effectively
**1.13 independent configurations, not 7**. The honest statement is that **−0.3349 Å is an
ORACLE ceiling on any per-target configuration-switching rule, with full leakage**, that a large
share of it is between-arm noise on nearly-identical arms, and that it is **not a selection
rule**. §3 pre-registers the only test that could convert it.

---

## 2. DELIVERABLE 2 — THE LANDSCAPE COMPARISON

Re-measured on the shipped K = 500 pool where the pool or the score can have moved the number;
**cited, not re-run**, where the quantity is a property of the energy function on the continuous
torsion manifold that no candidate pool can change.

| axis | Legacy | AMBER | source | physical reading |
|---|---|---|---|---|
| **energy scale** | median −10.9, IQR 9.7, **2.90 decades** | median 2.0e7, IQR 8.5e8, **15.26 decades**, worst point **1.87e29** | **[S25]** | Legacy is a smooth soft potential; AMBER on unrelaxed windows is a hard-core potential evaluated inside its own singularity. |
| **candidate-energy distribution** | 6.5% inside \|z_raw\|<0.1; top-10 variance share **0.183** | **97.0%** inside \|z_raw\|<0.1; top-10 variance share **0.9966** | **[S25]** | Ten of five hundred candidates carry 99.7% of AMBER's variance. Its mean and sd describe one clash, not the pool. |
| **ρ(E, distogram)** | **+0.3671, 95% CI [+0.285, +0.450]** (median +0.556, 26/126 negative) | **−0.0186, 95% CI [−0.0583, +0.0211] — CI INCLUDES ZERO, so ORTHOGONAL, not anti-correlated** (median −0.028, 70/126 negative) | **[S25]**, s24 +0.3875 / −0.0270 | Legacy is *partly the same instrument* as the deployed scorer. **AMBER's interval spans zero: the defensible word is ORTHOGONAL and a bare point estimate of −0.019 or s24's −0.027 invites over-reading as anti-correlation.** Orthogonality is the mechanistic reason AMBER is the one functional that buys a non-parallel direction — and §1.5 is why it still cannot be used. |
| **ρ(Legacy, AMBER)** | **−0.0896, 95% CI [−0.1254, −0.0537]**, median −0.0955, 79/126 negative | | **[S25]**; s24 retrieval −0.0829; s20 torsion −0.0886 | **Three instruments, three manifolds, the same number.** The two genuine Hamiltonians are mildly ANTI-correlated in rank — they are not two views of one potential. |
| **ρ(E, RMSD) whole pool** | +0.3073 ± 0.034 | **−0.0266 ± 0.021** | **[S25]** ORACLE DIAGNOSTIC | The whole-pool column is the documented trap: the **Rg-only control is +0.2794**, so essentially *all* of Legacy's apparent ranking skill here is radius of gyration. |
| **ρ(E, RMSD) partialling Rg** | +0.2250 ± 0.032 | **+0.0473 ± 0.017** | **[S25]** | With compactness removed Legacy retains real but modest skill; AMBER's sign flips from −0.027 to **+0.047** — its raw anti-correlation was an Rg artefact in the opposite direction. |
| **ρ(E, RMSD) in-band (<3 Å)** | +0.1086 ± 0.039 | **+0.0913 ± 0.028** | **[S25]** ORACLE DIAGNOSTIC | In-band, where compactness carries almost nothing (**Rg control +0.0213**), both energies are genuinely, weakly positive — `physics-ranks-real-geometry` reproduced, but at ⅓ of its quoted size and **beaten by the distogram's +0.2361**. |
| **compactness response** | top-75 is **−0.758 ± 0.023 Å** in Rg (MORE compact); ρ(E, Rg) **+0.6005** | top-75 is **+1.103 ± 0.042 Å** in Rg (MORE expanded); ρ(E, Rg) **−0.2697** | **[S25]**; s20 −0.4477 / prefers expanded | **The two energies point in opposite directions on the single axis that dominates this pool.** Legacy's "steric" term behaves as a compactness/typicality prior; AMBER's nonbonded core rejects compaction because unrelaxed windows cannot pack without clashing. |
| **steric sensitivity** | top-75 min-\|i−j\|≥3 CA–CA **−0.214 Å** vs pool | top-75 **+0.646 Å** vs pool; 58.6% of every pool above 1e4 kcal/mol | **[S25]**; s20 all-atom census | AMBER is measurably a **feasibility filter**: what it actually selects on is separation, not fold. Legacy tolerates tighter contacts because its steric term is soft. |
| gradient direction (joint) | | **cos(∇Legacy, ∇AMBER) = −0.3512 [−0.472, −0.222]**, 27/30 targets negative | **s20**, cited | The anti-correlation is not a rank artefact — the two forces genuinely push apart in torsion space. |
| gradient magnitude | moves 5.6× smaller | first steps are clash relief from `E_0 = +1.9e8` | **s20**, cited | AMBER's gradient is dominated by escaping its own singularity, which is why no minimisation is admissible here. |
| Hessian conditioning | baseline | **three orders larger**; still **473×** worse after log conditioning | **s20**, cited | Conditioning fixes the scale, not the curvature. |
| anisotropy / participation | 5.8 / 0.4221 | **18.8** (0W/30L) / **0.0745** (30W/0L) | **s20**, cited | ~7% of AMBER's modes carry all the curvature — a steric singularity, not a folding landscape. |
| ruggedness / negative curvature | baseline | **2×** local minima per 2π; **slightly LESS** negative curvature, −0.0346 [−0.071, +0.012] | **s20**, cited | s20's own prediction refuted: AMBER is rougher, not more saddle-ridden. |
| **locality** | full-register | full-register | `torsion-space-locality-theorem` | **"AMBER is less local" is a CATEGORY ERROR.** `d_ij` depends on exactly the `j−i−1` residues between i and j; all-atom supports are one residue wider; the union of supports is the whole chain for **both**. |
| **surrogate feasibility** | n/a | bonded subset tracks at **ρ 0.211** against a 0.7 bar | s16/s20, cited | **No defensible AMBER surrogate exists.** AMBER is used at full cost or not at all. |
| **trainability / CVaR concentration** | rank: 8.83 bits, ESS 60.1, m 74.5 | rank: 8.83 bits, ESS 59.4, m 73.2 | **[S25]** | Under rank standardisation all seven Hamiltonians train identically well — the differences in the endpoint are the *energies*, not the optimiser. |
| **trainability, raw moment** | 8.11 bits, ESS 24.2, m 35.9 | **9.00 bits, ESS 92.2, m 92.8** | **[S25]** | Raw-moment AMBER is a **9.00-of-9-bit uniform state**: 97% of the energies are numerically equal, so there is nothing for the circuit to train on. Un-normalised AMBER is not a hard Hamiltonian — it is an *empty* one. |
| **structural usefulness** | 3.7553 Å, **worse than random (3.4251)** | 3.8805 Å, **worse than random**, worse than its own noise | **[S25]** §1.4–1.5 | Both energies are complementary about GEOMETRY (opposite compactness signs, ρ = −0.09) and jointly useless about ACCURACY on this manifold. |

**Two standing corrections carried forward so they are not re-derived wrongly.** (i) "AMBER is
less local" is a category error — both energies are full-register in torsion space. (ii) "Legacy
is the less physical model" inverts the measurement: **Legacy is a compactness/typicality model
wearing a physics vocabulary** (ρ(E, Rg) = +0.60, and its whole-pool RMSD correlation of +0.307
sits on an Rg control of +0.279), which is the opposite of the naive prior given that its steric
term carries the largest weight.

---

## 3. THE NORMALISATION FORK — F4 FIRED AGAINST ME, AND THE MECHANISM IS WORSE THAN THE PREDICTION

I registered: under raw moment, configurations 2/5/6/7 move materially and get **worse**.
Measured, paired, VQE arm:

    #  configuration              moment    rank     d       x MDE   W/L      verdict
    2  AMBER                      3.7539  3.8805  -0.1266   -2.05   92/34    BETTER
    5  AMBER+Distogram            3.0869  3.1317  -0.0448   -0.31   64/62    NOT MEASURED
    3  Distogram                  3.0916  3.0580  +0.0336   +0.64   54/72    NOT MEASURED
    1  Legacy                     3.8517  3.7553  +0.0964   +0.26   55/68    NOT MEASURED
    7  Legacy+AMBER+Distogram     3.3466  3.2530  +0.0936   +0.41   53/73    NOT MEASURED
    4  Legacy+Distogram           3.3770  3.2151  +0.1619   +0.81   51/75    NOT MEASURED
    6  Legacy+AMBER               3.8969  3.6742  +0.2227   +0.54   57/69    NOT MEASURED

**C2 is BETTER under raw moment at 2.05× its own MDE. My falsifier fired.** So I chased it, and
the mechanism disqualifies the arm rather than the choice:

### 3.1 Raw-moment AMBER is not monotone in float64

`zmoment(x) = (x − mean)/sd` is monotone in exact arithmetic. It is **not monotone in double
precision** when the sd is set by a single 1e28 outlier: every candidate below ~1e12 maps to
**bitwise-identical** z, because their differences fall under `2.2e−16 × sd`.

    distinct values after zrank    :  mean 464.8 of 500   min 422
    distinct values after zmoment  :  mean 439.8 of 500   min  26
    largest EXACT tie block, zmoment: mean 18.4           max 462 of 500
    zmoment preserves raw AMBER's own argsort on           86 / 126 targets

### 3.2 Every Ångström of the apparent win lives on the 40 targets where the order was destroyed

    C2, classical top-75:
      86 targets where the argsort IS preserved :  max |moment − rank| = 0.00e+00   (identical)
      40 targets where the argsort is BROKEN    :  moment 3.8556   rank 3.9396   d −0.0840
      whole instrument                          :  d −0.0267, ALL of it on those 40

The VQE arm's larger −0.1266 Å decomposes exactly, and neither half is AMBER doing anything:

      ordering broken by float64 (classical top-75)        −0.0267
      m-ladder (moment tail 92.8 vs rank tail 73.2)        −0.1004
      sum                                                  −0.1271   vs measured −0.1266

The second term is a **larger** tail on a below-random objective, which moves the answer toward
the pool mean — i.e. toward the random null at 3.4251 — and therefore helps precisely because
C2 is worse than random. Neither channel of the "win" is information.

`np.argsort` on an exact tie block returns **array order**, which on this pool is the **BLOSUM
retrieval order**, and the retrieval order is *not* neutral: `ρ(pool index, ORACLE candidate
RMSD) = +0.0540 ± 0.0077`, seven standard errors from zero. So moment-C2's "advantage" is the
retrieval ranking leaking into a numerically-degenerate score. This is project memory's
documented trap — *tie-breaking leaks the pool order* — arriving through a normalisation choice
rather than through an `argmin`.

### 3.3 And where raw moment appears to help a combination, it helps by deleting the channel

    top-75 set overlap with C3 (Distogram alone), mean over 126 targets:
       rank   AMBER+Distogram vs C3 :  0.5266     <- AMBER genuinely changes half the set
       moment AMBER+Distogram vs C3 :  0.9930     <- AMBER changes nothing
       rank   Legacy+Distogram vs C3 :  0.5208
       moment Legacy+Distogram vs C3 :  0.4344
       moment vs rank, Legacy alone   :  1.0000    <- Legacy has no scale pathology (2.9 decades)

> **Under raw moment, "AMBER + Distogram" IS the distogram** — 99.3% the same 75 candidates. The
> moment arm is not a competing normalisation of AMBER; it is a normalisation that **annihilates**
> AMBER, and it scores well for exactly that reason. The rank choice is vindicated by a stronger
> argument than the one I pre-registered: raw moment is not merely uninformative here, it is
> **not order-preserving**, and an arm that is not order-preserving cannot be compared to one
> that is.

**Recorded as a fork I got wrong.** My pre-registration predicted the raw-moment arm would be
uniformly worse. It is better on one configuration, and I would have had to withdraw the
normalisation choice had I stopped at the table. What saved it was running the mechanism check
rather than the outcome — which is exactly the discipline s25 L2 used on the calibration arm.

---

## 4. WHAT THE SUITE SETTLES, AND WHAT MY REGISTERED PRIOR DID

**F3, my registered prior, is CONFIRMED in all four of its parts** and I state it plainly
because a confirmed pessimistic prior is weak evidence and should be labelled as such:

    (i)   C3 is the best single configuration                              CONFIRMED
    (ii)  no AMBER-containing configuration beats its AMBER-free twin
          PAST ITS OWN MDE                                                 CONFIRMED
            C5 vs C3  +0.0737 (0.58x MDE)   C7 vs C4  +0.0379 (0.39x MDE)  -- AMBER costs
            C6 vs C1  -0.0811 (0.91x MDE, TYPE-M ZONE)                     -- AMBER appears
              to help Legacy by 0.08 A, but at 0.91x its own MDE that is not a result and I
              do not claim it. All three are NOT MEASURED.
    (iii) C4 (Legacy+Distogram) is worse than C3                           SIGN CONFIRMED,
            +0.1571, 0.85x MDE, 79W/47L against, 4/5 folds -- NOT MEASURED, and I do not claim it
    (iv)  C1/C2/C6 are 0.5-1.5 A worse than C3                             CONFIRMED (0.63-0.83)

**What is genuinely new here rather than a reproduction:**

1. **The physics energies are worse than a matched random subset, measurably.** s24 established
   that AMBER's ordering buys nothing; this establishes that it *costs* — 0.46 Å against a null
   matched in the operator's own space, at 2.42× MDE, 30W/96L, 5/5 folds.
2. **The equal-weight parameterisation reaches the same place as the fitted one.** s24 L16
   closed the score form with an oracle `w` worth 0.0148 Å. This closes it again with **no free
   parameter at all**, which removes the "the fit was the problem" escape route.
3. **The rank-permuted control now fires in the harmful direction for Legacy** (+0.124 Å inside
   C4, 1.00× MDE) as well as being null for AMBER (+0.033, 0.43× MDE). s24 only had the AMBER
   half. Legacy is worse than noise of its own amplitude; AMBER is exactly noise.
4. **Normalisation changes the quantum selector even where it provably cannot change the
   classical one** — bit-identical top-75 for Legacy, tail 74.5 → 35.9, worth +0.0975 Å.
5. **The mechanism behind "AMBER selects a non-parallel direction"** (s24 L10-A, bias cos
   0.5693) is now named on this pool: it is **compactness with the opposite sign**. AMBER's
   top-75 is +1.10 Å more expanded than the pool; Legacy's is −0.76 Å more compact. They
   disagree about geometry (ρ = −0.09) and are jointly blind to accuracy.

---

## 5. THIRD DELIVERABLE — FORM 5, RUN AND CLOSED. THE FUNCTIONAL LEVER IS CLOSED IN ALL FIVE FORMS.

Pre-registered in full in `s25/PREREG_PHYS.md` **ADDENDUM A1–A7**, written to disk before the
run, with a loudly-stated expectation of failure. Authorised by the coordinator. Architecture
FROZEN throughout: nothing here was permitted to change what ships unless it cleared its
falsifier decisively. `s25/phys_form5.py` → `s25/results/phys_form5.json` (`complete: true`,
provenance-stamped), full log at `s25/results/phys_form5.log`.

### 5.1 The form, and why it is not one of the four already closed

The functional lever had failed as a **filter** (s24 D1-C, 5/5 harmful), a **partition** (s24
D1, every partition at its own floor), a **fitted score** (s24 L16, ORACLE ceiling 0.0148 Å) and
— §1.3–1.5 above — an **equal-weight score**. Form 5 is the one remaining class: the physics
never chooses a candidate, never splits the pool and never enters the ranking. It supplies a
**per-target SIGN at inference**, which is the only leverage `in-band-ordering-is-per-target`
leaves open.

    s_t = mean_{i in T_t} [ zrank(E_LEG)_i - zrank(E_AMB)_i ]     T_t = C3's own top-75, NATIVE-FREE

*Mechanism, stated before the run.* s25 L1 measured the distogram predicting distances
systematically **too long**, growing to −0.589 Å at long separation — its posterior favours
over-expanded structures. §2 measures the two physics channels as an **opposite-signed
compactness pair on this exact pool** (Legacy's top-75 −0.758 Å in Rg, AMBER's +1.103 Å; ρ(LEG,
Rg) = +0.60 against ρ(AMB, Rg) = −0.27). So `s_t` asks whether the distogram's own answer sits
where the physics says the pool's feasible region is — without either energy ranking anything.
The arm switches between C3 and C5; **threshold AND direction fitted on the other four folds
only**, with the two degenerate ends (always-C3, always-C5) on the grid so the fit is free to
decline to switch.

### 5.2 THE FALSIFIER DID NOT CLEAR. FORM 5 IS CLOSED.

    PRIMARY   switched arm (nested LOFO)  vs  C3
              3.0629 vs 3.0580     effect +0.0049   SE 0.0130   MDE 0.0363   **0.14x MDE**
              iid CI [-0.0203,+0.0290]   fold CI [-0.0113,+0.0242]   folds same sign 3/5
              14W / 15L / **97 TIES**    median +0.0000    VERDICT: NOT MEASURED

    N1 permuted-signal null (marginal kept, target correspondence destroyed)
                                   effect +0.0053   **0.11x MDE**   fold CI [-0.0140,+0.0204]

> **The real signal and a signal with all of its target-correspondence destroyed produce
> indistinguishable results: +0.0049 against +0.0053.** This is s24 L16's rank-permuted control
> firing again, in the one place the lever had not yet been tested, and it is the reason the
> closure is by absence of signal rather than by lack of power.

Three further diagnostics, each of which independently says the same thing:

* **The fitted direction does not survive resampling.** The nested fit chose `s > τ → C5` on
  folds 0, 2, 4 and `s < τ → C5` on folds 1, 3. A rule whose *sign* flips between training sets
  is fitting noise, and per-fold reporting is what surfaced it — the mean would have hidden it,
  exactly as in s24 L16.
* **97 of 126 targets are ties.** The fit switched only 3–9 held-out targets per fold: with the
  degenerate ends available, it very nearly chose "do nothing", which is the honest behaviour of
  a grid handed a signal with no content.
* **M1, ORACLE DIAGNOSTIC** (reads the native; can explain a result, can never produce one):
  `ρ(s_t, C5−C3) = +0.0506`. The signal carries essentially no information about which arm wins.
  C5 beats C3 on 58 of 126 targets, by −0.2565 Å when it does and losing +0.3553 Å when it does
  not — a real per-target spread that no native-free statistic tested here can point at.

### 5.3 The declared secondary, with its multiplicity honoured

    SECONDARY  r_t = Rg contrast of C3's top-75 vs the pool, IDENTICAL nested CV
               3.0568 vs 3.0580   effect -0.0012   MDE 0.0364   **0.03x MDE**
               fold CI [-0.0194,+0.0195]   folds same sign 2/5   9W/10L/107T   NOT MEASURED
               N1 permuted null +0.0524 (0.63x MDE) -- LARGER than the real arm

Registered as a second signal with a Bonferroni bar at α/2 **before** it ran, so it could not be
promoted after the fact. It does not clear the plain bar, let alone the corrected one. Its
permuted null being *larger* than its real arm is diagnostic in its own right: at this effect
size the fitting procedure's own noise dominates whatever the signal contributes.

### 5.4 THE ORACLE CEILING, AND WHY IT CANNOT BE PROMOTED

    per-target min(C3, C5)      2.9400 A     **-0.1181 A vs C3**, FULL LEAKAGE, ORACLE
    corr(C3, C5)                0.9536   ->  **k_eff = 1.02 independent arms, not 2**
    nested CV recovers          +0.0049 of -0.1181  -- it recovers none, and moves the wrong way

**A best-of-K null is INAPPLICABLE when the K arms are not independent, and `k_eff` is how you
know.** At corr 0.9536 these are effectively one arm sampled twice, so the −0.1181 Å is
overwhelmingly the noise between two nearly-identical draws rather than a switchable difference.
This is the same lesson as the suite's `share_accounted = 5.36` (§1.9) in its second instance,
and it generalises: *a null that returns an impossible number, or an arm count that collapses to
`k_eff ≈ 1`, is telling you about the null, not about the arm.*

### 5.5 Disposition

> **Form 5 fails at 0.14× its own MDE against a permuted-signal null it cannot be distinguished
> from. THE FUNCTIONAL LEVER IS NOW CLOSED IN ALL FIVE OF ITS FORMS — filter, partition, fitted
> score, equal-weight score, and per-target audit — and the fifth was closed by measurement, not
> by omission.**

The architecture is unchanged and nothing is brought back to the coordinator for the freeze.
**My pre-stated expectation was confirmed, which is weak evidence and is labelled as such**; the
load-bearing part is N1, not my prior.

---

## 6. WHAT DAMAGED MY OWN EXPECTATIONS

1. **My normalisation falsifier F4 fired against me** (§3). I predicted raw moment would be
   uniformly worse; it is better on AMBER alone at 2.05× MDE. Had I reported the table without
   the mechanism check I would have had to withdraw the primary normalisation. The rescue was
   not foresight — it was following the s25 L2 discipline of checking the mechanism, and the
   mechanism turned out to be a **float64 monotonicity failure**, which I did not anticipate
   despite `pauli-spectrum-delta-spike-artefact` recording the scale pathology that causes it.
   Project memory told me the moment z-score was *uninformative*; it did not tell me it was
   *not order-preserving*, and neither did I until I measured it.
2. **A best-of-K null is inapplicable when the K arms are correlated, and I nearly missed it
   twice.** `share_accounted = 5.36` on the best-of-seven null is my own harness misapplied by
   me. `ST.best_of_k_null` assumes i.i.d. draws; seven arms at mean pairwise correlation 0.86 are
   not that. I nearly reported "the best-of-7 ceiling is 536% accounted for by its own null",
   which would have been a confident statement produced by a null that does not apply. The
   correct reading is that the null is inapplicable and the ceiling stands as ORACLE. **A null
   that returns an impossible number is telling you about itself, not about the arm.** The same
   trap re-appeared in Form 5 at `corr(C3, C5) = 0.9536`, `k_eff = 1.02` (§5.4), which is why the
   general rule — **`k_eff` is how you know the null does not apply** — is now stated rather than
   rediscovered. This sits beside the audit lane's finding that an own-row null returns 83% on
   pure noise.
3. **My registered prior was right, and that is worth almost nothing.** F3 was confirmed in all
   four parts. A pessimistic prior that is confirmed is weak evidence; the load-bearing results
   here are the two *controls* — the operator-matched random null and the rank-permuted twin —
   not the fact that I guessed the ordering.
4. **α = 0.18 is a choice I made, on a native-free criterion, and it is a fork.** Median
   realised tail 75 is the production rung, and I pinned α there so the VQE and the classical
   control consume the same number of candidates. That is defensible and pre-registered, but it
   is *my* judgement, not a measurement: at α = 0.15 the tail is 64 and at 0.25 it is 104, and
   `operator-consumes-set-mean` says m matters. **Reviewed and approved as made; the suite is not
   re-run at another α**, and §1.0's paragraph is the standing answer to "why is α not the
   deployed value" — because the deployed value is tuned for a 128-candidate consensus-medoid
   readout, and this suite needs the realised tail to land on the production rung of 75 so that
   "same readout for all seven" and "CVaR-VQE for all seven" hold simultaneously.
5. **Form 5's failure is the weakest kind of confirmation and I say so.** I registered a loud
   expectation that it would fail and it failed. What makes §5 a result is not my prior but N1:
   the real signal and a signal stripped of all its target-correspondence give +0.0049 against
   +0.0053. Had I reported only the primary contrast I would have published an underpowered
   non-result as a closure; the permuted null is what makes it an absence of signal instead.

---

## 7. WHAT REMAINS LIMITING

Not a manufactured win, per BRIEF §6. **The closest defensible result in this lane is the
incumbent**: no configuration beats C3, and C3's VQE arm reproduces the incumbent to within
0.0097 Å (0.43× MDE). I propose no promotion.

What the suite establishes is *where the constraint is not*, with a sharper instrument than
Sprint 24 had:

* **Not in the choice of physical functional.** Two genuine, mutually anti-correlated
  Hamiltonians (ρ = −0.0896, reproduced on a third manifold) are both **worse than a random
  subset** of the same pool, and both worse than their own rank-permuted noise.
* **Not in any remaining USE of the functional.** With Form 5 measured and closed, the lever is
  closed in **all five** of its forms — filter, partition, fitted score, equal-weight score, and
  per-target audit at inference. There is no sixth form on the table and I do not propose one.
* **Not in how the two are combined.** Equal weight on exactly-matched marginals reaches the
  same null as s24's fitted `w`, and removes the "the fit was the problem" escape route.
* **Not in the selector.** All seven Hamiltonians train identically well under rank
  standardisation (8.82–8.84 bits, ESS 59–62), the realised tail spans 73.2–75.5 across all
  seven, and `VQE − classical top-m` is zero to machine precision on 1764 of 1764 cells.
* **Not in the normalisation** — provided it is rank. Raw moment is not a competing choice; it
  is numerically degenerate on AMBER and its apparent advantage is the retrieval order leaking
  through an exact tie.

**The architectural sentence for the final report, from this lane's evidence:**

> The functional lever is closed in all five of its available forms — filter, partition, fitted
> score, equal-weight score, and per-target audit at inference — the last of them by measurement.
> CVaR-VQE is a genuine, exactly simulated selector whose realised tail is provably the classical
> top-m of whatever Hamiltonian it is given, and which trains identically well on all seven.
> Legacy is a compactness/typicality prior wearing a physics vocabulary — ρ(E, Rg) = +0.60, and
> its apparent ranking skill sits almost entirely on an Rg control. AMBER is a genuine ff14SB/GBn2
> energy measurement and a feasibility diagnostic that is orthogonal to the deployed scorer
> (ρ = −0.019, CI includes zero) and selects for *expansion*, not for fold. On the shipped
> retrieval manifold both are measurably worse than choosing at random, and the distogram alone
> is the best of the seven configurations. The physics belongs in the system as a measurement and
> a feasibility check, which is what it is, and not as a selector, which it is not.
