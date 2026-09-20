# AUDIT V — adversarial audit of `s30/REPORT_S30.md` §1, §3, §4, Appendix A

Lane V, report adversary, contract rule 28. 2026-09-20.
Read-only audit. No project file was modified; this file is the only output.
Report state audited: 431 lines, including the 13:55 Appendix A edit (A.3 reversal + new §A.7).
Verifier state audited: `s30/s30_verify.py` at 22 checks.

Everything below was recomputed from artefacts unless marked otherwise. Two independent
recomputations were run from raw data (`s12/instrument.py` + `s8/generate_univ/*.npz`), not from
the lane JSONs, as requested.

---

## VERDICT (read this first)

**Not safe to publish as it stands.** Four defects are severe enough to change what the report
says, and two of them invert a headline:

- **§4.2's "5.2× / worth five times their face value" is a free parameter.** It is
  `I = −(d/2)log2(1−ρ²)` at `d = 3n−6 = 32.88`. At `d = 6` — the effective dimension the sprint's
  own S30-L15 established and which the lane's own ledger flags as the correction ("S30-L15
  corrects the `d` this should be evaluated at, and the corrected numbers are smaller") — the
  same arithmetic gives **6.69 bits, a 0.95× ratio**. The headline does not survive its own
  sprint's correction and the report carries no `d` at all.
- **§A.7's "the effect does not reach the tail" quotes the outcome-defined stratum and drops the
  filter-independent one, which reverses both legs.** Both rows are in `s30_G_disp2.json` and both
  are in S30-L26. The ledger printed both; the report printed one.
- **§4.1's "the 18 worst pools" is FAIL18, which is not the 18 worst.** On the genuinely worst 18
  the three quoted numbers are 2.5298 / 11 of 18 / 3.95, not 2.2842 / 13 of 18 / 3.54.
- **§4.1's "it replicates on all three tail definitions" is false as written.** The two
  filter-independent tails give ρ = 0.380 and 0.344, not 0.107.

Beyond those, §4.3's bolded "anti-aligned on the hard targets" is **0.67× its own MDE** — below
the sprint's own "not a result" line — and §6.1's "A control must match the operator's own space |
no instance this sprint" is falsified by §4.1's retrieval paragraph.

The report is *less* careful than its own LEDGER in at least four places (A.7 tail rows, the
shape-83% quotation condition, the F3c "statistic guessed" caveat, the S30-L15 `d` flag). Every
one of those is a case of the summary dropping a qualifier the lane attached on purpose. That is
the pattern to fix before §0 is written, because §0 will inherit it.

---

## CONFIRMED

Recomputed and matching. Where I recomputed from raw data rather than the lane's JSON it says so.

### §1 — the opening arithmetic: **all of it, exactly**
Recomputed from `s29/results/s29_O_chain_rows.jsonl` (`item == "prod"`, `rmsd_chain`), n = 126:

```
mean 3.2105   median 2.9661   max 8.2406
worst-18 mean 6.2758    other-108 mean 2.6997
cap worst 10 at 3.00 -> 2.9074 (-0.3031)
cap worst 18 at 3.00 -> 2.7426 (-0.4680)
cap worst 30 at 3.00 -> 2.5778 (-0.6328)
all 126 minus 0.20   -> 3.0105
```
Basis is correctly the **built chain**, and the stratum is correctly the worst-18-by-outcome.
§1 is clean. (One caveat about a line §1 *omitted* is in DEFECTS, D17.)

### §4.1 — 2.2842 Å, recomputed independently from raw
Priority item 1. Rebuilt from `s12.instrument` (`load_univ`, `pool_idx`, `shipped_record`) over
all 126 targets, not from `s30_F_stagegap.json`:

```
ORACLE best-of-pool-500   all 1.7108   FAIL18 2.2842   other108 1.6153
ORACLE best-of-universe   all 1.3134   FAIL18 1.6400   other108 1.2590
ORACLE best-of-top-75     all 2.3062   FAIL18 4.6774   other108 1.9109
FAIL18 with best-of-pool < 3.00:  13    worst: 3.5436
```
**2.2842, 13 of 18 and 3.54 all reproduce exactly — for the FAIL18 stratum.** The built-chain
counterpart 2.2845 is `s29_O_ladder_table.json` rung "1 best single member, K=500", `chain_fail18`
= 2.28453. The *number* is right and its basis pair is right. What is wrong is the label the
report puts on the set (D3).

### §4.4 — ΔR² −0.089 / +0.600, recomputed
Priority item 2. `s30_R_verdict.json → D1_locality`:

```
r2_m_only  0.31106   r2_local 0.22164   r2_global 0.91077
d_local   -0.08942   se 0.01264  fold CI [-0.1217,-0.0522]  2.53x MDE  5/5 folds
d_global  +0.59971   se 0.00839  fold CI [+0.5850,+0.6166] 25.52x MDE  5/5 folds
```
Both reproduce, both exclude zero on the fold CI, both 5/5. **The spine of "per-residue channels
are closed by theorem on this instrument" holds.** Two small corrections of wording are in D12.

### §3 — rows that reproduce
| row | check |
|---|---|
| ORACLE global ρ 0.169, LFO 0.012, best single 0.121 | `s30_D_gram.json`, verifier, MATCH |
| sparse 2-of-75 = 2.1683 vs argmin-128 = 2.1435 | `s29_O_ladder_table.json`, **both built chain**, MATCH |
| quadric +0.2059 (1.81× MDE) and +0.086 (2.95×) | S30-L12 / S30-L15, MATCH |
| `disp2` 3.3585 vs 3.0483 | S30-L12, **both CA cloud**, MATCH |
| m_eff = 1.4 of 75 | S30-L7 (224.18/160.36 = 1.398), MATCH |
| E2: 98% of clashes for +0.08 Å; 2/n = 15.4% at n=13 | S30-L13 / ledger 1661, 2505, MATCH |
| filter width: no knee, cross at ~55%/55% near k≈275 | S30-L9 / ledger 1405, MATCH |
| ORACLE global argmin over k **is** the shipped 75 | `s30_F_width_cloud.json phase2_cloud`: k75_m75 = 3.04834 is the minimum of all 41 cells. MATCH (cloud basis, see D10) |
| ρ = 0.358 needed for 3.00 Å | S30-L14 §6, MATCH |

### §4.2 — the ladder and the law
`s30_T_bits.json`, 126 rows, refitted from scratch:
```
D(R) = 4.1080 3.6059 3.1115 2.7289 2.4901 2.2973 2.1230 1.8978 1.8061 1.7108
my refit:  a 1.3176  c 2.7981  gamma 3.1917  R2 0.998573
reported:  a 1.3312  c 2.7859  gamma 3.1636  R2 0.998341
ln2/gamma = 0.2172 (mine) / 0.2191 (reported);  marginal at R=7 = 0.1329 / 0.1317
rho from D0=4.1080, D(7)=1.8978 : 0.88690  (reported 0.887)
```
The ladder, ρ, the functional form and the derivative coefficient all reproduce; the small `a`/`γ`
difference is a fitting-weight choice, not an error. **Candidate indexing 0.132 vs subset
cardinality 0.044 = 3.0×** reproduces.

### §4.3 — the radial measurement (the descriptive half)
`s30_D_radial.json`, 126 rows, recomputed:
```
radial_share      mean 0.5798  median 0.5879   16.81x MDE
cos_pc1_radial    mean 0.9472  median 0.9840   30.92x MDE
stable            mean 1.7050  ->  stable_after mean 2.6425
cos_u_radial      all mean -0.0675      FAIL18 mean -0.2524
```
All six printed values reproduce exactly. The *interpretation* of two of them is defective (D5, D8).

### §4.4 — the rest of the recognition verdict
`s30_R_verdict.json`: 43 channels; DIS `rho_A_part` +0.3468 with `rho_ANCHOR_part` +0.2132 →
contrast +0.1336; `A4_multiplicity.rho_A_part.p_max = 0.0`; RAMA `pref_near` 0.6401 and
`pref_pool` 0.7896; DIS `pref_near` 0.0584 → prefers production on 94.16%; instrument
`near_rmsd` 0.5553 ("a 0.55 Å structure"); LFO real `pref_near` 0.92982, `pref_pool` 1.0, margin
−0.070175, fold CI [−0.11009, −0.02830]; `LFO_garbage_check.real.pref_near` = 1.0 (the "3 Å rung
on 100%"); DIS `conc_near["0-0.25"]` = 0.52035. All MATCH.

### §4.5 — lane X
`s30_X_typicalgood.json`: `TERTIARY_avg_gain_vs_S` slope 0.414269, intercept −0.155912, r 0.885417
→ "avg_gain = 0.4143·S − 0.1559, r = 0.885" MATCH. `MOST_CONCENTRATED_SOURCE`: arm T0_helix,
S 0.56764, pool S 1.57700, endpoint 3.78916, `is_it_the_worst_endpoint` true. MATCH.
`set_mean² ≈ B² + S²` holds to ~2.5% on all four arms (3.503 vs 3.592 on T3_pool, 3.831 vs 3.866
on T0_helix) — the "≈" is doing real work but is honestly written.

### Appendix A — rows that reproduce
A.1's 0.3676/0.3688 pair (verifier); A.1's 1.86 pair-distance vs 3.4–3.6 coordinate (ledger 923–925,
2184–2189); A.2's 3.78-bit correction (S30-L15 §3); A.3's −0.128 [−0.316,+0.036] (`s30_F_score.json
→ predictors_of_rho_pool.F4_pool_rg_sd`: −0.12780, fold CI [−0.3163,+0.0362]); A.3's error×confidence
−0.729 vs −0.799 (`F3d`); A.4's 24.80% vs 25% (`PRIMARY_falsifier`: 0.248016 vs 0.25) and the
−0.4744/+0.0851 pair; A.5's F1c FAIL18 +1.7674 (`F1c.filter_fail18.effect` = 1.767381).

### §A.7 — the WRITHE self-kill: **sound, and the best-argued thing in the appendix**
`contrast_pref` +0.1641 [+0.1016,+0.2282], 2.17× MDE, 5/5 (ledger 3593); kind-matched rebuilt-native
percentiles 0.6061 / 0.6732 against DIS's 0.2876 (ledger 3608–3612); best chiral anchor contrast
+0.0405 against a max-over-3 null mean of +0.0408, p_max 0.430 (ledger 3553, 3563). The
`reflection_audit` in `s30_G_chiral_rows.jsonl` confirms the channels are genuinely chiral
(`flip_err` 0.0, `rot_err` ~1e-13). **No defect.** One mis-count of instances is in D13.

### §A.7 (a) — the reversal's stated reason: **SUPPORTED by the artefacts**
You asked me to attack this rather than accept it. It survives, with one omission.

- Lane F's −0.128 is `F4_pool_rg_sd` vs **`rho_pool`** — the score's Spearman with ORACLE in-pool
  RMSD, i.e. *filter ordering skill*. Source: `s30_F_score.py:96,154`, `s30_F_score.json →
  predictors_of_rho_pool`.
- Lane G's −0.316 is `DISP` vs **`d = k30 − proj`** — the AMBER relax gain. Source:
  `s30_G_disp.py:180`, `s30_G_disp.json → F_G1_continuous`.
- These are **different outcome variables**, exactly as you wrote. Lane L itself said so at the
  time (ledger 2630: *"Different outcome variable, so not a refutation, but…"*) and lowered the
  claim anyway. Lane G's own prior is documented in S30-L26 as resting on the common-mode /
  bias-repair argument, which is a different mechanism.

Your sentence is correct. **The omission**: it was also a *different predictor*. Lane F's artefact
carries two dispersion variables — `F4_pool_rg_sd` (−0.128, NOT MEASURED) and `F4_top75_rg_sd`
(−0.528, fold CI [−0.602,−0.452], well measured). Lane L withdrew on the pool-level one; lane G
confirmed using the top-75 one (`s30_G_disp.py:187` computes `rg.std(ddof=1)` over `W75`). So the
honest version is one clause longer: *different outcome, and the weaker of two dispersion
variables in the same artefact — the stronger of which is the one lane G later used.* That makes
the withdrawal look worse, not better, so it belongs in the appendix.

---

## NOT RECOMPUTABLE

| claim | why |
|---|---|
| §3 "generative structural spaces … at ρ≈0 the endpoint is near-unit-slope in the set mean, and width buys ceiling while costing mean" | Prose. No artefact key holds a slope for this; `s30_X_sourcelaw.json` was not cross-checked to this wording. Label it a synthesis or attach a number. |
| §3 "subset objective through an averaging readout — T1: the tail is *always* a prefix of the order induced by ∇V at the optimum" | A theorem statement. Not recomputable by construction; correctly presented as such ("Several are theorems rather than measurements"). |
| §3 "E1 — a prior on the shared bias's form | L | three independent ways" | Prose-only; no artefact enumerates the three. |
| §4.3 "a uniform 10% contraction of a 13-mer … moves it ~14 kT, while a size-matched reference moves by exactly 0.0000" | Lane L derivation; `s30/lit/s30_L_orthogonality.py` exists but I did not find a committed artefact carrying 14 kT. Prose-derivation, flag as such. |
| §4.3 radial-share CI `[+0.545, +0.603]` | Not in `s30_D_radial.json` (per-target rows only). It is in the LEDGER table (line 2339) as a fold CI. Recomputable in principle, not from the published artefact. The iid 95% band from the rows is [0.5555, 0.6041] — compatible. |
| §4.2 R² 0.9983 / a 1.3312 / γ 3.1636 exactly | The fit is not stored in `s30_T_bits.json` (only `rows`, `NS`, `NDRAW`). My independent refit gives 0.998573 / 1.3176 / 3.1917. Same conclusion, different third digit — the fitting procedure is not in the artefact. |
| §4.1 "the chain is fully mediated" | Ledger 1811 (+0.299, sign flips) — recomputable from the ledger, **not** from any `s30/results/*.json` key. The mediation numbers are not persisted. |
| §A.7 "97.2% of the gain in the divergent half" quoted next to "−0.0406, 3.56× MDE, length-matched" | Both reproduce but from **different splits**: 97.2% is `s30_G_disp.json → F_G1b_split.DISP_rmsd` (mean_hi −0.042913 / mean_lo −0.001234, the plain split), while −0.0406 / 3.56× is `s30_G_disp2.json → length_matched_split`. Two splits, one sentence. See D10. |

---

## DEFECTS

Ranked by severity.

### D1 (SEVERE) — §A.7: "the effect does not reach the tail" is argued on the outcome-defined stratum, and the filter-independent stratum in the same artefact reverses both legs

**Report says:** *"What survives the confirmation is narrower than the claim: **the effect does not
reach the tail.** FAIL18 dispersion against the 108 is NOT MEASURED (0.47× MDE), the relax gain on
FAIL18 is −0.0113 against −0.0239 on the other 108 — *the wrong direction* — and targeting the
divergent half buys −0.0215 against −0.0221, i.e. nothing."*

**Artefacts say** (`s30_G_disp2.json → is_divergence_the_tail`, and `s30_G_disp.json →
F_G1c_strata`):

```
                                        effect    xMDE   fold CI               folds
FAIL18 dispersion vs the 108           +0.4583   0.47   [-0.3845,+1.0472]     2/4   NOT MEASURED
ORACLE worst18-by-poolmean dispersion  +1.6116   1.44   [+0.5832,+2.0755]     4/5   EXCLUDES ZERO
relax gain: FAIL18 -0.0113 vs 108 -0.0239        0.22x MDE on the difference
relax gain: ORACLE worst18 -0.0386 vs rest -0.0193   0.31x MDE
FAIL18 in the high-dispersion half: 11 of 18 (9.0 expected if independent)
```

**Both legs reverse on the non-circular stratum.** On the ORACLE worst-18-by-pool-mean the tail
*is* measurably more divergent (1.44× MDE, CI excluding zero) and the relax gain is *twice* as
large there (−0.0386 vs −0.0193).

**My ruling on the circularity you asked me to make:** this is **not** the S30-L23 forced-arithmetic
circularity. Lane G is right and declared it in its prereg — FAIL18 is defined by the filter's
in-band recall, and the relax gain is not the filter's recall. **But the stratum is worse than
circular for this question: it is the wrong tail.** FAIL18 is barely enriched in high dispersion
(11 of 18 vs 9 expected) while the outcome-relevant tail is strongly enriched (+1.61 vs +0.46).
Using FAIL18 to ask "does a dispersion-graded effect reach the tail?" is biased toward the null by
construction of the stratum, and the report then reads that null as evidence of absence.

**Lane G's own artefact carries the warning the report dropped:** `"FAIL18 is defined by the
FILTER's in-band recall (S30-L23) … reported as a descriptive split"`. So does the ledger, which
printed **both** rows side by side. The report printed one.

**Corrected sentence:**
> What survives is narrower than the claim, and the tail question is **NOT MEASURED on either tail
> definition**, with the two definitions pointing opposite ways. On FAIL18 — the filter-defined
> stratum, reported by lane G as descriptive — the gain is smaller (−0.0113 vs −0.0239, 0.22× MDE)
> and dispersion is not elevated (0.47× MDE). On the ORACLE worst-18-by-pool-mean, dispersion **is**
> elevated (+1.61, 1.44× MDE, fold CI excluding zero) and the gain is **larger** (−0.0386 vs
> −0.0193, 0.31× MDE). What is established is that targeting buys nothing (−0.0215 against
> −0.0221) and that the whole arm is 0.69% of baseline against a counterfactual asking for −0.30 Å.
> E2 stays out of the "mechanism found" column on **size**, not on reach.

Note the last clause: the targeting-null leg (leg 3) is sound on its own and is sufficient for the
conclusion you want. You do not need legs 1 and 2, and they are the two that do not survive.

### D2 (SEVERE) — §4.2: the "5.2× / worth five times their face value" headline is a free parameter its own sprint corrected downward

**Report says:** *"Seven index bits move the ORACLE ladder 4.108 → 1.898 Å, which through the
displacement bound is ρ = 0.887 — **36.6 bits of displacement information out of 7 index bits, a
5.2× ratio**."* and the pull-quote *"the readout's 7 bits are worth five times their face value."*

**What the number is.** ρ is obtained by inverting S29's displacement bound,
`ρ = sqrt(1 − (D/D0)²)`; bits are then `I = −(d/2)log2(1−ρ²)` with `d = 3n−6 = 32.88` (S30-L14 §6).
I reproduce 36.633 at that `d`. **The report states neither the formula nor `d`.** Sensitivity:

```
d = 32.88 (3n-6, nominal)               I = 36.63 bits   5.23x   <- what the report quotes
d = 21    (k99, pair-distance)          I = 23.40        3.34x
d = 10                                  I = 11.14        1.59x
d =  6    (k90 / the pool's top-6)      I =  6.69        0.95x   <- headline inverts here
d =  2.72 (deviation-matrix stable rank) I =  3.03        0.43x
d =  2.06 (Gram stable rank, S30-L21)   I =  2.29        0.33x
```

**The sprint already made this correction and the report did not carry it.** S30-L14 §6 ends with
the line *"**S30-L15 corrects the `d` this should be evaluated at**, and the corrected numbers are
smaller."* S30-L15 §3 then evaluates the identical formula in the pool's top-6 subspace and
**withdraws a 12-bit figure down to 3.78** on exactly this ground — a withdrawal the report
records in A.2. So §4.2 reproduces a number whose `d` the appendix, two pages later, records as
corrected.

I am not claiming `d = 6` is the right `d` for *this* quantity — the top-6 correction is about
direction information and this is total displacement, and a case can be made either way. I am
claiming the report states a headline that moves from 5.2× to 0.95× under a change of assumption
its own ledger flags, without stating the assumption.

Secondary: this is also the §7-of-my-charge pattern. ρ = 0.887 is an inversion of the
`RMSD·sqrt(1−ρ²)` bound, which lane P measured as flattering (implied 3.0338, applied 3.0519,
production 3.0483). Here the inversion direction is conservative in ρ, so the bound is not the
problem — `d` is. But "36.6 bits" is still an implied conversion of an implied conversion, with
no applied measurement anywhere behind it.

**Corrected sentence:**
> Seven index bits move the ORACLE ladder 4.108 → 1.898 Å (CA point cloud), which through S29's
> displacement bound is ρ = 0.887. Converting that to bits requires an effective dimension: at the
> nominal `d = 3n−6 = 32.88` it is 36.6 bits, 5.2× the index budget; at the pool's own top-6
> effective dimension (S30-L15) it is 6.7 bits, 0.95×. **The qualitative claim — that bits are not
> conserved across an index, because the 500 deposited backbones carry the structure and the index
> only names one — does not depend on `d`. The 5.2× multiplier does, and should not be quoted
> without it.**

The qualitative claim in §4.2 is robust and is the part worth keeping. The multiplier is not.

### D3 (SEVERE) — §4.1: "the 18 worst pools" is FAIL18, and FAIL18 is not the 18 worst

**Report says:** *"The ORACLE best member of **the 18 worst pools** is 2.2842 Å (built chain
2.2845) … with **13 of the 18** holding a member under 3.00 and the worst tail pool bottoming out
at **3.54**."*

**The LEDGER says it correctly** (S30-L2 §2: *"the ORACLE best member of the **FAIL18** pools"*).
The report reworded it, and §1 four pages earlier defines "the worst 18 targets" as an order
statistic on production RMSD. A reader carries §1's definition into §4.1.

**They are different sets.** Recomputed:
```
worst-18-by-production(chain) vs FAIL18 overlap: 13 of 18
  in worst18 not FAIL18:  2LM8 2MAI 2MFV 2MSA 7YFS
  in FAIL18 not worst18:  1JBF 1LB7 2NB7 3BTB 9L1M

ORACLE best-of-pool on the genuinely worst 18:  2.5298
  under 3.00:  11 of 18       worst:  3.9523
```

The direction of the claim survives (2.5298 < 3.00), so **§4.1's headline "the tail is
selection-limited, not pool-limited" is safe**. Three numbers are not.

Also: *"Nothing in this claim is conditioned on the thing it measures"* is stated too strongly.
FAIL18's defining predicate (`band = rr <= rr.min() + 1.5` must miss `sub`, `s12/instrument.py:281–283`)
is a relation between the pool minimum and the top-75, and the claim is about the pool minimum.
The conditioning is weak and I could not sign it, but "nothing" is not the right word.

**Corrected sentence:**
> The ORACLE best member of the FAIL18 pools is 2.2842 Å point cloud (2.2845 built chain), with 13
> of the 18 holding a member under 3.00 and the worst bottoming out at 3.54. On the *outcome*-defined
> tail — the 18 worst targets by production RMSD, the set §1's arithmetic is about — the ORACLE
> best-of-pool is 2.5298 Å, 11 of 18 under 3.00, worst 3.95. Both are below the 3.00 cap in the
> mean. The material to fix the tail is already inside the candidate sets the pipeline is handed,
> on either definition.

### D4 (SEVERE) — §4.1: "it replicates on all three tail definitions" is false as written

**Report says:** *"**The strongest tail result, and it replicates on all three tail definitions:**
… +0.6446 on the easy 108 … and **+0.1066 on the hard 18 with the fold CI including zero.**"*

**`s30_F_score.json → F3a` has all three, and they are not the same number:**
```
other108              rho_pool 0.6446   fold CI [+0.5905,+0.7063]  5/5 folds
FAIL18                rho_pool 0.1066   fold CI [-0.0241,+0.1916]  3/5 folds
worst18_poolmean      rho_pool 0.3798   (no CI stored)
worst18_bestpool      rho_pool 0.3438   (no CI stored)
random-18 null:  mean 0.5682   CI95 [+0.4067,+0.7170]   p = 0.0 for FAIL18
```

The number the report quotes is the **filter-defined** one, and it is 3.2–3.6× lower than the two
filter-independent ones. The "fold CI including zero" exists **only** for FAIL18; no CI is stored
for the other two, so that half of the sentence cannot replicate at all. And 0.1066's 3/5
folds-same-sign is not reported.

What *does* replicate is weaker and still real: 0.3798 and 0.3438 both fall below the random-18
null's lower bound of 0.4067, so the score is genuinely worse at ordering hard pools on every
definition — just at roughly a third of the reported degradation, and much closer to the
all-target 0.5678 than to zero.

This is the S30-L23 failure mode in a softer form: the headline number comes from the stratum
defined by the filter's own recall, and the claim is about the filter's ordering.

**Corrected sentence:**
> The shipped score's Spearman with ORACLE in-pool RMSD is +0.6446 on the easy 108 (fold CI
> [+0.5905,+0.7063], 5/5) against a random-18 null of 0.5682 [+0.4067,+0.7170]. It degrades on
> every tail definition — **+0.1066 on FAIL18** (fold CI including zero, 3/5 folds; but FAIL18 is
> the filter's own recall predicate, so read it as a bound not an estimate), **+0.3798** on the
> worst 18 by pool mean and **+0.3438** on the worst 18 by best-of-pool, the latter two below the
> null's lower bound with no fold CI computed. *The score orders its own pool markedly worse on
> hard targets; on the non-circular definitions the effect is about a third of the FAIL18 figure.*

The italicised conclusion *"The score cannot order its own pool on hard targets"* should be
softened to "orders it markedly worse" — at ρ ≈ 0.35 it can still order it.

### D5 (MAJOR) — §4.3: the bolded "anti-aligned on the hard targets" is 0.67× its own MDE

**Report says (bold pull-quote):** *"…orthogonal to the answer in general and ***anti-aligned* on
the hard targets**"*, resting on `cos(direction to native, radial), FAIL18 = −0.2524`.

**Recomputed from `s30_D_radial.json`:**
```
all 126:  mean -0.0675   sd 0.3636   se 0.0324   MDE 0.0907   ->  0.74x MDE
FAIL18:   mean -0.2524   sd 0.5672   se 0.1337   MDE 0.3745   ->  0.67x MDE
          per-fold means  [+0.0048, -0.0050, -0.1667, -0.6491]   3 of 4 negative
          fold CI95 ~ [-0.5048, +0.0968]  -- includes zero, and one fold carries it
```

**0.67× MDE is below the sprint's own 0.7× "not a result" line** (§6.1, row 2). There is no CI for
this cell anywhere in the ledger table either — the FAIL18 column of S30-L22 is the only column
without one. The per-fold breakdown shows a single fold (−0.6491) driving it, which is the same
shape as lane X's withdrawn −0.4744 (A.4).

The "orthogonal in general" half is fine — 0.74× MDE supports a null and that is what is being
claimed. The "anti-aligned on the hard targets" half is a positive claim from a not-a-result.

**Corrected sentence:**
> The library spends the majority of its two available directions on a component that is orthogonal
> to the answer overall (cos −0.0675, 0.74× MDE — consistent with zero). On FAIL18 the point
> estimate is −0.2524, but at **0.67× its own MDE with one of four folds carrying it, that is NOT
> A RESULT** and must not be quoted as anti-alignment.

Losing it costs little: "orthogonal to the answer" plus lane F's shape-not-scale result already
carries §4.3's conclusion. But "Three lanes, one chain" then becomes two lanes.

### D6 (MAJOR) — §4.1: "Retrieval is exonerated" — the control is not in the operator's space, and §6.1 claims no such instance this sprint

**Report says:** *"Retrieval is exonerated at every stratum: BLOSUM's 500 against a random 500 of
the same universe is NOT MEASURED everywhere, most pointedly on the tail (0.01× its own MDE).
**This revises what the project had recorded** — the harm on hard targets was attributed to the
retrieval corpus; it is not retrieval."*

Three problems, in increasing order:

1. **NOT MEASURED is not exoneration.** `F1c.retrieval_fail18`: effect −0.0034, **MDE 0.2902**,
   power 0.050, type_m 71.3, folds_same_sign 1/4. That interval admits ±0.29 Å of retrieval effect
   on the tail. The other two strata are at 0.87× and 0.99× MDE — also in the NOT MEASURED band,
   with point estimates *favouring* BLOSUM. The correct reading is "the sprint could not measure
   it", not "retrieval is exonerated".
2. **The operator spaces do not match.** The record being revised (`sequence-conditioning-hurts-the-failures`)
   is about the **end-to-end pipeline endpoint** on FAIL18: blind 5.425 vs shipped 6.019, a 0.594 Å
   built-chain gap. The new measurement is the **ORACLE best-of-pool**, CA cloud. Retrieval's effect
   on a pool's ORACLE ceiling and its effect on what the pipeline emits from that pool are different
   quantities — the harm can be entirely downstream (the score orders BLOSUM pools worse) with the
   ceiling untouched. **This is the project's most repeated error, and §6.1 states "no instance this
   sprint".** That row is falsified by this paragraph.
3. Bases are mixed across the revision: 0.594 Å built chain vs an MDE of 0.2902 CA cloud.

**Corrected sentences:**
> At the ORACLE best-of-pool level, BLOSUM's 500 against a random 500 of the same universe is NOT
> MEASURED at every stratum — 0.87×, 0.99× and 0.01× of their own MDEs, with all three point
> estimates favouring BLOSUM. **On the tail the MDE is 0.29 Å, so this bounds nothing smaller than
> that and is not an exoneration.** It also does not speak to the record's claim that sequence
> conditioning hurts the failures (blind 5.425 vs shipped 6.019 built chain), which is a statement
> about the *emitted endpoint*, not about the pool's ceiling. That remains untested.

And §6.1's row should read *"the project's most repeated error; one instance this sprint, §4.1's
retrieval paragraph"* — or the paragraph should be fixed and the row kept.

### D7 (MAJOR) — §4.1: the shape-83% number is quoted against a condition its own lane attached, and without its "statistic guessed" caveat

**Report says:** *"**shape error is 83% of it — scale is refuted as the mechanism** (partial ρ =
−0.067, p = 0.46, against shape's −0.641 at p = 6.5e−16)."* and repeats "83% of the tail's error"
in §4.3.

The numbers reproduce (`F3c`: −0.067045 p 0.4557; −0.640711 p 6.49e−16; 83% = shape_err mean on
FAIL18 4.3670 → 3.6304/4.3670 = 0.8313). Three omissions:

1. **Lane F attached a quotation condition and the report broke it.** S30-L16 §3, verbatim:
   *"|scale error| is 3.23× on FAIL18 versus shape's 1.93×, but shape is 83.1% of the tail's total
   error (94.2% on the 108) … **Both facts belong in any quotation of this row**; the elevated
   ratio is real and the causal claim it invites is not."* The report quotes the 83.1% and drops
   the 3.23×.
2. **The artefact records `F3c.fires = False`** with `caveat: "statistic guessed by lane F; lane L
   was unreachable to specify it"`. An unregistered, self-chosen statistic is being used as the
   headline refutation of another lane's derivation. That caveat is not in the report.
3. 83% is the **FAIL18** figure; the all-target figure is 91.2% and the other-108 figure 94.2%. The
   report's §4.3 pull-quote "shape … is 83% of the tail's error" is correct; §4.1's "83% of it"
   with "it" ambiguous between rho_pool and the tail's error is not.

**Corrected sentence:**
> Shape error, not scale, carries it: partial ρ(rho_pool, shape | scale) = −0.641 (p = 6.5e−16)
> against ρ(rho_pool, |scale| | shape) = −0.067 (p = 0.46), and shape is 83.1% of the tail's total
> error. **Both halves of lane F's row belong here: |scale error| is elevated 3.23× on the tail
> against shape's 1.93×, so scale is elevated and still is not the mechanism. The partialling
> statistic was chosen by lane F, not pre-registered — its falsifier is recorded as not firing.**

### D8 (MODERATE) — §4.3: "fraction of lambda_1 that is radial 94.8%" is a ratio of two trace shares, not a projection

`94.8% = 0.5798 / 0.6117` = (radial share of the Gram trace) / (λ₁'s share of the trace),
ledger 2350. That quantity is **not** "the fraction of λ₁ that is radial" — it compares the radial
component's share of the *whole trace* to λ₁'s share of the *whole trace*, and the radial component
need not lie inside the λ₁ eigenspace. It is unbounded above in principle (it happens to be < 1 on
all 126 here). It is also a ratio of two means, where the per-target mean of the ratio is 0.9417.

The quantity the label describes is `E[cos²(pc1, radial)] = 0.9122`.

**Corrected line:** `radial share of lambda_1's own direction   0.912  (mean cos^2)` — or keep
0.5798/0.6117 and rename it *"radial share of the trace, relative to λ₁'s share"*.

### D9 (MODERATE) — §4.5: `corr(S,B) = +0.085` is the arm-dropped value, presented as the primary

`s30_X_typicalgood.json → SECONDARY_within_target_corr_S_B` = **−0.4744, n = 630**. The +0.0851
(n = 504) is what remains after dropping T0_helix, and A.4 discloses that. **§4.5 does not**, and
§4.5 is where a reader takes the number from.

Note the raw value *met* lane X's registered prediction ("near zero or negative"). Dropping the arm
moves it to the other side of zero and the report uses the post-drop value for a claim
("concentration is orthogonal to bias") that the pre-drop value also supports, differently
("concentration is anti-correlated with bias"). The conclusion is robust; the presentation is not.

**Corrected sentence:** *"…and within target `corr(S,B) = −0.474` over all 630 cells, +0.085 over
the 504 with the T0_helix arm dropped as an artefact (A.4) — **concentration does not buy bias on
either reading**."*

### D10 (MODERATE) — bases are omitted where they matter, and three different "production built chain" values are in circulation

The report's header pins the endpoint as built chain, and §1 uses 3.2105. Then:

| location | number | actual basis | stated? |
|---|---|---|---|
| §4.2 ladder 4.108 → 1.898 → 1.7108 | pool ORACLE `rr` | **CA point cloud** | no |
| §4.5 "the worst endpoint in the record, 3.789" | vs pool 3.0483 | **CA point cloud** | no |
| §3 filter-width row ("argmin over k IS 75") | `phase2_cloud`, artefact `basis = POINT CLOUD; built chain is phase 3` | **CA point cloud** | no |
| §3 sparse row 2.1683 / 2.1435 | ladder table | built chain | no |
| §3 `disp2` 3.3585 vs 3.0483 | | CA cloud | yes (via 3.0483) |
| §4.1 2.2842 / 2.2845 | | both, paired | yes |

**§3's table mixes bases row to row with no column for it.** A reader who has just read §1 will
read §4.5's 3.789 against 3.2105 and get a 0.58 Å error.

And three values for the same quantity are live in the sprint: **3.2105** (`s29_O_chain_rows`,
the report's endpoint), **3.2071** (`s30_D_meter`, §10.2, and `s30_R_verdict.instrument.prod_rmsd`),
**3.2052** (`s16/results/repair_A.json`, `proj` mean, the baseline lane G's relax arm is measured
against, while `endpoint_arithmetic` divides by 3.2105). Immaterial numerically; material for a
report that has published the wrong basis before. Pin one and footnote the other two.

Also in this family: §A.7 puts "−0.0406, 3.56× MDE, 5/5, length-matched" and "97.2% of the gain in
the divergent half" in one sentence; they come from two different splits (`length_matched_split`
vs `F_G1b_split.DISP_rmsd`, whose own ratio is 1.48× on the fold SE). Say which is which.

### D11 (MODERATE) — §4.2 / §3: "48 bits for four basins per residue at n = 12" is a factor of 2 off as worded

The ledger's arithmetic (S30-L14 §5) is `2n·log2(k)` over **2n = 25.9 torsions**: 25.9 at k=2,
41.1 at k=3, **51.8 at k=4**. The report's 48 = 2·12·log2(4) is the same formula evaluated at
n = 12 exactly rather than at the measured mean n = 12.96 — so 48 is arithmetically consistent with
*four basins per **torsion***. As written ("four basins per residue at n = 12") it computes to
**24**, not 48.

The conclusion is unaffected (24 or 48 or 51.8, all against 7 deployed). But this is a headline in
two places and a reader will check it.

**Corrected phrase:** *"arithmetically infeasible — naming one Ramachandran basin per torsion costs
2n·log2(k) bits: **51.8 at four basins and n = 12.96**, against **7** deployed."*

### D12 (MODERATE) — §4.4: three small mis-statements

1. **The LFO margin −0.070 is 1.04× MDE, 4/5 folds** (`LFO_combination.real.margin_ratio` =
   −1.0424, `margin_mde` 0.0673, `folds_same_sign` 4). The ledger states both; the report states
   neither. At 1.04× this sits one hundredth above the NOT MEASURED band, in a sprint whose §6.2
   says the surviving positives are the 2–10× ones. (The *conclusion* rests on `pref_pool = 1.0`,
   not on the margin, so it holds — but the margin must carry its ratio.) Also unmentioned: the
   shuffled control gives margin +0.114 with a fold CI **including zero**, which is the right
   control and passes.
2. **"reaching 0.822 only above 4 Å" is from the all-pairs table, not the near band.** DIS
   `conc_near` has bins 0–0.25 … with `n/a` at the top; 0.822 is `conc[">4"]` (ledger 2802 ALL
   PAIRS row: 0.525 0.566 0.618 0.698 0.771 0.822; the near-band row 2806 is 0.520 0.554 0.570
   0.624 n/a n/a). The sentence as written implies the near band reaches 0.822, which is empty by
   construction. Corrected: *"DIS separates two structures 0–0.25 Å apart at 0.520 inside the near
   band, and over all pairs only reaches 0.822 above 4 Å."*
3. **"The largest preference effect in the library"** is CONS at |pref_near − 0.5| = 0.4433
   (pref_near 0.0567); DIS is second at 0.4416. Say "one of the two largest" or name CONS. (CONS is
   a pool-reference channel, so quoting DIS is defensible — say why.)
4. Minor: "still adds nothing beyond knowing the perturbation budget" — ΔR² = −0.089 held out means
   the local block *subtracts*, i.e. overfits. "Adds nothing" understates a cleaner result.

### D13 (MODERATE) — §A.7's "third instance of the cross-kind confound" contradicts A.1 and A.5

A.7: *"**Third instance of the cross-kind confound this sprint**, after S28-L48 (withdrawn by lane
R) and the widening result (withdrawn by lane F)."*

A.1 and A.5 both describe the widening result's defect as **circularity / an outcome-defined
stratum**, not cross-kind: *"the widening-rescues-the-tail result does not replicate on either
filter-independent tail"* (A.1) and *"its widening result flagged by itself as circular"* (A.5).
Those are two different failure modes and A.6 lists them as two different checklist entries.

**Corrected:** *"Second instance of the cross-kind confound this sprint, after S28-L48."*

### D14 (MODERATE) — the sprint's "one cross-lane reproduction" is ρ = 1.0000 by construction

`s30_G_disp.json → cross_check_dispersion.DISP_rg_vs_laneF_top75_rg_sd = 1.0`, exactly.

`s30_G_disp.py:189` reads `F_top75_rg_sd` **directly out of lane F's artefact** into lane G's own
row, and `:187` computes `DISP_rg = rg.std(ddof=1)` over the same top-75 members. The two
quantities differ only by `sqrt(75/74)`, so their Spearman is 1.0000 **by construction**. It is a
genuine and useful check that lane G's top-75 membership matches lane F's — it would have caught a
mis-indexed pool — but it reproduces a *variable*, not a *result*, and it is not an independent
construction.

**Corrected:** *"lane G's top-75 membership is confirmed identical to lane F's (rank correlation
1.0000 between the two computations of top-75 Rg dispersion — exact by construction, since they
differ only in `ddof`). No independent reproduction of a result was performed this sprint."*

### D15 (MINOR) — §3 quotes a stable rank without its feature space, which is the exact error A.1 withdraws

§3: *"Gram stable rank **2.057** — the fields span ~2 directions, **not 11**"*.

- The project has at least seven stable ranks live: 1.681 (per-target Gram, `s30_D_gram.json`),
  1.705 (per-target, `s30_D_radial.json` — a *third* value for a per-target field rank, unreconciled
  with 1.681), 1.86 (pair-distance), 2.057 (aggregated Gram), 2.642 (radial removed), 2.717
  (deviation matrix), 3.4–3.6 (coordinate). A.1 withdraws a claim for exactly this omission.
- "**not 11**" appears to be a typo for **not 21** — there are 21 fields.
- S30-L21 also reports an **effective rank of 3.68** beside the stable rank of 2.057; the report
  quotes only the smaller one, which is the one that flatters the "closed" verdict.

### D16 (MINOR) — the verifier's coverage is narrower than its "22/22" implies, and one of its prints is wrong

`s30/s30_verify.py` now passes 22/22. **All 22 checked numbers come from lane D** — 7 from
`s30_D_gram.json`, 4 from `s30_D_radial.json`, 6 meter anchors, 4 draw contrasts, 1 sign-trap
assertion. Lanes T, X, R, F, Q, L, G and P are covered only by `show()` calls, which print and
never compare. **Every number in §4.1, §4.2, §4.4, §4.5 and every §3 row is unverified by the
verifier.** §10.1's "Verified by `s30/s30_verify.py`, 22/22 matched" is true of the anchors and
should not be read past them.

And one `show()` prints a wrong value silently:
```python
'B %s' % round(a.get('B', a.get('endpoint', 0)), 4)
```
The lane X artefact key is **`B_endpoint`**, so the verifier prints `B 0` for all four arms. The
set-mean decomposition — §4.5's entire algebraic claim — is displayed as zero and nobody noticed.
Easy fix; worth doing before anyone cites the verifier as coverage.

Cheap additions that would cover the load-bearing numbers: `F3a.other108/FAIL18/worst18_*` (§4.1),
`F1a` (§4.1), `D1_locality.d_local/d_global` (§4.4), `TERTIARY_avg_gain_vs_S` + `per_arm.*.B_endpoint`
(§4.5), the ladder `D(R)` row (§4.2), and `is_divergence_the_tail` **both** strata (§A.7).

### D17 (MINOR) — a wrong line in S30-L0 that §1 correctly dropped; do not reinstate it in §0 or §2

S30-L0's arithmetic block contains:
```
worst 18 all the way to 2.50 ->  mean 2.1334  (-1.0772)
```
That is **capping all 126 at 2.50** (I reproduce 2.1334 exactly for that operation). Capping the
worst 18 at 2.50 gives **2.6711 (−0.5394)** — the stated effect is inflated 2.0×. §1 does not carry
the line, correctly. Flagging it because §0 and §2 are still to be written from the same block, and
because the ledger is immutable so the error will sit there.

### D18 (MINOR) — §6.1's "no instance this sprint" row

Falsified by D6. Also worth noting that §6.1's row 2 ("Below 0.7× MDE is not a result; 0.7–1.0× is
NOT MEASURED") is violated by §4.3 (0.67×, D5) and by §A.7's leg 2 (0.22×, D1) in the same document.
A rule stated in §6 and broken in §4 is worse than a rule not stated.

---

## Things I checked and found clean

Recorded so you know the sweep was done, not assumed.

- **The `compare()` sign-convention class sweep you asked for.** I checked every rate, correlation,
  R², percentile, concordance and win-rate quoted in §1, §3, §4 and Appendix A against the sign of
  its own effect in the artefact. **No instance of the inverted reading leaked into a claim.**
  Specifically: lane R's `pref_near`/`pref_pool`/margin (−0.070 is correctly read as damning,
  because `pref_pool = 1.0 > pref_near = 0.930`); lane R's `conc`/`resolution`; lane G's WRITHE
  `contrast_pref` +0.1641 (correctly read as a *positive* before being killed on kind); lane G's
  rebuilt-native percentiles (0.6061 / 0.6732 correctly read as *worse* than the 0.5 line, where
  lower is better, and DIS's 0.2876 as better — that one is lower-is-better and is read
  lower-is-better); lane D's cosines (sign meaningful, read correctly); lane X's `corr(S,B)`
  (see D9 for a different problem); lane F's Spearmans (negative = more error, less skill,
  read correctly). The only direction error I found is §4.4's ΔR² wording in D12.4, which is not a
  sign inversion.
- **§4.4's "matched capacity"** — `r2_local` and `r2_global` are both increments over the same
  `r2_m_only` baseline (0.31106), and the local block is the ORACLE-advantaged one. The comparison
  is fair and the conclusion is the stronger direction.
- **Lane R's `LFO_note`** (the 0.500/0.500/0.500 null-input artefact) is in the artefact and in A.5.
  Good practice, correctly reported.
- **Lane G's `reflection_audit`** genuinely establishes the chiral channels are chiral
  (`flip_err` 0.0 against `rot_err` ~1e-13 to 1e-16). The chirality dichotomy's empirical leg is sound.
- **Lane G's mediation** (`rho(disp,d|contraction)` −0.178 vs `rho(contraction,d|disp)` +0.078)
  supports "divergence is primary, contraction is its shadow" as stated.
- **Artefact availability.** Every artefact I needed is now committed. At the start of this audit
  lane G's `s30_G_disp.json`, `s30_G_disp2.json`, `s30_G_chiral_rows.jsonl` and its four scripts
  were untracked; they were committed during the audit. Only `s30/results/s30_P_chain_rows.jsonl`
  remains untracked, and §§1/3/4/A do not depend on it.

---

## The one-line summary for §0

If §0 is written from §4 as it stands, it will inherit a bit-value multiplier that its own sprint
corrected (D2), a tail claim argued on the wrong tail (D1, D3, D4), and a "no instance this sprint"
that is false (D6, D18). **Fix D1–D6 and the report is publishable; the surviving science — the
tail is selection-limited, the pool is a codebook not a channel, the field library spends its rank
radially, ordering exists and preference does not, generation is closed jointly with the readout —
is not damaged by any of them.** Every defect above removes a multiplier, a stratum label or an
over-claim. None removes a result.
