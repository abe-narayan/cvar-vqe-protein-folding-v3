# PREREG S29-M-F2 -- CAN ANY NATIVE-FREE PREDICTOR OF THE SHELL PROFILE BREAK ASSUMPTION B2?
(lane M, written 2026-09-20 00:5x, BEFORE any F2 number exists; appended-to only, never edited
after a result exists -- contract rule 7. F1's 126-target run stays ahead of this; nothing here
has been executed.)

## 0. Why this is the live one, in the sprint's own currency

Lane T's bound (S29-L23) collapses every native-free operator onto one number: any output is
`A = c + u` with `u` native-free, and at the optimal step

    RMSD_achievable = RMSD_prod * sqrt(1 - rho_max^2),   rho = cos(u, t - c).

Assumption **B2** -- "rho_max <= 0.14 for every field constructible from the present information"
-- is the load-bearing one, and the bound says so itself. 3.00 A needs rho 0.358; 2.50 A needs
0.628.

The **shell profile** is the lowest-dimensional named quantity in the repository whose ORACLE
version clears that by a wide margin. `s12/obj_FINDINGS.md:188` (section 4's table, 126 targets,
coordinate-average arm): a scorer that knows **only the native's per-separation mean distance
profile** -- `n-2` numbers, 11 to 14 per target, no pair detail whatsoever -- emits **2.299 A**
against production's 3.048 on the point cloud. Inverting the bound's identity on that basis,
`sqrt(1 - (2.299/3.048)^2) = 0.657`: the ORACLE profile is a displacement field with a cosine of
about **0.66**, roughly **4.7x** B2's ceiling.

So F2 asks the only question that can break B2 cheaply:

> **Does a native-free predictor of the shell profile exist whose implied displacement clears
> rho = 0.14?**

If it does, the sprint's headline changes from a ceiling to an opening. If it does not, B2 is
confirmed on the one quantity that had the best a-priori chance of breaking it, which is worth more
than another arm that loses by 0.2 A.

## 1. RULE 10 DECLARATION: this is a re-opening, and the closure is TWO sections deep

Cited by line, as the contract requires.

**(a) The deployable profile scorer was measured and lost.** `s12/obj_FINDINGS.md:198`:
"shell-profile only (deployable) | 3.163 | 5.608 | 2.755 | 2.425 | 3.893 | 0.480 | -0.056 | 0.459 |
0.231" -- **+0.115 A worse than the shipped 3.048** (`:191`), with an in-band rho of **-0.056** and
the native at the **45.9th** percentile of its own pool. The ORACLE twin on the line above
(`:188`) is 2.299. **The entire gap between 2.299 and 3.163 is the profile predictor.**

**(b) And the "predict the profile better" question was asked directly, with eleven arms, and
every one lost.** `s12/obj_FINDINGS.md:306-331` (section 6, `obj_profile.json`, all leave-fold-out,
fitted on the other four folds and asserted in `fit_all`). The scored objective there swaps ONLY
the profile and keeps the distogram's pair-specific residual:
`target_p = q'_{sep(p)} + (E[d]_p - q_disto_{sep(p)})` (`s12/obj_profile.py:19`). Results:

    ORACLE true profile                                          emitted 2.402   (the ceiling)
    shipped distogram profile (baseline)                                  3.078
    LFO per-shell debias                                                  3.071
    LFO affine recal. with the pool profile as 2nd regressor  MAE 2.394   3.089   <- BEST MAE, WORSE emitted
    LFO ridge, ESM+composition+disto+pool profile                         3.309
    LFO ridge, ESM-2 PCA-32 pooled + composition                          3.339
    pool mean profile (native-free typicality)                            3.383

> "Nothing wins. The arm with the best profile MAE (cal2, 2.394) emits worse than the baseline --
> the project's MAE lesson again, now at the level of the profile."

**This is the closure F2 re-opens, and it is a strong one.** A leave-fold-out ridge over ESM-2
PCA-32 + composition + the distogram profile + the pool profile is already close to the feature set
a naive version of F2 would build, and it is **+0.23 to +0.26 A worse** than doing nothing.

### The new angle, stated narrowly and without inflation

Three things are new, and only three:

1. **The parameterisation is the RATIO, not the profile.** S12 fitted `q*` (or an affine recal of
   it) directly. F2 fits `r_s = q*_s / q_pool_s` -- the target's profile as a multiple of its own
   pool's typical profile, per shell. This matters because the pool profile already carries the
   length and the shell shape (it is the "typical peptide of this length" S18/S19 names), so the
   regression target is the *deviation from typical*, which is the quantity the bound says a field
   must correlate with. A model that predicts `r = 1` everywhere reproduces the pool profile
   exactly and is the registered zero-information control.
2. **Compactness proxies enter as named features, not as something a generic regressor might
   find.** Lane L (S29-L19) reports the in-band axis is compactness: per-target in-band skill
   correlates **+0.909** with the native's z-scored radius of gyration, and the record's own
   `in-band-ordering-is-per-target` puts achievable native-free compactness proxies at **0.24 to
   0.37**. S12's ridge had ESM means/stds and composition; it did not have an explicit predicted-Rg
   channel or the pool's own Rg dispersion.
3. **The primary readout is the COSINE, not the RMSD.** The bound's currency is rho, it is
   measurable per target without any projection, and it answers B2 directly. S12 judged these arms
   on emitted RMSD only, which is a 0.7x-MDE instrument for effects this size; a cosine at n = 126
   is far better powered and is the quantity the sprint needs.

**Registered prediction (lane T's, which I adopt): the deployable arm lands inside 0.05 A of
production.** I will be surprised by anything else. F2's value is the measurement of *how much*
of the ORACLE profile's 0.66 cosine survives leave-fold-out, which is a number nobody has.

## 2. THE TRAP THIS EXPERIMENT WALKS INTO, AND THE GUARD (contract rule 20)

**A ratio model regularised toward 1 IS the shrink the shrink rule warns about.** Rule 20 (contract
addendum 1, from lane T's theorem 2, S29-L7): shrinking a target map toward typicality,
`m -> tau + s(m - tau)`, drives the gradient cosine positive **with zero information added** while
moving the emitted structure toward the typical map. F2's ratio parameterisation is *literally*
that map with `tau = q_pool` and `s = ` the regression's shrinkage. A ridge fitted on a weak signal
shrinks toward the mean by construction, so **F2 is pre-loaded to produce a spurious cosine gain**.

Therefore, registered now, before any number: **every F2 arm reports the rule-20 signature in the
same table** -- (a) its implied shrink `s` (the sd of its predicted `r` about 1, divided by the sd
of the true `r` about 1), (b) meter number 3, the native's percentile in its own pool, (c) the
emitted structure's mean virtual bond and Rg. And the decisive control:

> **PROFILE-SHRINK CONTROL.** For each fitted arm, a zero-information twin that reproduces that
> arm's realised shrink `s` exactly, by scaling the TRUE-mean-free pool profile -- i.e.
> `r_ctrl = 1 + s * z`, with `z` a stable per-target random draw of the same sd, carrying no
> sequence information. If the fitted arm's cosine does not exceed its shrink twin's by more than
> 1.0x the MDE of that paired contrast, **the arm has bought nothing and F2 reports zero**, whatever
> the RMSD does.

This control is the point of the experiment. Without it a positive cosine here is unreadable, and
lane D would be right to veto it on sight.

## 3. THE ARMS

Pool, posterior, K = 500, M = 75, uniform average, production projection and RMSD all held at the
production setting (identical to F1's section 2). The objective is S12's own swap convention
(`s12/obj_profile.py:19`), so the comparison to `s12/obj_FINDINGS.md` section 6 is like-for-like:

    target_p = q'_{sep(p)} + ( E[d]_p - q_disto_{sep(p)} )      only the profile is swapped

| arm | `q'` | what it is |
|---|---|---|
| **PROD** | `q_disto` | the shipped distogram profile -- the comparator (section 6's 3.078 baseline; F2 re-derives it rather than quoting it) |
| **POOL** | `q_pool` | the pool's own mean profile: the zero-information typicality arm (section 6: 3.383) |
| **RATIO** (primary) | `q_pool * r_hat` | `r_hat` from a leave-fold-out model on native-free features (below) |
| **RATIO-SHRINK** (the rule-20 control) | `q_pool * (1 + s*z)` | RATIO's realised shrink with a random direction, no information |
| **ORACLE-RATIO** (diagnostic, labelled) | `q_pool * r_true` | the ceiling of this parameterisation; ORACLE, chooses nothing |
| **ORACLE-PROFILE** (diagnostic, labelled) | `q_true` | section 6's 2.402 / section 4's 2.299 reproduced as the instrument check |

**Features for `r_hat` (all native-free, all computable before any structure is chosen):** the
distogram's own profile ratio `q_disto/q_pool` per shell; the pool profile's own shape (its shell
values and their dispersion across the 500 members); the **compactness block** -- the pool's Rg
distribution (mean, sd, skew), the distogram's implied Rg from `E[d]`, and the ratio of the two
(memory: `prediction-pool-disagreement-is-a-native-free-signal`); length; composition and the five
`_PROP_TABLE` means/stds; ESM-2 PCA-32 pooled mean/std (as S12 had). Model: ridge with the penalty
chosen **inside** the training folds by nested CV, never on the 126. Fitted per shell index with
the shell as a feature so short targets are not thrown away.

**Leave-fold-out throughout**, on the 5 pinned folds read from the instrument (`ST.pinned_folds`),
with the NaN-poison test for the deployable arm (contract rule 7): replacing every native quantity
with NaN must leave `r_hat` and the emitted coordinates bit-identical.

## 4. THE FALSIFIERS, IN ORDER, AND WHAT EACH ONE DECIDES

**Stage 1 -- the meter first (contract rule 19).** The RATIO scorer goes through
`s29/s29_D_cost_audit.py` before any endpoint run: ladder Spearman, gradient cosine at production,
the native's percentile in the pool, the ORACLE-structure preference with its pool-member control.
Reported with the rule-20 signature attached. **A scorer that fails meter 3 (the native's
percentile no better than the shipped 0.368) does not proceed to the endpoint**, because that is
what the shell-profile deployable arm already did in S12 (0.459, `:198`).

**Stage 2 -- THE B2 TEST, which is the experiment (primary).** For each target, form the
displacement `u = A_RATIO - A_PROD` between the two arms' emitted point clouds (rigid body removed,
the bound's own construction) and measure `rho = cos(u, t - A_PROD)` against the native, ORACLE,
post hoc. Report the mean, the fold-clustered CI, the per-target sign fraction `q`, and the
signed-value `|rho|(2q-1)` the bound actually consumes (S29-L23, assumption B4).

> **B2 falls only if** the mean `rho` exceeds **0.140** with the fold-clustered CI excluding 0.140
> **AND** the RATIO arm's cosine exceeds its RATIO-SHRINK twin's by more than 1.0x the MDE of that
> paired contrast. Anything less is a confirmation of B2 on this quantity, and will be written as
> one.

**Stage 3 -- the endpoint, only if stage 2 clears.** Built-chain CA-RMSD, n = 126, paired against
PROD with `ST.compare` and the fold-clustered CI; a RESULT requires **0.7x MDE, the fold CI
excluding zero, and the rule-20 signature attached**. The point cloud is reported beside it as the
intermediate.

**Stage 0 -- the probe (contract rule 16).** 12 targets, the pinned `PROBE_TARGETS` set, before
any 126-target run; the probe is never quoted as evidence for the instrument.

## 5. WHAT WOULD MAKE ME REPORT ZERO EVEN WITH A POSITIVE NUMBER

- **The shrink twin matches** (section 2). This is the most likely way a positive appears.
- **`r_hat` collapses to 1.** If the fitted ratio's sd about 1 is a small fraction of the true
  ratio's, the arm IS the POOL arm with extra steps; report the shrink `s` and say so.
- **The cosine is carried by a handful of targets.** Report the median beside the mean and the
  drop-top-10 against a uniform-effect null (memory: `the median-vs-mean gap is the free early
  warning`; and a raw drop-top threshold is not a valid test on its own).
- **The gain is a length effect.** The profile is strongly length-dependent and the 126 span 9 to
  16 residues; a model that only learns length reproduces the pool profile. Registered control: the
  same fit with length as the ONLY feature.

## 6. MULTIPLICITY AND COST

Comparisons: stage 2 contributes **2** (RATIO vs 0.140, RATIO vs its shrink twin); stage 3, if
reached, contributes **2 endpoint comparisons** (RATIO vs PROD on the built chain, and the POOL arm
beside it). The ORACLE arms are diagnostics and are counted as such, never as endpoint comparisons.
The sprint total is the coordinator's to keep (`s29/STATE.md`).

Cost: profiles and pools are cached; the fit is a ridge over 126 rows x ~120 features (seconds);
the expensive part is the projection, ~3 s per target per arm. Probe (12 x 4 deployable arms)
~3 min; full (126 x 4) ~25 min, plus the two ORACLE diagnostics. All through
`python s26/jobrun.py --agent S29M --tag CPU`, checkpointed per target, resumable.

**F1's 126-target run has priority and is queued at the box's job cap; F2 does not start until F1's
rows are complete and its ledger entry is posted.**
