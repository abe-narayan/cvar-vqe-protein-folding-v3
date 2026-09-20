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

---

# ADDENDUM 1 (2026-09-20 02:2x, lane M) -- WHAT SUPPLIES THE PROFILE AT INFERENCE; A SUPPLY AUDIT THAT COSTS NO PROJECTIONS; AND A STOPPING RULE THAT SPENDS THE BOX ONLY IF IT CLEARS

Appended, not edited (contract rule 7). Written before any F2 number exists, in answer to the
coordinator's two bullets (00:5x): *what supplies those n-2 numbers at inference without the native,
and pre-register the falsifier before the run.*

## A1.1 The answer to "what supplies them" is NOT "nothing" -- and that matters, because the fallback closure the coordinator offered would be wrong as written

The coordinator's proposed fallback -- "if nothing native-free supplies them, report the class as an
ORACLE ceiling with no deployable instantiation" -- does not fit this class. **Native-free rules
that produce the shell profile exist, there are at least five of them, and every one is already
measured on this instrument** (`s12/obj_FINDINGS.md:306-331`, all leave-fold-out):

    q' supplied by                                              profile MAE   emitted
    the distogram's own profile        (THE INCUMBENT, native-free)   2.458     3.078
    LFO per-shell debias of it                                        2.458     3.071
    LFO affine recalibration + the pool profile as 2nd regressor      2.394     3.089
    LFO ridge over ESM + composition + disto + pool profile           2.960     3.309
    the pool's own mean profile        (pure typicality, no model)    2.717     3.383
    ------------------------------------------------------------------------------
    the TRUE profile                   (ORACLE)                       0.000     2.402

So the honest statement of the class is stronger and more specific than "no instantiation":

> **The shell profile HAS a native-free supply -- the shipped distogram is one, and it is what
> production already uses. The 0.75 A is not the gap between "an oracle" and "nothing"; it is the
> gap between the true profile and the best of five measured native-free estimates of it, none of
> which beats simply using the distogram's own.**

That reframing is what makes F2 a well-posed question rather than a fishing trip: **F2 is not
asking whether the profile can be supplied. It is asking whether it can be supplied BETTER, and by
how much in the bound's own currency.**

## A1.2 The supply audit: the decisive measurement costs minutes and NO projections

The expensive stage of every arm is the ideal-geometry projection (~3 s per target per arm). **The
cosine that decides B2 does not need it.** Emitting a point cloud is scoring 500 candidates, a
tie-safe top-75 and a 75-member coordinate average -- milliseconds. So the whole of stage 2 (the B2
test) can run on the point cloud for every arm over all 126 targets in minutes, and the built chain
is spent only on an arm that has already cleared.

**Stage 0, registered here, run before anything else:**

**(a) The supply gap, per predictor.** For each native-free candidate `r_hat` (the distogram's own
ratio; the pool profile, i.e. `r = 1`; length-only; the compactness block alone; the full LFO
ridge), measure leave-fold-out against the true ratio `r_true = q_true / q_pool`:
  - `corr(r_hat, r_true)` pooled over all (target, shell) cells and per target;
  - the realised shrink `s` = sd(`r_hat` - 1) / sd(`r_true` - 1) (rule 20's quantity (a));
  - the fraction of `r_true`'s variance the predictor explains, leave-fold-out.

**(b) The cosine ceiling this implies, and then the cosine MEASURED.** The implied ceiling is
reported first as arithmetic (a predictor capturing a fraction `c` of the true ratio's variation
supplies at most about `c * rho_oracle` of the ORACLE displacement's cosine, with
`rho_oracle ~ 0.66` from section 0) -- and then **measured directly**, not inferred: emit each arm's
cloud, form `u = A_arm - A_PROD`, and compute `rho = cos(u, t - A_PROD)` per target, ORACLE and post
hoc, with the fold-clustered CI, the per-target sign fraction `q`, and the signed value
`|rho|(2q-1)` the bound consumes (S29-L23 assumption B4).

**(c) The ORACLE arms in the same table**, as the instrument check: ORACLE-RATIO and ORACLE-PROFILE
must reproduce `s12/obj_FINDINGS.md`'s 2.402 / 2.299 through my code before any deployable number is
read. If they do not, F2 stops and the discrepancy is the finding.

## A1.3 THE STOPPING RULE -- pre-registered, and it is the point of this addendum

> **The deployable projection run does not happen unless stage 0(b) shows the best native-free arm's
> measured displacement cosine exceeding 0.140 with the fold-clustered CI excluding 0.140, AND
> exceeding its own rule-20 shrink twin by more than 1.0x the MDE of that paired contrast.**
>
> If it does not clear, F2 reports -- as its result -- **the ORACLE ceiling (cosine ~0.66, 2.299 A
> point cloud), the measured supply gap that closes it, and the five native-free instantiations that
> already exist and lose**, and spends no further box. That closure is written as a closure, in the
> heading, not as a null of an arm that was never run.

This is the coordinator's instruction implemented with one correction: the closure is *"the supply
is measured and insufficient"*, not *"there is no supply"*.

## A1.4 What I expect, stated before the numbers

The record's own arithmetic makes the direction fairly clear. `s12/obj_FINDINGS.md:306-331` reports
`r ~ 0.4-0.6` between the predicted and true profile, and the best-MAE arm emits worse than the
baseline. Lane L's compactness channel (S29-L19) is +0.909 for the ORACLE Rg and **0.24 to 0.37**
for achievable native-free proxies (`in-band-ordering-is-per-target`). Taking 0.37 as the optimistic
supply and 0.66 as the ORACLE cosine gives about **0.24**, which would clear 0.14 -- and taking the
record's more typical 0.24 gives **0.16**, which barely does, before the shrink twin subtracts its
share. **So the honest pre-run position is that this is genuinely close to the line, which is
exactly why it is worth the minutes and why the shrink twin decides it rather than the raw cosine.**
The registered endpoint prediction (T's, adopted at section 1) is unchanged: inside 0.05 A of
production.

One asymmetry worth stating now: a cosine that clears 0.14 but yields an endpoint inside 0.05 A is
**not a contradiction** -- the bound's identity converts rho to RMSD at the OPTIMAL step, and a
deployable arm takes the step its own score dictates, not the optimal one. If that is what happens,
the reportable quantity is the cosine (B2 falls, the route opens, the step is a separate problem),
and the entry must say so without inflating the Angstroms.

---

# ADDENDUM 2 (2026-09-20 02:3x, lane M) -- THE RULE-20 CITATION IS CORRECTED, THE SHRINK TWIN IS RE-JUSTIFIED ON A FOOTING THAT SURVIVES, AND THE SUPPLY GAP GAINS A SECOND CURRENCY: BITS OF THE READOUT'S 7

Appended, not edited (contract rule 7). Still before any F2 number exists.

## A2.1 Rule 20's cosine justification is withdrawn; my addendum 1 leaned on the withdrawn half

Contract addendum 5 (30), on lane D's S29-L37: the claim "a positive gradient cosine is purchasable
with zero information by shrinking the target map toward typicality" **is refuted on this
instrument** -- across the registered nine-point grid the cosine does not rise, it falls to a
minimum at s = 0.6 and returns, staying inside [-0.056, -0.033] and never approaching zero. **Do
not quote it.** What is confirmed, monotonically and 9 of 9 steps in the predicted direction, is the
other half: shrinking toward typicality **degrades the native's percentile in its own pool**,
0.3688 -> 0.4910.

**My addendum 1, section A1.3, justified the shrink twin by the withdrawn half** ("F2 is pre-loaded
to produce a spurious cosine gain"). That sentence is withdrawn here. The control itself stands, on
two footings that do survive:

1. **It is the matched zero-information control in the operator's own space** -- the project's most
   repeated error (`control-must-match-the-operators-space`, three instances in two sprints). A
   fitted ratio model with realised shrink `s` produces a displacement whose geometry is partly a
   property of `s` alone. The only way to separate "the fit knows something about this sequence"
   from "a map of that shrink displaces the answer that way" is a twin with the same `s` and no
   information. That argument never depended on rule 20.
2. **The percentile axis is confirmed and F2 sits exactly on it.** A shrunk ratio model is a map
   pulled toward the pool's typical profile, which is the move addendum 5 measures as degrading the
   native's percentile from 0.369 to 0.491. S12's deployable shell-profile arm already sits at
   **0.459** (`s12/obj_FINDINGS.md:198`) against the shipped 0.368 (`:191`) -- i.e. **the S12 arm's
   percentile is already 3/4 of the way along addendum 5's shrink grid.** That is a measured
   precedent that F2's arm will do the same thing, and it is now a pre-registered expectation.

**Consequently the F2 report gate is restated:** every F2 arm prints the shrink `s`, **the native's
percentile (the surviving axis)**, the emitted bond and Rg, beside the cosine -- and an arm whose
percentile drifts toward 0.49 while its cosine improves is reporting addendum 5's measured trade
(ladder rho -0.186 -> -0.054 and ORACLE-preference 0.198 -> 0.373 while recognition degrades), not a
new information channel. Four numbers, never one.

## A2.2 The MAE law, arriving from a fourth independent direction -- to be named in the entry

Inside section 6's own table the **incumbent wins**: the distogram's own profile emits **3.078**,
beating all four alternatives, including the arm with the **best profile MAE** (2.394 -> emits
3.089). That is `MAE does not price selected RMSD` (memory; r = 0.19 among achievable priors) at the
level of the *profile*, which is a fourth independent direction after the prior's MAE, the pool's
MAE and the distance-matrix MAE. The F2 entry names it as such, and it is also the reason F2's
primary readout is the cosine and the bits (below) rather than the profile's MAE or correlation.

## A2.3 THE SECOND CURRENCY: how many of the readout's 7 bits the best native-free profile estimate delivers

The coordinator's new architectural number: on the **built chain**, inside the identical top-128
candidate set, the single best member reaches **2.1549 A** while the prefix-average readout caps at
**2.9122 A** -- the readout costs **0.757 A at 4.0x MDE with the candidates held fixed**, and a rule
needs only **7 bits per target** (log2 128) to claim it.

**This falls out of F2 with one extra line and no new machinery**, because every F2 arm already
computes a score over the pool and the instrument already carries `oracle_rr`. Registered now:

    For each target, hold the candidate set FIXED at the production top-128 (the shipped score's
    own prefix -- the coordinator's set, so only the RANKING varies, never the set).
    Let r = the 1-based rank, under the arm's own score, of the ORACLE-best member of that set.
    bits_delivered(arm) = 7 - log2(r)        (7 when it is ranked first, 0 when it is ranked last)

Report, per arm: the mean over 126 targets, the median, the fold-clustered CI, and the **delta
against the shipped score's own value**, which is the supply gap in the readout's currency. The
ORACLE-RATIO and ORACLE-PROFILE arms give the ceiling of the class in the same units.

**Three caveats, registered before the number exists, because this currency is easy to over-read:**

1. **Bits about the BEST member are not bits the deployed readout can spend.** The terminal operator
   consumes the retained set's MEAN, not its best (`operator-consumes-set-mean`: d_out = 1.16 x set
   mean + 0.04 x set best, R2 0.89). So `bits_delivered` is a diagnostic in the coordinator's
   currency and **is not a claim that F2's arm would realise any part of the 0.757 A** through the
   average. An arm can gain bits and move the endpoint by nothing; that is the m-ladder's whole
   history.
2. **It is an ORACLE diagnostic** (it reads `oracle_rr` to identify the best member) and is labelled
   so in every sentence that quotes it. It chooses nothing and tunes nothing.
3. **7 - log2(r) is an order statistic**, so a per-target maximum over arms must be priced with
   `ST.best_of_k_within` (`grid-oracles-are-order-statistics`). Only the pre-specified primary arm
   is read as a result; the others are reported with the null beside them.

If `bits_delivered` for the best native-free profile estimate is at or below the shipped score's own
value, then the supply gap is total in this currency too -- the profile channel delivers **none** of
the 7 bits the readout leaves on the table -- and that sentence, with its number, is the single most
useful thing F2 can hand the report under the coordinator's framing.
