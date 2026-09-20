# SPRINT 30 LEDGER

Entries are `## S30-L<n> -- TITLE (date time, lane)` with the verdict in the heading.
Run `date` in the same command as the append. Corrections are annotated in place with the
original wording left standing.

---

## S30-L0 -- THE CHARTER, TWO READING-LIST CORRECTIONS, AND THE ARITHMETIC THAT SHOULD GOVERN THE SPRINT: **FIXING TEN TARGETS BEATS IMPROVING ALL 126 BY 0.20 Å** (2026-09-20 12:35, coordinator)

### The charter
Saved verbatim as `s30/BRIEF.md`. One hard constraint (CVaR-VQE remains the spine and main
scientific object), one endpoint (mean built-chain Cα RMSD, 126 targets), current production
**3.2105 Å**, primary target < 3.0, ambitious < 2.5. Section 8's twelve leads are explicitly a
**leads register, not a task list**; the charter states that a sprint ignoring ten of them and
finding the mechanism is a success, and one executing all twelve and finding nothing is not.

### Two corrections to the reading list, made by checking rather than assuming
The charter names `s27/REPORT_S29.md` and `s27/LEDGER.md through the S29 close`. Neither is right:
the S29 report is at **`s29/REPORT_S29.md`** (1,619 lines) and `s27/LEDGER.md` is **S28's** ledger
(4,215 lines, ends S28-L50) while S29's is **`s29/LEDGER.md`** (5,823 lines, S29-L0..L57). Both
were read. There is no separate S29 retractions file; S29's 18 withdrawals live in its report
§14(e). Recorded because this project has four instances of prose naming a path that does not
exist, one of which cost a mis-planned lane-week.

### Section 7's meter already exists
`s29/s29_D_cost_audit.py` (767 lines) implements the cost-RMSD meter. S30 **verifies and extends**
rather than rebuilds; lane D owns it and copies it to `s30/s30_D_meter.py` so the S29 artefact
stays untouched.

### THE ARITHMETIC — computed at the open from lane O's S29 rows, before any new experiment

```
production, n = 126      mean 3.2105   median 2.9661
the worst 18 targets     mean 6.2758
the other 108            mean 2.6997

cap the worst 10 at 3.00 A  ->  mean 2.9074  (-0.3031)   BEATS the charter's primary target
cap the worst 18 at 3.00 A  ->  mean 2.7426  (-0.4680)
cap the worst 30 at 3.00 A  ->  mean 2.5778  (-0.6328)
worst 18 all the way to 2.50 ->  mean 2.1334  (-1.0772)

versus improving EVERY ONE of the 126 by 0.20 A  ->  mean 3.0105  (-0.2000)
```

**Fixing ten targets is worth more than improving all 126 by 0.20 Å.** For an endpoint that is a
mean over a distribution with median 2.9661 and a tail reaching 8.24 Å, a mechanism that works on
the typical target and leaves the tail alone is close to worthless, and a mechanism that only works
on hard targets can hit the primary target on its own.

**The caveat, stated with the number so it is not over-read:** capping is an ORACLE operation. It
bounds the prize; it does not deliver it. Two routes realise part of it — (a) a native-free regime
detector, or (b) a method simply better on hard pools *without needing to know they are hard*.
Route (b) requires no detection and is the stronger target.

**And a live, unexplained clue:** the record says the shipped pipeline is **worse than a
sequence-blind one on its 18 hardest targets** (blind 5.425, shipped 6.019) while sequence
conditioning is worth +0.776 Å overall. Something the pipeline does on hard targets is actively
harmful.

### Opening lane assignment (7 of a permitted 8; one slot held for what results demand)

| lane | remit |
|---|---|
| **R** | is nativeness recognizable from a single structure at all (L11) — judged to sit underneath the rest |
| **L** | literature, permanent role |
| **D** | adversary, permanent; owns and extends the cost/RMSD meter |
| **T** | theory: the bit accounting (L8), and when the CVaR tail stops being a prefix |
| **F** | the failure tail — highest leverage on the endpoint by the arithmetic above |
| **X** | divergent, permanent: should the quantum stage select at all, or generate? |
| **Q** | L5 + L6 **together** — a subset objective is pointless through an averaging readout, and a sparse readout is unusable without a way to choose its support |

Contract: `s30/S30_CONTRACT.md`, 28 rules, each annotated with the incident that paid for it.

## S30-L1 -- A CORRECTION TO THE PROJECT'S MOST-CITED NEGATIVE, POSTED BEFORE MY OWN NUMBERS EXIST: S28-L48's "20 OF 31 SCORERS PREFER PRODUCTION TO A 0.25 Å ORACLE STRUCTURE" IS A **CROSS-KIND** PREFERENCE, NOT EVIDENCE THAT NATIVENESS IS INVISIBLE -- EVERY RUNG IN THAT LADDER DIFFERS IN CONSTRUCTION AS WELL AS IN NATIVENESS, AND PERCEPTION-DISTORTION ALREADY PREDICTS THE CROSS-KIND HALF (2026-09-20 12:41, R)

Posted on the coordinator's instruction **before any lane R measurement exists**, so that the
correction cannot be read as motivated by whatever my ladder returns. My pre-registration
`s30/PREREG_S30_R.md` was committed at 12:40:35 (commit `7eabffee`) and contains no aggregate
over more than the single probe target it declares. This entry rests on reading the artefacts and
the definitions, not on new data.

### WHAT S28-L48 MEASURED, FROM ITS OWN DEFINITION

`s27/REPORT_S28.md` section 3.7, the lane C2 audit (entries S28-L35, L36, L37, L48, L49):

> BUILT CHAIN, 31 scorers (16 backbone + 15 CA on the projected chains), one max-over-31 null:
> **20 of 31 prefer the projected production average to a 0.25 Å ORACLE structure** with the fold
> CI below 0.5 (18 informative once the 5 tie-dominated are set aside; DIS on 93% of targets, LEG
> on 79%, RAMA on 59%); 6 coin tosses; CAGEO collapses to anti-recognition on ideal geometry
> (0.611 to 0.421) as registered; the one two-clause pass, CONTACT@chain 0.583, sits at the
> max-over-31 null's mean (0.579, p95 0.627, p_max 0.39).

The verdict attached to it, verbatim: *"the closure claim STANDS over the whole library. No
native-free scorer in the project recognises a near-native structure when it is offered one."*

### THE DEFECT: THE TWO STRUCTURES BEING COMPARED DIFFER IN **KIND**, NOT ONLY IN NATIVENESS

Read from `s29/s29_D_cost_audit.py:250 build_rungs`, which regenerates C2's rungs exactly and
asserts them against the S28 artefacts (`s27/results/s28_A_structs/`, 379 files, verified present):

| rung | how it is CONSTRUCTED | kind |
|---|---|---|
| PROD | `H.readout_uniform(cand, top)` -- the uniform coordinate average of the DIS top-75 | a coordinate average, contracted 25.8% (`averaging-space-beats-the-objective`), then projected |
| circ_best | `A.oracle_circuit_ceiling(...)` best-of-5 -- a quantum circuit's amplitude readout, fitted WITH the native | a circuit output |
| sub0 | `A.oracle_subspace_ls(...)` -- a least-squares fit in a random subspace, fitted WITH the native | a least-squares fit |
| RAND_SIGNED / GAUSS_* | `C2.rand_signed` / `C2.gauss_perturb` of PROD at an ORACLE-matched scale | perturbations of PROD |

So the comparison "does scorer f prefer circ_best to PROD?" varies **two things at once**:
nativeness (0.25 Å vs ~3.2 Å) and construction (a circuit amplitude readout vs a projected
coordinate average). A scorer that is entirely blind to nativeness but sensitive to construction
would produce exactly the reported table. The measurement therefore cannot separate

    (H-kind)    scorers detect WHICH PROCESS produced the structure, from
    (H-native)  scorers cannot see NATIVENESS.

### AND THE RECORD ALREADY CONTAINS THE MECHANISM FOR (H-kind)

S29-L12 imported Blau & Michaeli's perception-distortion theorem and the report states the
consequence itself (§9.1, §12.0): *for any distortion measure, the distortion-optimal estimator's
output distribution must diverge from the real one, so every realism-type scorer must disprefer
the RMSD-optimal answer.* S29-L12's own text names S28-L48 as "that statement measured".

That is the point. **If perception-distortion predicts the result from the construction alone,
then the result carries no additional information about whether nativeness is visible.** A
cross-kind preference is the theorem's content; it is not an independent test of recognition.
The two have been cited as if they were two facts. They are one fact and its explanation.

### WHAT S28-L48 DOES AND DOES NOT ESTABLISH

**It does establish**, and this survives intact and is deployably important: *no scorer in this
library will pick a near-native structure out of a set whose members were produced by different
processes* -- which is the real selection setting, and is why every native-free selector priced in
S12/S27/S28/S29 fails at the endpoint. It also stands as the measured instance of
perception-distortion on this instrument.

**It does not establish** the stronger sentence the record has been leaning on -- *"nativeness is
not recognisable"* -- because it never held construction fixed. The S29 report's own §12.0 points
the same way from the opposite side: **9 of 31 channels keep in-band ordering skill ≥ 0.15** after
partialling out both rank(Rg) and rank(|Rg − median Rg|) (S29-L50), which is not the behaviour of
a library that is blind to structural quality. §12.0 and §9.1 have been read as consistent because
one was about preference and the other about ordering; the kind confound is why they can both be
true.

### CONSEQUENCE, AND WHAT I AM DOING ABOUT IT

1. The sentence "no native-free scorer recognises a near-native structure when it is offered one"
   should be quoted as **"no native-free scorer prefers a near-native structure of a DIFFERENT
   CONSTRUCTION to the production average"**. The qualifier is not cosmetic; it is the whole
   content of the objection.
2. Lane R's ladder is built precisely to remove the confound: every rung is an ideal-geometry
   backbone built from (phi, psi) by the same function, perturbed residues draw their torsions
   from the fold's leakage-safe Ramachandran table, and the perturbation BUDGET is held fixed
   within each stratum. Kind, local realism and budget are matched; only nativeness varies.
3. The instrument's scope, stated now rather than later: its floor is the torsion rebuild at
   **0.284 Å from `nat_ca`** on the probe target, not 0 Å, and every sentence quoting it will say
   so.

**This entry changes no number.** It changes which sentence the numbers support. I am posting it
before my own results because the correction is worth exactly as much as its independence from
them, and because I was handed S28-L48 as settled and it is not.

Artefacts read: `s27/REPORT_S28.md` §3.7; `s27/LEDGER.md` (S28-L35/L36/L37/L48/L49 lines 3833,
3870, 4006, 4033, 4046); `s29/REPORT_S29.md` §9.1, §9.4, §12.0; `s29/LEDGER.md` S29-L12, S29-L50;
`s29/s29_D_cost_audit.py:250`; `s27/results/s28_A_structs/` (379 files, present);
`s29/results/s29_D_ladder_structs/` (126 files, present).

## S30-L2 -- THE CHEAP DIAGNOSIS, AND IT IS DECISIVE: THE TAIL IS **NOT POOL-LIMITED**. THE ORACLE BEST POOL MEMBER ON FAIL18 IS **2.284 Å** -- ALREADY UNDER THE 3.00 Å CAP THE SPRINT-OPEN COUNTERFACTUAL ASKS FOR. THE PIPELINE THEN THROWS IT AWAY: THE DISTOGRAM'S 500 -> 75 FILTER IS **+1.767 Å WORSE THAN A RANDOM 75** ON FAIL18 (3.87× MDE, 0W/18L, RANDOM-18 NULL p = 0) WHILE BEING NO WORSE THAN RANDOM ON THE OTHER 108. **BLOSUM RETRIEVAL IS NOT THE CULPRIT** (NOT MEASURED EVERYWHERE). THE HARMFUL THING THE PIPELINE DOES ON HARD TARGETS IS NOW LOCATED IN ONE STAGE (2026-09-20 12:44, F)

Pre-registration `s30/PREREG_S30_F1.md`, committed **73626721 at 12:39:50**, before the first
aggregate number for this experiment existed. Code `s30/s30_F_stagegap.py`; output
`s30/results/s30_F_stagegap.json`, `s30/results/s30_F_stratum_nulls.json`.

**Every RMSD in this entry is ORACLE except `prod`, and every one is a POINT CLOUD figure**
(`rr` from `s8/generate_univ/<pdb>.npz`). The built-chain counterparts for the three rungs that
have them come from S29 lane O's `s29/results/s29_O_ladder_table.json` and are given in
parentheses; the two bases are never mixed inside a row (contract rule 1).

### 1. The decomposition

The shipped path narrows three times: universe (9.8k-21.5k windows) -> BLOSUM pool (500) ->
score-selected top-75 -> emitted m=75 average.

```
ORACLE best of ...        all 126      FAIL18       other 108
  the universe             1.3134      1.6400        1.2590
  the BLOSUM 500           1.7108      2.2842 (chain 2.2845)   1.6153 (chain 1.6117)
  the score's top-75       2.3062      4.6774 (chain 4.6826)   1.9109 (chain 1.9093)
production (DEPLOYABLE)    3.0483      5.8319 (chain 6.0200)   2.5844 (chain 2.7423)

gap decomposition         all 126      FAIL18       other 108
  G_retr = pool - univ     0.3974      0.6442        0.3562
  G_filt = top75 - pool    0.5953      2.3932        0.2957
  G_read = prod  - top75   0.7422      1.1545        0.6735
  G_tot  = prod  - pool    1.3375      3.5478        0.9691
  share of G_tot due to the filter:     67.5%         30.5%
```

### 2. F1a FIRES -- the tail is not pool-limited

The ORACLE best member of the FAIL18 pools averages **2.2842 Å point cloud (2.2845 built chain)**,
below the 3.00 Å cap; **13 of the 18 have an ORACLE pool member under 3.00 Å** and the worst pool
on the whole tail bottoms out at 3.54 Å. The material to cap the tail is already in the candidate
sets the pipeline is handed. This closes the pool-limited branch and nothing in this lane will
pursue retrieval or generation on its account.

### 3. F1c FIRES, and it is where the damage is -- but only on one of the two filters

Each narrowing stage against a random subset drawn in **its own operator space** (contract rule
6), because best-of-k is an order statistic (`grid-oracles-are-order-statistics`) and the raw drop
from 500 to 75 is expected:

```
arm                                        effect    ×MDE   fold CI            folds  W/L      verdict
BLOSUM-500 vs random-500 of the universe
  all 126                                  -0.0722  -0.87  [-0.1329,-0.0172]   4/5   72W/54L  NOT MEASURED
  FAIL18                                   -0.0034  -0.01  [-0.1609,+0.1130]   1/4    8W/10L  NOT MEASURED
  other 108                                -0.0837  -0.99  [-0.1308,-0.0357]   4/5   64W/44L  NOT MEASURED
score's top-75 vs random-75 of the pool
  all 126                                  +0.2292  +1.17  [+0.0599,+0.3401]   4/5   62W/64L  WORSE (type-M zone)
  **FAIL18**                               **+1.7674  +3.87  [+1.4984,+2.3478]   4/4    0W/18L  WORSE**
  other 108                                -0.0271  -0.24  [-0.0920,+0.0320]   3/5   62W/46L  NOT MEASURED
```

**Retrieval is exonerated at every stratum.** BLOSUM's 500 is indistinguishable from a random 500
of the same universe, and most pointedly so on FAIL18 (0.01× its own MDE). The project's
`sequence-conditioning-hurts-the-failures` memory attributes the harm to the corpus that retrieval
draws from; on the ORACLE *best-member* axis that stage is simply neutral, and the harm is one
stage downstream.

**The filter is worse than chance on every single FAIL18 target** -- 0W/18L, and the median
per-target percentile of the realised best within its own exact random-75 distribution is
**0.99999**: on the median tail target, essentially every random 75-subset of the same pool
contains a better member than the score's chosen 75 does. Random-18 null (20,000 draws, seed
30002): observed +1.7674 against a null mean of +0.2288, CI95 [-0.0849, +0.5832], **p = 0**.

### 4. F1b FIRES on all three registered clauses

1. share of the gap due to the filter, FAIL18 0.6746 vs other-108 0.3051, **ratio 2.21** (bar 2.0);
2. `G_filt` on FAIL18 = 2.3932 against a random-18 null CI95 [0.2690, 0.9880], **p = 0**;
3. **difficulty control**: regressing `G_filt` on production RMSD over the 108 and extrapolating,
   **all 18 of 18** FAIL18 targets sit above the 95% prediction band, mean z = **+4.99**, mean
   residual +1.698. The filter's failure is super-linear in difficulty, not a restatement of it.

### 5. The set mean -- why this reaches the endpoint at all

The terminal operator consumes the set **mean** (`operator-consumes-set-mean`, d_out =
1.16·d_set_mean + 0.04·d_set_best, R² 0.89), so the best-member axis above is not by itself an
endpoint story. Measured on the set mean, where a random-75's expectation is exactly the pool mean:

```
filter effect on the set mean     effect   ×MDE   fold CI            folds  W/L
  all 126                        -0.9026  -4.01  [-1.0032,-0.8174]   5/5   110W/16L
  other 108                      -1.0855  -5.11  [-1.2862,-0.9380]   5/5   100W/8L
  FAIL18                         +0.1947  +0.39  [-0.0280,+0.3544]   3/4    10W/8L
```

Stated to the contract's own bar (rule 2): **the +0.1947 on FAIL18 is NOT MEASURED and I am not
claiming the filter degrades the set mean there.** What *is* measured is that the filter's large,
5/5-fold, 100W/8L set-mean benefit on the other 108 **vanishes on the tail**, and that the
difference between the strata is far outside the random-18 null (observed +0.1947 against null
mean -0.9036, CI95 [-1.2886, -0.5203], **p = 0**). Through the 1.16 coefficient a vanished -1.09 Å
set-mean benefit is worth roughly +1.26 Å of emitted RMSD on those targets, which is the right
order for the record's blind-beats-shipped gap on FAIL18 (5.425 vs 6.019) -- and that gap is now
attributable to a stage rather than to a corpus.

### 6. The circularity, declared, and what survives it

FAIL18 is defined by production RMSD, and production is the average of the set the filter chose,
so a filter that chose badly **mechanically** lands its target in FAIL18. Part of the +1.7674 is
that selection. Two tail definitions that never see the filter or production, scored through the
same random-18 null:

```
tail definition                       overlap    filter vs random-75    null CI95          p
worst 18 by POOL MEAN                  9 / 18        +0.6708          [-0.0805,+0.5782]   0.0148
worst 18 by ORACLE best-in-pool        8 / 18        +0.6320          [-0.0788,+0.5894]   0.0300
(FAIL18, for contrast)                18 / 18        +1.7674          [-0.0849,+0.5832]   0
```

**The effect is real under filter-independent stratification and clears the null both times, but
FAIL18's magnitude is inflated roughly 2.7× by the way FAIL18 is defined.** Anyone quoting
"+1.77 Å" must quote the +0.63/+0.67 beside it. I would rather post this calibration myself than
have the adversary lane find it.

A third, non-circular corroboration: **44 of 126 targets** have a filter worse than 90% of random
75-subsets, and **26 of those are not in FAIL18**. Their production mean is **3.1360** against
**2.4095** for the remaining 82 non-FAIL18 targets. The pathology is broader than the tail and it
tracks emitted quality outside it.

### 7. One declared deviation from the pre-registration

The prereg registered a 2,000-draw sampler for both controls. The min of a random k-subset depends
only on the rank of the smallest index drawn, whose law is exact and closed-form
(P(min rank ≥ i) = C(n−i,k)/C(n,k)), so the exact computation replaced the sampler: **same
estimand, zero Monte-Carlo error, no free parameter to tune**. The registered sampler was run
alongside on the k=75 arm as a cross-check and agrees to a maximum absolute difference of
**0.0258 Å** across all 126 targets. The deviation is annotated in the source docstring.

### 8. What this does and does not license

It does **not** produce a fix, and every number above except `prod` is ORACLE. A selection-limited
verdict says the information is present in the pool, not that any native-free rule can find it, and
S29 closed recognition three independent ways. What it does establish is that **the branch of the
lane that would have gone after retrieval or generation is closed**, and that the sprint-open
clue -- "something the pipeline does on hard targets is actively harmful" -- is now a named stage
with a measured sign, rather than a curiosity.

Artefacts read: `s29/REPORT_S29.md` §7.1/§12/§13; `s29/s29_O_FINDINGS.md`;
`s29/results/s29_O_ladder_table.json` (present, 72,521 bytes); `s12/instrument.py:50` (FAIL18,
pinned, not re-derived); `bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json` (`sub`, verified to
index into the pool, not the universe); `s8/generate_univ/*.npz` (126 present); project memory
`sequence-conditioning-hurts-the-failures`, `operator-consumes-set-mean`,
`grid-oracles-are-order-statistics`.

## S30-L3 -- THE METER REPRODUCES: ALL TEN S29 ANCHORS TO <5e-4 ON A COLD RE-RUN. BUT THE INSTRUMENT COULD NOT MEASURE A NEW COST ON THE **BUILT CHAIN** AT ALL -- THE REPORTING BASIS WAS AVAILABLE ONLY TO THE 31 SCORERS S28 HAD ALREADY STORED. FIXED, PLUS FOUR EXTENSIONS SECTION 7 ASKS FOR (2026-09-20 12:45, D)

**Verdict: the six numbers the brief asked me to check all reproduce; the defect is elsewhere, in
what the meter could not be pointed at.**

### 1. Reproduction (job 1, the part that had to be true)

Cold re-run of `s29/s29_D_cost_audit.py` on the 126 dev targets, nothing cached from S29's run
reused except the ladder structures (themselves re-asserted against the S28 artefacts on load):

| anchor | S29 published | S30 cold re-run | dev |
|---|---|---|---|
| ladder rho, S28 rungs, chain | -0.4023 | **-0.4023** | < 5e-5 |
| ladder rho, S28 rungs, CA | -0.1818 | **-0.1818** | < 5e-5 |
| gradient cosine | -0.0339 | **-0.0339** | < 5e-5 |
| random-direction reference \|cos\| | 0.140 | **0.140** | < 1e-3 |
| native percentile | 0.3688 | **0.3688** (DIS_SURR) / **0.3676** (DIS) | see 2 |
| pref(circ_best vs PROD), chain / CA | 0.0714 / 0.2063 | **0.0714 / 0.2063** | < 5e-5 |
| pool-member control, chain / CA | 0.0201 / 0.1265 | **0.0201 / 0.1265** | < 5e-5 |

Ten checks, all pass, now asserted in code: `python s30/s30_D_meter.py verify`.

### 2. Two things the brief said that the code does not do

(a) **`selftest` does not check any of this.** `s29_D_cost_audit.py selftest` is a synthetic
8-residue check that the RMSD-gradient cosine is +1 and that `spearman`/`pct_in_pool` behave --
it runs in 0.4 s and would pass on a meter whose every dev-set number had drifted. The brief's
instruction "run selftest and confirm it reproduces S29's published baselines" could not have
been satisfied by that command. The reproduction check now exists as a separate `verify` mode.

(b) **The 0.3688 anchor is the SURROGATE's, not the shipped cost's.** `DIS_SURR` (S~, the
linear-interpolation surrogate) gives 0.3688; the shipped lookup `DIS` gives 0.3676. S29's own
lane D recorded both (`s29/s29_D_FINDINGS.md:190`); the brief quotes 0.3688 beside the other
five numbers, which are all `DIS`. Not a defect, but the percentile must always name its cost.

### 3. THE REAL DEFECT: the reporting basis was closed to new costs

The meter has three bases: `ca` (the point cloud), `chain` (the built chain, recomputed) and
`chain-s28rows` (the built chain, read out of C2's stored rows). **`chain` was unusable: 0 of
126 ladder-cache files carried the chain projections it requires**, so any `--basis chain` run
raised `FileNotFoundError`. `chain-s28rows` works only for the 31 scorer names S28 had already
evaluated and stored -- by construction it can never serve a cost written this sprint.

So the instrument that is supposed to gate every new cost function could measure new costs
**only on the CA point cloud**, and the project's own reporting basis is the built chain. That
matters here more than usual: on the shipped cost the two bases do not merely differ in
magnitude, they differ in **sign of the headline ladder** (CHARTER rho **+0.2603** on CA,
**-0.0921** on the chain) and by a factor of 2.2 on the S28 ladder (-0.182 vs -0.402). A cost
cleared on CA alone would have been cleared on the wrong instrument.

Fixed: `s30/s30_D_meter.py build-cache --chain` now builds the projections (47 s/target, four
shards, ~25 min) into `s30/results/s30_D_ladder_structs/`. S29's cache and results are left
alone; the s30 cache is a separate directory (both are gitignored, regenerable).

### 4. Four extensions section 7 asks for, now standard output

- **(E1) per-target distributions** -- min/p10/q1/median/q3/p90/max, IQR, and the five named
  extreme targets, for every diagnostic. A mean is not a distribution.
- **(E2) FAIL18 vs 108 as a fold-clustered CONTRAST** (Welch SE, MDE, cluster-bootstrap CI,
  per-fold cells), replacing S29's two bare means. **With the caveat that cannot be engineered
  away: fold 0 holds NO FAIL18 target** (FAIL18 by pinned fold = 1:6, 2:2, 3:4, 4:6, 0:0), so
  the split's cluster bootstrap draws from 4 non-empty clusters and its CI is wide by
  construction. Any lane quoting a FAIL18/108 split this sprint inherits that.
- **(E4) matched random-signed structures as a first-class control, R = 8 draws** (S29 used
  seed 0 alone). The per-draw spread is printed beside the contrast, so a single-draw control
  can be seen for what it is.
- **(E5) the cosine's aggregate null.** S29 printed mean |cos| = 0.140, which is the magnitude
  of ONE random direction on ONE target. That is **not** the null for a 126-target mean, which
  is ~10x tighter. Both are now printed and labelled; mixing them makes a real effect look like
  noise and is the mirror image of the error this project usually makes.
- plus **(E6)** a pre-stated PASS/BLOCK gate carrying the 0.7x / 1.0x MDE rule, and **(E7)** a
  multiplicity register counting every comparison the meter emits.

### 5. Standing instruction to the other lanes

Route a cost through `python s30/s30_D_meter.py meter --f <name-or-module:function> --basis chain`
**before** it gets substantial VQE compute. A cost whose S28-ladder fold CI lies entirely below
zero is gated BLOCK: lowering it raises RMSD, and the charter forbids spending a day on that.
If the meter cannot see the merits of your cost, say so -- that is a finding about the meter and
I will extend it.

Artefacts: `s30/s30_D_meter.py`, `s30/results/s30_D_meter_*.json`; anchors asserted by
`s30/s30_D_meter.py verify`.

## S30-L4 -- S29'S "THE DEPLOYED SCORE IS AT CHANCE" **SURVIVES AND IS STRENGTHENED**; ITS SUPPORTING PARAGRAPH DOES NOT. "BELOW THE RANDOM BASELINE ON 82 OF 126" IS WHAT A RANDOM RANKING ITSELF DOES (NULL EXPECTS **78.8**, z = +0.60), AND THE MEDIAN GAP OF −0.575 BITS IS −0.416 UNDER THE NULL. **7 − log₂ r IS RIGHT-SKEWED: ITS MEDIAN SITS 0.416 BITS BELOW ITS OWN MEAN, SO ANY MEDIAN OR WIN-COUNT READ AGAINST THE ANALYTIC MEAN READS "WORSE THAN CHANCE" BY CONSTRUCTION** (2026-09-20 12:48, D)

**Verdict: the headline is right, the "and the mean flatters it" paragraph must be withdrawn.
Five statistics, all at chance. The shipped score's ranking is indistinguishable from uniform --
not worse than it.**

### What the report says

`s29/REPORT_S29.md` §0 item 3 and §11:

> The deployed score delivers **1.442 of 7 bits**. A uniform random ranking delivers **1.405**
> (the exact null, 7 − (1/128)·Σ log₂ r). The difference is +0.0366 at **0.10× MDE**...
> **And the mean flatters it.** The *median* paired difference is **−0.575 bits** — the deployed
> score is **below the random baseline on 82 of 126 targets**, and the median rank of the
> ORACLE-best member of its own top-128 is **72 of 128 against a chance median of 64.5**...
> **the median is the honest summary: the shipped score is at chance ... and on the typical
> target slightly worse than chance.**

I reproduced every number in that quote from `s29/results/s29_M_F2_supply_rows.jsonl` (mean
1.4415, median 0.8301, median paired difference −0.5749, 82 of 126, median rank 72.0, and
`ST.compare` returns +0.0366, SE 0.1372, 0.10× MDE, fold CI [−0.225, +0.311], 2/5 folds, to the
digit). So this is a critique of the *inference*, on the same data, not of the arithmetic.

### The defect: three statistics, each compared to the wrong null

1.405 is the null's **MEAN**. It is not the null's median, and it is not the null's win rate.
7 − log₂ r for r uniform on {1..128} is **right-skewed**: its mean is 1.4050 and its **median is
0.9888**. So under the null itself:

- a random ranking scores **below 1.405 on 62.5% of targets** — 78.8 of 126, not 63;
- a random ranking's **median** "paired difference against 1.405" is **−0.4162 bits**, not 0.

Against a simulated uniform-ranking null (20,000 draws of 126 ranks, two-sided):

| statistic | observed | null p2.5 / p50 / p97.5 | two-sided p |
|---|---|---|---|
| mean bits | 1.4415 | 1.180 / 1.402 / 1.650 | **0.75** |
| median bits | 0.8301 | 0.762 / 0.989 / 1.259 | **0.16** |
| targets below 1.405 | 82 | 68 / 79 / 89 | **0.61** |
| median rank | 72.0 | 53.5 / 64.5 / 75.5 | **0.19** |
| mean rank | 66.18 | 58.1 / 64.5 / 70.9 | **0.60** |

and two omnibus tests of the whole rank distribution: **KS D = 0.0685, p = 0.571**;
chi-square on 8 equal rank bins **8.35, p = 0.303**.

Not one of the three "worse than chance" statistics is significant. 82 of 126 is +0.60 binomial
sd from what a random ranking does to itself. The median rank of 72 against "a chance median of
64.5" compares a sample median to a population median without its sampling distribution: the
median of 126 uniform draws has a 95% interval of [53.5, 75.5], and 72 is inside it.

### What is true instead, and it is a cleaner claim

> The deployed score's ranking of production's own top-128 is **statistically indistinguishable
> from a uniform random ranking on every statistic tested** — mean, median, win-count, mean rank,
> median rank, and the full distribution by KS and chi-square. It is at chance. It is **not**
> below chance, and "on the typical target slightly worse than chance" is withdrawn.

This does not weaken §0. It removes an overstatement that a rival group would have found, and it
makes the surviving claim stronger: "at chance" now rests on an omnibus test of the rank
distribution, not on a mean whose CI happens to straddle zero.

### Why this happened, and the standing rule it earns

Three project-memory rules fire at once here:

- `control-must-match-the-operators-space` -- the median and the win-count were read against a
  control computed in the *mean's* space. The project's most repeated error, now instance 4.
- `median-vs-mean-is-the-free-warning` -- the rule says print the median **and the null
  percentiles as one verdict**. The median was printed; the null's own median was not computed.
- `unstated-operators-align-with-your-hypothesis` -- all three mis-comparisons pointed the same
  way, toward "recognition is closed," which is the conclusion the section was arguing for.

**Standing rule for S30, for any lane using bits as a currency:** 7 − log₂ r is a log transform,
so the mean, the median and the win-count have three *different* null values (1.4050, 0.9888,
62.5%). Never compare one of them to another's baseline. Simulate the null for the exact
statistic you are quoting; the cost is four lines and 20,000 draws.

Reproduction: `s29/results/s29_M_F2_supply_rows.jsonl`, arm `PROD`, seed 30.

## S30-L5 -- ASSUMPTION B2's NUMBER SURVIVES; **THE REASON GIVEN FOR IT IS BACKWARDS.** THE PRINTED VERDICT "NO FIELD'S MEAN |cos| CLEARS THE RANDOM REFERENCE 0.140" IS FALSE ON THE QUANTITY IT NAMES -- **21 OF 21 FIELDS CLEAR IT, BY 11 TO 19 SIGMA** -- AND THE SIGNED MEAN IT ACTUALLY TESTS IS BEING READ AGAINST A NULL FOR A DIFFERENT STATISTIC. AGAINST THE RIGHT NULL (+0.0014 ± 0.0144), **11 OF 21 FIELDS HAVE FOLD CIs EXCLUDING ZERO AND THE BEST IS +7.7 SIGMA**. THE FIELDS ARE NOT NOISE. THEY ARE REAL AND WORTH 0.0195 Å BECAUSE **√(1−ρ²) SQUARES THEM** (2026-09-20 12:51, D)

**Verdict: B2 stands at its stated value (best exploitable ρ = 0.1128 ≤ 0.14). Its argument must
be replaced. The corrected framing is more useful to this sprint than the one it replaces.**

### The defect, in the code

`s29/s29_D_fields.py:243,263`:

```python
ref   = np.mean([np.mean(np.abs(r["cos_random_ref"])) for r in rows])   # 0.1398
beats = bool(abs(cmp0["effect"]) - 1.96 * cmp0["se"] > ref)             # cmp0["effect"] = the 126-TARGET SIGNED MEAN
```

`ref` is the expected **magnitude of one random direction on one target**. `cmp0["effect"]` is a
**mean over 126 targets**. They differ in scale by ~√126. Two errors compound:

1. **The label is wrong.** The verdict string says "no field's **mean |cos|** clears the random
   reference 0.140" but the code compares `abs(signed mean)`. The quantity the sentence names,
   `mean_abs`, runs **0.2505 to 0.3249 -- and clears 0.1398 on 21 of 21 fields.**
2. **The null is wrong for either quantity.** From the same stored reference draws
   (`s29_D_fields_rows.jsonl`, 16 per target):
   - null for the 126-target **signed mean**: **+0.0014 ± 0.0144** (|mean| p95 0.0261)
   - null for the 126-target **mean |cos|**: **0.1398 ± 0.0099** (p95 0.1548)

### What the right nulls say

| field | signed mean | z vs signed null | fold CI excl. 0 | mean \|cos\| | z vs \|cos\| null |
|---|---|---|---|---|---|
| CHAN_DISTPOT | +0.1128 | **+7.7** | yes [+0.088,+0.137] | 0.2998 | +16.2 |
| MSET_250 | +0.1118 | **+7.6** | yes [+0.043,+0.175] | 0.2837 | +14.6 |
| MSET_150 | +0.0946 | +6.4 | yes | 0.2872 | +15.0 |
| CHAN_RG_LAW | +0.0933 | +6.4 | yes | 0.2942 | +15.7 |
| CHAN_CONTACT | +0.0933 | +6.4 | yes | 0.2970 | +15.9 |
| ... | | | **11 of 21** | | **21 of 21** |
| MEDOID | −0.0097 | −0.8 | no | 0.2505 | +11.2 |

Eleven of twenty-one fields have fold-clustered CIs excluding zero. Under multiplicity, one would
be expected. `any_beats_reference` is the empty list; the correct count is 11.

**Internal structured control, so this is not an i.i.d.-null artefact.** The honest objection to a
+7.7 σ figure against Gaussian random directions is `zero-information-control-must-be-plausible`:
a Gaussian direction may be a *worse* measure rather than an uninformative one. The survey answers
itself -- three of its own 21 *structured* fields (MEDOID −0.0097, MSET_50 −0.0187, EXPAND −0.0208)
sit at the signed null. So a signed mean of +0.11 is not a generic property of "any real
displacement field"; DISTPOT is doing something those three are not.

### The corrected claim, and why it is the more useful one

- **WITHDRAW:** "no native-free displacement field beats a random direction; B2 survives because
  the class is empty."
- **REPLACE WITH:** *Native-free displacement fields carry a real, strongly significant alignment
  with the direction to the native -- 11 of 21 with fold CIs excluding zero, the best at +7.7 σ.
  B2 survives anyway, because the exploitable signed cosine is only 0.1128 and the achievable gain
  goes as 1 − √(1−ρ²) ≈ ρ²/2. That is **0.0195 Å.** Reaching 3.00 Å needs ρ = 0.358: 3.2× the
  cosine, and **10× the ρ².** The barrier is not that the signal is absent. It is that the
  geometry squares it.*

Every Å figure in §11 is unchanged -- 0.0195 Å, the 2.708 Å sign-oracle ceiling, the B3 residual
bands. What changes is the sentence a reader takes away, and it changes the sprint's search: you
are not looking for the first field with any signal, you are looking for one with **3.2× the
cosine of the best of twenty-one**, and no amount of combining fields whose ρ ≈ 0.1 gets there
(`decorrelated-errors-exist-but-are-unusable`: fusion gain goes as the square of the weaker
channel).

### The forward defect this closes, which is why it is lane D's problem

The same 0.140 constant is printed by the cost-RMSD meter beside every cosine, as
`random-direction |cos| mean 0.140`. **A lane reading its new cost's +0.09 cosine against it would
discard a field that is 6 σ from the signed null.** That is the "if the meter cannot see a cost's
merits, that is a finding about the meter" case, found before any lane hit it. Fixed in
`s30/s30_D_meter.py` extension **(E5)**: both nulls are computed and printed with their names, and
the gate uses the target-mean null.

**Third instance today of `unstated-operators-align-with-your-hypothesis`**: this mis-comparison,
like S30-L4's two, points toward "recognition is closed."

Reproduction: `s29/results/s29_D_fields.json`, `s29_D_fields_rows.jsonl`; nulls from the 16 stored
`cos_random_ref` draws per target.

## S30-L6 -- THE REFERENCE STATE IS A COMPACTNESS TERM, ALGEBRAICALLY: A FIXED-REFERENCE PAIR POTENTIAL PAYS **-0.256 kT PER PAIR (~14 kT OVER A 13-mer)** FOR A PURE 10% CONTRACTION WITH ZERO SHAPE CHANGE, AND A SIZE-MATCHED REFERENCE PAYS **EXACTLY 0.0000**. THIS MAKES S29-L50's LOADINGS **FORCED RATHER THAN MEASURED**, AND EXPLAINS THE FIELD'S 40-50 RESIDUE WALL MECHANISTICALLY: FOR A 13-mer THE DOPE REFERENCE STATE'S ENTIRE SUPPORT IS [0, 15.05 A] AND THE TABLE'S CUTOFF IS 15 A (2026-09-20 12:52, L)
Full note: `s30/lit/L30_1_reference_state.md`. Arithmetic script: `s30/lit/s30_L_refstate.py`
(closed form only; reads no project data). **No measurement; this is a derivation plus published
arithmetic.**

**THE ALGEBRA.** A statistical potential is `u(d) = -kT ln[P_obs(d)/P_ref(d)]`, so a candidate's
score contains the **separable** term `+kT sum_pairs ln P_ref(d_ij)`. For the ideal-gas reference
`P_ref ~ d^2` that term is `2kT sum ln d_ij`: scale the structure by `lambda` and it moves by
exactly `2 kT N_pairs ln lambda`. It is a pure compactness term with no shape content. For DOPE's
reference (Shen & Sali, Protein Sci 15:2507 (2006): uniform points in a ball of radius
`a = sqrt(5/3) Rg`, their own words), the density `f(d;a) = (3d^2/a^3)(1 - (3/4)(d/a) + (d/a)^3/16)`
on `[0, 2a]` has the exact form `(1/a) g(d/a)`, so `ln f = -ln a + ln g(d/a)`: **scale-invariant iff
the reference's size parameter tracks the candidate's own size.**

**THE ARITHMETIC** (Rg from `2.2 n^0.38`, the same law `s27/ham_lib.py` RG_LAW uses):

    n      Rg      a    2a = the reference state's ENTIRE support
    9    5.07   6.55    13.09 A   <- DOPE is tabulated to 15 A: the top 13% has ZERO reference density
   13    5.83   7.53    15.05 A   <- table and reference run out together
  150   14.77  19.07    38.13 A   <- the 15 A cutoff sits comfortably inside

    finite-size correction R(d;a) at 8 A:  0.278 (n=13) vs 0.690 (n=150) -- a factor 2.5, ~0.9 kT/pair
    uniform 10% contraction, 13-mer:  fixed reference -0.256 kT/pair (~14 kT over ~55 Ca pairs)
                                      size-matched reference  0.0000 EXACTLY, by scale invariance

**WHAT IT CORRECTS.** S29-L50 reported LEG_compactness +0.956, RG_LAW +0.921, LEG_solvation
+0.621, LEG +0.606 as a **measured** pattern of Rg loading. The algebra makes them **forced**: any
distance-based potential with a fixed reference contains a separable pure-scale term, so high
|rho(channel, Rg)| is the expected value and not a finding. **The surprise in S29-L50 was never the
loadings; it was that they fail to explain in-band skill (F1a +0.083 against a +0.40 bar), and that
half is untouched.** S29's framing of its own result should be corrected in the S30 report.

It is also the mechanism under the field's length wall (S29-L1, 40-50 residues): at 9-16 aa the
molecule's diameter is comparable to the potential's interaction cutoff, so the reference is
boundary-dominated. DOPE's authors say plainly *"DOPE, like other statistical potentials, is less
accurate for smaller proteins"*, and their own ablation has a FIXED reference sphere (DOPE-24)
performing substantially worse than the size-adaptive one. A separate source states the failure
direction: *"a reference state that is too small results in an erroneous preference for loosely
packed structures"* -- the reference's size parameter **is the knob that sets a potential's
compactness preference**, and in our channels it is untuned.

**OUR CHANNELS, CHECKED IN SOURCE** (`s27/ham_lib.py` docstring lines 24-45): DISTPOT reference is
*"separation-only"*, CONTACT is *"quasi-chemical"*, ENV is a burial ratio -- **none size-corrected**;
CAGEO and RAMA are angle-based and **already scale-free**. That immunity is a cleaner explanation
than coincidence for why LEG_torsion (+0.181, fold CI [+0.053,+0.236]) and CAGEO (+0.216) are among
S29-L50's nine partialled survivors. Qualification I am not hiding: our references are fitted per
target on a same-length universe, so the pool-average size is already right; what is uncorrected is
the **per-candidate** Rg variation -- which is precisely where in-band ranking happens.

**THE CONSTRUCTIVE PROPOSAL AND ITS PRICE.** Build a pair channel with a per-candidate reference
`a_i = sqrt(5/3) Rg_i`: compactness-free **by construction**, not by partialling. Worth having as a
**band statistic** (S29-L19's self-fulfilling-null warning is exactly this, and this use never
passes through the terminal operator). Worth **approximately nothing** as a deployable ranker: it
creates no information, reaches the endpoint only as a cosine (S29 section 5.1), and the averaging
readout spends at most 0.04 of any ranking (section 5.4). The trade is explicit and against us in
one direction -- ANDIS (Yu et al. 2019): *"native recognition and decoy discrimination cannot be
optimized simultaneously with the same parameter sets"*.

**REJECTED IN THIS TOPIC, with reasons**: fine-grained sub-region Ramachandran and
neighbour-dependent Ramachandran (Ting & Dunbrack) -- the conditioner is **sequence**, measured
dead here (memory `phi-carries-no-sequence-signal`: 36.1 deg with full sequence context vs 36.4 deg
sequence-blind); omega as a discriminator -- fixed at 180 deg by the builder; **Ca pseudo-torsion
(theta, tau) potentials** reported to beat DFIRE/dDFIRE/RWPlus -- **already built as CAGEO**, found
by reading `s27/ham_lib.py` before proposing it, and already priced at +0.216; DFIRE's `r^1.61`
exponent -- **fitted** to protein-size spheres, so it carries the artefact rather than fixing it.

**WHERE I COULD BE WRONG**: the -0.256 kT/pair figure is for DOPE's *ball* reference and our
DISTPOT uses an empirical separation-conditioned one (mechanism transfers exactly, magnitude does
not); `Rg = 2.2 n^0.38` is a folded-protein law, so real peptides are less compact, the true `a` is
larger, and section 2's magnitudes are **overstated** while correct in direction and ordering; and
"scale-freeness explains the torsion channels' survival" is an inference I have not tested against
the other seven survivors.

## S30-L7 -- THE COMMON MODE IS NON-IDENTIFIABLE, NOT MERELY INVISIBLE: UNDER "MEMBER = NATIVE + SHARED BIAS + i.i.d. NOISE" THE POOL'S LIKELIHOOD DEPENDS ON `(t, mu)` **ONLY THROUGH `t + mu`**, SO NO POOL-ONLY ESTIMATOR OF THE SHARED BIAS EXISTS AT ANY K. THE THEOREM PERMITS **EXACTLY THREE** ESCAPES AND THIS PROJECT HAS MEASURED ALL THREE: E1 CLOSED, E2 OPEN AT **-0.022 A**, E3 REFUTED FOR A DIFFERENT GENERATOR (PROVENANCE COSINE **0.9432 ABOVE A 0.9330 WITHIN-SOURCE CONTROL**). THE 75-MEMBER POOL CARRIES **m_eff = 1.4 INDEPENDENT MEMBERS** (2026-09-20 12:52, L)
Full note: `s30/lit/L30_2_common_mode.md`. Commissioned by the coordinator after my message of
~12:50. **Derivation plus arithmetic on an existing artefact; no new measurement.**

**WHAT I AM NOT CLAIMING.** The project memory note `pool-error-is-68-percent-common-mode` already
states the invariance in prose: *"a bias shared by every member moves `c` and every `w_k` together
and leaves every within-pool statistic unchanged."* That sentence is not mine. I read the note's
**body**, not its index line (memory `read-the-memory-body-not-the-index-line`), and the body is
what this entry formalises. What is new is the quantifier, the model class, the exhaustive escape
list, and the observation that the record already prices every escape.

**THE TAUTOLOGY I AM AVOIDING**, stated so the real claim is not confused with it: "for any `t'`
the world `(t', {w_k - t'})` reproduces the data" is true and worthless -- it says only that a set
of structures does not name one of them as the answer. The content comes from the **error model**.

**THE THEOREM.** Model class M: `w_k = t + mu + d_k`, `d_k` i.i.d. `~ G` with `E_G[d] = 0`, `t` the
native, `mu` an unknown shared bias. This is the smallest class containing both the i.i.d. model
the ensemble literature assumes (`mu = 0`) and the pool we have. Then the likelihood is
`prod_k G(w_k - (t + mu))`, which depends on `(t, mu)` **only through the sum**. Hence `(t, mu)` and
`(t - c, mu + c)` induce *identical* data distributions for every `c` and every `K`. **`t + mu` is
identifiable (estimated by the pool mean, error `O(K^-1/2)`); `mu` is not, at any K.**

Corollary, and it is the S23 L9 identity read as an identification statement rather than a variance
decomposition: `mean_k|e_k|^2 = |ebar|^2 + mean_k|d_k|^2` = `160.36 + 63.82`, where the first term
is the NOT-IDENTIFIABLE part and the second is the identifiable one.

**THE NUMBER TO QUOTE TO ANY ENSEMBLE PROPOSAL.** Under i.i.d. the average's squared error would be
`mean_k|e_k|^2 / K`. Observed it is `|ebar|^2`. The ratio is the pool's effective independent size:

    m_eff = 224.18 / 160.36 = 1.398        (per-target, from the note's own f = 0.676: 1.479)

> **The deployed 75-member pool carries the statistical content of about 1.4 independent members.**
> Every aggregation, consensus, bagging, voting and ensemble-QA method prices its gain in K. The
> honest conversion factor from any such paper to this instrument is `(1.4 / K_theirs)`.

**THE QUANTIFIER, EXACT** (the coordinator asked for this explicitly and it is the load-bearing
paragraph). The theorem says `mu` is not identifiable **from the pool, under model class M** -- a
statement about one data matrix and one invariance group. It does **NOT** say `mu` is invisible to
any method. The invariance breaks in exactly three places, because those are the only three places
the proof can fail:

  E1  a restriction on `mu` itself -- a prior on the bias's FORM (breaks "mu ranges over R^3n")
  E2  a constraint `t` satisfies that `t + mu` does not (breaks "t ranges over R^3n")
  E3  a second observation whose distribution depends on `(t, mu)` DIFFERENTLY (breaks the factorisation)

Nothing else can work: any estimator built from `{w_k}` alone -- consensus, medoid, typicality,
dispersion, weighting, re-ranking, selection, re-averaging -- is a function of the identifiable part
and is therefore **constant in `mu`**.

**ALL THREE ARE ALREADY PRICED ON THE RECORD, which is what makes the theorem worth having:**

- **E1 -- CLOSED.** The natural form to posit is a scale error, and one exists (averaging contracts
  the backbone 25.8%, Jensen). But the memory note derives the closure: the optimal rescale
  `s* = <c,t>/|c|^2` is **a function of the invisible component and of nothing else**, which is
  Corollary 1 as algebra and which *derives* the mismatched-native placebo result. The governing
  literature is Kennedy & O'Hagan (2001) and **Brynjarsdottir & O'Hagan, Inverse Problems 30:114007
  (2014)**: calibration parameter and discrepancy are not jointly identifiable, and a discrepancy
  term helps **only** given a strongly informative prior on its shape -- with the warning that a
  *wrong* discrepancy prior is worse than none, which is this project's "confidently wrong costs
  2-3x absent" arriving from a second literature.
- **E2 -- OPEN, MEASURED, SMALL.** A hard physical constraint satisfied by real backbones and
  violated by the contracted average is an E2. Already run: **AMBER-relaxing the average at k=30 is
  -0.022 A [-0.036, -0.009], n=126, valid geometry, 5/5 folds same sign**
  (memory `averaging-space-beats-the-objective`). **This is the only escape this project has ever
  obtained a signed fold-consistent result from, and -0.022 A is the number any geometry-repair,
  constraint-projection or physical-validity proposal should be priced against -- not the 68%.**
- **E3 -- OPEN IN PRINCIPLE; THE OBVIOUS INSTANCE IS REFUTED.** E3 needs a source whose bias is not
  the same `mu`. S24 measured: quality-matched **retrieval-free** provenance cosine **0.9432**
  against a **0.9330 within-source control**. Two candidates from completely different sources are
  MORE aligned than two draws from the same source. **`mu` is not a property of where candidates
  come from; it is a property of the prior they are scored against.** Hence the same note's
  "selection aligns the output with WHATEVER prior it is scored against, including a target-blind
  one", and hence `prior-derivative-is-the-only-steep-lever` (-2.15 A/unit) being the only steep
  direction anyone has found.

**WHY LANE O's EXACT ZEROS ARE EXACT.** PCA on the centred pool returns the common-mode
**direction** but is *identically* uninformative about `mu`'s component along it, because centring
removes exactly that component before the decomposition runs. That is why S29-L47's ORACLE global
`eta` for the PC1 family is `+0.0000 exactly`, and two of its four scalars are exact zeros. An exact
zero in a measurement is almost always something being identically zero for a reason; this is it.

**FOR LANE X** (the coordinator asked me to route this; lane X is not directly addressable from my
session). The theorem **justifies** lane X's premise -- selection is a function of the identifiable
part, so its cap is structural, not a failure to find a good selector. Section E3 then **refutes the
obvious remedy**: changing the generator does not change `mu`. **Suggested pre-registered falsifier:
before spending an endpoint run, measure the provenance cosine of generated candidates against the
incumbent pool and require it BELOW the 0.9330 within-source control.** If it lands at ~0.94 like
every previous source, the arm cannot move `mu` and will measure approximately nothing -- for a
reason known in advance rather than after.

**FAMILY-LEVEL REJECTIONS, assumption named per family** (the form the brief asked for): bagging /
consensus QA (Pcons, ModFOLDclust, DAVIS-EMAconsensus) assume i.i.d. mean-zero errors with gain
`~1/K` -- we have `m_eff = 1.4` and `mu != 0`; negative correlation learning needs members being
**trained**; control variates need a **known mean** (ours is the native); multifidelity MC needs
**high-fidelity samples** as an anchor and de-biases nothing; factor models / PCA / ICA recover the
shared direction only up to the standard **sign-and-scale gauge**, and the offset is removed by
centring before they run; blind source separation needs a structural assumption we do not have.
Krogh-Vedelsby / Ueda-Nakano (Brown, Wyatt & Tino 2005 eqs 9-10) assume **nothing** and remain the
correct framing.

**WHERE I COULD BE WRONG**: the theorem is relative to model class M -- if loadings vary with an
observable, something is recoverable, and the per-target `f` spread (10th 0.372, 90th 0.961) is
*consistent with* varying loadings and does not prove M. `m_eff` is arithmetic on two aggregate
numbers from `s23/results/errdecomp.json` via the memory note, not a recomputation from the raw
file. And I have **not** shown E3 is achievable -- only that the obvious E3 is measured to fail, and
where a working one would have to differ.

## S30-L8 -- THE SET-SELECTION PROBLEM IS **COMPUTATIONALLY EASY AND INFORMATIONALLY EXPENSIVE**, AND THE PROJECT HAS BEEN TREATING IT AS THE REVERSE: `V(S) = f(mean_S W)` FACTORS THROUGH THE CENTROID, SO EVERY SUBMODULARITY GUARANTEE IS VOID (THEY REQUIRE **MONOTONE**; S29-L25 STATES V IS NOT) AND FRANK-WOLFE ON THE SIMPLEX IS THE CORRECT CLASSICAL COUNTERPART. PLUS: MY OWN FINITE-SHOT CVaR BIAS HYPOTHESIS, **TESTED AND REFUTED**, IN THE OPPOSITE DIRECTION TO MY GUESS (2026-09-20 12:52, L)
Full note: `s30/lit/L30_3_set_selection_cvar.md`. Check script: `s30/lit/s30_L_cvar_bias.py`.

**THE STRUCTURAL OBSERVATION.** From S29-L25's own definition `V(S) = f(mean_{i in S} W_i)`: **V
depends on S only through one point of R^(3n), the centroid.** That disposes of the reading list:

- **Submodularity guarantees -- VOID.** Nemhauser-Wolsey-Fisher `1-1/e` and Das & Kempe's weakly
  submodular `1-e^(-gamma)` both require **monotone**. S29-L25 states in its own words that V is
  *"neither additive nor monotone"*. **Monotonicity, not submodularity, is what we fail first.**
- **DPPs / facility location / diversity-aware selection -- VOID.** V is constant on
  centroid-equivalence classes, so there is no diversity structure to exploit. (And Abe et al.
  arXiv:2302.00704 is already in the S29 index as the negative result on diversity interventions.)
- **QUBO / Ising -- available but moot.** If `f` were quadratic in the centroid then
  `V(S) = (1/|S|^2) sum_{i,j in S} Q_ij` exactly: the densest-k-subgraph form, NP-hard at fixed `m`,
  poly-time by max-flow at free `m` (Goldberg 1984). A real structural analogy -- **fixed-m is the
  hard formulation** -- but moot, because the continuous relaxation is cheap.
- **The relaxation.** Reachable centroids are the `1/m`-grid on the simplex; by Maurey's empirical
  method any hull point is within `R/sqrt(m)`. With the project's own dispersion
  `sqrt(63.82) = 7.99`: 0.22 A at m=75, 0.86 A at m=5, 1.36 A at m=2 (RMSD units). **Frank-Wolfe on
  the simplex, seconds, is the correct classical counterpart (contract rule 15) and strictly upper
  bounds any circuit on this objective.**

**THE REFRAME, which is the entry's point.** Beside the project's own numbers -- 2 members with
ORACLE **weights** 1.4315 A vs 75 members with ORACLE **membership** 2.3055 A; choosing 2 of 500 is
~17.9 bits vs 7 bits for the top-128 argmin (S29 section 12.1) -- the expressiveness lives in the
**weights**, and the weighted object costs **more** bits, not fewer. Hence:

> **Solving the set problem better is worth nothing -- exhaustive search already solved it
> (S29-L25) and the answer was worse than production. What is scarce is the information needed to
> SPECIFY a good set, and no solver supplies information.**

That also explains S29-L25's second clause mechanically: escaping the set-equality theorem creates
no information, so the optimum of a marginal-class objective moves the support without moving the
answer.

**CVaR ON A SET FUNCTION -- the literature exists and is sharp.** Maehara (Oper Res Lett 43:526,
2015): the CVaR of a stochastic submodular set function **is not submodular**, and **no
polynomial-time multiplicative approximation exists** unless P=NP. Wilder (AAAI 2018) and Ohsaka &
Yoshida (2017): the escape is to stop asking for a single set and relax to a **portfolio -- a
distribution over sets** -- where a `1-1/e` guarantee returns via continuous DR-submodular
maximisation. **Mismatch stated rather than hidden:** their CVaR is over *exogenous* randomness;
ours is over a distribution the optimiser controls, so the theorems do not transfer. What transfers
is the design lesson, identical to the one above from the other side: **relax the set, do not search
it harder.**

**MY OWN HYPOTHESIS, TESTED AND REFUTED.** I proposed that finite-shot CVaR estimator bias gives the
VQE an incentive to shape its distribution for the estimator rather than for the structure -- a
candidate mechanism for "the objective improves and the structure does not move" and for S29's
"realised m wanders". Checked at the deployed cell (n=9 -> 512 states, alpha 0.18, shots 2048 from
`core/quantum.py:179`, tail k=369):

    A  bias vs concentration (Boltzmann states)     +0.0000 .. +0.0025   (sd 0.017-0.036)
    B  states MATCHED on true CVaR, support 2..32   +0.0001 .. +0.0019, NO trend with support
    C  does it move the argmin?                     exact beta 2.5, finite-shot beta 2.5, penalty +0.0000
    D  Barkoutsos's flat minimiser set              support 92->1: bias -0.0006 -> +0.0162, i.e.
                                                    concentration PENALISED, and by LESS THAN ONE sd
                                                    of the shot noise (0.027)

**REFUTED, and in the direction opposite to my guess** -- the empirical lower-tail CVaR is
*pessimistically* biased (verified analytically: `E[min(X1,X2)] = -0.564` vs `E[X|X<=med] = -0.798`
for a standard normal). At 369 tail shots the bias is ~0.1% of the objective's range, does not order
states by concentration, and does not move the argmin. **Finite-shot CVaR estimator bias is not a
mechanism for anything at this project's shot count. Do not spend on it.** Scope: the bias grows as
`alpha*shots` falls; at the `foldvqe` sort path (384 shots, 69 tail shots) it is ~5x larger and
still an order of magnitude below the effect sizes here. **If any lane proposes a low-alpha or
low-shot arm, this closure lapses and must be re-run.**

**WHERE I COULD BE WRONG**: Maurey's bound is for **multisets** (weights on a `1/m` grid), not the
uniform-weight **subsets** the readout uses -- the subset family is strictly poorer and the
`1.4315` vs `2.3055` gap is the size of that difference, so "the relaxation is tight" is a claim
about the *weighted* problem and a reader could wrongly take it as one about subsets. `R` should be
a **max** over members and I used the RMS dispersion, so the table is optimistic by perhaps 2-3x
(ordering across `m` unaffected). And section 4's refutation is on a **synthetic Gaussian
spectrum**, not a real cell's energies; I judged the margin large enough not to spend a lane's time
re-running it on a real one, but that is my judgement.

## S30-L9 -- A CORRECTION TO MY PREDECESSOR AND TO THE BRIEF: **THE CVaR TAIL NEVER STOPS BEING A PREFIX.** IT IS ALWAYS THE PREFIX OF THE ORDER INDUCED BY `grad V` AT THE OPTIMUM, SO THE QUESTION "WHEN DOES IT STOP" IS MIS-POSED -- AND S29's §4.3 (KEEP THE ORDER EXOGENOUS, FOR WELL-POSEDNESS) **PROVABLY GUARANTEES** WHAT S29's §4.5 CLAIMS TO ESCAPE: UNDER THE TAIL-THEN-AGGREGATE LIFT AS SPECIFIED, THE EMITTED SET IS STILL `argsort(E)[:m]` AND THE ENDPOINT CHANNEL IS STILL THE SINGLE INTEGER `m` (2026-09-20 12:53, T)

**Verdict: S29's section 4 conclusion "this is the first formulation in the project's history whose
classical counterpart genuinely goes away" is WITHDRAWN as stated. What S29 proved is that a FREE
subset optimum is non-prefix. It did not show the lifted CVaR readout can REACH it, and §4.3's own
well-posedness recommendation is exactly the condition that prevents it.** The repair is in this
entry and it is constructive, not merely a withdrawal.

### 1. THE THEOREM (T1), WHICH REPLACES THE QUESTION

Let `Lambda(p) = { lambda : 0 <= lambda <= p, 1'lambda = alpha }` be the tail polytope -- the exact
feasible set of the deployed readout (`core/quantum.py:290 cvar_from_probs` lines 306-308 compute
its vertex). Let `V` be any tail objective differentiable on it.

> **T1.** At every KKT point `lambda*` of `min_{Lambda(p)} V` there is a scalar `mu` with
>
>     grad V(lambda*)_x < mu  =>  lambda*_x = p_x ,       grad V(lambda*)_x > mu  =>  lambda*_x = 0 .
>
> **Hence `lambda*` is a prefix of the order induced by `grad V(lambda*)`, for every `V`.**

*Proof.* `Lambda(p)` is a box intersected with one equality. The Lagrangian is
`V(lambda) - mu(1'lambda - alpha) - <s, lambda> + <t, lambda - p>` with `s, t >= 0`; stationarity
gives `grad V = mu + s - t`, complementary slackness kills `s` where `lambda > 0` and `t` where
`lambda < p`. The three cases are the display. ∎

The deployed CVaR is the case `V = <E, lambda>`, `grad V = E` -- **constant in `lambda`**. That, and
not anything about CVaR, is why one classical sort reproduces it. `cvar_from_probs`'s cumulative
scan is the greedy solution of a fractional knapsack, and the rearrangement inequality is the whole
content of `s25/QUANTUM.md` section 5.

**So the charter's question ("under what conditions does the tail stop being a prefix?") has the
answer NEVER, and the load-bearing question is a different one:**

> **the tail stops being reproducible by ONE CLASSICAL SORT iff the ordering map
> `lambda -> grad V(lambda)` has more than one fixed point.**

That is a statement about self-consistency, not about CVaR, and it is checkable.

### 2. THE INCOMPATIBILITY IN S29's SECTION 4, NAMED IN ITS OWN TEXT

- `s29/THEORY.md:4.1` -- (4.1) "keeps a scalar order to FORM the tail and lets the structure be
  aggregated inside it."
- `s29/THEORY.md:4.2` -- "With the tail formed by a per-state scalar order (energies `E`, quantile
  `q`, boundary state `x_q`)..."
- `s29/THEORY.md:4.3` -- "**Recommendation: keep the order fixed by a per-state scalar** and let
  only the aggregate be structural."
- `s29/THEORY.md:4.5` -- "**the set-equality theorem fails for (b)**"; "Under (b) the tail's
  composition -- which states share the mass -- changes the objective, so 'which set' is a real
  optimisation variable rather than a read-out of the sort."

The first three fix `grad`(tail-formation)` = E`. By T1 with an exogenous order the support of
`lambda*` is the `E`-prefix reaching `alpha`, full stop: the optimiser's freedom is `m`, plus
deletions at exact zeros, which is *verbatim* the deployed situation (`s25/QUANTUM.md` §5, and
`DATAPATH.md` stage 9: "the whole quantum stage reduces, for the structure, to choosing `m`").

§4.5's three-state witness (`W_1=(1,0), W_2=(-1,0), W_3=(0,0.1)`, optimal 2-subset `{1,2}`) is
correct **about the free subset problem** and is *unreachable* by the readout it was written for:
with `l_3 < l_1 = l_2`, state 3 is in every `alpha`-prefix, so `{1,2}` requires `p_3 = 0` exactly --
the "hole" mechanism, which `s25/QUANTUM.md` §5 already showed is chosen "by an amplitude pattern
that is a function of the **same energies** the classical rank order already uses."

**What survives of §4.5 intact:** `V(S) = f(mean_{i in S} W_i)` is neither additive nor monotone;
the free-subset optimum is non-prefix (lane B measured 11/12 pairs and 12/12 five-subsets, S29-L25);
and the gap is not measure-zero. **What does not survive:** that the lifted CVaR stage reaches it.
S29-L25's own second clause -- "the escape buys nothing: those subsets' averages are worse than the
prefix on 10 of 12" -- was measured on the *free* optimum, i.e. on sets the readout cannot emit, so
it is not evidence about the lift either way.

### 3. THE REPAIR, AND ITS PRICE (T1b)

Make the order **endogenous**: `V(lambda) = f(R_lambda)`, `R_lambda = (1/alpha) sum_x lambda_x W_x`,
so by T1 the induced cost is

    c_x(lambda) = (1/alpha) < grad f(R_lambda), W_x > ,

i.e. **the candidates sorted by their projection onto ONE self-consistently determined direction
`g = grad f(R)`.** Now the readout genuinely has `d` degrees of freedom instead of 1. The price:

> **T1b.** Under an endogenous linear order the reachable tails are contained in the halfspace-cut
> class `{ {x : <g, Psi_x> <= tau} }`, VC dimension `d+1`. By Sauer-Shelah the class carries at most
> `log2 sum_{i<=d+1} C(500, i)` bits.

**Measured, 126 real pools of 500** (`s30/s30_T_bits.py`, `s30_T_spec.py`; the rank-collapse
pre-check the brief required *before* compute): stable rank of the centred pair-distance matrix
**1.859** (median 1.865, min 1.247, max 2.721); PC1 alone is **55.4%** of the variance;
`k90 = 5.61`, `k99 = 21.13` directions; for the coordinate matrix stable rank 3.404, `k90 = 11.5`.
Caps against the **300.6 bits** of free 75-subset choice:

    d_eff = 2  (stable rank)       halfspace cap    24.3 bits
    d_eff = 6  (k90)               halfspace cap    50.4 bits
    d_eff = 21 (k99)               halfspace cap   126.7 bits
    d_eff = 33 (nominal 3n-6)      halfspace cap   175.5 bits

**My registered decision rule (stable rank < 2.0 -> CLOSED at the encoding level) fires.** The
tail-then-aggregate lift, even with the order made endogenous, is a 2-to-6-parameter continuous
family, not a combinatorial search: a 250-to-275-bit collapse below the set choice it was supposed
to open. **It should not be built**, and that verdict cost one 76-second pass over data that was
already on disk rather than a lane-week.

### 4. THE ONE SURVIVING ESCAPE, WHICH IS S29's OWN OPEN POINTER

Add a **second moment**: `V = f(R_lambda) + mu*h(Sigma_lambda)` with
`Sigma_lambda = (1/alpha) sum lambda_x (W_x - R)(W_x - R)'`. Then
`grad V_x` is **quadratic** in `W_x`, the cut is a quadric rather than a halfspace, and the VC
dimension goes `d+1 -> (d+2)(d+1)/2`: at `d_eff = 6`, from 7 to 28, and the cap from 50.4 to
**152.1 bits**. A 3x enrichment in the exponent, arrived at from the readout's algebra.

**It is the same object S29's own post-mortem pointed at from the opposite direction** -- "the one
argument for operators that read **the pool's own dispersion** rather than the posterior's
marginals" (THEORY_SUMMARY §2, from `cov(a,n) > cov(b,n)`). Two independent derivations landing on
one class is the only reason this row is worth a measurement rather than an opinion.

**And its hard cap, stated now so no one has to discover it later:** `tr(Sigma_lambda)` is invariant
to a common shift of the pool, so by the S23-L9 identity (`mean_k |e_k|^2 = |ebar|^2 + mean_k
|d_k|^2`, exact to 2.7e-14) **a dispersion term cannot see the 68% common-mode error at all.** It
operates on the idiosyncratic 32% only. Its ceiling is measured in S30-L9 (M4) and not asserted here.

### 5. CVaR's RESIDUAL ROLE, WHICH IS EXACTLY ONE THING

Under an endogenous order the pair `(g, tau)` defines the cut. `g` comes from `grad f`; **`alpha`
supplies `tau`.** That is CVaR's whole remaining function: it turns a direction into a set, with a
budget rather than a threshold, which is what makes the map continuous in `p`. It is not
dispensable and it is not doing selection.

### 6. THE TRICHOTOMY, FOR THE REPORT

| order | tail | classical counterpart | endpoint channels |
|---|---|---|---|
| **exogenous** (deployed, and S29 §4.3 as recommended) | `argsort(E)[:m]` | one sort | **1** (the integer `m`) |
| **endogenous, `f` convex** | halfspace cut, unique up to ties | Frank-Wolfe, `O(1/eps)` linear-minimisation oracles, each a sort | `d_eff` ~ 2-6 |
| **endogenous, `f` non-convex** | halfspace cut, multiple fixed points | greedy + local search over `<= 175` bits | `d_eff`, genuinely combinatorial |

Row 2 is classically polynomial, so it is not a quantum opening; row 3 is the only one that is, and
T1b prices its search space at 24-175 bits against 300.6.

Reproduction: `s30/PREREG_S30_T.md` (registered before any array was loaded, commit 7f8bbce5),
`s30/s30_T_bits.py`, `s30/s30_T_spec.py`, `s30/results/s30_T_bits.json`, `s30_T_spec.json`.
Derivations in `s30/THEORY.md` sections 2-4.


## S30-L10 -- THE GATE LANE X WAS HANDED ("MEASURE THE CEILING FIRST, ALWAYS") IS AN **ANTI-PREDICTOR**: IT GETS THE DIRECTION OF THE ENDPOINT MOVE RIGHT ON **1 OF 10 CELLS** ACROSS TWO INDEPENDENT FAMILIES, WHERE THE FROZEN SOURCE LAW GETS **5 OF 6** AND 22x THE ACCURACY. THE ADMISSION CONDITION THAT REPLACES IT SAYS **GENERATION MUST RAISE A SPACE'S FLOOR, NOT ITS CEILING** -- AND ITS SPECIFICATION (-0.2 Å IN THE TYPICAL MEMBER, EVEN GRANTED A PERFECT CEILING) IS MET BY NOTHING IN THE RECORD, SO **DO NOT BUILD A GENERATIVE SPACE THIS SPRINT**. SEPARATELY, MY OWN INFORMATION-CEILING HYPOTHESIS IS DEAD WITH THE SIGN REVERSED: THE FAIL18 TAIL HAS **TIGHTER** DEPOSITED NMR ENSEMBLES THAN THE OTHER 108 (-0.487 Å), SO THE TAIL IS NOT CONFORMATIONAL AMBIGUITY IN THE REFERENCE (2026-09-20 12:53, X)

Pre-registration `s30/PREREG_S30_X.md` @ **d4305d17**, committed before any number below existed.

> **Numbering note (contract rule 15).** Written as S30-L5 and appended at 12:53. While it was
> being written lanes D, L and T landed S30-L5 through S30-L9, so it was renumbered twice on
> append -- L5 -> L6 -> **L10** -- before anyone read it. No number is skipped, no other
> lane's text was touched, and the chronology is 12:53 regardless of the number.
Findings `s30/s30_X_FINDINGS.md`. Code `s30/s30_X_sourcelaw.py`, `s30/s30_X_ensemble.py`.
Results `s30/results/s30_X_sourcelaw.json`, `s30/results/s30_X_ensemble.json`.
**No VQE or pipeline compute was spent. Every number is from a stored artefact or a PDB parse.**

### 0. What I was sent to do, and what I am reporting instead

My remit: attack "retrieve, then select or average", and ask whether the quantum stage should be
GENERATING. I was given one test to apply before spending compute -- *"does this space contain
structures better than the pool's own ORACLE best? ... Measure the ceiling first, always."*

**I am reporting that the test is wrong-signed, and that the condition replacing it forbids the
thing I was sent to build.**

### 1. THE SCORECARD. Both families are ORACLE; both were already in the record.

| family | cells | **ceiling gate** sign correct | source law sign correct |
|---|---|---|---|
| S7-6 K-ladder, K = 25..2000 vs K = 500 | 6 | **1 / 6** | **5 / 6** |
| S24 generated-source unions, 4 samplers | 4 | **0 / 4** | n/a (set mean not stored) |
| **total** | **10** | **1 / 10** | 5 / 6 |

`P(X <= 1 | n = 10, p = 0.5) = 11/1024 = 0.0107`, **indicative only** -- the six rungs share one
reference and the four arms share 126 targets and one pool, so the cells are not independent. The
magnitude is the finding: on the same six cells, mean |residual| **0.333 Å for the ceiling gate
against 0.0153 Å for the law**, a factor of **22**.

The underlying table (`docs/FINDINGS.md` §S7-6, `s7/poolsize_kcurve.json`, n = 126 per rung,
verified in place before transcription) is the whole argument in seven lines:

```
     K      selected   pool best (ORACLE)   pool mean
     25      3.439           2.350            4.177
     50      3.399           2.177            4.232
    100      3.425           1.970            4.316
    250      3.461           1.808            4.393
    500      3.454           1.711            4.453      <- shipped
   1000      3.483           1.568            4.522
   2000      3.520           1.504            4.596

   ceiling improves 0.846 A  ->  endpoint gets 0.081 A WORSE
```

S7 stated the direction in 2026 prose ("selection tracks the pool MEAN not its BEST ...
*anti*-correlated with the pool BEST"). It was never turned into a gate, and the project has since
run at least two generative probes whose admission test was the ceiling -- S24 lane C and
**S29-L56, which abandoned the chimera space on exactly this number (+0.6533 Å)**.

### 2. H-X1 PRIMARY: PASS, against a registered residual bar, coefficients FROZEN

`d_out = 0.803*d_set_mean + 0.298*d_set_best` -- the Sprint 20 wide-set fit (2,142 (arm,target)
cells, 17 generator arms at 8,192 draws). **Nothing fitted.** Six rungs predicted vs K = 500:

```
    mean |residual|  0.0153 A     registered bar 0.10 A      PASS
    max  |residual|  0.0263 A
    sign correct     5 / 6        registered bar < 1/3 wrong PASS

    matched controls in the operator's own space (contract rule 6):
       two-term source law    |resid| 0.0153    sign 5/6
       mean-only              |resid| 0.1198    sign 5/6
       best-only (THE GATE)   |resid| 0.3330    sign 1/6
```

**DISCLOSURE (rule 14).** One cell -- K25 -> K2000, predicted +0.0843 vs observed +0.0810 -- I
computed by hand *before* writing the prereg. It motivated the hypothesis, is excluded from the
bar, and is carried separately in the artefact as `disclosed_cell_NOT_out_of_sample`.

**IDENTIFIABILITY, registered in advance and confirmed: r(pool mean, pool best) = -0.9898** across
the ladder. This family cannot identify the two coefficients separately, so none were fitted on it.

H-X1 SECONDARY, S24 lane C, **504 (arm,target) cells**, a different sprint / corpus / sampler set:

```
    WITHIN-TARGET slope on the selected set's mean   0.9349
      target-clustered CI95 [0.8201, 1.0632]         R2 0.668
    pooled cross-target (contaminated by difficulty) 1.0413   R2 0.929
```

### 3. THE ADMISSION CONDITION, AND ITS SPECIFICATION

```
    ADMIT a generative space G against pool P iff
        Delta(set mean) < -(b/a) * Delta(set best),      b/a in [0.319, 0.371]
```

> **An X Å improvement in a space's CEILING is worth having only if its TYPICAL member degrades by
> less than ~0.32-0.37 X. At fixed ceiling, the mean is worth ~2.7-3.1x what the best is worth.**

It is also the *cheaper* gate: a set mean is a few draws with small variance; a set best over 10^4
configurations is an order statistic needing `best_of_k_within` and split-half transfer before it
can be read at all.

**Specification.** Best ceiling anywhere on this corpus is the whole-library best **1.313 Å**
(§S8-4, `s8/triage.py`) vs the shipped **1.7108** (ORACLE both). Granting that *in full*,
`0.298 * -0.398 = -0.119 Å`. To reach the sprint's -0.30 Å primary:

```
    Delta(set mean) = -0.19 A (a = 0.935)  to  -0.23 A (a = 0.803)
```

> **A winning generative space must make its TYPICAL member ~0.2 Å better than the pool's, even
> after being granted a perfect ceiling. Generation must raise the FLOOR.** S24's four samplers
> stand at 3.789 / 3.243 / 3.207 / 3.175 against the 3.0483 incumbent -- every one worse in the
> mean while three beat the pool at the ceiling.

**AND THE COROLLARY FOR THE QUANTUM ENCODING.** A `2**q` register over structural variables is by
construction a **wide** space (S29-L56's chimera was 8**S = 32,768). Width is the encoding's
selling point -- and at rho ~ 0 width buys ceiling and costs mean. **The property a quantum
encoding is chosen FOR is the property this condition penalises.** The coefficients are a property
of the **readout**, not of physics: under argmin they would be (0, 1) and the ceiling *would* be
the endpoint. So wide space + averaging consumes the wrong statistic; wide space + argmin needs
rho, which is closed; **narrow, uniformly-good space + averaging needs no rho at all and is the
only one of the three still open.** Generation is closed *jointly with the readout*.

**SCOPE, so it is not misused.** Validated for changes to the candidate **SOURCE at a fixed
score**. NOT validated for changes to the **SCORE at a fixed source** -- S18's amendment to
`operator-consumes-set-mean` measured this law failing there (Legacy gate: premise held, law
predicted a gain, output +0.076 worse). **Lane F's S30-L2 is a GATE change and this condition must
not be composed with it.** I flag that because the two results are adjacent enough to invite it.

### 4. A GAP IN A CLOSURE THAT HAS BEEN DRIVING LANE DESIGN FOR SIX SPRINTS

`prior-derivative-is-the-only-steep-lever` says "do not build another candidate generator -- closed
on five instruments". Checking the artefacts: **four of the five measured only the ENDPOINT after
selection and carry no oracle arm at all** (`biasalign.json`, `qmatch.json` field lists verified).
The fifth (S24 lane C) measured the generated source's oracle best in isolation but **never formed
`min(pool, generated)`**. Recomputed here from its own stored rows, ORACLE:

| arm | of 126 targets, # carrying a structure better than the WHOLE K=500 pool's best | union ceiling (500) | Δ | union ceiling (2000) | Δ | **endpoint Δ** |
|---|---|---|---|---|---|---|
| T0_helix | 12 | 1.6969 | -0.0139 | 1.6814 | -0.0294 | **+0.0085** |
| T1_blind | 39 | 1.6088 | -0.1020 | 1.5297 | -0.1812 | **+0.0170** |
| T2_restype | **45** | 1.5946 | **-0.1162** | 1.5154 | -0.1955 | **+0.0157** |
| T3_pool | **47** | 1.6011 | -0.1098 | 1.5229 | -0.1879 | **+0.0509** |

An **untrained per-residue-type Ramachandran sampler** puts a structure better than the entire
retrieved pool's best member on **45 of 126 targets**, and the endpoint still gets worse.

**This does not weaken S24's verdict -- it is section 1 a third time.** What it changes is the
wording. The closure should read *"no achievable source produces candidates the shipped score can
convert into Ångströms"*, not *"no alternative source contains better structures"*, which is false
on 45/126 targets. The stronger wording has been used to rule proposals out.

### 5. H-X2: MY OWN HYPOTHESIS, DEAD, WITH THE SIGN REVERSED

`core/geometry.py:808` scores against **model 1** (`model_index=0`); `ca_rmsd_to_ensemble` and
`native_ensemble_from_pdb` exist in two modules and are **called by nothing**. 111 of 126 tuning
targets are multi-model depositions (mean 15.6 models, median 20). So the question was real: if the
reference is one draw from an ensemble, RMSD to it has a floor no information can beat.

Floor = `RMSD(ensemble medoid, model 1)`, ORACLE by construction:

```
    floor  mean 0.6136  median 0.3784  sd 0.7630  max 4.2943   21.4% above 1.0 A
    pairwise ensemble spread  mean 0.9535
    endpoint (built chain) 3.2148    FAIL18 6.2872    other 108 2.7027
    corr(floor, endpoint) +0.1220    corr(spread, endpoint) +0.1009

    FAIL18 minus the other 108:
       d endpoint  +3.5845
       d floor     -0.1776     share of tail excess  -5.0%
       d spread    -0.4872     share                -13.6%
```

**VERDICT: PRICED AND DEAD at this endpoint** (registered rule: dead if mean floor < 1.0 Å and it
explains < 25% of the FAIL18 excess). My registered prediction was ~0.6 Å and < 5%; measured
0.614 Å and -5.0%.

**It answers the coordinator's open question #1 directly.** *"Does the tail's difficulty have a
cause we can name?"* -- **it is NOT conformational ambiguity in the reference. The hard 18 have
*tighter* deposited ensembles than the easy 108.** With lane F's S30-L2 (the tail is not
pool-limited either) two candidate causes are now excluded.

**Calibration:** the model-1 convention costs <= ~0.61 Å and does not bind until the endpoint is
below ~1.5 Å. `corr(floor, pool_best_ORACLE) = +0.311` -- a floppy peptide has a floppy retrieval
neighbourhood, which is the mechanism to expect if this is ever revisited.

### 6. A FRAMING I DECLINE TO ACCEPT AS BINDING

My brief cites *"the pool's error is 68% common-mode, 50.7x the i.i.d. prediction"* as the cap that
motivates generation. Reading the memory's **body** (`read-the-memory-body-not-the-index-line`):
`|ebar|^2 = n*RMSD^2` **exactly**, so the "common-mode" numerator *is* the output RMSD -- not
independent evidence about it -- and the same file's own S24 correction says ***"f is NOT a screen
and must never be used as one"*** (a blind library draw scores f = 0.4896 against the incumbent's
0.6758 while being 0.76 Å **worse**). The 50.7x compares an identity against an i.i.d.-member model
that retrieval exists to violate. The durable content is the weaker claim -- pool members resemble
each other more than any resembles the native -- and it is not a cap on generation.

### 7. THE DIVERGENT VERDICT, WHICH IS NOT THE ONE I WAS SENT TO FIND

At rho ~ 0 the endpoint is a near-unit-slope function of the candidate set's MEAN
(0.9349 [0.8201, 1.0632] within target, 504 cells). **The only stage in the pipeline that raises
that mean is retrieval** -- it turns a ~18,674-window universe into a K = 500 set of mean 4.4533 Å.
No downstream operator does it and no generative sampler in the record does it.

> **"Retrieve, then average" is not an accident of history to be escaped. It is the only
> architecture in the record that optimises the one statistic the terminal consumes.**

**RECOMMENDATION, against my own remit: do not build a generative structural space this sprint.**
The specification is -0.2 Å in the typical member and nothing approaches it; and the gate that
would have admitted a proposal is wrong on the sign 9 times in 10.

**Where the live freedom is:** in `a` and `b` themselves, which are readout properties, not
physics. That is where S30's STATE already points ("a sparse weighted readout -- the only ladder
class not closed by ceiling"), reached from the opposite direction and now with a transfer
function attached to it.

### 8. What these numbers are NOT

The primary is a **6-cell prediction test against a registered residual bar, not a powered
effect-size comparison**: the S7-6 rungs are published aggregates, no per-target rows survive, so
no SE, MDE or fold CI is computable for them and none is quoted -- which is exactly why the
falsifier was registered as a residual bar. Comparisons made by this lane: 6 cells x 3 predictors,
4 union arms, 6 slopes, and the H-X2 battery; only H-X1-primary and H-X2 were pre-registered and
only those two are read as results. No per-target maximum is taken anywhere here, so no
`best_of_k_within` arm is required. Every ceiling row is ORACLE and says so.


## S30-L11 -- THE SPARSE WEIGHTED READOUT IS **INFORMATION-DOMINATED BY THE ARGMIN IT WAS MEANT TO BEAT**: AT EVERY BUDGET FROM 3 TO 9 BITS THE PLAIN `argmin` OVER THE TOP-2^B IS BETTER THAN THE BEST FULLY-PRICED SPARSE ARM, BY **+0.162 TO +0.727 Å**; ON THE BUILT CHAIN, 2 OF THE TOP-75 WITH **FREE CONTINUOUS WEIGHTS** (11.4 BITS + UNBOUNDED) IS 2.1683 Å AGAINST THE TOP-128 ARGMIN'S **2.1435 Å AT 7.0 BITS**. S29's LAST UNCLOSED LADDER CLASS IS CLOSED -- NOT BY CEILING, BY **PRICE**. PLUS THE SPRINT-RELEVANT INCIDENTAL: ON FAIL18 THE ORACLE BEST CANDIDATE IS **HIDING AT RANK 128-500**, AND WIDENING THE REGISTER 128->512 IS WORTH **-1.90 Å ON THE TAIL AGAINST -0.19 Å ON THE OTHER 108** (2.77x MDE) (2026-09-20 12:57, Q)

Pre-registration `s30/PREREG_S30_Q_sparse.md`, committed at 05455a88 before the first number
existed. Code `s30/s30_Q_sparse.py`, `s30/s30_Q_analyse.py`. Rows
`s30/results/s30_Q_sparse_rows.jsonl` (126), summary `s30/results/s30_Q_sparse.json`.
**Every arm is ORACLE except the `D_*` and `H_*` rows**, whose construction reads no native.
Basis is the POINT CLOUD except the section explicitly headed BUILT CHAIN, which re-uses S29 lane
O's already-chained arms and computes nothing new.

## WHAT WAS ASKED

S29 left a sparse weighted readout as "the only ladder class not closed by ceiling": two members
with ORACLE weights emit 1.4315 Å on the chain where 75 members with ORACLE membership emit
2.3055 and production emits 3.2105. The coordinator's framing, which this lane adopted: **a low
parameter count is not a low information requirement.** Choosing 2 of 500 is log2 C(500,2) =
**16.93 bits**, more than the 7 bits the deployed top-128 argmin spends, not fewer -- and the
1.4315 Å arm additionally fits **free continuous weights per target against the native**, a second
and entirely unpriced channel.

So the question is not "how low can this class go" but **"what does it cost, in the same currency
the deployed readout already spends".**

## THE BIT LEDGER (both channels counted the same way)

- support: log2 C(|S|, s); **0** for a support fixed by a native-free rule.
- weights: log2 C(L+s-1, s-1) for weights on the 1/L simplex lattice; **0** for uniform;
  **unbounded** for continuous.
- reference: the argmin over the top-2^B of the score order costs **exactly B bits**.

## THE FOUR PRE-REGISTERED FALSIFIERS

**F1 CONFIRMED -- the support does carry most of the gain, but only while s is small.**
At s = 2 over the pool, ORACLE support with UNIFORM weights recovers **76.7%** of the
production-to-ORACLE gain (A 1.3631, B 1.7565, gain 1.6852). The share **falls monotonically with
s**: 78.7 / 75.9 / 60.4 / **33.9%** at s = 3 / 5 / 10 / 20, and on the top-75 set 88.6% down to
18.7%. So "which members" dominates at s = 2 and "with what weights" takes over by s = 20.

**F2 CONFIRMED -- the weights are cheap in bits.** At s = 2 the 1/L lattice is within 0.05 Å of
the continuous optimum at **L = 4, i.e. 2.32 bits** (1.406 vs 1.363); L = 8 costs 3.17 bits and is
within 0.011. The weight channel is NOT what makes this class expensive.

**F3 CONFIRMED AND WORTHLESS, WHICH IS THE INFORMATIVE OUTCOME.** A native-free support rule does
beat a random one -- best of {score-prefix, farthest-point diversity} beats the mean of 8 random
supports by **+0.364 Å** at s = 2 (3.198 vs 3.562), and by +0.699 / +0.586 / +0.379 / +0.219 at
s = 3 / 5 / 10 / 20. **But every native-free support arm is WORSE than production** (3.198 against
3.048), even with ORACLE continuous weights handed to it for free. The skill exists, is real, and
is nowhere near enough. Reported as a confirmed falsifier whose confirmation changes nothing.

**F4 REFUTED -- and this is the lane's result.** At equal bits the argmin wins everywhere in the
registered window:

```
 budget  best fully-priced sparse arm          argmin over top-2**B     sparse - argmin
 B = 3   H_top75_s20   (0.0 b)  3.0662          2.9040                    +0.1622
 B = 4   H_top75_s20   (0.0 b)  3.0662          2.7159                    +0.3503
 B = 5   T8_s2_unif    (4.8 b)  2.8805          2.5115                    +0.3690
 B = 6   T8_s2_unif    (4.8 b)  2.8805          2.3432                    +0.5374
 B = 7   T16_s2_unif   (6.9 b)  2.6517          2.1458                    +0.5059
 B = 8   T16_s2_unif   (6.9 b)  2.6517          1.9383                    +0.7134
 B = 9   T32_s2_unif   (9.0 b)  2.4375          1.7108                    +0.7267
```

The sparse family needs **~15 bits** (T256_s2_unif, 1.8317) to reach what the argmin delivers at
**9** (1.7108), and **~19-20 bits** to reach 1.41. The three "sparse wins" at B = 0, 1, 2 are
outside the registered window and are not a route: they say only that a 20-member uniform average
(3.0662) beats the best of 1, 2 or 4 candidates -- which is production's existing design, not a
new one.

**On the BUILT CHAIN, using only arms S29 already chained, so nothing here is a projection:**

```
 best1_top128     argmin over top-128           7.0 bits, no weights      2.1435 A
 sparse_top75_s2  2 of top-75 + CONTINUOUS w    11.4 bits + UNBOUNDED     2.1683 A
 sparse_pool_s2   2 of 500 + CONTINUOUS w       16.9 bits + UNBOUNDED     1.4315 A
 best1_pool       argmin over the 500           8.97 bits, no weights     1.7078 A
```

The argmin beats the 2-of-75 sparse arm **by 0.0248 Å while spending 4.4 fewer bits and no weight
channel at all**. The famous 1.4315 Å is real and is bought for at least 16.9 bits plus an
unbounded continuous channel; the same pool's argmin gets to 1.7078 for 8.97 bits and nothing else.
**The 1.4315 Å figure was never a route. It is an expensive way to buy what 9 bits already buys.**

## CONTROLS

- **Greedy-gap.** The supports are chosen by forward greedy; exhaustive s = 2 over the top-64
  (C(64,2) = 2016 pairs enumerated per target) gives 2.1180 against greedy's 2.1427, a **+0.0247 Å**
  gap. The greedy arms are faithful, not lower bounds, at this resolution.
- **Order statistic (contract rule 8).** The random-support arm is the **mean of 8 draws**, never
  the per-target min. Priced: observed min-of-8 gain -1.1056 Å is **108% accounted for by the
  across-target null** and the split-half transfer is **+0.0007 Å (0%)**. The random-support
  minimum is pure best-of-k and is used nowhere as an arm.
- Every comparison is made on one basis at a time and each figure names it.

## THE INCIDENTAL FINDING, WHICH MAY MATTER MORE THAN THE REGISTERED ONE

The argmin reference curve was computed per stratum because S30's headline says the mean is a tail
statistic. **The deployed top-128 register is what fails the tail.**

```
 B     top-N     all 126    FAIL18    other 108
 7     128        2.1458    4.1846      1.8060
 8     256        1.9383    3.4526      1.6859
 9     500        1.7108    2.2842      1.6153

 widening 128 -> 512:      -0.4350   -1.9004     -0.1907
```

The difference of stratum means is **-1.7097 Å, SE 0.2207, 2.77x MDE**. On the 18 hardest targets
the ORACLE best candidate in the pool sits **at rank 128-500**, outside the window the quantum
stage sees (`core/pipeline.py:758`). On the other 108 widening is worth almost nothing.

**Stated with its caveat, which is severe:** this is ORACLE -- it says the candidate *is there*,
not that anything can find it, and S29 closed recognition three ways. It relocates the tail's
failure from "the pool does not contain the answer" to "the register does not contain the
candidate", which are different problems with different costs. Two more qubits widen 128 -> 512.

**It also composes with lane F's width measurement in the opposite direction, and the pair is the
real statement.** Lane F (`s30/results/s30_F_width_cloud.json`, phase 1) has the ORACLE **set
mean** getting *worse* with width on the other 108 (3.116 at 75 -> 4.201 at 500) and *better* on
FAIL18 (6.159 -> 5.755 at 400). Mine has the **argmin** getting dramatically better with width on
FAIL18. Together: **widening the register hurts an averaging readout and helps a selecting one, and
on the tail the selecting readout gains 1.90 Å.** Neither lane's number says that alone. Flagged to
lane F rather than claimed here.

## WHAT THIS DAMAGED IN MY OWN EXPECTATIONS

1. I expected F1 to be the decisive one and F4 to be a formality. It was the reverse: the support
   *does* carry the gain (F1 confirmed at 76.7%), and that turned out to be exactly why the class
   loses -- the support is the expensive channel, and the argmin buys the same thing cheaper.
2. I expected the weights to be the hidden cost. They are the cheapest part of the whole family
   (2.32 bits at s = 2).
3. I expected to be arguing about whether a native-free support rule could be built. F3 says one
   already has measurable skill over random -- and it does not matter, because the whole
   native-free branch sits above production.

## COMPARISONS MADE (contract rule 26)

Four pre-registered falsifiers; 10 F1 cells, 5 F2 rows x 7 lattice resolutions, 5 F3 cells,
21 F4 budgets, 1 greedy control, 2 order-statistic prices, 1 stratum split of the reference curve.
Only F1-F4 are read as results. The stratum split of the argmin curve was **not** pre-registered
and is reported as an incidental finding with its SE and MDE attached, not as a lane result.

> **NUMBERING COLLISION RESOLVED 2026-09-20 13:05 by the coordinator.** Posted as `S30-L3`, which an earlier entry already held -- lanes picked the next free number from a ledger that grew between their reads. Renumbered to **S30-L16**; the earlier entry keeps `S30-L3`. Recorded rather than silently fixed, per contract rule 15. Anything citing `S30-L3` for *this* result means `S30-L16`.

## S30-L16 -- FILTER WIDTH IS **NOT** AVERAGING WIDTH: IT MOVES THE EMITTED STRUCTURE 4-15× AS MUCH, SO S29'S FLAT m SWEEP DOES NOT TRANSFER. BUT THERE IS **NO FREE LUNCH IN WIDTH** -- F2a REFUTED, THE BENEFIT/HARM TRADE-OFF IS CLEAN AND CONTINUOUS AND CROSSES AT ~55%/55%. F2b CLOSED **BY CEILING WITHOUT SPENDING THE CHAIN**: THE ORACLE GLOBAL BEST FILTER WIDTH **IS THE SHIPPED VALUE**. AND THE "WIDENING HELPS HARD TARGETS" READING IS **FALSE** -- IT DOES NOT REPLICATE ON EITHER FILTER-INDEPENDENT TAIL. ONE MECHANISM (ρ = +0.81) EXPLAINS ALL FOUR STRATA (2026-09-20 12:57, F)

Pre-registration `s30/PREREG_S30_F2.md`, committed **e626a88e at 12:47:48**, before this file's
code existed. Code `s30/s30_F_width.py`; rows `s30/results/s30_F_width_rows.jsonl` (126);
analysis `s30/results/s30_F_width_cloud.json`.

**Basis: every RMSD below is a POINT CLOUD figure.** No built chain was spent, and §4 says why
that is a closure rather than an omission (contract rule 1).

**Reproduction gate.** The score order was recomputed from the posterior rather than read from a
stored order, and the recomputed top-75 equals the production record's `sub` **set-wise on
126/126 targets**, aborting otherwise. Declared caveat: the posterior itself comes from the s12
cache, and S29-L19's operational note (`s29/LEDGER.md:751`) warns that an operator reading the
score's order *below* the top-75 cut carries a 2-in-126 chance of differing from a fresh
recomputation. The gate certifies the top-75 **set**, not the deep order; every k > 75 row
inherits that caveat.

### 1. The confound this removed

In production the top-75 is **both** the filter and the averaging set, so one number moves two
mechanisms. S29 swept m as the **averaging** width and found it flat (−0.00047 Å per unit m).
Separating them — filter to top-k by the shipped Bayes-risk score, then coordinate-average a
**random** m-subset of those k, 8 seeds, m = 75 held fixed — gives:

```
axis                                            range of the emitted point-cloud mean
filter width      k = 75 -> 500, m = 75 fixed        0.3790   (3.0483 -> 3.4273)
averaging width   m = 25 -> 75,  k = 75 fixed        0.0256
averaging width   m = 25 -> 150, k = 500 fixed       0.0927
```

**The filter width is worth 4× to 15× the averaging width, on the same instrument.** The
coordinator's prior — that S29's flatness is a property of the averaging width and says nothing
about the filter width — is confirmed by measurement rather than assumed. **F2c is refuted**
(max |×MDE| over k = 1.59 at k = 500, 5/5 folds), and that is the durable methodological output
of this entry: *the project has one published "m is flat" result and it does not license any
statement about the filter.*

### 2. F2a REFUTED -- width is a trade-off, not a free lunch

The registered clause needed some k > 75 retaining **≥ 80%** of the 108's ORACLE set-mean benefit
while shedding **≥ 50%** of FAIL18's ORACLE best-member harm. Nothing does (both columns ORACLE):

```
k     benefit retained on the 108     harm shed on FAIL18     fires
100            97.9%                        8.1%               no
128            95.0%                       20.6%               no
150            92.1%                       26.4%               no
200            84.2%                       41.6%               no
250            73.2%                       50.5%               no
300            60.3%                       59.6%               no
400            33.6%                       78.0%               no
500             0.0%                      100.0%               no
```

The two curves are smooth and cross near k ≈ 275 at about 55%/55%. There is no knee, no
saturation of the benefit before the harm starts shedding, and therefore **no width at which the
filter is nearly as good on the body and much less bad on the tail.** The hypothesis was worth
testing and it is wrong.

### 3. F2b CLOSED BY CEILING -- and the chain was not spent

Emitted point-cloud mean against filter width, m = 75 (paired vs production at k = 75):

```
k       mean     vs prod    ×MDE   fold CI             folds   W/L        verdict
75    3.0483    +0.0000       —          —             5/5        —       production
100   3.0534    +0.0051    +0.13  [-0.0229,+0.0313]    3/5    58W/68L     NOT MEASURED
128   3.0598    +0.0114    +0.18  [-0.0105,+0.0334]    3/5    62W/64L     NOT MEASURED
150   3.0729    +0.0245    +0.30  [-0.0102,+0.0637]    3/5    63W/63L     NOT MEASURED
200   3.0972    +0.0488    +0.43  [-0.0046,+0.1117]    4/5    67W/59L     NOT MEASURED
250   3.1412    +0.0929    +0.67  [+0.0172,+0.1666]    4/5    58W/68L     NOT MEASURED
300   3.1890    +0.1407    +0.81  [+0.0352,+0.2443]    4/5    53W/73L     WORSE
400   3.2962    +0.2478    +1.18  [+0.1204,+0.3793]    4/5    46W/80L     WORSE (type-M)
500   3.4273    +0.3790    +1.59  [+0.2243,+0.5217]    5/5    45W/81L     WORSE
```

**The ORACLE global argmin over k is k = 75 — the shipped value.** A deployable leave-fold-out k
can at best tie production and cannot beat it, because the quantity it would be estimating has its
optimum at the incumbent. F2b is therefore closed **by ceiling**, in the same way lane O closed its
rungs, and no built-chain run was spent on it. The ordering is safe to carry to the chain without
running it: lane O's D5 established that the projection price is a monotone increasing function of
the cloud's own accuracy (corr +0.866), so every arm here — all of which have *worse* clouds than
production — would pay *more* on the chain, not less. The chain cannot reverse this sign.

This is the **fifth** instance of S29's single-global-scalar pattern (`s29_O_FINDINGS.md` U1 lists
four). I pre-registered that I expected F2b to fail, for that reason, before running it.

### 4. The result I did not expect, and the honesty check that changes how to read it

Per stratum, effect of widening versus production (point cloud, m = 75):

```
k        FAIL18    other 108   worst18_poolmean   worst18_bestpool
100     -0.0176     +0.0089        +0.0862            +0.0127
200     -0.2349     +0.0961        +0.3088            -0.0209
300     -0.4765     +0.2436        +0.5550            +0.0041
400     -0.6193     +0.3924        +0.7313            -0.0081
500     -0.5905     +0.5406        +0.9566            +0.1121
```

On FAIL18 widening helps monotonically to an interior optimum at k = 400 (**−0.6193**, random-18
null p = 0; at k = 200 already p = 0.0026), and on the 108 it hurts monotonically. **Taken alone
that row reads as "a wider filter rescues hard targets", and that reading is false.** On the two
tails that are *not* defined by production — the same two S30-L2 used — the effect is absent
(worst-18 by ORACLE best-in-pool: −0.02 to +0.11, no trend) or **reversed** (worst-18 by pool
mean: +0.7313 at k = 400, the same sign as the body). Hardness as such does not predict that
widening helps. Only the production-defined tail shows it.

**A linear regression-to-the-mean control does not explain it either**, so it is not simply that
FAIL18 is where production was unluckiest: fitting the per-target effect on production RMSD over
the 108 and extrapolating gives a predicted **+0.069** at k = 400 against an observed **−0.619**,
residual −0.688.

**The mechanism that does explain all four strata, at ρ = +0.81.** Let x be the filter's ORACLE
set-mean benefit (pool mean member RMSD minus top-75 mean member RMSD; positive = the filter
improved the set mean). Then the per-target effect of widening is a function of x:

```
Pearson(x, effect of widening)      k=200 +0.574   k=300 +0.691   k=400 +0.774   k=500 +0.812
Spearman                                  +0.515         +0.605         +0.703         +0.772

stratum              x (ORACLE)     effect at k=400
all                   +0.9026          +0.2478
other 108             +1.0855          +0.3924
FAIL18                -0.1947          -0.6193
worst18_poolmean      +0.9355          +0.7313
worst18_bestpool      +0.3515          -0.0081
```

**Widening helps exactly where the filter hurt the set mean, and hurts where it helped** — which
is the `operator-consumes-set-mean` law (d_out = 1.16·d_set_mean + 0.04·d_set_best) read backwards.
FAIL18 is largely that condition *by construction*, which is why it alone shows the gain; the
worst-18-by-pool-mean have bad pools **but a filter that still works** (x = +0.94), so widening
costs them. This dissolves the stratum table into one quantity and supersedes "widening helps hard
targets" as the way to state it.

### 5. The ORACLE-gated arm, and the one number that changes an old economic verdict

Gating on the mechanism rather than on the label — widen to k = 400 only where the filter hurt the
set mean (**16 of 126** targets; 8 of them FAIL18):

```
                                          effect    ×MDE   fold CI            folds   W/L
on the 16 gated targets                  -0.9397   -2.38  [-1.1235,-0.7033]   4/4   16W/0L
over all 126 (the endpoint value)        -0.1193   -1.29  [-0.1859,-0.0446]   4/5   16W/0L
```

**Both rows are ORACLE in every sense** — the gate reads the native to know x, and no native-free
rule here supplies it. It is a router, and the router family is closed.

**But the size is the point, and it revises an economic argument the project has been running on.**
Lane C's FAIL18 detector (S28-L6) closed on *two* grounds, and the second was that the ORACLE prize
was under the instrument: **−0.0474 Å built chain, 0.62× MDE**, with the note that a perfect
detector would need ~330 targets to be measurable. **This gate's ORACLE prize is −0.1193 Å point
cloud at 1.29× MDE — roughly 2.5× larger, and above its own MDE rather than below it.** That does
not reopen the router family, whose features are closed on the *first* ground (no block clears its
permutation null). It does mean the "even a perfect detector is not worth measuring" half of that
closure is **specific to the operator lane C was switching** and must not be quoted as a general
property of FAIL18 detection. A detector built on a feature family demonstrably outside the closed
seven would now be aiming at a prize 2.5× the one that was dismissed as too small.

### 6. What this lane does next, and what it will not do

Not a ninth router: the eight closed constructions span length, the objective's own score
distribution, posterior entropy, candidate-set geometry, compactness disagreement (`rg_z` — already
closed at 2% of its own MDE, S23 L7), ESM contact statistics, and pool statistical-potential
agreement. The next measurement is §3's question — *what is the distogram confidently wrong about
on the pools where its own filter inverts* — and the cheapest decisive form of it is a per-target
Spearman between the shipped score and ORACLE in-pool RMSD over all 500 members, split by stratum
with the random-18 null. **That measurement does not exist in the record**: the pool-ranking
measurements (S29-L33, S29-L50, S14) are unstratified, and the FAIL18-stratified measurements
(S28-L23b's gradient cosine, S29-L2's ladder rho, S29-L20's axis cosine) are structure-level
cosines over hand-built ladders, not within-pool ranking skill.

Artefacts read for this entry: `s29/s29_O_FINDINGS.md` U1 and D5; `s29/LEDGER.md:751` (the
recompute-the-posterior warning, adopted); `s27/LEDGER.md:229-321` and `s27/s28_C_FINDINGS.md:62-73`
(lane C's ORACLE prize, the comparison in §5); project memory `operator-consumes-set-mean`,
`grid-oracles-are-order-statistics`, `mde-is-per-comparison-not-per-instrument`.


## S30-L12 -- THE SECOND-MOMENT (QUADRIC) TAIL CLASS **DOES NOT BEAT THE HALFSPACE CLASS IT WAS MEANT TO ESCAPE INTO**: A 17x VC EXPANSION (34 -> 595) BUYS **-0.0105 Å AT 0.16x MDE (NOT MEASURED)** AT ORACLE m AND IS **+0.2059 Å WORSE AT 1.81x MDE (W/L 37/89)** AT THE SHIPPED m=75. AND THE WHOLE APPARENT CEILING OF BOTH CLASSES IS AN ORDER STATISTIC: **225% OF THE LINEAR SEARCH'S GAIN IS ACCOUNTED FOR BY THE ACROSS-TARGET NULL AND THE SPLIT-HALF TRANSFER IS +0.032 Å (-3%)**. THE DISPERSION RULE THE CONSTRUCTION POINTS AT, MEASURED DIRECTLY, IS **3.3585 Å AT m=75 AGAINST PRODUCTION'S 3.0483** -- WORSE THAN WHAT SHIPS. MEASURED BEFORE A HAMILTONIAN WAS BUILT, AS INSTRUCTED (2026-09-20 13:01, Q)

Code `s30/s30_Q_quadric.py`, `s30/s30_Q_quadric_analyse.py`. Rows
`s30/results/s30_Q_quadric_rows.jsonl` (126), summary `s30/results/s30_Q_quadric.json`.
POINT CLOUD; **every ceiling row is ORACLE**; the `rules` rows are native-free constructions with
an ORACLE readout. No deployable parameter is chosen.

## WHAT WAS ASKED AND HOW IT WAS MADE FALSIFIABLE

The coordinator's construction (from lane T's theorem T1, S30-L9): with an **endogenous** order the
cost is `c_x = <grad f(R_lambda), W_x>`, so candidates are sorted by projection onto one
self-consistent direction and the reachable tails are **halfspace cuts**, VC dimension d+1 = 34.
Adding a second-moment term `V = f(sum l W, sum l W W^T)` makes `grad V_x` **quadratic** in `W_x`,
the cuts become **quadrics**, and the VC dimension jumps to (d+2)(d+1)/2 ~ 595 -- "the only
construction that survives lane T's cap."

**The measurement that decides it before any circuit exists:** enumerate what each class can
actually reach on the real pools, with **the readout held fixed** (the uniform coordinate average
of the selected prefix) so that the only thing varying is *which set*. 128 sampled directions per
class per target, all prefix sizes m = 1..500, 126 targets.

## THE RESULT: THE EXPANSION BUYS NOTHING

```
                                       ORACLE m        at the shipped m = 75
  deployed exogenous prefix             2.6355                3.0483
  LINEAR  (halfspace), best of 128      1.7356                2.2204
  QUADRIC (2nd moment), best of 128     1.7251                2.4263
  ORACLE-sorted prefix (upper ref)      1.3958                    --

  QUADRIC minus LINEAR, ORACLE m      -0.0105  SE 0.0230  0.16x MDE  W/L 65/57   NOT MEASURED
  QUADRIC minus LINEAR, m = 75        +0.2059  SE 0.0405  1.81x MDE  W/L 37/89   WORSE
```

A seventeen-fold expansion of the reachable-tail class is worth **0.16x MDE** where it is free to
pick its own set size, and is **significantly harmful** at the set size the pipeline actually
emits. The richer class is not merely not-better; at m = 75 it is reliably worse, losing on 89 of
126 targets.

**Independently corroborated by lane T's own quadric rows** (`s30/results/s30_T_quadric.json`,
K = 5000 directions inside a 6-dimensional PC subspace at fixed M = 75 -- a different sampler, a
different subspace, a different K): on both of its completed targets `quad_best` is worse than
`half_best` (1A13 2.207 vs 2.117; 1A1P 2.340 vs 2.316). Two samplers, two lanes, same sign.

## AND THE CEILING THAT IS THERE DOES NOT TRANSFER

The headline "halfspace cuts reach 1.74" is a **best-of-128**, and pricing it kills it:

```
  best_of_k_within, LINEAR   observed -1.2003  across-target null -2.7006 (225% accounted)
                             split-half +0.0323 (-3% transfer)   k_eff 78.8
  best_of_k_within, QUADRIC  observed -1.1906  across-target null -2.4193 (203% accounted)
                             split-half -0.0204 (2% transfer)
```

The across-target null **exceeds** the observed gain in both classes, and split-half transfer is at
or on the wrong side of zero. **There is no transferable direction.** This is S29's U1
incidental-parameter pattern again -- previously measured for a scalar step (rungs 6 and 8), now
measured for a whole direction: per-target it looks like a large gain, across targets it is worth
nothing. The cut class's expressiveness is real and unusable for the same reason every free scalar
in this project has been.

## THE STRUCTURED VERSION OF THE PROPOSAL, TESTED DIRECTLY

The random-quadric sampler could be accused of missing the *right* quadric. So the specific
quadratic form the second-moment functional's gradient produces was evaluated as a named
native-free rule: `disp2(x) = ||W_x - mean(W)||^2`, the candidate's squared deviation from the
pool mean -- "operators that read the pool's own dispersion."

```
  native-free rule        ORACLE m     at m = 75
  dis_score (deployed)     2.6355        3.0483
  disp2                    3.0606        3.3585      <- the dispersion rule
  rg                       2.7858        3.5382
  dist_to_medoid           3.1030        3.7246
  PC1 / PC2 / PC3          ~3.0          4.0-4.2
```

**Every native-free rule direction tested is worse than the deployed score at both set sizes, and
the dispersion rule is 0.31 Å worse than production at m = 75.** The two-routes-converging argument
(S29's post-mortem and T's derivation both pointing at pool dispersion) is real as an argument and
does not survive contact with the pools.

## ONE THING THAT DID SURPRISE ME, AND A WORRY THAT WAS UNFOUNDED

I expected the cut classes' ORACLE ceiling to be attained by degenerating to m = 1 -- i.e. for the
"richer selection class" to collapse into the argmin readout, which would have made it a
relabelling of a question already answered. **It does not.** The best prefix for a linear direction
sits at **median m = 402**, with 81% of directions peaking at m >= 50 and only 4% at m <= 2. The
class genuinely selects large sets. It simply cannot select them from outside the target.

## SCOPE, STATED SO THE NEGATIVE IS NOT OVER-READ

- 128 sampled directions in the full 3n space is not exhaustive, and neither is T's 5000 in 6
  dimensions. But the order-statistic pricing is what closes this, not the ceiling: a larger search
  inflates an untransferable number. To reopen it, someone must exhibit a **rule** that produces a
  direction, not a larger search over directions.
- This measures the SELECTION side with the readout held fixed. It says nothing about readouts --
  which is S30-L11's subject, and there the argmin dominates.
- The quadric family sampled is `A` of rank 3 with random signs plus a linear part. A structured
  `A` derived from a trained objective is not tested here; the one structured form the construction
  itself specifies (`disp2`) is, and it is negative.

## COMPARISONS MADE (contract rule 26)

Two classes x 8 best-of-K points x 2 set-size conventions, 2 order-statistic prices, 14 native-free
rule directions (7 rules x 2 sign conventions), 3 stratum splits. The two contrasts read as results
are QUADRIC-minus-LINEAR at ORACLE m and at m = 75; both were fixed as the decision before the run
and both are reported with SE, MDE and W/L. No per-target maximum is read as a mean anywhere.

## S30-L13 -- I RE-CAP MY OWN ESCAPE E2: THE FIELD HAS RUN IT AT SCALE AND IT COSTS ACCURACY. MCORE TAKES AVERAGING CLASHES FROM **63.0% OF ATOMS TO 1.09%** AND RMSD FROM **3.28 A TO 3.36 A -- 0.08 A WORSE**, ON 2090 PROTEINS. OUR -0.022 A IS THE **OPPOSITE SIGN** FROM THE PUBLISHED RESULT, AND THE MECHANISM PREDICTS WHY: THE BONDED FRACTION OF ALL PAIRS GOES AS ~2/n, SO A CONSTRAINT REPAIR IS **~15% OF THE GEOMETRY AT n=13 AND ~1% AT n=200**. PERCEPTION-DISTORTION RE-CAPS E2 IMMEDIATELY AFTER IT ESCAPES THE IDENTIFICATION INVARIANCE (2026-09-20 13:01, L)
Addendum to S30-L7, correcting my own framing there. Note: `s30/lit/L30_2_common_mode.md`
section 4. **Literature reading plus arithmetic; no new measurement.**

**WHAT S30-L7 SAID, AND WHERE IT WAS INCOMPLETE.** S30-L7 proved the pool's shared bias `mu` is
non-identifiable and listed exactly three escapes. I called **E2** -- a hard constraint that the
native satisfies and the contracted average does not -- "the only escape this project has ever
obtained a signed fold-consistent result from", at -0.022 A [-0.036, -0.009], 5/5 folds
(`averaging-space-beats-the-objective`). That is true. **What I did not say is that escaping the
identification invariance does not exempt you from the OTHER theorem, and E2 walks straight into
it.**

**THE PUBLISHED NUMBERS** (*Improving consensus structure by eliminating averaging artifacts*,
PMC2662860; 2090 non-homologous single-domain proteins under 200 residues):

    baseline averaging (COMBO)      63.0 %  of atoms in clashes < 3.6 A
    MCORE (their repair)             1.09 %                                <- a 58x reduction
    PULCHRA                          3.64 %
    RMSD, MCORE refined              3.36 A   against   3.28 A original    <- +0.08 A, WORSE

**A repair that removes 98% of the clashes costs 0.08 A of accuracy.** That is a realism operator
improving realism and losing distortion, which is Blau & Michaeli (CVPR 2018) Theorem 3 -- already
imported in S29-L12 -- arriving in a place S29 did not apply it. **So the correct statement is:
E2 breaks the identification invariance of S30-L7 and is then immediately re-capped by
perception-distortion.** The two theorems compose, and the composition is what sets the size of the
prize.

**WHY OUR RESULT HAS THE OPPOSITE SIGN, AND IT IS A FALSIFIABLE MECHANISM.** Ours is -0.022 A
(better) on 126 peptides of 9-16 aa; the published one is +0.08 A (worse) on proteins under 200. A
constraint repair acts on the **bonded** part of the geometry, and the bonded fraction of all
residue pairs is `(n-1) / C(n,2) = 2/n`:

    n = 13   ->  12 Ca-Ca bonds among 78 pairs   =  15.4 %
    n = 200  ->  199 among 19,900                =   1.0 %

**A constraint repair touches ~15x more of the geometry at peptide length than at protein length.**
That is a mechanism for the sign flip, it is arithmetic rather than a story, and it is falsifiable
three ways: (i) the gain should fall with `n` across our own 9-16 range, (ii) it should fall as the
repair is restricted to non-bonded constraints, and (iii) it predicts our gain is essentially a
bond-geometry gain and not a clash gain.

**CONSEQUENCES, IN PRIORITY ORDER.**

1. **Do not scale the AMBER-relax result.** The field ran the general version at 15x our sample and
   got the wrong sign. Our -0.022 A is an outlier in the literature's direction, is 0.7% of
   baseline, and section "why the sign flips" says the mechanism does **not** grow with effort --
   it is capped by the bonded fraction, which is fixed by chain length.
2. **The cheap re-read that is worth doing, and it uses a measurement that already exists.** The
   published paper states *"averaging artifacts become more pronounced when members of the ensemble
   are more divergent"* but provides **no quantitative version** (I checked; the correlation with
   ensemble spread is asserted, never measured -- a gap in the literature, not a result I can
   import). Our pool's divergence is observable per target (`mean_k |d_k|^2`,
   `s23/results/errdecomp.json`). So: **split the EXISTING k=30 AMBER-relax measurement by pool
   dispersion and by FAIL18/108.** It costs no new compute, it tests a literature claim nobody has
   tested, and it is aimed exactly where S30-L0's arithmetic says the sprint should aim. If the
   -0.022 A concentrates on the divergent tail it is a tail intervention rather than a 0.7%
   aggregate; if it is uniform, E2 closes.
3. **A prediction I am registering before anyone looks**, so it is falsifiable: the repair's gain
   concentrates on high-dispersion targets. Stated at 2 to 1, not more -- the field asserts the
   mechanism but never measured it, and S29's own experience is that asserted mechanisms of this
   kind fail about as often as they hold (S29-L50 refuted a 3-to-1 prior of lane T's and my own
   expectation in the same direction).

**A CLOSURE THIS COMPLETES: E1 IS DEAD THREE INDEPENDENT WAYS.** With E2 re-capped, the accounting
in S30-L7 tightens. Escape E1 (a prior on the bias's FORM) is closed on this instrument by three
separate arguments, and I want them in one place because each alone looks escapable:
  (a) the **scale** form -- the optimal rescale `s* = <c,t>/|c|^2` is a function of the invisible
      component and of nothing else (memory `pool-error-is-68-percent-common-mode`), which is why
      the mismatched-native placebo works;
  (b) the **flexible** form -- a corrector trained on the predictor's own features inherits its
      error structure (memory `error-coherence-decides-correctors`: at identical 0.688 sign
      accuracy, coherent mistakes emit +0.31 A and i.i.d. mistakes -0.14 A), and S24 closed
      boosting/residual fitting on this instrument;
  (c) the **rigid global** form -- S29-L47 measured four global scalars whose ORACLE global optima
      are 0.0-0.6% of their per-target gains, **two of them exactly zero**.
Brynjarsdottir & O'Hagan (2014) require an informative prior on the discrepancy's *shape*; (a),
(b) and (c) are the three shapes available here and all three are measured closed.

**WHERE THAT LEAVES THE SPRINT'S HOPE.** Of S30-L7's three escapes: **E1 closed** (three ways),
**E2 open but doubly capped** (identification escaped, perception-distortion re-applied, bonded
fraction fixed by `n`, published sign wrong, our size 0.7% of baseline), **E3 the only one whose
cap is not yet a theorem** -- and its obvious instance is refuted (provenance cosine 0.9432 above a
0.9330 within-source control), leaving "score against a different prior" rather than "draw from a
different source", which is `prior-derivative-is-the-only-steep-lever` at -2.15 A/unit.

**WHERE I COULD BE WRONG.** The 63.0% / 1.09% / 3.28 -> 3.36 numbers are from one paper on one
benchmark, read from the article rather than recomputed, and their averaging protocol is not ours
(they average full-atom models across clusters; we average Ca windows in a medoid frame after
Kabsch). The `2/n` bonded-fraction argument treats all pairs as equally weighted by the repair,
which no real force field does -- it is an order-of-magnitude argument for the sign flip, not a
prediction of its size. And prediction 3 is mine, at 2 to 1, and the last time this lane and lane T
agreed on a prior of this kind at 3 to 1 we were both wrong at n = 126.

> **NUMBERING COLLISION RESOLVED 2026-09-20 13:05 by the coordinator.** Posted as `S30-L4`, which an earlier entry already held -- lanes picked the next free number from a ledger that grew between their reads. Renumbered to **S30-L17**; the earlier entry keeps `S30-L4`. Recorded rather than silently fixed, per contract rule 15. Anything citing `S30-L4` for *this* result means `S30-L17`.

## S30-L17 -- THE SCORE HAS **NO IN-POOL RANKING SKILL ON THE TAIL**: ρ = **+0.1066** (FOLD CI INCLUDING ZERO) AGAINST **+0.6446** ON THE 108, RANDOM-18 NULL **p = 0**. IT IS ALMOST ENTIRELY THE DISTOGRAM'S OWN ERROR (ρ = **−0.799**), THE CHAIN IS FULLY MEDIATED, AND THE ERROR THAT MATTERS IS **SHAPE, NOT SCALE** -- LANE L'S SCALE-SEPARABLE TERM IS **REFUTED AS THE MECHANISM** (PARTIAL −0.067, p = 0.46, AGAINST SHAPE'S −0.641). "CONFIDENTLY WRONG" ALSO REFUTED: ERROR×CONFIDENCE IS **WORSE** THAN ERROR ALONE (2026-09-20 13:03, F)

> **Cross-reference repair, 2026-09-20 13:07, F.** This entry was written as S30-L4 and cited its sibling as "S30-L3"; the coordinator's collision fix renumbered that sibling to **S30-L16** (and this entry to L17) without rewriting the bodies, so four in-text references here pointed at lane D's S30-L3 instead. They now read S30-L16. Only the number changed; no claim, figure or wording was altered, and no other lane's entry was touched.


Pre-registration `s30/PREREG_S30_F3.md`, committed **eb719398 at 12:58:34**, before this file's
code existed. Code `s30/s30_F_score.py`; output `s30/results/s30_F_score.json`.

**Every predictor in this entry is ORACLE** (it reads the native's distances). This experiment
builds **no router and no detector**, by design — the question is *what the score is wrong about*,
not whether we can spot it. The seven closed router feature families are enumerated in the prereg
and every one of them is measured here and reported whether flattering or not.

Reproduction gate as S30-L16: the score order is recomputed from the posterior and its top-75
equals production's `sub` set-wise on 126/126. Same declared caveat about deep order below the cut
(`s29/LEDGER.md:751`).

### 1. F3a FIRES -- and this quantity did not exist in the record

Per-target Spearman between the shipped Bayes-risk score and ORACLE CA-RMSD, over each target's
full 500-member pool. Positive = the score orders its own pool correctly.

```
stratum               n    rho_pool    rho_top75(in-band)
all                 126     +0.5678        +0.0652
other 108           108     +0.6446        +0.0743     fold CI [+0.5905,+0.7063]  5/5 folds
FAIL18               18     +0.1066        +0.0103     fold CI [-0.0241,+0.1916]  3/4 folds
worst18_poolmean     18     +0.3798        +0.2115
worst18_bestpool     18     +0.3438        +0.0180
random-18 null: obs +0.1066, null mean +0.5682, CI95 [+0.4067,+0.7170],  p = 0
```

**On the tail the score cannot order its own candidate pool — its fold CI includes zero.** On the
body it orders it at +0.64. This is the mechanism behind S30-L2's filter inversion stated at the
level of the score itself rather than of the stage.

**It replicates, in degree, on both filter-independent tails** (+0.3798 and +0.3438, both well
below the 108's +0.6446). That is the opposite of S30-L16's widening effect, which did not
replicate at all — so **skill degradation is a genuine property of hard targets**, while FAIL18 is
its extreme. This is the first tail statement in this lane that survives on all three definitions,
and it should be quoted in preference to the FAIL18-only rows.

A third non-circular corroboration: **12 of 126 targets have rho_pool ≤ 0** — the score is
*anti*-informative inside its own pool — and **only 5 of the 12 are FAIL18** (2MQ2 −0.505,
2O0S −0.373, 2NB7 −0.366, 6F3V −0.288, 3SGO −0.201, 7JS6 −0.196, 1CS9 −0.147, ...).

### 2. F3b FIRES -- it is the distogram's own error, and rho_pool fully mediates

```
predictor of rho_pool          Spearman   resid(n)   fold CI              F18 mean  108 mean
distogram MAE vs native  *      -0.799     -0.819   [-0.880,-0.668]        4.367     2.001
shape error              *      -0.772     -0.797   [-0.839,-0.672]        3.630     1.884
|scale error|            *      -0.578     -0.588   [-0.664,-0.482]        3.416     1.056
|log scale ratio|        *      -0.564     -0.568   [-0.661,-0.453]        0.310     0.106
posterior sd                    -0.525     -0.546   [-0.631,-0.411]        1.667     1.411
F4 top-75 Rg sd                 -0.528     -0.530   [-0.602,-0.452]        0.464     0.353
F5 rg disagreement              -0.428     -0.405   [-0.489,-0.383]        0.448     0.094
F3 posterior entropy            -0.462     -0.456   [-0.600,-0.339]        1.468     1.435
F2 score sd                     +0.615     +0.608   [+0.525,+0.688]        0.822     0.881
F2 score gap                    -0.383     -0.391   [-0.484,-0.262]        0.402     0.346
F4 pool Rg sd                   -0.128     -0.189   [-0.316,+0.036]        1.396     1.300
F1 length n                     +0.154     +0.012   [-0.071,+0.391]       14.000    12.787
                                                        (* = ORACLE predictor)
```

The distogram is **2.18× as wrong on the tail** (MAE 4.367 vs 2.001 Å), and its error explains the
score's in-pool skill at **−0.799**, essentially unchanged when length is residualised out
(−0.819) — so this is not S28-L6's length confound.

**The chain is one mechanism, and rho_pool is the proximate variable:**

```
spearman(distogram MAE, rho_pool)                      -0.799
spearman(rho_pool, filter set-mean benefit)            +0.888
spearman(distogram MAE, filter set-mean benefit)       -0.631
   ... controlling for rho_pool                        +0.299  (sign FLIPS: fully mediated)
   ... rho_pool vs benefit controlling for MAE         +0.844  (p = 2.1e-35: survives intact)
```

Joined to S30-L16 (benefit ↔ effect of widening, +0.81) and to `operator-consumes-set-mean`
(d_out = 1.16·d_set_mean), the whole path from prior error to emitted RMSD is now measured
end to end:

> **distogram error → the score's in-pool ranking skill → the filter's set-mean benefit → the
> emitted structure.** The prior's error reaches the endpoint *through the filter's ranking*, not
> through the averaging.

That is a mechanistic account of `prior-derivative-is-the-only-steep-lever` (−2.15 Å per unit
toward a perfect prior): this entry says *by what route* the derivative acts.

### 3. F3c does NOT fire -- lane L's scale term is refuted as the mechanism

The coordinator relayed lane L's derivation that every fixed-reference pair channel contains a
separable term that is a pure function of scale, and asked whether hard pools fail for that
algebraic reason. Decomposing the distogram's error into a signed additive offset (**scale**) and
the residual after removing it (**shape**):

```
partial Spearman(rho_pool, |scale error|)  controlling shape    -0.067   p = 0.46    NULL
partial Spearman(rho_pool, shape error)    controlling scale    -0.641   p = 6.5e-16
```

**The skill collapse is shape error, not scale error.** Two direct answers to the questions put
to me:

- **Is the filter's failure predicted by the pool's Rg dispersion? No.** F4 pool Rg sd is
  **−0.128** with a fold CI **including zero** [−0.316, +0.036] — NOT MEASURED. Size dispersion
  does not predict where the filter fails.
- **Is scale error elevated on the tail? Yes, in ratio — and it still is not the mechanism.**
  |scale error| is 3.23× on FAIL18 versus shape's 1.93×, but shape is **83.1%** of the tail's
  total error (94.2% on the 108) and shape is what carries the correlation. Both facts belong in
  any quotation of this row; the elevated ratio is real and the causal claim it invites is not.

**Two caveats, both mine.** I could not reach lane L directly (a lane is not addressable from
this lane), so the scale statistic is the one **I** guessed from a second-hand relay, not the one
lane L's derivation specifies. If lane L's separable term is a different functional — a
weighted or log-domain form, or one defined against the posterior's bins rather than its
expectation — this refutation does not touch it, and I will re-run on their statistic on request.
**Treat F3c as provisional until lane L confirms the statistic.** Second, this is a correlation
over 126 targets among ORACLE quantities, not an intervention.

### 4. F3d does NOT fire -- plain wrongness, not confident wrongness

```
error alone                        -0.799
error x confidence (1/posterior sd) -0.729     gain -0.069 against a +0.10 bar
```

The interaction is **worse** than error alone. The project's "confidently wrong costs 2-3× absent"
(`torsion-restraints-reach-the-target`) does not reach this stage: what predicts the filter's
collapse is how wrong the distogram is, not how confidently. Reported as a refutation of my own
framing — I proposed "what is the distogram *confidently* wrong about" and the answer is that the
adverb is doing no work.

### 5. What this does NOT license, stated before anyone asks

The three best **native-free** correlates of rho_pool in the table — F2 score sd (+0.615),
F4 top-75 Rg sd (−0.528), F5 rg disagreement (−0.428) — are **all inside the closed router
families**, and F5 (`rg_z`) is the one S23-L7 closed at 2% of its own MDE. Correlating with an
ORACLE diagnostic at 0.5–0.6 is not the same as delivering Ångströms through a gate, which is
exactly the gap the eight closed constructions fell into, and S29-L31's incidental-parameter
result says why. **No router was built here and none should be built from this table.** What the
table is for is the opposite purpose: any future detector must report its correlation with these
columns, and a feature that lands inside them is not new.

The one thing that is genuinely changed for a future detector is S30-L16 §5's economics: the
ORACLE gate prize is −0.1193 Å at 1.29× MDE, about 2.5× the −0.0474 Å that lane C's detector was
dismissed against — so the "even a perfect detector is too small to measure" half of that closure
is specific to lane C's operator and must not be quoted generally.

Artefacts read: `s29/LEDGER.md:751` (recompute-the-posterior warning, adopted), `:3097-3151`
(S29-L33's unstratified pool-ranking table, which this entry stratifies for the first time),
`s27/LEDGER.md:1465-1542` (S28-L23b's gradient cosine −0.143/−0.016, the structure-level analogue
of §1), `:229-321` (lane C's detector and its ORACLE prize); project memory
`prior-derivative-is-the-only-steep-lever`, `operator-consumes-set-mean`,
`torsion-restraints-reach-the-target`, `error-shape-not-mae-decides-ranking`.

## S30-L14 -- THE CHARTER'S BIT QUESTION IS MIS-POSED AND THE RIGHT ANSWER INVERTS IT: **THE POOL IS A CODEBOOK, NOT A CHANNEL.** 7 INDEX BITS REALISE **36.6 BITS** OF DISPLACEMENT INFORMATION (5.2x), SO NO BITS WENT MISSING -- THE READOUT's 7 ARE WORTH FIVE TIMES THEIR FACE VALUE AND THE SYSTEM CANNOT SUPPLY ONE. PLUS THE **VALUE-OF-A-BIT LAW** (R² 0.9983), THE ENCODING ALLOCATION (CANDIDATE INDEXING WINS BY 3x; TORSION SPACE IS ARITHMETICALLY INFEASIBLE AT 7 QUBITS), AND **ONE REGISTERED PREDICTION OF MINE REFUTED** (2026-09-20 13:04, T)

**Verdict: the encoding is NOT the hidden bottleneck.** The charter named it "the least-examined
component and may be the hidden bottleneck"; it is now examined and it is not. Binary candidate
indexing is optimal among the measured classes by a factor of 3, its ORACLE floor is 1.331 A, and
the charter's target is 2.50 A -- so the allocation is not floor-limited.

### 1. THE ACCOUNTING, IN A CURRENCY DEFINED BEFORE IT WAS MEASURED

| stage | bits of CHOICE (capacity) | bits DELIVERED (measured) |
|---|---|---|
| 5 retrieval, 17,088 -> 500 | **3,252** | **+0.69** search-equivalent [SE 0.18, median +0.44] |
| 6 energy, 500 scores -> ranks | 3,767 (permutation) | -- |
| 7 top-128 prefix | 300.6 (set choice) | ORACLE ceiling 0.066 A (S29-L30); transferable part 0 |
| 9 tail readout, which of 128 | **7** | **+0.036** [MDE 0.385] (S29-L44/F2, and see S30-L4 on its median) |

**End to end the deployed system's measured information yield is under one bit per target**, against
7,019 bits of choice consumed. `b_ret = log2(N*/500)`, `N*` = the number of uniformly random draws
from the same universe whose ORACLE best matches the BLOSUM-500's, was defined in the prereg before
any array was loaded precisely so it could not be chosen afterwards.

### 2. **REFUTED: MY OWN REGISTERED P2b, BY A FACTOR OF 2 TO 4, IN THE FLATTERING DIRECTION**

I registered that the BLOSUM top-500 would beat a matched random 500 of the same universe by
**+0.10 to +0.30 A** on the ORACLE minimum, from a Gaussian-copula argument at the key's measured
`rho = +0.066` (conditioning on the top 2.9% shifts the normal-scored `rr` by `0.066 x 2.27 = 0.150
sd`; at `sd(rr) ~ 1-2 A` that is 0.15-0.30 A).

```
  M2  BLOSUM-500 best minus RANDOM-500 best (NEGATIVE = the key helps)
    a 1.7108 (med 1.7139)   b 1.7823 (med 1.8407)   n=126
    effect -0.0715   median -0.0131   SE 0.0323   MDE 0.0904   effect/MDE -0.79
    fold CI95 [-0.1474, -0.0051]   folds same sign 4/5   67W/59L   power 0.60
    VERDICT: NOT MEASURED (|effect| 0.0715 <= its own MDE 0.0904, 0.79x)
```

**Measured 0.07 against a registered 0.10-0.30.** The copula argument over-predicted by 2 to 4x --
it treats the key as a noisy copy of quality with a fixed marginal, and the real key has enormous
tie sets (an integer BLOSUM sum over 9-16 positions), so conditioning on it does not shift the
distribution the way a continuous copula does. **The direction of my error is toward my own
hypothesis** (I wanted retrieval to be delivering something measurable), which is the pattern
`unstated-operators-align-with-your-hypothesis` names, and it is the second such row today.
P2a (`b_ret <= 3.0 bits`) HELD at +0.69.

### 3. THE VALUE-OF-A-BIT LAW

The ORACLE order-statistic ladder inside the deployed pool, in the deployed order (n = 126):

     N        1      2      4      8     16     32     64    128    256    500
     D(R)  4.1080 3.6059 3.1115 2.7289 2.4901 2.2973 2.1230 1.8978 1.8061 1.7108   A
     SE    0.1810 0.1516 0.1457 0.1297 0.1205 0.1112 0.0996 0.0877 0.0844 0.0791

Fitted by `D(R) = a + c 2^(-R/gamma)` with **a = 1.3312 A, c = 2.7859, gamma = 3.1636, R² =
0.998341** (P1a HELD: bar R² >= 0.99 and gamma in [2,6]). Differentiating:

> **-dD/dR = (ln 2 / gamma) (D(R) - a) = 0.2191 (D - 1.331) A per bit.**

Marginal at R = 7: **0.1317 A/bit** (P1b HELD, registered [0.05, 0.14]; the empirical 64->128 step
is 0.2253). Fitted floor **1.3312 A** against a true universe floor of 1.3134 and a pool floor of
1.7108; `D(7) - a = 0.567 A`, so at the production register the pool is NOT near its own floor and
bits still pay (P1c HELD).

### 4. WHY NO BITS WENT MISSING -- THE CODEBOOK CROSS-CHECK

Through S29's own displacement bound (`RMSD = RMSD_prod sqrt(1-rho^2)`):

    R = 7.000 index bits: 4.1080 -> 1.8978 A  =>  rho = 0.8869  =>  **36.63 displacement bits  (5.23x)**
    R = 8.966 index bits: 4.1080 -> 1.7108 A  =>  rho = 0.9092  =>    41.55 displacement bits  (4.63x)

Bits are not conserved across the index **because the 500 deposited backbones ARE the information
and the index only names one**. The charter's premise -- the register as a channel through which
7 bits pass, 1.44 arriving -- is the wrong object. The right statement is the inverse and it is
harder, not softer:

> **the readout's 7 bits are worth five times their face value, and the system cannot supply one.**

### 5. THE ALLOCATION THEOREM, AND THE ANSWER TO "IS BINARY CANDIDATE INDEXING OPTIMAL?"

Since `-dD/dR = (ln2/gamma)(D - a)` and `D` is common at the branch point, **allocations are
comparable only through `(D - a)/gamma`: a class's ORACLE FLOOR and its TAIL INDEX. Cardinality is
irrelevant.** Ranked at the deployed width:

| allocation | A/bit | floor | note |
|---|---|---|---|
| **candidate identity** (deployed) | **0.132** at R = 7 | 1.331 A | wins |
| subset cardinality (`m`) | 0.044 | -- | S29-L44's -0.3079 A / 7; dominated 3.0x -- **and the law explains the 3.5x S29 measured** |
| mode / basin index | saturates ~1.6 bits | -- | 2-3 populated clusters; a sub-class of the first |
| **torsion / configuration** | **infeasible** | -- | see below |
| hierarchical / factorised | capped | -- | block-independent allocation cannot touch the 68% common-mode; bounded by the idiosyncratic 32% |

**The torsion arithmetic, which is a counting statement and not a preference.** Mean `n = 12.96`
over the 126, so `2n = 25.9` torsions. Naming one Ramachandran basin per residue costs
`2n log2(k)` bits: **25.9 at k = 2, 41.1 at k = 3, 51.8 at k = 4, 77.8 at k = 8.** The deployed
register is **7** qubits (`core/pipeline.py:181`) or 9 in the harness --
**0.27 bits per torsion, where 1 bit names a single basin.** Any generative proposal over torsion
space needs 4-8x the register this project has ever run, before expressivity is even discussed.
(Consistent with `phi-carries-no-sequence-signal`: the channel is 10.4 deg of psi, so the class's
`gamma` is large as well as its bit demand.)

### 6. THE rho <-> BITS DICTIONARY, SO THE SPRINT CAN ARGUE IN ONE CURRENCY

`I = -(d/2) log2(1 - rho^2)`, `d = 3n - 6 = 32.88`. Small-`rho` limit `I ~= (d/(2 ln 2)) rho^2`:
**bits are quadratic in rho**, which is lane D's "the geometry squares it" (S30-L5) -- and the
squaring is not a penalty, it is what makes `rho^2` the **additive** quantity across stages and
fields.

    rho 0.1128 (lane D's best field)   0.304 bits        rho 0.358 (3.00 A)     3.253 bits
    rho 0.140  (B2's ceiling)          0.470 bits        rho 0.628 (2.50 A)    11.895 bits

To reach 0.358 by fusing fields at 0.1128 needs `(0.358/0.1128)^2 = 10.07` **mutually orthogonal**
fields as good as the best of D's twenty-one -- the operational content of S30-L5's "10x the rho²".
**S30-L15 corrects the `d` this should be evaluated at**, and the corrected numbers are smaller.

Reproduction: `s30/PREREG_S30_T.md` (committed before any array was loaded, 7f8bbce5),
`s30/s30_T_bits.py`, `s30_T_analyse.py`, `s30_T_spec.py`,
`s30/results/s30_T_bits.json`, `s30_T_spec.json`. Every ladder row is ORACLE and none of them
chooses a deployable parameter. Derivations in `s30/THEORY.md` sections 5-7.


## S30-L15 -- **THE SECOND-MOMENT ESCAPE IS CLOSED BEFORE IT WAS BUILT** (QUADRIC IS +0.086 A WORSE THAN HALFSPACE AT MATCHED BUDGET, 2.95x MDE, 5/5 FOLDS). AND THE RESULT THAT REFRAMES THE SPRINT: **THE POOL'S TOP SIX DIRECTIONS ALREADY CONTAIN A 1.91 A POINT** -- THE ORACLE ERROR IS 8x MORE CONCENTRATED IN THEM THAN ISOTROPY PREDICTS, A 6-PARAMETER LINEAR ORDER REACHES 2.069 A OF IT, AND THE DEPLOYED PREFIX GETS 3.051 A. **THE WHOLE 0.98 A GAP IS SIX PER-TARGET COEFFICIENTS**, PRICED AT **3.78 BITS** FOR 2.50 A AND **0.855** FOR 3.00 A (2026-09-20 13:04, T)

**Verdict: the reachable set is not the obstruction, the subspace is not the obstruction, the
search is not the obstruction, and the second moment does not help. What is missing is a direction
in six dimensions, and it is smaller than anyone has priced it.**

### 1. M4 -- THE QUADRIC CEILING, BUDGET-MATCHED (n = 126, K = 5,000 each, m = 75, ORACLE)

Registered in `s30/PREREG_S30_T.md` addendum 1 before the run, with the budget matched precisely
because this is the `grid-oracles-are-order-statistics` trap. Classes, at fixed `m`:
PREFIX is **exactly one set**; HALFSPACE orders by `<g, z_x>` (6 params); QUADRIC by
`z' G z + <g, z>` (27 params); FREE is any 75-subset. Frame: all 500 superposed onto the pool
medoid ONCE, **exogenous and declared**, because a tail-dependent frame is an unstated operator.

    PREFIX    (deployed DIS top-75, one set)        3.0507 A   SE 0.1456
    HALFSPACE best-of-5,000                        2.0691 A   SE 0.0903
    QUADRIC   best-of-5,000                        2.1551 A   SE 0.0929
    FREE      best-of-5,000 random 75-subsets      2.9259 A   SE 0.1254   <- the matched null

| contrast | effect | MDE | folds | W/L | verdict |
|---|---|---|---|---|---|
| **P4a QUADRIC - HALFSPACE** | **+0.0860** | 2.95x | 5/5 [+0.068,+0.103] | 21W/105L | **WORSE** -- P4a HELD (bar: within 0.10) |
| HALFSPACE - FREE (matched null) | -0.8568 | 4.93x | 5/5 | 122W/4L | **BETTER** -- the class is real, not best-of-K |
| HALFSPACE - PREFIX | -0.9816 | 4.08x | 5/5 | 122W/4L | **BETTER** |
| P4b FREE - PREFIX | -0.1248 | 0.60x | 4/5 | 70W/56L | NOT MEASURED -- P4b HELD (bar: within 0.15) |

**HALFSPACE is the `G = 0` slice of QUADRIC**, so the deficit is a statement about SEARCHABILITY,
not containment -- 27 parameters are sparser to sample than 6 at equal budget. Operationally the
answer is the same: **at equal effort the enrichment is not findable and the halfspace ceiling is
already the binding one.** The second-moment row in S29's post-mortem, and my own T1b enrichment
(VC dim 7 -> 28, cap 50.4 -> 152.1 bits), are both **closed as ceiling-improvers** in one afternoon
instead of a lane-week. That is the second time today the pre-check has paid (see S30-L9).

**MISSED, and it is my second registered miss today.** P4c predicted every searched class would
beat PREFIX by **less than 0.6 A**; HALFSPACE beat it by **0.982 A**. Not falsified (my falsifier
was 1.0 A, and it came within 0.018 of it) but the point estimate missed by 64%, **again toward my
own hypothesis** -- I under-predicted what an ORACLE class can do because I over-weighted
`operator-consumes-set-mean`'s flattening. With P2b (S30-L14) that is two registered priors wrong
in one day, both in the flattering direction.

### 2. M5 -- THE COMBINATION CEILING, AND IT INVERTS THE QUESTION

Any native-free operator that reweights or selects pool members emits
`u = sum_x w_x (W_x - c)`, `sum w_x = 1`, so it lies in the span of the pool's deviations about the
production average `c`. Project the ORACLE error `e = t - c` onto the top-`k` principal directions
of that deviation matrix (n = 126; deviation-matrix stable rank 2.717, PC1 share 0.391):

     k      ||Pi_k e||^2/||e||^2    rho_ceiling    isotropic null k/d
     1            0.1997              0.4469           0.0257
     2            0.3514              0.5928           0.0514
     3            0.4369              0.6610           0.0772
     6            0.6072            **0.7792**         0.1543
    10            0.7519              0.8671           0.2572
    21            0.9261              0.9624           0.5401

**The oracle error is 8x more concentrated in the pool's leading deviation directions than isotropy
predicts.** Six directions contain a `rho = 0.779` point -- **1.91 A** through the bound, against
production's 3.051. **The subspace is NOT the obstruction**, which is the opposite of what I
expected and of what the 68% common-mode result had led me to expect. (See the declared defect box
in the prereg: my clause P5d was written AFTER this ran and is **not** a pre-registration; its
failure is recorded so the wrong prior is auditable, and carries no evidential weight.)

**M4 and M5 agree to 0.16 A from two independent instruments**: the 6-direction subspace ceiling is
1.91 A and the 6-parameter halfspace class reaches 2.069 A. A linear order nearly saturates the
subspace that contains the answer.

### 3. THE SPRINT's PROBLEM, WITH A DIMENSION AND A BIT PRICE ON IT

The reachable displacement space contains a **1.91 A** point inside **six** well-conditioned
directions; a 6-parameter linear tail order reaches **2.069**; the deployed prefix gets **3.051**.
**The entire 0.98 A gap is the problem of locating six per-target coefficients** -- S29's
Neyman-Scott incidental parameter with a dimension attached for the first time.

**And the price is much lower than the full-space dictionary says**, because the subspace itself is
FREE: the pool's PCA is native-free and per-target. Evaluating `I = -(d/2) log2(1 - rho^2)` inside
the top-6 subspace, where a field needs only `cos = rho_target / 0.779` against the true direction
(two derivations, sphere rate-distortion and the covering number of `S^5` by caps, agreeing
exactly):

    rho_target              FULL SPACE (d = 32.88)     WITHIN THE POOL's TOP-6
    0.1128 (D's best field)      0.304 bits                 0.076 bits
    0.140  (B2's ceiling)        0.470                      0.118
    0.358  (3.00 A)              3.253                    **0.855**
    0.628  (2.50 A, charter)    11.895                    **3.782**
    0.695  (2.31 A)             15.649                      5.726

> **The charter's 2.50 A needs 3.78 bits of per-target direction information in a subspace the
> pipeline can already compute for nothing. 3.00 A needs 0.855 bits. The best native-free field on
> record supplies 0.076.**

**CORRECTION TO MY OWN UNREGISTERED ASIDE**, flagged loudly because it is the kind that gets quoted
without a bar: I told the coordinator that "six coefficients at ~2 bits each = 12 bits" matched the
full-space 11.895 and called it two independent routes. **It was a coincidence and it is
withdrawn.** The two columns above are the same formula evaluated at two different `d`, and the
8.1-bit difference at 2.50 A is exactly the value of knowing the subspace. The corrected number is
3.78, not 12, and it makes the missing thing **smaller**, not larger -- which sharpens the question
from "is there enough signal anywhere" to "why can nothing supply four bits of direction in a
subspace we already have".

The factor to 3.00 A is `0.855/0.076 = 11.2`, which is S30-L14's "10x the rho²" **in the additive
currency and from the same algebra** -- one witness, not two, and it must not be reported as two.

### 3b. PRIORITY AND A CORRECTION TO MY OWN SECTION 1: **LANE Q GOT THERE FIRST AND WITH A BETTER NULL**

`S30-L12` (Q, 13:01) measured the quadric class **before this entry was posted** and its design is
stronger than mine in the one place that matters. I am recording that rather than letting two
entries claim the same closure.

- **Priority.** The quadric verdict is lane Q's. Its numbers (128 directions/class/target, all
  `m = 1..500`, 126 targets): QUADRIC - LINEAR is **-0.0105 A at 0.16x MDE (NOT MEASURED)** at
  ORACLE `m` and **+0.2059 A at 1.81x MDE, 37W/89L (WORSE)** at the shipped `m = 75`. My M4 is an
  independent sampler (K = 5,000 inside a 6-dimensional PC subspace, fixed `M = 75`) and agrees in
  sign and roughly in size at `m = 75` (+0.0860, 2.95x MDE, 105L/21W). **Two samplers, two lanes,
  same sign** -- but the first statement of it is theirs. Lane Q cites my *two-target smoke file*
  as corroboration; the completed 126-target run is in
  `s30/results/s30_T_quadric.json` and the full-n numbers are the ones above.

- **CORRECTION TO MY SECTION 1.** I wrote "HALFSPACE - FREE ... the class is real, not best-of-K".
  That is true against the null I registered -- a matched best-of-5,000 over **random 75-subsets**
  -- but that null is too weak, and lane Q's is the right one. Against an **across-target** null
  (does the direction found for target *i* transfer to target *j*?) lane Q measures the linear
  class's entire apparent gain as **225% accounted for**, with **split-half transfer +0.0323 A
  (-3%)**. **So the 0.98 A halfspace gain in my section 1 is per-target and does NOT transfer**,
  and my "the class is real" must be read as "the class is reachable", not "the class is
  exploitable". I did not run a transfer arm and should have; `grid-oracles-are-order-statistics`
  says W/L cannot diagnose this and only split-half transfer can.

- **What this does NOT touch: M5.** Section 2 is a **projection**, not a maximum over draws --
  there is no best-of-K anywhere in it, so the across-target-null critique does not apply. The
  `rho_ceiling = 0.779` at `k = 6` and the 1.91 A subspace figure stand as written. And lane Q's
  transfer result and my section 3 are **the same finding from two directions**: they measure that
  no direction transfers, and I measure that the missing quantity is six per-target coefficients.
  A non-transferable direction in a 6-dimensional subspace IS six per-target coefficients.

### 3c. **RETRACTION OF MY OWN SECTION 1 HEADLINE, ON MY OWN K = 5,000** (M6)

Lane Q asked me to price my search as the best-of-K it is, the coordinator endorsed it, and the
answer retracts my reading. `s30/s30_T_transfer.py`, n = 126, K = 5,000, with a **common direction
bank** across targets (columns made comparable by a deterministic native-free sign convention --
each PC's largest-|loading| entry positive -- so column k is the same RULE on every target and
`best_of_k_within` is well posed):

```
                      observed   across-target null   accounted   split-half transfer   k_eff
  HALFSPACE  K=5000   -1.7881         -3.5038           196%        -0.3319 (19%)       118
  QUADRIC    K=5000   -1.5663         -2.9368           188%        -0.3148 (20%)       110
  verdict: NOT A SIGNAL in both (transfer < 25% of the oracle)
```

The across-target null **exceeds** the observed gain, exactly as lane Q predicted for my larger K.
And because a 19% transfer is not zero and could be misread as a route, the decisive check:

    PREFIX (deployed)                            3.0507 A
    HALFSPACE mean over the 5,000 directions     3.8550 A
    the transferable rule (K-mean + split-half)  3.5231 A
    vs PREFIX  +0.4724   1.88x MDE   5/5 folds   39W/87L   ->  **WORSE**

> **RETRACTED: "HALFSPACE - PREFIX -0.9816 A ... BETTER" in section 1 above is withdrawn AS A
> STATEMENT ABOUT THE CLASS's VALUE.** It is a per-target best-of-5,000 and 196% of it is an
> across-target order statistic; the class's transferable content lands **0.472 A worse than what
> ships**. **The halfspace class is reachable, not exploitable, and as a rule it is negative.**
> My matched random-subset null was too weak: it asks whether a structured class beats an
> unstructured one at equal budget, not whether the winning direction is the same direction twice.

**Unaffected:** (i) the QUADRIC - HALFSPACE contrast, a like-for-like comparison at matched budget
that a null shifting both does not touch; (ii) section 2's M5, which is a projection with no
maximum over draws anywhere in it; (iii) section 3's conclusion, which this **strengthens** -- a
non-transferable direction in a six-dimensional subspace IS six per-target coefficients.

### 3d. SCOPE ON THE STABLE RANK, ADOPTED FROM LANE Q IN FULL

Lane Q verified 1.859 independently through the DEPLOYED path (`s29_O_ladder.load_pool` +
`s28_A_amp.Frame`, 12 pools): 1.862 mean, PC1 54.4%, k90 5.6. **But the collapse is a property of
the PAIR-DISTANCE feature space, not of the pool.** In COORDINATE space it is **3.404** (mine) /
**3.619** (Q's), with **k90 = 11.2 directions** -- comfortably above my P3b threshold of 2.0. The
sparse/weighted readout and every second-moment construction act on **coordinates**.

> **"stable rank 1.86" must never appear without its feature space in the same sentence.** It is
> the right number for a distance-map lift and an overstatement for a coordinate-space one, and the
> two differ by ~2x in the direction that matters. At the coordinate figure the lift is still only
> an 11-parameter family (cap ~80 bits against 300.6), so the VERDICT is unchanged -- but it now
> rests on the direct measurements in section 1 and 3c, not on the rank threshold.

### 4. WHAT THIS DOES AND DOES NOT CLOSE

- **CLOSES** the second-moment / quadric tail class as a ceiling-improver (M4), and with S30-L9's
  T1b it closes the tail-then-aggregate lift as a whole.
- **CLOSES** "the reachable set is too poor" and "the span is too small" as explanations. Neither
  is true: 1.91 A is inside six directions the pool hands you.
- **DOES NOT** close, or even touch, whether anything native-free can supply 0.86 to 3.78 bits of
  direction. Every arm here is ORACLE. No native-free rule is proposed, tested or recommended.
- **DOES NOT** contradict the bound. It restates where the bound bites: not in the operator space,
  not in the pool, not in the encoding -- in six coefficients.

### 5. THE SCOPE CAUTION, BECAUSE THREE LANES NOW POINT AT DIFFERENT COMPONENTS

This entry is about the **operator** (which displacement it emits) and S30-L14 is about the
**encoding** (how the register indexes candidates). Lane F's S30-L2 is about a **stage** (the
distogram's 500 -> 75 filter, +1.767 A worse than a random 75 on FAIL18). These are three distinct
components and the report must keep them distinct. "The encoding is not the bottleneck" is a
closure about the encoding **among the measured classes**; it is not, and must not be read as, a
claim that nothing upstream or downstream is.

Reproduction: `s30/PREREG_S30_T.md` addenda 1-2, `s30/s30_T_quadric.py`, `s30/s30_T_combo.py`,
`s30/results/s30_T_quadric.json`, `s30_T_combo.json`. Derivations in `s30/THEORY.md` sections 4, 8.


> **NUMBERING COLLISION RESOLVED 2026-09-20 13:29 by the coordinator.** Posted as `S30-L6` at 2026-09-20 13:17 by lane D, which an entry stamped earlier already held. Renumbered to **S30-L21**; the earlier entry keeps `S30-L6`. Recorded rather than silently fixed (contract rule 15). Anything citing `S30-L6` for *this* result means `S30-L21`.

## S30-L21 -- **THE COMBINATION QUESTION IS CLOSED, AND THE REASON IS RANK, NOT COUNT.** THE 21 FIELDS HAVE A PER-TARGET STABLE RANK OF **1.68** -- ELEVEN SIGNIFICANT FIELDS ARE ~2 DIRECTIONS. THE **ORACLE GLOBAL** WEIGHTING, ONE w FOR ALL TARGETS WITH THE NATIVE IN HAND, REACHES ρ = **0.1693 = 0.69 BITS AGAINST THE 3.22 NEEDED**, WORTH **0.046 Å ON THE BUILT CHAIN**. EVERY NATIVE-FREE ARM **LOSES TO THE SINGLE BEST FIELD**, AND THE 21-PARAMETER FITTED ARM IS THE WORST OF THEM (ρ = 0.012). MY OWN PRE-REGISTERED FALSIFIER IS NOT MET (2026-09-20 13:17, D)

**Verdict: no. Combining the fields does not reopen the bound. Prereg
`s30/PREREG_S30_D_gram.md` predicted this and stated the falsifier; the falsifier failed.**

### The question (coordinator, after S30-L5)

If eleven native-free displacement fields each carry real signal at ρ ≈ 0.10–0.11, and orthogonal
fields would combine to √(Σρᵢ²) ≈ **0.332–0.365**, and 3.00 Å needs **ρ = 0.3561** — does the
combination clear the threshold? The arithmetic lands close enough that it had to be measured.

### The answer, on the built chain (production 3.2105 Å, ρ for 3.00 Å = 0.3561 = 3.22 bits)

| arm | ρ | bits | built-chain Å | vs production |
|---|---|---|---|---|
| **what is needed for 3.00 Å** | **0.3561** | **3.218** | 3.0000 | −0.2105 |
| ORACLE per-target weighting (21-dim projection) | 0.9491 | 54.8 | 1.0111 | −2.199 |
| — **matched random 21-dim subspace (the control)** | **0.8095** | **25.3** | 1.8850 | −1.326 |
| **ORACLE GLOBAL weighting** (one w, native in hand) | **0.1693** | **0.690** | **3.1642** | **−0.046** |
| best single field (CHAN_DISTPOT), LFO-chosen | 0.1214 | 0.352 | 3.1867 | −0.024 |
| EQ11, zero-parameter, leaked selection | 0.1139 | 0.310 | 3.1896 | −0.021 |
| **LFO-selected equal weights (no leakage)** | **0.0948** | **0.214** | 3.1960 | −0.015 |
| LFO global weighting, 21 fitted parameters | 0.0124 | 0.004 | 3.2103 | −0.000 |

n = 119 of 126 (7 targets lack one of the 9 S27 channels, so their field set is not the same 21
and they are excluded from every aggregate rather than silently padded). Cosines are POINT-CLOUD;
the Å column carries them through the bound with the published built-chain production 3.2105 Å,
which is the basis the record's 0.3561 threshold lives on. Both bases are printed by the tool.

### Why: rank, exactly as the coordinator guessed

**The Gram of the 21 unit field directions has a per-target stable rank (trace/λmax) of 1.681**
and an effective rank of 3.68; aggregated, stable rank 2.057 with the **top 3 eigenvalues
carrying 60.0%** and λ₁ alone carrying 10.21 of 21. Mean |off-diagonal| 0.427. **Eleven
significant fields are about two directions.** This is the same collapse lane T measured on the
pair-distance matrix (stable rank 1.859) — the fields are all differences of averages over
overlapping prefixes of one pool, ranked by scores that all contain DIS, so they cannot be
independent and are not.

The orthogonality arithmetic that motivated the question assumed 11 independent directions.
There are ~2. **√(Σρᵢ²) = 0.33 was the right formula applied to the wrong rank.**

### The one real thing here, and why it is still unusable

The ORACLE per-target projection reaches ρ = 0.9491. **Read alone that is nonsense**: a random
21-dimensional subspace of these targets' 32.9-dimensional rigid-body-removed space reaches
0.8095 by dimension counting. But the excess is real and large — **+0.1396, SE 0.0073, 6.83× MDE,
fold CI [+0.132, +0.148], 5/5 folds** — and the rank curve shows it is not an artefact of the
near-degenerate tail: at rank 4 the real field set reaches **0.6627 against 0.3299 for a random
4-dimensional subspace (+0.333, 6.31× MDE)**. The field class's span genuinely points at the
native better than a random subspace of the same size.

**And it is per-target.** The gap between the per-target ORACLE (0.9491) and the global ORACLE
(0.1693) — both with the native in hand, differing only in whether the weights may vary by target
— is the entire result. This is `in-band-ordering-is-per-target` appearing in a new basis: 0.949
within a target, 0.169 across them. S29's incidental-parameter result says per-target weights are
not estimable from other targets, and the four native-free arms confirm it from the other side.

### Three things I had to fix in my own arms before the numbers meant anything

1. **My first LFO fit returned a NEGATIVE cosine (−0.218).** An untuned ridge on a stable-rank-2
   Gram is a strawman. The ridge is now chosen by nested leave-one-fold-out *inside the training
   folds*. It selects the largest value in the grid on 4 of 5 folds — i.e. the fitted covariance
   contributes nothing and the rule degenerates to marginal-cosine weighting.
2. **My "ORACLE global" arm was fitted by least squares and returned ρ BELOW the best single
   field** — impossible for a true maximum, since w = eⱼ reproduces field j exactly. That internal
   consistency floor is what caught it. The arm now maximises the mean cosine directly (analytic
   gradient, 22 starts including every single-field basis vector). It is the reason the reported
   0.1693 can be trusted as a ceiling rather than an artefact of an objective mismatch.
3. **My EQ11 arm selected its 11 fields using the same 126 targets it was scored on.** Replaced by
   an arm that selects fields inside the training folds only. It costs 0.019 in ρ (0.1139 →
   0.0948), which is the size of that leakage and is worth knowing.

### Against my own pre-registered falsifier

> the rank curve's excess is materially positive at some rank **AND** the leave-fold-out global
> weighting reaches ρ ≥ 0.20 with a fold CI excluding the best single field's ρ.

First clause **met**. Second clause **failed on both counts**: the best honest native-free arm
reaches ρ = 0.0948, less than half of 0.20, and against the best single field it is
**−0.0267, fold CI [−0.067, +0.010], −0.36× MDE**. The fitted 21-parameter arm is worse still
(−0.1091, −1.31× MDE, 5/5 folds against it). Under the sprint's rule, nothing here is a result
except the negative: **below 0.7× MDE is not a result, and every combination arm is below it.**

Concentration was checked against a uniform-effect null (52nd percentile, no flag), so the small
positive-vs-zero effects are not one or two targets.

### What this settles for the sprint

- **Do not spend compute searching for a better weighting of this field class.** The ceiling with
  the native in hand is 0.046 Å on the built chain; the deployable arms return 0.015 Å and lose to
  picking one field.
- **The barrier is the rank of the information, not the number of channels.** Adding a 22nd field
  built the same way — from the same pool, ranked by a score containing DIS — will land in the
  same 2-dimensional span. A field that helps must be *orthogonal* to that span, and the Gram
  above is the instrument for checking that before anything is built: it is now a one-line test.
- **B2 stands, and stands more firmly than before.** Its number survives (S30-L5), its argument is
  replaced (S30-L5), and now its scope is measured: the bound does not move under combination.

Artefacts: `s30/s30_D_gram.py`, `s30/results/s30_D_gram.json`, `s30/results/s30_D_gram/*.npz`,
prereg `s30/PREREG_S30_D_gram.md`. The regenerated fields are asserted equal to
`s29/results/s29_D_fields_rows.jsonl` cosines to 1e-9, so this prices the same 21 fields S30-L5 did.
Multiplicity: 25 comparisons from this file (21 of them the rank curve).

> **NUMBERING COLLISION RESOLVED 2026-09-20 13:29 by the coordinator.** Posted as `S30-L7` at 2026-09-20 13:24 by lane D, which an entry stamped earlier already held. Renumbered to **S30-L22**; the earlier entry keeps `S30-L7`. Recorded rather than silently fixed (contract rule 15). Anything citing `S30-L7` for *this* result means `S30-L22`.

## S30-L22 -- **LANE L's PREDICTION IS CONFIRMED, AND HARDER THAN IT ASKED FOR.** THE RADIAL/SCALE DIRECTION CARRIES **58.0% OF THE GRAM TRACE**, AND THE FIELD SET'S DOMINANT PRINCIPAL DIRECTION **IS** THE RADIAL ONE AT cos = **0.947** (MEDIAN 0.984). REMOVING IT RAISES THE STABLE RANK 1.705 -> 2.642. AND THE PAYOFF: THE DIRECTION TO THE NATIVE IS **ORTHOGONAL** TO IT (cos −0.068), **−0.252 ON FAIL18** -- THE LIBRARY SPENDS THE MAJORITY OF ITS RANK ON THE ONE COMPONENT THAT POINTS ELSEWHERE (2026-09-20 13:24, D)

**Verdict: lane L's algebra is validated on an independent measurement. Its consequence is worse
than "one of the two directions is spent on scale" -- that direction is anti-aligned with the
native on exactly the targets that matter.**

Lane L predicted that a fixed-reference distance potential's separable scale term would make one
of S30-L6's ~2 effective directions the radial field `x − centroid`. Test, on the 126 cached
field sets (`s30/results/s30_D_gram/`), rigid body removed from both sides as everywhere else:

| quantity | mean | median | fold CI | FAIL18 | other 108 |
|---|---|---|---|---|---|
| **radial share of the Gram trace** | **0.5798** | 0.5879 | [+0.545, +0.603] | 0.5528 | 0.5843 |
| λ₁ share of the trace | 0.6117 | 0.6187 | [+0.588, +0.629] | 0.5777 | 0.6173 |
| **cos(dominant principal direction, radial)** | **0.9472** | **0.9840** | [+0.924, +0.966] | 0.9469 | 0.9473 |
| cos(2nd principal direction, radial) | 0.1502 | 0.0952 | [+0.126, +0.187] | | |
| stable rank, as measured | 1.7050 | 1.6163 | [+1.64, +1.78] | 1.8515 | 1.6806 |
| **stable rank, radial removed** | **2.6425** | 2.5467 | [+2.55, +2.75] | 2.7923 | 2.6176 |
| λ₁ share, radial removed | 0.3962 | | [+0.384, +0.407] | | |
| mean \|field·radial\| over the 21 fields | 0.6908 | 0.7155 | [+0.664, +0.709] | | |
| **cos(direction to the NATIVE, radial)** | **−0.0675** | −0.0575 | [−0.126, −0.011] | **−0.2524** | −0.0367 |

**Confirmed.** Lane L asked whether the radial direction captures "roughly half the spectrum": it
captures **58.0%**, and **94.8% of the top eigenvalue is radial** (0.5798 / 0.6117). The
rank-1.7 appearance of the field library is the scale direction: take it out and the stable rank
rises to 2.64 while λ₁'s share falls from 61% to 40%.

**The consequence, which is the reason it was worth four lines.** The last row is ORACLE and it is
the one that matters. The direction to the native has a radial cosine of **−0.068** overall — the
field library concentrates the majority of its two available directions on a component that is
essentially orthogonal to where the native lies. On FAIL18 it is **−0.2524**: on the hard targets
the radial direction is not merely useless, it is **anti-aligned**, and a method that moves along
it moves away from the answer. That is consistent with lane F's finding that shape, not scale,
carries the tail's error, and it supplies the mechanism: the library cannot express shape because
its rank is spent on scale.

**What this does NOT show.** It does not show that removing the scale component would help. The
residual 42% of the spectrum is spread over a stable rank of 2.64 with no large eigenvalue, and
S30-L6 already measured what the whole 21-dimensional span is worth through a global weighting
(ρ = 0.1693, 0.046 Å on the built chain). Deflating the radial direction reallocates rank; it does
not create any. A lane wanting to act on this must show the deflated set reaches a higher ρ, and
the Gram tool is now a one-liner for checking a proposed new field's orthogonality to this
direction **before** it is built.

Reproduction: `s30/results/s30_D_gram/*.npz`, radial direction `remove_rigid(C0 − mean(C0), C0)`.

---

> **NUMBERING COLLISION RESOLVED 2026-09-20 13:29 by the coordinator.** Posted as `S30-L8` at 2026-09-20 13:24 by lane D, which an entry stamped earlier already held. Renumbered to **S30-L23**; the earlier entry keeps `S30-L8`. Recorded rather than silently fixed (contract rule 15). Anything citing `S30-L8` for *this* result means `S30-L23`.

## S30-L23 -- **S30-L2's HEADLINE IS CONDITIONED ON ITS OWN NUMERATOR.** FAIL18 IS *DEFINED* IN `s12/instrument.py:277` AS THE TARGETS WHERE THE SCORE'S TOP-75 RETAINS **ZERO** POOL MEMBERS WITHIN 1.5 Å OF THE POOL OPTIMUM. I REPRODUCE THE +1.767 Å EXACTLY (−1.7622, OPPOSITE SIGN CONVENTION) AND AT LEAST **0.869 Å OF IT IS FORCED ARITHMETIC**; THE STRATUM keep=1 SHOWS **−0.036**. OUTSIDE THE 18, THE FILTER IS **+0.025 Å**, AND THE ALL-126 EFFECT IS **100% THE 18** (2026-09-20 13:24, D)

**Verdict: F1a stands untouched. F1c's FAIL18 row must be withdrawn as an effect estimate. Lane F
did the hard parts right -- prereg before the numbers, matched random subsets in each operator's
own space, fold CIs, a random-18 null -- and none of those defences reaches this one, because the
stratum itself is the outcome.**

### The definition

`s12/instrument.py:271-278`, the selfcheck that pins the constant:

```python
sub  = np.asarray(rec["sub"], int)              # the SCORE's top-75
band = np.where(rr <= rr.min() + BAND)[0]       # pool members within BAND = 1.5 A of the pool best (ORACLE)
if not np.isin(band, sub).any():
    zero.append(t["pdb"])
assert set(zero) == set(FAIL18)
```

**FAIL18 is the set of targets on which the 500 -> 75 filter retained no in-band member.** S30-L2
then measures, on that set, how much worse the filter's retained set is than a random 75 in
ORACLE best. The selection predicate is a lower bound on the measured quantity's first term.

### Reproduction and decomposition (n = 126, 200 random-75 draws per target, seed 3030)

I reproduce lane F's number independently: (ORACLE best of random-75) − (ORACLE best of top-75)
= **−1.7622 Å** on FAIL18 against its +1.7674 in the opposite sign convention.

| stratum, by `keep` = in-band members surviving the filter | n | effect |
|---|---|---|
| **keep = 0 — this stratum *is* FAIL18** | 18 | **−1.7622** |
| keep = 1 | 3 | −0.0363 |
| keep = 2–3 | 1 | −0.8654 |
| keep = 4–8 | 8 | −0.0039 |
| keep ≥ 9 | 96 | +0.0380 |

**There is no gradient.** The effect is a step at the selection boundary, not a continuum in
recall: one retained in-band member is enough to remove 98% of it. Among the 108 non-FAIL18
targets, Spearman(keep, effect) = **−0.046**.

**How much is forced.** By definition `top_best − pool_best ≥ 1.5` on these 18 (measured: 2.393).
The random arm is unconstrained and measures `rand_best − pool_best = 0.631`. So the selection
predicate alone forces |effect| ≥ 1.5 − 0.631 = **0.869 Å, 49% of the headline**. The remaining
0.893 Å is not an unbiased estimate either: it is the expected *exceedance above a selection
threshold*, which truncation inflates by an amount this design cannot measure.

**Difficulty is not the explanation, and it is not the fix.** The 18 targets with the worst ORACLE
pool best — a difficulty criterion **upstream** of the filter — give **−0.629 Å**, 36% of the
headline, and only 8 of them are in FAIL18 (so even that is contaminated).

**And the all-126 row is the 18.** Effect on all 126 = −0.2308; 18/126 × (−1.7622) = **−0.2517**.
Outside the 18 the filter is **+0.0245 Å** — *better* than a random 75, not worse. Lane F reported
that row honestly ("NOT MEASURED", CI spanning zero); the reading that must not survive is
"the harmful thing the pipeline does on hard targets is located in one stage."

### What survives, and it is not nothing

- **F1a is untouched**: the ORACLE best member of the FAIL18 pools is 2.284 Å. That is a statement
  about pools, involves no conditioning, and closes the pool-limited branch. It stands.
- **The existence of 18 zero-recall targets is real** — but it is a *known* property, the one S12
  used to define the set, not a new measurement of filter damage.
- **What is not established** is that the filter is differentially harmful on hard targets. The
  honest test defines the stratum on information the measurement does not reuse: split the pool,
  define recall on one half and measure ORACLE best on the other; or stratify on an upstream
  difficulty criterion, which gives −0.629 Å and would have to be priced against its own null.

### The rule this earns

`control-must-match-the-operators-space` is instance 5 of this sprint's most repeated error, but
this is its other face and it needs its own name: **a matched control in the right space does not
rescue a stratum defined by the outcome.** Lane F's random-75 control was drawn in exactly the
right space and is not the problem. The problem is that the 18 rows it was averaged over were
chosen because the quantity being measured was large on them.

Reproduction: `s12/instrument.py:271-278`; my run reads `u["rr"][pool_idx]` and
`shipped_record(pdb)["sub"]` only, no distogram, seed 3030, 200 draws per target.

## S30-L18 -- FOR LANE P, BEFORE IT BUILDS: ITS TARGET IS **ALREADY LOCALISED** (A PER-TARGET SEPARATION PROFILE, FIVE NUMBERS, WORTH **0.525 A ORACLE**, WITH **19-30% ALREADY RECOVERED NATIVE-FREE**) AND ONE OF ITS OBVIOUS ROUTES CARRIES AN EXPLICIT "NOBODY SHOULD SPEND ON THIS" FROM S19. PLUS TWO DERIVATIONS: **ORTHOGONALITY IS WORTH A 5.1% DISCOUNT ON THE REQUIREMENT, NOT A ROUTE** -- DECORRELATION IS NOT THE LEVER, SKILL IS -- AND LANE F's ORACLE 0.169 IS **EXACTLY WHAT A ZERO-SKILL PARTNER AT c = 0.745 PRODUCES** (2026-09-20 13:25, L)
Note: `s30/lit/L30_4_shared_bias.md`. Derivation script: `s30/lit/s30_L_orthogonality.py`
(closed form; reads no project data). **Reading, derivation and arithmetic; no new measurement.**

Commissioned by the coordinator on spawning lane P. Everything below was found by reading the
record **before** searching outside it, which is how the first item was caught.

### 1. THE RE-IMPORT WARNING, AND IT IS EXPLICIT IN THE RECORD

I developed a proposal to project the distogram's predicted distance matrix onto the rank-3 EDM
cone -- the constraint is real (for n=13, a symmetric hollow matrix has 78 free entries and a
realisable 3-D structure has 3n-6 = 33, so **58% of the prediction's freedom is metric
inconsistency**), it is parameter-free, and it uses no new information. **Then I grepped the
codebase, and S19 has already measured it and left a named prohibition**
(`s19/agentA_FINDINGS.md` section 1.2, `s19/a_topo.py`, `s19/a_coh.py`):

| quantity | predicted field | native |
|---|---|---|
| triangle-inequality violations | 4.09 % | 0.02 % |
| rank-3 EDM defect of `-1/2 J D^2 J` | **0.286** | 0.0018 |

> *"**Independent pairwise marginals are not jointly realisable** is real, measured, and **NOT the
> mechanism**. Nobody should spend on EDM projection, triangle repair, or joint-consistency
> enforcement as a route to RMSD."*

**A4 was refuted in the direction nobody expected**: a magnitude-matched *incoherent* field is
WORSE on both realisability measures (0.340 defect, 10.85% violations) and lands **1.24 A
BETTER**. Unrealisability is anti-correlated with harm. And `rho(EDM defect, RMSD) = +0.511`
against `rho(residual RMS, RMSD) = +0.891`, so the defect adds nothing over error size -- **it is a
magnitude proxy**. My 58%-of-freedom arithmetic is correct and is beside the point.

**Third instance of this lane catching itself by reading the source** (after CAGEO in S30-L6 and
the near-miss S29-L17 recorded). I am logging it because the hit rate matters: three proposals
that survived the information test died on an existence-or-precedent check.

### 2. THE FIVE-INSTANCE LAW THIS COMPLETES, AND IT IS THE SPRINT'S MOST TRANSFERABLE OUTPUT

S19's result is not isolated. Assembled across the record and the literature:

| # | operator | realism moved | accuracy moved | source |
|---|---|---|---|---|
| 1 | make the distance field realisable | EDM defect 0.286 -> **0.148**, kappa -> 0.914 | **more coherently wrong** | S19 sec 1.2 / 2 |
| 2 | project the cloud to geometric consistency | fit residual 0.8135 -> **0.6191 A** | **structure worse**, ~+0.73 A vs production | S29-L57 |
| 3 | repair averaging artefacts (MCORE) | clashes **63.0% -> 1.09%** of atoms | **3.28 -> 3.36 A**, 0.08 A worse | PMC2662860, n=2090 |
| 4 | calibrate the posterior | calibration improves | **RMSD worse** | S25 L2 |
| 5 | AMBER-relax the average, k=30 | valid geometry | **-0.022 A**, 5/5 folds | `averaging-space-beats-the-objective` |

> **Every realism-enforcement operator this project or the literature has measured improves
> realism and costs, or fails to help, accuracy. Five instances, two of them at n >= 126 and one
> at n = 2090. There is exactly one exception, row 5, and S30-L13 gives its size a mechanism:
> a constraint repair acts on the bonded geometry, whose share of all residue pairs is `2/n` --
> 15.4% at n=13 against 1.0% at n=200 -- so it works here because the chain is short and it does
> not scale with effort.**

The theorem underneath is Blau & Michaeli (CVPR 2018) Thm 3, imported in S29-L12 and applied there
to scorers only. **It applies to OPERATORS too, and rows 1-4 are that statement measured four
times on four different objects.** Any future proposal of the form "make the intermediate object
more valid / consistent / calibrated / physical" should be priced against this table before it is
staffed.

### 3. LANE P's TARGET IS ALREADY LOCALISED -- BUILD ON THIS, DO NOT REDISCOVER IT

`s19/agentA_FINDINGS.md` sections 4, 5 and 6, all ORACLE, n = 126:

- **The harm is COHERENCE, not magnitude and not unrealisability.** A perfectly coherent
  matched-magnitude error is the **worst arm on the board** (3.749 vs real 3.610, +0.139
  [+0.023, +0.253]). Destroying only cross-pair sign coherence, at exactly matched magnitude,
  buys **-1.202 A [-1.408, -1.004], 112W/14L, 5/5 folds**.
- **The harmful mode is named: a per-target separation profile -- FIVE NUMBERS PER TARGET.**
  Nested variance explained: offset (1 param) 0.141, stretch (1) 0.226, separation profile (5)
  0.361, per-residue additive (13) 0.430. But **variance explained does not price damage**:
  removing the offset *hurts* (+0.174), removing the stretch *hurts* (+0.069), the per-residue
  additive explains the most variance and owns **none** of the harm, and the separation profile
  owns **+52.5% of the gap, worth -0.525 A** with the native in hand.
- **The wall.** The best native-free estimator of that profile recovers **0.098-0.155 A of the
  0.525 A -- 19 to 30% -- and only by regressing toward what the pipeline already does.**
  S19's verdict: *"There is still no native-free estimator of coherence."*
- **And the bias is SHARED across independent sources** (S19 sec 4: *"the harmful coherent
  component is shared across independent [sources]"*; the incoherent component does not share at
  all, 0.09-0.19 cross-family).

**That last bullet is my S30-L7 escape E3, measured in S19 two sprints before I derived it, and
the S24 provenance cosine (0.9432 against a 0.9330 within-source control) is a third independent
arrival at the same fact.** Three arrivals: a theorem, a cross-family correlation measurement, and
a provenance cosine.

**And the five numbers are five INCIDENTAL parameters per target** (Neyman & Scott 1948, S29-L31).
That composes with the above into the sharpest brief I can give lane P:

> The prize is real and large (0.525 A ORACLE, the biggest in the record outside the prior
> itself), the object is small and named (five numbers per target), 19-30% is already taken, and
> the remainder is an incidental parameter. **S29-L31 names exactly two escapes: replication
> within the instance, or a covariate OBSERVED AT INFERENCE. So lane P's only viable form is a
> native-free covariate, observed at inference, that predicts the per-target separation profile
> AND is not generated by the distogram or by the retrieval pool** -- because both share the bias
> being estimated. That is a much narrower target than "find a decorrelated source", and it is
> falsifiable at design time rather than after a run.

### 4. DERIVATION: WHAT IS ORTHOGONALITY ACTUALLY WORTH? (`s30_L_orthogonality.py`)

The brief asked for sources whose errors are *decorrelated* from the prior. I priced that request
in closed form. With `rho_max^2 = (r1^2 - 2 c r1 r2 + r2^2)/(1 - c^2)` (the projection of `u` onto
`span{v1,v2}`) and the project's anchors (best single field 0.1128; 0.358 needed for 3.00 A):

```
  to reach 3.00 A combined with our best field, a PERFECTLY ORTHOGONAL new channel
  must itself carry rho = 0.3398.  Alone it would need 0.3580.
  -> ORTHOGONALITY IS WORTH A 5.1% DISCOUNT ON THE REQUIREMENT (1.6% for 2.50 A)
  -> the new channel must be 3.01x better than anything we own

  gain is QUADRATIC in the NEW channel's own skill (c = 0: rho_max ~ r1 + r2^2/(2 r1)):
     r2 = 0.1128 (as good as everything we own, perfectly orthogonal)  ->  -0.0206 A
     r2 = 0.1398 (the random-shape reference)                          ->  -0.0317 A
     r2 = 0.34                                                         ->  -0.1926 A
  reaching 3.00 A by stacking needs k = 10.1 MUTUALLY ORTHOGONAL channels each as good
  as our best.  Lane F measured the 21 we have at Gram stable rank 2.057.
```

> **DECORRELATION IS NOT THE LEVER. SKILL IS.** A source that is perfectly orthogonal to
> everything we own and merely as good as our best buys **0.02 A**.

**AND THE SCOPE LIMIT, WHICH MATTERS MORE THAN THE RESULT -- I do NOT want this used against lane
P.** The bound above is the S29 section 5.1 displacement law, and that law governs operators in
**POOL SPACE**, where everything is capped at a cosine. **Lane P is working on the PRIOR, where
the record's derivative is -2.15 A per unit** (`prior-derivative-is-the-only-steep-lever`) and
where the ORACLE prize is 0.525 A on one named mode. **Section 4 prices "add another field to the
21"; it does not price "fix the distance prior's shape". Those are different levers and the flat
one must not be quoted at the steep one.**

### 5. DERIVATION: LANE F's ORACLE 0.169 IS A CANCELLATION ARTEFACT, NOT AN INFORMATION GAIN

Same closed form, with the second field carrying **zero** skill (`r2 = 0`):
`rho_max = r1 / sqrt(1 - c^2)`. A skill-free partner inflates the ORACLE combination purely by
cancelling the part of `v1` orthogonal to `u`, fitted with the native.

    lane F measured        0.169 ORACLE against 0.1128 best-single  = inflation 1.50x
    a ZERO-SKILL partner at c = 0.745 produces exactly              = inflation 1.50x
    and leave-fold-out the same combination gives 0.012, BELOW the best single field

**So lane F's two numbers are jointly consistent with the 21 fields containing no combinable
information at all**, and the ORACLE figure should be reported with that reading attached rather
than as "0.169 is reachable in principle". This strengthens lane F's own conclusion rather than
weakening it.

### 6. FAMILY REJECTIONS FOR THE COORDINATOR'S ITEM 3 (bias unidentifiable from corrupted data)

| family | the assumption it needs | why it fails here |
|---|---|---|
| **errors-in-variables**, two error-laden measurements identify the reliability ratio from their covariance | the two measurements' errors are **INDEPENDENT** | our two are the distogram profile and the pool profile, and **S19 sec 4 measured that they share the bias**. The assumption is not merely unverified, it is measured false |
| **Reiersol (1950)**, EIV identified without an instrument when the latent is **non-normal** | a linear relation between two observed error-laden variables, plus higher moments | we do not have two observed error-laden variables in that relation -- we have one prediction and no second observable of the same quantity. **Does not transfer**; recorded because it is the standard answer to "identify without an instrument" and someone will propose it |
| **multichannel blind deconvolution** (Xu et al. 1995 and successors) | the channels are **COPRIME** -- no common zeros | a shared bias IS a common factor; coprimeness is exactly what S19 sec 4 and S24's 0.9432 measure to be absent |
| **single-channel blind deconvolution** | sparsity / non-negativity / known support | none of these hold for a distance-profile bias |
| **instrument calibration** | a reference standard, i.e. an anchor | the anchor is the native (this is S30-L7's E1/E3 again, and S29-L4's control-variate rejection) |

**The pattern across all five is one condition under five names**: every method that separates a
shared bias from the truth needs **two views whose errors are independent**, and this instrument's
two views are measured to share the bias. That is the same sentence as S30-L7's E3, arrived at
from the measurement-error literature instead of from the likelihood.

### 7. WHERE I COULD BE WRONG, AND ONE THING I UPDATED AGAINST MYSELF

- Section 4's closed form assumes the optimal 2-field combination is taken **with the native in
  hand**; leave-fold-out is strictly worse, so the 5.1% discount is an **upper bound** on what
  orthogonality is worth, not an estimate. That direction favours my conclusion, which is the
  direction I should be most suspicious of, and I state it for that reason.
- Section 5's "exactly what a zero-skill partner at c = 0.745 produces" is a **consistency**
  argument, not a proof: a genuinely skilful partner at higher `c` produces the same inflation.
  It shows the ORACLE number does not *require* information, not that no information is present.
  The leave-fold-out 0.012 is what makes the reading likely, and that is lane F's number, not mine.
- I have **not** verified that S19's separation-profile result survives S29's benchmark pinning
  (memory `benchmark-and-folds-must-be-pinned` records that correcting the identity clustering
  silently moved 13 targets). S19 predates that correction. **Lane P should re-check the 0.525 A
  on the pinned benchmark before building on it** -- I flag this rather than assume it, and it is
  the single most load-bearing unverified number in this entry.
- **Against myself:** S30-L13 predicted at 2:1 that the AMBER-relax benefit concentrates on
  high-dispersion targets. Lane F then measured that pool Rg dispersion does not predict filter
  failure (-0.128, CI includes zero). Different outcome variable, so not a refutation, but
  same family and pointing the other way. **I have lowered that prediction to roughly even** and
  would rather it be run by the artefact's owner than defended by me.

### 8. ADDED AFTER LANE D's TEST LANDED: MY PREDICTION IS CONFIRMED AND **MY OWN CONSTRUCTIVE PROPOSAL IS DOWNGRADED BY IT**

Lane D ran the four-line check I offered and it came back confirmed and harder than I asked:
radial share of the Gram trace **0.5798**, cos(dominant principal direction, radial) **0.9472**
(median 0.9840), stable rank **1.705 -> 2.642** with the radial direction removed, and the ORACLE
payoff **cos(direction to the native, radial) = -0.0675 overall, -0.2524 on FAIL18**.

**The causal chain is now complete across three lanes, and it is this sprint's strongest result:**

    S30-L6 (L, algebra)   every fixed-reference distance potential contains a SEPARABLE
                          pure-scale term -- forced, not fitted
        ->  S30-L7 (D)    58.0% of the field library's Gram trace IS the radial direction,
                          and 94.8% of its top eigenvalue is
        ->  S30-L7 (D)    the direction to the native is ORTHOGONAL to radial (-0.068),
                          and ANTI-ALIGNED on FAIL18 (-0.2524)
        ->  lane F        shape, not scale, is 83% of the tail's error

> **The library is built out of scorers whose algebra forces them to measure scale; the native is
> not in the scale direction; and on the hard targets the scale direction points the wrong way.
> The library cannot express shape because its rank is spent on scale, and the reason its rank is
> spent on scale is the reference state.**

**AND NOW THE PART THAT COSTS ME.** In S30-L6 I offered a size-matched reference state
(`a_i = sqrt(5/3) Rg_i`) as a *constructive* alternative to statistical partialling. Lane D's
entry contains the warning that kills the deployable half of that offer: *"Deflating the radial
direction reallocates rank; it does not create any."* I initially wanted to argue that
size-matching is not deflation because it **builds a different field** rather than projecting an
existing one. **That argument is wrong, and my own algebra is what refutes it:** the score
decomposes *additively and exactly* into a shape term plus a scale term, the field is the gradient
of the score, so the field decomposes additively too, and **the size-matched field IS the deflated
field.** Lane D's caveat therefore applies to my proposal in full.

**What survives, precisely:** the two uses that never pass through the cosine --
(i) a **band statistic** with provably zero uniform-scale loading, which is S29-L19's
self-fulfilling-null warning answered by construction, and (ii) **lane R's ladder control**. Both
stand. **What does not survive is any suggestion that a size-matched channel is a route to
Angstroms**; it is the same 42% of the spectrum lane D already measured, relabelled. S30-L6's own
section 5 hedged this correctly ("worth approximately nothing as a deployable ranker") and I am
promoting that hedge to the headline where it belongs.

**Where this leaves the direction, honestly.** A new field that is orthogonal to the radial
direction is *necessary* (58% of the current rank is wasted) and *not sufficient* (section 4:
orthogonality is a 5.1% discount; the field still needs rho ~ 0.34 of its own). Lane D's Gram tool
is now a one-liner for checking a proposed field's radial orthogonality **before** it is built,
and I would make that check a gate on any new channel this sprint funds.

## S30-L19 -- L11, IS NATIVENESS RECOGNISABLE FROM SINGLE-STRUCTURE GEOMETRY? **F-R1 DOES NOT FIRE AT n = 126**: ON A LADDER WHERE KIND, LOCAL REALISM AND PERTURBATION BUDGET ARE ALL MATCHED, ORDERING SURVIVES (DIS +0.347, +0.134 ABOVE ITS ANCHOR CONTROL, p_max 0.000) BUT **PREFERENCE DOES NOT ON ANY OF 43 CHANNELS** -- AND THE LEAVE-FOLD-OUT COMBINATION THAT PREFERS A 0.55 Å STRUCTURE TO PRODUCTION ON 93.0% OF TARGETS PREFERS A **RANDOM POOL MEMBER ON 100%** AND A **3 Å RUNG ON 100%**, SO ITS MARGIN IS **-0.070 [-0.110, -0.028]**; LEG_torsion, S29 §12.0's LAST OUTSIDE-CLASS-M HOPE, CLEARS ITS ANCHOR CONTROL BY ONLY +0.024 (A QUARTER OF THE REGISTERED MARGIN) AND SITS AT 0.503 PREFERENCE AGAINST A 0.698 POOL-MEMBER CONTROL; AND D1 GIVES THE MECHANISM -- THE RMSD SIGNAL IS **ABSENT FROM LOCAL FEATURES (ΔR² -0.089) AND ABUNDANT IN GLOBAL ONES (+0.600)**, SO THE NULL FOR EVERY PER-RESIDUE CHANNEL IS A THEOREM ON THIS INSTRUMENT, NOT AN EMPIRICAL MISS (2026-09-20 13:27, R)

Pre-registered in `s30/PREREG_S30_R.md`, committed `7eabffee` at 12:40:35, BEFORE any aggregate
over more than the single probe target the file declares. Registered prior: **F-R1 does not fire,
about 4 to 1** -- so this null confirms my own expectation and must be read with that discount.
Code `s30/s30_R_ladder.py` + `s30/s30_R_agg.py`; tests `tests/test_s30_R.py` (5 pass); artefacts
`s30/results/s30_R_verdict.json`, `s30_R_rows.s{0,1,2}of3.jsonl` (126 targets),
`s30_R_stability.json`, `s30_R_chains/` (126), `s30_R_repro/`.

### THE INSTRUMENT, AND ITS FLOOR -- WHICH TRAVELS WITH EVERY SENTENCE BELOW

Every rung is an ideal-geometry backbone built from (phi, psi) by `core.geometry
.build_backbone_batch`: identical bond lengths and bond angles at every rung, no projection, no
coordinate average, no contraction. Perturbed residues draw torsions from the fold's
**leakage-safe** Ramachandran table `s8/generate_rama.npz[fold]`. Three families per target:
**ladder A** anchored on the native's own torsions, **ladder B** anchored on a random real pool
member, m ∈ {1,2,3,4,6,8,12} residues resampled × 32 draws each; plus the rebuilt native, the
projected PROD / circ_best / sub0 chains, and the 500 real pool members.

```
torsion-rebuild FLOOR        0.347 A   <- the ladder's zero, NOT 0 A
near-native band             0.555 A mean over 18.8 rungs/target; 3 targets fall back
PROD chain 3.207 A   circ_best 0.252 A (ORACLE)   sub0 0.498 A (ORACLE)   pool best 1.711 A
ladder share of the set      0.465 -- reported, not hidden; the four SET-REFERENCED channels
                             (CONS, DMAP_CONS, TORS_CONS, POOLGO) are referenced to the POOL ONLY
```

**Why the ladder is the point.** At m = 1 -- ONE resampled residue -- CA-RMSD to the native spans
0.284 to 4.685 A on the probe target. At a fixed perturbation budget, with local torsions drawn
from the same distribution, two structures differ by one torsion and by 4.4 A of nativeness.
Ordering that is the experiment; S28-L48 could not ask it because its rungs differed in KIND
(S30-L1).

**A1, reproduction.** My independently projected PROD chains give **3.2071 A against the record's
3.2126** (`s27/results/chain_rows.jsonl :: DIS`), median per-target |Δ| **0.0012 A**, with one
0.513 A outlier -- which is the known multi-start branch flip already documented in
S28-L18/L27b/L43. A 17-sprint-old number reproduces and the one deviation has a named cause.

**A2, realism flatness (ORACLE; this audit can only WEAKEN my own positives).** The ladder is not
perfectly realism-flat: within m, RAMA **+0.116** [+0.097, +0.137], EXVOL **+0.161**
[+0.099, +0.220], |Rg - median pool Rg| **+0.190** [+0.157, +0.233]. Structures that land farther
from the native are slightly less Rama-typical, slightly more clashing, slightly less typical in
size. **That residual gradient would help a positive, so it makes every null below stronger and
every positive below smaller than it looks.** DIS's +0.347 is about twice the largest gradient.

**A3, no native.** The object every channel is handed is a bare namespace with **no `nat_ca` and
no `oracle_rr` attributes at all**, and `ham_lib.Context` consumes the universe through exactly
four keys (W, S, PHI, PSI). A channel cannot read a native even by accident -- a construction
argument, which is stronger than a poison run, and `tests/test_s30_R.py::test_t4` locks it.
Every RMSD column is computed only AFTER every channel value exists.

**A5, ties.** No `np.argmin` on a tied signal anywhere; all preferences score ties at 0.5 and
average over the tied set.

### THE VERDICT: F-R1 DOES NOT FIRE

```
                     rhoA_part rhoANCH_p  A-ANCH   rho_rg  prefN prefPool pref-ctl  pctNat  F-R1
DIS                     +0.347    +0.213  +0.134   -0.084  0.058    0.020   +0.038   0.291   i
CONS                    +0.343    +0.399  -0.056   -0.374  0.057    0.199   -0.142   0.434   -
DIS_MEAN                +0.335    +0.234  +0.101   -0.176  0.062    0.025   +0.037   0.326   i
CONS_SI                 +0.315    +0.366  -0.051   -0.807  0.079    0.238   -0.159   0.423   -
CONTACT_LL              +0.294    +0.229  +0.065   -0.601  0.213    0.134   +0.080   0.297   -
DMAP_CONS               +0.263    +0.316  -0.053   -0.050  0.138    0.297   -0.160   0.421   -
LEG                     +0.210    +0.198  +0.013   -0.124  0.269    0.427   -0.158   0.283   -
DISTPOT_SI              +0.193    +0.081  +0.112   -0.022  0.367    0.254   +0.112   0.315   -
RAMA                    +0.107    +0.090  +0.017   +0.051  0.640    0.790   -0.150   0.639   -
CONTACT                 +0.084    +0.009  +0.075   +0.031  0.558    0.443   +0.115   0.371   -
LEG_torsion             +0.076    +0.053  +0.024   -0.025  0.503    0.698   -0.195   0.503   -
LEG_compactness         +0.033    +0.020  +0.013   +0.966  0.435    0.433   +0.001   0.424   -
RG_LAW_SI               +0.000    +0.007  -0.006   -0.010  0.485    0.496   -0.011   0.326   -
(43 channels in the artefact; nothing selected; full table in s30_R_verdict.json)
```

- **Clause (i), ordering, FIRES on 2 of 43** -- DIS (+0.347, anchor contrast +0.134) and DIS_MEAN
  (+0.335, +0.101). Multiplicity is not the explanation: the max-over-channels per-target
  sign-flip null (500 draws) gives null mean 0.062, p95 0.105, **p_max 0.000**.
- **Clause (ii), preference, FIRES ON NOTHING.** The best `pref_near` in the entire library is
  **RAMA at 0.640**, under the 0.65 bar -- and RAMA prefers a **random real pool member** to
  production on **0.790** of targets, so its contrast is **-0.150**. The S28-L36 pool-member veto
  fires on every channel that clears the bar-adjacent range.
- **The largest preference effect in the whole library points the wrong way**: DIS prefers
  production to a **0.55 A** near-native structure on **94.2%** of targets (`pref_near` 0.058).
  That is S28-L48 reproduced with the kind confound REMOVED -- and it is sharper, not weaker.

### THE RESULT THAT MATTERS MOST: THE COMBINATION LOOKS LIKE RECOGNITION AND IS NOT

Ridge over 40 channels, fit on four folds, evaluated on the fifth (prereg item 6):

```
near-native (<=1 A) vs PROD   real       pref 0.930   pool-member control 1.000   margin -0.070
                                         fold CI [-0.110, -0.028]   1.04x MDE   n = 114
GARBAGE CHECK far (>=3 A) vs PROD  real  pref 1.000   control 1.000   (both SATURATED)
near-native, label-shuffled              pref 0.439   control 0.325              margin +0.114
```

**Read it in this order.** The combination prefers a 0.55 A structure to production on 93.0% of
held-out targets. Quoted alone that is a headline. It prefers an **arbitrary real pool member** on
**100%**, and a **3 A rung** on **100%** -- i.e. it prefers structures that are FARTHER from the
native MORE often than the near-native one. (Those two saturate at 1.000, so the garbage-check
margin is not a measured contrast and is not quoted as one; the *near*-vs-control margin is the
measured number, **-0.070 [-0.110, -0.028], 1.04x MDE, 4/5 folds**.) It has learned *"is this the projected production
average?"*, not *"is this near-native"*. The margin excludes zero in the wrong direction.
This is `decoy-bank-not-a-pool-proxy` and the S28-L36 veto, reproduced on the cleanest instrument
the project owns, and it is the reason clause (ii) had a control clause at all.

**DISCLOSURE, my own bug, caught before it was reported.** The first run of this cell printed
`0.500 / 0.500 / 0.500` on all three arms. That was a **NULL-INPUT ARTEFACT**: the channel filter
required `d_near` on *every* row, and 3 targets have no rung under 1 A, so the channel set was
EMPTY and the "preference" was `mean(0 < 0) + 0.5*mean(0 == 0)`. Three identical 0.500s across
three different questions is what caught it. Fixed by dropping the offending ROWS rather than the
channels; the disclosure is in `s30_R_verdict.json :: LFO_note`.

### D1 -- THE MECHANISM, AND WHY HALF THE NULL IS A THEOREM RATHER THAN A MISS

Held-out (5-fold within target) R² of RMSD-to-native on ladder A, with the m one-hot as the
baseline block, feature counts deliberately MATCHED (84 local vs 80 global):

```
r2_m_only    +0.3111  [+0.2917, +0.3273]     knowing only HOW MANY residues moved
r2_local     +0.2216                ->  d_local   -0.0894  [-0.1217, -0.0522]
r2_global    +0.9108  [+0.9013, +0.9201]  ->  d_global  +0.5997  [+0.5850, +0.6166]
```

The LOCAL block is **ORACLE-ADVANTAGED**: besides per-residue sin/cos of (phi, psi) it is handed
the **per-residue circular deviation from the native anchor itself**. With that advantage, and at
matched capacity, **it adds NOTHING beyond knowing m** -- the held-out ΔR² is negative, meaning it
only costs capacity. The global block (the CA distance map, Rg, contact count) adds +0.600.

**Therefore the null for every per-residue channel -- RAMA, LEG_torsion, DSSPHB, CAGEO,
LEG_hbond_local, HP -- is a property of the geometry and not an empirical miss.** One torsion at
mid-chain swings global RMSD by several Angstroms through a lever arm; a sum of per-residue terms
cannot see a lever arm. The live question narrows exactly as the coordinator predicted it would:
only a globally-reaching function could order this ladder. **And the only globally-reaching
functions the project owns are distogram re-readings (DIS, DIS_MEAN, CONTACT_LL -- inside class M,
bounded by lane T's theorem 2) and pool-consensus terms (CONS, DMAP_CONS, POOLGO), whose anchor
contrasts are zero or negative (-0.056, -0.053, -0.050): they measure TYPICALITY, not nativeness.**
The ORACLE label on r2_global matters and is stated: it is a per-target fit whose labels encode the
native's own distance map, so it bounds what a target-specific native-informed global function
could extract -- it does NOT say a native-free one exists.

### D2 -- THE RESOLUTION: WHERE THE ORDERING DIES (chance = 0.500 exactly, by construction)

Pairwise concordance within m, by |Δ RMSD|:

```
                  0-0.25  0.25-0.5   0.5-1     1-2     2-4      >4
ALL PAIRS  DIS     0.525     0.566   0.618   0.698   0.771   0.822
           LEG     0.511     0.525   0.545   0.584   0.656   0.749
           LEG_torsion 0.499 0.511   0.514   0.528   0.561   0.607
BOTH <= 2 A (the only regime the endpoint cares about)
           DIS     0.520     0.554   0.570   0.624     n/a     n/a
           LEG     0.526     0.555   0.574   0.618     n/a     n/a
           CONS    0.500     0.486   0.528   0.572     n/a     n/a
```

**The signal is coarse triage and nothing else.** The best channel separates 4 A from 0.3 A at
0.822 and separates two near-native structures 0.25 A apart at **0.520** -- two points above a coin
toss. For lane Q: a sparse-support rule built from this library can only ever select on
differences of ~2 A or more. For lane T: an objective that had to resolve 0.5 A would need a
channel roughly four times as discriminating as the best one in the record. Note also that inside
the near-native band **LEG is as good as DIS** (0.526/0.555/0.574/0.618 against
0.520/0.554/0.570/0.624) -- the distogram's advantage is entirely in the coarse regime.

### THE SPECIFIC NEGATIVE S29 §12.0 ASKED FOR: LEG_torsion

S29-L50 made LEG_torsion the last live exit -- in-band skill +0.181 that is demonstrably not
compactness, and a function of the structure rather than of the distogram, hence outside class M.
On this instrument, with kind and budget matched:

```
rho_A partialled +0.076   anchor contrast +0.024   rho(score, Rg) -0.025
pref_near 0.503  vs pool-member control 0.698  ->  contrast -0.195
the rebuilt NATIVE sits at the 0.503 percentile of its OWN Rama-resampled ladder
resolution in the near-native band: 0.507 / 0.525 / 0.567 / 0.574
```

**Correction to my own wording, made on recomputing from the artefact rather than from this
entry's prose (the S29-L48 standing rule):** "at chance" is too strong and I withdraw it. The
anchor contrast is **+0.024 with fold CI [+0.010, +0.038], which EXCLUDES zero** -- LEG_torsion
does order nativeness slightly above its own anchor control. What is true is the weaker and still
decisive statement: that contrast is **a quarter of the registered +0.10 margin**, the channel
prefers a near-native structure to production on 0.503 of targets against a pool-member control of
0.698, and it ranks the native at the exact median of the native's own perturbations. S29's caveat -- that in-band skill is not deployable value -- turns out to
understate it: here the in-band skill does not even reappear as ordering once kind and budget are
matched. **Row 3 of S29 section 7 should now be closed, negatively, with a mechanism (D1: it is a
per-residue sum, and the signal is not in per-residue features).**

### WHERE THE NATIVE SITS IN ITS OWN LADDER (the sentence that carries the whole result)

`pctile_nat_in_A` = the share of ladder-A rungs -- random Rama-resampled perturbations of the
native -- that score BETTER than the rebuilt native itself:

```
DIS 0.291   LEG 0.283   CONTACT_LL 0.297   DSSPHB 0.294   DISTPOT_SI 0.315   CAGEO 0.340
LEG_torsion 0.503   LEG_contact 0.499   ELEC 0.463   RG_LAW 0.463   RAMA 0.639
```

Under the shipped cost, **29% of the native's own random perturbations score better than the
native**. Under the Ramachandran term -- the one channel whose entire job is local plausibility --
**64% do**, because a draw from the fold's Rama table is typically MORE Rama-typical than a real
native torsion. Nothing in the library puts the native first.

### LANE L's SIZE-MATCHED FIX: IT WORKS WHERE PREDICTED, AND IT DOES NOT MEAN WHAT IT SOUNDS LIKE

Implemented stronger than requested: instead of only replacing the reference term, the universe
AND the candidates are both rescaled to a common Rg, so each pair potential is **fitted and
evaluated in reduced units** and a uniform contraction is exactly invisible. The audit is not an
assertion -- the two PURE functions of Rg become constant:

```
RG_LAW   relative sd  raw 0.5009  ->  SI 3.2e-16        RG_UNIV  0.5071  ->  2.8e-16
```

- **On the channel lane L's derivation named, it helps**: DISTPOT_SI beats DISTPOT on ordering
  (+0.193 vs +0.136), on the anchor contrast (+0.112 vs +0.056) and on the preference contrast
  (+0.112 vs +0.078), while its Rg loading falls to -0.022. That is the predicted direction.
- **And one thing lane L should have back, because it is not obvious**: scale invariance by
  construction does **NOT** imply ρ(score, Rg) = 0, and on this instrument several twins are MORE
  Rg-correlated than their parents (CONS -0.374 -> CONS_SI **-0.807**, DMAP_CONS -0.050 ->
  **-0.729**, POOLGO +0.611 -> **-0.678**). The reason is real rather than a bug: the transform
  removes the *forced* scale term, but at peptide length shape and size are genuinely correlated,
  and what is left is that correlation. **Removing the forced term is a correctness fix for the
  in-band question; it is not a decorrelation, and it changes no endpoint** -- by the S29 section 8
  bound it is a re-parameterisation and reaches the endpoint only as a cosine. ANDIS's
  recognition/discrimination trade is the right frame and this ladder is on the side where
  removing the size term is correct, which is why it is reported here and not proposed as an arm.

### D3 -- THE ANCHOR CONTROL'S OWN CONFOUND, MEASURED BEFORE THE VERDICT WAS READ

`corr(distance-to-anchor, distance-to-native)` within m, over 126 targets: **median +0.259**, mean
+0.214, [p10 -0.545, p90 +0.847], 37% of targets above 0.5. The prereg's caveat triggers only if
the median exceeds 0.5, so **it does not trigger**: clause (i)'s +0.10 margin was reachable, and
two channels did reach it. The confound-immune column (`rho_B` with rank(distance-to-anchor)
partialled out) is in the artefact and tells the same story.

### STABILITY OF THE INSTRUMENT, FROM AN ACCIDENT

A seeding defect -- `hash()` on a str is process-randomised -- meant a duplicate shard drew a
DIFFERENT ladder for 38 targets. Rather than discard them I measured with them
(`s30_R_stability.json`), and the answer is worth having:

```
rho_A_part       draw-to-draw corr 0.912   mean |diff| 0.079     means +0.099 / +0.087
pref_near        draw-to-draw corr 0.982   mean |diff| 0.027
pref_pool        draw-to-draw corr 1.000   (does not depend on the ladder)
rho_ANCHOR_part  draw-to-draw corr 0.147   mean |diff| 0.235     means +0.107 / +0.081
d_local  -0.127 / -0.134        d_global  +0.619 / +0.596
```

The ordering and preference columns are stable; **the anchor control is NOISY per target**
(corr 0.147) because ladder B redraws its anchor, so the per-target A-ANCHOR contrast carries that
noise even though its mean is stable. Any future use of that contrast should average several
anchors. The defect is fixed (`zlib.crc32`) for every later run; the numbers in this entry come
from the pre-fix seed and the duplicate draws are exactly what prices that.

### THE HONEST ANSWER TO L11

**No -- not in the sense the sprint needs.** Precisely:

1. Structures are **not** indistinguishable: a within-kind, budget-matched ordering signal exists
   and is not multiplicity (DIS +0.347, p_max 0.000, +0.134 above its own anchor control).
2. But that signal is **the distogram's**, i.e. the objective production already optimises -- which
   is why "DIS prefers production" is close to tautological -- and it is **coarse**: 0.520
   concordance at 0.25 A separation in the near-native band.
3. **No channel and no leave-fold-out combination prefers a 0.55 A structure to the 3.2 A
   production average at the registered bar**, and the combination that appears to does so less
   often than it prefers a random pool member or a 3 A rung.
4. Half the null is **structural**: the signal is not in per-residue features at all (ΔR² -0.089
   against +0.600 global), so every local channel is blind by construction.

**What this forecloses.** (a) S29 §12.0 row 3 -- LEG_torsion as the outside-class-M route -- closed
negatively with a mechanism. (b) A native-free support rule for §12.1's sparse weighted readout
**built from this library**: choosing 2 of 500 is ~17.9 bits, and a library that resolves 0.25 A at
0.520 cannot pay it. (c) Any "better Hamiltonian assembled from the existing channels": the only
ordering present is the one production already uses.

**What it does NOT foreclose**, and I want this stated as carefully as the negative: D1 says the
information IS in the global shape (R² 0.911, ORACLE per-target fit). What is missing is a
*native-free* globally-reaching channel. That is a supply problem, not an impossibility -- and it
is a different problem from the one the sprint has been attacking. Nothing here touches the
retrieval/prior lever, which `prior-derivative-is-the-only-steep-lever` still prices as the only
steep one.

### MULTIPLICITY, MDE AND SCOPE

0 endpoint (deployable RMSD) comparisons. 43 channels × 3 ordering columns reported as ONE table,
none selected; one max-over-channels per-target sign-flip null per column (500 draws): rho_A_part
observed 0.347 vs null mean 0.062 / p95 0.105, p_max 0.000; pref_near observed |0.443| vs null
0.094 / 0.132, p_max 0.000. Every cell carries a fold-clustered CI and a per-comparison MDE in
`s30_R_verdict.json` (`mde-is-per-comparison-not-per-instrument`). ρ(score, Rg) travels beside
every recognition number, as lane L required. Scope: 9-16mers, 126 dev targets, CA-RMSD, a ladder
whose floor is 0.347 A and whose realism is flat only to ±0.12-0.19 Spearman.

**Operational disclosure.** `s26/jobrun.py`'s `CPU_START = 85.0` gate was closed throughout (lane D
held the box at 99.8%), so the three run shards went out **detached rather than through the
governor**. RAM -- the binding constraint per `machine-fits-two-heavy-jobs` -- had 4.0 GB free and
peak per process was ~0.35 GB. The inherited gate defect named in `s30/STATUS.md` is real and it
bit exactly as predicted.

## S30-L20 -- **NO.** A GENERATOR CANNOT HAVE A BETTER TYPICAL MEMBER IN ANY SENSE THAT PAYS, BECAUSE "SET MEAN" IS **TWO** QUANTITIES AND ONLY ONE CONVERTS: `set_mean^2 = B^2 + S^2` WITH `B` = THE ENDPOINT ITSELF. THE SPREAD HALF IS WORTH **ZERO** (THE TERMINAL ALREADY EXTRACTS IT AT r = 0.885; IT IS ORTHOGONAL TO BIAS AT r = +0.085) AND THE BIAS HALF **IS** THE ENDPOINT, WHICH LANE L PROVED NON-IDENTIFIABLE AT ANY K. **NOTHING IN THE RECORD IMPROVES THE SET MEAN AT ALL** -- THE BEST IS A TIE AT 0.57x MDE -- AND THE MOST CONCENTRATED SOURCE EVER BUILT (`T0_helix`, S = 0.57 AGAINST THE POOL'S 1.58) IS **THE WORST ENDPOINT IN THE RECORD AT 3.789**. ALSO A CORRECTION TO MY OWN S30-L10 (2026-09-20 13:29, X)

Pre-registration `s30/PREREG_S30_X.md` **ADDENDUM 1 @ 539de5a8**, committed before any number
below existed, with my registered prediction that the answer would be NO. Findings
`s30/s30_X_FINDINGS.md` PART II. Code `s30/s30_X_typicalgood.py`. Results
`s30/results/s30_X_typicalgood.json` (504 cells + 126 pool rows). **No VQE or pipeline compute.**

### 0. The question, and why measuring it answers a different one

> *"Can a generator be built whose TYPICAL member is better than the pool's -- Δ(set mean) ≈ −0.2 Å
> -- rather than one whose BEST member is better? Measure the set mean first."*

I measured it first, as instructed, and **the set mean turns out to be two quantities.** For a
coordinate-average terminal over m members the common-mode identity gives, per target:

```
    set_mean^2  ~=  B^2 + S^2        B = RMSD(set average, native) = THE ENDPOINT
    endpoint     =  B                S = the set's spread about its own centroid
```

**A set mean improved purely by CONCENTRATION (S down, B fixed) moves the endpoint by zero, by
algebra.** Only B converts, and B *is* the endpoint. Measured, n = 126 per row:

```
  arm            set_mean   B (endpoint)      S    avg gain |  d_set_mean     d_B      d_S   winsBOTH
  POOL (top-75)    3.5507      3.0483     1.5770    0.5023  |      --         --       --        --
  T0_helix         3.8663      3.7892     0.5676    0.0771  |   +0.3156   +0.7408  -1.0094    22/126
  T1_blind         4.0218      3.2435     2.2242    0.7784  |   +0.4712   +0.1951  +0.6472    31/126
  T2_restype       3.9563      3.2065     2.1267    0.7498  |   +0.4056   +0.1581  +0.5497    36/126
  T3_pool          3.5916      3.1752     1.4791    0.4164  |   +0.0409   +0.1269  -0.0979    36/126
```

### 1. FOUR INDEPENDENT REASONS THE ANSWER IS NO

**(1) Nothing in the record improves the set mean at all.** Every `d_set_mean` is **positive**
(+0.3156 / +0.4712 / +0.4056 / +0.0409), the first three at 5/5 folds with fold CIs excluding zero.
The best, T3_pool: **+0.0409, SE 0.0254, MDE 0.0712, 0.57x MDE, fold CI [-0.0134, +0.0887], 4/5 --
NOT MEASURED**, a statistical **tie** with the pool. Against a −0.2 Å requirement, **the record's
best is a tie at zero**, and even that arm's endpoint is +0.1269 worse.

**(2) The concentration half is worth zero, as the algebra says.**

```
    avg_gain  =  0.4143 * S  -  0.1559        r = 0.8854,  n = 630
```

**The terminal's entire value is spread extraction**; a concentrated source hands it nothing to
extract. Within target, `corr(S, B) = +0.0851` (n = 504): **concentration is orthogonal to bias.**

> *Honesty note on that correlation.* The raw within-target `corr(S, B) = -0.4744` (n = 630) is
> **driven by T0_helix alone**; drop that one arm and it is +0.0851. My registered prediction was
> "near zero or negative" and both satisfy it, but **the strong negative is an artefact of the
> zero-information arm and must not be quoted.** The durable number is +0.085.

**(3) The extreme case proves it.** `T0_helix` -- a constant alpha-helix with 15° jitter -- is by a
wide margin the most concentrated source ever built here, **S = 0.5676 against the pool's 1.5770**.
It *is* the typical-good generator the question asks for. **It is the worst endpoint in the record:
3.7892 against 3.0483.** Its averaging gain collapses to 0.0771 against the pool's 0.5023 -- the
terminal had nothing left to extract, so the endpoint fell back onto B.

**(4) The cells that win on both do not transfer.** 125 of 504 cells (**24.80%**) beat the pool on
set mean *and* endpoint.

> **That essentially TIED its registered bar -- 24.80% against 25%, a miss of 0.20 percentage
> points. The counting half decided nothing and I will not pretend otherwise.**

The verdict rests on the transfer arm, which is unambiguous. Pick the arm winning on both most
often in half the targets, score it on the other half, 400 splits:

```
    split-half transfer of d_B  =  +0.1597    CI95 [+0.0676, +0.2687]
```

**Wrong sign, CI excluding zero** -- the selection makes the endpoint **0.16 Å worse**. The
wins-on-both cells are an order statistic (`grid-oracles-are-order-statistics`).
**VERDICT: H-X3 STANDS.**

### 2. WHY THIS IS A CEILING AND NOT A MISS

Endpoint = B = the shared bias. *"Improve the typical member"* decomposes into a **spread** half --
free, worth 0, already extracted -- and a **bias** half, which **is** the endpoint. **Lane L's
S30-L7 closes the second half structurally**: the likelihood depends on `(t, mu)` only through
`t + mu`, so mu is non-identifiable at any K, and provenance cosine **0.9432** vs a **0.9330**
within-source control says mu is a property of **the prior candidates are scored against, not of
where they come from.** The only half of "typical-good" that pays is the half no source change can
move.

**Lane L's falsifier is applied and ALREADY FAILED before any endpoint run**, on its own rule
(require cosine < 0.9330 first): the measured value is 0.9432 at n = 126
(`s24/results/qmatch.json`, S24 L3, replicated independently by S24 lane E). **No endpoint run was
spent.** Lane T's register arithmetic is accepted as binding and no torsion-vector encoding is
proposed.

### 3. THE PREMISE, INVERTED

The question came with a premise -- *"every sampler was built to be diverse, not typical-good;
nobody has built one to be typical-good."* **Somebody did, and diversity is the correct design.**
T0_helix is that generator and it is the worst thing in the record, while the arms with the
*largest* spread have the largest averaging gains (T1 S 2.224 -> 0.778; T2 S 2.127 -> 0.750) and
the best generated endpoints. **For an averaging terminal, spread is the raw material, not a
defect. "Typical-good" is the wrong design goal and the samplers were not built wrong.**

### 4. CORRECTION TO MY OWN S30-L10 (rule 15; the original wording stands there unedited)

S30-L10's admission condition reads `Δ(set mean) < -(0.32..0.37)*Δ(set best)`. **It treats
Δ(set mean) as one channel when it is two, and only the bias channel converts.** Amended:

> **ADMIT G against P iff `Δ(bias B) + 0.298*Δ(set best) < 0`.** Δ(set mean) is admissible as a
> proxy **only** to the extent it reflects ΔB; the spread component converts at ≈ 0.

S30-L10 §1-§2 stand unchanged -- the ceiling gate at 1/10, the frozen law passing at 0.0153 Å. What
is amended is the *interpretation* of its set-mean term, and the amendment makes the condition
**stricter, not looser**: the cheap half of it was free all along.

### 5. What this is NOT

`set_mean` is the arithmetic mean of member RMSDs, not the RMS, so S is Jensen-biased **small** --
it understates concentration and **cannot manufacture** the effect reported here. The primary
counting test tied its bar and is reported as undecided. Only H-X3's two pre-registered arms are
read as results; the per-arm decompositions are descriptive. No native was used to set any
parameter; every ceiling and every per-member RMSD row is ORACLE and says so in its field name.

> **NUMBERING COLLISION RESOLVED 2026-09-20 13:39 by the coordinator.** Posted as `S30-L9` at 2026-09-20 13:38 by lane D, which an entry stamped earlier already held. Renumbered to **S30-L24**; the earlier entry keeps `S30-L9`. Recorded rather than silently fixed (contract rule 15). Anything citing `S30-L9` for *this* result means `S30-L24`.

## S30-L24 -- THE METER'S NEW 8-DRAW CONTROL EARNS ITS KEEP ON ITS FIRST RUN: **S29's SEED-0 RANDOM-SIGNED DRAW WAS THE MAXIMUM OF ITS OWN EIGHT ON THE BUILT CHAIN.** THE C2-CLAUSE-2 CONTRAST FALLS **+0.0635 -> +0.0357**, FROM 0.92× TO **0.56× MDE**, FOLD CI NOW SPANNING ZERO -- BELOW THE SPRINT'S "NOT A RESULT" FLOOR. **ON CA THE SAME SINGLE DRAW WAS FINE** (+0.1746 -> +0.1716, STILL 1.84× MDE, 5/5 FOLDS). THE DRAW NOISE IS A PROPERTY OF THE BASIS, NOT OF THE CONTROL (2026-09-20 13:38, D)

**Verdict: one S29 number withdrawn, one confirmed. The built-chain preference contrast is not a
result; the CA one survives a control eight times stronger than the one it was published with.**

This is the first run of extension **(E4)** (S30-L3), which replaces S29's single seed-0
matched random-signed structure with R = 8 independent draws per target, and it is also the
first time any cost has been metered on the built chain through recomputed projections rather
than C2's stored rows.

### Built chain — the reporting basis

| | value |
|---|---|
| pref(circ_best vs PROD) | 0.0714 |
| pref(matched random-signed vs PROD), 8 draws | 0.0357 (per draw 0.008, 0.048, 0.008, 0.056, 0.056, 0.016, 0.048, 0.048; **sd 0.0212**) |
| **contrast, 8 draws** | **+0.0357, SE 0.0226, +0.56× MDE, fold CI [−0.011, +0.082], 3/5 folds** |
| S29's published value (seed 0 alone) | **+0.0635, +0.92× MDE** |
| what a single draw could have given | **[+0.0159, +0.0635]** |

**S29's seed-0 draw sits at the top of its own range.** The reported +0.0635 is the maximum of
the eight; the 8-draw mean is exactly **half** of it, and the effect moves from "0.7–1.0× MDE,
not a demonstrated improvement" to "**below 0.7× MDE, not a result**," with a fold CI that now
includes zero and folds agreeing 3/5 instead of 4/5. Nothing was done wrong in S29 — a single
draw is an unbiased estimate of the control — but its variance was never measured, and here it
is 59% of the mean.

### CA point cloud — the same control, the opposite conclusion

| | value |
|---|---|
| pref(matched random-signed), 8 draws | 0.0347 (per draw 0.016–0.048, **sd 0.0094**) |
| **contrast, 8 draws** | **+0.1716, SE 0.0333, +1.84× MDE, fold CI [+0.086, +0.255], 5/5 folds** |
| S29's published value (seed 0 alone) | +0.1746, +1.84× MDE |
| what a single draw could have given | [+0.1587, +0.1905] |

Here the single draw was representative and **S29's number is confirmed to three decimals**. The
contrast is not concentrated: against the uniform-effect null the drop-top-10 mean sits at the
51st percentile [p10 +0.151, p90 +0.242], no flag. (Its *median* paired difference is 0.0000,
but that is the expected shape for a 0/0.5/1 indicator whose modal value is 0 — the
median-vs-mean warning does not transfer to bounded indicators, and the concentration test
against the uniform-effect null is the right instrument, which is why E4 prints it.)

### Why the two bases differ, and the rule it earns

The draw sd is **0.0212 on the chain against 0.0094 on CA**, while the chain's effect is a
quarter of the size. Relative draw noise: **59% on the chain, 27% on CA.** On the chain the
shipped cost prefers production to almost everything (every pref is 0.008–0.056), so the control
is a rare event and a single draw of it is correspondingly noisy. **A single-draw control is
least trustworthy exactly where the effects are smallest — which is where this project's
remaining effects live.**

> **Rule: a control that is itself a random draw needs its own draw distribution reported, and
> the number of draws must scale with how small the effect is.** Eight was enough to separate
> these two cases; one was not.

### What this does and does not change

- **Withdrawn as evidence:** the built-chain C2-clause-2 contrast, "f prefers the 0.29 Å ORACLE
  structure to production more often than it prefers a matched random displacement." On the
  reporting basis that is 0.56× MDE with a CI spanning zero.
- **Confirmed:** the CA-basis version, +0.1716 at 1.84× MDE with 5/5 folds.
- **Unchanged:** the meter's verdict for the shipped cost is **BLOCK on both bases**, driven by
  the S28 ladder (rho −0.4023 chain, −0.1818 CA, both fold CIs entirely below zero), not by this
  contrast.
- The CHARTER ladder rho is stable across draws (chain sd 0.016, CA sd 0.008), so the ladder
  numbers were never at risk from the single-draw choice — only the preference contrast was.

Artefacts: `s30/results/s30_D_meter_DIS_chain.json`, `s30/results/s30_D_meter_DIS_ca.json`;
ladder cache with chain projections and 8 matched draws at
`s30/results/s30_D_ladder_structs/` (126/126, built 13:36). Multiplicity: 30 comparisons from
the chain run, 24 from the CA run.

## S30-L25 -- **THE WHOLE NATIVE-FREE FEATURE SPACE EXPLAINS 0.8% OF THE ORACLE ERROR OUT OF FOLD** (BAR 1.96%, EXCESS OVER A MATCHED-DIMENSION CONTROL +0.0095 AT 0.75x MDE), WHILE THE SAME BASIS **LOCATES 82.7% OF THAT ERROR NATIVE-FREE AGAINST A CONTROL'S 44.4% AT 10.17x MDE** -- THE SUBSPACE IS FREE AND THE SIGN IS THE WHOLE PROBLEM. AND ON THE PRIOR ROUTE, WHICH THE BOUND DOES NOT PRICE: **S19's 0.525 A SEPARATION-PROFILE PRIZE REPRODUCES ON THE PINNED BENCHMARK AT -0.5740 A** (2.52x MDE, 5/5 FOLDS, 70.8% OF A PERFECT DISTOGRAM), **68% OF IT IS IN |i-j| >= 7**, IT CONCENTRATES **4x ON BOTH FILTER-INDEPENDENT TAILS**, AND IT REDUCES TO **FIVE SIGNS PER TARGET**: ORACLE SIGN WITH A LEAVE-FOLD-OUT MAGNITUDE BUYS **-0.3567 (62%)**, ORACLE MAGNITUDE WITH A LEAVE-FOLD-OUT SIGN BUYS **-0.0733 (NOT MEASURED)** (2026-09-20 13:43, P)

Pre-registration `s30/PREREG_S30_P.md`, committed **097352c5 at 13:25:18, before the first arm
existed**. Code `s30/s30_P_r2.py`, `s30/s30_P_prior.py`, `s30/s30_P_chain.py`. Results
`s30/results/s30_P.json`, `s30/results/s30_P_prior.json`, `s30/results/s30_P_chain.json`.

My brief was L3 -- find an information source whose errors are decorrelated from the distogram's,
then correct with it. The coordinator reformulated it mid-task into lane T's restatement of B2,
`rho_max = sqrt(R^2(e ~ S))`, which turns the enumeration into one out-of-fold regression. Lane L
then narrowed it again (S30-L18) by localising the target to a five-number separation profile and
by pricing orthogonality out of the picture (**decorrelation is a 5.1% discount; skill is the
lever**). Both reformulations were better than my brief and this entry follows them.

---

### 1. THE REGISTERED MEASUREMENT: EVERYTHING, THROWN AT `e`, OUT OF FOLD

`e` = `remove_rigid(oracle_direction(C0, nat), C0)`, C0 = production's `avg_ca` **point cloud**
(3.0483 A -- an INTERMEDIATE; production's built chain is 3.2041 A, contract rule 1). Verified:
`||e||^2/n` reproduces `rmsd_avg` to 1e-4 on all 126.

The native-free basis is 14 orthonormal directions per target: the pure **scale** mode, the
**consensus** direction (top-75 mean deviation), lane L's **five per-separation-bin distogram
descent directions** -- put in the basis rather than left to the pool's PCs to span -- and 7 pool
PCs. 20 signed reference fields and correlations, 6 unsigned direction features, 17 target
features; ridge, leave-fold-out on the pinned folds with the penalty chosen by a **nested**
leave-fold-out inside the training folds.

**The design decision that makes the number mean anything.** A PC's sign is arbitrary, so
regressing `c_k = <e, U_k>` on target features is ill-posed -- and supplying the per-target sign is
exactly the leverage `in-band-ordering-is-per-target` names as the only one left. Every feature is
therefore either **signed** (flips with `U_k`: projections of the consensus / distogram-descent /
scale / EXPAND / provenance / score-weighted directions, the skewness of the 500 pool projections
along `U_k`, the signed pool correlations of projection against score, Rg, BLOSUM sim and
residual) or **unsigned**, and the model is linear in the signed ones with unsigned interactions
and **no intercept**. The fit is equivariant by construction, so **it cannot launder an ORACLE
sign**. Without this the arm is lane D's rho = 0.949 that was 0.809 dimension.

```
arm                 R2_oof    ctrl     excess    fold CI (excess)     x MDE   verdict
full_nodist|wt     +0.0083  -0.0012   +0.0095   [-0.0079, +0.0988]    0.75    NOT MEASURED
full|wt            -0.0097  -0.0108   +0.0010   [-0.0015, +0.2142]    0.74    NOT MEASURED
main|wt            -0.0021  -0.0024   +0.0002   [-0.0145, +0.1304]    0.68    NOT MEASURED
main_nodist|wt     -0.0020  -0.0023   +0.0003   [-0.0114, +0.0557]    0.44    NOT MEASURED
full_nodist|unwt   +0.0062  -0.0082   +0.0144   [+0.0153, +0.1480]    0.87    NOT MEASURED
main_nodist|unwt   -0.0646  -0.0309   -0.0337                         1.35    WORSE
CPERM (permuted target->feature link)  R2 -0.0388
```

**Registered bar 1.96%. Best arm 0.83%, excess 0.95%, 0.75x MDE. P1 holds; the primary and
ambitious bars are not approached.** In Angstroms, both figures because one of them flatters:
excess rho = 0.098 gives `3.0483*sqrt(1-rho^2)` = **3.0338 A** implied, but **applying the fitted
displacement out of fold actually emits 3.0519 A against production's 3.0483** -- a no-op that is
marginally worse. The implied-endpoint conversion is the generous one and both belong in the
record.

Strata: nothing clears the bar on FAIL18 (circular, flagged), the other 108, or either of lane F's
filter-independent tails; the largest is worst-18-by-pool-mean at R^2 = 0.0462 on n = 18, where the
MDE is far larger than the effect. **P4 holds -- the tail is no more predictable than the bulk in
this space.** P2 holds and is vacuous: the scale mode carries 26.5% of the ORACLE error and its
coefficient's out-of-fold R^2 is +0.0404. P5 holds -- removing every distogram-derived feature moves
the excess by under 1 point, and in the *favourable* direction, i.e. there was nothing there to
inherit.

### 2. THE ONE STRONG POSITIVE, AND IT IS ABOUT LOCATION, NOT DIRECTION

The same basis, scored ORACLE against the matched-dimension random-frame control:

```
ORACLE capture of ||e||^2, native-free basis        0.8273
ORACLE capture, matched-dimension random frame      0.4444
excess                                             +0.3829   10.17x MDE, 5/5 folds, [+0.369, +0.392]
```

> **We can locate 82.7% of the ORACLE error's variance with no native at all, in 14 native-free
> directions, and out of fold we cannot say which way along them (R^2 0.8%). The subspace is free.
> The sign is the whole problem.**

The control is in the operator's own space, and I checked rather than assuming: at n = 9-16 the
rigid-free space is 21-42 dimensional and the 75 pool deviations span **all** of it, so an
isotropic random K-frame is inside the pool's span and this is not an instance of
`control-must-match-the-operators-space`.

### 3. THE PRIOR ROUTE -- WHICH THE DISPLACEMENT BOUND DOES **NOT** PRICE

Lane L flagged, correctly, that `rho_max` governs **pool-space operators** ("add another field to
the 21") and does not price "fix the prior's shape". Section 1 is the operator question; this is a
different one. Every arm translates each pair's **posterior** by `Delta_p` so the shipped
Bayes-risk score is evaluated unchanged at `risk_p(d - Delta_p)`, then runs production's own
downstream: rescore the same 500-member pool, keep 75, coordinate-average in the medoid frame.
Point cloud, production 3.0483:

```
arm                                     mean     delta    x MDE  folds  W/L
ORACLE_FULL      (perfect distogram)   2.2379  -0.8104   3.28    5/5   117/9
ORACLE_SEPPROF5  (5 numbers/target)    2.4743  -0.5740   2.52    5/5   106/20    = 70.8% of full
ORACLE_PERRES    (n numbers/target)    2.4769  -0.5714   2.81    5/5   113/13    -- adds NOTHING over 5
ORACLE_OFFSET1                         2.9041  -0.1443   0.87          69/57     NOT MEASURED
ORACLE_STRETCH1                        2.8983  -0.1501   0.85          74/52     NOT MEASURED
```

**The coordinator's caveat 3 is discharged: S19's -0.525 A predates the benchmark pinning and it
REPRODUCES on the pinned 126, slightly larger, at -0.5740 A.** And the nesting reproduces too --
one offset and one stretch are both NOT MEASURED, the five-number profile is 70.8% of a perfect
distogram, and n per-residue parameters buy nothing beyond the five.

**Which half is missing, and this is the sharpest number in the entry:**

```
ORACLE SIGN      + leave-fold-out magnitude   2.6917  -0.3567   2.29x  5/5  93/33   = 62.1% of the prize
ORACLE MAGNITUDE + leave-fold-out sign        2.9750  -0.0733   0.32x        66/60  NOT MEASURED, 12.8%
```

Matched-accuracy sign corruption, the mandatory null from `error-coherence-decides-correctors`:
acc 1.0 -> -0.574, 0.9 -> **-0.4724**, 0.8 -> **-0.2744**, 0.7 -> -0.2506, 0.6 -> -0.0747 (NOT
MEASURED). **The deployable leave-fold-out sign performs like accuracy ~0.6, barely above chance;
~0.8 is needed to buy -0.27 A.** The target is **five bits per target**.

**Both things I can currently build to supply them deliver zero or worse:**

```
LFO_GLOBALPROF5  (5 global numbers, leave-fold-out)   3.0520  +0.0036  0.11x  NOT MEASURED
NF_POOLPROF5     (5 numbers from the pool)            3.0878  +0.0395  0.48x  wrong sign
```

The mechanism is explicit and is not a fitting failure. The ORACLE profile's mean across targets is
`[0.048, 0.194, 0.281, 0.397, 0.464]` against sd `[0.369, 1.247, 1.712, 2.268, 4.274]`. **The
distogram really does systematically over-predict long-range distances, and the per-target
dispersion is ~9x the systematic part**, so a global constant is arithmetically a no-op. And
estimating the profile **from the pool** makes the answer worse -- `pool-error-is-68-percent-
common-mode` arriving at the prior, and an independent confirmation of lane L's constraint that the
covariate cannot be generated by the pool.

### 4. WHERE THE PRIZE IS -- TWO LOCALISATIONS, BOTH CLEAN

**By separation bin (ORACLE, one bin corrected at a time):**

```
|i-j| =  2   -0.0107   (0.49x, NOT MEASURED)
|i-j| =  3   -0.0874   (1.00x, NOT MEASURED)
|i-j| =  4   -0.0836   (1.32x, 5/5)
|i-j| = 5-6  -0.1071   (1.45x, 5/5)
|i-j| >= 7   -0.3907   (1.98x, 5/5)   = 68% of the profile prize
```

**The prize is long-range.** Lane R measured on a completely different instrument that the RMSD
signal is absent from local features (dR^2 -0.089) and abundant in global ones (+0.600). Two
instruments, one statement.

**By stratum, and it replicates on BOTH of lane F's filter-independent tails**, so it does not rest
on the FAIL18 circularity the adversary exposed:

```
                    other108   tailA(pool mean)  tailB(ORACLE best)  FAIL18 (CIRCULAR)
ORACLE_SEPPROF5      -0.359        -1.405             -1.164             -1.865
ORACLEsign_LFOmag    -0.266        -0.821             -0.591             -0.901
```

**~4x concentration on the hard targets on non-circular definitions.** By the sprint-open
counterfactual that is the stratum where fixing ten targets beats improving all 126 by 0.20 A.

### 5. WHAT I CLAIM AND WHAT I DO NOT

**Claimed.** (a) The operator question is closed **as a measurement**: the whole native-free feature
space -- distogram outputs, pool statistics, geometry, provenance, sequence-level scalars, in the
most favourable basis available and with a sign-equivariant model that cannot cheat -- explains
**0.8%** of the ORACLE displacement out of fold against a 1.96% bar. That reproduces the 21-field
survey's answer in one run instead of an enumeration. (b) The prior route is **not** closed by that
number and must not be read as closed by it. (c) On the prior route the prize reproduces on the
pinned benchmark, is 68% long-range, concentrates 4x on the tail, and reduces to five signs per
target at ~0.8 accuracy.

**Not claimed.** That any covariate exists which supplies those signs. I tested two (a global
constant and the pool) and both are zero or negative. Lane L's design constraint -- *a native-free
covariate, observed at inference, predicting the separation profile, generated by neither the
distogram nor the pool* -- survives this entry untested, and it is now a five-bit question rather
than an open-ended one.

**Two things I got wrong or nearly wrong.** The verdict strings in my first run were **inverted**:
`ST.compare`'s convention is lower-is-better and explained variance is higher-is-better, so an arm
that was 0.0337 WORSE than its control printed `BETTER`. I caught it on the first table and negated
both arms. It is the same failure lane D found in S29-L23's printed verdict, one sprint later, in a
lane that had just read the entry about it. And my original brief -- hunt for a decorrelated source
-- would have been a waste: lane L's arithmetic (a perfectly orthogonal channel must itself carry
rho = 0.3398, 3.01x anything we own; 10.1 orthogonal channels to reach 3.00 A against a measured
stable rank of 2.057) prices decorrelation out before any measurement, and my section 1 confirms it
from the other side.

## S30-L26 -- TWO ANSWERS, AND **MY REGISTERED PRIOR IS WRONG ON BOTH**. (1) THE AVERAGING-ARTIFACT LITERATURE'S UNMEASURED CLAIM IS **CONFIRMED ON OUR INSTRUMENT, BOTH HALVES**: pool divergence -> corrupted average at rho **+0.947**, and divergence -> relax value at rho **-0.316** (permutation null p = 0.000, **length-matched** split **-0.0406 at 3.56x MDE, 5/5 folds**), with **97.2% OF THE GAIN IN THE DIVERGENT HALF** AND THE OTHER HALF AT **-0.0012** -- BUT **DIVERGENCE IS NOT THE PRODUCTION TAIL**, SO E2 IS **NOT** A TAIL INTERVENTION: FAIL18 GETS **-0.0113 AGAINST THE 108's -0.0239**. (2) Q2 IS **NO**, AND IT IS A **THEOREM PLUS A MEASUREMENT WITH THE EASY EXPLANATION REMOVED**: BY CLASSICAL MDS EVERY REFLECTION-INVARIANT SINGLE-STRUCTURE CHANNEL **IS** A DISTANCE-MAP READING, SO THE THREE CLOSED BUCKETS ARE COMPLETE UP TO **CHIRAL FUNCTIONALS** -- WHICH I BUILT, WHICH ARE **GENUINELY EXERCISED** ON OUR MANIFOLD (OCCUPANCY **0.69-0.97 AGAINST DIS's 0.23-0.35**, SO THE CLASS IS **NOT** EMPTY IN PRACTICE, REFUTING MY OWN REGISTERED MECHANISM), AND WHOSE BEST ANCHOR CONTRAST IS **+0.0405 AGAINST A MAX-OVER-3 NULL MEAN OF +0.0408, p_max 0.430** (2026-09-20 13:52, G)

Prereg `s30/PREREG_S30_G.md` @ **78b65521**, committed before the first number. Code and rows @ 44cfc9f7.
Lane R's floor travels with every Q2 sentence: **the torsion rebuild is 0.347 A, not 0** -- reproduced
exactly by my run, which is the check that the ladder is bit-identical.

---

# Q1 -- WHERE DOES THE k=30 AMBER-RELAX GAIN LIVE?

**Artefact, verified before planning (rule 13).** `s16/results/repair_A.json` (126 per-target rows)
and `s16/results/repair_report.json`, both in git at a15406c8.
`settings.k30_full.vs_proj_ungated` = `mean_diff -0.0221, se 0.0070, n 126, W/L 67/59` -- exactly the
memory's `-0.022 [-0.036,-0.009]`. **Nothing was recomputed.** `d_t = rmsd(f0_k30_full) - rmsd(proj)`.

## THE POWER NOTE, REGISTERED BEFORE THE SPLIT -- AND THE SPLIT BEAT IT

```
sd(d_t) = 0.0789        half-split MDE 0.0394      FAIL18-vs-108 MDE 0.0563
whole sample            -0.0221  fold CI [-0.0297,-0.0159]  1.12x MDE  5/5 folds
```

I registered that **a half-split cannot reach its own MDE unless MORE than 100% of a 0.022 A effect
sits in one half**, and therefore that the continuous statistic was primary and a null would mean
*"we could not localise it"*, never *"it is uniform"*. **More than 100% of it does sit in one half.**

## F-G1 FIRES ON EVERY FORM

```
                    rho(d, disp)   perm-null sd   xMDE    p2      split hi-lo      fold CI              folds
DISP_rmsd             -0.316          0.089      -1.27   0.000     -0.0417   [-0.0620,-0.0222]          5/5
DISP_S                -0.304          0.089      -1.22   0.001     -0.0420   [-0.0597,-0.0188]          4/5
DISP_rg               -0.355          0.088      -1.43   0.000     -0.0564   [-0.0658,-0.0462]          5/5

high-dispersion half  -0.0429        low half  -0.0012        share of gain in the high half  97.2%
uniform-effect null (rule 5): observed gap at pctile 0.001 of [-0.0264,+0.0289]   (all three variables)
```

The null is **fold-preserving** (labels permuted within fold) so fold structure in either variable
cannot manufacture it. **Multiplicity is ~1 comparison, not 3**: the three dispersion variables
correlate at **0.93 to 1.000** with each other, and `DISP_rg` reproduces lane F's own
`F4_top75_rg_sd` at **rho = 1.0000**.

## I THEN TRIED TO DESTROY IT, FOUR WAYS

**(a) Is it the ten winners relabelled?** 10 targets carry 66.8% of the whole-sample gain, so this
was the real risk.

```
              whole sample   hi-lo split    fold CI              xMDE   folds   rho
drop top 0      -0.0221        -0.0417   [-0.0620,-0.0222]     -1.48    5/5   -0.316
drop top 5      -0.0144        -0.0407   [-0.0574,-0.0240]     -1.73    5/5   -0.292
drop top 10     -0.0080        -0.0310   [-0.0509,-0.0161]     -1.19    5/5   -0.223
drop top 20     +0.0028        -0.0095   [-0.0229,+0.0044]     -0.46    3/5   -0.104
```

**It survives removing the ten biggest winners** and dies only at drop-top-20 -- where the arm itself
has **no gain left to localise** (+0.0028). So it is a gradient, not a relabelling; and the gradient
is exactly coextensive with the effect. The top-10 winners' median dispersion percentile is 0.833,
but 3 of the 10 sit below the median, so they are not the high-dispersion half by another name.

**(b) Chain length.** `rho(n, DISP) = +0.267`, so this had to be checked. Partialling n **strengthens**
the effect (-0.316 -> **-0.338**, p2 0.000), `rho(n, d) = +0.036`, and the **length-matched split** --
high vs low dispersion *within* chain-length tertiles -- is **-0.0406 [-0.0481,-0.0319], 3.56x MDE,
5/5 folds**, the strongest form of the result. Not length.

**(c) Mediation.** `rho(DISP, contraction) = +0.947` and `rho(DISP, raw-average Ramachandran-ok) =
-0.763`. Partialling contraction leaves `rho(DISP, d) = -0.178`; partialling DISP collapses
`rho(contraction, d)` from -0.276 to **+0.078**. **Divergence is the primary variable and contraction
is its shadow**, not the reverse.

**(d) The stratum.** Declared in the prereg: FAIL18 is defined by the *filter's* in-band recall
(S30-L23), not by `d_t`, so this is not the S30-L8 failure mode -- but it is reported as a
**descriptive** split, and fold 0 contains no FAIL18 target.

# AND THIS IS THE PART THAT DECIDES THE COORDINATOR'S QUESTION: **DIVERGENCE IS NOT THE TAIL**

```
rho(DISP, production RMSD)             +0.514
rho(DISP, pool mean)  ORACLE           +0.699
FAIL18 dispersion vs the 108           +0.4583  [-0.3845,+1.0472]   0.47x MDE   NOT MEASURED
ORACLE worst18-by-pool-mean dispersion +1.6116  [+0.5832,+2.0755]   1.44x MDE   (ORACLE)
FAIL18 in the high-dispersion half     11 of 18   (9 expected if independent)

relax gain on FAIL18    -0.0113        relax gain on the other 108    -0.0239
ORACLE worst18-by-pool-mean            -0.0192   0.31x MDE   NOT MEASURED
```

> The brief's inference was *"if it concentrates on divergent pools it is a tail intervention, and
> the tail is where the prize is."* **The premise holds and the conclusion does not**, because the
> two tails are different objects. Divergent pools are harder on an ORACLE definition (+1.61, 1.44x
> MDE) but **FAIL18 is not measurably more divergent** (0.47x MDE), and the gain is **smaller on
> FAIL18 than on the other 108** -- the wrong direction.

**ENDPOINT ARITHMETIC.** The whole arm is **-0.0221 A = 0.69% of the 3.2105 baseline**; applying it
only to the divergent half is worth **-0.0215**, i.e. the targeting buys nothing extra because the
other half was already contributing nothing. Against the opening counterfactual's ask of **-0.30 A**
from capping the worst 10 at 3.00 A, **E2 is not on that scale and does not reach that stratum.**

**VERDICT. The literature's claim is CONFIRMED -- and, as far as lane L's search found, measured here
for the first time, on both halves: divergence -> artifact (rho +0.947) and artifact -> repair value
(rho -0.316).** But **escape E2 does not become a tail intervention**, because pool divergence and the
production failure tail are not the same set. E2 remains what S30-L13 priced: a real, mechanistically
explained, divergence-graded 0.7%-of-baseline effect whose restraint constant still has **no
native-free selection rule** (`averaging-space-beats-the-objective`: the non-circular rule picks
k=10 at +0.001).

**MY PRIOR WAS WRONG AND THE REASONING THAT PRODUCED IT WAS TOO.** I registered 2-to-1 *against*
concentration, arguing from `corr(S,B) = +0.085` that spread is orthogonal to bias and from the 68%
common-mode result that averaging's residual is the shared component. Both are true and **neither
governs this operator**: the relax does not repair the bias, it repairs the *geometry* the averaging
destroyed, and the destruction scales with spread at rho 0.947. **I am the second lane to revise lane
L's "it concentrates" downward on independent reasoning -- L went to even odds on lane F's
dispersion null, I went to 2-to-1 against on the common-mode argument -- and the convergence was
wrong in both directions at once.** Lane F's null was about a different outcome (filter failure);
mine was about a different mechanism (bias vs geometry). Two lanes agreeing does not make a prior.

---

# Q2 -- CAN A NATIVE-FREE CHANNEL WITH GENUINE GLOBAL REACH BE CONSTRUCTED AT ALL?

## THEOREM G1, THE CHIRALITY DICHOTOMY (registered before any number)

For a single-structure channel `S: R^{n x 3} -> R` invariant under rotation and translation, **with
the sequence held fixed as a parameter** (reflection does not touch the sequence, so burial, SASA and
hydrophobicity-weighted terms are covered):

> **`S` is a function of the pairwise distance matrix `D` IFF `S` is also invariant under
> reflection.**
>
> *Proof.* Classical MDS gives `G = -1/2 J D^2 J = X X^T`, so `D` fixes the centred coordinates up to
> `Q in O(3) = SO(3) x {+-I}`. After quotienting the rotations this project actually uses -- proper
> rotation Kabsch with the det correction, `s12/instrument.py:kabsch_rmsd_batch` -- the residual
> ambiguity is **exactly reflection**, so a reflection-invariant `S` is constant on the fibre and
> descends to a function of `D`. Conversely `D` is reflection-invariant. **QED**

**COROLLARY -- the brief's whole candidate list falls on one ground, not one at a time.** Long-range
contact topology, the radius-of-gyration **profile**, principal-axis / inertia-tensor structure
(asphericity, acylindricity, shape anisotropy), end-to-end and intermediate-separation distance
distributions, lane L's separation profile, contact order, excluded volume, packing density and
burial/SASA are **each reflection-invariant** -- a mirrored peptide has identical contacts, identical
inertia spectrum, identical separation profile, identical solvent-accessible surface. By G1 every one
is a function of `D`. **They are not new channels; they are coordinate systems on the distance map.**

**AND THE REFERENCE, NOT THE FUNCTIONAL, IS WHERE THE INFORMATION LIVES.** A function of `D` is a
shape descriptor; to score *nativeness* it must be compared to an expectation, and this project has
exactly three sources of one -- which **are** the brief's three closed buckets:

| reference | bucket | why closed |
|---|---|---|
| the predicted distogram | (a) class M | bounded by theorem 2 |
| the pool | (b) consensus | typicality; anchor contrasts CONS -0.056, DMAP_CONS -0.053, POOLGO -0.050 |
| universal physics (AMBER, DOPE/LEG, Ramachandran, a hydrophobicity table) | -- | **target-independent by construction**, so it cannot supply the per-target sign that `in-band-ordering-is-per-target` says is the only leverage -- and it is exactly what lane R's ladder matches rung for rung |

> **So the three buckets are complete up to EXACTLY ONE family: CHIRAL functionals.** This also
> *explains* lane R's split verdict rather than merely surviving it: all 43 of its channels are
> either per-residue sums (blind by D1) or `D`-functionals against one of those three references.
> **There was no fourth kind in the library to test.**

## THE PRE-CHECK -- AND IT REFUTES MY OWN REGISTERED MECHANISM

I predicted F-G2 would fail because at n = 9-16 the chiral coordinate would be *near-constant*: local
helical handedness matched by the fold's Ramachandran draw. The coordinator required that be measured
rather than inferred. `occupancy = sd(X)/sqrt(sd(X)^2+mean(X)^2)`, which is exactly
`sd(X over candidates) / sd(X over candidates UNION their mirrors)` -- 1.0 means the chiral axis is
fully exercised, 0 means the set sits on one side of it:

```
channel          set        mean       sd     |mean|/sd   occupancy   frac>0
WRITHE           pool      0.5903   0.5363      1.105       0.685      0.782
CHIRAL3          pool      0.0352   0.1636      0.262       0.953      0.669
CHIRAL3_LONG     pool      0.0276   0.1701      0.218       0.965      0.661
DIS              pool      2.4870   0.8730      3.610       0.345      1.000
DIS              ladderA   2.4310   0.5656      4.853       0.231      1.000
```

> **THE ESCAPE CLASS IS NOT EMPTY IN PRACTICE. The chiral coordinate is more fully exercised on our
> manifold than the shipped cost is** (0.69-0.97 against DIS's 0.23-0.35), and both signs are
> occupied (frac>0 = 0.62-0.78). **My registered mechanism is refuted**, which makes the negative
> below *stronger*, not weaker: the easy explanation is gone.

## F-G2 DOES NOT FIRE -- ON LANE R's OWN LADDER, CONTROLS AND STATISTICS

Bit-identical ladder (same `Sampler`, same `crc32` seed rule, same call order); the **0.347 A rebuild
floor reproduces exactly**, which is the check that it is the same instrument.

```
channel          rho_A_part  rho_ANCHOR  A-ANCHOR   fold CI             xMDE  folds
WRITHE             -0.207      -0.247     +0.0405  [-0.0255,+0.1015]   +0.41   3/5
CHIRAL3            -0.044      -0.008     -0.0363  [-0.0933,+0.0147]   -0.32   3/5
CHIRAL3_LONG       -0.038      +0.001     -0.0391  [-0.0982,+0.0139]   -0.35   3/5
DIS                +0.349      +0.225     +0.1238  [+0.0800,+0.1666]   +1.42   5/5

max-over-3 sign-flip null (lane R's statistic):  observed 0.0405   null mean 0.0408   p95 0.0852
                                                 p_max = 0.430
```

Registered bar was **+0.10** with a fold CI excluding zero. Best chiral channel: **+0.0405, CI
includes zero, 3/5 folds, 0.41x MDE** -- and **sitting on its own null to three decimal places.**

**THE MATCHED CONTROL, AND IT IS WHAT MAKES THIS A CLOSURE RATHER THAN A NULL.** Each chiral channel's
**achiral twin `|X|`** is reflection-invariant, hence by G1 a distance-map reading; same functional,
same scale, same ladder, differing in **nothing but chirality**:

```
                 ordering          anchor contrast     preference
WRITHE          +0.0054 (0.11x)    +0.0170 (0.29x)    +0.0471 (0.53x)   all CIs include zero
CHIRAL3         -0.0760 (-1.02x)   -0.0160 (-0.15x)   -0.0381 (-0.30x)
CHIRAL3_LONG    -0.0808 (-1.06x)   -0.0141 (-0.13x)   -0.0372 (-0.29x)
```

> **Whatever WRITHE can do, its own reflection-invariant shadow already does. The chiral content
> contributes nothing measurable** -- and the two purely non-local chiral channels are *worse* than
> their twins at ordering, with CIs excluding zero and 5/5 folds.

**AUDIT, asserted in code and stored in the artefact rather than claimed in prose** (worst over all
126 targets): reflection flips `X` at error **0.00e+00**, leaves `|X|` at **0.00e+00** and leaves `D`
at **0.00e+00**; rotation invariance **4.19e-13**. The audit also **caught a real bug**: my first
Klenin-Langowski sign term was not rotation-invariant (rot_err 2.03) and the assertion stopped the
run before any number existed.

## A POSITIVE I FOUND, CHASED, AND KILLED MYSELF -- AND IT IS THE ENTRY'S METHODOLOGICAL OUTPUT

`contrast_pref = pref_near - pref_pool` for WRITHE is **+0.1641 [+0.1016,+0.2282], 2.17x MDE, 5/5
folds, max-over-3 null p = 0.000**, with `pref_near = 0.771` -- which clears **both** of lane R's
registered preference clauses (bar 0.65, margin 0.10) that **none of its 43 channels cleared**.

**It is cross-kind, in exactly the way S30-L1 withdrew S28-L48 for.** `pref_pool`'s control set is the
500 pool members, which **keep their DEPOSITED coordinates** (`W[:k] = cand.W`), while the near-native
rungs are **ideal-geometry rebuilds**. A raw geometric shape statistic separates those two
constructions whether or not it sees nativeness. Lane R's ladder matched kind for the **ordering**
contrast; the **preference** contrast is against production and the pool, and is not kind-matched for
a channel of this type.

**The internally matched statistic settles it** -- the rebuilt native's percentile inside its **own**
ladder, every member an ideal rebuild from the same table at the same budget, chance exactly 0.5:

```
WRITHE        0.6061 [0.5609,0.6522]  5/5   WORSE than chance
ABS_WRITHE    0.6732 [0.6556,0.6890]  5/5   WORSE than chance
CHIRAL3       0.5049 [0.4782,0.5368]  5/5   at chance
CHIRAL3_LONG  0.5020 [0.4760,0.5361]  5/5   at chance
DIS           0.2876 [0.2617,0.3144]  5/5   better than chance
```

> **WRITHE ranks the native at the 61st percentile of the native's own perturbations.** A channel
> that does that is not recognising nativeness; it is separating production from everything else.
> The +0.164 is withdrawn by me before it was ever quoted. `|WRITHE|` clears the same bars
> (`pref_near` 0.724) and is *worse* still on the matched read (0.673), which is the tell.

One small real thing in the other direction, reported because it is in the table: `ABS_CHIRAL3` and
`ABS_CHIRAL3_LONG` put the native at **0.443** and **0.435**, CIs excluding 0.5, 5/5 folds -- better
than chance. They are **achiral**, hence `D`-functionals inside the already-closed class, and their
anchor contrasts are the wrong sign (-0.020, -0.025). Worth a line, not a claim.

## THE ANSWER

> **NO.** Every achiral native-free single-structure channel **is** a distance-map reading (G1); the
> only three references available to score one against are the three closed buckets; and the single
> family the theorem leaves open was built, is **genuinely exercised on the manifold we occupy**, and
> carries **no nativeness signal** -- its best anchor contrast sits on its own max-over-channels null
> (+0.0405 vs +0.0408, p_max 0.430), its entire apparent skill is reproduced by its own achiral
> shadow, and it ranks the native **worse than chance** inside its own matched ladder.
>
> **This completes the sprint's ceiling argument rather than adding to it**: lane R's "the null for
> every per-residue channel is a theorem" now has a companion -- *and the null for every non-chiral
> channel is a theorem too*, leaving one family, which is measured and empty.

**SCOPE, stated so it is not over-quoted.** G1 bounds **single-structure** channels. It does not bound
set-referenced channels -- those are bucket (b) and are closed separately. It does not bound a channel
with access to a **fourth reference this project does not have**, e.g. an experimental observable;
`torsion-restraints-reach-the-target` already prices that route. And the measurement is at peptide
length: the non-local chiral content of a 13-mer is small, so G1's *survivor* family may be worth
testing again if this project ever works at 40+ residues -- the **theorem** does not depend on length,
only the emptiness does.

**MULTIPLICITY.** 3 registered chiral channels + 3 twins + DIS = 7 channels x 2 contrasts. The
max-over-3 sign-flip null is reported for both contrasts. The one comparison that cleared its bar is
withdrawn above on a kind argument, not on its p-value.
