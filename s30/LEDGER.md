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

## S30-L3 -- FILTER WIDTH IS **NOT** AVERAGING WIDTH: IT MOVES THE EMITTED STRUCTURE 4-15× AS MUCH, SO S29'S FLAT m SWEEP DOES NOT TRANSFER. BUT THERE IS **NO FREE LUNCH IN WIDTH** -- F2a REFUTED, THE BENEFIT/HARM TRADE-OFF IS CLEAN AND CONTINUOUS AND CROSSES AT ~55%/55%. F2b CLOSED **BY CEILING WITHOUT SPENDING THE CHAIN**: THE ORACLE GLOBAL BEST FILTER WIDTH **IS THE SHIPPED VALUE**. AND THE "WIDENING HELPS HARD TARGETS" READING IS **FALSE** -- IT DOES NOT REPLICATE ON EITHER FILTER-INDEPENDENT TAIL. ONE MECHANISM (ρ = +0.81) EXPLAINS ALL FOUR STRATA (2026-09-20 12:57, F)

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
