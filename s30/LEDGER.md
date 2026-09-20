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
