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
