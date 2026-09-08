# SPRINT 17 — AUDIT workstream findings

Modules: `s17/audit_ordernull.py`, `s17/audit_redundancy.py`, `s17/audit_repro.py`,
`s17/audit_frame.py`, `s17/audit_strat.py`.
Artefacts: `s17/results/audit_ordernull.json`, `audit_dupe.json`, `audit_clusters.json`,
`audit_repro.json`, `audit_repro_meta.json`, `audit_frame.json`, `audit_exactrand.json`,
`audit_strat.json`.

Every number below was **re-derived from `s8/generate_univ/*.npz`**, not read out of another
workstream's artefact. Where a coordinator number is quoted as reproducing, it was recomputed
by a second route (different dtype, different code path) and the max |Δ| is stated.

---

## 0. VERDICT TABLE — worst first

| # | claim under attack | verdict |
|---|---|---|
| **A1** | "the 1.711 Å K=500 oracle is a truncation artefact; the real ceiling is 1.313 Å" — read as *there is better structure out there that wider retrieval would find* | **NUMBERS EXACT, READING REFUTED.** The 0.397 Å gain is *below* the 0.473 Å a random ordering of the same universe delivers: obs − matched null **−0.075 Å [−0.132, −0.022]**. It is the order statistic of more draws, and BLOSUM under-delivers against it. |
| **A2** | oracle_map §C, `raw coordinate average, K = 500 → 3.396`, printed under `incumbent (deployed) 3.204` | **INVERTS THE SPRINT'S CENTRAL ARCHITECTURAL COMPARISON.** The programme's standing "raw coordinate average 3.050" is the average over the **shipped top-75** (recomputed: **3.0483**). The map's row averages the **BLOSUM top-500** — an operator nobody deploys, at a width the pipeline never uses. Averaging reads 0.19 Å *worse* than the incumbent where it is in fact 0.156 Å *better* (the incumbent is that average **plus** the +0.155 Å repair tax). |
| **A2b** | oracle_map §D, "widening K raises the ceiling; ranking is the blocker" | **THE NUMBER IS RIGHT AND STRONGER THAN THE SENTENCE.** Recovered is not "near zero", it is **−8% to −18%**: the realized answer gets **worse** with K (+0.141 Å at full vs K=75, fold-clustered CI [+0.046,+0.253]). Priced against the *exact* random control the selector is net −0.388 Å [−0.486, −0.268] — both halves must be stated. |
| **A3** | oracle_map table A prints `mean`, `diversity`, `coordavg` side by side | **LIVE SPRINT-16 TRAP.** `mean` is free-superposition, `diversity` is common-frame. The sprint's central identity is **exact** with the common-frame member error (residual 7×10⁻¹⁴ Å²) and off by **−2.10 Å²** with the printed column. |
| **A4** | oracle_map's matched random control | **UNDERPOWERED CONTROL.** `N_RAND = 3` ⇒ per-target SE **0.763 Å** (max 1.237). Its W/L counts and table-B interval are noise. The exact expectation is free (pool mean). |
| **A5** | "K = 500 pool" is a BLOSUM ranking | **PARTLY FALSE, HARMLESS.** 11.4% of the pool is chosen by corpus iteration order inside a tie group — but random tie-breaks give 1.708 Å vs the shipped 1.711 Å, per-target sd 0.018 Å. **CLEAN SURVIVAL.** |
| **A6** | "a universe of 13,000–27,000 windows per target" | **A BOOKKEEPING QUANTITY.** corr(U, chain length) = **−0.988**. U is 1/n. The exemplar 1CB3 is the shortest target with the biggest U. |
| **A7** | this audit's own hypothesis that the deep universe is near-duplicates | **MY H0 FAILS, HONESTLY REPORTED.** d* = 1.355 Å [1.272, 1.458]; the deep windows are genuinely distinct structures. The universe is a dense continuum, not copies. |
| **A8** | do the map's core numbers reproduce? | **YES, EXACTLY.** 3.4540 / 1.7108 / 1.3134; two RMSD implementations agree to 4.8×10⁻⁷ Å over ~2.1 M windows. |
| **A9** | leakage | **CLEAN.** No `s17/` module reads benchmark60 or `bench_results`; no bare `hash()`. |
| **A10** | `s17/results/oracle_map.json` as a shared file | **HAZARD.** No completion flag, unlike `I.write`'s `complete`/`n_expected`. Three workstreams read it while it held 60/126 rows with fold 3 under-represented 2.7×; nothing in the file said so. Now complete, hazard recurs on rerun. |
| **A14** | oracle_map §C, `ORACLE best in top-75 → 2.104` vs BRIEF §7's `shipped top-75 2.306` | **TWO SETS, ONE NAME.** They share 13.2/75 members. "The top-75 truncation costs 0.791 Å" attributes to *retrieval truncation* a loss the pipeline incurs by *selection* — the problem-A/problem-B distinction BRIEF §3 makes mandatory. |
| **A11** | power of the 126-target instrument | SE **0.030 Å**, MDE at 80% power **0.084 Å**. This instrument can detect what the sprint is denying. |
| **A12** | `selection_theory.rho_from_accuracy` | **SURVIVES.** ρ = sin(π(2·acc−1)/2) is the correct Gaussian-copula τ↔ρ map and s14's 0.600/0.638 really are pairwise ordering accuracies. One caveat: s14's matched-pair null is 0.505, not 0.500. |
| **A13** | `readout.py`'s error/diversity decomposition | **SURVIVES.** Its `err` is the common-frame quantity; `rmsd² = err² − div²` closes to 10⁻¹⁵ on its own artefact. It does correctly what oracle_map table A does not. |

---

## A1 — THE SPRINT'S OPENING PREMISE. Exact arithmetic, inverted meaning.

`s17/BRIEF.md` §2 and `oracle_map.py`'s docstring both state that the programme's 1.711 Å
ceiling "is not the retrieval ceiling — it is the ceiling of the first 500 entries of a
BLOSUM62 ranking over a universe of 13,000–27,000 windows per target," and offer 1CB3
(2.028 Å → 1.096 Å) as the exhibit. **Both numbers are right.** I reproduce
best(500) = 1.7108 Å and best(full) = 1.3134 Å independently (§A8).

The reading the sprint has taken from them is that truncation is *costing* the programme
0.397 Å of attainable accuracy. That reading requires the extra windows to be *better*, not
merely *more*. **The minimum of N draws falls with N even when nothing new is present**, so
the statistic needs a matched null: the same best-of-K statistic computed under a **random
permutation** of the identical universe (`audit_ordernull.py`, 200 permutations/target,
`stable_rng`).

```
                                              mean      95% CI (fold-clustered)   median
observed   best(500) − best(full)             0.397     [0.362, 0.429]            0.311
MATCHED NULL, same statistic, random order    0.473     [0.430, 0.511]            0.426
observed − null                              −0.075     [−0.132, −0.022]         −0.021
                                              W/L 56/70
```

The ceiling gain from widening K is **not larger than the order statistic — it is
significantly smaller.** The interpretation is not "retrieval is broken"; it is the
opposite, and it is the sharper statement: **BLOSUM retrieval has already front-loaded the
good structures into the first 500, so the marginal return on the next 12,500 is below the
i.i.d.-draw rate.** Widening K therefore buys the *tail* of a distribution the pool is
already sampling from the good end of.

**The retrieval ordering does have skill, and it is small.** `audit_ordernull.py` §A,
skill(K) = best-of-K under a random order minus best-of-K under the BLOSUM order:

```
      K     BLOSUM    random     skill    [95% CI]           W/L
     75      2.104     2.212     0.108   [+0.031,+0.179]    71/55
    500      1.711     1.786     0.075   [+0.022,+0.135]    70/56
   2000      1.504     1.553     0.049   [+0.019,+0.074]    84/42
   5000      1.410     1.429     0.020   [+0.001,+0.036]    83/43
```

Skill is positive at every K (so my own H0 — "the ordering is uninformative" — is falsified
as stated) but it **decays monotonically to zero**, which is exactly why the incremental
ceiling under-delivers against the null.

In recall terms, with the matched-K control (BLOSUM top-500 vs a random 500 of the same
universe, 200 draws/target — `audit_strat.recall_ctrl`):

```
  t = 1.5 Å :  BLOSUM-500 53/126   random-500 48.2/126   retrieval buys +4.8 targets
  t = 2.0 Å :  BLOSUM-500 73/126   random-500 71.5/126   retrieval buys +1.5 targets
  t = 2.5 Å :  BLOSUM-500 99/126   random-500 98.1/126   retrieval buys +0.9 targets
```

**The whole sequence-retrieval apparatus is worth 1.5 targets of 2.0 Å recall over drawing
500 windows at random from the same library.**

### The stratification that kills the operational reading (`audit_strat.py` §C)

```
                          n     observed gain   null gain   obs − null
  n 9–11                 33         0.496         0.529       −0.033
  n 12–13                42         0.396         0.364       +0.032
  n 14–16                51         0.335         0.526       −0.191
  K=500 best  > 2.0 Å    53         0.638         0.592       +0.046
  K=500 best <= 2.0 Å    73         0.222         0.386       −0.164
  FAIL18                 18         0.644         0.657       −0.013
```

**On the 53 targets that need it and on all 18 FAIL18 targets, the widening delivers exactly
the order-statistic amount and nothing more.** The apparent aggregate structure comes from
the *easy* targets, where widening under-delivers. This is the "a mean improvement carried by
easy targets is not a result" rule firing in reverse.

### Held to the sprint's own reporting standard (`audit_strat.py` §B)

```
mean −0.075 [−0.132,−0.022]   median −0.021   W/L 70/56   sd 0.335
drop the 10 most negative targets: −0.006     drop 20: +0.032
UNIFORM-EFFECT NULL for drop-top-10: [−0.084, +0.035]  → observed −0.006 is INSIDE
per fold: f0 −0.038(25)  f1 +0.016(23)  f2 −0.070(25)  f3 −0.188(23)  f4 −0.095(30)
```

Near-even W/L with a CI excluding zero is the programme's own early-warning signature, so the
drop-top statistic was compared against a **uniform-effect null** rather than a raw threshold
(the standing memory records that a raw threshold misfired here once). It lands **inside** the
null: the effect is uniform, not concentrated. Fold 1 has the opposite sign, and that is
reported rather than smoothed.

**Claim label: SUPPORTED (n=126, target-level, fold-clustered CI, matched null, uniform).**
**Consequence for the sprint: §2 of the BRIEF should not be read as "the ceiling is 1.313 Å
and truncation is costing us 0.397 Å". It should be read as "the K=500 pool is already
near the efficient frontier of this library's sampling, and 1.313 Å is what min-of-13,000
returns from a distribution the pool samples well."**

---

## A2 — THE SPRINT'S BUDGET TABLE PRINTS AN AVERAGING OPERATOR NOBODY DEPLOYS, AND THE
## SIGN OF THE SPRINT'S CENTRAL ARCHITECTURAL COMPARISON FLIPS

`oracle_map.report()` §C is titled "WHERE THE ANGSTROMS ARE" and is the table that sets the
sprint's budget. Its first two rows are:

```
  incumbent (deployed)                                3.204
  raw coordinate average, K = 500                     3.396
```

`BRIEF.md` §7's standing table carries **"raw coordinate average, before repair | 3.050"**,
and `BRIEF.md` §1 makes the whole sprint an argument for **switching from an averaging
readout to a selection readout**, with "coordinate averaging remains an explicit baseline,
never a silent partner."

Recomputed over all 126 targets:

```
coordinate average over the SHIPPED top-75   3.0483     <- this is the standing 3.050
coordinate average over the BLOSUM top-75    3.2928
coordinate average over the BLOSUM top-500   3.396      <- this is the map's row
```

The map's row is a **different operator on a different set**: the mean of the first 500
entries of a sequence ranking, at a width the deployed pipeline never uses and over a set the
deployed shortlist shares 13/75 members with (§A14). Placed directly under the incumbent, it
makes coordinate averaging look **0.19 Å worse than the incumbent**.

The deployed averaging readout is **3.048 Å, i.e. 0.156 Å BETTER than the incumbent** — and
the two reconcile exactly through the programme's own repair tax: 3.048 + 0.155 = 3.203 ≈ the
incumbent 3.204. **The incumbent *is* the coordinate average plus repair.**

**Consequence.** The sprint is about to argue that the averaging readout is capped and should
be replaced by selection. That argument must be made against 3.048, not against 3.396. As
printed, §C's table understates the baseline the sprint has to beat by 0.35 Å, and states the
sign of "averaging vs the incumbent" backwards. **This is the single most consequential
mislabel I found, because it sits in the table that sets the sprint's target.**

Recommendation: print both rows with their operators named —
`coordinate average, shipped top-75 (the deployed readout)  3.048` and
`coordinate average, BLOSUM top-500 (diagnostic)  3.396`.

---

## A2b — oracle_map §D's NUMBER IS STRONGER THAN ITS SENTENCE, AND THE RIGHT CONTROL
## CHANGES BOTH HALVES

`oracle_map.report()` §D prints "recovered", "the fraction of the extra CEILING that the
distance selector actually converts into realized accuracy", and warns that "**near zero**
means widening K buys a bigger ceiling and nothing else."

It is not near zero. At n = 126 the map's own recovered column runs **−8.4% to −17.9%**, and
the realized column is positive at every K. Recomputed from the map's artefact and priced
against the **exact, noise-free** random-selection expectation (`audit_frame.py` §C —
E[uniform pick from top-K] is the mean of `rr` over the top-K in closed form, so no draws are
needed at all), n = 126, fold-clustered:

```
      K   ceiling Δ   realized Δ        [95% CI]      random drift   net of drift
    150     −0.230      +0.031   [−0.002,+0.067]        +0.062   −0.031 [−0.074,+0.019]
    300     −0.322      +0.052   [+0.020,+0.086]        +0.125   −0.074 [−0.137,+0.002]
    500     −0.393      +0.033   [−0.008,+0.073]        +0.167   −0.134 [−0.200,−0.080]
   1000     −0.536      +0.062   [−0.009,+0.137]        +0.235   −0.173 [−0.262,−0.066]
   2000     −0.600      +0.099   [+0.033,+0.181]        +0.310   −0.211 [−0.305,−0.077]
   5000     −0.695      +0.114   [−0.010,+0.215]        +0.413   −0.300 [−0.419,−0.174]
   full     −0.791      +0.141   [+0.046,+0.253]        +0.529   −0.388 [−0.486,−0.268]
```

Two corrections, and they point in opposite directions — both must be stated:

1. **Widening K makes the deployed answer WORSE, significantly.** +0.141 Å at the full
   universe against K = 75. **Truncation is a filter that is helping**, which is the reverse
   of the framing BRIEF §2 sets up. *Interval caveat, stated rather than hidden:* the
   fold-clustered interval [+0.046, +0.253] excludes zero while `oracle_map`'s own plain
   target-level bootstrap gives [−0.008, +0.293], which does not. The two schemes disagree at
   the boundary, so the honest label is **SUPPORTED under fold clustering, SUGGESTIVE under
   the plain target bootstrap** — and the direction is consistent at every one of the seven K.
2. **But the selector is not idle, and §D cannot see that**, because it has no baseline for
   the candidate set's composition drift. Widening from 75 to the full universe degrades the
   *expected random pick* by **+0.529 Å**; the distance selector gives back only +0.141 Å of
   that. Net of drift the selector gains **−0.388 Å [−0.486, −0.268]** from the wider set. A
   "recovered %" computed without this baseline is BRIEF §44's error — an effect measured
   against a baseline that is itself moving, and here it moves 3.7× the size of the effect.

---

## A3 — THE FRAME ERROR IS LIVE IN THE SHARED INSTRUMENT

`oracle_map.report()` table A prints on one row:

| column | what it is |
|---|---|
| `mean` | mean of `I.kabsch_rmsd_batch(W, nat)` — **free superposition**, every window separately optimally rotated onto the native |
| `diversity` | RMS deviation of members about the coordinate average, in the medoid frame — a **common-frame** term |
| `coordavg` | RMSD of the average to the native |

The sprint's central identity (`BRIEF.md` §1) is `readout² = mean member error² − diversity²`.
That is the parallel-axis theorem, and it holds **only** when the member errors are measured
in the same fixed frame as the readout. Sprint 16 paid ≈2× attribution error for exactly this
substitution; the shared instrument of Sprint 17 puts the wrong column next to the right one.

`audit_frame.py` §F, n = 30 targets, residual of the identity in Å²:

```
    K   readout    div   err_common  err_free   resid common   resid free
   75     3.370   3.204     4.712      4.506      6.3e-17       −2.101
  500     3.436   3.353     4.881      4.671     −1.3e-15       −2.156

  K=75  : max |resid common| 7.1e-14 Å² (EXACT)   max |resid free| 7.852 Å²
  K=500 : max |resid common| 7.8e-14 Å² (EXACT)   max |resid free| 6.740 Å²
  using the ARITHMETIC mean of the free column instead of its quadratic mean:
          max |resid| 8.858 Å² (K=75), 9.017 Å² (K=500)
```

The identity is exact to machine precision with the common-frame member error and off by
−2.1 Å² — about **18% of readout²** — with the column the map prints. Anyone who reads
table A's `mean` and `diversity` together will conclude readout = 3.17 Å where it is 3.37 Å.
Both classic variants of the Sprint-16 error are reproduced here: the frame substitution
(−2.1 Å²) and, if the arithmetic rather than quadratic mean is used, a further error to
8.9 Å².

**This is not an error in a published number — the map's numbers are each correct. It is a
loaded gun in a shared instrument that three workstreams read.** Recommendation: rename the
column `mean_free`, and add `mean_common` next to it.

`readout.py` already does this correctly (§A13) — the fix exists inside the sprint.

---

## A4 — THE MAP'S OWN RANDOM CONTROL IS NOISE AT THE TARGET LEVEL

`oracle_map.py` sets `N_RAND = 3` and reports `sel_rand` as the matched random control that
BRIEF §4 makes mandatory. `audit_repro.r3_control_noise`, K = 500, all 126 targets:

```
per-target SE of sel_rand                  0.763 Å   (max 1.237 Å)
SE it contributes to the 126-target MEAN   0.071 Å
```

The **mean** column of table B is usable (0.071 Å of added SE on a 126-target average). The
**W/L counts and the `dist − random [95% CI]` column are not**: a per-target control carrying
0.76 Å of noise makes a win/loss tally on a ~0.5 Å effect close to a coin flip, and inflates
the interval. The correct control needs no draws at all — E[uniform pick from top-K] is the
mean of `rr` over the top-K in closed form (`audit_frame.exact_random`, all 126 targets, all
11 K, zero noise). Against it, the map's 3-draw estimate is biased low by 0.02–0.17 Å
depending on K, which is the sampling noise not averaging out at 60 rows.

---

## A5 — THE "K = 500 BLOSUM POOL" IS 11.4% CORPUS ORDER — AND IT DOES NOT MATTER

`s8/generate.py` builds `order = np.argsort(-sim, kind="stable")` where `sim` is an integer
BLOSUM62 sum over 9–16 residues. That gives **~57 distinct scores over ~17,000 windows**, so
the rank-500 cut always lands inside a large tie group, and `kind="stable"` resolves it by the
iteration order of `_windows_all` — peptide database first, then fold fragments, then
chain order, then sliding-window position. `audit_ordernull.py` §C, all 126 targets:

```
distinct BLOSUM scores in the universe                     mean 57  (min 42, max 76)
windows strictly above the rank-500 score                  mean 443
size of the tie group AT the rank-500 score                mean 117
pool members taken from that tie group by corpus order     mean 57  (11.4% of the pool)
```

The programme's standing memory records that reading a tie-break order as a signal
("tie-breaking leaks the pool order") once invented a 1.386 Å winner, so this needed measuring.
**It is benign here.** Averaging best-of-500 over 200 random tie-breaks per target:

```
shipped best-of-500                                1.711 Å
best-of-500 over RANDOM tie-breaks (n = 126)       1.708 Å   (mean per-target sd 0.018 Å)
```

**The 1.711 Å pool oracle is robust to the arbitrary 11.4%.** Clean survival — recorded
because a negative result on a plausible attack is worth as much as a kill.

One residual confound worth naming: the tie-break favours the peptide database, so the K=500
pool is **27.3% peptide** against **20.8%** in the universe. Since the programme's memory
records that the peptide corpus carries ~7× the sequence–structure channel of the protein
fragments, *widening K also dilutes the corpus the pool is enriched in*. Any K-effect is
confounded with that mix change, and no workstream has controlled for it.

---

## A6 — U IS CHAIN LENGTH, NOT STRUCTURAL RICHNESS

```
corr(U, chain length n) = −0.988
corr(U, best(full))     = −0.274
corr(U, best(500))      = −0.137
```

"13,000–27,000 windows per target" is a restatement of "9–16 residues": a library chain of
length m yields m − n + 1 windows, so U ≈ (constant)/n to a correlation of −0.988. The
universe is not richer for some targets than others in any structural sense.

This matters for the exhibit. **1CB3 — the four-target probe's headline, 2.028 Å → 1.096 Å —
is n = 11 with U = 26,746, i.e. the shortest chains and the largest number of draws.** It is
the target where min-of-N has the most N to work with. An exhibit chosen on the variable that
mechanically drives the effect is not an exhibit.

---

## A7 — MY OWN REDUNDANCY HYPOTHESIS FAILS, AND THE FAILURE SHARPENS A1

I pre-registered (`audit_redundancy.py`) that a "ceiling over 13,000 near-copies of 500
things" would be the same ceiling with a longer tail, and predicted the full-universe best
would be a near-duplicate of something already in the pool. **It is not.**

```
d* = RMSD of the full-universe best window to its NEAREST NEIGHBOUR in the K=500 pool
     mean 1.355 [1.272, 1.458]   median 1.363
RMSD of the full-best to the K=500 BEST window     mean 1.516   median 1.572
the ceiling gain it delivers                       mean 0.397   median 0.311
MATCHED REFERENCE: NN-to-pool distance of RANDOM out-of-pool windows   1.882 [1.872,1.895]
d* − reference   −0.527 [−0.619, −0.414]   W/L 75/51
within 1.0 Å of something already in the pool: 49/126;  within 1.5 Å: 70/126
```

The deep windows are genuine, structurally distinct conformations — the universe is a **dense
continuum**, not a set of copies. **H0(b) is REFUTED and I am recording it as such.**

That does not rescue A1; it explains it. Sampling a dense continuum more finely is precisely
the situation in which min-of-N improves for no informational reason, and the matched null in
A1 already prices that improvement. What the continuum picture adds is why the gain is
*unusable*:

```
BLOSUM rank of the full-universe best window: median 6176 of ~17,000; mean 8535
it is inside the top-500 on 18/126 targets
```

**The globally best window sits at the median of the retrieval ordering.** (18/126 is enriched
over the ~4/126 a uniform ordering would give — consistent with the small positive skill in
A1 — but it means that on 86% of targets the best structure carries no retrieval preference
at all.) A wider K makes it *present*; nothing in the sprint makes it *findable*, and A2 shows
the deployed selector goes backwards when it is.

### Distinct-conformation counts: novelty is strongly sublinear in K

Greedy leader clustering in BLOSUM rank order, exact CA-RMSD, 6 targets stratified by chain
length so U (which is 1/n, §A6) is spanned (`audit_clusters.json`; the 1.0 Å threshold was
dropped because ~9,000 leaders makes the O(U × leaders) sweep unaffordable, and that omission
is stated rather than hidden):

```
  threshold 1.5 A        K=75    K=500   K=2000   K=5000     full        U
    distinct conformations 40.3   222.3    764.2   1724.0   4249.8    18452
    C(full)/C(500) = 20.6x   while K grows 36.9x
    cluster density: 0.445 per window at K=500  ->  0.297 at full

  threshold 2.0 A        K=75    K=500   K=2000   K=5000     full        U
    distinct conformations 35.3   172.8    546.5   1127.7   2067.3    18452
    C(full)/C(500) = 12.2x   while K grows 36.9x
    cluster density: 0.346 per window at K=500  ->  0.179 at full
```

**At the 2.0 Å scale the sprint is actually targeting, the candidate set grows 36.9× and the
number of distinct conformations grows 12.2× — cluster density halves.** So the honest
statement is neither of the two available slogans: the deep universe is not near-copies
(§A7 proper), and it is not proportionally new either. It is a dense continuum whose novelty
per window falls by ~2× over the widening, which is precisely the regime in which min-of-N
keeps improving, ever more slowly, for no informational reason — and it is what produces the
decaying skill curve and the below-null marginal return of §A1.

---

## A8 — REPRODUCTION: EXACT, BY A SECOND ROUTE

`audit_repro.py` recomputes the map's core quantities without reading `oracle_map.json`.

* **R1, geometry.** `oracle_map.py` recomputes per-window RMSD as
  `I.kabsch_rmsd_batch(W, nat)` in float64; the npz already carries `rr`, computed
  independently in `s8/generate.py` by `audit.kabsch_rmsd_batch` on float32.
  Max |Δ| over **all ~2.1 M windows of all 126 targets: 4.77×10⁻⁷ Å** (worst target 1I6Y).
* **R2, the shipped anchor.** Distance-selector argmin at K = 500: **3.4540 Å**, against the
  programme's standing 3.4540. `I.selfcheck` rounds pair distances through float32 and
  `oracle_map` does not, so this is a real check, not a tautology. ORACLE best at K = 500:
  **1.7108** (standing 1.7108). ORACLE best full universe: **1.3134**.
* **Fresh interpreter** (`audit_strat.independent_headline`, separate process, re-derived from
  the npz files without importing `audit_ordernull`): max |Δ| on the observed gain **0.00**;
  the permutation-null arm differs by design (different seed stream) and lands at
  −0.072 [−0.131, −0.016] against the first run's −0.075 [−0.132, −0.022].

The coordinator's arithmetic is not the problem anywhere in this audit. As in all four
previous audits, **every number reproduced and the reading is what moved.**

---

## A9 — LEAKAGE AND SEEDING: CLEAN

* **Sealed benchmark.** A walk of every `.py` under `s17/` for `benchmark60|bench60|
  benchmark_60` and for `bench_results` returns **only this audit's own docstrings and its
  own search regex**. No sibling module references the sealed set, directly or by path.
* **Ordering leakage.** `u["order"]` was verified to equal `np.argsort(-sim, kind="stable")`
  on every target, and `sim = B62[S, encode(seq)] .sum(1)` is a function of the **window
  sequence and the target sequence only**. No native quantity enters the ordering. The one
  non-sequence input to the ordering is the corpus iteration order (§A5), which is
  native-free but not sequence-derived.
* **Seeding.** No bare `hash()` in `s17/`. All randomness in this audit's own modules goes
  through `s15.seed.stable_rng`.
* **Fold isolation.** The universes are built per target from out-of-fold peptides plus that
  fold's fragments (`s8/generate.py`), and identity-thresholded; nothing in this audit widened
  that.

---

## A10 — `oracle_map.json` IS A SHARED FILE WITH NO COMPLETION FLAG

`s12/instrument.write()` carries a documented HAZARD note — it has no config key, so a
partial run silently overwrites a complete one — and takes `n_expected` to stamp an explicit
`complete` flag. `oracle_map.py` does not use it: it checkpoints with a raw
`json.dump({"rows": rows}, ...)` every 20 targets to the same path it writes the final result
to.

While this audit was running the file passed through **60 of 126 rows**, and because targets
are processed in pinned pdb order the partial set was **fold-imbalanced**: `{fold 4: 16,
fold 2: 13, fold 0: 13, fold 1: 12, fold 3: 6}` — fold 3 under-represented 2.7×. Three sibling
workstreams read this file as "the shared instrument of this sprint" during that window, and
**nothing in it said it was incomplete**: any fold-clustered interval computed from it then is
wrong, and §A2b's own first pass (n = 60) is an instance — its `full` interval was
[−0.116, +0.456] and at n = 126 it is [+0.046, +0.253]. The run has since completed (126/126),
but the file still carries no flag, so the hazard recurs on the next rerun. One line
(`I.write("oracle_map", {...}, n_expected=126)`) fixes it.

---

## A11 — POWER OF THE INSTRUMENTS IN USE

BRIEF §5 records five Sprint-16 nulls quoted on an instrument that could not have detected the
effects being denied. For the 126-target universe instrument (`audit_repro.r5_power`):

| statistic | sd | SE | MDE at 80% power, α = 0.05 |
|---|---|---|---|
| paired order-null difference | 0.335 | 0.030 | **0.084 Å** |
| the K-widening ceiling gain | 0.359 | 0.032 | **0.090 Å** |

Contrast the 9-target enumerated instrument (SE 0.051 Å, MDE ≈ 0.16 Å). **A null declared on
the 126-target universe instrument is meaningful down to ≈0.09 Å and no further.** The A1
effect (−0.075 Å) sits just below that MDE and is nevertheless significant because the CI is
computed on the effect itself rather than by a power argument — but any workstream declaring
"no effect" on this instrument must quote 0.09 Å, not zero.

---

## A12 — `selection_theory.rho_from_accuracy` SURVIVES THE FUNCTIONAL CHECK

The module maps a pairwise ordering accuracy to a Gaussian-copula ρ by
`ρ = sin(π(2·acc − 1)/2)`, i.e. τ = 2·acc − 1 followed by the Gaussian-copula identity
τ = (2/π)·arcsin(ρ). Both steps are correct, and the module's own docstring records that the
coordinator initially read the recorded 0.600 as ρ = 0.600 when it is ρ = 0.309. I checked the
**axis** rather than the algebra, which is the check that matters here: the programme's 0.600 /
0.638 / 0.986 figures originate in `s14/obj_FINDINGS.md` and `s14/coord_FINDINGS.md`, where they
are explicitly **in-band ordering accuracy with null 0.5**, not correlations. The conversion
therefore consumes the right functional.

**One caveat the module should carry:** `s14/figures/README.md` records that because the pairs
are matched on sequence separation, the null for that accuracy is **0.505, not 0.500**. τ =
2·acc − 1 assumes 0.500, so every ρ derived this way is high by ≈0.016. Immaterial at the
precision reported; it should still be stated, because the sprint is about to price an
architecture on the difference between ρ = 0.309 and a requirement.

---

## A13 — `readout.py` DOES CORRECTLY WHAT `oracle_map.py` TABLE A DOES NOT

`s17/results/readout.json` reports per (target, K, m) an `err`, a `div`, an `rmsd` and a `cf`.
Checked on its own artefact: for 1A13, K = 500, m = 20, `sqrt(err² − div²) =
sqrt(2.68717² − 0.90588²) = 2.52988` against the recorded `rmsd = 2.529876` — the identity
closes to 10⁻¹⁵, and `cf` equals `rmsd` at every cell inspected. Its `err` is therefore the
**common-frame** member error, which is the functional the identity consumes. This is the
correct treatment and it is already inside the sprint; A3's fix is to copy it into the map.

Two smaller notes on the SELECT instrument (`sel_bench.py`), neither an error:

* `random` and `mean` are numerically identical by construction (the exact expectation of a
  uniform pick is the set mean). They are two names for one column and should not be counted
  as two pieces of evidence.
* `bestof_m` at m = k is the oracle by construction (`bestof75 = oracle = 2.243894` at K = 75).
  Any arm compared against `bestof75` at K = 75 is being compared against the ceiling, not
  against a control.

---

## A14 — TWO DIFFERENT "TOP-75" ARE PRINTED WITH THE SAME NAME, AND ONE MECHANISM IS
## ATTRIBUTED TO THE OTHER

`oracle_map.report()` §C prints

```
  ORACLE best in top-75      [truncation ceiling]    2.104
  ...
  the top-75 truncation costs  0.791 A of ceiling
```

and `BRIEF.md` §7's standing table carries **"ORACLE best in shipped top-75 | 2.306"**. A
reader will match 2.104 against 2.306 and conclude the map has improved the top-75 ceiling by
0.2 Å. **These are different sets.** Recomputed over all 126 targets:

```
BLOSUM top-75 (the map's k75)  oracle best   2.1041
shipped top-75 (rec["sub"])    oracle best   2.3062
difference −0.2020,  W/L 62/53
mean overlap of the two 75-member sets:      13.2 / 75
```

The two shortlists share **13 of 75 members**. The map's `k75` is the first 75 entries of the
BLOSUM ranking; the deployed top-75 is the distance objective's shortlist drawn from the whole
K = 500 pool. They are produced by different mechanisms and neither is a truncation of the
other.

The consequence is the mis-attribution, not the number. **"The top-75 truncation costs
0.791 Å of ceiling"** attributes to *truncation* a loss the deployed system incurs by
*selection*. The deployed shortlisting costs 2.3062 − 1.7108 = **0.595 Å relative to its own
input pool**, and it is a ranking loss, not a retrieval loss — which is the distinction BRIEF
§3 makes mandatory (problem **A** vs problem **B**). Recommendation: rename the map's row
`ORACLE best in BLOSUM top-75` and keep `shipped top-75` as a separate row.

---

## A15 — WHAT SURVIVED EVERYTHING

A clean survival is as useful as a kill, so these are stated with the same weight.

1. **Every published number reproduced, by a second route.** 3.4540 / 1.7108 / 1.3134 /
   2.3062, and two independent RMSD implementations at two dtypes agree to 4.8×10⁻⁷ Å over
   ~2.1 M windows. This is the fifth consecutive audit in this programme to find no
   transcription error, and the fifth to move conclusions anyway.
2. **`u["order"]` is native-free.** Verified equal to `np.argsort(-sim, kind="stable")` on
   all 126 targets, with `sim` a pure function of the window sequence and the target
   sequence. No native, no RMSD, no structural quantity enters the ordering.
3. **The 1.711 Å pool oracle is robust to its own arbitrariness.** 11.4% of the pool is
   corpus-order tie-break, and randomising it moves the number by 0.003 Å (§A5).
4. **No leakage of the sealed benchmark**, no bare `hash()`, in any `s17/` module (§A9).
5. **`selection_theory.rho_from_accuracy` consumes the right functional** (§A12) — the one
   place in the sprint where a law's input axis was checked before use, and it holds.
6. **`readout.py` measures the member error in the frame the identity consumes** (§A13).
7. **`sel_bench.py` carries a matched random control for every gate and every consensus
   arm**, which is what BRIEF §4 requires and what Sprint 16 lacked.
8. **The BLOSUM ordering does have retrieval skill** — my own H0 that it was structurally
   uninformative is falsified: skill is positive at every K with CIs excluding zero (§A1).
   It is just small (+1.5 targets of 2.0 Å recall over a random 500) and decaying.
9. **The deep universe is real structure, not near-copies** — my own redundancy hypothesis
   failed (§A7).

## A16 — WHAT THIS AUDIT DID NOT SETTLE

* **The distinct-conformation counts are n = 6 targets.** The 1.0 Å threshold was not run at
  all. The sublinearity in §A7 is a direction, not an interval.
* **The frame check is n = 30 targets.** The identity residual is exact/inexact by
  construction, so n does not affect the qualitative finding, but the −2.10 Å² magnitude is
  a 30-target mean.
* **§A2b's interval disagrees between the fold-clustered and plain target bootstraps** at the
  boundary. I have not established which is right for this statistic; the direction is
  consistent at all seven K and that is what I am relying on.
* **I did not audit the physics, quantum or repair workstreams' outputs**, which landed
  during this run (`phys_ident_leg.json`, `quantum_pareto.json`, `q_theory.py`,
  `sel_hard.py`, `sel_inband.py`). None of their claims are covered by anything above.
* **The peptide/protein corpus-mix confound (§A5) is named and not measured.** Widening K
  takes the pool from 27.3% peptide to 20.8% peptide, and the programme's own memory records
  that the peptide corpus carries ~7× the sequence–structure channel. Every K-effect in this
  sprint is confounded with that, and nobody has controlled it. **This is the single most
  valuable unrun experiment I can name**: rerun §A1's ceiling curve within each corpus
  separately.

## A17 — RECOMMENDATIONS, in the order I would act on them

1. **Fix `oracle_map` §C's two mislabelled rows** (§A2, §A14) before any workstream quotes
   the budget table. `raw coordinate average` must say which set, and `top-75` must say
   BLOSUM or shipped. As printed, the sprint's central architectural comparison has the
   wrong sign.
2. **Rename table A's `mean` to `mean_free` and add `mean_common`** (§A3), or delete the
   column. It is the Sprint-16 error sitting in a shared instrument.
3. **Replace `sel_rand` with the closed-form expectation** (§A4). It is one line, exact, and
   removes 0.76 Å of per-target control noise.
4. **Stop reading BRIEF §2 as "truncation is costing 0.397 Å"** (§A1). The correct sentence
   is: *the K = 500 pool already sits near the efficient frontier of this library's sampling;
   the extra ceiling at large K is the order statistic of more draws, it under-delivers
   against a random-order null, and on the 53 hard targets and all 18 FAIL18 targets it
   delivers exactly the null amount and nothing more.* Widening K is not a route to 2.0 Å.
5. **Add `n_expected` to the map's writer** (§A10).
6. **Run the corpus-mix control** (§A16).

## A18 — PRE-REGISTRATION RULES: DID THEY FIRE?

| module | rule | fired? |
|---|---|---|
| `audit_ordernull` | H0: BLOSUM order carries no structural skill | **NO — falsified.** Skill positive at every K. Reported as such. |
| `audit_ordernull` | success criterion for H1 (coordinator's reading): observed gain must exceed the matched null with a CI excluding zero | **NO.** It is *below* the null, CI [−0.132, −0.022]. H1 rejected. |
| `audit_redundancy` | H0(b): the full-universe best is a near-duplicate of a pool member | **NO — falsified.** d* = 1.355 Å. Reported as a failure of my hypothesis. |
| `audit_redundancy` | H0(a): distinct conformations grow far more slowly than K | **YES.** 12.2× against 36.9× at the 2.0 Å scale. |
| `audit_frame` | the identity is exact in the common frame and violated by the printed column | **YES.** 7×10⁻¹⁴ vs −2.10 Å². |
| `audit_strat` | concentration check against a uniform-effect null | **NO concentration** — observed drop-top-10 inside the null band. |

## A19 — ERRORS THIS AUDIT MADE, PRESERVED IN PLACE

**`audit_strat.py` §D, first version.** I printed the recall of the *random-order full
universe* as a control for the *BLOSUM full universe* and it read 77/77, 102/102, 118/118 —
"identical, therefore the recall gain is entirely an order-statistic effect". That is a
tautology: at K = full the candidate set is the whole universe under **every** permutation, so
the two columns are equal by construction and the "control" carries no information. It is
precisely the shape of mistake this audit exists to catch — a quantity computed correctly and
read as a different quantity — and I made it while writing the section that catches it.
The comment block in `audit_strat.py` preserves the error, and the replacement
(`recall_ctrl`, a **matched-K** control: BLOSUM-500 vs random-500) is what §A1 now quotes.
It changes the finding: the honest number is +1.5 targets of retrieval benefit, not zero.

---

*Modules: `s17/audit_ordernull.py`, `s17/audit_redundancy.py`, `s17/audit_repro.py`,
`s17/audit_frame.py`, `s17/audit_strat.py`. All randomness through `s15.seed.stable_rng`.
No module reads benchmark60.*
