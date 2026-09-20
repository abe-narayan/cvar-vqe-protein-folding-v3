# S30 — RUNNING STATE

Last update: sprint open.
Endpoint: mean built-chain Cα RMSD, 126 targets. Production incumbent **3.2105 Å**.

---

## HEADLINE: THE MEAN IS A TAIL STATISTIC, AND THE TAIL IS WHERE THE PRIZE IS

Computed at the sprint open from lane O's S29 rows, before any new experiment:

```
production, n = 126        mean 3.2105    median 2.9661
the worst 18 targets       mean 6.2758
the other 108              mean 2.6997
```

**Counterfactuals on the charter's own endpoint:**

```
cap the worst 10 at 3.00 Å   ->  mean 2.9074   (-0.3031)   BEATS the primary target
cap the worst 18 at 3.00 Å   ->  mean 2.7426   (-0.4680)
cap the worst 30 at 3.00 Å   ->  mean 2.5778   (-0.6328)   approaches the ambitious target
worst 18 all the way to 2.50 ->  mean 2.1334   (-1.0772)

versus: improve EVERY ONE of the 126 by 0.20 Å  ->  mean 3.0105  (-0.2000)
```

> **Fixing ten targets is worth more than improving all 126 by 0.20 Å.**

That single comparison should govern how this sprint allocates effort. A mechanism that works on
the typical target and leaves the tail alone is close to worthless *for this endpoint*; a mechanism
that only works on hard targets can hit the primary target on its own.

**The caveat, stated up front so the number is not over-read.** Capping is an ORACLE operation — it
requires knowing which targets are bad. It bounds the prize; it does not deliver it. Realising any
of it needs either (a) a native-free regime detector, or (b) a method that is simply better on hard
pools without needing to know they are hard. **(b) does not require detection and is therefore the
stronger target.**

**And there is a live, unexploited clue.** The record says the shipped pipeline is **worse than a
sequence-blind one on its 18 hardest targets** — 5.425 blind against 6.019 shipped — while
sequence conditioning is worth 0.776 Å overall. Something the pipeline does on hard targets is
actively harmful. That is a mechanism to find, not a curiosity.

---


## NOTE 1 (2026-09-20 12:40, lane R, before its numbers): A FINDING I HANDED OVER AS SETTLED IS CONFOUNDED

I briefed lane R with S28-L48 as established: *20 of 31 native-free scorers prefer the contracted
production average to a 0.25 Å ORACLE structure.* It is the most-cited negative in the project's
recent record — it is in the S29 report, in the published artifact, and it is why I called
recognition "the crux."

**Lane R's objection, which I accept: every rung in that ladder differs in KIND as well as in
nativeness.** Production is a contracted average; the ORACLE rungs are circuit outputs; the
controls are perturbations. So the measurement is a **cross-kind preference**, and the
perception–distortion theorem (S29-L12) already predicts exactly that without any claim about
whether nativeness is visible. **S28-L48 therefore establishes something weaker than I said: that
these scorers prefer one kind of object to another, not that nativeness is unrecognisable.**

**What this does and does not touch.** It does *not* undo S29's "recognition closed three ways" —
that rests on different evidence (the across-band identity, F2 failing on 58/70 *within* band, and
the incidental-parameter theorem). It removes one pillar from the specific claim that *single-
structure geometry cannot see nativeness*, which is precisely the question lane R is now testing on
an instrument where kind, local geometry and perturbation budget are matched and only nativeness
varies.

**The instrument's one-line case**, from lane R's own probe on 1A13: resample **one** residue's
torsions from the fold's leakage-safe Ramachandran table and RMSD-to-native spans **0.284 to
4.685 Å**. Same kind, same local statistics, one torsion, 4.4 Å of nativeness. If nothing can order
that, nothing can see nativeness.

**A derivation I have asked lane R to make explicit, because it may predict the answer.** If a
single torsion swings global RMSD by 4.4 Å while every local statistic is matched by construction,
then any purely *local* channel is blind to that variation **by construction**, and its ceiling is
set by how much RMSD variance is non-local. If that fraction is near zero for local channels, the
negative result becomes a theorem on this instrument rather than an empirical miss — and the live
question narrows to channels with genuinely global reach.


## NOTE 2 (2026-09-20 12:43, lane L, derived not cited): THE COMPACTNESS LOADING IS **ALGEBRAICALLY FORCED**, AND THERE IS A CONSTRUCTIVE FIX

A statistical potential is `u(d) = -kT ln[P_obs(d)/P_ref(d)]`, so its total score contains a
separable term `+kT * sum_pairs ln P_ref(d_ij)`. For DOPE's reference state -- non-interacting
points in a ball of radius `a = sqrt(5/3)*Rg` -- `P_ref(d;a) = (1/a) g(d/a)` **exactly**. That term
is a pure function of the structure's **scale**, forced by the construction rather than fitted.

**The decisive number.** Contract a 13-mer uniformly by 10% -- pure scale, **zero shape change** --
and a fixed-reference score moves **-0.256 kT per pair**. Over ~55 scored Ca pairs that is
**~14 kT of reward for being smaller with no shape content at all.** A size-matched reference gives
**exactly 0.000**, by scale invariance rather than by fitting.

**A second effect specific to our length.** At n = 13 the reference state's entire support is
[0, 15.05 A] while DOPE is tabulated to 15 A; at n = 9 the support is 13.09 A, so the top 13% of the
tabulated range has **zero reference density**. This is the mechanism behind the authors' own
"less accurate for smaller proteins", and behind the 40-50 residue wall S29-L1 found empirically.
We are not approaching that wall -- we are inside it.

**WHAT THIS CHANGES.** S29-L50 reported the compactness loadings (LEG_compactness +0.956,
RG_LAW +0.921, LEG_solvation +0.621, LEG +0.606) as a measured pattern. They are **algebraically
forced**. That is a stronger statement and a correction to how S29 framed its own result.

**AND IT GIVES A CONSTRUCTIVE ALTERNATIVE TO PARTIALLING.** Rebuild a pair channel with a
**per-candidate** reference `a_i = sqrt(5/3)*Rg_i`: provably scale-invariant, compactness-free
**by construction**, no partialling needed. Nothing in `s27/ham_lib.py` does this -- every reference
there is separation-only or quasi-chemical, which absorb chain connectivity and composition but
**not size**. This matters for a verdict: a channel that survives statistical partialling has had a
correlate removed post hoc and a sceptic can argue about what else went with it; a channel that is
scale-invariant by construction leaves nothing to argue about. Relayed to lane R mid-build.

**THE CONSTRAINT LANE L STATED UNPROMPTED, AND I AM HOLDING IT TO:** this creates **no new
information**. It is a re-parameterisation of an existing structure->score map, so by the bound it
reaches the endpoint only as a cosine, and the averaging readout spends at most 0.04 of any
ranking. It is a **correctness fix and a band statistic**, not an Angstrom route. The published
trade is explicit and runs against us in one direction -- ANDIS: native recognition and decoy
discrimination cannot be optimised simultaneously, which *is* this trade.

**PENDING FROM LANE L, AND I RATE IT HIGHER THAN THE ABOVE:** an exact non-identifiability proof
that the pool's common-mode error is invisible **from within the pool**. If it lands it converts
"68% common-mode" from a measured pathology into a structural impossibility for an operator class,
makes lane O's *exactly zero* PC1 result a corollary, and either justifies lane X's whole premise
or shows the escape must come from outside the pool. I have asked for the quantifier to be exact:
invisible *from within the pool* is very different from invisible *to any method*, and the
difference is where the sprint's remaining hope lives.


## NOTE 3 (2026-09-20 12:46, lane F, S30-L2): **THE TAIL IS SELECTION-LIMITED, AND THE DAMAGE IS IN ONE STAGE**

The sprint's first real mechanism, and it redirects the lane.

**The tail is NOT pool-limited.** The ORACLE best member of the FAIL18 pools is **2.2842 Å**
(chain 2.2845) -- already below the 3.00 Å cap the opening counterfactual asks for. **13 of the 18
have an ORACLE pool member under 3.00 Å**, and the worst tail pool bottoms out at 3.54. The
material to cap the tail is already inside the candidate sets we hand the pipeline. The
retrieval/generation branch of lane F is closed on its own evidence.

**The gap decomposes, and one stage owns the tail** (ORACLE, point cloud):

```
                          FAIL18    other 108
G_retr  pool - universe    0.644      0.356
G_filt  top75 - pool       2.393      0.296     <- 67.5% of the tail gap vs 30.5%
G_read  prod - top75       1.155      0.673
G_tot                      3.548      0.969
```

**Retrieval is exonerated at every stratum.** BLOSUM's 500 against a random 500 of the same
universe is NOT MEASURED everywhere, and most pointedly on FAIL18 (0.01× its own MDE). This
**revises the mechanism recorded in project memory** as `sequence-conditioning-hurts-the-failures`:
the harm is the distogram **filter**, not the retrieval corpus.

**The filter is worse than chance on every single tail target.** Score's top-75 against a random 75
of the same pool: other-108 NOT MEASURED (−0.027, 0.24×); **FAIL18 +1.7674 at 3.87× MDE, 4/4 folds,
0W/18L**, median per-target percentile **0.99999** -- on the median tail target essentially every
random 75-subset of the same pool contains a better member than the score's 75 does. Random-18
null, 20,000 draws: **p = 0**.

**Why it reaches the endpoint.** The terminal operator consumes the set mean, so best-member is not
an endpoint story by itself. On the set mean the filter buys **−1.0855 on the 108** (5.11× MDE,
5/5 folds, 100W/8L) and that benefit **vanishes on the tail** (+0.1947, 0.39× -- NOT MEASURED, and
lane F does not claim degradation). Through the operator's 1.16 coefficient a vanished −1.09 is
~+1.26 Å of emitted RMSD -- the right order for the record's blind-beats-shipped gap on FAIL18
(5.425 vs 6.019).

**CIRCULARITY, declared and calibrated by lane F before anyone asked.** FAIL18 is defined by
production, which is the filtered set's average, so a bad filter mechanically lands a target in
FAIL18. Two filter-independent tails through the same null: worst-18 by pool mean **+0.6708**
(p = 0.0148), worst-18 by ORACLE best-in-pool **+0.6320** (p = 0.0300). Both clear; FAIL18's +1.77
is inflated **~2.7×** by its own definition and the smaller figures are quoted beside it.
**The load-bearing, non-circular corroboration:** 44/126 targets have a filter worse than 90% of
random 75-subsets; **26 of them are outside FAIL18**, and their production mean is 3.1360 against
2.4095 for the other 82.

**WHAT I HAVE ASKED FOR NEXT.** (1) Price the deployable no-filter/soften-filter family on the
**endpoint** -- the only route-(b) member on the table, expected to lose net because the 108's
−1.09 cannot be given up. (2) **The variable is filter WIDTH, not filter on/off** -- if the 108's
benefit saturates early in k while the tail harm grows with aggressiveness, a wider filter keeps
the benefit and sheds the harm **with no detection at all**. Production's top-75 is *both* filter
and averaging set, so these must be separated: filter to k, average over m ≤ k. S29's flat m-ladder
(−0.00047 Å per unit) is a property of the **averaging** width and says nothing about the filter
width. (3) Then open the score: *what is the distogram confidently wrong about on those pools* --
not a ninth router, which the charter closes explicitly.

**Two joins to lane L.** If hard pools differ in Rg dispersion, a scale-loaded score behaves
differently on them for **algebraic** reasons (note 2). And if common-mode error is unidentifiable
from within the pool, a score computed over the pool is structurally blind to the shared component
-- which would say why confident wrongness concentrates where it does.


## NOTE 4 (2026-09-20 12:49, lane T): TWO THEOREMS. THE TAIL IS **ALWAYS** A PREFIX, AND **THE POOL IS THE CODEBOOK**

### (b) The tail condition — the question was mis-posed, mine as much as S29's

**T1, the endogenous-order prefix theorem.** For any tail objective V differentiable on the
box-simplex {0 ≤ λ ≤ p, Σλ = α}, every KKT point has λ*_x = p_x where ∇V(λ*)_x < μ and 0 where
> μ. **The tail is always a prefix — of the order induced by ∇V at the optimum.** Prefix-hood is
universal. What varies is whether that order is **exogenous** (known before solving) or
**endogenous** (a fixed point). The real condition: *the tail stops being reproducible by one
classical sort iff the ordering map λ → ∇V(λ) has more than one fixed point.*

**A correction to S29 that the sprint needed.** S29 §4.3 recommends fixing the tail's order by a
per-state scalar for well-posedness; §4.5 claims set-equality then fails. **Those are
incompatible.** With an exogenous order the emitted set is still `argsort(E)[:m]` and the endpoint
channel is still the single integer *m*. S29's counterexample proves a **free** subset optimum is
non-prefix; it does **not** show the lifted CVaR readout can reach it. **The tail-then-aggregate
lift as specified does not remove the bottleneck** — which is also why S29's endpoint arm came back
refuted, and now we have the reason and not only the null.

**T1b, the reachability cap — the pre-check fired before any box was spent.** With an endogenous
order the cost is `c_x = ⟨∇f(R_λ), W_x⟩`: candidates sorted by projection onto **one**
self-consistent direction. Reachable tails are **halfspace cuts**, VC dim d+1. Measured on all 126
real pools: **stable rank of the centred pair-distance matrix = 1.859** (median 1.865, max 2.72),
PC1 = 55.4% of variance, k90 = 5.6. The class carries ~31 effective bits against
**log₂ C(500,75) = 300.6** of free set choice — **a 269-bit collapse.** Lane T's pre-registered
rule (stable rank < 2.0) fires and **lane Q's direction is closed at the encoding level.** Same
failure mode that killed S29's first non-diagonal Hamiltonian, caught this time *before* the spend.

**THE ONE ESCAPE THAT SURVIVES.** Add a **second-moment (dispersion)** term,
`V = f(Σλ W, Σλ W Wᵀ)`. Then ∇V_x is **quadratic** in W_x, the cut is a **quadric not a
halfspace**, and VC dim jumps from 34 to **≈595**. It is the *same object* S29's own post-mortem
pointed at — "operators that read the pool's own dispersion rather than the posterior's marginals"
— reached by an independent derivation. **Two routes converging on one class is the strongest
signal this sprint has produced about where to build.** And it gives CVaR a precise residual role:
**α supplies the threshold that turns a direction into a set.**

### (a) The bit accounting — the charter's question was not well-posed either

```
stage 5 retrieval  3,252 bits of CHOICE  ->  DELIVERED +0.69 bits   (Å contrast NOT MEASURED, 0.79×)
stage 9 readout        7 bits of CHOICE  ->  DELIVERED +0.036       (median -0.575)
end to end: the deployed system's measured information yield is UNDER ONE BIT PER TARGET
```

**Why the 5.56 bits are not missing: THE POOL IS THE CODEBOOK.** The 500 deposited backbones carry
the structure; the index only names it. Bits are not conserved across an index. The cross-check
settles it: 7 index bits move 4.108 → 1.898 Å, which through the displacement bound is ρ = 0.887,
i.e. **36.6 bits of displacement information out of 7 index bits — a 5.2× ratio.** So the honest
statement is the reverse of the charter's premise: **the readout's 7 bits are worth five times
their face value and the system cannot supply even one of them.**

**THE VALUE-OF-A-BIT LAW** (ORACLE ladder in the deployed pool: 4.1080 → 1.8978 → 1.7108, fitted
`D(R) = a + c·2^(−R/γ)`, a = 1.3312 Å, γ = 3.1636, **R² = 0.9983**):

> **−dD/dR = (ln2/γ)·(D − a) = 0.219·(D − 1.331) Å per bit.** Marginal at R = 7: **0.132 Å/bit.**

**ALLOCATION, ranked — and the encoding is NOT the bottleneck.** Candidate identity 0.132 Å/bit
(floor 1.331, γ 3.16); subset cardinality 0.044 — **dominated 3.0×, which explains the 3.5× S29
measured**; mode/basin saturates at ~1.6 bits; **torsion/configuration space is arithmetically
infeasible at the deployed width** — 2n·log₂(k) = **48 bits** for 4 basins per residue at n = 12,
against 7 deployed or 9 in the harness, i.e. 0.29 bits per torsion where 1 bit names a single
Ramachandran basin. **Binary candidate indexing is optimal among the measured classes by 3×**, and
capped at 1.331 Å. The charter called the encoding the least-examined component and the likely
hidden bottleneck; lane T examined it and **it is not**.

**THE DICTIONARY THAT PUTS THE WHOLE SPRINT IN ONE CURRENCY** (d = 3n−6 = 32.9):

```
rho <= 0.04  (every field ever built)  ->   0.038 bits
rho  = 0.140 (B2's ceiling)            ->   0.470 bits
rho  = 0.358 (3.00 A)                  ->   3.25  bits
rho  = 0.628 (2.50 A, the charter)     ->  11.9   bits
```

**The charter's target needs ~12 bits of displacement information per target and the entire
native-free operator space supplies 1/26th of one bit.**

**Lane T's own registered prediction failed in the flattering direction and it said so**: +0.10 to
+0.30 Å predicted for retrieval, **0.07 measured**, copula reasoning overestimating 2–4×. Second
time this sprint a lane has reported its own prior failing that way.

**A caution I have given lane T about its own framing.** "The encoding is not the bottleneck" is
established *among the measured classes*. Its own dispersion escape is a statement about the
**operator**, and lane F has just located a large loss in a **stage**. Three lanes now point at
three different components and the report must keep them distinct.


## NOTE 5 (2026-09-20 12:53, lane D, S30-L5): **THE FIELDS ARE NOT NOISE.** B2's NUMBER STANDS; ITS ARGUMENT IS BACKWARDS

The most consequential correction of the sprint, and it **inverts the strategy** rather than
adjusting a number.

**The defect, found by reading the code and not the prose** (`s29/s29_D_fields.py:243,263`).
Two errors compound:
1. The verdict string says "no field's **mean |cos|** clears the random reference 0.140" but the
   code compares `abs(signed mean)`. The quantity the sentence names runs **0.2505-0.3249 and
   clears 0.1398 on 21 of 21 fields.**
2. The null is wrong for either quantity. `ref = 0.1398` is the magnitude of **one random direction
   on one target**; `cmp0["effect"]` is a **mean over 126 targets**. They differ in scale by ~sqrt(126).

**Against the right nulls** (from S29's own stored draws): null for the 126-target signed mean is
**+0.0014 +- 0.0144**; for the mean |cos| it is 0.1398 +- 0.0099. Result: **11 of 21 fields have
fold-clustered CIs excluding zero, best +7.7 sigma** (CHAN_DISTPOT +0.1128), where multiplicity
predicts one. S29 recorded `any_beats_reference = []`; the correct count is **11**.

**The objection answered from inside the survey.** A Gaussian random direction might be a *worse*
measure rather than an uninformative one (`zero-information-control-must-be-plausible`). The survey
answers itself: three of its own **structured** fields -- MEDOID -0.0097, MSET_50 -0.0187,
EXPAND -0.0208 -- sit *at* the signed null. So +0.11 is not a generic property of structured
directions.

> **B2's number stands (best exploitable rho = 0.1128 <= 0.14). Its argument does not. The fields
> are real and they are worth 0.0195 A because sqrt(1-rho^2) SQUARES them.**

### WHY THIS CHANGES THE SPRINT

Old reading: *nothing carries signal, there is nothing to build on.*
New reading: **eleven directions carry highly significant signal, and each is individually
worthless because the transfer function is quadratic near zero.**

Those imply different next moves, and the second has a question nobody here has asked:

> **What rho does the best COMBINATION of the 21 fields achieve?**

If they were orthogonal, a combination reaches rho = sqrt(sum rho_i^2); eleven fields at ~0.10-0.11
gives **~0.33-0.36**, against the **0.358 the bound requires for 3.00 A**. They are almost
certainly not orthogonal -- all built from the same distogram and pool, and project memory
(`decorrelated-errors-exist-but-are-unusable`) records fusion gains going as the **square** of the
weaker channel's skill. But the arithmetic lands close enough to the threshold that it must be
**measured, not assumed**, and the per-target displacement vectors already exist.

**Requested of lane D, ahead of its queue:** (1) the **Gram matrix** of the 21 directions and its
effective rank -- if it collapses to ~2 like the pair-distance matrix (stable rank 1.859) the
combination buys nothing and this closes in an afternoon; (2) the **ORACLE-optimal combination's
rho**, the ceiling of the whole field class, labelled ORACLE; (3) only if that ceiling clears ~0.3,
whether a **native-free, global, leave-fold-out** weighting recovers a useful fraction -- per-target
weights are not estimable (incidental parameter).

In lane T's currency: rho 0.14 = 0.470 bits, rho 0.358 = 3.25 bits. **The question is whether
eleven 0.04-bit channels can be combined into 3.25 bits, and the answer is governed by their RANK,
not their count.**

**The meta-lesson, one level up from the four dangling artefact paths:** this error survived a full
sprint, a published report and a published artifact because everyone downstream -- me included --
quoted the **verdict string**. A verdict string is a claim about a computation and is worth exactly
what reading the computation is worth.


## NOTE 18 (2026-09-20 14:12, lane V, `s30/AUDIT_V.md`): **THE REPORT ADVERSARY FOUND TEN DEFECTS IN MY REPORT, TWO OF WHICH INVERTED A HEADLINE**

Contract rule 28 earned its place. Lane V audited §1, §3, §4 and Appendix A against the artefacts
and returned **NOT SAFE TO PUBLISH**. All ten are now fixed in place with the original error stated.

**The two that inverted a headline:**

1. **D1 — §A.7 argued "the effect does not reach the tail" on `FAIL18` and dropped the
   filter-independent stratum sitting beside it in the same artefact**, where both legs reverse
   (relax gain **−0.0386 on the tail against −0.0193**, not −0.0113 against −0.0239; tail dispersion
   **+1.6116 at 1.44× MDE with the CI excluding zero**, not +0.4583 at 0.47×). And the FAIL18
   difference I called "the wrong direction" is **0.22× MDE** — *I used a number below my own
   not-a-result line as evidence of absence.* E2 now stays out on **size**, not on reach.
2. **D2 — the "36.6 bits / 5.2× / worth five times their face value" headline is a free parameter
   this sprint corrected downward.** It is `I = −(d/2)log₂(1−ρ²)` at `d = 3n−6 = 32.88`, published
   with **neither the formula nor the `d`**. At the corrected `d = 6` it is 6.69 bits and **0.95× —
   the multiplier inverts.** Lane T's own S30-L14 §6 says the corrected numbers are smaller.

**The rest:** "the 18 worst pools" was FAIL18, which is not the 18 worst (overlap 13/18; real
figures 2.5298 / 11 of 18 / 3.9523); "replicates on all three tail definitions" was false
(0.1066 vs 0.3798 vs 0.3438, and only the filter-defined one has a CI spanning zero); the bolded
"anti-aligned on the hard targets" is **0.67× MDE**, withdrawn; "retrieval is exonerated" is a
**control-space mismatch** — which falsified my own §6.1 row claiming no such instance this sprint.

> **My two priority recomputations both survived** — 2.2842 Å recomputed from raw rather than from
> lane F's JSON, and ΔR² −0.08942 / +0.59971 at 2.53× and 25.52× MDE. **The spine held; every
> defect removed a multiplier, a stratum label or an over-claim, and none removed a result.**

**Two things lane V corrected that were mine and not in the report:** the "cross-lane reproduction"
at ρ = 1.0000 is **exact by construction** (`s30_G_disp.py:189` reads lane F's variable straight in;
`:187` recomputes it with `ddof=1`, differing by √(75/74)) — a genuine membership check, *not* an
independent reproduction, and I had called it the sprint's only one. And the verifier's "22/22" was
narrower than it sounded: **all 22 checks came from lane D**, with every other lane `show()`-only.
It is now **36 checks** including the endpoint itself.

---

## NOTE 17 (2026-09-20 14:10, lane W, `s30/QUANTUM_W.md`): **"NO VQE" AND "NO QUANTUM COMPUTE" ARE DIFFERENT CLAIMS, AND ONLY THE FIRST IS TRUE**

My prior — that S30 ran no VQE — is right in substance and wrong in one checkable way that changes
how the report must be worded.

**S30 executed a 9-qubit depth-3 statevector circuit on all 126 targets**, exact parameter-shift
Jacobians, 300 Adam iterations each, to regenerate the meter's ORACLE rungs
(`s30/results/s30_D_ladder_structs/`, 126 `.npz`, asserted to 1e-6 against S28). It is **not** a VQE
— the objective is ORACLE point-cloud RMSD — but *"no quantum compute was spent this sprint"* would
have been false as written, and I would have written it.

**The row that answers charter items 12 and 16 at once, from this sprint's own cache:**

```
circ_best   same circuit, ORACLE objective, best of 5     0.2516 chain / 0.2884 CA
circ_s0     ORACLE objective, 1 start, regenerated S30    0.3175 / 0.3854
PROD        the deployed uniform average                  3.2071 / 3.0483
circ_opt    SAME CIRCUIT, DEPLOYED native-free score      3.4330 / 3.3850
```

> **The same circuit reaches 0.2516 Å with the native as objective and 3.4330 Å with the shipped
> one — worse than the classical average it was meant to improve. The circuit is not the problem.**

**And a defect in shipped code:** `core/pipeline.py:821` still asserts *"the CVaR tail is worth
+0.113 Å … that is the component's measured role."* **S25-L5 withdrew that number** and replaced it
with −0.1405 Å at 0.68× MDE. Anyone answering item 12 from the source file gets a **withdrawn
positive**. **Sixth instance of prose asserting a state that does not hold, and the first in shipped
code rather than a report.** Not fixed — `core/` was read-only — and listed in the report's §11.

---

## NOTE 16 (2026-09-20 14:05, lane P closing, S30-L27): **THE BUILT CHAIN, AND TWO CORRECTIONS TO ME**

**The endpoint numbers.** Five ORACLE signs are worth **−0.3259 Å on the built chain: 3.2126 →
2.8867**, 2.05× MDE, 5/5 folds, 92W/34L. `ORACLE_SEPPROF5` gives **2.6791, −0.5335**, 2.29× MDE.
**The prize clears the charter's primary target.**

**Correction 1, and it is a memory I have been misusing.** I transferred a cloud delta to the chain
with the 1.16 coefficient from `operator-consumes-set-mean`. **That law does not map cloud to
chain.** Its `d_set_mean` is the **mean RMSD of the retained 75** (3.5507), not the RMSD of their
average (3.0483) — substituting its own inputs gives **4.21 Å against production's actual 3.21**.
**The measured cloud→chain transfer for a prior correction is 0.92** (0.929 and 0.914 on two arms):
*the projection absorbs ~8% of a cloud gain rather than amplifying it.* Memory updated.

**Correction 2 — my tautology argument is not tight.** I argued that R²(e ~ S) ≈ 0 is near-tautological
for sequence-derived features because a well-fit model's residual is by construction unpredictable
from its own inputs. Lane P: that holds exactly only for a **Bayes-optimal** predictor, and lane M's
audit says **the distogram memorises by 8×** — it is demonstrably *not* well-fit, so its residual
**is** in principle recoverable from its own inputs. **The tautology does not do the work; the null
does.** Accepted.

**And my prediction was falsified on the raw statistic.** I put 3:1 on long-range R² also coming
back ≈ 0. It is **+0.1959**. Nested properly it is **+0.0758** genuinely new (the rest is
calibration and the shared-referent floor — a trap the random-feature control does not catch, and
one I have a written memory about and did not apply).

**Lane P's repair to lane G's enumeration, which is tighter than G1:** *"physics is target-independent"*
is a statement about the **function**, not its **value** — a universal function on a target-specific
argument yields target-specific output. **Physics is an operator, not a source.** The enumeration
collapses from three sources to **two: sequence and library.**

> **The closure: within {sequence, library}, every predictable component of the prior's error is
> common-mode, and common-mode correction is non-identifiable from within the pool. Not "we searched
> and found nothing" — "the set has a structural property that guarantees the search finds the wrong
> thing."**

**The named gap, honestly stated by the lane that produced it:** mutual information is *not* shown
to be zero. What is shown is that the extractable part is the wrong part, **by a theorem about which
part is identifiable.**

**And lane P made, and caught, the same error lane D found in S29:** its first run's verdict strings
were inverted, because `ST.compare` is lower-is-better and explained variance is higher-is-better.
Caught on the first table. *In a lane that had just read the entry about it.* Fifth instance of the
statistic/convention mismatch this sprint.

## NOTE 15 (2026-09-20 13:52, lane G, S30-L26): **BOTH MY PRIORS WERE DIRECTIONALLY RIGHT AND BOTH MY REASONS WERE WRONG** — E2 IS REAL AND DOES NOT REACH THE TAIL; THE CHIRAL ESCAPE IS EXERCISED AND EMPTY

Lane G closed both open questions. Prereg at `78b65521`, committed before the first number.

### (a) Q1 — the divergence claim is CONFIRMED, and it beat the power note registered against it

Lane G registered that a half-split **could not** reach its own MDE unless more than 100% of the
0.022 Å effect sat in one half. **More than 100% does.**

```
whole sample                 -0.0221  [-0.0297,-0.0159]  1.12x MDE  5/5   (reproduced from s16)
rho(d, DISP_rmsd)            -0.316   fold-preserving permutation null p2 = 0.000
high-dispersion half         -0.0429      low half  -0.0012      97.2% of the gain
LENGTH-MATCHED split         -0.0406  [-0.0481,-0.0319]  3.56x MDE  5/5 folds
uniform-effect null          observed gap at pctile 0.001 of [-0.0264,+0.0289]
```

Multiplicity is ~1 comparison, not 3 — the three dispersion variables correlate 0.93–1.000, and
lane G's reproduces **lane F's `F4_top75_rg_sd` at ρ = 1.0000**. It survives dropping the ten
biggest winners (−0.0310, 5/5) and dies only at drop-top-20, where there is no gain left to
localise — a gradient, not ten winners relabelled. Partialling chain length **strengthens** it.

**The mechanism is measured, not assumed:** `ρ(DISP, contraction) = +0.947`, and partialling DISP
collapses contraction's effect from −0.276 to **+0.078**. **Divergence is the primary variable;
contraction is its shadow.** The relax does not repair the bias — *it repairs the geometry the
averaging destroyed*, and the destruction scales with spread.

**But the inference I drew does not go through, and that is the finding.** Divergence and the
production tail are **different objects**:

```
FAIL18 dispersion vs the 108     +0.4583 [-0.3845,+1.0472]  0.47x MDE  NOT MEASURED
relax gain on FAIL18  -0.0113    vs the other 108  -0.0239   <- the WRONG direction
targeting the divergent half     -0.0215  against  -0.0221   <- buys nothing
```

> **E2 is a real, mechanistically explained, divergence-graded effect that does not reach the
> stratum the opening arithmetic cares about.** 0.69% of baseline, and its restraint constant still
> has no native-free selection rule.

### (b) This withdraws one of MY Appendix A entries — a withdrawal that is itself withdrawn

Appendix A recorded lane L lowering "the benefit concentrates on divergent pools" from 2:1 to
roughly even, on lane F's dispersion null. **That downgrade was wrong**, and so was lane G's own
2:1-against. The correction that matters is *why*:

> **Two lanes agreeing did not make a prior.** Lane L moved on lane F's null, which was about a
> **different outcome** (filter failure, not relax gain). Lane G moved on the common-mode argument,
> which was about a **different mechanism** (bias repair, not geometry repair). Two independent
> revisions in the same direction, each resting on an object that was not the one under test, read
> as convergent evidence and were not evidence at all.

### (c) Q2 — my pre-check refuted lane G's own mechanism, and made the negative stronger

I asked for the variance of the chiral coordinate across the manifold *before* any contrast, on the
grounds that a provably-chiral functional can still be near-constant on a set containing no mirrors.
Lane G predicted it would be. **It is not:**

```
occupancy (sd over candidates / sd over candidates + their mirrors)
WRITHE 0.685   CHIRAL3 0.953   CHIRAL3_LONG 0.965   |   DIS 0.345 (pool), 0.231 (ladder)
```

**The chiral axis is more fully exercised on our manifold than the shipped cost is.** The escape
class is not empty in practice — so the null below is a real negative, not a degenerate coordinate.

```
WRITHE anchor contrast   +0.0405 [-0.0255,+0.1015]  0.41x MDE  3/5 folds   (bar +0.10)
CHIRAL3 -0.0363   CHIRAL3_LONG -0.0394   (wrong sign)
max-over-3 sign-flip null:  observed 0.0405   null MEAN 0.0408   p_max 0.430
WRITHE minus its achiral twin |WRITHE|:  +0.0170  0.29x MDE  CI includes zero
```

The best chiral channel sits **on its own null to three decimal places**, and whatever WRITHE does,
its reflection-invariant shadow — a `D`-functional by G1 — already does. The in-code audit
(reflection flips X at 0.00e+00, leaves |X| and D at 0.00e+00, rotation invariance 4.19e-13)
**caught a real bug**: the first Klenin–Langowski sign term was not rotation-invariant (rot_err
2.03) and the assertion stopped the run before any number existed.

### (d) The positive lane G found, chased, and killed itself — the third cross-kind confound

WRITHE's preference contrast is **+0.1641 [+0.1016,+0.2282], 2.17× MDE, 5/5 folds, max-null
p = 0.000**, with `pref_near = 0.771` — **clearing both of lane R's registered preference clauses
that none of its 43 channels cleared.** It is **cross-kind**, in exactly the way S30-L1 withdrew
S28-L48: `pref_pool`'s control keeps **deposited** coordinates while the near rungs are **ideal
rebuilds**, and a raw shape statistic separates those constructions whether or not it sees
nativeness. The kind-matched statistic settles it — the rebuilt native's percentile inside its own
ladder, chance exactly 0.5:

```
WRITHE 0.6061   |WRITHE| 0.6732   both WORSE than chance, 5/5 folds   |   DIS 0.2876
```

A channel ranking the native at the 61st percentile of the native's **own** perturbations is
separating production from everything else, not recognising nativeness. **Withdrawn by lane G before
it was quoted anywhere.** Third instance of the cross-kind confound this sprint.

### (e) What it closes

> Lane R proved the null for every **per-residue** channel is a theorem. Lane G proves **the null
> for every non-chiral channel is a theorem too** — leaving exactly one family, which was built,
> is genuinely exercised, and is empty.

**Scope, so it is not over-quoted:** G1 bounds **single-structure** channels only — not
set-referenced ones (bucket b, closed separately), and not a fourth reference this project does not
have. **The theorem does not depend on chain length; only the emptiness does.** The survivor family
is worth retesting at 40+ residues, where a chain can actually cross itself.

## NOTE 14 (2026-09-20 13:43, lane P, S30-L25): **THE ONE MEASUREMENT CAME BACK NULL — AND THE SECOND HALF REPRICES THE SPRINT INTO A FIVE-BIT QUESTION**

The whole B2 bound reduced this morning to a single estimable number, `rho_max = sqrt(R2(e ~ S))`,
with bars pre-registered at 1.96% / 12.82% / 39.44% before the regression existed. It is measured.

### (a) The registered answer is the null

**R² = 0.0083** out of fold — `rho_max = 0.091` against the **0.358** that 3.00 Å requires and the
**0.628** that 2.50 requires. Best of 8 arms; matched-dimension control −0.0012; excess **+0.0095 at
0.75× MDE, NOT MEASURED**. Nothing in {main, full, no-distogram} × {weighted, unweighted} ×
{126, FAIL18, other-108, tailA, tailB} clears the bar. P1–P5 all hold.

And the line that matters more than the R²: **applying the fitted displacement out of fold emits
3.0519 Å against production's 3.0483 — it is very slightly WORSE.** The implied-endpoint conversion
(ρ → 3.0338) flatters a no-op by 0.018 Å. That is the sprint's **fourth checklist entry**:

> **An implied-endpoint conversion from ρ is an upper bound attained only by a perfectly calibrated
> correction. Always also APPLY the correction and measure the endpoint.**

Every lane that quoted a √(1−ρ²) number this sprint, me first, was quoting the flattering side.

### (b) Why the null is *more* conclusive than a failed search

A trained model's residual is by construction what it could not predict from its own inputs. So for
any S that is a function of the **sequence**, R²(e ~ S) ≈ 0 is near-tautological for a well-fit
predictor — and ESM-2 650M is one. The sequence arms were never the informative ones. **The only
genuinely non-tautological arm in the design was the POOL**, which the distogram never saw — and it
came back at **+0.0395 Å with the wrong sign**. Regressing the prior toward the pool *hurts*. That is
`pool-error-is-68-percent-common-mode` arriving at the prior, a third independent site.

### (c) The sprint's cleanest positive, at 10.17× MDE

The native-free basis captures **82.73%** of the ORACLE error in 14 directions against a
matched-dimension random frame's **44.44%** — excess **+0.3829, 10.17× MDE, 5/5 folds, CI
[+0.369, +0.392]**.

> **The subspace is free. The sign is the whole problem.**

### (d) The reframe — and it says the target is REACHABLE, just not from here

Lane L's 0.525 Å separation-profile prize reproduces on the pinned benchmark, slightly larger:

```
ORACLE_FULL (perfect distogram)   2.2379   -0.8104   3.28x MDE  5/5  117W/9L
ORACLE_SEPPROF5  (5 numbers)      2.4743   -0.5740   2.52x MDE  5/5  106W/20L   70.8% of the full gap
ORACLE_PERRES    (n numbers)      2.4769   -0.5714                              extra params add NOTHING

ORACLE SIGN + leave-fold-out magnitude   -0.3567   2.29x MDE  5/5  93W/33L    62.1%
ORACLE MAGNITUDE + leave-fold-out sign   -0.0733   0.32x MDE  NOT MEASURED    12.8%
```

Sign-corruption ladder (the mandatory null from `error-coherence-decides-correctors`):
acc 1.0 → −0.574, 0.9 → −0.472, 0.8 → −0.274, 0.7 → −0.251, 0.6 → −0.075. **The deployable
leave-fold-out sign performs like accuracy ≈ 0.6 — barely above chance — and ≈ 0.8 is needed.**

Both deployable arms are zero and the mechanism is explicit: 5 **global** numbers buy +0.0036 Å
(the profile's per-target dispersion is **9× its mean**, so the real systematic bias is swamped), and
5 numbers estimated **from the pool** buy +0.0395 Å, wrong sign.

**Where the prize is, on non-circular strata:** |i−j| ≥ 7 carries **68%** of it — lane R's
"absent from local (ΔR² −0.089), abundant in global (+0.600)" arriving from a different instrument —
and it concentrates **~4×** on both of lane F's filter-independent tails, not just on the circular
FAIL18.

**So the sprint does not close with "unreachable".** Five ORACLE signs are worth −0.3567 Å on the
cloud: 44% of a perfect distogram, production 3.0483 → 2.69. The closing statement is:

> **The target is reachable with five bits per target on long-range pairs, concentrated on the tail.
> Nothing in the distogram-plus-pool information set supplies them at the required accuracy.**

### (e) My registered prediction, made before lane P's built-chain number arrives

The terminal slope (`d_out = 1.16·d_set_mean + 0.04·d_set_best`) transfers −0.3567 to **≈ −0.41 on
the built chain, 3.2105 → ~2.80**, which would clear the charter's primary target. **I predict the
measured chain delta lands in −0.35 to −0.45.** If it lands materially short, that is *not* noise:
correcting the prior changes the score, which changes selection, which changes the set the terminal
averages, so the 1.16 slope is not clean here and a shortfall is the selection-feedback effect —
its own finding. I will quote lane P's measured number, not the transfer.

### (f) Lane P and lane G have met, and jointly they may close it

Lane G proved (G1, classical MDS) that target-specific information has **exactly three** sources:
the predicted distogram, the pool, and universal physics — the third **target-independent by
construction**, hence unable to supply a per-target sign. Lane P has now measured the first two at
**zero** and **negative**. If G1's reference argument transfers from *scoring a structure* to
*correcting the prior* — different operators, so the transfer is not free — then the missing
covariate provably does not exist inside this project's information set, and the only escape is a new
observable, which S27 already priced and closed (shifts on 54/126; ORACLE-perfect torsions still
leave 2.021 Å against the 55 needed). **I have asked lane P to attack the transfer rather than
accept it.**

### (g) An unregistered comparison, flagged as such

Lane T's law gives **0.376 Å per INDEX bit** at D = 3.05; lane P's −0.574/5 gives **0.115 Å per
prior-sign bit**. Roughly 3×, and lane T independently measured candidate indexing beating subset
cardinality by 3×. Two different bit currencies, so it may be coincidence — but if it is not, **the
readout is where bits are worth most**, which is load-bearing for the next sprint.

## NOTE 13 (2026-09-20 13:38, lane D, S30-L24): **THE METER'S NEW CONTROL WITHDRAWS AN S29 NUMBER ON ITS FIRST RUN** — S29 DREW THE MAXIMUM OF ITS OWN EIGHT

I asked lane D to extend the meter with a `verify` verb and an R-draw control. The control paid for
itself immediately, and what it withdrew is *ours*.

S29 published its built-chain preference contrast against **one** matched random-signed structure,
seed 0. Lane D ran the same contrast against **eight** draws:

| basis | S29 published (1 draw) | S30 (8 draws) | single-draw range | verdict |
|---|---|---|---|---|
| **built chain** | **+0.0635** (0.92× MDE) | **+0.0357** (0.56× MDE) | [+0.0159, +0.0635] | **WITHDRAWN** — below the 0.7× floor, fold CI spans zero, folds 3/5 |
| **CA point cloud** | +0.1746 (1.87× MDE) | **+0.1716** (1.84× MDE) | [+0.1587, +0.1905] | **CONFIRMED**, 5/5 folds |

**S29's seed-0 draw was the maximum of its own eight.** Not near it — the maximum. The published
number was the single most flattering value the control could have taken.

### Why the two bases split, which is the transferable part

Relative draw noise is **59% on the chain against 27% on CA**. On the chain the shipped cost prefers
production to nearly everything — every preference in the sweep sat between 0.008 and 0.056 — so the
control *is a rare event*, and a rare event's single realisation is mostly variance.

> **A single-draw control is least trustworthy exactly where the effects are smallest — which is
> where this project's remaining effects live.**

That is the whole of it. It is not a fact about S29's carelessness; it is a fact about where we now
operate. Every contrast still standing in this project is in the 0.02–0.15 band.

### The rule, in lane D's words, and it goes in the checklist

> *"A control that is itself a random draw needs its own draw distribution reported, and the number
> of draws must scale inversely with the effect size. Eight separated these two cases; one did not."*

### A nuance that stops a false alarm

Lane D also flagged that the preference statistic is a 0/0.5/1 indicator with modal value 0, so its
**median is uninformative by construction**. A median paired difference of 0.0000 there is *not* the
median-vs-mean early warning firing — it is what that statistic's median always does. Use the
uniform-effect-null concentration test instead. This is the third time this sprint that a statistic's
transform and its null have had to be matched by hand (see NOTE 10(a), my own `7 − log₂r` error).

### What it does and does not touch

It does **not** touch the sprint's own results: every S30 contrast was run against its own control in
its own space, and lane R's 43-channel verdict used a max-over-channels sign-flip null. It touches
one **S29** number, and it means the S29 report's built-chain preference row must be annotated rather
than cited. Per contract rule 13 the historical artefact stays untouched; the annotation lives here,
in S30-L24, and in the report's §3 closures table.

## NOTE 12 (2026-09-20 13:31, lane X, S30-L20): **NO GENERATOR CAN IMPROVE THE TYPICAL MEMBER IN THE SENSE THAT PAYS** — AND MY PREMISE WAS BACKWARDS

I asked lane X to test the one construction its own law would admit: a generator whose **typical**
member is better, rather than whose **best** member is better. It measured the set mean first, as
instructed, and found **the set mean is two quantities and only one of them converts.**

For a coordinate-average terminal over m members, the common-mode identity gives per target:

```
set_mean^2  ~=  B^2 + S^2        B = RMSD(set average, native) = THE ENDPOINT
                                 S = the set's spread about its own centroid
```

**A set mean improved purely by concentration is worth zero, by algebra.** Only B converts, and
**B is the endpoint itself.**

```
arm            set_mean   B(endpoint)     S    avg_gain |  d_set_mean     d_B      d_S
POOL (top-75)    3.5507      3.0483    1.5770   0.5023  |      --         --       --
T0_helix         3.8663      3.7892    0.5676   0.0771  |   +0.3156   +0.7408  -1.0094
T1_blind         4.0218      3.2435    2.2242   0.7784  |   +0.4712   +0.1951  +0.6472
T2_restype       3.9563      3.2065    2.1267   0.7498  |   +0.4056   +0.1581  +0.5497
T3_pool          3.5916      3.1752    1.4791   0.4164  |   +0.0409   +0.1269  -0.0979
```

**Four reasons, independent:**

1. **Nothing in the record improves the set mean at all** — every `d_set_mean` is positive, three
   at 5/5 folds with CIs excluding zero. The best is a **statistical tie** (T3_pool +0.0409,
   0.57x MDE) against my stated requirement of **-0.2 A**.
2. **The concentration half is worth zero, as the algebra says.** `avg_gain = 0.4143*S - 0.1559,
   r = 0.8854, n = 630` — **the terminal's entire value is spread extraction**. Within target,
   `corr(S, B) = +0.0851`: concentration is **orthogonal** to bias.
3. **The extreme case proves it.** `T0_helix` — constant alpha-helix plus jitter — **is** the
   typical-good generator I described, and the most concentrated source ever built here
   (S 0.568 against the pool's 1.577). **It is the worst endpoint in the record: 3.7892.** Its
   averaging gain collapsed to 0.077 and the endpoint fell back onto B.
4. **The cells that win on both do not transfer.** Split-half **+0.1597 [+0.0676, +0.2687]** —
   wrong sign, CI excluding zero.

> **MY PREMISE WAS BACKWARDS.** I wrote: *"every sampler was built to be diverse; nobody has built
> one to be typical-good."* Somebody did, and **diversity is the correct design**: the
> largest-spread arms have the largest averaging gains and the best generated endpoints.
> **For an averaging terminal, spread is the raw material, not a defect.** The samplers were not
> built wrong.

**Two honesty items lane X volunteered, and both matter for quoting it.** Its counting falsifier
**essentially tied its bar — 24.80% against 25%, a 0.20-point miss — and decided nothing**; the
verdict rests entirely on the transfer arm. And the raw within-target `corr(S,B) = -0.4744` is
**driven by T0_helix alone**; drop that arm and it is +0.0851. Both satisfy its registered
prediction, but **the strong negative is an artefact and must not be quoted** — the durable number
is the orthogonality.

**And it amended its own admission condition to be STRICTER**, not looser: S30-L10's
`d(set mean) < -(0.32..0.37)*d(set best)` treats one channel as two, so it becomes
**`ADMIT G iff d(bias B) + 0.298*d(set best) < 0`**. The cheap half was free all along.

**Why this is a ceiling and not a miss.** "Improve the typical member" decomposes into a **spread**
half (free, already extracted, orthogonal to the endpoint) and a **bias** half (worth everything,
but it *is* the endpoint). **Lane L's non-identifiability closes the second half structurally.**
The only half of "typical-good" that pays is the half no source change can move — and lane L's
falsifier had already failed before any endpoint run, at a provenance cosine of 0.9432 against a
0.9330 within-source control. **No endpoint compute was spent on this closure.**

## NOTE 11 (2026-09-20 13:28, lane R, S30-L19): **L11 IS ANSWERED. NATIVENESS IS NOT RECOGNISABLE FROM SINGLE-STRUCTURE GEOMETRY** — AND THE NULL IS A THEOREM, NOT A MISS

The question I judged to sit underneath the sprint, answered on an instrument built to remove the
confound that invalidated the project's previous attempt. **F-R1 does not fire at n = 126.**

**The instrument.** Every rung is an ideal-geometry backbone built from (phi, psi) — identical bond
lengths and angles at every rung, no projection, no coordinate average, no contraction. Perturbed
residues draw torsions from the fold's **leakage-safe** Ramachandran table. Two ladders: anchored
on the native's own torsions, and on a random real pool member. **Kind, local realism and
perturbation budget are all matched; only nativeness varies.** Floor is the torsion rebuild at
**0.347 A**, not 0, and it travels with every sentence.

**The split verdict, which is more informative than a flat null.**

- **ORDERING survives.** DIS reaches **+0.347**, **+0.134 above its own anchor control**,
  p_max 0.000. Channels can order structures by nativeness at fixed budget.
- **PREFERENCE does not, on ANY of 43 channels.** And the leave-fold-out combination exposes why in
  one line: it prefers a **0.55 A structure to production on 93.0% of targets** — and prefers a
  **random pool member on 100%** and a **3 A rung on 100%**. Its margin over the controls is
  **−0.070 [−0.110, −0.028]**. It is not detecting nativeness; it is detecting *not-production*.
- **LEG_torsion does not survive — but "at chance" was MY wording and lane R withdrew it**
  (corrected 2026-09-20 13:30). Its anchor contrast is **+0.024 with fold CI [+0.010, +0.038],
  excluding zero**: it *does* order nativeness slightly above its own anchor control. The correct
  and still-decisive statement is that this is **a quarter of the registered +0.10 margin**, its
  preference is **0.503 against a 0.698 pool-member control**, and it ranks the native at the
  **exact median (0.503) of the native's own perturbations**. S29 section 12.0's last
  outside-class-M hope closes negatively, but not by being noise.

**AND THE MECHANISM, WHICH IS THE PART THAT MAKES THIS DURABLE.** The derivation I asked lane R to
make explicit came back measured: the RMSD signal is **absent from local features (delta-R^2
−0.089) and abundant in global ones (+0.600)**.

> **So the null for every per-residue channel is a THEOREM on this instrument, not an empirical
> miss.** A structure can be locally perfect everywhere and globally wrong; matched local
> statistics make local channels blind by construction, and the measurement confirms it.

**Lane R discounts its own result correctly**: its registered prior was "F-R1 does not fire, about
4 to 1", so the null confirms its own expectation and it said so rather than presenting a
confirmation as a discovery.


**Three more numbers from lane R's close that sharpen this.**

- **The largest preference effect in the whole library points the wrong way:** DIS prefers
  production to a 0.55 A structure on **94.2%** of targets. That is S28-L48 with the kind confound
  removed -- **sharper, not weaker**. My note 1 said removing the confound might weaken it; it did
  the opposite.
- **The local block is ORACLE-advantaged and still adds nothing.** It is handed the per-residue
  circular deviation *from the native anchor* and still gives delta-R^2 **-0.089** against global's
  **+0.600**. So RAMA, LEG_torsion, DSSPHB, CAGEO and HP are **blind by construction** -- a sum of
  per-residue terms cannot see a lever arm.
- **Resolution: coarse triage only.** DIS concordance inside the near band is **0.520** at
  |delta| = 0-0.25 A, rising to 0.822 only above 4 A. **Resolving 0.5 A needs a channel ~4x the best
  in the record**, and inside the near band LEG matches DIS -- the distogram's edge is entirely
  coarse.

**AND WHAT IT DOES NOT FORECLOSE, which lane R was careful to separate.** D1 says the information
**is** in the global shape -- held-out R^2 **0.911**, ORACLE and native-informed. What is missing is
a **native-free globally-reaching channel**. That is a **supply problem, not an impossibility**, and
it is a different problem from the one this sprint has been attacking. Nothing in lane R's result
touches the retrieval or prior lever.

### WHAT THIS CLOSES

The charter said a negative here "would be one of the most important results the project could
produce", and it is the one that composes with everything else closed today. If nativeness cannot
be *preferred* from a single structure's geometry, then no Hamiltonian, no encoding, no ansatz and
no amount of search can be pointed at it — the objective would be optimising toward a target it
cannot see. Combined with:

- the field library spending **58% of its two directions on a radial component that is
  anti-aligned with the answer on hard targets** (D, confirming L's derivation),
- the combination of all 21 fields reaching **rho = 0.169 with the native in hand**,
- **ordering** skill existing while **preference** does not,

the picture is consistent: **there is orderable signal and it cannot be converted into a
preference, and the preference is what a cost function needs.**

**What is still open is exactly one thing:** lane P's out-of-fold `R^2(e ~ S)`. Lane R's D1 says the
answer must live in **global** features if it lives anywhere, which is a constraint lane P should
have — and lane P's basis is already built from the pool's principal directions, which are global by
construction.

## NOTE 10 (2026-09-20 13:28, lane D closing): A CLAIM I REPORTED IS WITHDRAWN, LANE L's DERIVATION IS CONFIRMED HARDER THAN IT ASKED, AND THE METER WAS BLIND ON THE REPORTING BASIS

### (a) I reported "the median target is worse than chance" and it is wrong

S29's "the deployed score is at chance" **survives and is strengthened**. Its *supporting paragraph*
— which I quoted to the user as "the mean flatters it; the median target is worse than chance" —
**is withdrawn.** The statistic `7 − log2 r` is **right-skewed**: null mean 1.4050, null **median
0.9888**, null **win rate 62.5%**. So "below the random baseline on 82 of 126" is simply *what a
random ranking does to itself* — the null expects 78.8, z = +0.60 — and the −0.575 median gap is
−0.416 under the null. Against a simulated null **all five statistics are at chance** (p = 0.16 to
0.75), and the rank distribution passes KS (p = 0.571) and chi-square (p = 0.303).

**The durable lesson, lane D's words:** *when a statistic is a nonlinear transform, its mean,
median and win-rate have three different nulls — never read one against another's baseline.* I read
the median against the mean's baseline and reported it twice.

### (b) Lane L's radial prediction is CONFIRMED, and the payoff is better than the prediction

Lane L derived that one of the Gram's two effective directions must be the scale/radial one, and
offered the four-line test with its concession pre-stated. Measured:

- the radial direction carries **58.0% of the Gram trace** [+0.545, +0.603]
- the field set's **dominant principal direction IS the radial one** — cos **0.947**, median 0.984;
  **94.8% of lambda_1 is radial**
- removing it raises stable rank 1.705 -> 2.642 and drops lambda_1's share 61% -> 40%

**And the payoff, which lane L did not predict:**

> **cos(direction to the native, radial) = −0.0675 overall and −0.2524 on FAIL18.**

**The field library spends the majority of its two available directions on a component that is
orthogonal to the answer in general and ANTI-ALIGNED on the hard targets.** Lane D is careful that
this does not license deflating scale: the residual 42% has no large eigenvalue, and S30-L6 already
priced the whole span.

### (c) The meter was blind on the reporting basis, and my brief asked for an impossible check

Two defects in what I handed lane D:

1. I told it to run `selftest` and confirm it reproduces S29's baselines. **`selftest` is a
   synthetic 8-residue check** (gradient cosine ~ +1, tie conventions) that runs in 0.4 s and
   **would pass on a meter whose every dev-set number had drifted.** The instruction could not have
   been satisfied by that command. The reproduction check now exists separately as
   `python s30/s30_D_meter.py verify`, asserted in code.
2. The 0.3688 anchor I listed is **DIS_SURR's**, not the shipped cost's — `DIS` is 0.3676. I put it
   beside five `DIS` numbers.

**And the material one, which lane D found rather than inherited:** `--basis chain` was **unusable**
— 0 of 126 ladder-cache files carried chain projections, so any chain run raised
`FileNotFoundError`, and `chain-s28rows` only ever serves the 31 scorer names S28 already
evaluated. **The instrument that gates every new cost could measure new costs only on the CA point
cloud.** That is not cosmetic here: on the shipped cost the two bases differ in the **sign** of the
headline charter ladder (**+0.2603 CA against −0.0921 chain**) and by 2.2x on the S28 ladder. Fixed.

**All ten anchors reproduce to < 5e-4 on a cold re-run**, and the extended meter now carries
per-target distributions, a fold-clustered FAIL18/108 contrast, matched random-signed structures as
a first-class control, the aggregate cosine null, a PASS/BLOCK gate at the 0.7x/1.0x rule, and a
multiplicity register. **The gate BLOCKs the shipped cost on its own ladder** — the honest verdict.

One caveat every lane inherits: **fold 0 contains no FAIL18 target** (1:6, 2:2, 3:4, 4:6, 0:0), so
every FAIL18/108 CI draws from four clusters and is wide by construction.

A first CI on a split S29 printed as bare means: **under the shipped cost the native sits at the
78.4th percentile of its own pool on FAIL18 against the 29.8th on the other 108** — +0.485,
3.48x MDE, 4/4 folds.

### (d) Two checklist entries earned today

- **A matched control in the right space does not rescue a stratum defined by the outcome** (lane D,
  from S30-L8).
- **When a statistic is a nonlinear transform, its mean, median and win-rate have three different
  nulls** (lane D, from S30-L4).

## NOTE 9 (2026-09-20 13:26, lane D, S30-L8): **THE SPRINT'S FIRST MECHANISM HAS ITS HEADLINE WITHDRAWN.** THE STRATUM *IS* THE OUTCOME

The adversary has gone after lane F's S30-L2 — the entry I reported as the sprint's first real
mechanism — and it is right. **F1c's FAIL18 row must be withdrawn as an effect estimate.**

**The definition, from `s12/instrument.py:271-278`:**

```python
sub  = np.asarray(rec["sub"], int)              # the SCORE's top-75
band = np.where(rr <= rr.min() + BAND)[0]       # pool members within 1.5 A of the pool best (ORACLE)
if not np.isin(band, sub).any():
    zero.append(t["pdb"])
assert set(zero) == set(FAIL18)
```

> **FAIL18 is DEFINED as the set of targets on which the 500 -> 75 filter retained no in-band
> member.** S30-L2 then measured, *on that set*, how much worse the filter's retained set is than a
> random 75 in ORACLE best. **The selection predicate is a lower bound on the measured quantity's
> first term.**

**Lane D reproduced the number independently** — −1.7622 A against lane F's +1.7674 in the opposite
sign convention — and then decomposed it:

```
stratum by `keep` = in-band members surviving the filter     n      effect
keep = 0   <- this stratum IS FAIL18                         18    -1.7622
keep = 1                                                      3    -0.0363
keep = 2-3                                                    1    -0.8654
keep = 4-8                                                    8    -0.0039
keep >= 9                                                    96    +0.0380
```

**There is no gradient.** The effect is a **step at the selection boundary**, not a continuum in
recall — one retained in-band member removes 98% of it. Among the 108 non-FAIL18 targets,
Spearman(keep, effect) = **−0.046**. Outside the 18 the filter is **+0.025 A**, and the all-126
effect is **100% the 18**.

**How much is forced arithmetic.** By definition `top_best − pool_best >= 1.5` on these 18
(measured 2.393); the random arm is unconstrained at `rand_best − pool_best = 0.631`. So the
predicate alone forces **|effect| >= 0.869 A — 49% of the headline.** The remaining 0.893 is not an
unbiased estimate either: it is expected *exceedance above a selection threshold*, which truncation
inflates by an amount this design cannot measure.

**Lane D's verdict, and the sentence is the durable output:** *"Lane F did the hard parts right --
prereg before the numbers, matched random subsets in each operator's own space, fold CIs, a
random-18 null -- and none of those defences reaches this one, because **the stratum itself is the
outcome.**"*

That is a failure mode distinct from every one on the project's list. Lane F **did** check
circularity — it noticed FAIL18 is defined by *production* and produced two filter-independent
tails. What it did not catch is that FAIL18 is defined by **the filter's own recall**, which is
precisely the quantity being measured.

### WHAT SURVIVES, AND IT IS NOT NOTHING

- **F1a is untouched**: the tail is **not pool-limited**. ORACLE best of the FAIL18 pools is
  2.2842 A, 13 of 18 under 3.00. That is a property of the pools, not of the filter's recall.
- **The filter effect on filter-independent tails survives**: +0.6708 (p = 0.0148) and +0.6320
  (p = 0.0300) on worst-18-by-pool-mean and worst-18-by-ORACLE-best. Lane F reported these itself
  and told me to quote them beside the headline; **they are now the headline.**
- **S30-L17 survives and is the strongest tail result**: the score's Spearman with ORACLE in-pool
  RMSD is +0.6446 on the 108 and **+0.1066 on the tail with the CI including zero**, and it
  **replicates on all three tail definitions** (+0.38, +0.34). Nothing about that rests on the
  FAIL18 predicate.
- The **mediated chain** — distogram shape error -> in-pool ranking skill -> filter set-mean
  benefit -> emitted RMSD — rests on S30-L17, not on the withdrawn row.

### WHAT I GOT WRONG IN REPORTING IT

I told the user the filter is "worse than chance on every single tail target, 0W/18L, random-18
null p = 0" as the sprint's first mechanism. **Roughly half that number is forced by the definition
of the set it was measured on**, and the rest is truncation-inflated by an unmeasurable amount. The
correct statement is the +0.63 to +0.67 pair on tails that are not defined by the filter.

## NOTE 8 (2026-09-20 13:20, lane T's closing line): **THE WHOLE BOUND COLLAPSES TO ONE OUT-OF-FOLD REGRESSION**

Lane T closed with a reformulation that turns the sprint's central question into a single cheap
measurement. B2 -- "no native-free field reaches rho > 0.14" -- is restatable as one estimable
number:

> **rho_max = sqrt( R^2( e ~ S ) )**

where `e` is the oracle error (production's answer to the native) and `S` is any native-free signal.
The bound becomes a variance-explained threshold:

```
rho = 0.14   B2's ceiling, the best field ever built   <->   R^2 =  1.96%
rho = 0.358  3.00 A, the charter's primary target      <->   R^2 = 12.8%
rho = 0.628  2.50 A, the charter's ambitious target    <->   R^2 = 39.4%
```

**So the question the whole sprint has circled -- is there enough information anywhere -- becomes:
can any native-free feature set explain 39.4% of the variance of the oracle error, out of fold?**
That is **falsifiable by one regression rather than by a field survey**, and it subsumes the survey:
throw everything available at `e` at once and read the number.

**And the bar is lower than it looks, because the subspace is free.** Lane T measured that the error
is **8x more concentrated in the pool's leading deviation directions than isotropy predicts** --
six directions contain a rho = 0.779 point (1.91 A). The pool's own PCA is native-free and
per-target, so regressing *within* that subspace costs nothing and drops the requirement to
**0.855 bits for 3.00 A and 3.78 bits for 2.50 A**, against **0.076** from the best field on record.

This is now the highest-value measurement left in the sprint and it is lane P's. Registered with it:
the bars above, a **matched-dimension control R^2 beside every number** (lane D's per-target ORACLE
combination hit rho = 0.949 where a random 21-dimensional subspace hit 0.809 -- almost all
dimension, no signal), a re-run with distogram-derived features removed if anything clears, and a
split by **shape versus scale** and by stratum, because lane F localised **83% of the tail's error
to shape** and refuted scale as the mechanism.

---

## LANE T CLOSED. Final ledger: S30-L9, L14, L15. `s30/THEORY.md`, 550 lines, 10 sections.

Durable, and nothing about them is an order statistic: **T1** and **T1b** (algebra); the
**value-of-a-bit law** at R^2 0.9983; the **codebook 5.2x**; the **torsion 48-vs-7** arithmetic;
**M5**'s subspace projection; the corrected **within-subspace dictionary**; and the combination
formula **rho_comb = rho_0/sqrt(s_1)**, registered in advance and accurate to **0.008**.

Withdrawn by its own hand: four claims, three toward its own hypothesis.

## NOTE 7 (2026-09-20 13:19, lane T): **IT RETRACTED ITS OWN HEADLINE**, ON A NULL LANE Q TOLD IT TO RUN

Lane Q asked lane T to price its K = 5,000 direction search with `best_of_k_within` before
reporting. Lane T ran it, using a **common direction bank across targets** (columns made comparable
by a deterministic native-free sign convention, so column k is the same *rule* everywhere and the
statistic is well posed). The result destroys its own result:

```
                    observed   across-target null   accounted   split-half transfer   k_eff
HALFSPACE K=5000    -1.7881        -3.5038            196%        -0.3319 (19%)        118
QUADRIC   K=5000    -1.5663        -2.9368            188%        -0.3148 (20%)        110
```

**The null exceeds the observed in both.** And the decisive follow-up, because 19% transfer is not
zero and could be misread as a route: the halfspace direction bank's own mean is 3.8550 A, so the
**transferable rule** (K-mean + split-half) gives **3.5231 A against PREFIX's 3.0507 — +0.4724,
1.88x MDE, 5/5 folds, 39W/87L, WORSE.**

> **"HALFSPACE - PREFIX = -0.9816 A, BETTER" is WITHDRAWN as a statement about the class's value.
> The halfspace class is REACHABLE, not EXPLOITABLE, and as a RULE it is negative.**

Lane T identified its own error precisely: its matched random-subset null asked whether a
*structured* class beats an *unstructured* one at equal budget, not whether **the winning direction
is the same direction twice**. That is the project's most repeated error arriving in a new costume,
caught by another lane and confirmed by the lane that made it.

**What is unaffected, each checked rather than assumed:** the QUADRIC-minus-HALFSPACE contrast
(like-for-like at matched budget; a null shifting both does not touch it), M5's subspace projection
(no maximum over draws anywhere), and **the six-coefficient conclusion, which this strengthens** --
a non-transferable direction in a six-dimensional subspace *is* six per-target coefficients.

**AND THE COMBINATION FORMULA HELD TO 0.008.** Lane T registered `rho_comb = rho_0 / sqrt(s_1)`
with the field count cancelling. Against lane D's independent measurement: the formula gives
**0.174** from D's best single field and **0.162** from the 0.1128 T registered with, against a
measured **0.1693**. P5a held on both clauses (s_1 = 0.486, per-target stable rank 1.681); P5b held
(0.1693 inside [0.11,0.22], point estimate 0.15). **P5c did not hold and lane T says its mechanism
was wrong** -- it predicted the in-sample 21-parameter arm would inflate to 0.25-0.45; lane D tuned
its ridge by nested CV *inside* the training folds, so the arm degenerated to 0.0124 instead.
Lane D's arm was better built than lane T anticipated.

> **The combination ceiling is `rho_0/sqrt(s_1)` -- set by the Gram's RANK, with the field count
> cancelling exactly. Eleven significant fields spanning two directions combine like two.
> The Gram is a one-line test to run before anyone builds a twenty-second field.**

**The stable-rank scope is adopted in full**: every appearance now carries the feature space in the
same sentence -- pair-distance 1.859/1.862, **coordinate 3.404/3.619, k90 = 11.2 not 5.6**. Lane T's
2.0 threshold fires for a distance-map lift and **does not** fire for a coordinate-space one; its
verdict is unchanged but now rests on direct measurement rather than on the threshold, and it said
so rather than letting the threshold carry weight it cannot bear.

**Lane T's running total of its own claims that came back wrong today: four** -- P2b (2-4x,
flattering), P4c (64%, flattering), the unregistered "12 bits" aside (withdrawn, corrected *down*
to 3.78), and now the halfspace headline (retracted on a null it should have run itself). Three of
the four ran toward its own hypothesis. What it still stands behind is what no null touches: T1 and
T1b, the value-of-a-bit law (R2 0.9983), the codebook 5.2x, the torsion 48-vs-7 arithmetic, M5, and
the corrected within-subspace dictionary.

## NOTE 6 (2026-09-20 13:17, lane D + lane T): **THE COMBINATION QUESTION IS ANSWERED. NO.** AND LANE T'S REGISTERED PREDICTION HELD

The measurement I called the sprint's central one is in, and it closes the direction.

**Lane T registered, before the numbers existed:** `rho_comb = rho_0 / sqrt(s_1)` — the combination's
cosine is the best single field's divided by the root of the Gram's PC1 share, **and the count k
cancels**. Its P5a put the 21x21 Gram's stable rank in [1.3, 3.5] and PC1 share in [0.45, 0.80];
its P5b put the **leave-fold-out ORACLE-optimal combination at [0.11, 0.22], point estimate 0.15**,
falsified above 0.25; its P5c warned the in-sample figure would inflate to 0.25-0.45 and must not
be read as a ceiling.

**Measured** (`s30/results/s30_D_gram.json`, n = 119 of 126, ORACLE):

```
Gram stable rank                    2.057      (T predicted [1.3, 3.5])
Gram top-3 share                    0.600      (T predicted PC1 [0.45, 0.80])
per-target stable rank              1.681

best SINGLE field                   rho 0.1214
ORACLE GLOBAL combination           rho 0.1693   implied 3.0253 A
LEAVE-FOLD-OUT combination          rho 0.0124   implied 3.0694 A
LFO equal-weight over the 11 sig.   rho 0.0948   implied 3.0558 A
```

**Three readings, in order of how much they matter.**

1. **Even with the native in hand and one global weighting, the 21 fields combine to rho = 0.169.**
   The bound needs **0.358 for 3.00 A**. The combination is not merely short, it is short by more
   than a factor of two on a quantity that enters as its square.
2. **Deployably it is worse than the best single field** — LFO 0.0124 against 0.1214. Fitting 21
   coefficients leave-fold-out destroys the signal rather than combining it. That is the
   incidental-parameter result arriving at the level of *weights*.
3. **The orthogonality hope is dead on arithmetic.** Eleven fields at 0.11 would need a Gram of
   full rank 11 (PC1 share 1/11 = 0.09) to reach 0.36. Measured PC1-and-friends: top-3 share
   **0.600**, stable rank **2.057**. The fields span about two effective directions, not eleven.

**And the per-target figure is a trap that the control caught.** The ORACLE *per-target* combination
reaches **rho = 0.949** (implied 0.967 A) — which looks spectacular and is meaningless: the matched
random 21-dimensional subspace control reaches **0.809** on the same targets. Fitting 21
coefficients per target to a ~33-dimensional vector gets you most of the way there whatever the
fields contain. **Excess over control: +0.1396.** Without that control this would have gone into
the record as a 0.97 A result.

**So lane D's own reframe — "the fields are real and worth 0.0195 A because sqrt(1-rho^2) squares
them" — survives, and the follow-up I built on it does not.** The fields are real, they are
significant at up to 7.7 sigma, they span two dimensions, and combining them buys nothing. I
proposed the combination as the sprint's central measurement; it took an afternoon and the answer
is no.

**Lane T's prediction landing at 0.169 against a registered point estimate of 0.15, with the count
cancelling as derived, is the strongest methodological result of the sprint so far** — a formula
registered in advance that predicted a measurement to within 0.02 in a quantity nobody had computed.

## CORRECTION TO THE SYNTHESIS BELOW (2026-09-20 13:08, lane F): **THE WIDENING HALF IS CIRCULAR AND THE FILTER-WIDTH HALF IS CLOSED**

Written within the hour of the synthesis it corrects. Both halves of my "two qubits" mechanism took
damage from lane F's own measurements, and I am recording it before the synthesis is quoted.

**1. Filter width is closed by ceiling, not open.** I hypothesised that if the 108's benefit
saturates early in k while the tail's harm grows, a wider filter keeps the benefit and sheds the
harm with no detection. Lane F tested it:
- **F2a refuted — there is no free lunch.** Benefit-retained against harm-shed crosses smoothly at
  ~55%/55% near k ≈ 275. **No knee.**
- **F2b closed by ceiling without spending the chain:** the **ORACLE global argmin over k IS the
  shipped 75**, so no leave-fold-out k arm can beat production. Fifth instance of the global-scalar
  pattern, and lane F pre-registered that it expected exactly this.

What survives from that experiment is a **methodological** result I should have had before
proposing the hypothesis: **filter width moves the emitted structure 0.3790 Å over its range against
0.0256–0.0927 Å for averaging width — 4 to 15×.** My claim that S29's flat m sweep says nothing
about filter width was right; the inference I drew from it was not.

**2. The widening-rescues-the-tail result is CIRCULAR, and lane F flagged it rather than selling
it.** Widening appears to rescue FAIL18 (−0.6193 Å at k = 400, null p = 0) — but **it does not
replicate on either filter-independent tail**, coming in at ≈0 and *reversed* at **+0.7313** on
worst-18-by-pool-mean. One quantity explains all four strata at **ρ = +0.81**: widening helps
exactly where the filter *hurt the set mean*. That is `operator-consumes-set-mean` read backwards,
not a tail rescue.

**3. And the same circularity has to be applied to lane Q's §5, which I put at the centre of the
synthesis.** Lane Q's finding is that on FAIL18 the ORACLE best candidate sits at rank 128–500,
outside the register. **But FAIL18 is defined by production, production is the filtered set's
average, and a bad filter is exactly what pushes good candidates down the ranking.** So "the good
candidate is outside the window" may be *produced by* the thing it is offered as evidence about.
Lane Q's number is on a different readout from lane F's (argmin versus set mean) so the two do not
directly conflict — but it has not been run against a filter-independent tail, and **until it is,
the two-qubit mechanism rests on a statistic with a known circularity.** That check is cheap and it
is now the gate on that whole direction.

**What the synthesis below keeps.** Items 1, 2, 3 and 6 stand — the tail is selection-limited, the
score has no in-pool skill there and that **replicates on all three tail definitions**, the filter
is worse than random on the tail, and the argmin dominates the sparse weighted family by price.
What falls is the *intervention* I built on top of them.

## SYNTHESIS (2026-09-20 13:07): SIX LANES CONVERGED ON ONE PLACE, AND IT IS TWO QUBITS WIDE

Everything closed today points at the same intervention, reached independently. Stating it as a
hypothesis with its falsifier, not as a result.

### THE CHAIN

1. **The tail is selection-limited, not pool-limited** (F). ORACLE best of the FAIL18 pools is
   **2.2842 Å** — already under the 3.00 cap. 13 of 18 have a member under 3.00.
2. **The score has no ranking skill there.** In-pool ρ is **+0.1066 on the tail (fold CI includes
   zero)** against **+0.6446 on the 108** (F).
3. **And its filter is actively worse than random on the tail** — 0W/18L, median per-target
   percentile 0.99999, random-18 null p = 0 (F).
4. **The tail's answer is OUTSIDE the register.** On FAIL18 the ORACLE best candidate sits at rank
   **128–500**, beyond the top-128 window the quantum stage sees (`core/pipeline.py:758`).
   Widening 128 → 512 is worth **−1.9004 Å on FAIL18** and −0.1907 on the 108; difference of
   stratum means −1.7097, SE 0.2207, **2.77× MDE** (Q).
5. **Widening helps a SELECTING readout and hurts an AVERAGING one**, and neither lane's number
   says this alone. Lane F has the ORACLE *set mean* worsening with width on the 108 (3.116 →
   4.201) and improving on FAIL18 (6.159 → 5.755); lane Q has the *argmin* improving sharply with
   width on FAIL18. **The pair is the statement.**
6. **And the argmin is the cheap operator, not the expensive one.** Lane Q priced support and
   weights in the readout's own currency: the argmin over top-2^B **dominates the sparse weighted
   family at every budget** (B=9: 1.7108 vs 2.4375). On the built chain, 2-of-75 with free
   continuous weights gives 2.1683 Å at 11.4+ bits; **argmin over top-128 gives 2.1435 Å at 7.0
   bits and no weight channel.** S29's last unclosed ladder class is closed **by price**, not by
   ceiling — the 1.4315 Å figure was never a route.

### SO THE CANDIDATE MECHANISM IS

> **A selecting readout over a wider register, aimed at the tail.** Two qubits take 128 → 512.
> It is quantum-relevant, it is in the stratum the opening arithmetic says holds the prize, and
> every closed direction today points away from the alternatives.

### WHAT MUST BE TRUE FOR IT TO PAY, AND IT IS THE HARD PART

**All of step 4–5 is ORACLE.** It says the candidate *is there*, not that anything can find it.
The mechanism needs a **native-free selector** that works on hard pools — and lane F measured that
the deployed score has no in-pool skill there at all. So this does not escape the recognition
problem; **it relocates it and prices it.** The register question is cheap and settled; the
selector question is lane R's and lane D's and is open.

Stated as the falsifier: *if no native-free rule can order the top-512 on FAIL18 better than the
top-128, widening buys nothing and the mechanism is dead.*

### WHAT IS CLOSED, AND BY WHAT

| direction | closed by | how |
|---|---|---|
| sparse weighted readout (L6) | Q, S30-L11 | **by price** — argmin dominates at every budget |
| second-moment / quadric escape (L5) | Q **and** T, independently | +0.086 Å (T, 2.95× MDE) and +0.2059 (Q, 1.81×); the structured `disp2` form scores 3.3585 vs production 3.0483 |
| subset objective through an averaging readout | T, S30-L9 | the tail is **always** a prefix — of the order induced by ∇V at the optimum; S29 §4.3 and §4.5 are incompatible |
| generative structural spaces | X, S30-L10 | closed **jointly with the readout**: at ρ≈0 the endpoint is near-unit-slope in the set **mean**, and width buys ceiling while costing mean |
| torsion/configuration encodings | T | arithmetically infeasible — 48 bits needed against 7 deployed |
| common-mode correction from pool data | L, S30-L7 | **non-identifiable**: the likelihood depends on (t, μ) only through t+μ at any K; m_eff = **1.4** of 75 |
| E1 (prior on μ's form) | L, S30-L13 | three independent ways |
| E2 (constraint repair) | L, S30-L13 | field-scale result: 98% of clashes removed costs **+0.08 Å**; works here only because 2/n is 15.4% at n=13 |

### CORRECTIONS TO MY OWN INSTRUCTIONS, BOTH FROM LANE X

- **The gate I handed lane X — "measure the ceiling first, always" — is an ANTI-PREDICTOR.**
  Scored on predicting the direction the endpoint moves: **1 of 10** cells correct, against 5/6 for
  the set-mean law, mean |residual| 0.333 Å vs 0.0153 — **22×**. Enlarging the source improves the
  ORACLE ceiling monotonically by 0.846 Å while the endpoint gets monotonically 0.081 Å **worse**.
  The replacement, and it is usable: **admit G against P iff Δ(set mean) < −(0.32…0.37)·Δ(set
  best)** — at fixed ceiling the mean is worth **2.7–3.1×** what the best is worth.
- **"Candidate generation closed on five instruments" is supported by one.** Four of the five
  measured only the endpoint and carry no oracle field at all. Recomputed: an untrained
  per-residue-type Ramachandran sampler **beats the entire K=500 pool's best member on 45 of 126
  targets** and the endpoint still worsens. The verdict should read *"no achievable source produces
  candidates the shipped score can convert into Ångströms"*, not *"no alternative source contains
  better structures"* — which is false and has been used to rule proposals out.

### THE SCOPE CAVEAT ON A NUMBER I HAVE BEEN QUOTING

Stable rank **1.86 is the PAIR-DISTANCE feature space**. In **coordinate** space it is **3.619**
(k90 = 11.2 directions, not 5.6) — above lane T's own 2.0 threshold. Sparse/weighted readouts and
second-moment constructions act on **coordinates**. I wrote the number without a space attached in
note 4; both values must always carry one.

## WHAT I AM TREATING AS BINDING FROM S29 (until a lane breaks it)

- Every native-free operator is a displacement; its whole value is one cosine. 3.00 Å needs
  ρ = 0.358, 2.50 needs 0.628; everything measured is ≤ 0.04 **signed** — but per-target |cos| is
  0.25–0.33, so **the magnitude exists and only the sign is missing**.
- The deployed architecture's ORACLE ceiling is 2.9027 Å built chain. 2.5 Å is unreachable through
  it *with the native in hand*.
- The terminal operator consumes the set **mean**, so a perfect rank-1 is worth ≈ −0.03 Å through
  the shipped m = 75 average. **A better channel needs a different readout to be worth anything.**
- The deployed score supplies 1.442 of 7 bits against 1.405 for random; the median target is
  *worse* than chance.
- Recognition is closed across bands (theorem), within bands (F2 fails 58/70) and per target
  (incidental parameter) — **but not in the class-B2-is-empty sense**: 9 of 31 channels keep skill
  after partialling out compactness both ways.

## WHAT S29 LEFT GENUINELY OPEN

1. **LEG_torsion and 8 sibling channels** — in-band skill that survives removing monotone *and*
   V-shaped Rg dependence; LEG_torsion is a function of the structure, not the distogram, so it is
   outside the class theorem 2 bounds.
2. **A sparse weighted readout** — the only ladder class not closed by ceiling (2 members with
   oracle weights reach 1.43 Å). Needs a native-free support rule and ~18 bits, not 7.
3. **A free-energy stage** — the only native-free selector class with peptide-length precedent. It
   does not exist on disk; a rebuild, not a resume.

## LANES

Seven of a permitted eight; one slot held for what the results demand.

| lane | remit | status |
|---|---|---|
| **R** | is nativeness recognizable from one structure at all (L11) | running — plan approved, building the matched-kind ladder |
| **F** | the failure tail; highest leverage on the endpoint by the opening arithmetic | running — first task is the pool-limited vs selection-limited diagnosis |
| **Q** | L5 **and** L6 together: a subset objective is pointless through an averaging readout, a sparse readout is unusable without a support rule | running |
| **X** | divergent, permanent: should the quantum stage select at all, or generate? | running |
| **T** | theory: the bit accounting (L8) and when the CVaR tail stops being a prefix | running |
| **D** | adversary, permanent; owns and extends the cost/RMSD meter | running |
| **L** | literature, permanent | running |

## OPEN QUESTIONS I AM HOLDING

- Does the tail's difficulty have a *cause* we can name, or is it just harder sequences?
- Is there a formulation that dissolves several leads at once rather than patching each?
- The launch-gate defect (`jobrun.py` CPU_START 85 vs `governor.py` CPU_CEILING 101) is still open
  and unresolved from S29; RAM, not CPU, is currently the binding constraint so it has not bitten
  yet this sprint.
