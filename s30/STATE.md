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
