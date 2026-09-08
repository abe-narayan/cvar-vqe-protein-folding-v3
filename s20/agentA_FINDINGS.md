# SPRINT 20 — AGENT A (structural inference / RMSD)

**Question owned.** *Can any candidate generator produce structures whose COHERENT error
decorrelates from the retrieval pool's?* Metric: `rho(e_coherent, e_pool)`, not member RMSD.

Pre-registration: `s20/PREREG_A.md`, written before any arm ran and **not edited since**.
Code: `s20/a_src.py` (sources), `a_sel.py` (corpus × selector), `a_gate.py` (selector identity +
partial correlation), `a_proj.py` (built basis), `a_report.py`.
Artefacts: `s20/results/a_{src,sel,gate,fuseceil,proj}.json`, **all five COMPLETE at n = 126**,
plus the rendered `a_report.json`. One superseded artefact is retained unread and unquoted as
`_SUPERSEDED_a_proj_norho.json`.

**Instrument.** 126 cluster-disjoint tuning targets, pinned folds, `s12/instrument.py`, full-chain
Cα-RMSD, frozen implementation. Sealed 60-target benchmark **untouched**;
`results/benchmark_manifest.json` never opened.

**Reproduction gates.** Two, both passed before anything was interpreted.

1. My reconstruction of the operator — `pool_idx` → `pair_dists` → `shipped_score` →
   `argsort[:75]` — is **set-identical to `shipped_record(pdb)["sub"]` on 126/126 targets**, and
   the coordinate average of that set reproduces 3.048 Å exactly.
2. My `pool` arm projected onto the ideal-geometry manifold reproduces the **shipped production
   structure** `shipped_record(pdb)["fit_ca"]` to **max |Δ| = 0.0004 Å** (a per-target identity,
   valid at any n). This independently confirms the coordinator's 2026-09-07 correction that
   projecting the coordinate average *gives* the incumbent, and it means the built basis in §7 is
   the incumbent's own basis and not a re-implementation of it.

I also reproduce the coordinate average's contraction independently at n = 126: mean virtual
Cα–Cα bond **2.961 Å against 3.804 Å physical — 22.2%**, matching the coordinator's figure exactly,
while the candidate **members** the average consumes are all physical (3.804–3.814 Å).

Every arm below runs the same operator with only the candidate index set changed.

---

# 0. LEAD WITH WHAT DAMAGES MY OWN HYPOTHESIS MOST

> ## No generator on my board decorrelates enough to matter, and the one stage that does decorrelate is not the generator. **A candidate set with NO sequence retrieval at all — 500 windows drawn uniformly at random from the library — emits a structure 0.712 Å from the incumbent's, with an error aligned 0.957 (0.807 after partialling out the zero-information mode), and lands at 3.157 Å against the incumbent's 3.048 Å.**

Swapping the **corpus** moves the emitted error alignment from 1.000 to 0.951–0.968.
Swapping the **selector** moves it to 0.794.
Removing the target's sequence from retrieval entirely moves it to **0.957**.

*(That headline is the outcome statement. The mechanistic half needs the qualification in §0.1,
which I pre-committed to in §7.1 before the deciding number existed: the corpus **does** carry a
measurable minority of the shared direction. It is worth nothing, which is why the outcome
statement stands and the mechanistic one is narrowed rather than kept.)*

**F-A1a FIRES** — `rho(e_pep, e_pool) = 0.968`, against a pre-registered bar of 0.85.
**F-A1b FIRES** — the corpus partition reaches **0.963** of the matched random-partition ceiling,
against a bar of 0.85.

**Candidate sources 1 and 2 of my brief (peptide-first corpus; independent structural corpus) are
CLOSED on outcome**, and source 3 (classical continuous-torsion generator) with them at 0.898.

**On the incumbent's own basis, nothing on my board beats it.** The best arm of the 13 built is
`fuse(pool, pep)` at **3.182 Å against 3.205 Å — −0.023 Å, a quarter of the MDE, 63W/63L, a dead-even
coin flip.** The peptide corpus alone is −0.017 [−0.096, +0.056], 62W/64L, median +0.006.

## 0.1 AND THE PART THAT DAMAGES MY OWN CLOSURE, PRE-COMMITTED IN §7.1

The sentence "the corpus is not the carrier" is **too strong as a mechanistic claim, and I said so
before I saw the number that shows it.** On a stricter post-hoc statistic — the alignment with the
zero-information "typical peptide of this length" mode partialled out — the corpus partition
decorrelates **materially** more than a matched random partition: **0.679 vs 0.870, share 0.780
against the 0.85 bar, −0.1917 [−0.2313, −0.1564], 115W/11L, 5/5 folds.**

**The corpus carries a real minority of the shared error direction.** What it does not do is
convert: it is worth −0.017 Å on the incumbent's basis with a coin-flip win rate, and §4.1b shows
that at target level the fusion gain is explained by the partner's *quality* at Spearman 0.96–0.99
and by its decorrelation at |ρ| ≤ 0.29 with an inconsistent sign.

> **The correct statement is: a decorrelating candidate source exists, the corpus is part of what
> makes it decorrelate, the selector is the larger part — and it is worth nothing.**

---

# 1. THE SOURCE TABLE

`s20/a_src.py` · `results/a_src.json` (COMPLETE, n = 126) · **BASIS: POINT CLOUD**
(`I.coordinate_average`; reference `pool` = 3.048 Å). Every source runs the identical operator:
500 candidates → shipped Bayes-risk distogram score → lowest 75 → coordinate average.
`rho`, `member`, ceilings and coverage are **ORACLE** diagnostics.

| source | fitted to | realised | gen ceil | sel ceil | member | cov<2Å | divGEO | divERR | bond | ρ\|pool | memρ\|pool |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `pool` (incumbent) | target seq (BLOSUM) + library | **3.048** | 1.711 | 2.306 | 3.551 | 0.229 | 2.486 | 0.303 | 2.961 | — | 0.697 |
| `pep` peptide corpus | target seq; **peptide DB only** | 3.018 | 1.671 | 2.173 | 3.545 | 0.235 | 2.618 | 0.329 | 2.925 | **0.968** | 0.673 |
| `prot` protein fragments | target seq; **fragments only** | 3.127 | 1.857 | 2.380 | 3.600 | 0.222 | 2.455 | 0.286 | 2.960 | 0.986 | 0.704 |
| `halfA` random half | target seq; random half of library | 3.071 | 1.684 | 2.302 | 3.574 | 0.227 | 2.527 | 0.304 | 2.934 | 0.990 | — |
| `halfB` random half | target seq; the complement | 3.066 | 1.685 | 2.277 | 3.567 | 0.229 | 2.532 | 0.305 | 2.933 | 0.990 | — |
| `rand500` | **nothing** (500 uniform windows) | 3.142 | 1.820 | 2.276 | 3.667 | 0.217 | 2.685 | 0.318 | 2.853 | 0.954 | 0.675 |
| `tors` Ramachandran gen. | **nothing** — literature constants | 3.480 | 2.046 | 2.426 | 3.963 | 0.171 | 2.910 | 0.330 | 2.858 | 0.898 | 0.630 |
| `unsel` 75 random, no score | **nothing**, no selector at all | 3.768 | 1.313 | 2.245 | 4.816 | 0.061 | 4.505 | 0.500 | 2.175 | **0.726** | 0.475 |
| `helix` zero-information | **nothing** | 4.065 | — | — | — | — | — | — | 3.804 | 0.814 | — |
| native | — | 0.000 | — | — | — | — | — | — | 3.812 | — | — |

**Set widths, because the operator law is width-dependent** (coordinator, 2026-09-07): every
retrieval arm is *m*=75 selected from *K*=500, the narrow regime where the recorded law puts 0.04
on `set_best`. `unsel`'s generation set is the **whole universe** (mean 18,674 windows), which is
why its generation ceiling (1.313 Å) is the best on the board and its realised number the worst.

**`tors` is fitted to nothing.** Its parameters (`_MU`, `_SG`, `_TRANS`, `_INIT` in `a_src.py`) are
hardcoded literature backbone-basin constants: a 3-state H/E/C first-order Markov chain with
per-state Gaussian (φ,ψ). It never sees the pool's marginals. This is the clean case the
coordinator asked for, and it aligns at **0.898**.

---

# 2. THE PRE-REGISTERED FALSIFIER

`rho` is the mean-centred Pearson correlation of `e = d(X) − d_true` across a target's pairs,
averaged over the 126 targets — the same semantics as `s19/a_source._corr`.
`X` is a real structure, so `e` **is** its coherent error and the Sprint-19 split `r = r_coh +
r_inc` is degenerate with `r_inc = 0`. **That is an identity, not a finding.**

| quantity | value | bar |
|---|---|---|
| `rho(e_pep, e_pool)` | **0.9684** (median 0.9865) | F-A1a: ≥ 0.85 → **FIRES** |
| `rho(e_pep, e_prot)` — corpus partition | 0.9510 (median 0.9758) | |
| `rho(e_halfA, e_halfB)` — **random partition, the matched ceiling** | 0.9881 (median 0.9940) | |
| share of ceiling | **0.9625** | F-A1b: ≥ 0.85 → **FIRES** |
| corpus − random partition | **−0.0371 [−0.0484, −0.0272]**, 118W/8L, **5/5 folds** | |
| long-range only (sep ≥ 5) | corpus 0.9497, random 0.9883 | same picture |

**The corpus effect is real and it is tiny.** Partitioning the universe by corpus decorrelates
significantly more than partitioning it at random — the CI excludes zero and every fold agrees —
but it moves the alignment from 0.988 to 0.951. It is a 3.7% effect where the falsifier needed a
15% one.

**At the population level the same statement is sharper.** Mean member-error alignment between the
two **disjoint** corpora is **0.673**, against **0.697** *inside the incumbent pool itself* and
0.695 between two random halves. **Two disjoint corpora are wrong in the same way that one corpus
is internally wrong.**

### 2.1 Stratified by length and by difficulty — the falsifier fires in every stratum

| stratum | n | ρ(pep,pool) | ρ(pep,prot) | ρ(halfA,halfB) | share of ceiling | `pep`−`pool` |
|---|---|---|---|---|---|---|
| ALL | 126 | 0.9684 | 0.9510 | 0.9881 | **0.9625** | −0.031 |
| len ≤ 10 | 20 | 0.9455 | 0.9269 | 0.9900 | 0.9363 | +0.113 |
| len 11–12 | 32 | 0.9658 | 0.9451 | 0.9897 | 0.9550 | −0.060 |
| len 13–14 | 37 | 0.9710 | 0.9572 | 0.9873 | 0.9695 | −0.122 |
| len ≥ 15 | 37 | 0.9804 | 0.9629 | 0.9864 | 0.9761 | +0.008 |
| **FAIL18** | 18 | **0.9820** | 0.9787 | 0.9963 | 0.9824 | −0.209 |
| not FAIL18 | 108 | 0.9661 | 0.9464 | 0.9867 | 0.9591 | −0.001 |

**Every stratum clears the 0.85 bar by a wide margin**, and the alignment is *highest* on the 18
hardest targets — the ones with zero in-band recall — which is what an identifiability wall should
look like.

**And the peptide corpus's entire apparent gain lives on FAIL18 and is not measured even there:**

    ALL         n=126   pep 3.018  pool 3.048   -0.031 [-0.104,+0.036]  68W/58L
    FAIL18      n= 18   pep 5.623  pool 5.832   -0.209 [-0.559,+0.039]  12W/6L
    not FAIL18  n=108   pep 2.584  pool 2.584   -0.001 [-0.062,+0.059]  56W/52L

On the 108 targets that are not FAIL18 the peptide corpus is worth **exactly nothing** (−0.001).
Project memory records that sequence conditioning *hurts* on FAIL18 and that the peptide corpus is
the one carrying the sequence–structure channel; this is the targeted check, and at n = 18 the CI
spans zero. **NOT MEASURED, and flagged as a subgroup that would need its own instrument.**

**Structurally, they are barely different objects.** Each emitted structure is ~3 Å from the
native; the emitted structures are **0.34–0.83 Å from each other**:

    pool ~ prot   0.344      halfA ~ halfB  0.385      pool ~ pep   0.661
    pool ~ halfA  0.363      pool ~ rand500 0.717      pep ~ prot   0.829
    pool ~ tors   1.419      pool ~ unsel   2.498      pool ~ helix 2.554

---

# 3. IT IS THE SELECTOR, NOT THE CORPUS — and this was NOT pre-registered

`s20/a_sel.py`, `a_gate.py` · `results/a_sel.json`, `a_gate.json` (both COMPLETE, n = 126) ·
**BASIS: POINT CLOUD** · **POST-HOC decomposition, reported as such.**

The coordinator's constraint — *any Q2 test whose candidate space is fitted to the retrieval pool
is asking a foregone question* — exposed a confound I had not pre-registered: **every arm in §1
swaps the corpus but keeps the selector**, and the selector is the shipped leave-fold-out
distogram, trained on the same library. So a null in §1 cannot separate "the corpus is not the
carrier" from "the corpus is not the carrier *because the selector is*". I crossed them 4 × 2.

| contrast | ρ | ρ partialled on the zero-info helix mode | RMSD(X_A, X_B) |
|---|---|---|---|
| **CORPUS swap**, selector held (`blosum.disto ~ pep.disto`) | 0.968 | — | 0.661 |
| **CORPUS partition**, selector held (`pep.disto ~ prot.disto`) | 0.951 | — | 0.829 |
| **corpus → uninformative**, selector held (`blosum.disto ~ univ.disto`) | **0.957** | **0.807** | **0.712** |
| **SELECTOR swap**, corpus held (`blosum.disto ~ blosum.rand`) | **0.794** | **0.396** | 1.864 |
| **SELECTOR swap**, corpus held (`pep.disto ~ pep.rand`) | 0.810 | — | 1.795 |
| both removed (`blosum.disto ~ univ.rand`) | 0.726 | 0.284 | 2.524 |
| corpus swap with **no selector anywhere** (`blosum.rand ~ univ.rand`) | 0.935 | 0.813 | 1.306 |

> **The distogram selector, applied to a candidate set that carries no information about the target
> at all, reproduces the incumbent's error direction at 0.957 raw and 0.807 after the generic
> peptide mode is partialled out. The pool's own corpus without that selector reproduces it at
> 0.794 raw and 0.396 partialled.**

**And it is not "any score".** A zero-information constant-α-helix gate over the same 500
candidates aligns at 0.819 raw but only **0.268 partialled**, and emits 4.040 Å (+0.991 [+0.710,
+1.294] vs the incumbent). So the shared error *direction* is carried by the **distogram
specifically**, not by score-ordering in general.

> **This does NOT contradict Sprint 19 C4/Z7 and must not be quoted against it.** C4/Z7 measures
> *gate damage to the emitted mean* — how much a score-ordered selection costs in RMSD — and finds
> the helix gate the worst arm. I measure something else: *which direction the resulting error
> points*. The helix gate is simultaneously the most damaging gate (+0.991) and the least
> pool-aligned one (partial 0.268). Both are true, they are different statistics, and the
> resolution is that the helix gate is wrong in its **own** direction rather than the pool's. The
> two arms are also not magnitude-matched, so no ranking of them on a common axis is claimed.

**Realised, paired vs the incumbent (point cloud):**

    pep.disto    -0.031 [-0.104,+0.036]  68W/58L     NOT MEASURED (below MDE, CI spans zero)
    prot.disto   +0.079 [+0.047,+0.113]  50W/76L
    univ.disto   +0.109 [+0.046,+0.175]  53W/73L     <- NO sequence retrieval at all
    blosum.rand  +0.387 [+0.221,+0.560]  40W/86L
    pep.rand     +0.416 [+0.219,+0.618]  48W/78L
    univ.rand    +0.779 [+0.556,+1.014]  38W/88L
    *.helixgate  +0.99 to +1.00

**BLOSUM retrieval of the target's own sequence is worth 0.109 Å** through this operator. The
distogram selector is worth 0.387 Å. Neither buys a different error *direction*.

---

# 4. FUSION — the payoff test, and it fails on every arm

**BASIS: POINT CLOUD**, paired against the `pool` arm, n = 126.

| arm | realised | vs pool | W/L | bond |
|---|---|---|---|---|
| `fuse(pool, pep)` | 3.019 | **−0.030 [−0.067, +0.005]** | 73W/53L | 2.912 |
| `merge(pool, pep)` (150-member set) | 3.023 | −0.026 [−0.064, +0.009] | 69W/57L | 2.923 |
| `fuse(halfA, halfB)` — **matched control** | 3.062 | +0.014 [−0.010, +0.038] | 57W/69L | 2.919 |
| `fuse(pool, prot)` | 3.082 | +0.033 [+0.017, +0.050] | 53W/73L | 2.949 |
| `fuse(pool, rand500)` | 3.074 | +0.026 [−0.006, +0.059] | 59W/67L | 2.868 |
| `fuse(pool, tors)` | 3.197 | +0.148 [+0.073, +0.228] | 50W/76L | 2.810 |
| `fuse(pool, unsel)` — the most decorrelated | 3.259 | +0.210 [+0.098, +0.326] | 46W/80L | 2.381 |

**The best arm on the board is −0.030 Å: below the 0.084 Å MDE, with a CI spanning zero. NOT
MEASURED, and never "an improvement".**

**F-A3 FIRES.** Difference of differences, `[fuse(pool,pep) − pool] − [fuse(halfA,halfB) − halfA]`
= **−0.0212 [−0.0597, +0.0153]**, 69W/57L, folds split 2/3. Any fusion gain is *"averaging two
structures"*, not *"averaging two decorrelated structures"*.

**F-A2's antecedent was never satisfied and I will not pretend otherwise.** I pre-registered
F-A2 at `rho ≤ 0.6 × ceiling` = 0.593. The most decorrelated source reached 0.726
(0.735 × ceiling). **The regime "strong decorrelation, measured for conversion" is therefore
UNTESTED**, and what I measured instead is the regime 0.73–0.99, where conversion is nil.

### 4.1 The variance prediction fails, and the way it fails is the finding

If two structures had equal error magnitude and correlation ρ, their average would carry
`sqrt((1+ρ)/2)` of the error. That predicts the gain for the *aligned* arms and inverts for the
decorrelated ones:

| fusion partner | ρ | **measured** gain | predicted gain |
|---|---|---|---|
| `pep` | 0.968 | −0.030 | −0.024 |
| `rand500` | 0.954 | +0.026 | −0.035 |
| `tors` | 0.898 | +0.148 | −0.079 |
| `unsel` | 0.726 | **+0.210** | **−0.216** |

**The decorrelated partner is worse by exactly as much as the equal-magnitude model says it should
be better.** The missing term is magnitude: `unsel`'s member error RMS is 4.681 against the pool's
2.986. **Decorrelation in this system is only purchasable by degrading quality, and the exchange
rate is unfavourable at every point measured.**

### 4.1b It dissolves under target-level analysis, and what replaces it is quality

The mechanism predicts that on targets where the partner is *more decorrelated*, fusion should gain
*more*. Per-target Spearman of the fusion gain against (a) the decorrelation and (b) the partner's
own quality on that target:

| partner | mean ρ | mean gain | ρ(gain, **decorrelation**) | ρ(gain, **partner quality**) |
|---|---|---|---|---|
| `pep` | 0.968 | −0.030 | **+0.172** *(wrong sign)* | **+0.976** |
| `prot` | 0.986 | +0.033 | −0.289 | **+0.983** |
| `rand500` | 0.954 | +0.026 | −0.095 | **+0.990** |
| `tors` | 0.898 | +0.148 | −0.281 | **+0.961** |
| `unsel` | 0.726 | +0.210 | −0.113 | **+0.963** |

> **Fusion gain is explained at Spearman 0.96–0.99 by "was this partner better or worse than the
> pool on this target", and at |ρ| ≤ 0.29 with an inconsistent sign by "was this partner
> decorrelated".** The decorrelation channel does not survive target-level analysis at all. This is
> the cleanest refutation in my lane and it does not depend on any basis, any ceiling, or any
> ORACLE weight.

### 4.2 The ORACLE ceiling closes it with a ceiling, not a null

`results/a_fuseceil.json` (COMPLETE, n = 126). Per target, the **best of 21 fusion weights chosen
with the native in hand**.

**EXACT, stated before the numbers**: this column *cannot* exceed `pool`, because w = 0 is in the
grid. Its W/L is xxW/**0**L by construction. It must therefore be read against its
zero-information and same-corpus partners, never against `pool`.

    ORACLE(pep)   - ORACLE(helix)  = +0.0226 [-0.0387,+0.1043]  57W/36L   NOT MEASURED
    ORACLE(tors)  - ORACLE(helix)  = +0.0293 [-0.0282,+0.0964]  44W/30L   NOT MEASURED
    ORACLE(pep)   - ORACLE(halfB)  = -0.0820 [-0.1377,-0.0364]  65W/37L   ~ at the MDE
    ORACLE(unsel) - ORACLE(helix)  = -0.0944 [-0.1653,-0.0258]  53W/27L

> **With a per-target ORACLE weight, the peptide corpus is statistically indistinguishable from a
> constant α-helix as a fusion partner, and so is the generic Ramachandran generator.** The only
> source that beats the zero-information partner is `unsel` — whose fused point cloud has a mean
> virtual bond of **2.650 Å**, i.e. it wins by contracting.

Equal-weight fusion of **all eight** sources: 3.185 Å, **+0.136 [+0.052, +0.229]**, 55W/71L, bond
2.682. More sources is worse.

---

# 5. DIVERSITY, SPLIT AS THE BRIEF REQUIRES

Geometric diversity = mean pairwise Cα-RMSD inside the 75-member set. Error diversity =
`1 − mean pairwise Pearson correlation of the members' own error fields` (ORACLE).

| source | divGEO | divERR | member-error alignment inside the set | member err RMS |
|---|---|---|---|---|
| `pool` | 2.486 | 0.303 | 0.697 | 2.986 |
| `pep` | 2.618 | 0.329 | 0.671 | 2.967 |
| `prot` | 2.455 | 0.286 | 0.714 | 3.032 |
| `tors` | 2.910 | 0.330 | 0.670 | 3.355 |
| `unsel` | 4.505 | 0.500 | 0.500 | 4.681 |

**Across the eight sources the two axes are almost perfectly coupled: Pearson +0.992** (n = 8
sources, descriptive, not inferential). Within a source across targets they decouple only for the
two non-retrieval sources (`tors` −0.08, `unsel` +0.05).

> **The brief's worry — "two populations can be geometrically different and make the same
> structural mistake" — is measured and is the *general* case here, not a corner case.** `pep` and
> `prot` are disjoint corpora of near-equal geometric diversity (2.618 vs 2.455) whose members'
> errors align at **0.673**, the same as the pool's internal 0.697. **No source in this study
> bought error diversity without buying geometric diversity, and every one that bought geometric
> diversity paid for it in quality.**

---

# 5b. THE CEILING LADDER, and where the in-band best is actually lost

The brief requires generation / selection / repair / realised kept separate. Mean over 126, ORACLE:

    whole universe (mean 18,674 windows)          1.313
    BLOSUM top-500  (generation ceiling)          1.711
    random 75 drawn from the whole universe       2.245
    random 75 drawn from the BLOSUM-500           2.061
    SHIPPED top-75, distogram  (selection ceiling)2.306
    realised, point cloud                         3.048
    realised, built chain                         3.204   (the incumbent)

**Almost the whole in-band ceiling is spent between the 500-pool and the 75-set** (1.711 → 2.306),
and the distogram filter is not what recovers it. At **matched N** (75 chosen from the same 500):

| corpus | best\|distogram | best\|random 75 | distogram − random | median | W/L |
|---|---|---|---|---|---|
| `blosum` | 2.306 | **2.061** | **+0.246 [+0.098, +0.398]** | **+0.000** | 50W/62L |
| `pep` | 2.173 | 2.115 | +0.058 [−0.077, +0.193] | +0.000 | 61W/56L |
| `prot` | 2.380 | 2.170 | +0.210 [+0.073, +0.358] | −0.008 | 64W/52L |
| `univ` | 2.288 | 2.200 | +0.087 [−0.075, +0.253] | −0.061 | 75W/43L |

> **Read the mean and the median together, or this result is misreported.** The median is
> **exactly zero** and the win/loss count favours the distogram on three corpora of four.

The distribution of `d = best(distogram 75) − best(random 75)` on `blosum`, since the mean alone
would misdescribe it (memory: `median-vs-mean-is-the-free-warning`; a drop-top threshold is not a
valid test on its own, so the decomposition is given in full):

    quantiles 0/5/10/25/50/75/90/95/100 :  -2.77  -0.55  -0.33  -0.09  +0.00  +0.48  +1.28  +1.94  +4.71
    mean +0.246   median +0.000   sd 0.858   14/126 targets exactly 0
    d < -0.5 : 29 targets are NOT it -- only  9 targets, contributing -0.084 to the mean
    d > +0.5 :                                29 targets, contributing +0.318 to the mean
    |d| <= 0.5:                              contributing +0.012 to the mean

**Against the uniform-effect null** — which would put every quantile at +0.246 and no target at
exactly 0 — the observed shape is decisively non-uniform: 14 exact ties, a heavy tail in *both*
directions, and essentially the whole mean carried by **29 targets (23%) where the filter discards
a much better member**, partly offset by 9 where it finds one.

> **The honest statement: the distogram filter is neutral on the in-band best for the typical
> target, finds a much better member on 7% of targets, and discards one on 23%** — while improving
> the set **mean** by 0.90 Å (3.551 vs 4.454).

That trade is correct for the incumbent terminal operator, which consumes the set mean far more
heavily than the set best. It is the wrong trade for any operator that consumes the best, and it is
where a re-priced `set_best` coefficient would bite. **Set widths: every row here is m=75 from
K=500**, the narrow regime.

---

# 6. BASIS — a correction to my own pre-registration

**DEVIATION, recorded rather than edited (brief §7).** `PREREG_A.md` §0 states that the raw
coordinate average is "the programme's best built object, 3.048 Å" and makes it the primary emitted
object. **That was wrong**, and the coordinator corrected it mid-run:
`I.coordinate_average` returns a **point cloud** whose mean virtual Cα–Cα bond is 2.961 Å against
the physical 3.804 Å — a 22.2% contraction. The pre-registration is left unedited.

**What survives and what does not.** All §1–§5 comparisons are **within one basis** — every arm is
a coordinate average of 75 candidates — so the *cross-arm* contrasts, the falsifier and the
diversity statistics are unaffected. What was wrong is the *reference*: 3.048 Å is the best point
cloud, not a structure any arm should be compared to across bases.

**The contraction is a function of the population's diversity, and it is the mechanism behind the
apparent gain of the decorrelated arms:**

    source     divGEO   emitted bond   min bond
    pool        2.486        2.961       2.260
    tors        2.910        2.858       2.092
    unsel       4.505        2.175       1.364      <- a 43% contraction
    helix          --        3.804       3.804
    native         --        3.812

Candidate **members** are all physical (mean member bond 3.808–3.818 Å for every retrieval source,
3.804 for `tors`). **The contraction is created entirely by the averaging operator**, and a
decorrelated population contracts harder. `unsel` at ρ = 0.726 buys its decorrelation and its
apparent point-cloud RMSD partly by collapsing the chain.

**§7 puts every arm on the built basis, where this cannot happen.**

---

# 7. BUILT-CHAIN BASIS — the decisive one. **COMPLETE, n = 126**

`s20/a_proj.py` · `results/a_proj.json` (COMPLETE, n = 126). Every arm's point cloud projected onto
the pipeline's own ideal-geometry manifold (`I.project(...)["fit_ca"]`, λ = 0) — the **incumbent's
own basis**. My `pool` arm lands at **3.205 Å against the shipped 3.2041** and reproduces the
shipped structure per target to max |Δ| = 0.0004 Å. **Every built chain has mean virtual bond
exactly 3.804 Å**, so nothing here can be a contraction artefact.

| arm | cloud | **BUILT** | repair tax | vs `pool` (built) | W/L | folds |
|---|---|---|---|---|---|---|
| `pool` (incumbent) | 3.048 | **3.205** | +0.157 | — | | |
| `pep` | 3.018 | **3.188** | +0.170 | **−0.017 [−0.096, +0.056]** | 62W/64L | 3/5 |
| `prot` | 3.127 | 3.302 | +0.174 | +0.097 [+0.052, +0.144] | 50W/76L | 5/5 |
| `halfA` | 3.071 | 3.238 | +0.167 | +0.033 [−0.006, +0.075] | 57W/69L | 4/5 |
| `halfB` | 3.066 | 3.235 | +0.169 | +0.030 [−0.008, +0.070] | 67W/59L | 3/5 |
| `rand500` | 3.142 | 3.332 | +0.190 | +0.127 [+0.050, +0.210] | 52W/74L | |
| `tors` | 3.480 | 3.676 | +0.196 | +0.471 [+0.301, +0.659] | 36W/90L | |
| `unsel` | 3.768 | 3.965 | +0.197 | +0.760 [+0.512, +1.006] | 35W/91L | |
| `fuse(pool,pep)` | 3.019 | **3.182** | +0.163 | **−0.023 [−0.064, +0.013]** | 63W/63L | 4/5 |
| `merge(pool,pep)` | 3.023 | 3.182 | +0.160 | −0.023 [−0.066, +0.015] | 64W/62L | 3/5 |
| `fuse(halfA,halfB)` | 3.062 | 3.224 | +0.162 | +0.019 [−0.010, +0.048] | 67W/59L | |
| `fuse(pool,tors)` | 3.197 | 3.384 | +0.187 | +0.178 [+0.091, +0.270] | 44W/82L | |
| `fuse(pool,unsel)` | 3.259 | 3.444 | +0.186 | +0.239 [+0.121, +0.360] | 44W/82L | |

> **Nothing beats the incumbent.** The best arm on the whole board is `fuse(pool,pep)` at
> **−0.023 Å**, a quarter of the MDE, CI spanning zero, **63W/63L — a dead-even coin flip**, and its
> median is +0.001. Under §9 of my pre-registration it fails all five promotion criteria.

**The repair tax rises monotonically with population diversity**: +0.157 for the incumbent (matching
the coordinator's 0.156) to +0.197 for the most diverse source. **A decorrelated population pays a
larger stereochemical-realisability tax**, which is the second axis the coordinator asked me to
measure, and it moves *against* decorrelation.

**F-A3 fires again on this basis**: difference of differences −0.0094 [−0.0536, +0.0320], 63W/63L.

## 7.1 THE COMMITMENT I MADE BEFORE SEEING THIS, AND IT LANDED AGAINST ME

I recorded in advance (unedited above this line in the file's history) that the helix-partialled
statistic was pointing against my closure at partial *n*, and that I would report it either way.
**It lands against me, and here it is.**

| statistic | pep ~ pool | corpus partition | random partition (ceiling) | share | bar |
|---|---|---|---|---|---|
| **RAW `rho`** — *pre-registered* | **0.923** | 0.893 | 0.958 | **0.9320** | 0.85 → **FIRES** |
| **helix-partialled** — *post-hoc* | 0.756 | 0.679 | 0.870 | **0.7798** | 0.85 → **does NOT fire** |

`corpus − random`, partialled: **−0.1917 [−0.2313, −0.1564], 115W/11L, 5/5 folds.**

**So the honest verdict is two-part, and I will not report only the half that suits me:**

1. **On the pre-registered statistic, F-A1a and F-A1b FIRE on both bases.** `PREREG_A.md` §6
   defines them on raw `rho`, which is also the scale Sprint 19's 0.85 bar was expressed in.
2. **On a stricter post-hoc statistic that removes the generic "typical peptide of this length"
   mode, they do not.** The corpus partition decorrelates materially more than a random partition
   of the same universe — 0.679 vs 0.870, 5/5 folds, a CI nowhere near zero. **The corpus carries a
   real minority of the shared error direction.** My §0 headline is therefore too strong as a
   *mechanistic* claim and is corrected in §0.1.

**What is unchanged is the part that decides the mission.** The conversion evidence never used
`rho` at all: `pep` is **−0.017 [−0.096, +0.056], 62W/64L, median +0.006** on the incumbent's own
basis, every fusion arm fails, F-A3 fires on both bases, and §4.1b shows the fusion gain is
explained by partner *quality* at Spearman 0.96–0.99 and by decorrelation at |ρ| ≤ 0.29 with an
inconsistent sign. **A decorrelating source exists; it is worth nothing.**

An earlier revision of this pass that lacked the built-basis `rho` was killed at 10/126 and is
retained as `_SUPERSEDED_a_proj_norho.json`, unread and unquoted.

---

# 7b. MY OWN DEFECTS, VOLUNTEERED

**1. A `.COMPLETE` flag that fired on the wrong condition.** `a_proj.run(tg)` wrote its completion
flag whenever `len(rows) == len(tg)` — the subset it was *called with*. A 2-target smoke test
therefore left a valid-looking `a_proj.COMPLETE` on disk, and a later wait-loop keyed on that flag
**exited immediately on a run that was 8/126 done**, reporting "finished". Nothing was
mis-published (the JSON's own `complete: false` caught it on the next line), but the flag would have
certified a partial as complete to any reader. Fixed: the flag now requires the full 126-target
instrument. **This is brief Z6 in a new place — not a gate that never fired, but a gate that fired
on the wrong condition, which is worse because it looks like a pass.**

**2. A basis error in my own pre-registration**, corrected by the coordinator and recorded in §6
rather than edited away.

**3. A confound I did not pre-register** — corpus and selector were never crossed in
`PREREG_A.md` — surfaced by the coordinator's constraint on pool-fitted candidate spaces, and
answered in §3, labelled post-hoc throughout.

**Index-construction self-audit passed** on the risky operations: `pep`/`prot` disjoint and
correctly typed, `halfA`/`halfB` disjoint, BLOSUM rank order preserved under every subset filter,
all cached top-75 sets contained in their parent pools, `stable_rng` reproducible across processes.

---

# 8. THE LEDGER

| # | branch | mechanism | result | status |
|---|---|---|---|---|
| A1 | **peptide-first corpus** (brief candidate 1) | change the corpus, not its consumption | raw `rho(e_pep,e_pool)` **0.968** cloud / **0.923** built, 0.93–0.96 of the matched ceiling; **BUILT realised −0.017 [−0.096,+0.056], 62W/64L, median +0.006** | **CLOSED ON OUTCOME — F-A1a/b FIRE on the pre-registered raw statistic; see A10** |
| A2 | **independent structural corpus** (candidate 2) | protein-fragment-only, disjoint from A1 | ρ = 0.986 cloud / 0.956 built; **BUILT realised +0.097 [+0.052,+0.144]**, 5/5 folds | **CLOSED** |
| A3 | **classical continuous-torsion generator** (candidate 3) | sample, don't retrieve; **fitted to nothing** | ρ = 0.898 cloud / 0.830 built; **BUILT realised +0.471 [+0.301,+0.659]**; ORACLE-fused indistinguishable from a constant helix | **CLOSED** |
| A4 | learned generative model (candidate 4) | — | **NOT RUN.** Gated on A1–A3 showing the effect is real and corpus-carried. It is neither. | **CLOSED by gate** |
| A5 | fusion of decorrelated sources | coherent-error cancellation | best −0.030 [−0.067,+0.005]; F-A3 fires; ORACLE ceiling ties a zero-information helix | **REFUTED** |
| A6 | **corpus vs selector** (post-hoc) | which stage imposes the shared direction | corpus swap 0.951–0.968; selector swap 0.794; uninformative corpus + distogram **0.957 (0.807 partialled)** | **ESTABLISHED — the carrier is the selector** |
| A7 | is it *any* score? (post-hoc) | zero-information helix gate | raw 0.819 but **partialled 0.268**, realised +0.991 | **the distogram specifically, not score-ordering in general** |
| A5b | fusion on the INCUMBENT'S basis | coherent-error cancellation, built chains | best arm `fuse(pool,pep)` **−0.023 [−0.064,+0.013], 63W/63L**, median +0.001; F-A3 fires again | **REFUTED** |
| A8 | strong decorrelation (ρ ≤ 0.6 × ceiling) | — | no source reached it on the raw statistic; lowest 0.726 cloud / 0.662 built | **NOT MEASURED / OPEN** |
| A9 | the repair (stereochemical-realisability) tax | project onto the ideal-geometry manifold | **+0.157 for the incumbent, rising monotonically to +0.197 for the most diverse source** | **ESTABLISHED — the tax moves AGAINST decorrelation** |
| A10 | **does the corpus carry ANY of the direction?** (post-hoc) | raw ρ vs ρ partialled on the zero-information mode | raw share 0.932 (fires), **partialled share 0.780 (does not)**; corpus−random partialled **−0.1917 [−0.2313,−0.1564], 115W/11L, 5/5 folds** | **SUPPORTED — the corpus carries a real MINORITY, and it does not convert** |

## What I would tell the next sprint

0. **The mission number did not move.** Incumbent 3.204 A (my reproduction 3.205); my best built arm
   3.182 A at −0.023 [−0.064, +0.013], 63W/63L. Nothing here is promotable and nothing here
   justifies unsealing the benchmark.
1. **Stop looking for a decorrelated generator.** Four sources spanning a full sequence-retrieval
   corpus, two disjoint sub-corpora, a random half, an unfitted Ramachandran sampler and an
   unselected random draw all emit structures 0.34–1.42 Å apart with errors aligned 0.90–0.99.
   The generator is not where the shared error lives.
2. **The selector is.** A candidate set with zero target information, scored by the distogram,
   reproduces the incumbent's error at 0.807 *after* the generic peptide mode is removed. If
   anyone reopens Q2, the object to attack is the **distogram-as-selector**, not the pool.
3. **Report the partial correlation, not the raw one.** A constant α-helix aligns with the
   incumbent at 0.814 raw. Any raw alignment below ~0.85 in this system is mostly the statement
   "both are peptides of this length". The helix-partialled figure separated the distogram
   (0.807) from a zero-information gate (0.268) and nothing else did.
4. **Check the virtual bond on anything that looks unexpectedly good.** The most decorrelated
   source in this study wins on the point-cloud metric partly by contracting the chain 43%.
5. **Decorrelation and quality did not dissociate, on any axis measured.** Geometric and error
   diversity are coupled at +0.992 across sources; the fusion gain predicted by the equal-magnitude
   variance model inverts in sign exactly where decorrelation is largest.
