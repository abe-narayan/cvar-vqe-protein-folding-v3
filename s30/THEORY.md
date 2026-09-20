# S30 THEORY (lane T)

The sprint's charter (§20) asks: *what is the deepest reason this system cannot turn optimisation
into structural accuracy, and what mathematical reformulation could remove that reason?*

This document answers it in three movements.

1. **The CVaR question is mis-posed.** The tail never stops being a prefix; it is always the prefix
   of the order induced by `grad V` at the optimum (§2). The real dichotomy is whether that order
   is exogenous, and the reachable class when it is not is capped (§3).
2. **The bit question is mis-posed.** The register is a codebook index, not a channel, so no bits
   went missing; the readout's 7 bits realise 36.6 (§5). The value of a bit obeys a measured law
   (§6) which settles the allocation question (§7).
3. **What is actually missing has a dimension and a price.** Six per-target coefficients, worth
   0.855 bits for 3.00 Å and 3.78 for 2.50 Å, in a subspace the pipeline can already compute for
   nothing (§8).

Ledger: **S30-L9** (§2-§4), **S30-L14** (§5-§7), **S30-L15** (§4, §8). Pre-registration:
`s30/PREREG_S30_T.md`, committed before each run, with one declared defect (P5d, §8.4).
Every ORACLE quantity is labelled. Nothing here proposes or prices a deployable operator.

---

## 1. NOTATION, AND THE ONE OBJECT EVERYTHING IS ABOUT

`D = 500` retrieved candidates per target; `W_x ∈ R^(n×3)` the CA coordinates of candidate `x`,
superposed onto the pool medoid in a frame fixed **once and exogenously** (§4.1 says why this is
not a detail); `c` the production average (the uniform coordinate mean of the tie-safe DIS top-75);
`t` the native in the same frame; `e = t - c` the **oracle error**, `d = 3n - 6 = 32.88` on average
over the 126. `p ∈ Δ^(D-1)` a readout distribution, `α` the CVaR level,

    Λ(p) = { λ : 0 ≤ λ ≤ p,  1'λ = α }

the **tail polytope** — the exact feasible set of `core/quantum.py:290 cvar_from_probs`, whose
lines 306-308 compute one of its vertices by a cumulative scan.

S29's displacement bound, which this document uses as a transfer function and does not re-derive:
for a native-free operator emitting `c + u`, with `ρ = cos(u, e)`,
**`RMSD_achievable = RMSD_prod · sqrt(1 - ρ²)`**.

---

## 2. THEOREM T1 — THE TAIL IS ALWAYS A PREFIX

> **T1.** Let `V` be differentiable on `Λ(p)`. At every KKT point `λ*` of `min_{Λ(p)} V` there is a
> scalar `μ` with
>
>     grad V(λ*)_x < μ  ⟹  λ*_x = p_x ,        grad V(λ*)_x > μ  ⟹  λ*_x = 0 .
>
> Hence `λ*` is a prefix of the order induced by `grad V(λ*)` — **for every `V`**.

*Proof.* `Λ(p)` is a box `0 ≤ λ ≤ p` intersected with one equality. The Lagrangian is
`V(λ) - μ(1'λ - α) - ⟨s, λ⟩ + ⟨r, λ - p⟩` with `s, r ≥ 0`. Stationarity gives
`grad V = μ·1 + s - r`. Complementary slackness kills `s_x` where `λ*_x > 0` and `r_x` where
`λ*_x < p_x`. If `grad V_x < μ` then `r_x > 0`, so `λ*_x = p_x`; if `grad V_x > μ` then `s_x > 0`,
so `λ*_x = 0`. ∎

### 2.1 What this says about the deployed stage

The deployed CVaR is `V(λ) = ⟨E, λ⟩`, so `grad V = E`, **constant in `λ`**. That — and nothing
about CVaR, quantiles, order statistics or the circuit — is why one classical sort reproduces it.
`cvar_from_probs`'s cumulative scan is the greedy solution of a fractional knapsack, and the
rearrangement inequality is the entire content of the set-equality theorem
(`s25/QUANTUM.md` §5). The theorem is not a fact about this pool, this Hamiltonian or this ansatz;
it is a fact about linear programs over `Λ(p)`.

### 2.2 The question that replaces the charter's

"Under what conditions does the tail stop being a prefix?" has the answer **never**. The
load-bearing question is:

> **the tail stops being reproducible by one classical sort iff the ordering map
> `λ ↦ grad V(λ)` has more than one fixed point.**

Define the **order operator** `Φ(λ) = prefix_α(grad V(λ))`. The deployed case is `Φ` constant — one
step, one sort. A general `V` makes `Φ` a genuine self-consistency problem, and the number of its
fixed points is the combinatorial content.

### 2.3 A correction to S29's section 4, and to the brief

S29 §4.1-4.3 specify the tail-then-aggregate lift with the tail's order **fixed by a per-state
scalar**, on well-posedness grounds: "**Recommendation: keep the order fixed by a per-state scalar**
and let only the aggregate be structural" (§4.3), because a self-consistent structural order can
make `p ↦ R` multivalued. S29 §4.5 then claims "the set-equality theorem fails for (b)".

**These are incompatible, and T1 says which one wins.** With an exogenous order, `grad`(tail
formation)` = E` and by T1 the emitted support is the `E`-prefix reaching `α`. The optimiser's
freedom is `m`, plus deletions where `p` has exact zeros — which is *verbatim* the deployed
situation (`s29/DATAPATH.md` stage 9: "the whole quantum stage reduces, for the structure, to
choosing `m`").

S29 §4.5's three-state witness (`W_1=(1,0), W_2=(-1,0), W_3=(0,0.1)`; optimal 2-subset `{1,2}`)
is **correct about the free subset problem and unreachable by the readout it was written for**:
with `l_3 < l_1 = l_2`, state 3 is in every α-prefix, so emitting `{1,2}` requires `p_3 = 0`
exactly — the hole mechanism, which `s25/QUANTUM.md` §5 already showed is driven by "an amplitude
pattern that is a function of the **same energies** the classical rank order already uses."

| S29 §4 claim | status |
|---|---|
| `V(S) = f(mean_{i∈S} W_i)` is neither additive nor monotone | **stands** |
| the free-subset optimum is non-prefix (11/12 pairs, 12/12 five-subsets, S29-L25) | **stands** |
| the gap is not measure-zero | **stands** |
| **the lifted CVaR readout reaches it** | **WITHDRAWN** — §4.3's own recommendation forbids it |
| "the first formulation whose classical counterpart genuinely goes away" | **WITHDRAWN as stated** |

S29-L25's second clause — "the escape buys nothing: those subsets' averages are worse than the
prefix on 10 of 12" — was measured on the **free** optimum, i.e. on sets the readout cannot emit.
It is therefore not evidence about the lift in either direction, and should not be cited as such.

---

## 3. THEOREM T1b — THE PRICE OF MAKING THE ORDER ENDOGENOUS

Take `V(λ) = f(R_λ)` with `R_λ = (1/α) Σ_x λ_x W_x`. Then

    grad V(λ)_x = (1/α) ⟨ grad f(R_λ), W_x ⟩ ,

so by T1 the tail is the α-prefix of the candidates **projected onto one self-consistently
determined direction `g = grad f(R)`**. The readout now has `d` degrees of freedom instead of 1.
The price:

> **T1b.** Under an endogenous linear order the reachable tails lie in the halfspace-cut class
> `H = { {x : ⟨g, Ψ_x⟩ ≤ τ} : g ∈ R^d, τ ∈ R }`, of VC dimension `d+1`. By Sauer–Shelah the class
> carries at most `log₂ Σ_{i≤d+1} C(D, i)` bits.

### 3.1 Rank collapse, measured before any compute (P3, HELD) — **AND THE FEATURE SPACE IS PART OF THE NUMBER**

The cap is governed by the **effective** dimension of the candidate feature matrix, not its nominal
one: directions past the leading few separate candidates by amounts below the objective's own
resolution. Measured on **126 real pools of 500** (`s30_T_bits.py`, `s30_T_spec.py`), and verified
independently by lane Q on 12 pools through the **deployed** path (`s29_O_ladder.load_pool` +
`s28_A_amp.Frame`, not `s8/generate_univ`):

| feature space | `r_stable` (T / Q) | λ₁/Σλ | k90 | k99 |
|---|---|---|---|---|
| **pair-distance** `A_dist` (500 × npairs) | **1.859** / 1.862 | 0.554 / 0.544 | 5.61 / 5.6 | 21.13 |
| **coordinate** `A_coord` (500 × 3n) | **3.404** / 3.619 | 0.300 / 0.286 | 11.52 / 11.2 | 27.02 |

`A_dist` variance share PC1..PC6: 0.554, 0.167, 0.107, 0.046, 0.025, 0.019.

> **SCOPE, at lane Q's request and adopted in full: the collapse is a property of the
> PAIR-DISTANCE feature space, not of the pool.** In coordinate space `r_stable` is 3.4–3.6 with
> `k90 = 11.2` directions — comfortably above my P3b threshold of 2.0. **The sparse/weighted
> readout and every second-moment construction act on COORDINATES.** So "stable rank 1.86, the lift
> is a relabelled one-dimensional sort" is correct for a **distance-map** lift and an
> **overstatement** for a **coordinate-space** one. The number must never appear without its
> feature space in the same sentence; it was already circulating without one.

Caps against the **300.6 bits** of free 75-subset choice:

| `d_eff` | source | halfspace VC | cap (bits) |
|---|---|---|---|
| 2 | `A_dist` stable rank | 3 | **24.3** |
| 6 | `A_dist` k90 | 7 | **50.4** |
| 11 | **`A_coord` k90 — the relevant one for a coordinate lift** | 12 | **≈ 80** |
| 21 | `A_dist` k99 | 22 | 126.7 |
| 33 | nominal `3n-6` | 34 | 175.5 |

**My registered decision rule (stable rank < 2.0 ⟹ closed at the encoding level) fires for a
distance-map lift and does NOT fire for a coordinate-space one.** Even at the coordinate figure the
lift is an 11-parameter continuous family — a ~220-bit collapse below the set choice it was
designed to open — so the *verdict* is unchanged, but it now rests on §4.2–§4.3's direct
measurement rather than on the rank threshold alone. The check cost 76 seconds over data already on
disk; S29's first non-diagonal Hamiltonian failed for the same reason discovered *after* the
attempt.

### 3.2 The trichotomy

| order | tail | classical counterpart | endpoint channels |
|---|---|---|---|
| **exogenous** (deployed; S29 §4.3 as recommended) | `argsort(E)[:m]` | one sort | **1** (the integer `m`) |
| **endogenous, `f` convex** | halfspace cut, unique up to ties | Frank–Wolfe, `O(1/ε)` linear-minimisation oracles, each a sort | `d_eff` ≈ 2–6 |
| **endogenous, `f` non-convex** | halfspace cut, several fixed points | greedy + local search over ≤ 175 bits | `d_eff`, genuinely combinatorial |

Row 2 is classically polynomial, so it is not a quantum opening. Row 3 is the only one that is, and
T1b prices its search space at 24–175 bits against 300.6.

### 3.3 CVaR's residual role, which is exactly one thing

Under an endogenous order the cut is a pair `(g, τ)`. `g` comes from `grad f`; **`α` supplies `τ`**.
That is CVaR's whole remaining function — it turns a direction into a set with a *budget* rather
than a threshold, which is what keeps the map continuous in `p`. It is not dispensable and it is
not doing selection.

---

## 4. THE SECOND-MOMENT CLASS: THE ONLY SURVIVING ESCAPE, AND ITS MEASURED VERDICT

Add a dispersion term: `V = f(R_λ) + μ·h(Σ_λ)` with
`Σ_λ = (1/α) Σ λ_x (W_x - R)(W_x - R)'`. Then `grad V_x` is **quadratic** in `W_x`, the cut is a
quadric, and VC dimension goes `d+1 → (d+2)(d+1)/2`: at `d_eff = 6`, 7 → 28, and the cap 50.4 →
**152.1 bits**. The same class is S29's own independent pointer — "operators that read **the pool's
own dispersion** rather than the posterior's marginals" (THEORY_SUMMARY §2, from
`cov(a,n) > cov(b,n)`).

### 4.1 The full design analysis, with a named failure mode per line

Not equations for appearance: each line below exists because it kills a specific way this
construction fails.

| property | analysis | failure mode addressed |
|---|---|---|
| **well-definedness** | `R`, `Σ` well defined on `Λ(p)`, `Σ ⪰ 0`. The superposition frame must be **exogenous** (the pool medoid), never drawn from the tail | a tail-dependent frame is an *unstated operator* (`unstated-operators-align-with-your-hypothesis`) |
| **differentiability** | `∂Σ/∂λ_y = [(W_y-R)(W_y-R)' - Σ]/α`, smooth. `f` kinks at the 17 bin centres; Clarke subdifferential non-empty, subgradient method well posed. Quantile crossings are measure zero | non-smoothness masquerading as non-convergence |
| **normalisation** | `h` must be **scale-free** (a ratio to the pool's own dispersion), else it is minimised by shrinking the set | **contraction artefact** — the average already contracts 22–26% |
| **trivial solutions** | `h` must *match* a target dispersion, not minimise one; otherwise the argmin is a singleton with `Σ = 0`, and the α-budget does not prevent it because `p_x → 1` is reachable | **entropy-induced trivial solution** |
| **identifiability** | the correction's direction is `sign(tr Σ_λ - σ̂²)` — one native-free scalar per target. This is the *only* place the construction can supply the missing sign | the sign is what S29 showed is missing |
| **finite pools** | `Σ` is `33×33` from 75 samples; but the *effective* rank is measured at `k90 = 5.61`, so at most ~6 scalars of `Σ` are estimable | over-parameterised covariance, silently fitted to noise |
| **multimodality** | preserved — a dispersion-matching term does not force unimodality, and can hold 2–3 clusters where a mean-only objective averages the poles back to the pool mean (S29 §3d) | mode collapse in the readout |
| **shot noise** | `tr Σ` from `S` shots in the tail has relative s.e. `≈ sqrt(2/(αS·r_eff))` — 3.3% at `αS = 1800`, `r_eff = 2` | hardware non-realisability |
| **order-statistic artefacts** | any per-target max over `K` variants of `h` is best-of-`K`; `h` must be single, pre-registered, and validated by **split-half transfer** | `grid-oracles-are-order-statistics` |
| **common-mode** | **hard cap**: `tr Σ_λ` is invariant to a common shift, so by the S23-L9 identity (`mean_k|e_k|² = |ē|² + mean_k|d_k|²`, exact to 2.7e-14) a dispersion term **cannot see the 68% common-mode error at all**. It operates on the idiosyncratic 32% only | the dominant error component, unobservable to this class by construction |
| **rank collapse** | measured: `r_stable = 1.859`, `k90 = 5.61`. Quadric VC at `d_eff = 6` is 28, not the nominal 595 | the failure that killed S29's first non-diagonal Hamiltonian |

### 4.2 The measured verdict: closed

**Lane Q's S30-L12 has priority and a better null.** At 128 directions/class/target over all
`m = 1..500`: QUADRIC − LINEAR is **−0.0105 Å at 0.16× MDE (NOT MEASURED)** at ORACLE `m` and
**+0.2059 Å at 1.81× MDE, 37W/89L (WORSE)** at the shipped `m = 75`. My independent M4 (K = 5,000
directions inside the 6-dimensional PC subspace, fixed `M = 75`, `s30_T_quadric.py`, n = 126):

    PREFIX    (deployed DIS top-75, one set)       3.0507 Å   SE 0.1456
    HALFSPACE best-of-5,000                        2.0691 Å   SE 0.0903
    QUADRIC   best-of-5,000                        2.1551 Å   SE 0.0929
    FREE      best-of-5,000 random 75-subsets      2.9259 Å   SE 0.1254   ← matched null

    P4a  QUADRIC − HALFSPACE   +0.0860  2.95× MDE  5/5 folds  21W/105L   WORSE   (HELD)
    P4b  FREE − PREFIX         −0.1248  0.60× MDE  4/5 folds              NOT MEASURED (HELD)

Two samplers, two lanes, same sign. `HALFSPACE ⊂ QUADRIC`, so the deficit is a **searchability**
statement (27 parameters sample more sparsely than 6 at equal budget) rather than a containment
one — but operationally the answer is identical: the enrichment is not findable and the halfspace
ceiling already binds. **The second-moment escape is closed as a ceiling-improver.**

### 4.3 RETRACTED: the halfspace "ceiling" is an order statistic, and priced honestly the class is NEGATIVE

My registered null — best-of-5,000 over *random 75-subsets* — is too weak. It asks whether a
structured class beats an unstructured one at equal budget; it does not ask whether the winning
direction is the same direction twice. Lane Q asked for the right pricing and the coordinator
endorsed it. **Run on my own K = 5,000** (M6, `s30_T_transfer.py`, n = 126), with a **common
direction bank** across targets made comparable by a deterministic native-free sign convention
(each PC's largest-|loading| entry positive), so that column `k` is the same *rule* on every target
and `best_of_k_within` is well posed:

```
                      observed   across-target null   accounted   split-half transfer   k_eff
  HALFSPACE  K=5000   -1.7881         -3.5038           196%        -0.3319 (19%)       118
  QUADRIC    K=5000   -1.5663         -2.9368           188%        -0.3148 (20%)       110
  verdict: NOT A SIGNAL in both (transfer < 25% of the oracle)
```

**The across-target null exceeds the observed gain in both classes**, exactly as lane Q predicted
for my larger `K`. And the decisive check — because a transfer of 19% is not zero and could be
misread as a route:

    PREFIX (deployed)                                3.0507 Å
    HALFSPACE mean over the 5,000 directions         3.8550 Å
    the transferable rule (K-mean + split-half)      3.5231 Å
    vs PREFIX: +0.4724, 1.88x MDE, 5/5 folds, 39W/87L   ->  **WORSE**

> **RETRACTION.** S30-L15 §1 and the earlier reading "HALFSPACE − PREFIX −0.9816 Å, BETTER" are
> **withdrawn as a statement about the class's value.** The 0.98 Å is a per-target best-of-5,000
> and 196% of it is accounted for by the across-target null. The class's *transferable* content
> lands **0.472 Å worse than what ships**, because the direction bank's own mean (3.855) is far
> worse than production and a 19% transfer does not close that. **The halfspace class is reachable,
> not exploitable, and as a rule it is negative.**

What survives unchanged: the *reachability* statement (a 6-parameter order can reach 2.07 Å), which
is what §8 uses, and every quadric contrast, which is a like-for-like comparison of two classes at
matched budget and is not affected by a null that shifts both. `W/L cannot diagnose this` —
only split-half transfer can, and I should have run it in M4 rather than after being asked.

---

## 5. THE BIT ACCOUNTING: THE REGISTER IS A CODEBOOK, NOT A CHANNEL

### 5.1 The accounting, in a currency defined before it was measured

`b_ret = log₂(N*/500)`, where `N*` is the number of uniformly random draws from the same universe
whose ORACLE best-of-`N` matches the BLOSUM-500's. Defined in the prereg precisely so it could not
be chosen afterwards.

| stage | bits of CHOICE | bits DELIVERED |
|---|---|---|
| 5 retrieval, 17,088 → 500 | **3,252** | **+0.69** [SE 0.18, median +0.44] |
| 6 energy, scores → ranks | 3,767 | — |
| 7 top-128 prefix | 300.6 | ORACLE ceiling 0.066 Å (S29-L30); transferable part 0 |
| 9 tail readout, which of 128 | **7** | **+0.036** [MDE 0.385] |

**Under one delivered bit per target, against 7,019 bits of choice consumed.**

### 5.2 Why no bits went missing

Through the displacement bound:

    R = 7.000 index bits:  4.1080 → 1.8978 Å  ⟹  ρ = 0.8869  ⟹  **36.63 displacement bits  (5.23×)**
    R = 8.966 index bits:  4.1080 → 1.7108 Å  ⟹  ρ = 0.9092  ⟹    41.55 displacement bits  (4.63×)

Bits are not conserved across an index because **the 500 deposited backbones are the information
and the index only names one of them**. The charter's premise — the register as a channel through
which 7 bits pass and 1.44 arrive — is the wrong object. The correct statement is the inverse, and
it is harder rather than softer:

> **the readout's 7 bits are worth five times their face value, and the system cannot supply one.**

### 5.3 Two things `bits_delivered` is not, stated because they look the same

- It is **not** a mutual information. `7 − E[log₂ r]` is a realised *search cost* of a particular
  ordering. A score can carry information and rank badly; the standard bounds relating the two
  carry a slack of `log₂(1 + ln 128) = 2.55` bits, which is too large to bite at `K = 128`. So a
  measured 1.44 does **not** upper-bound `I(B; score)`.
- It is **not** an Ångström claim. `operator-consumes-set-mean` (`d_out = 1.16·mean + 0.04·best`)
  means a perfect rank-1 is worth −1.74 Å through argmin and −0.03 Å through the `m = 75` average.
  Nothing in §5 claims any part of S29-L44's 0.757 Å is recoverable by ranking.

---

## 6. THE VALUE-OF-A-BIT LAW

The ORACLE order-statistic ladder inside the deployed pool, in the deployed order (n = 126):

     N       1      2      4      8     16     32     64    128    256    500
     D(R) 4.1080 3.6059 3.1115 2.7289 2.4901 2.2973 2.1230 1.8978 1.8061 1.7108  Å
     SE   0.1810 0.1516 0.1457 0.1297 0.1205 0.1112 0.0996 0.0877 0.0844 0.0791

If candidate quality has CDF `F(x) ≈ κ(x-a)^γ` near its lower endpoint `a` (the Weibull domain),
then `E[min of N] ≈ a + (κN)^(-1/γ)`, i.e. `D(R) = a + c·2^(-R/γ)`. Fitted:
**`a = 1.3312 Å, c = 2.7859, γ = 3.1636, R² = 0.998341`** (P1a HELD). Differentiating:

> **`−dD/dR = (ln 2 / γ)·(D(R) − a) = 0.2191·(D − 1.331)` Å per bit.**

Marginal at `R = 7`: **0.1317 Å/bit** (P1b HELD, registered [0.05, 0.14]). Fitted floor 1.3312 Å
against a true universe floor of 1.3134 and a pool floor of 1.7108; `D(7) − a = 0.567 Å`, so at the
production width the pool is not near its floor and bits still pay (P1c HELD).

**Why this law and not a curve fit.** It has one structural consequence, and that is its whole
purpose: the value of a bit is proportional to the **distance above the class's floor**, divided by
the class's **tail index**. Both are properties of the codebook, neither is a property of the
register.

---

## 7. THE ALLOCATION THEOREM

Since `−dD/dR = (ln2/γ)(D − a)` and `D` is common at the branch point, **allocations are comparable
only through `(D − a)/γ` — a class's ORACLE floor and its tail index. Cardinality is irrelevant.**

| allocation | Å/bit | floor | note |
|---|---|---|---|
| **candidate identity** (deployed) | **0.132** at R = 7 | 1.331 Å | wins |
| subset cardinality `m` | 0.044 | — | S29-L44's −0.3079/7; dominated 3.0× — **and the law explains S29's measured 3.5×** |
| mode / basin index | saturates ≈ 1.6 bits | — | 2–3 populated clusters; a sub-class of the first |
| **torsion / configuration** | **infeasible** | — | below |
| hierarchical / factorised | capped | — | block-independent allocation cannot touch the 68% common mode; bounded by the idiosyncratic 32% |

**The torsion arithmetic — a counting statement, not a preference.** Mean `n = 12.96`, so
`2n = 25.9` torsions. Naming one Ramachandran basin per residue costs `2n·log₂ k` bits:
**25.9 at k = 2, 41.1 at k = 3, 51.8 at k = 4, 77.8 at k = 8.** The deployed register is **7**
qubits (`core/pipeline.py:181`) or 9 in the harness — **0.27 bits per torsion, where 1 bit names a
single basin**. Any generative proposal over torsion space needs 4–8× the register this project has
ever run, before expressivity is discussed. Consistent with `phi-carries-no-sequence-signal`: the
channel is 10.4° of ψ, so the class's `γ` is large as well as its bit demand.

> **Answer to the charter.** Binary candidate indexing **is** optimal among the measured classes, by
> 3×, and it is capped at 1.331 Å however many bits are added. Since the charter's target is 2.50 Å,
> the allocation is not floor-limited. **The encoding is not the hidden bottleneck.**

**Scope.** That is a closure about the encoding *among the measured classes*. It is not a claim
about the operator (§4, §8) or about any stage (lane F's S30-L2 locates a +1.767 Å loss in the
distogram's 500→75 filter on FAIL18). Three lanes now point at three different components and the
report must keep them distinct.

---

## 8. WHAT IS ACTUALLY MISSING: SIX COEFFICIENTS

### 8.1 The ρ ↔ bits dictionary

`I = −(d/2)·log₂(1 − ρ²)`. Its small-ρ limit is `I ≈ (d / 2ln2)·ρ²`: **bits are quadratic in ρ**.
That is lane D's "the geometry squares it" (S30-L5) — and the squaring is not a penalty, it is what
makes `ρ²` the **additive** quantity across stages and fields, which is what lets the sprint argue
in one currency.

### 8.2 The subspace, and the result that inverts the question

Any native-free operator that reweights or selects pool members emits `u = Σ_x w_x (W_x − c)` with
`Σ w_x = 1`, so it lies in the span of the pool's deviations about `c`. Project `e` onto the top-`k`
principal directions of that deviation matrix (n = 126, ORACLE; `s30_T_combo.py`):

     k       ‖Π_k e‖²/‖e‖²     ρ_ceiling      isotropic null k/d
     1           0.1997         0.4469            0.0257
     2           0.3514         0.5928            0.0514
     6           0.6072       **0.7792**          0.1543
    21           0.9261         0.9624            0.5401

**The oracle error is 8× more concentrated in the pool's leading deviation directions than isotropy
predicts.** Six directions contain a `ρ = 0.779` point — **1.91 Å** through the bound, against
production's 3.051. The subspace is **not** the obstruction, which is the opposite of what the 68%
common-mode result had led me to expect. (There is no maximum over draws anywhere in §8.2, so the
across-target-null critique of §4.3 does not apply to it.)

**Two instruments agree to 0.16 Å**: the 6-direction subspace ceiling is 1.91 Å and the
6-parameter halfspace class reaches 2.069 Å. A linear order nearly saturates the subspace that
holds the answer.

**Both figures are REACHABILITY, not value.** §4.3 retracts any reading of the 2.069 Å as an
achievable gain: 196% of it is an across-target order statistic and its transferable part is
0.472 Å *worse* than production. The two numbers say that the answer is **inside** a
six-dimensional space the pool hands you — they say nothing about anyone being able to point at it.
That distinction is the whole of §8.4.

### 8.3 The price, evaluated in the right space

The subspace is **free** — the pool's PCA is native-free and per-target — so the dictionary should
be evaluated inside it, where a field needs only `cos = ρ_target/0.779` against the true direction.
Two derivations agree exactly (sphere rate–distortion; and the covering number of `S^5` by caps of
angular radius `arccos c`, leading term `(d−1)log₂(1/sin θ)`):

    ρ_target                 FULL SPACE (d = 32.88)    WITHIN THE POOL'S TOP-6
    0.1128 (best field, S30-L5)   0.304 bits                0.076 bits
    0.140  (B2's ceiling)         0.470                     0.118
    0.358  (3.00 Å)               3.253                   **0.855**
    0.628  (2.50 Å, charter)     11.895                   **3.782**

> **The charter's 2.50 Å needs 3.78 bits of per-target direction information in a subspace the
> pipeline already computes for nothing. 3.00 Å needs 0.855 bits. The best native-free field on
> record supplies 0.076.**

The factor to 3.00 Å is `0.855/0.076 = 11.2`, which is S30-L5's "10× the ρ²" **from the same
algebra** — one witness, not two.

**Withdrawn:** I told the coordinator that "six coefficients at ~2 bits each = 12 bits" matched the
full-space 11.895 and called it two independent routes. It was a coincidence. The two columns are
one formula at two values of `d`, and the 8.1-bit gap at 2.50 Å is exactly the value of knowing the
subspace. The corrected figure makes the missing thing **smaller**.

### 8.4 The deepest reason, stated plainly

The reachable displacement space contains a **1.91 Å** point inside **six** well-conditioned
directions the pool hands you for nothing. A 6-parameter linear tail order reaches **2.069 Å** of
it. The deployed prefix gets **3.051 Å**. And the direction does **not** transfer: lane Q measures
split-half **+0.032 Å (−3%)** at K = 128, and my own K = 5,000 with a common direction bank gives
**−0.332 Å (19%)** against an across-target null accounting for **196%** — whose transferable rule
lands **0.472 Å worse than production** (§4.3).

> **The entire 0.98 Å gap is the problem of locating six per-target coefficients** — S29's
> Neyman–Scott incidental parameter with a dimension attached for the first time. A
> non-transferable direction in a 6-dimensional subspace *is* six per-target coefficients, so lane
> Q's transfer result and this section are the same finding approached from two sides.

It is not the encoding (§7), not the objective's functional form (S29 Theorem 2), not the
Hamiltonian's off-diagonal structure (S29 §3), not the ansatz, not the optimiser, not the pool, and
not the search (§4.2, and S30-L8). It is under four bits of direction per target that nothing
native-free supplies.

**A declared defect.** Clause P5d of my prereg addendum 2 — predicting `‖Π_6 e‖²/‖e‖² ∈ [0.15,0.40]`
— was written **after** M5 had run. It is not a pre-registration, its falsification carries no
evidential weight, and it is retained only so that a wrong prior of mine stays auditable. P5a–P5c
(lane D's combination) are prospective and are the only clauses to be scored.

---

## 9. B2's SCOPE, WHICH THE BRIEF ASKED ME TO ESTABLISH

S29's assumption **B2** quantifies over "every field constructible from the present information".
The brief was right to doubt it as a universal statement, and right that 21 measured fields do not
establish a supremum. But the fix is not more fields — it is noticing that **B2 is not a claim about
all functions; it is a claim about all functions of one statistic.**

Every native-free field in the record is a function of the same `S = (posterior, pool)`. The
supremum over `φ` of `ρ(φ(S))` is attained by the Bayes-optimal `φ*(S) = E[e | S]`, and by the
tower property `Cov(E[e|S], e) = Var(E[e|S])`, so

> **`ρ_max = sqrt( R²(e ~ S) )` exactly** — the square root of the population variance of `e`
> explained by the statistic.

So B2's entire load-bearing content is **one estimable number**: `R² ≤ 0.0196` at `ρ = 0.14`, and
`R² ≤ 0.0127` at lane D's measured best of 0.1128. The charter's 2.50 Å needs `R² = 0.394`.

Three consequences.

1. **The falsifier is a single out-of-fold regression**, not a field survey: fit `E[e|S]` with any
   regressor on the pinned folds and report its out-of-fold `R²`. Above 1.96% and B2 falls. (It
   requires ORACLE labels — which is fine: B2 is a claim about what is *achievable*, and you only
   need the labels once, to decide whether to build.)
2. **S30-L50's nine in-band channels do not contradict B2.** In-band skill is a *scalar* rank
   statistic; `R²` is the *joint* one over 33 dimensions. Nine channels with in-band skill can
   coexist with `R² < 2%`, and `decorrelated-errors-exist-but-are-unusable` (fusion gain going as
   the square of the weaker channel) is exactly the mechanism by which they do.
3. **Combination is governed by rank, not count.** For `k` fields of cosine `ρ_0` with exchangeable
   correlation `c`, `ρ_comb² = kρ_0²/(1+(k−1)c)` and `(1+(k−1)c)/k = s_1`, the PC1 share of the
   field Gram. Hence **`ρ_comb = ρ_0 / sqrt(s_1)`** — the count cancels. Registered as P5a–P5c
   before lane D reports: Gram stable rank in [1.3, 3.5], `s_1 ∈ [0.45, 0.80]`, leave-fold-out
   combination in **[0.11, 0.22]** with point estimate 0.15, falsified above 0.25; and the
   in-sample all-21 figure inflated to 0.25–0.45 and not to be read as a ceiling. The orthogonal
   estimate of 0.33–0.36 corresponds to `s_1 = 1/11`, i.e. a Gram of full rank 11 — which nothing
   built from one distogram and one pool has.

**What B2 covers, then:** all functions of `(posterior, pool)`, with a one-number characterisation
and a one-measurement falsifier. **What it does not cover:** a genuinely new information source,
which changes `S` and therefore changes the number. That is unchanged from S29 and it is the only
door §8's six coefficients could come through.

---

## 10. REGISTERED PREDICTIONS AND THEIR OUTCOMES

| clause | registered | measured | verdict |
|---|---|---|---|
| P1a ladder form | R² ≥ 0.99, γ ∈ [2,6] | R² 0.9983, γ 3.164 | **HELD** |
| P1b marginal at R=7 | [0.05, 0.14] Å/bit | 0.1317 | **HELD** |
| P1c floor / headroom | a ≤ 1.6, D(7)−a ≥ 0.4 | 1.331, 0.567 | **HELD** |
| P2a retrieval yield | `b_ret` ≤ 3.0 bits | +0.69 | **HELD** |
| **P2b retrieval Å gain** | **+0.10 to +0.30 Å** | **0.0715 (0.79× MDE)** | **REFUTED, 2–4×** |
| P3a rank collapse | rs_dist ∈ [1.2,4] | 1.859 | **HELD** |
| P4a quadric vs halfspace | within 0.10 Å | +0.086 | **HELD** |
| P4b free vs prefix | within 0.15 Å | −0.125, not measured | **HELD** |
| **P4c class gain over prefix** | **< 0.6 Å** | **0.982 Å** | **MISSED by 64%** (falsifier was 1.0) |
| P5a–c lane D's combination | see §9 | *pending* | prospective |
| P5d subspace concentration | — | 0.607 | **NOT A PREREGISTRATION** (§8.4) |
| M6 halfspace transfer (lane Q's request, unregistered by me) | — | 196% accounted, transfer 19%, rule +0.472 Å worse | **RETRACTS my §4.2 reading** |

**The retraction is the largest correction in this document and it is to my own headline.** I
reported the halfspace class as BETTER than the deployed prefix by 0.98 Å against a matched
random-subset null, and that null was the wrong one. Lane Q named the right one, the coordinator
endorsed it, I ran it on my own K = 5,000, and it says the class is negative as a rule. The
reachability statement survives and is what §8 rests on; the value statement does not.

**Both misses ran toward my own hypothesis** — P2b because I wanted retrieval to be delivering
something measurable, P4c because I over-weighted `operator-consumes-set-mean`'s flattening and
under-predicted what an ORACLE class can do. That is the pattern
`unstated-operators-align-with-your-hypothesis` names, and recording it is the reason the rest
should be trusted. A third claim — the "12 bits ≈ six coefficients" aside of §8.3 — was wrong and
was **unregistered**, which is the most dangerous kind; it is flagged in §8.3 rather than deleted.
