# WORKSTREAM D — FINDINGS, SPRINT 22

Role: adversary, mathematician, librarian. Everything below is independently reproduced from
source (scripts in `s22/results/`) or independently derived (no scripts, pure math, checked
numerically where a check was cheap). Nothing here is quoted from a lane's table without having
recomputed it or derived it myself. Findings that damage a coordinator claim are marked **DAMAGE**.

---

## 1. AUDIT OF s22 LEDGER L1 — THE ROUTING CEILING **[sent to coordinator, formalised here]**

**Claim under audit.** Best fixed arm 3.048, ORACLE per-target routing 2.566 over 13 native-free
arms, headroom 0.482 Å, concentration 50%/21 targets, 80%/50 targets, median 0.317, corr(headroom,
incumbent)=+0.564, 12/126 incumbent-already-best.

**Method.** `s22/results/d_route_audit.py`. Sources: `s21/results/poolgap.json` (8 candidate arms:
`avg_500,150,75,20,5,1`, `medoid_all`, `medoid_75`, all on the pool's own WINDOW coordinates) and
`s21/results/latentsel.json` (6 candidate arms: `lat_rand1, lat_argmin, lat_avg75, lat_medoid`, on
the generative latent's BUILT-CHAIN coordinates, plus `pool_argmin`/`pool_avg75` REBUILT onto that
same manifold — a third, distinct basis for "the pool's argmin/average"). The coordinator's prose
("m-ladder(6) + medoids(2) + four latent arms") parses to 12, not 13; I built every plausible
12/13/14-arm completion.

**DAMAGE (minor, process not conclusion).** Two identity hazards inside `poolgap.json` itself:
`avg75_medoid75 ≡ avg_75` (max|diff| = 0.0, literally the same assignment in source) and
`pool_argmin ≡ avg_1` (max|diff| = 3e-14, float noise only) — BRIEF §9's "identity counted as a
trial" hazard, present in the data regardless of whether the coordinator's 13-arm list used them.

**Result.** All five reconstructions (12, 13, and 14 arms) give headroom **0.5047-0.5063 Å** —
tightly clustered, i.e. the exact arm-set ambiguity does not matter for the headroom's SIZE. This is
**0.02-0.03 Å larger than the coordinator's 0.482**, an unresolved discrepancy of unknown origin
(possibly a genuinely different arm list, possibly a different random draw of the stochastic
`lat_rand1`/`lat_avg75_rand` controls — the static `.json` used here is presumably identical to
whatever the coordinator read, so a differing arm menu is the more likely explanation). Concentration
diagnostics replicate in SHAPE, not to the digit: top21 51.1-51.2% (claim: 50%), top50 81.5-81.8%
(80%), median per-target headroom 0.313 (0.317), corr(headroom, incumbent RMSD) +0.60 (+0.564),
incumbent-already-best on 10-13 targets (12).

**Verdict: the ceiling is real, its size is robust to the stated construction ambiguity, and its
qualitative shape (broad, not needle-in-haystack, larger where the incumbent is worse) replicates
independently. The EXACT number 0.482 is not independently certified and should not be quoted to
three decimals until the coordinator supplies the literal 13-item list.**

**Independent bonus measurement.** The aggregate window-vs-rebuild basis price, computed directly
here (not quoted from Sprint 21): pool_argmin rebuild−window `+0.0528 [MDE 0.117]`; pool_avg75
rebuild−window `+0.0155 [MDE 0.029]`. The second reproduces Sprint 21 L20's `+0.016` almost to the
decimal — an independent cross-check that L20's number is sound, obtained as a side effect of this
audit rather than by re-reading L20.

---

## 2. THE CORRECT NULL FOR A PER-TARGET ROUTER OVER CORRELATED ARMS **[Duty 3 derivation; DAMAGE to the coordinator's method, not his conclusion]**

The coordinator's question: *"What IS the correct null for a per-target router scored across
strongly correlated arms? Is my diagnosis right? Is 'the ceiling is legitimate, the estimate needs
held-out folds' the correct disposition, or is there a bias in the ceiling too that I am waving
away?"*

### 2.1 The ceiling needs no null at all — this is the thing to fix first

Write `X_{a,t}` for the (deterministic) RMSD of arm `a` on target `t`. Both quantities in question
are then deterministic functions of the `13×126` table `X`:

    best_fixed      = min_a  (1/T) Σ_t X_{a,t}
    oracle_routing  = (1/T) Σ_t  min_a X_{a,t}

`min_a X_{a,t} ≤ X_{a*,t}` for every `t` and in particular for `a* = argmin_a mean_t X_{a,t}`, the
arm achieving `best_fixed`. Averaging a pointwise inequality preserves it:

    oracle_routing = (1/T) Σ_t min_a X_{a,t}  ≤  (1/T) Σ_t X_{a*,t} = best_fixed.

**`headroom = best_fixed - oracle_routing ≥ 0` is a theorem, not a measurement.** There is no
sampling distribution anywhere in this inequality — every `X_{a,t}` is a specific real number (an
RMSD between a specific candidate structure and a specific native structure), so there is nothing
for a null hypothesis to be a null hypothesis ABOUT. Asking "is 2.566 [or 2.542, whichever] a
statistically surprising minimum" is a category error of exactly the kind the project's own
`latent_oracle` entries were careful never to make (L14/L17 label it a ceiling and stop). The
coordinator's within-arm shuffle is not a wrong ANSWER to this question; it is an answer to a
DIFFERENT, well-posed question — see 2.2 — mistakenly aimed at the ceiling.

### 2.2 What the within-arm shuffle actually measures, and why it goes the "wrong" way

Shuffling `t → σ_a(t)` independently per arm `a`, then recomputing `E_t[min_a X_{a,σ_a(t)}]`,
answers: *"if the same 13 marginal per-arm distributions existed but arms were mutually independent
across targets (no shared difficulty factor), how low would the per-target minimum be?"* Because
target difficulty is a dominant factor common to ALL arms (every arm is worse on hard targets — this
is visible directly in the data: e.g. the incumbent's own RMSD spans 0.7 to 6.4 Å across deciles,
§3.2 below), decorrelating the arms lets a single easy-target's low value from ANY one of the 13
independently-shuffled arms fill the "slot" for what was originally a hard target. This is *easier*
to achieve than requiring the SAME target, under correlated arms that are all simultaneously bad on
it, to have even one of its 13 arms be good — so the decorrelated minimum runs LOWER (here, 1.115
vs 2.566/2.54), exactly as observed. **The coordinator's mechanistic diagnosis is correct.** But the
right conclusion is not "the shuffle is mis-specified so use a different shuffle" — it is that **no
shuffle of this table is the right instrument**, because the object being investigated (the ceiling)
is not a sample statistic.

### 2.3 The construction that IS correctly specified, derived from first principles

The real question — is the 0.48-0.51 Å gap *exploitable* — is a decision-theoretic one, and its
null has an exact, non-simulated answer. Let `f` be any (possibly randomized) rule mapping
native-free, pre-decision features of target `t` to an arm choice, and let its achieved value be
`V(f) = E_t[X_{f(t),t}]`.

**Proposition.** If `f` is measurable with respect to information that is *statistically independent
of the joint law of `(X_{1,t},...,X_{13,t})`* — i.e. `f` carries zero real signal — then
`E[V(f)] = best_fixed` exactly, not merely "no better than." *Proof sketch:* condition on the value
`f(t)=a`; independence of `f` from the outcome table means the conditional law of `X_{a,t}` given
`f(t)=a` equals its marginal law, so `E[X_{f(t),t}] = Σ_a P(f(t)=a) E[X_{a,t}]`, a probability-weighted
average of the 13 arms' own means — which is minimized (over all possible weightings, including the
degenerate ones a deterministic rule uses) by putting all weight on `argmin_a E[X_{a,t}]`, i.e. by
`best_fixed`. A constant rule achieves this exactly and cannot be beaten by a noisier zero-information
rule (a mixture of means is never below its own minimum term). So **the correct "no-skill" floor for
an achievable router is `best_fixed = 3.048` by construction — this is arithmetic, not something a
permutation procedure needs to manufacture**, and it is already what BRIEF states ("the correct
baseline for a router is the best FIXED arm"). It requires no null of its own for the same reason the
ceiling doesn't: it's the other endpoint of a deterministic sandwich, `best_fixed ≥ V(f) ≥
oracle_routing` for every rule `f`, with equality on the left iff `f` carries no real information and
equality on the right iff `f` is the oracle.

### 2.4 What DOES need a null, and it is one level up from where the coordinator was looking

An **achieved, held-out-fold-estimated** router (fit `f` on training folds, score on the held-out
fold) is a genuine random variable across folds/resamples even though the underlying `X` table is
fixed, because `f` itself is fit from noisy finite data (see §1's demonstration: my length-quartile
router, fit exactly this way, returned `+0.031 [SE 0.025, MDE 0.070]` relative to `best_fixed` — a
real, CI-bearing quantity, unlike the ceiling). **This is where a null belongs**, and BRIEF's own
Rule 0 clause 3 already supplies the correct one: if the router-fitting PROCEDURE itself searched
over multiple candidate constructions (feature choices, bin counts, candidate-arm subsets) and the
result reported is the best of those, its null is the distribution of the **best of K such
held-out estimates under a procedure with no real signal**, not the distribution of a single one.
This is a second, INDEPENDENT layer of best-of-K bias on top of the per-target `min_a` inside the
ceiling, and it is the multiple-comparisons correction the coordinator's OWN M4 memory ("a best-of-K
result's null is the distribution of the maximum") already requires but had not yet been pointed at
the router-selection step itself. **Recommendation, adopted as a rule for Workstream C**: report how
many router constructions were tried before the one quoted, and either pre-register the single
construction or apply a nested/second-fold correction for the search over constructions.

### 2.5 Disposition

**The coordinator's diagnosis is right and his final disposition is right, but for a cleaner reason
than the one given.** The ceiling needs no null because it is not a statistic. The estimate of an
achievable router's value needs held-out folds not merely as damage control for a bad permutation
test, but because it is the ONLY well-posed random quantity in this whole analysis, and it further
needs a multiple-comparisons correction if more than one router construction was tried before the
reported one. There is no hidden bias in the ceiling being "waved away" — the ceiling is exactly what
it claims to be, an upper bound on achievable value, and needs no defence beyond the one-line
inequality in §2.1.

---

## 3. THE RETRIEVAL/LATENT DIFFICULTY INTERACTION (L2/L3) **[DAMAGE — sent to coordinator, formalised here]**

Full numbers in `s22/PREREG_D.md` D-AUDIT-3 and `s22/results/d_l2_audit.py`. Restated briefly:
`corr(ship_avg75, lat_avg75) = +0.819`, OLS slope of `lat` on `inc` = `0.580` (R²=0.671). Because
`Cov(inc, lat-inc) = Cov(inc,lat) - Var(inc)`, ANY two variables related by a slope-<1 line produce a
mechanically decreasing `d = lat - inc` as a function of `inc`, with slope exactly `beta-1 = -0.420`,
with **zero requirement that `lat` carry genuine difficulty-regime information**. The linear-mechanical
component alone reproduces 81-202% of each of L2's four quartile means (overshooting on 3 of 4); the
residual left over is small, non-monotone, and 2 of 4 CIs include zero. A placebo that shuffles
`lat_avg75` within length-deciles (destroying true pairing, preserving marginals and the length
association) reproduces an EVEN LARGER version of the reported sign-flip shape purely from the same
algebra (+2.25 to -1.70 vs the real +1.47 to -0.46) — confirming the shape, on its own, is not
diagnostic of anything beyond "these two variables are correlated with slope < 1," which is the
uncontroversial and already-known fact that both methods track true difficulty.

**Consequence.** L2 downgrades from "a routing signal, cleaner than anything in Sprint 21" to a
single global fact (the latent's error scales more weakly with difficulty than the incumbent's) with
a small, non-monotone, partially-inside-CI residual beyond that fact. This directly explains L3
without needing to invoke proxy imprecision: there may not be much regime-dependent structure to find
past the linear-slope fact, and that fact itself is not native-free-actionable (the quartile split
consumed the incumbent's own RMSD to the native). It also matches §1's own achievable-router result
(the length-quartile router captures ≈0%, if anything slightly negative, of the L1 ceiling) — both
audits land on the same practical conclusion by different routes, which is a point in favour of both.

**FAIL18 not independently checked.** The coordinator's "72% vs 24%" win-rate claim on a specific
named 18-target failure set was not verifiable from the artefacts pulled here (the set membership was
not identified in `latentsel.json`); flagged OPEN rather than asserted either way.

---

## 4. THE COMPACTNESS CHANNEL'S Å CONVERSION (BRIEF's "central conversion") **[null result, pre-registered before running]**

Full pre-registration in `s22/PREREG_D.md` D-PREREG-1. L27/L28 (Sprint 21) demonstrated `rg_z`/
`rg_gap` — a native-free compactness-disagreement channel already latent in the shipped distogram —
correlates (partial ρ ≈ 0.34-0.38) with an ORACLE rank-percentile label, on the n≤13 panel (75
targets) only, and explicitly left the Å value unmeasured. Extended the channel to the full n=126
pool (`s22/results/d_rg_full126.py`) and built the most direct honest test of the conversion
question: a 5-fold held-out router choosing among ladder arms (top-75 vs a wider/narrower average or
the medoid), conditioned on `rg_z`/`rg_gap`/`|rg_z|` via a 3-bin quantile split fit on training folds
only (`s22/results/d_rg_route.py`).

**Result: every construction is null.** Best case `-0.021 [SE 0.019, MDE 0.053]`, NOT MEASURED; two
of six candidate-arm pairings (`avg_75` vs `avg_500`, `avg_75` vs `medoid_75`) never route away from
`avg_75` in ANY training fold — there is no in-sample crossover for those pairs at all, let alone an
out-of-sample one, at this bin resolution.

**Disposition.** The falsifier did not fire (no construction beat the incumbent past its own MDE).
This does **not** close the L27/L28 channel — a coarse 3-bin threshold on a single feature is a weak
instrument for a ρ≈0.35 correlation, and this is one specific, reasonably thorough but not
exhaustive attempt (a continuous/regression-based combiner, or a joint model with other features,
is untested; both would need their own pre-registration before running so as not to inherit a
best-of-K search bias, per §2.4). **What this DOES establish**: nobody should currently plan a
pipeline change on the strength of the L27/L28 signal without a stronger router than the one tested
here, and BRIEF's "central conversion" is not merely unattempted — one honest, pre-registered attempt
at it has failed, and that failure is itself now part of the record.

---

## 5. THE CVaR FACE AND DEGENERACY-BREAKING **[Duty 3 derivation]**

**Setup.** Sprint 21's Q6 (adopted, `s21/LEDGER.md` L1c) states: *for a pool-restricted selector
whose training objective and readout energy are the same H, CVaR-VQE and argmin have the SAME
optimal answer; CVaR changes the sampling distribution, not the selected point.* This section derives
that fact precisely and then answers BRIEF §3's specific question: what does the degeneracy of the
CVaR-optimal set imply for an attempt to break it with an added structural term?

**5.1 The optimal face, precisely.** Fix a finite pool of candidates with energies
`H_{(1)} ≤ H_{(2)} ≤ ... ≤ H_{(m)}` (ties allowed) under Hamiltonian `H`, and let `p` range over all
probability distributions on the pool (the fully-expressive limit; a real ansatz only reaches a
subset of this simplex, which can only shrink the achievable set below what is derived here). For a
tail fraction `α`, `CVaR_α(H;p)` is the α-quantile-conditional expectation of `H` under `p`. For ANY
`p`, `CVaR_α(H;p) ≥ H_{(1)}` (it is an average over a subset of `H`'s support, all of whose values are
`≥ H_{(1)}`), with **equality if and only if `p` places at least `α` total probability mass on the
level set `{i : H_i = H_{(1)}\}`** (the argmin, or the tied argmin set). Call this set of
optimal laws `F_α(H)`. It is a genuine FACE of the simplex: fixing `p(\text{argmin level set}) ≥ α`,
the remaining `1-p(\text{argmin})` mass may be placed **arbitrarily** over every other state in the
pool, including the single WORST state, without changing `CVaR_α` at all, and (if there are ties for
the minimum) the mass within the argmin level set may also be split arbitrarily among the tied states.
`F_α(H)` therefore has dimension `m - 2` in general (one degree of freedom fixed by the `≥α`
constraint being tight, one by normalisation) — a large, non-trivial degenerate optimum, not a
point.

**5.2 The readout kills the degeneracy question before it starts, when H is shared.** The deployed
readouts (`vqe_bitstring`, `vqe_modal_bitstring`, `best_seen_bitstring`; verified in
`core/quantum.py`) are order statistics of `H` over states with non-negligible probability under `p`
(or over everything sampled) — **not functionals of `p`'s shape**. For any `p ∈ F_α(H)`, the argmin
level set has probability `≥ α > 0`, so it is certainly represented among sampled/seen states, and
`argmin_i H_i` returns the same value `H_{(1)}` and the same optimal state(s) **regardless of which
member of `F_α(H)` training happened to converge to**. Adding a structural penalty term to the
training LOSS in an attempt to "select a nicer point within the degenerate face" therefore cannot
change a same-H argmin readout's output **at all**, as long as it does not change WHERE the global
minimum of `H` itself sits. If it is strong enough to move the location of the minimum, it is not
degeneracy-breaking — it is training on a **different, combined Hamiltonian** `H' = H + λ·H_struct`,
which is precisely the mandated `H(λ)` continuation architecture already built and measured (Sprint
21 C6/C7: degenerate in raw units, and the specific instance tried — Legacy as `H_struct` added to
AMBER — makes AMBER's landscape measurably worse, not better, `log₁₀E +1.49`). **Corollary: "break
CVaR's degeneracy with a structural term" is not a new lever distinct from Hamiltonian mixing — it is
the same lever under a different name, and the programme has already measured it and found it does
not help under a same-H readout.**

**5.3 Where a structural term COULD matter, and it is not through degeneracy-breaking.** If the
readout Hamiltonian differs from the training Hamiltonian (`H_readout ≠ H_train`, the one non-null
lever per B6/B11/B12), a `best_seen`-style readout is bounded by the **support of the training
distribution**, since it can only rank states that were actually explored. Here a structural term
added to `H_train` *does* matter, but as an **exploration-shaping** device (which states get visited
at all), not as degeneracy-breaking within a fixed optimal face. However, L14/L17 already establish
that **exploration/coverage is not the binding constraint on this instrument** — the exhaustive
enumeration of the entire latent (`2**n` configurations, all `n≤16`) ties a zero-evaluation pool, so
there is no coverage gap for a structural exploration term to close. **Both routes by which a
structural term inside CVaR training could matter are therefore closed by results the programme
already has**: same-H readout makes it provably inert (§5.2); different-H readout makes it an
exploration lever, and exploration is not the constraint (L14/L17).

**5.4 What is left, and it is already adopted.** The only place a "structural term" usefully enters
is not inside the degenerate face at all — it is as the **content of `H_readout` itself**, i.e.
building a better ranking energy to argmin against, which is exactly B12's adopted architecture
instruction (train on the physics, read out with something that ranks). This derivation supplies the
missing "why": every alternative way of using a structural term (regularizing the training loss to
prefer one member of the degenerate face) is provably inert or already-tested-and-null, which leaves
"change what the readout consults" as the only remaining door — not by elimination of imagination,
but by elimination of the other two formally available options.

---

## 6. TYPE-M / TYPE-S EXAGGERATION AS A FUNCTION OF EFFECT/MDE — a lookup table, not a constant

BRIEF §1 already forbids a constant MDE; the same applies to Type-M. Numerically derived (not
asserted) via `z ~ N(θ,1)`, `θ = 2.8016·k` where `k` = true-effect/MDE, significance at the
conventional two-sided 5% threshold (`|z|>1.96`), 4×10⁶ draws per row:

    k = true effect / MDE    power    Type-S      Type-M (|d_hat|/true, given "significant")
    0.3                      0.134    0.0194      2.93x
    0.5                      0.288    0.0013      1.85x
    0.7                      0.500    0.0001      1.41x
    0.9                      0.713    0.0000      1.19x
    1.0  (design power=80%) 0.800    0.0000      1.13x
    1.1                      0.869    0.0000      1.08x
    1.3                      0.954    0.0000      1.03x
    1.5                      0.988    0.0000      1.01x
    2.0                      1.000    0.0000      1.00x

**Reading for this programme's own rule ("0.7-1.3x MDE is the Type-M zone").** Across exactly that
zone, the exaggeration factor itself ranges **1.03x-1.41x**, i.e. a "just significant" result at the
bottom of the zone (`k=0.7`) should be read as **up to 41% inflated**, while one at the top (`k=1.3`)
is close to unbiased (`~3%`). Below the zone (`k=0.3-0.5`, comparisons far under-powered relative to
their own MDE that nonetheless clear a naive threshold — possible via multiple looks, a favourable
sub-analysis, or an unaccounted best-of-K selection per §2.4) the exaggeration is **1.85-2.93x**, and
non-trivial Type-S risk appears (1.9% sign errors at `k=0.3`). **Recommendation, adopted as a
standing rule**: any effect reported with `k = |effect|/MDE` between 0.7 and 1.3 should be labelled
with ITS OWN Type-M multiplier from this table (or a freshly computed one at its exact `k`), not a
blanket caveat; any effect at `k<0.5` reported as "significant" should be treated as a Type-S
candidate and re-examined for a multiple-comparisons or look-elsewhere explanation before being
trusted even in sign.

---

## 7. LITERATURE — searched via WebSearch/WebFetch, 2025-2026 emphasis, errata checked where feasible

### 7.1 CVaR-VQE — the project's Pillar 1 is a KNOWN METHOD applied to a KNOWN COMBINATION, not a new one

Barkoutsos, Nannicini, Robert, Tavernelli & Woerner, **"Improving Variational Quantum Optimization
Using CVaR"**, *Quantum* 4, 256 (2020) — the foundational paper. It already establishes that the
CVaR objective works with only the lowest-energy fraction `ρ` of sampled bitstrings and that this
beats the plain expectation value on combinatorial landscapes, including under realistic
(shot-noise-limited) measurement. The general fact this project's Q6 formalises — that a CVaR-optimal
law is supported on the tail and, for a fully expressive ansatz, coincides with a point mass at the
global minimum — is a direct, essentially immediate corollary of the Rockafellar–Uryasev (2000)
representation of CVaR as a tail-conditional expectation, and is implicit in how Barkoutsos et al.
already describe the objective. **Q6 is a correct and USEFULLY STATED derivation in this project's
own context (§5 above formalises it further), but it is not new mathematics; its genuine
contribution is the operational observation that a REAL deployed VQE's readout is `argmin`/`mode`,
not a tail average, which makes the corollary bite architecturally** (same-H CVaR training cannot
outreach argmin) — that architectural point, and the readout-vs-training-H separation it motivates
(B6/B11/B12), is this project's own and appears to be a genuine, if narrow, contribution.

CVaR-VQE applied specifically to **peptide/protein folding is NOT new as of this campaign**, and this
is the most important novelty correction in this review:

- **QuPepFold** (Roget et al., *PLoS ONE*, Feb 2026; PMC12893577) applies CVaR-optimised VQE to
  6-10-residue peptides (1,224 sequences / 21,600 conformers) on a **tetrahedral lattice** with
  **Miyazawa-Jernigan** contact energies, geometric/chirality constraints, and steric penalties.
  Reports energies only — **no RMSD against experimental structures at all** — and no best-of-N or
  untrained-circuit control is stated in what could be retrieved. CVaR is reported to reach the
  ground state "~30% faster" than plain expectation-value VQE.
- **"Quantum synergy in peptide folding: a comparative study of CVaR-VQE and MD simulation"**
  (*Chemical Physics Letters* or similar, ScienceDirect, 2024) and its companion
  **"A comparative insight into peptide folding with quantum CVaR-VQE algorithm, MD simulations and
  structural alphabet analysis"** (*Quantum Information Processing*, 2024) already benchmark
  CVaR-VQE against classical MD on ~50 seven-residue peptides, predating this campaign by roughly two
  years.

**Consequence for this project's novelty claims.** "Genuine CVaR-VQE as the selector for peptide
structure" (BRIEF Pillar 1) is a **known combination**, publicly demonstrated on peptides multiple
times since 2024, none of it under this project's name. What appears NOT yet published elsewhere,
based on this search: (i) a **continuous torsion-angle** representation combined with CVaR-VQE and a
**physically detailed all-atom force field (AMBER ff14SB/GBn2)** rather than a lattice + coarse
contact potential — the closest external work found (below) uses continuous torsions with VQE but a
different training objective, and the lattice+CVaR papers above use a coarse potential, not AMBER;
(ii) the explicit **readout-Hamiltonian-vs-training-Hamiltonian separation** as a measured,
architecture-level design lever (B6/B11/B12); (iii) the **RMSD-primary, honest-failure** empirical
programme (below) as a sustained, multi-sprint negative-results campaign — most external CVaR-peptide
papers report energy-only metrics and do not test whether the energy landscape ranks structures by
correctness at all.

### 7.2 Continuous-torsion VQE for off-lattice protein structure — a close, concurrent, external result

`arXiv:2609.02113` (already logged, Sprint 21 L6) remains the closest external match: a continuous
torsional-angle VQE for off-lattice structure prediction, dated 2 Sep 2026 — five days before this
document, effectively contemporaneous with this campaign and too recent for errata to exist yet. Its
independent finding that "energy-ranking imbalances persisted across sampled landscapes for ALL
functions" is a third external corroboration of this project's central negative
(`nothing-ranks-within-the-pool` / `the-objective-does-not-rank-the-native`), on different peptides
with different energies, and its own best-of-set-vs-final-model gap (0.623→1.199 Å chignolin,
2.501→3.512 Å Trp-cage) independently reproduces this project's basis-mismatch/selection-gap problem
in one sentence — worth restating here because it means **at least two independent groups, using
different systems, have now hit this project's own two central negatives**, which raises confidence
they are real physics/methodology facts and not artefacts of this codebase.

A third, less directly comparable, external result found in this pass: **"Assessing Cost Hamiltonian
Reliability in Quantum Protein Structure Prediction"** (arXiv:2606.21241, Roget/Damour/Cadet/Wang,
June 2026) is explicitly framed around whether a cost Hamiltonian correctly ranks candidate structures
by correctness — directly on-topic for this project's central finding — but full-text extraction
failed in this pass (only metadata was retrievable); **flagged for a follow-up read, not assessed**,
rather than summarised from a title.

### 7.3 Barren plateaus — status essentially unchanged from Sprint 21's L4, now with a 2025 Nature Reviews Physics treatment

A comprehensive review appeared in *Nature Reviews Physics* (2025) treating barren plateaus as "the
main issue at the moment" in variational quantum computing, and a 2025 *Quantum Information
Processing* systematic literature review taxonomises mitigation strategies into five families
(initialisation, ansatz design, cost-function design, warm-starting/classical pre-training, and
noise-aware strategies). Nothing found in this pass contradicts or supersedes Sprint 21 L4's own
finding (Cerezo et al., arXiv:2312.09121 — barren-plateau-free structure and classical simulability
are the same structure), which remains current best evidence and is **directly reinforced** rather
than merely uncontradicted: this project's own register sizes (9-16 qubits, one qubit per residue)
sit far below where barren plateaus or their absence would be diagnostic of anything quantum at all,
which is the closed half of L4 and is unaffected by anything found here.

### 7.4 Miyazawa-Jernigan — a known, coarse, mostly-one-body potential; relevant to Rule 0's HAMILTONIAN fork

The MJ contact potential (used by QuPepFold above and widely elsewhere in the quantum-folding
literature) is documented as being dominated by one-body (solvation-like) terms rather than genuine
pairwise contact information — a standing critique in the structural-bioinformatics literature
predating this search and not specific to 2024-2026. This is relevant to BRIEF's Rule 0 HAMILTONIAN
axis: any external comparison against an MJ-based competitor is not comparing like-for-like physics
against this project's AMBER ff14SB/GBn2, and that difference should be named explicitly whenever
this project's results are compared to the MJ-based literature above (§7.1), rather than treated as
"another CVaR-VQE peptide result" without qualification.

### 7.5 A claim requiring active scepticism: quantum beating AlphaFold3 on RMSD

A 2026 *Advanced Science* paper (Zhang et al., arXiv:2506.22677) reports a VQE-on-tetrahedral-lattice
framework ("QDockBank") for protein **binding-site** structure on a 127-qubit processor, claimed to
beat AlphaFold2/AlphaFold3 on RMSD and docking affinity. This is an extraordinary claim relative to
everything else in this literature pass (and relative to this project's own multi-sprint experience
that basis/operator mismatches routinely manufacture apparent wins of exactly this size — Sprint 21
L14's `+0.215→+0.025` operator-ladder collapse is the project's own cautionary tale for precisely this
pattern). **Recorded as a claim to be sceptical of, not evaluated in depth here** (full-text
methodology was not retrieved in this pass) — flagged specifically because BRIEF's Rule 0 exists to
catch exactly this shape of claim, and a professor reading this literature section should know it
exists and know why this project does not treat it as settled.

### 7.6 Learning-to-defer / selective prediction / MoE routing — the coordinator's specific request on L2/L3

This search is now less urgent than when it was requested, since §3 above shows L2's "regime
complementarity" is substantially a mechanical artefact rather than a genuine crossover — but the
literature is still worth recording, because L1's ceiling (a real, if unrouted, headroom) is
precisely a **learning-to-defer / model-selection** problem in the field's own terms, and the field's
answer to "can the deferral signal be estimated well enough to matter" is directly relevant to
whether Workstream C should keep investing here.

- **Two-stage learning-to-defer** (Mao, Mohri, Zhong et al., NeurIPS 2023 and follow-ups through 2025,
  e.g. arXiv:2603.14324 "Learning-to-Defer with Expert-Conditional Advice") establishes
  Bayes-consistency and finite-sample excess-risk bounds (via Rademacher complexity) for a learned
  rejector/router — i.e. the field has formal tools for exactly "how much data does the deferral rule
  need," which is the right frame for asking whether this project's 126-target instrument is even
  large enough in principle to learn a router, independent of whether a good feature exists.
- A recurring, explicitly stated limitation across this literature (2024-2025 survey material) is that
  L2D methods require **expert-annotated deferral labels** and are **tightly coupled to the specific
  experts they were trained with** — of direct relevance here, since this project's "experts" (the 13
  arms) are themselves highly correlated (§2.2), which is exactly the regime in which a learned
  router's effective sample size is smaller than its nominal one.
- A 2025 preprint (arXiv:2509.12573, "No Need for 'Learning' to Defer? A Training-Free Deferral
  Framework... through Conformal Prediction") argues a fitted deferral rule can be UNNECESSARY or
  unreliable relative to a conformal-prediction-based, training-free alternative — suggestive that
  this project's own null result (§4, and the length-quartile router in §2) is not an anomaly but
  consistent with a live methodological debate about whether LEARNED deferral rules are worth their
  estimation cost at all in exactly this kind of low-N, high-correlation regime.
- **Mixture-of-experts routing** (2024-2025 LLM literature) independently documents **representation
  collapse**: small external routers "lack sufficient capacity... forcing them to learn their
  prediction via indirect trial-and-error optimisation, which inevitably leads to suboptimal
  routing" — the same qualitative failure mode as this project's null router results, in a completely
  unrelated domain (large language models), which is circumstantial but genuine cross-domain support
  for treating "the router itself is hard to learn, even when the underlying experts are genuinely
  complementary" as a real, general phenomenon rather than a peculiarity of this instrument.

**Net read for the coordinator's question**: the literature does NOT say "the deferral signal is
provably unlearnable here" — it says the field has repeatedly found learned routers hard to get right
in exactly this correlated-experts, limited-data regime, and has produced formal tools (finite-sample
bounds, training-free alternatives) that this project has not yet applied. **Recommendation**: before
building more router infrastructure, compute the Rademacher-style finite-sample bound for a router
over 13 (correlated) arms at n≈100 training targets, and compare it to the 0.48-0.51 Å ceiling — if
the bound already exceeds the ceiling, the field's own theory says not to expect a learned router to
help here, independent of any further empirical search. Not computed in this pass; registered as the
natural next derivation.

### 7.7 Summary novelty table

| claim | classification |
|---|---|
| CVaR-VQE as a selector, in general | **known method** (Barkoutsos et al. 2020) |
| CVaR-VQE applied to peptide/protein folding | **known combination** (QuPepFold 2026; two 2024 papers) — NOT novel, contra any implicit assumption otherwise |
| Continuous-torsion representation for VQE structure prediction | **known combination**, concurrent (arXiv:2609.02113, 5 days old) |
| Continuous torsion + CVaR-VQE + all-atom AMBER ff14SB/GBn2 jointly | **plausible implementation novelty** — not found combined elsewhere in this pass, but the search was not exhaustive |
| Q6 (CVaR face supported on tail; pool-restricted argmin equivalence) | **not new math** (immediate corollary of Rockafellar-Uryasev + Barkoutsos); **its architectural application here (readout-vs-training-H separation, B6/B11/B12) is this project's own** |
| The readout-H / training-H separation as a measured design lever | **plausible empirical + implementation novelty** — not found stated this way elsewhere |
| "Energy does not rank structures / optimising harder can hurt" | **empirical finding with growing external corroboration** — now 3 independent groups (this project, arXiv:2609.02113, and by implication the still-unread arXiv:2606.21241) — increasingly looks like a real phenomenon, not a codebase artefact, but is not this project's alone to claim |
| The routing-ceiling / discrimination-not-search framing | **not found stated this way in the quantum-folding literature searched**; closest analogue is the learning-to-defer field's general framing, applied here to a new domain — **plausible empirical novelty**, contingent on the Rademacher-bound check in §7.6 not closing it first |

**Scope limit on this whole section, stated rather than hidden**: this is one focused literature
pass (roughly a dozen targeted searches and two full-text fetches), not a systematic review. Absence
of a match in these searches is evidence, not proof, that a given combination is unpublished —
several of the closest matches (QuPepFold, the 2024 peptide-CVaR papers) were found only because the
coordinator specifically asked about CVaR-peptide combinations; a claim of novelty for any item in the
table above should be re-checked before it appears in a final report exactly as thoroughly as it was
checked here, not assumed to remain true.

---

## 7.8 ADDENDUM (same session) — L2 WITHDRAWN BY THE COORDINATOR; THE ARM-LIST DISCREPANCY RESOLVED; A NEW LOAD-BEARING POSITIVE (mreal.py) AUDITED AND SURVIVES; A NEW LOAD-BEARING NEGATIVE (routercv.json) FOUND INDEPENDENTLY

**L2 formally withdrawn.** The coordinator accepted §3's decomposition in full, including that his
own attempted defence (re-stratifying on the ORACLE-independent `pool_oracle` rather than the
incumbent's realised RMSD) fails for the same algebraic reason:
`Cov(stratifier, lat) - Cov(stratifier, inc)` is negative whenever the latent is less
difficulty-sensitive than the incumbent, for ANY stratifier correlated with true difficulty, not only
the incumbent's own realisation. His "9.6-SE durable interaction" is withdrawn along with the
original claim. What survives is exactly §3's restatement: a single global slope fact, not a
complementary-regime router target.

**The L1 arm-list discrepancy (§1) is resolved** — see the dated addendum in `s22/PREREG_D.md`
D-AUDIT-1. Bookkeeping error on the coordinator's side, not a different construction; his 0.482 is
now independently confirmed to five significant figures via `s22/results/routerdata.json`'s
`A_plus_latent` configuration (an artefact from a concurrent Workstream C run that appeared in
`s22/results/` mid-session).

**A new positive claim, audited: `s22/mreal.py`, "65% of the m-ladder headroom transfers across
split halves."** Reproduced directly from `s22/results/mreal.json`
(`s22/results/d_mreal_audit.py`): primary `-0.2394 [SE 0.0353, MDE 0.0990]`, matching the
coordinator's report. Four checks:

1. **Fold-clustered CI, which the original `report()` does not compute** (a real, if minor, process
   gap against this project's own standing practice — e.g. `d_cvarop.py` always prints both):
   `[-0.292,-0.180]`, essentially identical to the i.i.d. CI `[-0.312,-0.173]`. Not fragile to
   clustering.
2. **Jackknife** on the primary mean: leave-one-out range `[-0.244,-0.228]` against a full mean of
   `-0.239`. Not driven by a handful of targets.
3. **The coordinator's own named weakness, tested directly**: the half-pool `m=500` rung is
   degenerate with "the whole 250-member half" (`min(500,250)=250` in source). Reran with that rung
   dropped from the selectable ladder, using the IDENTICAL RNG stream/seeding
   (`s22/results/d_mreal_notrunc.py`, so the same 8×126 half-splits are reproduced exactly, not
   resampled): `-0.1958 [SE 0.0325, MDE 0.0912]`, W/L 77/39. **82% of the original effect survives**
   the removal of the flagged weak point, still 2.15× its own MDE. The truncation is not carrying
   the result.
4. **A downgrade, not a refutation**: "median tied rungs 1.00 of 6" is computed via an EXACT
   floating-point tie test (`≤ min+1e-12`) on a continuous outcome (a coordinate-average RMSD), which
   will read ≈1 (no ties) essentially by construction regardless of whether the underlying ladder is
   statistically flat. This is a different and much weaker question than Sprint 21 L23's "rungs
   within 1 SE of the target's own minimum" (median 2.0/6) — the two numbers are not in tension and
   `mreal.py`'s statistic should not be read as "sharper than L23."

**Verdict: `mreal.py` survives every attack run against it, including the one the coordinator
pre-identified as its own weakest point.** The per-target optimal averaging breadth is a real,
stable property that transfers across two independent random draws from the same target's pool — a
genuine positive control, distinct from (and not undermined by) the router-estimation failures below.

**A new negative claim, found independently: `s22/results/routercv.json` and `routercv2.json`
(Workstream C, in progress, artefacts appeared mid-session).** Five achieved-router constructions,
using native-free features far richer than anything tested in §2/§4 of this document (`rg_z`,
`rg_gap`, `score_mean/sd/skew/gap`, `pool_spread`, `sim_mean_pool/top75`), 5-fold held out, tercile
binning. **Every single one is WORSE than the fixed incumbent**: `m_ladder_only` `+0.048 [MDE
0.081]`, `pool_all` `+0.059 [MDE 0.098]`, `A_plus_blend` `+0.032 [MDE 0.104]`, `ALL` (16 arms)
`+0.043 [MDE 0.104]` — all NOT MEASURED but uniformly wrong-signed (capture fractions -7.6% to
-15.4%) — and **`A_plus_latent`, the exact configuration that reproduces the headline 0.482 Å
ceiling, comes in at `+0.070 [+0.0048,+0.1459]`, a CI EXCLUDING ZERO on the HARMFUL side**: this
achieved router is significantly worse than doing nothing. `routercv2.json`'s two follow-up variants
(`in_fold` `-0.046 [MDE 0.069]`, `gated` `+0.028 [MDE 0.070]`) are both NOT MEASURED. This
independently reproduces, with a far richer feature set, the same wrong-signed null this document's
own §2/§4 length- and rg-based routers already found, and adds a genuinely new, sharper result: at
least one specific, well-defined construction is SIGNIFICANTLY harmful, not merely unhelpful.

**This needs a specific methodological check before it is read as "no native-free signal exists,"
per §2.4's multiple-comparisons rule**: was the tercile-binning procedure itself fit by trying
several bin counts/feature subsets before landing on these five, and if so, has that search's own
best-of-K exposure been accounted for? A router-construction search fit on small per-bin training
samples, choosing among near-tied arms (recall the m-ladder's six rungs are known-flat, L23), is
exactly the regime where fitting noise costs more than the ceiling can pay back — which would show up
as a uniformly wrong-signed, occasionally-significant-in-the-harmful-direction pattern, precisely
what is observed. **Recommended next step for Workstream C, not yet done here**: report how many
constructions preceded these five, and re-score exactly one pre-registered construction against a
fold that played no part in choosing it.

**The combined diagnosis, now resting on two independently-audited artefacts from two different
lanes**: the per-target routing opportunity is REAL and STABLE (mreal.py, survives attack) and, so
far, UNFINDABLE by every native-free feature construction tried, including a rich one (routercv.json,
independently discovered and audited here) — the textbook learning-to-defer shape named in §7.6.

## 9. THE FINITE-SAMPLE BOUND FOR ROUTER LEARNABILITY AT n≈126 OVER ~13 CORRELATED ARMS

**Question posed.** Is a native-free router learnable AT ALL at this project's sample size, or is `n`
itself the binding constraint — independent of which features are tried?

### 9.1 The bound

Standard uniform-convergence theory (Hoeffding + a union bound over a finite/discretised hypothesis
class `F`, tightened with the empirical-Bernstein variance term) gives: with probability `≥ 1-δ`,
for every router `f ∈ F` trained on `n` targets,

    R(f) ≤ R̂(f) + sqrt( 2σ² (log|F| + log(1/δ)) / n ) + (2M(log|F| + log(1/δ))) / (3n)

where `σ²` is the per-target variance of the achieved-minus-fixed outcome, `M` is the full range of a
single target's outcome (the Hoeffding/Bernstein range term), and `|F|` is the size of the
(discretised) class of candidate routers being searched over — this is Massart's lemma applied to a
finite class, which is the right instrument here because every router construction actually tried in
this campaign (tercile bins, a fixed feature list, a fitted Ridge/RF with a bounded parameter norm) is
equivalent, for generalisation purposes, to searching a finite or VC-bounded class, and `log|F|` (or
the VC-dimension `d`, which enters the analogous bound in the same way up to log factors) is exactly
the SAME complexity term that governs the multiple-comparisons correction derived in §2.4 — this
section is that derivation made quantitative rather than qualitative.

### 9.2 Plugging in this project's own measured numbers

`σ ≈ 0.41 Å`, measured TWICE independently and in agreement: `sqrt(n)·SE` from `mreal.json`
(`0.0353·√126 = 0.396`) and from `routercv.json`'s `ALL` construction (`0.0371·√126 = 0.417`). `M = 8`
Å, a conservative round envelope for a single target's worst-case outcome (individual arms in this
project's own tables reach 6-8 Å; this is a modelling choice, stated as such, not a measured
constant — see §9.4 for its effect on the conclusion). `δ = 0.2`, matching this project's own 80%
design-power convention. `n = 100`, approximating one training fold of the 5-fold CV this project
uses at `n=126`.

    router class                                          log|F| or VC-d   gap at n=100   n needed for gap=0.24   n needed for gap=0.48
    1 global threshold, K=13 arms                                 3.3           0.39                 180                    78
    3-bin tercile, 1 feature, K=13 (routercv's own design)        12.3          0.96                 513                   222
    5-bin quintile, 1 feature, K=13                               22.0          1.54                 874                   377
    linear router, p=5 features                                    6.0          0.57                 281                   121
    linear router, p=15 features (routercv's full feature list)   16.0          1.18                 650                   281
    shallow tree, depth 3, over 15 features                       26.6          1.82                1045                   451

(`s22/agentD_FINDINGS.md` derivation script not separately saved — this is a closed-form calculation,
reproducible from the formula in §9.1 and the numbers in this paragraph; anyone can re-derive it in
one line of Python.)

### 9.3 Reading the table against what actually happened

At the training-fold size this project actually has (`n≈100`), **even the single simplest possible
router (one global threshold choosing between two candidate arms) already has a generalisation gap of
0.39 Å — comparable to the ENTIRE 0.482 Å ceiling.** To reliably resolve even HALF the ceiling
(`ε=0.24`) with that simplest class needs `n≈180` — more targets than exist in the WHOLE 126-target
benchmark, let alone one training fold. Every router construction actually tried in this campaign sits
at or above the "linear, p=5-15 features" or "tercile bin" rows, which require `n` in the
**280-1045** range to resolve even half the ceiling — **1.3× to 3.6× the entire benchmark's size**.
**This is not a loose analogy to what was observed — routercv.json's Ridge/RF over 15 features is
almost exactly the "linear router, p=15" row, and its measured MDE (`0.098-0.104` for a SINGLE
pre-specified construction) already sits close to the WHOLE ceiling; the bound above says that once
you account for having SEARCHED over several such constructions (as every lane did), the effective
bar is 2-4× higher still.**

### 9.4 Honest limits on what this proves, stated before it is over-read

- This is a **distribution-free, worst-case** bound. Real learning algorithms with implicit
  regularisation, cross-validated model selection, or favourable structure in the true function
  routinely beat such bounds in practice. **The correct reading is "no proof of learnability exists at
  this n for anything beyond the simplest router," not "a router is provably impossible."** The
  distinction matters and must be kept in any final report.
- `M=8` is a modelling choice. Using the Rademacher/variance term ALONE (dropping the Bernstein range
  term, which dominates at small `n`) still gives `gap ≈ sqrt(2σ² log|F| / n)`: at `n=100`,
  `log|F|=12.3` (the tercile-bin row), this alone is `0.31` — still on the same order as the whole
  ceiling, so the conclusion is not an artefact of the `M=8` choice; it survives dropping the more
  arbitrary term entirely.
- Targets are **not i.i.d.** — they share fold/family structure (this project's own 5-fold clustering
  exists precisely because of this), so the EFFECTIVE `n` for a Rademacher-style bound is smaller than
  the nominal one, meaning this section's numbers are, if anything, **optimistic** (favourable to
  learnability) relative to reality.
- Arm correlation (target difficulty as a dominant shared factor, §2.2) does not enter this bound
  directly — it affects the ceiling's own size and the achievable ORACLE, not `|F|`, since `|F|` here
  indexes candidate ROUTING RULES, not candidate arms. A more refined analysis could shrink `|F|` by
  the arms' effective (rather than nominal) count, which would make the bound somewhat more
  favourable to learnability; not attempted here, flagged as a refinement.

### 9.5 Disposition

**Sample size, not feature quality, is the more defensible primary explanation for the four
independent router failures** (this document's length- and rg_z-routers, and Workstream C's Ridge,
RF, and candidate-geometry constructions). At the project's actual training-fold size, standard
learning theory already predicts a generalisation gap comparable to or exceeding the entire ceiling
for any router more complex than a single global threshold — and a single global threshold is
precisely what `s22/mreal2.py` already tested and found captures `+0.002`, nothing. **This converts
"we tried several routers and none worked" into a structural statement**: at this `n`, with this
`σ`, no router of realistic complexity carries a learning-theoretic guarantee of capturing the
ceiling, and the four empirical failures are consistent with, rather than surprising relative to,
that guarantee's absence. **This does not mean routing is closed forever** — it means the lever
available to open it is more data (more cluster-disjoint targets) or a variance-reduction device
(a much lower-noise feature, or a hierarchical/shrinkage estimator that trades bias for variance
explicitly, rather than a flat per-bin or per-leaf estimate) — not another pass through the same
class of router constructions at the same `n`.

---

## 10. arXiv:2606.21241, RECOVERED — A THIRD, VERY CLOSE, INDEPENDENT CORROBORATION, WITH ONE POINT OF GENUINE TENSION

Full text recovered via the HTML mirror (`arxiv.org/html/2606.21241v1`) after the PDF route failed.
**"Assessing Cost Hamiltonian Reliability in Quantum Protein Structure Prediction"** (Roget, Damour,
Cadet, Wang; June 2026).

**What it does.** Lattice-based (tetrahedral, Cα-only) quantum protein structure prediction with a
**Miyazawa-Jernigan contact-energy** cost Hamiltonian. Asks exactly this project's own question:
does the cost Hamiltonian's landscape correlate with structural correctness (RMSD)?

**What it finds, quoted directly.** *"For small peptides and on average, the energy landscape of the
considered cost Hamiltonian is not correlated well enough to the actual error to provide meaningful
predictions."* More specifically, for peptides **≤15 amino acids** (almost exactly this project's own
9-16-residue instrument): *"the RMSD of the solution with minimum cost is much higher at around 5 Å"*
against an achievable ~1.5 Å, and — the sentence that matters most — ***"the average RMSD across all
solutions is, on average across all peptides, better than the optimal solution according to the cost
Hamiltonian."*** Spearman correlation between cost and RMSD-error is reported **negative for small
peptides, null at length 50, and positive only for larger proteins**.

**This is a very close match to this project's central negative, independently obtained, on a
different energy (MJ contact potential, not AMBER/Legacy), a different representation (discrete
tetrahedral lattice, not continuous torsions), and a different group.** "The cost-optimal structure
is worse than the average of all candidates" is this project's own `C10`/`search-saturates-
discrimination-binds` shape (both physics energies worse than a matched random tail; optimising the
objective can make things worse, not better) arriving from a completely independent implementation.
Combined with `arXiv:2609.02113` (Sprint 21 L6, continuous-torsion, independently found "energy-
ranking imbalances... for ALL functions"), **this project's central finding — that a coarse or
simplified cost Hamiltonian for small-peptide structure prediction fails to rank candidates by
correctness, sometimes anti-ranking them — now has THREE independent external sightings**, spanning
lattice and continuous representations and at least two different energy families. **This changes the
finding's status from "a property of this codebase" to "a property of the class of methods,"
which is the stronger and more publishable claim, and it should be reported as such.**

**The point of genuine tension, not smoothed over.** The paper also reports that correlation
**IMPROVES with more interaction shells / larger instances** — i.e., it frames the failure as (at
least partly) a **resolution** problem: a coarser, fewer-shell contact potential ranks worse, and a
richer one (or a larger, more cooperative protein) ranks better. **This project's own result does not
follow that pattern**: AMBER (the far richer, all-atom, physically detailed energy) is not shown to
rank better than Legacy (the coarser compactness model) — Sprint 21's own matrix (C8, C10) has AMBER
performing WORSE than Legacy in most configurations, and both worse than a matched random tail. If
"more physical detail improves ranking" were the operative variable, AMBER should have been
the one energy in this project that ranks well, and it is not. **This is worth stating plainly rather
than quietly dropping**: either (a) "interaction shells" in a lattice contact potential is not the
same axis as "all-atom force-field detail" in a continuous model, and the two findings are
compatible but about different variables, or (b) there is a real disagreement between this project's
result and the external paper's proposed resolution-dependent mechanism, and only one of them can be
the operative explanation for THIS project's own AMBER-vs-Legacy ordering. **Recommend this tension be
named explicitly in any final report rather than cited as pure corroboration** — the corroboration is
real and strong on the headline (energy doesn't rank small-peptide structures), and it is NOT
corroboration on the proposed mechanism (resolution/detail fixes it), where this project's own AMBER
result is, if anything, evidence against that mechanism.

---

## 11. LEGACY/AMBER ANTI-CORRELATION ON THE COMPACTNESS AXIS (Workstream B's perturbation probe) — NOVELTY CHECK

**Claim under check.** B's controlled perturbation probe found Legacy and AMBER **anti-correlated on
a compactness axis at matched magnitude** (Spearman `+0.454` vs `-0.444`), with the STERIC response
shared rather than specialised between the two energies.

**Search conducted.** Targeted searches for knowledge-based-vs-physics-based potential disagreement,
contact-potential-vs-physics-based anti-correlation, and compactness-axis-specific perturbation
studies. **No paper matching this specific, controlled-perturbation result was found.** What IS well
documented, and predates this by roughly two decades, is the broader observation that knowledge-based
(statistical/contact) and physics-based (force-field) potentials for protein structure are built on
different, sometimes disputed, physical premises and are known to disagree in general (the
"knowledge-based potentials: physical or statistical?" debate — Ben-Naim and others have long
questioned the physical interpretability of potentials of mean force derived from structural
statistics, as distinct from genuine physics-based energies). **What was not found anywhere in this
pass is the SPECIFIC, sharper claim**: that the two disagree specifically and oppositely on a
COMPACTNESS axis, AT MATCHED MAGNITUDE, under a CONTROLLED PERTURBATION (rather than an observational
correlation across a decoy set), while SHARING their steric response.

**Classification: plausible empirical/mechanistic novelty.** The general phenomenon class
(knowledge-based vs. physics-based potential disagreement) is old and known; the SPECIFIC mechanism
(matched-magnitude anti-correlation on compactness, shared steric channel, obtained by perturbation
rather than correlation) is not found in this search and — pending a more exhaustive check than one
literature agent can give in one pass — is plausibly a genuine, useful, and reasonably narrow
empirical/mechanistic contribution. **Caveat identical to §7.7's scope limit**: absence from this
search is evidence, not proof, of novelty; the standard search terms for this specific mechanism are
not obvious (it may be filed under docking, decoy discrimination, or force-field validation
literature using different vocabulary), and a more determined search — ideally by someone who can
also check the CASP-era decoy-discrimination literature directly — would be needed before this claim
is put in a final report as confirmed-novel rather than plausibly-novel.

---

## 12. OUTSTANDING / NOT YET DONE, STATED PLAINLY

- ~~The Rademacher/finite-sample bound~~ — **DONE, §9.**
- ~~`arXiv:2606.21241` full text~~ — **DONE, §10.**
- ~~Legacy/AMBER compactness anti-correlation novelty check~~ — **DONE, §11.**
- FAIL18 membership was not recovered from available artefacts (§3); the 72%/24% win-rate claim is
  neither confirmed nor refuted here.
- A continuous/regression-form router for `rg_z`/`rg_gap` (§4) — the binned router tested is a weak
  instrument for a ρ≈0.35 signal; §9's bound now suggests this is LOW PRIORITY relative to the sample-
  size constraint (a regression router over even 2-3 features already sits in the "n needed: 121-281"
  rows of §9.2's table, well beyond what exists) — attempt only if the coordinator still wants it after
  reading §9, and pre-register it first per §2.4's rule about multiple router constructions.
- §9's refinement noted but not attempted: shrinking `|F|`'s effective size using the arms'
  correlation structure (rather than their nominal count) could tighten the bound; not done here.
- §11's search was one pass with fairly narrow search terms; a more exhaustive check (CASP-era decoy
  discrimination literature specifically) would strengthen or weaken the novelty classification.
