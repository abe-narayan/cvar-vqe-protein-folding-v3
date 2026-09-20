# S30 LANE L, TOPIC 3 -- SET-FUNCTION OPTIMISATION AND CVaR ON A SET

Brief priorities 3 and 4. The starting fact is S29-L25, which I re-read in full before reading
anything external:

> the f-optimal 2-subset is NOT the energy-order prefix on 11/12 targets (exhaustive over all
> 124,750 pairs), the f-optimal m=5 subset is non-prefix on 12/12, mean objective gap +0.152,
> **and the escape buys nothing**: those subsets' averages lose to the prefix on 10/12 and to
> production on 8/12.

So a genuine combinatorial problem exists. The brief asks whether the literature has machinery for
it. **My answer is that the machinery does not apply and, more usefully, is not needed: the problem
is computationally easy and informationally expensive, and this project has been treating it as the
reverse.**

---

## 1. THE OBJECTIVE FACTORS THROUGH THE CENTROID, AND THAT SETTLES MOST OF IT

From S29-L25's own definition:

        V(S) = f( mean_{i in S} W_i )

**V depends on S only through one point of R^(3n) -- the centroid of the selected members.** That
single observation disposes of most of the reading list:

**Submodularity and its guarantees: DO NOT APPLY, and the reason is a stated hypothesis, not a
technicality.** Nemhauser-Wolsey-Fisher's `1 - 1/e` and Das & Kempe's weakly-submodular extension
`1 - e^(-gamma)` (submodularity ratio) both require the set function to be **monotone**. S29-L25
states in its own words that V is *"neither additive nor monotone"*. Every guarantee in that family
is void here. This is a family-level rejection with the violated assumption named: **monotonicity,
not submodularity, is what we fail first.**

**DPPs, facility location, diversity-aware selection: DO NOT APPLY.** These parameterise a
*preference over sets* (repulsion, coverage). V has no such preference -- two sets with the same
centroid have identical value, so V is constant on enormous equivalence classes and carries no
diversity structure at all. Separately, Abe et al. (arXiv:2302.00704) is already in the S29 index
as the negative result that diversity interventions harm good ensembles.

**QUBO / Ising reductions: AVAILABLE BUT POINTLESS.** If `f` were quadratic in the centroid, then
`V(S) = (1/|S|^2) sum_{i,j in S} Q_ij` exactly -- the densest-k-subgraph form, NP-hard at fixed
`|S|`, poly-time by max-flow when `|S|` is free (Goldberg 1984). That is a real structural
observation and it says the *fixed-m* formulation is the hard one. But it is moot, because:

**The continuous relaxation is tight and cheap.** The reachable centroids are
`{ sum_i lambda_i W_i : lambda in the 1/m-grid on the simplex }`. By Maurey's empirical method,
for any hull point `c` there is an `m`-multiset whose centroid is within `R/sqrt(m)` of it, with
`R = max_i ||W_i - c||`. Using the project's own dispersion,
`sqrt(mean_k |d_k|^2) = sqrt(63.82) = 7.99` in the summed-coordinate norm of
`s23/results/errdecomp.json` (per-residue RMSD-equivalent 1.92 A), the grid is fine at deployed
`m` and coarse only at the small `m` S29-L25 actually searched:

```
   m       R/sqrt(m), RMSD units   (R taken as the RMS dispersion; a max-based R is ~2-3x larger)
   75              0.22 A
    5              0.86 A
    2              1.36 A
```

**So: minimise `f` over the simplex by Frank-Wolfe or projected gradient -- a convex-domain
problem, seconds -- and round. That is the correct classical counterpart (S29 contract rule 15)
for any quantum arm on this objective, and it strictly upper-bounds what a circuit could find.**

## 2. THE REFRAME THAT I THINK MATTERS MOST

Put section 1 beside two numbers the project already owns (S29 section 12.1, section 0 item 3):

```
  2 members with ORACLE WEIGHTS                      1.4315 A
 75 members with ORACLE MEMBERSHIP (uniform weights) 2.3055 A
 choosing 2 of 500                                  ~17.9 bits
 choosing the argmin of the top-128                   7   bits
```

The expressiveness lives in the **weights**, not in the membership; and the weighted object costs
*more* bits, not fewer. Therefore:

> **The set-selection problem here is computationally easy and informationally expensive.**
> Solving it better is worth nothing -- an exhaustive search already solved it (S29-L25) and the
> answer was worse than production. What is scarce is the information needed to *specify* a good
> set, and no solver supplies information.

This is the single sentence I would put in front of any future proposal to bring optimisation
machinery -- submodular, quantum, or otherwise -- to a selection problem in this pipeline. It also
explains S29-L25's second clause mechanically: escaping the set-equality theorem creates no
information, so the optimum of a marginal-class objective moves the support without moving the
answer.

## 3. CVaR ON A SET FUNCTION -- THE LITERATURE EXISTS, AND IT SAYS SOMETHING SHARP

The brief asked whether anyone has put CVaR on a *set function* rather than a sample. They have.

- **Maehara T, Oper Res Lett 43:526-529 (2015), "Risk averse submodular utility maximization".**
  The CVaR of a stochastic submodular set function **is not itself submodular**, and -- the result
  worth having -- **no polynomial-time algorithm attains any multiplicative approximation to the
  optimal CVaR** (under stated assumptions on the risk level), unless P = NP. Adding CVaR to a set
  function is not a decoration; it changes the complexity class of the problem.
- **Wilder B, AAAI 2018, "Risk-Sensitive Submodular Optimization".** The escape from Maehara's
  hardness is to stop asking for a single set: relax to a **portfolio (a distribution over sets)**.
  Wilder gives a black-box reduction from the discrete portfolio problem to CVaR maximisation of a
  *continuous* DR-submodular function over a down-closed convex set, with a `1 - 1/e` guarantee.
  Ohsaka & Yoshida (2017) reached the same relaxation independently for influence maximisation.

**The structural reading, and I am flagging the mismatch rather than hiding it.** The pattern is
the same as section 1: *the single-set formulation is the hard one; the fractional/portfolio
formulation is the tractable one.* A quantum state over candidate-index bitstrings **is** a
portfolio over sets, and the readout collapses it to one set. That resonance is suggestive and I
want to be honest that it is only that -- Maehara and Wilder take CVaR over **exogenous** randomness
in the objective, whereas our CVaR is over a distribution the optimiser itself controls. The
theorems do not transfer as stated. What transfers is the design lesson, and it is the same one
section 2 reaches from the other side: **relax the set, do not search it harder.**

## 4. A HYPOTHESIS OF MINE, TESTED AND REFUTED -- FINITE-SHOT CVaR BIAS

**The hypothesis.** The empirical CVaR is the mean of the `k = ceil(alpha * shots)` lowest of
`shots` draws. Estimators of tail means are biased at finite sample, the bias depends on the shape
of the distribution near the VaR cut, and the optimiser minimises *the estimate*. So the VQE should
have a shot-noise-driven incentive to shape its distribution in ways that flatter the estimator --
which would be a mechanism for the project's most repeated phenomenon, "the objective is optimised
but the structure does not move", and for S29's measured "the deployed arm is 0.0097 A worse than a
fixed 75-prefix because its realised m wanders".

**The check** (`s30/lit/s30_L_cvar_bias.py`), at the deployed cell's own parameters -- `n = 9`
(512 basis states), `alpha = 0.18`, `shots = 2048` from `core/quantum.py:179`, tail `k = 369`:

```
A. bias vs concentration (Boltzmann states)      bias ranges +0.0000 .. +0.0025, sd 0.017-0.036
B. states MATCHED on true CVaR, support 2..32    bias +0.0001 .. +0.0019, NO trend with support
C. does the bias move the argmin?                exact-minimising beta 2.5; finite-shot 2.5;
                                                 penalty +0.0000
D. Barkoutsos's flat minimiser set, finite shots  support 92 -> 1: bias -0.0006 -> +0.0162,
                                                 i.e. concentration is PENALISED, not rewarded,
                                                 and by less than one sd of the shot noise
```

**REFUTED, by my own check, in the direction opposite to my guess.** The empirical lower-tail CVaR
is *pessimistically* biased (verified analytically on a two-sample case: `E[min(X1,X2)] = -0.564`
against `E[X | X <= median] = -0.798` for a standard normal), the bias at 369 tail shots is ~0.1%
of the objective's range, it does not order states by concentration, and it does not move the
argmin. **Finite-shot CVaR estimator bias is not a mechanism for anything at this project's shot
count. Do not spend on it.**

Scope of the refutation, so it is not over-read: the bias grows as `alpha * shots` falls. At the
`foldvqe` sort path (384 shots, `alpha` 0.18 -> 69 tail shots) it is ~5x larger and still an order
of magnitude below the effect sizes this project cares about. **If any lane proposes a low-alpha or
low-shot arm, this closure lapses and must be re-run.**

What does survive from this reading is Barkoutsos's own falsifier, already in the S29 index and
worth repeating because section 4D is consistent with it: CVaR's global-minimiser set is
`{theta : overlap >= alpha}`, **large and flat**, so any accuracy change attributed to CVaR
optimisation must be shown not to be a tie-break inside that set. Finite shots do not break the
degeneracy usefully -- test D's degeneracy-breaking (+0.016) is **smaller than one standard
deviation of the shot noise (0.027)**, so it is not even reliably visible, let alone exploitable.

## 5. WHAT I REJECTED IN THIS TOPIC

| source / family | why rejected |
|---|---|
| Nemhauser-Wolsey-Fisher `1 - 1/e`; Das & Kempe submodularity ratio; randomised greedy for weakly submodular functions | all require **monotone**; V is stated non-monotone in S29-L25 |
| Buchbinder et al. double-greedy for unconstrained non-monotone submodular | V is not submodular either; it factors through a centroid |
| determinantal point processes; facility location; diversity-aware subset selection | V is constant on centroid-equivalence classes; no diversity structure exists to exploit |
| QUBO / Ising reductions of subset selection | available; moot, because the continuous relaxation is tight and cheap (section 1) |
| densest-k-subgraph hardness (Bhaskara et al.) and max-density-subgraph tractability (Goldberg 1984) | NOTED as the correct structural analogy -- fixed-m is the hard formulation, free-m is the easy one -- but not actionable, for the same reason |
| Maehara 2015; Wilder 2018; Ohsaka & Yoshida 2017 | KEPT as the design lesson (relax the set, do not search it); NOT transferred as theorems, because their CVaR is over exogenous randomness and ours is over a distribution the optimiser controls |

## 6. WHERE I COULD BE WRONG

- Section 1's Maurey bound is for **multisets** (weights on a `1/m` grid), not for the uniform-weight
  **subsets** the deployed readout uses. The subset family is strictly poorer, and the project's own
  `1.4315` vs `2.3055` gap is the size of that difference. My "the relaxation is tight" claim is
  therefore about the *weighted* problem, and I have stated it that way -- but a reader could take
  it as a claim about subsets, and it is not.
- `R` in the Maurey bound should be a **max** over members, and I used the RMS dispersion, which
  understates it by perhaps 2-3x. The table's numbers are therefore optimistic by that factor;
  the ordering across `m` is unaffected.
- Section 4's refutation is on a **synthetic Gaussian spectrum**, not on a real cell's energies. A
  real spectrum with heavy degeneracy or large gaps could behave differently. I judged the margin
  (0.1% of range, argmin unmoved) large enough that a real-spectrum re-run is not worth a lane's
  time, but that is a judgement and it is mine.
