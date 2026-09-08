# SPRINT 22 — DECISION LEDGER

## L1 — THE ROUTING CEILING: 0.482 Å, AND IT IS NOT CONCENTRATED

Measured before briefing any lane, over **13 deployable native-free arms** (the m-ladder
500/150/75/20/5/1, consensus medoids whole-pool and top-75, four latent arms), n=126:

    best single fixed arm (pool top-75 average)   3.048
    ORACLE per-target routing                     2.566
    HEADROOM                                      0.482 A

**Every arm wins on some target** (winner counts 1–17, none dominant). **The headroom is broadly
distributed, not a few outliers**: 50% of it comes from the top 21 targets (17% of the panel), 80%
from the top 50 (40%), median per-target headroom **0.317 Å**, and only **12/126** targets have the
incumbent already best. *A router does not need to find a needle.*

`corr(headroom, incumbent RMSD) = +0.564` — **the headroom is largest where the incumbent is worst.**

> **A NULL I GOT WRONG, RECORDED SO NO LANE REPEATS IT.** I ran a permutation null shuffling targets
> *within* each arm. It returned **1.115 — better than the observed oracle.** It is **mis-specified**:
> arms are strongly correlated through **target difficulty**, and shuffling destroys that shared
> factor, letting the min draw from easy targets. The oracle routing value is a legitimate ceiling
> (like `latent_oracle`); the min-of-K concern applies to **estimating a router's performance**, which
> is why any router must be scored on **held-out folds**.

## L2 — THE LATENT AND THE POOL ARE COMPLEMENTARY, AND THE CROSSOVER IS DIFFICULTY

L20 closed the generative source on an **aggregate**: `lat_avg75 - incumbent = +0.390`. **That
aggregate hides a monotone interaction with a SIGN FLIP.** Stratified by true difficulty (n=126):

| quartile | incumbent | lat_avg75 − incumbent | W/L |
|---|---|---|---|
| easiest | 1.17 | **+1.469** [+1.224,+1.699] | 0W/32L |
| Q2 | 2.40 | +0.448 [+0.232,+0.677] | 8W/23L |
| Q3 | 3.28 | +0.152 [−0.045,+0.372] | 17W/14L |
| **HARDEST** | 5.33 | **−0.455 [−0.717,−0.217]** | **22W/10L** |

**Monotone across four quartiles, sign-flipping, and significant at both ends.** On the known FAIL18
set a **latent arm wins 13 of 18 (72%) against 24% elsewhere**; a latent arm wins **66%** of the
hardest quartile against **6%** of the easiest.

**Retrieval is better where retrieval works; the generative latent gets relatively better exactly
where retrieval struggles.**

> **CORRECTED WITHIN THE HOUR BY MY OWN OBJECTION, AND IT PARTIALLY FIRES.** The stratifier above is
> **the incumbent's own RMSD**, so the easiest quartile is selected on the incumbent having been
> good — **regression to the mean can manufacture part of the gradient.** I raised that objection
> against myself and tested it by re-stratifying on **`pool_oracle`**, which is still ORACLE but is
> **independent of the incumbent's realisation**:
>
>     stratifier              easiest      Q2         Q3         HARDEST
>     incumbent RMSD          +1.469     +0.448     +0.152     -0.455 [-0.704,-0.203]  SIG
>     pool_oracle             +1.499     +0.262     -0.002     -0.155 [-0.388,+0.062]  n.s.
>
> **THE INTERACTION SURVIVES AND IS OVERWHELMING.** Continuous, unbinned, against the independent
> stratifier: `d = a + b·z(pool_oracle)` gives **b = −0.625, SE 0.065 — 9.6 SE.** The latent gains
> relative to retrieval as targets get harder, and that is not an artefact.
>
> **THE SIGN FLIP DOES NOT SURVIVE.** *"The latent significantly WINS on the hardest quartile"* was
> **partly regression to the mean** and is **withdrawn**: on the independent stratifier it is
> **−0.155 [−0.388,+0.062]**, a null. The intercept says the latent is still **+0.405 worse at mean
> difficulty**.
>
> **What is established: a strong monotone INTERACTION. What is not: that the latent beats the
> incumbent anywhere.** A router exploiting the interaction must therefore beat the incumbent by
> arm-*blending* or by finding the specific targets — not by a difficulty threshold, which L3 also
> rules out independently.

**This is the FIFTH arrival of *an aggregate hides a structure that decides the disposition*** (after
the MDE, the operator forks, the bimodal ordering, and the step response).

## L3 — BUT TWO NATIVE-FREE DIFFICULTY PROXIES FAIL TO REPRODUCE IT, WITH THE WRONG SIGN

**The L2 split is ORACLE** — difficulty is the incumbent's own RMSD, which consumes the native. So
the decisive question is whether a native-free estimate reproduces it. **Two do not:**

    proxy                r with true difficulty    hardest-quartile effect
    ORACLE (true)              1.000               -0.455 [-0.717,-0.217]   LATENT WINS
    score_mean                +0.433               +0.385 [+0.149,+0.640]   latent LOSES
    pool_spread               +0.360               +0.405 [+0.053,+0.765]   latent LOSES

`score_mean`'s hardest quartile **overlaps the true-hardest on only 13 of 32 targets.**

> **The interaction is REAL and the deployable version is NOT DEMONSTRATED.** The router's task is
> therefore harder than "estimate difficulty": difficulty alone, at the accuracy a native-free proxy
> achieves (r ≈ 0.36–0.43), **identifies the wrong targets**. This does not close routing — a learned
> multi-feature router is a different and untested object — but it closes the obvious one-feature
> version, and it closes it *before* a lane spent the night on it.

---

## L4 — MY OWN FALSIFIER FIRED: PER-TARGET m IS REAL AND TRANSFERABLE, NOT NOISE-FITTING

`s22/mreal.py`, **126/126 complete**, 8 split-half repeats, half-pools of 250.

I pre-registered the hypothesis that **per-target m does NOT transfer** — that the 0.350 Å m-ladder
headroom is an artefact of taking an argmin over six flat, correlated arms. **It is refuted.**

    FIXED m=75, held out                     3.072
    SELECTED m (on half A) -> half B         2.833
    IN-SAMPLE oracle m (the suspect)         2.705

    PRIMARY  selected - fixed, HELD OUT   -0.239  SE 0.035  MDE 0.099  [-0.309,-0.173]  79W/47L
    -> TRANSFERS, at 2.4x its own MDE
    65% OF THE APPARENT HEADROOM ACTUALLY TRANSFERS
    median tied rungs at the half-ladder minimum: 1.00 of 6

**The per-target optimal m is a STABLE PROPERTY OF THE TARGET.** Selecting it on one independent
half of the pool and applying it to the disjoint other half beats the fixed incumbent m by
**0.239 Å**. The noise-fitting explanation is dead.

> ### THE DECISIVE QUALIFICATION, AND IT IS THE WHOLE FINDING
>
> **The selection on half A uses half A's RMSD TO NATIVE. It is ORACLE.** What is established is
> that the per-target optimum is **real and stable**, *not* that it is reachable native-free. This is
> a **ceiling with a transferability guarantee**, not an achievable arm.

### Why that matters more than a deployable arm would

Three independent native-free routers were run this campaign and **all captured ~0%**: the audit
lane's held-out length-quartile router (+0.031 [MDE 0.070], wrong-signed); Workstream C's Ridge over
16 arms (lost on **all five** arm-set configurations, one significantly at +0.070 [+0.005,+0.146]);
and C's Random Forest, which produced the textbook signature — **in-fold −0.276 Å, more than half the
ceiling, reversing to +0.079 Å under honest CV.** C's nested-CV gated router mostly learned to
disable itself.

**Those two facts together are far sharper than either alone:**

> **THE PER-TARGET SIGNAL EXISTS AND IS STABLE ACROSS INDEPENDENT DATA FROM THE SAME TARGET, AND NO
> NATIVE-FREE FEATURE SET THIS PROGRAMME HAS TRIED CAN SEE IT.**

**The diagnosis moves from "the ceiling is an artefact" to "the ceiling is real and our FEATURES are
blind to it."** That is a different research problem with a different next step, and it is the
opposite of what I expected when I wrote the file.

### What it does not license

The transfer is measured on **half-pools of 250**, so `avg_500` and `avg_250` are the same arm there
and the ladder's top rung is truncated — both sides face the truncation identically, so the
comparison is valid, but the 0.239 Å is a **half-pool** number and is not the full-pool headroom.
And the incumbent's own half-pool value (3.072) is **not** its full-pool value (3.048).

**Nothing is promoted.** An ORACLE selection that transfers is a statement about the target, not a
pipeline.

---

## L5 — THE GAIN IS GENUINELY PER-TARGET, NOT A GLOBAL m CORRECTION

The obvious alternative explanation for L4: **maybe half A merely says "m=150 beats m=75" globally**,
and the transfer is a global correction wearing per-target clothes. `s22/mreal2.py` tests it
directly by adding the two missing controls. **FINAL, n=126, complete:**

    half-A ladder means:  m=500 3.280 | 150 3.062 | 75 2.986 | 20 3.025 | 5 3.120 | 1 3.369
    -> ONE global best m chosen on half A is m=75 -- EXACTLY THE INCUMBENT'S OWN CHOICE

    FIXED m=75 (incumbent)                2.988
    RANDOM m per target                   3.153
    GLOBAL best m from half A (m=75)      2.991
    PER-TARGET selected m                 2.744

    per-target sel - FIXED m=75      -0.244  iid[-0.314,-0.176]  fold[-0.290,-0.191]   77W/49L
    per-target sel - RANDOM m        -0.412  iid[-0.482,-0.346]  fold[-0.473,-0.361]  119W/ 7L
    per-target sel - GLOBAL best m   -0.243  iid[-0.317,-0.174]  fold[-0.289,-0.187]   82W/44L
    GLOBAL best m  - FIXED m=75      -0.001  iid[-0.007,+0.004]  fold[-0.009,+0.006]   61W/65L

**The global correction is worth −0.001 Å — nothing, with a CI through zero.** The per-target
selection is worth **−0.243** against the *global* best m, essentially identical to its −0.244 against
the fixed incumbent.
**The entire effect is per-target.** The alternative explanation is dead, and m=75 is confirmed as
the correct *global* choice — the incumbent's m was never the problem.

Fold-clustered CIs are now printed beside the iid ones, per Workstream D's correct observation that
this project's standing practice requires both and `mreal.py` printed only one.

## L6 — WORKSTREAM D's INDEPENDENT REPRODUCTION AND ATTACK ON L4: IT SURVIVES

D reproduced from the artefact (`−0.2394 SE 0.0353 MDE 0.0990`) and attacked it four ways.

1. **Fold-clustered CI** `[−0.292,−0.180]` — essentially identical to iid. Not fragile to clustering.
2. **Jackknife** on the primary mean: leave-one-out range `[−0.244,−0.228]` against −0.239. **Not
   driven by a handful of targets.**
3. **The weakness I named myself** — the half-pool `m=500` rung being degenerate with "the whole
   250-member half" — tested by dropping that rung, same seeding so the identical 8×126 half-splits
   reproduce: **−0.1958 SE 0.0325, MDE 0.0912, 82% of the effect survives at 2.15× its own MDE.**
   **The truncation is not carrying it.**
4. **A statistic of mine correctly downgraded.** "Median tied rungs 1.00 of 6" uses an *exact
   floating-point* tie test on a **continuous** outcome, so it reads ~1 essentially by construction.
   It is **near-vacuous** and is **not** in tension with L23's within-1-SE tied set of 2.0/6 — it
   answers a different and much weaker question. **Do not read it as "sharper than L23".**

## L7 — FOUR INDEPENDENT ROUTER CONSTRUCTIONS, ALL FAILING, ONE SIGNIFICANTLY HARMFUL

| construction | lane | result |
|---|---|---|
| length-quartile, 5-fold held out | D | +0.031 [MDE 0.070], wrong-signed |
| Ridge over 16 arms, 5 configurations | C | **lost on all five**; `A_plus_latent` **+0.070 [+0.005,+0.146], CI EXCLUDING ZERO — significantly WORSE than doing nothing** |
| Random Forest, original features | C | in-fold **−0.276** (>half the ceiling) → **+0.079 held out, sign reversed** |
| Ridge + Random Forest, NEW candidate-set geometry features | C | +0.014 [−0.051,+0.079] (−3.8% capture, the flattest yet); RF again **−0.185 in-fold → +0.060 [+0.004,+0.117] held out, significantly worse** |

C built the geometry features on my suggestion — `spread(m)` = mean pairwise RMSD inside the top-m,
and its shape — because that is mechanistically the quantity m controls. **It is the least harmful
feature family tried and it still does not beat the incumbent.**

> ### THE CAMPAIGN'S CENTRAL RESULT
>
> **The per-target optimal m is REAL, STABLE across independent draws from the same target, GENUINELY
> PER-TARGET rather than a global correction, and worth 0.24 Å. And it is INVISIBLE to every
> native-free feature family tried: distogram confidence, score-distribution shape, and
> candidate-set geometry — linear and nonlinear, across four independent router constructions.**
>
> This is a **positive control paired with a negative result**, which is far stronger than either
> alone: the target signal is *proven present*, so the routers' failure cannot be dismissed as
> "there was nothing there." **The blindness is the finding.**

**D's caution, adopted:** a router-*construction* search carries its own best-of-K exposure separate
from the per-target min-of-K inside the ceiling. C's five configurations and any bin-count/feature
choices must be counted; a single pre-registered construction on a genuinely fresh fold is the clean
form. **"Not yet found" is the correct disposition, not "does not exist."**

---

## L8 — WORKSTREAM B: THE DIRECT PERTURBATION PROBE ANSWERS THE MANDATED PHYSICS QUESTION

30 targets, **3,780 trials**, complete with config-hashed artefacts. This is the experiment the
directive asked for at §43 — *which physical variable does each Hamiltonian actually respond to* —
and it answers it by direct controlled perturbation rather than by correlation.

**(1) ON COMPACTNESS, LEGACY AND AMBER ARE ANTI-CORRELATED AT MATCHED MAGNITUDE.**

    Spearman(dE, dRg):   Legacy  +0.454      AMBER  -0.444

**An independent, direct-perturbation replication of the mechanism Sprint 21 inferred from staged
preconditioning.** Legacy prefers compact structures; AMBER's steric singularity punishes exactly
that. The two energies are not merely differently scaled — **on the compactness axis they point in
opposite directions.** That is the fundamental reason their rankings disagree, and it is now measured
directly rather than inferred.

**(2) A HYPOTHESIS B GOT WRONG AND CAUGHT ITSELF.** The pre-registered concentration test came back
refuted *backwards*, because **torsion-index concentration is not Cartesian concentration** — NeRF
chain-building means one residue's torsion move cascades down the whole downstream chain. B
diagnosed it, registered a same-day addendum, and reran against **each potential's own top Hessian
eigenvector**. Corrected: **AMBER is ~29× [CI 4.5–543×] more disproportionately sensitive to its own
dominant curvature direction than Legacy is.** The mechanism is confirmed once properly
operationalised, and the original registration is left unedited.

**(3) STERIC RESPONSE IS SHARED, NOT SPECIALISED.** Both Hamiltonians respond to steric perturbation
similarly — so the difference between them is **not** that "AMBER sees sterics and Legacy does not."
It is the **compactness axis**, where they oppose.

**(4) NO DEFENSIBLE AMBER SURROGATE EXISTS.** The bonded-force-group subset tracks true AMBER
response at **ρ = 0.211** against a pre-registered bar of 0.7 — **refuted**, and for a principled
reason: **the steric singularity lives in the nonbonded and solvation terms a bonded subset cannot
see.** Directive §18's surrogate route is closed, with a mechanism.

**Self-caught defect:** a bootstrap RNG re-seeded inside the resampling closure produced a degenerate
CI on the first pass — found, fixed, and documented rather than silently corrected.

**Process risk B flagged, and it is a real one:** `s22/results/` shows fresh concurrently-written
artefacts from other lanes throughout B's run. **Multi-lane AMBER/OpenMM activity happened in
parallel without a central serialisation gate.** B's own runs completed cleanly and its artefacts
verify, but the brief's serialisation rule was not enforceable by the lanes themselves — that is a
coordinator failure, not B's.

---

## L9 — HEDGING ACROSS m IS REFUTED. MY FALSIFIER FIRED.

If the per-target optimum cannot be **predicted**, it can in principle be **hedged**: averaging the
structures emitted at several widths is a **fixed, global, native-free** operator — same thing done
to every target, no feature consulted, no fold fit, nothing to overfit. `s22/mhedge.py`, n=126:

    arm            RMSD    virt bond   vs fixed m=75
    fixed_75     3.0483      2.961     (incumbent)
    hedge_all    3.0482      2.778     -0.0002  MDE 0.056  iid[-0.038,+0.039] fold[-0.035,+0.041]
    hedge_core   3.0367      2.930     -0.0117  MDE 0.027  iid[-0.030,+0.007] fold[-0.021,-0.002]
    hedge_wide   3.0922      2.628     +0.0438  MDE 0.099
    hedge_narrow 3.0713      3.054     +0.0230  MDE 0.068
    oracle_m     2.7268        -       <- ORACLE ceiling

**The pre-registered falsifier — "if `hedge_all` AND `hedge_core` both fail to beat `fixed_75` past
their own MDE with a CI excluding zero" — FIRES.** `hedge_all` is exactly zero. `hedge_core` is
**−0.0117 against its own MDE of 0.027**: below its detection threshold, and its iid CI includes
zero even though the fold-clustered CI does not. **Not measured, not promotable, 0.4% of baseline.**

**And the contraction column is the reason to be glad it failed.** Averaging contracts the backbone
(physical virtual Cα–Cα is 3.805 Å; the incumbent already sits at 2.961). **Averaging averages
contracts further** — `hedge_all` reaches **2.778** and `hedge_wide` **2.628**. The one arm that
*improved* the geometry, `hedge_narrow` at 3.054, is also the worst on RMSD. **Any Å the wide hedges
bought would have been bought by making the structure less physical.**

**So the 0.243 Å per-target signal is unreachable by SELECTION (L7, four constructions) and by
COMBINATION (here).** Both doors tried, both closed.

## L10 — THE BOUND: SAMPLE SIZE, NOT FEATURES, IS THE DEFENSIBLE EXPLANATION

Workstream D derived a finite-class/Bernstein uniform-convergence bound using **this project's own
measured** σ ≈ 0.41 Å per-target achieved-vs-fixed SD (cross-checked independently from `mreal.json`
and `routercv.json`), at n=100 — one training fold:

    router class                                      gap at n=100   n for 0.24   n for 0.48
    1 global threshold, K=13 arms                         0.39          180           78
    3-bin tercile, 1 feature (routercv's own design)      0.96          513          222
    linear router, p=5 features                           0.57          281          121
    linear router, p=15 features (routercv's full set)    1.18          650          281
    shallow tree depth=3, 15 features                     1.82         1045          451

**Even the simplest possible router — one global threshold — has a generalisation gap comparable to
the ENTIRE 0.482 Å ceiling at the actual training-fold size, and needs n≈180 to resolve half of it,
more than the whole 126-target benchmark.** Everything Workstream C tried sits in rows needing
**280–1045 targets, 2–8× the instrument.** Dropping the range term and using variance alone still
gives 0.31 at the tercile row — the same order, so it is not an artefact of that choice.

**And `mreal2`'s own result is exactly what the bound predicts**: the one router class simple enough
to possibly work at this n — a single global m — is worth **−0.001 Å, nothing.**

> **L7 IS RESTATED. It is not "our features are blind." It is: AT n=126, NO ROUTER CLASS RICH ENOUGH
> TO EXPRESS THE SIGNAL IS LEARNABLE, AND NO ROUTER CLASS SIMPLE ENOUGH TO BE LEARNABLE CAN EXPRESS
> IT.** Four empirical failures become one structural explanation.

**D's two caveats, kept:** this is a **worst-case, distribution-free** bound — *"no proof of
learnability exists at this n"*, **not** *"a router is provably impossible"*; real algorithms
sometimes beat such bounds. And **targets are not i.i.d.** (fold/family structure), so the bound is
**optimistic** relative to reality, not pessimistic.

**The forward implication is concrete:** more cluster-disjoint targets, or an explicit
variance-reduction estimator (shrinkage/hierarchical, trading bias for variance deliberately) — **not
another pass through the same router class at the same n.**

---

## L11 — EXTERNAL CORROBORATION FROM A THIRD INDEPENDENT GROUP, AND THE TENSION WE MUST NOT HIDE

**arXiv:2606.21241**, *"Assessing Cost Hamiltonian Reliability in Quantum Protein Structure
Prediction"* — recovered via the HTML mirror after the PDF route failed. Verbatim:

> *"the average RMSD across all solutions is, on average across all peptides, better than the optimal
> solution according to the cost Hamiltonian"*

**For peptides ≤15 residues — almost exactly this instrument's range.** Their Spearman(cost, error)
is **negative for small peptides**, null at length 50, and positive only for larger proteins.

**That is this project's own shape — the physics energy anti-ranks, and optimising it hurts —
obtained independently on a DIFFERENT energy (Miyazawa–Jennings contact, not AMBER or Legacy) and a
DIFFERENT representation (tetrahedral lattice, not continuous torsion).** With arXiv:2609.02113 that
is **three independent groups**.

> **UPGRADE ADOPTED: from "our pipeline has a problem" to "this is a property of the method class."**
> A negative result reproduced across three groups, two energies and two representations is a
> different and much stronger object than a local finding.

### The tension, stated rather than quietly dropped

**Their paper attributes the failure partly to RESOLUTION** — correlation improves with more
interaction shells, larger and more detailed instances. **Our data cuts against that specific
mechanism.** AMBER is far richer than Legacy and **does not rank better**: it is worse in most
configurations (C8/C10 in Sprint 21's register), and **both lose to a matched random tail**. If
"more physical detail fixes it" were the operative variable, **AMBER should be our best-ranking
energy, and it is not.**

**Two readings, and I do not know which is right:** either "interaction shells" in a lattice contact
model is simply a different axis from "force-field detail" in a continuous all-atom model (compatible
claims about different variables), or there is a **genuine disagreement about mechanism**, in which
case our matrix is **evidence against the resolution-fixes-it story**. **The corroboration on the
headline must not quietly import their proposed mechanism** — D flagged exactly this and it is
adopted.

### Novelty, checked rather than assumed

- **CVaR-VQE for peptide folding is a KNOWN combination** — QuPepFold (PLoS ONE, Feb 2026;
  tetrahedral lattice, Miyazawa–Jennings, energy-only, no RMSD) plus two 2024 CVaR-VQE-vs-MD peptide
  papers. **Nothing in this campaign may claim it as novel.** *(Caught before it could be
  misreported.)*
- **Plausibly still open:** continuous-torsion + CVaR-VQE + all-atom AMBER *jointly*; the
  **readout-H / training-H separation** as a measured architecture lever; and the candidate-identity
  register results below.
- **The Legacy/AMBER compactness anti-correlation** (B's direct perturbation probe): **no specific
  match found.** The broader knowledge-based-vs-physics-based potential debate is ~20 years old
  (Ben-Naim on the physical interpretability of potentials of mean force), but the **specific claim**
  — matched-magnitude anti-correlation on a compactness axis, by **controlled perturbation** rather
  than observational correlation, with a **shared** steric channel — was not found. Classified
  **plausible empirical/mechanistic novelty**, with D's caveat that one search pass is *evidence of
  absence, not proof*, and a CASP-era decoy-discrimination search is the natural next check.

---

## L12 — THE EXTENDED ABLATION (DIRECTIVE §57) AT n=126: DISTANCE ALONE WINS, AND THE TORSION PRIOR IS THE WORST TERM

`s22/ablate.py`, **126/126 complete**. Pool instrument, shipped top-75 tail readout, terms
rank-normalised per target (declared before any RMSD was read). **This is a DIFFERENT INSTRUMENT
from Sprint 21's n=12 CVaR-VQE matrix and the two must never be quoted as one table** — it exists
because the directive asks *which terms lower RMSD*, and the pool answers that at 10× the targets.

| cell | avg (top-75) | argmin | vs D alone, fold-clustered | MDE |
|---|---|---|---|---|
| **D** | **3.0483** | 3.4540 | — | — |
| D+G | 3.1107 | 3.4971 | +0.0623 [+0.000,+0.122] | 0.110 |
| D+L | 3.2129 | 3.5253 | +0.1646 [+0.052,+0.242] | 0.160 |
| D+T | 3.3435 | 3.4956 | +0.2951 [+0.150,+0.448] | 0.198 |
| D+T+G | 3.4158 | 3.5043 | +0.3675 [+0.229,+0.505] | 0.221 |
| G | 3.4240 | 4.0912 | +0.3756 [+0.212,+0.533] | 0.260 |
| **RANDOM tail** | **3.4259** | — | **+0.3776 [+0.219,+0.511]** | 0.242 |
| D+L+T+G | 3.6533 | 3.7983 | +0.6050 [+0.418,+0.746] | 0.294 |
| D+L+T | 3.6772 | 3.6832 | +0.6288 [+0.446,+0.764] | 0.301 |
| L | 3.7543 | 4.4897 | +0.7059 [+0.493,+0.840] | 0.331 |
| L+T | 3.8965 | 4.0499 | +0.8481 [+0.633,+1.012] | 0.371 |
| **T** | **3.9470** | 4.0268 | **+0.8986 [+0.677,+1.085]** | 0.391 |

**The falsifier did not fire: no combination beats Distance alone.** `D+G` is the closest and sits
**below its own MDE** with a CI touching zero — a tie, not a win, and exactly as predicted, since the
pool's members are real windows and therefore already valid chains.

> **THE TORSION PRIOR IS THE SINGLE MOST HARMFUL TERM TESTED — +0.899 Å alone, and +0.295 even when
> added to Distance.** The directive named `H_T = −Σ log p(φ,ψ|S)` as a candidate (§25); measured on
> this instrument it is **worse than a matched random tail by 0.52 Å.**

**Six cells are at or worse than a matched RANDOM tail:** G, D+L+T+G, D+L+T, L, L+T, T. **This
extends Sprint 21's C10 — "both physics energies are worse than chance" — to the torsional and
geometric terms**, at full n with every CI excluding zero. **Adding physically-motivated terms to
the structural objective makes it worse, monotonically in the number of terms added.**

**Scope, stated:** `H_AMBER` was **not** re-run here — its deployed form is `E ∘ Relax_50`, costing
an OpenMM context per candidate per target, and its cells already exist at n=12 (VQE) and n=126
(tail). **Stated, not silently skipped.**

---

## L13 — THE QUANTUM LANE CLOSES: THE ENTROPY-REGULARISED TAIL DOES NOT BEAT THE CLASSICAL SELECTOR

*Workstream A produced these artefacts and then stalled in a wait loop on an already-written file —
the exact "absence of signal read as a healthy state" hazard another lane named tonight. The
analysis below is the coordinator's, on A's data, and is labelled as such.*

`a2b_entropy.json`, n=16 targets, 200 iterations, α=0.15. Workstream A itself added the correct bar:
**not a size-matched random draw, but the size-matched CLASSICAL TOP-m selector** — because Sprint 21
priced the averaging operator at −0.28 to −0.31 Å **independent of the energy**, so *"a bigger tail
averages better"* is the null, not the finding.

    classical top-m tail-average (THE BAR)  3.2998
    classical score argmin                  3.8923
    pool ORACLE                             1.7802

    T        VQE tail   rand ctrl    ESS    tail   vs THE BAR                       vs random
    0          3.7232      4.6989   1.02     2.1   +0.4234 MDE 0.545 [+0.040,+0.769]  -0.976
    0.05       3.5559      4.3021   5.50    10.9   +0.2561 MDE 0.496 [-0.100,+0.578]  -0.746
    0.2        3.3125      3.9219  33.70    44.5   +0.0127 MDE 0.220 [-0.149,+0.149]  -0.609
    0.5        3.2711      3.8522  62.08    68.5   -0.0287 MDE 0.105 [-0.110,+0.021]  -0.581

**The entropy regulariser works exactly as designed** — ESS climbs 1.02 → 62.08 and the tail widens
2.1 → 68.5 as T rises, so the degeneracy lever is real and controllable. **It beats a size-matched
random draw at every temperature (−0.58 to −0.98). It never beats the classical top-m selector.**
The best cell (T=0.5) is **−0.029 against its own MDE of 0.105** — a null — and it gets there by
converging on classical top-m behaviour, which is what a flattening distribution must do.

> **Item 2 closes NEGATIVELY and CLEANLY. It extends Sprint 21's nine-arm list of native-free
> selectors that fail to beat the incumbent rather than escaping it.** And it is a textbook
> illustration of why the bar matters: **against the wrong control this would have read as a
> 0.6–1.0 Å quantum win at every temperature.**

## L14 — THE READOUT-H EFFECT REPLICATES ON A COMPLETELY DIFFERENT SUBSTRATE

`a_readouth.json`, n=16, α=0.15. Swapping **which Hamiltonian is read out** on the
candidate-identity register:

    R_score_argmin_dist   3.8923      R_score_argmin_leg   5.3998
    R_score_tailavg_dist  3.2998      R_score_tailavg_leg  4.4770

    readout swap (distance vs Legacy), tail-average  -1.177  SE 0.420  MDE 1.177  [-2.007,-0.410]  11W/5L
    readout swap (distance vs Legacy), argmin        -1.508  SE 0.668  MDE 1.873  [-2.820,-0.312]  11W/5L

**Sprint 21 measured the readout-H lever at −0.697 Å (Legacy) and −1.299 Å (AMBER) on the torsion
basin latent. It reproduces here at −1.18/−1.51 on the candidate-identity register — a different
encoding, a different state space, a different experiment.** Both CIs exclude zero.

**Both magnitudes sit AT or BELOW their own MDEs at n=16**, so **directions ESTABLISHED, magnitudes
NOT MEASURED** — the same disposition Sprint 21's matrix carried, for the same reason.

**This is the campaign's only architectural lever that has now survived two independent substrates.**

---

## L15 — MULTI-STAGE VQE: STAGING HELPS THE TAIL A LITTLE, NEITHER ARM BEATS THE CLASSICAL BAR, AND THE A-PRIORI CLOSURE HOLDS THROUGHOUT

`a3_multistage.json`, n=16, stage 1 = `H_Legacy` (60 iters) → stage 2 = `H_Distance` (140 iters) —
the directive's §36/§37 architecture, where the scalar being minimised **changes between stages**.

    classical top-m tail-average (THE BAR)   3.2998        classical score argmin  3.8923

    arm       tail-avg   argmin   tail size   ESS
    single      3.8188   3.8923         1.2   1.00
    staged      3.6461   3.8923         3.9   1.09

    staged - single (does staging help?)   -0.1727 SE 0.0795 MDE 0.223 [-0.334,-0.033]  8W/ 8L
    staged - classical top-m BAR           +0.3463 SE 0.1726 MDE 0.484 [+0.009,+0.666]  4W/12L
    single - classical top-m BAR           +0.5190 SE 0.2056 MDE 0.576 [+0.099,+0.886]  3W/13L

**Staging helps the tail readout by −0.173**, with a CI excluding zero — but **below its own MDE of
0.223 and at 8W/8L**, so **NOT MEASURED**. It widens the tail (1.2 → 3.9) and lifts ESS (1.00 →
1.09), so the mechanism does something. **Neither arm beats the classical top-m bar.**

> **AND THE A-PRIORI CLOSURE HOLDS UNDER STAGING.** The argmin readout is tied to the classical
> argmin on **100% of cells for BOTH the single-stage and the staged arm.** Changing the Hamiltonian
> between stages cannot move an argmin over an exactly-enumerable register — exactly as D1's theorem
> requires. **A third confirmation, from an experiment designed to test something else.**

**Directive §36/§37 is answered: a staged Hamiltonian progression is measurable, does move the tail,
and does not reach the classical selector.**

---

## L16 — THE MECHANISM PROOF: A2b's NULL AND THE READOUT-H "LEVER" ARE THE SAME FACT

Workstream A's final report traces both results to a single property of the implementation, verified
by **set equality** rather than by RMSD closeness. This supersedes the empirical framing in L13/L14
and in the published report's §04.

> **`core.quantum.cvar_from_probs` orders states by ENERGY, not by trained probability, and assigns
> mass to a PREFIX of that order.** Therefore the **tail-support SET is provably always the classical
> top-m set** for whatever m the training realises — **regardless of ansatz, entangler, or
> temperature.**

**Verified directly, not inferred:** trained-tail candidate indices are **element-for-element
identical** to `argsort(scores)[:m]` in every spot-check (18 cells, 7 targets); and on the argmin,
**max |diff| = 0.0 across all 64 cells** of the readout-H experiment, *independent of which
Hamiltonian trained it*.

### The consequence, and it collapses two findings into one

**L13's null and L14's "architectural lever" are not two results. They are one.** The
−1.177 / −1.508 readout-swap effect is simply **"Distance ranks better than Legacy," made visible
through a readout that training can never move away from either answer.**

- It **still replicates Sprint 21's sign and rough magnitude on a genuinely different substrate**,
  which is worth keeping.
- But **it has an a-priori mechanism, not an empirical one** — and it is **not a quantum lever**. It
  is the classical ranking quality of the two Hamiltonians, surfaced through a readout the VQE cannot
  influence.
- **A3 is a third instance of the same fact**: realised tail size 1.25 (single) vs 3.89 (staged), so
  the Legacy warm-start leaves stage 2 less collapsed and lands a marginally better m. The −0.173
  [−0.335,−0.036] is that, not a staging benefit.

**So the campaign's "only lever surviving two substrates" survives — but as a statement about the two
Hamiltonians, not about the quantum selector.** That is a weaker and more accurate claim, and it must
replace the stronger one everywhere.

### The one door the theorem does not close

**A probability-WEIGHTED tail average** — never tried here, and not what this project's shipped
operator or this workstream use — **is not covered by the set-equality argument.** Every result above
concerns an *unweighted* average over a tail whose *membership* is classical. **That is the one place
a genuinely different quantum-shaped signal could still live on this substrate**, and A named it
rather than leaving the closure looking total.

### Reconciliation of the two lanes' numbers — no discrepancy

A confirms L13's table exactly (T=0.5: cnot 3.264 + none 3.278 → 3.271 against my 3.2711; ESS
52.63/71.54 → 62.09 against my 62.08). My pooled readout figure is the per-target average of A's two
separately-registered comparisons: **train=Legacy fixed −1.929, independently significant, 14W/2L**;
**train=Distance fixed +0.523, same sign, not independently significant, 6W/10L.** That explains the
11W/5L sitting between A's two tallies. **Same conclusion, two compatible slices.**
