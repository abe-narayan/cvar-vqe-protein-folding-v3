# Sprint 21 — CLAIMS REGISTER

Every row carries its **n**, its **effect with a CI**, its **per-comparison MDE** (`2.8016 × SE`, an
operator and not the programme's old 0.084 Å constant — L8), and its **disposition**. ORACLE
quantities are labelled ORACLE and are **ceilings, never achievable arms**. Point-cloud and
built-chain RMSDs are never compared across bases.

**Instrument**: 126 cluster-disjoint targets, 5 folds, full-chain Cα-RMSD after Kabsch. The 60-target
benchmark remained **sealed** throughout — not inspected, probed, tuned against, or used for any
architecture decision.

**Incumbent**: 3.048 Å (point cloud) / 3.204 Å (built chain), score-filter → top-75 coordinate
average over the K=500 BLOSUM pool.

---

## A. THE SEARCH HALF — CLOSED, THREE INDEPENDENT WAYS

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| A1 | The VQE latent is exactly `2**n` (one qubit/residue, two von Mises basins): **median 4096**, and the 8192-evaluation budget **equals or exceeds the entire space on 75 of 126 targets**. The deployed sampler was **resampling a space it could have enumerated**. | 126 | — | **ESTABLISHED** (L11) |
| A2 | The **EXACT argmin over the entire enumerated latent** ties the argmin of a zero-evaluation retrieval pool. | 126 | +0.085 [-0.018,+0.169], MDE 0.197, 53W/73L | falsifier **DID NOT FIRE**; magnitude **NOT MEASURED** (L14/L17) |
| A3 | Given **exhaustive** access to every one of `2**n` configurations, the objective still misses the best one by 1.7–1.8 Å. **Zero reversals in 126 targets.** | 126 | +1.815 [+1.613,+2.059] squared, 0W/**122L**; +1.720 [+1.499,+1.975] Bayes, 0W/**119L** | **ESTABLISHED** (L17) |
| A4 | The objective is **not noise** — it beats the mean of the very space it searches. | 126 | -1.367 [-1.533,-1.210], 109W/17L | **ESTABLISHED** (L17) |
| A5 | Budget-versus-latent-size was **not** the operative variable: on the 51 targets where the budget was *smaller* than the latent — the only ones genuinely searching — exhaustive search still does not beat the pool, if anything less. | 51 | +0.174 [-0.008,+0.372], 19W/32L | **ESTABLISHED** (L17) |
| A6 | The objective's argmin over **512 random draws** is indistinguishable from its argmin over the **entire 8192-configuration latent**. Exhaustive enumeration buys nothing over 512 draws for the deployed readout. | 45 *(interim)* | -0.045 [-0.177,+0.137], MDE 0.153 | **INTERIM** — Workstream D, n=75 pending |
| A7 | The **mode restriction** used to enumerate deterministically costs nothing; the stochastic von Mises decoder behaves the same. | 75 | +0.083 [-0.326,+0.156], MDE 0.255 | **null**, the control that validates A2–A4 (L14) |

> **A2–A4 together convert `search-saturates-discrimination-binds` from an inference over six
> instruments into a property of the instrument.** No ansatz, optimiser, budget, or quantum resource
> can recover the 1.72 Å, because all of them are ways of arriving at an argmin that has now been
> computed exactly.

---

## B. THE QUANTUM PILLAR — GENUINE, CORRECT, AND MEASURED AT NULL ON RMSD

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| B1 | Measured bond dimension **reproduces the analytic `2^layers` exactly** on all seven MPS variants; a ring circuit reaches **χ = 48**. | 7 variants | exact | **ESTABLISHED** |
| B2 | The **entire ansatz ladder χ=1 → 48 spans 0.37 Å**, and gradient norms *rise* with depth — so it is **not** a trainability null. | — | 0.37 Å total span | **CLOSED** |
| B3 | Adam beats every alternative **on the CVaR loss** at 0W/20L on both Hamiltonians, 5/5 folds — yet **no optimiser beats `best_of_N` on RMSD**. | 20 | — | **the two columns disagree, and that is the answer** |
| B4 | **Training the circuit at all** is worth nothing. | 20 | +0.003 [-0.136,+0.144] | **null** |
| B5 | QNG fails because the **MPS Fisher is rank-deficient by construction** — a property of the *parameterisation*, not the landscape — and is significantly **worse** than no optimisation. | 20 | +0.097 [+0.006,+0.185] | **ESTABLISHED**, mechanism identified |
| B6 | **Changing the READOUT Hamiltonian is the one non-null lever in the selector**, at zero extra budget on the identical evaluated set. | 20 | **-0.697 [-1.059,-0.352]**, 5/5 folds; replicated in two blocks to **0.011 Å**; holds on `best_of_N` where no training occurs | **ESTABLISHED** — prices a design choice, does **not** beat the incumbent |
| B7 | Changing the **TRAINING** Hamiltonian is worth nothing. | 20 | +0.056 [-0.020,+0.129] | **null** |
| B8 | QNG on **AMBER** — the badly-conditioned landscape where its case was strongest a priori. | 3/20 targets | — | **NOT MEASURED, OPEN** (compute-bound) |

> **B6 vs B7 is the sprint's cleanest architectural statement: the Hamiltonian inside the training
> loop is nearly free to change; the one in the readout is not.** It also means the mandated
> `H(λ) = (1−λ)H_L + λH_A` continuation operates in the dimension measured as null (B7), so any
> effect along λ must be claimed as an **optimisation-path** effect, with B7 quoted beside it.

---

## C. LEGACY vs AMBER, AND THE READOUT

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| C1 | Picking the **lowest-energy** structure is worse than picking at random, **for both physics energies**. | 20, subset(20) | AMBER 4.990 \| Legacy 5.487 \| distogram 3.676, pool mean 4.739 | **ESTABLISHED** (Sprint 20; provenance corrected — see M4) |
| C2 | The two distance functionals (**squared** vs **shipped Bayes risk**) are interchangeable at the argmin on pool-like structures. | 126 | +0.050 [-0.019,+0.156], MDE 0.124; median ρ **0.973** | **null** — and D's own opposite hypothesis **REFUTED** |
| C3 | That interchangeability is a property of the **pool's** structure distribution and does **not** transfer unchanged to the latent's. | 75 | +0.114 SE 0.069 on the latent, vs +0.050 on pool windows | **ESTABLISHED** — a measured null does not automatically transfer |
| C4 | Legacy's partial rank skill beyond the distogram. | 42 | **-0.0076** | **≈ zero** |
| C5 | The **mandatory 7-cell matrix** (Legacy \| AMBER \| Legacy+AMBER \| Distance \| Distance+Legacy \| Distance+AMBER \| Distance+Legacy+AMBER), separably evaluable. | 4/12 targets | — | **INCOMPLETE — Workstream A, ~90 min outstanding** |

---

## D. THE ENCODING LEVER — NOT SUPPORTED

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| D1 | At **exactly matched step counts** the θ vs (sin θ, cos θ) encoding effect is **zero**, and on DIST/spsa that null is **properly powered** (own MDE 0.052, below the 0.084 bar). | 60 | -0.015 [-0.050,+0.022] | **NOT SUPPORTED** |
| D2 | The Sprint-20 result was a **step-count confound**: at fixed budget the embedded parameterisation has twice the dimension, so adam takes 4.6 steps against 9.6. | 60 | AMBc/spsa shrinks from -0.649 [-1.205,-0.165] "sig" to **-0.295 [-0.636,+0.022]** | **ESTABLISHED** |
| D3 | **The symmetry argument.** At matched budget the effect **flips sign with the objective** — LEG -0.165 [-0.273,-0.063], DIST +0.100 [+0.012,+0.189], both CIs excluding zero. *A change of coordinates that provably does not change the physics cannot help on one energy and hurt on another.* | 60 | — | **ESTABLISHED** |
| D4 | **And the step-count response PREDICTS the effect's size.** DIST residual **+0.004 against a predicted +0.087**. | 60 | LEG residual -0.039 [-0.122,+0.025]; DIST +0.004 [-0.023,+0.027] | **ESTABLISHED** — sufficiency, not merely impossibility |
| D5 | "Optimising harder hurts" is a property of the **LEGACY objective**, not of the instrument. The deployed distogram objective goes the **other way**. | 60 | LEG **+0.091 [+0.078,+0.109]**; DIST **-0.087 [-0.168,-0.018]** | **ESTABLISHED** — corrects how `search-saturates-discrimination-binds` was being quoted, including by the coordinator |
| D6 | A real **gauge defect** in the embedded arm: `‖u‖` grows 1.000 → 2.5 (max 4.03) along a direction whose true gradient is **exactly zero** (arctan2 is scale-invariant), so every unit of radial motion is injected finite-difference noise; `adam_fd` gives that pure-noise direction a full-size step. **`_EmbField` is running a different OPTIMISER, not a different set of COORDINATES.** | — | — | **BUG, unfixed by choice** — the durable output of the lane |
| D7 | Whether **any** encoding effect exists. The AMBc panel's own MDE is ~0.53 Å and cannot exclude an effect of the originally reported size. | 10 | — | **NOT MEASURED**, with a 55%-smaller point estimate — *not* "no effect exists" |

---

## E. DISCRIMINATION — WHERE EVERY REMAINING ÅNGSTRÖM IS

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| E1 | The **ORACLE-best latent configuration sits at median percentile 5.62** of the objective's own ordering (mean 20.27, max 99.91) against a **50.0** no-skill null. | 126 | <1% on 21/75, <10% on 42/75 | **ESTABLISHED** |
| E2 | **The objective's ordering skill is -1.184 Å at M=1 and gone by M=64.** The whole fall of the top-M ladder is **min-of-M**, not ordering. | 75 | M=1 **-1.184 [-1.448,-0.918]** 63W/12L; M=8 -0.202; M=64 +0.039; M=512 **+0.107 [-0.051,+0.300]** | **ESTABLISHED** — the D9 matched control; **overturned the coordinator's first reading** |
| E3 | Unconditionally, the ordered window's excess over the global best exceeds a random window's — recovering E2 from a construction with **no split and no selection effect**. | 126 | ordered +0.507 [+0.359,+0.673]; random +0.388 [+0.317,+0.460]; difference +0.119 [-0.038,+0.290] | **ESTABLISHED** (two routes to one number) |
| E4 | **The aggregate null is bimodal.** The objective contains the answer in its top-512 on **59%** of targets (random gives 6%) and is a **concentrated wrong region** on the other 41% — **3.26× further from the ORACLE than an arbitrary window.** | 126 | contained 74/126; failure cell excess 1.228 Å vs random 0.377 Å | **MECHANISM** — the failure cell's **sign is FORCED** by conditioning; only the **3.26× magnitude** is free |
| E5 | An **averaging** readout extracts 0.45 Å from the same ordering an **argmin** readout finds worthless — because `d_out = 1.16·d_set_mean + 0.04·d_set_best`. | 126 | -0.449 [-0.520,-0.348], 89W/31L, 5/5 folds | **ESTABLISHED** — **no statement about an objective's usefulness is well-posed until the readout is fixed** |
| E6 | Whether a **native-free** signal can separate E4's two regimes (worth ~0.4 Å if yes). | 45 *(interim)* | best-of-K permutation null: best-of-10 \|AUC\| = **0.625 is the MEDIAN of pure noise** at n=75 | **PENDING** — anything in 0.55–0.65 will be reported **NOT MEASURED / NOT DEMONSTRATED** |

---

## F. THE GENERATIVE SOURCE — CLOSED BY ARITHMETIC

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| F1 | Swapping the retrieval pool for the generative latent, through the **identical shipped operator** and the **same basis**, **LOSES**. | 126 | **+0.390 [+0.224,+0.558]**, MDE 0.242, 47W/79L, **5/5 folds** | pre-registered falsifier fired **against** the latent |
| F2 | The latent is also **poorer**: at matched set size its ORACLE is worse than the pool's. | 126 | +0.558 [+0.413,+0.701], 35W/91L | **ESTABLISHED** |
| F3 | **Under the most generous possible reading** — both sources, perfect selector, ORACLE of the union — the latent adds **0.103 Å**, supplying the better structure on 35/126. | 126 | -0.103 [-0.146,-0.067] | **CLOSED by arithmetic**, not by extraction failure |
| F4 | The objective *does* do work in a 75-member selection, contradicting a naive reading of E2. | 126 | -0.429 [-0.561,-0.298], 91W/35L | reconciled by E5 (min vs mean) |
| F5 | **Basis price** on the averaging operator, so the deployment footnote is decomposable rather than trusted. | 126 | +0.016 [-0.004,+0.037] — the composite is **96.5% source, 3.5% basis** | **ESTABLISHED** |
| F6 | Sensitivity: dropping the 9 targets where "512 draws" is the *entire* latent makes the latent look **worse**, monotonically. | 126→117→93 | +0.390 → +0.427 → +0.474 | **ESTABLISHED** |
| F7 | *Generative sets average worse than retrieved sets, as a property of the generator* (error coherence). | 45 | latent top-M averaging gain **-0.478** at M=64 vs the pool's -0.377 | **REFUTED** — withdrawn by its proposer; **must not appear as a next-sprint direction** |

---

## G. THE ACHIEVED COLUMN — THE 1.338 Å GAP IS ENTIRELY UNCAPTURED

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| G1 | A **perfect in-pool selector** is worth 1.338 Å and would land the instrument at **1.711 Å — past the mission target — using the pool that already ships.** | 126 | **-1.338 [-1.558,-1.138]**, 123/126 | **ESTABLISHED** (ORACLE) |
| G2 | **Not one of nine native-free arms beats the incumbent**, and the two nearest are flat against it. **The incumbent is already at the top of everything this programme knows how to do.** | 126 | top-150 +0.023; top-20 +0.042; top-5 +0.183; argmin +0.406; whole-pool medoid +0.658; medoid-of-top-75 +0.234 | registered falsifier ("more than half") **did not fire and could not have** |
| G2a | **SCOPE, and it is narrower than G2 invites.** What was measured is that **nine specific selectors** capture none of the gap — **eight of them variations on consensus/typicality**, a family project memory already prices at -0.172 Å. The claim is therefore **"the CONSENSUS FAMILY is exhausted"**, *not* "in-pool selection is closed". The 1.338 Å remains an **unclaimed** ceiling, not a **proven-unreachable** one. | 126 | — | **binding scope limit** (Workstream D) |
| G3 | The incumbent's value is a **conjunction**: the argmin alone (-30% of the gap) discards averaging; the whole-pool medoid alone (-49%) discards the score. **Each half loses more than either is worth.** | 126 | — | **ESTABLISHED** |
| G4 | **m\* has not moved from S8-11's 75** — a registered forward prediction (`operator-consumes-set-mean`: m\* shrinks as the objective improves), resolved **NEGATIVE**. | 126 | m=150 and m=20 both **FLAT** against 75 | **the objective has not improved** — as four other instruments say independently |
| G5 | m\* is **INTERIOR**, and the ladder **cannot distinguish m ∈ [20,150]**. | 126 | median tied-set size 2.0 of 6 rungs | **ESTABLISHED** — the non-vacuous form of G4, per D's caveat that an interior minimum on a nested ladder is near-vacuous |

---

## M. METHODOLOGY — THE ERRORS, AND THE RULES THEY PRODUCED

*The single most reproduced fact of this sprint: **an aggregate hides a structure that decides the
disposition.** It appears as M1 (the MDE), M2 (the operator forks), E2/E4 (the bimodal ordering), and
D5 (the step response) — four independent arrivals.*

| # | finding | consequence |
|---|---|---|
| M1 | **MDE is per-comparison, not per-instrument.** Measured across 26 comparisons it ranges **0.11× to 9.54×** the constant quoted for four sprints. The coordinator's own brief called a **6.09-SE, power-1.000** effect "a quarter of the MDE." | **BRIEF §1 now carries the operator `2.8016 × SE`.** |
| M2 | **Unstated operator differences compound AND align with the author's expectation.** The coordinator's primary carried **three** — argmin vs average, squared vs Bayes, built chain vs window — worth **+0.215 Å against a true effect of +0.025**. **89% was operator, and all three pointed the way the author expected.** | **BRIEF §7 RULE 0**: enumerate the five operator forks (functional, basis, readout, normalisation, null) in the docstring and **name the alternative not taken**. |
| M3 | **Rule 0 does not say who writes the forks.** Two of the three were caught by a lane that *did not know which way the coordinator wanted the answer*, so every fork looked symmetric. **That is a property of the ROLE.** | **Rule 0, clause 2**: where one person both designs a comparison and has a stake in its direction, the forks are enumerated by someone who does not. *(First applied the same hour: the fork list for `latentsel.py` was sent before its numbers existed, and returned two corrections.)* |
| M4 | **A best-of-K result's null is the distribution of the MAXIMUM.** Simulated at this panel's size: a single predictor's folded \|AUC\| has a 95th percentile of **0.633** under pure noise, and the **best of K=10 has a MEDIAN of 0.625**. A "nothing below 0.65 is usable" threshold would have licensed a value noise produces half the time. | **Rule 0, clause 3.** A **pre-specified ordering** (a monotone prediction across arms) is immune in a way a pre-specified *set* of arms is not — reporting all K does not help, because the reader's eye goes to the winner. Only a claim that can **fail across the arms** is immune. |
| M5 | **Atomicity protects readers from a half-written file. It does nothing against a second writer.** Two *atomic* writers to one path produce a **perfectly well-formed chimera** of two configurations. | Fix is a **lock or a per-run output path**, not stronger writes — derive the filename from the config. Caught only because a **completion flag demanded the full key set**, not just the row count. |
| M6 | **A completion flag can pass vacuously.** `tailprice.py`'s counted *skipped* rows, so 3/3 failures wrote `complete: true`. | Flags require **every row to carry every declared key**. |
| M7 | **Provenance.** The brief's opening row (`AMBER 4.990 \| Legacy 5.487 \| distogram 3.676`) was quoted without its operator: **n=20, subset(20), rebuild basis, squared functional** — three differences at once. subset(20) is simply harder (pool mean 4.732 vs 4.453). | **Never benchmark an n=126 arm against 3.676 or 4.739.** Sprint 20's own conclusion is untouched — its four numbers share one n, one basis, one pool. |
| M8 | **Never correct an ARGMIN comparison with a mean-shift constant measured per member.** Rebuilding changes *which* item wins; measured here at **7× the per-member figure**. | Same class as quoting an MDE as a constant. |
| M9 | **Hypotheses closed against their proposers.** Workstream D registered and then refuted **six** of its own: R1 (set diversity), D7 (distance functionals), D1's mechanism, D8's mechanism, Q1's interval, D9's magnitude, and T2 (error coherence). **In every case the instrument the hypothesis motivated produced the result.** | Raising a decisive hazard and being wrong about its size are different things, and only the first changes the answer. |

---

## THE SPRINT'S ARITHMETIC, IN ONE PLACE

| avenue | measured value | status |
|---|---|---|
| search / optimiser / budget | exhaustive argmin ties the pool (+0.085) | **CLOSED by exhaustion** |
| ansatz, χ = 1 → 48 | spans 0.37 Å, gradients rise with depth | **CLOSED** |
| training the circuit at all | +0.003 [-0.136,+0.144] | **null** |
| training Hamiltonian | +0.056 [-0.020,+0.129] | **null** |
| torsion encoding (θ vs sin/cos) | 0 at matched step count | **NOT SUPPORTED** |
| generative latent as a source | **-0.103 Å** at the ORACLE of the union | **CLOSED by arithmetic** |
| **readout Hamiltonian** | **-0.697 [-1.059,-0.352]**, 5/5 folds | **the one non-null lever** |
| **in-pool selection** | **-1.338 [-1.558,-1.138]** ceiling; the **consensus family** captures **0%** | **where the Ångströms are — unclaimed, not proven unreachable** |

**Every avenue this sprint was asked to test is measured at or below 0.1 Å except the two that are
selection.** The 3.2 → 2.5 → 2.0 Å ladder does not require a better sampler, a better generator, or
more quantum resource. **It requires a better discriminator over candidates the pipeline already
retrieves** — and E4 locates the one place that discrimination demonstrably breaks, without yet
supplying a native-free way to detect it.
