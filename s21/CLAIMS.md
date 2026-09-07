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
| B8 | QNG on **AMBER** — **ANSWERED, and two-sided.** `qng_diag` is the **one** arm on the **one** Hamiltonian where geometry-aware preconditioning **closes the objective gap to Adam** (+0.037 ns vs +0.214 / +0.291) — exactly where its a-priori case was strongest — **and the reward is +0.496 Å of harm.** | 20 | — | **ESTABLISHED.** The method works; the working is the damage. |
| B9 | **On AMBER, optimising at all is significantly HARMFUL**, and `best_of_N` (5.053) is itself worse than the zero-evaluation pool mean (4.739). | 20 | **+0.523 [+0.248,+0.791]**, 4W/16L, **5/5 folds** | **ESTABLISHED** — Adam optimises the loss best and yields the worst structure |
| B10 | Harm orders monotonically with how well the energy tracks structure, across three Hamiltonians through identical machinery. | 20 | ρ(E,RMSD) −0.071 / +0.191 / +0.497 → harm +0.523 / +0.113 / +0.091 | **SUPPORTED by mechanism only.** *The third clause of `search-saturates` does NOT reproduce — nothing helps on any of the three* — but see B10a before calling it refuted. |
| B10a | Reconciliation with D5. | 60 / 20 | step-count **within** an optimiser on the distogram: **−0.087**; optimiser endpoint **vs best-of-N**: **+0.091** | **BOTH TRUE, different contrasts.** An optimiser can improve with steps and still converge worse than the best of N draws. Third clause carried **UNCONFIRMED**, not refuted. *Fifth reconciliation-by-naming-the-operator this sprint.* |
| **B11** | **THE READOUT EFFECT ON AMBER IS NEARLY DOUBLE LEGACY'S, and its size tracks how badly the readout Hamiltonian RANKS** (AMBER −1.30, Legacy −0.70). | 20 | **−1.299 [−2.002,−0.550]**, 17W/3L, **5/5 folds**, all seven arm-level CIs excluding zero | **ESTABLISHED** |
| **B12** | **THE SPRINT'S ONE PRESCRIPTIVE STATEMENT.** With the training-H null (B7) beside the readout-H effect (B6/B11): **if AMBER is to enter a CVaR-VQE at all, put it in the TRAINING role and take the argmin with something that ranks.** | 20 | two independent measurements, not one arm | **ADOPTED.** Preserves both mandated pillars — CVaR-VQE stays the selector; Legacy and AMBER stay separably evaluable, in **different roles**, which is what separability was for. |

> **B6 vs B7 is the sprint's cleanest architectural statement: the Hamiltonian inside the training
> loop is nearly free to change; the one in the readout is not.** It also means the mandated
> `H(λ) = (1−λ)H_L + λH_A` continuation operates in the dimension measured as null (B7), so any
> effect along λ must be claimed as an **optimisation-path** effect, with B7 quoted beside it.
> **C6 then shows its natural parameterisation is degenerate in raw units.** Both of the user's named
> Hamiltonian-design hypotheses are therefore closed or blocked — and **neither closure rests on a
> null alone; each arrives with a mechanism.**

---

## C. LEGACY vs AMBER, AND THE READOUT

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| C1 | Picking the **lowest-energy** structure is worse than picking at random, **for both physics energies**. | 20, subset(20) | AMBER 4.990 \| Legacy 5.487 \| distogram 3.676, pool mean 4.739 | **ESTABLISHED** (Sprint 20; provenance corrected — see M4) |
| C2 | The two distance functionals (**squared** vs **shipped Bayes risk**) are interchangeable at the argmin on pool-like structures. | 126 | +0.050 [-0.019,+0.156], MDE 0.124; median ρ **0.973** | **null** — and D's own opposite hypothesis **REFUTED** |
| C3 | That interchangeability is a property of the **pool's** structure distribution and does **not** transfer unchanged to the latent's. | 75 | +0.114 SE 0.069 on the latent, vs +0.050 on pool windows | **ESTABLISHED** — a measured null does not automatically transfer |
| C4 | Legacy's partial rank skill beyond the distogram. | 42 | **-0.0076** | **≈ zero** |
| C5 | **AMBER is genuinely a harder optimisation landscape**, on metrics invariant to the energy's units: condition number **6.26 → 7,756**, participation ratio **0.497 → 0.067**, anisotropy **4.18 → 18.19**, near-zero eigenvalue fraction **0 → 0.481**. Raw distribution: skew 8.13, kurtosis 68.9, **16% of pool members above 1e6 kcal/mol**, max **1.8e11**. | 30 | — | **ESTABLISHED** — Sprint 20's open question, answered |
| C6 | **The mandated `H(λ) = (1−λ)H_L + λH_A` is DEGENERATE in raw units.** The crossover sits at **λ\* = 1.99e-05** (median) with a **3.7e7-fold spread across targets** — so a uniform grid is AMBER-dominated above λ≈2e-5, every λ a person would pick (0.1, 0.25, 0.5) is *the same Hamiltonian*, and **no single grid is comparable between two targets**. | 30 | — | **must be run in NORMALISED units**; unknowable without measuring components separately first |
| C7 | **AMBER preconditioning by Legacy is NOT SUPPORTED, and its mechanism is refuted by its own diagnostics.** AMBER energy reached **1.55e5 → 1.74e8**, condition number **7,756 → 10,840**, near-zero fraction **0.481 → 0.503**, ‖grad‖ **2.0e7 → 2.4e10** — all worse — after moving θ by only **0.063 rad**. | 30 | RMSD +0.154 [-0.032,+0.431], MDE 0.348, 14W/16L | **NOT MEASURED** on RMSD; **mechanism REFUTED** |
| C8 | **THE MANDATORY 7-CELL MATRIX — DELIVERED.** Genuine CVaR-VQE, 4 seeds, matched candidate-evaluation budget, every component separably evaluable throughout. **Distance 3.845/3.686 wins on mean AND median; all six Legacy/AMBER cells are worse — 6 of 6 harmful** (Dist+Amb +0.273, Dist+Leg +0.274, D+L+A +0.561, Legacy +0.922, Leg+Amb +0.944, AMBER +1.471). | 12 | 4 CIs exclude zero; **every effect below its own MDE** | **DIRECTION ESTABLISHED; EVERY MAGNITUDE NOT MEASURED.** F-A1 does not fire. |
| C8a | **The matrix is SEED-NOISE-LIMITED**, and this is the binding limit on it. | 12 | ansatz-seed sd **0.51–1.06 Å** (2–5× Sprint 20's 0.200); best-of-4-seeds beats the 4-seed mean by **0.65–1.36 Å**, exceeding most row-to-row gaps | **Only the 6/6 DIRECTION is safe.** An operator fork the lane did **not** enumerate — found in coordinator fork review. |
| C8b | **FORK REVIEW (rule 0 clause 2, first invocation).** The lane declared six forks, named the direction each pushes — **including the one pointing AGAINST its own conclusion** (wall-clock budget) — and requested independent review. | — | — | **Forks 2/3/5/6 PASS** (basis measured not assumed; readout declines `tail_avg`, where Legacy looks *worst*; null avoids uniform-torus and the initialisation mean; budget matched in the selector's space). |
| C8c | **The AMBER row is the BARE single point, not the deployed `E ∘ Relax_50`** (Sprint 20 L7c: relaxation is constitutive). | 12 | deployed AMBER costs **+0.021 [+0.014,+0.028]**, and at k=30 **beats** the projection by −0.022 | **+1.471 MUST NOT be quoted as AMBER's deployed cost.** The lane flagged this fork itself. |
| C8c′ | **How much of the direction is independent of that fork — CORRECTED.** The coordinator's review asserted *"the four Legacy rows never use the AMBER functional"*. **They do:** `Leg+Amb` and `Dist+Leg+Amb` both contain AMBER. | 12 | AMBER-free rows: **only Legacy +0.922 and Dist+Leg +0.274** | **Direction is 2/2 independent of the AMBER fork and 4/4 CONDITIONAL on it** — materially weaker than the review credited. Asserted, not counted. |
| C8d | **A pre-registered rule printed and not applied.** The lane declared *"if a non-declared normalisation wins, the declared choice is recorded as WRONG"*. | 12 | Matched at seed 0: **declared wins 3 of 4**; only `Dist+Leg+Amb` goes to `rank` (4.224 vs 4.363) | **DECLARED CHOICE RECORDED WRONG on that one cell** — not swapped, and 0.139 Å is inside seed noise, so nothing is promoted either. *Invoking low power after the fact does not release a rule declared before the run.* |
| C8d′ | **The coordinator's own load-bearing check was mis-computed — CORRECTED by the lane.** The audit table compared **seed-0 audit arms to a 4-SEED MEAN** on an instrument with 0.51–1.06 Å seed noise; the review repeated that arithmetic without checking the comparison was matched. | 12 | Matched at seed 0: best-normalised hybrid **3.866** vs Distance-only **3.807** = **+0.059 Å, NOT MEASURED** (the review claimed ~0.24 Å) | **Distance is NOT SHOWN to beat the best-normalised hybrid; it is merely not beaten by it.** C8's verdict rests on the **declared** comparison (0.27–0.94 Å), not on this reconstruction. |
| C8e | **Cross-readout — the priority cell — is a NEGATIVE.** Best off-diagonal `Distance→Dist+Amb` 3.922 vs best diagonal `Distance→Distance` 3.845. | 12 | **+0.077 [−0.187,+0.348]**, 5W/7L | **NOT MEASURED, worse.** Reconciles with B6/B11: **B swaps against a BAD readout; A measures against the BEST one.** |
| **C8f** | **THE PRESCRIPTION, CORRECTED BY THE CONJUNCTION.** B12 said *if* AMBER enters a CVaR-VQE, put it in the training role. C8 answers the prior question. | 12 | Distance alone beats all seven cells; no cross-readout recombination reaches it | **NEITHER PHYSICS ENERGY SHOULD ENTER THE SELECTOR AT ALL.** Quoting B12 without C8 turns a conditional into a recommendation it is not entitled to be. |
| C8g | **F-A2 is NOT passed**: VQE minus best-of-N from the **untrained** circuit. | 12 | −0.361 to +0.108 on all seven | **NOT MEASURED / OPEN.** A **scope restriction** on Sprint 20 Q4 (8192 evals, n=126), **not a refutation** — the lane's own distinction, and correct. F-A3's kill condition also holds. |
| C9 | ~~The mandatory 7-cell matrix~~ (Legacy \| AMBER \| Legacy+AMBER \| Distance \| Distance+Legacy \| Distance+AMBER \| Distance+Legacy+AMBER), separably evaluable. | 4/12 targets | — | **INCOMPLETE — Workstream A, ~90 min outstanding** |

---

## D. THE ENCODING LEVER — NOT SUPPORTED

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| D1 | At **exactly matched step counts** the θ vs (sin θ, cos θ) encoding effect is **zero**, and on DIST/spsa that null is **properly powered** (own MDE 0.052, below the 0.084 bar). | 60 | -0.015 [-0.050,+0.022] | **NOT SUPPORTED** |
| D2 | The Sprint-20 result was a **step-count confound**: at fixed budget the embedded parameterisation has twice the dimension, so adam takes 4.6 steps against 9.6. | 60 | AMBc/spsa shrinks from -0.649 [-1.205,-0.165] "sig" to **-0.295 [-0.636,+0.022]** | **ESTABLISHED** |
| D3 | **The symmetry argument.** At matched budget the effect **flips sign with the objective** — LEG -0.165 [-0.273,-0.063], DIST +0.100 [+0.012,+0.189], both CIs excluding zero. *A change of coordinates that provably does not change the physics cannot help on one energy and hurt on another.* | 60 | — | **ESTABLISHED** |
| D4 | **And the step-count response PREDICTS the effect's size.** DIST residual **+0.004 against a predicted +0.087**. | 60 | LEG residual -0.039 [-0.122,+0.025]; DIST +0.004 [-0.023,+0.027] | **ESTABLISHED** — sufficiency, not merely impossibility |
| D5 | "Optimising harder hurts" is a property of the **LEGACY objective**, not of the instrument. The deployed distogram objective goes the **other way**. | 60 | LEG **+0.091 [+0.078,+0.109]**; DIST **-0.087 [-0.168,-0.018]** | **ESTABLISHED** — corrects how `search-saturates-discrimination-binds` was being quoted, including by the coordinator |
| D6b | **THE MECHANISM, identified analytically.** A gauge sweep maps the two parameterisations onto each other through angular displacement: `emb_r05 ≈ th_s2` (‖Δθ‖ 3.68 vs 3.37, obj −1.485 vs −1.600) and `emb_r2 ≈ th_s05` (0.661 vs 0.649, −0.906 vs −0.893). **The circle chart's RADIUS and the angle chart's SCALE are the same knob** — the "encoding" is a **step-size reparameterisation**. | 20 | — | **ESTABLISHED** — the third and deepest route to D1 |
| D6c | Wrapping the torsion angle is **exactly the identity** for a gradient-state optimiser (Adam's state depends only on gradients; gradients at θ and θ+2π are identical). | 20 | **−0.000 [0.000, 0.000]**, wrap firing 1.5×/run | **REFUTED ANALYTICALLY**; only residue is float drift on SPSA (+0.039, ns) |
| D6 | A real **gauge defect** in the embedded arm: `‖u‖` grows 1.000 → 2.5 (max 4.03) along a direction whose true gradient is **exactly zero** (arctan2 is scale-invariant), so every unit of radial motion is injected finite-difference noise; `adam_fd` gives that pure-noise direction a full-size step. **`_EmbField` is running a different OPTIMISER, not a different set of COORDINATES.** | — | — | **BUG, real, and NOT load-bearing** — see D6a |
| D6a | **The cost of that bug, MEASURED rather than asserted.** The prescribed fix (retract ‖u‖ each step) was implemented. | 20 | lands within **0.015 Å** of the unfixed arm | **CORRECTION.** This register first carried *"fix before any embedded arm means anything"*. The diagnosis is right; the consequence is 0.015 Å. **Hygiene, not a blocker.** Sixth instance of a magnitude asserted rather than measured. |
| D7 | Whether **any** encoding effect exists. The AMBc panel's own MDE is ~0.53 Å and cannot exclude an effect of the originally reported size. | 10 | — | **NOT MEASURED**, with a 55%-smaller point estimate — *not* "no effect exists" |

---

## E. DISCRIMINATION — WHERE EVERY REMAINING ÅNGSTRÖM IS

| # | claim | n | effect | disposition |
|---|---|---|---|---|
| E1 | The **ORACLE-best latent configuration sits at median percentile 5.62** of the objective's own ordering (mean 20.27, max 99.91) against a **50.0** no-skill null. | 126 | <1% on 21/75, <10% on 42/75 | **ESTABLISHED** |
| E2 | **The objective's ordering skill is -1.184 Å at M=1 and gone by M=64.** The whole fall of the top-M ladder is **min-of-M**, not ordering. | 75 | M=1 **-1.184 [-1.448,-0.918]** 63W/12L; M=8 -0.202; M=64 +0.039; M=512 **+0.107 [-0.051,+0.300]** | **ESTABLISHED** — the D9 matched control; **overturned the coordinator's first reading** |
| E3 | Unconditionally, the ordered window's excess over the global best exceeds a random window's — recovering E2 from a construction with **no split and no selection effect**. | 126 | ordered +0.507 [+0.359,+0.673]; random +0.388 [+0.317,+0.460]; difference +0.119 [-0.038,+0.290] | **ESTABLISHED** (two routes to one number) |
| E4 | **The aggregate null is bimodal.** The objective contains the answer in its top-512 on **59%** of targets (random gives 6%) and is a **concentrated wrong region** on the other 41% — **3.26× further from the ORACLE than an arbitrary window.** | 126 | contained 74/126; failure cell excess 1.228 Å vs random 0.377 Å | **MECHANISM** — the failure cell's **sign is FORCED** by conditioning; only the **3.26× magnitude** is free |
| E5 | An **averaging** readout extracts 0.45 Å from the same ordering an **argmin** readout finds worthless — because the terminal is dominated by the set MEAN (`d_out = 0.803·d_set_mean + 0.298·d_set_best` on WIDE sets like this latent, R² 0.960; the narrow m=75-over-K=500 fit 1.16/0.04 understates set-best ~7× here). | 126 | -0.449 [-0.520,-0.348], 89W/31L, 5/5 folds | **ESTABLISHED** — **no statement about an objective's usefulness is well-posed until the readout is fixed** |
| E6 | Whether a **native-free** signal from the objective's own SCORE DISTRIBUTION can separate E4's two regimes. Eleven predictors. | 75 | +0.2775 against a permuted best-of-11 bar of **0.3170**; the raw-AUC winner (0.776) was **length predicting length** — a binarised label whose threshold is `512/2**n` | **NOT DEMONSTRATED.** *The objective cannot audit itself.* |
| **E7** | **A COMPACTNESS channel — `Rg² = (1/N²) Σ d̂²` in closed form from the shipped distogram — DOES separate them.** Six predictors, D's exact protocol, partialled on length **and** difficulty. | 75 | **`rg_z` +0.3788, CI [+0.156,+0.557]**; bar 0.2899. Pooled **sceptic's family of all 17 tried tonight**: bar 0.3341 — **still clears**. Top **four** arms all Rg (0.379/0.372/0.371/0.306); best non-Rg 0.278; the partly-circular arm ranks 7th of 17. | **DEMONSTRATED as a SIGNAL** — falsifier fired, registered estimate (`\|ρ\|<0.30`) beaten. Difficulty control **consumes the native**, i.e. stronger than deployable. |
| E7b | **E7 attacked from both sides and SURVIVES; the mechanism is DISAGREEMENT, not extension.** | 75 | Length artefact **impossible by construction** (stratified effect is *larger*, +0.462 vs +0.379; ablating the i,i+1 term changes ρ by 0.003, r=0.9999). Native-free difficulty control (`pool_spread`) clears at 0.343/bar 0.286; both controls 0.360/0.287. Partialling **the pool's realised extension**: scale-free arms survive (`rg_z` 0.335, `rg_gap` 0.326, bar 0.280), **every scale-carrying arm dies**. | **PRIMARY IS NOW `rg_gap`/`rg_z`; `rg_disto` DEMOTED** — it falls to 0.281 under the proper control |
| E7c | Leave-one-length-stratum-out on `rg_z`. | 75 | +0.367 / +0.370 / **+0.289** / +0.513 / +0.387 | **HONEST FLAG**: all positive, none decisive, but the weakest sits **on** the bar (0.286–0.290). One 13-target stratum does real work. Fourth arrival of *an aggregate hides a structure* — this time it **prices** the finding. |
| E7a | **The Ångström value of E7.** | — | — | **NOT MEASURED.** ρ=0.38 is a correlation with *position in an ordering*, **not 0.4 Å**; the ρ→RMSD chain is exactly the composition the source memory calls unmeasured. Magnitude is a **max over the family, upward-biased**. Channel is independent of the objective's *scores*, **not** of the distogram. n≤13 panel. **Promotion: none.** |

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
| G4 | ~~m\* has not moved from 75, a forward prediction resolved negative~~ | 126 | — | **RETRACTED.** m\*'s argument is the **prior's MAE**, not the selector's quality (`m*=500 @ MAE 3.5 → 75-110 @ 2.34 → 20 @ 1.2`); the prior's MAE did not change this sprint, and the source already states *"m=75 is correct today"*. The law is also explicitly scoped **out of score-ordered gates**, which the m-ladder is. **The conclusion it supported is independently established by A2–A4, E2 and F1.** |
| G4a | What the ladder *does* show: **m=75 remains correct at the current prior MAE.** | 126 | m=150 +0.023, m=20 +0.042, both flat | **CONFIRMATION**, unsurprising |
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
