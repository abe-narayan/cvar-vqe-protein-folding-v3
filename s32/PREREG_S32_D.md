# PREREG — S32 lane D — PHYSICAL RESPONSE AND DYNAMICS

**Registered before any number of this lane exists.** Commit this file, then cite its hash in
`s32/LEDGER.md` beside every lane-D row. Contract rules 8 and 19 bind.

Endpoint: **mean built-chain Cα RMSD, `tuning126`, n = 126. Production = 3.2105 Å.** Every number
below names its basis. Cloud (3.0483) and set mean (3.5507) are different objects and are never
differenced against the chain.

---

## PART 0 — D0: THE G1 ESCAPE, DERIVED BEFORE ANY COMPUTE

The lane may not spend compute on an observable until it has named the hypothesis of **G1** that the
observable breaks. This part does that, and in doing so it **closes more of charter §§23–24 than it
opens** — which is the point of doing it first.

### 0.1 G1, stated exactly as S30 proved it

From `s30/REPORT_S30.md` line 690, quoted verbatim:

> A **non-chiral single-structure channel** can carry information the distance map does not …
> **G1** (classical MDS): `G = −½JD²J = XXᵀ`, so `D` fixes the centred coordinates up to `O(3)`;
> after proper-rotation Kabsch the residual ambiguity is **exactly reflection**. Hence `S` is a
> function of `D` **iff** `S` is reflection-invariant.
> Scope: *"**Single-structure** channels only, sequence held fixed."*

And its companion, `s30/REPORT_S30.md` line 692 (Corollary G1b):

> this project has exactly three [references] — the predicted distogram …, the pool (typicality) …,
> and universal physics (**target-independent by construction**, so it cannot supply a per-target
> sign)

G1 therefore rests on **four** hypotheses, and an observable escapes only by breaking one by name:

| | hypothesis | how an observable could break it |
|---|---|---|
| **H-a** | the observable is **reflection-invariant** (achiral) | build a chiral functional |
| **H-b** | the observable is a function of **ONE structure** | use a pair, a path, or the pool |
| **H-c** | the observable is **SE(3)-invariant and scalar** | emit an **equivariant vector** instead |
| **H-d** | *"the distance map"* is the map the pipeline **already consumes** | use a **finer** map than the one consumed |

**H-d is not part of G1's proof; it is part of its *application*.** G1 proves `S = g∘D`. The step
from there to *"so it adds nothing"* requires that the pipeline already has `D`. That step is sound
for the Cα distance map of a candidate. It is **not** sound if the observable reads a strictly finer
map than `D_Cα` — which is exactly what an all-atom force field does.

### 0.2 What I close, before spending anything — Theorem D-E (the ensemble extension of G1)

Charter §23 lists: molecular dynamics, short trajectory ensembles, conformational covariance, basin
transitions, state populations, metastability, transition rates, autocorrelation, dynamic modes,
transition-path structure, free-energy landscapes. Charter §24 adds ensemble free energy, ensemble
reweighting and temperature-dependent response. **S31 closed only the single-structure half of this
by derivation and explicitly recorded that "no thermodynamics was actually computed."** The
following closes the *stochastic* half by the same route, and it is new.

> **Theorem D-E.** Let `U` be an achiral potential — `U(Rx + t) = U(x)` for every `R ∈ O(3)`,
> `t ∈ R³`. (S31 Lemma B1 **verified** this for the shipped `amber14/protein.ff14SB.xml +
> implicit/gbn2.xml`: all 1340 `PeriodicTorsionForce` phases are at distance 0.000e+00 from {0, π}
> and there is no CMAP.) Let the dynamics seeded at `x` be Langevin, Brownian, Monte-Carlo or
> Hamiltonian with an isotropic Maxwell–Boltzmann momentum draw, i.e. any propagator built from
> `−∇U` plus noise whose law is `O(3)`-invariant. Let `A` be any achiral SE(3)-invariant observable
> of a configuration (or of a whole trajectory, applied configuration-wise and reduced by an
> invariant functional). Then
>
> `p_t(R y + t | R x + t) = p_t(y | x)`  for every `R ∈ O(3)`,
>
> because `−∇U` is `O(3)`-**equivariant** and the noise law is `O(3)`-**invariant**. Hence
>
> `⟨A⟩_{x,T,t} := ∫ A(y) p_t(y|x) dy`  satisfies  `⟨A⟩_{Rx+t,T,t} = ⟨A⟩_{x,T,t}`,
>
> so `⟨A⟩` is itself an **achiral, SE(3)-invariant, single-structure** functional of `x`.
> **By G1, `⟨A⟩ = g(D(x))`.** ∎

**Corollary D-E1 — stochasticity is not an escape.** A finite-sample MD estimate is
`Â = ⟨A⟩ + ε` with `ε` a sampling error whose law depends on `x` only through the same
G1-closed channel. **Noise is not information**: `Â` can only lose skill relative to `⟨A⟩`, which is
already a distance-map reading. *Running the thermostat does not create a fourth reference; it adds
variance to a third one.*

**Corollary D-E2 — the §23 list falls on one ground, not one at a time.** Ensemble averages, basin
populations, state populations, metastability lifetimes, escape/transition rates *out of a basin
seeded at `x`*, autocorrelation times, conformational-covariance spectra, dynamic-mode frequencies
and every `d⟨A⟩/dT` or `d⟨A⟩/dλ` susceptibility (derivatives in a **universal** parameter of a
G1-closed function are G1-closed) are **achiral invariant scalars of `⟨·⟩_x`** and hence functions of
`D_Cα(x)` whenever they are computed on the Cα representation. **Charter §23 is closed as a source of
new scalar ranking information, and §24's "ensemble free energy / reweighting / temperature-dependent
response" closes with it.** This extends S31 §8 from `F`, `E`, `S` to the entire stochastic family.

**Corollary D-E3 — a deterministic minimiser is closed twice over.** A trajectory from a
deterministic minimiser started at `x` is a deterministic function of `x`; any achiral invariant
reduction of it is an achiral invariant functional of `x`. (`LocalEnergyMinimizer` is deterministic
given the same build.) *"Run a minimiser and score the result"* does **not** escape G1. **The lane
must not build it.**

**Corollary D-E4 — the barrier survives, exactly as S31 left it.** `F‡(x_a → x_b)` is a function of
the **pair**, not of `D(x_a)` and `D(x_b)` separately. It breaks **H-b**. It is not closed.

### 0.2b Theorem D-F — contracting two equivariant fields falls back into G1

The natural way to build a *directional* physics observable is to contract the force with something
target-specific: *"does the AMBER force agree with the distogram's gradient?"* **It does not escape,
and the one-line proof is worth more than the experiment it prevents.**

> **Theorem D-F.** Let `F, G : R^{n×3} → R^{n×3}` both be `O(3)`-equivariant
> (`F(Rx + t) = R F(x)`, likewise `G`). Then `⟨F(x), G(x)⟩` is `O(3)`-**invariant**:
> `⟨F(Rx+t), G(Rx+t)⟩ = ⟨R F(x), R G(x)⟩ = ⟨F(x), G(x)⟩` for every `R ∈ O(3)`, reflections included.
> It is therefore an **achiral invariant scalar** and **by G1 a function of the distance map.** ∎

`−∇U_amber` is `O(3)`-equivariant because `U` is `O(3)`-invariant (Lemma B1). `−∇` of the shipped
distogram cost is `O(3)`-equivariant because that cost is a function of pair distances. So **force
alignment, gradient agreement, physics-vs-prior consistency, force-covariance scalars, and every
other scalar built by contracting equivariant objects are closed.**

> **Corollary D-G — the only surviving escape for physics is to emit a displacement and APPLY it.**
> Every *scalar* a force field can produce about one structure — energy, component energies, ensemble
> averages (D-E), curvature and normal-mode spectra (S31 B1′), free energy and entropy (S31 B1),
> susceptibilities, and now every equivariant contraction (D-F) — is a distance-map reading. What is
> **not** a distance-map reading is the **vector itself**, used as a correction rather than as a
> score. *Physics is closed as a ranker and open as a mover.* This is also exactly the shape of the
> prize: S31 §20.1 prices a **direction**, not a quality estimate, and says it *"does not have to be
> accurate — it has to point the right way."*

### 0.2c D0-C, answered from the code: **door 4 is CLOSED. AMBER has no resolution advantage here.**

Read rather than assumed (`core/amber.py`, `s8/generate_univ/<PDB>.npz`, `s24/d_harness.py`):

- A pool candidate stores **only** `W` (the Cα trace, `(k, n, 3)`) and `PHI`/`PSI` `(k, n)`.
  **There are no deposited all-atom, and no deposited backbone N/C/O, coordinates anywhere in the
  pool.**
- The AMBER path rebuilds everything: `geo.build_backbone_batch(φ, ψ)` → ideal-geometry `N,CA,C,O,CB`;
  `build_full_structure(seq, backbone)` adds **all sidechains at Engh–Huber ideal geometry with a
  modal rotamer**; hydrogens are **not** re-added per candidate but reconstructed from **frozen local
  frames** calibrated once per sequence. **There is no PDBFixer in the repository** (0 hits).
- Therefore the all-atom structure is `Ψ(seq, φ, ψ)` for a **deterministic** `Ψ`, and `(φ, ψ)`
  determine the ideal-geometry backbone exactly.

> **Every atom AMBER sees beyond the ideal backbone is a fixed function of the sequence and the
> backbone torsions. The sidechains carry no per-candidate information; the modal rotamer is the
> same rotamer for every candidate of a given sequence.** So `E_amber = g(D_ideal(φ,ψ), seq)` —
> an achiral (Lemma B1) invariant scalar of a structure that is itself a deterministic function of
> `(φ, ψ)`. **H-d is not broken. Door 4 is closed, and S31's "all-atom supports are one residue
> wider" is a locality statement about the same information, not an extra channel.**

This is the single cheapest closure in the lane and it removes the main stated reason to spend AMBER
CPU on scoring. **Registered consequence: lane D will not build an all-atom static score.**

### 0.3 The three doors that are actually open, and their prices

| door | breaks | object | prior art | status |
|---|---|---|---|---|
| **1. chirality** | H-a | a chiral scalar | S30-L26 (WRITHE +0.0405, 0.41× MDE, beaten by its own achiral twin); S31 §9 (`XTWIST` −0.0104, its achiral twin wins) | **measured empty at 9–16 residues, twice.** Not reopened by this lane. |
| **2. multi-structure / relational** | H-b | `f(x_a, x_b)`, a path, or a pool-referenced quantity | S31 §8 left `F‡` open **on price, not theory**; consensus/typicality (a pool-referenced *static* scalar) is in-band **−0.2837, the wrong sign** | **open.** The static pool reference is measured and bad; the *path* reference is unmeasured. |
| **3. equivariance** | H-c | a **vector field** `Δ(x) ∈ R^{n×3}`, not a scalar | none in the G1 family — all 43 S30 channels and every S31 channel are **scalars used for ranking** | **open, and this is the one the endpoint asks for.** |
| **4. resolution** | H-d | a functional of `D_all-atom`, strictly finer than `D_Cα` | S31's locality theorem notes all-atom supports are **one residue wider** | **open iff** candidate all-atom coordinates are not a deterministic function of the Cα trace. **D0-C below decides this in code, for free.** |

### 0.3b The chirality door, addressed directly — because two things circulating about it are wrong

The coordinator proposed chirality as *"the cheapest escape"* and offered two supports. **Both are
checked here against the code and the record before any compute, and both fail.** Registered before
the measurement so the correction is not retrofitted.

**(i) "AMBER is chiral (it distinguishes L from D)" — FALSE for the shipped force field, and S31
verified it by computation, not by assertion.** `s31/REPORT_S31.md` §8, Lemma B1, verbatim:

> the shipped `amber14/protein.ff14SB.xml + implicit/gbn2.xml` is **reflection-invariant** —
> ff14SB's **1340 `PeriodicTorsionForce` phases sit at distance 0.000e+00 from {0, π}**, and there
> is no CMAP.

The derivation behind that measurement: bonds, angles, Coulomb, Lennard-Jones and GBn2 are functions
of interatomic **distances**; a periodic torsion `k(1 + cos(nφ − δ))` is **even in `φ`** exactly when
`δ ∈ {0, π}`, and a mirror sends `φ → −φ`. With no CMAP there is no odd term anywhere in the
Hamiltonian. **`E(x*) = E(x)` identically**, so AMBER energy is an achiral invariant scalar and is
G1-closed — S31's Corollary B1, which is already in the record. S31's registered falsifier fired at
3.673e-06 and its **matched rotation control gave 5.377e-06, larger**: the residual is floating point,
not chirality. *AMBER is not the project's chiral scorer. It has none.*
**A mover, however, survives**: `Δ(x) = M(x) − x` is `O(3)`-equivariant, which is door 3, not door 1.

**(ii) "the distance objective is exactly mirror-blind, so a distance-based multi-start selects
enantiomers" — a real sentence in `core/project.py`, but the next sentence of the same docstring says
the shipped objective is not that one.** Verbatim, `core/project.py`:

> Related: the distance objective is exactly mirror-blind, so a distance-based multi-start selects
> enantiomers. **This objective is a COORDINATE distance and is chirality-sensitive, which is what
> makes lowest-objective selection safe here.** `test_project.py` asserts every emitted structure is
> L-handed anyway, because that safety is a property of the objective and would be silently lost if
> the objective were ever reformulated.

The first clause is a **warning about a counterfactual** objective. The deployed projection consumes
**coordinates**, is chirality-sensitive, and has an assertion in its test suite that the emitted
chain is L-handed. **The pipeline does not have a live chirality leak at the projection stage**, and
lane D will not build on the claim that it does. (The degeneracy that *is* live there is the
*Ramachandran-plausible vs implausible* torsion branch — a different object, and lane R's.)

**(iii) What is nevertheless right in the framing, and the cheap decisive test it implies.** It is
true that *"every native-free ranking signal this project has tested is a distance-map function and
therefore achiral."* But the project has measured the chiral escape twice and found it empty
(S30-L26 `WRITHE` +0.0405, 0.41× MDE, **beaten by its own achiral twin**; S31 §9 `XTWIST` −0.0104,
**its own achiral twin beats it** at +0.1124), and the standing explanation is the scope note
*"empty at 9–16 residues."* **That is a scope, not a mechanism.** The mechanism is testable, free,
and is registered here as **D0-X**:

> **A channel can only carry in-band information about a degree of freedom that VARIES in band.**
> Pool candidates are real protein fragments: L-amino acids, right-handed α, α_L rare. If the chiral
> degree of freedom is essentially **constant across the members of a target's top-75 band**, then no
> chiral functional — writhe, torsion-sign, Ramachandran region, or any other — can have in-band
> skill, however informative chirality is about proteins in general.

**D0-X test (free, no AMBER, no MD):** per target, over the top-75 band, compute the Cα-trace
pseudo-torsion sign profile and a scalar chiral summary; report (a) its **between-candidate sd within
the band** against (b) its **between-target sd**, and (c) the fraction of the band's total chiral
variance that is in-band. **Registered prediction: the in-band share is below 15%**, and if so,
*chirality's in-band emptiness is explained by degeneracy of the degree of freedom, not by chain
length* — which converts two sprints' scope note into a mechanism. **Falsifier:** in-band share
≥ 15%, in which case door 1 is genuinely open at this length and must be re-measured.

**Registered constraint carried from S31 §9 and from memory
`zero-information-control-must-be-plausible`:** any Ramachandran/torsion-region arm must beat the
**constant α-helix at (−57°, −47°)**, which S31 measured **beating every torsion channel on both
bands** (ρ +0.3413 vs `LEG_torsion`'s +0.1676; in band −0.0901 at 1.81× MDE, 5/5 folds). *"Beats
random" proves nothing here.* Lane D will not build a Ramachandran term.

### 0.4 Door 3 is the door, and the reason is S31 §20.1

S31's most consequential result says the missing quantity is **`mu`, the pool's common mode** — *"how
all of them are wrong together"* — priced at **−0.8102 Å ORACLE / NOT DEPLOYABLE**, and that
*"a future common-mode channel does not have to be accurate — it has to point the right way"*
(the shrink curve saturates at `c = 0.75`, 97% of the benefit).

**`mu` is a direction in `R^{n×3}`. It is an equivariant vector, not an invariant scalar.** Every
channel G1 closes is an invariant scalar. **G1 has never been applied to, and does not bound, an
equivariant vector-valued functional.** A relaxation displacement `Δ(x) = M(x) − x`, where `M` is any
`O(3)`-equivariant map (every force-field minimiser is one), is precisely such an object, and it is
**target-specific**, because `U` depends on the sequence through atom typing.

**The honest limit, stated before the measurement so it cannot be forgotten after it.** None of
doors 2–4 adds **Shannon information** relative to the candidate's own coordinates: everything here
is computable from data the pipeline holds. `Δ(x)` is, formally, a function of `D(x)` expressed in
`x`'s own frame. **What doors 2–4 escape is the *function class*, not the information.** G1's real
force is that the achiral-invariant-scalar class has been searched (43 channels, S30-L26) and is
empty. *The claim this lane may make is that a different class has positive in-band directional
skill — never that it found a new information source.* Any lane-D sentence claiming a "new
information channel" is, by this paragraph, wrong when written.

### 0.5 The reality check, registered as a falsifier rather than an excuse

These are **9–16 residue peptides**. Conformational basins, metastability and folding cooperativity
may not be well defined at that length. **That is measurable, and it is the lane's first experiment**
(D1-R below), because either horn kills the dynamics hypothesis *with a named mechanism*:

- if physics **contracts** the pool (candidates relax toward a common attractor), the dynamics
  destroys the very discrimination it was invoked to supply;
- if physics **does not move** the pool (candidates sit where they started), the trajectory carries
  no bit that the starting structure did not already carry.

Only an intermediate regime — candidates move *differently*, and stay distinct — leaves room for a
response observable. **D1-R measures which regime this instrument is in, before any ranking number
exists.**

---

## PART 1 — THE LADDER (charter §23: cheap first, escalate only on a derived reason)

Rungs run in order. **A rung is only run if the rung before it left it open.** One AMBER process at
a time (charter §32, contract).

| rung | name | door | cost | run condition |
|---|---|---|---|---|
| **D0-C** | code audit: is a candidate's all-atom detail a function of its Cα trace? | 4 | **free** | **DONE — §0.2c. Door 4 closed.** |
| **D1-E** | **AMBER** single-point **in-band** ρ on the real 126 × top-75 band, `Rg`-partialled | — | **free** (`s24/cache_amber`, 126×500 single points already on disk) | always. Charter §32's *"do not assume AMBER is a good ranker"*, measured in band on the shipped instrument. |
| **D1-L** | **Legacy**, all 11 terms, same basis, same band | — | ~1 min batched | always. Charter §33. |
| **D0-X** | the chiral degree of freedom's **in-band variance share** | 1 | free | always. Converts two sprints' *"empty at 9–16 residues"* scope note into a mechanism. |
| **D2-R** | the **basin-reality check**: does the pool contract, freeze, or neither under relaxation? | — | ~1–2 targets, minutes, 1 AMBER process | always |
| **D3-M** | **PHYSICS AS MOVER** — AMBER-relax the production structure; (a) **ORACLE** diagnostic `cos(Δd, mu)`, `cos(Δd, y)`; (b) native-free **built-chain endpoint**, paired, same job | 3 | 126 relaxations, ~40–60 min, 1 process | always — this is the lane's primary arm |
| **D4-P** | pairwise interpolation barrier `F‡(x_a → x_b)` between pool members | 2 | expensive | only if D3-M shows a positive `cos` AND D2-R leaves an intermediate regime |

**Rungs not on the ladder, and why — registered so they are not quietly added later:** any
**scalar** built from a relaxed structure, an MD ensemble, a normal-mode spectrum, a basin
population, a temperature derivative, **or a contraction of two equivariant fields** is closed by
**Theorem D-E, Corollaries D-E1–D-E3, Theorem D-F** and **will not be built**. Any all-atom static
score is closed by **§0.2c** and will not be built. Any Ramachandran term is closed by S31 §9's
constant-α-helix control and will not be built.

### D3-M, stated precisely, because it is the arm that counts

`mu` and `y` are **pair-distance** objects, not coordinate objects — `s31/s31_E_lib.py`:
`y = exp − d_nat`, `mu = pool75_mean − d_nat`, both in `R^{npairs}` over pairs with `min_sep = 2`.
The mover's output is converted into the same space: `Δd = d(x_relaxed) − d(x_production)`, one vector
per target in `R^{npairs}`. Then

- **(a) ORACLE / NOT DEPLOYABLE diagnostic.** `cos(Δd, −mu)` and `cos(Δd, −y)` per target. (`y` and
  `mu` are *errors*; a helpful move points **against** them, hence the minus signs, fixed here before
  the numbers exist.) **Registered prediction: median `cos(Δd, −y) ≥ +0.10`.** *This prices the only
  open door directly against the priced prize (−0.8102 Å along `mu`, ORACLE).*
- **(b) native-free endpoint.** Built-chain Cα RMSD of the relaxed arm against production, **projected
  in the same job**, paired, `folds`, `s24/stats_lib.compare` (LOWER IS BETTER). Prior art to beat and
  not to repeat: memory `averaging-space-beats-the-objective` records AMBER-relaxing the average at
  k=30 as **−0.022 Å [−0.036, −0.009], n=126, 5/5 folds** — real, right sign, **0.7% of baseline**.
- **CONTROLS.** (i) a **displacement-magnitude-matched random direction** in the same pair-distance
  space, matched to `‖Δd‖` per target — contract rule 6, matched to *this* arm's own norm, not
  another's; (ii) the **G1-closed twin** `‖Δd‖` used as a score, which must NOT beat the directional
  arm; (iii) `Rg` partialled out of every ρ (memory `physics-ranks-real-geometry-not-lattice`: the
  project's most expensive physics confound).
- **FALSIFIER.** `median cos(Δd, −y) ≤ 0`, or the endpoint arm at `< 0.7×` MDE. If `cos ≈ 0`, the
  mechanism is named: **the physical force is orthogonal to the pool's error, so physics-as-mover
  cannot reach the prize however it is scaled** — and door 3 closes with the rest.

---

## PART 2 — HYPOTHESIS / MECHANISM / PREDICTION / FALSIFIER / CONTROL / TEST / RULE / DEPLOYMENT

### H-D1 (the lane's hypothesis, charter §22)

**HYPOTHESIS.** A candidate's **response to a physical perturbation**, read as an **equivariant
displacement** rather than as an invariant scalar, carries in-band directional information about
candidate quality that static geometry does not.

**MECHANISM.** A candidate drawn from a fragment library carries local backbone geometry that was
realised for a *different* sequence. Under the target's own `U_seq`, a candidate in the correct
basin is already near a local minimum of the sequence's own potential, so its relaxation
displacement is **small and locally distributed**; a candidate in a wrong basin carries sequence-
incompatible strain, and relieving it produces a **large, collectively organised** displacement.
The *magnitude* of that displacement is an achiral invariant scalar and is **closed** (D-E3); the
**direction**, contrasted against a target-specific external field the candidate itself did not
generate, is not.

**PREDICTION (registered, directional per contract rule 18).** The in-band Spearman ρ between the
response observable and true Cα-RMSD, computed **within each target's top-75 band and then averaged
over targets**, is **positive and ≥ +0.10**, against the shipped `DIS` in-band **−0.0262** and
`CONS` in-band **−0.2837**.

**FALSIFIER.** In-band ρ ≤ 0 on the screen set, **or** |in-band ρ| < 0.7× its own paired MDE, **or**
the observable's in-band ρ does not exceed that of its **magnitude-only twin** (the same relaxation,
scored by `‖Δ‖` alone, which is G1-closed) by at least 0.7× the paired MDE. *If the direction adds
nothing over the closed magnitude, the door-3 escape is decorative and the channel is closed.*

**CONTROLS** (contract rule 6 — each matched to the operator's own space):

1. **G1-closed twin**: `‖Δ(x)‖` — same relaxation, same cost, scalar-invariant reduction. This is
   the control that prices the *escape*, not the physics.
2. **Zero-information plausible control**: a **constant α-helix at (−57°, −47°)** as the external
   reference field, never uniform-on-the-torus (memory: `zero-information-control-must-be-plausible`).
   S31 §9 records this control **beating every torsion channel**, so it is a real bar.
3. **Compactness control**: partial `Rg` out. Memory `physics-ranks-real-geometry-not-lattice` records
   the project's single most expensive physics confound — AMBER's non-bonded + solvation reward
   compactness, so any physics-vs-RMSD correlation runs through `Rg` unless `Rg` is partialled.
   **Every lane-D in-band ρ ships with its `Rg`-partialled twin in the same row.**
4. **Tie handling**: where a signal ties, **average the outcome over the tied argmin set**; never
   `np.argmin` on a tied signal (memory: tie-breaking leaks the pool order).

**TEST.** In-band Spearman ρ against true Cα-RMSD inside each target's top-75 band, paired per
target, `s24/stats_lib.compare(..., folds=...)` on the frozen 5 folds.

**STATISTICAL RULE.** `MDE = 2.8016 × SE`, per comparison. **< 0.7× → NOT A RESULT. 0.7–1.0× → NOT
MEASURED. ≥ 1.0× with fold CI excluding zero and ≥ 4/5 folds agreeing → RESULT.** Median beside mean,
W/L beside both, SE beside every mean.

**DEPLOYMENT CONDITION.** An in-band ρ result is a **diagnostic**. It counts as a lane-D result only
when carried to the **built chain**, projected **in the same job** as its paired production rows, at
**≥ 1.0× MDE** with fold CI excluding zero and ≥ 4/5 folds — and with the geometry secondaries
(virtual-bond mean/sd, displacement from the structure it claims to repair) on **row one** if the arm
constructs rather than selects (contract rule 15).

### H-D2 (the reality check, charter §22's scope question)

**HYPOTHESIS.** Conformational-basin physics is **not operative** at 9–16 residues on this pool:
relaxation either contracts the pool toward a common attractor or leaves it essentially fixed.

**PREDICTION.** One of: (contraction) median pairwise Cα-RMSD among relaxed candidates ≤ 0.7× its
pre-relaxation value; or (frozen) median per-candidate Cα displacement ≤ 0.3 Å. **Registered before
the run, with both thresholds stated, so neither horn can be chosen after seeing the outcome.**

**FALSIFIER of H-D2** (= the condition that keeps the lane alive): neither horn fires — the pool
neither contracts below 0.7× nor freezes below 0.3 Å.

### H-D3 (the resolution door, D0-C)

**HYPOTHESIS.** A pool candidate's non-Cα atoms are **not** a deterministic function of its Cα trace,
so an all-atom observable reads a strictly finer map than `D_Cα` and breaks H-d.

**TEST.** Read the pool loader and the AMBER build path. If non-Cα atoms are **rebuilt** from the Cα
trace by a deterministic builder, H-D3 is **FALSE** and door 4 is closed by composition: the all-atom
structure is a function of the Cα trace, hence (achiral case) of `D_Cα`. If they are carried from the
source fragment's deposited coordinates, H-D3 is **TRUE** and door 4 is open.

**This test costs nothing and it decides whether AMBER is worth a single CPU-second in this lane.**

---

## PART 3 — LOGISTICS

- Results as JSON under `s32/results/`, sharded `.jsonl` rows, launched with `s26/jobrun.py`.
- Every emitted comparison appended to `s32/MULTIPLICITY.md` **as emitted**, marked
  registered/exploratory.
- ORACLE arms labelled **ORACLE / NOT DEPLOYABLE** on every occurrence.
- `benchmark60`, the frozen folds, the cluster assignments and the sealed splits are never opened.
- Cross-job chain comparison only with the printed bit-identity check on the shared `PROD` rows
  (contract rule 3).
- Lane D stages only its own files.
