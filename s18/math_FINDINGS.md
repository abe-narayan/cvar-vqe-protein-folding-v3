# SPRINT 18 — MATH WORKSTREAM FINDINGS

Owner: the mathematical definition of the degree-1 objective.

**Status: the branch is closed. This is a write-up, not a campaign.** On the coordinator's
instruction (2026-09-06) all further compute in this lane was stopped: the 126-target
coefficient build was killed at 40/126, the argmin-sensitivity arm at n = 29, the
lattice→continuum bridge check at n = 2. Those three are reported below **as smoke, with no
directional claim**, and are flagged everywhere they appear. The completed, sealed work is the
lattice algebra, the gauge test, the reference-measure sensitivity and the mandatory controls.

Everything is native-free except rows explicitly labelled ORACLE, which are post-hoc scoring.

| artefact | contents | status |
|---|---|---|
| `s18/results/math_lattice.json` | L1–L5, 19 enumerated targets, exact | sealed `4cb0a8dfe9b5376b` |
| `s18/results/math_mu.json` | L6, reference-measure sensitivity | sealed `4cb0a8dfe9b5376b` |
| `s18/results/math_control.json` | matched-random projection control | sealed |
| `s18/results/math_report.json` | aggregation | sealed |
| `s18/results/math_anova_selfcheck.json` | Control D + the locality theorem | sealed |
| `s18/results/math_hamil_decomp.json`, `math_hamil_prior.json` | §6, what `hamil`'s separability actually is | complete, 19/19 |
| `s18/results/math_bridge.json` | continuum↔lattice check | **PARTIAL, n = 2, SMOKE** |
| `s18/results/math_sens.json` | ∂θ*/∂d̂ | **PARTIAL, n = 29, SMOKE** |
| `s18/cache/anova_<pdb>_pool_512_16.npz` | persisted objective coefficients | **PARTIAL, 40/126** |

Modules: `math_lib.py` · `math_lattice.py` · `math_report.py` · `math_anova.py` ·
`math_iface.py` · `math_bridge.py` · `math_sens.py`. Interface contract with EXPERIMENT:
`s18/MATH_to_EXP.md`.

---

## 0. The answer, in the order it damages the hypothesis

1. **The object the brief's bridge derives is not the object that produced 2.411 Å, and it is
   null.** The residue-additive (RA) ANOVA has a certified argmin of **2.657 Å against the full
   objective's 2.661** — **−0.004 [−0.380, +0.307], 3W / 5L / 11 ties, folds 1/5**. That is
   **5 % of the instrument's 0.084 Å minimum detectable effect**. The 2.411 Å belongs to the
   strict Walsh weight-≤1 (W1) object, which the bridge does not produce.
2. **W1 is a statement about a bit encoding, not about the objective.** Relabelling which
   2-bit code names which torsion state — pure bookkeeping — moves W1's argmin to a 19-target
   mean of **2.862 ± 0.135 Å** over 24 relabellings, worse than the full objective on 23 of
   24, with the shipped labelling at the **0th percentile** of its own null. RA is **exactly
   invariant** to the same relabelling on 19/19 targets.
3. **The separability that motivated the branch belongs to a term the production objective does
   not contain.** `hamil` is 75 % an exactly-additive torsion prior. Measured here: the prior
   term is **1.0000** residue-additive, the distogram term **0.404**, the blend **0.9495** —
   and **ρ(RA(hamil), the pure torsion prior) = 0.986**. The degree-1 truncation of `hamil`
   *is*, to three nines, the retrieval-pool torsion prior. `s17/refine.py`'s objective has no
   prior term at all.

Any one of the three closes the transfer. They are independent.

**Independent reproduction.** The ADVERSARIAL lane reached these numbers separately, from a
different implementation: their gauge mean is **2.852** against my **2.862**, their
RA − full **−0.004 [−0.390, +0.318]** against my **−0.004 [−0.380, +0.307]**, and RA gauge-
invariance holds 19/19 in both lanes. The QUANTUM lane's `s18/results/q_report_anova.json`
agrees with the RA and W1 point estimates to 15 digits. Three implementations, one answer.

---

## 1. THE THEOREM (claim: **EXACT** — this is a theorem, not a discovery)

### 1.1 Statement

Let a configuration be `x = (x_1,…,x_n)` with `x_r ∈ {0,…,k−1}`, `k = 2^m`, encoded in `nm`
bits with residue `r` owning the bit block `bits(r) = {mr, …, mr+m−1}`. Let μ be uniform.
Write the Walsh basis `χ_S(x) = ∏_{q∈S} (−1)^{x_q}`, `S ⊆ [nm]`; it is orthonormal in L²(μ).

> **Theorem.** Under uniform μ,
>
> **(a)** the order-≤1 functional ANOVA in the **qubit** factorisation is the orthogonal
> projection onto `span{χ_S : |S| ≤ 1}` — the Walsh weight-≤1 projection, **W1**;
>
> **(b)** the order-≤1 functional ANOVA in the **residue** factorisation is the orthogonal
> projection onto `span{χ_S : S ⊆ bits(r) for some single r} ∪ {χ_∅}` — that is, weight 0,
> **all** weight 1, and **exactly those weight-2 coefficients whose two qubits lie in the same
> residue**. Call it **RA**.
>
> For `k = 4` (`m = 2`), `W1 ⊊ RA`, and the gap is precisely the intra-residue weight-2 mass.

### 1.2 Proof

For a product measure `μ = ⊗_r μ_r`, the order-≤1 ANOVA `E₀ + Σ_r (E_μ[E|x_r] − E₀)` is by
construction the L²(μ)-orthogonal projection onto
`V₁ = {c + Σ_r g_r(x_r) : g_r ∈ L²(μ_r)}` (standard; the centring `E_{μ_r}[g_r] = 0` makes the
summands mutually orthogonal, and the residual is orthogonal to every element of `V₁`).

Under uniform μ, `{χ_S : S ⊆ bits(r)}` is an orthonormal set of `2^m = k` functions, each of
which depends on `x_r` only. The space of all functions of `x_r` alone has dimension exactly
`k`. Therefore

    {functions of x_r alone} = span{χ_S : S ⊆ bits(r)}   —   exactly, not approximately.

Hence `V₁ = span{χ_∅} ⊕ ⊕_r span{χ_S : ∅ ≠ S ⊆ bits(r)}`, a direct sum of mutually orthogonal
Walsh spans (the index sets `S` are distinct). The orthogonal projection onto a span of basis
elements keeps exactly those coefficients, so **P_RA zeroes precisely the coefficients whose
support is not contained in a single residue block**. For `m = 2` a subset of one block has
size 0, 1 or 2, giving (b). Taking the residue to *be* the qubit (`m = 1`) gives (a). ∎

### 1.3 Verification procedure and residual

`s18/math_lattice.py`, all 19 exhaustively enumerated targets, run on the rank-uniformised
objective **and** the raw one, comparing the ANOVA computed by tensor contraction
(`math_lib.anova1_discrete`) against the coefficient-zeroing projection computed by fast
Walsh–Hadamard transform (`math_lib.walsh_truncate`, `walsh_residue_additive`):

| identity | worst residual, relative to max\|E\| |
|---|---|
| qubit-ANOVA ≡ W1 | **8.9e-16** |
| residue-ANOVA ≡ weight-0 + weight-1 + intra-residue weight-2 | **1.1e-15** |
| orthogonality `⟨E − PE, g⟩ = 0` for every additive `g` | **5.5e-17** |
| variance budget `Var(E) = Var(PE) + Var(E − PE)` | **5.0e-16** |

(The ADVERSARIAL lane verifies the same identity at 1e-12 by a different route.)

### 1.4 The correction to `s18/BRIEF.md` §4

The brief's *"When μ is uniform on the enumerated lattice this **is** the Walsh weight-≤1
projection"* is **false as stated**, and the coordinator has corrected it in place (ledger L5).
The true statement is §1.1(b). Nothing in this lane was built on the false version — the two
objects were built separately and compared from the first run, which is what produced §0.1.

### 1.5 How much the two objects differ

| | mean retained variance fraction, 19 targets |
|---|---|
| W1 | **0.613** |
| RA | **0.913** |

**W1 discards 33 % of the per-residue field itself.** The discarded part is the intra-residue
weight-2 mass Sprint 17 measured at 95.7 % of all weight-2 mass and then, in the same report,
truncated away.

---

## 2. WHY RA IS GAUGE-INVARIANT AND W1 IS NOT

This is the load-bearing point of the sprint and it is not visible unaided.

### 2.1 What the gauge is

The `k = 4` torsion states of a residue are **unsorted circular-k-means cluster centres**
(`torsion_lib2._class_library`, seed 1). For the GENERAL residue class they come out as

    state 0 = (+58°, +30°)   left-handed helix
    state 1 = (−102°, +135°) beta
    state 2 = (−68°, −38°)   alpha
    state 3 = (−120°, −9°)   bridge

and the binary encoding then declares that the **high bit** means "{L-helix, β} versus
{α, bridge}" and the **low bit** means "{L-helix, α} versus {β, bridge}". Nobody chose that.
It is the order k-means happened to emit. Re-running the clustering with a different seed, or
sorting the states by φ, or by population, gives a different assignment of the same four
physical conformations to the same four two-bit codes.

Formally: let `π = (π_1,…,π_n)` with each `π_r` a permutation of `{0,…,k−1}`, and define
`(T_π E)(x) = E(π(x))`. This is exactly a relabelling of the per-residue lookup table — the
lattice point `x` in the new labelling names the conformation that `π(x)` named in the old one.
**`T_π` changes no structure, no energy, no ordering of configurations, and no argmin
conformation.** Any object that is a property of the objective must therefore satisfy
`P(T_π E) = T_π(P E)`, i.e. its argmin must be the same *conformation*.

### 2.2 RA commutes with the gauge — proof

Two facts:
1. **`V₁` is `T_π`-invariant as a subspace.** If `g_r` is any function of `x_r`, then
   `g_r ∘ π_r` is again an arbitrary function of `x_r`. So `T_π V₁ = V₁`.
2. **Uniform μ is `T_π`-invariant** (a permutation of a finite set preserves counting measure);
   so is any μ whose per-residue law is permuted along with the labels.

An orthogonal projection onto a `T_π`-invariant subspace, with respect to a `T_π`-invariant
inner product, commutes with `T_π`. Hence `P_RA(T_π E) = T_π(P_RA E)` and the RA argmin is the
same conformation in every labelling. **EXACT.**

*Verified:* per-target standard deviation of the RA argmin RMSD across 24 stable-seeded
relabellings is **0.000e+00 on 19 of 19 targets**.

### 2.3 W1 does not — proof, and the concrete failure

`span{χ_S : |S| ≤ 1}` is **not** `T_π`-invariant. Take `k = 4` and residue `r` with bits
`b₁b₀`. The three weight-≥1 Walsh functions supported inside that residue are the three ways of
splitting four states into two pairs:

    chi_{b1}      : {0,1} vs {2,3}      weight 1   — KEPT by W1
    chi_{b0}      : {0,2} vs {1,3}      weight 1   — KEPT by W1
    chi_{b1 b0}   : {0,3} vs {1,2}      weight 2   — DISCARDED by W1

The three splits are on an equal footing physically; the encoding privileges two of them.
A permutation `π_r` that maps the third pairing onto one of the first two — for example the
transposition (1 2) — carries a weight-2 function to a weight-1 function and vice versa, so
`T_π` does not preserve the subspace and `P_W1(T_π E) ≠ T_π(P_W1 E)`. Only the 8-element
subgroup of `S₄` generated by the two bit flips and the bit swap preserves the split; the other
**16 of 24** permutations move it. **W1 is a function of the encoding, not of the objective.**

### 2.4 What that costs, measured

| | 19-target mean argmin RMSD (ORACLE scoring) |
|---|---|
| full objective | **2.661** |
| W1, the **shipped** labelling | **2.411**  ← the published number |
| W1, over 24 random relabellings | **2.862 ± 0.135**, range [2.600, 3.221] |
| relabellings worse than the full objective | **23 / 24** |
| gauge percentile of the shipped labelling | **0.000** |
| per-target sd of the W1 argmin across relabellings | **0.392 Å** (worst 1.46 Å) |
| RA, all 24 relabellings | **invariant, sd 0.000** |

> **The expected value of "W1 argmin" over the encoding is 2.862 Å — 0.20 Å WORSE than the
> objective it truncates. 2.411 is the best of twenty-five encodings.** ADVERSARIAL's
> independent gauge mean is 2.852, and their encoding-luck term W1(identity) − W1(gauge mean)
> is −0.440 [−0.844, −0.133], folds 5/5 — the only interval in the whole degree-1 story that
> excludes zero, and what it measures is the arbitrariness of a bit code.

### 2.5 The gauge obstruction also blocks any informative reference measure

A functional ANOVA — with orthogonal components and a closing variance budget — requires μ to
be a **product over the chosen factors**. `_factorises` confirms that neither the retrieval-pool
marginal nor the Ramachandran prior factorises over qubits on any of the 19 targets (a general
law on four states is not a product of two Bernoullis). So under either informative μ the qubit
decomposition is not an ANOVA at all: the per-qubit terms stop being mutually orthogonal, the
variance budget stops closing, and "zero the coefficients above weight 1" stops being the
projection. **W1 exists as an ANOVA only under the uniform lattice measure — the least
informative choice available.** RA is defined under all three.

---

## 3. The two objects, side by side

Certified argmins on the fully enumerated instrument — no search error, no Monte Carlo.
ORACLE scoring, post hoc.

| object | argmin RMSD | median | vs full, paired [95 % CI] | W/L/ties | folds | ρ(obj, RMSD) |
|---|---|---|---|---|---|---|
| full objective | 2.661 | 2.782 | — | — | — | +0.264 |
| **W1 (strict Walsh weight-≤1)** | **2.411** | 2.184 | −0.249 [−0.650, +0.093] | 7/5/7 | 4/5 | +0.153 |
| **RA (residue-additive ANOVA)** | **2.657** | 2.339 | **−0.004 [−0.380, +0.307]** | 3/5/11 | 1/5 | +0.172 |

Intervals are plain paired bootstraps. Fold-clustered intervals were also computed and are
**narrower** here (−0.152/+0.179 for RA) — with five clusters that bootstrap is not
trustworthy, so the **wider** interval is quoted throughout. Both are in the artefact.

**Rank-conditioning is not the confound.** Sprint 17 truncated `uniformise(hamil)` and
production rank-transforms nothing. It makes no difference: raw and rank-uniformised objectives
give **identical certified argmins on all 19 targets for both truncations** (their retained-
variance fractions do differ: 0.613/0.637 and 0.913/0.950). One candidate artefact eliminated.

---

## 4. Controls (both mandatory controls, lattice arm)

| | 19-target mean argmin RMSD |
|---|---|
| **zero-information reference** — uniform draw from the same space | 4.003 |
| **matched-random** — random Walsh subspace at the *same* retained variance (0.613), 24 draws | **3.552** |
| W1 vs matched-random | **−1.140 [−1.714, −0.561]**, 14W/5L |
| ORACLE ceiling — best structure in the space | 1.030 |

**This control cuts against the closure and is reported as prominently as the rest.**
Discarding 39 % of the objective's variance *at random* costs 0.89 Å; discarding the
*high-degree* 39 % costs nothing. **Low degree genuinely carries most of this objective's
useful content — degree is not a decoy.** What kills the branch is not that truncation is
damaging; it is that (i) the truncation the brief derives buys nothing, and (ii) the extra
0.25 Å the strict truncation appears to buy is encoding luck. Neither statement should be
quoted without the other.

---

## 5. Reference-measure sensitivity (kept, and framed as sensitivity)

The branch is closed; this section exists so a future reader asking *"was the null an artefact
of the reference measure?"* gets an answer rather than an assumption. Three native-free product
measures, RA under each, same 19 targets, certified argmins:

| μ | what it is | argmin RMSD | vs full, paired [95 % CI] | W/L/T | μ-weighted var frac | factorises over qubits? |
|---|---|---|---|---|---|---|
| **uniform** | counting measure on the torsion library — zero-information | 2.657 | −0.004 [−0.380, +0.307] | 3/5/11 | 0.913 | **yes** |
| **pool** | the target's own shipped top-75 per-residue torsion marginal — native-free, target-conditioned, the brief's "most informative legitimate choice" | 2.595 | −0.066 [−0.426, +0.211] | 3/4/12 | 0.443 | no |
| **rama** | the shipped leave-fold-out Ramachandran density of each residue's class — generic, sequence-conditioned, native-free | 2.640 | −0.021 [−0.372, +0.261] | 3/4/12 | 0.762 | no |

Every interval covers zero. The **argmin is identical across all three μ on 14 of 19 targets**,
even though the retained-variance fraction swings from 0.08 to 0.95 across targets and measures.
The largest point estimate any μ produces (−0.066 Å under `pool`) is **still below the
instrument's 0.084 Å minimum detectable effect**.

> **The null is not an artefact of the reference measure.** μ changes how much of the
> objective's variance the additive component captures; it does not change where the additive
> component's optimum sits. The honest object is the same under every defensible native-free
> choice.

---

## 6. Where the "separability" actually came from (claim: **ESTABLISHED**)

Sprint 17's motivating sentence — *"the deployed objective puts 61.3 % of its variance at Pauli
weight 1 and 93.0 % at weight ≤ 2; it is very nearly a separable per-residue field"* — is about
`hamil = 0.75·z(torsion prior) + 0.25·z(distogram Bayes risk)`. The torsion prior is
`−Σ_r log P_r(s_r)`: **exactly residue-additive by construction.**

Measured here, exact ANOVA on the full enumeration, 19/19 targets
(`s18/results/math_hamil_decomp.json`):

| term | residue-additive variance fraction |
|---|---|
| `z(prior)` alone | **1.0000** (exactly, on every target) |
| `z(distogram Bayes risk)` alone | **0.404** |
| `hamil`, the 0.75/0.25 blend | **0.9495** |

And the truncation is the prior (`s18/results/math_hamil_prior.json`):

| | mean over 19 targets |
|---|---|
| ρ(RA(hamil), the pure torsion prior) | **0.986** |
| ρ(W1(hamil), the pure torsion prior) | 0.799 |
| ρ(full `hamil`, the pure torsion prior) | 0.959 |
| certified argmin of the pure torsion prior (ORACLE) | **2.615** |
| certified argmin of RA(hamil) (ORACLE) | 2.657 |
| certified argmin of the pure distogram term (ORACLE) | 3.077 |

> **The degree-1 object on the 19-target instrument is, to ρ = 0.986, the retrieval-pool
> torsion prior.** Its argmin is the prior's argmin. It was never a distilled distance
> objective.

**And the production objective has no prior term.** `s17/refine.py::_obj` is
`Σ_p ((d_p − d̂_p)/sd_p)²` — the distance term alone, whose additive fraction is 0.404, not
0.95. The continuum degree-1 truncation therefore has nothing to reduce *to*. This is a third
mechanism, independent of the gauge argument and of the RA null, and it explains all three.

*(A note against my own number: `s14/vqe_hamil.py` prints "torsion-prior certified argmin
3.636". That is a **different object** — `s14.ladder.prior`, the held-out-database state
occupancy. The term inside `hamil` is `s14.retprior.state_prior(pdb, "top75")`, the retrieval
pool's marginal, which is what I measured at 2.615 over 19 targets / 2.878 over the nine of
`TARGETS9`. The record is not wrong; the two priors are not the same prior.)*

---

## 7. The continuous object

`s18/math_anova.py` + `s18/math_iface.py`; contract in `s18/MATH_to_EXP.md`. Delivered,
verified, and now parked with its coefficient cache at 40/126 targets.

### 7.1 Locality, and a corollary nobody had written down (claim: **EXACT**)

`d_ab` depends on **exactly** residues `{a+1,…,b−1}`, both torsions each, and nothing else —
verified to 1e-9 at n = 9 and n = 12, every pair, every residue, both angles. `ψ_a` rotates the
chain about an axis through CA_a, and `φ_b` only moves atoms placed after CA_b, so the boundary
half-torsions drop out. The programme's recorded locality theorem is confirmed at residue
resolution.

> **Corollary, EXACT.** The entire distance objective is a function of the **n − 2 interior
> residues**. `f_r ≡ 0` for `r = 0` and `r = n−1`, for the full objective and for both
> truncations. Neither the degree-1 argmin nor a full refinement can place the terminal
> residues at all — and the frozen metric scores them. On a 9-mer that is 2/9 of the chain.

It is identical across arms so it confounds no comparison, but any terminal-residue choice is a
modelling patch and must be labelled one. (`align_lib.raw_jacobian` records the same fact from
the geometry side as "the four inert torsions"; this is its objective-side statement.)

### 7.2 The questions the brief asked, answered

**Which physical terms contribute at each order.**

| order | content | fate |
|---|---|---|
| 0 | `E₀`, the mean distogram penalty of the reference ensemble | kept |
| 1 | per residue, the full 2-torus profile of how its torsions shift the *ensemble-mean* penalty of every pair spanning it | kept |
| ≥ 2 | every **co-operative** effect — that two residues' torsions must be chosen *together* to place a distant pair. A hairpin set by a compensating pair of turns has its entire signal here | **removed** |

**Do the coefficients depend on residue identity and on pair separation?** Yes, and in a way
that is forced rather than fitted. Residue identity enters only through `μ_r` when μ is
target-conditioned; the functional form is identical for every residue. Separation enters
twice — through `w_p = 1/sd_p²`, and through support counting, since residue `r` is spanned by
exactly `r(n−1−r)` pairs. **The field is therefore automatically strongest mid-chain and exactly
zero at the termini: a shape imposed by chain topology, not by data.**

**Marginalisation or conditional expectation?** Conditional expectation, strictly.
`E_μ[E | θ_r]` averages the *energy* over the other residues; nothing is normalised and no
probability is marginalised. Equivalently it is the L²(μ)-orthogonal projection onto the
additive subspace. (The genuine *marginalisation* alternative — the free energy
`F_r(t) = −β⁻¹ log E_μ[e^{−βE} | θ_r = t]`, a soft-min rather than a mean — is a different and
defensible object: it keeps the low-energy envelope instead of the average, and it is **not** an
orthogonal projection, has no variance budget, and introduces a temperature. It is one extra
accumulator inside `Target.fit` and was **not run**, because the branch closed. Recorded as
**OPEN** and costed, not as a rescue.)

**Does it change the landscape or merely discard correlations?** Those are the same statement.
`E_{≤1}` is a different function with a different minimiser; the correlations discarded are
exactly the ones that made the minimiser different. It also removes the search problem entirely:
an additive objective's global optimum is coordinate-wise and closed-form. Sprint 17's *"there
is no search problem here for a quantum device"* is, restated, *"the object being optimised was
almost a separable field"* — and §6 says which field.

**Is it still meaningful as a continuous function of torsion angles?** Yes. `f_r` is a genuine
function on the 2-torus, tabulated on a G×G mesh and carried by **exact trigonometric
interpolation** (2-D DFT) — periodic by construction, analytic, exact at the nodes,
differentiable in closed form. A polynomial basis would have been wrong; the variables are
angles. Verified: analytic gradients match central differences to rel. err ≤ 2.9e-9 on all
three arms.

**Interpretability, in one sentence:**

> `f_r(φ, ψ)` is the mean distogram penalty this target pays for residue r adopting torsions
> (φ, ψ), when every other residue is drawn from the reference ensemble.

### 7.3 Alternatives considered, and why they were not chosen

| alternative | verdict |
|---|---|
| **circular harmonics truncated at degree 1** (only cos θ, sin θ) | a *different* truncation — "degree-1 in angles" is not "order-1 in residues", and it is strictly coarser. Built as `harmonic_truncate(order)` and as the angle-additive arm `E_ang`; reported separately, never conflated |
| **angle-additive ANOVA** (order-1 over the 2n individual angles) | built and delivered. It is the honest continuous analogue of "finer than per-residue" — but it is **not** the continuous image of W1: the two lattice qubits index k-means clusters of the *joint* (φ,ψ) library, not φ and ψ. **There is no continuous image of W1.** |
| **per-residue marginalisation** (free-energy / soft-min) | genuinely different and defensible; not an orthogonal projection; not run (see §7.2). OPEN |
| **low-rank tensor decomposition** of the conditional tables | rejected: not an orthogonal projection onto an interpretable subspace, so no one-sentence physical reading of a factor. The brief's interpretability test decides it |
| **learned additive torsion potential** | rejected outright: fitting needs a target, and every native-free target available is a signal already closed |

### 7.4 Computability, and what makes the object cheap

By §7.1 only `r(n−1−r)` of the ~n²/2 pairs vary with `θ_r`. With **common random numbers**
across the whole mesh, the constant pairs cancel *to the last bit* — not merely in expectation —
so quadrature error enters only through pairs that genuinely depend on `r`.

And because `E = Σ_p w_p(d_p² − 2 d̂_p d_p + d̂_p²)` is **affine in (d, d²) at fixed θ**, the
cache stores conditional *distance* moments `E_μ[d_p|θ_r]`, `E_μ[d_p²|θ_r]` rather than energies.
Both truncations therefore rebuild for an **arbitrary** distance vector with no new chain builds
(`o.with_dhat(d)`, milliseconds). That would have made the coordinator's α-ladder free on the
degree-1 arms, and it is what made `∂θ*/∂d̂` analytic in §8.

### 7.5 The estimator, and the bug Control D caught

    h_r(t) = (1/S) Σ_s E(Θ_s with residue r set to t),   Ē = (1/S) Σ_s E(Θ_s)
    f_r = h_r − Ē,   E₀ = Ē    ⟹   E_{≤1}(θ) = Σ_r h_r(θ_r) − (n−1)·Ē

For an exactly additive `E = Σ_r g_r` this is **identically** `E`, at any S, because
`h_r = g_r + (Ē − ḡ_r^S)` and the `Ē` terms telescope. Re-centring each `f_r` on its own mesh
mean breaks the telescoping and leaves an O(σ/√S) constant behind. **That was the first
implementation and Control D caught it**: fields exact to 7e-15, `E_le1` off by 1.08 absolute.

### 7.6 MANDATORY CONTROL D (claim: **EXACT**, PASS)

A random band-limited additive function on `(T²)^n` (bandwidth 3, mesh 8–24, so interpolation
error is impossible too), with `g_0 = g_{n−1} = 0` to match the real object's terminal structure.

| check | result |
|---|---|
| fields recovered (centred) | **7.1e-15** absolute, 3.3e-16 relative |
| `E_le1` at 64 random **off-mesh** torsion vectors | **4.8e-14** absolute, **1.0e-15** relative |
| `E_ge2` at the same points (must be 0) | **4.8e-14** |
| `PASS_machine_precision` | **True** at n = 6/8/10, S = 64/256/512, G = 8/12/16/24 |
| the angle-additive object on the same residue-coupled synthetic | 32.2 — correctly **fails** |

The last row is the control's own control: an implementation that recovered a residue-coupled
function with an angle-additive model would be broken in the other direction.

### 7.7 Other verification

| check | result |
|---|---|
| `E_full` vs `s17/refine.py::_obj` at the same torsions | **2.3e-13** |
| gradients vs central differences (`E_full`, `E_res`, `E_ang`) | rel. err ≤ **2.9e-9** |
| `E_lambda(·, 1) − E_full`, `E_lambda(·, 0) − E_res` | **0.0 exactly** |

### 7.8 Lattice ↔ continuum bridge — **SMOKE, n = 2, stopped**

Enumerating the *production* objective over all `k^n` lattice configurations and comparing its
exact RA projection against the continuous deliverable built under the identical reference
measure (`μ = "lattice"`):

| target | Spearman(continuous, exact) over all configs | Pearson | RMS abs. error / sd | certified RA argmin: exact vs continuous |
|---|---|---|---|---|
| 1CS9 | 0.9991 | 0.9991 | 0.111 | 5.195 / **5.195** |
| 2MK7 | 0.9975 | 0.9977 | 0.078 | 0.944 / **0.944** |

The implementation reproduces the exact enumerated ANOVA's **ordering to ρ ≈ 0.998 and its
certified argmin exactly** on both targets; the residual 8–11 % of a standard deviation is
combined Monte-Carlo (S = 512) and mesh (G = 16) error, dominated by a small constant offset
(Pearson 0.999). **Two targets is smoke and no conclusion rests on it.** Also visible here and
consistent with §6: the RA projection of the *production* objective retains only **0.41–0.48**
of its variance, against 0.913 for `hamil`.

---

## 8. ∂θ*/∂d̂ — the coordinator's robustification question — **SMOKE, n = 29, stopped**

The hypothesis (`s18/COORD_FINDING.md`): if the distogram's harmful error acts through the
*pairwise* structure the objective induces, the additive projection should have an optimum that
moves *less* when d̂ is perturbed. That is `∂θ*/∂d̂`, exactly computable here and reading no
native:

    ∂θ*/∂d̂ = −H⁻¹ ∂²E/∂θ∂d̂ ,   H = ∇²_θ E(θ*)
    full:  J = 2 H⁻¹ Gᵀ diag(w),  G = align_lib.pair_jacobian
    RA:    H block-diagonal 2×2 per residue;  ∂²f_r/∂t∂d̂_p = −2 w_p ∇_t E_μ[d_p|θ_r=t]

pushed through `align_lib.sup_jacobian`, so the statistic `S = ‖J_sup J‖_F/√n` is in Å of Cα
RMS motion of the optimum per Å of perturbation of d̂, superposition-invariant. Cutoffs
1e-8/1e-6/1e-4 agree (the only modes dropped are the four inert torsions).

At **n = 29 of 126 — a prefix of the alphabetical target list, not a random subset**:

| perturbation scaling | full | RA | angle-additive | RA − full, paired [95 % CI] | W/L | folds |
|---|---|---|---|---|---|---|
| unweighted | 1.409 | 1.817 | 1.935 | **+0.407 [+0.135, +0.684]** | 8/21 | 5/5 |
| sd-scaled (matched to the distogram's own sd) | 1.040 | 1.252 | 1.239 | **+0.212 [−0.040, +0.479]** | 12/17 | 3/5 |

Zero-information control (the same functional with d̂ shuffled across pairs): 1.333.

**No directional claim is made.** The two weightings disagree on whether the interval excludes
zero, the subset is not random, and n = 29 of 126. What can be said is that **at n = 29 there is
no sign of the predicted *reduction* in either weighting** — both point estimates go the other
way. The mechanism that would explain it, if it survived: the full objective is massively
over-determined (`n(n−3)/2` restraints on `2n−4` free torsions) and therefore stiff, so `H⁻¹` is
small and its optimum barely moves; the additive projection deletes the cross-residue coupling
and leaves each residue's two torsions set by a conditional *average*, which is flatter — median
Hessian spectra bear this out. **Averaging over the other residues damps the signal at least as
much as the error.** Labelled **OPEN**; finishing it would cost ≈ 90 min of the coefficient
build plus ≈ 18 min of fits, and it was stopped rather than finished because the branch closed.

---

## 9. Falsifier verdicts, in my own words

**F1 — degree-1 lands ≈3.6 Å on the real instrument.** *Not tested in this lane.* My workstream
never ran the 126-target refinement; the deliverable was built for EXPERIMENT and its
coefficient cache stopped at 40/126. I decline to report a verdict on someone else's
measurement.

**F2 — degree-1 is no better than the full objective. FIRED,** on the object the brief's bridge
derives. RA − full = **−0.004 Å [−0.380, +0.307]**, 3W/5L/11 ties, folds 1/5 — five per cent of
the instrument's minimum detectable effect, on the very instrument where the effect was born.

**F3 — the advantage dies under matched controls. FIRED, but only on the right control, and the
distinction matters.** Against a *matched-variance random projection* the strict object wins
decisively (−1.140 [−1.714, −0.561]) — truncating by degree is genuinely not the same as
truncating at random, and I will not let that be quoted away. Against the *gauge* control the
advantage evaporates and reverses: mean 2.862 over relabellings against the full objective's
2.661, shipped labelling at the 0th percentile.

**F4 — the effect exists only on the 19-target instrument. FIRED, and the mechanism is now
named rather than suspected.** It is not merely that 19 targets are few. The effect is carried
by (i) an arbitrary k-means cluster ordering and (ii) a torsion-prior term that constitutes 75 %
of `hamil`, is exactly additive by construction, and does not appear in the production objective
at all.

**F5 — the continuous mapping is mathematically invalid. SPLIT, and both halves should be
stated.** The *construction* is sound: the ANOVA bridge is a theorem verified to 1e-15, and its
continuous implementation reproduces the exact enumerated ANOVA (ρ ≈ 0.998, identical certified
argmin, n = 2 smoke). The *transfer of the published result* is invalid: there is no continuous
object corresponding to W1, because the lattice has no φ/ψ factorisation and the continuum has
no qubits, and under any informative reference measure W1 is not even an ANOVA. **F5 fires on
the result, not on the machinery.**

---

## 10. Claim register

| claim | label |
|---|---|
| RA = projection onto within-residue-supported Walsh coefficients (weight 0 + 1 + intra-residue 2); W1 = qubit-ANOVA order 1 | **EXACT** (§1, residual 8.9e-16) |
| the projection is orthogonal and the variance budget closes | **EXACT** |
| RA commutes with any per-residue state relabelling; W1 does not | **EXACT** (§2.2–2.3), and **ESTABLISHED** empirically (sd 0.000 vs 0.392 Å) |
| `d_ab` depends on exactly residues {a+1..b−1}; the objective is a function of the n−2 interior residues; `f_0 = f_{n−1} = 0` | **EXACT** |
| Control D: an exactly additive objective is recovered to machine precision | **EXACT** |
| `E_lambda(·,1) ≡ E_full`; `E_full ≡ s17/refine.py::_obj` | **EXACT** |
| the 2.411 Å result is W1, not RA | **ESTABLISHED** |
| the torsion-prior term of `hamil` is exactly additive (1.0000); the distogram term is 0.404; ρ(RA(hamil), prior) = 0.986 | **ESTABLISHED** |
| RA is null against the full objective on the 19 enumerated targets | **SUPPORTED** — −0.004 [−0.380, +0.307], n = 19 with 11 ties; an interval this wide cannot establish a null, only fail to find an effect |
| the null is not an artefact of the reference measure | **SUPPORTED** |
| discarding variance *at random* is far worse than discarding it by degree | **ESTABLISHED** |
| the degree-1 optimum is *more* sensitive to d̂ than the full optimum | **OPEN** — smoke, n = 29, non-random subset, weighting-dependent |
| the continuous implementation reproduces the exact enumerated ANOVA | **OPEN** — smoke, n = 2 |
| free-energy (soft-min) per-residue marginalisation | **OPEN**, unrun, costed |
| all argmin RMSDs, `S_dir_ORACLE`, `resid_rms` | **ORACLE** |

---

## 11. What I got wrong, kept here rather than deleted

1. **The first estimator re-centred each `f_r` on its own mesh mean.** It looks harmless — the
   ANOVA centring convention only moves constants between `E₀` and the fields — but it destroys
   the telescoping identity that makes the estimator exact on additive functions and left
   `E_le1` biased by an O(σ/√S) constant. **Control D failed on it.** That is the mandatory
   control earning its cost in the first hour, and it is the reason the control exists.
2. **I first reported only a fold-clustered bootstrap interval**, which on 19 targets in five
   folds came out *narrower* than the plain paired bootstrap (−0.152/+0.179 against
   −0.380/+0.307 for RA). With five clusters that bootstrap is not trustworthy, and quoting the
   narrower of the two would have overstated the precision of a null. Both are now computed and
   the **wider** is quoted throughout.
3. **I suspected rank-uniformisation was the confound** — Sprint 17 truncated
   `uniformise(hamil)`, production rank-transforms nothing. It is not: raw and uniformised give
   identical argmins on all 19 targets. Recorded because it was a live alternative explanation
   that had to be closed rather than assumed away.
4. **I sized the continuous build without measuring first.** 126 targets at S = 512, G = 16 is
   ≈ 2 hours on a five-process box; I launched it and then discovered the rate. Had I measured
   one target before launching, I would have delivered the sensitivity arm complete instead of
   at n = 29. The lattice work — which is what closed the branch — needed none of it.
