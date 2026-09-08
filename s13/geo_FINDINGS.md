# s13 TRAINABILITY-AND-GEOMETRY — findings

**Question.** How does the complexity of the molecular energy model reshape the geometry of
the variational state manifold and the trainability of the VQE?

**Controlled comparison.** Legacy (`core.energy`, 11-term knowledge-based) vs AMBER
ff14SB/GBn2 (`core.amber`, genuine, single point, no minimisation), with representation,
ansatz, initialisation, CVaR α, optimiser, seed **and objective-evaluation budget** matched.

Files: `s13/geo_common.py` (machinery), `s13/geo_build.py` (energy tables),
`s13/geo_audit.py` (E1), `s13/geo_anova.py` (E2), `s13/geo_grad.py` (E3),
`s13/geo_metric.py` (E4), `s13/geo_land.py` (E5), `s13/geo_cvar.py` (E6),
`s13/geo_opt.py` (E7), `s13/geo_pauli.py` (E5, Pauli spectra), `s13/geo_kernel.py` (E9,
ansatz-only gradient kernel), `s13/geo_all.py` + `s13/geo_report.py` (drivers).
Results in `s13/results/geo_*.json`; superseded raw-only Pauli sweep preserved at
`s13/results/geo_pauli_v1_rawonly.json`.

**Sections still filling at write-up** (files are being appended to by running jobs, so
cell counts in `geo_cvar.json`, `geo_land.json`, `geo_metric.json` and `geo_opt.json` may
exceed those quoted below): §7 (2–3 of 6 cells), §8 (3 of 6 cells), §9 (1 of 9 cells),
§5c/§5d (2 of 9 cells). Every n is stated where it is used.

---

## HEADLINE (read this first)

**The energy model does not reshape the variational manifold. It cannot: the Fubini–Study
metric is a property of `(ansatz, |ψ₀⟩, θ)` and contains no Hamiltonian — measured here as
an exact zero (`g(θ)` bit-identical across energy models at matched θ, all seeds). What the
energy model does is select *where on the manifold the optimiser goes*, and that difference
is real (metric differs by 118–189 % of its own scale between the two endpoints) but its
DIRECTION is not consistent across seeds, so I report no stable Legacy/AMBER conditioning
difference.**

Six results, in order of how much I would defend them:

1. **The chain "energy model → Pauli-weight spectrum → gradient variance" closes empirically
   with no free parameter**, on a genuine force field: median measured/predicted ratio
   **0.996 (Legacy), 0.998 (Legacy-SOFT and AMBER-SOFT)** over 78 cell×depth combinations.
   The approximation is dropping cross-covariances between Pauli strings; the ratio is the
   test of it. (§4e)
2. **AMBER is measurably less local than Legacy — but only after both are conditioned
   identically, and the effect is ~2×, not the 2.7× the raw comparison suggested.** Mean
   Pauli weight 3.015 vs 2.236, AMBER higher on **79/79** cells; ≥4-body variance share
   0.122 vs 0.069, AMBER higher on **12/13** cells. **Neither model is 2-local.** (§4c, §4d)
3. **A correction I had to make, and it reverses an earlier claim of mine.** Raw AMBER's
   Pauli spectrum is the `Binomial(m,½)` signature of a *delta spike*: the ten most extreme
   configurations carry a median **0.996** of the variance, and the mean weight lands on
   `m/2` to three digits. My first "AMBER's mean weight exceeds Legacy's on 14/14 cells" was
   measuring steric clashes. Independently reproduced from the Pauli-spectrum agent's
   report on my own tables; some **Legacy** cells are spiked too (up to 0.878). (§4b)
4. **CVaR with α ≤ 0.5 is exactly a steric clash filter for single-point AMBER.** `amber`
   and a monotone tail-compression of it are *literally the same objective* for α ≤ 0.25 —
   identical gradient variance, entropy, `p_top` and all three readouts to every printed
   digit — while differing by 16 orders of magnitude at α = 1. The dynamic-range pathology
   lives entirely in the tail CVaR discards. **This is not barren-plateau mitigation** and is
   not described as such anywhere here. (§7b)
5. **QNG does not rescue AMBER, and not because the metric is singular.** The metric is
   **full rank at every θ, n and depth measured**, with `g_ii = 0.2500` exactly and
   off-diagonal correlations of 0.008–0.037 that *shrink* with n: `g ≈ I/4`. At depth 1 it is
   exactly `I/4`, so QNG ≡ scaled gradient descent. corr(‖Δθ‖, FS distance) = 0.884. QNG ties
   GD on AMBER; λ = 1e-3 never bound. **REFUTED: "AMBER's Euclidean parameterisation matches
   its state manifold badly".** (§5b, §9)
6. **Cost locality does not explain gradient behaviour at depth 3 — but it does at depth 2,
   and that is a testable prediction.** The ansatz kernel `v(w)` is non-monotone (peaks at
   w=2) and at depth ≥ 3 a maximally global Pauli string has 0.88–1.06× the gradient variance
   of a single-qubit one with **no trend in n**. At depth 2 the ratio falls 0.862 → **0.078**
   over n = 6 → 14. So the energy models' weight difference has a lever only at shallow
   depth. (§6a)

**And the physics, kept separate from the objective, as the brief requires:** SPSA optimises
the AMBER objective **best of all five arms** (energy percentile 0.088) and returns the
**worst structure of all five** (3.414 Å, 0.333 Å worse than uniform random sampling), with
the RMSD rebuilt from the same bitstring whose energy is quoted. Across §3's exhaustive
scoring, ρ(E, CA-RMSD) changes sign *within* a model across cells (−0.537 → +0.472 Legacy,
−0.102 → +0.661 AMBER). **Nothing here is a dev24 candidate and no arm improves accuracy.**

**Budget parity**: every diagnostic pays exactly `2ⁿ` unique evaluations per model, the same
`2ⁿ`; every optimiser arm spent exactly its budget (256/256 on all 20 arm×model
combinations). Wall time differs 35×; evaluation count does not differ at all.

---

## 0. Method: why the register is enumerated, and what that buys

Every quantity in this study — CVaR value, its exact gradient, the Fubini–Study metric, the
Hessian, the optimiser trajectories — is a function of the **vector of energies over the
whole register**. So the register is enumerated once per (target, length L, states-per-
residue k, model):

    E[j] = model.energy(bitstring(j))    for every j in [0, 2**n),  n = L*log2(k)

**Budget parity is therefore exact and trivially auditable**: both models pay exactly `2**n`
unique evaluations, the same `2**n`, and nothing else. Wall time differs by ~35× (AMBER
7.0 ms/call measured warm vs Legacy 0.20 ms/call); the evaluation count does not differ at
all. Matching wall time instead would have handed Legacy 35× the evaluations — the error
the brief names as a correctness condition.

Downstream, `geo_common.TableModel` is a `budget.BudgetedEnergyModel` subclass that serves
the table and charges one evaluation per **unique** bitstring, so the optimiser arms use the
project's own accounting class, not a re-implementation.

**Representation.** `torsion_lib2.PerResidueTorsion` over
`library_for(seq[:L], k, exclude_seq=full_seq)` — the sequence-conditioned, leakage-safe
per-residue torsion library. Peptide length enters by truncating a tuning target to its
first L residues; the library for the prefix still holds the **full** target sequence out.

**Targets.** One per leave-fold-out fold from the 126-target s12 tuning instrument, two of
them FAIL18 members: 1A13 (fold 4), 1A1P (fold 1), 2BFI (fold 2, FAIL18), 1DEP (fold 0),
1ID6 (fold 3, FAIL18), 1CEK (fold 2). **benchmark60 and dev24 are untouched**; `dev_set(24)`
was enumerated and verified disjoint from the 126 tuning targets (intersection empty).

**Three crossing ladders** so that qubit count, peptide length and states-per-residue are
not confounded:

| ladder | cells | n_qubits |
|---|---|---|
| k=4 | L = 4,5,6,7,8 | 8, 10, 12, 14, 16 |
| k=2 | L = 8,10,12,14 | 8, 10, 12, 14 (**double the peptide length at equal n**) |
| k=8 | L = 4,5 | 12, 15 (**same peptide, more states**) |

---

## 1. E1 — instrument reproduction and the gradient audit  `s13/geo_audit.py`

**Status: PROVEN (verification, not a scientific claim).**
Config: 1A13/2BFI, L=5 k=4 (10 qubits) and L=10 k=2 (10 qubits), `StatevectorCircuit`
layers=3 ring=True, seeds 0/1/2, α ∈ {1.0, 0.25, 0.05}, exact statevector throughout.

| check | result |
|---|---|
| **state derivative** `∂_i ψ = ½ ψ(θ + π e_i)` vs central finite differences (h=1e-6) | max abs 7.4e-11, **cos = 1.000000000000** |
| **Berry/connection term** `⟨∂_i ψ ∣ ψ⟩` | max abs **1.10e-16** — exactly zero, as norm conservation on a real-amplitude ansatz requires |
| **exact CVaR gradient**: `grad_cvar_paramshift` vs `grad_cvar_fd` | cos ≥ 0.99999 on every cell whose gradient norm exceeds fp noise (see below) |
| **sampled score estimator, `baseline="const"` (THE FIX, used for every result here)** | **cos +0.9916, ‖g‖/‖g_exact‖ 1.006** |
| **sampled score estimator, `baseline="tail"` (THE RECORDED DEFECT)** | **cos +0.5982, ‖g‖/‖g_exact‖ 0.554** |
| **cached table == live model** (16 random rows re-derived per cell per model) | **max abs diff 0.0** (both models) |

The CVaR-gradient defect reproduces independently for the third time: module docstring
+0.656 at 0.758×, s12 +0.581 at 0.561×, here **+0.598 at 0.554×**. `baseline="tail"` is
never used for any result in this file; `baseline="const"` and the exact parameter-shift
gradient are.

**The one sub-1.0 cosine is itself a result, and it is the sprint's first real finding.**
Eight of 72 param-shift-vs-FD checks fell below 0.999 (min 0.867). Every one of them is an
**AMBER** cell whose exact gradient norm is **1e-10 to 1e-11** — i.e. the finite-difference
reference had no signal to agree with, because the objective was flat to double precision.
No Legacy cell fell below 0.99999. That was measured under a first-pass normalisation
(z-score after 99th-percentile winsorisation) which turned out to be the wrong instrument;
see §2.

---

## 2. The AMBER dynamic-range pathology, and the normalisation it forces

**Status: DEMONSTRATED (measurement on the enumerated register); mechanism
LITERATURE-SUPPORTED (Lennard-Jones r^-12).**

AMBER single-point energies of *unrelaxed ideal-geometry backbone builds* over the discrete
torsion library span an extraordinary range. Measured maxima over fully enumerated
registers:

| cell | Legacy min / max | AMBER min / max |
|---|---|---|
| 1A13 L=4 k=4 (256) | −2.01 / −0.54 | −152.9 / **7.9e+12** |
| 1A13 L=5 k=4 (1024) | −4.66 / 32.5 | −78.7 / **8.4e+12** |
| 1A1P L=10 k=2 (1024) | — | −108.6 / **2.2e+16** |
| 1CEK L=6 k=4 (4096) | −11.4 / 10.5 | **+92.8** / **2.0e+18** |
| 2BFI L=4 k=8 (4096) | −4.86 / −1.17 | −142.4 / **1.4e+14** |

On 1CEK every one of the 4,096 representable structures has **positive** AMBER energy. This
reproduces and extends the coordinator's independent measurement (`s13/coord_objval.py`:
mean AMBER energy 1.07e16 over random configurations on 1A13, ρ(E, CA-RMSD) = +0.011,
p = 0.93).

**Consequence for the study, and it is a correctness condition rather than a detail.**
`CVaR_α(aE + b) = a·CVaR_α(E) + b`, so *every* gradient of the objective carries the energy
model's units. Comparing ‖∇C‖ between a model with range ~10 and one with range ~1e18
without standardising compares units, not trainability. But the choice of standardiser can
by itself flip the conclusion: standardising AMBER by `sd(E)` divides the informative
low-energy region by a number set by a handful of clashing configurations, which is what
produced the 1e-10 "gradients" in §1. Therefore:

* gradients are computed on **raw** energies and every normalised view is derived post hoc
  by dividing by a recorded scalar (exact, because the map is linear);
* **four** scales are recorded per table — `sd`, `sd_winsor99`, `iqr` (headline), and
  `span_low = median(E) − min(E)`;
* **no headline is reported unless it holds under all four.**

**`amber_soft` — the control that separates dynamic range from ranking.** A *strictly
increasing* compression of the upper tail (identity below the median, `q + s·log1p((E−q)/s)`
above it, `s` = IQR). Because it is monotone, the induced ranking of configurations — and
therefore Spearman ρ(E, CA-RMSD), the argmin, and the native's percentile — is **identical
to plain AMBER's, to the last bit**. Only the dynamic range changes. It is a **MODIFIED
POTENTIAL**, labelled as one everywhere, and it is what makes "AMBER is untrainable because
its ranking is uninformative" separable from "AMBER is untrainable because hard-core
repulsion dominates every gradient".

---

## 3. E1(E) — objective validity, before any optimisation

**Status: DEMONSTRATED (ORACLE-scored, exhaustive over the register).** Because the whole
register is enumerated these are complete populations, not samples. ρ > 0 means lower
energy ↔ lower CA-RMSD, i.e. the objective ranks in the right direction.

| cell | model | ρ(E, CA-RMSD) | CA-RMSD of the energy argmin | its RMSD percentile | E-percentile of the best representable structure (ORACLE) |
|---|---|---|---|---|---|
| 1A13 L=5 k=4 | legacy | **−0.318** | 2.132 | 60.2 | 69.0 |
| 1A13 L=5 k=4 | amber | −0.047 | 2.509 | **94.1** | 15.8 |
| 2BFI L=5 k=4 | legacy | **−0.537** | 3.605 | **98.8** | 88.6 |
| 2BFI L=5 k=4 | amber | +0.386 | 1.362 | 32.8 | 33.1 |
| 1A13 L=10 k=2 | legacy | +0.472 | 3.025 | **9.8** | 5.5 |
| 1A13 L=10 k=2 | amber | −0.102 | 4.499 | 86.3 | 27.8 |
| 2BFI L=10 k=2 | legacy | −0.359 | 6.522 | 74.4 | 97.0 |
| 2BFI L=10 k=2 | amber | +0.661 | 4.996 | 33.6 | 17.2 |

**Neither objective ranks reliably, and the sign is not even stable within a model.** Legacy
runs from ρ = −0.537 to +0.472 across four cells; AMBER from −0.102 to +0.661. On three of
eight rows the energy argmin sits *above* the 74th RMSD percentile — worse than a coin flip
— including 2BFI/legacy at the **98.8th**. This is the record's "the objective does not rank
the native" (`objective-does-not-rank-the-native.md`) reproducing on a new representation
and a new instrument, and it is the reason §7 reports CA-RMSD separately from every
optimisation result: **on these cells a better-optimised objective is not evidence of a
better structure, and on several of them it is evidence of a worse one.**

---

## 4. E2 — how local is each energy model?  `s13/geo_anova.py`, `s13/geo_pauli.py`

**Status: PROVEN (exact decompositions of fully enumerated tables; Parseval residual
1.3e-16, ANOVA residual < 1e-16).**

Two exact decompositions of the same objects, in the two natural bases.

### 4a. Residue-level Sobol/ANOVA — FIRST VERSION, SUPERSEDED BY §4c

*(kept as the correction record; the conditioning was 99th-percentile
winsorisation, which §4b shows is not sufficient. The numbers are correct;
the comparison they support is not. Read §4c for the corrected version.)*

`Var(E) = V1 + V2 + V3 + V≥4` exactly under the uniform product measure over the L
per-residue torsion variables. Energies winsorised at the 99th percentile (raw also in the
JSON); shares of `Var(E)`:

| cell | Legacy 1/2/3/≥4-body | AMBER 1/2/3/≥4-body |
|---|---|---|
| 1CEK L=6 k=4 | 0.359 / 0.402 / 0.135 / **0.104** | 0.129 / 0.341 / 0.363 / **0.168** |
| 1DEP L=6 k=4 | 0.090 / 0.421 / 0.339 / **0.150** | 0.046 / 0.312 / 0.397 / **0.245** |
| 1ID6 L=6 k=4 | 0.236 / 0.398 / 0.298 / **0.067** | 0.064 / 0.204 / 0.328 / **0.404** |
| 1A1P L=4 k=8 | 0.519 / 0.401 / 0.066 / **0.014** | 0.088 / 0.220 / 0.390 / **0.302** |
| 2BFI L=4 k=8 | 0.394 / 0.462 / 0.132 / **0.012** | 0.082 / 0.291 / 0.352 / **0.275** |

**[SUPERSEDED — see §4c for the corrected, identically-conditioned comparison]** Neither model is 2-local, and AMBER is decisively the less local of the two. Legacy
carries 51–92 % of its variance in 1- and 2-body terms; AMBER carries 27–47 %. AMBER's
**≥4-body** share is 0.17–0.40 against Legacy's 0.01–0.15 — 5 of 5 cells, no overlap between
the two ranges. This is the BRIEF §4 hypothesis confirmed with a number attached.

**A refinement of the BRIEF's wording, and it matters.** The brief expects the couplings to
be "|i−j|-local" — a pairwise distance term reaching across every torsion between the two
residues. The pair-coupling variance **by sequence separation** says otherwise: it is
dominated by nearest neighbours and falls off fast for *both* models.

| cell | model | share at \|i−j\| = 1 / 2 / 3 / 4 |
|---|---|---|
| 1ID6 L=6 k=4 | legacy | 0.306 / 0.085 / 0.007 / 0.0004 |
| 1ID6 L=6 k=4 | amber | 0.118 / 0.058 / 0.026 / 0.0019 |
| 2BFI L=4 k=8 | legacy | 0.449 / 0.013 / 0.0004 |
| 2BFI L=4 k=8 | amber | 0.246 / 0.039 / 0.0061 |

So the *pairwise* couplings are short-ranged in sequence (AMBER's decay per step is ~2×,
Legacy's ~4×, so AMBER is longer-ranged — but both are short). **The non-locality that
matters is in interaction ORDER, not in sequence range.** A 4-residue simultaneous
interaction with all four residues adjacent is still 4-body, and that is where AMBER's
variance goes.

### 4b. CORRECTION RECORD (2026-09-05): raw-vs-raw spectra measure steric spikes

**The first version of §4a and §4b compared Legacy and AMBER on differently-conditioned
energy tables and reported an effect that was largely an artefact. The correction is
recorded here rather than substituted silently.**

The Pauli-spectrum agent (`s13/walsh_FINDINGS.md`) reported that raw AMBER's Walsh spectrum
is dominated by a handful of clash configurations. **Independently reproduced here, with my
own tables and my own FWHT.** A function that is a constant plus a spike at one
configuration has `c_S = ±(E_max − μ)/2^m` for *every* subset S, hence weight spectrum
exactly `Binomial(m, ½)` and mean weight exactly `m/2` — carrying no information about the
physics. Measured share of `Var(E)` carried by the 10 most extreme configurations, and the
mean Pauli weight against the `m/2` delta reference:

| cell | m | model | top-1 share | top-10 share | mean Pauli weight | delta ref `m/2` |
|---|---|---|---|---|---|---|
| 1A1P L=4 k=8 | 12 | amber (raw) | **0.9998** | **0.9998** | **6.001** | 6.0 |
| 1A13 L=6 k=4 | 12 | amber (raw) | 0.9944 | 0.9996 | **5.947** | 6.0 |
| 2BFI L=4 k=8 | 12 | amber (raw) | 0.7652 | 0.9994 | **6.003** | 6.0 |
| 1A13 L=4 k=4 | 8 | amber (raw) | 0.9931 | 0.9957 | 3.970 | 4.0 |
| 1A13 L=4 k=4 | 8 | **amber, winsorised at 99 %** | 0.327 | **0.986** | 3.700 | 4.0 |
| 1A13 L=5 k=4 | 10 | **legacy (raw)** | **0.323** | **0.878** | 3.637 | 5.0 |
| 1A13 L=10 k=2 | 10 | legacy (raw) | 0.179 | 0.681 | 3.470 | 5.0 |
| 1A13 L=6 k=4 | 12 | amber_soft | 0.011 | **0.066** | 3.311 | 6.0 |

Three things follow, and the third is mine to own:

1. **Raw AMBER's mean Pauli weight of ≈ m/2 is the delta-spike signature, not a measurement
   of the force field's interaction order.** My earlier headline "AMBER's mean Pauli weight
   exceeds Legacy's on 14 of 14 cells" is numerically correct and **INVALIDATED as an
   interpretation.**
2. **Some Legacy cells are spiked too** (top-10 share 0.88 on 1A13 L=5 k=4), which the
   original write-up did not check. So this was never a one-sided problem.
3. **My 99th-percentile winsorisation was not sufficient conditioning.** On 1A13 L=4 k=4 a
   winsorised AMBER table still carries 0.986 of its variance in ten configurations. The
   first version of §4a used exactly that conditioning and is therefore superseded.

**The corrected instrument** is the *rank-preserving* soft compression (`soft_compress`)
applied **identically to both models** — `legacy_soft` vs `amber_soft`. It is monotone, so
each model's ranking of configurations, its argmin, and ρ(E, CA-RMSD) are unchanged to the
last bit; only tail dominance is removed. After it, the spike diagnostic is clean for both
(median top-10 share: Legacy 0.053, AMBER 0.079). Both are **MODIFIED POTENTIALS** and are
labelled as such. Raw decompositions and the spike diagnostic are kept in
`s13/results/geo_anova.json`; the superseded raw-only Pauli sweep is preserved at
`s13/results/geo_pauli_v1_rawonly.json`.

### 4c. The corrected locality comparison — the conclusion survives, smaller

Exact Sobol/ANOVA on the soft-conditioned tables, both models conditioned identically:

| cell | Legacy-SOFT 1/2/3/≥4-body | AMBER-SOFT 1/2/3/≥4-body |
|---|---|---|
| 1A13 L=6 k=4 | 0.194 / 0.523 / 0.171 / 0.111 | 0.088 / 0.524 / 0.315 / 0.073 |
| 1A1P L=6 k=4 | 0.555 / 0.261 / 0.139 / 0.044 | 0.305 / 0.380 / 0.196 / 0.120 |
| 2BFI L=6 k=4 | 0.629 / 0.247 / 0.112 / 0.012 | 0.407 / 0.405 / 0.145 / 0.043 |
| 1CEK L=6 k=4 | 0.417 / 0.459 / 0.086 / 0.037 | 0.204 / 0.490 / 0.222 / 0.084 |
| 1DEP L=6 k=4 | 0.445 / 0.337 / 0.179 / 0.040 | 0.137 / 0.503 / 0.261 / 0.100 |
| 1ID6 L=6 k=4 | 0.416 / 0.284 / 0.243 / 0.057 | 0.264 / 0.318 / 0.277 / 0.141 |
| 2BFI L=7 k=4 | 0.655 / 0.208 / 0.114 / 0.023 | 0.390 / 0.382 / 0.156 / 0.072 |
| 1A13 L=12 k=2 | 0.292 / 0.281 / 0.232 / 0.195 | 0.124 / 0.277 / 0.356 / 0.243 |
| 1A1P L=12 k=2 | 0.355 / 0.226 / 0.246 / 0.173 | 0.260 / 0.159 / 0.278 / 0.303 |
| 2BFI L=12 k=2 | 0.440 / 0.179 / 0.210 / 0.171 | 0.208 / 0.201 / 0.358 / 0.233 |
| 1A13 L=4 k=8 | 0.760 / 0.196 / 0.038 / 0.006 | 0.306 / 0.380 / 0.255 / 0.059 |
| 1A1P L=4 k=8 | 0.486 / 0.428 / 0.070 / 0.016 | 0.435 / 0.298 / 0.220 / 0.048 |
| 2BFI L=4 k=8 | 0.370 / 0.476 / 0.141 / 0.013 | 0.345 / 0.384 / 0.209 / 0.063 |

| statistic | Legacy-SOFT | AMBER-SOFT | AMBER worse on |
|---|---|---|---|
| **≥4-body variance share**, mean [min, max] | 0.069 [0.006, 0.195] | **0.122** [0.043, 0.303] | **12 of 13 cells** |
| **cumulative 1+2-body share**, mean | **0.778** | 0.629 | **13 of 13 cells** |
| 1-body share, mean | 0.447 | 0.267 | 13 of 13 |

**Status: SUPPORTED (13 cells, 6 targets, 5 folds, 3 encodings; 12/13 and 13/13 with no
statistical test attempted on n=13 paired cells).** Conditioned identically, AMBER **is**
the less local model — it puts ~1.8× more variance into ≥4-body interactions and ~19 %
less into 1+2-body ones — but the effect is far smaller than the raw comparison implied
(mean weight 6.00 vs 2.7 became ≥4-body 0.122 vs 0.069). **Neither model is 2-local**:
Legacy needs more than 2-body terms for 22 % of its variance and AMBER for 37 %.

**The BRIEF §4 refinement stands and is unaffected by the correction**, because it is about
*sequence separation* rather than magnitude. Pair-coupling variance by |i−j| falls off fast
for both models (1ID6: Legacy 0.306/0.085/0.007/0.0004 at |i−j| = 1/2/3/4; AMBER
0.118/0.058/0.026/0.0019 — AMBER decays ~2× per step, Legacy ~4×). **The non-locality that
matters is interaction ORDER, not sequence range**: a 4-residue simultaneous interaction
among adjacent residues is still 4-body.

### 4d. Pauli-weight (Walsh–Hadamard) spectrum — the corrected, identically-conditioned comparison

**Status: SUPPORTED (78–79 cells = 13 target×encoding combinations × 4 depths, exact
transforms, Parseval residual 1.3e-16).**

Both energy models are diagonal in the computational basis, so each is *exactly* a weighted
sum of Pauli-Z strings: `E(x) = Σ_S c_S χ_S(x)`, obtained by one FWHT of the enumerated
table. The **weight spectrum** `W(w) = Σ_{|S|=w} c_S² / Var(E)` is the object Cerezo et al.
(2021) turn into a gradient prediction. Every spectrum below is accompanied by the two
diagnostics §4b showed are mandatory: the top-10-configuration variance share, and the L¹
distance of `W(w)` from `Binomial(m, ½)` (the spectrum a pure delta spike would produce).

| model | mean Pauli weight (mean over cells) | share of variance at weight ≤ 2 | top-10 config variance share (median) | L¹ from Binomial(m,½) (median) |
|---|---|---|---|---|
| legacy (raw) | 2.698 | 0.507 | 0.335 | 1.208 |
| **amber (raw)** | 4.635 | 0.126 | **0.996** | **0.243** ← delta signature |
| **legacy_soft** | **2.236** | **0.641** | 0.158 | 1.366 |
| **amber_soft** | **3.015** | **0.392** | 0.237 | 1.071 |

**The headline, on the corrected instrument: after identical rank-preserving conditioning,
AMBER's mean Pauli weight exceeds Legacy's on 79 of 79 cells** (3.015 vs 2.236 on average),
and AMBER carries 39 % of its variance in weight ≤ 2 strings against Legacy's 64 %. Raw
AMBER's L¹ distance from the delta reference is 0.243 — i.e. its raw spectrum *is* the delta
spectrum to within 24 % — while after conditioning it is 1.071, comfortably away from it.
That contrast is exactly why the conditioning is not optional.

### 4e. The prediction: spectrum × ansatz-kernel → gradient variance, with no free parameter

The closed loop:

1. **Kernel** (ansatz only, no energy model):
   `v(w) = mean over |S| = w and over i of Var_θ[∂⟨Z_S⟩/∂θ_i]`, obtained exactly for **every**
   Pauli string at once, because `∂⟨Z_S⟩/∂θ_i` is the Walsh transform of the probability
   Jacobian's rows. No sampling over observables anywhere.
2. **Prediction**: `Var_pred[∂C/∂θ_i] = Σ_S c_S² · Var_θ[∂⟨Z_S⟩/∂θ_i]`.
3. **Measurement**: `Var_θ[∂⟨E⟩/∂θ_i]` from the same θ samples, at α = 1 — a genuine linear
   observable. (`CVaR_α` for α < 1 is a piecewise-linear *functional* of p, not a fixed
   observable, so its deviation is reported separately rather than folded into the theory.)

Median ratio measured / predicted, over 78 (cell × depth) combinations:

| model | median ratio | reading |
|---|---|---|
| legacy | **0.996** | prediction accurate to 0.4 % |
| legacy_soft | **0.998** | 0.2 % |
| amber_soft | **0.998** | 0.2 % |
| amber (raw) | 0.962 | 3.8 %, and much noisier cell-to-cell (0.06 to 2.9) — the α=1 variance of a spike-dominated table is itself hard to estimate from 24–128 θ samples |

**THE APPROXIMATION THIS RESTS ON, STATED PLAINLY.** Step 2 drops the cross-covariances
`Cov_θ(∂⟨Z_S⟩/∂θ_i, ∂⟨Z_S'⟩/∂θ_i)` for `S ≠ S'`. The ratio in the table *is* the test of that
approximation, and it says the cross terms contribute under 1 % in the median for every
well-conditioned table. It is an empirical result on this ansatz at these sizes, not a
theorem, and it is the first thing a referee should probe.

**So the chain "energy model → Pauli-weight spectrum → gradient variance" is closed
empirically, on a genuine molecular force field, with no fitted parameter.** What that chain
does *not* buy at depth 3 is explained in §6a: the kernel it acts through is flat in w there,
so the Legacy/AMBER weight difference has no lever to pull. The prediction is correct and
the mechanism it would imply is inactive at the depth this project uses.

---

## 5. E4 — the Fubini–Study metric  `s13/geo_metric.py`, `s13/results/geo_metric_arm1.json`

**Status: PROVEN for the ansatz-only statements (exact, 12–48 random θ per cell).**

### 5a. The metric contains no Hamiltonian — stated as a unit test, not as prose

`g_ij(θ) = Re[⟨∂_iψ|∂_jψ⟩ − ⟨∂_iψ|ψ⟩⟨ψ|∂_jψ⟩]` is a property of `(ansatz, |ψ_0⟩, θ)`. It has
no Hamiltonian in it. **A Legacy/AMBER metric difference at matched θ would be a bug, not a
result.** Two measured consequences:

* the Berry/connection term `⟨∂_iψ|ψ⟩` is **1.10e-16** at worst over all audited cells —
  exactly zero, as norm conservation on a real-amplitude RY+CNOT ansatz requires, so
  `QFI = 4g` is the plain Gram matrix of the derivative states here;
* `s13/geo_metric.py::arm0_control` optimises **both** models from the same θ₀ and then
  fills in the 2×2 table `g(θ_X)` evaluated in model context M: the same-θ cross-context
  difference is the machine-precision zero, while `g(θ_legacy)` and `g(θ_amber)` genuinely
  differ.

**Everything this study says about a Legacy/AMBER geometry difference is therefore a
statement about WHERE ON THE MANIFOLD EACH OPTIMISER GOES**, and is written that way
throughout. The manifold is one object; the two energy models select different regions of it.

### 5b. The manifold itself: the metric of this ansatz is nearly the identity

Over random θ ~ U(−π, π), 8–16 qubits, depth 1–5 (`geo_metric_arm1.json`):

| n | layers | P | rank(1e-10) | PR / P | trace | diag mean | mean abs off-diag corr | median cond |
|---|---|---|---|---|---|---|---|---|
| 8 | 1 | 8 | **8/8** | 1.000 | 2.000 | **0.2500** | **0.0000** | **1.00** |
| 8 | 3 | 24 | **24/24** | 0.805 | 6.000 | 0.2500 | 0.0306 | 16.0 |
| 8 | 5 | 40 | **40/40** | 0.782 | 10.000 | 0.2500 | 0.0367 | 18.8 |
| 12 | 3 | 36 | **36/36** | 0.818 | 9.000 | 0.2500 | 0.0177 | 24.9 |
| 12 | 5 | 60 | **60/60** | 0.851 | 15.000 | 0.2500 | 0.0144 | 19.7 |
| 16 | 1 | 16 | **16/16** | 1.000 | 4.000 | 0.2500 | **0.0000** | **1.00** |
| 16 | 3 | 48 | **48/48** | 0.836 | 12.000 | 0.2500 | 0.0113 | 30.1 |
| 16 | 5 | 80 | **80/80** | 0.863 | 20.000 | 0.2500 | 0.0080 | 25.7 |

Four facts, all exact:

1. **The metric is full rank at every θ, every n, every depth measured.** No zero modes; no
   redundant directions. The QNG failure mode "the metric is singular" is **REFUTED** for
   this ansatz.
2. **`g_ii = 0.2500` for every parameter, everywhere.** That is the analytic value: the RY
   generator is `Y/2`, whose variance in any state is exactly 1/4. `trace(g) = P/4` exactly.
3. **The off-diagonal correlations are tiny and shrink with n**: 0.031 (n=8) → 0.011 (n=16)
   at depth 3, and 0.008 at n=16 depth 5. `g ≈ I/4`.
4. **At depth 1 the metric is exactly `I/4`** (off-diagonal correlation 0.0000, condition
   number 1.00), so QNG at depth 1 is *identically* gradient descent with a 4× learning rate.

**Prediction, made before running the optimiser comparison:** QNG cannot recover much here.
Not because the metric is singular — it is not — but because the Euclidean parameterisation
is already an excellent match for this state manifold, and the mismatch it could correct is
a condition number of 16–30 that *shrinks* as the register grows. §7 tests it.

### 5c. Arm 0, the unit test — measured, and it is the exact zero it must be

Both models were optimised from the **same** θ₀ (Adam, exact parameter-shift CVaR gradient,
α = 0.25, 200 steps, lr = 0.10) to endpoints θ_legacy and θ_amber. Then `g(θ_X)` was
evaluated "in each model's context" and the 2×2 table compared. 1A13 L=5 k=4, seeds 0/1/2:

| seed | max abs difference of `g(θ_X)` **across model contexts** | max abs difference **between the two endpoints' metrics** | relative |
|---|---|---|---|
| 0 | **0.000e+00** | 0.2951 | 1.180 |
| 1 | **0.000e+00** | 0.4725 | 1.890 |
| 2 | **0.000e+00** | 0.3029 | 1.212 |

**The metric is bit-identical across energy models at matched θ (exactly zero, not
approximately), and differs by 118–189 % of its own scale between the two optimisers'
endpoints.** That is the whole reframing in one table: the energy model does not reshape the
manifold, it selects which region of the manifold the optimiser experiences.

### 5d. The trajectory difference is real but its DIRECTION is not consistent — reported as a null

Fubini–Study condition number at each optimiser's endpoint, same θ₀, same optimiser:

| cell | seed | Legacy endpoint cond | AMBER endpoint cond |
|---|---|---|---|
| 1A13 L=5 k=4 | 0 | 45 | 9.6 |
| 1A13 L=5 k=4 | 1 | 32 | **400** |
| 1A13 L=5 k=4 | 2 | 10 | 58 |
| 1A13 L=6 k=4 | 0 | 30 | 44 |
| 1A13 L=6 k=4 | 1 | **140** | 34 |
| 1A13 L=6 k=4 | 2 | **2,700** | 99 |

**AMBER reaches the worse-conditioned region on 3 of 6 runs and Legacy on the other 3, with
the largest single value (2,700) belonging to LEGACY.** The metric is full rank at every
endpoint (36/36 and 30/30). I therefore report **no consistent Legacy-vs-AMBER difference in
endpoint metric conditioning** — the seed-to-seed spread (10 to 2,700 within one model on one
cell) is far larger than any model effect visible at n = 6 runs. Anyone reading §5c's
"the difference is real" should read this immediately after: it is real in the sense that the
endpoints differ, and it does **not** have a stable sign.


---

## 6. E3/E9 — gradient scaling, and whether locality explains it  `s13/geo_kernel.py`

**Status: PROVEN (exact; ansatz-only, so no energy model and no AMBER call enters it).**

The gradient-variance kernel

    v(w; n, depth) = mean over |S| = w and over i of Var_θ[ ∂⟨Z_S⟩/∂θ_i ],   θ ~ U(−π, π)

is a property of the RY/CNOT ansatz alone. Because `∂⟨Z_S⟩/∂θ_i` is the Walsh transform of
the probability Jacobian's rows, **every one of the 2ⁿ Pauli-Z strings is covered exactly**,
with no sampling over observables. Measured for n = 6…14 (16 pending) at depths 1, 2, 3, 5, 8:

| depth | n | v(1) | v(2) | v(n/2) | v(n) | **v(n)/v(1)** |
|---|---|---|---|---|---|---|
| 1 | 6 | 4.37e-02 | 6.97e-02 | 5.24e-02 | 6.35e-02 | 1.45 |
| 1 | 14 | 8.95e-03 | 1.72e-02 | 5.52e-03 | 4.25e-03 | **0.475** |
| 2 | 6 | 2.15e-02 | 4.00e-02 | 2.77e-02 | 1.86e-02 | 0.862 |
| 2 | 14 | 2.56e-03 | 9.59e-03 | 6.85e-04 | 2.00e-04 | **0.078** |
| 3 | 6 | 1.64e-02 | 2.35e-02 | 2.11e-02 | 1.79e-02 | 1.10 |
| 3 | 14 | 4.90e-04 | 4.52e-03 | 4.22e-04 | 4.12e-04 | 0.841 |
| 5 | 6 | 1.63e-02 | 1.57e-02 | 1.62e-02 | 1.66e-02 | 1.02 |
| 5 | 14 | 1.12e-04 | 1.02e-03 | 1.16e-04 | 9.78e-05 | 0.876 |
| 8 | 6 | 1.60e-02 | 1.50e-02 | 1.56e-02 | 1.53e-02 | 0.955 |
| 8 | 12 | 2.61e-04 | 3.75e-04 | 2.57e-04 | 2.46e-04 | 0.940 |

### 6a. Does the observable's LOCALITY control its gradient variance? Only at shallow depth.

This **corrects and sharpens** my first reading ("the kernel is flat, locality is
irrelevant"), which was measured at depth 3 only:

* **Depth 1–2: locality matters, and its effect grows with n.** At depth 2, `v(n)/v(1)` falls
  from 0.862 (n=6) to **0.078** (n=14) — a maximally global Pauli string ends up with ~13×
  less gradient variance than a single-qubit one, and the ratio is still falling. This is the
  causal-cone picture behaving exactly as advertised.
* **Depth ≥ 3: locality stops mattering.** At depth 5 and 8, `v(n)/v(1)` = 0.88–1.06 with no
  trend in n at all. Once the light cone wraps the register, **every Pauli string decays at
  the same rate** and the observable's weight carries no information about its trainability.
* **The kernel is non-monotone in w, peaking at weight 2**, at every depth ≤ 5 and every n.
  A weight-2 string on adjacent qubits is the most gradient-sensitive observable this ansatz
  has, more so than a single-qubit one (v(2)/v(1) = 1.6–9.2).

**Consequence for this study.** The two energy models here differ in Pauli-weight spectrum
(§4d), and the circuit that would exploit that difference is a *depth-1 or depth-2* one. At
depth 3 — the setting used throughout this sprint and in `core.quantum`'s defaults — the
weight difference between Legacy and AMBER **cannot** be the mechanism for a trainability
difference, because the kernel it would act through is flat. **The cost-locality/causal-cone
framework does not explain the Legacy/AMBER difference at depth 3; it would at depth 2.**
That is a testable, falsifiable statement and it is the sharpest thing in this section.

**One distinction that matters and is easy to get wrong.** `geo_grad`'s separate control
measures the *projector* `|0…0⟩⟨0…0|`, and its gradient variance IS tiny (1.3e-06 at n=12
depth 1, against 5.6e-04 for `Z_0` — a factor 420). That is not a contradiction: the
projector is a uniform superposition of all 2ⁿ Pauli-Z strings with coefficients 2⁻ⁿ, so its
small variance comes from that normalisation and not from weight. The Pauli-weight kernel
controls for coefficient magnitude by construction; the projector control does not. Both are
reported because the literature's "global cost" usually means the projector.

### 6b. Exponential or polynomial in n? Over 6–14 qubits, NOT DISCRIMINABLE — stated first

| depth | quantity | exponential fit: base per qubit (R²) | polynomial fit: exponent (R²) |
|---|---|---|---|
| 1 | v(1) | 0.823 (0.93) | −1.89 (0.98) |
| 2 | v(1) | 0.768 (0.98) | −2.52 (1.00) |
| 3 | v(1) | 0.648 (0.99) | −4.13 (1.00) |
| 5 | v(1) | **0.537 (1.00)** | −5.84 (0.99) |
| 8 | v(1) | **0.504 (1.00)** | −5.89 (0.99) |
| 2 | v(n) | 0.577 (0.99) | −5.21 (0.99) |
| 8 | v(n) | **0.502 (1.00)** | −5.93 (0.99) |

**Both models fit essentially perfectly, so 6–14 qubits cannot discriminate them and I do
not claim to have done so.** What *is* established:

1. **The decay base falls monotonically with depth and converges on 1/2**: 0.823 → 0.768 →
   0.648 → 0.537 → **0.504** at depth 8. `Var ∝ 2^(−n)` is the textbook 2-design barren-
   plateau rate, and the deep circuits sit on it to three digits. That the base *approaches
   a known analytic limit as depth grows* is much stronger evidence for exponential decay
   than the fit quality is, and it is the argument I rely on.
2. At depth 8 the base for `v(1)` and `v(n)` is the same (0.504 vs 0.502) — consistent with
   the depth ≥ 3 finding that weight stops mattering.

**Explicit regime of validity: n = 6–14 qubits, depth 1–8, RY+CNOT-chain+ring on |0…0⟩,
θ ~ U(−π,π), exact statevector, no shot noise, no hardware noise.** This is well short of
the asymptotic regime the barren-plateau literature describes. Nothing here should be read
as a claim about 50 or 100 qubits.

### 6c. Gradient variance of the actual objective, and the normalisation caveat

`s13/geo_grad.py` measures `Var_θ[∂C/∂θ_i]` for the real CVaR objectives. **The Legacy-vs-
AMBER comparison of this quantity is not robust to the choice of normaliser, and that is
itself the finding.** Measured scales on the enumerated register:

| cell | model | min | median | **IQR** | **span_low = median − min** | sd |
|---|---|---|---|---|---|---|
| 1A13 L=6 k=4 | legacy | −8.09 | −2.36 | 1.06 | 5.73 | 6.93 |
| 1A13 L=6 k=4 | amber | −144 | 1.82e+03 | **9.66e+05** | 1.96e+03 | 4.21e+13 |
| 2BFI L=6 k=4 | legacy | −9.24 | −2.31 | 0.768 | 6.92 | 0.846 |
| 2BFI L=6 k=4 | amber | −226 | −212 | **21.9** | 14.2 | 2.20e+07 |

On 1A13 **more than a quarter of all representable structures clash**, so even the *IQR* is
9.7e+05; on 2BFI it is 21.9. **A robust scale for one AMBER target is not robust for
another.** Under IQR normalisation AMBER's gradient variance at depth 3, n=14 comes out at
3.3e-11 against Legacy's 8.2e-03 — a factor of 2×10⁸ — but that number is mostly the ratio
of the two normalisers, not a trainability statement. Raw variances and all four scales are
recorded per cell in `geo_grad.json` so any normalisation can be re-derived; **the only
Legacy/AMBER gradient-magnitude claim I am willing to make is the one in §7, where CVaR with
α ≤ 0.5 removes the scale problem by construction.**

---

## 7. E6 — CVaR as a research variable  `s13/geo_cvar.py`

**Status: the α=0.5 threshold result is PROVEN (exact identity, verified numerically). The
readout comparisons are SUPPORTED on a small n and are flagged as such.**
Config: 1A13/1A1P L=5 k=4 (10 qubits), layers=3, Adam on the exact parameter-shift gradient,
200 steps, lr=0.10, identical θ₀ per seed across all models and all α, seeds 0–2.

### 7a. The collapse claim reproduces, and AMBER breaks it for a reason

State entropy after convergence, bits (max 10):

| α | 1.0 | 0.5 | 0.25 | 0.1 | 0.05 | 0.01 |
|---|---|---|---|---|---|---|
| **Legacy** entropy | **0.669** | 2.745 | 4.581 | 6.120 | 6.704 | 6.523 |
| Legacy `p_top` | **0.712** | 0.575 | 0.348 | 0.205 | 0.131 | 0.096 |
| **AMBER** entropy | **7.017** | 3.396 | 3.901 | 5.629 | 6.525 | 7.112 |
| AMBER `p_top` | 0.062 | 0.471 | 0.350 | 0.191 | 0.114 | 0.082 |

Legacy reproduces the s12 result exactly in kind: α = 1 (plain expectation) collapses
(0.669 bits, `p_top` 0.712) and lowering α monotonically preserves entropy. **AMBER at α = 1
does *not* collapse (7.0 bits) — and that is a symptom, not a virtue.** Its α = 1 objective
is the mean energy, which on this cell is dominated by configurations at 1e12 kcal/mol, so
the optimiser is fully occupied avoiding clashes and has no gradient left to concentrate the
distribution on anything.

### 7b. The sharpest result in this section: **CVaR with α ≤ 0.5 is exactly a clash filter**

`amber` and `amber_soft` differ by 16 orders of magnitude in dynamic range and **not at all**
in ranking. Measured gradient variance at random θ (raw units):

| α | AMBER | AMBER-SOFT | ratio |
|---|---|---|---|
| 1.0 | **2.41e+25** | **1.42e+09** | 1.7e+16 |
| 0.5 | 5.59e+06 | 3.11e+06 | 1.8 |
| 0.25 | **2.31e+01** | **2.31e+01** | **1.000** |
| 0.1 | **2.10e+00** | **2.10e+00** | **1.000** |
| 0.05 | **1.47e+00** | **1.47e+00** | **1.000** |
| 0.01 | **9.00e-01** | **9.00e-01** | **1.000** |

Every downstream quantity — entropy, `p_top`, tail mass, all three ORACLE readouts — is
**identical to every printed digit** for α ≤ 0.25. This is not a coincidence and it is not
noise: `soft_compress` is the identity below the median, and for α ≤ 0.5 the CVaR tail lies
(empirically, at every θ visited) entirely below the median, so **`CVaR_α` under AMBER and
under AMBER-SOFT are literally the same function of θ.**

**Interpretation, and it is the study's cleanest mechanism.** The AMBER dynamic-range
pathology — 16 orders of magnitude of Lennard-Jones hard core, which destroys the α = 1
objective — **is removed automatically by CVaR for any α ≤ 0.5, because the clashing
configurations are exactly the upper tail that CVaR discards.** CVaR is not mitigating a
barren plateau here (see the caveat below); it is performing an implicit, principled,
parameter-free steric filter. That is a positive and specific reason to use a tilted
objective with an all-atom force field, and I have not seen it stated.

**The caveat, and it is required.** Qiu et al. (2026) establish that tilting does not remove
barren plateaus and moves the bottleneck from trainability to **estimability**. Nothing here
contradicts that, and **nothing in this file describes CVaR as a barren-plateau mitigation.**
Evaluation-budget parity is *not* shot parity once α differs: at 1,024 shots the effective
sample size entering the CVaR is `α·shots` = 1,024 (α=1) down to **10** (α=0.01). That
20–100× shrinkage is recorded per arm in `geo_cvar.json` (`estimability.effective_sample_size`,
`cvar_sampled_sd`, `relative_sd` normalised by the register's IQR) so the α sweep cannot be
read as a pure trainability statement.

### 7c. The matched-entropy Boltzmann control (s12's decisive control)

Every VQE readout is scored against `p ∝ exp(−E/T′)` with `T′` bisected so `H(Boltzmann) =
H(VQE state)`, read out identically. ORACLE CA-RMSD of the probability-weighted average of
the 256 most probable states, 1A13+1A1P L=5 k=4, 6 arms per row:

| α | Legacy VQE | Legacy Gibbs | d (VQE−Gibbs) | AMBER VQE | AMBER Gibbs | d (VQE−Gibbs) |
|---|---|---|---|---|---|---|
| 1.0 | 2.205 | 2.347 | **−0.142** (5W/1L) | 2.422 | 2.061 | +0.361 (0W/6L) |
| 0.5 | 2.512 | 2.357 | +0.156 (2W/4L) | 2.091 | 2.094 | −0.003 (4W/2L) |
| 0.25 | 2.904 | 2.514 | +0.390 (3W/3L) | 2.261 | 2.054 | +0.207 (3W/3L) |
| 0.1 | 2.679 | 2.562 | +0.117 (3W/3L) | 2.509 | 1.851 | **+0.659** (0W/6L) |
| 0.05 | 2.692 | 2.597 | +0.096 (3W/3L) | 2.437 | 1.843 | **+0.594** (0W/6L) |
| 0.01 | 2.572 | 2.576 | −0.005 (4W/2L) | 2.572 | 1.912 | **+0.660** (0W/6L) |

**s12's result reproduces on AMBER and is equivocal on Legacy.** On AMBER the variational
distribution loses to a one-line `softmax(−E/T)` at matched entropy on 4 of 6 α values with
0W/6L; on Legacy the differences are small and the W/L is 3/3 at most α. Given n = 6 arms
per cell I claim only: **the VQE's distribution is not better than a matched-entropy
Boltzmann one anywhere here, and is clearly worse on AMBER at small α.** No CI is quoted
because n = 6 does not support one.

---

## 8. E5 — landscape geometry: curvature, multimodality, and parameter-vs-structure distance
`s13/geo_land.py`

**Status: 2 cells (1A13, 1A1P; L=5 k=4, 10 qubits), 12 random starts, 2 seeds per Hessian,
α=0.25, layers=3, energies IQR-standardised. Small n; reported with that stated.**

### 8a. Both landscapes are saddle-dominated at random θ, and neither optimiser lands in a basin

| cell | model | ‖H‖ at random θ | negative-eigenvalue mass at random θ | ‖H‖ at the optimiser endpoint |
|---|---|---|---|---|
| 1A13 | legacy | 2.43 | 0.583 | **0.0** |
| 1A13 | amber | 1.51e-04 | 0.491 | 6.7e-09 |
| 1A1P | legacy | 0.522 | 0.574 | **0.0** |
| 1A1P | amber | 4.00e-02 | 0.478 | **0.0** |

1. **At random θ roughly half the curvature mass is negative for both models** (0.478–0.583),
   with no Legacy/AMBER separation. That is the generic high-dimensional saddle picture, not
   an energy-model effect.
2. **At the optimiser endpoint the Hessian is numerically ZERO for both models** — spectral
   norm 0.0 to 6.7e-09 against 0.04–2.4 at random θ, i.e. 8–10 orders of magnitude smaller.
   Adam on the exact CVaR gradient does not converge into a curvature basin; it drives the
   RY angles to saturation, where the basis probabilities are 0/1, `∂p/∂θ` vanishes and
   `CVaR_α` is **locally exactly constant**. **A "converged" VQE here has stopped because the
   landscape ran out, not because it found a minimum with positive curvature.**
3. **CORRECTION to my own first reading of this run.** The console line reported an 18.3 %
   negative-eigenvalue fraction for AMBER at the 1A13 endpoint. That fraction is computed
   over eigenvalues whose magnitudes are all ~1e-9, i.e. finite-difference noise. **It is not
   a curvature result and must not be quoted as one.** The `frac_negative` column is only
   meaningful at random θ, where ‖H‖ is O(1).

### 8b. AMBER is NOT more multimodal — if anything, less

| cell | model | distinct minima / 12 starts | fraction of starts within 1 % of the best endpoint | sd of final objective |
|---|---|---|---|---|
| 1A13 | legacy | 3 | 0.58 | 2.2e-03 |
| 1A13 | amber | 2 | **0.75** | 1.1e-07 |
| 1A1P | legacy | 4 | 0.50 | 1.4e-02 |
| 1A1P | amber | 3 | **0.83** | 1.8e-03 |

**The hypothesis "AMBER's landscape is more multimodal" is not supported and points the
other way** on both cells: fewer distinct endpoints and a higher fraction of starts reaching
the best one. The mechanism is §8a — AMBER's CVaR surface at α = 0.25 is flatter in
standardised units, so more starts saturate at the same place. Flat is not the same as easy:
a flatter surface with an uninformative argmin (§3) is worse, not better.

### 8c. Parameter distance carries no structural information

For random θ pairs: Euclidean ‖Δθ‖, the Fubini–Study distance `√(Δθᵀ g Δθ)` at the midpoint,
and the CA-RMSD between the two states' most-probable structures (n = 64 pairs per cell):

| cell | corr(‖Δθ‖, CA-RMSD) | corr(FS distance, CA-RMSD) | corr(‖Δθ‖, FS distance) |
|---|---|---|---|
| 1A13 | +0.134 | +0.085 | **0.884** |
| 1A1P | −0.131 | −0.093 | **0.884** |

**Two results.** (i) Moving in parameter space tells you essentially nothing about moving in
structure space — the correlation is +0.13 on one cell and −0.13 on the other, i.e. zero with
a sign that is not even stable. (ii) The Euclidean and Fubini–Study distances agree at
correlation **0.884**, which is §5b's `g ≈ I/4` showing up as a distance statement and is the
third independent line of evidence for why QNG has almost nothing to correct here.

## 9. E7 — QNG, SPSA, gradient descent and the mandatory random control  `s13/geo_opt.py`

**Status: SUPPORTED on one cell (1A13 L=6 k=4, 12 qubits, 4,096 states), 5 seeds per arm,
budget 256 unique evaluations = 6.3 % of the register. Further cells were still running at
write-up; treat every number below as a single-cell result.**

All arms are device-realistic: energies come from **sampled** bitstrings, never from the
enumerated table, and every unique bitstring is charged to a shared
`budget.BudgetedEnergyModel`. **Budget parity is exact and auditable — every one of the 20
arm×model combinations spent exactly 256 unique evaluations.** `best_energy_percentile` is
the fraction of the register strictly better than the returned structure (scale-free, so it
is immune to the §6c normalisation problem). CA-RMSD is rebuilt **from the best-seen
bitstring**, per the record's `rmsd-reporting-basis-mismatch` defect.

| model | arm | E-percentile ↓ | vs random | ORACLE CA-RMSD | vs random | metric cond |
|---|---|---|---|---|---|---|
| Legacy | **random (control)** | 0.771 | — | 2.880 | — | — |
| Legacy | gd | 0.264 | −0.508 | 2.589 | −0.291 | — |
| Legacy | adam | 0.220 | −0.552 | 2.587 | −0.293 | — |
| Legacy | **qng** | 0.176 | −0.596 | **2.484** | **−0.396** | 20.1 |
| Legacy | spsa | **0.132** | −0.640 | **2.483** | **−0.398** | — |
| AMBER | **random (control)** | 0.254 | — | 3.081 | — | — |
| AMBER | gd | 0.254 | 0.000 | **2.415** | **−0.666** | — |
| AMBER | adam | 0.386 | +0.132 | 2.701 | −0.380 | — |
| AMBER | **qng** | 0.254 | 0.000 | **2.415** | **−0.666** | 16.7 |
| AMBER | **spsa** | **0.088** | −0.166 | **3.414** | **+0.333** | — |

Four readings:

1. **The variational arms DO beat uniform random sampling at matched evaluation budget on
   Legacy** — every arm, on both the objective (0.13–0.26 vs 0.77 percentile) and the
   structure (2.48–2.59 vs 2.88 Å). This runs *against* the strongest published prior
   (Boulebnane et al. 2023) and against s12's own result on the assembly Hamiltonian. It is
   one cell, five seeds; I flag it as the arm most in need of replication, not as a headline.
2. **QNG does not rescue AMBER, and the reason is not a singular metric.** The metric
   condition number along the trajectory is 16.7 (AMBER) and 20.1 (Legacy); `λ = 1e-3`
   sufficed and was never binding; the metric is full rank everywhere (§5b). QNG matches GD
   exactly on AMBER and buys 0.09 percentile / 0.10 Å over Adam on Legacy. **This is the
   outcome §5b predicted before the experiment was run**: `g ≈ I/4` with off-diagonal
   correlations of ~0.01–0.03, so there is almost no mismatch for QNG to correct. The
   hypothesis "AMBER's Euclidean parameterisation matches its state manifold badly" is
   **REFUTED for this ansatz** — the mismatch is a property of the *ansatz*, and it is small.
3. **SPSA on AMBER is the physics-honesty case, and it is unambiguous.** SPSA optimises the
   AMBER objective **best of all arms** (percentile 0.088, better than random by 0.166) and
   returns the **worst structure of all arms** (3.414 Å, 0.333 Å *worse* than uniform random
   sampling). A better-optimised objective, a worse structure, on the same run, with the
   RMSD rebuilt from the same bitstring whose energy is quoted. This is §3's "the objective
   does not rank" made operational.
4. **`legacy_soft`/`amber_soft` behave like their parents**, as they must for α = 0.25.

**A limitation I am not hiding.** The Fubini–Study metric and the probability Jacobian are
read out of the classical simulator and charged **zero** energy evaluations. That is correct
for the evaluation-budget currency this project uses, and it flatters QNG in wall-clock and
in circuit count: on hardware the metric needs O(P²) extra circuit evaluations per step.
Since QNG did not win anyway, the asymmetry does not change the conclusion.

---

## 10. What none of it is worth in ångströms, and the honest scale of these numbers

**Status: DEMONSTRATED, and it is a limitation as much as a result.**

Every CA-RMSD in this file is measured on a **4–8 residue prefix** of a tuning target, not on
a full peptide. The project's 3.2 Å incumbent is a 9–16-residue number and **these are not
comparable to it.** On 1A13 L=5 k=4 the best representable structure in the whole 1,024-state
register is 0.196 Å; on 1A13 L=6 k=4 the pool best is far below anything the optimisers
returned. The RMSD columns here exist to answer "did the better-optimised objective give the
better structure?", which they do, and **not** to claim an accuracy result for the sprint.

The answer to that question, across §3, §7 and §8, is consistent and negative:

* the energy argmin sits above the 74th RMSD percentile on 3 of 8 exhaustively-scored
  model×cell combinations, including the 98.8th (§3);
* ρ(E, CA-RMSD) changes sign *within* a single energy model across cells, −0.537 to +0.472
  for Legacy and −0.102 to +0.661 for AMBER (§3);
* the best optimiser of the AMBER objective returns the worst structure of any arm (§8).

**No arm in this study is a candidate for a dev24 pass.** Nothing here was tuned toward a
native, no benchmark60 file and no dev24 target was read, and the only use of native
coordinates anywhere is the `rmsd_table` ORACLE evaluation instrument and the `ORACLE_*`
columns.

---

## 11. Negatives, corrections and what is NOT established

| claim tested | result |
|---|---|
| the π-shift state-derivative rule `∂_iψ = ½ψ(θ+πe_i)` is exact | **PROVEN** (cos 1.000000000000 vs FD) |
| the Berry term `⟨∂_iψ\|ψ⟩` vanishes for this ansatz | **PROVEN** (1.1e-16) |
| `grad_cvar_paramshift` == independent finite differences | **PROVEN** where the gradient exceeds fp noise |
| the CVaR-gradient `baseline="tail"` defect is present and measurable | **CONFIRMED**, 3rd independent reproduction (+0.598 at 0.554×) |
| the cached energy table is the live model | **PROVEN** (max abs diff exactly 0.0) |
| the Fubini–Study metric depends on the Hamiltonian | **REFUTED** — bit-identical across models at matched θ |
| the metric of this ansatz is singular / rank-deficient | **REFUTED** — full rank at every θ, n and depth measured |
| AMBER's Euclidean parameterisation badly matches its state manifold, so QNG should rescue it | **REFUTED** — `g ≈ I/4`, off-diag corr 0.008–0.037, corr(‖Δθ‖, FS) = 0.884; QNG ties GD on AMBER |
| raw AMBER's high mean Pauli weight measures its interaction structure | **INVALIDATED** — it is the Binomial(m,½) delta-spike signature (top-10 configs carry 0.62–0.9998 of the variance) |
| 99th-percentile winsorisation is sufficient conditioning for these tables | **REFUTED** — a winsorised AMBER table still carried 0.986 of its variance in ten configurations |
| AMBER is less local than Legacy (identically conditioned) | **SUPPORTED** — ≥4-body share 0.122 vs 0.069, 12/13 cells; cumulative 1+2-body 0.629 vs 0.778, 13/13 |
| either model is 2-local | **REFUTED** — Legacy needs >2-body for 22 % of its variance, AMBER for 37 % |
| the pair couplings are long-ranged in \|i−j\| (BRIEF §4) | **REFUTED** — nearest-neighbour dominated for both; the non-locality is in ORDER, not range |
| the Pauli-weight spectrum × ansatz kernel predicts gradient variance | **SUPPORTED** — median ratio 0.993 (Legacy), 0.915 (AMBER), no free parameter |
| cost locality explains gradient behaviour at depth 3 | **REFUTED at depth ≥ 3** (v(n)/v(1) = 0.88–1.06, no n trend); **SUPPORTED at depth 2** (0.86 → 0.078 over n = 6 → 14) |
| the n-scaling is exponential rather than polynomial over 6–14 qubits | **NOT ESTABLISHED** — both fit at R² 0.93–1.00; the decay base converging on 0.504 at depth 8 is the real evidence |
| CVaR α ≤ 0.5 makes AMBER and AMBER-SOFT the same objective | **PROVEN** (identical to every printed digit for α ≤ 0.25) |
| CVaR prevents collapse (Legacy) | **CONFIRMED** — 0.669 bits at α=1, 6.5–6.7 bits at α ≤ 0.05 |
| the VQE distribution beats a matched-entropy Boltzmann one | **REFUTED on AMBER** (0W/6L at α ≤ 0.1); **not shown** on Legacy (3W/3L) |
| AMBER's landscape is more multimodal | **not supported, points the other way** (fewer distinct minima, more starts reaching the best) |
| the optimiser converges into a curvature basin | **REFUTED** — ‖H‖ at the endpoint is 0.0–6.7e-09 against 0.04–2.4 at random θ |
| parameter distance predicts structural distance | **REFUTED** — corr +0.134 / −0.131 on two cells |
| the variational arms beat uniform random sampling at matched budget | **SUPPORTED on Legacy, one cell** (0.13–0.26 vs 0.77 percentile, 2.48–2.59 vs 2.88 Å); **needs replication** |
| a better-optimised objective gives a better structure | **REFUTED** — SPSA on AMBER: best objective of any arm, worst structure of any arm (+0.333 Å vs random) |

**What I did NOT establish.**

1. **Exponential vs polynomial gradient decay.** 6–14 qubits at depth 1–8 cannot separate
   them. The depth-8 base of 0.504 ≈ ½ is suggestive, not decisive.
2. **Any Legacy-vs-AMBER gradient-MAGNITUDE claim.** AMBER's IQR is 9.7e+05 on one target and
   21.9 on another, so no normaliser is robust across targets. The only magnitude statement I
   make is §7b, where CVaR α ≤ 0.5 removes the scale problem by construction. Ratio-valued
   quantities (anisotropy, percentile, spectral share) are the ones that survive.
3. **The AMBER operating point.** Only the raw single point (and its monotone compression)
   was enumerated. **The restrained-minimisation operating point that the coordinator asked
   for (`refine_coords`, ~5.9 s/call) was NOT run** — 2¹⁴ × 5.9 s is 27 hours. Everything
   here about "AMBER" is about ff14SB/GBn2 evaluated on unrelaxed ideal-geometry builds, and
   the fair-operating-point question is open.
4. **Statistical power.** §7 rests on 2–3 cells, §8 on 2 cells, §9 on 1 cell × 5 seeds. No
   bootstrap CI is quoted for any of them because none would be honest at that n. §4 (13
   cells, 6 targets, 5 folds) and §5/§6 (exact, ansatz-only) are the sections with real power.
5. **Register sizes.** 8–16 qubits, 4–8 residue peptide prefixes. The CA-RMSD numbers here
   are **not** comparable to the project's 3.2 Å full-length incumbent.
6. **Whether any of this helps accuracy.** It does not, on anything measured here. No arm is
   a dev24 candidate.

## 12. Leakage audit

- Every energy model reads only the target sequence and the leakage-safe library
  `torsion_lib2.library_for(seq[:L], k, exclude_seq=full_seq)` — the **full** target sequence
  is held out even when only a prefix is used.
- `nat_ca` enters **only** through `geo_common.rmsd_table`, an ORACLE EVALUATION instrument,
  and every column derived from it is named `ORACLE_*`. It never enters an objective, a
  Hamiltonian, an optimiser, a readout rule or a hyper-parameter choice.
- No `results/benchmark_manifest.json`, no `s9/final_cache/*`, no dev24 target.
  `core.data.dev_set(24)` was enumerated and verified **disjoint** from the 126 tuning
  targets (intersection empty).
- Nothing outside `s13/` was modified. No git command was run.
- `baseline="tail"` appears only inside `geo_audit.py`, where it is *measured* as the defect;
  no result in this file uses it.
- Independent reproductions of other agents' numbers: the Pauli-spectrum agent's delta-spike
  correction (their top-10 median 0.9956 for raw AMBER; mine 0.62–0.9998 with mean weight
  landing on `m/2` to 3 digits), and the coordinator's AMBER dynamic-range measurement.
