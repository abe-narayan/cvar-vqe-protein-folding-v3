# GEOMETRY

The conformational representation, the manifold the pipeline emits onto, the exact gradient that
makes solving possible, and the two opposed geometric pathologies that the incumbent pipeline
happens to contain. This is the document that establishes what is *not* the barrier, which is the
precondition for the rest of the programme's attribution.

---

## 1. THE REPRESENTATION

Backbone conformations are continuous torsions `(φ, ψ) ∈ R^{2n}` under **ideal bond geometry**:
fixed bond lengths and angles, with only the dihedrals free. The Cα trace is built by
`core.project.build_ca_exact`, which is the **bit-exact builder the production projection uses**, so
every structure fitted in this programme lies on the same manifold the shipped pipeline emits onto.
Nothing here is compared across manifolds.

Degrees of freedom: `2n` nominal, **`2n − 4` live** (§3).

---

## 2. THE EXACT ANALYTIC GRADIENT

`core.project.frames` together with `_torsion_grad` supplies an exact **O(n) reverse-mode** gradient
of any Cα-function with respect to `(φ, ψ)`.

**Verification.** Against central differences on real targets: **cosine 1.000000000000000**, with
the only disagreeing entries being those that are numerically zero.

This is what makes the generative frame possible at all. Every objective in the programme — weighted
least squares, maximum likelihood under a binned or continuous density, robust and redescending
losses, torsion priors — is differentiated through this single path, so all fitting is gradient-based
rather than sampled, and every arm is exactly comparable to every other because they differ only in
the coefficient handed to the same chain rule.

---

## 3. FOUR TORSIONS PER CHAIN ARE INERT

`phi[0]`, `psi[0]`, `phi[n−1]` and `psi[n−1]` do not move the Cα trace. Established by two
independent routes:

| route | measurement |
|---|---|
| direct perturbation of ±0.7 rad | Kabsch RMSD moves by **0 to 1.4e-07 Å**, across six targets |
| the analytic gradient at those entries | `|g|` = 0, **2.98e-14**, 0, 0 against a maximum of 3.31e+02 and an interior median of 6.76e+01 — a ratio of **9e-17** |

**Consequences.** The live parameter count is `2n − 4`, not `2n`. In the discrete encoding the live
qubit register is `(n − 2)·log₂k`, not `n·log₂k` — any qubit accounting that ignores this
**overstates the register by four**. And any optimiser that spends budget on those four coordinates
is spending it on nothing; any conditioning analysis that includes them is measuring padding.

---

## 4. THE REPRESENTATION IS NOT THE BARRIER

The single most important geometric result of the programme, because it licenses attributing the
entire remaining deficit to the restraints.

Fitting continuous torsions to the **true** distance matrix, by L-BFGS with the exact gradient, from
native-free starts, selecting among starts by the objective only:

| | all 126 targets |
|---|---|
| mean | **0.611 Å**† |
| median | **0.055 Å** |
| fraction under 2 Å | **0.86** |
| vs incumbent (paired) | **−2.593 [−2.901, −2.295]** |

A median of **0.055 Å** means that on half these targets, ideal-geometry torsions reproduce the
native Cα trace from its own distance matrix to within a rounding error. **The conformational
representation, the ideal-geometry constraint and the search are definitively not the barrier.**

This is a strictly stronger statement than the previous best ceiling — a k = 4 discrete torsion
descent bottoming out at 1.982 Å — because it is continuous, *solved* rather than searched, and does
not depend on a library.

### 4.1 The ceiling this reframes

Two standing project results said perfect distance knowledge caps the answer at **≈1.95–2.0 Å**.
Both were measured with the retrieval library in the loop, using distances to **select, filter and
weight** retrieved windows. That is the ceiling of *selection through a finite library*, not of the
distance channel.

> **Perfect distance knowledge is worth 0.611 Å, not 1.95 Å. The library was carrying the ceiling.**

Every prior estimate of the value of improving distance prediction was computed against the wrong
floor and understated the floor by **3.19×** and the addressable headroom by **2.07×**.

**The counterweight, which must always travel with it.** From *predicted* distances the same
machinery gives **3.644 Å and loses to the incumbent by +0.440 [+0.290, +0.592]**. The headroom is
real; the means to reach it is not demonstrated. The 0.611 Å is an ORACLE ceiling and is never quoted
without the 3.644 Å beside it.

---

## 5. TWO OPPOSED GEOMETRIC PATHOLOGIES

The incumbent pipeline contains two stages whose geometric errors point in **opposite directions**,
which is a striking observation in its own right and may explain part of why it works as well as it
does.

**Coordinate averaging CONTRACTS.** Superposing many structures that disagree and taking the mean
pulls every atom toward the centroid; the more the members disagree, the harder it pulls. Measured
contraction of the production consensus against the true distances: **3.5%**, and **9.5%**
against the predicted ones. *(An earlier draft of this document said 25.8%, reused from a
Sprint 14 record without being checked. It is WITHDRAWN — see `coord_FINDINGS.md` K11. The
direction stands; the magnitude was overstated sevenfold.)*

**The distogram EXPANDS.** It systematically over-predicts distance, and the over-prediction grows
monotonically with sequence separation:

| separation | n pairs | MAE | **bias** | RMSE | z-sd |
|---|---|---|---|---|---|
| 2–3 | 2636 | 0.974 | **+0.113** | 1.463 | 2.186 |
| 4–5 | 2132 | 2.148 | **+0.323** | 2.963 | 2.983 |
| 6–7 | 1628 | 2.818 | **+0.615** | 4.001 | 2.699 |
| 8–10 | 1506 | 3.745 | **+0.931** | 5.216 | 2.701 |
| 11–15 | 647 | 4.666 | **+1.492** | 6.328 | 2.671 |
| **global** | 8549 | 2.386 | **+0.509** | 3.704 | **2.633** |

So a least-squares fit to the distogram produces a systematically **over-extended** structure, while
the aggregation stage produces a systematically **contracted** one. Correcting either in isolation
can therefore make the composite *worse*, and any paper reporting one correction without the other
is reporting half a system.

### 5.1 Why recalibrating the uncertainty cannot help, and why the loss can

The `sd` is **2.633× over-confident** — a genuine reporting defect that must be declared. It is
tempting to treat recalibration as obvious low-hanging fruit. It is not:

The standardised residual's spread is nearly **constant** across separation (2.19, 2.98, 2.70, 2.70,
2.67). So the `sd` captures the error's **shape** correctly and is wrong only by a near-uniform
factor — and **a weighted least-squares argmin is invariant under uniform rescaling of all weights.**
Recalibration alone changes not one emitted structure.

The same argument says what *can* work. A uniform rescaling cannot distinguish a mildly wrong
restraint from a catastrophically wrong one, and the violations are **concentrated**: median
`z = |d_native − d̂|/sd` is **1.336** while the maximum is **7.237**. Only the **loss** can make that
distinction, which is why robust and redescending potentials — the standard tool of NMR restrained
refinement for exactly this failure — are the reformulation the geometry calls for.

---

## 6. THE PROJECTION, AND WHAT IT ACTUALLY COSTS

The pipeline's final geometric step projects the consensus onto the ideal-geometry manifold
(`ramah` penalty at λ = 0.3, multi-start, exact gradient). In the generative cascade this step costs
**+0.170 Å**, which after aggregation is the largest single remaining gap in the pipeline.

**The attribution matters and was initially got wrong.** Sprint 14 recorded this as "the price of
ideal geometry." The standing correction is that it is **the price of re-expanding a contracted
structure** — a different claim with a different fix. Ideal geometry is not expensive; §4 shows a
structure on that manifold can sit 0.055 Å from the native. What is expensive is handing the
projector a consensus that is the wrong size -- and the size error is worth a measured **0.275 A
[-0.397, -0.174]** if corrected with the true scale.

That correction implies a specific, cheap, native-free intervention: rescale the consensus so its
distance scale matches the predicted distances before projecting, using the closed form
`s* = Σ w d d̂ / Σ w d²`. Because the **production** pipeline performs the identical
average-then-project sequence, this is a candidate improvement to the shipped system rather than to
an experimental arm.

---

## 7. THE LOCALITY THEOREM, AND ITS SCOPE

`d_ij` depends on **exactly** the `j − i − 1` residues strictly between `i` and `j` — agreement
1.0000. This is exact for the Cα trace.

**Its scope is narrower than it looks, and the narrowing is the useful part.** All-atom supports are
one residue wider. And both energies in this programme — Legacy and AMBER — are **full-register**:
neither decomposes over a bounded neighbourhood. So the statement "AMBER is less local than Legacy"
is a **category error**; they are equally non-local, and the locality result applies to the geometry,
not to the physics.

---

## 8. RMSD, FROZEN

Independently reimplemented with **Horn's quaternion method** (no SVD; a proper rotation by
construction). Agreement with the production `kabsch_rmsd_batch`: **2.04e-13 Å over 63,000 real
structures**. All 15 edge cases pass, including the mirror test — both implementations return
2.9602 Å where the same code without the determinant fix returns **0.0**.

> **Frozen definition.** Cα only. Every residue including both termini. No trimming. Correspondence
> by residue index. Uniform weights. Proper rotations only. Ångström. Model 1 of the native. Chain
> breaks scored, not rejected.

Full-chain versus `[1:-1]` differs by up to **2.4 Å** on one structure in this set, so the choice is
not cosmetic. Any peptide-RMSD number reported without stating which convention was used is not
comparable — an ambiguity the closest recent preprint to this work leaves unstated.

---

## 9. WHAT THE GEOMETRY DOES NOT EXPLAIN

- It does not explain the accuracy deficit. §4 puts the entire 3.03 Å between the ORACLE and
  predicted fits in the **restraints**, since both arms share optimiser, starts, selection rule and
  manifold.
- It does not explain the selection failure. The objective ranks the native at the 34.7th percentile
  of its own retrieval pool while ordering that pool at ρ = +0.562 — a discrimination problem, not a
  geometric one.
- It does not license `min_sep = 2`, which is undocumented and has **never been swept**, yet is
  inherited by every pair statistic in this document.

---

† **Absolute constants carry a start-draw error bar.** The multi-start initialisation was seeded with
`hash(pdb)`, which Python salts per process, so a run's absolute mean varies with the draw: measured
**sd 0.132 Å** across four independent draws, with an independent workstream's draw returning
**0.569 Å** for this same arm. **Quote this ceiling as ≈ 0.6 Å with a start-draw sd of 0.13 Å, not as
0.611.** Paired *differences* are far more stable (the same quantity measured at 0.381 / 0.416 / 0.486
across draws), which is why the comparative claims in this document survive while its absolute
constants need the bar. Seeding is now stable (`s15/seed.py`), and the frozen protocol is re-run under
it before the benchmark is touched.

‡ **An empirical false-positive floor of ≈ 0.08 Å applies to single-draw paired comparisons.**
Demonstrated directly: a comparison that is **exactly zero by construction** (maximum likelihood
against a constant-width Gaussian *is* least squares) returned **+0.081 [+0.014, +0.169]** on one
start draw and −0.003 on another — a confidence interval excluding zero on a null effect. The
null-calibrated concentration test flags that row, and only that row, as FAIL (p = 0.024). **Effects
at or below 0.08 Å are not resolvable here without replication across draws**, and are marked ‡
wherever they appear.
