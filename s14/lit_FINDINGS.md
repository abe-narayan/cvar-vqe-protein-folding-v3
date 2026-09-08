# Sprint 14 — LITERATURE AND NOVELTY FINDINGS

Agent: `lit`. Date: 2026-09-05. Method: external literature only (WebSearch / WebFetch / local
PDF text extraction of a fetched preprint). **No repository compute run. No repository file
outside `s14/lit_FINDINGS.md` written or modified. Benchmark untouched.**

Tiering per BRIEF §5: **LITERATURE-SUPPORTED** (primary source fetched and read),
**LITERATURE-REPORTED (UNVERIFIED)** (search snippet or secondary only), **HYPOTHESIS** (my
inference, not stated by the source), **REFUTED**.

---

## 0. HEADLINE — THE SEPTEMBER 2026 PREPRINT WAS LOCATED AND READ IN FULL

**The paper exists and I have the complete primary text.** LITERATURE-SUPPORTED.

> Cumbo F, Raubenolt B, Puram V, Katzenmeyer N, Joshi J, Blankenberg D.
> **"Logarithmic-scale variational quantum eigensolver for off-lattice protein structure
> prediction in continuous torsional angle space."**
> arXiv:2609.02113 [quant-ph; physics.chem-ph]. Submitted 2 September 2026, 05:05:46 UTC.
> v1 only. 50 pages, 8 figures. CC BY 4.0. DOI 10.48550/arXiv.2609.02113.
> Affiliation: Computational Life Sciences, Cleveland Clinic Research, Cleveland OH; and
> Dept. Molecular Medicine, Cleveland Clinic Lerner College of Medicine, CWRU.
> Corresponding author: Daniel Blankenberg, blanked2@ccf.org.
> Method name: **QTF (Quantum Torsion Folder)**.
> Code: https://github.com/cumbof/qtf (MIT); support library PHEAT at
> https://github.com/BlankenbergLab/pheat. Data: https://doi.org/10.5281/zenodo.22088098.

I fetched `https://arxiv.org/pdf/2609.02113` (8.8 MB), extracted the text with `pdftotext
-layout`, and read the entire body: abstract, introduction, Materials and Methods, Results,
Discussion, Limitations, Availability, and all 50 references. Every quotation below is from
that text. **No number in this section is invented.** Where the PDF's table layout was
mangled by extraction I say so explicitly and fall back on the prose and abstract, which are
unambiguous.

The one caveat: **arXiv serves no HTML version, and I did not obtain the supplementary
material.** The full mathematical form of the custom energy function is deferred to the
supplement ("For a detailed mathematical description of the full potential and each of these
terms, please see the supplementary material"), as is Figure S2 (side-chain DOF map) and
Figure S1 (ibm_miami panels). So the *term-by-term* algebra of their custom Hamiltonian is
**UNVERIFIED**; its term list and weights quoted in the body are LITERATURE-SUPPORTED.

### The three sentences that matter for Sprint 14

1. **Their headline accuracy numbers are oracle-selected minima over ~7.3 million structures,
   ranked by RMSD to the experimental structure.** The predictive, energy-selected number is
   the *median final model*: **2.90 Å for chignolin** and **5.39 Å for Trp-cage** (custom
   energy, the best of their three backends). Our 3.213 Å across 126 targets is in the same
   regime, on ~13× more targets, with none of the oracle selection.
2. **They independently reproduce our central negative result.** They report Spearman
   energy–RMSD correlations, find them weak, find OpenMM/AMBER total potential energy
   *negatively* correlated with RMSD in both proteins, and state that "model selection remains
   the main bottleneck." Our "neither energy ranks the native" is not a local artefact of our
   pipeline. It is now corroborated by an independent group, on different targets, with three
   different energy functions, over 7.3 M structures.
3. **They have no classical control.** There is no arm anywhere in the paper that optimises
   the torsions directly with the same optimiser, the same energy and the same curriculum,
   without the circuit. The causal contribution of the quantum component is therefore
   unmeasured in this paper. That is precisely the control our brief makes mandatory, and it
   is our clearest scientific differentiator.

---

## 1. DEEP SECTION — arXiv:2609.02113 (QTF), points 1–12

### 1.1 What the paper claims is novel, and what "first" does and does not mean

Two "first" claims are made, both hedged with "to the best of our knowledge":

> "our work provides to the best of our knowledge the first quantum algorithm for protein
> structure prediction in all-atom continuous space" (Abstract/Conclusions)

> "Our team developed a novel method for encoding protein angles into qubit phases,
> presenting, to the best of our knowledge, the first ever QPSP approach modeling protein
> ensembles in continuous space at a full-atom resolution." (Discussion)

The four claimed improvements over the state of the art are enumerated verbatim in the
Introduction: "1) moving away from discretized space and towards a more continuous
representation, 2) modeling all heavy-atoms with full side chains, 3) incorporating
sophisticated molecular mechanics energy functions that capture interactions far beyond the
alpha trace level, 4) logarithmic scaling with respect to the number of protein angles and the
corresponding qubits." LITERATURE-SUPPORTED.

**What "first" specifically means:** the conjunction *(quantum generator) ∧ (continuous, i.e.
non-binned, torsions) ∧ (all heavy atoms with side chains) ∧ (log-qubit register)*. It is a
claim about a *combination*, and the combination does appear to be new.

**What "first" does NOT mean.** Each component is individually prior art, and the paper's own
reference list concedes most of it:

* Not the first quantum torsion-space PSP. Pamidimukkala et al., *JCTC* 20:10223 (2024) —
  their own ref [19] — is "Protein Structure Prediction with High Degrees of Freedom in a
  Gate-Based Quantum Computer", up to 114 qubits on IBM hardware. Torsion/DOF quantum
  approaches predate QTF; QFold (Casares, Campos, Martín-Delgado, 2022) encoded φ/ψ in a few
  bits per angle in 2021–22; Marchand et al. (2018, D-Wave) and Mato et al. (2022, Quantum
  Molecular Unfolding) did torsion-space quantum optimisation over rotatable bonds with real
  force fields. **QTF cites none of QFold, Marchand or Mato.**
* Not the first all-atom quantum treatment of a protein — refs [21]–[27] in their own list
  (Shajan, Kaliakin, Merz, SQD/DMET work) are all-atom quantum *chemistry* on proteins, but
  they score conformers rather than predict folds. The paper distinguishes these correctly.
* Not the first quantum PSP with a real molecular-mechanics force field: QTF's own OpenMM
  ff14SB arm is by construction the same class of object we already run.
* Not a claim of quantum advantage, speedup, or better accuracy than classical methods. **No
  such claim appears anywhere in the paper.** It is a feasibility and resource-scaling claim.
* Not a claim that the structures were produced on quantum hardware. See §1.10.

**Novelty assessment (HYPOTHESIS):** the genuinely new object is the *phase encoding* — using
gauge-fixed relative phases of statevector amplitudes as continuous real-valued variables. I
have not found that specific device elsewhere in the QPSP literature. It is, however, closer
to a *classical reparameterisation with a quantum-shaped nonlinearity* than to a quantum
optimisation (§1.7).

### 1.2 The encoding, written down

**Setup.** An n-qubit register spans a 2ⁿ-dimensional Hilbert space; a normalised state is

    |ψ(θ)⟩ = Σ_{j=0}^{2ⁿ-1} r_j e^{i γ_j} |j⟩ ,   with   n = ⌈log₂ N⌉ ,  N = degrees of freedom

with r_j ≥ 0 the amplitude magnitude and γ_j ∈ (−π, π] the principal phase of the j-th
basis-state amplitude, subject to Σ r_j² = 1. Global phase is unobservable, so 2ⁿ amplitudes
carry at most 2ⁿ−1 independent *relative* phases. QTF therefore uses gauge-fixed relative
phases as the torsion variables. LITERATURE-SUPPORTED, verbatim structure.

**The simulator decoder (the one that produced every headline number).** With
c_j(θ) = ⟨j|ψ(θ)⟩,

    φ̃_j(θ) = wrap_{[−π,π)} ( arg c_j(θ) − arg c_0(θ) )
            = Arg_{[−π,π)} ( c_j(θ) · c_0*(θ) )

The reference c_0 (the all-zeros basis state) is gauge-fixed to phase 0, which makes the
decoded torsion vector — and hence the reconstructed conformation — invariant under global
phase. That leaves **N−1 independent relative phases among N indexed amplitudes.** The zero
reference slot is assigned to the **N-terminal φ**, which is undefined anyway (no preceding
peptide bond) and does not affect reconstruction. The remaining N−1 relative phases "are
mapped to the allowed physical ranges of the corresponding backbone and side-chain torsions."
LITERATURE-SUPPORTED, verbatim.

**Qubit count.** `n = max(2, ⌈log₂ N_torsion⌉)`. Stated concretely in the paper: "For proteins
with up to N ≤ 1,024 degrees of freedom, only n = 10 qubits are needed."

**Instantiated (LITERATURE-SUPPORTED for N and n; the derived d and P are arithmetic from
their stated formulas, HYPOTHESIS-tier only in that the paper does not print them):**

| target | residues | N_torsion (their Table 2) | n = ⌈log₂N⌉ | d = ⌈N/n⌉+2 | P = 2n(d+1) |
|---|---|---|---|---|---|
| chignolin 5AWL | 10 | **44** | **6** (confirmed in Table 4) | 10 | **132** |
| Trp-cage 2JOF | 20 | **88** | **7** | 15 | **224** |

Table 4 independently confirms "Logical QTF qubits 6" for chignolin on both QPUs.

**Angular resolution.** *There is none in the usual sense, and this is the actual novelty.*
The torsions are real numbers read off continuous phases; there is no binning and no
discretisation grid. The paper is explicit: "QTF does not sample torsions in fixed bins, but
instead optimizes continuous torsional degrees of freedom." Their Table 2 discretisation (3,
6, 12 bins per torsion → 9.85e20 to 3.05e47 states for chignolin) is labelled "used only as an
illustrative reference." So the honest answer to "how many qubits for how many torsions at
what angular resolution" is: **6 qubits for 44 torsions at machine (float64) angular
resolution in simulation**, and at *shot-limited* resolution on hardware — where the resolution
is set by 1/n_shots through the CDF, not by the register width. They concede this: "the
accuracy of the extracted molecular geometry is heavily dependent on the statistical resolution
of the measurement vector."

**Periodicity / 2π handling.** Explicit and clean on the simulator side: `wrap_{[−π,π)}` shifts
by integer multiples of 2π to select the representative in [−π, π); `Arg_{[−π,π)}` is the
principal branch. Because the decoder reads a *phase*, periodicity is native — there is no
wraparound discontinuity to patch, which is a genuine and real merit of the representation over
a binned encoding. On the hardware side the CDF map "guarantees that every extracted angle
falls within the physically valid torsion range without any additional clipping or wrapping",
with θ ∈ [−π, π] "with the two endpoints representing the same angular direction."
LITERATURE-SUPPORTED.

**ω handling.** Three modes: excluded and fixed trans; included but affinely mapped into a
170–190° window; or fully free over [−π, π). In windowed mode "the circuit parameter is mapped
into the 170-190 degree interval, avoiding wasted optimizer effort on highly unfavorable
cis-like ω values." LITERATURE-SUPPORTED.

### 1.3 THE CRITICAL STRUCTURAL DEFECT I FOUND: the hardware CDF decoder is monotone

This is my most important technical finding about the paper and it is **derived directly from
their own two printed equations**, so it is LITERATURE-SUPPORTED as arithmetic, with the
consequences marked HYPOTHESIS.

Their hardware decoder is, verbatim:

    C_i = Σ_{j=0}^{i} P_j  ∈ [0,1]        (cumulative measured basis-state probability)
    θ_i = 2π C_i − π

Because P_j ≥ 0, the partial sums C_i are **monotonically non-decreasing in i**. Therefore

    θ_0 ≤ θ_1 ≤ ... ≤ θ_{2ⁿ−1} = π        and        θ_i − θ_{i−1} = 2π P_i ≥ 0

**Consequence 1 — the decoded torsion vector is a sorted staircase.** "The first N values of
the resulting angle vector are used as the molecular degrees of freedom." So on hardware the
torsion vector is *always monotonically non-decreasing in the DOF index ordering*, for every
circuit, every parameter setting, and every noise realisation. The reachable set is a
measure-zero (order-constrained) sliver of the torsion torus — of the order of 1/N! of it if
the DOFs were exchangeable. **The paper never states this and never analyses it.**

**Consequence 2 — the total angular variation is capped at 2π.** Telescoping,
Σ_{i=1}^{2ⁿ−1} (θ_i − θ_{i−1}) = 2π(1 − P_0) ≤ 2π. Across chignolin's 44 used slots in a
64-slot register, the *entire* backbone-plus-side-chain torsion vector has total variation
bounded by 360°, i.e. a mean consecutive step of ≤ 8.2°. Unless the measured probability
distribution is extremely spiky, the decoded angles form a near-linear ramp that is almost
independent of the circuit parameters.

**Consequence 3 (HYPOTHESIS, and it explains their headline asymmetry).** A decoder whose
output is monotone with small total variation is *structurally biased toward torsion profiles
that are nearly constant along the chain* — that is, toward **regular secondary structure**.
Chignolin is a β-hairpin: φ and ψ are close to constant within each strand. Trp-cage is a
helix packed against a polyproline-like segment through a turn: it requires *large jumps*
between three distinct torsion regimes, which is exactly what a bounded-total-variation
monotone decoder cannot produce. **This predicts the paper's central asymmetry — chignolin
works, Trp-cage never crosses 2.0 Å — from the decoder alone, without any appeal to sampling
difficulty or fold complexity.** The paper attributes the gap instead to Trp-cage's "more
heterogeneous fold". Both explanations are consistent with their data; theirs is untested and
mine is testable.

**Why this matters to us directly.** Our own Sprint 12/13 result is that a **zero-information
constant α-helix (φ=−63°, ψ=−42°) emits 4.065 Å and beats uniform random sampling by 0.457 Å
[−0.822, −0.090]**. We have already measured that near-constant torsion profiles are
deceptively strong baselines. A decoder biased toward near-constant profiles will look like it
is working for exactly the reason our brief warns about. The QTF hardware arm has no constant-
torsion control. **The constant-secondary-structure control that our brief mandates is the
single missing experiment that would settle their hardware claim.** HYPOTHESIS.

Note the scope carefully: this applies to the **hardware CDF path only**. The simulator
relative-phase decoder has no monotonicity constraint, and all of 0.623 / 1.199 / 2.501 /
3.512 Å come from the simulator path.

### 1.4 VQE, ansatz, optimiser, decoder — how they actually interact

* **Ansatz:** Qiskit `EfficientSU2` by default (they cite Kandala et al., *Nature* 549:242,
  2017). Each repetition: independent single-qubit Ry and Rz on every qubit, then a **circular
  nearest-neighbour CX entangling layer**; a final Ry/Rz block after the last entangler.
  Repetitions `d = ⌈N_torsion/n⌉ + 2`, chosen "heuristically". Parameter count
  **P = 2n(d+1)**. Their own worked example: "a 6-qubit circuit with five repetition layers
  contains 72 trainable parameters acting on a 64-dimensional Hilbert space." An optional
  **brickwork** ansatz (even-pair then odd-pair CX layers) is offered for hardware.
  LITERATURE-SUPPORTED.
* **Optimisers:** **COBYLA** (derivative-free) in stage 1; **SLSQP** (numerically
  differentiated) in stages 2 and 3. No gradient of the circuit is ever taken analytically —
  no parameter-shift rule, no adjoint, no natural gradient. LITERATURE-SUPPORTED.
  *Note the internal tension:* the paper stresses that NERF is "mathematically smooth and fully
  differentiable... ensuring that spatial gradients... can be effectively backpropagated
  through the geometrical reconstruction and mapped to the quantum circuit parameters", but
  the implementation then uses a derivative-free optimiser and a finite-difference one. The
  differentiability is architectural, not exploited. LITERATURE-SUPPORTED (both statements are
  theirs); the tension is my reading, HYPOTHESIS.
* **Initialisation:** deterministic. RNG seed = SHA-256 of the amino-acid sequence, offset by
  replica index / SLURM array task id. Each replica then runs a **basin-hopping scout**: "By
  default, 50 random parameter sets were evaluated, and the parameter set producing the lowest
  initial cost-function value was selected as the starting point for full optimization."
  LITERATURE-SUPPORTED. **This is a best-of-50 classical pre-selection before the VQE begins.**
* **Shot budget:** *none in the optimisation loop.* All optimisation is exact statevector.
  Hardware shots appear only in the final parameter-transfer jobs: **8,192 shots** per job via
  `SamplerV2`. LITERATURE-SUPPORTED.
* **The loop:** circuit → statevector → gauge-fixed relative phases → torsion vector → NERF →
  all-heavy-atom Cartesian coordinates → classical scalar energy → COBYLA/SLSQP updates θ.
  The energy is **never** measured as a qubit observable. Their own words: "the total scalar
  score is not decomposed into quantum observables and measured as the expectation value of a
  qubit Hamiltonian; it is computed classically from the reconstructed Cartesian structure
  after the circuit readout and minimized by the classical optimizer." LITERATURE-SUPPORTED,
  verbatim. **This is the most important admission in the methods.**

### 1.5 The energy functions — is any of it a real force field?

Three interchangeable backends, 400 replicas each per protein:

1. **"Custom" (their own).** A linear combination, term list quoted from the text:
   E_total = E_clash + E_vdW + E_elec + E_hbond + E_solv/burial + E_disulfide + E_rama +
   E_rotamer (+ E_ω when ω is a DOF, split into ω-trans and ω-window) + E_e2e (end-to-end
   hairpin bias) + E_geom (L-chirality, peptide-plane twist, proline ring closure).
   Atoms carry "coarse physicochemical parameters, including effective partial charges and van
   der Waals radii, **partly inspired by** AMBER-family molecular mechanics force fields."
   Lennard-Jones is **softened**, with the repulsive branch **logarithmically capped** —
   they give the example form `E = 50 + log(E − 49)`. Solvation is a "soft neighbor-counting"
   burial proxy, i.e. a sigmoid-weighted SASA approximation, not GB/PB. H-bonds are a geometric
   distance-and-angle term with virtual amide hydrogens.
   **Verdict: this is NOT a force field. It is a hand-weighted empirical scoring function with
   deliberate structural biases baked in** (Gaussian wells placed on the α-helix and β-sheet
   Ramachandran basins, rotamer wells on gauche+/trans, plus an end-to-end collapse bias).
   LITERATURE-SUPPORTED.
2. **Rosetta** all-atom and centroid energy via PyRosetta (Alford et al. 2017). A real,
   published, benchmarked scoring function — but a knowledge-based one, not a force field.
3. **OpenMM with AMBER ff14SB** (Eastman 2017; Maier 2015). **This is a genuine all-atom force
   field and it is the same object our project runs** (our `core/amber.py`, ff14SB/GBn2 via
   OpenMM). The paper does not state a solvent model for the OpenMM scorer.

Post-processing: heavy-atom models are protonated, parameterised and **energy-minimised in
GROMACS with amber99sb-ildn** (steepest descent). LITERATURE-SUPPORTED.

**The E_e2e term is the one to watch.** It is an "end-to-end hairpin bias" with weight 50.0 in
stage 1 and tapered 8.0 → 1.5 in stage 3, and Table 3 reports an "Exp. E2E" of 5.51 Å for
chignolin. The paper does not state where the E_e2e target distance comes from. **If it is set
from the experimental structure, it is a native-information leak into the objective.** I could
not resolve this from the main text — the term's algebra is in the unavailable supplement.
**UNVERIFIED, and it is the single most important thing a follow-up should check in the
github.com/cumbof/qtf source.** Flagging it because under our BRIEF §1.2 this would be
disqualifying, and because a strong collapse bias toward a known end-to-end distance would by
itself explain much of the chignolin result.

### 1.6 CVaR — is it used?

**No. CVaR is not used, not mentioned, and not cited anywhere in the paper.** I grepped the
complete extracted text for "CVaR", "conditional value at risk", "quantile", and "tail": zero
substantive hits (the only "tail" hits are "high-RMSD tails" in the results prose).
Barkoutsos et al. is not in the 50-reference list. LITERATURE-SUPPORTED (by exhaustive
negative search of the primary text).

**The closest equivalents in QTF, and they are all classical best-of-N with an oracle:**

| QTF mechanism | what it is | selector |
|---|---|---|
| Scout initialisation | best of 50 random parameter vectors | **energy** (legitimate) |
| Snapshot retention | up to 5,000 lowest-energy structures per replica, deduplicated by ≥0.25 raw energy units | **energy** (legitimate) |
| Replica ensemble | 400 replicas per energy function, 1,200 per protein | mixed |
| **"Best retained snapshot"** | argmin over the whole snapshot pool | **RMSD to the native** — ORACLE |
| **"Best final model"** | argmin over 1,157 final models | **RMSD to the native** — ORACLE |
| Hardware best-of-ten | lowest RMSD among 10 repeats of a circuit | **RMSD to the native** — ORACLE |

**Our CVaR is doing something structurally different from all six.** CVaR is a *tail objective
optimised in the loop* — the α-tail of the energy distribution under the ansatz enters the cost
that the optimiser minimises. QTF's snapshot mining is *post-hoc extraction from a completed
trajectory*, and the headline extraction is by RMSD, not by energy. These are not the same
mechanism and are not in competition. **Genuine in-loop CVaR remains unclaimed in QPSP torsion
space.** LITERATURE-SUPPORTED (absence) + HYPOTHESIS (the novelty inference).

### 1.7 Quantum versus classical — ruthless separation

**What the quantum register actually does:** it implements a fixed, parameterised, nonlinear
map

    R^P  ──EfficientSU2──▶  C^{2ⁿ}  ──gauge-fixed relative phase──▶  [−π,π)^{N−1}
    θ (P real parameters)      statevector                            torsion vector

and nothing else. Then:

* The **energy is computed classically**, on classical Cartesian coordinates, by a classical
  scoring function. Stated by the authors (§1.4 quote). There is no qubit Hamiltonian, no
  Pauli decomposition, no expectation value, no measurement of an observable.
* The **coordinates are built classically** by NERF.
* The **optimisation is classical** (COBYLA, SLSQP), on P real parameters, with no quantum
  gradient.
* The **initialisation is a classical best-of-50 random search.**
* The **snapshot mining, ranking, filtering and selection are classical.**
* The **final structures are further minimised classically** in GROMACS/amber99sb-ildn.
* On hardware, the quantum device is used **only as a sampler for a parameter-bound circuit,
  once, with no optimisation** — 8,192 shots, no feedback loop.

**The decisive arithmetic.** P = 2n(⌈N/n⌉ + 3) ≈ **2N + 6n**. So:

* chignolin: **132 classical parameters** to produce **43 torsions** using **6 qubits**;
* Trp-cage: **224 classical parameters** to produce **87 torsions** using **7 qubits**.

**The classical optimiser's search space is ~3× LARGER than the torsion vector it is used to
produce.** The "exponential reduction in qubits" does not reduce the optimisation problem at
all — it *inflates* it, and then hides the inflation in the circuit-parameter vector. This is
the crux of the paper. It is a **reparameterisation**, not a compression: a classical
optimiser searching R^132 through a quantum-shaped nonlinearity onto R^43. Nothing about the
result depends on the map being unitary, entangling, or quantum. LITERATURE-SUPPORTED
(the formulas are theirs); the reparameterisation reading is **HYPOTHESIS**, but it is a very
hard one to argue against.

The authors partly concede the resource point: "By converting physical qubit constraints into
circuit depth constraints..." — width is traded for depth, and depth is the scarcer NISQ
resource. Their measured chignolin circuit is **455 median transpiled depth with 204 CZ gates
on ibm_cleveland**, for a 6-qubit, 10-residue problem.

**The missing experiment.** Nowhere in 50 pages is there a control that runs COBYLA/SLSQP
directly on the 44 torsions, under the same energy, the same three-stage curriculum, the same
scout initialisation, and the same snapshot retention. Without it, the paper cannot attribute
any part of its result to the quantum component — and it does not claim to. **This is exactly
the causal VQE-vs-classical control at matched budget that BRIEF §1.4 requires of us, and it is
our single strongest differentiator.** LITERATURE-SUPPORTED (absence), HYPOTHESIS (import).

### 1.8 Every reported RMSD, and whether it is a retained snapshot or a final model

This is the table the sprint coordinator needs. **The distinction is load-bearing and the
paper is, to its credit, scrupulous about labelling it — but the abstract is not.**

| value | target | backend | class | selector | tier |
|---|---|---|---|---|---|
| **0.623 Å** | chignolin | Rosetta snapshot pool | **best retained snapshot**, global argmin over ~2.68 M structures | **RMSD (oracle)** | oracle |
| **1.199 Å** | chignolin | "any backend" | **best final model**, argmin over 1,157 finals | **RMSD (oracle)** | oracle |
| **2.90 Å** | chignolin | custom | **median final model** | energy (in-loop) | **predictive** |
| 4.20 Å | chignolin | Rosetta | median final model | energy | predictive |
| 3.84 Å | chignolin | OpenMM | median final model | energy | predictive |
| ~1.2–1.4 Å | chignolin | all three | **median of per-replica best snapshots** | RMSD (oracle) | oracle |
| **2.501 Å** | Trp-cage | — | best retained snapshot over ~4.65 M | **RMSD (oracle)** | oracle |
| **3.512 Å** | Trp-cage | — | best final model | **RMSD (oracle)** | oracle |
| **5.39 Å** | Trp-cage | custom | **median final model** | energy | **predictive** |
| 6.59 / 6.21 Å | Trp-cage | Rosetta / OpenMM | median final model | energy | predictive |
| >3 Å | Trp-cage | all | median per-replica best snapshot | RMSD (oracle) | oracle |
| **1.758 Å** | chignolin | OpenMM start | best of 300 ibm_cleveland jobs | **RMSD (oracle)** | oracle, hardware |
| **1.782 Å** | chignolin | OpenMM start | best of 300 ibm_miami jobs | **RMSD (oracle)** | oracle, hardware |
| 2.233 / 2.316 / 2.234 Å | chignolin | custom / Rosetta / OpenMM | **median** of 100 ibm_cleveland jobs each | none | **predictive-ish** |
| 2.301 / 2.475 / 2.293 Å | chignolin | custom / Rosetta / OpenMM | median of 100 ibm_miami jobs each | none | predictive-ish |
| 1.784 Å | chignolin | — | median of the 30 hand-picked warm-start simulation models | RMSD (they were chosen as the top-10 per backend) | oracle |

**Table 3 in full, recovered and cross-checked.** My first extraction (`pdftotext -layout`)
interleaved the rows; re-extracting pages 19–20 in raw reading order (`pdftotext -f 19 -l 20`)
recovers the table cleanly, column by column. All values in Å.

| | 5AWL custom | 5AWL Rosetta | 5AWL OpenMM | 2JOF custom | 2JOF Rosetta | 2JOF OpenMM |
|---|---|---|---|---|---|---|
| **Final RMSD med.** | **2.90** | 4.20 | 3.84 | **5.39** | 6.59 | 6.21 |
| Final RMSD best | 1.38 | **1.20** | 1.62 | 3.51 | 3.60 | 3.53 |
| Snapshot RMSD med. | 1.22 | 1.28 | 1.36 | 3.24 | 3.32 | 3.49 |
| Snapshot RMSD best | 0.64 | **0.62** | 0.68 | **2.50** | 2.62 | 2.77 |
| Exp. E2E | 5.51 | 5.51 | 5.51 | 11.25 | 11.25 | 11.25 |
| Final E2E med. | 7.84 | **16.78** | 14.86 | 14.11 | **24.66** | 21.61 |
| Best-snapshot E2E med. | 7.36 | 7.88 | 7.85 | 10.09 | 9.97 | 9.89 |
| Exp. Rg | 4.85 | 4.85 | 4.85 | 6.91 | 6.91 | 6.91 |
| Final Rg med. | 5.16 | 6.57 | 6.22 | 7.91 | 10.05 | 9.15 |
| Best-snapshot Rg med. | 5.02 | 5.04 | 5.03 | 7.46 | 7.48 | 7.50 |

**Four independent cross-checks confirm the recovery**, since the table's two-decimal cells must
match the abstract's three-decimal numbers: Rosetta 5AWL snapshot best **0.62 = 0.623**;
Rosetta 5AWL final best **1.20 = 1.199**; custom 2JOF snapshot best **2.50 = 2.501**; custom
2JOF final best **3.51 = 3.512**. The table is therefore **LITERATURE-SUPPORTED without
caveat.** Note also that this locates the two headline chignolin numbers precisely: both the
0.623 Å snapshot and the 1.199 Å final model came from the **Rosetta** backend, not the custom
one the paper otherwise favours.

**And Table 3 contains a result the paper does not draw out.** Look at the E2E row. Rosetta's
chignolin final models have a **median end-to-end distance of 16.78 Å against an experimental
5.51 Å** — three times over-extended — and OpenMM's is 14.86 Å, while the custom energy's is
7.84 Å. The same pattern holds on Trp-cage (24.66 and 21.61 vs 14.11, experimental 11.25). The
paper reads this as the custom function having better "biological fidelity". **The simpler
reading is that the custom energy contains an explicit end-to-end collapse bias term (E_e2e,
weight 8.0 tapering to 1.5) and the other two do not.** Its advantage on the one structural
metric its own bias term directly controls is close to tautological. This sharpens the
UNVERIFIED question in §1.5 considerably: if E_e2e's target distance is native-derived, the
custom energy's entire compactness advantage — and much of its 2.90 Å vs 4.20 Å final-model
edge — is a leak. HYPOTHESIS, and checkable in their source.

**Native-like hit rates (LITERATURE-SUPPORTED, prose, unambiguous), chignolin, fraction of
replicas producing at least one structure below 2.0 Å:**

| | snapshot mining | **final model selection** |
|---|---|---|
| custom | 100 % | **6 %** |
| OpenMM | 98 % | **4 %** |
| Rosetta | 99.3 % | **0.3 %** |

**Read that table again.** Their own energy functions convert a ~99 % sampling success rate
into a 0.3–6 % selection success rate. The 0.3 % for Rosetta — the most mature, most
benchmarked scoring function of the three — is a *20-fold worse than chance* outcome relative
to what was in the pool. This is the strongest quantitative statement of the energy-ranking
problem I have seen anywhere in the QPSP literature, and it is stronger than our own numbers.

**Trp-cage:** no final model and no retained snapshot crossed 2.0 Å, across 4.65 M structures.
Median snapshot-mining improvement: 2.18 Å (custom), 2.66 Å (OpenMM), 3.26 Å (Rosetta).

**One protocol caveat that changes comparability.** The paper states, in a hardware aside:
"the reported RMSD **excludes the terminal residues** and evaluates the structural core." I
found **no definition of the RMSD protocol in the Materials and Methods section at all** — no
statement of atom selection (Cα vs all-heavy), no alignment method, no residue range. The
abstract says "Cα RMSD". So chignolin's numbers are, apparently, Cα RMSD over roughly 8 of 10
residues. **Our `I.ca_rmsd` is full-chain Kabsch Cα RMSD over all residues.** Excluding termini
on a 10-mer is a materially easier metric; on our targets our own Sprint 12 work found
terminal dropout is 0.40–0.50 Å cheaper than uniform. **Do not compare their numbers to ours
without this adjustment.** LITERATURE-SUPPORTED (the exclusion statement is verbatim);
the magnitude of the effect is HYPOTHESIS.

### 1.9 Chignolin and Trp-cage separately, with protocol

**Shared protocol.** 1,200 replicas per protein, 400 per energy function, one SLURM array task
each. Seed = SHA-256(sequence) ⊕ task id. Scout: 50 random θ, keep lowest energy. Then the
curriculum: **Stage 1** global hydrophobic collapse, COBYLA, maximal bias weights
(w = 50.0, λ = 15.0); **Stage 2** local refinement, SLSQP; **Stage 3** natural relaxation
(custom energy only — Rosetta and OpenMM run two stages because they lack the bias terms),
burial coefficient 15.0 → 2.5, end-to-end strength 8.0 → 1.5, SLSQP. Up to 5,000 low-energy
snapshots retained per replica from the final stage, deduplicated at ≥0.25 raw energy units.
Optional GROMACS amber99sb-ildn steepest-descent minimisation. LITERATURE-SUPPORTED.

**Chignolin (PDB 5AWL, 10 residues, 166 atoms, compact antiparallel β-hairpin, Tyr–Trp
packing).** 1,157 completed final models; 2,681,469 retained snapshots; 2,682,626 total.
44 torsional DOF → 6 qubits. Best snapshot 0.623 Å (Rosetta pool, oracle); best final 1.199 Å
(oracle); median final 2.90 Å (custom). Hit rates as above. Median best-snapshot ~1.2–1.4 Å.
Best-snapshot Rg close to experimental (Exp. Rg 4.85 Å); best-snapshot E2E "somewhat larger"
than the experimental 5.51 Å.

**Trp-cage (PDB 2JOF, 20 residues, 284 atoms, N-terminal α-helix packed against a
polyproline-like C-terminal segment around a buried Trp).** 1,176 completed final models;
4,650,631 retained snapshots. 88 torsional DOF → 7 qubits. **Nothing below 2.0 Å anywhere.**
Best snapshot 2.501 Å, best final 3.512 Å, median final 5.39 Å (custom). The paper attributes
the difficulty to fold heterogeneity; see §1.3 for a competing decoder-level explanation of
the chignolin/Trp-cage gap — though note that gap appears in the *simulator* arm too, where
the monotonicity argument does not apply, so fold complexity is likely also real.

**Failure mode disclosed:** ~67 OpenMM replicas failed at stage 2 with "numerical instability
in the optimizer" / "numerical parameter-binding instabilities."

### 1.10 IBM hardware — what it actually demonstrates

**Devices.** `ibm_cleveland`, a **156-qubit IBM Heron R2**, heavy-hex. `ibm_miami`, a
**120-qubit IBM Nighthawk R1**, square topology. Native gates RZ, SX, X; CZ entangler.
**300 jobs on each.** LITERATURE-SUPPORTED.

**Resources used (Table 4, LITERATURE-SUPPORTED):** **6 logical qubits** on both. On
ibm_cleveland, transpiled depth 443–506 (median 455), total gates 975–1067 (median 1048), CZ
168–208 (median 204), 2Q depth 143–160. On ibm_miami, depth exactly 206 on all 300 jobs, 449–450
gates, exactly 60 CZ, 2Q depth 60. 8,192 shots per job. Error suppression: 32 twirled
randomisations, XY4 dynamical decoupling, gate and measurement twirling. ibm_cleveland job
≈5 s; ibm_miami ≈35 s (a 7× runtime penalty at *lower* depth, unexplained).

**What it does NOT demonstrate — and the paper says so plainly.** Four things:

1. **No optimisation ran on hardware.** "Each job comprised a single execution and sampling of
   the parameter-bound circuit using the SamplerV2 primitive with 8,192 shots, **rather than
   continued iterative optimization**." The QPU is a one-shot sampler.
2. **It is a warm start from already-good classical answers.** The 30 circuits came from "the
   top 10 replica final models... for each energy function" — simulation models already at
   **1.199 to 2.270 Å, median 1.784 Å**. The hardware was handed a near-native answer.
3. **A different decoder was used, so hardware never reproduced the simulator structures.**
   Their words: "the calculations from the hardware experiments constitute a **parameter-
   transfer experiment** using a measurement-compatible probability decoder, **rather than
   hardware recovery of the simulator-derived torsions**." And: "the empirical CDF mapping does
   not reconstruct or approximate the relative phases used by the statevector decoder."
4. **They state the confound is unresolvable in their design.** "Because a matched noiseless
   CDF-decoded baseline was not included, these contributions cannot be separated in the
   present analysis." So the difference between simulation and hardware structures mixes the
   decoder change with all hardware noise, and neither can be attributed. **This is the same
   class of missing control as the missing classical arm.**

**What it does demonstrate, honestly stated:** running a 6-qubit EfficientSU2 circuit
parameter-bound to good classical parameters, and decoding its measured probabilities through
a monotone CDF map, yields a structure that is *usually worse than the classical input but
retains recognisable β-hairpin character*. Below 2.0 Å in **88/300** ibm_cleveland jobs
(29/100 custom, 28/100 Rosetta, 31/100 OpenMM) and **45/300** ibm_miami (21/100, 12/100,
12/100). Best 1.758 Å (cleveland) / 1.782 Å (miami), both from OpenMM starts, both oracle-
selected minima.

**The most telling number in the hardware section:** **0/10 custom-energy starting models
improved on hardware**, on either device (5/10 Rosetta and 3/10 OpenMM improved on cleveland,
3/10 and 1/10 on miami — and the paper explains this correctly as regression from worse
starts). The distribution of hardware-minus-simulation RMSD for the custom energy function is
"almost entirely positive." **Hardware execution monotonically degraded the best inputs.**
LITERATURE-SUPPORTED.

Also disclosed: on ibm_miami the N-to-C terminal Cα distance ranged 4.482–13.285 Å against an
experimental 5.509 Å — i.e. the chain ends wander by up to 2.4× the native separation while the
"core" RMSD stays moderate. That is the terminal-exclusion caveat of §1.8 doing visible work.

### 1.11 Limitations — stated and unstated

**Stated by the authors (LITERATURE-SUPPORTED, and creditably candid):**

* The statevector is not observable on NISQ hardware; the phase encoding's "core limitation
  lies precisely in its internal representation."
* Shot cost: "as the number of degrees of freedom increases, the state space grows
  exponentially. This means the shot count must also scale up dramatically to accurately sample
  the increasingly diluted amplitude probabilities... this intensive shot requirement and the
  resulting noisy readouts effectively **trade the problem of qubit quantity for a problem of
  measurement complexity**." That is a self-refutation of the resource claim on NISQ hardware,
  and they say it themselves.
* Direct phase recovery would need interference measurements or tomography, "scales
  exponentially and is prohibitive for large registers."
* Energy-function generalisability and landscape imbalance; "no universal force field exists
  that reliably ranks structural accuracy across all contexts."
* Model selection is "the main bottleneck."
* Physical artefacts: low-RMSD snapshots with "physically suspect GROMACS energies... such as
  inferred terminal bond formation occurring when the N- and C-termini were placed in
  abnormally close covalent proximity"; they call for filtering on terminal separation,
  non-local contacts and ring penetration.
* Backend-dependent variability; transpiled gate count does not predict accuracy.
* ~67 OpenMM replicas crashed.

**Unstated, and material (HYPOTHESIS unless noted):**

* **No classical control at matched budget.** LITERATURE-SUPPORTED as an absence. Fatal to any
  causal quantum attribution.
* **The CDF decoder's monotonicity** (§1.3). Unanalysed.
* **P > N.** The parameter inflation is never remarked on; "logarithmic scaling" is stated
  about qubits only, and the reader is left to infer the classical cost.
* **n = 2 targets.** Two proteins, both α/β miniprotein benchmarks that essentially every PSP
  method is tuned against, both with abundant published structures. No held-out set, no fold
  diversity, no statistical comparison across targets, no confidence intervals anywhere in the
  paper. **There is not a single CI or p-value in 50 pages.**
* **The custom energy function's weights are hand-set** (50.0, 15.0, 8.0, 2.5, 1.5) and it is
  the backend that wins. No ablation, no cross-validation of the weights, no statement of how
  they were chosen. Chignolin is a hairpin and the custom energy contains an explicit
  "end-to-end **hairpin** bias." **UNVERIFIED whether the E_e2e target is native-derived**
  (§1.5); if it is, the custom-energy chignolin result is leaked.
* **Ramachandran Gaussian wells on α and β basins** are a strong structural prior deliberately
  built into the objective. Our own measurement says that generic Ramachandran accounts for
  **88 %** of what a sequence-conditioned torsion library buys — so a large fraction of QTF's
  apparent skill may be that prior, not the search. HYPOTHESIS with our own supporting number.
* **No RMSD protocol in the Methods.** Terminal exclusion is revealed only in passing.
* **No random-sampling control**, no constant-secondary-structure control, no
  scrambled-sequence control.
* **The scout (best-of-50) is never ablated**, so we cannot tell how much of the final quality
  is the VQE and how much is the classical random search that precedes it.

### 1.12 THE ENERGY-VERSUS-STRUCTURE MISMATCH — they measure it, and they find what we found

This is point 12 of the remit and the answer is emphatic. **They do report it. They report it
three separate ways. It reproduces our result.**

**(a) Spearman rank correlations, energy terms vs RMSD, final models (their Figure 5).**
Verbatim: "we computed the Spearman rank correlations between structural accuracy (RMSD) and
the calculated energies of the final selected models." Findings, verbatim:

* Custom: "a more balanced mixture of positive and negative correlations, indicating that no
  single class of terms dominates the relationship with RMSD... the current weighting does not
  yet provide a clean monotonic ranking of the best structures."
* Rosetta: "a more penalty-like pattern... many of the strongest centroid and full-atom terms
  are positively correlated with RMSD" (i.e. the correct direction) — "however... Rosetta can
  still produce broad low-energy regions containing both near-native and higher-RMSD
  structures, so these correlations are useful but not sufficient."
* **OpenMM (AMBER ff14SB): "OpenMM shows the opposite tendency: total OpenMM potential energy
  is negatively correlated with RMSD in both proteins. This means that lower OpenMM energies
  are not consistently associated with lower RMSD final models in this dataset, and may even
  favor structures that are physically favorable under the force field but topologically
  farther from the experimental fold."**

**That last sentence is our Sprint 13 finding, written by a different group.** Our numbers:
raw AMBER rank correlation with RMSD **−0.088**; native at the **32nd–40th percentile**;
Legacy's certified global optimum over 262,144 enumerated configurations is **+0.139 Å worse
than random sampling**. Theirs: AMBER total potential energy negatively correlated with RMSD on
both chignolin and Trp-cage, over millions of structures. Two independent instruments, two
independent codebases, two disjoint target sets, same sign, same conclusion.
LITERATURE-SUPPORTED.

**(b) Rescoring on a shared scale.** "When sampled structures were rescored on a shared
minimized-potential scale (GROMACS using amber99sb-ildn), low-RMSD conformations were often
present but were **not consistently assigned the lowest energies**, and snapshot energy-RMSD
Spearman correlations were **weak**." So a second all-atom AMBER variant, applied as a
post-hoc rescorer to 7.3 M structures, also fails to rank. LITERATURE-SUPPORTED.

**(c) Funnel plots (their Figure 4).** The custom energy "heavily concentrates native-like
5AWL snapshots into a distinct low-RMSD region, **yet it still fails to cleanly rank these
highest-accuracy structures as the lowest-energy structures among the sampled snapshots.**"
Rosetta and OpenMM "exhibit broad low-energy bands that encompass both near-native and much
higher-RMSD conformations, presenting a relatively flat selection gradient at the bottom of
the optimization basin." LITERATURE-SUPPORTED.

**Do they report the correlation between energy and RMSD?** Yes — qualitatively, by sign and
by figure. **What they do not report is a single numeric ρ.** Not one Spearman coefficient
appears as a number in the text; Figure 5 is described but its values are never tabulated, and
the snapshot correlations are called "weak" without a number. **That absence is itself a
finding, and it is where our project is quantitatively ahead of the published state of the
art:** we have ρ = +0.043 for Legacy inside the low-energy decile, −0.088 for raw AMBER, a
native percentile, a certified global optimum on 9 fully enumerated targets (262,144
configurations each), and a demonstration that monotone log compression leaves ρ unchanged
while collapsing AMBER's range from 16.1 to 1.8 decades — i.e. that the damage is in the
ordering, not the scale. **QTF establishes the phenomenon; we have measured it.**

**Is their energy's optimum structurally good?** Their evidence says no, and it agrees with
ours. The 0.3 % Rosetta final-model hit rate against a 99.3 % pool hit rate is the cleanest
statement of it. Their explanation is the standard field defence — "a useful scoring function
can still guide sampling toward native-like regions... while failing to rank the most accurate
structure as the absolute minimum" — which is a weaker claim than the one our enumeration
licenses. **We have shown the certified global optimum is worse than random; they have shown
the argmin is not the best of the pool. Ours is the stronger result and it is unpublished.**

---

## 2. LITERATURE LEDGER

Verification key: **[P]** = I fetched and read a primary source in this sprint;
**[P13]** = primary source fetched and read by the `lit` agent in Sprint 13, carried forward
from `s13/lit_FINDINGS.md`; **[S]** = secondary — search snippet, abstract-only fetch, or
citation appearing in another paper's reference list. Anything marked **[S]** must not be
quoted as fact without re-checking.

### 2A. Quantum protein structure prediction

| # | Citation | Date | Method | Representation | Objective | Quantum component | Classical component | Dataset | RMSD / result | Limitations | Relation to us | Ver |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **Cumbo, Raubenolt, Puram, Katzenmeyer, Joshi, Blankenberg. arXiv:2609.02113 (QTF)** | Sep 2026 | Log-scale VQE, "Quantum Torsion Folder" | **Continuous torsions** from gauge-fixed statevector relative phases; all heavy atoms + side chains | Custom empirical E; Rosetta REF15; OpenMM ff14SB | EfficientSU2, n=⌈log₂N⌉ qubits (6 chignolin / 7 Trp-cage); statevector phases; hardware = CDF-of-probabilities decoder | NERF build, all energy evaluation, COBYLA+SLSQP, best-of-50 scout, snapshot mining, GROMACS minimisation | chignolin 5AWL, Trp-cage 2JOF (**n=2**) | oracle best snapshot 0.623 / 2.501 Å; oracle best final 1.199 / 3.512 Å; **median final 2.90 / 5.39 Å**; hardware best 1.758 / 1.782 Å | No classical control; no CI/p-value anywhere; n=2; RMSD excludes termini and is undefined in Methods; CDF decoder is monotone; P≈2N+6n > N | **The direct comparator.** Corroborates our energy-ranking negative; lacks every control we have | **[P]** |
| 2 | Perdomo-Ortiz, Dickson, Drew-Brook, Rose, Aspuru-Guzik. *Sci Rep* 2:571 | 2012 | Quantum annealing | 2D/3D lattice HP | HP contact energy | D-Wave annealer | Embedding, energy construction | 6-mer lattice model | lattice ground state found | Toy lattice; no 3D structure | Origin of the field; the lattice tradition QTF departs from | **[S]** |
| 3 | Babbush, Perdomo-Ortiz, O'Gorman, Macready, Aspuru-Guzik. arXiv:1211.3422 | 2012 | AQC, energy-function construction | Lattice heteropolymer | Constraint-satisfaction | Adiabatic QC | QUBO construction | Lattice models | — | Constraint terms dominate the Hamiltonian | The constraint-overhead problem QTF cites as motivation | **[S]** |
| 4 | Fingerhuth, Babej, Ing. arXiv:1810.13411 | 2018 | QAOA with hard+soft constraints | Lattice turns | HP/MJ | Quantum alternating operator ansatz | Constraint encoding | Lattice peptides | — | Lattice only | Prior QAOA-on-lattice | **[S]** |
| 5 | **Robert, Barkoutsos, Woerner, Tavernelli. *npj Quantum Inf* 7:38** | 2021 | Resource-efficient VQE (CVaR) | **Turn encoding**, tetrahedral lattice, 2 qubits per turn | MJ-style contact + interaction qubits | VQE, ~2(N−3)+interaction qubits, on 20-qubit IBM | Hamiltonian construction, decoding | 7–10-mer peptides | native lattice folds recovered | Lattice; coarse-grained; no all-atom | The canonical turn encoding; **uses CVaR**, so CVaR-on-lattice-QPSP is taken | **[P13]** |
| 6 | **Casares, Campos, Martín-Delgado. QFold. arXiv:2101.10279; *Quantum Sci Technol* 7:025013** | 2021/22 | Quantum Metropolis via quantum walks | **Torsion angles, off-lattice** — "does not require a lattice model simplification... parameterization in terms of torsion angles" | Psi4 real energies | Quantum walk Metropolis; minimal realisation on IBMQ Casablanca | Minifold deep-learning dihedral initialiser; energy enumeration | small peptides | polynomial quantum advantage vs classical Metropolis annealing schedules | Angles binned into few bits; tiny systems | **Kills the unqualified "first off-lattice quantum PSP" reading of QTF.** QTF does not cite it | **[S]** + **[P13]** |
| 7 | Chandarana, Hegade, Montalban, Solano, Chen. *Phys Rev Appl* 20:014024 | 2023 | Digitised counterdiabatic QAOA | Lattice turns | HP/MJ | DC-QAOA | — | Lattice peptides | shallower circuits than QAOA | Lattice | Depth-reduction line; not torsion space | **[S]** |
| 8 | Wong & Chang. *J Parallel Distrib Comput* 164:178 | 2022 | Grover-based | HP lattice | HP | Grover search | — | HP model | quadratic speedup claim | HP toy model | Not comparable | **[S]** |
| 9 | **Doga, Raubenolt, Cumbo, Joshi, DiFilippo, Qin et al. *JCTC* 20:3359** | 2024 | Perspective / review | survey of lattice + emerging | survey | survey | survey | — | — | — | The field's own framing document; same group as QTF | **[S]** |
| 10 | **Pamidimukkala, Bopardikar, Dakshinamoorthy, Kannan, Dasgupta, Senapati. *JCTC* 20:10223** | Nov 2024 | Gate-based VQE, "high degrees of freedom" turn encoding | Turn-based with **greater DOF**; HP alphabet | HP | VQE, **up to 114 qubits on IBM hardware** | Encoding, decoding | proteins of varied lengths | folded structures predicted | HP simplification; qubit count linear and large | Closest prior "more DOF than lattice turns" work; cited by QTF | **[S]** |
| 11 | **Li, Doga, Raubenolt, Mostame, DiSanto, Cumbo, ... Blankenberg. arXiv:2507.08955** | Jul 2025 | PolyFit + VQE-with-constraints (Lagrangian duality) | **FCC lattice** | contact potential | VQE-C on IBM Eagle R3 (ibm_cleveland) and Heron R2 (ibm_kingston) | Constraint handling | KLVFFA (6-mer) | FCC "improved capability for modeling realistic secondary structures" vs other lattices | Single 6-mer; lattice | Same group, same QPU, one year earlier — QTF's own predecessor | **[P]** abstract |
| 12 | **Linn, Li, Holden, Saki, DiFilippo, Radivoyevitch, Blankenberg, García-Álvarez, Johansson. arXiv:2509.18263** | Sep 2025 | Problem-agnostic ansatz QPSP | Tetrahedral / BCC / FCC lattices, 2nd-nearest-neighbour | classically computable energy cost function | Hardware-efficient problem-agnostic ansatz; IBM Kingston + noiseless simulator | Cost evaluation | peptides **up to 26 aa** | — (not in abstract) | Lattice; results not in abstract | Same group; establishes that "problem-agnostic ansatz + classically-computed cost" is their house style — exactly QTF's architecture | **[P]** abstract |
| 13 | **Zhang, Yang, Chen, Lu, Saeidi, Volchenboum, Zhao, Chen, Jiang, Guan. arXiv:2510.06413** | Oct 2025 | Hybrid quantum-AI energy fusion | VQE conformations + **NSP3 (NetSurfP-3) secondary-structure and dihedral distributions as statistical potentials** | fused quantum + learned potential | VQE on IBM 127-qubit | NSP3 network, fusion, ranking | **375 conformations from 75 protein fragments** | **mean RMSD 4.9 Å**, p<0.001 vs AlphaFold3, ColabFold, quantum-only | Coarse quantum landscape by their own account | **The nearest analogue of our "sequence-conditioned torsion prior + VQE" idea, published Oct 2025.** Our 3.213 Å beats their 4.9 Å on a comparable fragment benchmark | **[P]** abstract |
| 14 | **Kannan, Pamidimukkala, Dakshinamoorthy, Bopardikar, Dasgupta, Senapati. arXiv:2510.15316** | Oct 2025 | **CVaR-VQE** turn encoding | FCC lattice, 5(N−1) qubits | Miyazawa–Jernigan | CVaR-VQE on IBM 133-qubit Heron | Decoding, comparison | 13 peptides 6–10 residues, chains to 20 | **RMSD 1.22–3.11 Å** (their classical MD control 0.79–5.32 Å) | Lattice; CG; MJ potential | **The number our <2.0 Å target will be read against.** Establishes CVaR-VQE on lattice QPSP as taken | **[P13]** |
| 15 | Zhang, Yang, Martin, Lin, Wang, Lu et al. *Adv Sci* 13:e13641 (arXiv:2506.22677) | 2026 | Hardware-executable quantum framework | binding-site structure | — | utility-level QPU | — | protein binding sites | — | — | Binding-site niche; not backbone folding | **[S]** |
| 16 | **Zhang, Fang, Yang, Cheng, Chen, Fang, Chen, Zhao, Guan. QSAD. arXiv:2607.06971** | Jul 2026 | Non-iterative Hamiltonian evolution sampling | amino-acid-level Hamiltonian sampling | — | IBM Heron R2 | — | **101 binding-pocket peptides, 5–18 residues** | "27–71 % improvement" over AI and quantum baselines | Binding-pocket regime; improvement is relative | **The only quantum peptide benchmark at our scale (n=101, 5–18 res).** Watch it | **[P]** abstract |
| 17 | Feng, Zhang, Huang, Fang, Xu, Guan, Zhang. Graph-VQE. arXiv:2607.02749 | Jul 2026 | Louvain graph partitioning of the Hamiltonian into weakly coupled clusters; multi-QPU CUDA-Q | lattice-style interaction graph | — | simulated multi-QPU | Louvain community detection | not specified in abstract | "lower final energies" than baselines | **Reports energy, not RMSD** | A live example of the failure our project names: optimising the objective and never pricing the structure | **[P]** abstract |
| 18 | Agathangelou, Manawadu, Tavernelli. arXiv:2507.19383 | 2025 | Quantum side-chain optimisation vs classical | fixed backbone, rotamer choice | rotamer energy | VQE/QAOA | comparison | side-chain packing | quantum vs classical comparison | Fixed backbone | **Has the classical comparison QTF lacks** — a model for how to do it | **[S]** |
| 19 | **Boulebnane, Lucas, Meyder, Adaszewski, Montanaro. *npj Quantum Inf* (2023)** | 2023 | QAOA, noiseless, 20 qubits | lattice protein | folding objective | QAOA | random-sampling control | simplified folding | *"the performance of QAOA can be matched by random sampling up to a small overhead... these results cast serious doubt on the ability of QAOA to address the protein folding problem in the near term, even in an extremely simplified setting"* | — | **The mandatory control.** Neither QTF nor most of this table runs it | **[P13]** |
| 20 | Marchand et al. | 2018 | Quantum annealing in torsion space | rotatable-bond torsions | *"flexible with respect to the choice of molecular force field"* | D-Wave | — | small molecules | — | Small molecules, not peptides | Prior torsion-space quantum optimisation **with a real force field** — predates QTF by 8 years | **[P13]** |
| 21 | Mato et al. Quantum Molecular Unfolding | 2022 | HUBO over rotatable bonds | torsions | force field | annealing / QAOA | — | ligands | — | Ligand-scale | Torsion-space HUBO prior art | **[P13]** |
| 22 | Shajan, Kaliakin, Liang, Merz et al. (SQD / SQD+DMET / implicit-solvent SQD); Merz et al. arXiv:2605.01138 | 2025–26 | Sample-based quantum diagonalization | all-atom quantum chemistry on proteins | electronic energy | genuine quantum chemistry on QPU | DMET embedding | protein–ligand complexes, to 12,000 atoms | electronic energies | **Scores conformers; does not predict folds** | The other meaning of "quantum protein" — orthogonal to us | **[S]** |

### 2B. Variational quantum algorithms — objective, trainability, controls

| # | Citation | Date | Contribution | Relation to us | Ver |
|---|---|---|---|---|
| 23 | **Barkoutsos, Nannicini, Robert, Tavernelli, Woerner. *Quantum* 4:256** | 2020 | **CVaR as a VQE objective** for combinatorial optimisation | The origin of our CVaR arm. **Not cited by QTF.** CVaR-VQE is established prior art in general and in lattice QPSP (row 5, 14) — but not in continuous torsion space | **[P13]** |
| 24 | Kandala, Mezzacapo, Temme, Takita, Brink, Chow, Gambetta. *Nature* 549:242 | 2017 | Hardware-efficient VQE; the EfficientSU2 lineage | QTF's ansatz; also the ansatz family our `StatevectorCircuit` sits in | **[S]** |
| 25 | McClean, Boixo, Smelyanskiy, Babbush, Neven. *Nat Commun* 9:4812 | 2018 | Barren plateaus from 2-design expressibility | Regime of validity: needs the circuit ensemble to 2-design. Our MPS/shallow arms are outside it | **[P13]** |
| 26 | Cerezo, Sone, Volkoff, Cincio, Coles. *Nat Commun* 12:1791 | 2021 | Cost-function-dependent BP: global observables → exponential decay even at shallow depth | **Our measurement refutes its relevance at 6–18 qubits**: the measured ansatz kernel is flat in Pauli weight; what matters is n | **[P13]** |
| 27 | Anschuetz & Kiani. *Nat Commun* (2022) | 2022 | Traps, not plateaus: BP-free shallow models still have superpolynomially few good local minima | **Absence of a plateau is not evidence of trainability.** Solution quality must be its own axis — which is exactly our "energy reached vs RMSD, ρ=−0.006" | **[P13]** |
| 28 | **Qiu, Lumbreras, Li, Rebentrost. arXiv:2605.02850** | 2026 | Tilted-loss family (containing CVaR) *"does not remove the barren plateau problem by itself"*; bottleneck moves *"from trainability to estimability"* | **CVaR must not be justified as BP mitigation in Sprint 14.** Justify it as a tail objective | **[P13]** |
| 29 | Cerezo et al. *Nat Commun* (2025), simulability critique | 2025 | Many provably BP-free models admit classical simulation after a classical data-acquisition phase — "a soft form of dequantization" | Live against our exact-statevector and exact-MPS arms; unaddressed by us and by QTF | **[P13]** |
| 30 | Ragone et al., Lie-algebraic variance unification | 2024 | Var = P_M(ρ)·P_M(O)/dim(M) | The language for our no-free-parameter spectrum→gradient-variance result (measured/predicted 1.006, 1.001) | **[P13]** |

### 2C. Chemical-shift → torsion / structure channel (for the parallel shift agent)

| # | Citation | Date | What it maps | Accuracy | Coverage / caveat | Ver |
|---|---|---|---|---|---|---|
| 31 | Cornilescu, Delaglio, Bax. TALOS. *J Biomol NMR* 13:289 | 1999 | shifts → (φ,ψ) | — | ~65 % coverage, ~3 % error, 20-protein DB | **[P13]** |
| 32 | Shen, Delaglio, Cornilescu, Bax. **TALOS+**. *J Biomol NMR* 44:213 | 2009 | shifts → (φ,ψ) | **±13°** excl. ~2.5 % gross errors | 88.5 % "Good"; 2.46 % error; **reference DB excludes chains ≤20 residues** | **[P13]** |
| 33 | **Shen & Bax. TALOS-N.** *J Biomol NMR* 56:227 | 2013 | shifts → (φ,ψ) **and χ1** | **~12° RMS** vs crystallographic (φ,ψ); **χ1 rotameric state for ~50 % of residues at 89 % consistency** | **≥90 % of residues**, ~3.5 % error rate. Post-filter figure. Gating: ≥3 shifts on ≥2 of {i−1,i,i+1}; heptapeptide (±3) matching unit; **RCI-S² ≤ 0.6 → `Dyn`, no restraint** | **[P13]** + χ1 detail **[S]** |
| 34 | **Shen & Bax. SPARTA+.** *J Biomol NMR* 48:13 | 2010 | **coordinates → shifts** (the inverse map) | RMS 2.45 / 1.09 / 0.94 / 1.14 / 0.25 / 0.49 ppm for ¹⁵N/¹³C′/¹³Cα/¹³Cβ/¹Hα/¹HN | Trained on 580 high-res X-ray structures. **Not a torsion predictor** | **[P13]** |
| 35 | Li, Brüschweiler. UCBShift | 2020 | coordinates(+sequence) → shifts | RMSE 0.31 (HN), 0.19 (Hα), 0.84 (C′), 0.81 (Cα), 1.00 (Cβ), 1.81 (N) ppm | Transfer module uses sequence + structural alignment to reference candidates | **[S]** |
| 36 | **LEGOLAS.** *JCTC* (2025), doi 10.1021/acs.jctc.5c00026 | 2025 | coordinates → shifts, TorchANI GNN | RMSE **2.53 (N), 0.91 (Cα), 1.14 (Cβ), 1.02 (C′), 0.49 (HN), 0.27 (Hα) ppm** | Order of magnitude faster than SHIFTX2; open source. **Also an inverse map — not a torsion predictor** | **[S]** |
| 37 | **Klukowski, Riek, Güntert. ARTINA.** *Nat Commun* 13:6151 (2022); and *Sci Adv* 9:adi9323 (2023) | 2022–23 | raw NMR spectra → assignments → restraints → structure | **1.44 Å median RMSD to PDB reference**; **91.4 %** correct chemical-shift assignments | **100-protein benchmark, 1,329 spectra, monomeric globular proteins.** 2023 integrative version with AlphaFold + UCBShift: 92.59 % from five 3D spectra | **[S]** |
| 38 | SPOT-1D-LM. *Sci Rep* (2022) | 2022 | sequence (LM) → (φ,ψ) | φ MAE **15.99°** / ψ **23.74°** (TEST2018); 20.67 / 36.57 (TEST2020). Single-seq SPOT-1D-Single: 22.16/40.58 and 22.92/44.25 | Full crystallised domains only; **no 9–16-mers in any test set**; a 13-mer is the Neff≈1 regime | **[P13]** |

### 2D. Geometry, force fields, and reconstruction (QTF's classical stack, and ours)

| # | Citation | Date | What | Relation to us | Ver |
|---|---|---|---|---|---|
| 39 | Parsons, Holmes, Rojas, Tsai, Strauss. *J Comput Chem* 26:1063 | 2005 | **NERF** — torsion → Cartesian | QTF's reconstructor. Ours is `core/geometry.build_backbone_batch` on ideal geometry — the same class of object | **[S]** |
| 40 | AlQuraishi. pNERF. *J Comput Chem* 40:885 | 2019 | parallelised NERF | The batched form; ours is batched too | **[S]** |
| 41 | Engh & Huber. *Acta Cryst* A47:392 | 1991 | ideal bond lengths/angles | QTF fixes geometry to these; our builder does the same. **Shared assumption, shared ceiling** | **[S]** |
| 42 | Alford, Leaver-Fay, Jeliazkov, O'Meara, DiMaio, Park et al. *JCTC* 13:3031 | 2017 | Rosetta REF15 all-atom energy | QTF's second backend. **Its final-model native-like hit rate on chignolin was 0.3 %** | **[S]** |
| 43 | Maier, Martinez, Kasavajhala, Wickstrom, Hauser, Simmerling. *JCTC* 11:3696 | 2015 | **AMBER ff14SB** | **The identical force field we run** via `core/amber.py`. QTF's OpenMM backend | **[S]** |
| 44 | Eastman, Swails, Chodera, McGibbon, Zhao, Beauchamp et al. *PLoS Comput Biol* 13:e1005659 | 2017 | OpenMM 7 | **The identical engine we run** | **[S]** |
| 45 | Lindorff-Larsen, Piana, Palmo, Maragakis, Klepeis, Dror, Shaw. *Proteins* 78:1950 | 2010 | amber99sb-ildn | QTF's GROMACS post-minimisation force field | **[S]** |
| 46 | Honda, Akiba, Kato, Sawada, Sekijima, Ishimura et al. *JACS* 130:15327; PDB 5AWL | 2008/15 | chignolin crystal structure | QTF target 1 | **[S]** |
| 47 | Barua, Lin, Williams, Kummler, Neidigh, Andersen. *PEDS* 21:171; PDB 2JOF | 2008 | Trp-cage | QTF target 2 | **[S]** |
| 48 | Miyazawa & Jernigan. *J Mol Biol* 256:623 (1996); *Macromolecules* 18:534 (1985) | 1985/96 | MJ contact potentials | The energy of nearly all lattice QPSP, including rows 5, 14 | **[S]** |
| 49 | Dill. *Biochemistry* 24:1501 (1985); Lau & Dill. *Macromolecules* 22:3986 (1989) | 1985/89 | HP model | The other lattice energy | **[S]** |
| 50 | Zhou & Zhou. DFIRE. *Protein Sci* 11:2714 | 2002 | distance-scaled statistical potential | QTF's PHEAT can generate DFIRE-like potentials. **Closest published analogue of our Legacy energy** | **[S]** |

**Count: 50 rows, ≥25 required.** 23 rows are primary-verified (this sprint or Sprint 13).

---

## 3. NOVELTY BOUNDARY

### 3.1 What is now definitively TAKEN

Each of these was, at some point, part of this project's implicit novelty story. Each is
published. Name the paper and stop claiming it.

| Claim | Taken by | Verdict |
|---|---|---|
| Quantum torsion-space (off-lattice) protein folding | QFold, arXiv:2101.10279 / *QST* 7:025013 (2021–22); Marchand 2018; Mato 2022; Pamidimukkala *JCTC* 2024 | **TAKEN.** Do not claim "first off-lattice / torsion-space quantum folding" in any form |
| Per-residue binned φ/ψ + binary index encoding + VQE at peptide length | arXiv:2510.06413 (Oct 2025), IBM 127-qubit, 75 fragments | **TAKEN** (Sprint 13 finding C1, re-confirmed). Our `PerResidueTorsion` with `n_bits = n_res·log₂k` is that object |
| **CVaR-VQE applied to protein folding** | Robert et al. *npj QI* 7:38 (2021); Kannan et al. arXiv:2510.15316 (Oct 2025) | **TAKEN on lattice representations.** See 3.2 for what survives |
| Torsion prior from a neural network fused into a quantum folding objective | arXiv:2510.06413 (NSP3 dihedral distributions as statistical potentials); QFold's Minifold initialiser (2021) | **TAKEN** |
| Real all-atom AMBER force field as the scorer for quantum-generated peptide conformations | QTF's OpenMM ff14SB backend, arXiv:2609.02113 | **TAKEN, as of three days ago** |
| Continuous (non-binned) torsions from a log-scaled quantum register | QTF, arXiv:2609.02113 | **TAKEN.** This is QTF's genuine contribution |
| All-heavy-atom side-chain reconstruction from a quantum-generated torsion vector | QTF, arXiv:2609.02113 | **TAKEN** |
| "Energy and RMSD are decoupled in quantum peptide folding" as a *qualitative* observation | QTF, arXiv:2609.02113 §Results/Discussion; and it is folklore in scoring-function literature | **TAKEN qualitatively.** See 3.2 |
| Snapshot/trajectory mining to beat endpoint selection | QTF, arXiv:2609.02113 (7.3 M snapshots) | **TAKEN** |
| Quantum peptide structure prediction at n≈100 targets of 5–18 residues | QSAD, arXiv:2607.06971 (101 binding-pocket peptides) | **TAKEN as a benchmark scale** |

### 3.2 What genuinely survives the audit

These are the claims I could not find published anywhere after aggressive search. Each is
stated at the tier the evidence supports.

1. **An exact locality theorem for torsion-space molecular objectives.** Our result — under an
   ideal-geometry backbone builder, `d_ij` depends on *exactly* the `j−i−1` residues strictly
   between i and j, agreement 1.0000, zero counterexamples at machine precision, holding at
   termini, on GLY/PRO and under cis-ω; with all-atom supports one residue wider (CA `i<m<j`,
   N `i≤m<j`, C/O `i<m≤j`, CB `i≤m≤j`); and the corollary that a separation-8 pair is a
   14-qubit interaction at k=4, so **no 2-local Ising form of a distance-based molecular
   objective exists in this encoding** — has no counterpart in this literature. QTF does not
   analyse locality at all; it does not need to, because it never forms a qubit Hamiltonian.
   **SURVIVES, and it is our strongest formal claim.** HYPOTHESIS→ our own DEMONSTRATED.
2. **A no-free-parameter chain from the Pauli spectrum to the gradient variance**
   (measured/predicted 1.006 and 1.001), together with the **delta-spike artefact correction**
   (raw AMBER's top-10 of 4,096 configurations carry a median 99.6 % of its Walsh variance,
   landing on Binomial(m,½) exactly; 99th-percentile winsorisation is insufficient). Nothing
   in the QPSP literature computes a Pauli/Walsh spectrum of a molecular energy at all.
   **SURVIVES.**
3. **A quantitative, certified statement of the energy–structure mismatch.** QTF establishes
   the phenomenon and reports *no numeric correlation coefficient anywhere*. We have
   ρ = +0.043 (Legacy, low-energy decile), −0.088 (raw AMBER), the native at the 32nd–40th
   percentile, and — uniquely — **a certified global optimum on 9 fully enumerated targets
   (262,144 configurations each, 2.36 M labelled structures) that is +0.139 Å WORSE than random
   sampling**, plus the demonstration that monotone log compression leaves ρ unchanged while
   collapsing AMBER's range from 16.1 to 1.8 decades. **The certified-optimum result is the one
   nobody else can make**, because nobody else has an enumerable state space. **SURVIVES, and
   it is publishable on its own.**
4. **In-loop CVaR in a continuous/high-resolution torsion representation.** CVaR-VQE is taken
   on lattices; QTF has no CVaR at all and its tail mechanism is post-hoc oracle-selected
   best-of-N. **SURVIVES narrowly** — and only if we state the boundary honestly as "CVaR in
   torsion space", never as "CVaR in protein folding".
5. **A causal VQE-vs-classical control at matched evaluation budget, on a molecular objective.**
   Row 18 (Agathangelou/Tavernelli, side-chain packing) is the only quantum-vs-classical
   comparison I found in the whole QPSP table, and it is on a fixed backbone. **QTF has none.**
   Doing this properly, and reporting it whichever way it falls, **SURVIVES** — and per BRIEF
   §1.6 a rigorous negative here is a success.
6. **A protected, pre-registered, held-out predictive benchmark with pinned folds and identity
   clusters, paired bootstrap CIs, W/L, per-fold values and drop-top-10/20 concentration
   checks.** **There is not a single confidence interval or p-value in QTF's 50 pages**, and
   n=2 targets. Our statistical discipline is not merely better; it is of a different kind.
   **SURVIVES, and it is the cheapest differentiator to defend.**
7. **The "optimising harder makes structures worse" curve** — 3.764 Å at 10 evaluations,
   3.667 at 300, **3.920 at the certified global optimum**, governed by the fraction of space
   explored rather than the evaluation count. **SURVIVES**, and it is the sharpest available
   critique of every energy-minimising quantum folding paper in §2A.
8. **The measured decomposition of what a sequence-conditioned torsion library actually buys**
   — 88 % generic Ramachandran, 12 % sequence-related; a wrong target's library costs only
   +0.135 Å; the entire measurable sequence channel is 10.4° of ψ; the incumbent pipeline is
   equivalent to σ≈29°. **SURVIVES**, and it directly prices the Ramachandran wells QTF builds
   into its custom energy.

### 3.3 What does NOT survive, and must be removed from the story

* Any framing of this project as "the first quantum method to fold peptides off-lattice in
  torsion space." Three separate lines of prior art (QFold 2021, Marchand 2018, Pamidimukkala
  2024), plus QTF itself.
* "Torsion-constrained VQE for peptides" as a novel object (Sprint 13 finding C1 stands).
* CVaR as a barren-plateau mitigation (Qiu et al. 2026; Sprint 13 finding C2 stands).
* Any claim that we are the first to notice energy does not rank structure. QTF says it in
  print, in September 2026, and so do decades of decoy-discrimination papers. **What we own is
  the certified, quantitative form.**

---

## 4. DIFFERENTIATION — ranked by strength

Ranked by how hard each would be for a referee, or for the QTF group, to dismiss.

**Tier 1 — differences of scientific kind, not degree.**

1. **The causal control.** We run VQE against classical optimisers at *matched evaluation
   budget* on the *same* objective and *same* representation, and report it whichever way it
   falls. QTF cannot attribute any of its result to the quantum component and does not try.
   This is the difference between a demonstration and an experiment. Concretely: their
   architecture is a classical COBYLA/SLSQP search over **P ≈ 2N + 6n = 132** parameters that
   produces **43** torsions; the obvious control — the same optimiser on the 43 torsions
   directly — is absent from 50 pages. **Do this control in Sprint 14 and name it in the arm
   name** (BRIEF §1.4).
2. **The enumerated state space.** `s13/results/qarch_enum_<PDB>.npz` — nine targets, 262,144
   configurations each, 2.36 M structures with true RMSD, Legacy, prior, eleven Legacy
   components, plus AMBER on a 2,955-config subset. This gives us *certified* global optima and
   exact spectra. No paper in §2A has an enumerable space. Every claim we make about "the
   optimum of this objective" is a theorem about a finite set; every claim QTF makes about its
   landscape is an inference from a biased sample of 7.3 M points drawn *by the optimiser it is
   trying to characterise*.
3. **Statistical discipline.** Paired bootstrap CIs, W/L, per-fold values, drop-top-10/20
   concentration checks, null controls (permutation, label shuffle, matched random,
   sequence-blind twin) on every attractive correlation, pre-declared selection logic, a
   pinned protected benchmark. Against n=2 targets and zero CIs.
4. **Explicit separation of objective quality from structural accuracy as a reporting axis.**
   We have the number: within-target ρ(energy reached, RMSD) across six optimisers is
   **−0.006**, and SPSA optimises the AMBER objective best of five arms while returning the
   **worst** structure. QTF observes the phenomenon; we have made it a measurement protocol.

**Tier 2 — differences of instrument.**

5. **Scale and honesty of the benchmark.** 126 development targets with five pinned
   leave-fold-out folds, a cluster-disjoint dev24, and a protected held-out benchmark evaluated
   once. QTF: two of the most-optimised miniproteins in structural biology, both with the
   method's own bias terms (a "hairpin" bias, α/β Ramachandran wells) plausibly matched to
   target 1. Our failure taxonomy alone (10/18 catastrophic targets are fibril segments or
   lasso peptides, Fisher p = 1.2e-6) is more diagnostic information than QTF reports in total.
6. **Two genuine energies with a measured relationship, not three interchangeable backends.**
   We have Legacy (non-all-atom) and AMBER ff14SB/GBn2 via OpenMM, with their Pauli-weight
   spectra measured (2.236 vs 3.015 after monotone conditioning, AMBER higher on 79/79 cells),
   their locality settled (**both full-register — "AMBER is less local than Legacy" is a
   category error and is refuted**), and their component decompositions (Legacy's entire skill
   is clash rejection; remove the steric term and it inverts to +0.395 Å worse than random).
   QTF's custom energy has hand-set weights, no ablation, and its algebra is in an unavailable
   supplement.
7. **A representation whose search space we can characterise.** Our k-state torsion library at
   k=4/8/16 has known live-qubit counts (25.9 / 38.9 / 51.8), known ORACLE descent ceilings
   (1.594 / 1.184 / 0.876) and a known random-start gap (1.982 from random vs 1.594 from the
   oracle start — 0.388 Å of the ceiling is the privileged start). We know what our
   representation can and cannot reach. QTF reports a hypothetical bin count in a table
   labelled "illustrative reference."

**Tier 3 — differences of care that referees notice.**

8. **A defined, documented RMSD protocol.** Full-chain Kabsch Cα RMSD over all residues via
   `I.ca_rmsd`, reproducible to 16 digits (`python -m s12.instrument` must emit
   `shipped 3.4540004952559396`, `pool_best 1.7108244199364904`,
   `synthesis_fit 3.2040761603809194`). QTF's Methods never define the RMSD; terminal exclusion
   surfaces only in a hardware aside.
9. **Known-defect disclosure.** We document that `core/project.py` leaves φ[0], ψ[n−1] and
   φ[n−1] inert, that residue n−1 is entirely invisible to Cα-RMSD, that log₂(k) qubits per
   chain are dead, and that every nominal `n_res·log₂k` count is an overcount — and we report
   live qubits. QTF reports 6 logical qubits for 44 torsions with no such audit.
10. **The corrections ledger.** BRIEF §1.5 requires corrections to be preserved next to the
    finding they correct. Our memory index carries a dozen explicit OVERTURNED / REINSTATED /
    CORRECTED entries. That is a research-integrity artefact no single paper can match.

**What we should learn from QTF without copying it.**

* **Report the median, not the minimum — and report the sampled-vs-selected gap explicitly.**
  Their 99 % pool hit rate against a 0.3–6 % selection hit rate is the single most informative
  table in the paper *because* they computed both. We have the analogous numbers
  (`pool_best 1.711`, `top75_best 2.306`, `synthesis_fit 3.204`) and should present them the
  same way: **sampling is solved, selection is not.**
* **Trajectory snapshot retention is cheap and we do not do it.** They kept up to 5,000
  low-energy structures per replica, deduplicated at ≥0.25 energy units, and it converted a
  1.199 Å result into a 0.623 Å one. We have `BestSeenTracker` in `core/quantum.py`; a
  **deduplicated low-energy snapshot pool per VQE run**, retained and reported separately from
  the returned structure, is directly transferable. **But it is an ORACLE DIAGNOSTIC unless the
  selection back out of the pool is native-free** — and our Sprint 13 result says nothing ranks
  within the pool (dev pool 2.355 Å vs returned 3.324 Å; all-atom AMBER reranking worth
  +0.004 Å). So adopt the *reporting*, not the claim.
* **The two-decoder discipline.** They kept the simulator decoder and the hardware decoder
  mathematically distinct and *said so*, repeatedly, rather than pretending the hardware
  reproduced the simulation. That is exactly the honesty BRIEF §1.4 demands, and it is worth
  imitating in how we describe any hardware or MPS arm.
* **The ω treatment.** Constraining ω to a 170–190° window instead of spending optimiser effort
  on cis is a cheap, information-free, defensible restriction. Worth checking whether our
  torsion library already implies it.
* **Don't adopt:** the phase encoding (see §5), the hand-weighted custom energy, the
  end-to-end bias term, or oracle-selected headline numbers.

---

## 5. RISK — does any of this invalidate, weaken, or change Sprint 14?

**Direct answer: nothing here invalidates the Sprint 14 direction. Two things change the
framing, one thing changes a priority, and one thing is a genuine threat to a claim we might
otherwise have made.**

### R1 — LOW RISK, and actually a strengthening. The energy-ranking negative is now corroborated.
QTF independently reproduces our central finding, on different targets, with three energy
functions, over 7.3 M structures, and reports OpenMM/AMBER energy *negatively* correlated with
RMSD. **This de-risks Sprint 14's most uncomfortable result.** It is no longer plausible that
"neither energy ranks the native" is a bug in our pipeline. It also means BRIEF §2's statement
should now cite arXiv:2609.02113 as external corroboration. **Action: cite it; do not re-run
anything.**

### R2 — MEDIUM. The novelty story must be rewritten, but the science does not change.
QTF lands three days before this sprint and takes "continuous torsions from a log-scaled
quantum register, all-atom, with a real force field backend." Combined with QFold (2021),
arXiv:2510.06413 (Oct 2025) and Kannan (Oct 2025), the *representational* novelty of this
project is now essentially zero. **What survives is entirely in §3.2 — the locality theorem,
the spectrum→variance chain, the certified global optimum, the causal control, and the
statistical instrument.** These are stronger claims than the representational ones ever were.
**Action: reposition Sprint 14's contribution from "a new quantum encoding" to "the first
properly controlled measurement of whether quantum conformational optimisation helps, with a
certified landscape." No experimental change.**

### R3 — MEDIUM, and this is the priority change. Their result raises the bar on *reporting*, not on accuracy.
Naive reading: "QTF gets 0.623 Å, we get 3.213 Å, we are far behind." **This reading is wrong
and the sprint coordinator must not act on it.** Their 0.623 Å is an oracle argmin over 2.68 M
structures on one 10-residue hairpin with terminal residues excluded from the metric. Their
predictive number is a **2.90 Å median final model on chignolin and 5.39 Å on Trp-cage**, and
their published hybrid-quantum-AI comparator (arXiv:2510.06413) reports **mean 4.9 Å across 75
fragments**. Our 3.213 Å across 126 targets, full-chain Cα, no oracle selection, is
**competitive with or better than every predictive number in this literature.** The risk is
that we chase an oracle number. **Action: state the sampled-vs-selected gap in our own results
the way they do, and never quote a pool minimum as performance — BRIEF §1.2 already forbids it.**

### R4 — HIGH, and it is a threat to a specific claim. The log-encoding is a reparameterisation, and if we borrow it we inherit the flaw.
If anyone in Sprint 14 is tempted by the phase encoding as a way to cut qubit counts, the
arithmetic in §1.7 must be read first: **P = 2n(⌈N/n⌉+3) ≈ 2N + 6n > N.** The classical
parameter vector is larger than the torsion vector. There is no compression of the search
problem, only of the register. Worse, the objective is never a qubit observable, so **there is
no VQE in the technical sense** — no Hamiltonian, no expectation value, no variational
principle over an operator. Adopting it would violate BRIEF §1.4 ("Do not call a classical
search 'VQE'") unless described with great care. **Action: do not adopt the phase encoding.
If it is evaluated at all, evaluate it as a classical reparameterisation arm and name it as
one.**

### R5 — MEDIUM-LOW, an unpriced opportunity and a trap. Snapshot mining.
Their 99 % → 0.3–6 % gap says the whole field's bottleneck is selection, which is also our
Sprint 13 conclusion ("nothing ranks within the pool"; "consensus is outlier avoidance, not a
nativeness signal"; "in-band discrimination is signal-limited, not sample-limited"). Retaining
a deduplicated low-energy snapshot pool per VQE run is cheap and improves our *diagnostic*
resolution. **The trap is that the improvement is oracle-only** — our own instrument already
says the selection back out of the pool buys +0.004 Å with all-atom AMBER. **Action: add
snapshot retention as instrumentation; label any pool minimum ORACLE DIAGNOSTIC.**

### R6 — LOW but watch it. Benchmark scale is being matched.
QSAD (arXiv:2607.06971, Jul 2026) reports 101 binding-pocket peptides of 5–18 residues on IBM
Heron R2. Our "n=126 vs their n=2" differentiator is safe against QTF but not against the whole
field for much longer. **Action: none this sprint; the protected held-out benchmark is the
durable asset, not the target count.**

### R7 — LOW. An unresolved question about their work that would matter if it went the other way.
Whether QTF's `E_e2e` end-to-end bias target is derived from the experimental structure is
**UNVERIFIED** (§1.5). If it is, their custom-energy chignolin results are leaked and the
comparison shifts further in our favour. If it is not, nothing changes for us. It is
checkable in `github.com/cumbof/qtf`. **Action: optional; a 20-minute source read, not a
sprint dependency.**

### Things that specifically do NOT threaten us
* Their hardware results. 6 qubits, no optimisation on device, warm-started from good classical
  answers, a different decoder, 0/10 improvement on the best inputs, and an explicit
  acknowledgement that the confound is unresolvable in their design.
* Their qubit-count advantage. It is width traded for depth (455 median transpiled depth, 204
  CZ gates for a 10-residue problem) plus an explicit admission that shot cost "effectively
  trade[s] the problem of qubit quantity for a problem of measurement complexity."
* Their claim of being "fully realized, scalable... not a preliminary proof-of-concept." Two
  targets, no controls, no statistics.

---

## 6. VERIFICATION STATUS

### Confirmed from a primary source I fetched and read in full this sprint
* **arXiv:2609.02113 exists, and I read the complete 50-page text.** Retrieved as PDF from
  `https://arxiv.org/pdf/2609.02113`, extracted with `pdftotext -layout` (124,283 characters,
  1,801 lines), and read end to end. Title, all six authors, Cleveland Clinic affiliations,
  1 Sep→2 Sep 2026 submission, v1 only, 50 pages/8 figures, CC BY 4.0, corresponding author
  blanked2@ccf.org — all confirmed against the arXiv abstract page independently.
* Everything in §1.1–§1.12 marked LITERATURE-SUPPORTED: the encoding equations, the CDF decoder
  equations, `n = max(2,⌈log₂N⌉)`, `d = ⌈N/n⌉+2`, `P = 2n(d+1)`, the 72-parameter worked
  example, EfficientSU2 + brickwork, COBYLA/SLSQP, the three-stage curriculum and its weights
  (50.0, 15.0, 8.0→1.5, 15.0→2.5), SHA-256 seeding, best-of-50 scouting, 5,000-snapshot
  retention at ≥0.25 energy units, the term list of the custom energy, the log-capped repulsion
  example, all three energy backends, GROMACS amber99sb-ildn post-minimisation, Table 1
  (1,157/2,681,469 and 1,176/4,650,631; total 7,334,433), Table 2 (44 and 88 torsions), Table 4
  (6 logical qubits; depths 443–506/206; CZ 168–208/60), 8,192 shots, 32 twirls + XY4, the
  device identities (ibm_cleveland = 156-qubit Heron R2 heavy-hex; ibm_miami = 120-qubit
  Nighthawk R1 square), 88/300 and 45/300, the per-backend hardware medians, 1.758/1.782 Å, the
  0/10 custom improvement, the hit rates 100/98/99.3 % vs 6/4/0.3 %, 0.623/1.199/2.501/3.512 Å,
  the median snapshot improvements 2.18/2.66/3.26 Å, the ~67 OpenMM failures, and all quoted
  sentences.
* **CVaR is absent.** Verified by exhaustive case-insensitive search of the complete extracted
  text for "cvar", "conditional value", "quantile", "tail". Zero substantive hits. Barkoutsos
  et al. is absent from the 50-reference list, which I read in full.
* **No classical control and no CI/p-value.** Verified by full read of Methods, Results and
  Discussion.
* **They do report Spearman energy–RMSD correlations** (Figure 5) and state the OpenMM
  negative-correlation finding verbatim; **and they report no numeric ρ anywhere.** Both
  verified by full read.
* **The RMSD protocol is undefined in Methods**; terminal exclusion appears only in the
  ibm_miami discussion. Verified by targeted search of the whole text for RMSD-definition
  language.
* arXiv:2510.06413 — title, ten authors, 7 Oct 2025, verbatim abstract, VQE on IBM 127-qubit,
  NSP3 statistical potentials, 375 conformations from 75 fragments, **mean RMSD 4.9 Å**.
  Fetched from the arXiv abstract page. **Correction to Sprint 13:** its abstract says **NSP3**,
  not "NetSurfP", and **does not mention CVaR**; Sprint 13's row implying CVaR for this paper
  should be treated as unconfirmed.
* arXiv:2510.15316 — title "Capturing Protein Free Energy Landscape using Efficient Quantum
  Encoding", six authors (Kannan, Pamidimukkala, Dakshinamoorthy, Bopardikar, Dasgupta,
  Senapati), 17 Oct 2025, FCC lattice, turn-based encoding, VQE on IBM 133-qubit, compared
  against classical simulated annealing and MD. **The 1.22–3.11 Å figure and the 5(N−1) qubit
  formula are carried from Sprint 13, not re-verified this sprint.**
* arXiv:2509.18263 — title, nine authors, 22 Sep 2025, tetrahedral/BCC/FCC lattices with
  second-nearest-neighbour interactions, problem-agnostic hardware-efficient ansatz, IBM
  Kingston + noiseless simulator, peptides up to 26 aa, 18 pages/8 figures.
* arXiv:2507.08955 — title, twenty authors, 11 Jul 2025, FCC lattice, PolyFit and VQEC via
  Lagrangian duality, IBM Eagle R3 (ibm_cleveland) and Heron R2 (ibm_kingston), KLVFFA 6-mer.
* arXiv:2607.06971 (QSAD) — title, nine authors, 8 Jul 2026, IBM Heron R2, 101 binding-pocket
  peptides of 5–18 residues, "27–71 % improvement". **Abstract-level only.**
* arXiv:2607.02749 (Graph-VQE) — title, seven authors, 2 Jul 2026, Louvain partitioning,
  CUDA-Q multi-QPU simulation, reports lower final energies. **Abstract-level only.**
* The 50-item QTF reference list, read in full, is the source for rows 2–4, 7–10, 15, 22, 24,
  39–49 of the ledger — **bibliographically reliable, content-wise secondary.**

### Carried forward from Sprint 13 as primary-verified there, not re-fetched here
TALOS/TALOS+/TALOS-N figures and gating mechanisms; SPARTA+ RMS values; SPOT-1D-LM MAEs;
Robert et al.; Kannan's RMSD range and qubit formula; Boulebnane et al.'s quotation; Barkoutsos
CVaR-VQE; Qiu et al. 2026; McClean/Cerezo/Anschuetz-Kiani/Ragone; Marchand 2018; Mato 2022;
QFold's bits-per-angle. These are marked **[P13]** in the ledger. I re-read
`s13/lit_FINDINGS.md` §0–§1 to confirm what it claims, but did not re-fetch its sources.

### Secondary or unverified — do not quote as fact
* **The algebra of QTF's custom energy function.** Deferred to a supplement I could not obtain.
  Term *names* and *stage weights* are primary; the functional forms are **UNVERIFIED**.
* **Whether `E_e2e`'s target distance is native-derived.** **UNVERIFIED.** The single most
  consequential open question about the paper.
* **Figures S1 and S2** (ibm_miami panels; side-chain DOF map): not obtained.
* ~~Table 3's exact cell values are uncertain due to PDF layout interleaving.~~
  **CORRECTED, same session:** re-extracting pages 19–20 in raw reading order recovers Table 3
  cleanly, and four cells cross-check exactly against the abstract's three-decimal numbers
  (0.62=0.623, 1.20=1.199, 2.50=2.501, 3.51=3.512). **Table 3 is now
  LITERATURE-SUPPORTED without caveat**, including the E2E/Rg rows. Correction retained per
  BRIEF §1.5.
* **QFold's exact qubit/bit-per-angle counts:** search-snippet level only this sprint.
* **Pamidimukkala *JCTC* 2024's** encoding details beyond "turn-based, higher DOF, up to 114
  qubits, HP alphabet": secondary.
* **LEGOLAS, UCBShift, ARTINA numbers:** search-snippet level. The ARTINA figures
  (1.44 Å median RMSD, 91.4 % assignment) come from a PMC/Science summary, not a fetched paper.
  **The parallel shift agent should re-fetch these before pricing the channel.**
* **TALOS-N's χ1 detail** (~50 % of residues, 89 % consistency): snippet-level.
* Ledger rows marked **[S]**: bibliographic details are reliable (most come from QTF's own
  reference list); method and result descriptions are not independently verified.

### My own inferences, clearly not the papers' claims
* **HYPOTHESIS:** the CDF decoder's monotonicity (§1.3) and its consequences — the sorted
  staircase, the ≤2π total-variation bound, and the resulting bias toward near-constant torsion
  profiles / regular secondary structure. *The monotonicity itself is arithmetic from their two
  printed equations and is not in serious doubt; the structural consequence is my inference and
  the paper neither states nor tests it.* This is testable in ~an hour against
  `github.com/cumbof/qtf` and is the highest-value follow-up in this document.
* **HYPOTHESIS:** P ≈ 2N + 6n > N, and the reading of the architecture as a classical
  reparameterisation rather than a compression (§1.7). *The formulas are theirs; the arithmetic
  and the interpretation are mine.*
* **HYPOTHESIS:** that a material fraction of QTF's apparent skill is the Ramachandran wells and
  the collapse bias rather than the search, by analogy with our own measurement that 88 % of
  what our sequence-conditioned library buys is generic Ramachandran and that a constant
  α-helix beats uniform random by 0.457 Å.
* **HYPOTHESIS:** that terminal exclusion makes their chignolin RMSD materially easier than our
  full-chain metric, by analogy with our Sprint 12 measurement that terminal dropout is
  0.40–0.50 Å cheaper than uniform dropout.

### What I did not establish
* I did not obtain the supplementary material, so I cannot write down their custom Hamiltonian.
* I did not read the QTF or PHEAT source, so `E_e2e` remains open.
* I did not download the Zenodo dataset (10.5281/zenodo.22088098). It contains the simulation
  and hardware data for both proteins and would permit an **independent recomputation of their
  energy–RMSD Spearman coefficients** — the numbers they never printed. That is a genuinely
  attractive, purely analytical, zero-native-leak opportunity for a future sprint, but it is
  not a Sprint 14 dependency and would breach the "no heavy compute" limit of this remit.
* I did not verify QFold's, Marchand's or Mato's primary texts this sprint.
* I ran no repository code and touched no benchmark artefact.

