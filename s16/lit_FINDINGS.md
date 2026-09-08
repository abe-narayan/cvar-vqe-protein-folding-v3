# Sprint 16 — LITERATURE FINDINGS: the continuous-torsion preprint, and the flagship's novelty boundary

Agent: `lit`. Date: 2026-09-06. Method: WebSearch / WebFetch, one PDF downloaded and extracted
locally with `pdftotext -layout`, one numerical check run in the scratchpad. **No repository
compute. No benchmark artefact touched. No repository file outside this one written.**

Extends `s15/LITERATURE.md` (961 lines, 67 papers) and `s14/lit_FINDINGS.md`. **It does not repeat
them.** Facts established there and not re-verified here are cited `[S14]` / `[S15]` and carry
their original tier. New work this sprint is marked with the tiers below.

---

## 0. TIERING, AND ONE HONESTY NOTE ABOUT THE TOOL

| tier | meaning |
|---|---|
| **[F]** | **FULL PRIMARY TEXT, read by me.** The PDF/HTML was downloaded to the scratchpad and I read the extracted text myself, including exhaustive greps. Only QTF qualifies this sprint. |
| **[A]** | **ABSTRACT / ARTICLE PAGE.** Obtained via WebFetch. |
| **[C]** | **CODE.** A source file or README in the authors' public repository, obtained via WebFetch. |
| **[S]** | **SNIPPET.** Search result only. **Never quotable as established fact.** |
| **[D]** | **MY OWN DERIVATION OR COMPUTATION**, run here and reproducible. |
| **[I]** | **MY INFERENCE**, not stated by any source. |

**The honesty note, and it matters for how much weight to put on the [A] and [C] rows.** `WebFetch`
does not hand me the page; it runs a small model over the page and hands me *that model's* rendering
of it, including its "verbatim" quotes. So an **[A]** or **[C]** row is *one summarising model
removed from the raw text*. It is much stronger than a snippet and much weaker than **[F]**.
Sprint 15's `[P15]` tier conflated these two things. **Where a claim in this document is
load-bearing, I say which tier carries it, and every "TAKEN" verdict below names the weakest link
in its own evidence.**

---

## PART A — THE SEPTEMBER 2026 CONTINUOUS-TORSION VQE PREPRINT, DISSECTED

### A.0 Identity, and what I did

> Cumbo F, Raubenolt B, Puram V, Katzenmeyer N, Joshi J, Blankenberg D.
> **"Logarithmic-scale variational quantum eigensolver for off-lattice protein structure prediction
> in continuous torsional angle space."** arXiv:2609.02113 [quant-ph; physics.chem-ph].
> Submitted 2 September 2026, 05:05:46 UTC. **v1 only — no v2 exists as of 2026-09-06.** 50 pages.
> CC BY 4.0. Cleveland Clinic (Computational Life Sciences; Lerner College of Medicine, CWRU).
> Method name **QTF (Quantum Torsion Folder)**. Code `github.com/cumbof/qtf` (MIT); support library
> PHEAT at `github.com/BlankenbergLab/pheat`. Data DOI 10.5281/zenodo.22088098.

**[F]** I downloaded the 9.18 MB PDF, extracted 1,801 lines with `pdftotext -layout`, and read it.
**[A]** I separately confirmed the arXiv abstract verbatim and the single-version status.
**[C]** I fetched the public repository — `README.md`, `qtf/scoring.py`, `qtf/metrics.py`,
`qtf/recipes.py`, `qtf/engines/qtf.py`, `qtf/assets/recipes/default.yaml`, and the PHEAT
`src/pheat/metrics.py` — **which Sprint 14 flagged as its single most important undone check.**

Sprint 14's dissection (`s14/lit_FINDINGS.md` §1.1–1.12) is **confirmed on every point I re-tested**
and is not repeated here. **This section reports only what is new, what I verified myself, and what
the code changed.** Three things are new: an exhaustive negative result on the RMSD protocol, the
resolution (and partial re-opening) of the end-to-end leak question from the source code, and a
barren-plateau claim that is asserted three times and never measured.

### A.1 What it actually does, and what is genuinely novel in it

QTF is **a classical optimiser searching a circuit-parameter vector through a fixed quantum-shaped
nonlinearity onto a torsion vector, scored classically.** That is not a hostile paraphrase; the
authors state the decisive half of it themselves.

Its claimed novelty is a **conjunction**: (quantum generator) ∧ (continuous, non-binned torsions) ∧
(all heavy atoms with side chains) ∧ (logarithmic qubit register). Each conjunct is prior art —
QFold (Casares et al.) binned torsions off-lattice in 2021–22; Marchand (2018) and Mato (2022) did
torsion-space quantum optimisation with real force fields; Pamidimukkala et al. (*JCTC* 20:10223,
2024, their own ref [19]) did high-DOF gate-based PSP — **and QTF cites none of QFold, Marchand or
Mato.** The combination does appear new [S14, re-confirmed].

**The genuinely novel object is the phase encoding**: using gauge-fixed relative phases of
statevector amplitudes as continuous real-valued torsions. I found nothing else like it in QPSP.
**[I]**, on a negative search.

**What it does not claim, and this is important for our positioning:** no quantum advantage, no
speedup, no accuracy claim against classical methods. It is a feasibility and qubit-resource paper.

### A.2 The logarithmic-scale encoding, exactly

**[F], verbatim structure from pp. 8–11.**

An n-qubit register with `n = max(2, ceil(log2 N_torsion))` spans `2^n` amplitudes carrying at most
`2^n − 1` independent *relative* phases after the unobservable global phase is fixed. QTF gauge-fixes
`c_0` (the all-zeros basis state) to phase 0 and reads

    φ̃_j(θ) = wrap_{[−π,π)} ( arg c_j(θ) − arg c_0(θ) )

as the torsion variables, mapping the `N−1` remaining relative phases onto the physical torsion
ranges. The zero-reference slot is assigned to the N-terminal φ, which is undefined anyway.
"For proteins with up to N ≤ 1,024 degrees of freedom, only n = 10 qubits are needed."

**Ansatz and parameter count, verified by my own arithmetic from their printed formulas [F]+[D]:**
Qiskit `EfficientSU2`, `d = ceil(N/n) + 2` repetitions, `P = 2n(d + 1)`, hence

    P = 2n·(ceil(N/n) + 3) ≈ 2N + 6n

| target | residues | N_torsion | n = ceil(log2 N) | d | **P** |
|---|---|---|---|---|---|
| chignolin 5AWL | 10 | 44 | **6** (Table 4 confirms) | 10 | **132** |
| Trp-cage 2JOF | 20 | 88 | 7 | 15 | **224** |

> **The classical optimiser's search space is ~3× larger than the torsion vector it produces.**
> The log-qubit claim is about register *width* only. Width is traded for depth (median transpiled
> depth 455, 204 CZ on ibm_cleveland for a 6-qubit, 10-residue problem) and for classical parameter
> count. **The paper never remarks on `P > N`.**

**Angular resolution:** none in the binned sense. In simulation the torsions are float64 reals read
off phases; on hardware the resolution is set by `1/n_shots` through the CDF, not by register width.
The authors concede this trades "the problem of qubit quantity for a problem of measurement
complexity" — a self-refutation of the resource claim on NISQ hardware, in their own words.

**Periodicity** is native and clean, and this is a real merit over binned encodings: the decoder
reads a phase, so there is no wraparound discontinuity to patch. **ω** is handled three ways:
fixed trans, affinely windowed into 170–190°, or free.

### A.3 What is quantum and what is classical — the ruthless split

**[F], and the authors say the decisive part themselves (p. 11, verbatim):**

> "the total scalar score is not decomposed into quantum observables and measured as the expectation
> value of a qubit Hamiltonian; it is computed classically from the reconstructed Cartesian structure
> after the circuit readout and minimized by the classical optimizer."

| component | quantum or classical |
|---|---|
| the map `R^P → C^{2^n} → [−π,π)^{N−1}` | **quantum** (a parameterised unitary; simulated exactly) |
| torsion → all-heavy-atom Cartesian (NERF) | classical |
| **the energy** | **classical** — no qubit Hamiltonian, no Pauli decomposition, no observable, no expectation value |
| the optimisation (COBYLA stage 1, SLSQP stages 2–3) | classical, derivative-free / finite-difference; **no parameter-shift, no adjoint, no natural gradient** |
| initialisation | classical **best-of-50 random scout** |
| snapshot mining, ranking, selection | classical |
| final minimisation | classical (GROMACS, amber99sb-ildn) |
| **hardware** | a **one-shot sampler**: 8,192 shots on a parameter-bound circuit, **no optimisation loop**, with a *different* (CDF) decoder |

**By Sprint 15's Rule 1 (`s15/LITERATURE.md` §4), QTF is NOT genuine VQE.** It is a
reparameterisation. Nothing in the reported result depends on the map being unitary or entangling.
**[F]** for every row; **[I]** for the reparameterisation reading — but it is a hard one to argue
against, because the authors supply the decisive sentence.

**The barren-plateau claim, which is new here and which Sprint 14 did not flag.** The abstract says
"We use an EfficientSU2 ansatz and multi-stage relaxation to mitigate barren plateaus", and the body
invokes the barren plateau twice more (ll. 54, 257, 645) as the motivation for the curriculum.
**I grepped the complete text for "gradient variance" and "variance of the gradient": zero hits.**
No gradient variance is ever measured, no scaling with `n` is shown, and the circuit is 6–7 qubits —
far outside any regime where a barren plateau has been demonstrated. **The BP claim is decorative.**
**[F]**, by exhaustive negative grep. This matters to us: `s15` V6 (Qiu et al. 2026) establishes that
tail objectives do not remove barren plateaus; a *curriculum* certainly does not, and no one has
claimed it does. **We must not make the symmetric mistake with CVaR.**

### A.4 How it extracts torsions, and continuous variables on a discrete register

Two decoders, and **they are not the same map** — a point the authors state plainly and which
invalidates any simulation↔hardware comparison in the paper:

1. **Simulator (produced every headline number):** gauge-fixed relative phases, above. Unconstrained.
2. **Hardware:** `C_i = Σ_{j≤i} P_j`, `θ_i = 2π C_i − π`, first `N` slots used.

**Sprint 14's structural defect stands and I re-derived it [F]+[D].** Because `P_j ≥ 0`, the `C_i`
are non-decreasing, so **the hardware-decoded torsion vector is a sorted staircase** with total
variation `2π(1 − P_0) ≤ 2π` across all `2^n − 1` steps — a mean consecutive step ≤ 8.2° over
chignolin's 44 used slots. The reachable set is an order-constrained sliver of the torus, and the
decoder is structurally biased toward *near-constant torsion profiles*, i.e. toward regular
secondary structure. **The paper never states this and never analyses it.**

**Why this is ours to say, not theirs.** Our own instrument measures that **a zero-information
constant α-helix (φ=−63°, ψ=−42°) emits 4.065 Å and beats uniform random by 0.457 Å**
[project record, `MEMORY: torsion-space-neither-energy-ranks`]. A decoder biased toward
near-constant profiles will look like it is working for exactly that reason, and **QTF's hardware arm
has no constant-torsion control.** **[I]**, and it is the single cheapest experiment that would settle
their hardware claim.

The authors themselves concede the confound in a different form: *"Because a matched noiseless
CDF-decoded baseline was not included, these contributions cannot be separated in the present
analysis."* **[F]**

### A.5 Energy evaluation and force field

Three interchangeable backends, 400 replicas each per protein (1,200 per protein):

1. **"Custom"** — a hand-weighted empirical scoring function, **not a force field**: clash, softened
   LJ with a *logarithmically capped* repulsive branch, Coulomb at ε≈4, geometric H-bonds, a
   sigmoid neighbour-counting burial proxy (not GB/PB), disulfide, **Gaussian Ramachandran wells on
   the α and β basins**, rotamer wells, aromatic stacking, Huber geometry-integrity penalties, and
   an **end-to-end constraint**. Parameters "partly inspired by" AMBER. Weights hand-set
   (50.0 / 15.0 / 8.0 → 2.5 / 1.5), never ablated, never cross-validated.
2. **Rosetta** REF15 + centroid via PyRosetta.
3. **OpenMM with AMBER ff14SB** — a genuine all-atom force field, **the same object we run**.
   No solvent model is stated.

Post-processing: GROMACS amber99sb-ildn steepest descent.

**A caution we should carry into our own Discussion.** Our record says generic Ramachandran accounts
for **88%** of what a sequence-conditioned torsion library buys. QTF's custom energy — the backend
that wins — contains explicit Gaussian wells on the α and β basins. **A large fraction of its
apparent skill may be that prior rather than the search.** **[I]**, with our own number behind it.

### A.6 Structural performance, targets, lengths — **and the RMSD convention**

Targets: **n = 2**, chignolin (5AWL, 10 res) and Trp-cage (2JOF, 20 res) — both among the most
optimised miniproteins in structural biology. 7,334,433 structures total.

| value | target | class | selector | tier |
|---|---|---|---|---|
| **0.623 Å** | chignolin | best retained snapshot (Rosetta pool, ~2.68 M) | **RMSD to native — ORACLE** | oracle |
| **1.199 Å** | chignolin | best final model (of 1,157) | **RMSD — ORACLE** | oracle |
| **2.90 Å** | chignolin | **median final model, custom** | energy | **predictive** |
| 4.20 / 3.84 Å | chignolin | median final, Rosetta / OpenMM | energy | predictive |
| **2.501 Å** | Trp-cage | best snapshot (~4.65 M) | **RMSD — ORACLE** | oracle |
| **5.39 Å** | Trp-cage | **median final model, custom** | energy | **predictive** |
| 6.59 / 6.21 Å | Trp-cage | median final, Rosetta / OpenMM | energy | predictive |
| 1.758 / 1.782 Å | chignolin | best of 300 hardware jobs | **RMSD — ORACLE**, warm-started | oracle |

Native-like hit rate, chignolin, fraction of replicas producing anything below 2.0 Å:
**snapshot mining 100 / 98 / 99.3 %** (custom / OpenMM / Rosetta) versus **final-model selection
6 / 4 / 0.3 %.** Trp-cage: **nothing below 2.0 Å anywhere, across 4.65 M structures.**

#### The RMSD convention — and this is a **new, exhaustive negative result**

**[F], by grep over the complete extracted text.** I searched for `terminal residues`,
`excludes the terminal`, `structural core`, `RMSD is/was computed`, `Kabsch`, `superimpos`,
`superpos`, `align`, `C-alpha RMSD`, `alpha-carbon`.

* **There is no RMSD definition in Materials and Methods. None.** No atom set, no superposition
  method, no residue range. The word `Kabsch` does not occur in the paper.
* The **abstract** says "Cα RMSD".
* The **only** protocol statement in 50 pages is a hardware aside on p. 31, verbatim:
  *"Because this end-to-end metric is defined by the two terminal Cα atoms, whereas the reported RMSD
  **excludes the terminal residues and evaluates the structural core**, a model can retain a
  moderately native-like core while exhibiting pronounced terminal displacement."*
* The paper repeatedly uses the undefined term **"effective RMSD"** (e.g. "the lowest-effective-RMSD
  snapshot from each replica") and **never defines it.**

**And the code does not agree with the paper. [C]:**
* `qtf/metrics.py` sets `PRIMARY_RMSD_ATOM_SET = "all-heavy"` and
  `DEFAULT_RMSD_ALIGNMENT_ATOM_SET = "same-as-rmsd"`, delegating to `pheat.metrics.structure_rmsd`.
* `pheat`'s `structure_rmsd` takes `atom_set` defaulting to `"all-heavy"`, offers `"ca"`, uses
  `kabsch_rmsd`, and **applies no residue-range filtering and no terminal exclusion whatsoever**.

> **Verdict: QTF's RMSD is under-determined in the paper and inconsistent with its own code.**
> The abstract says Cα; the shipped library's primary atom set is all-heavy; the paper says terminals
> are excluded; the library excludes nothing. **On a 10-mer, dropping both termini is 20% of the
> chain**, and our Sprint 12 measurement puts terminal dropout at **0.40–0.50 Å cheaper than uniform
> dropout**; BRIEF §3.3 records a **2.4 Å** full-chain-vs-`[1:-1]` gap on one structure in our set.
> **Every QTF number is therefore incomparable to ours until this is pinned, and it cannot be pinned
> from the paper.** This is a stronger statement than Sprint 14 could make, because Sprint 14 had
> only the aside; I have the aside, the exhaustive absence, *and* the contradicting source.

### A.7 The end-to-end bias: Sprint 14's open leak question, resolved — and a different one opened

Sprint 14 wrote: *"If [E_e2e's target] is set from the experimental structure, it is a native
information leak into the objective... the single most important thing a follow-up should check in
the github.com/cumbof/qtf source."* I checked. **[C]**

**Finding 1 — the README documents a purely length-aware target, with no native input:**

    target = 4.5 + 0.40 · max(0, N − 5)   (Å)
    slack  = 1.5 + 0.05 · N               (Å)
    E = λ · max(0, |d(Cα_first, Cα_last) − target| − slack)²

For N=10 this is a band of 4.5–8.5 Å (chignolin experimental 5.51 Å); for N=20, 8.0–13.0 Å
(Trp-cage experimental 11.25 Å). **Both natives fall inside a band set by sequence length alone.**
No PDB or experimental structure is read. **On this evidence the leak concern is discharged.**

**Finding 2 — but `qtf/scoring.py` at `main` does not implement that formula.** It implements:

    end_to_end_target = float(score_options.get("end_to_end_target", 5.5))
    constraint_strength = 5.0 if folder.current_stage == 3 else 50.0
    terms["end_to_end"] = constraint_strength * (dist_ends - end_to_end_target) ** 2

— **a hard-coded default of 5.5 Å with no length dependence, no slack band, and a plain harmonic.**
`qtf/assets/recipes/default.yaml` sets only *weights* (`end_to_end_weight` 8.0 → 1.5); a grep of the
recipe file returned **no occurrence of `end_to_end_target` and no per-target recipe for 5AWL or
2JOF.** `qtf/engines/qtf.py` computes the end-to-end distance as a reported metric but contains no
target formula.

**The paper reports chignolin's experimental end-to-end distance as 5.51 Å (Table 3, "Exp. E2E").**
The code's hard-coded default is **5.5 Å**.

> **[I], and it must be stated as an inference, not a finding.** I cannot establish which code
> version produced the runs — the repository has evolved, the README and `scoring.py` already
> disagree with each other and with the paper's stated λ values (paper 50.0 → 8.0 → 1.5; README
> 50.0/5.0; YAML 8.0/1.5). **Two readings are open and both are consequential:**
> (a) the runs used a length-aware target, in which case there is no leak and the prior merely
> happens to bracket both natives; or
> (b) the runs used the 5.5 Å default, in which case **chignolin was optimised toward its own
> experimental end-to-end distance to two decimal places**, and **Trp-cage was pulled toward 5.5 Å
> against a native 11.25 Å** — which would supply a mechanical, testable explanation for the paper's
> central asymmetry that is simpler than "Trp-cage has a more heterogeneous fold".
>
> **Either way, the paper never states the target distance anywhere.** For a term with weight 50.0
> in stage 1 that is a material omission, and it is the correct thing to raise if we ever review or
> cite this work. It is **not** something we may assert as a leak.

### A.8 Where it fails, and what it does not control for

**Stated by the authors, creditably [F]:** the statevector is unobservable on NISQ hardware and
"the core limitation lies precisely in its internal representation"; shot cost grows exponentially
with DOF, trading "qubit quantity for measurement complexity"; direct phase recovery needs
tomography; **"model selection remains the main bottleneck"**; no universal force field ranks
accuracy; physically suspect low-RMSD snapshots with "inferred terminal bond formation... when the
N- and C-termini were placed in abnormally close covalent proximity"; ~67 OpenMM replicas crashed
with numerical instability; and the hardware decoder confound quoted in A.4.

**Not controlled for, verified by exhaustive grep [F]:**

* **No classical control.** Grep for `baseline`, `classical control`, `random sampling`,
  `ablation`: the only hits are "a more favorable baseline for preserving chain dimensions"
  (a metric comparison) and the CDF-decoder admission. **Nowhere does anyone run COBYLA/SLSQP
  directly on the 44 torsions under the same energy and curriculum.** The causal contribution of
  the quantum component is unmeasured. They do not claim it, but the paper's framing invites the
  reader to.
* **No statistics.** `p-value` and `confidence interval`: **zero hits in 50 pages.** n = 2 targets.
* **No scout ablation.** Best-of-50 classical random search precedes every VQE run and is never
  removed, so the split between "the VQE" and "the classical pre-search" is unknown.
* **No random-sampling, constant-secondary-structure, or scrambled-sequence control.**
* **No gradient variance**, despite a barren-plateau claim in the abstract (A.3).
* **The custom energy's weights are hand-set and it is the backend that wins.** No ablation.
* **The hardware arm is warm-started from 1.199–2.270 Å classical answers** (median 1.784 Å) and
  **hardware execution monotonically degraded the best inputs**: 0/10 custom-energy starting models
  improved on either device.
* **P > N** is never remarked on.
* **The RMSD protocol** (A.6).

### A.9 QTF versus this project — the genuine differentiators, both directions

**What they have that we do not:**

| | |
|---|---|
| **Real quantum hardware execution** | 600 jobs across two IBM devices (Heron R2, Nighthawk R1), full transpilation/error-suppression reporting. We are simulator-only. |
| **All heavy atoms with side chains inside the generator** | Their torsion vector includes χ angles; our generative torsion arm is backbone φ/ψ with AMBER applied downstream. |
| **Scale of sampling per target** | 7.3 M structures across 2 targets; 1,200 replicas each. |
| **A continuous, un-binned angular representation with native periodicity** | Our library arm is discrete at k=4 (though `s15/distgeo` and `align_lib.fit` are continuous). |
| **A published, installable, MIT-licensed package with recipes and a support library** | Ours is a research repository. |
| **A log-qubit register** | 6 qubits for a 10-mer. Ours is `(n−2)·log2(k)`. |

**What we have that they do not — and every one of these is a control or a discipline, not a
capability:**

| | |
|---|---|
| **A defined, frozen RMSD protocol** | Cα, all residues, no trimming, proper rotations, model 1, chain breaks scored. QTF has none in Methods and its code contradicts its abstract (A.6). |
| **126 cluster-disjoint targets, 5 pinned folds, a disjoint dev24, a sealed 60-target benchmark** | QTF has 2 hand-picked miniproteins. |
| **Paired bootstrap CIs, per-fold detail, W/L, null-calibrated concentration checks** | QTF has zero CIs and zero p-values. |
| **A genuine qubit Hamiltonian** | Our objective *is* measured as an observable; QTF's is explicitly not. **By our own Rule 1, QTF is not genuine VQE and we are.** |
| **Genuine in-loop CVaR** | QTF has none — "CVaR" does not appear in the paper. |
| **The classical control at matched budget** | best-of-N from the untrained circuit, matched diversity, matched cost in evaluations *and* wall-clock, Boltzmann reweighting of the same samples, annealing, greedy, exact enumeration. QTF runs none. |
| **A certified global optimum over an enumerated space** | 19 targets, 1.28e7 exactly labelled structures. |
| **Leakage discipline** | Every ORACLE arm labelled in every table row. QTF's abstract leads with an argmin-against-the-native. |
| **A stated, audited objective** | Our AMBER/Legacy terms are identifiable and separately measurable; QTF's winning backend is a hand-weighted function whose end-to-end target is not stated anywhere (A.7). |

**The one-line positioning.** *QTF is the field's most ambitious continuous-torsion quantum
generator and has no controls; we have every control and a smaller generator. The two papers are
not competitors — they are the two halves of one missing experiment.* Comparability verdict is
unchanged from Sprint 15: **[I-oracle][I-metric][I-scale]** on the headline; only the 2.90 / 5.39 Å
medians are worth comparing, and even those are terminal-trimmed on n = 2.

---

## PART B — THE FLAGSHIP'S NOVELTY BOUNDARY

The hypothesis under test:

> *native-free estimation of structural error **DIRECTION**, combined with native-free estimation of
> the structurally **QUIET** (low-response) Jacobian subspace, used to steer generation.*

I searched each ingredient and the combination. **Two of the three ingredients are taken, one of
them decisively. The combination appears open. And the fusion law is taken outright.**

### B.1 Ingredient A — native-free estimation of structural error DIRECTION. **TAKEN. Twice.**

#### B.1.1 ATOMRefine — the decisive one

> Wu T, Guo Z, Cheng J. **"Atomic protein structure refinement using all-atom graph representations
> and SE(3)-equivariant graph transformer."** *Bioinformatics* **39**(5):btad298 (2023).
> doi:10.1093/bioinformatics/btad298. Code: `github.com/BioinfoMachineLearning/ATOMRefine`.

**[A]**, PMC article page fetched. Quoted from that page:

> *"A SE(3) transformer is used to predict the **coordinate shifts between the initial model and the
> native structure**... the final refined structure is simply equal to the initial model plus the
> predicted shift."*

Native information is used **only in training**; at inference the network sees the model and emits a
per-atom 3-D displacement vector toward the native. **That is native-free estimation of the
structural error vector — direction and magnitude — applied directly to steer the structure.**

**Our claim "native-free estimation of structural error direction" cannot be presented as new.**
It is published, in the strongest possible form (a full 3-D displacement field), three years ago,
in *Bioinformatics*.

#### B.1.2 DeepAccNet — the signed-error-as-restraint form

> Hiranuma N, Park H, Baek M, Anishchenko I, Dauparas J, Baker D. **"Improved protein structure
> refinement guided by deep learning based accuracy estimation."** *Nat Commun* **12**:1340 (2021).

**[A]**, PMC full text fetched. It predicts *"per-residue-pair distributions of **signed** Cβ–Cβ
distance error from the corresponding native structures (referred to as 'estograms')"*, with
*"Red and blue indicat[ing] that the pair of residues are too far apart and too close"* — and then
*"estograms were converted to residue–residue interaction potentials... added to the Rosetta energy
function as restraints."* Trained on ~1 M alternative minima over 7,510 proteins.

So the *sign* of every pairwise distance error is estimated native-free and used to steer
generation. **This is our S15 K2/K9 channel, in the pair-distance basis, published in 2021.**

#### B.1.3 What survives for us, stated narrowly

1. **The basis.** Ours is a **torsion-space** error direction (`θ_fit − θ_pool`); theirs are
   Cartesian displacements (ATOMRefine) and pair-distance signs (DeepAccNet). I found nothing
   estimating the error direction *in torsion coordinates*.
2. **The mechanism.** Ours is **native-free by disagreement between two channels**, not by
   supervised regression on native offsets. Both prior works are supervised on the native.
   For a 9–16-mer regime with ~126 targets, a disagreement surrogate is arguably the only
   available route, and that is worth saying.
3. **The use.** Both prior works use the estimate to **cancel** the error. Our flagship uses it to
   decide **which subspace to rotate the residual error into**, at *fixed magnitude*. That is a
   different intervention (§B.4).
4. **The measurement.** `|cos| 0.390` against a 0.157 null, alignment corr `+0.499 [+0.33, +0.65]`.
   **Neither prior work reports the accuracy of its error-direction estimate at all** — DeepAccNet
   reports cross-entropy, not a correlation; ATOMRefine reports no directional analysis. Reporting
   the cosine of a native-free error-direction estimate against the true error, with a null,
   appears unclaimed. **[I]**, on a negative search.

#### B.1.4 **The calibration warning, and it is the most useful thing in this section**

ATOMRefine — a supervised SE(3)-equivariant transformer trained on the exact quantity, from the
Cheng lab, published in *Bioinformatics* — improves **GDT-HA from 69.84 to 70.04 on 193 AlphaFoldDB
targets, and 63.91 to 64.42 on 69 CASP14 targets.** **[A]**

> **The published state of the art at this exact ingredient buys about 0.2–0.5 GDT-HA points.**
> Our best native-free alignment arm of 22 tried bought **−0.007 Å [−0.074, +0.055]**. Those two
> facts belong in the same sentence. **The flagship's prior should be set by ATOMRefine's margin,
> not by the 1.855 Å ORACLE.** A null result here is not a failure of our implementation; it is
> the field's result too, and we would be the first to *report* it as a bounded quantity with a CI
> rather than as a small positive delta.

### B.2 Ingredient B — the QUIET (low-response) Jacobian subspace. **TAKEN IN SUBSTANCE.**

The object is the SVD of `J = ∂(Cartesian)/∂(torsions)`, and its low-singular-value directions.

| # | prior art | what it establishes | tier |
|---|---|---|---|
| **J1** | **Manipulability ellipsoid** (Yoshikawa 1985) and standard robot kinematics | The SVD of the manipulator Jacobian is *the* textbook tool; the low-singular-value directions are the ones "difficult to move in", and the inverse condition number is the standard dexterity index. **Our spectrum, participation ratio and "15.7× more RMS angle in the quiet half" are a manipulability analysis of a polypeptide.** | **[S]** (Wikipedia/Modern Robotics level; the underlying result is textbook and uncontroversial) |
| **J2** | **Internal-coordinate NMA**: López-Blanco & Chacón, `iMod` (*Bioinformatics* 27:2843, 2011) and `iMODS` (*NAR* 42:W271, 2014); Noguti & Gō lineage | NMA in torsion space builds the torsion→Cartesian Jacobian explicitly to project mass and stiffness. Low-frequency modes are the standard steering directions for flexible fitting. | **[S]** |
| **J3** | **Protein inverse kinematics / self-motion manifolds**: Coutsias & Seok, "Inverse Kinematics in Biology: The Protein Loop Closure Problem"; Liu & Latombe, "On the Structure of the Inverse Kinematics Map of a Fragment of Protein Backbone"; *R Soc Open Sci* 11:240873 (2024) on the topology of protein self-motions | Null spaces of the torsion→endpoint map, their dimension and topology, are a mature literature. | **[S]** |
| **J4** | **Nullspace / kino-geometric sampling**: Fonseca, Budday, van den Bedem et al., "Nullspace Sampling with Holonomic Constraints Reveals Molecular Mechanisms of Protein Gαs", *PLoS Comput Biol* 11:e1004361 (2015); "Kinematic Flexibility Analysis" (arXiv:1802.08683) | **Sampling explicitly in the null space of a constraint Jacobian** to move a protein while preserving specified geometry. This is *steering in the quiet subspace*, named and implemented. | **[A]**/[S] |
| **J5** | **"Wriggling" and concerted rotation**: "Proteins Wriggle" (arXiv:cond-mat/0108218); biased Gaussian steps in torsional space (arXiv:cond-mat/0103580); Dodd/Theodorou concerted rotation; Ulmschneider & Jorgensen | Verbatim from the search return: *"small angular or torsional changes may result in large displacements far from the site of modification"*, and proteins *"alter their shapes by means of local motions that minimize the displacement of distant atoms."* **This is exactly the loud/quiet distinction, stated as the design principle of a Monte Carlo move set.** | **[S]** |

**Verdict.** *"Some torsion directions move the structure much less than others, and you should
prefer them"* is **established, in at least three separate literatures** — robot manipulability,
protein IK/null-space sampling, and polymer Monte Carlo move design. Our AL1.1 measurement
(participation ratio 3.47 of ~25.9; 90% of `trace(JᵀJ)` in 4.5 directions; 15.7× angular tolerance
in the quiet half) is **a quantification of a known phenomenon on a new object**, not a discovery.

**What plausibly survives, narrowly:**
* **The superposed Jacobian.** Ours is `∂(Kabsch-superposed CA)/∂(φ,ψ)` — the Jacobian *after*
  optimal superposition. The IK/NMA literature works in a body frame or with explicit constraints;
  the non-Eckart-frame issue is known (*J Math Chem* 2012, "Normal mode analysis... on a non-Eckart
  body-frame: an application to protein torsion dynamics" **[S]**) but I found no one computing and
  spectrally analysing the superposition-projected torsion Jacobian for a free peptide.
* **The count.** "Exactly five null directions, `2n−4` parameters onto a `2n−5`-dimensional shape
  space, the fifth a distributed whole-chain crank with `Σδφ ≈ −Σδψ`" is a clean structural result
  I did not find stated. **[I]**, on a negative search, and it is elementary enough that it may
  exist in the kinematics literature — S15 flagged one unrun robotics search pass for the locality
  theorem and **that pass is still not properly run for the null-direction count either.**
* **Nothing else.** The spectrum, the participation ratio and the tolerance ratio should be
  presented as measurements, with J1–J5 cited.

**And the native-free half of ingredient B is weaker than it reads.** Estimating `J(native)`'s quiet
subspace by `J(emitted)`'s (cos² 0.827 vs a 0.499 null, spectra correlated +0.980) is **the
Gauss–Newton assumption** — that the Jacobian at the current iterate approximates the Jacobian at
the solution. Every Gauss–Newton, Levenberg–Marquardt and NMA-based flexible-fitting method in
existence already assumes it. **Our contribution is measuring how well it holds here, which is
genuinely worth reporting, but it must not be framed as "we found that the quiet subspace is
identifiable native-free."** **[I]**, and I am confident in it.

### B.3 The MECHANISM — "error in quiet directions costs less" is the discrete Picard condition

AL1.2's ORACLE result (rotate the error out of the loud half at fixed magnitude → 1.855 Å,
−1.822 [−2.037, −1.613], W/L 121/5) is a statement that `‖J e‖` depends on the alignment of `e` with
`J`'s singular vectors. **That is the content of classical regularisation theory for discrete
ill-posed problems:**

> Hansen PC. **"The discrete Picard condition for discrete ill-posed problems."** *BIT* 30:658 (1990);
> and **"The truncated SVD as a method for regularization."** *BIT* 27:534 (1987).

**[S]**, and the underlying mathematics is textbook. The discrete Picard condition is precisely the
statement about *how the error's spectral coefficients decay relative to the singular values*, and
TSVD/Tikhonov are the standard responses. Our AL1.3 result — **isotropic Tikhonov beats both
geometry-aware metrics at every λ**, +0.124 [+0.025, +0.228] at the selected λ — is **the expected
outcome in that theory**, and we should present it that way rather than as a surprise.

> **Actionable.** Framing the flagship as *spectral regularisation of an ill-posed inverse problem
> under a structured, correlated noise model* connects it to sixty years of numerical analysis,
> gives us the right null (isotropic Tikhonov, already measured), and makes the honest claim —
> *we are measuring whether a native-free estimate of the noise's spectral orientation beats the
> isotropic default* — which is both new and modest.

### B.4 The COMBINATION — **appears OPEN**

I searched for work combining an estimated model-error direction with a Jacobian sensitivity /
null-space analysis to choose a search or steering direction, in structure prediction, refinement,
restraint-based modelling, robotics, and inverse problems. **I found nothing that does both.**
Robotics does task-priority null-space projection (secondary objective in the null space of the
primary), but the *primary* there is a commanded task, not an estimated error; DeepAccNet and
ATOMRefine estimate the error and apply it isotropically, with no sensitivity weighting; KGS
samples the null space with no error estimate at all.

**Verdict: the COMBINATION is the flagship's only defensible novelty claim. [I]**, on a negative
search across four literatures, which is the weakest kind of evidence — but the components are so
well-known that anyone who had combined them would have said so.

**And it must be claimed at the right size.** The honest statement is:

> *Both an error-direction estimate and a sensitivity-subspace estimate are separately available
> without the native. Whether combining them steers generation better than isotropic regularisation
> is, as far as we can find, untested — and the published precedent for each ingredient alone
> (ATOMRefine, +0.2 GDT-HA) says the effect will be small.*

### B.5 The FUSION LAW — **TAKEN. This is the clearest claim-taker in the document.**

Our record calls `d_avg ≈ √(r² − (s/2)²)` a **parameter-free law**, "elementary geometry, exact when
the native is equidistant from both and the three points are coplanar", predicting to 0.135–0.162 Å
on 24 and then 126 targets.

**It is the Krogh–Vedelsby ambiguity decomposition (1995), specialised to two members.** For
members `f_i`, mean `f̄`, target `y`, under squared loss:

    (f̄ − y)² = (1/B) Σ_b (f_b − y)² − (1/B) Σ_b (f_b − f̄)²
             = average error − ambiguity (diversity)

**[S]** for the citation trail (Krogh & Vedelsby, NIPS 1995; Ueda & Nakano 1996 bias–variance–
covariance; Brown et al., *JMLR* 6:1621, 2005), which is standard ensemble-learning material.
For `B = 2`, `(1/2)Σ(f_b − f̄)² = (s/2)²` exactly, giving our law.

**I verified the algebra numerically myself. [D]**, `scratchpad`, five random 40-point trials:

| | |
|---|---|
| `d_avg²` measured | equals `(r₁² + r₂²)/2 − (s/2)²` **to machine precision, every trial** |
| coplanarity | **not required** |
| equidistance | **not required** |
| using `r = (r₁+r₂)/2` instead of the quadratic mean | under-predicts by **0.07–0.15 Å** at our error scale |

> **Two consequences, and the second is a live correction to the sprint's own record.**
>
> 1. **The law is not ours and is not new.** It is an exact identity in any inner-product space,
>    known in ensemble learning for thirty-one years, and it has a structural-biology twin (the
>    standard relation between ensemble-average pairwise RMSD and RMSF about the mean structure —
>    *Biophys J*, "Determination of ensemble-average pairwise RMSD from experimental B-factors",
>    2010 **[S]**). **It must be cited, not claimed.** What survives is *empirical*: that it holds
>    to 0.162 Å on 126 real targets through a **Kabsch superposition**, which is a nonlinear
>    operation the identity does not cover, and that **`s` is native-free**, so it can be used as
>    an *a priori* screen for whether fusion is worth running. That last use is genuinely
>    practical and I found no one doing it. **[I]**
> 2. **Most of the "prediction error" may be an arithmetic artefact, not physics.** Our record
>    computes `r` as the *mean* of the two channels' RMSDs (3.270 and 3.469 → 3.370). The identity
>    needs the *quadratic* mean. With `r₁=3.270, r₂=3.469, s=2.980`: quadratic-mean prediction
>    **2.981 Å**, arithmetic-mean prediction **2.970 Å**, measured **3.105 Å**. Here the two agree
>    closely because `r₁ ≈ r₂`, so this particular row is not affected much — **but the substitution
>    is wrong in general and will bite on any channel pair with unequal accuracy.** **[D]**
>    **Recommendation: recompute the fusion-law residual with the quadratic mean before it appears
>    in any write-up, and report the Kabsch-nonlinearity as the residual's stated cause.**

### B.6 Realizability of the predicted distance matrix (K9) — **largely OPEN, weakly adjacent**

I searched EDM-cone arguments, projection onto the EDM cone, triangle-inequality/metric violation in
predicted distance maps, and distogram self-consistency.

* EDM/cone formulations of the *conformation* problem exist (e.g. *J Glob Optim* 2019,
  "A Euclidean distance matrix model for protein molecular conformation" — rank-constrained least
  squares with a cone constraint) **[S]**, but they are *solvers*, not diagnostics of a predictor's
  error.
* Distance-prediction assessment work (`DISTEVAL`, *BMC Bioinformatics* 2021; "Toward the assessment
  of predicted inter-residue distance", *Bioinformatics* 38:962) evaluates accuracy, not
  realizability **[S]**.
* The closest structural statement I found is that **"most false contact predictions are shifted by
  one or two residues with respect to true native contacts"** (contact-prediction literature,
  BMC Struct Biol 8:36, 2008) **[S]** — a coherent-error observation, but about registration shift,
  not about the errors describing a realizable alternative structure.

> **Verdict: "the predictor's errors describe a consistent WRONG STRUCTURE, and are 2.4× closer to
> realizable than matched-magnitude noise (ratio 0.413, paired −1.391 [−1.623, −1.184])" appears
> unclaimed. [I]** It is, with the CVaR defect triple and the untrained-circuit control, one of the
> three most defensible things the programme holds. **Two cautions:** (i) the *general* fact that
> structure predictors make systematic rather than i.i.d. errors is not new and must not be claimed
> — it is the premise of every consensus/QA method and of DeepAccNet's estogram; (ii) the specific,
> quantified, null-controlled realizability measurement is what is ours.

### B.7 Error SHAPE versus error MAGNITUDE — **partially adjacent; scope it**

The K12 result (outlier-shaped error at 4.12 Å RMS emits 1.993 Å; i.i.d. at 3.00 Å emits 2.561 Å)
is a shape-beats-magnitude statement. The nearest literature is the NMR restraint-error work —
"depending on the error bounds used, these distance restraints can seriously distort the structure,
leading to deviations even in rigid regions", and the error-distribution-derived restraint potentials
of Wang & Nilges (*J Biomol NMR* 2006) **[S]** — which establishes that the *shape of the restraint
error distribution* matters more than its width. **That is the same claim, in the same domain,
twenty years earlier, qualitatively.** Ours is a controlled dose–response surface with a matched
null. **SURVIVES AS A MEASUREMENT; the qualitative claim must be cited to the NMR restraint
literature.** **[I]**, on **[S]** evidence — this row deserves a primary fetch before it is written up.

### B.8 Novelty boundary — summary table

| ingredient | verdict | authority | what survives |
|---|---|---|---|
| native-free error **DIRECTION** estimation | **TAKEN** | **ATOMRefine, *Bioinformatics* 39:btad298 (2023)**; **DeepAccNet, *Nat Commun* 12:1340 (2021)** | the torsion-space basis; the disagreement (unsupervised) surrogate; reporting `\|cos\|` against a null |
| native-free **QUIET SUBSPACE** estimation | **TAKEN IN SUBSTANCE** | manipulability ellipsoids; iMod/iMODS; protein IK & self-motion manifolds; **KGS nullspace sampling, *PLoS CB* 11:e1004361**; "Proteins Wriggle" / concerted rotation | the *superposed* Jacobian; the five-null-direction count; the quantified spectrum. The native-free part is **the Gauss–Newton assumption** and must be labelled as such |
| the **MECHANISM** (quiet-direction error is cheap) | **TAKEN** (textbook) | Hansen, discrete Picard condition; TSVD/Tikhonov | the magnitude in this domain; the ORACLE ceiling |
| **the COMBINATION**, used to steer | **OPEN** | negative search across 4 literatures | this is the flagship's novelty, and it should be claimed alone |
| the **FUSION LAW** `√(r²−(s/2)²)` | **TAKEN** | **Krogh & Vedelsby (1995) ambiguity decomposition**; Brown et al. *JMLR* 6:1621 (2005); ensemble-RMSD/RMSF identity | that it holds through Kabsch on 126 real targets; that `s` is native-free, hence an *a priori* fusion screen |
| **realizability** of the error (K9) | **OPEN** | negative search | the quantified, null-controlled measurement — not the general "errors are systematic" |
| error **SHAPE > MAGNITUDE** | **PARTIALLY TAKEN** | NMR restraint-error literature | the controlled dose–response surface |

---

## PART C — WHAT MUST BE REMOVED OR RE-SCOPED

Sprint 15 removed four claims. **These are the new ones. The first three are the important ones.**

### C.1 REMOVE: "a parameter-free fusion law"

**It is the Krogh–Vedelsby ambiguity decomposition.** See §B.5. Replace with: *"the two-member case
of the standard ensemble ambiguity decomposition, verified to hold through Kabsch superposition on
126 targets, and used as an* a priori *native-free screen for whether fusion is worth running."*
**And fix the arithmetic mean → quadratic mean before quoting the residual.**

### C.2 REMOVE: "native-free estimation of structural error direction is new"

**ATOMRefine predicts the coordinate shift to the native and applies it.** See §B.1. Replace with
the torsion-basis / disagreement-surrogate / fixed-magnitude-rotation framing, and **calibrate the
expected effect against ATOMRefine's +0.2 GDT-HA.**

### C.3 REINSTATE THE THREAT: the AMBER ff14SB / GB-Neck2 peptide claim

Sprint 15 §2.4 downgraded the threat to the AMBER anti-ranking claim because Maffucci & Contini
(2016) tests ff96 and the ff99SB series, **not** ff14SB — and flagged the *PCCP* 2018 study as
"the highest-value unfetched source in this document". **I did not obtain its full text (RSC 403,
PubMed cookie wall), but I obtained enough to reverse the downgrade. [S], and I flag the tier
loudly:**

> Shao Q, Zhu W. **"Assessing AMBER force fields for protein folding in an implicit solvent."**
> *Phys. Chem. Chem. Phys.* **20**(10):7206–7216 (2018). doi:10.1039/C7CP08010G. PMID 29480910.

The search return states that the study runs enhanced-sampling MD over **six AMBER force fields —
FF99SBildn, FF99SBnmr, FF12SB, FF14ipq, **FF14SB**, and **FF14SBonlysc** — each coupled with the
**GB-Neck2** model**, on two helical and two β-sheet peptides, and concludes that *"a combination of
the force fields and GB-Neck2 implicit model able to describe all aspects of the folding transitions
towards the native structures of all the considered peptides was not identified"*, with
FF14SBonlysc/GB-Neck2 "a reasonably balanced combination".

> **This is our exact pair — ff14SB + GB-Neck2 — tested on peptides and found not to describe
> folding to the native.** Sprint 15's stated condition for the downgrade ("it does **not** test
> ff14SB") no longer holds.
>
> **What still separates us, and it is real:** their measurement is *folding thermodynamics from
> enhanced-sampling MD*, not **conformer ranking of a fixed pool**, and certainly not a **certified
> global optimum over an enumerated space with a paired CI and a filtering/ordering decomposition**.
> **Re-scope, do not withdraw:** cite Shao & Zhu as prior context establishing that ff14SB/GB-Neck2
> does not place peptide natives at the global free-energy minimum, then state exactly what
> enumeration + certification + CI adds.
>
> **ACTION: this source is still only [S]. Fetch the full text before any paper draft. It is now
> the highest-value unfetched source in the programme, promoted from Sprint 15's list — not
> demoted.**

### C.4 RE-SCOPE: "searching harder on a bad objective hurts" is a named phenomenon in two other fields

Our standing ledger entry reads *"the budget trap is a property of BAD OBJECTIVES: searching harder
hurts on a bad one, is neutral on a mediocre one, and HELPS on a good one."* **This is discovered
work elsewhere, twice, with mechanisms:**

* **Game-tree search pathology.** Nau (1979–1983): *"minimaxing amplifies the error of the heuristic
  evaluations and consequently deeper searches produce worse evaluations"*, with explicit conditions
  (branching factor, independence of sibling values) under which pathology occurs — and later work
  (Luštrek & Gams; Zuckerman & Nau, "Avoiding game-tree pathology in 2-player adversarial search")
  identifying when it does not. **[S]**
* **Proxy-objective over-optimisation.** Gao L, Schulman J, Hilton J. **"Scaling Laws for Reward
  Model Overoptimization."** ICML 2023 (PMLR 202), arXiv:2210.10760. A gold objective is optimised
  through a proxy; **for best-of-n sampling the gold score follows `R(d) = d(α − βd)` with
  `d = √KL` — a downward parabola** — and the coefficients scale smoothly with proxy quality
  (reward-model size). **[S]**, but the functional form was returned consistently.
* **The optimizer's curse.** Smith JE, Winkler RL, *Management Science* 52:311 (2006). **[S]**

> **Consequence for Sprint 16.** The phase boundary we intend to *measure* is a known phenomenon
> whose *shape* has been fitted in another field, using **best-of-n as the optimiser** — the same
> control family our BRIEF §3.5 mandates. **Cite Nau and Gao et al.; do not present the phenomenon
> as a discovery.** What is genuinely open is doing it for a *variational quantum sampler* on a
> *molecular* objective with `blend_objective`'s rank-space quality knob — i.e. the same
> re-scoping that saved claim 3 in Sprint 15.

### C.5 RE-SCOPE: the quiet-subspace measurements

Present AL1.1's spectrum, participation ratio and 15.7× tolerance as **a manipulability analysis of
a peptide**, citing J1–J5 (§B.2). Do not present "some torsion directions are quiet" as a finding.
Do not present the native-free quiet-subspace estimate as anything other than **a measurement of
how well the Gauss–Newton assumption holds** (cos² 0.827 vs 0.499).

### C.6 Sprint 15's four removals stand, reinforced

Certified-optimum-worse-than-random (Roget et al. 2606.21241); in-loop CVaR for peptide folding
(Uttarkar et al., *PLoS One* 21(2):e0342012); tail-restricted discrimination as a framework
(Truchon & Bayly, BEDROC); decoy-bank ranking non-transfer (Handl, Knowles & Lovell 2009). Nothing
this sprint weakens any of them.

---

## PART D — THE TWO QUESTIONS THE ARCHITECTURE NEEDS ANSWERED

### D.1 Does anyone run a control equivalent to best-of-N from an untrained circuit? **STILL NO — but the search found two near-misses Sprint 15 did not have, and they must be cited.**

I attacked this claim under six formulations: untrained-circuit baselines, random-parameter circuit
baselines, best-of-N sampling baselines for VQAs, trivial baselines in quantum optimisation
benchmarking, "does optimisation help", and matched-shot-budget comparisons.

| # | near-miss | what it actually does | why it is not our control | tier |
|---|---|---|---|---|
| **N-1** | **Chernyavskiy AY, Bantysh BI, Bogdanov YI, "Entropic property of randomized QAOA circuits", arXiv:2308.01807 (2023)** | Studies **QAOA circuits with random (untrained) angles** as an object; proves higher energy-distribution entropy than uniform bitstring sampling and shows *"the probability to obtain the global optima... appears to be higher on average than for random sampling."* | **It compares the untrained circuit only to uniform random bitstrings, never to the optimised circuit at matched shots.** It is the first paper I have found to treat the untrained-circuit distribution as a baseline object at all, and it establishes that the untrained circuit is *already* better than uniform — which is the premise our control rests on. **We should cite it as the justification for why uniform random is the wrong control.** | **[A]**, abstract verbatim |
| **N-2** | **Scriva G, Astrakhantsev N, Pilati S, Mazzola G, "Challenges of variational quantum optimization with measurement shot noise", *Phys Rev A* 109:032408 (2024), arXiv:2308.00044** | Scales the **shot cost** of VQE/QAOA against brute force. Concludes: *"standard VQE scales comparably to brute-force search"*, and *"when the parameters are optimized from random guesses, also the scaling of QAOA implies problematically long absolute runtimes"*, becoming practical only with physics-informed initialisation. | A **cost-scaling** argument, not a paired matched-budget comparison of an optimised circuit against its own untrained samples. But it is the strongest published statement that **optimising from a random start may not repay its shot cost** — adjacent to our result and it must be cited. | **[A]** |
| **N-3** | Boulebnane, Lucas, Meyder, Adaszewski, Montanaro, *npj QI* 9:70 (2023) | QAOA matched by **uniform random sampling** "up to a small overhead". | Sprint 15's finding; weaker control (does not hold the ansatz's induced distribution fixed), and a *parity* result, not a loss. | [S15] |
| **N-4** | Kesiku & Garcia-Zapirain, "How Quantum Circuits Actually Learn: A Causal Identification of Genuine Quantum Contributions", arXiv:2603.16321 (Mar 2026) | **Counterfactual causal mediation** decomposing architectural performance differences into direct effects and effects mediated by entanglement entropy / quantum mutual information; 43 configurations, 5 topologies. **"Direct architectural contributions exceed quantum-mediated effects at a 13.1:1 ratio"**, and current VQCs "operate substantially below their quantum potential." | A **different** causal methodology, on ML benchmark datasets, comparing *architectures* — not a same-circuit untrained baseline. **But it is the field beginning to ask our question**, and it is a competing method for causal attribution that a referee may raise. **Fetch it in full before drafting.** | **[A]** |
| **N-5** | Benchmarking practice generally (arXiv:2408.03073 and the QMI 2026 version; the quantum-optimisation benchmarking working groups) | Standard baselines are **uniform random sampling at matched shot count** and **greedy**. | The field's default is exactly the weaker control. | **[A]**/[S] |

> **VERDICT: the untrained-circuit best-of-N control SURVIVES as our lead novelty claim.** But the
> claim must now be stated with N-1 and N-2 named, and phrased precisely: *a paired comparison of
> the optimised variational state against best-of-the-same-number-of-shots drawn from the **same
> ansatz at its initial parameters**, at matched shot budget, reported with a CI.* Chernyavskiy et
> al. supply exactly the reason it is the right control (the untrained circuit already beats
> uniform); Scriva et al. supply the cost framing. **Citing both makes the claim stronger, not
> weaker** — it stops reading as "we found a control nobody thought of" and starts reading as
> "the field has been circling this and here is the measurement."
>
> **Standing caveat, unchanged:** this is an absence claim, and absence searches are weaker than
> presence searches.

### D.2 Is there published work on the objective-quality condition under which a variational sampler beats classical search? **NOT FOR QUANTUM. YES, TWICE, ELSEWHERE — AND WE MUST CITE IT.**

* **No quantum instance found.** The nearest is instance-level ("when does VQE beat greedy on *this*
  instance"), e.g. the benchmarking line in N-5, which correlates algorithms instance-by-instance
  but does not vary *objective fidelity* as a controlled axis. **[S]/[A]**
* **Game-tree search:** Nau's minimax pathology gives the *condition* on the evaluation function
  under which more search degrades decisions. **[S]**
* **RLHF:** Gao, Schulman & Hilton give the *functional form* — `R(d) = d(α − βd)` for best-of-n —
  and show its coefficients scale smoothly with proxy quality. **[S]**

> **VERDICT: the phase boundary exists in the literature, in two fields, with mechanisms and one
> fitted functional form — but never for a variational quantum sampler and never on a molecular
> objective.** Sprint 16 should **pre-register a comparison of its measured boundary against Gao et
> al.'s downward-parabola form**, since `s14/vqe_lib.blend_objective` varies objective quality in
> rank space with no scale or tail-shape confound and is therefore a *cleaner* instrument for the
> same question than reward-model size. **That reframes the experiment from "we discover a phase
> boundary" to "we test a published functional form in a new regime", which is a stronger paper and
> an honest one.**

---

## PART E — THE 2024–2026 SWEEP, TIERED

Only rows that are **new since `s15/LITERATURE.md`** or that **change a scoping decision**. The 64
rows of Sprint 15's §4 stand and are not repeated.

### E.1 Quantum PSP and VQAs

| # | source | what it adds | tier |
|---|---|---|---|
| E1 | **Cumbo et al., arXiv:2609.02113 (QTF)** | Part A. **v1 only, no v2.** Code fetched; RMSD convention undetermined and code-contradicted; E2E target unstated in the paper and hard-coded 5.5 Å at repo `main`; barren-plateau claim never measured. | **[F]** + **[C]** |
| E2 | Chernyavskiy, Bantysh, Bogdanov, arXiv:2308.01807 (2023) | Untrained-circuit sampling as a baseline object; random-angle QAOA beats uniform bitstrings. **Cite as the reason uniform random is the wrong control.** | **[A]** |
| E3 | Scriva, Astrakhantsev, Pilati, Mazzola, *Phys Rev A* 109:032408 (2024) | Shot-cost scaling: standard VQE ≈ brute force; optimisation from random starts is cost-problematic. | **[A]** |
| E4 | Kesiku & Garcia-Zapirain, arXiv:2603.16321 (2026) | Causal-mediation attribution of "genuine quantum contribution"; 13.1:1 direct-vs-quantum-mediated. **Fetch in full.** | **[A]** |
| E5 | **Larocca et al., "Barren plateaus in variational quantum computing", *Nature Reviews Physics* (2025)** | The canonical current review — what does and does not avoid BPs. **Cite instead of McClean 2018 alone.** | **[S]** |
| E6 | **Cerezo et al., "Does provable absence of barren plateaus imply classical simulability?", *Nat Commun* (2025); arXiv:2312.09121** | BP-free landscapes are broadly classically simulable — "a soft form of dequantization". **Live against our exact-statevector and exact-MPS arms; still unaddressed by us.** Carried from S15 V7 and now published. | **[S]** |
| E7 | CV-DV / bosonic variational line: CV-ADAPT-VQE (arXiv:2606.05297); CV quantum annealing (arXiv:2608.04601); variational quantum continuous optimization (arXiv:2210.03136); high-precision qubit-efficient continuous optimisation via amplitude estimation (arXiv:2606.08690) | Bosonic modes carry continuous variables **directly**, "avoiding the qubit overhead associated with discretization and binary encoding". **This is the principled contrast to QTF's phase encoding and to our binary register**, and it is unapplied to PSP. | **[S]** |
| E8 | "De Novo Design of Protein-Binding Peptides by Quantum Computing", arXiv:2503.05458 | Peptide *design*, not folding. Different task; note only. | **[S]** |
| E9 | Zhang et al., *Advanced Science* (2026), doi 10.1002/advs.202513641 — binding-site structure prediction on utility-scale QPUs | The published form of S15's Q11/Q12 line. **Still the only quantum peptide benchmark near our scale (101 peptides, 5–18 res); watch it.** | **[S]** |

### E.2 Structure prediction, refinement, error estimation

| # | source | what it adds | tier |
|---|---|---|---|
| E10 | **Wu, Guo, Cheng, ATOMRefine, *Bioinformatics* 39(5):btad298 (2023)** | **Takes ingredient A.** Predicts coordinate shifts to native; GDT-HA 69.84→70.04 (193 AlphaFoldDB) and 63.91→64.42 (69 CASP14). **The calibration for our flagship's expected effect size.** | **[A]** |
| E11 | **Hiranuma et al., DeepAccNet, *Nat Commun* 12:1340 (2021)** | **Takes ingredient A in the pair-distance basis.** Signed Cβ–Cβ error histograms ("estograms") → Rosetta restraints. Does **not** analyse direction vs magnitude, and does **not** check realizability — the two gaps we occupy. | **[A]** |
| E12 | AlphaFold3-like methods for protein–peptide complexes (bioRxiv 2025.03.09.642277): AF3, Protenix, Chai-1, Boltz-1 | Success rates ~89–91% by Fnat — but on **complexes**, scored by native-contact fraction, not on monomeric 9–16-mer Cα-RMSD. **Not a comparator; useful only to show the field's peptide metric is not ours.** | **[S]** |
| E13 | Gulsevin & Meiler, *Structure* 31(1) (2023) | Unchanged from S15. **The journal-version discrepancy (588 peptides / 10–40 aa / length-normalised) remains UNVERIFIED — ScienceDirect still 403.** | [S15] |
| E14 | 3D-Jury (Ginalski et al.) and the consensus-QA lineage | The pairwise-similarity-predicts-quality principle. Relevant to §B.5 and to our consensus-medoid result. | **[S]** |
| E15 | Ensemble-average pairwise RMSD ↔ RMSF relation (*Biophys J*, 2010, from B-factors) | The structural-biology twin of the ambiguity decomposition. §B.5. | **[S]** |

### E.3 Geometry, kinematics, inverse problems

| # | source | what it adds | tier |
|---|---|---|---|
| E16 | Manipulability ellipsoid (Yoshikawa 1985) and standard robot kinematics | Jacobian SVD, low-singular-value directions, inverse condition number. **The textbook frame for AL1.1.** | **[S]** |
| E17 | López-Blanco & Chacón, iMod (*Bioinformatics* 27:2843, 2011); iMODS (*NAR* 42:W271, 2014); Gō/Noguti torsional NMA; "Optimized torsion-angle normal modes reproduce conformational changes more accurately than Cartesian modes" (*Biophys J* 2011) | Torsion-space NMA builds the same Jacobian; ICS modes beat Cartesian modes with fewer modes. | **[S]** |
| E18 | Coutsias & Seok (protein loop closure IK); Liu & Latombe (structure of the protein-backbone IK map); *R Soc Open Sci* 11:240873 (2024) on self-motion topology | Null spaces of the torsion→endpoint map, dimension and topology. **The literature Sprint 15 flagged as an unrun search pass for the locality theorem — it is real and it is where a competing statement of our null-direction count would live.** | **[S]** |
| E19 | Fonseca, Budday, van den Bedem et al., *PLoS Comput Biol* 11:e1004361 (2015); Kinematic Flexibility Analysis (arXiv:1802.08683); KGSrna | **Nullspace sampling** of a constraint Jacobian: steering in the quiet subspace, named and implemented. | **[A]**/[S] |
| E20 | "Proteins Wriggle" (cond-mat/0108218); biased Gaussian steps in torsional space (cond-mat/0103580); concerted rotation (Dodd/Theodorou; Ulmschneider & Jorgensen) | Move sets designed **to minimise distant-atom displacement per torsion change.** The loud/quiet distinction as MC design principle. | **[S]** |
| E21 | Hansen, discrete Picard condition (*BIT* 30:658, 1990); TSVD as regularisation (*BIT* 27:534, 1987) | The mechanism behind AL1.2, and the reason AL1.3's isotropic-Tikhonov result is the expected one. | **[S]** |
| E22 | Krogh & Vedelsby (NIPS 1995); Ueda & Nakano (1996); Brown et al. *JMLR* 6:1621 (2005) | **The fusion law.** §B.5. | **[S]** + **[D]** |
| E23 | Non-Eckart body-frame torsional NMA (*J Math Chem*, 2012) | The frame/superposition subtlety in torsion-space mode analysis. Relevant to our superposed Jacobian. | **[S]** |

### E.4 Restraints, NMR, chemical shifts

| # | source | what it adds | tier |
|---|---|---|---|
| E24 | NMR restraint-error literature (error-distribution-derived NOE restraints, *J Biomol NMR* 2006; restraint-bound distortion effects) | **Shape-of-error-distribution beats width** — the qualitative twin of K12. §B.7. **Needs a primary fetch.** | **[S]** |
| E25 | Chemical-shift → torsion, 2025–26 sweep | **Re-confirms Sprint 15's absence result: no transformer-based successor to TALOS-N for shift→torsion exists.** Activity in 2025–26 is on the inverse map (structure→shift: NMRNet/SE(3), UCBShift 2.0, sequence-only shift predictors, LEGOLAS) and on small-molecule elucidation. **The channel stays closed on coverage arithmetic; nothing new reopens it.** | **[S]** |
| E26 | Distance geometry underdetermination for short extended chains | Unchanged from S15 C23; supports "the distance prior is the ceiling". | [S15] |

---

## PART F — VERIFICATION STATUS

### Read in full by me this sprint **[F]**
* **arXiv:2609.02113 (QTF).** 9.18 MB PDF downloaded to the scratchpad, `pdftotext -layout`,
  1,801 lines, read plus exhaustive greps for: RMSD protocol terms, `Kabsch`, `superpos`, `align`,
  `terminal residues`, `structural core`, `effective RMSD`, `end-to-end`, `baseline`,
  `classical control`, `random sampling`, `ablation`, `p-value`, `confidence interval`,
  `barren plateau`, `gradient variance`, `degrees of freedom`, `qubits`. Every negative result in
  Part A is a grep over the complete extracted text. Sprint 14's dissection is confirmed wherever
  I re-tested it.
  **Not obtained:** the supplementary material (so the custom Hamiltonian's term-by-term algebra
  remains UNVERIFIED, as in S14) and the Zenodo dataset.

### Verified by my own computation **[D]**
* The ambiguity-decomposition identity `d_avg² = (r₁²+r₂²)/2 − (s/2)²`, five random 40-point trials,
  exact to machine precision, no coplanarity or equidistance needed; and the size of the error
  introduced by substituting the arithmetic for the quadratic mean of `r`. §B.5.
* `P = 2n(ceil(N/n)+3) ≈ 2N + 6n`, hence 132 parameters for 44 torsions (chignolin) and 224 for 88
  (Trp-cage), from QTF's own printed formulas.

### Source code fetched **[C]**
`cumbof/qtf`: `README.md`, `qtf/scoring.py`, `qtf/metrics.py`, `qtf/recipes.py`,
`qtf/engines/qtf.py`, `qtf/assets/recipes/default.yaml`, plus directory listings.
`BlankenbergLab/pheat`: repository description and `src/pheat/metrics.py`.
**Caveat: all of this is WebFetch's summariser rendering the file, not the raw bytes.** The
`end_to_end_target = 5.5` line was returned identically on two independent fetches with different
prompts, which is the strongest corroboration available through this tool. **A follow-up should
clone the repo and read the raw file before this appears in any write-up.**

### Abstract / article-page level **[A]**
arXiv:2609.02113 abstract; ATOMRefine (PMC10191610); DeepAccNet (PMC7910447); arXiv:2308.01807;
arXiv:2308.00044; arXiv:2603.16321.

### Snippet only — **DO NOT QUOTE AS FACT** **[S]**
**Shao & Zhu, *PCCP* 20:7206 (2018)** — the force-field list including **ff14SB and FF14SBonlysc
with GB-Neck2** and the "not identified" conclusion come from a search return; RSC returned 403 and
PubMed a cookie wall. **This is the single most consequential [S] row in the document (§C.3) and it
must be fetched before any AMBER novelty claim is drafted.**
Also snippet-only: Krogh & Vedelsby and the ensemble-diversity line; Nau's pathology papers;
Gao/Schulman/Hilton; Hansen's Picard papers; Yoshikawa/manipulability; iMod/iMODS; Coutsias & Seok;
Liu & Latombe; KGS; "Proteins Wriggle"; the NMR restraint-error line; 3D-Jury; the RMSF/pairwise
identity; Larocca et al.; Cerezo et al. 2025; the CV-DV line; all E-section rows marked [S].

### My own inferences, explicitly not any source's claim **[I]**
* That the QTF hardware CDF decoder's monotonicity biases it toward near-constant torsion profiles
  and therefore toward regular secondary structure (Sprint 14's argument, re-derived).
* That QTF's 5.5 Å hard-coded end-to-end default *may* be chignolin's experimental value; **two
  readings are open and I assert neither.**
* That the native-free quiet-subspace estimate is the Gauss–Newton assumption.
* That the combination of error-direction and quiet-subspace steering is unclaimed — a negative
  search across four literatures.
* That the realizability result (K9) is unclaimed — a negative search.
* Every "SURVIVES"/"OPEN" verdict resting on absence. **Absence searches are weaker than presence
  searches, and I say so at each one.**

### What I did NOT establish
* The *PCCP* 2018 full text (§C.3) — **the highest-priority remaining fetch in the programme.**
* Kang, arXiv:2605.01319 — still unfetched; still needed to finalise the spectrum-to-variance claim
  (S15 claim 2). I did not re-attempt it this sprint.
* Kesiku & Garcia-Zapirain, arXiv:2603.16321, full text — needed before drafting the causal-control
  claim (D.1, N-4).
* The published *Structure* version of Gulsevin & Meiler (ScienceDirect 403, unchanged).
* Krogh & Vedelsby's original NIPS text — the identity is verified by my own computation **[D]**,
  but the citation itself is **[S]** and should be checked against the primary before printing.
* The QTF supplementary material and Zenodo dataset.
* A raw (non-summarised) read of the `qtf` and `pheat` source files.
* The robotics/kinematics search pass for the **exact locality theorem** (S15 claim 1) — I searched
  the adjacent null-space and manipulability literature and found the *quiet-subspace* prior art
  reported in §B.2, but **I did not run a dedicated pass for the `j−i−1` support result.** It
  remains unrun, as it was after Sprint 15.
* I ran no repository code and touched no benchmark artefact.

---

## STATUS: COMPLETE.
