# Sprint 13 — LITERATURE AND NOVELTY FINDINGS

Agent: `lit`. Date: 2026-09-05. Method: external literature only (WebSearch / WebFetch).
**No compute run.** No repository code executed, no files outside `s13/` touched.

Structured data: `s13/results/lit_methods.json` (50 rows, four strands).

Every claim below is tiered per BRIEF §1.4. Literature claims are **LITERATURE-SUPPORTED**
where the primary source was fetched and **LITERATURE-REPORTED (UNVERIFIED)** where only a
search snippet was seen — the JSON's `verified` field marks each row `primary` or `secondary`.
Inferences I drew from those sources are marked **HYPOTHESIS** or **PROPOSED**.

---

## 0. THINGS THE PROJECT CURRENTLY ASSUMES THAT THE LITERATURE CONTRADICTS

These lead because they are the most valuable output. Details and sources follow in §1–§4.

**C1. The torsion encoding this sprint is built on was published in October 2025.**
`arXiv:2510.06413` uses per-residue φ/ψ bins, binary index encoding, qubit count linear in
residue count, VQE + CVaR on a 127-qubit IBM processor, on **75 fragments of 10–14 residues**,
fused with a NetSurfP dihedral/secondary-structure prior. That is `PerResidueTorsion` with
`n_bits = n_res·log₂k`, at this length, with a learned torsion prior, one year ago.
QFold (Casares/Campos/Martín-Delgado) encoded φ/ψ in 2–5 bits per angle with a distilled-AlphaFold
torsion prior in 2021. Torsion-space quantum optimisation *with a real molecular force field*
dates to Marchand et al. 2018 (D-Wave, "flexible with respect to the choice of molecular force
field") and Mato et al. 2022 (Quantum Molecular Unfolding, HUBO over rotatable bonds).
**"Torsion-constrained VQE for peptides" is not new and must be removed from the novelty claim.**
LITERATURE-SUPPORTED (primary).

**C2. CVaR is not a barren-plateau mitigation, and the sprint should stop treating it as one.**
Qiu, Lumbreras, Li & Rebentrost (`arXiv:2605.02850`, 2026) prove the tilted-loss family that
contains CVaR "does not remove the barren plateau problem by itself": fixed tilt retains
exponential-in-n gradient-variance decay identical to the standard objective. What tilting does
is move the bottleneck "from trainability to **estimability**: from a flat landscape to a larger
sampling overhead" — as α falls, the effective sample size shrinks proportionally and the total
measurement budget must strictly increase. LITERATURE-SUPPORTED (primary).

**C3. The strongest published result in this exact problem area is negative, and the sprint has
no control for it.** Boulebnane, Lucas, Meyder, Adaszewski & Montanaro (npj QI 2023, 20 qubits,
noiseless): *"the performance of QAOA can be matched by random sampling up to a small overhead …
these results cast serious doubt on the ability of QAOA to address the protein folding problem in
the near term, even in an extremely simplified setting."* A random-sampling arm at matched
objective-evaluation budget is therefore **mandatory**, per energy model. Without it, every
positive VQE result in this sprint is presumed refuted by default. LITERATURE-SUPPORTED (primary).

**C4. TALOS-N's 90% / 12° is real, but both of its gating mechanisms fail on this target class,
and the gaps are clustered.** The 12° figure is an RMS against **crystallographic** (φ,ψ) on
**folded globular proteins**, computed over the residues that survived two filters:
(i) predictability requires ≥3 assigned shifts on ≥2 of {i−1, i, i+1}; (ii) the matching unit is a
**heptapeptide** (±3 residues) and confidence is graded by whether the top 25 / top 10 matches
agree; (iii) residues with **RCI-S2 ≤ 0.6 are labelled `Dyn`** and carry no usable restraint.
On a 13-mer, 6 of 13 residues sit in a truncated ±3 window, and an isolated peptide in water is by
construction low-S2. Both failure modes are **spatially clustered** (termini; flexible segments) —
the worst case for a sequentially built chain, and exactly the case Sprint 12's `fill_pool`
experiment showed is worse than random gaps. TALOS+'s reference database excludes chains ≤ 20
residues, so nothing in the TALOS family's corpus resembles this project's targets.
LITERATURE-SUPPORTED (primary), transfer argument HYPOTHESIS.

**C5. SPARTA+ is not a torsion predictor — it is the inverse map.** The brief lists it among
"chemical-shift-to-torsion methods". SPARTA+ (Shen & Bax 2010) takes **coordinates** and returns
**chemical shifts**, RMS 2.45 / 1.09 / 0.94 / 1.14 / 0.25 / 0.49 ppm for ¹⁵N / ¹³C′ / ¹³Cα /
¹³Cβ / ¹Hα / ¹HN, trained on 580 high-resolution X-ray structures. It cannot supply a φ/ψ prior.
It *can* supply something the sprint needs more: a native-coordinate-free objective-validity
instrument (§1.4). LITERATURE-SUPPORTED (primary).

**C6. Sequence-only torsion prediction cannot reach 2.0 Å, and this is now measured rather than
guessed.** SPOT-1D-LM (Sci Rep 2022), the language-model state of the art: **φ MAE 15.99°,
ψ MAE 23.74°** on TEST2018; **20.67° / 36.57°** on TEST2020; single-sequence SPOT-1D-Single
**22.16° / 40.58°** and **22.92° / 44.25°**. These are on full crystallised domains
(resolution < 2.5 Å, R-free < 0.25, ≤ 25 % identity); no test set contains 9–16-mers, and a 13-mer
has an MSA of depth ≈ 1, i.e. the Neff1 regime. Sprint 12's own coverage table gives
**2.408 Å at σ = 20° with full coverage**. A deployable sequence-only torsion prior therefore
lands around **2.4–2.9 Å** — better than the 3.213 Å incumbent, short of the target.
LITERATURE-SUPPORTED (primary) + arithmetic against s12's measured curve.

**C7. A lattice CVaR-VQE already reports 1.22–3.11 Å on peptides of comparable length.**
Kannan et al., `arXiv:2510.15316` (Oct 2025): FCC turn encoding, 5(N−1) qubits, MJ potential,
CVaR-VQE on IBM Heron, 13 peptides of 6–10 residues plus chains to 20, RMSD **1.22–3.11 Å**
against experimental references (their classical MD control: 0.79–5.32 Å). The sprint's < 2.0 Å
target will be read against this number. LITERATURE-SUPPORTED (primary).

**C8. The most likely fatal error in the geometry strand — the Hamiltonian cannot change the
metric.** For |ψ(θ)⟩ = U(θ)|ψ₀⟩ the Fubini–Study metric
g_ij(θ) = Re[⟨∂_iψ|∂_jψ⟩ − ⟨∂_iψ|ψ⟩⟨ψ|∂_jψ⟩] is a function of the ansatz, the input state and θ
**only**. Legacy and AMBER cannot produce different metrics at the same θ. They produce different
gradients, hence different trajectories θ(t), hence different metric values *visited*.
The one genuine exception is a **Hamiltonian-generated ansatz** (HVA / QAOA), where H defines the
generators — with a problem-independent hardware-efficient or MPS ansatz that exception does not
apply. Note also that "Hamiltonian-aware QNG" in the literature (`arXiv:2511.14511`) means
Hamiltonian-aware *metric estimation* — choosing which blocks are worth measuring — not a
Hamiltonian-dependent metric. LITERATURE-SUPPORTED (definitional + primary).

---

## 1. STRAND 1 — TORSION INFORMATION FROM REALISTIC CHANNELS

### 1.1 The TALOS lineage — verified claim and its evaluation set

| method | year | coverage | accuracy | error rate | evaluation |
|---|---|---|---|---|---|
| TALOS (orig.) | 1999 | ~65 % | — | ~3 % | 20-protein database |
| TALOS+ | 2009 | **88.5 % "Good"** | **±13°** (excl. ~2.5 % gross errors) | 2.46 % | folded proteins; DB = X-ray < 2.0 Å, chains **> 20 residues** |
| **TALOS-N** | 2013 | **≥ 90 %** | **~12° RMS** vs crystallographic (φ,ψ) | ~3.5 % | "an independent set of proteins"; matching DB of 9,523 proteins |

Verbatim (Shen & Bax 2013, *J Biomol NMR* 56:227): *"Validation on an independent set of proteins
indicates that backbone torsion angles can be predicted for a larger, ≥90 % fraction of the
residues, with an error rate smaller than ca 3.5 %, and a root mean square difference between
predicted and crystallographically observed (φ, ψ) torsion angles of ca 12°."*

**The figure is verified.** Three properties of it matter more than its value:

1. **It is post-filter.** The 12° is computed over predicted residues, after the ≥10 % declined
   and the ~3.5 % erroneous are removed. It is not a 12° error over all residues.
2. **The reference is a rigid crystal structure.** It is a rigid-body accuracy, not an accuracy
   against a solution ensemble. For a peptide whose NMR observable *is* an ensemble average,
   the object being predicted is not the same object.
3. **The corpus contains nothing of this length.** TALOS+'s database explicitly excludes chains
   ≤ 20 residues.

### 1.2 What decides coverage, and why the gaps cluster

Three gates, all from the primary sources (Shen & Bax 2013; protocol chapter *Methods Mol Biol*
1260:17, PMC4319698):

- **Shift-completeness gate.** *"If at least two of the three residues have at least three chemical
  shifts, the center residue is considered to be predictable."* Missing assignments propagate to
  neighbours, so an unassigned residue costs up to three predictions.
- **Heptapeptide window.** Confidence is graded by whether *"all 25 best matching heptapeptides
  locate in a consistent φ/ψ region"* (Strong) or only the top 10 (Generous). The window is ±3.
  On a 13-mer, residues 1–3 and 11–13 — **6 of 13, 46 %** — have a truncated window.
- **Dynamics gate.** *"Residues below the threshold RCI-S2 ≤ 0.6 are assigned as dynamic"* (`Dyn`)
  and carry no usable torsion restraint. Standard practice in restraint generation is stricter
  still: exclude S2 ≤ 0.7.

The protocol chapter's worked example (DinI, 81 residues, near-complete assignment, well folded):
**71 unambiguous, 2 ambiguous, 6 dynamic, M1 unpredicted** — 87.7 % on the easy case, and the
residual is termini plus dynamic stretches, i.e. **contiguous**.

**HYPOTHESIS (high confidence).** On isolated 9–16-mers the realised coverage will be materially
below 90 % and the gaps will be clustered at the termini and in the flexible core. Sprint 12
measured the consequence: the channel is worth **+0.078 Å [−0.147, +0.298]** — nothing, and
possibly negative — at 50 % coverage, and only reaches 2.0 Å at ≥ 90 % coverage and ≤ 12°.

**This is measurable now, from shifts alone.** RCI-S2 needs only chemical shifts and sequence.
The 79 BMRB-linked targets in `s12/results/lit_bmrb_coverage.json` are enough to compute, for each:
the fraction of residues with ≥3 assigned shifts, the fraction whose ±1 neighbourhood passes the
predictability gate, the fraction with RCI-S2 > 0.6, and — critically — **the run-length
distribution of the gaps**. No structure pipeline, no TALOS-N install, no heavy compute.
**PROPOSED, and it is the cheapest kill-or-continue decision in the sprint.**

### 1.3 Modern sequence-only φ/ψ prediction — the deployable ceiling

| method | φ MAE | ψ MAE | test set |
|---|---|---|---|
| SPOT-1D-LM (ProtTrans + ESM-1b) | 15.99° | 23.74° | TEST2018 (250 proteins) |
| SPOT-1D (MSA profile) | 16.88° | 24.87° | TEST2018 |
| SPOT-1D-Single (single sequence) | 22.16° | 40.58° | TEST2018 |
| SPOT-1D-LM | 20.67° | 36.57° | TEST2020 (671 proteins) |
| SPOT-1D-Single | 22.92° | 44.25° | TEST2020 |

Test sets: TEST2018 = 250 proteins released Jan–Jun 2018, resolution < 2.5 Å, R-free < 0.25,
< 25 % identity to pre-2018 structures. TEST2020 = 671 proteins, HMM-based homolog removal.
Neff1-2020 = 46 proteins with no homologs — the closest analogue to a peptide with no MSA.
No test set contains short peptides. 2024–2026 successors (DCMA, OPUS-TASS2) were checked and
none shows a regime change; treat ~16° φ / ~24° ψ as the state of the art.
LITERATURE-SUPPORTED (primary for SPOT-1D-LM; secondary for successors).

Cross-referenced against Sprint 12's own coverage curve, this closes the sequence-only route to
< 2.0 Å (C6 above). It is also consistent with the sprint's own leave-fold-out 4-state torsion-bin
predictor (0.690 overall, 0.517 on FAIL18) — that result is not an anomaly, it is what the field
delivers on out-of-distribution short sequences.

### 1.4 The one thing SPARTA+ is actually good for here

SPARTA+ / SHIFTX2 are **forward** maps: coordinates → shifts. That makes them a
**native-coordinate-free objective-validity instrument**, which BRIEF §5 demands and the project
does not currently have. For each candidate structure, predict shifts and score agreement against
the deposited BMRB shifts; the resulting score uses no native coordinates and no native torsions.
It directly answers §5's four questions (native percentile, Spearman ρ vs CA-RMSD, argmin quality,
scatter) with an experimental observable rather than an oracle.

Caveat: SHIFTX2's homology-transfer module looks up structurally homologous proteins with assigned
shifts, which is a memorisation path — use the ML module only, and label the arm. Leakage caveat
inherited from Sprint 12 travels with everything in this section: the deposited coordinates were
determined *using* these shifts plus NOEs, so any shift-driven arm is **NMR-restrained structure
determination reported in its own column**, never an improvement to the sequence-only system.
PROPOSED.

---

## 2. STRAND 2 — BARREN PLATEAUS, WITH REGIMES OF VALIDITY

The taxonomy below follows the Nature Reviews Physics review (Larocca et al. 2025,
`arXiv:2405.00781`). Its framing sentence is the one to hold onto: *"all the moving pieces of an
algorithm — choices of ansatz, initial state, observable, loss function and hardware noise — can
lead to BPs when ill-suited."* The **observable** is on that list, which is precisely why this
sprint's experiment is well posed.

| mechanism | claim | **regime — outside this it does not apply** | applies here? |
|---|---|---|---|
| **Expressibility / 2-design** (McClean 2018) | E[∂C] = 0, Var[∂C] exponentially small in n | random init **and** the circuit ensemble approximating a unitary 2-design (i.e. deep/expressive enough) | **Only if the ansatz 2-designs.** An MPS ansatz with bond dimension capped at 2^layers, or a shallow chain, is *not* in this regime. Do not cite McClean for flatness seen there. |
| **Cost globality** (Cerezo 2021) | global observable → exp. vanishing gradients **even at shallow depth**; local observable → at worst polynomial decay **for depth O(log n)** | *"assuming V(θ) is an alternating layered ansatz composed of blocks forming local 2-designs"*; local-cost guarantee only up to O(log n) depth | **YES — the central theorem for this sprint.** Both energies are diagonal in the torsion bitstring basis, hence weighted sums of Pauli-Z strings; "globality" is their Pauli-weight distribution. |
| **Lie-algebraic unification** (Ragone 2024) | Var = P_M(ρ)·P_M(O)/dim(M) | circuit 2-designs over its dynamical Lie algebra | **YES conceptually.** With ρ, U(θ) and dim(M) fixed across arms, any variance difference is attributable to **P_M(O)** — the observable factor. This is the language to state the hypothesis in. |
| **Entanglement-induced** (Ortiz Marrero 2021) | volume-law entanglement → deterministic exponentially flat landscape | volume law between measured and traced-out subsystems, or O(1/2ⁿ) module projection of ρ | **Probably not** — the bond-dimension cap is an explicit entanglement bound. Report entanglement entropy across the chain cut to close this alternative explanation. |
| **Noise-induced** (Wang 2021) | gradient vanishes exponentially in n if depth grows linearly with n, under local Pauli noise; unital noise concentrates the loss deterministically toward Tr[O]/2ⁿ | requires a **noisy device** and depth growing with n; deterministic, so init tricks and error mitigation do not fix it | **NO** for this sprint's numbers — statevector / MPS / lightning are noiseless. Say so explicitly. |
| **Traps, not plateaus** (Anschuetz & Kiani 2022) | shallow models with **no BP** still have only a superpolynomially small fraction of local minima within any constant energy of the global minimum | random init; the specific model class analysed | **YES.** Absence of a plateau is *not* evidence of trainability. Solution quality must be a separate reported axis. |
| **CVaR / tilted loss** (Qiu 2026) | tilting *"does not remove the barren plateau problem by itself"*; fixed tilt → exponential-in-n decay; bottleneck moves *"from trainability to estimability"* | fixed tilt is general; the polynomial lower bound needs structured settings and a Θ(−n) schedule | **YES.** See C2. |
| **Simulability critique** (Cerezo et al., Nat Commun 2025) | many models with provably BP-free landscapes also admit classical simulation given an initial classical-data-acquisition phase — *"a soft form of dequantization"*; mechanism is that BP-avoidance confines the evolved observable to a poly-size subspace | case-by-case evidence, not a theorem; attacks the *information-processing* advantage | **YES, and unaddressed.** The repo runs exact statevector and exact MPS of its own circuits — it is a live instance of the critique. |

Two further review points that are routinely misstated and should not be:

- **"Log-depth + local cost avoids barren plateaus" is false as stated.** The review is explicit:
  *"a logarithmic-depth hardware efficient ansatz can still exhibit a BP for a non-local
  measurement."* The observable must lie in a polynomially large module.
- **Narrow gorges.** Probabilistic BPs come with minima *"always exponentially narrow (meaning
  their relative volume in parameter space is exponentially small)"*. Finding a good point by luck
  does not contradict a plateau, and a plateau does not mean the landscape is featureless — it
  means the features are unfindable from random initialisation.

### 2.1 The measurement this literature implies

**PROPOSED — highest-value experiment in the sprint.** Both Legacy and AMBER are diagonal
functions of the torsion bitstring, so each is *exactly* a weighted sum of Pauli-Z strings, and
the weight distribution is computable: exactly by Walsh–Hadamard transform for n_bits ≲ 20
(k = 4 with ≤ 10 residues; k = 2 with ≤ 20), by sampling above that. Cerezo 2021 then converts
that spectrum into a prediction about gradient-variance scaling.

This is cheap, classical, needs no VQE run, and states BRIEF §4's intuition ("a distance between
residues i and j is |i−j|-local, not 2-local") as a number the audience's own theorems are written
in. It also gives the sprint a mechanism rather than a correlation: *AMBER is harder to train
because its spectral mass sits at higher Pauli weight* is a falsifiable claim; *AMBER is harder to
train* is not.

Independent support for the same physical point comes from the annealing literature: Quantum
Molecular Unfolding formulates torsion-space conformational optimisation as a **HUBO** — high
order, not quadratic — and quadratises to QUBO with ancillas; Marchand et al. restrict to
torsion *neighbourhoods* because only a subset can be jointly optimised. Two independent groups
hit the same |i−j|-locality wall.

---

## 3. STRAND 3 — QUANTUM NATURAL GRADIENT AND INFORMATION GEOMETRY

**The definitional fact the geometry agent must not get wrong** (C8): the Fubini–Study metric /
quantum geometric tensor of |ψ(θ)⟩ = U(θ)|ψ₀⟩ depends on the ansatz, the input state and θ, and
on nothing else. Stokes et al. (*Quantum* 4:269, 2020) define QNG as steepest descent with respect
to *"the real part of the Quantum Geometric Tensor (QGT), also known as the Fubini-Study metric
tensor"* — an object of the state family. The observable enters the *gradient*, never the metric.

Consequences for this sprint:

- Legacy and AMBER **must** give identical metrics at matched θ. Use that as a **unit test** on the
  geometry code: sample random θ, compute g under both arms, assert equality to numerical
  precision. A discrepancy means either a bug or an unintended Hamiltonian-dependent ansatz.
- Every geometry sentence should read *"the energy model selects a different region of a fixed
  state manifold"*, never *"the energy model reshapes the geometry"*.
- The exception is real and must be stated: with a **Hamiltonian-variational ansatz** (HVA/QAOA)
  the generators *are* the problem terms, so the manifold itself is energy-model-dependent and the
  metrics legitimately differ. Choose the ansatz deliberately and say which regime you are in.
- Beware the terminology trap: "Hamiltonian-aware quantum natural gradient" (`arXiv:2511.14511`)
  means Hamiltonian-aware *estimation* — selecting which metric blocks are worth measuring — not a
  Hamiltonian-dependent tensor.

**Approximations.** QNG is usually run with a block-diagonal (per-layer) or diagonal
approximation. Diagonal QNG is a per-parameter rescaled gradient — a preconditioned SGD — and must
not be described as natural-gradient descent on the state manifold. Report which one was used.

**Singularity and regularisation.** The FS metric is *generically singular* under redundant
parameterisation, so QNG updates are ill-posed; practice uses the Moore–Penrose pseudo-inverse or
Tikhonov g + λI. Reported curvature can depend on the regulariser — for Tikhonov the value at a
point can depend on how λ → 0 is taken — and trajectories can pass through genuinely singular
points. **Report λ, show stability under a λ sweep, and prefer regularisation-robust summaries
(rank, condition number) over raw eigenvalues.** LITERATURE-REPORTED (secondary).

**The right ansatz-side diagnostic.** Haug, Bharti & Kim (PRX Quantum 2021) define the *effective
quantum dimension* as the **rank of the QFI/QGT** — the number of independent directions the
circuit can actually move the state in; a rank deficit is parameter redundancy. Because it is
observable-independent, it should be reported **once per ansatz, not once per energy model**.
Doing so is the cleanest demonstration that the comparison is controlled, and it pre-empts the
reviewer's first question: did both arms have the same expressive power?

**Which Fisher information?** Abbas et al. (Nat Comput Sci 2021) use the **classical** Fisher
information of the model's output distribution to define an effective dimension. That is a
different object from the QFI of the state, and the two are routinely conflated. For a CVaR
objective over sampled bitstrings the natural object is the classical Fisher information of
p_θ(x) — which depends on the measurement basis but still not on the energy model. Name the
object explicitly; this audience will notice.

**When geometry-aware optimisation helps.** The honest summary of the literature is: QNG helps
when the parameterisation is badly conditioned relative to the state manifold — i.e. when Euclidean
steps in θ produce wildly unequal steps in state space — and it does not help when the metric is
near-isotropic, when the metric estimate is noisy, or when the failure is a trap rather than
curvature. It does **not** cure barren plateaus: an exponentially small gradient divided by an
exponentially small metric is still not measurable at finite shots.

---

## 4. STRAND 4 — NOVELTY AUDIT

### 4.1 What already exists (be specific)

**Quantum protein folding on lattices — a crowded, 14-year-old field.**
Perdomo-Ortiz et al. (Sci Rep 2012, D-Wave, HP lattice, up to 81 qubits) →
Robert, Barkoutsos, Woerner & Tavernelli (npj QI 2021: tetrahedral turn encoding, O(N⁴) scaling,
**10-aa Angiotensin on 22 qubits**, 7-aa neuropeptide on 9 qubits, CVaR) →
Chandarana et al. (dCD-VQE, 9 aa on 17 qubits) →
FCC-lattice encodings (`arXiv:2507.08955`, 2025) →
trapped-ion hardware demonstrations at ~64 qubits (2025–26) →
penalty-free formulations (`arXiv:2606.02104`, 2026).
Resource line from the field's own perspective (Doga et al., JCTC 2024): 118 qubits for 22 aa,
41 aa on 433 qubits, 67 aa on 1121 qubits, ~10,000 qubits for a 141-residue protein.

**CVaR-VQE for peptide folding — taken, repeatedly.**
Barkoutsos et al. (*Quantum* 2020) introduced CVaR aggregation; Robert et al. 2021 applied it to
folding; Uttarkar & Niranjan (QIP 2024 / IJBM 2024) ran CVaR-VQE against MD on 7-mers with
structural-alphabet analysis; Kannan et al. (2025) ran CVaR-VQE on IBM Heron reaching
**1.22–3.11 Å** on 6–20-residue peptides; QuPepFold (PLOS One 2026) packages it
(α = 0.025 noiseless / 0.05 on hardware, 1,224 sequences of 6–10 residues, ~1.76 M conformers).

**Torsion / dihedral quantum encodings — taken, and recently.**
Marchand et al. 2018 (annealer, torsion neighbourhoods, arbitrary force field) →
Adaszewski & Montanaro et al. 2022 (Frontiers: *"a specific encoding of molecular dihedral angles
into registers of qubits and a method for implementing, in quantum superposition, a Markov chain
Monte Carlo update step based on a classical all-atom force field"* — fault-tolerant resource
estimate) → QFold 2021 (φ/ψ in 2–5 bits, Minifold von-Mises torsion prior, Szegedy quantum walk,
di- to tetrapeptides, TTS only, no RMSD) → Mato et al. 2022 (torsion HUBO on D-Wave) →
**`arXiv:2510.06413` (Oct 2025): per-residue φ/ψ bins, binary indices, linear qubit scaling,
VQE + CVaR on 127 qubits, 75 fragments of 10–14 residues, NetSurfP dihedral prior fused into the
objective, mean 4.89 Å.**

**NMR-restrained quantum structure determination — nothing found.** No paper was found combining
chemical-shift-derived torsion restraints with a quantum optimiser. This appears genuinely open —
but it is open partly because it is a small idea (restraints are restraints, whatever optimises
them) and partly because of C4.

**All-atom force fields inside a quantum optimisation — essentially nothing.** Doga et al. (JCTC
2024) survey the field and report that all-atom force fields (FF99SB) appear only in **classical
post-processing refinement**, never inside the quantum computation, and that MJ potentials are
*"the gold standard for determining interaction energies"* in quantum PSP. Adaszewski/Montanaro
proposed all-atom-inside-the-loop but as a fault-tolerant resource estimate, not a demonstration.

### 4.2 Testing the coordinator's hypothesis

> *"the novelty is not 'using torsions' and not 'using VQE for folding', but the controlled
> comparison of two molecular energy models of differing physical realism under matched
> variational conditions, with geometry and trainability as the measured observables."*

**The hypothesis survives, and the first two clauses are confirmed emphatically.** Torsions: taken
(C1). VQE for folding: taken (§4.1). Verdict on the third clause:

- **No study was found that varies the energy model as the independent variable under matched
  variational conditions.** The existing comparison literature varies the *optimiser* with the
  energy model fixed (Uttarkar & Niranjan: CVaR-VQE vs MD; QuPepFold: CVaR-VQE vs expectation-value
  VQE vs MD; Kannan: VQE vs simulated annealing). Every one of them holds MJ fixed.
- **The field's own perspective names the gap.** Doga et al. (JCTC 2024) state that no empirical
  comparison exists showing whether coarse-grained results correlate with all-atom predictions,
  and that a rigorous resource estimate for a given PSP problem was beyond their scope.
- **The nearest methodological neighbour is on the wrong system.** `arXiv:2609.00235` (2026)
  benchmarks eight optimisers for exact-statevector VQE on **frustrated spin models** under
  matched function-evaluation budgets, characterising Hamiltonian–ansatz landscapes via local
  minima, gradients, curvature and reachability, with restart-based basin coverage as a control.
  Matched-budget landscape characterisation is therefore **established practice** — the sprint gets
  no novelty credit for the protocol, but it gets a citable template and a defence of BRIEF rule 9.
  *(Secondary source — open it before citing.)*

**Verdict: the contribution is real but narrower than the sprint brief implies, and it is a
trainability paper, not a folding paper.** The defensible claim is roughly:

> *A controlled study, at matched representation / ansatz / initialisation / CVaR α / optimiser /
> seed / objective-evaluation budget, of how the Pauli-weight structure of a molecular energy model
> governs the gradient-variance scaling and the state-manifold trajectory of a variational
> optimisation — instantiated as a knowledge-based potential versus a genuine all-atom
> ff14SB/GBn2 force field on a torsion-space peptide encoding.*

Three things make that non-trivial rather than an exercise: (i) an all-atom force field has never
been the objective *inside* a variational loop in this domain, and `core.amber.single_point` at
6 ms warm is what makes it possible; (ii) the two energies are natural, physically motivated
objects rather than synthetic Hamiltonians constructed to differ in locality; (iii) the sprint has
a 126-target instrument, which nothing in §4.1 has — every prior-art peptide set is between 13 and
75 hand-assembled targets.

**If it is to be sharpened, sharpen it here.** The single most valuable addition is the
**Walsh/Pauli-weight spectrum** (§2.1). It converts "AMBER is more complex" into a measured
spectral statement, connects directly to Cerezo 2021 and Ragone 2024, and would be, as far as this
survey found, the first time a real molecular force field's Pauli-weight structure has been
measured and tied to VQA trainability. That is a cleaner and more citable contribution than the
folding accuracy.

**What is definitively NOT novel and must not be claimed:** the torsion encoding; CVaR-VQE;
CVaR-VQE for peptides; peptide folding on quantum hardware; discretised torsion libraries; the
matched-evaluation-budget protocol; comparing a quantum optimiser with MD.

---

## 5. ANALYSIS — actionable implications for the sprint design

Each states what should change and why. Ordered by expected value.

**A1. Add a random-sampling control at matched objective-evaluation budget, per energy model, before
any VQE arm is reported.** *Why:* Boulebnane et al. (npj QI 2023) found QAOA on a simplified folding
potential is matched by random sampling up to a small overhead at 20 qubits, noiseless. That is the
strongest published prior in this exact area and it is negative. Sprint 12's own record points the
same way (a certified exact optimum emitted 0.169 Å *worse* than its anchor). Without this control
every positive result is presumed refuted. It is also nearly free — the same `BudgetedEnergyModel`
call count with uniform bitstrings.

**A2. Make the Walsh/Pauli-weight spectrum of E(bitstring) the first experiment, before any
optimisation.** *Why:* Cerezo et al. 2021 proves globality (Pauli weight), not "complexity",
controls gradient scaling; Ragone et al. 2024 localises the observable's contribution to a single
factor P_M(O). Both energies are diagonal, so their spectra are exactly computable at n_bits ≲ 20
and estimable above. This gives the sprint a *mechanism* instead of a correlation, and it is
classical, cheap, and cannot fail to produce a result. It also directly discharges BRIEF §4 — the
|i−j|-locality claim becomes a measured weight distribution.

**A3. Any barren-plateau claim requires a scaling sweep in n, not a variance measurement.**
*Why:* a BP is a statement about Var[∂C] as n → ∞ under random initialisation. A single-size
variance difference between Legacy and AMBER is not evidence of one and will be read as a category
error by this audience. The representation gives two clean knobs: n_bits = n_res·log₂k. Pre-register
a sweep (e.g. n_bits ≈ 8 → 28 via k ∈ {2,4,8} × n_res), ≥ 200 random parameter draws per point,
fit log Var vs n, and report the fitted slope with a CI. Report the ansatz's **effective quantum
dimension** (QFI rank) once per ansatz to show expressibility was held fixed.

**A4. Rewrite every geometry claim as trajectory selection, and add the metric-equality unit test.**
*Why:* the Fubini–Study metric of U(θ)|ψ₀⟩ cannot depend on the observable. Claiming AMBER
"reshapes the metric" is the single error most likely to end the conversation with this professor.
Report g(θ) along both trajectories **and** at matched random θ, where it must be identical — an
assertion in the test suite. State the HVA exception and which side of it the chosen ansatz sits on.
Report the regulariser λ and show stability under a λ sweep; prefer rank and condition number to
raw eigenvalues.

**A5. Stop describing CVaR as a trainability aid, and decouple evaluation-budget parity from shot
parity.** *Why:* Qiu et al. 2026 — tilting "does not remove the barren plateau problem by itself";
fixed tilt keeps exponential-in-n decay; the bottleneck becomes estimability, with effective sample
size shrinking as α falls. Practical consequence for BRIEF rule 9: if α is held equal across arms
this is neutral, but **any α sweep must report shots-to-a-fixed-gradient-SNR alongside evaluation
counts**, or budget parity is silently violated in the direction that flatters small α.

**A6. Measure realised TALOS-N-style coverage on the 79 BMRB-linked targets from shifts alone, and
pre-register a kill threshold, before committing to the NMR channel.** *Why:* Sprint 12 measured
that this system needs ≥ 90 % coverage at ≤ 12° to reach 2.0 Å and that the channel is worse than
useless below 50 %. TALOS-N's three gates (≥3 shifts on ≥2 of {i−1,i,i+1}; ±3 heptapeptide window;
RCI-S2 ≤ 0.6 → `Dyn`) all bite hardest on isolated 9–16-mers, and the resulting gaps are clustered
at termini and flexible segments — the worst case for a sequentially built chain. Report per-target
coverage **and the gap run-length distribution**, not just the mean. Suggested pre-registered kill
threshold: median coverage < 75 % ⇒ close the direction. Cost: shift parsing only, no structure
pipeline, no heavy compute.

**A7. Reclassify sequence-only torsion prediction as a ~2.4 Å route, not a < 2.0 Å route, and stop
funding accuracy work against that target.** *Why:* SPOT-1D-LM delivers φ MAE 16.0°/20.7° and
ψ MAE 23.7°/36.6° on *full crystallised domains*; single-sequence (the 13-mer regime) is 22.2°/40.6°
and 22.9°/44.3°. Sprint 12's own curve gives 2.408 Å at σ = 20° with full coverage. The arithmetic
closes the route. The < 2.0 Å claim requires the NMR channel and must be reported in its own column
as NMR-restrained structure determination, never folded into the sequence-only mean.

**A8. Move SPARTA+ from the torsion-prior list to the objective-validity instrument list.**
*Why:* it is the inverse map (coordinates → shifts, RMS 0.25–2.45 ppm depending on nucleus). As a
forward model it gives BRIEF §5 something it currently lacks: a **native-coordinate-free** score
for ranking candidates against an experimental observable. Use the ML module only if SHIFTX2 is
substituted (its homology-transfer module is a memorisation path).

**A9. Cite `arXiv:2510.06413`, QFold, Marchand 2018 and Mato 2022 explicitly, and remove the
torsion encoding from the novelty claim.** *Why:* the per-residue-bin binary torsion encoding with
linear qubit scaling on 10–14-residue fragments was published in October 2025 with a learned torsion
prior. Being scooped is survivable; being scooped and not citing it is not, and it is exactly the
kind of thing a VQA professor's group will already know. Replace the encoding claim with the
positioning argument the literature actually supports: the field is spending qubits on finer
lattices (tetrahedral → FCC), the torsion representation spends them on real backbone geometry, and
it is **natively penalty-free** — every bitstring is a valid conformation, so there are no
self-avoidance or chirality penalty terms and no constraint-weight hyperparameters, which lattice
encodings cannot say.

**A10. Add one paragraph engaging the classical-simulability critique, and state plainly that no
quantum advantage is claimed.** *Why:* Cerezo et al. (Nat Commun 2025) argue that models provably
free of barren plateaus tend to admit classical simulation — "a soft form of dequantization" — and
this repository is already exact-simulating its own circuits with `StatevectorCircuit` and
`MPSAnsatz`. Get there first: every number in the sprint is classically simulated; the contribution
is a controlled measurement of how problem structure shapes VQA landscapes; here is where a genuine
advantage could live and why it is not demonstrated here. Volunteering this is worth far more than
having it extracted.

**A11. Report solution quality on a separate axis from gradient statistics, and never infer one from
the other.** *Why:* Anschuetz & Kiani (Nat Commun 2022) prove shallow models with **no** barren
plateau can still have only a superpolynomially small fraction of good local minima. "AMBER's
gradient variance is no worse than Legacy's" would therefore not license "AMBER is trainable". Keep
BRIEF §5's four objective-validity measurements (native percentile, Spearman ρ vs CA-RMSD, argmin
vs random, scatter) as a first-class result table alongside the variance scaling.

**A12. Benchmark against Kannan et al.'s numbers explicitly, and lead with the instrument.**
*Why:* a lattice MJ CVaR-VQE reports 1.22–3.11 Å on 6–20-residue peptides. The sprint's < 2.0 Å
target sits inside that band, and the obvious question — why does an off-lattice all-atom system not
beat a lattice statistical potential? — needs an answer prepared. The honest and strong answer is
the instrument: their 13 hand-picked targets versus this project's 126-target set with pinned folds,
FAIL18 stratification, bootstrap CIs and drop-top-k. If bandwidth allows, run one or two of their
peptides through the pipeline for a direct anchor.

---

## 6. WHAT I DID NOT ESTABLISH

- **No realised coverage number for TALOS-N on short peptides.** No published study reports
  TALOS-N coverage on isolated 9–16-mers. C4 is a mechanism-based prediction, tiered HYPOTHESIS.
  A6 is the experiment that would settle it, and it is cheap.
- **I did not run TALOS-N, RCI, SPARTA+ or any predictor.** No compute was used.
- **I did not read four primary sources whose rows are in the JSON**: `arXiv:2606.01611` (CD-QAOA
  peptides), `arXiv:2606.02104` (penalty-free lattice folding), `arXiv:2506.07866` /
  `arXiv:2604.26861` (trapped-ion folding), and `arXiv:2609.00235` (VQE landscape geometry — the
  methodological template). All are flagged `verified: secondary — NOT READ` in the JSON. Open them
  before citing; `2609.00235` in particular is the closest methodological neighbour and could
  narrow the novelty claim further.
- **I did not verify the TALOS-N validation set's composition** (how many proteins, which). The
  paper says "an independent set of proteins"; the PDF would not render in this environment
  (no poppler). The gates and the worked example are verified; the set size is not.
- **I did not do an exhaustive patent or non-English search**, and I did not search conference
  proceedings (QIP, TQC) where an unpublished version of this comparison could exist.
- **I did not price the Walsh-spectrum experiment (A2) in compute.** It is classical and cheap at
  n_bits ≲ 20 but I did not measure it; the compute agent should.
- **I did not establish whether the 2.4 Å estimate in C6/A7 survives** the difference between MAE
  (what SPOT-1D-LM reports) and the σ of the Gaussian corruption Sprint 12 used. Circular error
  distributions are heavy-tailed and MAE understates RMS; the true figure is probably *worse* than
  2.4 Å, not better, but I did not compute the conversion.
