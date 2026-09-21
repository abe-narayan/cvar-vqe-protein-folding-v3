<!--
PROVENANCE: the S31 charter, saved verbatim at sprint open, 2026-09-20.
This comment is the only addition; everything below is the user's message as sent.
S30-L0 claimed to have done this and had not (77 minutes, five prior instances of the same
failure). `s31/s31_verify.py` asserts this path exists.
-->

You are the autonomous lead research scientist and engineering lead for the next sprint of the Protein Folding with CVaR-VQE project.

Repository:
https://github.com/abe-narayan/cvar-vqe-protein-folding-v3

Local repository:
C:\Users\abena\Protein-Folding-Algorithm

You are operating inside the repo with Claude Code. This is not a small coding task. Treat this as a serious computational physics and quantum algorithms research sprint.

Your objective is to discover and, if scientifically justified, deploy the next genuine accuracy breakthrough while preserving the project's scientific integrity.

============================================================

1. CORE OBJECTIVE
   ============================================================

The central scientific objective is:

Improve the canonical 126-target development endpoint from the current approximately 3.21 Å mean built-chain Cα RMSD toward:

PRIMARY TARGET: < 3.00 Å
AMBITION: < 2.50 Å

The main intellectual focus MUST be:

CVaR-VQE.

CVaR-VQE is not a decorative component, not an optional add-on, and not something to replace simply because classical approaches are easier.

Treat the quantum variational optimization problem itself as a central research object:

* objective design
* Hamiltonian design
* CVaR behavior
* non-diagonal structure
* candidate-state representation
* quantum readout
* ansatz geometry
* gradients
* optimization landscape
* expressivity
* information allocation across qubits
* subset selection
* structural observables
* physically meaningful quantum objectives
* connections between the classical candidate distribution and the VQE state
* whether the VQE is solving the wrong optimization problem
* whether the Hamiltonian encodes the wrong notion of usefulness
* whether the readout is extracting the wrong information from a successful optimization
* whether a new quantum formulation can convert already-existing candidate quality into endpoint accuracy

Broad architecture redesign is explicitly allowed.

Do NOT preserve the current architecture merely because it exists.

Do NOT preserve a component merely because it is called "CVaR-VQE."

Preserve the scientific identity of the project, but redesign anything necessary around the CVaR-VQE spine.

The goal is genuine predictive improvement, not cosmetic quantumization.

============================================================
2. CURRENT S30 SCIENTIFIC STATE
===============================

You MUST read the complete S30 report and relevant prior sprint reports/artifacts before doing substantive work.

S30 established:

Production:

* 126 dev targets
* 9–16 amino acids
* canonical production built-chain mean = 3.2105 Å
* median = 2.9661 Å
* current primary target < 3.00 Å has not been reached
* ambitious < 2.50 Å has not been reached

Important distinction:

* 3.0483 Å is the canonical CA point-cloud number
* approximately 3.21 Å is the built-chain endpoint
* future endpoint claims must use the correct built-chain basis

The tail dominates the mean:

* worst 18 mean = 6.2758 Å
* other 108 mean = 2.6997 Å
* capping worst 10 at 3 Å would imply 2.9074 Å
* capping worst 18 at 3 Å would imply 2.7426 Å
* capping worst 30 at 3 Å would imply 2.5778 Å

Therefore tail recovery is disproportionately important.

The current production candidate pool already contains useful structures.

On the worst 18:

* oracle best pool-member mean = 2.5298 Å
* 11/18 contain at least one member under 3 Å

Therefore the problem is not simply that the pool lacks good structures.

The tail is substantially selection-limited.

The current scoring objective fails to reliably order the good candidates inside the pool.

============================================================
3. THE FIVE-BIT / INFORMATION BOTTLENECK
========================================

S30 found a particularly important result:

Five ORACLE signs on selected long-range pairs would move the built-chain endpoint from approximately 3.2126 Å to 2.8867 Å.

This is an ORACLE result.

It is NOT deployable.

Treat it as evidence that the target is reachable if the correct target-specific information can be obtained.

The useful information is concentrated in a very small number of decisions.

The report estimates that approximately five useful bits per target are enough to cross the primary endpoint.

However:

The source enumeration appears to collapse to:

1. sequence
2. candidate library / pool

Physics acts as an operator on those inputs rather than being a third independent target-specific source.

The predictable part of the prior error was found to be strongly associated with pool common-mode structure.

The pool common mode is non-identifiable from pool observations at arbitrary K under the established model.

This means:

Do NOT casually claim "there is no information."

Do NOT claim an information-theoretic impossibility theorem unless you actually prove one.

The correct scientific statement is narrower:

The extractable component identified so far is coherent with the pool common mode and is therefore the wrong component for improving the endpoint.

What remains open is whether an observable exists whose error is incoherent with that common mode and which therefore supplies new target-specific information.

This is the major scientific question of the sprint.

============================================================
4. WHAT S30 CLOSED
==================

Do NOT casually spend the sprint redoing variants of the following unless you develop a materially different theoretical formulation that escapes the exact reason each was closed.

Closed directions:

A. Combining the 21 existing native-free displacement fields

* ORACLE global rho approximately 0.169
* leave-fold-out rho approximately 0.012
* worse than the best individual field
* Gram stable rank approximately 2.06
* the fields effectively span only about two directions

B. Sparse weighted readout in its tested form

* plain argmin dominated every measured bit budget
* 2-of-75 free-weight readout approximately 2.1683 Å at 11.4+ bits
* argmin over 128 approximately 2.1435 Å at 7 bits

BUT:
S30 explicitly leaves open a more powerful sparse-readout formulation if it can get native-free support.

C. Second-moment / quadric escape

* measured effects +0.2059 Å and +0.086 Å
* two lanes
* two samplers
* same harmful direction

D. Subset optimization through the current averaging readout
Theorem T1:
the CVaR tail is always a prefix of the order induced by the local objective gradient at the optimum.

Therefore simply asking current CVaR to magically become a free non-prefix subset selector does not solve the underlying problem.

E. Generative structural-space escape coupled to current averaging
Theorem:
set_mean² approximately B² + S²
where B is the endpoint and spread is free.

F. Large configuration / torsion encodings tested only by brute bit pricing
48 bits versus 7 deployed was judged too expensive in the tested implementation.

BUT:
Do not confuse "expensive encoding" with "scientifically impossible encoding."

A new quantum representation may be worth revisiting if the information is shown to be structurally valuable and the qubit economy is materially improved.

G. Common-mode correction from pool data
Non-identifiable under the established model.

H. Recognizing nativeness from one structure using the audited achiral/local channel family
Ordering survived on only 2 of 43 channels and preference failed on all 43.

I. G1 achiral theorem
Any rotation/translation-invariant reflection-invariant single-structure observable is a function of the distance map.

Therefore do not waste time independently renaming distance-map re-readings as fundamentally new information:

* contact topology
* Rg
* inertia
* burial/SASA
* related achiral geometry summaries

These require a materially different mathematical object to escape the theorem.

J. Chiral single-structure channels as tested
The easy explanation was falsified.
Chiral-axis occupancy was not useful enough.
Best measured contrast was approximately +0.0405 versus null approximately +0.0408.

K. Native-free prior correction as currently formulated
OOF R² approximately 0.83%, below pre-registered 1.96% bar.
Application slightly worsened point-cloud endpoint.

L. Simply making a Hamiltonian off-diagonal
The first graph-Hamiltonian experiment produced a nearly rank-one hopping term with extremely small gradient contribution.
Spread-spectrum graph construction reduced decay but did not produce endpoint gain.

Therefore "make H non-diagonal" is not itself a research result.

You need a physically and mathematically meaningful non-diagonal Hamiltonian if you use one.

============================================================
5. WHAT S30 ESTABLISHED POSITIVELY
==================================

A. The tail is selection-limited, not pool-limited.

B. The native error subspace can be approximately located without the native.
The native-free basis captured approximately 82.73% of ORACLE error in 14 dimensions versus approximately 44.44% for a matched-dimension random frame.

The sign is the bottleneck.

C. The pool is a codebook rather than a new information source.
Seven candidate-index bits can move the ORACLE ladder dramatically.

The index identifies a stored structure.

Therefore candidate indexing should be treated as a serious information-allocation problem.

D. The statistical-potential displacement library wastes many of its principal directions on scale/radial structure.
A radial direction carried approximately 58% of Gram trace.
Dominant PC cosine with radial direction approximately 0.947.
Native direction was approximately orthogonal.

This suggests that useful structural displacement information is present, but the existing field basis spends rank inefficiently.

E. The current quantum circuit is not obviously the expressivity bottleneck.

A 9-qubit depth-3 circuit was able to reach:

* approximately 0.2516 Å built chain
  with an ORACLE objective.

The same broad circuit family with the deployed native-free objective gave:

* approximately 3.4330 Å

Thus:
THE CIRCUIT CAN EXPRESS GOOD SOLUTIONS.
THE CURRENT OBJECTIVE DOES NOT DIRECT THE OPTIMIZER TOWARD THEM.

This is one of the most important conclusions in the entire project.

============================================================
6. CVaR-VQE MUST BE THE MAIN RESEARCH SPINE
===========================================

Treat this sprint primarily as a CVaR-VQE research sprint.

The central question should be:

"What objective and Hamiltonian should a CVaR-VQE optimize so that quantum optimization actually spends its limited information capacity on endpoint-relevant candidate selection?"

You should investigate at least:

1. CVaR objective mathematics
2. Hamiltonian construction
3. objective geometry
4. gradient alignment with endpoint-improving directions
5. state-space geometry
6. candidate indexing
7. readout
8. subset selection
9. structural observables
10. quantum/classical information allocation
11. qubit efficiency
12. whether the quantum state should represent candidates, structural modes, pair decisions, latent coordinates, or combinations thereof
13. whether CVaR should operate over candidate states, structural perturbations, pairwise decisions, or hierarchical choices
14. whether the tail should represent low-energy states, low-cost states, high-information states, or another physically meaningful subset
15. whether the Hamiltonian should encode energies, pair compatibility, long-range geometry, torsion physics, free-energy-like quantities, or a learned-but-scientifically-controlled physical observable
16. whether a hybrid classical-quantum architecture can preserve genuine CVaR-VQE optimization while moving the useful information bottleneck into a better basis

Do not assume that the present 7-qubit candidate-index formulation is optimal.

Do not assume that more qubits automatically help.

Do not assume that deeper circuits automatically help.

Do not assume that CVaR alpha values automatically improve selection.

Do not assume that SPSA/Adam is the problem.

Prove which layer is limiting the result.

============================================================
7. THREE EXPLICIT S30 RECOMMENDATIONS
=====================================

You MUST investigate all three unless a rigorous scientific result proves that one is redundant under a stronger formulation.

---

## 7A. FREE-ENERGY STAGE

Investigate a genuinely meaningful free-energy or thermodynamic stage.

This should not be a cosmetic scalar added to the existing score.

Ask:

Can a physically motivated free-energy-like observable supply target-specific information not reducible to sequence plus the static candidate library?

Possible directions include, but are not limited to:

* conformational free-energy differences
* relative free energies between candidate basins
* entropy/enthalpy decomposition
* local configurational entropy
* effective conformational entropy
* basin populations
* physically grounded ensemble observables
* restrained versus unrestrained free-energy-like quantities
* approximate thermodynamic integration where computationally sensible
* variational free-energy objectives
* minimum-free-energy rather than minimum-static-score formulations
* temperature-dependent candidate ranking
* candidate ensemble reweighting
* state occupancy
* free-energy differences induced by structural interactions

But:

Do not invent physical justification after seeing results.

Pre-register why the quantity should contain information not already contained in the pool common mode.

Determine whether it truly escapes the previous identifiability bottleneck.

Connect it directly to CVaR-VQE wherever possible.

A strong formulation would turn free-energy estimates into a physically meaningful Hamiltonian whose low-energy/CVaR sector corresponds to structurally useful states.

Test cheap proxies before expensive calculations.

Use hierarchy:
cheap theory -> analytical pre-check -> small-scale exact experiment -> fold-aware test -> full 126-target test only if justified.

---

## 7B. BACKBONE-TORSION CHANNEL

Investigate the backbone-torsion channel.

This remains open because it survived the earlier closure logic after controlling for simple size/compactness effects.

This is important because torsional structure may encode information not captured by the audited distance-map-only family.

Test:

* phi/psi distributions
* torsional pair structure
* local torsional strain
* torsional coupling
* Ramachandran basin occupancy
* sequence-conditioned torsion structure
* long-range torsional correlations
* chirality-sensitive torsional quantities
* torsion-space distances
* torsion-space geodesics
* conditional torsion entropy
* torsion transition structure
* backbone-state compatibility
* candidate-to-candidate torsional topology

But do not merely produce another correlated scalar.

The channel must:

1. survive the appropriate band/partialling tests,
2. show information not already explained by the current cost,
3. demonstrate some relationship to the in-pool selection problem,
4. ideally produce a direction that is useful for CVaR-VQE optimization.

Find a terminal operator that can actually spend this information.

Possible quantum formulations include:

* torsion-state qubits
* pairwise torsion compatibility Hamiltonians
* torsion-basin occupation variables
* local/global torsion constraints
* torsion-energy terms
* QUBO/Ising encodings
* non-diagonal transitions between torsion states
* CVaR over torsion-compatible candidate subsets
* hierarchical candidate -> torsion basin -> structural state optimization

Do not force the representation into the old candidate-index architecture if a better quantum formulation emerges.

---

## 7C. SPARSE WEIGHTED READOUT

Revisit sparse weighted readout, but solve the actual problem identified by S30.

Do NOT simply repeat the previous 2-of-75 experiment.

The question is:

Can a native-free mechanism identify a small number of candidate members and weights whose weighted structural combination produces a substantially better structure?

The ORACLE result shows the potential is large.

Investigate:

* sparse support selection
* learned support under strict native-free training
* group sparsity
* entropy-regularized support
* differentiable top-k
* Gumbel-like relaxation if scientifically justified
* sparse quantum measurement
* amplitude-based support
* candidate mixture models
* nonuniform CVaR readout
* posterior-weighted candidate averaging
* candidate clustering followed by quantum selection
* adaptive candidate cardinality
* support/weight co-optimization
* quantum-selected subset followed by deterministic geometric refinement

The key requirement is:
native-free support.

No deployable support selection may use target-native RMSD.

Treat oracle support results strictly as upper-bound diagnostics.

============================================================
8. GO BEYOND THE THREE RECOMMENDATIONS
======================================

You have broad freedom.

Investigate anything that can plausibly create a real CVaR-VQE breakthrough, including but not limited to:

* alternative Hamiltonian families
* physically motivated Ising/QUBO formulations
* sparse non-diagonal Hamiltonians
* pairwise compatibility Hamiltonians
* graph Hamiltonians with meaningful spectra
* hopping between structurally compatible candidates
* structural-mode Hamiltonians
* torsional Hamiltonians
* free-energy Hamiltonians
* hierarchical VQE
* multi-stage CVaR
* nested CVaR
* adaptive alpha schedules
* candidate-dependent alpha
* temperature-dependent CVaR
* CVaR with physically meaningful regularization
* robust variational objectives
* risk-sensitive variational optimization
* distributional objectives
* quantile matching
* coherent candidate superpositions
* quantum kernel ideas where genuinely useful
* low-rank Hamiltonian factorization
* operator compression
* basis redesign
* qubit-efficient encodings
* symmetry reductions
* adaptive operator pools
* ADAPT-VQE
* qubit-ADAPT-VQE
* problem-inspired ansätze
* data-reuploading only if it has a scientific purpose
* candidate graph spectral bases
* latent structural mode bases
* quantum PCA-like approaches where meaningful
* information-theoretic qubit allocation
* mutual-information-guided qubit assignment
* bit-value optimization
* readout optimization
* hybrid quantum/classical selection
* quantum optimization over candidate indices
* quantum optimization over pair decisions
* quantum optimization over structural modes
* quantum optimization over torsional states

You may invent new formulations.

But every new idea must answer:

"What target-specific information does this contain that the current pipeline does not?"

If you cannot answer that, deprioritize it.

============================================================
9. REQUIRED VQE OBJECTIVE DIAGNOSTICS
=====================================

For every serious VQE formulation, calculate more than final RMSD.

At minimum measure:

A. Objective versus RMSD correlation.

B. Local gradient alignment:
cosine between objective gradient/displacement and native-improving direction.

C. Objective-step experiment:
small steps in objective descent direction and multiple step sizes.

D. Gradient variance versus system size.

E. Optimization stability.

F. CVaR alpha dependence.

G. Hamiltonian spectrum.

H. spectral gap.

I. degeneracy.

J. effective rank.

K. participation ratio / inverse participation ratio.

L. entropy.

M. probability mass in useful candidate states.

N. candidate support concentration.

O. rank of the selected native candidate.

P. bit efficiency.

Q. information per qubit.

R. readout sensitivity.

S. trainability.

T. parameter-shift versus finite difference where appropriate.

U. statevector versus MPS agreement where applicable.

V. optimizer dependence.

W. random-objective controls.

X. shuffled-label controls.

Y. matched-information nulls.

Z. fold-aware uncertainty.

The key test is:

Does changing the VQE objective change the direction of optimization toward endpoint-relevant structures?

Do not settle for "the optimizer converges."

============================================================
10. REQUIRED ATTENTION TO THE OBJECTIVE BOTTLENECK
==================================================

The current evidence says:

Circuit expressivity is not the primary barrier.

Therefore do NOT spend the majority of the sprint merely:

* increasing circuit depth
* increasing optimizer iterations
* changing SPSA seeds
* tuning Adam learning rate
* adding arbitrary layers
* increasing restarts

Those are secondary until an objective with endpoint-aligned geometry exists.

For every architecture, explicitly answer:

1. What is the Hamiltonian?
2. What physical or structural quantity does it represent?
3. Why is it not just a reparameterized version of the old classical rank ladder?
4. Why should minimizing it correlate with endpoint quality?
5. What information does each qubit encode?
6. What does CVaR mean physically in this formulation?
7. What state does the optimized wavefunction represent?
8. What does measurement/readout mean?
9. How does the readout turn quantum information into a 3D structure?
10. What mechanism moves the answer toward the native-free useful subspace?
11. Why should the gradient point in a useful direction?
12. What is the falsifier?

============================================================
11. IMPORTANT THEOREM: DO NOT FOOL YOURSELF WITH CLASSICAL CVaR
===============================================================

The prior formulation had a diagonal Hamiltonian based on rank-standardized candidate scores.

For a diagonal H, CVaR over the measured distribution can reduce to classical tail selection.

S28/S29/S30 established the prefix behavior and related theorems.

Therefore:

Do not claim quantum advantage merely because a quantum circuit is optimizing a classical diagonal ordering.

For any new formulation, explicitly determine whether:

* CVaR is doing something genuinely quantum/nonclassical,
* or merely implementing a classical prefix/rank operation through a quantum state.

If it remains classically reducible, that can still be useful computationally, but do not describe it as a fundamentally new quantum optimization mechanism.

The burden is to find a quantum formulation where the state, Hamiltonian, interference, spectrum, or readout creates useful structure that the current classical procedure cannot cheaply reproduce.

============================================================
12. CANDIDATE-INDEX INFORMATION ALLOCATION
==========================================

Treat the 7-bit candidate index as a serious bottleneck.

The project has evidence that candidate-index bits are potentially much more valuable than bits spent on other choices.

Do not assume the binary encoding is optimal.

Investigate:

* Gray-code candidate encoding
* learned binary encodings
* locality-preserving index encodings
* hierarchical indexing
* cluster-first encoding
* coarse-to-fine indexing
* structure-aware indexing
* candidate graph embeddings
* spectral indexing
* bit-balanced candidate partitions
* adaptive register allocation
* mixed-radix encodings
* compressed candidate labels
* quantum-address-like formulations
* variable candidate cardinality

The goal is not just to represent 128 candidates.

The goal is to make the qubits correspond to scientifically meaningful candidate distinctions.

The 128->512 register widening has an existing undischarged gate.

Treat it as potentially useful, but do not run it blindly.

First resolve the circularity and preregister the test.

============================================================
13. FAIL18 MUST RECEIVE SPECIAL ATTENTION
=========================================

The current endpoint is heavily tail-dominated.

The FAIL18 targets are:

1ID6
1JBF
1LB7
2BFI
2BP4
2JN5
2MQ2
2N5C
2NB7
2NDM
3BTB
3SGO
5W52
7JS6
7LCW
8T63
9KAR
9L1M

Use FAIL18 as a diagnostic subset.

But do not tune a deployable method directly to these targets.

Use them to understand mechanism.

Ask:

* What distinguishes FAIL18?
* Is candidate selection failing differently?
* Does the quantum objective fail more strongly?
* Are useful candidates present but badly ranked?
* Is the native-improving subspace different?
* Does torsional structure differ?
* Does long-range structure differ?
* Does free-energy-like structure differ?
* Do candidate clusters behave differently?
* Does candidate index information have different value?
* Are there specific graph/spectral structures?

Do not make claims from outcome-defined subsets without appropriate controls.

============================================================
14. OPEN SCIENTIFIC QUESTION: NEW PHYSICAL OBSERVABLE
=====================================================

The most important open question from S30 is:

Can you find an observable of the molecule itself that is:

* not reducible to sequence,
* not merely a re-reading of the deposited pool,
* not simply a distance-map scalar,
* not just another common-mode predictor,
* and whose error is incoherent with the pool common mode?

This is the key requirement for a genuinely new channel.

Investigate physically meaningful observables.

Potential directions include:

* torsional thermodynamics
* conformational entropy
* free-energy differences
* dynamical/transition observables
* hydrogen-bond rearrangement patterns
* electrostatic interaction structure
* solvent-mediated effects
* local energetic frustration
* pairwise coupling observables
* correlated torsion transitions
* nonlocal structural compatibility
* state occupancy
* energy-landscape curvature
* basin connectivity
* graph spectral properties
* kinetic surrogates
* elastic/network observables
* many-body structural quantities

But:

You must establish whether the observable genuinely escapes existing equivalence classes and whether it can be computed without native information at inference.

============================================================
15. INFORMATION THEORY
======================

Use serious information-theoretic analysis.

Investigate:

* mutual information
* conditional mutual information
* information gain per qubit
* candidate-rank information
* entropy reduction
* effective candidate count
* posterior concentration
* common-mode decomposition
* independent-member effective sample size
* Fisher information
* rate-distortion style analysis
* channel capacity where appropriate
* information bottleneck formulations
* bit-value curves

Do not turn a fitted correlation into a grand information-theoretic claim.

A theorem is a theorem.

An empirical estimate is an empirical estimate.

An ORACLE quantity is an ORACLE quantity.

Keep those categories separate.

============================================================
16. PHYSICS AND MATHEMATICS STANDARD
====================================

Use maximum mathematical and physical depth.

You are encouraged to derive new results analytically before running experiments.

Use:

* variational calculus
* perturbation theory
* spectral analysis
* matrix inequalities
* convexity/concavity
* information geometry
* Fubini-Study geometry
* Riemannian geometry
* Lie algebra structure
* statistical mechanics
* thermodynamics
* graphical models
* Bayesian decision theory
* quantum measurement theory
* risk-sensitive optimization
* random matrix intuition
* compressed sensing ideas
* manifold geometry
* low-rank structure
* subspace estimation
* perturbation bounds
* optimization landscape analysis

Where useful, derive theorems before touching the full 126-target benchmark.

The strongest result is not "we tried 400 configurations."

The strongest result is:

"Here is the mathematical reason this architecture should or should not contain the information we need, and here is the smallest experiment that falsifies it."

============================================================
17. RESEARCH STYLE
==================

Operate like a highly skeptical research group.

For each idea:

HYPOTHESIS
-> MECHANISM
-> CHEAP PRE-CHECK
-> FALSIFIER
-> SMALL TEST
-> FOLD-AWARE TEST
-> 126-TARGET TEST
-> DEPLOYMENT ONLY IF JUSTIFIED

Do not spend large compute on an idea that can be cheaply falsified.

But do not prematurely kill an idea based on weak/noisy evidence.

Respect the project's MDE standard.

MDE = 2.8016 × SE.

Rules:

* below 0.7× MDE = not a result
* 0.7–1.0× MDE = NOT MEASURED
* > 1.0× MDE requires appropriate uncertainty interpretation
* use paired target-level comparisons
* use fold-clustered confidence intervals
* use concentration-aware nulls where appropriate
* report Type-M where useful
* distinguish exploratory measurements from registered tests

Never create a W/L scorecard that hides uncertainty.

============================================================
18. ORACLE DISCIPLINE
=====================

Oracle experiments are allowed and encouraged for diagnosis.

They MUST ALWAYS be explicitly labeled:

ORACLE
NOT DEPLOYABLE

Native information may be used to establish:

* expressivity ceilings
* theoretical readout ceilings
* attainable subspaces
* ideal candidate selections
* best possible signs
* upper bounds
* failure decomposition

Native information MUST NOT be used to tune any deployable parameter.

Do not quietly turn an oracle construction into a proposed method.

Do not report oracle results as achievements of the model.

============================================================
19. REPRODUCIBILITY RULES
=========================

The canonical benchmark is sacred.

DO NOT:

* reopen benchmark60
* regenerate its folds
* regenerate cluster assignments
* regenerate sealed benchmark splits
* change target identities
* retune production on native RMSD

Preserve all existing sprint artifacts.

Do not erase failed experiments.

Do not overwrite prior reports.

Do not silently recompute existing numbers under new definitions.

Every meaningful number needs:

* artifact path
* experiment ID
* configuration
* seed if applicable
* endpoint basis
* whether cloud or built-chain
* whether native-free or ORACLE
* comparison baseline

Pin the cloud-to-built-chain projection seed.

Resolve the current reproducibility discrepancy:

s29 canonical production:
3.2105 built-chain
3.0483 cloud

Other reprojections differed by as much as approximately 0.0107 Å due to an unpinned multi-start projection seed.

Any future sub-0.01 Å claim is invalid unless this is fixed.

============================================================
20. REQUIRED OPEN BUG FIXES
===========================

Before or during the research, inspect and fix the known integrity defects.

A. Withdrawn positive in:
core/pipeline.py:821

It still asserts the old CVaR-tail +0.113 Å measured role.

S25 withdrew that result.

Current measured role was approximately -0.1405 Å at 0.68× MDE.

Remove/correct the stale claim.

B. Projection seed
Pin the cloud->chain projection seed and make endpoint reproduction deterministic.

C. Governor / launcher inconsistency
Inspect:
s26/jobrun.py
s26/governor.py

Current behavior reportedly has:

* jobrun refusing launches above 85% CPU
* governor tolerating/suspending around 101%

This conflicts with the desired sustained approximately 94–95% resource utilization regime.

Resolve this cleanly and safely.

Do not make the machine unstable just to satisfy utilization numbers.

D. 128->512 widening
Resolve the existing circularity concern and either:

* formally register a valid experiment,
* run it,
* or document a rigorous reason to close it.

============================================================
21. COMPUTE RESOURCE POLICY
===========================

Use aggressive parallel research.

Maintain approximately:

CPU utilization: 94–95%
RAM utilization: 94–95%

Treat this as the operating target.

BUT:
do not exceed safe machine capacity,
do not intentionally trigger swapping,
do not freeze the OS,
do not corrupt processes,
do not create uncontrolled process explosions.

The goal is approximately 94–95% sustained utilization, not reckless oversubscription.

Always keep at least enough headroom for:

* OS stability
* Claude Code
* file operations
* logging
* process cleanup
* emergency termination

If the current governor conflicts with the target utilization, fix the governor/launcher architecture rather than bypassing safety blindly.

Use existing resource-control infrastructure where sensible.

Prefer one AMBER/OpenMM process at a time unless you have strong evidence the machine can safely support more.

============================================================
22. MULTI-AGENT POLICY
======================

Maintain at least 4 active research agents/threads at all times during substantive exploration.

Never use more than 8 simultaneously.

You decide whether the active count should be:
4
5
6
7
or 8

based on workload.

Default structure:

Agent 1:
Quantum/CVaR theory lead

* CVaR mathematics
* Hamiltonians
* spectra
* objective geometry
* VQE formulations
* quantum information

Agent 2:
Physical-model lead

* free energy
* torsions
* structural physics
* physical observables
* candidate compatibility
* thermodynamics

Agent 3:
Encoding/readout lead

* candidate indexing
* sparse readout
* qubit allocation
* subset formulations
* information theory
* quantum measurement

Agent 4:
Experimental verification lead

* preregistration
* controls
* artifacts
* endpoint calculations
* statistical validity
* reproducibility

Spawn additional agents when useful, up to 8 total.

Claude decides:

* whether additional agents are useful
* how tasks should be split
* which hypotheses deserve parallel attacks
* when agents should be reassigned
* which lanes should be merged

Do not create agents merely to create activity.

Every agent must have a scientific purpose.

============================================================
23. RESEARCH PRIORITIZATION
===========================

Use this order conceptually, but adapt dynamically if evidence strongly suggests another path:

Priority 1:
Find a genuinely better CVaR-VQE objective/Hamiltonian whose local geometry points toward useful structures.

Priority 2:
Free-energy stage.

Priority 3:
Backbone-torsion channel.

Priority 4:
Native-free sparse readout.

Priority 5:
Candidate-index redesign / bit allocation.

Priority 6:
Physically meaningful non-diagonal Hamiltonian.

Priority 7:
Hybrid quantum/classical formulations preserving genuine CVaR-VQE optimization.

Priority 8:
Only then, secondary ansatz/optimizer tuning.

The ordering is not absolute.

A surprising scientific result can reorder the sprint.

============================================================
24. REQUIRED EXPERIMENTAL LADDER
================================

Every candidate architecture should move through:

LEVEL 0:
mathematical derivation

LEVEL 1:
synthetic toy model

LEVEL 2:
small real subset

LEVEL 3:
fold-held-out subset

LEVEL 4:
126-target development instrument

LEVEL 5:
deployment candidate

Do not jump directly from a toy demonstration to the final endpoint.

============================================================
25. REQUIRED CONTROLS
=====================

For any promising new method, compare against:

* current production
* random candidate ranking
* plain argmin
* current distogram ranking
* relevant classical structural score
* matched-information random controls
* shuffled-label control
* native-free null
* ORACLE upper bound where relevant

Do not use "best found" as a baseline unless it is a properly defined baseline.

============================================================
26. IMPORTANT QUESTION ABOUT CVaR
=================================

Investigate deeply whether the current CVaR use is conceptually correct.

Questions:

Should CVaR operate on:

* candidate energies?
* candidate posterior losses?
* structural compatibility?
* pairwise decisions?
* free energies?
* torsional states?
* candidate subsets?
* latent structural modes?
* an ensemble of physically realizable structures?

Should the tail represent:

* low energy?
* high posterior probability?
* robustly good states?
* uncertainty-sensitive states?
* rare but physically meaningful conformations?

Could a two-level CVaR be useful?

For example:

CVaR over candidate families
followed by
CVaR within family

or:

CVaR over structural modes
followed by
candidate selection

or:

classical posterior
-> quantum state
-> CVaR over quantum structural choices
-> geometric reconstruction

Explore such formulations.

============================================================
27. DO NOT ASSUME AVERAGING IS THE RIGHT READOUT
================================================

The existing uniform structural average is not sacred.

Investigate whether VQE should instead produce:

* one candidate
* sparse mixture
* weighted mixture
* coherent structural mode combination
* mode coefficients
* pair-selection decisions
* torsion choices
* local structural corrections
* displacement coefficients
* candidate cluster selection
* coarse-to-fine reconstruction

But every readout must preserve physical validity.

Do not produce impossible coordinates just to reduce a cloud RMSD.

Built-chain endpoint is the governing metric.

============================================================
28. PHYSICAL VALIDITY
=====================

Preserve or strengthen:

* peptide bond geometry
* chirality
* torsional validity
* structural continuity
* steric sanity
* physically meaningful conformations

AMBER/GBn2 remains valid as an all-atom validity/refinement tool.

Do not automatically assume it is an accuracy ranker.

S26/S27 evidence showed that physics-only rankers can be harmful for the current endpoint.

Use AMBER primarily where it has a justified scientific role.

============================================================
29. REPORTING REQUIREMENTS
==========================

At the end, produce a complete sprint package.

Required:

1. REPORT_S31.md

2. LEDGER.md

3. machine-readable experiment records

4. verification script

5. before/after endpoint summary

6. complete hypothesis table

7. closed directions

8. surviving directions

9. open questions

10. reproducibility information

11. exact git commit(s)

12. list of modified files

13. resource-utilization summary

14. agent/task summary

15. mathematical derivations used

16. statistical methods

17. explicit native-free/deployable versus ORACLE separation

18. exact endpoint basis

19. exact experimental seeds

20. all negative results

Do not hide failed experiments.

Do not delete old artifacts.

============================================================
30. REPORT STRUCTURE
====================

Use a structure approximately like:

§0 Answer up front

§1 Endpoint and reproducibility

§2 Scientific hypothesis

§3 Mathematical mechanism

§4 CVaR-VQE architecture

§5 Hamiltonian

§6 Quantum encoding

§7 Objective geometry

§8 Free-energy work

§9 Torsion work

§10 Sparse/readout work

§11 Candidate-index information

§12 FAIL18 analysis

§13 Controls and nulls

§14 126-target endpoint

§15 Statistical analysis

§16 What was closed

§17 What improved

§18 What did not improve

§19 ORACLE ceilings

§20 Remaining information bottleneck

§21 Next-sprint recommendation

============================================================
31. SUCCESS CRITERIA
====================

A genuine success is NOT:

* lower training loss
* lower Hamiltonian energy
* better optimizer convergence
* prettier plots
* deeper circuits
* more qubits
* higher entropy
* lower cloud RMSD alone
* an ORACLE result
* a one-off target improvement
* a sub-MDE trend

A genuine success means a statistically and scientifically credible improvement on the canonical 126-target built-chain endpoint.

Primary:
mean < 3.00 Å

Ambitious:
mean < 2.50 Å

A result in the 2.99 Å range is meaningful only if it survives the project's statistical standards and is reproducible.

Do not engineer a fragile 2.98 result just to cross a threshold.

============================================================
32. IF THE PRIMARY TARGET IS NOT REACHED
========================================

A rigorous ceiling result is acceptable.

But only if it materially advances the science.

If you fail to beat production, determine:

* exactly where information is lost
* exactly where the VQE objective fails
* exactly what information the quantum representation can/cannot encode
* whether a theorem can explain the bottleneck
* whether the remaining barrier is identifiable
* whether the remaining barrier is information-theoretic, computational, representational, or merely architectural
* what experiment would actually distinguish those possibilities

A rigorous negative result is valuable.

A pile of weak negative experiments is not.

============================================================
33. ABSOLUTE SCIENTIFIC INTEGRITY RULES
=======================================

NEVER:

* fabricate a result
* alter a result to fit a hypothesis
* tune deployable parameters with native RMSD
* reopen sealed benchmark60
* regenerate frozen splits
* hide failed runs
* silently change endpoint definitions
* confuse cloud and built-chain metrics
* call ORACLE results deployable
* call correlation causation
* call 0.5× MDE a null result with high confidence
* claim "information-theoretically impossible" without an actual theorem
* call a classical diagonal rank operation a quantum advantage
* assume a complex Hamiltonian is useful merely because it is non-diagonal
* optimize harder on a scientifically invalid objective

NEVER trust prose over artifacts.

If a report says a job finished:
check the artifact.

If a file is being written:
do not interpret its partial contents as a completed result.

If a number disagrees with an artifact:
the artifact wins, and the discrepancy must be recorded.

============================================================
34. MAXIMUM CREATIVE FREEDOM
============================

You have broad freedom to:

* redesign data representations
* redesign the Hamiltonian
* redesign the quantum state
* redesign the readout
* redesign the candidate-selection problem
* redesign the VQE objective
* introduce new mathematically justified observables
* introduce new classical preprocessing
* introduce physically meaningful intermediate stages
* combine classical and quantum components differently
* invent new experiments
* derive new theorems
* create new diagnostics
* refactor code substantially
* replace obsolete components
* simplify where appropriate

But every architectural decision must be justified by information flow and physics.

Prefer elegant mechanisms over piles of heuristics.

============================================================
35. RESEARCH DEPTH
==================

Use maximum:

scientific reasoning
mathematical reasoning
quantum physics
statistical mechanics
optimization theory
information theory
protein structural physics
computational experimentation

Research the relevant literature when it can materially improve the design.

Use current authoritative sources for modern methods/software details when necessary.

Do not blindly copy papers.

Translate literature ideas into this exact problem.

============================================================
36. FINAL DECISION POLICY
=========================

Claude, you are responsible for deciding:

* which experiments deserve compute
* which hypotheses survive
* whether to use 4, 5, 6, 7, or 8 active agents
* which agents should work together
* whether to introduce new mathematical machinery
* whether to redesign the VQE entirely
* whether to deploy a result
* whether a result is scientifically strong enough
* whether a negative result should close a direction
* where the next compute budget should go

Do not ask me for permission for ordinary scientific decisions.

Make the best research decisions autonomously.

============================================================
37. STARTUP PROCEDURE
=====================

Before coding:

1. inspect git state
2. inspect current branch
3. inspect current untracked/modified files
4. read S30 report
5. read relevant S29/S28 reports
6. inspect current production pipeline
7. inspect s26 governor/jobrun infrastructure
8. inspect current VQE code
9. inspect Hamiltonian code
10. inspect candidate generation/indexing code
11. inspect diagnostics
12. inspect existing experiment ledger
13. verify that benchmark60 remains untouched and sealed

Then produce an internal research map.

Do not stop and wait for me.

Immediately launch the multi-agent research program.

============================================================
38. FIRST RESEARCH MOVES
========================

Start with parallel work on:

Lane A:
derive and test a genuinely non-classical CVaR-VQE formulation whose objective is structurally aligned with candidate selection

Lane B:
free-energy observables and thermodynamic formulations

Lane C:
backbone torsion information and a quantum terminal operator

Lane D:
native-free sparse readout plus candidate-index/qubit allocation

Then allow additional agents to attack whichever directions appear strongest.

In parallel:

* fix the stale CVaR prose
* pin projection reproducibility
* resolve governor/launcher inconsistency
* investigate 128->512 gate
* preserve all current artifacts

Do not let engineering cleanup consume the entire sprint.

============================================================
39. FINAL MANDATE
=================

Treat this as a serious attempt to solve the actual bottleneck.

The current evidence strongly suggests:

The candidate pool contains useful structures.
The existing scoring/selection objective cannot reliably select them.
The current CVaR-VQE circuit can express excellent states when given the right objective.
Therefore the central scientific task is to construct a physically meaningful, information-efficient, native-free CVaR-VQE objective/Hamiltonian/readout whose optimization geometry points toward the useful part of the candidate distribution.

Find that mechanism.

Do the mathematics.

Do the physics.

Do the experiments.

Use the compute aggressively.

Use 4 active agents at all times, never more than 8, with Claude deciding whether additional agents are useful.

Target approximately 94–95% CPU and RAM utilization while preserving machine stability.

Be maximally creative.

Be maximally skeptical.

Be maximally rigorous.

Preserve every useful result and every negative result.

Do not optimize for a pretty sprint report.

Optimize for a real scientific breakthrough.

BEGIN.

Lots of freedom, creativity, ambition. Main goal is RSMD. Report only after everything is done. keep cpu/ram usage ~95%. I believe in you! You got this!
