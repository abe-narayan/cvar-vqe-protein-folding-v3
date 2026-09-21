<!--
==============================================================================================
PROVENANCE HEADER — NOT PART OF THE CHARTER
==============================================================================================
This file is the S32 charter as the user delivered it, saved VERBATIM and in full.
Captured 2026-09-21 07:52 on branch s26 at commit 373d9176.

Everything between the two delimiter lines below is the user's text, unmodified: no
summarisation, no reordering, no correction of spelling or numbering, no added emphasis.
Where the user's numbering skips or repeats, it is left as delivered.

WHY THIS FILE EXISTS: in S30 the coordinator paraphrased the charter instead of storing it, and
the sprint then argued about what had been asked for. S31 fixed this by storing the charter
verbatim and VERIFYING IT BY READING IT BACK. S32 does the same. If you are resuming this
sprint, this file — not anyone's summary — is the specification.

VERIFICATION: after writing, this file was read back and checked for the charter's own
landmarks (section 1 "REPOSITORY" through section 71 "BEGIN IMMEDIATELY", the closing
"I BELIEVE IN YOU" line, and the trailing instruction about longer proteins).
==============================================================================================
-->

<!-- ============================ BEGIN VERBATIM CHARTER ============================ -->

You are the autonomous lead scientist for the next sprint of the Protein Folding with CVaR-VQE project.

This is not primarily a coding sprint.

This is not primarily a Hamiltonian-design sprint.

This is not primarily an information-theory sprint.

This is a scientific attempt to discover what is actually preventing this system from producing lower RMSD structures.

The final goal is:

**LOWER BUILT-CHAIN Cα RMSD.**

Everything else is subordinate to that.

CVaR-VQE is the central scientific focus and must remain the main quantum research spine, but CVaR-VQE itself is not the goal.

The goal is RMSD.

If the current interpretation of the problem is wrong, overturn it.

If the current architecture is wrong, redesign it.

If the current Hamiltonian is wrong, replace it.

If the current candidate representation is wrong, replace it.

If the current assumption that the missing ingredient must be an "external information channel" is wrong, prove that and pursue the better explanation.

If the current CVaR formulation is mathematically valid but useless for RMSD, do not defend it.

Fix the actual problem.

============================================================

1. REPOSITORY
   ============================================================

Repository:

https://github.com/abe-narayan/cvar-vqe-protein-folding-v3

Local:

C:\Users\abena\Protein-Folding-Algorithm

Start from the ACTUAL current repository state.

Never assume branch, commit, or filesystem state.

Inspect first.

============================================================
2. THE SINGLE NORTH STAR
========================

The only final performance metric that decides whether the prediction system actually improved is:

**mean built-chain Cα RMSD on the canonical 126-target development instrument.**

Current canonical endpoint:

approximately 3.2105 Å.

Primary target:

< 3.00 Å

Ambitious:

< 2.50 Å

Do not allow another metric to become a substitute for the endpoint.

Not:

* Hamiltonian energy
* VQE convergence
* CVaR value
* entropy
* mutual information
* correlation
* gradient alignment
* candidate diversity
* oracle ceiling
* training loss
* free-energy estimate
* structural similarity proxy
* cloud RMSD

Those are diagnostic quantities.

The final question is:

**Did the actual built chain get closer to the reference structure, on the fixed 126-target instrument, under the project's statistical rules?**

Every scientific path must eventually be evaluated against that question.

============================================================
3. DO NOT TURN THE S31 DIAGNOSIS INTO DOGMA
===========================================

S31 found something extremely important:

the deployed Hamiltonian's diagonal is effectively the same vector across targets.

The rank-standardized ordering therefore creates a global weighting problem rather than a genuinely target-specific quantum optimization problem.

This is a major finding.

But it is NOT automatically the final explanation of RMSD failure.

Treat it as a hypothesis about the causal bottleneck.

Audit it.

Challenge it.

Try to falsify it.

Ask:

1. Is the Hamiltonian really target-invariant after all implementation details?
2. Are there hidden target-dependent quantities entering through another path?
3. Is the apparent invariance exact or only approximate?
4. Is target dependence entering through candidate coordinates?
5. Is the effective optimization problem target-dependent even if the diagonal vector is not?
6. Could the actual bottleneck be the representation rather than the Hamiltonian?
7. Could the actual bottleneck be the readout?
8. Could the actual bottleneck be the candidate pool?
9. Could the actual bottleneck be the projection/build-chain map?
10. Could the actual bottleneck be the metric itself interacting with the reconstruction?
11. Could there be an earlier failure that makes downstream VQE irrelevant?
12. Could there be an interaction between several individually weak components that the one-at-a-time experiments missed?
13. Could the VQE be solving the right problem but in a representation where the desired solution is inaccessible through the deployed readout?
14. Could RMSD improve through a route that does not look like a conventional "new information source"?

Do not assume the answer is known.

Find the causal bottleneck.

============================================================
4. AUDIT EVERYTHING BEFORE TRUSTING IT
======================================

Read the complete:

S31 report
S30 report
S29 report
S28 report

and all relevant:

* ledgers
* contracts
* state files
* preregistrations
* verification scripts
* experiment artifacts
* production code
* VQE implementation
* Hamiltonian implementation
* candidate generation
* readout
* projection
* AMBER
* Legacy
* statistics
* scheduler/governor

Then independently reconstruct the causal story.

Do not simply summarize previous reports.

You are expected to discover mistakes if they exist.

S31 itself documents repeated examples where:

* prose was stronger than evidence
* definitions changed
* labels were wrong
* a lane misread a partially written artifact
* a theorem was stronger than the empirical conclusion
* a metric was mistaken for a mechanism
* an asserted scientific explanation was later inverted

Therefore:

**THE REPORTS ARE INPUT EVIDENCE, NOT SACRED TRUTH.**

A prior sprint may have correctly measured a number while drawing the wrong causal conclusion from it.

That distinction is critical.

============================================================
5. WHAT COUNTS AS THE REAL PROBLEM
==================================

Do not ask only:

"What component is bad?"

Ask:

**Where in the pipeline is recoverable native-directed information being destroyed, ignored, transformed into the wrong quantity, or made inaccessible to the final RMSD objective?**

Trace the complete chain:

sequence
→ candidate generation
→ candidate pool
→ candidate representation
→ structural prior
→ Hamiltonian
→ quantum state
→ CVaR objective
→ optimization
→ measurement
→ readout
→ structural reconstruction
→ geometric projection
→ physical refinement
→ built chain
→ RMSD

At every stage ask:

* What information enters?
* What information leaves?
* What information is destroyed?
* What information is compressed?
* What information is duplicated?
* What information is never represented?
* What information is target-specific?
* What information is global?
* What information is common-mode?
* What information is recoverable?
* What information is lost irreversibly?

Then determine where the loss of endpoint performance happens.

============================================================
6. PERMIT A COMPLETE ARCHITECTURAL REDESIGN
===========================================

You have broad freedom.

You may:

* redesign candidate generation
* change candidate count
* change candidate representation
* redesign the Hamiltonian
* redesign the quantum register
* redesign the ansatz
* redesign CVaR
* redesign the measurement process
* redesign the readout
* redesign the reconstruction
* introduce new physical observables
* introduce new structural coordinates
* introduce new latent variables
* introduce new thermodynamic variables
* introduce dynamics
* introduce trajectory information
* introduce additional physical simulations
* introduce target-conditioned operators
* introduce hierarchical models
* introduce multiscale formulations
* introduce hybrid quantum/classical stages
* introduce multiple quantum stages
* introduce adaptive quantum circuits
* introduce new optimization variables
* change the candidate-to-qubit mapping
* change what a basis state means
* change what the quantum state represents

You may even redesign the architecture around a completely different quantum formulation while keeping genuine CVaR-VQE as the central optimization mechanism.

Do not preserve bad architecture for historical reasons.

============================================================
7. BUT DO NOT THROW AWAY USEFUL INFRASTRUCTURE
==============================================

Broad freedom does NOT mean careless destruction.

Preserve:

* benchmark integrity
* frozen folds
* existing artifact history
* Legacy Hamiltonian
* AMBER pathway
* production baseline
* statistical infrastructure
* experiment records
* negative results
* reproducibility artifacts

You may replace their role.

You may experimentally supersede them.

You may wrap them.

But do not erase scientific history.

============================================================
8. THE REPORT'S CURRENT SCIENTIFIC STATE
========================================

Current endpoint:

~3.2105 Å built chain.

Cloud:

~3.0483 Å

These are NOT the same metric.

Set mean:

~3.5507 Å

This is also a different object.

Always label the basis.

The primary endpoint is built-chain RMSD.

The cloud is diagnostic.

============================================================
9. THE TAIL
===========

The mean is heavily influenced by the tail.

The project has repeatedly found that repairing a small number of bad targets can move the mean much more than making tiny improvements everywhere.

Therefore analyze:

* mean
* median
* full distribution
* upper quantiles
* worst targets
* FAIL18
* easy 108
* target-level error modes
* candidate availability
* candidate rank failures

But do not optimize directly on a post-hoc target-defined tail.

The final system must work on the whole instrument.

The tail is a mechanism diagnostic.

============================================================
10. THE CANDIDATE POOL
======================

S31 shows that useful candidates can exist inside bad target pools.

The problem may therefore be selection.

But do not assume selection is the only problem.

Explicitly decompose:

POOL QUALITY
vs
POOL DIVERSITY
vs
POOL COMMON-MODE BIAS
vs
WITHIN-POOL RANKING
vs
READOUT
vs
RECONSTRUCTION

Measure each independently.

Do not call a bad candidate pool a selection problem just because the best candidate exists somewhere inside it.

Do not call a selection problem a pool-generation problem merely because better candidates could hypothetically be generated.

Find which stage is actually causing the RMSD loss.

============================================================
11. COMMON MODE
===============

S30/S31 established a strong common-mode interpretation.

It may be correct.

It may be incomplete.

It may also be mathematically correct under the chosen decomposition while being insufficient to explain the ultimate RMSD bottleneck.

Investigate all three possibilities.

Ask:

* Does common-mode error really control endpoint RMSD?
* Is common-mode error stable across candidate-generation methods?
* Does it survive reconstruction?
* Does it survive projection?
* Does it survive the averaging operator?
* Does it differ between easy and hard targets?
* Is it spatially structured?
* Is it low-rank?
* Is it torsional?
* Is it long-range?
* Is it radial?
* Is it handedness-related?
* Is it dynamical?
* Is it an artifact of the current coordinate representation?
* Is it an artifact of comparing against a single reference?
* Is it actually a geometric bias rather than an information deficit?

Do not accept "common mode" as an explanation until you connect it quantitatively to RMSD.

============================================================
12. THE CRITICAL ORACLE RESULTS
===============================

Use all ORACLE results aggressively as diagnostic probes.

But label:

ORACLE
NOT DEPLOYABLE

The known oracle results are valuable because they tell you:

* the target is not obviously outside the available structural region
* useful candidate choices exist
* useful correction directions exist
* the circuit family may have much larger structural capability than the deployed objective exploits
* the tail can contain recoverable error

Use oracle information to identify:

**WHAT INFORMATION WOULD HAVE TO EXIST FOR RMSD TO DROP?**

Then work backward toward a native-free way of obtaining it.

Do not stop at "the answer exists."

Identify the minimal information required.

============================================================
13. REOPEN CLOSED DIRECTIONS WHEN JUSTIFIED
===========================================

A major change from previous prompts:

"closed" does NOT mean "forbidden forever."

A direction is closed only relative to the exact mechanism that was tested or proved.

You may reopen a direction if:

1. its previous theorem has a hidden assumption that does not hold,
2. the implementation differed from the stated mathematical object,
3. the endpoint mapping invalidated the original closure,
4. a genuinely different formulation escapes the theorem,
5. a new information source changes the problem,
6. multiple weak effects interact,
7. new code inspection reveals the previous test was not testing the intended object,
8. a stronger mathematical formulation makes the previous experiment irrelevant.

But do NOT reopen something merely because you dislike its negative result.

The burden for reopening is a concrete mechanism.

============================================================
14. CVaR-VQE IS THE CENTRAL QUANTUM OBJECT
==========================================

The main research question is:

**Can a genuine CVaR-VQE be designed whose quantum state and objective actually contain information that can improve RMSD?**

Do not assume the current candidate-index CVaR formulation is the final formulation.

Consider quantum states representing:

* candidates
* candidate subsets
* structural modes
* torsional states
* long-range pair configurations
* conformational basins
* transition states
* physical response states
* correction directions
* low-dimensional structural coordinates
* candidate mixtures
* hybrid states

The computational basis does not have to mean "candidate number."

============================================================
15. TARGET-CONDITIONAL CVaR
===========================

The current major limitation may be that the deployed Hamiltonian is effectively global.

Therefore aggressively investigate target-conditioned formulations.

For a target t:

H_t

should depend on a legitimate target-specific quantity.

Potential sources:

* sequence-conditioned physical calculations
* candidate-specific physical response
* thermodynamics
* dynamics
* external measurements
* experimentally measured quantities
* sequence-derived but genuinely target-specific predicted observables
* environment-dependent physics
* structural ensemble response
* learned target representation, PROVIDED it is trained strictly out-of-fold and does not use native leakage

But do not merely insert a target embedding into the Hamiltonian.

The target representation must have a causal reason for improving structural prediction.

============================================================
16. QUANTUM OBJECTIVE DESIGN
============================

For every serious candidate objective derive:

1. mathematical definition
2. physical interpretation
3. state interpretation
4. basis interpretation
5. Hamiltonian construction
6. CVaR meaning
7. gradient
8. optimizer landscape
9. global/local minima
10. relation to candidate quality
11. relation to structural RMSD
12. information content
13. classical equivalent
14. failure mode
15. deployment condition

If you cannot explain the mechanism, it is not ready for expensive testing.

============================================================
17. DO NOT ASSUME DIAGONAL HAMILTONIANS
=======================================

The old diagonal rank Hamiltonian is largely a classical ordering encoded into quantum notation.

You may use diagonal Hamiltonians if they are genuinely target-dependent and scientifically useful.

But investigate whether the correct Hamiltonian lives in:

* Ising form
* QUBO
* pairwise compatibility space
* graph Laplacian
* transition operator
* torsion-space operator
* free-energy operator
* response operator
* block Hamiltonian
* sparse low-rank operator
* physically coupled operator

Non-diagonal terms must have meaning.

Not just:

"add J"

or

"add hopping."

============================================================
18. DEEPLY RECONSIDER WHAT CVaR IS SELECTING
============================================

Maybe CVaR should not select candidates directly.

Maybe it should select:

* structural modes
* pair assignments
* torsion states
* candidate clusters
* conformational basins
* transition paths
* physical perturbation responses
* candidate mixtures
* correction coordinates

Investigate.

The question is not:

"How can CVaR rank these 128 candidates?"

The deeper question is:

**What discrete or low-dimensional decision is actually responsible for lower RMSD, and can CVaR-VQE solve that decision?**

============================================================
19. INVESTIGATE A QUANTUM-CENTRIC REPRESENTATION
================================================

Try to discover a representation where quantum optimization has something meaningful to optimize.

Candidate examples:

STRUCTURAL MODE BASIS:
Each basis state corresponds to a structural mode or displacement pattern.

PAIRWISE BASIS:
Qubits encode long-range pair decisions.

TORSION BASIS:
Qubits encode meaningful torsional states.

BASIN BASIS:
Basis states correspond to conformational basins.

CANDIDATE-MIXTURE BASIS:
Amplitudes encode candidate mixture coefficients.

TRANSITION BASIS:
Basis states encode transitions between conformational states.

HYBRID BASIS:
Part of the register encodes candidates, part physical response variables.

These are examples.

Invent better ones.

============================================================
20. INFORMATION BOTTLENECK ANALYSIS
===================================

Treat every representation as an information channel.

Measure:

* entropy
* mutual information
* conditional mutual information
* effective candidate count
* rank information
* bits per decision
* continuous parameter dimension
* information per qubit
* useful information per computation
* information lost at each pipeline stage

But do not let information theory become detached from RMSD.

A quantity matters only if information in it can actually be converted into lower RMSD.

============================================================
21. DIRECTION MATTERS MORE THAN CORRELATION
===========================================

S31 demonstrates an important distinction:

A predictor can capture a real and substantial component of error while still pointing in the wrong direction for RMSD improvement.

Therefore always ask:

Does the observable point toward a structure with lower RMSD?

Not merely:

Does it predict something about the error?

Use oracle directional diagnostics where appropriate.

Build native-free tests that approximate directional utility.

============================================================
22. INVESTIGATE PHYSICAL RESPONSE
=================================

Aggressively investigate observables based on how candidate structures respond to physical perturbations.

Possible perturbations:

* torsional perturbation
* restrained relaxation
* local force perturbation
* temperature changes
* solvent/environment changes
* normal modes
* elasticity
* conformational transitions
* trajectory evolution
* basin stability
* relaxation response
* energy curvature
* susceptibility-like observables

The scientific hypothesis is:

A candidate's response to a physical perturbation may contain information about whether it lies in the correct structural basin that is invisible to static geometry alone.

Test it.

Do not assume it works.

============================================================
23. DYNAMICS
============

Investigate:

* molecular dynamics
* short trajectory ensembles
* conformational covariance
* basin transitions
* state populations
* metastability
* transition rates
* autocorrelation
* dynamic modes
* transition-path structure
* free-energy landscapes

Start with cheap simulations.

Escalate only after deriving a reason the observable could add information.

============================================================
24. FREE ENERGY
===============

Investigate free-energy methods only in ways that genuinely escape previous static-energy closure arguments.

Potential directions:

* relative basin free energy
* ensemble free energy
* conformational entropy
* thermodynamic response
* transition free energy
* path-dependent free energy
* restrained/unrestrained differences
* temperature-dependent response
* ensemble reweighting

Do not simply rename another single-structure scalar "free energy."

============================================================
25. EXPERIMENTAL INFORMATION
============================

Investigate legitimate target-specific experimental observables when accessible.

Examples can include:

* NMR observables
* chemical shifts
* residual dipolar couplings
* NOE-derived information
* exchange-related observables
* FRET-like constraints
* other experimentally measured physical quantities

But aggressively audit leakage.

For every proposed experimental channel ask:

* Is it genuinely measured independently?
* Is it available at inference?
* Does it depend on the native structure?
* Is it derived from deposited coordinates?
* Does it provide new target information?
* Is it available for enough targets?
* What is its cost?
* Can it be converted into a CVaR-VQE Hamiltonian?

Reference-derived restraints are ORACLE.

============================================================
26. SEQUENCE MODELS
===================

Do not blindly ban learned models.

You may use:

* sequence models
* protein language models
* embeddings
* distograms
* structure predictors
* neural physical surrogates

provided:

* training is properly out-of-fold,
* native leakage is prevented,
* the information source is documented,
* its output has a clear role,
* it is not merely another arbitrary scalar score,
* it is evaluated against the fixed endpoint.

A learned model may be:

* candidate generator
* target encoder
* observable predictor
* Hamiltonian parameter generator
* physical surrogate
* auxiliary channel

But the final scientific question remains whether it lowers RMSD.

============================================================
27. REVISIT CANDIDATE GENERATION
================================

The candidate pool may still be part of the real problem.

You have permission to redesign:

* retrieval
* diversity
* candidate source mixture
* sequence profile search
* embedding retrieval
* structural motif retrieval
* distogram-conditioned retrieval
* torsion-conditioned retrieval
* generative candidates
* physics-generated candidates
* perturbation-generated candidates
* dynamically generated candidates
* local refinement candidates
* hierarchical candidates

But do not optimize only pool diversity.

A diverse pool can still be wrong.

Measure:

POOL QUALITY
and
POOL ERROR CORRELATION
and
IN-BAND SELECTION SKILL.

============================================================
28. GENERATIVE STRUCTURES
=========================

You are allowed to investigate generative or optimization-based candidate generation.

Potentially:

* torsion perturbations
* local structural optimization
* physics-guided proposals
* normal-mode moves
* learned generative structures
* diffusion-derived proposals
* fragment recombination
* constrained structural sampling
* candidate interpolation where physically valid

But a generated structure only matters if it improves the final RMSD after the complete built-chain endpoint.

============================================================
29. DO NOT CONFUSE "BETTER CANDIDATE" WITH "BETTER PREDICTION"
==============================================================

A candidate can be extremely good while:

* being unreachable by the selector,
* being destroyed by averaging,
* being projected badly,
* being rejected by the Hamiltonian,
* being inaccessible to the quantum state,
* being inaccessible through the readout.

Trace the entire path.

============================================================
30. READOUT REMAINS AVAILABLE IF A NEW MECHANISM EXISTS
=======================================================

Previous readout families are closed under their tested assumptions.

They may be reopened only with a new upstream information source or a genuinely different mathematical formulation.

Investigate readouts that produce:

* a candidate
* sparse candidate mixture
* weighted candidate mixture
* mode coefficients
* pair assignments
* torsion configurations
* correction vectors
* local structural edits
* quantum amplitudes mapped to structural weights

Do not assume uniform averaging is sacred.

Do not assume argmin is sacred.

Do not assume the convex reconstruction is sacred.

Change it if RMSD evidence says it is the bottleneck.

============================================================
31. PROJECTION AND RECONSTRUCTION
=================================

The built-chain conversion is part of the causal system.

Do not assume it is a harmless post-processing step.

Study:

* sensitivity
* conditioning
* branch selection
* local minima
* initialization
* multi-start
* geometry constraints
* chain continuity
* chirality
* torsions
* projection-induced displacement
* projection-induced ranking changes

Determine whether the reconstruction itself is destroying structural improvements.

A perfect cloud solution is useless if the built chain becomes worse.

============================================================
32. AMBER
=========

Keep AMBER ff14SB/GBn2 available as a real physical Hamiltonian.

Do not assume it is a good ranker.

Investigate whether AMBER is more useful for:

* response
* dynamics
* perturbation
* curvature
* physical filtering
* ensemble generation
* structural refinement

rather than static candidate ranking.

One AMBER process at a time by default.

============================================================
33. LEGACY
==========

Keep Legacy as a real independent Hamiltonian.

Investigate its relationship to:

* candidate structure
* structural modes
* CVaR landscape
* AMBER
* torsion
* long-range interactions
* compactness
* optimization geometry

Do not turn it into a renamed heuristic.

============================================================
34. CVaR + PHYSICS
==================

Consider whether CVaR can be used as a risk-sensitive physical variational principle.

Examples:

CVaR of conformational free energy.

CVaR of physical response.

CVaR over basin energies.

CVaR over transition costs.

CVaR over candidate physical stability.

CVaR over structural compatibility.

CVaR over a hybrid physical/combinatorial Hamiltonian.

Invent better formulations.

Derive them.

Test them.

============================================================
35. CVaR + ENSEMBLES
====================

Investigate whether the wavefunction should represent uncertainty over structures rather than simply selecting a candidate.

Possible idea:

quantum amplitudes represent a structured posterior over conformational hypotheses.

Then CVaR operates on physically meaningful risk.

Then the readout produces a posterior structure.

Or:

quantum state represents competing correction hypotheses.

Or:

quantum state represents structural basin occupancy.

These are examples only.

You are allowed to invent the formulation.

============================================================
36. VQE OPTIMIZER
=================

Optimizer tuning is secondary.

Do not spend the majority of the sprint on:

* Adam learning rates
* SPSA hyperparameters
* restarts
* depth
* arbitrary ansätze
* iteration counts

until the objective itself has evidence of endpoint relevance.

Once a genuinely useful Hamiltonian exists, aggressively optimize the VQE implementation.

============================================================
37. ANSATZ FREEDOM
==================

You may redesign the ansatz completely.

Potential options:

* hardware-efficient
* problem-inspired
* symmetry-preserving
* ADAPT-VQE
* qubit-ADAPT-VQE
* pair-correlated
* excitation-like operators
* graph-based operators
* physically motivated generators
* shallow structured ansatz
* tensor-network-inspired ansatz

But:

ansatz complexity is justified only if it unlocks a meaningful state that produces lower endpoint RMSD.

============================================================
38. MATHEMATICAL DEPTH
======================

Use maximum mathematical depth.

Explore:

* variational calculus
* convex analysis
* nonconvex optimization
* perturbation theory
* spectral theory
* matrix inequalities
* low-rank approximations
* information geometry
* Fubini-Study geometry
* Riemannian geometry
* Lie algebras
* statistical mechanics
* Bayesian decision theory
* rate-distortion ideas
* inverse problems
* identifiability
* causal decomposition
* sufficient statistics
* manifold learning
* operator theory
* quantum measurement theory
* variational principles

A theorem that closes 1,000 experiments is extremely valuable.

A theorem that tells you where the RMSD gain must live is even more valuable.

============================================================
39. PHYSICS DEPTH
=================

Use maximum physical reasoning.

Think about:

* energy landscapes
* frustration
* conformational entropy
* torsional coupling
* hydrogen bonding
* long-range interactions
* solvent
* electrostatics
* sterics
* chirality
* topology
* basin stability
* metastability
* kinetic trapping
* thermodynamic selection
* collective coordinates
* structural response

Do not use physics terminology merely to decorate a heuristic.

============================================================
40. INFORMATION THEORY DEPTH
============================

Use information theory only when it connects to prediction.

For every information argument ask:

"How does this ultimately lower RMSD?"

Potential quantities:

* mutual information with useful structural choices
* conditional mutual information
* information gain
* entropy reduction
* effective number of candidates
* channel capacity
* information per qubit
* information per physical simulation
* information per unit compute
* sufficient statistics

Do not infer a missing information source merely because a correlation is small.

============================================================
41. THE FIVE-BIT RESULT
=======================

Treat the previous five-bit oracle result as an information-pricing clue, not as a hard limit.

Ask:

* Are five bits really sufficient?
* Are they necessary?
* Is the binary representation responsible?
* Could a continuous low-dimensional signal be more efficient?
* Could one physical observable encode many useful bits?
* Could a quantum state compress those bits?
* Could the required information be encoded in amplitudes rather than basis states?
* Could a nonlinear observable amplify the useful information?
* Is the apparent five-bit structure an artifact of the chosen oracle construction?

Try to derive a genuine minimal-information theorem where possible.

============================================================
42. DO NOT OVERFIT TO THE CURRENT PROBLEM FORMULATION
=====================================================

The current formulation is:

protein sequence
→ retrieved candidate pool
→ score candidates
→ CVaR/VQE
→ tail
→ average
→ project
→ RMSD

But perhaps the correct formulation is:

sequence
→ target-conditioned physical representation
→ quantum state
→ CVaR
→ structural reconstruction

Or:

sequence
→ physical candidate ensemble
→ quantum basin selection
→ reconstruction

Or:

sequence
→ candidate pool
→ target-specific physical response
→ CVaR-VQE
→ sparse readout
→ final structure

Or something entirely different.

Do not constrain yourself to one pipeline graph.

============================================================
43. FINAL DEPLOYMENT CAN BE HYBRID
==================================

You may use:

classical
→ quantum
→ classical
→ quantum
→ physical
→ geometric

if that improves RMSD.

"Quantum" does not mean every operation must be quantum.

The final selector must retain a genuine CVaR-VQE role.

============================================================
44. SCIENTIFIC FALSIFICATION
============================

For every major hypothesis create:

HYPOTHESIS
MECHANISM
PREDICTION
FALSIFIER
CONTROL
TEST
STATISTICAL RULE
DEPLOYMENT CONDITION

Never wait until after seeing the outcome to decide what the experiment meant.

============================================================
45. MULTIPLE TESTS AND FALSE DISCOVERY
======================================

The project has had enough experiments that multiple testing is a serious concern.

Track:

* number of hypotheses
* number of exploratory arms
* preregistered arms
* confirmatory arms
* independent replication
* fold structure
* number of comparisons

Do not turn the best result from a massive search into a headline without accounting for search.

Independent replication is especially important for anything surprising.

============================================================
46. ORACLE DISCIPLINE
=====================

Oracle work is encouraged.

Use it to answer:

"What would be possible if the missing information were available?"

But label every oracle result:

ORACLE
NOT DEPLOYABLE

Do not accidentally allow native coordinates to influence:

* feature selection
* hyperparameters
* architecture choice
* deployment thresholds
* training
* candidate generation
* readout parameters

============================================================
47. EXPERIMENTAL CONTROL
========================

Use strong controls.

At minimum consider:

* current production
* random
* matched random
* shuffled
* magnitude-matched
* information-matched
* compute-matched
* classical control
* oracle ceiling
* no-op
* projection-only
* physics-only
* quantum-only

The exact controls should match the mechanism being tested.

============================================================
48. DO NOT CHASE P-VALUES
=========================

A tiny p-value on the wrong quantity is worthless.

A large effect with weak measurement may still be unresolved.

Use the project's MDE standard:

MDE = 2.8016 × SE

Use:

* paired target-level analysis
* fold-clustered CIs
* concentration-aware nulls
* Type-M where useful

Below 0.7× MDE:
not a result.

0.7–1.0×:
not measured.

Do not inflate weak results into wins.

============================================================
49. FAIL18
==========

Study FAIL18 deeply.

But treat it as a diagnostic mechanism.

Ask:

Why does the existing score lose within-pool skill there?

Why does the best candidate hide?

Why is the score near or below an uninformative ranking null?

Why is the recoverable prize concentrated there?

Does the new mechanism specifically restore useful information there?

Then verify that it also improves the full 126-target endpoint.

============================================================
50. LOOK FOR INTERACTIONS
=========================

Most previous work appears to have tested components separately.

Do not assume:

weak + weak = useless.

Investigate carefully selected interactions such as:

new physical observable
×
new Hamiltonian

new candidate pool
×
CVaR

new target representation
×
new VQE

new quantum basis
×
new physical score

new physical response
×
sparse readout

candidate generation
×
target-conditioned Hamiltonian

But use mechanistic reasoning first.

Do not launch combinatorial parameter sweeps blindly.

============================================================
51. COMPUTE STRATEGY
====================

Use the machine aggressively.

Target approximately:

94–95% CPU utilization
94–95% RAM utilization

while maintaining system stability.

Do not:

* intentionally swap into unusability
* spawn uncontrolled processes
* crash the machine
* corrupt artifacts
* starve Claude Code
* overwhelm Windows process management

The target is sustained high utilization with controlled headroom.

============================================================
52. FOUR AGENTS AT ALL TIMES
============================

Maintain at least 4 active Claude research agents during substantive research.

Never exceed 8.

Claude chooses:

4
5
6
7
or 8

based on the scientific workload.

Default roles:

AGENT 1:
CVaR-VQE / quantum theory

AGENT 2:
physical chemistry / structural biology / dynamics

AGENT 3:
information / candidate generation / prediction theory

AGENT 4:
experimental verification / statistics / adversarial audit

Additional agents can be created for:

* literature research
* mathematical theorem proving
* independent implementation
* adversarial replication
* computational acceleration
* visualization/diagnostics

Claude decides dynamically.

Do not create agents just for activity.

Every active agent needs a real scientific question.

============================================================
53. ADVERSARIAL SCIENCE
=======================

At least one active agent should aggressively try to disprove every promising result.

Require adversarial attacks on:

* leakage
* circularity
* target dependence
* hidden native information
* projection artifacts
* selection artifacts
* optimizer artifacts
* order-statistics
* multiple testing
* fold leakage
* compute mismatch
* metric mismatch
* bad null
* implementation mismatch

A result should survive an attempt to kill it.

============================================================
54. LITERATURE RESEARCH
=======================

Perform serious literature research where useful.

Investigate:

* CVaR-VQE
* variational quantum algorithms
* risk-sensitive optimization
* quantum optimization
* quantum encodings
* free-energy methods
* molecular dynamics
* protein conformational landscapes
* Markov state models
* enhanced sampling
* transition path theory
* torsional landscapes
* structural ensembles
* NMR observables
* protein language models
* physical response functions
* inverse problems
* information theory

Do not copy literature blindly.

Find mechanisms that matter for THIS endpoint.

============================================================
55. START WITH A CAUSAL AUDIT
=============================

Before major experiments, construct a causal map:

SOURCE
↓
REPRESENTATION
↓
INFORMATION
↓
HAMILTONIAN
↓
VQE
↓
CVaR
↓
READOUT
↓
STRUCTURE
↓
BUILT CHAIN
↓
RMSD

For each arrow:

What mathematical transformation occurs?

What information is preserved?

What information is destroyed?

What assumption is introduced?

What error can arise?

Where could a lower RMSD solution disappear?

This map is the backbone of the sprint.

============================================================
56. IDENTIFY THE EARLIEST IRREVERSIBLE LOSS
===========================================

A key research goal:

Find the earliest stage at which the information required for a lower RMSD structure becomes unavailable.

Possible locations:

* candidate generation
* candidate pool
* scoring
* Hamiltonian
* quantum encoding
* optimization
* CVaR tail
* measurement
* readout
* averaging
* projection

If an earlier stage destroys the necessary information, later optimization cannot recover it.

If later stages destroy it, the upstream pool may already be sufficient.

This distinction should determine where compute goes.

============================================================
57. FIND THE MOST LEVERAGED INTERVENTION
========================================

Once the causal map exists, estimate:

Potential RMSD gain
/
compute cost
/
information required
/
experimental uncertainty

Use this only for prioritization.

Do not turn it into a fake objective.

The intervention with the largest expected scientific leverage should get the most resources.

============================================================
58. FINAL GOAL OVERRIDES PREVIOUS CHARTER WORDING
=================================================

The original charter and reports are historical constraints for scientific integrity.

They are NOT commands to preserve every design decision.

The final goal is:

**produce lower-RMSD structures.**

The CVaR-VQE spine is important because this is a quantum research project.

But do not force a useless quantum component to remain in a place where it prevents the scientific problem from being solved.

Instead:

find the formulation in which CVaR-VQE is genuinely useful.

============================================================
59. SCIENTIFIC BREAKTHROUGH THRESHOLD
=====================================

A genuine breakthrough should ideally satisfy:

* mean built-chain RMSD improves by >= 1 MDE
* result survives paired statistics
* at least 4/5 folds agree where applicable
* mechanism survives adversarial audit
* no native leakage
* reproducible
* deployable
* explanation connects intervention to RMSD

A result below MDE is not a breakthrough.

An oracle result is not a breakthrough.

A training metric is not a breakthrough.

============================================================
60. NEGATIVE RESULTS
====================

Negative results are valuable when they identify mechanism.

Do not report:

"didn't work."

Report:

"didn't work because..."

Examples:

* no new information
* wrong target dependence
* wrong direction
* representation mismatch
* insufficient quantum capacity
* classical equivalence
* readout loss
* reconstruction loss
* projection loss
* physical irrelevance
* statistics insufficient
* theorem-closed
* computationally impossible

============================================================
61. REPRODUCIBILITY
===================

Fix reproducibility problems that affect interpretation.

In particular:

* pin projection/build-chain behavior
* verify cloud versus chain
* verify seeds
* verify branch state
* verify folds
* verify target set
* verify artifacts
* verify experiment completeness
* verify verifier correctness

Do not claim tiny improvements beyond the measurement floor.

============================================================
62. BENCHMARK INTEGRITY
=======================

Never open or alter benchmark60.

Never regenerate the frozen folds.

Never tune deployable parameters using target-native RMSD.

Native RMSD is for evaluation and explicitly marked oracle diagnostics only.

============================================================
63. ENGINEERING FREEDOM
=======================

You may refactor the repo substantially.

You may:

* reorganize modules
* remove obsolete code
* create cleaner abstractions
* build experiment registries
* create reusable quantum objects
* create physical-observable interfaces
* improve the scheduler
* improve artifact tracking
* improve verification

Do not destroy scientific history.

============================================================
64. RESOURCE-SAFE AUTONOMY
==========================

Do not ask the user for routine confirmation.

You decide:

* experiment ordering
* agent allocation
* compute allocation
* architecture
* stopping
* continuation
* reopening
* closure
* deployment candidate

Pause only for genuinely destructive external actions that require user authorization.

============================================================
65. WHAT YOU SHOULD NOT DO
==========================

Do not spend the sprint blindly repeating:

* Gray code
* arbitrary index permutation
* optimizer roulette
* deeper circuits without a useful objective
* more iterations without objective evidence
* arbitrary hopping terms
* old sparse-readout variants
* old single-structure achiral channels
* the same short torsion scalar families
* another generic native-free residual regressor
* another global scalar cost
* another cosmetically different rank transform

unless you have identified a materially different mechanism.

============================================================
66. WHAT YOU SHOULD DO IF YOU FIND SOMETHING SURPRISING
=======================================================

Stop and understand it.

Do not immediately build a deployment around it.

Ask:

* Is it real?
* Is it reproducible?
* Is it causal?
* Is it native-free?
* Is it target-specific?
* Is it RMSD-relevant?
* Does the effect survive controls?
* Is the effect concentrated in one fold?
* Is it a projection artifact?
* Is it an oracle artifact?
* Is it a metric artifact?
* Is it actually new?

Then decide whether to scale.

============================================================
67. THE ULTIMATE RESEARCH QUESTION
==================================

Do not assume this sentence has already been answered:

> **Why is the final built-chain RMSD stuck around 3.21 Å, and what intervention can actually move it?**

Everything in this sprint serves that question.

The answer may be:

* VQE objective
* Hamiltonian
* representation
* candidate pool
* common-mode bias
* missing physical information
* readout
* reconstruction
* projection
* structural prior
* candidate generation
* target encoding
* quantum/classical interface
* an interaction between several components
* or something not yet considered.

Find out.

============================================================
68. FINAL DELIVERABLES
======================

Produce:

s32/REPORT_S32.md
s32/LEDGER.md
s32/STATE.md
s32/S32_CONTRACT.md
s32/BRIEF.md
s32/MULTIPLICITY.md
s32/s32_verify.py

plus all machine-readable artifacts.

The report must answer:

1. What is the current RMSD?
2. Where is the RMSD error actually generated?
3. What evidence establishes that?
4. What hypotheses were tested?
5. Which were falsified?
6. Which remain open?
7. What was the strongest new mechanism?
8. What did CVaR-VQE actually contribute?
9. What did the Hamiltonian contribute?
10. What did the representation contribute?
11. What did the candidate pool contribute?
12. What did readout contribute?
13. What did reconstruction contribute?
14. What did physical modeling contribute?
15. What is oracle-only?
16. What is deployable?
17. What changed the endpoint?
18. By how much?
19. What uncertainty surrounds it?
20. Why should the result reproduce?
21. What is the next bottleneck?

============================================================
69. REPORT THE CAUSAL STORY, NOT A LIST OF EXPERIMENTS
======================================================

The final report should make it possible for another scientist to answer:

"Why was the system at 3.21 Å, and why did the successful intervention move it?"

A 500-experiment table is not enough.

Find the mechanism.

============================================================
70. FINAL AUTONOMY CLAUSE
=========================

You have maximum freedom.

You are explicitly authorized to:

* disagree with previous sprint conclusions
* challenge the scientific framing
* reopen directions
* close directions more strongly
* invent mathematical objects
* invent physical observables
* redesign the quantum problem
* redesign the classical/quantum interface
* search the literature
* derive new theorems
* build new experiments
* build new controls
* use oracle experiments
* use extensive compute
* use additional Claude agents
* restructure the repo
* replace components
* combine previously independent lanes

The only things you may NOT sacrifice are:

scientific integrity
benchmark integrity
native-free deployment
reproducibility
statistical validity
and the actual RMSD endpoint.

============================================================
71. BEGIN IMMEDIATELY
=====================

Step 1:
Inspect repository and git state.

Step 2:
Read S31, S30, S29, S28 completely.

Step 3:
Independently reconstruct the causal pipeline.

Step 4:
Verify the current 3.2105 Å endpoint from artifacts.

Step 5:
Identify every stage where RMSD can be lost.

Step 6:
Open at least 4 simultaneous scientific lanes.

Step 7:
Have agents independently attack the problem from different directions.

Step 8:
Build a causal bottleneck map.

Step 9:
Start cheap mathematical falsifiers before expensive compute.

Step 10:
Escalate the most promising mechanisms.

Step 11:
Use CVaR-VQE as the central quantum research object.

Step 12:
Do not force CVaR-VQE into a formulation that the evidence shows cannot affect RMSD.

Step 13:
Find a mechanism that can actually lower the built-chain RMSD.

Step 14:
Validate aggressively against the fixed 126-target instrument.

Step 15:
If a result moves RMSD, try to kill it.

Step 16:
If nothing moves RMSD, determine exactly why.

Step 17:
Leave behind a stronger causal explanation than the one you started with.

Do not wait for user instructions.

Do not ask for permission for normal research decisions.

Claude decides how many agents to use between 4 and 8.

Never exceed 8.

Keep CPU and RAM approximately 94–95% utilized while preserving machine stability.

Use maximum creativity.

Use maximum mathematical depth.

Use maximum physics depth.

Use maximum quantum depth.

Use maximum structural-biology depth.

Use maximum computational experimentation.

Use maximum skepticism.

And above all:

**THE FINAL SCORE IS RMSD.**

**FIND THE PROBLEM.**

**THEN FIX IT.**

BEGIN.
 final report after everything is done. Say if its possible to lower rsmd at the end, maybe test on longer proteins etc. I BELIEVE IN YOU, you got this!

<!-- ============================ END VERBATIM CHARTER ============================ -->
