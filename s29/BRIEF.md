# SPRINT 29 — CVaR-VQE PROTEIN FOLDING: FIND THE PATH BELOW 2.5 Å

(The user's charter, verbatim, received 2026-09-19. This file is the authority for S29; `s29/S29_CONTRACT.md` is its operational form.)

## 0. HOW TO READ THIS PROMPT

This is not a task list. It is a research charter.

There is exactly one hard constraint (CVaR-VQE remains the central scientific object) and exactly one primary endpoint (mean built-chain Cα RMSD on the 126-target instrument). Everything else in this document — every architecture, encoding, Hamiltonian, objective, readout, and representation currently in the repository — is open for replacement.

Do not treat the sections below as sequential instructions. Read the evidence first, form your own scientific position about where the largest remaining leverage lives, and then attack that position with everything you have. If the strongest thing you can do contradicts a suggestion in this prompt, do the strongest thing and explain why.

This sprint is expected to take significantly longer than previous sprints. Depth is the point. Do not converge early. Do not ship the first plausible idea. A sprint that spends most of its budget on reading, deriving, and falsifying before committing to a build is a correctly run sprint.

## 1. READ FIRST, IN ORDER

Before doing any work:

1. Read `s27/REPORT_S28.md` end to end.
2. Read `s27/LEDGER.md` through the S28 close.
3. Read `s27/RESUME_S28.md`.
4. Read the relevant S28 findings and retractions, including the retractions specifically — a retracted claim tells you where the project's instincts were wrong, which is more informative than a confirmed one.
5. Inspect the current repository architecture and every artifact needed to understand the deployed CVaR-VQE end to end: candidate generation, pool construction, the posterior, the Hamiltonian build, the encoding, the ansatz, the objective, the optimizer configuration, the readout, structure construction, geometry handling, and the evaluation harness.
6. Reconstruct, in your own words, the full data path from sequence to final built chain, and write down every place where information is created, transformed, compressed, or destroyed. You will need this map to reason about where the missing information lives.
7. Identify every place in the current system where a modeling choice was made for convenience rather than for a stated scientific reason. These are prime redesign targets.
8. Search the scientific literature deeply before committing to a direction.

The S28 report is the current state of the science. Do not re-derive closed experiments unnecessarily. Build past them.

But also: do not inherit S28's framing uncritically. Prior sprints closed questions under specific formulations. A closure is a statement about a formulation, not about the universe.

## 2. PROJECT GOAL

The main objective of this sprint is now extremely clear:

> **Reduce the real built-chain mean Cα RMSD as much as scientifically possible, with the CVaR-VQE remaining the central component of the system.**

Current built-chain production is about **3.21 Å** on the 126-target instrument.

The aspirational target is: **< 2.5 Å mean built-chain RMSD**

- A result below **3.0 Å** would already be meaningful movement.
- A result below **2.5 Å** would be a major result.
- A sub-**2 Å** result would be exceptional, but do not optimize around a number you cannot honestly support.

Do not assume the target is reachable from the present candidate pool, the present objective, the present representation, or the present information content of the system. Determine experimentally what is possible, and if the answer is "not from here," determine precisely what must change and whether that change is achievable inside this sprint.

An honest, well-evidenced "this ceiling is imposed by X, and here is the measurement proving it" is worth more than a fragile 2.9 Å.

## 3. CORE DESIGN PRINCIPLE

The CVaR-VQE is the spine.

Do not redesign the project into a classical model with a token quantum component bolted on at the end. Do not produce a system where the quantum stage could be deleted without changing the output.

**There is deliberately no prescribed architecture in this prompt.**

Previous sprints supplied a pipeline sketch. That sketch is withdrawn — not because it was wrong, but because it has become an anchor. The shape of the system is now part of the research question. You are to decide what the architecture should be, including whether it should resemble anything the project has built so far.

The only central constraint is: **CVaR-VQE remains the main scientific object, not a decorative final stage.**

What that means operationally:
- There is a parameterized quantum state.
- There is a Hamiltonian or set of observables with a defensible meaning.
- There is a CVaR (or rigorously justified CVaR-family) objective over a measured quantity.
- The optimization of that objective is causally connected to the quality of the final structure.
- Removing or randomizing the quantum stage measurably degrades the result.

Everything else — what the qubits mean, what a basis state is, what the Hamiltonian encodes, how many quantum stages exist, whether the quantum stage runs once or iteratively, whether it selects, samples, scores, generates, constrains, or refines — is yours to determine.

## 4. SCIENTIFIC STARTING POINT FROM S28

Internalize these findings before choosing the next direction. For each one, also ask the second-order question in italics.

1. **The current deployed circuit family can express structures near 0.25 Å mean RMSD in ORACLE diagnostics.** Circuit expressivity is not the demonstrated bottleneck. *Does expressivity-in-principle imply reachability under a trainable objective? Expressible is not the same as findable.*
2. **The optimizer can reduce the current objective.** Optimization of the current objective is not the demonstrated bottleneck. *Is the optimizer reducing the objective along directions that matter, or only along cheap ones?*
3. **The current objective is the demonstrated bottleneck.** *Is it the objective's functional form, its inputs, its normalization, or the information available to construct it?*
4. **Along the important structural ladder, reducing the shipped cost can correspond to increasing RMSD.** S28 measured Spearman ≈ −0.40. *An anti-correlated objective is not merely uninformative — it is actively adversarial. Is the sign stable across the ladder, or does it flip regionally?*
5. **The local gradient of the shipped objective at the production point is essentially uninformative about the native direction.** Mean cosine ≈ −0.03, random-direction reference ≈ 0.14. *A gradient worse than random directions is a strong statement. What structure in the objective produces this?*
6. **Signed-amplitude freedom did not solve the problem** because the current objective actively drives the emitted structure away from the native. *Additional representational freedom cannot repair a misaimed objective. What class of fix can?*
7. **The first non-diagonal Hamiltonian did not solve the problem.** Its Gaussian similarity graph was nearly rank one, and the spread-spectrum follow-up altered gradient scaling without producing an accuracy result. *Rank-one similarity means the off-diagonal structure carried almost no information. Was the failure "non-diagonal Hamiltonians don't help" or "that particular similarity measure was degenerate"? These are very different closures.*
8. **The current native-free scorer library does not reliably recognize the extremely good ORACLE structures** that the circuit family can express. *This may be the deepest finding in the list. If nothing available can recognize a good structure, no amount of search or expressivity will help. Recognition may be the true bottleneck behind the objective bottleneck.*
9. **Alternative Hamiltonians alone have not solved the problem.** *Hamiltonian variation within a fixed encoding and fixed information set may be a low-dimensional search over the wrong space.*
10. **The candidate pool still contains far better structures than production.** ORACLE top-75 ceiling ≈ 2.31 Å; pool-best ORACLE ceiling ≈ 1.71 Å. These are ORACLE diagnostics, not achieved results. *There is real headroom inside the existing pool. The gap between 3.21 and 1.71 is a recognition and aggregation gap, not necessarily a generation gap.*
11. **The pool has substantial common-mode error, around 68 percent.** *Consensus and averaging methods inherit common-mode error rather than cancelling it. Any method whose power comes from agreement across candidates is structurally limited by this number.*
12. **Therefore the central unresolved problem is not merely "search harder."** It is: **How can the system construct and optimize an objective whose reduction actually corresponds to structural improvement, and how can CVaR-VQE exploit that objective without discarding the useful information?**
13. **The ideal architecture must eventually connect:** better structural information → a meaningful Hamiltonian → genuine CVaR-VQE behavior → useful quantum output → improved final built-chain RMSD. *Do not assume these components have to be connected in the same way they are today, in the same order, in the same number, or as a single pass.*

Before you commit to a direction, state explicitly which of findings 1–11 your hypothesis attacks, and which ones it accepts as binding. A hypothesis that does not clearly engage the recognition failure (8) or the common-mode error (11) should have a stated reason for ignoring them.

## 5. RESEARCH FREEDOM

You have broad freedom. Use it. Do not interpret this prompt as a list of experiments to run one by one. Instead: **Read the evidence, form your own scientific hypothesis about where the largest remaining leverage is, and attack it aggressively.**

You have full authority to redesign any of the following, individually or jointly: the candidate representation; the candidate pool and how it is generated; whether a candidate pool is the right state space at all; the Hamiltonian, or the entire notion of a single Hamiltonian; the encoding and what a qubit corresponds to; the meaning of a basis state; the ansatz, including adaptive and problem-derived ansätze; the CVaR objective, its tail parameter, and what quantity it acts on; entropy, temperature, and free-energy terms; the observables and what is measured; the meaning of the quantum amplitudes; the readout and how the quantum state is converted to structure; the structural representation itself; how distogram or other structural information enters, and at what stage; how physics enters, and whether it should enter as a term, a constraint, a prior, or a projection; how candidates interact, or whether "candidates" persist as a concept; whether the VQE acts over candidates, structural modes, fragments, torsions, constraints, pair relations, graph states, latent coordinates, basins, or something else entirely; whether multiple Hamiltonian components should act simultaneously, hierarchically, or at different stages; whether a multi-stage, iterative, or adaptive CVaR-VQE architecture is better; whether candidate generation and the VQE should be coupled rather than sequential; whether the system should run once or as a loop; whether an entirely different quantum-compatible formulation is required.

You may reject the current binary candidate-index encoding outright. You may replace it with something that has no notion of candidate identity. You may build a non-diagonal Hamiltonian whose off-diagonal structure has genuine physical, geometric, or statistical meaning rather than generic similarity. You may redesign the meaning of a basis state. You may make the quantum state represent a distribution over structures, a superposition of structural hypotheses, or a coherent basin rather than a ranking over candidate IDs. You may make the Hamiltonian encode pairwise compatibility, structural consistency, transition cost, geometric coherence, constraint satisfaction, residual error structure, or another principled relationship. You may use multiple Hamiltonian terms where each has a defensible scientific reason. You may introduce new information sources into the system if they are legitimate and native-free.

You may investigate whether CVaR should act on: candidate energies; structural observables; a composite operator; a constrained energy; a free-energy landscape; multiple correlated observables; a risk measure over a posterior predictive; another mathematically justified quantity. You may investigate whether the entropy term should remain, change form, be replaced by an explicit temperature, or interact differently with CVaR. You may investigate whether the quantum output should be sampled, read through amplitudes, converted into a subset, clustered, interpreted as modes, used to build a coherent structural ensemble, used to select a basin, used as a search distribution, or combined with a classical estimator.

Do not make any of these choices in advance, and do not inherit them silently from the existing code. Every inherited choice that survives into your final design should survive because you examined it and kept it, and you should be able to say why.

## 6. DEEP LITERATURE RESEARCH

Do serious research before settling on the final design. This sprint expects substantially more literature work than previous ones.

Search deeply across at least the following areas where relevant:

**Quantum algorithms and variational methods:** CVaR-VQE and risk-sensitive variational objectives; variational quantum algorithms generally; non-diagonal and interacting Hamiltonians in variational optimization; quantum optimization and quantum-assisted combinatorial optimization; variational quantum simulation; QAOA/VQE relationships, and what QAOA's structure implies about problem encoding; quantum Boltzmann machines, Gibbs-state preparation, and variational thermal states; quantum optimal control and its relationship to ansatz reachability; adaptive ansätze, ADAPT-VQE, qubit-ADAPT-VQE, operator-pool selection; dynamical Lie algebras and controllability; barren plateaus, gradient scaling, cost-function-dependent plateaus, locality of observables; entanglement, expressibility, and trainability trade-offs; quantum landscape geometry, Fisher information, and natural-gradient methods; measurement strategies, shot noise, and estimator variance for tail objectives; classical simulability boundaries and what distinguishes a genuinely quantum advantage claim.

**Structure, physics, and statistics:** statistical protein potentials and knowledge-based scoring functions; coarse-grained protein energy functions; backbone geometry potentials and torsional statistics; distance-based and distogram-driven structure prediction; structural ensembles and ensemble scoring; model quality assessment / native-free structure scoring (this is directly relevant to finding 8); multimodal posterior decision theory; Bayesian risk, proper scoring rules, and decision-theoretic point estimation; uncertainty calibration and miscalibration diagnostics; energy-based models and their training pathologies; graph-based molecular representations and geometric deep learning formulations; protein fragment assembly and fragment libraries; consensus, medoid, and centroid methods, and their failure under correlated error; structure-conditioned selection and reranking; physics-informed structure prediction; common-mode / correlated error in ensembles, and methods that break it; residual modeling and systematic-bias correction.

Prefer primary papers and mathematically explicit formulations. Read the equations, not the abstracts. Where a paper's claim depends on a specific assumption, note whether that assumption holds here. Do not copy a method because it is fashionable.

For every method you consider importing, answer: **"What information does this method contain that the current project does not?"** If the answer is "none, but it is formulated more elegantly," that is not a reason to introduce complexity. If the answer is "it contains a source of signal orthogonal to what we have," that is exactly the reason to introduce it.

Keep a short written record of what you read and what you rejected, including why. A negative literature result ("this entire family assumes X, which our pool violates") is a real finding and should survive into the report.

## 7. MAXIMUM SCIENTIFIC DEPTH

Think at the level of: mathematical formulation; objective geometry and conditioning; optimization landscape structure; operator structure and commutation relations; spectral properties and gap structure; Lie algebra structure and reachable sets; representability and expressivity; information content and information destruction; statistical decision theory; structural biology and physical plausibility; computational complexity and scaling.

Derive important relationships when useful. Where a derivation is short, do it explicitly. Where it is long, state the result and the argument sketch precisely enough to be checked.

Analyze whether every proposed objective is: well-defined; differentiable or subdifferentiable, and where it is not; stable under resampling and shot noise; properly normalized and scale-invariant where it should be; compatible with CVaR's order-statistic structure; identifiable; meaningful under finite candidate pools; likely to preserve or destroy multimodality; vulnerable to order-statistic artifacts; vulnerable to entropy-induced trivial solutions; vulnerable to contraction artifacts; vulnerable to common-mode error; vulnerable to degenerate spectra or rank collapse (see S28 finding 7); monotone with respect to anything you actually care about.

Do not use sophisticated mathematics as decoration. Every mathematical idea should connect to a measurable failure mode or a potential accuracy improvement. If you cannot name the failure mode a piece of math addresses, cut it. Conversely: do not avoid hard mathematics when the problem calls for it.

## 8. AGENT ORCHESTRATION

Use **4 active agents at all times whenever computationally possible.** You may use **up to 8 agents concurrently.** Never exceed 8. Never intentionally drop below 4 while meaningful work is available.

You decide how many agents are needed, which questions they attack, which implement, which derive theory, which research literature, which run experiments, which reproduce results, which act as adversaries, and when an agent should be reassigned. You are responsible for coordination, for keeping agents from duplicating work, and for integrating their findings into a single coherent scientific position.

Standing requirements:
- **Do not force every agent onto the same hypothesis.**
- **At least one active agent must always be trying to falsify or destroy the strongest current hypothesis.** This is a permanent role, not a phase.
- **At least one agent must always be exploring a substantially different direction** rather than incrementally modifying the leading design.
- **At least one agent should be reading literature for most of the sprint**, not just at the start.
- Reuse finished agents on new work immediately. Do not wait for every agent to finish before making progress. Rotate the adversary role.

Maintain a running written state of: current leading hypothesis, current strongest objection to it, what would falsify it, what each agent is doing, and what has been closed. Update this as results arrive.

## 9. EXPERIMENTAL FREEDOM

Run whatever experiments are necessary to distinguish: objective failure; recognition failure; representation failure; Hamiltonian failure; encoding failure; ansatz failure; optimization failure; information failure; readout failure; candidate-pool failure; aggregation failure; physical-validity failure; evaluation-harness failure (yes, check this too).

Do not assume the bottleneck from prior sprints is immutable. Try to falsify it under a genuinely new formulation. However, do not blindly reopen closed questions. A reopened question must have: a genuinely new formulation; an explicitly stated reason the old closure may not apply; a preregistered falsifier; a meaningful control. No cosmetic reruns.

Design diagnostic experiments that are cheap and decisive before expensive experiments that are broad and ambiguous. Where possible, build **measurements of mechanism**, not just measurements of outcome.

## 10. ACCURACY IS THE PRIMARY OUTCOME

The main scientific endpoint is: **mean built-chain Cα RMSD on the 126-target instrument.** Point-cloud RMSD is diagnostic only.

Track: mean; SE; MDE; fold-clustered CI; paired differences; per-target effects; FAIL18 / remaining 108; concentration; W/L as descriptive only; relevant controls; computational cost.

A lower objective value is not success unless it produces a lower structural RMSD. A lower energy is not success unless it improves the structure. A better gradient is not success unless it translates into meaningful optimization or structure. A more expressive circuit is not success unless it helps recognition or the final structure. A better correlation is not success unless it survives to the built chain. A beautiful theoretical result is not an accuracy result. Report the built-chain number for every serious configuration, including the ones that failed.

## 11. QUANTUM RESULT REQUIREMENT

For every serious quantum formulation, separately establish: (1) what the basis states mean; (2) what the Hamiltonian means; (3) why the Hamiltonian should correlate with useful structural information; (4) what the CVaR objective is optimizing, and over what distribution; (5) what the ansatz can represent, and what it provably cannot; (6) whether the optimizer can actually reach the relevant states; (7) whether the quantum output contains information unavailable to the classical baseline; (8) whether a classical control can reproduce the effect; (9) whether the result changes actual built-chain RMSD.

If a quantum formulation produces an accuracy improvement, test it against: a classical equivalent; a diagonalized equivalent; a permuted control; a random control; a matched-budget control; an untrained circuit; a simpler ansatz; appropriate order-statistic controls; a product-state / separable-state restriction.

Do not call something quantum merely because a quantum circuit produced it. If a classical formulation turns out to be required to make CVaR useful, preserve the CVaR-VQE spine and study precisely what the quantum component contributes — and say honestly if the answer is "little." A clean, well-controlled demonstration that the quantum stage contributes a specific, measurable, classically-irreproducible piece of information is a major result on its own, even at unchanged RMSD.

## 12. POSSIBLE HIGH-UPSIDE QUESTIONS

Prompts for creativity, not required experiments: a coherent structural basin instead of a candidate ID; a non-diagonal Hamiltonian encoding structural compatibility that avoids rank collapse; off-diagonal terms as transitions, overlaps, constraints, coherence; the Hamiltonian encoding the disagreement between the prior and the pool; preserving multimodality the average destroys; CVaR over a structural observable; a free-energy objective balancing correctness and diversity without contraction; a proper uncertainty-aware structural objective from the posterior; learning the systematic residual independently enough to break the 68% common-mode error; a new objective from off-pool structures able to recognise the ORACLE family; "which structural configuration?" instead of "which candidate?"; coupling VQE and generation; the quantum state as a search distribution; ADAPT under a genuinely correlated ground state; protein geometry on a graph with a non-classical operator algebra; hierarchical Hamiltonian terms; CVaR over basins; fixing the gradient locally without fixing the ranking; exploiting the anti-correlation; phase carrying structural information; partitioning the 126 by regime with the quantum stage deciding.

## 13. NO ARTIFICIAL CONSTRAINT ON ARCHITECTURE

Architecture redesign is allowed and expected. Constraints: do not delete historical experiments; do not erase negative results; do not silently replace old evidence with new; keep the production implementation and historical records recoverable; if a new design is promising, create a clean new experimental architecture rather than contaminating the historical baseline; any new architecture must be reproducible from its own artifacts, with seeds and configuration recorded.

## 14. STATISTICAL DISCIPLINE

Every major experiment must have its falsifier stated before running. Use the project's rules: paired per-target analysis; comparison-specific MDE; fold-clustered CIs; iid interval secondary; concentration checks; permutation controls; matched random controls; order-statistic controls; leakage poisoning; implementation audits; replication of positives. Below 0.7x MDE is not a result; 0.7 to 1x is not a demonstrated improvement; do not turn suggestive noise into a mechanism claim; replicate strong positives. Do not spend the sealed benchmark; do not regenerate `peptide_folds.json` or `peptide_clusters.json`; do not use native RMSD to tune deployable parameters; ORACLE analyses only when labeled and kept separate. Track multiplicity across the sprint.

## 15. RESOURCE MANAGEMENT

The machine is memory constrained. Respect the governor and job scheduling; checkpoint; never lose a long run; do not recompute completed work; small preregistered probe before scaling to 126; a cheap probe is not evidence for the full instrument. Budget deliberately between reading and theory, probes, and full runs. Keep RAM and CPU usage at 94-95 percent; try not to go over 95.

## 16. ADVERSARIAL SCIENCE

For every promising result: is it an order statistic; concentrated in a few targets; a tie-break; a change in entropy; contraction; a projection artifact; effective temperature; a pool change; simply classical; leakage; does it disappear under permutation, at matched budget; does it survive both seeds, fold clustering, a simpler control, the built chain; is the mechanism independently measured; would I believe this from someone else with this evidence? At least one agent must actively try to kill every serious positive. A positive that has not survived a genuine attempt to destroy it is provisional.

## 17. SUCCESS CRITERIA

Ideal: a genuine CVaR-VQE-centred architecture with a meaningful Hamiltonian and useful structural information, substantially below 3.21 Å, ideally below 2.5 Å. Strong partial: below 3.0 Å with a clear causal explanation. Strong scientific success without accuracy: a new CVaR-VQE formulation where the quantum state contributes information a matched classical method cannot reproduce, with a clear mechanism. Strong negative: a rigorous demonstration that a major route cannot work. Also valuable: a decisive measurement of what imposes the ceiling (recognition, information, common-mode error, representation).

## 18. THE STANDARD TO AIM FOR

Think beyond incremental tuning, another Hamiltonian, another optimizer, more qubits. Think about the mathematical object being optimized; what the quantum state represents; where the missing information lives; whether the encoding makes the quantum problem artificially classical; whether the pool is the wrong state space; whether CVaR is on the wrong observable; whether the readout destroys the state's structure; whether the posterior is miscalibrated or biased; whether a residual prior can break common-mode error; whether the Hamiltonian can encode relationships; whether a correlated quantum state could represent a family of compatible structures; whether the problem as posed is the problem to solve. Take the time this deserves. Explore, test, falsify, redesign, test again, then once more.

## 19. FINAL REPORT

Cover: what was tested; what changed; what failed; what succeeded; why, at the level of mechanism; what the VQE contributed, with controls; the final built-chain RMSD with full statistics; uncertainty, MDE, fold CI, concentration; every control and what it ruled out; the literature relied on and rejected; the best architecture and why; what remains unresolved; the next question and the evidence that would settle it; every hypothesis entertained and killed, with the reason. No hand-waving, no unsupported claims, no hidden tuning, no deletion of inconvenient results. Final report with artifact after everything is finished.
