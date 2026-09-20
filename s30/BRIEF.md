<!--
PROVENANCE NOTE -- added 2026-09-20 13:52 by the coordinator; NOT part of the charter.

S30-L0 (2026-09-20 12:35) states "Saved verbatim as `s30/BRIEF.md`". THE FILE WAS NEVER
WRITTEN. The claim stood unchecked for 77 minutes while eight lanes worked against a charter
that had no canonical on-disk copy.

This is the FIFTH instance in this project of prose naming a path that does not exist -- and
S30-L0 is the entry that RECORDED THE OTHER FOUR, in the same paragraph, as a reason to check
rather than assume. Writing the sentence was mistaken for doing the thing.

Everything below this comment is the charter EXACTLY as the user sent it, recovered from the
session transcript at
  ~/.claude/projects/C--Users-abena-Protein-Folding-Algorithm/2b059c1c-...-46bf.jsonl
line 13069 (39,834 characters). It is the message itself, not a reconstruction from the
paraphrases in LEDGER.md, STATE.md or S30_CONTRACT.md. Nothing has been edited, reordered or
elided; this comment is the only addition and is enclosed so it cannot be mistaken for charter
text. S30-L0 is annotated in place per contract rule 15, with its original wording left standing.
-->

# SPRINT 30 — FIND THE FIRST REAL ACCURACY BREAKTHROUGH

---

## 0. HOW TO READ THIS PROMPT

This is a research charter, not a work order.

There is one hard conceptual constraint (CVaR-VQE remains the spine and main scientific object) and one primary endpoint (mean built-chain Cα RMSD on the 126-target instrument). Everything else — every architecture, encoding, Hamiltonian, objective, readout, representation, and every recommendation listed below — is open for replacement.

Section 8 contains twelve leads carried forward from S29. **They are a leads register, not a task list.** They exist so that accumulated evidence is not lost, not to constrain you. You may pursue some, reorder them, merge them, split them, or reject them outright. If your reading of the evidence says the largest leverage is somewhere none of them point, go there and say why. A sprint that ignores ten of the twelve and finds the mechanism is a success. A sprint that dutifully executes all twelve and finds nothing is not.

This sprint is expected to run long. Depth is the point. Do not converge early. Spend real time reading, deriving, and falsifying before committing to a build, and expect to return to that phase at least once mid-sprint when results come in.

---

## 1. STARTING STATE

Read these first, completely and in this order:

1. `s27/REPORT_S28.md`
2. `s27/REPORT_S29.md`
3. `s27/LEDGER.md` through the S29 close
4. `s27/RETRACTIONS_S28.md`
5. the current architecture and all S29 artifacts needed to reproduce the final state

Treat S29 as the current state of the science.

- Do not restart old experiments from scratch.
- Do not recompute something merely because a report mentions it.
- Use existing artifacts wherever they already answer a question.

Two additional reading tasks before you start building:

6. **Reconstruct the full data path** from sequence to final built chain in your own words, and annotate every point where information is created, transformed, compressed, thrown away, or silently assumed. You will need this map constantly. Where did the 1.44 usable bits go, and where did the other 5.56 get spent?

7. **Identify every modeling choice in the current system that exists for convenience rather than for a stated scientific reason.** Inherited conveniences are the most likely hiding place for the bottleneck, precisely because nobody has had to defend them.

Read the retractions carefully. A retracted claim tells you where the project's instincts were wrong, which is more informative than a confirmed one.

---

## 2. PROJECT OBJECTIVE

The main objective is now:

> **Produce a genuinely better protein structure using a CVaR-VQE-centered architecture, while driving the mean built-chain Cα RMSD as low as possible.**

Current mean built-chain RMSD is approximately **3.21 Å** on the 126-target instrument.

Targets:

- primary practical target: **< 3.0 Å**
- ambitious target: **< 2.5 Å**
- stretch target: substantially below 2.5 Å

Do not optimize around the target number itself. Optimize around discovering a mechanism that actually reduces RMSD.

The CVaR-VQE must remain the **main scientific component**. You may redesign almost everything around it, and a good deal of what is inside it.

Do not assume the target is reachable from the present pool, objective, encoding, or information content. Determine experimentally what is possible. If the answer is "not from here," determine precisely what must change, measure it, and say so with evidence. A well-evidenced ceiling argument that redirects the next three sprints is worth more than a fragile 2.98 Å.

---

## 3. CORE CONCLUSION FROM S29

Internalize this before choosing experiments. For each finding, the italic line is the second-order question you should also be holding.

S29 established:

- the deployed circuit family can express structures around **0.25 Å** in ORACLE diagnostics
  *Expressible is not findable. Does expressivity-in-principle survive contact with a trainable objective?*
- the optimizer can optimize its objective
  *Along directions that matter, or only along cheap ones?*
- **the current objective is the barrier**
- the current objective's cost can decrease while RMSD increases
  *An anti-correlated objective is adversarial, not merely uninformative. Is the sign stable across the ladder, or does it flip regionally? A regional flip is exploitable; a global one is not.*
- the production-point objective gradient contains essentially no useful native direction
- the best displacement fields measured so far are at or below random-direction information
  *Worse than random is a strong structural statement about the objective, not a weak one. What produces it?*
- better ranking alone is not enough because the terminal operator consumes the ranking inefficiently
  *Two failures in series. Fixing either alone may show nothing. Are you measuring them separately?*
- the pool itself contains dramatically better structures
- the current 7-bit/top-128 information budget is poorly spent
- the current diagonal candidate-index Hamiltonian is too classically structured
- the first non-diagonal Hamiltonian was not enough because its graph spectrum made the hopping term effectively a typicality projector
  *Was that a closure on non-diagonal Hamiltonians, or on that similarity measure? These are very different statements.*
- the spread-spectrum follow-up changed trainability but not endpoint accuracy
- **current native-free scorers mostly cannot recognize the structures the VQE family can express**
  *Possibly the deepest item here. If nothing available can recognize a good structure, no amount of search, expressivity, or Hamiltonian design will help.*
- the mean is dominated by the difficult targets
- FAIL18 remains scientifically important
- a backbone-torsion signal remains open and is outside the previously bounded scorer class
- **a free subset-search objective produces non-prefix optima**, meaning there is a genuine combinatorial selection problem available that the deployed VQE has not yet solved
- the current architecture spends its information budget on the wrong question

The central question is now:

> **How do we make CVaR-VQE optimize a structurally meaningful problem whose solutions actually correspond to better protein structures?**

Do not assume the answer is one particular Hamiltonian.
Do not assume the answer is better calibration.
Do not assume the answer is a better ansatz.
Do not assume the answer is on the list in Section 8.

Find it.

Before committing to a direction, state explicitly which of these findings your hypothesis attacks and which it accepts as binding. A hypothesis that engages neither the recognition failure nor the common-mode error should carry a stated reason for bypassing both.

---

## 4. MAIN SCIENTIFIC PRIORITY

Priority order:

1. genuine reduction in built-chain RMSD
2. mechanistic understanding of why the reduction occurs
3. genuinely meaningful CVaR-VQE contribution
4. robustness and replication
5. computational efficiency

A beautiful quantum result that does not improve structure is scientifically useful, but secondary to a real accuracy improvement.

A sixth, unranked item: **a decisive measurement of what imposes the ceiling** — recognition, information, common-mode error, encoding, or readout — with evidence strong enough to redirect future work. This can be the most valuable output of the sprint if the accuracy result does not materialize.

---

## 5. SCIENTIFIC FREEDOM

You have broad freedom. Use it. This is deliberately **not** a prescribed experiment list.

You may:

- redesign the architecture, in whole or in part
- introduce entirely new Hamiltonians
- introduce new non-diagonal Hamiltonians
- redesign the basis-state encoding and what a qubit means
- redesign the candidate representation
- move from candidate IDs to structural variables
- abandon candidate identity as a concept entirely
- introduce interacting candidate variables
- introduce subset-selection Hamiltonians
- introduce graph-based structural couplings
- introduce constraint Hamiltonians
- introduce physically motivated off-diagonal operators
- use multiple Hamiltonians hierarchically
- use different Hamiltonians at different stages
- make the Hamiltonian adaptive
- redesign the CVaR objective, including its tail parameter and the quantity it acts on
- alter how entropy is used, or remove it
- alter temperature, or make it explicit and learned
- redesign the ansatz
- use ADAPT-VQE or related adaptive constructions
- redesign the readout
- preserve multimodality rather than collapsing it
- cluster quantum output
- interpret amplitudes as structural hypotheses
- use phase, not only amplitude, to carry information
- couple generation and selection rather than running them in sequence
- generate candidate structures from the variational state
- run the quantum stage iteratively rather than once
- change the information budget
- increase the number of qubits if scientifically justified
- reduce the number of qubits if that creates a better formulation
- investigate a larger-than-7-bit selection space
- introduce learned native-free statistical potentials
- introduce residual priors
- introduce uncertainty-aware objectives
- use torsion geometry
- use structural consistency
- use sequence-structure information
- use a free-energy formulation
- use a hybrid classical-quantum formulation
- use a completely new formulation you discover in the literature
- use a formulation nobody has published, if you can justify it

The only hard conceptual constraint:

> **CVaR-VQE remains the spine and main scientific object.**

Operationally that means: there is a parameterized quantum state; there is a Hamiltonian or observable set with defensible meaning; there is a CVaR or rigorously justified CVaR-family objective over a measured quantity; optimizing it is causally connected to final structure quality; and removing or randomizing the quantum stage measurably degrades the result.

Do not turn this into a generic classical protein-prediction project with a quantum afterthought.

Every inherited design choice that survives into your final system should survive because you examined it and kept it, not because it was already there. Be able to say why for each one.

---

## 6. LITERATURE RESEARCH

Before committing heavily to an approach, perform deep literature research. This sprint expects substantially more reading than previous ones, and expects you to return to the literature mid-sprint when results reshape the question.

Search primary and high-quality secondary sources on:

**Quantum optimization and variational methods**
CVaR variational quantum algorithms · CVaR-VQE · risk-sensitive quantum optimization · non-diagonal VQE Hamiltonians · variational quantum simulation · Hamiltonian variational ansätze · QAOA and VQE connections · quantum combinatorial optimization · subset-selection Hamiltonians · Ising/QUBO formulations · graph Hamiltonians · spectral graph theory · quantum walks · stoquastic and non-stoquastic Hamiltonians · free-energy variational principles · entropy regularization · quantum Boltzmann and Gibbs-like states · adaptive ansätze · ADAPT-VQE · qubit-ADAPT-VQE · Lie algebraic controllability · barren plateaus · gradient troughs · objective geometry · variational landscape geometry · quantum natural gradients · measurement strategies and estimator variance for tail objectives · shot-noise behavior of order statistics · classical simulability boundaries and what separates a real advantage claim from a simulable one

**Structure, statistics, and decision theory**
protein energy functions · structural priors · protein distance geometry · torsion potentials · Ramachandran distributions · statistical protein potentials · structural ensembles · multimodal protein representations · uncertainty-aware structure prediction · Bayesian decision theory · proper scoring rules · energy-based models · peptide structure prediction · fragment assembly · conformer selection · consensus structure methods · clustering of conformational ensembles · model quality assessment and native-free structure scoring · correlated / common-mode error in ensembles and methods that break it · residual and systematic-bias modeling · calibration and miscalibration diagnostics · submodularity and set-function optimization · determinantal point processes and diversity-aware selection · geometric deep learning representations of molecules

Read the equations, not the abstracts. Where a paper's claim rests on an assumption, check whether that assumption holds here — several of these literatures assume ensemble errors are independent, which finding 11 says is false for this pool.

Do not copy a method merely because it is established or fashionable.

For each imported idea, ask:

> **What information or mathematical structure does this contribute that the current pipeline does not already contain?**

That is the standard for adding complexity. "Formulated more elegantly" is not a reason. "Contains a source of signal orthogonal to what we have" is exactly the reason.

Keep a written record of what you read and what you rejected, with reasons. "This entire family assumes X, which our pool violates" is a real finding and belongs in the report.

---

## 7. THE COST/RMSD METER — BUILD THIS FIRST

Of everything carried forward from S29, this one is not optional, because it is infrastructure rather than a hypothesis.

Implement and formalize the cost-audit machinery. For any candidate cost function f(structure), produce a standard diagnostic measuring:

- ladder Spearman correlation between cost and RMSD
- local gradient cosine at production
- native percentile under the cost
- preference for ORACLE structure versus production
- preference versus matched random signed structures
- pool-member controls
- fold-clustered confidence intervals
- per-target distributions
- FAIL18 versus 108 split

Current shipped-cost baselines:

- ladder correlation ≈ **−0.40**
- production-point cosine ≈ **−0.03**
- native percentile ≈ **36.9**
- ORACLE preference signal ≈ **0.07**

Every new cost passes through this meter before receiving substantial VQE compute. Do not let an objective that moves in the wrong structural direction consume a day of compute.

Treat the meter as a permanent research instrument. Extend it when you find a diagnostic it is missing. If you design a cost whose merits the meter cannot see, that is itself worth reporting — it means the meter has a blind spot, and blind spots in your instruments are how sprints get lost.

---

## 8. LEADS REGISTER — CARRIED FORWARD FROM S29

**Read this section as evidence, not instruction.** These are the directions S29 left open, with what it knew about each. Pursue what your judgment supports. Reject what it does not, and say why in the report. Inventing a thirteenth lead that beats all of these is an entirely acceptable outcome.

### L1 — Calibrate the posterior
The distogram is roughly 2× over-confident. Possible approaches: temperature calibration, per-bin calibration, fold-safe recalibration, variance-aware weighting, alternative proper scoring constructions, other justified posterior transformations.

Requirements if pursued: leave-fold-out; native-free at deployment; evaluated with the cost/RMSD meter; evaluated on built-chain endpoints; strictly separated from native-dependent diagnostics.

Do not assume calibration works. S29's registered prior was skeptical because common-mode error may dominate. Test that hypothesis directly — and note that if common-mode error does dominate, that is a finding about the pool, not just about calibration.

### L2 — Remove the contraction blind spot
Production average is ~22% contracted. Investigate whether the cost is rewarding contraction rather than correctness.

Principled options: scale-adjusted structure evaluation, predicted-Rg normalization, geometry-aware cost normalization, de-contracted coordinates, torsion-space scoring, bond-length-aware distance-map construction, alternative structural coordinates.

Do not globally rescale the output and call it a result. The question is:

> **Can the objective distinguish a physically valid structure from an artificially contracted one without destroying legitimate structural variation?**

Meter first, then VQE.

### L3 — Residual prior
S29 rated this among the highest-upside accuracy directions. Pool error is strongly common-mode (~68%). Investigate whether a model can predict the systematic residual of the current distance prior using information sufficiently independent of the original prediction.

Possible sources: sequence-derived features, ESM embeddings, structural-context features, torsion information, candidate-pool statistics, geometry, uncertainty, complementary statistical potentials, anything else you find.

Critical: do **not** train a second model on the same feature representation and claim independent residual correction. Measure error correlation explicitly and require demonstrably different information. Use the pinned fold structure, train leave-fold-out, evaluate held-out cost geometry, then evaluate inside the CVaR-VQE.

Measure: cost/RMSD correlation, local gradient alignment, candidate recognition, built-chain RMSD, FAIL18 behavior, common-mode error reduction.

### L4 — Backbone-torsion channel
The major newly open direction. S29 found a torsion signal that survives removal of size/compactness effects, placing it outside the previously bounded scorer class.

Determine: what information it contains; whether it predicts native-like local conformations; whether it is complementary to the distogram; whether it improves ranking; whether it improves the meter; whether it gives a useful local gradient; whether it survives permutation and matched-random controls; whether it helps FAIL18; whether it is auxiliary-only or can support a full Hamiltonian.

Explore formulations, not one arbitrary scalar: local torsion compatibility, torsion distributions, Ramachandran energy, sequence-conditioned torsion consistency, torsion transitions, local pairwise couplings, torsion smoothness, secondary-structure-compatible sectors, non-additive interactions.

Do not assume it is useful because an audit left it open. Prove it or kill it.

### L5 — Free-energy / subset-selection CVaR-VQE
S29 flagged this as the highest-priority quantum reformulation, because a genuine subset objective has optima that are **not prefixes of the energy ordering**. That means a nontrivial combinatorial problem exists which is not equivalent to ranking candidate energies.

The question is whether CVaR-VQE can actually solve that problem.

Explore: subset selection, pairwise compatibility, structural coherence, subset energy, diversity/coverage, mutual compatibility, representative medoids, cluster-aware selection, basin selection, constrained ensembles, set functions, submodular and supermodular structure, graph cuts, QUBO/Ising reductions, nonlinear set objectives, variational free-energy formulations.

The quantum state should represent something more meaningful than "candidate 17 versus candidate 53." It might instead represent inclusion/exclusion variables, structural modes, compatible subsets, latent basin assignments, pairwise relations, fragment choices, torsion sectors, or structural constraints.

Ideal properties of a good formulation: the optimum is not trivially the classical energy prefix; the Hamiltonian is genuinely interacting; the state can represent correlated choices; the problem has real structural semantics; the classical equivalent is clearly defined for control; CVaR has a nontrivial role; the output maps to a valid structure.

Do not force a quantum advantage. If the classical problem is better, say so clearly — but determine whether VQE changes the optimization frontier.

### L6 — Sparse weighted readout
S29 found that oracle weights over only a few members can produce very low RMSD. Investigate whether a native-free sparse readout can identify a small number of structurally complementary members.

Options: sparse support selection, structural clustering, Bayesian support inference, topological clustering, subset VQE, medoid-support selection, low-rank reconstruction, sparse variational weights, entropy-constrained sparse distributions.

The system currently spends its budget on *how many to average* rather than *which structures*. S29 estimated the relevant selection task may need substantially more than 7 bits. Larger registers are permitted where justified — never merely to raise the qubit count.

### L7 — FAIL18
The mean is dominated by the failure tail. Treat this as a structural regime problem, not outlier removal.

Investigate: what makes FAIL18 different, fibril/lasso structure, candidate-pool structure, torsion behavior, score uncertainty, retrieval geometry, distance-map disagreement, structural diversity, sequence composition, alternative priors, torsion signatures, pool statistics.

A native-free regime detector is allowed but must be genuinely different from the closed router family. Do not build increasingly elaborate routers from the same posterior-derived features. The target is:

> **discover a scientifically meaningful latent regime distinction** — and then determine whether CVaR-VQE should behave differently in that regime.

### L8 — Information-budget study
S29 showed that changing how many candidates are averaged pays little, choosing the correct member pays a lot, and current ranking provides only about **1.44 of 7** useful bits.

Treat the quantum register as an information budget and ask:

> **How should the available bits be allocated to maximize structural information?**

Directions: candidate identity bits, subset bits, mode bits, torsion bits, structural basin bits, hierarchical encoding, factorized encoding, mixed encodings.

Do not assume binary candidate indexing is optimal. The encoding may be the hidden bottleneck, and it is the least-examined component in the system.

### L9 — Non-diagonal Hamiltonians, only if actually new
Do not spend a sprint on superficial graph variants. The first graph failed because its spectrum was essentially a typicality projector.

A new non-diagonal Hamiltonian needs a scientifically motivated interaction structure: structural compatibility, conformer transitions, torsional adjacency, pairwise consistency, graph Laplacians, constraint couplings, overlap operators, physically derived transition amplitudes, non-stoquastic interactions, Hamiltonian variational structures.

The question is not "is H off-diagonal?" It is:

> **"Does the off-diagonal structure encode information that materially changes the optimization problem?"**

Analyze: spectrum, eigenvectors, spectral participation, locality, commutator structure, DLA, representability, gradient variance, optimization behavior, classical counterpart, endpoint RMSD. Check for rank collapse before spending compute — that check is cheap and would have saved the previous attempt.

### L10 — ADAPT-VQE
Determine whether adaptive ansätze become meaningful under a new formulation. S26 found product-circuit behavior at α = 1 in the old diagonal problem; the question is whether a genuinely correlated structural Hamiltonian creates a meaningful adaptive growth process.

Investigate: operator-pool design, qubit-ADAPT, problem-inspired pools, graph-local pools, torsion-inspired pools, pairwise-compatibility pools, commutator-generated directions, gradient troughs, DLA growth, circuit depth, gradient variance.

Do not read large or small gradients alone as proof of trainability or barren plateaus. Distinguish product-state triviality, optimization difficulty, gradient troughs, genuine concentration, finite-size scaling, and DLA limitations.

### L11 — Off-pool recognition
The scorer library failed to recognize very good structures. Take this seriously — it may sit underneath several of the other leads.

Determine whether nativeness is recognizable from the geometry of an individual structure **at all**. Build careful leave-fold-out experiments using ORACLE-generated structures, random matched structures, production structures, projected structures, structurally perturbed structures, torsion-space perturbations, and other meaningful ladders.

Candidate predictors: torsion features, local geometry, secondary-structure patterns, distance-map consistency, contact structure, structural entropy, local smoothness, backbone geometry, physically motivated terms, statistical structural descriptors.

A successful cost here could become the Hamiltonian objective. This is a fundamentally different question from ranking the original pool, and a negative answer would be one of the most important results the project could produce.

### L12 — Combinations and anything not listed
The leads above are not independent. L3 and L11 may be the same problem viewed from two sides. L5 and L8 both concern what the register should encode. L2 may be a special case of L11. Look for the formulation that dissolves several at once rather than the one that patches each in turn — and be open to the possibility that the right move is none of these.

---

## 9. AGENT ORCHESTRATION

Maintain **at least 4 active agents whenever meaningful work is available.** Never exceed **8 concurrent agents.** You decide whether to run 4, 5, 6, 7, or 8.

You decide how many agents belong on each direction, which do literature research, which implement, which derive mathematics, which run experiments, which reproduce, which inspect artifacts, which act as adversaries, and when to terminate and reassign.

Standing requirements:

- At least one active agent must **continuously attack the strongest positive or leading hypothesis.** This is a permanent role, not a phase.
- At least one active agent must pursue a **substantially different mechanism**, not incremental variations of the current favorite. If the leading idea is a new Hamiltonian, this agent questions whether the state space is right at all.
- At least one agent should be **reading literature for most of the sprint**, not only at the start.
- Do not let all agents converge prematurely on one idea. Parallel agents on the same hypothesis is wasted parallelism.
- Reassign finished agents immediately. No artificial idle time.
- Rotate the adversary role so no agent is only ever attacking its own prior conclusions.

Maintain a running written state: current leading hypothesis, strongest live objection to it, its preregistered falsifier, what each agent is doing, what has been closed and why. Update it as results land. Do not hold it only in your head — with 8 agents, an unwritten research state is how contradictory conclusions coexist unnoticed.

---

## 10. RESOURCE USE

The machine is intentionally being used aggressively.

Target roughly **94–95% overall CPU/RAM utilization** when meaningful work is available, while respecting the existing governor, memory safety, AMBER serialization, jobrun limits, checkpointing, and OS stability.

Do **not** exceed safe limits to produce a utilization number.

- Use the governor intelligently.
- Prefer many medium concurrent jobs over one enormous fragile job when that yields more scientific throughput.
- Never let memory pressure destroy long runs or lose checkpoint state.
- Every expensive experiment must be resumable.
- Do not throw away work because of a resource kill.

Budget deliberately: decide roughly how much of the sprint goes to reading and theory, how much to cheap probes, and how much to full-instrument runs — and revise that split consciously rather than letting it drift. High utilization on the wrong hypothesis is not throughput.

---

## 11. STATISTICAL REQUIREMENTS

For every major experiment:

- pre-register the hypothesis
- pre-register the falsifier
- use paired target-level comparisons
- calculate MDE separately for each comparison
- report SE
- report fold-clustered CI
- retain i.i.d. CI as secondary
- run concentration tests
- run permutation controls where appropriate
- use matched-random controls
- use order-statistic corrections
- use leakage poisoning
- replicate promising findings
- separate ORACLE from deployable rows

Hard rules:

- Never use the native to tune a deployable parameter.
- Never let a native-dependent ceiling become a production result.
- Do not change benchmark60. Do not regenerate the pinned folds or clusters.
- Anything below **0.7× MDE is not a result.** Between 0.7× and 1× MDE is not a demonstrated improvement.
- Do not call a borderline numerical movement a breakthrough.
- Every strong positive must survive the Adversary.

Also track your own multiplicity. Across a long sprint with 8 agents running many comparisons, a spurious 1× MDE positive is expected rather than surprising. Count the comparisons run and hold positives to a standard that accounts for the count.

---

## 12. ADVERSARY AGENT

One active agent functions as a hostile reviewer at all times.

Its job: attack every positive; look for order-statistic effects; look for concentration; detect hidden native leakage; compare against matched classical controls; inspect normalization; inspect effective-temperature changes; inspect projection artifacts; inspect tie conventions; inspect basis changes; inspect implementation drift; inspect whether the effect is simply a classical selection effect; inspect whether the quantum effect is merely entropy.

Add to the standing checklist:

- Does it survive both seeds?
- Does it survive fold clustering?
- Does it survive a simpler control?
- Does it survive at matched budget?
- Does it remain on the built chain?
- Is the mechanism independently measured, or only inferred from the endpoint?
- Would I believe this result if a rival group reported it with exactly this evidence?

The main team may not approve its own positive. A positive that has not survived a genuine attempt to destroy it is provisional and must be labeled as such in every document it appears in.

---

## 13. ARCHITECTURE POLICY

- Production remains recoverable.
- Historical experiments remain untouched.
- New work lives in a clearly separated experimental area.
- You may redesign the experimental architecture aggressively.
- Do not delete inconvenient results.
- Do not overwrite historical ledgers.
- Do not change benchmark60.
- Do not regenerate the pinned folds or clusters.
- Do not mutate the historical production commit.
- Any new architecture must be reproducible from its own artifacts, with seeds and configuration recorded.
- If a new formulation proves superior, document exactly how it differs.

---

## 14. THE QUANTUM SPINE — NO PRESCRIBED ARCHITECTURE

Previous sprints supplied a conceptual pipeline template. **That template is withdrawn.** Not because it was wrong, but because it has become an anchor, and the shape of the system is now part of the research question.

There is no diagram in this prompt to work toward. Design the architecture you can defend.

The only criterion for whether a stage belongs:

> **Does this stage help CVaR-VQE solve a structurally meaningful optimization problem and reduce final built-chain RMSD?**

Do not preserve a stage because it existed before. Do not preserve an ordering because it was the ordering. Do not assume the system is a feed-forward pipeline at all — it may be a loop, a multi-stage process, a coupled generator-selector, or something with no clean stage boundaries.

If your final system does resemble the old template, that is a legitimate outcome, but it should be a conclusion you reached, not a starting assumption you never tested.

---

## 15. DEEP MATHEMATICS

Where useful, derive: exact objective gradients · envelope derivatives · CVaR subgradients · entropy/free-energy identities · Hamiltonian spectra · spectral decompositions · graph Laplacians · DLA closures · parameter-count and expressivity bounds · information-theoretic bounds · subset objective formulations · representability bounds · structural distance geometry · bias/variance decompositions · error-vector decompositions · decision-theoretic relationships between ranking and structural reconstruction · estimator bias under finite shots and finite pools.

For every proposed objective, analyze whether it is: well-defined · differentiable or subdifferentiable, and where not · stable under resampling and shot noise · properly normalized and scale-invariant where it should be · compatible with CVaR's order-statistic structure · identifiable · meaningful under finite candidate pools · likely to preserve or destroy multimodality · vulnerable to order-statistic artifacts · vulnerable to entropy-induced trivial solutions · vulnerable to contraction artifacts · vulnerable to common-mode error · vulnerable to degenerate spectra or rank collapse · monotone with respect to anything you actually care about.

Use mathematics to answer the project's actual bottlenecks. Do not add equations for appearance — if you cannot name the failure mode a piece of math addresses, cut it. Conversely, do not avoid hard mathematics when the problem demands it. If the right answer requires the spectrum of an operator family, the reachable set of an ansatz, or the bias of an estimator, work it out properly.

---

## 16. PHYSICS DEPTH

Investigate whether the Hamiltonian has a real physical interpretation: steric compatibility · torsion energetics · hydrogen bonding · electrostatics · solvation · local backbone physics · conformational transitions · structural compatibility · free-energy landscapes · coarse-grained interactions · conformational basin connectivity.

Do not assume a physical force field is automatically a good inference objective. Keep these three explicitly distinct at all times:

- **physical validity** — is this structure possible?
- **energetic plausibility** — is this structure favorable?
- **native discrimination** — is this structure the right one?

Conflating them is a standing risk in this project, and a cost function that optimizes the first two while ignoring the third would reproduce exactly the failure mode S29 measured.

---

## 17. COMPUTATIONAL DEPTH

Run: probes · full 126-target evaluations · replications · ablations · classical controls · exact small-n quantum checks · MPS checks · larger-n studies only where justified.

- Do not burn compute on a hypothesis that already fails its low-cost falsifier.
- Design the cheap decisive measurement before the expensive ambiguous one. A 20-minute experiment that kills a direction beats a day of compute that confirms one.
- Use checkpointing everywhere. Preserve intermediate outputs. Never recompute complete experiments unnecessarily.
- Build **measurements of mechanism**, not only measurements of outcome. If something helps, show the intermediate quantity that moved.

---

## 18. WHAT COUNTS AS A BREAKTHROUGH

A breakthrough is **not**: lower energy · lower point-cloud RMSD only · a more expressive circuit · a larger gradient · lower gradient variance · more qubits · a more complicated Hamiltonian · more candidates · more parameters · a better correlation that does not survive to the built chain.

A breakthrough **is**:

> **a reproducible reduction in built-chain RMSD with a clear causal mechanism, surviving the project's controls.**

A quantum breakthrough additionally requires:

> **the CVaR-VQE to contribute information or optimization behavior that a matched classical control does not reproduce.**

For any serious quantum formulation, separately establish: what the basis states mean · what the Hamiltonian means · why it should correlate with useful structural information · what CVaR is optimizing and over what distribution · what the ansatz can and provably cannot represent · whether the optimizer can reach the relevant states · whether the quantum output holds information unavailable classically · whether a classical control reproduces it · whether it changes built-chain RMSD.

Controls for any quantum positive: classical equivalent · diagonalized equivalent · permuted · random · matched-budget · untrained circuit · simpler ansatz · order-statistic controls · a product-state/separable restriction to test whether entanglement is doing anything at all.

Do not call something quantum merely because a quantum circuit produced it.

---

## 19. FINAL REPORT

Produce a report in the same scientific style as S28/S29, answering:

1. What was the strongest new architecture?
2. What was the final mean built-chain RMSD?
3. Did it beat 3.21 Å?
4. Did it beat 3.0 Å?
5. Did it beat 2.5 Å?
6. What was the exact paired effect versus production?
7. What was the fold-clustered CI?
8. What was the MDE?
9. What happened on FAIL18?
10. What happened on the other 108?
11. What information source produced the gain?
12. What did CVaR-VQE specifically contribute?
13. What Hamiltonian was ultimately used?
14. Why is it meaningfully different from the earlier diagonal rank-ladder?
15. What did the quantum state represent?
16. Could a classical control reproduce the result?
17. What happened to trainability?
18. What happened to gradient variance?
19. What happened to DLA / ansatz structure?
20. Which hypotheses were killed?
21. Which remain open?
22. What is the next highest-value scientific question?

Additionally include:

23. Which of the twelve S29 leads you pursued, which you rejected, and why.
24. What literature you relied on, and what you read and rejected.
25. Every control you ran and what each one ruled out.
26. A full list of hypotheses entertained and killed, with reasons. The dead ends are part of the result.
27. Where the cost/RMSD meter's baselines ended up for every serious cost function tested.
28. A final architecture diagram of whatever you actually built.

No hand-waving. No unsupported claims. No hidden tuning. No deletion of inconvenient results.

---

## 20. FINAL MINDSET

Do not think:

> "What experiment should I run next?"

Think:

> **"What is the deepest reason this system cannot currently turn optimization into structural accuracy, and what mathematical reformulation could remove that reason?"**

Do not optimize harder on a bad objective.
Do not keep adding Hamiltonians that contain no new information.
Do not keep increasing qubits without changing the problem.
Do not mistake an ORACLE ceiling for an achievable result.
Do not settle for a tiny null because it is easy.
Do not force a positive because the target is ambitious.
Do not let a list of twelve recommendations substitute for your own scientific judgment.

You have permission to take large conceptual swings.
You have permission to redesign the architecture.
You have permission to invent a new Hamiltonian.
You have permission to invent a new encoding.
You have permission to formulate a new CVaR-VQE objective.
You have permission to combine protein physics, statistical decision theory, graph theory, quantum optimization, information theory, and structural biology.
You have permission to ignore this document's suggestions when the evidence points elsewhere.
You have permission to abandon an idea quickly when its falsifier fires.
You have permission to spend substantial compute when the hypothesis survives.
You have permission to spend a large fraction of the sprint thinking before building.

Use maximum research depth.
Use maximum mathematical depth.
Use maximum physics depth.
Use maximum experimental depth.
Use maximum adversarial testing.
Use maximum creativity.
Use maximum scientific ambition.

But above all:

> **Make the CVaR-VQE genuinely useful for structure selection, and make the final built-chain RMSD move.**

Go find the mechanism. Final report through atrifact skill after EVERYTHING is done. Broad freedom final goal is RSMD. I believe in you! You got this!