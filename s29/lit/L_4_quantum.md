# L_4 -- THE QUANTUM LITERATURE, AND THE CONDITION FOR NON-CLASSICALITY (charter section 6 group 1; brief topic 4)

Lane L, 2026-09-20. The brief's question, which is the charter's question in section 11: under
what condition does the quantum object differ from its classical counterpart, and what
Hamiltonian structure would make the variational state NOT classically reproducible by a sort or
an eigensolver? Name the papers that show a variational state doing something an eigensolver on
the same H cannot.

The answer, stated first, in three parts.

1. For a DIAGONAL Hamiltonian, CVaR is a classical order statistic BY DEFINITION, not by an
   accident of our implementation. Barkoutsos et al. define CVaR on the SORTED samples, and
   state that their target is "classical optimization problems, which yield diagonal
   Hamiltonians". The project's set-equality theorem (S25: the CVaR tail's support is always a
   subset of an initial prefix of the energy order, 2,592 adversarial cells, 0 violations) is not
   a discovery about our pipeline -- it is eq (12) of the source paper. The only quantum content
   in CVaR-VQE on a diagonal H is the SAMPLING DISTRIBUTION the trial state induces over
   bitstrings. Anything else has to come from a non-diagonal H.
2. There is a second, sharper obstruction that the charter should know about, and it is recent:
   Cerezo et al. (Nat Commun 16:7907, 2025; arXiv:2312.09121) collect case-by-case evidence that
   the very structure which makes a variational model free of barren plateaus also tends to make
   its loss classically simulable, because "barren plateaus result from a curse of dimensionality,
   and ... current approaches for solving them end up encoding the problem into some small,
   classically simulable, subspaces". Trainability and classical hardness pull against each
   other. Any S29 quantum claim must be positioned against this, because "our circuit trains
   well" is now evidence FOR classical simulability, not against it.
3. Therefore the honest target for this sprint is not "beat an eigensolver". It is a formulation
   where the object being prepared is NOT an eigenstate of the Hamiltonian at all -- a thermal /
   Gibbs-like state, or a distribution used as a generator -- built on a Hamiltonian whose terms
   do not commute. That is the only structure in which "sort the diagonal" and "diagonalise H"
   are both the wrong classical counterpart. The project's own record already tells us the cost:
   S28's first non-diagonal Hamiltonian had a near-rank-one similarity graph, so its
   off-diagonal carried almost no information (S28-L8b/L11).

---

## 1. CVaR: what it is, and two properties the project should use

### 1.1 Barkoutsos PKl, Nannicini G, Robert A, Tavernelli I, Woerner S. "Improving Variational Quantum Optimization using CVaR." Quantum 4, 256 (2020); arXiv:1907.04769

Construction, verbatim. "CVaR_alpha(X) = E[X | X <= F_X^{-1}(alpha)]" (their eq 11), and for K
samples sorted in nondecreasing order H_{k+1} >= H_k,

    CVaR_alpha = (1 / ceil(alpha K)) * sum_{k=0}^{ceil(alpha K)} H_k                     (12)

"the limit alpha -> 0 corresponds to the minimum, and alpha = 1 corresponds to the expected
value of X. In this sense, CVaR is a generalization of both the sample mean and the best observed
sample." The motivation, verbatim from the abstract: "In the case of classical optimization
problems, which yield diagonal Hamiltonians, we argue that aggregating the samples in a different
way than the expected value is more natural."

Check against this instrument. Eq (12) IS the set-equality theorem. The tail is a prefix of the
sorted samples; nothing else is available to it. The project's S25 result therefore has an
external, definitional source, and should be cited as such rather than presented as a surprise:
on a diagonal H, CVaR-VQE's aggregation is classical by construction, and the quantum stage's
only possible contribution is which bitstrings get sampled and with what probability.

### 1.2 The two properties of CVaR the project should be using, and one of them is a warning

(a) DEGENERACY OF THE OPTIMUM. The paper states: "for any problem (1) and parameters theta* such
that |psi(theta*)> has overlap rho > 0 with the ground state, theta* is a global minimum of
CVaR_alpha(X(theta*)) for alpha <= rho". So the global-minimiser SET of CVaR_alpha is
{theta : overlap with the best state >= alpha} -- an enormous, flat set. CVaR at small alpha
cannot distinguish any two states that both place at least alpha of their mass on the best
candidate.

This is a mechanism for one of the sprint's central puzzles, and I believe it has not been
stated in the record in this form. The project observes that the optimiser reaches the
objective's optimum on 126/126 while the emitted structure does not move (S28-L18b, S28-L26b).
With a readout that AVERAGES over the tail, the objective is indifferent to exactly the degree of
freedom the readout consumes: CVaR fixes the tail's mean energy, but the structure depends on
WHICH other candidates populate the tail, and over the flat optimum set that is unconstrained.
"Optimising the objective does not determine the output" is then not a pathology of our
Hamiltonian; it is a property of CVaR plus an averaging readout. FLAGGED for lane T to check
and for lane D's meter: an objective whose argmin is a large flat set cannot, by itself, be
blamed for the emitted structure -- and conversely, any accuracy change attributed to CVaR
optimisation must be shown not to be a tie-break within that flat set.

(b) ESTIMATOR VARIANCE. "The variance of the empirical CVaR_alpha estimator using K samples is
O(1/(K alpha^2)) ... implying that the resulting standard error increases as 1/alpha. Thus, for
a fixed number of samples K, to achieve the same accuracy as for the expected value we need to
increase the number of samples to K/alpha." Not binding for this project (the simulator is exact
and the statevector is available, S25), but it is the correct price to quote if any arm ever
moves to sampling, and it is the reason alpha cannot be taken arbitrarily small in practice.

(c) Proposition 5.1: "A local minimum of (1) does not necessarily correspond to a local minimum
of (13), and vice versa" -- the expectation-value landscape and the CVaR landscape are genuinely
different objects, with a two-qubit example where the expectation is constant in theta (no
optimisation possible) while CVaR_0.5 = sin^2(theta/2) is informative. This is the paper's real
argument for CVaR and it is about LANDSCAPE, not about solution quality.

Information test for the CVaR family: it contains a landscape repair for flat expectation-value
objectives. This project's objective is not flat in that way (the optimiser reduces it on
126/126), so the repair is not the binding issue here.
KEPT as the definitional source (it closes "is our set-equality theorem a bug?" -- no, it is the
definition) and for properties (a) and (b).

---

## 2. Trainability, and the result that reframes the whole charter

### 2.1 Cerezo M, Sone A, Volkoff T, Cincio L, Coles PJ. "Cost function dependent barren plateaus in shallow parametrized quantum circuits." Nat Commun 12:1791 (2021); arXiv:2001.00550

Verbatim from the abstract: "Our first result states that defining C in terms of global
observables leads to exponentially vanishing gradients (i.e., barren plateaus) even when V(theta)
is shallow ... our second result states that defining C with local observables leads to at worst
a polynomially vanishing gradient, so long as the depth of V(theta) is O(log n). Our results
establish a connection between locality and trainability." Assumption: an alternating layered
ansatz whose blocks form local 2-designs.

Check. The assumption does not hold here and the project already knows it: S13's ansatz is
"nowhere near a 2-design" (`docs/STATE_BRIEF` 5.7 item 4), and the project measured no barren
plateau at any width (S25). Contract rule 9 forbids calling any slope a barren plateau or its
absence, and this paper is the reason the distinction matters: locality of the OBSERVABLE, not
depth alone, is the lever. The project's observable is the pool-energy ladder, a global object
over the candidate register. NOTED, not actionable: the project is not in the 2-design regime
where the theorems bite.

### 2.2 Larocca M, Czarnik P, Sharma K, Muraleedharan G, Coles PJ, Cerezo M. "Diagnosing barren plateaus with tools from quantum optimal control." Quantum 6, 824 (2022); arXiv:2105.14377

Construction. Build the DYNAMICAL LIE ALGEBRA g generated by the ansatz's Hermitian generators
under commutation; the reachable set of the ansatz is the group exp(g), and controllability
(g = su(2^n)) versus uncontrollability is the diagnostic. The paper links gradient scaling to
the scaling of dim(g), and proves no-go results for obtaining ground states with variational
ansatze in the controllable case.

Check, and this is a concrete gap in the project's record: `docs/STATE_BRIEF` 5.7 item 4 lists
the ansatz's DLA as NOT MEASURED. This is the one quantum diagnostic the charter names
(section 6) that the project has never computed, it is purely classical linear algebra on the
generators, it costs minutes, and it determines both the reachable set (charter question 5:
"what the ansatz can represent, and what it provably cannot") and the expected gradient scaling.
RECOMMENDED to lane T or X as a cheap, decisive, non-endpoint measurement. I do not run it.

### 2.3 Ragone M, Bakalov BN, Sauvage F, Kemper AF, Ortiz Marrero C, Larocca M, Cerezo M. "A Lie algebraic theory of barren plateaus for deep parameterized quantum circuits." Nat Commun 15 (2024); arXiv:2309.09342

Verbatim from the abstract: "we present a general Lie algebraic theory that provides an exact
expression for the variance of the loss function of sufficiently deep parametrized quantum
circuits ... This theoretical leap resolves a standing conjecture about a connection between loss
concentration and the dimension of the Lie algebra of the circuit's generators."

The variance is expressed through g-purities of the input state and the observable divided by
dim(g); the direction is INVERSE -- a larger DLA gives a smaller variance, i.e. worse
trainability. CAVEAT ON MY SOURCING: I read the abstract and the surrounding text, not the
displayed equation, so I quote the direction and the structure (P_g(rho) P_g(O) / dim g) and not
a verbatim formula; anyone using the exact expression should pull the equation from the paper.
The direction matters because it sets up section 2.4.

### 2.4 Cerezo M, Larocca M, Garcia-Martin D, Diaz NL, Braccia P, Fontana E, Rudolph MS, Bermejo P, Ijaz A, Thanasilp S, Anschuetz ER, Holmes Z. "Does provable absence of barren plateaus imply classical simulability?" Nat Commun 16:7907 (2025); arXiv:2312.09121

Verbatim from the abstract: "Can the structure that allows one to avoid barren plateaus also be
leveraged to efficiently simulate the loss classically? We collect evidence-on a case-by-case
basis-that many commonly used models whose loss landscapes avoid barren plateaus can also admit
classical simulation, provided that one can collect some classical data from quantum devices
during an initial data acquisition phase. This follows from the observation that barren plateaus
result from a curse of dimensionality, and that current approaches for solving them end up
encoding the problem into some small, classically simulable, subspaces."

The notion of simulability is conditional: efficient classical computation of the loss AFTER an
initial data-acquisition phase on a quantum device, not universal simulation for all inputs. The
authors list their own caveats: average-case arguments, smart initialisations, models outside
their assumptions, the potential for provably superpolynomial advantages, and heuristic behaviour
at larger scale.

Check against this instrument, and this is the most important framing item in topic 4. The
project's trainability result -- the one genuine quantum positive in the record (S13/S25: the
optimiser trains, beats best-of-200 from the untrained circuit, closes 78-89% of the free-energy
gap, no barren plateau at any measured width) -- sits exactly in the regime this paper warns
about. Small, structured, trainable, low-dimensional: the paper's argument says such a model is
a candidate for classical simulation of its loss. This does NOT retract the project's result
(which is a trainability claim, correctly scoped, and the project already states that an exact
simulator is never a cause, contract rule 9). What it does is tell the sprint how to position
any quantum claim: showing "we train well" is no longer evidence of quantum content, and a
"classical equivalent" control (charter section 11's first control) is not a formality but the
central test. The charter's ten controls are, in effect, this paper's programme.
KEPT as the framing result for every quantum claim in S29.

---

## 3. What would make the variational state non-classical here

Putting the pieces together, for the charter's question. The classical counterpart of our current
formulation is a SORT (because H is diagonal and CVaR reads a prefix of it, section 1.1). An
eigensolver is the classical counterpart as soon as H is non-diagonal but the target is its
ground state. So a formulation is only non-classical in the relevant sense if BOTH of the
following hold:

(C1) the Hamiltonian's terms do not commute, so the eigenbasis is not the computational basis
     and "sort the diagonal" is not available; AND
(C2) the object being prepared is NOT an eigenvector of that Hamiltonian, so "diagonalise H" is
     not the classical counterpart either.

(C2) is the part the project has never used. It is satisfied by a THERMAL / GIBBS state
rho = e^{-beta H}/Z, by a state whose role is to be a SAMPLING distribution rather than a
minimiser, and by any objective that is a free energy (an energy MINUS an entropy) rather than an
energy. The project's own record is suggestive here from an unexpected direction: S25 measured
that the trained state sits 0.902 nats / 45% of its mass away from its own analytic Gibbs
optimum, and that RMSD tracks READOUT ENTROPY (rho -0.74) rather than alpha or T. The entropy is
where the output actually lives, and a Gibbs-type target is the formulation in which entropy is
the object rather than a side effect.

Two caveats, stated so nobody over-reads this:
- For a DIAGONAL H, the Gibbs state is a classical Boltzmann distribution over candidates, which
  is samplable classically whenever the partition function over 2^n candidates is tractable --
  and S21 exhaustively enumerated the 2^n latent on 75/126 targets. So (C2) alone does not buy
  non-classicality; (C1) and (C2) are both required.
- Even with both, section 2.4 says a trainable small-DLA model is a simulability candidate. The
  defensible claim available to this sprint is therefore NOT "classically impossible" but "a
  specific quantity that the matched classical control in our ten-control set does not
  reproduce", which is exactly what charter section 17 calls a strong scientific success.

### 3.1 Papers where a variational state does something an eigensolver on the same H cannot

The brief asks for these by name. The honest answer is that the class is narrow, and it is not
the optimisation literature:
- Variational Gibbs-state preparation and variational thermal states: the prepared object is
  e^{-beta H}/Z, which no eigensolver on H returns; the standard constructions minimise a free
  energy F = <H> - T S, with the entropy estimated by an auxiliary register or a truncated
  expansion. This is a genuine "not an eigenvector" target.
- Quantum Boltzmann machines (Amin, Andriyash, Rolfe, Kulchytskyy, Melko, Phys Rev X 8:021050,
  2018): the model distribution is the measurement distribution of a thermal state of a
  TRANSVERSE-FIELD Hamiltonian; the off-diagonal terms make the model's distribution not equal
  to any classical Boltzmann distribution over the same variables, which is the precise sense in
  which a quantum Boltzmann machine differs from a classical one. That is the cleanest existing
  statement of (C1)+(C2) together, and it is the reference the charter's "quantum Boltzmann
  machines" bullet is pointing at.
- QAOA's mixer, as the meaning of an off-diagonal term: the mixer B = sum_i X_i is what makes the
  dynamics non-trivial; with the cost Hamiltonian alone (diagonal) the state never leaves its
  initial computational-basis distribution. This is the structural statement the charter wants
  about "what an off-diagonal term means": it is the only thing that MOVES AMPLITUDE BETWEEN
  CANDIDATES. S28's non-diagonal Hamiltonian failed not because off-diagonals are useless but
  because its particular similarity graph was near rank one (S28-L8b/L11) -- a degenerate mixer
  moves amplitude along one direction only.
- ADAPT-VQE (Grimsley HR, Economou SE, Barnes E, Mayhall NJ, Nat Commun 10:3007, 2019;
  arXiv:1812.11173): "instead of fixing an ansatz upfront, this algorithm grows it systematically
  one operator at a time in a way dictated by the molecule being simulated. This generates an
  ansatz with a small number of parameters, leading to shallow-depth circuits." The selection
  rule is the energy gradient with respect to adding operator A_i at zero angle, i.e. the
  commutator expectation |<psi|[H, A_i]|psi>|, maximised over the pool. It is an ansatz-
  construction method, not a source of non-classicality: it changes the reachable set, and its
  relevance here is that the operator pool is where problem structure enters. If a lane builds a
  non-diagonal Hamiltonian with genuine structure, ADAPT's selection rule is the principled way
  to grow a circuit matched to it, and the DLA of the resulting generator set (2.2) is what
  should be reported.

---

## 4. Verdict table (topic 4)

| # | paper | what it adds | verdict |
|---|---|---|---|
| 1.1 | Barkoutsos et al., Quantum 4:256 (2020), eqs (11), (12) | the CVaR definition ON SORTED SAMPLES, and the explicit scoping to diagonal Hamiltonians | KEPT -- the definitional source of the project's set-equality theorem; cite it instead of presenting S25's theorem as a surprise |
| 1.2 | the same paper's Prop 5.1 and the overlap result | CVaR's global optimum set is {theta : overlap >= alpha}, i.e. flat and large; estimator SE grows as 1/alpha | KEPT and FLAGGED -- a candidate mechanism for "the objective is optimised but the structure does not move"; lane T to check |
| 2.1 | Cerezo et al., Nat Commun 12:1791 (2021) | locality of the OBSERVABLE, not just depth, controls gradient decay (under a 2-design assumption we do not satisfy) | NOTED, not actionable -- our ansatz is nowhere near a 2-design and rule 9 applies |
| 2.2 | Larocca et al., Quantum 6:824 (2022) | the dynamical Lie algebra as the diagnostic for reachability and gradient scaling | KEPT as a RECOMMENDED cheap measurement: the project's DLA is listed as never measured, it is classical linear algebra, and it answers charter question 5 |
| 2.3 | Ragone et al., Nat Commun 15 (2024) | an exact variance expression through the DLA; variance falls as dim(g) grows | KEPT with a sourcing caveat (abstract read, equation not quoted verbatim) |
| 2.4 | Cerezo et al., Nat Commun 16:7907 (2025), arXiv:2312.09121 | trainable, barren-plateau-free models tend to be classically simulable | KEPT -- the framing result: "it trains" is now evidence FOR simulability; the charter's ten controls are this paper's programme |
| 3.1 | Amin et al., Phys Rev X 8:021050 (2018), quantum Boltzmann machines | a thermal state of a non-commuting H whose measurement distribution is not a classical Boltzmann distribution | KEPT as the cleanest existing instance of (C1)+(C2) |
| 3.1 | QAOA's mixer | an off-diagonal term is the only thing that moves amplitude between candidates | KEPT as the structural reading of S28's rank-one failure |
| 3.1 | Grimsley et al., ADAPT-VQE, Nat Commun 10:3007 (2019) | grow the ansatz by the commutator gradient over an operator pool | KEPT conditionally -- the right way to build a circuit for a structured non-diagonal H; not itself a source of non-classicality |

KEPT 7, NOTED 1, none importable as an accuracy operator. Topic 4 yields no RMSD lever; it yields
a definition, a warning, a cheap measurement, and a condition.

## 5. What topic 4 tells the sprint

1. Stop treating the set-equality theorem as a finding about our pipeline. It is eq (12) of
   Barkoutsos et al.: on a diagonal H, CVaR IS a classical prefix, by definition. Cite it.
2. The condition for non-classicality is (C1) non-commuting terms AND (C2) a target that is not
   an eigenvector. Gibbs/thermal targets and generator-style sampling satisfy (C2); S28's
   non-diagonal attempt satisfied (C1) with a degenerate (rank-one) off-diagonal. Neither has
   been tried with both.
3. CVaR's global optimum is a large flat set ({overlap >= alpha}). Any claim that CVaR
   optimisation caused an accuracy change must rule out a tie-break within that set. This is a
   new, cheap falsifier the sprint can apply to its own quantum arms.
4. The DLA of the deployed ansatz has never been measured, is classical linear algebra, and
   answers the charter's "what can the ansatz provably not represent". Cheapest open quantum
   question in the project.
5. Position every quantum claim against arXiv:2312.09121: trainability is no longer evidence of
   quantum content. The reachable major result is "a specific quantity the matched classical
   control does not reproduce", not "classically impossible".
