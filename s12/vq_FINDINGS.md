# s12 QUANTUM-ROLE — findings

**Question.** Is there a scientifically real combinatorial optimisation problem in this
pipeline where the value of a choice depends on the other choices — and how does CVaR-VQE
compare, honestly, against classical solvers on it?

**Headline.** Yes — and it is **anchored multi-segment fragment assembly**, whose objective
(the shipped distogram Bayes risk) decomposes **exactly, bit-for-bit** into a fully-connected
2-local Ising Hamiltonian over one-hot segment choices. That is a genuine textbook QUBO over
the real production objective, not a proxy. **The interaction is real but small (2.4 % of the
objective's variance) and the landscape is nearly unimodal, so a 70-microsecond greedy +
1-opt local search reaches the certified optimum on 92–98 % of instances.** At the project's
own shared unique-evaluation budget the CVaR-VQE loses to that local search on **0 of 126
targets** (d = +0.014 [+0.009, +0.019]) and to annealing on **0 of 126** in a second problem
family (d = +0.020 [+0.014, +0.026]); it also loses to uniform random sampling on 4 of 5
budgets, because its state entropy never falls below 9.1 of 12 bits. **CVaR does deliver the
one thing claimed for it — it prevents the collapse to the argmin (0.000 bits at α=1 vs 6.1
bits at α=0.05, converged under 12× the optimisation, matching an analytic bound) — and that
breadth is worthless: a matched-entropy Boltzmann distribution beats the VQE's readout on
20 of 20 α×T arms, and the VQE's most-probable states never cover the near-native band better
than the lowest-energy states.** End to end, every arm is worse than doing nothing (exact
optimum +0.169 Å [+0.103, +0.236] 43W/83L). **Recommendation: keep the component, re-site it
from the top-128 selector (whose energy is a per-candidate score, i.e. no interaction
structure at all) onto this Hamiltonian, and publish the classical comparison beside it.
A genuine quantum method on a genuine instance, benchmarked and losing, is a defensible
scientific role; a decorative selector is not.** One incidental non-quantum positive: inside
an *anchored* family the shipped objective has real ranking skill (corr +0.449, optimum at
the 24.8th RMSD percentile, shuffled-objective and rg nulls both passed) — the record's
"the objective does not rank" was measured on unanchored sets.

---

## 0. Instrument validation (run first, before any claim)

`python -m s12.instrument`:

| quantity | expected | observed |
|---|---|---|
| shipped distogram argmin | 3.4540 | 3.4540004952559396 |
| K=500 pool best (oracle) | 1.7108 | 1.7108244199364904 |
| top-75 best (oracle) | 2.3062 | 2.3061526409453816 |
| synthesis (avg → project, λ=0) | 3.2041 | 3.2040761603809194 |
| zero-recall targets | 18 | 18 (== FAIL18) |

Code: `s12/vq_lib.py` (instances + Hamiltonian), `s12/vq_classical.py` (solvers),
`s12/vq_quantum.py` (CVaR-VQE arm), `s12/vq_run.py` (drivers).

---

## 1. Choosing the problem

Two candidates were handed to me.

| | Problem A — multi-piece assembly | Problem B — subset selection for the terminal operator |
|---|---|---|
| variables | one library fragment per segment | which of the shipped candidates to average |
| interaction | cross-segment pair distances | risk of an average ≠ average of risks |
| exact form | **exactly 2-local** in the segment choices (proved below) | quadratic only after an `sd`-weighted surrogate; the true `F` carries a `1/k` and `1/k²` normalisation, so it is a *ratio* of quadratics, not a QUBO |
| instances | built here (`vq_lib.make_instance`) | already dumped by the aggregation agent (`s12/results/agg_subset_instances.json`), 126 × 16 binary vars, certified exhaustive optimum |
| oracle headroom | pool best 1.711 → assembly best 1.175 Å (asm agent) | 3.048 → 2.281 Å on the 16-hypothesis family (agg agent) |

**Problem A is the primary formulation** because its Hamiltonian is *exactly* the objective
with no surrogate anywhere, and because it is the textbook one-hot 2-local QUBO. Problem B
is run as a **second family** on the aggregation agent's dumped instances, where a certified
exhaustive optimum already exists for all 126 targets, so the quantum arm is scored against
ground truth without re-deriving it.

Both are honest combinatorial problems in the sense the brief asks for. Neither is a
generation lever: the assembly agent measured deployable assembly at **+0.218 Å worse** than
the incumbent synthesis, and the aggregation agent measured its subset objective's exact
minimiser at the **57.8th** RMSD percentile. That is stated up front and is not walked back
anywhere below.

---

## 2. The Hamiltonian, and its verification

### 2a. Formulation (VQ-A: anchored multi-segment assembly)

Variables. A target of length `n` is cut into `k` contiguous segments `s = 0..k-1`, each of
length ≥ 4 (`vq_lib.segments`). Segment `s` chooses one of `m` library pieces,
`x_s ∈ {0..m-1}`.

Candidates. For each segment, the top-200 library windows of that length by BLOSUM62
against the target **sub**-sequence (deployable retrieval, out-of-fold banks from
`s12/asm_lib.bank`), agglomeratively clustered on their `200×200` CA-RMSD matrix into `m`
clusters, keeping each cluster medoid. Clustering is what makes the `m` choices genuinely
distinct shapes rather than 8 copies of the top hit.

Placement. Each piece is superposed (Kabsch) onto the corresponding residues of a
**deployable anchor** — the shipped emitted structure `fit_ca`. No native is read anywhere
in the construction.

Emission. `X(x)` = concatenation of the placed pieces.

Objective. The **shipped** leave-fold-out distogram Bayes risk of `X(x)`, i.e. exactly
`I.shipped_score`:

    F(x) = (1/npairs) Σ_{p=(i,j), j-i≥2} risk_p( ‖X_i(x) − X_j(x)‖ )

### 2b. The decomposition claim

A placed piece's coordinates depend only on its own segment's choice. Every pair `(i,j)`
therefore falls in exactly one block: internal to a segment `s`, or straddling `(s,t)`.
Hence, with **no approximation of any kind**,

    npairs · F(x) = Σ_s h_s(x_s) + Σ_{s<t} J_st(x_s, x_t)

`h_s ∈ R^m`, `J_st ∈ R^{m×m}` are tabulated once per instance by summing the shipped risk
over the pairs of each block.

One-hot QUBO over `N = k·m` binary variables `x_{s,a}`:

    E(x) = Σ_{s,a} h_s(a)/npairs · x_{s,a}
         + Σ_{s<t} Σ_{a,b} J_st(a,b)/npairs · x_{s,a} x_{t,b}
         + P · Σ_s ( Σ_a x_{s,a} − 1 )²      ,   P = 2 · (max F − min F)

and the Ising form `x = (1−z)/2` is `vq_lib.ising_from_qubo`. An **index encoding**
(`⌈log₂ m⌉` qubits per segment, indices taken mod `m`, no constraint, whole register
feasible) is also provided; it is *not* 2-local but it uses the register `k·m / (k·log₂ m)`
times more efficiently, and both are run.

### 2c. Verification — the Hamiltonian equals the real objective

For each instance: 64 random configurations, plus the full `m^k` enumeration, scored two
ways — through the tables, and by rebuilding the assembled structure from scratch and
calling the shipped scorer.

| check | n=9 (1CS9) | n=14 (1A13) | n=16 (1KWE) |
|---|---|---|---|
| max abs(tables − shipped scorer, float32 accumulation) | 2.75e-07 | 4.85e-07 | 1.06e-06 |
| **max abs(tables − truth, float64 accumulation)** | **0.0** | **0.0** | **0.0** |
| **risk-grid bin of every pair identical under both routes** | **True** | **True** | **True** |
| max abs(one-hot QUBO − tables), feasible states | 2.2e-16 | 8.9e-16 | 1.1e-15 |
| max abs(Ising − QUBO), 500 random bitstrings incl. infeasible | 1.6e-13 | 1.7e-13 | 2.3e-13 |

The decomposition is **exact, not approximate**: the float64 route agrees bit-for-bit, and
the ≈5e-07 residual against `I.shipped_score` is entirely that function's float32
accumulation of the 78×760 risk lookup (the *bins addressed are identical*, which is the
whole content of the claim). This is a genuine 2-local Ising Hamiltonian for the real
production objective, not a proxy for it.

---

## 3. Is the interaction structure real, and how hard is the problem? (`s12/vq_run.py`)

Three instance families, every tuning target, **certified exhaustive optima** throughout.
`A2_8` = 2 segments × 8 pieces (64 assemblies), `A3_5` = 3 × 5 (125), `A2_64` = 2 × 64
(4,096). Objective values are the shipped distogram Bayes risk; gaps are to the exhaustive
optimum of the same instance.

### 3a. The interaction is exactly real — and it is small

An exact functional ANOVA of `F` under the uniform product measure on the choices splits
`Var(F)` into one-body and pure-interaction parts. The split is exact by construction
(orthogonality of the centred terms); the residual `|Var(F) − V1 − V2|` measures only
float error.

| family | configs | interaction share of Var(F), mean | median | max | ANOVA residual (max) | mean # 1-change local minima | % instances multimodal |
|---|---|---|---|---|---|---|---|
| A2_8   | 64    | **0.0253** | 0.0203 | 0.115 | 9.6e-17 | 1.10 | 9.5 % |
| A3_5   | 125   | **0.0225** | 0.0191 | 0.099 | 6.1e-17 | 1.06 | 6.5 % |
| A2_64  | 4,096 | **0.0245** | 0.0213 | 0.085 | 8.1e-17 | 1.24 | 19.0 % |

**The value of a choice genuinely does depend on the other choices — but only ≈2.4 % of the
objective's variance is interaction.** The remaining 97.6 % is separable. And the landscape
is nearly unimodal: on 81–94 % of instances there is exactly one 1-change local minimum, so
steepest descent from anywhere is the global optimum.

The interaction is not zero, and the right control proves it: `separable_anova`, the
*optimum of the best additive approximation to F* (not merely `argmin h`, which is a much
weaker control), finds the true optimum on only 51–71 % of instances and carries a mean
objective gap of 0.012–0.018. So the couplings do real work — just not much of it.

### 3b. Classical solvers, quality vs cost (mean over all targets)

| family | solver | mean gap to optimum | % exactly optimal | mean objective evals | median wall (ms) |
|---|---|---|---|---|---|
| **A2_64** (4,096 cfgs) | exhaustive | 0 | 100 % | 4,096 | 0.22 |
| | MILP (Glover linearisation, HiGHS) | 0 | 100 % | — | 1,681 |
| | SA, 1,000 steps | 0.00024 | **98.4 %** | 1,065 | 19.4 |
| | **greedy + 1-opt local search** | **0.00142** | **92.1 %** | **368** | **0.069** |
| | beam-16 | 0.00389 | 92.9 % | 1,088 | 2.16 |
| | SA at matched evals | 0.01504 | 71.4 % | 433 | 7.9 |
| | separable_anova (ignore couplings) | 0.01836 | 50.8 % | 128 | 0.000 |
| | beam-4 | 0.02041 | 69.0 % | 320 | 0.74 |
| | greedy (no local search) | 0.06901 | 34.9 % | 128 | 0.077 |
| | random at matched evals | 0.10737 | 10.3 % | 368 | 0.14 |
| | separable_h (`argmin h` only) | 0.12763 | 19.0 % | 128 | 0.000 |
| A2_8 | greedy + 1-opt | 0.00056 | 97.6 % | 37 | 0.061 |
| A3_5 | greedy + 1-opt | 0.00003 | 97.8 % | 41 | 0.123 |

**This problem is solved by a 70-microsecond local search.** Greedy + 1-opt reaches the
certified optimum on 92 % of the 4,096-assembly instances in 368 evaluations (9 % of the
space) and 0.07 ms; annealing closes the rest. The MILP proves optimality on every instance
but is 24,000× slower than the local search that already found the same answer. Random
search at a matched evaluation count is 75× worse in objective gap — so the problem is not
trivial in the sense of "anything works", it is trivial in the sense of "the right cheap
classical method works".

### 3c. What the optimum is worth in ångströms (evaluation only)

| family | quantity | mean | FAIL18 | other-108 | frac < 2 Å |
|---|---|---|---|---|---|
| A2_64 | ORACLE best assembly in the family | 2.735 | 5.400 | 2.290 | 0.373 |
| | **anchor (`fit_ca`, the shipped synthesis)** | **3.204** | 6.026 | 2.734 | 0.278 |
| | greedy + 1-opt | 3.372 | 6.196 | 2.901 | 0.246 |
| | **exhaustive optimum of F** | **3.378** | 6.193 | 2.909 | 0.254 |
| | separable_anova | 3.416 | 6.213 | 2.949 | 0.222 |
| | random assembly | 3.529 | 6.191 | 3.086 | 0.167 |
| | mean over all 4,096 assemblies | 4.029 | 6.256 | 3.658 | 0.008 |

Paired, exact optimum vs the anchor it perturbs: **d = +0.174 Å, CI [+0.108, +0.242],
42W/84L, drop-top-10 +0.216** (A2_8: +0.166 [+0.103, +0.236] 40/86; A3_5: +0.058
[+0.027, +0.093] 29/64). Every CI excludes zero on the wrong side.

Two honest readings, and they point opposite ways:

1. **The objective is positively aligned here, unlike elsewhere in the record.** The exact
   optimum sits at the **24.8th** percentile of achievable RMSD among the 4,096 assemblies
   (`corr(F, RMSD) = +0.449` overall), against the aggregation agent's 57.8th percentile on
   Problem B and the assembly agent's 74th percentile on FAIL18 in the unanchored setting.
   Anchoring the search on the incumbent structure is what buys that: the family is a local
   neighbourhood of a decent structure, so the score's coarse skill is enough to rank inside
   it. On FAIL18 the alignment collapses to the 45.1st percentile and `corr = +0.119`, which
   is the record's own FAIL18 signature.
2. **Optimising it still loses accuracy.** The whole family lies *above* the anchor: even
   its ORACLE best is 2.735 Å against the anchor's 3.204 Å only on 37 % of targets, and the
   deployable optimum is +0.174 Å worse than doing nothing. Solving the problem harder is
   marginally worse than solving it well (greedy+1-opt 3.372 vs exhaustive 3.378), the
   record's "optimise harder, get worse" law again — here at only +0.006 Å, because the
   objective is not badly anti-aligned in this anchored family, merely useless.

So: **a well-posed, exactly-2-local, genuinely-coupled combinatorial problem that a
microsecond of classical local search solves, and whose solution costs 0.17 Å of accuracy.**

### 3d. It stays easy as the instance grows

| family | configs | % instances with >1 local minimum | greedy+1-opt: mean gap / % optimal / evals / wall | exhaustive wall | MILP wall |
|---|---|---|---|---|---|
| A2_8 | 64 | 9.5 % | 0.00056 / 97.6 % / 37 / 0.06 ms | 0.18 ms | 29 ms |
| A3_5 | 125 | 6.5 % | 0.00003 / 97.8 % / 41 / 0.12 ms | 0.17 ms | 99 ms |
| A2_64 | 4,096 | 19.0 % | 0.00142 / 92.1 % / 368 / 0.07 ms | 0.22 ms | 1,681 ms |
| A3_32 | 32,768 | 19.4 % | 0.00072 / 94.6 % / 309 / 0.14 ms | 1.3 ms | 19,265 ms |
| A2_256 | 65,536 | 27.8 % | 0.00147 / 88.1 % / 1,597 / 0.09 ms | 2.1 ms | (not run) |

Multimodality grows with `m` (9.5 % → 27.8 %) but the search cost does not: greedy + 1-opt
still finds the certified optimum on 88 % of 65,536-assembly instances after touching
**2.4 %** of the space, in 0.09 ms. The exhaustive optimum's accuracy is +0.173 Å worse
than the anchor at m=256 (CI [+0.113, +0.237], 41W/85L) and +0.067 Å at k=3
(CI [+0.041, +0.095], 27W/66L), i.e. the conclusion of §3c is stable across the family.

---

## 4. The quantum arm at matched budget (`s12/vq_qrun.py`, `s12/vq_quantum.py`)

### 4a. Protocol — the project's own shared-budget accounting

The device-realistic arm is `core.quantum.run_global_cvar_vqe`: `lightning.qubit`, the
`layers × (RY on every wire, CNOT chain + ring)` circuit, bitstrings **sampled** from the
measured distribution, the **genuine sample CVaR** as the objective, SPSA as the optimiser.
It and every classical search are handed the SAME `budget.BudgetedEnergyModel`, which
charges one evaluation per **unique** bitstring and caches thereafter — the sprint's own
protocol, not one invented here. The metric is the best energy any arm ever saw, against
the certified exhaustive optimum. Family `A2_64`: 12 qubits, 4,096 assemblies, all 126
targets, alpha = 0.25, layers = 3, shots = 16, restarts = 1 (the setting that lets SPSA take
100+ iterations rather than tripping the module's own budget-starvation warning).

### 4b. Result: the classical solvers win at every budget

| budget (unique evals) | VQE gap | random gap | anneal gap | **greedy+1-opt gap** | VQE % optimal | greedy % optimal | VQE median wall | greedy median wall |
|---|---|---|---|---|---|---|---|---|
| 256   | 0.13213 | 0.10308 | 0.11132 | **0.05463** | 2.4 % | 28.6 % | 76 ms | 1.9 ms |
| 512   | 0.08725 | 0.07102 | 0.07391 | **0.02906** | 10.3 % | 47.6 % | 186 ms | 4.3 ms |
| 1,024 | 0.05734 | 0.04105 | 0.03938 | **0.00973** | 23.8 % | 71.4 % | 636 ms | 10.8 ms |
| 2,048 | 0.02390 | 0.01578 | 0.01851 | **0.00291** | 46.8 % | 88.1 % | 3,566 ms | 46.7 ms |
| 3,200 | 0.01536 | 0.00635 | 0.00631 | **0.00138** | 61.1 % | **97.6 %** | 9,203 ms | 144 ms |

Paired, VQE minus each classical arm (positive = VQE worse), n = 126:

| budget | vs random | vs annealing | vs greedy+1-opt |
|---|---|---|---|
| 256   | +0.02905 [+0.01132, +0.04693] 54W/71L | +0.02081 [+0.00017, +0.04191] 55W/68L | **+0.07751 [+0.06106, +0.09435] 18W/104L** |
| 512   | +0.01623 [+0.00262, +0.03023] 54W/66L | +0.01334 [−0.00238, +0.03012] 55W/60L | **+0.05819 [+0.04333, +0.07360] 20W/92L** |
| 1,024 | +0.01629 [+0.00245, +0.02992] 48W/66L | +0.01796 [+0.00471, +0.03174] 51W/59L | **+0.04761 [+0.03625, +0.05983] 13W/86L** |
| 2,048 | +0.00812 [+0.00064, +0.01595] 43W/55L | +0.00539 [−0.00325, +0.01304] 20W/54L | **+0.02099 [+0.01403, +0.02850] 7W/63L** |
| 3,200 | +0.00901 [+0.00284, +0.01530] 21W/47L | +0.00904 [+0.00278, +0.01570] 21W/45L | **+0.01397 [+0.00929, +0.01916] 0W/63L*** |

*(at 3,200 the VQE never once beat greedy+1-opt on any of the 126 targets; the 48 recorded
losses are the targets where the two differ at all — the rest are exact ties at the optimum.)*

**The classical solvers win, cleanly, at every budget, with every CI excluding zero.**
Greedy + 1-opt reaches the certified optimum on 97.6 % of instances at 3,200 evaluations
and 144 ms; the VQE reaches it on 61.1 % using the same 3,200 unique evaluations and 9.2 s
— **64× the wall time for a 11× larger objective gap**. The VQE also loses to uniform
random sampling at four of five budgets. That is not a surprise and it is not a defect of
this implementation: `core/quantum.py`'s own module docstring records the same fact
("at these problem sizes the VQE ties uniform random sampling and loses to annealing").
It reproduces here, on a different and better-posed Hamiltonian, over 126 instances.

**Why**: the state entropy stays at **9.1 bits of a possible 12** at every budget, from 13
SPSA iterations to 1,800. The circuit distribution never concentrates, so the search *is*
biased random sampling — which is exactly what its gap being between random's and
annealing's says.

### 4c. Gradient audit, re-measured on these Hamiltonians

Not quoted from the record — measured on the VQ-A energy vectors (`vqe_gradient_audit`):

| estimator | cos with exact parameter-shift gradient | ‖g‖ / ‖g_exact‖ |
|---|---|---|
| parameter shift vs. independent finite differences | **+1.0000000** | — |
| sampled score function, `baseline="const"` (the fix, used everywhere here) | **+0.978** | 0.996 |
| sampled score function, `baseline="tail"` (the recorded DEFECT) | **+0.581** | 0.561 |

The defect reproduces independently (record: +0.656 at 0.758×; here +0.581 at 0.561× on a
different Hamiltonian). Every result in this file uses `baseline="const"`, and
`run_cvar_vqe` differentiates the exact parameter-shift gradient, so no number here is
contaminated by it. Separately, `lightning.qubit` and the exact statevector agree to
**1.4e-16** max abs on the same angles (`pennylane_probs_match`), so the fast path and the
pinned device are the same circuit.

---

## 5. Null controls for the one attractive correlation (`s12/vq_null.py`)

§3c found `corr(F, CA-RMSD) = +0.449` and the exact optimum at the 24.8th RMSD percentile —
the shipped objective appearing to rank *inside* this anchored family, which is not what the
record says it does elsewhere. Two nulls, all 126 targets, family A2_64.

**Null 1 — shuffled objective.** The per-pair risk rows of the shipped distogram are
permuted across pairs (identical 17-bin marginals, wrong pair assignment) and the whole
2-local Hamiltonian is rebuilt from the same placed pieces. 3 shuffles per target.

**Null 2 — radius of gyration.** `rg` drives both risk and RMSD; the partial correlation
`corr(F, RMSD | rg)` is the control that overturned "physics ranks real geometry" in the
record. Also reported: what the `|rg − rg_native|` argmin alone achieves (an ORACLE control).

| quantity | real objective | shuffled null | rg-only (oracle) |
|---|---|---|---|
| corr(F, RMSD), mean / FAIL18 / other-108 | **+0.449** / +0.119 / +0.504 | +0.205 / +0.187 / +0.208 | — |
| **partial corr(F, RMSD ∣ rg)** | **+0.333** / +0.217 / +0.352 | **−0.017** / −0.064 / −0.009 | — |
| corr(F, rg) | +0.275 | — | — |
| corr(RMSD, rg) | +0.163 | — | — |
| RMSD percentile of the optimum | **24.8th** | 53.3rd (chance) | 35.2nd |
| mean RMSD of the optimum | **3.378** | 4.055 | 3.759 |
| (best / mean assembly in the family) | 2.735 / 4.029 | | |

Paired, real optimum vs shuffled-null optimum: **d = −0.678 Å, CI [−0.835, −0.531],
104W/22L, drop-top-10 −0.505**, negative in all five folds (−0.40 to −0.99).

**Both nulls are passed, and this is a genuinely positive and slightly surprising result.**
The shipped distogram's ranking skill inside an *anchored* assembly family is real: it is
not the geometry of the family (shuffled null is at chance, 53.3rd percentile), and it is not
compactness (partialling out `rg` leaves +0.333 of the +0.449, while the shuffled null drops
to −0.017). The record's finding that this objective ranks the native near chance was
measured on *unanchored* candidate sets (the assembly agent's best assembly at the 74th score
percentile on FAIL18, the aggregation agent's exact optimum at the 57.8th RMSD percentile).
**Anchoring the search on the incumbent structure restores the objective's discrimination**,
and the FAIL18 split shows why it does not help: on FAIL18 the correlation is +0.119 and the
optimum sits at the 45.1st percentile — chance again, exactly the targets that matter.

This does not rescue the accuracy: the entire anchored family lies *above* the anchor, so a
better rank inside it still buys −0.17 Å of harm (§3c). It does say the objective is not
useless when the hypothesis set is a tight neighbourhood of a decent structure — which is
information the coordinator should have, and is a *generation*-side statement, not a
quantum one.

---

## 6. Does CVaR contribute anything a classical optimiser does not? (`s12/vq_qrun.py` stage 2)

This is the question that matters once the classical solvers have won. The one honest
argument on the record is that **CVaR prevents collapse to the argmin and returns a
DISTRIBUTION over hypotheses rather than a point**, which is the right output shape when
the near-native band holds several structural modes. It is tested directly.

Setup: the EXACT-gradient driver `core.quantum.run_cvar_vqe` (Adam on the audited
parameter-shift gradient of `F = CVaR_alpha − T·H`), 12 qubits, 4,096 assemblies, 50
iterations, 10 stratified targets (FAIL18-proportional), full alpha × T grid. Reference on
this subset: anchor 3.081 Å, exact optimum 3.105 Å, oracle best assembly 2.567 Å.

### 6a. The collapse claim is TRUE

State entropy in bits (max 12) at T = 0, i.e. with no entropy term doing the work:

| alpha | 1.00 | 0.50 | 0.25 | 0.10 | 0.05 |
|---|---|---|---|---|---|
| **entropy at T=0** | **0.51** | 2.37 | 4.56 | 5.73 | **7.27** |
| entropy at T=0.02 | 1.95 | 5.57 | 7.76 | 9.10 | 9.74 |
| entropy at T=0.30 | 11.51 | 11.53 | 11.46 | 11.32 | 11.26 |

`alpha = 1` (the plain expectation) collapses to **0.51 bits** — the record's 0.076 bits
reproduces in kind. Lowering alpha monotonically preserves entropy, reaching 7.27 bits at
alpha = 0.05 with **no entropy term at all**. So the claim is confirmed: CVaR does prevent
the collapse, and it does so on its own.

The mechanism is worth stating precisely, because it is *indifference*, not preference.
`CVaR_alpha` is minimised by ANY distribution with at least `alpha` mass on the argmin —
the upper `1 − alpha` of the distribution is unconstrained. The maximum-entropy exact
minimiser therefore carries
`H_max(alpha) = −α log₂ α − (1−α) log₂((1−α)/(dim−1))` bits: 11.68 at α=0.05, 9.81 at
α=0.25, 7.00 at α=0.50, 0 at α=1 (dim = 4,096). The measured entropies sit *below* that
ceiling but far above zero. CVaR does not seek diversity; it declines to penalise it.

### 6b. …and the surviving breadth carries no usable information

Every arm below is the same 10 targets. `rm_pw` = the distribution readout (p-weighted
coordinate average of all 4,096 assemblies). `rm_gibbs` = **the decisive control**: a
classical Boltzmann distribution `p ∝ exp(−E/T′)` with `T′` chosen so its entropy equals
the VQE state's, read out the same way. `cov8` = best true CA-RMSD among the 8
most-probable states (VQE / Gibbs) or the 8 lowest-energy states.

| alpha | T | H bits | rmsd argmax | rmsd best-in-support | **rmsd p-weighted** | **rmsd Gibbs @ matched H** | cov8 VQE | cov8 Gibbs = cov8 lowest-E | paired (pw − gibbs) |
|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 0.00 | 7.27 | 3.251 | 3.117 | 3.272 | **3.039** | 3.071 | **2.910** | +0.233 [+0.06,+0.42] 3W/7L |
| 0.05 | 0.02 | 9.74 | 3.242 | 3.111 | 3.292 | **3.080** | 3.034 | **2.910** | +0.212 [+0.05,+0.39] 2W/8L |
| 0.05 | 0.10 | 10.65 | 3.194 | 3.065 | 3.294 | **3.086** | 2.990 | **2.910** | +0.208 [+0.09,+0.33] 2W/8L |
| 0.10 | 0.00 | 5.73 | 3.322 | 3.295 | 3.249 | **3.024** | 3.197 | **2.910** | +0.225 [+0.05,+0.42] 3W/7L |
| 0.10 | 0.10 | 10.40 | 3.180 | 3.120 | 3.302 | **3.071** | 2.986 | **2.910** | +0.231 [+0.09,+0.38] 2W/8L |
| 0.25 | 0.00 | 4.56 | 3.363 | 3.175 | 3.202 | **3.001** | 3.230 | **2.910** | +0.201 [+0.05,+0.38] 3W/7L |
| 0.25 | 0.10 | 9.52 | 3.327 | 3.105 | 3.191 | **3.017** | 3.171 | **2.910** | +0.174 [+0.06,+0.31] 3W/7L |
| 0.50 | 0.10 | 9.32 | 3.393 | 3.200 | 3.132 | **3.030** | 3.151 | **2.910** | +0.101 [−0.01,+0.23] 4W/6L |
| 1.00 | 0.00 | 0.51 | 3.535 | 3.583 | 3.507 | **3.098** | 3.184 | **2.910** | +0.408 [+0.08,+0.75] 3W/7L |
| 1.00 | 0.30 | 11.51 | 3.813 | 3.120 | 3.285 | **3.168** | 3.311 | **2.910** | +0.118 [+0.05,+0.19] 1W/9L |

*(all 20 alpha × T arms are in `s12/results/vq_stage2_A2_64.json`; the 10 shown span the grid)*

Three results, and they all point one way.

1. **The distribution readout is worse than the matched-entropy classical one on all 20
   arms**, by +0.10 to +0.41 Å, with the bootstrap CI excluding zero on 16 of 20 and a
   typical 2W/8L. A Boltzmann soft-min — one line of numpy, no circuit — dominates the
   variational state at every entropy the VQE reaches.
2. **The VQE's most-probable states never cover the near-native band better than the lowest
   energy states.** `cov8_VQE` is 2.94–3.31 across all 20 arms; `cov8_Gibbs` and
   `cov8_lowest_energy` are both exactly **2.910** — identical, because a Boltzmann
   distribution's most-probable states *are* the lowest-energy ones. The extra breadth the
   VQE carries points at structures that are neither low-energy nor near-native.
3. **The distribution readout is also worse than the point readout** (`rm_pw` 3.13–3.51 vs
   `rm_best_in_support` 3.07–3.58) and worse than doing nothing (anchor 3.081) on every arm.

**Answer: CVaR does contribute the thing it is claimed to contribute — a non-collapsed
distribution — and that distribution is not worth having here.** The right output shape does
not compensate for content: what the state is broad *over* is chosen by an objective whose
in-band skill is +0.13 (record) and whose optimum is +0.17 Å worse than the incumbent (§3c).
Breadth over an uninformative ranking is breadth over noise, and a classical Gibbs
distribution provides the same breadth with a strictly better ordering.

---

## 7. Problem B — subset selection for the terminal operator (`s12/vq_brun.py`)

The aggregation agent's dumped instances: 16 binary variables (one per structurally-distinct
hypothesis drawn from the shipped top-75), all 126 targets, certified exhaustive optimum.

    H(x) = x'Mx/k² − 2b'x/k + const ,  k = Σx ≥ 2

This is **not a QUBO** — the `1/k` and `1/k²` make it a ratio of quadratics, which is
precisely what makes it non-separable (the risk of an average is not the average of risks).
It is kept exactly; the `k < 2` states are given a penalty of `max + span`. 16 qubits, whole
register otherwise feasible.

**Emission validated first**: the uniform coordinate average of the selected
medoid-superposed candidates reproduces the aggregation agent's own dumped RMSDs to
**max 1.9e-07 Å** over 24 checks, so the operator here is theirs, not a lookalike.
Independent reproductions of their headline numbers: RMSD percentile of the exact optimum
**57.80th** (they report 57.8th), `corr(H, RMSD) = +0.0774` (they report +0.077),
avg75 3.0483, uniform-16 3.0724, oracle best subset 2.2810 — all exact.

| budget | VQE gap | random gap | anneal gap | greedy+1-opt gap | VQE % opt | anneal % opt | VQE wall | anneal wall |
|---|---|---|---|---|---|---|---|---|
| 256 (0.4 % of the space) | 0.13640 | 0.14302 | **0.00365** | 0.00755 | 3.2 % | **84.9 %** | 0.21 s | **0.05 s** |
| 1,024 | 0.09151 | 0.09732 | **0.00066** | 0.00044 | 5.6 % | **96.8 %** | 0.88 s | 0.10 s |
| 4,096 | 0.05044 | 0.06452 | **0.00000** | 0.00001 | 15.1 % | **100 %** | 5.79 s | 0.30 s |
| 16,384 (25 % of the space) | 0.01975 | 0.03089 | **0.00000** | 0.00000 | 32.5 % | **100 %** | 51.45 s | 0.96 s |

Paired, VQE minus classical (positive = VQE worse), n = 126:

| budget | vs random | vs annealing | vs greedy+1-opt |
|---|---|---|---|
| 256 | −0.00662 [−0.0366, +0.0163] 56W/70L | **+0.13275 [+0.1005, +0.1738] 1W/122L** | +0.12885 [+0.0995, +0.1652] 2W/121L |
| 1,024 | −0.00581 [−0.0432, +0.0348] 61W/63L | **+0.09085 [+0.0634, +0.1298] 0W/119L** | +0.09107 [+0.0635, +0.1300] 0W/119L |
| 4,096 | −0.01407 [−0.0372, +0.0033] 62W/58L | **+0.05044 [+0.0374, +0.0651] 0W/107L** | +0.05044 [+0.0374, +0.0651] 0W/107L |
| 16,384 | −0.01115 [−0.0257, −0.0002] 49W/52L | **+0.01975 [+0.0141, +0.0263] 0W/85L** | +0.01975 [+0.0141, +0.0263] 0W/85L |

**Annealing beats the VQE on 122 of 126 targets at the smallest budget and on 100 % of the
targets where they differ at every larger one — it never once wins.** Annealing finds the
certified optimum on 84.9 % of instances after seeing **0.4 %** of the space in 50 ms; the
VQE finds it on 32.5 % after seeing **25 %** of the space in 51 s (a 1,000× wall-time
disadvantage). Against uniform random search the VQE is **statistically indistinguishable**
(three of four CIs cross zero) — the module docstring's own "ties uniform random sampling
and loses to annealing", now measured on 126 instances of a second, independent problem.
State entropy again pins the reason: **12.2–12.3 bits of a possible 16 at every budget**.

The accuracy side confirms the aggregation agent in full and adds nothing: every solver's
emitted RMSD is 3.146–3.199 against avg75's **3.048** and uniform-16's 3.072, and the better
the optimiser the *worse* the structure (VQE-at-256 3.146 < anneal 3.192 = exact optimum
3.199). The exact optimum sits at the 57.8th RMSD percentile — worse than a coin flip. **On
Problem B a better optimiser is measurably harmful.**

---

## 8. End-to-end through the real path (`s12/vq_e2e.py`)

Every arm's assembled point cloud is pushed through the production projection
(`I.project`, λ=0 arm = `fit_ca`), so the numbers sit on the same scale as the shipped
3.2041. Family A2_64, budget 3,200, all 126 targets. `vqe_cvar_tail_avg` is the honest
"return a distribution, not a point" readout: the uniform average of the CVaR-α tail of
everything the budgeted search actually visited.

| arm | raw | **projected** | FAIL18 | other-108 | frac < 2 Å | paired vs anchor | CI95 | W/L | drop-10 | per-fold |
|---|---|---|---|---|---|---|---|---|---|---|
| **anchor (shipped `fit_ca`)** | 3.2041 | **3.2041** | 6.0258 | 2.7338 | 0.278 | — | — | — | — | — |
| vqe_cvar_tail_avg | 3.2471 | 3.3499 | 6.1617 | 2.8813 | 0.238 | **+0.1459** | [+0.1014, +0.1932] | 31/95 | +0.1812 | +0.07…+0.23 |
| greedy + 1-opt | 3.3718 | 3.3687 | 6.1921 | 2.8981 | 0.246 | **+0.1646** | [+0.1003, +0.2305] | 43/83 | +0.2048 | +0.07…+0.28 |
| exact optimum | 3.3779 | 3.3727 | 6.1904 | 2.9031 | 0.254 | **+0.1686** | [+0.1033, +0.2360] | 43/83 | +0.2112 | +0.06…+0.29 |
| VQE (point) | 3.3962 | 3.3885 | 6.1896 | 2.9216 | 0.238 | **+0.1844** | [+0.1166, +0.2538] | 42/84 | +0.2334 | +0.07…+0.29 |

Head-to-head:

| comparison | d | CI95 | W/L |
|---|---|---|---|
| VQE − exact optimum | +0.0158 | [−0.0133, +0.0462] | 18/31 |
| VQE − greedy+1-opt | +0.0198 | [−0.0148, +0.0558] | 20/34 |
| VQE CVaR-tail average − VQE point | −0.0386 | [−0.0942, +0.0156] | 75/51 |

**Every arm is worse than doing nothing, every CI excludes zero, every fold agrees, and
none of it survives dropping the top ten targets (it gets worse).** The quantum answer is
statistically indistinguishable from the exact optimum and from greedy+1-opt in accuracy
(+0.016 / +0.020, both CIs crossing zero) — the optimiser choice does not matter because
the objective's optimum is the wrong place to be. The distribution readout is the least-bad
arm (+0.146 vs +0.184 for the point readout, 75W/51L between them, CI crossing zero), which
is a real but not significant point in favour of the "return a distribution" argument.

Nothing here is a candidate for a dev24 pass: the one arm that is not a null is negative.

---

## 9. The encoding matters, and one-hot is the expensive one

Two encodings of the same instance were built and both are exercised.

| family | encoding | qubits | register states | feasible | feasible fraction | non-zero couplings | mean problem coupling | mean penalty coupling | penalty / problem |
|---|---|---|---|---|---|---|---|---|---|
| A2_5 | one-hot | 10 | 1,024 | 25 | 2.44 % | 45 (all) | 0.324 | 1.478 | **4.8×** |
| A2_8 | one-hot | 16 | 65,536 | 64 | **0.098 %** | 120 (all) | 0.325 | 1.584 | **5.1×** |
| A3_5 | one-hot | 15 | 32,768 | 125 | 0.38 % | 105 (all) | 0.125 | 1.202 | **10.1×** |
| A2_64 | index | 12 | 4,096 | 4,096 | 100 % | (not 2-local) | — | — | — |
| A2_256 | index | 16 | 65,536 | 65,536 | 100 % | (not 2-local) | — | — | — |

The Ising is **fully connected** in every one-hot case (all `N(N−1)/2` couplings non-zero),
which is the genuine article — but the one-hot penalty couplings are **5–10× larger** than
the couplings that encode the actual problem, and 99.9 % of the register is infeasible at
A2_8. The encoding overhead is measured, on 40 stratified targets (`vq_stage3_A2_8`):

| budget | VQE gap | random gap | anneal gap | greedy+1-opt gap | VQE % opt | anneal % opt | VQE wall | anneal wall |
|---|---|---|---|---|---|---|---|---|
| 1,024 | **2.30994** | 0.70884 | **0.04304** | 0.20274 | 0 % | 75 % | 441 ms | 81 ms |
| 4,096 | 0.40826 | 0.36529 | **0.00000** | 0.06271 | 5 % | 100 % | 1,796 ms | 297 ms |
| 16,384 | 0.15231 | 0.16625 | **0.00000** | 0.03031 | 35 % | 100 % | 7,820 ms | 1,213 ms |

Paired at 1,024: VQE − annealing = **+2.267 [+1.637, +2.932], 0W/40L**; VQE − random =
+1.601 [+0.983, +2.255], 13W/26L. **In the textbook one-hot QUBO form the VQE is worse
than uniform random sampling**, because a near-uniform 16-qubit state lands in the feasible
0.098 % about once every thousand shots while the penalty term dominates everything it
does see. Its state entropy is 12.3–12.4 bits of 16 at every budget — it never learns the
constraint. The index encoding (§4) is the one where the VQE is merely mediocre rather than
catastrophic, and it is precisely the encoding that is *not* 2-local. That tension — the
faithful 2-local form is the one the VQE handles worst — is a real result about this
component's placement, not an artefact of tuning.

---

## 10. The collapse control, converged (`s12/vq_collapse.py`)

§6a's entropy-vs-alpha reading is only worth anything if the surviving entropy is a
converged property of the CVaR objective rather than an under-optimised circuit. 4
stratified targets, **T = 0** (no entropy term anywhere), 12 qubits, iterations 50 → 200 →
600, exact parameter-shift gradient throughout. `H_max` is the analytic ceiling
`−α log₂α − (1−α) log₂((1−α)/(dim−1))` — the entropy of the *maximum-entropy exact
minimiser* of CVaR_α.

| alpha | H bits @ 50 it | @ 200 it | @ 600 it | **H_max (analytic)** | p(argmin) @600 | p_top @600 | gap @600 |
|---|---|---|---|---|---|---|---|
| 0.05 | 7.146 | 6.721 | **6.091** | 11.686 | 0.188 | 0.196 | 0.0573 |
| 0.10 | 5.455 | 4.565 | **4.565** | 11.269 | 0.200 | 0.334 | 0.2912 |
| 0.25 | 4.015 | 3.958 | **4.216** | 9.811 | 0.191 | 0.379 | 0.0907 |
| 0.50 | 2.507 | 1.520 | **1.184** | 7.000 | 0.203 | 0.725 | 0.1064 |
| **1.00** | 0.375 | **0.000** | **0.000** | **0.000** | 0.000 | **1.0000** | 0.3061 |

**Confirmed, with the control.** At `alpha = 1` (plain expectation) the state collapses to
**exactly 0.000 bits, p_top = 1.0000**, and stays there under 12× the optimisation — the
analytic prediction. At `alpha ≤ 0.25` the entropy is stable at 4–6 bits under the same 12×;
it is a converged fixed point, not an unfinished run. The record's claim that CVaR prevents
collapse is **true, is a property of the objective and not of the optimiser, and needs no
entropy term to hold**.

Two caveats that belong with it. First, the measured entropy is far *below* `H_max`, so the
variational circuit is not exploiting the degeneracy — it lands somewhere arbitrary between
collapse and the maximum-entropy minimiser. Second, and more damaging: at `alpha = 1` the
collapsed state has `gap = 0.306`, i.e. **it collapses onto a state that is not the optimum**.
Collapse and correctness are separate failures and the VQE has both.

---

## 11. Answers to the questions I was set

**Is there a scientifically real combinatorial optimisation problem in this pipeline where
the value of a choice depends on the other choices?**
**Yes — anchored multi-segment fragment assembly, and it is exactly 2-local.** The shipped
distogram Bayes risk of an assembly decomposes *without approximation* into one-body tables
per segment and pairwise coupling matrices between segments; the float64 route agrees with a
from-scratch recomputation **bit for bit** and addresses identical risk-grid bins. It is the
textbook one-hot Ising/QUBO a VQE is built for, fully connected, over the real production
objective. Problem B (subset selection) is also genuinely non-separable but is a ratio of
quadratics rather than a QUBO, so Problem A is the cleaner formulation.

**How strong is the interaction?** Real but small: **2.4 % of the objective's variance**, by
an exact functional ANOVA (residual 1e-16). The best additive surrogate still misses the
optimum on 30–50 % of instances, so the couplings do work — but 81–94 % of instances have a
single 1-change local minimum.

**How does CVaR-VQE compare against classical solvers at matched budget?**
**It loses, everywhere, with every CI excluding zero.** Two independent problem families,
several budgets each, 126 targets each, the project's own shared unique-evaluation accounting:

| | Problem A (12 q, 4,096) | Problem B (16 q, 65,536) | Problem A one-hot (16 q, 0.098 % feasible) |
|---|---|---|---|
| best classical | greedy+1-opt: **97.6 % optimal** at 3,200 evals, 144 ms | annealing: **84.9 % optimal** at 256 evals (0.4 % of space), 50 ms | annealing: **75 % optimal** at 1,024 evals |
| CVaR-VQE at the same budget | 61.1 % optimal, 9.2 s | 32.5 % optimal at 16,384 evals, 51 s | 0 % optimal, 441 ms |
| paired (VQE − best classical) | +0.0140 [+0.0093, +0.0192] **0W/63L** | +0.0198 [+0.0141, +0.0263] **0W/85L** | +2.267 [+1.637, +2.932] **0W/40L** |
| vs uniform random | VQE **loses** at 4 of 5 budgets | statistically tied | VQE **loses** (+1.60) |
| state entropy | 9.1 / 12 bits at every budget | 12.3 / 16 bits at every budget | 12.4 / 16 bits at every budget |

The mechanism is one number: the circuit distribution never concentrates, so the search is
biased random sampling. This is not an implementation failure — `core/quantum.py`'s own
docstring records "the VQE ties uniform random sampling and loses to annealing", and it
reproduces here on two new problems and 126×2 instances.

**Is there anything CVaR specifically contributes?**
**Yes, exactly one thing, and it is not worth having.** CVaR does prevent collapse: entropy
at T=0 runs 0.000 bits at alpha=1 to 6.1 bits at alpha=0.05, converged under 12x more
optimisation and matching the analytic degeneracy bound. But the surviving breadth is
*indifference*, not preference, and it carries no information:

- the distribution readout loses to a **matched-entropy Boltzmann distribution on all 20
  alpha x T arms** (+0.10 to +0.41 A, CI excluding zero on 16 of 20);
- the VQE's 8 most-probable states **never** cover the near-native band better than the 8
  lowest-energy states (2.94–3.31 A vs 2.910 A, on every arm);
- the distribution readout is worse than the point readout and worse than doing nothing.

**End-to-end?** Through the real projection path, all 126 targets: every arm is worse than
the incumbent — exact optimum **+0.169 A [+0.103, +0.236] 43W/83L**, VQE **+0.184 A
[+0.117, +0.254] 42W/84L**, VQE CVaR-tail distribution **+0.146 A [+0.101, +0.193] 31W/95L**
— every CI excluding zero, every fold agreeing, all worsening on drop-top-10. The VQE and the
exact optimum are statistically indistinguishable in accuracy (+0.016, CI crossing zero),
because the optimiser is not what is wrong.

## 12. Recommendation for the component's role in the final architecture

Do not remove it, do not claim it buys accuracy, and **move it to where the problem is real**.

1. **Re-site the mandated CVaR-VQE from the top-128 selector to the anchored assembly
   Hamiltonian of section 2.** The present placement is indefensible on its own terms: its
   energy is a per-candidate score, so there is no interaction structure and the
   "optimisation" re-expresses what an argmin already has. The assembly Hamiltonian is a
   genuine, verified, fully-connected 2-local Ising over the real production objective. That
   change makes the component scientifically well-posed **at zero accuracy cost either way**
   — both placements measure null, and this one is at least a real instance of the problem
   class the method is for.
2. **Report the classical comparison alongside it, permanently.** The honest headline is
   "greedy + 1-opt solves this in 0.07 ms and 97.6 % of the time; the CVaR-VQE reaches
   61.1 % in 9.2 s under the same evaluation budget." That is a defensible scientific role —
   a genuine quantum method on a genuine instance, benchmarked and losing — and it is far
   stronger than a decorative component nobody can defend.
3. **Keep the CVaR, drop the claim for its distribution.** CVaR earns its place as the thing
   that keeps the state from collapsing (measured, converged, analytically explained), and
   that is a real property worth documenting. Do not claim the distribution is a better
   multi-hypothesis output: a matched-entropy Boltzmann distribution beats it on 20 of 20
   arms, which means the readout can be replaced by `softmax(-E/T)` with no loss and a large
   speed-up. If a multi-hypothesis output is wanted, take it from the Gibbs distribution.
4. **Use the index encoding, not one-hot, if the component must produce an answer**; use
   one-hot when the point is the 2-local Hamiltonian. They are different jobs: one-hot is the
   faithful Ising (and the VQE is catastrophic on it — worse than random); the index encoding
   is where it is merely mediocre, and it is not 2-local. State which one is being shown.
5. **The one genuinely new result here is not quantum** (see section 5): inside an *anchored*
   candidate family the shipped distogram has real ranking skill — `corr = +0.449`, the
   optimum at the 24.8th RMSD percentile, both nulls passed (shuffled objective at chance,
   `+0.333` surviving an `rg` partial). The record's "the objective does not rank" was
   measured on unanchored sets. This is a **generation-side** lead: build hypothesis sets as
   tight perturbations of the incumbent rather than as wide pools. It does not pay here
   because this particular anchored family lies entirely above its anchor — but the ranking
   mechanism works, and on FAIL18 it does not (+0.119, 45th percentile), which localises the
   failure exactly where the record already puts it.
6. **Nothing here earns a dev24 pass.** Every deployable arm is negative or null.

## 13. Negatives, in full

| claim tested | result |
|---|---|
| the assembly objective decomposes exactly into 2-local terms | **CONFIRMED** (bit-exact, identical risk bins) |
| the interaction is a large part of the objective | **refuted** — 2.4 % of the variance |
| the problem is hard for classical solvers | **refuted** — greedy+1-opt, 97.6 % optimal, 0.07 ms |
| CVaR-VQE beats greedy+1-opt at matched evaluations | **refuted**, 0W/63L |
| CVaR-VQE beats annealing (Problem B) at matched evaluations | **refuted**, 0W/85L, never once wins |
| CVaR-VQE beats uniform random sampling | **refuted** on Problem A and on one-hot; tied on B |
| the one-hot (2-local) encoding is workable for the VQE | **refuted** — worse than random, +1.60 at 1,024 evals |
| CVaR prevents collapse to the argmin | **CONFIRMED**, converged, with an analytic bound |
| the CVaR distribution covers the modes better than a point answer | **refuted** — never beats the 8 lowest-energy states |
| a distribution readout beats a point readout | **not shown** — -0.039 A, CI [-0.094, +0.016], 75W/51L |
| the CVaR distribution beats a matched-entropy Boltzmann one | **refuted**, 20 of 20 arms |
| solving the assembly Hamiltonian improves the emitted structure | **refuted** — +0.169 A, CI excludes 0, 43W/83L |
| solving Problem B improves the emitted structure | **refuted** — optimum at the 57.8th RMSD percentile |
| the objective's ranking skill inside the anchored family is an artefact | **refuted** — both nulls passed |

## 14. Leakage audit

- **Deployable inputs only** in every instance and every objective: the target sequence,
  BLOSUM62, the fold's out-of-fold library (`s12/asm_lib.bank`, keyed by (fold, L); the
  assembly agent verified it bit-identical to the shipped universe at L = n), the shipped
  leave-fold-out distogram (`I.distogram`), and the production pipeline's own emitted
  structure `fit_ca` as the anchor. `fit_ca` is an output of the deployable pipeline, so
  anchoring on it introduces no label.
- `nat_ca` enters **only** as an evaluation label (`vq_run.native`, the RMSD columns, the
  ORACLE rows, and the `rg`-control's `rg_native`), never in a solver, an objective, a
  Hamiltonian, a readout or a hyperparameter choice. Every oracle row is labelled ORACLE.
- No `results/benchmark_manifest.json`, no `s9/final_cache/*`, no dev24 target, no
  `esm_cache.npz`. A grep over `s12/vq_*.py` for those paths returns nothing.
- Nothing outside `s12/` was modified; no `git` write command was run. `OMP_NUM_THREADS=2`
  throughout; `I.free_gb()` gates every heavy loop.
- Three independent reproductions of shipped / other-agent numbers: the instrument selfcheck
  (3.4540 / 1.7108 / 2.3062 / 3.2041 / 18); the aggregation agent's emitted RMSDs to
  **1.9e-07 A** and its Problem-B headline numbers (57.8th percentile, corr +0.077, avg75
  3.0483) exactly; and `lightning.qubit` vs the exact statevector to **1.4e-16**.

## 15. Files

| file | what |
|---|---|
| `s12/vq_lib.py` | instance construction, the exact 2-local tables, both encodings, verification |
| `s12/vq_build.py` | builds and caches the six instance families |
| `s12/vq_classical.py` | exhaustive / greedy / 1-opt / beam / SA / MILP, categorical and binary |
| `s12/vq_quantum.py` | CVaR-VQE wrappers, gradient + device audits, shared-budget classical searches |
| `s12/vq_run.py` | stage 1: ANOVA, local-minima count, solver ladder, per-arm RMSD |
| `s12/vq_qrun.py` | stage 2 (alpha/T sweep, Gibbs control) and stage 3 (matched budget) |
| `s12/vq_brun.py` | Problem B on the aggregation agent's dumped instances |
| `s12/vq_null.py` | shuffled-objective and rg null controls |
| `s12/vq_collapse.py` | the collapse/convergence control and the analytic entropy bound |
| `s12/vq_e2e.py` | end-to-end through `I.project` |
| `s12/vq_report.py` | pooled tables and paired statistics |
| `s12/results/vq_stage1_{A2_8,A3_5,A2_64,A3_32,A2_256}.json` | problem structure + classical ladder |
| `s12/results/vq_stage2_A2_64.json` | alpha/T sweep, entropy, readouts, Gibbs control |
| `s12/results/vq_stage3_{A2_64,A2_8}.json`, `vq_stage3_B.json` | matched-budget comparisons |
| `s12/results/vq_null_A2_64.json`, `vq_collapse_A2_64.json`, `vq_e2e_A2_64.json` | controls, end-to-end |
| `s12/results/vq_encoding_stats.json`, `vq_report_*.json`, `vq_*.log` | encoding table, pooled reports, logs |
