# S31 THEORY — lane A (quantum / CVaR)

What the sprint asked lane A: **Q1** what must the state specify to be worth more than one
integer; **Q2** is there a physically meaningful non-diagonal Hamiltonian; **Q3** discharge the
classical-reducibility burden. Mid-sprint the coordinator added **R1** (a readout-capacity
claim to attack) and then redirected to **(a)** realised capacity and **(b)** the ansatz's
inductive bias.

The answers, up front:

1. **Q1.** T1's conclusion ("the tail is a prefix") never fails. What can fail is "*of a fixed
   order*", and it fails only when `grad V` depends on `λ`. Making `H` non-diagonal does **not**
   do that under either definition of CVaR that is well posed for a non-diagonal `H` (§2,
   Theorem A2). Adding the entropy term does not do it either: it replaces "one integer" with
   **one monotone curve from a two-parameter family**, pinned by **one scalar** (§1, A1) —
   which is the same answer lane L reached by duality.
2. **Q2.** There **is** a forced non-diagonal object, and it is not the one the brief proposed
   and not a Hamiltonian. The exact error of the weighted-average readout is
   `⟨w,a⟩ − ½w'Bw` (§3, A3), so the native-free half of the objective is **pairwise with a
   minus sign** — repulsive, where the brief proposed attractive. The corresponding operator
   `H[w] = diag(â) − B` is **mean-field**, quartic in ψ, with no per-shot eigenvalue: it cannot
   be a CVaR-VQE, which is lane L's closure reached along a second road.
3. **Q3.** Everything in this lane is classically reducible and the costs are small (§5).
   Stated plainly and without hedging, because §11 of the charter requires it.
4. **R1** half-survives and the surviving half is stronger than stated (§4).
5. **(b)** The ansatz is an **MPS Born machine of bond dimension exactly 2^layers** over the
   *bits of the DIS rank*, while the objective's optimum is a **hinge in the rank's value**
   (§6).

Endpoint basis: **CA point cloud** everywhere in this document (production **3.0483 Å**).
Nothing here is measured on the built chain (production 3.2105 Å). ORACLE quantities are
labelled ORACLE and are **NOT DEPLOYABLE**.

---

## 1. THEOREM A1 — the deployed state is one scalar

> **A1.** Let `F(p) = CVaR_α(E;p) − T·H(p)` on `Δ^(D−1)`, `T > 0`, `α ∈ (0,1]`. Then `F` is
> strictly convex and its unique minimiser is
>
>     p*_x  ∝  exp( (μ − E_x)_+ / (αT) )
>
> with `μ` the unique scalar satisfying `Σ_{x : E_x < μ} p*_x = α`. At `α = 1`, `μ > max E`
> and the law collapses to the plain Boltzmann `p*_x ∝ exp(−E_x/T)`.

*Proof.* By Rockafellar–Uryasev,
`CVaR_α(E;p) = max_μ [ μ − (1/α) Σ_x p_x (μ − E_x)_+ ]`, a pointwise maximum of functions
**affine in `p`**, hence convex in `p`; `−T·H` is strictly convex; so `F` is strictly convex and
the minimiser is unique. The inner maximiser is the α-quantile, and by the envelope theorem
`∂CVaR/∂p_x = −(μ − E_x)_+/α`. Simplex stationarity `∂F/∂p_x = const` gives
`T log p_x = const + (μ − E_x)_+`, i.e. the stated law; `μ` is then fixed by the quantile
condition, which is one monotone scalar equation. ∎

**Why the scan over the atoms of `E` fails, which is worth recording because I made the
mistake.** The optimum generically sits at the **kink** where the mass strictly below the
threshold equals `α` *exactly*; there the inner maximiser `μ` ranges over an interval between
two atoms and is a genuinely continuous parameter. Restricting `μ` to values of `E` — the
natural first guess, and what `cvar_exact`'s `q` returns — lands **1.5e-3 high in `p`** and
strictly above the optimum in `F`. The bisection on `μ` is exact.

### 1.1 What A1 says about the architecture

* **The deployed quantum state's entire content is a two-parameter monotone reweighting of one
  classical sort.** `p*` is a non-increasing function of `E` alone, so every level set of `p*`
  is a prefix of the `E`-order. T1's "one integer `m`" is the `T → 0` limit; `T > 0` buys a
  curve, not a direction.
* **The irony belongs in the record.** The entropy term exists (`core/quantum.py:993-1002`) to
  break the degeneracy T1 creates. It does break it — and the unique point it selects is still
  a closed-form function of the same scalar energies.
* **Three of the five pinned folds run at `α = 1`** (`VQE_LFO`, `core/pipeline.py:113`), where
  `CVaR_1 = ⟨E,p⟩` identically. On those folds **the CVaR is doing nothing at all** and the
  objective is a plain Boltzmann free energy. Agreement with a mirror-descent solver there:
  **1.4e-15**.

### 1.2 Measured against the circuit (n = 126, `s31/results/s31_A_r1.json`)

| quantity | value |
|---|---|
| `F(closed form) − F(mirror descent)`, max over 126 | **+1.8e-15** (the closed form is never worse) |
| `KL(p_θ ‖ p*)` mean / median / p90, bits | **0.930 / 1.298 / 1.303** |
| total variation `p_θ` vs `p*` | **0.378** |
| KKT residual of `p*` / of `p_θ` (median) | 0.137 (kink artefact, one state) / **0.416** |
| cloud RMSD, `p_θ`-weighted synthesis − `p*`-weighted | **−0.0044 Å, SE 0.0177, 0.088× MDE — NOT A RESULT** |
| the same, per target in absolute value | 0.102 Å |

> **The circuit misses the objective's optimum by ~1 bit of KL and the readout cannot tell.**
> Registered clause A1-e predicted `KL < 0.10` bits: **FALSIFIED**, by roughly 10×. The second
> half of the same clause — that the two readouts agree to 0.02 Å — **held**, at 0.088× MDE.
> Both halves are informative and they point in opposite directions, which is the finding.

**A registered clause I got wrong, recorded rather than quietly rewritten.** A1-v was
registered against `max|p_closed − p_numeric| < 1e-8`. That statistic measures my **reference
solver's** convergence, not the closed form's correctness, and the measured 3.6e-4 refutes the
clause as written while the closed form is provably the better point. The verification of
record is the `F`-comparison above.

---

## 2. THEOREM A2 — non-diagonal `H` does not make the order endogenous

T1 (S30) concludes that at any KKT point of `min_{Λ(p)} V` the tail is a prefix of the order
induced by `grad V(λ*)`. **Its conclusion never fails.** The load-bearing hypothesis is not
"`H` is diagonal" — it is that `grad V` is *constant in `λ`*, so that the order is exogenous.
The question is therefore whether a non-diagonal `H` makes `grad V` depend on `λ`. There are
exactly three ways to define CVaR for a non-diagonal `H`, and two of them do not.

**(N1) Computational-basis energies, `E_x = ⟨x|H|x⟩ = H_xx`.** The off-diagonal is invisible to
the objective. Reduces to the diagonal case exactly.

**(N2) Spectral CVaR.** Measure `H` projectively; the outcome is `ε_k` with probability
`q_k = ⟨ψ|Π_k|ψ⟩`. Then:

> **A2.** (i) `CVaR_α(H;ψ)` depends on `ψ` only through the spectral weights `q`; (ii) the tail
> is the α-prefix of the **fixed** eigenvalue order, so A1 applies verbatim with `E → spec(H)`;
> (iii) `min_ψ CVaR_α(H;ψ) = ε_min` for **every** `α`, so CVaR cannot change *what* is optimal,
> only the landscape — but its **argmin set is strictly larger**, namely
> `{ψ : ⟨ψ|Π_min|ψ⟩ ≥ α}`; (iv) the map `ψ ↦ q` has fibres of dimension ≥ D−1 (relative
> phases), along which the computational-basis readout `p_x = |⟨x|ψ⟩|²` varies freely.

(iii) is why CVaR helps a shallow ansatz at all — it relaxes "be the ground state" to "have α of
the mass there". (iv) is the damaging part: **the extra degrees of freedom a non-diagonal `H`
creates are exactly the ones the objective is blind to and the readout consumes.** At `T = 0`
the emitted structure on the argmin set is decided by the optimiser's path, not the objective.
This is a mechanism for S30 item L's "nearly rank-one hopping term with negligible gradient
contribution" that does not require the term to have been badly chosen.

**(N3) Local energy, `E_loc(x) = ⟨x|H|ψ⟩/⟨x|ψ⟩`.** This is the only definition under which the
off-diagonal enters a per-outcome energy, and it **does** make the order endogenous:
`E_loc` depends on `ψ`. It is also not a measurement of an observable — it is the variational
Monte-Carlo local energy, and `⟨H⟩ = Σ_x p_x E_loc(x)` is an identity, not a sampling rule.
A CVaR of `E_loc` is a well-defined classical functional of `ψ` and an ill-defined quantum
one: **no projective measurement produces `E_loc(x)` as an outcome**, so there is no CVaR-VQE
here, only a nonlinear classical optimisation wearing a wavefunction.

> **Q1, answered.** The escape from "one fixed order" requires a `λ`- or `ψ`-dependent gradient.
> A non-diagonal `H` supplies one **only** under (N3), which is not a quantum measurement. Under
> (N1) and (N2) a **generalised prefix theorem** applies and the direction is closed. This is
> the same closure lane L reached from the per-shot-eigenvalue side (Barkoutsos et al.,
> Quantum 4:256) and from dimension counting; two independent arguments, same verdict.

---

## 3. THEOREM A3 — the readout's objective is an exact identity, and it forces the sign

> **A3.** Let `W_x` be candidate coordinates in any fixed frame, `t` the native in the same
> frame, `a_x = ‖W_x − t‖²_F`, `B_xy = ‖W_x − W_y‖²_F`. Then for **every** `w` with
> `Σ_x w_x = 1` — non-negativity is **not** required —
>
>     ‖ Σ_x w_x W_x − t ‖²_F   =   ⟨w, a⟩  −  ½ w' B w .

*Proof.* `‖Σ w_x e_x‖² = Σ_{xy} w_x w_y ⟨e_x,e_y⟩` with `e_x = W_x − t`; substitute
`⟨e_x,e_y⟩ = (‖e_x‖² + ‖e_y‖² − ‖e_x − e_y‖²)/2` and use `e_x − e_y = W_x − W_y` and
`Σ w = 1`. ∎ Measured max relative error over 126 targets × 80 draws each, on simplex draws
and on affine draws with negative weights: **7.7e-14**.

Note `½ w'Bw = tr Σ_w`, the weighted dispersion about the weighted mean, so the identity is the
general-`w` form of the project's S23-L9 identity (`mean|e|² = |ē|² + mean|d|²`), which is its
uniform-weight special case.

### 3.1 Four consequences, in decreasing order of how much they hurt

1. **The objective is fully known and half of it is free.** `B` is native-free and exact; `a` is
   the only unknown. **The readout problem and the ranking problem are the same problem** — not
   two components to be designed separately, which is how the pipeline treats them.
2. **The native-free half is pairwise and carries a MINUS sign.** At fixed quality the readout
   should **maximise** weighted mutual spread. The brief's candidate
   `H = diag(zrank) − λ·W(similarity)` with `λ > 0` is **attractive**; the derivation says
   **repulsive**. The sign is wrong.
3. **The attractive branch cannot produce a distribution at all.** `B` is a squared-distance
   matrix, hence conditionally negative definite, so `w ↦ w'Bw` is **concave on the simplex**.
   Therefore `⟨w,â⟩ − γ·½w'Bw` is convex for `γ > 0` and **concave for `γ < 0`**, and a concave
   function on a polytope attains its minimum at a **vertex**. The consensus-attraction branch
   degenerates to "select the single best-scoring candidate" — the shipped argmin. A theorem,
   not a measurement, and it is why no amount of tuning could have rescued that sign.
4. **The forced operator is mean-field, not a Hamiltonian.** `∂/∂w_x [⟨w,â⟩ − ½w'Bw] = â_x −
   (Bw)_x`, which is the VMC local energy of `H[w] = diag(â) − B`. The energy functional is
   **quartic in ψ** (`w = |ψ|²`), so it is a Gross–Pitaevskii-type nonlinear problem with no
   per-shot eigenvalue. Independent confirmation of §2's (N3) closure from the readout side.

### 3.2 The class ladder, which reprices charter §7C

`min_w ‖Σ w_x W_x − t‖²` is a convex QP (ORACLE), and its value depends entirely on the
constraint set:

| readout class | constraint | ORACLE ceiling |
|---|---|---|
| selection (readout 1) | `w` a vertex | best pool member, 2.1458 Å |
| convex (uniform average, `p_θ`-weighted) | `w ∈ Δ`, min-norm point of `conv{e_x}` | §7 table |
| **affine (readout 2, `s27/s28_A_amp.py:105-117`)** | `Σw = 1`, **negatives allowed** | **0 Å whenever `rank(aff{W_x}) ≥ 3n−3`** |

The last row is a rank argument, not an experiment: with `D = 128` candidates spanning a
`(3n−3)`-dimensional centred configuration space (`3n−3 = 36-45` here), the affine hull
generically contains the native, so the class **interpolates exactly**. Measured on 1A13:
`rank = 36`, ORACLE affine residual **0.0000 Å**, with total negative weight mass 24.3.

> **So "0.2516 Å with an ORACLE objective" is not a statement about readout 2's ceiling — its
> ceiling is zero.** Readout 2 is an **over-parameterised interpolator**: 127 free weights
> against ~39 residual dimensions. Its behaviour is decided entirely by the regulariser, and
> the simplex constraint *is* the regulariser the uniform average enjoys. That is the exact
> mechanism behind the coordinator's S31-L2(2) — "expressivity without an aligned objective is
> harmful" — and it predicts that result rather than observing it.

---

## 4. R1, attacked

**Survives, strengthened.** `argmin(P e_j) = j` requires `P[i,j] > 0` for `i ≠ j`. Duplicate
candidates give `P[i,j] = 0` off-diagonal, and `np.argmin` returns the first index, so those
vertices are unreachable. The number of reachable vertices equals the number of byte-distinct
structures **on every target**, so the mechanism is confirmed rather than inferred.

**Falsified for the emitted structure.** The cited lines 869-871 are `quantum_stage`'s
**selector**. `average_weighted` (`core/pipeline.py:880-895`) uses the medoid only as the
superposition **frame** (line 890) and emits a continuous convex combination at line **894**.
`p` therefore enters the structure through a ≤7-bit piecewise-constant frame **and** a
continuous `(D−1)`-dimensional weight vector; only the frame is capped.

Numbers in §7.

---

## 5. THE CLASSICAL-REDUCIBILITY BURDEN (charter §11), discharged without hedging

| construction | classical algorithm that reproduces it | cost |
|---|---|---|
| deployed `CVaR_α − T·H` optimum | A1's closed form: one bisection on `μ`, one `exp` | `O(D log(1/ε))`, microseconds |
| deployed CVaR at `T = 0` | one sort (T1) | `O(D log D)` |
| any diagonal `H` on a candidate register | sort the diagonal | `O(D log D)` |
| any **non-diagonal** `H` on a candidate register, spectral CVaR (§2 N2) | `eigh` of a `D×D` matrix | `O(D³)`, ≈1 ms at `D = 512` |
| `H[w] = diag(â) − B` mean field (§3) | Frank–Wolfe / FISTA on a convex QP | `O(D²)` per iteration |
| MEB weights (`argmax_Δ ½w'Bw`) | minimum-enclosing-ball dual, Frank–Wolfe | `O(D²)` per iteration |
| the readout's ORACLE ceiling | convex QP (simplex) or `lstsq` (affine) | `O(D³)` once |

> **Nothing in this lane is a quantum mechanism.** The reason is structural and it is lane L's,
> stated here in my own terms because it is the single most important architectural fact:
> **a candidate-index register has Hilbert dimension equal to the number of candidates**, so
> every operator on it is a `D×D` matrix and every spectral question about it is a call to
> `eigh`. Index encoding cannot be classically hard *by construction*, whatever the Hamiltonian.
> This is not a criticism of the Hamiltonian; it is a criticism of the encoding, and it applies
> to every Hamiltonian the encoding admits.
>
> **Classically reducible is not the same as useless.** A1's closed form is a strict improvement
> over the shipped solver — same objective, exactly optimal, microseconds instead of seconds.
> §3's identity is a genuinely better objective for the continuous readout than anything the
> pipeline currently minimises. Both should be judged on the endpoint, and neither should ever
> be described as quantum.

---

## 6. (b) THE INDUCTIVE BIAS OF THE REACHABLE SET

The circuit is `layers × (RY on every wire, CNOT chain + ring)` on a real statevector
(`core/quantum.py:818-895`), `n = 7`, `layers = 3`, **21 parameters** against a 127-dimensional
simplex. Three structural facts:

1. **It is an MPS Born machine of bond dimension exactly `2^layers`.** Measured Schmidt ranks of
   `ψ(θ)` across the six contiguous cuts: **2, 4, 8, 8, 4, 2** — capped by `2^L = 8`, which is
   the light-cone bound for a depth-3 nearest-neighbour circuit. This pins lane L's Han et al.
   (PRX 8:031012) framing to this specific ansatz rather than to the family.
2. **The index is the DIS rank.** `o = top[:dim]` with `top` the DIS-sorted order, so basis
   state `x` *is* rank `x` and qubit `b` is bit `b` of that rank.
3. **The objective's optimum is a hinge in the rank's VALUE** (A1: `log p* = (μ−E_x)_+/(αT)`),
   while a product state's log-probability is `Σ_b log q_b(x_b)`, i.e. **additive over bits** —
   and with equal per-qubit odds it is a function of `popcount(x)` alone. The two live in
   different bases.

Measurements in §7. The registered intervention — relabel candidates so rank `i` goes to the
`i`-th bitstring in `(popcount, value)` order, which is **native-free** (it uses only `E`) —
tests fact 3 directly by putting the objective into the ansatz's basis.

---

## 7. MEASUREMENTS

*(filled from `s31/results/s31_A_r1.json`, `s31_A_cap.json`, `s31_A_readout.json`)*

---

## 8. REGISTERED PREDICTIONS AND THEIR OUTCOMES

*(scored in §8 of the ledger entry)*
