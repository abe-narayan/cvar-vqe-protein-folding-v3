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

n = 126, CA point-cloud basis, fold-clustered SE on the pinned 5 folds, MDE = 2.8016 × SE.
Artefacts: `s31/results/s31_A_r1.json`, `s31_A_cap.json`, `s31_A_readout.json` and the three
`*_rows.jsonl` beside them. Data path validated: my reconstruction of production's uniform
DIS-top-75 average gives **3.048338 Å** against the canonical **3.0483**, and the shipped
score argmin gives **3.4540** against the 3.454 quoted in `core/quantum.py`.

### 7.1 The deployed quantum stage, on the cloud basis

| comparison | mean | SE | × MDE | verdict |
|---|---|---|---|---|
| quantum synthesis − PROD75 | +0.0178 | 0.0170 | 0.37 | **NOT A RESULT** |
| quantum synthesis − uniform-128 (matched set) | +0.0125 | 0.0165 | 0.27 | **NOT A RESULT** |
| uniform-128 − PROD75 | +0.0053 | 0.0152 | 0.12 | **NOT A RESULT** |
| sel(p_θ) − sel(uniform medoid) | −0.0308 | 0.0321 | 0.34 | **NOT A RESULT** |
| quantum synthesis − `p*` synthesis | −0.0044 | 0.0177 | 0.09 | **NOT A RESULT** |

> **On the CA cloud at n = 126, every comparison involving the deployed quantum stage is below
> 0.4× MDE.** It is not measurably better or worse than the classical alternatives it competes
> with, including the closed form that replaces it. I earlier quoted the +0.0178 and +0.0125 to
> the coordinator without their MDEs, which made them look like measured costs. They are not.

### 7.2 R1 and realised capacity

| quantity | value |
|---|---|
| byte-distinct candidates in the top-128 | **118.45 mean, 94 min** |
| reachable vertices of `argmin(P e_j)` | **118.45 mean, 94 min** — equal on every target |
| alphabet capacity | **6.886 bits mean, 6.555 worst** (not 7) |
| emitted structure to nearest pool member | **1.1144 Å mean, 0.0730 Å min** |
| `sel(p_θ)` agrees with `sel(uniform)` | 22.2 %, i.e. the weighting moves it on **77.8 %** |

ORACLE bits delivered (S30 currency, `k − E[log₂ rank]`), **ORACLE — NOT DEPLOYABLE**:

| selector | bits |
|---|---|
| `p_θ`-weighted medoid | **1.805** |
| `p*`-weighted medoid | 1.876 |
| score argmin (no quantum stage) | 1.693 |
| uniform medoid | 1.659 |

> The quantum stage moves the selection on 78 % of targets and delivers **0.112 bits more than
> the plain score argmin**, against an alphabet of 6.886 — **1.6 % of the register's capacity**.
> And the closed form `p*` delivers more of them than the circuit does.

### 7.3 The readout identity and its arms

Identity max relative error over 126 × 80 draws (simplex **and** affine with negatives):
**1.66e-11**. Deployed `Pt` substituted for `B`: **2.39 % median relative error** — the
distinction is real and `Pt` is not licensed by the derivation.

| arm | mean cloud RMSD | status |
|---|---|---|
| shipped score argmin | 3.4540 | deployable |
| `GAM(LFO)` | 3.1942 | deployable, one LFO scalar |
| **`MEB` = `argmax_Δ ½w'Bw`, quality-blind** | **3.1919** | **deployable, ZERO parameters, no score** |
| `CAL` (LFO-calibrated `â`, γ = 1) | 3.1817 | deployable, one LFO scalar |
| best cell of the whole γ grid (γ = 2) | 3.1687 | not a deployable selection rule |
| **`PROD75`** | **3.0483** | production |
| ORACLE convex QP over the simplex | **1.8290** (support 6.54 / 128) | **ORACLE — NOT DEPLOYABLE** |
| ORACLE affine (readout 2) | **0.0000** (rank 32.9) | **ORACLE — NOT DEPLOYABLE** |

| comparison | mean | SE | × MDE | folds | verdict |
|---|---|---|---|---|---|
| **PRIMARY `CAL − PROD75`** | **+0.1334** | 0.0409 | **1.17** | 5/5 | **WORSE** |
| `GAM(LFO) − PROD75` | +0.1459 | 0.0471 | 1.10 | 5/5 | WORSE |
| `MEB − PROD75` | +0.1436 | 0.0420 | 1.22 | 5/5 | WORSE |
| **`MEB − argmin(score)`** | **−0.2621** | 0.0342 | **2.74** | 5/5, 81W/45L | **BETTER** |
| `γ=1 − shuffled-B (8 draws)` | −0.0250 | 0.0165 | 0.54 | — | **NOT A RESULT** |
| `γ=1 − shuffled-score` | −0.0500 | 0.0306 | 0.58 | — | **NOT A RESULT** |
| ORACLE convex QP − PROD75 | −1.2193 | 0.0454 | 9.59 | 5/5, 126W | ORACLE |
| ORACLE affine − PROD75 | −3.0483 | 0.0907 | 12.00 | 5/5, 126W | ORACLE |

> **The one positive, and its own control demolishing the reason for it.** Pure dispersion
> maximisation with **no score at all** beats the shipped argmin selector by **−0.2621 Å at
> 2.74× MDE, 5/5 folds**. But the **shuffled-B control fires**: replacing `B` with a random
> relabelling of itself costs only 0.0250 Å at **0.54× MDE — NOT A RESULT**. So the mechanism
> is **"spread the weights over many candidates"**, not **"spread them along the real
> geometry"**. `B`'s *content* is not being used, only its effect on the support's size. That
> is `operator-consumes-set-mean` arriving once more, and **no part of the −0.2621 Å may be
> attributed to the pairwise structure the derivation is about.**

Within the derived family the score is worth 0.050 Å (0.58× MDE, NOT A RESULT), and the whole
family — quality-blind or not — sits **0.12 to 0.15 Å below production**, measured.

### 7.4 The price of the quality estimate (ORACLE sweep — NOT DEPLOYABLE)

`â` interpolated from the shipped DIS z-rank toward the true `a`, γ at the derived value 1 in
Å²:

| ρ(`â`, `a`) | 0.118 | 0.330 | 0.604 | 0.853 | 0.974 | 1.000 |
|---|---|---|---|---|---|---|
| cloud RMSD | 3.170 | 2.893 | 2.614 | 2.332 | 2.027 | **1.829** |

> **The shipped DIS score supplies ρ = 0.1176.** The derived readout crosses production at
> **ρ = 0.211**, reaches **3.00 Å cloud at ρ = 0.248**, and **2.50 Å cloud at ρ = 0.705**.
> Beating production costs a factor of **1.79 in ρ (3.2× in ρ²)**; 3.00 Å cloud costs
> **2.11× in ρ (4.4× in ρ²)**.

This is a **different `ρ`** from S30 THEORY §8.3's `cos(u, e)`. The numerical proximity of
0.1176 to that document's 0.1128 is a coincidence of two different quantities and must not be
quoted as agreement between two instruments.

### 7.5 The inductive bias (the coordinator's redirect (b))

| quantity | value |
|---|---|
| Schmidt ranks of `ψ(θ)` across the six contiguous cuts | **2, 4, 8, 8, 4, 2** — capped at `2^layers` |
| R² of `log p*` on the hinge `(μ − E)₊` | **1.000** (by construction — the yardstick) |
| R² of `log p_θ` on the hinge | **0.467** |
| R² of `log p_θ` on popcount(x) | 0.122 |
| R² of `log p_θ` on hinge + popcount jointly | 0.483 |
| R² of `log p` on the hinge, **random θ from the init law**, mean of 400 | **0.0071** |
| the same, **best of 400** | **0.0538** |
| entropy at random θ | 4.889 bits |

> **In one sentence: the reachable set is an MPS Born machine of bond dimension `2^layers = 8`
> over the bits of the DIS rank; it starts essentially orthogonal to the objective's shape
> (R² 0.007, best-of-400 0.054), optimisation carries it to R² 0.467, and it stops there.**
> That is not a bias *toward* anything structural — it is a **ceiling at about half the right
> shape**, which is the same statement as lane L's expressivity floor in a different currency.

**A hypothesis of mine, refuted.** I predicted the residual was a *basis* mismatch: a product
state's log-probability is additive over bits and, at equal per-qubit odds, a function of
`popcount(x)` alone, whereas `log p*` is a hinge in the rank's *value*. The registered
native-free intervention — relabel so rank `i` goes to the `i`-th bitstring in
`(popcount, value)` order, putting the objective into the ansatz's own basis — **made it
slightly worse**: KL to `p*` 0.967 against 0.930, cloud RMSD **−0.0031 Å at 0.060× MDE, NOT A
RESULT**. Popcount adds only 0.016 of R² beyond the hinge. **The residual is not a
popcount-versus-value basis mismatch.** I do not have a replacement explanation and would
rather say so than fit one after the fact.

### 7.6 The `core/quantum.py` entropy-collapse claim, corrected on the real instrument

The docstring asserted that "for ANY alpha the minimiser concentrates p on the lowest-energy
basis states … state entropy 0.01 bits at alpha=1 … and it is a property of CVaR, not of the
optimiser." Measured at `T = 0`, 8 seeds per target, at the deployed per-fold α:

| | n | entropy at `T = 0` | seed sd of cloud RMSD |
|---|---|---|---|
| `α = 1` (folds 0, 3, 4) | 78 | **0.258 bits** | 0.1634 Å |
| `α = 0.25` (folds 1, 2) | 48 | **3.596 bits** | 0.1599 Å |

At `α = 1`, `CVaR₁ = ⟨E,p⟩` and the minimiser is the unique argmin vertex, so the collapse is
real and **is** a property of CVaR. At `α < 1` the argmin set is `{p : p_x₀ ≥ α}`, a face of
**positive volume**, so the objective does not determine `p` at all and the optimiser's path
picks the point — **3.596 bits, no collapse, and the entropy that exists belongs to the
optimiser**, which is the opposite of the sentence. The quantifier "for ANY alpha" is false by
theorem. Corrected in place, with the original quoted.

---

## 8. REGISTERED PREDICTIONS AND THEIR OUTCOMES

| clause | registered | measured | verdict |
|---|---|---|---|
| A1-v | `max abs(closed − numeric) < 1e-8` | 3.6e-4 | **REFUTED AS WRITTEN** — the statistic measures my reference *solver*, not the closed form. The correct check, `F(closed) ≤ F(numeric)`, holds at **+1.8e-15**. My clause, my error. |
| A1-e (KL) | `KL(p_θ ‖ p*) < 0.10` bits | **0.930** | **FALSIFIED, ~10×** |
| A1-e (readout) | the two readouts within 0.02 Å | −0.0044, 0.09× MDE | **HELD** |
| A1-d (`T = 0` seed sd) | `> 0.15 Å` | 0.162 | **HELD** |
| A1-d (deployed `T` seed sd) | `< 0.05 Å` | 0.096 | **REFUTED, 1.9×** |
| A3-i | identity max rel err `< 1e-10` | 1.66e-11 | **HELD** |
| A3-frame | deployed `Pt` breaks it by `> 1 %` | 2.39 % | **HELD** |
| **A3 PRIMARY** | `CAL − PROD75` in `[−0.15, +0.10]` | **+0.1334, 1.17× MDE, WORSE** | **MISSED** by 0.033 — and it is the direction I said in §6 of the prereg I would report plainly rather than search past |
| A3-sign | LFO γ `> 0` on `≥ 4` of 5 folds | **5 of 5** (5.0, 2.0, 3.0, 1.5, 1.0) | **HELD**, and the negative branch is closed by theorem rather than by the grid |
| A3-meb | `MEB − PROD75` in `[+0.3, +1.5]` WORSE | +0.1436 WORSE | **direction HELD, magnitude MISSED** (2× too small) |
| A3-oracle | `QPORACLE < 1.2 Å` | **1.829** | **REFUTED** by 0.63 Å |
| (b) popcount relabel | the residual bias is popcount-versus-value | −0.0031 Å, 0.060× MDE; KL worse | **REFUTED** |

**Direction of my misses, recorded because it is the only thing that makes the rest
trustworthy.** Two ran *against* my hypothesis (the primary landed worse than my band; A1-e's
KL was 10× my prediction). Two ran *toward* it: I put the convex ORACLE ceiling 0.63 Å too low,
and I over-predicted how badly the quality-blind MEB arm would do, which flattered the
derivation by making its failure look inevitable. `A3-frame` is the only clause where being
right cost me something — it forbids substituting the deployed `Pt` for `B`, which would have
been convenient.

---

## 9. WHAT THIS LANE CLOSES, AND THE ONE THING IT OPENS

**Closed.**

1. **Q1** — the escape from a fixed order needs a `λ`- or `ψ`-dependent gradient; non-diagonal
   `H` supplies one only under the VMC local energy, which is not a measurement (§2).
2. **Q2 as a quantum question** — the forced operator is mean-field and quartic in `ψ`, with no
   per-shot eigenvalue (§3.1 item 4); independently, a candidate-index register's Hilbert
   dimension *is* the candidate count, so every operator on it is an `eigh` away (§5).
3. **Q2's proposed sign** — the attractive/consensus branch is concave on the simplex and
   therefore degenerates to the shipped argmin, by theorem (§3.1 item 3).
4. **The derived objective as a deployable readout** — `CAL`, `GAM` and `MEB` are all
   **0.12–0.15 Å worse than production at 1.1–1.2× MDE, 5/5 folds**, and the shuffled-B control
   shows the pairwise content is not what produces even the part that works (§7.3).
5. **My own popcount explanation of the ansatz gap** (§7.5).

**Open, and it is the only thing in this lane worth compute next sprint.**

> The exact objective for every averaging readout is `⟨w,a⟩ − ½w'Bw`. **Half of it is free and
> exact. The entire deficit is `a`.** The convex readout's ORACLE ceiling is **1.829 Å** against
> production's 3.048, and the crossing price is **ρ(â,a) = 0.211** against the shipped score's
> **0.1176**.

Two consequences the sprint should carry:

* **The readout question and the ranking question are one problem, with an equals sign.** Any
  future work on "sparse weighted readout" (charter §7C) that does not improve `â` is spending
  effort on the half that is already exact. And the exact optimum is *already* sparse —
  **6.54 of 128 members, with no sparsity penalty imposed** — so sparsity is an output of the
  correct objective, not a design choice to be tuned.
* **Readout 2's ORACLE ceiling is exactly 0 Å**, by rank (`rank(aff{W_x}) = 32.9 ≥ 3n−3` on all
  126, measured residual 0.0000). It is an over-parameterised interpolator: 127 weights against
  ~33 residual dimensions. So "0.2516 Å under an ORACLE objective" is a statement about a
  *regulariser*, not about a class ceiling — and the simplex constraint is precisely the
  regulariser the uniform average enjoys for free. That is the mechanism behind S31-L2(2)
  ("expressivity without an aligned objective is harmful"), derived rather than observed.

---

## 10. ADDENDUM — the `â` arms, and a correction to §7.4's own headline

Requested by the coordinator after §9 was written, because the identity makes any quality
estimate convert into an endpoint for free. Ledger **S31-L19**; code `s31/s31_A_ahat.py`,
`s31/s31_A_setweights.py`; artefacts `s31/results/s31_A_{ahat,setweights,consensus_is_gradient}.json`.
n = 126, CA cloud, fold-clustered, every `â` fitted **leave-fold-out** on the pinned folds.

### 10.1 The arms — ρ first

| arm | ρ global | **ρ in band** (top-24 by true `a`, ORACLE) | amplitude | cloud RMSD | vs PROD75 |
|---|---|---|---|---|---|
| 1 `DIS` z-rank | 0.1176 | **−0.0262** | 0.142 | 3.1846 | +0.1363, 1.20× **WORSE** |
| 2 `CONS` (lane B's `P.mean(1)`) | **0.4278** | **−0.2837** | 0.601 | 3.1826 | +0.1342, 1.71× **WORSE** |
| 3 **`DIS + CONS` LFO** | 0.4281 | −0.2827 | 0.602 | **3.1830** | **+0.1347, 1.71×, 5/5, WORSE** |
| 4 `const` (= MEB) | 0 | — | 0 | 3.1919 | +0.1436, 1.22× WORSE |

`3 − 1`, i.e. **what consensus adds on top of DIS: −0.0016 Å at 0.02× MDE — NOT A RESULT.**
The LFO regression puts essentially all its weight on consensus (β_CONS 0.30–0.36 against
β_DIS 0.004–0.015) and the endpoint does not move.

### 10.2 **Correction to §7.4: ρ_global is not a sufficient statistic, and my crossing price is withdrawn as a target**

Consensus reaches **ρ_global = 0.4278**, double §7.4's crossing price of 0.211, and its
full-amplitude ORACLE-scaled twin lands at **3.1614** where the §7.4 curve predicts ≈ 2.80 at
that ρ. The resolution:

> **Consensus's ρ_global = +0.428 is entirely outlier detection. In band it is −0.284 — the
> wrong sign among the candidates that matter.** `DIS` in band is −0.026, also negative but
> nearly null. §7.4's rungs interpolate toward the truth, so their in-band ρ rises with λ;
> a real feature's need not. **The binding axis is in-band ρ.**

So ρ_global = 0.211 is the crossing price **along the interpolation path only** —
necessary-not-sufficient, and it must not be quoted as a target for a real feature. Amplitude is
not the escape either: the full-amplitude twins move the endpoint by 0.02 Å and stay 2.1× MDE
worse. This is `in-band-is-the-only-ranking-metric` and `consensus-is-outlier-avoidance`
arriving together on the exact objective.

### 10.3 Why consensus can never be the `â` — by identity

The dispersion term's mean-field gradient at uniform weights is `(B·1/D)_x = mean_y ‖W_x−W_y‖²`;
lane B's medoid criterion is `mean_y ‖W_x−W_y‖_RMSD`. Over **all 126 targets**:

    corr( medoid criterion , B @ uniform ) = 0.9667 Pearson (min 0.809), 0.9659 Spearman (min 0.740)

> **The consensus criterion IS the free half of the objective, read at uniform weights.** The QP
> gradient is `â_x − (Bw)_x`, so `â ∝ +consensus` **cancels** the term it was meant to
> complement. The whole family is one scalar κ = (consensus coefficient)/(dispersion
> coefficient): κ < 1 → MEB (3.192), κ → ∞ → the uniform medoid, a single candidate (3.344).
> **Production's uniform top-75 average (3.048) is not on that axis at all.**

That also explains lane B's C2 result — every node-level graph observable is absorbed by
consensus because they are all absorbed by the *same object*, which the readout already
contains.

### 10.4 The 2×2: the deficit is the WEIGHTING RULE, not the set

`{top-75, top-128} × {uniform, MEB}`, with `â ≡ const` so no quality model can confound it.

| | uniform | MEB (derived weights) |
|---|---|---|
| **top-75** | **3.0483** (= production) | 3.1585 |
| top-128 | 3.0532 | 3.1919 |

**SET effect** (top-128 − top-75, at uniform): **+0.0049, 0.14× MDE — NOT A RESULT.**
**WEIGHT effect** (MEB − uniform, on production's own top-75): **+0.1102, 2.68× MDE, 5/5 folds,
49W/77L — WORSE.**

> The entire +0.144 Å is the weighting rule. The MEB support is 6.33 of 75 and loads the extreme
> points; the extreme points really are worse. That is `consensus-is-outlier-avoidance` measured
> *through* the exact objective instead of beside it.

### 10.5 The closing sentence this earns

> The readout's objective is exact and half of it is free; the conversion from any quality
> estimate to an endpoint is a tuning-free convex program; that program's ORACLE ceiling is
> **1.829 Å** against production's **3.048**. **What is missing is a per-candidate quality
> estimate with positive IN-BAND skill.** Every native-free candidate now measured — the shipped
> distogram (−0.026), consensus (−0.284), and every graph observable consensus absorbs — has
> in-band skill that is zero or negative. Not small: **the wrong sign**.
