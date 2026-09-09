# THE QUANTUM COMPONENT — A SPECIFICATION

**Protein-Folding-Algorithm, Sprint 25.** Written for a reader who works on variational
quantum algorithms and will read the source. Every claim below carries a file and line
reference or a measurement; where the project has not measured something it says **NOT
MEASURED** rather than reasoning to a conclusion.

Everything here was independently re-verified in Sprint 25 (`s25/q_verify.py`,
`s25/q_alpha.py`, `s25/q_gibbs.py`, `s25/q_plateau.py` → `s25/results/`), not inherited from
earlier sprints. Two claims in this document's own first draft failed their own registered
falsifiers and were cut; both cuts are marked where they occur.

---

## 0. READ THIS FIRST — WHAT THE COMPONENT DOES AND DOES NOT CONTRIBUTE

**The short version, stated up front because you would find it in ten minutes anyway and the
work is more credible with it stated than hidden.**

**(1) What is genuinely and verifiably quantum-algorithmic here.**

**The deployed selector** is a parameterised real-amplitude circuit — 3 layers of (RY on every
wire, CNOT chain, ring closure) on 7 qubits, 21 parameters — simulated as an **exact dense
statevector**, all 128 amplitudes carried. Its objective is Hamiltonian-derived: `H = diag(E)`
on a candidate-identity register, with a genuine CVaR of the state's own probability
distribution. Its optimiser is Adam on the **exact parameter-shift gradient**. Re-verified this
sprint: the statevector agrees with a dense simulator written independently from the gate list
to **5.6e-17**, and the gradient agrees with central finite differences at cosine
**1.000000000**, relative error **4.6e-10**, for the CVaR term and the entropy term together.

**Separately, and it is not the selector:** the *generation lane* (s19–s21) uses an exact
matrix-product-state ansatz whose bond dimension is `2^layers` **by construction** — no SVD and
no truncation anywhere in the class — verified this sprint against the same dense simulator to
**3.3e-16**, with the Schmidt rank saturating χ exactly. **`core/pipeline.py` never
instantiates `MPSAnsatz`.** The MPS credit belongs to the generation lane and is stated here so
it cannot be carried onto the selector by a fast reader; §1 and §2 keep them apart throughout.

**(2) What that machinery provably cannot do.** The realised CVaR tail's support is always a
**subset of an initial prefix of the energy order**, and equals that prefix exactly when every
state in the prefix carries positive probability. `p_theta` can *delete* a member of the
classical top-*m*; it can never *add* one from outside it. This is a theorem about the
algorithm (§5), not a property of any candidate pool, and it means **the selection is the
classical top-*m* set of whatever energy the component is given.**

**(3) What it demonstrably does not do here.** The one channel the theorem leaves open is
*where the α-mass cuts* — i.e. how spread out the readout weights are. Sprint 25 measured that
channel and it is not a quantum effect:

* Across 18 arms sharing one readout operator, mean RMSD tracks the **entropy of the readout
  weights** at ρ = **−0.7423**, against α at **+0.2700** and T at **−0.0234**. α contributes
  **3.2%** of the variance remaining after H and H².
* Fit that curve on the **nine no-circuit arms only** — arms the circuit is not in — and the
  nine circuit arms sit on it at **+0.0090 Å**, against the fit's own residual sd of 0.0268 Å.
* At the **deployed** temperature a genuine tail constraint (α < 1) is worth
  **−0.0011 Å, 0.01× MDE, NULL.**
* The deployed leave-fold-out table (`core/pipeline.py:118`) runs **α = 1.0 on three of five
  folds**, i.e. **78 of 126 targets, 61.9%, carry no tail constraint at all.** At α = 1 the
  CVaR *is* the mean; there is no tail.
* The whole selector against the shipped argmin is **−0.1405 Å at 0.68× its own MDE —
  below its own MDE and not a result by this project's rule** — although 5 of 5 folds agree in
  sign. Against an exact classical Boltzmann weighting at the same temperature it is
  **−0.0002 Å**.

**(4) What it does do — and this is where an earlier draft of this document was wrong.**

At α = 1 the objective is `mean_p(E) − T·H(p)`, whose **unconstrained minimiser over the
probability simplex is exactly the Gibbs distribution `p*(x) = exp(−E(x)/T)/Z`**. So the
classical `boltz_T` arm is not an analogy — it is *the exact optimum the 21-parameter state is
approximating*, with a closed form. That makes two separate claims measurable, and they came
out differently.

* **The optimiser genuinely optimises, against the control that matters.** Differenced against
  **best-of-200 draws from the same untrained circuit** — never against an initialisation mean,
  per this project's standing rule — the trained state wins at all three temperatures by a
  margin of ≈0.70–0.73 in `F`, and **closes 78–89% of the free-energy gap** between the
  untrained mean and the analytic optimum.
* **It does NOT reach that optimum, and the draft that said it did has been cut.** At the
  **deployed** T = 0.3 the residual is **KL = 0.902 nats (1.30 bits), total variation 0.453** —
  the trained distribution and its own analytic optimum **disagree on 45% of their mass**, and
  the trained state sits **broader** than optimal (5.67 bits against 4.91), so this is not
  under-convergence toward a collapse.
* **And the endpoint cannot tell the difference.** The circuit arm and the exact-Boltzmann arm
  differ by **−0.0302 Å, 0.24× MDE, NULL**.

**(5) THE HEADLINE READING.**

> ## **The readout is insensitive to a distributional difference of nearly half the mass.**

That is the finding. It is the mechanism behind (3): it is *why* only the coarse entropy of the
weights shows up in the endpoint, and why **an exact classical Gibbs state, a trained quantum
state 45% away from it, and a flat average over the top 64 candidates all land within a few
hundredths of an Ångström of each other.** The scale that makes "45%" mean something is the
comparator ladder — a divergence has no units without one:

    T = 0.30, the DEPLOYED temperature (folds 0, 3, 4)
      state                                    F        KL nats   KL bits      TV   H bits
      random theta (untrained)             -1.275443    3.927427   5.66608  0.78052  4.4562
      uniform over 2**n                    -1.455609    3.326874   4.79966  0.70154  7.0000
      point mass at argmin (the collapse)  -1.718572    2.450332   3.53508  0.91374  0.0000
      5 Adam steps                         -1.644194    2.698258   3.89276  0.70994  5.1019
      15 Adam steps                        -1.908150    1.818405   2.62340  0.63160  5.4014
      50 Adam steps (DEPLOYED)             -2.183096    0.901916   1.30119  0.45308  5.6706
      the Gibbs optimum itself             -2.453671    0.000000   0.00000  0.00000  4.9135

**Measured, not inferred — and that distinction is the reason to believe it.** For
`F(p) = E_p[E] − T·H(p)`, the Gibbs state is the exact simplex minimiser and

    F(p) − F(p*)  =  T · KL(p ‖ p*)      **an identity, for any p**

so the free-energy suboptimality **is** the divergence. Both sides come out of the same exact
statevector the deployed driver uses; **no sampling and no estimation enters anywhere**, and
`s25/q_gibbs.py` **asserts the identity at runtime to < 1e-9 on every row**, so it cross-checks
both computations rather than assuming either.

**(6) THE DEEPER POINT, WHICH CORRECTS A CLAIM THIS PROJECT HAD ALREADY MADE.** The standing
explanation was *"the target distribution is classical and cheap to compute, so the circuit
buys nothing"*. **That is wrong, and it was wrong for a reason nobody had checked: the circuit
does not reach the cheap target either.** Cheapness is not the mechanism. **Readout slack is**
— and readout slack would bind identically on a landscape where the Gibbs state were expensive.
That is at once strictly more honest about this system and strictly more favourable to
variational methods in general, and it converts the result from *"our quantum component is
redundant here"* into a transferable statement about a class of pipelines:

> **This readout cannot resolve distributions at this scale. No improvement to the state — by a
> better ansatz, a better optimiser, or better hardware — can reach the endpoint through it.**

**(7) One line.** *The implementation is genuine and correct and the optimiser genuinely
optimises; the accuracy contribution is not established; and the sharpest scientific content is
negative — the downstream readout cannot distinguish the trained state from its own analytic
optimum, from which it differs by 45% of its mass.*

**Two results that cut against us, kept here rather than in a footnote.** At **T = 0.1** the
trained state's free energy (−1.7176) is marginally **worse** than the plain point mass at the
argmin (−1.7186), and its entropy is **0.075 bits**: the entropy term is too weak at that
temperature to hold the state open, and the circuit converges to what is effectively the
collapse. At **T = 1.0** the trained state (KL 0.373) is only modestly better than the
**uniform** distribution (KL 0.458) — the objective is nearly flat there.

**A retraction that belongs here rather than in a footnote.** Earlier project documents state
that "the CVaR tail is worth +0.113 Å by preventing the collapse". The arithmetic is exact and
reproduces to four decimals. It is nonetheless **withdrawn on two independent grounds**: it is
measured at T = 0.1, which is not the deployed temperature, and as a *paired* contrast it is
**0.51× its own minimum detectable effect with a fold-clustered CI spanning zero and a median
of exactly 0.0000 — NULL by this project's own rule.** It was labelled off marginal means where
a paired statistic was required. Recorded in `s25/LEDGER.md` L5.

---

## 1. WHERE THE COMPONENT SITS

    sequence -> distance prior (distogram) -> retrieval pool (k=500)
             -> score filter -> [2^n = 128 candidates]
             -> CVaR-VQE  ->  p_theta  ->  weighted consensus readout  ->  structure

`core/pipeline.py:818-865`, `quantum_stage`. There are **two** distinct quantum paths in this
codebase and conflating them is the first mistake available:

| | **the selector** (this document's subject) | **the generation lane** |
|---|---|---|
| ansatz | `StatevectorCircuit` (`core/quantum.py:885`) | `MPSAnsatz` (`core/quantum.py:459`) |
| register | candidate identity, n = 7, dim = 128 | torsion/basin latent, n up to 96 |
| deployed by | `core/pipeline.quantum_stage` | s19–s21 lanes |
| simulation | exact dense statevector | exact MPS, χ = 2^layers |

`core/pipeline.py` does **not** use `MPSAnsatz`. Statements about χ apply to the generation
lane; statements about the deployed selector apply to `StatevectorCircuit`.

---

## 2. THE ANSATZ, EXACTLY

### 2.1 The deployed selector — `StatevectorCircuit` (`core/quantum.py:885-982`)

    |psi(theta)>  =  prod_{l=1..L} [ U_ent · (X)_{q=0}^{n-1} RY(theta_{l,q}) ]  |0...0>

with `RY(t) = exp(-i t Y / 2)` and the entangler a **CNOT chain plus a ring closure**:

    U_ent = CNOT(n-1, 0) · CNOT(n-2, n-1) · ... · CNOT(1,2) · CNOT(0,1)

Deployed configuration (`core/pipeline.py:186-189`, `Config`):

    n       = 7 qubits          dim = 2^7 = 128 basis states = 128 candidate hypotheses
    layers  = 3                 P = layers * n = 21 parameters
    RY gates = 21               CNOT gates = 21   (6 chain + 1 ring, per layer)
    two-qubit depth = 21        total depth ~ 24
    iters   = 50 Adam steps     restarts = 1      seed = 0

**Why the amplitudes are real.** RY and CNOT are both real matrices and the initial state is
real, so `psi(theta) ∈ R^{2^n}` exactly — the imaginary part is a guaranteed zero, not a small
number. The reachable manifold therefore lies in **SO(2^n)**, not SU(2^n). This matters for any
Lie-algebraic statement about the ansatz and is stated for that reason. `dim so(2^7) = 8128`;
the circuit has 21 parameters.

**Simulation is exact.** All `2^n` amplitudes are carried. Two efficiency choices are pinned by
tests as arithmetic-preserving: the whole CNOT chain plus ring is composed into a single
basis-index permutation once in `__init__` (a permutation is exact), and `probs_batch`
simulates all `2P` shifted parameter vectors in one pass (elementwise-identical arithmetic).
`t_statevector_probs_bit_identical_to_legacy` asserts max |diff| = 0.0. **Verified
independently this sprint** against a dense simulator built from Kronecker products in
`s25/q_verify.py`: max |p − p_dense| = **5.6e-17**.

### 2.2 The generation-lane ansatz — `MPSAnsatz` (`core/quantum.py:459-776`)

    layers x (RY on every wire, CNOT chain)  [+ one trailing RY layer if final_ry]
    no ring closure

**In what precise sense this is an MPS, and why χ is exact rather than truncated.** A CNOT is
`|0><0| ⊗ I + |1><1| ⊗ X`, an exact bond-dimension-2 MPO. Applying the nearest-neighbour chain
therefore **doubles one bond and discards nothing**: after `L` entangling layers every interior
bond is exactly `2^L`, and the two boundary bonds are 1. There is no SVD, no truncation
threshold, and no approximation to trade off. The state is represented as one padded array
`A` of shape `(n, χ, 2, χ)`, and

    amp(x) = e_0^T ( prod_q A[q][:, x_q, :] ) e_0 ,     chi = 2^layers.

Cost is `O(n χ^3)` — **linear in the number of qubits** — so there is no 30-qubit statevector
wall; 64- and 96-qubit registers are routine. Sampling is exact ancestral sampling from the
left-to-right conditional decomposition with right environments, `O(n χ^2)` per shot, so the
`2^n` probability vector is never materialised.

**Verified this sprint** (`s25/q_verify.py`, falsifier: any deviation above 1e-12):

    truncation primitives inside the MPSAnsatz class            NONE
    every "svd" / "truncat" in core/quantum.py                  PROSE ONLY (lines 25, 465, 1591)
    chi by layers                {1: 2, 2: 4, 3: 8, 4: 16}      == 2**layers
    chi at n = 64, layers = 2                                   4   (independent of n)
    built tensor shape at n=8, layers=2                         (8, 4, 2, 4) = (n, chi, 2, chi)
    max |p_MPS - p_dense| over 24 configurations                3.331e-16
    max | <psi|psi> - 1 |                                       6.661e-16

**And χ is saturated, not merely bounded.** The Schmidt spectrum across the middle cut at
n = 6, layers = 2 is `[0.879145, 0.461538, 0.105625, 0.054141, 0, 0, 0, 0]` — **rank exactly
4 = χ**, four non-zero values and four exact zeros. The entangler is doing work: the same
angles through `entangler="none"` (the product-state control) move a probability by up to
0.334.

**Two design points worth knowing.** `final_ry` appends one RY layer with no entangler after
it and is not cosmetic: without it, `b_q` is a prefix XOR of independent Bernoullis, so
`|1 − 2 P(b_q = 1)|` is a running product and must be non-increasing along the wire — 53% of
chain steps in the real problem violate that ordering, and matching a random marginal vector
leaves a residual of 0.21 for the plain chain against 0.000 with one trailing RY layer. And
`entangler="none"` gives a product distribution: it is **the control that says whether the
entanglement is doing anything** (see §7.4).

---

## 3. THE OBJECTIVE, EXACTLY

### 3.1 How the energies are formed

`core/pipeline.py:852`. The top `2^n` candidates by the shipped distogram Bayes-risk score are
taken, and

    E = zrank( score[ top[:2^n] ] )          H = diag(E),   H|i> = E_i |i>

`zrank` is the standardised rank: rank the scores, subtract the mean, divide by the sd. **Two
consequences a reader should draw.**

* `zrank` is **strictly monotone**, so it changes no ordering, no argmin, no level set — and
  therefore, by §5's theorem, **no tail membership**. It changes only the energy *gaps*, which
  is what the entropy term trades against.
* Rank standardisation is not cosmetic. Raw energies on unrelaxed windows are
  outlier-dominated (a moment z-score of raw AMBER puts 98–99.7% of a pool inside |z| < 0.1 and
  measures which candidate has the worst steric clash, not which candidate the potential
  prefers). Conditioning must be monotone; 99th-percentile winsorisation is not enough.

### 3.2 THE HAMILTONIAN BARELY CHANGES BETWEEN TARGETS. **THIS EXPLAINS A LOT.**

`top = argsort(sc)` (`core/pipeline.py:768`), so `sc[top[:128]]` is **already ascending**, and
`_zrank` of an ascending vector returns the standardised ranks 1..128. **E is therefore very
nearly the same vector for every one of the 126 targets.**

It is *not exactly* the same, and the difference matters enough to measure rather than assert:
`rankdata` **averages ties**, and the shipped score does tie. On 8 real targets:

    worst |E - standardised ranks|   4.06e-02   against an E range of 3.4371  =  1.18% of range
    targets where E is EXACTLY the standardised ranks                            0 of 8

> **The spectrum of H is target-independent to within tie-averaging — not exactly constant, but
> within about 1% of its own range. ALL of the per-target information enters through WHICH
> CANDIDATE OCCUPIES WHICH RANK, and essentially none of it through the spectrum of H.**

**Three consequences a reader should draw, and the third is the interesting one.**

1. The deployed selector solves very nearly the **same variational problem** on all 126
   targets. There are effectively **two trained states in the whole deployment** — one per
   `(alpha, T)` cell in `VQE_LFO` — not 126.
2. It is therefore no surprise that the quantum stage is **insensitive to the target**: the
   object it optimises barely knows which protein it is looking at.
3. **It independently explains a result this project measured five separate times and never
   had a mechanism for: deeper ansätze and larger χ order nothing.** A more expressive state
   can only pay if there is target-specific structure in `H` for it to capture. There is almost
   none — `H` is a fixed ladder of standardised ranks. Extra expressivity has nothing to be
   expressive *about*. That is a mechanism, not a restatement, and it was available from four
   lines of source the whole time.

### 3.3 CVaR from the state's probabilities

`p_theta(x) = |<x|psi(theta)>|^2`. For a confidence level `alpha ∈ (0,1]`, with `q` the
`alpha`-quantile of E under `p_theta`,

    CVaR_alpha(theta) = (1/alpha) [ sum_{E(x) < q} p(x) E(x)  +  (alpha - P(E < q)) q ]

i.e. the mean of the lowest-`alpha` mass of the energy distribution, with the boundary state
contributing **fractionally**. Two implementations, both exact and mutually checked:

* `cvar_exact` (`core/quantum.py:985-1005`) — used by the deployed driver. Returns
  `(value, q, dCVaR/dp)`.
* `cvar_from_probs` (`core/quantum.py:336-359`) — the vectorised form that additionally
  returns the **per-state mass the tail used**, which is what makes the value a differentiable
  function of `p` and what the finite-difference tests differentiate.

**The clip in `cvar_from_probs`, line by line, because it is where the theorem lives:**

    order = np.argsort(e, kind="stable")                       # 352   the ENERGY order
    cum   = np.cumsum(p[order])                                # 353   inclusive prefix sum
    take  = np.clip(alpha - (cum - p[order]), 0.0, p[order])   # 354

`take[j] > 0` requires **both** conjuncts of that clip:

* **(i)** `alpha - (cum[j] - p[order][j]) > 0` — the mass strictly *before* position `j` **in
  the energy order** is below `alpha`. `cum - p[order]` is the exclusive prefix sum, which is
  non-decreasing because `p >= 0`, so it crosses `alpha` **exactly once**. Conjunct (i) alone
  defines an initial **prefix** of the energy order.
* **(ii)** `p[order][j] > 0` — a zero-probability state inside that prefix is **punched out**.

So the support is a prefix **with holes**, and every hole is exactly a zero-probability state.
That is §5.

### 3.4 The deployed objective is a free energy, not CVaR alone

`free_energy` (`core/quantum.py:1072-1100`):

    F(theta) = CVaR_alpha(E; p_theta)  -  T * H(p_theta),     H(p) = -sum_x p(x) log p(x)

**Minimising CVaR alone is degenerate for a *selection* task, and this was measured, not
assumed.** For any `alpha` the minimiser concentrates `p` on the lowest-energy basis states, so
a consensus readout collapses to the argmin. Measured: at `alpha = 1, T = 0.1` the state
carries **0.0761 bits** of a possible 7 and the arm is **identical, target by target**, to the
plain argmin selector. The entropy term is what makes the optimum interior.

The deployed `(alpha, T)` come from a **leave-fold-out table** (`core/pipeline.py:113-118`) —
fold *f*'s cell was chosen on the other four folds, so reading it is not tuning on the
instrument:

    VQE_LFO = {0: (1.0, 0.3), 1: (0.25, 0.3), 2: (0.25, 0.3), 3: (1.0, 0.3), 4: (1.0, 0.3)}

**Read that table.** T = 0.3 on all five folds; **α = 1.0 on three of them**, where CVaR is the
full mean and there is no tail constraint at all. See §0(3) and §6.

### 3.5 What is and is not differentiable

* **`CVaR_alpha` as a function of `p`** is concave and **piecewise linear**. By the envelope
  theorem (the quantile `q` is itself the maximiser of `max_t t − (1/alpha) E_p (t − E)_+`, so
  its own dependence on `p` cancels),

      dCVaR/dp(x) = (E(x) − q)/alpha   on the strict tail,   0 elsewhere,

  with the boundary state carrying its partial mass. It is **non-differentiable exactly where
  the quantile crosses a state or where energies tie** — a measure-zero set in `p`, and the
  reason the finite-difference cross-check in §4 is a real check rather than a restatement.
* **`H(p)`** is smooth on the interior; `dH/dp = −(log p + 1)`. The implementation clips
  `p` at 1e-15 inside the log, so the entropy gradient is exact except at states with
  effectively zero amplitude.
* **`p_theta` as a function of `theta`** is analytic and obeys an exact two-term shift rule
  (§4).
* **The composition** `F(theta)` is therefore differentiable almost everywhere, and its
  gradient is exact wherever it exists.
* **What is NOT differentiable at all, and is not treated as if it were:** the *readout*. The
  consensus medoid is an `argmin` over candidates; the emitted structure is a discontinuous
  function of `p`. The optimisation never differentiates through it.

---

## 4. THE GRADIENT

### 4.1 The shift rule on this ansatz

Each parameter enters as `RY(theta) = exp(-i theta Y/2)`. The generator `Y/2` has eigenvalues
`±1/2`, so **every expectation value obeys the two-term parameter-shift rule with shift
`pi/2` and coefficient `1/2`**. The observable here is a **projector**: `p(x) = <psi|Pi_x|psi>`
with `Pi_x = |x><x|`. Therefore *every basis probability individually* satisfies

    dp(x)/dtheta_k  =  [ p(x; theta_k + pi/2)  -  p(x; theta_k - pi/2) ] / 2

exactly — no finite-difference step, no approximation. Chaining with `dCVaR/dp` from §3.5:

    dCVaR/dtheta_k = sum_x  dCVaR/dp(x) · dp(x)/dtheta_k
    dF/dtheta_k    = sum_x  [ dCVaR/dp(x) − T·dH/dp(x) ] · dp(x)/dtheta_k

`grad_cvar_paramshift` (`core/quantum.py:1008-1023`) and `free_energy`
(`core/quantum.py:1072-1100`). Cost: `2P` circuit evaluations per gradient — 42 at the deployed
`P = 21` — batched into one simulation pass.

**Verified this sprint** against central finite differences on the exact objective, which is
independent machinery (no projector algebra, no envelope theorem):

    cos(parameter shift, exact FD)      1.000000000    min over 12 cells, alpha in {0.1,0.25,1}
    rel |g_ps - g_fd| / |g_fd|          4.597e-10      max over 12
    cos(free-energy grad, FD)           1.000000000    min over 4 (alpha, T), incl. the deployed
    rel err, free-energy gradient       4.663e-10      max over 4 -- BOTH terms

`run_cvar_vqe` (`core/quantum.py:1103-1130`) is Adam (β₁=0.9, β₂=0.999, lr=0.15) on that
gradient. **It samples nothing** — verified by source inspection: the only RNG use is drawing
the initial angles `theta ~ N(0, 0.6^2)`.

### 4.2 The recorded defect in the score-function estimator's baseline. **DO NOT REINTRODUCE IT.**

The sampled (device-realisable) estimator uses the score-function form. At the optimal `t` the
envelope theorem leaves

    grad CVaR = E_p[ f(x) grad log p_theta(x) ],     f(x) = -(q - E(x))_+ / alpha

A baseline `b` may be subtracted from `f` **if and only if it is constant in x**, because the
correction term is `b · E_p[grad log p] = 0`. **The historically shipped `qansatz.cvar_gradient`
subtracted the TAIL MEAN from the tail entries and left the rest at zero** — that is
`b(x) = m · 1[x ∈ tail]`, a **function of x**, and the identity does not hold on a
data-dependent subset. It is a **bias**, not extra variance: it does not shrink with shots.

Both `cvar_gradient` (`core/quantum.py:816-863`) and `grad_cvar_score`
(`core/quantum.py:1040-1069`) now default to `baseline="const"`, and `baseline="tail"`
reproduces the defect **verbatim** so the regression test measures its absence rather than
asserting it. Confirmed this sprint: both defaults are `"const"`.

**The magnitude is instrument-dependent and must always be cited with its instrument.** Three
independent measurements, all with **zero sampling noise** (exact expectation values), which is
what proves the disagreement is bias:

| instrument | cos with the exact gradient | ‖g‖ / ‖g_exact‖ |
|---|---|---|
| s9, 10 qubits, 1024 amplitudes enumerated (`core/quantum.py:66`) | **+0.655634** | 0.758 |
| s8, 7 qubits × 3 layers, 36 checks over 12 targets × 3 α (project memory) | **+0.524** | 0.534 |
| **s25 re-verification**, n=7, layers=3, α=0.15, deployed energy shape | **+0.566586** | 0.519 |
| any of the above, CONSTANT baseline | **+1.000000** | 1.000 |
| sampled estimator, constant baseline | +0.994 | ~1 |

**There is no universal constant here and quoting one would be wrong.** What reproduces across
three instruments is the sign, the order of magnitude and the mechanism. A matching constant
would have been a coincidence.

---

## 5. THE SET-EQUALITY THEOREM, AND ITS SCOPE

> **The realised CVaR tail's support is always a SUBSET of an initial prefix of the energy
> order, and equals that prefix exactly when every state in the prefix carries positive
> probability. `p_theta` can delete a member; it can never add one outside the classical
> top-*m*.**

**Proof sketch, from the source (§3.3).** `take[j] > 0` needs both conjuncts of the clip on
line 354. Conjunct (i) is a condition on the exclusive prefix sum in the **energy** order,
which is non-decreasing (because `p >= 0`) and therefore crosses `alpha` exactly once — so (i)
alone defines an **initial prefix**. Conjunct (ii) `p[order][j] > 0` removes zero-probability
states from that prefix. Hence: prefix with holes, every hole a zero-probability state. On the
other path, `tail_indices` (`core/quantum.py:229-272`) takes **no probability vector and no
circuit parameters at all** — `k = ceil(alpha·n)` depends on `(n, alpha)` only and the mask is
a pure function of the energy array.

**Consequences, stated flatly.** The trained state can move exactly two things: **where the
prefix cuts** (a rung `m` on the classical top-*m* ladder, plus which prefix members it deletes
where it has exact zeros), and **the weights inside the prefix**. It can never **add** a
candidate the classical energy order excluded. Membership is bounded above by `argsort(E)[:m]`
as an *identity*, independent of the landscape's shape, degeneracy, multimodality and scale,
and therefore independent of the candidate manifold that produced it. **Changing the pool
cannot break it.** Deletion is not a hidden information channel either: the deleted member is
chosen by an amplitude pattern that is a function of the **same energies** the classical rank
order already uses.

**Scope conditions.** *Subset-hood* is the theorem. *Equality* is the empirical regime: a
trained RY/CNOT state has generic angles and therefore full support, which is why equality
holds in every cell anyone has run. Do not conflate them — an earlier draft of this statement
asserted prefix-hood without conjunct (ii) and was wrong.

**Verification, s25, independent of the s24 test:**

    families rebuilt from scratch WITH EXACT ZEROS (p = 0.0, never +1e-30)
    6 register sizes x 9 energy structures x 8 probability structures x 6 alphas

    n_cells                                              2592
    cells containing an EXACT zero probability           1620      <- the assertion CAN fire
    SUBSET-hood violations   <- THE THEOREM                 0
    holes that were NOT exactly zero-probability            0      <- the mechanism, directly
    prefix-hood violations (the wrong, looser claim)      1424     = 54.9%
    full-support cells                                    972
    full-support cells with EXACT value equality      972 / 972    = 100%

s24's independent run: 3,888 cells, 0 subset violations, 29.9% hole rate. The coordinator's:
17,574 trials, 0 violations, 58.8%. **The hole rates differ because the family mixes differ;
the three assertions agree exactly.**

**AND THE HONEST STATEMENT ABOUT GATE 1.** Sprint 22 reported a set-equality gate passing on
**2016 of 2016 cells** and this was, at the time, read as evidence about the candidate pool.
**It was not. It was the code being read back.** The same applies to the harness test
`t_matched_control_reproduces_the_quantum_arm`, which passes to 1e-9: a reader who sees the
quantum arm exactly matching the classical control should read this section, not celebrate.
A 100% pass rate on a property that is an algebraic identity carries no information about the
problem, and reporting it as if it did is the failure mode this section exists to prevent.

**Two additional guards found this sprint, both in `s24/d_harness.py`, both reported and
neither yet applied** (see `s25/agentQ_FINDINGS.md` §5): `aggregate` computes
`gate_equality_rate` as `r.get("gate_equality", True)`, so a row that never measured equality
is counted as having passed it; and `write` skips the `_COMPLETE` sidecar entirely when
`required_keys=None`, leaving a stale sidecar asserting completeness over new content.

---

## 6. THE HONEST CONTRIBUTION, WITH THE MEASUREMENTS

**BASIS NOTICE.** Two instruments appear and they never share a column.

| instrument | operator | numbers |
|---|---|---|
| **s8** (`s8/integrate_vqe.json`) | consensus-medoid **selection** of one member from a 128-candidate filtered set | 3.4540 / 3.3414 / 3.2835 / 3.3135 |
| **126-target dev** | score-filtered uniform **top-75 coordinate average** | 3.0483 point cloud; **3.2148 built chain** |

The 3.0483 Å arm is a **point cloud**, not a structure: its virtual Cα–Cα bonds are **22.3%
short, worst case 0.649 Å — shorter than a covalent C–C bond** — and `core/bench.py:682` labels
it *"raw average (illegal)"*. **The system's structural number is the built chain, 3.2148 Å**,
independently re-derived this sprint across all 126 targets through a different RMSD
implementation than the one that wrote the records, with zero disagreement. The synthesis arm's
projection gap is **+0.1664 Å** (not +0.156, which is `rmsd_fit`, a different arm). Everything
in this section is on the s8 instrument and is not comparable to either.

All contrasts are paired at n = 126 over 5 pinned folds, with **MDE = 2.8016 × SE per
comparison**; 0.7–1.3× MDE is a Type-M zone and is not a result; below 0.7× is underpowered.
`s25/q_alpha.py` → `s25/results/q_alpha.json`.

### 6.1 The α effect does not survive the temperature

    T = 0.1     a=0.10 - a=1.0   -0.1126  SE 0.0792  0.51x MDE  55W/44L  med +0.0000  NULL
                a=0.25 - a=1.0   -0.1067  SE 0.0644  0.59x MDE  36W/25L  UNDERPOWERED
    T = 0.3     a=0.10 - a=1.0   +0.0279  SE 0.0518  0.19x MDE  50W/50L  NULL
                a=0.25 - a=1.0   -0.0011  SE 0.0477  0.01x MDE  44W/47L  NULL     <- DEPLOYED T
    T = 1.0     a=0.10 - a=1.0   +0.0126  SE 0.0268  0.17x MDE  28W/19L  NULL
                a=0.25 - a=1.0   +0.0039  SE 0.0240  0.06x MDE  18W/24L  NULL

**The effect appears only at T = 0.1, and it changes sign at the other two temperatures.**

### 6.2 One curve: the readout's entropy is the variable

Every one of the 18 arms in the artefact goes through the **same** operator,
`consensus_medoid(D, o, w)` — the member of the same 128-candidate set minimising a `w`-weighted
mean distance. The arms differ **only in `w`**. So the entropy of `w` is exactly comparable
across circuit and non-circuit arms, and is computable in closed form for the non-circuit ones
(`argmin` → 0 bits, `topfrac_f` → log₂ k, `boltz_T` → H(exp(−E/T)), `medoid128` → 7 bits).

    corr(H_readout, mean RMSD)   over the 9 VQE cells      -0.7423
    corr(alpha,     mean RMSD)                             +0.2700
    corr(T,         mean RMSD)                             -0.0234
    alpha's marginal share of the variance left after H, H^2   0.032
    quadratic RMSD ~ H + H^2 over all 18 arms              R^2 = 0.7045

**The construction that makes this airtight:** fit the curve on the **nine no-circuit arms
only** — arms the circuit is not in — then score the nine circuit arms against it.

    mean residual, circuit arms       +0.0090 A   (sd 0.0342, n = 9)
    residual sd of the fit itself      0.0268 A

**The circuit sits on a curve fitted without it, one third of a residual sd above.**

### 6.3 The circuit against its own analytic optimum — at the endpoint AND in distribution

At `alpha = 1` the objective is `mean_p(E) − T·H(p)`, whose unconstrained minimiser over the
simplex is **exactly** `p*(x) = exp(−E(x)/T)/Z`. So `boltz_T` is not an analogy; it is *the
optimum the 21-parameter state is approximating*. **Both the endpoint comparison and the
distributional one are reported, because they say different things and only one of them is a
statement about the ansatz.**

**At the endpoint** (`s25/q_alpha.py`, n=126, s8 instrument):

    circuit - exact Boltzmann @T=0.1   +0.0499  SE 0.0427  0.42x MDE  24W/37L  UNDERPOWERED
    circuit - exact Boltzmann @T=0.3   -0.0302  SE 0.0450  0.24x MDE  43W/44L  NULL
    circuit - exact Boltzmann @T=1.0   +0.0445  SE 0.0511  0.31x MDE  38W/47L  NULL

**In distribution** (`s25/q_gibbs.py`), using the exact identity `F(p) − F(p*) = T·KL(p‖p*)`,
asserted at runtime to < 1e-9 on every row so it is a cross-check rather than an assumption.
The comparator ladder is the point: a divergence has no units without one.

    T = 0.30   (DEPLOYED on folds 0, 3, 4)
      state                                    F        KL nats   KL bits      TV   H bits
      random theta (untrained)             -1.275443    3.927427   5.66608  0.78052  4.4562
      uniform over 2**n                    -1.455609    3.326874   4.79966  0.70154  7.0000
      point mass at argmin (the collapse)  -1.718572    2.450332   3.53508  0.91374  0.0000
      5 Adam steps                         -1.644194    2.698258   3.89276  0.70994  5.1019
      15 Adam steps                        -1.908150    1.818405   2.62340  0.63160  5.4014
      50 Adam steps (DEPLOYED)             -2.183096    0.901916   1.30119  0.45308  5.6706
      the Gibbs optimum itself             -2.453671    0.000000   0.00000  0.00000  4.9135

    KL(trained || optimum), nats:   T=0.1  1.449   T=0.3  0.902   T=1.0  0.373
    total variation:                T=0.1  0.761   T=0.3  0.453   T=1.0  0.351

**And the mandatory control** — best-of-200 from the *untrained* circuit, never an
initialisation mean:

      T      F_trained   F_init_mean   F_init_best    F_gibbs   gap closed   beats best-of-N
    0.10     -1.717571     -0.526843     -0.988733  -1.862495        0.891        True
    0.30     -2.183096     -1.207063     -1.483528  -2.453671        0.783        True
    1.00     -4.936715     -3.587831     -4.214344  -5.309822        0.783        True

> **Three readings, and the order matters.** (i) **The optimiser genuinely trains**: it beats
> best-of-200 from the untrained circuit at every temperature and closes 78–89% of the
> available free-energy gap. (ii) **It does not reach the optimum**: at the deployed
> temperature it is 0.902 nats away and disagrees with it on 45% of its mass, sitting *broader*
> than optimal (5.67 bits against 4.91), not collapsed. (iii) **The endpoint cannot tell**:
> 0.24× MDE. **The readout is insensitive to a distributional difference of nearly half the
> mass**, which is the mechanism for §6.2 and the sharpest thing in this document.

**Two details a careful reader will want.** At T = 0.1 the trained state's free energy
(−1.7176) is marginally *worse* than the plain point mass at the argmin (−1.7186), and its
entropy is 0.075 bits: at that temperature the entropy term is too weak to hold the state open
and the circuit converges to what is effectively the collapse. And at T = 1.0 the trained state
(KL 0.373) is only modestly better than the **uniform** distribution (KL 0.458) — the objective
is nearly flat there.

### 6.4 Against the no-circuit arms

    VQE_LFO - argmin (shipped)    -0.1405  SE 0.0732  0.68x MDE  66W/48L  5/5 folds
                                                                 UNDERPOWERED, NOT A RESULT
    VQE_LFO - Boltzmann T=0.3     -0.0002  SE 0.0457  0.00x MDE  39W/43L  NULL
    VQE_LFO - uniform top-64      +0.0262  SE 0.0279  0.34x MDE  26W/37L  NULL
    VQE_LFO - uniform top-128     -0.0308  SE 0.0590  0.19x MDE  49W/49L  NULL

The headline direction is consistent across all 5 folds, and the effect is still **below its
own minimum detectable effect**. It is reported as underpowered, not as a result.

A concentration check was run because the medians are zero: `VQE_LFO − argmin` has median
−0.0081 against a mean of −0.1405, 12 of 126 targets tie exactly, and 5 targets carry 45.3% of
the effect — but the drop-top-5 statistic sits at the **64.7th percentile of a uniform-effect
null**, so concentration is *suggested* by the median–mean gap and **not established** by the
drop-top test.

### 6.5 The α = 1 folds are a claims problem, not a performance problem

Forcing α < 1 on every fold (an **ORACLE** counterfactual, not a proposal, and not applied) is
worth **−0.0311 Å at 0.27× MDE** (α=0.25) or −0.0021 Å (α=0.10). The leave-fold-out table is
not costing accuracy. **Nothing about the deployed configuration should change on the strength
of this section; what changes is what is claimed about it.**

---

## 7. WHAT A VQE RESEARCHER WILL ACTUALLY ASK

### 7.1 Barren plateaus on this ansatz at this width

**Measured this sprint** (`s25/q_plateau.py` → `s25/results/q_plateau.json`), exact
parameter-shift gradients so there is **no shot noise**: the variance reported is variance over
`theta` alone, which is what a barren-plateau statement is about. `theta ~ N(0, 0.6^2)`, the
deployed initialisation law. `layers = 3` throughout the width sweep.

**Width sweep, `Var_theta[ dF/dtheta_0 ]`.** The deployed cell is the `n = 7` row.

    n     dim     P   |  alpha=1,T=0    alpha=.25,T=0   alpha=.10,T=0   a=1,T=.3    a=.25,T=.3
                      |  (LINEAR COST)  (non-linear)    (non-linear)    (deployed form)
    4      16    12   |  5.358e-02      1.453e-02       4.699e-03       5.491e-02   2.088e-02
    6      64    18   |  1.699e-02      1.140e-02       3.214e-03       1.188e-02   1.559e-02
    7     128    21   |  2.062e-02      1.064e-02       1.922e-03       1.526e-02   1.121e-02
    8     256    24   |  4.358e-03      6.340e-03       2.284e-03       8.389e-03   7.837e-03
    10   1024    30   |  2.713e-03      5.517e-03       3.175e-03       5.731e-03   6.121e-03
    12   4096    36   |  1.218e-03      3.222e-03       2.524e-03       7.741e-03   6.376e-03
    13   8192    39   |  1.072e-03      3.600e-03       2.869e-03       4.802e-03   4.113e-03

    fitted log2 Var per qubit    -0.6492        -0.2522        -0.0472    -0.3105     -0.2429

**No exponential plateau is present at any width measured, in any column.** The steepest decay
is the **linear** cost at −0.649 log₂ per qubit — i.e. Var falls by about 1.6× per qubit added,
against the ≈ −1 per qubit (a factor of 2) a 2-design would give. That is a statement about
*this* shallow ansatz at these widths, not a general one; see the scope conditions below.
Draws per row are 250 (n ≤ 8), 200 (n = 10), 120 (n = 12), 80 (n = 13); the relative SE of a
variance estimate is `sqrt(2/(m−1))`, so ≈ 9% at m = 250 and ≈ 16% at m = 80. Individual rows
should be read with those error bars; the fitted slope spans seven points.

**Depth sweep at the deployed width, n = 7, α = 0.25, T = 0.3, 250 draws each.**

    L      P    Var[dF/dth_0]     E|g|^2 / P
    1      7    2.381e-02         2.491e-02
    2     14    1.658e-02         1.905e-02
    3     21    9.197e-03         1.152e-02      <- DEPLOYED
    4     28    7.999e-03         9.383e-03
    6     42    7.644e-03         8.073e-03
    8     56    7.440e-03         8.046e-03
    12    84    8.213e-03         7.713e-03

**Depth costs a factor of ~3 in gradient variance and then saturates.** Between L = 4 and
L = 12 the variance is flat to within its own sampling error — there is no exponential
attenuation with depth in the measured range, and the deployed L = 3 already sits close to the
plateau. (This row's n = 7 entry, 9.20e-3, uses a different seed from the width sweep's
1.12e-2; the two are about two sampling-SEs apart and are consistent.)

**Scope conditions, so this is not over-read.**

* The **deployed register is n = 7**. Every larger `n` here is an extrapolation instrument, not
  a deployed configuration.
* `P = 3n` at `layers = 3`, against `dim so(2^n) = 2^{n−1}(2^n − 1)` — at n = 7 that is **21
  parameters against 8128**. **This circuit is nowhere near a 2-design at any width measured**,
  so an observed decay rate is a property of *this shallow, structured ansatz* and must not be
  quoted as a 2-design result. The classic exponential-in-`n` plateau is a statement about
  Haar-random or 2-design circuits; this is not one.
* The ansatz is real-amplitude, so the relevant algebra is a subalgebra of **so(2^n)**, not
  su(2^n). **NOT MEASURED:** the dynamical Lie algebra actually generated by
  `{Y_q}` and the CNOT-chain conjugations at fixed depth, and therefore whether this ansatz is
  controllable or lies in a polynomially-sized DLA. This project has not computed it and this
  document will not guess it. It is a well-posed question and the honest answer is that nobody
  here has done the work.
* **NOT MEASURED:** gradient variance at the deployed *initialisation* versus along the
  optimisation trajectory. The numbers above are at random initialisations.

### 7.2 Does the CVaR non-linearity change the gradient-variance picture?

**A methodological point first, because it is the part of this section worth reusing.** This
question normally requires comparing a CVaR objective against a *separately chosen* linear
objective, which confounds the non-linearity with the choice of observable, its spectrum, its
norm and its optimiser. **Here no such choice is needed.** At `alpha = 1` the CVaR is not
*analogous* to the linear cost — it **reduces to it identically**:

    CVaR_1(theta) = E_{p_theta}[E] = <psi(theta)| diag(E) |psi(theta)>

So setting `alpha = 1` gives the linear-cost control **for free, from the same expression**:
identical circuit, identical Hamiltonian, identical parameter distribution, identical
estimator, identical code path — with **only the non-linearity switched on or off**. Every
confound is held fixed by construction rather than by matching. The question then reduces to a
ratio at matched `n`, and is answered by measurement rather than by argument.

    ratio  Var[ grad CVaR_alpha ] / Var[ grad MEAN ],  identical circuit, identical theta law
      n      alpha = 0.25     alpha = 0.10
      4         0.2712           0.0877
      6         0.6707           0.1892
      7         0.5159           0.0932        <- DEPLOYED WIDTH
      8         1.4549           0.5241
     10         2.0336           1.1704
     12         2.6447           2.0715
     13         3.3585           2.6765

**Yes, and in two opposite directions depending on what you ask.**

* **In magnitude, at the deployed width, the non-linearity SHRINKS the gradient.** At n = 7 the
  CVaR gradient's variance is **0.52×** the mean's at α = 0.25 and **0.093×** at α = 0.10.
  This is unsurprising and is not a plateau: a tail objective's gradient weight
  `(q − E)_+/alpha` is supported on an `alpha` fraction of the states, so it is a smaller
  quantity.
* **In SCALING, the non-linearity FLATTENS the decay, and this is the interesting half.** The
  linear cost decays at −0.649 log₂ per qubit; α = 0.25 at −0.252; α = 0.10 at **−0.047, i.e.
  essentially flat over 4 → 13 qubits.** The ratio therefore crosses 1 near n ≈ 8 and reaches
  **3.36× (α=0.25) and 2.68× (α=0.10) by n = 13**. On this ansatz, the CVaR objective's
  gradient is better conditioned at width than the linear one it reduces to.
* The two deployed-form columns (`T = 0.3`) sit between: −0.311 at α = 1 and −0.243 at
  α = 0.25. **The entropy term also flattens the decay** — expected, since `−T·H(p)` grows in
  importance as `p` spreads over a larger register.

**A reading, offered as a reading and not as a result.** The mean cost averages the gradient
weight over all `2^n` states, so its per-parameter signal dilutes with the register. The CVaR
weight is supported on a fixed *fraction* rather than a fixed *number*, and the tail's spread
in E does not shrink with `n` under the rank-standardised spectrum of §3.1 — so the dilution
that hits the mean does not hit it the same way. **This project has not tested that mechanism
and it should not be quoted as established.**

**What this does and does not say.** It is a measurement on this ansatz, at these widths, at
this initialisation law, on this Hamiltonian spectrum (standardised ranks, §3.1). It is **not**
a general theorem about CVaR objectives, and no such theorem is claimed. **NOT MEASURED:**
whether the same ordering holds for a different spectral shape, or at depths where the circuit
approaches a 2-design.

### 7.3 Shot noise, and what would change on hardware

**Everything reported in this document is computed exactly and is therefore
shot-count-independent.** Specifically: `p_theta` comes from an exact simulation of all `2^n`
amplitudes; the CVaR and its tail are read off `cvar_exact`/`cvar_from_probs` rather than
estimated from draws; the gradient is the exact parameter-shift gradient. There is no shot
noise anywhere on the deployed selector path, which is what makes the gradient audit in §4 a
verification rather than a noise comparison.

**On hardware, at the deployed size:**

    qubits                    7
    RY gates                 21          CNOT gates  21     (6 chain + 1 ring per layer)
    two-qubit depth          21          total depth ~24
    parameters               21
    circuits per gradient    2P = 42
    circuits per target      50 Adam steps x 42 = 2100 circuit settings, x shots each

**And one thing that would change qualitatively, which is the honest answer to "why not run it
on hardware".** The deployed gradient is assembled from `dCVaR/dp(x)` for *every basis state*,
so it needs the **full 128-outcome distribution** resolved to useful precision. That is not a
Pauli-expectation estimation problem with `O(1/eps^2)` cost on a few observables; it is
distribution reconstruction, and the shot cost scales with the size of the support. Two
routes exist for a device and both are in this codebase:

* **Sampled CVaR with a derivative-free optimiser.** `cvar_from_samples` estimates the CVaR
  from the empirical `alpha`-quantile of measured energies — the standard CVaR-VQE construction
  — and `run_global_cvar_vqe` (`core/quantum.py:1448`) drives it with SPSA, which costs
  2 evaluations per iteration regardless of dimension and whose convergence theory is *for*
  noisy evaluations. The plug-in quantile carries an `O(1/(alpha·S))` bias.
* **The score-function gradient.** `cvar_gradient` (`core/quantum.py:816`) with
  `baseline="const"` and an ansatz whose `grad log p` is available analytically (the
  `OneLayerAnsatz` and `MPSAnsatz` families both provide it). **`baseline="tail"` must never be
  used** — §4.2.

Neither is what produced the numbers in §6, and this document does not claim hardware results.

### 7.4 Is the entanglement structure doing anything?

**Two separate questions, and the answers differ.**

* *Does the entangler change the state?* **Yes, provably and measurably.** χ saturates at
  `2^layers` (Schmidt rank exactly 4 at layers = 2, §2.2), and the product-state control
  (`entangler="none"`) at identical angles moves probabilities by up to 0.334.
* *Does it change the answer?* **NOT MEASURED, and specifically not refuted.** The standing
  project result, from deleting the CNOTs and changing nothing else, is
  **−0.013 Å [−0.095, +0.077]** — an interval containing zero and containing effects in both
  directions. The project also has a standing result that **χ orders nothing**: increasing the
  bond dimension does not reorder candidates. This lane did not re-derive either and does not
  claim more than they say.

Given §6.3, there is a coherent reading available and it is worth stating as a *reading*, not
a result: the target distribution `exp(−E/T)/Z` over a **diagonal** H is a product-like object
with no correlations to represent, so there may be nothing for the entanglement to do *on this
problem*. That is consistent with everything measured and is not itself measured.

### 7.5 The questions this document cannot answer

* **NOT MEASURED:** whether the circuit contributes on a landscape where the classical Gibbs
  optimum is *not* cheap. Every result in §6 rests on a 128-dimensional **diagonal** H whose
  exact Gibbs state is one `exp` call. That is a property of this deployment, not of CVaR-VQE.
  **But §0(6) and §6.3 have narrowed why this matters, and the narrowing goes the other way
  from what you would expect:** the circuit does *not* reach that cheap optimum either, so
  "the target is classical" is **not** the mechanism. The binding constraint is **readout
  slack** — a trained state and its own optimum, 45% of their mass apart, land 0.24× MDE apart
  at the endpoint — **and readout slack would bind identically on a landscape where the Gibbs
  state were expensive.** So this particular caveat is weaker than it looks: the limitation
  identified here is not special to a cheap-target deployment.
* **NOT MEASURED:** the dynamical Lie algebra of this ansatz (§7.1).
* **NOT MEASURED:** anything on hardware.
* **A standing negative worth knowing:** at 8192 objective evaluations, no sampler — quantum or
  classical — beat the zero-evaluation shipped retrieval pool on realised RMSD. In this project
  **search is not the binding constraint; discrimination is.** Concentrating a distribution is
  the wrong move when discrimination binds, which is the same fact §3.4's entropy term is
  working around from the other side.

---

## 8. HOW TO CHECK ALL OF THIS YOURSELF

    python s25/q_verify.py    # MPS exactness and chi, parameter shift, the theorem, the audit
    python s25/q_alpha.py     # the alpha/temperature/entropy analysis, n=126
    python s25/q_gibbs.py     # KL and TV to the analytic optimum, plus the training control
    python s25/q_plateau.py   # gradient variance vs width and depth
    python -m pytest tests/test_quantum.py

`tests/test_quantum.py` uses the five pre-consolidation modules (`qansatz.py`, `vqe.py`,
`foldvqe.py`, `objective.py`, `hamiltonian.py`, all still on disk) as an **equivalence oracle**
for every routine that was rewritten — they are the reference, not the implementation. The
single most important test in that file asserts the §4.2 defect is gone **by measuring it**
rather than by reading the source.

Artefacts, all provenance-stamped via `s24/stats_lib.save_atomic`:
`s25/results/{q_verify,q_alpha,q_gibbs,q_plateau}.json`.
Full working notes and the five-defect methodological list: `s25/agentQ_FINDINGS.md`.
Pre-registration and fork lists: `s25/PREREG_Q.md`.
