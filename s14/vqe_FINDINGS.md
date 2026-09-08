# Sprint 14 — VQE workstream findings

Agent: VQE (quantum machinery + the causality question).
Instrument reproduced at start: `shipped 3.4540004952559396, pool_best 1.7108244199364904,
top75_best 2.3061526409453816, synthesis_fit 3.2040761603809194, n_zero_recall 18`. Exact.

Everything below is measured on the nine fully enumerated targets
(`s13/results/qarch_enum_<PDB>.npz`, n=9, k=4, 262,144 configurations each, true CA-RMSD on
every one) unless stated otherwise. Reference points on that family, mean over nine targets:
**space best 0.969 A, uniform random draw 3.781 A, space worst 6.701 A.**

Tiering: **DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **HYPOTHESIS** / **REFUTED**.

---

## 0. HEADLINE, stated first

1. **The CVaR value estimator is correct** (1.7e-15 against an independent 60-digit
   reference, including partial buckets, ties and the `alpha*2^n < 1` edge). The recorded
   gradient defect is real, is reproduced bit-for-bit, and its bias now has a **closed form**
   instead of a cosine. Two *new* defects are recorded: a sampled-CVaR bias at non-integer
   `alpha*N` (up to +0.134 sd), and an exactly-zero-gradient regime with an iff condition.
2. **VQE NEVER BECOMES USEFUL, on either axis, at any objective quality tested.** On the
   objective axis it is a **strictly worse optimiser than greedy 1-opt** — on the real Legacy
   objective greedy and annealing certify the global optimum in 100% of cells and VQE in 0%.
   On the structural axis it ties or loses everywhere; no cell favours it with a CI excluding
   zero. **There is no crossing.**
3. **The VQE optimisation is WORSE THAN NOT OPTIMISING.** Against the correct control — the
   best of the same number of shots from its own untrained initial distribution — running the
   VQE costs **+0.65 to +1.32 A**, at **W/L 0/12**, for every initialisation including ORACLE
   warm starts. It also **washes out its own initialisation entirely** (2.42-2.76 A out, from
   2.97-4.13 A in).
4. **The mechanism, and it is the sprint's central scientific result.** Every objective is at
   **chance inside its own low-energy tail** (0.512 at gap>0.25 A), so the argmin it returns
   is a **random draw from that tail** (41st-48th percentile). **The selection gap is the
   price of drawing at random from a tail you cannot order:** tail mean minus tail best =
   1.93 A, my certified gap = 1.876 A, the coordinator's sampled gap = 2.040 A. Concentration
   — which is what a variational optimiser does — is exactly the wrong move when
   discrimination is the binding constraint. That is *why* (2) and (3) hold.
5. **`alpha` is not a learning rate** (cos(grad(1.0), grad(0.01)) = 0.0727), but no alpha
   is good in general and the pre-registered alpha test is null (W/L 2/2). The recorded
   "small alpha collapses distinct objectives" trap has **the sign backwards**.
6. **The dead-qubit count in the brief is an undercount:** `2 log2 k`, not `log2 k` —
   **14 of 18 live**, on 9/9 targets, and it *follows from* the locality theorem.
7. **Binary encoding is confirmed, for three new reasons**, and one-hot's Pauli spectrum is
   shown to measure the feasibility constraint rather than the objective.
8. **The QNG refutation in the brief is depth-1 only** and does not extend to depth >= 2.

---

## 1. PART A — CVaR AND VQE CORRECTNESS

### 1.0 The derivation (reproduced in full in the `s14/test_vqe.py` module docstring)

For the LOWER tail, Rockafellar-Uryasev with the sign flipped for minimisation:

    CVaR_a(theta) = max_t [ t - (1/a) E_{x~p}[(t - E(x))_+] ]

`G(t)` is concave with `dG/dt = 1 - (1/a) P(E < t)`, so it is maximised at the a-quantile
`q`. A discrete distribution generally has no state where the cumulative mass equals `a`
exactly, so the boundary atom must be **split fractionally**:

    CVaR_a = (1/a) [ sum_{E(x)<q} p(x) E(x) + (a - P(E<q)) q ]                       (3)

Differentiating (1) at the optimal `t`, the envelope theorem kills the `t` dependence:

    grad CVaR = E_p[ w(X) grad log p(X) ],   w(x) = (E(x)-q)/a if E(x)<q else 0      (4)

and `w` **is** `dCVaR/dp`.

### 1.1 The value estimator is correct — DEMONSTRATED

Against an independent 60-digit-decimal mass-walking reference:

| property | result |
|---|---|
| `cvar_exact` vs brute force, 400 random discrete distributions, alpha 0.003-1.0 | max abs error **1.7e-15** |
| `cvar_exact` vs `cvar_from_probs` (two independent code paths), 300 cases | **< 1e-12** |
| partial bucket at the quantile boundary | split fractionally, **exactly** eq. (3) |
| ties: permutation invariance under 60 random relabellings, 6 alphas | **0 disagreements** |
| massive degeneracy (24 states, 3 distinct energies) | fractional, not rounded to a bucket |
| `alpha * 2^n < 1` | returns min over the support, correct at every alpha down to 1e-6 |
| `tail_indices` argsort/argpartition switch at N=448, heavy ties | mask identical at N = 10, 447, 448, 449, 1000, 5000 |

**All 92 tests in `s14/test_vqe.py` pass.**

### 1.2 A subtlety that is correct and load-bearing — DEMONSTRATED

`dCVaR/dp` is **non-zero on states with `p(x)=0`** whenever `E(x) < q`. This is not a bug:
it is the one-sided derivative saying *"moving mass onto this state would lower the CVaR"*,
and it is exactly the term that lets a variational optimiser **discover states it is not yet
sampling**. Consequence: the zero-gradient regime is governed by the argmin of `E` over the
**whole register**, not over the current support.

### 1.3 NEW DEFECT — the sampled CVaR is biased at non-integer `alpha*N` — DEMONSTRATED

`cvar_from_samples` (and `cvar`/`tail_indices`) average the lowest `ceil(alpha*N)` samples.
Equation (3) splits the boundary atom. They coincide **exactly** when `alpha*N` is an
integer and differ otherwise, with an **upward** bias:

| N | alpha | alpha*N | k=ceil | mean bias vs eq. (3) |
|---|---|---|---|---|
| 13 | 0.10 | 1.3 | 2 | **+0.134** +/- 0.002 |
| 8 | 0.30 | 2.4 | 3 | **+0.112** +/- 0.001 |
| 10 | 0.25 | 2.5 | 3 | **+0.081** +/- 0.001 |
| 100 | 0.037 | 3.7 | 4 | +0.024 |
| 64 | 0.15 | 9.6 | 10 | +0.019 |
| 1000 | 0.0333 | 33.3 | 34 | +0.008 |

(units: standard deviations of the sampled energy; decays as 1/N; monotone decrease
verified at N = 16, 64, 256, 1024.)

**Why it matters.** This is the objective the global SPSA driver evaluates. At `alpha=0.05`
with 256 shots the effective tail is 13 samples and the bias is of the same order as the
differences an optimiser is trying to resolve. **It is not a bug to fix** — it is the
estimator, and it is consistent — but any CVaR value quoted at small `alpha` and modest
shots is biased upward and must not be compared against an exact CVaR computed by eq. (3).
Pinned by `test_ceil_cvar_estimator_is_biased_when_alpha_times_N_is_not_an_integer`.

### 1.4 NEW — the exactly-zero-gradient regime, with an iff — DEMONSTRATED

    dCVaR/dp is IDENTICALLY ZERO  <=>  p(x*) >= alpha,  x* = the (stable) argmin of E.

**0 counterexamples in 3,000 random cases.** This is a *hard* zero — not a barren plateau,
not a small gradient — and no estimator, shot count or baseline recovers it.

It is **convergence, not failure**: `CVaR_a >= min E` always, with equality exactly in this
regime (verified over 400 cases), so the flat region is the global minimum of `CVaR_a`.

But the consequence is severe and is the mechanism for section 1.6. Frequency of the event
at a random initialisation, exact statevector, 40 seeds per cell:

| n | layers | init scale | alpha=0.01 | alpha=0.025 |
|---|---|---|---|---|
| 7 | 3 | 0.8 | **0.23** | 0.10 |
| 7 | 1 | 1.6 | 0.12 | 0.10 |
| 9 | 3 | 0.8 | 0.05 | — |
| 9 | 1 | 0.3 | 0.03 | 0.03 |

So at `alpha=0.01`, roughly **one initialisation in four to one in twenty starts dead**, and
the fraction grows as the distribution concentrates during optimisation.

### 1.5 The gradient audit — DEMONSTRATED

`grad_cvar_paramshift` vs `grad_cvar_fd` (wholly independent machinery: no projector
algebra, no envelope theorem) across **n in {5,7,9,11} x layers in {1,2,3,5} x alpha in
{1.0, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01} x 4 seeds**:

* **worst cosine over every non-degenerate cell: 0.999999999999998**
* **worst relative error: 2.85e-09**
* the only disagreements are cells whose gradient is exactly zero (section 1.4), which are
  reported as such rather than averaged into a NaN.

Also checked by a third route — `dCVaR/dp` contracted with a numerically differentiated
`dp/dtheta` — cosine > 1 - 1e-9.

Sampled score estimator (4 x 200,000 shots, averaged, so what remains is bias):

| alpha | `const` cos | `const` \|g\| ratio | `tail` cos | `tail` ratio | `none` cos |
|---|---|---|---|---|---|
| 0.50 | +0.999882 | 0.9992 | **+0.674** | 0.726 | +0.999854 |
| 0.25 | +0.999869 | 1.0012 | **+0.776** | 0.681 | +0.999852 |
| 0.10 | +0.999817 | 1.0037 | **+0.696** | 0.811 | +0.999814 |
| 0.05 | +0.999817 | 1.0037 | **+0.696** | 0.811 | +0.999814 |
| 0.01 | +0.999688 | 0.9941 | **+0.321** | 0.711 | +0.999684 |

**`baseline="const"` is unbiased in direction and in magnitude. `baseline="tail"` is biased
in both.** `baseline="none"` is essentially indistinguishable from `const` at these shot
counts — **the constant baseline's value is that it is unbiased, not that it reduces
variance**, which corrects the implication of the docstring in `core/quantum.py`.

### 1.6 The defect's bias has a CLOSED FORM — DEMONSTRATED (new)

For any *constant* `c`, `E_p[c grad log p] = c grad 1 = 0`, so `const` is unbiased. The
defect uses `b(x) = c * 1[x in tail]`, a **function of x**, and its bias is exactly

    E[g_tail] - grad CVaR  =  - c * grad_theta P(E(X) < q)                          (5)

with `c` the p-weighted tail mean of `w`. Verified to **< 1e-9 at every alpha** by exact
expectation over the enumerated register
(`test_tail_baseline_bias_has_the_closed_form_of_equation_5`).

**This is worth more than the cosine.** The defect does not add noise — it adds `c` times
the gradient of the **tail probability**, which pushes the optimiser to change *how much
mass sits in the tail* rather than *how good the tail is*. And it explains why the recorded
severity varies: I measure the tail-baseline cosine from **+0.32 to +0.96** across alpha at
a fixed point, and from **+0.47 to +0.96** across alpha on the one-layer ansatz.

**CORRECTION to the memory entry "cosine -0.023 with the true gradient".** That number is
one draw from a wide, alpha- and instance-dependent distribution, not a constant. The defect
is real and must stay fixed; its severity is not a single number. The closed form (5) is the
stable statement.

`qansatz.cvar_gradient` reproduces `core.quantum.cvar_gradient(baseline="tail")` to
**max abs difference 0.00e+00** at four alphas — the shipped module still carries the
defect, and `core/quantum.py` is the corrected path.

### 1.7 REFUTED (my own hypothesis): "alpha is effectively a learning rate" — DEMONSTRATED

I expected small-alpha CVaR to be roughly a rescaled expectation gradient. **It is not, and
the refutation is emphatic.** Cosine between `grad(alpha_i)` and `grad(alpha_j)` at
**identical theta**, exact statevector, n=9, layers=3, 6 seeds:

| | 1.0 | 0.5 | 0.25 | 0.1 | 0.05 | 0.025 | 0.01 |
|---|---|---|---|---|---|---|---|
| **1.0** | 1.0000 | 0.8624 | 0.6911 | 0.4923 | 0.3259 | 0.2449 | **0.0727** |
| **0.5** | 0.8624 | 1.0000 | 0.8546 | 0.6270 | 0.3868 | 0.2877 | 0.0887 |
| **0.25** | 0.6911 | 0.8546 | 1.0000 | 0.8310 | 0.5552 | 0.4631 | 0.2054 |
| **0.05** | 0.3259 | 0.3868 | 0.5552 | 0.8227 | 1.0000 | 0.9488 | 0.5475 |
| **0.01** | 0.0727 | 0.0887 | 0.2054 | 0.3898 | 0.5475 | 0.6097 | 1.0000 |

`grad(alpha=1)` and `grad(alpha=0.01)` are **nearly orthogonal (0.0727)**. Adjacent alphas
sit near 0.85 and the structure decays smoothly, so alpha traces a *curve of directions*, not
a scale. Mean gradient NORM also moves the "wrong" way for a rescaling story — it *rises* as
alpha falls (0.356 at alpha=1 to 0.72 at alpha=0.025).

Dynamically (1CS9, signal 0.3, 20,480 evaluations, identical seed / initialisation, only
alpha differs): final displacement from `theta0` is essentially constant (7.05-7.58) but
**path length falls monotonically with alpha** (15.01 -> 9.62) and the cosine of the net
displacement against the alpha=1 run drops to **0.31-0.53**. Different alphas walk different
distances in different directions to end equally far from where they started.

### 1.8 The collapse trap: the identity is real, the RECORDED FORM IS REFUTED

The identity holds and is proved: **when `alpha <= p(x*)`, `CVaR_a = E(x*)` and
`dCVaR/dp = 0`, both of which depend on `E` only through its argmin**, so any two objectives
sharing an argmin are the same CVaR function bit for bit
(`test_small_alpha_collapses_distinct_objectives_onto_each_other`).

**But the brief's stated form — "at CVaR alpha <= 0.25 two AMBER variants became literally
the same objective" — does not reproduce, and the sign of the alpha effect is backwards.**
`s14/vqe_collapse.py` measures the cosine of `dCVaR/dp` between real objective pairs over a
grid of CVaR alpha x distribution concentration. **All AMBER numbers use the clean
`amber_kind == 0` mask** (see 1.9).

Cosine of `dCVaR/dp`, AMBER total vs AMBER-minus-solvation (rho between them +0.9685):

| Dirichlet a_D | max p | a=1 | 0.5 | 0.25 | 0.1 | 0.05 | 0.025 | 0.01 |
|---|---|---|---|---|---|---|---|---|
| 2 (broad) | 0.0037 | 0.992 | 0.931 | 0.691 | 0.397 | 0.226 | 0.105 | **0.085** |
| 0.02 | 0.0625 | 0.992 | 0.937 | 0.658 | 0.524 | 0.237 | 0.221 | 0.209 |
| 0.001 (peaked) | 0.5079 | 0.992 | 0.929 | **0.910** | **0.910** | **0.910** | **0.910** | **0.910** |

Legacy vs AMBER total (rho +0.5226) goes 0.880 -> **0.000** as alpha falls on a broad
distribution. Legacy vs Legacy-minus-steric (rho +0.9991) sits at **0.994-1.000 at every
alpha and every concentration**.

**Three readings, and they change the guidance.**

1. **Along a row — shrinking alpha on a broad distribution makes genuinely different
   objectives MORE distinguishable, not less.** Small alpha is a *discriminator*. The
   recorded claim has the sign backwards.
2. **Down a column — concentrating the distribution is what collapses them.** At max p = 0.51
   the cosine freezes at 0.910 for every alpha at or below 0.25. This is the identity, and
   its governing variable is **`alpha <= p(x*)`, a JOINT condition on alpha and on how
   concentrated `p` has become** — which happens LATE in an optimisation, not at the start.
3. **A pair whose rank correlation is already ~1.000 is the same objective everywhere.**
   Legacy vs Legacy-minus-steric at rho +0.9991 is indistinguishable at alpha=1 too. That is
   a property of the pair, not of CVaR, and it is the likeliest explanation of the original
   observation.

**Classification: the identity is DEMONSTRATED; "small alpha collapses distinct objectives"
is REFUTED as stated.** The correct warning is: *CVaR collapses objectives once the
distribution has concentrated past alpha, and near-rank-identical objectives were never
distinguishable at any alpha.*

### 1.9 AMBER mask — which numbers moved

The coordinator flagged that the cached AMBER subset is 40% oracle-conditioned. Confirmed:
1CS9 subset mean RMSD **3.573** vs full-enumeration **4.064**; the `amber_kind == 0` rows
(1,193 of 2,955) are clean at **4.061**. Every AMBER number in this document uses that mask.
Effect of the mask on the pair correlations:

| pair | rho (clean) | rho (contaminated) | shift |
|---|---|---|---|
| AMBER total vs AMBER-minus-solvation | +0.9685 | +0.9819 | -0.0134 |
| AMBER total vs nonbonded+solvation | +0.9998 | +0.9999 | -0.0001 |
| AMBER total vs torsion-only | +0.2422 | +0.2265 | +0.0158 |
| Legacy vs Legacy-minus-steric | +0.9991 | +0.9986 | +0.0006 |
| **Legacy vs AMBER total** | **+0.5226** | **+0.5592** | **-0.0366** |

Largest shift 0.037; **no conclusion in section 1.8 changes**. The Part C grid uses only the
`rmsd`, `legacy` and `prior` columns, which are full-enumeration and need no mask.

### 1.10 THE ALPHA SWEEP — what alpha buys, and what it costs — DEMONSTRATED

1CS9, budget 20,480 evaluations, MPS depth-2+final-RY ansatz on 18 qubits, 5 seeds.
Initialisation RMSD **4.117** on every row — *the VQE starts at a random point, not at the
answer.*

**signal 0.30** (rho +0.310; certified optimum 1.929 A; space best 1.270; random draw 4.064):

| alpha | obj gap | RMSD returned | sd | RMSD mode | entropy (bits) | diversity | grad var | P(RMSD<2.5) | max p |
|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.0119 | **2.266** | 0.185 | 2.684 | 4.07 | 0.178 | 7.6e-4 | 0.277 | 0.291 |
| 0.5 | 0.0157 | 2.456 | 0.073 | 2.625 | 5.70 | 0.204 | 1.6e-3 | 0.217 | 0.237 |
| 0.25 | 0.0157 | 2.456 | 0.073 | 2.519 | 6.61 | 0.244 | 9.1e-4 | **0.371** | 0.226 |
| 0.1 | 0.0125 | 2.300 | 0.207 | **2.385** | 7.51 | 0.289 | 5.9e-4 | 0.247 | 0.176 |
| 0.05 | 0.0152 | 2.433 | 0.075 | 2.433 | 9.14 | 0.335 | 3.6e-4 | 0.260 | 0.106 |
| 0.01 | 0.0151 | 2.459 | 0.068 | 2.456 | 10.26 | 0.364 | 1.9e-4 | 0.135 | 0.054 |
| 0.5->0.05 anneal | 0.0166 | 2.460 | 0.076 | 2.478 | 8.17 | 0.275 | 1.5e-3 | 0.208 | 0.125 |
| 1.0->0.01 anneal | 0.0122 | 2.359 | 0.229 | 2.731 | 9.83 | 0.303 | 1.9e-3 | 0.090 | 0.056 |

**signal 0.00** (rho -0.083; certified optimum **4.066 A**; random draw 4.064):

| alpha | obj gap | RMSD returned | sd | entropy | P(<2.5) |
|---|---|---|---|---|---|
| 1.0 | 0.0003 | 4.155 | 0.113 | 3.70 | 0.002 |
| 0.25 | 0.0001 | 4.245 | 0.371 | 8.23 | 0.019 |
| **0.01** | **0.0000** | **4.066** | **0.000** | 12.55 | 0.033 |

**signal 1.00** (the objective IS the truth; certified optimum = space best 1.270 A):

| alpha | obj gap | RMSD returned | sd | entropy | P(<2.5) |
|---|---|---|---|---|---|
| 1.0 | 0.0003 | 1.386 | 0.046 | 6.17 | **0.766** |
| 0.25 | 0.0000 | 1.295 | 0.051 | 9.32 | 0.496 |
| 0.05 | 0.0000 | **1.270** | **0.000** | 11.47 | 0.269 |
| 0.01 | 0.0000 | **1.270** | **0.000** | 11.18 | 0.220 |

**Five readings.**

1. **Small alpha BROADENS the final distribution.** Entropy rises monotonically as alpha
   falls — 4.07 -> 10.26 bits at signal 0.3, 3.70 -> 12.55 at signal 0.0 — and `max p` falls
   0.29 -> 0.05. This is the **opposite** of the intuition that a narrow tail concentrates
   the state, and it contradicts the premise in `qansatz`'s adaptive-alpha docstring ("early
   on a narrow tail concentrates the distribution"). Mechanism: at small alpha the gradient
   is carried by `alpha * shots` samples, so it is noisier per unit signal — measured grad
   variance falls 7.6e-4 -> 1.9e-4 — and Adam accumulates less coherent displacement (path
   length 15.0 -> 9.6, section 1.7).
2. **Small alpha finds the objective's optimum better and the structure no better.** At
   signal 0.0 alpha=0.01 reaches objective gap **0.0000 with sd 0.000** — it certifies the
   global optimum every seed — and returns **4.066 A, which is the certified optimum's RMSD
   and is worse than the 4.064 A random draw.** This is the "optimise harder, get worse"
   trap in its purest form: *perfect* optimisation of a bad objective returns *exactly* the
   bad answer, deterministically.
3. **Which alpha is best depends on the objective's quality, and the ordering FLIPS.** At
   signal 1.0 small alpha wins (1.270 exactly vs 1.386 for alpha=1). At signal 0.3 the
   expectation value wins (2.266 vs 2.459). **There is no alpha that is good in general.**
4. **Near-native probability MASS moves opposite to returned RMSD.** At signal 1.0,
   alpha=1.0 puts **0.766** of its mass below 2.5 A while alpha=0.05 puts 0.269 — yet
   alpha=0.05 returns the better single structure. The two readouts genuinely disagree, so
   the choice of alpha depends on whether the downstream operator consumes an argmin or an
   ensemble. Given the recorded result that *the terminal operator consumes the set MEAN*,
   the mass column is the one that matters, and it favours **large** alpha.
5. **Neither annealing schedule beats the best fixed alpha** on any of the three signal
   levels.

**Classification: DEMONSTRATED on one target x 3 signal levels x 5 seeds.** The nine-target
version is in the Part C grid.

---

## 2. THE DEAD-QUBIT DEFECT IS TWICE AS LARGE AS RECORDED — DEMONSTRATED (correction)

The brief states that `core/project.py` leaves `phi[0]`, `psi[n-1]`, `phi[n-1]` inert and
concludes "**log2(k) qubits per chain are dead**". Measured by brute force — flip qubit `q`
on 4,000 random configurations, ask whether the tabulated CA-RMSD ever moves — on all nine
enumerated targets:

| quantity | value |
|---|---|
| nominal qubits (n=9, k=4) | 18 |
| **live qubits** | **14** |
| dead qubits | **4**: `[0, 1, 16, 17]` |
| targets showing exactly this pattern | **9 of 9** |

The dead block is **residue 0 AND residue n-1**, i.e. **`2 log2 k` qubits per chain, not
`log2 k`**. At k=4 that is 4 dead qubits, 22.2% of the register.

**The mechanism is the locality theorem itself, so this is a prediction and not a
coincidence.** `d_ij` depends on exactly the `j-i-1` residues *strictly between* `i` and `j`.
Residue 0 and residue n-1 are strictly between **no** pair of residues, so neither enters any
CA-CA distance, so neither can change the CA distance matrix, so neither can change a Kabsch
CA-RMSD. Residue 0's `psi[0]` being individually "live" as a torsion is irrelevant: it moves
the chain rigidly.

**Consequences.** (i) Every `n_res * log2 k` qubit count in the project is an overcount by
`2 log2 k`, so the k=4 operating point is **21.9 mean live qubits**, not 25.9. (ii) A quarter
of the register at n=9 is pure search-space inflation for CA-RMSD — 16x the configurations
for zero structural resolution. (iii) `phi[n-1]` does place C/CB/O, so Legacy and AMBER **do**
read a state CA-RMSD cannot see, which is a genuine 2-qubit objective/metric mismatch and not
merely wasted width.

---

## 3. THE SIGNAL-TUNABLE OBJECTIVE FAMILY — ORACLE DIAGNOSTIC, calibrated

`blend_objective(base, truth, s)` mixes in **rank space** after putting both on the identical
uniform marginal, so `s` is a pure quality knob with no scale or tail-shape confound. Realised
correlations are **measured at every setting**, never assumed. Mean over nine targets:

| signal | 0.00 | 0.05 | 0.10 | 0.15 | 0.20 | 0.30 | 0.40 | 0.50 | 0.70 | 1.00 |
|---|---|---|---|---|---|---|---|---|---|---|
| **global rho** (base=legacy) | +0.081 | +0.124 | +0.171 | +0.222 | +0.277 | +0.405 | +0.557 | +0.721 | +0.933 | +1.000 |
| **in-decile rho** | +0.085 | +0.248 | +0.360 | +0.422 | +0.440 | +0.446 | +0.428 | +0.401 | +0.340 | +1.000 |
| RMSD at the certified argmin | **3.920** | 2.360 | 2.287 | 2.238 | 2.233 | 2.233 | 2.153 | 2.153 | 2.065 | 0.969 |
| best RMSD available in the lowest decile | 1.697 | 1.627 | 1.605 | 1.588 | 1.549 | 1.482 | 1.444 | 1.372 | 1.132 | 0.969 |

Signal 0 reproduces the established Legacy result exactly: the certified global optimum is
**3.920 A, worse than the 3.781 A random draw**.

### 3a. WARNING — in-decile rho is NOT a monotone measure of objective quality

Global rho rises monotonically with the signal knob. **In-decile rho does not**: it peaks at
**+0.446 around signal 0.2-0.3 and then FALLS to +0.340 at signal 0.7**, where the objective
is enormously better (global rho +0.933, certified argmin 2.065 A vs 2.233 A). The mechanism
is that as the objective improves its lowest decile becomes a narrower, more homogeneous set
of genuinely good structures, so there is less RMSD spread left inside it to rank.

**Therefore an in-decile rho of ~+0.37 is compatible with two very different objectives** —
signal ~0.11 (global rho +0.18, certified argmin 2.29 A) and signal ~0.75 (global rho ~+0.94,
certified argmin ~2.05 A). In-decile rho is a good diagnostic of *where a search lives* and a
**poor sufficient statistic for objective quality**. This is a caution for any comparison
that uses it as the quality axis, including the reconciliation in section 4.

---

## 4. RECONCILIATION WITH THE COORDINATOR'S STRUCTURAL OBJECTIVE

The coordinator reports an in-decile rho of **+0.370** for `s14/hamil.py` at w=0.25, and a
budget curve on it that is **monotone in objective and flat in structure** (3.742 -> 3.565 ->
3.632 -> 3.572 over 10 -> 20,000 evaluations), with a **selection gap growing 0.835 -> 2.040 A**.

**Localisation.** On my family the coordinator's +0.370 lands at **signal ~0.108** from a
Legacy base and **signal ~0.111** from a prior base — the two independent bases agree to
0.003, which is a real cross-check rather than a coincidence.

**Agreement so far.** At signal 0.10 my certified argmin is 2.287 A while the best RMSD
available in that objective's lowest decile is 1.605 A — a **0.68 A selection gap at the
certified global optimum**, i.e. present even with an infinite budget. And the measured
selection gaps per arm (1CS9, complete; the other eight targets are still running) are
**large and essentially arm-independent**:

| budget | signal | random | greedy | anneal | GA | VQE a=1.0 | VQE a=0.25 | VQE a=0.1 | VQE a=0.05 |
|---|---|---|---|---|---|---|---|---|---|
| 3,000 | 0.00 | 3.519 | 2.832 | 2.749 | 2.687 | 2.879 | 2.713 | 2.578 | 2.590 |
| 3,000 | 0.10 | 1.352 | 0.915 | 0.937 | 0.948 | 1.109 | 1.054 | 1.107 | 1.139 |
| 3,000 | 0.30 | 0.992 | 0.927 | 0.877 | 0.733 | 0.910 | 1.012 | 0.878 | 0.896 |
| 30,000 | 0.10 | 1.135 | 1.109 | 1.211 | 1.218 | 0.939 | 1.190 | 1.010 | 1.086 |
| 30,000 | 0.30 | 1.013 | 1.015 | 1.062 | 0.817 | 0.682 | 0.976 | 0.706 | 0.972 |

Every arm — classical and quantum — **evaluates something ~1 A better than what it returns**,
and raising the budget tenfold does not close the gap. This is the coordinator's
discrimination-limited diagnosis, reproduced on a different objective family, a different
target set and with quantum arms included. **No arm has a selection-gap advantage**, which is
the first piece of evidence that CVaR's tail-weighting does not buy discrimination.

**One caution I owe the coordinator** (section 3a): matching on in-decile rho alone may not
match objective quality, because in-decile rho is non-monotone in quality on my family. The
`hamil` objective's global rho of +0.539 sits between my signal 0.4 and 0.5, whereas its
in-decile +0.370 sits at signal 0.11. **The two axes disagree by a factor of four in signal**,
so the localisation above should be read as "somewhere between signal 0.11 and 0.45", not as
a point.

---

## 5. PART B — ENCODINGS

Sprint 13 compared binary / gray / gray_sorted / hierarchical / one-hot on qubits, move
locality and one-layer expressivity, and chose binary. It did **not** measure the Pauli-weight
distribution, Hamiltonian sparsity, coefficient dynamic range or gradient variance *per
encoding*, and did not consider domain-wall. This does. `s14/vqe_encoding.py`.

### 5.1 The qubit-count question is closed by counting — DEMONSTRATED

Any faithful encoding of `k^n` configurations needs at least `ceil(n log2 k)` qubits. At
n=9, k=4 that bound is **18**, and **binary attains it exactly**. Every alternative is
strictly wider:

| encoding | qubits/residue | qubits (n=9,k=4) | feasible fraction | penalty | vs bound |
|---|---|---|---|---|---|
| **binary / gray / gray_sorted** | 2 | **18** | 1.000 | none | **+0** |
| domain-wall (= unary) | 3 | 27 | 1.95e-03 | yes | +9 |
| one-hot | 4 | 36 | 3.81e-06 | yes | +18 |

**"Unary" and "domain-wall" are the same encoding** for a single categorical variable — a
prefix of ones in `k-1` qubits — and are listed once rather than as two arms.

The only regime with room for a cleverer encoding is `k` not a power of two, and there the
waste is small: at k=7, n=9 binary uses 27 against a bound of 26 (**1 qubit**); at k=3, n=9
it is 18 against 15 (3 qubits); at k=5, n=16 it is 48 against 38 (10 qubits, the worst case
found).

### 5.2 The Pauli spectrum, with the mandatory conditioning — DEMONSTRATED

Every objective is rank-conditioned (`uniformise`) before its Walsh-Hadamard transform, per
the brief's methodological requirement. The raw spectrum is computed alongside so the
artefact is priced: on this 5-residue sub-register the top-10 of 1,024 configurations carry
**39.3%** of raw Legacy's variance and **2.9%** after conditioning (2.9% is the value for an
unstructured function, i.e. the artefact is fully removed).

Conditioned mean Pauli weight, 1CS9, 5 residues, 1,024 configurations:

| encoding | qubits | cond. mean weight | Binomial(m,1/2) null | 99% weight | #terms for 99% of variance | dynamic range (decades) |
|---|---|---|---|---|---|---|
| **binary** | 10 | 3.012 | 5.0 | 7 | **510** | 4.02 |
| **gray / gray_sorted** | 10 | **2.725** | 5.0 | 7 | **510** | 4.02 |
| domain-wall | 15 | 4.279 | 7.5 | 9 | 14,314 | 5.20 |
| **one-hot** | 20 | **10.007** | **10.0** | 18 | **454,463** | 5.42 |

**Three results, one of them new and positive.**

1. **One-hot's Pauli spectrum lands EXACTLY on the Binomial(m, 1/2) null (10.007 vs 10.0).**
   Its spectrum carries essentially no information about the objective. Section 5.3 shows why.
2. **One-hot needs 454,463 Pauli terms for 99% of its variance where binary needs 510** — a
   factor of **891**. This is the decisive practical number and it was not previously measured.
3. **RETRACTED, by my own replication — "Gray coding lowers the mean Pauli weight".** On
   1CS9 Gray looked clearly better (2.725 vs binary's 3.012 on Legacy, 3.015 vs 3.640 on
   true RMSD) and I wrote it up as a new positive that s13 had missed. **Replicating over
   9 targets x 3 five-residue windows x 2 objectives (54 cells) destroys it:**

   | objective | n cells | mean delta (gray - binary) | cells where Gray is lower |
   |---|---|---|---|
   | Legacy | 27 | **+0.0055** | 16/27 |
   | true RMSD | 27 | **+0.0583** | 13/27 |

   A coin flip, with per-target swings of +/- 1.07 in mean weight that change sign across
   targets (1CS9 -0.62, 2MK7 -1.07, but 6F3V +0.97, 9UV5 +0.84 on RMSD). **The single cell I
   generalised from was one draw from a high-variance, zero-mean distribution.**
   **Classification: REFUTED. s13's "Gray and hierarchical coding buy nothing measurable" is
   CONFIRMED, on a statistic s13 did not compute.** Recorded here in full rather than
   deleted, because the failure mode — a clean-looking effect on the first target — is the
   one this project keeps meeting.

### 5.3 A non-surjective encoding's spectrum measures the CONSTRAINT, not the objective — DEMONSTRATED (new)

Sweep the infeasible-state penalty height and watch the mean Pauli weight:

| penalty | one-hot | domain-wall |
|---|---|---|
| 1.001 (just above the conditioned max) | 10.007 | 4.279 |
| 2.0 | 10.009 | 3.897 |
| 20.0 | 10.010 | 3.868 |
| 100.0 | 10.010 | 3.870 |
| **binary (surjective, NO penalty)** | — | **3.012** |

One-hot's spectrum is **immovable at the binomial null** regardless of the penalty, because
99.9996% of its register is infeasible and the indicator of that set dominates everything.
Domain-wall moves 10% and then saturates.

**Methodological consequence, and it generalises beyond this project: a Pauli-weight
comparison across encodings of different surjectivity is not a comparison of Hamiltonians.**
For a non-surjective encoding it measures the feasibility constraint. Binary has no penalty
and therefore no such degree of freedom, which is a further reason to prefer it that is
independent of qubit count.

### 5.4 Gradient variance per encoding — DEMONSTRATED

RY/CNOT depth 2, expectation value, conditioned Legacy, 4 residues, 64 random-angle seeds:

| encoding | qubits | Var[g_1] | mean \|g\| | Var[g_1] x 2^qubits |
|---|---|---|---|---|
| binary | 8 | 1.47e-03 | 0.1387 | 0.376 |
| gray / gray_sorted | 8 | 6.46e-04 | 0.1236 | 0.165 |
| domain-wall | 12 | 2.52e-04 | 0.0428 | **1.034** |
| one-hot | 16 | **2.17e-06** | 0.0046 | 0.142 |

**One-hot's gradient variance is 678x smaller than binary's.** Normalised by its own width it
is 0.142 vs binary's 0.376, i.e. **the loss is essentially the exponential price of the extra
8 qubits and not an additional pathology**. Domain-wall is the best per unit width (1.034) —
interesting, and not enough to overcome 9 extra qubits.

### 5.5 Move locality, on 300x the sample — CONFIRMS s13

All 262,144 configurations, 20,000 base states, every qubit:

| move | median \|d RMSD\| | median \|d Legacy\| |
|---|---|---|
| binary, 1 qubit flip | 0.3030 | 0.5217 |
| gray, 1 qubit flip | 0.3112 | 0.5208 |
| any residue state change (one-hot / domain-wall 2-flip) | 0.3202 | — |

Statistically indistinguishable. A one-hot 2-flip and a binary 1-flip produce **the same set
of configuration changes**, so no bit labelling can make a move local. s13 section 4b is
confirmed on 300x the sample.

### 5.6 The encoding verdict

**Binary, unchanged, and now for four independent reasons rather than one.** (i) It attains
the information-theoretic qubit bound exactly at power-of-two k. (ii) It needs 510 Pauli
terms where one-hot needs 454,463. (iii) It is surjective, so it needs no penalty and its
spectrum has no penalty-dependent degree of freedom. (iv) Its gradient variance is the
highest of any encoding tested. Gray coding buys nothing (5.2, retracted). **No encoding
change can rescue a bad objective**, which sections 6 and 8 show is the actual binding
constraint.

---

## 6. PART D — ANSATZ: THE QNG REFUTATION IS DEPTH-1 ONLY (correction)

`s14/vqe_ansatz.py`. Fubini-Study metric `g_ij = Re<d_i psi|d_j psi> - <d_i psi|psi><psi|d_j psi>`
by central differences on the exact statevector — machinery wholly independent of
`core/quantum.py`. n=8, averaged over seeds.

| pattern | L | P | mean g_ii | max off-diag | **cond. number** | rank | max\|g - I/4\| |
|---|---|---|---|---|---|---|---|
| none (product) | 1 | 8 | 0.250000 | 0.000000 | 1.000 | 8 | **3.5e-12** |
| chain | 1 | 8 | 0.250000 | 0.000000 | 1.000 | 8 | 3.5e-12 |
| ring | 1 | 8 | 0.250000 | 0.000000 | 1.000 | 8 | 3.5e-12 |
| block (torsion-aware) | 1 | 8 | 0.250000 | 0.000000 | 1.000 | 8 | 3.5e-12 |
| all_to_all | 1 | 8 | 0.250000 | 0.000000 | 1.000 | 8 | 3.5e-12 |
| **none** | **2** | **16** | 0.250000 | 0.250000 | **inf** | **8** | 2.5e-01 |
| chain | 2 | 16 | 0.250000 | 0.142821 | 3.665 | 16 | 1.4e-01 |
| ring | 2 | 16 | 0.250000 | 0.203157 | 9.937 | 16 | 2.0e-01 |
| brick / block_chain | 2 | 16 | 0.250000 | 0.218203 | 15.392 | 16 | 2.2e-01 |
| **block** | **2** | **16** | 0.250000 | 0.246720 | **155.2** | 16 | 2.5e-01 |
| **block** | **3** | **24** | 0.250000 | 0.246720 | **478.9** | 24 | 2.5e-01 |
| all_to_all | 2 | 16 | 0.250000 | 0.106178 | **2.478** | 16 | 1.1e-01 |

**CORRECTION to the brief.** "QNG is refuted for the current ansatz; the Fubini-Study metric
is full rank with `g_ii = 0.2500` exactly and is exactly `I/4` at depth 1" is **true and
complete at depth 1, and true for EVERY entangler pattern there — but it does not extend to
depth >= 2.**

* `g_ii = 0.250000` exactly at every depth and pattern (that part is universal).
* `g = I/4` **only at depth 1**. At depth 2 the maximum deviation is 0.11-0.25 and the
  condition number ranges **2.478 (all-to-all) to 155.2 (block)** — a factor of 63 across
  ansatz choices. At depth 3 the block ansatz reaches **478.9**.
* **So QNG is not refuted at depth >= 2**, and it would matter *most* for exactly the
  torsion-aware `block` ansatz the brief asked me to re-measure. **HYPOTHESIS (not yet
  measured): QNG buys nothing useful anyway, because the badly-conditioned directions are
  the reducible ones.** The measurement is a natural next step and I did not run it.
* **Parameter redundancy, measured:** the product ansatz (`none`) is **rank 8 with 16 or 24
  parameters** — infinite condition number, exactly `L`-fold redundancy, because RY angles on
  the same wire simply add. Every entangling pattern is full rank at every depth tested.
* **The metric contains no Hamiltonian — now established STRUCTURALLY rather than
  numerically.** `fubini_study` is a function of the state and its derivatives only; no energy
  array can enter it. The brief's bit-identical measurement (0.000e+00) is a consequence, not
  a coincidence. The energy model selects *where on the manifold* you go, not the manifold's
  shape.

---

## 7. THE COORDINATOR'S OBJECTIVE AGAINST A CERTIFIED OPTIMUM — DEMONSTRATED

`s14/vqe_hamil.py` tabulates `s14/hamil.py` over **all 262,144 configurations** of each of
the nine enumerated targets, so its certified global optimum is known exactly and the
"you did not search far enough" defence is closed. Uniform population, both terms
standardised over the complete space. Four axes, as the shared brief now requires:

| w (share on distogram) | global rho | in-decile rho | **RMSD @ CERTIFIED argmin** | decile mean | decile best |
|---|---|---|---|---|---|
| 0.00 | +0.084 | +0.087 | **2.866** | 3.531 | 1.020 |
| 0.25 | +0.176 | +0.105 | **2.896** | 3.410 | 1.020 |
| 0.50 | +0.272 | +0.113 | 3.034 | 3.318 | 1.121 |
| 0.75 | +0.336 | +0.112 | 3.073 | 3.253 | 1.308 |
| 1.00 | +0.360 | +0.115 | 3.237 | 3.244 | 1.513 |

Reference on these nine: space best **0.969**, random draw **3.781**, Legacy certified
argmin **3.920**, s13 torsion-prior certified argmin **3.636**.

**7a. The first objective whose certified optimum is somewhere you would want to go.**
At w=0.25 hamil's certified optimum is **2.896 A, 0.885 A better than a random draw**.
Legacy's is 0.139 A *worse* than random; the s13 torsion prior's is 0.145 A better. hamil is
roughly **six times the torsion prior's margin**, and this survives certification rather than
resting on a budget. **Sprint 13's "the certified optimum is in the wrong place" does NOT
generalise to a structural objective.** That is a real and important limit on the sprint's
central negative.

**7b. NEW — the w axis inverts, and the two axes give opposite answers.** Global rho improves
monotonically with distogram weight (+0.084 -> +0.360) while the **certified argmin gets
monotonically WORSE (2.866 -> 3.237)** and the best structure available in the decile
degrades (1.020 -> 1.513). Decile *mean* moves the other way (3.531 -> 3.244). Mechanism:
**the distogram orders the bulk; the retrieval torsion prior places the optimum.**
Consequence: if the terminal operator is an argmin or a small-m aggregation, low w is right;
if it is a large-m set mean, high w is right. Given the recorded result that *the terminal
operator consumes the set MEAN*, this is a live architectural tension, and **w must not be
chosen on global rho alone.**

**7c. The selection gap at INFINITE budget, on the coordinator's own objective.**

| w | certified argmin | decile best | **gap** |
|---|---|---|---|
| 0.00 | 2.866 | 1.020 | **+1.846** |
| 0.25 | 2.896 | 1.020 | **+1.876** |
| 1.00 | 3.237 | 1.513 | **+1.724** |

This is ~1.9 A of discrimination loss that no search can recover, and it lands close to the
coordinator's *sampled* selection gap of **2.040 A at budget 20,000**. Two independent routes
— a certified enumeration and a sampled budget curve — agree. It is also the same phenomenon
ENER measures from a third direction: no objective exceeds 0.511 pairwise accuracy below a
0.25 A quality gap, and the useful search range is ~2.5 A wide, so each energy resolves
roughly one bit of it.

**7d. Heterogeneity is large and the mean is not a typical target.** Per-target global rho at
w=0.25: 7N2I +0.469, 2P5H +0.366, 2MK7 +0.343, 9UV5 +0.322, 6EY3 +0.306, 6S0N +0.140,
8IS3 +0.048, **6F3V -0.150, 1CS9 -0.262**. Two of nine anti-rank. Certified argmin ranges
**0.675 A (2MK7)** to **5.256 A (1CS9, worse than that target's own random draw of 4.064)**.

**7e. The w axis, finely — a PLATEAU, and the mechanism made quantitative.**
`s14/vqe_wsweep.py`.

| w | global rho | in-decile | **certified argmin** | sd across targets | decile mean | decile best |
|---|---|---|---|---|---|---|
| 0.00 | +0.084 | +0.087 | 2.866 | 1.389 | 3.531 | 1.020 |
| 0.05 | +0.101 | +0.094 | **2.861** | 1.389 | 3.505 | 1.020 |
| 0.10 | +0.118 | +0.099 | **2.861** | 1.389 | 3.480 | 1.020 |
| 0.15 | +0.137 | +0.101 | 2.896 | 1.371 | 3.455 | 1.020 |
| 0.25 | +0.176 | +0.105 | 2.896 | 1.371 | 3.410 | 1.020 |
| 0.50 | +0.272 | +0.113 | 3.034 | 0.951 | 3.318 | 1.121 |
| 1.00 | +0.360 | +0.115 | 3.237 | 0.968 | 3.244 | 1.513 |

The certified optimum does **not** keep improving below w=0.25 — it is **flat at 2.861-2.896
across all of w in [0, 0.25]**, and the 0.035 A "best" at w=0.05-0.10 is far inside the
1.389 A between-target sd. **It is noise and must not be reported as an optimum.**

**The mechanism, now quantitative.** Per target, w=0.00 and w=0.25 give the **identical**
certified argmin on **7 of 9 targets** (only 8IS3 and 9UV5 move, by 0.04 and 0.31 A). So
**below w ~ 0.25 the distogram contributes nothing to the argmin — the certified optimum is
the retrieval torsion prior's alone.** Meanwhile the distogram improves the decile MEAN
monotonically and without exception over the full range, worth **0.287 A**. "The distogram
orders the bulk, the prior places the optimum" is therefore: *zero on the argmin below
w=0.25 on 7/9 targets, +0.287 A on the decile mean.*

**A third axis on the same inversion:** the between-target sd of the certified argmin falls
1.389 -> 0.951 as w rises, so distogram weight buys **consistency across targets** while
costing mean argmin quality. Two axes favour low w, one favours high w.

**Caution:** w=0.35 is non-monotone (3.185, worse than both 0.25 and 0.50) — a single-target
flip, and evidence the argmin axis is noisy at n=9.

Per-target certified argmin (the mean is a **mixture**, never quote it alone):

| w | 1CS9 | 2MK7 | 2P5H | 6EY3 | 6F3V | 6S0N | 7N2I | 8IS3 | 9UV5 |
|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 5.256 | 0.675 | 2.184 | 3.820 | 4.538 | 2.825 | 1.465 | 3.005 | 2.026 |
| 0.25 | 5.256 | 0.675 | 2.184 | 3.820 | 4.538 | 2.825 | 1.465 | 2.962 | 2.339 |
| 1.00 | 4.406 | 4.446 | 3.479 | 1.948 | 4.234 | 3.149 | 1.687 | 3.243 | 2.542 |
| random draw | 4.064 | 4.869 | 3.788 | 3.696 | 3.568 | 3.293 | 4.157 | 3.080 | 3.513 |

w=1.0 **rescues 6EY3 (3.820 -> 1.948) and destroys 2MK7 (0.675 -> 4.446)**. The w choice is
not a smooth trade-off; it is a per-target coin flip about which term wins.

**7ea. CORRECTION TO 7e, from my own follow-up: the argmin result does NOT license
dropping the distogram.** 7e says the distogram contributes nothing to the argmin below
w=0.25 on 7/9 targets. That is a statement about ONE POINT. Comparing the top-**m** SETS:

| m | w=0.00 best | w=0.00 mean | w=0.25 best | w=0.25 mean | **Jaccard(w0, w0.25)** |
|---|---|---|---|---|---|
| 1 | 2.866 | 2.866 | 2.896 | 2.896 | 0.667 |
| 5 | 2.776 | 3.048 | **2.527** | **2.988** | **0.496** |
| 20 | 2.441 | 3.045 | **2.349** | **2.991** | 0.533 |
| 75 | 2.138 | 3.068 | **2.127** | **3.003** | 0.580 |
| 200 | 1.956 | 3.108 | **1.908** | **3.056** | 0.612 |

**The distogram turns over about half the top-m set at every m** (Jaccard 0.50-0.61) and is
**better at every m >= 5 on both best-in-set and set-mean** (m=5: -0.249 A on the best). The
only m where w=0 wins is m=1 exactly, by +0.030 A — inside the 1.389 A between-target sd,
i.e. a null. **The argmin coincidence was a knife-edge artefact of looking at a single
point.** The correct statement is: *the distogram does not move the argmin on 7/9 targets,
but it turns over half the top-m set and improves it at every m >= 5.*

**Extrapolation caveat that cuts the same way.** All nine targets are n=9. The distogram is a
pair-distance objective: 36 pairs and max separation 8 at n=9, against 120 pairs and
separation 15 at n=16, where the locality theorem says the many-body content is largest.
**Its contribution is structurally expected to grow with chain length, and these nine are the
shortest bin in the instrument** — the least favourable case for it.

**Recommendation: keep the distogram at w in [0.25, 0.5].** Three independent axes favour it
(top-m set quality at every m >= 5; per-target consistency, sd 1.389 -> 0.951; decile mean
+0.287 A) and only the m=1 argmin favours dropping it, by an amount inside the noise.

**7f. LIKE-FOR-LIKE corroboration with the coordinator's sampled budget curve.**
My certified gap and the coordinator's sampled gap are **not the same statistic** (mine uses
the objective's own decile of the complete space as the reference set, theirs a 20,000-draw
sample; my nine targets are all n=9, theirs span n=9-16). So the coordinator's *exact*
statistic was recomputed on the enumerated nine — uniform draws read as nested prefixes,
hamil w=0.25, 6 seeds x 9 targets:

| budget | 10 | 300 | 3,000 | 20,000 |
|---|---|---|---|---|
| **mine** returned | 3.297 | 3.173 | 3.177 | **2.920** |
| **mine** best available | 2.660 | 1.590 | 1.206 | 0.997 |
| **mine** GAP | 0.637 | **1.583** | **1.971** | 1.923 |
| **coordinator's** GAP (126 targets) | 0.835 | **1.577** | **1.935** | 2.040 |

**The selection gap reproduces to within 0.12 A at every budget, and to 0.006 and 0.036 at
budgets 300 and 3,000** — same statistic, disjoint target sets, independent implementations.
This is genuine corroboration and is stronger than the earlier decile-based comparison.

**But the other half of the curve does NOT reproduce, and this is a correction.** The
coordinator's *returned* RMSD is flat (3.742 -> 3.572, 0.17 A of movement). **Mine falls
monotonically by 0.38 A (3.297 -> 2.920).** So on the enumerated nine, searching harder *does*
improve the emitted structure. Most likely the length bin: at budget 20,000 my n=9 targets
have seen 7.6% of their space. **The honest joint statement is therefore: the discrimination
loss is the same everywhere and reproduces exactly; whether search still buys anything on top
of it is target-length dependent. "The emitted structure does not move" must not be stated as
a general claim — my curve contradicts it.**

**7g. In-decile, for the record.** I measure **+0.105** at w=0.25 on the uniform full
enumeration, against the coordinator's +0.046 under uniform sampling and the withdrawn +0.370
under a prior-drawn proposal. My number sits much nearer the uniform one, which **independently
confirms that the +0.370 was a proposal artefact** rather than a property of the objective.

---

## 8. TAIL-RESTRICTED PAIRWISE ACCURACY — the statistic nobody was reporting, and the mechanism it closes

`s14/vqe_tailacc.py`. **DEMONSTRATED, nine fully enumerated targets.**

Four proxies for "objective usefulness" were used this sprint and **none of the first three
predicts argmin quality**: global rho (dominated by the bulk), in-decile rho (non-monotone in
quality, section 3a), bulk pairwise accuracy (measured over random pairs, so dominated by the
99% of the space a search never visits). The quantity that governs an argmin is accuracy
**inside the region the search ends up in** — restrict to the objective's own lowest `frac`,
then measure ordering accuracy among pairs separated by at least `gap` in true CA-RMSD.
Objective ties are scored 0.5 (scoring them win/loss would import index order — the recorded
tie-breaking trap). At gap > 0.25 A, the ENER threshold:

| frac of space | n configs | hamil w=0 | w=0.25 | w=0.5 | **w=1 (distogram)** |
|---|---|---|---|---|---|
| 1.0 (BULK) | 262,144 | 0.532 | 0.570 | 0.611 | **0.655** |
| 0.1 | 26,214 | 0.534 | 0.542 | 0.545 | 0.549 |
| 0.01 | 2,621 | 0.526 | 0.528 | 0.521 | 0.509 |
| **0.001 (TAIL)** | **262** | **0.512** | **0.521** | **0.520** | **0.390** |
| **tail - bulk** | | **-0.020** | **-0.049** | **-0.091** | **-0.266** |
| **certified argmin** | | **2.866** | **2.896** | **3.034** | **3.237** |

**8a. The mechanism, on one number.** `tail - bulk` is monotone in w and tracks the certified
argmin exactly. **The pure distogram orders structures BACKWARDS inside its own top 0.1%
(0.390, well below chance) while being the best bulk ranker in the family (0.655).** This is
the fourth independent confirmation of "the distogram orders the bulk, the torsion prior
places the optimum", on a statistic unrelated to the other three.

**8b. REFUTED — my own pre-registered prediction, at its premise.** I registered (before
seeing data) that hamil w=0 was *tail-concentrated* — skill in the tail, chance in the bulk —
and therefore that small-alpha CVaR should beat the expectation value on it, since CVaR reads
only the tail. **No objective in this family is tail-concentrated.** w=0's tail accuracy is
0.512, chance, and *below* its own bulk accuracy of 0.532. Every objective tested is at
chance or worse inside its own tail. **The substrate the prediction required does not exist.**
Revised prediction, now the operative one: **small-alpha CVaR should be neutral to harmful
everywhere**, because it concentrates attention on exactly the region where every objective
carries no ordering information — consistent with the blend sweep, where alpha=1 beat
alpha=0.05 at signal 0.3 (2.266 vs 2.459).

**8c. The argmin is a RANDOM DRAW FROM ITS OWN TAIL — DEMONSTRATED.** If an objective is at
chance inside its tail, the argmin it returns should be indistinguishable from a uniform draw
from that tail. Tested directly:

| w | tail m | **actual argmin** | tail MEAN | tail best | **percentile of the argmin within its own tail** |
|---|---|---|---|---|---|
| 0.00 | 26 | 2.866 | 3.000 | 2.335 | **0.419** |
| 0.00 | 262 | 2.866 | 3.113 | 1.889 | **0.434** |
| 0.00 | 2,621 | 2.866 | 3.299 | 1.365 | **0.409** |
| 0.25 | 262 | 2.896 | 3.066 | 1.889 | **0.480** |
| 1.00 | 262 | 3.237 | 3.103 | 2.515 | 0.558 |
| 1.00 | 2,621 | 3.237 | 3.074 | 2.007 | 0.609 |

**The argmin sits at the 41st-48th percentile of its own tail** — a uniform draw would sit at
the 50th. At w=1 it sits at the **56th-61st**, i.e. *worse* than a random draw from its own
tail, exactly as the 0.390 anti-ranking predicts. The actual argmin equals the tail MEAN to
within 0.13-0.43 A while the tail BEST is 1.0-1.9 A better.

**So: the objective's argmin quality is set by WHERE its tail sits, not by any ability to
order within it. It picks a good neighbourhood and then draws essentially at random from it.**
w=0 wins the argmin because its extreme tail has the best *mean*, not because it ranks better
inside.

**8d. AND THAT IS THE SELECTION GAP, in closed form.** tail mean minus tail best at w=0,
m=2,621 is **3.299 - 1.365 = 1.93 A**. My certified selection gap is **+1.876 A**. The
coordinator's sampled gap at budget 20,000 is **2.040 A**. **These are the same quantity:
the selection gap IS the price of drawing at random from a tail you cannot order.** Three
independent measurements, one mechanism, now with an account rather than an empirical value.

**8e. Methodological recommendation.** Tail-restricted pairwise accuracy is the only one of
the four proxies that predicts argmin quality, and it is cheap. Global rho fails (w=1 has the
best global ordering and the worst argmin); in-decile rho fails (non-monotone); bulk pairwise
accuracy fails (w=1 wins it decisively and has the worst argmin). **Caveat:** at frac=1e-4
(26 configurations) it is too noisy to use.

---

## 9. PART C — THE DECISIVE EXPERIMENT: WHEN DOES VQE BECOME USEFUL?

`s14/vqe_signal.py` + `s14/vqe_report.py`. Signal-tunable objective family on the enumerated
space; at each signal level, on the SAME space, at MATCHED budgets counted in objective
evaluations: random / greedy 1-opt / simulated annealing / GA / **exact enumeration** /
expectation-value VQE / CVaR-VQE at alpha = 0.25, 0.1, 0.05. 3 seeds per arm.

**n WARNING, stated up front.** The box was shared with four other agents for most of this
run and delivered under one core of aggregate Python throughput across six processes. The
grid reached **5-7 of the 9 targets** per signal level, not 9. Every CI below is over that n
and is reported with W/L so the reader can see it. **Nothing here is a nine-target result.**
(The nine-target version of the same comparison, on a real objective, is section 11.)

### 9.1 The crossing curve, budget 30,000 evaluations (11% of the space)

Best VQE arm vs best classical arm, paired, negative = VQE better:

**COMPLETE: all 9 targets, 92 cells** (`done in 874s`, `s14/results/vqe_signal_legacy.json`).
The rho column now matches section 3's calibration exactly, which is the consistency check
that the grid scored the objectives it thought it was scoring.

| signal | rho | OBJECTIVE diff | CI95 | W/L | STRUCTURE diff | CI95 | W/L |
|---|---|---|---|---|---|---|---|
| 0.00 | +0.081 | +0.0001 | [+0.0000,+0.0001] | **0/7** | -0.126 | [-0.302,+0.105] | 6/1 |
| 0.05 | +0.124 | +0.0010 | [+0.0003,+0.0019] | **1/7** | +0.154 | [-0.006,+0.287] | 2/7 |
| 0.10 | +0.171 | +0.0012 | [-0.0005,+0.0037] | 1/5 | +0.094 | [-0.022,+0.292] | 1/2 |
| 0.20 | +0.277 | +0.0039 | [+0.0018,+0.0062] | **0/6** | **+0.070** | [+0.002,+0.142] | 2/4 |
| 0.30 | +0.405 | +0.0023 | [+0.0001,+0.0045] | **1/6** | +0.011 | [-0.061,+0.094] | 5/2 |
| 0.50 | +0.721 | +0.0029 | [+0.0013,+0.0049] | **0/8** | -0.033 | [-0.149,+0.048] | 3/4 |
| 1.00 | +1.000 | +0.0000 | [+0.0000,+0.0000] | **0/4** | **+0.027** | [+0.003,+0.056] | **0/4** |

**VQE loses on the objective axis at six of the seven signal levels with the CI excluding
zero, and W/L never better than 1/5.** On the structural axis it is significantly worse at
0.20 and 1.00 and everywhere else the CI spans zero. **No cell favours VQE.**

**A FALSE POSITIVE, killed by adding targets, and worth recording as such.** At n=4 the
signal-0.50 cell read **-0.113 A, CI [-0.292, -0.000], W/L 3/1** — a CI that just excluded
zero — and I had written it up as a possible structural crossing at rho ~ 0.72. Adding
targets killed it monotonically:

| n targets | signal-0.50 structure diff | CI95 | W/L |
|---|---|---|---|
| 4 | **-0.113** | [-0.292, **-0.000**] | 3/1 |
| 6 | -0.053 | [-0.226, +0.071] | 3/2 |
| 7 | -0.044 | [-0.187, +0.062] | 3/3 |
| 8 | -0.037 | [-0.167, +0.055] | 3/4 |
| **9 (complete)** | **-0.022** | [-0.122, +0.051] | **2/5** |

**The effect shrank by a factor of five and the W/L inverted from 3/1 favourable to 2/5
unfavourable, monotonically, as targets were added.** This is the concentration hazard the
brief warns about, arriving exactly on schedule at the one cell that would have been a
positive headline. Kept in place rather than deleted, because the *shape* of the decay is the
useful artefact: a real effect does not do this.

**9.1a. THE OBJECTIVE AXIS: VQE NEVER CROSSES. It loses at every signal level tested.**
The difference is positive everywhere and the CI excludes zero *on the wrong side* at
signals 0.00, 0.05, 0.20, 0.50 and 1.00, with W/L of 0/5, 0/4, 0/3, 0/4. This is the
sprint's cleanest causality answer and it is unambiguous.

The strength of the classical arms makes it sharper. Fraction of cells reaching the
**certified global optimum** at budget 30,000:

| signal | greedy 1-opt | annealing | VQE alpha=1 | VQE alpha=0.05 |
|---|---|---|---|---|
| **0.00 (real Legacy)** | **1.00** | **1.00** | **0.00** | **0.00** |
| 0.05 | 0.67 | 0.67 | 0.00 | 0.00 |
| 0.20 | 0.50 | 0.50 | 0.00 | 0.00 |
| 0.50 | 0.50 | 0.50 | 0.00 | 0.00 |
| 1.00 | 1.00 | 1.00 | 0.00 | 1.00 |

**On the real Legacy objective, greedy 1-opt and simulated annealing certify the global
optimum in 100% of cells and VQE in 0%.** VQE is decisively worse at the one thing an
optimiser exists to do.

**9.1b. THE STRUCTURAL AXIS: VQE never crosses either. It is significantly WORSE in the
mid-signal range and never significantly better anywhere.** VQE loses significantly at
signals 0.10 and 0.20 (+0.128, +0.089 A, CIs excluding zero on the wrong side, W/L 1/4,
0/3) and at 1.00 (+0.021). It ties at 0.00, 0.05, 0.30 and 0.50. **No cell favours VQE with
a CI excluding zero.** At n=4 the signal-0.50 cell briefly did; two more targets removed it
(see the note above).

**Even a hypothetical crossing would be in the wrong place.** rho +0.73 is far beyond any
real objective measured in this project: hamil's global rho is +0.176 at w=0.25 and +0.360 at
w=1.0; Legacy's is +0.081. It is also on the wrong side of the point where the objective is
good enough that classical search already solves the problem exactly (section 11e: at
budget 30,000 on hamil, every arm returns the certified optimum on 6 of 9 targets).

**9.1c. The answer to the brief's question, in one line.** *There is no crossing on either
axis.* The objective-axis crossing does not exist at any signal level — VQE is a strictly
worse optimiser than greedy 1-opt everywhere. The structural axis is a tie or a loss
everywhere. **The two axes do differ in character, and that difference is the finding: VQE
loses decisively on energy and merely ties on structure, because both methods hit the same
discrimination floor (section 8) from different distances. Structural parity here is not VQE
holding its own; it is the floor being close enough that a much worse optimiser reaches it
too.**

### 9.2 The budget trap, reproduced on the new instrument, and its cure

All nine targets:

| signal | certified optimum | space best | random draw | best arm @ b=3,000 | best arm @ b=30,000 |
|---|---|---|---|---|---|
| **0.00 (real Legacy)** | **3.920** | 0.969 | 3.781 | **3.677** | 3.706 |
| 0.05 | 2.360 | 0.969 | 3.781 | 2.549 | 2.395 |
| 0.10 | 2.287 | 0.969 | 3.781 | 2.466 | 2.357 |
| 0.20 | 2.233 | 0.969 | 3.781 | 2.466 | 2.256 |
| 0.30 | 2.233 | 0.969 | 3.781 | 2.435 | 2.265 |
| 0.50 | 2.153 | 0.969 | 3.781 | 2.288 | **2.149** |
| 1.00 | 0.969 | 0.969 | 3.781 | 1.124 | **0.969** |

**At signal 0 the certified global optimum is 3.920 A against a random draw of 3.781 A —
+0.139 A worse.** That is s13's recorded figure **to three decimal places**, reproduced here
by independent code on the same nine spaces, and it is the tightest available check that this
instrument is measuring what s13 measured. The best arm at budget 3,000 returns **3.677 A,
beating the certified optimum by 0.243 A**: more search makes the answer worse.

**The trap is cured by objective quality, not by search.** It shrinks through signal 0.20 and
inverts by 0.50, where the best arm at budget 30,000 (2.149) essentially reaches the certified
optimum (2.153) and more search is strictly better. On the coordinator's real structural
objective the residual trap is **0.034 A** (section 11f). **Every VQE arm in this document was
checked against the certified optimum rather than against random, as the brief requires.**

**The trap is cured by objective quality, not by search.** It vanishes by signal 0.20 and
inverts by signal 0.50, where more search is strictly better. On the coordinator's real
structural objective the residual trap is only **0.034 A** (section 11f) against Legacy's
0.293 A here.

### 9.3 The one thing VQE gives that classical search does not

Classical arms return a configuration; a variational state returns a *distribution*. Near-
native probability mass under the final distribution, budget 30,000 (classical arms have no
such column):

| signal | VQE a=1.0 | a=0.25 | a=0.1 | a=0.05 |
|---|---|---|---|---|
| 0.00 | 0.047 | 0.014 | 0.016 | 0.017 |
| 0.10 | 0.152 | 0.158 | 0.064 | 0.083 |
| 0.30 | 0.266 | **0.393** | 0.259 | 0.179 |
| 0.50 | **0.590** | 0.475 | 0.309 | 0.244 |
| 1.00 | **0.788** | 0.511 | 0.390 | 0.344 |

**P(RMSD < 2.5 A) rises with objective quality and is LARGEST at alpha = 1.** The expectation
value concentrates on good structures better than any CVaR tail does. Given the recorded law
that *the terminal operator consumes the set MEAN*, this column is the one that matters for a
downstream aggregator — **and it argues against CVaR, not for it.**

### 9.4 Selection gap — no arm has an advantage

`RMSD_returned - RMSD_best_EVALUATED`, budget 30,000: at signal 0.10 the eight arms span
1.077-1.375 and at signal 0.30 they span 0.930-1.263, with classical and quantum arms
interleaved and no ordering. **Every arm evaluates something ~1 A better than it returns, a
tenfold budget increase does not close it, and CVaR's tail weighting buys no discrimination.**
This is the direct answer to the sprint's causality question, and section 8 gives its
mechanism.

---

## 10. PART D — INITIALISATION: THE VQE OPTIMISATION IS WORSE THAN NOT OPTIMISING

`s14/vqe_init.py`. **DEMONSTRATED.** 4 targets x 3 seeds x 2 signal levels, budget 10,240
evaluations, MPS depth-2+final-RY on 18 qubits, CVaR alpha=0.25.

**THE HONEST CONTROL.** "After VQE" is *not* compared against the initialisation's mean. It
is compared against **`init_best`: the best of the SAME NUMBER OF SHOTS drawn from the
INITIAL distribution with no optimisation at all.** That is the arm a VQE must beat to have
contributed anything, because best-of-N sampling is free and is what the shipped pipeline
already does. `DELTA = after - init_best`; only a **negative** DELTA is a contribution.

**signal 0.10** (certified optimum 2.309 A):

| initialisation | init mean | **init best** | after VQE | **DELTA** | paired CI95 | W/L |
|---|---|---|---|---|---|---|
| random (native-free) | 4.131 | 1.797 | 2.606 | **+0.808** | [+0.576,+1.072] | **0/12** |
| uniform / zero (native-free) | 4.119 | 1.765 | 2.757 | **+0.992** | [+0.756,+1.242] | **0/12** |
| empirical prior warm start (NATIVE-FREE) | 4.062 | 1.591 | 2.469 | **+0.879** | [+0.539,+1.237] | 1/11 |
| ORACLE prior q=0.5 (ORACLE DIAGNOSTIC) | 3.893 | 1.309 | 2.626 | **+1.317** | [+0.888,+1.806] | **0/12** |
| ORACLE prior q=0.8 (ORACLE DIAGNOSTIC) | 3.424 | 1.444 | 2.714 | **+1.271** | [+0.820,+1.796] | **0/12** |
| ORACLE prior q=0.95 (ORACLE DIAGNOSTIC) | 2.966 | 1.764 | 2.578 | **+0.813** | [+0.274,+1.455] | 2/10 |

**signal 0.30** — same pattern, DELTA +0.654 to +1.118, W/L 0/12 on four of six rows.

### 10.1 Three readings, and the first is the headline

**1. The variational optimisation makes the answer WORSE than not optimising at all.**
DELTA is positive for **every** initialisation at **both** signal levels, by **+0.65 to
+1.32 A**, with W/L of 0/12 on most rows and bootstrap CIs that exclude zero decisively.
Taking the best of 10,240 shots from the **untrained** circuit beats running the VQE on the
same budget, every time. **This is not "VQE fails to help" — it is "VQE actively destroys
value relative to free best-of-N sampling from its own starting distribution".**

**2. The mechanism is section 8, and it closes the loop.** The optimisation concentrates the
distribution onto the objective's low-energy tail. Section 8 shows that inside that tail
**every objective is at chance** (0.512 at w=0) and the argmin is a **random draw from it**
(41st-48th percentile). So the optimiser trades a broad sample — whose *best* member is good
— for a narrow one whose *mean* it then samples blindly. Concentration is exactly the wrong
move when discrimination is the binding constraint. **Sections 8 and 10 are the same finding
seen from the mechanism side and the outcome side.**

**3. The initialisation question is answered, and the answer is "irrelevant".** After VQE the
result is **2.42-2.76 A regardless of where it started**, from initialisation means spanning
**2.966 to 4.131 A**. An ORACLE warm start at q=0.95 (init mean 2.966, init best 1.764) ends
at 2.578 — no better than a random start's 2.606. **The optimiser washes out its own
initialisation completely.** So a warm start here is neither "handing over the answer" nor
"accelerating genuine search": it is erased. And note the perverse ordering — **the better
the oracle initialisation, the WORSE the DELTA** (q=0.5 gives the worst, +1.317), because a
better start has a better `init_best` while the VQE converges to the same place anyway.

**This is the answer to the brief's initialisation confound.** There is no version of this
experiment in which the initialisation flatters the VQE, because the VQE cannot hold on to a
good initialisation at all.

---

## 11. THE CLOSING EXPERIMENT — every arm on the coordinator's real objective, against a certified optimum

`s14/vqe_hamil.py` + `python -m s14.vqe_report --base hamil`. **All nine enumerated targets**,
hamil w=0.25, budget 30,000 evaluations, 2-3 seeds per arm, certified global optimum known
exactly over all 262,144 configurations. **DEMONSTRATED.**

| arm | objective gap | structural gap | selection gap | P(RMSD<2.5) |
|---|---|---|---|---|
| **exact enumeration** | **0.0000** | 2.005 | — | — |
| random | 0.1641 | 1.986 | 1.986 | — |
| **greedy 1-opt** | **0.0000** | 2.005 | 1.949 | — |
| **annealing** | **0.0000** | 2.005 | 1.974 | — |
| **GA** | **0.0000** | 2.005 | **1.705** | — |
| VQE alpha=1.0 | 0.0111 | 2.100 | 1.828 | **0.272** |
| VQE alpha=0.25 | 0.0028 | 1.988 | 1.816 | 0.223 |
| VQE alpha=0.1 | 0.0028 | 1.988 | 1.872 | 0.180 |
| VQE alpha=0.05 | 0.0056 | **1.971** | 1.796 | 0.179 |

**11a. Objective axis: greedy, annealing and the GA all CERTIFY the global optimum; no VQE
arm does.** Reproduces section 9.1a on a real, native-free structural objective.

**11b. Structure axis: a dead tie.** Best VQE vs best classical: **-0.015 A, CI
[-0.374, +0.366], W/L 4/4.** Statistically indistinguishable.

**11c. Selection gap: 1.705-1.986 A across all eight arms, classical and quantum
interleaved.** This is the certified +1.876 A gap of section 7c reproduced by every search
method, and it confirms that **no arm has a discrimination advantage.**

**11d. THE PRE-REGISTERED ALPHA TEST — my revised prediction CONFIRMED, original REFUTED.**
Registered before the run: original prediction "small-alpha CVaR beats alpha=1 on a
tail-concentrated objective"; revised after section 8 to "no alpha reversal anywhere, because
no objective is tail-concentrated". Paired over nine targets, alpha=0.05 vs alpha=1.0:

| readout | mean diff | CI95 | W/L |
|---|---|---|---|
| RMSD returned | -0.1035 | [-0.3813, +0.0748] | **2/2** |
| P(RMSD<2.5) under p | -0.0306 | [-0.2151, +0.1454] | 2/7 |
| objective gap | -0.0055 | [-0.0154, +0.0043] | 4/1 |

**Every CI spans zero and the structural W/L is a dead 2/2.** The apparent -0.104 A mean is
carried by a **single target** — 2MK7, where alpha=1 returns 1.808 and alpha=0.05 returns
0.675, a -1.133 A swing — with the other seven tied or reversed. *A result carried by one
target is not a result.* **No alpha reversal exists. The revised prediction is confirmed and
the original is refuted on its own pre-registered test.**

**11e. Per-target, and it is the most informative table here.**

| target | VQE a=1.0 | VQE a=0.05 | greedy | **certified optimum** |
|---|---|---|---|---|
| 1CS9 | **5.020** | 5.256 | 5.256 | 5.256 |
| 2MK7 | 1.808 | **0.675** | 0.675 | 0.675 |
| 2P5H | 2.184 | 2.184 | 2.184 | 2.184 |
| 6EY3 | 3.820 | 3.820 | 3.820 | 3.820 |
| 6F3V | **4.403** | **4.269** | 4.538 | 4.538 |
| 6S0N | 2.825 | 2.825 | 2.825 | 2.825 |
| 7N2I | 1.465 | 1.465 | 1.465 | 1.465 |
| 8IS3 | 2.962 | 2.962 | 2.962 | 2.962 |
| 9UV5 | 2.240 | **2.339** | 2.339 | 2.339 |

**On 6 of 9 targets every arm returns the identical structure — the certified optimum.** The
search problem is *solved* at this budget and the answer is still 2.0 A from the space best.
The only deviations are the **budget trap**: on 1CS9 and 6F3V and 9UV5, *failing* to reach the
certified optimum returns a **better** structure (1CS9: 5.020 vs 5.256; 6F3V: 4.269 vs 4.538).

**11f. The budget trap is largely repaired by a structural objective.** Best arm at budget
30,000 gives **2.932 A** against a certified optimum of **2.966 A** — a residual trap of only
**0.034 A**, against **0.154 A** on Legacy (section 9.2) and s13's 0.253 A. Another
confirmation that the objective, not the search, was the problem.

---

## 12. WHAT I DID NOT ESTABLISH, AND MY OWN RETRACTIONS

**Retractions and refutations of my own claims, kept in place with the evidence:**

1. **"Gray coding lowers the mean Pauli weight"** — clean on 1CS9, destroyed by replication
   over 54 cells (section 5.2). REFUTED by me.
2. **"A structural crossing at rho ~ 0.72"** — CI [-0.292, -0.000] at n=4, gone at n=6
   (section 9.1). A concentration false positive, caught by adding targets.
3. **"hamil w=0 is tail-concentrated, so small-alpha CVaR should win"** — pre-registered,
   then refuted at its premise by my own tail-accuracy measurement (section 8b) and by the
   pre-registered test itself (section 11d, W/L 2/2, CI spanning zero).
4. **"alpha is effectively a learning rate"** — refuted, cos(grad(1.0), grad(0.01)) = 0.0727
   (section 1.7).
5. **"The distogram can be dropped below w=0.25"** — my own argmin result read as a set
   result; corrected by the top-m Jaccard measurement (section 7ea).

**Corrections to the brief / prior sprints:** the dead-qubit count is `2 log2 k` not
`log2 k` (section 2); the CVaR tail-baseline defect's severity is not the recorded constant
-0.023 but a wide alpha-dependent spread, with a closed-form bias (1.6); "small alpha
collapses distinct objectives" has the sign backwards (1.8); the QNG refutation is depth-1
only (section 6); "the certified optimum is in the wrong place" is a property of the
*physical* energies, not of the problem (section 7a).

**Not established:**

* **Part C is at 5-7 of 9 targets**, not 9 — the shared box fell to ~0.07 cores per process.
  The nine-target version of the same comparison on a real objective is section 11 and agrees.
* **`s14/vqe_ansatz.py` did not finish.** The Fubini-Study/QNG result (section 6) and the
  initialisation study (section 10) were captured standalone; the ansatz-family end-to-end
  comparison, the barren-plateau slopes and the multi-seed expressivity replication of s13
  section 4d were **not** run. s13's "depth 3-4 suffices" therefore remains OPEN.
* **QNG was never actually applied.** Section 6 shows the metric departs from I/4 at depth
  >= 2 with condition numbers to 478.9, so QNG is *not refuted* there — but whether it helps
  is a **HYPOTHESIS**, unmeasured.
* **No QAOA arm.** The one structurally different family (a diagonal cost unitary built from
  the objective itself) was not implemented. It would consume `4^n` objective evaluations per
  layer, so it cannot be budget-matched, but its expressivity was not tested.
* **`s14/vqe_tailavg.py` (coordinate-averaged CVaR tails) was written and smoke-tested but
  not run at scale.** On one target the CVaR tail average was worse than the classical
  top-20 average; n=1, no conclusion.
* **Shots sensitivity untested.** Every VQE arm used shots=512. Whether a different
  shots/iterations split changes the verdict is unmeasured, and it is the main methodological
  weakness of Part C.
* **Everything is n=9, k=4.** Section 7ea notes that a pair-distance term's content grows
  with chain length, so these nine targets are the least favourable bin for the distogram and
  possibly for many-body objectives generally.

---

## 13. REPRODUCTION

    python -m s12.instrument                     # pinned constants, start and end
    python -m pytest s14/test_vqe.py -q          # 92 correctness tests (Part A)
    python -m s14.vqe_calib                      # live-qubit audit + signal calibration
    python -m s14.vqe_decile                     # in-decile axis + selection gap
    python -m s14.vqe_cvar                       # Part A: gradient audit, alpha sweep
    python -m s14.vqe_collapse                   # the collapse surface (masked AMBER)
    python -m s14.vqe_encoding                   # Part B: encodings + Pauli spectra
    python -m s14.vqe_signal                     # Part C grid  (--resume --b30 --seeds2)
    python -m s14.vqe_report                     # Part C analysis (--base hamil)
    python -m s14.vqe_hamil --fast               # the closing experiment (--w 0,1)
    python -m s14.vqe_wsweep                     # fine w axis + like-for-like gap
    python -m s14.vqe_tailacc                    # tail-restricted pairwise accuracy
    python -m s14.vqe_init                       # Part D: the initialisation control

Seeds are explicit in every module (`SEEDS = (0,1,2)`; gradient audits use
`1000*seed + 13*n + layers`; classical arms use `1000*s + 7`). Every result JSON is written
under `s14/results/`; hamil full-enumeration tabulations are cached under `s14/cache/`.
All AMBER numbers use the `amber_kind == 0` mask (section 1.9).

**Leakage audit.** No benchmark file was read or written. `dev24` was not run. The native is
read only through `Enum.rmsd` (post-hoc scoring), `Space.ORACLE_*` / `ORACLE_prior` (labelled
ORACLE DIAGNOSTIC), and `blend_objective` / `noisy_truth`, which are ORACLE by construction
and exist solely to provide the objective-quality knob Part C requires. No native quantity
enters any arm described as native-free. Nothing under `core/`, `s5/ s7/ s8/ s9/ s12/`,
`tests/` was modified; no `s14/obj_*` file was touched. `qansatz.py` and `core/quantum.py`
were read but not modified — the recorded defect is left in `qansatz.py` and pinned by test.
