# Sprint 15 — QGEOM workstream findings

Agent: QGEOM (the quantum-mathematical thread: metric, curvature, trainability, CVaR).
Instrument reproduced at start, exactly:
`shipped 3.4540004952559396, pool_best 1.7108244199364904, top75_best 2.3061526409453816,
synthesis_fit 3.2040761603809194, n_zero_recall 18`.

Tiering: **DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **HYPOTHESIS** / **REFUTED**.

Everything is measured on the nine fully enumerated targets
(`s13/results/qarch_enum_<PDB>.npz`, n=9, k=4, 262,144 configurations, a true CA-RMSD on
every one) unless stated otherwise. Native quantities are read only for post-hoc scoring.

**Status: IN PROGRESS.** Written continuously; every module checkpoints to
`s15/results/qgeom_*.json` after every cell.

---

## 0. HEADLINE

**The paper's proposed thesis is refuted, and what replaces it is sharper.**

1. **CONDITIONING DOES RESHAPE THE VISITED GEOMETRY — AND ONLY ITS ENTROPY DOES.** Preparing
   a state for the retrieval-conditioned torsion prior raises the Fubini-Study condition
   number by **3.85x** (`block` L=2) to **36x** (`ring` L=2), cuts the metric volume by 1.8 to
   3.7 log units and the effective-rank FRACTION by 4.1 to 7.5 percentage points, at
   **W/L 27/0**. Against an
   **entropy-matched scramble of its own prior** — every per-residue entropy identical to the
   bit, all target information destroyed — it is **NULL on every one of seven statistics, on
   both ansatz, 54 paired cells**. The scramble reproduces the effect to within 3-7%. The one
   statistic whose CI excluded zero fails the null-calibrated concentration check
   (p = 0.000-0.001). **The chain `target-specific information -> manifold -> metric ->
   trainability` does not exist: the middle link is real and information-blind.**
2. **QNG'S ANSWER IS TWO-SIDED, AND THE BUDGET CONVENTION DECIDES ONE SIDE.** On a well-posed
   distribution-fitting cost QNG beats plain gradient descent **5 of 5** in the ill-conditioned
   cells at equal iterations — Sprint 14's depth-1 refutation genuinely does not extend. **At
   equal hardware cost the full geometric tensor wins 0 of 7 cells**, with two clean sign
   flips: reporting only the equal-iteration convention would have produced a false positive.
   And **on the cost this project actually optimises, QNG is statistically indistinguishable
   from plain gradient descent where the condition number is 688-2,677 (W/L 7/10, CI spanning
   zero) and loses to plain Adam 17 of 17.**
3. **THE METRIC A SAMPLING VQE LIVES ON IS `4 g` EXACTLY** (max abs diff **3.3e-16**, plus a
   proof). For a real-amplitude ansatz with a diagonal cost in the computational basis, QNG
   and classical natural gradient are the same preconditioner. There is no second metric to
   try, and no "wrong metric" escape from a negative result.
4. **THE GRADIENT LIVES IN THE RANGE OF THE METRIC (theorem, confirmed to 1e-30) AND AVOIDS
   ITS SMALL DIRECTIONS** — the bottom eigen-decile carries 0.0016-0.0043 of the gradient's
   squared norm against 0.10 for uniform. The rescaling QNG actually applies along the
   directions the gradient occupies is only **1.4-3.1x** even at condition number 2,000.
   **The ill-conditioning is real and lives where the gradient is not.**
5. **NO BARREN PLATEAU ANYWHERE.** `Var[dC/dtheta]` decays at **2^(-0.18 to -0.36 n)** with
   every bootstrap CI excluding the textbook -1.0; the retrieval prior decays at **+0.006**,
   i.e. not at all. What does cost trainability is **depth** (`L^-1` at 8 qubits, `L^-2.6` at
   16) and **cost locality** (6.05 log2 across the weight axis at n=14, against 2.8 log2
   across the whole 6-to-18-qubit width axis) — the axis Sprint 14 recorded as irrelevant.
6. **CVaR HAS NO DEFENSIBLE ROLE AS AN OBJECTIVE, A TAIL-SHAPER OR A CONSTRAINT** — including
   at the set level, and including as a risk measure on a clash term where a plain expectation
   penalty drives the violating mass to exactly zero and CVaR leaves 2.8-5.6%. Its one real
   effect is a **diversity dial**: 2.7 distinct structures in 2,048 draws at alpha=1 against
   408 at alpha=0.01. **Concentration buys the set MEAN (-0.85 A) and pays the set BEST
   (+1.09 A) in near-equal measure, and `alpha` scales both.**
7. **THE ONE POSITIVE VQE RESULT IN THIS PROJECT, hedged.** Consumed as an **ensemble with no
   ranker anywhere** — the project's actual situation — a CVaR-VQE beats best-of-N from its
   own untrained start by **0.36-0.57 A** on a coordinate-average readout, CIs excluding zero,
   at alpha=1 and at alpha=0.05 where the ensemble is genuinely diverse. The concentration
   check passes but with **LOW POWER** (mean/sd = -0.40); replicate before building on it.
8. **CONVERGENCE IS NOT A BASIN, WITH EXACT CURVATURE.** The gradient norm falls **313x** and
   the mode RMSD moves **0.049 A**. Only **25%** of runs reach a positive-semidefinite
   Hessian; the rest stop on saddle-adjacent plateaus with 27% of the gradient in
   negative-curvature directions. **The ansatz that optimises best (`block`, cost down 143x)
   reaches a minimum in 0 of 18 runs and improves the structure by 0.46 A; the ansatz that
   reaches a minimum (`ring`, PSD in 50%) makes the structure 0.37 A worse.**
9. **SPRINT 14'S METRIC TABLE IS AN AVERAGING ARTEFACT.** Its condition numbers and ranks are
   statistics of a **seed-averaged** metric, are not properties of any point, and do not
   reproduce from the recorded code at either width. Per point the true spread is far larger
   (`chain` L=3 median **735**, max **1.5e6**) and the `block` ansatz is **exactly
   rank-deficient** at depth >= 2, which the seed-average hides. `g_ii = 0.250000` survives.

---

## V. VERIFICATION — three independent checks before any claim

`s15/qgeom_lib.py` computes every derivative **analytically and exactly**, via the identity

    d/dt RY(t) = (1/2) RY(t + pi)

so `d|psi>/dtheta_j` is the circuit state with `theta_j` shifted by `pi`, halved. There is no
step size. Sprint 14's `fubini_study` uses central differences at `h = 1e-5`. The two are
mathematically independent routes to the same object.

| check | result |
|---|---|
| exact analytic metric vs Sprint 14's central differences, 7 patterns x 3 depths x 4 seeds | max abs diff **2.02e-11** (the `h^2` truncation) |
| **classical Fisher information of `p` vs `4 g`** | max abs diff **3.33e-16** |
| `uniform_theta` gives the exactly uniform distribution | max abs dev **8.67e-19** |
| exact gradient of `<E>` vs central differences | cosine **1.000000**, max rel err **3.25e-09** |
| exact gradient of `CVaR_0.25` vs central differences | cosine **1.000000**, max rel err **5.78e-09** |

### V1. `F = 4g` is a theorem in this setting — DEMONSTRATED (new)

RY + CNOT produces **real** amplitudes. Then `p(x) = psi(x)^2` and
`dp_j = 2 psi (d_j psi)`, so the classical Fisher information of the measured distribution is

    F_ij = sum_x (dp_i)(dp_j)/p = 4 sum_x (d_i psi)(d_j psi) = 4 (D D^T)

while normalisation gives `<d_j psi|psi> = 0` identically, so the Fubini-Study metric is

    g_ij = <d_i psi|d_j psi> - <d_i psi|psi><psi|d_j psi> = (D D^T)_ij .

Hence **`F = 4g` exactly**, which the table confirms to machine precision.

**Why it matters.** The literature distinguishes quantum natural gradient (QFI/Fubini-Study)
from classical natural gradient (Fisher information of the sampled outcomes) and treats the
choice as a design decision. **For a real-amplitude ansatz with a diagonal cost read out in
the computational basis — which is exactly this project's setting — they are the same
preconditioner up to a factor of four.** So there is no "we used the wrong metric" escape
hatch for a negative QNG result, and no second arm to try. The distinction only reappears
when amplitudes are complex (the Berry curvature term) or the observable is non-diagonal.

**Classification: DEMONSTRATED**, proof plus a 3.3e-16 numerical confirmation over 84 cells.

---

## A1. THE SPRINT 14 CONDITION NUMBERS ARE AN AVERAGING ARTEFACT — DEMONSTRATED (correction)

`s14/vqe_ansatz.metric_study` computes `g` at six different random `theta`, **averages the six
matrices**, and reports the condition number and rank of the average. `g` is a function of
`theta`; the average of six metrics at six different points of the manifold is not the metric
at any point, and neither its condition number nor its rank is a property of the ansatz.

Recomputed **per point**, 12 seeds, n=10, condition number over the non-null eigendirections:

| pattern | L | P | rank | cond median | cond min | cond max | s14-style cond(mean g) |
|---|---|---|---|---|---|---|---|
| none | 1 | 10 | 10 | 1.00 | 1.00 | 1.0 | 1.000 |
| none | 2 | 20 | **10** | 1.00 | 1.00 | 1.0 | **inf** |
| none | 3 | 30 | **10** | 1.00 | 1.00 | 1.0 | **inf** |
| chain | 1 | 10 | 10 | 1.00 | 1.00 | 1.0 | 1.000 |
| chain | 2 | 20 | 20 | **20.85** | 8.87 | 88,546 | 1.990 |
| chain | 3 | 30 | 30 | **735.40** | 22.98 | 1,528,176 | 2.167 |
| ring | 2 | 20 | 20 | 8.82 | 4.39 | 1,366 | 1.350 |
| ring | 3 | 30 | 30 | 10.25 | 4.61 | 1,861 | 1.491 |
| brick / block_chain | 2 | 20 | 20 | **202.39** | 5.36 | 1,091,653 | 3.806 |
| brick / block_chain | 3 | 30 | 30 | **634.10** | 33.76 | 7,474,642 | 4.564 |
| **block** | 2 | 20 | **15** | 16.24 | 3.87 | 1,598 | 15.301 |
| **block** | 3 | 30 | **15** | 6.04 | 2.68 | 9.7 | 20.174 |
| all_to_all | 2 | 20 | 20 | 12.46 | 3.56 | 81.2 | 3.769 |
| all_to_all | 3 | 30 | 30 | 40.72 | 8.34 | 278.5 | 7.701 |

**Three corrections.**

1. **The recorded numbers do not reproduce from the recorded code.** Running
   `s14.vqe_ansatz.metric_study`'s own function with its own seeds gives `block` L=2 = **39.3**
   at n=8 and **39.3** at n=10, and `block` L=3 = 59.7 at n=8 / 446.3 at n=10 — against the
   recorded **155.2** and **478.9**. Neither n reproduces the recorded pair. `chain` L=2 gives
   2.04 against the recorded 3.665, `ring` L=2 gives 3.18/2.04 against the recorded 9.937, and
   `all_to_all` L=2 gives 3.5/3.8 against the recorded 2.478. *The specific values in the
   Sprint 14 table should not be cited.* The qualitative claim they support — the metric leaves
   `I/4` at depth ≥ 2 and becomes ill-conditioned — is **confirmed and is in fact much
   stronger** than recorded.
2. **The true conditioning is far worse and far more variable than either number.** Per point,
   `chain` at depth 3 has a median condition number of **735** and reaches **1.5e6**;
   `brick`/`block_chain` reach **7.5e6**. Averaging across seeds destroys this by mixing
   independent random matrices, which drives the average toward its own mean and *lowers* the
   condition number.
3. **The rank column is the artefact in reverse.** The seed-average of six rank-15 matrices
   with different null spaces has rank 20, so `s14` records the `block` ansatz as full rank at
   depth 2 when **it is exactly rank-deficient at every point**.

`g_ii = 0.250000` exactly at every pattern and depth is **confirmed** — that part of the
Sprint 14 record is universal and survives.

Reproduce: `python -m s15.qgeom_metric`; the n=8 side-by-side is
`s15/results/qgeom_metric_n8_repro.json`.

### A1b. The rank of a block-factorised ansatz is capped by the block — DEMONSTRATED (new)

A `block` ansatz's state is a **product over blocks**, so its reachable manifold has at most
`2^b - 1` real dimensions per block regardless of how many parameters are spent on it.
Predicted rank `min(P, n_blocks * (2^b - 1))`, measured over 6 seeds:

| pattern | b | L | P | rank | predicted | match |
|---|---|---|---|---|---|---|
| none | 2 | 2 | 24 | 12 | 12 | yes |
| block | 2 | 2 | 24 | 18 | 18 | yes |
| block | 2 | 3 | 36 | 18 | 18 | yes |
| block | 2 | 4 | 48 | 18 | 18 | yes |
| block | 3 | 2 | 24 | 24 | 24 | yes |
| block | 3 | 3 | 36 | 28 | 28 | yes |
| block | 3 | 4 | 48 | 28 | 28 | yes |

**16 of 16 cells match the prediction exactly.** The product ansatz's infinite condition
number that Sprint 14 identified is the `b = 1` case of the same law. Any "torsion-aware"
block entangler saturates at `n_res * (2^{log2 k} - 1) = n_res * (k-1)` parameters — which is
**exactly the number of free parameters of a per-residue categorical distribution over k
states**, so the block ansatz is *precisely* expressive enough for a product torsion prior and
not one parameter more. Depth beyond 2 (at k=4) buys nothing at all.

---

## A2. THE GRADIENT LIES IN THE RANGE OF THE METRIC, AND AVOIDS ITS SMALL DIRECTIONS — DEMONSTRATED (new)

**Theorem.** For any cost depending on `theta` only through the state,
`grad_j C = <d_j psi | dC/dpsi> = (D u)_j`, so `grad C in range(D)`. And
`g = D D^T` (using `D psi = 0`), so `range(g) = range(D)`. Therefore

    grad C  lies in  range(g)   ALWAYS,   for every cost, every ansatz, every point.

A rank-deficient metric's infinite condition number is therefore **free**: the flat directions
are exactly the directions the gradient never has a component along.

Measured, n=10, 8 seeds, exact gradient of a native-free conditioned Legacy objective:

| pattern | L | rank | grad share in null space | **bottom eigen-decile** | top eigen-decile | cond\|range |
|---|---|---|---|---|---|---|
| none | 2 | 10 | 6.3e-31 | 0.2119 | 0.0885 | 1.00 |
| chain | 2 | 20 | 0.0 | **0.0175** | 0.1640 | 86.8 |
| chain | 3 | 30 | 0.0 | **0.0043** | 0.2954 | 1652.7 |
| ring | 2 | 20 | 0.0 | 0.0324 | 0.0712 | 9.6 |
| ring | 3 | 30 | 0.0 | **0.0192** | 0.1199 | 12.1 |
| brick / block_chain | 2 | 20 | 0.0 | **0.0040** | 0.1122 | 225.9 |
| brick / block_chain | 3 | 30 | 0.0 | **0.0020** | 0.2029 | 634.1 |
| block | 2 | 15 | 2.1e-30 | **0.0016** | 0.0064 | 24.3 |
| all_to_all | 2 | 20 | 0.0 | **0.0043** | 0.0833 | 12.5 |
| all_to_all | 3 | 30 | 0.0 | **0.0018** | 0.2761 | 31.1 |
| *(depth 1, any pattern)* | 1 | 10 | 0.0 | 0.09–0.23 | 0.07–0.18 | 1.00 |

**Readings.**

1. **The theorem holds numerically**: the gradient's share in the metric's exact null space is
   `0.0` or `~1e-30` in every cell. Sprint 14's HYPOTHESIS is now proved for the exactly-null
   directions.
2. **Where the metric is well conditioned (depth 1, `g = I/4`) the gradient is spread evenly**
   across the spectrum (0.09–0.23 per decile against 0.10 for uniform) — as it must be, since
   there is no spectrum to align with.
3. **Where the metric is ill-conditioned the gradient systematically AVOIDS the small
   eigendirections.** The bottom decile carries 0.0016–0.0043 of the gradient's squared norm
   in every ill-conditioned cell — a 25x to 60x depletion against uniform — while the top
   decile is enriched to 0.11–0.30.
4. **This predicts that QNG must hurt, and says why.** The natural gradient divides each
   component by its eigenvalue. In the bottom decile that multiplies a component carrying
   ~0.2% of the gradient by up to `cond ~ 1e3`–`1e6`. Those components are the ones most
   contaminated by shot noise and by finite-difference error, and they are the least
   informative. **Preconditioning by `g` amplifies precisely the directions the exact gradient
   has already found to be useless.** A pre-registered prediction, recorded here before A3/A4
   were read: **QNG will be neutral at best and will need heavy Tikhonov regularisation
   (i.e. to be turned back into ordinary gradient descent) to avoid being worse.**

---



---

## B. DOES TARGET-SPECIFIC CONDITIONING RESHAPE THE MANIFOLD? — NO. Only the ENTROPY of the conditioning is visible to the geometry.

`s15/qgeom_cond.py`. **Pre-registered before the run** (the hypotheses are in the module
docstring, written first):

* **H0** the conditioned point is geometrically indistinguishable from a random point.
* **H1** the geometry tracks the **entropy** of whatever distribution the state was prepared
  for, and not its target-specificity.
* **H2** target-specific conditioning is geometrically special — the paper's thesis.

**The discriminating control is an entropy-matched scramble**: the target's own retrieval
posterior with the `k` state labels permuted independently within each residue. Every
per-residue entropy is *identical to the bit*; only the assignment of probability to states —
i.e. all of the target information — is destroyed.

Design: 9 enumerated targets x 3 seeds x 7 arms, 12-qubit sub-register (residues 1-6),
torsion-aware `block` ansatz at depth 2 — which A1b proves is *exactly* expressive for a
product torsion prior (`n_res (k-1)` free parameters, no more and no less). Each arm's
`theta` is obtained by minimising `KL(target || p_theta)` with the identical optimiser,
budget and initialisation seed. The retrieval prior is
`s14/retprior.state_prior(pdb, "top75", k=4)`, the best channel in the project
(phi 33.6 / psi 59.2 deg).

**Fit-quality control — the confound that would invalidate everything.**
`KL(COND) - KL(scram_state) = +0.00458 [-0.02117, +0.03542]`, W/L 11/16, **NULL**, median
`1.5e-12`. Both arms fit essentially exactly, as A1b's rank law predicts. The comparison is
clean.

### B2. The paired table, n = 27 cells — DEMONSTRATED

Negative = the first arm is lower. `cond` is compared on `log10`, so its effect is a ratio.

| statistic | COND - uniform | COND - **scram_state** | scram_state - uniform |
|---|---|---|---|
| **log10 condition number** | **+0.585** [+0.508,+0.663] 0/27 **SIG** | +0.013 [-0.093,+0.119] 16/11 **NULL** | **+0.573** [+0.497,+0.656] 0/27 **SIG** |
| **effective rank fraction** | **-0.0408** [-0.0475,-0.0345] 27/0 **SIG** | -0.0020 [-0.0119,+0.0065] 13/14 **NULL** | **-0.0388** [-0.0454,-0.0321] 27/0 **SIG** |
| **log metric volume** | **-1.803** [-2.162,-1.473] 27/0 **SIG** | -0.092 [-0.597,+0.359] 12/15 **NULL** | **-1.711** [-2.058,-1.371] 27/0 **SIG** |
| **max abs(g - I/4)** | **-0.0037** [-0.0050,-0.0026] 27/0 **SIG** | -0.0010 [-0.0024,+0.0003] 16/11 **NULL** | **-0.0027** [-0.0037,-0.0019] 27/0 **SIG** |
| gradient variance | +0.0000 [-0.0003,+0.0003] 12/15 NULL | -0.0005 [-0.0010,-0.0001] 17/10 SIG | +0.0005 [+0.0002,+0.0009] 9/18 SIG |
| **state entropy (bits)** | **-4.621** [-5.605,-3.688] 27/0 **SIG** | +0.011 [-0.034,+0.062] 14/13 **NULL** | **-4.632** [-5.610,-3.697] 27/0 **SIG** |
| grad share, bottom eigen-decile | **-0.0234** [-0.0498,-0.0029] 17/10 **SIG** | -0.0033 [-0.0209,+0.0148] 11/16 **NULL** | -0.0202 [-0.0559,+0.0116] 16/11 NULL |

**Three readings, and they settle the question.**

1. **H0 is REFUTED.** Conditioning genuinely moves the state to a geometrically different
   region: the condition number rises by a factor of **3.85** (`10^0.585`), the effective
   rank falls 4.1 points of a fraction, the metric volume falls **1.80 log units** (a factor
   of 6.1) — all at W/L 27/0 with CIs nowhere near zero. **The geometry at the point you
   actually visit is not a constant, and preparing an informative state is a move to a
   worse-conditioned, lower-volume, lower-effective-rank neighbourhood.**
2. **H2 is REFUTED, on every statistic.** Against the entropy-matched scramble of its *own*
   prior, the target-conditioned point is **NULL on all seven statistics**, W/L within one or
   two of even, medians at or below `1e-4`. The single exception — gradient variance at
   -0.0005 — is contradicted by the `scram_state - uniform` row on the same axis (+0.0005,
   opposite sign) and is unsupported by any other statistic.
3. **H1 is CONFIRMED, quantitatively.** `scram_state - uniform` reproduces `COND - uniform`
   almost exactly on every axis: +0.573 vs +0.585 on log condition number, -1.711 vs -1.803
   on log volume, -4.632 vs -4.621 on entropy, -0.0388 vs -0.0408 on effective rank.
   **The whole effect is carried by the entropy of the prepared distribution; none of it is
   carried by whether that distribution was the right one for the target.**

**The ordering across controls is itself evidence for H1.** `scram_state` and `scram_res`
match the prior's per-residue entropies *exactly* and are null; `cross` (a different real
target's prior) is null; `dirichlet`, which matches only the **total** entropy and not the
per-residue profile, is the one control that separates (`COND - dirichlet` = -0.259 on log
cond, SIG). **The closer the entropy match, the closer the geometry** — which is what H1 says.

### B2b. REPLICATED on a structurally different ansatz — DEMONSTRATED

The whole of B2 was re-run on `ring` L=2 — **full rank 24/24**, where `block` L=2 is exactly
rank-deficient 18/24 — with the same 9 targets x 3 seeds and the same seven statistics.

| statistic | COND - uniform | COND - **scram_state** | scram_state - uniform |
|---|---|---|---|
| log10 condition number | **+1.555** [+1.142,+2.001] 0/27 SIG | -0.042 [-0.440,+0.360] 13/14 **NULL** | **+1.597** [+1.281,+1.913] 0/27 SIG |
| effective rank fraction | **-0.0753** [-0.0896,-0.0604] 27/0 SIG | +0.0087 [-0.0065,+0.0237] 11/16 **NULL** | **-0.0840** [-0.0993,-0.0686] 27/0 SIG |
| log metric volume | **-3.731** [-4.858,-2.636] 27/0 SIG | +0.270 [-0.672,+1.204] 11/16 **NULL** | **-4.000** [-4.979,-3.039] 27/0 SIG |
| max abs(g - I/4) | **+0.1865** [+0.1611,+0.2087] 0/27 SIG | -0.0107 [-0.0354,+0.0119] 15/12 **NULL** | **+0.1972** [+0.1746,+0.2174] 0/27 SIG |
| gradient variance | -0.0000 [-0.0002,+0.0001] 17/10 NULL | -0.0001 [-0.0003,+0.0001] 15/12 **NULL** | +0.0000 [-0.0001,+0.0002] 15/12 NULL |
| state entropy (bits) | **-3.088** [-3.883,-2.387] 27/0 SIG | +0.100 [-0.202,+0.405] 13/14 **NULL** | (matched by construction) |

**Identical structure, larger effect.** `COND - scram_state` is NULL on every statistic
again, and `scram_state - uniform` reproduces `COND - uniform` to within 3-7% on every axis
(+1.597 vs +1.555; -0.0840 vs -0.0753; -4.000 vs -3.731; +0.1972 vs +0.1865). The
conditioning effect is **larger** on the full-rank ansatz — the condition number rises **36x**
(`10^1.555`) against `block`'s 3.85x — so the result is not an artefact of the rank-deficient
torsion-aware ansatz; if anything `block` understates it.

**Two ansatz, 54 paired cells, seven statistics, one conclusion: conditioning moves the
geometry, and only its entropy does.**

### B3. The positive form: after entropy, the arm label carries nothing — DEMONSTRATED

`s15/qgeom_condlaw.py`, 189 points = 27 cells x 7 arms. `R2 arm | H resid` is the variance a
one-way ANOVA on the ARM LABEL still explains after the prepared state's entropy has been
regressed out linearly. The last column restricts to the **five fitted arms** (COND,
scram_state, scram_res, cross, dirichlet) — all produced by the identical optimiser on a
product prior, differing **only in whose prior it was**. For 5 groups and 135 points the
chance expectation of that R2 is `(5-1)/(135-1) = 0.030`.

| statistic | rho(stat, H) | R2 linear in H | R2 arm label | R2 arm \| H resid | **R2 arm \| H, FITTED ARMS ONLY** |
|---|---|---|---|---|---|
| log10 condition number | -0.422 | 0.192 | 0.488 | 0.258 | **0.099** |
| effective rank fraction | +0.435 | 0.227 | 0.459 | 0.205 | **0.008** |
| log metric volume | +0.438 | 0.180 | 0.370 | 0.172 | **0.032** |
| max abs(g - I/4) | +0.521 | 0.162 | 0.193 | 0.130 | **0.052** |
| gradient variance | +0.184 | 0.003 | 0.072 | 0.085 | **0.045** |
| grad share, bottom decile | +0.256 | 0.034 | 0.082 | 0.036 | **0.038** |
| max p | **-0.989** | 0.791 | 0.141 | 0.356 | **0.053** |

**Six of the seven statistics sit at or within 2x of the 0.030 chance level**, and the
seventh (condition number, 0.099) is 3.3x chance — which is exactly the statistic B2's paired
test also finds NULL (`COND - scram_state = +0.013 [-0.093, +0.119]`, W/L 16/11). *Which of
five priors the state was prepared for explains essentially none of its geometry once you
know how concentrated the prepared distribution was.*

**One caution against over-claiming.** The relationship to entropy is **monotone but not
linear**: `R2 linear in H` is only 0.16-0.23 for the metric statistics, while `rho` is
+0.42 to +0.52; the best single predictor of the geometry in the panel is `max p`
(`rho = -0.989` against H, `R2 = 0.791`), i.e. **concentration**, of which entropy is one
imperfect summary. So the correct statement is *"the geometry is a monotone function of how
concentrated the prepared distribution is"*, not *"the geometry is a linear function of its
entropy"*. The decisive evidence remains B2's entropy-matched paired test, which is
non-parametric and does not depend on the functional form at all.

### B1. How this sits beside the Sprint 14 result

Sprint 14: the Fubini-Study metric is bit-identical across energy models at matched
parameters (0.000e+00), so **the energy model does not reshape the manifold**.
`qgeom_metric.V` upgrades that from a measurement to a structural fact — `g` is built from
the state and its derivatives alone, and no energy array can enter it.

Sprint 15 supplies the other half. **Conditioning does move you to a different geometry, but
the map from information to geometry factors entirely through one scalar: the entropy of the
distribution you prepared.** Two states prepared to the same entropy have the same geometry
whether one is the target's true posterior and the other a scrambled fake.

**So the chain `target-specific information -> state manifold -> metric/curvature ->
trainability` does not exist as a causal chain. The middle link is real, and it is
information-blind: `entropy -> metric`, with target identity factored out.** Anything hoping
to use manifold geometry as a channel for target-specific information is looking where,
measurably, that information is not.

**Classification: H0 REFUTED, H1 DEMONSTRATED, H2 REFUTED** — 27 paired cells, seven
statistics, fit-quality confound checked and null.

---

## C. GRADIENT VARIANCE AND SCALING — IT IS NOT A BARREN PLATEAU (part 1 of 2; the depth, alpha, locality and energy-model axes are in "C (continued)" below)

`s15/qgeom_grad.py`. Exact gradients by **adjoint (reverse-mode) differentiation**, verified
against the `O(P)` derivative-state route to **2.43e-15**. Haar-random angles (the textbook
barren-plateau protocol). Objectives rank-conditioned so their variance is `1/12` at every
width, removing the scale confound. 48 seeds per cell; `Var[g_1]`, the textbook statistic.

| pattern | L | n=6 | 8 | 10 | 12 | 14 | 16 | 18 | **slope log2/qubit** | CI95 | excludes -1? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ring | 1 | 3.68e-3 | 1.53e-3 | 6.88e-4 | 7.22e-4 | 6.77e-4 | 6.34e-4 | 2.99e-4 | **-0.240** | [-0.348, -0.095] | **yes** |
| ring | 2 | 2.55e-3 | 6.98e-4 | 7.51e-4 | 5.49e-4 | 3.45e-4 | 4.08e-4 | 2.48e-4 | **-0.228** | [-0.353, -0.117] | **yes** |
| ring | 3 | 1.91e-3 | 4.62e-4 | 3.34e-4 | 3.88e-4 | 2.02e-4 | 2.57e-4 | 4.16e-5 | **-0.339** | [-0.457, -0.103] | **yes** |
| block | 1 | 5.38e-3 | 6.22e-4 | 4.09e-4 | 5.75e-4 | 4.64e-4 | 4.64e-4 | 2.67e-4 | **-0.244** | [-0.459, -0.014] | **yes** |

**The measured decay is `2^(-0.23 to -0.34 n)` and every bootstrap CI excludes the textbook
`-1.0`.** Over 6 to 18 qubits the variance falls by a factor of **12**; a barren plateau
demands a factor of **4,096**. **This is not a barren plateau and must not be called one.**
The brief's rule is met by evidence, not assertion: the scaling was measured over a 12-qubit
range with CIs, and it rejects the barren-plateau exponent.

Consistent with Sprint 14: the standard theorem needs local 2-design blocks this RY+CNOT
ansatz does not form, so it licenses no prediction here; and rank-conditioning removes the
delta-spike artefact that inflates a raw molecular energy's spectrum.

*The depth, alpha, locality and energy-model axes are in the C-continued section below.*

---

## D. CVaR, EXTENDED — THE CONCENTRATION IT BUYS IS DESTRUCTIVE, AND `alpha` IS A DIVERSITY DIAL

`s15/qgeom_cvar.py`. Nothing from Sprint 14's correctness audit is re-verified. Every arm is
given the **exact statevector gradient** (no shot noise, no SPSA) — deliberately the most
favourable setting a CVaR-VQE can have — and the budget is charged only for the draws taken
for the readout, with the control getting the same draws. 12-qubit sub-registers of the nine
enumerated targets, `hamil` w=0.25, 200 iterations, 3 seeds, 2,048 draws.

### D2. The mandatory control, at the SET level — DEMONSTRATED

Paired over 27 cells per alpha. Negative = VQE better than **best-of-N from its own untrained
initial distribution at the same budget**.

| readout | alpha=1.0 | alpha=0.25 | alpha=0.05 | alpha=0.01 |
|---|---|---|---|---|
| **argmin (by objective)** | -0.036 NULL | **+0.251 SIG** | +0.068 NULL | +0.048 NULL |
| top20 mean | -0.177 NULL | +0.111 NULL | -0.070 NULL | -0.091 NULL |
| **top20 COORDINATE AVERAGE** | +0.068 NULL | **+0.372 SIG** | +0.197 NULL | +0.177 NULL |
| **top75 COORDINATE AVERAGE** | +0.147 NULL | **+0.392 SIG** | **+0.276 SIG** | +0.164 NULL |
| **drawn-set MEAN** | **-0.847 SIG** [22/5] | **-0.372 SIG** | **-0.191 SIG** | -0.066 NULL |
| **drawn-set BEST** | **+1.092 SIG** [1/25] | **+0.345 SIG** | **+0.155 SIG** | +0.087 NULL |

**D2a. THE TRADE, in one pair of rows.** At alpha=1 the VQE's drawn set has a mean RMSD
**0.847 A better** than the control's and a best member **1.092 A worse** — two significant
effects of opposite sign and near-equal size, both shrinking monotonically to zero as alpha
falls. **That is what concentration is, measured directly: it buys the mean and pays the
best, and `alpha` scales both.**

**D2b. A TRAP THE PROJECT SHOULD RECORD.** The recorded law is that the terminal operator
consumes the **set mean** (`d_out = 1.16 d_set_mean + 0.04 d_set_best`). Read naively, the
`drawn-set MEAN` row says CVaR-VQE wins by 0.85 A. **The actual coordinate average says it
loses.** The resolution is the diversity column:

| alpha | distinct configurations among 2,048 VQE draws | control | final entropy (bits) | max p |
|---|---|---|---|---|
| 1.00 | **2.7** (median **2**) | 510.3 | 0.73 | 0.702 |
| 0.25 | 149.4 | 510.3 | 3.93 | 0.319 |
| 0.05 | 314.9 | 510.3 | 5.79 | 0.155 |
| 0.01 | 408.0 | 510.3 | 6.58 | 0.101 |

**At alpha=1 the VQE returns two distinct structures out of 2,048 draws.** Its "top-20 set"
is twenty copies of one structure, so its set *mean* is that structure's RMSD and there is
nothing to average. The set-mean law was fitted on sets with real diversity and is **out of
domain on a collapsed set**. *Set mean is a valid proxy for the terminal operator only at
fixed set diversity; a method that concentrates improves the proxy and degrades the
operator.* This is the project's recurring "better matrix, worse ranking" shape on a new axis.

### D3. The regime sweep — no regime, and the same artefact twice — DEMONSTRATED

Objective quality (ORACLE blend knob) x alpha x readout, 6 targets x 2 seeds.
Cells are `VQE / control (W/L)`; W = VQE better.

| readout | signal | alpha=0.05 | alpha=1.0 |
|---|---|---|---|
| top20 mean | 0.1 | 2.893 / 3.252 (**9/1**) | 3.515 / 2.970 (3/7) |
| top20 mean | 0.3 | 2.631 / 3.085 (**11/0**) | 3.071 / 3.085 (6/5) |
| **top20 COORDAVG** | 0.1 | 2.893 / 2.970 (7/3) | 3.515 / 2.970 (1/9) |
| **top20 COORDAVG** | 0.3 | 2.631 / 2.680 (6/5) | 3.054 / 2.680 (2/9) |
| **top75 COORDAVG** | 0.7 | 2.394 / **2.097** (**1/10**) | 2.709 / 2.097 (2/9) |
| **top75 COORDAVG** | 0.0 | 4.256 / **3.696** (**2/10**) | 4.132 / 3.696 (4/8) |

Small-alpha CVaR looks like a clear winner on top-m **mean** at every signal >= 0.1 (W/L to
11/0) — **and the same cells on the coordinate average are a tie or a loss**. The apparent win
is the collapse artefact of D2b, seen again on an independent axis. **On the readout the
pipeline actually uses, there is no signal level and no alpha at which CVaR-VQE beats
best-of-N from its own start.**

### D4. CVaR as a CONSTRAINT — the last defensible role, and it fails too — DEMONSTRATED

`minimise E[distogram] + mu * RISK(steric clash)`, with `RISK` either an expectation (the
incumbent) or a CVaR of the clash term's **upper** tail — the use a risk measure was invented
for, and one Sprint 14 never tested. 6 targets x 3 seeds.

| mu | risk | violating mass | argmin | top20 coordavg | vs control (coordavg) |
|---|---|---|---|---|---|
| 0.5 | **expectation** | **0.0000** | 3.028 | 3.045 | +0.486 [+0.170,+0.799] 3/15 SIG |
| 0.5 | CVaR 0.10 | 0.0278 | 3.195 | 3.168 | +0.609 [+0.339,+0.903] 3/15 SIG |
| 0.5 | CVaR 0.02 | 0.0278 | 2.996 | 2.967 | +0.407 [+0.093,+0.737] 5/13 SIG |
| 2.0 | **expectation** | **0.0000** | 3.116 | 3.116 | +0.557 [+0.229,+0.834] 3/15 SIG |
| 2.0 | CVaR 0.10 | 0.0556 | 3.299 | 3.259 | +0.700 [+0.319,+1.071] 3/15 SIG |
| 2.0 | CVaR 0.02 | 0.0278 | 3.102 | 3.067 | +0.508 [+0.097,+0.947] 5/13 SIG |

**A plain expectation penalty drives the violating mass to exactly zero; the CVaR risk
measure leaves 2.8-5.6% of the mass in violation and returns worse structures.** CVaR is
worse than an expectation at the one job risk measures exist for — because once the tail is
clean the CVaR term stops caring, while an expectation keeps pushing. **And every arm, CVaR
or not, loses to best-of-N from its own initial distribution by +0.41 to +0.70 A with CIs
excluding zero.**

### D. Verdict (extended by E1/E2 below)

**CVaR has no role as an objective, no role as a tail-shaper, and no role as a constraint
device in this problem.** The one thing it demonstrably *is* is a **diversity dial with an
exactly known mechanism** — Sprint 14's own iff condition `dCVaR/dp == 0 <=> p(argmin) >=
alpha`, read as a **stopping rule on concentration**. That mechanism, and the last remaining
steelman — the VQE consumed as an **ensemble generator with no ranker anywhere**, which is
the project's actual situation given that nothing ranks within a pool — are tested in
`s15/qgeom_ens.py`.

---

## C (continued). THE FULL SCALING PICTURE — DEPTH AND COST LOCALITY, NOT WIDTH

`s15/qgeom_grad.py`, `s15/qgeom_loc.py`. All numbers `Var[g_1]` over 48-64 Haar-random
initialisations, exact adjoint gradients, rank-conditioned objectives.

### C1. Width — DEMONSTRATED

Nine pattern x depth cells, n = 6..18. **Every slope is between -0.18 and -0.36 log2 per
qubit and every bootstrap CI excludes the textbook -1.0.**

| pattern | L=1 | L=2 | L=3 |
|---|---|---|---|
| ring | -0.240 [-0.348,-0.095] | -0.228 [-0.353,-0.117] | -0.339 [-0.457,-0.103] |
| block | -0.244 [-0.459,-0.014] | -0.181 [-0.394,+0.073] | -0.205 [-0.402,+0.095] |
| all_to_all | -0.363 [-0.544,-0.008] | -0.204 [-0.352,-0.024] | -0.175 [-0.315,+0.091] |

Over 6 -> 18 qubits the gradient variance falls by a factor of ~12. A barren plateau demands
a factor of 4,096. **This is not a barren plateau.**

### C2. DEPTH IS THE DOMINANT AXIS, AND ITS PENALTY GROWS WITH WIDTH — DEMONSTRATED (new)

`ring`, L = 1, 2, 3, 4, 6, 8 at fixed width:

| n | L=1 | L=2 | L=3 | L=4 | L=6 | L=8 | **slope in log2 L** |
|---|---|---|---|---|---|---|---|
| 8 | 1.53e-3 | 6.98e-4 | 4.62e-4 | 2.95e-4 | 2.56e-4 | 2.02e-4 | **-0.98** |
| 12 | 7.22e-4 | 5.49e-4 | 3.88e-4 | 1.15e-4 | 8.10e-5 | 2.78e-5 | **-1.57** |
| 16 | 6.34e-4 | 4.08e-4 | 2.57e-4 | 5.25e-5 | 1.50e-5 | 2.76e-6 | **-2.61** |

**`Var ~ L^-1` at 8 qubits and `L^-2.6` at 16 qubits.** Going from L=1 to L=8 at n=16 costs a
factor of **230**; going from 6 to 18 qubits at fixed depth costs a factor of 12. **Depth,
not width, is what removes the gradient in this ansatz, and the depth exponent itself grows
with width** — the two axes interact, which a pure `2^-n` framing cannot express. Practical
consequence: the shipped depth-2/3 operating point is the right side of this curve, and the
usual "add layers for expressivity" move is the expensive one here.

### C3. CVaR alpha lowers the LEVEL and leaves the EXPONENT alone — DEMONSTRATED

`ring` L=2, n = 8, 10, 12, 14:

| alpha | n=8 | n=10 | n=12 | n=14 | slope/qubit |
|---|---|---|---|---|---|
| 1.00 | 6.98e-4 | 7.51e-4 | 5.49e-4 | 3.45e-4 | -0.175 [-0.334,+0.053] |
| 0.25 | 7.03e-4 | 6.82e-4 | 5.34e-4 | 3.53e-4 | -0.167 [-0.299,-0.022] |
| 0.10 | 5.00e-4 | 3.66e-4 | 1.96e-4 | 1.40e-4 | -0.321 [-0.451,-0.225] |
| 0.05 | 2.83e-4 | 1.55e-4 | 1.26e-4 | 8.09e-5 | -0.286 [-0.437,-0.146] |
| 0.01 | 9.28e-5 | 3.51e-5 | 2.91e-5 | 3.00e-5 | -0.258 [-0.702,+0.022] |

Small alpha costs about **an order of magnitude of gradient variance** (7.0e-4 -> 9.3e-5 at
n=8) but the width exponent is statistically unchanged. **A shrinking tail is a smaller
effective sample, not a plateau**, and the distinction matters: the fix for a small effective
sample is more shots, and the fix for a plateau is a different ansatz.

### C4. COST LOCALITY — MY FIRST MEASUREMENT WAS CONFOUNDED; CORRECTED IT IS LARGE — DEMONSTRATED

**Retraction of my own cell.** `qgeom_grad.locality` built each objective from
`terms = 32` random Pauli-Z strings of exact weight `w`. At `w=1` only `n` supports exist and
at `w=n` only **one**, so those columns are built from far fewer distinct terms than the
middle of the range. The resulting profile is a deep U (4.9e-1 at w=1, 9.9e-3 at w=6,
2.2e+0 at w=8 for n=8) and the slope fitted through it is meaningless. **That is a defect in
my own first measurement.** `s15/qgeom_loc.py` fixes it: every weight uses exactly 28
*distinct* strings, variance is 1 by construction at every weight, and only weights with
`C(n,w) >= 28` enter the trend.

| n | trend weights | log2 range over the trend | slope per unit weight |
|---|---|---|---|
| 8 | 2..6 | 1.72 | -0.341 |
| 10 | 2..8 | 3.07 | -0.487 |
| 12 | 2..10 | 4.49 | -0.581 |
| 14 | 2..12 | **6.05** | -0.426 |

**Corrected verdict, and it disagrees with the Sprint 14 record.** At n=14 the gradient
variance falls **6.05 log2 units (a factor of 66) across the cost-locality axis**, against
**2.8 log2 (a factor of 7)** across the entire 6-to-18-qubit width axis measured on the same
machinery. **Cost locality is the LARGER of the two axes here, not the smaller.** The
per-unit-weight slope is roughly constant at **-0.45**, so a maximally global cost would decay
at `2^(-0.45 n)` — twice the exponent of a real structural objective, and still not `-1`.

**Scope, carefully.** Sprint 14's statement ("the measured ansatz kernel is flat in Pauli
weight; what matters is n") is about a different quantity — the ansatz's response kernel, not
the gradient variance of a weight-controlled cost family — and the two need not agree. What
is refuted is the *operational* reading, that cost locality can be ignored at 6-18 qubits in
favour of width. **It cannot: it is the bigger lever.**

### C5. The energy model changes the EXPONENT, and it does so through locality — DEMONSTRATED

`ring` L=2, n = 8, 10, 12, 14, all objectives rank-conditioned to the identical marginal:

| objective | n=8 | n=10 | n=12 | n=14 | slope/qubit |
|---|---|---|---|---|---|
| random (unstructured, maximally global) | 8.21e-4 | 4.56e-4 | 1.41e-4 | 1.14e-4 | **-0.511** [-0.845,-0.152] |
| Legacy | 6.98e-4 | 7.51e-4 | 5.49e-4 | 3.45e-4 | -0.175 [-0.334,+0.053] |
| **retrieval torsion prior (1-local)** | 1.25e-3 | 1.59e-3 | 1.49e-3 | 1.31e-3 | **+0.006** [-0.096,+0.177] |
| true RMSD (ORACLE DIAGNOSTIC) | 2.47e-3 | 5.25e-4 | 4.37e-4 | 8.90e-4 | -0.234 [-1.116,+0.514] |

**The 1-local retrieval prior has NO measurable decay with width at all** (slope +0.006, CI
straddling zero), while an unstructured random objective decays at -0.511. This is C4's
locality law seen on real objectives and is an independent confirmation of it: *the width
scaling of the gradient is set by how local the cost is, and the project's best information
channel is the most local objective it has.*

### C. Verdict

**No barren plateau at any width, depth, alpha or energy model tested.** The trainability
budget at this scale is spent on **depth** (up to `L^-2.6`) and **cost locality**
(up to `2^-0.45n`), not on width (`2^-0.23n` for a real objective, `2^0.0n` for the retrieval
prior). Calling any of this a barren plateau would be a category error, and the brief's
warning is met by measurement.

---

## A5. WHY QNG HELPS ON ONE COST AND NOT ANOTHER — MY MECHANISM IS REFUTED TOO

`s15/qgeom_align.py`. A2 measured the gradient's spectral profile for the **expectation
value** and I pre-registered from it that **QNG would be neutral at best**. A3 (below) finds
QNG clearly helping on the KL fitting task in exactly the ill-conditioned cells, so the
prediction is wrong there, and I proposed a mechanism: the KL gradient carries a `t(x)/p(x)`
factor that blows up on low-probability states, which is where the small eigendirections
live, so KL should load the bottom of the spectrum where the expectation does not.

**Measured on identical circuits and identical points, that mechanism is REFUTED.**

| pattern | L | cond | expectation bot10 / amplif | CVaR 0.25 bot10 / amplif | KL bot10 / amplif |
|---|---|---|---|---|---|
| ring | 1 | 1.00 | 0.170 / 1.0 | 0.071 / 1.0 | 0.093 / 1.0 |
| all_to_all | 2 | 7.5 | 0.0212 / 1.4 | 0.0204 / 1.6 | 0.0160 / 3.0 |
| ring | 2 | 9.6 | 0.0425 / 2.1 | 0.0228 / 2.1 | 0.0277 / 1.9 |
| block | 2 | 24.3 | 0.0229 / 1.4 | 0.0194 / 1.4 | **0.0002** / 1.9 |
| chain | 2 | 85.9 | 0.0498 / 2.3 | 0.0344 / 2.1 | 0.0239 / 2.1 |
| brick | 2 | 148.9 | 0.0529 / 2.1 | 0.0381 / 2.3 | **0.0025** / 2.0 |
| chain | 3 | 2003.7 | 0.0053 / 2.6 | 0.0073 / 2.6 | 0.0058 / 2.6 |
| brick | 3 | 1477.2 | 0.0060 / 3.1 | 0.0032 / 2.6 | 0.0016 / 2.6 |

The KL gradient loads the bottom eigen-decile **less**, not more, than the expectation
gradient (0.0002 vs 0.0229 at `block` L=2; 0.0025 vs 0.0529 at `brick` L=2). **The spectral
profile is a property of the ANSATZ AND THE POINT, essentially not of the cost**, and it does
not explain the difference in QNG's usefulness.

**What the table does establish**, and it is the operative number for QNG: the
**amplification factor** `lam_max / lam_eff` — the ratio between the stiffest direction and
the share-weighted eigenvalue the natural gradient actually divides by — is only **1.4 to
3.1** in every ill-conditioned cell, on every cost. **So even where the metric's condition
number is 2,000, the rescaling QNG applies along the directions the gradient actually
occupies is a factor of two or three.** Whatever QNG buys here, it is not the rescue of a
`10^3` ill-conditioning; the ill-conditioning lives in directions the gradient does not use.

**Two of my own predictions are refuted and both are kept in place**: "QNG will be neutral at
best" (wrong on the KL task at equal iterations, see A3) and "the cost's spectral profile
explains where it helps" (wrong, this table).

---

## Applicability of the coordinator's Phase 0 corrections to this workstream

Checked mechanically, not asserted:

* **`s15/distml.py`'s `LogPTable` grid bug.** `grep -rn "LogPTable|ShiftedLogP|SumLogP|distml"
  s15/qgeom_*.py` returns **nothing**. No QGEOM module imports it; no number here is affected.
* **The `amber_kind == 0 AND amber_idx != snap_index` rule.**
  `grep -rn "amber|AMBER" s15/qgeom_*.py` returns **nothing**. This workstream uses only
  `Enum.rmsd` (post-hoc), `Enum.legacy`, `Enum.leg["steric"]`, `Enum.prior`, the `hamil`
  prior/distogram tabulation and the retrieval torsion prior. **No AMBER quantity enters any
  QGEOM measurement**, so neither the tail rule nor the transposed 21.7%/40% share touches
  anything reported here.

---

## A3 / A4. DOES QNG HELP? — YES ON A DISTRIBUTION-FITTING COST, NO ON THE ACTUAL VQE

`s15/qgeom_qng.py`, analysed by `s15/qgeom_qngreport.py`. Six optimisers on identical
circuits, objectives, initialisations and seeds; **every arm gets its own (lr, lam) sweep on
a fixed tuning seed and is then evaluated on fresh seeds**, so the comparison is
best-tuned against best-tuned and the reported number is not the best of a sweep over the
evaluation seeds. The baselines (`sgd`, `adam`) keep the widest grids throughout — a negative
result about QNG must never rest on an under-tuned baseline.

`F = 4g` (section V1) means **there is no second metric to try**: `qng` and a classical
natural gradient are the same preconditioner here up to the learning rate. `ng_shots`
estimates the same object from 512 finite shots, which is the hardware-available version.

### A3. On a WELL-POSED distribution-fitting cost, QNG helps exactly where the metric is ill-conditioned — DEMONSTRATED

Target: `exp(-6 * uniformise(CA-RMSD))` over a 12-qubit sub-register — deliberately
multimodal. Minimise `KL(target || p_theta)`, 150 iterations, 3 seeds, equal ITERATIONS.
Lower is better.

| pattern | L | cond | sgd | adam | **qng** | qng_diag | qng_adam | ng_shots | winner |
|---|---|---|---|---|---|---|---|---|---|
| ring | 1 | **1.00** | 0.6378 | 0.6488 | 0.6378 | 0.6378 | 0.6488 | 0.6378 | tie |
| all_to_all | 2 | 2.97 | **0.7976** | 0.8146 | 0.8042 | 5.4590 | 0.9961 | 0.8082 | sgd |
| ring | 2 | 9.69 | 0.7916 | 1.1905 | 0.6614 | 0.6876 | 1.7159 | **0.6484** | ng_shots |
| ring | 3 | 9.73 | 1.9790 | 2.3319 | **1.0955** | 2.5285 | 2.3701 | 2.5559 | qng |
| block | 2 | 29.81 | 0.6100 | 0.6092 | **0.6092** | 0.6094 | 0.6092 | 0.6092 | qng |
| chain | 2 | 162.30 | 0.5647 | 0.9177 | **0.5415** | 0.7798 | 0.7552 | 0.6302 | qng |
| brick | 2 | 207.91 | 1.0406 | 1.0826 | **0.7948** | 0.8803 | 1.6348 | 0.8225 | qng |

**A clean threshold.** Where `cond > 5` (5 cells) **QNG beats plain gradient descent 5 out of
5**, by a mean of **0.257 in KL** and by as much as 0.883 (`ring` L=3: 1.096 vs 1.979). Where
`cond < 5` (2 cells) it wins **0 of 2** — and at `cond = 1.00`, where `g = I/4` exactly, it is
bit-for-bit the same optimiser as SGD, as the depth-1 theorem requires. **Sprint 14's
depth-1 refutation is confirmed at depth 1 and does not extend: the metric earns its keep at
depth >= 2 on this cost.**

**This refutes my own pre-registered prediction** (recorded in A2 before A3 was read): "QNG
will be neutral at best and will need heavy Tikhonov regularisation to avoid being worse."
It is wrong on this cost, and I record it as refuted rather than deleting it.

### A4. On the REAL VQE cost, QNG helps nowhere — DEMONSTRATED

Native-free Legacy sub-objective on a 12-qubit sub-register of three enumerated targets,
certified global optimum known by enumeration, alpha in {1.0, 0.25}, 4 conditioning rungs,
2 seeds, 90 iterations, equal ITERATIONS. **24 cells.** Objective gap to the certified
optimum; lower is better; negative paired difference = the first arm is better.

| conditioning band | n | sgd | adam | **qng** | qng_diag | qng_adam | ng_shots |
|---|---|---|---|---|---|---|---|
| `cond ~ 1` (g = I/4) | 6 | 0.19904 | **0.05162** | 0.22273 | 0.22273 | 0.07639 | 0.22914 |
| mild (2-100) | 1 | 0.00324 | **0.00025** | 0.00239 | 0.00321 | 0.00532 | 0.02614 |
| **ill-conditioned (>100)** | **17** | 0.14843 | **0.08346** | 0.15003 | 0.16325 | 0.09308 | 0.16556 |

Paired, in the ill-conditioned band (17 cells, condition numbers 688 to 2,677):

| comparison | mean diff | CI95 | W/L | mean/sd | verdict |
|---|---|---|---|---|---|
| **qng - sgd** | **+0.00160** | [-0.01001, +0.01353] | **7/10** | +0.06 | **NULL** |
| **qng - adam** | **+0.06658** | [+0.03924, +0.09639] | **0/17** | +1.11 | **SIGNIFICANT (worse)** |
| qng_adam - adam | +0.00962 | [-0.01067, +0.03343] | 6/11 | +0.20 | NULL |
| qng_diag - sgd | +0.01482 | [+0.00874, +0.02052] | 2/15 | +1.17 | SIGNIFICANT (worse) |
| ng_shots - sgd | +0.01713 | [+0.00604, +0.02891] | 7/10 | +0.67 | SIGNIFICANT (worse) |
| **adam - sgd** | **-0.06498** | [-0.09297, -0.03882] | **16/1** | -1.12 | **SIGNIFICANT (better)** |

Null-calibrated concentration on `qng - adam`: the top 10% of cells carry 0.165 of the
absolute effect against a uniform expectation of 0.059 — **DIFFUSE**, so the effect is not
carried by a handful of cells.

**A4a. The answer to the brief's question.** *On the cost this project actually optimises,
QNG does not help anywhere on the conditioning range.* It is **statistically
indistinguishable from plain gradient descent** even where the Fubini-Study condition number
is 688-2,677 (W/L 7/10, CI spanning zero), and it **loses to plain Adam in 17 of 17 cells**.
The hardware-available sampled version (`ng_shots`, 512 shots) is significantly worse than
SGD. Preconditioning Adam's step with the metric (`qng_adam`) adds nothing (NULL).

**A4b. Why, with the measurement that supports it.** A5 shows the amplification factor
`lam_max / lam_eff` — the rescaling the natural gradient actually applies *along the
directions the gradient occupies* — is only **1.4 to 3.1** even where the condition number is
2,000, because the gradient sits on the stiff directions (A2). **The ill-conditioning is real
and lives where the gradient is not.** A per-parameter adaptive method (Adam) captures more of
the available rescaling, for free, than the exact metric does.

**A4c. Why the two costs differ.** The KL cost demands the state match a target on all `2^n`
amplitudes, so it constrains many independent directions and the metric is the right currency
for trading them off; and its `t/p` weighting gives it a dynamic range that destabilises
Adam's second-moment estimate (Adam is the WORST arm in three A3 cells). The VQE cost is a
single linear functional of `p`: one descent direction suffices, and Adam's diagonal scaling
is enough. **HYPOTHESIS** — stated as such; the direct test (QNG's advantage as a function of
the number of modes the target constrains) was not run.

### A4d. The structural axis: the optimiser is structurally irrelevant, and all six beat the certified optimum — DEMONSTRATED

Mean mode-RMSD over the same 24 cells:

| arm | sgd | adam | qng | qng_diag | qng_adam | ng_shots | **certified optimum** | random draw |
|---|---|---|---|---|---|---|---|---|
| mode RMSD (A) | 4.349 | 4.390 | 4.332 | 4.398 | 4.372 | 4.342 | **4.684** | 4.381 |

**Two readings.**

1. **Six optimisers spanning a threefold range in objective gap land within 0.07 A of each
   other.** The optimiser choice does not move the structure. Whatever the differences in
   A4's objective column mean, they do not reach the answer.
2. **The certified global optimum of this native-free Legacy sub-objective is 4.684 A,
   0.303 A WORSE than a uniform random draw at 4.381 A — and every arm beats it** by
   0.29-0.35 A, purely by failing to reach it. This is Sprint 13/14's central negative
   (Legacy's certified argmin is worse than random: +0.139 A on the full n=9 register)
   reproduced on a different instrument, a different register width and with six optimisers
   instead of one. *The budget trap is not about budget: at fixed budget, the arms that
   optimise less well return the better structures.*
   Per-arm `rho(objective gap, mode RMSD)` runs -0.32 to +0.13 across the six arms — weak and
   mixed, so the *correlation* is not a reliable effect and is reported here only to say so;
   **the reliable statement is the level comparison against the certified optimum.**

### A3b / A4b. AT EQUAL HARDWARE COST — and the conclusion CHANGES on one of the two costs

The brief requires both budget conventions and requires any conflict to be reported.
Convention B charges a parameter-shift gradient `2P` circuit evaluations, the full quantum
geometric tensor a further `P(P+1)/2`, and the diagonal a further `P`. At `P = 24` that is
**348 evaluations per QNG step against 48 for a plain gradient step**, so QNG gets 7.25x
fewer iterations for the same hardware; at `P = 36` it is 10x.

**A4b, the real VQE, equal cost (8,000 circuit evaluations), 1CS9 across the ladder:**

| cell | P | iters (sgd / adam / qng / qng_diag) | sgd | adam | **qng** | qng_diag |
|---|---|---|---|---|---|---|
| a=1.0 ring L1 | 12 | 333 / 333 / 78 / 222 | 0.10381 | **0.00464** | 0.25599 | 0.22804 |
| a=1.0 ring L2 | 24 | 166 / 166 / 22 / 111 | 0.03236 | **0.00926** | 0.21307 | 0.04750 |
| a=1.0 chain L2 | 24 | 166 / 166 / 22 / 111 | 0.04309 | **0.00681** | 0.16516 | 0.06417 |
| a=1.0 chain L3 | 36 | 111 / 111 / 10 / 74 | 0.01268 | **0.00977** | 0.26178 | 0.01512 |
| a=0.25 ring L1 | 12 | 333 / 333 / 78 / 222 | 0.31709 | **0.06356** | 0.37276 | 0.36767 |
| a=0.25 ring L2 | 24 | 166 / 166 / 22 / 111 | 0.35898 | **0.15342** | 0.40209 | 0.36519 |
| a=0.25 chain L2 | 24 | 166 / 166 / 22 / 111 | 0.28784 | **0.16519** | 0.34341 | 0.33604 |

**No conflict on the VQE cost: QNG loses at equal iterations and loses far worse at equal
cost** — by up to 27x on the objective gap (`chain` L3: 0.262 vs Adam's 0.0098). 11 cells run,
7 shown; every one has Adam first and `qng` last or next-to-last.

**A3b, the KL cost, equal cost (6,000 circuit evaluations) — HERE THE CONCLUSION FLIPS:**

| pattern | L | cond | sgd | adam | **qng** | **qng_diag** | ng_shots | winner | A3 winner (equal iters) |
|---|---|---|---|---|---|---|---|---|---|
| ring | 1 | 1.00 | **0.6378** | 0.6433 | 0.6378 | 0.6378 | 0.6378 | sgd | tie |
| all_to_all | 2 | 2.97 | **0.7976** | 0.8151 | 1.2761 | 1.7504 | 0.8183 | sgd | sgd |
| ring | 2 | 9.69 | 0.7918 | 1.9632 | 1.4694 | **0.6977** | 1.0389 | **qng_diag** | qng (0.661) |
| ring | 3 | 9.73 | 2.6546 | 2.0002 | 2.0863 | 2.2212 | **1.8959** | ng_shots | qng (1.096) |
| block | 2 | 29.81 | 0.6102 | 0.6092 | 0.6092 | 0.6103 | **0.6092** | tie | qng |
| chain | 2 | 162.30 | **0.5655** | 0.9645 | 0.5953 | 0.7966 | 1.3694 | **sgd** | qng (0.542) |
| brick | 2 | 207.91 | 1.0407 | 1.0840 | 1.0045 | **0.8805** | 1.2206 | **qng_diag** | qng (0.795) |

**THE CONFLICT, stated plainly.** On the KL cost the full quantum geometric tensor beats
plain gradient descent **5 of 5** in the ill-conditioned cells at equal iterations and **2 of
5, marginally, at equal cost**, with two clean sign flips (`ring` L2: 0.661 vs SGD's 0.792
becomes 1.469 vs 0.792; `chain` L2: 0.542 vs 0.565 becomes 0.595 vs 0.566). **The full QGT
wins 0 of 7 cells outright at equal cost.** Reporting only the equal-iteration convention
would have produced a positive QNG headline that the cost accounting removes.

**And the useful residue.** The arm that survives equal-cost accounting is the **diagonal**
preconditioner, which pays only `P` extra evaluations rather than `P(P+1)/2`: `qng_diag` wins
2 of 7 cells outright and is never catastrophic, where the full QGT is. *If any metric-aware
step is worth taking in this setting it is the diagonal one, and the full geometric tensor is
not worth its price on either cost.*

---

## D5 / E1. THE ENTROPY FLOOR — REAL, BUT MY PROPOSED LAW IS REFUTED

`s15/qgeom_ens.py`. Sprint 14 measured that the final distribution BROADENS as alpha falls
and attributed it to the gradient being carried by `alpha * shots` samples, i.e. to shot
noise. **With an exact statevector gradient there is no shot noise at all and the broadening
is still there**, so the recorded mechanism cannot be the whole story. I pre-registered the
alternative — Sprint 14's own iff condition `dCVaR/dp == 0 <=> p(argmin E) >= alpha` read as a
stopping rule, predicting `p(argmin) -> alpha from above`.

9 targets x 2 seeds, 400 iterations, exact gradient, `hamil` w=0.25 on a 12-qubit sub-register:

| alpha | n | **p(argmin E)** | **p/alpha** | max p | final entropy (bits) | ESS fraction | E_p[RMSD] |
|---|---|---|---|---|---|---|---|
| 1.000 | 18 | **0.00000** | — | 0.7184 | 0.698 | 0.00045 | **2.802** |
| 0.500 | 18 | 0.14261 | 0.29 | 0.5333 | 1.942 | 0.00068 | 3.032 |
| 0.250 | 18 | 0.14223 | 0.57 | 0.3007 | 3.849 | 0.00165 | 3.322 |
| 0.100 | 18 | 0.15118 | 1.51 | 0.2012 | 4.904 | 0.00315 | 3.501 |
| 0.050 | 18 | 0.15812 | 3.16 | 0.1774 | 5.408 | 0.00426 | 3.471 |
| 0.025 | 18 | 0.11177 | 4.47 | 0.1551 | 5.727 | 0.00590 | 3.526 |
| 0.010 | 18 | 0.08053 | 8.05 | 0.1134 | 6.278 | 0.00896 | 3.645 |
| 0.005 | 18 | 0.05370 | 10.74 | 0.1062 | 6.470 | 0.00950 | 3.638 |

**E1a. My prediction is REFUTED.** `p(argmin)/alpha` runs from 0.29 to **10.74** and is never
near 1 except by accident at alpha=0.1. `p(argmin)` is essentially **flat at 0.05-0.16 across
two orders of magnitude in alpha** — it does not track alpha at all. Recorded as my fourth
refuted own-hypothesis.

**E1b. What survives, and it is still worth stating.** The stopping *condition* `p(x*) >=
alpha` is satisfied in **every** row with alpha <= 0.1 and in **no** row with alpha >= 0.25.
So for alpha <= 0.1 the optimiser really has reached the exact global minimum of `CVaR_alpha`
and its gradient is exactly zero — but it **overshoots** the condition by 1.5x to 10.7x rather
than stopping at it, because Adam's momentum carries it past the boundary. **`alpha` is a
monotone diversity dial (entropy 0.70 -> 6.47 bits) whose stopping condition is exact and
whose stopping POINT is not.**

**E1c. A separate result hiding in the first row.** At `alpha = 1` the expectation-value VQE
converges to a point mass (`max p = 0.718`, entropy 0.698 bits) on a state that carries
**`p(argmin E) = 0.00000`** — *it is not the objective's global minimum, in 18 of 18 runs.*
The variational landscape traps the optimiser away from the certified optimum, which is the
curvature-side statement of Sprint 14's "VQE reaches the certified optimum in 0% of cells" and
of A4's finding that no arm closes the objective gap.

**E1d. Diversity costs expected quality, monotonically.** `E_p[RMSD]` rises from 2.802 at
alpha=1 to 3.645 at alpha=0.005. So on the *expected* structural axis, concentration is
strictly good — the exact opposite of what the coordinate-average readout says in E2. **The
two readouts disagree over the whole alpha range, and which one is right is a property of the
downstream operator, not of CVaR.**

---

## E2. THE ONE PLACE A VQE WINS: AN ENSEMBLE CONSUMED WITHOUT A RANKER — DEMONSTRATED, hedged

This is the last steelman, and it is the project's actual situation: *nothing ranks within a
pool*. So the question is not "does the VQE find a good structure you can select" but **"is
the distribution it produces better than the one it started from, when consumed with no
objective-based selection anywhere?"** 9 targets x 4 alphas x 3 seeds, 2,048 draws, matched
budget, control = best-of-N from the arm's own untrained initial distribution.

| readout | alpha | VQE | control | diff | CI95 | W/L | verdict |
|---|---|---|---|---|---|---|---|
| **rand-5 coordinate average** | 1.00 | 2.722 | 3.289 | **-0.567** | [-1.080, -0.068] | 18/9 | **SIG** |
| rand-5 coordinate average | **0.05** | 2.930 | 3.289 | **-0.359** | [-0.708, -0.033] | 16/11 | **SIG** |
| rand-20 coordinate average | 1.00 | 2.693 | 3.072 | -0.380 | [-0.798, +0.018] | 17/10 | NULL |
| **rand-75 coordinate average** | 1.00 | 2.673 | 3.110 | **-0.437** | [-0.836, -0.032] | 15/12 | **SIG** |
| whole-set coordinate average | 1.00 | 2.602 | 3.015 | -0.414 | [-0.843, +0.012] | 16/11 | NULL |
| drawn-set mean | 1.00 | 2.863 | 3.710 | **-0.847** | [-1.240, -0.450] | 22/5 | **SIG** |
| **drawn-set best** | 1.00 | 2.459 | 1.413 | **+1.046** | [+0.725, +1.395] | **1/24** | **SIG (worse)** |
| mean pairwise RMSD (diversity) | 1.00 | 1.542 | 3.074 | **-1.532** | [-2.053, -1.024] | 24/3 | SIG (less diverse) |
| mean pairwise RMSD (diversity) | 0.05 | 3.008 | 3.074 | -0.067 | [-0.169, +0.032] | 16/11 | NULL |
| distinct configurations | 0.05 | 315 | 457 | -142 | | 24/3 | SIG |

**E2a. The result.** Consumed with **no ranker**, a CVaR-VQE ensemble beats best-of-N from its
own untrained start by **0.36 to 0.57 A** on a coordinate-average readout, with CIs excluding
zero, at both alpha=1 and alpha=0.05. **This is the only positive VQE result in this project.**

**E2b. Three hedges, all of which must be read with it.**

1. **At alpha=1 it is not an ensemble.** The final distribution has 2.7 distinct members in
   2,048 draws and a mean pairwise RMSD of 1.54 A against the control's 3.07. Its "5-member
   coordinate average" is one structure. The honest statement at alpha=1 is *a single
   converged structure beats the coordinate average of an equal number of unranked random
   draws.* **At alpha=0.05 the hedge does not apply**: the ensemble keeps 315 distinct members
   and a mean pairwise RMSD statistically indistinguishable from the control's (-0.067, NULL),
   and it still wins the rand-5 coordinate average by 0.359 A. *That cell is a genuine
   ensemble result.*
2. **It buys the mean and pays the best, as always.** `set_best` is **1.046 A worse** at
   alpha=1 with W/L 1/24. The gain is collectible only because no ranker exists to find the
   control's better member; the moment a ranker exists the sign reverses. This is the same
   trade as D2a, priced on the same targets.
3. **The null-calibrated concentration check PASSES but with LOW POWER.** `s15/qgeom_nullconc.py`
   simulates a uniform-effect-plus-observed-noise null: the observed top-10% share of 0.224
   (rand-5, alpha=1) against a null median of 0.192, `p = 0.145`; `p_drop = 0.458`. **Verdict
   DIFFUSE (PASS)** — but `mean/sd = -0.40`, below the brief's 0.5 power threshold, so the
   PASS is weak evidence and not strong evidence. **Tier: DEMONSTRATED at n=27 with a
   low-power concentration PASS; it should be replicated before it carries any architecture
   decision.**

---

## CORRECTION TO MY OWN SECTION D2, FROM THE MANDATED NULL CHECK

Applying `qgeom_nullconc` to the D2 table retrospectively **downgrades one of my own
"SIGNIFICANT" rows**:

| D2 claim | mean | bootstrap verdict | **null-calibrated concentration** | corrected status |
|---|---|---|---|---|
| argmin, alpha=0.25 | +0.251 | SIGNIFICANT | share 0.277 vs null 0.191, **p = 0.006** | **CONCENTRATED (FAIL)** — downgrade |
| argmin, alpha=0.05 | +0.068 | NULL | p = 0.004 | CONCENTRATED (FAIL) |
| argmin, alpha=0.01 | +0.048 | NULL | p = 0.009 | CONCENTRATED (FAIL) |
| top-20 coordavg, alpha=0.25 | +0.372 | SIGNIFICANT | p = 0.069, mean/sd +0.81 | **PASS, holds** |
| top-75 coordavg, alpha=0.25 | +0.392 | SIGNIFICANT | p = 0.458, mean/sd +0.69 | **PASS, holds** |
| top-20 coordavg, alpha=0.05 | +0.197 | NULL | p = 0.050 | CONCENTRATED (FAIL) |
| drawn-set mean, alpha=1.0 | -0.847 | SIGNIFICANT | p = 0.477, mean/sd -0.79 | **PASS, holds** |
| drawn-set best, alpha=1.0 | +1.046 | SIGNIFICANT | p = 0.163, mean/sd +1.17 | **PASS, holds** |

**So the corrected D2 verdict is:** the *argmin* comparisons at alpha <= 0.25 are carried by a
few targets and must not be quoted as general effects; the *coordinate-average* comparisons at
alpha=0.25 and the *set mean / set best* trade at alpha=1 survive the null calibration with
adequate power and stand. **The direction of every one is unchanged — CVaR-VQE is worse than
best-of-N on every ranked readout — but the argmin numbers should be quoted as
target-heterogeneous rather than as a clean mean effect.**

---

## E. THE LANDSCAPE — CONVERGENCE IS NOT A BASIN, MEASURED WITH EXACT CURVATURE

`s15/qgeom_land.py`. **Exact Hessians**, not finite differences: for
`C = sum_x E(x) psi(x)^2` with real amplitudes, `d_i d_j psi = (1/4) psi_{i+j}` for `i != j`
and `d_i^2 psi = -(1/4) psi` follow from the same `RY(t+pi)` identity used throughout, so

    H_ij = 2 sum_x E(x) [ (d_i d_j psi)(x) psi(x) + (d_i psi)(x)(d_j psi)(x) ]

is exact in `P(P-1)/2 + P` statevector builds. **Verified against central differences of the
exact adjoint gradient: max relative error 5.22e-11.**

9 enumerated targets x 2 ansatz x 2 seeds = **36 runs**, 12-qubit sub-register, native-free
Legacy objective, 200 Adam iterations, curvature taken at iteration 0, 100 and 200.

| point | \|grad\| | lam_max | lam_min | **# negative dirs** (of 24) | # flat | \|cond\| | eff. rank | **grad share in negative curvature** | cost | **mode RMSD** |
|---|---|---|---|---|---|---|---|---|---|---|
| init | 0.18160 | +0.2170 | -0.1572 | **10.92** | 0.00 | 327 | 15.92 | **0.473** | 0.4198 | 3.784 |
| mid | 0.00277 | +0.4923 | -0.0013 | 3.81 | 0.00 | 1.9e7 | 11.80 | 0.080 | 0.0187 | 3.774 |
| **final** | **0.00058** | +0.5443 | +0.0010 | **2.28** | 1.22 | **9.8e7** | 11.54 | **0.272** | 0.0157 | **3.735** |

**E1. Convergence is real and the structure does not follow it.** The gradient norm falls
**313-fold** and the cost falls **27-fold**, and the mode RMSD moves **3.784 -> 3.735 A, a
total of 0.049 A over the whole optimisation.** The brief's warning is met with a
measurement: this optimiser converges, and converging does not take it anywhere structurally.

**E2. Only a quarter of the runs reach a genuine local minimum.** At the final point the
Hessian is positive semi-definite in **25% of cells**; the other 75% still carry a mean of
**2.28 negative-curvature directions** and **27.2% of the gradient's squared norm lies along
them**. The optimiser stops on a saddle-adjacent plateau because the gradient has become
small, not because the curvature says it has arrived. (At initialisation 10.92 of 24
directions are negative and PSD holds in **0%** of cells, so the starting point is deep in
saddle country.)

**E3. The two ansatz go in OPPOSITE structural directions while both converge.** This is the
most informative split in the table:

| ansatz | metric rank | point | \|grad\| | # neg | # flat | PSD | cost | **mode RMSD** |
|---|---|---|---|---|---|---|---|---|
| **ring** L2 | 24 / 24 | init | 0.13763 | 10.06 | 0.00 | 0% | 0.4381 | 3.566 |
| | | mid | 0.00329 | 1.06 | 0.00 | 22% | 0.0332 | 3.988 |
| | | final | 0.00084 | **0.61** | 0.06 | **50%** | 0.0287 | **3.932** |
| **block** L2 | **18 / 24** | init | 0.22557 | 11.78 | 0.00 | 0% | 0.4015 | 4.002 |
| | | mid | 0.00225 | 6.56 | 0.00 | 0% | 0.0043 | 3.560 |
| | | final | 0.00032 | **3.94** | **2.39** | **0%** | 0.0028 | **3.539** |

* `ring` converges to a clean minimum (PSD in half its runs, 0.61 negative directions) and
  its structure gets **0.37 A WORSE** while its cost falls 15-fold.
* `block` never reaches a minimum in **0 of 18** runs, drives the cost down **143-fold** — the
  best optimisation in the study — and its structure gets **0.46 A BETTER**.
* **The ansatz that optimises best reaches no minimum, and the ansatz that reaches a minimum
  makes the answer worse.** The population mean (0.049 A of movement) is a *mixture of two
  opposite effects* and must never be quoted alone. This is the budget trap seen through
  curvature: reaching the objective's basin is what costs the structure.

**E4. The metric's rank deficiency shows up in the curvature.** `block`'s metric is exactly
rank 18 of 24 (A1b), and it is the only arm whose final Hessian carries flat directions:
**2.39 against `ring`'s 0.06, a 40x difference**, and it is the only arm that never reaches
PSD. The two objects are not the same — the Hessian's flat count (2.39) is smaller than the
metric's exact redundancy (6), so the redundant parameter directions are not all numerically
flat in the Hessian at the tolerance used — but the qualitative link is present and is the
only place in this study where a metric property predicts a landscape property.

**E5. Hessian conditioning and metric conditioning behave differently.** The Hessian's
absolute condition number **explodes** during optimisation, 327 -> 1.9e7 -> 9.8e7, and its
effective rank falls from 15.9 to 11.5 of 24 — the landscape becomes effectively
half-dimensional as the optimiser converges. The **metric's** rank is constant at 24 (`ring`)
and 18 (`block`) throughout. *A well-conditioned metric is no guarantee of a well-conditioned
landscape, and the object QNG preconditions with is not the object that makes the landscape
hard.* This is an independent reason, on top of A2 and A5, why preconditioning by `g` does
not fix the optimisation.

---

## C6-CORRECTED, AND THE FINAL CLOSURE OF H2

**Retraction of my own cell.** `qgeom_grad.conditioning_effect` reported `Var[g_1]` over seeds
for the `uniform` and `COND` arms, but those arms are a *single point with a 0.01 jitter*, so
the across-seed variance is small by construction and is not comparable to the Haar column
(which varies over the whole parameter space). Its numbers (1e-7 against Haar's 3e-4) mean
nothing and are withdrawn. `s15/qgeom_c6.py` reports instead the statistics of the gradient
**at** the point — its norm and its per-component variance — which is what "trainability at
this point" means. 9 targets x 3 seeds, `block` L=2.

| comparison | grad norm | grad var (components) |
|---|---|---|
| COND - uniform | -0.0186 [-0.0416, +0.0027] 14/13 **NULL** | +0.00002 [-0.00029, +0.00030] 12/15 **NULL** |
| COND - haar | -0.0249 [-0.0532, +0.0019] 16/11 **NULL** | -0.00038 [-0.00083, +0.00005] 16/11 **NULL** |
| **COND - scram_state** | -0.0328 [-0.0628, -0.0061] 16/11 *SIG* | -0.00053 [-0.00105, -0.00009] 17/10 *SIG* |

**Trainability is not degraded by conditioning**: the gradient at the target-conditioned point
is the same size as at a Haar-random point and at the uniform point (both NULL). *Informative
state preparation costs geometry (B2) and does not cost gradient.*

**The one surviving candidate for a target-specific effect, and it dies to the mandated
check.** `COND - scram_state` on gradient magnitude is the single statistic in this entire
workstream whose bootstrap CI excluded zero on the target-specificity axis — and its W/L is a
near-even 16/11 and 17/10, which the brief flags as the signature of a concentrated effect.
Run through the null-calibrated concentration test (`s15/qgeom_nullconc.py`, uniform-effect
plus observed-noise simulation):

| statistic | ansatz | mean | mean/sd | top-10% share | null median | p | **verdict** |
|---|---|---|---|---|---|---|---|
| grad norm, COND - scram_state | block | -0.0328 | -0.42 | 0.329 | 0.191 | **0.001** | **CONCENTRATED (FAIL)** |
| grad var, COND - scram_state | block | -0.00053 | -0.39 | 0.374 | 0.192 | **0.000** | **CONCENTRATED (FAIL)** |
| grad var, COND - scram_state | ring | -0.000062 | -0.13 | 0.310 | 0.192 | **0.001** | **CONCENTRATED (FAIL)** |
| grad norm, COND - haar | block | -0.0249 | -0.34 | 0.237 | 0.192 | 0.080 | DIFFUSE (PASS) |
| grad var, COND - uniform | block | +0.00002 | +0.02 | 0.232 | 0.192 | 0.106 | DIFFUSE (PASS) |

**It is carried by two or three of twenty-seven cells on both ansatz.** With that removed,
**H2 has no surviving evidence anywhere in this workstream: 54 paired cells across two
ansatz, seven geometric statistics, three trainability statistics, and every comparison of
the target-conditioned point against an entropy-matched scramble of its own prior is null or
concentration-driven.**

---

## WHAT I REFUTED, INCLUDING MY OWN HYPOTHESES

**My own, kept in place with the evidence:**

1. **"QNG will be neutral at best."** Pre-registered in A2 from the gradient's spectral
   profile. **Wrong on the KL cost**, where QNG beats plain gradient descent 5/5 in the
   ill-conditioned cells at equal iterations (A3). Right on the real VQE cost (A4).
2. **"The KL gradient loads the metric's small eigendirections, which is why QNG helps
   there."** **Refuted by direct measurement** (A5): the KL gradient loads the bottom
   eigen-decile *less* than the expectation gradient does (0.0002 vs 0.0229 on `block` L=2).
   The spectral profile is a property of the ansatz and the point, not of the cost.
3. **"`p(argmin E) -> alpha` at convergence"** — my proposed exact form of the entropy-floor
   law. **Refuted**: the ratio runs 0.29 to 10.74 and `p(argmin)` is flat at 0.05-0.16 across
   two orders of magnitude in alpha (E1a). The stopping *condition* holds for alpha <= 0.1;
   the stopping *point* is set by momentum overshoot, not by alpha.
4. **My own cost-locality measurement was confounded** by an uncontrolled term count at
   `w = 1` and `w = n`; the corrected version (C4) reverses the practical conclusion.
5. **My own C6 cell** measured across-seed variance at a single jittered point and was
   meaningless; withdrawn and replaced.
6. **My own D2 "argmin at alpha=0.25 is significantly worse"** fails the null-calibrated
   concentration check (p = 0.006) and is downgraded to target-heterogeneous.

**Of the record I inherited:**

7. **Sprint 14's metric condition numbers and rank column** are statistics of a seed-averaged
   metric, are not properties of the ansatz, and do not reproduce from the recorded code at
   either width. The underlying qualitative claim survives and is stronger than recorded.
8. **"Cost locality does not explain trainability at 6-18 qubits; what matters is n"** —
   refuted in its operational reading. On a term-count-controlled family, locality is the
   **larger** axis (6.05 log2 across weight at n=14, against 2.8 log2 across the whole
   6-to-18-qubit width range), and the 1-local retrieval prior has **no measurable width
   decay at all** (slope +0.006).
9. **The paper's own thesis (H2)** — that target-specific conditioning reshapes the manifold.
   Refuted on 54 paired cells across two ansatz, with the only apparent exception failing the
   concentration check.

## WHAT REMAINS OPEN

* **Why QNG helps on a distribution-fitting cost and not on a VQE cost.** A4c offers a
  framing (the KL cost constrains many independent directions of the state and has a dynamic
  range that destabilises Adam; the VQE cost is one linear functional of `p`) but it is a
  **HYPOTHESIS**. The direct test — QNG's advantage as a function of how many modes the target
  constrains — was designed and not run.
* **A3's remaining two ladder cells** (`chain` L=3, `brick` L=3, the two most ill-conditioned
  points) at equal iterations, and A4b beyond 11 cells. The direction is established in both;
  the cells are missing for completeness, not for the conclusion.
* **Everything is n=9 targets, k=4, 12-qubit sub-registers** (with A5, A1, A2, C1-C5 running
  to 18 qubits). Whether the conditioning result survives at chain lengths where the distogram
  carries more many-body content is untested.
* **The E2 ensemble-generator positive is low-powered** (mean/sd = -0.40) and should be
  replicated on the full 126-target instrument before it carries an architecture decision.

---

## REPRODUCTION

    python -m s12.instrument          # pinned constants, run at start and at end

    python -m s15.qgeom_metric        # V verification; A1 metric per point; A1b rank law;
                                      #   A2 gradient-in-range theorem + spectral profile
    python -m s15.qgeom_align         # A5 the profile is a property of the ansatz, not the cost
    python -m s15.qgeom_qng           # A3 / A3b / A4 / A4b -- does QNG help, both budget
                                      #   conventions.  RESUMABLE: rerun to continue.
    python -m s15.qgeom_cond          # B1/B2 conditioning vs entropy-matched scramble
    python -m s15.qgeom_condlaw       # B3 the positive form: geometry as a function of H(p)
    python -m s15.qgeom_grad          # C1/C2/C3/C5 scaling in width, depth, alpha, energy
    python -m s15.qgeom_loc           # C4-CORRECTED cost locality, term count controlled
    python -m s15.qgeom_c6            # C6-CORRECTED trainability AT the conditioned point
    python -m s15.qgeom_cvar          # D1 trajectory, D2 control test, D3 regimes, D4 constraint
    python -m s15.qgeom_ens           # E1 entropy-floor law, E2 ensemble generator
    python -m s15.qgeom_land          # Part E exact Hessians at init / mid / final

Every module writes `s15/results/qgeom_<tag>.json` **after every cell** via
`qgeom_lib.ck`, which merges with the file on disk before writing and retries on the
Windows `os.replace` lock, so a partial run is never lost and a concurrent reader cannot
kill it. `qgeom_qng`, `qgeom_cond`, `qgeom_cvar` and `qgeom_ens` skip cells already in the
checkpoint, so an interrupted run resumes by rerunning the same command.

Seeds are explicit everywhere: metric/ansatz points use `numpy.random.default_rng(100*s+7)`,
optimiser evaluation seeds are `s in {0,1,2}` with the tuning seed fixed at `0`, and every
`fit_kl` uses `random_theta(circ, seed)`.

Set `OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=2` (the modules also set it at
import) and run **one at a time** on this box.

### Leakage audit

No benchmark file was read or written; `results/benchmark_manifest.json`,
`peptide_folds.json` and `peptide_clusters.json` were never opened. `dev24` was not run.
Native quantities enter only as: `Enum.rmsd` and `s12.instrument.load_univ(...)["nat_ca"]`
for **post-hoc scoring**, and `vqe_lib.blend_objective` in the D3 regime sweep, which is an
**ORACLE DIAGNOSTIC** objective-quality knob and is labelled as such everywhere it appears.
Every conditioning prior is `s14.retprior.state_prior`, a BLOSUM retrieval output. Nothing
under `core/`, `s5/ s7/ s8/ s9/ s12/ s13/ s14/` or `tests/` was modified; `qansatz.py` and
`core/quantum.py` were imported and not changed.
