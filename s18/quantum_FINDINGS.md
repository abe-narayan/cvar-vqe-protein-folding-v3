# Sprint 18 — QUANTUM / ADVERSARIAL workstream findings

**Agent: QUANTUM/ADVERSARIAL.** Two jobs, the second more important: prepare the genuine
CVaR-VQE test on the new objective, and **actively try to falsify the degree-1 hypothesis**.

Pre-registration: `s18/PREREG_quantum.md`, written before any variational arm ran.
Modules: `s18/q_attack.py`, `q_anova.py`, `q_gauge.py`, `q_cond.py`, `q_report.py`,
`q_ceilattack.py`, `q_ceilattack2.py`, `q_temper.py`. Artefacts:
`s18/results/q_{attack_a1,anova,gauge,qcond,qvqe_*,ceilattack,ceilattack2,temper}.json`
and `s18/results/quantum_report.txt`, which is the regenerated output of every report mode.
**No number in this document is typed by hand.**

Tiering: **EXACT** (a theorem) · **ESTABLISHED** · **SUPPORTED** · **PLAUSIBLE** · **OPEN** ·
**REFUTED** · **ORACLE** (needs the native; never predictive).

---

## 0. HEADLINE — the degree-1 result is an artefact of which bit pattern names which torsion state

> **The lattice encodes k = 4 torsion states of a residue in 2 qubits. Which 2-bit code names
> which state is an arbitrary labelling. The deployed objective, its argmin, every RMSD and the
> residue-additive ANOVA object are all invariant under relabelling. The strict Walsh weight-≤1
> object is not — and the entire 2.411 Å result lives in that non-invariance.**

| | mean certified argmin RMSD (ORACLE), n = 19 | vs full |
|---|---|---|
| full deployed objective | **2.661** | — |
| **W1 — strict Walsh weight-≤1, the codebase's labelling** | **2.411** | −0.249 [−0.667, +0.147] |
| **W1 — under a random relabelling of the same states** | **2.852** | **+0.191** [−0.072, +0.530] |
| **RA — residue-additive ANOVA, the object the brief defines** | **2.657** | **−0.004** [−0.393, +0.317] |
| matched-random control (same construction, information deleted) | 4.002 | |
| zero-information control (space mean) | 4.003 | |
| ORACLE space best | 1.030 | |

### 0a. The gauge orbit, written out — this is the load-bearing evidence

**Why `W1` moves and `RA` does not, exactly.** Write `mᵢ[v] = E_μ[E | residue i = v]` for the
per-residue marginal of the objective under uniform μ, v = 0…3. Because χ_q depends only on
residue i's own state, the weight-1 Walsh coefficient at qubit q of residue i is
`c_q = (1/k) Σ_v mᵢ[v] χ_q(v)` — **a function of `mᵢ` alone**. Hence

* `RA` minimises `mᵢ` directly: `v*ᵢ = argmin_v mᵢ[v]`. **Label-free.**
* `W1` keeps only the two single-qubit coefficients of the same 4-vector and picks the two bits
  **independently** by their signs. It discards `mᵢ`'s intra-residue weight-2 coefficient — the
  one that says whether the 4-vector's minimum sits at a corner two independent bit choices can
  reach. **Label-dependent.**

Relabelling the 4 states permutes which value sits at which corner. Nothing physical moves;
`W1`'s argmin does.

**The sampling procedure.** The gauge group is (S₄)ⁿ — 10¹²–10¹³ elements at n = 9–10, so it is
sampled, but the *per-residue* orbit is enumerated **exhaustively**:

1. For each residue i and each of the **24** permutations σ of {0,1,2,3}, form the relabelled
   4-vector `f[c] = mᵢ[σ⁻¹(c)]`, apply `W1`'s independent-bit rule to it, and map the winning
   code back through σ⁻¹. This gives `chosen[i, g]`, the state `W1` picks for residue i under
   labelling g — a complete 9×24 or 10×24 table, no sampling.
2. Because the labelling of each residue is an independent coordinate of the group and `W1`'s
   choice at residue i depends only on residue i's labelling, drawing `g` uniformly per residue
   and reading `chosen[i, gᵢ]` samples the induced distribution over **configurations** exactly.
   4,000 draws per target, `stable_rng(pdb, "gauge", salt="s18quantum")`.
3. Each drawn configuration index is scored by the ORACLE `rmsd` table. No refit, no search.

**Two verifications before any claim.** The weight-1 coefficients computed from the marginals
agree with those from the full 2ⁿ�q fast Walsh–Hadamard transform to **≤ 1.4e−13** on every target,
so step 1's shortcut is exact, not approximate. And the identity labelling's own row reproduces
Sprint 17's `W1` argmin **target by target** (2.411 Å mean), so the object being permuted is the
object that produced the headline.

**The null.** Over 4,000 draws the statistic `deg1 − full` has null mean **+0.191**, median
+0.196, 95 % band **[−0.074, +0.433]**, and the codebase's own labelling sits at **gauge
percentile 0.00** (one-sided p < 0.0001, 0/4000).

**The tie trap, handled.** Five of nineteen targets have a **degenerate orbit** — every labelling
yields the identical configuration — and a strict `<` percentile reads **0.000** on all five,
which would have manufactured a much stronger result out of nothing. Ties are split (mid-rank), so
a degenerate orbit correctly scores **0.500**. Both columns are persisted
(`pct_of_identity_strict`, `pct_of_identity_midrank`). On the 14 non-degenerate targets the
identity labelling's mid-rank percentile is mean 0.356, median 0.288, below 0.5 on 9/14.

**The one interval in the whole degree-1 story that excludes zero measures the arbitrariness of
the encoding:**

```
W1(codebase labelling) - W1(gauge mean)   mean -0.440   median +0.000
                                          CI95 [-0.846, -0.130]   W/L 9/5   folds 5/5
```

— **ESTABLISHED**, `q_report gauge`. **Falsifier F5 fires.** The lattice-to-continuous mapping
is invalid not because the arithmetic is wrong but because **the object carrying the result has
no continuous-torsion analogue at all**: there is no "which bit is which" in continuous torsion
space, so there is nothing to port.

**Independently confirmed by the MATH workstream**, which reached the same conclusion by a
different implementation (`s18/math_lattice.py` L5, 24 random relabellings per target):
their W1 gauge mean **2.862** / **+0.202** against my 2.852 / +0.191, and their
`residue_add_invariant` is **true on 19/19**. Two implementations, agreement to 0.01 Å.

---

## 0b. SECOND HEADLINE — the shipped confidence weighting is VALIDATED, three independent ways

The lane's other result points the opposite way from everything else in the programme's record,
and it is the rarer kind.

> **Deleting the distogram's `1/sd²` confidence weighting costs +0.153 Å [+0.066, +0.247],
> CI excluding zero, 50W/76L, all five folds the same sign — and the full `sd^-p` curve is
> monotone up to the deployed exponent and flat above it. The shipped choice is at the optimum
> and is not a lever.** — **SUPPORTED**, `q_ceilattack2`, `q_temper`, n = 126, **NATIVE-FREE**.

Three independent confirmations:

1. the point comparison — uniform weights lose by **+0.153 [+0.066, +0.247]** (§5d);
2. the whole curve — monotone from p = 0 to p = 2 with **5/5 fold sign consistency at every rung
   below the deployed value**, turning over inside the noise floor above it (§5e);
3. EXPERIMENT's independently derived **ρ(1/sd², |residual|) = −0.476** — the confidences really
   are rank-ordered in the correct direction, which is *why* removing them costs 0.153 Å.

Against a standing record in which AMBER cannot rank, Legacy cannot rank, selection is closed at
five levels and refinement is harmful, this is **the first measured native-free vindication of a
deployed design choice** in this programme. It is stated as a result in its own right, not as an
appendix to the tempering null, because a curve saying *the shipped value is right* is a stronger
statement than a point test — and because it protects the incumbent from a future arbitrary
retune.

---

## 1. THE THREE OBJECTS CALLED "DEGREE-1" ARE NOT THE SAME OBJECT — and the brief's §4 asserts an equivalence that is false

`s18/BRIEF.md` §4: *"When μ is uniform on the enumerated lattice this **is** the Walsh weight-≤1
projection — that equivalence must be verified numerically, not asserted."*

It was verified, and **it is false as stated**. `q_report anova`:

```
max over 19 targets |RA - (weight<=1 + INTRA-RESIDUE weight-2)|   = 1.07e-12   <- EXACT
max over 19 targets |RA - W1| / sd(objective)                     = 2.238
mean over 19 targets |RA - W1| / sd(objective)                    = 1.435
mean Spearman(RA, W1)                                             = 0.804
```

The residue-additive ANOVA object `E₀ + Σᵢ(E_μ[E|θᵢ] − E₀)` equals the projection onto Walsh
coefficients **supported inside one residue** — weight 0, weight 1 **and intra-residue weight 2**
— not the weight-≤1 projection. With k = 4 states in 2 qubits a per-residue field *needs* its
intra-residue weight-2 term, and Sprint 17 measured **95.7 %** of all weight-2 mass as
intra-residue. This is not a rounding difference; it is most of the second-order mass, and it is
the difference between the two objects' argmins. — **EXACT**.

**The consequence for the headline.** The object the brief tells the sprint to port is `RA`, and

```
RA - full   n = 19   mean -0.004   median +0.000   CI95 [-0.393, +0.317]   W/L 3/5 (11 ties)
                     folds same sign 1/5
```

**−0.004 Å is 5 % of the 126-target instrument's 0.084 Å minimum detectable effect.** The brief's
own §6 says a null below that magnitude is uninformative — this one is twenty times below it. The
entire "degree-1 advantage" lives in `W1`, and §0 shows `W1` is a property of the encoding.

**BRIEF §4's caution was right and its equivalence claim was wrong at the same time.** It
correctly flagged "per-qubit ≠ per-residue" in point 1 and then asserted the equivalence in the
sentence above it. Both objects were built, as it instructed, and building both is what caught it.

---

## 2. THE PHASE-0 TENSION IS RESOLVED, AND IT RESOLVES AGAINST DEGREE-1

Phase 0 recorded a tension: `deg1` has a **better argmin** (2.411 vs 2.661) but is a **worse
global correlate of RMSD** (ρ +0.153 vs +0.264). One of the two had to give.

**The decomposition that settles it.** An objective's argmin RMSD is exactly

```
argmin_RMSD  =  [ mean RMSD of the objective's own top 1% ]     <- BAND QUALITY (alignment)
              + [ argmin_RMSD - that band mean ]                <- WITHIN-BAND DRAW (min-of-N)
```

`q_attack a1`, n = 19, fold-clustered paired bootstrap, identity residual 0.00e+00:

| term | mean | median | CI95 | W/L | folds |
|---|---|---|---|---|---|
| TOTAL `deg1 − full` | −0.249 | +0.000 | [−0.667, +0.147] | 7/5 (7 ties) | 4/5 |
| **(a) BAND QUALITY** | **+0.108** | +0.072 | **[+0.011, +0.250]** | **4/15** | **5/5** |
| (b) WITHIN-BAND DRAW | −0.357 | −0.198 | [−0.730, +0.009] | 15/4 | 4/5 |

**143 % of the total mean is carried by the within-band draw, and −43 % by band quality.** The
only component whose interval excludes zero is band quality, **and it points against degree-1**:
the truncation's own top 1 % is 0.108 Å *worse* on average, on 15 of 19 targets, with all five
folds agreeing. That is the same fact the global ρ reports, measured where it matters.

So "better argmin, worse correlation" is **not** a paradox and neither number is an artefact.
Under any monotone-copula model the expected minimum over an objective's own band decreases in
that objective's correlation with RMSD, so the model predicts degree-1's argmin should be
*worse*. It is better because it took a luckier draw from a worse band:

```
deg1 argmin - its own band mean   -0.972  [-1.573, -0.278]
full argmin - its own band mean   -0.615  [-1.229, +0.130]
```

**The min-of-N null, with its band named.** The matched-random control is a uniform draw from the
same top-1 % band the argmin comes from (expectation = the band mean): **deg1 3.383 Å, full
3.276 Å.** The zero-information control is the space mean, 4.003 Å. Both objectives beat both
controls; degree-1 beats neither of them by more than the full objective does.

**And the per-target coupling is the wrong sign for an alignment story:** Spearman(per-target ρ
advantage of deg1, per-target argmin advantage of deg1) = **−0.345 (p = 0.147, n = 19)**. If the
argmin gain came from better ordering this would be strongly negative *and* significant; it is
neither.

### The tie structure, and the three carrying targets

**All 19 targets have `argmin_ties = 1` on the degree-1 objective** — the seven ties in Phase 0's
W/L are *deg1 and full landing on the identical configuration*, not a tied argmin inside either.
So the ledger's tie-breaking trap (`np.argmin` reading the oracle sort order) **did not fire
here**, and I record that as a check that passed rather than a defect found. The seven agreeing
targets are the ones with the highest weight-1 variance fraction (2MJQ 0.784, 2MD2 0.727,
1TOR 0.740) — where the objective is most nearly its own separable part.

The three carriers are **7VI4 (−2.539), 1CS9 (−1.601), 7T3H (−1.365)**, and their common
property is the one that damages the hypothesis: in each the difference is **almost entirely the
within-band draw term** (−2.594, −1.400, −1.330) with band quality contributing +0.055, −0.201,
−0.035. On 7VI4 the degree-1 argmin is the **global best structure of a 1,048,576-configuration
register**, from a band whose mean is *worse* than the full objective's — and it sits at gauge
percentile 0.167, i.e. it is one labelling's luck. On 1CS9 the deg1 argmin sits at the **31st
percentile** of RMSD: it does not find a good structure, the full objective's argmin is simply
bad. Removing the three carriers leaves +0.048 [−0.315, +0.342] against a **uniform-effect null**
of −0.021 [−0.411, +0.366] — at the null's 63.6th percentile, so the drop-3 statistic is
**uninformative on its own**, exactly as the ledger's median-vs-mean law requires it to be
reported.

**Verdict: F4 fires — the effect exists only on the 19-target exhaustive instrument, and on that
instrument it is a min-of-N draw from a band the objective orders slightly worse than the full
objective does.** Combined with §0 and §1: F5 fires too.

---

## 3. THE QUANTUM CONDITIONAL — pre-registered, measured, and it cannot fire on this objective

`s18/PREREG_quantum.md` §2, fixed before any arm ran:

> If the new objective contains useful higher-order correlations, then **(a)** greedy 1-opt
> should cease to certify its optimum at tiny budget, and **(b)** the CNOT-free product ansatz
> should become measurably worse than the entangled VQE at matched budget, with a target-level
> interval excluding zero.

### 3a. Separability — the coordinator's requested measurement, and it is decisive on its own

`q_report anova`, Walsh variance fractions on the full enumerated registers, 19 targets, of the
**rank-uniformised** objective (the ledger's mandatory monotone conditioning):

| object | weight-1 | inside one residue | **INTER-RESIDUE** | mean Pauli weight |
|---|---|---|---|---|
| **full (deployed)** | 0.6133 | 0.9132 | **0.0868** | 1.538 |
| Walsh weight-≤2 | 0.6588 | 0.9823 | 0.0177 | 1.341 |
| **RA (residue-additive)** | 0.6706 | **1.0000** | **0.0000** | 1.329 |
| **W1 (strict weight-≤1)** | **1.0000** | 1.0000 | **0.0000** | 1.000 |

**The truncation raises the weight-1 variance fraction from 0.613 to 0.671 (RA) or 1.000 (W1)
and drives the inter-residue coupling variance to EXACTLY ZERO.** — **EXACT**.

> **The conditional cannot fire on this objective. A degree-1 objective has no correlations for
> an entangling ansatz to represent, because it has no correlations at all. The truncation makes
> the quantum case strictly worse than the full objective already was, and that is a result, not
> an omission.**

Both truncations' global optima are closed forms needing **no search**: `RA`'s is a per-residue
argmin at n·k ≈ 38 table reads (verified against the enumerated argmin on **19/19** targets);
`W1`'s is a per-qubit argmin at 2·n_qubits ≈ 38 reads. — **Q3 fires.**

### 3b. Clause (a) — greedy certifies the new objective FASTER, not slower

`q_cond budget`, 19 targets × 4 seeds, % of cells reaching the **certified** global optimum:

| budget (obj. evals) | **RA** | **W1** | full |
|---|---|---|---|
| 20 | 6.6 % | 0.0 % | 0.0 % |
| **36** (one coordinate pass = n(k−1)+1 = 30) | **100.0 %** | **100.0 %** | 63.2 % |
| 64 | 100.0 % | 100.0 % | 76.3 % |
| 1,024 | 100.0 % | 100.0 % | 100.0 % |
| 8,192 *(the CVaR-VQE budget)* | 100.0 % | 100.0 % | 100.0 % |

A cold Metropolis chain shows the same ordering far more sharply — at 8,192 evaluations it
certifies **100 %** on RA and W1 and only **47.4 %** on the full objective. The zero-information
control (uniform random draw) certifies **0–1.3 %** everywhere.

Paired, target as the unit:

```
greedy cert @   36 evals   RA - full   mean +0.368  med +0.500  CI[+0.197,+0.529]  W/L 0/11  folds 5/5
greedy cert @   64 evals   RA - full   mean +0.237  med +0.000  CI[+0.083,+0.412]  W/L 0/7   folds 5/5
greedy cert @  128 evals   RA - full   mean +0.145  med +0.000  CI[+0.026,+0.292]  W/L 0/4   folds 4/5
```
(identical for W1 at 36/64/128). **Positive = the new objective is EASIER.** — **Q1 fires**,
**ESTABLISHED**.

### 3c. Clause (b) — the entanglement control, run anyway, and null

`q_cond vqe`, 19 targets, seed 0. `mps2fn` is the **identical circuit with the CNOTs removed** —
a product Bernoulli model trained by the identical CVaR score-function estimator, the identical
Adam optimiser, the identical 8,192-evaluation budget, shots (512), learning rate (0.15),
`const` baseline and **the same seed**. It is a classical algorithm.

`mps2f − mps2fn` on the objective `RA`, paired over targets (positive = the entangled circuit is
worse on M and on the readout; D is a *coordinate*, not a quality):

| α | selM75 | selreadout75 | selset_best75 |
|---|---|---|---|
| 0.05 | −0.133 [−0.340, +0.006] | −0.102 [−0.304, +0.039] | +0.056 [−0.094, +0.238] |
| 0.25 | −0.081 [−0.236, +0.092] | −0.043 [−0.271, +0.177] | −0.048 [−0.286, +0.293] |
| 1.00 | −0.017 [−0.195, +0.196] | −0.025 [−0.282, +0.232] | +0.046 [−0.217, +0.413] |

**Every interval on the readout the pipeline consumes spans zero, at every α, on `RA`, on `W1`
and on `full`** (all three tables in `s18/results/quantum_report.txt`). — **Q2 fires**;
Sprint 17's §7 reproduces on the new objectives.

The one systematic difference is on the *uniform-draw* coordinates, and it is the collapse
trajectory Sprint 17's representability theorem predicts, not a quality difference: on `RA` the
entangled circuit is less converged at every α (D75 +0.260 [+0.162, +0.367] at α = 0.05, folds
5/5) and correspondingly higher in member error (M75 +0.310 [+0.045, +0.596]). **Deleting the
CNOTs makes the sampler converge faster, not better.**

For context, on `RA` the plain classical **greedy** arm reaches selreadout75 **2.462** against the
best VQE arm's 2.403 and the product ansatz's 2.448, at 8,192 evaluations it does not need — and
its certified optimum is available at **36**.

### 3d. Verdict on the quantum pillar

**Q1, Q2 and Q3 all fire.** By the pre-registered rule the pillar earns no role: it reaches no
Pareto point the classical controls cannot, the entangled ansatz does not beat its CNOT-free
control, and there is no optimisation problem left to attack. **The quantum branch is closed for
this sprint, and the new objective closes it harder than the old one did.** No quantum advantage
is claimed and none was manufactured.

---

## 4. AUDIT OF THE MATH WORKSTREAM — PASS, and the defect is in the brief, not in their code

Independently re-derived (`q_anova.py`) rather than read.

| check | their result | my independent result | verdict |
|---|---|---|---|
| qubit-factorisation ANOVA ≡ Walsh weight-≤1 | `L1` 8.9e−16 | — | **PASS (EXACT)** |
| residue-factorisation ANOVA ≡ residue-block projection | `L2` ≤ 1.3e−15 | ≤ 1.07e−12 | **PASS (EXACT)** |
| the two objects differ | `argmin` 3.655 vs 5.471 on 1CS9 | identical | **PASS** |
| relabelling invariance (their L5, my gauge orbit) | W1 gauge 2.862 / RA invariant 19/19 | 2.852 / invariant | **PASS, agreement 0.01 Å** |
| variance budget / orthogonality of the ANOVA residual | ≤ 5.0e−16 | — | **PASS** |
| off-mesh continuous reconstruction (`control_D`) | `PASS_machine_precision: false` | — | **flagged below** |

**The `math_anova` self-check does not pass at machine precision and should not be read as if it
did.** `control_D` reports `field_max_abs_err` 7.1e−15 (the field itself is exact) but
`E0_rel_err` **0.875** at S = 512 and **0.905** at S = 2048 — the Monte-Carlo constant `E₀` is
not converging, and `offmesh_E_ge2_max_abs` equals `offmesh_E_le1_max_abs_err` exactly, i.e. the
whole off-mesh error is the constant. `E₀` is an additive constant and therefore changes **no
ranking and no argmin**, so this does not invalidate any degree-1 comparison — but it does mean
`E_le1`'s *absolute value* is not trustworthy, and anything that mixes `E_le1` with another term
at a fixed λ **is** affected. That should be stated where the λ ladder is reported.

**The defect the brief carries and MATH does not.** MATH built both objects and measured them
separately, correctly. The false equivalence is in `BRIEF.md` §4's own sentence. **F5's most
likely quiet death was averted by the brief's own instruction to build both objects.**

---

## 5. ATTACKS ON THE COORDINATOR'S `objceil` HEADLINE — n = 126, complete

Attacks 1 and 2 (the separation-stratified shuffle and the weight–residual pairing) were taken
by the coordinator and are not duplicated. Attacks 3, 4 and 5 are here. `q_ceilattack.py`,
`q_ceilattack2.py`, `q_temper.py`. **Every α > 0 arm and every arm reading `d_true` is an ORACLE
DIAGNOSTIC**; the §5c arms read no native distance at all and are labelled NATIVE-FREE.

### 5a. Attack 3 — α = 1's 1.152 Å is a BASIN, not a ceiling. The ceiling is 0.083 Å.

**The parameterisation floor.** Project the **native itself** into the same ideal-geometry torsion
parameterisation the refinement optimises in: mean **0.083 Å**, median 0.060, sd 0.084, max 0.435.
The parameterisation can represent the native essentially exactly, so it is not the constraint on
any target.

**Multi-start at α = 1** (perfect distances), identical machinery and weights, with the objective
printed beside the RMSD because the objective is what decides convergence:

| start | start RMSD | end RMSD | median | objective at end | W/L vs avgproj |
|---|---|---|---|---|---|
| `avgproj` (the coordinator's start) | 3.213 | **1.152** | 0.307 | **3.04** | — |
| `member` (a random pool member) | 3.571 | 1.253 | 0.554 | 3.86 | 55/71 |
| **`helix` (ZERO-INFORMATION control)** | 4.070 | **1.028** | 0.412 | 3.68 | 61/65 |
| `randtors` (matched-random start) | 5.258 | 2.080 | 2.278 | 18.38 | 37/89 |
| **`nativeproj` (ORACLE start)** | 0.083 | **0.082** | 0.033 | **0.30** | **92/33** |

```
nativeproj_ORACLE - avgproj :  -1.070 [-1.321, -0.829]   median -0.094   92W/33L
helix             - avgproj :  -0.124 [-0.366, +0.109]   median +0.000   61W/65L
randtors          - avgproj :  +0.928 [+0.614, +1.241]   median +0.892   37W/89L
```

**The objective column settles it.** From the ORACLE start the fit reaches objective **0.30**;
from the coordinate average it stops at **3.04** — ten times higher, on the same objective with
the same weights. The optimum is not at 1.152 Å; the optimiser stops there. — **SUPPORTED**.

**And the median must be read beside it.** `avgproj`'s median is **0.307** and the paired median
difference is only **−0.094**. On more than half the targets the coordinator's start already
reaches the native; **the 1.152 Å mean is carried by a minority that fall into a distant basin**.
This is the ledger's median-vs-mean warning firing in the coordinator's favour on the substance
and against the figure: *α = 1 reaches the objective's optimum on most targets and a distant basin
on a minority, and the mean reports the minority.*

**Consequences, stated separately because they point in opposite directions.**

* **The "functional form is sound" headline SURVIVES and is understated.** With perfect distances
  this functional's reachable optimum is **0.083 Å**, not 1.152.
* **The α ladder cannot be read as a requirement curve.** Every rung is one refinement from one
  start, so each rung prices *distance quality* and *basin reachability* jointly. Nothing in the
  ladder isolates the former.
* **The coordinate average is NOT a privileged start.** A constant ideal α-helix — zero
  information, starting 0.86 Å *worse* — lands at 1.028, indistinguishable from the average's
  1.152 (−0.124 [−0.366, +0.109], 61W/65L). This is the zero-information control the arm lacked,
  and it bears on every arm in the programme that begins from the coordinate average.

### 5b. Attack 4 — "halve the residual → 2.5 Å" is one point in a 0.55 Å band that straddles 2.5

**First, a correction that dissolves the stated concern and sharpens it.** The α ladder is

```
d_alpha = (1-alpha)*dhat + alpha*d_true  =  d_true + (1-alpha)*(dhat - d_true)
```

— **exactly a uniform rescaling of the residual vector**, which already preserves its direction
and its entire correlation structure. "Shrink the residual preserving correlation structure" is
what the ladder does. The live objection is different: **a real predictor does not improve
uniformly across pairs — it gets better somewhere.**

A family of error models **all at the identical residual RMS** (0.5 × the deployed, the rung
quoted as 2.5 Å), differing only in *where* the improvement lands, n = 126:

| error model | resid RMS | RMSD | median | vs `unif` | W/L |
|---|---|---|---|---|---|
| **`fix_confident`** (fixes the pairs `sd` called reliable) | 1.604 | **2.145** | 1.865 | **−0.308 [−0.462, −0.150]** | 91/35 |
| **`fix_short`** (fixes |i−j| ≤ 5) | 1.604 | **2.159** | 1.969 | **−0.294 [−0.424, −0.167]** | 87/39 |
| `fix_big` (fixes the largest residuals) | 1.604 | 2.324 | 2.157 | −0.129 [−0.267, +0.004] | 72/54 |
| `debias_sep` (removes the per-separation bias) | 1.604 | 2.349 | 2.177 | −0.104 [−0.207, +0.002] | 64/62 |
| **`unif`** (= the coordinator's α = 0.5) | 1.604 | **2.453** | 2.509 | — | — |
| `fix_small` | 1.604 | 2.487 | 2.319 | +0.034 [−0.065, +0.143] | 66/60 |
| **`fix_long`** (fixes |i−j| > 5) | 1.604 | **2.618** | 2.643 | **+0.165 [+0.028, +0.303]** | 48/78 |
| **`fix_unconfident`** | 1.604 | **2.697** | 2.578 | **+0.244 [+0.097, +0.386]** | 35/91 |

**Spread 2.145 → 2.697 Å at identical residual RMS: a range of 0.552 Å straddling the 2.5 Å
target.** Residual RMS does not determine the outcome, so the requirement must be quoted as a
band, not a point. — **SUPPORTED**.

**And the direction is the opposite of the intuitive one, which makes it actionable.** Improving
the pairs the predictor was *already confident about* (−0.308) and the *short-range* pairs
(−0.294) buys more than a uniform improvement; improving the long-range, low-confidence pairs —
where the errors actually are — is significantly *worse* than uniform. A predictor that reduces
its RMS by fixing its worst pairs would land at 2.62–2.70, not 2.45.

### 5c. Attack 5 — `shuf_paired` prices a WEIGHTING change, not an error-assignment one

`s18/objceil.py` builds its strongest control as `A.fit(pw, sd[pi], i, j, phi0, psi0)` — note
`sd[pi]`. The arm permutes the residual **and** passes a permuted weight vector into the fit, so
it is no longer the deployed functional `Σ_p (d_p − d̂_p)²/sd_p²`. `shuffled` (r permuted, weights
kept) reaches 2.609 Å; `shuf_paired` (same permutation, weights permuted too) reaches 2.072 Å.
The two differ **only** in whether `sd` moves, so the whole **0.537 Å** between them prices the
weight permutation, not where the errors land. **Reported to the coordinator before it shipped;
accepted, and the "assignment to pairs" claim withdrawn.**

### 5d. Attack 5, isolated — the distogram's confidences ARE informative, and deleting them costs 0.153 Å

`q_ceilattack2.py`, n = 126, same start and machinery. The first three arms **read no native
distance at all** and are therefore NATIVE-FREE, not ceilings.

| arm | reads native? | RMSD | median | vs a0 = 3.610 | W/L |
|---|---|---|---|---|---|
| `a0_wkeep` — the deployed objective | no | 3.610 | 3.370 | — | — |
| `wperm_only` — real d̂, weights permuted | **no** | 3.701 | 3.547 | +0.091 [−0.017, +0.197] | 49/77 |
| **`wflat` — real d̂, uniform weights** | **no** | **3.763** | 3.654 | **+0.153 [+0.066, +0.247]** | 50/76 |
| `shufr_wkeep` — r permuted, weights kept | ORACLE | 2.635 | 2.570 | −0.975 [−1.170, −0.783] | 102/24 |
| `shufr_wperm` — r and weights permuted together | ORACLE | 2.167 | 2.055 | −1.442 [−1.680, −1.216] | 113/13 |
| `shufr_wflat` — r permuted, uniform weights | ORACLE | 2.440 | 2.301 | −1.169 [−1.396, −0.942] | 106/20 |
| `a1_wkeep` — perfect distances | ORACLE | 1.152 | 0.307 | −2.458 [−2.775, −2.149] | 117/9 |
| `a1_wperm` — perfect distances, permuted weights | ORACLE | 1.162 | 0.388 | −2.448 [−2.743, −2.153] | 117/9 |

> **A NATIVE-FREE POSITIVE, and they are scarce here. Deleting the distogram's confidence
> estimates costs +0.153 Å [+0.066, +0.247], CI excluding zero, 50W/76L. The confidences are
> informative and `1/sd²` is a defensible functional form for consuming them.** — **SUPPORTED**.

In a programme whose standing record is that AMBER cannot rank, Legacy cannot rank, selection is
closed at five levels and refinement is harmful, a measured native-free vindication of a shipped
design choice is worth stating as a result rather than as an appendix.

**Independently corroborated from a different direction.** EXPERIMENT's H8 leverage diagnostic
measures **ρ(1/sd², |residual|) = −0.476** (n = 74, partial, sign held from n ≈ 45): the pairs the
distogram says it is confident about really do carry smaller errors. That is *why* removing the
weights costs 0.153 Å, measured by a different lane on a different statistic.

**An internal consistency check that validates the whole arm family.** `a1_wperm` (perfect
distances, permuted weights) = **1.162** against `a1_wkeep`'s **1.152**. With perfect distances
the weighting *must* be irrelevant, and it is. Any arm family in which it were not would be
broken.

**And the resolution of §5c's confound is the reverse of what I predicted — recorded as my error,
in place.** I told the coordinator the likely reading was "the distogram's confidences are
anti-informative inside the fit". The data refutes that. What the `shuffled` → `shuf_paired` gap
actually prices is **residual–weight alignment**:

```
shufr_wperm - shufr_wkeep  =  -0.467 [-0.598, -0.338]
```

`shufr_wkeep` gives pair p the residual `r[π(p)]` but the weight `sd[p]` — residual and weight
**mismatched**. `shufr_wperm` gives it both `r[π(p)]` and `sd[π(p)]` — they **travel together**.
So the 0.537 Å between the coordinator's two arms means *orphaning a residual from its own weight
is bad*, confirmed on real data by `wperm_only` (+0.091). **The confound identification was
correct and load-bearing; my prediction of its direction was wrong. Those are separate acts and
both are recorded.**

**The coordinator's core headline is untouched, and I reproduced it independently.**
`shufr_wkeep` = **2.635** on my start against their `shuffled`'s 2.609 on theirs, −0.975
[−1.170, −0.783], 102W/24L — a different start, a different implementation, the same conclusion.
*The distogram's errors are worse than random errors of the same magnitude* has now survived
separation-stratified permutation, isotropic replacement, weight permutation, weight flattening
and an independent reimplementation.

### 5e. Attack 5b — the tempering curve is flat where it matters, and `p = 2` is fine

`q_temper.py`, n = 126, **every arm NATIVE-FREE** (real d̂, real `sd`, no native distance).
Weight = `sd^-p`, realised exactly by passing `s = sd^(p/2)` into the same fit. **The whole curve
is reported, not its argmin.**

| arm | weight | RMSD | median | vs deployed p = 2 | W/L | folds |
|---|---|---|---|---|---|---|
| `p0.0` (zero-information control) | uniform | 3.763 | 3.654 | +0.153 [+0.045, +0.272] | 50/76 | 5/5 |
| `p0.5` | sd⁻⁰·⁵ | 3.732 | 3.567 | +0.122 [+0.038, +0.219] | 47/79 | 5/5 |
| `p1.0` | sd⁻¹ | 3.725 | 3.522 | +0.116 [+0.042, +0.202] | 47/79 | 5/5 |
| `p1.5` | sd⁻¹·⁵ | 3.664 | 3.481 | +0.055 [+0.009, +0.115] | 44/82 | 5/5 |
| **`p2.0` (deployed)** | **1/sd²** | **3.610** | 3.370 | — | — | — |
| `p3.0` | sd⁻³ | **3.584** | 3.506 | **−0.026 [−0.108, +0.056]** | 66/60 | 3/5 |
| `p4.0` | sd⁻⁴ | 3.678 | 3.663 | +0.068 [−0.053, +0.209] | 61/65 | 4/5 |
| **`wperm` (MATCHED-RANDOM)** | 1/sd² permuted | 3.752 | 3.606 | +0.142 [−0.005, +0.288] | 49/77 | 4/5 |

**Neither Reading 1 nor Reading 2. The exponent does not matter once it is at or above ≈ 1.5.**
The curve falls monotonically from p = 0 to p = 2, turns over at p = 3 and rises again at p = 4.
The nominal argmin is p = 3 at 3.584 Å, but:

* it is worth **−0.026 Å**, which is **31 % of the instrument's 0.084 Å minimum detectable
  effect** — by the brief's own §6 that is uninformative, not an improvement;
* its interval **[−0.108, +0.056] spans zero**, W/L is a near-even 66/60 and fold sign
  consistency is 3/5 — three of the ledger's concentration warnings at once;
* and it is an argmin selected on the 126-target **tuning** instrument. It is a **hyperparameter
  chosen on the tuning set, not a validated setting**, and the sealed benchmark stays sealed.

**The honest reading, and it is a perfectly good result: the deployed exponent is fine, and the
exponent is not a lever.** That protects the incumbent from a future arbitrary retune. Given
EXPERIMENT's ρ(1/sd², |residual|) = −0.476 — informative but far from perfectly rank-ordered —
this says the ordering is what carries the signal and the sharpness of its consumption does not.
— **SUPPORTED** (a null above the noise floor at the ends, below it in the middle).

**And the matched-random control gives the sharpest form of §5d.** Permuting the `1/sd²` weights
across pairs — identical multiset, identical entropy, only the assignment destroyed — reaches
**3.752**, and uniform weighting reaches **3.763**:

```
uniform (p = 0) - permuted 1/sd^2  =  +0.011 [-0.120, +0.134]
```

**Scrambling the weights is indistinguishable from deleting them.** So the distogram's confidence
information lives *entirely* in **which pair gets which weight**, not in the spread of the weight
values — and that assignment is worth +0.153 Å. — **SUPPORTED**.

**A LOGIC SLIP IN THIS MODULE'S OWN PRE-REGISTERED READING, caught by the coordinator and
preserved in place.** `s18/q_temper.py`'s pre-registration said *"if uniform merely MATCHES the
permuted weights, `sd` carries no usable weighting information either way"*. **That inference is
invalid.** Uniform and permuted are **two ways of destroying the same information**, so their
agreeing with each other is a **consistency check on the null**, not evidence of absence — and
both are significantly *worse* than the deployed weighting (+0.153 and +0.142). Read as written,
the sentence would have inverted this lane's best positive result. The module now prints the
correction beside the table rather than the original sentence.

---

## 6. WHAT I REFUTED, INCLUDING MY OWN EXPECTATIONS

1. **REFUTED — that "degree-1" names one object.** Three objects were being called that; the two
   that matter differ by 1.4 objective sd and have different argmins.
2. **REFUTED — the degree-1 hypothesis, at its root.** The object the brief defines is worth
   −0.004 Å, twenty times below the instrument's minimum detectable effect. **F2 and F4 fire.**
3. **REFUTED — that the 2.411 Å number is a property of the objective.** It is a property of
   which 2-bit code names which torsion state, at gauge percentile 0.00 with p < 0.0001, and the
   encoding-luck term is the only interval in the story that excludes zero. **F5 fires.**
4. **RESOLVED, against the hypothesis — Phase 0's "better argmin, worse ρ" tension.** Both
   numbers are real; the argmin gain is 143 % within-band draw and the only significant component
   (band quality, +0.108 [+0.011, +0.250], folds 5/5) points the other way.
5. **CHECKED AND CLEAN — the tie-breaking trap.** All 19 `argmin_ties = 1`; Phase 0's seven
   "ties" are the two objectives agreeing, not a tied argmin. Recorded as a check that passed.
6. **REFUTED — my own framing that "the new objective might have higher-order structure".** It
   has exactly zero inter-residue Walsh variance, by construction, on both truncations.
7. **REFUTED — the pre-registered quantum conditional, on both clauses.** Greedy certifies the
   new objective at **36** evaluations in 100 % of cells (against 63.2 % on the full objective),
   and the CNOT-free control is null at every α on all three objectives.
8. **NOT CLAIMED.** No barren-plateau claim, no scaling claim, no quantum win of any kind. The
   VQE arms were run because a forced answer that is never measured is an assertion — and the
   **pre-registration was written, and records the separability measurement as a prior, before
   any variational arm ran.** That order is deliberate: it is the difference between a null and a
   null someone can accuse us of arranging.
9. **REFUTED — a defect in the coordinator's `shuf_paired` control, before it shipped.** The arm
   passed a permuted weight vector into the fit, so it was not the deployed functional; the
   0.537 Å attributed to "where the errors land" prices residual–weight alignment instead.
   Reported, accepted, the claim withdrawn and a downstream lane redirected.
10. **REFUTED — MY OWN prediction of how that confound would resolve, recorded in place.** I said
   the likely reading was "the confidences are anti-informative inside the fit". It is the
   opposite: deleting them costs **+0.153 [+0.066, +0.247]**. Identifying a confound and
   predicting its direction are separate acts; the first was right and load-bearing, the second
   was wrong, and both are on the record.
11. **REFUTED — α = 1's 1.152 Å as the functional form's ceiling.** The ceiling is the
   parameterisation floor at **0.083 Å**; 1.152 is a basin, and the objective column (0.30 from
   the ORACLE start against 3.04 from the coordinate average) is what shows it. The
   coordinator's conclusion survives and is understated; the figure does not.
12. **REFUTED — the α ladder as a requirement curve, and "halve the residual → 2.5 Å" as a point
   figure.** Each rung prices distance quality and basin reachability jointly; and at identical
   residual RMS the outcome spans **2.145–2.697 Å** depending on which pairs improve.
13. **REFUTED — that the coordinate average is a privileged refinement start.** A constant ideal
   α-helix, starting 0.86 Å worse, ends indistinguishable (−0.124 [−0.366, +0.109], 61W/65L).
14. **NOT A LEVER — the weight exponent.** The `sd^-p` curve's nominal argmin (p = 3) is worth
   −0.026 Å with an interval spanning zero, below the instrument's 0.084 Å MDE, chosen on the
   tuning instrument. The deployed p = 2 is fine and should not be retuned.
15. **A LOGIC SLIP IN MY OWN PRE-REGISTERED READING, caught by the coordinator, preserved.**
   `q_temper`'s pre-registration said "if uniform merely matches the permuted weights, `sd`
   carries no usable weighting information either way". Invalid: the two are both ways of
   destroying the same information, so their agreeing is a null consistency check, and what
   decides it is that both are worse than the deployed weighting. As written it would have
   inverted this lane's best positive result.
16. **A PROCESS ERROR OF MY OWN, preserved.** The first `q_temper` launch ran p ∈ {0…2} only; my
   `pkill` did not take, the old process completed, and the sweep past the deployed exponent —
   the one direction in which a monotone curve could still have an optimum — was missing from
   the first report. It was caught by reading the printed table against the module's own `PS`
   constant, and added by `q_temper extend` rather than by rerunning. Had it not been caught, the
   curve would have been reported as monotone with the deployed value at its edge.

---

## 7. HONEST LIMITATIONS

* **n = 19 exhaustively enumerated targets, k = 4, n ≤ 10.** Everything in §§0–3 is on that
  instrument. The gauge argument, however, is **structural** and does not depend on n: there is
  no bit-labelling in continuous torsion space at all, so `W1` has no continuous analogue for any
  instrument size.
* **The 4,000-sample gauge orbit is exact but Monte-Carlo over (S₄)ⁿ** (10¹²–10¹³ elements). Five
  of nineteen targets have a **degenerate orbit** — every labelling gives the same configuration
  — and are scored at mid-rank 0.500, not 0.000; scoring them by strict `<` would have
  manufactured a stronger result and is the tie trap in a new costume. The strict and mid-rank
  columns are both persisted.
* **`q_cond vqe` is one seed (0) per target**, not the 8 seeds of the budget arm. It is a null on
  intervals that are wide enough to have shown a real effect, but it is not a precision
  measurement of a null.
* **The box was contended throughout** (CPU at 97–100 % from three sibling workstreams). Every
  wall-clock number is contended and none is used for a cost comparison; the objective-evaluation
  counts are.
* **`E₀`'s non-convergence in `math_anova`'s self-check** (§4) is inherited by anything that
  consumes `E_le1`'s absolute value. Nothing in this document does.
* **`q_temper`'s p = 3 and p = 4 rungs were added by `extend` after the main sweep**, on the
  identical rows, start, seeds and machinery. Two `extend` processes raced once (a `pkill` that
  did not take); the operation is deterministic and idempotent and the final artefact was
  verified complete with both rungs present on all 126 rows before it was reported.
* **The §5 attacks are on the 126-target TUNING instrument.** Every RMSD there is a tuning-set
  number. Nothing in §5 was validated on the sealed benchmark and nothing in it should be.
* **`q_ceilattack`'s multi-start uses one random start per target per family**, not a multi-start
  distribution; it separates basin-limited from objective-limited, which is what it was built
  for, but it does not price how many basins there are.
* **RA vs the pool-marginal μ.** Everything here uses μ = **uniform on the enumerated lattice**,
  which is what makes the Walsh comparison exact. The brief's alternative μ (the retrieval pool's
  empirical per-residue marginal) is **not** tested here; it is a different, legitimate object and
  the gauge argument does not speak to it. That sensitivity is open and belongs to MATH/EXPERIMENT.

---

## 8. REPRODUCTION

```
python -m s18.q_attack a1            # A1: band quality vs within-band draw, the min-of-N null
python -m s18.q_anova  run           # A2: build W1, RA, RA_walsh, W2 on all 19 registers
python -m s18.q_gauge  run           # A3: the exact relabelling orbit
python -m s18.q_cond   budget        # A4a: greedy/metropolis/random certification vs budget
python -m s18.q_cond   vqe RA        # A4b: entangled vs CNOT-free  (also `W1`, `full`)
python -m s18.q_ceilattack  run      # attacks 3 and 4 on the coordinator's objceil
python -m s18.q_ceilattack2 run      # attack 5: the `shuf_paired` weighting confound
python -m s18.q_temper      run      # attack 5b: the sd^-p tempering curve
python -m s18.q_temper      extend 3.0 4.0   # the rungs past the deployed exponent
python -m s18.q_report {anova,gauge,quantum,vqe}
```

`s18/results/quantum_report.txt` is the regenerated output of every report mode; every table
above comes from it.

### Leakage audit

`Inst.rmsd`, `nat`, `u_truth`, every `*_ORACLE` column, the band-quality and within-band-draw
terms, the gauge orbit's RMSD values and every arm's M/D/readout are **ORACLE** quantities used
for post-hoc scoring only. No objective, truncation, projection, labelling, temperature, α,
budget, threshold or stopping rule reads a native quantity. Every α > 0 arm in §5 is a labelled
ORACLE DIAGNOSTIC by construction. Seeds are `s15/seed.py::stable_rng` with salts `s18quantum`
and `s18qceil`; **`hash()` appears nowhere in `s18/q_*.py`**. **The 60-target sealed benchmark
was not read, not probed and not derived.**
