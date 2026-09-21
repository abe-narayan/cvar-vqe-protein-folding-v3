# SPRINT 31 LEDGER

Entries are `## S31-L<n> -- TITLE (date time, lane)` with the verdict in the heading.
Run `date` in the same command as the append. Corrections are annotated in place with the
original wording left standing.

Ledger numbers are a read-modify-write with no lock and **collided three times in S30**. Take the
next free number and re-check it immediately before writing; if you collide, renumber the
later-stamped entry and annotate in place.

---

## S31-L0 -- THE CHARTER (SAVED AND **VERIFIED**), THE BENCHMARK SEAL, AND THE ARCHITECTURAL READING THAT DETERMINED THE LANE MAP: **A 512-DIMENSIONAL STATE IS BEING USED TO SPECIFY ONE INTEGER** (2026-09-20 23:44, coordinator)

### The charter

Saved verbatim as `s31/BRIEF.md` — **1,660 lines, 47,348 bytes, existence and size verified by
reading the file back, not by asserting it.** S30-L0 made this same claim and the file was never
written; it stood false for 77 minutes and was the *fifth* instance in this project of prose naming
a path that does not exist — recorded, as it happens, by the very entry that committed it.
`s31/s31_verify.py` (lane D) asserts every path this ledger claims to have written.

One hard constraint (CVaR-VQE is the spine and the main scientific object), one endpoint
(**mean built-chain Cα RMSD on 126 targets, currently 3.2105 Å**), primary target < 3.00,
ambition < 2.50.

### Benchmark seal — verified, not assumed

`benchmark60` is sealed. `core/pipeline.py:386-388`: *"ITS SINGLE PRE-REGISTERED PASS IS SPENT — the
harness refuses to run it without `--i-am-spending-the-benchmark`, and it is never a tuning
instrument."* Two independent guard sites: `core/bench.py:1074` and `core/pipeline.py:1664`. The
working instrument this sprint is **`tuning126`**. Folds and clusters are pinned and will not be
regenerated.

### THE ARITHMETIC (carried from S30, unchanged — the mean is a tail statistic)

```
production, n = 126        mean 3.2105   median 2.9661
the worst 18 targets       mean 6.2758
the other 108              mean 2.6997

cap the worst 10 at 3.00 A  ->  2.9074  (-0.3031)   BEATS the primary target
cap the worst 18 at 3.00 A  ->  2.7426  (-0.4680)
cap the worst 30 at 3.00 A  ->  2.5778  (-0.6328)
```

And the tail is **selection-limited, not pool-limited**: on the genuinely worst 18 the ORACLE best
pool member is **2.5298 Å** with 11 of 18 under 3.00. The material is already in the candidate sets.

### THE ARCHITECTURAL READING, made from the code before any lane was briefed

`core/pipeline.py` ~838:

```python
E = _zrank(np.asarray(pool["sc"], float)[o])     # DIAGONAL: rank-standardised candidate scores
p, cvar, H, _circ = qm.run_cvar_vqe(E, alpha, T, n=cfg.vqe_qubits, ...)
block = np.asarray(Pt, float)[:dim, :dim]        # PAIRWISE distances between candidates
local = consensus_medoid(block, p)               # readout: p-weighted 1-median over `block`
```

**Two structural facts follow, and between them they set the sprint:**

1. **Theorem T1 (S30): for a diagonal `H`, the CVaR tail is always a prefix of the order induced by
   the objective.** The optimised state's entire structural output is therefore **one integer** —
   the prefix length `m`.

   > **A 9-qubit, 512-dimensional state is being used to specify a single number.** That is the
   > information-allocation defect at the centre of the architecture, and it is not a defect of
   > expressivity, optimisation or circuit depth.

2. **The Hamiltonian and the readout are decoupled.** `H` encodes *how good is each candidate
   alone* (diagonal, scores). The readout consumes `block` — *how structurally similar are
   candidates to each other* (off-diagonal). **The optimiser never sees the geometry it will be
   read out through.** This matters because the project's one established in-band discriminator is
   consensus/typicality (S12: score-filter + consensus medoid, −0.172 Å, the project's first CI
   excluding zero), and consensus is **irreducibly pairwise** — it cannot be written as a diagonal
   function of candidate index without precomputing it, which collapses it back to a diagonal.

And the fact that makes the whole thing worth attacking rather than abandoning: the same circuit
family reaches **0.2516 Å** built chain under an ORACLE objective and **3.4330 Å** under the
deployed native-free one, against a classical average's 3.2071. **The circuit can express excellent
states. The objective does not point at them.**

### THE LANE MAP — five launched, three slots held

| lane | remit | first question |
|---|---|---|
| **A** | quantum/CVaR theory | Characterise the objectives whose CVaR argmin is *not* a prefix. Which hypothesis of T1 fails off-diagonal — and does a generalised prefix theorem still apply? Then: is `H = diag(zrank) − λ·W(block)` physically meaningful, and is it classically reducible (charter §11)? |
| **B** | physical model | §7A free energy and §7B torsion. **Theoretical pre-check first:** `F = E − TS`; the enthalpic half is already closed (AMBER puts the native at the 51st percentile and is *worse* than chance as a ranker), so the novelty must live in `S` — and `S` must have a target argument that (sequence, pool) does not already carry |
| **C** | encoding / readout | §7C native-free sparse support; §12 candidate-index allocation; and the **128→512 widening's undischarged circularity gate**, which is cheap and either revives or closes a direction |
| **D** | verification / integrity | The four defects: the withdrawn positive in shipped code at `core/pipeline.py:821`; **the unpinned projection seed** (gates every sub-0.01 Å claim); the launcher/governor contradiction; and the verifier + sprint-wide multiplicity register |
| **E** | the incoherence hypothesis | See below — the sprint's sharpest new idea |

### THE NEW HYPOTHESIS, AND WHY I THINK IT IS THE SPRINT'S BEST SHOT (lane E)

S30's deepest result is that **what can be predicted is coherent with the pool's common mode and
therefore harmful**: at identical out-of-fold R² = 0.2355 a fitted corrector emits **+0.0554 Å
worse** while a synthetic i.i.d. one emits **−0.2466 Å better**. Every corrector the project owns
*raises* the residual's coherence with the common mode (0.6931 → 0.786 / 0.783 / 0.917).

S30 tested *"fit a corrector and apply it."* **It never tested "apply only the component orthogonal
to the common mode."**

And the enabling observation, which I believe S30 missed: since `mu = pool75_mean − d_nat` and the
distogram's `expected` is an estimate of `d_nat`, the quantity

```
mu_hat = pool75_mean - expected
```

**is computable at inference with no native at all** — and both terms are already features in
`s30_P_lr.py`. This is the project's standing "prediction/pool disagreement is a native-free signal"
arriving where it can actually be spent.

**My registered prior: roughly 2:1 against a deployable gain.** If the fitted corrector's residual
has coh 0.9172 it is *almost entirely* in the `mu` direction, so the orthogonal complement may be
nearly pure noise — giving the i.i.d. arm's *magnitude* without its *information*. **That is still a
clean result**, because it would show the i.i.d. benefit comes from genuine orthogonal signal rather
than from incoherence per se, which sharpens what a fourth observable must supply.

### WHAT I AM CARRYING FORWARD AS BINDING (until a lane breaks it)

- **Basis discipline.** Built chain **3.2105**; CA cloud **3.0483**; *set mean* 3.5507. Three
  different objects. The measured cloud→chain transfer for a correction is **0.92** — the 1.16
  coefficient in `operator-consumes-set-mean` maps *set mean* to output and **must not** be used to
  transfer a cloud delta (its own inputs give 4.21 Å against production's actual 3.21).
- **MDE = 2.8016 × SE, per comparison.** Below 0.7× is not a result; 0.7–1.0× is NOT MEASURED.
- **G1 (S30):** every achiral rotation/translation-invariant single-structure observable is a
  distance-map reading. Renaming one is not a new channel.
- **ORACLE labels travel inside the sentence carrying the number**, not a paragraph later.

---

## S31-L1 -- **R1, THE READOUT CEILING: THE ENTIRE QUANTUM STAGE CARRIES AT MOST k BITS AND CANNOT EMIT ANYTHING OUTSIDE THE POOL -- AND THIS IS A PROPERTY OF THE READOUT, NOT THE HAMILTONIAN** (2026-09-20 23:50, coordinator)

Derived from the code, not from prose, after the lanes were briefed. Handed to A (verify + measure
realised capacity), C (design a readout that escapes it) and F (the medoid's ceiling). **Posted
before any of them reported, so it is falsifiable by them rather than confirmed by them.**

### The statement

The quantum stage's output is exactly

```
v = P @ (p / sum(p));   return argmin_i v_i         # core/pipeline.py:795-803, :869-871
```

where `P` is the **pairwise Kabsch CA-RMSD matrix** between the top candidates
(`Pt[a] = aud.kabsch_rmsd_batch(W64top, W64top[a])`, `core/pipeline.py:771-777`) and `p` is the
state's measurement distribution. Four consequences:

1. **The output is always one of the `2^k` pool members.** The stage cannot construct a structure;
   it can only *name* a deposited one. `dim = min(1 << vqe_qubits, len(top))` with a hard `raise`
   if the pool is smaller, so `k` qubits means exactly `2^k` candidates.
2. **The state enters only through a linear map followed by an argmin.** `p` influences the answer
   solely through which cell of the hyperplane arrangement `{(P_i - P_j) . p = 0}` it lands in --
   a piecewise-constant map from the `(2^k - 1)`-simplex onto at most `2^k` outcomes.
3. **Therefore the whole quantum stage carries at most `k` bits.** Seven at the deployed
   `vqe_qubits = 7`.
4. All `2^k` outcomes are reachable: at a vertex `p = e_j`, `(P e_j)_i = P[i,j]` and `P[j,j] = 0`,
   so the argmin is `j`.

### Why it matters more than T1

T1 (S30) says a **diagonal** `H` makes the CVaR tail a prefix, so the state specifies one integer.
R1 says something strictly stronger and **Hamiltonian-independent**:

> **Even a perfect non-diagonal Hamiltonian cannot make this stage express more than `k` bits, and
> cannot make it emit anything outside the pool.** The cap lives in the readout.

This reframes the sprint. The charter's §6 asks what Hamiltonian the CVaR-VQE should optimise; R1
says that question is **downstream of a readout question nobody had stated**. It also explains the
Hamiltonian/readout decoupling recorded in S31-L0 from the other side: the readout already consumes
the pairwise matrix `P` and is the thing actually choosing, so a Hamiltonian built to "see the
pairwise geometry" is competing with the operator that consumes it.

### What it does NOT say

**It is a cap on information, not on value.** ORACLE argmin over 128 is **2.1435 A on the CA point
cloud** against production's 3.0483 -- roughly **0.90 A of genuine headroom**, so this readout is
*not* capped below production. **ORACLE / NOT DEPLOYABLE**, cloud basis; measured cloud->chain
transfer for a correction is 0.92.

It also does not say index redesign is worthless: re-labelling cannot raise capacity above `k`
bits, but it can make a low-entropy or partially-measured state land in a better cell, and make
each qubit correspond to a meaningful distinction. It *does* say that anyone hoping index redesign
buys **information** should stop.

### The open sub-question, which is sharper than the theorem

How many outcomes does the deployed stage **actually** produce? R1 bounds capacity at `k` bits;
the *realised* capacity could be far lower if the ansatz cannot reach the cells. Lane A owns this.
**A realised capacity strictly below `k` bits would be a worse and more interesting result than R1
itself.**

### Falsifier

R1 fails if `argmin(P p)` over the reachable `p` is not confined to the index set -- i.e. if any
code path lets the stage emit a structure that is not a pool member -- or if point 2's arrangement
argument is wrong about the reachable set. Both are checkable in an hour and lane A was asked to
attack them.

---

## S31-L2 -- **THERE ARE TWO QUANTUM READOUTS AND THE DEPLOYABLE PATH USES THE WEAK ONE.** THE CHARTER'S "THE CIRCUIT CAN EXPRESS GOOD SOLUTIONS" IS ABOUT THE OTHER ONE -- AND THE REAL FINDING IS THE INVERSE: **EXPRESSIVITY WITHOUT AN ALIGNED OBJECTIVE IS HARMFUL** (2026-09-20 23:53, coordinator)

Read from the code. Handed to C (whose remit it reframes) and A (whose premise it qualifies).

### The two readouts

**Readout 1 -- deployable** (`core/pipeline.py:795-803`, used by `quantum_stage`):

```python
v = Dblock @ (w / sum(w));  return int(np.argmin(v))     # w = measurement PROBABILITIES
```

Returns **one pool member**. **7 bits** at the deployed `vqe_qubits = 7` (R1, S31-L1).
ORACLE ceiling, argmin over 128: **2.1435 A** (CA cloud).

**Readout 2 -- every ceiling experiment** (`s27/s28_A_amp.py:105-117`):

```python
w = psi / psi.sum();  C = (w @ frame.Wf).reshape(n, 3)   # psi = real AMPLITUDES
```

`psi` are **amplitudes, not probabilities**, so `w` sums to 1 but **individual weights may be
negative**. This is an **affine** combination and it can leave the candidates' convex hull. Hence
`circ_best` = **0.2516 A** built chain while the best single pool member averages 1.7108 A.

### What this does to the charter's §5E

The charter concludes, from 0.2516 A (ORACLE objective) against 3.4330 A (deployed objective),
that *"the circuit can express good solutions; the objective does not point at them."* **Both
numbers use readout 2** (`circ_opt` is `circ_l1_i80` through the same amplitude readout), so the
comparison is fair *on readout* and **the conclusion holds -- for readout 2.** It says nothing
about readout 1, which is what a deployable quantum stage would use. **The report must not inherit
that conflation.**

### AND THE FINDING, WHICH IS THE INVERSE OF THE FRAMING

Readout 2 is enormously more expressive than readout 1. Under the deployed objective it emits
**3.4330 A on the built chain -- WORSE than the uniform average's 3.2071**, and the uniform average
is simply the **heavily regularised special case** of the same affine readout (all weights `1/k`,
no negatives).

> **Expressivity without an aligned objective is not neutral, it is harmful.** More degrees of
> freedom let a coherent-error objective do more damage. This is S30's coherence result --
> *what can be predicted is coherent and therefore harmful* -- arriving at the readout.

**Corollary, and it is why this is useful rather than merely deflating:** a sparse or norm-bounded
readout is not a bit-budget compromise, it is **the regulariser standing between a degenerate
affine class and a usable one.** That converts charter §7C from *"can sparsity buy capacity"* into
*"what is the right constraint set on a readout that is otherwise too expressive to be safe"* --
and makes non-negativity, an `ess` floor, an `||w||_1` bound and group sparsity **principled**
rather than heuristic.

### Two cautions handed to lane C with the work

1. **Readout 2's ORACLE ceiling is probably vacuous.** For n = 9-16 the coordinate space is
   3n = 27-48 dimensional and the affine hull of 128-500 candidates generically **spans it** --
   `oracle_affine_ls` reportedly reaches the native to 1e-7, which is what a degenerate ceiling
   looks like. Coordinate-space stable rank is ~3.4-3.6 with k90 = 11.2, so severe
   ill-conditioning is expected: reaching the native should require enormous cancelling weights.
   `weight_diag` already returns `frac_neg`, `neg_mass`, `ess`. **Measure them.** If the ceiling
   needs a tiny `ess`, it is a numerical artefact of an ill-conditioned basis, not a reachable
   structure.
2. **The published ladder is non-monotonic because the candidate SET changes.** 2-of-75 with free
   weights (2.1683 A, 11.4+ bits) is *worse* than argmin-over-128 (2.1435 A, 7 bits) -- but 75 != 128.
   **Hold the candidate set fixed across the ladder or the comparison is not about the readout.**
   This alone may explain S30's "argmin dominates at every bit budget".

### The deficit this sprint actually has to close, stated plainly

`quantum: bool = False` in `PROD = Config()` (`core/pipeline.py:179, :241`): **the 3.2105 A endpoint
is produced with the quantum stage OFF.** Turning it on, as built, costs **+0.2260 A on the built
chain** (`circ_opt` 3.4330 against PROD 3.2071 in the S30 meter cache). So a quantum-spine
architecture must recover **0.226 A to reach parity with the classical path it replaces**, and then
**0.211 A more** to take 3.2105 below 3.00 -- **~0.437 A in total**, not the 0.21 A that
"3.21 -> 3.00" suggests.

## S31-L3 -- **THE COMMON-MODE DIRECTION IS NOT NATIVE-FREE ESTIMABLE, AND THE CLOSED FORM SAYS EXACTLY WHY**: `mu_hat = mu - y` EXACTLY, SO THE ESTIMATOR'S ERROR *IS* THE PRIOR ERROR IT EXISTS TO PREDICT -- AND THE FITTED CORRECTOR TURNS OUT TO PREDICT THE COMPONENT **ORTHOGONAL** TO THE COMMON MODE, NOT THE COMMON MODE (2026-09-20 23:57, E)

Registered in `s31/PREREG_S31_E.md`, committed at `593bdd2e` **before the first number existed**;
AMENDMENT 1 at `f246eaba`, timestamped and marked as decided *after* seeing E1.
Artefacts: `s31/results/s31_E1_direction.json`, `s31/results/s31_E2_deltas.json`,
`s31/results/s31_E3_decomp.json`. Code: `s31/s31_E_lib.py`, `s31/s31_E1_direction.py`,
`s31/s31_E2_deltas.py`, `s31/s31_E3_decomp.py`.

### The closed form, which is the durable output

With `y = expected - d_nat` (prior error, ORACLE), `mu = pool75_mean - d_nat` (the pool's
common-mode pair error, ORACLE, the S30-L7 quantity) and `mu_hat = pool75_mean - expected`
(**native-free**, and the trace is *proved in code*, not asserted -- `mu_hat` is bit-identical when
`nat_ca` is replaced by NaN, `s31_E_lib.native_free_trace`):

```
mu_hat = mu - y          EXACTLY.   max |dev| = 1.8e-15 over all 8549 pairs, 126/126 targets

=>  corr(mu_hat, mu) = (1 - coh0 * r) / sqrt(1 + r^2 - 2 * coh0 * r),    r = sd(y)/sd(mu)
    reproduces the measured per-target correlation to 1.1e-15.
```

**`mu_hat` is not a noisy estimate of `mu`. Its error is exactly `-y`** -- the quantity the corrector
exists to predict. To use the estimator you would already need the answer.

### MEASURED: r = 1.5997, fold95 [1.4421, 1.7549], median 1.305, **100% of targets above 1**

The prior's within-target error dispersion is **1.6x** the pool's common-mode dispersion, so the
"estimate" is dominated by its own error. At those values the closed form is *negative*; the small
positive number we measure is Jensen curvature across targets, not signal.

**Not a contradiction of `pool-error-is-68-percent-common-mode`.** That 68% splits *pool members'*
coordinate error into a shared bias and an idiosyncratic part. `r` compares the **distogram's** error
against the pool's shared part. Different objects; both can be true and are.

### E1, the registered gate -- DEAD on both clauses

Within-target cosine(`mu_hat`, `mu`), 126 targets, all pairs `min_sep = 2`:

```
cos                 mean +0.0495  median +0.0614  fold95 [+0.0143, +0.0779]
cos (perm control)  mean +0.0213      <- mu_hat permuted within target: same marginal, same norm
excess over control       +0.0282  0.51x MDE  fold95 [-0.0120, +0.0820]  2/5 folds  NOT MEASURED
centered corr       mean +0.0408; excess over the same control 0.50x MDE, also NOT MEASURED
long range (sep>=7, 106 targets)  cos +0.0205; excess 1.20x MDE but the LEVEL is 0.02
only 1% of targets reach the registered USABLE bar of 0.577
```

Registered bars were STRONG 0.707 / USABLE 0.577 / WEAK 0.30 / **DEAD** below 0.30 *or* not
separated from the matched control at 1.0x MDE. It fails **both**.

### The mechanism finding, which corrects S30 §12 and contract rule 29

S30's `N3_plus_pool` corrector is reproduced **bit-exactly** here (R² 0.235453 vs S30's 0.235453;
`coh` 0.693096 vs S30's 0.693096), so this is the same object, not a re-implementation.

```
within-target corr(yhat, mu)       -0.0295      pooled -0.0058     <- essentially UNCORRELATED
within-target corr(yhat, mu_hat)   -0.9329      pooled -0.9078     <- it IS the pool-disagreement feature
within-target corr(y - yhat, mu)   +0.9172      pooled +0.9658     <- S30's number, reproduced

orthogonal SSE decomposition against mu (ORACLE / NOT DEPLOYABLE):
  energy of y      54.70% along mu / 45.30% perp   (per-target mean; POOLED 71.10 / 28.90)
  the corrector's sum-of-squares reduction:
        79.84% from the PERP component  (per-target mean)
       100.79% from the PERP component, -0.79% along  (POOLED)
  yhat's own energy along mu:  8.94% mean, 4.12% median
```

> **S30 §12 says "the predictable part of the prior's error *is* the common mode." Measured on
> S30's own corrector, the opposite is true: essentially all of what it predicts is ORTHOGONAL to
> the common mode, and it leaves the common mode untouched.** `coh` rises 0.6931 -> 0.9172 for the
> *opposite* reason to the one published -- not because the corrector captures the common mode, but
> because it strips everything *except* the common mode, leaving a residual that is nearly pure
> common mode. The observation S30 reported is exact; its explanation is inverted.

The reason is algebra again: the corrector's dominant feature is `exp_minus_pool75 = -mu_hat`, and
`mu_hat = mu - y` contains `-y`. Regressing `y` on `mu - y` recovers `y` through the `-y` term. That
is legitimate (native-free at inference) but it means the pool's "information about the prior's
error", S30's +0.0758 out-of-fold at long range, is the **non-common-mode** part of that error.

**Contract rule 29** ("correlating with the error is not the test; incoherence with the pool's
common mode is") rests on S30's inverted mechanism and should be re-read in light of this. The
*measurement* behind it -- +0.0554 coherent vs -0.2466 i.i.d. at matched R² -- stands untouched.

### Status of the lane's own prediction

The coordinator registered 2:1 against E2/E3 and I registered 4:1 / 3:1, both on the reasoning that
the corrector's residual sits at `coh` 0.9172 so its orthogonal complement would be noise. **The
direction was right and the mechanism was wrong**: the failure is one level earlier, at `mu_hat` not
being an estimate of `mu` at all, and the corrector is only 8.9% along `mu` rather than "almost
entirely" along it. A prediction that got the sign right for the wrong reason is worth less than a
correct one, and it is recorded that way.

### What is closed, what is open, and the specification this leaves

**CLOSED:** "the common-mode direction is native-free estimable as `pool75_mean - expected`". It is
*computable* and it is not an *estimate*. S30's `coh` admission test keeps its ORACLE label and loses
its only proposed native-free surrogate.

**THE SPECIFICATION, which is the point of the entry:** the closed form says exactly what a future
common-mode estimator must satisfy -- **its error must be small relative to the common mode itself,
and the distogram's error is 1.6x too large.** That is a quantitative target, not an exhausted search.

**OPEN, and running now:** the ORACLE class ceiling registered in AMENDMENT 1 -- over all corrections
orthogonal to the true `mu`, what is the best achievable endpoint? Early indications from the
orthogonal decomposition are that the prize is **along** the common mode, which would invert the
sprint's open question rather than answer it. Endpoint arms pending; nothing is claimed from them here.
