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

> **ANNOTATION (2026-09-20 23:58, coordinator, on lane C's refutation) — R1's THEOREM IS RIGHT AND ITS
> QUANTIFIER IS WRONG. MY OWN FALSIFIER FIRED, IN SHIPPED CODE.**
>
> I wrote the falsifier as *"R1 fails if any code path lets the stage emit a structure that is not
> a pool member."* Lane C took it seriously and found that path **in `core/pipeline.py`**:
>
> ```python
> def average_weighted(Wo, block, w, clk):            # core/pipeline.py:880-895
>     b   = consensus_medoid(block, w)                #   medoid only as the superposition FRAME
>     Sup = cc.superpose_batch(Wo, Wo[b])
>     ww  = w / w.sum()
>     C   = np.tensordot(ww, Sup, axes=(0, 0))        #   a NEW structure, generically not a member
> ```
>
> called at `:1110` (p-weighted) and `:1116` (uniform), scored at `:1179-1182`, registered as
> comparison arms at `:1611` and `:1613`. **`p` enters continuously.** So there are **three**
> readouts, not the one R1 describes and not the two of S31-L2:
>
> | | where | what `p` does | reachable set | capacity |
> |---|---|---|---|---|
> | **R1-sel** | `core/pipeline.py:795-803` | picks a cell of an argmin arrangement | the `2^k` pool members | **<= k bits — R1 is exactly right here** |
> | **R1-convex** | `core/pipeline.py:880-895` **SHIPPED** | enters continuously as weights | the **convex hull** of the posed candidates | not k bits |
> | **R2-affine** | `s27/s28_A_amp.py:105-117`, harness only | signed amplitudes | **outside** the hull | not k bits |
>
> **Points 1 and 3 hold for R1-sel only.** The sentence *"the entire quantum stage carries at most
> k bits"* is false as written, and S31-L2 compounded it by merging the convex and affine readouts —
> a distinction that is load-bearing, because **the convex one is shipped and the affine one is
> not**, and they have different ceilings and opposite verdicts. The original wording of both
> entries stands above.
>
> **What survives, and it is still the useful part:** the *selection* readout — the one the charter's
> §6 Hamiltonian question is really about — is capped at k bits and cannot leave the pool, and that
> cap is a property of the readout rather than of any Hamiltonian.

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

> **ANNOTATION (2026-09-20 23:58, coordinator, on lane C's refutation) — "TWO READOUTS" IS WRONG; THERE ARE
> THREE, AND THE ONE I OMITTED IS THE SHIPPED CONTINUOUS ONE.** See the annotation on S31-L1. The
> merge of *convex* (`core/pipeline.py:880-895`, shipped, `p >= 0`, reachable set = the convex hull)
> with *affine* (`s27/s28_A_amp.py`, harness only, signed amplitudes, leaves the hull) is the error;
> they have different ceilings and opposite verdicts. Original wording stands.
>
> **And one rung of the ladder I proposed already has a measured answer, which lane C supplied:**
> the convex optimum under the deployed objective, started from production and converged, is
> **3.0522 A (CA cloud) — i.e. the best convex reweighting the objective can find over that set IS
> the uniform average it already emits.** Two independent constructions agree (S23-L8's
> probability-weighted readout at +0.0213, 0.65x MDE, a NULL; S23-L5's leakage-ORACLE grid selecting
> the uniform point as optimal). **The convex rung is closed by ceiling, not by a failed fit** —
> nobody should re-derive it. The open rung is the **norm- or sparsity-bounded middle.**

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

---

## S31-L4 -- **THE DEPLOYED CVaR FREE ENERGY IS A CONVEX PROGRAM WITH A CLOSED-FORM GLOBAL MINIMISER.** THE SHIPPED CIRCUIT IS STRICTLY WORSE THAN IT IN 12/12 CELLS, AND THE GAP IS AN **EXPRESSIVITY** FLOOR, NOT AN OPTIMISER ONE -- SO CHARTER §11 IS CLOSED FOR THE DEPLOYED OBJECTIVE (2026-09-21 00:01, L)

Full working, tables and caveats: `s31/LIT_L.md` §L1.1-L1.4. Scripts:
`scratchpad/lit_L_dequant_check.py`, `scratchpad/lit_L_gap_budget.py`.

### THEOREM (derived here from Rockafellar-Uryasev 2000 + the Gibbs variational principle)

For `F(p) = CVaR_alpha(E;p) - T*H(p)` -- the objective at `core/quantum.py:993` --
the Rockafellar-Uryasev lower-tail form `CVaR^low = max_s { s - (1/alpha) sum_i p_i (s-E_i)_+ }`
is **affine in p** inside the max, so `F` is **convex in p**, strictly for `T>0`, and concave
in `s`. Sion's minimax swaps the order:

```
min_p F = max_s { s - T*log sum_i exp( (s-E_i)_+ / (alpha*T) ) }      [1-D, concave]
p*_i   proportional to   exp( (s* - E_i)_+ / (alpha*T) )
```

`p*` is a **hinged Gibbs distribution**: uniform on every candidate at or above the VaR level
`s*`, exponentially tilted only below it. **The whole 2**n-dimensional optimisation is pinned
by one scalar, and that scalar has a closed form.** This is strictly stronger than T1: the
entropy term was added to break T1's degeneracy, and the result is still one number.

### OURS -- verified against the shipped code, 12 cells, n in {7,9}, alpha in {0.1,0.25,1.0}, T in {0.1,0.05}

- **Strong duality holds to 1e-9 ... 1e-16 in all 12 cells.** An independent mirror descent
  over the full simplex reproduces the same optimum to ~1e-5. The closed form IS the optimum.
- **`run_cvar_vqe` at its deployed settings is strictly worse in 12/12**, gap +0.0186 to
  +0.2368 in F units, total-variation distance **0.246 to 0.962** from `p*`.
- The closed form is always the **more entropic** distribution (n=9, alpha=1, T=0.1:
  `H* = 3.693` bits vs the circuit's `1.181`).

### OURS -- the gap is EXPRESSIVITY, not optimisation

80 -> 2000 Adam iterations at fixed depth moves the gap by **nothing** (0.025764 -> 0.025944
at alpha=0.1; 0.233580 -> 0.233116 at alpha=1.0). Only **depth** and **restarts** close it,
and even at layers=12/400 iters/8 restarts (~60 s) the gap is +0.0042/+0.0176 with TV ~0.10 --
against **microseconds** for the closed form. 21 parameters (`n*layers`) cannot cover a
127-dimensional simplex.

### THE CORRECTION THIS FORCES ON OUR OWN SHIPPED DOCSTRING

`core/quantum.py:997-1002` states the entropy collapse "**is a property of CVaR, not of the
optimiser**."

- At **T = 0** this is CORRECT, and now has a one-line proof: any `p` with mass `>= alpha` on
  the argmin attains `CVaR = E_min` exactly, so the minimiser set is a **positive-volume flat
  face** and the optimiser's path -- not the objective -- picks the point in it. **Anything the
  readout reads at T=0 is an Adam artefact.**
- At the **deployed T = 0.1** it is **WRONG**. The objective's own optimum at `alpha=1, T=0.1`
  carries **2.310 bits** (n=7); the circuit delivers **0.671**. The collapse is a property of
  **the ansatz** -- neither CVaR nor the optimiser.

Re-pricing the inherited alpha reading (`pipeline.py:825`, 0.076 -> 6.36 bits): on a shuffled
z-rank ladder the circuit spans 3.88 bits while the **objective** asks for 2.375, so ~60% of the
alpha-effect on state entropy is real and ~40% is the ansatz amplifying it. **CAVEAT CARRIED
WITH THE NUMBER: synthetic ladder at n=7, not our measured pool `E`. Indicative, not measured
on the instrument.**

### WHAT IS CLOSED, AND THE ONE THING THAT IS NOT

**CLOSED -- charter §11 for the deployed objective.** The shipped CVaR-VQE is a *lossy
approximate solver for a convex program with an analytic solution*. No interference, no
spectrum, no entanglement is doing anything. Whatever it contributes, it contributes **by
failing to optimise**.

**NOT CLOSED.** The circuit's reachable set is a 21-parameter manifold inside the simplex; that
constraint is an inductive bias, and S20's law says optimising a bad objective harder makes
things worse. So the failure may be *useful* -- but if it is, the bias is a **Born machine**
`|psi|^2` from a shallow 1-D RY+CNOT circuit, i.e. an **MPS Born machine** (Han et al., PRX
8:031012, 2018), which is a classical tensor-network model. **Even the escape hatch is
classical.**

### THE DECISIVE EXPERIMENT, WHICH COSTS NOTHING (handed to lanes A and C)

Substitute the closed-form `p*` for `run_cvar_vqe`'s `p` and rerun tuning126.
- endpoint **improves** -> the quantum layer is a softmax plus a root-find;
- endpoint **worsens** -> the circuit's *inability* to optimise is the active ingredient, which
  is S20 arriving from the objective side and is a real publishable negative.

There is no third outcome. **CAVEAT: the readout is `consensus_medoid(block, p)`, so a large TV
in `p` need not move the medoid. Score at the ENDPOINT, never on `p`.**

### TWO FURTHER CLOSURES FROM THE LITERATURE (details in `LIT_L.md`)

1. **`H = diag(zrank) - lambda*W(block)` is not a CVaR-VQE and cannot be made into one.** CVaR
   needs an energy per *shot*; only a diagonal `H` gives every bitstring a definite eigenvalue
   (Barkoutsos et al., *Quantum* 4, 256). Off-diagonal terms cannot be assigned to a single
   outcome. **And dimension counting kills it independently:** a candidate-index register has
   Hilbert dimension = number of candidates, so *any* operator on it is a 128x128 or 512x512
   matrix and its spectrum is a microsecond `eigh`. **No Hamiltonian on a candidate-index
   register can be classically hard.** The clustering literature that does this properly
   (*Front. Phys.* 2025, arXiv:2502.06542) uses **one qubit per data point** -- 128 qubits for
   our pool, not 7 -- and demonstrates **no** quantum advantage at 175 points on D-Wave.
2. **Our `Var ~ 16/D = 16*2^(-n)` gradient law is the predicted global-cost barren plateau**
   (Cerezo et al., *Nat. Commun.* 12, 1791, 2021: global costs vanish exponentially **even at
   O(1) depth**; depth does not fix it). Quantitatively this prices lane C's widening at **4x
   less gradient signal per two added qubits**: n=7 -> 0.125, n=9 -> 0.031, n=12 -> 0.0039.
   128->512 is affordable; past ~12 qubits is the wall. **CAVEAT THAT CUTS THE OTHER WAY: at
   n=7-9 we are NOT yet gradient-limited** -- depth and restarts still move the objective a
   lot, so this is a bound on where the architecture can go, not a diagnosis of where it is
   stuck now.

### ALPHA SCHEDULING -- the literature's alpha and ours are different objects

Kolotouros & Wallden (*PRResearch* 4, 023225): Ascending-CVaR, `alpha_{t+1} = alpha_t + lambda`,
`lambda in [0.025,0.045]`. Their Prop. 1 -- **all `CVaR_alpha` with `alpha <= kappa` share the
same ground state** -- holds at `T=0`, so their alpha is a **landscape homotopy**, and the
schedule runs **upward to alpha=1**. At our `T=0.1` the optimum is genuinely alpha-dependent
(`H*` 4.685 vs 2.310 bits), so our alpha is a **distribution-shape** knob and their endpoint is
our *narrowest* ensemble -- the direction the consensus readout least wants. **The usable
translation**: schedule alpha upward while solving `T(alpha)` to hold `H*` fixed, which the
closed form gives analytically. Ranked **below** the `p*` substitution, because it fixes the
landscape and L1.1 says the landscape is not the binding constraint.

## S31-L5 -- **THE 128->512 WIDENING SURVIVES ITS CIRCULARITY GATE BUT IS RE-PRICED FROM -1.9004 TO -1.1811 A**, AND THE CIRCULARITY IS A **LITERAL IDENTITY**: THE 18 TARGETS WITH THE LARGEST WIDENING GAIN FROM THE TOP-75 **ARE** `FAIL18`, 18 OF 18. ON THE CLEAN FILTER-INDEPENDENT TAIL THE ORACLE EFFECT IS **-1.1811 A CLOUD / -1.1762 A CHAIN, p = 0.0001** AGAINST ITS OWN RANDOM-18 NULL -- AND A LEVEL CONTROL REMOVES TWO THIRDS OF `FAIL18`'s EXCESS (2026-09-21 00:02, C)

Pre-registration `s31/PREREG_S31_C.md`, committed at **8ce5e1a0** (swept into the coordinator's R1
commit) **before the first lane-C number existed**; amendment 1 at **2e12e02c**, also before.
Code `s31/s31_C_cache.py` (shared pool table, 126 npz) and `s31/s31_C_widen.py`.
Artefacts `s31/results/s31_C_widen.json` and `s31/results/s31_C_widen_rows.jsonl` (126 rows).

**EVERY NUMBER IN THIS ENTRY IS ORACLE AND NOT DEPLOYABLE.** `best1(N)` is the minimum native
RMSD over the first `N` of the shipped DIS order: it says the good candidate **is there**, not that
anything native-free can find it. S29 closed recognition three ways and S30 confirmed it.

## WHAT WAS ASKED

S30-L11 published, as an incidental finding, that widening the quantum register 128 -> 512 is worth
**-1.9004 A on FAIL18** against -0.1907 on the other 108 (2.77x MDE), and the S30 report carried it
forward as **"OPEN, with a known circularity the sprint never discharged"**: `FAIL18` is defined by
the filter's own recall, so *"the good candidate is outside the window"* may be **produced by** the
thing it is offered as evidence about. The coordinator called the filter-independent-tail check
"the gate on that whole direction". No S30 ledger entry ran it. Lane F ran the analogous check for
the **averaging** operator (S30-L16) and found the effect absent or reversed; that is a different
operator and its own entry says so.

## THE CIRCULARITY IS NOT A RISK, IT IS AN IDENTITY

The pre-registration defined a **definition-matched stratum** `defn18` -- the 18 targets with the
largest `best1(75) - best1(500)`, i.e. selected explicitly for the quantity `FAIL18`'s rule
thresholds -- to price how much of the published figure is forced by its own construction.

```
defn18 INTERSECT FAIL18  =  18 of 18          (s31_C_widen.json: overlap["FAIL18&defn18"])
defn18 effect            =  -1.9004 A          identical to FAIL18's, to four decimals
```

**`FAIL18` *is* the top-18 by the widening gain from the top-75.** Choosing the 18 targets whose
best candidate hides below rank 75 and then reporting that their best candidate also hides below
rank 128 is one statement, not two. The published -1.9004 is a near-tautology and **is not
quotable again.**

The rank diagnostic says the same thing from the other side. Under a score with no within-pool
skill the ORACLE-best member's rank is uniform on 1..500, so `P(rank > 128) = 0.744`:

```
stratum              frac(best member below rank 128)
other 108                       0.417        the score genuinely pulls the answer into the window
worst18_poolmean                0.667        near the uninformative null
worst18_bestpool                0.722        near the uninformative null
FAIL18                          1.000        ABOVE the null; p = 0.744^18 = 0.005 by chance
```

`FAIL18` is the only stratum **worse than an uninformative score**, which cannot happen by sampling
and is the definition showing through.

## THE GATE'S ANSWER: THE DIRECTION SURVIVES, AT 62% OF THE PUBLISHED SIZE

Effect of widening = `best1(500) - best1(128)`, negative = widening helps. **ORACLE / NOT
DEPLOYABLE.** CA point cloud; the chain column is s29 lane O's already-built chains, nothing was
rebuilt.

```
stratum               n   effect(cloud)  effect(chain)   best1_128   best1_500   p vs random-18
FAIL18               18      -1.9004        -1.8996        4.1846      2.2842       0.0000
defn18 (= FAIL18)    18      -1.9004        -1.8996        4.1846      2.2842       0.0000
worst18_poolmean     18      -1.1811        -1.1762        3.4097      2.2286       0.0001   <- headline
worst18_bestpool     18      -0.8260        -0.8260        3.9190      3.0930       0.0168
other 108           108      -0.1907        -0.1917        1.8060      1.6153          --
all 126             126      -0.4350        -0.4357        2.1458      1.7108          --

random-18 null on the stratum mean: mean -0.4363, sd 0.1652, 2.5th percentile -0.7953
```

**F-C3a fires as registered**: the headline filter-independent stratum is at **-1.1811 A**, below
the registered -1.00 bar, and below the random-18 null's 2.5th percentile at **p = 0.0001**.
`worst18_bestpool` -- which the pre-registration flagged in advance as biased *against* the effect,
because `best1(500)` is the subtrahend -- clears its null at p = 0.0168 but not the -1.00 bar.

**The registered verdict is REVIVE, and the registered consequence is that the quotable number is
now -1.1811 A (cloud) / -1.1762 A (chain), ORACLE, not -1.9004.**

## THE LEVEL CONTROL, WHICH IS THE PART THAT SHOULD TRAVEL WITH THE NUMBER

The effect is mechanically bounded by how bad `best1(128)` already is. Regressing the per-target
effect on `best1(128)` over all 126:

```
slope -0.4356   R2 0.5564

residual (the part NOT explained by the level of best1_128)
  FAIL18            -0.5774
  defn18            -0.5774
  worst18_poolmean  -0.1956
  worst18_bestpool  +0.3813        <- positive: widening helps this stratum LESS than its level predicts
  other 108         +0.0962
```

**Once the level is removed, only -0.196 A of the clean tail's -1.181 A is unexplained, against
FAIL18's -0.577 A.** So roughly two thirds of `FAIL18`'s excess over the general level-effect is
the stratum definition, and the honest mechanism is not "hard targets hide their answer deeper" but
**"the score has no within-pool skill on hard targets, so the best member lands roughly uniformly
and a 128-window misses it ~74% of the time"** -- which is S30 section 4.1's within-pool rho
degradation (+0.645 on the easy 108 against +0.107 / +0.380 / +0.344 on the three tails) arriving
at the register.

## A SIGN ERROR I MADE AND CAUGHT BEFORE ANY NUMBER LEFT THE LANE

The first version of `s31_C_widen.py` defined the contrast as `best1(128) - best1(500)`, which is
non-negative by construction (a running minimum cannot rise), and then read the falsifiers and the
null's **2.5th** percentile as though more negative were better -- the wrong tail of the right null,
which is S31 contract rule 8's exact failure mode, committed by the lane that quoted rule 8 in its
own pre-registration. The first run therefore printed `F-C3b FIRES -- CLOSE`, the opposite of the
truth. Corrected in place with the original error stated in the file's own comment at
`s31/s31_C_widen.py:49-54`; no number from the wrong-tail run was reported anywhere.

## WHAT THIS DOES AND DOES NOT LICENCE

- It **does** relocate the tail's failure from "the pool does not contain the answer" to "the
  register does not contain the candidate", at a re-priced -1.18 A rather than -1.90 A, on a
  stratum that is not defined by the filter.
- It **does not** say two more qubits are worth 1.18 A. The measurement is the ORACLE argmin, and
  production does not achieve the argmin at *any* width -- production is 3.2105 on the chain where
  `best1_top128` is 2.1435. Widening only pays if a selector exists, and recognition is closed.
- It **does** compose with the set-matched readout ladder (S31-L6, next): the convex readout over
  the **same** 128 reaches 1.8538 A on the built chain where the argmin over those 128 reaches
  2.1435, so the register width and the readout class are separable levers and the second one is
  larger at the deployed width.

## COMPARISONS MADE (contract rule 23)

5 strata x 1 contrast on 2 bases = 10; 1 random-18 null (20,000 draws) reused for all strata;
2 diagnostics (rank fraction, definition-match); 1 regression level control; 2 fold-clustered
paired comparisons. **13 read as results, all pre-registered.** No grid, no per-target maximum,
no split-half arm needed.
