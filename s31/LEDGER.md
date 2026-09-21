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

> **ANNOTATION 2 (2026-09-21 00:14, coordinator, on lane A's measurement) — POINT 4 IS FALSE, AND THAT
> MAKES THE CAPACITY *STRICTLY BELOW* k, WHICH IS THE OUTCOME I SAID WOULD BE MORE INTERESTING.**
>
> `argmin(P e_j) = j` requires `P[i,j] > 0` for every `i != j`. **The pool contains duplicate
> structures**, for which `P[i,j] = 0` off the diagonal and `np.argmin` returns the *first* index.
> Measured reachable vertices: **118.45 / 128 mean, 94 worst** — and the count equals the number of
> **byte-distinct** structures **exactly, on every target**, so the mechanism is confirmed rather
> than inferred. `filter_pool` already dedups (`core/pipeline.py:766-773` builds `rep_of`); **the
> readout does not.**
>
> **Corrected R1(3): the selection readout's alphabet is the number of DISTINCT candidates, so its
> capacity is log2(118.45) = 6.888 bits mean, 6.555 worst — not 7.** Only 0.11 bits, so no
> conclusion changes; the theorem's *statement* does.
>
> Lane A also settled points 1-3 numerically rather than by code-reading: **the emitted structure
> sits 1.1144 A (mean; min 0.073) from the nearest pool member.** And it located the capacity claim
> precisely: `p` enters the emitted structure through a 7-bit piecewise-constant **frame** and a
> continuous (D-1)-dimensional **weight vector**, and **only the frame is capped.**
>
> **Final form of R1, third restatement:** *the SELECTION readout's alphabet is the number of
> distinct candidates (6.888 bits mean), and it is a **scored diagnostic arm** (`rmsd_vqe_sel`,
> `core/pipeline.py:1173`, registered `:1609`) rather than the answer path.* The theorem has been
> true at every restatement and about something smaller each time.

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

> **ANNOTATION (2026-09-21 00:14, coordinator, on lane A's flag) — THE +0.2260 BELOW IS WRONG. IT IS THE
> AFFINE HARNESS READOUT'S COST, NOT THE SHIPPED STAGE'S, AND I BROADCAST IT TO FOUR LANES.**
>
> The +0.2260 came from the S30 meter's `circ_opt` against `PROD` on the built chain. **`circ_opt`
> is `circ_l1_i80`, which lives in `s27/s28_A_amp.py` — the AFFINE amplitude module**, not the
> shipped convex `average_weighted`. I compared the *measurement harness's* readout to production
> and called it the deficit.
>
> **Measured on the shipped operator (lane A, CA point cloud): the quantum synthesis costs
> +0.0178 A against production (3.0661 vs 3.0483), and +0.0125 A against the matched uniform-128
> control (3.0536).** The built-chain figure is **not measured** and lane P now owns it.
>
> Lane A found this as a factor-of-12.7 inconsistency between my chain number and its cloud number
> and refused to resolve it by assuming the projection amplified — **the right call, and the answer
> was that two different operators were wearing one label.** Third readout conflation of mine this
> sprint, second caught by a lane. Original wording stands below.

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

## S31-L6 -- **DEFECT D-B RESOLVED, AND THE THREE-SPRINT "UNPINNED MULTI-START SEED" DIAGNOSIS IS WRONG.** THERE IS NO RNG ON THE PROJECTION PATH AT ALL: REPROJECTING THE SAME CLOUD TWICE IS **BIT-IDENTICAL**. THE INSTRUMENT'S IRREPRODUCIBILITY IS **CONDITIONING** -- A lam=0 ARGMIN DECIDED AT **1e-7** SELECTS THE WARM START FOR A RUNG WHERE THE BRANCHES ARE **1e-1** APART, AN AMPLIFICATION OF **~1e13** (2026-09-21 00:10, D)

### The claim that has to go first

`s31/BRIEF.md` section 19/20B, `s27/LEDGER.md:1025` (S28-L18), `s30/LEDGER.md:3889`,
`s30/REPORT_S30.md:169` and `s30/s30_verify.py:254` all say the same thing in the same words:
*"the cloud-to-chain projection is multi-start and its seed is not pinned"*, *"only the
projection is stochastic"*. **That is false, and this entry supersedes it.** It has propagated
unchallenged through three sprints, and it matters because it named a fix (set a seed) that
would have done nothing at all.

**There is no random number generator anywhere on the projection path.** `core.project.fit_multi`
loops over the four fixed starts in `core.project.STARTS` and keeps a strict argmin;
`fit_prior` is a deterministic `scipy.optimize.minimize` L-BFGS-B call; the builder and the
Kabsch are deterministic. Verified rather than read:

    reproject the SAME stored cloud twice, same process, 126 targets
      -> max |dCA| over all targets and all atoms   EXACTLY 0.0     (bit-identical)
      -> max |d rmsd|                                EXACTLY 0.0
      -> unchanged under OMP/MKL/OPENBLAS_NUM_THREADS = 1 vs default

### The actual mechanism, which is worse than a seed

`lam_path` solves at **lam = 0** from four generic starts and takes the argmin. Those four
converge to objectives that agree to ~1e-7 -- **numerical noise** -- while sitting on
**different torsion branches**. 1A13:

    lam=0   0.5091924865   0.5091925170   0.5091924488   0.5091924188   <- spread 1e-7
    lam=0.3 0.9470961958   0.5190590069   0.5615765828   2.3960369367   <- spread 1e+0

The winner of that noise-level argmin becomes the **warm start** for the lam=0.3 rung. So a
decision taken where the objective **cannot discriminate** determines an outcome where it
**can**, and the two differ by seven orders of magnitude. **The lam=0 branch selection is a
selection made on noise.**

### Demonstrated, not argued

Perturb the input cloud by **1e-14 relative** -- not a stress test, that is the scale at which
two code paths that build the same average actually disagree (S28-L27b measured the production
cloud against S27's rows at 5.7e-14) -- and reproject, same job, same code path:

    input moved:  max |dcoord| 1.39e-12 A   (the cloud is essentially unchanged)
    output moved: n = 126 targets

      |d chain|   p50 1.63e-03   p75 8.29e-03   p90 1.75e-02   p95 4.22e-02   p99 1.86e-01
                  mean 1.35e-02   max 5.11e-01
      counts      >1e-6: 122/126    >1e-4: 97/126    >1e-3: 71/126
                  >0.01: 28/126     >0.05:  6/126    >0.1:    4/126
      worst       2LNG 0.5113   1RSW 0.1970   2NDN 0.1519   9BFL 0.1309   8HVS 0.0789

**NOT ONE TARGET of 126 is unchanged to 1e-9.** The median target's emitted chain moves
1.6e-3 A -- nine orders of magnitude more than its input did.

**2LNG's emitted chain moves +0.511 A.** That reproduces the historical 2LNG discrepancy
(0.517 A between `s29_O_chain_rows.jsonl::prod` and `chain_rows.jsonl::DIS`) essentially
exactly, from an input perturbation of ~1e-13 A. Amplification ~1e13.

### This is a property of the instrument, not one pathological target

The lam=0 branch margin, all 126 targets:

    p1 4.56e-10   p10 5.39e-09   p25 2.68e-08   p50 3.11e-07   p75 2.79e-04   p90 1.42e-02
    below 1e-9:  3/126     below 1e-7: 52/126     below 1e-6: 73/126
    below 1e-5: 82/126     below 1e-4: 91/126

**58% of the benchmark (73/126) has its branch chosen at a margin below 1e-6**, and a quarter
of it below 2.7e-8. This is not one pathological target; it is the normal condition of the
operator. The median is 3.1e-7 -- i.e. **the typical target's branch is selected by numerical
noise**, and it is only because most branches happen to lead to similar structures that the
endpoint is stable at all.

### Consequences, and they bind on every lane

1. **Reprojection is reproducible only from bit-identical input clouds.** Both sides of any
   built-chain contrast must be projected **in the same job from the same stored clouds**, and
   the entry must say so. Two sides from different code paths carry a per-target floor with a
   0.5 A tail.
2. **The canonical endpoint is 3.2105 A**, `s29/results/s29_O_chain_rows.jsonl :: item=prod`.
   It is canonical because it is the run the endpoint was declared from **and** the only one
   whose input clouds are persisted target-by-target (`s29/results/s29_O_structs/*.npz`), so it
   is the only one that can be reprojected from its own bits. **Not changed, not retro-fitted.**
   This lane's independent reprojection reproduces it **bit-for-bit, per target, max |diff|
   0.0** -- so the endpoint IS deterministically reproducible today, from the pinned inputs.
3. **The 0.0107 A spread is verified and its composition is now known.** Five values circulate
   for "production, built chain" (S30 AUDIT_Z row 6): 3.2105 (s29 O prod), 3.2126
   (`chain_rows.jsonl::DIS`, and independently `s30_P_chain_rows.jsonl::PROD_chain`), 3.2071
   (s30 lane R / the meter), 3.2148 (s30 lane X), 3.2041 (`s30_P_chain_rows.jsonl::rec_fit`).
   Range 3.2148 - 3.2041 = **0.0107**. I recomputed 3.2105, 3.2126 and 3.2041 from their
   artefacts; 3.2071 and 3.2148 are taken from the S30 record and not recomputed here.
   **Against this project's one confirmed effect of 0.0221 A, any built-chain claim below
   0.0107 A is inside the instrument's own reprojection noise.** I will flag any that appears.
4. **A heavier consequence nobody has drawn.** The measured cloud->chain transfer coefficient of
   0.92 (S30-L25 annotation) is an **average over a map that is locally chaotic on a substantial
   minority of targets**. Any lane proposing a cloud-level improvement must know the chain
   response is not smooth there: a small cloud gain can be erased or reversed by a branch flip.

### What "pinning" actually is here, since there is no seed

Two halves, both in `s31/results/s31_D_projection_pin.json`:
(a) the **operator pin** -- every constant that determines the output (STARTS, gradient mode,
FD_EPS, lam, maxiter, penalty, tie-break rule, and an explicit `"rng": "NONE"`), digested, so a
future run can assert it is the same operator; (b) the **input pin** -- the sha256 of each
canonical cloud's float64 bytes, because the input is where the variation actually entered and
without this half the operator pin is worthless.

### Not done, deliberately

I did **not** change the tie-break, the start set, or any default in `core/project.py`. Every
number in the record stands. The follow-on question -- whether deciding the branch at lam=0.3,
where the objective separates by 1e-1, is worth anything in **accuracy** -- is a separate,
pre-registered experiment (`s31/PREREG_S31_D_branch.md`), registered before any of its arms was
computed, with H0 (no accuracy gain) as the expected outcome.

Artefacts: `s31/s31_pin_projection.py`, `s31/results/s31_D_projection_pin.json`,
`s31/results/s31_D_projection_pin_rows_shard{0,1,2,3}of4.jsonl`.
Registered in `s31/MULTIPLICITY.md` rows 1-3 (instrument characterisation, no alpha spent).


## S31-L7 -- **DEFECTS D-A AND D-C FIXED IN PLACE.** A WITHDRAWN POSITIVE HAD BEEN ASSERTED IN SHIPPED CODE FOR SIX SPRINTS; AND THE LAUNCHER HAD BEEN CONTRADICTING THE GOVERNOR SINCE S29 BECAUSE ONLY ONE HALF OF A PAIRED THRESHOLD WAS EVER MOVED (2026-09-21 00:10, D)

### D-A: `core/pipeline.py`, the CVaR docstring

`quantum_stage`'s docstring asserted, under the heading **"WHY CVaR, MEASURED"**:

> *"So the CVaR tail is worth +0.113 A **by preventing the collapse**, and that is the
> component's measured role."*

**S25-L5 (`s25/LEDGER.md:233`) withdrew that number.** It was read off MARGINAL MEANS where a
PAIRED statistic was required. Paired, `s25/results/q_alpha.json`, n = 126:

    vqe_a0.1_T0.1 - vqe_a1.0_T0.1   -0.1126  SE 0.0792  MDE 0.2220  0.51x MDE
                                    fold CI spans zero, 55W/44L, median EXACTLY 0.0000  -> NULL
    VQE_LFO - argmin  (the replacement figure, s25/LEDGER.md:280)
                                    -0.1405  SE 0.0732  MDE 0.2051  0.68x MDE
                                    66W/48L, 5/5 folds same sign  -> NOT MEASURED

Sign convention stated because it has to be: lower RMSD is better, so both contrasts are
**nominally in the CVaR arm's favour and neither clears its MDE**. The defensible statement is
that alpha < 1 prevents the entropy collapse -- a fact about the **state**, visible in 0.076 vs
6.36 bits -- and that **no Angstrom effect of that mechanism has been measured.**

Corrected in place, with the old wording **quoted verbatim in the new text** so the correction
is visible to anyone who read the old version rather than silently erased, and with a direct
instruction not to quote a CVaR contribution from that docstring.

**This is the sixth instance in this project of prose asserting a state that does not hold, and
the first in shipped code** -- which is the part that matters: `s26/pr_changes.py:79` recorded
in S26 that this docstring carried a withdrawn claim, and the docstring was still wrong four
sprints later. Recording a defect in a table is not fixing it. A reader who greps the source
rather than the ledger got a withdrawn positive presented as an established role.

### D-C: `s26/jobrun.py` v3 -- the launch gate is now DERIVED from the governor

**The defect.** `jobrun.py` last changed in **S26** (v2.3). `governor.py` changed **four times
in S29** -- v2.5, v2.6, v2.6b, v2.6c -- and every one of those changes moved in the same
direction, away from treating CPU as a safety signal, ending at
`CPU_CEILING = 101.0  # CPU-triggered suspension is DISABLED ... A saturated CPU has no failure
mode; RAM does (OOM)`. **The launcher never followed.** It still refused to launch while the
governor's 15 s CPU mean exceeded **85%** -- below the **94-95%** band the charter instructs.
So in exactly the regime the charter asks for, `jobrun` blocked every launch while the governor
was content, and lanes responded by launching **detached**, which registers nothing in
`s26/jobs/` and leaves the governor unable to suspend, resume, kill or even see the process.
**That is strictly worse than either policy alone.**

Caught live while writing the fix: `s26/governor_state.json` at 23:52:33 read `cpu_smooth 96.8,
n_jobs 0`. The 96.8% was four of **this lane's own** shards, launched detached for this exact
reason. I am the instance.

**The fix, and why this one.** RAM is the binding constraint on this box: 16.75 GB total, ~10-11
GB of it the user's own baseline load, and an OOM is the only way a job here takes the machine
down. A saturated CPU makes everything slower and nothing unsafe. So the launch gate is now
**imported from `governor.py`** rather than duplicated:

    CEILING   = governor.CEILING - 1.0   = 93.0   RAM, UNCHANGED IN VALUE, one point of
                                                  hysteresis below the suspend ceiling
    CPU_START = governor.CPU_CEILING     = 101.0  i.e. CPU does not gate launches
    _cap()    = min(launch_cap.json, governor.MAX_JOBS)

**The structural point is the import, not the constant.** If a future sprint re-enables
CPU-triggered suspension by lowering `governor.CPU_CEILING`, the launcher follows automatically
and nobody has to remember to grep for the partner. This is the **third** instance of the
standing `paired-thresholds-move-together` failure, and the first fix that makes the next one
impossible rather than merely unlikely.

**Nothing safety-relevant was relaxed.** The per-job headroom test (`avail < est_ram + 0.5 GB`),
the refusal to start unsupervised above 0.5 GB with no live governor, the AMBER cap of 2, and
the concurrency cap all stand -- and the cap is now **clamped** to the governor's `MAX_JOBS`, so
an over-large `launch_cap.json` cannot register more jobs than the governor will manage. The
fallback when `launch_cap.json` is unreadable stays the conservative 4, deliberately not the
governor's 8: an unreadable cap file is a fault, and a fault should launch fewer jobs, not more.
`--cpu-gate` restores a CPU launch gate for any caller that wants one.

**Verified live, not asserted.** A job launched successfully at `cpu_smooth 100.0%` -- which the
old gate would have blocked indefinitely -- and the four shards of `s31_D_branch` are registered
in `s26/jobs/` and visible to the governor as I write this.

Two smaller things fixed in the same file while I was in it:
- the registration and completion records were written through a **shared** `<name>.tmp`;
  `os.replace` is atomic but a shared temp path is not, so two writers can publish an
  interleaved file atomically. Now `<name>.<pid>.tmp`.
- the wait message said only `box at ram X% cpu Y% jobs n/m`; it now names **which** gate is
  actually blocking, because "waiting" with no reason is how a lane decides to go detached.

**Correction to a standing memory:** the note *"s26/jobrun.py does not dedupe by --name"* is no
longer true of the current file -- `main()` reads `s26/jobs/<name>.json` and refuses with exit 3
if that pid is alive. The S29 incident of four copies of one job is fully explained by
**detached launches bypassing jobrun entirely**, which is the defect fixed above. The operative
half of that lesson stands and is now load-bearing: **check the pid, not the name** -- and a
detached process has no registration to check at all.

### The rule this puts on every lane

**Stop launching detached.** `nohup python s31/foo.py &` is invisible to the governor: it cannot
be suspended when RAM climbs, cannot be killed before an OOM, and does not appear in
`governor_state.json`. Launch through `s26/jobrun.py`, which now works in the band the charter
asks for. If `jobrun` still blocks you, it will now tell you which gate and why -- send me that
line rather than going around it.

Artefacts: `core/pipeline.py` (docstring), `s26/jobrun.py` (v3), `s31/MULTIPLICITY.md`.


## S31-L8 -- **THE BRANCH SET CONTAINS 0.0938 A OF REAL ACCURACY AND THE OBJECTIVE CANNOT FIND IT.** DECIDING THE BRANCH WHERE THE OBJECTIVE DISCRIMINATES (1e-4 INSTEAD OF 1e-7) FIXES THE CONDITIONING BY **1082x** AND IS WORTH **-0.0055 A AT 0.44x MDE -- NOT A RESULT.** H0 AS PRE-REGISTERED (2026-09-21 00:10, D)

Pre-registered in `s31/PREREG_S31_D_branch.md` **before any arm was computed**, including the
decision rule, the expected outcome (H0), and the structural caveat below. Proposed by the
coordinator off S31-L<D-B>'s mechanism.

### The arms, all projected in ONE job from the SAME stored clouds

`s29/results/s29_O_structs/<pdb>.npz :: prod`, n = 126, built chain, `grad="exact"`, lam 0.3.

    PROD        lam_path(multi=True) replicated exactly -- the incumbent
    B4          continue EACH of the four lam=0 branches to lam=0.3, union the four fresh
                lam=0.3 starts, argmin of the LAM=0.3 OBJECTIVE over all eight.
                Native-free, deterministic, DEPLOYABLE.
    ORACLE_B4   the same eight candidates, argmin of CA-RMSD TO THE NATIVE.
                **ORACLE / NOT DEPLOYABLE.**

**Identity gate passed before anything was read:** PROD reproduces
`s29_O_chain_rows.jsonl :: prod` **bit-for-bit on all 126 targets, max |diff| exactly 0.0**. So
the arms are comparable to the canonical endpoint and both sides of every contrast below come
from the same code path, as this lane's own D-B rule requires.

### The ceiling first, as pre-registered

    ORACLE_B4 - PROD   -0.0938 A   SE 0.0114   MDE 0.0321   2.92x MDE
                       fold CI [-0.1097, -0.0786]   114W / 0L / 12 tied   median -0.0279
                       best single target -0.573
                       **ORACLE / NOT DEPLOYABLE -- this is a ceiling, not a method.**

**The eight candidates the projection already computes contain 0.0938 A of accuracy, and it
clears its MDE at 2.92x with 114 wins and zero losses.** That is four times the project's one
confirmed effect (0.0221 A) and it is sitting inside an operator that is run on every target of
every arm of every experiment. It is ORACLE and it is NOT DEPLOYABLE.

### The native-free arm: the primary

    B4 - PROD          -0.0055 A   SE 0.0044   MDE 0.0124   **0.44x MDE**
                       34W / 21L / **71 tied**   median EXACTLY +0.0000
                       fold CI [-0.0071, -0.0031]   iid CI [-0.0144, +0.0031]   5/5 folds
                       -> **NOT A RESULT** by this project's fixed rule.

**I am stating this as a null and I want to be explicit about why, because the temptation here
is real.** The fold CI excludes zero and all five folds agree in sign. That is not enough. The
rule is MDE = 2.8016 x SE per comparison, and **below 0.7x is not a result** -- this is 0.44x.
The iid CI spans zero. Reaching for the fold CI because it is the one that excludes zero,
having seen both, is the exact move this project has a standing rule against. **B4 is a null on
accuracy.**

And a second disqualification, from my own instrument work an hour earlier: **-0.0055 A is
half the instrument's own 0.0107 A reprojection spread.** By the rule I set for every other
lane this sprint, a built-chain claim of that size is inside the noise of the thing measuring
it. The rule applies to me.

### Concentration: suggested by the median, NOT established by the null

    mean -0.0055   median EXACTLY 0.0000   71 of 126 tied
    drop-top-10 mean +0.0031 (the SIGN FLIPS)
    uniform-effect null: p10 -0.00038, p50 +0.0027, p90 +0.0069
    -> observed drop-top-10 sits at the **55.9th percentile of the null.  flag = FALSE.**

The median-vs-mean gap is exactly the free early warning the standing rule describes, and the
sign flip on drop-top-10 looks alarming -- **and the matched uniform-effect null says it is
ordinary.** A raw drop-top threshold is not a valid test; the null is. Reported as SUGGESTED
and NOT ESTABLISHED, which is what the rule requires and what it would have required if the
answer had gone the other way.

### Strata, reported whichever way they fell

    FAIL18     +0.0049   0.08x MDE   4W/3L    (NOT MEASURED, and a null at n=18)
    other 108  -0.0072   0.69x MDE   30W/18L  (below 0.7x: NOT A RESULT)

No switch, no stratum rescue.

### What DID move, and it is the thing the exercise was for

    decision margin, median:   PROD (at lam=0)  3.11e-07  ->  B4 (at lam=0.3)  3.36e-04
    targets with margin < 1e-6:            73/126        ->             21/126
    improvement in the median margin:                 **1082x**

**The conditioning fix is real and it is large.** B4 decides the branch where the objective
separates by ~1e-4 instead of ~1e-7, which collapses most of the 1e13 amplification that makes
the instrument irreproducible. That is a property of the operator and it does **not** depend on
the accuracy null above.

### The structural caveat, registered before the result and repeated here

**B4's candidate set is a strict superset of PROD's and both minimise the same objective, so
B4's objective is <= PROD's on every target BY CONSTRUCTION** -- verified, 126/126. "B4 reaches
a lower objective" is therefore not evidence of anything and is not offered as any. B4 is
strictly a harder search on the same objective, and it bought **-0.0055 A at 0.44x MDE** while
changing the emitted branch on **55 of 126 targets**, with a per-target range from -0.348 to
+0.321. It moves a great deal and nets nothing.

### The reading, and it is the project's own central finding again

The branch set contains **0.0938 A** (ORACLE, NOT DEPLOYABLE, 2.92x MDE, 114W/0L). The
objective recovers **0.0055 A of it at 0.44x MDE**, i.e. **within noise of nothing** -- about
6% of the ceiling, and not distinguishable from zero. So: **search is not the barrier here;
discrimination is.** A 1082x better-conditioned argmin over a set demonstrably containing 0.094
A finds essentially none of it, because it is an argmin of an objective that does not rank the
native. This is the same wall as S15, S29's certified optimum and S30, arriving from a new
direction -- the numerical conditioning of stage 3b -- and it is worth recording precisely
because the direction was new and the wall was in the same place.

### Recommendation, and it is not mine to take

1. **Do not adopt B4 for accuracy.** It is a null and it is inside the instrument spread.
2. **B4 is defensible purely as a CONDITIONING fix** (1082x, 73/126 -> 21/126 below 1e-6), and
   if the coordinator wants the endpoint reproducible against re-implementation rather than
   only against re-running the same bits, that is the argument for it -- **stated as a
   conditioning change with a null accuracy effect, never as an improvement.** It would move
   the canonical 3.2105 to 3.2050 and that is a decision above this lane.
3. **The 0.0938 A ORACLE ceiling is the interesting object and it is native-free-adjacent**:
   the candidates exist, are already computed, and cost nothing. Whether any native-free signal
   orders them better than the objective does is a real question for a future sprint -- and the
   honest prior, from S12/S29/S30, is that nothing will.

Artefacts: `s31/s31_D_branch.py`, `s31/results/s31_D_branch.json`,
`s31/results/s31_D_branch_rows*.jsonl`, prereg `s31/PREREG_S31_D_branch.md`.
Registered in `s31/MULTIPLICITY.md` row 7 (k = 2 primary, one basis).

---

## S31-L9 -- **92.9% OF tuning126 IS NMR-DETERMINED, SO THE NEW-OBSERVABLE QUESTION IS CLOSED BY PROVENANCE, NOT BY GEOMETRY** -- AND THE ALARMING FOLLOW-UP IS A NULL: THE DEPOSITED ENSEMBLE IS 1.08 A WIDE BUT **UNCORRELATED WITH FAIL18** (0.20x MDE) (2026-09-21 00:13, L)

Full working: `s31/LIT_L.md` §L1.5, §L2, §L2.1. Scripts: `scratchpad/lit_L_expmethod.py`,
`lit_L_ensemble_spread.py`, `lit_L_tail_vs_spread.py`.

### OURS -- what actually determined our reference coordinates (RCSB GraphQL, all 126)

```
  115  ( 91.3%)  SOLUTION NMR          4  (  3.2%)  X-RAY DIFFRACTION
    5  (  4.0%)  ELECTRON CRYSTALLOGRAPHY    2  (  1.6%)  SOLID-STATE NMR
  -> NMR-determined 117/126 = 92.9%
```

### THE L2 VERDICT: THE BINDING CONSTRAINT IS INFORMATION PROVENANCE, NOT G1

1. **Any observable computed at inference from `(sequence, pool)` adds NO information**,
   whatever equivalence class it occupies -- data-processing inequality. It can only be a
   better estimator, and that ceiling is already measured (`in-band-ordering-is-per-target`:
   0.600 across targets vs the 0.638 needed for 2.0 A). **G1 is not the binding constraint.**
2. **Any NMR observable of these targets IS the data that determined the reference.** For
   117/126 the deposited coordinates are a fit to deposited NOEs, J-couplings and torsion
   restraints. Feeding them back is **ORACLE through a different door**.
3. **This re-prices the chemical-shift direction retroactively.** TALOS+/TALOS-N dihedral
   restraints derived from shifts are standard practice in NMR peptide structure
   determination, so on NMR-determined targets "ORACLE-perfect torsions" was close to a
   **tautology**. *Stated as standard practice, NOT verified per-target here.* The direction
   is already closed; the **reason** should be recorded correctly.
4. **The only genuine escape is an observable measured on the molecule and NOT used in its
   structure determination.** Best physical candidate: **VCD / ROA** -- genuinely chiral, so
   outside G1 by construction, and they work **in our 9-16 band** (ROA run on ~6-residue
   peptides; Keiderling, *Chem. Rev.* 2020, 120(7):3381-3419). **But there is no repository of measured VCD/ROA
   spectra keyed to PDB entries** -- I searched and found only method papers. Measured-spectrum
   count for these 126 targets: **zero**, against 54/126 for shifts which (2) disqualifies.

> **L2 is CLOSED for this benchmark, by DATA AVAILABILITY rather than physics.** The physics
> leaves chiral and many-body channels open; the benchmark supplies no measurement to put in
> them. Same wall as `no-fresh-benchmark-exists`, arriving from the observable side.

### THE ALARMING FOLLOW-UP, AND THE NULL THAT KILLS IT

The manifest reference is **"deposited coordinates, MODEL 1"** -- an arbitrary member of an
NMR ensemble. ORACLE DIAGNOSTIC of the benchmark, never an inference-time signal; 111/126
ensembles resolved.

```
mean pairwise CA-RMSD BETWEEN DEPOSITED MODELS : mean 1.0823  median 0.9929  max 4.2648
model 1 -> ensemble medoid                     : mean 0.6965  median 0.4565  max 4.2943
spread >1.0 A : 55/111     >2.0 A : 15/111     >3.0 A : 2/111
widest: 3BTB 4.265, 6CEJ 4.200, 6GIJ 2.958, 6EY3 2.674, 2MIG 2.629
```

**The reference really is uncertain. It does NOT explain the tail.** Production arm, n=111:

```
corr(production RMSD, ensemble spread) = +0.1118  95% CI [-0.0762,+0.2921]  SPANS ZERO
corr(production RMSD, model1->medoid)  = +0.1341  95% CI [-0.0536,+0.3127]  SPANS ZERO

WORST 18 by production RMSD : mean RMSD 6.0291   mean ensemble spread 0.9672
the other 93                : mean RMSD 2.6842   mean ensemble spread 1.1046
difference tail-minus-rest  : -0.1374  SE 0.2421  MDE 0.6783  -> 0.20x MDE = NULL
```

The tail's targets have if anything **slightly narrower** deposited ensembles. **FAIL18 is
real failure against a reference no worse determined than any other target's.** Same null in
all four other arms (0.08x-0.90x MDE, every CI spanning zero).

**Reported prominently BECAUSE it is a null.** Stopping at "the reference is uncertain by
1.08 A" would have been quotable, alarming and wrong, and would have redirected the sprint's
tail work. `control-at-the-decisive-step`: the cheap null sits exactly where the wrong answer
would first have become quotable.

### WHAT SURVIVES -- an interpretive caveat on the ABSOLUTE number only

The endpoint carries a **uniform** ~0.70 A reference term: a perfect predictor aiming at the
ensemble medoid still scores ~0.70 A against model 1. In quadrature that is
`sqrt(3.21^2 - 0.70^2) = 3.13` vs 3.21 -- **~0.08 A now**, and far more material near 1 A.
Because it is uniform it **cancels in every arm-to-arm delta**, which is the project's actual
currency.

> **NOT a reason to re-score against the medoid.** `benchmark-and-folds-must-be-pinned`: the
> reference is part of the sealed instrument. Report the uncertainty; change nothing.

### ADAPT-VQE / qubit-ADAPT -- ranked last, with reasons

qubit-ADAPT (arXiv:1911.10205), ADAPT-QAOA (arXiv:2005.10258). (a) ADAPT is an **expressivity**
fix and S31-L4 shows our expressivity gap has a **free** alternative -- ADAPT would be an
efficient way to approximate a quantity with a closed form. (b) **ADAPT's selection rule is
undefined for CVaR**: the criterion `|<psi|[H,A]|psi>|` presumes the cost is `<H>`, a *linear*
functional of the state, and CVaR is not the expectation of any observable. Not worth a lane.

### A SMALL INTEGRITY NOTE FOR LANE D

`results/summary/results.csv` gives the production mean as **3.2126** (n=126) against the
charter's **3.2105** -- a 0.0021 discrepancy. Flagged, not reconciled: the charter already
lists **the unpinned projection seed** as an open defect gating every sub-0.01 A claim, and
this is consistent with exactly that. Lane D owns it.
