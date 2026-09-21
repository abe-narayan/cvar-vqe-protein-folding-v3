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
`s31/lit_L/lit_L_dequant_check.py`, `s31/lit_L/lit_L_gap_budget.py`.

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

Full working: `s31/LIT_L.md` §L1.5, §L2, §L2.1. Scripts: `s31/lit_L/lit_L_expmethod.py`,
`s31/lit_L/lit_L_ensemble_spread.py`, `s31/lit_L/lit_L_tail_vs_spread.py`.

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

## S31-L10 -- **THE FREE ENERGY IS ACHIRAL FOR THE SAME REASON THE ENERGY IS, SO G1 CLOSES BOTH HALVES OF `F = E - TS` -- AND THE TORSION CHANNEL'S ESCAPE FROM G1 IS REAL, MEASURABLE, AND WORKS ONLY ON THE PROBLEM THE PIPELINE DOES NOT HAVE.** THE 32-COST SWEEP REPRODUCES EXACTLY AND ITS TOP ROWS MEASURE TRIAGE, NOT NATIVENESS: **`pref(circ_best vs PROD)` IS AT OR BELOW A COIN FLIP FOR EVERY ONE OF THEM** (2026-09-21 00:15, B)

Pre-registered in `s31/PREREG_S31_B.md`, committed at **e2fc6257 before the first number**.
Artefacts: `s31/results/s31_B1_achirality.json`, `s31_B2_sweep_parity.json`,
`s31_B2_inpool.json` + `s31_B2_inpool_rows.jsonl`, `s31_B3_graph.json`, and the three `.log`s.
Code: `s31/s31_B1_achirality.py`, `s31_B2_sweep_parity.py`, `s31_B2_inpool.py`, `s31_B3_graph.py`.
**Basis:** B2's correlations are **per-candidate CA point-cloud RMSD to native**, which is neither
the 3.2105 Å built-chain endpoint nor the 3.0483 Å cloud endpoint. The sweep contrasts are on the
**built-chain** rungs. **ORACLE / NOT DEPLOYABLE**: every correlation here is labelled by `rr`, and
the coherence quantity needs the native in both arguments. Lane B emitted **146 comparisons**.

---

### 1. B1 — the free-energy stage, CLOSED BY DERIVATION, no thermodynamics computed

**What `S` is a function of.** `S(x) = Φ_{β,B}[U_seq](x)`. Its arguments are the candidate `x` (a
pool member), the sequence (only through the atom typing that instantiates `U`), and `β` and the
basin map `B`, which are universal. **There is no fourth argument**, so by lane P's repaired source
enumeration `S` is an *operator* on (sequence, pool), not a third source — exactly like `E`.

**Lemma B1, verified not assumed.** The shipped potential
`app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")` (`core/amber.py:927`) is
reflection-invariant, `U(Rx) = U(x)`.

**Corollary B1.** `F_β`, `E`, `S` and every temperature derivative of `F` are rotation-,
translation- **and reflection**-invariant single-structure observables, so **by G1 each is a
function of the candidate's distance map** and lies in the class S30 closed by theorem and measured
empty on 43 channels. The enthalpic half was already closed by measurement; the entropic half is
closed by the same theorem that closed contact topology and Rg, because **`S` is achiral for exactly
the reason `E` is.**

> ### Corollary B1′ — the part nobody had connected, and the more valuable half
> An **ANM or GNM Hessian is built from pairwise distances**, so its spectrum, its log-determinant
> — which *is* the harmonic configurational entropy, `S_harm = const − ½ ln det′ H` — and **every
> spectral-graph observable derived from it are distance-map functions.** This closes the whole
> *elastic/network — normal-mode — energy-landscape-curvature* family of charter §14 by the same
> theorem, without computing one of them. Basin populations, ensemble reweighting and
> temperature-dependent ranking are monotone transforms of single-structure free energies and close
> with it; local/per-residue configurational entropy closes **twice**, by B1′ and by S30-R's
> locality result.

**F2, registered falsifier — did not fire.** ff14SB's `PeriodicTorsionForce`, **1340 torsions**,
maximum distance of any phase to `{0, π}` = **0.000e+00 exactly**; no `CMAPTorsionForce`. The
algebra `cos(−nφ − γ) = cos(nφ − γ) ⟺ γ ∈ {0, π}` therefore holds term by term.

**F1, registered falsifier — FIRED AGAINST ITS THRESHOLD, AND THE MATCHED CONTROL INVERTS THE
READING.** `max |U(x) − U(Rx)| / (|U(x)|+1)` over 10 real target structures (208–390 atoms) is
**3.673e-06** against my pre-registered **1e-6** bar. A **pure proper rotation**, equally
analytically exact for this potential, gives **5.377e-06 — larger**; reflection/rotation = **0.68×**.
Per force group, reflect vs rotate: `PeriodicTorsionForce` **2.0e-16** vs 3.6e-16, `HarmonicBond`
1.9e-15, `HarmonicAngle` 1.3e-15 — *the only terms that could carry chirality are exact to machine
epsilon*. The 1e-6-scale residual is confined to `NonbondedForce` (5.5e-07 vs 3.2e-06) and
`CustomGBForce` (6.8e-07 vs 7.9e-07), which are pairwise-distance functions and analytically
invariant under **any** isometry; their residual is summation-order noise.

> **A registered threshold set at an absolute value, on a quantity only a matched control can read,
> is a mis-set threshold, not a falsification.** The original wording stands; F1 is never quoted
> without the rotation control in the same sentence. *(14 of 24 PDBs failed to build — non-standard
> residues `CGU`/`DAL`, waters, missing atoms — unrelated to the test.)*

**Left OPEN, on price and not on theory.** A **multi-structure** free energy — a barrier
`F‡(x_a → x_b)` — depends on the *path* and so is not a function of `D(x_a)` and `D(x_b)`. G1 does
not bind it. 10³–10⁵ trajectories per target is out of reach here; it is neither closed nor run.

**Registered directional prediction, recorded before measurement (§B1.5), for whoever runs it.** A
free energy **restrained toward the distogram prior is coherent with the prior's error by
construction**, hence with the pool common mode, hence **harmful rather than neutral**
(+0.0554 coherent against −0.2466 i.i.d. at matched R²).

---

### 2. The 32-cost meter sweep — REPRODUCED EXACTLY, and what its top rows actually measure

Lane B rebuilt the published statistic independently and, once the meter's own tie rule
(`s29_D_cost_audit.py:454`, **ties count 0.5**) is used, reproduces it to four decimals:
`LEG_steric` **+0.2515**, `CAGEO` **+0.2341**, `LEG_torsion` **+0.2212**, `DIS` **+0.0357**,
`DSSPHB` **−0.1835**. **The rows are correctly computed. Three things about them are new.**

**(a) It is NOT an artefact of the `RAND_SIGNED` construction.** The same contrast against
`GAUSS_MATCHED_0` — a Gaussian perturbation at matched distance, not a signed affine combination of
the pool — is **the same size**: `LEG_torsion` +0.2222 vs +0.2212, `LEG_steric` +0.2460 vs +0.2515,
`CAGEO` +0.2778 vs +0.2341. **The rows stand as measurements of what they measure.**

**(b) What they measure is that the CONTROL is worse, not that the near-native is better.**
Splitting the contrast into its two terms:

```
cost              pref(circ_best vs PROD)   pref(RAND_SIGNED vs PROD)   contrast
LEG_steric                  0.5079                    0.2564            +0.2515
CAGEO                       0.4206                    0.1865            +0.2341
LEG_torsion                 0.4048                    0.1835            +0.2212
RAMA                        0.4127                    0.2629            +0.1498
```

**Not one of the top rows reaches 0.5 on `pref(circ_best vs PROD)`** — `LEG_steric` is a coin flip
at 0.508 and `LEG_torsion` is **below** one at 0.405, i.e. it prefers production to the ORACLE
near-native structure on 60% of targets. A cost that recognised nativeness would score high there.
**The contrast is coarse triage against a displaced structure.**

**(c) The torsion-family contrasts are carried by the CHIRAL part; `LEG_steric`'s is not.**
Decomposing each cost into `f_even = ½[f(x)+f(Rx)]` and `f_odd = ½[f(x)−f(Rx)]` under the exact
point reflection and recomputing the *same* contrast:

```
cost              contrast_tot   contrast_even   contrast_odd    odd/tot
LEG_torsion         +0.2212        +0.1220         +0.2698  2.78x  +1.22
RAMA                +0.1498        +0.1002         +0.2242  2.25x  +1.50
LEG                 +0.0546        +0.0317         +0.2133  1.87x  +3.91
TORS_CONS_POOL      +0.0159        -0.0714         +0.1429  1.40x  +9.00
LEG_steric          +0.2515        +0.1523         +0.1225         +0.49
CAGEO               +0.2341        +0.1845         +0.0863         +0.37
DIS / CONTACT       +0.0357        +0.0357         +0.0000          0.00   <- self-check
```

**`LEG_steric` splits roughly evenly, so my own opening hypothesis about it was half wrong and is
corrected here**: its *variance* is dominated by chirality (odd share 13.0) but its *contrast* is
not. **`CAGEO` is not achiral either** (mirror self-check 6.57) though its contrast is mostly even.
The two achiral CA costs come out at **exactly 0.0000** odd, which is the self-check that the mirror
is being applied correctly.

**The construction-level mechanism, not a correlation.** `core/geometry.build_backbone_batch` places
`CB = −0.58273431·cross(b, d) + …`. `cross` is a **pseudovector**, so every CB-dependent Legacy term
inherits chirality. Measured odd variance share on 40 targets × 300 real pool members: steric
**13.02**, torsion **0.301**, electrostatic 0.230, contact 0.105, solvation 0.100, aromatic 0.089 —
and **exactly 0.000** for `hbond_local`, `hbond_longrange`, `coop_helix`, `coop_sheet`,
`compactness`, which are the achiral terms.

---

### 3. B2 — the torsion channel, CLOSED. The escape from G1 is real and works on the wrong problem.

**The decomposition theorem.** For an ideal-geometry backbone the point reflection acts on torsions
**exactly** as `(φ,ψ) → (−φ,−ψ)`. So `T_even` is a reflection-invariant single-structure observable
and **by G1 a function of the distance map**. **All of a torsion channel's escape from G1 lives in
`T_odd`.** Crossed with separability (S30-R: a sum of per-residue terms cannot see a lever arm), the
only cell no existing theorem closes is **chiral AND non-separable**. It was built (`XTWIST`) and
measured.

**The result, in one table.** `LEG_torsion`'s odd variance share over the pool is **0.3001** (median
0.2844) — so the channel is genuinely chiral, and the registered `G1-check` (fires below 10%) does
not fire. But the skill sits in different halves on the two problems:

```
                                 EVEN half (G1-CLOSED)      ODD half (G1-ESCAPING)
coarse triage (sweep contrast)   +0.1220   1.35x MDE        +0.2698   2.78x MDE
in-pool selection, 500 band      +0.2428   2.51x MDE        +0.0183   0.28x MDE
in-pool selection, top-75 band   +0.0288   0.40x MDE        +0.0423   0.66x MDE
```

> **The chiral escape from G1 is real, measurable and not degenerate — and it only works on a
> problem the pipeline does not have.** On the problem it does have, the surviving skill is in the
> G1-closed half, and even that dies in band.

**G2 fires exactly as registered.** In-pool ρ on the shipped **top-75** band: `LEG_torsion` +0.0444
at **0.65× MDE**, `RAMA` +0.0368 at 0.51×, `LEG_tors_odd` +0.0423 at 0.66×. **Below 0.7× is not a
result.** For scale, the *shipped* cost `DIS` scores +0.0652 at **0.83× MDE — NOT MEASURED** in its
own band. **G3 does not fire**: no channel's partial ρ on `DIS` clears 1.0× MDE (best +0.0401,
0.61×). Partialling on Rg changes nothing.

**G5 fires, and this is the sharpest single number.** A **constant α-helix** at (−57°, −47°) — a
plausible zero-*target*-information reference, never uniform-on-the-torus (contract rule 9) —
**beats every torsion channel on both bands**: on the 500 band `HELIX_CONST` ρ = **+0.3413 (2.89×
MDE, 5/5 folds)** against `LEG_torsion`'s +0.1676, a paired gap of **−0.1722 at 3.17× MDE with the
fold CI [−0.2081, −0.1453]**; in band the gap is −0.0901 at 1.81× MDE, 5/5 folds. `HELIX_CONST` also
beats Rg (+0.2794). **The Ramachandran channel's in-pool skill is a constant-prior effect — "this
benchmark's peptides are helical" — and a single constant does it better than the fold-conditioned
channel.**

**M6, the open cell, is empty.** `XTWIST` = mean `sin(crossing dihedral)` over spatially contacting
segment pairs at sequence separation ≥ s — chiral by construction (`χ → −χ` under reflection) and
non-separable by construction. In-pool ρ at the 500 band: `XTWIST_4` **−0.0104 (−0.23×)**,
`XTWIST_8` −0.0853 (−1.02×, i.e. *backwards*); in band +0.0320 (0.36×) and +0.0091 (0.09×). **Its
own achiral twin `XTWABS` beats it** (+0.1124 and +0.1486 at the 500 band) — the matched control in
the operator's own space. The α-helix control beats `XTWIST_4` by −0.3429 at 4.06× MDE, 5/5 folds.

**M2 — why triage and selection come apart.** `LEG_torsion` of the `RAND_SIGNED` rung sits at
**z = +2.363** in the pool's *own* torsion distribution (median +2.069), and **70.6% of targets have
the control above the pool's 95th percentile**. `PROD` sits at +0.844 and the ORACLE `circ_best` at
**+1.129 — worse than production**. The channel's dynamic range is spent outside the pool; inside it
there is nothing left to spend.

**M0 — chiral dynamic range does exist.** Mirroring a candidate changes its RMSD by **0.5168 Å**
against a pool spread of 1.3206 Å, a ratio of **0.402**. *Reported as a diagnostic and explicitly
not as a ceiling*: `D(x)` determines `x` up to reflection and the real pool holds no mirror pairs,
so on **this** pool an achiral observable is not information-limited. **G1 forbids a new channel,
not a good function of the old one** — and that is what the charter asked about.

**My own registered prior failed and is recorded as failed.** I gave **4 : 1 that `LEG_torsion` is
predominantly odd**. It is **70% even**. The three other priors (no in-pool skill at 1.0× MDE;
nothing passes the coherence bar; B2 ends closed) all held.

---

### 4. M5 — WHY NO RANKER WILL EVER PASS S30'S COHERENCE ADMISSION TEST

**The derivation.** Let the readout be an affine combination of pool members, `C = Σ a_m x_m` with
`Σ a_m = 1` (argmin is `a = δ`; the uniform average is `a = 1/75`). To first order in pair-distance
space, `e_p = Σ_m a_m (d_{m,p} − d_nat,p) = μ_p + Σ_m a_m η_{m,p}`. **`Σ a_m = 1` passes the common
mode `μ_p` through with coefficient exactly one, whatever the weights are.** So
`coh = corr(e, μ)` depends on the **concentration of `a`**, not on the ranker that produced it.

**Measured, n = 126, ORACLE / NOT DEPLOYABLE (both arguments need the native):**

```
uniform mean in pair space       coh = 1.0000   sd 0.0000     <- EXACTLY, as derived
coordinate average (production)        0.9780   sd 0.0322     <- leaves the hull, barely
argmin by HELIX_CONST                  0.8684
argmin by RAMA                         0.8552
argmin by LEG_torsion                  0.8467
argmin by DIS  (the shipped cost)      0.8288
argmin by a RANDOM pool member         0.8244   <- the shipped cost and a coin are the same
argmin by LEG_tors_odd                 0.7861
ORACLE best member                     0.6708   <- the ONLY arm under the 0.6931 bar
```

> **The shipped cost and a random pick differ by 0.0044 in coherence. The entire 43-channel ranking
> search was searching a dimension along which the admission test does not vary** — and the ORACLE
> *ceiling* of all in-pool ranking is 0.6708 against a bar of 0.6931, so even perfect selection
> barely clears it. **The admission test cannot be passed by choosing better inside the pool; it can
> only be passed by leaving the pool's affine hull.** The project's one confirmed positive — the
> AMBER relax at k = 30, −0.0221 Å — is exactly such an operator.

---

### 5. The coordinator's multi-structure proxy — pre-checked before anything was built

Graph observables on the pairwise Kabsch-RMSD matrix `P` of the top-75 are functions of the **set**
of distance maps, so they escape G1 on the same argument that leaves barriers open, at the cost of
one matrix operation. Both failure modes were checked first (n = 126, kNN = 8).

**C1 — collapse to the direct distance: DOES NOT FIRE.** Within-target Spearman against the direct
`P` over all 75·74/2 pairs: **geodesic 0.867**, **commute time 0.792** (p10 0.787 / 0.662). These
are **not** relabellings of `P`; they carry real independent rank variance.

**C2 — reduction to consensus: FIRES, decisively.** In-band ρ with the ORACLE label, then the same
partialled on the **medoid criterion** `mean_j P_ij`:

```
node statistic     rho      xMDE          rho | medoid      xMDE
medoid_crit      +0.2237    2.27   5/5     +0.2094 (|DIS)   2.22   5/5
geo_cent         +0.2099    2.10   5/5     +0.0377          0.53   4/5
commute_cent     +0.1697    1.88   5/5     +0.0121          0.19   3/5
degree           +0.1519    1.87   5/5     -0.0225         -0.49   4/5
fiedler_abs      +0.0738    0.88   5/5     +0.0458          0.54   3/5
```

**Every node-level graph observable's in-band skill is fully absorbed by consensus** — which the
project already owns and already priced at −0.172 Å. Note the consensus criterion itself is
**+0.2094 at 2.22× MDE, 5/5 folds, partialled on the shipped cost**, i.e. it is the one thing in
this entry that beats `DIS` in band (`DIS` is 0.83×, NOT MEASURED).

> **The escape is real at the PAIR level and empty at the NODE level.** Anyone spending on it needs
> a **pair-level terminal operator**, not a candidate ranking. That is the handoff.

---

### 6. What lane B closed, and what it did not

| direction | status | closed by |
|---|---|---|
| single-candidate free energy, entropy, `E`/`S` decomposition, basin populations, ensemble reweighting, temperature-dependent ranking | **CLOSED** | Corollary B1 (G1 + Lemma B1, verified) |
| ANM/GNM spectra, normal-mode/harmonic entropy, landscape curvature, spectral-graph-on-one-structure | **CLOSED** | Corollary B1′ |
| local/per-residue configurational entropy | **CLOSED twice** | B1′ and S30-R locality |
| the backbone-torsion channel for **in-pool selection** | **CLOSED** | G2 fires (0.65× MDE), G3 does not fire, **G5 fires at 1.81–3.85× MDE** — a constant α-helix beats it |
| the chiral **and** non-separable cell (the only one G1 left) | **built and empty at this length** | `XTWIST` −0.23× / 0.36×; its own achiral twin beats it |
| node-level candidate-graph observables | **CLOSED** | C2: fully absorbed by consensus |
| passing S30's coherence bar **by ranking** | **CLOSED by derivation** | M5: `Σa=1` passes μ through with coefficient 1 |
| **multi-structure free energy (true barriers)** | **OPEN, on price only** | 10³–10⁵ trajectories/target |
| **pair-level** graph quantities | **OPEN** | C1 did not fire; needs a pair-level operator |
| chiral functionals at 40+ residues | **still open, still unaskable here** | the theorem is length-free, the emptiness is not |

## S31-L11 -- **`bestm128 = 2.9027` IS AN ORDER STATISTIC, AND THE PREFIX AXIS IS A *WORSE-THAN-ARBITRARY* 7-BIT INDEX.** A MATCHED RANDOM-SUBSET FAMILY OVER THE SAME TOP-128, SAME OPERATOR, SAME K, SAME SIZE DISTRIBUTION, BUYS **-0.4191 A** AGAINST THE PREFIX FAMILY'S **-0.2879 A** -- **146% OF IT.** BOTH PRE-REGISTERED BARS FIRE (2026-09-21 00:17, F)

**Pre-registered** in `s31/PREREG_S31_F.md` §11 (commit `de852bef`), **before any aggregate of this
family existed**, with two bars written so they could fire. Both fired.

Artefacts: `s31/s31_F3_prefix.py`, `s31/results/s31_F3_prefix.json`,
`s31/results/s31_F3_randfamily_rows.jsonl` (126 rows). Seed 31007. Source of the 126x128 matrix:
`s29/results/s29_O_p128_rows.jsonl`.

### 0. Provenance gate, passed exactly

My recomputation of the score-ordered prefix curve reproduces S29's stored `curve128`
**bit-for-bit on all 126 targets: max |curve - stored| = 0.0e+00.** So everything below is about
the same object S29 measured, not a lookalike.

### 1. The number is real on the endpoint basis, and I confirm it first

| arm | basis | mean | vs production | MDE | x MDE | fold CI | W/L | verdict |
|---|---|---|---|---|---|---|---|---|
| `bestm128` vs production | **BUILT CHAIN** | **2.9027** | **-0.3079** | 0.0982 | **-3.14** | [-0.3625, -0.2496] | 117/9 | **MEASURED** |
| `bestm128` vs m=75 | CA point cloud | 2.7605 | -0.2879 | 0.0875 | -3.29 | [-0.3312, -0.2389] | 126/0 | MEASURED |

**ORACLE / NOT DEPLOYABLE**, both rows. And one structural fact that must travel with every
future quotation of 2.9027: **it is a CLOUD-SELECTED oracle, projected** --
`s29/s29_O_ladder.py:542` picks `m_best` on the cloud curve and projects that one structure.

### 2. THE MATCHED CONTROL, AND IT INVERTS THE READING (contract rule 7)

The registered control holds **everything** fixed except the one thing under test: the same
top-128, the same operator (superpose on the subset medoid, uniform coordinate mean), the same
K = 128 variants, the same variant-size distribution 1..128 -- but each variant is a **random
subset of that size** instead of the **score-ordered prefix**. The variant index therefore carries
no score information and nothing else changes. Four independent draws, and the **mean** over draws
is reported, never the maximum (contract rule 10).

```
per-target min over the 128 score-ordered PREFIXES     -0.2879 A   (this is bestm128)
per-target min over 128 RANDOM SUBSETS, same sizes     -0.4191 A   sd over 4 draws 0.0065
                                                       -------------------------------
share of the prefix gain reached by the null family     146%
```

**Registered bar: "the m axis is an order statistic" fires if the random family reaches >= 80%.
It reached 146%.** The prefix ordering does not merely fail to beat an arbitrary 7-bit index over
the same set -- **it loses to one, by 0.13 A.**

### 3. The mechanism, measured, so this is not just a null

| family | within-target sd | lag-1 autocorrelation of the curve VALUES | local minima per curve |
|---|---|---|---|
| score-ordered prefixes | 0.1812 | **0.9172** | 25.3 |
| random subsets | 0.1559 | **0.1259** | 42.0 |

The random family has **lower** dispersion and still a **larger** per-target minimum. The
difference is **correlation**: prefix variants are *nested* -- `curve[m]` and `curve[m+1]` share m
members -- so the 128 of them are ~0.92-autocorrelated and contribute far fewer effectively
independent draws than 128. Decorrelate the family and the minimum grows. **That is the
order-statistic mechanism itself, measured directly.**

And the growth curve settles it. Minimum over a random k-subset of the columns, gain vs m = 75,
CA point cloud, 200 repetitions per cell:

```
k          1        2        4        8       16       32       64      128
prefix  +0.0266  -0.0618  -0.1282  -0.1782  -0.2118  -0.2414  -0.2640  -0.2879
random  +0.0348  -0.0377  -0.0824  -0.1317  -0.1854  -0.2547  -0.3325  -0.4241
```

Monotone, unsaturated, in both families. **A quantity that keeps growing with K and has not begun
to flatten at K = 128 is a best-of-K, not a ceiling.** Note also `k = 1`: a *randomly chosen* m is
**worse** than the shipped m = 75 (+0.0266), which is the same statement from the other end.

### 4. And nothing of it transfers -- three independent ways

| arm | basis | effect vs m=75 | MDE | x MDE | verdict |
|---|---|---|---|---|---|
| ORACLE **global** m (m* = 72, one m for all 126) | cloud | -0.0018 | 0.0082 | -0.22 | NOT A RESULT |
| **leave-fold-out** global m | cloud | **+0.0079** (worse) | 0.0260 | +0.30 | NOT A RESULT |
| `stats_lib.split_half_transfer` over the 128 columns | cloud | transfer **-0.0039** of a -0.3179 oracle = **1.2%**, CI [-0.0261, +0.0390] spanning zero | | | |

Best **fold-held-out native-free per-target m-rule**, over four native-free features with a
two-parameter monotone rule fitted per fold:

```
n_distinct   -0.0221   0.51x MDE   65W/61L   fold CI [-0.0477, +0.0064]     NOT A RESULT
DISP128      -0.0039   0.09x MDE                                            NOT A RESULT
n            -0.0054   0.11x MDE                                            NOT A RESULT
rg_sd128     +0.0301   0.61x MDE   (wrong sign)                             NOT A RESULT
```

**Registered bar: "no part of it is deployable" fires if the best rule is > -0.7 x its own MDE.
The best rule is -0.51x. It fired.** And note the *fit itself* touches the native, so even that
-0.0221 is an upper bound on a deployable rule, not a clean one.

### 5. What this changes

**`2.9027 A` must never again be quoted as "the architectural ceiling" or as "what the one integer
T1 says the state can specify is worth".** The defensible sentence is:

> **ORACLE / NOT DEPLOYABLE: a per-target minimum over 128 nested prefix averages of the top-128
> reaches 2.9027 A built chain. A per-target minimum over 128 *arbitrary* subsets of the same set
> reaches 46% further, the gain is still growing at K = 128, and the transferable content of the
> m axis is 1.2% of it with every deployable rule below 0.7x its own MDE.**

The corollary is the useful part. S29-L44 observed that `best1_top128` (choose a member, 7 bits,
-1.0657) beats `bestm128` (choose m, 7 bits, -0.3084) by 3.5x at equal information cost, and read
it as *the architecture spends its bits on the wrong question*. **This entry sharpens that: it is
not that m is a less valuable question than membership. It is that m is not a question at all** --
an arbitrary 7-bit index into the same set outperforms it. The whole -0.3084 is the price of the
minimum, and the axis contributes nothing.

Two prior odds are settled on the record. The coordinator registered ~2:1 that a substantial
fraction is best-of-128; I registered ~4:1 and 8:1-against-deployability, in `PREREG §11.4`, before
measuring. **Both of us were directionally right and both of us were too generous: the surviving
fraction is not "40%", it is negative.**

### 5b. TWO ANNOTATIONS APPENDED IN PLACE AT 00:25, ORIGINAL WORDING LEFT STANDING (rule 13)

**(i) A seeding defect in my own control, found and repaired by me.** The first version of
`run_randfamily` seeded its draws with `hash(pdb)`, which Python salts per process
(`PYTHONHASHSEED`), so the *specific* random subsets could not be regenerated in another process.
The statistic was unaffected -- it is a mean over four draws of a matched family -- but the artefact
was not reproducible, which is not acceptable for a control that carries a headline. Fixed to
`crc32(pdb.encode())` and **the whole control was re-run from zero**. The original rows are kept at
`s31/results/s31_F3_randfamily_rows_hashseed.jsonl`.

```
                                     hash-seeded (original)   crc32-seeded (reproducible)
gain, prefix family                        -0.2879                    -0.2879
gain, matched random family                -0.4191                    -0.4279
  per draw                    -0.4241 -0.4246 -0.4110 -0.4168   -0.4232 -0.4016 -0.4206 -0.4662
share of the prefix gain                     146%                       149%
prefix curve vs S29's stored curve         0.0e+00                    0.0e+00
```

**Both versions fire the 80% bar by more than 60 points.** The headline is quoted from the
reproducible run: **149%**, draw sd 0.0273 over 4 draws.

**(ii) The coordinator's correction, and it makes this entry SMALLER and the record BETTER.** I was
briefed on the premise that the split-half transfer had never been run on `bestm128`. **It had
been, in the entry that produced the number.** S29-L30's own heading reads *"ITS TRANSFERABLE PART
IS ZERO -- THE ORACLE GLOBAL PREFIX IS m = 72 (WORTH -0.0018 A) AND THE LEAVE-FOLD-OUT PREFIX IS
+0.0079 A WORSE THAN PRODUCTION AT 0.30x MDE"*, and S29's report lists *"A transferable prefix
length m"* with verdict **FALSIFIED**. `grid-oracles-are-order-statistics` was **honoured, not
violated** -- and my section 4 above independently re-derives S29's own numbers to the digit, which
I did not realise at the time I wrote it.

So this entry is **not** the discovery of a missed audit. It is a **second, independent support for
a conclusion S29 already reached by a different route**: S29 showed the gain *does not transfer*;
this entry shows the *axis is worse than an arbitrary 7-bit index at the same budget*. Different
evidence, same conclusion, and the honest framing is one finding with two supports.

One further correction to my own section 5, which over-reached: **order-statistic inflation makes an
ORACLE number optimistically biased, and an optimistically biased upper bound is still a valid
upper bound.** So S29's *"2.5 A is unreachable through this architecture"* is **safe, and safer than
it was stated.** What was never licensed is reading 2.9027 as a -0.3079 A lead available to a better
selector. My sentence *"must never again be quoted as the architectural ceiling"* is **too strong
and is withdrawn**; the correct restriction is the one in the display block above it, which
restricts what the number licenses rather than forbidding its use as a bound.

### 5c. OPEN, AND OWNED

The **cloud-to-chain price on this rung is +0.1422** (2.7605 cloud, 2.9027 chain), and S29-L30's
transfer arms are **point-cloud**, so *"the transferable part is zero"* is currently quoted across
two bases -- which the charter forbids and which this sprint has been bitten by three times.
`s31/s31_F3_chain.py` runs the ORACLE-global-m and leave-fold-out-m arms, plus the matched random
family, **on the built chain in one job from the same stored clouds**. Queued behind `AVG_SEP`.

### 6. What is NOT claimed

* **The random-family control is CA point cloud only.** `F3-a` is confirmed on the built chain
  (-0.3079, 3.14x MDE) but the 146% figure is a within-basis cloud comparison. A same-job chain
  version of the control is the one open piece and is queued.
* This says nothing about whether *some other* per-target readout decision is valuable.
  `best1_top128` at -1.0657 is a different arm and is untouched here.
* Nothing here is a deployable gain. Lane F emitted **12 comparisons** in this family; the
  registered family count is 7. Appended to `s31/MULTIPLICITY.md`.

---

## S31-L12 -- **THE COLLAPSE TO ONE SOURCE IS TRUE AND, AS AN INFORMATION STATEMENT, VACUOUS**: ANFINSEN MAKES `I(N;sequence) = H(N)`, SO THE DPI BOUND READS "THE POOL CONTAINS AT MOST EVERYTHING" AND CANNOT LIMIT ACCURACY -- AND **`bestm128 = 2.9027` SURVIVES AS A BOUND BUT ITS "LEAD" FRAMING WAS FALSIFIED IN S29 ITSELF** (2026-09-21 00:19, L)

Requested by the coordinator, who asked to be attacked rather than agreed with. Full working:
`s31/LIT_L.md` §L2.2, §L3.1.

### PART 1 -- THE PREMISE SURVIVES, THE CONCLUSION DOES NOT

**Premise, upheld.** `pool = f(sequence; library)` with a universal library is processing, so
`I(N;pool) <= I(N;sequence)`. A codebook indexed by a sequence-derived key carries no target
argument. **I could not break it.** Enumerating sources is closed: no field, energy, graph
statistic or free energy is a second source.

**Caveat the premise needs:** the library is NOT independent of the native for every target.
`containment-threshold-is-at-the-null` records **4/126 targets carrying a verbatim copy in
their own distogram's training set**. For those `N -> S -> pool` is broken and the DPI does not
apply. That is **LEAKAGE, NOT A CHANNEL** -- it argues for excluding those 4, never for
counting the library as a source.

**Conclusion, REFUTED.** By Anfinsen the native is a function of the sequence, so

```
I(N ; sequence) = H(N)      -- the sequence already contains ALL of it
```

The DPI bound therefore reads **"the pool contains at most everything"**: true, and empty. It
places **no constraint whatever on achievable accuracy**.

> **If "all target-specific information is the sequence" implied a ceiling, AlphaFold would be
> impossible.** A different computation on the same single source extracts far more than ours.

**And this project has its own internal counterexample.** `esm-adds-nothing-for-short-peptides`
was OVERTURNED: on SELECTION, ESM buys **0.288 A over one-hot (p = 0.005, n = 126)** -- a pure
computation on the sequence beating a weaker computation on the same sequence, no new source,
no new measurement. Under the strong reading that gain could not exist. It does.

### WHAT THE ARGUMENT ACTUALLY LICENSES

**NOT** "you need a measurement of the molecule, not a computation" -- that is unlicensed and
adopting it would wrongly close the two escapes the field actually used. What is licensed:

> **No re-combination of the objects we currently hold can help. The remaining moves are a
> better ESTIMATOR or a better LIBRARY -- not a further source.**

### THE DISTINCTION THE COLLAPSE INVITES US TO LOSE

The library plays **two roles** and the enumeration conflates them:

| role | verdict |
|---|---|
| as an information SOURCE | a constant, no target argument -- **not a source**. Coordinator is right |
| as an architectural CEILING | **binding and measured**: ORACLE distances still give only ~1.95-2.0 A through this library |

**"Not a source" does not mean "not a constraint."** The ~2.0 A figure is a LIBRARY statement,
and it survives perfect information -- which is exactly what an informational argument can
never produce. The two constraints that actually bind are both estimator/architecture facts:
(1) in-band ordering **0.600** across targets vs **0.638** needed for 2.0 A; (2) ORACLE
distances capping at **~1.95-2.0 A** through this library.

### PART 2 -- AUDIT OF `bestm128 = 2.9027 A` (coordinator's §5)

**The transfer arm EXISTS, in the entry that produced the number.** `s29/LEDGER.md:2777`,
S29-L30 heading: *"ITS TRANSFERABLE PART IS ZERO -- THE ORACLE GLOBAL PREFIX IS m = 72 (WORTH
-0.0018 A, i.e. THE SHIPPED 75) AND THE LEAVE-FOLD-OUT PREFIX IS +0.0079 A WORSE THAN
PRODUCTION AT 0.30x MDE."* The 33-arm ladder (`s29/LEDGER.md:4288-4290`) carries a dedicated
**`transferable`** column saying the same on all three prefix-m rungs, and
`s29/REPORT_S29.md:252` lists *"A transferable prefix length m"* with verdict **FALSIFIED**.

So the per-target ORACLE gain of **-0.3079 A is ~99.4% order statistics**, and S29 measured
that and said so. **`grid-oracles-are-order-statistics` was HONOURED here, not violated.**

**THE NUMBER DOES NOT DEFLATE. ITS USE DOES.** The intuition runs the wrong way:

> Order-statistics inflation makes an ORACLE number **optimistically biased** -- and an
> optimistically biased **UPPER BOUND IS STILL A VALID UPPER BOUND**.

So S29's actual inference -- *"2.5 A is unreachable through this architecture"* -- is **safe,
and safer than stated**: the true ceiling is WORSE than 2.9027, which strengthens
unreachability. **The "architectural ceiling" framing was licensed.**

What is **NOT** licensed is reading 2.9027 as *"a -0.3079 A lead available to a better
selector."* S29 priced that at **-0.0018 A global / +0.0079 A LFO** and marked it FALSIFIED. If
this sprint treats "close the gap to 2.9027" as its most promising remaining lead, that lead
was measured at **0.6% of face value two sprints ago** and the caveat was dropped in
re-quotation -- exactly the L3 pattern.

### THE ONE THING GENUINELY STILL OPEN, FOR LANE F

The transfer numbers were measured on the **POINT CLOUD** (S29-L30's arm is 2.7605 vs
production 3.0483, both cloud). **2.9027 is the BUILT CHAIN**, and the cloud->chain price on
that rung is **+0.1422**. So *"the transferable part is zero"* is quoted **ACROSS BASES**.
Very likely fine, but not literally measured on the chain, and the charter's basis discipline
forbids silent cross-basis transfer.

> **Lane F: do not re-derive the ladder -- S29 already did. Pin the one open thing: run the
> global-m and leave-fold-out-m arms ON THE BUILT CHAIN.** A two-arm measurement, not a
> sprint. It either closes the last gap in a two-sprint-old claim or finds the only place it
> could move.

## S31-L13 -- **THE EXACT OBJECTIVE OF EVERY AVERAGING READOUT IS AN IDENTITY: ||C-t||^2 = <w,a> - 1/2 w'Bw, WITH THE NATIVE-FREE HALF PAIRWISE AND **REPULSIVE**. SO THE READOUT AND RANKING QUESTIONS ARE ONE PROBLEM AND THE WHOLE DEFICIT IS THE PER-CANDIDATE QUALITY `a`: CROSSING PRICE rho = 0.211 AGAINST THE SHIPPED SCORE'S 0.1176, ORACLE CEILING 1.829 A AGAINST PRODUCTION 3.048. THE DEPLOYABLE ARMS ARE **WORSE** (+0.133 A, 1.17x MDE, 5/5) AND THE SHUFFLED-B CONTROL FIRES ON THE ONE POSITIVE. R1 IS HALF FALSIFIED AND THE SURVIVING HALF CAPS CAPACITY AT **6.886 BITS, NOT 7** (2026-09-21 00:21, A)

**Basis: CA POINT CLOUD throughout (production 3.0483 Å). Nothing here is on the built chain
(production 3.2105 Å).** Fold-clustered SE on the pinned 5 folds, MDE = 2.8016 × SE.
Pre-registration `s31/PREREG_S31_A.md`, committed inside 8ce5e1a0 before any number existed.
Full working `s31/THEORY_A.md`. Code `s31/s31_A_r1.py`, `s31_A_cap.py`, `s31_A_readout.py`.
Artefacts `s31/results/s31_A_{r1,cap,readout}.json` + `*_rows.jsonl`.

**Instrument validated before anything was claimed:** my reconstruction of production's uniform
DIS-top-75 average gives **3.048338 Å** against the canonical **3.0483**, and the shipped score
argmin gives **3.4540** against the 3.454 asserted in `core/quantum.py`.

---

### 1. THE IDENTITY — the exact objective of every averaging readout, and it forces the sign

For **every** `w` with `Σ_x w_x = 1` (non-negativity NOT required, so it covers the selection,
convex and affine readouts in one formula) and any fixed frame:

```
|| Σ_x w_x W_x − t ||²_F   =   ⟨w, a⟩  −  ½ w' B w        EXACT
   a_x  = ||W_x − t||²_F     ORACLE       (per-candidate squared error)
   B_xy = ||W_x − W_y||²_F   NATIVE-FREE  (and ½w'Bw = tr Σ_w, the weighted dispersion)
```

Verified to **1.66e-11 max relative error** over 126 targets × 80 draws, on simplex draws *and*
affine draws with negative weights. Four consequences, all derived, not measured:

1. **The readout question and the ranking question are the same problem, with an equals sign.**
   `B` is free and exact; `a` is the entire deficit.
2. **The native-free half carries a MINUS sign** — at fixed quality the readout should
   *maximise* weighted mutual spread. The brief's `H = diag(zrank) − λ·W(similarity)` with
   `λ > 0` is attractive; the derivation says repulsive.
3. **The attractive branch cannot produce a distribution at all.** `B` is a squared-distance
   matrix, hence conditionally negative definite, so `w ↦ w'Bw` is **concave on the simplex**;
   the objective is convex for `γ > 0` and **concave for `γ < 0`**, and a concave function on a
   polytope is minimised at a **vertex** — i.e. the negative branch degenerates to the shipped
   argmin, **by theorem**, and no tuning could have rescued it.
4. **The forced operator is mean-field, not a Hamiltonian.** `∂/∂w_x = â_x − (Bw)_x` is the VMC
   local energy of `H[w] = diag(â) − B`; the energy is **quartic in ψ**, so there is no per-shot
   eigenvalue and no CVaR-VQE. Independent confirmation of lane L's S31-L4 closure from the
   readout side.

> **Caveat I am enforcing on myself: the deployed `Pt` is NOT `B`.** `Pt` is pairwise RMSD with
> *per-pair* optimal superposition; `B` needs one *common* frame. Substituting `n·Pt²` breaks
> the identity at **2.39 % median relative error**. Small, real, and not licensed.

### 2. THE DERIVED READOUT IS MEASURED AND IT IS WORSE THAN PRODUCTION

| arm (cloud, n = 126) | mean | status |
|---|---|---|
| shipped score argmin | 3.4540 | deployable |
| `GAM(LFO)` — one LFO scalar γ | 3.1942 | deployable |
| **`MEB` = `argmax_Δ ½w'Bw`** — quality-blind, **zero parameters, no score** | **3.1919** | deployable |
| `CAL` — LFO-calibrated `â`, γ = 1 exactly | 3.1817 | deployable |
| best cell of the whole γ grid (γ = 2) | 3.1687 | not a selection rule |
| **`PROD75`** | **3.0483** | production |
| ORACLE convex QP over the simplex | **1.8290**, support 6.54/128 | **ORACLE — NOT DEPLOYABLE** |
| ORACLE affine (readout 2) | **0.0000**, rank 32.9 | **ORACLE — NOT DEPLOYABLE** |

| comparison | mean | SE | × MDE | folds | verdict |
|---|---|---|---|---|---|
| **PRIMARY `CAL − PROD75`** | **+0.1334** | 0.0409 | 1.17 | 5/5 | **WORSE** |
| `GAM(LFO) − PROD75` | +0.1459 | 0.0471 | 1.10 | 5/5 | WORSE |
| `MEB − PROD75` | +0.1436 | 0.0420 | 1.22 | 5/5 | WORSE |
| **`MEB − argmin(score)`** | **−0.2621** | 0.0342 | 2.74 | 5/5, 81W/45L | **BETTER** |
| `γ=1 − shuffled-B` (8 draws) | −0.0250 | 0.0165 | **0.54** | — | **NOT A RESULT** |
| `γ=1 − shuffled-score` | −0.0500 | 0.0306 | 0.58 | — | **NOT A RESULT** |

> **The one positive, demolished by its own control.** Quality-blind dispersion maximisation
> beats the shipped argmin by **−0.2621 Å at 2.74× MDE**. But replacing `B` with a **random
> relabelling of itself** costs only **0.0250 Å at 0.54× MDE — NOT A RESULT**. The mechanism is
> *"spread the weights over many candidates"*, not *"spread them along the real geometry"*.
> **No part of the −0.2621 Å may be attributed to the pairwise structure the derivation is
> about.** That is `operator-consumes-set-mean` again, and the control was placed at the
> decisive step precisely so it could fire.

### 3. THE PRICE OF `a` (ORACLE sweep — NOT DEPLOYABLE)

`â` interpolated from the shipped DIS z-rank toward the true `a`, γ at the derived value 1:

| ρ(`â`,`a`) | 0.118 | 0.330 | 0.604 | 0.853 | 0.974 | 1.000 |
|---|---|---|---|---|---|---|
| cloud RMSD | 3.170 | 2.893 | 2.614 | 2.332 | 2.027 | **1.829** |

> **The shipped DIS score supplies ρ = 0.1176.** The derived readout crosses production at
> **ρ = 0.211**, reaches **3.00 Å cloud at ρ = 0.248**, and 2.50 Å cloud at ρ = 0.705. Beating
> production costs **1.79× in ρ (3.2× in ρ²)**; 3.00 Å cloud costs **2.11× in ρ (4.4× in ρ²)**.

This ρ is **a different quantity** from S30 THEORY §8.3's `cos(u,e)`; the closeness of 0.1176 to
that document's 0.1128 is a coincidence of two different objects and must not be quoted as two
instruments agreeing.

### 4. R1, ATTACKED — half survives and is stronger; half is falsified

* **Falsified (point 4, and it makes capacity *lower*).** `argmin(P e_j) = j` needs `P[i,j] > 0`
  for `i ≠ j`; duplicate candidates give `P[i,j] = 0` and `np.argmin` returns the first.
  Reachable vertices **118.45 mean, 94 min**, and that count **equals the byte-distinct
  structure count on every one of the 126 targets** — mechanism confirmed, not inferred.
  `filter_pool` already dedups (`core/pipeline.py:766-773`); the readout does not.
  **Corrected capacity: log₂(118.45) = 6.886 bits mean, 6.555 worst — not 7.**
* **Falsified (points 1-3) for the arm that emits structure.** The cited `:869-871` are
  `quantum_stage`'s *selector*. `average_weighted` (`:880-895`) uses the medoid only as the
  superposition **frame** (`:890`) and emits a continuous convex combination at **`:894`**. The
  emitted structure is **1.1144 Å (min 0.0730) from the nearest pool member**. `p` enters
  through a ≤7-bit piecewise-constant frame **and** a continuous (D−1)-dimensional weight
  vector; only the frame is capped.
* **Realised capacity, ORACLE bits delivered (S30 currency) — ORACLE, NOT DEPLOYABLE:**
  `p_θ`-weighted medoid **1.805**, `p*`-weighted **1.876**, score argmin **1.693**, uniform
  medoid **1.659**. The stage moves the selection on **77.8 %** of targets and delivers
  **0.112 bits** over the plain score argmin — **1.6 % of the 6.886-bit alphabet** — and the
  closed form delivers more of them than the circuit does.

### 5. EVERY DEPLOYED-QUANTUM CLOUD COMPARISON AT n = 126 IS BELOW 0.4× MDE

| comparison | mean | SE | × MDE | verdict |
|---|---|---|---|---|
| quantum synthesis − PROD75 | +0.0178 | 0.0170 | 0.37 | NOT A RESULT |
| quantum synthesis − uniform-128 (matched set) | +0.0125 | 0.0165 | 0.27 | NOT A RESULT |
| sel(p_θ) − sel(uniform medoid) | −0.0308 | 0.0321 | 0.34 | NOT A RESULT |
| quantum synthesis − `p*` synthesis | −0.0044 | 0.0177 | 0.09 | NOT A RESULT |

**Correction to my own message to the coordinator:** I quoted the +0.0178 and +0.0125 without
their MDEs, which made them read as measured costs. They are **NOT A RESULT**. The honest
statement is that on the cloud basis the deployed quantum stage is **indistinguishable** from
every classical alternative it competes with, including the closed form that replaces it.

### 6. A1 — THE DEPLOYED OBJECTIVE'S OPTIMUM IS ONE SCALAR (independent of lane L's S31-L4)

`p*_x ∝ exp((μ − E_x)_+/(αT))`, μ fixed by `Σ_{E_x<μ} p*_x = α`. Derived by simplex KKT +
Rockafellar–Uryasev; lane L reached the same object by Sion duality. Two independent routes.
At **α = 1 — three of the five pinned folds — `CVaR₁ = ⟨E,p⟩` and the CVaR is doing nothing at
all**; the law collapses to plain Boltzmann and agrees with mirror descent to **1.4e-15**.

Measured against the circuit: `KL(p_θ‖p*) = 0.930 bits` mean (median 1.298), TV 0.378, and
`F(closed) ≤ F(numeric)` at **+1.8e-15**. Yet the emitted structures differ by **−0.0044 Å at
0.088× MDE** while differing **0.102 Å per target in absolute value**. **The optimisation gap is
real and cancels in the mean** — the sharpest available statement of "the readout discards what
the objective computes".

### 7. THE ANSATZ'S INDUCTIVE BIAS (coordinator's redirect (b))

Schmidt ranks of `ψ(θ)` across the six cuts: **2, 4, 8, 8, 4, 2** — an **MPS Born machine of
bond dimension exactly `2^layers = 8`**, which pins lane L's Han-et-al. framing to this ansatz
rather than to Born machines in general. R² of `log p` on the objective's own shape (the hinge
`(μ−E)₊`, on which `log p*` scores 1.000 by construction): **random θ 0.0071, best of 400
draws 0.0538, optimised `p_θ` 0.467.**

> **One sentence: the reachable set starts essentially orthogonal to the objective's shape,
> optimisation carries it to half of it, and it stops there.** That is not a bias *toward*
> anything structural — it is a ceiling at about half the right shape, which is lane L's
> expressivity floor in a different currency.

**My own explanation of the residual, refuted.** I predicted a *basis* mismatch (a product
state's log-probability is additive over bits, hence a function of `popcount(x)`; `log p*` is a
hinge in the rank's *value*). The registered native-free intervention — relabel so rank `i` goes
to the `i`-th bitstring in `(popcount, value)` order — made it **worse**: KL 0.967 vs 0.930,
cloud **−0.0031 Å at 0.060× MDE, NOT A RESULT**; popcount adds 0.016 of R² beyond the hinge.
**I have no replacement explanation and am not fitting one after the fact.**

### 8. DEFECT FIXED: `core/quantum.py` free_energy docstring, on the real instrument

Asserted: *"for ANY alpha the minimiser concentrates p on the lowest-energy basis states … state
entropy 0.01 bits at alpha=1 … and it is a property of CVaR, not of the optimiser."* Measured at
`T = 0`, 8 seeds/target, deployed per-fold α:

| | n | entropy at `T = 0` |
|---|---|---|
| α = 1 (folds 0, 3, 4) | 78 | **0.258 bits** |
| α = 0.25 (folds 1, 2) | 48 | **3.596 bits** |

At α = 1 the minimiser is the unique argmin vertex and the collapse **is** a property of CVaR.
At α < 1 the argmin set is `{p : p_x₀ ≥ α}`, a face of **positive volume**, so the objective
does not determine `p` at all and the optimiser's path picks the point — **no collapse, and the
entropy belongs to the optimiser**, the opposite of the sentence. The quantifier "for ANY alpha"
is **false by theorem**. Corrected in place with the original quoted, as lane D did for
`pipeline.py:821`.

### 9. WHAT THIS CLOSES

1. **Q1** — escaping a fixed order needs a `λ`/`ψ`-dependent gradient; non-diagonal `H` supplies
   one only under the VMC local energy, which is not a measurement. Under computational-basis or
   spectral CVaR a **generalised prefix theorem** applies (A2).
2. **Q2 as a quantum question** — mean-field, quartic in ψ, no per-shot eigenvalue; and a
   candidate-index register's Hilbert dimension *is* the candidate count, so every operator on
   it is one `eigh` away. Agrees with lane L by two independent arguments.
3. **Q2's proposed sign** — the consensus/attractive branch degenerates to the shipped argmin,
   by theorem.
4. **The derived objective as a deployable readout** — 0.12–0.15 Å worse than production at
   1.1–1.2× MDE, 5/5 folds, with the shuffled-B control showing the pairwise content is not
   producing even the part that works.
5. **Q3** — every construction in this lane is classically reducible: A1's optimum is one
   bisection (`O(D log 1/ε)`), any diagonal `H` is one sort, any non-diagonal `H` on a candidate
   register is one `eigh` (`O(D³)` ≈ 1 ms at D = 512), the mean-field objective is a convex QP.
   **Nothing here is a quantum mechanism, and the reason is the encoding, not the Hamiltonian.**

### 10. WHAT IT OPENS — the only thing in this lane worth compute

> The exact objective is `⟨w,a⟩ − ½w'Bw`. **Half is free and exact; the entire deficit is `a`.**
> ORACLE ceiling of the convex readout **1.829 Å** against production 3.048, crossing price
> **ρ(â,a) = 0.211** against the shipped **0.1176**.

Two consequences for other lanes:

* **Charter §7C (sparse weighted readout): the exact optimum is already sparse — 6.54 of 128
  members with no sparsity penalty imposed.** Sparsity is an *output* of the correct objective,
  not a design choice; and any sparse-readout work that does not improve `â` is spending effort
  on the half that is already exact.
* **Readout 2's ORACLE ceiling is exactly 0 Å**, by rank (`rank(aff{W_x}) = 32.9 ≥ 3n−3` on all
  126; measured residual 0.0000 on every target). It is an over-parameterised interpolator: 127
  weights against ~33 residual dimensions. So **"0.2516 Å under an ORACLE objective" is a
  statement about a regulariser, not a class ceiling** — and the simplex constraint is exactly
  the regulariser the uniform average enjoys for free. That *derives* S31-L2(2) ("expressivity
  without an aligned objective is harmful") rather than observing it.

### 11. REGISTERED PREDICTIONS, SCORED (full table in `s31/THEORY_A.md` §8)

**HELD:** A3-i (identity 1.66e-11); A3-frame (2.39 %); A3-sign (LFO γ > 0 on **5 of 5** folds);
A1-e readout half (0.09× MDE); A1-d `T=0` seed sd (0.162 > 0.15).
**REFUTED / MISSED:** A1-v as written (the statistic measured my reference solver, not the
closed form — my error); A1-e KL (0.930 vs `< 0.10`, **10×**); A1-d deployed-`T` seed sd
(0.096 vs `< 0.05`); **A3 PRIMARY** (+0.1334 vs the registered `[−0.15,+0.10]`); A3-meb
magnitude (+0.144 vs `[+0.3,+1.5]`); **A3-oracle** (1.829 vs `< 1.2`, by 0.63 Å); the popcount
relabel hypothesis.

**Direction of the misses:** two ran *against* my hypothesis (the primary landed worse than my
band; the KL was 10× my prediction) and two ran *toward* it (I put the ORACLE ceiling 0.63 Å too
low, and over-predicted how badly quality-blind MEB would do, which flattered the derivation by
making its failure look inevitable). `A3-frame` is the only clause where being right cost me
something — it forbade substituting the deployed `Pt` for `B`.


## S31-L14 -- **S30's "THE ARGMIN DOMINATES AT EVERY BIT BUDGET" IS A SET-MISMATCH ARTEFACT: ON A FIXED TOP-128 THE CONVEX READOUT IS -0.3450 A BETTER THAN THE ARGMIN (4.38x MDE, 5/5 FOLDS, 119W/4L).** AND THE WHOLE NATIVE-FREE CLASS IS CLOSED BY AN **EXACT IDENTITY**: A SUM-TO-ONE READOUT CAN ONLY ACT THROUGH THE CROSS TERM `2<c,Dw>`, AND EVERY NATIVE-FREE RULE MEASURED CAPTURES **3.96%** OF THE ORACLE'S WHILE PAYING DISPERSION THAT CANCELS IT TO WITHIN **0.01 A^2** (2026-09-21 00:20, C)

Pre-registration `s31/PREREG_S31_C.md` at **8ce5e1a0**, amendment 1 at **2e12e02c**, both before the
first lane-C number. Code `s31/s31_C_ladder.py`, `s31/s31_C_ident.py`. Artefacts
`s31/results/s31_C_ladder.json`, `s31_C_ladder_rows.jsonl` (126), `s31_C_ident.json`.
**Basis: CA POINT CLOUD** except the block headed BUILT CHAIN, which re-uses s29 lane O's already
built chains and projects nothing (S31 operational rule: both sides of a chain contrast must be
projected in the same job).

## 1. THE SET-MISMATCH, WHICH INVERTS A PUBLISHED CLOSURE

S30-L11 closed the sparse weighted readout "by price -- a plain argmin dominates at every bit
budget", and supported it on the chain with **2-of-75 (2.1683) against 1-of-128 (2.1435)**. Those
are different candidate sets, so the comparison is not about the readout. Holding the set fixed at
the deployed top-128 (**ORACLE / NOT DEPLOYABLE**, CA point cloud, n = 126):

```
arm                                rmsd     vs argmin over the SAME 128
prod_top75_uniform               3.0483        +0.9025   4.65x MDE  5/5 folds    3W/123L
unif_prefix_m128 (native-free)   3.0435        +0.8977
ORACLE argmin over 128           2.1458             --
ORACLE convex over 128           1.8008        -0.3450  -4.38x MDE  5/5 folds  119W/4L
ORACLE affine over 128           0.0000        -2.1458   (a tautology -- see section 4)
```

**BUILT CHAIN, from `s29/results/s29_O_chain_rows*.jsonl`, nothing rebuilt:**

```
production (uniform top-75 average)                3.2105
bestm128  -- ORACLE prefix-average over the 128    2.9027    <- T1's ENTIRE reach, 0.308 below prod
best1_top128 -- ORACLE argmin over 128 (7 bits)    2.1435
hull_top128  -- ORACLE CONVEX weights over the 128 1.8538    <- 0.290 better than naming one
best1_pool   -- ORACLE argmin over 500 (9 bits)    1.7078
```

**The honest two-sentence form, and it is the correction:** *as classes at a fixed candidate set,
combining beats naming -- the convex readout is 0.290 A better on the built chain than naming the
single best of the same 128. As an A-per-bit question across sets, naming wins -- the argmin over
the whole 500 reaches 1.7078 for 9 bits where 2-of-128 with free continuous weights needs 12.99
support bits plus an unbounded weight channel to reach 1.9138.* S30 published only the second and
supported it with a set-mismatched pair.

**One correction I owe the coordinator, made against my own message.** I told him S30's
`T128_s2_unif` = 2.0700 was "2-of-128 with uniform weights, so the 0.076 A over the argmin is pure
error cancellation and costs zero weight bits". The weights are indeed free, but that arm's
**support is ORACLE-chosen** (greedy against the native, 12.99 bits, charged by S30). The
increment from 1 to 2 members at matched ORACLE support is error cancellation; the arm as a whole
is not free, and my sentence implied it was.

## 2. THE IDENTITY THAT CLOSES THE NATIVE-FREE CLASS

For **any** weight vector with `sum(w) = 1`, in the deployed common frame:

```
||X(w) - nat||^2 / n  =  c2  +  2<c, Dw>/n  +  ||Dw||^2 / n
        c = Xbar - nat   (Xbar the unweighted set mean)      D = W - Xbar
```

`c2` is a property of the candidate set and **no readout of this form can touch it**: the shared
component is preserved exactly by every convex *and* every affine combination. A readout beats the
set mean only by making the **cross term** negative faster than it pays in `||Dw||^2`. Measured
(A^2, mean over 126, same fixed top-128):

```
                          cross      ||Dw||^2    cross + ||Dw||^2     emitted rmsd
unif_prefix_m128         +0.0000       0.0000        +0.0000            3.0435
soft_T4.00 (native-free) -0.0353       0.0257        -0.0096            3.0391
soft_T2.00               -0.0846       0.0969        +0.0123            3.0383
soft_T1.00               -0.2025       0.3097        +0.1072            3.0442
soft_T0.50               -0.4122       0.7096        +0.2974            3.0638
soft_T0.25               -0.6615       1.2107        +0.5492            3.0983
typ_T1.00 (typicality)   +0.1844       0.2708        +0.4552            3.0957   <- WRONG SIGN
ORACLE argmin128        -10.8999       6.6853        -4.2146            2.1458
ORACLE convex128        -10.4214       3.9984        -6.4230            1.8008
ORACLE affine128        -24.1357      12.0678         0.0000            0.0000
```

**A sum-to-one readout helps if and only if `|cross| > ||Dw||^2`.** The ORACLE convex optimum runs
at `|cross| = 2.61 x ||Dw||^2`. **Every native-free rule runs at approximately 1** -- the score's
orientation toward the native is real and is cancelled, to within **0.01 A^2**, by the dispersion
the concentration costs. Sharpening the softmax buys cross monotonically (-0.035 -> -0.662 as
T falls 4 -> 0.25) and pays `||Dw||^2` faster every time.

**The best native-free rule captures 3.96% of the ORACLE convex cross term.** That single number is
the size of the whole native-free readout channel.

*Frame caveat, stated because the decomposition needs it:* the identity is exact in the deployed
common frame while every RMSD elsewhere is Kabsch-optimal. On the `m = 128` arm the two differ by
**0.0216 A** (3.0651 in-frame against 3.0435 optimally superposed), and that gap is the size of the
approximation in reading the table above as A rather than as A^2 of the frame.

## 3. THE SATURATION LAW THE COORDINATOR PREDICTED -- CONFIRMED AT R^2 = 0.997

Uniform averages over **random** `m`-subsets of the 128 (the exchangeable arm; MEAN of 8 draws,
never the per-target minimum, contract rule 10), fitted as mean-square across targets:

```
MS(m) = c2 + v2/m      c2 = 11.8495    v2 = 4.2120    R^2 = 0.99700
common-mode share at m = 1:  0.7378        (the coordinator predicted 0.68)
fitted floor  c_rms = 3.4423 A            (RMS across targets; the mean-RMSD reading is 3.0435)
```

The law holds, and the floor is the set mean's own error -- which by section 2's identity **no
sum-to-one readout, convex or affine, can go below with native-free weights.** The measured share
0.7378 sits above the project's published 67.6%; the two are measured on different sets (this one
is the top-128, the published one the top-75) and the difference is not interpreted here.

## 4. F-C1d REFUTED -- THE AFFINE CEILING IS A DIMENSION-COUNTING TAUTOLOGY

Registered pre-check (contract rule 22): the ORACLE affine ceiling is informative only if
`ess >= 5` and `neg_mass <= 1` on a majority of targets.

```
ORACLE affine over the 128:  rmsd 0.000000
median ess 1.13 of 128 | median neg_mass 3.38 | median ||w||_1 7.76 | median 3n 39
targets with ess >= 5 AND neg_mass <= 1:  1.6%
```

**REFUTED, as I registered I expected.** 128 generic windows span a 39-dimensional coordinate
space, so the affine hull contains the native exactly, reached by cancelling weights with an
effective sample size near one. `s27/s28_A_FINDINGS.md:85-90` already said not to quote an
affine-hull ceiling as a bound; this is that statement measured on this set.

## 5. THE NORM-BOUNDED MIDDLE HAS NO KNEE, SO THE CONSTRAINT CANNOT BE CHOSEN BY CEILING

The registered ladder rung: `w = argmin ||Aw - nat||^2 + lam*||w - u||^2` s.t. `sum(w) = 1`,
closed form on the zero-sum subspace. **ORACLE / NOT DEPLOYABLE.**

```
lam      1e4     3162    1000     316     100      31.6     10      3.16     1      0.01    1e-6
rmsd   2.7318  2.4044  1.9957  1.5816  1.1977  0.8598  0.5855  0.3881  0.2555  0.0250  0.0000
||w||_1  1.01    1.13    1.47    2.09    3.08    4.49    6.29    8.37   10.82   24.80   28.32
ess    112.3    86.3    56.2    30.9    15.4     7.8     4.5     3.0     2.2     1.5     1.5
neg_m   0.007   0.066   0.234   0.545   1.040   1.746   2.647   3.683   4.911  11.899  13.662
```

**Smooth and monotone; there is no interior structure to find.** At `||w||_1 = 2.09` -- 0.545 of
negative mass, ess 31 -- the ORACLE already reaches **1.5816 A, better than the argmin over the
whole 500 (1.7108)**. So the coordinator's reframe ("what is the right constraint set on an
otherwise too-expressive affine readout") **has no ORACLE answer**: the ceiling is a continuous
function of the norm budget. A constraint here must be chosen by *generalisation*, which is
section 6.

## 6. F-C1c REFUTED -- THE ORACLE WEIGHTS ARE NOT IDENTIFIABLE FROM NATIVE-FREE FEATURES

Target: the ORACLE convex weight vector minus uniform. Features: 7 native-free per-candidate
quantities (z-scored DIS, rank, typicality, distance to the set medoid, Rg, virtual-bond length,
score gap) plus a constant. Ridge, **leave-fold-out on the pinned folds**.

```
out-of-fold R^2  0.02103   against a registered 5% bar   ->  REFUTED
APPLIED (contract rule 21, never quote an implied conversion alone):
   fitted weights   3.0792      against the score-prefix-75 on the same set   3.0483
   effect           +0.0309     0.37x MDE, 60W/66L -- i.e. WORSE, and below its own MDE
   ORACLE convex on the same set                                             1.8262
```

The S30 lane-P pattern exactly: the predictable part of the weight vector is not the part that
pays. **F-C1a is REFUTED** (best native-free arm `soft_T2.00`, **-0.0101 A at 0.56x MDE, NOT A
RESULT**, and the record already had -0.0038 from a 30-arm sweep and a leakage-ORACLE grid that
selects uniform). **F-C1b is CONFIRMED**: no native-free arm moves the bias by 0.10 A; all of the
movement is dispersion, or a cross term two orders of magnitude below the ORACLE's.

## 7. WHAT WOULD IT TAKE -- the axis priced, ORACLE / NOT DEPLOYABLE

Stepping a fraction `f` along the ORACLE convex direction, `w(f) = u + f*(w_or - u)`:

```
f      0.00   0.10   0.20   0.30   0.40   0.50   0.70   1.00
rmsd  3.0651 2.9232 2.7845 2.6495 2.5189 2.3934 2.1631 1.8994
f needed for 3.00 A: 0.06      f needed for 2.50 A: 0.42
```

**This is not comparable to the 3.96% cross capture and must not be read as "4% against 6%".** The
native-free rules obtain their cross term in a *different direction* that costs far more in
`||Dw||^2`; `f` prices the ORACLE direction only. The row is here to give the axis a scale, not to
suggest the two are 2 percentage points apart.

## CONVERGENCE WITH TWO LANES THAT LANDED WHILE THIS RAN

**Lane A (S31-L13)** derives the same object from the other side -- `||C-t||^2 = <w,a> - (1/2)w'Bw`
with the native-free half pairwise and **repulsive** -- and reports its deployable arms **worse**
(+0.133 A, 1.17x MDE, 5/5 folds) with the shuffled-`B` control firing on its one positive. Section
2 here is the same conclusion in the `c2 / cross / ||Dw||^2` basis, reached from a different
construction and with a different control. **Two independent derivations agreeing is worth stating,
and per contract rule 25 it is worth stating that they are the same object and not two.**

**Lane F (S31-L11)** finds `bestm128 = 2.9027` is an order statistic and the *prefix* axis is a
worse-than-arbitrary 7-bit index, a matched random-subset family buying 146% of the prefix
family's gain. S31-L15's independent finding that the best index map for a **deployable** partial
readout is the **random permutation** is the same phenomenon at the index rather than at the
prefix length.

## COMPARISONS MADE (contract rule 23)

Ladder: 11 named arms x 1 aggregate = 11, plus a 21-point m-curve on two supports (prefix and
random) = 42 descriptive, 8 fold-clustered contrasts against the argmin, 33 against production.
Identifiability: 1 out-of-fold fit + 1 applied contrast. Cross accounting: 8 native-free + 3
ORACLE cross terms. Ridge path: 21 lam values, descriptive. **Read as results: F-C1a, F-C1b,
F-C1c, F-C1d and the set-matched inversion -- 5.** No per-target maximum is reported anywhere in
this entry.



### CORRECTION 1 to S31-L14, made by the lane's own self-audit at 2026-09-21 00:31, original numbers left standing above

**The defect.** The `ORACLE_convex128` rung was solved by two iterative solvers (FISTA with simplex
projection, and s29 lane O's alternating convex-NNLS) and the better of the two was taken. **On 4
of 126 targets both returned a point WORSE than the best simplex VERTEX** — which is feasible, so
the returned value could not have been the convex optimum. The targets are **1ID6 (+0.1821),
2NDM (+0.1794), 9L1M (+0.0267), 2BP4 (+0.0031)**, three of them `FAIL18`, i.e. exactly the
ill-conditioned candidate sets.

**The fix.** The best vertex is now a third candidate and the bound is `assert`ed rather than
assumed (`s31/s31_C_ladder.py`, the `ORACLE_convex128` block). Solver used across the 126:
**hull 115, fista 6, vertex 5.**

**Every affected number, old -> new.** The correction moves in the direction that STRENGTHENS the
entry's conclusion, which is why it is stated here in full rather than quietly applied:

```
                                              published        corrected
ORACLE convex over the 128 (cloud)              1.8008           1.7977
   vs the ORACLE argmin over the same 128      -0.3450          -0.3481
   xMDE / fold CI / folds / W-L        -4.38x [-0.4070,-0.2938]  -4.48x [-0.4072,-0.2984]
                                        5/5  119W/4L             5/5  119W/0L
ORACLE convex cross term (A^2)                -10.4214         -10.3193
ORACLE convex ||Dw||^2 (A^2)                    3.9984           4.2054
|cross| / ||Dw||^2                               2.61             2.45
best native-free capture of the ORACLE cross     3.96%            3.99%
f needed for 2.50 A on the ORACLE axis           0.42             0.44
f needed for 3.00 A on the ORACLE axis           0.06             0.06
```

**Nothing qualitative changes**, and the W/L improves from 119W/4L to **119W/0L** because the four
losses *were* the four solver failures. The BUILT CHAIN row `hull_top128 = 1.8538` is s29 lane O's
own artefact and is untouched by this; it carries whatever solver tolerance S29 used, and the
0.290 A chain contrast against `best1_top128` is quoted from S29's rows unchanged.

**Why this is recorded rather than fixed silently.** A ceiling that violates its own feasibility
bound is a solver failure that flatters nothing — but it would have been invisible to any reader,
and the audit that caught it was a one-line check (`is the convex optimum ever worse than a
vertex?`) that should have been in the file from the start. It is now an `assert`.

## S31-L15 -- **RE-INDEXING THE 7-BIT CANDIDATE REGISTER IS WORTH -0.0655 A ON AN ORACLE CEILING AND NOTHING DEPLOYABLE**, AND THE MECHANISM IS A **CONFLICT BY CONSTRUCTION**: STRUCTURE-AWARE INDEXING RAISES THE ORACLE VALUE OF A PARTIAL MEASUREMENT (2.4216 -> 2.2712) AND **LOWERS** ITS DEPLOYABLE VALUE (3.2369 -> 3.3721) -- THE BEST MAP FOR A DEPLOYABLE PARTIAL READOUT IS THE **RANDOM PERMUTATION**. PLUS: GRAY CODING IS A PROVEN NO-OP HERE, AND NO INDEX BIT CARRIES 0.07 BITS ABOUT CANDIDATE QUALITY (2026-09-21 00:20, C)

Pre-registration as above. Code `s31/s31_C_index.py`; artefacts `s31/results/s31_C_index.json`,
`s31_C_index_rows.jsonl` (126). **CA point cloud, n = 126.** The candidate set is held at the same
deployed top-128 for every map, so only the assignment of basis states to candidates varies
(contract rule 7, the matched-space control).

## THE SIX MAPS AND THE THREE MEASUREMENTS

```
map            rho(Hamming, RMSD)   product-state ORACLE ceiling   total MI (bits)
score (deployed)     +0.0487                  2.0646                   0.0824
gray                 +0.0494                  2.0631                   0.0865
bisect               +0.1782                  2.0092                   0.1584
bisect_score         +0.1843                  1.9991                   0.1650
spectral             +0.1758                  2.0141                   0.1421
perm (control)       -0.0005                  2.0846                   0.0395

anchors: production top-75 uniform 3.0483 | uniform over the 128 3.0435
         ORACLE argmin over the 128 2.1458 | ORACLE convex over the 128 1.8008
```

**The matched zero-information control lands at rho = -0.0005, which is where it must land**, and
it is the reason the +0.0487 of the deployed map can be read as a real if tiny amount of
structural locality in the score order.

## F-C2c REFUTED AS REGISTERED, AND THE DIRECTION IS RIGHT

The registered pre-check needed `bisect` and `spectral` to raise rho by **>= 0.15** over `score`.
They raise it by **+0.130 and +0.127** -- a 3.7x multiple of a very small number, and short of the
bar. Reported as refuted, with the direction stated, because a near-miss on a near-tautological
pre-check is not evidence for the family.

## GRAY CODING IS A NO-OP, PROVED AND THEN MEASURED

Registered in advance as an analytic expectation: a `j`-bit prefix cell under the reflected Gray
code is the same block as under binary, up to reflection, so Gray coding can change only Hamming
geometry. **Asserted in code on 126/126 targets: the j-bit prefix PARTITION is identical for
j = 1..6.** Measured consequence: rho +0.0494 against +0.0487, product-state ceiling **-0.0014 A at
0.03x MDE**, partial-measurement cells identical to five decimals. **Gray-code candidate encoding
is closed.**

## F-C2a REFUTED -- the ORACLE ceiling moves, by less than the registered bar

```
map vs score, ORACLE product-state ceiling    effect    xMDE   fold CI            folds   W/L
bisect_score                                  -0.0655  -1.34   [-0.0818,-0.0480]   5/5   82W/44L
bisect                                        -0.0553  -1.11   [-0.0660,-0.0465]   5/5   80W/46L
spectral                                      -0.0504  -0.84   [-0.0591,-0.0365]   5/5   73W/53L
gray                                          -0.0014  -0.03   [-0.0262,+0.0313]   3/5   60W/66L
perm                                          +0.0200  +0.34   [-0.0172,+0.0535]   4/5   61W/65L
```

Real, above its own MDE, 5/5 folds -- and **-0.0655 A against a registered -0.10 bar, on an ORACLE
ceiling that deploys nothing.** REFUTED.

## THE ROW THAT IS WORTH MORE THAN THE FALSIFIER

The product-state ceilings are **2.0646 (score map) and 1.9991 (bisect_score)**, and the ORACLE
argmin over the same 128 is **2.1458**. **ORACLE / NOT DEPLOYABLE, CA point cloud:**

> **Seven qubits read out as a product-state weighted average beat the same seven qubits read out
> as a selection, by 0.081 A on the deployed index map and 0.147 A on the best one.**

And the product restriction costs **0.26 A** against the unrestricted convex optimum over the same
128 (1.8008). So the ansatz's product structure is a real but second-order loss, and the readout
class is the first-order one. This is the same ordering S31-L14 found from the classical side.

## F-C2b REFUTED, AND ITS MECHANISM IS THE C2 FINDING

Measuring `j` of the 7 bits localises a cell of `2^(7-j)` candidates, which is then averaged.
ORACLE-best cell against two native-free cell rules (best mean DIS; most typical):

```
map            j    ORACLE   by_score   by_consensus   mean cell
score          5    2.4216     3.2369      3.3126       3.2441
bisect         5    2.2712     3.3721      3.3541       3.4093
spectral       5    2.3022     3.2430      3.2840       3.3370
perm           5    2.5170     3.1753      3.1932       3.1912
```

**Structure-aware indexing raises the ORACLE value of a partial measurement and lowers its
deployable value, and the two move in opposite directions BY CONSTRUCTION.** Differentiating the
cells makes the best cell better (2.4216 -> 2.2712) and the chosen cell worse (3.2369 -> 3.3721),
because a coherent cell is a *concentrated* set and concentration is worth zero through an
averaging terminal -- S31-L14's identity, arriving at the index.

The best native-free triple over the whole 72-cell grid is **`perm`, j = 2, by_score -> 3.0323,
i.e. -0.0160 A against production**, which is REFUTED against the -0.10 bar. **The best index map
for a deployable partial readout is the random permutation.** Order-statistic price of the grid
(contract rule 11): per-target best-of-72 gain -0.6466, across-target null -1.2522, **194%
accounted**; split-half transfer **-0.1486, 23% of the ORACLE gain** -- so a per-target choice of
(map, j, rule) is mostly best-of-k and the transferable part is not the arm reported above.

## THE ANSWER TO "MAKE EACH QUBIT A MEANINGFUL DISTINCTION", IN BITS

Mutual information between each index bit and "this candidate is in the ORACLE-best decile of the
128" (ORACLE, and a diagnostic only):

```
bit             1       2       3       4       5       6       7     total
score      0.0343  0.0176  0.0082  0.0067  0.0064  0.0044  0.0047   0.0824
bisect     0.0658  0.0388  0.0236  0.0134  0.0091  0.0045  0.0032   0.1584
perm       0.0053  0.0048  0.0065  0.0069  0.0045  0.0066  0.0050   0.0395
```

**The best available map doubles a quantity that is 0.034 bits.** No single index bit carries even
0.07 bits about candidate quality, against the 7 bits the register nominally holds. The charter's
goal -- make each qubit correspond to a scientifically meaningful candidate distinction -- is
achievable in *direction* (bisect roughly doubles bits 1-3, and the permutation control sits flat
at 0.005) and is **two orders of magnitude short in size.**

## COMPARISONS MADE (contract rule 23)

6 maps x (1 rho + 1 ceiling + 6 j-levels x 4 cell statistics + 7 bit MIs) = **222 descriptive**;
5 fold-clustered contrasts for F-C2a; 1 grid maximum priced two ways. **Read as results: F-C2a,
F-C2b, F-C2c and the Gray no-op -- 4.** The 72-cell grid is never quoted without its
order-statistic price.

---

## S31-L16 -- **CROSS-LANE AUDIT, RUN AHEAD OF THE REQUEST BECAUSE CROSS-LANE SYNTHESIS IS WHERE THIS SPRINT'S ERRORS LIVE (4/4).** LANE C's ENTRY IS CLEAN AND SELF-GUARDED; THE RISK IS ENTIRELY IN RE-QUOTATION -- **A BROKEN CROSS-REFERENCE, A LEVER COMPARISON AT UNEQUAL INFORMATION COST, FIVE DIFFERENT "WORST 18" STRATA, AND A WIDENING WHOSE DEPLOYABLE SIGN IS PREDICTED *NEGATIVE* RATHER THAN NEUTRAL** (2026-09-21 00:24, L)

The coordinator reports his cross-lane syntheses are **4 for 4** on carrying the sprint's errors
(this lane's two, plus lane C's falsification of his R1 quantifier, lane B's correction to his
sweep reading, and lane A's catch of his wrong-operator deficit). **That is a measurement, and the
implication is that cross-lane claims need the same adversary single-lane claims get.** Serving the
standing request now rather than at report time, because a synthesis is cheapest to kill before it
is written down. Working: `s31/LIT_L.md` §L3.2.

### FIRST, THE PART THAT IS NOT A FAULT

**Lane C's S31-L5 is among the best-disciplined entries in the sprint** and I want that on the
record, because an auditor who only ever finds faults is useless. It labels **every number ORACLE
AND NOT DEPLOYABLE** in its own heading block; it found its circularity was a **literal identity**
(`defn18 INTERSECT FAIL18 = 18 of 18`) and withdrew the published -1.9004 as a near-tautology; it
carries a **level control** showing only -0.196 A of the clean tail's -1.181 A survives the
`best1(128)` level; it **self-caught a sign error** before any number left the lane; and its
"WHAT THIS DOES AND DOES NOT LICENCE" section says in its own words *"It does not say two more
qubits are worth 1.18 A... Widening only pays if a selector exists, and recognition is closed."*
**Nothing below is a criticism of the lane's own work.** All four items are re-quotation hazards.

### 1. A BROKEN CROSS-REFERENCE -- `ledger-numbers-collide-under-parallel-lanes`, FOURTH INSTANCE

S31-L5 says the result *"composes with the set-matched readout ladder (**S31-L6, next**)."*
**S31-L6 is lane D's projection-conditioning defect entry.** I checked: it contains neither
`1.8538`, nor `2.1435`, nor a readout ladder. Lane C reserved the next number in prose and lane D
took it first -- the same read-modify-write with no lock that collided three times in S30. A reader
following the pointer lands on an unrelated entry.

**The numbers are also not S31 lane work.** `1.8538` and `2.1435` originate in **S29**
(`s29/LEDGER.md:4170`), so the sentence re-quotes a two-sprint-old ORACLE ladder in a voice that
reads as this sprint's finding. Fix the pointer to the S29 origin, not to any S31 number.

### 2. THE LEVER COMPARISON IS MADE AT UNEQUAL INFORMATION COST -- AND S30 ALREADY NAMED THIS ERROR

S31-L5 concludes: *"the convex readout over the same 128 reaches 1.8538 A ... where the argmin over
those 128 reaches 2.1435, so the register width and the readout class are separable levers and
**the second one is larger** at the deployed width."*

The two ceilings do not cost the same:

```
best1_top128 = 2.1435   ORACLE argmin      -- costs 7 BITS (naming one of 128)
hull_top128  = 1.8538   ORACLE CONVEX weights over 128 -- costs 128 REAL NUMBERS, unbounded bits
```

S29 labels that rung **`[EXPRESSIVENESS]`** with **free weights** (`s29/LEDGER.md:4249`), and lane
F's own S31 line states it correctly: *"hull_top128 -- ORACLE CONVEX weights over the 128."* The
0.290 A is bought with unbounded oracle real numbers, so **"larger lever" is not established --
only "larger ceiling at unpriced cost."**

**S30 §9.4 already identified exactly this comparison as the project's characteristic error**: 2
members with ORACLE *weights* reach 1.4315 A against 75 with ORACLE *membership* at 2.3055, while
choosing 2 of 500 costs ~17.9 bits against 7 for the top-128 argmin -- from which S30 drew
*"solving the set problem better is worth nothing; what is scarce is the information needed to
SPECIFY a good set."* **Any lever comparison in this project must carry its bit cost in the same
sentence as its Angstroms.**

### 3. THE WIDENING'S DEPLOYABLE SIGN IS PREDICTED **NEGATIVE**, NOT NEUTRAL

Lane C's caveat is *"widening only pays if a selector exists."* Three facts compose to something
sharper, and no single lane holds all three:

1. **Lane C's own rank diagnostic is the strongest anti-deployment evidence in its entry, and is
   not flagged as such.** FAIL18's ORACLE-best member sits below rank 128 with frequency **1.000**
   against an uninformative null of **0.744** -- the score is **worse than uninformative** on
   exactly the targets widening is meant to help.
2. **My S31-L4 / §L1.4 gradient law**: widening 128 -> 512 costs **4x in gradient variance**
   (`Var ~ 16/D = 16*2^(-n)`; n=7 -> 0.125, n=9 -> 0.031). Widening *degrades the very selector it
   would require.*
3. **`operator-consumes-set-mean`** (`d_out = 1.16*d_set_mean + 0.04*d_set_best`, R2 0.89):
   admitting 384 more candidates to a set whose selector cannot order them **raises the set mean**,
   which the terminal operator consumes at coefficient 1.16.

> **Prediction, falsifiable and cheap: a DEPLOYABLE 128 -> 512 arm should come out WORSE than
> production, not merely flat.** That is a stronger and more useful statement than "it does not
> pay," and it should be measured before anyone widens anything. **Caveat that cuts the other way:
> my gradient law is measured at n=7-9 where we are NOT yet gradient-limited (§L1.4), so item 2 is
> a bound on the direction, not a quantified Angstrom cost.**

### 4. FIVE DIFFERENT "WORST 18" STRATA ARE NOW IN PLAY, WITH MATERIALLY DIFFERENT NUMBERS

This is the quotation hazard most likely to fire, because all five are called *"the worst 18"* in
prose somewhere:

```
stratum                          defined by                          a headline number
FAIL18                           the filter's own recall             best1_500 = 2.2842
defn18                           top-18 by best1(75)-best1(500)      = FAIL18, 18 of 18 (identity)
worst18_poolmean                 pool mean                           best1_500 = 2.2286  <- lane C's headline
worst18_bestpool                 best pool member                    best1_500 = 3.0930
"the genuinely worst 18" (S31-L0) charter's arithmetic                ORACLE best pool member 2.5298
worst 18 by production RMSD (L)  production endpoint                 mean RMSD 6.0291 (n=111 matched)
```

`worst18_bestpool` at **3.0930** and `FAIL18` at **2.2842** differ by **0.81 A** on the same-named
quantity. **Every "worst 18" number must name its stratum in the same sentence**, exactly as ORACLE
labels must (S31-L0's own rule). I recommend the report carry this table once, as a key.

### WHAT I AM *NOT* CLAIMING

I have not re-run lane C's arms and am not disputing a single number in S31-L5. Items 1 and 4 are
bookkeeping defects; item 2 is a framing defect with a named precedent; item 3 is a **prediction**,
not a measurement, and is labelled as one.

---

## S31-L17 -- **THE CVaR-VQE's HAMILTONIAN IS A TARGET-INDEPENDENT CONSTANT. THE QUANTUM STAGE CARRIES ZERO TARGET-SPECIFIC INFORMATION** -- STRONGER THAN T1 AND R1 TOGETHER (2026-09-21 00:36, lane P, verified independently by the coordinator)

Found by lane P at n = 36 of 126. **I verified it myself before recording it**, because it is the
sprint's headline and because every cross-lane claim I have made this sprint has been wrong.

### The fact

`core/pipeline.py` ~838 builds the Hamiltonian diagonal as `E = _zrank(pool["sc"][o])`, and
`_zrank` (`core/pipeline.py:788-792`) is `rankdata` followed by standardisation. **The ranks of any
128 distinct values are 1..128**, so the standardised vector is a **constant**. Verified on three
independent random score vectors:

```
trial 0  first5: [-1.718572 -1.691507 -1.664443 -1.637379 -1.610315]   last3: [1.664443 1.691507 1.718572]
trial 1  first5: [-1.718572 -1.691507 -1.664443 -1.637379 -1.610315]   last3: [1.664443 1.691507 1.718572]
trial 2  first5: [-1.718572 -1.691507 -1.664443 -1.637379 -1.610315]   last3: [1.664443 1.691507 1.718572]
```

**Identical to six decimals.** Up to the pool's tie structure, `E` is **the same vector on every
target in the benchmark.**

### The consequence

`p*` is a closed-form function of `(E, alpha, T)` alone (S31-L4), and `run_cvar_vqe` is seeded at 0.
Therefore **both the closed-form optimum and the circuit's output are ONE FIXED WEIGHTING CURVE PER
alpha**, identical across targets. Lane P's measurement:

```
H(p*) = 4.9136 bits at alpha = 1     sd across 126 targets = 1.4e-4
H(p*) = 6.6392 bits at alpha = 0.25  sd across 126 targets = 4.5e-3
ESS(p*) = 22.19 / 128  and  58.32 / 128
```

> **The quantum stage carries ZERO target-specific information.** The target enters the answer only
> through the readout's own `P` and `W` -- never through the objective, the Hamiltonian, the CVaR,
> or the state.

**This is strictly stronger than both of the sprint's earlier capacity theorems.** T1 said the state
specifies **one integer**. R1 said the selection readout's alphabet is **6.886 bits**. This says the
state specifies **nothing at all**: the stage answers *"what fixed weight should rank k receive?"*,
which is a **128-number global hyperparameter, not a per-target computation.** The rank -> candidate
*mapping* is target-specific; the *weight on each rank* is not.

### Why it explains every other negative in the sprint

- **Why the objective does not point at good solutions** (charter §5E): it is not a function of the
  target. There is nothing for it to point *at*.
- **Why alpha does nothing beyond reshaping a fixed curve**, and why alpha = 1 on three of five
  pinned folds makes CVaR inactive by construction (lane A).
- **Why the circuit's 0.930-bit optimisation gap cancels in the mean** (lane A: readouts differ at
  0.088x MDE while differing per-target by 0.102 A) -- both arms are the same fixed curve, perturbed.
- **Why S31-L2's convex optimum over that family converged to the uniform average** (3.0522 from
  production): the family's whole design space had already been searched, and the answer is uniform.

### The deployed temperature, which corrects S31-L4

`core/pipeline.py:113`: `VQE_LFO = {0:(1.0,0.3), 1:(0.25,0.3), 2:(0.25,0.3), 3:(1.0,0.3),
4:(1.0,0.3)}` -- **T = 0.3 on every fold.** Lane L's 12 verification cells were run at T = 0.1 and
0.05, so **none was at the deployed temperature**, and its statement that the closed form is always
the more entropic distribution **reverses at the deployed T**: at alpha = 1 (folds 0/3/4, three of
five) `H* = 4.91` against `H_vqe = 5.67` -- the circuit is *more* entropic than the optimum. Lane L
has closed; recorded here against its entry with the original standing. Lane P's own registered
prediction rests on the reversed premise and it is reporting that prediction as **falsified** rather
than dropping it.

### A control nobody asked for, which calibrates every built-chain comparison this project makes

Lane P found that arm A uses `I.coordinate_average` while arms D/E/F use
`core.pipeline.average_weighted` -- **the same operator in two implementations, agreeing only to
~1e-14**, fed through a projection with **~1e13 amplification** (S31-L6). So `D - A` and `F - A`
contain an unknown amount of pure implementation noise. It is running `A2 = average_weighted(uniform,
top-75)` so that `A2 - A` **measures that noise on the built chain** and `F - A2` is the clean
effect. **Without it the sprint would have published a contaminated number on a row I promoted.**
