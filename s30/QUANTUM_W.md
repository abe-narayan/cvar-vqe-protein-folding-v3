# S30 — THE QUANTUM SPINE, ITEMS 12–19 (lane W)

Answers to charter §19 items 12–19 (`s30/BRIEF.md:563-570`). Investigation only: nothing here was
computed for this document except two reads of artefacts already on disk (the `circ_*` rung means
in §12, recomputed from `s30/results/s30_D_ladder_structs/*.npz`, n = 126). No lane result is
restated without its file, line or artefact.

---

## 0. THE HONEST HEADLINE

**No CVaR-VQE was trained in Sprint 30.** No new Hamiltonian was built, no ansatz was designed, no
trainability, gradient-variance or DLA measurement was taken. S30 was a measurement-and-theory
sprint, and its lanes concluded — with theorems in two cases and pre-build measurements in two more
— that the binding constraint sits **upstream of the quantum stage**, in what the system can
*recognise*, not in what it can *express*, *optimise* or *encode*.

Five of the eight items therefore read **NOT APPLICABLE THIS SPRINT**. That is not eight blanks. It
is the sprint's position, and three of the five were closed by S30 work that makes them *stay*
closed:

- **T1** (`s30/THEORY.md:43-80`) proves the CVaR tail is *always* a prefix — of the order induced by
  `∇V` at the optimum — so item 12's question has a closed-form answer rather than an empirical one,
  and the charter's own L5 lead was mis-posed.
- **T1b + the quadric measurement** (`s30/THEORY.md:112-275`, S30-L12, S30-L15) priced and closed the
  two Hamiltonian classes that *would* have made item 14 answerable — **before either was built.**
  This is the first time in the project's history the cheap pre-check fired before the spend rather
  than after (S29's first non-diagonal Hamiltonian failed for exactly this reason, discovered
  afterwards: `s30/THEORY.md:161-166`).
- **S30-L19** (`s30/LEDGER.md:2680`) closed single-structure nativeness *preference* on 43 channels
  with a matched-kind ladder, which is the statement that makes the whole quantum stage not the
  binding constraint: *"If nativeness cannot be preferred from a single structure's geometry, then
  no Hamiltonian, no encoding, no ansatz and no amount of search can be pointed at it — the
  objective would be optimising toward a target it cannot see."* (`s30/STATE.md:722-726`).

**Three things the checking turned up that a reasonable prior would have got wrong**, stated here
because they are the part that cannot be guessed:

1. **"No VQE" and "no quantum circuit" are different claims, and only the first is true.** S30 *did*
   execute a parameterised 9-qubit statevector circuit on all 126 targets, with exact
   parameter-shift Jacobians and 300 Adam iterations per target — see §12.1. It was not a VQE (its
   objective is ORACLE RMSD to the native, not CVaR of a Hamiltonian) and it was a regeneration
   check, asserted to 1e-6 against S28's stored value. But the circuit ran, the outputs are on disk
   in `s30/results/`, and a report sentence saying "no quantum compute was spent this sprint" would
   be false as written.
2. **The most damaging single number about the quantum stage in the whole record is an S30 cache
   entry.** `circ_opt` — the same circuit trained on the **deployed native-free objective** — sits at
   **3.4330 Å** on the built chain against production's **3.2071 Å**. Same circuit, ORACLE objective:
   **0.2516 Å**. S30 recomputed and re-cached both on 126 targets and used the pair as the meter's
   rungs (§12.1, §16).
3. **`core/pipeline.py:821` still asserts a number withdrawn in S25.** The `quantum_stage` docstring
   states *"the CVaR tail is worth +0.113 Å by preventing the collapse, and that is the component's
   measured role."* S25-L5 (`s25/LEDGER.md:233-244`) withdrew it: *"was never a measured effect …
   labelled off MARGINAL MEANS … by this project's own fixed rule that is a NULL."* The replacement
   figure is **−0.1405 Å at 0.68× MDE, not a result** (`s25/LEDGER.md:280`). Anyone answering item 12
   from the source file rather than the ledger gets a withdrawn positive. Flagged, not touched
   (`core/` is read-only for this lane). This is a **sixth** instance of the
   `findings-prose-is-not-evidence-of-code` pattern, and the first one located in shipped code
   rather than in a report.

### How this was verified, so the negative is auditable

| check | result |
|---|---|
| `grep -rniE "cvar\|vqe\|ansatz\|qiskit\|pennylane\|qubit\|hamiltonian\|statevector\|SparsePauli\|quantum"` over all 24 `s30/s30_*.py` (16,557 lines) | **2 hits, both comments**: `s30_P_prior.py:6` (prose), `s30_Q_quadric.py:3` ("measured BEFORE any Hamiltonian is built") |
| `grep -rn "quantum=True"` over the repository | 8 hits, **all in `s26/`** (S26 lane E/Q), **none in `s30/`** |
| `grep -rniE "\bDLA\b\|ansatz\|trainab\|barren\|gradient variance\|parameter.shift"` over all of `s30/` excluding `BRIEF.md` | **3 hits**, all prose in `STATE.md`/`THEORY.md` saying the ansatz is *not* the barrier |
| scan of all 30 `s30/results/*.json` for keys `cvar`/`theta`/`p_theta`/`vqe`/`qubits`/`entropy_bits`/`ansatz` | **zero** |
| lane assignment (`s30/LEDGER.md:85-95`, `s30/STATE.md:1194-1203`) | 7 lanes: R, F, Q, X, T, D, L. **No VQE lane was staffed**; the 8th slot was held open and never used for one |
| 28 ledger entries S30-L0 … S30-L27 | none reports a VQE run |
| `s26/jobs_done/s30*` | 12 job records, none a VQE job |

---

## ITEM 12 — What did CVaR-VQE specifically contribute?

**NOT APPLICABLE AS AN ENDPOINT CONTRIBUTION THIS SPRINT, because no CVaR-VQE was run and the
production anchor does not pass through the quantum stage.** What S30 contributed instead is a
*theorem* about what CVaR can contribute, and a *measurement* of what the circuit does when pointed
at the deployed objective.

### 12.0 The production anchor does not pass through the quantum stage — from the code

- `core/pipeline.py:241` — `PROD = Config()`, i.e. all defaults.
- `core/pipeline.py:179` — `quantum: bool = False   #: CVaR-VQE over the discrete hypothesis set`.
- `core/pipeline.py:173-178`, the comment governing that default, verbatim: *"the four mandated
  components. `s9/final.py` is retrieve -> filter -> synthesise -> AMBER: AMBER participates, Legacy
  is only imported for scoring, and **VQE/CVaR do not participate at all**."*
- `core/pipeline.py:1076` — `if cfg.quantum:` is the only gate on `quantum_stage`, and
  `core/pipeline.py:758` — the top-128 widening `(1 << cfg.vqe_qubits) if cfg.quantum else 0` is the
  only other place the flag is read.

So the sprint's 3.2105 Å anchor and every S30 comparison against it are computed with the quantum
stage off. **Stated as the file:line the charter asked for: `core/pipeline.py:179` plus
`core/pipeline.py:241`.**

### 12.1 What DID run: a circuit, on 126 targets, not a VQE

`s30/s30_D_meter.py build-cache --chain` (four shards; job records `s26/jobs_done/s30D_cache_0..3.json`)
calls `M.build_rungs` at `s30/s30_D_meter.py:145`, which is `s29/s29_D_cost_audit.py:250`. That
function, at lines 255 and 265-266:

```
from core import quantum as Q
circ = Q.StatevectorCircuit(A.N_QUBITS, A.LAYERS)          # 9 qubits, depth 3
oc = A.oracle_circuit_ceiling(circ, frame, cand.nat_ca, starts=1)      # ORACLE
```

`s27/s28_A_amp.py:422-449 oracle_circuit_ceiling` runs Adam for `ORACLE_ITERS = 300`
(`s27/s28_A_amp.py:60`) on the exact parameter-shift Jacobian
(`circ.states_batch(th[None,:] + math.pi * np.eye(P))`, line 439), minimising **point-cloud RMSD to
the native**. The result is asserted against S28's stored `per_start[0]` to `dev_s0 < 1e-6`
(`s29/s29_D_cost_audit.py:271`).

Artefact: **`s30/results/s30_D_ladder_structs/` — 126 `.npz` files**, keys including `circ_s0`,
`circ_best`, `circ_opt`, with CA and built-chain coordinates and ORACLE RMSD labels.

Recomputed from those 126 files for this document:

```
rung          what it is                                            CA cloud        built chain
NATIVE        the rebuilt native (projection floor)                 0.0000          0.0837  (SE 0.0079)
circ_best     9q circuit, ORACLE objective, best of 5 starts        0.2884 (0.0173) 0.2516  (0.0152)
circ_s0       9q circuit, ORACLE objective, 1 start (regenerated
              IN S30 on all 126)                                    0.3854 (0.0221) 0.3175  (0.0188)
sub0          least squares in a random subspace, ORACLE            0.6199 (0.0426) 0.4992  (0.0358)
PROD          the deployed uniform average of the DIS top-75        3.0483 (0.1467) 3.2071  (0.1541)
circ_opt      THE SAME CIRCUIT, trained on the DEPLOYED
              NATIVE-FREE OBJECTIVE (lam=1, 80 iters)               3.3850 (0.1420) 3.4330  (0.1454)
```

`circ_opt`'s provenance: `s29/s29_D_cost_audit.py:108`, `CIRC_OPT_KEY = "circ_l1_i80"` —
*"lane A's native-free recognition optimum, lam = 1, 80 iters"*, loaded at line 260 from
`s27/results/s28_A_structs/<pdb>_recog.npz`.

**This three-line contrast is the honest answer to item 12.** The circuit's amplitude readout can
represent a 0.25 Å structure and an optimiser can find it *when the objective is the native*. Point
the identical circuit at the deployed native-free objective and it lands **0.2259 Å worse than the
classical average it was meant to improve** (3.4330 vs 3.2071, built chain, n = 126). The gap
between 0.2516 and 3.4330 is the objective, not the circuit.

### 12.2 The theorem that replaces the empirical question — T1

`s30/THEORY.md:43-56` (lane T, S30-L9, `s30/LEDGER.md:840`):

> **T1.** Let `V` be differentiable on `Λ(p) = {λ : 0 ≤ λ ≤ p, 1'λ = α}`. At every KKT point `λ*` of
> `min_{Λ(p)} V` there is a scalar `μ` with `∇V(λ*)_x < μ ⟹ λ*_x = p_x` and `∇V(λ*)_x > μ ⟹ λ*_x = 0`.
> Hence `λ*` is a prefix of the order induced by `∇V(λ*)` — **for every `V`**.

`Λ(p)` is the exact feasible set of `core/quantum.py:290 cvar_from_probs`, whose lines 306-308
compute one of its vertices by a cumulative scan. The deployed case is `V(λ) = ⟨E, λ⟩`, so
`∇V = E`, **constant in λ** — *"That — and nothing about CVaR, quantiles, order statistics or the
circuit — is why one classical sort reproduces it"* (`s30/THEORY.md:60-62`).

**CVaR's residual contribution, stated exactly once and precisely** (`s30/THEORY.md:178-185`):

> Under an endogenous order the cut is a pair `(g, τ)`. `g` comes from `∇f`; **`α` supplies `τ`.**
> That is CVaR's whole remaining function — it turns a direction into a set with a *budget* rather
> than a threshold, which is what keeps the map continuous in `p`. **It is not doing selection.**

### 12.3 And it corrects the charter's own lead

`s30/BRIEF.md:111` lists as an S29 finding: *"a free subset-search objective produces non-prefix
optima, meaning there is a genuine combinatorial selection problem available that the deployed VQE
has not yet solved"*, and L5 (`s30/BRIEF.md:304-315`) builds on it. T1 says the free optimum is
real and **the lifted CVaR readout cannot reach it**, because S29 §4.3's own well-posedness
recommendation (fix the order by a per-state scalar) *provably guarantees* the prefix. Lane T's
disposition table (`s30/THEORY.md:97-105`): the non-prefix optimum **stands**; *"the lifted CVaR
readout reaches it"* is **WITHDRAWN**; *"the first formulation whose classical counterpart genuinely
goes away"* is **WITHDRAWN as stated**.

Independently, lane L (S30-L8, `s30/LEDGER.md:764-838`) tested its own CVaR hypothesis and refuted
it: finite-shot CVaR estimator bias at the deployed cell (9 qubits, α = 0.18, 2048 shots, tail
k = 369; `s30/lit/s30_L_cvar_bias.py`) is **+0.0001 to +0.0019 with no trend in support**, and in
the direction *opposite* to the guess — *"Finite-shot CVaR estimator bias is not a mechanism"*
(`s30/LEDGER.md:825`).

---

## ITEM 13 — What Hamiltonian was ultimately used?

**NOT APPLICABLE THIS SPRINT: none. No Hamiltonian was built, and the charter withdrew the
obligation to build one.**

`s30/BRIEF.md:481-495`, §14 "THE QUANTUM SPINE — NO PRESCRIBED ARCHITECTURE", verbatim:

> Previous sprints supplied a conceptual pipeline template. **That template is withdrawn.** … There
> is no diagram in this prompt to work toward. Design the architecture you can defend. The only
> criterion for whether a stage belongs: **Does this stage help CVaR-VQE solve a structurally
> meaningful optimization problem and reduce final built-chain RMSD?**

No stage was proposed that met that criterion, so nothing was built. For completeness, the
Hamiltonian that *would* be used if `quantum=True` is unchanged from S25/S26:

```
core/pipeline.py:838    E = _zrank(np.asarray(pool["sc"], float)[o])
core/pipeline.py:788    _zrank -> standardised rank (scipy rankdata, centred, unit sd)
core/pipeline.py:839    alpha, T = VQE_LFO[int(fold) % len(VQE_LFO)]    # leave-fold-out (alpha, T)
core/pipeline.py:113    VQE_LFO = {0:(1.0,0.3), 1:(0.25,0.3), 2:(0.25,0.3), 3:(1.0,0.3), 4:(1.0,0.3)}
```

`H = diag(zrank(top-2^n distogram scores))`, objective `F = CVaR_α(E; p_θ) − T·H(p_θ)`, 7 qubits
(`core/pipeline.py:181`), 3 RY/CNOT layers, 50 Adam iterations, exact statevector. Deployed **off**.

**Two Hamiltonian classes were considered and both were closed before any was built** — which is
S30's actual contribution to this item:

| class | what it is | closed by | verdict |
|---|---|---|---|
| **endogenous-order linear** (T1b) | `V(λ) = f(R_λ)` ⟹ `∇V_x = (1/α)⟨∇f(R_λ), W_x⟩`; tail = halfspace cut, VC dim `d+1` | `s30/THEORY.md:112-166`; rank collapse measured on 126 real pools by T (`s30_T_bits.py`, `s30_T_spec.py`) and independently by Q through the **deployed** path | Reachable, **not exploitable**. `r_stable = 1.859` (pair-distance) / 3.404 (coordinate); class carries 24–80 bits against 300.6 of free 75-subset choice. Row 2 of the trichotomy is classically polynomial (Frank–Wolfe), **so it is not a quantum opening** (`s30/THEORY.md:167-176`) |
| **second-moment / quadric** | `V = f(R_λ) + μ·h(Σ_λ)` ⟹ `∇V_x` quadratic in `W_x`; cut is a quadric, VC dim 34 → 595 | S30-L12 (Q, `s30/LEDGER.md:1531`) **and** S30-L15 (T, `s30/LEDGER.md:2012`), independently, two samplers | **CLOSED.** `+0.2059 Å` worse at the shipped m = 75 (1.81× MDE, 37W/89L) and `+0.086 Å` worse than the halfspace class it was meant to escape into (2.95× MDE, 5/5 folds). The dispersion rule it points at measures **3.3585 Å against production's 3.0483** |

Lane Q's S30-L12 heading carries the sprint's methodological point on this item verbatim:
**"MEASURED BEFORE A HAMILTONIAN WAS BUILT, AS INSTRUCTED"**.

One retraction belongs here, because the class's *apparent* value was the sprint's largest
self-correction (`s30/THEORY.md:235-274`, §4.3): lane T first reported the halfspace class at
−0.9816 Å BETTER than the deployed prefix. Lane Q named the right null, the coordinator endorsed it,
lane T ran it on its own K = 5,000, and **196% of the gain is accounted for by an across-target
null**; the transferable rule lands **+0.472 Å worse than what ships** (1.88× MDE, 5/5 folds,
39W/87L). *"The halfspace class is reachable, not exploitable, and as a rule it is negative."*

---

## ITEM 14 — Why is it meaningfully different from the earlier diagonal rank-ladder?

**NOT APPLICABLE THIS SPRINT: there is no new Hamiltonian to be different, and S30 established
that the diagonal rank-ladder is not the thing that is wrong.**

### 14.1 The referent, made concrete

The "earlier diagonal rank-ladder" is `H = diag(zrank(sorted top-128 scores))`
(`core/pipeline.py:838`). Its defining property, measured in S25 and never overturned:

- **It is the same ladder on every target to 1.18% of its range** — `s25/results/q_gibbs.json ::
  results/spectrum_target_independence/{worst, E_range, worst_frac_of_range}` = `0.04063`, `3.43714`,
  `0.011821`; reproduced at 0.394% on 1S9Z and 9KAR in `s26/EXAMINATION.md:221`; also
  `ARCHITECTURE.md:176`, `s25/LEDGER.md:1165`.
- Consequence (S25-L17): *"the circuit trains toward one of two states regardless of target"*
  (`docs/REPORT_S26_SUMMARY.md:17`).
- Its tail is provably a subset of the classical prefix — **0 violations on 2,592 cells**
  (`s25/results/q_verify.json`), and T1 now explains *why* without reference to CVaR at all (§12.2).

### 14.2 What S30 established about it

The charter's item 14 presupposes that the rank-ladder's *diagonality* is the defect to escape. S30
says that presupposition is wrong in a specific, measured way:

1. **The off-diagonal escape was already closed in S28 and S30 explains the closure.** S28's
   `H = diag(E) − J·A` (`s27/REPORT_S28.md:104-121`) changed the emitted structure by nothing
   measurable at every `J`, graph and seed. S29-L11 gave the law
   (`Var[∂F/∂θ] ≈ r_stable(A)/D²`, `s29/THEORY.md:394`) and S30's independent rank-collapse
   measurement on **126 real pools** (`s30/THEORY.md:126-166`) confirms the mechanism generalises:
   every similarity structure this pool admits is near-rank-one in the feature space that matters.
2. **The scoped correction S30 added, which must travel with the number.** At lane Q's request and
   adopted in full (`s30/THEORY.md:140-152`): *"the collapse is a property of the PAIR-DISTANCE
   feature space, not of the pool."* In **coordinate** space `r_stable` is 3.404/3.619 with
   `k90 = 11.2` — above lane T's own 2.0 threshold. So "stable rank 1.86, the lift is a relabelled
   one-dimensional sort" is correct for a distance-map lift and an **overstatement** for a
   coordinate-space one. The verdict is unchanged (an 11-parameter continuous family is still a
   ~220-bit collapse below the set choice) but it now rests on direct measurement, not the rank
   threshold. *"The number must never appear without its feature space in the same sentence; it was
   already circulating without one."*
3. **The barrier is not the Hamiltonian's structure at all.** `s30/THEORY.md:450-453`, listing what
   the residual 0.98 Å gap is *not*: *"It is not the encoding (§7), not the objective's functional
   form (S29 Theorem 2), not the Hamiltonian's off-diagonal structure (S29 §3), not the ansatz, not
   the optimiser, not the pool, and not the search (§4.2, and S30-L8). It is under four bits of
   direction per target that nothing native-free supplies."*

**So the answer to item 14 is: a new Hamiltonian would be meaningfully different only if it changed
the order the tail is a prefix of, and S30 measured the two candidate ways of doing that and found
both reachable-but-not-transferable. The rank-ladder is not the barrier; the thing that would
replace it is not available.**

---

## ITEM 15 — What did the quantum state represent?

**Partly applicable.** No state was prepared for a protein objective this sprint. But S30 answered,
in closed form and with a measurement, *what it can represent and what that is worth* — which is
the durable half of the question.

### 15.1 What it would represent, if run

`p_θ` over `2^n` basis states = the `2^n` filtered candidates, index `i` at bit-pattern `i`
(`core/pipeline.py:836-838`; the labelling is the identity permutation —
`s22/qcand_lib.py :: Encoding` exposes a free relabelling that the deployment does not use,
`s29/THEORY.md:1105-1112`). The state enters the answer two ways: as a `p_θ`-weighted consensus
medoid (`core/pipeline.py:847`) and, in the ablation arm, as a `p_θ`-weighted coordinate average
(`core/pipeline.py:855 average_weighted`).

### 15.2 T1 says what the state can determine: one integer

Under the deployed (exogenous-order) formulation, T1 collapses the state's entire structural
influence to the tail's cardinality. `s30/THEORY.md:167-176`, the trichotomy:

| order | tail | classical counterpart | **endpoint channels** |
|---|---|---|---|
| **exogenous** (deployed) | `argsort(E)[:m]` | one sort | **1** (the integer `m`) |
| endogenous, `f` convex | halfspace cut | Frank–Wolfe, `O(1/ε)` sorts | `d_eff` ≈ 2–6 |
| endogenous, `f` non-convex | halfspace cut, several fixed points | greedy + local search over ≤ 175 bits | `d_eff`, genuinely combinatorial |

Row 1 is *verbatim* `s29/DATAPATH.md` stage 9: *"the whole quantum stage reduces, for the structure,
to choosing `m`."*

### 15.3 The register is a codebook index, not a channel — S30-L14

Lane T's bit accounting (`s30/LEDGER.md:1893`, `s30/THEORY.md:275-317`), in a currency defined in
the prereg before it was measured:

| stage | bits of CHOICE | bits DELIVERED |
|---|---|---|
| retrieval, 17,088 → 500 | 3,252 | **+0.69** [SE 0.18, median +0.44] |
| energy, scores → ranks | 3,767 | — |
| top-128 prefix | 300.6 | ORACLE ceiling 0.066 Å; transferable part **0** |
| **tail readout, which of 128** | **7** | **+0.036** [MDE 0.385] |

> **Under one delivered bit per target, against 7,019 bits of choice consumed.**

And the inversion that matters for item 15: **7 index bits realise 36.6 bits of displacement
information (5.2×)**, because the basis state names a real 33-dimensional structure rather than
carrying a payload. *"No bits went missing — the readout's 7 are worth five times their face value
and the system cannot supply one."*

The **value-of-a-bit law** prices the state's remaining headroom (`s30/THEORY.md:318-341`, P1a/P1b/P1c
all HELD): `D(R) = a + c·2^(-R/γ)` with `a = 1.3312 Å, c = 2.7859, γ = 3.1636, R² = 0.9983`, so
`−dD/dR = 0.2191·(D − 1.331)` = **0.1317 Å per bit at R = 7**. The allocation theorem
(`s30/THEORY.md:343-374`) then settles what the state *should* index: candidate identity **wins by
3×** over subset cardinality (0.132 vs 0.044 Å/bit); basin/mode indexing saturates at ~1.6 bits;
and **torsion encodings are arithmetically infeasible** — 25.9 bits at 2 basins/residue, 51.8 at 4,
against the deployed 7 (`s30/LEDGER.md:1985-1990`).

### 15.4 The corollary lane X added about a quantum encoding specifically

`s30/LEDGER.md:1093-1100` (lane X, S30-L10/L20; `s30/s30_X_FINDINGS.md:129-135`):

> A `2**q` register over structural variables is by construction a **wide** space … Width is the
> encoding's selling point — and at ρ ≈ 0 width buys ceiling and costs mean. **The property a
> quantum encoding is chosen FOR is the property this condition penalises.** … **Generation is
> closed jointly with the readout.**

The coefficients are a property of the **readout**, not of physics: under argmin they would be
(0, 1) and the ceiling *would* be the endpoint. Lane Q's S30-L11 then shows the argmin already wins
on price (§16.3), which is what leaves the whole wide-encoding direction with no consumer.

---

## ITEM 16 — Could a classical control reproduce the result?

**THE HONEST FORM OF THE QUESTION THIS SPRINT IS: IS THERE ANY POSITIVE FOR A CLASSICAL CONTROL TO
REPRODUCE? THE ANSWER IS NO.** S30 produced no quantum positive, and no endpoint improvement of any
kind that runs through the quantum stage. The sprint's results are overwhelmingly negative or
theoretical; §3 of the report (`s30/REPORT_S30.md:55-73`) is a table of eleven *closures*, and §0,
§11, §12 and §13 are still `[PENDING]`.

Reading item 16 as written would invite a sentence of the form "yes, a classical control reproduces
it", which would imply a quantum result existed to be reproduced. **It did not.** The defensible
statement is the stronger one below.

### 16.1 Where a classical control is decisive, it is decisive against the quantum stage, not for it

| claim | classical control | result |
|---|---|---|
| The deployed CVaR tail is a quantum selection | **one classical sort** | Exact. T1 (`s30/THEORY.md:43-64`): `∇V = E` is constant in `λ`, so `cvar_from_probs`'s cumulative scan *is* greedy fractional knapsack and the rearrangement inequality is the whole content. Empirically 0 violations on 2,592 cells (`s25/results/q_verify.json`); S28 re-verified `gate_set_equality` at 4,914/4,914 cells to **1.1e-13 Å** at every `J` (`s27/REPORT_S28.md:116`) |
| The endogenous-order lift is a quantum opening | **Frank–Wolfe on the simplex** | `s30/THEORY.md:167-176` row 2 — *"classically polynomial, so it is not a quantum opening"*. Lane L (S30-L8, `s30/LEDGER.md:764`) independently: `V(S) = f(mean_S W)` **factors through the centroid**, so every submodularity guarantee is void (they require monotone; S29-L25 says `V` is not) and Frank–Wolfe is the correct classical counterpart |
| A sparse weighted quantum readout beats the classical argmin | **plain `argmin` over top-2^B** | S30-L11 (`s30/LEDGER.md:1206`): argmin wins **at every budget from 3 to 9 bits**, by +0.162 to +0.727 Å. On the built chain, 2-of-75 with **free continuous weights** (11.4 bits + unbounded) is 2.1683 Å against top-128 argmin's **2.1435 Å at 7.0 bits**. *"S29's last unclosed ladder class is closed — not by ceiling, by PRICE"* |
| The circuit's amplitude readout beats the classical average | **the deployed uniform average (PROD)** | §12.1: `circ_opt` **3.4330 Å** vs PROD **3.2071 Å** on the built chain, n = 126. The circuit optimising the deployed objective is **worse than not running it** — the `concentration-is-wrong-when-discrimination-binds` pattern, exactly as recorded |

### 16.2 The one control that would matter cannot be run, because its subject does not exist

The charter's bar (`s30/BRIEF.md:540-542`) is: *"the CVaR-VQE to contribute information or
optimization behavior that a matched classical control does not reproduce."* No S30 arm reached the
stage where that control would be informative. The report should say so in exactly those terms
rather than leaving item 16 to imply a null test was run and passed.

### 16.3 One control discipline S30 did exercise, worth carrying

S30-L24 (`s30/LEDGER.md:3104`) is the sprint's cleanest control result and it is about the *meter*,
not the circuit: the meter's new 8-draw random-signed control caught that **S29's seed-0 draw was
the maximum of its own eight on the built chain**. The C2 clause-2 contrast fell +0.0635 → +0.0357,
0.92× → **0.56× MDE**, fold CI now spanning zero — below the sprint's "not a result" floor. On CA
the same single draw was fine. *"The draw noise is a property of the basis, not of the control."*
This is the project memory `draw-controls-need-their-own-distribution` earning its keep on first
run, and it is the reason the `circ_best`-vs-PROD preference numbers in §12.1 can be quoted at all.

---

## ITEM 17 — What happened to trainability?

**NOT APPLICABLE THIS SPRINT: no circuit was trained on a protein objective, so no trainability was
measured.** The standing position is S28/S29's and S30 neither tested nor moved it. Recorded here so
the report can cite it rather than leave a blank.

| finding | source |
|---|---|
| The deployed diagonal problem is the **trivial regime**. At α = 1 the state the optimiser heads for is a **basis state** — a product state that seven parameters represent exactly — so an ADAPT-grown circuit has nothing to train; grown circuits hold 1–3 distinct operators, the rest repeats of one single-qubit Y. *"A large gradient from a product circuit is not trainability; it is the trivial regime."* | S26 lane Q A4; `s26/results/q_var.json`; `s26/LEDGER.md` L35; `docs/REPORT_S26.md:332` |
| At α = 0.25, T = 0.3 the L2 pool grows 4–21 distinct 2-local strings and decays at −0.302/qubit, with the fixed ansatz's −0.243 inside the error of a 7-point slope | same |
| S28's `J = 3` ground state is **not representable** (max overlap 0.80–0.88 over 16 starts, saturated at 80 iterations), but a representable state has `F` below the ground state's on 10/12 targets and below every reached `F` on 24/24 cells: **the gap is OPTIMISATION by basin selection, not expressivity.** The trained state is **bimodal on all 126** (43/32 cells by seed fully coherent, the rest at the `J = 0` state, none between) | S28-L41, L43; `s27/results/s28_B_represent.json`; `s27/REPORT_S28.md:120` |
| A spread-spectrum kNN graph halves the gradient decay and the circuit then finds the coherent basin on 25–39% of (target, seed) cells — **and the emitted structure does not move on either class** | S28-L25, L29, L46, L47; `s27/REPORT_S28.md:124-130` |
| S29's summary: *"the deployed circuit family can express structures around 0.25 Å in ORACLE diagnostics; the optimizer can optimize its objective; **the current objective is the barrier**"* | `s30/BRIEF.md:91-95` |

S30's only addition is the *reason* the trainability question is no longer load-bearing:
`s30/THEORY.md:450-453` removes the ansatz and the optimiser from the list of things the residual
gap could be, and §12.1's `circ_opt` number is the direct demonstration — the optimiser succeeds at
its objective and lands 0.23 Å worse than the classical baseline.

**One caveat the report must carry with any trainability claim** (S29-L11, `s29/THEORY.md:1558`,
listed in S29's own withdrawal table): *"Centre the coupling matrix to fix trainability"* was the
coordinator's own premise, lane B was spawned on it, and lane T's derivation killed it — centring
makes the decay **worse** (−2.305 vs −1.830).

---

## ITEM 18 — What happened to gradient variance?

**NOT APPLICABLE THIS SPRINT: not measured, because nothing was differentiated against a protein
objective.** The standing law is S29's and S30 did not revisit it.

**The law** (`s29/THEORY.md:380-398`, S29-L11), for an observable of unit spectral norm on a
`D`-dimensional register:

```
Var[∂F/∂θ] ≈ tr(A_0²)/D² = r_stable(A)/D² ,     r_stable = ||A||_F² / ||A||_2²
```

Three checks with no free parameter (`s29/THEORY.md:415-431`): S28's Gaussian graph at
`r_stable = 1.036`, `n = 9`; the kNN graphs at 1.591/1.675; measured hop-only variance 4.06e-3 at
`n = 4` → 4.16e-6 at `n = 9` (`s27/REPORT_S28.md:112`).

**The design rule, and the wall it hits** (`s29/THEORY.md:427-437`):

> `Var_hop / Var_diag = r_stable / (D² Var_diag)`, so parity at `n = 9` needs `r_stable ~ 8000`,
> i.e. `r_stable ~ 16 D`. **No dense similarity kernel can do this, centered or not**; a `k`-regular
> sparse graph reaches `r_stable = M/k` (parity would need `k ~ 0.06`); and any Gram of structural
> features has `r_stable ≤ rank ≤ 3 N_res − 6 ≤ 42` on this instrument.

The escapes it permits: a **sum of local Pauli terms**, `r_stable(Σ_q X_q) = D/n`, giving
`Var ≈ 1/(nD)`, slope exactly −1 per qubit and `J* = 11.8` at `n = 9` (`s29/THEORY.md:756`); or a
smaller register. Neither was built.

**What S30 contributed to this item:** it measured the *input* to the law on 126 real pools rather
than on a constructed kernel — `r_stable = 1.859` (pair-distance) and `3.404` (coordinate),
`k90 = 5.61` / `11.52`, PC1 carrying 55.4% / 30.0% (`s30/THEORY.md:133-139`, verified independently
by lane Q through the deployed `s29_O_ladder.load_pool` + `s28_A_amp.Frame` path). **The check cost
76 seconds over data already on disk** (`s30/THEORY.md:163-166`), against S29's first non-diagonal
Hamiltonian, which failed for the same reason discovered *after* the attempt. That is the durable
output: the gradient-variance law's governing quantity is now a measured property of the real pools,
so the next off-diagonal proposal can be priced before it is written.

---

## ITEM 19 — What happened to DLA / ansatz structure?

**NOT APPLICABLE THIS SPRINT: no ansatz was designed, grown, or analysed. Zero DLA computations in
`s30/`.** The census is S29's, complete, and unchanged.

**The ansatz** (`s29/THEORY.md:1037-1039`): RY on every wire, CNOT chain plus ring closure, depth 3,
real amplitudes; `P = nL = 27` parameters at `n = 9` (`core/quantum.py:818-903`).

**The DLA census** (`s29/s29_T_reach.py --dla` → `s29/results/s29_T_reach.json :: dla`, using
`s26/q_dla.py`'s exact Pauli-set closure, cross-checked against an independent dense
nested-commutator closure with SVD rank at every `n ≤ 5`, 12/12 agree):

- **The obstruction occurs exactly when `3 | n`**, over `n = 3..11` — S26's unexplained observation
  turned into a rule by S29's new `n = 3` row. In those cases the closure climbs a chain of proper
  subalgebras one rung per layer. Labelled a **CONJECTURE**; its next test, `n = 12`, is not runnable
  at `s26/q_dla.py`'s `CAP = 2^21`.
- **Scope correction to `s26/REPORT.md` V.9** (`s29/THEORY.md:1069-1073`): *"The algebra at the
  deployed cell is maximal"* is true at `n = 7, L = 3` (8128 = so(128)) and **false at the S27/S28
  register `n = 9, L = 3`, where it is 65535 = dim su(256), exactly half of dim so(512) = 130816.**

**The algebra is not the obstruction; the parameter count is** (`s29/THEORY.md:1075-1086`):
`su(N)` acts transitively on the unit sphere, so even `su(256) ⊂ so(512)` generates a group whose
orbit of `|0…0⟩` is the whole real sphere `S^511`. **Nothing is unreachable by the algebra.** What is
unreachable is everything outside the image of a smooth map `R^27 → S^511` — a semialgebraic set of
dimension ≤ 27 in 511, so a *generic* target has best squared overlap of order `P/D = 27/512 = 0.05`.
S28 measured 0.80–0.88, so the `J = 3` ground state is very far from generic. *"The correct 'provably
not reachable' statement for this ansatz is a dimension statement, not a Lie-algebra statement, and
any argument that reaches for controllability here is answering the wrong question."*

**The cut-rank result, which is the ansatz-side twin of the set-equality theorem**
(`s29/THEORY.md:1088-1104`): `L` entangling layers give bond dimension exactly `2^L` (verified to
3.3e-16); measured max-cut rank at `n = 9` over 200 random subsets per size —

```
|S|                              8    16    32    75   128
random subset  (median)          7    10    14    16    16
PREFIX {0..|S|-1}                1     1     1     2     1
```

> **The deployed encoding makes exactly the classical prefixes cheap and every other set expensive.**
> A prefix of the energy order is a rank-2 object; a generic 75-subset saturates the register's
> maximum. This is an **ansatz-side** mechanism for "the circuit reproduces the classical top-m",
> entirely independent of the CVaR clip the set-equality theorem argues from — **two separate
> reasons for the same fact.**

**The design variable nobody has used**, recorded so it is not lost: `s22/qcand_lib.py :: Encoding`
takes a free permutation `label` of the `2^n` basis indices and the deployment sets it to the
identity. *"The permutation decides which FAMILY of candidate sets is cheap. If a formulation wants
the circuit to select non-prefix sets, relabelling is the zero-cost lever, and it comes with the
matching control (a random relabelling) built in."* (`s29/THEORY.md:1108-1112`.) S30 did not test it.

---

## 20. WHAT THE REPORT'S §5 SHOULD SAY, IN ONE PARAGRAPH

Sprint 30 ran no CVaR-VQE. It ran a 9-qubit statevector circuit on all 126 targets, but only to
regenerate the meter's ORACLE rungs (`s30/results/s30_D_ladder_structs/`, 126 files, asserted to
1e-6 against S28), and the production anchor is computed with the quantum stage off
(`core/pipeline.py:179`, `:241`). What the sprint contributed to the quantum spine is three closures
and a theorem: **T1**, which shows the CVaR tail is always a prefix of the order induced by `∇V` at
the optimum, so the deployed stage's structural output is the single integer `m` and the charter's
L5 lead was mis-posed; **T1b plus two independent quadric measurements**, which price the only two
Hamiltonian classes that could have changed that order and find both reachable but
non-transferable — closed *before* either was built, which is the first time this project's cheap
pre-check has fired ahead of the spend; and **the bit accounting**, which inverts the charter's
question by showing the 7-bit register realises 36.6 bits and delivers under one, so the encoding is
not where the loss is. The reason none of this was carried into a circuit is S30-L19: on a
kind-matched, budget-matched ladder, ordering survives on 2 of 43 native-free channels and
**preference fails on all 43**, with a mechanism (local ΔR² −0.089 against global +0.600) that makes
the null a theorem on this instrument. An objective cannot be pointed at a target nothing can see.
The one number that shows this is not an abstraction is in the sprint's own cache: the identical
circuit reaches **0.2516 Å** on the built chain when its objective is the native and **3.4330 Å**
when its objective is the deployed native-free score, against a classical average's **3.2071 Å**.

---

*Lane W, 2026-09-20. Read-only on `s29/`, `s28/`, `s27/`, `s16/`, `s12/`, `core/`, `s30/LEDGER.md`,
`s30/STATE.md`, `s30/REPORT_S30.md`, `s30/BRIEF.md`. No pipeline or VQE compute spent; the only
computation performed for this document is a mean over 126 cached `.npz` rung labels.*
