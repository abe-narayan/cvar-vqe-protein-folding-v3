# S32 LANE Q — the derivations

Prereg `s32/PREREG_S32_Q.md` @ **`a8f9d6a7`**. Endpoint = **built-chain Cα RMSD, n = 126,
production 3.2105 Å**. Everything in this file that carries an Ångström is on the **CA point
cloud** unless it says otherwise, and every arm that reads the native is labelled
**ORACLE / NOT DEPLOYABLE**.

---

## Q0 — the target-invariance claim survives a test that could fail

S31 §5.1's own verification was run on **random floats, which never tie** — structurally incapable
of testing its own caveat (contract rule 5). This one runs on the real 126 pools plus two synthetic
inputs built to break the claim.

**Self-test (must break the ramp, and does):** a deliberately tied vector gives max deviation
**0.1221**; a deliberately unsorted vector gives **3.3830**; a sorted distinct vector gives
**0.0** exactly. The test can fail.

**Q0-H1 — the residual is EXACTLY the tie pattern.** On all 126 real targets:

```
sc[o] non-decreasing                        126 / 126
targets with >= 1 tie in the top-128        125 / 126
targets whose E is EXACTLY the ramp           1 / 126     (the one with no ties)
max |E - ramp|,  max over targets          0.0406840      S31 quoted 0.0407
max |E - ramp|,  mean over targets         0.0218
tied positions in excess of blocks, mean       9.556      S31 quoted 9.56
BLOCK-MODEL error  max |E - E(block sizes)|  0.0e+00      exact on 126/126
```

The block model reproduces `E` **bit-for-bit on every target**, and `E ≠ ramp` **iff** the target
has a tie, on 126/126. **H1 holds.** The Hamiltonian's entire target-dependent content is the
*composition of 128 into runs of equal score*, which is the duplicate-structure pattern.

**Q0-H2 — no other path.** `len(order) = 500` on 126/126 so `dim = 128` always; `(α,T)` takes
**exactly two values** over the whole instrument, selected by **fold**, not by target. The
optimisation's complete argument list is `(E, α, T, n, layers, iters, seed)` and nothing else in it
is target-specific.

**Q0-H3 — the residual, priced.** At the deployed `(α,T)`, comparing the real `E` against the exact
ramp, CA cloud, n = 126:

| arm | effect | SE | ×MDE | W/L/T | verdict |
|---|---|---|---|---|---|
| `Q0_TIE − Q0_RAMP`, convex readout | +0.0000 | 0.0000 | **0.26×** | 70/55/1 | NOT A RESULT |
| `Q0_TIE_SEL − Q0_RAMP_SEL`, selection readout | 0.0000 | — | — | **0/0/126** | **identically tied on every target** |
| `Q0_PERM − Q0_RAMP`, blocks relocated (8 draws, sizes matched) | +0.0004 | 0.0003 | 0.38× | 68/57/1 | NOT A RESULT |

`TV(p*_tie, p*_ramp)` is **3.3e-03 mean, 1.8e-02 max**; the selected medoid is **identical on
126/126**.

> ### VERDICT: S31 §5.1 is NOT falsified, and is now sharper.
> The control matched to the operator's own space (`Q0_PERM`, which keeps the block sizes and
> relocates them) is **the same size as the real pattern**. So the residual channel is the
> *existence* of duplicates, **not which candidates are duplicates** — it carries the integer
> "how many distinct structures are in the top-128" and nothing else. The deployed Hamiltonian is a
> **global 128-number hyperparameter plus one integer per target**, and the integer is worth
> 0.26× MDE through the readout that ships and **exactly nothing** through the readout that selects.

---

## Q1-T1 — `â` and `μ` are the same object. S31 §20.3 is resolved.

S31 §20.3 states as *"the sharpest question this sprint produces, and it is stated as open"*:
whether the per-candidate quality `â` and the pool common mode `μ` are one requirement or two.
They are one.

**Setup.** Candidates `W_1..W_K` in a common frame, flattened to `R^d` (`d = 3n`); native `t`;
weights `w` with `Σ_x w_x = 1` (non-negativity not required). `U` is the `K × d` matrix of rows
`W_x`, so the emitted point is `x(w) = Uᵀw`.

**(I) The readout identity.** With `a_x = ‖W_x − t‖²` and `B_xy = ‖W_x − W_y‖²`,

```
||U'w - t||^2  =  <w,a> - (1/2) w'Bw
```

*Proof.* `Uᵀw − t = Σ_x w_x (W_x − t)` by `Σw = 1`; expand and substitute
`⟨W_x−t, W_y−t⟩ = ½(a_x + a_y − B_xy)`. ∎ **Verified to 4.1e-14 relative** over 126 targets × 80
draws, simplex and signed-affine.

**(II) `a` is affine in `t`.** `a_x = ‖W_x‖² − 2⟨W_x, t⟩ + ‖t‖²`. The first term is native-free; the
third is **constant in `x`**, and `Σw = 1` makes an `x`-constant an additive constant of the
objective, which cannot move an argmin over the simplex. **Verified to 3.0e-15 relative.**

> **Therefore the readout-relevant content of the quality vector `a ∈ R^K` is the linear functional
> `g_x = ⟨W_x, t⟩`, i.e. `P_aff{W} t` — the native's coordinates in the candidate set's own affine
> hull.**

**(III) The bijection with the common mode.** `μ = X̄ − t` with `X̄` the pool mean (native-free), so
`t = X̄ − μ`, and by (II) `a` is an affine image of `t`, hence of `μ`. Conversely, `a` on any
`rank(aff{W}) + 1` candidates in general position determines `P_aff{W} t` by least squares.

```
              NATIVE-FREE AFFINE BIJECTION
   mu  <----------------------------------->  t  <----------------->  a
        t = Xbar - mu                              a_x = c_x - 2<W_x,t> + ||t||^2
```

**Numerical falsifier, run:** replace `t` by `P_aff{W} t` everywhere and re-solve the readout's
convex program. Predicted identical; measured **max |Δw| = 2.5e-14**, **max Δx = 6.1e-14 Å**.
`rank(aff{W_x})` is **≈ 32–33** against `d ≈ 38`, and `‖t − P_aff t‖ = 2.0e-14` — the native already
lies in the candidates' affine span, which is S31 §10.5's observation arriving from the other side.

> ### THEOREM Q1-T1. The per-candidate quality vector and the pool common mode are the same
> ### object up to a known native-free affine bijection, and the sufficient statistic for the
> ### entire sum-to-one readout problem is ≈ 33 real numbers, not 128.
>
> **Consequence for the charter.** A per-candidate quality estimator with real in-band skill **is**
> a structure predictor, and a common-mode corrector **is** a per-candidate quality estimator.
> There is one missing channel this project has been describing in two vocabularies, and no
> arrangement of the pipeline supplies one from the other. *This is a closure, not a lead.*

---

## Q1-T2 — the sensitivity question, answered analytically

*(This is the question the coordinator asked: **how precisely does `a` have to be known for the
optimum `w*` to be useful?** It is answerable on the convex program with no circuit.)*

By Q1-T1 the program is `min_{w ∈ Δ} ‖Uᵀw − t‖²`, i.e. **the Euclidean projection of the native
onto the convex hull of the candidates**, `x(t) = P_C(t)` with `C = conv{W_x}`.

### (a) Global: gain ≤ 1, always

Projection onto a convex set is firmly non-expansive, so `‖P_C(t̂) − P_C(t)‖ ≤ ‖t̂ − t‖`. **The
readout cannot amplify an error in the quality/structure estimate. It also cannot attenuate one
that lies in the hull.** With `d = dist(t, C)` (the ORACLE hull floor):

```
d  <=  || P_C(that) - t ||  <=  d + eps        where eps = || that - t ||
```

### (b) Local: the window is `|S| − 1` real numbers wide

Let `S = supp(w*)`, `A = U_Sᵀ`, `G_SS = A'A`. Stationarity on the support plus `Σw = 1` gives, for
a perturbation that does not change `S`,

```
dw_S = M dg_S ,     M = G_SS^-1 - (G_SS^-1 11' G_SS^-1)/(1' G_SS^-1 1)
dx   = A M A' dt  =  P_aff{W_x : x in S} dt
```

*(the rank-one correction is exactly `uu'/u'u` with `u = A G_SS^{-1}1` the unique vector in
`span(A)` with `⟨W_x,u⟩ = 1` on `S`, i.e. the affine-hull normal — so `P_span(A) − P_u` is the
projector onto the **direction space of the active candidates' affine hull**, of dimension
`|S| − 1`.)*

> ### THEOREM Q1-T2. `∂x/∂t` is an orthogonal projector: gain **exactly 1** inside the active
> ### candidates' affine hull (dimension `|S| − 1`), **exactly 0** outside. In terms of the quality
> ### vector, `∂x/∂a_x = 0` for every candidate outside the support.

**Verified by finite differences, n = 126 × 4 draws each:**

```
|| dx - P_aff(S) dt || / || P_aff(S) dt ||     7.9e-08   (mean)     <- the identity
gain along a PURE in-hull direction            1.0000000  sd 3.1e-09
gain along a PURE orthogonal direction         4.3e-08               <- exactly blind
support size |S|, unconstrained convex optimum   see s32_Q1_sufficiency.json
```

### (c) The three consequences that decide the route

1. **How precisely must `a` be known: to first order, only on the support, and only along
   `|S| − 1` directions.** Everything else in `a` is exactly invisible. That is the minimal
   information requirement charter §12 asks for, and it is **~5 real numbers locally, ~33
   globally** (across active-set changes).
2. **The readout is dominated by its own input.** `‖P_C(t̂) − t‖ ≤ ε + d`, while emitting `t̂`
   directly costs `ε`. So the projection readout beats direct emission **iff** the estimate's error
   exceeds `d` *and* lies outside the hull's affine span — measured crossover at `ε ≈ 2.2 Å`
   against a hull floor `d ≈ 1.96 Å` (CA cloud, **ORACLE / NOT DEPLOYABLE**). *A structure
   estimate good enough to make the readout worth solving is already good enough to emit.*
3. **Gain exactly 1 means no noise suppression.** There is no regime in which a noisy `â` is
   cleaned up by the convex program. This is why S31 measured "solving the objective exactly
   reshuffles the answer everywhere and buys nothing": the program is a faithful, non-contracting
   transcription of whatever quality estimate it is handed.

---

## Q2 — the fifteen items, for each objective that was seriously considered

Charter §16. **Item 13 (classical equivalent) is decided first**, because it is the item that
killed S31 and no compute is spent on an objective that fails it.

### Objective A — CVaR over a posterior on candidate quality (lead (a))

| # | item | |
|---|---|---|
| 1 | definition | `min_{w∈Δ} CVaR_α( ⟨w,a⟩ − ½w'Bw )` with `a ~ π` a posterior |
| 2 | physical | risk-averse projection of an uncertain native onto the fragment hull |
| 3 | state | amplitudes ↦ a distribution over quality hypotheses |
| 4 | basis | candidate index (or a discretised quality vector) |
| 5 | Hamiltonian | diagonal in a sampled `a`; **not** an operator whose expectation is the objective |
| 6 | CVaR meaning | risk over *estimator error*, not over candidates |
| 7 | gradient | `∇_w = â_tail − Bw`, available in closed form |
| 8 | landscape | **convex on the simplex** (Q1-T1 makes it `‖Uᵀw − t̂‖²` plus a norm) |
| 9 | minima | unique global; no local minima |
| 10 | candidate quality | `a` *is* candidate quality, by definition |
| 11 | RMSD | exact — it **is** the squared cloud error |
| 12 | information | the whole of `a` is `P_aff t`: ≈ 33 reals (Q1-T1) |
| **13** | **classical equivalent** | **Gaussian `π` ⇒ `⟨w,â⟩ − ½w'Bw + φ(α)√(w'Σw)` = a second-order cone program: convex, unique optimum, milliseconds. CLOSED.** |
| 14 | failure mode | `a` unknown; setting `a = const` reduces it to "maximise spread", which S31 measured at **+0.1436 Å, 1.22× MDE, WORSE** |
| 15 | deployment condition | a native-free `â` with in-band `ρ > 0` — which by Q1-T1 is a structure predictor |

**VERDICT: closed by classical equivalence, and the closure is stronger than "a circuit is
unnecessary" — Q1-T2 says the program does not even need solving accurately, because its output
inherits its input's error at gain exactly 1.**

### Objective C — sparse convex combinations, `s`-of-`K` (lead (c))

| # | item | |
|---|---|---|
| 1 | definition | `min ‖Uᵀw − t‖²` s.t. `w ∈ Δ_K`, `|supp(w)| ≤ s` |
| 2–4 | physical / state / basis | a bitstring is a **subset**; register width `K`, dimension `2^K` |
| 5 | Hamiltonian | `E(x) = min_{supp(w)⊆x} ‖Uᵀw − t‖²` — **diagonal in the subset basis** |
| 6 | CVaR | risk over subsets; a genuine lower-tail objective |
| 7–9 | gradient / landscape / minima | non-convex in `x`; best-subset selection |
| 10–11 | quality / RMSD | exact, as A |
| 12 | information | as A: the objective needs `P_aff t` |
| **13** | **classical equivalent** | **see the theorem below — the constraint is SLACK** |
| 14 | failure mode | the objective needs `a`; and the constraint does not bind |
| 15 | deployment condition | same missing channel as A |

> ### THEOREM Q1-T3 (the cardinality constraint does not bind).
> Let `w°` be the unconstrained minimiser over `Δ_K` and `s° = |supp(w°)|`. The `s`-cardinality
> feasible set is a **subset** of `Δ_K` containing `w°` whenever `s ≥ s°`; therefore `w°` is the
> exact global minimiser of the constrained problem. **For `s ≥ s°` the "sparse" problem IS the
> convex program**, solvable in milliseconds with a KKT certificate.

Measured `s°` (unconstrained convex optimum, **ORACLE / NOT DEPLOYABLE**): see
`s32_Q1_sufficiency.json :: Q1C_cardinality_ORACLE`. **This is the item that decides whether the
charter's 2.10 Å sparse headroom is a combinatorial prize or a convex one.**

**VERDICT: the combinatorial hardness is an artefact of quoting `s` below the solution's own
sparsity. It escapes S31 §5.2's obstructions 1 and 3 — the subset basis gives every bitstring a
definite energy, and its Hilbert dimension is `2^K`, not `K` — and is then closed by a different
argument: an inactive constraint plus an unknowable objective.**

### Objective B — reconstruction branch selection (lead (b), owned by lane R)

`core/project.py:100-118`: *"A CA trace admits **two** ideal-geometry torsion solutions at
near-equal objective distance, one Ramachandran-plausible and one not … the structures this stage
returns on those targets are **not determined by the objective; they are determined by the
arithmetic**."*

Q's contribution is the **formulation**, not the pricing:

- **If the branch variable is per-target with a handful of values**, the decision is exhaustively
  enumerable and CVaR over it is an argmin: at `T = 0, α = 1` the minimiser of `CVaR_α(E;p)` is the
  argmin vertex, and at `α < 1` the argmin set is the face `{p : p_{x₀} ≥ α}` of positive volume, so
  the objective does not determine `p` at all (S31 §7, measured). **CVaR adds nothing over
  argmax for a fixed, known energy vector.** Saying otherwise would be S31 §5's objection
  reappearing in a new basis, which is precisely what this lane exists to prevent.
- **CVaR becomes meaningful only if the branch SCORE is itself uncertain**, and then it is the
  1-D classical rule *"pick the branch minimising score + λ·(score uncertainty)"*. Still classical.
- **The formulation earns a circuit iff the branch variable is PER-RESIDUE and the energy is
  non-separable across residues** (chain closure couples them), giving `2^{n_res}` with a diagonal
  per-shot energy. That breaks all three of S31 §5.2's obstructions at once.

**Concrete diagnostic handed to lane R, not duplicated here:** when two multi-start branches differ,
**do they differ at a sparse set of residues or globally?** Sparse ⇒ the decision is per-residue and
combinatorial; global ⇒ it is a binary label and CVaR is an argmax. This is one line on artefacts
lane R is already generating, and it decides whether objective B is a quantum problem.

**Whatever the answer, the scoring channel is the interesting half**, and it is *not* Q's to price:
the branch decision is **chiral**, and every native-free ranker this project has tested is a
distance-map function and therefore achiral — so it is the first decision in the pipeline that the
existing score family is structurally unable to see.

---

## Q3 — non-diagonal Hamiltonians: which obstruction each basis breaks

Charter §13 permits reopening with a **named mechanism**, quoting the original closure.

| S31 §5.2 obstruction | quoted | does a new basis break it? |
|---|---|---|
| **1. CVaR needs an energy per shot** — *"Only a diagonal `H` gives every measured bitstring a definite eigenvalue"* | kills `diag(zrank) − λ·W(block)` on a **candidate-index** register | **Not an obstruction to changing what the bitstring labels.** A subset basis or a per-residue branch basis is diagonal **by construction** and has genuine combinatorial structure. Obstruction 1 forbids adding hopping; it does not forbid a QUBO. |
| **2. The forced operator is quartic in ψ** — *"`H[w] = diag(â) − B`, whose VMC local energy is … a mean-field operator, quartic in ψ"* | arises because the readout weights **are** the state's probabilities, so the objective is quadratic in `p = |ψ|²` | **Broken by decoupling the weights from the probabilities.** Let the bitstring select a *discrete object* and solve the continuous weights classically inside `E(x)`; then `⟨E⟩ = Σ_x p_x E(x)` is linear in `p` and is an ordinary diagonal expectation. This is the only one of the three with a clean escape. |
| **3. Dimension counting** — *"a candidate-index register has Hilbert dimension equal to the candidate count … No Hamiltonian on a candidate-index register can be classically hard"* | correct, and it is about the **register**, not the Hamiltonian | **Broken by the subset basis** (`2^K`, and `C(500,10) ≈ 2.5e20` is not enumerable) and by a per-residue basis (`2^{n_res}`) — though at `n_res ≈ 12–16` the latter is enumerable on this instrument, which is the memory entry *"exhaustive enumeration closes the search half"* arriving again. |

> **So obstructions 1 and 3 are properties of the candidate-index register, not of CVaR-VQE, and
> both are escaped by changing what a basis state means. Obstruction 2 is escaped by solving the
> continuous part classically inside a diagonal `E(x)`.** The family is **reopened** on that
> mechanism — and then closed again, one level down, by Q1-T1 (the objective needs `P_aff t`) and
> Q1-T3 (the cardinality constraint that would make it hard does not bind). *The closure moved from
> "a circuit cannot be built" to "a circuit can be built and there is nothing for it to compute",
> which is a different and more useful statement.*

---

## Q4 — the five-bit result is priced in the wrong currency

Charter §41 asks whether five bits are necessary, sufficient, and an artefact of the oracle
construction.

**By Q1-T2 the decision is not a selection at all.** The readout reads the native through an
orthogonal projector of rank `|S| − 1` at gain exactly 1 — a **continuous** window. A bit count
prices a *selection alphabet* (S31 §6: 6.886 bits mean, the distinct-candidate count); it does not
price this decision, and the two must not be differenced.

**The registered pricing curve, in real numbers rather than bits.** Truncate the native's
representation to the top `r` directions of the pool's **own** spread (the ordering is native-free;
only the coefficients are ORACLE), re-solve the convex readout, and emit. **ORACLE / NOT
DEPLOYABLE**, CA cloud — see `s32_Q1_sufficiency.json :: Q4_curve_ORACLE`. `r = 0` is the uniform
mean of the 128 and is the curve's own baseline; `r = rank(aff{W}) ≈ 33` is the hull floor.

**Answers to §41's questions, from the derivation:**

- *Could the information live in amplitudes rather than basis states?* It already does — the
  quantity needed is `≈33` continuous coefficients, not a label. **But representational capacity was
  never the binding constraint:** the numbers required are the native's coordinates in the pool's
  basis, and no encoding manufactures them.
- *Could a nonlinear observable amplify it?* **No.** Q1-T2 gives gain exactly 1; a nonlinear
  readout of a state that does not contain `P_aff t` cannot create it, and by consequence (c) there
  is no noise-suppression regime to exploit.
- *Is the five-bit structure an artefact of the oracle construction?* **Yes, in the specific sense
  that it prices an alphabet and the decision is continuous** — and the honest replacement is the
  `r`-reals curve above, which is a different currency and is never differenced against a bit count.

---

## The charter's own question, answered plainly

**Charter §14: *"Can a genuine CVaR-VQE be designed whose quantum state and objective actually
contain information that can improve RMSD?"*** Charter §58 frees this lane from defending the
spine. So: **on this instrument, no — and the reason is a property of the instrument, not of
CVaR-VQE.** Here is the derivation rather than the opinion.

A CVaR-VQE does real work on a decision only if **all five** of these hold:

| | condition | why |
|---|---|---|
| **A** | the decision space is discrete and **too large to enumerate** | otherwise `argmin` by brute force |
| **B** | the energy `E(x)` is **per-shot computable, target-dependent, native-free** | otherwise there is no diagonal Hamiltonian, or no target in it |
| **C** | choosing better in the space **lowers built-chain RMSD** | charter §2 |
| **D** | **no cheap exact classical algorithm** — not convex, not separable, not greedy-optimal | charter §16 item 13 |
| **E** | `E` is a genuine **random variable**, so the lower tail differs from the minimum | otherwise CVaR is `argmin` in risk notation |

Every discrete decision this architecture contains, scored against those five:

| decision | A | B | C | D | E | verdict |
|---|---|---|---|---|---|---|
| **candidate index** (deployed) | ✗ 128 | ✗ target-independent (Q0) | — | ✗ argmin | ✗ | dead 4 ways |
| **subset / sparse `s`-of-`K`** | ✓ `C(500,10)≈2.5e20` | ✗ `E(x)` needs `a` = ORACLE | ✓ | **✗ Q1-T3: the constraint is slack** | ✗ | dead |
| **reconstruction branch, per target** | ✗ a handful | ✓ | ✓ open (lane R) | ✗ enumerable | ✗ | argmax in quantum notation |
| **reconstruction branch, per residue** | **✗ `2^n_res ≤ 65536` on 126/126** | ✓ | ✓ open (lane R) | ✗ enumerable | ✗ | **blocked only by chain length** |
| **fragment assembly** (which window where) | ✓ `K^n_res` | ✓ | ✓ | ✓ genuinely NP-hard | ✗ | **does not exist at this length** |

Measured on the instrument: **`n_res` is 9–16, mean 12.96, so `2^n_res ≤ 65536` on 126 of 126
targets.** Contract-adjacent memory: *exhaustive enumeration closes the search half — the budget
exceeds the `2**n` latent on 75/126.*

> ### The two properties a problem would need, stated so they can be checked rather than argued.
>
> **P1 — a decision space that GROWS WITH THE TARGET and outruns enumeration.** On 9–16-residue
> peptides every decision here is a total ordering over ≤ 500 objects, a convex program in
> disguise, or a search of ≤ 2^16. **Fragment assembly and per-residue branch selection are the
> two places where a genuine combinatorial problem appears — and both only exist at chain lengths
> where one retrieved fragment no longer spans the target.** That is the charter's own closing
> instruction (*"maybe test on longer proteins"*) arriving as a derived requirement rather than a
> suggestion.
>
> **P2 — an energy that is GENUINELY STOCHASTIC, so the lower tail is not the minimum.** For every
> decision in the deployed pipeline `E(x)` is a deterministic function of the bitstring, so
> `CVaR_α` reduces to a reweighting of a fixed vector whose minimiser is a face of the argmin set
> (S31 §7, measured). **CVaR earns its name only where the energy is a sampled quantity** — a free
> energy from a finite MD sample, a physically noisy observable. Charter §34 (*"CVaR of
> conformational free energy"*) is the one place in the charter where P2 could be met, and it is
> met by the *sampling*, not by the physics vocabulary.
>
> **P1 and P2 must hold together.** P1 alone gives a quantum optimiser with a deterministic
> objective — QAOA, not CVaR-VQE. P2 alone gives risk-sensitive selection over a small set — a
> classical one-dimensional rule. *This project has never had either.*

---

## Scope of Q1-T1 and Q1-T2 — what the theorems do NOT cover

Stated explicitly, because S31's G1 was over-applied and this derivation is now load-bearing.

1. **`Σw = 1` is assumed.** Identity (I) and Q1-T1's bijection hold for **signed affine** weights
   too (verified on signed draws, 4.1e-14). Q1-T2's *active-set* formula does **not**: with no
   non-negativity there is no active set, the window is `rank(aff{W}) ≈ 33` everywhere, and by
   S31 §10.5 the affine hull already spans the residual space — so the theorem is **vacuous**, not
   false, on the affine readout.
2. **No cardinality constraint is assumed.** With a **binding** `|supp(w)| ≤ s` the feasible set is
   a non-convex union of faces; the global 1-Lipschitz statement fails and only the local,
   fixed-support projector formula survives. Q1-T3 establishes that the constraint is **slack** at
   the sparsity the solution chooses for itself, which is what makes the convex statement apply
   here — it is **not** a claim about sparse readouts in general.
3. **The emitted object is assumed to lie in the candidates' affine hull.** Arms that **re-embed**
   are outside the hypotheses entirely — S31's `AVG_SEP` left the hull and measured **0.045 Å**
   from it, and nothing here applies to it.
4. **A fixed common frame is assumed.** The endpoint does one further Kabsch; both fixed-frame and
   re-superposed values are reported and the gap is small but not zero.
5. **"Gain exactly 1" is a statement about the CA POINT CLOUD and must not be carried to the
   built chain.** The projection to a chain is **not** 1-Lipschitz — S32-L4 measures it as
   discontinuous in its input at one float64 ULP, ~1e13 amplification. *The contraction argument
   stops at the cloud, and saying otherwise would be this project's signature defect: a number
   re-used across a boundary its definition does not cross.*
