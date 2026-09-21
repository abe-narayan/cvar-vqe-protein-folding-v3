# Sprint 31 — CVaR-VQE Protein Folding: the objective, the readout, and where the information is not

**Status: DRAFT IN PROGRESS.** Sections marked `[PENDING]` await lanes still running. Assembled as
results land so nothing is reconstructed from memory at the end.

Branch `s26` · instrument: `tuning126`, 126 targets, 9–16 aa · endpoint: **mean built-chain Cα RMSD**
Charter: `s31/BRIEF.md` (verbatim, 1,660 lines) · Ledger: `s31/LEDGER.md` · State: `s31/STATE.md`
Contract: `s31/S31_CONTRACT.md` (32 rules) · Verifier: `s31/s31_verify.py` · Multiplicity:
`s31/MULTIPLICITY.md`

---

## 0. The answer, up front

**The endpoint did not move. It is 3.2105 Å.** Nothing was deployed and nothing earned deployment:
**not one arm in the sprint reached its own MDE in the helpful direction.** The charter's primary
target of < 3.00 Å and its ambition of < 2.50 Å are both unmet.

What the sprint did instead was **dissolve the question it was given and replace it with a sharper
one that points the opposite way.**

### 0.1 The charter's central question does not have an answer, because it is not well posed

The charter asks what objective and Hamiltonian a CVaR-VQE should optimise so that quantum
optimisation spends its capacity on endpoint-relevant selection. Reading the code rather than the
record moved the diagnosis twice, and the second move ended the question:

> **The Hamiltonian's diagonal is `E = _zrank(pool["sc"][o])` — the standardised *rank* of 128
> sorted values — so it is THE SAME VECTOR ON EVERY TARGET**, to within 0.0407 in max-norm.
> `p*` and the circuit's output are therefore **one fixed weighting curve per α**, two values
> selected by fold, with `sd(H(p*)) = 1.4e-4` bits across 126 targets. **The stage answers *"what
> fixed weight should rank `k` receive?"* — a 128-number global hyperparameter, not a per-target
> computation.**

That is strictly stronger than either capacity theorem the sprint began with, and it explains every
other negative at once: why the objective does not point at good solutions (*it is not a function of
the target*), why α is inert on three of five folds by construction, why the circuit's 0.930-bit
optimisation gap **cancels in the mean**, and why the convex optimum over that family converges to
the uniform average production already emits.

**And the deployed objective is not a quantum problem.** `F(p) = CVaR_α(E;p) − T·H(p)` is a **convex
program with a closed-form global minimiser** pinned by one scalar (Rockafellar–Uryasev + Sion,
derived twice independently). `run_cvar_vqe` is **strictly worse in 126/126 targets** at the deployed
settings; `p*` costs **0.0012 s against 0.0985 s**. Substituting it changes the selected candidate on
**66 of 126 targets** and is worth **−0.0112 Å at 0.19× MDE — a null.** *Solving the objective
exactly reshuffles the answer everywhere and buys nothing.*

### 0.2 The bottleneck was stated backwards, and the correction is the sprint's most consequential result

S30 concluded — and this sprint's own opening contract encoded — that the next channel must supply an
observable whose error is **incoherent** with the pool's common mode. Measured with ORACLE magnitude
**and** a perfect direction, an incoherent correction is worth **+0.0747 Å: harmful.**

```
correcting ONLY the along-mu half    -0.8102 A   3.29x MDE, 5/5 folds, 121W/5L    <- the whole prize by SIZE
correcting ONLY the orthogonal half  +0.0747 A   0.87x MDE, 57W/69L               <- worth less than nothing
      ^ the SIZE is real; that it is a separate MECHANISM is NOT -- the along arm is NOT MEASURED against the shrink curve
```

And the reason is an **identity, not a fit**: `mu_hat = mu − y` forces `y_perp = −P_perp(mu_hat; mu)`
to 9.3e-15, so out-of-fold `R²` is **+0.9403** on the orthogonal half and **−0.0051** on the
along-`mu` half. **The predictable half is not predicted — it is observed.**

**And `mu` survives its own falsifier.** The asymmetry could have been magnitude or generic
geometry; an **energy-matched** arbitrary direction — same captured fraction of `y`, otherwise
arbitrary — gives an along-half worth **0.13 Å less** and a perp-half worth **0.38 Å more** than
`mu`'s (2.07× and 2.53× MDE, both 5/5 folds). **Two qualifications belong here and not only in
§20.1.** *First*, that same arbitrary direction still shows along-beating-perp by −0.3734, so **≈42%
of the raw −0.8848 asymmetry is generic geometry** and `mu` supplies 58% — ***`mu` doubles an
asymmetry that is already there rather than creating one.*** *Second, and larger:* `y_along(mu)`
carries **84% of `y`'s energy**, and against its own magnitude control — the shrink curve — it is
**NOT MEASURED** (0.94× and 0.64×, 65W/61L at c = 0.90). ***The controlled half of this result is
the negative one — that the orthogonal complement is specifically worthless. The positive half
restates, at matched magnitude, that most of a perfect correction is most of the prize.*** And that
negative half is exactly what charter §14 and contract rule 29 proposed to build on. And **three
quarters of a perfect correction buys 97% of its benefit**, so the channel §21 asks for must point
the right way, not be accurate.

> ### What can be predicted is the component ORTHOGONAL to the common mode, and it is harmful.
> ### What would help is the common mode itself, and it is unpredictable.

**S30's slogan survives verbatim; its mechanism was inverted.** And two lanes that never spoke
measured the same fact from opposite ends: the tail's pools are **coherently wrong, not diversely
wrong** — `S/B` falls 0.78 → 0.42 while the distinct-candidate count is unchanged, *they are
displaced together* — and **the entire recoverable prize lies along that displacement**, 4×
concentrated on the tail.

**So the next sprint's question is the inverse of this one's:** find information about the pool's
**common mode** that does not come from the pool. It is **priced** (−0.8102 Å ORACLE, which would put
the chain at **2.4025 Å**, past the *ambitious* target), the **conversion is free and tuning-free**
(an exact convex program, §3(ii)), and the **failure mode is named** (non-identifiable from pool data
at any K, so a channel that reads the pool cannot supply it).

### 0.3 What closed, and how much of it by derivation

Thirteen directions, **most by theorem or by derivation before compute was spent** — the free-energy
stage and the entire elastic-network/normal-mode family (achiral, hence distance-map readings by G1);
non-diagonal Hamiltonians (three independent obstructions); ADAPT-VQE (its selection rule is
*undefined* for CVaR); the torsion channel (the chiral escape is real **and only works on a problem
the pipeline does not have**); charter §14 by **provenance** — 92.9% of the benchmark is
NMR-determined and for **117/126** the deposited coordinates *are* a fit to the deposited restraints,
so any NMR observable is **ORACLE through a different door**; index redesign; the per-target prefix
length (a matched random-subset family reaches **149%** of it); and the constrained-affine readout,
which **cannot be selected by ceiling and cannot be selected by fit**.

**The readout is now exact and closed.** For any `Σw = 1`,
`‖Σ w_x W_x − t‖² = ⟨w,a⟩ − ½ w'Bw`, verified to 1.66e-11 — **half of it is free, the only unknown is
per-candidate quality, and the native-free channel is 3.99% of the ORACLE cross term.** Every
native-free quality estimate measured has **in-band skill that is zero or the wrong sign.**

### 0.4 How much of this report is negative, and the pattern in the errors

Nine lanes, 21 ledger entries, ten pre-registrations. **Six registered predictions fired against the
lane that wrote them, three of those six mine.** The sprint's only live deployable candidate,
`AVG_SEP`, is **refuted at +0.4609 Å, 2.34× MDE** — and the *"0.4–0.6 Å better"* that circulated
mid-sprint was **the magnitude of that deficit with the sign inverted**, read off partial rows the
lane itself had refused to average.

**Fourteen claims of mine were withdrawn, and the distribution is the finding:**

> **Every single-lane result held. Every cross-lane synthesis of mine failed** — four for four before
> an adversary was assigned, then eleven more defects, four severe. **A single-lane claim is audited
> by the lane that owns the data; a cross-lane claim is audited by nobody**, because each lane sees
> only its own half and assumes the other was checked.

And the shape is one shape. Every defect the adversary found — and every one of mine — is **a
quantity transplanted across a boundary its definition does not cross**: an object, a search size, a
bit budget, a basis, a data regime, a moment, a key, a sign.

> **A number carries its definition, not just its value.**

With one instrument corollary that cost the sprint a headline: **an audit that checks for an absent
label cannot catch a wrong one.** The verifier passed 97/97 on *"is a basis named?"* while the
headline carried a **misnamed** basis. It now asserts each number against the value on the basis it
names, and ships with a self-test on the defect that motivated it.

---

## 1. Endpoint and reproducibility

### 1.1 The endpoint did not move

**Production remains 3.2105 Å** (`s29/results/s29_O_chain_rows.jsonl :: item=prod`), mean built-chain
Cα RMSD over the 126 `tuning126` targets. The CA point cloud is **3.0483 Å** and the *set mean* is
**3.5507 Å** — three different objects, and this report names the basis in the same sentence as every
number.

`benchmark60` is **sealed**; its single pre-registered pass is spent and two guard sites enforce it
(`core/pipeline.py:386-388`, `core/bench.py:1074`, `core/pipeline.py:1664`). Folds and clusters were
not regenerated.

### 1.2 The projection is deterministic and chaotic — there is no seed, and the record said otherwise for three sprints

S28, S30's ledger, S30's report, S30's verifier and this sprint's own charter all describe *"an
unpinned multi-start projection seed."* **There is no RNG anywhere on the projection path.**
`core.project.fit_multi` loops over four fixed starts and takes a strict argmin; `fit_prior` is a
deterministic L-BFGS-B call. Reprojecting the same cloud twice is **bit-identical — max |ΔCA|
exactly 0.0**, unchanged under BLAS thread count.

**The real defect is conditioning, and it is worse than a seed would be:**

```
1A13's four lam=0 objectives   0.5091924865 / 0.5091925170 / 0.5091924488 / 0.5091924188
                               -> spread 1e-7, i.e. NOISE
the same four branches at lam=0.3   0.947 / 0.519 / 0.562 / 2.396
                               -> spread 1e-1
```

A decision taken where the objective **cannot** discriminate fixes an outcome where it **can**.
Measured across the benchmark: the median λ=0 branch margin is **3.1e-7**, and **73 of 126 targets
(58%)** have their branch chosen below a 1e-6 margin. A **1e-14 relative** perturbation of the input
cloud moves the emitted chain by a median of **1.6e-3 Å with a 0.511 Å tail** (2LNG), reproducing the
historical 0.517 Å discrepancy exactly. **Amplification ~1e13. Not one target of 126 is unchanged to
1e-9.**

**Three rules follow and were binding on every lane from the moment they were measured:**

1. Reprojection is reproducible **only from bit-identical clouds**. Both sides of any built-chain
   contrast must be projected **in the same job from the same stored clouds**.
2. **Canonical stays 3.2105 Å** — the run the endpoint was declared from, and the only one whose
   clouds are persisted per target. Not retro-fitted.
3. **No built-chain MEAN claim below 0.0107 Å** — the spread across the five circulating
   values. **And the per-target floor is far higher.** Lane P measured it directly (S31-L18)
   by running the *same operator in two implementations* — `average_weighted(uniform,
   top-75)` against `coordinate_average(top-75)`, whose **clouds agree to 3.6e-14**:

   ```
   built-chain MEAN     3.2108 vs 3.2105, paired +0.0003, 0.03x MDE  -- a null
   PER TARGET  |d|      mean 0.0134   median 0.0026   p90 0.0329   MAX 0.2285 A
                        exactly zero on 0 of 126; above 0.0107 on 28 of 126 (22%)
   ```

   > **A 1e-14 input difference — the same method written twice — moves the built chain by
   > up to 0.23 Å on individual targets while cancelling to nothing in the mean.** Any
   > per-target chain statement below **~0.03 Å** is indistinguishable from re-running the
   > same computation in a different implementation. This is also a **matched null** for
   > the chaos result above, arriving from shipped code rather than an injected
   > perturbation, and with a distribution rather than a single worst case.

**A consequence nobody had drawn:** the measured 0.92 cloud→chain transfer is **an average over a map
that is locally chaotic on a substantial minority of targets.** A lane proposing a cloud-level gain
must expect a non-smooth chain response there.

### 1.3 The reference itself is uncertain — and it does not explain the tail

**92.9% of `tuning126` is NMR-determined** (115 solution NMR, 5 electron crystallography, 4 X-ray, 2
solid-state NMR; RCSB GraphQL over all 126). The manifest reference is *"deposited coordinates, model
1"* — an arbitrary member of an ensemble. The ensembles are genuinely wide: mean pairwise CA-RMSD
between deposited models **1.0823 Å**, with model 1 sitting **0.6965 Å** from its own ensemble medoid
on average (max 4.294); 55/111 resolved targets have spread > 1.0 Å and 15/111 > 2.0 Å.

**It does not explain the tail:**

```
corr(production RMSD, ensemble spread)  +0.1118   95% CI [-0.0762,+0.2921]   SPANS ZERO
worst 18    mean RMSD 6.0291   mean ensemble spread 0.9672
other 93    mean RMSD 2.6842   mean ensemble spread 1.1046
tail-minus-rest  -0.1374   SE 0.2421   MDE 0.6783  ->  0.20x MDE   NULL
```

**The tail's targets have, if anything, narrower deposited ensembles.** The same null holds in all
four other arms. FAIL18 is real failure against a reference no worse determined than any other
target's.

**What survives is a caveat on the absolute number only** — and it is **~2× larger than first
stated, for two reasons lane V found.** The reference term is **not uniform** (mean 0.6965 but
**median 0.4565, sd 0.7767, CV 1.12, max 4.2943, and 13/111 exactly 0**), and the quadrature needs
the **RMS** of the per-target reference RMSDs — **1.0407**, not the mean — because
`d_obs² = d_true² + r²` is an identity **in squares** applied to **first** moments. So the endpoint
inflation is **+0.1645 Å, not ~0.08**, and the attenuation at `d = 1.0` is **30.7%, not 17.9%**.

**And paired deltas are protected for a different reason than the one first given.** They do not
survive because the term is uniform — it is not. They survive because **both arms share the same
reference on the same target**, which holds regardless of uniformity, and the effect on a delta is
**attenuation, not cancellation**: `observed = true × d/√(d²+r²)`, i.e. **2.3% at the endpoint**
(0.0005 Å on this project's confirmed 0.0221 Å effect — immaterial) rising to **30.7% at d = 1.0**.

**A positive that could have been claimed and was not:** **13 of 111 deposited `model 1`s are
exactly their own ensemble medoid (11.7%)** against ~5% expected by chance at ~20 models — depositors
commonly order NMR models by agreement. **Direct evidence that the reference is a better-than-random
ensemble member**, which strengthens the "change nothing" recommendation on its own terms.

*The original framing, for the record:* a **uniform ~0.70 Å reference term** — a
perfect predictor aiming at the ensemble medoid still scores ~0.70 Å against model 1. In quadrature
`sqrt(3.21² − 0.70²) = 3.13` against 3.21, i.e. **~0.08 Å today and material only near 1 Å**. Because
it is uniform it **cancels in every arm-to-arm delta**, which is this project's actual currency.
**This is explicitly not a reason to re-score against the medoid** — the reference is part of the
sealed instrument. Reported as an uncertainty line; nothing changed.

---

## 2. Scientific hypothesis

The charter's premise was that **the objective is the barrier** — justified by S30's finding that
the same circuit family reaches **0.2516 Å** under an ORACLE objective and **3.4330 Å** under the
deployed native-free one, against a classical average's 3.2071. *"The circuit can express good
solutions; the objective does not point at them."*

The sprint's opening move was to read the code rather than the record, and that moved the diagnosis
twice — first one stage later, then to a place no one had looked:

1. **The readout, not the Hamiltonian, caps the stage.** (§6)
2. **The objective is not a function of the target at all.** (§5)

Both are stronger than the premise they replace, and the second dissolves the charter's central
question rather than answering it.

---

## 3. Mathematical mechanism

Four exact results, each verified numerically, together accounting for every negative in the sprint.

**(i) The deployed CVaR free energy is a convex program with a closed-form global minimiser.**
For `F(p) = CVaR_α(E;p) − T·H(p)` (`core/quantum.py:993`), the Rockafellar–Uryasev lower-tail form
is affine in `p` inside a max, so `F` is convex in `p` and concave in `s`; Sion's minimax gives

```
min_p F = max_s { s - T*log sum_i exp( (s - E_i)_+ / (alpha*T) ) }     [1-D concave]
p*_i  proportional to  exp( (s* - E_i)_+ / (alpha*T) )                  [hinged Gibbs]
```

**The whole 2^n optimisation is pinned by one scalar.** Strictly stronger than T1 — *the entropy
term was added to break T1's degeneracy and the answer is still one number.* Derived twice
independently: by Sion duality (lane L) and by KKT on the simplex (lane A), agreeing exactly.

**(ii) The readout's objective is an exact identity, and half of it is free.** For **any** weights
with `Σ_x w_x = 1` — non-negativity **not** required, so this covers the selection, convex and
affine readouts *and* the uniform average in one formula:

```
|| sum_x w_x W_x - t ||^2_F  =  <w, a>  -  (1/2) w' B w
    a_x  = ||W_x - t||^2_F      ORACLE       (per-candidate quality)
    B_xy = ||W_x - W_y||^2_F    NATIVE-FREE  ((1/2) w'Bw = the weighted dispersion)
```

Verified to **1.66e-11** over 126 × 80 draws, simplex and affine. **The only unknown is `a`.
The readout question and the ranking question are the same problem.**

**(iii) The same identity from the other end, with a measured constant.** In the deployed common
frame, `||X(w) − nat||²/n = c2 + 2⟨c,Dw⟩/n + ||Dw||²/n` with `c = X̄ − nat`, `D = W − X̄`. **`c2` is
a property of the candidate set and no sum-to-one readout can touch it**; a readout beats the set
mean only through the cross term. Measured:

> **A sum-to-one readout helps iff `|cross| > ||Dw||²`. ORACLE convex runs at 2.45×. Every
> native-free rule runs at ≈ 1** — the score's orientation toward the native is real and is
> cancelled, to within **0.01 Å²**, by the dispersion that concentration costs. **The best
> native-free rule captures 3.99% of the ORACLE cross term. That is the size of the entire
> native-free readout channel.**

**(iv) No ranker can pass the coherence admission test.** With `Σa = 1`, in pair-distance space
`e_p = μ_p + Σ_m a_m η_{m,p}` — so **`Σa = 1` passes the pool's common mode through with
coefficient exactly one, whatever the weights are.** `coh` is therefore a function of the readout's
*concentration*, **not of the ranker**:

```
uniform mean in pair space        coh = 1.0000  sd 0.0000   <- exactly, as derived
coordinate average (shipped)            0.9780
argmin by the SHIPPED DIS               0.8288
argmin by a RANDOM pool member          0.8244   <- the shipped cost and a coin, 0.0044 apart
ORACLE best member                      0.6708
```

> **THE 0.6931 BAR DOES NOT BELONG IN THIS TABLE, and an earlier draft of this report put it
> there.** S30's `coh` grades a **corrector's residual** — the distogram's own prediction
> error, an *input* to scoring. These arms are **emitted readout errors**, an *output*. Same
> pipeline, same 126 targets, same `mu`, and **production reads 0.6931 in S30's table and
> 0.9780 in this one**, because they are two different errors. **No readout-space admission
> bar has been established**, and against the matched reference (production's own readout at
> 0.9780) *every* arm in this table is lower.
>
> **What survives is the theorem and the within-space comparisons**, which need no imported
> reference: `Σa = 1` passes the common mode with coefficient exactly one, `coh` tracks
> concentration rather than the ranker, and **the shipped cost sits 0.0044 from a random pool
> member.** The claim *"no ranker can ever pass the coherence bar"* is withdrawn and replaced
> by *"coherence does not vary along the ranking axis at all."*

**The entire 43-channel ranking search of S30 was searching a dimension along which the admission
test does not vary**, and the ORACLE ceiling of *all* in-pool ranking barely clears the bar. **The
test can only be passed by leaving the pool's affine hull** — which is exactly what this project's
one confirmed positive, the AMBER relax at k = 30, does.

---

## 4. CVaR-VQE architecture

The deployed stage (`core/pipeline.py:806-877`) is: take the top `2^k` candidates by score, build a
**diagonal** Hamiltonian from their standardised score ranks, minimise a CVaR free energy over a
9-qubit (deployed: **7**-qubit) RY/CNOT statevector circuit, and read the resulting distribution out
in **three** different ways, only one of which emits the answer:

| readout | where | what `p` does | reachable set | capacity |
|---|---|---|---|---|
| **selection** | `core/pipeline.py:795-803` | picks a cell of an argmin arrangement | the distinct pool members | **6.886 bits mean, 6.555 worst** *(the mean of `log2(distinct)` per target; note `log2(118.45)` = 6.888 is the log of the mean, a different quantity by Jensen)* |
| **convex** | `core/pipeline.py:880-895` **SHIPPED, emits the answer** | enters continuously as weights | the **convex hull** | not bit-capped |
| **affine** | `s27/s28_A_amp.py:105-117`, **harness only** | signed amplitudes | **outside** the hull | not bit-capped |

**Two facts about this that the record did not contain.**

**`quantum: bool = False` in `PROD = Config()`** — the 3.2105 Å endpoint is produced with the stage
**off**. And **all three readout classes were conflated in the record**: the charter's *"the circuit
can express good solutions"* rests on numbers from the **affine harness** readout, while a
deployable quantum stage would use the **convex** one. This report distinguishes them everywhere.

---

## 5. Hamiltonian

### 5.1 The finding: it is a target-independent constant

`E = _zrank(pool["sc"][o])`, and `_zrank` (`core/pipeline.py:788-792`) is `rankdata` followed by
standardisation. **The ranks of any 128 distinct values are 1..128, so the standardised vector is a
constant.** Verified on three independent random score vectors:

```
trial 0  first5: [-1.718572 -1.691507 -1.664443 -1.637379 -1.610315]   last3: [1.664443 1.691507 1.718572]
trial 1  first5: [-1.718572 -1.691507 -1.664443 -1.637379 -1.610315]   last3: [1.664443 1.691507 1.718572]
trial 2  first5: [-1.718572 -1.691507 -1.664443 -1.637379 -1.610315]   last3: [1.664443 1.691507 1.718572]
```

> **CORRECTION — the algebra is right and MY VERIFICATION WAS NOT.** Those three trials were fed an
> **already-sorted** vector, so they could only print a ramp — and **random floats never tie**, which
> is precisely the case the caveat is about. **My check was structurally incapable of testing its own
> caveat.** The claim survives only because `sc[o]` *is* sorted in the real pipeline
> (`core/pipeline.py:757` `order = argsort(sc)`, `:863` `o = top[:dim]`, verified non-decreasing on
> all 126 real targets). Measured on the **real** scores instead of synthetic ones:
>
> ```
> targets with >=1 tie in the top-128      125 of 126
> targets whose E is EXACTLY the ramp        1 of 126
> tied positions   mean 9.56  median 8  p90 18  max 34
> max |E - ramp|                           0.0407   (the vector spans +/-1.7186)
> max |E_i - E_j| between any two targets  0.0812
> distinct structures in the top-128       mean 118.45, min 94   <- the source of the ties
> ```
>
> **`118.45` reproduces the independently measured distinct-candidate count exactly** — the ties
> **are** duplicate structures carrying identical scores. So `E` does carry the tie pattern, which
> is a target-specific fact, and **"zero information" overstates a quantity that is now measured.**
> A measured bound beats an absolute; **the conclusion is unchanged.**

Since `p*` is a closed-form function of `(E, α, T)` alone and the circuit is seeded at 0, **both the
closed-form optimum and the circuit's output are one fixed weighting curve per α**, identical across
targets — the ties in a real pool leave a residual of **max deviation 4.07e-02** across the 126, but
the vector is otherwise fixed. Measured: **`H(p*) = 4.9137` bits at α = 1 with standard deviation
`3.1e-04` across the 126
targets**; 6.6392 bits at α = 0.25 with sd 4.5e-3.

> ### The quantum stage is target-independent to within a measured bound.
> **0.0407 max-norm on `E`; 0.0812 between any two targets; downstream `sd(H(p*)) = 1.4e-4`
> bits.** And `(α,T) = VQE_LFO[fold]`, so `p` takes **exactly two values, selected by fold** —
> *one fixed weighting curve per α.*
> The target enters the answer **only** through the readout's own `P` and `W` — never through the
> objective, the Hamiltonian, the CVaR, or the state. The stage answers *"what fixed weight should
> rank `k` receive?"*, which is a **128-number global hyperparameter, not a per-target
> computation.**

**This is strictly stronger than both of the sprint's earlier capacity theorems** — T1 said the
state specifies one integer, R1 said the selection readout's alphabet is 6.886 bits; this says the
state specifies **almost nothing — 0.0407 in max-norm, and two curves selected by fold** — and it
explains every other negative at once: why the objective does
not point at good solutions (*it is not a function of the target*), why α is inert beyond reshaping
a fixed curve, why the circuit's 0.930-bit optimisation gap cancels in the mean, and why the convex
optimum over that family converges to the uniform average.

### 5.2 A non-diagonal Hamiltonian is closed by three independent obstructions

1. **CVaR needs an energy per shot.** Only a diagonal `H` gives every measured bitstring a definite
   eigenvalue; an off-diagonal term cannot be attributed to a single outcome. So
   `H = diag(zrank) − λ·W(block)` **is not a CVaR-VQE and cannot be made into one.**
2. **The forced operator is not even linear in the state.** The identity of §3(ii) forces
   `H[w] = diag(â) − B`, whose VMC local energy is the exact gradient of the readout error — a
   **mean-field operator, quartic in ψ**, with no per-shot eigenvalue.
3. **Dimension counting.** A candidate-index register has Hilbert dimension equal to the candidate
   count, so **any** operator on it is 128×128 and its spectrum is a microsecond `eigh`. **No
   Hamiltonian on a candidate-index register can be classically hard.** The clustering literature
   that does this properly uses **one qubit per data point** — 128, not 7 — and shows no advantage
   at 175 points on hardware.

**And the sign was backwards in the proposal I gave lane A.** The identity's native-free half is
**repulsive**: at fixed candidate quality the readout should *maximise* weighted dispersion. Worse,
the attractive branch **cannot produce a distribution at all** — `B` is conditionally negative
definite, so `w ↦ w'Bw` is concave on the simplex and a concave function on a polytope is minimised
**at a vertex**. Consensus-attraction degenerates to *"pick the single best-scoring candidate"*,
**which is the shipped argmin. A theorem, not a measurement.**

---

## 6. Quantum encoding

**R1, in its final form after two falsifications.** The selection readout is
`v = P·(p/Σp); return argmin_i v_i`, with `P` the pairwise Kabsch CA-RMSD matrix. The state enters
only through a linear map then an argmin, so it matters solely through which cell of the arrangement
`{(P_i − P_j)·p = 0}` it lands in.

- **Falsified once (lane C):** the quantifier *"the entire quantum stage"* is false — the **shipped**
  `average_weighted` emits a continuous convex combination. Measured **1.1144 Å (min 0.073) from the
  nearest pool member**.
- **Falsified twice (lane A):** `argmin(P e_j) = j` fails because the pool contains **duplicate
  structures**, for which `P[i,j] = 0` off-diagonal. Reachable vertices **118.45 / 128 mean, 94
  worst** — *equal to the byte-distinct count on every target*. `filter_pool` already dedups
  (`:766-773`); the readout does not.

> **Final form: the SELECTION readout's alphabet is the number of distinct candidates — 6.886 bits
> mean, 6.555 worst — and it is a scored diagnostic arm (`rmsd_vqe_sel`, `:1173`, registered
> `:1609`), not the answer path.**

**Realised value, ORACLE bits delivered:** `p_θ` medoid 1.805, `p*` medoid 1.876, score argmin
1.693, uniform medoid 1.659. The stage moves the selection on **77.8%** of targets and delivers
**0.112 bits over the plain argmin — 1.6% of the alphabet.**

**Index redesign is closed.** Re-labelling cannot raise capacity above `k` bits, and it does not
help in practice either: **Gray coding is a proven no-op** (identical prefix partition on 126/126,
asserted in code; −0.0014 Å at 0.03× MDE); the best index map raises an ORACLE ceiling by −0.0655 Å
and **lowers** the deployable value; **the best map for a deployable partial readout is the random
permutation**; and no index bit carries more than **0.0658 bits** about candidate quality (deployed
map 0.0343, permutation control 0.005). **The best available map doubles a quantity that is 0.03
bits.** The mechanism is a conflict by construction: differentiating the cells makes the best cell
better and the chosen cell worse, because *a coherent cell is a concentrated set* — §3(iii)
arriving at the index.

---

## 7. Objective geometry

**The circuit does not reach the closed form, and the readout cannot tell.** `KL(p_θ‖p*) = 0.930
bits mean, TV 0.378`; `run_cvar_vqe` is **strictly worse in 12/12** verification cells and in
**36/36** targets on the real instrument. **The gap is expressivity, not optimisation**: 80 → 2000
Adam iterations move it by nothing (0.025764 → 0.025944); only depth and restarts close it, and at
layers = 12 / 400 iters / 8 restarts (~60 s) it is still short — against **microseconds** for the
closed form. *21 parameters cannot cover a 127-dimensional simplex.*

**And yet the two readouts differ by −0.0044 Å at 0.088× MDE** while differing **per-target by
0.102 Å in absolute value**. *The optimisation gap is real and cancels in the mean.*

**The inductive bias, characterised.** The reachable set is an **MPS Born machine of bond dimension
exactly `2^layers = 8`** (measured Schmidt ranks 2, 4, 8, 8, 4, 2). Against the objective's own hinge
shape, `R²` is **0.0071** for a random θ, **0.0538** for the best of 400 draws, and **0.467** for the
optimised `p_θ` — *a ceiling at about half the right shape, not a bias toward anything structural.*

**α is inert on most of the benchmark.** The deployed table is
`VQE_LFO = {0:(1.0,0.3), 1:(0.25,0.3), 2:(0.25,0.3), 3:(1.0,0.3), 4:(1.0,0.3)}` — **T = 0.3 on every
fold**, and **α = 1 on three of the five pinned folds, where `CVaR_1 = ⟨E,p⟩` and CVaR is inactive by
construction.**

**A shipped docstring corrected.** `core/quantum.py:997-1002` asserted the entropy collapse *"is a
property of CVaR, not of the optimiser."* Measured at T = 0 on the real pool: **0.258 bits at α = 1
(78 targets), 3.596 bits at α = 0.25 (48 targets)**. At α = 1 the minimiser is the unique argmin
vertex, so the collapse **is** a property of CVaR; at α < 1 the argmin set is `{p : p_x₀ ≥ α}`, a
face of **positive volume**, so the objective does not determine `p` at all and **the entropy
belongs to the optimiser** — the opposite of the sentence. *"For ANY alpha" is false by theorem.*
Corrected in place with the original quoted.

---

## 8. Free-energy work

**Closed by derivation. No thermodynamics was computed**, which is the charter's §7A hierarchy
working as intended.

**The gating question — what is `S` a function of?** `S(x) = Φ_{β,B}[U_seq](x)`. Its arguments are
the candidate, the sequence (only through atom typing in `U`), and `β` plus the basin map, which are
universal. **There is no fourth argument**, so `S` is an *operator* on (sequence, pool), not a third
source.

**Lemma B1, verified rather than assumed:** the shipped `amber14/protein.ff14SB.xml +
implicit/gbn2.xml` is **reflection-invariant** — ff14SB's **1340 `PeriodicTorsionForce` phases sit at
distance 0.000e+00 from {0, π}**, and there is no CMAP.

**Corollary B1:** `F`, `E`, `S` and every `dF/dT` are rotation-, translation- **and
reflection**-invariant single-structure observables ⇒ by G1 each is a **function of the distance
map**, inside the class S30 closed by theorem and measured empty on 43 channels. The enthalpic half
was already closed by measurement; **the entropic half is closed by the same theorem that closed
contact topology and Rg, because `S` is achiral for exactly the reason `E` is.**

**Corollary B1′, which nobody had connected:** an ANM/GNM Hessian is built from **pairwise
distances**, so its spectrum, its log-determinant — *which is the harmonic configurational entropy
`S_harm = const − ½ ln det′ H`* — and every spectral-graph observable derived from it are
distance-map functions. **This closes the elastic-network, normal-mode and energy-landscape-curvature
family of charter §14 without computing one of them.** Basin populations, ensemble reweighting and
temperature-dependent ranking are monotone transforms of single-structure free energies and close
with it; local configurational entropy closes twice.

**A registered falsifier fired and its matched control inverted the reading.** F1's 1e-6 threshold
fired at **3.673e-06** — but a pure **rotation**, equally analytically exact, gives **5.377e-06,
larger** (ratio 0.68×), and the only force groups that could carry chirality are exact to
**2.0e-16**; the residual lives entirely in two terms invariant under *any* isometry. **A registered
threshold set at an absolute value, on a quantity only a matched control can read, is a mis-set
threshold — not a falsification.** Recorded with the original standing; F1 is never quoted without
the rotation control in the same sentence.

**Left open on price, not theory:** a genuine **multi-structure** barrier `F‡(x_a → x_b)` depends on
the *path*, not on `D(x_a)` and `D(x_b)`, so **G1 does not bind it.** Not closed, not run — 10³–10⁵
trajectories per target against a one-AMBER-process box.

---

## 9. Torsion work

**Closed, and the reason is worth more than the closure.**

**The decomposition theorem:** for an ideal-geometry backbone, the point reflection acts on torsions
**exactly** as `(φ,ψ) → (−φ,−ψ)`, so `T_even` is reflection-invariant and by G1 a distance-map
function. **All of a torsion channel's escape from G1 lives in `T_odd`.** Crossed with separability,
the only cell no theorem closes is **chiral AND non-separable**.

`LEG_torsion` is **30.0% odd** — genuinely chiral, so the registered G1 check does not fire. But the
skill sits in **different halves on the two problems**:

| | EVEN half (G1-closed) | ODD half (G1-escaping) |
|---|---|---|
| coarse triage (meter contrast) | +0.1220, 1.35× MDE | **+0.2698, 2.78× MDE** |
| in-pool selection, 500 band | **+0.2428, 2.51× MDE** | +0.0183, 0.28× MDE |
| in-pool selection, top-75 band | +0.0288, 0.40× | +0.0423, 0.66× |

> **The chiral escape from G1 is real, measurable and not degenerate — and it only works on a
> problem the pipeline does not have.** On the problem it does have, the surviving skill is in the
> **G1-closed** half, and even that dies in band.

- **G2 fired as registered:** in-band `LEG_torsion` is **0.65× MDE**, below the 0.7× bar — where the
  **shipped** cost `DIS` is itself **0.83×, NOT MEASURED in its own band**.
- **G5 fired and is the sharpest number:** a **constant α-helix** at (−57, −47) — plausible, never
  uniform-on-the-torus — **beats every torsion channel on both bands**: ρ +0.3413 against
  `LEG_torsion`'s +0.1676 on the 500 band (paired gap **−0.1722 at 3.17× MDE**), and **−0.0901 at
  1.81× MDE, 5/5 folds** in band. **The Ramachandran channel's in-pool skill is a constant-prior
  effect** — *"this benchmark's peptides are helical"* — and one constant does it better than the
  fold-conditioned channel.
- **M6, the only open cell, is empty:** `XTWIST` (chiral and non-separable by construction) scores
  −0.0104 / +0.0320 and **its own achiral twin beats it** (+0.1124). Scoped honestly: **empty at
  9–16 residues**, which is where S30 left chiral functionals too.

**Mechanism, measured:** `LEG_torsion` of the `RAND_SIGNED` rung sits at **z = +2.363** in the pool's
own torsion distribution, with 70.6% of targets above the pool's 95th percentile, against PROD at
+0.844 and ORACLE `circ_best` at +1.129. **The channel's dynamic range is spent outside the pool.**

---

## 10. Sparse / readout work

### 10.1 The readout class question, answered — with its information cost attached

Built chain, n = 126, **ORACLE / NOT DEPLOYABLE**, from `s29/results/s29_O_chain_rows*.jsonl`
(nothing reprojected):

```
production (uniform top-75 average)                3.2105
bestm128   ORACLE prefix-average over the 128      2.9027    T1's entire reach
best1_128  ORACLE argmin over 128                  2.1435    7 BITS
hull_128   ORACLE CONVEX weights over the SAME 128 1.8538    128 REAL NUMBERS
best1_pool ORACLE argmin over 500                  1.7078    9 BITS
```

**Combining reaches further than naming at a fixed candidate set. Naming wins per bit.** Both
sentences are needed and only the pair is honest — and the second is the one that governs, because
the project's scarce resource is information, not expressiveness.

> **Every lever comparison in this report carries its information cost in the same sentence as its
> Ångströms.** Reading `hull_128 = 1.8538` against `best1_128 = 2.1435` as *"the readout class is
> the larger lever"* compares **7 bits against 128 free reals**, and that is the error S30 §9.4
> already named as this project's characteristic one. It is a **larger ceiling at unpriced cost**,
> not a larger lever.

**And S30-L11's conclusion survives its own mismatched support.** S30 claimed *"a plain argmin
dominates at every bit budget"* and supported it with 2-of-**75** against 1-of-**128**. The support
was mismatched — but **S30 computed the correct reference curve in the same file**
(`s30_Q_sparse.json :: argmin_ref`: 1.9383 at 8 bits, **1.7108 at 9**), and the sparse arm prices
out at ~**9.1 bits**, where naming gives 1.7108 — **0.36 Å better at matched information. Nothing
in the ladder beats the argmin at a matched budget.** *An earlier draft of this report stated that
the set-matched ladder "inverts" S30-L11. That was a withdrawal that was itself wrong; it is
retracted here.*

### 10.2 The matched-search-size correction, which is a finding

`2-of-128` is a per-target minimum over **C(128,2) = 8128** ORACLE supports; `argmin over 128` is a
minimum over **128**. Comparing them directly is a **64× oracle-search mismatch**. Measured at
matched `K` (**ORACLE / NOT DEPLOYABLE**, CA cloud):

```
ORACLE argmin over 128 singletons   (K=128)     2.1458        --
2-of-128 uniform, 128 RANDOM pairs  (K=128)     2.2229   +0.0771  WORSE  1.11x MDE, 5/5, 44W/82L
2-of-128 uniform, EXHAUSTIVE        (K=8128)    1.9070   -0.2388  better (3.81x MDE)
```

**And it decomposes exactly:**

```
family LEVEL   singles 3.5847  pairs 3.3469   -0.2378   <- error cancellation IS real
dispersion loss (pairing compresses the family)        +0.3149
net at matched K                                       +0.0771  pairs WORSE
```

> **Error cancellation is real and worth −0.238 Å in level — and is then more than consumed by
> averaging compressing the family's dispersion.** This is §11's prefix law, **second independent
> instance, opposite operator**: one lane decorrelated a family and its minimum **grew**; another
> contracted one and its minimum **shrank**.

### 10.3 The native-free channel, measured exactly

§3(iii)'s identity gives the closure. A sum-to-one readout beats the set mean **only** through the
cross term, and:

```
                        cross    ||Dw||^2   emitted rmsd (CA cloud)
unif_prefix_m128       +0.0000    0.0000        3.0435
soft_T4.00             -0.0353    0.0257        3.0391
soft_T1.00             -0.2025    0.3097        3.0442
typicality T1.00       +0.1844    0.2708        3.0957   <- WRONG SIGN
ORACLE convex128      -10.3193    4.2054        1.7977   <- 2.45x
```

> **ORACLE convex runs at 2.45×; every native-free rule runs at ≈ 1.** The score's orientation
> toward the native is **real** and is cancelled, to within **0.01 Å²**, by the dispersion that
> concentration costs. Sharpening buys cross monotonically and pays dispersion faster every time.
> **The best native-free rule captures 3.99% of the ORACLE convex cross term — that number is the
> size of the entire native-free readout channel.**

### 10.4 The constrained-affine middle has no selection rule

- **No knee by ceiling.** The ridge path from uniform to unconstrained affine is smooth and
  monotone (`||w||₁` 1.01 → 28.3 gives ORACLE 2.7318 → 0.0000), and at `||w||₁ = 2.09` the ORACLE
  already beats the argmin over the whole 500. **There is no ORACLE answer to "what is the right
  constraint".**
- **No answer by fit either.** ORACLE weights are not identifiable: out-of-fold **R² = 0.021**
  against a 5% bar, and **applied it is +0.0309 Å worse**.

**So a constrained-affine readout cannot be selected by ceiling and cannot be selected by fit.**

**And sparsity is an output, not a design choice:** the exact objective's optimum is **already
sparse — 6.54 of 128 members with no sparsity penalty imposed.** Sparse-readout work that does not
improve `â` is spending effort on the half of the objective that is already exact.

### 10.5 The affine readout's ceiling is vacuous by rank

`rank(aff{W_x}) = 32.9 ≥ 3n−3`, so the affine hull of 128 candidates **spans the residual space**
and the ORACLE affine ceiling is **exactly 0.0000 Å on all 126 targets** — 127 weights against ~33
dimensions, an over-parameterised interpolator (median `ess` **1.13 of 128**, `neg_mass` 3.38).

> **Therefore "0.2516 Å under an ORACLE objective" is a statement about a REGULARISER, not a class
> ceiling** — and the simplex constraint is precisely the regulariser the uniform average enjoys for
> free. That *derives* "expressivity without an aligned objective is harmful" instead of observing
> it: the uniform average is the heavily regularised special case of the same affine readout, and
> under the deployed objective it **beats** the unconstrained optimisation.

---

## 11. Candidate-index information

**Closed. Re-labelling cannot buy information, and in practice it does not buy anything else
either.**

- **Gray coding is a proven no-op** — the j-bit prefix partition is **identical to binary's on
  126/126 targets** (asserted in code); measured ceiling difference **−0.0014 Å at 0.03× MDE**.
- **Structure-aware indexing helps the ORACLE and hurts the deployable**, by construction: it
  raises the ORACLE value of a partial measurement (2.4216 → 2.2712) and **lowers** its deployable
  value (3.2369 → 3.3721). *Differentiating the cells makes the best cell better and the chosen
  cell worse, because a coherent cell is a concentrated set* — §3(iii) arriving at the index.
- **The best index map for a deployable partial readout is the random permutation** — the same
  phenomenon §12 finds at the prefix length.
- **No index bit carries meaningful information about candidate quality.** Mutual information
  between an index bit and "in the ORACLE-best decile": deployed map's best bit **0.0343 bits**,
  best map **0.0658**, permutation control 0.005. **The best available map doubles a quantity that
  is 0.03 bits.**

**The 128 → 512 register widening** survives its circularity gate at **−1.1811 Å cloud / −1.1762
chain** on a filter-independent tail (**ORACLE / NOT DEPLOYABLE**, p = 0.0001), **62% of the
published −1.9004** which is retired as a near-tautology. But **it only pays if a selector exists**,
and three facts held by three different lanes compose into a prediction that it would not:

1. **On the tail, the score is worse than uninformative:** FAIL18's ORACLE-best member sits below
   rank 128 at frequency **1.000** against an uninformative-score null of **0.744**.
2. **Widening degrades the very selector it requires:** `Var ~ 16/D` is the predicted global-cost
   barren plateau; 128 → 512 costs **4×** in gradient variance.
3. **`operator-consumes-set-mean`:** admitting 384 more candidates to a set the selector cannot
   order **raises the set mean**, consumed at coefficient 1.16.

> **PREDICTION, labelled as one and not run: a DEPLOYABLE 128 → 512 arm should come out WORSE than
> production, not merely flat.** Caveat carried because it cuts the other way: at n = 7–9 we are
> **not yet gradient-limited**, so fact 2 bounds the direction without quantifying the Ångströms.

---

## 12. FAIL18 analysis

**`FAIL18` is not a stratum. It is an identity, and this sprint proved it.**

A pre-registered definition-matched stratum `defn18` — the 18 targets with the largest widening gain
from the top-75 — turned out to be **`defn18 ∩ FAIL18 = 18 of 18. The same targets, effect identical
to four decimals.** Choosing the 18 targets whose best candidate hides below rank 75, then reporting
that their best candidate also hides below rank 128, **is one statement and not two.**

**The rank diagnostic confirms it independently and says something sharper.** Under a score with no
within-pool skill the ORACLE best member's rank is uniform on 1..500, so `P(rank > 128) = 0.744`:

```
other 108           0.417     the score genuinely pulls the answer into the window
worst18_poolmean    0.667     near the uninformative null
worst18_bestpool    0.722     near the uninformative null
FAIL18              1.000     ABOVE the null -- p = 0.744^18 = 0.005 by chance
```

> **`FAIL18` is the only stratum *worse than an uninformative score*, which cannot happen by
> sampling.** That is the tail's mechanism: **not "hard targets hide their answer deeper" but "the
> score has no within-pool skill on hard targets, so the best member lands roughly uniformly and a
> 128-window misses it ~74% of the time."**

**And the excess is mostly definitional.** Regressing per-target widening effect on `best1(128)`
gives slope −0.4356, R² 0.5564; the residuals are **−0.577 for `FAIL18` against −0.196 for the clean
filter-independent tail**. **Two thirds of `FAIL18`'s excess is the stratum definition.**

**Consequences enforced throughout this report:** filter-independent tails are primary (`worst18` by
pool mean and by best-in-pool), `FAIL18` appears only as a labelled diagnostic, and **Appendix B is
the key** — five different strata are called "the worst 18" somewhere in this project's record and
they differ by up to **0.81 Å on the same-named quantity**.

**The reference is not the explanation** (§1.3): the tail's targets have, if anything, **narrower**
deposited NMR ensembles, and the tail-minus-rest difference in ensemble spread is **0.20× MDE — a
null**. **FAIL18 is real failure.**

**Where the tail's prize is, ORACLE / NOT DEPLOYABLE:** the along-`mu` correction (§20.1) is **4×
concentrated on the tail** — FAIL18 6.0195 → 3.7392, **−2.2803 Å**, against −0.5651 on the other 108.
The tail is where the endpoint's mean is decided and where the only priced route would pay most.

---

## 13. Controls and nulls

Charter item 25. **Several results in this sprint were killed by their own lane's control**, which
is the point of the table.

| control | what it was matched to | what it ruled out |
|---|---|---|
| **Shuffled-`B`** (random relabelling of the pairwise matrix) | quality-blind dispersion maximisation, which beat the shipped argmin by −0.2621 Å at 2.74× MDE | **Killed the reason for its own lane's only positive.** Costs 0.0250 Å at **0.54× MDE — NOT A RESULT** — so the mechanism is *"spread the weights over many candidates"*, **not** *"spread them along the real geometry"*. No part of the gain may be attributed to the pairwise structure |
| **Matched random subsets at fixed `K`** | the per-target ORACLE prefix length | **Destroyed `bestm128` as a lead.** Random subsets reach **149%** of the prefix family's gain: the prefix axis is *worse than an arbitrary 7-bit index* |
| **Matched search size `K`** | 2-of-128 against argmin-over-128 | **Inverted the sign.** A minimum over C(128,2) = 8128 supports read against one over 128; at matched `K` the 0.076 Å gain is a **0.077 Å penalty** |
| **`A2` = the same operator written twice** | every built-chain arm-to-arm comparison | **Recalibrated the whole instrument.** Clouds agree to 3.6e-14, chain means at 0.03× MDE — but per-target \|Δ\| is 0.0134 mean / 0.0329 p90 / **0.2285 max**. **The 0.0107 Å floor is a MEAN floor** |
| **Pure rotation** (analytically exact, like reflection) | a reflection-invariance test whose 1e-6 threshold had fired at 3.673e-06 | **Inverted a registered falsifier.** Rotation gives **5.377e-06, larger** — the threshold was mis-set at an absolute value only a matched control could read |
| **Constant α-helix** (plausible, never uniform-on-the-torus) | every torsion channel | **Beat all of them on both bands** (−0.1722 at 3.17× MDE) — the Ramachandran channel's in-pool skill is a **constant-prior effect** |
| **Achiral twin** `\|X\|` | chiral functionals | Whatever the chiral channel does, its reflection-invariant shadow already does |
| **Deposited-ensemble spread** | the tail's difficulty | **Ruled out the reference as the explanation** — 0.20× MDE, and the tail's ensembles are *narrower* |
| **Random-feature control** | the predictability split | `y_perp` R² **+0.9403** against −0.0026; `y_along` −0.0051 against −0.0039 |
| **Magnitude-matched shrinkage** | the orthogonal-correction arm | **The perp arm is +0.6169 Å worse than its own magnitude-matched control** (3.25× MDE, 9W/117L) — magnitude is not the explanation |
| **Norm-matched shrinkage** | the projection operation itself | 0.26× MDE — projecting is indistinguishable from making the corrector smaller |
| **Uniform-weight ablation** (`:1116`, already in the code) | the p-weighted convex readout | The shipped ablation, used rather than rebuilt |
| **Shuffled-feature / shuffled-label** | the `â` combination arms | 0.13× and 0.71× MDE — NOT MEASURED |
| **Random permutation index map** | structure-aware indexing | **It is the best map** for a deployable partial readout |
| **Random-18 stratum null** (20,000 draws) | tail claims | Rules out *"an ordinary 18-subset"* — and **does not** rule out *"a stratum defined by the quantity measured on it"*, which is why `FAIL18` needed §12 |
| **Vertex feasibility assert** | the ORACLE convex solver | Caught 4 of 126 targets where both iterative solvers returned a point worse than a feasible vertex; **119W/4L → 119W/0L**, the four losses *were* the four failures |

**Two controls that changed a lane's own conclusion rather than confirming it** deserve naming:
the **shuffled-`B`** control (which demolished the reason for its lane's only positive) and the
**rotation** control (which inverted a registered falsifier that had fired). Both were run by the
lane that stood to lose from them.

---

## 14. The 126-target endpoint

**Production is 3.2105 Å and nothing this sprint moved it.** No arm was deployed, and no arm earned
deployment.

### 14.1 The quantum stage, measured end-to-end for the first time

All six arms **built and projected in one process from one stored cache**, with arm A re-projected
to **3.2105 Å against the canonical 3.2105 Å, worst per-target deviation 0.0000 Å** — so the
comparisons are internally valid under §1.2's rules. Built chain ± SE; **CA point cloud in brackets,
a different object**:

```
A   quantum OFF (production)      3.2105 +/- 0.1540   [3.0483]
A2  the same operator, recoded    3.2108 +/- 0.1542   [3.0483]   <- the implementation-noise null
B   selection readout, VQE p      3.3117 +/- 0.1648   [3.3135]
C   selection readout, p*         3.2838 +/- 0.1586   [3.2852]
D   convex readout, VQE p         3.2281 +/- 0.1558   [3.0661]   <- the SHIPPED quantum arm
E   convex readout, p*            3.2169 +/- 0.1526   [3.0605]
F   convex readout, uniform p     3.2415 +/- 0.1528   [3.0532]
```

| contrast | effect | × MDE | verdict |
|---|---|---|---|
| **`P1` = E − D** — solving the objective *exactly* instead of with the circuit, shipped readout | **−0.0112** | **0.19×** | **NULL** |
| `P2` = C − B — the same, selection readout | −0.0280 | 0.29× | NULL |
| **`S3` = D − A** — the built-chain cost of the **shipped quantum stage** | **+0.0175** | **0.35×** | **NULL** |
| B − A — selection readout against production, **CA cloud** | +0.2652 | 2.15× | **MEASURED, worse** |
| B − A — the same contrast, **built chain** | — | 0.93× | NOT MEASURED |

**`S3` decomposes exactly:** `+0.0175 = +0.0003` (code path) `+ 0.0307` (prefix widening 75 → 128)
`− 0.0134` (weights). **The largest component is the widening that exists only to fill the `2^7`
register** (`core/pipeline.py:757`) and is pointless when `quantum = False`, which is production.
**NOT MEASURED at 0.42×** — *and its fold CI [+0.0028, +0.0557] excludes zero, which is the exact
shape `stats_lib._verdict` was hardened against; the MDE gate binds.*

### 14.2 The null is not "no effect"

```
                  |d| mean  median    p90     MAX    exactly 0   > 0.0107
P1  E - D          0.1432  0.0494  0.3706  1.1457    0 / 126    100 / 126
P2  C - B          0.1802  0.0132  0.6288  1.8766   60 / 126     64 / 126
X1  A2 - A (null)  0.0134  0.0026  0.0329  0.2285    0 / 126     28 / 126
```

**`P1`'s per-target spread is 10.7× the implementation-noise null, and 74 of 126 targets move by
more than that null's own p90** — while the mean is −0.0112 at 0.19×. The honest sentence is
***"solving the objective exactly reshuffles the answer everywhere and buys nothing."*** The two
selection readouts choose a **different candidate on 66 of 126 targets** (54 of 78 at α = 1) and are
identically tied on the other 60.

### 14.3 The caveat that limits what §14.1 can mean

**Both substituted arms are target-independent rank curves (§5).** `P1` therefore compares two
points inside a family **none of whose members knows which target it is** — so it could not have
found an aggregate effect even if one existed in the objective. **This experiment does not test
whether a better objective would help, and must not be recorded as evidence either way.** What it
establishes is what the charter needed: **the optimisation quality of the deployed CVaR-VQE is
invisible at the endpoint**, with certificates at the deployed settings (duality gap 3.56e-09, KKT
2.67e-15, **circuit strictly worse 126/126**, `p*` at 0.0012 s against 0.0985 s).

### 14.4 The terminal-operator arms — the answer is keep the uniform top-75 average

n = 126, built chain, comparator = **production re-projected in the same job** (3.2126, not the
canonical 3.2105 — the per-target deviation from the record reaches **0.3256 Å**, above lane P's
0.2285 same-operator maximum, which is exactly why the comparator must be re-projected rather than
quoted):

| arm | mean | vs production | × MDE | verdict |
|---|---|---|---|---|
| `MED` — the consensus medoid, a real deposited member | 3.2814 | +0.0688 | 0.81× | NOT MEASURED |
| `AVG_RG` — uniform Rg restore | 3.2602 | +0.0475 | 1.11× | **MEASURED WORSE** |
| `AVG_SEP` — per-separation restore + MDS re-embed | 3.6736 | **+0.4609** | **2.34×** | **MEASURED MUCH WORSE** |

**`AVG_SEP` was the sprint's only live deployable candidate and it is dead.** Its full paired
distribution: mean +0.4609, **median +0.2045**, sd 0.788, **41W/85L**, p10 −0.2283, **p90 +1.7176**,
best −0.7328, **worst +2.9878**. *The early 0.4–0.6 Å wins seen at ~50/126 were the left tail of a
distribution whose right tail is three times longer* — **and no mean was computed before 126/126.**

**Every gate failed in the registered direction** (dispersion→MED +0.0591, move→MED +0.0600,
move→AVG_SEP +0.1412, dispersion→AVG_SEP +0.1457), and the ORACLE per-target `min(AVG, MED)` is
worth only **−0.0782 Å** with a split-half transfer of −0.0323 whose CI spans zero. **Choosing
per-target between the two best operators, with the native in hand, is worth 0.078 Å.**

**The lane's pre-registration predicted its own endpoint.** §2, written before any number existed,
derived `d_MED² ≈ B² + s_b²` and therefore a chain cost of `+0.2339 − 0.1622 ≈ +0.0717`. **Measured
+0.0688.**

**And the mechanism I proposed is refuted, monotonically.** The registered falsifier required the
high-minus-low dispersion contrast of `chain(MED) − chain(AVG)` to be **negative**; it is **+0.0990**,
rising across dispersion tertiles (+0.0049 / +0.0604 / +0.1410). Both halves of the trade-off scale
with dispersion — but **the medoid's own cost scales faster, so there is no operating point where it
turns over.** *That is why no gate can work, not merely why these gates did not.*

**And the reading that makes this a result rather than only a retraction.** `AVG_RG` was demoted to
a **registered negative control** — the correction matched to the *withdrawn* uniform-scale
mechanism. `AVG_SEP` was the correction matched to the *corrected* per-separation mechanism, and was
expected to work precisely where `AVG_RG` would not. **Both fail, and `AVG_SEP` fails hard.**

> **The per-separation distortion is real and measured — and correcting it makes things worse. So
> the distortion is a SYMPTOM of averaging, not the MECHANISM of its error**, which closes the
> entire *"repair the average's shape"* family rather than just this arm.

**Two process notes that belong with the number.** First, `AVG_SEP` is the **only arm carrying none
of the four registered geometry secondaries** (`bond_cloud / bond_chain / rg_cloud / rg_chain`) —
and it is the only arm that builds a *new* structure by rescaling pair distances and re-embedding by
classical MDS, where invalid geometry is a live possibility rather than a formality. Its
`move_AVG_SEP` is **0.792 Å mean, 1.746 max**: stage-3b has to drag its output a long way to make it
a chain at all, which is consistent with the refutation and **would have been visible earlier**.
Second, the interim *"0.4–0.6 Å better"* that circulated mid-sprint was **the magnitude of the
deficit with the sign inverted** — `stats_lib.compare` is lower-is-better, and the final effect is
+0.4609 with fold CI [0.378, 0.565], *exactly that band*. **The lane published no aggregate before
126/126; the inverted reading was the coordinator's summary of partial rows**, which is precisely
what a partial reading gets wrong.

> **A correction to my corrected mechanism, and it is mine to carry.** After the 25.8% contraction
> figure was withdrawn I offered *"averaging contracts short-range and expands long-range, crossing
> unity near \|i−j\| = 8"*. **That is wrong as stated.** The 75 members' own profile also rises (to
> 1.448 at s = 14), and **the average is below the members at every single separation.** Averaging
> **contracts everywhere and contracts least at long range**, against a pool that is already long
> there — *the crossing of 1.00 near s = 8 is an artefact of comparing the average to the NATIVE
> instead of to the members that produced it.* The band falsifier fired too: the medoid is **less**
> flat at long range (0.2918 against 0.2824). What the projection actually does is **repair the
> short band** (0.1609 → 0.1033) and leave the long band alone.

### 14.5 What distinguishes the tail, on filter-independent strata

`FAIL18` used only as a cross-check (§12). CA point cloud throughout.

- **Useful candidates are present and badly ranked.** The pool's best member sits at **rank 391 of
  500** under the shipped score on the outcome tail against **134** elsewhere (286 against 151 on
  the pool-mean tail).
- **Whether the tail is a filter failure or a readout failure is NOT MEASURED** — and the attempt
  to say otherwise is instructive. A table of *means* showed filter loss growing 2.4–3.4× on the
  tails against readout loss at 1.9×, which reads as *"disproportionately a filter failure."* **The
  statistic that claim requires** — per target, `(filter loss) − (readout loss)`, tail minus rest,
  fold-clustered — **does not clear MDE on either filter-independent stratum:**

  ```
  T_POOL  (filter-independent)      +0.1925   0.16x MDE   fold CI [-0.401,+1.256]   2/5 folds
  T_BEST  (filter-independent)      +0.4298   0.38x MDE   fold CI [-0.271,+1.270]   2/5 folds
  T_CHAIN (defined by the OUTCOME)  +0.8781   0.84x MDE   fold CI [+0.312,+1.559]   4/5 folds
  FAIL18  (the filter's OWN zero-recall set)
                                    +1.6165   1.24x MDE   fold CI [+1.026,+2.933]   4/4 folds
  ```

  > **The effect size rises monotonically with how circular the stratum is — 0.16×, 0.38×, 0.84×,
  > 1.24×. That gradient *is* the signature of the stratum's definition doing the work**, and it is
  > S30-L2's failure mode exactly.

  **The defensible statement:** on a filter-independent tail, **the filter's loss and the readout's
  loss are both roughly doubled, and which is hurt more is NOT MEASURED.** Likewise
  `n_top75_under3 = 0.0` is outcome-defined; the clean version is **5.7 against 33.9** on the
  pool-mean tail. *(An earlier draft of this report carried the stronger claim, and the lane
  retracted it with its own registered control before anything was published.)*
- **The tail's pool is COHERENTLY wrong, not diversely wrong.** `S/B` falls 0.78 → 0.42 while
  `n_distinct` is unchanged (70.3 against 69.0). ***The pools are not smaller or less varied — they
  are displaced together.***
- Independent confirmation of the long-range concentration: chain `|ratio−1|` at `s ≥ 7` is **0.915**
  on the tail against 0.177 elsewhere (**5.2×**), against 0.200/0.087 (2.3×) at `s < 7`.
- **The readout is not failing in a way that favours the medoid**: `chain(MED) − chain(AVG)` is
  +0.087 / +0.203 / +0.173 on the three tails against +0.066 / +0.046 / +0.051 elsewhere.

**And the consequence that retires the readout direction by derivation rather than by exhaustion.**
Every readout tested in this sprint is an affine combination with `Σa = 1`, and §3(iv)'s theorem says
such a readout passes the common mode through **with coefficient exactly one, whatever the weights**.

> **An operator cannot remove the component that carries the prize.** `AVG_SEP` was the attempt to
> leave that hull and measured **0.045 Å from it**. So the readout is not where the tail is lost —
> and this lane's four refutations are **four instances of one theorem**, which is a better reason to
> stop spending on the readout than its ORACLE ceiling of 0.078 Å, *because it is a derivation
> rather than an exhausted search.*

> **The coherently-wrong finding above is the sprint's strongest cross-lane agreement, and neither
> half is a synthesis.** Lane F measured that the tail's pools are *displaced together*; lane E measured,
> independently and on a different object, that **the entire recoverable prize lies along the pool's
> common mode** (−0.8102 Å) while the orthogonal component is harmful. **"The tail's pools are
> displaced together" and "the prize is the common mode" are the same fact** — and it is why the
> along-`mu` correction pays **4× on the tail**: *the tail is where the common mode is large.*

### 14.6 The prefix-length arms on the built chain — S29-L30's cross-basis quotation, repaired

S29-L30's transfer arms were quoted on the **CA point cloud**; the cloud→chain price on this rung is
**+0.14 to +0.16 Å**, so they were never endpoint statements. All six arms below are
**projected in one process from the same stored clouds**, every variant selected on the cloud so the
two oracles see the same thing, comparator **M75 = 3.2126 re-projected in that same job**:

| arm | chain | vs M75 | × MDE | folds | W/L | verdict |
|---|---|---|---|---|---|---|
| `PREFIX` — per-target ORACLE `m` (`bestm128`) | 2.9027 | **−0.3100** | **3.17×** | 5/5 | 118/8 | **MEASURED** — ORACLE |
| `M_GLOBAL` — one ORACLE global `m* = 72` | 3.2083 | −0.0044 | 0.19× | 3/5 | 69/57 | **NULL** — still ORACLE |
| `M_LFO` — leave-fold-out `m`, the deployable arm | 3.2201 | **+0.0075** | 0.22× | 2/5 | **63/63** | **NULL, and the wrong sign** |
| `RANDOM` — matched random-subset family, 2 draws | 2.7770 | −0.4356 | — | 5/5 | — | **ORACLE, and it BEATS the prefix** |

**The deployable arm is a coin flip on the endpoint: 63W/63L, +0.0075 Å, 2/5 folds.** On the cloud
the same transfer read −0.0039 (**1.2%** of the ORACLE gain, CI spanning zero); on the built chain it
is **−2.4%**. ***Both are nulls, and the sign difference between them is noise inside the null — not
an inversion.***

> **What the repair actually buys, stated so it is not oversold.** S29-L30's transfer arms were
> **cloud** numbers quoted as endpoint statements, which was unlicensed *at the time*. Measured on
> the endpoint, **the conclusion is unchanged** — *"the transferable part of the prefix axis is
> zero"* is now true on the basis it is stated on, in one job. **And lane F found the reason it
> carries, which generalises past these four arms:** the cloud→chain price is **flat across
> prefixes** — M75 **+0.1643**, PREFIX **+0.1422**, global-`m` **+0.1618**, LFO-`m` **+0.1639**.
> *The projection charges every prefix nearly the same, so conclusions about `m` transfer from cloud
> to chain by construction.* **That is a licence for the `m` axis specifically, and it is not a
> licence anywhere else** — §14.4's `AVG_SEP` pays a completely different price because it
> re-embeds.

**And S31-L11's refutation carries to the endpoint intact.** A matched random-subset family — same
top-128, same operator, same `K`, same size distribution, but an arbitrary subset instead of the
score-ordered prefix — reaches **141% of the prefix gain on the built chain** (149% on the cloud),
beating `PREFIX` by **−0.1368 (1.13× MDE, 5/5 folds, 90W/36L)** on draw 0 and −0.1146 (0.93×, 5/5,
83W/43L) on draw 1, with a **draw-to-draw sd of 0.0157** — 11% of the effect, so not a lucky draw.
***Draw 1 is at 0.93× and is NOT MEASURED on its own; the claim rests on the mean of the draw
distribution, not on the better draw*** (contract rule 10, and the lane flagged it rather than
quoting draw 0). The lane's registered bar (*"no part of it is deployable"*) **fires on the chain as
it fired on the cloud.**

> **The mechanism is an order statistic, and the lane measured it rather than asserting it.** Prefix
> variants are **nested** — `curve[m]` and `curve[m+1]` share `m` members — so they are far more
> correlated than random subsets of the same sizes (lag-1 autocorrelation **0.917 against 0.112**;
> 25.3 local minima against 42.0). **A family of less-correlated variants has a larger per-target
> minimum.** *The prefix ordering does not merely fail to beat an arbitrary index — it loses to
> one*, and `2.9027` is therefore a **best-of-K order statistic**, not a property of the score axis.
> Its `K_eff` is well under 128: only **6.3** values of `m` lie within 0.01 Å of the minimum and
> **30.1** within 0.05, so the argmin is not sharply identified.

**Geometry is valid on every arm** (virtual bond 3.80395 Å, sd ~1e-15), so none of this is an
invalid-structure artefact — the check that `AVG_SEP` failed in §14.4.

---

## 15. Statistical analysis

**MDE = 2.8016 × SE, per comparison.** Below **0.7×** is *not a result*; **0.7–1.0×** is *NOT
MEASURED*. Fold-clustered CIs on the pinned folds; `s24.stats_lib.compare` **refuses a verdict** if
`folds` is not passed.

**Multiplicity was tracked sprint-wide and written to as comparisons were emitted**, not
reconstructed at the end — the build item S30 carried forward. `s31/MULTIPLICITY.md` holds the
register; the running total is in the hundreds with a Bonferroni multiplier table computed up front
(k = 32 → 1.43×, k = 912 → 1.74×). **One lane's first total was wrong by 145 and was corrected in
place with the error stated.**

**The verifier** (`s31/s31_verify.py`) recomputes **129 numbers from artefacts, 0 mismatched, 0
missing, 0 flagged.** It **parses the sprint's own documents for every path they name** (28 found, 28
exist) rather than using a hand-kept list — *because a hand-kept list is exactly what fails*. It
checks both reporting bases, asserts `MDE == 2.8016 × SE` from the SEs, and carries forward S30's
sign-convention trap. **It caught two defects in its own author's work.**

### Three places the rules bound against the person applying them

1. **A fold CI that excludes zero does not rescue a sub-MDE effect.** The branch-carry arm came in at
   **0.44× MDE** with a fold CI excluding zero and 5/5 folds agreeing — *and was reported as NOT A
   RESULT*. The widening arm is the same shape at **0.42× MDE** with CI [+0.0028, +0.0557]. This is
   the documented sibling of the underpowered bug (`s24/stats_lib.py:145`), and **the MDE gate
   binds.**
2. **A null in the mean is not "no effect".** `P1`'s per-target spread is **10.7× the implementation
   noise null**, with 74 of 126 targets moving more than that null's p90, while the mean is −0.0112
   at 0.19×. The honest sentence is *"it reshuffles the answer everywhere and buys nothing."*
3. **The per-target floor is not the mean floor.** 0.0107 Å bounds *mean* built-chain claims; the
   per-target floor is **~0.03 Å with a 0.23 Å tail** (§1.2).

### Registered predictions, scored

Ten pre-registrations were committed before their first number. **Predictions that failed are
recorded as failures, including the coordinator's.**

| prediction | by | outcome |
|---|---|---|
| The orthogonal complement is noise (2:1 against E2/E3) | coordinator | **Right for the wrong reason** — it is 94% predictable and *actively harmful*, not noise |
| Long-range R² ≈ 0 (3:1) | coordinator | **Falsified** on the raw statistic (+0.1959) |
| A substantial fraction of `bestm128` survives (2:1) | coordinator | **Falsified** — the surviving fraction is *negative* |
| `LEG_torsion` is predominantly odd (4:1) | lane B | **Failed** — it is 70% even |
| `F-G2` does not fire | lane G (S30) | Held, and its *mechanism* was refuted by its own pre-check |
| The derived readout beats production ([−0.15, +0.10]) | lane A | **Missed at +0.1334**, on the side that says the direction fails |
| Disagreement count 20–45 | lane P | **Wrong — 66** |
| Arm E lies between D and F | lane P | **Falsified in premise and outcome** |
| `p*` does not improve on the VQE (2:1) | lane P | Held |
| A deployable prefix-`m` rule transfers | lane F | **Falsified**, three ways |

> **A pre-registration that only ever confirms is decoration.** Six of these went against the lane
> that wrote them, and three of those six are the coordinator's.

---

## 16. What was closed

Several of these are theorems. Where a direction closed by **derivation before compute was spent**,
that is noted — it is the charter's §24 ladder working as intended, and it is most of this table.

| direction | closed by | how |
|---|---|---|
| **The deployed CVaR objective as a quantum problem** | **theorem** | It is a **convex program with a closed-form global minimiser** pinned by one scalar. `run_cvar_vqe` is strictly worse in **126/126** targets at the deployed settings; `p*` costs 0.0012 s against 0.0985 s. Charter §11 closed for the deployed objective |
| **Any target-specific role for the quantum state** | **derivation + measurement** | `E = _zrank(...)` is a **target-independent constant** (§5). `H(p*)` has sd **3.1e-04** across 126 targets. The stage answers a global hyperparameter question |
| **A non-diagonal Hamiltonian** | **three independent theorems** | CVaR needs a per-shot eigenvalue; the forced operator is mean-field and quartic in ψ; and a candidate-index register's Hilbert dimension *is* the candidate count, so no operator on it can be classically hard |
| **ADAPT-VQE / qubit-ADAPT** | **theorem** | Its selection rule `\|⟨ψ\|[H,A]\|ψ⟩\|` presumes the cost is `⟨H⟩`, a linear functional; **CVaR is not the expectation of any observable**, so the criterion is *undefined*, not merely unhelpful |
| **The free-energy stage (§7A)** | **derivation, no compute** | `S` has no target argument, and the shipped potential is **reflection-invariant** (1340 torsion phases at distance 0.000e+00 from {0,π}), so `F`, `E`, `S` and every `dF/dT` are distance-map functions by G1 |
| **The elastic-network / normal-mode / landscape-curvature family** | **derivation, no compute** | **Corollary B1′** — an ANM/GNM Hessian is built from pairwise distances, so its spectrum and log-determinant (*which is the harmonic entropy*) are distance-map functions |
| **The backbone-torsion channel (§7B)** | measurement | The chiral escape from G1 is real (odd half 2.78× MDE on coarse triage) **and only works on a problem the pipeline does not have** (0.28× in-pool). A **constant α-helix beats every torsion channel on both bands** |
| **Charter §14, the new-observable question** | **provenance** | 92.9% of the benchmark is NMR-determined; for **117/126** the deposited coordinates *are* a fit to the deposited restraints, so any NMR observable is **ORACLE through a different door**. The one genuine escape (VCD/ROA, chiral, works at 9–16 residues) has **zero measured spectra** for these targets |
| **Candidate-index redesign (§12)** | measurement | Gray coding is a **proven no-op**; the best map raises an ORACLE ceiling and **lowers** the deployable value; the best map for a deployable partial readout is the **random permutation**; no index bit carries more than **0.066 bits** about candidate quality |
| **The per-target prefix length `m`** | **matched control** | A per-target minimum over 128 **random subsets** reaches **149%** of the prefix family's gain — the prefix axis is *worse than an arbitrary 7-bit index*. And the transfer was already **FALSIFIED in S29** (−0.0018 global, +0.0079 leave-fold-out) |
| **The constrained-affine readout** | **by exhaustion of selection rules** | No knee by ceiling (the ridge path is smooth and monotone) and no answer by fit (out-of-fold R² 0.021, applied +0.0309 worse). **It cannot be selected by ceiling and cannot be selected by fit** |
| **Consensus as a quality estimate** | **identity** | The consensus criterion **is** the free half of the objective read at uniform weights — `corr(medoid criterion, B·1/D) = 0.9667` over all 126. Setting `â ∝ +consensus` **cancels** the term it was meant to complement |
| **My own opening hypothesis (E1)** | **exact identity** | `mu_hat = mu − y`, so the proposed common-mode estimator's error **is** the prior error it was meant to help predict. Circular |
| **The orthogonal-correction family (E2/E3)** | **ORACLE ceiling** | With perfect magnitude *and* a perfect direction, correcting only the component orthogonal to the common mode is worth **+0.0747 Å — harmful** (§20) |

---

## 17. What improved

**Nothing at the endpoint.** No arm was deployed; every arm tested is NOT MEASURED or worse.

**What did improve is the instrument and the record**, and several of these gate future measurement:

| improvement | why it matters |
|---|---|
| **The projection's real defect diagnosed** — there is no RNG; the λ=0 argmin is decided at a 1e-7 spread among branches 1e-1 apart, amplification ~1e13 | A **seed** story had propagated through four sprints' documents and a verifier. Now corrected in `core/project.py` itself, where the next person will look |
| **The per-target built-chain floor measured** — 0.0134 mean / 0.0329 p90 / **0.2285 max**, from the *same operator written twice* | The 0.0107 Å figure had been quoted as bounding per-target statements. **It does not.** This calibrates every arm-to-arm comparison the project makes |
| **A withdrawn positive removed from shipped code** (`core/pipeline.py`, the +0.113 Å CVaR-tail claim) | Recorded as a defect in S26 and **still wrong four sprints later**. *Recording a defect in a table is not fixing it* |
| **The launcher now imports its gate from the governor** rather than duplicating it | The governor changed four times in S29 and the launcher never followed. **The pair can no longer drift** — the first fix that makes the next instance impossible |
| **A cross-basis audit with a self-test** | The verifier passed 97/97 on *"is a basis named?"* while the headline carried a **misnamed** basis. **An audit that checks for an absent label cannot catch a wrong one** |
| **The verifier grew to 136 checks** and parses the sprint's documents for every path they name | A hand-kept list is exactly what fails |
| **A sprint-wide multiplicity register**, written to as comparisons were emitted | The build item S30 carried forward |

---

## 18. What did not improve

**The endpoint, and every arm proposed to move it.**

| arm | basis | result |
|---|---|---|
| `p*` substituted for the VQE's `p`, shipped readout | built chain | **−0.0112, 0.19× — NULL** |
| `p*`, selection readout | built chain | −0.0280, 0.29× — NULL |
| The shipped quantum stage, on | built chain | +0.0175, 0.35× — NULL (and a **cost**) |
| The derived convex readout at γ = 1 | CA cloud | **+0.1334, 1.17× — WORSE** |
| Quality-blind dispersion maximisation | CA cloud | +0.1436, 1.22× — WORSE *(and its one positive against the argmin was killed by its own shuffled-`B` control)* |
| `DIS + CONS` leave-fold-out through the convex readout | CA cloud | +0.1347, 1.71× — WORSE |
| `AVG_SEP` — per-separation restore | built chain | **+0.4609, 2.34× — MUCH WORSE** |
| `AVG_RG` — uniform Rg restore | built chain | +0.0475, 1.11× — WORSE |
| The consensus medoid | built chain | +0.0688, 0.81× — NOT MEASURED |
| Branch-carry at λ = 0.3 | built chain | −0.0055, 0.44× — NOT A RESULT *(a conditioning fix, 1082×, with a null accuracy effect)* |
| The per-target prefix length `m` (leave-fold-out, the deployable arm) | built chain | **+0.0075, 0.22× MDE, 2/5 folds, 63W/63L — a coin flip** (§14.6). *The **1.2%** figure quoted before is **CA cloud**; the endpoint reads −2.4%, and **both are nulls** — the conclusion is unchanged, it is now measured on the basis it is stated on* |
| Orthogonal-only prior correction (E2/E3) | built chain | +0.0414 / +0.0323 — **FALSIFIED as registered** |
| Native-free index maps | CA cloud | Best −0.0160; **the best map is the random permutation** |

**Not one deployable arm in the sprint reached its own MDE in the helpful direction.**

---

## 19. ORACLE ceilings

Every number here is **ORACLE / NOT DEPLOYABLE**. They bound what an operator could reach with
perfect information; none is achievable, and several are vacuous for reasons worth stating.

| ceiling | basis | value | note |
|---|---|---|---|
| affine readout over 128 | CA cloud | **0.0000** | **vacuous by rank** — 127 weights against ~33 residual dimensions |
| convex readout over 128 | CA cloud | **1.7977** | support 6.54/128 *without* a sparsity penalty; 9.59× MDE, 119W/0L |
| convex readout over 128 | built chain | 1.8538 | S29's rows, quoted unchanged |
| argmin over 500 (9 bits) | built chain | 1.7078 | |
| argmin over 128 (7 bits) | built chain | 2.1435 | **1.067 Å of headroom below production** |
| per-target prefix length `m` | built chain | 2.9027 | **an order statistic — see §12**; confirmed on the chain in §14.6, where a matched random family reaches **141%** of it |
| best pool member, genuinely worst 18 | built chain | 2.5298 | the tail is selection-limited |
| branch choice among the 8 already computed | built chain | **−0.0938 vs production** | 2.92× MDE, 114W/0L; free to compute |

**Two of these deserve their caveat repeated in any quotation.**

**The affine ceiling is not a bound.** `rank(aff{W_x}) = 32.9 ≥ 3n−3`, so the hull spans the space
and the "ceiling" is an interpolation artefact. The source file's own instruction is *"do not quote
any affine-hull ceiling as a bound."*

**`bestm128 = 2.9027` is an order statistic, and the framing that was licensed is narrower than the
one that was used.** See §12. The *ceiling* interpretation stands — **an optimistically biased upper
bound is still a valid upper bound**, so "2.5 Å is unreachable through this architecture" is safe
and **safer than stated**. The *lead* interpretation — "0.308 Å available to a better selector" —
was priced at **−0.0018 global / +0.0079 leave-fold-out** and marked **FALSIFIED in S29**, and the
caveat was dropped in re-quotation.

---

## 20. The remaining information bottleneck

### 20.1 The requirement was stated backwards, and the correction is the sprint's most consequential result

S30 concluded — and this report's own contract rule 29 encoded — that the next channel must supply
an observable whose error is **incoherent** with the pool's common mode. **Measured with ORACLE
magnitude and a perfect direction, that is worth `+0.0747 Å`: harmful.** The decomposition, splitting
the ORACLE ideal correction `y` against the ORACLE common mode `mu` per target (built chain, all
**ORACLE / NOT DEPLOYABLE**):

```
correcting ONLY the along-mu half    -0.8102 A   3.29x MDE, 5/5 folds, 121W/5L   <- the whole prize by SIZE
correcting ONLY the orthogonal half  +0.0747 A   0.87x MDE, 57W/69L              <- worth less than nothing
      ^ the SIZE is real; that it is a separate MECHANISM is NOT -- the along arm is NOT MEASURED against the shrink curve
along vs perp, paired                -0.8848 A   3.50x MDE, 5/5 folds, 118W/8L
```

And the reason is an **identity, not a fit**. Refitting the same ridge on each half:

```
y_PERP_mu    out-of-fold R2  +0.9403    (random-feature control -0.0026)
y_ALONG_mu   out-of-fold R2  -0.0051    (random-feature control -0.0039)
```

`mu_hat = mu − y` forces `y_perp = −P_perp(mu_hat; mu)` **exactly to 9.3e-15**, and with
`cos(mu_hat, mu) = 0.05` that makes `y_perp ≈ −mu_hat` at corr 0.971. **The predictable half is not
predicted — it is observed.** The unpredictable half is exactly the common mode that S30-L7 proved
non-identifiable from pool data at any K.

> ### What can be predicted is the component ORTHOGONAL to the common mode, and it is harmful.
> ### What would help is the common mode itself, and it is unpredictable.

**S30's slogan survives verbatim; its mechanism was inverted.** The measurement behind it (+0.0554
coherent against −0.2466 i.i.d. at matched R²) is **untouched and still stands** — it was explained
backwards, not measured wrongly.

**Controls, all registered before the numbers:** magnitude is not the explanation (the perfect
correction shrunk to the perp arm's exact norm gives −0.5422, so the perp arm is **+0.6169 worse
than its own magnitude-matched control**, 3.25× MDE, **9W/117L**); the projection is not the
mechanism (against a norm-matched shrinkage, 0.26× MDE); and projecting against the *true* `mu` buys
nothing over not projecting (0.18× / 0.08×).

**The last control, and it is the one that could have killed the result.** The asymmetry could be
*magnitude* (the along component is simply bigger) or *geometry* (splitting any vector against any
direction produces this pattern). Both are now tested at n = 126, **ORACLE / NOT DEPLOYABLE**, built
chain. An **energy-matched** direction `u` is built per target to capture *exactly the fraction of
`y`'s energy that `mu` captures* — `cos²(t) = ‖y_along(mu)‖² / ‖y‖²` — but is otherwise arbitrary:

```
                                        effect   xMDE  folds    W/L     reading
energy-matched ALONG  vs true along-mu  +0.1317  2.07   5/5    25/101   true mu is BETTER
energy-matched PERP   vs true perp-mu   -0.3797  2.53   5/5   114/12    true mu is WORSE
```

**`mu` is special on both sides, both at 5/5 folds and both past MDE.** Its along-component is worth
0.13 Å *more* than an arbitrary direction carrying identical energy; its perp-component is worth
0.38 Å *less*. The control was built so that a null result would have refuted the section, and it
did not fire.

**Which half of the finding that actually controls, stated plainly, because the two halves are not
equally earned.** `y_along(mu)` carries **84% of `y`'s energy** (rms 3.12 of 3.70), and against the
shrink curve — *its own* magnitude control — it is **NOT MEASURED in both directions**:

```
ALONG-mu vs SHRINK_Y_075   -0.0614   0.94x MDE   81W/45L   NOT MEASURED
ALONG-mu vs SHRINK_Y_090   -0.0405   0.64x MDE   65W/61L   NOT MEASURED   <- a coin flip
```

***"Correcting along `mu` is the whole prize" is, at matched magnitude, the same statement as "most
of a perfect correction is the whole prize."*** The along arm is not a separate mechanism; it is 84%
of the right answer, delivered by a direction that happens to select it. **The controlled half of
this section is the NEGATIVE one** — that the orthogonal complement is *specifically* worthless,
worse than its own norm-matched shrinkage (**+0.6169, 3.25×**) and worse than an energy-matched
arbitrary direction (**+0.3797, 2.53×**) — **and the negative half is exactly what charter §14 and
contract rule 29 proposed to build on.**

**`mu`'s direction still does the selecting, though, which is why the section is not vacuous.** An
*energy-matched arbitrary* direction's along-half **is** distinguishable from the same shrink
(**+0.0912, 1.43× MDE**) while `mu`'s is not. *Any direction capturing 84% of `y` gets 84% of `y`'s
magnitude; only `mu` gets 84% of `y`'s* **value**.

> **A trap worth recording, because I walked up to it.** `CTRL_SHRINK_ORACLE_Y` is the perfect
> correction shrunk **to the perp arm's norm (1.99)**, not the along arm's (3.12). Differencing the
> along arm against it reads **−0.2680 at 2.06× — an apparently controlled positive** — and it is
> meaningless, because the control is matched to the *other* arm. That is this sprint's signature
> defect exactly (**a number carries its definition**), one boundary from being quoted as the
> headline. The along arm's magnitude control is the shrink curve, and it says NOT MEASURED.

**But the asymmetry is not entirely `mu`'s, and that qualification is the honest one.** At matched
energy an arbitrary direction *still* shows along-beating-perp, by **−0.3734 (1.98× MDE, 5/5,
84W/42L)**. Against the true split's **−0.8848**, roughly **42% of the asymmetry is generic** and `mu`
supplies the remaining 58%. ***`mu` more than doubles an asymmetry that is already there — it does
not create one.*** The mu-specific part is where the prize is: **−1.4697 on FAIL18 against −0.1907
elsewhere.**

**The random-direction null behaves exactly as pure magnitude predicts, which is what makes it a
scale.** An isotropic direction captures almost none of `y`, so its along-arm is **−0.0229 at 0.66×
MDE — NOT MEASURED** and its perp-arm is nearly the entire correction (−0.7884). *Given a meaningless
direction the labels carry no meaning, and the split collapses to how much of `y` survived it.*

**The shrink curve prices accuracy against direction** (ORACLE, built chain, paired against the
same-job production):

```
c        0.25     0.50     0.75     0.90     1.00
delta  -0.2978  -0.5857  -0.7488  -0.7697  -0.7756
frac      38%      76%      97%      99%     100%
```

**Three quarters of a perfect correction buys 97% of its benefit; a quarter still buys 38%.** The
curve is steeply concave and effectively saturated by `c = 0.75`. ***A future common-mode channel
does not have to be accurate — it has to point the right way.*** That is the most encouraging number
in this report, and it is why §21 prices the *direction* and not the magnitude.

> **Instrument fact, and it sharpens when cross-job chain comparison is legal.** These control arms
> were computed in a **different process** from the `mu`-split they are differenced against — which
> lane D's ~1e13 amplification result would normally forbid. It is legal here, and checked rather
> than assumed: the two jobs' production rows are **bit-identical on all 126 targets in both bases**
> (max |diff| exactly `0.000e+00`, 126/126 exact ties). **The projection is deterministic; what
> lane D measured is sensitivity to *differing* inputs, not nondeterminism.** So the rule is sharper
> than *"always project in one job"*: **identical clouds give identical chains to the last bit, and
> the check costs one line.**

**Sizing, ORACLE / NOT DEPLOYABLE:** the along-`mu` correction alone puts the built chain at
**2.4025 Å** — which would clear the charter's *ambitious* 2.50 target — and it is **4× concentrated
on the tail** (FAIL18 6.0195 → 3.7392, −2.2803, against −0.5651 on the other 108). It does **not**
beat a perfect prior: ALONG against FULL is −0.0345 at **0.51× MDE, 3/5 folds — a tie, NOT
MEASURED**, and it is reported as a tie.

### 20.2 The second bottleneck, which is a different object

Lane A's closing statement, from the exact readout identity:

> The readout's objective is exact and half of it is free; the conversion from any quality estimate
> to an endpoint is a tuning-free convex program; that program's ORACLE ceiling is **1.7977 Å**
> against production's **3.0483** (CA cloud). **What is missing is a per-candidate quality estimate
> with positive IN-BAND skill.**

Every native-free candidate now measured has in-band skill that is **zero or the wrong sign**:

```
DIS z-rank (the shipped cost)    rho global +0.1176    IN BAND  -0.0262
CONS (the medoid criterion)      rho global +0.4278    IN BAND  -0.2837
DIS + CONS, leave-fold-out       rho global +0.4281    IN BAND  -0.2827
```

**Consensus reaches double the crossing price on the global axis and the wrong sign in band** — its
`ρ_global = +0.428` is **entirely outlier detection**. *The binding axis is in-band `ρ`, and the
global crossing price of 0.211 is necessary-not-sufficient: it holds along an interpolation path and
must not be quoted as a target for a real feature.*

### 20.3 The two bottlenecks are different objects, and whether they are one requirement is open

`â` is *which candidate is better*. `mu` is *how all of them are wrong together*. **Both are things
the pool cannot tell you about itself**, and both would have to come from outside it. Whether they
are two faces of one requirement or two separate ones is **the sharpest question this sprint
produces, and it is stated as open rather than resolved by assertion.**

### 20.4 What is NOT the bottleneck, so the next sprint does not pay for it again

The circuit's expressivity (`p*` is free and closed-form, and substituting it is worth **−0.0112 Å at
0.19× MDE**); the optimiser (80 → 2000 iterations move the gap by nothing); the CVaR α (inactive by
construction on three of five folds); the ansatz depth; the register width (the 75 → 128 widening is
a **+0.0307 Å COST** at **0.42× MDE, NOT MEASURED** — paid only to fill the `2**7` register
(`core/pipeline.py:757`) and pointless when `quantum = False`, which is production; *its fold CI
[+0.0028, +0.0557] excludes zero while the effect is 0.42×, the exact shape `_verdict` was hardened
against, so the MDE gate binds and the CI does not rescue it*); the index encoding; the readout
class (exact, and its
native-free channel is **3.99%** of ORACLE); and the candidate set (the set effect is **+0.0049 Å at
0.14× MDE — NOT A RESULT**, while the weighting rule costs **+0.1102 Å at 2.68× MDE**).

---

## 21. Next-sprint recommendation

**One question, stated as a requirement rather than a direction.**

> ### Find a source of information about the POOL'S COMMON MODE — the shared component of the retrieved candidates' error — that does not come from the pool.

This is the inverse of what S30 and this sprint's own opening contract asked for, and the inversion
is measured (§20.1), not argued. Three things make it a well-posed target rather than a wish:

1. **It is priced.** Perfect knowledge of the along-`mu` component alone is worth **−0.8102 Å on the
   built chain (ORACLE / NOT DEPLOYABLE)**, taking the endpoint to **2.4025 Å** — past the charter's
   *ambitious* target — with **4× concentration on the tail**.
2. **The conversion is free and tuning-free.** The readout identity makes any quality estimate
   convert to an endpoint through a **convex program with no hyperparameters** (§3(ii)).
3. **The failure mode is named.** The common mode is **non-identifiable from pool data at any K**
   (S30-L7), so a channel that reads the pool — however cleverly — cannot supply it. *That is the
   one thing we now know for certain about where it must come from.*

**How big the ask is, stated so the framing cannot flatter it.** `y_along(mu)` is **84% of the
perfect correction's energy**, and applying it is **NOT MEASURED** against simply applying 75–90% of
that perfect correction (§20.1). **So "supply the common mode" is not a clever projection that
extracts a large gain from a small piece of information — it is a request for most of a perfect
prior, expressed in the coordinates that make the *remaining* structure legible.** The reason to
state it in those coordinates anyway is the part that *is* controlled: **the 94%-predictable
orthogonal complement is worth +0.0747 Å, which is worse than doing nothing.** *The project's
corrector is not weak — it is aimed at the wrong component*, and no amount of improvement to it
changes that, because it is an **identity** (`mu_hat = mu − y`) and not a fit. **That is the
durable result: it rules out a whole family of next steps rather than endorsing one.**

**And the parallel requirement, which may or may not be the same one:** a per-candidate quality
estimate with **positive in-band skill**. Every native-free candidate measured has zero or negative
in-band skill (§20.2). **State which of the two any proposal addresses.**

### What NOT to spend the next sprint on

- **Any further readout, index, ansatz, optimiser, α-schedule or register-width work.** All are
  closed above, most by theorem, and §20.4 lists them so the case does not have to be re-made.
- **A deployable 128 → 512 widening**, on the prediction in §11 — three facts held by three
  different lanes compose to say it should come out **worse** than production, not flat.
- **Chiral single-structure functionals at 9–16 residues.** Empty here, twice. **The theorem is
  length-free and only the emptiness is length-dependent**, so this is worth retesting at 40+
  residues on an instrument that does not exist yet.

### Two engineering items that gate future measurement

- **The per-target built-chain floor is ~0.03 Å, not 0.0107** (§1.2). Any S32 claim below that on
  individual targets is indistinguishable from re-running the same code twice. **The branch-carry
  fix (B4) is available, improves the conditioning 1082×, and is a NULL for accuracy at 0.44× MDE**
  — adopt it deliberately at a sprint boundary as a *conditioning change*, never as an improvement,
  and note it moves canonical 3.2105 → 3.2050.
- **The verifier should assert each number against the value on the basis it names**, not merely
  check that a basis is named. An audit for an *unstated* basis does not catch a *misstated* one,
  and this sprint had exactly one of the latter — in its own headline.

---

## Appendix A — every claim withdrawn this sprint

Recorded in full because the charter requires it and because it is the reason the surviving results
are worth anything. **The distribution is the finding**: every single-lane result held; **every
cross-lane synthesis by the coordinator failed.**

### A.1 The coordinator's

| claim | how it died |
|---|---|
| **R1's quantifier** — *"the entire quantum stage carries at most k bits"* | **My own falsifier fired in shipped code.** I wrote *"R1 fails if any code path lets the stage emit a structure that is not a pool member"*; lane C found `average_weighted` (`core/pipeline.py:880-895`), **shipped**, called at `:1110`/`:1116`, emitting a continuous convex combination — measured **1.1144 Å from the nearest pool member**. The theorem holds for the *selection* readout only |
| **R1 point 4** — all `2^k` vertices reachable | Duplicate structures make `P[i,j] = 0` off-diagonal, so `argmin` returns the first index. **Capacity is 6.886 bits, not 7** — and the reachable count equals the byte-distinct count on every target |
| **"Two readouts"** (S31-L2) | There are **three**, and I merged the two that matter most to distinguish: the **convex** one ships, the **affine** one does not |
| **The +0.2260 Å deficit**, broadcast to four lanes | It is the **affine harness readout's** cost (`circ_l1_i80` lives in `s28_A_amp.py`), not the shipped stage's. Lane A caught it as a **12.7× inconsistency** with its own cloud measurement and refused to resolve it by assuming the projection amplified. The shipped stage costs **+0.0178 Å cloud / +0.0175 built chain, both NULL** |
| **"Averaging contracts the backbone 25.8%"**, which I built a lane's entire mechanism on | **A number withdrawn two sprints earlier** (`s15/coord_FINDINGS.md:914-921`; correct: **3.5%**, and it is a *separation-dependent distortion* crossing 1.00 near \|i−j\| = 8, not a contraction). Still stale in S30's ledger and in project memory when I quoted it |
| **My opening hypothesis (E1)** — estimate the common mode as `mu_hat = pool75_mean − expected` | **Circular, killed by one line of algebra.** `mu_hat = mu − y` exactly (1.8e-15), so the estimator's error **is** the prior error it was meant to help predict. To use it you would already need the answer |
| **"You need a measurement of the molecule, not a computation"** — nearly the sprint's headline | **True premise, vacuous conclusion.** By Anfinsen `I(N; sequence) = H(N)`, so the bound reads *"the pool contains at most everything"*. **If it implied a ceiling, AlphaFold would be impossible** — and our own ESM result (+0.288 Å over one-hot on the same sequence) could not exist |
| **"The set-matched ladder inverts S30-L11"** | **A withdrawal that was itself wrong** (contract rule 25). S30 computed the correct reference curve **in the same file**: argmin reaches **1.7108 at 9 bits** against the sparse arm's ~9.1. **S30's support was mismatched; its conclusion survives** |
| **"0.076 Å of pure error cancellation at zero weight bits"** | Two errors at once. The support is **ORACLE-greedy at 12.99 bits**, which S30 *did* charge; and at **matched search size** it is a **0.077 Å penalty** — a minimum over C(128,2) = 8128 read against one over 128 |
| **The headline's basis** — *"2.1435 Å (CA cloud) … ~0.90 Å of headroom"* | **2.1435 is the built chain** (cloud is 2.1458), differenced against the cloud production — so the headline **understated its own headroom by 0.16 Å** in the sentence carrying the framing claim |
| **"If `coh(AVG_SEP) < 0.6931` it is the first native-free operator to pass"** | **The bar is measured on a different object.** S30's `coh` grades a *corrector's* residual (the distogram's prediction error, an **input**); lane B's and lane F's grades the *emitted structure's* error (the **output**). Same pipeline, same `mu`, **0.6931 against 0.9780** — two different errors. It reached **shipped code in two lanes** before anyone read the definition at source |
| **The lever comparison** — *"the convex readout is 0.290 Å better, so the readout class is the larger lever"* | **Unequal information cost**: 7 bits against **128 free reals**. A larger *ceiling at unpriced cost*, not a larger lever — **the error S30 §9.4 had already named as this project's characteristic one** |
| **"`bestm128` may deflate"** | **Inverted.** The transfer arm existed *in the entry that produced the number*, and S29 marked it **FALSIFIED**. Order-statistic inflation makes an oracle number **optimistically biased, and an optimistically biased upper bound is still a valid upper bound** — so the ceiling reading was licensed. **My use of it as a "0.308 Å lead" was not**, and I had dropped S29's caveat in re-quotation |
| **The widening null as 0.35×** | That is `D − A`; the widening is `F − A2` at **0.42×** |
| **My 2:1 reasoning on E2/E3** | Right in direction, **wrong in mechanism, twice**: I predicted the orthogonal complement would be *noise*. It is **94% predictable and actively harmful** |

### A.2 The lanes'

| lane | claim | how it died |
|---|---|---|
| **A** | The derived readout as a deployable improvement — its **registered primary** | **+0.1334 Å, 1.17× MDE, WORSE**, outside a registered band of [−0.15, +0.10], *missed on the side that says the direction fails* |
| **A** | Its one positive — quality-blind dispersion beating the shipped argmin by −0.2621 Å at 2.74× | **Demolished by its own shuffled-`B` control at 0.54× MDE.** The mechanism is *"spread the weights"*, not *"spread along the real geometry"* |
| **A** | `ρ = 0.211` as a crossing **target** | Consensus reaches **double** it and is the **wrong sign in band**. The global price is necessary-not-sufficient and holds only along an interpolation path |
| **A** | Its popcount explanation of the optimisation residual | Registered and **refuted** by its own relabel intervention; **no post-hoc replacement offered** |
| **A** | +0.0178 / +0.0125 quoted without MDEs | **NOT A RESULT** (0.37× / 0.27×) — self-corrected and propagated |
| **B** | 4:1 that `LEG_torsion` is predominantly odd | **It is 70% even** |
| **B** | That the 32-cost sweep's contrast was a `RAND_SIGNED` artefact | **Refuted by its own control** — `GAUSS_MATCHED` gives the same contrast. The construction is exonerated |
| **B** | That `LEG_steric`'s contrast is chirality-dominated | Its **variance** is (odd share 13.02); its **contrast** is not (odd/tot 0.49). Self-corrected |
| **C** | A wrong-tail sign error in C3 | Defined a contrast that could only rise, then read the null's lower tail. **The first run printed the opposite verdict; caught before any number left the lane** |
| **C** | The ORACLE convex ceiling | **Violated its own feasibility bound on 4 of 126** (3 of them FAIL18). Self-audited, vertex added as a third candidate, bound asserted: **119W/4L → 119W/0L — the four losses *were* the four failures** |
| **C** | *"must never again be quoted as the architectural ceiling"* | Withdrawn by lane C as **too strong** |
| **D** | Two defects in its own work | Caught by its own verifier |
| **E** | *"It must carry orthogonal INFORMATION"* | **Half withdrawn by its author.** Orthogonal information is exactly what we have, and it is worth **+0.075 Å** |
| **F** | Its `hash()`-seeded random family | Python salts `hash()` per process, so the subsets were not regenerable. **Repaired to `crc32` and the control re-run from zero** — 146% → **149%**, conclusion strengthened, **draw sd 4.2× larger** |
| **F** | A missed grep — `medoid75` was already on disk | **Declared as its own defect in its prereg before anyone could find it**, and its derivation had predicted the sign and margin first |
| **L** | *"Consistent with the unpinned projection seed"* | Withdrawn — **there is no RNG on that path.** *"I inferred from the charter's defect list rather than from the code, the same error this lane exists to catch in others"* |
| **L** | Its 12 verification cells | Run at **T = 0.1 / 0.05**; the deployed `T` is **0.3 on every fold**. Its entropy-direction claim **reverses** there |
| **P** | Disagreement count 20–45 | **66** |
| **P** | *"Arm E lies between D and F"* | **Falsified in premise and in outcome**, and the premise failure was flagged by its author before the outcome was known |
| **P** | Lane L's *"there is no third outcome"* | **There was, and it happened** — a null at 0.19× MDE |

### A.3 The rules this sprint earned

1. **A claim combining two lanes' numbers needs a named owner who holds both**, and must carry both
   caveats in the sentence carrying the number. *Every single-lane result held; every cross-lane
   synthesis failed — because a single-lane claim is audited by the lane that owns the data and a
   cross-lane claim is audited by nobody.*
2. **Every lever comparison carries its information cost in the same sentence as its Ångströms.**
3. **A bar imported from another sprint must be checked against the object it graded there.** The
   `0.6931` bar crossed from corrector-space into readout-space with no owner holding both, and
   reached shipped code in two lanes.
4. **A fold CI excluding zero does not rescue a sub-MDE effect** — twice this sprint, at 0.44× and
   0.42×.
5. **A null in the mean is not "no effect"**: report the per-target distribution against a matched
   implementation-noise null.
6. **An optimistically biased upper bound is still a valid upper bound.** Order-statistic inflation
   deflates a *lead*, not a *ceiling*.

---

## Appendix B — the five "worst 18" strata, as a key

Five different strata are called "the worst 18" somewhere in this project's record. They differ by up
to **0.81 Å on the same-named quantity**. Every such number in this report names its stratum in the
same sentence.

| stratum | defined by | `best1_500` |
|---|---|---|
| `FAIL18` | the filter's own recall (`s12/instrument.py:271-278`) | 2.2842 |
| `defn18` | top-18 by widening gain | **= FAIL18, 18 of 18 — an identity** |
| `worst18_poolmean` | pool mean | 2.2286 |
| `worst18_bestpool` | best pool member | 3.0930 |
| "the genuinely worst 18" | production built-chain RMSD | ORACLE best pool member 2.5298 |
| worst 18 by production RMSD, n = 111 matched | ensemble arm | mean RMSD 6.0291 |
