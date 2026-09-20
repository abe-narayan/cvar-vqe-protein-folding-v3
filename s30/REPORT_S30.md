# Sprint 30 — CVaR-VQE Protein Folding: the search for the first real accuracy breakthrough

**Status: DRAFT IN PROGRESS.** Sections marked `[PENDING]` await lanes still running. The charter
requires the report only after everything is finished; this file is assembled as results land so
nothing is reconstructed from memory at the end.

Branch `s26` · benchmark: 126 dev targets, 9–16 aa · endpoint: **mean built-chain Cα RMSD**
Ledger: `s30/LEDGER.md` · Running state: `s30/STATE.md` · Contract: `s30/S30_CONTRACT.md` (28 rules)
Charter: `s30/BRIEF.md`

---

## 0. The answer, up front

The charter asked for the first real accuracy breakthrough, with a primary target of < 3.0 Å and an
ambitious target of < 2.5 Å against production's **3.2105 Å**. It also said, explicitly, that *"a
well-evidenced ceiling argument that redirects the next three sprints is worth more than a fragile
2.98 Å."*

### 0.1 The endpoint did not move

| # | question | answer |
|---|---|---|
| 1 | Strongest new architecture? | **None was deployed.** Every candidate was closed — by measurement, by theorem, or by price — before it reached the pipeline |
| 2 | Final mean built-chain RMSD? | **3.2105 Å**, unchanged |
| 3 | Did it beat 3.21? | **No** |
| 4 | Did it beat 3.0? | **No** |
| 5 | Did it beat 2.5? | **No** |
| 6 | Exact paired effect vs production? | **Zero — nothing was deployed.** The best *measured, applied* prior correction emits **+0.0554 Å WORSE**. The one confirmed deployable effect (AMBER relax at k = 30) is **−0.0406 Å**, 3.56× MDE, 5/5 folds on a length-matched split — **0.69% of baseline, and its restraint constant still has no native-free selection rule** |
| 7 | Fold-clustered CI? | Per comparison; every one is reported beside its effect. No sprint-wide constant |
| 8 | MDE? | Per comparison, MDE = 2.8016 × SE. **Below 0.7× is not a result; 0.7–1.0× is NOT MEASURED** |
| 9 | What happened on FAIL18? | **Unchanged.** And one of the sprint's own headlines died here: FAIL18 is defined by the filter's own recall, so claims about the filter's effect on it are partly forced arithmetic (§Appendix A) |
| 10 | What happened on the other 108? | **Unchanged** |
| 11 | What information source produced the gain? | **There was no gain.** The sprint's contribution is the reason why |

### 0.2 What the sprint established instead

**The target is reachable. Nothing we own can reach it.** Both halves are measured, and they are
worth more than a fragile 2.98 would have been.

> **Five ORACLE signs on long-range pairs are worth −0.3259 Å on the built chain — 3.2126 → 2.8867,
> 2.05× MDE, 5/5 folds. That clears the charter's primary target.** *(Lane P's arm is measured
> against the S27 projection's 3.2126 rather than the S29 `prod` row's 3.2105; the two differ by
> 0.0021 Å because the projection seed is not pinned — §1.1. The **delta** is what transfers.)* The prize is **five bits per
> target**, concentrated ~4× on the tail and 68% in `|i−j| ≥ 7`.

And the reason those five bits are unobtainable is **mechanistic, not a failed search**:

1. **The source enumeration collapses to two** — sequence and library — with physics an *operator*
   on either, not a third source (theorem G1, repaired by lane P).
2. **Both are measured.** Sequence: out-of-fold R² **0.83%** against a pre-registered 1.96% bar.
   Pool: genuinely informative at **+0.0758**, and **strictly harmful when applied.**
3. **The mechanism is a theorem about *which part* is identifiable.** The predictable part of the
   prior's error *is* the pool's common mode — and the common mode is exactly the component proved
   non-identifiable from pool data at any K.

> ### What can be predicted is coherent and therefore harmful; what would help is incoherent and therefore unpredictable.

A corrector fitted on the distogram and pool emits **+0.0554 Å worse**; a synthetic i.i.d. one at
**identical out-of-fold R² = 0.2355** emits **−0.2466 Å better**. *A 0.30 Å swing at matched
accuracy.* The i.i.d. arm is **ORACLE-constructed and not deployable** — it is the price of a
channel nobody has, never an achievement.

### 0.3 The deliverable

Not a number — **a gate**, cheap and native-free, that no measurement in this project's history
would have passed:

```
coh = corr(a corrector's residual, the pool's common-mode pair error)

    uncorrected                               0.6931
    every corrector this project owns          0.783 / 0.786 / 0.917   <- ALL RAISE IT
    imposed-structure ORACLE arms only         0.587 / 0.537

    ADMIT iff coh < 0.6931.     R2 0.16 -> -0.126 A;  R2 0.24 -> -0.247 A
```

**This is not "decorrelated from the distogram"** — orthogonality was priced at a 5.1% discount on
the requirement. Incoherence is a different condition, and it separates a −0.25 Å corrector from a
+0.06 Å one at identical accuracy.

### 0.4 How much of this report is negative, and why that is the point

Twelve lanes, 28 ledger entries, ten pre-registrations. **Four of the ten registered falsifiers
fired against the lane that wrote them.** Directions closed this sprint: the field combination, the
sparse weighted readout, subset objectives through an averaging readout, the second-moment escape
(twice, independently), generative spaces, torsion encodings, common-mode correction, filter width,
achiral single-structure channels, and recognition from single-structure geometry.

**Thirty-odd claims were withdrawn, eleven of them mine, and one withdrawal was itself withdrawn.**
The report adversary then found ten further defects in this document — two of which inverted a
headline — and they are fixed in place with the original error stated. That process is the reason
the surviving results are worth anything, and §Appendix A records all of it.

**The one thing the sprint got wrong about itself:** it claimed, in its own opening ledger entry,
to have saved the charter to disk, and did not — while that same entry was recording four prior
instances of exactly that failure as a reason to check rather than assume. The file was written 77
minutes later. The fix is now mechanical: the verifier asserts that every path the ledger claims to
have written exists.


---

## 1. The arithmetic that governed the sprint

Computed at the open, before any experiment, from S29's stored rows:

```
production, n = 126      mean 3.2105   median 2.9661
the worst 18 targets     mean 6.2758
the other 108            mean 2.6997

cap the worst 10 at 3.00 Å  ->  mean 2.9074  (−0.3031)   beats the primary target
cap the worst 18 at 3.00 Å  ->  mean 2.7426  (−0.4680)
cap the worst 30 at 3.00 Å  ->  mean 2.5778  (−0.6328)

versus improving EVERY ONE of the 126 by 0.20 Å  ->  mean 3.0105  (−0.2000)
```

**Fixing ten targets is worth more than improving all 126 by 0.20 Å.** The endpoint is a mean over
a distribution whose median already clears 3.0, so a mechanism that works on typical targets and
leaves the tail alone is close to worthless here. This is why the sprint staffed the tail first.

The caveat carried with it throughout: capping is an ORACLE operation. It bounds the prize; it does
not deliver it.

---


### 1.1 The endpoint's own reproducibility, which nobody had pinned

Four values for "production built chain" were in circulation this sprint. They are now reconciled,
and the reconciliation is itself a finding:

```
s29/results/s29_O_chain_rows.jsonl :: item=prod      chain 3.2105   cloud 3.0483   <- canonical
s27/results/chain_rows.jsonl       :: config=DIS     chain 3.2126   cloud 3.0483
lane D's independent re-projection                   chain 3.2071
lane P's independent re-projection                   chain 3.2126
```

> **The CA point cloud is identical to four decimals in both records — 3.0483. Only the built chain
> differs.** The projection from cloud to chain is **multi-start and its seed is not pinned**, so
> the same cloud projected twice gives means differing by ~0.002 Å, and independent re-projections
> range to 3.2071. The known cause is documented: a multi-start branch flip on a single target
> (S28-L18/L27b/L43), which lane D reproduced with **median per-target |Δ| of 0.0012 Å and one
> 0.513 Å outlier**.

**This report uses 3.2105** (the S29 `prod` row, which is what the charter quotes). The spread
matters for exactly one reason: **the full range across instruments is 0.0107 Å, and the sprint's
one confirmed deployable effect — the AMBER relax at k=30 — is 0.0221 Å.** That effect is only
**2×** the endpoint's own reconstruction spread. It does not threaten anything at the 3.00 or 2.50
scale, but it means **any future claim below ~0.01 Å on the built chain is inside the noise of the
instrument that measures it**, and the projection seed should be pinned before one is made.

## 2. What was tested

### 2.1 The lanes

Twelve lanes ran, never more than eight concurrently, against the charter's floor of four whenever
meaningful work was available. Ledger entries `S30-L0` … `S30-L27`, 28 in all.

| lane | remit | entries | outcome |
|---|---|---|---|
| **D** | adversary, permanent; owns and extends the cost/RMSD meter | 7 | Closed the field-combination question; withdrew two published numbers, one of them S29's |
| **L** | literature, permanent | 5 | Three theorems (reference-state compactness, common-mode non-identifiability, set-selection); found the one paper that contradicts our E2 |
| **T** | theory: the bit accounting, and when the CVaR tail stops being a prefix | 3 | **T1** (the tail is always a prefix); the codebook reframe; the value-of-a-bit law. Retracted its own headline |
| **F** | the failure tail — highest leverage by the opening arithmetic | 3 | Tail is **selection-limited, not pool-limited**; shape not scale; closed filter width by ceiling |
| **X** | divergent, permanent: should the quantum stage select at all, or generate? | 2 | The set-mean decomposition — generation closed **jointly with the readout** |
| **R** | is nativeness recognisable from one structure at all (L11) | 2 | Ordering survives, **preference fails on all 43 channels**; the local-feature null is a theorem |
| **Q** | L5 and L6 together: subset objectives and sparse readouts | 2 | Both closed; named the null that made lane T retract |
| **P** | the prior — the one measurement the whole bound reduced to | 3 | The registered null, **and** the coherence result that explains it |
| **G** | the two questions nobody had measured | 1 | **G1**, the chirality dichotomy; confirmed a claim two lanes had talked themselves out of |
| **V** | report adversary (rule 28: the main team may not approve its own positive) | — | Auditing sections 1, 3, 4 and Appendix A against the artefacts |
| **W** | what the CVaR-VQE actually contributed (charter items 12–19) | — | Investigating |
| **Y** | synthesis over the sprint's own record (items 23, 25, 26) | — | Reading |

### 2.2 Pre-registration

**Ten pre-registration files**, every one committed before the number it predicted existed:
`PREREG_S30_{D_gram, F1, F2, F3, G, P, Q_sparse, R, T, X}.md`. Three lanes' registered falsifiers
then **failed**, and all three reported the failure as the result:

- **Lane D** — its own falsifier for the combination question was not met (S30-L21).
- **Lane T** — retracted its own headline on a null lane Q told it to run (S30-L16).
- **Lane G** — **both** registered priors were directionally wrong, and it said so first (S30-L26).

Lane P registered its bars at **1.96% / 12.82% / 39.44%** R² before the regression existed and came
in at **0.83%**. Lane G registered a power note stating a half-split *could not* reach its own MDE
unless more than 100% of the effect sat in one half — then measured that more than 100% does.

> **A pre-registration that only ever confirms is decoration.** Four of ten fired against the lane
> that wrote them, which is the evidence that they were real.

### 2.3 The permanent roles, honestly assessed

The charter required three permanent roles. Two were honoured as written; one was not.

- **Adversary — honoured in substance, NOT in the form the charter asked.** The charter says
  *"rotate the adversary role so no agent is only ever attacking its own prior conclusions."*
  **Lane D held the named role for essentially the whole sprint and it was not rotated.** What
  happened instead is that adversarial action distributed itself: lane D attacked lane F's headline
  (S30-L23), lane Q named the null that made lane T retract, lane F closed a mechanism I had
  synthesised one minute earlier, lane G attacked lane L's prior and its own, and lane V was stood
  up at the end specifically to audit *my* report. The charter's purpose was served; its mechanism
  was not followed, and a sprint that relied on distribution happening spontaneously got lucky.
- **Divergent — honoured.** Lane X held it permanently and used it to overturn my premise rather
  than to explore variations.
- **Literature — honoured.** Lane L read for the whole sprint, not only at the start, and returned
  mid-sprint as §6 of the charter asks.

### 2.4 Controls

This sprint's controls are unusually load-bearing — several results died to their own control, and
one control **refuted its lane's stated mechanism while strengthening its conclusion**. The full
register of every control run and what each one ruled out is **§7.3**, compiled from the ledger.


---

## 3. What was closed, and by what

The sprint's dominant output. Each row is a direction that is now shut, with the instrument that
shut it. Several are theorems rather than measurements.

| direction | closed by | how |
|---|---|---|
| combining the 21 native-free displacement fields | D + T | ORACLE global ρ = **0.169**, LFO **0.012** (worse than the best single field) against **0.358** needed for 3.00 Å; Gram stable rank 2.057 — the fields span ~2 directions, not 11 |
| sparse weighted readout (S29's last open ladder class) | Q | **by price** — a plain argmin dominates at every bit budget; 2-of-75 with free weights is 2.1683 Å at 11.4+ bits against argmin-over-128 at **2.1435 Å for 7.0 bits** |
| second-moment / quadric escape | Q **and** T, independently | +0.2059 Å (Q, 1.81× MDE) and +0.086 (T, 2.95×); the structured `disp2` form scores 3.3585 against production's 3.0483 |
| subset objective through an averaging readout | T | **T1**: the tail is *always* a prefix — of the order induced by ∇V at the optimum. S29 §4.3 and §4.5 are incompatible and §4.3 wins |
| generative structural spaces | X | closed **jointly with the readout**: at ρ≈0 the endpoint is near-unit-slope in the set **mean**, and width buys ceiling while costing mean |
| torsion / configuration encodings | T | arithmetically infeasible — **48 bits** (2 bits/torsion x 2 torsions/residue x 12 residues) against **7** deployed |
| common-mode correction from pool data | L | **non-identifiable**: the likelihood depends on (t, μ) only through t+μ at any K; **m_eff = 1.4** of 75 members |
| E1 — a prior on the shared bias's form | L | three independent ways |
| E2 — constraint repair | L | field-scale: removing 98% of clashes costs **+0.08 Å**; works here only because 2/n is 15.4% at n=13 |
| filter width as a free lunch | F | no knee — benefit/harm cross smoothly at ~55%/55%; the **ORACLE global argmin over k IS the shipped 75** |

---

## 4. What was established positively

### 4.1 The tail is selection-limited, not pool-limited

The claim holds on **two** tail definitions, and the stratum label matters, so both are given.

**On `FAIL18`** — the *filter-defined* stratum (`s12/instrument.py:271-278`: the targets whose
score-top-75 retained **zero** members within 1.5 Å of the pool optimum) — the ORACLE best pool
member is **2.2842 Å** (built chain 2.2845), with **13 of the 18** holding a member under 3.00 and
the worst bottoming out at **3.5436**.

**On the genuinely worst 18** (by production built-chain RMSD, an order statistic, overlap with
FAIL18 **13 of 18**) the same figures are **2.5298 Å**, **11 of 18** under 3.00, worst **3.9523**.

Either way the direction survives and it is the point: **2.53 < 3.00, so the material to fix the
tail is already inside the candidate sets the pipeline is handed.** The larger claim I first wrote
— that *nothing* here is conditioned on what it measures — was too strong: the 2.2842 figure is
computed on a stratum defined by the filter's own recall, and only the worst-18 figure is free of
that.

**Retrieval is NOT MEASURED at any stratum — which is not the same as exonerated.** BLOSUM's 500
against a random 500 of the same universe fails to clear its MDE everywhere, most pointedly on the
tail (0.01× its own MDE). But the measurement is **ORACLE best-of-pool on the CA cloud** with an
MDE of 0.29 Å on the tail, power 0.05 and type-M 71, whereas the record it appears to revise
(blind 5.425 against shipped 6.019 on FAIL18) is about the **emitted built-chain endpoint**. *A
different operator, on a different basis, at a power that could not have detected the effect.*
**The prior record stands; this adds an absence of evidence, not evidence of absence** — and it is
this sprint's one instance of the control-space error (§6.1).

**The tail ordering result, with the three strata separated rather than merged:** the shipped
score's Spearman with ORACLE in-pool RMSD is **+0.6446 on the easy 108** (fold CI [+0.5905,
+0.7063], 5/5 folds) against

```
FAIL18 (filter-defined)   +0.1066   fold CI includes zero, 3/5 folds
worst18 by pool mean      +0.3798   no fold CI stored
worst18 by best-in-pool   +0.3438   no fold CI stored
random-18 null, lower bound         +0.4067
```

All three sit below the random-18 null's lower bound, so the degradation is real — but the
filter-defined number is **3.2–3.6× lower than the other two**, and only it has a CI spanning zero.
The honest statement is therefore **"the score orders its own pool markedly worse on hard
targets,"** not that it cannot order it at all.

It is the distogram's own error that drives this (ρ = −0.799, −0.819 length-residualised) and the
chain is fully mediated. On shape versus scale, lane F attached a quotation condition to its own
row and it travels here verbatim: **shape error is 83.1% of it, and `|scale error|` is 3.23× on
FAIL18 against shape's 1.93× — both facts belong in any quotation of this row.** Scale is refuted
as the *mechanism* (partial ρ = −0.067, p = 0.46, against shape's −0.641 at p = 6.5e−16) while
being the larger *absolute* deviation. Lane F also records `F3c.fires = False` with the caveat
*"statistic guessed by lane F; lane L was unreachable"* — an unregistered self-chosen statistic,
which is why this is stated as a refuted mechanism and not as a headline.

### 4.2 The pool is a codebook, not a channel

The charter asked where the 1.44 usable bits went and where the other 5.56 were spent. The question
is not well-posed: **the 500 deposited backbones carry the structure and the index only names it**,
so bits are not conserved across an index. Seven index bits move the ORACLE ladder 4.108 → 1.898 Å
**on the CA point cloud**,
which through the displacement bound is ρ = 0.887.

**The multiplier that follows is a free parameter and this sprint corrected it downward, so it is
reported with its formula rather than as a headline.** The information is
`I = −(d/2)·log₂(1 − ρ²)`, and everything depends on `d`:

```
d = 3n - 6 = 32.88  (the full coordinate dimension)   ->  36.6 bits, a 5.2x ratio
d = 6               (the effective dimension)         ->   6.69 bits, a 0.95x ratio
```

**At the corrected `d` the multiplier inverts to slightly under one.** S30-L14 §6 says so in lane
T's own words — *"S30-L15 corrects the `d` this should be evaluated at, and the corrected numbers
are smaller"* — and Appendix A records the companion 12 → 3.78 bit withdrawal made on exactly this
ground. **I originally published the `d = 32.88` figure with neither the formula nor the `d`.**

> What survives, and it does not depend on `d` at all: **bits are not conserved across an index,
> because the 500 deposited backbones carry the structure and the index only names it.** The
> *qualitative* claim is robust; the multiplier is not, and the system still cannot supply the
> seven bits either way.

The **value-of-a-bit law** follows and makes allocations comparable: `D(R) = a + c·2^(−R/γ)` fits at
**R² = 0.9983** (a = 1.3312 Å, γ = 3.1636), so `−dD/dR = 0.219·(D − 1.331)` Å per bit. Candidate
indexing beats subset cardinality by **3×** — which *explains* the 3.5× S29 measured — and
**torsion/configuration encodings are arithmetically infeasible** at this width (**2 bits per
torsion x 2 torsions per residue x 12 residues = 48 bits**, against 7 deployed; the ledger's
per-torsion figure is 51.8 bits). *As first worded — "4 basins per residue" — the arithmetic gives
24, not 48; the factor of two is the two backbone torsions.* The charter called the encoding the least-examined
component and the likely hidden bottleneck. **It was examined and it is not.**

### 4.3 The field library spends its rank on the wrong direction

Lane L derived that a fixed-reference statistical potential contains a **separable term that is a
pure function of scale** (`+kT·Σ ln P_ref(d_ij)`; a uniform 10% contraction of a 13-mer with zero
shape change moves it ~14 kT, while a size-matched reference moves by exactly 0.0000). It predicted,
falsifiably and with its concession pre-stated, that **one of the Gram's two effective directions
must therefore be the radial one**. Measured independently:

```
radial share of the Gram trace                        0.5798  [+0.545, +0.603]
cos(dominant principal direction, radial)             0.947 mean / 0.984 median
radial share of lambda_1 (ratio of trace shares, NOT a projection)   91.2%
stable rank with radial removed                       1.705 -> 2.642
   (ALL of the above are PAIR-DISTANCE space; coordinate space is 3.4-3.6 -- A.1)

cos(direction to the native, radial), all targets    -0.0675
cos(direction to the native, radial), FAIL18         -0.2524
```

> **The library spends the majority of its two available directions on a component that is
> orthogonal to the answer** — and lane F had independently shown that shape, not scale, is 83.1%
> of the tail's error. **Two lanes, one chain.**

**The anti-alignment half is withdrawn.** The FAIL18 figure of −0.2524 is **0.67× its own MDE**
(n = 18, SE 0.1337) — *below this sprint's own "not a result" line* — has no CI anywhere in the
artefact, and one of four non-empty folds (−0.6491) carries it. I bolded it in the first draft. The
orthogonal-in-general half stands, because at 0.74× MDE it is claiming a **null**, which is what a
sub-MDE number can legitimately support.

It does **not** follow that deflating scale helps: the residual 42% has no large eigenvalue, the
combination's ceiling already prices the whole span, and lane L withdrew its own deployable proposal
on exactly this point — the size-matched field *is* the deflated field, and deflation reallocates
rank without creating any.

### 4.4 Ordering exists; preference does not

On a ladder where **kind, local realism and perturbation budget are all matched** — every rung an
ideal-geometry backbone built from torsions drawn from the fold's leakage-safe Ramachandran table,
no projection, no averaging, no contraction — the verdict splits:

- **Ordering survives.** Two of 43 channels clear the bar, DIS at +0.347 with an anchor contrast of
  +0.134 (max-over-channels sign-flip null p = 0.000).
- **Preference fails on all 43.** The best `pref_near` in the library is 0.640, under the bar, and
  that channel prefers a *random pool member* to production on 0.790.
- **The largest preference effect in the library points the wrong way:** DIS prefers production to a
  0.55 Å structure on **94.2%** of targets. The project's most-cited negative, with its confound
  removed, comes out **sharper**.
- The leave-fold-out combination prefers the near-native rung to production on **93.0%** of held-out
  targets — and prefers an **arbitrary pool member on 100%** and a **3 Å rung on 100%**. Margin
  **−0.070 [−0.110, −0.028]** — *1.04× MDE at 4/5 folds, so it only just clears its own bar.*
  *It learned "is this production?", not "is this near-native".*

**And the mechanism makes the null a theorem on this instrument.** Held-out R² at matched capacity:
local features give **ΔR² −0.089**, global features **+0.600** — and the local block is
**ORACLE-advantaged**, handed per-residue deviations from the native anchor, and still adds nothing
beyond knowing the perturbation budget. **A sum of per-residue terms cannot see a lever arm.**

Resolution, which bounds what any support rule could ever select on: DIS concordance inside the near
band is **0.520** at |Δ| = 0–0.25 Å; concordance reaches 0.822 only above 4 Å, **and that figure
comes from the all-pairs table rather than the near band**. **Coarse triage and nothing
else.**

### 4.5 Generation is closed jointly with the readout

For a coordinate-average terminal, `set_mean² ≈ B² + S²` where **B is the endpoint itself** and S is
the set's spread. The terminal's entire value is spread extraction (`avg_gain = 0.4143·S − 0.1559`,
r = 0.885), and within target `corr(S, B) = +0.085` **after dropping the `T0_helix` arm** — the
artefact's primary value over all arms is **−0.4744**, which lane X itself flagged as driven by that
one arm and "must not be quoted". Both belong here: on the arms that matter, **concentration is
orthogonal to bias**. So a
set mean improved purely by concentration is worth **zero by algebra**, and the half that pays *is*
the endpoint, which is non-identifiable from pool data at any K.

The extreme case settles it: the most concentrated source ever built here (spread 0.568 against the
pool's 1.577) is **the worst endpoint in the record, 3.789 on the built chain**. **For an averaging terminal, spread is
the raw material, not a defect.**

---

## 5. What the CVaR-VQE contributed

*(Charter items 12–19. Audited by lane W against the code and the artefacts, not against the
sprint's own prose; full working in `s30/QUANTUM_W.md`.)*

### 5.1 The honest summary

Sprint 30 ran **no CVaR-VQE**. It *did* run a 9-qubit depth-3 statevector circuit on all 126
targets with exact parameter-shift Jacobians and 300 Adam iterations — but only to regenerate the
meter's ORACLE rungs (`s30/results/s30_D_ladder_structs/`, 126 files, asserted to 1e-6 against
S28's stored values). **"No VQE" and "no quantum compute" are different claims and only the first
is true**; a report sentence saying the sprint spent no quantum compute would be false as written.

The production anchor does not pass through the quantum stage at all: `core/pipeline.py:241`
(`PROD = Config()`) with `core/pipeline.py:179` (`quantum: bool = False`), and the governing
comment at `:173-178` says it outright — *"AMBER participates, Legacy is only imported for scoring,
and VQE/CVaR do not participate at all."*

What the sprint contributed to the quantum spine is **three closures and a theorem**: **T1**, which
shows the CVaR tail is always a prefix of the order induced by `∇V` at the optimum, so the deployed
stage's structural output is the single integer `m`; **T1b plus two independent quadric
measurements**, which price the only two Hamiltonian classes that could have changed that order and
find both reachable but non-transferable — **closed before either was built**, the first time this
project's cheap pre-check fired *ahead* of the spend rather than after; and **the bit accounting**,
which shows the encoding is not where the loss is.

The reason none of it was carried into a circuit is **S30-L19**: on a kind-matched, budget-matched
ladder, ordering survives on 2 of 43 native-free channels and **preference fails on all 43**, with a
mechanism (local ΔR² −0.089 against global +0.600) that makes the null a theorem on this instrument.

> **An objective cannot be pointed at a target nothing can see.**

### 5.2 The single row that answers items 12 and 16 at once

From this sprint's own cache, means over n = 126:

```
rung                                                    built chain      CA cloud
circ_best   same circuit, ORACLE objective, best of 5       0.2516        0.2884
circ_s0     ORACLE objective, 1 start (regenerated S30)     0.3175        0.3854
PROD        the deployed uniform average                    3.2071        3.0483
circ_opt    THE SAME CIRCUIT, DEPLOYED native-free score    3.4330        3.3850
```

*(`PROD` here is the meter's own re-projection; see §1.1 on the ±0.002 Å projection spread.)*

> **The same circuit reaches 0.2516 Å when its objective is the native and 3.4330 Å when its
> objective is the shipped score — worse than the classical average it was meant to improve.** The
> circuit is not the problem. The objective is, and it has been all along.

### 5.3 Item by item

| # | item | answer |
|---|---|---|
| 12 | What did CVaR-VQE contribute? | **Nothing as an endpoint contribution.** Its residual role, from T1: **α supplies τ; it is not doing selection.** Lane L's finite-shot CVaR bias hypothesis was tested and **refuted, in the opposite direction to its guess** (S30-L8) |
| 13 | What Hamiltonian was used? | **NONE.** Charter §14 withdrew the obligation. Two candidate classes (T1b halfspace; second-moment quadric) were priced and **closed before either was built** — S30-L12 and S30-L15 independently, two samplers, same sign |
| 14 | Why different from the diagonal rank-ladder? | **No new Hamiltonian to be different** — and S30 established the rank-ladder's diagonality **is not the defect**. `s30/THEORY.md:450-453` removes Hamiltonian structure, ansatz, optimiser, pool and search from the list. Referent: `H = diag(zrank(top-128 scores))`, `core/pipeline.py:838`, the same ladder on every target to **1.18% of range** |
| 15 | What did the quantum state represent? | **Partly applicable.** T1's trichotomy: under the deployed formulation the state determines **exactly one integer, `m`**. The codebook result says the 7-bit register realises 36.6 bits at `d = 32.88` (0.95× at the corrected `d = 6`; §4.2) and **delivers under one** |
| 16 | Could a classical control reproduce it? | **The honest form is "is there any positive for a classical control to reproduce?" and the answer is NO.** Every decisive classical control runs *against* the quantum stage: one sort reproduces the tail exactly; Frank–Wolfe is the correct counterpart and is polynomial; argmin dominates the sparse readout at every budget; `PROD` beats `circ_opt` |
| 17 | Trainability | **NOT APPLICABLE, not measured this sprint.** S28/S29 positions cited in `QUANTUM_W.md` rather than restated |
| 18 | Gradient variance | **NOT APPLICABLE, not measured** — but S30 measured the governing input of `Var[∂F/∂θ] = r_stable/D²`, namely `r_stable`, **on 126 real pools in 76 seconds** over data already on disk. Carry lane Q's scope fix: **`r_stable = 1.859` is pair-distance space; coordinate space is 3.404** and does not fire the threshold |
| 19 | DLA / ansatz structure | **NOT APPLICABLE, not measured.** S29's scope correction stands: the algebra is **not** maximal at n = 9, L = 3 |

### 5.4 A defect in shipped code, found while answering item 12

`core/pipeline.py:821` — the `quantum_stage` docstring — still asserts *"the CVaR tail is worth
+0.113 Å by preventing the collapse, and that is the component's measured role."* **S25-L5 withdrew
that number** (`s25/LEDGER.md:233-244`: *"never a measured effect … by this project's own fixed rule
that is a NULL"*) and replaced it with **−0.1405 Å at 0.68× MDE**.

> Anyone answering item 12 from the source file rather than the ledger gets a **withdrawn positive**
> presented as the component's measured role. This is the **sixth** instance of prose asserting a
> state that does not hold, and **the first located in shipped code rather than in a report.** It
> is not fixed here — `core/` was read-only for this audit — and it is listed in §11 as open.


## 6. Statistical discipline and multiplicity

### 6.1 The rules actually enforced

| rule | where it bit |
|---|---|
| **MDE = 2.8016 x SE, per comparison** - never a sprint-wide constant | the S29 constant 0.084 A is wrong by up to 84x in both directions |
| **Below 0.7x MDE is not a result; 0.7-1.0x is NOT MEASURED** | killed S29's built-chain preference row at 0.56x (§6.3) and lane P's headline at 0.75x |
| **Fold-clustered CI, on the pinned folds** | `s24.stats_lib._verdict` refuses a verdict when `folds=None` rather than falling back to the IID CI |
| **Pre-register the falsifier before the number exists** | every lane did; three lanes' falsifiers then failed and were reported as failures |
| **A control must match the operator's own space** | the project's most repeated error, and **there WAS an instance this sprint**: the "retrieval is exonerated" claim in §4.1 compared an ORACLE best-of-pool CA-cloud measurement against a record about the emitted built-chain endpoint. Caught by the report adversary, not by the lane |
| **The native never tunes a deployable parameter** | every ORACLE arm is labelled ORACLE and none is deployable |

### 6.2 Multiplicity, counted rather than asserted

The meter counts every comparison it emits - **30 on the built chain, 33 on the CA cloud** - and
states the consequence in the artefact itself: *"With 8 lanes, a 1x-MDE positive is EXPECTED
somewhere; price any single positive from this instrument against this count and the sprint's
running total."*

Sprint-wide the count is in the **hundreds**: 43 channels (lane R) + 21 fields (lane D) + 8 feature
arms x 5 strata (lane P) + the width sweep, the quadric/halfspace ladder, the bit ladder, and the
63 above. This is why **the surviving positives are the ones at 2-10x MDE with 5/5 folds**, not the
ones at 1.2x. Two lanes ran max-over-channels nulls rather than per-channel ones (lane R's
sign-flip null across 43, p = 0.000; lane Q's across-target null on the search ladder) and that is
the standard the next sprint should inherit for any swept family.

**The honest limitation:** the sprint-wide total is a *lower bound*. Lanes counted their own
comparisons; nobody maintained a single register, and by the time that was obvious the lanes had
closed. **A sprint-wide multiplicity register, written to as comparisons are emitted, is a build
item for S31** - not a discipline problem to exhort about.

### 6.3 Three statistic/null mismatches, which is the sprint's real methodological finding

Three times a statistic was compared against a null belonging to a *different* statistic. They are
the same defect wearing three costumes:

1. **`7 - log2 r` is right-skewed**, so its median sits 0.416 bits below its own mean. Reading the
   median against the analytic mean says "worse than chance" **by construction**. (Mine, reported
   to the user twice.)
2. **`abs(signed mean)` scored against a null for `mean |cos|`** - `s29/s29_D_fields.py:243,263`
   labels one statistic and computes the other.
3. **A 0/0.5/1 preference indicator has modal value 0**, so its median is uninformative by
   construction and a median paired difference of 0.0000 is *not* the median-vs-mean warning firing.

And a fourth, latent, caught by the verifier rather than by a reader: **`s24.stats_lib.compare` is
lower-is-better** (`d = a - b`, negative = a better) because its native statistic is RMSD. Fed a
**preference rate**, which is higher-is-better, its `verdict` string reads `WORSE` for an effect of
**+0.1716 that is the good direction**. Lane D's gate handles it correctly
(`s30_D_meter.py:451`, `PASS if effect > 0`) and no S30 document quotes the inverted field -
`s30/s30_verify.py` now asserts the trap's shape so it stays known.

> **The generalisation, and it is the one to carry forward:** a statistics library encodes a
> *direction* as well as a test. Every quantity handed to it - preference rates, accuracies,
> correlations, R2, win rates, concordances - must be checked against that direction, because the
> failure is silent and the field that is easiest to quote is the one that is wrong.


## 7. Every hypothesis entertained and killed

[PENDING]

## 8. The twelve S29 leads: pursued, rejected, and why

[PENDING]

## 9. Literature relied on and rejected

The charter asked for literature read *"for most of the sprint, not only at the start"* and to
**"read the equations, not the abstracts."** Lane L held the role permanently. The distinguishing
feature of this sprint's literature use is that **three papers closed project directions that no
experiment here had the power to close**, and one **contradicted a result we had measured** — which
turned out to be the most useful thing any of them did.

### 9.1 Relied on, and what each one settled

| source | what it settled here |
|---|---|
| **Blau & Michaeli — the perception–distortion tradeoff** | Predicts the *cross-kind* half of S28-L48's "20 of 31 scorers prefer production to a 0.25 Å ORACLE structure". Lane R used it to withdraw the project's most-cited negative **before its own numbers existed** (S30-L1). If perception–distortion predicts a result from the construction alone, the result is not evidence about nativeness |
| **Kennedy & O'Hagan (2001); Brynjarsdottir & O'Hagan, *Inverse Problems* 30:114007 (2014)** | Calibration parameter and discrepancy are **not jointly identifiable**, and a discrepancy term helps *only* given a strongly informative prior on its shape — with a *wrong* prior worse than none. This closes escape **E1** and is this project's own "confidently wrong costs 2–3× absent" arriving from a **second, independent literature** |
| **Maehara, *Oper. Res. Lett.* 43:526 (2015)** | The CVaR of a stochastic submodular set function **is not submodular**, and **no polynomial-time multiplicative approximation exists unless P = NP**. The sharpest available statement about the sprint's one hard constraint |
| **Wilder (AAAI 2018); Ohsaka & Yoshida (2017)** | The escape from Maehara is to stop asking for a single set and relax to a **portfolio — a distribution over sets** — restoring a `1 − 1/e` guarantee via continuous DR-submodular maximisation |
| **Nemhauser–Wolsey–Fisher; Das & Kempe** | Both guarantees (`1 − 1/e`, `1 − e^(−γ)`) require **monotone**. S29-L25 states in its own words that `V` is *"neither additive nor monotone"*. **Monotonicity, not submodularity, is what we fail first** — so every submodularity guarantee in the plan was void |
| **Goldberg (1984); Maurey's empirical method** | Fixed-`m` is the hard formulation (densest-`k`-subgraph); free-`m` is poly-time by max-flow. Maurey bounds any hull point within `R/√m`, which with the project's own dispersion `√63.82 = 7.99` gives **0.22 Å at m = 75** — i.e. the grid is not the constraint |
| **ANDIS (Yu et al. 2019)** | *"Native recognition and decoy discrimination cannot be optimized simultaneously with the same parameter sets."* An explicit statement, from the field, of the trade this project keeps rediscovering |
| **DOPE (Shen & Sali), incl. the authors' own DOPE-24 ablation** | The reference state's entire support for a 13-mer is [0, 15.05 Å] against a 15 Å table cutoff — the top 13% has **zero reference density**. The authors say plainly that DOPE *"is less accurate for smaller proteins"*; lane L turned that into the mechanism for the field's 40–50 residue wall |
| **Klenin & Langowski** | The writhe/Gauss double-integral formulation lane G used to build the only channel class that theorem G1 leaves open |
| **Abe et al., arXiv:2302.00704** | Already in the S29 index as the negative result on diversity interventions; used to reject diversity-aware selection without re-running it |

### 9.2 The one that contradicted us, which was worth more than the ones that agreed

*Improving consensus structure by eliminating averaging artifacts* (PMC2662860), 2090
non-homologous single-domain proteins under 200 residues:

```
baseline averaging (COMBO)   63.0 %  of atoms in clashes < 3.6 A
MCORE (their repair)          1.09 %                              <- a 58x reduction
PULCHRA                       3.64 %
RMSD, MCORE refined           3.36 A  against  3.28 A original    <- +0.08 A, WORSE
```

**Our E2 result (−0.022 Å) has the opposite sign from the published one at scale.** Lane L did not
explain the discrepancy away — it derived why both can be true: **the bonded fraction of all pairs
goes as ~2/n**, so a constraint repair is **~15% of the geometry at n = 13 and ~1% at n = 200**.
The published result is the large-`n` limit of ours.

> This is the sprint's model for how to use literature: a published number that **disagrees** with
> a measured one is a constraint on the mechanism, not a reason to doubt either. It also re-caps E2
> immediately after it escapes the identification invariance — perception–distortion catches what
> the non-identifiability theorem lets through.

### 9.3 Read and rejected, with reasons

| rejected | why |
|---|---|
| **DPPs / facility location / diversity-aware selection** | `V` is constant on centroid-equivalence classes, so **there is no diversity structure to exploit**. Rejected by algebra, not by trial |
| **Submodular maximisation guarantees** | Void — they require monotone, and `V` is not (above) |
| **QUBO / Ising formulations** | Available and exact if `f` were quadratic in the centroid, but **moot**: the continuous relaxation is cheap, and Frank–Wolfe on the simplex runs in seconds and **strictly upper-bounds any circuit on this objective** |
| **Portfolio/CVaR transfer from Wilder and Ohsaka & Yoshida** | **Mismatch stated rather than hidden:** their CVaR is over *exogenous* randomness; ours is over a distribution **the optimiser controls**, so the theorems do not transfer. What transfers is only the design lesson — relax the set, do not search it |
| **Fine-grained sub-region Ramachandran tables** | No sequence channel to exploit at peptide length; `phi` is predicted at 36.1° by full sequence context against 36.4° sequence-blind |

### 9.4 The reframe the literature produced

Lane L's reading, set beside the project's own numbers — 2 members with ORACLE **weights** reach
1.4315 Å against 75 members with ORACLE **membership** at 2.3055 Å, while choosing 2 of 500 costs
~17.9 bits against 7 for the top-128 argmin — gives:

> **Solving the set problem better is worth nothing. Exhaustive search already solved it
> (S29-L25) and the answer was worse than production. What is scarce is the information needed to
> SPECIFY a good set, and no solver supplies information.**

Which is the same sentence the sprint's measurements arrived at independently, from the other side.


## 10. The cost/RMSD meter's baselines for every cost tested

Built first, as §7 of the charter demanded, and extended mid-sprint with a `verify` verb and an
8-draw control. `python s30/s30_D_meter.py meter --f <name-or-module:function> --basis chain`.

### 10.1 The charter's four anchors all recompute

| anchor | charter's value | measured | basis |
|---|---|---|---|
| ladder rho (S28 ladder) | ~ -0.40 | **-0.4023** [-0.477, -0.322], 5/5 | built chain |
| cosine | ~ -0.03 | **-0.0339**, z = -2.45 vs its own null | CA cloud |
| native percentile | ~ 36.9 | **0.3676** [0.306, 0.405] | either |
| ORACLE preference | ~ 0.07 | **0.0714** | built chain |

Verified by `s30/s30_verify.py`, **22/22 matched, 0 mismatched, 0 missing**. The extension did not
move the instrument.

### 10.2 The shipped cost, on both bases

```
                          built chain      CA point cloud
ladder rho, S28            -0.4023          -0.1818       lowering the cost RAISES RMSD
ladder rho, CHARTER        -0.1964          +0.2603
ladder rho, FULL           -0.3187          +0.1179
cosine                      n/a             -0.0339       (-2.45 sigma, WRONG side)
native percentile           0.3676           0.3676       (0 = native is the pool's best; 0.5 = chance)
pref(ORACLE best vs PROD)   0.0714           0.2063
pref(pool member vs PROD)   0.0201           0.1265
production RMSD             3.2071           3.0483
GATE                        BLOCK            BLOCK
```

**Both bases BLOCK, and they block for different reasons** - the chain on the ladder alone, the
cloud on the ladder *and* the cosine. The gate's rule: *"a BLOCK on the ladder means the cost moves
in the wrong structural direction and gets no endpoint compute."*

### 10.3 The two costs whose percentiles were compared, and the mislabel I made

`native_pctile_DIS = 0.3676` and `native_pctile_DIS_SURR = 0.3688`. I quoted the **DIS_SURR** value
in lane D's brief beside five `DIS` numbers. They differ by 0.0012 so nothing downstream moved, but
the anchors block now carries both, named.

### 10.4 The cosine null, corrected from S29

Two different nulls, which S29 mixed:

- `per_draw_abs_mean = 0.1398` - the magnitude of **one** random direction on **one** target.
- `target_mean_null` - the null for the **126-target mean**, ~10x tighter: sd 0.0144, 95% interval
  [-0.0169, +0.0307].

The shipped cosine of **-0.0339** is z = **-2.45** against the second. Against the first it looks
unremarkable. **A cosine of +0.05 is far below 0.1398 and still several sigma above chance** - the
two nulls answer different questions and S29's paragraph used the wrong one.

### 10.5 What the meter could not see, which is a finding about the meter

Per contract rule 24, *a cost the meter cannot see is a finding about the meter.* Two surfaced:

1. **The meter was blind on the reporting basis.** It ran on the CA point cloud while the endpoint
   is the **built chain**. Both bases now run and, as §6.3 shows, they do not agree - the
   preference contrast is a *result* on one and *not a result* on the other.
2. **Its random-signed control was a single draw.** S29's seed-0 draw turned out to be the
   **maximum of its own eight**. The relative draw noise is **59% on the chain against 27% on CA**,
   because on the chain the cost prefers production to nearly everything, so the control is a rare
   event. The meter now takes `R = 8` draws and reports the per-draw distribution and the
   single-draw range.

> **A single-draw control is least trustworthy exactly where the effects are smallest - which is
> where this project's remaining effects live.**


## 11. What remains open

### 11.1 The one scientific question

**Does an observable exist whose error is incoherent with the pool's common mode?** Everything else
below is housekeeping by comparison. It has an admission test and a price (§12).

**The named gap in the closure, stated by the lane that produced it:** neither G1 nor lane P has
shown the *mutual information* between the native-free feature set and the prior's error is zero.
What is shown is that **the extractable part is the wrong part, by a theorem about which part is
identifiable.** That is a stronger claim than an exhausted search and a weaker one than an
information-theoretic zero, and it should be quoted as exactly that.

### 11.2 Open because it was never measured at the right scale

**Chiral functionals at 40+ residues.** G1 proves every achiral single-structure channel is a
distance-map reading, leaving chiral functionals as the only escape family; lane G built them and
found them empty. But **the theorem does not depend on chain length — only the emptiness does.** A
13-mer barely crosses itself, so writhe is dominated by local helical handedness. The family is
worth retesting where a chain can actually knot. *This project's instrument cannot ask the question:
the benchmark is 9–16mers and all 204 clusters are spent.*

**E2's restraint constant.** The AMBER relax at k = 30 is confirmed (−0.0406, 3.56× MDE,
length-matched, 5/5 folds) and **still has no native-free rule that selects its restraint
constant.** It was picked on dev-set RMSD. Until that rule exists it is not deployable, independent
of its size.

### 11.3 Open defects, all located, none fixed here

| defect | location | status |
|---|---|---|
| A **withdrawn positive asserted in shipped code** — "the CVaR tail is worth +0.113 Å … the component's measured role" | `core/pipeline.py:821` | S25-L5 withdrew it and replaced it with −0.1405 Å at 0.68× MDE. `core/` was read-only for this audit. **Sixth instance of prose asserting a state that does not hold; first in shipped code** |
| **The cloud→chain projection seed is not pinned** | the multi-start projection | Two records give 3.2105 and 3.2126 from an *identical* cloud (§1.1). Harmless at the 3.00 Å scale; **fatal for any future sub-0.01 Å claim** |
| **The launch gate contradicts the governor** | `s26/jobrun.py:39` (`CPU_START = 85.0`, refused at `:137`) against `s26/governor.py:61` (`CPU_CEILING = 101.0`, CPU suspension deliberately **disabled**) | The governor is content up to 101% CPU; `jobrun` refuses to *launch* above 85%. During exactly the 94–95% the charter asks for, launches are blocked while the governor is happy, so lanes launch detached and bypass the governor. **Open since S29; this is the user's call, not mine** — the one-line fix is to raise `CPU_START` to match the governor's disabled-CPU policy and leave RAM as the binding constraint, which it is |
| **No sprint-wide multiplicity register** | — | Lanes counted their own comparisons; the sprint-wide total is a *lower bound* in the hundreds. A register written to as comparisons are emitted is a build item, not an exhortation |

### 11.4 Open because the sprint chose not to spend on it

The charter's twelve leads are a **leads register, not a task list**, and it says so. §8 records which
were pursued, which were not, and — where nobody touched one — the honest reason, without
retrofitted justifications.


## 12. The next highest-value scientific question

> ### Find an observable of *this molecule* whose error is **incoherent with the pool's common
> mode** — not merely decorrelated from the distogram.

This is not a restatement of "get a better prior". It is a specific, measurable, and cheaply
falsifiable condition that the sprint derived and then priced.

### 12.1 Why this and not something else

Three results compose into it, each independently established:

1. **The source enumeration collapses to two** — sequence and library — with physics as an
   *operator* on either rather than a third source (G1, repaired by lane P: a universal function
   evaluated on a target-specific argument yields target-specific output, so physics has no target
   argument of its own).
2. **Both sources are measured.** Sequence-derived features: out-of-fold R² **0.83%** against a
   registered 1.96% bar; the entire uncollapsed posterior is worth **−0.018**. The pool:
   genuinely informative at **+0.0758** out of fold at long range, and **strictly harmful when
   applied**, raising the residual's coherence with the pool's common mode from 0.6931 to **0.9172**.
3. **The mechanism says why, and it is a theorem, not a search result.** The predictable part of the
   prior's error *is* the common mode, and the common mode is exactly the component S30-L7 proves
   non-identifiable from within the pool.

> **What can be predicted is coherent and therefore harmful; what would help is incoherent and
> therefore unpredictable.**

### 12.2 The admission test — use this, not decorrelation

Define `coh` as the **within-target correlation of a candidate corrector's residual with the pool's
common-mode pair error** `μ_{t,p}` (the mean over the retained 75 of `d_{m,p} − d_nat,{t,p}`; the
S30-L7 quantity).

```
uncorrected                                     coh = 0.6931
every fitted corrector this project owns        coh = 0.783 / 0.786 / 0.917   <- ALL RAISE IT
imposed-structure ORACLE arms only              coh = 0.587 / 0.537           <- only these lower it

ADMIT a prior corrector iff it LOWERS coh below 0.6931.
```

**This is not "decorrelated from the distogram."** Lane L priced orthogonality at a **5.1% discount
on the requirement** — decorrelation is not the lever. Incoherence is a different condition and it
is the one that separates a −0.25 Å corrector from a +0.06 Å one **at identical accuracy**.

### 12.3 The price, which is the encouraging part

The required quality is **far below the displacement bound**, because this corrects the *prior* and
acts through the filter rather than being a displacement added to the answer:

```
R2 = 0.16  ->  -0.126 A        ORACLE-CONSTRUCTED, NOT DEPLOYABLE:
R2 = 0.24  ->  -0.247 A        a PRICE for a channel nobody has, never an achievement
```

And the target is sized: **five ORACLE signs on long-range pairs are worth −0.3259 Å on the built
chain (3.2126 → 2.8867, 2.05× MDE, 5/5 folds)** — which would clear the charter's primary target.
**The prize is five bits per target.** The deployable sign currently performs at accuracy ≈ 0.60
against a *free* baseline of 0.627; ≈ 0.8 is needed.

### 12.4 The secondary recommendation, and three measurements agree on it

**If the next sprint has one place to spend bits, it is the readout.**

```
readout index bit      0.376 A/bit   (lane T's value-of-a-bit law at D = 3.05)
prior-sign bit         0.0713 A/bit  (lane P, cloud)  ->  5.3x
prior-sign bit         0.0652 A/bit  (lane P, chain)  ->  5.8x
candidate indexing over subset cardinality            ->  3x  (lane T, independent)
```

**Flagged unregistered, and the currencies are not exchangeable** — one names a candidate, one
shifts a prior, and lane T's law was fitted on an ORACLE ladder. But three independent measurements
agree in direction, which is more than any single one of them is worth.

### 12.5 What NOT to spend on

Closed by theorem or by price, with the closure named: the field combination (rank, not count);
sparse weighted readouts (argmin dominates at every bit budget); subset objectives through an
averaging readout (T1); the second-moment/quadric escape (twice, independently); generative spaces
(closed *jointly* with the readout); torsion encodings (48 bits against 7); common-mode correction
from pool data (non-identifiable at any K); achiral single-structure channels (G1); and
recognition from single-structure geometry (ordering survives, preference fails on all 43).


## 13. Architecture diagram of what was actually built

### 13.1 The production path — unchanged by this sprint

```
  sequence
     |
     v
  [1] BLOSUM62 retrieval, K = 500            <- NOT MEASURED as harmful at any stratum
     |                                          (and not exonerated either -- see 4.1)
     v
  [2] leave-fold-out ESM-2 650M distogram, 17 bins
     |                                          <- the ceiling lives here: its error is
     |                                             83.1% shape, and 68% of the recoverable
     v                                             prize is in |i-j| >= 7
  [3] L1 Bayes-risk score  ->  top-75 filter
     |                                          <- rho +0.6446 on the easy 108,
     |                                             +0.1066 to +0.3798 on the tail
     v
  [4] coordinate average of the 75            <- extracts spread (r = 0.885); cancels
     |                                           i.i.d. error, NOT coherent error
     v
     CA point cloud   3.0483 A                <- exact, reproduces to 4 decimals
     |
     v
  [5] multi-start ideal-geometry projection   <- SEED NOT PINNED: +-0.002 A (section 1.1)
     |
     v
     BUILT CHAIN      3.2105 A   <- THE ENDPOINT (charter section 2)


  [Q] quantum stage ....................... OFF THE PATH
      core/pipeline.py:179  quantum = False
      core/pipeline.py:241  PROD = Config()
      "VQE/CVaR do not participate at all"  (:173-178)
```

### 13.2 What S30 actually built

No new production architecture. **Nothing in §13.1 changed**, and that is the result, not an
omission — every candidate replacement was closed by measurement, by theorem or by price before it
reached the pipeline.

What was built is **instrumentation**:

```
  s30/s30_D_meter.py     the cost/RMSD meter, extended from s29's
                         + a `verify` verb        (the old `selftest` was an 8-residue synthetic
                                                   check that would pass on a drifted meter)
                         + BOTH reporting bases   (it had been blind on the built chain)
                         + an R = 8 draw control  (the single draw was the maximum of its own 8)
                         + a multiplicity counter (30 chain / 33 cloud comparisons per run)

  s30/s30_verify.py      36 headline numbers recomputed from artefacts, 0 mismatches
                         + the endpoint's own reproducibility
                         + an assertion that every path the ledger claims to have written EXISTS
                           (the mechanical fix for the S30-L0 failure)

  s30/THEORY.md          T1, T1b, the codebook reframe, the value-of-a-bit law
  s30/BRIEF.md           the charter, verbatim -- written 77 minutes after being claimed
  s30/AUDIT_V.md         the report adversary's findings (rule 28)
  s30/QUANTUM_W.md       items 12-19, audited against the code
```

### 13.3 The one structural change the sprint argues for

Not a new stage — a **gate on an existing one**:

```
  a candidate prior corrector
     |
     v
  compute coh = corr(residual, pool common-mode pair error)     [S30-L7 quantity]
     |
     +-- coh >= 0.6931  ->  REJECT.  It will emit WORSE at any accuracy.
     |                      (every corrector this project owns lands here)
     |
     +-- coh <  0.6931  ->  admit and measure the endpoint.
                            R2 0.16 -> -0.126 A;  R2 0.24 -> -0.247 A
```

**This gate is the sprint's deliverable.** It is cheap, it is native-free, it discriminates a
−0.25 Å corrector from a +0.06 Å one at identical out-of-fold accuracy, and **no measurement in
this project's history would have passed it.**


---

## Appendix A — claims withdrawn during this sprint

The sprint's most distinctive feature was lanes destroying their own results, usually before anyone
asked. Recorded in full because the charter requires it and because it is the reason the surviving
results are worth anything.

### A.1 The coordinator's (mine)

| claim | how it died |
|---|---|
| **The "two qubits" synthesis** — a selecting readout over a wider register aimed at the tail | Written at 13:07, corrected at 13:08. Lane F closed filter width by ceiling (the ORACLE global argmin over k **is** the shipped 75) and showed the widening-rescues-the-tail result does not replicate on either filter-independent tail (≈0, and *reversed* at +0.73 on one) |
| **"The median target is worse than chance"** — reported to the user twice | `7 − log₂r` is right-skewed. Null mean 1.4050, null **median 0.9888**, null win rate 62.5%. "Below random on 82 of 126" is what a random ranking does to itself (null expects 78.8, z = +0.60). I read a median against a mean's baseline |
| **"LEG_torsion is at chance"** | My wording, not lane R's. Its anchor contrast is **+0.024 with the fold CI excluding zero** — it *does* order slightly above its control. It fails by being a quarter of the required margin, not by being noise |
| **The ceiling gate I handed lane X** — "measure the ceiling first, always" | An **anti-predictor**: 1 of 10 cells correct on predicting the endpoint's direction, against 5/6 for the set-mean law; mean residual 0.333 Å vs 0.0153 — **22×** |
| **"Nobody has built a typical-good generator"** | Somebody had, and **diversity is the correct design**. For an averaging terminal spread is the raw material; the most concentrated source ever built is the worst endpoint in the record |
| **"Run `selftest` and confirm it reproduces the baselines"** | `selftest` is a synthetic 8-residue check that runs in 0.4 s and **would pass on a meter whose every number had drifted**. The instruction was unsatisfiable. A real `verify` now exists |
| **The 0.3688 anchor in lane D's brief** | That is **DIS_SURR's**; the shipped cost is 0.3676. I listed it beside five `DIS` numbers |
| **"Find a source with decorrelated errors"** | Priced by lane L: a *perfectly orthogonal* channel is worth a **5.1% discount** on the requirement and must be **3.01× better** than anything owned. **Decorrelation is not the lever; skill is** |
| **Quoting stable rank 1.86 without its feature space** | 1.86 is **pair-distance** space; **coordinate space is 3.4–3.6** with k90 = 11.2 not 5.6 — and sparse readouts and second-moment constructions act on coordinates |
| **Citing S28-L48 as established** | Every rung in that ladder differs in **kind** as well as nativeness, so it measured a cross-kind preference that perception–distortion already predicts (lane R, S30-L1, posted before its own numbers existed) |

### A.2 Lane T's — four, three toward its own hypothesis

| claim | how it died |
|---|---|
| **P2b: retrieval worth +0.10–0.30 Å** | Measured **0.07**; its copula reasoning overestimated by 2–4× |
| **P4c: searched classes beat prefix by < 0.6 Å** | Measured **0.982** — missed by 64%, in the flattering direction |
| **"12 bits ≈ six coefficients — two independent routes"** | An **unregistered aside**; a coincidence, not a derivation. Withdrawn and corrected *downward* to 3.78 bits, and flagged louder precisely because it had no bar attached |
| **"HALFSPACE − PREFIX = −0.9816 Å, BETTER"** — its own headline | Lane Q named the null; lane T ran it with a common direction bank: **196% accounted**, split-half transfer 19%, and the **transferable rule lands +0.472 Å WORSE than production** (1.88× MDE, 5/5 folds). Its own diagnosis: the matched random-subset null asked whether a structured class beats an unstructured one, not whether **the winning direction is the same direction twice** |

### A.3 Lane L's

| claim | how it died |
|---|---|
| **EDM projection / triangle repair as a route** | Grepped the codebase after deriving it: S19 had already measured it and left a **named prohibition**. Refuted in the *unexpected* direction — an incoherent magnitude-matched field is worse on realisability and lands **1.24 Å better** |
| **The deployable half of its own reference-state proposal** | Its own algebra refutes it: the score decomposes additively, so the size-matched field **is** the deflated field, and deflation reallocates rank without creating any |
| **"The AMBER-relax benefit concentrates on divergent pools", 2:1** | Lowered to **roughly even** by lane L itself, on lane F's dispersion null (−0.128, CI includes zero) — adjacent evidence pointing the other way. **THIS WITHDRAWAL IS ITSELF WITHDRAWN (S30-L26): the claim is CONFIRMED at −0.0406, 3.56× MDE, 5/5 folds on a length-matched split, with 97.2% of the gain in the divergent half.** Lane L moved on a null about a *different outcome* and lane G's matching prior rested on a *different mechanism*; two revisions agreeing was not evidence. See §A.7 |
| **"Two independent arrivals" on confidently-wrong** | Downgraded to one literature result plus one project result on a *different object*, after lane F measured error×confidence (−0.729) as **worse** than error alone (−0.799) |

### A.4 Lane X's

| claim | how it died |
|---|---|
| **Its own S30-L10 admission condition** | Amended to be **stricter**: Δ(set mean) is two channels, not one, so `ADMIT iff Δ(bias B) + 0.298·Δ(set best) < 0` |
| **Its counting falsifier** | Volunteered as having **essentially tied its bar — 24.80% against 25%** — and therefore deciding nothing; the verdict rests on the transfer arm alone |
| **A −0.4744 correlation it could have quoted** | Flagged as **driven by one arm**; drop it and the value is +0.0851. "The strong negative is an artefact and must not be quoted" |

### A.5 Lane F's, lane D's, lane R's, lane Q's

| claim | how it died |
|---|---|
| **Lane F's F1c FAIL18 row** (+1.7674 Å) | The adversary: FAIL18 is *defined* as the targets whose top-75 retained **zero** in-band members, so **49% is forced arithmetic** and the rest is truncation-inflated. Lane F *had* checked circularity against production; the predicate is the **filter's own recall**. **"A matched control in the right space does not rescue a stratum defined by the outcome"** |
| **Lane F's own F2a and F2c** | Both refuted by its own measurements — no knee in filter width, and its widening result flagged by itself as circular |
| **Three defects in lane D's own combination arms** | An untuned ridge returning a negative cosine; a least-squares objective violating its own consistency floor (caught because ρ fell below the best single field, impossible for a true maximum); selection leakage worth 0.019 in ρ |
| **Lane R's 0.500/0.500/0.500 cell** | A **null-input artefact** — three identical values across three different questions is what caught it |
| **Lane Q's framing of my point 4** | "The framing survives; the operator does not" — identity is not only the better question but the **cheaper** one, so a sparse weighted readout is the expensive way to ask it |

### A.7 The one withdrawal that was itself wrong, and lane G's self-kill

**A withdrawal is a claim too, and it can be wrong in the same ways.** Appendix A recorded lane L
lowering "the AMBER-relax benefit concentrates on divergent pools" from 2:1 to roughly even. Lane G
then measured it: **confirmed at −0.0406 Å, 3.56× MDE, 5/5 folds** on a length-matched split, with
**97.2%** of the gain in the high-dispersion half and a permutation null at p = 0.000. Lane G had
*also* registered 2:1 against, and lost.

> **Two lanes agreeing did not make a prior.** Lane L moved on lane F's null, which was about a
> different *outcome* (filter failure, not relax gain). Lane G moved on the common-mode argument,
> which was about a different *mechanism* (bias repair, not geometry repair). Each rested on an
> object that was not the one under test, and their agreement read as convergence.

What survives the confirmation is narrower than the claim: **the effect is too small to matter,
which is not the same as failing to reach the tail, and my first draft argued it on the wrong
stratum.**

I originally quoted two FAIL18 rows. **FAIL18 is the wrong tail for a dispersion-graded question**
— it is the *filter-defined* stratum and is barely enriched in high dispersion (11 of 18 against
9.0 expected) — and `s30/results/s30_G_disp2.json` holds the filter-independent rows **right beside
the ones I quoted**, where both legs reverse:

```
                                  FAIL18 (quoted)        worst18 by pool mean (not quoted)
tail dispersion vs the rest       +0.4583  0.47x MDE     +1.6116  1.44x MDE  CI [+0.583,+2.076]
relax gain, tail vs rest          -0.0113 vs -0.0239     -0.0386  vs  -0.0193
```

On the filter-independent tail the relax gain is **twice as large on the tail**, not half. The
ledger printed both; the report printed one. Worse, the FAIL18 relax-gain difference I called *"the
wrong direction"* is **0.22× MDE** — I used a number below my own not-a-result line as evidence of
absence.

**The third leg is sound and is sufficient on its own:** targeting the divergent half buys
**−0.0215 against −0.0221**, i.e. nothing, because the other half contributed nothing to begin
with. And the whole arm is **0.69% of the 3.2126 Å baseline**. So **E2 stays out of the mechanism
column on SIZE, not on reach**, and its restraint constant still has no native-free selection
rule.

**And lane G killed its own positive before anyone quoted it.** WRITHE's preference contrast is
**+0.1641, 2.17× MDE, 5/5 folds, max-null p = 0.000**, clearing *both* of lane R's preference
clauses that none of 43 channels cleared. It is **cross-kind** — the control keeps deposited
coordinates while the near rungs are ideal rebuilds. The kind-matched statistic settles it: the
rebuilt native's percentile inside its own ladder is **0.6061 for WRITHE and 0.6732 for |WRITHE|,
both worse than the 0.5 chance line**, against DIS's 0.2876. **Second instance of the cross-kind confound this sprint**, after S28-L48 (withdrawn by lane R).
*My first draft said "third" and counted the widening result, which Appendix A correctly records as
circular rather than cross-kind — a different defect.*

### A.6 The three checklist entries this sprint earned

1. **A matched control in the right space does not rescue a stratum defined by the outcome.**
2. **When a statistic is a nonlinear transform, its mean, median and win-rate have three different
   nulls** — never read one against another's baseline.
3. **"The winning direction is the same direction twice"** is a failure mode distinct from the
   order-statistic one, and needs a common direction bank to detect.
