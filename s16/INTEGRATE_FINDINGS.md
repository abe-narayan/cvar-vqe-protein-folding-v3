# SPRINT 16 / THE INTEGRATED ARCHITECTURE — AND LEGACY vs AMBER INSIDE IT

**Modules** `s16/integrate.py`, `s16/diversity.py` · **artefacts** `s16/results/integrate.json`,
`s16/results/diversity.json` · **logs** the matching `.log` files
**Instrument**: the 9 fully enumerated targets, 3 seeds, budget 6144 objective evaluations

> ## CORRECTED 2026-09-06 AFTER VERIFY RECOMPUTED WITH THE TARGET AS THE UNIT
>
> The caveat below was written before VERIFY reported. **VERIFY has now reported and the correction
> is worse than the caveat implied.** Both are kept.
>
> - **The AMBER filter does not survive.** With the target as the unit: **-0.054 [-0.148, +0.040]**,
>   CI **2.79x wider**, **4 targets worse / 5 better**, **dropping any of 8 of 9 targets puts the CI
>   over zero**, and **seed 2 alone gives -0.0003**. No per-generator CI excludes zero. The
>   description "the only positive accuracy effect in Sprint 16" is **WITHDRAWN** for this arm. The
>   Legacy-vs-AMBER *contrast* is unaffected in direction but inherits the same width.
>   **Note added after the ENERGY workstream reported:** AMBER's effect on the SHIPPED 126-target
>   ensemble is a different measurement and it does survive in aggregate -- **-0.0234 [-0.0375,
>   -0.0099], n = 123 gated** (`s16/energy_ablate.py`, report section 5.7). That is the sprint's one
>   surviving accuracy effect. It is not this one: that arm is AMBER as a **relaxation operator** on
>   the shipped ensemble, this arm is AMBER as a **ranking filter** on VQE-generated ensembles, and
>   the two must not be conflated.
> - **The diversity attribution in Section 4 is wrong by about a factor of two.** The table prints
>   the **free-superposition** member error (spread 0.056 A) while the identity consumes the
>   **common-frame** term (spread **0.122 A**). Exact per-target attribution: **cvar0.25 is 52%
>   member error, cvar1.0 51%, cvar0.5 42%.** It is roughly half diversity and half conformer
>   quality, not diversity alone.
> - **The tie trap does NOT apply** — censused, zero ties in Legacy or AMBER.
> - **Section 5 survives at target level**: coordinate-over-torsion readout **+0.446 [+0.193, +0.706],
>   8/9 targets**. It is the only accuracy claim in this file that holds.
>
> **Original caveat, as written before VERIFY reported:** the bootstrap in `s16/integrate.py`
> resamples **rows**, and a row is one (target, seed, generator) cell — so 162 rows sit on only
> **9 targets**. The programme's own rule is that the TARGET is the unit of analysis. These
> intervals are therefore **too narrow**.

---

## 1. What this experiment is for

It is the only place in the sprint where all four mandated pillars meet on the same targets, the
same ensembles and the same readout. Legacy and AMBER are made **explicitly comparable** by being
applied as *the same operation* — drop the worst quartile of the ensemble by that energy — so the
only thing that differs between them is which energy does the ranking. Comparing an energy to an
*operation* (a filter against a minimisation) would not be a comparison of the energies at all.

    CVaR-VQE ensemble → structure-level readout → Legacy detection → AMBER repair

- **Genuine VQE and genuine CVaR-VQE** are arms of one ladder: `α = 1.0` is plain VQE, `α < 1` is
  CVaR-VQE. Real MPS ansatz, shot sampling, Adam, a real qubit Hamiltonian, a hard
  objective-evaluation budget.
- **Genuine Legacy**: the eleven-component model at `DEFAULT_WEIGHTS`, never fitted.
- **Genuine AMBER**: ff14SB/GBn2 single points computed **on demand**, because the precomputed
  table covers 1.1% of the enumerated space and a sampled ensemble essentially never intersects it.

---

## 2. Legacy versus AMBER, the same operation on the same ensembles

Each energy priced against **dropping the same number of ensemble members at random** — the control
a filter has to beat.

| m | contrast | mean | [95% CI]* | median | W/L |
|---|---|---|---|---|---|
| 5 | AMBER − random | −0.024 | [−0.074, +0.027] | −0.009 | 85/75 |
| 5 | **Legacy − random** | **+0.068** | **[+0.012, +0.125]** | +0.036 | 70/90 |
| 5 | AMBER − Legacy | −0.092 | [−0.134, −0.050] | −0.068 | 92/52 |
| 20 | **AMBER − random** | **−0.051** | **[−0.090, −0.011]** | −0.045 | 101/61 |
| 20 | Legacy − random | −0.003 | [−0.048, +0.044] | −0.002 | 81/81 |
| 75 | **AMBER − random** | **−0.054** | **[−0.088, −0.020]** | −0.061 | 101/61 |
| 75 | Legacy − random | +0.005 | [−0.034, +0.044] | −0.025 | 85/77 |
| 75 | **AMBER − Legacy** | **−0.059** | **[−0.089, −0.030]** | −0.043 | 97/65 |
| 75 | Legacy→AMBER − random | −0.053 | [−0.090, −0.017] | −0.051 | 99/63 |

\* provisional — see the caveat at the head of this file.

**The answers the sprint was asked for:**

- ~~**AMBER contributes.** As a ranker inside this architecture it is worth about **−0.05 Å**
  against a matched random drop at m = 20 and m = 75, with intervals excluding zero and win/loss
  above 60%. It is the only positive accuracy effect anywhere in Sprint 16.~~
  **[CORRECTED 2026-09-06.]** At target level the point estimate is unchanged (−0.054) but the
  interval is **[−0.148, +0.040]**, 4 targets worse and 5 better, and the effect is destroyed by
  dropping any of 8 of the 9 targets. **AMBER's contribution AS A RANKING FILTER on these ensembles
  is not established.** Its contribution as a **relaxation operator** on the shipped 126-target
  ensemble is a separate measurement and does survive in aggregate (−0.0234 [−0.0375, −0.0099]);
  see report section 5.7.
- **Legacy does not.** At m = 20 and m = 75 it is indistinguishable from dropping members at
  random. At m = 5 it is **significantly worse than random** (+0.068 [+0.012, +0.125], 70W/90L) —
  on a small ensemble, ranking by Legacy is worse than not ranking at all.
- **AMBER beats Legacy directly**, −0.059 [−0.089, −0.030] at m = 75, 97W/65L.
- **The composition adds nothing over AMBER alone.** Legacy→AMBER is −0.053 against random where
  AMBER alone is −0.054. The interaction term is −0.004 [−0.042, +0.032] at m = 75. Legacy
  contributes nothing to the pipeline either alone or in front of AMBER.
- This is consistent with, and sharpens, the recorded result that Legacy in GENERATION hurts while
  Legacy as a corrector of the selector was only −0.16 Å at p = 0.23.

**Both energies at every m are dwarfed by the generator choice** (§3), which moves 0.30 Å.

---

## 3. What the quantum component did

Every CVaR-VQE arm at m = 75, against its two mandatory controls at matched budget:

| arm | minus untrained-circuit best-of-N | minus uniform random |
|---|---|---|
| `cvar0.1` | +0.115 [−0.054, +0.287], 12W/15L | **+0.212 [+0.059, +0.365]**, 8W/19L |
| `cvar0.25` | **+0.204 [+0.038, +0.366]**, 9W/18L | **+0.301 [+0.163, +0.437]**, 5W/22L |
| `cvar0.5` | +0.137 [−0.050, +0.317], 11W/16L | **+0.234 [+0.087, +0.382]**, 8W/19L |
| `cvar1.0` (plain VQE) | +0.112 [−0.062, +0.281], 9W/18L | **+0.209 [+0.046, +0.357]**, 8W/19L |

Positive means the VQE is **worse**. Every arm loses to uniform random sampling at the same budget,
with intervals excluding zero; against the untrained-circuit control the losses are the same size
but only `cvar0.25`'s interval excludes zero.

**This reconfirms the standing result — running the VQE is worse than not running it — on a
readout it had never been tested against.** The earlier finding was measured through an argmin;
this is an ensemble → coordinate-average readout, a different terminal operator, and the conclusion
survives the change. **α does not rescue it**: no monotone trend in α, and plain VQE (α = 1.0) is
among the better arms.

---

## 4. WHY — and this is the part that is new

`s16/diversity.py`. At m = 75 the six generators have **almost identical mean member error**:

| generator | readout | mean member (quadratic) | diversity | mean pairwise RMSD |
|---|---|---|---|---|
| `uniform` | **3.214** | 3.966 | **2.626** | 3.346 |
| `bestofN` | 3.311 | 3.966 | 2.501 | 3.217 |
| `cvar1.0` | 3.423 | 3.983 | 2.475 | 3.075 |
| `cvar0.1` | 3.426 | 3.950 | 2.393 | 2.947 |
| `cvar0.5` | 3.448 | 3.971 | 2.425 | 3.001 |
| `cvar0.25` | 3.515 | 4.006 | 2.418 | 2.966 |

The mean member error spans **0.056 Å**. The readout spans **0.301 Å**. The difference is
**diversity**, and the ordering of the diversity column reproduces the ordering of the readout
column almost exactly.

That is not a coincidence, and it is not a hypothesis either — it is arithmetic. The
Krogh–Vedelsby ambiguity decomposition,

    ‖avg − native‖² = mean_i ‖x_i − native‖² − mean_i ‖x_i − avg‖²

holds through this operator to **machine precision**: residual −0.0000 [−0.0000, +0.0000] with
max |residual| 0.0000 at m = 5, 20 and 75, measured in the frame the operator itself builds in.

**Stated honestly: panel A is a theorem check, not a discovery.** It confirms the implementation is
correct and nothing more — the RETRACT workstream established independently, and at machine
epsilon, that this identity is exact rather than predictive, and that Sprint 15's "validated to
0.162 Å" residual was 54% an arithmetic-mean error and 46% a frame-convention mismatch, with 0%
attributable to the identity failing.

**The content is the measurement, not the identity.** Because the decomposition is exact, saying
which of its two terms the generators differ in is a complete causal account.

~~And the answer is: **they differ in diversity, not in conformer quality.** The VQE does not produce
worse conformers — it produces *less various* ones, and the terminal operator's entire gain is the
diversity term it destroys.~~ **[CORRECTED 2026-09-06 by VERIFY — wrong by about a factor of two.]**

The `mean member (quadratic)` column above is the **free-superposition** member error, which is not
the term the identity consumes; the identity uses the **common-frame** term, whose spread across
generators is **0.122 Å**, not 0.056 Å. Exact per-target attribution of each generator's readout
penalty:

| generator | share attributable to member error | share attributable to diversity |
|---|---|---|
| `cvar0.25` | **52%** | 48% |
| `cvar1.0` | **51%** | 49% |
| `cvar0.5` | **42%** | 58% |

**The corrected statement: the VQE's readout penalty is roughly half concentration and half worse
conformers.** Diversity is the larger single term on one arm and the smaller on two. Concentration
is a real and substantial part of the mechanism; it is not the whole of it, and the original sentence
overstated it about twofold. This is the same failure mode VERIFY names throughout the sprint —
a quantity measured correctly, then read as if it were the quantity a law consumes.

~~This also **refines a recorded project law.** The record says the terminal operator consumes the
set MEAN rather than the set BEST. Here the set means are equal to 0.056 Å and the readouts differ
by 0.301 Å, so the mean does not account for it; the missing term is diversity, and with it the
account is exact.~~ **[CORRECTED 2026-09-06.]** The 0.056 Å was the wrong term. In the frame the
identity uses, the member-error spread is 0.122 Å and accounts for roughly half the readout spread.
**The recorded law is not overturned — it is incomplete**: the operator consumes the set mean *and*
the set diversity, in equal and opposite weights, and the earlier claim that the mean plays no part
here was an artefact of reading the wrong frame.

**And it names the mechanism behind "concentration is wrong when discrimination binds":**
concentration is precisely the destruction of the diversity term. That is why searching harder on
this objective loses, and why CVaR's α reads as a diversity dial rather than a tail parameter.

**One negative to keep.** Diversity is native-free, so it is a candidate generator-selection
criterion — but its rank correlation with the readout is only −0.145 to −0.270 within target, far
weaker than the exact identity might suggest, because the readout depends on *both* terms and the
member-error term varies too. It orders the six generators; it does not order individual runs.
RETRACT independently found the related native-free screen (structural disagreement `s` as an
a-priori fusion screen) **refuted** at AUC 0.401 against a 0.500 null. Treat native-free diversity
as a weak generator-level heuristic, not as a selection signal.

---

## 5. The architecture's own design claim, priced

The pivot to a coordinate-space readout is a claim, so it is measured rather than assumed:

| m | torsion readout minus coordinate readout | W/L for coordinate |
|---|---|---|
| 5 | +0.214 [+0.149, +0.284] | 104/58 |
| 20 | +0.297 [+0.218, +0.376] | 118/44 |
| 75 | **+0.446 [+0.356, +0.537]** | **124/38** |

Coordinate readout wins at every ensemble size and **the advantage grows with ensemble size** —
which is exactly what the ambiguity decomposition predicts, since more members means more diversity
to cancel, and torsion averaging cannot cash it. This is the same geometry the flagship measured
from the opposite direction, and it is the one architectural choice in Sprint 16 that survives.

---

## 6. The one constructive experiment the mechanism suggested — and its clean null

**Module** `s16/divselect.py` · **artefact** `s16/results/divselect.json` · **n = 126, the real
instrument** (not the 9 enumerated targets §1–5 use)

If the terminal operator's gain *is* the diversity term, then the incumbent — which selects its 75
windows by distogram score alone — is optimising one of the identity's two terms and ignoring the
other, which enters with equal weight and opposite sign. So: select for score **and** diversity.

The falsifier was written before the run: *if no (m, λ) beats the shipped top-75 out of fold with an
interval excluding zero, the mechanism is exact arithmetic without a lever.*

**Instrument check first.** The greedy selector at λ = 0 reproduces the shipped top-75 exactly
(3.048 Å, +0.000), so the arms differ only in selection, never in readout.

**The trade is real, measurable, and cancels.** As λ rises at m = 75:

| λ | 0.0 | 0.05 | 0.2 | 0.35 | 0.5 | 0.8 | 1.5 |
|---|---|---|---|---|---|---|---|
| diversity (native-free) | 1.902 | 1.916 | 1.959 | 1.996 | 2.039 | 2.126 | **2.316** |
| member error (ORACLE) | 3.615 | 3.618 | 3.637 | 3.649 | 3.669 | 3.694 | **3.778** |
| readout RMSD | 3.048 | 3.045 | 3.049 | 3.040 | 3.049 | 3.040 | 3.042 |

Diversity rises 22%, member error rises 4.5%, and **the readout moves by at most 0.008 Å**. The two
terms rise together and cancel, exactly as an identity with equal and opposite weights requires.
The experiment is informative *because* both terms are printed: this is not "diversity failed to
rise", it is "diversity rose and was paid for".

**The deployable arm is a null.** (m, λ) chosen leave-fold-out: **+0.015 [−0.019, +0.049]**, median
+0.000, **61W/65L**. And the module's own pre-registered trap fires — three of five folds picked
λ = 1.5, the **end of the ladder**. A selector pinned at a boundary is choosing an endpoint rather
than a trade, which is the same failure `s16/csteer.py` caught; it was flagged in advance here and
it still delivered nothing out of fold.

**One zero-information result worth keeping, with its caveat stated in the same breath.** Selecting
75 windows by **pure diversity, ignoring the distogram score entirely** (`maxdiv75`) gives 3.051 Å
against the shipped score-based 3.048 Å — indistinguishable — while a matched random selection gives
3.090 Å. **Caveat that must travel with it:** the candidate set is the top-200 of the K = 500 pool
*by score*, so `maxdiv` is not score-blind — it inherits that pre-filter. The correct statement is
narrow and is the only one supportable: **within the top-200 by score, ranking the remaining
candidates by score adds nothing over ranking them by diversity.** The score's selective work is
done in the 500 → 200 cut, not in the 200 → 75 one. That is consistent with the recorded result that
top-24 filtering caps any reranker, and with the shrinking m\* the record already tracks.

**Verdict.** The Krogh–Vedelsby account of the terminal operator is exact and it explains the
quantum result completely (§4). It does **not** yield a selection lever on the real instrument. The
sprint reports a mechanism without a lever, which is the honest description.
