# PREREG — S31 LANE F — THE TERMINAL OPERATOR

**Committed before the first number exists.** Contract rule 19. Lane F, 2026-09-20.

Endpoint basis for every falsifier below: **mean built-chain Cα RMSD over the 126 `tuning126`
targets**, production incumbent **3.2105 Å**, re-projected in the *same process* as the test arm so
that contract rule 31 (unpinned projection seed, 0.0107 Å instrument spread) cannot decide a
comparison. The CA point cloud (production 3.0483 Å) and the *set mean* (3.5507 Å) are named
explicitly wherever they appear and are never the endpoint.

---

## 1. The two operators, stated exactly

Both consume the **same** input: the shipped top-75 of the 500-member BLOSUM pool, ordered by the
shipped distogram Bayes-risk score (`s12.instrument.shipped_score`), with the production
reproduction gate (recomputed top-75 must equal the record's `sub` set-wise on all 126 targets).

* **AVG** — production's terminal. `core/pipeline.py:951`: superpose the 75 on their
  unweighted consensus medoid, take the uniform coordinate mean → `C` (a point cloud, not a
  structure) → stage-3b projection `I.project(C, seq, fold)` → built chain.
* **MED** — the quantum stage's readout operator (`core/pipeline.py:795` `consensus_medoid`,
  `w=None`): return the member of the 75 that minimises its mean pairwise RMSD to the other 74.
  This is **one real deposited backbone**. It is passed through the *same* stage-3b projection so
  the two arms differ in exactly one place.

`MED`'s cloud RMSD is a pool member's own `rr` entry — no new geometry is created, so `MED` cannot
violate physical validity by construction. That is the asymmetry the hypothesis rests on.

---

## 2. The theorem the coordinator asked me to establish, DERIVED HERE BEFORE MEASURING

For an **averaging** terminal the project already has `set_mean² = B² + S²` with `B` the endpoint,
so pure concentration is worth zero by algebra. The medoid is not an average, so that identity does
not describe it. **The analogous statement, derived:**

Superpose the 75 members `x_1..x_75` into a common frame; let `μ` be their centroid (= the
average operator's output cloud), `t` the native, `B = |μ − t|` (this is `AVG`'s cloud error), and
`S² = mean_i |x_i − μ|²` (the set's spread). Let `b` be the medoid index and
`s_b = |x_b − μ|`. Then, with the residual direction `x_b − μ` approximately orthogonal to
`μ − t` (high-dimensional, 3n ≈ 39 coordinates):

```
d_MED²  ≈  B² + s_b²   ≥   B²  =  d_AVG²          (on the CLOUD)
```

and since the medoid is the member nearest the set's own 1-median, `s_b > 0` strictly whenever the
set is not degenerate. **So on the point cloud the average dominates the medoid by construction,
and the margin grows with the spread `S`.** Concentration is *not* worth zero for the medoid — but
it works in the direction that makes the medoid *worse*, not better, on every target.

**Therefore the medoid's only possible route to a win is the projection.** Write
`P(op) = chain(op) − cloud(op)` for the stage-3b penalty. `AVG` emits a contracted non-structure
and pays `P(AVG) = 3.2105 − 3.0483 = +0.1622 Å` on average; `MED` emits a real backbone and should
pay `P(MED) ≈ 0`. The whole hypothesis is the single inequality

```
    [ d_MED_cloud − d_AVG_cloud ]     <     [ P(AVG) − P(MED) ]
      the medoid's structural cost          the projection penalty it avoids
```

with `P(AVG) ≈ 0.1622` the budget. **This localises the lane to one measurable number and one
budget, and it does not close it** — but it is a *prior-lowering* derivation, because the left-hand
side is a spread-scaled quantity and the coordinator's gate (high dispersion → prefer MED) pushes
on the variable that *inflates* the left-hand side as well as the right.

**Registered consequence, which I will check first as a cheap pre-check (contract rule 22):**
if `d_MED_cloud − d_AVG_cloud` is already larger than 0.1622 Å *on the whole sample*, the whole-
sample F1 is dead on arrival and the lane reduces to the gated and tail questions.

---

## 3. F1 — the gate, and that it is native-free

**Gate variable `DISP`** = mean of the off-diagonal of the top-75 pairwise-RMSD block
(`block.mean()`), i.e. **the very criterion the medoid minimises**. It uses only the retained
candidates. It requires no native, no label and no held-out fit. `s30_G_disp.json` records it under
`DISP_rmsd`; lane F's S30 `F4_top75_rg_sd` is an `rg`-based sibling at ρ = 1.0000 with lane G's
`DISP_rg` **by construction** (one reads the other's column) — a membership check, not a
reproduction.

**Registered gating rules, in order of strength. I commit to reporting all three.**

| id | rule | free parameters | deployable? |
|---|---|---|---|
| **G0** | `DISP > median(DISP over the 126)` → MED, else AVG | 0 | deployable (a rank rule, no native) |
| **G1** | `DISP > τ`, τ fitted leave-fold-out on native chain RMSD | 1, held out | **ORACLE-ADJACENT**: τ touches the native. Reported as a diagnostic, labelled in the same sentence as its number. |
| **G2** | per-target `min(AVG, MED)` | — | **ORACLE / NOT DEPLOYABLE — best-of-2, not skill.** Reported only with a split-half transfer arm (fit the choice on half the targets, apply to the other half, both directions). |

Registered *direction*: G0 sends the **high**-dispersion half to MED. This direction is fixed now by
the mechanism (averaging damage concentrates in the divergent half; S30 `−0.0406 Å` hi-vs-lo
contrast, 3.56× MDE, 97.2% of the k=30 repair gain in the high-dispersion half). **If the measured
effect is the other way round, G0 has fired against me and I report it as a refutation, not as a
sign flip.**

### Falsifier, in Å on the built chain

Let `Δ = mean_chain(G0) − mean_chain(AVG)` over the 126, paired, fold-clustered CI on the pinned
folds, MDE = 2.8016 × SE **per comparison**.

* **F1 is CONFIRMED** only if `Δ ≤ −1.0 × MDE(Δ)` and the fold CI excludes zero.
* `−1.0 × MDE < Δ < −0.7 × MDE`: **NOT MEASURED**.
* `Δ > −0.7 × MDE`: **F1 REFUTED**; I will say so and quote the number.

**Mechanism falsifier, registered separately** (a gain with the wrong mechanism is not this lane's
claim): the hi-minus-lo `DISP`-half contrast of `chain(MED) − chain(AVG)` must be **negative**. If
the whole-sample `Δ` is favourable but the contrast is not negative, the gate is working through
something other than dispersion and I report the gain as *unexplained*, not as confirmation.

### Registered secondary: geometry vs accuracy, split

The coordinator's stated worry is that MED helps geometry and hurts accuracy. I register that I
will report, for both operators and both DISP halves, **separately**: (a) chain Cα RMSD (accuracy);
(b) the projection penalty `P` and the radius-of-gyration contraction `Rg(out)/Rg(native)`
(geometry). **If MED wins on (b) and loses on (a), both halves are quoted, and the accuracy half
governs the verdict.**

---

## 4. F2 — what distinguishes the tail (charter §13)

**Primary tail definitions are filter-independent** (contract rule 12: a stratum defined by the
outcome cannot measure the thing that defines it):

* **T-POOL** — worst 18 by `pool_mean` (a property of retrieval, not of the filter).
* **T-BEST** — worst 18 by `pool_best` (ORACLE, but *not* a function of the filter's recall).
* **T-CHAIN** — worst 18 by production built-chain RMSD (the outcome; overlap with FAIL18 is 13/18).

`FAIL18` appears **only** as a cross-check and is labelled *defined by the filter's own recall
(`s12/instrument.py:271-278`) and therefore unable to measure filter recall* every time it appears.

Registered F2 questions and the statistic each is answered with:

1. *Are useful candidates present but badly ranked?* → `pool_best`, `top_m_best`, and the **rank of
   the pool's best member under the shipped score**, tail vs rest.
2. *Is the terminal operator failing differently?* → `chain(MED) − chain(AVG)`, `P(AVG)`, and
   contraction, tail vs rest.
3. *Do candidate clusters behave differently?* → `DISP`, the number of distinct structures
   (`n_distinct`), and the top-75's spread `S` relative to `B`.
4. *Is the native-improving subspace different?* → the fraction of the 75 that beat `AVG`'s own
   cloud error, tail vs rest.

No F2 number tunes anything. F2 is diagnostic and is labelled ORACLE wherever it uses `rr`.

---

## 5. Honest prior odds

* **F1 (G0 beats production at ≥1.0× MDE on the built chain): ~1:4 against.** The coordinator
  registered ~1:1; I am registering *lower* than my coordinator, and the reason is §2 above, which I
  wrote before seeing a number: the medoid's cloud cost and the projection penalty it avoids are
  **both** increasing functions of the same spread, so the gate's premise does not obviously
  separate them, and the projection budget is only 0.1622 Å.
* **The mechanism contrast (hi-minus-lo negative): ~1:1.** This one I genuinely do not know, and it
  is the part worth the compute even if F1 fails.
* **F2 producing a durable statement about the tail: ~3:1 in favour**, because the tail statistics
  are cheap, filter-independent and have never been put side by side.
* **That the lane closes by derivation:** ~1:3 against — §2 constrains but does not close.

## 6. Multiplicity

Registered comparison count for lane F before any number: **3 operator arms × (whole sample + 2
DISP halves) + 3 gate rules + 4 tail definitions × 4 F2 statistics = 9 + 3 + 16 = 28 comparisons.**
Any comparison beyond these is declared as unregistered when it is reported.

## 7. Artefacts

`s31/s31_F_terminal.py` → `s31/results/s31_F_terminal_rows.jsonl`, `s31/results/s31_F_terminal.json`.
Seed 31006. Config: top-75, POOL_K=500, MIN_SEP=2, `I.project(lam=0.3, multi=True, maxiter=300)`.

---

## 8. ADDENDUM — a third operator arm, registered before the main run produces any number

**Appended 2026-09-20 (see the commit timestamp). No number from the main run exists at this
point; the only numbers seen so far are the three-target scout (`s31/s31_F_scout.py`), whose role
was the reproduction gate and the projection timing and whose values are reported in the ledger as
a scout, not as a result.**

§2's derivation says `AVG` pays `P(AVG) ≈ +0.1622 Å` at stage 3b because it emits a *contracted
non-structure* (the scout saw a CA–CA bond mean of **1.97 Å** against a real backbone's 3.79 Å on
one target). The project's standing measurement is that coordinate averaging **contracts the
backbone by 25.8%**. That suggests a third terminal operator which is neither the average nor the
medoid:

* **AVG_RG** — the uniform coordinate average, rescaled about its own centroid so that its radius
  of gyration equals the **mean Rg of the 75 members** that produced it, then projected through the
  same stage 3b.

The rescaling target is a property of the retained candidates only: **native-free, zero free
parameters, deployable.** Falsifier identical in form to F1: `mean_chain(AVG_RG) − mean_chain(AVG)`
must reach `−1.0 × MDE` with a fold CI excluding zero to count; `> −0.7 × MDE` refutes it.

Registered prior: **~1:2 against**, because stage-3b projection restores ideal bond geometry and is
therefore *already* a partial de-contraction, so the correction may be double-counted.

This raises lane F's registered comparison count from 28 to **32**.
