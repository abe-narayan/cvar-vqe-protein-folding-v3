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

> **CORRECTION, appended 2026-09-21 00:05 — the 25.8% in the paragraph above is WITHDRAWN and the
> original wording is left standing (contract rule 13).** The coordinator supplied it from project
> memory `averaging-space-beats-the-objective`; it is carried stale there and at
> `s30/LEDGER.md:129`, and it was withdrawn at `s15/coord_FINDINGS.md:914-921`. **This was the
> coordinator's error, not this lane's, and it is recorded that way at the coordinator's own
> request.** The corrected measurement is *contraction against true distances = 3.5%*, and — far
> more importantly — the distortion is **not a scale effect at all**. See §10.
>
> `AVG_RG` is therefore a correction of the **wrong shape**, and it is retained in the run only as
> a **registered negative control**: it is the uniform-scale fix that the corrected mechanism says
> should *not* work. Its ~1:2-against prior is now closer to ~1:5 against. Retaining it is
> deliberate — a pre-registered arm is not deleted because its motivation was withdrawn, and an
> arm that fails for the newly-understood reason is evidence for the new mechanism.

* **AVG_RG** — the uniform coordinate average, rescaled about its own centroid so that its radius
  of gyration equals the **mean Rg of the 75 members** that produced it, then projected through the
  same stage 3b.

The rescaling target is a property of the retained candidates only: **native-free, zero free
parameters, deployable.** Falsifier identical in form to F1: `mean_chain(AVG_RG) − mean_chain(AVG)`
must reach `−1.0 × MDE` with a fold CI excluding zero to count; `> −0.7 × MDE` refutes it.

Registered prior: **~1:2 against**, because stage-3b projection restores ideal bond geometry and is
therefore *already* a partial de-contraction, so the correction may be double-counted.

This raises lane F's registered comparison count from 28 to **32**.

---

## 9. SECOND ADDENDUM — a native-free gate variable I did not think of first, and exactly what I
had seen when I added it

**Appended 2026-09-20 23:58. FULL DISCLOSURE, because this is registered mid-run and the honest
thing is to say what was on the screen.** At the moment of writing, the main run had emitted
**8 of 126 rows**, and the two things visible in them were: (i) `P(MED)` is ≈ 0 to ±0.003 Å on
every row, exactly as §2 derived; (ii) `P(AVG)` is *not* one-signed — on `1D6X` it is **−0.309 Å**,
i.e. stage-3b projection sometimes *improves* the average. **No aggregate, no gate, no contrast and
no F1 number existed.**

Fact (ii) is what prompted this. `P(AVG)` is ORACLE (it needs the native), but the closely related

```
MOVE(op)  =  CA-RMSD( input cloud ,  its stage-3b projection )
```

is **fully native-free** — both arguments are available at inference. It measures *how far the
projection has to drag the operator's output to make it a chain at all*, which is the mechanism
§2 identified, measured directly instead of through the dispersion proxy `DISP`.

* **G3** — `MOVE(AVG) > median(MOVE(AVG))` → MED, else AVG. 0 free parameters, native-free.
* **G3τ** — leave-fold-out threshold on `MOVE(AVG)`: **ORACLE-ADJACENT**, diagnostic only.

**Falsifier, identical in form to G0:** `mean_chain(G3) − mean_chain(AVG)` must reach `−1.0 × MDE`
with a fold CI excluding zero. `> −0.7 × MDE` refutes it. **G3 is declared EXPLORATORY** — it was
registered after the run began — and if it is the only arm that clears, it is reported as
*exploratory, requiring confirmation*, never as the lane's confirmed result, and it carries a
split-half transfer arm.

Registered prior for G3: **~1:2 against**, better than G0's ~1:4 because it measures the mechanism
instead of a proxy for it, but still against, because §2's inequality has the medoid's cost rising
with the same quantity.

Lane F's registered comparison count rises from 32 to **36**. The run was restarted from zero so
that every row carries the `MOVE` columns; the 8 rows seen above were discarded and recomputed.

---

## 10. THIRD AMENDMENT — the withdrawn 25.8%, the corrected mechanism, and a defect of my own

**Appended 2026-09-21 00:10, on the coordinator's course correction, with the run stopped at 8 rows
and restarted from zero. No F1 aggregate exists at this point.**

### 10.1 The withdrawal (the coordinator's error, recorded as such at their request)

The brief that opened this lane asserted *"averaging contracts the backbone by 25.8%."* **That
number is withdrawn** at `s15/coord_FINDINGS.md:914-921`, in its author's own words: *"the magnitude
was overstated by a factor of seven."* I verified the withdrawal at that path rather than accepting
the relay (contract rule 14). It is still carried stale at `s30/LEDGER.md:129` and in project memory
`averaging-space-beats-the-objective`. The corrected figures:

```
contraction against TRUE distances          3.5%   (not 25.8%)
contraction against PREDICTED distances     9.5%
separation profile, monotone                0.773 at the virtual bond -> 1.10 at |i-j| = 13,
                                            crossing 1.00 near |i-j| = 8
```

Every sentence of this prereg that rested on the 25.8% is annotated in place above with the
original wording left standing (contract rule 13). §2's derivation does **not** rest on it — it
rests on `P(AVG) = chain − cloud = +0.1622 Å`, which is measured on the endpoint basis and is
unaffected.

### 10.2 The corrected mechanism, and it is sharper than the one it replaces

The distortion is **not a scale effect**. It is a **separation-dependent shape distortion**:
averaging *contracts* short-range CA–CA distances and *expands* long-range ones, crossing unity
near |i−j| = 8. Beside two results already on disk:

* the tail's error is **83.1% shape**, with scale **refuted** as the mechanism (partial
  ρ = −0.067, p = 0.46, against shape's −0.641);
* **68% of the recoverable prize is in |i−j| ≥ 7**, concentrating ~4× on both filter-independent
  tails.

**The averaging operator's distortion crosses unity almost exactly where the prize lives.**

### 10.3 What I am adding, and its falsifier

1. **Instrumentation (ORACLE diagnostic, no parameter):** every arm's **separation profile** —
   mean CA–CA distance at each s = 1..n−1 — for the native, the 75 members, and every operator
   output on **both** the cloud and the chain basis. Registered questions: *is the medoid's profile
   flatter than the average's in the s ≥ 7 band?* and *does `chain(MED) − chain(AVG)` correlate
   with the average's s ≥ 7 distortion?*
   **Mechanism falsifier, band version:** the medoid's advantage must be larger in the s ≥ 7 band's
   distortion than in the s < 7 band's. If the medoid's profile is *not* flatter at long range, the
   separation-band framing is refuted for this operator pair and I say so.
2. **A fourth operator, AVG_SEP (native-free, 0 free parameters):** rescale the average's pair
   distances so that its per-separation mean matches **the 75 members' own per-separation mean**,
   re-embed by classical MDS, project through the same stage 3b. This is the correction *matched to
   the corrected mechanism*, where `AVG_RG` is the correction matched to the *withdrawn* one.
   **Falsifier:** `mean_chain(AVG_SEP) − mean_chain(AVG) ≤ −1.0 × MDE` with a fold CI excluding
   zero; `> −0.7 × MDE` refutes it. Registered prior **~1:3 against**, because §10.4's record shows
   every *scale* correction of the average is already closed on the endpoint basis and a
   per-separation correction is a richer member of the same family.
3. **`AVG_RG` is demoted to a registered negative control** and to a reproduction check against the
   closed `pool`-scale arm (+0.095 Å, CI excluding zero on the bad side).

### 10.4 A DEFECT OF MY OWN, declared before it can be found for me

**I wrote §§1–3 of this prereg without searching the record for a prior medoid-as-terminal
measurement, and there is one.** `s12/agg_FINDINGS.md:43-70`:

```
medoid75      3.2822 (CA POINT CLOUD)   +0.2339 vs avg75   CI [+0.162, +0.305]   34W/92L
              FAIL18 6.100 vs avg75's 5.832   other108 2.813 vs 2.584
```

That is contract rule 14 in the other direction — I should have run the `grep` before writing the
registration, and lane C had already found it. **Two things follow and I record both:**

* **My §2 derivation predicted this before I saw it** — it says the medoid must lose on the cloud,
  by a margin that grows with the spread, and the record's +0.2339 is *larger* than the entire
  projection budget of +0.1622. **That is my own registered pre-check firing against my lane's
  hypothesis, and it is the best thing in this registration.** The whole-sample F1 arm is therefore
  expected to be **refuted**, with an implied chain-basis value near +0.2339 − 0.1622 ≈ **+0.07 Å
  worse**, and I am writing that predicted number down now so it can be checked against the
  measurement.
* **What is NOT on disk, and is what this lane still measures:** the medoid on the **built-chain**
  basis (the record is cloud-only and the projection is exactly the medoid's advantage); the
  **dispersion- and MOVE-gated** per-target choice; the **separation-band mechanism**; and the
  tail behaviour on the two **filter-independent** tails (the record's tail column is FAIL18, which
  cannot measure the filter that defines it).

I am **not** re-deriving the convex-reweighting rung: lane C reports the convex optimum under the
deployed objective converges to **3.0522 CA cloud**, i.e. to the uniform average production already
emits, on two independent constructions. The medoid/average *choice* is a different object.

### 10.5 The three readouts, named so my arms cannot be confused

| code | what it emits | this lane's arm |
|---|---|---|
| `core/pipeline.py:951` `average_struct`, called at `:1116` | uniform coordinate mean of the top-75 → a **new** point cloud | **AVG** (the deployed terminal) |
| `core/pipeline.py:795` `consensus_medoid` | an **index** — one real deposited pool member | **MED** |
| `core/pipeline.py:880-895` `average_weighted`, called at `:1110` | p_θ-weighted convex mean → a **new** point cloud | **not run here** — it needs the VQE's `p`, and its convex ceiling is closed at 3.0522 (lane C) |

Lane F's registered comparison count rises from 36 to **44**.

---

## 11. FOURTH AMENDMENT — F3: IS `bestm128 = 2.9027` AN ORDER STATISTIC?

**Appended 2026-09-21 00:33 on the coordinator's new highest-priority item, before any aggregate
from the `bestm128` family has been computed.** The only row I have opened is `1A13`
(`s29/results/s29_O_p128_rows.jsonl`), to learn the schema.

### 11.1 What is being audited, and the structural fact that decides how

`s29/LEDGER.md:3943` and `s29/results/s29_O_p128_rows.jsonl` hold, per target, the full
`curve128` — the CA-point-cloud RMSD of the prefix average of the first `m` of the score-ordered
top-128, for **every** `m = 1..128` — and `bestm128 = min_m curve128[m]`. The **built-chain**
value 2.9027 comes from projecting *one* structure per target: the prefix average at the
**cloud-argmin** `m_best` (`s29/s29_O_ladder.py:542`). **So the chain-basis 2.9027 is a
cloud-selected oracle, projected** — not a chain-selected one. That is stated in every sentence
this lane writes about it.

`bestm128` is a **per-target minimum over K = 128 variants**, which is the exact shape contract
rule 11 exists for.

### 11.2 The registered arms

All on the **CA point cloud** where the full 126×128 matrix exists (exact, K = 128, no
projections), with the **built chain** carried for the arms that survive, projected **in the same
job from the same stored clouds** (lane D's 1e13 branch amplification).

| id | arm | native? |
|---|---|---|
| **F3-a** | `min_m curve128[m]` vs `curve128[75]` — the ORACLE gain being audited | **ORACLE / NOT DEPLOYABLE** |
| **F3-b** | random `m` per target — the **zero-skill** control | native-free |
| **F3-c** | ORACLE **global** `m` (one m for all 126) and its **leave-fold-out** twin — the transferable part of the axis | LFO = held out |
| **F3-d** | `stats_lib.split_half_transfer` over the 126×128 matrix | — |
| **F3-e** | **the matched random-variant-family control** (contract rule 7): per-target min over K = 128 *random* subsets of the same top-128, sizes drawn to match the prefix-size distribution. Same operator, same set, same K, **index carries no score order** | ORACLE min over a null family |
| **F3-f** | native-free per-target `m` rules (n, DISP, MOVE, pool-score spread, `rg` sd, `n_distinct`), **fold-held-out** | native-free, deployable if it clears |
| **F3-g** | effective K: curve smoothness / count of distinct local minima, and the best-of-K inflation re-priced at K_eff | — |

### 11.3 REGISTERED BARS — written before the numbers, so they can fire

1. **"The m axis is an order statistic"** fires if **F3-e reaches within 20% of F3-a's gain**, i.e.
   if a min over 128 *arbitrary* subsets of the same 128 buys ≥ 0.80 × 0.3084 = **0.247 Å**. In
   that case 2.9027 must be reported as *"7 bits of ORACLE index over the top-128"* and never as
   *"the value of choosing the prefix length"*.
2. **"No part of it is deployable"** fires if the best **F3-f** rule's fold-held-out gain is
   `> −0.7 × MDE` of its own comparison. Then the number may never again be quoted as a reachable
   ceiling without that sentence attached.
3. **F3 is a POSITIVE** only if some F3-f rule reaches `≤ −1.0 × MDE` with a fold CI excluding
   zero **on the built chain**.

### 11.4 My registered prior

* **That a substantial fraction of the 0.3084 Å is best-of-128: ~4:1 in favour** (the coordinator
  registered ~2:1). My reason is the arithmetic already in the record: `best1_top128` spends the
  **same 7 bits** and is worth **−1.0657**, i.e. 3.5×. If 7 bits of *arbitrary* oracle index over
  this set is worth ≥ 1.07, then 0.31 for one particular 7-bit index is not evidence that *that
  index* is special — it is a lower bound on what any 7-bit oracle buys.
* **That a native-free per-target m-rule clears 1.0× MDE on the built chain: ~8:1 against.**
  S29 already measured the ORACLE **global** m at **−0.0018** and its LFO twin at **+0.0079
  worse**, and S30 closed filter width globally. The per-target question is genuinely open, but
  every per-target native-free rule this project has built has landed in the hundredths.
* **I do not expect F3 to be the sprint's best shot, and I am saying so before I measure it.** The
  coordinator rates it highest; I rate `AVG_SEP` higher, and both of us should be on record.

### 11.5 Multiplicity

F3 adds 7 registered comparison families. Lane F's registered count rises from 44 to **51**, and
the family is appended to `s31/MULTIPLICITY.md` as it is emitted.

---

## 12. FIFTH AMENDMENT — two gates for `AVG_SEP`, registered before any `AVG_SEP` aggregate exists

**Appended 2026-09-21 00:47. The main run stands at 93 of 126 rows. I have NOT computed any mean,
median, W/L or contrast for `AVG_SEP` and will not until 126/126 — that commitment is on the record
in two messages to the coordinator.** What I have seen is individual log lines, and they are mixed:
`1G89` −0.409 and `1I6Y` −0.578 in `AVG_SEP`'s favour, `2MID` **+0.884** against it. A mixed
per-target distribution is exactly the shape a gate exists for, so I register the gates now rather
than after the aggregate tells me whether I need one.

* **G4** — `MOVE(AVG) > median(MOVE(AVG))` → `AVG_SEP`, else `AVG`. Native-free, 0 free parameters.
* **G5** — `DISP > median(DISP)` → `AVG_SEP`, else `AVG`. Native-free, 0 free parameters.
* **G6** — per-target `min(AVG, AVG_SEP)`: **ORACLE / NOT DEPLOYABLE, best-of-2 not skill**,
  reported only with a split-half transfer arm.

Both G4 and G5 are **EXPLORATORY** — registered after the run began — and each carries a
split-half transfer arm. Falsifier identical in form to F1: `−1.0 × MDE` with a fold CI excluding
zero to count, `> −0.7 × MDE` refutes.

**Registered directions, fixed now.** G4 sends **high**-MOVE targets to `AVG_SEP`, and G5 sends
**high**-DISP targets to `AVG_SEP`, because both variables index *how badly the average is
distorted* and `AVG_SEP` exists to undo that distortion. **If the measured effect runs the other
way, these fire against me and are reported as refutations, not as sign flips.**

**One thing I am registering so it cannot be claimed afterwards.** `s31/results/s31_F_coh.json`
already shows `coh(AVG_SEP) = 0.9689` against `AVG`'s 0.9780 and a bar of 0.6931, with **0 of 126**
targets under the bar and an affine-hull residual of 0.045 Å RMS per coordinate. So **whatever
`AVG_SEP` does at the endpoint, it does NOT do it by leaving the pool's affine hull.** If it wins, I
will say the mechanism is unexplained rather than attach lane B's theorem to it.

Lane F's registered comparison count rises from 51 to **57**.
