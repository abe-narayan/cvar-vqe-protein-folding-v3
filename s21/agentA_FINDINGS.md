# SPRINT 21 — WORKSTREAM A FINDINGS
## The mandatory Hamiltonian matrix, priced before it was built

Pre-registration: `s21/PREREG_A.md` (falsifier first, unedited) plus the coordinator's
pre-registration in the `s21/tailprice.py` docstring and the coordinator-issued amendment logged
beneath it. Artefacts: `s21/results/tailprice.json`, `s21/results/a_matrix.json`.

**Basis is stated on every structural row, and so is the READOUT** (§1.5b: a Hamiltonian's sign is
a function of its readout).

**MDE IS PER COMPARISON.** Every paired contrast here carries its own standard error and its own
`2.8016 × SE` minimum detectable effect. The project's pooled **0.084 Å constant is not used as a
threshold** — it is wrong by up to 84× in both directions on an individual contrast, and a "null"
judged against the wrong MDE is not a null.

---

## 0. THE ONE-LINE RESULT, and it damages my own lane

**The mandatory matrix was bounded before a single circuit was built, and then the circuits agreed
with the bound.** Over the shipped K=500 pool, the α-tail of **AMBER is significantly WORSE than a
matched-count random subset of the same size** — on every readout, at every tail width up to
α=0.30, on both geometric bases, with 0/5 folds on the good side. **No Hamiltonian containing
Legacy or AMBER beats the structural-only arm.** The pre-registered falsifier F-A1 fires: a
CVaR-VQE selecting over this candidate set is bounded at chance on those two Hamiltonians *no
matter how well it optimises*.

**The mandatory matrix (§2), run with genuine CVaR-VQE on 12 targets, 4 seeds, matched budget,
delivered exactly that.** Distance-only **3.845 / 3.686** (mean/median) wins; every one of the six
Hamiltonians containing Legacy or AMBER is worse, all six point estimates on the harmful side,
AMBER worst at **+1.471** vs Distance with **2W/10L**. **AMBER's VQE does not beat a constant
α-helix with noise** (+0.131 [−0.681, +1.035], 6W/6L).

> **Two scope limits on that paragraph, both found by the coordinator's fork review and both
> narrowing my own claim.** (i) `H_AMBER` is the **bare single point**, not the deployed
> `E ∘ Relax₅₀` whose own cost Sprint 20 puts at +0.021 — **+1.471 Å is not the pipeline's AMBER
> cost**, and only **two** of the six rows (Legacy, Dist+Leg) are AMBER-free, so the direction is
> 2/2 independent of that fork and 4/4 conditional on it (§2.8, OPEN O7). (ii) **Seed count was a
> fork I failed to enumerate** and it is the one that most constrains the table (§2.3).

**And the honest counterweight: at this budget the matrix cannot resolve its own magnitudes.**
Ansatz-seed sd is **0.51–1.06 Å**, two to five times Sprint 20's figure, and the VQE does **not**
beat best-of-N from the untrained circuit on any Hamiltonian (F-A2 **NOT MEASURED**, a scope
restriction on Sprint 20 Q4, not a refutation). Every §2.2 magnitude is **NOT MEASURED**; only the
**direction** — six of six against physics — is established.

**Three qualifications, each of which cost me a claim I would otherwise have made.**

1. **A Hamiltonian's sign is a function of the readout, not of the Hamiltonian.** Legacy *beats*
   its matched control at the `member` readout (−0.41 [−0.56, −0.26], n=126) and *loses* to it at
   `medoid` (+0.30) and `avg` (+0.33). Both are true. **"AMBER and Legacy are worse than chance"
   is only a complete sentence with a readout attached**, and the mechanism is error coherence:
   averaging preserves Legacy's systematic compactness bias where a single random member escapes
   it. AMBER is the one that fails on *every* readout.
2. **The tail operator does help AMBER, a lot, and it does not matter.** With the averaging
   confound removed, `medoid(tail) − argmin` on AMBER is **−1.862 [−2.555, −1.147]** at α=0.30, 5/5 folds.
   The coordinator's original hypothesis was right. The level reached is still worse than random.
3. **My headline number for the distogram is the shipped pipeline, not a discovery.** At α=0.15
   the tail-average selector under the distogram *is* the production operator — it reproduces
   `rmsd_avg` to `max |diff| = 0.000000` on 126/126. That is a soundness gate passed.

So the informative content of the mandatory matrix is **how much damage each physics Hamiltonian
does to a good selector, per readout** — not a search for a winner.

---

## 1. THE TAIL-PRICING BOUND — `s21/tailprice.py`

### 1.1 What was fixed before it ran, and why the fixes are not cosmetic

The file was handed over with two wrong API calls. Both are recorded verbatim in the module
docstring so the diff is auditable; neither changed an arm, a null or an endpoint.

1. `legacy_components_of_windows(pdb, W, seq, fold)` → **`(seq, PHI, PSI)`**. The function takes
   torsions, not coordinates. The stale artefact in `s21/results/tailprice.json` proves the call
   never worked: three rows, all `skipped`, all
   `ValueError("could not convert string to float: 'INWKGIAAMAKKLL'")` — the sequence being fed
   to a float cast. It is retired as `_SUPERSEDED_tailprice_coordinator_smoke.json`.
2. `s13.qarch_lib.amber_energies(space, S)` is the **lattice** path (a discrete `Space` and
   integer states) and cannot score an arbitrary retrieved window. Replaced with the genuine
   ff14SB/GBn2 **single point** on the coordinate path, exactly as `s18/phys_down.py` does it and
   as Sprint 20's gate GC20a certified:
   `ConstrainedBox(seq, rep).energy_point(build_backbone(phi, psi))`. **No minimisation:
   H_AMBER = E, not E∘Relax₅₀.**

**A third defect was mine to catch and it would have cost the run.** The resume path loaded the
stale artefact's three `skipped` rows, so the restarted job would have **silently dropped three
targets from n=126** while reporting a full run. `run()` now discards any row failing
`row_complete` before resuming, and says how many it discarded.

**A fourth is the completion flag**, corrected on the coordinator's instruction: `complete` now
requires n_expected rows **and zero skips and all 946 Hamiltonian × α × readout cells per row**.
The old flag counted skipped rows and had certified a run of total failures as complete.

### 1.2 Basis — verified, not assumed

Legacy and AMBER exist only on the **ideal-geometry rebuild** of a member's torsions; the shipped
distogram filter scores the **real retrieved window**. I checked rather than assumed: the shipped
top-75 index set is reproduced **75/75** by the distogram score on `W`, and only **54/75** on the
rebuild. So `disto` scores `W`, and Legacy/AMBER score the rebuild.

Because that asymmetry could by itself manufacture the result, **every arm is emitted and scored
twice** — once as the coordinate average / medoid / member of the selected **real windows** (`W`),
once of the selected **rebuilds** (`RB`).

> **The basis is not the explanation.** Per-member rebuild − window shift is **+0.010 Å**,
> Spearman(d_W, d_RB) = **0.988**, and — the check that actually counts — **the RB table
> recomputes the medoid, the superposition frame and the average entirely in the rebuild basis**
> and reproduces the W table to within ~0.06 Å on every cell. AMBER's failure survives being
> scored entirely inside its own geometry.

> **Two different quantities, and only one of them is the operator's price.** My +0.010 Å is a
> **per-member** shift; the coordinator's **operator-level** measurement is
> `pool_avg75(rebuilt) − ship_avg75(window) = +0.016 [−0.004, +0.037]` at n=126 — negligible on a
> 75-member coordinate average. **The per-member number could not have been borrowed for the
> operator**, because the rebuild can change *which member is the medoid*, and the medoid sets the
> superposition frame. That is precisely why every arm here is emitted on both bases rather than
> corrected by a scalar.

### 1.3 Normalisation, declared before RMSD (BRIEF §6)

Measured first, per target, over each target's own K=500 pool (medians of the per-target
statistic, **n=42**):

| component | mean | sd | skew | q50 | q99 | q100 | frac > 1e4 | finite |
|---|---|---|---|---|---|---|---|---|
| disto | 2.28 | 0.91 | +1.55 | 2.05 | 5.29 | 6.14 | 0.00 | 1.00 |
| legacy | −10.9 | 7.6 | +0.64 | −12.0 | +6.8 | +16.2 | 0.00 | 1.00 |
| **amber** | **1.0e14** | **1.5e15** | **+22.1** | **3.6e4** | **3.2e12** | **3.1e16** | **0.582** | **1.00** |

AMBER's per-target range is ~16 decades and a median **58.2%** of pool members sit above
1e4 kcal/mol — Sprint 20's L5c (53.5% on the top-75) reproduced and *exceeded* on the full K=500
pool, which is the right direction since the top-75 is already distogram-filtered. It is
**finite everywhere** (finite fraction 1.000): *finite-but-meaningless, nothing throws.*

**And the reason neither can select — but only after a correction I had to be given.** ORACLE
rank skill against true RMSD, per target, n=42, **fold-clustered CI**:

| | Spearman(E, d_true) MARGINAL | fold CI | right sign |
|---|---|---|---|
| disto | **+0.618** | [+0.523, +0.707] | 39/42 |
| legacy | **+0.328** | [+0.301, +0.355] | 34/42 |
| **amber** | **−0.042** | **[−0.116, +0.032]** | **19/42** |

> **I first reported that middle row as "Legacy carries about half the distogram's in-band rank
> skill". That reading is wrong, and my own artefact contained the refutation three lines above
> it.** `rho(disto, legacy)` is **+0.478** (median +0.624) — Legacy is strongly rank-correlated
> with *the score already in production*. **A marginal correlation with the truth, for a score
> correlated with the deployed selector, is not that score's contribution.** Partialling the
> distogram out costs nothing; all three coefficients were already recorded per target:

| | PARTIAL given disto | fold CI | negative |
|---|---|---|---|
| legacy \| disto | **−0.0076** | [−0.0616, +0.0681] | 25/42 |
| amber \| disto (control) | −0.0290 | [−0.0593, −0.0036] | 25/42 |

**Legacy retains −2% of its marginal rank skill: it carries NO in-band rank information beyond
the distogram.** The AMBER row is the control that shows the partialling is not manufacturing
structure — its marginal was already ~0 and its partial stays ~0. **And my own tables already
contained the consequence: `d+l` never beats `d`, at any α or any readout. This is the mechanism
for that.** Same class of error as Sprint 20's Q11, where every circuit-side landscape metric
collapsed once target difficulty was partialled out.

**The mandatory matrix's Legacy row is therefore reported as a PARTIAL, not a marginal.**

**AMBER carries no rank information about accuracy at all** — and its *tail* is nonetheless
significantly worse than random, so the relationship is not merely absent but **non-monotone**:
the members AMBER scores best are systematically the wrong ones. Sprint 20 L2c already named the
mechanism — AMBER-preferred candidates are *expanded* with open sterics while Legacy-preferred are
0.45 Å more compact — and an expanded chain is further from a compact 9–16-mer native. That
mechanism is **cited, not re-measured here.**

A raw sum of these three is an AMBER-outlier detector wearing a hybrid's name. The declared form
is therefore the **within-target rank-to-normal** transform — scale-free, invariant to any
monotone reparameterisation of either energy, outlier-robust by construction, no fitted constant.
Its cost is declared too: it discards magnitude, deliberately, because Sprint 20 L5c shows the
magnitude channel is a steric singularity and not an accuracy signal.

**Audit** (raw sum and per-target z-score), computed alongside and reported whatever they say —
see §1.7, which records the one cell where a non-declared form wins.

**A legend in my own report was inverted and is fixed.** `tailprice.py` annotated the rank-skill
block *"(negative = lower energy is better)"*. The sign is backwards — `disto`, the score that
works, is **+0.618 with 39/42 targets positive** — so as printed the legend said the arm that
works is the arm that fails. Corrected to *"POSITIVE = lower energy goes with lower RMSD =
CORRECT ordering"*.

### 1.4 THE PRIMARY (amended) — arm minus its matched-count random control

The coordinator re-primaried this before the table was read, and the reason is worth keeping:
the original primary `tail(0.15) − argmin` compares an **average of 75 structures** to **one
structure**, so it differs by the averaging operator (~1.0 Å, Sprint 19) and would have
"confirmed" the >0.3 Å prediction on *every* Hamiltonian including a random one. The endpoint is
now `tail{α} − random{α}` at matched count, with each readout compared to **the same readout of a
random subset of the same size**.

**n = 42 of 126.** The run was paused to hand the shared AMBER context to §2's mandatory matrix;
the artefact's `complete` flag correctly reads **false**, `n_expected` stays 126, and the achieved
n is printed on every table. Mean Å, negative = the Hamiltonian's tail carries information over
chance. CIs are paired target-level bootstrap, **i.i.d. and fold-clustered both quoted**, with
**SE and the comparison's own `2.8016 × SE` MDE**, W/L and per-fold sign on every row.

**α = 0.15, all three readouts side by side — because the readout decides the sign:**

| readout | disto | legacy | **amber** | d+l | d+a | l+a | d+l+a |
|---|---|---|---|---|---|---|---|
| **member** (1 structure) | **−0.821** ✓ | **−0.348** ✓ | **+0.605** ✗ | **−0.814** ✓ | **−0.690** ✓ | **−0.289** ✓ | **−0.677** ✓ |
| **medoid** (1 structure) | **−0.469** ✓ | **+0.323** ✗ | **+0.541** ✗ | **−0.310** ✓ | **−0.318** ✓ | **+0.165** ✗ | **−0.235** ✓ |
| **AVERAGED** (point cloud) | **−0.431** ✓ | +0.232 | **+0.485** ✗ | −0.248 | **−0.332** ✓ | +0.082 | **−0.266** ✓ |

✓ = i.i.d. CI excludes zero on the better side; ✗ = on the **worse** side.

**Legacy flips sign between the first and second rows** (−0.348 [−0.601, −0.098], 29W/13L, **5/5
folds** vs +0.323 [+0.141, +0.531], 11W/31L, **0/5 folds**) — my n=42 reproducing the coordinator's
n=126 to the third decimal. **AMBER is the only Hamiltonian that fails on every readout.**

*Medoid readout across α (the family the pillar emits):*

| α | disto | legacy | **amber** | d+l | d+a | l+a | d+l+a |
|---|---|---|---|---|---|---|---|
| 0.01 | **−0.648** ✓ | −0.008 | **+1.895** ✗ | **−0.627** ✓ | **−0.629** ✓ | −0.243 | **−0.683** ✓ |
| 0.05 | **−0.502** ✓ | +0.155 | **+1.383** ✗ | **−0.372** ✓ | **−0.443** ✓ | −0.244 | **−0.320** ✓ |
| 0.15 | **−0.469** ✓ | **+0.323** ✗ | **+0.541** ✗ | **−0.310** ✓ | **−0.318** ✓ | **+0.165** ✗ | **−0.235** ✓ |

AMBER at α=0.01: **+1.895 [+1.250, +2.512]**, se 0.321 (own MDE 0.899), median **+1.948**,
**10W/32L, 0/5 folds negative** — the effect is more than twice its own MDE. **Mean and median
agree**, so this is not a handful of catastrophic targets carrying an average; the median-vs-mean
early warning does not fire.

**And the cell the mandatory matrix actually asks about — does physics add to the structural
arm?** Each Hamiltonian minus `disto`, same α, same readout (AVERAGED, α=0.15): legacy
**+0.666 [+0.315, +1.037]** se 0.186 (13W/29L, **0/5 folds**), amber **+0.920 [+0.573, +1.289]**
se 0.185 (6W/36L, **0/5 folds**), d+l +0.189, d+a +0.092, l+a +0.497, d+l+a +0.168 — **every one
positive (worse).** The hybrids' intervals graze zero and sit below their own MDEs, so they are
**NOT MEASURED**, not "matched" — but not one of them has a point estimate on the helpful side.
**FALSIFIER F-A1 fires on the pool-restricted instrument.**

**Read.** AMBER is not merely uninformative — its low-energy tail is **actively anti-selective on
every readout**, and the effect is largest exactly where a CVaR selector concentrates (α → 0).
**Legacy is readout-conditional**: helpful at `member`, harmful at `medoid` and `AVERAGED`, and
worth nothing beyond the distogram in either case (§1.3's partial). This is the same object as
Sprint 20's argmin result, but it now covers the *whole tail family* at *every readout*, which is
what closes the door the sprint was opened to check: **the argmin→CVaR distinction does not
rescue AMBER, and it does not create a Legacy contribution that survives partialling.**

### 1.4b OPERATOR DECOMPOSITION — and the coordinator's original hypothesis, judged properly

Separating the two operators (n=42, W basis) settles the confounded primary in both
directions at once.

**(i) The TAIL, averaging held out — `medoid(tail) − argmin`** (n=42, with each contrast's own SE
and own MDE):

| α | disto | legacy | **amber** |
|---|---|---|---|
| 0.05 | −0.090 [−0.285,+0.089] (MDE 0.27) | −0.395 [−1.236,+0.376] (MDE 1.18) | −0.432 [−1.077,+0.175] (MDE 0.92) |
| 0.15 | −0.141 [−0.348,+0.058] (MDE 0.29) | −0.337 [−1.183,+0.445] (MDE 1.20) | **−1.344 [−1.979,−0.689]** (MDE 0.93), 32W/9L, 5/5 folds |
| 0.30 | +0.015 [−0.207,+0.227] (MDE 0.31) | −0.372 [−1.218,+0.408] (MDE 1.19) | **−1.862 [−2.555,−1.147]** (MDE 1.01), 31W/11L, 5/5 folds |

Only the AMBER rows at α ≥ 0.15 clear both their CI and their own MDE. **Legacy's rows are NOT
MEASURED**: an SE of 0.42 makes the smallest detectable effect 1.2 Å, far larger than anything
present, so the −0.34/−0.39 point estimates carry no weight in either direction.

**(ii) The AVERAGING, tail held fixed — `tail_avg − medoid(tail)`:** −0.18 to −0.35 Å on **every**
Hamiltonian at α=0.15, every CI excluding zero, every effect above its own MDE, 4–5/5 folds,
up to 37W/5L. It is a Hamiltonian-independent lever of about a third of an Ångström, and it is
the whole of the confounded original primary that was not the tail.

**So the coordinator's original prediction is CONFIRMED for AMBER, NOT MEASURED for Legacy — and
it did not matter.** A tail mean genuinely does beat a minimum on an energy that ranks badly, by up to
**1.9 Å** on AMBER, exactly as the "a minimum chases the worst-ranked outlier" argument said. It
buys **nothing** on the distogram (CI spans zero and the effect is under its own MDE at every
width), which the prediction also said.
But **the level it reaches is still below chance**: after the tail has rescued it,
AMBER's medoid readout is still **+0.541 [+0.244, +0.857]** *worse* than a matched-count random
subset at α=0.15 (11W/31L, 1/5 folds negative). **The tail operator converts a catastrophe into a
merely-worse-than-random selector.** That is why the matched random control, not argmin, had to be
the endpoint.

### 1.5 The trap the file named, and it did not spring

A tail arm that beats argmin proves nothing, because as α→1 the tail becomes the whole-pool
average — a zero-information operator that is already decent (**3.275 Å** at n=42, against a pool
mean of **4.357**). Every number in §1.4 is against the **matched-count random control at the same
count**, which is the only comparison that means anything, and the α=1 row is printed on every
table so the degeneracy is visible: at α=1 all seven Hamiltonians, the random control and the
ORACLE tail collapse to the identical 3.214.

### 1.5b THE READOUT DECIDES A HAMILTONIAN'S SIGN — coordinator, n=126

Measured by the coordinator at the full n=126 against matched-count random controls (α=0.15; same
pattern at 0.05 and 0.30, all CIs excluding zero):

| readout | disto | legacy |
|---|---|---|
| tail_member | −0.90 [−1.08, −0.73] | **−0.41 [−0.56, −0.26]** — BEATS its control |
| tail_medoid | −0.45 [−0.67, −0.24] | **+0.30 [+0.20, +0.40]** — WORSE, 35W/91L |
| tail_avg | −0.38 [−0.55, −0.22] | **+0.33 [+0.21, +0.45]** — WORSE, 41W/85L |

**Legacy's sign flips with the readout, and both readings are true.** A matrix quoted at
`tail_avg` reports Legacy as harmful; the same matrix at `tail_member` reports it as helpful.
So **no row in this document is quoted without its readout**, and "which Hamiltonian is best" is
reported *per readout*, never unqualified. My n=42 pool-restricted medoid numbers agree closely
with the coordinator's n=126 (legacy medoid +0.305 [+0.097, +0.530] vs +0.30 [+0.20, +0.40]),
which is the cross-check that lets me stand my §1.4 table down at n=42 rather than spend the
shared AMBER context finishing it.

**The mechanism is ERROR COHERENCE, not diversity.** Averaging cancels i.i.d. error and
*preserves* systematic error; Legacy's tail carries a coherent compactness bias (0.45 Å more
compact, 124W/2L — Sprint 20 L2c), so averaging preserves it and a single random member escapes
it. The first explanation offered — that Legacy's tail is low-*diversity* — was refuted by its own
sign control: **maximising set diversity is dead on all three bands and *minimising* it helps**
(−0.214 [−0.397, −0.040]). **No diversity term is built into any arm of my matrix**, and none
should be built on the strength of the earlier account.

**BANKED — this is a two-lane result.** The readout sign-flip reproduces across independent
implementations, different n, different code, to the third decimal:

| | audit lane (n=126) | this lane (n=42) |
|---|---|---|
| tail_member \| legacy | −0.407 [−0.557, −0.259] | **−0.406 [−0.682, −0.122]** |
| tail_medoid \| legacy | +0.302 [+0.204, +0.404] | **+0.305 [+0.097, +0.530]** |
| tail_member \| disto | −0.904 [−1.078, −0.729] | **−0.874 [−1.166, −0.572]** |

**And the pipeline identity now holds in three places** — my `tail0.15|disto`, the audit lane's
`tail_avg0.15|disto`, and production's `rmsd_avg` — at **max |diff| = 0.000000 on all 42 targets
this lane finished.** It is this lane's standing soundness gate.

### 1.6 The cross-readout matrix, and an identity that is *not* a discovery

Workstream D's near-theorem: CVaR's minimiser is a face supported on the α-tail, and the α-tail
of H **contains H's global pool minimum** — so for a pool-restricted selector whose training and
readout energy are the same H, CVaR-VQE and argmin have the **same optimal answer**.

**Verified, exactly**: `max |diagonal − argmin(train)| = 0.00e+00` at every α and on both bases.
That is an identity of the construction, printed as a check, and it is stated here as a check
rather than as a finding.

The informative cells are therefore the ones breaking its scope: **readout H ≠ training H**. All
7×7 (train, read) pairs are emitted at every α, single structure throughout (mean Å, n=42, W
basis, α=0.05):

| train ＼ read | disto | legacy | amber | d+l | d+a | l+a | d+l+a |
|---|---|---|---|---|---|---|---|
| disto | **3.246** | 3.186 | 3.365 | 3.178 | 3.235 | 3.209 | 3.233 |
| legacy | 3.260 | 4.234 | 4.011 | 3.293 | 3.378 | 3.989 | 3.387 |
| amber | 3.428 | 3.962 | **5.458** | 3.590 | 3.643 | 4.099 | 3.562 |
| d+l | 3.253 | 3.557 | 3.336 | 3.280 | **3.136** | 3.429 | 3.260 |
| d+l+a | 3.174 | 3.399 | 3.510 | 3.285 | 3.255 | 3.437 | **3.234** |

**What this says, and it is not what it looks like.** Reading out with a distogram-containing
energy recovers most of the damage a physics training Hamiltonian does (`amber→disto` 3.428 vs
`amber→amber` 5.458). But the best off-diagonal cell does **not** beat the best diagonal by
anything that survives its own selection: `d+l→d+a` 3.136 against `d+l+a→d+l+a` 3.234 is
**−0.098 [−0.243, +0.047]**, median **+0.000**, se 0.077, **own MDE 0.215** — the interval spans
zero *and* the effect is less than half the smallest effect this comparison could detect, so it is
**NOT MEASURED**. And it is a **max over 42 off-diagonal and 7 diagonal cells chosen on this
instrument**, which is a selection this test does not correct for. **Reported as a CEILING, never
as a method.**

**The durable statement is the simpler one: the training Hamiltonian's only effect is to restrict
the candidate set the readout sees, and every restriction by a physics energy is a net loss
against no restriction.**

### 1.7 Audit normalisations

At α=0.15, W basis, n=42 (mean Å, AVERAGED readout): declared **rank** vs **raw sum** vs
**z-score** — d+l 3.081 / 3.530 / 3.161; d+a 2.984 / 3.810 / **2.893**; l+a 3.389 / 3.789 / 3.556;
d+l+a 3.060 / 3.784 / 3.160. The declared rank transform beats the raw sum on all four by
0.40–0.83 Å and beats the z-score on three of four. **It loses to the z-score on `d+a` by
0.091 Å**, and that is recorded here rather than used to change the declaration: the pre-registered
choice stands, and this is the one cell where a non-declared form wins.

### 1.8 Ceilings, kept separate (BRIEF §6)

| ceiling | value (n=42, W basis, point cloud) |
|---|---|
| pool mean (zero-information reference) | 4.357 |
| whole-pool average (α=1, zero-information operator) | 3.275 |
| **best achieved by any native-free arm here** | **2.889** (disto, AVERAGED α=0.05) |
| ORACLE selector ceiling at α=0.05 (**not achieved**) | 1.618 |
| ORACLE pool best member (**generation ceiling, not achieved**) | 1.790 |

**Do not read the AVERAGED disto row as a CVaR result.** The coordinator's n=126 anchor is
decisive here: `tail_avg(0.15)|disto` = **3.0483**, matching the shipped `rmsd_avg` at
**max |diff| = 0.000000 on 126/126**, and `argmin|disto` = **3.454**, matching the instrument's
pinned constant. **At α=0.15 under the distogram, the tail-average selector *is* the production
pipeline** — the arm reproduces the shipped number because it is the shipped operator. That is a
soundness gate passed, not a result produced.

---

## 2. THE MANDATORY MATRIX — `s21/a_matrix.py`

> **COMPLETE.** `s21/results/a_matrix.json` — **12/12 targets, 12/12 complete rows, 106/106 cells
> per row, `complete=true`**, 9,385 s wall on one heavy process holding the AMBER context (now
> released). Report: `s21/results/a_matrix_report.txt`. The completion flag requires the **full**
> configuration — 7 Hamiltonians × 5 arms × their declared seeds × 4 readouts × 3 α, the 7×7
> cross-readout on all four variational seeds, and both audit normalisations — so it cannot
> certify a subset.

Seven Hamiltonians, genuine CVaR-VQE as the selector, every component separably evaluable.

> **SCOPE, stated before the numbers so it cannot be read as an excuse afterwards.** This matrix
> is a **comparison among Hamiltonians at a fixed small budget**, not an attempt on the incumbent.
> Sprint 20 Q3 established that **no sampler, quantum or classical, at 8192 evaluations beat the
> zero-evaluation retrieval pool (3.230)**; at 512 evaluations every arm here will sit well above
> the 3.204 Å incumbent, and the smoke confirms it (3.3–4.5 Å). **No number in this section is a
> candidate for the ladder.** What the matrix measures is the *relative* standing of seven
> Hamiltonians under an identical selector, identical budget and identical seeds — which is
> exactly what §4 of the brief asks for and nothing more.

**One definition differs between §1 and §2 and must not be conflated.** The structural
Hamiltonian is the shipped **Bayes-risk** distogram score in §1 (`I.shipped_score`, which is the
production top-75 filter, verified 75/75) and Sprint 20's **inverse-variance weighted squared
deviation** `Σ (d − d̂)²/σ²` in §2 (`s20.qb2_lib.Ham("DIST")`). Both are functions of the same
leave-fold-out distogram; §2 uses Sprint 20's so the VQE stack and every Sprint-20 comparison
transfer unmodified. Declared, not discovered.

- **Selector**: `s20.qb2_opt.arm_vqe` — `core.quantum.MPSAnsatz` (RY + CNOT chain + final RY)
  simulated exactly, with the analytic CVaR score-function gradient at `baseline="const"` (the
  corrected estimator, never the recorded `tail` defect). The bitstring picks a per-residue
  conformer basin and **the angle is drawn continuously from that basin's von Mises component** —
  a continuous density with full support on the torus. **No lattice.**
- **Readout**: four of them, at each of α ∈ {0.05, 0.15, 0.30}, because §1.5b shows a
  Hamiltonian's sign is a function of the readout — `tail_min` (= `Field.best_z`, the argmin of
  the training H over everything seen, and the operator the deployed pillar actually has),
  `tail_member`, `tail_medoid`, and `tail_avg` (**AVERAGED**, labelled everywhere). **CVaR is the
  training objective; the readout is an argmin.** Each readout carries **its own matched-count
  random control drawn from that arm's own generated set.**
- **Budget**: matched in **candidate evaluations**, which is the space the selector works in — a
  three-component hybrid gets no extra look at the space, only a different opinion about what it
  saw. Wall-clock is reported separately and is *not* matched.
- **Seeds**: 4 on every variational arm (Sprint 20: within-target seed sd 0.200 Å = 2.4× MDE).
- **Controls**: `vqe_untrained` (best-of-N from the **untrained** circuit, θ never stepped — never
  an initialisation mean), `best_of_N` (i.i.d. draws from the same basin mixtures, zero
  optimisation), `metro` (matched-move-class continuous Metropolis), `helix` (zero-information but
  plausible; uniform-on-the-torus is *not* a zero-information control).

**Normalisation for the hybrid cells, declared before RMSD**: per component,
`z = (cond(e) − median_pool) / (IQR_pool/1.349)`, with `cond` = signed log about the pool median
for AMBER and the identity for Distance and Legacy, all fitted free and native-free on the
target's own K=500 pool. **The property that makes this defensible is an identity, not a result:**
`cond` is strictly monotone, so for a **single-component** Hamiltonian every argmin readout is
*exactly* unchanged — conditioning cannot flip a single-Hamiltonian ranking. The normalisation is
therefore a decision about the **hybrid cells only**. Audits: `robustz` (no conditioning, Sprint
20's primary) and `rank` (the form tailprice declared).

**Readout column, per the coordinator's §1.5b finding.** Every arm emits four readouts at each of
α ∈ {0.05, 0.15, 0.30} — `tail_min` (≡ argmin, by the identity), `tail_member`, `tail_medoid`,
`tail_avg` (AVERAGED, labelled) — **each with its own matched-count random control drawn from that
arm's own generated set**, because a control must be matched in the space the operator works in.
No row is reported without its readout.

**Cross-readout is where the AMBER budget goes**, on the coordinator's instruction and for a
reason my own §1.6 confirms: the pool-restricted, same-H, order-based case is closed analytically
and verified at `tail_min == argmin`, max difference exactly 0, so the only cells that can be
informative are those breaking that scope — **readout H ≠ training H**. All 7×7 (train, read)
pairs are emitted, **on all four variational seeds**, at every α.

**OPERATOR FORKS (BRIEF §7 rule 0), with the alternative not taken named.** This lane has a
directional hypothesis (F-A1 expects physics to fail), so the forks are enumerated in the module
docstring and repeated here, each with the direction it would have pushed:

| fork | taken | **not taken** | direction it would have pushed |
|---|---|---|---|
| functional | AMBER as a **bare single point** | `E ∘ Relax₅₀`, the deployed object | **toward my conclusion** — relaxation would make AMBER less catastrophic. The fork to distrust most. |
| basis | energies on the rebuild, RMSD on `build_ca_exact` | scoring the distogram on the rebuild too (54/75 vs 75/75) | unknown — so it was **measured, not chosen**: the pool run carries both bases and they agree to ~0.06 Å |
| readout | **four readouts, none privileged** | quoting `tail_avg` alone | **toward my conclusion** — `tail_avg` is exactly where Legacy looks worst |
| normalisation | signed-log + robust z on the target's own pool | raw sum; pool-quantile rank | **toward my conclusion** — a raw sum is an AMBER-outlier detector. Both computed as declared audits |
| null | matched-count random from the **arm's own** generated set | uniform-on-the-torus; the initialisation **mean** | both forbidden by BRIEF §7 rule 4 and the untrained-circuit rule |
| budget | matched in **candidate evaluations** | matching **wall clock** | **AGAINST my conclusion** — it would give single-component Hamiltonians ~3× more candidates and favour Legacy-only/AMBER-only |
| **seed count** *(FORK I OMITTED — added on review)* | **4 seeds** on variational arms, 2 on controls | more seeds; or quoting best-of-seeds | **the fork that most constrains the table.** Seed sd is 0.51–1.06 Å and best-of-4 beats the 4-seed mean by 0.65–1.36 Å — larger than most row-to-row gaps. I disclosed this in §2.3 but did **not** enumerate it as a fork, which is the omission. Consequence: only the DIRECTION is safe; **every row-to-row gap sits inside seed noise** |

I designed these comparisons and have a stake in their direction, so this list is offered for
**independent fork review** rather than treated as sufficient on its own.

**Gates, each reporting how many times it FIRED** — all three PASS:

| gate | what it asserts | result | times fired |
|---|---|---|---|
| **G-A1** | `AmberSP` bit-exact against `core.amber.refine_coords(k=0, steps=−1)` | max relative difference **0.000e+00** | **12** |
| **G-A2** | the conditioning is exactly order-preserving on real pool energies | **0** rank breaks | **576** |
| **G-A3** | a single-component `MultiHam`'s argmin equals the genuine energy's argmin under all three normalisations | **0** breaks | **27** |
| **identity** | `tail_min` ≡ argmin over everything seen, on the VQE's own generated sets | max \|diff\| **0.000e+00** | **56 cells** |

G-A3 is the gate the normalisation argument leans on, and it is not vacuous: it fired 27 times
across three components × three normalisations × three targets and never once found a break.

**One guard DID pass vacuously and is reported as such.** `Field`'s non-finite sentinel
substitution fired **0 times over 2,048 evaluations per Hamiltonian** — AMBER's single point is
finite on every continuous torsion sample the VQE and its controls drew. That is consistent with
§1.3 (finite fraction 1.000 on the K=500 pool) and with Sprint 20 L5c: on plausible geometry AMBER
is *finite-but-meaningless*, and nothing throws. **A guard that never fires is not evidence that
the hazard is absent** — it is evidence that this sampler never visited it.

---

### 2.1 THE MANDATORY MATRIX — final Cα-RMSD, the endpoint

n = 12 targets, BUILT-CHAIN basis (single structure from continuous torsions via
`build_ca_exact`), **readout = argmin of the training Hamiltonian over everything seen**, every
variational cell the mean over **4 ansatz seeds**. ORACLE references, never achieved: retrieval
pool mean 4.919, pool best member 1.535.

| Hamiltonian | **VQE** mean/median | vqe_untrained | best_of_N | metro | helix (zero-info) |
|---|---|---|---|---|---|
| Legacy | 4.767 / 4.775 | 4.817 | 4.856 | 4.832 | 5.176 |
| **AMBER** | **5.316 / 5.161** | 5.524 | 5.687 | 5.378 | 5.185 |
| Leg+Amb | 4.789 / 4.504 | 4.721 | 4.703 | 4.918 | 5.181 |
| **Distance** | **3.845 / 3.686** | 4.079 | 3.890 | 3.751 | 4.295 |
| Dist+Leg | 4.119 / 3.790 | 4.480 | 4.301 | 4.095 | 4.956 |
| Dist+Amb | 4.118 / 4.067 | 4.238 | 4.252 | 3.924 | 4.496 |
| Dist+Leg+Amb | 4.406 / 4.265 | 4.298 | 4.257 | 4.269 | 5.089 |

**Distance-only wins the matrix on both mean and median. All six Hamiltonians containing Legacy
or AMBER are worse.** Mean and median agree throughout, so no row is an outlier artefact.

### 2.2 Target-level paired differences, CI, W/L, folds — before any diagnostic

**vs the structural-only reference (Distance, VQE).** Negative would mean physics helps:

| Hamiltonian | mean | median | se | own MDE | i.i.d. CI | fold CI | W/L |
|---|---|---|---|---|---|---|---|
| Legacy | **+0.922** | +0.949 | 0.396 | 1.110 | [+0.180, +1.642] | [+0.010, +1.834] | 4W/8L |
| **AMBER** | **+1.471** | +1.616 | 0.630 | 1.764 | [+0.266, +2.625] | [+0.378, +2.601] | **2W/10L** |
| Leg+Amb | +0.944 | +0.859 | 0.528 | 1.478 | [−0.026, +1.951] | [−0.249, +1.969] | 4W/8L |
| Dist+Leg | +0.274 | +0.179 | 0.253 | 0.708 | [−0.190, +0.748] | [−0.379, +0.854] | 5W/7L |
| Dist+Amb | +0.273 | +0.149 | 0.140 | 0.391 | [+0.021, +0.552] | [−0.021, +0.561] | **2W/10L** |
| Dist+Leg+Amb | +0.561 | +0.313 | 0.292 | 0.819 | [+0.013, +1.114] | [−0.144, +1.278] | 3W/9L |

**All six point estimates are on the harmful side; four have i.i.d. CIs excluding zero.** But
**each of those four has an effect smaller than its own MDE** — the winner's-curse flag:
significant *and* underpowered means the point estimate is very likely inflated. The honest
reading is **consistently harmful in direction, magnitude NOT MEASURED at n=12**. Per-fold means
carry the same sign structure: AMBER positive on 4/5 folds (+3.13, +1.25, +1.68, −0.14, +1.41),
Legacy on 4/5.

**FALSIFIER F-A1 fires on the mandatory matrix as well as on the pool.**

**MANDATORY CONTROL — VQE minus best-of-N from the UNTRAINED circuit** (θ never stepped; not an
initialisation mean). Every cell **NOT MEASURED**: Legacy −0.050, AMBER −0.208, Leg+Amb +0.068,
Distance −0.234 [−0.512, +0.033], Dist+Amb −0.120, Dist+Leg+Amb +0.108. Only Dist+Leg reaches
−0.361 [−0.663, −0.093] (9W/3L, fold CI excludes zero), and it is still under its own 0.423 MDE.

> **FALSIFIER F-A2 is not passed at this budget, and that is a negative against my own pillar.**
> Sprint 20 Q4 (−0.210 [−0.339, −0.082], 5/5 folds) was measured at **8192** evaluations; at
> **512** the same comparison is NOT MEASURED on every Hamiltonian. This is a **scope restriction
> on Q4, not a refutation of it** — 16× less budget and n=12 rather than 126. Recorded as
> **OPEN**, never as a contradiction.

**MATCHED CLASSICAL CONTROLS.** VQE − `best_of_N`: all seven NOT MEASURED (−0.371 to +0.148).
VQE − `metro`: all seven NOT MEASURED (−0.129 to +0.194). **F-A3's condition holds — at this
budget the quantum arm is reproduced by simpler classical baselines.**

**ZERO-INFORMATION reference** (constant α-helix + matched torsion noise, same H, same budget).
As a *reference* it is a null; as a *gate* it is the worst arm measured. The VQE beats it on six
of seven, but only `Dist+Leg` (−0.837 [−1.424, −0.272], 9W/3L) has a CI excluding zero — and
**AMBER does not beat the helix at all** (+0.131 [−0.681, +1.035], 6W/6L). **A genuine CVaR-VQE
on genuine AMBER is indistinguishable from a constant α-helix with noise.**

### 2.3 Ansatz-seed sensitivity — why 4 seeds was the floor, and why it was not enough

| Hamiltonian | within-target seed sd | untrained | best-of-4-seeds | 4-seed mean |
|---|---|---|---|---|
| Legacy | 0.649 | 0.630 | 3.873 | 4.767 |
| **AMBER** | **1.062** | 1.147 | 3.952 | 5.316 |
| Distance | 0.539 | 0.460 | 3.192 | 3.845 |
| Dist+Leg+Amb | 0.511 | 0.571 | 3.698 | 4.406 |

**Seed sd here is 0.51–1.06 Å — two to five times Sprint 20's 0.200 Å**, because the budget is
16× smaller. **The whole matrix is seed-noise-limited**: the gap between best-of-4-seeds and the
4-seed mean (0.65–1.36 Å) exceeds most of the Hamiltonian differences it is being used to
resolve. That is the strongest single reason the §2.2 magnitudes are NOT MEASURED, and it is why
quoting a best-seed number here would be indefensible.

### 2.4 THE READOUT COLUMN — the pool result replicates on generated candidates

α = 0.15, VQE arm, each readout against **its own** matched-count random control drawn from
**that arm's own generated set**:

| Hamiltonian | member | medoid | AVERAGED |
|---|---|---|---|
| **Legacy** | **−0.368 [−0.554, −0.196], 11W/1L — BETTER** | −0.321 [−0.571, −0.091], under own MDE | −0.133 [−0.326, +0.038] NOT MEASURED |
| **AMBER** | +0.087 [−0.279, +0.420], 3W/9L | **+0.335** fold[+0.023, +0.626], 4W/8L | **+0.368** fold[+0.014, +0.722], 3W/9L |
| Distance | **−0.794**, 10W/2L — BETTER | −0.618, 9W/3L | −0.582, 8W/4L |
| Dist+Leg | **−0.628**, 9W/3L — BETTER | −0.509, 10W/2L | −0.411, 9W/3L |
| Dist+Amb | **−0.677**, 11W/1L — BETTER | **−0.502** — BETTER | **−0.556** — BETTER |
| Dist+Leg+Amb | **−0.558**, 10W/2L — BETTER | **−0.518** — BETTER | −0.410, 8W/4L |

**A third independent instrument for the readout effect, on a candidate distribution that is not
the pool.** Legacy's advantage is largest at `member` (−0.368, above its own MDE), halves at
`medoid`, and vanishes at `AVERAGED` — the same monotone ordering the coordinator measured at
n=126 and I measured at n=42, now reproduced on VQE-generated candidates. **AMBER is the mirror
image: positive on all three readouts, with fold CIs excluding zero on the worse side at medoid
and averaged.**

### 2.5 CROSS-READOUT — the coordinator's priority cell, answered with a negative

Rows = training H (defines the candidate set), columns = readout H (argmin inside it), single
structure, all 4 variational seeds:

| train ＼ read | Legacy | AMBER | Leg+Amb | **Distance** | Dist+Leg | Dist+Amb | D+L+A |
|---|---|---|---|---|---|---|---|
| Legacy | 4.767 | 5.385 | 4.869 | 4.070 | 4.369 | 4.036 | 4.354 |
| AMBER | 5.073 | 5.316 | 4.924 | 4.202 | 4.403 | 4.213 | 4.229 |
| Leg+Amb | 4.916 | 5.504 | 4.789 | 4.198 | 4.326 | 4.121 | 4.115 |
| **Distance** | 4.918 | 5.118 | 4.841 | **3.845** | 4.239 | 3.922 | 4.248 |
| Dist+Leg | 4.761 | 5.120 | 4.930 | 4.049 | 4.119 | 3.942 | 4.226 |
| Dist+Amb | 4.966 | 5.256 | 4.786 | 4.045 | 4.243 | 4.118 | 4.297 |
| Dist+Leg+Amb | 4.851 | 5.359 | 4.849 | 4.031 | 4.371 | 3.989 | 4.406 |

**IDENTITY CHECK** (a check, never a discovery): `max |diagonal − argmin(train)| = 0.000e+00`,
confirming on generated sets what §1.6 confirmed on the pool.

**Best off-diagonal vs best diagonal: `Distance→Dist+Amb` (3.922) against `Distance→Distance`
(3.845) is +0.077 [−0.187, +0.348], 5W/7L — the off-diagonal is WORSE, and NOT MEASURED.**
Changing the readout Hamiltonian away from the best single Hamiltonian buys nothing here. The
column structure is the whole story: **every training Hamiltonian read out with Distance lands at
4.03–4.20, and Distance→Distance alone reaches 3.845.** Restricting the candidate set by physics
and then reading with the distogram is worse than the distogram alone.

> **Scope, stated plainly, because this looks like a contradiction and is not one.** Workstream B
> measures changing the READOUT Hamiltonian at **−0.697 [−1.059, −0.352], 5/5 folds**, against a
> null for changing the TRAINING Hamiltonian. My matrix does not contradict that: B compares
> readouts *across* Hamiltonians on a fixed generated set, whereas my off-diagonal is measured
> against **the best diagonal**, which already carries the best readout available. Both are true
> — a readout swap is worth a great deal relative to a *bad* readout and nothing relative to the
> *best* one. Neither result displaces the other.

### 2.6 Explanation — objective convergence, gradient, CVaR behaviour

`d_cvar` is negative on **all seven** Hamiltonians (−0.149 to −0.412): **the VQE genuinely
optimises its objective on every one.** Gradient norms 0.26–1.47, ESS/shots 0.41–0.60,
cos(g_α, g_α=1) 0.21–0.64, |Δθ| 4.3–4.6 — the machinery moves and the CVaR tail is not degenerate.

**And it does not convert.** ρ(objective gained, RMSD gained), trained minus untrained: Distance
**+0.035**, Dist+Amb **−0.476**, Dist+Leg+Amb −0.112, Legacy +0.182, Leg+Amb +0.441, Dist+Leg
+0.413, AMBER +0.762. **On the Hamiltonian that actually works, optimising the objective harder is
uncorrelated with getting a better structure.** A sixth independent instrument for *search
saturates, discrimination binds*.

**The sentinel guard fired 0 times over 24,576 evaluations per Hamiltonian** — reported as
vacuous, per BRIEF §7.3.

### 2.7 Normalisation audit — a defect in my own table, then the rule applied

**First, the defect, because it changes the answer.** The audit arms run at **seed 0 only**. The
version of this table I first published put the declared form's **four-seed mean** beside them —
an **unmatched comparison across the seed dimension**, on an instrument whose seed sd is
0.51–1.06 Å. Corrected, with all three forms read on the identical circuit:

| Hamiltonian | **declared** @s0 | robustz @s0 | rank @s0 | winner | *(declared 4-seed, not the comparator)* |
|---|---|---|---|---|---|
| Leg+Amb | **4.757** | 4.880 | 5.194 | declared | 4.789 |
| Dist+Leg | **3.866** | 4.151 | 4.494 | declared | 4.119 |
| Dist+Amb | **4.080** | 4.086 | 4.198 | declared | 4.118 |
| Dist+Leg+Amb | 4.363 | 4.503 | **4.224** | **rank** | 4.406 |

**The unmatched table cost a cell.** It showed `robustz` beating the declared form on `Dist+Amb`
(4.086 vs 4.118); matched at seed 0 the declared form wins that cell (4.080 vs 4.086). So the
count is **declared 3 of 4, not 2 of 4** — a correction that happens to favour my own declared
choice, which is exactly why it is stated first and in full.

**Now the rule, applied.** I wrote: *"If a non-declared form wins, the declared choice is recorded
as WRONG, not quietly swapped."* On **`Dist+Leg+Amb` the `rank` transform wins by 0.139 Å, so the
declared signed-log + robust-z normalisation is recorded as WRONG on that cell.** It is not
swapped: the pre-registered form stands for every other cell and no alternative is promoted, since
0.139 Å is far inside the 0.51–1.06 Å seed noise and is therefore **NOT MEASURED**. Invoking low
power after the fact does not release me from a rule I declared before the run — the cell is
labelled wrong *and* the difference is labelled unresolvable, and both are true.

**The fork is not load-bearing, and this is the strong form of the answer.** Taking the **best**
normalisation on every hybrid cell — 4.757, 3.866, 4.080, 4.224, minimum **3.866** — still does
not reach **Distance-only at the same seed, 3.807**, which carries no normalisation at all.

> **The margin is 0.059 Å, not the ~0.24 Å the fork review inferred from my unmatched table.** The
> review's conclusion survives; its margin does not. 0.059 Å is deep inside seed noise, so the
> correct label is **NOT MEASURED in either direction** — Distance is not shown to beat the
> best-normalised hybrid at this budget, it is merely not beaten by it. The verdict of §2.1 rests
> on Distance beating the hybrids **as declared** (0.27–0.94 Å), not on this reconstruction.


---

### 2.8 SCOPE OF THE AMBER MAGNITUDE — the fork I flagged was worse than I said

**`H_AMBER` here is the BARE single point. My +1.471 Å is NOT the deployed Hamiltonian's cost.**

The deployed object is `E ∘ Relax₅₀` — restrain at k=100, minimise, report the *unrestrained*
energy at the minimised coordinates — and Sprint 20 L7c established that the relaxation is
**constitutive**, not cosmetic. Sprint 20 measured deployed AMBER's own cost at **+0.021
[+0.014, +0.028]**, and at k=30 it *beats* the projection by **−0.022**. **The deployed AMBER is
roughly two orders of magnitude less catastrophic than the bare one.**

So every AMBER magnitude in §2 — +1.471 vs Distance, +0.335/+0.368 vs its readout controls, the
5.316 cell — is a property of the **bare single point** and **does not transfer to the deployed
Hamiltonian**. Anyone quoting +1.471 as "AMBER's cost in the pipeline" is quoting the wrong
object. This was the fork I named as the one to distrust most, and naming it is what stopped it
being quoted wrongly.

**And a correction to how far the direction survives it, which cuts against me.** The fork review
states that "the four Legacy-containing rows never touch the AMBER functional". That is not right,
and the artefact settles it: `Leg+Amb = (LEG, AMB)` and `Dist+Leg+Amb = (DIST, LEG, AMB)` both
contain AMBER. Of the six non-reference rows, only **two are AMBER-free**:

| row | components | AMBER-free? | vs Distance |
|---|---|---|---|
| **Legacy** | (LEG,) | **yes** | **+0.922** |
| **Dist+Leg** | (DIST, LEG) | **yes** | **+0.274** |
| AMBER | (AMB,) | no | +1.471 |
| Leg+Amb | (LEG, AMB) | no | +0.944 |
| Dist+Amb | (DIST, AMB) | no | +0.273 |
| Dist+Leg+Amb | (DIST, LEG, AMB) | no | +0.561 |

**So the 6/6 direction is not 6/6 independent of the AMBER functional fork — it is 2/2 on the
AMBER-free rows and 4/4 on rows that carry it.** The honest statement is:

- **Independent of the AMBER fork:** Legacy and Dist+Leg are both harmful (+0.922, +0.274). Two
  rows, both on the harmful side, and Legacy's is one of the two largest effects in the table.
- **Conditional on the bare-AMBER fork:** the remaining four rows, whose *magnitudes* would
  shrink substantially under `E ∘ Relax₅₀` and whose *signs* are untested under it.

That is a weaker claim than the review credited me with, and it is the correct one. **The untested
regime — the whole matrix re-run with `H_AMBER = E ∘ Relax₅₀` — is recorded as OPEN (O7).**


---

## 3. PROCESS NOTES

**The shared box is a measurement condition, not a nuisance.** The n=126 tail-pricing run was
**killed at target 13/126** by `core.amber.memory_guard` ("physical memory at 96% exceeds the 92%
ceiling") because the box was full. The guard is correct and was not bypassed; the caller was made
resilient (hold for memory, retry with backoff). Per project memory the box fits **two** heavy
jobs; with three or more, OpenMM context creation is refused. Any timing in this document is from
a contended box and is **not** a performance measurement (BRIEF §10).

**And a share of that contention was mine.** When the coordinator reported four OpenMM contexts
competing — one lane measuring **52 ms per AMBER single point against a 6 ms nominal, a 9×
slowdown** — I checked my own processes and found **three of mine alive**: two detached
`s21.tailprice` runs that a `pkill -f` had failed to reap, plus the matrix. I had believed I was
running one. **`pkill -f` returning success is not evidence that the process died** — the check is
`ps`, and I did not make it until asked. The strays were killed and the lane is now one process.
Nine-fold contention on a shared box is not a fair trade for a second lane of my own results.

**Self-corrections recorded**: the stale-artefact resume defect (§1.1) was mine to catch and would
have silently reduced n; the completion-flag defect was flagged by the coordinator and fixed to
require the full configuration; the confounded primary was re-primaried by the coordinator before
any table was read, and the original pre-registration is kept unedited with the amendment logged
beneath it.

---

---

## 3b. UNTESTED REGIMES, recorded as OPEN

Per BRIEF §8: a pre-registration that could not be executed in full says so, keeps its text
unedited, and records what it did not reach.

| # | regime | why it is open | what would close it |
|---|---|---|---|
| O1 | **Tail pricing at the pre-registered n=126** | paused at n=42 to hand the shared AMBER context to the mandatory matrix | `python -m s21.tailprice` — resumes; all 42 rows verified to survive the completeness filter |
| O2 | **The matrix at n > 12** | 52 ms/AMBER point under contention vs a 6 ms nominal | more box time; the per-cell configuration needs no change |
| O3 | **Budget > 512 evaluations per arm** | not a scientific choice — the matrix is a Hamiltonian comparison, not an attempt on the incumbent | a budget ladder; Sprint 20 Q3 predicts no arm crosses the pool at any budget |
| O4 | **The PARTIAL of each component on the VQE's *generated* set** | §1.3's partial is measured on the K=500 **pool**; the matrix records each component's marginal on generated candidates but not the cross-component correlations needed to partial | one extra recorded coefficient per cell; the matrix's Legacy rank-skill column is therefore reported as a **marginal only** |
| O5 | **Legacy's `member`-readout advantage, causally** | measured twice (n=126 and n=42) and reproduced to the third decimal, but the error-coherence mechanism is **cited from Sprint 20 L2c, not re-measured here** | a compactness/coherence decomposition of the member vs medoid gap |
| O7 | **The whole matrix with `H_AMBER = E ∘ Relax₅₀`** | §2 uses the BARE single point; the deployed object relaxes first, and Sprint 20 prices its own cost at +0.021 [+0.014, +0.028] — two orders of magnitude smaller | re-run `a_matrix.py` with the relaxed functional. **Until then the four AMBER-containing rows' magnitudes and signs are untested under the deployed Hamiltonian** |
| O6 | **AMBER's non-monotonicity, directly** | inferred from ρ≈0 marginal *plus* a tail significantly worse than random; the radius-of-gyration decomposition was not recorded | add `rg` to the tail-pricing row and regress tail membership on it |

---

## 4. CLAIMS

Every structural row states its **readout** and its **basis**, because §1.5b shows a Hamiltonian's
sign is a function of the readout.

| # | claim | label | evidence |
|---|---|---|---|
| A1 | **AMBER's α-tail is significantly WORSE than a matched-count random subset of the same size**, on every readout, at every α ≤ 0.15, on both bases. Medoid readout: **+1.895 [+1.250, +2.512]** se 0.321 (median +1.948, 10W/32L, **0/5 folds negative**) at α=0.01; **+1.383** at α=0.05; **+0.541 [+0.244, +0.857]** (11W/31L) at α=0.15. AVERAGED readout: +1.530 / +0.998 / +0.485. MEMBER readout: **+0.605 [+0.364, +0.866]** (11W/31L, 0/5 folds) — it fails on all three. Mean and median agree, so it is not a few catastrophic targets. | ESTABLISHED (n=42 of 126, W and RB bases) | `_PARTIAL_tailprice_n42.json`, `_PARTIAL_tailprice_report_n42.txt` |
| A2 | **Legacy's α-tail is at or below chance at the medoid and averaged readouts** — CI spans zero at α ≤ 0.05 with the point estimate on the wrong side, and excludes zero on the **worse** side at α=0.15–0.30 (medoid +0.305 [+0.097, +0.530]; averaged +0.243 [+0.000, +0.485], 1/5 folds negative). **But its sign flips at the `member` readout**, where the coordinator measures −0.41 [−0.56, −0.26] at n=126. **Both are true; the readout decides.** | ESTABLISHED, readout-conditional | ibid.; coordinator n=126 |
| A3 | **No Hamiltonian containing Legacy or AMBER beats the structural-only arm.** Minus `disto` at α=0.05, AVERAGED: legacy +0.646 [+0.258, +1.035] (8W/27L, 0/5 folds), amber +1.511 [+0.984, +2.060] (3W/32L, 0/5 folds), and every hybrid positive with 1/5 folds negative. | ESTABLISHED (pool-restricted) | ibid. |
| A4 | **The pre-registered falsifier F-A1 fires on the pool-restricted instrument.** A CVaR-VQE whose selection is over this pool is bounded at chance on Legacy and AMBER regardless of optimisation quality. | SUPPORTED | PREREG_A stage 0 |
| A5 | **The argmin→CVaR distinction genuinely helps AMBER and still does not rescue it.** With the averaging confound removed, `medoid(tail) − argmin` on AMBER is **−1.597 [−2.304, −0.850]** at α=0.15 and **−2.094 [−2.856, −1.306]** at α=0.30 (5/5 folds) — the coordinator's original prediction, confirmed. The level reached is still **+0.518 [+0.204, +0.868]** worse than matched random. **The tail converts a catastrophe into a merely-worse-than-random selector.** | ESTABLISHED | §1.4b |
| A6 | **For a selector whose training and readout energy are the same H, the α-tail's minimum IS the global argmin** — `max |tail_min − argmin| = 0.00e+00`, on the pool (all α, both bases, n=42) and on the VQE's own generated sets (56/56 smoke cells, and 0.000e+00 again in the live matrix). An identity of the construction, printed as a check, never as a discovery. | **EXACT** (verified twice) | §1.6; `_SMOKE_a_matrix.json` |
| A7 | **The basis asymmetry is not the explanation.** Rebuild − window shift **+0.011 Å**, ρ(d_W, d_RB) = **0.989**, and the RB table reproduces the W table on every cell. AMBER's failure survives being scored entirely inside its own geometry. | ESTABLISHED | ibid. |
| A8 | **AMBER carries no rank information about accuracy**: marginal Spearman(E_AMBER, d_true) = **−0.042**, fold CI [−0.116, +0.032], right sign on 19/42 targets. Its tail being *worse* than random therefore means the relation is **non-monotone**, not merely absent. | ESTABLISHED | ibid. |
| A8a | **Legacy carries no in-band rank information BEYOND the distogram.** Marginal +0.328 fold[+0.301, +0.355] collapses to a PARTIAL of **−0.0076 fold[−0.0616, +0.0681]**, 25/42 negative, once the deployed distogram is partialled out — **−2% of its marginal retained**. AMBER's partial (−0.029) is the control showing the partialling manufactures nothing. This is the mechanism for `d+l` never beating `d`. | ESTABLISHED | §1.3 |
| A8b | **I reported Legacy's marginal as its contribution, and my own artefact contained the refutation three lines above it** (`rho(disto, legacy)` = +0.478). A marginal correlation with truth, for a score correlated with the score in production, is not that score's contribution. Same class as Sprint 20 Q11 and the shared-referent floor. | RETRACTED (my reading) / recorded as process | §1.3 |
| A8c | **A legend in my own report was inverted** — "negative = lower energy is better" on a block where the working score is +0.618 with 39/42 positive. As printed it said the arm that works is the arm that fails. Fixed. | recorded as process | §1.3 |
| A8d | **The readout sign-flip replicates across two independent lanes to the third decimal** (tail_member\|legacy −0.407 vs −0.406; tail_medoid\|legacy +0.302 vs +0.305; tail_member\|disto −0.904 vs −0.874), and the pipeline identity `tail_avg(0.15)\|disto ≡ rmsd_avg` holds at max\|diff\| = 0.000000 in three places. | ESTABLISHED, two-lane | §1.5b |
| A9 | **AMBER is finite everywhere on the K=500 pool** (finite fraction 1.000) with a ~16-decade per-target range and a median **58.2%** of members above 1e4 kcal/mol — Sprint 20 L5c (53.5% on the already-filtered top-75) reproduced and exceeded on the unfiltered pool. | ESTABLISHED | ibid. |
| A10 | **A completion flag that counts skipped rows certifies failure as success**, and a resume path that trusts it silently shrinks n. Both were live in this lane's own artefact; both fixed. | recorded as process | §1.1 |
| A11 | **`pkill -f` returning success is not evidence a process died.** Three of my processes were alive while I believed one was; the check is `ps`. This lane contributed to a measured 9× AMBER slowdown on the shared box. | recorded as process | §3 |
| A12 | **THE MANDATORY MATRIX, n=12, genuine CVaR-VQE, 4 seeds, matched budget, argmin readout: Distance-only 3.845/3.686 wins on mean and median; all six Hamiltonians containing Legacy or AMBER are worse** (Legacy +0.922, AMBER +1.471 with 2W/10L, Leg+Amb +0.944, Dist+Leg +0.274, Dist+Amb +0.273 with 2W/10L, D+L+A +0.561). Six of six on the harmful side; four i.i.d. CIs exclude zero but **each effect is under its own MDE**, so the DIRECTION is established and every MAGNITUDE is NOT MEASURED. | ESTABLISHED (direction) / NOT MEASURED (magnitude) | `a_matrix.json` complete=true, `a_matrix_report.txt` §1–2 |
| A12a | **SCOPE ON A12, and it is narrower than the fork review credited.** `H_AMBER` is the **bare single point**, not the deployed `E ∘ Relax₅₀` (Sprint 20: deployed cost +0.021 [+0.014, +0.028]). **+1.471 Å is not the pipeline's AMBER cost.** Only **two** of the six non-reference rows are AMBER-free — Legacy +0.922 and Dist+Leg +0.274 — not four. **The direction is 2/2 independent of the AMBER functional fork and 4/4 conditional on it.** | ESTABLISHED (the scope) / OPEN (O7) | §2.8 |
| A13 | **A genuine CVaR-VQE on genuine AMBER is indistinguishable from a constant α-helix with noise** — +0.131 [−0.681, +1.035], 6W/6L, at matched budget on the same Hamiltonian. Every other Hamiltonian beats the zero-information reference in direction. | ESTABLISHED | §2.2 |
| A14 | **FALSIFIER F-A2 is not passed at 512 evaluations**: VQE minus best-of-N from the UNTRAINED circuit is NOT MEASURED on all seven Hamiltonians (−0.361 to +0.108). This is a **scope restriction on Sprint 20 Q4** (measured at 8192 evaluations, n=126), **not a refutation of it**, and the untested regime is OPEN. | NOT MEASURED | §2.2 |
| A15 | **F-A3's condition holds at this budget**: VQE minus `best_of_N` and VQE minus `metro` are NOT MEASURED on all seven Hamiltonians. The quantum arm is reproduced by simpler classical baselines here. | NOT MEASURED / kill-rule condition met | §2.2 |
| A16 | **The matrix is seed-noise-limited.** Within-target ansatz-seed sd **0.51–1.06 Å**, two to five times Sprint 20's 0.200 Å at 16× the budget; best-of-4-seeds beats the 4-seed mean by 0.65–1.36 Å, which exceeds most Hamiltonian differences. Four seeds was the declared floor and it was **not enough**. | ESTABLISHED | §2.3 |
| A17 | **The readout effect replicates on VQE-GENERATED candidates**, a distribution that is not the pool: Legacy beats its matched control at `member` (−0.368 [−0.554, −0.196], 11W/1L, above its own MDE), is under-MDE at `medoid`, and NOT MEASURED at `AVERAGED`. **AMBER is the mirror image — positive on all three readouts**, fold CIs excluding zero on the worse side at medoid (+0.335) and averaged (+0.368). Third independent instrument. | ESTABLISHED | §2.4 |
| A18 | **The cross-readout cell — the coordinator's priority — is a NEGATIVE.** Best off-diagonal `Distance→Dist+Amb` (3.922) vs best diagonal `Distance→Distance` (3.845) = **+0.077 [−0.187, +0.348], 5W/7L**: the off-diagonal is *worse*, and NOT MEASURED. Every training Hamiltonian read out with Distance lands 4.03–4.20 while Distance→Distance alone reaches 3.845. **Does not contradict Workstream B's −0.697 readout effect** — B measures a readout swap against a *bad* readout, this measures it against the *best* one. | ESTABLISHED (negative) | §2.5 |
| A19 | **The VQE optimises every one of the seven objectives** (`d_cvar` negative on all seven, −0.149 to −0.412) **and it does not convert**: ρ(objective gained, RMSD gained) is **+0.035** on Distance and **−0.476** on Dist+Amb. Sixth independent instrument for *search saturates, discrimination binds*. | ESTABLISHED | §2.6 |
| A20 | **RETRACTED AND RESTATED — my audit table was unmatched across the seed dimension.** It compared seed-0 audit arms to the declared form's 4-seed mean. Matched at seed 0 the declared form wins **3 of 4** hybrids, not 2 of 4 (`robustz` does *not* win Dist+Amb: 4.080 vs 4.086). | RETRACTED (the 2-of-4 count) | §2.7 |
| A20a | **The declared normalisation is recorded as WRONG on `Dist+Leg+Amb`**, where `rank` wins by 0.139 Å — my own pre-registered rule, applied rather than escaped by invoking low power. It is **not swapped**: the declared form stands elsewhere and no alternative is promoted, since 0.139 Å is inside the 0.51–1.06 Å seed noise. | recorded as declared-WRONG on one cell / NOT MEASURED (the difference) | §2.7 |
| A20b | **The normalisation fork is not load-bearing, but the margin is thinner than the review inferred from my broken table.** Best-per-cell normalisation reaches **3.866** against Distance-only at **3.807** — a **+0.059 Å** margin, not ~0.24 Å, and itself deep inside seed noise. Distance is **not shown** to beat the best-normalised hybrid; it is merely not beaten by it. §2.1's verdict rests on the declared comparison (0.27–0.94 Å), not on this reconstruction. | NOT MEASURED | §2.7 |
| A20c | **SEED COUNT was a fork I failed to enumerate**, and it is the one that most constrains the matrix — disclosed in §2.3 but absent from my fork list until the review. It is the reason every row-to-row magnitude stays NOT MEASURED. | recorded as process | fork table, §2.3 |
| A21 | **The identity `tail_min ≡ argmin` holds on VQE-generated sets too** — `max |diag − argmin(train)| = 0.000e+00` across the full 7×7 × 4-seed cross-readout. Verified now on the pool (n=42, both bases) and on generated candidates. | **EXACT** (verified three times) | §2.5 |
