# SPRINT 19 — AGENT C (LEGACY / AMBER / PHYSICS) — PRE-REGISTRATION

**Written before any Sprint-19 physics number existed.** The only numbers that existed when this
was written are the two verification gates (`GC0`, `GC1`, §0), which are identity checks and carry
no scientific direction, and the Sprint-18 record.

Per `s19/BRIEF.md` §14: *do not edit a pre-registration after seeing results.* If something below
is mis-specified I will say so, keep the text, and record the untested regime as OPEN.

My honest prior on each hypothesis is stated with it, before the run.

---

## 0. GATES — passed before any scientific number is quoted

`python -m s19.agentC_lib` → `s19/results/agentC_gates.json`

| gate | what it checks | result |
|---|---|---|
| **GC0** | `agentC_lib.common_frame_members` **is** the deployed averaging operator (`s15.phys_repl.averaged_backbone_from`), its mean **is** the metric's argument, and the Krogh–Vedelsby identity holds | **PASS** — mean delta **0.00e+00 Å**, readout-vs-metric 8.1e-15 Å, KV residual **3.4e-14 Å²**, min FRAME² +5.2e-04 Å² (non-negative, as the theorem requires) |
| **GC1** | my Legacy scores **are** the Sprint-18 artefact's Legacy scores, candidate by candidate, and my ORACLE labels are its ORACLE labels | **PASS** — max \|Δ Legacy term\| **0.00e+00**, max \|Δ d\| **0.00e+00 Å**, 8 targets × 75 candidates × 11 terms |

**A defect GC0 caught in its own first draft, recorded rather than repaired silently.** My first
(KV3) used the operator law's `d` — the free-superposition RMSD of the **raw window CA trace** —
as the member error, and FRAME² came out **negative (−0.059 Å²)**, which is impossible for a
quantity defined as a minimum-over-rigid-motions gap. The cause: the operator averages the
**ideal-geometry rebuild** of each window's (φ, ψ), not the window. The module now carries both
(`d_reb` inside the identity, `d_win` beside it) and never one in the other's place. This is
exactly the error class §7 of the brief warns about (a prior sprint broke the identity by 18% the
same way), caught by the gate before it reached a result.

---

## 1. WHAT I AM ATTACKING, AND WHAT I AM NOT

Sprint 18 closed, and I do **not** reopen: Legacy or AMBER as generic rankers; `leg_contact` as an
objective term; `leg_torsion` as a gate; fitted Legacy weights; a stronger AMBER minimiser.

I attack exactly three things:

1. **Q1 — the mechanism** behind the programme's hardest negative: *why* score-ordering damages an
   averaged set beyond what its mean and its best predict.
2. **Q2 — AMBER as a constrained terminal operator**: the Cα-accuracy / physical-validity Pareto
   frontier, reported as two axes, never one scalar.
3. **Q3 — the one Legacy role not yet refuted**: hard rejection of a specific class of impossible
   candidate.

---

## 2. THE INSTRUMENT, FIXED FOR EVERY ARM BELOW

* **Targets**: all 126 cluster-disjoint tuning targets, `s12.instrument.targets()`, pinned order.
  Sealed benchmark **not read**.
* **Candidate set**: the shipped top-75 windows, `s14.avgspace.top75_windows` — the identical 75
  Sprint 18's Phase 8 used, so every arm here is comparable to the record. No filter is handed a
  set built by its own criterion.
* **Primary readout**: full-chain Cα-RMSD of the **all-atom coordinate average** of the survivor
  set (`avg_rmsd`, the 3.048 Å object), frozen implementation `s12.instrument.ca_rmsd`. The
  ideal-geometry projection and the AMBER repair are reported as secondary readouts on the arms
  that reach them.
* **Statistics**: `s18.phys_lib.paired` unchanged — target-level paired bootstrap, **fold-clustered
  CI quoted as primary**, median and W/L beside every mean. n = 126 unless stated. MDE 0.084 Å.
* **Seeding**: `s15.seed.stable_rng` only.
* **Legacy**: `core.energy`, eleven components, `DEFAULT_WEIGHTS`, never fitted.
* **AMBER**: genuine ff14SB/GBn2 through `core.amber`. Convergence gate **declared here, before
  use**: finite final energy ≤ `CONVERGE_MAX_KCAL` = 1000 kcal/mol with the restraint switched
  off, reported with its exclusion count, every gated arm compared to **its own gated input**.
* **ORACLE**: the native enters evaluation and labelled ORACLE diagnostics only. `D` (diversity)
  and every gate below are native-free; `readout`, `E_mem`, `FRAME` are ORACLE diagnostics.

---

## 3. Q1 — THE AMBIGUITY DECOMPOSITION OF GATE DAMAGE

### 3.1 The identity (EXACT — a theorem, not a finding)

With `Y_i` the survivors in the operator's **own common frame**, `Ybar` their mean, `T` the native
superposed onto `Ybar`, and `||·||²` the per-residue mean square:

    readout²  =  E_mem²  −  D²                                                          (KV)
    readout²  =  <d_reb²>  +  FRAME²  −  D²                                             (KV3)

`E_mem` = common-frame member error, `D` = ambiguity/diversity, `<d_reb²>` = mean squared
free-superposition member RMSD, `FRAME²` = `E_mem² − <d_reb²>` ≥ 0. Asserted to 1e-9 Å² on every
call. **This is algebra. Nothing in §3.1 is a result.**

### 3.2 The hypothesis under test (the open question Sprint 18 left)

> **H-C1.** A physics score gate selects survivors whose errors are *correlated* — it collapses the
> ensemble's ambiguity `D` — and that is why the output degrades even when the set mean improves.
> A matched-random gate leaves the errors decorrelated, which is why the operator law
> `d_out = 1.16·d_set_mean + 0.04·d_set_best` holds for it (miss +0.006 [−0.003, +0.014]) and
> fails for every physics gate (Legacy +0.094, `leg_torsion` +0.079).

Only `<d²>` (and the best) is visible to that law. `FRAME²` and `D²` are invisible to it. H-C1
says the missing angstroms live in `D²`.

### 3.3 Arms

All at **m = 37 survivors of K = 75** (f = 0.50, the Sprint-18 primary), same 75 for every arm.

| arm | rule | role |
|---|---|---|
| `none` | no gate, m = 75 | the incumbent |
| `legacy` | 37 lowest by the eleven-term total at `DEFAULT_WEIGHTS` | score gate |
| `leg_torsion`, `leg_contact`, `leg_steric` | 37 lowest by that component | score gates |
| `amber_sp` | 37 lowest by genuine ff14SB/GBn2 single point | score gate |
| `helix` | 37 closest to a constant ideal α-helix | **ZERO-INFORMATION control** (§8 of the brief) |
| `rand` | 37 at random, **R = 30 `stable_rng` draws**, averaged | **MATCHED-RANDOM control** |

### 3.4 PRIMARY PRE-REGISTERED OUTCOME

For each score gate `g`, per target, against the **mean over the 30 random draws** at the same m:

    Δreadout²(g) = readout²(g) − readout²(rand)
    ΔE_mem²(g)   = E_mem²(g)   − E_mem²(rand)
    ΔD²(g)       = D²(g)       − D²(rand)          with  Δreadout² ≡ ΔE_mem² − ΔD²  (exact)

**The primary outcome is the Legacy arm's diversity channel `−ΔD²`**, at n = 126, with a
fold-clustered paired CI, and its **share** of the total damage,
`S_D = mean(−ΔD²) / mean(Δreadout²)`.

**H-C1 is SUPPORTED iff, for the `legacy` arm:**

1. `mean(Δreadout²) > 0` with a fold-clustered CI excluding zero (the damage reproduces), **and**
2. `mean(−ΔD²) > 0` with a fold-clustered CI excluding zero (the gate really does collapse
   ambiguity relative to random), **and**
3. `S_D ≥ 0.50` (the diversity channel carries at least half the damage).

**FALSIFIER F-C1 — declared now.** If (2) fails (CI spans zero) **or** `S_D < 0.50`, then
**H-C1 is REFUTED as the dominant mechanism** and the damage is attributed to whichever of
`ΔE_mem²` / `ΔFRAME²` carries it. I will report the attribution either way and will **not**
re-cut the arms, the m, or the readout to rescue it.

**Secondary, reported whatever happens**: the same three channels for every other gate; the full
(KV3) attribution `Δ<d_reb²> + ΔFRAME² − ΔD²`; the mean pairwise error-vector cosine and the
shared fraction `readout²/E_mem²`; the operator-law miss recomputed on my own arms so it is
verified rather than imported.

**My honest prior**: 0.6 that (1) and (2) hold, but only **0.35** that `S_D ≥ 0.50`. I expect the
`FRAME²` channel — the ensemble's internal rigid misalignment, which no one has measured — to take
a large share, because a score gate that keeps structurally similar candidates should shrink
`FRAME²` and `D²` together, and those two enter (KV3) with **opposite signs**.

### 3.5 The design-rule test that follows from H-C1 (pre-registered as a conditional)

If H-C1 predicts anything operational it is: **a gate must preserve error decorrelation.** Two
native-free, diversity-preserving gates of **exactly the same count** are therefore run in the same
module, whatever §3.4 returns:

| arm | rule |
|---|---|
| `legacy_spread` | reject on Legacy to 2m survivors, then take m by greedy farthest-point on CA-RMSD |
| `legacy_clust` | partition the 75 into m CA-RMSD clusters, keep each cluster's best-Legacy member |

**Pre-registered secondary outcome H-C2**: `legacy_clust` (and/or `legacy_spread`) beats the plain
`legacy` gate at the same m, with a fold-clustered CI excluding zero. **And the harder bar the
brief demands**: it must also not lose to `rand` at the same m. A diversity-preserving gate that
merely ties `rand` is a **mechanism confirmation, not a deployable gain**, and will be reported as
such.

**My honest prior**: 0.55 that it beats plain `legacy`; **0.2** that it beats `rand`.

---

## 4. Q2 — AMBER AS A CONSTRAINED TERMINAL OPERATOR: THE PARETO FRONTIER

**This is not an attempt to make AMBER a stronger minimiser.** The standing role — stereochemical
repair at an accuracy cost of +0.133 Å [+0.112, +0.165] at k = 30 — is not under attack. The
question is whether a *differently constrained* AMBER buys the same repair for less displacement.

### 4.1 Input, identical for every arm

The all-atom coordinate average of the **ungated** top-75 (`none`), i.e. the 3.048 Å object. Every
arm therefore has literally the same gated input, and each arm is compared to **that input** as
well as to the incumbent k = 30 arm.

### 4.2 THE TWO AXES — never collapsed into one number

* **Accuracy axis**: `ΔCα-RMSD` against the arm's own input (and the absolute Cα-RMSD).
* **Validity axis**, all from `s16.energy_lib.panel` on the same emitted structure, reported as a
  vector, never averaged into a scalar: `n_clash_2A`, `n_clash_2p6A`, `min_heavy`,
  `bond_strain`, `angle_strain`, `rama_favoured`, `rama_outlier`, `cis_frac`, `chirality_L_frac`.

A Pareto point is **useful** iff it is not dominated on both axes by k = 30, i.e. its clash and
strain statistics are no worse than k = 30's within noise **and** its ΔCα-RMSD is smaller.

### 4.3 Arms (declared in full, before the run)

| family | arm | protocol |
|---|---|---|
| restraint strength | `k1 k3 k10 k30 k100 k300 k1000` | restrained set {N, CA, C} (incumbent), `steps=0`, `tolerance=1.0` |
| restrained atom set | `caonly@k30`, `caonly@k100`, `caonly@k300` | **Cα pinned, N/C/O/CB and side chains free** |
| restrained atom set | `heavybb@k30` | {N, CA, C, O, CB} restrained |
| staged | `stage_300_30` | k = 300 to convergence, then k = 30 from that structure |
| staged | `pullback_30_1000` | k = 30, then re-minimise with the restraint reference reset to the **original input** at k = 1000 |
| projection | `blend@α` for α ∈ {0.25, 0.5, 0.75} | linear step-back of the k = 30 output toward its input — **free**, no extra AMBER; a validity-preserving projection *candidate*, and expected to fail the validity axis |
| adaptive | `adaptive_val` | **native-free**: per target, the *softest* k in the ladder whose output has `n_clash_2p6A == 0` and `bond_strain ≤` the k = 30 arm's; falls back to k = 300 |
| null | `frame` | the input rigidly rotated and translated before k = 30 — **exactly zero by construction**; reported with its **MAXIMUM** over targets, per §7 of the brief, never its mean |

`caonly` and `heavybb` change **only the positional-restraint atom list**, which is not part of
ff14SB; the force field, the solvent model, the topology and the minimiser are untouched. The
builder is constructed with a patched `core.amber.RESTRAINED_BACKBONE` and kept in a private
cache, because `core.amber.builder_for`'s cache key does not include the restraint set.

### 4.4 PRIMARY PRE-REGISTERED OUTCOME

**A Pareto point exists iff some arm has `n_clash_2p6A` and `bond_strain` statistically
indistinguishable from (or better than) `k30`'s, while its ΔCα-RMSD is smaller than `k30`'s by
more than the 0.084 Å MDE, with a fold-clustered CI excluding zero, at n = 126, gated.**

**FALSIFIER F-C2**: if no arm clears that bar, the answer is *"the +0.133 Å is the price of the
repair and the frontier has no better point"* — a clean negative, reported as the result.

**My honest prior**: 0.45 that `caonly` clears it (Cα is the axis the metric scores, and freeing
N/C/O is where most of the ideal-geometry strain lives), 0.15 for anything else, 0.35 that nothing
does.

---

## 5. Q3 — CAN LEGACY REJECT A CLASS OF IMPOSSIBLE CANDIDATE?

### 5.1 The class

Not "bad candidates" — **impossible** ones: hard steric violation. The detector is Legacy's own
`steric` component (never fitted, `DEFAULT_WEIGHTS`) and, as an independent physical criterion, the
candidate's minimum heavy-atom distance from `s16.energy_lib.panel`.

### 5.2 Arms — small rejections, unlike the m = 37 halving

Reject `r ∈ {2, 5, 10, 19}` of the 75 by:

| arm | rule |
|---|---|
| `steric_r` | the `r` worst by Legacy `steric` |
| `minheavy_r` | the `r` worst by minimum heavy-atom distance |
| `amber_r` | the `r` worst by genuine AMBER single point |
| `rand_r` | **MATCHED-RANDOM: `r` at random, 30 draws** |
| `randdiv_r` | **DIVERSITY-PRESERVING RANDOM: `r` at random, drawn to preserve the ensemble's `D`** — the extra control the brief demands given Q1 |

### 5.3 PRIMARY PRE-REGISTERED OUTCOME

**Legacy's steric rejection is a live role iff `steric_r` beats BOTH `rand_r` AND `randdiv_r`, at
some declared `r`, with a fold-clustered CI excluding zero at n = 126.** The primary `r` is
**r = 10**; the others are shape and are reported as shape.

**FALSIFIER F-C3**: if it does not beat both at r = 10, **this role is CLOSED** and I say so. No
sweep over `r` will be used to find a winner post hoc; r = 10 is the pre-registered primary and any
other `r` that looks better will be labelled a hyperparameter chosen on the tuning instrument.

**My honest prior**: 0.25. The Sprint-18 result is that the ordering is what costs, and a small
rejection is a small dose of the same ordering. The one thing that could differ is that the *tail*
of `steric` may be genuinely non-physical rather than merely low-ranked — which is precisely the
distinction this test exists to draw.

---

## 6. WHAT WOULD MAKE ME REPORT A NEGATIVE, LOUDLY

* F-C1 fires → the correlated-survivor hypothesis, which is the most interesting thing in my lane,
  is wrong, and I lead the findings with that.
* F-C2 fires → AMBER's displacement tax has no better Pareto point; the incumbent stands.
* F-C3 fires → the last unrefuted Legacy role closes and Legacy is a demonstrated non-contributor
  in every architectural role this programme has tested.

All three firing is a completely acceptable sprint outcome and I will not soften it.

---

## 7. ARTEFACTS

| file | contents |
|---|---|
| `s19/results/agentC_gates.json` | GC0, GC1 |
| `s19/results/agentC_kv.json` | Q1 + §3.5, per target, per arm, all three channels |
| `s19/results/agentC_pareto.json` | Q2, per target, per arm, both axes, convergence flags |
| `s19/results/agentC_reject.json` | Q3, per target, per arm, per `r` |
| `s19/agentC_FINDINGS.md` | the report |

Every artefact carries `complete`, a row count and a config hash, written through
`s18.phys_lib.write`; `read_complete` refuses to read a partial as a result. Nothing overwrites a
Sprint-18 file. Smokes are written to `_SMOKE_*` names and are never quoted as results.
