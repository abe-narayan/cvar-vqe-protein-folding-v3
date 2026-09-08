# SPRINT 15 — QRESTRAINT findings

**Does Sprint 14's decisive quantum negative survive a change of OBJECTIVE CLASS, from
selective to generative?**

Tiering: **DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **HYPOTHESIS** / **REFUTED**.

Module: `s15/qrestraint.py`. Artefacts: `s15/results/qrestraint.json` (raw),
`s15/results/qrestraint_analysed.json` (all tables), `s15/results/qrestraint_alpha.json`,
`s15/results/qrestraint_mode.json`, `s15/results/qrestraint_untrained_enrich.json`,
`s12/results/s15_qrestraint.json`.

**Scale.** 19 fully enumerated targets (nine n=9 at 4^9, ten n=10 at 4^10) =
**1.31e7 exactly-labelled structures**; 6 objectives; 3 seeds; 3 budgets; 7 arms; matched
hard budget enforced by `Counter` and verified cell by cell.

**Instrument, at the end of the workstream.** `python -m s12.instrument` returns
`shipped 3.4540004952559396, pool_best 1.7108244199364904, top75_best 2.3061526409453816,
synthesis_fit 3.2040761603809194, n_zero_recall 18` — unchanged. Per the Phase 0 audit this is
a **cache read** that asserts only three of those five constants, so it is reported as
"nothing moved", not as a reproduction of the pipeline.

---

## THE FOUR ANSWERS

**(a) Does CVaR-VQE beat best-of-N from the untrained circuit on a generative restraint
objective, at matched budget, with a CI excluding zero?** — **Yes, but the finding is void as
a statement about objective class**: at 8192 evaluations VQE beats the untrained control by
−0.202 Å [−0.400, −0.035] on `E_ml_pred` (13/19) and −0.184 [−0.368, −0.010] on `E_ls_pool`
(14/19), **and by −0.192 Å [−0.320, −0.076] on `S14_disto_bayes` (15/19), which is the
SELECTIVE control** — the objective with the worst in-tail ordering in the whole study. The
effect is the same size on every objective regardless of ordering skill, so it is a property
of the sampler and the budget, not of the objective.

**(b) Does the answer change with budget?** — **Yes, and it is non-monotone**: **1 of 6**
objectives has a CI excluding zero at 2048 (means −0.018 to −0.162), **3 of 6** at 8192 (means
−0.147 to −0.205), and **0 of 6** at 32768 (means −0.094 to +0.041). The gap **widens then
narrows**, peaking at 3.1% coverage of an n=9 space and gone by 12.5% — and §6.4 gives the
mechanism: more budget means more concentration, which costs the distinct samples a best-of-N
readout is paid in.

**(c) Is the generative objective's in-tail ordering skill materially better than the
selective energies'?** — **Yes on the rho statistic, marginally on Sprint 14's own statistic,
and NOT in the direction Sprint 14 measured**: in-tail rho is +0.127 (z = 13.1) for
`E_combined` against +0.048 for `S14_disto_bayes`, −0.042 for Legacy and −0.030 for AMBER; but
Sprint 14's `tail − bulk` degradation is **negative for every objective including all four
generative ones** (−0.134 to −0.237), so "every objective gets worse inside its own tail"
reproduces exactly and is **not** a pathology of the selective class.

**(d) If VQE still loses despite better in-tail ordering, what is the mechanism instead?** —
VQE does not lose; it modestly wins against the untrained control and **ties or loses against
classical search at the same budget** — **0 of 16 cells where VQE beats greedy with a CI
excluding zero and 3 of 16 where it loses; 0 of 16 against annealing and 3 of 16 where it
loses**. The mechanism is that **the VQE-over-control gain is
objective-independent** — six objectives spanning in-tail rho from −0.075 to +0.379 and global
rho from 0.316 to 0.797 all give the same −0.145 to −0.205 Å at 8192 — so it is not
discrimination at all: it is a moderately concentrated sampler beating a fixed one at finding
low-E points, which greedy does better. Underneath, **the selection gap (1.47–2.14 Å) is
identical for VQE and the control and does not move with budget**, and collapses to 0.29–0.44 Å
under ORACLE restraints. The binding quantity is restraint error, not the sampler.

**The one-sentence version.** Sprint 14's negative does not reproduce in its strong form —
running VQE is not *worse* than not running it on a generative objective — but its *mechanism*
reproduces intact, and the replacement finding is stronger and more awkward: **the choice of
sampler is worth ~0.2 Å at one budget and nothing at the others, against a selection gap of
~1.9 Å that no sampler and no objective in this study moves at all.**

**The finding I did not expect and would have missed without the defect arm.** *Less
optimisation is better*, confirmed on three causally independent axes — raising the budget,
raising alpha, and **correcting the CVaR gradient** each increase concentration, each cost
distinct samples, and each **reduce** the gain. The recorded `baseline="tail"` gradient defect,
whose bias has a closed form, turns out to be a de-facto step-size reduction and is
**better** on returned RMSD (−0.266 vs −0.145 on the ORACLE objective at α = 0.25). §6.3.

**The one unambiguously positive quantum-side result, and why it does not show up above.** The
trained variational state multiplies the probability of drawing a sub-2 Å structure by
**5.2–6.7× relative to its own untrained initialisation** (18.6× under ORACLE restraints), with
a clean internal control — the weakest objective gives 0.98×. That is a real, exactly-enumerated
distribution property. It converts into ≤ 0.2 Å through the arm table because the readout is an
argmin over samples, **and an argmin does not care about probability mass**. Realising it needs
an operator that consumes the ensemble — `coord_FINDINGS` Family C. §7.

---

## 0. WHAT WAS BROKEN, WHAT WAS NOT, AND WHAT IS NOW VERIFIED

**DEMONSTRATED.** The module was handed over untested. It ran. Recording the negative too:
**I found no arithmetic bug in the handed-over code** — the defects were omissions of controls,
not errors.

### 0.1 `all_ca` reproduces `Enum.rmsd`, and the residual is the LABEL's dtype, not ours

`_check_consistency` passes at `tol = 1e-6` on all 19 targets. Maxima over the whole set:

| quantity | max over 19 targets |
|---|---|
| &#124;rebuilt − `Enum.rmsd`&#124;, CA stored **float32** (as written) | **7.178e-07** |
| &#124;rebuilt − `Enum.rmsd`&#124;, CA rebuilt entirely in **float64** | **7.119e-07** |

The brief was explicit that the tolerance must not be weakened to make it pass, so I traced
the residual instead of adjusting the threshold. The float64 rebuild — which shares no storage
rounding with our array — disagrees by essentially the same amount. **The cached `rmsd` label
is itself stored as float32** (`s14/obj_enum.py` writes `np.empty(B, np.float32)`; confirmed
`d["rmsd"].dtype == float32`), whose spacing at 3 Å is 2.4e-07 and at 5 Å is 4.8e-07. So ~99%
of the residual is the label's own dtype and the rest is our float32 CA storage. Both numbers
are returned by `_check_consistency` on every target and written into the artefact.

The consequence that matters: **the objective and the RMSD index the same structure.** The
Sprint 12 defect that cost a whole table is not present.

This is a real check, not a tautology: `Enum.rmsd` is built by `s13.qarch_lib.Space.rmsd` as
`I.kabsch_rmsd_batch(I.build_ca(PHI[rows, S], PSI[rows, S]), nat_ca)`, so agreement confirms
the big-endian base-k index decoding in `states()` as well as the geometry.

### 0.2 The corrected `LogPTable` / `SumLogP` are in use, and the shapes are right

`rowsum` on a `(B, npairs)` block returns `(B,)`; `SumLogP.rowsum` adds the two tables'
rowsums. All five objectives are finite over the full 262,144 / 1,048,576 configuration spaces
with no NaN and no infinity, on all 19 targets. `s15/distml.py` on disk is the `searchsorted`
version, so nothing here inherits K4-RETRACTED.

One check the handover lacked, now added and passing on every target:

```
assert np.array_equal(i, dg["i"]) and np.array_equal(j, dg["j"])
```

`I.pair_index(n)` and the distogram's own pair order must agree or `dhat`, `sd`, `prob` and the
computed distances would be permuted relative to one another.

### 0.3 `V.summarise` and `R.run` are called with what they expect

Verified by running them. `R.run` returns every key the caller reads whenever `rmsd` is passed.
The handed-over `run_target` took a positional `CA` argument it never used; removed.

### 0.4 Memory

Never a constraint, contrary to the handover note. Largest array is the n=10 CA block at
**126 MB** float32 with a **151 MB** pair-distance block beside it, both freed before any arm
runs. `V.wait_for_memory(1.6)` gated every target and never fired.

### 0.5 THE CONTROL IS THE VQE'S OWN INITIALISATION — bit-identical, verified

The load-bearing property of the whole experiment. `untrained_best_of_n` builds `theta0` with
`R.init_theta(an, np.random.default_rng(seed), 0.8, "random")`; `R.run` builds its start the
same way from `np.random.default_rng(seed)`. Measured at seeds 0, 1, 2:
**`theta0` identical, max difference exactly 0.0.**

The control is therefore not "some untrained circuit" — it is **the same circuit at the same
parameters the VQE starts from**, and the arms differ in exactly one thing: whether the
parameters are subsequently updated. (The sampling RNG streams differ, since `R.run`'s
generator is advanced by `init_theta`; both draw from the same `p_theta0`, so this is noise,
not bias.)

### 0.6 STATUS OF THE THREE RECORDED CVaR DEFECTS

| Sprint 14 defect | active in the main table? | how handled |
|---|---|---|
| **1. `baseline="tail"` gradient bias** | **NO** | `R.run` defaults to `baseline="const"`, the corrected control variate (cos +1.000000 against the exact-expectation reference). `alpha_sweep` additionally runs the **defective** estimator as its own arm at matched everything, so the defect is measured on this objective class rather than asserted away. |
| **2. sampled-CVaR upward bias at non-integer `alpha·N`** | **NO** | shots = 512, alpha = 0.25 → `alpha·N = 128`, an integer. Active only at the sweep's 0.02 / 0.05 / 0.1, flagged per row. |
| **3. `dCVaR/dp` identically zero iff `p(argmin E) ≥ alpha`** | measured | every VQE run records `grad_norm_first` and a `dead_start` flag. |

**So the negative-vs-positive verdict below is not attributable to a defective gradient
estimator.**

---

## 1. WHAT I CHANGED IN THE HANDED-OVER DESIGN

1. **Nineteen targets, not nine.** `s14/cache/obj_enum_<PDB>.npz` holds ten more n=10 targets
   on the identical schema (verified: same key set, index algebra
   `snap_states @ k^(n-1-i) == snap_index`). This mattered: see §5, where n = 9 and n = 1
   both produced conclusions that n = 19 overturned.
2. **A budget ladder**: 2048 / 8192 / 32768, a 16× span — which is 0.78 / 3.12 / 12.50% of an
   n=9 space and 0.20 / 0.78 / 3.12% of an n=10 space. The ladder is only interpretable in
   those units and both are printed.
3. **The comparators are recomputed here, not quoted.** Sprint 14's in-tail numbers are
   *tail-restricted pairwise accuracy*; the handover's `ordering_skill` is *Spearman rho over
   the tail*. Comparing one to the other compares statistics, not objectives. Both are
   computed, on the same targets and tails, for the generative objectives **and** for
   `legacy`, `prior`, `amber_total`.
4. **`S14_disto_bayes` added — the sharpest control available, and it changed the answer.**
   Sprint 14's `w = 1` objective (in-tail accuracy 0.390, its worst anti-ranking) is
   `s14/hamil.py::e_disto` = `I.shipped_score`, the distogram's **Bayes-risk** score. It reads
   *the same distogram* as `E_ls_pred` and `E_ml_pred` and differs only in how it **consumes**
   it: an expected-loss lookup scoring a candidate for nativeness (SELECTIVE) versus a
   residual against predicted distances (GENERATIVE). Same information, same targets, same
   tails, same optimiser, same budget. **Without it I would have reported a positive quantum
   result caused by objective class.**
5. **A random-tail null on every in-tail statistic**, per objective per target.
6. **Concentration measured, not presumed** — every VQE run records the mean global
   E-percentile of its first and last quarter of samples.
7. **`wait_for_cpu`**, since `wait_for_memory` gates RAM and nothing gated CPU.

AMBER rows obey the Phase 0 binding rule `amber_kind == 0 AND amber_idx != snap_index`.
Surviving n per target: **1190–1195 of ~2,900** on the n=9 targets and **698–700 of ~1,780** on
the n=10 targets — i.e. **39–41%**, which independently confirms the Phase 0 correction that
40% is the unbiased share (the brief's original §2 one-liner had it transposed). The 0.1% tail
of a ~1,100-row subsample is ~1 structure, below `tail_accuracy`'s own `min_tail = 30` guard,
so AMBER's q = 0.01 and q = 0.001 cells return **NaN by design**. They are left visible.

---

## 1b. INDEPENDENT VERIFICATION — Sprint 14's flagship table reproduces exactly

**DEMONSTRATED.** The brief requires an independent check for every headline claim, and the
claim here is a *comparison against Sprint 14*, so what had to be verified is that this
instrument reproduces Sprint 14's own numbers.

`S14_disto_bayes` reaches Sprint 14's `w = 1` objective by a different code path —
`I.shipped_score` on distances from `all_ca`-built traces over the **full enumerated space**,
where `Terms.e_disto` used `Space`-built traces over a **4,000-configuration sample** — and is
scored with `s14.ener_lib.tail_accuracy` at Sprint 14's own settings. On Sprint 14's own nine
targets:

| tail-restricted pairwise accuracy, gap > 0.25 Å | q=1.0 | q=0.1 | q=0.01 | q=0.001 | tail − bulk | certified argmin |
|---|---|---|---|---|---|---|
| **Sprint 14 recorded** (`s14/results/vqe_tailacc.json`, w=1.0) | 0.655 | 0.549 | 0.509 | **0.390** | −0.266 | 3.237 |
| **this module** (`S14_disto_bayes`) | 0.655 | 0.549 | 0.509 | **0.391** | −0.264 | 3.237 |

Identical to three decimals on all four fractions; the certified argmin agrees exactly. Two
consequences. First, every comparison below is like-for-like rather than across instruments.
Second, **Sprint 14's `w=1` number is confirmed by re-derivation** — it was not a sampling
artefact, since this is an exhaustive enumeration of the same space and theirs was a 4,000-point
sample.

---

## 2. RESULTS — all nineteen targets

### 2.1 ORDERING SKILL, Spearman rho, against a matched random-tail null

**DEMONSTRATED.** Mean over 19 targets. `z` = rho / (sd of the same statistic over random
tails of the same size).

| objective | class | rho global | rho best-1% | z | rho best-0.1% | z | tail gap 1% | argmin RMSD |
|---|---|---|---|---|---|---|---|---|
| `E_ls_pred` | GENERATIVE | 0.449 | +0.076 | 7.2 | +0.030 | 0.9 | 1.234 | 3.123 |
| `E_ml_pred` | GENERATIVE | 0.513 | +0.104 | 10.6 | −0.026 | −0.8 | 1.274 | 2.855 |
| `E_ls_pool` | GENERATIVE | 0.316 | +0.044 | 3.5 | +0.001 | 0.0 | 1.352 | 3.269 |
| **`E_combined`** | GENERATIVE | 0.523 | **+0.127** | **13.1** | +0.050 | 1.6 | 1.287 | **2.767** |
| `E_ORACLE_true` | **ORACLE** | 0.797 | +0.379 | 68.8 | +0.294 | 16.3 | 1.208 | 1.439 |
| `S14_disto_bayes` | SELECTIVE | 0.485 | +0.048 | 4.7 | **−0.075** | **−2.4** | 1.256 | 3.077 |
| `legacy` | SELECTIVE | 0.105 | **−0.042** | **−2.9** | −0.041 | −0.9 | 1.676 | 4.085 |
| `prior` | SELECTIVE | 0.005 | −0.009 | −0.6 | +0.023 | 0.5 | 2.515 | 3.463 |
| `amber_total` | SELECTIVE | 0.148 | −0.030 | −0.2 | −0.030 | −0.2 | 1.484 | 3.931 |

*(space best RMSD 1.030 Å for the non-AMBER rows, 1.557 Å on AMBER's masked subsample.)*

**The generative objectives are the only ones with positive, significant in-tail rho.** At the
1% tail all four are positive and three clear their null by z ≥ 3.5, while Legacy is
*negative* at z = −2.9, `prior` and AMBER are at zero, and the distogram consumed as a Bayes
risk is only +0.048 at 1% and **−0.075 at 0.1%**. At the 0.1% tail only the ORACLE is
convincingly positive.

### 2.2 SPRINT 14's OWN STATISTIC — the LEVEL differs, the DIRECTION does not

**DEMONSTRATED.** Tail-restricted pairwise accuracy at gap > 0.25 Å, each read against **its
own** matched random-tail null (never against 0.500). The null column equals the q=1.0 column
to within 0.001 throughout, which is the correct behaviour and is printed in the artefact.

| objective | class | q=1.0 | q=0.1 | q=0.01 | q=0.001 | **tail − bulk** |
|---|---|---|---|---|---|---|
| `E_ls_pred` | GENERATIVE | 0.692 | 0.558 | 0.533 | 0.507 | −0.185 |
| `E_ml_pred` | GENERATIVE | 0.721 | 0.574 | 0.546 | 0.484 | −0.237 |
| `E_ls_pool` | GENERATIVE | 0.636 | 0.537 | 0.520 | 0.502 | −0.134 |
| **`E_combined`** | GENERATIVE | 0.727 | 0.580 | **0.558** | **0.532** | −0.196 |
| `E_ORACLE_true` | **ORACLE** | 0.867 | 0.726 | 0.678 | 0.643 | −0.224 |
| `S14_disto_bayes` | SELECTIVE | 0.708 | 0.559 | 0.521 | **0.466** | **−0.243** |
| `legacy` | SELECTIVE | 0.545 | 0.544 | 0.484 | 0.475 | −0.070 |
| `prior` | SELECTIVE | 0.504 | 0.495 | 0.496 | 0.512 | +0.008 |
| `amber_total` | SELECTIVE | 0.557 | 0.486 | NaN | NaN | — |

Two statements that must be kept apart, because conflating them is exactly how §5's retracted
reads went wrong:

- **The DIRECTION reproduces.** `tail − bulk` is negative for every objective, all four
  generative ones included, and the ORACLE's is −0.224 — as negative as anything else.
  Sprint 14's *"every objective is at chance or worse inside its own tail"* is **not a
  pathology of the selective class**. It survives the change of objective class and it
  survives perfect restraints.
- **The LEVEL is better, but only modestly.** At q = 0.001 the generative objectives sit at
  0.484–0.532 against 0.466 for the same distogram consumed selectively and 0.475 for Legacy.
  `E_combined` 0.532 versus `S14_disto_bayes` 0.466 is **+0.066** on identical information.
  Real, and an order of magnitude smaller than the 1.9 Å selection gap it would have to close.

### 2.3 THE SELECTION DECOMPOSITION

**DEMONSTRATED.** `sel = pool_mean + FILTERING + ORDERING` at q = 0.01
(`s14.ener_lib.selection_decomposition`). Sprint 14: *"the ORDERING term's sign IS
optimise-harder-get-worse, measured directly."*

| objective | class | pool mean | tail mean | argmin | FILTERING | **ORDERING** |
|---|---|---|---|---|---|---|
| `E_ls_pred` | GENERATIVE | 4.003 | 3.210 | 3.123 | −0.793 | −0.087 |
| `E_ml_pred` | GENERATIVE | 4.003 | 3.032 | 2.855 | −0.971 | −0.177 |
| `E_ls_pool` | GENERATIVE | 4.003 | 3.544 | 3.269 | −0.459 | −0.275 |
| `E_combined` | GENERATIVE | 4.003 | 3.032 | 2.767 | −0.971 | −0.264 |
| `E_ORACLE_true` | **ORACLE** | 4.003 | 2.238 | 1.439 | −1.765 | −0.799 |
| `S14_disto_bayes` | SELECTIVE | 4.003 | 3.154 | 3.077 | −0.849 | **−0.077** |
| `legacy` | SELECTIVE | 4.003 | 3.845 | 4.085 | −0.158 | **+0.240** |
| `prior` | SELECTIVE | 4.003 | 4.117 | 3.463 | +0.114 | −0.654 |
| `amber_total` | SELECTIVE | 4.006 | 3.782 | 3.931 | −0.223 | **+0.149** |

**Legacy's ordering term is +0.240 here against Sprint 14's recorded +0.263** — an independent
replication of the single number Sprint 14 called *"'optimise harder, get worse', measured
directly"*. AMBER's is +0.149. Both physical energies still have the pathology.

> **CORRECTION, in place.** On the nine n=9 targets alone this table showed a perfectly clean
> split — every generative objective negative, **every** selective objective positive, with
> `S14_disto_bayes` at **+0.163**. At n = 19 `S14_disto_bayes` is **−0.077** and `prior` is
> **−0.654**, so *"selective objectives have a positive ordering term"* is **REFUTED as a
> general claim**. What survives is narrower and is what is written above: the two physical
> energies (Legacy, AMBER) have a positive ordering term; the distance-derived objectives do
> not, however they consume the distogram.

The generative objectives' real advantage in this table is **FILTERING**, not ordering:
−0.79 to −0.97 Å against Legacy's −0.16 and AMBER's −0.22. They are much better at putting
the tail in a good neighbourhood, which is a bulk property.

### 2.4 THE ARM TABLE — the pre-declared primary comparison

**DEMONSTRATED.** `vqe_cvar_a0.25` vs `untrained_bestofN`, paired over 19 targets, negative =
VQE better, `*` = CI excludes zero. The n=9 / n=10 split is printed because the two register
sizes are different instruments.

| objective | 2048 | 8192 | 32768 | n=9 / n=10 at 8192 |
|---|---|---|---|---|
| `E_ls_pred` | −0.018 [−0.174,+0.130] 7/19 | −0.147 [−0.362,+0.062] 12/19 | +0.041 [−0.124,+0.201] 6/19 | −0.215 / −0.085 |
| `E_ml_pred` | −0.062 [−0.245,+0.124] 12/19 | **\*−0.202 [−0.400,−0.035]** 13/19 | +0.005 [−0.108,+0.109] 9/19 | −0.211 / −0.193 |
| `E_ls_pool` | **\*−0.162 [−0.304,−0.022]** 12/19 | **\*−0.184 [−0.368,−0.010]** 14/19 | −0.036 [−0.156,+0.072] 9/19 | −0.099 / −0.261 |
| `E_combined` | −0.109 [−0.306,+0.074] 10/19 | −0.205 [−0.424,+0.003] 13/19 | −0.043 [−0.172,+0.094] 10/19 | −0.038 / −0.356 |
| `E_ORACLE_true` | −0.092 [−0.266,+0.052] 9/19 | −0.145 [−0.388,+0.050] 12/19 | −0.094 [−0.252,+0.052] 11/19 | +0.032 / −0.305 |
| **`S14_disto_bayes`** (SELECTIVE) | — | **\*−0.192 [−0.320,−0.076]** 15/19 | — | −0.133 / −0.244 |

**THE RESULT THAT VOIDS THE TIDY STORY.** The strongest cell in the table, by W/L and by CI
width, is the **selective** control — the objective with the *worst* in-tail ordering measured
anywhere in this study (0.466 at q=0.001; in-tail rho −0.075 at z = −2.4). An objective that
anti-ranks its own tail produces the same VQE-over-control gain as one that ranks it at
z = +13.1. **The gain cannot be attributed to the generative objective's in-tail ordering,
because it does not depend on in-tail ordering at all.**

Every concentration check is reported as one verdict: at 8192 the four significant cells all
read **PASS** against the simulated uniform-effect null (the effect is spread across targets,
not carried by a few); the non-significant cells read **NO POWER** (|mean/sd| < 0.3), which is
the honest statement, not PASS.

### 2.5 VQE vs CLASSICAL SEARCH — the control a quantum claim actually needs

**DEMONSTRATED.** Same matched budget. Negative = VQE better.

| objective | budget | vs greedy | vs anneal |
|---|---|---|---|
| `E_ls_pred` | 2048 / 8192 / 32768 | +0.147 / −0.012 / +0.041 | +0.079 / −0.013 / +0.041 |
| `E_ml_pred` | 2048 / 8192 / 32768 | −0.008 / −0.077 / +0.023 | −0.012 / −0.065 / +0.023 |
| `E_ls_pool` | 2048 / 8192 / 32768 | **+0.191\*** / **+0.152\*** / +0.073 | +0.154 / **+0.148\*** / +0.073 |
| `E_combined` | 2048 / 8192 / 32768 | +0.044 / +0.009 / −0.015 | +0.091 / +0.006 / −0.015 |
| `E_ORACLE_true` | 2048 / 8192 / 32768 | **+0.270\*** / +0.155 / +0.131 | **+0.264\*** / **+0.193\*** / +0.128 |
| `S14_disto_bayes` | 8192 | −0.011 | −0.063 |

> **0 of 16 cells where VQE beats greedy with a CI excluding zero; 3 of 16 where it LOSES
> (`E_ls_pool` at 2048 and 8192, `E_ORACLE_true` at 2048). Against annealing: 0 of 16 wins,
> 3 of 16 losses (`E_ls_pool` at 8192, `E_ORACLE_true` at 2048 and 8192).**

And on the **objective axis** the two arms are not close: at 8192 on `E_ls_pred`, greedy's mean
objective gap is **0.165** and the VQE's is **3.441**; at 32768 greedy and anneal reach
**0.000** — the certified global optimum — on most cells while the VQE is still at 2.1–4.2.
The two axes come apart exactly as Sprint 14 said. **There is no quantum advantage here, on
either axis.**

### 2.6 THE MECHANISM — question (d) answered

**DEMONSTRATED.** `vqe_cvar_a0.25` vs `untrained_bestofN`, 19 targets. `conc` = mean global
E-percentile of the last quarter of samples minus the first (negative = concentrated);
`dist_frac` = fraction of the budget spent on distinct configurations; `selgap` = returned
RMSD minus the best RMSD the arm actually sampled.

| budget | objective | conc | dist_frac V | dist_frac C | d_objgap | **d_rmsd** | selgap V | selgap C |
|---|---|---|---|---|---|---|---|---|
| 8192 | `E_ls_pred` | −0.194 | 0.508 | 0.610 | −9.364 | **−0.147** | 1.765 | 1.797 |
| 8192 | `E_ml_pred` | −0.205 | 0.508 | 0.610 | −4.743 | **−0.202** | 1.472 | 1.504 |
| 8192 | `E_ls_pool` | −0.234 | 0.495 | 0.610 | −2.065 | **−0.184** | 2.027 | 2.137 |
| 8192 | `E_combined` | −0.198 | 0.505 | 0.610 | −7.280 | **−0.205** | 1.481 | 1.512 |
| 8192 | `E_ORACLE_true` | −0.203 | 0.520 | 0.610 | −17.716 | **−0.145** | 0.413 | 0.292 |
| 8192 | `S14_disto_bayes` | −0.205 | 0.499 | 0.610 | −0.058 | **−0.192** | 1.721 | 1.793 |
| 2048 | (range over 6) | −0.20 to −0.24 | ~0.785 | 0.779 | — | −0.018 to −0.162 | 0.38–1.91 | 0.30–2.04 |
| 32768 | (range over 6) | −0.12 to −0.14 | ~0.228 | 0.418 | — | −0.094 to +0.041 | 0.44–2.02 | 0.37–2.08 |

Four findings, in order of importance:

**(i) The gain is OBJECTIVE-INDEPENDENT, and that is the mechanism.** At 8192 the six
objectives span in-tail rho from −0.075 to +0.379, global rho from 0.316 to 0.797, in-tail
accuracy from 0.466 to 0.643, and an objective-gap improvement from −0.058 to −17.7 in their
own units — and every one of them yields **−0.145 to −0.205 Å**. A quantity that is constant
across a 5-fold range of objective quality is not caused by objective quality. The VQE at
intermediate budget is simply a **moderately concentrated sampler** beating a fixed one at
finding low-E points; any objective whose global rho is positive converts that into slightly
lower RMSD. Greedy does the same thing far better on the objective axis and no better on the
structural one.

**(ii) The premise of the Sprint 14 mechanism holds — and is not what decides the outcome.**
`conc` is −0.12 to −0.24: the optimiser genuinely moves its distribution onto the low-E tail.
So tail-blindness *is* applicable in principle. It just does not govern the result, because
the result does not depend on the tail.

**(iii) Concentration costs diversity, and a best-of-N readout is paid for in independent
draws.** At 32768 the VQE covers **0.230** of its budget in distinct configurations against the
control's **0.418** — a 45% loss. This is exactly why the gain **vanishes at the top budget**:
the concentration advantage and the diversity penalty cross over. It is a loss channel wholly
independent of ordering skill, and it is the one that survives when tail-blindness does not.

**(iv) THE BINDING QUANTITY IS THE SELECTION GAP, and nothing here moves it.** Measured two
ways, which are related but not the same quantity and are therefore quoted separately:

- *Arm-level* (returned RMSD minus the best RMSD the arm actually sampled): **1.47–2.14 Å** on
  the predicted objectives, **the same for VQE and the control to within 0.11 Å** on every
  cell, and flat across the 16× budget range (1.37–2.04 at 2048, 1.47–2.03 at 8192, 1.51–2.02
  at 32768). Under **ORACLE** restraints it collapses to **0.29–0.44 Å**.
- *Objective-level* (`tail mean − tail best` at q = 0.01, the direct analogue of Sprint 14's
  1.93 Å): **1.23–1.35 Å** for the four generative objectives, 1.26 for the distogram Bayes
  risk, 1.68 for Legacy, 1.48 for AMBER, 2.52 for the prior, and **1.21 even for the ORACLE**.

Both say the same thing. The gap is a property of **restraint error and of the k = 4 space** —
not of the objective class, not of the optimiser, not of the sampler, and not of the budget.
Against ~1.5–1.9 Å of selection gap, the entire quantum-versus-classical question is worth
0.2 Å at one budget and nothing at the others.

### 2.7 THE BUDGET-TRAP LAW — REFUTED at n = 19, on this objective class

The brief's standing law: *"the budget trap is a property of BAD OBJECTIVES. Searching harder
makes structures worse on a bad objective, is neutral on a mediocre one, and HELPS
monotonically on a good one."* Tested per (objective, target) cell with x = in-tail accuracy at
q = 0.01 minus its own matched null, y = RMSD(32768) − RMSD(2048). The law predicts
rho(x, y) < 0.

| arm | n = 9 cells (45) | **n = 19 cells (95)** | ORACLE cells removed | mean y | frac(y<0) |
|---|---|---|---|---|---|
| `greedy` | **−0.491** | **−0.073** | −0.105 (n=76) | −0.029 Å | 0.22 |
| `vqe_cvar_a0.25` | **−0.360** | **−0.223** | −0.179 (n=76) | −0.107 Å | 0.58 |

> **REFUTED for `greedy` at n = 19, and my own n = 9 reading of it is retracted.** On the nine
> n=9 targets I measured rho = −0.491 and wrote that "the law holds on this objective class".
> On nineteen it is **−0.073** — indistinguishable from no relationship. The law is not
> confirmed on real native-free objectives by this experiment. It survives weakly for the VQE
> arm (−0.223, −0.179 without the ORACLE) and that is the most that can be said.

The confound is reported rather than hidden: `x` is not a clean quality axis, because a better
objective has a *narrower* tail (the brief's own trap). `E_ORACLE_true`, by far the best
objective here, has the most negative `x` in the table (−0.189). Both the pooled and the
ORACLE-removed correlations are given for that reason.

---

## 3. WHAT THIS MEANS FOR THE SPRINT 14 NEGATIVE

**Sprint 14's strong claim does NOT reproduce.** *"Running VQE is worse than not running it"*
— 0 wins in 12 at +0.65 to +1.32 Å — does not hold on a generative restraint objective. At
matched budget the VQE is at parity or modestly ahead of best-of-N from its own untrained
initialisation (0/6 objectives worse with a CI excluding zero at any budget; 4 cells better).
Anyone quoting the 0/12 as a general property of CVaR-VQE would be overreaching, and the
scope condition is now measured: **it was a property of Legacy and AMBER**, whose ordering
terms are +0.240 and +0.149 here and whose in-tail rho is −0.042 and −0.030.

**Sprint 14's MECHANISM reproduces intact, and generalises further than Sprint 14 claimed.**
Every objective still degrades inside its own tail (`tail − bulk` −0.134 to −0.243), including
all four generative ones **and the ORACLE**. The arm-level selection gap reproduces at
1.47–2.03 Å against the recorded 1.93. Legacy's ordering term reproduces at **+0.240** against
the recorded **+0.263**. The distogram Bayes-risk objective's in-tail accuracy reproduces at
**0.391** against the recorded **0.390**. Three independent replications of Sprint 14 numbers,
on a larger target set and by different code.

**And the replacement claim is stronger than either.** The experiment was designed to ask
whether objective class decides the quantum result. The answer is that **objective class does
not decide it, because nothing about the objective decides it**: a 5-fold range of measured
ordering skill produces an identical −0.145 to −0.205 Å. What decides the returned structure is
the selection gap, which is set by restraint error and which no sampler in this study moves.

This is a cleaner statement of the project's standing position than "concentration is wrong
when discrimination binds", because it survives the case where discrimination does *not* bind:
even with ORACLE restraints (in-tail rho +0.379, z = 68.8; in-tail accuracy 0.643) the VQE
still fails to beat greedy, still loses 0.29–0.44 Å to its own selection gap, and still shows
the same objective-independent −0.145 Å against the untrained control.

---

## 3b. THE CASCADE, AND WHERE THE LARGEST GAP IS

The brief requires the cascade for every architecture, never a single number. This experiment
has no aggregation stage (it returns a single enumerated configuration), so the chain is
G → S → F with the enumerated space as G:

| stage | quantity | best generative objective (`E_combined`) | ORACLE restraints |
|---|---|---|---|
| **G** generation | mean RMSD over the whole k=4 space | 4.003 | 4.003 |
| | *best member of the space* | *1.030* | *1.030* |
| **S** selection, filtering | mean RMSD of the objective's best 1% | 3.032 | 2.238 |
| **S** selection, ordering | RMSD of the certified argmin | **2.767** | **1.439** |
| **F** what arms actually return at 8192 | `vqe_cvar_a0.25` / `greedy` / `untrained_bestofN` | **2.776** / 2.767 / 2.981 | **1.615** / 1.460 / 1.760 |

Gaps: **G→S(filter) = −0.97**, **S(filter)→S(order) = −0.26**, and
**S(order) → space best = 1.74 Å**, which is by a wide margin the largest and is not a search
gap at all — greedy reaches the certified argmin exactly (objective gap 0.000) and it is still
1.74 Å above the best structure the space contains. Under ORACLE restraints the same gap is
**0.41 Å**.

> **The largest gap is between the objective's own optimum and the space's best member, and it
> is 4× larger under predicted restraints than under true ones. That is the next research
> target, and it is a restraint-accuracy target, not a search or sampler target.**

This is the same conclusion `s15/coord_FINDINGS.md` K1-CORRECTED and K5 reach by a completely
different route (continuous L-BFGS on 126 targets), which is worth stating: the discrete
enumerated instrument and the continuous fitting instrument agree that the binding constraint
is distogram error.

---

## 4. HONEST LIMITATIONS

- **The 32768 budget is partly degenerate at n = 9**: it is 12.5% of a 262,144 space, and
  greedy, anneal, random and VQE all converge on the certified optimum there. The n=10 targets
  at the same budget sit at 3.1% and are the informative half of that row. This is why the
  n=9 / n=10 split is printed on every primary cell.
- **One alpha in the main ladder** (0.25). The alpha dimension is in §6.
- **`E_ORACLE_true` is an ORACLE ceiling and a diagnostic**, labelled as such in every row and
  in every sentence above. It is never a predictive result.
- **k = 4 discretisation floors everything.** The best structure in the whole enumerated space
  averages **1.030 Å** and the ORACLE objective's certified argmin is **1.439 Å**, against
  K1-CORRECTED's continuous oracle fit at 0.611 Å. Nothing measured here can beat 1.030 Å, and
  no arm came within 1.3 Å of it.
- **Only `rmsd_returned` (best-seen-by-E) is the headline readout**, which is a best-of-N
  selector for both arms. The `mode` readout — the argmax of the trained distribution, the only
  readout for which "the VQE did something" is meaningful — is in §7.
- **Machine contention.** The main ladder ran on a box at 90–100% CPU from other agents' jobs;
  `wait_for_cpu` gated the start. This affects wall-clock only, not any number.

---

## 5. MY OWN ERRORS, PRESERVED IN PLACE

Three claims of mine were overturned by more data **within this same workstream**, and the
brief's warning — *"n ≤ 4 samples from the enumerated nine have reversed a conclusion four
times"* — was vindicated twice more, once at n = 1 and once at n = 9.

**5a. REFUTED — "generative objectives get BETTER inside their own tail" (n = 1).**
On 1CS9 alone I measured `tail − bulk` at **+0.290 / +0.473 / +0.411** for three generative
objectives against Sprint 14's −0.266, cleared the matched random-tail null by 0.28, and wrote
it up as a probable sign reversal on the exact quantity Sprint 14's negative rests on. At
n = 9 the column is negative for every objective; at n = 19 it is −0.134 to −0.237. **1CS9 is
an outlier and I quoted it.** Note precisely how it fooled me: the numbers were internally
consistent, cleared their own null, and told a mechanistically satisfying story that inverted
a prior result. Nothing about the observation flagged it. Only n did.

**5b. REFUTED — "every selective objective has a positive ordering term" (n = 9).**
True on Sprint 14's own nine, where `S14_disto_bayes` read **+0.163**. At n = 19 it reads
**−0.077**. The surviving claim is narrower: the two *physical energies* have a positive
ordering term. See the correction inside §2.3.

**5c. REFUTED — "the budget-trap law holds on this objective class" (n = 9).**
rho = −0.491 on 45 cells became **−0.073** on 95. See §2.7.

**5d. NOT an error, but the design decision that saved the headline.** The handed-over design
had no selective control on the same information. With only the five handed-over objectives,
the correct reading of §2.4 would have been *"CVaR-VQE beats the untrained-circuit control at
matched budget on generative restraint objectives, CI excluding zero, 13–14 of 19"* — the first
positive quantum result in the project's history, and wrong. `S14_disto_bayes` cost about
15 minutes of compute and refuted it.

---

## 6. THE CVaR DIMENSION — and the mechanism confirmed on two more axes

**DEMONSTRATED.** `s15/results/qrestraint_alpha.json`. All 19 targets, 3 objectives, 3 seeds,
budget 8192, shots 512, matched. Six alphas on the **corrected** `baseline="const"` estimator,
plus two arms on the **recorded defective** `baseline="tail"` estimator at matched everything.
Paired against `untrained_bestofN`; negative = VQE better.

| arm | α·shots | `E_ls_pred` | `E_combined` | `E_ORACLE_true` |
|---|---|---|---|---|
| `vqe_const_a0.02` | 10.24 | −0.180 [−0.375,+0.001] | −0.167 [−0.388,+0.044] | **−0.217 [−0.415,−0.020]** |
| `vqe_const_a0.05` | 25.60 | −0.061 [−0.307,+0.173] | −0.155 [−0.376,+0.051] | −0.210 [−0.458,+0.002] |
| `vqe_const_a0.1` | 51.20 | −0.146 [−0.340,+0.059] | −0.154 [−0.329,+0.016] | −0.177 [−0.437,+0.046] |
| `vqe_const_a0.25` | 128 | −0.147 [−0.362,+0.062] | −0.205 [−0.424,+0.003] | −0.145 [−0.388,+0.050] |
| `vqe_const_a0.5` | 256 | −0.128 [−0.302,+0.023] | −0.144 [−0.295,+0.010] | −0.092 [−0.301,+0.090] |
| `vqe_const_a1.0` | 512 | −0.051 [−0.258,+0.159] | −0.167 [−0.383,+0.024] | −0.068 [−0.340,+0.146] |
| **`vqe_TAILDEFECT_a0.1`** | 51.20 | −0.113 [−0.333,+0.106] | **−0.246 [−0.456,−0.052]** | **−0.193 [−0.417,−0.011]** |
| **`vqe_TAILDEFECT_a0.25`** | 128 | −0.103 [−0.287,+0.086] | **−0.219 [−0.421,−0.029]** | **−0.266 [−0.541,−0.024]** |

### 6.1 How the best alpha relates to the objective's tail structure

**On the ORACLE objective — the only one with real in-tail ordering (in-tail rho +0.379,
z = 68.8; in-tail accuracy 0.643) — alpha matters, monotonically, and small alpha wins:**
−0.217 at α = 0.02 falling to −0.068 at α = 1.0, the only clean monotone trend in the table.
**On the predicted objectives, whose in-tail rho is +0.076 and +0.127, alpha barely matters**:
`E_combined` runs −0.144 to −0.205 across a 50× range of alpha with no trend, and `E_ls_pred`
is noisy. So the alpha dimension pays out **in proportion to how much orderable structure the
tail actually contains** — which is the tail-structure relationship the sprint asked for, and
it is a small payout: the entire alpha axis is worth ~0.15 Å against a ~1.9 Å selection gap.

This also **corrects the direction** of an older record. The brief notes that "the recorded
'small alpha collapses two AMBER variants' has the sign backwards: small alpha discriminates
better; concentration collapses them." That is confirmed here and given a mechanism (§6.3).

### 6.2 THE THREE RECORDED DEFECTS, measured rather than asserted

- **Defect 3 — `dCVaR/dp` identically zero iff `p(argmin E) ≥ alpha`, "one initialisation in
  four to twenty starts dead at α = 0.01".** Measured on every run:
  **`dead_start` = 0.000 across all 19 targets × 3 seeds × 3 objectives × 8 arms**, with
  first-iteration gradient norms of 9.8–230. **No run in this study started dead.** At α ≥ 0.02
  with this ansatz and these objectives the defect is inert. That is a scope condition on the
  Sprint 14 record, not a contradiction of it.
- **Defect 2 — sampled-CVaR upward bias at non-integer `alpha·N`.** Active at α = 0.02, 0.05,
  0.1 (α·N = 10.24, 25.6, 51.2) and inactive at 0.25, 0.5, 1.0. **Caveat that must travel with
  the result:** the single best cell in the whole sweep (`E_ORACLE_true` at α = 0.02, −0.217,
  CI excluding zero) is one where this bias **is** active. The α = 0.25 cells, which are clean,
  are weaker. I am not able to separate the two on this data and say so.
- **Defect 1 — the `baseline="tail"` gradient bias.** See below; the result is the surprise of
  the section.

### 6.3 REFUTED — the gradient-baseline defect does not cost accuracy. It HELPS.

**The recorded defective estimator is not worse on the structural axis. It is better**, and on
`E_combined` and `E_ORACLE_true` it produces the only α = 0.1 / 0.25 cells whose CI excludes
zero (−0.246, −0.219, −0.193, −0.266) where the *corrected* estimator at the same alpha gives
−0.154, −0.205, −0.177, −0.145 with CIs crossing zero.

The diagnostics say exactly why, and they close the mechanism:

| arm | first-iteration &#124;grad&#124; | `conc` (concentration) | distinct fraction |
|---|---|---|---|
| `vqe_const_a0.02` | 20.8 – 137 | −0.096 to −0.121 | 0.551 – 0.583 |
| `vqe_const_a0.25` | 20.2 – 141 | −0.194 to −0.203 | 0.505 – 0.520 |
| `vqe_const_a1.0` | 24.2 – 230 | **−0.288 to −0.291** | **0.458 – 0.472** |
| **`vqe_TAILDEFECT_a0.25`** | **9.8 – 64** | **−0.106 to −0.124** | **0.584 – 0.591** |

The biased estimator's gradient norm is **2.2–2.5× smaller** than the corrected one's at the
same alpha, so it moves the parameters less, concentrates less, and keeps more distinct
samples — and it returns better structures. Likewise alpha: larger alpha gives a larger
gradient, more concentration (−0.096 → −0.291), fewer distinct samples (0.583 → 0.472), and a
**smaller** gain.

> **The CVaR gradient-baseline defect is a de-facto step-size reduction, and on this problem
> a weaker optimiser is a better sampler.**

### 6.4 THE MECHANISM, NOW CONFIRMED ON THREE INDEPENDENT AXES

Everything in this workstream points the same way, and the three axes are causally independent
of one another:

| axis | more optimisation ⇒ | outcome |
|---|---|---|
| **budget** 2048 → 8192 → 32768 | distinct fraction 0.78 → 0.61 → 0.42 | gain −0.15/−0.20 at 8192, **gone** at 32768 |
| **alpha** 0.02 → 1.0 | conc −0.10 → −0.29, distinct 0.58 → 0.47 | gain −0.217 → −0.068 (ORACLE) |
| **gradient quality** biased → corrected | &#124;grad&#124; ×2.4, conc −0.11 → −0.20 | gain −0.266 → −0.145 (ORACLE) |

**In all three, less optimisation is better.** This is Sprint 14's *"concentration is the wrong
move"* vindicated as a mechanism on a new objective class, a new statistic and two new axes —
even though Sprint 14's *headline* (running VQE is worse than not running it) does not
reproduce. The two are consistent: concentration is harmful at the margin, but starting from
an untrained circuit a *little* concentration is still net positive, and the optimum amount is
small and objective-independent.

**And the negative cannot be blamed on a defective estimator.** The main table uses the
corrected `baseline="const"`; the defective one was run side by side and does *better*, so
neither the parity-with-classical result nor the small win over the untrained control is an
artefact of gradient bias.

---

## 7. THE `mode` READOUT — the one place the VQE clearly does something

**DEMONSTRATED.** `s15/results/qrestraint_mode.json` and
`s15/results/qrestraint_untrained_enrich.json`. Nine n=9 targets, 3 seeds, budget 8192,
α = 0.25, `exact_dist=True` so the whole 2^18 distribution is enumerated exactly rather than
sampled. `s14/vqe_run.py` states the reason this readout exists: `best_seen` is a best-of-N
selector, so it "makes the VQE a (bad) random search dressed in a circuit unless the
distribution actually concentrates", while `mode` and the distribution "are the only readout
for which 'the VQE did something' is a meaningful claim".

**The control matters more than the number, so it was measured rather than skipped.** Quoting
near-native mass against a *uniform* base rate would credit the optimiser with whatever bias
the ansatz already has at `theta0` — and `theta0` is the control this entire workstream rests
on. So the identical quantity was computed at `theta0`, on the same targets and seeds.

| distribution | entropy (bits) | effective support | P(<2.0 Å) | P(<2.5 Å) | × uniform @2.0 | **× UNTRAINED @2.0** |
|---|---|---|---|---|---|---|
| uniform over 4^9 | 18.00 | 262,144 | 0.0149 | 0.0712 | 1.00 | — |
| **UNTRAINED (`theta0`)** | 12.74 | 7,418 | 0.0126 | 0.0594 | **0.85** | 1.00 |
| trained: `E_ls_pred` | 9.46 | 931 | 0.0657 | 0.1817 | 4.41 | **5.20** |
| trained: `E_ml_pred` | 9.39 | 854 | 0.0778 | 0.2052 | 5.22 | **6.16** |
| trained: `E_ls_pool` | 9.15 | 790 | 0.0123 | 0.0699 | 0.83 | **0.98** |
| trained: `E_combined` | 9.68 | 1,025 | 0.0847 | 0.1881 | 5.68 | **6.70** |
| trained: `E_ORACLE_true` | 9.93 | 1,157 | 0.2356 | 0.3996 | 15.81 | **18.64** |

**The control strengthens the result rather than weakening it.** The untrained circuit is
*slightly worse than uniform* (0.85× at 2.0 Å) — it is not a near-native-biased prior — so the
enrichment measured against it is **larger** than against uniform: **5.2×, 6.2× and 6.7× at
2.0 Å** on the three good generative objectives, and 18.6× on the ORACLE.

**And it has a clean internal control.** `E_ls_pool`, the weakest objective in the study
(global rho 0.316, in-tail rho +0.044), gives **0.98×** — no enrichment at all. The enrichment
tracks objective quality, so it is not an artefact of optimisation per se.

**So the variational state genuinely learns something.** Training narrows the distribution from
an effective support of 7,418 to ~1,000 configurations (12.74 → 9.4–9.9 bits) and multiplies
the probability of drawing a sub-2 Å structure by 5–7×. This is the one unambiguously positive
quantum-side finding in the workstream, and it is a **distribution** property.

**And it is worth almost nothing through the readout the pipeline actually uses.** The mode
lands at 2.76–3.53 Å against a certified argmin of 2.71–3.32 and a best-seen of 1.24–1.33: the
distribution concentrates on a region whose *modal* member is ~1.5 Å worse than the best member
it drew. The arm table (§2.4) converts all of that into ≤ 0.2 Å, because `best_seen` is an
argmin over samples and an argmin does not care about mass.

> **A 5–7× enrichment in near-native probability mass is invisible to an argmin readout and
> would only be realisable by an operator that consumes the ENSEMBLE.** That is exactly
> `s15/coord_FINDINGS.md`'s **Family C** (conditional ensemble generation), and this is the
> first measurement in the project that prices it.

**The control this does NOT have, stated plainly.** Whether a *classical* sampler at matched
budget — e.g. a Boltzmann or tempered distribution fitted to the same objective — reaches the
same enrichment is **untested here**. Without it, this is a measurement of what the VQE does,
not a claim that only a VQE does it. That is the obvious next experiment and it is cheap.

---

## 8. WHAT IS OPEN

1. **A classical ensemble control for §7.** Does a tempered/Boltzmann classical sampler at
   matched budget reach the same 5–7× near-native enrichment? Until that is run, §7 measures
   what the VQE does, not what only a VQE can do. Cheapest high-value experiment remaining.
2. **An ensemble-consuming operator.** §7's enrichment is invisible to an argmin. Feeding the
   trained distribution's samples to `I.coordinate_average` / the Family C machinery is the
   only route by which any of this becomes accuracy.
3. **The 5.2–6.7× enrichment is measured at α = 0.25 only.** §6 shows concentration is monotone
   in alpha, so the enrichment/diversity trade-off almost certainly has an interior optimum that
   is *not* the one that optimises `best_seen`. Untested.
4. **Nothing here escapes k = 4.** Best structure in the whole space averages 1.030 Å;
   K1-CORRECTED's continuous oracle fit reaches 0.611 Å. The discrete instrument cannot test
   anything below its own floor.
5. **`E_ls_pool` behaves unlike the other three generative objectives** (global rho 0.316, no
   enrichment, VQE loses to greedy with a CI excluding zero at two budgets). Pool-median
   distances are a genuinely different channel and were not diagnosed here.

---

## 9. REPRODUCTION

```
python -m s15.qrestraint main           # 19 targets x 6 objectives x 3 seeds x 3 budgets
python -m s15.qrestraint alpha          # CVaR alpha sweep + the tail-baseline defect arm
python -m s15.qrestraint mode           # the `mode` readout, n=9 targets, exact_dist=True
python -m s15.qrestraint enrich         # the UNTRAINED-circuit control for the mode probe
python -m s15.qrestraint report         # summarise the checkpoint, complete or partial
python -m s15.qrestraint skill          # the mechanistic readout alone
```

Every mode checkpoints after each target (and after each budget in `main`), so an interrupted
run is readable with `report` rather than lost — the failure mode that cost Sprint 14 an entire
ansatz study. `report` regenerates Tables 0–7 from the artefact, so no number in this document
is transcribed from a console: they are all in
`s15/results/qrestraint_analysed.json`.

Wall clock on a contended box (90–100% CPU from other agents throughout): main ladder ~2.5 h,
alpha sweep ~12 min, mode probe ~4 min, untrained-enrichment control ~3 min. `wait_for_cpu`
gated every launch; the main ladder waited ~35 min and the alpha sweep ~30 min for headroom.
