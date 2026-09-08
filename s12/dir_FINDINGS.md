# S12 — ERROR DIRECTION

**Question.** Can the SIGN of the shipped distogram's per-pair error be predicted from
deployable features, and is a sign-corrected objective worth what the ORACLE says it is?

**Instrument.** tuning126 only. benchmark60 never read; `s9/final_cache` never opened;
dev24 never run. `rr` / `nat_ca` / peptide-database native CA are used for EVALUATION and
as leave-fold-out TRAINING LABELS only. Every arm that consumes them at inference time is
named `o_*` and labelled ORACLE.

Code: `s12/dir_common.py`, `dir_repro.py`, `dir_gate.py`, `dir_gate2.py`, `dir_gatepick.py`,
`dir_esmsub.py`, `dir_corpus.py`, `dir_feats.py`, `dir_head.py`, `dir_headanal.py`,
`dir_head2.py`, `dir_emit.py`, `dir_emitavg.py`, `dir_emitreport.py`, `dir_report.py`,
`dir_conf.py`, `dir_why.py`. Results: `s12/results/dir_*.json` and `dir_*.log`.

**Status.** Complete on all 126 targets, including the projected end-to-end run
(`dir_emit.py` → `dir_emit_summary.json`), whose `bayes` arm returns 3.205 against the
production record's 3.204.

---

## D0. What the forensics agent's `o_sign` arm actually does

Read from `s12/fail_pairs.py` lines 59–81, not from the summary.

The baseline arm scores the K=500 pool with the shipped **Bayes risk** table
`risk0` (npairs × 760 grid), i.e. the L1 Bayes risk of the full 17-bin predicted
distribution. The `o_sign` arm **replaces that whole table** with

```python
tgt  = exp + 2.0 * np.sign(dtrue - exp)          # DELTA = 2.0, fixed magnitude
risk = np.abs(grid[None, :] - tgt[:, None])      # L1 risk of a POINT estimate
```

So `o_sign` is a *fixed-magnitude* shift in the correct direction — a genuine sign channel
in its second factor — but it changes **two** things relative to the baseline at once:

| factor | baseline | `o_sign` |
|---|---|---|
| (a) objective form | Bayes risk over the 17-bin distribution | L1 against a single point |
| (b) the point | — | `E[d]` moved 2.0 Å the right way |

Only (b) is direction information. The published `o_magonly` control
(`exp + |err| * random_sign`, −1.226 Å on FAIL18) shares factor (a) and carries **zero**
direction information, and it reproduces the entire `o_sign` effect (−1.217 Å) on FAIL18.
That is the first reason to suspect the FAIL18 half of this lead is factor (a).

The correct reference for a sign channel is therefore the **`pt` arm**: the same L1
point-estimate objective at `E[d]` with **no** shift. D1 measures it.

*(Note on scoring convention: `score_risk(D, |grid − tgt|)` equals
`mean_p |D_p − tgt_p|` up to the 0.05 Å grid quantisation and clipping to [2, 40] Å, so
the sign arms are plain L1 objectives.)*

---

## D1. Reproduction and de-confounding (`dir_repro.py`, `dir_repro.json`)

Instrument selfcheck first: `python -m s12.instrument` → 3.4540 / 1.7108 / 2.3062 /
3.2041 / 18 zero-recall = FAIL18. Reproduced.

FAIL18 + the 18 length-matched controls of `fail_contrast.json`, emitted through the real
path (λ = 0 synthesis arm, the 3.204 convention; λ = 0.3 stored in the same JSON).

| arm (ORACLE except `base`/`pt`) | FAIL18 | Δ vs base | Δ vs **pt** | W/L vs pt | MATCH18 | Δ vs base | Δ vs pt |
|---|---|---|---|---|---|---|---|
| `base` shipped Bayes risk | 6.034 | 0 | +0.077 | 8/10 | 2.256 | 0 | −0.083 |
| **`pt` L1 point at E[d], no shift** | **5.956** | **−0.077** | **0** | — | **2.339** | **+0.083** | **0** |
| `o_sign` Δ = 0.5 | 5.730 | −0.303 | −0.226 | 16/2 | 2.235 | −0.020 | −0.103 |
| `o_sign` Δ = 1.0 | 5.355 | −0.678 | −0.601 | 18/0 | 2.163 | −0.093 | −0.176 |
| **`o_sign` Δ = 2.0** (the published arm) | **4.784** | **−1.250** | −1.172 | 18/0 | **2.058** | **−0.198** | −0.281 |
| `o_sign` Δ = 3.0 | 4.340 | −1.694 | −1.617 | 18/0 | 2.022 | −0.234 | −0.317 |
| `r_sign` Δ = 2.0 (random sign, same Δ) | 5.926 | −0.108 | −0.030 | 12/6 | 2.448 | +0.193 | +0.109 |
| `r_sign` Δ = 3.0 | 5.866 | −0.168 | −0.091 | 13/5 | 2.613 | +0.357 | +0.274 |
| `o_magonly` (true magnitude, random sign) | 4.788 | −1.245 | −1.168 | 13/5 | 2.202 | −0.054 | −0.137 |
| `o_true` (perfect distogram) | 3.641 | −2.392 | −2.315 | 18/0 | 1.932 | −0.323 | −0.406 |

**Reproduction is exact.** The forensics agent published base 6.019, `o_sign` 4.803
(Δ −1.217, 18/0), `o_magonly` 4.793 (Δ −1.226), `o_true` 3.640 (Δ −2.379), MATCH18
`o_sign` Δ −0.199. Every one of those is matched to ≤ 0.02 Å; the residual is the λ = 0 vs
λ = 0.3 readout convention.

**The confound is real but small, and the lead survives it.** The objective-form change
(factor (a) of D0) is worth only −0.077 Å on FAIL18 and **+0.083 Å (i.e. a loss)** on
MATCH18. So of the published −1.217 Å, ~0.08 Å is the point-estimate form and ~1.17 Å is
direction. Against the honest `pt` reference the sign channel is **−1.172 Å on FAIL18
(18/0) and −0.281 Å on MATCH18**, and its own null at the same magnitude (`r_sign` 2.0) is
−0.030 / +0.109 — dead, and slightly harmful on the controls. The channel is real.

Two things the published table did not show:

* **The magnitude was not optimised.** Δ = 3.0 is worth 40 % more than Δ = 2.0 on FAIL18
  (−1.617 vs −1.172 against `pt`). D2 sweeps this properly.
* **`o_magonly` is not a direction null.** It reproduces `o_sign` on FAIL18 (−1.168)
  because the true magnitude is huge exactly on the pairs that matter, so a random-signed
  shift of that size destroys the misleading pair either way — but it is 13/5 rather than
  18/0 and it collapses on the controls (−0.137 vs −0.281). The correct direction null is
  `r_sign` at the same fixed Δ, and that one is dead.

---

## D2. The PRE-REGISTERED accuracy gate (`dir_gate.py`, `dir_gate2.py`, `dir_gatepick.py`)

All 126 targets, ORACLE sign corrupted at a fixed flip rate (0 → 0.5; 0.5 = random = the
null), 3 seeds per cell, magnitudes 0.25 → 6.0 plus the Δ → ∞ **linear limit**

```
|D_p − (E[d]_p + Δ s_p)| → Δ − s_p (D_p − E[d]_p)     ⇒  score = −mean_p s_p (D_p − E[d]_p)
```

Reported on the coordinate average (obj_FINDINGS 0b: r = 0.995 with the projected
structure, near-constant +0.16 Å), against `pt` = 3.082 (all126) / 5.747 (FAIL18) /
2.638 (other-108); `bayes` = 3.048 / 5.832 / 2.584; `o_true` = 2.214 / 3.365 / 2.023.

**Gain over `pt` (Å, all126, coordinate average):**

| Δ \ flip rate | 0.00 | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 |
|---|---|---|---|---|---|---|
| 0.25 | −0.086 | −0.069 | −0.049 | −0.038 | −0.019 | −0.001 |
| 1.0 | −0.306 | −0.244 | −0.187 | −0.122 | −0.061 | +0.001 |
| 2.0 | −0.511 | −0.442 | −0.334 | −0.222 | −0.104 | +0.032 |
| 3.0 | −0.630 | −0.562 | −0.436 | −0.295 | −0.120 | +0.081 |
| **4.0** | **−0.647** | −0.586 | −0.466 | **−0.335** | −0.111 | +0.171 |
| 6.0 | −0.561 | −0.506 | −0.405 | −0.264 | +0.033 | +0.380 |
| linear limit | −0.138 | −0.088 | +0.046 | +0.140 | +0.472 | +0.973 |

FAIL18 is monotone all the way to the linear limit (Δ = 4 −1.887, Δ = 6 −2.099,
lin −2.229); other-108 peaks at Δ ≈ 3 (−0.457) and the linear limit is a **disaster**
there (+0.210). The optimum magnitude is finite and ≈ 3–4 Å because the shift has to stay
small enough that E[d] still anchors the objective; a pure sign-matching objective throws
the anchor away and helps only the targets whose objective was already broken.

**Random sign (rate 0.5) is null everywhere** (+0.032 at Δ = 2, +0.171 at Δ = 4 on
all126) — the sweep is measuring direction, not perturbation.

### The gate

Sign accuracy required to reach a given gain over `pt` (linear interpolation in accuracy):

| all126, Δ | break-even | −0.05 Å | −0.10 Å | −0.20 Å | −0.30 Å | −0.50 Å |
|---|---|---|---|---|---|---|
| 1.0 | 0.502 | 0.583 | 0.664 | 0.823 | 0.991 | — |
| 2.0 | 0.524 | 0.561 | 0.597 | 0.682 | 0.770 | 0.984 |
| 3.0 | 0.540 | 0.565 | 0.590 | 0.646 | 0.703 | 0.851 |
| **4.0** | **0.561** | 0.578 | **0.596** | **0.640** | **0.684** | 0.828 |
| 6.0 | 0.611 | 0.628 | 0.645 | 0.678 | 0.725 | 0.894 |

On FAIL18 the gate is far lower — Δ = 4 clears −0.20 Å at chance and −0.50 Å at 56.6 %
accuracy — because there the objective is so wrong that almost any correct-direction
pressure helps. On other-108 Δ = 4 needs 59.0 % to break even and 70.6 % for −0.20 Å.

**Pre-registered decision, taken before any trained model's emitted RMSD was inspected:
Δ = 4.0**, the argmin of the ORACLE gate curve evaluated at the head's measured accuracy.
Disclosure: Δ was chosen on the all-126 ORACLE curve, i.e. in sample on the evaluation
instrument; the choice is between magnitudes whose curves are within 0.06 Å of each other
at the relevant accuracy, so it cannot manufacture the effect, but it is not free.

**So the gate is: ~56 % to break even, ~60 % for 0.10 Å, ~64 % for 0.20 Å, ~68 % for
0.30 Å (all126, Δ = 4).**

---

## D3. The sign head (`dir_corpus.py`, `dir_feats.py`, `dir_head.py`, `dir_headanal.py`)

**Training corpus.** Not the 126 tuning targets — every one of the **787 peptides** in
`peptide_db`, each given its own leave-fold-out distogram (`fold_model(fold(P))`, which
never saw fold(P)). **103,876 labelled pairs** from 787 peptides, vs ~8,500 from the 126.
ESM is served from a 787-sequence extract of `s12/esm_bank.py`; `esm_cache.npz` was never
opened. Validation: the recomputed distogram reproduces the cached shipped distogram on
all 126 tuning targets to `max|ΔE[d]| = 0.0000 Å`.

Label `y = 1` iff `d_true > E[d]`. Corpus base rate 0.4743.
124 deployable features: separation and position; E[d], sd, sd/E[d]; the shape of the
17-bin predicted distribution (entropy, mode, quartiles, IQR, skew, off-mode mass,
P(d<6/8/10/12)); within-target shell residuals and per-residue row means; the implied
radius of gyration and its ratio to a length expectation; ESM-2 contact probability, its
2×2 neighbourhood, per-residue degree and separation-shell contact mass; residue
propensities; and the shipped PCA-32 ESM embedding of residues i and j.

**Leave-fold-out accuracy** (the head scoring fold f is trained on folds ≠ f, so no
target's own errors and no error from any sequence in its identity-clustered fold ever
trains the head that judges it):

| group | n pairs | **full** | sep-only | shell-majority | constant | shuffled | AUC full | AUC sep |
|---|---|---|---|---|---|---|---|---|
| corpus (787) | 103,876 | **0.6760** | 0.5246 | 0.5414 | 0.5257 | 0.5238 | 0.741 | 0.518 |
| tuning126 | 8,549 | **0.6865** | 0.5595 | 0.5563 | 0.5562 | 0.5534 | 0.751 | 0.518 |
| **FAIL18** | 1,429 | **0.6816** | 0.6137 | 0.5479 | 0.5472 | 0.5444 | 0.747 | 0.669 |
| other-108 | 7,120 | **0.6875** | 0.5486 | 0.5580 | 0.5580 | 0.5552 | 0.752 | 0.487 |
| non-tuning (661) | 95,327 | 0.6750 | 0.5215 | 0.5401 | 0.5230 | 0.5212 | 0.740 | 0.517 |

The head is **not** chance-level on FAIL18 (0.682, AUC 0.747) — the brief's precondition
holds. Null (b) — separation alone — reaches only 0.525 corpus-wide and 0.559 on
tuning126, i.e. essentially the constant baseline; it does carry real content on FAIL18
(0.614, AUC 0.669), which is consistent with the FAIL18 compactness defect being
separation-structured, but it is 7 points below the full head there. Null (a) — shuffled
labels — collapses to the constant baseline exactly as it must.

### What the head learned (`dir_headanal.json`)

| | tuning126 | FAIL18 | other-108 |
|---|---|---|---|
| head accuracy | 0.6879 | 0.6676 | 0.6912 |
| ORACLE per-target majority sign (**one oracle bit per target**) | 0.6637 | 0.7445 | 0.6503 |
| constant (always "over-predicted") | 0.5520 | 0.5133 | 0.5584 |
| **head accuracy on MINORITY pairs** (true sign ≠ target majority) | **0.7096** | **0.7326** | **0.7062** |
| minority fraction | 0.3363 | 0.2555 | 0.3497 |
| head predictions agreeing with the head's own per-target majority | 0.644 | 0.659 | 0.642 |
| head's per-target majority equals the ORACLE target majority | 0.651 | 0.667 | 0.648 |

**This is pair-specific, not one bit per target.** A head that had only learned the global
scaling mode would score 0 % on the minority pairs; this one scores **71 %** there,
*higher* than its overall accuracy, and only 64 % of its predictions follow its own
per-target majority. It also beats the ORACLE single-bit ceiling on tuning126 (0.688 vs
0.664), which a target-level model cannot do.

**Accuracy rises monotonically with the size of the error it is predicting** — the shape
the channel needs, since the value is front-loaded onto a handful of badly-wrong pairs:

| \|error\| decile | 0.00–0.17 | 0.36–0.59 | 0.88–1.30 | 1.91–2.77 | 2.77–3.97 | 3.97–6.09 | 6.09–21.9 |
|---|---|---|---|---|---|---|---|
| accuracy | 0.519 | 0.613 | 0.703 | 0.704 | 0.752 | 0.786 | **0.835** |

Per-separation accuracy is flat at 0.65–0.73 for shells 2–13 (0.57 at sep 15, n = 21) —
the head is not a separation prior in disguise.

### How this differs from the S7-3 de-biasing that failed

S7-3 fitted a *magnitude* correction: regress the true distance on the predicted one and
shrink. It died because the residual is not a scaled version of the truth (calibration
slope +0.376, the predictor is uncorrelated rather than shrunk), so there is no scalar
to apply. The head here never estimates how wrong a pair is, never touches calibration,
and is scored by a metric that has no scale in it at all — a binary label whose base rate
is 0.47. The fixed Δ is not a fitted magnitude either: it is chosen on the ORACLE-sign
gate, i.e. it is a property of the *operator's* sensitivity, not of the error distribution.
Whether that distinction buys anything downstream is D4.

## D4/D5. End to end, with all three nulls (`dir_emitavg.py`, `dir_report.py`)

All 126 targets, every arm at every magnitude, coordinate-average proxy, paired against
the **shipped Bayes incumbent** (`bayes` = 3.048 avg / 3.204 projected). `pt` = 3.082.
Negative = better. `rand_acc` = ORACLE sign corrupted to the head's own per-target
accuracy (3 seeds) — **null (c)**. `sep` = separation-only head — **null (b)**.
`shuf` = head trained on shuffled labels — **null (a)**. `const` = the training folds'
majority sign applied everywhere. `o_tmaj` = ORACLE one-bit-per-target majority sign.

### all126 (n = 126), Δ emitted vs the incumbent

| Δ | `o_sign` ORACLE | **`head`** | `head_soft` | **`rand_acc` (c)** | `sep` (b) | `shellmaj` | `const` | `shuf` (a) | `o_tmaj` | `head_tmaj` |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.0 | −0.272 | +0.028 | +0.021 | −0.078 | −0.001 | +0.074 | +0.073 | +0.070 | −0.084 | +0.079 |
| 2.0 | −0.478 | +0.066 | +0.022 | −0.138 | +0.075 | +0.188 | +0.178 | +0.166 | −0.076 | +0.217 |
| 3.0 | −0.597 | +0.157 | +0.042 | −0.169 | +0.164 | +0.298 | +0.289 | +0.278 | +0.009 | +0.398 |
| **4.0** (pre-registered) | **−0.614** | **+0.291** | +0.074 | **−0.142** | +0.297 | +0.415 | +0.409 | +0.398 | +0.185 | +0.623 |
| 6.0 | −0.528 | +0.526 | +0.170 | +0.016 | +0.491 | +0.563 | +0.563 | +0.524 | +0.427 | +0.896 |
| **W/L at Δ = 4** | 98/28 | **45/81** | 55/71 | 61/65 | 40/86 | 39/87 | 41/85 | 41/85 | 50/76 | 34/92 |
| **drop-top-10 at Δ = 4** | −0.422 | +0.456 | +0.170 | +0.005 | +0.454 | +0.545 | +0.540 | +0.525 | +0.365 | +0.776 |

95 % CI at Δ = 4: `head` **[+0.124, +0.455]**, `rand_acc` [−0.270, −0.025],
`o_sign` [−0.781, −0.453]. Per-fold `head` @4: {0: +0.43, 1: +0.16, 2: −0.04, 3: +0.51,
4: +0.38} — harmful in 4 of 5 folds.

### FAIL18 (n = 18) — `bayes` 5.832, `pt` 5.747

| Δ | `o_sign` | **`head`** | `head_soft` | **`rand_acc`** | `sep` | `shellmaj` | `const` | `shuf` | `o_tmaj` |
|---|---|---|---|---|---|---|---|---|---|
| 2.0 | −1.248 | −0.379 | −0.239 | −0.457 | −0.436 | −0.243 | −0.236 | −0.234 | −0.691 |
| **4.0** | **−1.972** | **−0.741** | −0.338 | **−0.961** | −0.698 | −0.383 | −0.385 | −0.367 | −0.884 |
| 6.0 | −2.185 | −1.016 | −0.479 | −1.202 | −0.891 | −0.544 | −0.527 | −0.512 | −0.935 |
| W/L @4 | 18/0 | 14/4 | 10/8 | 14/4 | 15/3 | 15/3 | 15/3 | 14/4 | 17/1 |
| drop-top-10 @4 | −1.191 | **+0.165** | +0.193 | −0.112 | +0.008 | +0.115 | +0.127 | +0.141 | −0.196 |

CI @4: `head` [−1.278, −0.256], `rand_acc` [−1.495, −0.505], `o_sign` [−2.425, −1.558].

### other-108 (n = 108) — `bayes` 2.584, `pt` 2.638

| Δ | `o_sign` | **`head`** | `head_soft` | **`rand_acc`** | `sep` | `shuf` |
|---|---|---|---|---|---|---|
| 2.0 | −0.349 | +0.140 | +0.065 | −0.085 | +0.160 | +0.232 |
| **4.0** | **−0.387** | **+0.463** | +0.143 | **−0.006** | +0.462 | +0.525 |
| 6.0 | −0.252 | +0.783 | +0.278 | +0.219 | +0.721 | +0.697 |
| W/L @4 | 80/28 | 31/77 | 45/63 | 47/61 | 25/83 | 27/81 |

### Reading

1. **The gate was cleared on accuracy and failed on structure.** The head measured 0.6865
   accuracy; the pre-registered gate said ≥ 0.684 buys −0.30 Å at Δ = 4. It emits
   **+0.291 Å [+0.124, +0.455], 45W/81L** — 0.60 Å on the wrong side of the prediction,
   with the sign reversed. This is the project's standing "better predictor, worse
   structure" pattern, reproduced on a channel specifically designed to avoid it.
2. **Null (c) is the whole story and it fires against the head.** ORACLE signs corrupted
   to the head's *own accuracy* emit **−0.142 Å [−0.270, −0.025]**. Same accuracy, 0.43 Å
   apart, opposite signs. It is not that a 68 %-accurate sign channel is worthless — it is
   that *this* head's 32 % of mistakes are the damaging kind. Read the other way: the only
   arm at the head's accuracy that helps is the one that cannot be built.
3. **Null (b) is indistinguishable from the head downstream.** `sep` (55.9 % accuracy)
   emits +0.297 vs the head's +0.291 on all126, and +0.462 vs +0.463 on other-108. Thirteen
   points of per-pair accuracy buy **0.00 Å**. The r_sep = 0.196 result generalises: what
   the objective consumes is what survives after separation, and the head added nothing
   there that the operator can use.
4. **Null (a) does not vanish, and that is informative.** `shuf` collapses to the constant
   sign, and the constant sign is not a no-op: +0.398 on all126, −0.367 on FAIL18. On
   FAIL18 *shrinking every predicted distance by 4 Å regardless of anything* is worth
   −0.37 Å. So of the head's −0.741 on FAIL18, half is the shuffled-label null. Any
   FAIL18 claim must be net of `const`/`shuf`, and this one nearly isn't.
5. **The FAIL18 gain is concentrated and dies to drop-top-10** (−0.741 → **+0.165**),
   while `o_sign`'s survives (−1.972 → −1.191). Brief rule 8 disposes of it.
6. `head_soft` (shift Δ(2p−1), confidence-weighted) is uniformly less harmful than the
   hard head and never helpful (+0.074 all126, −0.338 FAIL18): shrinking the correction
   towards zero recovers `pt`, which is the same thing as saying the correction has no
   value.

### Why the head is worse than random at its own accuracy (`dir_why.py`)

Per target, the head's *signed mistake field* vs a random field of the same rate (5 draws):

| group | arm | err rate | coherence¹ | top-1 eigen share | eff. rank | global-scale share | **MAE of the corrected objective** |
|---|---|---|---|---|---|---|---|
| all126 | head | 0.3121 | 0.999 | 0.3158 | 6.62 | **0.0889** | **3.658** |
| all126 | random | 0.3121 | 0.951 | 0.2741 | 7.28 | 0.0527 | 3.916 |
| FAIL18 | head | 0.3324 | 1.011 | 0.3259 | 6.82 | 0.1271 | 4.721 |
| FAIL18 | random | 0.3324 | 0.978 | 0.2843 | 7.57 | 0.1051 | 4.992 |

¹ P(both wrong | pairs share a residue) / rate², chance = 1.0.

The head's mistakes are **more coherent than i.i.d.**: 1.7× more of the squared mistake
sits in the single global "expand/contract everything" mode (0.089 vs 0.053), the leading
eigenvector takes a larger share (0.316 vs 0.274) and the effective rank is lower
(6.62 vs 7.28). That is the same defect obj_FINDINGS 1 measured in the distogram itself,
and the same reason a displacement-field corruption is far more damaging than i.i.d. noise
at equal MAE: a coherent wrong field moves the whole coordinate average, while independent
flips cancel in it. A sign head trained on ESM and distogram features inherits the error
structure of the very model it is correcting.

**And the divergence is explicit in the metrics.** The head's corrected objective has a
**lower** MAE than the random-sign one (3.658 vs 3.916 Å) and emits **0.43 Å worse**
structures. Note also that both corrected objectives are much worse distance predictors
than the uncorrected `pt` (MAE 2.339) while the ORACLE-sign version of the same
construction emits −0.61 Å: **the value of this channel has nothing to do with the
accuracy of the resulting distance objective, only with how it orders the pool.**

### D4b. Confirmation through the REAL projection (`dir_emit.py`, `dir_emitreport.py`)

Full path: score K=500 → top-75 → superpose on medoid → coordinate average → L-BFGS
projection. `fit` is the λ = 0 synthesis arm (the 3.204 convention); `lam` is λ = 0.3.
Δ = 4.0, paired against the shipped Bayes incumbent, **all 126 targets**. The `bayes` arm
returns **3.205** through this path against the production record's **3.204** — an
independent end-to-end reproduction of the incumbent.

| arm | avg | **fit** | lam | **Δ fit vs bayes** | 95 % CI | W/L | drop-10 | **MAE** | **r_sep** | top-75 best |
|---|---|---|---|---|---|---|---|---|---|---|
| `bayes` | 3.048 | **3.205** | 3.213 | 0 | — | — | — | 2.339 | 0.196 | 2.306 |
| `pt` | 3.082 | 3.246 | 3.264 | +0.041 | [−0.005,+0.088] | 57/69 | +0.082 | 2.339 | 0.196 | 2.314 |
| **`o_sign` ORACLE** | 2.435 | **2.577** | 2.597 | **−0.628** | [−0.800,−0.466] | 99/27 | −0.428 | 2.792 | **0.537** | 1.821 |
| **`head`** | 3.339 | **3.511** | 3.513 | **+0.306** | [+0.128,+0.481] | 42/84 | +0.476 | 3.658 | **0.078** | 2.554 |
| `head_soft` | 3.122 | 3.293 | 3.305 | +0.088 | [−0.014,+0.188] | 55/71 | +0.184 | **2.199** | 0.197 | 2.332 |
| `sep` (null b) | 3.345 | 3.517 | 3.519 | +0.312 | [+0.150,+0.467] | 43/83 | +0.479 | 4.264 | 0.196 | 2.419 |
| `o_tmaj` ORACLE 1 bit | 3.233 | 3.443 | 3.442 | +0.238 | [+0.054,+0.420] | 55/71 | +0.421 | 3.686 | 0.195 | 2.238 |
| `head_tmaj` | 3.671 | 3.871 | 3.863 | +0.666 | [+0.489,+0.845] | 37/89 | +0.821 | 4.268 | 0.196 | 2.733 |
| **`rand_acc` (null c)** | 2.907 | **3.069** | 3.080 | **−0.136** | [−0.283,+0.008] | 63/63 | +0.030 | 3.910 | 0.227 | 2.117 |

**FAIL18, complete (n = 18)** — `bayes` 6.034, exactly the D1 value:

| arm | fit | Δ vs bayes | 95 % CI | W/L | drop-10 | r_sep |
|---|---|---|---|---|---|---|
| `o_sign` ORACLE | **4.119** | **−1.915** | [−2.448,−1.430] | **18/0** | **−0.978** | 0.449 |
| `rand_acc` (null c) | 5.120 | −0.914 | [−1.571,−0.320] | 12/6 | +0.250 | 0.160 |
| `o_tmaj` ORACLE 1 bit | 5.251 | −0.783 | [−1.228,−0.395] | 16/2 | −0.069 | 0.089 |
| **`head`** | **5.312** | **−0.721** | [−1.344,−0.209] | 13/5 | **+0.224** | 0.080 |
| `sep` (null b) | 5.397 | −0.637 | [−1.077,−0.246] | 15/3 | +0.080 | 0.088 |
| `head_soft` | 5.713 | −0.321 | [−0.640,−0.030] | 10/8 | +0.182 | 0.140 |
| `pt` | 5.956 | −0.077 | [−0.202,+0.019] | 10/8 | +0.075 | 0.089 |

On the subgroup the whole lead was built for, the head is **beaten by one ORACLE bit per
target** (−0.783), **beaten by its own matched-accuracy random control** (−0.914), clears
the separation-only null by **0.084 Å**, and **does not survive drop-top-10** (+0.224).
Other-108: `bayes` 2.734, `o_sign` 2.320 (−0.414, 81/27), `rand_acc` 2.727 (−0.007),
**`head` 3.210 (+0.477, 29/79)**, `sep` 3.204 (+0.470) — head and null within 0.007 Å.

The projection reproduces the coordinate-average verdict exactly, including the sign
reversal between `head` and its matched-accuracy random control.

**MAE and r_sep make the divergence explicit, in three different directions:**

* `o_sign` has a **worse** MAE than the incumbent (2.780 vs 2.324) and emits **−0.641 Å**.
  What it improves is r_sep: 0.176 → **0.546**. Correcting the direction restores the
  *within-shell ordering* that obj_FINDINGS 1 showed the distogram has almost none of.
* **The head does the opposite: it DESTROYS within-shell information, r_sep 0.176 →
  0.079**, despite 0.688 per-pair sign accuracy. Its sign predictions are largely a
  function of quantities the objective already contains, so applying them subtracts
  pair-specific signal instead of adding it. That is the mechanism, stated in the
  objective's own currency.
* And r_sep is not a sufficient statistic either: `head_soft` has a **better** MAE (2.159)
  *and* a better r_sep (0.203) than the incumbent and still emits **+0.127 Å**, and `sep`
  has exactly the incumbent's r_sep (0.176) and emits +0.329 Å. No intermediate metric
  measured here prices the emitted structure; only the ORACLE arm does.

### D5d. The remedy implied by the mechanism also fails (`dir_head2.py`)

If the head hurts because its mistakes are coherent, deny it every target-global and
row-context feature (rg, rg_rel, global means, shell residuals, row means, ESM degree and
shell mass — 104 of 124 features kept) and retrain leave-fold-out.

| | accuracy | global-scale share | top-1 | eff. rank | Δ@2 | Δ@3 | **Δ@4** | rand@4 |
|---|---|---|---|---|---|---|---|---|
| full head, all126 | 0.6879 | 0.0889 | 0.316 | 6.62 | +0.066 | +0.157 | **+0.291** | −0.142 |
| **local head, all126** | 0.6868 | 0.0941 | 0.317 | 6.55 | +0.062 | +0.174 | **+0.356** | −0.190 |
| local head, FAIL18 | 0.6787 | 0.1088 | 0.319 | 6.90 | −0.370 | −0.693 | −0.848 | −1.144 |
| local head, other-108 | 0.6882 | 0.0917 | 0.317 | 6.49 | +0.135 | +0.319 | +0.557 | −0.031 |

Same accuracy, same coherence, same harm. **The coherence is not carried by the global
features** — it is intrinsic to predicting a pair's error sign from the pair's own local
distogram and ESM features, because adjacent pairs share residues and therefore share
those features. There is no feature-pruning fix.

### D5c. Confidence gating does not rescue it (`dir_conf.py`)

Apply the head's sign only to the top-q fraction of pairs by |p − 0.5|; the matched null
is ORACLE signs corrupted to the head's *conditional* accuracy **inside the same gate**.
Δ emitted vs the incumbent, coordinate average.

| q | head acc. in gate | head @3 | rand @3 | head @4 | rand @4 |
|---|---|---|---|---|---|
| **all126** | | | | | |
| 0.10 | **0.890** | +0.028 | +0.021 | +0.023 | +0.010 |
| 0.25 | 0.834 | +0.017 | −0.015 | +0.040 | −0.024 |
| 0.50 | 0.785 | +0.060 | −0.095 | +0.116 | −0.105 |
| 1.00 | 0.688 | +0.157 | −0.175 | +0.291 | −0.142 |
| **FAIL18** | | | | | |
| 0.10 | 0.790 | −0.167 | −0.172 | −0.166 | −0.169 |
| 0.50 | 0.752 | −0.348 | −0.488 | −0.423 | −0.658 |
| 1.00 | 0.668 | −0.590 | −0.770 | −0.741 | −1.003 |
| **other-108** | | | | | |
| 0.10 | 0.907 | +0.061 | +0.054 | +0.054 | +0.040 |
| 1.00 | 0.691 | +0.282 | −0.076 | +0.463 | +0.002 |

The head reaches **89 % accuracy on its most confident decile** — well past every gate in
D2 — and that arm is worth **+0.023 Å**, i.e. nothing, because ten per cent of the pairs
is not enough coverage to move the objective (its matched random control is +0.010, also
nothing). As coverage grows the channel acquires value (the random control goes to
−0.142) and the head simultaneously acquires its coherent-error penalty, and the two cross
so that the head is never better than its matched control at any coverage or any Δ.
**There is no operating point at which this head is deployable.**

## D6. The ceiling, and what fraction of it the head captured

Coordinate-average scale (add ≈ +0.16 Å for the projected structure; the projected run
`dir_emit.py` confirms the headline rows separately).

| | all126 | FAIL18 | other-108 |
|---|---|---|---|
| shipped Bayes incumbent | 3.048 | 5.832 | 2.584 |
| `pt` (L1 at E[d]) | 3.082 | 5.747 | 2.638 |
| **ORACLE sign, Δ = 4** | **2.434** | **3.860** | **2.197** |
| ORACLE sign, Δ → best (6.0 / lin on FAIL18) | 2.520 | 3.647 | 2.181 |
| **perfect distogram (`o_true`)** | **2.214** | **3.365** | **2.023** |
| K=500 pool best member (ORACLE selection) | 1.711 | 2.284 | 1.616 |
| the head | 3.339 | 5.091 | 3.047 |

**The channel's own ceiling.** Distance-channel headroom from the incumbent to a perfect
distogram is 0.834 Å on all126 and 2.467 Å on FAIL18. A **perfect** error-direction head
captures **73.6 %** of it on all126 (0.614/0.834) and **79.9 %** on FAIL18
(1.972/2.467). So direction alone is nearly the whole distance channel — the forensics
agent's framing is confirmed and, at the optimised magnitude, strengthened (they measured
51 % at Δ = 2).

**What the trained head captured: −35 % on all126** (it moved the wrong way), and on
FAIL18 30.0 % gross, **15.2 % net of the `shuf`/`const` null**, which does not survive
drop-top-10.

**And the channel is capped below what would matter anyway.** Even a perfect distogram —
strictly above any sign head — emits 2.214 Å (avg) / 2.395 Å (projected, obj_FINDINGS 0)
on all126 and 3.365 / 3.641 Å on FAIL18, against a 2.284 Å FAIL18 pool best and a 2.285 Å
ORACLE-single-best emitted. The remaining 1.36 Å on FAIL18 is the operator/selection gap
and no distance model closes it. So the honest bound on this line of work is:

* a **perfect** sign head is worth ≈ −0.61 Å on tuning126 and cannot reach 2.0 Å;
* a **perfect distogram** is worth ≈ −0.83 Å and cannot reach 2.0 Å either;
* the gap that decides < 2.0 Å is selection, which this channel does not touch.

---

## D7. Leakage audit

* benchmark60 never touched: no `results/benchmark_manifest.json`, no `s9/final_cache/*`,
  no `bench_results/cache` read other than `I.shipped_record` inside `instrument.selfcheck`
  (the production tuning record). **dev24 was not run.**
* `esm_cache.npz` was **never loaded**. ESM comes from `s12/cache/dir_esm_pep.npz`, a
  787-sequence extract of the compact `s12/esm_bank.py` bank built in a throwaway process.
  Peak RSS of every process here stayed under ~0.5 GB; `free_gb()` was checked before every
  target and the job waited whenever it was under 1.5 GB.
* `rr` / `nat_ca` / peptide-database native CA appear only in (a) evaluation, (b) arms named
  `o_*` (`o_sign`, `o_true`, `o_tmaj`, `o_magonly`, `rand_acc`), (c) leave-fold-out training
  **labels**. `rand_acc` is a DIAGNOSTIC, not a deployable arm: it starts from the oracle
  sign. The deployable arms (`head`, `head_soft`, `sep`, `shellmaj`, `const`, `shuf`) read
  only the target sequence, the shipped LFO distogram and ESM-2 features of the target.
* Fold discipline: the head scoring a target in fold f is trained only on peptides in folds
  ≠ f (folds pinned by `peptide_folds.json`, keyed by sequence, verified to agree with
  `instrument.targets()` on all 126). Second-order disclosure: a training peptide P in fold
  g ≠ f is labelled using `fold_model(g)`, which *did* see fold f. That leaks nothing about
  the scored target's structure — the target's native never enters — but the head is
  learning the error pattern of a distogram trained on a slightly different set than the one
  it will correct. Note it; it cannot manufacture the negative result reported here.
* Selection-bias disclosure: Δ = 4.0 was chosen on the all-126 ORACLE gate curve at the
  head's measured accuracy, i.e. in sample on the evaluation instrument, before the head's
  emitted RMSD was looked at. The full Δ ladder is reported, and the head is harmful at
  every Δ ≥ 1 on all126, so the choice does not drive the conclusion.
* Tie-breaking: no `argmin` on a tied signal is taken anywhere; `np.argsort(..., kind=
  "stable")` on continuous L1 scores, as in the shipped pipeline.

---

## D8. Answer, and what the coordinator should do next

**Q: Can the sign of the distogram's per-pair error be predicted?** **Yes, clearly.**
0.688 leave-fold-out accuracy on tuning126 against a 0.552 constant baseline, AUC 0.751,
0.682 on FAIL18, 0.835 on the largest-error decile, and it is genuinely pair-specific
(0.710 accuracy on the pairs whose true sign contradicts their own target's majority).
Trained on 103,876 pairs from all 787 corpus peptides under strict fold discipline.

**Q: Is a sign-corrected objective worth what the oracle says it is?** **No.** The ORACLE
channel is worth even more than published once the magnitude is optimised (−0.614 Å on
tuning126 at Δ = 4, 98W/28L; −1.972 Å on FAIL18, 18/0; 74–80 % of the entire
perfect-distogram headroom). The trained head, at an accuracy that clears the
pre-registered gate, emits **+0.291 Å [+0.124, +0.455], 45W/81L** — worse than the
incumbent, and 0.43 Å worse than ORACLE signs corrupted to its own accuracy.

**Verdict: not deployable, at any magnitude, at any confidence coverage, on any subgroup
that survives drop-top-10.** Recommend closing the error-direction head.

### Why this is not the S7-3 failure repeated, and what it is instead

S7-3 died on *calibration* — there was no scalar shrinkage that helped because the
predictor is uncorrelated with the truth rather than shrunk. This head has no calibration
in it: the label is binary, the base rate is 0.47, and its accuracy is a real 0.688. It
died on a different and, for this project, more important mechanism: **error coherence.**
A model trained on the distogram's own features (E[d], sd, distribution shape, ESM) makes
mistakes that are *correlated in the same low-rank, global-scaling way the distogram's own
errors are* — 0.089 of the squared mistake in the single global expand/contract mode
against 0.053 for i.i.d. flips at the same rate, top-1 eigen share 0.316 vs 0.274,
effective rank 6.62 vs 7.28. Coherent wrongness moves the coordinate average bodily;
independent wrongness cancels in it. The corrected objective is *more accurate by MAE*
(3.658 vs 3.916) and emits *worse structures*, which is the project's standing pattern
in its sharpest form yet.

The generalisable lesson: **the correlation structure of a corrector's errors is a
first-class design constraint, and a corrector trained on the features of the thing it
corrects will inherit that thing's error structure.** Any future per-pair corrector needs
an explicit decoherence objective or an independent feature channel, and it needs to be
evaluated against random-corruption-at-matched-accuracy, not against chance.

### Ranked recommendations

1. **Close the per-pair error-direction head** (forensics recommendation 1). It is the
   only item in `fail_FINDINGS.md` that was flagged as untested-and-promising, and it is
   now tested and negative. The channel is real (ORACLE −0.61 Å) and unreachable by a
   head trained on distogram/ESM features.
2. **Retire "predictor accuracy" as an intermediate target for this pipeline, explicitly.**
   Three independent measurements now say it does not price emitted RMSD: MAE (S7-3,
   S7-11), in-band rank correlation (S11), and now per-pair sign accuracy — where 13
   points of accuracy over the separation-only null bought exactly 0.00 Å. The only
   intermediate metric that has ever tracked emitted RMSD here is the ORACLE arm itself.
3. **If a corrector is attempted again, add an error-decoherence constraint and the
   matched-accuracy random control.** The matched control is cheap (it is just a corrupted
   oracle), it is the only null that caught this, and it should be mandatory for every
   future "learned correction" arm in this project. Note that feature pruning is *not* a
   decoherence mechanism: D5d removed every global and row-context feature and the
   coherence statistics did not move at all.
4. **The quantity a distance corrector must improve is r_sep, not MAE or accuracy.** The
   ORACLE sign arm has a *worse* MAE than the incumbent and emits −0.64 Å; it is the only
   arm that lifts within-shell correlation (0.176 → 0.546). The head, at 0.688 sign
   accuracy, *lowers* it to 0.079. If a future corrector is built, gate it on r_sep before
   paying for a projection run — while remembering (`sep`, `head_soft`) that r_sep is
   necessary, not sufficient.
5. **Do not spend dev24 on this.** Nothing here warrants confirmation: the deployable arm
   is harmful on the tuning instrument with a CI excluding zero on the wrong side.
6. **The standing ceiling is unchanged and should be quoted as such:** a perfect distance
   objective emits 2.395 Å projected on tuning126 and 3.641 Å on FAIL18 against a 2.285 Å
   ORACLE-selection floor. Selection, not the distance objective, is what stands between
   this system and 2.0 Å. Every distance-side channel — including the best-shaped one
   anyone had — is bounded by that.

