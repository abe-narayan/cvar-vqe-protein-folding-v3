# s25 — RMSD LANE PRE-REGISTRATION. THE LOCATION QUESTION.

**Endpoint: mean full-chain Cα-RMSD, 126 cluster-disjoint dev targets, POINT CLOUD basis, shipped
uniform top-75 coordinate average in the medoid frame. Incumbent 3.0483 Å.**
Only the SELECTION varies in every arm below. Nothing downstream of selection is touched.

Written before any number is read. Two modules, run in order, with a stated prediction crossing
between them.

---

## 0. WHAT THE MANDATE NARROWS TO, AND THE ONE THING I AM ADDING

s25 L2 established that the shipped score is `risk(t) = Σ_c p_c |t − C_c|` — an **L1 functional whose
minimiser is the posterior MEDIAN**. Width is inert (tempering moves calibration 3× and the endpoint
0.003 Å); location is not (the prior ladder moves the endpoint at −2.15 Å per unit γ).

Every closed direction in the DO-NOT-REDO table reshapes each pair's posterior **independently**, or
reweights pairs. So does every arm below except family D. **The one structurally new thing available
is that a distance vector is not a free object: it must be realisable by a point set in R³.** That
constraint couples pairs and is the only source of location information in this system that is
neither the sequence nor the native.

### The location statistic is the whole lever, so I enumerate it as one axis

The incumbent's effective per-pair target is the risk minimiser on the grid — the posterior median.
Three one-parameter families move that statistic, each containing the incumbent as an **exact**
interior point (asserted numerically, not argued):

| family | loss `ρ(u)`, `u = t − C_c` | minimiser | identity point |
|---|---|---|---|
| **A1 POWER** | `\|u\|^q` | q→0 mode, q=1 **median**, q=2 mean | **q = 1, bit-exact** |
| **A2 TRUNC** | `min(\|u\|, δ)` | δ→0 mode, δ→∞ **median** | **δ = 40 ≥ max\|grid−C\| = 35.95, bit-exact** |
| **B QUANT** | `2·[τ·u⁺ + (1−τ)·u⁻]` | τ-quantile | **τ = 0.5, bit-exact** |
| **C SHELL** | `\|t + α·off(sep) − C_c\|` | median shifted by α·off | **α = 0, bit-exact** |
| **D METRIC** | `\|t − C_c\|` on mass moved to `(1−η)L_med + η·L_proj` | projected location | **η = 0, bit-exact** |

A1/A2 are the **shape** axis (mechanisms 2 and 3). B/C are the **shift** axis (mechanisms 1 and 5),
in two different normalisations. D is the **constraint** axis (mechanism 4).

---

## 1. MODULE 1 — `s25/loc.py`. THE LOCATION DIAGNOSTIC. NO RMSD IS COMPUTED.

Natives are read for **diagnosis only**. Nothing is fitted, nothing is selected, no production object
is modified. Same standing as `s25/calib.py`.

For every candidate location statistic `L` above, per pair, measure the move `ΔL = L − L_med` against
the direction to truth `ΔT = D_native − L_med`:

    gam_eff  =  <ΔL, ΔT> / |ΔT|²        the fraction of the way to truth actually travelled
    cos      =  <ΔL, ΔT> / (|ΔL|·|ΔT|)  how much of the move is aligned at all
    |ΔL|/|ΔT|                            how far it moved

`gam_eff` is the **prior ladder's own currency**. The ladder's TILT arm is location-only and tracks
MASS to 0.024 Å at γ=0.1, so γ is a legitimate unit for a pure location move. The ladder's γ has
`cos = 1` by construction; any real functional has `cos < 1` and its orthogonal component adds error
the ladder never carried. **Therefore −2.1496 × gam_eff is an UPPER BOUND on the endpoint gain, not a
prediction of it.** I will state it as a bound and say so at every appearance.

Three further free measurements, each of which can kill a family before it costs a run:

**(D1) THE CANCELLATION TEST — the coordinator's stated warning, made falsifiable.** s24 Workstream C
measured the distogram's signed offset at +0.4062 and the whole legal universe of **real protein
windows** at +0.4419 against the same natives, i.e. it is a corpus/native scale mismatch the
CANDIDATES SHARE. If the pool's own signed offset per shell matches the prior's, then shifting only
the prior moves the target *away* from where the candidates live and breaks a cancellation.
**PREDICTION, STATED BEFORE THE RUN: `pool_signed[shell] ≈ prior_signed[shell]` to within ~0.1 Å at
every shell, and families B and C will therefore FAIL end-to-end despite positive `gam_eff`.**
If that is what happens it is a clean result, and it is the one I expect.

**(D2) THE MULTIMODAL SUB-POPULATION.** On the 24.1% of pairs with ≥2 peaks above 0.02: is the native
closer to the NEAREST MODE than to the median? If not, mechanism 2 has no headroom and A2 is dead
before it runs.

**(D3) ZERO-INFORMATION CONTROLS, matched in the operator's space** (`control-must-match-the-
operators-space`, the project's most repeated error). Two: (i) a random location move of matched
per-pair magnitude `|ΔL|`, which is a plausible location move rather than a uniform-nonsense one;
(ii) the per-shell offset profile of a **randomly permuted other target**, which preserves the
correction's shape and destroys only its target correspondence.

---

## 2. MODULE 2 — `s25/locrun.py`. END-TO-END, NESTED CV, THE ENDPOINT.

Every arm: shipped K=500 BLOSUM pool → modified risk table → top-75 → shipped uniform coordinate
average in the medoid frame → Cα-RMSD to native, point cloud. **Only the risk table changes.**

**The risk table is never re-implemented.** A genuine `core.predict.Distogram` is constructed from the
posterior and its own `w` and `grid` are used; the kernel `|grid − CENTRES|` is replaced by `ρ`. At
the identity parameter the resulting float32 table is asserted **`np.array_equal` to the object's own
`_risk`** — already verified at max abs diff 0.0 on 1A13 — per target, before any number is read.

**One global parameter per family**, fitted by nested leave-one-fold-out CV over the 5 pinned folds.
No per-target routers (`n≈100/fold` carries a 0.39 Å generalisation gap for even a one-threshold
router). The per-shell arm C-5 (5 parameters) is reported ONLY as the optimism, never as a result.

**Every per-target minimum over the grid is reported beside `ST.best_of_k_null` with `share_accounted`,
residual and `k_eff`.** W/L is never used to diagnose one.

### The scale confound in family A, named and controlled

`|u|^q` and `min(|u|,δ)` change the location statistic AND the relative weight of one pair against
another (a pair with large typical `|u|` gains weight at q=2, loses it at q<1). So each A arm is run
twice: **raw**, and **scale-matched** — each pair's risk row rescaled so its mean over the grid equals
the q=1 row's, which is exactly the identity at the identity parameter and isolates SHAPE from
WEIGHT. Same device as `priorladder`'s MASSFIXW and `temper`'s SDFIXW.

---

## 3. RULE 0 — SIX OPERATOR FORKS, EACH NAMING THE ALTERNATIVE NOT TAKEN

Enumerated by me, who has a stake in the answer. The audit lane should re-enumerate.

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| **functional** | The shipped Bayes-risk form with `w` and `grid` taken from a genuine `core.predict.Distogram` and only the loss KERNEL replaced; bit-exact identity asserted per target. | Re-implementing the risk table or the weight `w` by hand; and modifying `w` itself, which is confidence reweighting and is CLOSED (s24 L6, s25 L2). |
| **basis** | Point cloud throughout, medoid frame, on **both sides of every contrast**. | The built-chain basis, 0.156 Å away and not comparable to the pinned 3.0483. |
| **readout** | Shipped uniform top-75 coordinate average, unchanged in every arm, so the only thing that varies is WHICH 75 candidates are selected. | Re-optimising `m`, weighting members, or clustering — s23 L5 shows uniform is the in-sample optimum and that removing geometric outliers COSTS +0.142 Å. |
| **normalisation** | Family A raw AND scale-matched, reported side by side; family B self-scaled by the posterior's own quantile; family C in absolute Ångströms off a per-shell profile. The two normalisations of the shift are run as separate arms rather than one being chosen silently. | Choosing one normalisation for the shift and reporting it as "the" location correction. |
| **null** | The identity parameter of each family (q=1, δ=40, τ=0.5, α=0, η=0), which reproduces the incumbent **bit-for-bit through this module's own code path** and is therefore also the module's self-check. Plus the two matched zero-information location moves of D3. | Comparing against the pinned constant 3.0483 alone, which would not catch a defect in my own scoring path. |
| **THE LABEL** | Continuous Cα-RMSD **and** `gam_eff` **and** the held-out location MAE, reported together. A family that improves `gam_eff` and not RMSD, or RMSD and not `gam_eff`, refutes the stated mechanism either way and I will say which. | RMSD alone, which would let a lucky ranking pass as a mechanism; and any binarised win rate, which cannot diagnose a best-of-K arm. |

---

## 4. HYPOTHESES AND FALSIFIERS, ONE PER FAMILY

| # | H | falsifier | my prior |
|---|---|---|---|
| A1/A2 | The median is the wrong location statistic for a posterior that is multimodal on 24% of pairs; a mode-seeking or mean-seeking loss ranks better. | The nested-CV arm failing to beat the identity parameter past its own MDE with a fold CI excluding zero. | **Genuinely open.** This is the only untested axis that is neither width nor a global shift. |
| B/C | The posterior's location is systematically too long and a single global shift toward shorter distances improves the ranking. | Same. **Plus the D1 cancellation test, which I predict will pre-refute it.** | **Expect FAILURE**, for the reason in D1. Stated before the run. |
| D | Projecting the location field onto the metrically realisable set is a location correction that uses no new information, only the constraint. | Same, plus `gam_eff ≤ 0` in module 1 would refute it before it runs. | **Expect small or negative.** `core/predict.py`'s own `realize` docstring already records "projecting the matrix and scoring against the projection measurably hurts", and `better-matrix-worse-ranking` item 1 quantifies it (MAE 2.13→1.96, in-band ρ 0.379→0.334). That measurement was made on the **decoy bank the project has since retracted** (`decoy-bank-not-a-pool-proxy`) and with **no interpolation parameter**, so a retest on the valid instrument with η=0 as an exact interior point is legitimate — but it is a RETEST OF A CLOSED DIRECTION and I label it that way. |

**Sprint-level falsifier.** If no family clears its own MDE, the honest output is that the location of
this posterior cannot be improved from information already available by any reshaping of the
consuming functional, and that Phase I's open question is answered NO on the functional side. That is
a result and it will be reported as one, not padded.

---

## 5. DISCIPLINE ATTACHED TO THIS LANE

- `ST.compare` / `ST.fmt` for every contrast; MDE = 2.8016×SE per comparison; fold CI beside iid;
  0.7–1.3× MDE is Type-M and is not a result.
- `ST.save_atomic(..., complete_keys=NEED, rows=rows, n_expected=126, module_file=__file__)`.
- `s25/results/LOCK_TRAIN` via `os.open(..., O_CREAT|O_EXCL)` for the end-to-end run; released
  promptly. The AMBER/OpenMM lock is the physics lane's and is not taken.
- No module is edited while a job launched from it is running.
- ORACLE / ACHIEVABLE / PRODUCTION labelled at every appearance. `gam_eff` is a DIAGNOSTIC computed
  from natives and is never a selection criterion.

---

## APPENDIX — OUTCOME, APPENDED AFTER BOTH RUNS. THE PRE-REGISTRATION ABOVE IS UNEDITED.

Full write-up in `s25/agentRMSD_FINDINGS.md`; artefacts `s25/results/loc.json`,
`s25/results/locrun.json`, both `complete: true`, n=126, provenance-stamped.

**Amendments made BEFORE the runs they affect, at the coordinator's direction, each recorded here:**
the cosine beside `gam_eff` per family and per shell; the multimodal stratification pre-registered
before run 1 rather than requested after; the grid-quantisation diagnostic (`d2centre`/`on_ctr`);
the pool's signed offset AFTER score selection as well as over the full K=500; a SIXTH family
DEQUANT added after run 1's quantisation finding, with a SEVENTH operator fork for the bin-support
convention; and families QUANT and SHELL DROPPED before running, closed by run 1's own null and by
the audit lane's A3/A8. `s25/loc.py`'s opening motivation cites `s25/LEDGER.md` L7, which retracts
the L2 framing this lane was commissioned under; the correction is carried in the module rather than
edited away.

**Every family is NOT MEASURED against its own bit-exact identity parameter, and for METRIC, POWER
and TRUNC the FULL-LEAKAGE best global parameter IS the identity.** Across 51 non-identity arms
`corr(gam_eff, endpoint delta) = +0.054`. The falsifier stated in §4 — "every family's nested-CV arm
failing to beat its own identity past its own MDE with a fold CI excluding zero" — **FIRED, on all of
them.** Phase I's open question is answered NO on the functional side.

**Prediction scorecard (four of six wrong or half wrong):** D1 half wrong; METRIC small-or-negative
confirmed; family A open on multimodality **refuted by my own run 1**; DEQUANT null; P1 half
confirmed; P2 refuted. The two that landed both damaged directions I had argued for.
