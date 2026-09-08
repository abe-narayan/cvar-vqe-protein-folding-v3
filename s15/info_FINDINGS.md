# SPRINT 15 — INFO workstream findings

Owner: INFO (Phase 1, information half). Written continuously; newest sections appended.

**Instrument selfcheck, run at the start AND at the end of the workstream**
(`python -m s12.instrument`, `s15/results/instrument_end.log`) reproduces exactly:

    shipped 3.4540004952559396   pool_best 1.7108244199364904
    top75_best 2.3061526409453816   synthesis_fit 3.2040761603809194   n_zero_recall 18

**Reproduction commands** (each writes JSON to `s15/results/` and caches to `s15/cache/`):

    python -m s15.info_sign        # B.1  the skill-definition decomposition, 19 targets
    python -m s15.info_chain       # B.2  proxy -> sign -> emitted RMSD, enumerated space
    python -m s15.info_signpool    # B.3  the same chain on the 126 real pools
    python -m s15.info_channels    # A    the channel table + the `sd` column
    python -m s15.info_phase       # C.1/C.2  the phase surfaces (~6 min)
    python -m s15.info_regime      # C.3  where the real channels sit on the surface
    python -m s15.info_errstruct   # C.3  surrogate destruction of the error structure
    python -m s15.info_null        # C.4  alignment with the RMSD-quiet subspace
    python -m s15.info_mi          # D    information content + higher-order structure
    python -m s15.info_figure      # -> s15/figures/info_phase_diagram.png

Shared helpers: `s15/info_lib.py` (geometry caches off the cached enumerations, circular and
rank statistics, and `concentration_verdict`, which emits ONE null-calibrated PASS/FAIL
against a simulated uniform-effect null instead of fields a reader can choose between).

All threads capped at 2 (`OMP/MKL/OPENBLAS_NUM_THREADS`); no job here allocates more than a
few hundred MB, and `esm_cache.npz` (1.5 GB) was never loaded — only `s12/esm_bank.py`.

---

## B. THE SKILL-DEFINITION DISCREPANCY — RESOLVED

`python -m s15.info_sign` → `s15/results/info_sign_grid.json`

### B.0 What was actually in conflict

Sprint 14 carried one number in two states: OBJ's **+0.909** (per-target skill vs native Rg
z-scored against the enumerated space, 12 targets, leave-one-**target**-out, skill =
`rho_global` on the full space) and the coordinator's **+0.416, p = 0.177** (skill =
`cross_rho` from the ceiling experiment, band-trained and band-scored, compactness
length-residualised). **Five** things differ, not one. I varied them one at a time on a
common target set, using the geometry caches rebuilt directly from the enumerations'
own `PHI`/`PSI` tables.

**Bit-level check that the reconstruction is the same instrument.** My independently
recomputed space statistics reproduce OBJ's §2.3 table exactly on all 12 shared targets
(e.g. 1CS9 natRg 5.257 / spaceRg 5.433 / z −0.20 against the recorded 5.26 / 5.43 / −0.21;
2MK7 8.375 / 5.428 / +3.91 against 8.37 / 5.43 / +3.92).

### B.1 The decomposition — DEMONSTRATED

n = 19 enumerated targets (9 at n=9 from `s13/results`, 10 at n=10 from `s14/cache`),
Spearman, 10,000-permutation p.

| per-target skill statistic | vs native-Rg z (space) | vs length-residual | vs raw | perm p (z) |
|---|---:|---:|---:|---:|
| **model-free `rho(Rg, RMSD)` over the space** | **−0.935** | −0.914 | −0.946 | **0.000** |
| model-free `rho(Rg, RMSD)` inside the band | −0.900 | −0.914 | −0.861 | 0.000 |
| **Legacy `rho_global`** | **−0.775** | −0.823 | −0.747 | **0.000** |
| LFO learned model, `rho_global` (full space) | +0.154 | +0.125 | +0.288 | 0.531 |
| LFO learned model, in-band pair acc (<1.5 A) | **−0.409** | −0.489 | −0.307 | 0.084 |
| LFO learned model, **−d_top100 (emitted RMSD)** | −0.028 | −0.079 | +0.119 | 0.915 |
| CEILING `cross_rho` (12 targets) | +0.406 | +0.343 | +0.392 | 0.193 |
| CEILING cross-target in-band pair acc (12) | +0.294 | +0.273 | +0.350 | 0.358 |
| CEILING −d_top100 in band (12) | +0.056 | −0.098 | −0.077 | 0.872 |

**The resolution, in three statements.**

1. **The coordinator's +0.416 replicates exactly.** My independent implementation returns
   Spearman **+0.406** / Pearson **+0.421** (z-scored) and **+0.416** (length-residualised)
   for `cross_rho`, p = 0.19–0.28. The number was not an implementation error.
2. **The +0.909 is a different quantity, and it is very close to a tautology.** The
   magnitude only appears on a *structural* skill axis: `rho(Rg, RMSD)` over the space
   correlates **−0.935** with the native Rg z-score (p = 0.000, n = 19). That link is
   near-definitional — if the native is more extended than the space mean, extended
   structures are nearer the native, so `rho(Rg, RMSD)` must be negative. OBJ's +0.909 is
   this tautology inherited through the intermediate claim "the learned model ≈ a
   compactness detector". It is a statement about geometry, not about a learned model's
   transferability, and **it should not be quoted as a ceiling on a conditioning signal.**
3. **On the actual leave-fold-out learned model over 19 targets the link is absent**:
   `rho_global` vs z = **+0.154, p = 0.531**. The +0.909 does not survive the change from
   leave-one-target-out on 12 to leave-fold-out on 19.

**Two further facts the discrepancy was hiding.**

- **The two skill statistics have OPPOSITE signs on the same arm.** For the LFO learned
  model, `rho_global` correlates **+0.154** with extendedness while in-band pair accuracy
  correlates **−0.409** (−0.489 length-residualised, p = 0.032). Global rank correlation is
  better on extended natives; near-native ordering is better on compact ones. "The
  per-target sign" is therefore **not a single well-defined quantity** — it depends on which
  ordering statistic you ask about, and the two that matter most anti-correlate. Any future
  claim about "the sign" must name the statistic.
- **Compactness does not predict the emitted structure at all.** Against `−d_top100`, the
  quantity a selection operator actually emits, every compactness normalisation gives
  |rho| ≤ 0.12 with p ≥ 0.63 (learned arm, n = 19) and |rho| ≤ 0.10 with p ≥ 0.77 (ceiling
  arm, n = 12). **The chain breaks at the last link, not the first.**

**Tier: DEMONSTRATED (a negative).** Recorded beside, not over, the sprint-14 entries.
`s14/LEDGER.md` N1 is answered: the common definition is **the emitted quantity**, and on it
neither +0.909 nor +0.416 is the operative number — the correlation is ~0.

### B.2 The direct chain, measured end to end — proxy → sign → emitted RMSD

`python -m s15.info_chain` → `s15/results/info_chain.json`. 19 enumerated targets, operator
= top-100 of a 60,000-configuration uniform sample (or of the RMSD ≤ 2.5 A band), objective
`E_s = s·Rg` and `E_s = s·Legacy`; the sign rule's threshold **and** its orientation are
fitted leave-one-target-out, so nothing about the scored target is used.

| region / axis | random | always + | always − | coin | fixed (LOO) | **ORACLE sign** | sign channel |
|---|---:|---:|---:|---:|---:|---:|---:|
| full space, ±Rg | 4.002 | 5.899 | 3.938 | 4.918 | 3.938 | **3.486** | −0.452 [−0.838,−0.143] W/L 6/0 |
| full space, ±Legacy | 4.002 | 3.887 | 6.002 | 4.945 | 3.887 | **3.853** | −0.035 [−0.094,+0.000] W/L 2/0 |
| **in-band (≤2.5 A), ±Rg** | 2.199 | 2.371 | 2.146 | 2.258 | 2.428 | **2.089** | −0.339 [−0.564,−0.159] W/L 19/0 |
| **in-band (≤2.5 A), ±Legacy** | 2.199 | 2.251 | 2.107 | 2.179 | 2.251 | **1.982** | −0.269 [−0.496,−0.087] W/L 9/0 |

Native-free sign rules, in-band, against the LOO fixed sign:

| proxy | ±Rg emitted | sign acc | vs fixed | ±Legacy emitted | sign acc | vs fixed |
|---|---:|---:|---|---:|---:|---|
| distogram-predicted Rg | **2.090** | **0.947** | −0.338 [−0.563,−0.158] W/L 18/0 | **2.006** | 0.789 | −0.245 [−0.475,−0.059] W/L 7/2 |
| retrieval-pool mean Rg | 2.163 | 0.842 | −0.265 [−0.465,−0.111] W/L 16/0 | 2.106 | 0.684 | −0.145 [−0.352,+0.022] |
| incumbent emitted Rg | 2.090 | 0.947 | −0.338 [−0.563,−0.158] | 2.006 | 0.789 | −0.245 [−0.475,−0.059] |
| ORACLE native Rg z (ceiling) | 2.098 | 0.947 | −0.330 [−0.555,−0.147] | 2.005 | 0.895 | −0.246 [−0.476,−0.061] |

On the enumerated band this is the **first native-free per-target sign signal the program
has measured with a CI excluding zero**: the distogram-predicted Rg captures **99.7%** of the
ORACLE sign channel on the compactness axis and **91%** on Legacy, and it *equals* the
ORACLE compactness variable. On the full space it captures nothing (−0.130 [−0.554,+0.286]),
because there the fixed sign is already right on 68% of targets.

### B.3 …and it does NOT transfer to the real pools — I refute my own result

`python -m s15.info_signpool` → `s15/results/info_signpool.json`. Identical logic, run on
the 126-target shipped K=500 BLOSUM pools; band = `pool_best + 1.5 A`; 2,000-permutation
null on the sign rule.

| region / operator | fixed (LOO) | ORACLE sign | sign channel | distogram proxy | acc | perm p |
|---|---:|---:|---|---|---:|---:|
| full pool, top-75 | 4.378 | 4.038 | −0.340 [−0.519,−0.185] | **+0.431 [+0.136,+0.725]** (harmful) | 0.603 | 0.004 |
| full pool, top-24 | 4.566 | 4.124 | −0.441 [−0.676,−0.237] | **+0.677 [+0.280,+1.077]** (harmful) | 0.603 | 0.004 |
| in-band, top-75 | 2.607 | 2.569 | −0.038 [−0.053,−0.025] | −0.008 [−0.026,+0.009] | 0.587 | 0.050 |
| in-band, top-24 | 2.668 | 2.585 | −0.083 [−0.108,−0.060] | −0.020 [−0.051,+0.010] | 0.587 | 0.067 |

**The whole ORACLE sign channel on the real in-band problem is 0.038–0.083 A**, an order of
magnitude below the 0.339 A it is worth on the enumerated band, and the distogram proxy
recovers 22–24% of it with a **CI that crosses zero**. Even the ORACLE native Rg only buys
−0.056 A [−0.084,−0.030] at top-24. On the full pool the proxy is actively **harmful**,
because the fixed sign is already correct on 83% of targets and the proxy overrides it.

**This is `decoy-bank-not-a-pool-proxy.md` happening to my own headline.** The enumerated
k=4 space has a different compactness distribution from a BLOSUM-retrieved window pool
(pool members are real fragments and already near-native in shape), so the sign has much
less room to matter there. Reported as a refutation of B.2's apparent significance, not as a
positive.

**Bug-grade caveat, stated rather than buried.** The `retrieval pool mean Rg` proxy is
**degenerate on the pool instrument**: `z_rg_pool = (pool mean Rg − pool mean Rg)/sd ≡ 0`, so
its LOO rule collapses to the majority sign and it reproduces the fixed arm exactly (`vs
fixed` = 0.000 with W/L 0/0). It is not a working proxy there and its 0.833 "accuracy" is
just the majority-class rate. On the enumerated space it is a real proxy (the space mean Rg
differs from the pool mean Rg).

### B.4 Verdict on the open thread the sprint inherited

> **Is there a native-free quantity, computable at inference, that predicts the per-target
> ordering sign?** — **Yes in principle, no in practice.** The distogram-predicted radius of
> gyration predicts the sign at 0.947 accuracy on the enumerated band and captures
> essentially all of an ORACLE sign channel there; but on the instrument that decides the
> paper, that entire channel is worth **0.083 A at most**, the proxy recovers **0.020 A with
> a CI crossing zero**, and applied outside the band it **costs 0.43–0.68 A**. Tier:
> **REFUTED as a route to 2.5 A.** The per-target sign is not where the missing factor of two
> lives.

---

## A. THE CHANNEL TABLE — every target-specific signal on one yardstick

`python -m s15.info_channels` → `s15/results/info_channels.json`, `info_channels.log`

Instrument: the shipped K=500 BLOSUM pool of each of the 126 tuning targets. Band =
`pool_best + 1.5 A`. Metric = **in-band pairwise ordering accuracy** (ties at 0.5), the axis
the sprint-14 arithmetic is stated on (**0.638 required for 2.0 A through a top-100
operator**; chance 0.500; sprint-14 best measured 0.539 for Legacy).

| channel | availability | in-band acc | CI95 | global acc | emit@24 | emit@75 |
|---|---|---:|---|---:|---:|---:|
| **distogram Bayes risk, low-sd half** (NEW) | 126/126 | **0.566** | [0.547,0.586] | 0.704 | 3.536 | 3.599 |
| distogram L2 on the low-sd quartile (NEW) | 126/126 | 0.560 | [0.539,0.579] | 0.678 | 3.653 | 3.701 |
| distogram L2 on the low-sd half (NEW) | 126/126 | 0.559 | [0.538,0.579] | 0.705 | 3.549 | 3.604 |
| distogram L2 weighted by 1/sd² (NEW) | 126/126 | 0.548 | [0.525,0.570] | 0.714 | 3.520 | 3.587 |
| **distogram Bayes risk (the incumbent)** | 126/126 | 0.545 | [0.522,0.567] | 0.718 | **3.501** | 3.551 |
| distogram plain L2 on `expected` | 126/126 | 0.519 | [0.495,0.542] | 0.713 | 3.589 | 3.619 |
| ESM-2 contact map agreement | 126/126 | 0.513 | [0.494,0.531] | 0.626 | 4.384 | 4.118 |
| pool-consensus secondary structure | 126/126 | 0.513 | [0.498,0.527] | 0.611 | 4.070 | 4.080 |
| radius of gyration (compactness) | 126/126 | 0.497 | [0.475,0.517] | 0.598 | 4.566 | 4.378 |
| pool typicality (consensus/medoid criterion) | 126/126 | 0.495 | [0.465,0.526] | 0.672 | 3.723 | 3.790 |
| random score (calibration null) | — | 0.498 | [0.490,0.507] | 0.499 | 4.492 | 4.462 |
| *reference*: Legacy / learned obj. / 1-local prior / AMBER (sprint 14) | | 0.539 / 0.523 / 0.501 / 0.463 | | | | |

One target has a band too small to score (125/126 usable).

### A.1 THE `sd` COLUMN — used for the first time, and it works — DEMONSTRATED

The shipped distogram has carried a per-pair `sd` column since Sprint 8 that no experiment
has ever read.

**It ranks its own error very well and its scale is wrong by a factor of 2.6.**

| sd quintile | n | mean sd | rms error | z sd |
|---|---:|---:|---:|---:|
| [0.08, 0.58] | 1710 | 0.44 | **1.34** | 3.17 |
| [0.58, 0.96] | 1710 | 0.76 | 2.23 | 2.90 |
| [0.96, 1.46] | 1709 | 1.20 | 3.09 | 2.58 |
| [1.46, 2.23] | 1710 | 1.82 | 4.13 | 2.24 |
| [2.23, 8.62] | 1710 | 3.09 | **5.93** | 1.89 |

`corr(sd, |error|) = +0.473`, Spearman **+0.540**, over 8,549 pairs on 126 targets. Pooled
`z = (true − expected)/sd` has sd **2.633**, not 1.000, and the mis-calibration is itself
sd-dependent (3.17 → 1.89 across the quintiles), so **`sd` is a usable RANKING of pair
reliability and is not usable as a likelihood without recalibration**.

**Conditioning the shipped score on it is free and it works on the metric that binds:**
restricting the shipped Bayes-risk score to the confident half of the pairs improves in-band
ordering accuracy by **+0.021 [+0.010, +0.034]**, 74 targets better / 51 worse, median
+0.011, concentration **PASS** (p = 0.79 against a simulated uniform-effect null,
mean/sd −0.32). At 0.566 it is **the best in-band ordering number this project has
measured**, above Legacy's 0.539 and the incumbent's 0.545.

**And it does not convert.** The same arm's emitted RMSD at top-24 is **+0.033
[−0.051, +0.116] WORSE**, W/L 55/69. The whole family behaves this way — every low-sd arm
gains in-band accuracy and loses (or ties) emitted RMSD. This is `in-band-is-the-only-
ranking-metric.md` and `nothing-ranks-within-the-pool.md` reproduced on a new channel: **the
distogram `sd` column is a real discriminator that the terminal operator does not consume.**
Tier: **DEMONSTRATED** for the accuracy gain, **REFUTED** as a route to lower RMSD through
the incumbent operator.

### A.2 The other channels, with their limits stated

- **The retrieval pool (500 BLOSUM windows, coordinates AND torsions).** 126/126, no
  missingness. Best member 1.711 A, top-75 best 2.306 A, one random top-75 window emits
  3.594 A. Per-torsion error of its circular mean: **RMS phi 48.1 / psi 82.6 deg**
  (MAE 45.3). Its *discriminative* content is nil: pool typicality is at chance in-band
  (0.495 [0.465, 0.526]) — consistent with `consensus-is-outlier-avoidance.md`, the
  consensus operator is an outlier filter, not an in-band ranker.
- **ESM-2 (`s12/esm_bank.py`, 100 MB; the 1.5 GB cache was NOT loaded).** Contact map
  in-band 0.513 [0.494, 0.531] — the CI touches chance; global 0.626. As a *probability*
  it is catastrophic: see D.4.
- **Predicted secondary structure** (pool-consensus H/E/C): in-band 0.513 [0.498, 0.527],
  emit@24 4.070 A. Marginal.
- **Compactness / radius of gyration**: in-band 0.497, i.e. nothing, and emit@24 4.566 A —
  worse than the random-score null. See section B for the per-target-sign version.
- **Chemical shifts**: closed by sprint-14 arithmetic (54/126 targets; ORACLE-perfect
  torsions on all 54 leave 2.021 A; 55 are needed). Not re-opened.
- **Sequence composition / physicochemical descriptors**: no per-structure score exists, so
  they cannot be scored on this yardstick; their generative content is measured in D.1
  (the retrieval channel is the sequence channel's best realisation) and by sprint 14's
  `phi carries no sequence signal at peptide length`.

### A.3 The one-line verdict for Part A

**No available channel reaches 0.638.** The best is 0.566 and it does not convert to RMSD.
The channel table's ceiling is **0.072 accuracy short** of what 2.0 A needs through a top-100
operator, and the gap has not moved by adding the two channels nobody had read (`sd` and ESM
contacts).

---

## C. THE PHASE DIAGRAM — coverage x uncertainty -> predictive RMSD

`python -m s15.info_phase` → `s15/results/info_phase.json` + `.log`;
figure `s15/figures/info_phase_diagram.png` (`python -m s15.info_figure`).

**Construction (ORACLE DIAGNOSTIC).** Covered residues get the native torsion plus noise
from the stated error model; uncovered residues are filled from the **native-free** top-75
retrieval pool, so coverage 0 *is* the retrieval channel and coverage 1 with sigma 0 is the
ideal-geometry build floor. 126 targets, 24 reps, common random numbers across the sigma
axis. Structures built with `I.build_ca`, the bit-exact builder the pipeline projects onto.

### C.1 The surface (i.i.d. noise, uniform dropout, pool fallback)

Mean CA-RMSD, 126 targets:

| sigma \ coverage | 0.0 | 0.3 | 0.5 | 0.7 | 0.8 | 0.9 | 1.0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | 4.025 | 3.829 | 3.523 | 3.084 | 2.598 | 1.664 | **0.347** |
| 10 | 4.025 | 3.865 | 3.607 | 3.242 | 2.841 | 2.155 | 1.302 |
| 15 | 4.025 | 3.910 | 3.703 | 3.391 | 3.078 | 2.503 | 1.862 |
| 20 | 4.025 | 3.963 | 3.790 | 3.581 | 3.293 | 2.851 | 2.410 |
| 30 | 4.025 | 4.100 | 4.025 | 3.930 | 3.781 | 3.526 | 3.314 |
| 50 | 4.025 | 4.337 | 4.462 | 4.501 | 4.452 | 4.437 | 4.422 |
| 80 | 4.025 | 4.599 | 4.761 | 4.869 | 4.905 | 4.916 | 4.921 |

**Contours — the largest per-torsion sigma that still reaches each level:**

| level \ coverage | ≤0.7 | 0.8 | 0.9 | 1.0 |
|---|---|---:|---:|---:|
| 3.0 A | unreachable | 13.3 deg | 22.0 | **26.1** |
| 2.5 A | unreachable | unreachable | 15.0 | **20.9** |
| 2.0 A | unreachable | unreachable | 7.6 | **16.3** |

**Four anchors that verify the surface against numbers measured elsewhere.**

- The ideal-geometry build floor is **0.347 A** (mean; median 0.272, max 1.474). No
  torsion-space generator on this manifold can beat it.
- Coverage 0 with a **constant alpha-helix** fill emits **4.065 A** — the brief's control,
  reproduced to three decimals by an independent construction.
- The incumbent's 3.204 A sits at **sigma-equivalent 28.6 deg** on this surface against the
  sprint-14 record of 27.1 deg (5% agreement, independent implementation).
- 2.0 A at full coverage needs **16.3 deg** here against sprint 14's **15.1 deg** (8%
  agreement, independent implementation). Both are i.i.d. figures. **Both are the wrong
  instrument for real channels**; see C.3.

### C.2 WHERE the gaps fall beats HOW MANY there are, by 1.79 A — DEMONSTRATED

At **sigma = 0** (perfect torsions wherever the channel reaches), 126 targets:

| coverage | gaps at the chain ENDS | one contiguous interior gap | gaps scattered uniformly | gap centred mid-chain |
|---|---:|---:|---:|---:|
| 0.5 | **2.765** | 3.549 | 3.523 | 3.758 |
| 0.6 | **2.260** | 3.213 | 3.342 | 3.595 |
| 0.7 | **1.716** | 2.859 | 3.084 | **3.505** |
| 0.8 | **1.014** | 2.332 | 2.598 | 3.212 |
| 0.9 | **0.393** | 1.728 | 1.664 | 2.624 |

**1.79 A of spread at coverage 0.7 between the cheapest and the most expensive missingness
geometry.** And it changes what is reachable at all: with terminal gaps, **2.0 A is
reachable at coverage 0.7** (sigma 7.9 deg) and 2.5 A at coverage 0.6; with a mid-chain gap,
2.5 A is unreachable at **any** coverage below 1.0 at **any** sigma.

Sprint 14 recorded terminal dropout as "0.40-0.50 A cheaper than uniform". That was measured
at high coverage. **Across the coverage axis the effect is 3-4x larger than recorded**
(1.37 A at coverage 0.7, 1.58 A at 0.8), and it is the single largest lever in the diagram.
Recorded beside, not over, the sprint-14 entry.

*Methodological note recorded rather than buried:* an earlier "contiguous" mask in this
module was mathematically identical to the terminal mask (a contiguous COVERED block puts the
gaps at the ends). It was returning identical numbers under two names; it was replaced by a
genuinely distinct contiguous-interior-gap model rather than double-counted.

Fill choice barely matters next to gap geometry: at coverage 0.5, pool fill 3.523 A,
constant-helix fill 3.477 A, uniform-random fill 4.768 A. **A constant alpha-helix fill is as
good as a retrieval-pool fill** — 0.05 A apart across the coverage axis.

### C.3 THE PHASE DIAGRAM OVER-PRICES EVERY REAL CHANNEL — stated explicitly

`python -m s15.info_regime`, `s15.info_errstruct`, `s15.info_null`.

Read each real channel's measured per-torsion RMS angular error off the i.i.d. surface at
coverage 1 and compare with the RMSD it really emits:

| channel | RMS phi | RMS psi | RMS | surface says | actually emits | gap |
|---|---:|---:|---:|---:|---:|---:|
| retrieval pool-500 circular mean | 48.8 | 86.3 | 70.1 | 4.793 | 4.175 | −0.618 |
| retrieval top-75 circular mean | 48.1 | 82.6 | 67.6 | 4.761 | 4.072 | −0.688 |
| retrieval top-75 sim-weighted | 55.5 | 85.8 | 72.2 | 4.820 | 3.514 | −1.306 |
| top-75 medoid window | 55.1 | 85.7 | 72.0 | 4.817 | 3.314 | −1.503 |
| constant alpha-helix (control) | 52.5 | 95.4 | 77.0 | 4.882 | 4.065 | −0.817 |
| incumbent projected torsions | 61.7 | 92.0 | 78.4 | 4.899 | 3.215 | −1.685 |
| **ORACLE best pool window** | 56.5 | 71.4 | 64.4 | 4.714 | **1.773** | **−2.941** |

**Per-torsion angular error is not a sufficient statistic for CA-RMSD.** A fragment with
64 deg RMS torsion error emits 1.77 A where an i.i.d. field of the same magnitude emits
4.71 A. Therefore **sprint 14's "per-torsion sigma required for 2.0 A = 15.1 deg", and my
own 16.3 deg, are SUFFICIENT conditions, not necessary ones.** A channel does not have to
reach 15 deg to reach 2.0 A; it has to get the error *direction* right.

**What the structure is — surrogate destruction, 126 targets, 24 reps.** Take a channel's
real error field, destroy one property at a time, re-add it to the native torsions, rebuild:

| channel | real | rot_res | sign_flip (paired) | perm_res | i.i.d. Gaussian matched |
|---|---:|---:|---:|---:|---:|
| retrieval top-75 circular mean | 4.072 | 4.208 | 3.922 | 4.379 | 4.817 |
| incumbent projected torsions | 3.215 | 4.319 | 4.093 | 4.659 | 4.937 |
| **ORACLE best pool window** | **1.773** | 3.883 | 3.369 | 4.243 | **4.751** |

Paired costs vs the real field (bootstrap CI): for the ORACLE window, sign flipping alone —
which preserves every |error| at every position exactly — costs **+1.596 [+1.385, +1.813]**;
cyclic rotation **+2.110 [+1.918, +2.311]**; residue permutation **+2.470 [+2.305, +2.646]**;
an i.i.d. Gaussian of matched sd **+2.978 [+2.790, +3.174]**. For the incumbent projection:
+0.878 / +1.104 / +1.444 / +1.722, every CI excluding zero. For the retrieval circular mean
the sign-flip term is **not significant** (−0.150 [−0.366, +0.070]) — circular averaging has
already destroyed the directional coherence, which is the mechanism behind sprint 14's
"coordinate averaging beats torsion averaging by 1.024 A".

**And none of it is visible in a pairwise correlation.** On the same fields:
`r(dphi_i, dpsi_i)` ranges −0.089 to +0.007; the peptide-plane compensating pair
`r(dpsi_i, dphi_i+1)` −0.072 to −0.026; chain lag-1 +0.002 to +0.155. Sprint 14's statement
that every real emitter sits in a near-i.i.d. band is **confirmed and is non-informative**:
the dependence worth up to 2.9 A is invisible to second-order statistics. Any future error
model fitted to measured correlations will be wrong by that much.

**The mechanism, named and measured.** Linearise the CA trace about the native torsions,
`J = d(superposed CA)/d(theta)` by central differences at 1 deg, and define
alignment(e) = ||J e|| / (||e|| · sqrt(trace(JᵀJ)/2n)) = damage per unit angular error
relative to a random direction.

| channel | alignment | CI95 | frac of targets < 1 |
|---|---:|---|---:|
| **ORACLE best pool window** | **0.565** | [0.494, 0.645] | 0.841 |
| incumbent projected torsions | 0.665 | [0.596, 0.739] | 0.794 |
| top-75 medoid window | 0.693 | [0.629, 0.757] | 0.817 |
| constant alpha-helix (control) | 0.709 | [0.646, 0.771] | 0.794 |
| retrieval pool-500 circular mean | 0.727 | [0.661, 0.796] | 0.810 |
| retrieval top-75 circular mean | 0.738 | [0.670, 0.810] | 0.746 |
| retrieval top-75 sim-weighted | 0.774 | [0.701, 0.846] | 0.722 |
| one random top-75 window | 0.800 | [0.728, 0.872] | 0.746 |
| *null*: random Gaussian direction | **0.945** | (sd 0.319) | — |
| *null*: the ORACLE window error, signs flipped | 0.726 | (sd 0.393) | — |

Every CI excludes the random-direction null. **Real torsion errors concentrate in the
RMSD-quiet directions of torsion space.** For the best fragment that is a first-order
discount of 0.565; the observed damage ratio is 0.373 (1.773 / 4.751), so **about two-thirds
of the discount is linear subspace alignment and one-third is nonlinear.** Stated as a
limitation rather than smoothed over. (The random-direction null returns 0.945 rather than
1.000 because the statistic is a ratio of norms and Jensen's inequality biases it downward at
2n = 18–32 dimensions; every channel is compared to that same 0.945, not to 1.)

Note the control: a **constant alpha-helix** already scores 0.709. Alignment below 1 is
therefore not by itself evidence of a good channel — it must clear the helix, and only the
ORACLE window (0.565) and the incumbent projection (0.665) do so decisively.

**Synthetic coherence is not the same thing and points the other way.** AR(1) and class-bias
noise — the usual "correlated error" models — need a *smaller* sigma than i.i.d. to reach a
level (2.0 A: i.i.d. 16.3, AR(1) rho=0.5 14.5, AR(1) rho=0.8 13.6, residue-class bias 14.9,
**SS-class bias 13.1**, the most damaging of all). **Correlated is not aligned**; no
synthetic model in this study produces alignment, and only alignment helps. Sprint 14's
"error coherence is worth nearly a factor of two in required sigma" is directionally right
(coherent error is more expensive) but the factor I measure is **1.24**, not 2.

**Sparse gross outliers.** With the remaining residues *exact*, a fraction f of residues
carrying a uniform-random torsion emits, at full coverage: f=0.05 → 1.513 A, 0.10 → 2.303,
0.15 → 2.929, 0.20 → 3.398, 0.30 → 4.093. **One bad residue in ten costs 1.96 A** against the
0.347 A floor; 2.5 A tolerates about **11%** gross outliers and 2.0 A about **8%**.

### C.4 THE PHASE DIAGRAM'S VERDICT ON 3.0 / 2.5 / 2.0 A

- **3.0 A — supported.** Needs coverage ≥ 0.8 at sigma ≤ 13 deg, or coverage 1.0 at
  sigma ≤ 26 deg, or **coverage 0.5 with terminal gaps at sigma ≤ 11 deg**. The incumbent
  sits at 28.6 deg-equivalent, so 3.0 A is roughly a 10% improvement in the effective
  channel. **Plausible with what exists.**
- **2.5 A — hard, not closed.** Reachable only at coverage ≥ 0.9 **and** sigma ≤ 15 deg
  i.i.d., or coverage ≥ 0.6 if the gaps are terminal. The best channel measured is 67.6 deg
  RMS at alignment 0.74. The honest route is **not** a smaller sigma; it is a better-aligned
  error, which fragment retrieval already supplies and which circular averaging destroys.
- **2.0 A — not supported by any i.i.d.-class channel.** It needs coverage 1.0 at
  sigma ≤ 16.3 deg or coverage 0.9 at 7.6 deg, and no available channel is within a factor
  of four of that. The only measured object that reaches 1.77 A is the ORACLE best pool
  window — a *selection* result, not a generation one, and numerically the 1.711 A pool best
  the program has had since Sprint 10. **A 2.0 A architecture built on driving sigma down is
  not supported by the information. A 2.0 A architecture built on SELECTING well-aligned
  fragments is arithmetically possible and is bounded by the 0.638 in-band ordering accuracy
  that Part A's channel table tops out 0.072 short of (best 0.566).**
- **A model that reaches 2 A only under unrealistically perfect constraints is not a 2 A
  solution.** Two arms in this study are exactly that and are labelled: (a) the sigma = 0
  column of C.2 — 2.0 A at coverage 0.7 needs *exact* torsions on 70% of residues *and* the
  gaps at the chain ends, which no channel supplies; (b) the ORACLE best-pool-window row of
  C.3, which presumes the selection problem already solved.

---

## D. INFORMATION CONTENT — with the estimator checked and a structural measure beside it

`python -m s15.info_mi` → `s15/results/info_mi.json` + `.log`

The primary estimator is **predictive**: `I_pred = H0 − CE(channel, truth)` in bits, where
`H0` is the cross-entropy of a target-blind marginal built from the OTHER 125 targets
(leave-one-target-out) and `CE` is the channel's own held-out cross-entropy on the native
value. It is a valid lower bound on the mutual information, it cannot be inflated by sampling
noise, and it measures the quantity the project needs — *bits the target-specific channel adds
over generic Ramachandran* — rather than raw association.

### D.1 The retrieval pool adds ZERO bits over generic Ramachandran — DEMONSTRATED

| bins | bin width | I_pred (bits/torsion) | CI95 | target-shuffled null |
|---:|---:|---:|---|---:|
| 6 | 60 deg | **−0.112** | [−0.229, −0.002] | −0.698 |
| 8 | 45 deg | −0.110 | [−0.231, +0.006] | −0.747 |
| 12 | 30 deg | −0.191 | [−0.318, −0.066] | −0.803 |
| 18 | 20 deg | −0.176 | [−0.295, −0.060] | −0.769 |
| 24 | 15 deg | −0.219 | [−0.336, −0.101] | −0.780 |
| 36 | 10 deg | −0.225 | [−0.329, −0.120] | −0.699 |

**Negative at every bin count**, with CIs excluding zero at four of six. The target-shuffled
null is far worse (−0.70 to −0.80), so the pool distribution IS target-specific — it is just
not *better* than a generic peptide Ramachandran prior at predicting the native torsion.

Independent cross-check in a second metric: at 30 deg bins the pool's modal bin hits the
native bin **0.371** of the time against the generic modal bin's **0.360** (chance 0.083);
at 10 deg bins **0.168** against **0.176**. Point accuracy is a wash; as a *distribution* the
pool is 0.18–0.23 bits worse because it is more diffuse.

**Paired structural measure, and it agrees exactly:** the top-75 circular mean emits
**4.072 A** and a zero-information constant alpha-helix emits **4.065 A**. Two independent
measurements, one conclusion: **at the per-residue marginal level the retrieval channel
carries no usable target-specific torsion information.** Its value is entirely in the joint
structure (D.5) and in the coordinates it supplies for selection.

### D.2 The plug-in estimator says the opposite — a measured demonstration of the trap

| bins | plug-in MI | Miller–Madow corrected | shuffle null | null sd |
|---:|---:|---:|---:|---:|
| 6 | 0.403 | 0.408 | 0.006 | 0.002 |
| 12 | 0.500 | 0.518 | 0.024 | 0.003 |
| 24 | 0.679 | 0.730 | 0.104 | 0.006 |
| 36 | 0.791 | 0.881 | 0.204 | 0.007 |

The plug-in MI between the pool's modal bin and the native bin is 0.40–0.88 bits and is
enormously "significant" against its shuffle null — while the predictive estimator on the
same data says **negative**. The plug-in is measuring the *deterministic association* of two
variables that both concentrate in the alpha basin: it re-measures generic Ramachandran and
labels it target-specific information. Note also that the plug-in rises monotonically with
bin count while its own null rises 34-fold — the signature of a bias-dominated estimator.
**Reported as the honest answer to "where is mutual information unreliable here": exactly
this setting.**

### D.3 The distogram is NEGATIVE-information as a probability, positive as a ranking

`I_pred = −0.818 bits/pair`, CI **[−1.265, −0.419]**, median −0.138, positive on only
**52/126** targets, against the pooled distance marginal of the other 125 targets. Split by
its own uncertainty: **low-sd half −0.570 bits, high-sd half −1.070 bits.**

The baseline here is *unconditional* on |i−j|, which makes it weaker than a fair baseline
would be; the distogram loses to it anyway, so the negative is conservative.

**Paired structural measure:** the same channel is nevertheless the best in-band *ranker*
available (0.545 shipped, 0.566 on its own low-sd half; A.1). **A channel can be a good
ranking function and a bad probability at the same time**, and this one is. It is the
mechanism behind `error-shape-not-mae-decides-ranking.md` (calibration slope +0.376) stated
in bits, and it is the reason `sd` recalibration — not a better point prediction — is the
cheap Phase-2 move.

### D.4 ESM-2 contact maps are catastrophically over-confident on 9–16mers

`I_pred = −2.570 bits/pair`, CI [−2.710, −2.438], positive on **0/126** targets, against a
base contact rate of 0.407. Paired structural measure: in-band ordering accuracy 0.513
[0.494, 0.531] — at chance. Tier: **REFUTED** as a usable channel at this chain length.
This is about the ESM *contact map*; it does not touch
`esm-adds-nothing-for-short-peptides.md`'s overturned status for ESM *embeddings* in
selection, which is a different object.

### D.5 The pool's higher-order structure — worth −0.457 A on the mean, +0.245 A on the best

500 draws per target, 126 targets:

| | mean | best-of-500 |
|---|---:|---:|
| whole windows (the full joint) | **3.564 A** | 2.293 A |
| independent per-residue resampling (marginals only) | 4.021 A | **2.048 A** |
| paired difference (joint − marginal) | **−0.457 [−0.604, −0.322]**, W/L 95/31, conc PASS | **+0.245 [+0.131, +0.357]**, W/L 24/102, conc PASS |

**The pool's joint structure is worth 0.457 A on the average member and costs 0.245 A on the
best attainable member.** Per-residue recombination of the *same* pool produces a strictly
better achievable pool (2.048 A vs 2.293 A, and vs the 2.306 A top-75 pool best) at the same
draw budget, and a strictly worse typical member.

Read through `operator-consumes-set-mean.md` (`d_out = 1.16·d_set_mean + 0.04·d_set_best`),
recombination is predicted to cost `1.16(+0.457) + 0.04(−0.245) = +0.52 A` through the
incumbent m=75 averaging operator and to gain `−0.245 A` through a sharp argmin. **So
per-residue recombination is a Phase-2 option only if the selector is good, and is harmful
otherwise** — the same conditional as `The budget trap is a property of BAD OBJECTIVES`.

---

## THE ONE JUDGEMENT PART D WAS ASKED FOR

**Does the available target-specific signal contain enough information to reach 3.0 / 2.5 /
2.0 A?**

- **3.0 A: YES**, with what exists. Effective channel needs ~10% improvement over the
  incumbent's 28.6 deg-equivalent, and three independent levers of that size are measured
  here (gap placement, `sd` conditioning, recombination-with-a-selector).
- **2.5 A: MARGINAL, and only through selection.** Generation is 4x short on the i.i.d.
  axis and the honest lever is error *alignment*, not error size. The selection side tops
  out at 0.566 in-band accuracy against 0.638 required, and the best in-band gain measured
  in this workstream does not convert to RMSD.
- **2.0 A: NOT SUPPORTED** by any enumerated channel. Every route that reaches it in this
  study requires an oracle — perfect torsions on 70% of residues with terminal gaps, or the
  best pool window already identified. The information is present in the pool (1.711 A pool
  best, 2.048 A recombined best-of-500); **what is absent is any measured signal able to
  find it**, and Part A shows the gap is 0.072 of ordering accuracy wide with nothing left
  unread to close it.

## WHAT I REFUTED, INCLUDING MY OWN HYPOTHESES

1. **My own B.2 headline.** The native-free per-target sign signal (0.947 accuracy, 99.7% of
   the ORACLE channel on the enumerated band) **does not transfer** to the 126 real pools,
   where the whole ORACLE channel is 0.038–0.083 A and the proxy recovers 0.020 A with a CI
   crossing zero — and is +0.43 to +0.68 A **harmful** outside the band.
2. **That per-torsion sigma prices a generation channel.** Refuted: a 64 deg-RMS fragment
   emits 1.77 A where a 64 deg i.i.d. field emits 4.71 A.
3. **That "near-i.i.d. lag-1 autocorrelation" licenses i.i.d. error modelling.** Confirmed
   the autocorrelations and refuted the inference: the dependence worth 2.9 A is invisible at
   second order.
4. **That the retrieval pool's per-residue torsion prior carries target-specific bits.**
   Refuted at every bin count (−0.11 to −0.23 bits vs generic Ramachandran), and the
   structural twin agrees (4.072 A vs a 4.065 A constant helix).
5. **That plug-in mutual information is safe here.** Refuted with the counter-example in the
   same table: +0.88 bits plug-in where the predictive bound is −0.23.
6. **That the distogram `sd` column, once used, buys RMSD.** It buys the best in-band
   ordering accuracy the project has measured (0.566) and **not** RMSD (+0.033 at top-24).
7. **That ESM-2 contact maps are a usable structural channel at peptide length.** Refuted:
   −2.57 bits/pair, 0/126 targets positive, in-band accuracy at chance.
8. **That "+0.909" was a measurement of a conditioning ceiling.** It is a near-tautology
   between two structural quantities; the leave-fold-out learned model's own skill correlates
   +0.154 (p = 0.53) with the same variable.

## OPEN AFTER THIS WORKSTREAM

- **Recalibrating the distogram `sd`** (z sd = 2.633, sd-dependent) into a proper likelihood.
  Untested; it is the only channel here whose failure is a *calibration* failure rather than
  a signal failure.
- **Whether error ALIGNMENT can be engineered rather than inherited.** Every aligned emitter
  in this study is aligned because it is a real fragment. Nothing measured says a generator
  can be steered into the quiet subspace.
- **Terminal-gap-shaped coverage.** The largest lever in the phase diagram (1.79 A) and no
  channel in Part A has a coverage *pattern* that exploits it.
