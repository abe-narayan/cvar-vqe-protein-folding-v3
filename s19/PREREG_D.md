# PRE-REGISTRATION — AGENT D (theoretical / adversarial / research auditor), Sprint 19

Written **before** any result was inspected. Per `s19/BRIEF.md` §14 this file is **not edited
after results are seen**. Where a pre-registration turns out to be mis-specified, the defect is
recorded below it and the untested regime is left OPEN.

Unit of analysis is the **target** throughout (§9). Pair-level statistics are diagnostics and are
reported with a target-level aggregation beside them, never alone.

---

## Block P1 — ATTACK: is the coordinator's multimodality census real?

**Target of attack.** `s19/BRIEF.md` §3: "fraction of pairs with two or more modes (>2% mass) =
**0.217**", measured on 12 targets, offered as evidence that the moment collapse in
`s15/distcal.gather` aims the objective at values the model considers unlikely.

**Mechanism under test.** `core/predict.BIN_EDGES` is **non-uniform**: bin widths run
0.5 Å (4.5–8 Å) → 1.0 → 1.5 → 2 → 3 → 4 Å, plus two unbounded catch-all bins. A mode detector that
looks for local maxima of the **probability-mass** vector is therefore biased: under any smooth
density, mass ∝ width, so a wide bin can be a local mass maximum with no density maximum at all. A
long right tail crossing into the 3–4 Å-wide bins will manufacture a "second mode".

**Statistics computed** on all 126 cached distograms (`s12/cache/disto_*.npz`), 
so the census is at n=126 not n=12:

- `mass_modes` — strict local maxima of `p`, each ≥ 0.02 mass. (Reconstruction of the
  coordinator's definition; the exact code is his, so an exact match is not expected and a
  reconstruction is labelled as one.)
- `dens_modes` — same detector on `p / w`, `w` = bin width (catch-all bins assigned their
  nominal 1.0 / 4.0 Å).
- `prom_modes` — density modes surviving a prominence gate: the valley between two accepted
  modes must fall below 0.5 × the smaller peak, and each mode must carry ≥ 0.05 mass.
- `rebin2_a`, `rebin2_b` — density modes after merging adjacent bins pairwise at two different
  offsets (a coarsening that any real mode must survive).

**Primary outcome.** Fraction of pairs multimodal under each detector, and the target-level mean
of that fraction with a bootstrap CI.

**F-D1 (falsifier of the census).** If `dens_modes ≥ 2` on **< half** the pairs that
`mass_modes ≥ 2` flags, the 21.7% figure is substantially an artefact of the non-uniform binning
and must be withdrawn as stated. If the density-corrected and re-binned fractions all stay within
±25% of the mass-based fraction, the census stands and F-D1 does not fire.

**Note recorded in advance.** F-D1 firing would damage the *stated evidence*, not necessarily the
*hypothesis*. Blocks P2/P3/P5 test the hypothesis itself and are run regardless of P1's outcome.

---

## Block P2 — ATTACK: is multimodality HARMFUL, or merely PRESENT? (ORACLE diagnostic)

**Claim under test.** "On 21.7% of pairs the deployed objective aims at a value the model itself
considers unlikely… at a distance the geometry may never realise."

**The confound that must be removed first.** Multimodal pairs are wider and longer-range. Any
naive comparison of |d̂ − d_true| between modal classes restates "wide distributions have big
errors", which is not a mechanism. All comparisons are therefore **matched within
(separation bin × sd quintile)** cells, and the cell-matched difference is reported.

**Arms, per pair (ORACLE — the native is read only to score):**

| symbol | aim point | nature |
|---|---|---|
| `e_mean` | the deployed `expected` (first moment) | the incumbent |
| `e_med` | the distribution's median under `CENTRES` | native-free alternative |
| `e_mode1` | the highest-density mode | native-free alternative |
| `e_best` | the mode **nearest the truth** | **ORACLE ceiling of any mode-selection scheme** |
| `e_flip` | mean + ε·δ, δ = \|best mode − mean\|, ε = ±1 by `stable_rng` | **matched-magnitude random control** |

**F-D2a.** If, inside matched (sep × sd) cells, `e_mean` on multimodal pairs is **no larger** than
on unimodal pairs, then multimodality is present but not harmful at the aim point, and the
proposed mechanism loses its premise.

**F-D2b (the control that decides it).** If `e_best − e_mean` is not better than
`e_flip − e_mean` by a clear margin, then "the modes" carry no information about where the truth
is beyond the size of the move, and mode-selection is dead even with an oracle.

---

## Block P3 — ATTACK: is the multimodality COHERENT across pairs?

**Claim under test.** "…and it does so **coherently across correlated pairs**, because a whole
region flipping between two conformer families moves many pairs the same way at once." The brief
asserts this without evidence.

- **P3a, clustering.** Fraction of multimodal pairs sharing a residue with another multimodal
  pair, against a **separation-stratified label permutation** null (permuting the multimodal flag
  within separation bins, which preserves the fact that long-range pairs are both more multimodal
  and differently connected). 1000 permutations, per target, aggregated at target level.
- **P3b, signed coherence (ORACLE).** δ_pq = (nearest-mode-to-truth − mean), signed. If a region
  flips as a unit, δ must agree in sign between pairs sharing a residue more often than between
  disjoint pairs. Report sign-agreement adjacent vs non-adjacent, and the leading-eigenvalue share
  of the pair-δ structure against a sign-permutation null.

**F-D3.** If adjacency-conditioned agreement matches the stratified permutation null on both P3a
and P3b, "coherent across correlated pairs" is **REFUTED** as stated, and the distinction the
brief draws from a permuted residual (which is incoherent by construction) collapses.

---

## Block P4 — ATTACK: is the multimodality manufactured by the ENSEMBLE?

`Distogram.for_target` sets `prob = mean_m predict_proba_m(X)` over a dropout ensemble.
**A mixture of unimodal members is multimodal whenever the members disagree.** If so, the modes
are not a belief any member holds; they are a display of ensemble disagreement, and the correct
reading of the diagnosis changes.

**Measurement.** For a fold's models, per-member modality vs the ensemble-mean modality on the
same pairs; the fraction of ensemble-multimodal pairs whose members are individually unimodal;
and member-mean-spread as a predictor of ensemble multimodality.

**F-D4.** If ≥ 2/3 of ensemble-multimodal pairs have all members unimodal, the modes are an
ensembling artefact, and *that* — not the moment collapse — is the finding. Recorded as a
separate claim either way.

---

## Block P5 — THE HEADROOM TEST (decisive; run even if P1–P4 all fire)

**The question that decides whether the branch is worth a sprint arm at all:** with a *perfect*
oracle mode-picker, does aiming the deployed objective at modes instead of means lower Cα-RMSD?

**Instrument.** Exactly `s17/refine.py`'s `refine_full`: start = ideal-geometry projection of the
coordinate average of the shipped top-75; objective `Σ ((d − target)/sd)²`; optimiser
`s15/align_lib.fit`; leave-fold-out separation debias fitted on the **full** `I.targets()` per
§7's subset trap; `target = max(·, 2.0)`. Only the **target vector** changes between arms.

| arm | target vector |
|---|---|
| `raw` | deployed d̂ (debiased) — must reproduce `s17` refine_full ≈ 3.610 as a validity gate |
| `med` | distribution median, same debias offset applied |
| `mode1` | highest-density mode, same debias offset |
| `ORACLE_bestmode` | mode nearest the truth — the ceiling of all mode selection |
| `rand_flip` | mean ± δ with δ = \|bestmode − mean\| and random sign — matched-magnitude null |
| `ORACLE_true` | the true distances — the parameterisation floor, for scale |

**Primary outcome.** Paired target-level mean difference `ORACLE_bestmode − raw`, n=126,
fold-aware bootstrap CI, median and W/L beside it.

**F-D5 (the kill shot).** If `ORACLE_bestmode − raw` does not reach **−0.084 Å** (the instrument's
MDE) with a CI excluding zero, then **no mode-selection scheme of any kind can move this metric**,
and the coordinator's mechanism — whatever its truth as a description — has **no headroom** as a
lever. Report immediately.

**F-D5b.** If `ORACLE_bestmode` does beat `raw` but not `rand_flip`, the gain is the size of the
move, not the modes.

**Validity gates, declared in advance.** (i) `raw` must reproduce `s17/refine.py`'s `refine_full`
per-target to |Δ| < 0.01 Å on at least 120/126, else the arm is not the deployed functional and no
number from it is reported. (ii) `sd` is held **identical across every arm** — this sprint's
predecessor retracted a result (G5, `objceil.py` line 163) precisely because a control permuted the
weights along with the residuals. (iii) Seeds by `s15/seed.stable_rng`; no bare `hash()`.

---

## Block P6 — THEORY (stated in advance so it cannot be back-fitted)

**T1 (to be checked as an identity, not measured as a discovery — §10).** The selection path's
Bayes risk `_risk[·,x] = Σ_b p_b |x − c_b| · w` is a **sum of absolute values with non-negative
weights**, hence **convex in x**, with minimiser the weighted median of `CENTRES`. If so, the
brief's "the selection path already consumes the full distribution" is true in letter but the
selection path is *incapable of expressing multimodality*: it aims at a single robust location
estimate. This would **remove** the brief's own counter-argument ("why did it not show there?")
rather than support it, and is reported as EXACT if the algebra holds.

**T2.** Any per-pair-separable objective Σ_pq f_pq(d_pq) — which includes both the deployed
squared-error form and any full-distribution −log p form — cannot represent a *joint* constraint
that pairs choose the same conformer. Whether a non-convex per-pair f can nevertheless let the
geometry select a self-consistent set of modes is an empirical question and is exactly what P5
bounds from above.

---

## Deliverables and artefacts

`s19/agentD_FINDINGS.md`, this file, and immutable per-experiment directories under
`s19/results/` each carrying `config.json`, the result JSON, and a `COMPLETE` flag written only
after the last target. No file of any prior sprint is overwritten. The sealed 60-target benchmark
is not read, probed, or derived from; `results/benchmark_manifest.json` is not opened.
