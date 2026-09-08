# Sprint 12 — FAILURE-TARGET FORENSICS (FAIL18)

**Question.** Why do the 18 zero-recall targets fail, what do they have in common, and
what information available to a real predictor would have let the system recognise the
correct structure on them?

**Instrument validation.** `python -m s12.instrument` conventions reproduced through
`instrument.selfcheck()` paths used by every script here: pool best 1.711, top-75 best
2.306, 18 zero-recall targets identical to `I.FAIL18`. All experiments run on
tuning126 only. benchmark60 and dev24 untouched. Natives (`rr`, `nat_ca`, deposited
coordinates, PDB headers) are used **as labels / diagnostics only**; every arm that uses
them is labelled ORACLE or DIAGNOSTIC.

Code: `s12/fail_headers.py`, `s12/fail_dossier.py`, `s12/fail_contrast.py`,
`s12/fail_ensemble.py`, `s12/fail_recognise.py`, `s12/fail_oracle.py`.
Results: `s12/results/fail_*.json`.

---

## E1. Provenance audit of the deposited entries (DIAGNOSTIC — headers only)

`s12/fail_headers.py` parses HEADER/TITLE/COMPND/KEYWDS/EXPDTA/NUMMDL/SSBOND/LINK/HETATM
of `pdbs_ext/<PDB>.pdb` for all 126 tuning targets and flags provenance classes.

**Hypothesis.** The 18 are enriched in entries whose deposited conformation is not
determined by the isolated sequence (fibril/steric-zipper segments, partner-bound
fragments, micelle-bound peptides) or whose backbone topology is not linear
(lasso peptides, head-to-tail cyclics).

| flag | FAIL18 | other-108 | Fisher p |
|---|---|---|---|
| **fibril / amyloid / steric zipper** | **6/18** | 4/108 | 5.5e-04 |
| **lasso peptide** | **4/18** | 2/108 | 3.8e-03 |
| **fibril OR lasso** | **10/18 (56%)** | 6/108 (5.6%) | **1.2e-06** |
| cyclic | 1/18 | 6/108 | 1.00 |
| membrane / micelle / TFE | 5/18 | 45/108 | 0.31 (n.s., *depleted*) |
| partner-bound in title | 2/18 | 9/108 | 0.66 |
| X-ray / MicroED | 3/18 | 6/108 | 0.12 |
| SSBOND/LINK/cyclic ("constrained topology") | 8/18 | 52/108 | 0.80 |

Negative controls matter here: **membrane/micelle peptides are NOT enriched** (they are
the single largest class in the database and are *depleted* among the failures), nor are
disulfide/LINK-constrained peptides, nor "FRAGMENT:" entries. The enrichment is specific
to two structural classes.

### The lasso class is a complete failure class
There are exactly **6 lasso peptides** in tuning126. Four are in FAIL18; the other two are
not zero-recall but are still failures by RMSD:

| pdb | in FAIL18 | emitted CA-RMSD | pool best | distogram MAE on native | native score pct |
|---|---|---|---|---|---|
| 2N5C | yes | 6.39 | 2.46 | 3.38 | 0.55 |
| 7JS6 | yes | 7.18 | 3.13 | 5.36 | 0.95 |
| 7LCW | yes | 6.61 | 3.19 | 4.20 | 0.86 |
| 9KAR | yes | 7.54 | 1.40 | 4.95 | 0.88 |
| 2MAI | no | 6.41 | 3.04 | 3.58 | 0.76 |
| 2MFV | no | 6.10 | 3.95 | 3.72 | 0.81 |

**6/6 lasso peptides emit ≥ 6.1 Å.** Mean 6.71 Å against a 3.236 Å system mean. A lasso
peptide's backbone is threaded through an 8–9-residue isopeptide macrolactam ring; the
conformation is enforced by a covalent constraint that (i) no linear window in the library
possesses and (ii) the sequence-only distogram has never seen. This is not a ranking
failure — it is a **representation/coverage** failure.

### The fibril class splits by length
| pdb | n | emitted | pool best | note |
|---|---|---|---|---|
| 7N2I | 9 | 1.28 | 0.74 | MicroED single β-strand — works |
| 5V5B | 10 | 3.44 | 1.07 | works-ish |
| 7VI4 | 10 | 2.85 | 0.88 | works-ish |
| 1S9Z | 16 | 0.32 | 0.16 | works (β-hairpin-like) |
| 3SGO | 11 | 6.28 | 1.60 | FAIL |
| 5W52 | 11 | 5.67 | 1.12 | FAIL |
| 2BFI | 12 | 7.05 | 2.22 | FAIL |
| 9L1M | 12 | 5.00 | 2.96 | FAIL |
| 2JN5 | 12 | 5.57 | 2.69 | FAIL (also partner-bound) |
| 2BP4 | 16 | 5.41 | 0.63 | FAIL (also TFE cosolvent) |

The failing amyloid segments have pool bests of 0.63–2.96 Å: **the right conformation is
in the pool**, the objective refuses to pick it. See E2.

---

## E2. Length- and fold-matched control group (18 vs 18)

`s12/fail_contrast.py` builds a control by matching each FAIL18 target to the
best-emitting unused other-108 target of the *same length* (fold-preferred). Means are
identical in length (14.000 vs 14.000), so nothing below is a length artefact.

| quantity (ORACLE where `o_`) | FAIL18 | MATCH18 | other-108 | diff (F−M) 95% CI |
|---|---|---|---|---|
| emitted CA-RMSD | 6.056 | 2.300 | 2.765 | +3.756 [3.03, 4.52] |
| pool best (K=500) | 2.284 | 1.449 | 1.615 | +0.835 [0.26, 1.41] |
| universe best | 1.640 | 1.120 | 1.259 | +0.520 [0.05, 0.99] |
| top-75 best | 4.677 | 1.639 | 1.911 | +3.039 [2.42, 3.67] |
| near-native band size in pool | 43.4 | 167.0 | 149.2 | −123.6 [−159, −86] |
| **band's score percentile** | **0.687** | 0.251 | 0.287 | +0.436 [0.34, 0.53] |
| **native's own score percentile** | **0.784** | 0.321 | 0.298 | +0.463 [0.33, 0.59] |
| ρ(score, true RMSD) in pool | 0.107 | 0.734 | 0.645 | −0.627 [−0.80, −0.45] |
| ρ(score, true RMSD) in band | −0.219 | 0.245 | 0.184 | −0.464 [−0.68, −0.27] |
| **distogram MAE on the native** | **4.367** | 1.764 | 2.001 | +2.603 [1.81, 3.47] |
| distogram r vs true CA dists | 0.445 | 0.801 | 0.738 | −0.356 [−0.54, −0.17] |
| distogram calibration slope | 0.436 | 0.874 | 0.775 | −0.438 [−0.63, −0.23] |
| MAE, separation 2–4 | 1.885 | 1.113 | 1.208 | +0.772 [0.33, 1.25] |
| MAE, separation 5–8 | 5.130 | 2.033 | 2.490 | +3.097 [1.92, 4.39] |
| bias (pred−true), sep 2–4 | −0.474 | +0.515 | +0.274 | −0.989 [−1.87, −0.12] |
| native helix fraction | 0.133 | 0.523 | 0.346 | −0.390 [−0.60, −0.16] |
| native strand fraction | 0.109 | 0.046 | 0.036 | +0.063 [−0.02, 0.14] |
| native coil fraction | 0.758 | 0.431 | 0.619 | +0.326 [0.11, 0.52] |
| pool fraction from peptide DB | 0.253 | 0.294 | 0.277 | −0.042 [−0.09, 0.01] |
| top-75 fraction from peptide DB | 0.177 | 0.280 | 0.276 | −0.103 [−0.21, −0.00] |
| band fraction from peptide DB | 0.389 | 0.385 | 0.374 | +0.004 (n.s.) |
| distinct window sequences in pool | 385 | 380 | 366 | n.s. |

Reading of the table:

1. **The pool is not the problem.** Pool best 2.284 Å; on 15/18 it is under 3.2 Å, on 6/18
   under 1.7 Å. Retrieval diversity (distinct sequences, peptide fraction) is
   indistinguishable from the control.
2. **The objective is the problem.** The native itself sits at the **78th percentile** of
   the shipped distogram score over its own pool (32nd on matched controls): the score
   actively prefers ~4 out of 5 pool members over the true answer. ρ(score, RMSD) is 0.107
   vs 0.734, and *negative* inside the band.
3. **The distogram is 2.5× worse on these natives, mostly at medium separation** (5–8
   residues: MAE 5.13 vs 2.03). Short-range it is *too compact* (bias −0.47 vs +0.52),
   i.e. it predicts a turn/helix where the native runs straight.
4. The band is thin (43 vs 167 members): a scorer with the control's in-band skill would
   still have ~4× less mass to hit. But that is a consequence of the native being far from
   the pool mode, not a separate cause.
5. The top-75 under-selects peptide-database windows (0.177 vs 0.276) while the band is
   *enriched* in them (0.389): the distogram is systematically steering the filter towards
   protein fragments on exactly the targets where real peptides hold the answer.

### Failure taxonomy (ORACLE diagnostic; same classifier run on MATCH18 as the null)

| class | rule | FAIL18 | MATCH18 (null) |
|---|---|---|---|
| (c) objective | native score pct > 0.5 **or** ρ(score,RMSD) < 0.25 | **17/18** | 7/18 |
| (a) query | pool best − universe best > 0.5 Å | 12/18 | 6/18 |
| (d) label/context | fibril, lasso, cyclic or partner-bound header | 12/18 | 3/18 |
| (b) filter only | band in pool, excluded, without (c) | **1/18** (1JBF) | 0/18 |

**The dominant mechanism is (c) objective failure**, present on 17/18 and at 2.4× the
control rate. Pure filter failure — the objective is fine but the top-75 cut misses — is
one target. Query loss is real but secondary (0.52 Å of the 0.84 Å pool-best gap).

---

## E3. NMR ensemble analysis — the label-convention hypothesis is REFUTED

`s12/fail_ensemble.py` parses every deposited model of all 126 entries
(`protein_geometry.native_ensemble_from_pdb`). Model-1 in the cached universes is
verified identical to model 1 on disk (max discrepancy 0.0000 Å over 126 targets).
The project's declared endpoint (CA-RMSD to model 1) is reported first and unchanged;
the ensemble-aware numbers are a secondary diagnostic.

| quantity | FAIL18 | MATCH18 | other-108 |
|---|---|---|---|
| deposited models | 14.7 | 13.3 | 15.8 |
| **ensemble spread** (mean pairwise CA-RMSD) | **0.942** | 0.946 | 0.955 |
| model-1 → ensemble medoid | 0.722 | 0.659 | 0.596 |
| model-1 centrality percentile | 0.276 | 0.346 | 0.365 |
| emitted vs model-1 (**declared metric**) | **6.056** | 2.300 | 2.765 |
| emitted vs best deposited model (diagnostic) | 5.878 | 2.083 | 2.453 |
| pool best vs model-1 | 2.284 | 1.449 | 1.615 |
| pool best vs best deposited model | 2.046 | 1.238 | 1.383 |

**The 18 are not floppier than the controls.** Ensemble spread is 0.942 vs 0.946 Å —
indistinguishable. Switching from the model-1 convention to a best-of-ensemble convention
recovers 0.178 Å on FAIL18 and 0.217 Å on MATCH18; the paired difference is
−0.039 Å [−0.194, +0.113], 10W/7L — **null**. Model-1 is not a systematically worse
representative on the failures.

Three individual exceptions worth naming (they do not carry the group effect):

| pdb | spread | m1→medoid | note |
|---|---|---|---|
| 3BTB | **4.265** | **4.294** | band-3 peptide in fast exchange with GAPDH; the deposited ensemble is genuinely scattered and model-1 is not representative. Genuine (d) label problem. |
| 8T63 | 1.555 | 2.551 | model-1 at the 90th percentile of non-centrality |
| 1ID6 | 1.669 | 1.791 | model-1 at the 91st percentile of non-centrality |

**Conclusion: (d) accounts for ~1–3 of 18, not the mechanism.** Three of the 18 (2BFI,
3SGO, 5W52) have a single deposited model (crystal/MicroED), so ensemble spread is not
even defined for them; their failure has a different cause (E4).

---

## E4. The recognition question — every conventional signal is ANTI-correlated on the failures

`s12/fail_recognise.py`. Inside each K=500 pool, AUC for separating the near-native band
(ORACLE label, evaluation only) from the rest, using only deployable signals. 0.5 = no
skill; **< 0.5 = the signal points the wrong way.**

| signal | FAIL18 AUC | MATCH18 AUC | FAIL18 drop-top-4 | FAIL18 rg-matched | z(band − top-75) | r(AUC, pool best) |
|---|---|---|---|---|---|---|
| shipped distogram score | **0.301** | 0.869 | 0.249 | 0.308 | −1.74 | −0.46 |
| ... restricted to sep 2–4 | 0.522 | 0.769 | 0.412 | 0.472 | −0.88 | −0.33 |
| ... restricted to sep ≥ 5 | **0.273** | 0.876 | 0.195 | 0.258 | −1.70 | −0.54 |
| consensus / centrality | **0.252** | 0.793 | 0.138 | 0.224 | −1.31 | −0.59 |
| **ESM-2 contact agreement** | **0.593** | 0.631 | 0.495 | 0.619 | **+0.16** | **+0.055** |
| Legacy 11-term energy | 0.380 | 0.743 | 0.285 | 0.352 | −0.75 | −0.53 |
| rg vs distogram-implied rg | 0.293 | 0.687 | 0.205 | 0.393 | −1.60 | −0.36 |
| SS vs propensity SS guess | 0.332 | 0.753 | 0.232 | 0.307 | −0.77 | −0.65 |

Two conclusions, both important.

1. **On FAIL18 the near-native band is the pool's OUTLIER set.** Consensus AUC 0.252
   means the band members are 3:1 *less* central than random pool members. This is a
   direct, per-target confirmation of the standing memory result that consensus is
   outlier-*avoidance* rather than a nativeness signal: on the 108 the native happens to
   sit near the pool mode and centrality reads as skill (AUC 0.793); on the 18 it does
   not and the same signal inverts. Every centrality-flavoured method (medoid, typicality,
   the coordinate average itself) is structurally guaranteed to fail on this subgroup.
   Null control N1 confirms it: AUC correlates −0.59 with pool best, i.e. these signals
   are largely reporting difficulty.

2. **ESM-2 contact agreement is the only signal that survives.** It is the only one with
   AUC > 0.5 on FAIL18, the only one with a positive z-gap (the band scores *higher* than
   what the shipped filter actually picks), the only one that survives rg-matching
   (0.619 — it gets *stronger*), and the only one that passes N1 (r = +0.055 with
   difficulty: it is not just describing easy targets). Note it is the *weakest* signal on
   the controls (0.631 vs 0.869 for the distogram) — this is a genuinely decorrelated
   channel, strongest exactly where the shipped one dies.
   Concentration caveat: dropping the 4 best targets takes it to 0.495, so the FAIL18 mean
   rests on ~8 targets (1JBF 0.98, 2NDM 0.97, 2N5C 0.95, 7LCW 0.84, 3SGO 0.83, 5W52 0.81,
   8T63 0.80, 2MQ2 0.77) and it is inverted on four (2BP4 0.02, 9KAR 0.13, 2BFI 0.21,
   2NB7 0.25).

### The unifying quantitative defect: the compactness channel is broken, in both directions

The distogram's own expected distances imply a radius of gyration
(rg² = Σd²ᵢⱼ/2n², with |i−j|=1 filled at 3.81 Å). Against the native rg:

| | bias | MAE |
|---|---|---|
| MATCH18 | +0.519 | **0.824** |
| FAIL18 | +0.255 | **2.440** |

The *bias* is the same; the *magnitude* is 3× worse, because FAIL18 splits into two
opposite classes that cancel in the mean:

| class | targets | rg_pred | rg_native | error |
|---|---|---|---|---|
| **too compact predicted, native EXTENDED** | 2BFI, 3SGO, 5W52, 8T63, 2MLQ-like | 5.75–7.36 | 10.0–12.1 | −3.4 to −4.7 |
| **too extended predicted, native COMPACT/knotted** | 2MQ2, 7JS6, 9KAR, 7LCW, 2NB7, 2N5C | 7.1–10.1 | 5.1–6.3 | +1.8 to +4.7 |

The extended class is exactly the amyloid steric-zipper set (a fibril strand is a
straight β-ribbon, ~12 Å rg at n=11–12, a conformation no globular fragment library and
no compactness-regularised objective will ever prefer). The compact class is exactly the
lasso peptides plus the disulfide-free protegrin hairpin (a knotted macrolactam is far
more compact than a linear 15-mer). On the extended targets the rg-fit AUC is 0.002–0.023
— literally every band member is on the wrong side of the objective's compactness
preference.

---

## E5. Can the ESM contact channel be USED? — a positive that dies to its own control

`s12/fail_esmrescore.py` (all 126, no projection). Filter arm:
`rescore = z(−shipped score) + w·z(ESM-2 contact agreement)`, take top-75. Evaluated by
`top75_best` (best true CA-RMSD inside the selected 75 — what synthesis can reach) and
band recall. `w = 0` reproduces the shipped filter exactly. **Null control**: the same
statistic computed against a *different, length-matched target's* ESM contact map.

| w | top75_best all-126 | FAIL18 | other-108 | recall FAIL18 | recall other | **null all-126** | **null FAIL18** |
|---|---|---|---|---|---|---|---|
| 0.00 (shipped) | 2.306 | 4.677 | 1.911 | 0.00 | 1.00 | 2.306 | 4.677 |
| 0.25 | **2.181** | 4.038 | 1.871 | 0.44 | 1.00 | 2.247 | 4.196 |
| 0.50 | 2.186 | 3.858 | 1.907 | 0.44 | 0.99 | 2.229 | 3.892 |
| 0.75 | 2.198 | 3.679 | 1.951 | 0.56 | 0.97 | 2.308 | 3.768 |
| 1.00 | 2.219 | 3.614 | 1.986 | 0.56 | 0.96 | 2.394 | 3.644 |
| 3.00 | 2.317 | 3.560 | 1.910 | 0.61 | 0.90 | 2.665 | **3.519** |

Paired on all 126 at w = 0.25: −0.126 Å [−0.206, −0.052], 37W/20L — but
**drop-top-10 = −0.024 and drop-top-20 = +0.016**, and the *null* map alone gives
−0.060 [−0.152, +0.034]. On FAIL18 the null tracks the real map at every w and **beats it
at w = 3**. The apparent ESM gain is not sequence-specific contact information.

---

## E6. What the ESM arm was actually doing: the top-75 is a degenerate set

`s12/fail_diversity.py` (all 126). Arms all take 75 of the 500; w = 0.75.

**Redundancy of the shipped top-75** (mean pairwise CA-RMSD inside the selected set):

| | shipped top-75 | random 75 from the same pool |
|---|---|---|
| FAIL18 | **2.879** | 4.608 |
| other-108 | 2.421 | 4.210 |

The filter collapses onto ~one basin. On the 108 that basin contains the native; on the
18 it does not, and there is no second opinion in the set for the coordinate average to
recover from.

| arm | top75_best all-126 | FAIL18 | other-108 | **recall FAIL18** | Δ all-126 vs shipped | 95% CI | drop-top-10 |
|---|---|---|---|---|---|---|---|
| shipped | 2.306 | 4.677 | 1.911 | **0/18** | 0 | — | — |
| esm (own contacts) | 2.198 | 3.679 | 1.951 | 10/18 | −0.108 | [−0.231, +0.010] | +0.032 |
| esmnull (foreign contacts) | 2.308 | 3.768 | 2.065 | 7/18 | +0.002 | [−0.134, +0.138] | +0.151 |
| ncontact (contact count only) | 2.230 | 3.904 | 1.951 | 7/18 | −0.076 | [−0.199, +0.037] | +0.071 |
| **random N(0,1)** | 2.120 | 3.764 | 1.846 | 8/18 | **−0.186** | [−0.289, −0.094] | −0.068 |
| **divmax (greedy max-min diversity)** | **2.094** | **3.306** | 1.892 | **14/18** | **−0.212** | [−0.330, −0.102] | −0.057 |

In-band AUC against its own null:

| signal | FAIL18 | other-108 |
|---|---|---|
| shipped score | 0.301 | 0.812 |
| ESM contact agreement | 0.593 | 0.629 |
| **foreign-map agreement (null)** | **0.574** | 0.480 |
| contact count | 0.541 | 0.643 |

On the 108, ESM content is real (0.629 vs a 0.480 null). **On FAIL18 it is not**
(0.593 vs 0.574): the surviving skill in E4 was almost entirely the *shape* of the
statistic — how many long-range contacts the candidate has — not what ESM says about
*this* sequence.

**The corrected reading.** Zero recall on FAIL18 is not primarily "no signal can see the
band". It is that the shipped filter selects 75 near-duplicates of one wrong basin.
Adding *pure Gaussian noise* to the score recovers 8/18 recall and beats the ESM arm on
all 126 (−0.186 vs −0.108, and unlike ESM it survives drop-top-10). Explicit diversity
does better still: 14/18 recall, top-75 best 4.677 → 3.306 Å on FAIL18 and 2.306 → 2.094 Å
on all 126, drop-top-10 −0.057.

**Necessary caveat, stated before the emitted test.** `top75_best` is a *minimum* over
the selected set and a more dispersed set lowers a minimum mechanically. Diversity buys
reachability and *costs* set quality: mean true RMSD of the selected 75 goes 3.551 → 4.281
on all 126 (3.116 → 4.006 on the other-108; on FAIL18 it *improves*, 6.159 → 5.932).
Since the coordinate average tracks the set mean, the emitted-structure test in E8 is the
one that decides whether this is a real lever or a bookkeeping artefact.

---

## E7. (coordinator request) Which FORM of ESM contact agreement carries content?

`s12/fail_contactforms.py`, all 126. Seven forms of candidate-vs-ESM contact agreement,
each scored **twice**: against the target's own ESM-2 contact map and against a
length-matched **foreign** target's map. `net = own − null` is the part of the signal that
is about *this* sequence. In-band AUC, judged before any emitted arm as asked.

| form | FAIL18 own | FAIL18 null | **FAIL18 net** | 95% CI | other-108 own | other-108 null | other-108 net |
|---|---|---|---|---|---|---|---|
| corr3 (the E4 form) | 0.593 | 0.567 | +0.026 | [−0.076, +0.157] | 0.634 | 0.485 | **+0.149** |
| corr5 (|i−j| ≥ 5) | 0.567 | 0.527 | +0.040 | [−0.051, +0.164] | 0.580 | 0.523 | +0.057 |
| cb5 (CB–CB) | 0.548 | 0.515 | +0.033 | [−0.061, +0.161] | 0.584 | 0.524 | +0.060 |
| soft5 (sigmoid) | 0.527 | 0.530 | −0.003 | [−0.119, +0.135] | 0.582 | 0.502 | +0.080 |
| loglik5 | 0.348 | 0.374 | −0.025 | [−0.074, +0.011] | 0.594 | 0.543 | +0.052 |
| rank5 (Spearman) | 0.469 | 0.516 | −0.047 | [−0.170, +0.107] | 0.632 | 0.516 | **+0.116** |
| **topk** (top-n ESM pairs) | **0.676** | 0.634 | +0.043 | [+0.006, +0.095] | 0.512 | 0.481 | +0.030 |

`topk` has the **highest raw in-band skill of any native-free signal measured anywhere in
this study (0.676 on FAIL18)** — but 0.634 of it is reproduced by a foreign contact map.
Only one of seven forms has a net CI excluding zero, marginally, on n = 18, after seven
comparisons. **On the other-108 ESM content is real (corr3 net +0.149); on FAIL18 it is
not.** What ESM contributes on the failures is the generic shape of the statistic —
"does this candidate have long-range contacts where any ESM map puts its mass" — i.e. a
compactness/contact-order axis, not sequence-specific structure.

## E8. (coordinator request) Fusion sweep and the routing question

`s12/fail_fusion.py`, all 126, `rescore = z(−score) + w·z(signal)`, top-75.

| arm | top75_best all-126 | FAIL18 | other-108 | recall FAIL18 | Δ vs shipped | 95% CI | **drop-top-10** |
|---|---|---|---|---|---|---|---|
| shipped (w=0) | 2.306 | 4.677 | 1.911 | 0/18 | 0 | — | — |
| topk @0.25 | 2.195 | 3.950 | 1.902 | 9/18 | −0.111 | [−0.198, −0.034] | **−0.000** |
| topk @0.5 | 2.169 | 3.519 | 1.944 | 12/18 | −0.138 | [−0.265, −0.016] | **+0.033** |
| topk-NULL @0.25 | 2.173 | 3.887 | 1.888 | 7/18 | −0.133 | [−0.226, −0.048] | −0.015 |
| topk-NULL @0.5 | 2.199 | 3.572 | 1.970 | 10/18 | −0.107 | [−0.232, +0.015] | +0.053 |
| corr3 @0.5 | 2.173 | 3.852 | 1.893 | 9/18 | −0.133 | [−0.237, −0.037] | −0.008 |
| corr3-NULL @0.5 | 2.224 | 3.904 | 1.943 | 6/18 | −0.083 | [−0.195, +0.028] | +0.047 |
| **divmax (diversity only, no ESM)** | **2.094** | **3.306** | 1.892 | **14/18** | **−0.212** | [−0.330, −0.102] | **−0.057** |

Every ESM arm is matched to within noise by its own foreign-map null, and every ESM arm's
effect **vanishes or reverses under drop-top-10** — rule 8 of the brief says ignore it.
The diversity arm, which uses no ESM at all, is larger, has a CI excluding zero, keeps
14/18 recall, and is the only arm that survives drop-top-10.

### The router is null
LFO logistic regression (5 pinned folds) on 13 native-free features — distogram entropy
and mean sd, distogram-implied rg, pool/top-75 median rg and their gaps, top-75 spread,
ESM contact sharpness and mass, score gap and sd, n:

| | value |
|---|---|
| LFO AUC, FAIL-regime vs rest | **0.558** |
| label-permutation null (200 draws) | 0.511 [0.314, 0.680] |
| p | **0.295** |
| precision@18 | 0.278 (base rate 0.143) |

Strongest single feature: distogram mean sd (r = +0.245), then distogram-implied rg
(+0.241) and n (+0.198) — nothing usable. **There is no native-free router for the FAIL
regime**, so a regime-dependent blend is not available and a fixed setting is the honest
option. This is the third independent confirmation in the record that weak routers fail
here.

## E9. (coordinator request) Are the lasso and fibril classes reachable in principle?

`s12/fail_reach.py`: how much of the *whole leakage-safe universe* (not just the K=500
pool) lies within 2.0 Å of the native, and where BLOSUM ranks it.

| class | n | universe best | pool best | emitted | windows ≤ 2.0 Å | % of universe | of which in K=500 | BLOSUM rank of the universe best (median) |
|---|---|---|---|---|---|---|---|---|
| all 126 | 126 | 1.31 | 1.71 | 3.24 | 1171 | 6.01% | 54.5 | 6176 |
| other-108 | 108 | 1.26 | 1.62 | 2.77 | 1342 | 6.79% | 62.3 | 6439 |
| FAIL18 | 18 | 1.64 | 2.28 | 6.06 | 146 | 1.34% | 7.8 | 5492 |
| **lasso (all 6)** | 6 | **2.01** | 2.86 | 6.70 | **2.3** | **0.02%** | **0.3** | 3629 |
| **fibril ∈ FAIL18** | 6 | **1.20** | 1.87 | 5.83 | **371** | **3.51%** | 16.7 | **12175** |
| fibril ∉ FAIL18 | 4 | 0.55 | 0.71 | 1.97 | 1565 | 7.87% | 78.2 | 5023 |
| FAIL18, neither | 8 | 1.91 | 2.47 | 5.79 | 49 | 0.36% | 4.8 | 5058 |

**Verdict on reachability.** The two classes are *opposite* problems and must not be
lumped together.

*Lasso peptides are not reachable by this architecture.* Of ~10,000 candidate windows per
target, a mean of **2.3** lie within 2.0 Å of the native; three of the six (2MAI, 2MFV,
7LCW) have **zero**, and the best single window in the entire universe averages 2.01 Å —
worse than the *pool* best of an average target. A lasso backbone is threaded through an
8–9-residue isopeptide macrolactam; that topology is a covalent constraint, and a library
of linear windows cut from linear chains does not contain it at any depth. The library-
saturation result (3.2× more fragments moves universe best by 0.049 Å) says more of the
same corpus will not fix it. This class should be **declared out of scope** for a
retrieval architecture, or handled by a constrained builder that is given the macrolactam
as an explicit bond — not by better ranking. It costs the tuning mean ≈ 6 targets × ~3.5 Å
excess / 126 ≈ 0.17 Å.

*Fibril/amyloid segments ARE reachable and are being lost by the query and the objective.*
371 windows per target (3.5% of the universe) are within 2.0 Å and the universe best is
1.20 Å — the extended β-conformation is abundant in a fragment library cut from folded
proteins. But BLOSUM ranks the best one at position **12,175** (vs 5,023 for the fibril
targets that succeed), so only ~17 of those 371 enter K=500, and the compactness-preferring
objective then discards them (rg-fit AUC 0.002–0.023, E4). This is a **retrieval-key**
problem, and S6-6 left exactly this hypothesis open: SS as a retrieval key (changing which
windows enter the pool) has never been tested, only SS as a post-hoc filter.

## E10. SS (and rg) as a RETRIEVAL KEY — the open S6-6 hypothesis, now closed

`s12/fail_retrieve.py`, all 126. Each key re-ranks the **whole universe**
(7k–40k windows), takes K=500, then the UNCHANGED shipped distogram filter selects 75.

| key | pool best 126 | Δ | 95% CI | pool best FAIL18 | pool best fibril | top-75 best 126 | Δ | top-75 best FAIL18 | recall FAIL18 | drop-10 |
|---|---|---|---|---|---|---|---|---|---|---|
| blosum (shipped) | 1.711 | 0 | — | 2.284 | 1.407 | 2.306 | 0 | 4.677 | 0.00 | — |
| shortD @0.5 (distogram i,i+2..4) | 1.684 | −0.027 | [−0.083, +0.040] | 2.465 | 1.714 | 2.370 | +0.064 | 4.770 | 0.11 | +0.110 |
| ssdisto @0.5 (SS read off the distogram) | 1.708 | −0.003 | [−0.056, +0.063] | 2.568 | 1.734 | 2.356 | +0.050 | 4.678 | 0.22 | +0.096 |
| ssprop @0.5 (propensity SS) | 1.817 | +0.107 | [+0.049, +0.172] | 2.757 | 1.918 | 2.294 | −0.012 | 4.535 | 0.33 | +0.065 |
| **ORACLE_ss @1.0** | 1.607 | **−0.104** | [−0.167, −0.044] | 2.128 | 1.528 | 2.196 | −0.110 | 4.114 | 0.33 | −0.005 |
| **ORACLE_rg @2.0** | 1.533 | **−0.177** | [−0.229, −0.128] | 1.866 | **1.062** | 2.099 | −0.207 | 3.677 | 0.33 | −0.072 |

**Negative, and it closes an open question.** No deployable SS or shape retrieval key
improves the pool: all three make the FAIL18 pool *worse* (2.284 → 2.47–2.90), because
re-ranking on a wrong SS/shape guess costs more sequence signal than it buys. More
importantly, **even a perfect native SS string used as a retrieval key is worth only
−0.104 Å of pool best and −0.110 Å of top-75 best** — comparable to the 0.21 Å S6-6
measured for oracle SS as a post-hoc filter, and nowhere near enough to matter. S6-6's
open hypothesis ("SS as a retrieval key is a different hypothesis and is open") should be
recorded as **closed and negative**. Oracle rg is the better key (−0.177 Å pool best, and
it does halve the fibril class's pool best from 1.407 to 1.062) but still buys only
0.207 Å of top-75 best. **The query is not where the missing information lives.**

---

## E11. The emitted-structure test of the diversity lever — it does NOT survive

`s12/fail_emit.py`: FAIL18 + MATCH18, each arm's 75 (or 25) through
`I.coordinate_average` → `I.project(λ=0.3)`. All arms deployable.

| arm | FAIL18 | Δ | W/L | MATCH18 | Δ | W/L | mean of the 36 |
|---|---|---|---|---|---|---|---|
| shipped75 (control) | 6.019 | 0 | — | 2.259 | 0 | — | 4.139 |
| score25 (cardinality control) | 6.010 | −0.010 | 10/8 | 2.298 | +0.040 | 8/10 | 4.154 |
| **divmax75** | **5.566** | **−0.453** | **15/3** | 3.103 | **+0.844** | 3/15 | 4.334 |
| divmax25 | 5.565 | −0.455 | 13/5 | 3.242 | +0.983 | 4/14 | 4.403 |
| random75 | 5.863 | −0.156 | 14/4 | 2.375 | +0.116 | 7/11 | 4.119 |

**The reachability gain does not reach the emitted structure.** divmax turns 0/18 recall
into 14/18 and drops the top-75 best from 4.677 to 3.306 Å, yet the emitted RMSD only
improves 0.45 Å on FAIL18 and gets 0.84 Å *worse* on the matched controls. The coordinate
average tracks the set *mean*, and a diversified set has a worse mean — so the operator
throws away exactly the reachability that diversity bought. This is consistent with C3 and
with the aggregation agent's result that the terminal operator is already at its
native-free optimum. **Diversifying the filter is only useful if paired with an in-band
selector that can pick out of the diversified set — and the record says no such selector
exists.** Reported as a negative for the deployable pipeline and as a *decomposition*: the
filter can be made to contain the answer; nothing native-free can then find it.

---

## E12. ORACLE counterfactuals — ranking the missing channels in emitted Å

`s12/fail_oracle.py`. Every arm gives the system exactly ONE piece of native information,
then runs the real path (subset → top-75 by the shipped score → coordinate average →
project λ=0.3). **All ORACLE / DIAGNOSTIC; none is deployable.**

| channel supplied (ORACLE) | FAIL18 emitted | **Δ vs base** | 95% CI | W/L | MATCH18 Δ |
|---|---|---|---|---|---|
| base (shipped top-75) | 6.019 | 0 | — | — | 0 |
| native secondary-structure string | 5.232 | **−0.788** | [−1.436, −0.218] | 14/4 | +0.015 |
| native radius of gyration | 5.179 | **−0.841** | [−1.283, −0.442] | 14/4 | −0.008 |
| native 8 Å contact map | 5.186 | **−0.833** | [−1.426, −0.337] | 14/4 | +0.058 |
| **true distance for the 10 most-wrong pairs** | 4.527 | **−1.493** | [−2.135, −0.890] | 16/2 | −0.265 |
| which 75 pool members are best | 3.331 | −2.689 | [−3.186, −2.206] | 18/0 | −0.481 |
| the single best pool member | **2.285** | −3.735 | [−4.244, −3.208] | 18/0 | −0.814 |

Ranked missing channels (FAIL18, emitted Å):
**selection (3.74) ≫ 10 correct distances (1.49) > rg (0.84) ≈ contact map (0.83) ≈ SS (0.79)**.

Three readings:

1. **The projection path is essentially lossless.** ORACLE single-best emits 2.285 Å
   against a pool best of 2.284 Å. Every one of the 3.73 Å of FAIL18 error is created
   between the pool and the selection — none of it by synthesis or projection.
2. **Ten pair distances out of ~85 buy 1.49 Å.** The distogram does not need to be
   globally better; it needs its *worst ten pairs* fixed. Combined with the E2 shell
   decomposition (MAE 5.13 at separation 5–8 vs 2.03 on controls) and the standing
   error-structure clue (native + random-sign residuals beats the real distogram at equal
   MAE), the missing channel is **the sign/location structure of the distogram's medium-
   separation errors**, not its average accuracy.
3. **Every "one global property" oracle is worth ~0.8 Å and only on the failures** (SS,
   rg and the contact map buy 0.00–0.06 Å on the matched controls). They are a failure-mode
   repair, not a general improvement — and none of the three has a deployable proxy that
   works (E4: propensity SS AUC 0.332; distogram-implied rg MAE 2.44 Å on FAIL18; ESM
   contacts net ≈ 0 over their null on FAIL18).

Per-target highlights: on 2BP4 the ORACLE SS alone takes 5.28 → 0.95 Å and the ORACLE
contact map 5.28 → 0.84 Å; on 2NB7 SS gives 4.56 → 1.97. On the lasso peptides no global
oracle helps much (9KAR 7.43 → 6.28/6.36/6.56; 7JS6 6.97 → 5.40 at best) because the pool
does not contain the answer (E9) — only the selection oracles move them, and only to their
pool best.

---

## E13. Leakage audit

- No benchmark60 file, `results/benchmark_manifest.json` or `s9/final_cache/*` was opened
  by any script here; all pools come from `s8/generate_univ/`. dev24 not touched.
- `rr` / `nat_ca` / deposited PDB coordinates / PDB headers appear only in (a) evaluation
  metrics, (b) the failure taxonomy, (c) arms explicitly named `o_*` / `ORACLE_*` /
  DIAGNOSTIC. No arm labelled deployable reads them: `divmax`, `random`, `esm`, `topk`,
  `corr3`, `shortD`, `ssprop`, `ssdisto` and the router use only window coordinates and
  torsions, the shipped LFO distogram, ESM-2 features of the target sequence, and Legacy
  energies.
- The router is trained leave-fold-out on the 5 pinned folds with a label-permutation null.
- The ESM contact bank is a subset (`s12/cache/esm_con_targets.npz`, 126 target sequences)
  of the shipped `esm_bank`; `esm_cache.npz` was never loaded. Peak process RSS stayed
  under ~0.7 GB; `free_gb()` never dropped below 1.1 GB.
- Selection-bias disclosure: the mixing weights in E5/E8 and the diversity λ = 1.0 in E6
  were chosen by looking at all-126 `top75_best`, i.e. **in sample**. That is a reason to
  discount those effect sizes, and it makes the E11 emitted negative *more* trustworthy,
  not less. Nothing here is proposed for a dev24 pass.

---

## E14. Anatomy of the pair channel — how few distances, and is it sign or magnitude?

`s12/fail_pairs.py`, FAIL18 + MATCH18, emitted through the real path. All ORACLE.

| arm (ORACLE) | FAIL18 emitted | Δ | W/L | MATCH18 Δ | top-75 best FAIL18 |
|---|---|---|---|---|---|
| base | 6.019 | 0 | — | 0 | 4.677 |
| true distance for the **3** worst pairs | 5.368 | −0.651 | 16/2 | −0.124 | 3.652 |
| ... 5 worst | 5.077 | −0.942 | 16/2 | −0.153 | 3.527 |
| ... 10 worst | 4.527 | −1.493 | 16/2 | −0.265 | 2.755 |
| ... 20 worst | 3.832 | −2.187 | 17/1 | −0.327 | 2.408 |
| **... ALL pairs (a perfect distogram)** | **3.640** | −2.379 | 18/0 | −0.318 | 2.287 |
| **o_sign**: every pair nudged 2.0 Å in the correct *direction* only | 4.803 | −1.217 | **18/0** | −0.199 | 3.223 |
| o_magonly: correct *magnitude*, random sign | 4.793 | −1.226 | 15/3 | −0.006 | 2.773 |

Three findings.

1. **Three distances are worth 0.65 Å; ten are worth 1.49 Å.** The value is extremely
   front-loaded — a handful of specific pairs, mean separation 9.1 (vs 5.7 for all pairs),
   63% of them over-predicted.
2. **A perfect distogram still emits 3.640 Å on FAIL18**, against a pool best of 2.284 Å
   and an ORACLE-single-best emitted of 2.285 Å. So even with every CA-CA distance exactly
   right, top-75-by-L1-risk → coordinate average → project loses 1.36 Å. This is a
   quantitative version of the standing result that the distance objective does not rank
   the native: **fixing the distogram completely does not fix FAIL18.** The remaining
   1.36 Å is the operator/argmin gap, and it is larger than every non-selection oracle
   channel except the full-distance one.
3. **Direction alone carries half the value** (o_sign −1.217 of the −2.379 total), with no
   magnitude information whatsoever. On the *controls* direction is worth −0.199 and
   magnitude-with-random-sign is worth −0.006 — direction is the whole story there. This
   is the formulable target: a learned per-pair **error-direction** classifier
   (is this predicted distance too long or too short?) is worth ≈ 1.2 Å on FAIL18 and
   ≈ 0.2 Å elsewhere, and needs no calibration at all.

### Q1: the distogram partially knows which pairs it got wrong

| | FAIL18 | MATCH18 |
|---|---|---|
| r(per-pair sd, \|error\|) | **0.492** | 0.479 |
| r(entropy, \|error\|) | −0.004 | 0.109 |
| r(separation, \|error\|) | **0.672** | 0.447 |
| sd-percentile of the 10 worst pairs | 0.713 | 0.749 |
| mean separation of the 10 worst pairs | 9.12 | 8.16 |
| fraction of the 10 worst that are over-predicted | 0.633 | 0.667 |

The model's own sd is an honest, moderately informative error flag (r ≈ 0.49, equally on
both groups) — the FAIL18 errors are bigger, not less anticipated. Entropy is useless.

---

## E15. The deployable version of the pair channel: sd-weighted scoring

`s12/fail_sdweight.py`, all 126, no new model, no native information — just reweight the
existing Bayes-risk score by the distogram's own per-pair sd.

| arm | top-75 best all-126 | Δ | 95% CI | FAIL18 | other-108 | recall FAIL18 | drop-top-10 |
|---|---|---|---|---|---|---|---|
| shipped (uniform) | 2.306 | 0 | — | 4.677 | 1.911 | 0/18 | — |
| weight ∝ sd^−0.5 | 2.301 | −0.005 | [−0.042, +0.033] | 4.563 | 1.924 | 1/18 | +0.033 |
| weight ∝ sd^−1 | 2.258 | −0.049 | [−0.137, +0.026] | 4.275 | 1.921 | 5/18 | +0.048 |
| **weight ∝ sd^−2** | **2.192** | **−0.115** | **[−0.230, −0.012]** | **3.868** | **1.912** | **9/18** | +0.016 |
| drop the 20 highest-sd pairs | 2.283 | −0.023 | [−0.088, +0.050] | 4.438 | 1.924 | 4/18 | +0.047 |
| drop every pair with \|i−j\| ≥ 9 | 2.235 | −0.072 | [−0.163, +0.018] | 4.101 | 1.924 | 7/18 | +0.026 |

`invsd2` has the profile the sprint is looking for and the others did not: **FAIL18 top-75
best 4.677 → 3.868 Å and recall 0/18 → 9/18, at literally zero cost to the other 108
(1.911 → 1.912)**, unlike the diversity lever which paid 0.84 Å on the controls (E11).
Its drop-top-10 is +0.016 — but that is *expected and not disqualifying here*: the effect
is by construction concentrated on 18 of 126 targets, so removing the ten biggest gainers
removes the mechanism. The honest statement is: this is a FAIL18-subgroup intervention
worth 0.8 Å of reachability there and nothing elsewhere, not a general improvement.

## E16. …and the sd-weighted objective does not survive the emitted test either

`s12/fail_sdemit.py`, same 36 targets, the identical selection projected through
`coordinate_average → project(λ=0.3)`; the baseline is re-used from E11 (same selection).

| group | base | invsd2 | Δ | W/L | 95% CI |
|---|---|---|---|---|---|
| FAIL18 | 6.019 | 5.850 | **−0.169** | 12/6 | [−0.429, +0.092] |
| MATCH18 | 2.259 | 2.491 | **+0.232** | 6/12 | [+0.041, +0.450] |
| all 36 | 4.139 | 4.171 | +0.031 | 18/18 | [−0.154, +0.208] |

FAIL18 wins are real where they happen (1JBF −1.11, 3SGO −1.09, 9KAR −1.02, 7LCW −0.90)
but the group CI includes zero, and the controls are significantly harmed.

### The unified negative result (E11 + E16), which is the most useful thing here

Two structurally unrelated filter modifications — geometric diversification and per-pair
uncertainty weighting — both:

| | FAIL18 recall | FAIL18 top-75 best | other-108 top-75 best | FAIL18 emitted | MATCH18 emitted |
|---|---|---|---|---|---|
| divmax | 0/18 → **14/18** | 4.677 → **3.306** | 1.911 → 1.892 | −0.453 | **+0.844** |
| invsd2 | 0/18 → **9/18** | 4.677 → **3.868** | 1.911 → 1.912 | −0.169 | **+0.232** |

They make the answer *reachable* and then lose most of it at the terminal operator,
because the coordinate average tracks the set MEAN and every intervention that widens the
set to include the band also widens it to include worse members. **The filter is not the
binding constraint on its own.** The binding constraint is the pair
(filter that contains the answer, selector that can find it) — and E4 shows every
native-free in-band signal on FAIL18 is at or below chance once its null is subtracted.
This is the same wall the record hits from the other side ("nothing ranks within the
pool", "consensus is the only in-band discriminator, and it is outlier avoidance"). Our
contribution is to show it is *not* a filter-capacity problem: the filter can be made to
hold the band on 14/18 of the failures with a one-line change.

---

# Conclusions and what to build next

## Answers to the three questions asked

**Why do the 18 fail?** Not because the answer is missing. Their K=500 pool best is
2.284 Å and the projection path is lossless (ORACLE single-best emits 2.285 Å). They fail
because the shipped distogram objective is 2.5× more wrong on their natives
(MAE 4.37 vs 1.76 length-matched; separation 5–8 shell 5.13 vs 2.03), which puts the
native itself at the **78th percentile** of its own pool's score, makes ρ(score, RMSD)
0.107 instead of 0.734, and collapses the top-75 onto one wrong basin
(intra-set spread 2.879 Å vs 4.608 Å for a random 75). **17/18 are objective failures**
(vs 7/18 when the same classifier is run on the matched controls); 12/18 also lose ~0.5 Å
at the query; exactly **1/18** is a pure filter-cut failure; and the label-convention
explanation is **refuted** — their NMR ensembles are no floppier than the controls'
(spread 0.942 vs 0.946 Å) and the ensemble-aware metric buys the same 0.18–0.22 Å in both
groups (3BTB is the one genuine label problem, spread 4.27 Å).

**What do they have in common?** Two structural classes, at 10× and 12× enrichment
(fibril-or-lasso: 10/18 vs 6/108, Fisher p = 1.2e-06), and one shared physical defect: the
compactness channel is broken in *opposite directions* for the two classes — the distogram
implies rg 5.8–7.4 Å for amyloid segments whose natives are 10.0–12.1 Å, and rg 7.6–10.1 Å
for lasso peptides whose knotted natives are 5.1–5.3 Å. Membrane/micelle peptides, the
largest class in the database, are *not* enriched (they are depleted).

**What information would have let a real predictor recognise the correct structure?**
Ranked by measured emitted-Å value on FAIL18 (E12/E14, all ORACLE):

| rank | channel | FAIL18 Δ emitted | MATCH18 Δ | deployable proxy tested? |
|---|---|---|---|---|
| 1 | which pool member is best (selection) | **−3.735** | −0.814 | none exists (E4: all in-band signals ≤ chance after their null) |
| 2 | which 75 pool members are best | −2.689 | −0.481 | — |
| 3 | every CA-CA distance exactly right | −2.379 | −0.318 | — (and still only reaches 3.640 Å) |
| 4 | the 20 most-wrong distances | −2.187 | −0.327 | sd-weighting: 0.8 Å of reachability, lost at the operator (E15/E16) |
| 5 | the 10 most-wrong distances | −1.493 | −0.265 | " |
| 6 | **the DIRECTION of every pair error (no magnitude)** | **−1.217** | −0.199 | **untested — the best-shaped open target** |
| 7 | the 3 most-wrong distances | −0.651 | −0.124 | " |
| 8 | native radius of gyration | −0.841 | −0.008 | distogram-implied rg: MAE 2.44 Å on FAIL18 — fails |
| 9 | native 8 Å contact map | −0.833 | +0.058 | ESM-2 contacts: net ≈ 0 over their null on FAIL18 (E7) |
| 10 | native secondary structure | −0.788 | +0.015 | propensity SS AUC 0.332; as a retrieval key even the ORACLE is worth 0.11 Å (E10) |

## Ranked recommendations

1. **Build a per-pair error-DIRECTION head, not a better distogram.** Direction alone
   (a 2 Å nudge the right way, zero magnitude information) is worth −1.217 Å on FAIL18 and
   −0.199 Å on the controls, 18/0 W/L. It is a binary classification problem per pair, it
   is trainable leave-fold-out on `rr`/`nat_ca` as labels, and it needs no calibration —
   which sidesteps the standing "calibration slope 0.376, MAE is unusable" wall entirely.
   Feature the existing distogram's sd (r = 0.49 with |error|), separation (r = 0.67 on
   FAIL18), and the ESM contact map, which does have net content on the 108.
2. **Stop trying to fix FAIL18 from the filter alone; the operator eats it.** Two
   independent filter changes each recovered 9–14 of 18 recall and each lost it at the
   coordinate average (E11, E16). Any filter work must be paired with a selector, and the
   selector is the thing that does not exist. If a diversified filter is combined with an
   in-band ranker of even modest skill, the headroom is 4.677 → 3.306 Å of top-75 best.
3. **Declare the lasso class out of scope, explicitly.** All 6 lasso peptides in tuning126
   emit ≥ 6.10 Å (mean 6.71). Of ~10,000 windows per target, a mean of 2.3 are within
   2.0 Å and three of the six have none; the best window in the entire universe averages
   2.01 Å. This is a coverage limit of linear-window retrieval, not a ranking failure, and
   the library-saturation result says more fragments will not fix it. It costs ≈ 0.17 Å of
   the tuning mean. Say so rather than absorb it.
4. **Treat the amyloid/fibril class as a retrieval-key problem, but do not expect much.**
   371 windows per target are within 2.0 Å (universe best 1.20 Å) and BLOSUM ranks the best
   at position ~12,000. However E10 shows even a *perfect* SS retrieval key is worth only
   −0.104 Å of pool best; ORACLE rg is worth −0.177 Å. **S6-6's open hypothesis (SS as a
   retrieval key) should be recorded as closed and negative.**
5. **Retire two attractive leads.** (a) ESM-2 contact agreement as a scoring channel:
   real content on the other-108 (net AUC +0.149) but net ≈ +0.03 [−0.08, +0.16] on
   FAIL18, and every fused arm is matched by a foreign-contact-map null and dies to
   drop-top-10. (b) Consensus/centrality on hard targets: on FAIL18 the band is the pool's
   *outlier* set (AUC 0.252), so every centrality-flavoured method is guaranteed to fail
   exactly where it is needed.
6. **Do not expect a router.** LFO AUC 0.558 against a permutation null of 0.511
   [0.314, 0.680], p = 0.295, on 13 native-free features. Regime-dependent blending is not
   available; a fixed setting is the honest choice.
7. **A ceiling worth knowing.** Even a *perfect* distogram emits 3.640 Å on FAIL18 through
   the shipped top-75 + average + project path, against a 2.285 Å ORACLE-selection floor.
   1.36 Å of the FAIL18 error is the operator/argmin gap and survives any improvement to
   the distance objective. That number should replace any assumption that a better
   distogram alone can close this subgroup.

## Nothing here is proposed for a dev24 pass.
Every positive in this report either fails its own null (E5, E7, E8), fails the
emitted-structure test (E11, E16), or is an ORACLE diagnostic (E12, E14). The one item
that would justify a dev pass is the error-direction head of recommendation 1, and it has
not been built.
