# s12 — LITERATURE & FOUNDATION-MODEL agent — findings

**Question.** What does the field know in 2026 that this project is not using, that (a) runs on 8 CPU
cores / ~5 GB free, (b) supplies information ORTHOGONAL to BLOSUM retrieval + the ESM-2 distogram,
and (c) could rescue the 18 zero-recall targets?

**Answer in one line.** The field's standard fragment-retrieval key has been *predicted local
backbone conformation* since 2011 (Rosetta fragment picker) and *is the entire architecture* of the
only two peptide-specific predictors that beat this system (PEP-FOLD, APPTEST). This project retrieves
by BLOSUM62 sum alone. Measured here, swapping the key is worth **2.284 → 1.640 Å pool-best on FAIL18**
and takes near-native band recall from **16/500 to 92/500** — and it survives 30 % prediction error.
Second: **56 % of FAIL18 is two chemically identifiable classes** (steric-zipper amyloid segments,
lasso peptides) that no linear-window retrieval can ever represent, and both have dedicated 2024–2025
sequence-only tools.

Files: `s12/results/lit_methods.json` (24 method rows), `lit_ideal_templates.json`,
`lit_distogram_extension.json`, `lit_abego_retrieval.json`, `lit_abego_noise.json`,
`lit_abego_deployable.json`, `lit_bmrb_coverage.json`.

**Compute discipline.** No heavy compute. All local work is array arithmetic on already-cached
universes plus 126 RCSB REST calls; no projection, no AMBER, no model inference, `OMP_NUM_THREADS=2`,
single process < 1 GB. **No dev24 run.** Natives used only as ORACLE labels for diagnostics, every
such arm labelled.

**Baseline caveat.** `shipped_record['shipped']` is the **distogram-argmin arm (mean 3.4540,
FAIL18 6.008)**, which is how FAIL18 is defined. Deltas below are on that arm. They must be
re-derived on the synthesis arm (3.204/3.236) before any of them is claimed as a system improvement.

---

## Part 0 — instrument validation

```
python -m s12.instrument  ->  3.4540 / 1.7108 / 2.3062 / 3.2041 / 18   (exact match)
```
`fail_headers.json` (produced by another agent) reproduced independently for the class counts below.

---

## Part 1 — the FAIL18 is not 18 random hard targets; it is two chemistry classes plus 8 others

| class | FAIL18 | other-108 | odds ratio | Fisher p |
|---|---|---|---|---|
| fibril / amyloid | 6/18 (33 %) | 4/108 (3.7 %) | 13.0 | 5.5e-4 |
| lasso peptide | 4/18 (22 %) | 2/108 (1.9 %) | 15.1 | 3.8e-3 |
| **either** | **10/18 (56 %)** | **6/108 (5.6 %)** | **21.2** | **1.2e-6** |

FAIL18 lasso members: `2N5C` chaxapeptin, `7JS6` des-citrulassin F, `7LCW` lihuanodin, `9KAR`
streptacidin. All four carry the canonical lasso motif (Gly at position 1, Asp/Glu at 8–9 forming the
macrolactam), i.e. **the class is detectable from sequence alone**.

FAIL18 fibril members: `2BFI`, `2BP4`, `2JN5`, `3SGO`, `5W52`, `9L1M`.

### 1a. Literature check — is this a known failure mode?

| finding | number | source |
|---|---|---|
| AF2 / AF3 / ESMFold fail on monomeric **NMR** structures | 128/190 (67.3 %), 135/190 (71 %), 139/190 (73 %) vs 95 %/83 % success on X-ray/cryoEM | NAR Genom Bioinf 2026 `lqag002` |
| AF2 on 588 NMR peptides 10–40 aa: no rank/pLDDT ↔ RMSD correlation; fails on kinks, turns, extended flexible regions | RMSD/residue 0.068 (disulfide-rich) → 0.202 (mixed membrane) | Structure 2023, PMC9883802 |
| AF2/ESMFold **cannot produce the lariat topology at all**; AF3 fails 10/12 PDB lasso peptides, 8 % success on the RODEO set | 83 % / 92 % failure | Nat Commun 2025 `s41467-025-60544-4` |
| ESM-2's training distribution barely contains peptides | 2.8 % of UniRef50 < 50 residues | PMC12844563 |

So: the isolated 9–16-mer NMR peptide is a **documented blind spot of the entire AF3-like family**,
and the lasso class is a documented total failure. The project is not behind the field on these
targets; the field has no answer either. That is worth saying plainly in any write-up.

### 1b. ORACLE diagnostic — the amyloid class is solved by a two-parameter model

Built ideal CA traces from a *single constant* (φ,ψ) and measured CA-RMSD to native.
No retrieval, no distogram, no ESM, no AMBER. (`lit_ideal_templates.json`)

| target | class | fully-extended template (−180,180) | shipped (argmin arm) | pool best | universe best |
|---|---|---|---|---|---|
| 2BFI | zipper, X-ray | **0.27** | 7.03 | 2.22 | 1.35 |
| 5W52 | zipper, microED | **0.43** | 5.07 | 1.12 | 0.47 |
| 5V5B | zipper | **0.65** | 3.61 | 1.07 | 0.68 |
| 7VI4 | zipper | **0.72** | 3.90 | 0.89 | 0.61 |
| 7N2I | zipper | **1.23** | 2.50 | 0.74 | 0.74 |
| 3SGO | zipper, X-ray | **2.41** | 6.62 | 1.60 | 1.04 |
| 1S9Z | designed coiled coil (fibril-forming) | 8.94 | 0.30 | 0.16 | 0.16 |
| 2BP4 | Aβ in TFE (helical) | 9.07 | 5.91 | 0.63 | 0.50 |
| 2JN5 | α-syn bound (helical) | 10.32 | 5.61 | 2.69 | 1.84 |
| 9L1M | low-complexity region | 8.84 | 5.14 | 2.96 | 1.98 |

**On 2BFI the two-parameter model beats the best window in the entire 21,547-window universe** (0.27
vs 1.35). The information needed is not in the library and never was; it is a class label.

Ceiling if a perfect steric-zipper detector routed those 6 targets to an extended template:

```
126 mean   3.4540 -> 3.2712   (-0.183)
FAIL18     6.008  -> 5.141    (-0.867)   [3 of the 6 are in FAIL18]
```

Adding LassoPred at its *published* 3.2 ± 1.0 Å for the 4 lasso FAIL18 targets:

```
126 mean   3.4540 -> 3.1612   (-0.293)
FAIL18     6.008  -> 4.370    (-1.638)
```

Detector literature is mature and CPU-free: WALTZ-DB 2.0 (1,416 assayed hexapeptides), CORDAX
(Bioinformatics 2024, predicts the zipper topology class), and a 2025 proteome-scale ML model
replicating ZipperDB at 93.5 % classification overlap. **None of these has ever seen a PDB
coordinate of these targets — they are trained on in-vitro aggregation phenotypes. Zero leakage.**

**Honest limits.** The routing gain is bounded by detector precision: `2JN5` (α-syn NAC-adjacent,
plausibly zipper-positive by sequence) is a *helix* in its deposited structure, so a false positive
there costs +4.7 Å on one target. The experiment must report precision/recall of the detector, not
just the ceiling.

---

## Part 2 — the distogram systematically under-predicts chain extension

Normalised end-to-end distance `d(1,n) / (3.8·(n−1))`; 1.0 = fully extended. (`lit_distogram_extension.json`)

| quantity | value |
|---|---|
| corr(predicted, true) normalised end-to-end, n=126 | r = 0.499 |
| corr(**\|**true − predicted**\|**, shipped RMSD) | **r = 0.531** |
| n with gap > 0.25 | 12 |
| mean shipped RMSD, gap > 0.25 vs rest | **4.898 vs 3.302** |
| FAIL18 fraction, gap > 0.25 vs rest | 0.33 vs 0.12 |

Worst offenders (gap, predicted, true, shipped): 2MK7 (0.56, 0.27, 0.83, 3.92), 2BFI (0.47, 0.43,
0.90, 7.03), 3SGO (0.43, 0.36, 0.79, 6.62), 7VI4 (0.38), 5W52 (0.38), 8T63 (0.37), 2MSA (0.36).

This is **S7-2's distribution shift made numeric**: the distogram was trained on fragments cut out of
folded proteins, which are more compact than free peptides, and it regresses global extension toward
the compact mean. A single scalar — the chain's global dimension — is mis-estimated, and its error
carries r = 0.53 with the final RMSD.

**The field has a purpose-built, zero-leakage, CPU-instant predictor for exactly this scalar:**
ALBATROSS (Nat Methods 2024) predicts Rg, end-to-end distance, scaling exponent ν and asphericity
from sequence, trained **entirely on coarse-grained simulations validated against SAXS — no PDB
coordinates whatsoever**. Sub-million-parameter BiLSTM, milliseconds, pip-installable (`sparrow`).
The obvious caveat is that its domain is IDRs, not structured peptides; that is a one-hour falsifier
(does predicted Rg correlate with native Rg on the 126 at r > 0.4?).

---

## Part 3 — THE MAIN RESULT: the retrieval key is the wrong question

### 3a. What the field actually uses as a fragment key

| system | retrieval key |
|---|---|
| Rosetta fragment picker (PLoS One 2011) | weighted sum of **sequence profile + predicted secondary structure + `torsion_bin_probs`** |
| PEP-FOLD 1–4 (NAR 2009/2012/2016/2023) | **per-position 27-state structural-alphabet profile**, SVM over a PSI-BLAST PSSM. This *is* PEP-FOLD. |
| CS-Rosetta (Bax lab) | **chemical-shift-based fragment selection** (again local conformation, not sequence) |
| APPTEST (Brief Bioinform 2021) | no retrieval at all: predicts distances **and φ/ψ**, folds by restrained annealing |
| **this project** | **BLOSUM62 sum. Only.** |

Every peptide-specific method that outperforms this system picks fragments by predicted **local
conformation**. This project picks by raw sequence similarity — the weakest of the Rosetta picker's
three terms, and in its weakest form (a BLOSUM sum, not even a profile score).

Recorded in memory: *"structure and sequence are decoupled — best-matching window has 12 % identity;
sequence is the wrong retrieval key."* The literature agrees, and names the replacement.

### 3b. ORACLE experiment — torsion-bin (ABEGO) key vs BLOSUM key

4-state ABEGO bins from real parent torsions. Key = the ABEGO string of the true-best window
(**ORACLE**). Retrieve K=500 by bin-match count, same stable-argsort tie-breaking.
46 targets (all 18 FAIL18 + every 4th other). (`lit_abego_retrieval.json`)

| group | n | pool-best BLOSUM | pool-best ABEGO\* | pool-mean BLOSUM | pool-mean ABEGO\* | band members BLOSUM | band ABEGO\* |
|---|---|---|---|---|---|---|---|
| FAIL18 | 18 | 2.284 | **1.640** | 5.96 | 5.11 | 16.0 | **91.6** |
| other (sample) | 28 | 1.366 | **1.050** | 4.01 | 2.94 | 108.3 | 251.0 |
| all | 46 | 1.725 | **1.281** | 4.78 | 3.79 | 72.2 | 188.6 |

("band" = windows within pool-best + 1.5 Å, out of 500.)

The ABEGO key recovers the **universe best** on 17/18 FAIL18 targets. Per-target band recall:
2BFI 2→139, 3SGO 15→307, 5W52 3→133, 2NDM 8→167, 9L1M 12→121, 2BP4 91→420.

**This directly attacks the defining FAIL18 defect.** The stated problem is "the top-75 keeps NONE of
the near-native band". With 92/500 band members instead of 16/500, a filter with the *existing*
in-band skill (ρ ≈ +0.13) retains some by chance alone.

### 3c. How accurate must the torsion-bin predictor be? (`lit_abego_noise.json`, `lit_abego_deployable.json`)

Corrupt a fraction p of the oracle key. Two error models: uniform-random substitution, and the
harsher **biased** model (corrupted positions collapse to the string's modal bin — how a real
predictor fails).

| key | FAIL18 pool-best | other pool-best |
|---|---|---|
| BLOSUM (current) | 2.284 | 1.366 |
| ABEGO oracle (p=0) | 1.640 | 1.050 |
| uniform noise p=0.2 | 1.641 | 1.085 |
| uniform noise p=0.3 | 1.669 | 1.114 |
| uniform noise p=0.5 | 1.780 | 1.217 |
| **biased noise p=0.2** | **1.694** | 1.157 |
| **biased noise p=0.3** | **1.701** | 1.212 |
| **biased noise p=0.5** | **1.865** | 1.326 |
| random key (p=1.0) | 2.533 | 2.035 |
| **deployable: modal ABEGO of the shipped top-75** | **3.890** ✗ | 1.437 |
| soft version of the same | 4.020 ✗ | 1.483 |
| soft + BLOSUM z-sum | 2.946 ✗ | 1.337 |

Two conclusions, one positive and one negative, both important:

1. **A ~70 %-accurate torsion-bin predictor beats BLOSUM as a retrieval key**, even under the biased
   error model. Only a *fully random* key is worse than BLOSUM. This is a low bar: LOCUSTRA (2008,
   SVM) reaches Q16 = 61 % on the 16-state Protein Blocks alphabet and Q3 = 79 %; PBC5.0d reports
   Q16 ≈ 80 %; 4-state ABEGO is easier than 16-state PB.
2. **Bootstrapping the key from the pipeline's own output does NOT work.** The modal ABEGO of the
   shipped top-75 agrees with the oracle on only 44 % of FAIL18 positions and makes FAIL18 pool-best
   *worse* than BLOSUM (3.890 vs 2.284). This is `consensus-is-outlier-avoidance` reappearing exactly
   where predicted: a failed pool cannot diagnose itself. **The predictor must be INDEPENDENT of the
   pool** — trained on sequence/ESM features under leave-fold-out, never on the retrieved set.

The required accuracy therefore sits in a measurable window: **above ~0.44 (self-consensus, fails)
and at or above ~0.70 (biased-0.3, works)**. That is the number the falsifying experiment must produce.

**Caveat I will not hide.** Corrupting an oracle string is optimistic even under the biased model,
because it preserves the string's *length-scale structure*. A real predictor that outputs "all A"
for a helical-looking sequence would behave like the p=1.0 row on the extended targets. The gate
below is written to catch that.

---

## Part 4 — the terminal operator has no escape hatch

The pipeline can only emit a convex combination of retrieved windows, projected onto ideal geometry.
On FAIL18 the retrieved set contains nothing near-native, so no filter, reranker, weighting or
quantum selector can help — a fact the record has now established eleven different ways
(`nothing-ranks-within-the-pool`, `score-axis-does-not-transfer`, `in-band-is-the-only-ranking-metric`).

APPTEST is the existence proof that a different terminal operator works at this length:

| set | n | APPTEST | PEP-FOLD | PEPstrMOD |
|---|---|---|---|---|
| short, 9–25 aa | 42 | **1.96 Å** | 2.05 | 4.66 |
| long, 26–40 aa | 30 | 3.20 | 3.43 | — |
| cyclic, 10–40 aa | 34 | 1.99 | 2.71 | — |

Selection is by lowest XPLOR-NIH energy / CYANA target function — **deployable, not oracle**.
Architecture: small 1D-CNN + gated residual blocks → predicted CA–CA/CB–CB distances **and φ/ψ** →
restrained simulated annealing. Two things this project does not have: a **torsion head**, and a
**fold-from-restraints operator with no retrieval bottleneck**.

The project already has the expensive half (a trained distogram) and the machinery
(`core/project.py` does L-BFGS over (φ,ψ) with multi-start). What is missing is running that
optimiser against the **distogram** rather than against a retrieved coordinate average.

---

## Part 5 — the NMR channel (task-change, not leakage)

Measured (RCSB REST, `lit_bmrb_coverage.json`): **79/126 targets (63 %) and 10/18 FAIL18 have a
linked BMRB entry with deposited chemical shifts.**

TALOS-N converts shifts to φ/ψ at **~12° RMSD for ≥90 % of residues**. Over a 13-mer that is roughly
~1 Å CA-RMSD — it would essentially solve the shift-covered subset.

This is **not leakage** (shifts are experimental observables, not deposited coordinates) but it **is a
task change**: from *predict from sequence* to *determine from NMR data*. With 63 % coverage it cannot
be a system component. Its value is as a **ceiling instrument**: it answers the question the record
has been circling for three sprints — *is the residual 3.2 Å an information problem or an operator
problem?* If a shift-derived torsion prior takes the 79 covered targets to ~1.5 Å, every further
reranking effort is misdirected and the answer is "more information". If it does not, the operator is
the problem. Either way the sprint learns something decisive.

---

## Part 6 — what the field has that will NOT work here (funded negatives)

| candidate | why not |
|---|---|
| AF2 / AF3 / Boltz-1/2 / Chai-1 / Protenix | GPU-only; **and** they fail on 67–73 % of monomeric NMR structures; and MSA depth for a 13-mer is zero, which removes their only advantage. Leakage HIGH (PDB ≤ 2021-09-30). |
| ESMFold | needs ESM-2 **3B** (≈5.6 GB fp16) — over budget; and its errors share a source with the project's ESM-2 distogram features, so low orthogonality. |
| Foldseek 3Di / ProstT5 / SaProt | 3Di encodes each residue's relation to its nearest **tertiary** neighbour, and the Foldseek paper states its information density is *lowest* in coils. A 13-mer has almost no tertiary contacts → the alphabet degenerates. The right alphabet for peptides is a **backbone** alphabet (ABEGO / Protein Blocks / PEP-FOLD's 27-state SA). ProstT5 is a 3B T5 — over budget anyway. **Do not fund this.** |
| PepFlow / PPFlow / PepGLAD / THFlow | Every 2024–25 "peptide generative model" is **receptor-conditioned binder design**. None predicts the solution conformation of a free peptide. Task mismatch. |
| RFdiffusion / Genie2 / FrameFlow / Chroma | Trained to produce *designable, idealised, well-packed* backbones. This benchmark's hard cases are floppy, extended, membrane-bound and topologically exotic — the prior is actively wrong-signed. GPU-only. |
| MDGen | Transferable regime is **tetrapeptides**. Length mismatch is disqualifying. |
| BioEmu | Genuinely orthogonal (MD-derived) but needs the AF2 evoformer + GPU. Out of scope for this box. |
| SimpleFold-100M | The only end-to-end folder that might fit (100M params on 13 tokens ≈ 1 min CPU). **But** it conditions on ESM-2 embeddings; if the required size is 3B it is infeasible, and if it is 650M it shares the project's error source. Verify the ESM dependency *before* spending time. Leakage HIGH (PDB + AFDB distillation). |

---

## Part 7 — leakage audit for everything recommended

| component | training data | benchmark contamination |
|---|---|---|
| torsion-bin head (H1) trained **locally**, LFO | this project's own out-of-fold windows | **NONE** by construction (pinned folds) |
| WALTZ-DB / CORDAX zipper propensity | 1,416 in-vitro hexapeptide aggregation assays | **NONE** — no PDB coordinates in the training signal |
| ALBATROSS | Mpipi-GG coarse-grained simulations, SAXS-validated | **NONE** — has never seen a PDB structure |
| metapredict V3 | disorder consensus + AFDB pLDDT distillation | LOW |
| ProteinMPNN | PDB ≤ 2020, ~1.7M chains at 30 % identity | LOW–MEDIUM: it scores *sequence given backbone*, so target-specific memorisation is indirect — but it is a PDB-trained component and must be labelled |
| LassoPred | 47 lasso PDB structures, cutoff **2024-03-15** | **HIGH for 3 of 4**: `2N5C` (2015), `7JS6` (2020), `7LCW` (2021) are almost certainly in the 47. **`9KAR` (released 2025-04-16) is post-cutoff and is the only clean instance.** Must be reported target-by-target. |
| TALOS-N / UCBShift | shift↔structure pairs; UCBShift's *transfer* module explicitly looks up structural homologs with assigned shifts | **task change**, plus a direct memorisation path in UCBShift's transfer module — use the ML module only |
| any imported protein SS predictor (NetSurfP-3.0, S4PRED, Porter) | PDB chains, standard protocols **exclude chains < 20 residues** | MEDIUM, and accuracy on 9–16-mers is unmeasured. Prefer a local LFO retrain. |

---

## Part 8 — RANKED HYPOTHESES

Ranked by (expected gain × orthogonality × FAIL18 relevance) ÷ (compute × leakage).

---

### H1 — Retrieve by predicted torsion-bin profile, not by BLOSUM. ★★★★★

**Claim.** The retrieval key is the binding constraint on FAIL18, and the field has used a local-
conformation key for 15 years. An independent, leave-fold-out sequence→ABEGO-posterior head used as
the K=500 key takes FAIL18 pool-best from 2.284 to ≈1.70 Å and near-native band recall from 16/500 to
≈90/500.

**Evidence.** §3b, §3c. Oracle key 1.640; biased-noise-0.3 key 1.701; random key 2.533; BLOSUM 2.284.
Rosetta fragment picker `torsion_bin_probs` (2011); PEP-FOLD's SA profile *is* its architecture;
LOCUSTRA Q16 = 61 %, PBC5.0d Q16 ≈ 80 % — 4-state ABEGO is easier than either.

**Experiment (~4 h).** Train a small MLP/BiLSTM on the banked ESM-2 PCA-32 features (`s12/esm_bank.py`)
→ 4-way ABEGO posterior per residue, **leave-fold-out on the pinned folds, never touching the pool**.
Report per-residue accuracy on FAIL18 separately from the 108. Then retrieve K=500 by expected
bin-match and by (z(bin-match) + z(BLOSUM)), and report pool-best / pool-mean / band-recall, paired
with bootstrap CI, FAIL18 vs other-108, per-fold, drop-top-10/20.

**Gate.** If per-residue ABEGO accuracy on FAIL18 < 0.60, stop — §3c says the key will not clear
BLOSUM. Also require a **shuffled-label control** (train on permuted ABEGO targets): it must land at
the p=1.0 row (≈2.53), not near BLOSUM.

**Falsified if.** Predicted-key pool-best on FAIL18 ≥ BLOSUM's 2.284, or band recall does not rise.

**Risk.** The corruption model is optimistic; a real predictor's errors are correlated along the
chain. The gate and the shuffled control are there to catch that.

---

### H2 — Conformational-class routing: detect steric-zipper segments, emit an extended strand. ★★★★☆

**Claim.** 6/126 targets are cross-β steric-zipper segments whose native is reproduced to
0.27–2.41 Å by a *two-parameter* model, while the pipeline returns 2.50–7.03 Å. A sequence-only
amyloid detector routes them out of the retrieval path entirely.

**Evidence.** §1b. Ceiling −0.183 Å on the 126 mean, **−0.867 Å on the FAIL18 mean**. Detector:
WALTZ-DB 2.0 (1,416 assayed hexapeptides), CORDAX (Bioinformatics 2024), 2025 zipper-propensity ML
(93.5 % overlap with ZipperDB). **Zero leakage** — trained on aggregation assays, not coordinates.

**Experiment (~3 h).** Train a hexapeptide amyloid classifier on WALTZ-DB (external, no PDB
coordinates), scan each target's 8 hexapeptide windows, route positives to a strand template whose
(φ,ψ) is fit by the *existing* projection objective under the distogram — not to the oracle (−180,180).
Report **detector precision/recall on the 126** before reporting any RMSD delta.

**Gate.** Detector precision must exceed ~0.5 on the 126; one false positive on a helical target costs
+4.7 Å (see `2JN5`).

**Falsified if.** The detector cannot separate the 6 zippers from the 120 others at usable precision,
or the routed template (fit, not oracle) does not reach < 2 Å on the 6.

---

### H3 — A de-novo fold-from-restraints arm that bypasses retrieval entirely. ★★★★☆

**Claim.** On the 18 zero-recall targets the pipeline is structurally incapable of emitting a
near-native answer, because its output space is the convex hull of the retrieved windows. APPTEST
shows a restraint-driven annealing operator reaches 1.96 Å on 42 short peptides with deployable
(energy-based) selection.

**Evidence.** §4. The project already owns both halves: a trained distogram and an L-BFGS
(φ,ψ) optimiser with multi-start and a Ramachandran prior (`core/project.py`). They have simply never
been connected — the optimiser currently fits a *retrieved coordinate average*, not the distogram.

**Experiment (~6 h).** Multi-start L-BFGS over (φ,ψ) minimising the distogram's L1 Bayes risk
directly + the hinged Ramachandran prior, ~64 restarts (seeded from ABEGO-consistent torsion draws if
H1 lands). Select by objective value (deployable). Compare against the retrieval average, paired, on
FAIL18 vs other-108. Cost: the projection is 2–8 s/call, so this is minutes per target — the most
expensive hypothesis here, hence the rank.

**Falsified if.** The de-novo arm does not beat the retrieval arm on FAIL18. Note the *known* prior
risk: `objective-does-not-rank-the-native` says the native is the distance objective's argmin on
3/126 targets at the 37th percentile — so optimising the objective harder may well get worse. That is
precisely why this is a **hypothesis** and not a plan: it is the cleanest test of whether the
distogram objective is worth optimising at all.

---

### H4 — A global-dimension (Rg / end-to-end) prior corrects the distogram's compactness bias. ★★★★☆

**Claim.** The distogram's error in one scalar — global chain extension — carries **r = 0.531** with
the final RMSD across all 126 targets. ALBATROSS predicts that scalar from sequence in milliseconds
with literally zero leakage.

**Evidence.** §2. 12 targets with gap > 0.25 have mean shipped RMSD 4.898 vs 3.302. ALBATROSS:
Nat Methods 2024, R² = 0.994 vs Mpipi-GG Re, SAXS-validated, sub-million-parameter BiLSTM, CPU.

**Experiment (~2 h).** (1 h falsifier first) Does predicted Rg/Re/ν correlate with **native** Rg on
the 126 at r > 0.4? If yes: (a) re-rank the K=500 by |Rg(window) − Rg(predicted)|; (b) add a global
Rg term to the projection objective. Paired stats, FAIL18 vs 108.

**Gate.** ALBATROSS's domain is IDRs; it may be systematically wrong for the helical majority. The
1-hour correlation check kills it cheaply if so. A **length-matched null** (predict Rg from n alone)
is mandatory — most of the apparent correlation may be chain length.

---

### H5 — Lasso-peptide topology routing. ★★★☆☆

**Claim.** 4 FAIL18 targets are lasso peptides at 6.67 Å mean. No linear window from any library can
express a threaded macrolactam. LassoPred (Nat Commun 2025) reaches 3.2 ± 1.0 Å from sequence, on CPU.

**Evidence.** §1. AF2/ESMFold cannot make the topology; AF3 fails 10/12. Sequence motif (G/C at 1,
D/E at 8–9) makes detection near-trivial. Ceiling **−0.110 Å** on the 126 mean, **−0.770 Å** on FAIL18.

**Cost / risk.** Highest leakage of anything recommended: **3 of the 4 targets are almost certainly
inside LassoPred's 47-structure training set (cutoff 2024-03-15)**. `9KAR` (2025-04-16) is the only
clean instance. Licence is CC BY-NC-ND. So the gain is real but must be reported target-by-target
with the contamination stated, and it cannot be presented as a clean benchmark improvement.

**Cheapest version.** Skip the tool: just **detect** the class and report those 4 as
out-of-model-class. Honest scoping is worth more than a contaminated 0.11 Å.

---

### H6 — ProteinMPNN inverse-folding likelihood as an orthogonal in-band reranker. ★★★☆☆

**Claim.** Every in-band ranker tried so far (223 features, ridge/GBT/MLP → ρ ≈ +0.20, ≈0 Å) reads the
distogram or window statistics. ProteinMPNN asks a different question of the *geometry*: would this
sequence design this backbone? 1.7M parameters — the cheapest external model available, milliseconds
per candidate on CPU.

**Evidence.** Inverse-folding log-likelihoods correlate strongly with stability and mutational
fitness (arXiv 2506.05596). MIT-licensed, open weights.

**Experiment (~3 h).** Thread the target sequence onto each of the K=500 backbones, score by NLL,
report **in-band** Spearman on the 108 and on FAIL18 separately (`in-band-is-the-only-ranking-metric`).

**Gate.** Require standalone in-band ρ > 0.30 before attempting any fusion —
`decorrelated-errors-exist-but-are-unusable` shows fusion gain goes as the *square* of the weaker
channel's skill, so a ρ = 0.15 channel buys nothing. Tie-break by averaging over the tied argmin set
(`consensus-is-the-only-in-band-discriminator#the-methodological-trap`).

**Expected shape.** ProteinMPNN's designability prior is *wrong* for floppy/extended peptides — expect
it to help the 108 and hurt FAIL18. Report both.

---

### H7 — NMR-assisted ceiling arm (TALOS-N / shift-agreement reranking). ★★★☆☆ (diagnostic)

**Claim.** 63 % of the benchmark (79/126) and 10/18 FAIL18 have deposited BMRB shifts. TALOS-N gives
φ/ψ at ~12° RMSD for ≥90 % of residues — enough to essentially solve those targets.

**Value.** Not a system component (coverage too low, and it changes the task). It is the **decisive
diagnostic** for the question the last three sprints have circled: information-limited or
operator-limited? Run it, label it `NMR-ASSISTED / DIAGNOSTIC`, and let the answer redirect the sprint.

**Falsified if.** Even a shift-derived torsion prior fails to reach ~2 Å on the 79 — in which case the
operator, not the information, is the bottleneck, and H3 becomes the priority.

---

### H8 — Score against the deposited NMR ensemble, not against model 1. ★★☆☆☆ (nearly free)

**Claim.** The natives are NMR ensemble representatives (`fail_headers.json` records 14, 21, 40 models
for individual targets). Scoring one averaged structure against model 1 of a 40-model ensemble
conflates prediction error with ensemble spread. The 2024–25 peptide-ensemble literature
(AFsample2; AF2 peptide ensembles, RMSD < 2.5 Å over 557 peptides) treats the ensemble as the target.

**Experiment (~1 h).** Recompute the 126 against (a) the best deposited model and (b) the ensemble
mean, and report how much of the 3.2 Å is ensemble spread. Costs almost nothing and re-calibrates
every number in the record. **This changes the ruler, not the system — it can only be reported as
context, never as an improvement.**

---

### H9 — Full multi-term fragment score (profile + SS + torsion-bin), LFO-weighted. ★★☆☆☆

The generalisation of H1: reimplement the Rosetta picker's weighted score with learned LFO weights and
ablate the terms. Only worth doing **after** H1 shows the torsion term carries signal; otherwise it is
hyperparameter search on a key that does not work.

---

## Part 9 — what I would tell the coordinator to fund

1. **H1** — the retrieval key. Largest measured effect, directly on the defining defect, no leakage,
   cheap, and it is the one thing the whole peptide-prediction literature does that this project does not.
2. **H2** — steric-zipper routing. Zero-leakage detector, −0.87 Å on FAIL18, and it establishes
   class-conditional routing as an architecture.
3. **H4** — the 1-hour Rg falsifier. Cheapest possible test of an r = 0.53 error mode.
4. **H7** — the NMR ceiling arm, as a diagnostic, to settle information-vs-operator.
5. **H3** — the de-novo arm, if H1 lands (H1 supplies its torsion seeds).

**Do not fund:** 3Di/ProstT5/SaProt retrieval (wrong alphabet for peptides, over budget); any
AF3-like co-folding model (GPU, and they fail this task class); peptide generative design models
(receptor-conditioned, wrong task); backbone diffusion (designability prior is wrong-signed).

**The single most important structural claim in this report:** on FAIL18 the system is
*retrieval-limited, not ranking-limited*, and the entire 2025–26 reranking/selection literature —
like the project's own eleven negative reranking results — is aimed at the wrong stage.
