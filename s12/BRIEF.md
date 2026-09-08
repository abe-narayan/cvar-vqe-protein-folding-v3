# Sprint 12 — shared agent brief (READ FULLY BEFORE DOING ANYTHING)

You are one of several research agents on a maximum-depth sprint whose objective is to find
the MISSING INFORMATION CHANNEL that would let a 9–16-residue peptide structure predictor reach
<2.0 Å mean CA-RMSD on unseen structures, with full scientific integrity. This file is the
shared context. The coordinator (the parent session) integrates your report with the others'.

## 0. Hard rules (non-negotiable)

1. **benchmark60 is SPENT and PROTECTED.** Never read `results/benchmark_manifest.json`
   targets, never load `s9/final_cache/*` (those 60 npz are the benchmark pools), never run
   anything on the benchmark, never use benchmark numbers to make a decision.
2. **Development instrument = tuning126** (the 126 npz in `s8/generate_univ/`). **dev24** is
   the cluster-disjoint confirmation split: do NOT run on dev24 unless the coordinator asks.
   If you believe a result deserves a dev pass, say so in your report and stop.
3. **Natives are labels only.** `rr` (per-window true CA-RMSD) and `nat_ca` may be used for
   evaluation and as leave-fold-out TRAINING LABELS. They must never enter an inference-time
   decision. Every deployable rule must be computable from the deployable view (window
   coordinates/torsions/sequence codes/BLOSUM sim, the shipped distogram, ESM features of
   the TARGET sequence, energies). If you use an oracle, label the arm ORACLE/DIAGNOSTIC.
4. **Leave-fold-out discipline.** Folds are pinned (`peptide_folds.json`, 5 folds, keyed by
   sequence; each tuning target's fold is in `instrument.targets()`). Any learned component
   scoring target T must be trained WITHOUT fold(T). The library windows for T are already
   out-of-fold (that is what the universe files contain).
5. **Do NOT modify** `core/`, root modules, `s5/ s7/ s8/ s9/`, `tests/`, `verify/`, `pdbs/`,
   `prots/`, `pdbs_ext/`, `peptide_*.json`, `*.npz` at the root, `bench_results/`. Do not
   run `git` write commands. Put ALL new code in `s12/<yourname>_*.py`, results in
   `s12/results/<yourname>_*.json`, and your findings in `s12/<yourname>_FINDINGS.md`.
6. **Memory.** This box has 16 GB and other processes hold ~10 GB. Call
   `s12.instrument.free_gb()` before any heavy step; keep each Python process under ~1.2 GB;
   NEVER load `esm_cache.npz` (1.5 GB) — use `s12/esm_bank.py` (compact, 100 MB) instead.
   Set `torch.set_num_threads(2)` and `OMP_NUM_THREADS=2` (4 agents share 8 cores).
   If free_gb() < 1.5, wait/poll rather than launch.
7. **Every experiment**: explicit hypothesis, control, paired comparison on the 126 targets,
   bootstrap CI, W/L, per-fold effects, concentration (drop-top-10/20), and a null control
   for any attractive correlation. `instrument.paired()` does the bookkeeping. Report
   negatives as fully as positives. Report the FAIL18 subgroup and the other-108 subgroup
   separately, always.
8. Ignore effects < ~0.03 Å or that vanish when the top 10 targets are dropped, unless
   scientifically interesting for another reason.
9. Be honest. Do not manufacture improvements. A rigorous negative is a deliverable.

## 1. What the system is (verified from code, not from summaries)

Pipeline (`core/pipeline.py`, reference `s9/final.py`):
- **retrieve**: for target of length n, every length-n window of (out-of-fold peptides from
  the 787-peptide `peptide_db` + this fold's ~6000-fragment `fragment_db` from 13,751 PDB
  proteins in `prots/`) → BLOSUM62 sum vs target → stable argsort → **K=500** pool.
  Universe ≈ 7k–40k windows/target; ~20% from peptides (`org=True`), ~80% protein fragments.
- **filter**: the shipped distogram (leave-fold-out MLP over pair features incl. ESM-2 PCA-32
  per-residue embeddings and ESM contact map; 17 distance bins; L1 Bayes-risk score) →
  **top-75** by score.
- **synthesise**: superpose the 75 on their medoid, coordinate-average (not a peptide,
  CA-CA ≈ 2.96 Å), project onto the ideal-geometry manifold by L-BFGS over (phi,psi) with a
  hinged Ramachandran prior λ=0.3, multi-start. (`core/project.py`)
- **relax**: AMBER ff14SB/GBn2 restrained minimisation k=10 (validity stage, +0.02 Å cost).
- Optional mandated components: CVaR-VQE selector over top-128 hypotheses (energy = zrank of
  distogram score; readout = p-weighted medoid/average) and Legacy 11-term knowledge-based
  energy as late refiner. Both measured null (S11-1).

Data: `pdbs/`+`pdbs_ext/` (1,463 RCSB single-chain 8–26-mers → 787 usable after rebuild
gating and dedup), identity clusters at 0.6 (Needleman-Wunsch / longer length), 5 folds
over clusters. ESM-2 650M embeddings exist for all 787 peptides + 22,007 large fragments
(`s12/esm_bank.py: load() -> {seq: (pca32 (n,32), pca128 (n,128), contacts (n,n))}`).
The ESM-2 650M weights are on disk (~2.6 GB in RAM to run — only if essential and alone).
CPU only (8 cores, 6.4 fast-core equivalents), torch 2.13 CPU, OpenMM CPU, pennylane.

## 2. The empirical record you must internalise (tuning126 unless stated)

| quantity | value |
|---|---|
| shipped distogram argmin (baseline) | 3.454 |
| synthesis (avg → project, λ=0) | 3.204 |
| synthesis + λ=0.3 + AMBER (full system) | 3.236 |
| top-75 best member (oracle) | 2.306 |
| K=500 pool best member (oracle) | 1.711 |
| whole-universe best window (oracle) | 1.313 |
| oracle re-weighting of top-75 (emitted) | 2.044 |
| convex hull of K=500 (oracle, emitted) | 0.853 |
| **benchmark60**: baseline 2.9507 vs system 2.9610, Δ +0.010 [-0.160,+0.180], 31W/29L | **no validated gain** |

Distogram: MAE ≈ 2.34 Å, r ≈ 0.70 vs true CA distances; calibration slope +0.376 (errors are
correlated across pairs and only 28% aligned with truth). The native is the objective's
argmin on 3/126 targets (percentile ≈ 37). Within the near-native band (pool_best + 1.5 Å),
the score's rank correlation with true RMSD is ≈ +0.13; every feature-based in-band ranker
tried (223 features, ridge/GBT/MLP, LFO) reaches ≈ +0.20 and buys ≈ 0 Å. The consensus
medoid is the only in-band discriminator and it is outlier-avoidance capped at the pool mode.

**The 18 zero-recall targets (FAIL18)**: the top-75 keeps NONE of the near-native band;
they return ≈ 6.0 Å from a 2.28 Å pool best. On them ρ(score, true RMSD) ≈ +0.11 vs +0.65
on the other 108. `instrument.FAIL18` lists them. Fixing them alone would move the mean by
≈ -0.5 Å. A method that does not help them is probably not solving the central problem;
a method that helps them but hurts the 108 must be reported as such.

Corrections you must respect (earlier closures that broke under audit):
- C1: the "perfect filter caps at 2.406" pool cap was measured through argmin; averaging the
  same oracle pool gives 1.925 and oracle top-25 averaging 1.644. The terminal operator
  matters enormously. Do not cite 2.406 as a cap.
- C2: "in-band discrimination is informationally impossible" rested on a capacity curve at
  fixed sample size; a learning curve (8→+0.095, 16→+0.137, 32→+0.171, 50→+0.179,
  90→+0.188, 101→+0.19) says sample-limited, not proven impossible. The untested class is
  LEARNED AGGREGATION over the full n×n candidate-vs-objective deviation map (previous
  rankers saw scalar poolings).
- C3: "operator loss = 0.741 × set mean − 1.673" is an ecological regression across 12 sets;
  within-target slope 0.236. Descriptive only; cardinality confounds it.
- C4: the hull-capacity ladder is largely a linear-algebra capacity artefact (centered
  coordinate dimension 24–45 vs m ≥ 75 candidates); donor-pool control showed most of the
  hull advantage is generic. Do not use hull numbers as achievable bounds.
- C5: the "loop is a contraction to 3.2 Å" fixed point depends on the terminal operator.
- Library size: a 3.2× larger fragment library from the same corpus does not move pool mean
  or pool best (recon_library_saturation.json); universe best moves 0.049. DEPTH is
  saturated; BREADTH/coverage of conformational space is the open question.
- S7-2: training the distogram on 3–10× more protein-fragment pairs made dev MAE
  monotonically WORSE (distribution shift: fragments in folded proteins are more compact).
- S7-11: ESM features buy −0.288 Å on selection vs one-hot (MAE said tie — MAE is not the
  metric). S7-12: BLOSUM retrieval beats random windows at every K.
- S6-6: an ORACLE secondary-structure filter (native SS vector within τ) before the argmin
  was worth only 0.21 Å as a FILTER. SS as a RETRIEVAL KEY (i.e., changing which windows
  enter the pool) is a different hypothesis and is open.
- Error-structure clue: at the SAME MAE, native + random-sign residuals through the
  projection path performs far better than the real distogram. Error correlation/sign
  structure matters more than MAE. Projection needs ≈ MAE 0.87 / r 0.94 to reach 2 Å.
- Two-piece assembly oracle floor (reconnaissance, not yet in the repo record): ≈ 0.53 Å on
  the FAIL18 and ≈ 0.32 Å on the easier set; fragment placement floor ≈ 0.36 Å. Unverified —
  re-derive before relying on it.
- No fresh cluster-disjoint benchmark exists: ALL 204 identity clusters of 9–16-mers in the
  database are already in tuning126/dev24/benchmark60. A fresh benchmark needs NEW data.

## 3. The instrument (`s12/instrument.py`) — use it

```python
from s12 import instrument as I
tg = I.targets()                # 126 dicts: pdb, n, fold, seq (pinned order)
u  = I.load_univ(pdb)           # W (nw,n,3), PHI, PSI, S, org, sim, order, rr (ORACLE), nat_ca (ORACLE)
p  = I.pool_idx(u)              # the shipped K=500 pool (indices into the universe)
rec = I.shipped_record(pdb)     # production record: sub (top-75 idx INTO THE POOL), ca, fit_ca, amber_ca, rmsd_*
dg = I.distogram(pdb)           # shipped LFO distogram: prob (npairs,17), i, j, expected, sd, risk, centres, grid
sc = I.shipped_score(dg, I.pair_dists(u["W"][p], *I.pair_index(u["n"])))   # lower = better
C, b = I.coordinate_average(W)  # S8-11 operator;  out = I.project(C, seq, fold) -> ca/phi/psi/fit_ca
I.kabsch_rmsd_batch(W, T); I.ca_rmsd(a, b); I.superpose_batch(W, T); I.pairwise_rmsd(W)
I.ss_of(phi, psi)               # simplified DSSP string from torsions
I.paired(a, b, folds=..., names=...)   # paired stats with CI, W/L, per-fold, concentration
I.write("myagent_result", obj)  # -> s12/results/myagent_result.json
I.FAIL18; I.free_gb(); I.decode(codes)
```
`python -m s12.instrument` reproduces 3.4540 / 1.7108 / 2.3062 / 3.2041 / 18 exactly.
The projection costs ~2–8 s per call; the coordinate average is instant; AMBER ~6 s per
structure (use `core.amber.refine_coords` only if you must, and sparingly).
Legacy energy: `core.energy.components_batch` / `BatchLegacy` (0.3 ms per structure).
Torsions for any window are real parent torsions; `I.build_ca(phi, psi)` rebuilds ideal CA.

## 4. Reporting

Write `s12/<yourname>_FINDINGS.md` incrementally (so work survives interruption) with: the
question, the instrument validation you ran first, each experiment (hypothesis / control /
result table / CI / W/L / FAIL18 vs other-108 / per-fold / concentration / interpretation),
negatives, leakage audit, and a final ranked list of what the coordinator should do next.
Numbers in tables, not prose. Your final message to the coordinator must be a compact
summary (≤ 60 lines) of the findings and the exact file paths.
