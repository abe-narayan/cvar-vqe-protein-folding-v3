# Sprint 12 / agent "assembly" — multi-piece fragment assembly as a candidate generator

Question. Whole-window retrieval caps at a universe best of 1.313 Å (real coordinates) and a
K=500 pool best of 1.711 Å because a 9–16-mer rarely has a single library window matching
its whole conformation. Do 3–10-residue pieces (which almost always exist locally) assembled by
torsion concatenation give (a) a lower ORACLE floor and (b) a DEPLOYABLE generator whose pool
beats retrieval, on tuning126, FAIL18 and other-108 reported separately?

Written incrementally; sections are appended as experiments finish.

## 0. Instrument validation

`python -m s12.instrument` → shipped 3.4540 / pool best 1.7108 / top-75 best 2.3062 /
synthesis 3.2041 / 18 zero-recall = FAIL18. Reproduced exactly.

## 0b. Piece banks (`s12/assembly_bank.py`, cache `s12/cache/pieces_f{fold}_{res,L3..L10}.npz`)

Library for fold f = out-of-fold peptides (`peptide_db.folds(5)`) + `core.predict._fold_fragments(f, 5)`,
i.e. exactly the universe's chain set. Per fold: ≈ 6,630 chains (≈ 630 peptides + ≈ 6,000
fragments), ≈ 92,200 residues; pieces per L: 79k (L=3) … 32.5k (L=10). Each piece carries its
parent's real torsions, the ideal-geometry CA rebuilt from them (`instrument.build_ca`), and the
parent's real CA. Peptide chain-end torsions (phi[0], psi[-1]) are placeholders and are
flagged (`valid_phi/valid_psi`); fragment ends are real. Torsion states: circular k-means k=8
fitted per fold on library residues; bigram transition table per fold.

Geometry fact that frames everything below (control, not an experiment):
**native torsions → ideal-geometry rebuild → CA-RMSD to the native = 0.347 Å mean on the 126
(FAIL18 0.474, other-108 0.326; worst 1ID6 1.474 Å).** Any torsion-concatenation assembly
inherits this floor; the reconnaissance numbers (0.32 / 0.53) sit AT this floor, i.e. they say
"the pieces exist", not "assembly reaches 0.3 Å on top of retrieval's geometry".

CA-relevant torsions: the CA trace depends only on psi[0..n-2] and phi[1..n-1]; a junction
between pieces A=[0,c) and B=[c,n) is therefore exactly the pair (psi_{c-1} from A, phi_c from
B). Under the cut convention each residue's (phi,psi) pair stays intact from one parent; the
overlap-1 conventions (A / B / circular average) change the shared residue's pair.

## E1 — ORACLE floors of torsion-concatenation assembly (`s12/assembly_e1.py`, `s12/results/assembly_e1_floors.json`)

All arms ORACLE/DIAGNOSTIC (the native picks pieces / combinations). Ideal geometry throughout;
whole-window controls given on both real and ideal geometry. Mean CA-RMSD (Å), all / FAIL18 / other-108.

Local piece existence: mean over all positions of the best bank piece's local RMSD to the native sub-trace.

| L | best (ideal) | best (real CA) | 10th best (ideal) |
|---|---|---|---|
| 3 | 0.014 / 0.014 / 0.014 | 0.003 | 0.014 |
| 4 | 0.038 / 0.042 / 0.038 | 0.036 | 0.068 |
| 5 | 0.133 / 0.157 / 0.128 | 0.130 | 0.206 |
| 6 | 0.260 / 0.309 / 0.252 | 0.258 | 0.371 |
| 7 | 0.399 / 0.491 / 0.384 | 0.394 | 0.545 |
| 8 | 0.546 / 0.664 / 0.527 | 0.539 | 0.724 |
| 9 | 0.700 / 0.841 / 0.676 | 0.689 | 0.912 |
| 10 (FAIL18 only, n≥10) | – / 0.996 / – | 0.974 | 1.331 |

Pieces exist at every length; the FAIL18 sub-traces are ≈ 20–25 % harder to match locally at every L.
Ideal vs real piece CA differ by < 0.01 Å at every L (ideal geometry is not what limits short pieces).

Assembly floors (mid cut = cut nearest n/2; best cut = oracle over all cuts with both lengths in [3,10]):

| arm | all | FAIL18 | other-108 |
|---|---|---|---|
| native torsions rebuilt (ideal-geometry floor, control) | 0.347 | 0.474 | 0.326 |
| whole-window universe best, real CA (control) | 1.313 | 1.640 | 1.259 |
| whole-window universe best, ideal rebuilt (m=1 control) | 1.313 | 1.622 | 1.262 |
| K=500 retrieval pool best, real (control) | 1.711 | 2.284 | 1.615 |
| m=2 mid cut, per-piece ARGMIN concatenated | 3.267 | 3.871 | 3.166 |
| m=2 mid cut, argmin pieces + junction (psi_b, phi_b+1) optimised ("placement floor") | 2.164 | 2.734 | 2.069 |
| m=2 mid cut, COMBINATION oracle over top-50 × top-50 (2 500 builds) | 0.841 | 1.090 | 0.800 |
| m=2 mid cut, overlap-1 avg convention, combination top-30² | 0.777 | 1.081 | 0.727 |
| m=2 best cut, per-piece argmin | 1.595 | 2.164 | 1.500 |
| m=2 best cut, junction optimised | 1.934 | 2.809 | 1.788 |
| m=2 best cut, combination top-50² | 0.690 | 0.925 | 0.651 |
| m=2 best cut, overlap-1 avg, combination top-30² | 0.651 | 0.909 | 0.608 |
| m=3 equal split, per-piece argmin | 3.871 | 4.940 | 3.693 |
| m=3 equal split, junction optimised (4 torsions) | 2.267 | 2.913 | 2.159 |
| m=3 equal split, combination top-15³ (3 375 builds) | 1.306 | 2.040 | 1.183 |
| m=3 best composition (L∈[3,8]), combination top-15³ | 0.848 | 1.016 | 0.819 |
| m=4 equal split (93 targets n≥12), per-piece argmin | 4.602 | 5.359 | 4.444 |
| m=4 equal split, junction optimised (6 torsions) | 2.511 | 2.884 | 2.433 |
| m=4 equal split, combination top-8⁴ (4 096 builds) | 2.012 | 2.404 | 1.930 |

Paired (negative = first arm better): m2-mid combo50 vs universe-ideal best −0.472 [−0.542, −0.407],
118 W / 8 L, drop-top-10 −0.399, per fold −0.585/−0.390/−0.501/−0.513/−0.385. m2-mid combo50 vs
native rebuild +0.495 [+0.398, +0.585] (the combination floor sits 0.5 Å ABOVE the ideal-geometry floor,
not at it). Fraction of targets < 1.0 Å: universe best 0.373, m2-mid combo 0.603, m2-best-cut combo 0.738.

FAIL18 individually (natreb / universe / m2-mid combo / m2-best combo / m3-best combo): 1ID6 1.47/2.66/1.81/1.67/1.54;
1JBF 0.93/1.67/0.84/0.62/0.72; 2BFI 0.13/1.35/0.24/0.21/0.34; 2JN5 0.15/1.84/0.93/0.93/1.44; 2MQ2 0.05/2.47/1.53/1.26/1.50;
2N5C 0.88/1.40/1.46/1.27/1.17; 3BTB 0.99/2.77/1.56/1.42/1.44; 7JS6 0.12/1.92/1.70/1.41/1.26; 7LCW 1.28/2.32/2.38/1.44/1.12;
8T63 0.40/2.43/1.01/1.01/0.90; 9KAR 0.16/1.40/1.18/1.18/1.28; 9L1M 0.05/1.98/1.33/1.11/2.09 (others < 1 Å).

Interpretation (error decomposition for m=2 mid cut):
* local-piece error ≈ 0.35 Å (pieces of length n/2 exist to 0.35 Å; at n/3 to 0.08, at n/4 to 0.02);
* naive concatenation of the two LOCALLY best pieces: 3.27 Å — worse than the retrieval pool best.
  Freeing the two CA-relevant junction torsions recovers only 1.1 Å (→ 2.16). A junction has 2 DOF
  (psi_b, phi_b+1; omega and bond angles fixed) for a 3-DOF relative orientation, so the pieces'
  END FRAMES must already be right; the locally-best piece's end frame usually is not;
* choosing the COMBINATION (top-50 × top-50 on global RMSD) reaches 0.84 (0.69 with cut choice, 0.65
  with overlap-averaging), i.e. the 2 500–17 500 assembled candidates contain a member 0.47–0.62 Å
  better than the best of the ~7k–40k whole windows. This is the verified form of the reconnaissance
  claim; the 0.32/0.53 numbers were not reproduced (they would be BELOW the native-rebuild floor for
  many targets and are therefore not achievable by torsion concatenation);
* more pieces do not help at fixed beam: m=3 with k=15 (3 375 combos) is 1.31, m=4 with k=8 is 2.01 —
  the junction-frame problem compounds and the oracle beam is too narrow. With the best 3-composition
  m=3 matches m=2 (0.85). Assembly's value is in 2 (at most 3) pieces of 5–8 residues, and it is a
  COMBINATORIAL selection problem, not a local-matching problem;
* remaining gap to the ideal-geometry floor: 0.84 − 0.35 ≈ 0.5 Å (global topology / beam limit).

## E2 — DEPLOYABLE piece selection (`s12/assembly_e2.py`, `s12/results/assembly_e2_selection.json`)

Hypothesis: local structure at 4–8 residues is recognisable from sequence even though whole-window
structure is not, so a deployable key should put a near-native piece in its top-k.
Control: `random` (uniform over the 32k–79k-piece bank) — the decisive control, because the metric is
min-over-top-k and a large bank has a non-trivial base rate. ORACLE = the E1 local floor.
Metric: mean over all (target, position) of the minimum local CA-RMSD among the top-k, and hit@k(τ).

LFO torsion-state classifier (8 circular k-means states, ESM pca128 + ±2 one-hot neighbours,
per-fold, chain-held-out validation): library holdout accuracy 0.606–0.636 (top-2 0.77–0.79) vs
majority 0.39–0.43. On the tuning targets' NATIVE torsions (ORACLE labels): 0.29–0.50, mean 0.41
(top-2 0.48–0.66) vs majority 0.30–0.43 — i.e. the local-state signal is real on library residues
and roughly HALF as strong on the targets, barely above the majority-state rate on 2 of 5 folds.
Training on peptides only (8.4k residues) loses library accuracy (0.48–0.58) but gains on targets
(0.35–0.47): the fragment residues are a distribution shift for local state as they were for the
distogram (S7-2).

min local RMSD in top-20, all / FAIL18 (lower is better); hit@20(≤1.0 Å):

| key | L=4 | L=5 | L=6 | L=7 | L=8 |
|---|---|---|---|---|---|
| ORACLE (E1 floor) | 0.038 / 0.042 | 0.133 / 0.157 | 0.260 / 0.309 | 0.399 / 0.491 | 0.546 / 0.664 |
| **random (control)** | **0.343 / 0.421** (h1 .986) | **0.653 / 0.770** (.771) | **0.962 / 1.139** (.523) | **1.215 / 1.451** (.407) | **1.488 / 1.798** (.320) |
| esm | 0.348 / 0.457 (.959) | 0.655 / 0.843 (.775) | 0.908 / 1.205 (.585) | 1.164 / 1.535 (.469) | 1.391 / 1.913 (.408) |
| blosum | 0.402 / 0.555 (.941) | 0.731 / 0.962 (.713) | 1.006 / 1.316 (.510) | 1.248 / 1.588 (.424) | 1.450 / 2.006 (.373) |
| blosum+esm | 0.392 / 0.530 | 0.717 / 0.965 | 0.968 / 1.281 | 1.239 / 1.592 | 1.460 / 2.077 |
| lfo (state classifier) | 0.671 / 0.977 | 1.044 / 1.583 | 1.319 / 2.053 | 1.554 / 2.477 | 1.793 / 2.929 |
| lfo_pep | 0.671 / 1.072 | 1.057 / 1.660 | 1.359 / 2.109 | 1.617 / 2.507 | 1.838 / 2.834 |
| tors (torsion regression) | 0.756 / 1.122 | 1.127 / 1.689 | 1.405 / 2.135 | 1.650 / 2.557 | 1.868 / 2.907 |
| combo (z blosum+esm+lfo+tors) | 0.579 / 0.891 | 0.908 / 1.427 | 1.163 / 1.852 | 1.384 / 2.196 | 1.610 / 2.622 |

Paired vs the random control (positive = key WORSE than random), k=20:

| comparison | L=5 | L=6 | L=8 |
|---|---|---|---|
| blosum − random | +0.078 [+0.044,+0.113] 42W/84L | +0.043 [−0.003,+0.091] | −0.038 [−0.102,+0.028] |
| esm − random | +0.002 [−0.033,+0.040] | **−0.054 [−0.100,−0.006]** 78W/48L | **−0.097 [−0.166,−0.027]** 78W/48L |
| lfo_pep − random | +0.404 [+0.306,+0.505] | +0.397 [+0.285,+0.514] | +0.350 [+0.220,+0.485] |
| combo − random | +0.256 [+0.176,+0.338] | +0.201 [+0.109,+0.298] | +0.122 [−0.001,+0.247] |
| esm − blosum | −0.076 [−0.104,−0.049] 92W/34L | −0.098 [−0.130,−0.067] 93W/33L | −0.059 [−0.108,−0.012] |

**Result: negative, and the negative is the finding.** No deployable key beats a uniformly random
draw from the piece bank by more than 0.10 Å at any length. ESM per-residue similarity is the only
key with a CI excluding zero against random (−0.054 at L=6, −0.097 at L=8, ~0.03 after dropping the
top 10 targets), consistent with S7-11 (ESM > one-hot on selection). BLOSUM is at or WORSE than
random for local pieces — the opposite of its behaviour as a whole-window retrieval key (S7-12) —
and the learned local-structure predictors are much worse than random.

Mechanism of the "worse than random" result (not a bug): the metric is min-over-top-k, so it rewards
DIVERSITY as much as accuracy. The LFO/torsion keys rank by agreement with a single predicted state
sequence, so their top-20 are ~20 copies of one conformation; a random 20 spans the bank's
conformational distribution. Concretely: at L=6 the base rate of a random piece being within 1.0 Å is
0.122, so 20 random draws hit 0.521 of the time and that is already most of what any key achieves
(esm 0.585). The bank is dense enough at short length that coverage, not ranking, supplies the hits.
Corollary for E3: a top-k-by-key piece shortlist is nearly a random shortlist, so the ASSEMBLY
objective must do all the work, and diversity in the shortlist is worth more than key quality.

Leakage audit for E2: all keys use only (target sequence, target ESM from the bank of all 787
peptides, out-of-fold library torsions/states). The LFO is trained per fold on out-of-fold chains
only; native torsions appear solely in `native_state_accuracy` (labelled ORACLE) and in the
evaluation metric.

## E3 — DEPLOYABLE assembly (`s12/assembly_e3.py`, `s12/results/assembly_e3_combo_k20.json`)

Setup: key = `combo`, k = 20 pieces per position, all 2-piece cuts (L ∈ [4,10]) and all 3-piece
compositions (L ∈ [4,8]; [3,8] for n < 12) — mean 50,714 assembled candidates per target, each built
by torsion concatenation. Objectives: shipped distogram Bayes risk (`dist`); + junction bigram
−log P(state|state) (`junc`); + |rg − rg_pool_median| ; + Legacy total energy on the ideal backbone
(top-2000 by dist). Terminal operator = the shipped one (top-75 by dist → coordinate_average → project).

| quantity | all | FAIL18 | other-108 |
|---|---|---|---|
| assembled space best (ORACLE eval of the DEPLOYABLE space) | 1.854 | 3.030 | 1.657 |
| … 2-piece structures only | 2.073 | 3.457 | 1.842 |
| … 3-piece structures only | 1.970 | 3.175 | 1.769 |
| retrieval K=500 pool best (control) | 1.711 | 2.284 | 1.615 |
| whole-universe best (control) | 1.313 | 1.640 | 1.259 |
| assembled top-500-by-dist pool best | 2.689 | 5.188 | 2.272 |
| assembled argmin | 3.384 | 5.790 | 2.984 |
| shipped argmin (control) | 3.454 | 6.008 | 3.028 |
| assembled avg → project | 3.359 | 5.949 | 2.928 |
| hybrid (retrieval 500 ∪ assembled 500) → project | 3.330 | 5.949 | 2.893 |
| shipped synthesis (control) | 3.204 | 6.026 | 2.734 |

Paired (negative = assembly better): assembled space best vs universe best **+0.540 [+0.406,+0.684]**,
30W/96L. Assembled pool best vs retrieval pool best **+0.978 [+0.792,+1.178]**, 19W/107L.
Assembled argmin vs shipped argmin −0.070 [−0.187,+0.044] 64W/62L (drop-top-10 +0.059 — a wash).
Assembled projection vs shipped synthesis +0.155 [+0.047,+0.270] 47W/79L; hybrid +0.126 [+0.022,+0.233].
Subgroups: on the other-108 the assembled projection is WORSE (+0.194 [+0.080,+0.317], 39W/69L); on the
FAIL18 it is nominally better (−0.077 [−0.395,+0.195], 8W/10L — CI includes zero, no effect).
`dist+junc`, `+rg` and `+legacy` all move the argmin by < 0.13 Å and none beats `dist` on projection.
Error correlation between the assembled answer and the shipped answer across the 126: **r = 0.928**
(hybrid 0.936); an ORACLE pick of the better of the two per target gives 3.050 vs 3.204 — the two
systems fail on the same targets, so assembly is not a decorrelated second opinion.

**Result: the deployable generator is not better than retrieval, and the deployable answer is not
better than the shipped one.** Three separate losses stack up:
1. the deployable SHORTLIST costs the oracle advantage. With the oracle top-50 per piece the m=2
   combination floor is 0.841 (E1); with the deployable top-20 by `combo` the best in a 50k-candidate
   assembled space is 1.854 — worse than the 1.711 retrieval pool it was meant to replace. This is the
   direct consequence of E2 (keys ≈ random for local pieces) and it hits FAIL18 hardest
   (gen best vs retrieval pool best +0.746 [+0.249,+1.267] on FAIL18, +0.042 on the other 108);
2. the distogram cannot rank inside the assembled space: the assembled pool's top-500 best is 2.689
   against an in-space best of 1.854, i.e. selection loses 0.84 Å on top of generation;
3. **the assembled candidates GAME the objective**: in the hybrid pool, 91.7 % of the distogram top-75
   is assembled rather than retrieved, yet the hybrid answer is no better than the shipped one. Torsion
   concatenation can drive the Bayes-risk score below anything retrieval offers while moving away from
   the native — a concrete instance of the "optimise harder, get worse" record (objective-does-not-rank-the-native).

Solvers (1,402 sub-problems, mean 4,558 feasible assignments, objective `dist+junc`):
simulated annealing (2,000 steps) finds the exhaustive optimum on **90.5 %** of problems (mean rank gap
0.28); left-to-right beam with width 20 on a prefix-distogram surrogate finds it on only **44.2 %**
(mean rank gap 580). The combinatorial problem is easy for a classical stochastic solver at this size —
solver quality is NOT the bottleneck, and any quantum-advantage claim on this instance family must be
measured against SA at 2,000 steps, not against beam search.

## E5 — ranking inside the assembled space, and the DIVERSITY correction (`s12/assembly_e5.py`, `s12/results/assembly_e5_ranking.json`)

E2 predicted that a top-k-by-key shortlist is worse than a diverse one because the metric that matters
is min-over-shortlist. E5 tests that at the ASSEMBLY level: the same enumeration as E3 (mean 50,714
candidates/target) with three deployable shortlist rules — `combo`, `esm`, and `random` (a fixed
seeded draw from the bank; deployable, uses nothing about the target).

| | retrieval K=500 (control) | combo k=20 | esm k=20 | **random k=20** |
|---|---|---|---|---|
| space best (all) | 1.711 | 1.854 | **1.447** | **1.456** |
| space best (FAIL18) | 2.284 | 3.030 | 1.974 | **1.854** |
| space best (other-108) | 1.615 | 1.657 | **1.360** | 1.390 |
| best at k=4 / k=8 / k=20 | — | 2.383 / 2.158 / 1.854 | 2.134 / 1.813 / 1.447 | 2.098 / 1.765 / 1.456 |
| mean pairwise RMSD of the scored top-40 (diversity) | — | 1.074 | 1.536 | 1.903 |
| ρ(score, true RMSD) global | 0.568 | 0.375 | 0.531 | 0.562 |
| ρ(score, true RMSD) in-band | 0.126 | 0.079 | 0.075 | 0.185 |
| percentile of the NATIVE's own score | 36.5 | 56.8 | 37.5 | **23.6** |
| top-75 keeps a band member (recall) | 0.857 | 0.754 | 0.635 | 0.619 |
| top-75 best | 2.306 | 2.911 | 2.853 | 2.753 |
| argmin | 3.454 | 3.384 | 3.507 | 3.531 |

Paired vs the retrieval pool best: combo +0.143 [+0.020,+0.272] 65W/61L; **esm −0.263 [−0.344,−0.187]
95W/31L; random −0.255 [−0.340,−0.168] 95W/31L.**

**This is the sprint-relevant positive: two-and-three-piece assembly with a DIVERSE deployable
shortlist generates a candidate space whose best member is 0.26 Å better than the shipped K=500
retrieval pool (1.45 vs 1.71), and 0.43 Å better on the FAIL18 (1.85 vs 2.28), using no oracle.**
It is also 0.14 Å better than the whole-universe best window on the FAIL18 (1.854 vs 1.640 — no, WORSE
there; against the K=500 pool it wins, against the full universe it does not: universe best is 1.313).
So assembly beats the pool the pipeline actually uses, not the retrieval universe as a whole.

And the reason it does not turn into accuracy: the shipped score cannot rank inside it. Top-75 best
2.75 vs the retrieval pool's 2.306, band recall drops 0.857 → 0.619, and the argmin does not move.
The `combo`-shortlisted space is the worst on every ranking measure AND the one where the native sits
at the 56.8th percentile of the objective — the key that "looks" most sequence-informed produces the
space in which the objective is most misleading. Ranking skill in-band is ~0.13–0.19 everywhere, i.e.
unchanged from the record; assembly does not create a better-behaved objective landscape.

Union of the two diverse spaces (esm ∪ random, still deployable): best 1.294 mean
(FAIL18 1.713, other-108 1.224) — equal to the whole-universe best window (1.313) at 1/1000 of the
enumeration, and better than it on the other-108. Union with the retrieval pool: 1.262.

### E3 with the diverse shortlists (`assembly_e3_random_k20.json`, `assembly_e3_esm_k20.json`)

The full deployable chain re-run with the `random` and `esm` shortlists that E5 showed generate the
better space. Mean CA-RMSD, all / FAIL18 / other-108:

| arm | combo k=20 | random k=20 | esm k=20 | shipped control |
|---|---|---|---|---|
| assembled space best (ORACLE eval) | 1.854 / 3.030 / 1.657 | **1.438 / 1.839 / 1.371** | 1.447 / 1.974 / 1.360 | pool best 1.711 / 2.284 / 1.615 |
| … 2-piece only / 3-piece only | 2.073 / 1.970 | 1.741 / 1.495 | — | — |
| top-500-by-dist pool best | 2.689 | 2.321 | 2.502 | 1.711 |
| top-75-by-dist best | 2.911 | 2.684 | 2.853 | 2.306 |
| argmin | 3.384 | 3.551 | 3.507 | 3.454 |
| coordinate average → project | 3.359 | 3.397 | 3.410 | 3.204 |
| hybrid (∪ retrieval 500) → project | 3.330 | 3.385 | 3.406 | 3.204 |
| hybrid pool best | 1.686 | **1.631** | 1.636 | 1.711 |
| fraction of the hybrid top-75 that is assembled | 0.917 | 0.839 | 0.924 | — |

Paired vs shipped: argmin +0.097 [−0.013,+0.208] (random), +0.053 [−0.061,+0.168] (esm);
projection +0.193 [+0.087,+0.307] (random), +0.206 [+0.082,+0.331] (esm). Every arm is a wash or a
loss at the emitted structure while the space it draws from is 0.27 Å better. Adding the junction
bigram, rg and Legacy to the objective changes the pool best by up to 0.19 Å (random: 2.321 → 2.277
with `+junc`) but never the emitted answer.

Solver numbers are unchanged across shortlists (SA 88–90 % exact, beam-20 38–44 %).

## E6 (ORACLE/DIAGNOSTIC) — what the assembled space is worth if selection were solved (`s12/assembly_e6.py`, `s12/results/assembly_e6_oracle_random.json`)

The same DEPLOYABLE `random`-shortlist assembled space (mean 50,714 candidates), selected by an ORACLE
(true CA-RMSD) top-25 / top-75, then through the shipped terminal operator (coordinate_average →
project). Control: the identical procedure on the retrieval K=500 pool, which reproduces the C1 record
(oracle top-25 averaging 1.644, oracle pool averaging ≈1.93 → here 1.963 at m=75).

| | assembled (random k=20) | retrieval K=500 | paired diff |
|---|---|---|---|
| space best | 1.464 / 1.955 / 1.382 | 1.711 / 2.284 / 1.615 | — |
| ORACLE top-25 → coordinate average | 1.203 / 1.560 / 1.143 | 1.644 / 2.416 / 1.515 | **−0.441 [−0.529,−0.358] 112W/14L** |
| ORACLE top-25 → project (emitted) | **1.265 / 1.714 / 1.190** | 1.760 / 2.623 / 1.616 | **−0.495 [−0.594,−0.401] 107W/19L** |
| ORACLE top-75 → coordinate average | 1.249 / 1.629 / 1.186 | 1.963 / 3.069 / 1.779 | −0.714 [−0.826,−0.607] 122W/4L |
| ORACLE top-75 → project (emitted) | **1.309 / 1.782 / 1.230** | 2.102 / 3.283 / 1.905 | **−0.793 [−0.917,−0.673] 118W/8L** |
(all / FAIL18 / other-108)

**This is the strongest number in this report: a solved selector over the assembled space emits
1.27–1.31 Å — comfortably under the 2.0 Å objective, and on the FAIL18 1.71–1.78 Å where the same
solved selector over the retrieval pool reaches only 2.62–3.28 Å.** The assembled space raises the
achievable ceiling by 0.50–0.79 Å, with 107–122 wins out of 126, and the effect survives dropping the
top 10 targets. The generator is therefore NOT the bottleneck any more; selection is, and it is now
the only bottleneck: E3 shows the shipped objective converts a 1.46 Å space into a 3.40 Å answer.

## E4 — QUBO / Ising form for the VQE/CVaR agent (`s12/assembly_e4.py`, `s12/results/assembly_e4_qubo.json`)

One-hot encoding x[p,i] ∈ {0,1} over k candidate pieces at each of m positions, N = m·k binary
variables, index p·k + i. Energy

    E(x) = Σ_p Σ_i h[p,i] x[p,i] + Σ_{p<q} Σ_{i,j} J[(p,i),(q,j)] x[p,i] x[q,j] + A Σ_p (Σ_i x[p,i] − 1)²

* h[p,i] = shipped distogram Bayes risk summed over the residue pairs INSIDE piece i at position p
  (exact — depends only on that piece's torsions);
* J[(p,i),(q,j)] for ADJACENT q = p+1 = the risk summed over the cross pairs of the exactly built
  two-piece sub-assembly + W_J·(−log P(state_first(q,j) | state_last(p,i))) from the fold's bigram;
* J for NON-adjacent (p,q) (m = 3 only) is the MEAN over the middle piece's candidates — the cross
  distances genuinely depend on the middle piece, so the m=3 QUBO is a mean-field approximation and
  the m=2 QUBO is exact;
* A = 2(max|h| + max row-sum|J|), so every infeasible bitstring costs more than every feasible one.
Ising: x = (1 − s)/2. Stored per instance: `Q` (N×N upper-triangular), `const`, `A`, `comp`,
`piece_start`, `piece_bank_idx`, plus the full `energies_qubo`, `energies_true` and ORACLE `rr` vectors
over all k^m feasible assignments.

Files: `s12/cache/qubo_<pdb>_m{2,3}_k{4,8}.npz` (126 targets × 4 = 504 instances, `combo` shortlist)
and `..._random.npz` (the diverse shortlist E5/E6 showed is the better space).

| instance | N | QUBO opt = true opt | corr(E_qubo, E_true) | ORACLE rr of the QUBO optimum | best rr in the instance |
|---|---|---|---|---|---|
| m=2, k=4 | 8 | 100 % | 1.000 (exact) | 3.468 | 3.025 |
| m=2, k=8 | 16 | 100 % | 1.000 (exact) | 3.492 | 2.711 |
| m=3, k=4 | 12 | 46.0 % | 0.832 | 3.570 | 3.008 |
| m=3, k=8 | 24 | 32.5 % | 0.830 | 3.571 | 2.659 |

The classical optimum is given per instance (`opt_qubo_energy`, exhaustive over k^m). For every
N ≤ 16 instance a brute-force sweep over all 2^N bitstrings confirms the penalty weight: the global
minimum is feasible and equals the feasible optimum on 100 % of instances. **A quantum solver on this
family must be judged against the exhaustive optimum, and note that simulated annealing at 2,000 steps
already finds the optimum of the much larger full problem (mean 4,558 assignments) 90 % of the time
(E3).** Note also that reaching the QUBO optimum is not the same as reaching a good structure: the
optimum's ORACLE RMSD (≈3.47–3.57 Å) is worse than the best assignment in the same instance
(2.66–3.03 Å), i.e. the objective — not the solver — is what limits this family.

## Leakage audit

* Piece banks are built per fold from out-of-fold peptides (`peptide_db.folds(5)[seq] != fold`) plus
  `core.predict._fold_fragments(fold, 5)` — the same chain set the shipped universe uses, so the
  leakage rules are identical to the pipeline's. Verified: no chain with a target's own sequence
  appears in that target's bank (126/126); maximum NW identity of any peptide-bank chain to the target
  is 0.31–0.46 on a sample of 8 (the 0.6 clustering threshold is respected by construction).
* The LFO local-state classifier is trained per fold on out-of-fold chains only; the ESM embeddings are
  the shipped bank (all 787 peptides, target included by design — the target's own SEQUENCE embedding
  is a deployable input, as in the shipped distogram).
* Native torsions / `nat_ca` / `rr` enter only: E1 (all arms labelled ORACLE), the `native_state_accuracy`
  diagnostic, E6 (labelled ORACLE), and evaluation metrics everywhere else. No inference-time decision
  in E2/E3/E4/E5 reads a native. The `random` shortlist is a fixed seeded draw keyed on (pdb, s, L) and
  uses nothing about the structure.
* benchmark60 was never read; no `s9/final_cache` access; dev24 not touched.

The `random`-shortlist instance family (`qubo_<pdb>_m{2,3}_k{4,8}_random.npz`,
`assembly_e4_qubo_random.json`) behaves identically as an optimisation problem (m=2 exact, m=3
corr 0.865, penalty verified) but contains better structures (best rr in the m=3 k=8 instance 2.233 vs
2.659), and its QUBO optimum is *further* from the native (3.866) — the same objective/structure
divergence, sharper. This is the family to hand to a quantum solver, because on it the gap between
"finds the optimum" and "finds a good structure" is largest and therefore most informative.

## Ranked recommendations for the coordinator

1. **The generator question is answered and the answer is positive, but only with a DIVERSE shortlist.**
   Two/three-piece torsion-concatenation assembly with a random or ESM-similarity shortlist of 20
   pieces per position produces ~50k candidates per target whose best member is 1.44–1.46 Å
   (FAIL18 1.85–1.97) against the shipped pool's 1.711 (FAIL18 2.284), −0.26 [−0.34,−0.19] paired,
   95W/31L; and whose ORACLE top-25 → project emits **1.265 Å (FAIL18 1.714)** against retrieval's
   1.760 (FAIL18 2.623). Assembly moves the CEILING below 2.0 Å, including on the FAIL18.
2. **The single most valuable next experiment**: build an in-band ranker for the ASSEMBLED space
   specifically. Every existing selection result (in-band ρ ≈ 0.13–0.20, consensus medoid −0.172 Å)
   was measured on retrieval pools whose members are whole library windows; the assembled space is
   structurally different (mean pairwise RMSD of the scored top-40 is 1.90 Å vs retrieval's tighter
   cluster, band recall 0.62, native at the 23.6th score percentile — the BEST native percentile of
   any space measured in this sprint, better than retrieval's 36.5). Concretely: run the consensus
   medoid + score-filter arm (`consensus-is-the-only-in-band-discriminator`) and the learned-aggregation
   arm (C2's untested class, the n×n candidate-vs-objective deviation map) ON THE ASSEMBLED SPACE.
   The prize is the 1.27 Å ceiling; the current selector converts it to 3.40.
3. **Do not spend more effort on piece-selection keys.** BLOSUM, ESM, a learned torsion-state
   classifier and their combinations are all within 0.1 Å of a random draw for local pieces
   (ESM the only one significantly better, −0.05 to −0.10), and the ranked keys are actively harmful
   because they collapse shortlist diversity. Diversity is the design variable; the natural next step
   there is an explicit diversity-maximising shortlist (e.g. k-medoids over piece conformations)
   rather than a better score.
4. **Warning for any objective work**: assembled candidates GAME the distogram — they occupy 84–92 %
   of the hybrid distogram top-75 and drive the Bayes risk below anything retrieval offers while being
   no closer to the native, and in the `combo` space the native sits at the 56.8th percentile of the
   objective. Any future objective must be validated against a generator that can exploit it.
5. **For the VQE/CVaR agent**: 504 + 504 instances are cached (see E4). The m=2 QUBO is exact; the m=3
   is mean-field (corr 0.83–0.87). Classical baselines to beat: exhaustive optimum (given per instance)
   and SA-2000, which solves 90 % of the far larger full problems exactly. A quantum result on this
   family is only interesting as a solver benchmark — the objective, not the search, caps the accuracy.
6. **Not worth pursuing**: junction-torsion optimisation as a repair step (frees only 1.1 Å of the
   3.27 Å concatenation error and still lands at 2.16 Å, worse than just choosing the right combination);
   m ≥ 4 pieces (2.01 Å even with an oracle combination search); and the `+junc`/`+rg`/`+legacy`
   objective terms (≤ 0.19 Å on the pool best, 0 Å on the emitted answer).
