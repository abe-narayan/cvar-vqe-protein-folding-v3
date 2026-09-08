# forensics — FAIL18 forensics and structure-aware retrieval (Sprint 12)

Agent: forensics.  Instrument: tuning126 (`s8/generate_univ`).  benchmark60 and dev24 untouched.
Code: `s12/forensics_lib.py`, `s12/forensics_build.py`, `s12/forensics_part1.py`,
`s12/forensics_alphabet.py`, `s12/forensics_part2.py`.  Results: `s12/results/forensics_*.json`.

## 0. Question

(1) Why do the 18 filter-blind targets (`instrument.FAIL18`) fail — query, filter, multimodality,
SS misprediction, or the reference itself — and what is the common causal structure?
(2) Is the information needed to recognise the correct structure present in the library, and under
which retrieval KEY (BLOSUM, ESM, predicted SS, structural alphabet, distogram score) could it
have been retrieved?  Every key is compared at matched K=500 through the full production chain
(shipped-score top-75 -> coordinate average -> projection) with paired CIs, FAIL18 vs other-108.

## 1. Instrument validation

`python -m s12.instrument`: shipped 3.4540 / pool_best 1.7108 / top75_best 2.3062 / synthesis 3.2041 /
18 zero-recall = FAIL18.  Reproduced exactly.

`forensics_lib.chain` (pool -> shipped score -> top-75 -> coordinate_average -> project) was checked
against the production record on 1ID6 and 1A13: the top-75 index set is identical, `rmsd_avg` agrees to
4 decimals, `rmsd_arm`/`rmsd_fit` agree to <1e-3 (multi-start projection noise).

Parent map: every universe window was mapped back to `(parent, start)` by replaying the generation order
(out-of-fold peptides in `peptide_db.load()` order, then `core.data.fold_fragments(fold)`), verified by
exact equality of the decoded codes for EVERY window and CA coordinates on a 400-window sample per target.
The target's own sequence never appears in its universe (checked).

SS instrument: the torsion-rebuilt DSSP (`instrument.ss_of`) agrees with a CA-geometry assigner
(`forensics_lib.ca_ss`, P-SEA-like) on only 62% of native residues and calls the 1JBF hairpin and the
5W52 steric-zipper strand "coil"; `ca_ss` is used for all SS comparisons (natives and windows alike).

## 2. Part 1 — metadata of the 18 (PDB header text only: HEADER/TITLE/KEYWDS/EXPDTA/REMARK 210/SSBOND/LINK)

| pdb | n | fold | seq | expdta | title | flags |
|---|---|---|---|---|---|---|
| 1ID6 | 15 | 3 | SVQARWEAAFDLDLY | SOLUTION NMR | SOLUTION STRUCTURES OF SYR6 (100% DMSO) | cosolvent (DMSO in REMARK 210 SAMPLE) |
| 1JBF | 15 | 1 | NLPRCTEGPWGWVCM | SOLUTION NMR | HAIRPIN PEPTIDE THAT INHIBITS IGE ACTIVITY BY BINDING TO THE HIGH AFFI | bound, SSBOND 5-14 |
| 1LB7 | 16 | 4 | RNCFESVAALRRCMYG | SOLUTION NMR | IGF-F1-1, A PEPTIDE ANTAGONIST OF IGF-1 (loop-helix) | designed, SSBOND 3-13 |
| 2BFI | 12 | 2 | KFFEAAAKKFFE | X-RAY | MOLECULAR BASIS FOR AMYLOID FIBRIL FORMATION AND STABILITY | fibril, xray, designed |
| 2BP4 | 16 | 1 | DAEFRHDSGYEVHHQK | SOLUTION NMR | ZINC-BINDING DOMAIN OF ALZHEIMER'S A-BETA IN TFE-WATER (80-20) | fibril, cosolvent (80% TFE) |
| 2JN5 | 12 | 1 | MDVFMKGLSKAK | SOLUTION NMR | DODECAPEPTIDE FROM ALPHA-SYNUCLEIN BOUND WITH SYNPHILIN-1 | fibril, membrane, bound |
| 2MQ2 | 14 | 4 | RGGRLYRRRFVVGR | SOLUTION NMR | CYSTEINE DELETED PROTEGRIN-1 (CDP-1) | membrane (hairpin without its disulfides) |
| 2N5C | 15 | 4 | GFGSKPLDSFGLNFF | SOLUTION NMR | LASSO PEPTIDE CHAXAPEPTIN (100% DMSO) | lasso (LINK N-G1 : CG-D8), cosolvent |
| 2NB7 | 14 | 2 | MENTSITIEFSSKF | SOLUTION NMR | N-TERMINAL EXTRAMEMBRANE DOMAIN OF SH PROTEIN (bicelles) | membrane |
| 2NDM | 13 | 4 | GRPCYTLQSCFPD | SOLUTION NMR | PAWS DERIVED PEPTIDE 21 (PDP-21) | cyclic (LINK N-G1 : C-D13), SSBOND 4-10 |
| 3BTB | 15 | 3 | MEELQDDYEDMMEEN | SOLUTION NMR | BAND 3 PEPTIDE INHIBITOR BOUND TO G3PDH (exchange-transferred NOE) | membrane, bound |
| 3SGO | 11 | 1 | KVKVLGDVIEV | X-RAY | AMYLOID-RELATED SEGMENT OF ALPHAB-CRYSTALLIN 90-100 | fibril, xray |
| 5W52 | 11 | 1 | DLIIKGISVHI | ELECTRON CRYST | MICROED STRUCTURE OF THE SEGMENT DLIIKGISVHI FROM TDP-43 (steric zipper) | fibril, xray |
| 7JS6 | 15 | 4 | LLGRSGNDRLILSKN | SOLUTION NMR | DES-CITRULASSIN F | lasso (LINK N-L1 : CG-D8) |
| 7LCW | 15 | 4 | GSKYSDTADESSYRW | SOLUTION NMR | ASPARTIMIDYLATED LASSO PEPTIDE LIHUANODIN | lasso |
| 8T63 | 16 | 3 | WHMWNTVPNAKQVIAA | SOLUTION NMR | DESIGNED PEPTIDE PH1 | designed |
| 9KAR | 15 | 3 | GGWGTVPDWFFNMNW | SOLUTION NMR | LASSO PEPTIDE STREPTACIDIN (DMSO) | lasso, cosolvent |
| 9L1M | 12 | 1 | MMMKPDMMMKPD | SOLUTION NMR | MMMKPD2 - LOW COMPLEXITY REGION OF LMP (fibril) | fibril |

Structural alphabet codebook (K=8, k-means on fragment residues only), (phi, psi) centres in degrees:
(-143,160) beta; (-116,75); (-114,129) beta/PPII; (-96,-4) bridge; (-73,144) PPII; (-64,-39) alphaR;
(74,18) alphaL; (82,-160).

## 3. Part 1 — what separates the 18 from the 108 (`forensics_part1.json`, `forensics_causal.json`)

Mann-Whitney across all 126.  `o_` = oracle-derived (evaluation only), `dep_` = deployable,
`meta_` = PDB header metadata.  `auc_fail` = P(a FAIL18 target ranks above an other-108 target).

| quantity | FAIL18 | other-108 | p | auc_fail |
|---|---|---|---|---|
| o_band_score_pct (score percentile of the in-band windows, 0 = best) | 0.687 | 0.287 | 3.4e-10 | 0.963 |
| o_band_score_pct_min (BEST in-band window's percentile) | 0.391 | 0.009 | 3.6e-13 | — |
| rho(score, true RMSD) within the pool | +0.107 | +0.645 | 1.3e-07 | 0.111 |
| o_library_SS_agreement of the K=500 pool with the native SS | 0.319 | 0.443 | 5.4e-06 | 0.164 |
| **meta_lasso_or_fibril (header text only)** | **0.556** | **0.056** | **4.4e-09** | **0.750** |
| o_native strand fraction (ca_ss E) | 0.470 | 0.166 | 1.2e-04 | 0.742 |
| o_pool_best | 2.284 | 1.615 | 4.9e-03 | 0.708 |
| o_universe_best | 1.640 | 1.259 | 3.7e-02 | 0.654 |
| o_windows within 2.0 A (universe / pool) | 146 / 7.8 | 1341 / 62.3 | 1.6e-02 / 1.7e-02 | 0.322 |
| dep_pool_nclusters (avg-linkage 2.5 A) | 159.8 | 119.5 | 1.1e-02 | 0.688 |
| dep_disto_rg (distogram-implied radius of gyration) | 7.42 | 6.76 | 2.8e-02 | 0.662 |
| top-75 SS agreement with native | 0.242 | 0.503 | 1.3e-05 | — |
| top-75 members in a near-native pool cluster | 8.8 | 47.3 | 1.2e-07 | — |
| **near-native windows' BLOSUM rank (median)** | **6715** | **8463** | **0.169** | — |
| **near-native windows inside K=500 (fraction)** | **0.079** | **0.073** | **0.964** | — |
| near-native windows: sequence identity to target | 0.066 | 0.071 | 0.116 | — |
| near-native windows: fraction from peptides (`org`) | 0.346 | 0.269 | 0.591 | — |

### 3.0 The 18, target by target (`forensics_causal.json:fail18_detail`)

| pdb | n | native SS | distogram SS | pool best | univ best | <=2A univ/pool | band in pool | band score pct mean/best | band BLOSUM rank med | top75 SS agree | top75 SS mode | top75 best | rho(score,rr) | rg nat/top75/disto |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1ID6 | 15 | CCCCCCCCCCCCCCC | EEECHHHHHCCCCCC | 3.54 | 2.66 | 0/0 | 112 | 0.52/0.19 | 4948 | 0.06 | HHHHHHHHHHHHHHH | 5.21 | +0.36 | 5.7/6.9/7.2 |
| 1JBF | 15 | CCEEEEECCEEEEEC | EEECEEECCEEECCC | 1.67 | 1.67 | 9/3 | 24 | 0.52/0.23 | 2611 | 0.43 | CCCCCCCCCCCCCCC | 3.49 | +0.47 | 7.3/6.4/6.5 |
| 1LB7 | 16 | CCCCCHHHHHHHHHHH | CCCHHHHHHHHHHHCC | 2.01 | 0.77 | 6/0 | 53 | 0.54/0.33 | 2749 | 0.66 | HHHHHHHHHHHHHHHH | 4.45 | +0.27 | 6.1/7.3/7.4 |
| 2BFI | 12 | EEEEEEEEEECC | HHHHHHHHHHHH | 2.22 | 1.35 | 121/0 | 9 | 0.99/0.97 | 11624 | 0.01 | HHHHHHHHHHHH | 6.92 | +0.02 | 11.9/5.7/5.8 |
| 2BP4 | 16 | HHHHHHHHHHHHHHHC | CCCCEEECEEECCCCC | 0.63 | 0.50 | 1245/91 | 94 | 0.51/0.23 | 3228 | 0.24 | CCCCCCCCCCCCCCCC | 3.07 | +0.10 | 7.3/7.6/8.5 |
| 2JN5 | 12 | CCCCCCCCCCCC | HHHHHHCCCCCC | 2.69 | 1.84 | 1/0 | 94 | 0.67/0.35 | 11574 | 0.06 | HHHHHHHHHHHH | 5.39 | -0.06 | 5.3/5.7/6.5 |
| 2MQ2 | 14 | CCCCCCCCEEEEEC | EEEEEEEEECEEEC | 3.40 | 2.47 | 0/0 | 72 | 0.65/0.21 | 7200 | 0.47 | CCEEEEEEEEEECC | 5.49 | -0.51 | 5.4/10.3/10.1 |
| 2N5C | 15 | CCCCEEEEECEEEEE | CCCEEECCCCCCCCC | 2.46 | 1.40 | 6/0 | 7 | 0.54/0.28 | 2949 | 0.24 | CCCCCCCCCCCCCCC | 4.88 | +0.47 | 5.3/6.7/7.1 |
| 2NB7 | 14 | CCCCHHHHHHHHHH | CCEEEEEECCCCCC | 1.28 | 1.12 | 303/34 | 66 | 0.63/0.29 | 5268 | 0.20 | CCCEEEEECCCCCC | 3.72 | -0.37 | 6.3/9.4/9.3 |
| 2NDM | 13 | EEEEECCCEEEEE | EEECCCHHHHHCC | 1.65 | 1.39 | 76/1 | 16 | 0.63/0.39 | 9149 | 0.05 | HHHHHHHHHHHHH | 4.00 | +0.26 | 6.5/6.2/6.9 |
| 3BTB | 15 | CCEEEEECCEEEEEC | CHHHHHHHHHHHHHH | 2.97 | 2.77 | 0/0 | 24 | 0.66/0.29 | 4936 | 0.03 | HHHHHHHHHHHHHHH | 4.82 | +0.44 | 7.6/6.9/6.6 |
| 3SGO | 11 | EEEEEEEEEEC | EEEEEECCEEE | 1.60 | 1.04 | 474/6 | 34 | 0.84/0.35 | 13184 | 0.48 | EEEEEEEEEEC | 3.60 | -0.20 | 10.0/6.2/6.7 |
| 5W52 | 11 | EEEEEEEEEEC | CCEEEEEECCC | 1.12 | 0.47 | 385/3 | 7 | 0.95/0.88 | 13241 | 0.29 | CCCCCCCCCCC | 3.46 | +0.14 | 10.7/6.8/7.0 |
| 7JS6 | 15 | CEEEEECCCCCCCCC | CEEECEEECCEEECC | 3.13 | 1.92 | 1/0 | 11 | 0.83/0.29 | 5724 | 0.40 | CCCCCCEEEEECCCC | 5.40 | -0.20 | 5.2/8.4/8.7 |
| 7LCW | 15 | CCEEEEECCCEEEEE | CCCHHHHHCEEEEEE | 3.19 | 2.32 | 0/0 | 13 | 0.65/0.22 | 5298 | 0.25 | CCCCCCCCCCCCCCC | 4.81 | +0.23 | 5.2/7.2/7.6 |
| 8T63 | 16 | EEEEECEEEEEEEEEE | EEECCCCHHHHHCCCC | 3.20 | 2.43 | 0/0 | 27 | 0.95/0.74 | 3796 | 0.05 | HHHHHHHHHHHHHHHH | 5.64 | +0.02 | 12.1/6.9/7.4 |
| 9KAR | 15 | CCCCEEEEECCCCCC | EEEEEECCEEEEEEC | 1.40 | 1.40 | 7/2 | 5 | 0.81/0.64 | 980 | 0.40 | CCCCHHHHHHHHHHH | 5.34 | +0.14 | 5.1/7.6/8.1 |
| 9L1M | 12 | EEEEECCCCCCC | CHHHHHHHHHCC | 2.96 | 1.98 | 1/0 | 113 | 0.50/0.16 | 12408 | 0.05 | HHHHHHHHHHHH | 4.49 | +0.31 | 6.1/5.8/6.2 |

Every row: the near-native band IS in the pool (5-113 windows), and the score puts it in the worst
half (`band score pct mean` >= 0.50 on all 18).  `rho(score, rr)` is near zero or negative on 8 of 18.

### 3.1 It is NOT a query failure

The near-native band is exactly as (in)accessible to BLOSUM on the 18 as on the 108: identical
in-pool fraction (0.079 vs 0.073, p=0.96) and no worse median BLOSUM rank.  BLOSUM is
uninformative about nativeness on BOTH groups — the near-native band sits at median rank ~7,500
out of ~17,000 everywhere, i.e. **at chance**.  What differs is what happens to the band once it
is in the pool.

### 3.2 It IS a filter failure, and the filter is anti-ranked

On the 108, the best in-band window is at score percentile 0.009 — the score puts it essentially
first.  On the 18 it is at 0.391, and the in-band mean is at 0.687: the near-native windows are
scored in the **worst third**.  This is not a weak signal, it is a reversed one.  Seventeen of 18
have `band_score_pct_mean > 0.5` vs 11/108 (Fisher p < 1e-8).  Per target, `band_score_pct_min`
reaches 0.97 (2BFI), 0.88 (5W52), 0.74 (8T63), 0.64 (9KAR).

### 3.3 The mechanism: the top-75's shape is the distogram's shape, and it is wrong

Per-target (`forensics_causal.json:fail18_detail`), the top-75's mean radius of gyration tracks the
distogram's implied rg almost exactly, and BOTH miss the native — **in both directions**:

| target | native rg | top-75 rg | distogram rg | native SS | distogram SS |
|---|---|---|---|---|---|
| 2BFI | 11.9 | 5.7 | 5.8 | EEEEEEEEEECC | HHHHHHHHHHHH |
| 5W52 | 10.7 | 6.8 | 7.0 | EEEEEEEEEEC | CCEEEEEECCC |
| 3SGO | 10.0 | 6.2 | 6.7 | EEEEEEEEEEC | EEEEEECCEEE |
| 8T63 | 12.1 | 6.9 | 7.4 | EEEEECEEEEEEEEEE | EEECCCCHHHHHCCCC |
| 2MQ2 | 5.4 | 10.3 | 10.1 | CCCCCCCCEEEEEC | EEEEEEEEECEEEC |
| 2NB7 | 6.3 | 9.4 | 9.3 | CCCCHHHHHHHHHH | CCEEEEEECCCCCC |
| 7JS6 | 5.2 | 8.4 | 8.7 | CEEEEECCCCCCCCC | CEEECEEECCEEECC |
| 9KAR | 5.1 | 7.6 | 8.1 | CCCCEEEEECCCCCC | EEEEEECCEEEEEEC |
| 3BTB | 7.6 | 6.9 | 6.6 | CCEEEEECCEEEEEC | CHHHHHHHHHHHHHH |

The fibril segments are **extended** and get a compact answer; the lasso peptides are **compact
threaded rings** and get an extended answer.  The score is not broken — it is faithfully executing
a distogram whose implied global size and SS class are wrong for these targets.  This is a
misprediction failure that presents as a filter failure, not a compactness bias per se.
(`forensics_compact.py` quantifies the rg calibration and bounds an oracle rg-matched pool.)

### 3.4 Multimodality is present but secondary

`pool_nclusters` 160 vs 119 (p=0.011) and `top75_in_near_cluster` 8.8 vs 47.3 (p=1.2e-7): the pool
does contain a near-native basin for 12/18, and the score picks a different one.  But
`top75_in_near_cluster` is a consequence of 3.2, not an independent cause — the near basin is
exactly the one the wrong-shape prediction scores worst.

### 3.5 The common causal structure: an out-of-distribution REFERENCE

`meta_lasso_or_fibril` — readable from the PDB header before anything is computed — carries a
**10/16 failure rate vs 8/110** (p=4.4e-09) and mean answer 5.15 A vs 2.93 A.  Native strand
content is largely the same fact (E = 0.499 in lasso/fibril vs 0.168 elsewhere), but retains
independent signal within the non-lasso/fibril targets: 5/30 high-E fail vs 3/80 low-E
(Fisher p=0.034), mean answer 4.29 A (high-E) vs 2.71 A (low-E).

The unifying statement: **the 18 are peptides whose deposited conformation is not the free-solution
conformation of that sequence.**  A lasso peptide is threaded through a macrolactam ring; a fibril
segment is a single extended strand held by neighbouring molecules in the lattice; 2BP4 is 80% TFE,
3BTB is a transferred-NOE structure of a bound peptide.  The sequence does not determine this
conformation, so a sequence-conditioned distance predictor cannot predict it, and the score built
on that prediction anti-ranks the very windows that match it.

### 3.6 Deployable early warning is NULL

An LFO logistic regression over 10 deployable target-level features (distogram rg, distogram SS
composition, pool cluster count, top-75 spread, length, ESM-SS composition, alphabet log-likelihood)
predicts FAIL18 membership at **AUC 0.600**, inside the label-shuffled null (mean 0.456, p95 0.616).
Flagging the 18 highest-risk targets catches 3 of 18.  **Nothing we can compute at inference time
from this feature set knows the answer will be catastrophic.**  The only thing that does know is the
header metadata — which for a genuinely novel peptide would come from the experimentalist's intent
("this is a lasso peptide", "this is a fibril segment"), not from a structure.

## 4. Part 2a — what the near-native windows share: structural CLASS (`forensics_class.json`)

For each target, the universe's 20 nearest windows (ORACLE selection) vs 20 random windows and the
20 the shipped score ranks best.

| | near-native 20 | random 20 | score-best 20 |
|---|---|---|---|
| BLOSUM rank (median) | 7478 | 9281 | 219 |
| ESM-key rank (median) | 6868 | 9141 | 925 |
| sequence identity to target | 0.079 | 0.057 | 0.216 |
| CA-SS agreement with native | **0.654** | 0.393 | 0.490 (FAIL18: **0.261**, other: 0.528) |
| fraction from peptides | 0.324 | — | — |

Neither BLOSUM nor ESM ranks the near-native windows meaningfully above random (7478 / 6868 vs 9281 /
9141 out of ~17,000).  **Sequence identity is 7.9% — the correct template is not a homologue.**  What
the near-native windows DO share with the target is secondary structure (0.654 vs 0.393 background)
and, for the flagged classes, the class itself:

| target class | n targets | near-20 share from that class | random-20 | library base rate | enrichment | p |
|---|---|---|---|---|---|---|
| **lasso** | 6 | **0.789** | 0.079 | 0.042 | **18.8x** | 0.025 |
| **fibril** | 8 | **0.345** | 0.000 | 0.052 | **6.6x** | 0.032 |
| xray | 6 | 0.510 | 0.000 | 0.065 | 7.9x | 0.009 |
| membrane | 36 | 0.489 | 0.312 | 0.318 | 1.5x | 0.024 |
| cyclic | 10 | 0.146 | 0.075 | 0.078 | 1.9x | 0.376 |
| designed | 21 | 0.131 | 0.114 | 0.147 | 0.9x | 0.692 |

**The information needed to recognise the correct structure IS in the library, and its key is the
structural class, not the sequence.**  Concretely (rank B = BLOSUM, E = ESM-key, of ~10-25k windows):

- 2N5C (lasso) best windows: 9KAR @1.40 A (B826, **E64**), 7BZ9 @1.62 (B1018, **E170**), 6M19 @1.68 (B681, **E31**)
- 9KAR (lasso) best windows: 2N5C @1.40 (B487, **E52**), 7BW5 @1.59 (B92, **E115**), 7BZ9 @1.76 (B812, **E93**)
- 7JS6 (lasso) best: 7BZ9 @1.92 (B4832, E2133), 6XTI @2.16 (B4831, E1711)
- 5W52 (fibril) best: **2BFI @0.47 A** (B18149, E17391) — the other amyloid zipper in the database
- 3SGO (fibril) best: protein fragments only (B21377, E13495)
- 9L1M (fibril) best: 1OI0_26 @1.98 (B2142, E12443) — BLOSUM better than ESM here

So ESM-2 embedding similarity **does** retrieve the lasso class (ranks 31-170, inside K=500 where
BLOSUM puts them at 487-1018 and misses them), and **does not** retrieve the fibril class (both keys
are at rank 10,000+).  The lasso peptides are recognisably a sequence family; the amyloid segments
are not — they are ordinary-looking sequences whose conformation is imposed by the lattice.

## 5. Part 2c — the structural alphabet is not predictable from sequence (`forensics_alphabet.json`)

Codebook: k-means K=8 on (cos/sin phi, cos/sin psi) over **fragment residues only**.  Centres
(phi, psi in degrees): (-143,160) beta, (-116,75), (-114,129) PPII/beta, (-96,-4) bridge, (-73,144)
PPII, (-64,-39) alphaR, (74,18) alphaL, (82,-160).  Predictor: per-residue MLP (256 hidden) over
ESM-2 pca128 of residues i-2..i+2 plus position features, trained **leave-fold-out** on the same
library the retrieval universe uses (~92k residues, ~6.6k members per fold).

| per-residue accuracy | all 126 | FAIL18 | other-108 |
|---|---|---|---|
| predicted 8-state alphabet | 0.375 | 0.224 | 0.400 |
| **majority-class baseline (always alphaR)** | **0.393** | 0.242 | 0.419 |
| predicted 3-state SS from ESM | 0.554 | 0.417 | 0.577 |
| **majority-class SS baseline** | **0.575** | 0.400 | 0.604 |
| shipped distogram's implied SS | 0.447 | 0.383 | 0.457 |
| held-out training accuracy of the alphabet MLP (library residues) | 0.645 | — | — |

**The LFO structural-alphabet predictor does not beat the majority class on the targets** (0.375 vs
0.393), even though it reaches 0.645 on held-out library residues.  The same is true of the 3-state
SS predictor (0.554 vs 0.575).  This is the distribution shift of S7-2 appearing again in a second
model family: fragments-in-proteins are predictable from sequence; free 9-16-mers are not.  A
Foldseek-style structural-alphabet key is therefore built on a predicted string that carries
essentially no information about these targets, which caps the whole idea before retrieval starts.

## 5b. Part 2b — retrieval keys at matched K=500, full chain (`forensics_part2.json`)

Every key ranks the WHOLE universe, takes the top 500, and runs the shipped chain (distogram
top-75 -> coordinate average -> projection).  Outcome = post-projection `rmsd_arm`.  Control =
`blosum` (which reproduces the production record exactly).  n=126 paired.

| key | pool_best | pool_mean | in-band | rmsd_arm all/F18/O108 |
|---|---|---|---|---|
| **blosum** (shipped) | 1.711 | 4.453 | 134.1 | **3.213 / 6.019 / 2.745** |
| esm (pca128 mean cosine) | 1.685 | 4.354 | 142.2 | 3.215 / 5.946 / 2.760 |
| ss_pred (distogram-implied SS agreement) | 2.037 | 4.449 | 146.0 | 3.369 / 5.838 / 2.957 |
| hyb1 / hyb3 / hyb10 (BLOSUM + lam x SS agree) | 1.689 / 1.716 / 1.897 | 4.29 / 4.25 / 4.33 | 153 / 162 / 157 | 3.226 / 3.255 / 3.302 |
| ss_esm (LFO ESM SS predictor) | 1.899 | 4.296 | 145.8 | 3.305 / 6.059 / 2.846 |
| sshyb3 | 1.708 | 4.177 | 165.8 | 3.237 / 5.998 / 2.777 |
| alpha_pred (predicted 8-state alphabet) | 2.019 | 4.018 | 205.3 | 3.395 / 5.969 / 2.966 |
| alpha_ll (soft alphabet log-likelihood) | 2.067 | 3.870 | 229.1 | 3.307 / 5.662 / 2.915 |
| ahyb3 / allhyb3 | 1.797 / 1.913 | 4.02 / 3.87 | 196 / 224 | 3.225 / 3.237 |
| score_univ (distogram score over the whole universe) | 2.161 | **3.605** | 235.4 | 3.362 / 5.965 / 2.929 |
| *oss* ORACLE native-SS key | 1.587 | 4.247 | 150.1 | 3.033 / 5.429 / 2.634 |
| *ohyb3* ORACLE BLOSUM+3xSS | 1.558 | 4.160 | 169.6 | 3.090 / 5.533 / 2.683 |
| *oalpha* ORACLE native-alphabet key | 1.532 | 3.720 | 242.3 | 3.023 / 5.477 / 2.614 |
| *oahyb3* ORACLE BLOSUM+3xalphabet | 1.564 | 3.858 | 212.9 | 3.090 / 5.581 / 2.674 |

Paired vs `blosum` (negative = better):

| key | d rmsd_arm [CI95] | W/L | **drop-top-10** | d FAIL18 [CI95] | d other-108 |
|---|---|---|---|---|---|
| esm | +0.003 [-0.054,+0.057] | 61/65 | +0.066 | -0.074 [-0.194,+0.040] | +0.015 |
| hyb1 | +0.013 [-0.030,+0.056] | 63/63 | +0.057 | -0.019 | +0.018 |
| ahyb3 | +0.012 [-0.051,+0.071] | 57/69 | +0.085 | +0.025 | +0.010 |
| sshyb3 | +0.024 [-0.029,+0.076] | 58/68 | +0.079 | -0.022 | +0.032 |
| allhyb3 | +0.025 [-0.059,+0.104] | 61/65 | +0.106 | -0.096 | +0.045 |
| hyb3 | +0.042 [-0.025,+0.110] | 50/76 | +0.108 | -0.081 | +0.063 |
| hyb10 | +0.090 [-0.005,+0.193] | 50/76 | +0.179 | -0.138 | +0.128 |
| ss_esm | +0.092 [+0.017,+0.169] | 55/71 | +0.157 | +0.039 | +0.101 |
| alpha_ll | +0.095 [-0.047,+0.228] | 59/67 | +0.219 | **-0.358 [-0.882,-0.044]** | +0.170 |
| score_univ | +0.150 [+0.047,+0.264] | 53/73 | +0.233 | -0.055 | +0.184 |
| ss_pred | +0.156 [+0.039,+0.283] | 52/74 | +0.265 | -0.181 | +0.212 |
| alpha_pred | +0.183 [+0.069,+0.303] | 50/76 | +0.270 | -0.051 | +0.222 |
| *oss* ORACLE | **-0.179** [-0.327,-0.034] | 80/46 | **+0.004** | -0.590 [-1.348,+0.134] | -0.111 |
| *oalpha* ORACLE | **-0.190** [-0.314,-0.074] | 80/46 | **-0.044** | -0.543 [-1.186,+0.026] | -0.131 |
| *ohyb3* ORACLE | -0.123 [-0.236,-0.024] | 78/48 | +0.001 | -0.487 | -0.062 |
| *oahyb3* ORACLE | -0.123 [-0.221,-0.042] | 64/62 | -0.009 | -0.439 [-1.012,-0.031] | -0.070 |

**BLOSUM/ESM fusion** (`forensics_part2_fuse.json`, run separately, same protocol):

| key | pool_best d [CI] | rmsd_arm d [CI] | W/L | drop-10 | d FAIL18 | d other-108 |
|---|---|---|---|---|---|---|
| union (250 BLOSUM + 250 ESM) | **-0.039 [-0.076,-0.004]** | -0.020 [-0.062,+0.018] | 61/65 | +0.029 | -0.042 | -0.016 |
| rrf (reciprocal rank fusion, k=60) | **-0.047 [-0.083,-0.013]** | -0.004 [-0.044,+0.035] | 64/62 | +0.040 | -0.026 | +0.000 |

Both fusions **do** improve the pool (pool_best CI excludes zero, FAIL18 -0.107 / -0.125) and both
deliver **nothing** downstream (answer CI spans zero, drop-10 positive).  This is the cleanest single
demonstration in the report that the filter consumes whatever retrieval gains: a measurably better
pool produces an unmeasurably different answer.

**Not one deployable key beats BLOSUM.**  ESM is an exact tie (+0.003, CI spans zero, 61W/65L);
everything else is worse.  Three further points:

1. **Even the ORACLE SS key fails the concentration rule.**  `oss` is -0.179 [-0.327,-0.034], but
   dropping the 10 best-helped targets leaves **+0.004** — the entire benefit of a perfect
   secondary-structure retrieval key lives in 10 of 126 targets.  Only the 8-state oracle alphabet
   key survives at all (-0.044 after drop-10), and only because it also encodes local geometry
   beyond H/E/C.  Combined with section 7's `o_shape_only` (-0.348, drop-10 -0.149, which adds the
   SIZE match), the ordering is: SS alone < SS + size, and both are oracles.
2. **`alpha_ll` on FAIL18 (-0.358 [-0.882,-0.044]) is very likely multiplicity.**  Sixteen keys were
   tested against the FAIL18 subgroup; one 95% CI excluding zero is what chance produces.  It is
   worse overall (+0.095) and on the 108 (+0.170).  Do not act on it without a pre-registered replication.
3. **Correction to the record on `score_univ`.**  The record says a universe-wide score key gives a
   worse pool mean.  Measured here it gives a *much better* pool mean (3.605 vs 4.453) and a much
   worse pool best (2.161 vs 1.711), and on the FAIL18 the in-band count collapses from 43.4 to
   **0.167** — it retrieves windows the (wrong) distogram likes, so it amplifies exactly the error
   that causes the failure.  The answer is +0.150 worse.  The direction of the conclusion stands;
   the stated reason does not.

**Error correlation.**  Pearson r between each key's per-target answer-RMSD vector and BLOSUM's:
esm 0.983, hyb1 0.990, sshyb3 0.984, ahyb3 0.979, hyb3 0.975, ss_esm 0.969, allhyb3 0.964, hyb10
0.945, score_univ 0.934, alpha_pred 0.929, ss_pred 0.915, alpha_ll 0.899, and even the oracles
oss 0.872 / oalpha 0.916 / ohyb3 0.935 / oahyb3 0.954.  **Every key fails on the same targets.**
An oracle that picked the best of all 19 keys per target would reach 2.636 (FAIL18 4.893) versus
BLOSUM's 3.213 (6.019) — but with r >= 0.87 across the board there is no decorrelation to exploit,
which is the "decorrelated errors exist and cannot be exploited" result reappearing on the
retrieval axis: the shared component is the target, not the key.

## 6. The SIZE channel: the strongest oracle signal, and why it is not deployable

### 6.1 The distogram's implied size barely tracks the true size (`forensics_compact.json`)

| | slope vs native rg | pearson | MAE (A) | MAE FAIL18 | MAE other-108 |
|---|---|---|---|---|---|
| shipped distogram's implied rg | 0.238 | +0.363 | 0.979 | **2.440** | 0.736 |
| LFO ridge on ESM pca128 (target sequence only) | 0.213 | +0.341 | 1.080 | 2.049 | 0.918 |
| LFO GBT on ESM pca128 | 0.151 | +0.345 | 0.958 | **1.647** | 0.843 |
| constant (training mean) | 0.001 | +0.042 | 1.681 | 2.165 | 1.600 |

A dedicated leave-fold-out size regressor trained on 787-minus-fold peptides (excluding every
tuning-target sequence) is **no better than the shipped distogram**, and both have slope ~0.2: they
predict the library mean plus a small correction.  Size prediction fails by 1.6-2.4 A of rg on
exactly the 18 targets where size is the thing that is wrong.

### 6.2 Size as a retrieval GATE: oracle real, deployable null (`forensics_rggate.json`)

Gate the WHOLE universe by |rg_window - rg_hat| < w, then BLOSUM top-500, then the shipped chain.
Outcome = `rmsd_avg` (the coordinate average before projection; see the proxy note in section 9).
All paired against the production BLOSUM pool, n=126.

| arm | rg_hat | w | mean | FAIL18 | other-108 | d (all) [CI] | d FAIL18 [CI] | d other | drop-top-10 |
|---|---|---|---|---|---|---|---|---|---|
| base (production) | — | — | 3.048 | 5.832 | 2.584 | — | — | — | — |
| **o_gate10** ORACLE | native rg | 1.0 | 2.738 | **4.699** | 2.411 | **-0.310** [-0.455,-0.190] | **-1.133** [-1.859,-0.516] | -0.173 | -0.121 |
| o_gate05 ORACLE | native rg | 0.5 | 2.745 | 5.093 | 2.354 | -0.303 [-0.446,-0.183] | -0.739 [-1.477,-0.141] | -0.231 | -0.129 |
| gbt_gate10 DEPLOYABLE | LFO GBT | 1.0 | 3.026 | 5.256 | 2.655 | -0.022 [-0.121,+0.069] | -0.576 [-1.013,-0.222] | +0.070 | +0.086 |
| gbt_gate15 DEPLOYABLE | LFO GBT | 1.5 | 3.033 | 5.431 | 2.633 | -0.015 [-0.093,+0.053] | -0.401 [-0.790,-0.108] | +0.049 | +0.068 |
| disto_gate10 DEPLOYABLE | distogram rg | 1.0 | 3.111 | 5.745 | 2.672 | +0.062 [+0.018,+0.110] | -0.087 [-0.223,+0.042] | +0.087 | +0.098 |
| **shuf_gate10 NULL** | another target's rg | 1.0 | 3.365 | 5.198 | 3.059 | +0.316 [+0.124,+0.522] | **-0.634** [-1.216,-0.240] | +0.475 | +0.482 |
| const_gate10 CONTROL | library mean rg | 1.0 | 3.392 | 5.355 | 3.065 | +0.344 [+0.192,+0.498] | -0.477 [-0.846,-0.143] | +0.481 | +0.478 |

**The null control is the whole story.** A gate centred on *another target's* true rg — carrying zero
target-specific information — already buys **-0.634 A on the FAIL18**.  Subtracting the null:

| | oracle rg | GBT-predicted rg | interpretation |
|---|---|---|---|
| FAIL18, vs base | -1.133 | -0.576 | both "help" |
| FAIL18, vs the shuffled null | **-0.499** | **+0.058** | the deployable predictor has **zero** skill here |
| other-108, vs the shuffled null | -0.648 | -0.405 | the predictor has real skill where it wasn't needed |

Most of the apparent FAIL18 benefit of any size gate is **generic pool restriction**, i.e. outlier
avoidance — the same mechanism the record already priced as "consensus is outlier avoidance, not a
nativeness signal".  The genuinely informational part of the size channel is worth ~0.50 A on the 18
and ~0.65 A on the 108, and the deployable predictor captures 0% of it on the 18 and ~60% of it on
the 108.  This is `MAE does not price selected RMSD` in its sharpest form: the size predictor's mean
error is fine, and it is useless precisely where the mean is set.

## 7. Structure as a retrieval KEY: the S6-6 open hypothesis, answered (`forensics_classkey.json`)

S6-6 measured an ORACLE secondary-structure FILTER before the argmin at 0.21 A.  Whether SS as a
retrieval KEY (changing which windows enter the pool) is worth more was open.  Arms at matched
K=500; outcome `rmsd_avg`; paired vs the production BLOSUM pool.

| arm | pool_best | in-band | mean | FAIL18 | other-108 | d rmsd_avg [CI] | d FAIL18 | d other | drop-10 | W/L |
|---|---|---|---|---|---|---|---|---|---|---|
| blosum (production) | 1.711 | 134.1 | 3.048 | 5.832 | 2.584 | — | — | — | — | — |
| class_meta (BLOSUM + PDB-class bonus) | 1.672 | 138.3 | 3.045 | 5.841 | 2.579 | -0.004 [-0.013,+0.006] | +0.009 | -0.006 | +0.006 | 47/41 |
| class_only (PDB class as primary key) | 1.635 | 153.6 | 3.021 | 5.696 | 2.575 | -0.027 [-0.105,+0.050] | -0.136 | -0.009 | +0.062 | 44/49 |
| **o_shape** ORACLE (BLOSUM + 3x[native-SS agreement + size match]) | 1.504 | 200.8 | 2.875 | 5.373 | 2.458 | -0.174 [-0.284,-0.070] | -0.458 | -0.126 | -0.045 | 85/41 |
| **o_shape_only** ORACLE (shape agreement as the PRIMARY key) | 1.497 | 201.7 | 2.700 | **4.716** | 2.364 | **-0.348** [-0.509,-0.198] | **-1.116** | -0.220 | -0.149 | 91/35 |

Two conclusions:

1. **SS-as-a-key beats SS-as-a-filter, but only from 0.21 A to 0.35 A.**  Making the pool's shape
   right raises the in-band count by 68 windows (134 -> 202) and lowers pool_best by 0.21 A, and the
   answer still only moves 0.35 A.  Even with *perfect* structural knowledge steering retrieval, the
   FAIL18 stay at 4.72 A and the whole set at 2.70 A.  **Retrieval is not the bottleneck.**
2. **The structural class, though 18.8x enriched in the near-native windows, is worth nothing as a
   key** (-0.027 A, CI spans zero, drop-10 +0.062).  The class flags are too coarse (membrane alone
   covers 32% of the library) and the score still anti-ranks whatever the class brings in.

## 7b. Where the 18 would land if each stage were perfect (`forensics_ceiling.json`)

ORACLE/DIAGNOSTIC bounds, one component replaced at a time, outcome `rmsd_avg`:

| arm | all | **FAIL18** | other-108 | frac < 2.0 A |
|---|---|---|---|---|
| base (production BLOSUM-500 -> shipped score top-75 -> average) | 3.048 | **5.832** | 2.584 | 0.29 |
| **perfect FILTER, m=75** (the 75 lowest-rr members of the SAME pool) | 1.963 | **3.069** | 1.779 | 0.56 |
| **perfect FILTER, m=25** | **1.644** | **2.416** | 1.515 | 0.71 |
| perfect KEY + perfect FILTER (25 lowest-rr of the whole universe) | 1.106 | 1.470 | 1.045 | 0.95 |
| pool argmin (single best member, for reference) | 1.711 | 2.284 | 1.615 | 0.58 |
| universe best window | 1.313 | 1.640 | 1.259 | 0.81 |

(Instrument cross-check: `ofilt25` = **1.644** reproduces correction C1's "oracle top-25 averaging
1.644" exactly, independently derived here.)

**The filter is worth 1.40 A; the key, on top of a perfect filter, a further 0.54 A.**  A perfect
filter alone takes the FAIL18 from 5.83 to 2.42 and the whole set below 2.0 A — the sprint's target
— *without changing retrieval at all*.  This is the quantitative answer to Part 2's question: the
information needed to recognise the correct structure is **already inside the shipped K=500 pool**
for both groups; retrieval is not what is missing.  It also explains why the oracle shape key
(section 7) only bought 0.35 A: it improves a pool whose selection step then discards the gain.

## 8. Early warning, posed correctly (`forensics_warning.json`)

Section 3.6 asked "can we classify the FAIL18?" and got a null (AUC 0.600 vs null p95 0.616).  The
better-posed question is leave-fold-out REGRESSION of the production answer's CA-RMSD on 17
deployable features.  Single-feature rank correlations with the answer (permutation p over 2000
shuffles):

| deployable feature | rho vs answer RMSD | p_perm | rho on other-108 |
|---|---|---|---|
| top-75 cluster count (2.5 A avg linkage) | **+0.526** | <5e-4 | — |
| top-75 mean pairwise RMSD | +0.516 | <5e-4 | — |
| top-75 rg spread | +0.501 | <5e-4 | — |
| distogram helix fraction | -0.449 | <5e-4 | — |
| distogram strand fraction | +0.427 | <5e-4 | — |
| ESM-SS helix fraction | -0.415 | <5e-4 | — |
| ESM-SS strand fraction | +0.377 | <5e-4 | +0.376 |
| pool cluster count | +0.314 | 5e-4 | +0.238 |
| distogram implied rg | +0.249 | 0.006 | +0.168 |

LFO ridge over all 17: **Spearman +0.516** with the answer (null mean -0.002, p95 +0.192).  So a
genuine confidence signal exists.  But it is a *spread* signal, not a *nativeness* signal, and it
does not find the catastrophes:

| abstain on the k worst-flagged targets | k=0 | k=6 | k=12 | k=18 | k=25 |
|---|---|---|---|---|---|
| LFO ridge | 3.215 | 3.169 | 3.076 | 3.000 | 2.828 |
| ORACLE abstention | 3.215 | 3.015 | 2.847 | 2.703 | 2.550 |

Flagging 18 targets catches **6 of the 18** and moves the mean by 0.215 A where perfect abstention
would move it 0.512 A.  **You can rank targets by expected error; you cannot name the catastrophes.**
Every top feature measures how much the top-75 disagrees with itself — the ensemble knows it is
uncertain, but not that it is wrong in a particular direction.

## 9. Method notes

- **`rmsd_avg` as the outcome for the secondary arms.**  Sections 6-7 report `rmsd_avg` (the
  coordinate average before the L-BFGS projection) instead of `rmsd_arm`, because the projection
  costs ~5 s per arm per target and those experiments have 6-8 arms.  Fidelity, measured on 1,207
  arm-target pairs from Part 2 where both were computed: levels r=0.995, mean(arm - avg)=+0.171 A
  (a near-constant offset); **paired differences r=0.971, slope 1.084, mean |disagreement| 0.090 A**.
  Every Part 2 table (section 5) reports the true post-projection `rmsd_arm`.
- **SS assignment.**  `instrument.ss_of` rebuilds an ideal backbone from torsions and runs a
  simplified DSSP; it agrees with a CA-geometry assigner on only 62% of native residues and calls
  the 1JBF hairpin and the 5W52 steric zipper "coil".  All SS comparisons here use
  `forensics_lib.ca_ss`, a P-SEA-style CA-only assigner, applied identically to natives and windows.
- **Tie-breaking.**  Every key is ranked with `np.argsort(kind="stable")` after adding 1e-3 uniform
  jitter to integer-valued keys (SS agreement, alphabet agreement), so a tied key cannot inherit the
  universe's construction order the way `np.argmin` on a tied signal did in the S11 trap.

## 10. Leakage audit

- No file under `s9/final_cache/`, `results/benchmark_manifest.json`, `bench_results/` (except the
  read-only production record for tuning targets), or any dev24 target was read.  All 126 targets
  used are the `s8/generate_univ` tuning set.
- `rr` and `nat_ca` appear only in (a) evaluation fields, (b) arms whose names begin `o_`/`ORACLE`/
  `DIAGNOSTIC` and which are labelled as such in every table.  No deployable arm reads them.
- The alphabet/SS/rg predictors are trained leave-fold-out on library members only
  (`peptide_db.folds(5)` + `core.data.fold_fragments(fold)`), and the rg predictor additionally
  excludes every tuning-target sequence by exact match.  The alphabet codebook is fitted on
  **fragment** residues only.
- The library universes are already out-of-fold by construction (`s8/generate.stage_univ`), and the
  parent map was verified to reproduce every window's codes exactly, so no in-fold parent leaked in.
- The PDB header text (`header_class`) reads HEADER/TITLE/KEYWDS/EXPDTA/COMPND/REMARK 210/SSBOND/LINK
  records only — metadata, never coordinates.  It is nevertheless **not** a deployable inference-time
  feature for a novel peptide and every table says so.
- The FAIL18/other-108 split is the pinned `instrument.FAIL18`; no arm was selected on it.
- Gate widths (0.5/1.0/1.5) and the hybrid lambdas (1/3/10) are reported for every value tried;
  nothing was chosen after seeing the outcome.

## 11. Answers to the two questions

**Q1. Why do the 18 fail, and what is the common causal structure?**

Not a query failure (the near-native band is equally accessible to BLOSUM on both groups: in-pool
fraction 0.079 vs 0.073, p=0.96).  It is a **filter failure driven by a reference that the sequence
does not determine**.  The 18 are, at a rate of 10/16 vs 8/110 (p=4.4e-09), lasso peptides and
amyloid/fibril segments — conformations imposed by a covalent thread or by a crystal lattice, not by
the free-solution folding of that sequence.  The distogram therefore predicts the wrong global shape
(rg calibration slope 0.238; SS agreement 0.383), the Bayes-risk score faithfully executes that wrong
shape (top-75 rg tracks the distogram rg to within ~0.5 A in every case, missing the native by 2-5 A
in **both** directions), and the near-native windows land in the worst third of the score
(percentile 0.687 vs 0.287, p=3.4e-10; best in-band window at 0.391 vs 0.009).  Multimodality is a
consequence, not a separate cause.

**Q2. Is the information in the library, and under what key?**

The information is present *twice over*, and the key is not the problem:

- **In the pool already.**  Perfect filtering of the shipped K=500 pool at m=25 gives 1.644 A overall
  and 2.416 A on the FAIL18 — the sprint's target — with retrieval untouched (section 7b).
- **In the library, keyed by structural class.**  The near-native windows are 18.8x enriched in
  lasso parents for lasso targets and 6.6x in fibril parents for fibril targets; 5W52's best window
  (0.47 A) is a window of 2BFI, the other amyloid zipper in the database (section 4).
- **But no key exploits it.**  Sequence identity of the near-native windows is 7.9% — they are not
  homologues.  None of 12 deployable keys beats BLOSUM (best: ESM at +0.003, an exact tie).  The
  class as an explicit key buys -0.027 A (CI spans zero).  Even ORACLE keys are limited: perfect SS
  as a key is -0.179 but **+0.004 after drop-top-10**; perfect SS+size as a key is -0.348
  (drop-10 -0.149) and still leaves the FAIL18 at 4.72 A.  All 19 keys' error vectors correlate
  r >= 0.87 with BLOSUM's.

The structural-alphabet (Foldseek-in-miniature) idea is closed at its first step: an LFO per-residue
alphabet predictor reaches 0.645 on held-out library residues but **0.375 on the targets, below the
0.393 majority-class baseline**.  The string you would retrieve with carries no information.

## 12. Ranked recommendations for the coordinator

1. **Spend the sprint on the FILTER, not on retrieval.**  Perfect filtering of the existing pool is
   worth 1.40 A (3.048 -> 1.644 on `rmsd_avg`; FAIL18 5.83 -> 2.42); perfect retrieval on top of it
   only a further 0.54 A.  Every retrieval experiment in this report — 12 deployable keys, 4 oracle
   keys, 2 metadata keys, 6 size gates — returns null or negative, and even the oracles fail the
   drop-top-10 rule.  Retrieval is measured out.
2. **The single experiment I would run next: learned aggregation over the full candidate-vs-objective
   deviation map (correction C2's untested class), trained and evaluated with the FAIL18 held out of
   training.**  Rationale: it is the only in-band discrimination class not yet tried; the ceiling it
   is aiming at (1.644) is the largest unclaimed quantity in the system; and the FAIL18 forensics say
   exactly what the input must be — the per-pair signed deviation between a candidate's distances and
   the objective's, not a scalar pooling of it, because the failure is a *coherent directional*
   error (whole-structure too compact or too extended) that any scalar pooling averages away.  Hold
   the 18 out of training so the result cannot be read as fitting the failure mode.
3. **Do not build a structure-aware retrieval key.**  Sections 4, 5b, 7 and the alphabet result close
   this.  If any structural key is revisited, it must be SS **plus size** and it must be shown to
   survive drop-top-10 (SS alone does not: +0.004).
4. **Do not pursue a size/rg predictor.**  Oracle size is worth -1.13 A on FAIL18, but a shuffled-
   label null recovers -0.634 of that, and the deployable LFO predictor recovers **0%** of the
   remainder on the 18 (it does recover ~60% on the 108, where it is not needed).
5. **Consider scoping the claim rather than fixing the 18.**  10 of the 16 lasso/fibril targets fail;
   these are conformations a sequence-conditioned model has no mechanism to predict.  A defensible
   headline is "free-solution monomeric peptides", with the covalently-threaded and lattice-imposed
   entries reported separately.  This is a reporting decision, not a scientific fix, and needs the
   coordinator's judgement — I am not proposing to drop them from the instrument.
6. **A confidence signal is available if wanted** (LFO ridge, Spearman +0.516 with the answer, null
   p95 +0.192), driven by top-75 internal diversity.  It ranks error but does not find the
   catastrophes (6 of 18 in the flagged 18), so it supports abstention/uncertainty reporting, not
   triage.

**Nothing here has been run on dev24 or benchmark60.**  If the coordinator wants a dev pass on any
arm, the only one with a positive deployable point estimate is ESM-as-key (+0.003 overall, -0.074 on
FAIL18) — and on this instrument it is a tie, so I do not recommend spending dev24 on it.
