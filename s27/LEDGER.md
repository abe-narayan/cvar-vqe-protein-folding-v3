# S27 LEDGER -- ALTERNATIVE HAMILTONIANS (2026-09-14)

## L1 -- SCOPE AND GROUND
Question: can a different Hamiltonian or a principled combination beat the shipped distogram
Bayes-risk score on the 126-target instrument with the production readout and the genuine
CVaR-VQE intact. Pre-registration `s27/PREREG.md` (base + addenda 1, 2) before every wave.
Anchors reproduced: point cloud 3.048338 (DIS top-75, 126/126), built chain 3.2126 (rebuild
basis), the S25 VQE arm 3.0580. Nothing in production, the benchmark or an earlier sprint touched.

## L2 -- WAVE 1 (203 configurations, point cloud): NOTHING BEATS DIS; 45 WORSE BEYOND MDE
`s27/results/pool_summary.json`. Best effect -0.026 A (DIS+0.5*SS_MATCH, 0.35x MDE). Every new
single worse than DIS; most worse than a random 75-subset. DISTPOT (Sippl-style pair potential
fitted on the leakage-safe universe) is the best new single at 3.199 (-0.222 vs random, 0.9x MDE).

## L3 -- MECHANISM: RANKING INFORMATION IS ANTI-USEFUL ON THE AVERAGING READOUT
`s27/results/redundancy.json`, `mechanism.json`. Spearman(partial rho of X with the ORACLE RMSD
given DIS, effect of DIS+X on the top-75 average) = +0.48 over 18 channels; CONS has the most
information (+0.25) and does the most harm (+0.29 A, worse than its own permuted control at
1.7x MDE). The set-mean law (S12) and the outlier-removal cost (S23 L5) in one measurement.

## L4 -- WAVE 2 (m-ladder, in-band, non-additive, near-corpus): NULL; DIS+CONS FALSIFIED AT EVERY m
`s27/results/wave2_summary.json`. No arm beats DIS at any m in 3..150; the pre-registered
small-m prediction for DIS+CONS fired the other way (+0.159 at m = 3). In-band re-ranking null;
CONS worst in-band (+0.184). Per-arm minima over m transfer 19 to 29% in split halves.

## L5 -- WAVE 3: THE SS_MATCH BUMP IS A GRID ARTEFACT
`s27/results/wave3_summary.json`. With the four-weight grid the nested held-out effect is
-0.008 A (0.1x MDE), `best_of_k_within` "NOT A SIGNAL (split-half 10%)"; SS_MATCH2 (H-bond SS
call, universe propensities) no better.

## L6 -- THE GENUINE CVaR-VQE ARM (80 configurations x 126 x 2 seeds) CHANGES NO VERDICT
`s27/results/arms_summary.json :: vqe`. Set-equality gate 126/126 on every cell; residual
RMSD_VQE - RMSD_top-m <= 6.3e-4 A for the 47 tie-free configurations (tie ambiguity for the
33 tied ones; 99%-tied channels give the uniform state). Best vs DIS-VQE: -0.037 A (0.47x),
replicates on seed 1 (-0.029), not a result; 51 configurations worse with the fold CI above zero.

## L7 -- BUILT CHAIN (10 finalists): NOTHING BEATS DIS
`arms_summary.json :: chain`. DIS+0.5*SS_MATCH +0.0003; REJ[POOLGO]->DIS -0.017 (0.34x);
DISTPOT +0.171, DIS+CONS +0.249, CONS +0.553 worse with the fold CI above zero.

## L8 -- TRAINABILITY (H8): HAMILTONIAN-INDEPENDENT
`trainability.json`. Gradient variance at n = 9, alpha 0.18, T 0.5 within 0.86x to 1.14x of
DIS's under the bounded standardisation, within 0.87x to 1.09x under zrank. Falsifier (2x) silent.

## L9 -- REGIMES
`strata.json`. Every alternative helps on FAIL18 by 0.1 to 0.4 A and is null-to-harmful on the
108 others; the switch is not native-free (routers: S22 L7, S23 L7, S26 L115; H6 here null).

## L10 -- CLOSE
`s27/REPORT.md` (tables from `make_tables.py`). Verdict: the shipped score stands; no candidate
for the next architecture; the mechanism (L3) is the sprint's finding; follow-ups in §10.
Caches `s27/cache/*.npz` (16 MB, regenerable in 13 min by `run_pool.py`) are not committed.
