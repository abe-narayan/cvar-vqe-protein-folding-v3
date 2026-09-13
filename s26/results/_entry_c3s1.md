`s26/ph_c3.py stage1`, `s26/results/ph_c3_stage1.json` (complete 126/126), job
`s26/jobs_done/ph_c3_stage1.json` (exit 0, 10 s, peak RSS 0.041 GB, held 345 s at the job cap).
Pre-registered in `s26/PREREG_c3_control.md` sections 2 and 3 and addendum 1 (prediction:
AMBER worse than the matched random move by about +0.013 A, cosine negative). No AMBER compute:
the production coordinates `ca` and `amber_ca` from `bench_results/cache/1fc9f2dcf489e2fb`.

Bases, named on both sides. Arm: input the BUILT CHAIN (`rmsd_arm`, 3.2148), output the
RELAXED CHAIN (`rmsd_full`, 3.2355), both reproduced from the stored coordinates to 1e-6.
Controls: the built chain's CA trace displaced by a vector of the SAME per-atom RMS magnitude
as AMBER's displacement (0.220 A mean, measured after superposing `amber_ca` onto `ca`);
`rand` = S16's construction (isotropic Gaussian, six rigid components projected out on
`rigid_basis`, 16 draws, mean); `member` = the same magnitude along the straight line toward a
random member of the shipped K=500 pool superposed onto `ca` (16 draws, mean; 6 of 2,016 draws
overshot the member). A displaced CA trace is not an ideal-geometry chain: the controls match the
operator's SPACE (magnitude), not its validity, exactly as S16 did. Full-chain CA-RMSD to the
native, ORACLE evaluation of native-free operators. `ST.fmt`, verbatim, all seven contrasts:

    stage1 AMBER (relaxed) minus do-nothing (built chain)
      a 3.2355 (med 2.9757)   b 3.2148 (med 2.9661)   n=126
      effect +0.0207   median +0.0158   SE 0.0034   MDE 0.0096   effect/MDE +2.16
      iid  CI95 [+0.0141, +0.0274]
      fold CI95 [+0.0154, +0.0290]   folds same sign 5/5   per-fold 0:+0.021 1:+0.037 2:+0.017 3:+0.013 4:+0.017
      40W/86L/0T   worst degradation +0.1653 (2MID)   p90 +0.0640   power 1.00  Type-M 1.00
      concentration: drop-top10 +0.0258 vs uniform-effect null p10/p50/p90 +0.0213/+0.0256/+0.0302 -> pctile 0.529
      VERDICT: WORSE
    stage1 AMBER minus matched-magnitude RANDOM (16 draws)
      a 3.2355 (med 2.9757)   b 3.2244 (med 2.9669)   n=126
      effect +0.0111   median +0.0070   SE 0.0035   MDE 0.0099   effect/MDE +1.12
      iid  CI95 [+0.0044, +0.0178]
      fold CI95 [+0.0062, +0.0171]   folds same sign 5/5   per-fold 0:+0.011 1:+0.023 2:+0.009 3:+0.003 4:+0.010
      52W/74L/0T   worst degradation +0.1518 (2MID)   p90 +0.0539   power 0.88  Type-M 1.07
      concentration: drop-top10 +0.0168 vs uniform-effect null p10/p50/p90 +0.0122/+0.0168/+0.0216 -> pctile 0.497
      VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.07x]
    stage1 AMBER minus matched-magnitude TOWARD-MEMBER (16 draws)
      a 3.2355 (med 2.9757)   b 3.1969 (med 2.9755)   n=126
      effect +0.0385   median +0.0338   SE 0.0061   MDE 0.0170   effect/MDE +2.27
      iid  CI95 [+0.0272, +0.0506]
      fold CI95 [+0.0298, +0.0479]   folds same sign 5/5   per-fold 0:+0.039 1:+0.052 2:+0.049 3:+0.029 4:+0.027
      29W/97L/0T   worst degradation +0.3497 (2BP4)   p90 +0.1092   power 1.00  Type-M 1.00
      concentration: drop-top10 +0.0485 vs uniform-effect null p10/p50/p90 +0.0407/+0.0481/+0.0561 -> pctile 0.521
      VERDICT: WORSE
    stage1 RANDOM minus do-nothing
      a 3.2244 (med 2.9669)   b 3.2148 (med 2.9661)   n=126
      effect +0.0096   median +0.0071   SE 0.0015   MDE 0.0042   effect/MDE +2.28
      iid  CI95 [+0.0070, +0.0127]
      fold CI95 [+0.0077, +0.0122]   folds same sign 5/5   per-fold 0:+0.010 1:+0.015 2:+0.008 3:+0.010 4:+0.007
      29W/97L/0T   worst degradation +0.1147 (1I93)   p90 +0.0276   power 1.00  Type-M 1.00
      concentration: drop-top10 +0.0114 vs uniform-effect null p10/p50/p90 +0.0094/+0.0113/+0.0135 -> pctile 0.528
      VERDICT: WORSE
    stage1 TOWARD-MEMBER minus do-nothing
      a 3.1969 (med 2.9755)   b 3.2148 (med 2.9661)   n=126
      effect -0.0178   median -0.0118   SE 0.0043   MDE 0.0120   effect/MDE -1.49
      iid  CI95 [-0.0266, -0.0096]
      fold CI95 [-0.0255, -0.0124]   folds same sign 5/5   per-fold 0:-0.018 1:-0.015 2:-0.031 3:-0.016 4:-0.010
      90W/36L/0T   worst degradation +0.1204 (8T62)   p90 +0.0274   power 0.99  Type-M 1.01
      concentration: drop-top10 -0.0088 vs uniform-effect null p10/p50/p90 -0.0141/-0.0092/-0.0040 -> pctile 0.543
      VERDICT: BETTER
    stage1 ORACLE cos: AMBER minus RANDOM
      a -0.0491 (med -0.0498)   b 0.0013 (med 0.0028)   n=126
      effect -0.0503   median -0.0509   SE 0.0154   MDE 0.0432   effect/MDE -1.17
      iid  CI95 [-0.0799, -0.0211]
      fold CI95 [-0.0735, -0.0253]   folds same sign 5/5   per-fold 0:-0.050 1:-0.088 2:-0.066 3:-0.001 4:-0.046
      74W/52L/0T   worst degradation +0.3845 (5MXS)   p90 +0.1647   power 0.90  Type-M 1.06
      concentration: drop-top10 -0.0225 vs uniform-effect null p10/p50/p90 -0.0433/-0.0232/-0.0028 -> pctile 0.517
      VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.06x]
    stage1 ORACLE cos: AMBER minus TOWARD-MEMBER
      a -0.0491 (med -0.0498)   b 0.1229 (med 0.1090)   n=126
      effect -0.1719   median -0.1733   SE 0.0232   MDE 0.0649   effect/MDE -2.65
      iid  CI95 [-0.2156, -0.1266]
      fold CI95 [-0.2104, -0.1388]   folds same sign 5/5   per-fold 0:-0.168 1:-0.190 2:-0.244 3:-0.121 4:-0.140
      94W/32L/0T   worst degradation +0.5848 (2NDN)   p90 +0.1612   power 1.00  Type-M 1.00
      concentration: drop-top10 -0.1320 vs uniform-effect null p10/p50/p90 -0.1626/-0.1327/-0.1020 -> pctile 0.512
      VERDICT: BETTER

    ORACLE cos(AMBER displacement, true residual)   mean -0.049 (SE 0.015), median -0.050, positive on 36.5% of targets
    orthogonal-move cost predicted by the S16 identity from the magnitude alone   +0.0107 (SE 0.0011), median +0.0079
    identity n RMSD^2 = |r|^2 - 2 v.r + |v|^2 as an upper bound on the superposed RMSD   max violation 0.0062 A (it is a bound, not an equality, under re-superposition)

Reading. (1) The production step reproduces exactly: +0.0207 [+0.0154, +0.0290], 2.16x MDE.
(2) About half of that cost is the SIZE of the move: a random displacement of the same 0.220 A
costs +0.0096, which is what the S16 identity predicts for a move orthogonal to the residual
(+0.0107). The other half is DIRECTION: AMBER is worse than its own random twin by +0.0111
(fold CI excluding zero, 5/5 folds, but 1.12x MDE, so the magnitude is an upper bound), and its
displacement points away from the native, cos -0.049 (36.5% positive), which reproduces S16 L27
to the fourth decimal (S16: cos -0.0521, AMBER minus random -0.0491; here AMBER minus random on
the cosine -0.0503). The registered prediction (+0.013, negative cosine) held. (3) A move of the
same size TOWARD A RANDOM POOL MEMBER improves the built chain by -0.0178 [-0.0255, -0.0124],
90W/36L, 1.49x MDE, 5/5 folds, with ORACLE cos +0.123. This is a zero-information control, not
a proposal, and it is PROVISIONAL until the replication registered in `s26/PREREG_c3_control.md`
addendum 2 lands (job `ph_c3_stage1_rep`: new seeds, reversed order). Its mechanism is on the
record: the projection moved the chain 0.166 A away from the point cloud with a negative cosine
(S16 L27), and the pool members surround the cloud, so any move back toward one of them
recovers part of that displacement. It says nothing about physics and everything about the
projection's cost; it will be handed to the Adversary as such. (4) Power: SEs 0.0015 to 0.0061,
MDEs 0.004 to 0.017; every contrast except AMBER-vs-random and the AMBER-vs-random cosine is
above 1.3x its MDE; those two are in the Type-M zone, direction measured, magnitude inflated
about 1.07x.

DECISION RULE (`s26/PREREG_c3_control.md` section 3): AMBER does not beat either matched
control; it is worse than both. "Refine with physics" is DROPPED as an accuracy step and KEPT
as a validity step. `s26/C3_RESULT.md` carries the sentence with the numbers for lane P and the
presentation. Stage 2 repeats this on lane P's best C2 rung when it is delivered.
