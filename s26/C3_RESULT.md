# C3 RESULT -- AMBER REFINEMENT AGAINST A MATCHED-MAGNITUDE RANDOM DISPLACEMENT (lane PH, S26)

Stage 1, on the PRODUCTION built chain, no new AMBER compute. `s26/results/ph_c3_stage1.json`
(complete 126/126), pre-registered in `s26/PREREG_c3_control.md`, ledger entry "C3 STAGE 1"
(2026-09-13, PH). Stage 2 (lane P's best C2 rung, AMBER compute) will be appended below when the
rung file `s26/results/p_best_rung_chains.json` exists.

## The sentence, for lane P and the presentation

The production relaxation (restrained ff14SB/GBn2, k = 10 kcal/mol/A^2, converged) moves the
built chain 0.220 A and costs +0.0207 A of CA-RMSD [fold CI +0.0154, +0.0290; 40W/86L; n = 126;
built chain 3.2148 A to relaxed chain 3.2355 A]. A random displacement of exactly the same size
costs +0.0096 A [+0.0077, +0.0122]. The relaxation is therefore WORSE than a random move of its
own size by +0.0111 A [fold CI +0.0062, +0.0171; 5/5 folds; 1.12x its MDE of 0.0099, so the size
of that number is an upper bound and its sign is what is measured], and worse than a move of the
same size toward a random pool member by +0.0385 A [+0.0298, +0.0479; 2.27x MDE]. Its
displacement points slightly away from the native (ORACLE cosine -0.049, positive on 36.5% of
targets), reproducing S16 L27 (-0.052; paired against random -0.0491 there, -0.0503 here).

By the decision rule registered before the run, "refine with physics" is dropped as an accuracy
step and kept as a validity step. What the validity step buys, with numbers
(`s26/results/ph_c3_nativefree.json`, ledger L24): the built chain's own AMBER energy is above
1e4 kcal/mol on 58.7% of targets before relaxation (median 8.6e4); the relaxation brings 125 of
126 below +1000 kcal/mol (mean -560; 9KAR ends at +1262), at a cost of 0.220 A of movement, a
0.06 A stretch of the virtual CA-CA bond (3.804 to 3.867 A) and a broken virtual bond on two
targets (2BP4 5.38 A, 9KAR 4.86 A). The heavy-atom validity axis (clash count, bond deviation)
of the relaxed emission is measured in stage 2.

## The numbers (all `ST.compare`, paired per target, fold-clustered CI beside iid)

| contrast (basis a minus basis b) | effect | SE | MDE | x MDE | iid CI | fold CI | W/L | verdict |
|---|---|---|---|---|---|---|---|---|
| AMBER (relaxed chain) minus do-nothing (built chain) | +0.0207 | 0.0034 | 0.0096 | 2.16 | [+0.0141, +0.0274] | [+0.0154, +0.0290] | 40/86 | WORSE |
| AMBER minus random, same magnitude (16 draws) | +0.0111 | 0.0035 | 0.0099 | 1.12 | [+0.0044, +0.0178] | [+0.0062, +0.0171] | 52/74 | WORSE, Type-M zone |
| AMBER minus toward-member, same magnitude (16 draws) | +0.0385 | 0.0061 | 0.0170 | 2.27 | [+0.0272, +0.0506] | [+0.0298, +0.0479] | 29/97 | WORSE |
| random minus do-nothing | +0.0096 | 0.0015 | 0.0042 | 2.28 | [+0.0070, +0.0127] | [+0.0077, +0.0122] | 29/97 | WORSE |
| toward-member minus do-nothing | -0.0178 | 0.0043 | 0.0120 | 1.49 | [-0.0266, -0.0096] | [-0.0255, -0.0124] | 90/36 | BETTER, provisional until replicated |
| ORACLE cos: AMBER minus random | -0.0503 | 0.0154 | 0.0432 | 1.17 | [-0.0799, -0.0211] | [-0.0735, -0.0253] | 74/52 | AMBER points further from the native, Type-M zone |
| ORACLE cos: AMBER minus toward-member | -0.1719 | 0.0232 | 0.0649 | 2.65 | [-0.2156, -0.1266] | [-0.2104, -0.1388] | 94/32 | same, measured |

The controls displace the built chain's CA trace by AMBER's own per-atom RMS magnitude (0.220 A
mean, measured after superposition); they match the operator's space, not its validity (a
displaced trace is not an ideal-geometry chain), as S16 did. The toward-member line is a
zero-information control that improves the chain: the projection moved the chain away from the
point cloud (S16 L27) and the pool members surround the cloud, so any move back toward one
recovers part of that cost. It is not a proposal and it waits for its replication
(`ph_c3_stage1_rep`, new seeds, reversed order) before it is more than provisional.

## What stage 2 adds

The same seven contrasts on lane P's best C2 rung output (the production relaxation re-run under
`--tag AMBER`, one-target probe first, per-target cells), plus the heavy-atom validity panel
(`s16.energy_lib.panel`) of the input and the relaxed output.
