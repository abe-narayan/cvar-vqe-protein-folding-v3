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

---
## ADDENDUM 1 (2026-09-13 19:35) -- THE ADVERSARY'S CAVEATS ON L39 (ledger L46, STANDS WITH CAVEAT), FOR LANE P AND LANE PR

Quote the result this way and not otherwise.

1. The +0.0111 A "AMBER worse than a random move of its own size" is a TYPE-M number: 1.12x
   its own MDE, fold CI [+0.0062, +0.0171] excluding zero, 5/5 folds, so its SIGN is measured
   and its MAGNITUDE is inflated about 1.07x and is not a result as a magnitude. Any sentence
   that quotes +0.0111 carries the Type-M flag. The two clean contrasts, both well clear of the
   Type-M zone, are: AMBER minus do-nothing +0.0207 A (2.16x MDE, fold CI [+0.0154, +0.0290],
   40W/86L) and AMBER minus a same-size move toward a random pool member +0.0385 A (2.27x MDE,
   fold CI [+0.0298, +0.0479], 29W/97L). Lean on those two.

2. The decision rule holds regardless of the Type-M flag: "accuracy step" required AMBER to
   BEAT both matched controls, and it is worse than both. "Refine with physics" is not an
   accuracy step. Note also what the Adversary notes: a zero-information move of the same size
   TOWARD a random pool member IMPROVES the built chain (-0.0178, fold CI [-0.0255, -0.0124],
   90W/36L), where the physics move worsens it. That line is provisional until its registered
   replication (`ph_c3_stage1_rep`) lands, and it is a statement about the projection's cost
   (S16 L27), not about physics.

3. "Validity step" carries an exception: on 124 of 126 targets the relaxation converges with
   a sane virtual CA-CA bond (3.804 to 3.867 A on average); on 2BP4 it stretches a virtual bond
   to 5.38 A, and on 9KAR to 4.86 A while ending at +1262 kcal/mol, above the 1000 kcal/mol
   convergence gate. The sentence is: "a validity step on 124 of 126 targets; on 2 it breaks a
   virtual bond, one of which does not converge."

The presentation line, with these applied: the physics relaxation costs 0.021 A of accuracy
[+0.015, +0.029], more than a random displacement of the same size and 0.039 A more than a
same-size move toward any other pool member; it is kept because it turns a chain with an
energy above 1e4 kcal/mol on 59% of targets into one below 1000 kcal/mol on 125 of 126, at the
price of a broken virtual bond on two.

---
## ADDENDUM 2 (2026-09-13 23:05) -- THE REPLICATION LANDED

`s26/results/ph_c3_stage1_rep.json` (new draw seeds, reversed target order; ledger entry "C3
STAGE 1 REPLICATION"): every contrast inside the first run's fold CI. Toward-member minus
do-nothing -0.0207 [-0.0250, -0.0166], 94W/32L; AMBER minus random +0.0100 [+0.0059, +0.0175]
(still Type-M, 1.02x MDE); AMBER minus toward-member +0.0414 [+0.0331, +0.0501]. The
"provisional" label in the table above is lifted; nothing else changes.

---
## ADDENDUM 3 (2026-09-14 00:25) -- THE HEAVY-ATOM VALIDITY AXIS, THE NUMBER "VALIDITY STEP" RESTS ON

`s26/results/ph_validity.json` (126/126; the re-run reproduces the cache bit-for-bit on every
target), ledger entry "THE VALIDITY AXIS OF THE PRODUCTION RELAXATION". Built chain to relaxed
chain: heavy-atom pairs below 2.0 A 0.444 per target (34 targets) to 0.008 (1 target); contacts
below 2.6 A 3.45 to 0.12; closest pair 2.36 to 2.78 A; bond strain 0.0 to 1.3% (12% on 2BP4),
angle strain 0.0 to 2.5%, omega non-planarity 0.0 to 6.6 degrees (48 on 1D6X); Ramachandran
favoured 0.930 to 0.911 (not measured). The constant alpha-helix reference has zero clashes,
ideal geometry and rama 1.000, so the defensible claim is the conjunction: the relaxation
removes the builder's clashes while holding the torsions and moving the CA trace 0.220 A, at
a covalent price. The sentence for the presentation: "a validity step that turns 34 emissions
with a sub-2 A heavy-atom overlap into 1 and 125 of 126 energies above the 1000 kcal/mol gate
into converged ones, at the price of 0.021 A of accuracy, 1.3% bond and 2.5% angle strain,
6.6 degrees of peptide-bond non-planarity, and a broken virtual bond on 2BP4 and 9KAR."
